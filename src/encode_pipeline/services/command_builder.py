"""Pure command-spec construction boundary for planned workflow runs."""

from __future__ import annotations

from datetime import datetime, timezone
import os
from pathlib import Path
from typing import Any, Mapping, Protocol
from uuid import uuid4

from encode_pipeline.platform.adapters import (
    COMMAND_CAPABILITY,
    CommandSpec,
    WorkflowAdapter,
    WorkflowInputs,
    WorkspacePlan,
)
from encode_pipeline.platform.managed_containers import managed_container_scope
from encode_pipeline.platform.planning import (
    ExecutionPlan,
    PlanStatus,
    WorkspacePathPolicy,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.reference_profiles import BoundWorkflowReference
from encode_pipeline.platform.results import Issue, Result


class _ReferenceProfileRuntimeResolver(Protocol):
    def resolve_run(
        self,
        run_id: str,
        workflow_id: str,
        inputs: WorkflowInputs,
        *,
        require_enabled: bool,
    ) -> Result[BoundWorkflowReference | None]: ...


def _bundled_snakefile_path(project_root: Path | None = None) -> Path:
    """Return the controlled project-relative path to the bundled Snakefile."""
    root = Path(__file__).resolve().parents[3] if project_root is None else project_root
    return root / "workflow" / "Snakefile"


class CommandBuilder:
    """Pure command-spec construction boundary for planned workflow runs."""

    def __init__(
        self,
        registry: WorkflowRegistry,
        *,
        project_root: Path | None = None,
        reference_profile_resolver: _ReferenceProfileRuntimeResolver | None = None,
        snakemake_executable: Path | None = None,
        conda_prefix: Path | None = None,
    ) -> None:
        """Initialize with an adapter registry for engine validation."""
        if not isinstance(registry, WorkflowRegistry):
            raise ValueError("registry must be a WorkflowRegistry")
        root = (
            Path(__file__).resolve().parents[3]
            if project_root is None
            else project_root
        )
        if not isinstance(root, Path) or not root.is_absolute():
            raise ValueError("project_root must be an absolute pathlib.Path")
        if reference_profile_resolver is not None and not callable(
            getattr(reference_profile_resolver, "resolve_run", None)
        ):
            raise ValueError("reference_profile_resolver is invalid")
        if (snakemake_executable is None) != (conda_prefix is None):
            raise ValueError(
                "snakemake_executable and conda_prefix must be configured together"
            )
        if snakemake_executable is not None:
            if (
                not isinstance(snakemake_executable, Path)
                or not snakemake_executable.is_absolute()
                or snakemake_executable.name != "snakemake"
                or any(
                    character in str(snakemake_executable)
                    for character in ("\x00", "\n", "\r")
                )
            ):
                raise ValueError("snakemake_executable is invalid")
            if (
                not isinstance(conda_prefix, Path)
                or not conda_prefix.is_absolute()
                or any(
                    character in str(conda_prefix) for character in ("\x00", "\n", "\r")
                )
            ):
                raise ValueError("conda_prefix is invalid")
        self._registry = registry
        self._project_root = root
        self._reference_profile_resolver = reference_profile_resolver
        self._snakemake_executable = snakemake_executable
        self._conda_prefix = conda_prefix

    def build_command(
        self,
        plan: ExecutionPlan,
        base_dir: Path,
        *,
        require_reference_enabled: bool = True,
    ) -> Result[ExecutionPlan]:
        """Build a controlled CommandSpec for ``plan`` under ``base_dir``."""
        if not isinstance(plan, ExecutionPlan):
            return Result.failure(
                [
                    Issue(
                        code="COMMAND_BUILD_INVALID_PLAN",
                        message="plan must be an ExecutionPlan.",
                        severity="error",
                        path="plan",
                        source="command_builder",
                    )
                ]
            )

        if not isinstance(base_dir, Path) or not base_dir.is_absolute():
            return Result.failure(
                [
                    Issue(
                        code="COMMAND_BUILD_BASE_DIR_RELATIVE",
                        message="base_dir must be an absolute Path.",
                        severity="error",
                        path="base_dir",
                        source="command_builder",
                    )
                ]
            )

        if plan.status is not PlanStatus.PENDING:
            return Result.failure(
                [
                    Issue(
                        code="COMMAND_BUILD_INVALID_PLAN_STATUS",
                        message="plan.status must be PENDING.",
                        severity="error",
                        path="plan",
                        source="command_builder",
                    )
                ]
            )

        if plan.workspace_plan is None:
            return Result.failure(
                [
                    Issue(
                        code="COMMAND_BUILD_MISSING_WORKSPACE_PLAN",
                        message="plan.workspace_plan is required.",
                        severity="error",
                        path="workspace_plan",
                        source="command_builder",
                    )
                ]
            )

        try:
            adapter = self._registry.get(plan.workflow_id)
        except KeyError:
            return Result.failure(
                [
                    Issue(
                        code="COMMAND_BUILD_UNSUPPORTED_WORKFLOW",
                        message="Workflow is not supported.",
                        severity="error",
                        path="workflow",
                        source="command_builder",
                    )
                ]
            )
        registered_adapter = adapter

        resolved_adapter = self._resolve_reference_profile(
            plan,
            adapter,
            require_enabled=require_reference_enabled,
        )
        if resolved_adapter.is_failure:
            return Result.failure(resolved_adapter.issues)
        assert resolved_adapter.value is not None
        adapter = resolved_adapter.value

        if COMMAND_CAPABILITY in adapter.capabilities.supports:
            return self._build_adapter_command(
                adapter=adapter,
                plan=plan,
                base_dir=base_dir,
            )

        if not (
            self._registry.uses_encode_execution_fallback(adapter)
            or self._registry.uses_encode_execution_fallback(registered_adapter)
        ):
            return Result.failure(
                [
                    Issue(
                        code="COMMAND_BUILD_UNSUPPORTED_ENGINE",
                        message="Workflow engine is not supported.",
                        severity="error",
                        path="workflow",
                        source="command_builder",
                    )
                ]
            )

        snakefile = _bundled_snakefile_path(self._project_root)
        if not snakefile.is_file():
            return Result.failure(
                [
                    Issue(
                        code="COMMAND_BUILD_SNAKEFILE_NOT_FOUND",
                        message="Bundled Snakefile was not found.",
                        severity="error",
                        path="workflow",
                        source="command_builder",
                    )
                ]
            )

        config_path = self._resolve_config_path(base_dir, plan.workspace_plan)
        if isinstance(config_path, Result):
            return config_path

        cores_result = self._resolve_cores(plan.inputs_snapshot)
        if isinstance(cores_result, Result):
            return cores_result
        cores = cores_result

        runtime = self._scientific_runtime(base_dir)
        if isinstance(runtime, Result):
            return runtime
        executable, runtime_arguments, environment = runtime
        argv = (
            executable,
            "--snakefile",
            str(snakefile),
            "--directory",
            str(base_dir),
            "--configfile",
            str(config_path),
            "--cores",
            str(cores),
            *runtime_arguments,
        )
        command_spec = CommandSpec(
            argv=argv,
            cwd=None,
            env=environment,
            preflight_argv=argv + ("-n",),
        )

        return Result.success(self._planned_plan(plan, command_spec))

    def _scientific_runtime(
        self,
        workspace: Path,
    ) -> tuple[str, tuple[str, ...], dict[str, str]] | Result[ExecutionPlan]:
        executable = self._snakemake_executable
        prefix = self._conda_prefix
        if executable is None or prefix is None:
            return "snakemake", (), {}
        runtime_root = executable.parent.parent.parent
        mamba_root = runtime_root / "mamba-root"
        conda_executable = executable.parent / "conda"
        activate = mamba_root / "bin" / "activate"
        micromamba = runtime_root / "runner" / "libexec" / "micromamba"
        try:
            observed_files = tuple(
                (path, path.lstat())
                for path in (executable, conda_executable, activate, micromamba)
            )
            observed_directories = tuple(
                (path, path.lstat()) for path in (runtime_root, mamba_root, prefix)
            )
            if (
                executable.parent.parent != runtime_root / "runner"
                or prefix != runtime_root / "conda-envs"
                or any(
                    path.is_symlink()
                    or not path.is_file()
                    or not os.access(path, os.X_OK)
                    or witness.st_mode & 0o022
                    for path, witness in observed_files
                )
                or any(
                    path.is_symlink() or not path.is_dir() or witness.st_mode & 0o022
                    for path, witness in observed_directories
                )
            ):
                raise OSError
        except OSError:
            return Result.failure(
                [
                    Issue(
                        code="COMMAND_BUILD_SCIENTIFIC_RUNTIME_UNAVAILABLE",
                        message="The admitted scientific runtime is unavailable.",
                        severity="error",
                        path="workflow",
                        source="command_builder",
                    )
                ]
            )
        path = ":".join(
            (
                str(executable.parent),
                "/usr/sbin",
                "/usr/bin",
                "/sbin",
                "/bin",
            )
        )
        return (
            str(executable),
            (
                "--use-conda",
                "--conda-prefix",
                str(prefix),
                "--conda-base-path",
                str(mamba_root),
                "--conda-frontend",
                "conda",
            ),
            {
                "PATH": path,
                "CONDA_DEFAULT_ENV": "",
                "CONDA_EXE": str(conda_executable),
                "CONDA_PREFIX": "",
                "CONDA_SHLVL": "0",
                "HOME": str(workspace),
                "MAMBA_ROOT_PREFIX": str(mamba_root),
                "PYTHONDONTWRITEBYTECODE": "1",
                "PYTHONNOUSERSITE": "1",
                "TMPDIR": str(workspace),
                "XDG_CACHE_HOME": str(workspace / ".snakemake"),
                "_CONDA_EXE": str(conda_executable),
                "_CONDA_ROOT": str(mamba_root),
            },
        )

    def _resolve_reference_profile(
        self,
        plan: ExecutionPlan,
        adapter: WorkflowAdapter,
        *,
        require_enabled: bool,
    ) -> Result[WorkflowAdapter]:
        resolver = self._reference_profile_resolver
        if resolver is None:
            return Result.success(adapter)
        from encode_pipeline.services.planning import WorkspacePlanner

        inputs_result = WorkspacePlanner._reconstruct_inputs(plan.inputs_snapshot)
        if inputs_result.is_failure:
            return self._reference_failure()
        assert inputs_result.value is not None
        try:
            resolved = resolver.resolve_run(
                plan.run_id,
                plan.workflow_id,
                inputs_result.value,
                require_enabled=require_enabled,
            )
        except Exception:
            return self._reference_failure()
        if not isinstance(resolved, Result):
            return self._reference_failure()
        if resolved.is_failure:
            return Result.failure(resolved.issues)
        bound = resolved.value
        if bound is None:
            return Result.success(adapter)
        if (
            not isinstance(bound, BoundWorkflowReference)
            or not isinstance(bound.inputs, WorkflowInputs)
            or not isinstance(bound.adapter, WorkflowAdapter)
            or bound.identity.workflow_id != plan.workflow_id
            or bound.adapter.metadata.workflow_id != plan.workflow_id
        ):
            return self._reference_failure()
        return Result.success(bound.adapter)

    @staticmethod
    def _reference_failure() -> Result[WorkflowAdapter]:
        return Result.failure(
            [
                Issue(
                    code="REFERENCE_PROFILE_BINDING_INVALID",
                    message=(
                        "The selected Reference Profile binding could not be verified."
                    ),
                    severity="error",
                    path="reference_profile_revision_id",
                    source="command_builder",
                )
            ]
        )

    def _build_adapter_command(
        self,
        *,
        adapter: WorkflowAdapter,
        plan: ExecutionPlan,
        base_dir: Path,
    ) -> Result[ExecutionPlan]:
        """Delegate command construction without exposing adapter failures."""
        assert plan.workspace_plan is not None
        try:
            adapter_result = adapter.build_command(plan.workspace_plan, base_dir)
        except Exception:
            return self._adapter_failure()

        if not isinstance(adapter_result, Result) or adapter_result.is_failure:
            return self._adapter_failure()
        command_spec = adapter_result.value
        if not isinstance(command_spec, CommandSpec):
            return self._adapter_failure()
        if not self._command_workspace_is_safe(
            command_spec,
            base_dir,
            plan.workspace_plan,
        ):
            return self._adapter_failure()
        return Result.success(self._planned_plan(plan, command_spec))

    @staticmethod
    def _command_workspace_is_safe(
        command_spec: CommandSpec,
        base_dir: Path,
        workspace_plan: WorkspacePlan,
    ) -> bool:
        """Require adapter-selected runtime paths to remain in the workspace."""
        if command_spec.cwd is None:
            return False
        if (
            command_spec.preflight_argv is not None
            and command_spec.preflight_argv[0] != command_spec.argv[0]
        ):
            return False
        cwd = Path(command_spec.cwd)
        if not cwd.is_absolute() or any(part in {".", ".."} for part in cwd.parts):
            return False
        try:
            cwd.relative_to(base_dir)
        except ValueError:
            return False
        if (
            command_spec.managed_container_scope is not None
            and command_spec.managed_container_scope
            != managed_container_scope(base_dir)
        ):
            return False
        planned_files = {
            str(base_dir / relative_path)
            for relative_path, _contents in workspace_plan.files
        }
        managed_paths = tuple(
            path
            for _stream_name, path in (
                *command_spec.preflight_managed_logs,
                *command_spec.execution_managed_logs,
            )
        )
        for path_value in managed_paths:
            path = Path(path_value)
            try:
                relative = path.relative_to(base_dir)
            except ValueError:
                return False
            if not relative.parts or str(path) in planned_files:
                return False
        return True

    @staticmethod
    def _adapter_failure() -> Result[ExecutionPlan]:
        return Result.failure(
            [
                Issue(
                    code="COMMAND_BUILD_ADAPTER_FAILED",
                    message="Workflow command could not be built.",
                    severity="error",
                    path="command_spec",
                    source="command_builder",
                )
            ]
        )

    @staticmethod
    def _planned_plan(
        plan: ExecutionPlan,
        command_spec: CommandSpec,
    ) -> ExecutionPlan:
        return ExecutionPlan(
            plan_id=str(uuid4()),
            run_id=plan.run_id,
            workflow_id=plan.workflow_id,
            status=PlanStatus.PLANNED,
            inputs_snapshot=plan.inputs_snapshot,
            dag_preview=plan.dag_preview,
            workspace_plan=plan.workspace_plan,
            command_spec=command_spec,
            created_at=datetime.now(timezone.utc),
            issues=(
                *plan.issues,
                Issue(
                    code="COMMAND_BUILDING_COMPLETE",
                    message="CommandSpec built successfully.",
                    severity="info",
                    path="command_spec",
                    source="command_builder",
                ),
            ),
        )

    def _resolve_config_path(
        self,
        base_dir: Path,
        workspace_plan: WorkspacePlan,
    ) -> Path | Result[Any]:
        """Find config/config.yaml in the workspace plan and resolve it safely."""
        for index, (file_path, _) in enumerate(workspace_plan.files):
            if file_path == "config/config.yaml":
                policy = WorkspacePathPolicy(base_dir=base_dir)
                try:
                    return policy.resolve(file_path)
                except Exception:
                    return Result.failure(
                        [
                            Issue(
                                code="COMMAND_BUILD_MISSING_CONFIG",
                                message="Workspace config file could not be resolved.",
                                severity="error",
                                path=f"workspace_plan.files[{index}]",
                                source="command_builder",
                            )
                        ]
                    )

        return Result.failure(
            [
                Issue(
                    code="COMMAND_BUILD_MISSING_CONFIG",
                    message="Workspace plan must include config/config.yaml.",
                    severity="error",
                    path="workspace_plan.files",
                    source="command_builder",
                )
            ]
        )

    def _resolve_cores(
        self,
        inputs_snapshot: Mapping[str, Any],
    ) -> int | Result[Any]:
        """Read a positive integer cores value from inputs_snapshot options."""
        options = inputs_snapshot.get("options")
        if not isinstance(options, Mapping):
            return 1
        cores = options.get("cores")
        if cores is None:
            return 1
        if type(cores) is not int or cores <= 0:
            return Result.failure(
                [
                    Issue(
                        code="COMMAND_BUILD_INVALID_CORES",
                        message="cores must be a positive integer.",
                        severity="error",
                        path="plan.inputs_snapshot.options.cores",
                        source="command_builder",
                    )
                ]
            )
        return cores
