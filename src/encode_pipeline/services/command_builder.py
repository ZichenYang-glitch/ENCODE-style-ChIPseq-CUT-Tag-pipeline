"""Pure command-spec construction boundary for planned workflow runs."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping
from uuid import uuid4

from encode_pipeline.platform.adapters import CommandSpec, WorkspacePlan
from encode_pipeline.platform.planning import ExecutionPlan, PlanStatus, WorkspacePathPolicy
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result


def _bundled_snakefile_path() -> Path:
    """Return the controlled repo-relative path to the bundled Snakefile."""
    return Path(__file__).resolve().parents[3] / "workflow" / "Snakefile"


class CommandBuilder:
    """Pure command-spec construction boundary for planned workflow runs."""

    def __init__(self, registry: WorkflowRegistry) -> None:
        """Initialize with an adapter registry for engine validation."""
        if not isinstance(registry, WorkflowRegistry):
            raise ValueError("registry must be a WorkflowRegistry")
        self._registry = registry

    def build_command(
        self,
        plan: ExecutionPlan,
        base_dir: Path,
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

        if "snakemake" not in adapter.metadata.engines:
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

        snakefile = _bundled_snakefile_path()
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

        command_spec = CommandSpec(
            argv=(
                "snakemake",
                "--snakefile",
                str(snakefile),
                "--directory",
                str(base_dir),
                "--configfile",
                str(config_path),
                "--cores",
                str(cores),
            ),
            cwd=None,
            env={},
        )

        new_plan = ExecutionPlan(
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

        return Result.success(new_plan)

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
