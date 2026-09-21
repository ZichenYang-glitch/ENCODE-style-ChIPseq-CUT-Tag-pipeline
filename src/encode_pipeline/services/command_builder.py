"""Pure command-spec construction boundary for planned workflow runs."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import Protocol
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


class CommandBuilder:
    """Pure command-spec construction boundary for planned workflow runs."""

    def __init__(
        self,
        registry: WorkflowRegistry,
        *,
        reference_profile_resolver: _ReferenceProfileRuntimeResolver | None = None,
    ) -> None:
        """Initialize with an adapter registry for engine validation."""
        if not isinstance(registry, WorkflowRegistry):
            raise ValueError("registry must be a WorkflowRegistry")
        if reference_profile_resolver is not None and not callable(
            getattr(reference_profile_resolver, "resolve_run", None)
        ):
            raise ValueError("reference_profile_resolver is invalid")
        self._registry = registry
        self._reference_profile_resolver = reference_profile_resolver

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

        trusted = self._registry.uses_encode_execution_fallback(registered_adapter)
        same_contract = (
            type(adapter) is type(registered_adapter)
            and adapter.metadata == registered_adapter.metadata
            and adapter.capabilities == registered_adapter.capabilities
        )
        if trusted and same_contract:
            return self._build_trusted_adapter_command(adapter, plan, base_dir)
        if COMMAND_CAPABILITY in adapter.capabilities.supports:
            return self._build_adapter_command(
                adapter=adapter,
                plan=plan,
                base_dir=base_dir,
            )
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

    def _build_trusted_adapter_command(
        self,
        adapter: WorkflowAdapter,
        plan: ExecutionPlan,
        base_dir: Path,
    ) -> Result[ExecutionPlan]:
        """Delegate under existing exact-instance legacy authority (cwd may be None)."""
        assert plan.workspace_plan is not None
        try:
            result = adapter.build_command(plan.workspace_plan, base_dir)
        except Exception:
            return self._adapter_failure()
        if not isinstance(result, Result):
            return self._adapter_failure()
        if result.is_failure:
            return Result.failure(result.issues)
        if not isinstance(result.value, CommandSpec):
            return self._adapter_failure()
        return Result.success(self._planned_plan(plan, result.value))

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
