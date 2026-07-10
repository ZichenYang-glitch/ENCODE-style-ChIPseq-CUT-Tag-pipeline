"""Execution planning service boundary."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Any
from collections.abc import Mapping
from uuid import uuid4

from encode_pipeline.platform.adapters import WorkspacePlan
from encode_pipeline.platform.planning import (
    ExecutionPlan,
    PlanStatus,
    WorkspacePathError,
    WorkspacePathPolicy,
)
from encode_pipeline.platform.results import Issue, Result

if TYPE_CHECKING:
    from encode_pipeline.services.runs import RunService


class ExecutionPlanner:
    """Read-only planning boundary for workflow runs.

    ``ExecutionPlanner`` never mutates run lifecycle state, emits events,
    appends logs, records artifacts, or calls adapter planning methods. PR106
    returns an ``UNSUPPORTED`` plan for every run that exists.
    """

    def __init__(self, run_service: "RunService") -> None:
        from encode_pipeline.services.runs import RunService

        if not isinstance(run_service, RunService):
            raise ValueError("ExecutionPlanner requires a RunService instance")
        self._run_service = run_service

    def plan_run(self, run_id: str) -> Result[ExecutionPlan]:
        """Return an execution plan for the given run ID."""
        try:
            record = self._run_service.get_run(run_id)
        except KeyError:
            return Result.failure(
                [
                    Issue(
                        code="EXECUTION_RUN_NOT_FOUND",
                        message="Run not found.",
                        severity="error",
                        path=run_id,
                        source="execution_planner",
                    )
                ]
            )

        plan = ExecutionPlan(
            plan_id=str(uuid4()),
            run_id=record.run_id,
            workflow_id=record.workflow_id,
            status=PlanStatus.UNSUPPORTED,
            inputs_snapshot=record.inputs,
            created_at=datetime.now(timezone.utc),
            issues=(
                Issue(
                    code="EXECUTION_PLANNING_UNSUPPORTED",
                    message="Execution planning is not supported yet.",
                    severity="info",
                    path="execution_plan",
                    source="execution_planner",
                ),
            ),
        )
        return Result.success(plan, issues=plan.issues)


class WorkspacePlanner:
    """Platform-safe workspace planning boundary that delegates to adapters.

    ``WorkspacePlanner`` consumes an immutable ``ExecutionPlan``, reconstructs
    the original ``WorkflowInputs``, looks up the workflow adapter, and returns
    a fresh ``ExecutionPlan`` with a populated ``workspace_plan``. It validates
    adapter-returned paths with ``WorkspacePathPolicy`` but does not touch the
    filesystem, mutate run state, or build commands.
    """

    def __init__(self, registry: WorkflowRegistry) -> None:
        from encode_pipeline.platform.registry import WorkflowRegistry

        if not isinstance(registry, WorkflowRegistry):
            raise ValueError("WorkspacePlanner requires a WorkflowRegistry instance")
        self._registry = registry

    @staticmethod
    def _reconstruct_inputs(
        inputs_snapshot: Mapping[str, Any],
    ) -> Result[WorkflowInputs]:
        """Rebuild WorkflowInputs from an ExecutionPlan inputs snapshot."""
        from encode_pipeline.platform.adapters import WorkflowInputs
        from encode_pipeline.platform.results import Issue

        def _malformed(path: str) -> Result[WorkflowInputs]:
            return Result.failure(
                [
                    Issue(
                        code="WORKSPACE_PLAN_MALFORMED_INPUTS",
                        message="Execution plan inputs snapshot is malformed and cannot be reconstructed into workflow inputs.",
                        severity="error",
                        path=path,
                        source="workspace_planner",
                    )
                ]
            )

        if not isinstance(inputs_snapshot, Mapping):
            return _malformed("inputs_snapshot")

        config = inputs_snapshot.get("config")
        if not isinstance(config, Mapping):
            return _malformed("inputs_snapshot.config")

        samples = inputs_snapshot.get("samples")
        if samples is not None:
            if isinstance(samples, Path):
                samples = str(samples)
            if isinstance(samples, str):
                pass
            elif isinstance(samples, list):
                for index, row in enumerate(samples):
                    if not isinstance(row, Mapping):
                        return _malformed(f"inputs_snapshot.samples[{index}]")
                    if not all(
                        isinstance(k, str) and isinstance(v, str)
                        for k, v in row.items()
                    ):
                        return _malformed(f"inputs_snapshot.samples[{index}]")
            else:
                return _malformed("inputs_snapshot.samples")

        options = inputs_snapshot.get("options")
        if options is None:
            options = {}
        elif not isinstance(options, Mapping):
            return _malformed("inputs_snapshot.options")
        else:
            options = dict(options)

        try:
            return Result.success(
                WorkflowInputs(
                    config=config,
                    samples=samples,
                    options=options,
                )
            )
        except ValueError:
            return _malformed("inputs_snapshot")

    _WORKSPACE_PATH_ERROR_MESSAGES = {
        "WORKSPACE_BASE_DIR_RELATIVE": "Workspace base directory must be absolute.",
        "WORKSPACE_PATH_ABSOLUTE": "Planned workspace path must be relative.",
        "WORKSPACE_PATH_TRAVERSAL": "Planned workspace path must not contain parent references.",
        "WORKSPACE_PATH_INVALID": "Planned workspace path is invalid.",
        "WORKSPACE_PATH_ESCAPE": "Planned workspace path escapes the workspace directory.",
    }

    def plan_workspace(
        self,
        plan: ExecutionPlan,
        base_dir: Path,
    ) -> Result[ExecutionPlan]:
        """Return a new ExecutionPlan with a safe workspace plan."""
        if not isinstance(plan, ExecutionPlan):
            return Result.failure(
                [
                    Issue(
                        code="WORKSPACE_PLAN_INVALID",
                        message="plan must be an ExecutionPlan",
                        severity="error",
                        path="plan",
                        source="workspace_planner",
                    )
                ]
            )

        try:
            policy = WorkspacePathPolicy(base_dir=base_dir)
        except WorkspacePathError as exc:
            return Result.failure(
                [
                    Issue(
                        code=exc.code,
                        message=self._WORKSPACE_PATH_ERROR_MESSAGES.get(
                            exc.code, "Workspace path policy rejected the base directory."
                        ),
                        severity="error",
                        path="base_dir",
                        source="workspace_planner",
                    )
                ]
            )

        inputs_result = self._reconstruct_inputs(plan.inputs_snapshot)
        if inputs_result.is_failure:
            return inputs_result
        inputs = inputs_result.value

        try:
            adapter = self._registry.get(plan.workflow_id)
        except KeyError:
            return Result.failure(
                [
                    Issue(
                        code="WORKSPACE_PLAN_WORKFLOW_NOT_FOUND",
                        message="Workflow is not supported.",
                        severity="error",
                        path="workflow_id",
                        source="workspace_planner",
                    )
                ]
            )

        adapter_result = adapter.plan_workspace(inputs, base_dir)
        if adapter_result.is_failure:
            return adapter_result
        adapter_plan = adapter_result.value

        for index, directory in enumerate(adapter_plan.directories):
            try:
                policy.resolve(directory)
            except WorkspacePathError as exc:
                return Result.failure(
                    [
                        Issue(
                            code=exc.code,
                            message=self._WORKSPACE_PATH_ERROR_MESSAGES.get(
                                exc.code, "Planned workspace path violates workspace path policy."
                            ),
                            severity="error",
                            path=f"workspace_plan.directories[{index}]",
                            source="workspace_planner",
                        )
                    ]
                )

        for index, (file_path, _) in enumerate(adapter_plan.files):
            try:
                policy.resolve(file_path)
            except WorkspacePathError as exc:
                return Result.failure(
                    [
                        Issue(
                            code=exc.code,
                            message=self._WORKSPACE_PATH_ERROR_MESSAGES.get(
                                exc.code, "Planned workspace path violates workspace path policy."
                            ),
                            severity="error",
                            path=f"workspace_plan.files[{index}]",
                            source="workspace_planner",
                        )
                    ]
                )

        updated_plan = ExecutionPlan(
            plan_id=str(uuid4()),
            run_id=plan.run_id,
            workflow_id=plan.workflow_id,
            status=PlanStatus.PENDING,
            inputs_snapshot=plan.inputs_snapshot,
            dag_preview=plan.dag_preview,
            workspace_plan=adapter_plan,
            command_spec=None,
            created_at=datetime.now(timezone.utc),
            issues=plan.issues + adapter_result.issues,
        )
        return Result.success(updated_plan, issues=updated_plan.issues)
