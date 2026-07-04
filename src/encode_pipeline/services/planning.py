"""Execution planning service boundary."""

from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING
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
    """Platform-safe workspace planning boundary.

    ``WorkspacePlanner`` consumes an immutable ``ExecutionPlan`` and returns a
    fresh ``ExecutionPlan`` with a populated ``workspace_plan``. It does not
    touch the filesystem, mutate run state, call adapter planning methods, or
    build commands.
    """

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
                        message=str(exc),
                        severity="error",
                        path="base_dir",
                        source="workspace_planner",
                    )
                ]
            )

        directories = ("logs", "results")
        files: tuple[tuple[str, bytes], ...] = ()

        try:
            for directory in directories:
                policy.resolve(directory)
            for file_path, _ in files:
                policy.resolve(file_path)
        except WorkspacePathError as exc:
            return Result.failure(
                [
                    Issue(
                        code=exc.code,
                        message=str(exc),
                        severity="error",
                        path="workspace_plan",
                        source="workspace_planner",
                    )
                ]
            )

        workspace_plan = WorkspacePlan(
            directories=directories,
            files=files,
        )

        updated_plan = ExecutionPlan(
            plan_id=str(uuid4()),
            run_id=plan.run_id,
            workflow_id=plan.workflow_id,
            status=PlanStatus.PENDING,
            inputs_snapshot=plan.inputs_snapshot,
            dag_preview=plan.dag_preview,
            workspace_plan=workspace_plan,
            command_spec=None,
            created_at=datetime.now(timezone.utc),
            issues=plan.issues
            + (
                Issue(
                    code="WORKSPACE_PLANNING_COMPLETE",
                    message="Workspace planning completed.",
                    severity="info",
                    path="workspace_plan",
                    source="workspace_planner",
                ),
            ),
        )
        return Result.success(updated_plan, issues=updated_plan.issues)
