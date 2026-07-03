"""Execution planning service boundary."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import TYPE_CHECKING
from uuid import uuid4

from encode_pipeline.platform.planning import ExecutionPlan, PlanStatus
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
