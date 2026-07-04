"""Fail-closed local execution driver skeleton."""

from __future__ import annotations

from typing import TYPE_CHECKING

from encode_pipeline.platform.planning import PlanStatus
from encode_pipeline.platform.results import Issue, Result

if TYPE_CHECKING:
    from encode_pipeline.platform.planning import ExecutionPlan
    from encode_pipeline.platform.runs import RunRecord
    from encode_pipeline.services.runs import RunService


class LocalRunDriver:
    """Fail-closed local execution driver skeleton.

    PR108 validates the plan/run association and refuses execution with
    documented issue codes. It does not invoke subprocess, Snakemake, or
    filesystem operations.
    """

    def __init__(self, run_service: "RunService") -> None:
        from encode_pipeline.services.runs import RunService

        if not isinstance(run_service, RunService):
            raise ValueError("LocalRunDriver requires a RunService instance")
        self._run_service = run_service

    def run(self, run_id: str, plan: "ExecutionPlan") -> "Result[RunRecord]":
        """Validate plan and refuse execution."""
        record = self._run_service.get_run(run_id)

        if record.run_id != plan.run_id or record.workflow_id != plan.workflow_id:
            issue = Issue(
                code="LOCAL_RUN_PLAN_MISMATCH",
                message="Plan run_id or workflow_id does not match the run record.",
                severity="error",
                path="plan",
                source="local_run_driver",
            )
            return self._refuse(issue, plan, record.run_id)

        if not plan.can_execute:
            if plan.status is not PlanStatus.PLANNED:
                issue = Issue(
                    code="LOCAL_RUN_NOT_PLANNED",
                    message="Execution plan is not in the planned state.",
                    severity="error",
                    path="plan",
                    source="local_run_driver",
                )
            else:
                issue = Issue(
                    code="LOCAL_RUN_MISSING_COMMAND_SPEC",
                    message="Execution plan has no command spec.",
                    severity="error",
                    path="plan",
                    source="local_run_driver",
                )
            return self._refuse(issue, plan, record.run_id)

        issue = Issue(
            code="LOCAL_RUN_NOT_IMPLEMENTED",
            message="Local execution is not implemented yet.",
            severity="error",
            path="plan",
            source="local_run_driver",
        )
        return self._refuse(issue, plan, record.run_id)

    def _refuse(
        self,
        issue: Issue,
        plan: "ExecutionPlan",
        run_id: str,
    ) -> "Result[RunRecord]":
        """Record a runner_refused event and return a failure result."""
        self._run_service.add_event(
            run_id=run_id,
            event_type="runner_refused",
            message="Local execution refused.",
            status=None,
            context={
                "reason_code": issue.code,
                "plan_status": plan.status.value,
                "can_execute": plan.can_execute,
                "has_command_spec": plan.command_spec is not None,
            },
        )
        return Result.failure([issue])
