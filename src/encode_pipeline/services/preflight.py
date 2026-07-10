"""Local preflight orchestrator for workflow runs."""

from __future__ import annotations

from typing import TYPE_CHECKING

from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.runs import RunRecord, RunStatus

if TYPE_CHECKING:
    from encode_pipeline.services.local_run_driver import LocalRunDriver
    from encode_pipeline.services.planning import ExecutionPlanner, WorkspacePlanner
    from encode_pipeline.services.runs import RunService


class LocalPreflightService:
    """Compose planning, workspace authoring, and dry-run into a preflight path.

    The service is synchronous so it can be unit-tested without HTTP. FastAPI
    routes should call ``preflight()`` directly for the full
    ``CREATED -> VALIDATING -> PLANNED/FAILED`` path, or schedule
    ``run_preflight()`` via ``BackgroundTasks`` after they have already
    transitioned the run to ``VALIDATING``.
    """

    def __init__(
        self,
        run_service: "RunService",
        execution_planner: "ExecutionPlanner",
        workspace_planner: "WorkspacePlanner",
        local_run_driver: "LocalRunDriver",
    ) -> None:
        self._run_service = run_service
        self._execution_planner = execution_planner
        self._workspace_planner = workspace_planner
        self._local_run_driver = local_run_driver

    def preflight(self, run_id: str) -> Result[RunRecord]:
        """Run the full local preflight path for *run_id*.

        Accepts runs only in ``CREATED`` and transitions them to
        ``VALIDATING`` before executing. Any other status is treated as a
        duplicate trigger.
        """
        try:
            record = self._run_service.get_run(run_id)
        except KeyError:
            return Result.failure(
                [
                    Issue(
                        code="PREFLIGHT_RUN_NOT_FOUND",
                        message="Run was not found.",
                        severity="error",
                        path="run_id",
                        source="preflight_service",
                    )
                ]
            )

        if record.status is not RunStatus.CREATED:
            return Result.failure(
                [
                    Issue(
                        code="PREFLIGHT_ALREADY_TRIGGERED",
                        message="Preflight has already been triggered for this run.",
                        severity="error",
                        path="run_id",
                        source="preflight_service",
                        context={"current_status": record.status.value},
                    )
                ]
            )

        self._run_service.transition_run(
            run_id,
            RunStatus.VALIDATING,
            stage="preflight",
            message="Local preflight started.",
        )
        return self._run_preflight(run_id)

    def run_preflight(self, run_id: str) -> Result[RunRecord]:
        """Background-worker entry point.

        The caller must already have transitioned the run to ``VALIDATING``.
        """
        try:
            record = self._run_service.get_run(run_id)
        except KeyError:
            return Result.failure(
                [
                    Issue(
                        code="PREFLIGHT_RUN_NOT_FOUND",
                        message="Run was not found.",
                        severity="error",
                        path="run_id",
                        source="preflight_service",
                    )
                ]
            )

        if record.status is not RunStatus.VALIDATING:
            return Result.failure(
                [
                    Issue(
                        code="PREFLIGHT_ALREADY_TRIGGERED",
                        message="Preflight has already been triggered for this run.",
                        severity="error",
                        path="run_id",
                        source="preflight_service",
                        context={"current_status": record.status.value},
                    )
                ]
            )

        return self._run_preflight(run_id)

    def _run_preflight(self, run_id: str) -> Result[RunRecord]:
        try:
            plan_result = self._execution_planner.plan_run(run_id)
            if plan_result.is_failure:
                return self._fail(run_id, plan_result.issues)

            base_plan = plan_result.value
            workspace_dir = self._local_run_driver.derive_workspace_dir(run_id)
            workspace_result = self._workspace_planner.plan_workspace(
                base_plan,
                base_dir=workspace_dir,
            )
            if workspace_result.is_failure:
                return self._fail(run_id, workspace_result.issues)

            prepared_plan = workspace_result.value
            run_result = self._local_run_driver.run(run_id, prepared_plan)
            if run_result.is_failure:
                return self._fail(run_id, run_result.issues)

            final_plan = run_result.value
            current = self._run_service.get_run(run_id)
            if current.status is RunStatus.CANCELLED:
                return self._cancelled_result()

            updated = self._run_service.transition_run(
                run_id,
                RunStatus.PLANNED,
                stage="preflight",
                message="Local preflight completed; dry-run succeeded.",
                context={"reason_code": "PREFLIGHT_COMPLETED"},
            )
            self._run_service.add_event(
                run_id=run_id,
                event_type="preflight_completed",
                message="Run is planned and ready for execution.",
                status=RunStatus.PLANNED,
                stage="preflight",
                context={
                    "plan_status": final_plan.status.value,
                    "has_command_spec": final_plan.command_spec is not None,
                    "reason_code": "PREFLIGHT_COMPLETED",
                },
            )
            return Result.success(updated)
        except Exception:
            issue = Issue(
                code="PREFLIGHT_UNEXPECTED_ERROR",
                message="An unexpected error occurred during preflight.",
                severity="error",
                source="preflight_service",
                path="preflight",
            )
            return self._fail(run_id, [issue])

    def _fail(
        self,
        run_id: str,
        issues: tuple[Issue, ...] | list[Issue],
    ) -> Result[RunRecord]:
        issue_list = list(issues)
        reason_code = issue_list[0].code if issue_list else "PREFLIGHT_FAILED"
        current = self._run_service.get_run(run_id)

        if current.status is RunStatus.CANCELLED:
            return self._cancelled_result()

        record_issue = Issue(
            code="PREFLIGHT_FAILED",
            message="Local preflight failed.",
            severity="error",
            path="preflight",
            source="preflight_service",
            context={"reason_code": reason_code},
        )
        self._run_service.transition_run(
            run_id,
            RunStatus.FAILED,
            stage="preflight",
            message="Local preflight failed.",
            issue=record_issue,
            context={"reason_code": reason_code},
        )
        return Result.failure(issue_list)

    def _cancelled_result(self) -> Result[RunRecord]:
        return Result.failure(
            [
                Issue(
                    code="PREFLIGHT_CANCELLED",
                    message="Preflight was cancelled before it completed.",
                    severity="error",
                    path="run_id",
                    source="preflight_service",
                )
            ]
        )
