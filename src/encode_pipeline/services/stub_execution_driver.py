"""Synchronous stub execution driver for workflow runs."""

from __future__ import annotations

from encode_pipeline.platform.results import Issue
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.services.runs import RunService


class StubExecutionDriver:
    """Synchronous stub driver that advances a run to a terminal state."""

    def __init__(self, run_service: RunService) -> None:
        if not isinstance(run_service, RunService):
            raise ValueError("StubExecutionDriver requires a RunService instance")
        self._run_service = run_service

    def advance_to_terminal(self, run_id: str) -> RunRecord:
        """Advance an active run through its lifecycle to a terminal state."""
        record = self._run_service.get_run(run_id)
        if record.status.is_terminal:
            return record

        self._transition(run_id, RunStatus.VALIDATING, stage="validation")
        self._append_log(run_id, "[stub] validating inputs")

        if record.tags.get("stub_outcome") == "failure":
            self._transition(
                run_id,
                RunStatus.FAILED,
                issue=Issue(
                    code="STUB_FAILURE",
                    message="Stub driver simulated a failure during validation.",
                    severity="error",
                    path="stub",
                    source="stub_driver",
                ),
            )
            self._append_log(run_id, "[stub] marking failed state")
            return self._run_service.get_run(run_id)

        self._transition(run_id, RunStatus.PLANNED, stage="plan")
        self._append_log(run_id, "[stub] recording plan stage")
        self._transition(run_id, RunStatus.QUEUED, stage="queue")
        self._append_log(run_id, "[stub] marking queued state")
        self._transition(run_id, RunStatus.RUNNING, stage="run")
        self._append_log(run_id, "[stub] marking running state")
        self._transition(run_id, RunStatus.SUCCEEDED)
        self._append_log(run_id, "[stub] marking completed state")
        return self._run_service.get_run(run_id)

    def _transition(
        self,
        run_id: str,
        to_status: RunStatus,
        *,
        stage: str | None = None,
        issue: Issue | None = None,
    ) -> None:
        message = self._transition_message(to_status)
        self._run_service.transition_run(
            run_id,
            to_status,
            stage=stage,
            message=message,
            issue=issue,
        )

    def _append_log(self, run_id: str, line: str) -> None:
        self._run_service.append_log(run_id, "stdout", [line])

    @staticmethod
    def _transition_message(to_status: RunStatus) -> str:
        messages = {
            RunStatus.VALIDATING: "Validation stage started.",
            RunStatus.PLANNED: "Validation stage completed.",
            RunStatus.FAILED: "Validation stage failed.",
            RunStatus.QUEUED: "Run queued.",
            RunStatus.RUNNING: "Run progress marked running.",
            RunStatus.SUCCEEDED: "Run completed.",
        }
        return messages[to_status]
