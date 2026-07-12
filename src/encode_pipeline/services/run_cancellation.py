"""Truthful cancellation orchestration for durable workflow runs."""

from __future__ import annotations

from dataclasses import dataclass

from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.services.run_queue import RunQueueError, RunStopQueue
from encode_pipeline.services.runs import (
    RunCancellationNotAvailableError,
    RunService,
)


@dataclass(frozen=True)
class RunCancellationResult:
    """Canonical run snapshot after a cancellation request."""

    record: RunRecord
    stop_requested: bool

    def __post_init__(self) -> None:
        if not isinstance(self.record, RunRecord):
            raise ValueError("record must be a RunRecord")
        if not isinstance(self.stop_requested, bool):
            raise ValueError("stop_requested must be a bool")


class RunCancellationConflictError(RuntimeError):
    """Durable execution ownership cannot be matched safely."""

    def __init__(self, record: RunRecord) -> None:
        super().__init__("Run cancellation conflicts with durable execution state.")
        self.record = record


class RunCancellationUnavailableError(RuntimeError):
    """The stop command could not be confirmed and may be retried."""

    def __init__(self, record: RunRecord) -> None:
        super().__init__("The execution queue could not confirm cancellation.")
        self.record = record


class RunCancellationService:
    """Persist intent before asking the worker backend to stop an execution."""

    def __init__(self, run_service: RunService, run_queue: RunStopQueue) -> None:
        if not isinstance(run_service, RunService):
            raise ValueError("run_service must be a RunService")
        if not isinstance(run_queue, RunStopQueue):
            raise ValueError("run_queue must implement RunStopQueue")
        self._run_service = run_service
        self._run_queue = run_queue

    def cancel_run(
        self,
        run_id: str,
        *,
        reason: str,
    ) -> RunCancellationResult:
        """Cancel before execution, or request a truthful RQ stop while running."""
        if not isinstance(reason, str) or not reason.strip():
            raise ValueError("cancellation reason must be a non-empty string")

        current = self._run_service.get_run(run_id)
        if current.status.is_terminal:
            return RunCancellationResult(record=current, stop_requested=False)
        if current.status is not RunStatus.RUNNING:
            try:
                cancelled = self._run_service.cancel_run(run_id, reason=reason.strip())
            except RunCancellationNotAvailableError as exc:
                current = exc.record
            else:
                return RunCancellationResult(
                    record=cancelled,
                    stop_requested=False,
                )

        assignment = self._run_service.get_execution_assignment(run_id)
        if (
            assignment is None
            or assignment.run_id != run_id
            or assignment.backend != self._run_queue.backend
            or assignment.queue_name != self._run_queue.queue_name
            or assignment.claimed_at is None
        ):
            raise RunCancellationConflictError(self._run_service.get_run(run_id))

        try:
            requested = self._run_service.request_execution_cancellation(
                run_id,
                job_id=assignment.job_id,
                backend=assignment.backend,
                queue_name=assignment.queue_name,
                reason=reason.strip(),
            )
        except (KeyError, ValueError) as exc:
            raise RunCancellationConflictError(
                self._run_service.get_run(run_id)
            ) from exc

        if requested.record.status.is_terminal:
            return RunCancellationResult(
                record=requested.record,
                stop_requested=False,
            )

        try:
            self._run_queue.request_stop(requested.assignment)
        except RunQueueError as exc:
            canonical = self._run_service.get_run(run_id)
            if canonical.status.is_terminal:
                return RunCancellationResult(
                    record=canonical,
                    stop_requested=False,
                )
            raise RunCancellationUnavailableError(canonical) from exc

        canonical = self._run_service.get_run(run_id)
        return RunCancellationResult(
            record=canonical,
            stop_requested=not canonical.status.is_terminal,
        )
