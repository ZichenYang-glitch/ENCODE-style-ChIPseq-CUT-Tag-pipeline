"""Durable submission orchestration for planned workflow runs."""

from __future__ import annotations

from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.services.run_queue import (
    RunQueue,
    RunQueueError,
    RunQueueIdentityError,
    RunQueueJobUnavailableError,
)
from encode_pipeline.services.run_repositories import ConcurrentRunUpdateError
from encode_pipeline.services.runs import RunService


class RunSubmissionError(RuntimeError):
    """Base error for a rejected or unconfirmed run submission."""

    def __init__(self, message: str, *, record: RunRecord) -> None:
        super().__init__(message)
        self.record = record


class RunNotReadyError(RunSubmissionError):
    """The run has not completed preflight."""


class RunBuildIdentityMissingError(RunSubmissionError):
    """A legacy planned run has no durable workflow build identity."""


class RunStartConflictError(RunSubmissionError):
    """Durable run and queue identities cannot be reconciled."""


class RunSubmissionUnavailableError(RunSubmissionError):
    """The queue could not confirm submission and the run remains retryable."""


class RunSubmissionService:
    """Reserve, enqueue, and durably queue one workflow execution."""

    def __init__(self, run_service: RunService, run_queue: RunQueue) -> None:
        if not isinstance(run_service, RunService):
            raise ValueError("run_service must be a RunService instance")
        if not isinstance(run_queue, RunQueue):
            raise ValueError("run_queue must implement RunQueue")
        self._run_service = run_service
        self._run_queue = run_queue

    def start_run(self, run_id: str) -> RunRecord:
        """Submit a planned run, converging retries on one durable job ID.

        A successful return means the backend accepted the stable job identity
        and the database contains both a dispatch marker and ``QUEUED`` (or a
        later state). Queue failures never synthesize ``QUEUED`` state, so a
        durable reservation remains safe to retry unless a worker independently
        advanced the database.
        """
        current = self._run_service.get_run(run_id)
        if current.status in {RunStatus.CREATED, RunStatus.VALIDATING}:
            raise RunNotReadyError(
                "Run must complete preflight before it can start.",
                record=current,
            )
        if (
            current.status is RunStatus.PLANNED
            and self._run_service.get_workflow_build_identity(run_id) is None
        ):
            raise RunBuildIdentityMissingError(
                "Run must be preflighted with a durable workflow build identity.",
                record=current,
            )

        assignment = self._resolve_assignment(current)
        current = self._run_service.get_run(run_id)
        if current.status in {
            RunStatus.RUNNING,
            RunStatus.SUCCEEDED,
            RunStatus.FAILED,
            RunStatus.CANCELLED,
        }:
            latest_assignment = self._run_service.get_execution_assignment(run_id)
            if latest_assignment is not None and latest_assignment.dispatched_at:
                return current
            raise RunStartConflictError(
                "Run reached a later state without a confirmed dispatch.",
                record=current,
            )
        if current.status not in {RunStatus.PLANNED, RunStatus.QUEUED}:
            raise RunStartConflictError(
                "Run changed before execution submission could begin.",
                record=current,
            )

        try:
            backend_job_id = self._run_queue.enqueue_execution(assignment)
        except (RunQueueIdentityError, RunQueueJobUnavailableError) as exc:
            latest, confirmed = self._latest_submission_state(
                run_id,
                current,
                expected_assignment=assignment,
            )
            if confirmed is not None:
                return confirmed
            raise RunStartConflictError(
                "The durable execution job cannot be reused.",
                record=latest,
            ) from exc
        except RunQueueError as exc:
            latest, confirmed = self._latest_submission_state(
                run_id,
                current,
                expected_assignment=assignment,
            )
            if confirmed is not None:
                return confirmed
            raise RunSubmissionUnavailableError(
                "The execution queue could not confirm submission.",
                record=latest,
            ) from exc

        if backend_job_id != assignment.job_id:
            raise RunStartConflictError(
                "The queue returned an unexpected execution job identity.",
                record=self._latest(run_id, current),
            )

        try:
            self._run_service.mark_execution_dispatched(
                run_id,
                job_id=assignment.job_id,
            )
            return self._run_service.queue_dispatched_run(
                run_id,
                job_id=assignment.job_id,
                backend=self._run_queue.backend,
                queue_name=self._run_queue.queue_name,
            )
        except ConcurrentRunUpdateError:
            raced = self._run_service.get_run(run_id)
            raced_assignment = self._run_service.get_execution_assignment(run_id)
            if (
                raced.status
                in {
                    RunStatus.QUEUED,
                    RunStatus.RUNNING,
                    RunStatus.SUCCEEDED,
                    RunStatus.FAILED,
                    RunStatus.CANCELLED,
                }
                and raced_assignment is not None
                and raced_assignment.dispatched_at is not None
                and self._assignment_matches_queue(raced_assignment)
            ):
                return raced
            raise RunStartConflictError(
                "Run changed while execution submission was being persisted.",
                record=raced,
            ) from None
        except ValueError as exc:
            raise RunStartConflictError(
                "The durable execution assignment is inconsistent.",
                record=self._latest(run_id, current),
            ) from exc

    def _resolve_assignment(self, current: RunRecord) -> RunExecutionAssignment:
        assignment = self._run_service.get_execution_assignment(current.run_id)
        if current.status is RunStatus.PLANNED:
            try:
                return self._run_service.ensure_execution_assignment(
                    current.run_id,
                    backend=self._run_queue.backend,
                    queue_name=self._run_queue.queue_name,
                )
            except (ConcurrentRunUpdateError, ValueError) as exc:
                raced = self._run_service.get_run(current.run_id)
                assignment = self._run_service.get_execution_assignment(current.run_id)
                if (
                    raced.status
                    in {
                        RunStatus.QUEUED,
                        RunStatus.RUNNING,
                        RunStatus.SUCCEEDED,
                        RunStatus.FAILED,
                        RunStatus.CANCELLED,
                    }
                    and assignment is not None
                    and assignment.dispatched_at is not None
                    and self._assignment_matches_queue(assignment)
                ):
                    return assignment
                raise RunStartConflictError(
                    "Run changed while its execution assignment was reserved.",
                    record=raced,
                ) from exc

        if assignment is None or assignment.dispatched_at is None:
            raise RunStartConflictError(
                "Run has no confirmed durable execution dispatch.",
                record=current,
            )
        if not self._assignment_matches_queue(assignment):
            raise RunStartConflictError(
                "The durable execution assignment does not match the queue.",
                record=current,
            )
        return assignment

    def _latest_submission_state(
        self,
        run_id: str,
        fallback: RunRecord,
        *,
        expected_assignment: RunExecutionAssignment,
    ) -> tuple[RunRecord, RunRecord | None]:
        """Return one latest error snapshot and any confirmed worker progress."""
        assignment = self._run_service.get_execution_assignment(run_id)
        current = self._latest(run_id, fallback)
        if current.status in {
            RunStatus.RUNNING,
            RunStatus.SUCCEEDED,
            RunStatus.FAILED,
            RunStatus.CANCELLED,
        }:
            # The run may have advanced after the first assignment read. Re-read
            # the monotonic ownership marker before treating the queue exception
            # as an idempotent success.
            assignment = self._run_service.get_execution_assignment(run_id)
        if (
            current.status
            in {
                RunStatus.RUNNING,
                RunStatus.SUCCEEDED,
                RunStatus.FAILED,
                RunStatus.CANCELLED,
            }
            and assignment is not None
            and assignment.run_id == expected_assignment.run_id
            and assignment.job_id == expected_assignment.job_id
            and assignment.dispatched_at is not None
            and self._assignment_matches_queue(assignment)
        ):
            return current, current
        return current, None

    def _assignment_matches_queue(
        self,
        assignment: RunExecutionAssignment,
    ) -> bool:
        return (
            assignment.backend == self._run_queue.backend
            and assignment.queue_name == self._run_queue.queue_name
        )

    def _latest(self, run_id: str, fallback: RunRecord) -> RunRecord:
        try:
            return self._run_service.get_run(run_id)
        except KeyError:
            return fallback
