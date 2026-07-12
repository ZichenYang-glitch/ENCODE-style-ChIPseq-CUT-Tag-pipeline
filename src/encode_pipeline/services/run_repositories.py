"""Storage boundary for workflow run aggregates."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
from threading import RLock
from typing import Any, Protocol, TypeVar

from encode_pipeline.platform.execution import (
    RunExecutionAssignment,
    RunExecutionCancellationRequest,
    RunExecutionClaim,
    RunExecutionOwnership,
    RunExecutionStopAcknowledgement,
)
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.results import Issue
from encode_pipeline.platform.runs import (
    RunArtifactRef,
    RunEvent,
    RunLogChunk,
    RunRecord,
    RunStatus,
)


_T = TypeVar("_T")


class ConcurrentRunUpdateError(RuntimeError):
    """Raised when a persisted run changed after it was read."""


@dataclass(frozen=True)
class RunEventDraft:
    """Event data whose repository-owned identity has not been assigned yet."""

    event_type: str
    message: str
    status: RunStatus | None = None
    stage: str | None = None
    context: Mapping[str, Any] = field(default_factory=dict)
    issue: Issue | None = None


class RunRepository(Protocol):
    """Persistence contract consumed by :class:`RunService`."""

    def contains_run(self, run_id: str) -> bool: ...

    def create_run(self, record: RunRecord, event: RunEventDraft) -> RunEvent: ...

    def get_run(self, run_id: str) -> RunRecord: ...

    def list_runs(self) -> tuple[RunRecord, ...]: ...

    def update_run(
        self,
        record: RunRecord,
        *,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent: ...

    def complete_preflight(
        self,
        record: RunRecord,
        identity: WorkflowBuildIdentity,
        *,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent: ...

    def get_workflow_build_identity(
        self,
        run_id: str,
    ) -> WorkflowBuildIdentity | None: ...

    def fail_interrupted_run_if_unowned(
        self,
        record: RunRecord,
        *,
        expected_status: RunStatus,
        required_ownership: RunExecutionOwnership | None,
        event: RunEventDraft,
    ) -> bool: ...

    def add_event(self, run_id: str, event: RunEventDraft) -> RunEvent: ...

    def list_events(
        self,
        run_id: str,
        *,
        after: str | None = None,
        limit: int = 50,
    ) -> tuple[RunEvent, ...]: ...

    def append_log(
        self,
        run_id: str,
        stream_name: str,
        lines: Iterable[str],
    ) -> RunLogChunk: ...

    def list_logs(
        self,
        run_id: str,
        stream_name: str = "stdout",
        *,
        after: str | None = None,
        limit: int = 50,
    ) -> tuple[RunLogChunk, ...]: ...

    def record_artifact(
        self,
        run_id: str,
        artifact: RunArtifactRef,
    ) -> RunArtifactRef: ...

    def replace_artifacts(
        self,
        run_id: str,
        artifacts: tuple[RunArtifactRef, ...],
        *,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent | None: ...

    def list_artifacts(self, run_id: str) -> tuple[RunArtifactRef, ...]: ...

    def ensure_execution_assignment(
        self,
        assignment: RunExecutionAssignment,
        *,
        expected_status: RunStatus,
    ) -> RunExecutionAssignment: ...

    def get_execution_assignment(
        self,
        run_id: str,
    ) -> RunExecutionAssignment | None: ...

    def mark_execution_dispatched(
        self,
        run_id: str,
        *,
        job_id: str,
        dispatched_at: datetime,
        allowed_statuses: frozenset[RunStatus],
    ) -> RunExecutionAssignment: ...

    def queue_dispatched_run(
        self,
        record: RunRecord,
        *,
        expected_status: RunStatus,
        job_id: str,
        backend: str,
        queue_name: str,
        event: RunEventDraft,
    ) -> bool: ...

    def claim_execution_assignment(
        self,
        run_id: str,
        *,
        job_id: str,
        backend: str,
        queue_name: str,
        claimed_at: datetime,
        allowed_statuses: frozenset[RunStatus],
        event: RunEventDraft,
    ) -> RunExecutionClaim: ...

    def request_execution_cancellation(
        self,
        run_id: str,
        *,
        job_id: str,
        backend: str,
        queue_name: str,
        requested_at: datetime,
        reason: str,
        event: RunEventDraft,
    ) -> RunExecutionCancellationRequest: ...

    def acknowledge_execution_stop(
        self,
        run_id: str,
        *,
        job_id: str,
        backend: str,
        queue_name: str,
        acknowledged_at: datetime,
        cancellation_event: RunEventDraft,
        unexpected_stop_event: RunEventDraft,
    ) -> RunExecutionStopAcknowledgement: ...


class InMemoryRunRepository:
    """Thread-safe in-memory implementation used by unit tests and adapters."""

    def __init__(self) -> None:
        self._lock = RLock()
        self._runs: dict[str, RunRecord] = {}
        self._events: dict[str, list[RunEvent]] = {}
        self._logs: dict[str, dict[str, list[RunLogChunk]]] = {}
        self._artifacts: dict[str, dict[str, RunArtifactRef]] = {}
        self._execution_assignments: dict[str, RunExecutionAssignment] = {}
        self._execution_run_ids_by_job: dict[str, str] = {}
        self._workflow_build_identities: dict[str, WorkflowBuildIdentity] = {}

    def contains_run(self, run_id: str) -> bool:
        with self._lock:
            return run_id in self._runs

    def create_run(self, record: RunRecord, event: RunEventDraft) -> RunEvent:
        with self._lock:
            if record.run_id in self._runs:
                raise ValueError(f"Duplicate run_id: {record.run_id!r}")
            created_event = self._make_event(record.run_id, 1, event)
            self._runs[record.run_id] = record
            self._events[record.run_id] = [created_event]
            self._logs[record.run_id] = {}
            self._artifacts[record.run_id] = {}
            return created_event

    def get_run(self, run_id: str) -> RunRecord:
        with self._lock:
            return self._runs[run_id]

    def list_runs(self) -> tuple[RunRecord, ...]:
        with self._lock:
            return tuple(self._runs.values())

    def update_run(
        self,
        record: RunRecord,
        *,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent:
        with self._lock:
            current = self._runs[record.run_id]
            if current.status is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {record.run_id!r} changed while it was being updated."
                )
            updated_event = self._make_event(
                record.run_id,
                len(self._events[record.run_id]) + 1,
                event,
            )
            self._runs[record.run_id] = record
            self._events[record.run_id].append(updated_event)
            return updated_event

    def complete_preflight(
        self,
        record: RunRecord,
        identity: WorkflowBuildIdentity,
        *,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent:
        """Atomically bind a build identity to the PLANNED transition."""
        with self._lock:
            if record.status is not RunStatus.PLANNED:
                raise ValueError("completed preflight record must be planned")
            if identity.workflow_id != record.workflow_id:
                raise ValueError("workflow build identity does not match the run")
            current = self._runs[record.run_id]
            if current.status is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {record.run_id!r} changed while preflight completed."
                )
            if record.run_id in self._workflow_build_identities:
                raise ValueError("run already has a workflow build identity")
            completed_event = self._make_event(
                record.run_id,
                len(self._events[record.run_id]) + 1,
                event,
            )
            self._runs[record.run_id] = record
            self._workflow_build_identities[record.run_id] = identity
            self._events[record.run_id].append(completed_event)
            return completed_event

    def get_workflow_build_identity(
        self,
        run_id: str,
    ) -> WorkflowBuildIdentity | None:
        with self._lock:
            return self._workflow_build_identities.get(run_id)

    def fail_interrupted_run_if_unowned(
        self,
        record: RunRecord,
        *,
        expected_status: RunStatus,
        required_ownership: RunExecutionOwnership | None,
        event: RunEventDraft,
    ) -> bool:
        with self._lock:
            current = self._runs[record.run_id]
            if current.status is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {record.run_id!r} changed while it was being recovered."
                )
            assignment = self._execution_assignments.get(record.run_id)
            if _assignment_has_ownership(assignment, required_ownership):
                return False
            updated_event = self._make_event(
                record.run_id,
                len(self._events[record.run_id]) + 1,
                event,
            )
            self._runs[record.run_id] = record
            self._events[record.run_id].append(updated_event)
            return True

    def add_event(self, run_id: str, event: RunEventDraft) -> RunEvent:
        with self._lock:
            self._runs[run_id]
            return self._append_event(run_id, event)

    def list_events(
        self,
        run_id: str,
        *,
        after: str | None = None,
        limit: int = 50,
    ) -> tuple[RunEvent, ...]:
        with self._lock:
            events = self._events[run_id]
            start = (
                0 if after is None else self._index_after(events, after, "event_id") + 1
            )
            return tuple(events[start : start + limit])

    def append_log(
        self,
        run_id: str,
        stream_name: str,
        lines: Iterable[str],
    ) -> RunLogChunk:
        with self._lock:
            stream_logs = self._logs[run_id].setdefault(stream_name, [])
            sequence = len(stream_logs) + 1
            chunk = RunLogChunk(
                chunk_id=f"log-{sequence}",
                run_id=run_id,
                stream_name=stream_name,
                sequence=sequence,
                timestamp=datetime.now(timezone.utc),
                lines=lines,
            )
            stream_logs.append(chunk)
            return chunk

    def list_logs(
        self,
        run_id: str,
        stream_name: str = "stdout",
        *,
        after: str | None = None,
        limit: int = 50,
    ) -> tuple[RunLogChunk, ...]:
        with self._lock:
            stream_logs = self._logs[run_id].get(stream_name, [])
            start = (
                0
                if after is None
                else self._index_after(stream_logs, after, "chunk_id") + 1
            )
            return tuple(stream_logs[start : start + limit])

    def record_artifact(
        self,
        run_id: str,
        artifact: RunArtifactRef,
    ) -> RunArtifactRef:
        with self._lock:
            artifacts = self._artifacts[run_id]
            if artifact.artifact_id in artifacts:
                raise ValueError(f"Duplicate artifact_id: {artifact.artifact_id!r}")
            artifacts[artifact.artifact_id] = artifact
            return artifact

    def replace_artifacts(
        self,
        run_id: str,
        artifacts: tuple[RunArtifactRef, ...],
        *,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent | None:
        with self._lock:
            current = self._runs[run_id]
            if current.status is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer eligible for artifact replacement."
                )
            replacement: dict[str, RunArtifactRef] = {}
            for artifact in artifacts:
                if artifact.run_id != run_id:
                    raise ValueError("artifact run_id does not match the run")
                if artifact.artifact_id in replacement:
                    raise ValueError("duplicate artifact_id in replacement")
                replacement[artifact.artifact_id] = artifact
            existing = tuple(self._artifacts[run_id].values())
            indexed = any(
                item.event_type == "artifacts_indexed" for item in self._events[run_id]
            )
            if existing == artifacts and indexed:
                return None
            indexed_event = self._make_event(
                run_id,
                len(self._events[run_id]) + 1,
                event,
            )
            self._artifacts[run_id] = replacement
            self._events[run_id].append(indexed_event)
            return indexed_event

    def list_artifacts(self, run_id: str) -> tuple[RunArtifactRef, ...]:
        with self._lock:
            return tuple(self._artifacts[run_id].values())

    def ensure_execution_assignment(
        self,
        assignment: RunExecutionAssignment,
        *,
        expected_status: RunStatus,
    ) -> RunExecutionAssignment:
        with self._lock:
            current = self._runs[assignment.run_id]
            if current.status is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {assignment.run_id!r} is no longer assignable."
                )
            existing = self._execution_assignments.get(assignment.run_id)
            if existing is not None:
                return existing
            assigned_run_id = self._execution_run_ids_by_job.get(assignment.job_id)
            if assigned_run_id is not None:
                raise ValueError(
                    f"Execution job_id {assignment.job_id!r} is already assigned "
                    f"to run {assigned_run_id!r}."
                )
            self._execution_assignments[assignment.run_id] = assignment
            self._execution_run_ids_by_job[assignment.job_id] = assignment.run_id
            return assignment

    def get_execution_assignment(
        self,
        run_id: str,
    ) -> RunExecutionAssignment | None:
        with self._lock:
            return self._execution_assignments.get(run_id)

    def mark_execution_dispatched(
        self,
        run_id: str,
        *,
        job_id: str,
        dispatched_at: datetime,
        allowed_statuses: frozenset[RunStatus],
    ) -> RunExecutionAssignment:
        with self._lock:
            current = self._runs[run_id]
            assignment = self._execution_assignments[run_id]
            if assignment.job_id != job_id:
                raise ValueError("job_id does not match the execution assignment")
            if assignment.dispatched_at is not None:
                return assignment
            if current.status not in allowed_statuses:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer dispatchable."
                )
            assignment = replace(assignment, dispatched_at=dispatched_at)
            self._execution_assignments[run_id] = assignment
            return assignment

    def queue_dispatched_run(
        self,
        record: RunRecord,
        *,
        expected_status: RunStatus,
        job_id: str,
        backend: str,
        queue_name: str,
        event: RunEventDraft,
    ) -> bool:
        """Queue a dispatched planned run and append exactly one event."""
        with self._lock:
            if expected_status is not RunStatus.PLANNED:
                raise ValueError("expected_status must be planned")
            if record.status is not RunStatus.QUEUED:
                raise ValueError("queued record must have queued status")

            current = self._runs[record.run_id]
            assignment = self._execution_assignments[record.run_id]
            _require_assignment_identity(
                assignment,
                job_id=job_id,
                backend=backend,
                queue_name=queue_name,
            )
            if assignment.dispatched_at is None:
                raise ValueError("execution assignment has not been dispatched")
            if current.status is not expected_status:
                return False

            queued_event = self._make_event(
                record.run_id,
                len(self._events[record.run_id]) + 1,
                event,
            )
            self._runs[record.run_id] = record
            self._events[record.run_id].append(queued_event)
            return True

    def claim_execution_assignment(
        self,
        run_id: str,
        *,
        job_id: str,
        backend: str,
        queue_name: str,
        claimed_at: datetime,
        allowed_statuses: frozenset[RunStatus],
        event: RunEventDraft,
    ) -> RunExecutionClaim:
        with self._lock:
            current = self._runs[run_id]
            assignment = self._execution_assignments[run_id]
            _require_assignment_identity(
                assignment,
                job_id=job_id,
                backend=backend,
                queue_name=queue_name,
            )
            if assignment.claimed_at is not None:
                return RunExecutionClaim(assignment=assignment, acquired=False)
            if current.status not in allowed_statuses:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer claimable."
                )
            assignment = replace(
                assignment,
                dispatched_at=assignment.dispatched_at or claimed_at,
                claimed_at=claimed_at,
            )
            self._execution_assignments[run_id] = assignment
            self._append_event(
                run_id,
                replace(event, status=current.status),
            )
            return RunExecutionClaim(assignment=assignment, acquired=True)

    def request_execution_cancellation(
        self,
        run_id: str,
        *,
        job_id: str,
        backend: str,
        queue_name: str,
        requested_at: datetime,
        reason: str,
        event: RunEventDraft,
    ) -> RunExecutionCancellationRequest:
        with self._lock:
            current = self._runs[run_id]
            assignment = self._execution_assignments[run_id]
            _require_assignment_identity(
                assignment,
                job_id=job_id,
                backend=backend,
                queue_name=queue_name,
            )
            if assignment.claimed_at is None:
                raise ValueError("execution assignment has not been claimed")
            if current.status is not RunStatus.RUNNING:
                if current.status.is_terminal:
                    return RunExecutionCancellationRequest(
                        assignment=assignment,
                        record=current,
                        created=False,
                    )
                raise ConcurrentRunUpdateError(f"Run {run_id!r} is no longer running.")
            if assignment.cancellation_requested_at is not None:
                return RunExecutionCancellationRequest(
                    assignment=assignment,
                    record=current,
                    created=False,
                )

            updated_assignment = replace(
                assignment,
                cancellation_requested_at=requested_at,
                cancellation_reason=reason,
            )
            created_event = self._make_event(
                run_id,
                len(self._events[run_id]) + 1,
                event,
            )
            self._execution_assignments[run_id] = updated_assignment
            self._events[run_id].append(created_event)
            return RunExecutionCancellationRequest(
                assignment=updated_assignment,
                record=current,
                created=True,
            )

    def acknowledge_execution_stop(
        self,
        run_id: str,
        *,
        job_id: str,
        backend: str,
        queue_name: str,
        acknowledged_at: datetime,
        cancellation_event: RunEventDraft,
        unexpected_stop_event: RunEventDraft,
    ) -> RunExecutionStopAcknowledgement:
        with self._lock:
            current = self._runs[run_id]
            assignment = self._execution_assignments[run_id]
            _require_assignment_identity(
                assignment,
                job_id=job_id,
                backend=backend,
                queue_name=queue_name,
            )
            if current.status.is_terminal:
                return RunExecutionStopAcknowledgement(
                    assignment=assignment,
                    record=current,
                    transitioned=False,
                )

            if assignment.cancellation_requested_at is not None:
                if assignment.claimed_at is None:
                    raise ValueError("execution assignment has not been claimed")
                if current.status is not RunStatus.RUNNING:
                    raise ConcurrentRunUpdateError(
                        f"Run {run_id!r} is no longer running."
                    )
                reason = assignment.cancellation_reason
                assert reason is not None
                updated_assignment = replace(
                    assignment,
                    cancellation_acknowledged_at=acknowledged_at,
                )
                updated = _terminal_record(
                    current,
                    status=RunStatus.CANCELLED,
                    ended_at=acknowledged_at,
                    cancellation_reason=reason,
                    issue=None,
                )
                event = _canonical_cancellation_event(
                    cancellation_event,
                    reason=reason,
                )
            else:
                if current.status not in {
                    RunStatus.PLANNED,
                    RunStatus.QUEUED,
                    RunStatus.RUNNING,
                }:
                    raise ConcurrentRunUpdateError(
                        f"Run {run_id!r} is not worker-owned."
                    )
                updated_assignment = assignment
                updated = _terminal_record(
                    current,
                    status=RunStatus.FAILED,
                    ended_at=acknowledged_at,
                    cancellation_reason=None,
                    issue=unexpected_stop_event.issue,
                )
                event = _canonical_unexpected_stop_event(
                    unexpected_stop_event,
                    previous_status=current.status,
                )

            created_event = self._make_event(
                run_id,
                len(self._events[run_id]) + 1,
                event,
            )
            self._execution_assignments[run_id] = updated_assignment
            self._runs[run_id] = updated
            self._events[run_id].append(created_event)
            return RunExecutionStopAcknowledgement(
                assignment=updated_assignment,
                record=updated,
                transitioned=True,
            )

    def _append_event(self, run_id: str, draft: RunEventDraft) -> RunEvent:
        events = self._events[run_id]
        sequence = len(events) + 1
        event = self._make_event(run_id, sequence, draft)
        events.append(event)
        return event

    @staticmethod
    def _make_event(
        run_id: str,
        sequence: int,
        draft: RunEventDraft,
    ) -> RunEvent:
        event = RunEvent(
            event_id=f"evt-{sequence}",
            run_id=run_id,
            sequence=sequence,
            event_type=draft.event_type,
            timestamp=datetime.now(timezone.utc),
            status=draft.status,
            stage=draft.stage,
            message=draft.message,
            context=draft.context,
            issue=draft.issue,
        )
        return event

    @staticmethod
    def _index_after(items: list[_T], cursor: str, id_attr: str) -> int:
        for index, item in enumerate(items):
            if getattr(item, id_attr) == cursor:
                return index
        raise KeyError(cursor)


def _require_assignment_identity(
    assignment: RunExecutionAssignment,
    *,
    job_id: str,
    backend: str,
    queue_name: str,
) -> None:
    if (
        assignment.job_id != job_id
        or assignment.backend != backend
        or assignment.queue_name != queue_name
    ):
        raise ValueError("execution assignment identity does not match")


def _assignment_has_ownership(
    assignment: RunExecutionAssignment | None,
    required_ownership: RunExecutionOwnership | None,
) -> bool:
    if required_ownership is None:
        return False
    if assignment is None or assignment.backend != "rq":
        return False
    if required_ownership is RunExecutionOwnership.DISPATCHED:
        return assignment.dispatched_at is not None
    return assignment.claimed_at is not None


def _terminal_record(
    current: RunRecord,
    *,
    status: RunStatus,
    ended_at: datetime,
    cancellation_reason: str | None,
    issue: Issue | None,
) -> RunRecord:
    return replace(
        current,
        status=status,
        updated_at=ended_at,
        ended_at=ended_at,
        cancellation_reason=cancellation_reason,
        error=issue,
    )


def _canonical_cancellation_event(
    event: RunEventDraft,
    *,
    reason: str,
) -> RunEventDraft:
    context = dict(event.context)
    context.update(
        {
            "previous_status": RunStatus.RUNNING.value,
            "new_status": RunStatus.CANCELLED.value,
            "cancellation_reason": reason,
        }
    )
    return replace(event, context=context)


def _canonical_unexpected_stop_event(
    event: RunEventDraft,
    *,
    previous_status: RunStatus,
) -> RunEventDraft:
    context = dict(event.context)
    context.update(
        {
            "previous_status": previous_status.value,
            "new_status": RunStatus.FAILED.value,
        }
    )
    return replace(event, context=context)
