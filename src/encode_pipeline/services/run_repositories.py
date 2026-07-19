"""Storage boundary for workflow run aggregates."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field, replace
from datetime import datetime, timezone
from decimal import Decimal, InvalidOperation
import re
from threading import RLock
from typing import Any, Protocol, TypeVar

from encode_pipeline.platform.execution import (
    RunExecutionAssignment,
    RunExecutionCancellationRequest,
    RunExecutionClaim,
    RunExecutionOwnership,
    RunExecutionStopAcknowledgement,
)
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.results import Issue
from encode_pipeline.platform.result_generations import (
    RunResultState,
    artifact_manifest_digest,
    build_artifact_generation,
    build_qc_generation,
    qc_metric_manifest_digest,
    validate_artifact_generation,
    validate_qc_generation,
    validate_result_attempt_id,
)
from encode_pipeline.platform.run_history import RunHistoryCursor, RunSummary
from encode_pipeline.platform.runs import (
    RunArtifactRef,
    RunEvent,
    RunLogChunk,
    RunQcMetric,
    RunRecord,
    RunStatus,
    build_qc_metric_id,
    validate_qc_identifier_token,
)
from encode_pipeline.platform.snapshots import (
    ValidatedInputSnapshot,
    ValidatedSnapshotRunCreation,
    canonical_workflow_inputs_json,
)


_T = TypeVar("_T")
_CANONICAL_DECIMAL = re.compile(r"^-?(?:0|[1-9]\d{0,25})(?:\.\d{1,12})?$")
_QC_METRIC_ID = re.compile(r"^qcmetric-[0-9a-f]{64}$")
_QC_METRIC_KEY = re.compile(r"^[a-z][a-z0-9_]*(?:\.[a-z][a-z0-9_]*)*$")
_QC_SOURCE_ARTIFACT_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,127}$")
_QC_UNITS = frozenset({"count", "fraction", "ratio", "score"})
_QC_SCOPES = frozenset({"run", "sample", "experiment"})
_QC_FLAGS = frozenset({"pass", "warning", "fail"})
_QC_OUTCOME_TYPES = frozenset(
    {
        "qc_metrics_indexed",
        "qc_metrics_indexing_failed",
        "qc_metrics_invalidated",
    }
)


class ConcurrentRunUpdateError(RuntimeError):
    """Raised when a persisted run changed after it was read."""


class ResultGenerationChangedError(ConcurrentRunUpdateError):
    """Raised when a generation-bound result read or write became stale."""


class ValidatedSnapshotExpiredError(RuntimeError):
    """Raised when an unconsumed validation snapshot passed its first-use TTL."""


class ValidatedSnapshotBuildMismatchError(RuntimeError):
    """Raised when current workflow source differs from validated source."""


class ValidatedSnapshotReplayConflictError(RuntimeError):
    """Raised when a consumed snapshot is replayed with different metadata."""


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

    def create_validated_input_snapshot(
        self,
        snapshot: ValidatedInputSnapshot,
    ) -> ValidatedInputSnapshot: ...

    def get_validated_input_snapshot(
        self,
        snapshot_id: str,
    ) -> ValidatedInputSnapshot: ...

    def consume_validated_input_snapshot(
        self,
        snapshot_id: str,
        *,
        workflow_id: str,
        expected_build_identity: WorkflowBuildIdentity,
        record: RunRecord,
        consumed_at: datetime,
        event: RunEventDraft,
    ) -> ValidatedSnapshotRunCreation: ...

    def get_run(self, run_id: str) -> RunRecord: ...

    def list_runs(self) -> tuple[RunRecord, ...]: ...

    def list_run_summaries(
        self,
        *,
        after: RunHistoryCursor | None = None,
        limit: int = 50,
        workflow_id: str | None = None,
        status: RunStatus | None = None,
    ) -> tuple[RunSummary, ...]: ...

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

    def get_result_state(self, run_id: str) -> RunResultState: ...

    def begin_artifact_result_attempt(
        self,
        run_id: str,
        *,
        attempt_id: str,
        expected_status: RunStatus,
    ) -> RunResultState: ...

    def begin_qc_result_attempt(
        self,
        run_id: str,
        *,
        attempt_id: str,
        expected_artifact_generation: str,
        expected_artifacts: tuple[RunArtifactRef, ...],
        expected_status: RunStatus,
    ) -> RunResultState: ...

    def replace_artifacts(
        self,
        run_id: str,
        artifacts: tuple[RunArtifactRef, ...],
        *,
        attempt_id: str,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent | None: ...

    def record_artifact_failure(
        self,
        run_id: str,
        *,
        attempt_id: str,
        reason_code: str,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent | None: ...

    def list_artifacts(
        self,
        run_id: str,
        *,
        after: str | None = None,
        limit: int | None = None,
    ) -> tuple[RunArtifactRef, ...]: ...

    def list_artifacts_page(
        self,
        run_id: str,
        *,
        expected_generation: str | None,
        after: str | None = None,
        limit: int | None = None,
    ) -> tuple[str | None, tuple[RunArtifactRef, ...]]: ...

    def get_artifact(self, run_id: str, artifact_id: str) -> RunArtifactRef: ...

    def get_artifact_at_generation(
        self,
        run_id: str,
        artifact_id: str,
        *,
        expected_generation: str | None,
    ) -> tuple[str, RunArtifactRef]: ...

    def replace_qc_metrics(
        self,
        run_id: str,
        metrics: tuple[RunQcMetric, ...],
        *,
        attempt_id: str,
        expected_artifact_generation: str,
        expected_artifacts: tuple[RunArtifactRef, ...],
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent | None: ...

    def list_qc_metrics(
        self,
        run_id: str,
        *,
        after: str | None = None,
        limit: int | None = None,
    ) -> tuple[RunQcMetric, ...]: ...

    def list_qc_metrics_page(
        self,
        run_id: str,
        *,
        expected_generation: str | None,
        after: str | None = None,
        limit: int | None = None,
    ) -> tuple[str | None, tuple[RunQcMetric, ...]]: ...

    def record_qc_metrics_failure(
        self,
        run_id: str,
        *,
        attempt_id: str,
        expected_artifact_generation: str,
        reason_code: str,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent | None: ...

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
        self._qc_metrics: dict[str, dict[str, RunQcMetric]] = {}
        self._result_states: dict[str, RunResultState] = {}
        self._result_attempts: dict[str, tuple[str, str, str | None]] = {}
        self._execution_assignments: dict[str, RunExecutionAssignment] = {}
        self._execution_run_ids_by_job: dict[str, str] = {}
        self._workflow_build_identities: dict[str, WorkflowBuildIdentity] = {}
        self._validated_input_snapshots: dict[str, ValidatedInputSnapshot] = {}

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
            self._qc_metrics[record.run_id] = {}
            self._result_states[record.run_id] = RunResultState(run_id=record.run_id)
            return created_event

    def create_validated_input_snapshot(
        self,
        snapshot: ValidatedInputSnapshot,
    ) -> ValidatedInputSnapshot:
        with self._lock:
            validated = _validated_snapshot_copy(snapshot)
            if validated.snapshot_id in self._validated_input_snapshots:
                raise ValueError(
                    f"Duplicate validated snapshot ID: {validated.snapshot_id!r}"
                )
            self._validated_input_snapshots[validated.snapshot_id] = validated
            return _validated_snapshot_copy(validated)

    def get_validated_input_snapshot(
        self,
        snapshot_id: str,
    ) -> ValidatedInputSnapshot:
        with self._lock:
            return _validated_snapshot_copy(
                self._validated_input_snapshots[snapshot_id]
            )

    def consume_validated_input_snapshot(
        self,
        snapshot_id: str,
        *,
        workflow_id: str,
        expected_build_identity: WorkflowBuildIdentity,
        record: RunRecord,
        consumed_at: datetime,
        event: RunEventDraft,
    ) -> ValidatedSnapshotRunCreation:
        with self._lock:
            snapshot = _validated_snapshot_copy(
                self._validated_input_snapshots[snapshot_id]
            )
            if snapshot.workflow_id != workflow_id:
                raise KeyError(snapshot_id)
            _validate_snapshot_consumption_candidate(
                snapshot,
                workflow_id=workflow_id,
                expected_build_identity=expected_build_identity,
                record=record,
                consumed_at=consumed_at,
            )
            if snapshot.consumed_run_id is not None:
                current = self._runs.get(snapshot.consumed_run_id)
                if current is None:
                    raise ValueError("validated snapshot references a missing run")
                _validate_snapshot_linked_run(snapshot, current)
                if dict(current.tags) != dict(record.tags):
                    raise ValidatedSnapshotReplayConflictError(
                        "validated snapshot replay metadata differs"
                    )
                return ValidatedSnapshotRunCreation(record=current, created=False)
            if consumed_at >= snapshot.expires_at:
                raise ValidatedSnapshotExpiredError(
                    "validated snapshot expired before first use"
                )
            if record.run_id in self._runs:
                raise ValueError(f"Duplicate run_id: {record.run_id!r}")

            created_event = self._make_event(record.run_id, 1, event)
            consumed = snapshot.with_consumption(record.run_id, consumed_at)
            self._runs[record.run_id] = record
            self._events[record.run_id] = [created_event]
            self._logs[record.run_id] = {}
            self._artifacts[record.run_id] = {}
            self._qc_metrics[record.run_id] = {}
            self._result_states[record.run_id] = RunResultState(run_id=record.run_id)
            self._validated_input_snapshots[snapshot_id] = consumed
            return ValidatedSnapshotRunCreation(record=record, created=True)

    def get_run(self, run_id: str) -> RunRecord:
        with self._lock:
            return self._runs[run_id]

    def list_runs(self) -> tuple[RunRecord, ...]:
        with self._lock:
            return tuple(self._runs.values())

    def list_run_summaries(
        self,
        *,
        after: RunHistoryCursor | None = None,
        limit: int = 50,
        workflow_id: str | None = None,
        status: RunStatus | None = None,
    ) -> tuple[RunSummary, ...]:
        """Return one validated, descending keyset page of public summaries."""
        if isinstance(limit, bool) or not isinstance(limit, int) or limit <= 0:
            raise ValueError("run summary limit must be a positive integer")
        if workflow_id is not None and not isinstance(workflow_id, str):
            raise ValueError("run summary workflow_id filter must be a string")
        if status is not None and not isinstance(status, RunStatus):
            raise ValueError("run summary status filter must be a RunStatus")
        if after is not None and not isinstance(after, RunHistoryCursor):
            raise ValueError("run summary cursor is invalid")
        if after is not None and (
            after.workflow_id != workflow_id or after.status is not status
        ):
            raise ValueError("run summary cursor filters do not match")

        with self._lock:
            if after is not None:
                boundary_record = self._runs.get(after.run_id)
                if boundary_record is None:
                    raise KeyError(after.run_id)
                boundary = _summary_from_record(after.run_id, boundary_record)
                if boundary.created_at != after.created_at or not _summary_matches(
                    boundary,
                    workflow_id=workflow_id,
                    status=status,
                ):
                    raise KeyError(after.run_id)

            selected: list[RunSummary] = []
            for mapping_key, record in self._runs.items():
                if workflow_id is not None and record.workflow_id != workflow_id:
                    continue
                if status is not None and record.status is not status:
                    continue
                summary = _summary_from_record(mapping_key, record)
                if after is not None and not _summary_is_after(summary, after):
                    continue
                selected.append(summary)
            selected.sort(
                key=lambda summary: (summary.created_at, summary.run_id),
                reverse=True,
            )
            return tuple(selected[:limit])

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

    def get_result_state(self, run_id: str) -> RunResultState:
        with self._lock:
            self._runs[run_id]
            return self._result_states[run_id]

    def begin_artifact_result_attempt(
        self,
        run_id: str,
        *,
        attempt_id: str,
        expected_status: RunStatus,
    ) -> RunResultState:
        validate_result_attempt_id(attempt_id)
        with self._lock:
            if self._runs[run_id].status is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer eligible for artifact indexing."
                )
            state = self._result_states[run_id]
            previous = self._result_attempts.get(attempt_id)
            if previous is not None:
                if previous == (run_id, "artifact", None) and (
                    state.artifact_attempt_id == attempt_id
                ):
                    return state
                raise ConcurrentRunUpdateError(
                    "artifact result attempt was already superseded"
                )
            self._result_attempts[attempt_id] = (run_id, "artifact", None)
            state = replace(
                state,
                artifact_attempt_id=attempt_id,
                artifact_attempt_status="pending",
            )
            self._result_states[run_id] = state
            return state

    def begin_qc_result_attempt(
        self,
        run_id: str,
        *,
        attempt_id: str,
        expected_artifact_generation: str,
        expected_artifacts: tuple[RunArtifactRef, ...],
        expected_status: RunStatus,
    ) -> RunResultState:
        validate_result_attempt_id(attempt_id)
        validate_artifact_generation(expected_artifact_generation)
        with self._lock:
            if self._runs[run_id].status is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer eligible for QC indexing."
                )
            state = self._result_states[run_id]
            if state.artifact_generation != expected_artifact_generation:
                raise ResultGenerationChangedError(
                    f"Run {run_id!r} artifact generation changed before QC indexing."
                )
            current_artifacts = _sorted_artifacts(self._artifacts[run_id].values())
            if current_artifacts != _validated_expected_artifacts(
                run_id,
                expected_artifacts,
            ):
                raise ResultGenerationChangedError(
                    f"Run {run_id!r} artifact manifest changed before QC indexing."
                )
            previous = self._result_attempts.get(attempt_id)
            if previous is not None:
                if previous == (run_id, "qc", expected_artifact_generation) and (
                    state.qc_attempt_id == attempt_id
                    and state.qc_attempt_artifact_generation
                    == expected_artifact_generation
                ):
                    return state
                raise ConcurrentRunUpdateError(
                    "QC result attempt was already superseded"
                )
            self._result_attempts[attempt_id] = (
                run_id,
                "qc",
                expected_artifact_generation,
            )
            state = replace(
                state,
                qc_attempt_id=attempt_id,
                qc_attempt_status="pending",
                qc_attempt_artifact_generation=expected_artifact_generation,
            )
            self._result_states[run_id] = state
            return state

    def replace_artifacts(
        self,
        run_id: str,
        artifacts: tuple[RunArtifactRef, ...],
        *,
        attempt_id: str,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent | None:
        validate_result_attempt_id(attempt_id)
        with self._lock:
            current = self._runs[run_id]
            if current.status is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer eligible for artifact replacement."
                )
            state = self._result_states[run_id]
            _require_current_attempt(
                state.artifact_attempt_id,
                attempt_id,
                "artifact",
            )
            sorted_replacement = _validated_expected_artifacts(run_id, artifacts)
            manifest_digest = artifact_manifest_digest(sorted_replacement)
            if state.artifact_attempt_status == "succeeded":
                if state.artifact_manifest_digest == manifest_digest:
                    return None
                raise ConcurrentRunUpdateError("artifact attempt already committed")
            if state.artifact_attempt_status != "pending":
                raise ConcurrentRunUpdateError("artifact attempt is no longer current")
            changed = state.artifact_manifest_digest != manifest_digest
            replacement = {
                artifact.artifact_id: artifact for artifact in sorted_replacement
            }
            should_emit = changed or state.artifact_outcome != "succeeded"
            prepared_events: list[RunEvent] = []
            if changed:
                revision = state.artifact_revision + 1
                generation = build_artifact_generation(
                    run_id=run_id,
                    revision=revision,
                    artifacts=sorted_replacement,
                )
                had_qc_state = bool(self._qc_metrics[run_id]) or (
                    state.qc_generation is not None or state.qc_outcome is not None
                )
                state = _state_after_artifact_change(
                    state,
                    artifacts=sorted_replacement,
                    artifact_revision=revision,
                    artifact_generation=generation,
                    attempt_id=attempt_id,
                    attempt_status="succeeded",
                    outcome="succeeded",
                )
                if had_qc_state:
                    prepared_events.append(
                        self._make_event(
                            run_id,
                            len(self._events[run_id]) + len(prepared_events) + 1,
                            _event_with_context(
                                RunEventDraft(
                                    event_type="qc_metrics_invalidated",
                                    message="Workflow QC metrics invalidated.",
                                    status=RunStatus.SUCCEEDED,
                                    stage="qc_summary_indexing",
                                ),
                                artifact_generation=generation,
                                qc_generation=state.qc_generation,
                            ),
                        )
                    )
            else:
                state = replace(
                    state,
                    artifact_attempt_status="succeeded",
                    artifact_outcome="succeeded",
                    artifact_reason_code=None,
                )
            result_event: RunEvent | None = None
            if should_emit:
                result_event = self._make_event(
                    run_id,
                    len(self._events[run_id]) + len(prepared_events) + 1,
                    _event_with_context(
                        event,
                        attempt_id=attempt_id,
                        artifact_generation=state.artifact_generation,
                        artifact_count=len(sorted_replacement),
                    ),
                )
                prepared_events.append(result_event)

            if changed:
                self._artifacts[run_id] = replacement
                self._qc_metrics[run_id] = {}
            self._result_states[run_id] = state
            self._events[run_id].extend(prepared_events)
            return result_event

    def record_artifact_failure(
        self,
        run_id: str,
        *,
        attempt_id: str,
        reason_code: str,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent | None:
        validate_result_attempt_id(attempt_id)
        with self._lock:
            if self._runs[run_id].status is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer eligible for artifact failure."
                )
            state = self._result_states[run_id]
            _require_current_attempt(state.artifact_attempt_id, attempt_id, "artifact")
            if state.artifact_attempt_status == "succeeded":
                return None
            if state.artifact_attempt_status == "failed":
                if state.artifact_reason_code == reason_code:
                    return None
                raise ConcurrentRunUpdateError("artifact attempt already failed")
            next_state = replace(
                state,
                artifact_attempt_status="failed",
                artifact_outcome="failed",
                artifact_reason_code=reason_code,
            )
            result_event = self._make_event(
                run_id,
                len(self._events[run_id]) + 1,
                _event_with_context(
                    event, attempt_id=attempt_id, reason_code=reason_code
                ),
            )
            self._result_states[run_id] = next_state
            self._events[run_id].append(result_event)
            return result_event

    def list_artifacts(
        self,
        run_id: str,
        *,
        after: str | None = None,
        limit: int | None = None,
    ) -> tuple[RunArtifactRef, ...]:
        _, artifacts = self.list_artifacts_page(
            run_id,
            expected_generation=None,
            after=after,
            limit=limit,
        )
        return artifacts

    def list_artifacts_page(
        self,
        run_id: str,
        *,
        expected_generation: str | None,
        after: str | None = None,
        limit: int | None = None,
    ) -> tuple[str | None, tuple[RunArtifactRef, ...]]:
        if expected_generation is not None:
            validate_artifact_generation(expected_generation)
        with self._lock:
            self._runs[run_id]
            state = self._result_states[run_id]
            if (
                expected_generation is not None
                and state.artifact_generation != expected_generation
            ):
                raise ResultGenerationChangedError("artifact generation changed")
            if self._artifacts[run_id] and state.artifact_generation is None:
                raise ValueError("artifact generation is unbound")
            artifacts = tuple(
                sorted(
                    self._artifacts[run_id].values(),
                    key=lambda artifact: artifact.artifact_id,
                )
            )
            start = (
                0
                if after is None
                else self._index_after(artifacts, after, "artifact_id") + 1
            )
            if limit is not None and limit <= 0:
                raise ValueError("limit must be positive")
            return (
                state.artifact_generation,
                artifacts[start:]
                if limit is None
                else artifacts[start : start + limit],
            )

    def get_artifact(self, run_id: str, artifact_id: str) -> RunArtifactRef:
        _, artifact = self.get_artifact_at_generation(
            run_id,
            artifact_id,
            expected_generation=None,
        )
        return artifact

    def get_artifact_at_generation(
        self,
        run_id: str,
        artifact_id: str,
        *,
        expected_generation: str | None,
    ) -> tuple[str, RunArtifactRef]:
        if expected_generation is not None:
            validate_artifact_generation(expected_generation)
        with self._lock:
            self._runs[run_id]
            state = self._result_states[run_id]
            if (
                expected_generation is not None
                and state.artifact_generation != expected_generation
            ):
                raise ResultGenerationChangedError("artifact generation changed")
            try:
                artifact = self._artifacts[run_id][artifact_id]
            except KeyError:
                raise KeyError((run_id, artifact_id)) from None
            if state.artifact_generation is None:
                raise ValueError("artifact generation is unbound")
            return state.artifact_generation, artifact

    def replace_qc_metrics(
        self,
        run_id: str,
        metrics: tuple[RunQcMetric, ...],
        *,
        attempt_id: str,
        expected_artifact_generation: str,
        expected_artifacts: tuple[RunArtifactRef, ...],
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent | None:
        validate_result_attempt_id(attempt_id)
        validate_artifact_generation(expected_artifact_generation)
        with self._lock:
            current = self._runs[run_id]
            if current.status is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer eligible for QC replacement."
                )
            state = self._result_states[run_id]
            _require_current_attempt(state.qc_attempt_id, attempt_id, "QC")
            if (
                state.artifact_generation != expected_artifact_generation
                or state.qc_attempt_artifact_generation != expected_artifact_generation
            ):
                raise ResultGenerationChangedError(
                    f"Run {run_id!r} artifact generation changed during QC indexing."
                )
            current_artifacts = _sorted_artifacts(self._artifacts[run_id].values())
            if current_artifacts != _validated_expected_artifacts(
                run_id,
                expected_artifacts,
            ):
                raise ResultGenerationChangedError(
                    f"Run {run_id!r} artifact generation changed during QC indexing."
                )
            replacement = _validated_qc_replacement(
                run_id,
                metrics,
                {artifact.artifact_id for artifact in current_artifacts},
            )
            manifest_digest = qc_metric_manifest_digest(replacement)
            if state.qc_attempt_status == "succeeded":
                if (
                    state.qc_manifest_digest == manifest_digest
                    and state.qc_artifact_generation == expected_artifact_generation
                ):
                    return None
                raise ConcurrentRunUpdateError("QC attempt already committed")
            if state.qc_attempt_status != "pending":
                raise ConcurrentRunUpdateError("QC attempt is no longer current")
            changed = (
                state.qc_manifest_digest != manifest_digest
                or state.qc_artifact_generation != expected_artifact_generation
                or state.qc_outcome != "succeeded"
            )
            if changed:
                revision = state.qc_revision + 1
                generation = build_qc_generation(
                    run_id=run_id,
                    revision=revision,
                    artifact_generation=expected_artifact_generation,
                    metrics=replacement,
                )
                replacement_by_id = {metric.metric_id: metric for metric in replacement}
                state = replace(
                    state,
                    qc_revision=revision,
                    qc_generation=generation,
                    qc_manifest_digest=manifest_digest,
                    qc_attempt_status="succeeded",
                    qc_artifact_generation=expected_artifact_generation,
                    qc_outcome="succeeded",
                    qc_reason_code=None,
                )
            else:
                state = replace(state, qc_attempt_status="succeeded")
                self._result_states[run_id] = state
                return None
            result_event = self._make_event(
                run_id,
                len(self._events[run_id]) + 1,
                _event_with_context(
                    event,
                    attempt_id=attempt_id,
                    artifact_generation=expected_artifact_generation,
                    qc_generation=state.qc_generation,
                    metric_count=len(replacement),
                ),
            )
            self._qc_metrics[run_id] = replacement_by_id
            self._result_states[run_id] = state
            self._events[run_id].append(result_event)
            return result_event

    def list_qc_metrics(
        self,
        run_id: str,
        *,
        after: str | None = None,
        limit: int | None = None,
    ) -> tuple[RunQcMetric, ...]:
        _, metrics = self.list_qc_metrics_page(
            run_id,
            expected_generation=None,
            after=after,
            limit=limit,
        )
        return metrics

    def list_qc_metrics_page(
        self,
        run_id: str,
        *,
        expected_generation: str | None,
        after: str | None = None,
        limit: int | None = None,
    ) -> tuple[str | None, tuple[RunQcMetric, ...]]:
        if expected_generation is not None:
            validate_qc_generation(expected_generation)
        with self._lock:
            self._runs[run_id]
            state = self._result_states[run_id]
            if (
                expected_generation is not None
                and state.qc_generation != expected_generation
            ):
                raise ResultGenerationChangedError("QC generation changed")
            if self._qc_metrics[run_id] and (
                state.qc_generation is None
                or state.qc_artifact_generation != state.artifact_generation
            ):
                raise ValueError("QC generation is unbound")
            by_id = self._qc_metrics[run_id]
            ordered_ids = tuple(sorted(by_id))
            if after is not None:
                try:
                    cursor = by_id[after]
                except KeyError:
                    raise KeyError((run_id, after)) from None
                _validate_stored_qc_metric_identity(
                    cursor,
                    run_id=run_id,
                    storage_key=after,
                )
                start = ordered_ids.index(after) + 1
            else:
                start = 0
            selected_ids = (
                ordered_ids[start:]
                if limit is None
                else ordered_ids[start : start + limit]
            )
            selected: list[RunQcMetric] = []
            for storage_key in selected_ids:
                metric = by_id[storage_key]
                _validate_stored_qc_metric_identity(
                    metric,
                    run_id=run_id,
                    storage_key=storage_key,
                )
                selected.append(metric)
            return state.qc_generation, tuple(selected)

    def record_qc_metrics_failure(
        self,
        run_id: str,
        *,
        attempt_id: str,
        expected_artifact_generation: str,
        reason_code: str,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent | None:
        validate_result_attempt_id(attempt_id)
        validate_artifact_generation(expected_artifact_generation)
        with self._lock:
            if self._runs[run_id].status is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer eligible for QC failure."
                )
            state = self._result_states[run_id]
            _require_current_attempt(state.qc_attempt_id, attempt_id, "QC")
            if (
                state.artifact_generation != expected_artifact_generation
                or state.qc_attempt_artifact_generation != expected_artifact_generation
            ):
                raise ResultGenerationChangedError(
                    "QC attempt artifact generation changed"
                )
            if state.qc_attempt_status == "succeeded":
                return None
            if state.qc_attempt_status == "failed":
                if state.qc_reason_code == reason_code:
                    return None
                raise ConcurrentRunUpdateError("QC attempt already failed")
            empty: tuple[RunQcMetric, ...] = ()
            revision = state.qc_revision + 1
            generation = build_qc_generation(
                run_id=run_id,
                revision=revision,
                artifact_generation=expected_artifact_generation,
                metrics=empty,
            )
            next_state = replace(
                state,
                qc_revision=revision,
                qc_generation=generation,
                qc_manifest_digest=qc_metric_manifest_digest(empty),
                qc_attempt_status="failed",
                qc_artifact_generation=expected_artifact_generation,
                qc_outcome="failed",
                qc_reason_code=reason_code,
            )
            result_event = self._make_event(
                run_id,
                len(self._events[run_id]) + 1,
                _event_with_context(
                    event,
                    attempt_id=attempt_id,
                    artifact_generation=expected_artifact_generation,
                    qc_generation=generation,
                    reason_code=reason_code,
                ),
            )
            self._qc_metrics[run_id] = {}
            self._result_states[run_id] = next_state
            self._events[run_id].append(result_event)
            return result_event

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

    def _latest_event_type(
        self,
        run_id: str,
        event_types: frozenset[str],
    ) -> str | None:
        return next(
            (
                item.event_type
                for item in reversed(self._events[run_id])
                if item.event_type in event_types
            ),
            None,
        )

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


def _summary_from_record(mapping_key: str, record: RunRecord) -> RunSummary:
    if not isinstance(mapping_key, str) or not isinstance(record, RunRecord):
        raise ValueError("persisted run summary identity is invalid")
    if mapping_key != record.run_id:
        raise ValueError("persisted run summary identity is invalid")
    return RunSummary(
        run_id=record.run_id,
        workflow_id=record.workflow_id,
        status=record.status,
        created_at=record.created_at,
        updated_at=record.updated_at,
        started_at=record.started_at,
        ended_at=record.ended_at,
        current_stage=record.current_stage,
    )


def _summary_matches(
    summary: RunSummary,
    *,
    workflow_id: str | None,
    status: RunStatus | None,
) -> bool:
    return (workflow_id is None or summary.workflow_id == workflow_id) and (
        status is None or summary.status is status
    )


def _summary_is_after(summary: RunSummary, cursor: RunHistoryCursor) -> bool:
    return (summary.created_at, summary.run_id) < (
        cursor.created_at,
        cursor.run_id,
    )


def canonical_decimal_text(value: Decimal) -> str:
    """Serialize a bounded Decimal without SQLite numeric affinity."""
    if not isinstance(value, Decimal) or not value.is_finite():
        raise ValueError("QC value must be a finite Decimal")
    text = format(value, "f")
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    if text in {"-0", ""}:
        text = "0"
    if _CANONICAL_DECIMAL.fullmatch(text) is None:
        raise ValueError("QC value exceeds durable decimal bounds")
    return text


def decimal_from_canonical_text(value: str) -> Decimal:
    """Reconstruct a Decimal only from the durable canonical grammar."""
    if not isinstance(value, str) or _CANONICAL_DECIMAL.fullmatch(value) is None:
        raise ValueError("Persisted QC decimal is invalid")
    try:
        result = Decimal(value)
    except InvalidOperation as exc:
        raise ValueError("Persisted QC decimal is invalid") from exc
    if not result.is_finite() or canonical_decimal_text(result) != value:
        raise ValueError("Persisted QC decimal is not canonical")
    return result


def _require_current_attempt(
    current_attempt_id: str | None,
    expected_attempt_id: str,
    name: str,
) -> None:
    if current_attempt_id != expected_attempt_id:
        raise ConcurrentRunUpdateError(f"{name} result attempt is no longer current")


def _event_with_context(event: RunEventDraft, **context: object) -> RunEventDraft:
    merged = dict(event.context)
    merged.update({key: value for key, value in context.items() if value is not None})
    return replace(event, context=merged)


def _state_after_artifact_change(
    state: RunResultState,
    *,
    artifacts: tuple[RunArtifactRef, ...],
    artifact_revision: int,
    artifact_generation: str,
    attempt_id: str | None,
    attempt_status: str | None,
    outcome: str,
) -> RunResultState:
    had_qc_state = state.qc_generation is not None or state.qc_outcome is not None
    values: dict[str, object] = {
        "artifact_revision": artifact_revision,
        "artifact_generation": artifact_generation,
        "artifact_manifest_digest": artifact_manifest_digest(artifacts),
        "artifact_attempt_id": attempt_id,
        "artifact_attempt_status": attempt_status,
        "artifact_outcome": outcome,
        "artifact_reason_code": None,
        # A changed artifact generation always supersedes a pending QC attempt.
        "qc_attempt_id": None,
        "qc_attempt_status": None,
        "qc_attempt_artifact_generation": None,
    }
    if had_qc_state:
        qc_revision = state.qc_revision + 1
        values.update(
            {
                "qc_revision": qc_revision,
                "qc_generation": None,
                "qc_manifest_digest": None,
                "qc_artifact_generation": None,
                "qc_outcome": "invalidated",
                "qc_reason_code": None,
            }
        )
    return replace(state, **values)


def _sorted_artifacts(
    artifacts: Iterable[RunArtifactRef],
) -> tuple[RunArtifactRef, ...]:
    return tuple(sorted(artifacts, key=lambda artifact: artifact.artifact_id))


def _validated_expected_artifacts(
    run_id: str,
    artifacts: tuple[RunArtifactRef, ...],
) -> tuple[RunArtifactRef, ...]:
    if not isinstance(artifacts, tuple):
        raise ValueError("expected_artifacts must be a tuple")
    seen: set[str] = set()
    for artifact in artifacts:
        if not isinstance(artifact, RunArtifactRef) or artifact.run_id != run_id:
            raise ValueError("expected artifact does not match the run")
        if artifact.artifact_id in seen:
            raise ValueError("expected artifacts contain a duplicate")
        seen.add(artifact.artifact_id)
    return _sorted_artifacts(artifacts)


def _validated_qc_replacement(
    run_id: str,
    metrics: tuple[RunQcMetric, ...],
    source_artifact_ids: set[str],
) -> tuple[RunQcMetric, ...]:
    if not isinstance(metrics, tuple):
        raise ValueError("metrics must be a tuple")
    seen: set[str] = set()
    for metric in metrics:
        if not isinstance(metric, RunQcMetric) or metric.run_id != run_id:
            raise ValueError("QC metric does not match the run")
        if metric.metric_id in seen:
            raise ValueError("QC replacement contains a duplicate metric_id")
        _validate_qc_metric_fields(metric)
        if metric.source_artifact_id not in source_artifact_ids:
            raise ValueError("QC metric source is not in the artifact generation")
        seen.add(metric.metric_id)
    return tuple(sorted(metrics, key=lambda metric: metric.metric_id))


def _validate_qc_metric_fields(metric: RunQcMetric) -> None:
    if (
        not isinstance(metric.metric_id, str)
        or _QC_METRIC_ID.fullmatch(metric.metric_id) is None
    ):
        raise ValueError("QC metric_id is invalid")
    if (
        not isinstance(metric.metric_key, str)
        or len(metric.metric_key) > 128
        or _QC_METRIC_KEY.fullmatch(metric.metric_key) is None
    ):
        raise ValueError("QC metric key is invalid")
    if (
        not isinstance(metric.display_name, str)
        or not 1 <= len(metric.display_name) <= 255
        or not metric.display_name.isprintable()
        or "/" in metric.display_name
        or "\\" in metric.display_name
    ):
        raise ValueError("QC metric display name is invalid")
    canonical_decimal_text(metric.value)
    if not isinstance(metric.unit, str) or metric.unit not in _QC_UNITS:
        raise ValueError("QC metric unit is invalid")
    if not isinstance(metric.scope, str) or metric.scope not in _QC_SCOPES:
        raise ValueError("QC metric scope is invalid")
    for value in (metric.sample_id, metric.experiment_id, metric.assay):
        if value is not None:
            try:
                validate_qc_identifier_token(value)
            except ValueError:
                raise ValueError("QC metric identifier is invalid") from None
    if metric.scope == "run" and (
        metric.sample_id is not None or metric.experiment_id is not None
    ):
        raise ValueError("run QC scope cannot have sample or experiment IDs")
    if metric.scope == "sample" and metric.sample_id is None:
        raise ValueError("sample QC scope requires sample_id")
    if metric.scope == "experiment" and (
        metric.sample_id is not None or metric.experiment_id is None
    ):
        raise ValueError("experiment QC scope requires only experiment_id")
    if metric.qc_flag is not None and (
        not isinstance(metric.qc_flag, str) or metric.qc_flag not in _QC_FLAGS
    ):
        raise ValueError("QC metric flag is invalid")
    if (
        not isinstance(metric.source_artifact_id, str)
        or _QC_SOURCE_ARTIFACT_ID.fullmatch(metric.source_artifact_id) is None
    ):
        raise ValueError("QC metric source artifact ID is invalid")
    if (
        not isinstance(metric.produced_at, datetime)
        or metric.produced_at.tzinfo is None
        or metric.produced_at.utcoffset() is None
    ):
        raise ValueError("QC metric produced_at must be timezone-aware")
    expected_metric_id = build_qc_metric_id(
        metric.metric_key,
        metric.scope,
        metric.sample_id,
        metric.experiment_id,
    )
    if metric.metric_id != expected_metric_id:
        raise ValueError("QC metric_id does not match its semantic coordinates")


def _validate_stored_qc_metric_identity(
    metric: object,
    *,
    run_id: str,
    storage_key: str,
) -> None:
    if (
        not isinstance(metric, RunQcMetric)
        or metric.run_id != run_id
        or metric.metric_id != storage_key
    ):
        raise ValueError("Persisted QC metric identity is invalid")
    _validate_qc_metric_fields(metric)


def _validated_snapshot_copy(
    snapshot: ValidatedInputSnapshot,
) -> ValidatedInputSnapshot:
    if not isinstance(snapshot, ValidatedInputSnapshot):
        raise ValueError("snapshot must be a ValidatedInputSnapshot")
    return replace(snapshot)


def _validate_snapshot_consumption_candidate(
    snapshot: ValidatedInputSnapshot,
    *,
    workflow_id: str,
    expected_build_identity: WorkflowBuildIdentity,
    record: RunRecord,
    consumed_at: datetime,
) -> None:
    if not isinstance(expected_build_identity, WorkflowBuildIdentity) or not (
        snapshot.workflow_build_identity.matches(expected_build_identity)
    ):
        raise ValidatedSnapshotBuildMismatchError(
            "validated snapshot workflow build identity differs"
        )
    if not isinstance(record, RunRecord):
        raise ValueError("record must be a RunRecord")
    if record.workflow_id != workflow_id or record.workflow_id != snapshot.workflow_id:
        raise ValueError("run workflow does not match validated snapshot")
    if record.status is not RunStatus.CREATED:
        raise ValueError("validated snapshot can only create a CREATED run")
    if any(
        value is not None
        for value in (
            record.started_at,
            record.ended_at,
            record.current_stage,
            record.cancellation_reason,
            record.error,
        )
    ):
        raise ValueError("new validated run contains lifecycle evidence")
    if (
        not isinstance(consumed_at, datetime)
        or consumed_at.tzinfo is None
        or consumed_at.utcoffset() is None
    ):
        raise ValueError("consumed_at must be timezone-aware")
    if record.created_at != consumed_at or record.updated_at != consumed_at:
        raise ValueError("new validated run timestamps must match consumption time")
    candidate_payload = canonical_workflow_inputs_json(
        WorkflowInputs(
            config=record.inputs.get("config", {}),
            samples=record.inputs.get("samples"),
            options=record.inputs.get("options", {}),
        )
    )
    if set(record.inputs) != {"config", "samples", "options"} or (
        candidate_payload != snapshot.canonical_payload
    ):
        raise ValueError("run inputs do not match validated snapshot")


def _validate_snapshot_linked_run(
    snapshot: ValidatedInputSnapshot,
    record: RunRecord,
) -> None:
    if (
        snapshot.consumed_run_id is None
        or record.run_id != snapshot.consumed_run_id
        or record.workflow_id != snapshot.workflow_id
    ):
        raise ValueError("validated snapshot linked run identity is invalid")
    if (
        snapshot.consumed_at is None
        or snapshot.consumed_at >= snapshot.expires_at
        or record.created_at != snapshot.consumed_at
    ):
        raise ValueError("validated snapshot linked run consumption time is invalid")
    payload = canonical_workflow_inputs_json(
        WorkflowInputs(
            config=record.inputs.get("config", {}),
            samples=record.inputs.get("samples"),
            options=record.inputs.get("options", {}),
        )
    )
    if set(record.inputs) != {"config", "samples", "options"} or (
        payload != snapshot.canonical_payload
    ):
        raise ValueError("validated snapshot linked run inputs are invalid")


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
