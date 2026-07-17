"""SQLAlchemy implementation of the workflow run repository boundary."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import replace
from datetime import datetime, timezone
from threading import RLock

from sqlalchemy import and_, delete, func, or_, select, update
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session, sessionmaker

from encode_pipeline.persistence.models import (
    RunArtifactRow,
    RunEventRow,
    RunExecutionAssignmentRow,
    RunLogRow,
    RunQcMetricRow,
    RunResultAttemptRow,
    RunResultStateRow,
    RunRow,
    RunWorkflowBuildIdentityRow,
    ValidatedInputSnapshotRow,
)
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.execution import (
    RunExecutionAssignment,
    RunExecutionCancellationRequest,
    RunExecutionClaim,
    RunExecutionOwnership,
    RunExecutionStopAcknowledgement,
)
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
)
from encode_pipeline.platform.snapshots import (
    ValidatedInputSnapshot,
    ValidatedSnapshotRunCreation,
)
from encode_pipeline.services.run_repositories import (
    ConcurrentRunUpdateError,
    ResultGenerationChangedError,
    RunEventDraft,
    ValidatedSnapshotExpiredError,
    ValidatedSnapshotReplayConflictError,
    _event_with_context,
    _require_current_attempt,
    _sorted_artifacts,
    _state_after_artifact_change,
    _validate_qc_metric_fields,
    _validated_expected_artifacts,
    _validated_qc_replacement,
    _validate_snapshot_consumption_candidate,
    _validate_snapshot_linked_run,
    canonical_decimal_text,
    decimal_from_canonical_text,
)


class SqlAlchemyRunRepository:
    """Persist run aggregates without exposing ORM rows to service callers."""

    def __init__(self, session_factory: sessionmaker[Session]) -> None:
        if not isinstance(session_factory, sessionmaker):
            raise ValueError("session_factory must be a SQLAlchemy sessionmaker")
        self._session_factory = session_factory
        self._lock = RLock()

    def contains_run(self, run_id: str) -> bool:
        with self._session_factory() as session:
            return (
                session.scalar(select(RunRow.id).where(RunRow.run_id == run_id))
                is not None
            )

    def create_run(self, record: RunRecord, event: RunEventDraft) -> RunEvent:
        with self._lock:
            try:
                with self._session_factory.begin() as session:
                    _begin_write(session)
                    session.add(_run_row(record))
                    session.flush()
                    session.add(RunResultStateRow(run_id=record.run_id))
                    session.flush()
                    return self._insert_event(session, record.run_id, event)
            except IntegrityError as exc:
                raise ValueError(f"Duplicate run_id: {record.run_id!r}") from exc

    def create_validated_input_snapshot(
        self,
        snapshot: ValidatedInputSnapshot,
    ) -> ValidatedInputSnapshot:
        validated = _validated_input_snapshot_from_row(
            _validated_input_snapshot_row(snapshot)
        )
        try:
            with self._lock, self._session_factory.begin() as session:
                _begin_write(session)
                session.add(_validated_input_snapshot_row(validated))
                session.flush()
        except IntegrityError as exc:
            raise ValueError(
                f"Duplicate validated snapshot ID: {validated.snapshot_id!r}"
            ) from exc
        return validated

    def get_validated_input_snapshot(
        self,
        snapshot_id: str,
    ) -> ValidatedInputSnapshot:
        with self._session_factory() as session:
            row = session.get(ValidatedInputSnapshotRow, snapshot_id)
            if row is None:
                raise KeyError(snapshot_id)
            return _validated_input_snapshot_from_row(row)

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
        try:
            with self._lock, self._session_factory.begin() as session:
                _begin_write(session)
                row = session.get(ValidatedInputSnapshotRow, snapshot_id)
                if row is None:
                    raise KeyError(snapshot_id)
                snapshot = _validated_input_snapshot_from_row(row)
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
                    run_row = session.scalar(
                        select(RunRow).where(RunRow.run_id == snapshot.consumed_run_id)
                    )
                    if run_row is None:
                        raise ValueError("validated snapshot references a missing run")
                    current = _record_from_row(run_row)
                    _validate_snapshot_linked_run(snapshot, current)
                    if dict(current.tags) != dict(record.tags):
                        raise ValidatedSnapshotReplayConflictError(
                            "validated snapshot replay metadata differs"
                        )
                    return ValidatedSnapshotRunCreation(
                        record=current,
                        created=False,
                    )
                if consumed_at >= snapshot.expires_at:
                    raise ValidatedSnapshotExpiredError(
                        "validated snapshot expired before first use"
                    )

                session.add(_run_row(record))
                session.flush()
                session.add(RunResultStateRow(run_id=record.run_id))
                session.flush()
                created_event = self._insert_event(
                    session,
                    record.run_id,
                    event,
                )
                if created_event.sequence != 1:
                    raise ValueError("validated run initial event is invalid")
                row.consumed_run_id = record.run_id
                row.consumed_at = consumed_at
                session.flush()
                return ValidatedSnapshotRunCreation(record=record, created=True)
        except IntegrityError as exc:
            raise ValueError("validated snapshot could not create a run") from exc

    def get_run(self, run_id: str) -> RunRecord:
        with self._session_factory() as session:
            return _record_from_row(self._require_run(session, run_id))

    def list_runs(self) -> tuple[RunRecord, ...]:
        with self._session_factory() as session:
            rows = session.scalars(select(RunRow).order_by(RunRow.id)).all()
            return tuple(_record_from_row(row) for row in rows)

    def list_run_summaries(
        self,
        *,
        after: RunHistoryCursor | None = None,
        limit: int = 50,
        workflow_id: str | None = None,
        status: RunStatus | None = None,
    ) -> tuple[RunSummary, ...]:
        """Return a bounded summary-only keyset query from SQLite."""
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

        with self._session_factory() as session:
            if after is not None:
                boundary_row = session.execute(
                    _run_summary_select().where(RunRow.run_id == after.run_id)
                ).one_or_none()
                if boundary_row is None:
                    raise KeyError(after.run_id)
                boundary = _summary_from_selected_row(boundary_row)
                if boundary.created_at != after.created_at or not _summary_matches(
                    boundary,
                    workflow_id=workflow_id,
                    status=status,
                ):
                    raise KeyError(after.run_id)

            statement = _run_summary_select()
            if workflow_id is not None:
                statement = statement.where(RunRow.workflow_id == workflow_id)
            if status is not None:
                statement = statement.where(RunRow.status == status.value)
            if after is not None:
                statement = statement.where(
                    or_(
                        RunRow.created_at < after.created_at,
                        and_(
                            RunRow.created_at == after.created_at,
                            RunRow.run_id < after.run_id,
                        ),
                    )
                )
            rows = session.execute(
                statement.order_by(
                    RunRow.created_at.desc(),
                    RunRow.run_id.desc(),
                ).limit(limit)
            ).all()
            return tuple(_summary_from_selected_row(row) for row in rows)

    def update_run(
        self,
        record: RunRecord,
        *,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent:
        with self._lock, self._session_factory.begin() as session:
            _begin_write(session)
            result = session.execute(
                update(RunRow)
                .where(
                    RunRow.run_id == record.run_id,
                    RunRow.status == expected_status.value,
                )
                .values(**_run_values(record))
            )
            if result.rowcount != 1:
                if (
                    session.scalar(
                        select(RunRow.id).where(RunRow.run_id == record.run_id)
                    )
                    is None
                ):
                    raise KeyError(record.run_id)
                raise ConcurrentRunUpdateError(
                    f"Run {record.run_id!r} changed while it was being updated."
                )
            return self._insert_event(session, record.run_id, event)

    def complete_preflight(
        self,
        record: RunRecord,
        identity: WorkflowBuildIdentity,
        *,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent:
        """Atomically persist a build identity, PLANNED run, and event."""
        if record.status is not RunStatus.PLANNED:
            raise ValueError("completed preflight record must be planned")
        if identity.workflow_id != record.workflow_id:
            raise ValueError("workflow build identity does not match the run")

        try:
            with self._lock, self._session_factory.begin() as session:
                _begin_write(session)
                if session.get(RunWorkflowBuildIdentityRow, record.run_id) is not None:
                    raise ValueError("run already has a workflow build identity")
                result = session.execute(
                    update(RunRow)
                    .where(
                        RunRow.run_id == record.run_id,
                        RunRow.status == expected_status.value,
                    )
                    .values(**_run_values(record))
                )
                if result.rowcount != 1:
                    if (
                        session.scalar(
                            select(RunRow.id).where(RunRow.run_id == record.run_id)
                        )
                        is None
                    ):
                        raise KeyError(record.run_id)
                    raise ConcurrentRunUpdateError(
                        f"Run {record.run_id!r} changed while preflight completed."
                    )
                session.add(_workflow_build_identity_row(record.run_id, identity))
                session.flush()
                return self._insert_event(session, record.run_id, event)
        except IntegrityError as exc:
            raise ValueError("workflow build identity could not be persisted") from exc

    def get_workflow_build_identity(
        self,
        run_id: str,
    ) -> WorkflowBuildIdentity | None:
        with self._session_factory() as session:
            row = session.get(RunWorkflowBuildIdentityRow, run_id)
            return _workflow_build_identity_from_row(row) if row is not None else None

    def add_event(self, run_id: str, event: RunEventDraft) -> RunEvent:
        with self._lock, self._session_factory.begin() as session:
            _begin_write(session)
            self._require_run(session, run_id)
            return self._insert_event(session, run_id, event)

    def fail_interrupted_run_if_unowned(
        self,
        record: RunRecord,
        *,
        expected_status: RunStatus,
        required_ownership: RunExecutionOwnership | None,
        event: RunEventDraft,
    ) -> bool:
        with self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, record.run_id)
            if RunStatus(current.status) is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {record.run_id!r} changed while it was being recovered."
                )
            assignment = session.get(RunExecutionAssignmentRow, record.run_id)
            if _assignment_row_has_ownership(assignment, required_ownership):
                return False
            result = session.execute(
                update(RunRow)
                .where(
                    RunRow.run_id == record.run_id,
                    RunRow.status == expected_status.value,
                )
                .values(**_run_values(record))
            )
            if result.rowcount != 1:
                raise ConcurrentRunUpdateError(
                    f"Run {record.run_id!r} changed while it was being recovered."
                )
            self._insert_event(session, record.run_id, event)
            return True

    def list_events(
        self,
        run_id: str,
        *,
        after: str | None = None,
        limit: int = 50,
    ) -> tuple[RunEvent, ...]:
        with self._session_factory() as session:
            self._require_run(session, run_id)
            statement = select(RunEventRow).where(RunEventRow.run_id == run_id)
            if after is not None:
                cursor_sequence = session.scalar(
                    select(RunEventRow.sequence).where(
                        RunEventRow.run_id == run_id,
                        RunEventRow.event_id == after,
                    )
                )
                if cursor_sequence is None:
                    raise KeyError(after)
                statement = statement.where(RunEventRow.sequence > cursor_sequence)
            rows = session.scalars(
                statement.order_by(RunEventRow.sequence).limit(limit)
            ).all()
            return tuple(_event_from_row(row) for row in rows)

    def append_log(
        self,
        run_id: str,
        stream_name: str,
        lines: Iterable[str],
    ) -> RunLogChunk:
        normalized_lines = tuple(lines)
        with self._lock, self._session_factory.begin() as session:
            _begin_write(session)
            self._require_run(session, run_id)
            sequence = self._next_log_sequence(session, run_id, stream_name)
            timestamp = datetime.now(timezone.utc)
            row = RunLogRow(
                chunk_id=f"log-{sequence}",
                run_id=run_id,
                stream_name=stream_name,
                sequence=sequence,
                timestamp=timestamp,
                lines=list(normalized_lines),
            )
            session.add(row)
            session.flush()
            return RunLogChunk(
                chunk_id=row.chunk_id,
                run_id=run_id,
                stream_name=stream_name,
                sequence=sequence,
                timestamp=timestamp,
                lines=normalized_lines,
            )

    def list_logs(
        self,
        run_id: str,
        stream_name: str = "stdout",
        *,
        after: str | None = None,
        limit: int = 50,
    ) -> tuple[RunLogChunk, ...]:
        with self._session_factory() as session:
            self._require_run(session, run_id)
            statement = select(RunLogRow).where(
                RunLogRow.run_id == run_id,
                RunLogRow.stream_name == stream_name,
            )
            if after is not None:
                cursor_sequence = session.scalar(
                    select(RunLogRow.sequence).where(
                        RunLogRow.run_id == run_id,
                        RunLogRow.stream_name == stream_name,
                        RunLogRow.chunk_id == after,
                    )
                )
                if cursor_sequence is None:
                    raise KeyError(after)
                statement = statement.where(RunLogRow.sequence > cursor_sequence)
            rows = session.scalars(
                statement.order_by(RunLogRow.sequence).limit(limit)
            ).all()
            return tuple(_log_from_row(row) for row in rows)

    def get_result_state(self, run_id: str) -> RunResultState:
        with self._session_factory() as session:
            self._require_run(session, run_id)
            return _result_state_from_row(self._require_result_state(session, run_id))

    def begin_artifact_result_attempt(
        self,
        run_id: str,
        *,
        attempt_id: str,
        expected_status: RunStatus,
    ) -> RunResultState:
        validate_result_attempt_id(attempt_id)
        with self._lock, self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            if RunStatus(current.status) is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer eligible for artifact indexing."
                )
            row = self._require_result_state(session, run_id)
            state = _result_state_from_row(row)
            previous = session.get(RunResultAttemptRow, attempt_id)
            if previous is not None:
                if (
                    previous.run_id == run_id
                    and previous.result_kind == "artifact"
                    and previous.artifact_generation is None
                    and state.artifact_attempt_id == attempt_id
                ):
                    return state
                raise ConcurrentRunUpdateError(
                    "artifact result attempt was already superseded"
                )
            session.add(
                RunResultAttemptRow(
                    attempt_id=attempt_id,
                    run_id=run_id,
                    result_kind="artifact",
                    artifact_generation=None,
                )
            )
            session.flush()
            state = replace(
                state,
                artifact_attempt_id=attempt_id,
                artifact_attempt_status="pending",
            )
            _apply_result_state(row, state)
            session.flush()
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
        with self._lock, self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            if RunStatus(current.status) is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer eligible for QC indexing."
                )
            row = self._require_result_state(session, run_id)
            state = _result_state_from_row(row)
            if state.artifact_generation != expected_artifact_generation:
                raise ResultGenerationChangedError(
                    f"Run {run_id!r} artifact generation changed before QC indexing."
                )
            current_artifacts = _sorted_artifacts(
                _artifact_from_row(artifact_row)
                for artifact_row in session.scalars(
                    select(RunArtifactRow).where(RunArtifactRow.run_id == run_id)
                ).all()
            )
            if current_artifacts != _validated_expected_artifacts(
                run_id,
                expected_artifacts,
            ):
                raise ResultGenerationChangedError(
                    f"Run {run_id!r} artifact manifest changed before QC indexing."
                )
            previous = session.get(RunResultAttemptRow, attempt_id)
            if previous is not None:
                if (
                    previous.run_id == run_id
                    and previous.result_kind == "qc"
                    and previous.artifact_generation == expected_artifact_generation
                    and state.qc_attempt_id == attempt_id
                    and state.qc_attempt_artifact_generation
                    == expected_artifact_generation
                ):
                    return state
                raise ConcurrentRunUpdateError(
                    "QC result attempt was already superseded"
                )
            session.add(
                RunResultAttemptRow(
                    attempt_id=attempt_id,
                    run_id=run_id,
                    result_kind="qc",
                    artifact_generation=expected_artifact_generation,
                )
            )
            session.flush()
            state = replace(
                state,
                qc_attempt_id=attempt_id,
                qc_attempt_status="pending",
                qc_attempt_artifact_generation=expected_artifact_generation,
            )
            _apply_result_state(row, state)
            session.flush()
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
        try:
            with self._lock, self._session_factory.begin() as session:
                _begin_write(session)
                current = self._require_run(session, run_id)
                if RunStatus(current.status) is not expected_status:
                    raise ConcurrentRunUpdateError(
                        f"Run {run_id!r} is no longer eligible for artifact replacement."
                    )
                sorted_replacement = _validated_expected_artifacts(run_id, artifacts)
                state_row = self._require_result_state(session, run_id)
                state = _result_state_from_row(state_row)
                _require_current_attempt(
                    state.artifact_attempt_id, attempt_id, "artifact"
                )
                manifest_digest = artifact_manifest_digest(sorted_replacement)
                if state.artifact_attempt_status == "succeeded":
                    if state.artifact_manifest_digest == manifest_digest:
                        return None
                    raise ConcurrentRunUpdateError("artifact attempt already committed")
                if state.artifact_attempt_status != "pending":
                    raise ConcurrentRunUpdateError(
                        "artifact attempt is no longer current"
                    )
                changed = state.artifact_manifest_digest != manifest_digest
                should_emit = changed or state.artifact_outcome != "succeeded"
                had_qc_rows = (
                    session.scalar(
                        select(func.count())
                        .select_from(RunQcMetricRow)
                        .where(RunQcMetricRow.run_id == run_id)
                    )
                    or 0
                ) > 0
                had_qc_state = (
                    had_qc_rows
                    or state.qc_generation is not None
                    or (state.qc_outcome is not None)
                )
                if changed:
                    revision = state.artifact_revision + 1
                    generation = build_artifact_generation(
                        run_id=run_id,
                        revision=revision,
                        artifacts=sorted_replacement,
                    )
                    session.execute(
                        delete(RunArtifactRow).where(RunArtifactRow.run_id == run_id)
                    )
                    session.add_all(
                        [_artifact_row(artifact) for artifact in sorted_replacement]
                    )
                    session.execute(
                        delete(RunQcMetricRow).where(RunQcMetricRow.run_id == run_id)
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
                    _apply_result_state(state_row, state)
                    session.flush()
                    if had_qc_state:
                        self._insert_event(
                            session,
                            run_id,
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
                else:
                    state = replace(
                        state,
                        artifact_attempt_status="succeeded",
                        artifact_outcome="succeeded",
                        artifact_reason_code=None,
                    )
                    _apply_result_state(state_row, state)
                    session.flush()
                if not should_emit:
                    return None
                return self._insert_event(
                    session,
                    run_id,
                    _event_with_context(
                        event,
                        attempt_id=attempt_id,
                        artifact_generation=state.artifact_generation,
                        artifact_count=len(sorted_replacement),
                    ),
                )
        except IntegrityError as exc:
            raise ValueError("artifact replacement could not be persisted") from exc

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
        with self._lock, self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            if RunStatus(current.status) is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer eligible for artifact failure."
                )
            row = self._require_result_state(session, run_id)
            state = _result_state_from_row(row)
            _require_current_attempt(state.artifact_attempt_id, attempt_id, "artifact")
            if state.artifact_attempt_status == "succeeded":
                return None
            if state.artifact_attempt_status == "failed":
                if state.artifact_reason_code == reason_code:
                    return None
                raise ConcurrentRunUpdateError("artifact attempt already failed")
            state = replace(
                state,
                artifact_attempt_status="failed",
                artifact_outcome="failed",
                artifact_reason_code=reason_code,
            )
            _apply_result_state(row, state)
            session.flush()
            return self._insert_event(
                session,
                run_id,
                _event_with_context(
                    event, attempt_id=attempt_id, reason_code=reason_code
                ),
            )

    def list_artifacts(
        self,
        run_id: str,
        *,
        after: str | None = None,
        limit: int | None = None,
    ) -> tuple[RunArtifactRef, ...]:
        with self._session_factory() as session:
            self._require_run(session, run_id)
            state = _result_state_from_row(self._require_result_state(session, run_id))
            has_rows = session.scalar(
                select(func.count())
                .select_from(RunArtifactRow)
                .where(RunArtifactRow.run_id == run_id)
            )
            if has_rows and state.artifact_generation is None:
                raise ValueError("artifact generation is unbound")
            if after is not None:
                cursor = session.scalar(
                    select(RunArtifactRow.artifact_id).where(
                        RunArtifactRow.run_id == run_id,
                        RunArtifactRow.artifact_id == after,
                    )
                )
                if cursor is None:
                    raise KeyError((run_id, after))
            statement = (
                select(RunArtifactRow)
                .where(
                    RunArtifactRow.run_id == run_id,
                    *(
                        (RunArtifactRow.artifact_id > after,)
                        if after is not None
                        else ()
                    ),
                )
                .order_by(RunArtifactRow.artifact_id)
            )
            if limit is not None:
                statement = statement.limit(limit)
            rows = session.scalars(statement).all()
            return tuple(_artifact_from_row(row) for row in rows)

    def get_artifact(self, run_id: str, artifact_id: str) -> RunArtifactRef:
        with self._session_factory() as session:
            self._require_run(session, run_id)
            state = _result_state_from_row(self._require_result_state(session, run_id))
            row = session.scalar(
                select(RunArtifactRow).where(
                    RunArtifactRow.run_id == run_id,
                    RunArtifactRow.artifact_id == artifact_id,
                )
            )
            if row is None:
                raise KeyError((run_id, artifact_id))
            if state.artifact_generation is None:
                raise ValueError("artifact generation is unbound")
            return _artifact_from_row(row)

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
        try:
            with self._lock, self._session_factory.begin() as session:
                _begin_write(session)
                current = self._require_run(session, run_id)
                if RunStatus(current.status) is not expected_status:
                    raise ConcurrentRunUpdateError(
                        f"Run {run_id!r} is no longer eligible for QC replacement."
                    )
                state_row = self._require_result_state(session, run_id)
                state = _result_state_from_row(state_row)
                _require_current_attempt(state.qc_attempt_id, attempt_id, "QC")
                if (
                    state.artifact_generation != expected_artifact_generation
                    or state.qc_attempt_artifact_generation
                    != expected_artifact_generation
                ):
                    raise ResultGenerationChangedError(
                        f"Run {run_id!r} artifact generation changed during QC indexing."
                    )
                current_artifacts = _sorted_artifacts(
                    _artifact_from_row(row)
                    for row in session.scalars(
                        select(RunArtifactRow).where(RunArtifactRow.run_id == run_id)
                    ).all()
                )
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
                if not changed:
                    state = replace(state, qc_attempt_status="succeeded")
                    _apply_result_state(state_row, state)
                    session.flush()
                    return None
                revision = state.qc_revision + 1
                generation = build_qc_generation(
                    run_id=run_id,
                    revision=revision,
                    artifact_generation=expected_artifact_generation,
                    metrics=replacement,
                )
                session.execute(
                    delete(RunQcMetricRow).where(RunQcMetricRow.run_id == run_id)
                )
                session.add_all([_qc_metric_row(metric) for metric in replacement])
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
                _apply_result_state(state_row, state)
                session.flush()
                return self._insert_event(
                    session,
                    run_id,
                    _event_with_context(
                        event,
                        attempt_id=attempt_id,
                        artifact_generation=expected_artifact_generation,
                        qc_generation=generation,
                        metric_count=len(replacement),
                    ),
                )
        except IntegrityError as exc:
            raise ValueError("QC replacement could not be persisted") from exc

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
        with self._session_factory() as session:
            _begin_consistent_read(session)
            self._require_run(session, run_id)
            state = _result_state_from_row(self._require_result_state(session, run_id))
            if (
                expected_generation is not None
                and state.qc_generation != expected_generation
            ):
                raise ResultGenerationChangedError("QC generation changed")
            has_rows = session.scalar(
                select(func.count())
                .select_from(RunQcMetricRow)
                .where(RunQcMetricRow.run_id == run_id)
            )
            if has_rows and (
                state.qc_generation is None
                or state.qc_artifact_generation != state.artifact_generation
            ):
                raise ValueError("QC generation is unbound")
            if after is not None:
                cursor_row = session.scalar(
                    select(RunQcMetricRow).where(
                        RunQcMetricRow.run_id == run_id,
                        RunQcMetricRow.metric_id == after,
                    )
                )
                if cursor_row is None:
                    raise KeyError((run_id, after))
                _qc_metric_from_row(cursor_row)
            statement = (
                select(RunQcMetricRow)
                .where(
                    RunQcMetricRow.run_id == run_id,
                    *((RunQcMetricRow.metric_id > after,) if after is not None else ()),
                )
                .order_by(RunQcMetricRow.metric_id)
            )
            if limit is not None:
                statement = statement.limit(limit)
            rows = session.scalars(statement).all()
            return state.qc_generation, tuple(_qc_metric_from_row(row) for row in rows)

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
        with self._lock, self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            if RunStatus(current.status) is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer eligible for QC failure."
                )
            row = self._require_result_state(session, run_id)
            state = _result_state_from_row(row)
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
            session.execute(
                delete(RunQcMetricRow).where(RunQcMetricRow.run_id == run_id)
            )
            state = replace(
                state,
                qc_revision=revision,
                qc_generation=generation,
                qc_manifest_digest=qc_metric_manifest_digest(empty),
                qc_attempt_status="failed",
                qc_artifact_generation=expected_artifact_generation,
                qc_outcome="failed",
                qc_reason_code=reason_code,
            )
            _apply_result_state(row, state)
            session.flush()
            return self._insert_event(
                session,
                run_id,
                _event_with_context(
                    event,
                    attempt_id=attempt_id,
                    artifact_generation=expected_artifact_generation,
                    qc_generation=generation,
                    reason_code=reason_code,
                ),
            )

    def ensure_execution_assignment(
        self,
        assignment: RunExecutionAssignment,
        *,
        expected_status: RunStatus,
    ) -> RunExecutionAssignment:
        try:
            with self._session_factory.begin() as session:
                _begin_write(session)
                current = self._require_run(session, assignment.run_id)
                if RunStatus(current.status) is not expected_status:
                    raise ConcurrentRunUpdateError(
                        f"Run {assignment.run_id!r} is no longer assignable."
                    )
                row = RunExecutionAssignmentRow(
                    run_id=assignment.run_id,
                    job_id=assignment.job_id,
                    backend=assignment.backend,
                    queue_name=assignment.queue_name,
                    created_at=assignment.created_at,
                    dispatched_at=assignment.dispatched_at,
                    claimed_at=assignment.claimed_at,
                    cancellation_requested_at=(assignment.cancellation_requested_at),
                    cancellation_reason=assignment.cancellation_reason,
                    cancellation_acknowledged_at=(
                        assignment.cancellation_acknowledged_at
                    ),
                )
                session.add(row)
                session.flush()
                return _execution_assignment_from_row(row)
        except IntegrityError as exc:
            with self._session_factory() as session:
                existing = session.get(
                    RunExecutionAssignmentRow,
                    assignment.run_id,
                )
                if existing is not None:
                    return _execution_assignment_from_row(existing)
                if (
                    session.scalar(
                        select(RunRow.run_id).where(RunRow.run_id == assignment.run_id)
                    )
                    is None
                ):
                    raise KeyError(assignment.run_id) from exc
                assigned_run_id = session.scalar(
                    select(RunExecutionAssignmentRow.run_id).where(
                        RunExecutionAssignmentRow.job_id == assignment.job_id
                    )
                )
                if assigned_run_id is not None:
                    raise ValueError(
                        f"Execution job_id {assignment.job_id!r} is already "
                        f"assigned to run {assigned_run_id!r}."
                    ) from exc
            raise

    def get_execution_assignment(
        self,
        run_id: str,
    ) -> RunExecutionAssignment | None:
        with self._session_factory() as session:
            row = session.get(RunExecutionAssignmentRow, run_id)
            return _execution_assignment_from_row(row) if row is not None else None

    def mark_execution_dispatched(
        self,
        run_id: str,
        *,
        job_id: str,
        dispatched_at: datetime,
        allowed_statuses: frozenset[RunStatus],
    ) -> RunExecutionAssignment:
        with self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            row = session.get(RunExecutionAssignmentRow, run_id)
            if row is None:
                raise KeyError(run_id)
            if row.job_id != job_id:
                raise ValueError("job_id does not match the execution assignment")
            if row.dispatched_at is not None:
                return _execution_assignment_from_row(row)
            if RunStatus(current.status) not in allowed_statuses:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer dispatchable."
                )
            row.dispatched_at = dispatched_at
            session.flush()
            return _execution_assignment_from_row(row)

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
        if expected_status is not RunStatus.PLANNED:
            raise ValueError("expected_status must be planned")
        if record.status is not RunStatus.QUEUED:
            raise ValueError("queued record must have queued status")

        with self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, record.run_id)
            assignment = session.get(RunExecutionAssignmentRow, record.run_id)
            if assignment is None:
                raise KeyError(record.run_id)
            _require_assignment_row_identity(
                assignment,
                job_id=job_id,
                backend=backend,
                queue_name=queue_name,
            )
            if assignment.dispatched_at is None:
                raise ValueError("execution assignment has not been dispatched")
            if RunStatus(current.status) is not expected_status:
                return False

            result = session.execute(
                update(RunRow)
                .where(
                    RunRow.run_id == record.run_id,
                    RunRow.status == expected_status.value,
                )
                .values(**_run_values(record))
            )
            if result.rowcount != 1:
                raise ConcurrentRunUpdateError(
                    f"Run {record.run_id!r} changed while it was being queued."
                )
            self._insert_event(session, record.run_id, event)
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
        with self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            row = session.get(RunExecutionAssignmentRow, run_id)
            if row is None:
                raise KeyError(run_id)
            _require_assignment_row_identity(
                row,
                job_id=job_id,
                backend=backend,
                queue_name=queue_name,
            )
            assignment = _execution_assignment_from_row(row)
            if assignment.claimed_at is not None:
                return RunExecutionClaim(assignment=assignment, acquired=False)

            current_status = RunStatus(current.status)
            if current_status not in allowed_statuses:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer claimable."
                )
            row.dispatched_at = row.dispatched_at or claimed_at
            row.claimed_at = claimed_at
            session.flush()
            self._insert_event(
                session,
                run_id,
                replace(event, status=current_status),
            )
            return RunExecutionClaim(
                assignment=_execution_assignment_from_row(row),
                acquired=True,
            )

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
        with self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            row = session.get(RunExecutionAssignmentRow, run_id)
            if row is None:
                raise KeyError(run_id)
            _require_assignment_row_identity(
                row,
                job_id=job_id,
                backend=backend,
                queue_name=queue_name,
            )
            assignment = _execution_assignment_from_row(row)
            if assignment.claimed_at is None:
                raise ValueError("execution assignment has not been claimed")
            record = _record_from_row(current)
            if record.status is not RunStatus.RUNNING:
                if record.status.is_terminal:
                    return RunExecutionCancellationRequest(
                        assignment=assignment,
                        record=record,
                        created=False,
                    )
                raise ConcurrentRunUpdateError(f"Run {run_id!r} is no longer running.")
            if assignment.cancellation_requested_at is not None:
                return RunExecutionCancellationRequest(
                    assignment=assignment,
                    record=record,
                    created=False,
                )

            row.cancellation_requested_at = requested_at
            row.cancellation_reason = reason
            self._insert_event(session, run_id, event)
            session.flush()
            return RunExecutionCancellationRequest(
                assignment=_execution_assignment_from_row(row),
                record=record,
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
        with self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            row = session.get(RunExecutionAssignmentRow, run_id)
            if row is None:
                raise KeyError(run_id)
            _require_assignment_row_identity(
                row,
                job_id=job_id,
                backend=backend,
                queue_name=queue_name,
            )
            assignment = _execution_assignment_from_row(row)
            record = _record_from_row(current)
            if record.status.is_terminal:
                return RunExecutionStopAcknowledgement(
                    assignment=assignment,
                    record=record,
                    transitioned=False,
                )

            if assignment.cancellation_requested_at is not None:
                if assignment.claimed_at is None:
                    raise ValueError("execution assignment has not been claimed")
                if record.status is not RunStatus.RUNNING:
                    raise ConcurrentRunUpdateError(
                        f"Run {run_id!r} is no longer running."
                    )
                reason = assignment.cancellation_reason
                assert reason is not None
                row.cancellation_acknowledged_at = acknowledged_at
                updated = replace(
                    record,
                    status=RunStatus.CANCELLED,
                    updated_at=acknowledged_at,
                    ended_at=acknowledged_at,
                    cancellation_reason=reason,
                    error=None,
                )
                context = dict(cancellation_event.context)
                context.update(
                    {
                        "previous_status": RunStatus.RUNNING.value,
                        "new_status": RunStatus.CANCELLED.value,
                        "cancellation_reason": reason,
                    }
                )
                event = replace(cancellation_event, context=context)
            else:
                if record.status not in {
                    RunStatus.PLANNED,
                    RunStatus.QUEUED,
                    RunStatus.RUNNING,
                }:
                    raise ConcurrentRunUpdateError(
                        f"Run {run_id!r} is not worker-owned."
                    )
                updated = replace(
                    record,
                    status=RunStatus.FAILED,
                    updated_at=acknowledged_at,
                    ended_at=acknowledged_at,
                    cancellation_reason=None,
                    error=unexpected_stop_event.issue,
                )
                context = dict(unexpected_stop_event.context)
                context.update(
                    {
                        "previous_status": record.status.value,
                        "new_status": RunStatus.FAILED.value,
                    }
                )
                event = replace(unexpected_stop_event, context=context)

            result = session.execute(
                update(RunRow)
                .where(
                    RunRow.run_id == run_id,
                    RunRow.status == record.status.value,
                )
                .values(**_run_values(updated))
            )
            if result.rowcount != 1:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} changed while stop was acknowledged."
                )
            self._insert_event(session, run_id, event)
            session.flush()
            return RunExecutionStopAcknowledgement(
                assignment=_execution_assignment_from_row(row),
                record=updated,
                transitioned=True,
            )

    @staticmethod
    def _require_run(session: Session, run_id: str) -> RunRow:
        row = session.scalar(select(RunRow).where(RunRow.run_id == run_id))
        if row is None:
            raise KeyError(run_id)
        return row

    @staticmethod
    def _require_result_state(session: Session, run_id: str) -> RunResultStateRow:
        row = session.get(RunResultStateRow, run_id)
        if row is None:
            raise ValueError("run result generation state is missing")
        return row

    @staticmethod
    def _next_log_sequence(
        session: Session,
        run_id: str,
        stream_name: str,
    ) -> int:
        latest = session.scalar(
            select(func.max(RunLogRow.sequence)).where(
                RunLogRow.run_id == run_id,
                RunLogRow.stream_name == stream_name,
            )
        )
        return (latest or 0) + 1

    @staticmethod
    def _insert_event(
        session: Session,
        run_id: str,
        draft: RunEventDraft,
    ) -> RunEvent:
        latest = session.scalar(
            select(func.max(RunEventRow.sequence)).where(RunEventRow.run_id == run_id)
        )
        sequence = (latest or 0) + 1
        timestamp = datetime.now(timezone.utc)
        row = RunEventRow(
            event_id=f"evt-{sequence}",
            run_id=run_id,
            sequence=sequence,
            event_type=draft.event_type,
            timestamp=timestamp,
            status=draft.status.value if draft.status is not None else None,
            stage=draft.stage,
            message=draft.message,
            context=dict(draft.context),
            issue=_issue_to_json(draft.issue),
        )
        session.add(row)
        session.flush()
        return RunEvent(
            event_id=row.event_id,
            run_id=run_id,
            sequence=sequence,
            event_type=draft.event_type,
            timestamp=timestamp,
            status=draft.status,
            stage=draft.stage,
            message=draft.message,
            context=draft.context,
            issue=draft.issue,
        )


def _run_row(record: RunRecord) -> RunRow:
    return RunRow(run_id=record.run_id, **_run_values(record))


def _workflow_build_identity_row(
    run_id: str,
    identity: WorkflowBuildIdentity,
) -> RunWorkflowBuildIdentityRow:
    return RunWorkflowBuildIdentityRow(
        run_id=run_id,
        workflow_id=identity.workflow_id,
        adapter_version=identity.adapter_version,
        scheme=identity.scheme,
        logical_entrypoint=identity.logical_entrypoint,
        digest=identity.digest,
        captured_at=identity.captured_at,
    )


def _validated_input_snapshot_row(
    snapshot: ValidatedInputSnapshot,
) -> ValidatedInputSnapshotRow:
    if not isinstance(snapshot, ValidatedInputSnapshot):
        raise ValueError("snapshot must be a ValidatedInputSnapshot")
    identity = snapshot.workflow_build_identity
    return ValidatedInputSnapshotRow(
        snapshot_id=snapshot.snapshot_id,
        workflow_id=snapshot.workflow_id,
        adapter_version=snapshot.adapter_version,
        schema_version=snapshot.schema_version,
        schema_dialect=snapshot.schema_dialect,
        canonical_payload=snapshot.canonical_payload,
        payload_digest_scheme=snapshot.payload_digest_scheme,
        payload_digest=snapshot.payload_digest,
        validation_outcome=snapshot.validation_outcome,
        validation_issue_codes=list(snapshot.validation_issue_codes),
        validated_at=snapshot.validated_at,
        expires_at=snapshot.expires_at,
        build_adapter_version=identity.adapter_version,
        build_scheme=identity.scheme,
        build_logical_entrypoint=identity.logical_entrypoint,
        build_digest=identity.digest,
        build_captured_at=identity.captured_at,
        consumed_run_id=snapshot.consumed_run_id,
        consumed_at=snapshot.consumed_at,
    )


def _begin_write(session: Session) -> None:
    """Serialize SQLite writers before read-then-increment sequence allocation."""
    bind = session.get_bind()
    if bind.dialect.name == "sqlite":
        session.connection().exec_driver_sql("BEGIN IMMEDIATE")


def _begin_consistent_read(session: Session) -> None:
    """Establish one SQLite snapshot before a generation-bound multi-read."""
    bind = session.get_bind()
    if bind.dialect.name == "sqlite":
        # Python sqlite3 legacy transaction control does not emit BEGIN for
        # SELECT. Without an explicit read transaction, a writer can commit
        # between the generation and row queries and produce a mixed page.
        session.connection().exec_driver_sql("BEGIN")


def _run_values(record: RunRecord) -> dict[str, object]:
    return {
        "workflow_id": record.workflow_id,
        "inputs": dict(record.inputs),
        "status": record.status.value,
        "created_at": record.created_at,
        "updated_at": record.updated_at,
        "started_at": record.started_at,
        "ended_at": record.ended_at,
        "current_stage": record.current_stage,
        "cancellation_reason": record.cancellation_reason,
        "error": _issue_to_json(record.error),
        "tags": dict(record.tags),
    }


def _record_from_row(row: RunRow) -> RunRecord:
    return RunRecord(
        run_id=row.run_id,
        workflow_id=row.workflow_id,
        inputs=row.inputs,
        status=RunStatus(row.status),
        created_at=_as_utc(row.created_at),
        updated_at=_as_utc(row.updated_at),
        started_at=_optional_utc(row.started_at),
        ended_at=_optional_utc(row.ended_at),
        current_stage=row.current_stage,
        cancellation_reason=row.cancellation_reason,
        error=_issue_from_json(row.error),
        tags=row.tags,
    )


def _run_summary_select():
    return select(
        RunRow.run_id,
        RunRow.workflow_id,
        RunRow.status,
        RunRow.created_at,
        RunRow.updated_at,
        RunRow.started_at,
        RunRow.ended_at,
        RunRow.current_stage,
    )


def _summary_from_selected_row(row: object) -> RunSummary:
    mapping = getattr(row, "_mapping", None)
    if mapping is None:
        raise ValueError("persisted run summary row is invalid")
    return RunSummary(
        run_id=mapping["run_id"],
        workflow_id=mapping["workflow_id"],
        status=mapping["status"],
        created_at=_as_utc(mapping["created_at"]),
        updated_at=_as_utc(mapping["updated_at"]),
        started_at=_optional_utc(mapping["started_at"]),
        ended_at=_optional_utc(mapping["ended_at"]),
        current_stage=mapping["current_stage"],
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


def _validated_input_snapshot_from_row(
    row: ValidatedInputSnapshotRow,
) -> ValidatedInputSnapshot:
    issue_codes = row.validation_issue_codes
    if not isinstance(issue_codes, list):
        raise ValueError("validated snapshot issue evidence is invalid")
    identity = WorkflowBuildIdentity(
        workflow_id=row.workflow_id,
        adapter_version=row.build_adapter_version,
        scheme=row.build_scheme,
        logical_entrypoint=row.build_logical_entrypoint,
        digest=row.build_digest,
        captured_at=_as_utc(row.build_captured_at),
    )
    return ValidatedInputSnapshot(
        snapshot_id=row.snapshot_id,
        workflow_id=row.workflow_id,
        adapter_version=row.adapter_version,
        schema_version=row.schema_version,
        schema_dialect=row.schema_dialect,
        workflow_build_identity=identity,
        canonical_payload=row.canonical_payload,
        payload_digest_scheme=row.payload_digest_scheme,
        payload_digest=row.payload_digest,
        validation_outcome=row.validation_outcome,
        validation_issue_codes=tuple(issue_codes),
        validated_at=_as_utc(row.validated_at),
        expires_at=_as_utc(row.expires_at),
        consumed_run_id=row.consumed_run_id,
        consumed_at=_optional_utc(row.consumed_at),
    )


def _event_from_row(row: RunEventRow) -> RunEvent:
    return RunEvent(
        event_id=row.event_id,
        run_id=row.run_id,
        sequence=row.sequence,
        event_type=row.event_type,
        timestamp=_as_utc(row.timestamp),
        status=RunStatus(row.status) if row.status is not None else None,
        stage=row.stage,
        message=row.message,
        context=row.context,
        issue=_issue_from_json(row.issue),
    )


def _log_from_row(row: RunLogRow) -> RunLogChunk:
    return RunLogChunk(
        chunk_id=row.chunk_id,
        run_id=row.run_id,
        stream_name=row.stream_name,
        sequence=row.sequence,
        timestamp=_as_utc(row.timestamp),
        lines=tuple(row.lines),
    )


def _artifact_from_row(row: RunArtifactRow) -> RunArtifactRef:
    if row.revision is None:
        raise ValueError("persisted artifact generation is unbound")
    return RunArtifactRef(
        artifact_id=row.artifact_id,
        run_id=row.run_id,
        artifact_type=row.artifact_type,
        name=row.name,
        uri=row.uri,
        mime_type=row.mime_type,
        produced_at=_as_utc(row.produced_at),
        revision=row.revision,
        metadata=row.artifact_metadata,
    )


def _artifact_row(artifact: RunArtifactRef) -> RunArtifactRow:
    return RunArtifactRow(
        artifact_id=artifact.artifact_id,
        run_id=artifact.run_id,
        artifact_type=artifact.artifact_type,
        name=artifact.name,
        uri=artifact.uri,
        mime_type=artifact.mime_type,
        produced_at=artifact.produced_at,
        revision=artifact.revision,
        artifact_metadata=artifact.to_dict()["metadata"],
    )


def _result_state_from_row(row: RunResultStateRow) -> RunResultState:
    return RunResultState(
        run_id=row.run_id,
        artifact_revision=row.artifact_revision,
        artifact_generation=row.artifact_generation,
        artifact_manifest_digest=row.artifact_manifest_digest,
        artifact_attempt_id=row.artifact_attempt_id,
        artifact_attempt_status=row.artifact_attempt_status,
        artifact_outcome=row.artifact_outcome,
        artifact_reason_code=row.artifact_reason_code,
        qc_revision=row.qc_revision,
        qc_generation=row.qc_generation,
        qc_manifest_digest=row.qc_manifest_digest,
        qc_attempt_id=row.qc_attempt_id,
        qc_attempt_status=row.qc_attempt_status,
        qc_attempt_artifact_generation=row.qc_attempt_artifact_generation,
        qc_artifact_generation=row.qc_artifact_generation,
        qc_outcome=row.qc_outcome,
        qc_reason_code=row.qc_reason_code,
    )


def _apply_result_state(row: RunResultStateRow, state: RunResultState) -> None:
    if row.run_id != state.run_id:
        raise ValueError("result state run_id does not match")
    row.artifact_revision = state.artifact_revision
    row.artifact_generation = state.artifact_generation
    row.artifact_manifest_digest = state.artifact_manifest_digest
    row.artifact_attempt_id = state.artifact_attempt_id
    row.artifact_attempt_status = state.artifact_attempt_status
    row.artifact_outcome = state.artifact_outcome
    row.artifact_reason_code = state.artifact_reason_code
    row.qc_revision = state.qc_revision
    row.qc_generation = state.qc_generation
    row.qc_manifest_digest = state.qc_manifest_digest
    row.qc_attempt_id = state.qc_attempt_id
    row.qc_attempt_status = state.qc_attempt_status
    row.qc_attempt_artifact_generation = state.qc_attempt_artifact_generation
    row.qc_artifact_generation = state.qc_artifact_generation
    row.qc_outcome = state.qc_outcome
    row.qc_reason_code = state.qc_reason_code


def _qc_metric_row(metric: RunQcMetric) -> RunQcMetricRow:
    return RunQcMetricRow(
        metric_id=metric.metric_id,
        run_id=metric.run_id,
        metric_key=metric.metric_key,
        display_name=metric.display_name,
        value_text=canonical_decimal_text(metric.value),
        unit=metric.unit,
        scope=metric.scope,
        sample_id=metric.sample_id,
        experiment_id=metric.experiment_id,
        assay=metric.assay,
        qc_flag=metric.qc_flag,
        source_artifact_id=metric.source_artifact_id,
        produced_at=metric.produced_at,
    )


def _qc_metric_from_row(row: RunQcMetricRow) -> RunQcMetric:
    metric = RunQcMetric(
        metric_id=row.metric_id,
        run_id=row.run_id,
        metric_key=row.metric_key,
        display_name=row.display_name,
        value=decimal_from_canonical_text(row.value_text),
        unit=row.unit,
        scope=row.scope,
        sample_id=row.sample_id,
        experiment_id=row.experiment_id,
        assay=row.assay,
        qc_flag=row.qc_flag,
        source_artifact_id=row.source_artifact_id,
        produced_at=_as_utc(row.produced_at),
    )
    _validate_qc_metric_fields(metric)
    return metric


def _execution_assignment_from_row(
    row: RunExecutionAssignmentRow,
) -> RunExecutionAssignment:
    return RunExecutionAssignment(
        run_id=row.run_id,
        job_id=row.job_id,
        backend=row.backend,
        queue_name=row.queue_name,
        created_at=_as_utc(row.created_at),
        dispatched_at=_optional_utc(row.dispatched_at),
        claimed_at=_optional_utc(row.claimed_at),
        cancellation_requested_at=_optional_utc(row.cancellation_requested_at),
        cancellation_reason=row.cancellation_reason,
        cancellation_acknowledged_at=_optional_utc(row.cancellation_acknowledged_at),
    )


def _workflow_build_identity_from_row(
    row: RunWorkflowBuildIdentityRow,
) -> WorkflowBuildIdentity:
    return WorkflowBuildIdentity(
        workflow_id=row.workflow_id,
        adapter_version=row.adapter_version,
        scheme=row.scheme,
        logical_entrypoint=row.logical_entrypoint,
        digest=row.digest,
        captured_at=_as_utc(row.captured_at),
    )


def _require_assignment_row_identity(
    row: RunExecutionAssignmentRow,
    *,
    job_id: str,
    backend: str,
    queue_name: str,
) -> None:
    if row.job_id != job_id or row.backend != backend or row.queue_name != queue_name:
        raise ValueError("execution assignment identity does not match")


def _assignment_row_has_ownership(
    row: RunExecutionAssignmentRow | None,
    required_ownership: RunExecutionOwnership | None,
) -> bool:
    if required_ownership is None:
        return False
    if row is None or row.backend != "rq":
        return False
    if required_ownership is RunExecutionOwnership.DISPATCHED:
        return row.dispatched_at is not None
    return row.claimed_at is not None


def _issue_to_json(issue: Issue | None) -> dict[str, object] | None:
    return issue.to_dict() if issue is not None else None


def _issue_from_json(payload: dict[str, object] | None) -> Issue | None:
    return Issue(**payload) if payload is not None else None


def _as_utc(value: datetime) -> datetime:
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _optional_utc(value: datetime | None) -> datetime | None:
    return _as_utc(value) if value is not None else None
