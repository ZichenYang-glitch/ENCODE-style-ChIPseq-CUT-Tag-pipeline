"""SQLAlchemy implementation of the workflow run repository boundary."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import replace
from datetime import datetime, timezone
from threading import RLock

from sqlalchemy import delete, func, select, update
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session, sessionmaker

from encode_pipeline.persistence.models import (
    RunArtifactRow,
    RunEventRow,
    RunExecutionAssignmentRow,
    RunLogRow,
    RunRow,
    RunWorkflowBuildIdentityRow,
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
from encode_pipeline.platform.runs import (
    RunArtifactRef,
    RunEvent,
    RunLogChunk,
    RunRecord,
    RunStatus,
)
from encode_pipeline.services.run_repositories import (
    ConcurrentRunUpdateError,
    RunEventDraft,
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
                    return self._insert_event(session, record.run_id, event)
            except IntegrityError as exc:
                raise ValueError(f"Duplicate run_id: {record.run_id!r}") from exc

    def get_run(self, run_id: str) -> RunRecord:
        with self._session_factory() as session:
            return _record_from_row(self._require_run(session, run_id))

    def list_runs(self) -> tuple[RunRecord, ...]:
        with self._session_factory() as session:
            rows = session.scalars(select(RunRow).order_by(RunRow.id)).all()
            return tuple(_record_from_row(row) for row in rows)

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

    def record_artifact(
        self,
        run_id: str,
        artifact: RunArtifactRef,
    ) -> RunArtifactRef:
        with self._lock:
            try:
                with self._session_factory.begin() as session:
                    _begin_write(session)
                    self._require_run(session, run_id)
                    session.add(
                        RunArtifactRow(
                            artifact_id=artifact.artifact_id,
                            run_id=run_id,
                            artifact_type=artifact.artifact_type,
                            name=artifact.name,
                            uri=artifact.uri,
                            mime_type=artifact.mime_type,
                            produced_at=artifact.produced_at,
                            artifact_metadata=dict(artifact.metadata),
                        )
                    )
                    session.flush()
            except IntegrityError as exc:
                raise ValueError(
                    f"Duplicate artifact_id: {artifact.artifact_id!r}"
                ) from exc
            return artifact

    def replace_artifacts(
        self,
        run_id: str,
        artifacts: tuple[RunArtifactRef, ...],
        *,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent | None:
        try:
            with self._lock, self._session_factory.begin() as session:
                _begin_write(session)
                current = self._require_run(session, run_id)
                if RunStatus(current.status) is not expected_status:
                    raise ConcurrentRunUpdateError(
                        f"Run {run_id!r} is no longer eligible for artifact replacement."
                    )
                replacement_ids: set[str] = set()
                for artifact in artifacts:
                    if artifact.run_id != run_id:
                        raise ValueError("artifact run_id does not match the run")
                    if artifact.artifact_id in replacement_ids:
                        raise ValueError("duplicate artifact_id in replacement")
                    replacement_ids.add(artifact.artifact_id)
                existing_rows = session.scalars(
                    select(RunArtifactRow)
                    .where(RunArtifactRow.run_id == run_id)
                    .order_by(RunArtifactRow.id)
                ).all()
                existing = tuple(_artifact_from_row(row) for row in existing_rows)
                latest_outcome = session.scalar(
                    select(RunEventRow.event_type)
                    .where(
                        RunEventRow.run_id == run_id,
                        RunEventRow.event_type.in_(
                            ("artifacts_indexed", "artifact_extraction_failed")
                        ),
                    )
                    .order_by(RunEventRow.id.desc())
                    .limit(1)
                )
                if existing == artifacts and latest_outcome == "artifacts_indexed":
                    return None
                session.execute(
                    delete(RunArtifactRow).where(RunArtifactRow.run_id == run_id)
                )
                session.add_all(
                    [
                        RunArtifactRow(
                            artifact_id=artifact.artifact_id,
                            run_id=run_id,
                            artifact_type=artifact.artifact_type,
                            name=artifact.name,
                            uri=artifact.uri,
                            mime_type=artifact.mime_type,
                            produced_at=artifact.produced_at,
                            artifact_metadata=dict(artifact.metadata),
                        )
                        for artifact in artifacts
                    ]
                )
                session.flush()
                return self._insert_event(session, run_id, event)
        except IntegrityError as exc:
            raise ValueError("artifact replacement could not be persisted") from exc

    def list_artifacts(self, run_id: str) -> tuple[RunArtifactRef, ...]:
        with self._session_factory() as session:
            self._require_run(session, run_id)
            rows = session.scalars(
                select(RunArtifactRow)
                .where(RunArtifactRow.run_id == run_id)
                .order_by(RunArtifactRow.id)
            ).all()
            return tuple(_artifact_from_row(row) for row in rows)

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


def _begin_write(session: Session) -> None:
    """Serialize SQLite writers before read-then-increment sequence allocation."""
    bind = session.get_bind()
    if bind.dialect.name == "sqlite":
        session.connection().exec_driver_sql("BEGIN IMMEDIATE")


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
    return RunArtifactRef(
        artifact_id=row.artifact_id,
        run_id=row.run_id,
        artifact_type=row.artifact_type,
        name=row.name,
        uri=row.uri,
        mime_type=row.mime_type,
        produced_at=_as_utc(row.produced_at),
        metadata=row.artifact_metadata,
    )


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
