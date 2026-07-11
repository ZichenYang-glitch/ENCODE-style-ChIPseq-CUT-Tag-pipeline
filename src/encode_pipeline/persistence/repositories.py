"""SQLAlchemy implementation of the workflow run repository boundary."""

from __future__ import annotations

from collections.abc import Iterable
from datetime import datetime, timezone
from threading import RLock

from sqlalchemy import func, select, update
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session, sessionmaker

from encode_pipeline.persistence.models import (
    RunArtifactRow,
    RunEventRow,
    RunLogRow,
    RunRow,
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

    def add_event(self, run_id: str, event: RunEventDraft) -> RunEvent:
        with self._lock, self._session_factory.begin() as session:
            self._require_run(session, run_id)
            return self._insert_event(session, run_id, event)

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

    def list_artifacts(self, run_id: str) -> tuple[RunArtifactRef, ...]:
        with self._session_factory() as session:
            self._require_run(session, run_id)
            rows = session.scalars(
                select(RunArtifactRow)
                .where(RunArtifactRow.run_id == run_id)
                .order_by(RunArtifactRow.id)
            ).all()
            return tuple(_artifact_from_row(row) for row in rows)

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
