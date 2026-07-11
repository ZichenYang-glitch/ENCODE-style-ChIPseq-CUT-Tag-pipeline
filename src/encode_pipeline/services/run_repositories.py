"""Storage boundary for workflow run aggregates."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from datetime import datetime, timezone
from threading import RLock
from typing import Any, Protocol, TypeVar

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

    def list_artifacts(self, run_id: str) -> tuple[RunArtifactRef, ...]: ...


class InMemoryRunRepository:
    """Thread-safe in-memory implementation used by unit tests and adapters."""

    def __init__(self) -> None:
        self._lock = RLock()
        self._runs: dict[str, RunRecord] = {}
        self._events: dict[str, list[RunEvent]] = {}
        self._logs: dict[str, dict[str, list[RunLogChunk]]] = {}
        self._artifacts: dict[str, dict[str, RunArtifactRef]] = {}

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

    def list_artifacts(self, run_id: str) -> tuple[RunArtifactRef, ...]:
        with self._lock:
            return tuple(self._artifacts[run_id].values())

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
