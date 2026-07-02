from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import Protocol


@dataclass(frozen=True)
class AuditEvent:
    event_type: str
    workflow_id: str | None
    session_id: str | None
    issue_count: int
    tool_name: str | None
    filtered: bool
    redaction_counts: Mapping[str, int]
    issue_codes: tuple[str, ...]
    issue_sources: tuple[str, ...]
    issue_severities: tuple[str, ...]


class AgentAuditSink(Protocol):
    def record(self, event: AuditEvent) -> None: ...
    def events(self) -> Sequence[AuditEvent]: ...


class NoOpAuditSink:
    """Audit sink that records nothing.

    Use when AgentService must satisfy the ``AgentAuditSink`` protocol
    without retaining events.
    """

    def record(self, event: AuditEvent) -> None:
        return None

    def events(self) -> Sequence[AuditEvent]:
        return ()


class InMemoryAuditSink:
    """Bounded in-memory audit sink."""

    def __init__(self, bound: int = 1000) -> None:
        if bound <= 0:
            raise ValueError("bound must be positive")
        self._bound = bound
        self._events: list[AuditEvent] = []

    def record(self, event: AuditEvent) -> None:
        self._events.append(event)
        if len(self._events) > self._bound:
            self._events.pop(0)

    def events(self) -> Sequence[AuditEvent]:
        return tuple(self._events)
