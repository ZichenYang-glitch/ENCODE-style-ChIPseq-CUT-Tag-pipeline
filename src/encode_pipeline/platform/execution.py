"""Durable execution-assignment primitives for workflow runs."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from hashlib import sha256

from encode_pipeline.platform.runs import RunRecord


class RunExecutionOwnership(str, Enum):
    """Durable hand-off evidence required for an active lifecycle state."""

    DISPATCHED = "dispatched"
    CLAIMED = "claimed"


@dataclass(frozen=True)
class RunExecutionAssignment:
    """Stable identity linking one persisted run to an execution backend job."""

    run_id: str
    job_id: str
    backend: str
    queue_name: str
    created_at: datetime
    dispatched_at: datetime | None = None
    claimed_at: datetime | None = None
    cancellation_requested_at: datetime | None = None
    cancellation_reason: str | None = None
    cancellation_acknowledged_at: datetime | None = None

    def __post_init__(self) -> None:
        for field_name in ("run_id", "job_id", "backend", "queue_name"):
            value = getattr(self, field_name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{field_name} must be a non-empty string")
        if not isinstance(self.created_at, datetime):
            raise ValueError("created_at must be a datetime")
        for field_name in (
            "dispatched_at",
            "claimed_at",
            "cancellation_requested_at",
            "cancellation_acknowledged_at",
        ):
            value = getattr(self, field_name)
            if value is not None and not isinstance(value, datetime):
                raise ValueError(f"{field_name} must be a datetime or None")
        if self.claimed_at is not None and self.dispatched_at is None:
            raise ValueError("claimed_at requires dispatched_at")
        if (self.cancellation_requested_at is None) != (
            self.cancellation_reason is None
        ):
            raise ValueError(
                "cancellation request and reason must be recorded together"
            )
        if self.cancellation_reason is not None and (
            not isinstance(self.cancellation_reason, str)
            or not self.cancellation_reason.strip()
        ):
            raise ValueError("cancellation_reason must be a non-empty string")
        if self.cancellation_requested_at is not None and self.claimed_at is None:
            raise ValueError("cancellation request requires claimed_at")
        if (
            self.cancellation_acknowledged_at is not None
            and self.cancellation_requested_at is None
        ):
            raise ValueError("cancellation acknowledgement requires a request")


@dataclass(frozen=True)
class RunExecutionClaim:
    """Result of an atomic attempt to claim one durable execution job."""

    assignment: RunExecutionAssignment
    acquired: bool

    def __post_init__(self) -> None:
        if not isinstance(self.assignment, RunExecutionAssignment):
            raise ValueError("assignment must be a RunExecutionAssignment")
        if not isinstance(self.acquired, bool):
            raise ValueError("acquired must be a bool")


@dataclass(frozen=True)
class RunExecutionCancellationRequest:
    """Canonical result of atomically recording cancellation intent."""

    assignment: RunExecutionAssignment
    record: RunRecord
    created: bool

    def __post_init__(self) -> None:
        if not isinstance(self.assignment, RunExecutionAssignment):
            raise ValueError("assignment must be a RunExecutionAssignment")
        if not isinstance(self.record, RunRecord):
            raise ValueError("record must be a RunRecord")
        if not isinstance(self.created, bool):
            raise ValueError("created must be a bool")


@dataclass(frozen=True)
class RunExecutionStopAcknowledgement:
    """Canonical result after an RQ work horse has been reaped."""

    assignment: RunExecutionAssignment
    record: RunRecord
    transitioned: bool

    def __post_init__(self) -> None:
        if not isinstance(self.assignment, RunExecutionAssignment):
            raise ValueError("assignment must be a RunExecutionAssignment")
        if not isinstance(self.record, RunRecord):
            raise ValueError("record must be a RunRecord")
        if not isinstance(self.transitioned, bool):
            raise ValueError("transitioned must be a bool")


def build_execution_job_id(run_id: str) -> str:
    """Return a deterministic RQ-safe job identity for ``run_id``."""
    if not isinstance(run_id, str) or not run_id.strip():
        raise ValueError("run_id must be a non-empty string")
    digest = sha256(run_id.encode("utf-8")).hexdigest()
    return f"run-execution-{digest}"
