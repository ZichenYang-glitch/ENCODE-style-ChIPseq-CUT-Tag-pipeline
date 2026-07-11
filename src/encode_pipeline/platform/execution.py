"""Durable execution-assignment primitives for workflow runs."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from hashlib import sha256


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

    def __post_init__(self) -> None:
        for field_name in ("run_id", "job_id", "backend", "queue_name"):
            value = getattr(self, field_name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{field_name} must be a non-empty string")
        if not isinstance(self.created_at, datetime):
            raise ValueError("created_at must be a datetime")
        for field_name in ("dispatched_at", "claimed_at"):
            value = getattr(self, field_name)
            if value is not None and not isinstance(value, datetime):
                raise ValueError(f"{field_name} must be a datetime or None")
        if self.claimed_at is not None and self.dispatched_at is None:
            raise ValueError("claimed_at requires dispatched_at")


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


def build_execution_job_id(run_id: str) -> str:
    """Return a deterministic RQ-safe job identity for ``run_id``."""
    if not isinstance(run_id, str) or not run_id.strip():
        raise ValueError("run_id must be a non-empty string")
    digest = sha256(run_id.encode("utf-8")).hexdigest()
    return f"run-execution-{digest}"
