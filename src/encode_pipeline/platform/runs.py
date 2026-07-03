"""Workflow-platform run lifecycle primitives."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Any, Mapping

from encode_pipeline.platform.results import Issue


class RunStatus(str, Enum):
    """Finite lifecycle statuses for a workflow run."""

    CREATED = "created"
    VALIDATING = "validating"
    PLANNED = "planned"
    QUEUED = "queued"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    CANCELLED = "cancelled"

    @property
    def is_terminal(self) -> bool:
        """Return True for terminal, absorbing states."""
        return self in {
            RunStatus.SUCCEEDED,
            RunStatus.FAILED,
            RunStatus.CANCELLED,
        }


_VALID_TRANSITIONS: dict[RunStatus, set[RunStatus]] = {
    RunStatus.CREATED: {RunStatus.VALIDATING, RunStatus.CANCELLED},
    RunStatus.VALIDATING: {RunStatus.PLANNED, RunStatus.FAILED, RunStatus.CANCELLED},
    RunStatus.PLANNED: {RunStatus.QUEUED, RunStatus.CANCELLED},
    RunStatus.QUEUED: {RunStatus.RUNNING, RunStatus.CANCELLED},
    RunStatus.RUNNING: {RunStatus.SUCCEEDED, RunStatus.FAILED, RunStatus.CANCELLED},
    RunStatus.SUCCEEDED: set(),
    RunStatus.FAILED: set(),
    RunStatus.CANCELLED: set(),
}


def can_transition(from_status: RunStatus, to_status: RunStatus) -> bool:
    """Return True if ``to_status`` is a legal next state."""
    return to_status in _VALID_TRANSITIONS.get(from_status, set())


def require_transition(from_status: RunStatus, to_status: RunStatus) -> None:
    """Raise ValueError if the transition is illegal; otherwise return None."""
    if not can_transition(from_status, to_status):
        raise ValueError(
            f"Illegal transition from {from_status.value!r} to {to_status.value!r}."
        )
