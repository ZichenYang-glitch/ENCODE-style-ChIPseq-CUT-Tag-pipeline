"""Workflow-platform run lifecycle primitives."""

from __future__ import annotations

from copy import deepcopy
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


@dataclass(frozen=True)
class RunRecord:
    """Immutable snapshot of a workflow run."""

    run_id: str
    workflow_id: str
    inputs: Mapping[str, Any]
    status: RunStatus
    created_at: datetime
    updated_at: datetime
    started_at: datetime | None
    ended_at: datetime | None
    current_stage: str | None
    cancellation_reason: str | None
    error: Issue | None
    tags: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "inputs", _copy_mapping(self.inputs, "inputs")
        )
        object.__setattr__(self, "tags", _copy_string_mapping(self.tags, "tags"))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready dict with fresh copies of mutable fields."""
        return {
            "run_id": self.run_id,
            "workflow_id": self.workflow_id,
            "inputs": deepcopy(dict(self.inputs)),
            "status": self.status.value,
            "created_at": self.created_at,
            "updated_at": self.updated_at,
            "started_at": self.started_at,
            "ended_at": self.ended_at,
            "current_stage": self.current_stage,
            "cancellation_reason": self.cancellation_reason,
            "error": self.error.to_dict() if self.error is not None else None,
            "tags": dict(self.tags),
        }


def _copy_mapping(mapping: Mapping[str, object], name: str) -> dict[str, object]:
    if not isinstance(mapping, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return deepcopy(dict(mapping))


def _copy_string_mapping(mapping: Mapping[str, str], name: str) -> dict[str, str]:
    if not isinstance(mapping, Mapping):
        raise ValueError(f"{name} must be a mapping")
    copied = dict(mapping)
    for key, value in copied.items():
        if not isinstance(key, str) or not isinstance(value, str):
            raise ValueError(f"{name} keys and values must be strings")
    return copied
