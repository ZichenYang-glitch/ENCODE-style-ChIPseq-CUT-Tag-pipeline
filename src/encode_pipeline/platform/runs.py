"""Workflow-platform run lifecycle primitives."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from datetime import datetime
from decimal import Decimal
from enum import Enum
from hashlib import sha256
import re
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
    RunStatus.QUEUED: {
        RunStatus.RUNNING,
        RunStatus.FAILED,
        RunStatus.CANCELLED,
    },
    # RUNNING -> CANCELLED is guarded by the execution-stop acknowledgement
    # repository operation; RunService.transition_run cannot take this edge.
    RunStatus.RUNNING: {
        RunStatus.SUCCEEDED,
        RunStatus.FAILED,
        RunStatus.CANCELLED,
    },
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


def _normalize_status(value: RunStatus | str) -> RunStatus:
    """Return a RunStatus, accepting either a RunStatus or its string value."""
    if isinstance(value, RunStatus):
        return value
    if isinstance(value, str):
        try:
            return RunStatus(value)
        except ValueError as exc:
            raise ValueError(f"Invalid run status: {value!r}") from exc
    raise ValueError(f"Invalid run status: {value!r}")


@dataclass(frozen=True)
class RunRecord:
    """Immutable snapshot of a workflow run."""

    run_id: str
    workflow_id: str
    inputs: Mapping[str, Any]
    status: RunStatus | str
    created_at: datetime
    updated_at: datetime
    started_at: datetime | None
    ended_at: datetime | None
    current_stage: str | None
    cancellation_reason: str | None
    error: Issue | None
    tags: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "status", _normalize_status(self.status))
        object.__setattr__(self, "inputs", _copy_mapping(self.inputs, "inputs"))
        object.__setattr__(self, "tags", _copy_string_mapping(self.tags, "tags"))
        if self.error is not None and not isinstance(self.error, Issue):
            raise ValueError("RunRecord error must be an Issue or None")

    def to_dict(self) -> dict[str, Any]:
        """Return a dict representation with fresh copies of mutable fields."""
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


@dataclass(frozen=True)
class RunEvent:
    """Immutable record of a run lifecycle event."""

    event_id: str
    run_id: str
    sequence: int
    event_type: str
    timestamp: datetime
    status: RunStatus | None
    stage: str | None
    message: str
    context: Mapping[str, Any] = field(default_factory=dict)
    issue: Issue | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "context", _copy_mapping(self.context, "context"))
        if self.status is not None and not isinstance(self.status, RunStatus):
            raise ValueError("RunEvent status must be a RunStatus or None")
        if not isinstance(self.sequence, int) or isinstance(self.sequence, bool):
            raise ValueError("RunEvent sequence must be an integer")
        if self.issue is not None and not isinstance(self.issue, Issue):
            raise ValueError("RunEvent issue must be an Issue or None")

    def to_dict(self) -> dict[str, Any]:
        """Return a dict representation with fresh copies of mutable fields."""
        return {
            "event_id": self.event_id,
            "run_id": self.run_id,
            "sequence": self.sequence,
            "event_type": self.event_type,
            "timestamp": self.timestamp,
            "status": self.status.value if self.status is not None else None,
            "stage": self.stage,
            "message": self.message,
            "context": deepcopy(dict(self.context)),
            "issue": self.issue.to_dict() if self.issue is not None else None,
        }


@dataclass(frozen=True)
class RunLogChunk:
    """Immutable append-only chunk of a run log stream."""

    chunk_id: str
    run_id: str
    stream_name: str
    sequence: int
    timestamp: datetime
    lines: tuple[str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "lines", _normalize_string_tuple(self.lines, "lines"))
        if not isinstance(self.sequence, int) or isinstance(self.sequence, bool):
            raise ValueError("RunLogChunk sequence must be an integer")

    def to_dict(self) -> dict[str, Any]:
        """Return a dict representation."""
        return {
            "chunk_id": self.chunk_id,
            "run_id": self.run_id,
            "stream_name": self.stream_name,
            "sequence": self.sequence,
            "timestamp": self.timestamp,
            "lines": list(self.lines),
        }


def _normalize_string_tuple(values: Any, name: str) -> tuple[str, ...]:
    if isinstance(values, str):
        raw_values = (values,)
    else:
        try:
            raw_values = tuple(values)
        except TypeError as exc:
            raise ValueError(f"{name} must be an iterable of strings") from exc

    normalized = []
    for value in raw_values:
        if not isinstance(value, str):
            raise ValueError(f"{name} entries must be strings")
        normalized.append(value)
    return tuple(normalized)


@dataclass(frozen=True)
class RunArtifactRef:
    """Opaque reference to an artifact produced by a run."""

    artifact_id: str
    run_id: str
    artifact_type: str
    name: str
    uri: str
    mime_type: str | None
    produced_at: datetime
    metadata: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "metadata", _copy_mapping(self.metadata, "metadata"))

    def to_dict(self) -> dict[str, Any]:
        """Return a dict representation with fresh copies of mutable fields."""
        return {
            "artifact_id": self.artifact_id,
            "run_id": self.run_id,
            "artifact_type": self.artifact_type,
            "name": self.name,
            "uri": self.uri,
            "mime_type": self.mime_type,
            "produced_at": self.produced_at,
            "metadata": deepcopy(dict(self.metadata)),
        }


_QC_IDENTIFIER_TOKEN = re.compile(r"^[A-Za-z0-9_.-]{1,255}$")


def build_qc_metric_id(
    metric_key: str,
    scope: str,
    sample_id: str | None,
    experiment_id: str | None,
) -> str:
    """Build the stable identifier for one semantic QC metric coordinate."""
    if (
        not isinstance(metric_key, str)
        or not isinstance(scope, str)
        or (sample_id is not None and not isinstance(sample_id, str))
        or (experiment_id is not None and not isinstance(experiment_id, str))
    ):
        raise ValueError("QC metric identity coordinates must be strings")
    coordinates = (
        metric_key,
        scope,
        "" if sample_id is None else sample_id,
        "" if experiment_id is None else experiment_id,
    )
    digest = sha256()
    for value in coordinates:
        encoded = value.encode()
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return f"qcmetric-{digest.hexdigest()}"


def validate_qc_identifier_token(value: object) -> str:
    """Return a canonical ENCODE-safe identifier or reject it."""
    if (
        not isinstance(value, str)
        or value in {".", ".."}
        or _QC_IDENTIFIER_TOKEN.fullmatch(value) is None
    ):
        raise ValueError("QC identifier token is invalid")
    return value


@dataclass(frozen=True)
class RunQcMetric:
    """One durable workflow-neutral numeric QC metric."""

    metric_id: str
    run_id: str
    metric_key: str
    display_name: str
    value: Decimal
    unit: str
    scope: str
    source_artifact_id: str
    produced_at: datetime
    sample_id: str | None = None
    experiment_id: str | None = None
    assay: str | None = None
    qc_flag: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.value, Decimal):
            raise ValueError("value must be a Decimal")

    def to_dict(self) -> dict[str, Any]:
        """Return the explicit durable QC fields without an open metadata bag."""
        return {
            "metric_id": self.metric_id,
            "run_id": self.run_id,
            "metric_key": self.metric_key,
            "display_name": self.display_name,
            "value": self.value,
            "unit": self.unit,
            "scope": self.scope,
            "sample_id": self.sample_id,
            "experiment_id": self.experiment_id,
            "assay": self.assay,
            "qc_flag": self.qc_flag,
            "source_artifact_id": self.source_artifact_id,
            "produced_at": self.produced_at,
        }
