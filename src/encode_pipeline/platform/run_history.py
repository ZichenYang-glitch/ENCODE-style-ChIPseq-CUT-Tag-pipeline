"""Workflow-neutral, disclosure-safe run-history primitives."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import base64
import binascii
import json
import re

from encode_pipeline.platform.runs import RunStatus


RUN_HISTORY_CURSOR_PREFIX = "runhist_"
RUN_HISTORY_CURSOR_MAX_LENGTH = 1024
_CURSOR_BODY_PATTERN = re.compile(r"^[A-Za-z0-9_-]+$")
_RUN_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_WORKFLOW_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,254}$")
_STAGE_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,254}$")
_CURSOR_KEYS = frozenset({"v", "created_at", "run_id", "workflow_id", "status"})


def normalize_run_history_workflow_filter(value: str | None) -> str | None:
    """Normalize a public workflow filter or reject it as unsafe."""
    if value is None:
        return None
    if not isinstance(value, str):
        raise ValueError("workflow_id filter must be a string")
    normalized = value.strip()
    if _WORKFLOW_ID_PATTERN.fullmatch(normalized) is None:
        raise ValueError("workflow_id filter is invalid")
    return normalized


def normalize_run_history_status(value: RunStatus | str | None) -> RunStatus | None:
    """Normalize an optional exact lifecycle filter."""
    if value is None:
        return None
    if isinstance(value, RunStatus):
        return value
    if isinstance(value, str):
        try:
            return RunStatus(value)
        except ValueError as exc:
            raise ValueError("run status filter is invalid") from exc
    raise ValueError("run status filter is invalid")


@dataclass(frozen=True)
class RunSummary:
    """Minimal persisted run projection safe for history consumers."""

    run_id: str
    workflow_id: str
    status: RunStatus | str
    created_at: datetime
    updated_at: datetime
    started_at: datetime | None
    ended_at: datetime | None
    current_stage: str | None

    def __post_init__(self) -> None:
        if (
            not isinstance(self.run_id, str)
            or _RUN_ID_PATTERN.fullmatch(self.run_id) is None
        ):
            raise ValueError("run summary run_id is not public-safe")
        if (
            not isinstance(self.workflow_id, str)
            or _WORKFLOW_ID_PATTERN.fullmatch(self.workflow_id) is None
        ):
            raise ValueError("run summary workflow_id is not public-safe")
        object.__setattr__(self, "status", _required_status(self.status))
        object.__setattr__(
            self, "created_at", _required_utc(self.created_at, "created_at")
        )
        object.__setattr__(
            self, "updated_at", _required_utc(self.updated_at, "updated_at")
        )
        object.__setattr__(
            self, "started_at", _optional_utc(self.started_at, "started_at")
        )
        object.__setattr__(self, "ended_at", _optional_utc(self.ended_at, "ended_at"))
        if self.current_stage is not None and (
            not isinstance(self.current_stage, str)
            or _STAGE_PATTERN.fullmatch(self.current_stage) is None
        ):
            raise ValueError("run summary current_stage is not public-safe")
        _validate_timestamp_order(self)
        _validate_lifecycle_evidence(self)

    def to_dict(self) -> dict[str, object]:
        """Return only the public run-history fields."""
        return {
            "run_id": self.run_id,
            "workflow_id": self.workflow_id,
            "status": self.status.value,
            "created_at": self.created_at,
            "updated_at": self.updated_at,
            "started_at": self.started_at,
            "ended_at": self.ended_at,
            "current_stage": self.current_stage,
        }


@dataclass(frozen=True)
class RunHistoryCursor:
    """Validated keyset boundary and the exact filters it belongs to."""

    created_at: datetime
    run_id: str
    workflow_id: str | None
    status: RunStatus | str | None

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "created_at", _required_utc(self.created_at, "created_at")
        )
        if (
            not isinstance(self.run_id, str)
            or _RUN_ID_PATTERN.fullmatch(self.run_id) is None
        ):
            raise ValueError("run history cursor run_id is invalid")
        object.__setattr__(
            self,
            "workflow_id",
            normalize_run_history_workflow_filter(self.workflow_id),
        )
        object.__setattr__(self, "status", normalize_run_history_status(self.status))


@dataclass(frozen=True)
class RunHistoryPage:
    """One bounded run-history page and its optional next cursor."""

    runs: tuple[RunSummary, ...]
    next_cursor: str | None

    def __post_init__(self) -> None:
        if not isinstance(self.runs, tuple) or any(
            not isinstance(run, RunSummary) for run in self.runs
        ):
            raise ValueError("run history page runs must be RunSummary values")
        if self.next_cursor is not None and not isinstance(self.next_cursor, str):
            raise ValueError("run history next_cursor must be a string or None")


def encode_run_history_cursor(cursor: RunHistoryCursor) -> str:
    """Encode canonical protocol state; the result is not an auth credential."""
    if not isinstance(cursor, RunHistoryCursor):
        raise ValueError("run history cursor is invalid")
    payload = {
        "created_at": _canonical_datetime(cursor.created_at),
        "run_id": cursor.run_id,
        "status": cursor.status.value if cursor.status is not None else None,
        "v": 1,
        "workflow_id": cursor.workflow_id,
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    body = base64.urlsafe_b64encode(raw).decode("ascii").rstrip("=")
    encoded = f"{RUN_HISTORY_CURSOR_PREFIX}{body}"
    if len(encoded) > RUN_HISTORY_CURSOR_MAX_LENGTH:
        raise ValueError("run history cursor exceeds its public bound")
    return encoded


def decode_run_history_cursor(
    value: str,
    *,
    workflow_id: str | None,
    status: RunStatus | str | None,
) -> RunHistoryCursor:
    """Strictly decode a canonical cursor and bind it to current filters."""
    normalized_workflow_id = normalize_run_history_workflow_filter(workflow_id)
    normalized_status = normalize_run_history_status(status)
    if (
        not isinstance(value, str)
        or len(value) <= len(RUN_HISTORY_CURSOR_PREFIX)
        or len(value) > RUN_HISTORY_CURSOR_MAX_LENGTH
        or not value.startswith(RUN_HISTORY_CURSOR_PREFIX)
    ):
        raise ValueError("run history cursor is invalid")
    body = value[len(RUN_HISTORY_CURSOR_PREFIX) :]
    if _CURSOR_BODY_PATTERN.fullmatch(body) is None:
        raise ValueError("run history cursor is invalid")
    padding = "=" * (-len(body) % 4)
    try:
        raw = base64.b64decode(
            body + padding,
            altchars=b"-_",
            validate=True,
        )
        payload = json.loads(raw.decode("utf-8"))
    except (binascii.Error, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("run history cursor is invalid") from exc
    if not isinstance(payload, dict) or set(payload) != _CURSOR_KEYS:
        raise ValueError("run history cursor payload is invalid")
    if payload["v"] != 1 or isinstance(payload["v"], bool):
        raise ValueError("run history cursor version is invalid")
    if not isinstance(payload["created_at"], str):
        raise ValueError("run history cursor timestamp is invalid")
    try:
        created_at = datetime.fromisoformat(
            payload["created_at"].replace("Z", "+00:00")
        )
    except ValueError as exc:
        raise ValueError("run history cursor timestamp is invalid") from exc
    cursor = RunHistoryCursor(
        created_at=created_at,
        run_id=payload["run_id"],
        workflow_id=payload["workflow_id"],
        status=payload["status"],
    )
    if encode_run_history_cursor(cursor) != value:
        raise ValueError("run history cursor is not canonical")
    if (
        cursor.workflow_id != normalized_workflow_id
        or cursor.status is not normalized_status
    ):
        raise ValueError("run history cursor filters do not match")
    return cursor


def _required_status(value: RunStatus | str) -> RunStatus:
    normalized = normalize_run_history_status(value)
    if normalized is None:  # pragma: no cover - type guard
        raise ValueError("run summary status is required")
    return normalized


def _required_utc(value: datetime, name: str) -> datetime:
    if (
        not isinstance(value, datetime)
        or value.tzinfo is None
        or value.utcoffset() is None
    ):
        raise ValueError(f"run history {name} must be timezone-aware")
    return value.astimezone(timezone.utc)


def _optional_utc(value: datetime | None, name: str) -> datetime | None:
    return None if value is None else _required_utc(value, name)


def _canonical_datetime(value: datetime) -> str:
    return (
        _required_utc(value, "created_at")
        .isoformat(timespec="microseconds")
        .replace("+00:00", "Z")
    )


def _validate_timestamp_order(summary: RunSummary) -> None:
    if summary.updated_at < summary.created_at:
        raise ValueError("run summary timestamp order is invalid")
    for value in (summary.started_at, summary.ended_at):
        if value is not None and (
            value < summary.created_at or value > summary.updated_at
        ):
            raise ValueError("run summary timestamp order is invalid")


def _validate_lifecycle_evidence(summary: RunSummary) -> None:
    terminal_evidence = summary.ended_at is not None
    if summary.status.is_terminal is not terminal_evidence:
        raise ValueError("run summary lifecycle evidence is invalid")
    if (
        summary.status
        in {
            RunStatus.CREATED,
            RunStatus.VALIDATING,
            RunStatus.PLANNED,
            RunStatus.QUEUED,
        }
        and summary.started_at is not None
    ):
        raise ValueError("run summary lifecycle evidence is invalid")
    if summary.status in {RunStatus.RUNNING, RunStatus.SUCCEEDED} and (
        summary.started_at is None
    ):
        raise ValueError("run summary lifecycle evidence is invalid")
    if (
        summary.started_at is not None
        and summary.ended_at is not None
        and summary.ended_at < summary.started_at
    ):
        raise ValueError("run summary timestamp order is invalid")
