"""Workflow-neutral run-history read-model and cursor contract tests."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timedelta, timezone

import pytest

from encode_pipeline.platform.run_history import (
    RunHistoryCursor,
    RunSummary,
    decode_run_history_cursor,
    encode_run_history_cursor,
)
from encode_pipeline.platform.runs import RunStatus


NOW = datetime(2026, 7, 14, 8, 0, tzinfo=timezone.utc)


def _summary(**changes) -> RunSummary:
    values = {
        "run_id": "run-1",
        "workflow_id": "workflow-1",
        "status": RunStatus.SUCCEEDED,
        "created_at": NOW,
        "updated_at": NOW + timedelta(minutes=2),
        "started_at": NOW + timedelta(seconds=5),
        "ended_at": NOW + timedelta(minutes=2),
        "current_stage": "qc_summary_indexing",
    }
    values.update(changes)
    return RunSummary(**values)


def test_run_summary_normalizes_status_and_timestamps() -> None:
    summary = _summary(status="succeeded")

    assert summary.status is RunStatus.SUCCEEDED
    assert summary.created_at.tzinfo is timezone.utc


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("run_id", "/private/run"),
        ("workflow_id", "WORKFLOW_TOKEN=/private"),
        ("current_stage", "execution\nsecret"),
    ],
)
def test_run_summary_rejects_unsafe_public_text(field: str, value: str) -> None:
    with pytest.raises(ValueError):
        _summary(**{field: value})


def test_run_summary_rejects_inconsistent_timestamp_order() -> None:
    with pytest.raises(ValueError, match="timestamp"):
        _summary(updated_at=NOW - timedelta(seconds=1))


@pytest.mark.parametrize(
    ("status", "started_at", "ended_at"),
    [
        (RunStatus.SUCCEEDED, NOW + timedelta(seconds=5), None),
        (RunStatus.RUNNING, None, None),
        (RunStatus.PLANNED, NOW + timedelta(seconds=5), None),
        (
            RunStatus.RUNNING,
            NOW + timedelta(seconds=5),
            NOW + timedelta(minutes=2),
        ),
    ],
)
def test_run_summary_rejects_status_without_consistent_lifecycle_evidence(
    status: RunStatus,
    started_at: datetime | None,
    ended_at: datetime | None,
) -> None:
    with pytest.raises(ValueError, match="lifecycle"):
        _summary(status=status, started_at=started_at, ended_at=ended_at)


def test_run_history_cursor_round_trips_canonically_and_binds_filters() -> None:
    cursor = RunHistoryCursor(
        created_at=NOW,
        run_id="run-1",
        workflow_id="workflow-1",
        status=RunStatus.SUCCEEDED,
    )

    encoded = encode_run_history_cursor(cursor)

    assert encoded.startswith("runhist_")
    assert (
        decode_run_history_cursor(
            encoded,
            workflow_id="workflow-1",
            status=RunStatus.SUCCEEDED,
        )
        == cursor
    )
    with pytest.raises(ValueError, match="filters"):
        decode_run_history_cursor(
            encoded,
            workflow_id="workflow-2",
            status=RunStatus.SUCCEEDED,
        )


@pytest.mark.parametrize(
    "value",
    ["", "not-a-cursor", "runhist_!", "runhist_" + "a" * 1025],
)
def test_run_history_cursor_rejects_malformed_or_unbounded_tokens(value: str) -> None:
    with pytest.raises(ValueError):
        decode_run_history_cursor(value, workflow_id=None, status=None)


def test_run_history_cursor_rejects_noncanonical_payload() -> None:
    encoded = encode_run_history_cursor(
        RunHistoryCursor(
            created_at=NOW,
            run_id="run-1",
            workflow_id=None,
            status=None,
        )
    )
    mutated = replace(
        RunHistoryCursor(
            created_at=NOW,
            run_id="run-1",
            workflow_id=None,
            status=None,
        ),
        run_id="run-2",
    )

    assert encoded != encode_run_history_cursor(mutated)
