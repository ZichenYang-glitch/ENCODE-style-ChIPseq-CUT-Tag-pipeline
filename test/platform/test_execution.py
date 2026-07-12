"""Tests for durable workflow execution identities."""

from __future__ import annotations

from dataclasses import FrozenInstanceError
from datetime import datetime, timezone
import re

import pytest

from encode_pipeline.platform.execution import (
    RunExecutionAssignment,
    RunExecutionClaim,
    build_execution_job_id,
)


def test_execution_assignment_is_an_immutable_value_object():
    created_at = datetime.now(timezone.utc)
    assignment = RunExecutionAssignment(
        run_id="run-1",
        job_id="job-1",
        backend="rq",
        queue_name="default",
        created_at=created_at,
    )

    assert assignment.created_at is created_at
    assert assignment.dispatched_at is None
    assert assignment.claimed_at is None
    with pytest.raises(FrozenInstanceError):
        assignment.job_id = "replacement"  # type: ignore[misc]


def test_execution_assignment_requires_dispatch_before_claim():
    with pytest.raises(ValueError, match="claimed_at requires dispatched_at"):
        RunExecutionAssignment(
            run_id="run-1",
            job_id="job-1",
            backend="rq",
            queue_name="default",
            created_at=datetime.now(timezone.utc),
            claimed_at=datetime.now(timezone.utc),
        )


def test_execution_assignment_requires_a_claimed_job_for_cancellation_intent():
    now = datetime.now(timezone.utc)

    with pytest.raises(ValueError, match="cancellation request requires claimed_at"):
        RunExecutionAssignment(
            run_id="run-1",
            job_id="job-1",
            backend="rq",
            queue_name="default",
            created_at=now,
            cancellation_requested_at=now,
            cancellation_reason="User requested cancellation.",
        )


@pytest.mark.parametrize(
    ("requested_at", "reason", "acknowledged_at", "message"),
    [
        (None, "User requested cancellation.", None, "request and reason"),
        (datetime.now(timezone.utc), None, None, "request and reason"),
        (None, None, datetime.now(timezone.utc), "acknowledgement requires"),
        (datetime.now(timezone.utc), "   ", None, "non-empty"),
    ],
)
def test_execution_assignment_rejects_incomplete_cancellation_evidence(
    requested_at,
    reason,
    acknowledged_at,
    message,
):
    now = datetime.now(timezone.utc)

    with pytest.raises(ValueError, match=message):
        RunExecutionAssignment(
            run_id="run-1",
            job_id="job-1",
            backend="rq",
            queue_name="default",
            created_at=now,
            dispatched_at=now,
            claimed_at=now,
            cancellation_requested_at=requested_at,
            cancellation_reason=reason,
            cancellation_acknowledged_at=acknowledged_at,
        )


def test_execution_assignment_accepts_complete_cancellation_evidence():
    now = datetime.now(timezone.utc)

    assignment = RunExecutionAssignment(
        run_id="run-1",
        job_id="job-1",
        backend="rq",
        queue_name="default",
        created_at=now,
        dispatched_at=now,
        claimed_at=now,
        cancellation_requested_at=now,
        cancellation_reason="User requested cancellation.",
        cancellation_acknowledged_at=now,
    )

    assert assignment.cancellation_reason == "User requested cancellation."


def test_execution_claim_records_whether_this_worker_acquired_ownership():
    now = datetime.now(timezone.utc)
    assignment = RunExecutionAssignment(
        run_id="run-1",
        job_id="job-1",
        backend="rq",
        queue_name="default",
        created_at=now,
        dispatched_at=now,
        claimed_at=now,
    )

    assert RunExecutionClaim(assignment=assignment, acquired=True).acquired is True


def test_execution_job_id_is_deterministic_distinct_and_rq_safe():
    first = build_execution_job_id("run/with spaces:and-unicode-测试")

    assert first == build_execution_job_id("run/with spaces:and-unicode-测试")
    assert first != build_execution_job_id("another-run")
    assert re.fullmatch(r"[A-Za-z0-9_-]+", first)


@pytest.mark.parametrize("run_id", ["", "   ", None])
def test_execution_job_id_rejects_missing_run_identity(run_id):
    with pytest.raises(ValueError, match="run_id must be a non-empty string"):
        build_execution_job_id(run_id)
