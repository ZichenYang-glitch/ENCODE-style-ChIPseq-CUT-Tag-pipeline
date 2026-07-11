"""Contract tests shared by run repository implementations."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timezone

import pytest

from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.services.run_repositories import (
    InMemoryRunRepository,
    RunEventDraft,
)


def test_in_memory_create_is_atomic_when_event_is_invalid():
    repository = InMemoryRunRepository()

    with pytest.raises(ValueError, match="context must be a mapping"):
        repository.create_run(
            _record(),
            RunEventDraft(
                event_type="status_changed",
                message="Run created.",
                status=RunStatus.CREATED,
                context=object(),
            ),
        )

    assert repository.list_runs() == ()


def test_in_memory_update_is_atomic_when_event_is_invalid():
    repository = InMemoryRunRepository()
    record = _record()
    repository.create_run(record, _created_event())

    with pytest.raises(ValueError, match="context must be a mapping"):
        repository.update_run(
            replace(record, status=RunStatus.VALIDATING),
            expected_status=RunStatus.CREATED,
            event=RunEventDraft(
                event_type="status_changed",
                message="Run validating.",
                status=RunStatus.VALIDATING,
                context=object(),
            ),
        )

    assert repository.get_run(record.run_id) == record
    assert len(repository.list_events(record.run_id)) == 1


def test_in_memory_execution_assignment_is_idempotent_per_run():
    repository = InMemoryRunRepository()
    record = _record()
    repository.create_run(record, _created_event())
    original = _assignment(record.run_id, "job-original")
    replacement = _assignment(record.run_id, "job-replacement")

    assert repository.get_execution_assignment(record.run_id) is None
    assert (
        repository.ensure_execution_assignment(
            original, expected_status=RunStatus.CREATED
        )
        == original
    )
    assert (
        repository.ensure_execution_assignment(
            replacement, expected_status=RunStatus.CREATED
        )
        == original
    )
    assert repository.get_execution_assignment(record.run_id) == original


def test_in_memory_dispatch_mark_is_idempotent_after_status_changes():
    repository = InMemoryRunRepository()
    record = _record()
    repository.create_run(record, _created_event())
    assignment = repository.ensure_execution_assignment(
        _assignment(record.run_id, "job-1"),
        expected_status=RunStatus.CREATED,
    )
    dispatched_at = datetime.now(timezone.utc)
    dispatched = repository.mark_execution_dispatched(
        record.run_id,
        job_id=assignment.job_id,
        dispatched_at=dispatched_at,
        allowed_statuses=frozenset({RunStatus.CREATED}),
    )
    repository.update_run(
        replace(record, status=RunStatus.VALIDATING),
        expected_status=RunStatus.CREATED,
        event=RunEventDraft(
            event_type="status_changed",
            message="Run advanced.",
            status=RunStatus.VALIDATING,
        ),
    )

    retried = repository.mark_execution_dispatched(
        record.run_id,
        job_id=assignment.job_id,
        dispatched_at=dispatched_at,
        allowed_statuses=frozenset({RunStatus.CREATED}),
    )

    assert retried == dispatched
    assert retried.dispatched_at == dispatched_at


def test_in_memory_execution_assignment_rejects_cross_run_job_reuse():
    repository = InMemoryRunRepository()
    first = _record()
    second = replace(first, run_id="run-2")
    repository.create_run(first, _created_event())
    repository.create_run(second, _created_event())
    repository.ensure_execution_assignment(
        _assignment(first.run_id, "shared-job"),
        expected_status=RunStatus.CREATED,
    )

    with pytest.raises(ValueError, match="shared-job.*already assigned"):
        repository.ensure_execution_assignment(
            _assignment(second.run_id, "shared-job"),
            expected_status=RunStatus.CREATED,
        )


def test_in_memory_execution_assignment_requires_a_persisted_run():
    repository = InMemoryRunRepository()

    with pytest.raises(KeyError, match="missing"):
        repository.ensure_execution_assignment(
            _assignment("missing", "job-1"),
            expected_status=RunStatus.CREATED,
        )
    assert repository.get_execution_assignment("missing") is None


def _record() -> RunRecord:
    now = datetime.now(timezone.utc)
    return RunRecord(
        run_id="run-1",
        workflow_id="fake",
        inputs={"config": {}, "samples": None, "options": {}},
        status=RunStatus.CREATED,
        created_at=now,
        updated_at=now,
        started_at=None,
        ended_at=None,
        current_stage=None,
        cancellation_reason=None,
        error=None,
        tags={},
    )


def _created_event() -> RunEventDraft:
    return RunEventDraft(
        event_type="status_changed",
        message="Run created.",
        status=RunStatus.CREATED,
    )


def _assignment(run_id: str, job_id: str) -> RunExecutionAssignment:
    return RunExecutionAssignment(
        run_id=run_id,
        job_id=job_id,
        backend="rq",
        queue_name="default",
        created_at=datetime.now(timezone.utc),
    )
