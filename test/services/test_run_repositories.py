"""Contract tests shared by run repository implementations."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timezone

import pytest

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
