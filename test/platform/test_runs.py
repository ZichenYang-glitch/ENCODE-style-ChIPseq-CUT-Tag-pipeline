from dataclasses import FrozenInstanceError

import pytest


def test_run_status_values_and_terminal_detection():
    from encode_pipeline.platform.runs import RunStatus

    assert RunStatus.CREATED.value == "created"
    assert not RunStatus.CREATED.is_terminal
    assert RunStatus.SUCCEEDED.is_terminal
    assert RunStatus.FAILED.is_terminal
    assert RunStatus.CANCELLED.is_terminal


def test_can_transition_accepts_pr99_path():
    from encode_pipeline.platform.runs import RunStatus, can_transition

    assert can_transition(RunStatus.CREATED, RunStatus.VALIDATING)
    assert can_transition(RunStatus.VALIDATING, RunStatus.PLANNED)
    assert can_transition(RunStatus.VALIDATING, RunStatus.FAILED)
    assert can_transition(RunStatus.PLANNED, RunStatus.QUEUED)
    assert can_transition(RunStatus.QUEUED, RunStatus.RUNNING)
    assert can_transition(RunStatus.RUNNING, RunStatus.SUCCEEDED)
    assert can_transition(RunStatus.RUNNING, RunStatus.FAILED)
    assert can_transition(RunStatus.RUNNING, RunStatus.CANCELLED)


def test_can_transition_rejects_invalid_transitions():
    from encode_pipeline.platform.runs import RunStatus, can_transition

    assert not can_transition(RunStatus.CREATED, RunStatus.RUNNING)
    assert not can_transition(RunStatus.RUNNING, RunStatus.PLANNED)
    assert not can_transition(RunStatus.SUCCEEDED, RunStatus.FAILED)
    assert not can_transition(RunStatus.FAILED, RunStatus.CANCELLED)


def test_require_transition_raises_on_invalid():
    from encode_pipeline.platform.runs import RunStatus, require_transition

    with pytest.raises(ValueError, match="Illegal transition"):
        require_transition(RunStatus.SUCCEEDED, RunStatus.FAILED)


def test_run_record_stores_status_and_error():
    from datetime import datetime, timezone
    from encode_pipeline.platform.runs import RunRecord, RunStatus
    from encode_pipeline.platform.results import Issue

    now = datetime.now(timezone.utc)
    error = Issue(code="RUN_FAILED", message="Run failed")
    record = RunRecord(
        run_id="run-1",
        workflow_id="wf-1",
        inputs={"config": {"samples": "samples.tsv"}},
        status=RunStatus.FAILED,
        created_at=now,
        updated_at=now,
        started_at=None,
        ended_at=now,
        current_stage=None,
        cancellation_reason=None,
        error=error,
        tags={"env": "test"},
    )

    assert record.run_id == "run-1"
    assert record.status == RunStatus.FAILED
    assert record.error == error
    assert record.ended_at == now


def test_run_record_is_frozen():
    from datetime import datetime, timezone
    from encode_pipeline.platform.runs import RunRecord, RunStatus

    now = datetime.now(timezone.utc)
    record = RunRecord(
        run_id="run-1",
        workflow_id="wf-1",
        inputs={},
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

    with pytest.raises(FrozenInstanceError):
        record.status = RunStatus.RUNNING


def test_run_record_to_dict_is_json_ready():
    from datetime import datetime, timezone
    from encode_pipeline.platform.runs import RunRecord, RunStatus

    now = datetime.now(timezone.utc)
    record = RunRecord(
        run_id="run-1",
        workflow_id="wf-1",
        inputs={"config": {}},
        status=RunStatus.CREATED,
        created_at=now,
        updated_at=now,
        started_at=None,
        ended_at=None,
        current_stage=None,
        cancellation_reason=None,
        error=None,
        tags={"env": "test"},
    )

    data = record.to_dict()
    assert data["run_id"] == "run-1"
    assert data["status"] == "created"
    assert data["created_at"] is now
    assert data["inputs"] == {"config": {}}
    assert data["tags"] == {"env": "test"}


def test_run_event_stores_sequence_and_issue():
    from datetime import datetime, timezone
    from encode_pipeline.platform.runs import RunEvent, RunStatus
    from encode_pipeline.platform.results import Issue

    now = datetime.now(timezone.utc)
    issue = Issue(code="STAGE_WARN", message="Stage warning", severity="warning")
    event = RunEvent(
        event_id="evt-1",
        run_id="run-1",
        sequence=3,
        event_type="issue_added",
        timestamp=now,
        status=RunStatus.RUNNING,
        stage="alignment",
        message="Warning during alignment.",
        context={"sample": "S1"},
        issue=issue,
    )

    assert event.sequence == 3
    assert event.event_type == "issue_added"
    assert event.issue == issue
    assert event.to_dict()["sequence"] == 3
    assert event.to_dict()["issue"]["code"] == "STAGE_WARN"


def test_run_event_context_is_defensively_copied():
    from datetime import datetime, timezone
    from encode_pipeline.platform.runs import RunEvent, RunStatus

    context = {"sample": "S1"}
    event = RunEvent(
        event_id="evt-1",
        run_id="run-1",
        sequence=1,
        event_type="status_changed",
        timestamp=datetime.now(timezone.utc),
        status=RunStatus.CREATED,
        stage=None,
        message="Run created.",
        context=context,
        issue=None,
    )

    context["sample"] = "S2"
    assert event.context == {"sample": "S1"}
