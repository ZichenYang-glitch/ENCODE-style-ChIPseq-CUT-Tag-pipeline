"""Shared recovery transaction tests for SQL and in-memory repositories."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timedelta, timezone

import pytest

from encode_pipeline.persistence import (
    SqlAlchemyRunRepository,
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.results import Issue
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.services.run_repositories import (
    ConcurrentRunUpdateError,
    InMemoryRunRepository,
    RunEventDraft,
)


CREATED_AT = datetime(2026, 8, 9, 10, 0, tzinfo=timezone.utc)
DISPATCHED_AT = CREATED_AT + timedelta(minutes=1)
CLAIMED_AT = CREATED_AT + timedelta(minutes=2)
REQUESTED_AT = CREATED_AT + timedelta(minutes=3)
CONFIRMED_AT = CREATED_AT + timedelta(minutes=4)
FAILED_AT = CREATED_AT + timedelta(minutes=5)
RESULT_ATTEMPT_ID = "resultattempt-" + "a" * 64


@pytest.fixture(params=("memory", "sql"), ids=("in-memory", "sql"))
def repository_case(request, tmp_path):
    if request.param == "memory":
        yield InMemoryRunRepository(), request.param
        return

    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    try:
        yield SqlAlchemyRunRepository(create_session_factory(engine)), request.param
    finally:
        engine.dispose()


def test_sql_and_memory_prepare_and_confirm_are_atomic_and_idempotent(
    repository_case,
) -> None:
    repository, _kind = repository_case
    assignment = _create_active_run(repository, status=RunStatus.QUEUED)
    initial_result_state = repository.get_result_state("run-1")

    prepared = repository.prepare_execution_requeue(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=assignment,
        requested_at=REQUESTED_AT,
        event=_recovery_event("requeue_requested", RunStatus.QUEUED),
    )

    assert prepared.created is True
    assert prepared.assignment == replace(
        assignment,
        requeue_requested_at=REQUESTED_AT,
    )
    assert repository.get_execution_assignment("run-1") == prepared.assignment
    assert [event.event_type for event in repository.list_events("run-1")] == [
        "run_created",
        "requeue_requested",
    ]

    retry = repository.prepare_execution_requeue(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=prepared.assignment,
        requested_at=REQUESTED_AT + timedelta(minutes=10),
        event=_recovery_event("duplicate_request", RunStatus.QUEUED),
    )
    assert retry.created is False
    assert retry.assignment == prepared.assignment

    confirmed = repository.confirm_execution_requeue(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=prepared.assignment,
        confirmed_at=CONFIRMED_AT,
        event=_recovery_event("requeue_confirmed", RunStatus.QUEUED),
    )
    assert confirmed == replace(
        assignment,
        requeue_requested_at=REQUESTED_AT,
        requeue_confirmed_at=CONFIRMED_AT,
    )
    assert repository.get_execution_assignment("run-1") == confirmed

    duplicate_confirmation = repository.confirm_execution_requeue(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=confirmed,
        confirmed_at=CONFIRMED_AT + timedelta(minutes=10),
        event=_recovery_event("duplicate_confirmation", RunStatus.QUEUED),
    )
    assert duplicate_confirmation == confirmed
    assert [event.event_type for event in repository.list_events("run-1")] == [
        "run_created",
        "requeue_requested",
        "requeue_confirmed",
    ]
    assert repository.get_result_state("run-1") == initial_result_state


def test_sql_and_memory_failure_is_atomic_and_preserves_non_lifecycle_evidence(
    repository_case,
) -> None:
    repository, _kind = repository_case
    assignment = _create_active_run(repository, status=RunStatus.RUNNING)
    result_state = repository.begin_artifact_result_attempt(
        "run-1",
        attempt_id=RESULT_ATTEMPT_ID,
        expected_status=RunStatus.RUNNING,
    )
    current = repository.get_run("run-1")
    failed = _failed_record(current)
    failure_event = _recovery_event(
        "run_failed_by_recovery",
        RunStatus.FAILED,
        issue=failed.error,
    )

    assert repository.fail_run_by_recovery(
        failed,
        expected_status=RunStatus.RUNNING,
        expected_assignment=assignment,
        event=failure_event,
    )

    persisted = repository.get_run("run-1")
    assert persisted == failed
    assert repository.get_execution_assignment("run-1") == assignment
    assert repository.get_result_state("run-1") == result_state
    assert persisted.workflow_id == current.workflow_id
    assert persisted.inputs == current.inputs
    assert persisted.created_at == current.created_at
    assert persisted.started_at == current.started_at
    assert persisted.current_stage == current.current_stage
    assert persisted.cancellation_reason == current.cancellation_reason
    assert persisted.tags == current.tags
    events = repository.list_events("run-1")
    assert [event.event_type for event in events] == [
        "run_created",
        "run_failed_by_recovery",
    ]
    assert events[-1].issue == failed.error


@pytest.mark.parametrize("operation", ("prepare", "confirm", "fail"))
def test_sql_and_memory_recovery_rolls_back_when_audit_event_write_fails(
    repository_case,
    monkeypatch,
    operation,
) -> None:
    repository, kind = repository_case
    status = RunStatus.RUNNING if operation == "fail" else RunStatus.QUEUED
    expected_assignment = _create_active_run(repository, status=status)

    if operation == "confirm":
        expected_assignment = repository.prepare_execution_requeue(
            "run-1",
            expected_status=RunStatus.QUEUED,
            expected_assignment=expected_assignment,
            requested_at=REQUESTED_AT,
            event=_recovery_event("requeue_requested", RunStatus.QUEUED),
        ).assignment

    before_record = repository.get_run("run-1")
    before_assignment = repository.get_execution_assignment("run-1")
    before_events = repository.list_events("run-1")
    before_result_state = repository.get_result_state("run-1")

    def reject_event_write(*_args, **_kwargs):
        raise RuntimeError("injected audit event failure")

    event_method = "_make_event" if kind == "memory" else "_insert_event"
    monkeypatch.setattr(repository, event_method, reject_event_write)

    with pytest.raises(RuntimeError, match="injected audit event failure"):
        if operation == "prepare":
            repository.prepare_execution_requeue(
                "run-1",
                expected_status=RunStatus.QUEUED,
                expected_assignment=expected_assignment,
                requested_at=REQUESTED_AT,
                event=_recovery_event("requeue_requested", RunStatus.QUEUED),
            )
        elif operation == "confirm":
            repository.confirm_execution_requeue(
                "run-1",
                expected_status=RunStatus.QUEUED,
                expected_assignment=expected_assignment,
                confirmed_at=CONFIRMED_AT,
                event=_recovery_event("requeue_confirmed", RunStatus.QUEUED),
            )
        else:
            repository.fail_run_by_recovery(
                _failed_record(before_record),
                expected_status=RunStatus.RUNNING,
                expected_assignment=expected_assignment,
                event=_recovery_event("run_failed_by_recovery", RunStatus.FAILED),
            )

    assert repository.get_run("run-1") == before_record
    assert repository.get_execution_assignment("run-1") == before_assignment
    assert repository.list_events("run-1") == before_events
    assert repository.get_result_state("run-1") == before_result_state


def test_sql_and_memory_prepare_rejects_claim_races_and_stale_cas(
    repository_case,
) -> None:
    repository, _kind = repository_case
    stale_assignment = _create_active_run(repository, status=RunStatus.QUEUED)
    events_before_claim = repository.list_events("run-1")
    claimed = repository.claim_execution_assignment(
        "run-1",
        job_id=stale_assignment.job_id,
        backend=stale_assignment.backend,
        queue_name=stale_assignment.queue_name,
        claimed_at=CLAIMED_AT,
        allowed_statuses=frozenset({RunStatus.QUEUED}),
        event=_recovery_event("execution_claimed", RunStatus.QUEUED),
    ).assignment
    events_after_claim = repository.list_events("run-1")
    assert len(events_after_claim) == len(events_before_claim) + 1

    with pytest.raises(ConcurrentRunUpdateError, match="recovery was being committed"):
        repository.prepare_execution_requeue(
            "run-1",
            expected_status=RunStatus.QUEUED,
            expected_assignment=stale_assignment,
            requested_at=REQUESTED_AT,
            event=_recovery_event("requeue_requested", RunStatus.QUEUED),
        )
    with pytest.raises(ValueError, match="claimed.*cannot be requeued"):
        repository.prepare_execution_requeue(
            "run-1",
            expected_status=RunStatus.QUEUED,
            expected_assignment=claimed,
            requested_at=REQUESTED_AT,
            event=_recovery_event("requeue_requested", RunStatus.QUEUED),
        )

    assert repository.get_execution_assignment("run-1") == claimed
    assert repository.list_events("run-1") == events_after_claim


def test_sql_and_memory_confirmation_preserves_a_monotonic_worker_race(
    repository_case,
) -> None:
    repository, _kind = repository_case
    assignment = _create_active_run(repository, status=RunStatus.QUEUED)
    prepared = repository.prepare_execution_requeue(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=assignment,
        requested_at=REQUESTED_AT,
        event=_recovery_event("requeue_requested", RunStatus.QUEUED),
    ).assignment
    claimed = repository.claim_execution_assignment(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
        claimed_at=CLAIMED_AT,
        allowed_statuses=frozenset({RunStatus.QUEUED}),
        event=_recovery_event("execution_claimed", RunStatus.QUEUED),
    ).assignment
    queued = repository.get_run("run-1")
    repository.update_run(
        replace(
            queued,
            status=RunStatus.RUNNING,
            updated_at=CLAIMED_AT,
            started_at=CLAIMED_AT,
        ),
        expected_status=RunStatus.QUEUED,
        event=_recovery_event("run_started", RunStatus.RUNNING),
    )

    confirmed = repository.confirm_execution_requeue(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=prepared,
        confirmed_at=CONFIRMED_AT,
        event=_recovery_event("requeue_confirmed", RunStatus.RUNNING),
    )

    assert confirmed == replace(claimed, requeue_confirmed_at=CONFIRMED_AT)
    assert repository.get_execution_assignment("run-1") == confirmed
    assert repository.get_run("run-1").status is RunStatus.RUNNING
    assert [event.event_type for event in repository.list_events("run-1")][-3:] == [
        "execution_claimed",
        "run_started",
        "requeue_confirmed",
    ]


def test_sql_and_memory_failure_rejects_stale_assignment_cas(
    repository_case,
) -> None:
    repository, _kind = repository_case
    stale_assignment = _create_active_run(repository, status=RunStatus.RUNNING)
    repository.request_execution_cancellation(
        "run-1",
        job_id=stale_assignment.job_id,
        backend=stale_assignment.backend,
        queue_name=stale_assignment.queue_name,
        requested_at=REQUESTED_AT,
        reason="operator request",
        event=_recovery_event("cancellation_requested", RunStatus.RUNNING),
    )
    current = repository.get_run("run-1")
    assignment_after_race = repository.get_execution_assignment("run-1")
    events_after_race = repository.list_events("run-1")

    with pytest.raises(ConcurrentRunUpdateError, match="recovery was being committed"):
        repository.fail_run_by_recovery(
            _failed_record(current),
            expected_status=RunStatus.RUNNING,
            expected_assignment=stale_assignment,
            event=_recovery_event("run_failed_by_recovery", RunStatus.FAILED),
        )

    assert repository.get_run("run-1") == current
    assert repository.get_execution_assignment("run-1") == assignment_after_race
    assert repository.list_events("run-1") == events_after_race


def test_sql_and_memory_failure_rejects_immutable_evidence_changes(
    repository_case,
) -> None:
    repository, _kind = repository_case
    assignment = _create_active_run(repository, status=RunStatus.RUNNING)
    current = repository.get_run("run-1")
    valid_failure = _failed_record(current)
    invalid_failures = (
        replace(valid_failure, workflow_id="other-workflow"),
        replace(valid_failure, inputs={"tampered": True}),
        replace(valid_failure, created_at=CREATED_AT - timedelta(minutes=1)),
        replace(valid_failure, started_at=CLAIMED_AT + timedelta(seconds=1)),
        replace(valid_failure, current_stage="tampered-stage"),
        replace(valid_failure, cancellation_reason="tampered reason"),
        replace(valid_failure, tags={"owner": "tampered"}),
    )
    original_events = repository.list_events("run-1")

    for invalid_failure in invalid_failures:
        with pytest.raises(
            ValueError,
            match="changes immutable run evidence",
        ):
            repository.fail_run_by_recovery(
                invalid_failure,
                expected_status=RunStatus.RUNNING,
                expected_assignment=assignment,
                event=_recovery_event("run_failed_by_recovery", RunStatus.FAILED),
            )
        assert repository.get_run("run-1") == current
        assert repository.get_execution_assignment("run-1") == assignment
        assert repository.list_events("run-1") == original_events


def _create_active_run(
    repository,
    *,
    status: RunStatus,
) -> RunExecutionAssignment:
    record = RunRecord(
        run_id="run-1",
        workflow_id="workflow-a",
        inputs={"config": {"private_path": "/private/input"}},
        status=status,
        created_at=CREATED_AT,
        updated_at=CREATED_AT,
        started_at=CLAIMED_AT if status is RunStatus.RUNNING else None,
        ended_at=None,
        current_stage="execution",
        cancellation_reason=None,
        error=None,
        tags={"owner": "test"},
    )
    repository.create_run(
        record,
        _recovery_event("run_created", status),
    )
    assignment = RunExecutionAssignment(
        run_id=record.run_id,
        job_id="job-1",
        backend="rq",
        queue_name="default",
        created_at=CREATED_AT,
        dispatched_at=DISPATCHED_AT,
        claimed_at=CLAIMED_AT if status is RunStatus.RUNNING else None,
    )
    return repository.ensure_execution_assignment(
        assignment,
        expected_status=status,
    )


def _failed_record(current: RunRecord) -> RunRecord:
    return replace(
        current,
        status=RunStatus.FAILED,
        updated_at=FAILED_AT,
        ended_at=FAILED_AT,
        error=Issue(
            code="RUN_RECOVERY_FAILED",
            message="Run failed by administrator recovery.",
            source="platform",
        ),
    )


def _recovery_event(
    event_type: str,
    status: RunStatus,
    *,
    issue: Issue | None = None,
) -> RunEventDraft:
    return RunEventDraft(
        event_type=event_type,
        message=event_type.replace("_", " ").capitalize() + ".",
        status=status,
        stage="execution",
        context={"reason_code": event_type.upper()},
        issue=issue,
    )
