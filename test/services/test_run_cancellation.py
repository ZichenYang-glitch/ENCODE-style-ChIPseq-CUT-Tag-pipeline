"""Tests for durable RUNNING cancellation state ownership."""

from __future__ import annotations

from encode_pipeline.platform.adapters import WorkflowInputs
import pytest

from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.defaults import create_default_workflow_registry
from encode_pipeline.services.run_cancellation import (
    RunCancellationConflictError,
    RunCancellationService,
    RunCancellationUnavailableError,
)
from encode_pipeline.services.run_queue import RunQueueStopUnavailableError
from encode_pipeline.services.runs import RunService


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
REASON = "User requested cancellation."


class RecordingStopQueue:
    backend = "rq"
    queue_name = "runs"

    def __init__(self, run_service: RunService) -> None:
        self.run_service = run_service
        self.assignments: list[RunExecutionAssignment] = []
        self.error: Exception | None = None
        self.callback = None

    def request_stop(self, assignment: RunExecutionAssignment) -> None:
        persisted = self.run_service.get_execution_assignment(assignment.run_id)
        assert persisted == assignment
        assert persisted.cancellation_requested_at is not None
        assert self.run_service.get_run(assignment.run_id).status is RunStatus.RUNNING
        self.assignments.append(assignment)
        if self.callback is not None:
            self.callback(assignment)
        if self.error is not None:
            raise self.error


def test_running_cancellation_persists_intent_without_premature_terminal_state():
    service, assignment = _running_service()
    events_before = service.list_events("run-1", limit=100)

    requested = service.request_execution_cancellation(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
        reason=REASON,
    )

    assert requested.created is True
    assert requested.record.status is RunStatus.RUNNING
    assert requested.record.ended_at is None
    assert requested.record.cancellation_reason is None
    assert requested.assignment.cancellation_requested_at is not None
    assert requested.assignment.cancellation_reason == REASON
    assert requested.assignment.cancellation_acknowledged_at is None
    events = service.list_events("run-1", limit=100)
    assert events[:-1] == events_before
    assert events[-1].event_type == "cancellation_requested"
    assert events[-1].status is RunStatus.RUNNING
    assert events[-1].context == {
        "job_id": assignment.job_id,
        "backend": assignment.backend,
        "queue_name": assignment.queue_name,
    }


def test_repeated_running_cancellation_preserves_first_intent_and_event():
    service, assignment = _running_service()
    first = service.request_execution_cancellation(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
        reason=REASON,
    )
    events_after_first = service.list_events("run-1", limit=100)

    second = service.request_execution_cancellation(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
        reason="A later reason must not replace the durable intent.",
    )

    assert second.created is False
    assert second.assignment == first.assignment
    assert second.assignment.cancellation_reason == REASON
    assert service.list_events("run-1", limit=100) == events_after_first


def test_stop_acknowledgement_is_the_only_running_to_cancelled_path():
    service, assignment = _running_service()
    service.request_execution_cancellation(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
        reason=REASON,
    )

    acknowledged = service.acknowledge_execution_stop(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
    )

    assert acknowledged.transitioned is True
    assert acknowledged.record.status is RunStatus.CANCELLED
    assert acknowledged.record.ended_at is not None
    assert acknowledged.record.cancellation_reason == REASON
    assert acknowledged.assignment.cancellation_acknowledged_at is not None
    events = service.list_events("run-1", limit=100)
    assert events[-1].event_type == "cancellation_acknowledged"
    assert events[-1].status is RunStatus.CANCELLED
    assert events[-1].context["previous_status"] == "running"
    assert events[-1].context["new_status"] == "cancelled"
    assert events[-1].context["cancellation_reason"] == REASON

    events_after_first = events
    repeated = service.acknowledge_execution_stop(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
    )
    assert repeated.transitioned is False
    assert repeated.record == acknowledged.record
    assert repeated.assignment == acknowledged.assignment
    assert service.list_events("run-1", limit=100) == events_after_first


def test_generic_transition_cannot_bypass_stop_acknowledgement():
    service, _assignment = _running_service()

    try:
        service.transition_run("run-1", RunStatus.CANCELLED)
    except ValueError as exc:
        assert "acknowledgement" in str(exc)
    else:  # pragma: no cover - explicit guard against a lying lifecycle
        raise AssertionError("generic transition bypassed process-stop acknowledgement")

    assert service.get_run("run-1").status is RunStatus.RUNNING


def test_stop_without_user_intent_fails_instead_of_impersonating_cancellation():
    service, assignment = _running_service()

    acknowledged = service.acknowledge_execution_stop(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
    )

    assert acknowledged.transitioned is True
    assert acknowledged.record.status is RunStatus.FAILED
    assert acknowledged.record.cancellation_reason is None
    assert acknowledged.record.error is not None
    assert acknowledged.record.error.code == "RUN_WORKER_STOPPED_UNEXPECTEDLY"
    assert acknowledged.record.error.context == {
        "reason_code": "WORKER_STOP_WITHOUT_CANCELLATION"
    }
    assert acknowledged.assignment.cancellation_acknowledged_at is None
    assert service.list_events("run-1", limit=100)[-1].event_type == (
        "execution_stopped_unexpectedly"
    )


def test_natural_terminal_state_wins_before_stop_acknowledgement():
    service, assignment = _running_service()
    service.request_execution_cancellation(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
        reason=REASON,
    )
    succeeded = service.transition_run("run-1", RunStatus.SUCCEEDED)
    events_before_ack = service.list_events("run-1", limit=100)

    acknowledged = service.acknowledge_execution_stop(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
    )

    assert acknowledged.transitioned is False
    assert acknowledged.record == succeeded
    assert acknowledged.assignment.cancellation_acknowledged_at is None
    assert service.list_events("run-1", limit=100) == events_before_ack


def test_stop_acknowledgement_wins_before_natural_completion():
    service, assignment = _running_service()
    service.request_execution_cancellation(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
        reason=REASON,
    )
    acknowledged = service.acknowledge_execution_stop(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
    )

    try:
        service.transition_run("run-1", RunStatus.SUCCEEDED)
    except ValueError as exc:
        assert "Illegal transition" in str(exc)
    else:  # pragma: no cover - terminal states must be absorbing
        raise AssertionError("natural completion replaced acknowledged cancellation")

    assert service.get_run("run-1") == acknowledged.record
    terminal_events = [
        event
        for event in service.list_events("run-1", limit=100)
        if event.status is not None and event.status.is_terminal
    ]
    assert len(terminal_events) == 1


def test_wrong_execution_identity_cannot_request_or_acknowledge_cancellation():
    service, assignment = _running_service()
    record_before = service.get_run("run-1")
    events_before = service.list_events("run-1", limit=100)

    for method, kwargs in (
        (
            service.request_execution_cancellation,
            {
                "job_id": "wrong-job",
                "backend": assignment.backend,
                "queue_name": assignment.queue_name,
                "reason": REASON,
            },
        ),
        (
            service.acknowledge_execution_stop,
            {
                "job_id": assignment.job_id,
                "backend": assignment.backend,
                "queue_name": "wrong-queue",
            },
        ),
    ):
        try:
            method("run-1", **kwargs)
        except ValueError as exc:
            assert "identity" in str(exc)
        else:  # pragma: no cover - identity checks are safety critical
            raise AssertionError("wrong execution identity mutated the run")

    assert service.get_run("run-1") == record_before
    assert service.get_execution_assignment("run-1") == assignment
    assert service.list_events("run-1", limit=100) == events_before


def test_cancellation_service_sends_stop_only_after_durable_intent():
    run_service, assignment = _running_service()
    queue = RecordingStopQueue(run_service)
    cancellation = RunCancellationService(run_service, queue)

    result = cancellation.cancel_run("run-1", reason=REASON)

    assert result.stop_requested is True
    assert result.record.status is RunStatus.RUNNING
    assert result.record.ended_at is None
    assert result.record.cancellation_reason is None
    assert queue.assignments == [run_service.get_execution_assignment("run-1")]
    assert queue.assignments[0].job_id == assignment.job_id


def test_cancellation_service_retries_stop_without_duplicate_intent():
    run_service, _assignment = _running_service()
    queue = RecordingStopQueue(run_service)
    cancellation = RunCancellationService(run_service, queue)

    first = cancellation.cancel_run("run-1", reason=REASON)
    second = cancellation.cancel_run("run-1", reason="later reason")

    assert first.stop_requested is True
    assert second.stop_requested is True
    assert len(queue.assignments) == 2
    persisted = run_service.get_execution_assignment("run-1")
    assert persisted is not None
    assert persisted.cancellation_reason == REASON
    assert [
        event.event_type for event in run_service.list_events("run-1", limit=100)
    ].count("cancellation_requested") == 1


def test_cancellation_service_stop_failure_is_retryable_without_false_terminal():
    run_service, _assignment = _running_service()
    queue = RecordingStopQueue(run_service)
    queue.error = RunQueueStopUnavailableError(
        "redis://user:password@private-host:6379 is unavailable"
    )
    cancellation = RunCancellationService(run_service, queue)

    with pytest.raises(RunCancellationUnavailableError) as captured:
        cancellation.cancel_run("run-1", reason=REASON)

    assert captured.value.record.status is RunStatus.RUNNING
    assert run_service.get_run("run-1").status is RunStatus.RUNNING
    persisted = run_service.get_execution_assignment("run-1")
    assert persisted is not None
    assert persisted.cancellation_requested_at is not None
    assert persisted.cancellation_acknowledged_at is None


@pytest.mark.parametrize(
    "status",
    [RunStatus.CREATED, RunStatus.VALIDATING, RunStatus.PLANNED, RunStatus.QUEUED],
)
def test_cancellation_service_keeps_pre_running_immediate_semantics(status):
    run_service = RunService(
        create_default_workflow_registry(),
        id_factory=lambda: "run-1",
    )
    run_service.create_run(WORKFLOW_ID, WorkflowInputs(config={}))
    for next_status in (
        RunStatus.VALIDATING,
        RunStatus.PLANNED,
        RunStatus.QUEUED,
    ):
        if status is RunStatus.CREATED:
            break
        run_service.transition_run("run-1", next_status)
        if next_status is status:
            break
    queue = RecordingStopQueue(run_service)

    result = RunCancellationService(run_service, queue).cancel_run(
        "run-1",
        reason=REASON,
    )

    assert result.stop_requested is False
    assert result.record.status is RunStatus.CANCELLED
    assert queue.assignments == []


def test_cancellation_service_returns_existing_terminal_without_stop():
    run_service, _assignment = _running_service()
    succeeded = run_service.transition_run("run-1", RunStatus.SUCCEEDED)
    queue = RecordingStopQueue(run_service)

    result = RunCancellationService(run_service, queue).cancel_run(
        "run-1",
        reason=REASON,
    )

    assert result.stop_requested is False
    assert result.record == succeeded
    assert queue.assignments == []


def test_cancellation_service_rejects_assignment_for_another_queue():
    run_service, _assignment = _running_service()
    queue = RecordingStopQueue(run_service)
    queue.queue_name = "other"

    with pytest.raises(RunCancellationConflictError) as captured:
        RunCancellationService(run_service, queue).cancel_run(
            "run-1",
            reason=REASON,
        )

    assert captured.value.record.status is RunStatus.RUNNING
    assert queue.assignments == []
    persisted = run_service.get_execution_assignment("run-1")
    assert persisted is not None
    assert persisted.cancellation_requested_at is None


def test_cancellation_service_returns_callback_terminal_truthfully():
    run_service, _assignment = _running_service()
    queue = RecordingStopQueue(run_service)
    queue.callback = lambda assignment: run_service.acknowledge_execution_stop(
        assignment.run_id,
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
    )

    result = RunCancellationService(run_service, queue).cancel_run(
        "run-1",
        reason=REASON,
    )

    assert result.stop_requested is False
    assert result.record.status is RunStatus.CANCELLED


def test_cancellation_service_does_not_report_503_after_terminal_race():
    run_service, _assignment = _running_service()
    queue = RecordingStopQueue(run_service)

    def complete_then_fail(_assignment):
        run_service.transition_run("run-1", RunStatus.SUCCEEDED)

    queue.callback = complete_then_fail
    queue.error = RunQueueStopUnavailableError("job already completed")

    result = RunCancellationService(run_service, queue).cancel_run(
        "run-1",
        reason=REASON,
    )

    assert result.stop_requested is False
    assert result.record.status is RunStatus.SUCCEEDED


def _running_service() -> tuple[RunService, RunExecutionAssignment]:
    service = RunService(
        create_default_workflow_registry(),
        id_factory=lambda: "run-1",
    )
    service.create_run(WORKFLOW_ID, WorkflowInputs(config={}))
    service.transition_run("run-1", RunStatus.VALIDATING)
    service.transition_run("run-1", RunStatus.PLANNED)
    assignment = service.ensure_execution_assignment("run-1", queue_name="runs")
    service.mark_execution_dispatched("run-1", job_id=assignment.job_id)
    service.transition_run("run-1", RunStatus.QUEUED)
    service.claim_execution_assignment(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
    )
    service.transition_run("run-1", RunStatus.RUNNING)
    claimed = service.get_execution_assignment("run-1")
    assert claimed is not None
    return service, claimed
