"""Tests for durable planned-run submission orchestration."""

from __future__ import annotations

from collections.abc import Callable

import pytest

from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.defaults import (
    create_default_workflow_build_identity_provider,
    create_default_workflow_registry,
)
from encode_pipeline.services.run_queue import (
    RunQueueJobUnavailableError,
    RunQueueUnavailableError,
)
from encode_pipeline.services.run_submission import (
    RunBuildIdentityMissingError,
    RunNotReadyError,
    RunStartConflictError,
    RunSubmissionService,
    RunSubmissionUnavailableError,
)
from encode_pipeline.services.runs import RunService


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


class RecordingRunQueue:
    backend = "rq"
    queue_name = "test-runs"

    def __init__(self) -> None:
        self.assignments: list[RunExecutionAssignment] = []
        self.error: Exception | None = None
        self.callback: Callable[[RunExecutionAssignment], None] | None = None

    def enqueue_execution(self, assignment: RunExecutionAssignment) -> str:
        self.assignments.append(assignment)
        if self.callback is not None:
            self.callback(assignment)
        if self.error is not None:
            raise self.error
        return assignment.job_id


def _service() -> RunService:
    return RunService(
        registry=create_default_workflow_registry(),
        id_factory=lambda: "run-1",
    )


def _planned_run(service: RunService) -> None:
    service.create_run(WORKFLOW_ID, WorkflowInputs(config={}))
    service.transition_run("run-1", RunStatus.VALIDATING)
    identity_result = create_default_workflow_build_identity_provider(
        registry=create_default_workflow_registry()
    ).capture(WORKFLOW_ID)
    assert identity_result.is_success
    service.complete_preflight("run-1", identity_result.value)


def _legacy_planned_run(service: RunService) -> None:
    service.create_run(WORKFLOW_ID, WorkflowInputs(config={}))
    service.transition_run("run-1", RunStatus.VALIDATING)
    service.transition_run("run-1", RunStatus.PLANNED)


def test_start_run_dispatches_and_queues_exactly_once():
    run_service = _service()
    _planned_run(run_service)
    queue = RecordingRunQueue()
    submission = RunSubmissionService(run_service, queue)

    first = submission.start_run("run-1")
    events_after_first = run_service.list_events("run-1")
    second = submission.start_run("run-1")

    assert first.status is RunStatus.QUEUED
    assert second == first
    assert len(queue.assignments) == 2
    assert queue.assignments[0].run_id == queue.assignments[1].run_id
    assert queue.assignments[0].job_id == queue.assignments[1].job_id
    assert queue.assignments[0].backend == queue.assignments[1].backend
    assert queue.assignments[0].queue_name == queue.assignments[1].queue_name
    assignment = run_service.get_execution_assignment("run-1")
    assert assignment is not None
    assert assignment.dispatched_at is not None
    assert assignment.job_id == queue.assignments[0].job_id
    assert run_service.list_events("run-1") == events_after_first
    assert (
        len([event for event in events_after_first if event.status is RunStatus.QUEUED])
        == 1
    )


def test_start_rejects_legacy_planned_run_without_build_identity():
    run_service = _service()
    _legacy_planned_run(run_service)
    queue = RecordingRunQueue()

    with pytest.raises(RunBuildIdentityMissingError):
        RunSubmissionService(run_service, queue).start_run("run-1")

    assert run_service.get_run("run-1").status is RunStatus.PLANNED
    assert run_service.get_execution_assignment("run-1") is None
    assert queue.assignments == []


def test_queue_failure_preserves_planned_reservation_for_retry():
    run_service = _service()
    _planned_run(run_service)
    queue = RecordingRunQueue()
    queue.error = RunQueueUnavailableError("redis://secret-host is unavailable")
    submission = RunSubmissionService(run_service, queue)

    with pytest.raises(RunSubmissionUnavailableError) as raised:
        submission.start_run("run-1")

    assert raised.value.record.status is RunStatus.PLANNED
    assignment = run_service.get_execution_assignment("run-1")
    assert assignment is not None
    assert assignment.dispatched_at is None
    assert run_service.get_run("run-1").status is RunStatus.PLANNED

    queue.error = None
    retried = submission.start_run("run-1")
    assert retried.status is RunStatus.QUEUED
    assert queue.assignments[0].job_id == queue.assignments[1].job_id


def test_ambiguous_queue_error_keeps_confirmed_queued_state_retryable():
    run_service = _service()
    _planned_run(run_service)
    queue = RecordingRunQueue()

    def worker_repairs(assignment: RunExecutionAssignment) -> None:
        run_service.mark_execution_dispatched(
            assignment.run_id,
            job_id=assignment.job_id,
        )
        run_service.queue_dispatched_run(
            assignment.run_id,
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )

    queue.callback = worker_repairs
    queue.error = RunQueueUnavailableError("connection dropped after enqueue")
    submission = RunSubmissionService(run_service, queue)

    with pytest.raises(RunSubmissionUnavailableError) as raised:
        submission.start_run("run-1")

    assert raised.value.record.status is RunStatus.QUEUED
    assignment = run_service.get_execution_assignment("run-1")
    assert assignment is not None
    assert assignment.dispatched_at is not None


def test_backend_terminal_race_returns_confirmed_worker_progress():
    run_service = _service()
    _planned_run(run_service)
    queue = RecordingRunQueue()

    def worker_advances(assignment: RunExecutionAssignment) -> None:
        run_service.mark_execution_dispatched(
            assignment.run_id,
            job_id=assignment.job_id,
        )
        run_service.queue_dispatched_run(
            assignment.run_id,
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )
        run_service.claim_execution_assignment(
            assignment.run_id,
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )
        run_service.transition_run(assignment.run_id, RunStatus.RUNNING)

    queue.callback = worker_advances
    queue.error = RunQueueJobUnavailableError("job finished during retry")
    submission = RunSubmissionService(run_service, queue)

    result = submission.start_run("run-1")

    assert result.status is RunStatus.RUNNING


def test_backend_error_rechecks_progress_before_building_conflict(monkeypatch):
    run_service = _service()
    _planned_run(run_service)
    queue = RecordingRunQueue()

    def worker_repairs(assignment: RunExecutionAssignment) -> None:
        run_service.mark_execution_dispatched(
            assignment.run_id,
            job_id=assignment.job_id,
        )
        run_service.queue_dispatched_run(
            assignment.run_id,
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )

    queue.callback = worker_repairs
    queue.error = RunQueueJobUnavailableError("job completed during retry")
    submission = RunSubmissionService(run_service, queue)
    original_get_assignment = run_service.get_execution_assignment
    advanced = False

    def advance_after_assignment_read(run_id: str):
        nonlocal advanced
        assignment = original_get_assignment(run_id)
        if (
            not advanced
            and assignment is not None
            and assignment.dispatched_at is not None
        ):
            advanced = True
            run_service.claim_execution_assignment(
                run_id,
                job_id=assignment.job_id,
                backend=assignment.backend,
                queue_name=assignment.queue_name,
            )
            run_service.transition_run(run_id, RunStatus.RUNNING)
        return assignment

    monkeypatch.setattr(
        run_service,
        "get_execution_assignment",
        advance_after_assignment_read,
    )

    result = submission.start_run("run-1")

    assert advanced is True
    assert result.status is RunStatus.RUNNING


def test_unreusable_backend_job_is_a_start_conflict():
    run_service = _service()
    _planned_run(run_service)
    queue = RecordingRunQueue()
    queue.error = RunQueueJobUnavailableError("failed")
    submission = RunSubmissionService(run_service, queue)

    with pytest.raises(RunStartConflictError) as raised:
        submission.start_run("run-1")

    assert raised.value.record.status is RunStatus.PLANNED
    assert run_service.get_run("run-1").status is RunStatus.PLANNED


def test_start_requires_planned_run_and_confirmed_terminal_dispatch():
    run_service = _service()
    run_service.create_run(WORKFLOW_ID, WorkflowInputs(config={}))
    submission = RunSubmissionService(run_service, RecordingRunQueue())

    with pytest.raises(RunNotReadyError):
        submission.start_run("run-1")

    run_service.cancel_run("run-1")
    with pytest.raises(RunStartConflictError, match="confirmed durable"):
        submission.start_run("run-1")


def test_cancel_race_after_enqueue_persists_dispatch_without_queue_event():
    run_service = _service()
    _planned_run(run_service)
    queue = RecordingRunQueue()
    queue.callback = lambda _assignment: run_service.cancel_run("run-1", reason="race")
    submission = RunSubmissionService(run_service, queue)

    cancelled = submission.start_run("run-1")
    retried = submission.start_run("run-1")

    assert cancelled.status is RunStatus.CANCELLED
    assert retried == cancelled
    assert len(queue.assignments) == 1
    assignment = run_service.get_execution_assignment("run-1")
    assert assignment is not None
    assert assignment.dispatched_at is not None
    events = run_service.list_events("run-1")
    assert events[-1].status is RunStatus.CANCELLED
    assert not any(event.status is RunStatus.QUEUED for event in events)
