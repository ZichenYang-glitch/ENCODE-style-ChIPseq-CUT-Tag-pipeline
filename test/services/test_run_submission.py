"""Tests for durable planned-run submission orchestration."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import replace

import pytest

import encode_pipeline.services.run_submission as run_submission_module
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.registry import WorkflowRegistry
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
    RunExecutionUnavailableError,
    RunNotReadyError,
    RunStartConflictError,
    RunSubmissionService,
    RunSubmissionUnavailableError,
    RunWorkflowBuildChangedError,
)
from encode_pipeline.services.run_repositories import ConcurrentRunUpdateError
from encode_pipeline.services.runs import RunService


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


class RecordingRunQueue:
    backend = "rq"
    queue_name = "test-runs"

    def __init__(self) -> None:
        self.assignments: list[RunExecutionAssignment] = []
        self.error: Exception | None = None
        self.callback: Callable[[RunExecutionAssignment], None] | None = None
        self.returned_job_id: str | None = None

    def enqueue_execution(self, assignment: RunExecutionAssignment) -> str:
        self.assignments.append(assignment)
        if self.callback is not None:
            self.callback(assignment)
        if self.error is not None:
            raise self.error
        return (
            assignment.job_id if self.returned_job_id is None else self.returned_job_id
        )


def _service() -> RunService:
    return RunService(
        registry=create_default_workflow_registry(),
        id_factory=lambda: "run-1",
    )


def _submission(run_service: RunService, queue: RecordingRunQueue):
    return RunSubmissionService(
        run_service,
        queue,
        build_identity_provider=create_default_workflow_build_identity_provider(
            registry=run_service.registry
        ),
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
    submission = _submission(run_service, queue)

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
        _submission(run_service, queue).start_run("run-1")

    assert run_service.get_run("run-1").status is RunStatus.PLANNED
    assert run_service.get_execution_assignment("run-1") is None
    assert queue.assignments == []


def test_start_rechecks_current_build_before_reserving_or_enqueueing(monkeypatch):
    run_service = _service()
    _planned_run(run_service)
    queue = RecordingRunQueue()
    provider = create_default_workflow_build_identity_provider(
        registry=run_service.registry
    )
    monkeypatch.setattr(
        provider,
        "capture",
        lambda _workflow_id: Result.failure(
            [
                Issue(
                    code="PRIVATE_RUNTIME_FAILURE",
                    message="private runtime detail",
                )
            ]
        ),
    )
    submission = RunSubmissionService(
        run_service,
        queue,
        build_identity_provider=provider,
    )

    with pytest.raises(RunExecutionUnavailableError):
        submission.start_run("run-1")

    assert run_service.get_run("run-1").status is RunStatus.PLANNED
    assert run_service.get_execution_assignment("run-1") is None
    assert queue.assignments == []


def test_start_rejects_changed_build_before_reserving_or_enqueueing(monkeypatch):
    run_service = _service()
    _planned_run(run_service)
    queue = RecordingRunQueue()
    persisted = run_service.get_workflow_build_identity("run-1")
    assert persisted is not None
    provider = create_default_workflow_build_identity_provider(
        registry=run_service.registry
    )
    monkeypatch.setattr(
        provider,
        "capture",
        lambda _workflow_id: Result.success(replace(persisted, digest="f" * 64)),
    )
    submission = RunSubmissionService(
        run_service,
        queue,
        build_identity_provider=provider,
    )

    with pytest.raises(RunWorkflowBuildChangedError):
        submission.start_run("run-1")

    assert run_service.get_run("run-1").status is RunStatus.PLANNED
    assert run_service.get_execution_assignment("run-1") is None
    assert queue.assignments == []


def test_queue_failure_preserves_planned_reservation_for_retry():
    run_service = _service()
    _planned_run(run_service)
    queue = RecordingRunQueue()
    queue.error = RunQueueUnavailableError("redis://secret-host is unavailable")
    submission = _submission(run_service, queue)

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
    submission = _submission(run_service, queue)

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
    submission = _submission(run_service, queue)

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
    submission = _submission(run_service, queue)
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
    submission = _submission(run_service, queue)

    with pytest.raises(RunStartConflictError) as raised:
        submission.start_run("run-1")

    assert raised.value.record.status is RunStatus.PLANNED
    assert run_service.get_run("run-1").status is RunStatus.PLANNED


def test_start_requires_planned_run_and_confirmed_terminal_dispatch():
    run_service = _service()
    run_service.create_run(WORKFLOW_ID, WorkflowInputs(config={}))
    submission = _submission(run_service, RecordingRunQueue())

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
    submission = _submission(run_service, queue)

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


def test_start_fails_closed_when_planned_workflow_disappears(
    monkeypatch,
) -> None:
    run_service = _service()
    _planned_run(run_service)
    queue = RecordingRunQueue()
    submission = _submission(run_service, queue)

    monkeypatch.setattr(
        WorkflowRegistry,
        "get",
        lambda _registry, _workflow_id: (_ for _ in ()).throw(
            KeyError("workflow disappeared")
        ),
    )

    with pytest.raises(RunExecutionUnavailableError):
        submission.start_run("run-1")

    assert run_service.get_execution_assignment("run-1") is None
    assert queue.assignments == []


def test_start_fails_closed_when_execution_admission_becomes_unavailable(
    monkeypatch,
) -> None:
    run_service = _service()
    _planned_run(run_service)
    queue = RecordingRunQueue()
    submission = _submission(run_service, queue)

    monkeypatch.setattr(
        run_submission_module,
        "resolve_workflow_availability",
        lambda _adapter: type(
            "Unavailable",
            (),
            {"execution": "unavailable"},
        )(),
    )

    with pytest.raises(RunExecutionUnavailableError):
        submission.start_run("run-1")

    assert run_service.get_execution_assignment("run-1") is None
    assert queue.assignments == []


def test_terminal_race_without_dispatch_is_not_reported_as_submitted(
    monkeypatch,
) -> None:
    run_service = _service()
    _planned_run(run_service)
    submission = _submission(run_service, RecordingRunQueue())
    original_resolve = submission._resolve_assignment

    def reserve_then_cancel(record):
        assignment = original_resolve(record)
        run_service.cancel_run(record.run_id, reason="cancel won before dispatch")
        return assignment

    monkeypatch.setattr(submission, "_resolve_assignment", reserve_then_cancel)

    with pytest.raises(RunStartConflictError, match="confirmed dispatch"):
        submission.start_run("run-1")

    assignment = run_service.get_execution_assignment("run-1")
    assert assignment is not None
    assert assignment.dispatched_at is None


def test_nonterminal_snapshot_race_outside_submission_states_fails_closed(
    monkeypatch,
) -> None:
    run_service = _service()
    _planned_run(run_service)
    submission = _submission(run_service, RecordingRunQueue())
    original_resolve = submission._resolve_assignment
    original_get = run_service.get_run
    reserved = False

    def reserve_then_expose_race(record):
        nonlocal reserved
        assignment = original_resolve(record)
        reserved = True
        return assignment

    def raced_snapshot(run_id):
        record = original_get(run_id)
        if reserved:
            return replace(record, status=RunStatus.VALIDATING)
        return record

    monkeypatch.setattr(submission, "_resolve_assignment", reserve_then_expose_race)
    monkeypatch.setattr(run_service, "get_run", raced_snapshot)

    with pytest.raises(RunStartConflictError, match="changed before"):
        submission.start_run("run-1")


def test_generic_queue_error_returns_confirmed_terminal_worker_progress():
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
        run_service.transition_run(assignment.run_id, RunStatus.RUNNING)

    queue.callback = worker_advances
    queue.error = RunQueueUnavailableError("connection lost after worker claim")

    result = _submission(run_service, queue).start_run("run-1")

    assert result.status is RunStatus.RUNNING


def test_queue_returning_a_different_job_identity_is_a_conflict():
    run_service = _service()
    _planned_run(run_service)
    queue = RecordingRunQueue()
    queue.returned_job_id = "rq-unexpected-job"

    with pytest.raises(RunStartConflictError, match="unexpected execution job"):
        _submission(run_service, queue).start_run("run-1")

    assignment = run_service.get_execution_assignment("run-1")
    assert assignment is not None
    assert assignment.dispatched_at is None


def test_concurrent_queue_commit_reconciles_confirmed_matching_state(
    monkeypatch,
) -> None:
    run_service = _service()
    _planned_run(run_service)
    submission = _submission(run_service, RecordingRunQueue())
    original_queue = run_service.queue_dispatched_run

    def commit_then_report_race(*args, **kwargs):
        original_queue(*args, **kwargs)
        raise ConcurrentRunUpdateError("commit acknowledgement raced")

    monkeypatch.setattr(
        run_service,
        "queue_dispatched_run",
        commit_then_report_race,
    )

    result = submission.start_run("run-1")

    assert result.status is RunStatus.QUEUED


def test_concurrent_queue_commit_without_matching_state_is_a_conflict(
    monkeypatch,
) -> None:
    run_service = _service()
    _planned_run(run_service)
    submission = _submission(run_service, RecordingRunQueue())
    monkeypatch.setattr(
        run_service,
        "mark_execution_dispatched",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            ConcurrentRunUpdateError("dispatch marker lost")
        ),
    )

    with pytest.raises(RunStartConflictError, match="being persisted"):
        submission.start_run("run-1")


def test_invalid_queue_commit_is_redacted_as_a_start_conflict(
    monkeypatch,
) -> None:
    run_service = _service()
    _planned_run(run_service)
    submission = _submission(run_service, RecordingRunQueue())
    monkeypatch.setattr(
        run_service,
        "mark_execution_dispatched",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            ValueError("private assignment detail")
        ),
    )

    with pytest.raises(
        RunStartConflictError,
        match="assignment is inconsistent",
    ):
        submission.start_run("run-1")


def test_assignment_reservation_race_reuses_confirmed_matching_dispatch(
    monkeypatch,
) -> None:
    run_service = _service()
    _planned_run(run_service)
    queue = RecordingRunQueue()
    submission = _submission(run_service, queue)
    original_reserve = run_service.ensure_execution_assignment

    def reserve_and_dispatch_then_race(run_id, *, backend, queue_name):
        assignment = original_reserve(
            run_id,
            backend=backend,
            queue_name=queue_name,
        )
        run_service.mark_execution_dispatched(run_id, job_id=assignment.job_id)
        run_service.queue_dispatched_run(
            run_id,
            job_id=assignment.job_id,
            backend=backend,
            queue_name=queue_name,
        )
        raise ConcurrentRunUpdateError("reservation acknowledgement raced")

    monkeypatch.setattr(
        run_service,
        "ensure_execution_assignment",
        reserve_and_dispatch_then_race,
    )

    result = submission.start_run("run-1")

    assert result.status is RunStatus.QUEUED
    assert len(queue.assignments) == 1


def test_assignment_reservation_race_without_dispatch_is_a_conflict(
    monkeypatch,
) -> None:
    run_service = _service()
    _planned_run(run_service)
    submission = _submission(run_service, RecordingRunQueue())
    monkeypatch.setattr(
        run_service,
        "ensure_execution_assignment",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            ValueError("private reservation detail")
        ),
    )

    with pytest.raises(RunStartConflictError, match="assignment was reserved"):
        submission.start_run("run-1")


def test_queued_assignment_must_match_the_configured_queue():
    run_service = _service()
    _planned_run(run_service)
    assignment = run_service.ensure_execution_assignment(
        "run-1",
        backend="other-backend",
        queue_name="other-queue",
    )
    run_service.mark_execution_dispatched("run-1", job_id=assignment.job_id)
    run_service.queue_dispatched_run(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
    )

    with pytest.raises(RunStartConflictError, match="does not match"):
        _submission(run_service, RecordingRunQueue()).start_run("run-1")


def test_latest_snapshot_falls_back_if_run_disappears(monkeypatch):
    run_service = _service()
    _planned_run(run_service)
    fallback = run_service.get_run("run-1")
    submission = _submission(run_service, RecordingRunQueue())
    monkeypatch.setattr(
        run_service,
        "get_run",
        lambda _run_id: (_ for _ in ()).throw(KeyError("run disappeared")),
    )

    assert submission._latest("run-1", fallback) is fallback
