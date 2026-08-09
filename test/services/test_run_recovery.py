"""Behavior tests for explicit administrator run recovery."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timezone
import json

import pytest

from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.result_generations import RunResultState
from encode_pipeline.platform.results import Issue
from encode_pipeline.platform.run_recovery import (
    ExecutionQueueEvidenceState,
    RunExecutionQueueEvidence,
    RunRecoveryAction,
    RunRecoveryCleanupState,
    RunRecoveryGap,
    RunResultIndexingState,
)
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.services.run_recovery import (
    RUN_RECOVERY_FAIL_EVENT,
    RUN_RECOVERY_FAIL_REASON_CODE,
    RUN_RECOVERY_REQUEUE_CONFIRMED_EVENT,
    RUN_RECOVERY_REQUEUE_REQUESTED_EVENT,
    RunRecoveryError,
    RunRecoveryService,
)
from encode_pipeline.services.run_queue import (
    RunQueueError,
    RunQueueIdentityError,
    RunQueueJobUnavailableError,
    RunQueueUnavailableError,
)
from encode_pipeline.services.run_repositories import (
    ConcurrentRunUpdateError,
    InMemoryRunRepository,
    RunEventDraft,
)


NOW = datetime(2026, 8, 9, 12, 0, tzinfo=timezone.utc)


class _Queue:
    def __init__(
        self,
        state: ExecutionQueueEvidenceState,
        *,
        requeue_delivery_matches_request: bool = False,
    ) -> None:
        self.state = state
        self.requeue_delivery_matches_request = requeue_delivery_matches_request
        self.inspected: list[RunExecutionAssignment] = []
        self.requeued: list[RunExecutionAssignment] = []

    def inspect_execution(
        self,
        assignment: RunExecutionAssignment,
    ) -> RunExecutionQueueEvidence:
        self.inspected.append(assignment)
        return RunExecutionQueueEvidence(
            state=self.state,
            requeue_delivery_matches_request=(self.requeue_delivery_matches_request),
        )

    def requeue_execution(self, assignment: RunExecutionAssignment) -> str:
        self.requeued.append(assignment)
        self.state = ExecutionQueueEvidenceState.QUEUED
        self.requeue_delivery_matches_request = True
        return assignment.job_id


class _Cleanup:
    def __init__(self, succeeds: bool = True) -> None:
        self.succeeds = succeeds
        self.run_ids: list[str] = []

    def __call__(self, run_id: str) -> bool:
        self.run_ids.append(run_id)
        return self.succeeds


class _CallbackQueue(_Queue):
    def __init__(
        self,
        state: ExecutionQueueEvidenceState,
        *,
        on_second_inspection=None,
        on_requeue=None,
    ) -> None:
        super().__init__(state)
        self._on_second_inspection = on_second_inspection
        self._on_requeue = on_requeue

    def inspect_execution(
        self,
        assignment: RunExecutionAssignment,
    ) -> RunExecutionQueueEvidence:
        evidence = super().inspect_execution(assignment)
        if len(self.inspected) == 2 and self._on_second_inspection is not None:
            self._on_second_inspection()
        return evidence

    def requeue_execution(self, assignment: RunExecutionAssignment) -> str:
        backend_job_id = super().requeue_execution(assignment)
        if self._on_requeue is not None:
            self._on_requeue()
        return backend_job_id


class _UnavailableMutationQueue(_Queue):
    def requeue_execution(self, assignment: RunExecutionAssignment) -> str:
        self.requeued.append(assignment)
        raise RunQueueUnavailableError("redis://secret@private-host:6379")


class _MutationOutcomeQueue(_Queue):
    def __init__(
        self,
        state: ExecutionQueueEvidenceState,
        *,
        error: Exception | None = None,
        job_id: str | None = None,
    ) -> None:
        super().__init__(state)
        self._error = error
        self._job_id = job_id

    def requeue_execution(self, assignment: RunExecutionAssignment) -> str:
        self.requeued.append(assignment)
        if self._error is not None:
            raise self._error
        return self._job_id if self._job_id is not None else assignment.job_id


class _SequencedQueue(_Queue):
    def __init__(self, states: list[ExecutionQueueEvidenceState]) -> None:
        super().__init__(states[0])
        self._states = states

    def inspect_execution(
        self,
        assignment: RunExecutionAssignment,
    ) -> RunExecutionQueueEvidence:
        self.inspected.append(assignment)
        state = self._states[min(len(self.inspected) - 1, len(self._states) - 1)]
        return RunExecutionQueueEvidence(state=state)


class _InvalidQueue(_Queue):
    def inspect_execution(self, assignment: RunExecutionAssignment):
        self.inspected.append(assignment)
        return {"private_path": "/private/queue/state"}


class _ExplodingQueue(_Queue):
    def inspect_execution(self, assignment: RunExecutionAssignment):
        self.inspected.append(assignment)
        raise RuntimeError("redis://secret@private-host:6379/private")


class _RepositoryView:
    def __init__(self, repository, *, result_state=..., identity=...) -> None:
        self._repository = repository
        self._result_state = result_state
        self._identity = identity

    def __getattr__(self, name):
        return getattr(self._repository, name)

    def get_result_state(self, run_id: str):
        if self._result_state is ...:
            return self._repository.get_result_state(run_id)
        return self._result_state

    def get_workflow_build_identity(self, run_id: str):
        if self._identity is ...:
            return self._repository.get_workflow_build_identity(run_id)
        return self._identity


class _RepositoryMethodOverride:
    def __init__(
        self,
        repository,
        method_name: str,
        *,
        result=...,
        error: Exception | None = None,
    ) -> None:
        self._repository = repository
        self._method_name = method_name
        self._result = result
        self._error = error

    def __getattr__(self, name):
        if name != self._method_name:
            return getattr(self._repository, name)

        def override(*_args, **_kwargs):
            if self._error is not None:
                raise self._error
            return self._result

        return override


def _result_state_for(classification: RunResultIndexingState) -> RunResultState:
    artifact_generation = "artifactgen-" + "1" * 64
    artifact_fields = {
        "artifact_revision": 1,
        "artifact_generation": artifact_generation,
        "artifact_manifest_digest": "2" * 64,
        "artifact_outcome": "succeeded",
    }
    if classification is RunResultIndexingState.ARTIFACT_FAILED:
        return RunResultState(
            run_id="run-1",
            artifact_attempt_id="resultattempt-" + "3" * 64,
            artifact_attempt_status="failed",
            artifact_outcome="failed",
            artifact_reason_code="ARTIFACT_INDEXING_FAILED",
        )
    if classification is RunResultIndexingState.INCONSISTENT:
        return RunResultState(
            run_id="run-1",
            artifact_outcome="failed",
            artifact_reason_code="ARTIFACT_STATE_INCONSISTENT",
        )
    if classification is RunResultIndexingState.ARTIFACT_INDEXED:
        return RunResultState(run_id="run-1", **artifact_fields)
    if classification is RunResultIndexingState.QC_INDEXING:
        return RunResultState(
            run_id="run-1",
            **artifact_fields,
            qc_attempt_id="resultattempt-" + "4" * 64,
            qc_attempt_status="pending",
            qc_attempt_artifact_generation=artifact_generation,
        )
    if classification is RunResultIndexingState.QC_FAILED:
        return RunResultState(
            run_id="run-1",
            **artifact_fields,
            qc_attempt_id="resultattempt-" + "5" * 64,
            qc_attempt_status="failed",
            qc_attempt_artifact_generation=artifact_generation,
        )
    if classification is RunResultIndexingState.QC_INVALIDATED:
        return RunResultState(
            run_id="run-1",
            **artifact_fields,
            qc_outcome="invalidated",
        )
    raise AssertionError(f"unsupported classification: {classification}")


def _mark_run_succeeded(repository: InMemoryRunRepository) -> None:
    current = repository.get_run("run-1")
    repository.update_run(
        replace(
            current,
            status=RunStatus.SUCCEEDED,
            updated_at=NOW,
            ended_at=NOW,
        ),
        expected_status=RunStatus.RUNNING,
        event=RunEventDraft(
            event_type="status_changed",
            message="Scientific execution succeeded.",
            status=RunStatus.SUCCEEDED,
        ),
    )


def test_requeue_preserves_assignment_identity_and_records_request_and_confirmation():
    repository, assignment = _queued_repository()
    queue = _Queue(ExecutionQueueEvidenceState.MISSING)
    service = RunRecoveryService(repository, queue, clock=lambda: NOW)

    result = service.requeue_run(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=assignment,
    )

    persisted = repository.get_execution_assignment("run-1")
    assert persisted is not None
    assert persisted.job_id == assignment.job_id
    assert persisted.created_at == assignment.created_at
    assert persisted.dispatched_at == assignment.dispatched_at
    assert persisted.claimed_at is None
    assert persisted.requeue_requested_at == NOW
    assert persisted.requeue_confirmed_at == NOW
    assert result.assignment == persisted
    assert result.reason_code == "RUN_REQUEUED_BY_ADMIN_RECOVERY"
    assert queue.requeued == [replace(assignment, requeue_requested_at=NOW)]
    assert [event.event_type for event in repository.list_events("run-1")[-2:]] == [
        RUN_RECOVERY_REQUEUE_REQUESTED_EVENT,
        RUN_RECOVERY_REQUEUE_CONFIRMED_EVENT,
    ]


def test_requeue_retry_confirms_a_prepared_request_without_duplicate_enqueue():
    repository, assignment = _queued_repository()
    repository.prepare_execution_requeue(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=assignment,
        requested_at=NOW,
        event=RunEventDraft(
            event_type=RUN_RECOVERY_REQUEUE_REQUESTED_EVENT,
            message="Run requeue requested by an administrator.",
            status=RunStatus.QUEUED,
        ),
    )
    prepared = repository.get_execution_assignment("run-1")
    assert prepared is not None
    queue = _Queue(
        ExecutionQueueEvidenceState.QUEUED,
        requeue_delivery_matches_request=True,
    )
    service = RunRecoveryService(repository, queue, clock=lambda: NOW)

    result = service.requeue_run(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=prepared,
    )

    assert result.assignment.requeue_confirmed_at == NOW
    assert queue.requeued == []
    assert [event.event_type for event in repository.list_events("run-1")].count(
        RUN_RECOVERY_REQUEUE_REQUESTED_EVENT
    ) == 1
    assert [event.event_type for event in repository.list_events("run-1")].count(
        RUN_RECOVERY_REQUEUE_CONFIRMED_EVENT
    ) == 1


def test_requeue_retry_binds_an_active_original_job_to_the_pending_request():
    repository, assignment = _queued_repository()
    prepared = repository.prepare_execution_requeue(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=assignment,
        requested_at=NOW,
        event=RunEventDraft(
            event_type=RUN_RECOVERY_REQUEUE_REQUESTED_EVENT,
            message="Run requeue requested by an administrator.",
            status=RunStatus.QUEUED,
        ),
    ).assignment
    queue = _Queue(ExecutionQueueEvidenceState.QUEUED)
    service = RunRecoveryService(repository, queue, clock=lambda: NOW)

    result = service.requeue_run(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=prepared,
    )

    assert queue.requeued == [prepared]
    assert result.assignment.requeue_confirmed_at == NOW


@pytest.mark.parametrize(
    "queue_state",
    [
        ExecutionQueueEvidenceState.MISSING,
        ExecutionQueueEvidenceState.FINISHED,
        ExecutionQueueEvidenceState.FAILED,
        ExecutionQueueEvidenceState.STOPPED,
        ExecutionQueueEvidenceState.CANCELED,
    ],
)
def test_requeue_retry_uses_the_same_unclaimed_identity(queue_state):
    repository, assignment = _queued_repository()
    repository.prepare_execution_requeue(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=assignment,
        requested_at=NOW,
        event=RunEventDraft(
            event_type=RUN_RECOVERY_REQUEUE_REQUESTED_EVENT,
            message="Run requeue requested by an administrator.",
            status=RunStatus.QUEUED,
        ),
    )
    prepared = repository.get_execution_assignment("run-1")
    assert prepared is not None
    queue = _Queue(queue_state)
    service = RunRecoveryService(repository, queue, clock=lambda: NOW)

    diagnostic = service.diagnose("run-1")
    assert diagnostic.can_requeue is True

    result = service.requeue_run(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=prepared,
    )

    assert queue.requeued == [prepared]
    assert result.assignment.requeue_confirmed_at == NOW


def test_requeue_confirmation_preserves_a_worker_claim_and_status_race():
    repository, assignment = _queued_repository()

    def worker_claims_and_starts() -> None:
        repository.claim_execution_assignment(
            "run-1",
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
            claimed_at=NOW,
            allowed_statuses=frozenset({RunStatus.QUEUED}),
            event=RunEventDraft(
                event_type="worker_dependencies_rebuilt",
                message="Worker claimed the exact execution.",
                status=RunStatus.QUEUED,
            ),
        )
        current = repository.get_run("run-1")
        repository.update_run(
            replace(
                current,
                status=RunStatus.RUNNING,
                updated_at=NOW,
                started_at=NOW,
            ),
            expected_status=RunStatus.QUEUED,
            event=RunEventDraft(
                event_type="status_changed",
                message="Run started.",
                status=RunStatus.RUNNING,
            ),
        )

    queue = _CallbackQueue(
        ExecutionQueueEvidenceState.MISSING,
        on_requeue=worker_claims_and_starts,
    )
    service = RunRecoveryService(repository, queue, clock=lambda: NOW)

    result = service.requeue_run(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=assignment,
    )

    persisted = repository.get_execution_assignment("run-1")
    assert persisted is not None
    assert persisted.claimed_at == NOW
    assert persisted.requeue_requested_at == NOW
    assert persisted.requeue_confirmed_at == NOW
    assert repository.get_run("run-1").status is RunStatus.RUNNING
    assert result.status is RunStatus.RUNNING
    assert queue.requeued == [replace(assignment, requeue_requested_at=NOW)]


@pytest.mark.parametrize(
    (
        "run_kind",
        "queue_state",
        "expected_code",
        "expected_gaps",
        "expected_actions",
        "expected_cleanup",
    ),
    [
        (
            "queued",
            ExecutionQueueEvidenceState.QUEUED,
            "RUN_RECOVERY_OWNER_CONFIRMED",
            set(),
            set(),
            RunRecoveryCleanupState.NOT_REQUIRED,
        ),
        (
            "running",
            ExecutionQueueEvidenceState.STARTED_LIVE,
            "RUN_RECOVERY_OWNER_CONFIRMED",
            set(),
            set(),
            RunRecoveryCleanupState.NOT_REQUIRED,
        ),
        (
            "queued",
            ExecutionQueueEvidenceState.MISSING,
            "RUN_RECOVERY_JOB_MISSING_BEFORE_CLAIM",
            {RunRecoveryGap.QUEUE},
            {RunRecoveryAction.FAIL, RunRecoveryAction.REQUEUE},
            RunRecoveryCleanupState.NOT_REQUIRED,
        ),
        (
            "running",
            ExecutionQueueEvidenceState.MISSING,
            "RUN_RECOVERY_OWNER_UNPROVEN",
            {RunRecoveryGap.QUEUE, RunRecoveryGap.CLEANUP},
            set(),
            RunRecoveryCleanupState.BLOCKED,
        ),
        (
            "running",
            ExecutionQueueEvidenceState.STARTED_UNPROVEN,
            "RUN_RECOVERY_OWNER_UNPROVEN",
            {RunRecoveryGap.QUEUE, RunRecoveryGap.CLEANUP},
            set(),
            RunRecoveryCleanupState.BLOCKED,
        ),
        (
            "running",
            ExecutionQueueEvidenceState.IDENTITY_DRIFT,
            "RUN_RECOVERY_IDENTITY_DRIFT",
            {RunRecoveryGap.QUEUE, RunRecoveryGap.CLEANUP},
            set(),
            RunRecoveryCleanupState.BLOCKED,
        ),
        (
            "running",
            ExecutionQueueEvidenceState.UNAVAILABLE,
            "RUN_RECOVERY_QUEUE_UNAVAILABLE",
            {RunRecoveryGap.QUEUE, RunRecoveryGap.CLEANUP},
            set(),
            RunRecoveryCleanupState.BLOCKED,
        ),
        (
            "running",
            ExecutionQueueEvidenceState.FAILED,
            "RUN_RECOVERY_TERMINAL_CALLBACK_PENDING",
            {RunRecoveryGap.CALLBACK, RunRecoveryGap.CLEANUP},
            {RunRecoveryAction.FAIL},
            RunRecoveryCleanupState.PENDING,
        ),
    ],
)
def test_diagnosis_distinguishes_owner_and_coordination_states(
    run_kind,
    queue_state,
    expected_code,
    expected_gaps,
    expected_actions,
    expected_cleanup,
):
    repository, _assignment = (
        _queued_repository() if run_kind == "queued" else _running_repository()
    )
    cleanup = _Cleanup()
    service = RunRecoveryService(
        repository,
        _Queue(queue_state),
        cleanup=cleanup,
        clock=lambda: NOW,
    )

    diagnostic = service.diagnose("run-1")

    assert diagnostic.diagnosis_code == expected_code
    assert set(diagnostic.gaps) == expected_gaps
    assert set(diagnostic.allowed_actions) == expected_actions
    assert diagnostic.cleanup is expected_cleanup
    assert cleanup.run_ids == []


@pytest.mark.parametrize(
    (
        "queue_state",
        "expected_code",
        "expected_gaps",
        "expected_actions",
        "expected_cleanup",
    ),
    [
        (
            ExecutionQueueEvidenceState.STARTED_LIVE,
            "RUN_RECOVERY_OWNER_CONFIRMED",
            set(),
            set(),
            RunRecoveryCleanupState.NOT_REQUIRED,
        ),
        (
            ExecutionQueueEvidenceState.FAILED,
            "RUN_RECOVERY_TERMINAL_CALLBACK_PENDING",
            {RunRecoveryGap.CALLBACK, RunRecoveryGap.CLEANUP},
            {RunRecoveryAction.FAIL},
            RunRecoveryCleanupState.PENDING,
        ),
        (
            ExecutionQueueEvidenceState.MISSING,
            "RUN_RECOVERY_OWNER_UNPROVEN",
            {RunRecoveryGap.QUEUE, RunRecoveryGap.CLEANUP},
            set(),
            RunRecoveryCleanupState.BLOCKED,
        ),
    ],
)
def test_queued_claim_is_a_real_transient_owner_not_database_corruption(
    queue_state,
    expected_code,
    expected_gaps,
    expected_actions,
    expected_cleanup,
):
    repository, assignment = _queued_repository()
    claimed = repository.claim_execution_assignment(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
        claimed_at=NOW,
        allowed_statuses=frozenset({RunStatus.QUEUED}),
        event=RunEventDraft(
            event_type="worker_dependencies_rebuilt",
            message="Worker claimed the exact execution.",
            status=RunStatus.QUEUED,
        ),
    ).assignment
    service = RunRecoveryService(
        repository,
        _Queue(queue_state),
        cleanup=_Cleanup(),
        clock=lambda: NOW,
    )

    diagnostic = service.diagnose("run-1")

    assert diagnostic.assignment == claimed
    assert RunRecoveryGap.DATABASE not in diagnostic.gaps
    assert diagnostic.diagnosis_code == expected_code
    assert set(diagnostic.gaps) == expected_gaps
    assert set(diagnostic.allowed_actions) == expected_actions
    assert diagnostic.cleanup is expected_cleanup


@pytest.mark.parametrize("run_kind", ["queued", "running"])
@pytest.mark.parametrize(
    "queue_state",
    [
        ExecutionQueueEvidenceState.QUEUED,
        ExecutionQueueEvidenceState.DEFERRED,
        ExecutionQueueEvidenceState.SCHEDULED,
    ],
)
def test_claimed_assignment_with_schedulable_job_reports_status_queue_mismatch(
    run_kind,
    queue_state,
):
    if run_kind == "queued":
        repository, assignment = _queued_repository()
        repository.claim_execution_assignment(
            "run-1",
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
            claimed_at=NOW,
            allowed_statuses=frozenset({RunStatus.QUEUED}),
            event=RunEventDraft(
                event_type="worker_dependencies_rebuilt",
                message="Worker claimed the exact execution.",
                status=RunStatus.QUEUED,
            ),
        )
    else:
        repository, _assignment = _running_repository()
    service = RunRecoveryService(
        repository,
        _Queue(queue_state),
        cleanup=_Cleanup(),
        clock=lambda: NOW,
    )

    diagnostic = service.diagnose("run-1")

    assert diagnostic.diagnosis_code == "RUN_RECOVERY_STATUS_QUEUE_MISMATCH"
    assert RunRecoveryGap.QUEUE in diagnostic.gaps
    assert diagnostic.allowed_actions == ()


def test_terminal_worker_race_can_close_only_the_pending_requeue_audit():
    repository, assignment = _queued_repository()
    prepared = repository.prepare_execution_requeue(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=assignment,
        requested_at=NOW,
        event=RunEventDraft(
            event_type=RUN_RECOVERY_REQUEUE_REQUESTED_EVENT,
            message="Run requeue requested by an administrator.",
            status=RunStatus.QUEUED,
        ),
    ).assignment
    repository.claim_execution_assignment(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
        claimed_at=NOW,
        allowed_statuses=frozenset({RunStatus.QUEUED}),
        event=RunEventDraft(
            event_type="worker_dependencies_rebuilt",
            message="Worker claimed the exact execution.",
            status=RunStatus.QUEUED,
        ),
    )
    current = repository.get_run("run-1")
    repository.update_run(
        replace(
            current,
            status=RunStatus.RUNNING,
            updated_at=NOW,
            started_at=NOW,
        ),
        expected_status=RunStatus.QUEUED,
        event=RunEventDraft(
            event_type="status_changed",
            message="Run started.",
            status=RunStatus.RUNNING,
        ),
    )
    current = repository.get_run("run-1")
    repository.update_run(
        replace(
            current,
            status=RunStatus.SUCCEEDED,
            updated_at=NOW,
            ended_at=NOW,
        ),
        expected_status=RunStatus.RUNNING,
        event=RunEventDraft(
            event_type="status_changed",
            message="Run succeeded.",
            status=RunStatus.SUCCEEDED,
        ),
    )
    pending = repository.get_execution_assignment("run-1")
    assert pending is not None
    assert pending != prepared
    queue = _Queue(ExecutionQueueEvidenceState.FINISHED)
    service = RunRecoveryService(repository, queue, clock=lambda: NOW)

    diagnostic = service.diagnose("run-1")
    assert diagnostic.diagnosis_code == "RUN_RECOVERY_REQUEUE_CONFIRMATION_PENDING"
    assert diagnostic.allowed_actions == (RunRecoveryAction.REQUEUE,)

    result = service.requeue_run(
        "run-1",
        expected_status=RunStatus.SUCCEEDED,
        expected_assignment=pending,
    )

    assert result.status is RunStatus.SUCCEEDED
    assert result.assignment.requeue_confirmed_at == NOW
    assert repository.get_run("run-1").status is RunStatus.SUCCEEDED
    assert queue.requeued == []


def test_claimed_pending_requeue_confirms_after_active_queue_evidence_expires():
    repository, assignment = _queued_repository()
    prepared = repository.prepare_execution_requeue(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=assignment,
        requested_at=NOW,
        event=RunEventDraft(
            event_type=RUN_RECOVERY_REQUEUE_REQUESTED_EVENT,
            message="Run requeue requested by an administrator.",
            status=RunStatus.QUEUED,
        ),
    ).assignment
    claimed = repository.claim_execution_assignment(
        "run-1",
        job_id=prepared.job_id,
        backend=prepared.backend,
        queue_name=prepared.queue_name,
        claimed_at=NOW,
        allowed_statuses=frozenset({RunStatus.QUEUED}),
        event=RunEventDraft(
            event_type="worker_dependencies_rebuilt",
            message="Worker claimed the exact execution.",
            status=RunStatus.QUEUED,
        ),
    ).assignment
    queue = _Queue(ExecutionQueueEvidenceState.MISSING)
    service = RunRecoveryService(repository, queue, clock=lambda: NOW)

    diagnostic = service.diagnose("run-1")
    assert diagnostic.can_requeue is True
    assert RunRecoveryGap.QUEUE in diagnostic.gaps

    result = service.requeue_run(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=claimed,
    )

    assert result.status is RunStatus.QUEUED
    assert result.assignment.requeue_confirmed_at == NOW
    assert queue.requeued == []


@pytest.mark.parametrize(
    "queue_state",
    [ExecutionQueueEvidenceState.MISSING, ExecutionQueueEvidenceState.UNAVAILABLE],
)
def test_claimed_terminal_pending_requeue_confirms_without_redis_proof(queue_state):
    repository, assignment = _queued_repository()
    prepared = repository.prepare_execution_requeue(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=assignment,
        requested_at=NOW,
        event=RunEventDraft(
            event_type=RUN_RECOVERY_REQUEUE_REQUESTED_EVENT,
            message="Run requeue requested by an administrator.",
            status=RunStatus.QUEUED,
        ),
    ).assignment
    claimed = repository.claim_execution_assignment(
        "run-1",
        job_id=prepared.job_id,
        backend=prepared.backend,
        queue_name=prepared.queue_name,
        claimed_at=NOW,
        allowed_statuses=frozenset({RunStatus.QUEUED}),
        event=RunEventDraft(
            event_type="worker_dependencies_rebuilt",
            message="Worker claimed the exact execution.",
            status=RunStatus.QUEUED,
        ),
    ).assignment
    current = repository.get_run("run-1")
    repository.update_run(
        replace(
            current,
            status=RunStatus.FAILED,
            updated_at=NOW,
            ended_at=NOW,
            error=Issue(code="RUN_WORKER_FAILED", message="Worker failed."),
        ),
        expected_status=RunStatus.QUEUED,
        event=RunEventDraft(
            event_type="status_changed",
            message="Run failed.",
            status=RunStatus.FAILED,
        ),
    )
    queue = _Queue(queue_state)
    service = RunRecoveryService(repository, queue, clock=lambda: NOW)

    diagnostic = service.diagnose("run-1")
    assert diagnostic.can_requeue is True

    result = service.requeue_run(
        "run-1",
        expected_status=RunStatus.FAILED,
        expected_assignment=claimed,
    )

    assert result.status is RunStatus.FAILED
    assert result.assignment.requeue_confirmed_at == NOW
    assert queue.requeued == []


@pytest.mark.parametrize("delivery_matches", [False, True])
def test_unclaimed_terminal_pending_requeue_requires_request_bound_delivery_marker(
    delivery_matches,
):
    repository, assignment = _queued_repository()
    prepared = repository.prepare_execution_requeue(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=assignment,
        requested_at=NOW,
        event=RunEventDraft(
            event_type=RUN_RECOVERY_REQUEUE_REQUESTED_EVENT,
            message="Run requeue requested by an administrator.",
            status=RunStatus.QUEUED,
        ),
    ).assignment
    current = repository.get_run("run-1")
    repository.update_run(
        replace(
            current,
            status=RunStatus.FAILED,
            updated_at=NOW,
            ended_at=NOW,
            error=Issue(code="RUN_WORKER_FAILED", message="Worker failed."),
        ),
        expected_status=RunStatus.QUEUED,
        event=RunEventDraft(
            event_type="status_changed",
            message="Run failed.",
            status=RunStatus.FAILED,
        ),
    )
    queue = _Queue(
        ExecutionQueueEvidenceState.FAILED,
        requeue_delivery_matches_request=delivery_matches,
    )
    service = RunRecoveryService(repository, queue, clock=lambda: NOW)

    diagnostic = service.diagnose("run-1")
    assert diagnostic.can_requeue is delivery_matches

    if not delivery_matches:
        with pytest.raises(RunRecoveryError) as raised:
            service.requeue_run(
                "run-1",
                expected_status=RunStatus.FAILED,
                expected_assignment=prepared,
            )
        assert raised.value.code == "RUN_RECOVERY_NOT_SAFE"
        assert queue.requeued == []
        return

    result = service.requeue_run(
        "run-1",
        expected_status=RunStatus.FAILED,
        expected_assignment=prepared,
    )
    assert result.status is RunStatus.FAILED
    assert result.assignment.requeue_confirmed_at == NOW
    assert queue.requeued == []


def test_worker_confirmed_preclaim_failure_stays_closed_after_rq_ttl_expiry():
    repository, assignment = _queued_repository()
    prepared = repository.prepare_execution_requeue(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=assignment,
        requested_at=NOW,
        event=RunEventDraft(
            event_type=RUN_RECOVERY_REQUEUE_REQUESTED_EVENT,
            message="Run requeue requested by an administrator.",
            status=RunStatus.QUEUED,
        ),
    ).assignment
    confirmed = repository.confirm_execution_requeue(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=prepared,
        confirmed_at=NOW,
        event=RunEventDraft(
            event_type="run_requeue_delivery_observed_by_worker",
            message="Worker observed the exact administrator requeue.",
            status=None,
            stage="execution",
            context={
                "reason_code": "RUN_REQUEUED_BY_ADMIN_RECOVERY",
                "confirmation_source": "worker_entry",
            },
        ),
    )
    current = repository.get_run("run-1")
    repository.update_run(
        replace(
            current,
            status=RunStatus.FAILED,
            updated_at=NOW,
            ended_at=NOW,
            error=Issue(code="RUN_WORKER_FAILED", message="Worker failed."),
        ),
        expected_status=RunStatus.QUEUED,
        event=RunEventDraft(
            event_type="status_changed",
            message="Run failed before claim.",
            status=RunStatus.FAILED,
        ),
    )

    diagnostic = RunRecoveryService(
        repository,
        _Queue(ExecutionQueueEvidenceState.MISSING),
        clock=lambda: NOW,
    ).diagnose("run-1")

    assert confirmed.claimed_at is None
    assert confirmed.requeue_confirmed_at == NOW
    assert diagnostic.diagnosis_code == "RUN_RECOVERY_NOT_REQUIRED"
    assert diagnostic.gaps == ()
    assert diagnostic.allowed_actions == ()


@pytest.mark.parametrize("dispatched", [False, True])
def test_planned_assignment_reservation_is_normal_and_never_inspects_queue(
    dispatched,
):
    repository = InMemoryRunRepository()
    _create_execution_ready_run(repository, RunStatus.PLANNED)
    assignment = RunExecutionAssignment(
        run_id="run-1",
        job_id="job-1",
        backend="rq",
        queue_name="default",
        created_at=NOW,
        dispatched_at=NOW if dispatched else None,
    )
    repository.ensure_execution_assignment(
        assignment,
        expected_status=RunStatus.PLANNED,
    )
    queue = _Queue(ExecutionQueueEvidenceState.QUEUED)
    service = RunRecoveryService(repository, queue, clock=lambda: NOW)

    diagnostic = service.diagnose("run-1")

    assert diagnostic.diagnosis_code == "RUN_RECOVERY_NOT_REQUIRED"
    assert diagnostic.gaps == ()
    assert diagnostic.allowed_actions == ()
    assert queue.inspected == []


@pytest.mark.parametrize(
    ("queue_state", "expected_code"),
    [
        (
            ExecutionQueueEvidenceState.UNAVAILABLE,
            "RUN_RECOVERY_QUEUE_UNAVAILABLE",
        ),
        (ExecutionQueueEvidenceState.IDENTITY_DRIFT, "RUN_RECOVERY_NOT_SAFE"),
        (ExecutionQueueEvidenceState.QUEUED, "RUN_RECOVERY_NOT_SAFE"),
        (ExecutionQueueEvidenceState.STARTED_LIVE, "RUN_RECOVERY_NOT_SAFE"),
        (ExecutionQueueEvidenceState.STARTED_UNPROVEN, "RUN_RECOVERY_NOT_SAFE"),
        (ExecutionQueueEvidenceState.UNKNOWN, "RUN_RECOVERY_NOT_SAFE"),
    ],
)
def test_claimed_failure_refuses_unproven_or_live_queue_evidence(
    queue_state,
    expected_code,
):
    repository, assignment = _running_repository()
    cleanup = _Cleanup()
    service = RunRecoveryService(
        repository,
        _Queue(queue_state),
        cleanup=cleanup,
        clock=lambda: NOW,
    )

    with pytest.raises(RunRecoveryError) as error:
        service.fail_run(
            "run-1",
            expected_status=RunStatus.RUNNING,
            expected_assignment=assignment,
        )

    assert error.value.code == expected_code
    assert repository.get_run("run-1").status is RunStatus.RUNNING
    assert cleanup.run_ids == []


def test_claimed_terminal_failure_cleans_then_atomically_fails_without_result_writes():
    repository, assignment = _running_repository()
    cleanup = _Cleanup()
    initial_result_state = repository.get_result_state("run-1")
    service = RunRecoveryService(
        repository,
        _Queue(ExecutionQueueEvidenceState.FAILED),
        cleanup=cleanup,
        clock=lambda: NOW,
    )

    result = service.fail_run(
        "run-1",
        expected_status=RunStatus.RUNNING,
        expected_assignment=assignment,
    )

    record = repository.get_run("run-1")
    assert result.status is RunStatus.FAILED
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == RUN_RECOVERY_FAIL_REASON_CODE
    assert cleanup.run_ids == ["run-1"]
    assert repository.get_result_state("run-1") == initial_result_state
    event = repository.list_events("run-1")[-1]
    assert event.event_type == RUN_RECOVERY_FAIL_EVENT
    assert event.context == {
        "previous_status": "running",
        "new_status": "failed",
        "reason_code": RUN_RECOVERY_FAIL_REASON_CODE,
    }


def test_cancellation_intent_refuses_admin_failure_even_with_terminal_evidence():
    repository, assignment = _running_repository(cancellation_requested=True)
    service = RunRecoveryService(
        repository,
        _Queue(ExecutionQueueEvidenceState.STOPPED),
        cleanup=_Cleanup(),
        clock=lambda: NOW,
    )

    with pytest.raises(RunRecoveryError) as error:
        service.fail_run(
            "run-1",
            expected_status=RunStatus.RUNNING,
            expected_assignment=assignment,
        )

    assert error.value.code == "RUN_RECOVERY_NOT_SAFE"
    assert repository.get_run("run-1").status is RunStatus.RUNNING


def test_cancellation_race_wins_the_final_failure_cas():
    repository, assignment = _running_repository()

    def request_cancellation() -> None:
        repository.request_execution_cancellation(
            "run-1",
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
            requested_at=NOW,
            reason="private operator reason",
            event=RunEventDraft(
                event_type="cancellation_requested",
                message="Cancellation requested.",
                status=RunStatus.RUNNING,
            ),
        )

    queue = _CallbackQueue(
        ExecutionQueueEvidenceState.FAILED,
        on_second_inspection=request_cancellation,
    )
    cleanup = _Cleanup()
    service = RunRecoveryService(
        repository,
        queue,
        cleanup=cleanup,
        clock=lambda: NOW,
    )

    with pytest.raises(RunRecoveryError) as error:
        service.fail_run(
            "run-1",
            expected_status=RunStatus.RUNNING,
            expected_assignment=assignment,
        )

    assert error.value.code == "RUN_RECOVERY_CONFLICT"
    assert repository.get_run("run-1").status is RunStatus.RUNNING
    assert repository.get_execution_assignment("run-1").cancellation_requested_at == NOW
    assert RUN_RECOVERY_FAIL_EVENT not in {
        event.event_type for event in repository.list_events("run-1")
    }
    assert cleanup.run_ids == ["run-1"]


def test_terminal_callback_race_wins_before_failure_commit():
    repository, assignment = _running_repository()

    def complete_callback() -> None:
        current = repository.get_run("run-1")
        repository.update_run(
            replace(
                current,
                status=RunStatus.SUCCEEDED,
                updated_at=NOW,
                ended_at=NOW,
            ),
            expected_status=RunStatus.RUNNING,
            event=RunEventDraft(
                event_type="status_changed",
                message="Worker callback completed.",
                status=RunStatus.SUCCEEDED,
            ),
        )

    queue = _CallbackQueue(
        ExecutionQueueEvidenceState.FINISHED,
        on_second_inspection=complete_callback,
    )
    service = RunRecoveryService(
        repository,
        queue,
        cleanup=_Cleanup(),
        clock=lambda: NOW,
    )

    with pytest.raises(RunRecoveryError) as error:
        service.fail_run(
            "run-1",
            expected_status=RunStatus.RUNNING,
            expected_assignment=assignment,
        )

    assert error.value.code == "RUN_RECOVERY_CONFLICT"
    assert repository.get_run("run-1").status is RunStatus.SUCCEEDED
    assert RUN_RECOVERY_FAIL_EVENT not in {
        event.event_type for event in repository.list_events("run-1")
    }


def test_cleanup_failure_refuses_claimed_failure_without_lifecycle_mutation():
    repository, assignment = _running_repository()
    cleanup = _Cleanup(succeeds=False)
    service = RunRecoveryService(
        repository,
        _Queue(ExecutionQueueEvidenceState.STOPPED),
        cleanup=cleanup,
        clock=lambda: NOW,
    )

    with pytest.raises(RunRecoveryError) as error:
        service.fail_run(
            "run-1",
            expected_status=RunStatus.RUNNING,
            expected_assignment=assignment,
        )

    assert error.value.code == "RUN_RECOVERY_CLEANUP_FAILED"
    assert repository.get_run("run-1").status is RunStatus.RUNNING
    assert RUN_RECOVERY_FAIL_EVENT not in {
        event.event_type for event in repository.list_events("run-1")
    }


def test_queue_failure_after_prepare_keeps_auditable_retryable_intent():
    repository, assignment = _queued_repository()
    queue = _UnavailableMutationQueue(ExecutionQueueEvidenceState.MISSING)
    service = RunRecoveryService(repository, queue, clock=lambda: NOW)

    with pytest.raises(RunRecoveryError) as error:
        service.requeue_run(
            "run-1",
            expected_status=RunStatus.QUEUED,
            expected_assignment=assignment,
        )

    assert error.value.code == "RUN_RECOVERY_QUEUE_UNAVAILABLE"
    assert "private-host" not in str(error.value)
    prepared = repository.get_execution_assignment("run-1")
    assert prepared is not None
    assert prepared.requeue_requested_at == NOW
    assert prepared.requeue_confirmed_at is None
    assert [event.event_type for event in repository.list_events("run-1")].count(
        RUN_RECOVERY_REQUEUE_REQUESTED_EVENT
    ) == 1
    assert [event.event_type for event in repository.list_events("run-1")].count(
        RUN_RECOVERY_REQUEUE_CONFIRMED_EVENT
    ) == 0

    retry_queue = _Queue(ExecutionQueueEvidenceState.MISSING)
    retried = RunRecoveryService(
        repository,
        retry_queue,
        clock=lambda: NOW,
    ).requeue_run(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=prepared,
    )
    assert retried.assignment.requeue_confirmed_at == NOW
    assert retry_queue.requeued == [prepared]


def test_admin_fail_intentionally_closes_a_pending_requeue_handshake():
    repository, assignment = _queued_repository()
    prepared = repository.prepare_execution_requeue(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=assignment,
        requested_at=NOW,
        event=RunEventDraft(
            event_type=RUN_RECOVERY_REQUEUE_REQUESTED_EVENT,
            message="Run requeue requested by an administrator.",
            status=RunStatus.QUEUED,
        ),
    ).assignment
    queue = _Queue(ExecutionQueueEvidenceState.MISSING)
    service = RunRecoveryService(repository, queue, clock=lambda: NOW)

    failed = service.fail_run(
        "run-1",
        expected_status=RunStatus.QUEUED,
        expected_assignment=prepared,
    )
    diagnostic = service.diagnose("run-1")

    assert failed.status is RunStatus.FAILED
    assert failed.assignment.requeue_confirmed_at is None
    assert diagnostic.diagnosis_code == "RUN_RECOVERY_NOT_REQUIRED"
    assert diagnostic.gaps == ()
    assert diagnostic.allowed_actions == ()
    assert RUN_RECOVERY_FAIL_EVENT in {
        event.event_type for event in repository.list_events("run-1")
    }


def test_natural_terminal_winner_does_not_leave_a_false_cancellation_gap():
    repository, assignment = _running_repository()
    repository.request_execution_cancellation(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
        requested_at=NOW,
        reason="operator request",
        event=RunEventDraft(
            event_type="cancellation_requested",
            message="Cancellation requested.",
            status=RunStatus.RUNNING,
        ),
    )
    current = repository.get_run("run-1")
    repository.update_run(
        replace(
            current,
            status=RunStatus.FAILED,
            updated_at=NOW,
            ended_at=NOW,
            error=Issue(
                code="RUN_WORKER_FAILED",
                message="Run failed before cancellation acknowledgement.",
            ),
        ),
        expected_status=RunStatus.RUNNING,
        event=RunEventDraft(
            event_type="status_changed",
            message="Natural terminal state won.",
            status=RunStatus.FAILED,
        ),
    )
    service = RunRecoveryService(
        repository,
        _Queue(ExecutionQueueEvidenceState.STOPPED),
        cleanup=_Cleanup(),
        clock=lambda: NOW,
    )

    diagnostic = service.diagnose("run-1")

    assert diagnostic.diagnosis_code == "RUN_RECOVERY_NOT_REQUIRED"
    assert diagnostic.gaps == ()
    assert "RUN_RECOVERY_CANCELLATION_ACK_PENDING" not in diagnostic.issue_codes
    assert diagnostic.allowed_actions == ()


def test_terminal_cleanup_failure_remains_visible_without_a_mutation_action():
    repository, _assignment = _running_repository()
    repository.add_event(
        "run-1",
        RunEventDraft(
            event_type="execution_cleanup_failed",
            message="Managed workflow cleanup could not be confirmed.",
            status=None,
            stage="execution",
            context={"reason_code": "MANAGED_CONTAINER_CLEANUP_FAILED"},
            issue=Issue(
                code="MANAGED_CONTAINER_CLEANUP_FAILED",
                message="Managed workflow cleanup could not be confirmed.",
            ),
        ),
    )
    current = repository.get_run("run-1")
    repository.update_run(
        replace(
            current,
            status=RunStatus.FAILED,
            updated_at=NOW,
            ended_at=NOW,
            error=Issue(code="RUN_WORKER_FAILED", message="Worker failed."),
        ),
        expected_status=RunStatus.RUNNING,
        event=RunEventDraft(
            event_type="status_changed",
            message="Run failed.",
            status=RunStatus.FAILED,
        ),
    )
    service = RunRecoveryService(
        repository,
        _Queue(ExecutionQueueEvidenceState.FAILED),
        cleanup=_Cleanup(),
        clock=lambda: NOW,
    )

    diagnostic = service.diagnose("run-1")

    assert diagnostic.diagnosis_code == "RUN_RECOVERY_CLEANUP_NOT_CONFIRMED"
    assert diagnostic.cleanup is RunRecoveryCleanupState.BLOCKED
    assert diagnostic.gaps == (RunRecoveryGap.CLEANUP,)
    assert diagnostic.allowed_actions == ()


@pytest.mark.parametrize(
    ("queue_state", "expected_code", "expected_gaps"),
    [
        (
            ExecutionQueueEvidenceState.STARTED_LIVE,
            "RUN_RECOVERY_NOT_REQUIRED",
            set(),
        ),
        (
            ExecutionQueueEvidenceState.FINISHED,
            "RUN_RECOVERY_ARTIFACT_INDEXING",
            {RunRecoveryGap.CALLBACK, RunRecoveryGap.RESULT_INDEXING},
        ),
        (
            ExecutionQueueEvidenceState.MISSING,
            "RUN_RECOVERY_ARTIFACT_INDEXING",
            {RunRecoveryGap.QUEUE, RunRecoveryGap.RESULT_INDEXING},
        ),
        (
            ExecutionQueueEvidenceState.UNAVAILABLE,
            "RUN_RECOVERY_QUEUE_UNAVAILABLE",
            {RunRecoveryGap.QUEUE, RunRecoveryGap.RESULT_INDEXING},
        ),
    ],
)
def test_succeeded_result_indexing_distinguishes_live_worker_from_lost_owner(
    queue_state,
    expected_code,
    expected_gaps,
):
    repository, _assignment = _running_repository()
    current = repository.get_run("run-1")
    repository.update_run(
        replace(
            current,
            status=RunStatus.SUCCEEDED,
            updated_at=NOW,
            ended_at=NOW,
        ),
        expected_status=RunStatus.RUNNING,
        event=RunEventDraft(
            event_type="status_changed",
            message="Scientific execution succeeded.",
            status=RunStatus.SUCCEEDED,
        ),
    )
    repository.begin_artifact_result_attempt(
        "run-1",
        attempt_id="resultattempt-" + "1" * 64,
        expected_status=RunStatus.SUCCEEDED,
    )
    service = RunRecoveryService(
        repository,
        _Queue(queue_state),
        cleanup=_Cleanup(),
        clock=lambda: NOW,
    )

    diagnostic = service.diagnose("run-1")

    assert diagnostic.result_indexing.state is RunResultIndexingState.ARTIFACT_INDEXING
    assert diagnostic.diagnosis_code == expected_code
    assert set(diagnostic.gaps) == expected_gaps
    assert diagnostic.allowed_actions == ()


def test_qc_indexed_success_never_requires_historical_queue_access():
    repository, _assignment = _running_repository()
    current = repository.get_run("run-1")
    repository.update_run(
        replace(
            current,
            status=RunStatus.SUCCEEDED,
            updated_at=NOW,
            ended_at=NOW,
        ),
        expected_status=RunStatus.RUNNING,
        event=RunEventDraft(
            event_type="status_changed",
            message="Run and result indexing succeeded.",
            status=RunStatus.SUCCEEDED,
        ),
    )
    artifact_generation = "artifactgen-" + "2" * 64
    result_state = RunResultState(
        run_id="run-1",
        artifact_revision=1,
        artifact_generation=artifact_generation,
        artifact_manifest_digest="3" * 64,
        artifact_attempt_id="resultattempt-" + "4" * 64,
        artifact_attempt_status="succeeded",
        artifact_outcome="succeeded",
        qc_revision=1,
        qc_generation="qcgen-" + "5" * 64,
        qc_manifest_digest="6" * 64,
        qc_attempt_id="resultattempt-" + "7" * 64,
        qc_attempt_status="succeeded",
        qc_attempt_artifact_generation=artifact_generation,
        qc_artifact_generation=artifact_generation,
        qc_outcome="succeeded",
    )
    queue = _ExplodingQueue(ExecutionQueueEvidenceState.UNAVAILABLE)
    service = RunRecoveryService(
        _RepositoryView(repository, result_state=result_state),
        queue,
        clock=lambda: NOW,
    )

    diagnostic = service.diagnose("run-1")

    assert diagnostic.result_indexing.state is RunResultIndexingState.QC_INDEXED
    assert diagnostic.diagnosis_code == "RUN_RECOVERY_NOT_REQUIRED"
    assert diagnostic.gaps == ()
    assert queue.inspected == []


@pytest.mark.parametrize(
    ("repository_option", "invalid_value"),
    [
        ("result_state", {"private_path": "/private/result/state"}),
        ("identity", {"digest": "/private/build/identity"}),
    ],
)
def test_repository_contract_corruption_fails_closed_through_public_diagnose(
    repository_option,
    invalid_value,
):
    repository, _assignment = _queued_repository()
    options = {repository_option: invalid_value}
    service = RunRecoveryService(
        _RepositoryView(repository, **options),
        _Queue(ExecutionQueueEvidenceState.QUEUED),
        clock=lambda: NOW,
    )

    with pytest.raises(RunRecoveryError) as error:
        service.diagnose("run-1")

    assert error.value.code == "RUN_RECOVERY_DATA_INVALID"
    assert "/private" not in str(error.value)


def test_invalid_queue_projection_fails_closed_through_public_diagnose():
    repository, _assignment = _queued_repository()
    service = RunRecoveryService(
        repository,
        _InvalidQueue(ExecutionQueueEvidenceState.UNKNOWN),
        clock=lambda: NOW,
    )

    with pytest.raises(RunRecoveryError) as error:
        service.diagnose("run-1")

    assert error.value.code == "RUN_RECOVERY_DATA_INVALID"
    assert "/private" not in str(error.value)


def test_queue_exception_and_private_run_inputs_are_redacted_from_diagnosis():
    repository, _assignment = _running_repository()
    service = RunRecoveryService(
        repository,
        _ExplodingQueue(ExecutionQueueEvidenceState.UNKNOWN),
        clock=lambda: NOW,
    )

    diagnostic = service.diagnose("run-1")
    payload = json.dumps(diagnostic.to_dict(), sort_keys=True)

    assert diagnostic.diagnosis_code == "RUN_RECOVERY_QUEUE_UNAVAILABLE"
    for secret in (
        "/tmp/private/sample.tsv",
        "private-host",
        "redis://",
        "a" * 64,
    ):
        assert secret not in payload


def test_diagnosis_reports_build_assignment_and_result_marker_inconsistencies():
    repository = InMemoryRunRepository()
    repository.create_run(
        _record(RunStatus.RUNNING),
        RunEventDraft(
            event_type="status_changed",
            message="Run created.",
            status=RunStatus.RUNNING,
        ),
    )
    service = RunRecoveryService(
        repository,
        _Queue(ExecutionQueueEvidenceState.UNKNOWN),
        clock=lambda: NOW,
    )

    diagnostic = service.diagnose("run-1")

    assert diagnostic.status is RunStatus.RUNNING
    assert diagnostic.assignment is None
    assert diagnostic.queue_evidence.state is ExecutionQueueEvidenceState.UNKNOWN
    assert diagnostic.build_identity_present is False
    assert set(diagnostic.issue_codes) >= {
        "RUN_RECOVERY_BUILD_IDENTITY_MISSING",
        "RUN_RECOVERY_ASSIGNMENT_MISSING",
        "RUN_RECOVERY_CLAIM_MARKER_MISSING",
    }
    assert "/" not in str(diagnostic.to_dict())


@pytest.mark.parametrize(
    "invalid_field",
    ["repository", "queue", "cleanup", "clock"],
)
def test_recovery_service_rejects_invalid_dependencies(invalid_field):
    repository, _assignment = _queued_repository()
    queue = _Queue(ExecutionQueueEvidenceState.QUEUED)
    kwargs = {"cleanup": None, "clock": None}
    if invalid_field == "repository":
        repository = None
    elif invalid_field == "queue":
        queue = None
    else:
        kwargs[invalid_field] = 1

    with pytest.raises(ValueError):
        RunRecoveryService(repository, queue, **kwargs)


def test_recovery_error_rejects_unknown_public_codes():
    with pytest.raises(ValueError, match="error code is invalid"):
        RunRecoveryError("PRIVATE_INTERNAL_FAILURE")


@pytest.mark.parametrize(
    ("case", "expected_code"),
    [
        ("empty_run_id", "RUN_RECOVERY_DATA_INVALID"),
        ("status_string", "RUN_RECOVERY_DATA_INVALID"),
        ("assignment_object", "RUN_RECOVERY_DATA_INVALID"),
        ("assignment_wrong_run", "RUN_RECOVERY_DATA_INVALID"),
    ],
)
def test_recovery_mutations_require_exact_typed_identity(case, expected_code):
    repository, assignment = _queued_repository()
    run_id = "run-1"
    expected_status = RunStatus.QUEUED
    expected_assignment = assignment
    if case == "empty_run_id":
        run_id = " "
    elif case == "status_string":
        expected_status = "queued"
    elif case == "assignment_object":
        expected_assignment = object()
    else:
        expected_assignment = replace(assignment, run_id="other-run")
    service = RunRecoveryService(
        repository,
        _Queue(ExecutionQueueEvidenceState.MISSING),
        clock=lambda: NOW,
    )

    with pytest.raises(RunRecoveryError) as error:
        service.requeue_run(
            run_id,
            expected_status=expected_status,
            expected_assignment=expected_assignment,
        )

    assert error.value.code == expected_code
    assert repository.get_execution_assignment("run-1") == assignment


@pytest.mark.parametrize(
    ("method_name", "mode", "expected_code"),
    [
        ("get_run", "key_error", "RUN_RECOVERY_NOT_FOUND"),
        ("get_run", "exception", "RUN_RECOVERY_DATA_INVALID"),
        ("get_run", "invalid_type", "RUN_RECOVERY_DATA_INVALID"),
        ("get_run", "wrong_run", "RUN_RECOVERY_DATA_INVALID"),
        ("get_execution_assignment", "exception", "RUN_RECOVERY_DATA_INVALID"),
        (
            "get_execution_assignment",
            "invalid_type",
            "RUN_RECOVERY_DATA_INVALID",
        ),
        ("get_execution_assignment", "wrong_run", "RUN_RECOVERY_DATA_INVALID"),
        ("get_workflow_build_identity", "exception", "RUN_RECOVERY_DATA_INVALID"),
        (
            "get_workflow_build_identity",
            "invalid_type",
            "RUN_RECOVERY_DATA_INVALID",
        ),
        ("get_result_state", "exception", "RUN_RECOVERY_DATA_INVALID"),
        ("get_result_state", "invalid_type", "RUN_RECOVERY_DATA_INVALID"),
        ("get_result_state", "wrong_run", "RUN_RECOVERY_DATA_INVALID"),
        ("list_events", "exception", "RUN_RECOVERY_DATA_INVALID"),
        ("list_events", "invalid_type", "RUN_RECOVERY_DATA_INVALID"),
        ("list_events", "invalid_member", "RUN_RECOVERY_DATA_INVALID"),
    ],
)
def test_diagnose_fails_closed_for_each_corrupted_repository_read(
    method_name,
    mode,
    expected_code,
):
    repository, assignment = _queued_repository()
    result = ...
    error = None
    if mode == "key_error":
        error = KeyError("/private/missing/run")
    elif mode == "exception":
        error = RuntimeError("/private/repository/failure")
    elif mode == "invalid_type":
        result = (
            list(repository.list_events("run-1"))
            if method_name == "list_events"
            else {"private_path": "/private/repository/value"}
        )
    elif mode == "invalid_member":
        result = (object(),)
    elif method_name == "get_run":
        result = replace(repository.get_run("run-1"), run_id="other-run")
    elif method_name == "get_execution_assignment":
        result = replace(assignment, run_id="other-run")
    else:
        result = RunResultState(run_id="other-run")
    service = RunRecoveryService(
        _RepositoryMethodOverride(
            repository,
            method_name,
            result=result,
            error=error,
        ),
        _Queue(ExecutionQueueEvidenceState.QUEUED),
        clock=lambda: NOW,
    )

    with pytest.raises(RunRecoveryError) as raised:
        service.diagnose("run-1")

    assert raised.value.code == expected_code
    assert "/private" not in str(raised.value)


def test_recovery_summary_reports_all_bounded_gap_categories_and_queue_outage():
    repository, _assignment = _queued_repository()
    summary = RunRecoveryService(
        repository,
        _Queue(ExecutionQueueEvidenceState.UNAVAILABLE),
        clock=lambda: NOW,
    ).summarize()

    assert summary.to_dict() == {
        "runs_examined": 1,
        "database_gaps": 0,
        "queue_gaps": 1,
        "callback_gaps": 0,
        "result_indexing_gaps": 0,
        "cleanup_gaps": 0,
        "queue_unavailable": True,
    }


def test_recovery_summary_counts_a_fail_closed_diagnosis_as_database_gap():
    repository, _assignment = _queued_repository()
    summary = RunRecoveryService(
        _RepositoryMethodOverride(
            repository,
            "get_run",
            error=RuntimeError("/private/corrupted/run"),
        ),
        _Queue(ExecutionQueueEvidenceState.QUEUED),
        clock=lambda: NOW,
    ).summarize()

    assert summary.runs_examined == 1
    assert summary.database_gaps == 1
    assert summary.queue_gaps == 0
    assert summary.queue_unavailable is False


@pytest.mark.parametrize("mode", ["exception", "invalid_type"])
def test_recovery_summary_fails_closed_when_run_inventory_is_not_authoritative(mode):
    repository, _assignment = _queued_repository()
    view = _RepositoryMethodOverride(
        repository,
        "list_runs",
        result=[] if mode == "invalid_type" else ...,
        error=(
            RuntimeError("/private/inventory/failure") if mode == "exception" else None
        ),
    )

    with pytest.raises(RunRecoveryError) as error:
        RunRecoveryService(
            view,
            _Queue(ExecutionQueueEvidenceState.QUEUED),
            clock=lambda: NOW,
        ).summarize()

    assert error.value.code == "RUN_RECOVERY_DATA_INVALID"
    assert "/private" not in str(error.value)


@pytest.mark.parametrize(
    "classification",
    [
        RunResultIndexingState.ARTIFACT_FAILED,
        RunResultIndexingState.ARTIFACT_INDEXED,
        RunResultIndexingState.QC_INDEXING,
        RunResultIndexingState.QC_FAILED,
        RunResultIndexingState.QC_INVALIDATED,
        RunResultIndexingState.INCONSISTENT,
    ],
)
def test_diagnose_projects_each_durable_result_indexing_outcome(classification):
    repository, _assignment = _running_repository()
    _mark_run_succeeded(repository)
    service = RunRecoveryService(
        _RepositoryView(
            repository,
            result_state=_result_state_for(classification),
        ),
        _Queue(ExecutionQueueEvidenceState.STARTED_LIVE),
        clock=lambda: NOW,
    )

    diagnostic = service.diagnose("run-1")

    assert diagnostic.result_indexing.state is classification
    if classification in {
        RunResultIndexingState.ARTIFACT_INDEXED,
        RunResultIndexingState.QC_INDEXING,
    }:
        assert RunRecoveryGap.RESULT_INDEXING not in diagnostic.gaps
    else:
        assert RunRecoveryGap.RESULT_INDEXING in diagnostic.gaps


def test_invalid_clock_refuses_requeue_before_recording_an_intent():
    repository, assignment = _queued_repository()
    service = RunRecoveryService(
        repository,
        _Queue(ExecutionQueueEvidenceState.MISSING),
        clock=lambda: "/private/not-a-time",
    )

    with pytest.raises(RunRecoveryError) as error:
        service.requeue_run(
            "run-1",
            expected_status=RunStatus.QUEUED,
            expected_assignment=assignment,
        )

    assert error.value.code == "RUN_RECOVERY_DATA_INVALID"
    assert repository.get_execution_assignment("run-1") == assignment


def test_cleanup_exception_refuses_failure_without_exposing_operator_details():
    repository, assignment = _running_repository()

    def fail_cleanup(_run_id: str) -> bool:
        raise RuntimeError("/private/container/runtime")

    service = RunRecoveryService(
        repository,
        _Queue(ExecutionQueueEvidenceState.STOPPED),
        cleanup=fail_cleanup,
        clock=lambda: NOW,
    )

    with pytest.raises(RunRecoveryError) as error:
        service.fail_run(
            "run-1",
            expected_status=RunStatus.RUNNING,
            expected_assignment=assignment,
        )

    assert error.value.code == "RUN_RECOVERY_CLEANUP_FAILED"
    assert "/private" not in str(error.value)
    assert repository.get_run("run-1").status is RunStatus.RUNNING


@pytest.mark.parametrize(
    ("final_state", "expected_code"),
    [
        (
            ExecutionQueueEvidenceState.UNAVAILABLE,
            "RUN_RECOVERY_QUEUE_UNAVAILABLE",
        ),
        (ExecutionQueueEvidenceState.STARTED_LIVE, "RUN_RECOVERY_NOT_SAFE"),
    ],
)
def test_final_queue_ownership_race_refuses_admin_failure(
    final_state,
    expected_code,
):
    repository, assignment = _queued_repository()
    queue = _SequencedQueue([ExecutionQueueEvidenceState.MISSING, final_state])
    service = RunRecoveryService(repository, queue, clock=lambda: NOW)

    with pytest.raises(RunRecoveryError) as error:
        service.fail_run(
            "run-1",
            expected_status=RunStatus.QUEUED,
            expected_assignment=assignment,
        )

    assert error.value.code == expected_code
    assert repository.get_run("run-1").status is RunStatus.QUEUED


@pytest.mark.parametrize(
    ("queue_error", "expected_code"),
    [
        (RunQueueIdentityError("/private/identity"), "RUN_RECOVERY_NOT_SAFE"),
        (
            RunQueueJobUnavailableError("/private/job"),
            "RUN_RECOVERY_NOT_SAFE",
        ),
        (RunQueueError("/private/queue"), "RUN_RECOVERY_NOT_SAFE"),
        (RuntimeError("/private/unexpected"), "RUN_RECOVERY_QUEUE_UNAVAILABLE"),
    ],
)
def test_requeue_maps_queue_failures_to_bounded_errors(queue_error, expected_code):
    repository, assignment = _queued_repository()
    queue = _MutationOutcomeQueue(
        ExecutionQueueEvidenceState.MISSING,
        error=queue_error,
    )
    service = RunRecoveryService(repository, queue, clock=lambda: NOW)

    with pytest.raises(RunRecoveryError) as error:
        service.requeue_run(
            "run-1",
            expected_status=RunStatus.QUEUED,
            expected_assignment=assignment,
        )

    persisted = repository.get_execution_assignment("run-1")
    assert error.value.code == expected_code
    assert "/private" not in str(error.value)
    assert persisted is not None
    assert persisted.requeue_requested_at == NOW
    assert persisted.requeue_confirmed_at is None


def test_requeue_rejects_a_backend_job_identity_mismatch_after_durable_prepare():
    repository, assignment = _queued_repository()
    queue = _MutationOutcomeQueue(
        ExecutionQueueEvidenceState.MISSING,
        job_id="different-job",
    )
    service = RunRecoveryService(repository, queue, clock=lambda: NOW)

    with pytest.raises(RunRecoveryError) as error:
        service.requeue_run(
            "run-1",
            expected_status=RunStatus.QUEUED,
            expected_assignment=assignment,
        )

    persisted = repository.get_execution_assignment("run-1")
    assert error.value.code == "RUN_RECOVERY_NOT_SAFE"
    assert persisted is not None
    assert persisted.requeue_requested_at == NOW
    assert persisted.requeue_confirmed_at is None


@pytest.mark.parametrize(
    ("repository_error", "expected_code"),
    [
        (ConcurrentRunUpdateError("race"), "RUN_RECOVERY_CONFLICT"),
        (KeyError("missing"), "RUN_RECOVERY_CONFLICT"),
        (TypeError("/private/type"), "RUN_RECOVERY_DATA_INVALID"),
        (ValueError("/private/value"), "RUN_RECOVERY_DATA_INVALID"),
    ],
)
def test_requeue_prepare_maps_repository_failures_without_queue_delivery(
    repository_error,
    expected_code,
):
    repository, assignment = _queued_repository()
    queue = _Queue(ExecutionQueueEvidenceState.MISSING)
    service = RunRecoveryService(
        _RepositoryMethodOverride(
            repository,
            "prepare_execution_requeue",
            error=repository_error,
        ),
        queue,
        clock=lambda: NOW,
    )

    with pytest.raises(RunRecoveryError) as error:
        service.requeue_run(
            "run-1",
            expected_status=RunStatus.QUEUED,
            expected_assignment=assignment,
        )

    assert error.value.code == expected_code
    assert "/private" not in str(error.value)
    assert queue.requeued == []
    assert repository.get_execution_assignment("run-1") == assignment


@pytest.mark.parametrize(
    ("repository_error", "expected_code"),
    [
        (ConcurrentRunUpdateError("race"), "RUN_RECOVERY_CONFLICT"),
        (KeyError("missing"), "RUN_RECOVERY_CONFLICT"),
        (TypeError("/private/type"), "RUN_RECOVERY_DATA_INVALID"),
        (ValueError("/private/value"), "RUN_RECOVERY_DATA_INVALID"),
    ],
)
def test_requeue_confirmation_failure_preserves_auditable_unconfirmed_intent(
    repository_error,
    expected_code,
):
    repository, assignment = _queued_repository()
    queue = _Queue(ExecutionQueueEvidenceState.MISSING)
    service = RunRecoveryService(
        _RepositoryMethodOverride(
            repository,
            "confirm_execution_requeue",
            error=repository_error,
        ),
        queue,
        clock=lambda: NOW,
    )

    with pytest.raises(RunRecoveryError) as error:
        service.requeue_run(
            "run-1",
            expected_status=RunStatus.QUEUED,
            expected_assignment=assignment,
        )

    persisted = repository.get_execution_assignment("run-1")
    assert error.value.code == expected_code
    assert "/private" not in str(error.value)
    assert queue.requeued != []
    assert persisted is not None
    assert persisted.requeue_requested_at == NOW
    assert persisted.requeue_confirmed_at is None


@pytest.mark.parametrize(
    ("repository_error", "expected_code"),
    [
        (ConcurrentRunUpdateError("race"), "RUN_RECOVERY_CONFLICT"),
        (KeyError("missing"), "RUN_RECOVERY_CONFLICT"),
        (TypeError("/private/type"), "RUN_RECOVERY_DATA_INVALID"),
        (ValueError("/private/value"), "RUN_RECOVERY_DATA_INVALID"),
    ],
)
def test_admin_failure_maps_repository_cas_errors_without_lifecycle_mutation(
    repository_error,
    expected_code,
):
    repository, assignment = _queued_repository()
    service = RunRecoveryService(
        _RepositoryMethodOverride(
            repository,
            "fail_run_by_recovery",
            error=repository_error,
        ),
        _Queue(ExecutionQueueEvidenceState.MISSING),
        clock=lambda: NOW,
    )

    with pytest.raises(RunRecoveryError) as error:
        service.fail_run(
            "run-1",
            expected_status=RunStatus.QUEUED,
            expected_assignment=assignment,
        )

    assert error.value.code == expected_code
    assert "/private" not in str(error.value)
    assert repository.get_run("run-1").status is RunStatus.QUEUED


def _queued_repository() -> tuple[InMemoryRunRepository, RunExecutionAssignment]:
    repository = InMemoryRunRepository()
    _create_execution_ready_run(repository, RunStatus.QUEUED)
    assignment = RunExecutionAssignment(
        run_id="run-1",
        job_id="job-1",
        backend="rq",
        queue_name="default",
        created_at=NOW,
        dispatched_at=NOW,
    )
    repository.ensure_execution_assignment(
        assignment,
        expected_status=RunStatus.QUEUED,
    )
    return repository, assignment


def _running_repository(
    *,
    cancellation_requested: bool = False,
) -> tuple[InMemoryRunRepository, RunExecutionAssignment]:
    repository = InMemoryRunRepository()
    _create_execution_ready_run(repository, RunStatus.RUNNING)
    assignment = RunExecutionAssignment(
        run_id="run-1",
        job_id="job-1",
        backend="rq",
        queue_name="default",
        created_at=NOW,
        dispatched_at=NOW,
        claimed_at=NOW,
        cancellation_requested_at=NOW if cancellation_requested else None,
        cancellation_reason=(
            "User requested cancellation." if cancellation_requested else None
        ),
    )
    repository.ensure_execution_assignment(
        assignment,
        expected_status=RunStatus.RUNNING,
    )
    return repository, assignment


def _create_execution_ready_run(
    repository: InMemoryRunRepository,
    status: RunStatus,
) -> None:
    repository.create_run(
        _record(RunStatus.VALIDATING),
        RunEventDraft(
            event_type="status_changed",
            message="Run validating.",
            status=RunStatus.VALIDATING,
        ),
    )
    repository.complete_preflight(
        _record(RunStatus.PLANNED),
        WorkflowBuildIdentity(
            workflow_id="fake",
            adapter_version="1.0.0",
            scheme="sha256-tree-v1",
            logical_entrypoint="workflow/main.nf",
            digest="a" * 64,
            captured_at=NOW,
        ),
        expected_status=RunStatus.VALIDATING,
        event=RunEventDraft(
            event_type="preflight_completed",
            message="Preflight completed.",
            status=RunStatus.PLANNED,
        ),
    )
    if status is RunStatus.PLANNED:
        return
    repository.update_run(
        _record(status),
        expected_status=RunStatus.PLANNED,
        event=RunEventDraft(
            event_type="status_changed",
            message="Run entered execution.",
            status=status,
        ),
    )


def _record(status: RunStatus) -> RunRecord:
    return RunRecord(
        run_id="run-1",
        workflow_id="fake",
        inputs={"private_path": "/tmp/private/sample.tsv"},
        status=status,
        created_at=NOW,
        updated_at=NOW,
        started_at=NOW if status is RunStatus.RUNNING else None,
        ended_at=None,
        current_stage="execution",
        cancellation_reason=None,
        error=None,
        tags={},
    )
