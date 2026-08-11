"""Real-Redis coverage for administrator run-recovery handshakes."""

from __future__ import annotations

from dataclasses import dataclass, replace
from datetime import datetime, timezone
import os
from pathlib import Path
from typing import Iterator
from uuid import uuid4

import pytest
from redis import Redis
from rq.job import Job, JobStatus

from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.run_recovery import ExecutionQueueEvidenceState
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.services.run_queue import (
    RunQueueIdentityError,
    RunQueueJobUnavailableError,
)
from encode_pipeline.services.run_recovery import (
    RUN_RECOVERY_REQUEUE_CONFIRMED_EVENT,
    RUN_RECOVERY_REQUEUE_REQUESTED_EVENT,
    RunRecoveryError,
    RunRecoveryService,
)
from encode_pipeline.services.run_repositories import (
    InMemoryRunRepository,
    RunEventDraft,
)
from encode_pipeline.workers.rq_queue import REQUEUE_REQUEST_META_KEY, RqRunQueue
from encode_pipeline.workers.settings import WorkerSettings


pytestmark = pytest.mark.platform_real_execution


TEST_REDIS_URL_ENV = "ENCODE_PIPELINE_TEST_REDIS_URL"
RECOVERY_AT = datetime(2026, 8, 9, 12, 0, tzinfo=timezone.utc)
RECOVERY_MARKER = "2026-08-09T12:00:00.000000Z"


@dataclass
class _RedisHarness:
    run_queue: RqRunQueue
    namespace: str
    job_ids: set[str]


@pytest.fixture
def real_redis_harness(tmp_path: Path) -> Iterator[_RedisHarness]:
    redis_url = os.getenv(TEST_REDIS_URL_ENV)
    if redis_url is None:
        pytest.skip(f"{TEST_REDIS_URL_ENV} is not configured")

    namespace = uuid4().hex
    queue_name = f"run-recovery-{namespace}"
    settings = WorkerSettings(
        database_url=f"sqlite:///{tmp_path / 'unused.db'}",
        redis_url=redis_url,
        queue_name=queue_name,
        workspace_root=tmp_path / "workspaces",
    )
    connection = Redis.from_url(redis_url)
    run_queue = RqRunQueue(settings, connection=connection)
    job_ids: set[str] = set()
    connected = False
    try:
        assert connection.ping() is True
        connected = True
        yield _RedisHarness(
            run_queue=run_queue,
            namespace=namespace,
            job_ids=job_ids,
        )
    finally:
        if connected:
            run_queue._queue.delete(delete_jobs=True)
            for job_id in job_ids:
                job = run_queue._queue.fetch_job(job_id)
                if job is not None:
                    job.delete()
        connection.close()


def test_real_redis_diagnosis_does_not_prune_a_stale_queue_entry(
    real_redis_harness: _RedisHarness,
) -> None:
    repository, assignment = _repository_with_assignment(
        real_redis_harness,
        label="read-only-diagnosis",
    )
    run_queue = real_redis_harness.run_queue
    connection = run_queue._queue.connection
    connection.rpush(run_queue._queue.key, assignment.job_id)
    queue_before = connection.lrange(run_queue._queue.key, 0, -1)

    assert connection.exists(Job.key_for(assignment.job_id)) == 0
    diagnostic = RunRecoveryService(repository, run_queue).diagnose(assignment.run_id)

    assert diagnostic.queue_evidence.state is ExecutionQueueEvidenceState.MISSING
    assert connection.lrange(run_queue._queue.key, 0, -1) == queue_before
    assert connection.exists(Job.key_for(assignment.job_id)) == 0


@pytest.mark.parametrize("terminal_job", [False, True], ids=["missing", "terminal"])
def test_real_redis_recovery_requeues_and_confirms_the_same_durable_identity(
    real_redis_harness: _RedisHarness,
    terminal_job: bool,
) -> None:
    repository, assignment = _repository_with_assignment(
        real_redis_harness,
        label="candidate",
    )
    run_queue = real_redis_harness.run_queue
    if terminal_job:
        assert run_queue.enqueue_execution(assignment) == assignment.job_id
        old_job = run_queue._queue.fetch_job(assignment.job_id)
        assert old_job is not None
        old_job.meta["terminal-record"] = "must-be-replaced"
        old_job.save_meta()
        old_job.set_status(JobStatus.FAILED)

    expected_initial_state = (
        ExecutionQueueEvidenceState.FAILED
        if terminal_job
        else ExecutionQueueEvidenceState.MISSING
    )
    assert run_queue.inspect_execution(assignment).state is expected_initial_state

    result = RunRecoveryService(
        repository,
        run_queue,
        clock=lambda: RECOVERY_AT,
    ).requeue_run(
        assignment.run_id,
        expected_status=RunStatus.QUEUED,
        expected_assignment=assignment,
    )

    persisted = repository.get_execution_assignment(assignment.run_id)
    replacement = run_queue._queue.fetch_job(assignment.job_id)
    assert persisted == result.assignment
    assert persisted is not None
    assert persisted.job_id == assignment.job_id
    assert persisted.created_at == assignment.created_at
    assert persisted.dispatched_at == assignment.dispatched_at
    assert persisted.claimed_at is None
    assert persisted.requeue_requested_at == RECOVERY_AT
    assert persisted.requeue_confirmed_at == RECOVERY_AT
    assert result.changed is True
    assert replacement is not None
    assert replacement.id == assignment.job_id
    assert replacement.args == [assignment.run_id]
    assert replacement.get_status(refresh=True) is JobStatus.QUEUED
    assert replacement.meta == {REQUEUE_REQUEST_META_KEY: RECOVERY_MARKER}
    assert run_queue._queue.get_job_ids() == [assignment.job_id]
    assert run_queue.inspect_execution(persisted).state is (
        ExecutionQueueEvidenceState.QUEUED
    )
    assert [event.event_type for event in repository.list_events(assignment.run_id)][
        -2:
    ] == [
        RUN_RECOVERY_REQUEUE_REQUESTED_EVENT,
        RUN_RECOVERY_REQUEUE_CONFIRMED_EVENT,
    ]


def test_real_redis_recovery_retry_of_exact_active_job_is_idempotent(
    real_redis_harness: _RedisHarness,
) -> None:
    repository, assignment = _repository_with_assignment(
        real_redis_harness,
        label="retry",
    )
    preparation = repository.prepare_execution_requeue(
        assignment.run_id,
        expected_status=RunStatus.QUEUED,
        expected_assignment=assignment,
        requested_at=RECOVERY_AT,
        event=_requeue_event(RUN_RECOVERY_REQUEUE_REQUESTED_EVENT),
    )
    prepared = preparation.assignment
    run_queue = real_redis_harness.run_queue

    assert run_queue.requeue_execution(prepared) == assignment.job_id
    active_job = run_queue._queue.fetch_job(assignment.job_id)
    assert active_job is not None
    active_job.meta["survives-retry"] = True
    active_job.save_meta()
    assert run_queue.requeue_execution(prepared) == assignment.job_id

    result = RunRecoveryService(
        repository,
        run_queue,
        clock=lambda: RECOVERY_AT,
    ).requeue_run(
        assignment.run_id,
        expected_status=RunStatus.QUEUED,
        expected_assignment=prepared,
    )

    persisted = repository.get_execution_assignment(assignment.run_id)
    preserved = run_queue._queue.fetch_job(assignment.job_id)
    assert result.assignment == persisted
    assert persisted is not None
    assert persisted.requeue_requested_at == RECOVERY_AT
    assert persisted.requeue_confirmed_at == RECOVERY_AT
    assert preserved is not None
    assert preserved.meta == {
        REQUEUE_REQUEST_META_KEY: RECOVERY_MARKER,
        "survives-retry": True,
    }
    assert preserved.get_status(refresh=True) is JobStatus.QUEUED
    assert run_queue._queue.get_job_ids() == [assignment.job_id]
    event_types = [
        event.event_type for event in repository.list_events(assignment.run_id)
    ]
    assert event_types.count(RUN_RECOVERY_REQUEUE_REQUESTED_EVENT) == 1
    assert event_types.count(RUN_RECOVERY_REQUEUE_CONFIRMED_EVENT) == 1


def test_real_redis_recovery_refuses_identity_drift_and_claimed_execution(
    real_redis_harness: _RedisHarness,
) -> None:
    run_queue = real_redis_harness.run_queue
    repository, assignment = _repository_with_assignment(
        real_redis_harness,
        label="identity-drift",
    )
    foreign_assignment = replace(
        assignment,
        run_id=f"foreign-run-{real_redis_harness.namespace}",
    )
    assert run_queue.enqueue_execution(foreign_assignment) == assignment.job_id
    foreign_job = run_queue._queue.fetch_job(assignment.job_id)
    assert foreign_job is not None
    foreign_job.meta["foreign-record"] = "must-survive"
    foreign_job.save_meta()
    assert run_queue.inspect_execution(assignment).state is (
        ExecutionQueueEvidenceState.IDENTITY_DRIFT
    )

    with pytest.raises(RunRecoveryError) as drift_error:
        RunRecoveryService(
            repository,
            run_queue,
            clock=lambda: RECOVERY_AT,
        ).requeue_run(
            assignment.run_id,
            expected_status=RunStatus.QUEUED,
            expected_assignment=assignment,
        )
    assert drift_error.value.code == "RUN_RECOVERY_NOT_SAFE"
    with pytest.raises(RunQueueIdentityError):
        run_queue.requeue_execution(
            replace(assignment, requeue_requested_at=RECOVERY_AT)
        )

    preserved_foreign_job = run_queue._queue.fetch_job(assignment.job_id)
    assert repository.get_execution_assignment(assignment.run_id) == assignment
    assert preserved_foreign_job is not None
    assert preserved_foreign_job.args == [foreign_assignment.run_id]
    assert preserved_foreign_job.meta == {"foreign-record": "must-survive"}
    assert run_queue._queue.get_job_ids() == [assignment.job_id]

    claimed_repository, claimed = _repository_with_assignment(
        real_redis_harness,
        label="claimed",
        status=RunStatus.RUNNING,
        claimed=True,
    )
    assert run_queue.enqueue_execution(claimed) == claimed.job_id
    claimed_job = run_queue._queue.fetch_job(claimed.job_id)
    assert claimed_job is not None
    claimed_job.meta["claimed-record"] = "must-survive"
    claimed_job.save_meta()

    with pytest.raises(RunRecoveryError) as claimed_error:
        RunRecoveryService(
            claimed_repository,
            run_queue,
            clock=lambda: RECOVERY_AT,
        ).requeue_run(
            claimed.run_id,
            expected_status=RunStatus.RUNNING,
            expected_assignment=claimed,
        )
    assert claimed_error.value.code == "RUN_RECOVERY_NOT_SAFE"
    with pytest.raises(RunQueueJobUnavailableError):
        run_queue.requeue_execution(replace(claimed, requeue_requested_at=RECOVERY_AT))

    preserved_claimed_job = run_queue._queue.fetch_job(claimed.job_id)
    assert claimed_repository.get_execution_assignment(claimed.run_id) == claimed
    assert preserved_claimed_job is not None
    assert preserved_claimed_job.args == [claimed.run_id]
    assert preserved_claimed_job.meta == {"claimed-record": "must-survive"}
    assert set(run_queue._queue.get_job_ids()) == {
        assignment.job_id,
        claimed.job_id,
    }


def _repository_with_assignment(
    harness: _RedisHarness,
    *,
    label: str,
    status: RunStatus = RunStatus.QUEUED,
    claimed: bool = False,
) -> tuple[InMemoryRunRepository, RunExecutionAssignment]:
    run_id = f"run-recovery-{label}-{harness.namespace}"
    job_id = f"run-recovery-job-{label}-{harness.namespace}"
    harness.job_ids.add(job_id)
    repository = InMemoryRunRepository()
    validating = RunRecord(
        run_id=run_id,
        workflow_id="recovery-test-workflow",
        inputs={},
        status=RunStatus.VALIDATING,
        created_at=RECOVERY_AT,
        updated_at=RECOVERY_AT,
        started_at=None,
        ended_at=None,
        current_stage="preflight",
        cancellation_reason=None,
        error=None,
        tags={},
    )
    repository.create_run(
        validating,
        RunEventDraft(
            event_type="status_changed",
            message="Run created for real Redis recovery testing.",
            status=RunStatus.VALIDATING,
        ),
    )
    planned = replace(validating, status=RunStatus.PLANNED)
    repository.complete_preflight(
        planned,
        WorkflowBuildIdentity(
            workflow_id=planned.workflow_id,
            adapter_version="recovery-test-v1",
            scheme="sha256",
            logical_entrypoint="workflow/recovery-test",
            digest="0" * 64,
            captured_at=RECOVERY_AT,
        ),
        expected_status=RunStatus.VALIDATING,
        event=RunEventDraft(
            event_type="preflight_completed",
            message="Recovery test build identity bound.",
            status=RunStatus.PLANNED,
        ),
    )
    assignment = RunExecutionAssignment(
        run_id=run_id,
        job_id=job_id,
        backend="rq",
        queue_name=harness.run_queue.queue_name,
        created_at=RECOVERY_AT,
        dispatched_at=RECOVERY_AT,
    )
    repository.ensure_execution_assignment(
        assignment,
        expected_status=RunStatus.PLANNED,
    )
    queued = replace(
        planned,
        status=RunStatus.QUEUED,
        current_stage="execution",
    )
    assert repository.queue_dispatched_run(
        queued,
        expected_status=RunStatus.PLANNED,
        job_id=job_id,
        backend="rq",
        queue_name=harness.run_queue.queue_name,
        event=RunEventDraft(
            event_type="run_queued",
            message="Recovery test run queued.",
            status=RunStatus.QUEUED,
        ),
    )
    if status is RunStatus.RUNNING:
        claim = repository.claim_execution_assignment(
            run_id,
            job_id=job_id,
            backend="rq",
            queue_name=harness.run_queue.queue_name,
            claimed_at=RECOVERY_AT,
            allowed_statuses=frozenset({RunStatus.QUEUED}),
            event=RunEventDraft(
                event_type="execution_claimed",
                message="Recovery test execution claimed.",
                status=RunStatus.QUEUED,
            ),
        )
        assert claim.acquired is True
        assert repository.update_run(
            replace(
                queued,
                status=RunStatus.RUNNING,
                started_at=RECOVERY_AT,
            ),
            expected_status=RunStatus.QUEUED,
            event=RunEventDraft(
                event_type="status_changed",
                message="Recovery test run started.",
                status=RunStatus.RUNNING,
            ),
        )
        assignment = claim.assignment
    elif status is not RunStatus.QUEUED:
        raise ValueError("recovery test helper status is unsupported")
    if claimed is not (assignment.claimed_at is not None):
        raise ValueError("recovery test helper claim state is inconsistent")
    assert repository.get_execution_assignment(run_id) == assignment
    assert repository.get_workflow_build_identity(run_id) is not None
    assert repository.get_run(run_id).status is status
    return repository, assignment


def _requeue_event(event_type: str) -> RunEventDraft:
    return RunEventDraft(
        event_type=event_type,
        message="Run recovery handshake advanced.",
        status=RunStatus.QUEUED,
    )
