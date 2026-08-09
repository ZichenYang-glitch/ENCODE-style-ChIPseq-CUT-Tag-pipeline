"""Tests for the RQ queue adapter and its serialized payload contract."""

from __future__ import annotations

import os
from dataclasses import replace
from datetime import datetime, timezone

import fakeredis
import pytest
from redis.exceptions import ConnectionError as RedisConnectionError
from redis.exceptions import TimeoutError as RedisTimeoutError
from rq import Worker
from rq.exceptions import DeserializationError, DuplicateJobError, InvalidJobOperation
from rq.job import Job, JobStatus
from rq.registry import DeferredJobRegistry, ScheduledJobRegistry
from rq.serializers import JSONSerializer

from encode_pipeline.workers import rq_queue
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.run_recovery import (
    ExecutionQueueEvidenceState,
    RunExecutionQueueEvidence,
)
from encode_pipeline.services.run_queue import (
    RunQueue,
    RunQueueInspector,
    RunRecoveryQueue,
)
from encode_pipeline.workers.rq_queue import (
    FAILURE_TTL_SECONDS,
    RQ_JOB_CLEANUP_GRACE_SECONDS,
    RQ_JOB_STARTUP_ALLOWANCE_SECONDS,
    REQUEUE_REQUEST_META_KEY,
    RESULT_TTL_SECONDS,
    RqRunQueue,
    RunQueueIdentityError,
    RunQueueJobUnavailableError,
    RunQueueUnavailableError,
    RunQueueStopUnavailableError,
    STOPPED_CALLBACK_PATH,
    STOPPED_CALLBACK_TIMEOUT_SECONDS,
    create_api_redis_connection,
    create_worker_redis_connection,
    rq_job_timeout_seconds,
)

from .conftest import worker_settings


def test_rq_run_queue_enqueues_only_run_id_with_canonical_job_identity(tmp_path):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment(configured.queue_name)

    returned_job_id = run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)

    assert isinstance(run_queue, RunQueue)
    assert returned_job_id == assignment.job_id
    assert job is not None
    assert job.func_name == "encode_pipeline.workers.jobs.run_execution_job"
    assert job.args == [assignment.run_id]
    assert job.kwargs == {}
    assert job.id == assignment.job_id
    assert job.origin == configured.queue_name
    assert job.serializer is JSONSerializer
    assert job.timeout == rq_job_timeout_seconds(configured.job_timeout_seconds)
    assert job.timeout - configured.job_timeout_seconds == (
        RQ_JOB_STARTUP_ALLOWANCE_SECONDS + RQ_JOB_CLEANUP_GRACE_SECONDS
    )
    assert RQ_JOB_STARTUP_ALLOWANCE_SECONDS == 300
    assert job.timeout - configured.job_timeout_seconds == (
        rq_queue.RQ_JOB_STARTUP_ALLOWANCE_SECONDS + RQ_JOB_CLEANUP_GRACE_SECONDS
    )
    assert job.result_ttl == RESULT_TTL_SECONDS
    assert job.failure_ttl == FAILURE_TTL_SECONDS
    assert job._stopped_callback_name == STOPPED_CALLBACK_PATH
    assert job.stopped_callback_timeout == STOPPED_CALLBACK_TIMEOUT_SECONDS


def test_rq_timeout_keeps_fixed_cleanup_window_for_one_second_workflow():
    assert rq_job_timeout_seconds(1) == (
        RQ_JOB_STARTUP_ALLOWANCE_SECONDS + 1 + RQ_JOB_CLEANUP_GRACE_SECONDS
    )


def test_pre_spawn_work_cannot_consume_the_process_runner_or_cleanup_budgets():
    workflow_timeout_seconds = 20
    pre_spawn_seconds = RQ_JOB_CLEANUP_GRACE_SECONDS + 1
    process_runner_deadline = pre_spawn_seconds + workflow_timeout_seconds
    required_outer_deadline = process_runner_deadline + RQ_JOB_CLEANUP_GRACE_SECONDS

    assert rq_job_timeout_seconds(workflow_timeout_seconds) >= required_outer_deadline


def test_redis_connection_profiles_keep_api_commands_bounded_and_worker_reads_blocking(
    tmp_path,
):
    configured = replace(
        worker_settings(tmp_path),
        redis_url=(
            "redis://localhost:6379/0?socket_connect_timeout=99&socket_timeout=99"
            "&retry_on_timeout=true"
        ),
        redis_connect_timeout_seconds=1.25,
        redis_api_read_timeout_seconds=4.5,
    )

    api_connection = create_api_redis_connection(configured)
    worker_connection = create_worker_redis_connection(configured)
    try:
        api_options = api_connection.connection_pool.connection_kwargs
        worker_options = worker_connection.connection_pool.connection_kwargs

        assert api_options["socket_connect_timeout"] == 1.25
        assert api_options["socket_timeout"] == 4.5
        assert api_options["retry_on_timeout"] is False
        assert worker_options["socket_connect_timeout"] == 1.25
        assert worker_options["socket_timeout"] is None
        assert worker_options["retry_on_timeout"] is False
    finally:
        api_connection.close()
        worker_connection.close()


def test_rq_run_queue_requires_a_durable_assignment(tmp_path):
    run_queue = RqRunQueue(
        worker_settings(tmp_path),
        connection=fakeredis.FakeRedis(),
    )

    with pytest.raises(ValueError, match="RunExecutionAssignment"):
        run_queue.enqueue_execution(object())  # type: ignore[arg-type]


@pytest.mark.parametrize(
    "assignment_changes",
    [
        {"backend": "other"},
        {"queue_name": "other"},
    ],
)
def test_rq_run_queue_rejects_assignment_configuration_drift(
    tmp_path,
    assignment_changes,
):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())
    assignment = replace(
        _assignment(configured.queue_name),
        **assignment_changes,
    )

    with pytest.raises(RunQueueIdentityError, match="configured RQ queue"):
        run_queue.enqueue_execution(assignment)


def test_rq_run_queue_is_idempotent_for_the_same_durable_identity(tmp_path):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment(configured.queue_name)

    first = run_queue.enqueue_execution(assignment)
    second = run_queue.enqueue_execution(assignment)

    assert first == second == assignment.job_id
    assert len(run_queue._queue) == 1


def test_pending_requeue_duplicate_atomically_binds_the_request_marker(tmp_path):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    initial = _assignment(configured.queue_name)
    run_queue.enqueue_execution(initial)
    requested_at = datetime.now(timezone.utc)
    pending = replace(
        initial,
        dispatched_at=requested_at,
        requeue_requested_at=requested_at,
    )

    assert run_queue.enqueue_execution(pending) == pending.job_id
    job = run_queue._queue.fetch_job(pending.job_id)
    assert job is not None
    assert job.meta == {REQUEUE_REQUEST_META_KEY: _requeue_marker(pending)}
    job.set_status(JobStatus.FAILED)
    assert run_queue.inspect_execution(pending) == RunExecutionQueueEvidence(
        state=ExecutionQueueEvidenceState.FAILED,
        requeue_delivery_matches_request=True,
    )


@pytest.mark.parametrize(
    "ghost_status",
    [JobStatus.QUEUED, JobStatus.DEFERRED, JobStatus.SCHEDULED],
)
def test_rq_run_queue_rejects_exact_schedulable_job_without_queue_ownership(
    tmp_path,
    ghost_status,
):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    run_queue._queue.remove(job)
    job.set_status(ghost_status)
    job.save()

    with pytest.raises(RunQueueJobUnavailableError, match="scheduling state"):
        run_queue.enqueue_execution(assignment)

    preserved = run_queue._queue.fetch_job(assignment.job_id)
    assert preserved is not None
    assert preserved.get_status() is ghost_status
    assert len(run_queue._queue) == 0


def test_rq_run_queue_rejects_duplicate_started_job_without_live_owner(tmp_path):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.set_status(JobStatus.STARTED)
    job.worker_name = "absent-worker"
    job.save()

    with pytest.raises(RunQueueJobUnavailableError, match="scheduling state"):
        run_queue.enqueue_execution(assignment)

    assert job.get_status(refresh=True) is JobStatus.STARTED


def test_rq_run_queue_rejects_duplicate_without_truthful_stop_callback(tmp_path):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment(configured.queue_name)
    job = run_queue._queue.create_job(
        "encode_pipeline.workers.jobs.run_execution_job",
        args=(assignment.run_id,),
        kwargs={},
        job_id=assignment.job_id,
        status=JobStatus.QUEUED,
    )
    job.save()

    with pytest.raises(RunQueueIdentityError, match="durable execution identity"):
        run_queue.enqueue_execution(assignment)


def test_rq_run_queue_rejects_job_id_reuse_for_another_run(tmp_path):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)

    with pytest.raises(RunQueueIdentityError, match="durable execution identity"):
        run_queue.enqueue_execution(replace(assignment, run_id="run-other"))

    assert len(run_queue._queue) == 1


@pytest.mark.parametrize(
    "terminal_status",
    [JobStatus.FAILED, JobStatus.STOPPED, JobStatus.CANCELED],
)
def test_rq_run_queue_rejects_unsuccessful_terminal_duplicate(
    tmp_path,
    terminal_status,
):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.set_status(terminal_status)

    with pytest.raises(RunQueueJobUnavailableError, match="scheduling state"):
        run_queue.enqueue_execution(assignment)


def test_rq_run_queue_rejects_successful_terminal_duplicate(tmp_path):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.set_status(JobStatus.FINISHED)

    with pytest.raises(RunQueueJobUnavailableError, match="scheduling state"):
        run_queue.enqueue_execution(assignment)


def test_rq_run_queue_rejects_created_job_that_was_never_queued(tmp_path):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment(configured.queue_name)
    job = run_queue._queue.create_job(
        "encode_pipeline.workers.jobs.run_execution_job",
        args=(assignment.run_id,),
        kwargs={},
        job_id=assignment.job_id,
        status=JobStatus.CREATED,
    )
    job.save()
    assert len(run_queue._queue) == 0

    with pytest.raises(RunQueueIdentityError, match="durable execution identity"):
        run_queue.enqueue_execution(assignment)

    assert len(run_queue._queue) == 0


def test_rq_run_queue_maps_duplicate_job_deletion_race_to_unavailable(
    tmp_path,
    monkeypatch,
):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None

    def duplicate(*_args, **_kwargs):
        raise DuplicateJobError

    def deleted_during_refresh(*_args, **_kwargs):
        raise InvalidJobOperation("job was deleted")

    monkeypatch.setattr(run_queue._queue, "enqueue", duplicate)
    monkeypatch.setattr(run_queue._queue, "fetch_job", lambda _job_id: job)
    monkeypatch.setattr(job, "get_status", deleted_during_refresh)

    with pytest.raises(RunQueueJobUnavailableError, match="scheduling state"):
        run_queue.enqueue_execution(assignment)


@pytest.mark.parametrize("error_type", [RedisConnectionError, RedisTimeoutError])
def test_rq_run_queue_sanitizes_backend_connection_errors(
    tmp_path,
    monkeypatch,
    error_type,
):
    run_queue = RqRunQueue(
        worker_settings(tmp_path),
        connection=fakeredis.FakeRedis(),
    )
    assignment = _assignment(run_queue.queue_name)

    def unavailable(*_args, **_kwargs):
        raise error_type("redis://password@private-host:6379")

    monkeypatch.setattr(run_queue._queue, "enqueue", unavailable)

    with pytest.raises(RunQueueUnavailableError) as raised:
        run_queue.enqueue_execution(assignment)

    assert "private-host" not in str(raised.value)


def test_rq_run_queue_sends_public_stop_command_for_strict_started_identity(
    tmp_path,
    monkeypatch,
):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment_with_cancellation_intent(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    assert job.func_name == "encode_pipeline.workers.jobs.run_execution_job"
    job.set_status(JobStatus.STARTED)
    job.worker_name = "worker-1"
    job.save()
    calls = []

    def send_stop(connection, job_id, serializer=None):
        calls.append((connection, job_id, serializer))

    monkeypatch.setattr(rq_queue, "send_stop_job_command", send_stop)

    run_queue.request_stop(assignment)

    assert calls == [(connection, assignment.job_id, JSONSerializer)]


@pytest.mark.parametrize(
    "mutation",
    ["wrong_args", "wrong_origin", "missing_callback", "queued", "no_worker"],
)
def test_rq_run_queue_never_stops_a_job_with_mismatched_identity(
    tmp_path,
    monkeypatch,
    mutation,
):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment_with_cancellation_intent(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    assert job.func_name == "encode_pipeline.workers.jobs.run_execution_job"
    job.set_status(JobStatus.STARTED)
    job.worker_name = "worker-1"
    if mutation == "wrong_args":
        job.args = ["another-run"]
    elif mutation == "wrong_origin":
        job.origin = "another-queue"
    elif mutation == "missing_callback":
        job._stopped_callback_name = None
    elif mutation == "queued":
        job.set_status(JobStatus.QUEUED)
    elif mutation == "no_worker":
        job.worker_name = None
    job.save()
    calls = []
    monkeypatch.setattr(
        rq_queue,
        "send_stop_job_command",
        lambda *_args, **_kwargs: calls.append(True),
    )

    with pytest.raises(RunQueueStopUnavailableError):
        run_queue.request_stop(assignment)

    assert calls == []


@pytest.mark.parametrize("error_type", [RedisConnectionError, RedisTimeoutError])
def test_rq_run_queue_sanitizes_stop_backend_errors(
    tmp_path,
    monkeypatch,
    error_type,
):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment_with_cancellation_intent(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.set_status(JobStatus.STARTED)
    job.worker_name = "worker-1"
    job.save()

    def unavailable(*_args, **_kwargs):
        raise error_type("redis://password@private-host:6379")

    monkeypatch.setattr(rq_queue, "send_stop_job_command", unavailable)

    with pytest.raises(RunQueueStopUnavailableError) as raised:
        run_queue.request_stop(assignment)

    assert "private-host" not in str(raised.value)


@pytest.mark.parametrize(
    ("status", "expected_state"),
    [
        (JobStatus.CREATED, ExecutionQueueEvidenceState.UNKNOWN),
        (JobStatus.QUEUED, ExecutionQueueEvidenceState.QUEUED),
        (JobStatus.DEFERRED, ExecutionQueueEvidenceState.DEFERRED),
        (JobStatus.SCHEDULED, ExecutionQueueEvidenceState.SCHEDULED),
        (JobStatus.FINISHED, ExecutionQueueEvidenceState.FINISHED),
        (JobStatus.FAILED, ExecutionQueueEvidenceState.FAILED),
        (JobStatus.STOPPED, ExecutionQueueEvidenceState.STOPPED),
        (JobStatus.CANCELED, ExecutionQueueEvidenceState.CANCELED),
    ],
)
def test_rq_run_queue_projects_every_non_started_rq_status(
    tmp_path,
    status,
    expected_state,
):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.set_status(status)
    _register_scheduling_status(job, status)

    evidence = run_queue.inspect_execution(assignment)

    assert isinstance(run_queue, RunQueueInspector)
    assert evidence == RunExecutionQueueEvidence(state=expected_state)


@pytest.mark.parametrize(
    "status",
    [JobStatus.QUEUED, JobStatus.DEFERRED, JobStatus.SCHEDULED],
)
def test_rq_run_queue_requires_schedulable_status_to_be_queue_owned(
    tmp_path,
    status,
):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())
    assignment = _assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    run_queue._queue.remove(job)
    job.set_status(status)

    assert run_queue.inspect_execution(assignment) == RunExecutionQueueEvidence(
        state=ExecutionQueueEvidenceState.UNKNOWN
    )


def test_rq_run_queue_reports_missing_execution_without_backend_details(tmp_path):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())

    evidence = run_queue.inspect_execution(_assignment(configured.queue_name))

    assert evidence == RunExecutionQueueEvidence(
        state=ExecutionQueueEvidenceState.MISSING
    )
    assert vars(evidence) == {
        "state": ExecutionQueueEvidenceState.MISSING,
        "requeue_delivery_matches_request": False,
    }


def test_rq_run_queue_inspection_does_not_prune_a_stale_queue_entry(tmp_path):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment(configured.queue_name)
    connection.rpush(run_queue._queue.key, assignment.job_id)
    queue_before = connection.lrange(run_queue._queue.key, 0, -1)

    assert connection.exists(Job.key_for(assignment.job_id)) == 0
    assert run_queue.inspect_execution(assignment) == RunExecutionQueueEvidence(
        state=ExecutionQueueEvidenceState.MISSING
    )
    assert connection.lrange(run_queue._queue.key, 0, -1) == queue_before
    assert connection.exists(Job.key_for(assignment.job_id)) == 0


@pytest.mark.parametrize(
    ("attribute", "value"),
    [
        ("func_name", "private.module.wrong_job"),
        ("args", ["another-run"]),
        ("kwargs", {"secret": "value"}),
        ("origin", "another-queue"),
        ("_stopped_callback_name", None),
        ("_stopped_callback_timeout", 999),
    ],
)
def test_rq_run_queue_reports_complete_identity_drift(
    tmp_path,
    monkeypatch,
    attribute,
    value,
):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())
    assignment = _assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    setattr(job, attribute, value)
    monkeypatch.setattr(
        rq_queue,
        "_fetch_job_read_only",
        lambda _queue, _job_id: job,
    )

    evidence = run_queue.inspect_execution(assignment)

    assert evidence == RunExecutionQueueEvidence(
        state=ExecutionQueueEvidenceState.IDENTITY_DRIFT
    )
    assert "private" not in repr(evidence)


def test_rq_run_queue_reports_assignment_configuration_drift(tmp_path):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())
    assignment = replace(_assignment(configured.queue_name), queue_name="private-queue")

    evidence = run_queue.inspect_execution(assignment)

    assert evidence == RunExecutionQueueEvidence(
        state=ExecutionQueueEvidenceState.IDENTITY_DRIFT
    )
    assert "private-queue" not in repr(evidence)


@pytest.mark.parametrize(
    "error",
    [
        RedisConnectionError("redis://password@private-host:6379"),
        RedisTimeoutError("redis://password@private-host:6379"),
        DeserializationError(),
    ],
)
def test_rq_run_queue_maps_inspection_failures_to_redacted_unavailable_evidence(
    tmp_path,
    monkeypatch,
    error,
):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())

    def unavailable(_queue, _job_id):
        raise error

    monkeypatch.setattr(rq_queue, "_fetch_job_read_only", unavailable)

    evidence = run_queue.inspect_execution(_assignment(configured.queue_name))

    assert evidence == RunExecutionQueueEvidence(
        state=ExecutionQueueEvidenceState.UNAVAILABLE
    )
    assert "private-host" not in repr(evidence)


def test_rq_run_queue_fails_closed_when_job_is_deleted_during_status_refresh(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())
    assignment = _assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    monkeypatch.setattr(
        rq_queue,
        "_fetch_job_read_only",
        lambda _queue, _job_id: job,
    )

    def deleted(*_args, **_kwargs):
        raise InvalidJobOperation("private deletion race")

    monkeypatch.setattr(job, "get_status", deleted)

    assert run_queue.inspect_execution(assignment) == RunExecutionQueueEvidence(
        state=ExecutionQueueEvidenceState.UNAVAILABLE
    )


def test_rq_run_queue_proves_started_owner_only_from_exact_live_local_process(
    tmp_path,
    monkeypatch,
):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.set_status(JobStatus.STARTED)
    worker = Worker(
        [run_queue._queue],
        name="recovery-worker",
        connection=connection,
        serializer=JSONSerializer,
    )
    worker.register_birth()
    connection.hset(worker.key, "current_job", assignment.job_id)
    job.worker_name = worker.name
    job.save()
    kill_calls = []
    monkeypatch.setattr(
        rq_queue.os,
        "kill",
        lambda pid, signal: kill_calls.append((pid, signal)),
    )

    evidence = run_queue.inspect_execution(assignment)

    assert evidence == RunExecutionQueueEvidence(
        state=ExecutionQueueEvidenceState.STARTED_LIVE
    )
    assert kill_calls == [(os.getpid(), 0)]
    assert vars(evidence) == {
        "state": ExecutionQueueEvidenceState.STARTED_LIVE,
        "requeue_delivery_matches_request": False,
    }


@pytest.mark.parametrize(
    "owner_mutation",
    [
        "missing_name",
        "missing_worker",
        "wrong_job",
        "remote_host",
        "nonpositive_pid",
        "dead_pid",
    ],
)
def test_rq_run_queue_keeps_started_owner_unproven_without_complete_local_proof(
    tmp_path,
    monkeypatch,
    owner_mutation,
):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.set_status(JobStatus.STARTED)
    worker = Worker(
        [run_queue._queue],
        name="recovery-worker",
        connection=connection,
        serializer=JSONSerializer,
    )
    worker.register_birth()
    connection.hset(worker.key, "current_job", assignment.job_id)
    job.worker_name = worker.name
    if owner_mutation == "missing_name":
        job.worker_name = None
    elif owner_mutation == "missing_worker":
        job.worker_name = "absent-worker"
    elif owner_mutation == "wrong_job":
        connection.hset(worker.key, "current_job", "another-job")
    elif owner_mutation == "remote_host":
        connection.hset(worker.key, "hostname", "private-remote-host")
    elif owner_mutation == "nonpositive_pid":
        connection.hset(worker.key, "pid", "0")
    job.save()

    def check_process(pid, signal):
        assert signal == 0
        if owner_mutation == "dead_pid":
            raise ProcessLookupError(pid)

    monkeypatch.setattr(rq_queue.os, "kill", check_process)

    evidence = run_queue.inspect_execution(assignment)

    assert evidence == RunExecutionQueueEvidence(
        state=ExecutionQueueEvidenceState.STARTED_UNPROVEN
    )
    assert "recovery-worker" not in repr(evidence)
    assert "private-remote-host" not in repr(evidence)


def test_rq_run_queue_missing_worker_inspection_never_mutates_registry(
    tmp_path,
    monkeypatch,
):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.set_status(JobStatus.STARTED)
    job.worker_name = "absent-worker"
    job.save()

    def reject_registry_write(*_args, **_kwargs):
        raise AssertionError("read-only inspection mutated the worker registry")

    monkeypatch.setattr(connection, "srem", reject_registry_write)

    assert run_queue.inspect_execution(assignment) == RunExecutionQueueEvidence(
        state=ExecutionQueueEvidenceState.STARTED_UNPROVEN
    )


def test_rq_run_queue_live_owner_does_not_infer_death_from_elapsed_time(
    tmp_path,
    monkeypatch,
):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.set_status(JobStatus.STARTED)
    worker = Worker(
        [run_queue._queue],
        name="recovery-worker",
        connection=connection,
        serializer=JSONSerializer,
    )
    worker.register_birth()
    connection.hset(
        worker.key,
        mapping={
            "current_job": assignment.job_id,
            "last_heartbeat": "2000-01-01T00:00:00.000000Z",
        },
    )
    job.worker_name = worker.name
    job.save()
    monkeypatch.setattr(rq_queue.os, "kill", lambda _pid, _signal: None)

    assert run_queue.inspect_execution(assignment) == RunExecutionQueueEvidence(
        state=ExecutionQueueEvidenceState.STARTED_LIVE
    )


def test_rq_requeue_missing_job_uses_the_same_stable_identity(tmp_path):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())
    assignment = _requeue_assignment(configured.queue_name)

    job_id = run_queue.requeue_execution(assignment)
    repeated_job_id = run_queue.requeue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)

    assert isinstance(run_queue, RunRecoveryQueue)
    assert job_id == repeated_job_id == assignment.job_id
    assert job is not None
    assert job.args == [assignment.run_id]
    assert job.origin == assignment.queue_name
    assert job._stopped_callback_name == STOPPED_CALLBACK_PATH
    assert job.stopped_callback_timeout == STOPPED_CALLBACK_TIMEOUT_SECONDS
    assert job.meta == {REQUEUE_REQUEST_META_KEY: _requeue_marker(assignment)}


@pytest.mark.parametrize(
    "terminal_status",
    [JobStatus.FINISHED, JobStatus.FAILED, JobStatus.STOPPED, JobStatus.CANCELED],
)
def test_rq_requeue_replaces_only_an_exact_terminal_record(
    tmp_path,
    terminal_status,
):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())
    assignment = _requeue_assignment(configured.queue_name)
    run_queue.enqueue_execution(replace(assignment, requeue_requested_at=None))
    old_job = run_queue._queue.fetch_job(assignment.job_id)
    assert old_job is not None
    old_job.meta["old-terminal-record"] = True
    old_job.save_meta()
    old_job.set_status(terminal_status)

    assert run_queue.requeue_execution(assignment) == assignment.job_id
    replacement = run_queue._queue.fetch_job(assignment.job_id)

    assert replacement is not None
    assert replacement.get_status() is JobStatus.QUEUED
    assert replacement.meta == {REQUEUE_REQUEST_META_KEY: _requeue_marker(assignment)}
    assert len(run_queue._queue) == 1


@pytest.mark.parametrize(
    "active_status",
    [JobStatus.QUEUED, JobStatus.DEFERRED, JobStatus.SCHEDULED],
)
def test_rq_requeue_is_idempotent_for_exact_active_job_with_pending_request(
    tmp_path,
    active_status,
):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())
    assignment = _requeue_assignment(configured.queue_name)
    run_queue.enqueue_execution(replace(assignment, requeue_requested_at=None))
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.set_status(active_status)
    _register_scheduling_status(job, active_status)
    job.meta["must-survive"] = True
    job.save_meta()

    first = run_queue.requeue_execution(assignment)
    second = run_queue.requeue_execution(assignment)

    assert first == second == assignment.job_id
    preserved = run_queue._queue.fetch_job(assignment.job_id)
    assert preserved is not None
    assert preserved.meta == {
        REQUEUE_REQUEST_META_KEY: _requeue_marker(assignment),
        "must-survive": True,
    }


def test_rq_requeue_accepts_only_a_live_started_owner(tmp_path, monkeypatch):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _requeue_assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.set_status(JobStatus.STARTED)
    worker = Worker(
        [run_queue._queue],
        name="requeue-live-worker",
        connection=connection,
        serializer=JSONSerializer,
    )
    worker.register_birth()
    connection.hset(worker.key, "current_job", assignment.job_id)
    job.worker_name = worker.name
    job.save()
    monkeypatch.setattr(rq_queue.os, "kill", lambda _pid, _signal: None)

    assert run_queue.requeue_execution(assignment) == assignment.job_id


def test_rq_requeue_rejects_started_job_without_live_owner(tmp_path):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())
    assignment = _requeue_assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.set_status(JobStatus.STARTED)
    job.worker_name = "absent-worker"
    job.save()

    with pytest.raises(RunQueueJobUnavailableError):
        run_queue.requeue_execution(assignment)

    assert job.get_status(refresh=True) is JobStatus.STARTED


def test_rq_terminal_replacement_marker_is_request_bound_and_not_replaced(
    tmp_path,
):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())
    assignment = _requeue_assignment(configured.queue_name)
    run_queue.requeue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.set_status(JobStatus.FAILED)

    evidence = run_queue.inspect_execution(assignment)
    assert evidence == RunExecutionQueueEvidence(
        state=ExecutionQueueEvidenceState.FAILED,
        requeue_delivery_matches_request=True,
    )
    assert run_queue.requeue_execution(assignment) == assignment.job_id
    assert run_queue._queue.fetch_job(assignment.job_id) is not None


def test_rq_requeue_marker_mismatch_is_identity_drift(tmp_path):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())
    assignment = _requeue_assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.meta[REQUEUE_REQUEST_META_KEY] = "2000-01-01T00:00:00.000000Z"
    job.save_meta()

    assert run_queue.inspect_execution(assignment).state is (
        ExecutionQueueEvidenceState.IDENTITY_DRIFT
    )
    with pytest.raises(RunQueueIdentityError):
        run_queue.requeue_execution(assignment)


@pytest.mark.parametrize(
    "status",
    [JobStatus.QUEUED, JobStatus.DEFERRED, JobStatus.SCHEDULED],
)
def test_rq_requeue_never_confirms_unowned_schedulable_status(tmp_path, status):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())
    assignment = _requeue_assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    run_queue._queue.remove(job)
    job.set_status(status)
    job.meta["must-survive"] = True
    job.save_meta()

    with pytest.raises(RunQueueJobUnavailableError):
        run_queue.requeue_execution(assignment)

    preserved = run_queue._queue.fetch_job(assignment.job_id)
    assert preserved is not None
    assert preserved.meta == {
        REQUEUE_REQUEST_META_KEY: _requeue_marker(assignment),
        "must-survive": True,
    }


def test_rq_requeue_fails_closed_when_queue_ownership_changes_during_read(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())
    assignment = _requeue_assignment(configured.queue_name)
    run_queue.enqueue_execution(replace(assignment, requeue_requested_at=None))
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    original_get_job_ids = run_queue._queue.get_job_ids

    def remove_after_read(*args, **kwargs):
        job_ids = original_get_job_ids(*args, **kwargs)
        run_queue._queue.remove(job)
        return job_ids

    monkeypatch.setattr(run_queue._queue, "get_job_ids", remove_after_read)

    with pytest.raises(RunQueueUnavailableError):
        run_queue.requeue_execution(assignment)

    assert run_queue._queue.fetch_job(assignment.job_id) is not None


@pytest.mark.parametrize(
    "assignment",
    [
        pytest.param("undispatched", id="undispatched"),
        pytest.param("claimed", id="claimed"),
        pytest.param("canceling", id="canceling"),
        pytest.param("unrequested", id="unrequested"),
        pytest.param("confirmed", id="confirmed"),
    ],
)
def test_rq_requeue_requires_exact_pending_durable_permission(
    tmp_path,
    assignment,
):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())
    candidate = _requeue_assignment(configured.queue_name)
    now = datetime.now(timezone.utc)
    if assignment == "undispatched":
        candidate = _assignment(configured.queue_name)
    elif assignment == "claimed":
        candidate = replace(candidate, claimed_at=now)
    elif assignment == "canceling":
        candidate = replace(
            candidate,
            claimed_at=now,
            cancellation_requested_at=now,
            cancellation_reason="private cancellation reason",
        )
    elif assignment == "unrequested":
        candidate = replace(candidate, requeue_requested_at=None)
    elif assignment == "confirmed":
        candidate = replace(candidate, requeue_confirmed_at=now)

    with pytest.raises(RunQueueJobUnavailableError) as raised:
        run_queue.requeue_execution(candidate)

    assert "private" not in str(raised.value)
    assert run_queue._queue.fetch_job(candidate.job_id) is None


def test_rq_requeue_never_deletes_identity_drift_even_when_terminal(tmp_path):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())
    assignment = _requeue_assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.set_status(JobStatus.FAILED)
    job._stopped_callback_name = None
    job.save()

    with pytest.raises(RunQueueIdentityError):
        run_queue.requeue_execution(assignment)

    assert run_queue._queue.fetch_job(assignment.job_id) is not None


def test_rq_requeue_never_deletes_unknown_exact_job(tmp_path):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())
    assignment = _requeue_assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.set_status(JobStatus.CREATED)

    with pytest.raises(RunQueueJobUnavailableError):
        run_queue.requeue_execution(assignment)

    assert run_queue._queue.fetch_job(assignment.job_id) is not None


@pytest.mark.parametrize("error_type", [RedisConnectionError, RedisTimeoutError])
def test_rq_requeue_redacts_backend_errors(
    tmp_path,
    monkeypatch,
    error_type,
):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())
    assignment = _requeue_assignment(configured.queue_name)

    def unavailable(_job_id):
        raise error_type("redis://password@private-host:6379")

    monkeypatch.setattr(run_queue._queue, "fetch_job", unavailable)

    with pytest.raises(RunQueueUnavailableError) as raised:
        run_queue.requeue_execution(assignment)

    assert "private-host" not in str(raised.value)


def test_rq_requeue_terminal_delete_race_never_deletes_newly_active_job(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=fakeredis.FakeRedis())
    assignment = _requeue_assignment(configured.queue_name)
    run_queue.enqueue_execution(replace(assignment, requeue_requested_at=None))
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.set_status(JobStatus.FAILED)
    monkeypatch.setattr(run_queue._queue, "fetch_job", lambda _job_id: job)
    original_delete = job.delete

    def become_active_before_delete(*, pipeline=None, **kwargs):
        job.connection.hset(job.key, "status", JobStatus.STARTED.value)
        return original_delete(pipeline=pipeline, **kwargs)

    monkeypatch.setattr(job, "delete", become_active_before_delete)

    with pytest.raises(RunQueueUnavailableError):
        run_queue.requeue_execution(assignment)

    assert job.connection.exists(job.key)
    assert job.get_status() is JobStatus.STARTED


def test_rq_run_queue_closes_only_owned_connections(tmp_path, monkeypatch):
    class FakeConnection:
        def __init__(self):
            self.closed = False

        def close(self):
            self.closed = True

    class FakeQueue:
        def __init__(self, connection):
            self.connection = connection

    owned = FakeConnection()

    def create_fake_queue(_settings, *, connection=None):
        return FakeQueue(owned if connection is None else connection)

    monkeypatch.setattr(rq_queue, "create_rq_queue", create_fake_queue)
    configured = worker_settings(tmp_path)
    injected = FakeConnection()

    injected_queue = RqRunQueue(configured, connection=injected)
    injected_queue.close()
    assert injected.closed is False

    owned_queue = RqRunQueue(configured)
    owned_queue.close()
    assert owned.closed is True


def _assignment(queue_name: str) -> RunExecutionAssignment:
    return RunExecutionAssignment(
        run_id="run-123",
        job_id="job-456",
        backend="rq",
        queue_name=queue_name,
        created_at=datetime.now(timezone.utc),
    )


def _assignment_with_cancellation_intent(queue_name: str) -> RunExecutionAssignment:
    now = datetime.now(timezone.utc)
    return RunExecutionAssignment(
        run_id="run-123",
        job_id="job-456",
        backend="rq",
        queue_name=queue_name,
        created_at=now,
        dispatched_at=now,
        claimed_at=now,
        cancellation_requested_at=now,
        cancellation_reason="User requested cancellation.",
    )


def _requeue_assignment(queue_name: str) -> RunExecutionAssignment:
    now = datetime.now(timezone.utc)
    return RunExecutionAssignment(
        run_id="run-123",
        job_id="job-456",
        backend="rq",
        queue_name=queue_name,
        created_at=now,
        dispatched_at=now,
        requeue_requested_at=now,
    )


def _requeue_marker(assignment: RunExecutionAssignment) -> str:
    assert assignment.requeue_requested_at is not None
    return (
        assignment.requeue_requested_at.astimezone(timezone.utc)
        .isoformat(timespec="microseconds")
        .replace("+00:00", "Z")
    )


def _register_scheduling_status(job, status: JobStatus) -> None:
    registry_type = {
        JobStatus.DEFERRED: DeferredJobRegistry,
        JobStatus.SCHEDULED: ScheduledJobRegistry,
    }.get(status)
    if registry_type is None:
        return
    registry_type(
        name=job.origin,
        connection=job.connection,
        job_class=job.__class__,
        serializer=job.serializer,
    ).add(job)
