"""Tests for the RQ queue adapter and its serialized payload contract."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timezone

import fakeredis
import pytest
from redis.exceptions import ConnectionError as RedisConnectionError
from redis.exceptions import TimeoutError as RedisTimeoutError
from rq.exceptions import DuplicateJobError, InvalidJobOperation
from rq.job import JobStatus
from rq.serializers import JSONSerializer

from encode_pipeline.workers import rq_queue
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.services.run_queue import RunQueue
from encode_pipeline.workers.rq_queue import (
    FAILURE_TTL_SECONDS,
    RQ_JOB_CLEANUP_GRACE_SECONDS,
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
        RQ_JOB_CLEANUP_GRACE_SECONDS
    )
    assert job.result_ttl == RESULT_TTL_SECONDS
    assert job.failure_ttl == FAILURE_TTL_SECONDS
    assert job._stopped_callback_name == STOPPED_CALLBACK_PATH
    assert job.stopped_callback_timeout == STOPPED_CALLBACK_TIMEOUT_SECONDS


def test_rq_timeout_keeps_fixed_cleanup_window_for_one_second_workflow():
    assert rq_job_timeout_seconds(1) == 1 + RQ_JOB_CLEANUP_GRACE_SECONDS


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
