"""Tests for the RQ queue adapter and its serialized payload contract."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timezone

import fakeredis
import pytest
from rq.job import JobStatus
from rq.serializers import JSONSerializer

from encode_pipeline.workers import rq_queue
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.services.run_queue import RunQueue
from encode_pipeline.workers.rq_queue import (
    FAILURE_TTL_SECONDS,
    RESULT_TTL_SECONDS,
    RqRunQueue,
    RunQueueIdentityError,
    RunQueueJobUnavailableError,
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
    assert job.timeout == configured.job_timeout_seconds
    assert job.result_ttl == RESULT_TTL_SECONDS
    assert job.failure_ttl == FAILURE_TTL_SECONDS


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


def test_rq_run_queue_returns_successful_terminal_duplicate(tmp_path):
    connection = fakeredis.FakeRedis()
    configured = worker_settings(tmp_path)
    run_queue = RqRunQueue(configured, connection=connection)
    assignment = _assignment(configured.queue_name)
    run_queue.enqueue_execution(assignment)
    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.set_status(JobStatus.FINISHED)

    assert run_queue.enqueue_execution(assignment) == assignment.job_id


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

    with pytest.raises(RunQueueJobUnavailableError, match="scheduling state"):
        run_queue.enqueue_execution(assignment)

    assert len(run_queue._queue) == 0


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
