"""Integration tests for the PR123 worker ownership handshake."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timezone

import fakeredis
from rq import SimpleWorker
from rq.serializers import JSONSerializer

from encode_pipeline.persistence.runtime import open_run_persistence
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.defaults import (
    create_default_run_service,
    create_default_workflow_registry,
)
from encode_pipeline.workers.rq_queue import RqRunQueue
from encode_pipeline.workers.settings import (
    QUEUE_NAME_ENV,
    REDIS_URL_ENV,
    WORKSPACE_ROOT_ENV,
)

from .conftest import create_planned_run, worker_settings


def _configure_worker_environment(monkeypatch, configured):
    monkeypatch.setenv("ENCODE_PIPELINE_DATABASE_URL", configured.database_url)
    monkeypatch.setenv(REDIS_URL_ENV, configured.redis_url)
    monkeypatch.setenv(QUEUE_NAME_ENV, configured.queue_name)
    monkeypatch.setenv(WORKSPACE_ROOT_ENV, str(configured.workspace_root))


def _run_burst(connection, run_queue):
    worker = SimpleWorker(
        [run_queue._queue],
        connection=connection,
        serializer=JSONSerializer,
    )
    return worker.work(burst=True, logging_level="WARNING")


def _read_run(configured, run_id):
    persistence = open_run_persistence(configured.database_url)
    service = create_default_run_service(
        registry=create_default_workflow_registry(),
        repository=persistence.repository,
    )
    try:
        return (
            service.get_run(run_id),
            service.list_events(run_id, limit=100),
            service.get_execution_assignment(run_id),
        )
    finally:
        persistence.close()


def test_rq_worker_rebuilds_dependencies_and_persists_handshake_event(
    tmp_path, monkeypatch
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "handshake-run",
        assign_queue=configured.queue_name,
    )
    _configure_worker_environment(monkeypatch, configured)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    assert assignment is not None
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    record, events, persisted_assignment = _read_run(configured, "handshake-run")
    handshake = [
        event for event in events if event.event_type == "worker_dependencies_rebuilt"
    ]
    assert record.status is RunStatus.PLANNED
    assert persisted_assignment is not None
    assert persisted_assignment.dispatched_at is not None
    assert persisted_assignment.claimed_at is not None
    assert len(handshake) == 1
    assert handshake[0].context == {
        "backend": "rq",
        "job_id": assignment.job_id,
        "queue_name": configured.queue_name,
        "workflow_id": record.workflow_id,
    }

    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.delete()
    run_queue.enqueue_execution(assignment)
    assert _run_burst(connection, run_queue) is True

    _, events_after_retry, assignment_after_retry = _read_run(
        configured,
        "handshake-run",
    )
    assert assignment_after_retry == persisted_assignment
    assert [event.event_type for event in events_after_retry].count(
        "worker_dependencies_rebuilt"
    ) == 1


def test_rq_worker_rejects_stale_job_identity_without_writing_handshake(
    tmp_path, monkeypatch
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "stale-run",
        assign_queue=configured.queue_name,
    )
    _configure_worker_environment(monkeypatch, configured)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    assert assignment is not None
    stale_assignment = replace(assignment, job_id="stale-job-id")
    run_queue.enqueue_execution(stale_assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job("stale-job-id")
    _, events, _ = _read_run(configured, "stale-run")
    assert job is not None
    assert job.is_failed
    assert all(event.event_type != "worker_dependencies_rebuilt" for event in events)


def test_rq_worker_rejects_missing_durable_assignment(tmp_path, monkeypatch):
    configured = worker_settings(tmp_path)
    create_planned_run(configured, "unassigned-run")
    _configure_worker_environment(monkeypatch, configured)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    orphan_assignment = RunExecutionAssignment(
        run_id="unassigned-run",
        job_id="orphan-job-id",
        backend="rq",
        queue_name=configured.queue_name,
        created_at=datetime.now(timezone.utc),
    )
    run_queue.enqueue_execution(orphan_assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job("orphan-job-id")
    _, events, _ = _read_run(configured, "unassigned-run")
    assert job is not None
    assert job.is_failed
    assert all(event.event_type != "worker_dependencies_rebuilt" for event in events)


def test_rq_worker_rejects_job_after_run_is_cancelled(tmp_path, monkeypatch):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "cancelled-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    persistence = open_run_persistence(configured.database_url)
    try:
        run_service = create_default_run_service(
            registry=create_default_workflow_registry(),
            repository=persistence.repository,
        )
        run_service.cancel_run("cancelled-run", reason="Cancelled before worker claim.")
    finally:
        persistence.close()

    _configure_worker_environment(monkeypatch, configured)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, events, _ = _read_run(configured, "cancelled-run")
    assert job is not None
    assert job.is_failed
    assert record.status is RunStatus.CANCELLED
    assert all(event.event_type != "worker_dependencies_rebuilt" for event in events)
