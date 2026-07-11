"""Real Redis and independent-process integration for the worker boundary."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
from uuid import uuid4

import pytest
from redis import Redis

from encode_pipeline.persistence.runtime import open_run_persistence
from encode_pipeline.services.defaults import (
    create_default_run_service,
    create_default_workflow_registry,
)
from encode_pipeline.workers.rq_queue import RqRunQueue
from encode_pipeline.workers.settings import (
    QUEUE_NAME_ENV,
    REDIS_URL_ENV,
    WORKSPACE_ROOT_ENV,
    WorkerSettings,
)

from .conftest import create_planned_run


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
TEST_REDIS_URL_ENV = "ENCODE_PIPELINE_TEST_REDIS_URL"


def test_real_redis_worker_process_rebuilds_from_sqlite(tmp_path):
    """A fresh worker process receives only run_id and writes through SQLite."""
    redis_url = os.getenv(TEST_REDIS_URL_ENV)
    if redis_url is None:
        pytest.skip(f"{TEST_REDIS_URL_ENV} is not configured")

    run_id = f"redis-process-{uuid4().hex}"
    queue_name = f"encode-pipeline-test-{uuid4().hex}"
    configured = WorkerSettings(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        redis_url=redis_url,
        queue_name=queue_name,
        workspace_root=tmp_path / "workspaces",
    )
    assignment = create_planned_run(
        configured,
        run_id,
        assign_queue=queue_name,
    )
    assert assignment is not None

    connection = Redis.from_url(redis_url)
    run_queue = RqRunQueue(configured, connection=connection)
    connected = False
    try:
        assert connection.ping() is True
        connected = True
        assert run_queue.enqueue_execution(assignment) == assignment.job_id
        assert run_queue.enqueue_execution(assignment) == assignment.job_id
        assert len(run_queue._queue) == 1

        environment = dict(os.environ)
        environment.update(
            {
                "ENCODE_PIPELINE_DATABASE_URL": configured.database_url,
                REDIS_URL_ENV: configured.redis_url,
                QUEUE_NAME_ENV: configured.queue_name,
                WORKSPACE_ROOT_ENV: str(configured.workspace_root),
                "PYTHONDONTWRITEBYTECODE": "1",
                "PYTHONPATH": str(REPOSITORY_ROOT / "src"),
            }
        )
        completed = subprocess.run(
            [sys.executable, "-m", "encode_pipeline.workers.cli", "--burst"],
            cwd=REPOSITORY_ROOT,
            env=environment,
            capture_output=True,
            text=True,
            timeout=30,
        )

        assert completed.returncode == 0, completed.stderr
        job = run_queue._queue.fetch_job(assignment.job_id)
        assert job is not None
        assert job.is_finished
        assert job.args == [run_id]

        persistence = open_run_persistence(configured.database_url)
        try:
            run_service = create_default_run_service(
                registry=create_default_workflow_registry(),
                repository=persistence.repository,
            )
            events = run_service.list_events(run_id, limit=100)
            persisted_assignment = run_service.get_execution_assignment(run_id)
        finally:
            persistence.close()

        assert [event.event_type for event in events].count(
            "worker_dependencies_rebuilt"
        ) == 1
        assert persisted_assignment is not None
        assert persisted_assignment.dispatched_at is not None
        assert persisted_assignment.claimed_at is not None
    finally:
        if connected:
            job = run_queue._queue.fetch_job(assignment.job_id)
            if job is not None:
                job.delete()
            run_queue._queue.delete()
        connection.close()
