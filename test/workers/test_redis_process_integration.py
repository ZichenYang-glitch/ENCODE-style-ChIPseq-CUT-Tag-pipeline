"""Real Redis and independent-process integration for the worker boundary."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys
import time
from uuid import uuid4

import pytest
from redis import Redis
from rq import Queue
from rq.serializers import JSONSerializer

from encode_pipeline.persistence.runtime import open_run_persistence
from encode_pipeline.services.defaults import (
    create_default_run_service,
    create_default_workflow_registry,
)
from encode_pipeline.workers.rq_queue import RqRunQueue
from encode_pipeline.workers.settings import (
    JOB_TIMEOUT_SECONDS_ENV,
    QUEUE_NAME_ENV,
    REDIS_URL_ENV,
    WORKSPACE_ROOT_ENV,
    WorkerSettings,
)

from .conftest import create_planned_run
from .process_helpers import run_burst_worker, terminate_rq_worker
from .signal_timeout_helpers import CAUGHT_MARKER_ENV, ENTERED_MARKER_ENV


pytestmark = pytest.mark.platform_real_execution


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
TEST_REDIS_URL_ENV = "ENCODE_PIPELINE_TEST_REDIS_URL"


def _create_fake_snakemake(tmp_path: Path) -> Path:
    """Return a deterministic executable for the worker-boundary test."""
    executable_dir = tmp_path / "test-bin"
    executable_dir.mkdir()
    snakemake = executable_dir / "snakemake"
    snakemake.write_text(
        "#!/bin/sh\n"
        'while [ "$#" -gt 0 ]; do\n'
        '  if [ "$1" = "--directory" ]; then shift; cd "$1"; break; fi\n'
        "  shift\n"
        "done\n"
        "mkdir -p results/multiqc\n"
        "printf 'output_type\\tpath\\n' > results/multiqc/result_manifest.tsv\n"
        "printf 'worker stdout\\n'\n"
        "printf 'worker stderr\\n' >&2\n",
        encoding="utf-8",
    )
    snakemake.chmod(0o755)
    return snakemake


def test_real_redis_worker_process_rebuilds_from_sqlite(tmp_path):
    """A fresh worker receives only run_id and rebuilds from SQLite."""
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

        fake_snakemake = _create_fake_snakemake(tmp_path)
        environment = dict(os.environ)
        environment.update(
            {
                "ENCODE_PIPELINE_DATABASE_URL": configured.database_url,
                REDIS_URL_ENV: configured.redis_url,
                QUEUE_NAME_ENV: configured.queue_name,
                WORKSPACE_ROOT_ENV: str(configured.workspace_root),
                "PYTHONDONTWRITEBYTECODE": "1",
                "PYTHONPATH": str(REPOSITORY_ROOT / "src"),
                "PATH": os.pathsep.join(
                    [str(fake_snakemake.parent), environment.get("PATH", os.defpath)]
                ),
            }
        )
        completed = run_burst_worker(
            environment,
            cwd=REPOSITORY_ROOT,
            timeout_seconds=30,
        )

        assert completed.returncode == 0, completed.stderr
        job = run_queue._queue.fetch_job(assignment.job_id)
        assert job is not None
        failure = job.latest_result()
        assert job.is_finished, (
            job.get_status(refresh=True),
            failure.exc_string if failure is not None else None,
            completed.stdout,
            completed.stderr,
        )
        assert job.args == [run_id]

        persistence = open_run_persistence(configured.database_url)
        try:
            run_service = create_default_run_service(
                registry=create_default_workflow_registry(),
                repository=persistence.repository,
            )
            record = run_service.get_run(run_id)
            artifacts = run_service.list_artifacts(run_id)
            events = run_service.list_events(run_id, limit=100)
            persisted_assignment = run_service.get_execution_assignment(run_id)
            stdout = run_service.list_logs(run_id, "stdout", limit=100)
            stderr = run_service.list_logs(run_id, "stderr", limit=100)
        finally:
            persistence.close()

        assert record.status.value == "succeeded"
        assert len(artifacts) == 1
        assert artifacts[0].metadata["output_type"] == "result_manifest"
        assert [event.event_type for event in events].count(
            "worker_dependencies_rebuilt"
        ) == 1
        assert persisted_assignment is not None
        assert persisted_assignment.dispatched_at is not None
        assert persisted_assignment.claimed_at is not None
        assert stdout[0].lines == ("worker stdout",)
        assert stderr[0].lines == ("worker stderr",)
    finally:
        if connected:
            job = run_queue._queue.fetch_job(assignment.job_id)
            if job is not None:
                job.delete()
            run_queue._queue.delete()
        connection.close()


def test_real_rq_sigalrm_persists_timeout_and_reaps_worker(tmp_path):
    """A real DurableWorker SIGALRM cannot be swallowed by application code."""
    redis_url = os.getenv(TEST_REDIS_URL_ENV)
    if redis_url is None:
        pytest.skip(f"{TEST_REDIS_URL_ENV} is not configured")
    if os.name != "posix":
        pytest.skip("RQ SIGALRM integration requires a POSIX worker")

    run_id = f"redis-timeout-{uuid4().hex}"
    queue_name = f"encode-pipeline-timeout-{uuid4().hex}"
    configured = WorkerSettings(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        redis_url=redis_url,
        queue_name=queue_name,
        workspace_root=tmp_path / "workspaces",
        job_timeout_seconds=30,
    )
    assignment = create_planned_run(
        configured,
        run_id,
        assign_queue=queue_name,
    )
    assert assignment is not None

    entered_marker = tmp_path / "timeout-entered"
    caught_marker = tmp_path / "timeout-caught"
    connection = Redis.from_url(redis_url)
    queue = Queue(queue_name, connection=connection, serializer=JSONSerializer)
    worker_process: subprocess.Popen[str] | None = None
    connected = False
    try:
        assert connection.ping() is True
        connected = True
        queue.enqueue(
            "workers.signal_timeout_helpers.run_execution_with_blocking_exception_handler",
            args=(run_id,),
            kwargs={},
            job_id=assignment.job_id,
            job_timeout=3,
            result_ttl=60,
            failure_ttl=60,
        )

        environment = dict(os.environ)
        python_path = os.pathsep.join(
            [str(REPOSITORY_ROOT / "src"), str(REPOSITORY_ROOT / "test")]
        )
        environment.update(
            {
                "ENCODE_PIPELINE_DATABASE_URL": configured.database_url,
                REDIS_URL_ENV: configured.redis_url,
                QUEUE_NAME_ENV: configured.queue_name,
                WORKSPACE_ROOT_ENV: str(configured.workspace_root),
                JOB_TIMEOUT_SECONDS_ENV: str(configured.job_timeout_seconds),
                ENTERED_MARKER_ENV: str(entered_marker),
                CAUGHT_MARKER_ENV: str(caught_marker),
                "PYTHONDONTWRITEBYTECODE": "1",
                "PYTHONPATH": python_path,
            }
        )
        started_at = time.monotonic()
        worker_process = subprocess.Popen(
            [sys.executable, "-m", "encode_pipeline.workers.cli", "--burst"],
            cwd=REPOSITORY_ROOT,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True,
        )
        try:
            stdout, stderr = worker_process.communicate(timeout=15)
        except subprocess.TimeoutExpired as exc:
            stdout, stderr = terminate_rq_worker(worker_process)
            raise AssertionError(
                f"DurableWorker timeout hung; stdout={stdout!r}; stderr={stderr!r}"
            ) from exc
        elapsed = time.monotonic() - started_at

        assert worker_process.returncode == 0, stderr
        assert elapsed < 12
        assert entered_marker.read_text(encoding="utf-8") == "entered\n"
        assert not caught_marker.exists()

        job = queue.fetch_job(assignment.job_id)
        assert job is not None
        job.refresh()
        assert job.is_failed
        failure = job.latest_result()
        assert failure is not None
        assert "WorkerHardTimeout" in (failure.exc_string or "")

        persistence = open_run_persistence(configured.database_url)
        try:
            run_service = create_default_run_service(
                registry=create_default_workflow_registry(),
                repository=persistence.repository,
            )
            record = run_service.get_run(run_id)
        finally:
            persistence.close()

        assert record.status.value == "failed"
        assert record.error is not None
        assert record.error.code == "RUN_WORKER_FAILED"
        assert record.error.context == {"reason_code": "WORKER_JOB_TIMEOUT"}

        with pytest.raises(ProcessLookupError):
            os.killpg(worker_process.pid, 0)
    finally:
        if worker_process is not None and worker_process.poll() is None:
            terminate_rq_worker(worker_process)
        if connected:
            job = queue.fetch_job(assignment.job_id)
            if job is not None:
                job.delete()
            connection.delete(queue.key)
        connection.close()
