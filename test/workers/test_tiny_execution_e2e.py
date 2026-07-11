"""Real Redis, API submission, worker process, and tiny Snakemake execution."""

from __future__ import annotations

import asyncio
import os
from pathlib import Path
import shutil
from uuid import uuid4

import pytest
import httpx
from redis import Redis
import yaml

from encode_pipeline.api.main import create_app
from encode_pipeline.persistence.runtime import open_run_persistence
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
    WorkerSettings,
)

from .process_helpers import run_burst_worker


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
PROFILE_ROOT = REPOSITORY_ROOT / "test" / "profiles" / "platform_worker_tiny"
TEST_REDIS_URL_ENV = "ENCODE_PIPELINE_TEST_REDIS_URL"
WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


def _request(app, method: str, url: str, **kwargs) -> httpx.Response:
    async def send() -> httpx.Response:
        transport = httpx.ASGITransport(app=app)
        async with httpx.AsyncClient(
            transport=transport,
            base_url="http://testserver",
            follow_redirects=True,
        ) as client:
            response = await client.request(method, url, **kwargs)
            await response.aread()
            return response

    return asyncio.run(send())


def _cleanup_redis_queue(
    *,
    redis_url: str,
    queue_name: str,
    database_url: str,
    workspace_root: Path,
    run_identity: dict[str, str | None],
) -> None:
    """Remove this test's unique queue and durable job even after assertions fail."""
    job_id: str | None = None
    run_id = run_identity["run_id"]
    if run_id is not None:
        persistence = open_run_persistence(database_url)
        try:
            run_service = create_default_run_service(
                registry=create_default_workflow_registry(),
                repository=persistence.repository,
            )
            assignment = run_service.get_execution_assignment(run_id)
            if assignment is not None:
                job_id = assignment.job_id
        finally:
            persistence.close()

    settings = WorkerSettings(
        database_url=database_url,
        redis_url=redis_url,
        queue_name=queue_name,
        workspace_root=workspace_root,
    )
    connection = Redis.from_url(redis_url)
    run_queue = RqRunQueue(settings, connection=connection)
    try:
        if job_id is not None:
            job = run_queue._queue.fetch_job(job_id)
            if job is not None:
                job.delete()
        run_queue._queue.delete()
    finally:
        connection.close()


def test_api_to_worker_executes_tiny_snakemake_and_persists_lifecycle(
    tmp_path,
    monkeypatch,
    request,
):
    redis_url = os.getenv(TEST_REDIS_URL_ENV)
    if redis_url is None:
        pytest.skip(f"{TEST_REDIS_URL_ENV} is not configured")
    if shutil.which("snakemake") is None:
        pytest.skip("snakemake is not available on PATH")

    queue_name = f"encode-pipeline-e2e-{uuid4().hex}"
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    workspace_root = tmp_path / "workspaces"
    run_identity: dict[str, str | None] = {"run_id": None}
    request.addfinalizer(
        lambda: _cleanup_redis_queue(
            redis_url=redis_url,
            queue_name=queue_name,
            database_url=database_url,
            workspace_root=workspace_root,
            run_identity=run_identity,
        )
    )
    monkeypatch.setenv("ENCODE_PIPELINE_DATABASE_URL", database_url)
    monkeypatch.setenv(REDIS_URL_ENV, redis_url)
    monkeypatch.setenv(QUEUE_NAME_ENV, queue_name)
    monkeypatch.setenv(WORKSPACE_ROOT_ENV, str(workspace_root))

    config = yaml.safe_load((PROFILE_ROOT / "config.yaml").read_text(encoding="utf-8"))
    samples_path = (PROFILE_ROOT / "samples.tsv").resolve()
    config["samples"] = str(samples_path)

    app = create_app(database_url=database_url, workspace_root=workspace_root)
    try:
        created_response = _request(
            app,
            "POST",
            f"/api/v1/workflows/{WORKFLOW_ID}/runs",
            json={
                "config": config,
                "samples": str(samples_path),
                "options": {},
            },
        )
        assert created_response.status_code == 201
        run_id = created_response.json()["run"]["run_id"]
        run_identity["run_id"] = run_id

        preflight_response = _request(
            app,
            "POST",
            f"/api/v1/runs/{run_id}/preflight",
        )
        assert preflight_response.status_code == 202
        assert (
            _request(app, "GET", f"/api/v1/runs/{run_id}").json()["run"]["status"]
            == "planned"
        )

        start_response = _request(
            app,
            "POST",
            f"/api/v1/runs/{run_id}/start",
        )
        assert start_response.status_code == 202
        assert start_response.json()["run"]["status"] == "queued"
    finally:
        app.state.run_queue.close()
        app.state.persistence.close()

    restarted_app = create_app(
        database_url=database_url,
        workspace_root=workspace_root,
    )
    try:
        assert restarted_app.state.recovered_run_ids == ()
        restarted = _request(
            restarted_app,
            "GET",
            f"/api/v1/runs/{run_id}",
        )
        assert restarted.json()["run"]["status"] == "queued"
    finally:
        restarted_app.state.run_queue.close()
        restarted_app.state.persistence.close()

    environment = dict(os.environ)
    environment.update(
        {
            "ENCODE_PIPELINE_DATABASE_URL": database_url,
            REDIS_URL_ENV: redis_url,
            QUEUE_NAME_ENV: queue_name,
            WORKSPACE_ROOT_ENV: str(workspace_root),
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONPATH": str(REPOSITORY_ROOT / "src"),
        }
    )
    completed = run_burst_worker(
        environment,
        cwd=REPOSITORY_ROOT,
        timeout_seconds=60,
    )
    assert completed.returncode == 0, completed.stderr

    persistence = open_run_persistence(database_url)
    try:
        run_service = create_default_run_service(
            registry=create_default_workflow_registry(),
            repository=persistence.repository,
        )
        record = run_service.get_run(run_id)
        assignment = run_service.get_execution_assignment(run_id)
        events = run_service.list_events(run_id, limit=100)
        stdout = run_service.list_logs(run_id, "stdout", limit=100)
        stderr = run_service.list_logs(run_id, "stderr", limit=100)
    finally:
        persistence.close()

    assert record.status is RunStatus.SUCCEEDED
    assert record.started_at is not None
    assert record.ended_at is not None
    assert assignment is not None
    assert assignment.dispatched_at is not None
    assert assignment.claimed_at is not None
    assert [event.status for event in events if event.status is not None][-3:] == [
        RunStatus.QUEUED,
        RunStatus.RUNNING,
        RunStatus.SUCCEEDED,
    ]
    persisted_output = "\n".join(
        line for chunk in (*stdout, *stderr) for line in chunk.lines
    )
    assert "localrule all:" in persisted_output
    assert "1 of 1 steps (100%) done" in persisted_output

    connection = Redis.from_url(redis_url)
    run_queue = RqRunQueue(
        WorkerSettings(
            database_url=database_url,
            redis_url=redis_url,
            queue_name=queue_name,
            workspace_root=workspace_root,
        ),
        connection=connection,
    )
    try:
        job = run_queue._queue.fetch_job(assignment.job_id)
        assert job is not None
        assert job.is_finished
    finally:
        connection.close()
