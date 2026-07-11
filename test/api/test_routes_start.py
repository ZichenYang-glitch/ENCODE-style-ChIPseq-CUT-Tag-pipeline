"""Tests for explicit durable run submission routes."""

from __future__ import annotations

import asyncio
from collections.abc import Iterator
import socket
from threading import Event, Thread, Timer
import time

import fastapi.routing
import httpx
import pytest

from encode_pipeline.api.main import create_app
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.run_queue import RunQueueUnavailableError
from encode_pipeline.services.run_submission import RunSubmissionService
from encode_pipeline.workers.settings import (
    REDIS_API_READ_TIMEOUT_SECONDS_ENV,
    REDIS_CONNECT_TIMEOUT_SECONDS_ENV,
    REDIS_URL_ENV,
)
from api_test_client import ApiTestClient


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


async def _run_in_joined_test_thread(function, *args, **kwargs):
    """Exercise sync FastAPI routes without leaking Python 3.13 executor threads."""
    completed = Event()
    results: list[object] = []
    exceptions: list[BaseException] = []

    def invoke() -> None:
        try:
            results.append(function(*args, **kwargs))
        except BaseException as exc:
            exceptions.append(exc)
        finally:
            completed.set()

    thread = Thread(target=invoke)
    thread.start()
    try:
        while not completed.is_set():
            await asyncio.sleep(0.001)
    finally:
        thread.join(timeout=3)
        if thread.is_alive():  # pragma: no cover - test deadlock guard
            raise RuntimeError("test threadpool call did not terminate")
    if exceptions:
        raise exceptions[0]
    return results[0]


@pytest.fixture(autouse=True)
def joined_test_threadpool(monkeypatch):
    """Use a one-shot joined thread for this environment's ASGI test client."""
    monkeypatch.setattr(
        fastapi.routing,
        "run_in_threadpool",
        _run_in_joined_test_thread,
    )


class RecordingRunQueue:
    backend = "rq"
    queue_name = "api-tests"

    def __init__(self) -> None:
        self.assignments: list[RunExecutionAssignment] = []
        self.error: Exception | None = None

    def enqueue_execution(self, assignment: RunExecutionAssignment) -> str:
        self.assignments.append(assignment)
        if self.error is not None:
            raise self.error
        return assignment.job_id


@pytest.fixture
def client_and_queue(tmp_path) -> Iterator[tuple[ApiTestClient, RecordingRunQueue]]:
    app = create_app(database_url=f"sqlite:///{tmp_path / 'platform.db'}")
    queue = RecordingRunQueue()
    app.state.run_submission_service = RunSubmissionService(
        app.state.run_service,
        queue,
    )
    try:
        with ApiTestClient(app) as client:
            yield client, queue
    finally:
        try:
            app.state.run_queue.close()
        finally:
            app.state.persistence.close()


def _create_run(
    client: ApiTestClient,
    *,
    planned: bool,
    bind_build_identity: bool = True,
) -> str:
    service = client.app.state.run_service
    record = service.create_run(WORKFLOW_ID, WorkflowInputs(config={}))
    if planned:
        service.transition_run(record.run_id, RunStatus.VALIDATING)
        if bind_build_identity:
            identity_result = client.app.state.build_identity_provider.capture(
                WORKFLOW_ID
            )
            assert identity_result.is_success
            service.complete_preflight(record.run_id, identity_result.value)
        else:
            service.transition_run(record.run_id, RunStatus.PLANNED)
    return record.run_id


def test_start_run_returns_202_and_durable_queued_state(client_and_queue):
    client, queue = client_and_queue
    run_id = _create_run(client, planned=True)

    response = client.post(f"/api/v1/runs/{run_id}/start")
    repeated = client.post(f"/api/v1/runs/{run_id}/start")

    assert response.status_code == 202
    assert response.json()["run"]["status"] == "queued"
    assert repeated.status_code == 202
    assert repeated.json()["run"]["status"] == "queued"
    assert len(queue.assignments) == 2
    assert queue.assignments[0].job_id == queue.assignments[1].job_id
    events = client.get(f"/api/v1/runs/{run_id}/events").json()["events"]
    assert len([event for event in events if event["status"] == "queued"]) == 1


def test_start_run_returns_409_until_preflight_completes(client_and_queue):
    client, queue = client_and_queue
    run_id = _create_run(client, planned=False)

    response = client.post(f"/api/v1/runs/{run_id}/start")

    assert response.status_code == 409
    body = response.json()
    assert body["ok"] is False
    assert body["run"]["status"] == "created"
    assert body["issues"][0]["code"] == "RUN_NOT_READY"
    assert body["issues"][0]["context"] == {"current_status": "created"}
    assert queue.assignments == []


def test_start_run_rejects_legacy_planned_run_without_build_identity(
    client_and_queue,
):
    client, queue = client_and_queue
    run_id = _create_run(client, planned=True, bind_build_identity=False)

    response = client.post(f"/api/v1/runs/{run_id}/start")

    assert response.status_code == 409
    body = response.json()
    assert body["run"]["status"] == "planned"
    assert body["issues"][0]["code"] == ("RUN_WORKFLOW_BUILD_IDENTITY_MISSING")
    assert body["issues"][0]["context"] == {"current_status": "planned"}
    assert queue.assignments == []


def test_start_run_returns_404_for_missing_run(client_and_queue):
    client, _queue = client_and_queue

    response = client.post("/api/v1/runs/missing/start")

    assert response.status_code == 404
    assert response.json()["issues"][0]["code"] == "RUN_NOT_FOUND"


def test_start_run_returns_sanitized_503_and_keeps_planned(client_and_queue):
    client, queue = client_and_queue
    run_id = _create_run(client, planned=True)
    queue.error = RunQueueUnavailableError(
        "redis://user:password@private-host:6379 is unavailable"
    )

    response = client.post(f"/api/v1/runs/{run_id}/start")

    assert response.status_code == 503
    body = response.json()
    assert body["run"]["status"] == "planned"
    assert body["issues"][0]["code"] == "RUN_QUEUE_UNAVAILABLE"
    assert body["issues"][0]["context"] == {
        "current_status": "planned",
        "retryable": True,
    }
    assert "private-host" not in response.text
    assignment = client.app.state.run_service.get_execution_assignment(run_id)
    assert assignment is not None
    assert assignment.dispatched_at is None


def test_start_run_submission_does_not_block_the_api_event_loop(tmp_path):
    app = create_app(database_url=f"sqlite:///{tmp_path / 'platform.db'}")
    service = app.state.run_service
    record = service.create_run(WORKFLOW_ID, WorkflowInputs(config={}))
    service.transition_run(record.run_id, RunStatus.VALIDATING)
    planned = service.transition_run(record.run_id, RunStatus.PLANNED)
    entered = Event()
    release = Event()

    class BlockingSubmissionService:
        def start_run(self, run_id):
            assert run_id == planned.run_id
            entered.set()
            if not release.wait(timeout=2):  # pragma: no cover - test deadlock guard
                raise RuntimeError("submission test gate timed out")
            return planned

    app.state.run_submission_service = BlockingSubmissionService()

    @app.get("/event-loop-probe")
    async def event_loop_probe():
        return {"ok": True}

    async def exercise() -> tuple[httpx.Response, httpx.Response, float]:
        transport = httpx.ASGITransport(app=app)
        async with httpx.AsyncClient(
            transport=transport,
            base_url="http://testserver",
        ) as client:
            started_at = time.monotonic()
            start_task = asyncio.create_task(
                client.post(f"/api/v1/runs/{planned.run_id}/start")
            )
            entered_deadline = time.monotonic() + 2
            while not entered.is_set() and time.monotonic() < entered_deadline:
                await asyncio.sleep(0.001)
            assert entered.is_set()
            probe_response = await client.get("/event-loop-probe")
            probe_elapsed = time.monotonic() - started_at
            release.set()
            start_response = await start_task
            return start_response, probe_response, probe_elapsed

    # A regression to an async route doing synchronous submission blocks the
    # event loop until this independent thread releases the service gate.
    deadlock_guard = Timer(1.0, release.set)
    deadlock_guard.start()
    try:
        start_response, probe_response, probe_elapsed = asyncio.run(exercise())
    finally:
        release.set()
        deadlock_guard.cancel()
        app.state.run_queue.close()
        app.state.persistence.close()

    assert start_response.status_code == 202
    assert probe_response.json() == {"ok": True}
    assert probe_elapsed < 0.75


def test_start_run_real_redis_read_timeout_is_bounded_and_sanitized(
    tmp_path,
    monkeypatch,
):
    listener = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    listener.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    listener.bind(("127.0.0.1", 0))
    listener.listen()
    listener.settimeout(0.05)
    redis_host, redis_port = listener.getsockname()
    stop = Event()
    accepted = Event()

    def accept_without_replying() -> None:
        connections: list[socket.socket] = []
        try:
            while not stop.is_set():
                try:
                    connection, _address = listener.accept()
                except TimeoutError:
                    continue
                except OSError:
                    break
                connections.append(connection)
                accepted.set()
        finally:
            for connection in connections:
                connection.close()

    server_thread = Thread(target=accept_without_replying, daemon=True)
    server_thread.start()
    monkeypatch.setenv(REDIS_URL_ENV, f"redis://{redis_host}:{redis_port}/0")
    monkeypatch.setenv(REDIS_CONNECT_TIMEOUT_SECONDS_ENV, "0.2")
    monkeypatch.setenv(REDIS_API_READ_TIMEOUT_SECONDS_ENV, "0.15")
    app = create_app(database_url=f"sqlite:///{tmp_path / 'platform.db'}")
    try:
        with ApiTestClient(app) as client:
            run_id = _create_run(client, planned=True)
            started_at = time.monotonic()
            response = client.post(f"/api/v1/runs/{run_id}/start")
            elapsed = time.monotonic() - started_at

        assert accepted.is_set()
        assert response.status_code == 503
        assert response.json()["run"]["status"] == "planned"
        assert response.json()["issues"][0]["code"] == "RUN_QUEUE_UNAVAILABLE"
        assert redis_host not in response.text
        assert str(redis_port) not in response.text
        assert elapsed < 1.5
    finally:
        stop.set()
        listener.close()
        server_thread.join(timeout=1)
        app.state.run_queue.close()
        app.state.persistence.close()
