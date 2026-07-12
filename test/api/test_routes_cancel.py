"""API tests for truthful durable process cancellation."""

from __future__ import annotations

import asyncio
from collections.abc import Iterator
import inspect
import socket
from threading import Event, Thread
import time

import fastapi.routing
import pytest

from encode_pipeline.api.main import create_app
from encode_pipeline.api.routes.runs import cancel_run
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.run_cancellation import RunCancellationService
from encode_pipeline.services.run_queue import RunQueueStopUnavailableError
from encode_pipeline.workers.settings import (
    REDIS_API_READ_TIMEOUT_SECONDS_ENV,
    REDIS_CONNECT_TIMEOUT_SECONDS_ENV,
    REDIS_URL_ENV,
)
from api_test_client import ApiTestClient


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


async def _run_in_joined_test_thread(function, *args, **kwargs):
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
    monkeypatch.setattr(
        fastapi.routing,
        "run_in_threadpool",
        _run_in_joined_test_thread,
    )


class RecordingStopQueue:
    backend = "rq"

    def __init__(self, queue_name: str) -> None:
        self.queue_name = queue_name
        self.assignments: list[RunExecutionAssignment] = []
        self.error: Exception | None = None
        self.callback = None

    def request_stop(self, assignment: RunExecutionAssignment) -> None:
        self.assignments.append(assignment)
        if self.callback is not None:
            self.callback(assignment)
        if self.error is not None:
            raise self.error


@pytest.fixture
def client_and_queue(tmp_path) -> Iterator[tuple[ApiTestClient, RecordingStopQueue]]:
    app = create_app(database_url=f"sqlite:///{tmp_path / 'platform.db'}")
    queue = RecordingStopQueue(app.state.worker_settings.queue_name)
    app.state.run_cancellation_service = RunCancellationService(
        app.state.run_service,
        queue,
    )
    try:
        with ApiTestClient(app) as client:
            yield client, queue
    finally:
        app.state.run_queue.close()
        app.state.persistence.close()


def test_running_cancel_returns_202_but_stays_running_until_callback(
    client_and_queue,
):
    client, queue = client_and_queue
    assignment = _create_running_run(client, queue.queue_name)

    response = client.post(f"/api/v1/runs/{assignment.run_id}/cancel")

    assert response.status_code == 202
    body = response.json()
    assert body["ok"] is True
    assert body["run"]["status"] == "running"
    assert body["run"]["ended_at"] is None
    assert body["run"]["cancellation_reason"] is None
    assert queue.assignments == [
        client.app.state.run_service.get_execution_assignment(assignment.run_id)
    ]

    persisted = queue.assignments[0]
    acknowledged = client.app.state.run_service.acknowledge_execution_stop(
        persisted.run_id,
        job_id=persisted.job_id,
        backend=persisted.backend,
        queue_name=persisted.queue_name,
    )
    assert acknowledged.record.status is RunStatus.CANCELLED


def test_running_cancel_retries_stop_without_duplicate_events(client_and_queue):
    client, queue = client_and_queue
    assignment = _create_running_run(client, queue.queue_name)

    first = client.post(f"/api/v1/runs/{assignment.run_id}/cancel")
    second = client.post(f"/api/v1/runs/{assignment.run_id}/cancel")

    assert first.status_code == second.status_code == 202
    assert len(queue.assignments) == 2
    events = client.app.state.run_service.list_events(assignment.run_id, limit=100)
    assert [event.event_type for event in events].count("cancellation_requested") == 1


def test_running_cancel_returns_terminal_200_when_callback_wins_response_race(
    client_and_queue,
):
    client, queue = client_and_queue
    assignment = _create_running_run(client, queue.queue_name)
    service = client.app.state.run_service
    queue.callback = lambda current: service.acknowledge_execution_stop(
        current.run_id,
        job_id=current.job_id,
        backend=current.backend,
        queue_name=current.queue_name,
    )

    response = client.post(f"/api/v1/runs/{assignment.run_id}/cancel")

    assert response.status_code == 200
    assert response.json()["run"]["status"] == "cancelled"


def test_running_cancel_returns_sanitized_retryable_503_without_false_terminal(
    client_and_queue,
):
    client, queue = client_and_queue
    assignment = _create_running_run(client, queue.queue_name)
    queue.error = RunQueueStopUnavailableError(
        "redis://user:password@private-host:6379 is unavailable"
    )

    response = client.post(f"/api/v1/runs/{assignment.run_id}/cancel")

    assert response.status_code == 503
    body = response.json()
    assert body["ok"] is False
    assert body["run"]["status"] == "running"
    assert body["run"]["ended_at"] is None
    assert body["run"]["cancellation_reason"] is None
    assert body["issues"][0]["code"] == "RUN_CANCELLATION_UNAVAILABLE"
    assert body["issues"][0]["context"] == {
        "current_status": "running",
        "retryable": True,
    }
    assert "private-host" not in response.text
    persisted = client.app.state.run_service.get_execution_assignment(assignment.run_id)
    assert persisted is not None
    assert persisted.cancellation_requested_at is not None
    assert persisted.cancellation_acknowledged_at is None


def test_cancel_route_is_sync_so_queue_io_runs_off_event_loop(client_and_queue):
    _client, _queue = client_and_queue

    assert inspect.iscoroutinefunction(cancel_run) is False


def test_cancel_real_redis_read_timeout_is_bounded_and_sanitized(
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

    server = Thread(target=accept_without_replying, daemon=True)
    server.start()
    monkeypatch.setenv(REDIS_URL_ENV, f"redis://{redis_host}:{redis_port}/0")
    monkeypatch.setenv(REDIS_CONNECT_TIMEOUT_SECONDS_ENV, "0.2")
    monkeypatch.setenv(REDIS_API_READ_TIMEOUT_SECONDS_ENV, "0.15")
    app = create_app(database_url=f"sqlite:///{tmp_path / 'platform.db'}")
    try:
        with ApiTestClient(app) as client:
            assignment = _create_running_run(
                client,
                app.state.worker_settings.queue_name,
            )
            started_at = time.monotonic()
            response = client.post(f"/api/v1/runs/{assignment.run_id}/cancel")
            elapsed = time.monotonic() - started_at

        assert accepted.is_set()
        assert response.status_code == 503
        assert response.json()["run"]["status"] == "running"
        assert response.json()["issues"][0]["code"] == ("RUN_CANCELLATION_UNAVAILABLE")
        assert redis_host not in response.text
        assert str(redis_port) not in response.text
        assert elapsed < 1.5
    finally:
        stop.set()
        listener.close()
        server.join(timeout=1)
        app.state.run_queue.close()
        app.state.persistence.close()


def _create_running_run(
    client: ApiTestClient,
    queue_name: str,
) -> RunExecutionAssignment:
    service = client.app.state.run_service
    record = service.create_run(WORKFLOW_ID, WorkflowInputs(config={}))
    service.transition_run(record.run_id, RunStatus.VALIDATING)
    service.transition_run(record.run_id, RunStatus.PLANNED)
    assignment = service.ensure_execution_assignment(
        record.run_id,
        queue_name=queue_name,
    )
    service.mark_execution_dispatched(record.run_id, job_id=assignment.job_id)
    service.transition_run(record.run_id, RunStatus.QUEUED)
    service.claim_execution_assignment(
        record.run_id,
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
    )
    service.transition_run(record.run_id, RunStatus.RUNNING)
    claimed = service.get_execution_assignment(record.run_id)
    assert claimed is not None
    return claimed
