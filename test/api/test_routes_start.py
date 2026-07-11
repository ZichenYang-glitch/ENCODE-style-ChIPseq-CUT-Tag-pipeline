"""Tests for explicit durable run submission routes."""

from __future__ import annotations

from collections.abc import Iterator

import pytest

from encode_pipeline.api.main import create_app
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.run_queue import RunQueueUnavailableError
from encode_pipeline.services.run_submission import RunSubmissionService
from api_test_client import ApiTestClient


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


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


def _create_run(client: ApiTestClient, *, planned: bool) -> str:
    service = client.app.state.run_service
    record = service.create_run(WORKFLOW_ID, WorkflowInputs(config={}))
    if planned:
        service.transition_run(record.run_id, RunStatus.VALIDATING)
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
