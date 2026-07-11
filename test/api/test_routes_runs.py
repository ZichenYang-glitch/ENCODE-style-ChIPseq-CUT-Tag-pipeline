"""Tests for the FastAPI run lifecycle routes."""

from __future__ import annotations

import ast
from collections.abc import Iterator
from pathlib import Path

import pytest

from encode_pipeline.api.main import create_app
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.runs import RunStatus
from api_test_client import ApiTestClient

fastapi = pytest.importorskip("fastapi")


@pytest.fixture
def client() -> Iterator[ApiTestClient]:
    """Default app wired to the bundled ENCODE-style adapter."""
    app = create_app()
    with ApiTestClient(app) as tc:
        yield tc


@pytest.fixture
def workflow_id() -> str:
    return "encode-style-chipseq-cuttag-atac-mnase"


def _create_run(client: ApiTestClient, workflow_id: str) -> str:
    response = client.post(
        f"/api/v1/workflows/{workflow_id}/runs",
        json={"config": {"samples": "samples.tsv"}, "tags": {"env": "test"}},
    )
    assert response.status_code == 201
    return response.json()["run"]["run_id"]


def test_create_run_returns_201_with_run_response(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    response = client.post(
        f"/api/v1/workflows/{workflow_id}/runs",
        json={"config": {"samples": "samples.tsv"}, "tags": {"env": "test"}},
    )
    assert response.status_code == 201
    data = response.json()
    assert data["ok"] is True
    run = data["run"]
    assert run["workflow_id"] == workflow_id
    assert run["status"] == "created"
    assert run["tags"] == {"env": "test"}
    assert data["issues"] == []


def test_create_run_unknown_workflow_returns_404(client: ApiTestClient) -> None:
    response = client.post(
        "/api/v1/workflows/missing-workflow/runs",
        json={"config": {}},
    )
    assert response.status_code == 404
    data = response.json()
    assert data["ok"] is False
    assert data["run"] is None
    assert data["issues"][0]["code"] == "WORKFLOW_NOT_FOUND"


def test_create_run_malformed_workflow_id_returns_400(client: ApiTestClient) -> None:
    # FastAPI path parameters cannot be empty, so simulate malformed by relying
    # on registry normalization via a whitespace-only id that becomes empty.
    # This is the closest realistic path; the route catches ValueError.
    response = client.post(
        "/api/v1/workflows/%20/runs",
        json={"config": {}},
    )
    assert response.status_code == 400
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "API_REQUEST_INVALID"


def test_create_run_malformed_body_returns_400(
    client: ApiTestClient, workflow_id: str
) -> None:
    response = client.post(
        f"/api/v1/workflows/{workflow_id}/runs",
        content="not-json",
        headers={"content-type": "application/json"},
    )
    assert response.status_code == 400
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "API_REQUEST_INVALID"


def test_get_run_returns_200(client: ApiTestClient, workflow_id: str) -> None:
    run_id = _create_run(client, workflow_id)
    response = client.get(f"/api/v1/runs/{run_id}")
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert data["run"]["run_id"] == run_id


def test_get_run_unknown_returns_404(client: ApiTestClient) -> None:
    response = client.get("/api/v1/runs/run-missing")
    assert response.status_code == 404
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "RUN_NOT_FOUND"


def test_cancel_run_active_returns_200(client: ApiTestClient, workflow_id: str) -> None:
    service = client.app.state.run_service
    record = service.create_run(
        workflow_id, WorkflowInputs(config={"samples": "samples.tsv"})
    )
    run_id = record.run_id

    response = client.post(f"/api/v1/runs/{run_id}/cancel")
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert data["run"]["status"] == "cancelled"


def test_cancel_run_created_returns_cancelled_unchanged(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    run_id = _create_run(client, workflow_id)
    response = client.post(f"/api/v1/runs/{run_id}/cancel")
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert data["run"]["status"] == "cancelled"
    response = client.post(f"/api/v1/runs/{run_id}/cancel")
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert data["run"]["status"] == "cancelled"


def test_cancel_run_running_returns_stable_409_without_mutation(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    service = client.app.state.run_service
    record = service.create_run(
        workflow_id,
        WorkflowInputs(config={"samples": "samples.tsv"}),
    )
    service.transition_run(record.run_id, RunStatus.VALIDATING)
    service.transition_run(record.run_id, RunStatus.PLANNED)
    service.transition_run(record.run_id, RunStatus.QUEUED)
    running = service.transition_run(record.run_id, RunStatus.RUNNING)
    events_before = service.list_events(record.run_id, limit=100)

    response = client.post(f"/api/v1/runs/{record.run_id}/cancel")

    assert response.status_code == 409
    data = response.json()
    assert data["ok"] is False
    assert data["run"]["status"] == "running"
    assert data["run"]["ended_at"] is None
    assert data["run"]["cancellation_reason"] is None
    assert data["issues"][0] == {
        "code": "RUN_CANCELLATION_NOT_AVAILABLE",
        "message": "Run cancellation is not available in the current state.",
        "severity": "error",
        "path": "run_id",
        "source": "api",
        "technical_message": None,
        "hint": None,
        "context": {"current_status": "running"},
    }
    assert service.get_run(record.run_id) == running
    assert service.list_events(record.run_id, limit=100) == events_before


def test_cancel_run_openapi_declares_not_found_and_conflict_responses(
    client: ApiTestClient,
) -> None:
    responses = client.app.openapi()["paths"]["/api/v1/runs/{run_id}/cancel"]["post"][
        "responses"
    ]

    for status_code in ("404", "409"):
        assert responses[status_code]["content"]["application/json"]["schema"] == {
            "$ref": "#/components/schemas/RunResponse"
        }


def test_cancel_run_unknown_returns_404(client: ApiTestClient) -> None:
    response = client.post("/api/v1/runs/run-missing/cancel")
    assert response.status_code == 404
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "RUN_NOT_FOUND"


def test_list_events_returns_200(client: ApiTestClient, workflow_id: str) -> None:
    run_id = _create_run(client, workflow_id)
    response = client.get(f"/api/v1/runs/{run_id}/events")
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert data["run_id"] == run_id
    assert len(data["events"]) == 1
    assert data["next_cursor"] is None


def test_list_events_pagination_next_cursor(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    run_id = _create_run(client, workflow_id)
    service = client.app.state.run_service
    service.add_event(run_id, "stage_started", "Stage 1")
    service.add_event(run_id, "stage_completed", "Stage 2")

    response = client.get(f"/api/v1/runs/{run_id}/events?limit=1")
    assert response.status_code == 200
    data = response.json()
    assert len(data["events"]) == 1
    assert data["next_cursor"] == data["events"][-1]["event_id"]

    response = client.get(
        f"/api/v1/runs/{run_id}/events?after={data['next_cursor']}&limit=1"
    )
    assert response.status_code == 200
    data = response.json()
    assert len(data["events"]) == 1
    assert data["events"][0]["message"] == "Stage 1"


def test_list_events_unknown_run_returns_404(client: ApiTestClient) -> None:
    response = client.get("/api/v1/runs/run-missing/events")
    assert response.status_code == 404
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "RUN_NOT_FOUND"


def test_list_events_invalid_cursor_returns_400(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    run_id = _create_run(client, workflow_id)
    response = client.get(f"/api/v1/runs/{run_id}/events?after=evt-missing")
    assert response.status_code == 400
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "RUN_CURSOR_NOT_FOUND"


def test_list_logs_returns_200(client: ApiTestClient, workflow_id: str) -> None:
    run_id = _create_run(client, workflow_id)
    service = client.app.state.run_service
    service.append_log(run_id, "stdout", ["line 1"])

    response = client.get(f"/api/v1/runs/{run_id}/logs")
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert data["run_id"] == run_id
    assert data["stream_name"] == "stdout"
    assert len(data["chunks"]) == 1
    assert data["chunks"][-1]["lines"] == ["line 1"]
    assert data["next_cursor"] is None


def test_list_logs_pagination_next_cursor(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    run_id = _create_run(client, workflow_id)
    service = client.app.state.run_service
    service.append_log(run_id, "stdout", ["line 1"])
    service.append_log(run_id, "stdout", ["line 2"])

    response = client.get(f"/api/v1/runs/{run_id}/logs?limit=1")
    assert response.status_code == 200
    data = response.json()
    assert len(data["chunks"]) == 1
    assert data["next_cursor"] == data["chunks"][0]["chunk_id"]


def test_list_logs_unknown_run_returns_404(client: ApiTestClient) -> None:
    response = client.get("/api/v1/runs/run-missing/logs")
    assert response.status_code == 404
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "RUN_NOT_FOUND"


def test_list_logs_invalid_cursor_returns_400(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    run_id = _create_run(client, workflow_id)
    response = client.get(f"/api/v1/runs/{run_id}/logs?after=log-missing")
    assert response.status_code == 400
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "RUN_CURSOR_NOT_FOUND"


def test_list_events_invalid_limit_returns_400(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    run_id = _create_run(client, workflow_id)
    response = client.get(f"/api/v1/runs/{run_id}/events?limit=0")
    assert response.status_code == 400
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "API_REQUEST_INVALID"


def test_api_routes_runs_import_boundary() -> None:
    """Run routes must not import execution or engine internals."""
    source_path = (
        Path(__file__).resolve().parents[2] / "src/encode_pipeline/api/routes/runs.py"
    )
    source = source_path.read_text(encoding="utf-8")
    tree = ast.parse(source)

    forbidden_modules = {
        "encode_pipeline.config.validator",
        "encode_pipeline.samples",
        "snakemake",
        "subprocess",
        "openai",
    }
    imported_modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_modules.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imported_modules.add(node.module)

    assert not any(
        module == forbidden or module.startswith(f"{forbidden}.")
        for module in imported_modules
        for forbidden in forbidden_modules
    )
    assert "encode_pipeline.config.validator" not in source
    assert "encode_pipeline.samples.load" not in source


def test_create_run_returns_created_with_default_driver(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    response = client.post(
        f"/api/v1/workflows/{workflow_id}/runs",
        json={"config": {"samples": "samples.tsv"}},
    )
    assert response.status_code == 201
    data = response.json()
    assert data["ok"] is True
    assert data["run"]["status"] == "created"


def test_create_run_generates_created_event(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    response = client.post(
        f"/api/v1/workflows/{workflow_id}/runs",
        json={"config": {"samples": "samples.tsv"}},
    )
    run_id = response.json()["run"]["run_id"]

    response = client.get(f"/api/v1/runs/{run_id}/events")
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    statuses = [event["status"] for event in data["events"]]
    assert statuses == ["created"]
