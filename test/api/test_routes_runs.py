"""Tests for the FastAPI run lifecycle routes."""

from __future__ import annotations

import ast
from collections.abc import Iterator
from pathlib import Path

import pytest

from encode_pipeline.api.main import create_app
from encode_pipeline.platform.adapters import WorkflowInputs

fastapi = pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402


@pytest.fixture
def client() -> Iterator[TestClient]:
    """Default app wired to the bundled ENCODE-style adapter."""
    app = create_app()
    with TestClient(app) as tc:
        yield tc


@pytest.fixture
def workflow_id() -> str:
    return "encode-style-chipseq-cuttag-atac-mnase"


def _create_run(client: TestClient, workflow_id: str) -> str:
    response = client.post(
        f"/api/v1/workflows/{workflow_id}/runs",
        json={"config": {"samples": "samples.tsv"}, "tags": {"env": "test"}},
    )
    assert response.status_code == 201
    return response.json()["run"]["run_id"]


def test_create_run_returns_201_with_run_response(
    client: TestClient,
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
    assert run["status"] == "succeeded"
    assert run["tags"] == {"env": "test"}
    assert data["issues"] == []


def test_create_run_unknown_workflow_returns_404(client: TestClient) -> None:
    response = client.post(
        "/api/v1/workflows/missing-workflow/runs",
        json={"config": {}},
    )
    assert response.status_code == 404
    data = response.json()
    assert data["ok"] is False
    assert data["run"] is None
    assert data["issues"][0]["code"] == "WORKFLOW_NOT_FOUND"


def test_create_run_malformed_workflow_id_returns_400(client: TestClient) -> None:
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


def test_create_run_malformed_body_returns_400(client: TestClient, workflow_id: str) -> None:
    response = client.post(
        f"/api/v1/workflows/{workflow_id}/runs",
        data="not-json",
        headers={"content-type": "application/json"},
    )
    assert response.status_code == 400
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "API_REQUEST_INVALID"


def test_get_run_returns_200(client: TestClient, workflow_id: str) -> None:
    run_id = _create_run(client, workflow_id)
    response = client.get(f"/api/v1/runs/{run_id}")
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert data["run"]["run_id"] == run_id


def test_get_run_unknown_returns_404(client: TestClient) -> None:
    response = client.get("/api/v1/runs/run-missing")
    assert response.status_code == 404
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "RUN_NOT_FOUND"


def test_cancel_run_active_returns_200(client: TestClient, workflow_id: str) -> None:
    service = client.app.state.run_service
    record = service.create_run(workflow_id, WorkflowInputs(config={"samples": "samples.tsv"}))
    run_id = record.run_id

    response = client.post(f"/api/v1/runs/{run_id}/cancel")
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert data["run"]["status"] == "cancelled"


def test_cancel_run_terminal_returns_200_unchanged(
    client: TestClient,
    workflow_id: str,
) -> None:
    run_id = _create_run(client, workflow_id)
    response = client.post(f"/api/v1/runs/{run_id}/cancel")
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert data["run"]["status"] == "succeeded"
    response = client.post(f"/api/v1/runs/{run_id}/cancel")
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert data["run"]["status"] == "succeeded"


def test_cancel_run_unknown_returns_404(client: TestClient) -> None:
    response = client.post("/api/v1/runs/run-missing/cancel")
    assert response.status_code == 404
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "RUN_NOT_FOUND"


def test_list_events_returns_200(client: TestClient, workflow_id: str) -> None:
    run_id = _create_run(client, workflow_id)
    response = client.get(f"/api/v1/runs/{run_id}/events")
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert data["run_id"] == run_id
    assert len(data["events"]) == 6
    assert data["next_cursor"] is None


def test_list_events_pagination_next_cursor(
    client: TestClient,
    workflow_id: str,
) -> None:
    run_id = _create_run(client, workflow_id)
    service = client.app.state.run_service
    service.add_event(run_id, "stage_started", "Stage 1")
    service.add_event(run_id, "stage_completed", "Stage 2")

    response = client.get(f"/api/v1/runs/{run_id}/events?limit=6")
    assert response.status_code == 200
    data = response.json()
    assert len(data["events"]) == 6
    assert data["next_cursor"] == data["events"][-1]["event_id"]

    response = client.get(f"/api/v1/runs/{run_id}/events?after={data['next_cursor']}&limit=1")
    assert response.status_code == 200
    data = response.json()
    assert len(data["events"]) == 1
    assert data["events"][0]["message"] == "Stage 1"


def test_list_events_unknown_run_returns_404(client: TestClient) -> None:
    response = client.get("/api/v1/runs/run-missing/events")
    assert response.status_code == 404
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "RUN_NOT_FOUND"


def test_list_events_invalid_cursor_returns_400(
    client: TestClient,
    workflow_id: str,
) -> None:
    run_id = _create_run(client, workflow_id)
    response = client.get(f"/api/v1/runs/{run_id}/events?after=evt-missing")
    assert response.status_code == 400
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "RUN_CURSOR_NOT_FOUND"


def test_list_logs_returns_200(client: TestClient, workflow_id: str) -> None:
    run_id = _create_run(client, workflow_id)
    service = client.app.state.run_service
    service.append_log(run_id, "stdout", ["line 1"])

    response = client.get(f"/api/v1/runs/{run_id}/logs")
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert data["run_id"] == run_id
    assert data["stream_name"] == "stdout"
    assert len(data["chunks"]) == 6
    assert data["chunks"][-1]["lines"] == ["line 1"]
    assert data["next_cursor"] is None


def test_list_logs_pagination_next_cursor(
    client: TestClient,
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


def test_list_logs_unknown_run_returns_404(client: TestClient) -> None:
    response = client.get("/api/v1/runs/run-missing/logs")
    assert response.status_code == 404
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "RUN_NOT_FOUND"


def test_list_logs_invalid_cursor_returns_400(
    client: TestClient,
    workflow_id: str,
) -> None:
    run_id = _create_run(client, workflow_id)
    response = client.get(f"/api/v1/runs/{run_id}/logs?after=log-missing")
    assert response.status_code == 400
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "RUN_CURSOR_NOT_FOUND"


def test_list_events_invalid_limit_returns_400(
    client: TestClient,
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
    source_path = Path(__file__).resolve().parents[2] / "src/encode_pipeline/api/routes/runs.py"
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
