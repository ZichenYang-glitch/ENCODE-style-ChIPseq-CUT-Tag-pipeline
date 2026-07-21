"""Tests for the FastAPI run lifecycle routes."""

from __future__ import annotations

import ast
import asyncio
from collections.abc import Iterator
from dataclasses import replace
from datetime import timedelta
from pathlib import Path
from threading import Event, Thread

import fastapi.routing
import pytest
from sqlalchemy import text

from encode_pipeline.api.main import create_app
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.results import Result
from encode_pipeline.platform.runs import RunStatus
from api_test_client import ApiTestClient

fastapi = pytest.importorskip("fastapi")


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


@pytest.fixture
def client() -> Iterator[ApiTestClient]:
    """Default app wired to the bundled ENCODE-style adapter."""
    app = create_app()
    with ApiTestClient(app) as tc:
        yield tc


@pytest.fixture
def workflow_id() -> str:
    return "encode-style-chipseq-cuttag-atac-mnase"


def _valid_inline_inputs() -> dict:
    return {
        "config": {},
        "samples": [
            {
                "sample": "S1",
                "fastq_1": "/tmp/S1.R1.fastq.gz",
                "layout": "SE",
                "assay": "chipseq",
                "target": "CTCF",
                "peak_mode": "narrow",
                "genome": "hs",
                "bowtie2_index": "/tmp/indices/hs",
            }
        ],
        "options": {},
    }


def _create_snapshot(
    client: ApiTestClient,
    workflow_id: str,
    inputs: dict | None = None,
) -> str:
    response = client.post(
        f"/api/v1/workflows/{workflow_id}/validate",
        json=inputs or _valid_inline_inputs(),
    )
    assert response.status_code == 200
    assert response.json()["ok"] is True
    return response.json()["snapshot"]["snapshot_id"]


def _create_run_response(
    client: ApiTestClient,
    workflow_id: str,
    *,
    inputs: dict | None = None,
    tags: dict[str, str] | None = None,
):
    snapshot_id = _create_snapshot(client, workflow_id, inputs)
    return client.post(
        f"/api/v1/workflows/{workflow_id}/runs",
        json={"snapshot_id": snapshot_id, "tags": tags or {}},
    )


def _create_run(client: ApiTestClient, workflow_id: str) -> str:
    response = _create_run_response(
        client,
        workflow_id,
        tags={"env": "test"},
    )
    assert response.status_code == 201
    return response.json()["run"]["run_id"]


def test_create_run_returns_201_with_run_response(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    response = _create_run_response(
        client,
        workflow_id,
        tags={"env": "test"},
    )
    assert response.status_code == 201
    data = response.json()
    assert data["ok"] is True
    run = data["run"]
    assert run["workflow_id"] == workflow_id
    assert run["status"] == "created"
    assert run["tags"] == {"env": "test"}
    assert data["issues"] == []


def test_create_run_persists_inline_sample_rows_without_server_path(
    client: ApiTestClient,
    workflow_id: str,
    tmp_path: Path,
) -> None:
    row = {
        "sample": "S1",
        "fastq_1": str((tmp_path / "S1.R1.fastq.gz").resolve()),
        "layout": "SE",
        "assay": "chipseq",
        "target": "CTCF",
        "peak_mode": "narrow",
        "genome": "hs",
        "bowtie2_index": str((tmp_path / "indices/hs").resolve()),
    }

    response = _create_run_response(
        client,
        workflow_id,
        inputs={"config": {}, "samples": [row], "options": {}},
    )

    assert response.status_code == 201
    run = response.json()["run"]
    assert run["inputs"] == {"config": {}, "samples": [row], "options": {}}
    assert "encode-platform-inline-samples" not in response.text


def test_create_run_rejects_oversized_authoring_body_without_persisting_run(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    before = client.app.state.run_service.list_runs()
    oversized = "x" * (2 * 1024 * 1024 + 1)

    response = client.post(
        f"/api/v1/workflows/{workflow_id}/runs",
        json={"config": {"oversized": oversized}},
    )

    assert response.status_code == 413
    data = response.json()
    assert data["ok"] is False
    assert data["run"] is None
    assert data["issues"][0]["code"] == "API_REQUEST_TOO_LARGE"
    assert client.app.state.run_service.list_runs() == before
    assert oversized[:128] not in response.text


def test_create_run_unknown_workflow_returns_404(client: ApiTestClient) -> None:
    snapshot_id = _create_snapshot(
        client,
        "encode-style-chipseq-cuttag-atac-mnase",
    )
    response = client.post(
        "/api/v1/workflows/missing-workflow/runs",
        json={"snapshot_id": snapshot_id},
    )
    assert response.status_code == 404
    data = response.json()
    assert data["ok"] is False
    assert data["run"] is None
    assert data["issues"][0]["code"] == "VALIDATED_SNAPSHOT_NOT_FOUND"


def test_create_run_cross_workflow_path_does_not_disclose_snapshot(client) -> None:
    snapshot_id = _create_snapshot(
        client,
        "encode-style-chipseq-cuttag-atac-mnase",
    )
    response = client.post(
        "/api/v1/workflows/%20/runs",
        json={"snapshot_id": snapshot_id},
    )
    assert response.status_code == 404
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "VALIDATED_SNAPSHOT_NOT_FOUND"


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


def test_create_run_rejects_mutable_inputs_and_client_validation_claim(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    response = client.post(
        f"/api/v1/workflows/{workflow_id}/runs",
        json={"config": {}, "samples": [], "validated": True},
    )

    assert response.status_code == 400
    assert response.json()["issues"][0]["code"] == "API_REQUEST_INVALID"
    assert client.app.state.run_service.list_runs() == ()


def test_create_run_identical_snapshot_replay_returns_same_canonical_run(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    snapshot_id = _create_snapshot(client, workflow_id)
    request = {"snapshot_id": snapshot_id, "tags": {"owner": "lab"}}

    first = client.post(f"/api/v1/workflows/{workflow_id}/runs", json=request)
    replay = client.post(f"/api/v1/workflows/{workflow_id}/runs", json=request)

    assert first.status_code == 201
    assert replay.status_code == 200
    assert replay.json()["run"] == first.json()["run"]
    assert len(client.app.state.run_service.list_runs()) == 1


def test_create_run_snapshot_replay_with_different_tags_is_conflict(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    snapshot_id = _create_snapshot(client, workflow_id)
    first = client.post(
        f"/api/v1/workflows/{workflow_id}/runs",
        json={"snapshot_id": snapshot_id, "tags": {"owner": "lab"}},
    )
    conflict = client.post(
        f"/api/v1/workflows/{workflow_id}/runs",
        json={"snapshot_id": snapshot_id, "tags": {"owner": "other"}},
    )

    assert first.status_code == 201
    assert conflict.status_code == 409
    assert conflict.json()["issues"][0]["code"] == (
        "VALIDATED_SNAPSHOT_REPLAY_CONFLICT"
    )
    assert len(client.app.state.run_service.list_runs()) == 1


def test_create_run_expired_snapshot_fails_without_run(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    from encode_pipeline.services.validated_inputs import ValidatedRunCreationService

    snapshot_id = _create_snapshot(client, workflow_id)
    snapshot = client.app.state.run_service.get_validated_input_snapshot(snapshot_id)
    client.app.state.validated_run_creation_service = ValidatedRunCreationService(
        run_service=client.app.state.run_service,
        build_identity_provider=client.app.state.build_identity_provider,
        clock=lambda: snapshot.expires_at,
    )

    response = client.post(
        f"/api/v1/workflows/{workflow_id}/runs",
        json={"snapshot_id": snapshot_id},
    )

    assert response.status_code == 409
    assert response.json()["issues"][0]["code"] == "VALIDATED_SNAPSHOT_EXPIRED"
    assert client.app.state.run_service.list_runs() == ()


def test_create_run_stale_build_identity_fails_without_run(
    client: ApiTestClient,
    workflow_id: str,
    monkeypatch,
) -> None:
    snapshot_id = _create_snapshot(client, workflow_id)
    provider = client.app.state.build_identity_provider
    current = provider.capture(workflow_id).value
    assert current is not None
    monkeypatch.setattr(
        provider,
        "capture",
        lambda requested_workflow_id: Result.success(
            replace(
                current,
                workflow_id=requested_workflow_id,
                digest="b" * 64,
                captured_at=current.captured_at + timedelta(seconds=1),
            )
        ),
    )

    response = client.post(
        f"/api/v1/workflows/{workflow_id}/runs",
        json={"snapshot_id": snapshot_id},
    )

    assert response.status_code == 409
    assert response.json()["issues"][0]["code"] == "VALIDATED_SNAPSHOT_STALE"
    assert client.app.state.run_service.list_runs() == ()


def test_create_run_bulk_runtime_unavailable_fails_closed_without_run(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    source_id = _create_snapshot(client, workflow_id)
    source = client.app.state.run_service.get_validated_input_snapshot(source_id)
    bulk_adapter = client.app.state.registry.get("bulk-rnaseq")
    bulk_identity = replace(
        source.workflow_build_identity,
        workflow_id="bulk-rnaseq",
        adapter_version=bulk_adapter.metadata.version,
    )
    bulk_snapshot = replace(
        source,
        snapshot_id="vsnap_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
        workflow_id="bulk-rnaseq",
        adapter_version=bulk_adapter.metadata.version,
        schema_version=bulk_adapter.schema().schema_version,
        workflow_build_identity=bulk_identity,
    )
    client.app.state.persistence.repository.create_validated_input_snapshot(
        bulk_snapshot
    )

    response = client.post(
        "/api/v1/workflows/bulk-rnaseq/runs",
        json={"snapshot_id": bulk_snapshot.snapshot_id},
    )

    assert response.status_code == 503
    assert response.json()["issues"][0]["code"] == ("WORKFLOW_EXECUTION_UNAVAILABLE")
    assert client.app.state.run_service.list_runs() == ()


def test_create_run_corrupt_snapshot_returns_safe_data_invalid_without_run(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    snapshot_id = _create_snapshot(client, workflow_id)
    with client.app.state.persistence.engine.begin() as connection:
        connection.execute(
            text(
                "UPDATE validated_input_snapshots SET canonical_payload = :payload "
                "WHERE snapshot_id = :snapshot_id"
            ),
            {
                "payload": '{"config":{"secret":"/private/value"}}',
                "snapshot_id": snapshot_id,
            },
        )

    response = client.post(
        f"/api/v1/workflows/{workflow_id}/runs",
        json={"snapshot_id": snapshot_id},
    )

    assert response.status_code == 409
    body = response.json()
    assert body["issues"][0]["code"] == "VALIDATED_SNAPSHOT_DATA_INVALID"
    assert "/private/value" not in response.text
    assert body["issues"][0]["technical_message"] is None
    assert client.app.state.run_service.list_runs() == ()


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


def test_cancel_run_running_without_durable_assignment_returns_stable_409(
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
        "code": "RUN_CANCELLATION_CONFLICT",
        "message": "Run cancellation conflicts with durable execution state.",
        "severity": "error",
        "path": "run_id",
        "source": "api",
        "technical_message": None,
        "hint": None,
        "context": {"current_status": "running"},
    }
    assert service.get_run(record.run_id) == running
    assert service.list_events(record.run_id, limit=100) == events_before


def test_cancel_run_openapi_declares_async_and_error_responses(
    client: ApiTestClient,
) -> None:
    responses = client.app.openapi()["paths"]["/api/v1/runs/{run_id}/cancel"]["post"][
        "responses"
    ]

    for status_code in ("202", "404", "409", "503"):
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
    response = _create_run_response(client, workflow_id)
    assert response.status_code == 201
    data = response.json()
    assert data["ok"] is True
    assert data["run"]["status"] == "created"


def test_create_run_generates_created_event(
    client: ApiTestClient,
    workflow_id: str,
) -> None:
    response = _create_run_response(client, workflow_id)
    run_id = response.json()["run"]["run_id"]

    response = client.get(f"/api/v1/runs/{run_id}/events")
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    statuses = [event["status"] for event in data["events"]]
    assert statuses == ["created"]
