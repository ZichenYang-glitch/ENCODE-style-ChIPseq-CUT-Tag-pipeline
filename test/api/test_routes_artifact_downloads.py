"""Tests for descriptor-safe artifact download API routes."""

from __future__ import annotations

from collections.abc import Iterator
from datetime import datetime, timezone
import asyncio
import ast
import inspect
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from encode_pipeline.api.main import create_app
from encode_pipeline.api.main import _handle_internal_server_error
from encode_pipeline.api.routes.artifacts import download_run_artifact
from encode_pipeline.api.routes import artifacts as artifact_routes
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.runs import RunArtifactRef, RunStatus
from encode_pipeline.platform.results import Issue, Result
from fastapi import Request
from fastapi.responses import JSONResponse, StreamingResponse
from starlette.requests import ClientDisconnect
from api_test_client import ApiTestClient


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
ARTIFACT_REVISION = f"artifactrev-{'0' * 64}"


@pytest.fixture
def client(tmp_path: Path) -> Iterator[ApiTestClient]:
    app = create_app(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        workspace_root=tmp_path / "workspaces",
    )
    with ApiTestClient(app) as test_client:
        yield test_client


def _create_run(client: ApiTestClient) -> str:
    service = client.app.state.run_service
    run_id = service.create_run(
        WORKFLOW_ID,
        WorkflowInputs(config={}),
    ).run_id
    for status in (
        RunStatus.VALIDATING,
        RunStatus.PLANNED,
        RunStatus.QUEUED,
        RunStatus.RUNNING,
        RunStatus.SUCCEEDED,
    ):
        service.transition_run(run_id, status)
    return run_id


def _record_download(
    client: ApiTestClient,
    run_id: str,
    artifact_id: str,
    content: bytes,
    *,
    name: str = "result manifest.tsv",
    mime_type: str = "text/tab-separated-values",
    metadata_size: int | None = None,
) -> RunArtifactRef:
    relative_path = f"results/multiqc/{name}"
    source = client.app.state.workspace_root / run_id / relative_path
    source.parent.mkdir(parents=True, exist_ok=True)
    source.write_bytes(content)
    artifact = RunArtifactRef(
        artifact_id=artifact_id,
        run_id=run_id,
        artifact_type="file",
        name=name,
        uri=f"run://runs/{run_id}/artifacts/{artifact_id}",
        mime_type=mime_type,
        produced_at=datetime.now(timezone.utc),
        revision=ARTIFACT_REVISION,
        metadata={
            "relative_path": relative_path,
            "output_type": "result_manifest",
            "size_bytes": len(content) if metadata_size is None else metadata_size,
            "catalog_id": "result_manifest",
            "scope": "project",
        },
    )
    service = client.app.state.run_service
    service.replace_artifacts(run_id, (*service.list_artifacts(run_id), artifact))
    return artifact


def _download_url(run_id: str, artifact_id: str) -> str:
    return f"/api/v1/runs/{run_id}/artifacts/{artifact_id}/download"


def _invoke(client: ApiTestClient, run_id: str, artifact_id: str):
    return download_run_artifact(
        run_id,
        artifact_id,
        client.app.state.artifact_download_service,
    )


def _json(response: JSONResponse) -> dict:
    return json.loads(response.body)


def test_download_streams_exact_persisted_bytes_with_controlled_headers(client):
    run_id = _create_run(client)
    content = b"output_type\tstatus\tpath\nresult_manifest\tpresent\tresults/a.tsv\n"
    artifact = _record_download(client, run_id, "artifact-manifest", content)
    events_before = client.app.state.run_service.list_events(run_id)

    response = _invoke(client, run_id, artifact.artifact_id)

    assert response.status_code == 200
    assert isinstance(response, StreamingResponse)
    assert (
        response.headers["content-type"] == "text/tab-separated-values; charset=utf-8"
    )
    assert response.headers["content-length"] == str(len(content))
    assert response.headers["x-content-type-options"] == "nosniff"
    assert response.headers["cache-control"] == "private, no-store"
    disposition = response.headers["content-disposition"]
    assert disposition.startswith("attachment; filename=")
    assert "filename*=UTF-8''result%20manifest.tsv" in disposition
    assert str(client.app.state.workspace_root) not in disposition
    assert client.app.state.run_service.list_events(run_id) == events_before
    assert response.background is not None
    response.background.func()


def test_unknown_and_cross_run_artifact_share_same_redacted_404(client):
    first_run = _create_run(client)
    second_run = _create_run(client)
    artifact = _record_download(client, first_run, "artifact-private", b"secret")
    _record_download(client, second_run, "artifact-present", b"present")

    missing = _invoke(client, second_run, "artifact-missing")
    cross_run = _invoke(client, second_run, artifact.artifact_id)

    assert missing.status_code == cross_run.status_code == 404
    missing_body = _json(missing)
    cross_body = _json(cross_run)
    assert missing_body["issues"][0]["code"] == "RUN_ARTIFACT_NOT_FOUND"
    assert cross_body["issues"][0]["code"] == "RUN_ARTIFACT_NOT_FOUND"
    assert missing_body["issues"] == cross_body["issues"]
    serialized = cross_run.body.decode()
    assert str(client.app.state.workspace_root) not in serialized
    assert "secret" not in serialized


def test_unknown_run_uses_stable_redacted_404(client):
    response = _invoke(client, "run-missing", "artifact-missing")

    assert response.status_code == 404
    body = _json(response)
    assert body["issues"][0]["code"] == "RUN_NOT_FOUND"
    assert "technical_message" not in body["issues"][0]


@pytest.mark.parametrize("conflict", ["missing", "symlink", "size"])
def test_changed_or_unsafe_source_returns_redacted_conflict(client, conflict):
    run_id = _create_run(client)
    artifact = _record_download(client, run_id, "artifact-source", b"indexed")
    source = (
        client.app.state.workspace_root / run_id / artifact.metadata["relative_path"]
    )
    if conflict == "missing":
        source.unlink()
    elif conflict == "symlink":
        source.unlink()
        outside = source.parents[3] / "outside-secret.txt"
        outside.write_text("outside secret", encoding="utf-8")
        source.symlink_to(outside)
    else:
        source.write_bytes(b"different size")

    response = _invoke(client, run_id, artifact.artifact_id)

    assert response.status_code == 409
    body = _json(response)
    assert body["issues"][0]["code"] == "RUN_ARTIFACT_DOWNLOAD_CONFLICT"
    assert body["issues"][0]["context"]["reason_code"].startswith("ARTIFACT_DOWNLOAD_")
    serialized = response.body.decode()
    assert str(client.app.state.workspace_root) not in serialized
    assert "outside secret" not in serialized
    assert "technical_message" not in body["issues"][0]


def test_corrupt_persisted_metadata_returns_redacted_500(client, monkeypatch):
    run_id = _create_run(client)
    artifact = _record_download(client, run_id, "artifact-corrupt", b"indexed")
    corrupted = RunArtifactRef(
        **{
            **artifact.__dict__,
            "uri": "file:///private/workspace/secret.txt",
        }
    )
    monkeypatch.setattr(
        client.app.state.run_service,
        "get_artifact",
        lambda requested_run_id, requested_artifact_id: (
            corrupted
            if (requested_run_id, requested_artifact_id)
            == (run_id, artifact.artifact_id)
            else (_ for _ in ()).throw(KeyError(requested_artifact_id))
        ),
    )

    response = _invoke(client, run_id, artifact.artifact_id)

    assert response.status_code == 500
    body = _json(response)
    assert body["issues"][0]["code"] == ("RUN_ARTIFACT_DOWNLOAD_DATA_INVALID")
    serialized = response.body.decode()
    assert "/private/workspace" not in serialized
    assert "technical_message" not in body["issues"][0]


def test_unexpected_service_failure_uses_download_specific_redacted_envelope(client):
    run_id = _create_run(client)
    request = Request(
        {
            "type": "http",
            "method": "GET",
            "path": _download_url(run_id, "artifact-broken"),
            "headers": [],
            "query_string": b"",
            "path_params": {
                "run_id": run_id,
                "artifact_id": "artifact-broken",
            },
            "route": SimpleNamespace(operation_id="downloadRunArtifact"),
        }
    )
    response = asyncio.run(
        _handle_internal_server_error(
            request,
            RuntimeError("DATABASE_URL=sqlite:////private/platform.db"),
        )
    )

    assert response.status_code == 500
    body = _json(response)
    assert body["ok"] is False
    assert body["run_id"] == run_id
    assert body["artifact_id"] == "artifact-broken"
    assert body["issues"][0]["code"] == "INTERNAL_SERVER_ERROR"
    serialized = response.body.decode()
    assert "DATABASE_URL" not in serialized
    assert "/private" not in serialized


def test_download_survives_sqlite_repository_reopen(tmp_path: Path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    workspace_root = tmp_path / "workspaces"
    first_app = create_app(database_url=database_url, workspace_root=workspace_root)
    with ApiTestClient(first_app) as first_client:
        run_id = _create_run(first_client)
        artifact = _record_download(
            first_client,
            run_id,
            "artifact-restart",
            b"durable download",
        )
    first_app.state.persistence.close()
    first_app.state.run_queue.close()

    second_app = create_app(database_url=database_url, workspace_root=workspace_root)
    try:
        with ApiTestClient(second_app) as second_client:
            result = second_client.app.state.artifact_download_service.prepare(
                run_id, artifact.artifact_id
            )
        assert result.is_success
        assert result.value is not None
        assert b"".join(result.value.iter_bytes()) == b"durable download"
    finally:
        second_app.state.persistence.close()
        second_app.state.run_queue.close()


def test_download_route_never_resolves_or_opens_workspace_paths():
    tree = ast.parse(inspect.getsource(artifact_routes))
    function = next(
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == "download_run_artifact"
    )
    referenced_names = {
        node.id for node in ast.walk(function) if isinstance(node, ast.Name)
    }
    called_attributes = {
        node.attr for node in ast.walk(function) if isinstance(node, ast.Attribute)
    }

    assert {"open", "Path", "FileResponse", "os"}.isdisjoint(referenced_names)
    assert {"open", "read_bytes", "read_text", "resolve"}.isdisjoint(called_attributes)
    assert "prepare" in called_attributes


def test_streaming_response_closes_descriptors_when_transmission_raises(
    client,
    monkeypatch,
):
    run_id = _create_run(client)
    artifact = _record_download(client, run_id, "artifact-disconnect", b"content")
    response = _invoke(client, run_id, artifact.artifact_id)
    plan = response._close_callback.__self__

    async def disconnect(_response, _scope, _receive, _send):
        raise ClientDisconnect()

    monkeypatch.setattr(StreamingResponse, "__call__", disconnect)

    async def receive():
        return {"type": "http.disconnect"}

    async def send(_message):
        return None

    with pytest.raises(ClientDisconnect):
        asyncio.run(
            response(
                {"type": "http", "asgi": {"spec_version": "2.4"}},
                receive,
                send,
            )
        )
    assert plan.closed is True


def test_route_closes_prepared_plan_when_response_construction_raises(
    client,
    monkeypatch,
):
    run_id = _create_run(client)
    artifact = _record_download(client, run_id, "artifact-response", b"content")
    prepared = client.app.state.artifact_download_service.prepare(
        run_id, artifact.artifact_id
    )
    plan = prepared.value

    class PreparedService:
        def prepare(self, _run_id, _artifact_id):
            return prepared

    def fail_response(*_args, **_kwargs):
        raise MemoryError("private response construction failure")

    monkeypatch.setattr(artifact_routes, "_ArtifactStreamingResponse", fail_response)

    with pytest.raises(MemoryError):
        download_run_artifact(run_id, artifact.artifact_id, PreparedService())
    assert plan.closed is True


def test_download_error_projects_only_allowlisted_context(client):
    run_id = _create_run(client)

    class MaliciousService:
        def prepare(self, _run_id, _artifact_id):
            return Result.failure(
                [
                    Issue(
                        code="RUN_ARTIFACT_DOWNLOAD_CONFLICT",
                        message="Artifact content is no longer available as indexed.",
                        path="artifact_id",
                        source="artifact_download_service",
                        technical_message="DATABASE_URL=sqlite:////private/platform.db",
                        context={
                            "reason_code": "ARTIFACT_DOWNLOAD_SOURCE_CHANGED",
                            "path": "/private/workspace/result.tsv",
                            "secret": "DATABASE_PASSWORD=private",
                        },
                    )
                ]
            )

    response = download_run_artifact(
        run_id,
        "artifact-malicious",
        MaliciousService(),
    )

    assert response.status_code == 409
    body = _json(response)
    assert body["issues"][0]["context"] == {
        "reason_code": "ARTIFACT_DOWNLOAD_SOURCE_CHANGED"
    }
    serialized = response.body.decode()
    assert "private" not in serialized
    assert "DATABASE" not in serialized
    assert "technical_message" not in body["issues"][0]
