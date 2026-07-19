"""Tests for the read-only persisted artifact API."""

from __future__ import annotations

import asyncio
from collections.abc import Iterator
from dataclasses import replace
from datetime import datetime, timezone

import httpx
import pytest

from encode_pipeline.api.main import create_app
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.result_generations import (
    ArtifactCursor,
    encode_artifact_cursor,
)
from encode_pipeline.platform.runs import RunArtifactRef, RunStatus
from api_test_client import ApiTestClient


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
ARTIFACT_REVISION = f"artifactrev-{'a' * 64}"
ARTIFACT_REVISION_B = f"artifactrev-{'b' * 64}"


@pytest.fixture
def client(tmp_path) -> Iterator[ApiTestClient]:
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


def _artifact(
    run_id: str,
    artifact_id: str,
    *,
    metadata: dict[str, object] | None = None,
) -> RunArtifactRef:
    relative_path = f"results/{artifact_id}.txt"
    return RunArtifactRef(
        artifact_id=artifact_id,
        run_id=run_id,
        artifact_type="file",
        name=f"{artifact_id}.txt",
        uri=f"run://runs/{run_id}/artifacts/{artifact_id}",
        mime_type="text/plain",
        produced_at=datetime.now(timezone.utc),
        revision=ARTIFACT_REVISION,
        metadata=metadata
        or {
            "relative_path": relative_path,
            "output_type": "summary",
            "size_bytes": 7,
            "catalog_id": "summary",
            "scope": "project",
        },
    )


def _record(client: ApiTestClient, artifact: RunArtifactRef) -> None:
    service = client.app.state.run_service
    service.replace_artifacts(
        artifact.run_id,
        (*service.list_artifacts(artifact.run_id), artifact),
    )


def _generation(client: ApiTestClient, run_id: str) -> str:
    generation = client.app.state.run_service.get_result_state(
        run_id
    ).artifact_generation
    assert generation is not None
    return generation


def test_list_existing_run_without_artifacts_returns_empty_page(client):
    run_id = _create_run(client)

    response = client.get(f"/api/v1/runs/{run_id}/artifacts")

    assert response.status_code == 200
    assert response.json() == {
        "ok": True,
        "run_id": run_id,
        "artifact_generation": None,
        "artifacts": [],
        "issues": [],
    }


def test_list_is_stably_sorted_and_uses_bounded_artifact_cursor(client):
    run_id = _create_run(client)
    for artifact_id in ("artifact-z", "artifact-a", "artifact-m"):
        _record(client, _artifact(run_id, artifact_id))

    first = client.get(f"/api/v1/runs/{run_id}/artifacts?limit=2")
    generation = first.json()["artifact_generation"]
    cursor = first.json()["next_cursor"]
    second = client.get(
        f"/api/v1/runs/{run_id}/artifacts",
        params={"after": cursor, "generation": generation, "limit": 2},
    )

    assert first.status_code == 200
    assert [item["artifact_id"] for item in first.json()["artifacts"]] == [
        "artifact-a",
        "artifact-m",
    ]
    assert cursor.startswith("artifactcur_")
    assert [item["artifact_id"] for item in second.json()["artifacts"]] == [
        "artifact-z"
    ]
    assert "next_cursor" not in second.json()


def test_list_exposes_generation_and_uses_an_opaque_generation_bound_cursor(client):
    run_id = _create_run(client)
    for artifact_id in ("artifact-a", "artifact-b"):
        _record(client, _artifact(run_id, artifact_id))

    first = client.get(f"/api/v1/runs/{run_id}/artifacts?limit=1")
    body = first.json()
    generation = body["artifact_generation"]
    cursor = body["next_cursor"]

    assert first.status_code == 200
    assert generation.startswith("artifactgen-")
    assert cursor.startswith("artifactcur_")
    assert "artifact-a" not in cursor
    second = client.get(
        f"/api/v1/runs/{run_id}/artifacts",
        params={"after": cursor, "generation": generation, "limit": 1},
    )
    assert second.status_code == 200
    assert second.json()["artifact_generation"] == generation
    assert [item["artifact_id"] for item in second.json()["artifacts"]] == [
        "artifact-b"
    ]


def test_stale_artifact_page_returns_409_after_same_id_count_replacement(client):
    run_id = _create_run(client)
    original = (_artifact(run_id, "artifact-a"), _artifact(run_id, "artifact-b"))
    client.app.state.run_service.replace_artifacts(run_id, original)
    first = client.get(f"/api/v1/runs/{run_id}/artifacts?limit=1").json()

    replacement = tuple(
        replace(artifact, revision=ARTIFACT_REVISION_B) for artifact in original
    )
    client.app.state.run_service.replace_artifacts(run_id, replacement)
    response = client.get(
        f"/api/v1/runs/{run_id}/artifacts",
        params={
            "after": first["next_cursor"],
            "generation": first["artifact_generation"],
            "limit": 1,
        },
    )

    assert response.status_code == 409
    assert response.json()["artifact_generation"] is None
    assert response.json()["artifacts"] == []
    assert response.json()["issues"][0]["code"] == "RUN_ARTIFACT_GENERATION_CHANGED"


def test_stale_deleted_detail_returns_generation_conflict_before_not_found(client):
    run_id = _create_run(client)
    artifact = _artifact(run_id, "artifact-a")
    client.app.state.run_service.replace_artifacts(run_id, (artifact,))
    generation_a = client.app.state.run_service.get_result_state(
        run_id
    ).artifact_generation
    assert generation_a is not None
    client.app.state.run_service.replace_artifacts(run_id, ())
    generation_b = client.app.state.run_service.get_result_state(
        run_id
    ).artifact_generation
    assert generation_b is not None

    stale = client.get(
        f"/api/v1/runs/{run_id}/artifacts/{artifact.artifact_id}",
        params={"generation": generation_a},
    )
    current = client.get(
        f"/api/v1/runs/{run_id}/artifacts/{artifact.artifact_id}",
        params={"generation": generation_b},
    )

    assert stale.status_code == 409
    assert stale.json()["issues"][0]["code"] == "RUN_ARTIFACT_GENERATION_CHANGED"
    assert current.status_code == 404


def test_detail_returns_strict_public_projection(client):
    run_id = _create_run(client)
    artifact = _artifact(
        run_id,
        "artifact-summary",
        metadata={
            "relative_path": "results/artifact-summary.txt",
            "output_type": "summary",
            "size_bytes": 7,
            "catalog_id": "summary",
            "scope": "project",
            "private_extension": {"path": "/private/workspace"},
        },
    )
    _record(client, artifact)

    response = client.get(
        f"/api/v1/runs/{run_id}/artifacts/{artifact.artifact_id}",
        params={"generation": _generation(client, run_id)},
    )

    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert data["run_id"] == run_id
    assert data["artifact"] == {
        "artifact_id": "artifact-summary",
        "run_id": run_id,
        "artifact_type": "file",
        "name": "artifact-summary.txt",
        "uri": f"run://runs/{run_id}/artifacts/artifact-summary",
        "mime_type": "text/plain",
        "produced_at": artifact.produced_at.isoformat().replace("+00:00", "Z"),
        "revision": ARTIFACT_REVISION,
        "relative_path": "results/artifact-summary.txt",
        "output_type": "summary",
        "size_bytes": 7,
        "metadata": {"catalog_id": "summary", "scope": "project"},
    }
    assert "private" not in response.text


def test_detail_requires_the_generation_seen_by_the_caller(client):
    run_id = _create_run(client)
    artifact = _artifact(run_id, "artifact-summary")
    _record(client, artifact)

    response = client.get(f"/api/v1/runs/{run_id}/artifacts/{artifact.artifact_id}")

    assert response.status_code == 400
    assert response.json()["artifact_generation"] is None
    assert response.json()["artifact"] is None
    assert response.json()["issues"][0]["code"] == "API_REQUEST_INVALID"


@pytest.mark.parametrize("target", ["RNA polymerase II", "Pol II/III"])
def test_legal_scientific_metadata_text_is_returned(client, target):
    run_id = _create_run(client)
    artifact = _artifact(
        run_id,
        "artifact-polymerase",
        metadata={
            "relative_path": "results/artifact-polymerase.txt",
            "output_type": "summary",
            "size_bytes": 7,
            "target": target,
        },
    )
    _record(client, artifact)

    response = client.get(
        f"/api/v1/runs/{run_id}/artifacts/{artifact.artifact_id}",
        params={"generation": _generation(client, run_id)},
    )

    assert response.status_code == 200
    assert response.json()["artifact"]["metadata"]["target"] == target


def test_unknown_run_returns_stable_404(client):
    response = client.get("/api/v1/runs/run-missing/artifacts")

    assert response.status_code == 404
    assert response.json()["issues"][0]["code"] == "RUN_NOT_FOUND"


def test_unknown_and_cross_run_artifacts_share_the_same_404(client):
    first_run = _create_run(client)
    second_run = _create_run(client)
    _record(client, _artifact(first_run, "artifact-present"))
    other = _artifact(second_run, "artifact-other")
    _record(client, other)

    generation = _generation(client, first_run)
    missing = client.get(
        f"/api/v1/runs/{first_run}/artifacts/artifact-missing",
        params={"generation": generation},
    )
    cross_run = client.get(
        f"/api/v1/runs/{first_run}/artifacts/{other.artifact_id}",
        params={"generation": generation},
    )

    assert missing.status_code == cross_run.status_code == 404
    assert missing.json()["issues"][0]["code"] == "RUN_ARTIFACT_NOT_FOUND"
    assert cross_run.json()["issues"][0]["code"] == "RUN_ARTIFACT_NOT_FOUND"


def test_unknown_and_cross_run_cursors_share_the_same_400(client):
    first_run = _create_run(client)
    second_run = _create_run(client)
    other = _artifact(second_run, "artifact-other")
    _record(client, other)
    client.app.state.run_service.replace_artifacts(
        first_run, (_artifact(first_run, "artifact-present"),)
    )
    first_generation = client.app.state.run_service.get_result_state(
        first_run
    ).artifact_generation
    second_generation = client.app.state.run_service.get_result_state(
        second_run
    ).artifact_generation
    assert first_generation is not None and second_generation is not None
    missing_cursor = encode_artifact_cursor(
        ArtifactCursor(
            run_id=first_run,
            artifact_generation=first_generation,
            after_artifact_id="artifact-missing",
        )
    )
    cross_run_cursor = encode_artifact_cursor(
        ArtifactCursor(
            run_id=second_run,
            artifact_generation=second_generation,
            after_artifact_id=other.artifact_id,
        )
    )

    missing = client.get(
        f"/api/v1/runs/{first_run}/artifacts",
        params={"after": missing_cursor, "generation": first_generation},
    )
    cross_run = client.get(
        f"/api/v1/runs/{first_run}/artifacts",
        params={"after": cross_run_cursor, "generation": first_generation},
    )

    assert missing.status_code == cross_run.status_code == 400
    assert missing.json()["issues"][0]["code"] == ("RUN_ARTIFACT_CURSOR_NOT_FOUND")
    assert cross_run.json()["issues"][0]["code"] == ("RUN_ARTIFACT_CURSOR_NOT_FOUND")


@pytest.mark.parametrize("limit", [0, 101])
def test_list_rejects_out_of_bounds_limits(client, limit):
    run_id = _create_run(client)

    response = client.get(
        f"/api/v1/runs/{run_id}/artifacts",
        params={"limit": limit},
    )

    assert response.status_code == 400
    assert response.json() == {
        "ok": False,
        "run_id": run_id,
        "artifact_generation": None,
        "artifacts": [],
        "issues": [
            {
                "code": "API_REQUEST_INVALID",
                "message": "Artifact query parameters are invalid.",
                "severity": "error",
                "path": "limit",
                "source": "api",
                "hint": "Use an integer limit between 1 and 100.",
                "context": {},
            }
        ],
    }


@pytest.mark.parametrize("limit", [1, 100])
def test_list_accepts_bounded_limit_endpoints(client, limit):
    run_id = _create_run(client)

    response = client.get(
        f"/api/v1/runs/{run_id}/artifacts",
        params={"limit": limit},
    )

    assert response.status_code == 200
    assert response.json()["artifacts"] == []


@pytest.mark.parametrize(
    "mutation,private_value",
    [
        ("relative_path", "/private/workspace/output.txt"),
        ("uri", "file:///private/workspace/output.txt"),
        ("name", "private/output.txt"),
        ("target", "/private/target"),
        ("windows_target", "C:" + chr(92) + "private" + chr(92) + "target"),
        ("unc_target", chr(92) * 2 + "private" + chr(92) + "target"),
        ("environment_target", "WORKSPACE=/private/target"),
        ("prefixed_environment", "prefix WORKSPACE=/private/target"),
        ("json_path", '"path":"/private/target"'),
        ("colon_path", "see:/private/target"),
    ],
)
def test_unsafe_persisted_artifact_fails_closed_without_disclosure(
    client,
    mutation,
    private_value,
):
    run_id = _create_run(client)
    safe = _artifact(run_id, "artifact-a")
    unsafe = _artifact(run_id, "artifact-b")
    if mutation == "relative_path":
        unsafe = replace(
            unsafe,
            metadata={**unsafe.metadata, "relative_path": private_value},
        )
    elif mutation in {
        "target",
        "windows_target",
        "unc_target",
        "environment_target",
        "prefixed_environment",
        "json_path",
        "colon_path",
    }:
        unsafe = replace(
            unsafe,
            metadata={**unsafe.metadata, "target": private_value},
        )
    else:
        unsafe = replace(unsafe, **{mutation: private_value})
    _record(client, safe)
    _record(client, unsafe)

    response = client.get(f"/api/v1/runs/{run_id}/artifacts")

    assert response.status_code == 500
    assert response.json()["artifacts"] == []
    assert response.json()["issues"][0]["code"] == "RUN_ARTIFACT_DATA_INVALID"
    assert private_value not in response.text

    detail = client.get(
        f"/api/v1/runs/{run_id}/artifacts/{unsafe.artifact_id}",
        params={"generation": _generation(client, run_id)},
    )
    assert detail.status_code == 500
    assert detail.json()["artifact"] is None
    assert detail.json()["issues"][0]["code"] == "RUN_ARTIFACT_DATA_INVALID"
    assert private_value not in detail.text


@pytest.mark.parametrize(
    ("operation", "service_method", "expected"),
    [
        (
            "list",
            "list_artifacts_page",
            {
                "artifact_generation": None,
                "artifacts": [],
                "issues_path": "artifacts",
            },
        ),
        (
            "detail",
            "get_artifact_at_generation",
            {
                "artifact_generation": None,
                "artifact": None,
                "issues_path": "artifact_id",
            },
        ),
    ],
)
def test_unexpected_artifact_repository_failure_uses_declared_safe_envelope(
    client,
    monkeypatch,
    operation,
    service_method,
    expected,
):
    run_id = _create_run(client)
    persisted = _artifact(run_id, "artifact-a")
    _record(client, persisted)
    generation = _generation(client, run_id)
    private_text = "sqlite:///private/artifacts.db SECRET_TOKEN=value"

    def fail(*_args, **_kwargs):
        raise RuntimeError(private_text)

    monkeypatch.setattr(client.app.state.run_service, service_method, fail)
    if operation == "list":
        url = f"/api/v1/runs/{run_id}/artifacts"
        params = {"generation": generation}
    else:
        url = f"/api/v1/runs/{run_id}/artifacts/{persisted.artifact_id}"
        params = {"generation": generation}

    async def request() -> httpx.Response:
        transport = httpx.ASGITransport(
            app=client.app,
            raise_app_exceptions=False,
        )
        async with httpx.AsyncClient(
            transport=transport,
            base_url="http://testserver",
        ) as async_client:
            return await async_client.get(url, params=params)

    response = asyncio.run(request())
    body = response.json()

    assert response.status_code == 500
    assert body["ok"] is False
    assert body["run_id"] == run_id
    assert body["artifact_generation"] is expected["artifact_generation"]
    assert body.get("artifacts") == expected.get("artifacts")
    assert body.get("artifact") == expected.get("artifact")
    assert body["issues"] == [
        {
            "code": "INTERNAL_SERVER_ERROR",
            "message": "Run artifacts are temporarily unavailable.",
            "severity": "error",
            "path": expected["issues_path"],
            "source": "runtime",
            "context": {},
        }
    ]
    assert private_text not in response.text
    assert "technical_message" not in response.text


def test_artifacts_survive_sqlite_reopen(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    workspace_root = tmp_path / "workspaces"
    first_app = create_app(database_url=database_url, workspace_root=workspace_root)
    with ApiTestClient(first_app) as first_client:
        run_id = _create_run(first_client)
        artifact = _artifact(run_id, "artifact-persisted")
        _record(first_client, artifact)

    second_app = create_app(database_url=database_url, workspace_root=workspace_root)
    with ApiTestClient(second_app) as second_client:
        response = second_client.get(
            f"/api/v1/runs/{run_id}/artifacts/{artifact.artifact_id}",
            params={"generation": _generation(second_client, run_id)},
        )

    assert response.status_code == 200
    assert response.json()["artifact"]["artifact_id"] == "artifact-persisted"
