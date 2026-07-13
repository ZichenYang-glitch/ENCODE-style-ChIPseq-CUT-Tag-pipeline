"""Tests for the read-only persisted artifact API."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import replace
from datetime import datetime, timezone

import pytest

from encode_pipeline.api.main import create_app
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.runs import RunArtifactRef
from api_test_client import ApiTestClient


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


@pytest.fixture
def client(tmp_path) -> Iterator[ApiTestClient]:
    app = create_app(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        workspace_root=tmp_path / "workspaces",
    )
    with ApiTestClient(app) as test_client:
        yield test_client


def _create_run(client: ApiTestClient) -> str:
    return client.app.state.run_service.create_run(
        WORKFLOW_ID,
        WorkflowInputs(config={}),
    ).run_id


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
    client.app.state.run_service.record_artifact(artifact.run_id, artifact)


def test_list_existing_run_without_artifacts_returns_empty_page(client):
    run_id = _create_run(client)

    response = client.get(f"/api/v1/runs/{run_id}/artifacts")

    assert response.status_code == 200
    assert response.json() == {
        "ok": True,
        "run_id": run_id,
        "artifacts": [],
        "issues": [],
    }


def test_list_is_stably_sorted_and_uses_bounded_artifact_cursor(client):
    run_id = _create_run(client)
    for artifact_id in ("artifact-z", "artifact-a", "artifact-m"):
        _record(client, _artifact(run_id, artifact_id))

    first = client.get(f"/api/v1/runs/{run_id}/artifacts?limit=2")
    second = client.get(
        f"/api/v1/runs/{run_id}/artifacts",
        params={"after": "artifact-m", "limit": 2},
    )

    assert first.status_code == 200
    assert [item["artifact_id"] for item in first.json()["artifacts"]] == [
        "artifact-a",
        "artifact-m",
    ]
    assert first.json()["next_cursor"] == "artifact-m"
    assert [item["artifact_id"] for item in second.json()["artifacts"]] == [
        "artifact-z"
    ]
    assert "next_cursor" not in second.json()


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

    response = client.get(f"/api/v1/runs/{run_id}/artifacts/{artifact.artifact_id}")

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
        "relative_path": "results/artifact-summary.txt",
        "output_type": "summary",
        "size_bytes": 7,
        "metadata": {"catalog_id": "summary", "scope": "project"},
    }
    assert "private" not in response.text


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

    response = client.get(f"/api/v1/runs/{run_id}/artifacts/{artifact.artifact_id}")

    assert response.status_code == 200
    assert response.json()["artifact"]["metadata"]["target"] == target


def test_unknown_run_returns_stable_404(client):
    response = client.get("/api/v1/runs/run-missing/artifacts")

    assert response.status_code == 404
    assert response.json()["issues"][0]["code"] == "RUN_NOT_FOUND"


def test_unknown_and_cross_run_artifacts_share_the_same_404(client):
    first_run = _create_run(client)
    second_run = _create_run(client)
    other = _artifact(second_run, "artifact-other")
    _record(client, other)

    missing = client.get(f"/api/v1/runs/{first_run}/artifacts/artifact-missing")
    cross_run = client.get(f"/api/v1/runs/{first_run}/artifacts/{other.artifact_id}")

    assert missing.status_code == cross_run.status_code == 404
    assert missing.json()["issues"][0]["code"] == "RUN_ARTIFACT_NOT_FOUND"
    assert cross_run.json()["issues"][0]["code"] == "RUN_ARTIFACT_NOT_FOUND"


def test_unknown_and_cross_run_cursors_share_the_same_400(client):
    first_run = _create_run(client)
    second_run = _create_run(client)
    other = _artifact(second_run, "artifact-other")
    _record(client, other)

    missing = client.get(
        f"/api/v1/runs/{first_run}/artifacts",
        params={"after": "artifact-missing"},
    )
    cross_run = client.get(
        f"/api/v1/runs/{first_run}/artifacts",
        params={"after": other.artifact_id},
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

    detail = client.get(f"/api/v1/runs/{run_id}/artifacts/{unsafe.artifact_id}")
    assert detail.status_code == 500
    assert detail.json()["artifact"] is None
    assert detail.json()["issues"][0]["code"] == "RUN_ARTIFACT_DATA_INVALID"
    assert private_value not in detail.text


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
            f"/api/v1/runs/{run_id}/artifacts/{artifact.artifact_id}"
        )

    assert response.status_code == 200
    assert response.json()["artifact"]["artifact_id"] == "artifact-persisted"
