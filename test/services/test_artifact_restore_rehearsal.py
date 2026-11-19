"""IRENE-26: local single-artifact restore rehearsal on isolated copies."""

from __future__ import annotations

import asyncio
from hashlib import sha256
import importlib
from pathlib import Path
import shutil
import sqlite3
import sys

import httpx

from conftest import (
    TEST_AUTH_CSRF_TOKEN,
    TEST_AUTH_SESSION_TOKEN,
    seed_test_authentication,
)
from encode_pipeline.api.main import create_app
from encode_pipeline.platform.adapters import (
    ExtractedArtifactCandidate,
    WorkflowInputs,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.artifact_extraction import ArtifactExtractionService
from test_artifact_extraction import ArtifactAdapter, _project

_API_DIR = Path(__file__).resolve().parents[1] / "api"
if str(_API_DIR) not in sys.path:
    sys.path.insert(0, str(_API_DIR))
ApiTestClient = importlib.import_module("api_test_client").ApiTestClient


RELATIVE_PATH = "results/browser/sample_peaks.narrowPeak"
SOURCE_CONTENT = b"chr1\t100\t200\tpeak-a\t50\t+\nchr1\t300\t450\tpeak-b\t80\t+\n"


def _close(app) -> None:
    app.state.run_queue.close()
    app.state.persistence.close()


def _app(tmp_path: Path, registry: WorkflowRegistry):
    app = create_app(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        workspace_root=tmp_path / "workspaces",
        project_root=tmp_path / "project",
        registry=registry,
    )
    seed_test_authentication(app)
    run_service = app.state.run_service
    run_id = run_service.create_run("fake", WorkflowInputs(config={})).run_id
    run_service.transition_run(run_id, RunStatus.VALIDATING)
    identity = app.state.build_identity_provider.capture("fake").value
    run_service.complete_preflight(run_id, identity)
    for status in (RunStatus.QUEUED, RunStatus.RUNNING, RunStatus.SUCCEEDED):
        run_service.transition_run(run_id, status)
    source = app.state.workspace_root / run_id / RELATIVE_PATH
    source.parent.mkdir(parents=True, exist_ok=True)
    source.write_bytes(SOURCE_CONTENT)
    return app, run_id, source


def _extract(app, registry: WorkflowRegistry, run_id: str):
    service = ArtifactExtractionService(
        run_service=app.state.run_service,
        registry=registry,
        build_identity_provider=app.state.build_identity_provider,
        workspace_root=app.state.workspace_root,
    )
    result = service.extract(run_id)
    assert result.is_success, result.issues
    assert len(result.value) == 1
    artifact = result.value[0]
    generation = app.state.run_service.get_result_state(run_id).artifact_generation
    assert generation is not None
    return generation, artifact


def _download_url(run_id: str, artifact_id: str) -> str:
    return f"/api/v1/runs/{run_id}/artifacts/{artifact_id}/download"


def _download(client: ApiTestClient, run_id: str, artifact, generation: str):
    return client.get(
        _download_url(run_id, artifact.artifact_id),
        params={"generation": generation, "revision": artifact.revision},
    )


def _prepare_track(client: ApiTestClient, run_id: str, artifact, generation: str):
    response = client.post(
        f"/api/v1/runs/{run_id}/browser-tracks/{artifact.artifact_id}/derivative",
        json={"generation": generation, "revision": artifact.revision},
    )
    assert response.status_code == 200, response.text
    return response.json()["track"]


def _ranges(client: ApiTestClient, track) -> tuple[bytes, bytes]:
    data = client.get(track["content_path"], headers={"Range": "bytes=0-3"})
    index = client.get(track["index_path"], headers={"Range": "bytes=0-3"})
    assert data.status_code == index.status_code == 206
    assert data.headers["content-range"].startswith("bytes 0-3/")
    assert index.headers["content-range"].startswith("bytes 0-3/")
    assert data.content == index.content == b"\x1f\x8b\x08\x04"
    return data.content, index.content


def _backup_database(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with sqlite3.connect(f"file:{source}?mode=ro", uri=True) as original:
        with sqlite3.connect(destination) as backup:
            original.backup(backup)
    with sqlite3.connect(f"file:{destination}?mode=ro", uri=True) as backup:
        assert backup.execute("PRAGMA integrity_check").fetchone() == ("ok",)


def _backup_file(source: Path, destination: Path) -> str:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    digest = sha256(source.read_bytes()).hexdigest()
    assert sha256(destination.read_bytes()).hexdigest() == digest
    return digest


def _anonymous_get(app, url: str) -> httpx.Response:
    async def request() -> httpx.Response:
        transport = httpx.ASGITransport(app=app)
        async with httpx.AsyncClient(
            transport=transport,
            base_url="http://testserver",
        ) as client:
            return await client.get(url, headers={"Range": "bytes=0-3"})

    return asyncio.run(request())


def test_single_artifact_restore_requires_original_layout_and_reindexing(
    tmp_path: Path,
) -> None:
    candidate = ExtractedArtifactCandidate(
        output_type="macs3_peak",
        relative_path=RELATIVE_PATH,
        mime_type="text/tab-separated-values",
        metadata={"catalog_id": "macs3_peak", "scope": "project"},
    )
    adapter = ArtifactAdapter((candidate,))
    registry = WorkflowRegistry([adapter])
    live = tmp_path / "live"
    backup = tmp_path / "backup"
    _project(live / "project")
    adapter.identity_root = live / "project"
    apps = []
    try:
        app, run_id, source = _app(live, registry)
        apps.append(app)
        generation, artifact = _extract(app, registry, run_id)
        client = ApiTestClient(app)

        downloaded = _download(client, run_id, artifact, generation)
        assert downloaded.status_code == 200
        assert downloaded.content == SOURCE_CONTENT
        track = _prepare_track(client, run_id, artifact, generation)
        original_ranges = _ranges(client, track)
        derivative = Path(
            app.state.workspace_root,
            "browser-derivatives",
            run_id,
            artifact.artifact_id,
            artifact.revision,
        )
        backed_up_digests = {
            "source": sha256(SOURCE_CONTENT).hexdigest(),
            "track": sha256((derivative / "track.gz").read_bytes()).hexdigest(),
            "index": sha256((derivative / "track.gz.tbi").read_bytes()).hexdigest(),
        }

        _close(app)
        apps.remove(app)
        database = live / "platform.db"
        _backup_database(database, backup / "platform.db")
        for name, path in (
            ("source", source),
            ("track", derivative / "track.gz"),
            ("index", derivative / "track.gz.tbi"),
            ("provenance", derivative / "provenance.json"),
        ):
            digest = _backup_file(path, backup / name)
            if name in backed_up_digests:
                assert digest == backed_up_digests[name]

        database.unlink()
        source.unlink()
        shutil.copy2(backup / "platform.db", database)
        restored_app = create_app(
            database_url=f"sqlite:///{database}",
            workspace_root=live / "workspaces",
            project_root=live / "project",
            registry=registry,
        )
        apps.append(restored_app)
        restored_app.state.test_auth_tokens = (
            TEST_AUTH_SESSION_TOKEN,
            TEST_AUTH_CSRF_TOKEN,
        )
        restored_client = ApiTestClient(restored_app)
        restored_state = restored_app.state.run_service.get_result_state(run_id)
        restored_artifact = restored_app.state.run_service.get_artifact(
            run_id, artifact.artifact_id
        )
        assert restored_state.artifact_generation == generation
        assert restored_artifact.revision == artifact.revision

        missing = _download(restored_client, run_id, artifact, generation)
        assert missing.status_code == 409
        assert missing.json()["issues"][0]["context"] == {
            "reason_code": "ARTIFACT_DOWNLOAD_SOURCE_MISSING_OR_UNSAFE"
        }

        shutil.copy2(backup / "source", source)
        assert sha256(source.read_bytes()).hexdigest() == backed_up_digests["source"]
        copied = _download(restored_client, run_id, artifact, generation)
        assert copied.status_code == 409
        assert copied.json()["issues"][0]["context"] == {
            "reason_code": "ARTIFACT_DOWNLOAD_REVISION_MISMATCH"
        }

        restored_generation, indexed = _extract(restored_app, registry, run_id)
        assert indexed.run_id == artifact.run_id
        assert indexed.artifact_id == artifact.artifact_id
        assert restored_generation != generation
        assert indexed.revision != artifact.revision
        stale = _download(restored_client, run_id, artifact, generation)
        assert stale.status_code == 409
        assert stale.json()["issues"][0]["context"] == {
            "reason_code": "ARTIFACT_DOWNLOAD_GENERATION_CHANGED"
        }

        recovered = _download(restored_client, run_id, indexed, restored_generation)
        assert recovered.status_code == 200
        assert sha256(recovered.content).hexdigest() == backed_up_digests["source"]
        restored_track = _prepare_track(
            restored_client, run_id, indexed, restored_generation
        )
        assert restored_track["revision"] == indexed.revision
        assert _ranges(restored_client, restored_track) == original_ranges
        anonymous = _anonymous_get(restored_app, restored_track["content_path"])
        assert anonymous.status_code == 401
        assert str(live) not in anonymous.text
    finally:
        for running_app in reversed(apps):
            _close(running_app)
