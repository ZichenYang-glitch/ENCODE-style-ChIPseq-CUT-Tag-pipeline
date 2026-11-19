"""IRENE-33: resume a scratch-to-durable copy without changing the saved link."""

from __future__ import annotations

import asyncio
import errno
import importlib
import os
from hashlib import sha256
from pathlib import Path
import shutil
import stat
import sys
import threading

import httpx

from conftest import seed_test_authentication
from encode_pipeline.api.main import create_app
from encode_pipeline.platform.adapters import (
    ExtractedArtifactCandidate,
    WorkflowInputs,
)
from encode_pipeline.platform.result_generations import (
    ARTIFACT_PATH_IDENTITY_METADATA_KEY,
    build_artifact_descriptor_revision,
    build_artifact_path_identity,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.artifact_disk_relocation import (
    MAX_RELOCATION_BYTES,
    PUBLIC_BYTES_MISMATCH,
    PUBLIC_PERMISSION_MISMATCH,
    PUBLIC_UNAVAILABLE,
    ArtifactDiskRelocation,
    ArtifactDiskRelocationError,
    DownloadRelocation,
    RelocationJournal,
    lookup_switched_download,
)
from encode_pipeline.services import artifact_downloads as artifact_downloads_module
from encode_pipeline.services.artifact_extraction import ArtifactExtractionService
from test_artifact_extraction import ArtifactAdapter, _project

_API_DIR = Path(__file__).resolve().parents[1] / "api"
if str(_API_DIR) not in sys.path:
    sys.path.insert(0, str(_API_DIR))
ApiTestClient = importlib.import_module("api_test_client").ApiTestClient


RELATIVE_PATH = "results/browser/sample_peaks.narrowPeak"
SOURCE_CONTENT = b"chr1\t100\t200\tpeak-a\t50\t+\nchr1\t300\t450\tpeak-b\t80\t+\n"
REPLACEMENT_CONTENT = b"chr1\t100\t200\tpeak-Z\t99\t+\nchr1\t300\t450\tpeak-Y\t99\t+\n"
SOURCE_MODE = 0o640
OUTPUT_TYPE = "macs3_peak"


class _InterruptCommit(ArtifactDiskRelocation):
    """Leave a partial durable file and stop before the copy stage commits."""

    def _commit_payload(self, destination: Path, content: bytes, mode: int) -> None:
        del content, mode
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.with_name(destination.name + ".partial").write_bytes(b"PARTIAL")
        raise OSError(errno.EIO, "interrupted")


def _close(app) -> None:
    app.state.run_queue.close()
    app.state.persistence.close()


def _reason(response: httpx.Response) -> str | None:
    if response.status_code == 200:
        return None
    context = response.json()["issues"][0].get("context") or {}
    return context.get("reason_code")


def _assert_private(response: httpx.Response, *paths: Path) -> None:
    text = response.text
    for path in paths:
        assert str(path) not in text


def _anonymous(app, url: str, params: dict[str, str]) -> httpx.Response:
    async def request() -> httpx.Response:
        transport = httpx.ASGITransport(app=app)
        async with httpx.AsyncClient(
            transport=transport,
            base_url="http://testserver",
        ) as client:
            return await client.get(url, params=params)

    return asyncio.run(request())


def _current_binding(path: Path) -> tuple[str, str]:
    infos: list[os.stat_result] = []
    current = Path("/")
    for component in path.parts[1:]:
        current /= component
        infos.append(os.stat(current, follow_symlinks=False))
    file_info = infos[-1]
    file_identity = (
        file_info.st_dev,
        file_info.st_ino,
        file_info.st_mode,
        file_info.st_nlink,
        file_info.st_uid,
        file_info.st_gid,
        file_info.st_size,
        file_info.st_mtime_ns,
        file_info.st_ctime_ns,
    )
    revision = build_artifact_descriptor_revision(
        output_type=OUTPUT_TYPE,
        relative_path=RELATIVE_PATH,
        descriptor_identity=":".join(str(value) for value in file_identity),
    )
    path_identity = build_artifact_path_identity(
        parent_identities=tuple(
            (info.st_dev, info.st_ino, info.st_mode) for info in infos[:-1]
        ),
        file_identity=file_identity,
    )
    return revision, path_identity


def test_scratch_change_before_copy_does_not_publish(tmp_path: Path) -> None:
    source_root = tmp_path / "scratch"
    durable_root = tmp_path / "durable"
    target = source_root / "run-1" / RELATIVE_PATH
    target.parent.mkdir(parents=True)
    target.write_bytes(SOURCE_CONTENT)
    os.chmod(target, SOURCE_MODE)
    revision, path_identity = _current_binding(target)
    relocation = ArtifactDiskRelocation(
        source_root=source_root,
        durable_root=durable_root,
    )
    relocation.capture(
        run_id="run-1",
        artifact_id="artifact-1",
        generation="artifactgen-" + "ab" * 32,
        revision=revision,
        relative_path=RELATIVE_PATH,
        output_type=OUTPUT_TYPE,
        path_identity=path_identity,
        expected_size=len(SOURCE_CONTENT),
        content=SOURCE_CONTENT,
        mode=SOURCE_MODE,
    )
    target.write_bytes(b"changed\n")
    os.chmod(target, SOURCE_MODE)

    try:
        relocation.resume(run_id="run-1", artifact_id="artifact-1")
    except ArtifactDiskRelocationError as exc:
        assert exc.reason_code == "relocation_source_changed"
    else:
        raise AssertionError("changed scratch bytes were copied")

    assert relocation.load("run-1", "artifact-1").stage == "copying"
    payload = relocation.payload_file("run-1", "artifact-1")
    assert not payload.exists()
    payload.parent.mkdir(parents=True, exist_ok=True)
    payload.write_bytes(b"PARTIAL")
    decision = lookup_switched_download(
        durable_root,
        run_id="run-1",
        artifact_id="artifact-1",
        generation="artifactgen-" + "ab" * 32,
        revision=revision,
        relative_path=RELATIVE_PATH,
        size_bytes=len(SOURCE_CONTENT),
    )
    assert decision.kind == "absent"


def test_equal_length_edit_before_capture_keeps_revision_conflict(
    tmp_path: Path,
) -> None:
    assert len(REPLACEMENT_CONTENT) == len(SOURCE_CONTENT)
    assert REPLACEMENT_CONTENT != SOURCE_CONTENT
    scratch = tmp_path / "scratch"
    durable = tmp_path / "durable"
    project = scratch / "project"
    _project(project)
    candidate = ExtractedArtifactCandidate(
        output_type=OUTPUT_TYPE,
        relative_path=RELATIVE_PATH,
        mime_type="text/tab-separated-values",
        metadata={"catalog_id": "macs3_peak", "scope": "project"},
    )
    adapter = ArtifactAdapter((candidate,))
    adapter.identity_root = project
    registry = WorkflowRegistry([adapter])
    app = create_app(
        database_url=f"sqlite:///{scratch / 'platform.db'}",
        workspace_root=scratch / "workspaces",
        project_root=project,
        registry=registry,
        relocation_root=durable,
    )
    try:
        seed_test_authentication(app)
        run_service = app.state.run_service
        run_id = run_service.create_run("fake", WorkflowInputs(config={})).run_id
        run_service.transition_run(run_id, RunStatus.VALIDATING)
        identity = app.state.build_identity_provider.capture("fake").value
        run_service.complete_preflight(run_id, identity)
        for status in (RunStatus.QUEUED, RunStatus.RUNNING, RunStatus.SUCCEEDED):
            run_service.transition_run(run_id, status)
        source = app.state.workspace_root / run_id / RELATIVE_PATH
        source.parent.mkdir(parents=True)
        source.write_bytes(SOURCE_CONTENT)
        os.chmod(source, SOURCE_MODE)
        extraction = ArtifactExtractionService(
            run_service=run_service,
            registry=registry,
            build_identity_provider=app.state.build_identity_provider,
            workspace_root=app.state.workspace_root,
        )
        extracted = extraction.extract(run_id)
        assert extracted.is_success, extracted.issues
        artifact = extracted.value[0]
        generation = run_service.get_result_state(run_id).artifact_generation
        assert generation is not None
        saved_url = f"/api/v1/runs/{run_id}/artifacts/{artifact.artifact_id}/download"
        saved_params = {"generation": generation, "revision": artifact.revision}
        client = ApiTestClient(app)
        relocation = ArtifactDiskRelocation(
            source_root=app.state.workspace_root,
            durable_root=durable,
        )
        original = client.get(saved_url, params=saved_params)
        assert original.status_code == 200
        assert original.content == SOURCE_CONTENT

        source.write_bytes(REPLACEMENT_CONTENT)
        os.chmod(source, SOURCE_MODE)
        assert len(source.read_bytes()) == len(SOURCE_CONTENT)
        try:
            relocation.capture(
                run_id=run_id,
                artifact_id=artifact.artifact_id,
                generation=generation,
                revision=artifact.revision,
                relative_path=RELATIVE_PATH,
                output_type=str(artifact.metadata["output_type"]),
                path_identity=str(
                    artifact.metadata[ARTIFACT_PATH_IDENTITY_METADATA_KEY]
                ),
                expected_size=int(artifact.metadata["size_bytes"]),
                content=REPLACEMENT_CONTENT,
                mode=SOURCE_MODE,
            )
        except ArtifactDiskRelocationError as exc:
            assert exc.reason_code == "relocation_indexed_revision_mismatch"
        else:
            raise AssertionError("equal-length scratch edit was sealed")

        try:
            app.state.artifact_download_service.capture_relocation(
                relocation,
                run_id=run_id,
                artifact_id=artifact.artifact_id,
                expected_generation=generation,
                expected_revision=artifact.revision,
            )
        except ArtifactDiskRelocationError as exc:
            assert exc.reason_code == "relocation_indexed_revision_mismatch"
        else:
            raise AssertionError("indexed capture published the edited scratch file")

        journal = durable / "journals" / run_id / f"{artifact.artifact_id}.json"
        assert not journal.exists()
        decision = lookup_switched_download(
            durable,
            run_id=run_id,
            artifact_id=artifact.artifact_id,
            generation=generation,
            revision=artifact.revision,
            relative_path=RELATIVE_PATH,
            size_bytes=len(SOURCE_CONTENT),
        )
        assert decision.kind == "absent"
        conflict = client.get(saved_url, params=saved_params)
        assert conflict.status_code == 409
        assert _reason(conflict) == "ARTIFACT_DOWNLOAD_REVISION_MISMATCH"
        assert REPLACEMENT_CONTENT not in conflict.content
        _assert_private(conflict, scratch, durable)
        assert source.read_bytes() == REPLACEMENT_CONTENT
    finally:
        _close(app)


def test_saved_link_reads_durable_bytes_after_resumable_switch(tmp_path: Path) -> None:
    scratch = tmp_path / "scratch"
    durable = tmp_path / "durable"
    project = scratch / "project"
    _project(project)
    candidate = ExtractedArtifactCandidate(
        output_type="macs3_peak",
        relative_path=RELATIVE_PATH,
        mime_type="text/tab-separated-values",
        metadata={"catalog_id": "macs3_peak", "scope": "project"},
    )
    adapter = ArtifactAdapter((candidate,))
    adapter.identity_root = project
    registry = WorkflowRegistry([adapter])
    app = create_app(
        database_url=f"sqlite:///{scratch / 'platform.db'}",
        workspace_root=scratch / "workspaces",
        project_root=project,
        registry=registry,
        relocation_root=durable,
    )
    try:
        seed_test_authentication(app)
        run_service = app.state.run_service
        run_id = run_service.create_run("fake", WorkflowInputs(config={})).run_id
        run_service.transition_run(run_id, RunStatus.VALIDATING)
        identity = app.state.build_identity_provider.capture("fake").value
        run_service.complete_preflight(run_id, identity)
        for status in (RunStatus.QUEUED, RunStatus.RUNNING, RunStatus.SUCCEEDED):
            run_service.transition_run(run_id, status)
        source = app.state.workspace_root / run_id / RELATIVE_PATH
        source.parent.mkdir(parents=True)
        source.write_bytes(SOURCE_CONTENT)
        os.chmod(source, SOURCE_MODE)
        extraction = ArtifactExtractionService(
            run_service=run_service,
            registry=registry,
            build_identity_provider=app.state.build_identity_provider,
            workspace_root=app.state.workspace_root,
        )
        extracted = extraction.extract(run_id)
        assert extracted.is_success, extracted.issues
        artifact = extracted.value[0]
        generation = run_service.get_result_state(run_id).artifact_generation
        assert generation is not None
        saved_url = f"/api/v1/runs/{run_id}/artifacts/{artifact.artifact_id}/download"
        saved_params = {"generation": generation, "revision": artifact.revision}
        client = ApiTestClient(app)
        relocation = ArtifactDiskRelocation(
            source_root=app.state.workspace_root,
            durable_root=durable,
        )

        def read() -> httpx.Response:
            return client.get(saved_url, params=saved_params)

        original = read()
        assert original.status_code == 200
        assert original.content == SOURCE_CONTENT
        anonymous = _anonymous(app, saved_url, saved_params)
        assert anonymous.status_code == 401
        _assert_private(anonymous, scratch, durable)

        sealed = app.state.artifact_download_service.capture_relocation(
            relocation,
            run_id=run_id,
            artifact_id=artifact.artifact_id,
            expected_generation=generation,
            expected_revision=artifact.revision,
        )
        assert sealed.stage == "captured"
        assert sealed.mode == SOURCE_MODE
        assert sealed.content_sha256 == sha256(SOURCE_CONTENT).hexdigest()
        assert sealed.size_bytes == len(SOURCE_CONTENT)
        journal_text = (
            durable / "journals" / run_id / f"{artifact.artifact_id}.json"
        ).read_text(encoding="utf-8")
        assert str(scratch) not in journal_text
        assert str(durable) not in journal_text

        interrupted = _InterruptCommit(
            source_root=app.state.workspace_root,
            durable_root=durable,
        )
        try:
            interrupted.resume(run_id=run_id, artifact_id=artifact.artifact_id)
        except OSError as exc:
            assert exc.errno == errno.EIO
        else:
            raise AssertionError("copy committed through an interruption")
        assert relocation.load(run_id, artifact.artifact_id).stage == "copying"
        partial = relocation.payload_file(run_id, artifact.artifact_id).with_name(
            "sample_peaks.narrowPeak.partial"
        )
        assert partial.read_bytes() == b"PARTIAL"

        during_copy = read()
        assert during_copy.status_code == 200
        assert during_copy.content == SOURCE_CONTENT
        assert b"PARTIAL" not in during_copy.content

        copied = relocation.resume(run_id=run_id, artifact_id=artifact.artifact_id)
        assert copied.stage == "copied"
        payload = relocation.payload_file(run_id, artifact.artifact_id)
        assert payload.read_bytes() == SOURCE_CONTENT
        assert stat.S_IMODE(payload.stat().st_mode) == SOURCE_MODE
        payload.write_bytes(b"tampered-durable")
        os.chmod(payload, SOURCE_MODE)
        assert read().content == SOURCE_CONTENT
        try:
            relocation.resume(run_id=run_id, artifact_id=artifact.artifact_id)
        except ArtifactDiskRelocationError as exc:
            assert exc.reason_code == "relocation_payload_mismatch"
        else:
            raise AssertionError("verify accepted tampered durable bytes")
        assert relocation.load(run_id, artifact.artifact_id).stage == "copied"
        assert read().content == SOURCE_CONTENT

        restored = relocation.copy(run_id=run_id, artifact_id=artifact.artifact_id)
        assert restored.stage == "copied"
        verified = relocation.resume(run_id=run_id, artifact_id=artifact.artifact_id)
        assert verified.stage == "verified"
        assert read().content == SOURCE_CONTENT
        switched = relocation.resume(run_id=run_id, artifact_id=artifact.artifact_id)
        assert switched.stage == "switched"
        assert source.read_bytes() == SOURCE_CONTENT
        assert (
            relocation.resume(run_id=run_id, artifact_id=artifact.artifact_id)
            == switched
        )

        moved = read()
        assert moved.status_code == 200
        assert moved.content == SOURCE_CONTENT
        assert sha256(moved.content).hexdigest() == sealed.content_sha256
        state = run_service.get_result_state(run_id)
        current = run_service.get_artifact(run_id, artifact.artifact_id)
        assert state.artifact_generation == generation
        assert current.revision == artifact.revision

        source.write_bytes(b"scratch-changed\n")
        assert read().content == SOURCE_CONTENT
        os.chmod(payload, 0o660)
        permission = read()
        assert permission.status_code == 409
        assert _reason(permission) == PUBLIC_PERMISSION_MISMATCH
        _assert_private(permission, scratch, durable)
        os.chmod(payload, SOURCE_MODE)
        assert read().content == SOURCE_CONTENT

        payload.write_bytes(b"durable-changed\n")
        os.chmod(payload, SOURCE_MODE)
        mismatch = read()
        assert mismatch.status_code == 409
        assert _reason(mismatch) == PUBLIC_BYTES_MISMATCH
        _assert_private(mismatch, scratch, durable)
        payload.write_bytes(SOURCE_CONTENT)
        os.chmod(payload, SOURCE_MODE)
        assert read().content == SOURCE_CONTENT

        secret = durable / "secret.txt"
        secret.write_bytes(b"SECRET-BYTES")
        payload.unlink()
        payload.symlink_to(secret)
        unsafe = read()
        assert unsafe.status_code == 409
        assert _reason(unsafe) == PUBLIC_UNAVAILABLE
        assert b"SECRET-BYTES" not in unsafe.content
        _assert_private(unsafe, scratch, durable)
        payload.unlink()
        payload.write_bytes(SOURCE_CONTENT)
        os.chmod(payload, SOURCE_MODE)
        assert read().content == SOURCE_CONTENT

        source.unlink()
        without_scratch = read()
        assert without_scratch.status_code == 200
        assert without_scratch.content == SOURCE_CONTENT
        anonymous_after = _anonymous(app, saved_url, saved_params)
        assert anonymous_after.status_code == 401
        _assert_private(anonymous_after, scratch, durable)
        wrong_revision = dict(saved_params)
        wrong_revision["revision"] = "artifactrev-" + "ef" * 32
        rejected = client.get(saved_url, params=wrong_revision)
        assert rejected.status_code == 409
        assert _reason(rejected) == "ARTIFACT_DOWNLOAD_REVISION_CHANGED"
    finally:
        _close(app)


def _capture_small(tmp_path: Path) -> tuple[ArtifactDiskRelocation, RelocationJournal]:
    source_root = tmp_path / "scratch"
    durable_root = tmp_path / "durable"
    target = source_root / "run-1" / RELATIVE_PATH
    target.parent.mkdir(parents=True)
    target.write_bytes(SOURCE_CONTENT)
    os.chmod(target, SOURCE_MODE)
    revision, path_identity = _current_binding(target)
    relocation = ArtifactDiskRelocation(
        source_root=source_root,
        durable_root=durable_root,
    )
    journal = relocation.capture(
        run_id="run-1",
        artifact_id="artifact-1",
        generation="artifactgen-" + "ab" * 32,
        revision=revision,
        relative_path=RELATIVE_PATH,
        output_type=OUTPUT_TYPE,
        path_identity=path_identity,
        expected_size=len(SOURCE_CONTENT),
        content=SOURCE_CONTENT,
        mode=SOURCE_MODE,
    )
    return relocation, journal


def test_directory_layers_reject_symlinks_before_writing_outside(
    tmp_path: Path,
) -> None:
    outside = tmp_path / "outside"
    outside.mkdir()
    durable = tmp_path / "durable"
    durable.mkdir()
    (durable / "journals").symlink_to(outside, target_is_directory=True)
    source_root = tmp_path / "scratch"
    target = source_root / "run-1" / RELATIVE_PATH
    target.parent.mkdir(parents=True)
    target.write_bytes(SOURCE_CONTENT)
    os.chmod(target, SOURCE_MODE)
    revision, path_identity = _current_binding(target)
    relocation = ArtifactDiskRelocation(
        source_root=source_root,
        durable_root=durable,
    )
    try:
        relocation.capture(
            run_id="run-1",
            artifact_id="artifact-1",
            generation="artifactgen-" + "ab" * 32,
            revision=revision,
            relative_path=RELATIVE_PATH,
            output_type=OUTPUT_TYPE,
            path_identity=path_identity,
            expected_size=len(SOURCE_CONTENT),
            content=SOURCE_CONTENT,
            mode=SOURCE_MODE,
        )
    except ArtifactDiskRelocationError as exc:
        assert exc.reason_code == "relocation_source_unavailable"
    else:
        raise AssertionError("journal write followed a symlinked journals directory")
    assert list(outside.iterdir()) == []

    (durable / "journals").unlink()
    relocation, _journal = _capture_small(tmp_path / "payload-case")
    payload_root = tmp_path / "payload-case" / "durable" / "payload" / "run-1"
    payload_root.mkdir(parents=True)
    (payload_root / "results").symlink_to(outside, target_is_directory=True)
    try:
        relocation.resume(run_id="run-1", artifact_id="artifact-1")
    except ArtifactDiskRelocationError as exc:
        assert exc.reason_code == "relocation_source_unavailable"
    else:
        raise AssertionError("copy followed a symlinked results directory")
    assert list(outside.iterdir()) == []
    assert (
        tmp_path / "payload-case" / "scratch" / "run-1" / RELATIVE_PATH
    ).read_bytes() == (SOURCE_CONTENT)


def test_switched_fifo_fails_without_blocking(tmp_path: Path) -> None:
    relocation, journal = _capture_small(tmp_path)
    assert relocation.resume(run_id="run-1", artifact_id="artifact-1").stage == "copied"
    assert (
        relocation.resume(run_id="run-1", artifact_id="artifact-1").stage == "verified"
    )
    assert (
        relocation.resume(run_id="run-1", artifact_id="artifact-1").stage == "switched"
    )
    payload = relocation.payload_file("run-1", "artifact-1")
    payload.unlink()
    os.mkfifo(payload)
    holder: dict[str, object] = {}

    def read() -> None:
        try:
            holder["decision"] = lookup_switched_download(
                tmp_path / "durable",
                run_id="run-1",
                artifact_id="artifact-1",
                generation=journal.generation,
                revision=journal.revision,
                relative_path=RELATIVE_PATH,
                size_bytes=journal.size_bytes,
            )
        except Exception as exc:  # noqa: BLE001 - the thread must not hang the test
            holder["error"] = exc

    worker = threading.Thread(target=read, daemon=True)
    worker.start()
    worker.join(2)
    assert not worker.is_alive()
    assert "error" not in holder
    decision = holder["decision"]
    assert getattr(decision, "kind") == "reject"
    assert getattr(decision, "reason_code") == PUBLIC_UNAVAILABLE


def test_capture_rejects_oversized_file_before_reading_it(
    tmp_path: Path,
    monkeypatch,
) -> None:
    reads: list[int] = []

    def refuse_read(nodes, expected_size: int) -> bytes:
        del nodes
        reads.append(expected_size)
        raise AssertionError("relocation read past the size cap")

    monkeypatch.setattr(
        artifact_downloads_module,
        "_read_bounded_content",
        refuse_read,
    )
    scratch = tmp_path / "scratch"
    durable = tmp_path / "durable"
    project = scratch / "project"
    _project(project)
    candidate = ExtractedArtifactCandidate(
        output_type="macs3_peak",
        relative_path=RELATIVE_PATH,
        mime_type="text/tab-separated-values",
        metadata={"catalog_id": "macs3_peak", "scope": "project"},
    )
    adapter = ArtifactAdapter((candidate,))
    adapter.identity_root = project
    registry = WorkflowRegistry([adapter])
    app = create_app(
        database_url=f"sqlite:///{scratch / 'platform.db'}",
        workspace_root=scratch / "workspaces",
        project_root=project,
        registry=registry,
        relocation_root=durable,
    )
    try:
        seed_test_authentication(app)
        run_service = app.state.run_service
        run_id = run_service.create_run("fake", WorkflowInputs(config={})).run_id
        run_service.transition_run(run_id, RunStatus.VALIDATING)
        identity = app.state.build_identity_provider.capture("fake").value
        run_service.complete_preflight(run_id, identity)
        for status in (RunStatus.QUEUED, RunStatus.RUNNING, RunStatus.SUCCEEDED):
            run_service.transition_run(run_id, status)
        source = app.state.workspace_root / run_id / RELATIVE_PATH
        source.parent.mkdir(parents=True)
        source.write_bytes(b"")
        os.truncate(source, MAX_RELOCATION_BYTES + 1)
        os.chmod(source, SOURCE_MODE)
        extraction = ArtifactExtractionService(
            run_service=run_service,
            registry=registry,
            build_identity_provider=app.state.build_identity_provider,
            workspace_root=app.state.workspace_root,
        )
        extracted = extraction.extract(run_id)
        assert extracted.is_success, extracted.issues
        artifact = extracted.value[0]
        generation = run_service.get_result_state(run_id).artifact_generation
        assert generation is not None
        try:
            app.state.artifact_download_service.capture_relocation(
                ArtifactDiskRelocation(
                    source_root=app.state.workspace_root,
                    durable_root=durable,
                ),
                run_id=run_id,
                artifact_id=artifact.artifact_id,
                expected_generation=generation,
                expected_revision=artifact.revision,
            )
        except ArtifactDiskRelocationError as exc:
            assert exc.reason_code == "relocation_too_large"
        else:
            raise AssertionError("oversized artifact was copied")
        assert reads == []
        journal = durable / "journals" / run_id / f"{artifact.artifact_id}.json"
        assert not journal.exists()
    finally:
        _close(app)


def _switch_small(tmp_path: Path) -> tuple[ArtifactDiskRelocation, RelocationJournal]:
    relocation, journal = _capture_small(tmp_path)
    for stage in ("copied", "verified", "switched"):
        assert (
            relocation.resume(run_id="run-1", artifact_id="artifact-1").stage == stage
        )
    return relocation, journal


def _lookup_small(tmp_path: Path, journal: RelocationJournal) -> DownloadRelocation:
    return lookup_switched_download(
        tmp_path / "durable",
        run_id="run-1",
        artifact_id="artifact-1",
        generation=journal.generation,
        revision=journal.revision,
        relative_path=RELATIVE_PATH,
        size_bytes=journal.size_bytes,
    )


def test_verify_and_download_do_not_follow_managed_symlinks(tmp_path: Path) -> None:
    relocation, journal = _capture_small(tmp_path)
    assert relocation.resume(run_id="run-1", artifact_id="artifact-1").stage == "copied"
    durable = tmp_path / "durable"
    outside_payload = tmp_path / "outside-payload"
    copied = outside_payload / "run-1" / RELATIVE_PATH
    copied.parent.mkdir(parents=True)
    copied.write_bytes(SOURCE_CONTENT)
    os.chmod(copied, SOURCE_MODE)
    payload = durable / "payload"
    held = durable / "payload-real"
    payload.rename(held)
    payload.symlink_to(outside_payload, target_is_directory=True)
    try:
        relocation.resume(run_id="run-1", artifact_id="artifact-1")
    except ArtifactDiskRelocationError as exc:
        assert exc.reason_code == "relocation_source_unavailable"
    else:
        raise AssertionError("verify followed a symlinked payload directory")
    payload.unlink()
    held.rename(payload)
    assert (
        relocation.resume(run_id="run-1", artifact_id="artifact-1").stage == "verified"
    )
    assert (
        relocation.resume(run_id="run-1", artifact_id="artifact-1").stage == "switched"
    )

    leaked = tmp_path / "leaked-payload"
    leaked_file = leaked / "run-1" / RELATIVE_PATH
    leaked_file.parent.mkdir(parents=True)
    leaked_file.write_bytes(b"chr1\t1\t2\toutside\t1\t+\n")
    os.chmod(leaked_file, SOURCE_MODE)
    payload.rename(held)
    payload.symlink_to(leaked, target_is_directory=True)
    decision = _lookup_small(tmp_path, journal)
    assert decision.kind == "reject"
    assert decision.reason_code == PUBLIC_UNAVAILABLE
    assert b"outside" not in decision.content
    payload.unlink()
    held.rename(payload)
    assert _lookup_small(tmp_path, journal).kind == "serve"

    outside_journals = tmp_path / "outside-journals"
    shutil.copytree(durable / "journals", outside_journals)
    journals = durable / "journals"
    journals_held = durable / "journals-real"
    journals.rename(journals_held)
    journals.symlink_to(outside_journals, target_is_directory=True)
    followed = _lookup_small(tmp_path, journal)
    assert followed.kind == "reject"
    assert followed.reason_code == PUBLIC_UNAVAILABLE
    try:
        relocation.load("run-1", "artifact-1")
    except ArtifactDiskRelocationError as exc:
        assert exc.reason_code == "relocation_journal_invalid"
    else:
        raise AssertionError("resume read followed a symlinked journals directory")


def test_switched_journal_fifo_fails_without_blocking(tmp_path: Path) -> None:
    _relocation, journal = _switch_small(tmp_path)
    journal_path = tmp_path / "durable" / "journals" / "run-1" / "artifact-1.json"
    journal_path.unlink()
    os.mkfifo(journal_path)
    holder: dict[str, object] = {}

    def read() -> None:
        try:
            holder["decision"] = _lookup_small(tmp_path, journal)
        except Exception as exc:  # noqa: BLE001 - the thread must not hang the test
            holder["error"] = exc

    worker = threading.Thread(target=read, daemon=True)
    worker.start()
    worker.join(2)
    assert not worker.is_alive()
    assert "error" not in holder
    decision = holder["decision"]
    assert getattr(decision, "kind") == "reject"
    assert getattr(decision, "reason_code") == PUBLIC_UNAVAILABLE
