"""Security and descriptor-lifecycle tests for artifact downloads."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import replace
from datetime import datetime, timezone
import os
from pathlib import Path
from pathlib import PurePosixPath
import stat
from types import SimpleNamespace

import pytest

from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.result_generations import (
    ARTIFACT_PATH_IDENTITY_METADATA_KEY,
    build_artifact_content_revision,
    build_artifact_descriptor_revision,
    build_artifact_path_identity,
)
from encode_pipeline.platform.runs import RunArtifactRef, RunStatus
from encode_pipeline.services import artifact_downloads
from encode_pipeline.services.artifact_downloads import (
    ArtifactDownloadService,
    ArtifactDownloadStreamError,
)
from encode_pipeline.services.defaults import create_default_workflow_registry
from encode_pipeline.services.run_repositories import ResultGenerationChangedError
from encode_pipeline.services.runs import RunService


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
CHUNK_SIZE = 64 * 1024
WINDOWS_SEPARATOR = chr(92)
ARTIFACT_REVISION = f"artifactrev-{'0' * 64}"


def _runtime(tmp_path: Path, *run_ids: str):
    identifiers = iter(run_ids or ("run-1",))
    run_service = RunService(
        create_default_workflow_registry(),
        id_factory=lambda: next(identifiers),
    )
    for _run_id in run_ids or ("run-1",):
        record = run_service.create_run(WORKFLOW_ID, WorkflowInputs(config={}))
        assert record.run_id == _run_id
        for status in (
            RunStatus.VALIDATING,
            RunStatus.PLANNED,
            RunStatus.QUEUED,
            RunStatus.RUNNING,
            RunStatus.SUCCEEDED,
        ):
            run_service.transition_run(_run_id, status)
    workspace_root = tmp_path / "workspaces"
    workspace_root.mkdir(parents=True)
    service = ArtifactDownloadService(
        run_service=run_service,
        workspace_root=workspace_root,
    )
    return service, run_service, workspace_root


def _artifact(
    run_id: str,
    *,
    artifact_id: str = "artifact-file",
    relative_path: str = "results/file.bin",
    size_bytes: object = 7,
    name: str = "file.bin",
    mime_type: str | None = "application/octet-stream",
    artifact_type: str = "file",
    uri: str | None = None,
) -> RunArtifactRef:
    return RunArtifactRef(
        artifact_id=artifact_id,
        run_id=run_id,
        artifact_type=artifact_type,
        name=name,
        uri=uri or f"run://runs/{run_id}/artifacts/{artifact_id}",
        mime_type=mime_type,
        produced_at=datetime.now(timezone.utc),
        revision=ARTIFACT_REVISION,
        metadata={
            "relative_path": relative_path,
            "output_type": "summary",
            "size_bytes": size_bytes,
        },
    )


def _write(workspace_root: Path, run_id: str, relative_path: str, content: bytes):
    target = workspace_root / run_id / relative_path
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(content)
    return target


def _record(run_service: RunService, artifact: RunArtifactRef) -> None:
    run_service.replace_artifacts(
        artifact.run_id,
        (*run_service.list_artifacts(artifact.run_id), artifact),
    )


def _bind_artifact(
    artifact: RunArtifactRef,
    workspace_root: Path,
    *,
    content_bound: bool = True,
) -> RunArtifactRef:
    relative_path = str(artifact.metadata["relative_path"])
    target = workspace_root / artifact.run_id / relative_path
    components = (
        *(workspace_root / artifact.run_id).parts[1:],
        *PurePosixPath(relative_path).parts,
    )
    current = Path("/")
    infos: list[os.stat_result] = []
    for component in components:
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
    path_identity = build_artifact_path_identity(
        parent_identities=tuple(
            (info.st_dev, info.st_ino, info.st_mode) for info in infos[:-1]
        ),
        file_identity=file_identity,
    )
    output_type = str(artifact.metadata["output_type"])
    revision = (
        build_artifact_content_revision(
            output_type=output_type,
            relative_path=relative_path,
            content=target.read_bytes(),
        )
        if content_bound
        else build_artifact_descriptor_revision(
            output_type=output_type,
            relative_path=relative_path,
            descriptor_identity=":".join(str(value) for value in file_identity),
        )
    )
    return replace(
        artifact,
        revision=revision,
        metadata={
            **artifact.metadata,
            ARTIFACT_PATH_IDENTITY_METADATA_KEY: path_identity,
        },
    )


def _prepare(
    service: ArtifactDownloadService,
    run_service: RunService,
    run_id: str,
    artifact_id: str,
):
    generation = "artifactgen-" + "a" * 64
    revision = ARTIFACT_REVISION
    try:
        state = run_service.get_result_state(run_id)
        if state.artifact_generation is not None:
            generation = state.artifact_generation
        revision = run_service.get_artifact(run_id, artifact_id).revision
    except KeyError:
        pass
    return service.prepare(
        run_id,
        artifact_id,
        expected_generation=generation,
        expected_revision=revision,
    )


def _failure_code(result) -> str:
    assert result.is_failure
    assert result.value is None
    assert len(result.issues) == 1
    issue = result.issues[0]
    assert issue.technical_message is None
    return issue.code


def test_prepare_streams_only_persisted_bytes_and_closes_every_descriptor(
    tmp_path,
    monkeypatch,
):
    content = b"artifact-download\n"
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact(
        "run-1",
        relative_path="results/report data.tsv",
        size_bytes=len(content),
        name="report data.tsv",
        mime_type="text/tab-separated-values",
    )
    _write(workspace_root, "run-1", artifact.metadata["relative_path"], content)
    artifact = _bind_artifact(artifact, workspace_root)
    _record(run_service, artifact)
    original_open = os.open
    opened: list[int] = []

    def tracking_open(*args, **kwargs):
        descriptor = original_open(*args, **kwargs)
        opened.append(descriptor)
        return descriptor

    monkeypatch.setattr(artifact_downloads.os, "open", tracking_open)
    events_before = run_service.list_events("run-1")

    result = _prepare(service, run_service, "run-1", artifact.artifact_id)

    assert result.is_success
    plan = result.value
    assert plan.size_bytes == len(content)
    assert plan.media_type == "text/tab-separated-values"
    assert plan.filename == "report data.tsv"
    assert "\r" not in plan.content_disposition
    assert "\n" not in plan.content_disposition
    assert "filename*=UTF-8''report%20data.tsv" in plan.content_disposition
    assert b"".join(plan.iter_bytes()) == content
    assert plan.closed is True
    assert run_service.list_events("run-1") == events_before
    assert opened
    for descriptor in opened:
        with pytest.raises(OSError):
            os.fstat(descriptor)


def test_prepare_is_run_scoped_and_never_opens_cross_run_artifact(
    tmp_path, monkeypatch
):
    service, run_service, workspace_root = _runtime(tmp_path, "run-1", "run-2")
    other = _artifact("run-2", artifact_id="artifact-other")
    _record(run_service, _artifact("run-1", artifact_id="artifact-present"))
    _record(run_service, other)
    _write(workspace_root, "run-2", "results/file.bin", b"private")
    original_open = os.open
    open_calls: list[tuple[tuple[object, ...], dict[str, object]]] = []

    def tracking_open(*args, **kwargs):
        open_calls.append((args, kwargs))
        return original_open(*args, **kwargs)

    monkeypatch.setattr(artifact_downloads.os, "open", tracking_open)

    missing_run = _prepare(service, run_service, "run-missing", "artifact-other")
    missing_artifact = _prepare(service, run_service, "run-1", "artifact-missing")
    cross_run = _prepare(service, run_service, "run-1", "artifact-other")

    assert _failure_code(missing_run) == "RUN_NOT_FOUND"
    assert _failure_code(missing_artifact) == "RUN_ARTIFACT_NOT_FOUND"
    assert _failure_code(cross_run) == "RUN_ARTIFACT_NOT_FOUND"
    assert missing_artifact.issues == cross_run.issues
    assert open_calls == []


@pytest.mark.parametrize(
    "relative_path",
    [
        "",
        "/private/file.bin",
        "~/private/file.bin",
        f"C:{WINDOWS_SEPARATOR}private{WINDOWS_SEPARATOR}file.bin",
        "results\\file.bin",
        "results/../file.bin",
        "results/./file.bin",
        "results//file.bin",
        "other/file.bin",
        "results/file\x00.bin",
    ],
)
def test_prepare_rejects_persisted_path_attacks_without_disclosure(
    tmp_path,
    relative_path,
):
    service, run_service, _workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", relative_path=relative_path)
    _record(run_service, artifact)

    result = _prepare(service, run_service, "run-1", artifact.artifact_id)

    assert _failure_code(result) == "RUN_ARTIFACT_DOWNLOAD_DATA_INVALID"
    serialized = str(result.issues[0].to_dict())
    assert str(tmp_path) not in serialized
    if relative_path:
        assert relative_path not in serialized


@pytest.mark.parametrize(
    "mutation",
    [
        {"artifact_type": "directory"},
        {"name": "../file.bin"},
        {"name": "file\nprivate.bin"},
        {"mime_type": "text/plain\r\nX-Private: yes"},
        {"uri": "file:///private/workspace/file.bin"},
    ],
)
def test_prepare_rejects_corrupt_persisted_identity_and_headers(tmp_path, mutation):
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = replace(_artifact("run-1"), **mutation)
    _record(run_service, artifact)
    _write(workspace_root, "run-1", "results/file.bin", b"content")

    result = _prepare(service, run_service, "run-1", artifact.artifact_id)

    assert _failure_code(result) == "RUN_ARTIFACT_DOWNLOAD_DATA_INVALID"
    assert "private" not in str(result.issues[0].to_dict()).lower()


@pytest.mark.parametrize("size_bytes", [True, -1, 2**63])
def test_prepare_rejects_invalid_persisted_size(tmp_path, size_bytes):
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", size_bytes=size_bytes)
    _record(run_service, artifact)
    _write(workspace_root, "run-1", "results/file.bin", b"content")

    result = _prepare(service, run_service, "run-1", artifact.artifact_id)

    assert _failure_code(result) == "RUN_ARTIFACT_DOWNLOAD_DATA_INVALID"


@pytest.mark.parametrize("actual", [b"short", b"content-too-long"])
def test_prepare_rejects_larger_or_smaller_files_as_conflicts(tmp_path, actual):
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", size_bytes=7)
    source = _write(workspace_root, "run-1", "results/file.bin", b"content")
    artifact = _bind_artifact(artifact, workspace_root)
    _record(run_service, artifact)
    source.write_bytes(actual)

    result = _prepare(service, run_service, "run-1", artifact.artifact_id)

    assert _failure_code(result) == "RUN_ARTIFACT_DOWNLOAD_CONFLICT"
    assert result.issues[0].context == {
        "reason_code": "ARTIFACT_DOWNLOAD_SIZE_MISMATCH"
    }


def test_prepare_rejects_same_length_content_rewrite_against_persisted_revision(
    tmp_path,
):
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", size_bytes=7)
    source = _write(workspace_root, "run-1", "results/file.bin", b"indexed")
    artifact = _bind_artifact(artifact, workspace_root)
    _record(run_service, artifact)
    source.write_bytes(b"changed")

    result = _prepare(service, run_service, "run-1", artifact.artifact_id)

    assert _failure_code(result) == "RUN_ARTIFACT_DOWNLOAD_CONFLICT"
    assert result.issues[0].context == {
        "reason_code": "ARTIFACT_DOWNLOAD_REVISION_MISMATCH"
    }


def test_prepare_rejects_same_content_inode_replacement(tmp_path):
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", size_bytes=7)
    source = _write(workspace_root, "run-1", "results/file.bin", b"content")
    artifact = _bind_artifact(artifact, workspace_root)
    _record(run_service, artifact)
    replacement = source.with_name("replacement.bin")
    replacement.write_bytes(b"content")
    replacement.replace(source)

    result = _prepare(service, run_service, "run-1", artifact.artifact_id)

    assert _failure_code(result) == "RUN_ARTIFACT_DOWNLOAD_CONFLICT"
    assert result.issues[0].context == {
        "reason_code": "ARTIFACT_DOWNLOAD_PATH_IDENTITY_MISMATCH"
    }


def test_prepare_rejects_content_mutation_during_bounded_revision_read(
    tmp_path,
    monkeypatch,
):
    content = b"indexed"
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", size_bytes=len(content))
    target = _write(workspace_root, "run-1", "results/file.bin", content)
    artifact = _bind_artifact(artifact, workspace_root)
    _record(run_service, artifact)
    original_read = os.read
    mutated = False

    def racing_read(descriptor, size):
        nonlocal mutated
        value = original_read(descriptor, size)
        if value and not mutated:
            mutated = True
            target.write_bytes(b"changed")
        return value

    monkeypatch.setattr(artifact_downloads.os, "read", racing_read)

    result = _prepare(service, run_service, "run-1", artifact.artifact_id)

    assert _failure_code(result) == "RUN_ARTIFACT_DOWNLOAD_CONFLICT"
    assert result.issues[0].context == {
        "reason_code": "ARTIFACT_DOWNLOAD_SOURCE_CHANGED"
    }


@pytest.mark.parametrize("content_bound", [True, False])
def test_prepare_rejects_path_replacement_during_final_generation_guard(
    tmp_path,
    monkeypatch,
    content_bound,
):
    content = b"indexed"
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", size_bytes=len(content))
    target = _write(workspace_root, "run-1", "results/file.bin", content)
    artifact = _bind_artifact(
        artifact,
        workspace_root,
        content_bound=content_bound,
    )
    _record(run_service, artifact)
    original_get = run_service.get_artifact_at_generation
    calls = 0

    def replace_during_final_guard(*args, **kwargs):
        nonlocal calls
        result = original_get(*args, **kwargs)
        calls += 1
        if calls == 2:
            replacement = target.with_name("replacement.bin")
            replacement.write_bytes(b"changed")
            replacement.replace(target)
        return result

    monkeypatch.setattr(
        run_service,
        "get_artifact_at_generation",
        replace_during_final_guard,
    )

    result = _prepare(service, run_service, "run-1", artifact.artifact_id)

    assert _failure_code(result) == "RUN_ARTIFACT_DOWNLOAD_CONFLICT"
    assert result.issues[0].context == {
        "reason_code": "ARTIFACT_DOWNLOAD_SOURCE_CHANGED"
    }


@pytest.mark.parametrize(
    ("failure", "expected_code"),
    [
        (
            ResultGenerationChangedError("private generation failure"),
            "RUN_ARTIFACT_DOWNLOAD_CONFLICT",
        ),
        (KeyError("private missing row"), "RUN_ARTIFACT_DOWNLOAD_CONFLICT"),
        (ValueError("private corrupt row"), "RUN_ARTIFACT_DOWNLOAD_DATA_INVALID"),
        (
            RuntimeError("private repository failure"),
            "RUN_ARTIFACT_DOWNLOAD_UNAVAILABLE",
        ),
    ],
)
def test_prepare_closes_descriptor_chain_when_final_repository_guard_fails(
    tmp_path,
    monkeypatch,
    failure,
    expected_code,
):
    content = b"content"
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", size_bytes=len(content))
    _write(workspace_root, "run-1", "results/file.bin", content)
    artifact = _bind_artifact(artifact, workspace_root, content_bound=False)
    _record(run_service, artifact)
    original_get = run_service.get_artifact_at_generation
    original_open = os.open
    calls = 0
    opened: list[int] = []

    def fail_final_guard(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 2:
            raise failure
        return original_get(*args, **kwargs)

    def tracking_open(*args, **kwargs):
        descriptor = original_open(*args, **kwargs)
        opened.append(descriptor)
        return descriptor

    monkeypatch.setattr(run_service, "get_artifact_at_generation", fail_final_guard)
    monkeypatch.setattr(artifact_downloads.os, "open", tracking_open)

    result = _prepare(service, run_service, "run-1", artifact.artifact_id)

    assert _failure_code(result) == expected_code
    assert opened
    for descriptor in opened:
        with pytest.raises(OSError):
            os.fstat(descriptor)
    assert "private" not in str(result.issues[0].to_dict()).lower()


def test_prepare_rejects_missing_directory_fifo_and_symlink_sources(tmp_path):
    scenarios: list[tuple[str, Callable[[Path], object]]] = [
        ("missing", lambda target: None),
        ("directory", lambda target: target.mkdir(parents=True)),
        (
            "fifo",
            lambda target: (
                target.parent.mkdir(parents=True, exist_ok=True),
                os.mkfifo(target),
            ),
        ),
        (
            "symlink",
            lambda target: (
                target.parent.mkdir(parents=True, exist_ok=True),
                (target.parent / "real.bin").write_bytes(b"content"),
                target.symlink_to(target.parent / "real.bin"),
            ),
        ),
    ]
    for index, (_name, arrange) in enumerate(scenarios):
        case_root = tmp_path / f"case-{index}"
        service, run_service, workspace_root = _runtime(case_root, "run-1")
        artifact = _artifact("run-1", size_bytes=7)
        target = _write(workspace_root, "run-1", "results/file.bin", b"content")
        artifact = _bind_artifact(artifact, workspace_root)
        _record(run_service, artifact)
        target.unlink()
        arrange(target)

        result = _prepare(service, run_service, "run-1", artifact.artifact_id)

        assert _failure_code(result) == "RUN_ARTIFACT_DOWNLOAD_CONFLICT"
        assert set(result.issues[0].context) == {"reason_code"}
        assert str(case_root) not in str(result.issues[0].to_dict())


def test_prepare_rejects_device_mode_even_when_descriptor_open_succeeds(
    tmp_path,
    monkeypatch,
):
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", size_bytes=7)
    _write(workspace_root, "run-1", "results/file.bin", b"content")
    artifact = _bind_artifact(artifact, workspace_root)
    _record(run_service, artifact)
    original_fstat = os.fstat

    def device_fstat(descriptor):
        info = original_fstat(descriptor)
        if stat.S_ISREG(info.st_mode):
            return SimpleNamespace(
                st_dev=info.st_dev,
                st_ino=info.st_ino,
                st_mode=stat.S_IFCHR | 0o600,
                st_nlink=info.st_nlink,
                st_uid=info.st_uid,
                st_gid=info.st_gid,
                st_size=info.st_size,
                st_mtime_ns=info.st_mtime_ns,
                st_ctime_ns=info.st_ctime_ns,
            )
        return info

    monkeypatch.setattr(artifact_downloads.os, "fstat", device_fstat)

    result = _prepare(service, run_service, "run-1", artifact.artifact_id)

    assert _failure_code(result) == "RUN_ARTIFACT_DOWNLOAD_CONFLICT"
    assert result.issues[0].context == {
        "reason_code": "ARTIFACT_DOWNLOAD_SOURCE_TYPE_INVALID"
    }


def test_prepare_rejects_symlinked_workspace_run_and_parent_components(tmp_path):
    for component in ("workspace", "run", "parent"):
        case_root = tmp_path / component
        service, run_service, workspace_root = _runtime(case_root, "run-1")
        artifact = _artifact("run-1", size_bytes=7)
        _write(workspace_root, "run-1", "results/file.bin", b"content")
        artifact = _bind_artifact(artifact, workspace_root)
        _record(run_service, artifact)
        if component == "workspace":
            real_root = case_root / "real-workspaces"
            workspace_root.rename(real_root)
            workspace_root.symlink_to(real_root, target_is_directory=True)
        elif component == "run":
            original_run = workspace_root / "run-1"
            real_run = case_root / "real-run"
            original_run.rename(real_run)
            (workspace_root / "run-1").symlink_to(real_run, target_is_directory=True)
        else:
            original_results = workspace_root / "run-1/results"
            real_results = case_root / "real-results"
            original_results.rename(real_results)
            (workspace_root / "run-1/results").symlink_to(
                real_results,
                target_is_directory=True,
            )

        result = _prepare(service, run_service, "run-1", artifact.artifact_id)

        assert _failure_code(result) == "RUN_ARTIFACT_DOWNLOAD_CONFLICT"


def test_stream_rejects_same_size_path_replacement_before_first_byte(tmp_path):
    content = b"content"
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", size_bytes=len(content))
    target = _write(workspace_root, "run-1", "results/file.bin", content)
    artifact = _bind_artifact(artifact, workspace_root, content_bound=False)
    _record(run_service, artifact)
    plan = _prepare(service, run_service, "run-1", artifact.artifact_id).value

    target.unlink()
    target.write_bytes(b"changed")

    stream = plan.iter_bytes()
    with pytest.raises(ArtifactDownloadStreamError):
        next(stream)
    assert plan.closed is True


def test_stream_rejects_parent_path_replacement_before_first_byte(tmp_path):
    content = b"content"
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", size_bytes=len(content))
    target = _write(workspace_root, "run-1", "results/file.bin", content)
    artifact = _bind_artifact(artifact, workspace_root, content_bound=False)
    _record(run_service, artifact)
    plan = _prepare(service, run_service, "run-1", artifact.artifact_id).value

    old_results = target.parent.with_name("results-old")
    target.parent.rename(old_results)
    target.parent.mkdir()
    target.write_bytes(b"changed")

    with pytest.raises(ArtifactDownloadStreamError):
        next(plan.iter_bytes())
    assert plan.closed is True


def test_stream_does_not_yield_a_chunk_mutated_during_read(tmp_path, monkeypatch):
    content = b"a" * (CHUNK_SIZE + 32)
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", size_bytes=len(content))
    target = _write(workspace_root, "run-1", "results/file.bin", content)
    artifact = _bind_artifact(artifact, workspace_root, content_bound=False)
    _record(run_service, artifact)
    plan = _prepare(service, run_service, "run-1", artifact.artifact_id).value
    original_read = os.read
    mutated = False

    def racing_read(descriptor, size):
        nonlocal mutated
        value = original_read(descriptor, size)
        if value and not mutated:
            mutated = True
            target.write_bytes(b"b" * len(content))
            current = target.stat()
            os.utime(
                target,
                ns=(current.st_atime_ns, current.st_mtime_ns + 10_000_000),
            )
        return value

    monkeypatch.setattr(artifact_downloads.os, "read", racing_read)

    with pytest.raises(ArtifactDownloadStreamError):
        next(plan.iter_bytes())
    assert plan.closed is True


def test_stream_allows_unrelated_directory_content_changes(tmp_path):
    content = b"content"
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", size_bytes=len(content))
    target = _write(workspace_root, "run-1", "results/file.bin", content)
    artifact = _bind_artifact(artifact, workspace_root, content_bound=False)
    _record(run_service, artifact)
    plan = _prepare(service, run_service, "run-1", artifact.artifact_id).value

    (workspace_root / "unrelated-run").mkdir()
    (target.parent / "unrelated.txt").write_text("unrelated", encoding="utf-8")

    assert b"".join(plan.iter_bytes()) == content


def test_service_rejects_noncanonical_workspace_root(tmp_path):
    run_service = RunService(create_default_workflow_registry())

    with pytest.raises(ValueError, match="workspace_root"):
        ArtifactDownloadService(
            run_service=run_service,
            workspace_root=tmp_path / "parent" / ".." / "workspaces",
        )


@pytest.mark.parametrize("mutation", ["rewrite", "truncate", "grow"])
def test_stream_detects_file_changes_after_a_partial_chunk(tmp_path, mutation):
    content = b"a" * (CHUNK_SIZE + 32)
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", size_bytes=len(content))
    target = _write(workspace_root, "run-1", "results/file.bin", content)
    artifact = _bind_artifact(artifact, workspace_root, content_bound=False)
    _record(run_service, artifact)
    plan = _prepare(service, run_service, "run-1", artifact.artifact_id).value
    stream = plan.iter_bytes()

    assert next(stream) == content[:CHUNK_SIZE]
    if mutation == "rewrite":
        target.write_bytes(b"b" * len(content))
    elif mutation == "truncate":
        target.write_bytes(b"b")
    else:
        with target.open("ab") as handle:
            handle.write(b"growth")
    current = target.stat()
    os.utime(
        target,
        ns=(current.st_atime_ns, current.st_mtime_ns + 10_000_000),
    )

    with pytest.raises(ArtifactDownloadStreamError):
        next(stream)
    assert plan.closed is True


def test_plan_close_is_idempotent_and_prevents_late_reads(tmp_path, monkeypatch):
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", size_bytes=7)
    _write(workspace_root, "run-1", "results/file.bin", b"content")
    artifact = _bind_artifact(artifact, workspace_root, content_bound=False)
    _record(run_service, artifact)
    original_open = os.open
    opened: list[int] = []

    def tracking_open(*args, **kwargs):
        descriptor = original_open(*args, **kwargs)
        opened.append(descriptor)
        return descriptor

    monkeypatch.setattr(artifact_downloads.os, "open", tracking_open)
    plan = _prepare(service, run_service, "run-1", artifact.artifact_id).value

    plan.close()
    plan.close()

    assert plan.closed is True
    with pytest.raises(ArtifactDownloadStreamError):
        next(plan.iter_bytes())
    for descriptor in opened:
        with pytest.raises(OSError):
            os.fstat(descriptor)


def test_client_stopping_iteration_closes_every_descriptor(tmp_path):
    content = b"a" * (CHUNK_SIZE + 1)
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", size_bytes=len(content))
    _write(workspace_root, "run-1", "results/file.bin", content)
    artifact = _bind_artifact(artifact, workspace_root, content_bound=False)
    _record(run_service, artifact)
    plan = _prepare(service, run_service, "run-1", artifact.artifact_id).value
    stream = plan.iter_bytes()

    assert next(stream) == content[:CHUNK_SIZE]
    stream.close()

    assert plan.closed is True


@pytest.mark.parametrize("fail_at", [1, 2])
def test_prepare_closes_new_descriptor_when_fstat_fails(
    tmp_path,
    monkeypatch,
    fail_at,
):
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", size_bytes=7)
    _write(workspace_root, "run-1", "results/file.bin", b"content")
    artifact = _bind_artifact(artifact, workspace_root, content_bound=False)
    _record(run_service, artifact)
    original_open = os.open
    original_fstat = os.fstat
    opened: list[int] = []
    fstat_calls = 0

    def tracking_open(*args, **kwargs):
        descriptor = original_open(*args, **kwargs)
        opened.append(descriptor)
        return descriptor

    def failing_fstat(descriptor):
        nonlocal fstat_calls
        fstat_calls += 1
        if fstat_calls == fail_at:
            raise OSError("private fstat failure")
        return original_fstat(descriptor)

    monkeypatch.setattr(artifact_downloads.os, "open", tracking_open)
    monkeypatch.setattr(artifact_downloads.os, "fstat", failing_fstat)

    result = _prepare(service, run_service, "run-1", artifact.artifact_id)

    assert _failure_code(result) == "RUN_ARTIFACT_DOWNLOAD_UNAVAILABLE"
    assert opened
    for descriptor in opened:
        with pytest.raises(OSError):
            original_fstat(descriptor)


def test_stream_read_error_is_redacted_and_closes_every_descriptor(
    tmp_path,
    monkeypatch,
):
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", size_bytes=7)
    _write(workspace_root, "run-1", "results/file.bin", b"content")
    artifact = _bind_artifact(artifact, workspace_root, content_bound=False)
    _record(run_service, artifact)
    plan = _prepare(service, run_service, "run-1", artifact.artifact_id).value

    def fail_read(_descriptor, _size):
        raise OSError("private read failure")

    monkeypatch.setattr(artifact_downloads.os, "read", fail_read)

    with pytest.raises(ArtifactDownloadStreamError) as exc_info:
        next(plan.iter_bytes())

    assert "private" not in str(exc_info.value)
    assert plan.closed is True


def test_unexpected_open_failure_is_redacted_as_unavailable(tmp_path, monkeypatch):
    service, run_service, workspace_root = _runtime(tmp_path, "run-1")
    artifact = _artifact("run-1", size_bytes=7)
    _write(workspace_root, "run-1", "results/file.bin", b"content")
    artifact = _bind_artifact(artifact, workspace_root, content_bound=False)
    _record(run_service, artifact)

    def fail_open(*_args, **_kwargs):
        raise OSError(f"private failure in {workspace_root}")

    monkeypatch.setattr(artifact_downloads.os, "open", fail_open)

    result = _prepare(service, run_service, "run-1", artifact.artifact_id)

    assert _failure_code(result) == "RUN_ARTIFACT_DOWNLOAD_UNAVAILABLE"
    assert result.issues[0].context == {
        "reason_code": "ARTIFACT_DOWNLOAD_IO_UNAVAILABLE"
    }
    assert str(workspace_root) not in str(result.issues[0].to_dict())
