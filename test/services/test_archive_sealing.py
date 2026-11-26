"""Provenance proofs for sealing a registered revision or a run publication."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import replace
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
from pathlib import PurePosixPath
import tempfile
from types import SimpleNamespace

import pytest

from encode_pipeline.persistence.models import ExperimentRevisionRow, ProjectRow
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.experiment_catalog import catalog_not_found
from encode_pipeline.platform.registered_results import (
    RevisionRecord,
    build_revision_digest,
    parse_metadata,
)
from encode_pipeline.platform.result_generations import (
    ARTIFACT_PATH_IDENTITY_METADATA_KEY,
    build_artifact_content_revision,
    build_artifact_path_identity,
)
from encode_pipeline.platform.runs import RunArtifactRef, RunStatus
from encode_pipeline.services.archive_sealing import (
    ArchiveSealError,
    CatalogExperimentGate,
)
from encode_pipeline.services.archive_sealing import ArchiveSealingService
from encode_pipeline.services.artifact_downloads import ArtifactDownloadService
from encode_pipeline.services.browser_bgzf import write_indexed_bed
from encode_pipeline.services.defaults import create_default_workflow_registry
from encode_pipeline.services.private_storage_pools import (
    PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
)
from encode_pipeline.services.runs import RunService


PROJECT = "prj_" + "a" * 32
EXPERIMENT = "exp_" + "b" * 32
REVISION = "expr_" + "c" * 32
OTHER_REVISION = "expr_" + "d" * 32
RESULT = "rreg_" + "1" * 32
RESULT_REVISION = "rrev_" + "2" * 32
HEAD_REVISION = "rrev_" + "3" * 32
INDEX_RESULT = "rreg_" + "4" * 32
INDEX_REVISION = "rrev_" + "5" * 32
POOL = "stgp_" + "6" * 32
WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


class _Objects:
    def __init__(self) -> None:
        self.payloads: list[bytes] = []

    def seal_fd(self, fd: int, expected_sha256: str | None = None) -> SimpleNamespace:
        os.lseek(fd, 0, os.SEEK_SET)
        payload = bytearray()
        while True:
            chunk = os.read(fd, 65536)
            if not chunk:
                break
            payload.extend(chunk)
        digest = hashlib.sha256(payload).hexdigest()
        if expected_sha256 is not None and digest != expected_sha256:
            raise ArchiveSealError(
                "archive_integrity",
                "The sealed bytes do not match the source.",
            )
        self.payloads.append(bytes(payload))
        return SimpleNamespace(size_bytes=len(payload), sha256=digest)

    @contextmanager
    def open_member(self, member):
        payload = next(
            p for p in self.payloads if hashlib.sha256(p).hexdigest() == member.sha256
        )
        yield SimpleNamespace(read=lambda start, end: payload[start:end])


class _Repository:
    def __init__(self) -> None:
        self.calls: list[dict[str, object]] = []
        self._versions: dict[tuple[object, ...], SimpleNamespace] = {}

    def publish(self, **kwargs: object) -> SimpleNamespace:
        self.calls.append(kwargs)
        members = kwargs["members"]
        assert isinstance(members, dict)
        identity = (
            kwargs["project_id"],
            kwargs["experiment_id"],
            kwargs["track_key"],
            tuple(sorted(kwargs["source"].items())),  # type: ignore[attr-defined]
            tuple((role, member.sha256) for role, member in sorted(members.items())),
        )
        version = self._versions.get(identity)
        if version is None:
            version = SimpleNamespace(
                archive_version_id=f"arc-{len(self._versions) + 1}",
                project_id=kwargs["project_id"],
                experiment_id=kwargs["experiment_id"],
                experiment_revision_id=kwargs["experiment_revision_id"],
                track_key=kwargs["track_key"],
                version_number=len(self._versions) + 1,
                display_name=kwargs["display_name"],
                assembly=kwargs["assembly"],
                track_type=kwargs["track_type"],
                source=dict(kwargs["source"]),  # type: ignore[arg-type]
                has_index="index" in members,
                members=members,
            )
            self._versions[identity] = version
        return version


class _Revisions:
    def __init__(self, *records: RevisionRecord) -> None:
        self._records = {record.result_revision_id: record for record in records}

    def get_revision(self, result_revision_id: str) -> RevisionRecord | None:
        return self._records.get(result_revision_id)


class _Pools:
    def config_key_for(self, storage_pool_id: str) -> str:
        if storage_pool_id != POOL:
            raise ArchiveSealError(
                "archive_offline",
                "The registered storage pool is unavailable.",
            )
        return "primary-inputs"


class _Allow:
    def __init__(self) -> None:
        self.revisions: list[dict[str, str]] = []
        self.publications: list[dict[str, str]] = []

    def require_selected_revision(self, **kwargs: str) -> None:
        self.revisions.append(kwargs)

    def require_run_publication(self, **kwargs: str) -> None:
        self.publications.append(kwargs)


class _Rows:
    def __init__(self, rows: dict[tuple[object, str], object]) -> None:
        self._rows = rows

    def __enter__(self) -> "_Rows":
        return self

    def __exit__(self, *args: object) -> bool:
        return False

    def get(self, model: object, key: str) -> object | None:
        return self._rows.get((model, key))


def _record(
    *,
    result_id: str,
    result_revision_id: str,
    relative_path: str,
    content: bytes,
    head_result_revision_id: str,
    reference: dict[str, object] | None = None,
    result_format: str = "bigwig",
) -> RevisionRecord:
    payload = parse_metadata(
        {
            "display_name": "signal.bw",
            "format": result_format,
            "origin": {"kind": "external", "source": {"status": "unknown"}},
            "method": {"status": "unknown"},
            "reference": reference or {"status": "unknown"},
        }
    ).canonical_json()
    digest = hashlib.sha256(content).hexdigest()
    return RevisionRecord(
        project_id=PROJECT,
        project_archived_at=None,
        result_id=result_id,
        stable_key="signal",
        result_archived_at=None,
        result_revision_id=result_revision_id,
        revision_number=1,
        canonical_payload=payload,
        revision_digest=build_revision_digest(
            project_id=PROJECT,
            result_id=result_id,
            result_revision_id=result_revision_id,
            revision_number=1,
            canonical_payload=payload,
            storage_pool_id=POOL,
            relative_path=relative_path,
            size_bytes=len(content),
            content_sha256=digest,
            parents=(),
        ),
        created_at=datetime.now(timezone.utc),
        storage_pool_id=POOL,
        relative_path=relative_path,
        size_bytes=len(content),
        content_sha256=digest,
        head_result_revision_id=head_result_revision_id,
        parents=(),
    )


def _pool(tmp_path: Path) -> tuple[Path, Path]:
    root = tmp_path / "pool"
    root.mkdir()
    config = tmp_path / "pools.json"
    config.write_text(
        json.dumps(
            {
                "schema_version": PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
                "storage_pool_roots": {"primary-inputs": str(root)},
            }
        ),
        encoding="utf-8",
    )
    return root, config


def _write_pool(root: Path, relative_path: str, content: bytes) -> None:
    target = root.joinpath(*relative_path.split("/"))
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(content)


def _service(
    tmp_path: Path,
    *,
    records: tuple[RevisionRecord, ...],
    experiments: object | None = None,
    assemblies: object | None = None,
) -> tuple[ArchiveSealingService, _Objects, _Repository, Path]:
    root, config = _pool(tmp_path)
    archive = tmp_path / "archive"
    archive.mkdir()
    objects = _Objects()
    repository = _Repository()
    service = ArchiveSealingService(
        archive_root=archive,
        objects=objects,
        repository=repository,
        experiments=experiments or _Allow(),
        registered_revisions=_Revisions(*records),
        pool_lookup=_Pools(),
        pool_config_path=config,
        assemblies=assemblies,
    )
    return service, objects, repository, root


def _seal_registered(
    service: ArchiveSealingService, **overrides: str
) -> dict[str, object]:
    arguments = {
        "project_id": PROJECT,
        "experiment_id": EXPERIMENT,
        "experiment_revision_id": REVISION,
        "track_key": "signal",
        "display_name": "signal.bw",
        "assembly": "GRCh38",
        "track_type": "bigwig",
        "result_id": RESULT,
        "result_revision_id": RESULT_REVISION,
    }
    arguments.update(overrides)
    return service.seal_registered(**arguments)


def test_registered_seal_reads_the_revision_path_not_the_display_name(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    content = b"v1-signal"
    decoy = b"display-name"
    record = _record(
        result_id=RESULT,
        result_revision_id=RESULT_REVISION,
        relative_path="tracks/signal-v1.bw",
        content=content,
        head_result_revision_id=HEAD_REVISION,
    )
    service, objects, repository, root = _service(tmp_path, records=(record,))
    _write_pool(root, "tracks/signal-v1.bw", content)
    _write_pool(root, "signal.bw", decoy)
    created: list[str] = []
    real_mkdtemp = tempfile.mkdtemp

    def tracking(*args: object, **kwargs: object) -> str:
        path = real_mkdtemp(*args, **kwargs)  # type: ignore[arg-type]
        created.append(path)
        return path

    monkeypatch.setattr(tempfile, "mkdtemp", tracking)

    receipt = _seal_registered(service)

    assert objects.payloads == [content]
    assert receipt["assembly"] == "hg38"
    assert receipt["source"] == {
        "kind": "registered",
        "result_id": RESULT,
        "result_revision_id": RESULT_REVISION,
    }
    assert receipt["has_index"] is False
    assert receipt["members"] == {
        "data": {
            "sha256": hashlib.sha256(content).hexdigest(),
            "size_bytes": len(content),
        }
    }
    assert "tracks/signal-v1.bw" not in json.dumps(receipt)
    assert repository.calls[0]["source"]["result_revision_id"] == RESULT_REVISION
    assert created == []
    assert (root / "tracks" / "signal-v1.bw").read_bytes() == content


@pytest.mark.parametrize("declared", ["bed", "bam", "hic", "tbi"])
def test_declared_format_conflict_does_not_read_the_pool_file(
    tmp_path: Path, declared: str
) -> None:
    content = b"bed-not-bigwig"
    record = _record(
        result_id=RESULT,
        result_revision_id=RESULT_REVISION,
        relative_path="tracks/signal-v1.bw",
        content=content,
        head_result_revision_id=RESULT_REVISION,
        result_format=declared,
    )
    service, objects, _repository, root = _service(tmp_path, records=(record,))
    _write_pool(root, "tracks/signal-v1.bw", content)

    with pytest.raises(ArchiveSealError) as captured:
        _seal_registered(service)
    assert captured.value.code == "archive_invalid"
    assert objects.payloads == []
    assert (root / "tracks" / "signal-v1.bw").read_bytes() == content


def test_indexed_longrange_without_an_index_is_refused(tmp_path: Path) -> None:
    content = b"pairs"
    record = _record(
        result_id=RESULT,
        result_revision_id=RESULT_REVISION,
        relative_path="tracks/pairs",
        content=content,
        head_result_revision_id=RESULT_REVISION,
        result_format="indexed-longrange",
    )
    service, objects, _repository, root = _service(tmp_path, records=(record,))
    _write_pool(root, "tracks/pairs", content)

    with pytest.raises(ArchiveSealError) as captured:
        _seal_registered(service, track_type="longrange")
    assert captured.value.code == "archive_invalid"
    assert objects.payloads == []


def test_bound_reference_assembly_is_kept_and_a_conflict_is_refused(
    tmp_path: Path,
) -> None:
    content = b"bound-signal"
    record = _record(
        result_id=RESULT,
        result_revision_id=RESULT_REVISION,
        relative_path="tracks/signal-v1.bw",
        content=content,
        head_result_revision_id=RESULT_REVISION,
        reference={
            "status": "known",
            "kind": "reference_profile_revision",
            "reference_profile_revision_id": "refpr_" + "a" * 32,
            "public_identity_sha256": "ab" * 32,
        },
    )

    class _Bound:
        def registered_assembly(self, reference: dict[str, object]) -> str:
            assert reference["kind"] == "reference_profile_revision"
            return "GRCh38"

    service, objects, _repository, root = _service(
        tmp_path,
        records=(record,),
        assemblies=_Bound(),
    )
    _write_pool(root, "tracks/signal-v1.bw", content)

    receipt = _seal_registered(service)

    assert receipt["assembly"] == "hg38"
    assert objects.payloads == [content]
    with pytest.raises(ArchiveSealError) as captured:
        _seal_registered(service, assembly="mm10")
    assert captured.value.code == "archive_invalid"
    assert objects.payloads == [content]


def test_registered_retry_returns_the_repository_identity_for_the_old_revision(
    tmp_path: Path,
) -> None:
    old = b"v1-signal"
    latest = b"v2-signal-different"
    old_record = _record(
        result_id=RESULT,
        result_revision_id=RESULT_REVISION,
        relative_path="tracks/signal-v1.bw",
        content=old,
        head_result_revision_id=HEAD_REVISION,
    )
    head = _record(
        result_id=RESULT,
        result_revision_id=HEAD_REVISION,
        relative_path="tracks/signal-v2.bw",
        content=latest,
        head_result_revision_id=HEAD_REVISION,
    )
    service, objects, _repository, root = _service(
        tmp_path,
        records=(old_record, head),
    )
    _write_pool(root, "tracks/signal-v1.bw", old)
    _write_pool(root, "tracks/signal-v2.bw", latest)

    first = _seal_registered(service)
    second = _seal_registered(service)

    assert first["archive_version_id"] == second["archive_version_id"] == "arc-1"
    assert objects.payloads == [old, old]
    assert first["source"]["result_revision_id"] == RESULT_REVISION


@pytest.mark.parametrize(
    "index_case", ["matching", "other-contig", "wrong-bin", "truncated", "not-tabix"]
)
def test_registered_index_is_a_separate_revision(
    tmp_path: Path, index_case: str
) -> None:
    indexed = write_indexed_bed(("chr1\t10\t20\tchr1:40-50,1",))
    data, index = indexed.bgzf, indexed.tabix
    if index_case == "other-contig":
        index = write_indexed_bed(("chr2\t10\t20\tchr2:40-50,1",)).tabix
    if index_case == "wrong-bin":
        index = write_indexed_bed(("chr1\t500000\t500010\tchr1:50-60,1",)).tabix
    if index_case == "truncated":
        index = index[: len(index) // 2]
    data_record = _record(
        result_id=RESULT,
        result_revision_id=RESULT_REVISION,
        relative_path="tracks/pairs.bed",
        content=data,
        head_result_revision_id=HEAD_REVISION,
        result_format="longrange",
    )
    index_record = _record(
        result_id=INDEX_RESULT,
        result_revision_id=INDEX_REVISION,
        relative_path="tracks/pairs.bed.tbi",
        content=index,
        head_result_revision_id=INDEX_REVISION,
        result_format="bai" if index_case == "not-tabix" else "tbi",
    )
    service, objects, repository, root = _service(
        tmp_path,
        records=(data_record, index_record),
    )
    _write_pool(root, "tracks/pairs.bed", data)
    _write_pool(root, "tracks/pairs.bed.tbi", index)

    arguments = dict(
        project_id=PROJECT,
        experiment_id=EXPERIMENT,
        experiment_revision_id=REVISION,
        track_key="pairs",
        display_name="pairs.bed",
        assembly="mm10",
        track_type="longrange",
        result_id=RESULT,
        result_revision_id=RESULT_REVISION,
        index_result_id=INDEX_RESULT,
        index_result_revision_id=INDEX_REVISION,
    )
    if index_case != "matching":
        with pytest.raises(ArchiveSealError) as captured:
            service.seal_registered(**arguments)
        assert captured.value.code == "archive_invalid"
        assert repository.calls == []
        return
    receipt = service.seal_registered(**arguments)

    assert objects.payloads == [data, index]
    assert receipt["has_index"] is True
    source = repository.calls[0]["source"]
    assert source["index_result_id"] == INDEX_RESULT
    assert source["index_result_revision_id"] == INDEX_REVISION
    assert set(receipt["members"]) == {"data", "index"}


def test_changed_registered_bytes_are_not_published(tmp_path: Path) -> None:
    content = b"v1-signal"
    record = _record(
        result_id=RESULT,
        result_revision_id=RESULT_REVISION,
        relative_path="tracks/signal-v1.bw",
        content=content,
        head_result_revision_id=RESULT_REVISION,
    )
    service, objects, repository, root = _service(tmp_path, records=(record,))
    _write_pool(root, "tracks/signal-v1.bw", b"v1-OTHER!")

    with pytest.raises(ArchiveSealError) as captured:
        _seal_registered(service)

    assert captured.value.code == "archive_integrity"
    assert objects.payloads == []
    assert repository.calls == []


def test_wrong_experiment_does_not_open_registered_bytes(tmp_path: Path) -> None:
    content = b"v1-signal"
    record = _record(
        result_id=RESULT,
        result_revision_id=RESULT_REVISION,
        relative_path="tracks/missing.bw",
        content=content,
        head_result_revision_id=RESULT_REVISION,
    )

    class _Deny:
        def require_selected_revision(self, **kwargs: str) -> None:
            raise ArchiveSealError(
                "archive_version_not_found",
                "The experiment revision was not found.",
            )

    service, objects, repository, _root = _service(
        tmp_path,
        records=(record,),
        experiments=_Deny(),
    )

    with pytest.raises(ArchiveSealError) as captured:
        _seal_registered(service)

    assert captured.value.code == "archive_version_not_found"
    assert objects.payloads == []
    assert repository.calls == []


def test_filename_and_latest_are_rejected(tmp_path: Path) -> None:
    record = _record(
        result_id=RESULT,
        result_revision_id=RESULT_REVISION,
        relative_path="tracks/signal-v1.bw",
        content=b"v1-signal",
        head_result_revision_id=RESULT_REVISION,
    )
    service, objects, _repository, _root = _service(tmp_path, records=(record,))

    with pytest.raises(ArchiveSealError) as captured:
        _seal_registered(service, result_revision_id="latest")
    assert captured.value.code == "archive_invalid"

    with pytest.raises(TypeError):
        service.seal_registered(  # type: ignore[call-arg]
            project_id=PROJECT,
            experiment_id=EXPERIMENT,
            experiment_revision_id=REVISION,
            track_key="signal",
            display_name="signal.bw",
            assembly="hg38",
            track_type="bigwig",
            result_id=RESULT,
            result_revision_id=RESULT_REVISION,
            path="/tmp/signal.bw",
        )
    assert objects.payloads == []


def test_catalog_gate_uses_the_named_revision_and_publication() -> None:
    revision = SimpleNamespace(
        experiment_id=EXPERIMENT,
        project_id=PROJECT,
        experiment_revision_id=REVISION,
    )
    project = SimpleNamespace(project_id=PROJECT, archived_at=None, kind="user")
    rows = {
        (ExperimentRevisionRow, REVISION): revision,
        (ProjectRow, PROJECT): project,
    }

    class _Catalog:
        def __init__(self) -> None:
            self.traces: list[dict[str, str]] = []

        def get_experiment_revision(
            self,
            experiment_id: str,
            experiment_revision_id: str,
        ) -> dict[str, str]:
            return {
                "experiment_id": experiment_id,
                "experiment_revision_id": experiment_revision_id,
            }

        def get_publication_experiment_trace(self, **kwargs: str) -> dict[str, str]:
            self.traces.append(kwargs)
            if kwargs["artifact_generation"].endswith("b" * 64):
                return {
                    "trace_status": "unresolved",
                    "reason": "experiment_snapshot_binding_absent",
                }
            return {
                "trace_status": "experiment_snapshot",
                "result": {"artifact_revision": "artifactrev-" + "a" * 64},
                "project_id": kwargs["project_id"],
                "experiment_revision_id": REVISION,
                "run_id": kwargs["run_id"],
                "artifact_id": kwargs["artifact_id"],
                "artifact_generation": kwargs["artifact_generation"],
            }

    catalog = _Catalog()
    gate = CatalogExperimentGate(catalog, lambda: _Rows(rows))
    generation = "artifactgen-" + "a" * 64
    gate.require_run_publication(
        project_id=PROJECT,
        experiment_id=EXPERIMENT,
        experiment_revision_id=REVISION,
        run_id="run-1",
        artifact_id="signal",
        generation=generation,
        revision="artifactrev-" + "a" * 64,
    )
    assert catalog.traces[0]["artifact_generation"] == generation

    with pytest.raises(ArchiveSealError) as missing:
        gate.require_selected_revision(
            project_id=PROJECT,
            experiment_id=EXPERIMENT,
            experiment_revision_id=OTHER_REVISION,
        )
    assert missing.value.code == "archive_version_not_found"

    with pytest.raises(ArchiveSealError) as unbound:
        gate.require_run_publication(
            project_id=PROJECT,
            experiment_id=EXPERIMENT,
            experiment_revision_id=REVISION,
            run_id="run-1",
            artifact_id="signal",
            generation="artifactgen-" + "b" * 64,
            revision="artifactrev-" + "a" * 64,
        )
    assert unbound.value.code == "archive_invalid"

    class _MissingCatalog:
        def get_publication_experiment_trace(self, **kwargs: str) -> dict[str, str]:
            raise catalog_not_found()

        def get_experiment_revision(
            self,
            experiment_id: str,
            experiment_revision_id: str,
        ) -> dict[str, str]:
            return {
                "experiment_id": experiment_id,
                "experiment_revision_id": experiment_revision_id,
            }

    missing_gate = CatalogExperimentGate(_MissingCatalog(), lambda: _Rows(rows))
    with pytest.raises(ArchiveSealError) as not_found:
        missing_gate.require_run_publication(
            project_id=PROJECT,
            experiment_id=EXPERIMENT,
            experiment_revision_id=REVISION,
            run_id="run-1",
            artifact_id="signal",
            generation=generation,
            revision="artifactrev-" + "a" * 64,
        )
    assert not_found.value.code == "archive_version_not_found"


def _runtime(tmp_path: Path) -> tuple[ArtifactDownloadService, RunService, Path]:
    run_service = RunService(
        create_default_workflow_registry(),
        id_factory=lambda: "run-1",
    )
    record = run_service.create_run(WORKFLOW_ID, WorkflowInputs(config={}))
    assert record.run_id == "run-1"
    for status in (
        RunStatus.VALIDATING,
        RunStatus.PLANNED,
        RunStatus.QUEUED,
        RunStatus.RUNNING,
        RunStatus.SUCCEEDED,
    ):
        run_service.transition_run("run-1", status)
    workspace = tmp_path / "workspaces"
    workspace.mkdir()
    return (
        ArtifactDownloadService(run_service=run_service, workspace_root=workspace),
        run_service,
        workspace,
    )


def _artifact(relative_path: str, content: bytes, artifact_id: str) -> RunArtifactRef:
    return RunArtifactRef(
        artifact_id=artifact_id,
        run_id="run-1",
        artifact_type="file",
        name=PurePosixPath(relative_path).name,
        uri=f"run://runs/run-1/artifacts/{artifact_id}",
        mime_type="application/octet-stream",
        produced_at=datetime.now(timezone.utc),
        revision="artifactrev-" + "0" * 64,
        metadata={
            "relative_path": relative_path,
            "output_type": "file",
            "size_bytes": len(content),
        },
    )


def _bind(artifact: RunArtifactRef, workspace: Path, content: bytes) -> RunArtifactRef:
    relative_path = str(artifact.metadata["relative_path"])
    target = workspace / artifact.run_id / relative_path
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(content)
    components = (
        *(workspace / artifact.run_id).parts[1:],
        *PurePosixPath(relative_path).parts,
    )
    current = Path("/")
    infos = []
    for component in components:
        current /= component
        infos.append(os.stat(current, follow_symlinks=False))
    file_info = infos[-1]
    path_identity = build_artifact_path_identity(
        parent_identities=tuple(
            (info.st_dev, info.st_ino, info.st_mode) for info in infos[:-1]
        ),
        file_identity=(
            file_info.st_dev,
            file_info.st_ino,
            file_info.st_mode,
            file_info.st_nlink,
            file_info.st_uid,
            file_info.st_gid,
            file_info.st_size,
            file_info.st_mtime_ns,
            file_info.st_ctime_ns,
        ),
    )
    revision = build_artifact_content_revision(
        output_type=str(artifact.metadata["output_type"]),
        relative_path=relative_path,
        content=content,
    )
    return replace(
        artifact,
        revision=revision,
        metadata={
            **artifact.metadata,
            ARTIFACT_PATH_IDENTITY_METADATA_KEY: path_identity,
        },
    )


def _publish_run(
    run_service: RunService,
    workspace: Path,
    files: tuple[tuple[str, str, bytes], ...],
) -> tuple[str, dict[str, RunArtifactRef]]:
    artifacts = []
    for artifact_id, relative_path, content in files:
        artifacts.append(
            _bind(_artifact(relative_path, content, artifact_id), workspace, content)
        )
    run_service.replace_artifacts("run-1", tuple(artifacts))
    generation = run_service.get_result_state("run-1").artifact_generation
    assert generation is not None
    stored = {
        artifact.artifact_id: artifact
        for artifact in run_service.list_artifacts("run-1")
    }
    return generation, stored


def test_run_seal_reads_the_artifact_descriptor_not_a_filename(tmp_path: Path) -> None:
    downloads, runs, workspace = _runtime(tmp_path)
    content = b"run-signal"
    (workspace / "run-1").mkdir()
    (workspace / "run-1" / "signal.bw").write_bytes(b"decoy-bw")
    generation, stored = _publish_run(
        runs,
        workspace,
        (("signal", "results/tracked-signal.bw", content),),
    )
    archive = tmp_path / "archive"
    archive.mkdir()
    objects = _Objects()
    repository = _Repository()
    gate = _Allow()
    service = ArchiveSealingService(
        archive_root=archive,
        objects=objects,
        repository=repository,
        experiments=gate,
        downloads=downloads,
    )

    receipt = service.seal_run(
        project_id=PROJECT,
        experiment_id=EXPERIMENT,
        experiment_revision_id=REVISION,
        track_key="signal",
        display_name="signal.bw",
        assembly="hg38",
        track_type="bigwig",
        run_id="run-1",
        artifact_id="signal",
        generation=generation,
        revision=stored["signal"].revision,
    )

    assert objects.payloads == [content]
    assert receipt["source"]["generation"] == generation
    assert receipt["assembly"] == "hg38"
    assert receipt["source"]["revision"] == stored["signal"].revision
    assert gate.publications[0]["generation"] == generation
    assert (
        workspace / "run-1" / "results" / "tracked-signal.bw"
    ).read_bytes() == content
    with pytest.raises(ArchiveSealError) as captured:
        service.seal_run(
            project_id=PROJECT,
            experiment_id=EXPERIMENT,
            experiment_revision_id=REVISION,
            track_key="signal",
            display_name="signal.bw",
            assembly="hg38",
            track_type="bed",
            run_id="run-1",
            artifact_id="signal",
            generation=generation,
            revision=stored["signal"].revision,
        )
    assert captured.value.code == "archive_invalid"
    assert objects.payloads == [content]


def test_run_assembly_conflict_does_not_open_the_artifact(tmp_path: Path) -> None:
    downloads, runs, workspace = _runtime(tmp_path)
    content = b"run-signal"
    generation, stored = _publish_run(
        runs,
        workspace,
        (("signal", "results/tracked-signal.bw", content),),
    )
    opened = {"count": 0}

    class _Downloads:
        def prepare(self, *args: object, **kwargs: object) -> object:
            opened["count"] += 1
            return downloads.prepare(*args, **kwargs)  # type: ignore[arg-type]

    class _Bound:
        def __init__(self) -> None:
            self.assembly = "mm10"

        def run_assembly(self, run_id: str) -> str:
            assert run_id == "run-1"
            return self.assembly

    bound = _Bound()
    archive = tmp_path / "archive"
    archive.mkdir()
    objects = _Objects()
    service = ArchiveSealingService(
        archive_root=archive,
        objects=objects,
        repository=_Repository(),
        experiments=_Allow(),
        downloads=_Downloads(),
        assemblies=bound,
    )
    arguments = {
        "project_id": PROJECT,
        "experiment_id": EXPERIMENT,
        "experiment_revision_id": REVISION,
        "track_key": "signal",
        "display_name": "signal.bw",
        "assembly": "hg38",
        "track_type": "bigwig",
        "run_id": "run-1",
        "artifact_id": "signal",
        "generation": generation,
        "revision": stored["signal"].revision,
    }
    with pytest.raises(ArchiveSealError) as captured:
        service.seal_run(**arguments)
    assert captured.value.code == "archive_invalid"
    assert opened["count"] == 0
    assert objects.payloads == []

    bound.assembly = "GRCh38"
    receipt = service.seal_run(**arguments)
    assert receipt["assembly"] == "hg38"
    assert objects.payloads == [content]
    assert opened["count"] == 1


def test_run_index_uses_its_own_artifact(tmp_path: Path) -> None:
    downloads, runs, workspace = _runtime(tmp_path)
    indexed = write_indexed_bed(("chr1\t10\t20\tchr1:40-50,1",))
    data, index = indexed.bgzf, indexed.tabix
    generation, stored = _publish_run(
        runs,
        workspace,
        (
            ("pairs", "results/pairs.pairs", data),
            ("pairsIndex", "results/pairs.pairs.tbi", index),
        ),
    )
    archive = tmp_path / "archive"
    archive.mkdir()
    objects = _Objects()
    service = ArchiveSealingService(
        archive_root=archive,
        objects=objects,
        repository=_Repository(),
        experiments=_Allow(),
        downloads=downloads,
    )

    receipt = service.seal_run(
        project_id=PROJECT,
        experiment_id=EXPERIMENT,
        experiment_revision_id=REVISION,
        track_key="pairs",
        display_name="pairs.bed",
        assembly="hg38",
        track_type="longrange",
        run_id="run-1",
        artifact_id="pairs",
        generation=generation,
        revision=stored["pairs"].revision,
        index_artifact_id="pairsIndex",
        index_generation=generation,
        index_revision=stored["pairsIndex"].revision,
    )

    assert objects.payloads == [data, index]
    assert receipt["has_index"] is True
    assert receipt["source"]["index_artifact_id"] == "pairsIndex"
    assert receipt["source"]["index_generation"] == generation


def test_replaced_run_bytes_are_not_sealed(tmp_path: Path) -> None:
    downloads, runs, workspace = _runtime(tmp_path)
    content = b"run-signal"
    generation, stored = _publish_run(
        runs,
        workspace,
        (("signal", "results/tracked-signal.bw", content),),
    )
    target = workspace / "run-1" / "results" / "tracked-signal.bw"
    target.write_bytes(b"run-OTHER!")
    archive = tmp_path / "archive"
    archive.mkdir()
    objects = _Objects()
    repository = _Repository()
    service = ArchiveSealingService(
        archive_root=archive,
        objects=objects,
        repository=repository,
        experiments=_Allow(),
        downloads=downloads,
    )

    with pytest.raises(ArchiveSealError) as captured:
        service.seal_run(
            project_id=PROJECT,
            experiment_id=EXPERIMENT,
            experiment_revision_id=REVISION,
            track_key="signal",
            display_name="signal.bw",
            assembly="hg38",
            track_type="bigwig",
            run_id="run-1",
            artifact_id="signal",
            generation=generation,
            revision=stored["signal"].revision,
        )

    assert captured.value.code == "archive_integrity"
    assert objects.payloads == []
    assert repository.calls == []


def test_wrong_generation_does_not_seal_another_artifact(tmp_path: Path) -> None:
    downloads, runs, workspace = _runtime(tmp_path)
    generation, stored = _publish_run(
        runs,
        workspace,
        (("signal", "results/tracked-signal.bw", b"run-signal"),),
    )
    archive = tmp_path / "archive"
    archive.mkdir()
    objects = _Objects()
    service = ArchiveSealingService(
        archive_root=archive,
        objects=objects,
        repository=_Repository(),
        experiments=_Allow(),
        downloads=downloads,
    )
    del generation

    with pytest.raises(ArchiveSealError) as captured:
        service.seal_run(
            project_id=PROJECT,
            experiment_id=EXPERIMENT,
            experiment_revision_id=REVISION,
            track_key="signal",
            display_name="signal.bw",
            assembly="hg38",
            track_type="bigwig",
            run_id="run-1",
            artifact_id="signal",
            generation="artifactgen-" + "b" * 64,
            revision=stored["signal"].revision,
        )

    assert captured.value.code == "archive_version_conflict"
    assert objects.payloads == []


def _seed_bound_run(app):
    from conftest import seed_test_authentication
    from encode_pipeline.persistence.database import create_session_factory
    from encode_pipeline.persistence.models import (
        RunExperimentBindingRow,
        RunProjectBindingRow,
        RunResultStateRow,
        RunRow,
    )
    from encode_pipeline.services.data_registry import DataRegistryService

    persistence = app.state.persistence
    actor = seed_test_authentication(app)
    registry = DataRegistryService(repository=persistence.data_registry_repository)
    project = registry.create_project("Switched publication").project_id
    sample = registry.import_sample(
        project, stable_key="control", display_name="Control", attributes={}
    )
    unknown = {"status": "unknown"}
    experiment = app.state.experiment_catalog_service.create_experiment(
        actor_user_id=actor.user_id,
        idempotency_key="archive-experiment",
        raw_body=b"{}",
        project_id=project,
        stable_key="range",
        display_name="Range",
        assay=unknown,
        target=unknown,
        batch=unknown,
        samples=[
            {
                "sample_revision_id": sample.revision.sample_revision_id,
                "replicate_role": unknown,
            }
        ],
    )
    revision_id = experiment["experiment_revision_id"]
    now = datetime.now(timezone.utc)
    sessions = create_session_factory(persistence.engine)
    with sessions.begin() as db:
        db.add(
            RunRow(
                run_id="run-1",
                workflow_id=WORKFLOW_ID,
                inputs=WorkflowInputs(config={}).to_dict(),
                tags={},
                status="succeeded",
                created_at=now,
                updated_at=now,
                ended_at=now,
            )
        )
        db.flush()
        db.add_all(
            [
                RunResultStateRow(run_id="run-1"),
                RunProjectBindingRow(
                    run_id="run-1",
                    project_id=project,
                    binding_mode="bound_v1",
                    provenance="resolved",
                    workflow_inputs_digest="0" * 64,
                    binding_digest_scheme="sha256-framed-project-sample-binding-v1",
                    binding_digest="0" * 64,
                    created_at=now,
                ),
                RunExperimentBindingRow(
                    run_id="run-1",
                    experiment_revision_id=revision_id,
                    payload_digest=db.get(
                        ExperimentRevisionRow, revision_id
                    ).payload_digest,
                ),
            ]
        )
    return project, experiment["experiment_id"], revision_id


def test_installed_cli_seals_switched_generation_after_reindex(tmp_path: Path):
    import asyncio
    import subprocess
    import sys

    import httpx
    from encode_pipeline.api.main import create_app
    from encode_pipeline.services.artifact_disk_relocation import ArtifactDiskRelocation

    workspace, durable, archive = (
        tmp_path / name for name in ("workspaces", "relocation", "archive")
    )
    archive.mkdir()
    database = f"sqlite:///{tmp_path / 'platform.db'}"
    app = create_app(
        database_url=database, workspace_root=workspace, relocation_root=durable
    )
    persistence = app.state.persistence
    try:
        project, experiment, revision_id = _seed_bound_run(app)
        path = "results/peaks.bed"
        old, new = b"chr1\t10\t20\tA\n", b"chr1\t10\t20\tB\n"
        runs = app.state.run_service
        first_gen, first = _publish_run(runs, workspace, (("peaks", path, old),))
        old_rev = first["peaks"].revision
        relocation = ArtifactDiskRelocation(source_root=workspace, durable_root=durable)
        receipt = app.state.artifact_download_service.capture_relocation(
            relocation,
            run_id="run-1",
            artifact_id="peaks",
            expected_generation=first_gen,
            expected_revision=old_rev,
        )
        for _ in range(2):
            receipt = relocation.resume(run_id="run-1", artifact_id="peaks")
        assert receipt.stage == "verified"
        second_gen, second = _publish_run(runs, workspace, (("peaks", path, new),))
        assert second_gen != first_gen

        def seal(generation, revision, *extra, artifact="peaks", track_type="bed"):
            command = [
                str(Path(sys.executable).with_name("helixweave")),
                "archive",
                "seal",
                "run",
                "--database-url",
                database,
                "--archive-root",
                str(archive),
                "--workspace-root",
                str(workspace),
                "--project-id",
                project,
                "--experiment-id",
                experiment,
                "--experiment-revision-id",
                revision_id,
                "--track-key",
                artifact,
                "--display-name",
                artifact,
                "--assembly",
                "GRCh38",
                "--track-type",
                track_type,
                "--run-id",
                "run-1",
                "--artifact-id",
                artifact,
                "--generation",
                generation,
                "--revision",
                revision,
                *extra,
            ]
            env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
            proc = subprocess.run(
                command,
                cwd=tmp_path,
                env=env,
                capture_output=True,
                text=True,
                check=False,
            )
            return proc.returncode, json.loads(
                proc.stderr if proc.returncode else proc.stdout
            )

        async def download_old():
            async with httpx.AsyncClient(
                transport=httpx.ASGITransport(app=app),
                base_url="http://test",
                cookies={
                    app.state.auth_cookie_policy.session_cookie.name: app.state.test_auth_tokens[
                        0
                    ]
                },
            ) as client:
                return await client.get(
                    "/api/v1/runs/run-1/artifacts/peaks/download",
                    params={"generation": first_gen, "revision": old_rev},
                )

        assert (
            seal(first_gen, old_rev)[1]["error"]["code"] == "archive_version_conflict"
        )
        option = ("--relocation-root", str(durable))
        assert (
            seal(first_gen, old_rev, *option)[1]["error"]["code"]
            == "archive_version_conflict"
        )
        assert (
            relocation.resume(run_id="run-1", artifact_id="peaks").stage == "switched"
        )
        code, v1 = seal(first_gen, old_rev, *option)
        assert code == 0, v1
        code, v2 = seal(second_gen, second["peaks"].revision, *option)
        assert code == 0, v2
        assert [v1["version_number"], v2["version_number"]] == [1, 2]
        assert [v1["source"]["revision"], v2["source"]["revision"]] == [
            old_rev,
            second["peaks"].revision,
        ]
        assert [
            (archive / v["members"]["data"]["sha256"]).read_bytes() for v in (v1, v2)
        ] == [old, new]
        assert (
            seal(first_gen, old_rev, *option)[1]["archive_version_id"]
            == v1["archive_version_id"]
        )
        public = asyncio.run(download_old())
        assert public.status_code == 409
        assert (
            public.json()["issues"][0]["context"]["reason_code"]
            == "ARTIFACT_DOWNLOAD_GENERATION_CHANGED"
        )
        relocation.payload_file("run-1", "peaks").write_bytes(new)
        assert (
            seal(first_gen, old_rev, *option)[1]["error"]["code"] == "archive_integrity"
        )
        assert (
            seal(first_gen, second["peaks"].revision, *option)[1]["error"]["code"]
            == "archive_version_conflict"
        )

        pair = write_indexed_bed(("chr1\t10\t20\tchr1:40-50,1",))
        other = write_indexed_bed(("chr2\t10\t20\tchr2:40-50,1",))
        generation, artifacts = _publish_run(
            runs,
            workspace,
            (
                ("links", "results/links.longrange.gz", pair.bgzf),
                ("index", "results/links.longrange.gz.tbi", pair.tabix),
                ("other", "results/other/links.longrange.gz.tbi", other.tabix),
            ),
        )

        def pair_with(index):
            return seal(
                generation,
                artifacts["links"].revision,
                "--index-artifact-id",
                index,
                "--index-generation",
                generation,
                "--index-revision",
                artifacts[index].revision,
                artifact="links",
                track_type="longrange",
            )

        assert pair_with("other")[1]["error"]["code"] == "archive_invalid"
        code, indexed = pair_with("index")
        assert code == 0, indexed
        assert indexed["version_number"] == 1 and indexed["has_index"] is True
        assert (
            indexed["members"]["index"]["sha256"]
            == hashlib.sha256(pair.tabix).hexdigest()
        )
    finally:
        app.state.run_queue.close()
        persistence.close()
