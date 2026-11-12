"""HTTP contracts for WashU track identity, ranges, and display files."""

from __future__ import annotations

from datetime import datetime, timezone
import asyncio
import os
from pathlib import Path
from pathlib import PurePosixPath

import httpx
import pytest

from encode_pipeline.api.main import create_app
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.reference_profiles import (
    AdapterReferenceBindingIdentity,
    ReferenceProfile,
    build_reference_profile_revision,
    build_reference_profile_revision_binding,
)
from encode_pipeline.platform.reference_profiles import (
    ReferenceProfileWorkflowBinding,
)
from encode_pipeline.platform.result_generations import (
    ARTIFACT_PATH_IDENTITY_METADATA_KEY,
    build_artifact_descriptor_revision,
    build_artifact_path_identity,
)
from encode_pipeline.platform.runs import RunArtifactRef, RunStatus
from encode_pipeline.persistence.models import RunReferenceBindingRow
from encode_pipeline.services.artifact_downloads import _descriptor_identity_text
from api_test_client import ApiTestClient
from conftest import seed_test_authentication


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
FIXTURE = Path(__file__).parents[1] / "fixtures" / "browser" / "chr1-signal.bw"


@pytest.fixture
def client(tmp_path: Path):
    app = create_app(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        workspace_root=tmp_path / "workspaces",
    )
    seed_test_authentication(app)
    with ApiTestClient(app) as test_client:
        yield test_client


def _create_run(client: ApiTestClient) -> str:
    service = client.app.state.run_service
    run_id = service.create_run(WORKFLOW_ID, WorkflowInputs(config={})).run_id
    for status in (
        RunStatus.VALIDATING,
        RunStatus.PLANNED,
        RunStatus.QUEUED,
        RunStatus.RUNNING,
        RunStatus.SUCCEEDED,
    ):
        service.transition_run(run_id, status)
    return run_id


def _record(
    client: ApiTestClient,
    run_id: str,
    artifact_id: str,
    content: bytes,
    *,
    name: str,
    output_type: str,
    revision_kind: str = "descriptor",
    reference_assembly: str | None = None,
) -> RunArtifactRef:
    relative_path = f"results/browser/{name}"
    source = client.app.state.workspace_root / run_id / relative_path
    source.parent.mkdir(parents=True, exist_ok=True)
    source.write_bytes(content)
    components = (
        *(client.app.state.workspace_root / run_id).parts[1:],
        *PurePosixPath(relative_path).parts,
    )
    current = Path("/")
    infos = []
    for component in components:
        current /= component
        infos.append(current.stat())
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
    if revision_kind == "descriptor":
        revision = build_artifact_descriptor_revision(
            output_type=output_type,
            relative_path=relative_path,
            descriptor_identity=_descriptor_identity_text(file_info),
        )
    else:
        from encode_pipeline.platform.result_generations import (
            build_artifact_content_revision,
        )

        revision = build_artifact_content_revision(
            output_type=output_type,
            relative_path=relative_path,
            content=content,
        )
    metadata = {
        "relative_path": relative_path,
        "output_type": output_type,
        "size_bytes": len(content),
        "catalog_id": output_type,
        "scope": "project",
        ARTIFACT_PATH_IDENTITY_METADATA_KEY: path_identity,
    }
    if reference_assembly is not None:
        metadata["reference_assembly"] = reference_assembly
    artifact = RunArtifactRef(
        artifact_id=artifact_id,
        run_id=run_id,
        artifact_type="file",
        name=name,
        uri=f"run://runs/{run_id}/artifacts/{artifact_id}",
        mime_type="application/octet-stream",
        produced_at=datetime.now(timezone.utc),
        revision=revision,
        metadata=metadata,
    )
    service = client.app.state.run_service
    existing = tuple(
        item
        for item in service.list_artifacts(run_id)
        if item.artifact_id != artifact_id
    )
    service.replace_artifacts(run_id, (*existing, artifact))
    return artifact


def _bind_reference(client: ApiTestClient, run_id: str, *, assembly: str) -> None:
    now = datetime.now(timezone.utc)
    profile = ReferenceProfile(
        profile_id="refp_" + "ab" * 16,
        safe_key="browser-fixture",
        created_at=now,
    )
    identity = AdapterReferenceBindingIdentity(
        workflow_id=WORKFLOW_ID,
        contract_version="browser-fixture-v1",
        identity_sha256="c" * 64,
    )
    revision = build_reference_profile_revision(
        revision_id="refpr_" + "cd" * 16,
        profile_id=profile.profile_id,
        revision_number=1,
        display_name="Browser fixture",
        organism="Homo sapiens",
        assembly=assembly,
        config_key="browser-fixture",
        workflow_bindings=(
            ReferenceProfileWorkflowBinding.from_adapter_identity(identity),
        ),
        created_at=now,
    )
    client.app.state.reference_profile_service.repository.create_profile(
        profile,
        revision,
    )
    binding = build_reference_profile_revision_binding(
        profile_id=profile.profile_id,
        revision_id=revision.revision_id,
        workflow_id=WORKFLOW_ID,
        revision_public_identity_sha256=revision.public_identity_sha256,
        adapter_identity=identity,
        bound_at=now,
    )
    repository = client.app.state.run_service._repository
    with repository._session_factory.begin() as session:
        session.add(
            RunReferenceBindingRow(
                run_id=run_id,
                profile_id=binding.profile_id,
                revision_id=binding.revision_id,
                workflow_id=binding.workflow_id,
                revision_public_identity_scheme=(
                    binding.revision_public_identity_scheme
                ),
                revision_public_identity_sha256=(
                    binding.revision_public_identity_sha256
                ),
                adapter_contract_version=binding.adapter_contract_version,
                adapter_identity_scheme=binding.adapter_identity_scheme,
                adapter_identity_sha256=binding.adapter_identity_sha256,
                binding_digest_scheme=binding.binding_digest_scheme,
                binding_digest=binding.binding_digest,
                bound_at=binding.bound_at,
            )
        )


def _generation(client: ApiTestClient, run_id: str) -> str:
    return client.app.state.run_service.get_result_state(run_id).artifact_generation


def _byte_url(run_id: str, artifact_id: str, generation: str, revision: str) -> str:
    return (
        f"/api/v1/runs/{run_id}/browser-tracks/{artifact_id}"
        f"/bytes/{generation}/{revision}"
    )


def test_download_route_still_refuses_ranges(client: ApiTestClient) -> None:
    run_id = _create_run(client)
    artifact = _record(
        client,
        run_id,
        "artifact-signal",
        FIXTURE.read_bytes(),
        name="chr1-signal.bw",
        output_type="cpm_bigwig",
    )
    generation = _generation(client, run_id)
    response = client.get(
        f"/api/v1/runs/{run_id}/artifacts/{artifact.artifact_id}/download",
        params={"generation": generation, "revision": artifact.revision},
        headers={"Range": "bytes=0-3"},
    )
    assert response.status_code == 200
    assert response.headers["accept-ranges"] == "none"
    assert "content-range" not in response.headers


def test_original_bigwig_range_reads_only_the_requested_bytes(
    client: ApiTestClient,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    run_id = _create_run(client)
    content = FIXTURE.read_bytes()
    artifact = _record(
        client,
        run_id,
        "artifact-signal",
        content,
        name="chr1-signal.bw",
        output_type="cpm_bigwig",
    )
    generation = _generation(client, run_id)
    inode = (
        (client.app.state.workspace_root / run_id / "results/browser/chr1-signal.bw")
        .stat()
        .st_ino
    )
    reads: list[int] = []
    real_read = os.read

    def spy(fd: int, size: int) -> bytes:
        chunk = real_read(fd, size)
        try:
            if os.fstat(fd).st_ino == inode:
                reads.append(len(chunk))
        except OSError:
            pass
        return chunk

    monkeypatch.setattr(os, "read", spy)
    response = client.get(
        _byte_url(run_id, artifact.artifact_id, generation, artifact.revision),
        headers={"Range": "bytes=4-13"},
    )
    assert response.status_code == 206
    assert response.headers["accept-ranges"] == "bytes"
    assert response.headers["content-range"] == f"bytes 4-13/{len(content)}"
    assert response.headers["content-length"] == "10"
    assert response.content == content[4:14]
    assert sum(reads) == 10
    assert str(client.app.state.workspace_root) not in response.text


def test_head_range_and_rejected_ranges_release_descriptors(
    client: ApiTestClient,
) -> None:
    run_id = _create_run(client)
    content = b"abcdefghijklmnopqrstuvwxyz"
    artifact = _record(
        client,
        run_id,
        "artifact-signal",
        content,
        name="chr1-signal.bw",
        output_type="cpm_bigwig",
    )
    generation = _generation(client, run_id)
    url = _byte_url(run_id, artifact.artifact_id, generation, artifact.revision)
    before = len(list(Path("/proc/self/fd").iterdir()))
    head = client.request("HEAD", url, headers={"Range": "bytes=0-3"})
    assert head.status_code == 206
    assert head.headers["content-range"] == f"bytes 0-3/{len(content)}"
    assert head.content == b""
    missing = client.get(url, headers={"Range": "bytes=100-120"})
    assert missing.status_code == 416
    assert missing.headers["content-range"] == f"bytes */{len(content)}"
    multi = client.get(url, headers={"Range": "bytes=0-1,2-3"})
    assert multi.status_code == 416
    full = client.get(url)
    assert full.status_code == 200
    assert full.content == content
    after = len(list(Path("/proc/self/fd").iterdir()))
    assert after <= before + 2


def test_anonymous_and_stale_identity_do_not_return_bytes(
    client: ApiTestClient,
) -> None:
    run_id = _create_run(client)
    artifact = _record(
        client,
        run_id,
        "artifact-signal",
        b"signal-bytes",
        name="chr1-signal.bw",
        output_type="cpm_bigwig",
    )
    generation = _generation(client, run_id)
    url = _byte_url(run_id, artifact.artifact_id, generation, artifact.revision)

    async def anonymous() -> httpx.Response:
        transport = httpx.ASGITransport(app=client.app)
        async with httpx.AsyncClient(
            transport=transport,
            base_url="http://testserver",
        ) as http:
            return await http.get(url)

    anon = asyncio.run(anonymous())
    assert anon.status_code == 401
    assert str(client.app.state.workspace_root) not in anon.text

    missing = client.get(
        _byte_url(run_id, "missing-artifact", generation, artifact.revision)
    )
    assert missing.status_code == 404
    assert missing.json()["issues"][0]["code"] == "RUN_ARTIFACT_NOT_FOUND"

    stale = client.get(
        _byte_url(
            run_id,
            artifact.artifact_id,
            generation,
            "artifactrev-" + "d" * 64,
        )
    )
    assert stale.status_code == 409
    assert stale.json()["issues"][0]["context"]["reason_code"] == "revision_changed"
    changed = client.get(
        _byte_url(
            run_id,
            artifact.artifact_id,
            "artifactgen-" + "e" * 64,
            artifact.revision,
        )
    )
    assert changed.status_code == 409
    assert "generation" in changed.text
    assert str(client.app.state.workspace_root) not in changed.text


def test_pet_derivative_and_signal_share_a_reference_window(
    client: ApiTestClient,
) -> None:
    run_id = _create_run(client)
    _bind_reference(client, run_id, assembly="hg38")
    signal = _record(
        client,
        run_id,
        "artifact-signal",
        FIXTURE.read_bytes(),
        name="chr1-signal.bw",
        output_type="cpm_bigwig",
        reference_assembly="hg38",
    )
    pet = _record(
        client,
        run_id,
        "artifact-pet",
        b"chr1\t1200\t1400\tchr1\t1600\t1800\tread\t211\t+\t+\n",
        name="sample_all.bedpe.gz",
        output_type="hitrac_bedpe_all",
        revision_kind="content",
        reference_assembly="hg38",
    )
    _record(
        client,
        run_id,
        "artifact-cool",
        b"cooler",
        name="contacts.mcool",
        output_type="mcool",
    )
    generation = _generation(client, run_id)
    session = client.get(f"/api/v1/runs/{run_id}/browser-session")
    assert session.status_code == 200
    body = session.json()
    assert body["reference"]["status"] == "bound"
    assert body["reference"]["viewer_scaffold"] == "hg38"
    assert body["reference"]["assembly_label"] == "hg38"
    assert len(body["reference"]["identity_sha256"]) == 64
    assert body["interaction_status"] == "available"
    assert body["unavailable"][0]["reason"] == "unsupported_contact_matrix"
    prepared = client.post(
        f"/api/v1/runs/{run_id}/browser-tracks/artifact-pet/derivative",
        json={"generation": generation, "revision": pet.revision},
    )
    assert prepared.status_code == 200, prepared.text
    track = prepared.json()["track"]
    assert track["kind"] == "pet_records"
    assert track["provenance"]["pair_count"] == 1
    assert track["provenance"]["anchor_row_count"] == 2
    assert track["preparation"] == "ready"
    summary = client.get(
        f"/api/v1/runs/{run_id}/browser-tracks/artifact-pet/region-summary/"
        f"{generation}/{pet.revision}",
        params={"chrom": "chr1", "start": 1000, "end": 2000},
    )
    assert summary.status_code == 200
    assert summary.json()["records"] == 2
    empty = client.get(
        f"/api/v1/runs/{run_id}/browser-tracks/artifact-pet/region-summary/"
        f"{generation}/{pet.revision}",
        params={"chrom": "chr1", "start": 8000, "end": 9000},
    )
    assert empty.json()["records"] == 0
    assert "not a zero-signal" in empty.json()["detail"]
    signal_summary = client.get(
        f"/api/v1/runs/{run_id}/browser-tracks/{signal.artifact_id}/region-summary/"
        f"{generation}/{signal.revision}",
        params={"chrom": "chr1", "start": 1000, "end": 2000},
    )
    assert signal_summary.json()["records"] is None
    assert signal_summary.json()["kind"] == "signal"
    derivative = client.get(track["content_path"], headers={"Range": "bytes=0-3"})
    assert derivative.status_code == 206
    assert str(client.app.state.workspace_root) not in derivative.text
    index = client.get(track["index_path"], headers={"Range": "bytes=0-3"})
    assert index.status_code == 206
    assert index.content == b"\x1f\x8b\x08\x04"


def test_replaced_display_index_cannot_be_read(client: ApiTestClient) -> None:
    run_id = _create_run(client)
    pet = _record(
        client,
        run_id,
        "artifact-pet",
        b"chr1\t1200\t1400\tchr1\t1600\t1800\tread\t211\t+\t+\n",
        name="sample_all.bedpe.gz",
        output_type="hitrac_bedpe_all",
        revision_kind="content",
    )
    generation = _generation(client, run_id)
    prepared = client.post(
        f"/api/v1/runs/{run_id}/browser-tracks/artifact-pet/derivative",
        json={"generation": generation, "revision": pet.revision},
    )
    track = prepared.json()["track"]
    index_file = (
        client.app.state.workspace_root
        / "browser-derivatives"
        / run_id
        / pet.artifact_id
        / pet.revision
        / "track.gz.tbi"
    )
    index_file.write_bytes(index_file.read_bytes() + b"tamper")
    response = client.get(track["index_path"], headers={"Range": "bytes=0-3"})
    assert response.status_code == 409
    assert (
        response.json()["issues"][0]["context"]["reason_code"] == "derivative_changed"
    )
    assert str(client.app.state.workspace_root) not in response.text
    assert "track.gz" not in response.text
