"""Same-origin platform page with two hg38 runs and one mm10 run.

The browser proof selects local files itself. This server only supplies the
authorized runs those files are overlaid on, plus a control that changes one
run's artifact revision.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path
from pathlib import PurePosixPath
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "test" / "browser"))

import uvicorn  # noqa: E402
from encode_pipeline.api.main import create_app  # noqa: E402
from encode_pipeline.deployment.frontend import ManifestStaticApplication  # noqa: E402
from encode_pipeline.frontend_assets import load_packaged_frontend_assets  # noqa: E402
from encode_pipeline.platform.adapters import WorkflowInputs  # noqa: E402
from encode_pipeline.platform.reference_profiles import (  # noqa: E402
    AdapterReferenceBindingIdentity,
    ReferenceProfile,
    ReferenceProfileWorkflowBinding,
    build_reference_profile_revision,
    build_reference_profile_revision_binding,
)
from encode_pipeline.platform.result_generations import (  # noqa: E402
    ARTIFACT_PATH_IDENTITY_METADATA_KEY,
    build_artifact_content_revision,
    build_artifact_descriptor_revision,
    build_artifact_path_identity,
)
from encode_pipeline.platform.runs import RunArtifactRef, RunStatus  # noqa: E402
from encode_pipeline.persistence.models import RunReferenceBindingRow  # noqa: E402
from encode_pipeline.services.artifact_downloads import (  # noqa: E402
    _descriptor_identity_text,
)
from encode_pipeline.services.browser_tracks import BrowserTrackService  # noqa: E402
from local_overlay_fixtures import load_consensus, write_signal_bigwig  # noqa: E402

WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


# Same fixed session as test/conftest.py. Inlined so this server can boot in a
# worktree whose shared environment is editable-installed from another checkout.
_AUTH_PASSWORD_HASH = (
    "$argon2id$v=19$m=65536,t=3,p=4$"
    "pGLQVK/BOQLC0oPSA8RQTg$TdiQWUP9gwBAI8iAXUT7oEtjOPqKjJhqyW0JS8ye/Ag"
)
_AUTH_SESSION_TOKEN = "a" * 43
_AUTH_CSRF_TOKEN = "b" * 43
_AUTH_USER_ID = "usr_" + "0" * 32


def seed_test_authentication(app) -> None:
    from datetime import datetime, timedelta, timezone

    from encode_pipeline.platform.authentication import (
        UserAccount,
        UserRole,
        UserStatus,
    )
    from encode_pipeline.services.authentication import (
        SessionSecrets,
        new_session_record,
    )
    from encode_pipeline.services.authentication_repositories import (
        AuthenticationAccountConflictError,
    )

    repository = app.state.authentication_repository
    now = datetime.now(timezone.utc)
    account = UserAccount(
        user_id=_AUTH_USER_ID,
        username="test-admin",
        role=UserRole.ADMINISTRATOR,
        status=UserStatus.ENABLED,
        password_hash=_AUTH_PASSWORD_HASH,
        created_at=now,
        updated_at=now,
        password_changed_at=now,
    )
    try:
        repository.create_account(account)
    except AuthenticationAccountConflictError:
        pass
    secrets = SessionSecrets._from_generated(_AUTH_SESSION_TOKEN, _AUTH_CSRF_TOKEN)
    record = new_session_record(
        user_id=account.user_id,
        secrets=secrets,
        created_at=now,
        lifetime=timedelta(hours=8),
    )
    try:
        repository.get_session(record.session_digest)
    except KeyError:
        repository.create_session(record)
    app.state.test_auth_tokens = (_AUTH_SESSION_TOKEN, _AUTH_CSRF_TOKEN)


class _Client:
    def __init__(self, app) -> None:
        self.app = app


def _create_run(client: _Client) -> str:
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
    client: _Client,
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


def _generation(client: _Client, run_id: str) -> str:
    return client.app.state.run_service.get_result_state(run_id).artifact_generation


def _bind_reference(client: _Client, run_id: str, *, assembly: str, token: str) -> None:
    now = datetime.now(timezone.utc)
    profile = ReferenceProfile(
        profile_id="refp_" + token,
        safe_key=f"overlay-{assembly}-{token[:8]}",
        created_at=now,
    )
    identity = AdapterReferenceBindingIdentity(
        workflow_id=WORKFLOW_ID,
        contract_version="browser-fixture-v1",
        identity_sha256=token + token,
    )
    revision = build_reference_profile_revision(
        revision_id="refpr_" + token,
        profile_id=profile.profile_id,
        revision_number=1,
        display_name=f"Overlay {assembly}",
        organism="Homo sapiens" if assembly == "hg38" else "Mus musculus",
        assembly=assembly,
        config_key=f"overlay-{assembly}-{token[:8]}",
        workflow_bindings=(
            ReferenceProfileWorkflowBinding.from_adapter_identity(identity),
        ),
        created_at=now,
    )
    client.app.state.reference_profile_service.repository.create_profile(
        profile, revision
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
                revision_public_identity_scheme=binding.revision_public_identity_scheme,
                revision_public_identity_sha256=binding.revision_public_identity_sha256,
                adapter_contract_version=binding.adapter_contract_version,
                adapter_identity_scheme=binding.adapter_identity_scheme,
                adapter_identity_sha256=binding.adapter_identity_sha256,
                binding_digest_scheme=binding.binding_digest_scheme,
                binding_digest=binding.binding_digest,
                bound_at=binding.bound_at,
            )
        )


def _signal(
    path: Path, chrom: str, size: int, start: int, end: int, value: float
) -> bytes:
    write_signal_bigwig(path, [(chrom, size)], [(chrom, start, end, value)])
    return path.read_bytes()


def build_app(database: Path, workspace: Path):
    app = create_app(database_url=f"sqlite:///{database}", workspace_root=workspace)
    seed_test_authentication(app)
    client = _Client(app)
    spec = load_consensus()
    platform = spec["platform"]
    interval = platform["run_signal_interval"]
    signal_dir = workspace / "signals"
    signal_dir.mkdir(parents=True, exist_ok=True)
    run_a = _create_run(client)
    _bind_reference(client, run_a, assembly="hg38", token="a" * 32)
    signal_a = _record(
        client,
        run_a,
        "artifact-signal-a",
        _signal(
            signal_dir / "run-a.bw",
            interval["chrom"],
            spec["chrom_sizes_hg38"]["chr7"],
            interval["start"],
            interval["end"],
            float(platform["run_signal_value"]),
        ),
        name="run-a-signal.bw",
        output_type="cpm_bigwig",
        reference_assembly="hg38",
    )
    run_b = _create_run(client)
    _bind_reference(client, run_b, assembly="hg38", token="b" * 32)
    peak_b = _record(
        client,
        run_b,
        "artifact-peak-b",
        (platform["run_peak_line"] + "\n").encode("ascii"),
        name="run-b-peaks.bed",
        output_type="macs3_peak",
        revision_kind="content",
        reference_assembly="hg38",
    )
    run_mm10 = _create_run(client)
    _bind_reference(client, run_mm10, assembly="mm10", token="c" * 32)
    signal_mm10 = _record(
        client,
        run_mm10,
        "artifact-signal-mm10",
        _signal(
            signal_dir / "run-mm10.bw",
            "chr7",
            spec["mm10_chr7"],
            1000,
            1100,
            float(platform["mm10_signal_value"]),
        ),
        name="run-mm10-signal.bw",
        output_type="cpm_bigwig",
        reference_assembly="mm10",
    )
    generation_b = _generation(client, run_b)
    service = BrowserTrackService(
        run_service=app.state.run_service,
        download_service=app.state.artifact_download_service,
        workspace_root=workspace,
        reference_profile_service=app.state.reference_profile_service,
    )
    service.prepare(
        run_b,
        peak_b.artifact_id,
        generation=generation_b,
        revision=peak_b.revision,
    )
    info = {
        "run_a": run_a,
        "run_b": run_b,
        "run_mm10": run_mm10,
        "signal_a": signal_a.artifact_id,
        "peak_b": peak_b.artifact_id,
        "signal_mm10": signal_mm10.artifact_id,
        "generation_b": generation_b,
    }
    app.state.overlay_proof = info

    @app.post("/proof-control/revise-run-b")
    def _revise() -> dict[str, str]:
        generation = revise_run_b(app)
        return {"ok": "true", "generation": generation}

    app.mount(
        "/", ManifestStaticApplication(load_packaged_frontend_assets()), name="frontend"
    )
    return app, info


def revise_run_b(app) -> str:
    """Replace run B's peak so its revision and generation move."""

    info = app.state.overlay_proof
    spec = load_consensus()
    client = _Client(app)
    updated = _record(
        client,
        info["run_b"],
        info["peak_b"],
        (spec["platform"]["run_peak_revised_line"] + "\n").encode("ascii"),
        name="run-b-peaks.bed",
        output_type="macs3_peak",
        revision_kind="content",
        reference_assembly="hg38",
    )
    generation = _generation(client, info["run_b"])
    info["generation_b"] = generation
    info["peak_revision"] = updated.revision
    return generation


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    args = parser.parse_args()
    app, info = build_app(args.database, args.workspace)
    print(
        "LOCAL_OVERLAY_SERVER "
        f"run_a={info['run_a']} run_b={info['run_b']} run_mm10={info['run_mm10']} "
        f"signal_a={info['signal_a']} peak_b={info['peak_b']} "
        f"signal_mm10={info['signal_mm10']}",
        flush=True,
    )
    uvicorn.run(app, host="127.0.0.1", port=args.port, log_level="warning")


if __name__ == "__main__":
    main()
