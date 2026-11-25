"""Serve three synthetic runs and the packaged page on one origin.

The opened run and a second run are both hg38, with different reference
profile identities. A third run is mm10. Track bytes stay on the existing
session and artifact routes. POST /__washu_overlay__/invalidate replaces the
second run's peak so its previous generation and revision no longer read.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import secrets
import sys

import uvicorn
from fastapi import Request
from fastapi.responses import JSONResponse

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "test"))
sys.path.insert(0, str(ROOT / "test" / "api"))

from encode_pipeline.api.main import create_app  # noqa: E402
from encode_pipeline.deployment.frontend import ManifestStaticApplication  # noqa: E402
from encode_pipeline.frontend_assets import load_packaged_frontend_assets  # noqa: E402
from encode_pipeline.persistence.models import RunReferenceBindingRow  # noqa: E402
from encode_pipeline.platform.reference_profiles import (  # noqa: E402
    AdapterReferenceBindingIdentity,
    ReferenceProfile,
    ReferenceProfileWorkflowBinding,
    build_reference_profile_revision,
    build_reference_profile_revision_binding,
)
from encode_pipeline.services.browser_tracks import BrowserTrackService  # noqa: E402
from conftest import seed_test_authentication  # noqa: E402
from test_routes_browser_tracks import (  # noqa: E402
    WORKFLOW_ID,
    _bind_reference,
    _create_run,
    _generation,
    _record,
)

ARTIFACT_ID = "artifact-peak"
WINDOW = "chr1:4500-6700"
PRIMARY_PEAK = b"chr1\t5500\t5700\tpeak-primary\t61\t+\n"
SECOND_PEAK = b"chr1\t6000\t6200\tpeak-second\t73\t+\n"
FOREIGN_PEAK = b"chr1\t5500\t5700\tpeak-foreign\t11\t+\n"
PRIMARY_NAME = "primary_peaks.narrowPeak"
SECOND_NAME = "second_peaks.narrowPeak"
FOREIGN_NAME = "foreign_peaks.narrowPeak"


class _Client:
    def __init__(self, app) -> None:
        self.app = app


def _bind_distinct(
    client: _Client,
    run_id: str,
    *,
    assembly: str,
    organism: str,
    marker: str,
    display_name: str,
) -> None:
    """Bind one run to its own profile. The public identity is not shared."""

    now = datetime.now(timezone.utc)
    profile = ReferenceProfile(
        profile_id="refp_" + marker,
        safe_key=f"overlay-{assembly}-{marker[:8]}",
        created_at=now,
    )
    identity = AdapterReferenceBindingIdentity(
        workflow_id=WORKFLOW_ID,
        contract_version=f"overlay-{marker[:8]}-v1",
        identity_sha256=marker[0] * 64,
    )
    revision = build_reference_profile_revision(
        revision_id="refpr_" + marker,
        profile_id=profile.profile_id,
        revision_number=1,
        display_name=display_name,
        organism=organism,
        assembly=assembly,
        config_key=f"overlay-{assembly}-{marker[:8]}",
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


def _peak(client: _Client, run_id: str, name: str, content: bytes, assembly: str):
    return _record(
        client,
        run_id,
        ARTIFACT_ID,
        content,
        name=name,
        output_type="macs3_peak",
        revision_kind="content",
        reference_assembly=assembly,
    )


def _prepare(app, client: _Client, run_id: str, artifact) -> None:
    service = BrowserTrackService(
        run_service=app.state.run_service,
        download_service=app.state.artifact_download_service,
        workspace_root=app.state.workspace_root,
        reference_profile_service=app.state.reference_profile_service,
    )
    service.prepare(
        run_id,
        artifact.artifact_id,
        generation=_generation(client, run_id),
        revision=artifact.revision,
    )


def _install_invalidate(app, client: _Client, second_run_id: str, token: str) -> None:
    bumps = {"count": 0}

    async def invalidate(request: Request):
        supplied = request.headers.get("x-washu-overlay-token")
        if supplied != token:
            return JSONResponse({"ok": False}, status_code=404)
        try:
            payload = await request.json()
        except Exception:
            payload = None
        if not isinstance(payload, dict) or payload.get("run_id") != second_run_id:
            return JSONResponse(
                {"ok": False, "error": "unknown run"},
                status_code=400,
            )
        bumps["count"] += 1
        score = 80 + bumps["count"]
        content = f"chr1\t6000\t6200\tpeak-second\t{score}\t+\n".encode()
        try:
            artifact = _peak(client, second_run_id, SECOND_NAME, content, "hg38")
            generation = _generation(client, second_run_id)
        except Exception:
            return JSONResponse(
                {"ok": False, "error": "invalidate failed"},
                status_code=500,
            )
        return {
            "ok": True,
            "run_id": second_run_id,
            "artifact_id": artifact.artifact_id,
            "generation": generation,
            "revision": artifact.revision,
        }

    app.add_api_route(
        "/__washu_overlay__/invalidate",
        invalidate,
        methods=["POST"],
    )


def build_app(database: Path, workspace: Path):
    app = create_app(
        database_url=f"sqlite:///{database}",
        workspace_root=workspace,
    )
    seed_test_authentication(app)
    client = _Client(app)
    primary_run_id = _create_run(client)
    _bind_reference(client, primary_run_id, assembly="hg38")
    primary = _peak(client, primary_run_id, PRIMARY_NAME, PRIMARY_PEAK, "hg38")
    _prepare(app, client, primary_run_id, primary)

    second_run_id = _create_run(client)
    _bind_distinct(
        client,
        second_run_id,
        assembly="hg38",
        organism="Homo sapiens",
        marker="12" * 16,
        display_name="Overlay second hg38",
    )
    second = _peak(client, second_run_id, SECOND_NAME, SECOND_PEAK, "hg38")
    _prepare(app, client, second_run_id, second)

    foreign_run_id = _create_run(client)
    _bind_distinct(
        client,
        foreign_run_id,
        assembly="mm10",
        organism="Mus musculus",
        marker="34" * 16,
        display_name="Overlay mouse mm10",
    )
    foreign = _peak(client, foreign_run_id, FOREIGN_NAME, FOREIGN_PEAK, "mm10")
    _prepare(app, client, foreign_run_id, foreign)

    token = secrets.token_hex(24)
    _install_invalidate(app, client, second_run_id, token)
    app.mount(
        "/",
        ManifestStaticApplication(load_packaged_frontend_assets()),
        name="frontend",
    )
    manifest = {
        "primary_run_id": primary_run_id,
        "second_run_id": second_run_id,
        "foreign_run_id": foreign_run_id,
        "artifact_id": ARTIFACT_ID,
        "primary_generation": _generation(client, primary_run_id),
        "primary_revision": primary.revision,
        "second_generation": _generation(client, second_run_id),
        "second_revision": second.revision,
        "window": WINDOW,
        "local_bed": "overlay-peak-local.bed",
        "local_peak": "peak-local",
        "token": token,
    }
    return app, manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    args = parser.parse_args()
    _app, manifest = build_app(args.database, args.workspace)
    print(
        "WASHU_OVERLAY " + json.dumps(manifest, separators=(",", ":")),
        flush=True,
    )
    uvicorn.run(_app, host="127.0.0.1", port=args.port, log_level="warning")


if __name__ == "__main__":
    main()
