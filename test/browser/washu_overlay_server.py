"""Serve three synthetic runs and the packaged page on one origin.

By default the opened run and a second run are both hg38, with different
reference profile identities, and a third run is mm10. ``--assembly mm10``
or ``GRCm38`` makes those two runs the mouse target instead. The local peak
then declares mm10, and the third run is a different assembly (hg38 by
default, or mm39). Track bytes stay on the existing session and artifact
routes. POST /__washu_overlay__/invalidate replaces the second run's peak
so its previous generation and revision no longer read.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from dataclasses import dataclass
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
LOCAL_PEAK = "chr1\t5800\t6000\tpeak-local\n"
PRIMARY_NAME = "primary_peaks.narrowPeak"
SECOND_NAME = "second_peaks.narrowPeak"
FOREIGN_NAME = "foreign_peaks.narrowPeak"
DEFAULT_LOCAL_BED = "overlay-peak-local.bed"


@dataclass(frozen=True)
class _Assembly:
    label: str
    scaffold: str
    organism: str


@dataclass(frozen=True)
class AssemblyPlan:
    """Target the proof asked for, plus the one assembly that must be refused."""

    label: str
    scaffold: str
    organism: str
    local_declaration: str | None
    foreign_label: str
    foreign_scaffold: str
    foreign_organism: str


_TARGETS = {
    "hg38": _Assembly("hg38", "hg38", "Homo sapiens"),
    "mm10": _Assembly("mm10", "mm10", "Mus musculus"),
    "grcm38": _Assembly("GRCm38", "mm10", "Mus musculus"),
}
_FOREIGN = {
    "hg38": _Assembly("hg38", "hg38", "Homo sapiens"),
    "mm10": _Assembly("mm10", "mm10", "Mus musculus"),
    "mm39": _Assembly("mm39", "mm39", "Mus musculus"),
}


def assembly_plan(assembly: str = "hg38", foreign: str | None = None) -> AssemblyPlan:
    """Resolve a target label and a scaffold that cannot share its window.

    hg38 keeps the undeclared local fixture. mm10 and GRCm38 both navigate
    on the mm10 scaffold, and the synthetic local peak declares mm10 either way.
    """

    target = _TARGETS.get(assembly.strip().lower())
    if target is None:
        known = ", ".join(item.label for item in _TARGETS.values())
        raise ValueError(f"assembly must be one of {known}")
    if foreign is None:
        foreign_name = "mm10" if target.scaffold == "hg38" else "hg38"
    else:
        foreign_name = foreign
    other = _FOREIGN.get(foreign_name.strip().lower())
    if other is None:
        known = ", ".join(item.label for item in _FOREIGN.values())
        raise ValueError(f"foreign assembly must be one of {known}")
    if other.scaffold == target.scaffold:
        raise ValueError(
            f"{other.label} uses the same {target.scaffold} scaffold as {target.label}"
        )
    declaration = None if target.scaffold == "hg38" else "mm10"
    return AssemblyPlan(
        label=target.label,
        scaffold=target.scaffold,
        organism=target.organism,
        local_declaration=declaration,
        foreign_label=other.label,
        foreign_scaffold=other.scaffold,
        foreign_organism=other.organism,
    )


def _write_declared_bed(workspace: Path, declaration: str) -> str:
    workspace.mkdir(parents=True, exist_ok=True)
    path = workspace / "overlay-peak-local-mm10.bed"
    path.write_text(
        f"browser genome {declaration}\n{LOCAL_PEAK}",
        encoding="utf-8",
    )
    return str(path)


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


def _install_invalidate(
    app,
    client: _Client,
    second_run_id: str,
    token: str,
    assembly: str,
) -> None:
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
            artifact = _peak(client, second_run_id, SECOND_NAME, content, assembly)
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


def build_app(
    database: Path,
    workspace: Path,
    assembly: str = "hg38",
    foreign: str | None = None,
):
    plan = assembly_plan(assembly, foreign)
    app = create_app(
        database_url=f"sqlite:///{database}",
        workspace_root=workspace,
    )
    seed_test_authentication(app)
    client = _Client(app)
    primary_run_id = _create_run(client)
    if plan.label == "hg38":
        # The default primary profile is the shared browser fixture, so the
        # hg38 evidence keeps that identity, organism, and display name.
        _bind_reference(client, primary_run_id, assembly="hg38")
    else:
        _bind_distinct(
            client,
            primary_run_id,
            assembly=plan.label,
            organism=plan.organism,
            marker="56" * 16,
            display_name=f"Overlay primary {plan.label}",
        )
    primary = _peak(client, primary_run_id, PRIMARY_NAME, PRIMARY_PEAK, plan.label)
    _prepare(app, client, primary_run_id, primary)

    second_run_id = _create_run(client)
    second_name = "Overlay second hg38" if plan.label == "hg38" else (
        f"Overlay second {plan.label}"
    )
    _bind_distinct(
        client,
        second_run_id,
        assembly=plan.label,
        organism=plan.organism,
        marker="12" * 16,
        display_name=second_name,
    )
    second = _peak(client, second_run_id, SECOND_NAME, SECOND_PEAK, plan.label)
    _prepare(app, client, second_run_id, second)

    foreign_run_id = _create_run(client)
    foreign_name = (
        "Overlay mouse mm10"
        if plan.foreign_label == "mm10"
        else f"Overlay foreign {plan.foreign_label}"
    )
    _bind_distinct(
        client,
        foreign_run_id,
        assembly=plan.foreign_label,
        organism=plan.foreign_organism,
        marker="34" * 16,
        display_name=foreign_name,
    )
    foreign = _peak(
        client,
        foreign_run_id,
        FOREIGN_NAME,
        FOREIGN_PEAK,
        plan.foreign_label,
    )
    _prepare(app, client, foreign_run_id, foreign)

    if plan.local_declaration is None:
        local_bed = DEFAULT_LOCAL_BED
    else:
        local_bed = _write_declared_bed(workspace, plan.local_declaration)
    token = secrets.token_hex(24)
    _install_invalidate(app, client, second_run_id, token, plan.label)
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
        "assembly": plan.label,
        "scaffold": plan.scaffold,
        "organism": plan.organism,
        "foreign_assembly": plan.foreign_label,
        "foreign_scaffold": plan.foreign_scaffold,
        "foreign_organism": plan.foreign_organism,
        "local_declaration": plan.local_declaration,
        "local_bed": local_bed,
        "local_peak": "peak-local",
        "token": token,
    }
    return app, manifest


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--assembly", default="hg38")
    parser.add_argument("--foreign", default=None)
    args = parser.parse_args()
    try:
        _app, manifest = build_app(
            args.database,
            args.workspace,
            assembly=args.assembly,
            foreign=args.foreign,
        )
    except ValueError as exc:
        parser.error(str(exc))
    print(
        "WASHU_OVERLAY " + json.dumps(manifest, separators=(",", ":")),
        flush=True,
    )
    uvicorn.run(_app, host="127.0.0.1", port=args.port, log_level="warning")


if __name__ == "__main__":
    main()
