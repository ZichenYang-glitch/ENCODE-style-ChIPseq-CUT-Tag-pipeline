"""Assembly selection for the packaged multi-run overlay fixture.

The default remains two hg38 runs and an mm10 run that cannot share the
window. mm10 and GRCm38 put both authorized runs on the mouse scaffold, write
a local peak that declares mm10, and keep a different assembly (hg38 or mm39)
as the run that must be refused.
"""

from __future__ import annotations

import asyncio
from pathlib import Path
import sys

import httpx
import pytest

_API_TESTS = Path(__file__).resolve().parents[1] / "api"
if str(_API_TESTS) not in sys.path:
    sys.path.insert(0, str(_API_TESTS))

from api_test_client import ApiTestClient  # noqa: E402
from washu_overlay_server import assembly_plan, build_app  # noqa: E402


def test_hg38_plan_keeps_the_mm10_foreign_run() -> None:
    plan = assembly_plan()
    assert plan.label == "hg38"
    assert plan.scaffold == "hg38"
    assert plan.organism == "Homo sapiens"
    assert plan.local_declaration is None
    assert plan.foreign_label == "mm10"
    assert plan.foreign_scaffold == "mm10"


def test_mouse_plans_declare_mm10_and_refuse_a_different_scaffold() -> None:
    mm10 = assembly_plan("mm10")
    grcm38 = assembly_plan(" grcm38 ", foreign="MM39")
    assert mm10.label == "mm10"
    assert mm10.scaffold == "mm10"
    assert mm10.local_declaration == "mm10"
    assert mm10.foreign_label == "hg38"
    assert grcm38.label == "GRCm38"
    assert grcm38.scaffold == "mm10"
    assert grcm38.organism == "Mus musculus"
    assert grcm38.local_declaration == "mm10"
    assert grcm38.foreign_label == "mm39"
    assert grcm38.foreign_scaffold == "mm39"


def test_same_scaffold_and_unknown_labels_are_refused() -> None:
    with pytest.raises(ValueError, match="same mm10 scaffold"):
        assembly_plan("GRCm38", "mm10")
    with pytest.raises(ValueError, match="assembly must be one of"):
        assembly_plan("mm9")
    with pytest.raises(ValueError, match="foreign assembly must be one of"):
        assembly_plan("mm10", "rn6")


def _session(client: ApiTestClient, run_id: str) -> dict:
    response = client.get(f"/api/v1/runs/{run_id}/browser-session")
    assert response.status_code == 200, response.text
    body = response.json()
    assert body["ok"] is True
    return body


def _artifact_assembly(app, run_id: str) -> str:
    artifacts = app.state.run_service.list_artifacts(run_id)
    assert len(artifacts) == 1
    return artifacts[0].metadata["reference_assembly"]


def _organism(app, revision_id: str) -> str:
    summary = app.state.reference_profile_service.get_revision_summary(revision_id)
    return summary.organism


def _anonymous_status(app, path: str) -> int:
    async def request() -> int:
        transport = httpx.ASGITransport(app=app)
        async with httpx.AsyncClient(
            transport=transport,
            base_url="http://testserver",
        ) as client:
            response = await client.get(path)
            await response.aread()
            return response.status_code

    return asyncio.run(request())


def test_default_fixture_keeps_hg38_evidence(tmp_path: Path) -> None:
    app, manifest = build_app(tmp_path / "overlay.sqlite", tmp_path / "workspace")
    assert manifest["assembly"] == "hg38"
    assert manifest["scaffold"] == "hg38"
    assert manifest["foreign_assembly"] == "mm10"
    assert manifest["foreign_scaffold"] == "mm10"
    assert manifest["local_declaration"] is None
    assert manifest["local_bed"] == "overlay-peak-local.bed"
    assert manifest["window"] == "chr1:4500-6700"
    assert not (tmp_path / "workspace" / "overlay-peak-local-mm10.bed").exists()

    client = ApiTestClient(app)
    primary = _session(client, manifest["primary_run_id"])
    second = _session(client, manifest["second_run_id"])
    foreign = _session(client, manifest["foreign_run_id"])
    assert primary["reference"]["assembly_label"] == "hg38"
    assert second["reference"]["assembly_label"] == "hg38"
    assert primary["reference"]["viewer_scaffold"] == "hg38"
    assert second["reference"]["viewer_scaffold"] == "hg38"
    assert (
        primary["reference"]["identity_sha256"]
        != second["reference"]["identity_sha256"]
    )
    assert foreign["reference"]["assembly_label"] == "mm10"
    assert foreign["reference"]["viewer_scaffold"] == "mm10"
    assert _organism(app, primary["reference"]["revision_id"]) == "Homo sapiens"
    assert _organism(app, second["reference"]["revision_id"]) == "Homo sapiens"
    assert _organism(app, foreign["reference"]["revision_id"]) == "Mus musculus"
    assert _artifact_assembly(app, manifest["primary_run_id"]) == "hg38"
    assert _artifact_assembly(app, manifest["second_run_id"]) == "hg38"
    assert _artifact_assembly(app, manifest["foreign_run_id"]) == "mm10"
    for body, generation in (
        (primary, manifest["primary_generation"]),
        (second, manifest["second_generation"]),
    ):
        track = next(
            item
            for item in body["tracks"]
            if item["artifact_id"] == manifest["artifact_id"]
        )
        assert track["preparation"] == "ready"
        assert body["artifact_generation"] == generation
    anonymous = _anonymous_status(
        app,
        "/api/v1/runs/{run}/browser-tracks/{artifact}/bytes/{generation}/{revision}".format(
            run=manifest["primary_run_id"],
            artifact=manifest["artifact_id"],
            generation=manifest["primary_generation"],
            revision=manifest["primary_revision"],
        ),
    )
    assert anonymous == 401


@pytest.mark.parametrize(
    ("assembly", "foreign", "foreign_label", "foreign_organism"),
    [
        ("mm10", None, "hg38", "Homo sapiens"),
        ("GRCm38", "mm39", "mm39", "Mus musculus"),
    ],
)
def test_mouse_runs_share_mm10_and_the_local_peak_declares_it(
    tmp_path: Path,
    assembly: str,
    foreign: str | None,
    foreign_label: str,
    foreign_organism: str,
) -> None:
    workspace = tmp_path / "workspace"
    app, manifest = build_app(
        tmp_path / "overlay.sqlite",
        workspace,
        assembly=assembly,
        foreign=foreign,
    )
    assert manifest["assembly"] == assembly
    assert manifest["scaffold"] == "mm10"
    assert manifest["organism"] == "Mus musculus"
    assert manifest["local_declaration"] == "mm10"
    assert manifest["foreign_assembly"] == foreign_label
    assert manifest["foreign_organism"] == foreign_organism
    assert manifest["foreign_scaffold"] == foreign_label
    bed = Path(manifest["local_bed"])
    assert bed.is_file()
    assert not bed.is_relative_to(workspace)
    text = bed.read_text(encoding="utf-8")
    assert text.startswith("browser genome mm10\n")
    assert "chr1\t5800\t6000\tpeak-local\n" in text

    client = ApiTestClient(app)
    primary = _session(client, manifest["primary_run_id"])
    second = _session(client, manifest["second_run_id"])
    foreign_body = _session(client, manifest["foreign_run_id"])
    assert primary["reference"]["assembly_label"] == assembly
    assert second["reference"]["assembly_label"] == assembly
    assert primary["reference"]["viewer_scaffold"] == "mm10"
    assert second["reference"]["viewer_scaffold"] == "mm10"
    assert (
        primary["reference"]["identity_sha256"]
        != second["reference"]["identity_sha256"]
    )
    assert _organism(app, primary["reference"]["revision_id"]) == "Mus musculus"
    assert _organism(app, second["reference"]["revision_id"]) == "Mus musculus"
    assert _artifact_assembly(app, manifest["primary_run_id"]) == assembly
    assert _artifact_assembly(app, manifest["second_run_id"]) == assembly
    assert foreign_body["reference"]["assembly_label"] == foreign_label
    assert foreign_body["reference"]["viewer_scaffold"] == foreign_label
    assert _organism(app, foreign_body["reference"]["revision_id"]) == foreign_organism
    assert _artifact_assembly(app, manifest["foreign_run_id"]) == foreign_label
    for body in (primary, second):
        track = next(
            item
            for item in body["tracks"]
            if item["artifact_id"] == manifest["artifact_id"]
        )
        assert track["preparation"] == "ready"

    replaced = client.post(
        "/__washu_overlay__/invalidate",
        headers={"x-washu-overlay-token": manifest["token"]},
        json={"run_id": manifest["second_run_id"]},
    )
    assert replaced.status_code == 200, replaced.text
    payload = replaced.json()
    assert payload["ok"] is True
    assert payload["generation"] != manifest["second_generation"]
    assert _artifact_assembly(app, manifest["second_run_id"]) == assembly
    refreshed = _session(client, manifest["second_run_id"])
    assert refreshed["artifact_generation"] == payload["generation"]
    assert refreshed["reference"]["assembly_label"] == assembly
