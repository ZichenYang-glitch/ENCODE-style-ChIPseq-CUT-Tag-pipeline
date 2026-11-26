"""Black-box acceptance tests for strict and per_group_fallback 2/3 routing.

Covers public CLI (`validate_samples` / `encode-validate`, `make_manifest` /
`encode-manifest`), `EncodeStyleWorkflowAdapter` (file-path, inline, workspace
planning, artifact extraction, and WashU browser descriptions), real Snakemake
dry-runs, run-identity isolation and resumption, and lightweight real consensus
execution without real sequencing reads.
"""

from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pytest
import yaml

from conftest import seed_run_identity
from encode_pipeline.adapters.browser.catalog import describe_browser_artifact
from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.config.scientific_policy import SCIENTIFIC_POLICY_REVISION
from encode_pipeline.errors import ValidationError
from encode_pipeline.manifest.make import build_manifest_rows
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.services.materialization import WorkspaceMaterializer


SAMPLE_COLUMNS = (
    "sample",
    "layout",
    "fastq_1",
    "fastq_2",
    "assay",
    "peak_mode",
    "genome",
    "bowtie2_index",
    "chrom_sizes",
    "target",
    "experiment",
    "role",
    "biological_replicate",
    "technical_replicate",
)

WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
EXPECTED_MANIFEST_COLUMNS = [
    "sample_id",
    "experiment_id",
    "assay",
    "target",
    "genome",
    "output_type",
    "method",
    "path",
    "status",
    "qc_flag",
]


def _write_synthetic_inputs(root: Path) -> dict[str, Path]:
    inp = root / "inputs"
    inp.mkdir(parents=True, exist_ok=True)
    r1 = inp / "r1.fq.gz"
    r2 = inp / "r2.fq.gz"
    r1.write_bytes(b"")
    r2.write_bytes(b"")
    bt2 = inp / "hg38"
    for ext in (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"):
        (inp / f"hg38{ext}").write_bytes(b"")
    cs = inp / "hg38.chrom.sizes"
    cs.write_text("chr1\t248956422\n", encoding="utf-8")
    return {"r1": r1, "r2": r2, "bt2": bt2, "cs": cs}


def _render_samples_tsv(
    inputs: dict[str, Path],
    groups: list[tuple[str, str, str, tuple[int, ...], str]],
) -> str:
    lines = ["\t".join(SAMPLE_COLUMNS)]
    for exp, assay, peak_mode, bioreps, layout in groups:
        target = "H3K27me3" if peak_mode == "broad" else "CTCF"
        technical = Counter()
        for br in bioreps:
            technical[br] += 1
            sample = f"{exp}_r{br}"
            if technical[br] > 1:
                sample += f"_t{technical[br]}"
            fq2 = str(inputs["r2"]) if layout == "PE" else ""
            row = [
                sample,
                layout,
                str(inputs["r1"]),
                fq2,
                assay,
                peak_mode,
                "hs",
                str(inputs["bt2"]),
                str(inputs["cs"]),
                target,
                exp,
                "treatment",
                str(br),
                str(technical[br]),
            ]
            lines.append("\t".join(row))
    return "\n".join(lines) + "\n"


def _build_fixture_workspace(
    tmp_path: Path,
    tmp_config,
    *,
    name: str,
    groups: list[tuple[str, str, str, tuple[int, ...], str]],
    extra_config: dict,
) -> tuple[Path, Path, Path, dict, list[dict[str, str]]]:
    inputs = _write_synthetic_inputs(tmp_path / f"inp_{name}")
    samples_text = _render_samples_tsv(inputs, groups)
    inline_samples = list(csv.DictReader(samples_text.splitlines(), delimiter="\t"))
    for row in inline_samples:
        # The adapter takes chrom_sizes from the config for inline rows.
        del row["chrom_sizes"]
    cfg = {
        "outdir": str(tmp_path / "results"),
        "genome": "hs",
        "bowtie2_index": str(inputs["bt2"]),
        "chrom_sizes": str(inputs["cs"]),
        "use_control": False,
        "replicate_analysis": True,
        "threads": 1,
        "chipseq_idr": False,
        "blacklist_mode": "off",
    }
    cfg.update(extra_config)
    workdir, config_path, samples_path = tmp_config(config=cfg, samples=samples_text)
    # Keep --configfile as the only source for Snakemake, including on resume.
    # Stage the adapter's config/config.yaml only when its consumer needs it.
    loaded_cfg = yaml.safe_load(Path(config_path).read_text(encoding="utf-8"))
    return workdir, Path(config_path), Path(samples_path), loaded_cfg, inline_samples


def _parse_summary(stdout: str) -> dict[str, dict[str, str]]:
    parsed = {}
    for line in stdout.splitlines():
        cols = line.split("\t")
        if len(cols) >= 6 and cols[0] != "output_file":
            parsed[cols[0]] = {
                "date": cols[1],
                "rule": cols[2],
                "log": cols[3],
                "status": cols[4],
                "plan": cols[5],
            }
    return parsed


def _assert_three_replicate_command(run_snakemake, config_path, target, *, minimum):
    # A summary lists the default DAG's files. A separate dry-run prints the
    # command for exactly this experiment and caller.
    sm = run_snakemake(config_path, extra_args=["--printshellcmds", str(target)])
    output = sm.stdout + sm.stderr
    assert sm.returncode == 0, output
    assert "--bioreps 1 2 3" in output
    assert f"--min-replicates {minimum}" in output


def _extract_workspace_artifacts(adapter, inputs, workdir):
    plan = adapter.plan_workspace(inputs, workdir)
    assert plan.is_success, plan.issues
    materialized = WorkspaceMaterializer().materialize(plan.value, workdir)
    assert materialized.is_success, materialized.issues
    return adapter.extract_artifacts(inputs, workdir)


def _assert_final_track(item):
    desc = describe_browser_artifact(
        workflow_id=WORKFLOW_ID,
        output_type=item.output_type,
        name=Path(item.relative_path).name,
    )
    assert desc is not None and desc.available
    assert (desc.wash_u_type, desc.role) == ("bed", "interval")
    reading = f"{desc.detail} {desc.score_meaning}".lower()
    if "idr" in item.output_type:
        assert "idr" in reading or "rank" in reading
    else:
        assert "support" in reading and "replicate" in reading


def _materialize_present_outputs(
    manifest_rows: list[dict],
    *,
    output_types: set[str],
) -> None:
    for row in manifest_rows:
        if row["output_type"] not in output_types:
            continue
        path = Path(row["path"])
        if not row["path"] or path.suffix == "":
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        peak = (
            "chr1\t10\t20\n"
            if path.suffix == ".bed"
            else "chr1\t10\t20\tp1\t1000\t.\t1.0\t-1\t-1\t5\n"
        )
        path.write_text(peak, encoding="utf-8")


def _run_validate_cli(config_path: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "encode_pipeline.cli.validate",
            "--config",
            str(config_path),
        ],
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )


def _run_manifest_cli(
    config_path: Path, output_path: Path
) -> tuple[subprocess.CompletedProcess[str], list[dict[str, str]]]:
    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "encode_pipeline.cli.manifest",
            "--config",
            str(config_path),
            "--output",
            str(output_path),
        ],
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    rows: list[dict[str, str]] = []
    if proc.returncode == 0 and output_path.exists():
        with output_path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            assert reader.fieldnames == EXPECTED_MANIFEST_COLUMNS
            rows = list(reader)
    return proc, rows


def _run_snakemake_identity_targets(
    workdir: Path,
    snakefile: str,
    config_path: Path,
    marker: Path,
    eff: Path,
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "-s",
            snakefile,
            "--configfile",
            str(config_path),
            "--profile",
            "none",
            "--workflow-profile",
            "none",
            "--cores",
            "1",
            str(marker),
            str(eff),
        ],
        cwd=workdir,
        capture_output=True,
        text=True,
        env={
            **os.environ,
            "PYTHONDONTWRITEBYTECODE": "1",
            "XDG_CACHE_HOME": str(workdir / ".cache"),
        },
    )


def _assert_all_public_entrypoints_reject(
    *,
    workdir: Path,
    config_path: Path,
    samples_path: Path,
    cfg: dict,
    inline_samples: list[dict[str, str]],
    run_validator,
    run_snakemake,
) -> None:
    """Assert every public entrypoint rejects the configuration pre-execution."""
    val = run_validator(config_path)
    assert val.returncode != 0
    assert val.stderr.strip()

    cli_val = _run_validate_cli(config_path)
    assert cli_val.returncode != 0
    assert cli_val.stderr.strip()

    sm = run_snakemake(config_path, extra_args=["--summary"])
    assert sm.returncode != 0

    adapter = EncodeStyleWorkflowAdapter()
    path_res = adapter.validate(WorkflowInputs(config=cfg))
    assert not path_res.is_success
    assert path_res.issues

    inline_cfg = dict(cfg)
    inline_cfg.pop("samples", None)
    inline_res = adapter.validate(
        WorkflowInputs(config=inline_cfg, samples=inline_samples)
    )
    assert not inline_res.is_success
    assert inline_res.issues

    with pytest.raises((ValidationError, ValueError)):
        build_manifest_rows(cfg, samples_path=str(samples_path))

    manifest_proc, _ = _run_manifest_cli(config_path, workdir / "rejected_manifest.tsv")
    assert manifest_proc.returncode != 0


# ---------------------------------------------------------------------------
# 1. Strict-policy baseline & cross-mode identity isolation
# ---------------------------------------------------------------------------


def test_strict_policy_pure_two_chipseq_idr_and_non_rewrite_identity(
    tmp_path, tmp_config, run_validator, run_snakemake, snakefile
):
    assert SCIENTIFIC_POLICY_REVISION == 7
    workdir, config_path, samples_path, cfg, _ = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name="pure2_chip",
        groups=[("GRP2", "chipseq", "narrow", (1, 2), "PE")],
        extra_config={
            "chipseq_idr": True,
            "reproducibility": {"enabled": True, "consensus": {"enabled": True}},
        },
    )

    assert run_validator(config_path).returncode == 0
    assert _run_validate_cli(config_path).returncode == 0
    adapter = EncodeStyleWorkflowAdapter()
    assert adapter.validate(WorkflowInputs(config=cfg)).is_success

    sm = run_snakemake(config_path, extra_args=["--summary"])
    assert sm.returncode == 0, sm.stderr
    summary = _parse_summary(sm.stdout)
    outdir = Path(cfg["outdir"])
    assert (
        str(outdir / "experiments/GRP2/06_idr/final/conservative.narrowPeak") in summary
    )
    assert str(outdir / "experiments/GRP2/06_idr/final/optimal.narrowPeak") in summary

    rows, _, _ = build_manifest_rows(cfg, samples_path=str(samples_path))
    pairs = {(r["experiment_id"], r["output_type"]) for r in rows}
    assert ("GRP2", "idr_conservative") in pairs
    assert ("GRP2", "idr_optimal") in pairs
    assert ("GRP2", "chipseq_macs3_narrow_consensus_peak") in pairs

    # I-01 / I-01b: marker + effective_config non-rewrite vs threads 1->2 control
    marker = outdir / ".helixweave/analysis_identity.json"
    eff = outdir / ".helixweave/effective_config.json"
    create = _run_snakemake_identity_targets(
        workdir, snakefile, config_path, marker, eff
    )
    assert create.returncode == 0, create.stderr
    marker_bytes = marker.read_bytes()
    eff_bytes = eff.read_bytes()
    assert json.loads(marker_bytes)["scientific_policy_revision"] == 7

    sm_reuse = run_snakemake(
        config_path, extra_args=["--summary", str(marker), str(eff)]
    )
    assert sm_reuse.returncode == 0, sm_reuse.stderr
    reuse_summary = _parse_summary(sm_reuse.stdout)
    for rec in (marker, eff):
        assert reuse_summary[str(rec)]["status"] == "ok"
        assert reuse_summary[str(rec)]["plan"] == "no update"

    reexec = _run_snakemake_identity_targets(
        workdir, snakefile, config_path, marker, eff
    )
    assert reexec.returncode == 0, reexec.stderr
    assert marker.read_bytes() == marker_bytes
    assert eff.read_bytes() == eff_bytes

    # Explicit strict resumes the omitted-key directory with identical identity.
    cfg_explicit_strict = dict(cfg)
    cfg_explicit_strict["reproducibility"] = {
        "enabled": True,
        "mixed_replicate_policy": "strict",
        "consensus": {"enabled": True},
    }
    cfg_strict_path = workdir / "config_explicit_strict.yaml"
    cfg_strict_path.write_text(
        yaml.safe_dump(cfg_explicit_strict, sort_keys=False), encoding="utf-8"
    )
    sm_strict_reuse = run_snakemake(
        cfg_strict_path, extra_args=["--summary", str(marker), str(eff)]
    )
    assert sm_strict_reuse.returncode == 0, sm_strict_reuse.stderr
    strict_summary = _parse_summary(sm_strict_reuse.stdout)
    for rec in (marker, eff):
        assert strict_summary[str(rec)]["status"] == "ok"
        assert strict_summary[str(rec)]["plan"] == "no update"
    reexec_strict = _run_snakemake_identity_targets(
        workdir, snakefile, cfg_strict_path, marker, eff
    )
    assert reexec_strict.returncode == 0, reexec_strict.stderr
    assert marker.read_bytes() == marker_bytes
    assert eff.read_bytes() == eff_bytes

    # Cross-mode isolation: per_group_fallback must refuse to write into strict outdir.
    cfg_fallback = dict(cfg)
    cfg_fallback["reproducibility"] = {
        "enabled": True,
        "mixed_replicate_policy": "per_group_fallback",
        "consensus": {"enabled": True},
    }
    cfg_fallback_path = workdir / "config_fallback.yaml"
    cfg_fallback_path.write_text(
        yaml.safe_dump(cfg_fallback, sort_keys=False), encoding="utf-8"
    )
    cross_mode = _run_snakemake_identity_targets(
        workdir, snakefile, cfg_fallback_path, marker, eff
    )
    assert cross_mode.returncode != 0
    assert marker.read_bytes() == marker_bytes
    assert eff.read_bytes() == eff_bytes

    # Negative control: changing threads=2 triggers params changed / update pending
    cfg_t2 = dict(cfg)
    cfg_t2["threads"] = 2
    cfg_t2_path = workdir / "config_t2.yaml"
    cfg_t2_path.write_text(yaml.safe_dump(cfg_t2, sort_keys=False), encoding="utf-8")
    sm_t2 = run_snakemake(cfg_t2_path, extra_args=["--summary", str(marker), str(eff)])
    assert sm_t2.returncode == 0, sm_t2.stderr
    t2_summary = _parse_summary(sm_t2.stdout)
    assert t2_summary[str(marker)]["status"] == "params changed"
    assert t2_summary[str(marker)]["plan"] == "update pending"
    assert t2_summary[str(eff)]["status"] == "params changed"
    assert t2_summary[str(eff)]["plan"] == "update pending"


def test_strict_policy_mixed_2_3_with_idr_off_routes_both_to_consensus(
    tmp_path, tmp_config, run_validator, run_snakemake
):
    _, config_path, samples_path, cfg, _ = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name="mix23_noidr",
        groups=[
            ("GRP2", "cuttag", "narrow", (1, 2), "PE"),
            ("GRP3", "cuttag", "narrow", (1, 2, 3), "PE"),
        ],
        extra_config={
            "cuttag": {"seacr": {"enabled": True, "mode": "stringent"}},
            "reproducibility": {"enabled": True, "consensus": {"enabled": True}},
        },
    )

    assert run_validator(config_path).returncode == 0
    adapter = EncodeStyleWorkflowAdapter()
    assert adapter.validate(WorkflowInputs(config=cfg)).is_success

    sm = run_snakemake(config_path, extra_args=["--summary"])
    assert sm.returncode == 0, sm.stderr
    summary = _parse_summary(sm.stdout)
    outdir = Path(cfg["outdir"])
    for exp in ("GRP2", "GRP3"):
        assert (
            str(
                outdir / f"experiments/{exp}/06_reproducibility/final/"
                f"{exp}.cuttag.macs3.narrow.replicate_validated.consensus.narrowPeak"
            )
            in summary
        )
        assert (
            str(
                outdir / f"experiments/{exp}/06_reproducibility/final/"
                f"{exp}.cuttag.seacr.stringent.replicate_validated.consensus.bed"
            )
            in summary
        )

    rows, _, _ = build_manifest_rows(cfg, samples_path=str(samples_path))
    pairs = {(r["experiment_id"], r["output_type"]) for r in rows}
    for exp in ("GRP2", "GRP3"):
        assert (exp, "cuttag_macs3_narrow_consensus_final_peak") in pairs
        assert (exp, "cuttag_seacr_consensus_final_peak") in pairs


def test_strict_policy_cross_mode_mixed_2_3_passes_end_to_end(
    tmp_path, tmp_config, run_validator, run_snakemake
):
    _, config_path, samples_path, cfg, _ = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name="cross_mode",
        groups=[
            ("GRP2_N", "chipseq", "narrow", (1, 2), "PE"),
            ("GRP3_B", "chipseq", "broad", (1, 2, 3), "PE"),
        ],
        extra_config={
            "chipseq_idr": True,
            "reproducibility": {"enabled": True, "consensus": {"enabled": True}},
        },
    )

    assert run_validator(config_path).returncode == 0
    adapter = EncodeStyleWorkflowAdapter()
    assert adapter.validate(WorkflowInputs(config=cfg)).is_success

    sm = run_snakemake(config_path, extra_args=["--summary"])
    assert sm.returncode == 0, sm.stderr
    summary = _parse_summary(sm.stdout)
    outdir = Path(cfg["outdir"])
    assert (
        str(outdir / "experiments/GRP2_N/06_idr/final/conservative.narrowPeak")
        in summary
    )
    assert (
        str(
            outdir / "experiments/GRP3_B/06_reproducibility/final/"
            "GRP3_B.chipseq.macs3.broad.replicate_validated.consensus.broadPeak"
        )
        in summary
    )

    rows, _, _ = build_manifest_rows(cfg, samples_path=str(samples_path))
    pairs = {(r["experiment_id"], r["output_type"]) for r in rows}
    assert ("GRP2_N", "idr_conservative") in pairs
    assert ("GRP3_B", "chipseq_macs3_broad_consensus_final_peak") in pairs


@pytest.mark.parametrize("policy_override", [None, "strict"], ids=["omitted", "strict"])
@pytest.mark.parametrize(
    ("assay", "peak_mode", "idr_key"),
    [
        ("chipseq", "narrow", "chipseq_narrow"),
        ("atac", "narrow", "atac_narrow"),
        ("cuttag", "narrow", "cuttag_narrow"),
        ("chipseq", "broad", "chipseq_broad_experimental"),
        ("cuttag", "broad", "cuttag_broad_experimental"),
    ],
)
def test_strict_policy_same_mode_mixed_2_3_rejected_consistently(
    tmp_path,
    tmp_config,
    run_validator,
    run_snakemake,
    policy_override,
    assay,
    peak_mode,
    idr_key,
):
    repro = {
        "enabled": True,
        "consensus": {"enabled": True},
        "idr": {idr_key: True},
    }
    if policy_override is not None:
        repro["mixed_replicate_policy"] = policy_override

    workdir, config_path, samples_path, cfg, inline_samples = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name=f"reject_{policy_override or 'omitted'}_{idr_key}",
        groups=[
            ("GRP2", assay, peak_mode, (1, 2), "PE"),
            ("GRP3", assay, peak_mode, (1, 2, 3), "PE"),
        ],
        extra_config={"reproducibility": repro},
    )

    _assert_all_public_entrypoints_reject(
        workdir=workdir,
        config_path=config_path,
        samples_path=samples_path,
        cfg=cfg,
        inline_samples=inline_samples,
        run_validator=run_validator,
        run_snakemake=run_snakemake,
    )


def test_cuttag_independent_seacr_coexists_with_macs3_idr_d00_m00(
    tmp_path, tmp_config, run_validator, run_snakemake
):
    workdir, config_path, samples_path, cfg, _ = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name="d00_m00",
        groups=[("CT2", "cuttag", "narrow", (1, 2), "PE")],
        extra_config={
            "cuttag": {
                "seacr": {
                    "enabled": True,
                    "mode": "stringent",
                    "threshold": 0.01,
                    "normalization": "non",
                }
            },
            "reproducibility": {
                "enabled": True,
                "consensus": {"enabled": True, "min_replicates": 2},
                "idr": {"cuttag_narrow": True},
            },
        },
    )

    assert run_validator(config_path).returncode == 0
    sm = run_snakemake(config_path, extra_args=["--summary"])
    assert sm.returncode == 0, sm.stderr
    summary = _parse_summary(sm.stdout)
    outdir = Path(cfg["outdir"])
    idr_final = (
        outdir / "experiments/CT2/06_reproducibility/final/"
        "CT2.cuttag.macs3.narrow.replicate_validated.idr.narrowPeak"
    )
    seacr_final = (
        outdir / "experiments/CT2/06_reproducibility/final/"
        "CT2.cuttag.seacr.stringent.replicate_validated.consensus.bed"
    )
    macs3_cons_final = (
        outdir / "experiments/CT2/06_reproducibility/final/"
        "CT2.cuttag.macs3.narrow.replicate_validated.consensus.narrowPeak"
    )
    assert str(idr_final) in summary
    assert str(seacr_final) in summary
    assert str(macs3_cons_final) not in summary

    rows, _, _ = build_manifest_rows(cfg, samples_path=str(samples_path))
    pairs = {(r["experiment_id"], r["output_type"]) for r in rows}
    assert ("CT2", "cuttag_macs3_narrow_idr_final_peak") in pairs
    assert ("CT2", "cuttag_seacr_consensus_final_peak") in pairs
    assert ("CT2", "cuttag_macs3_narrow_consensus_final_peak") not in pairs

    _materialize_present_outputs(
        rows,
        output_types={
            "cuttag_macs3_narrow_idr_final_peak",
            "cuttag_seacr_consensus_final_peak",
        },
    )
    adapter = EncodeStyleWorkflowAdapter()
    extracted = _extract_workspace_artifacts(
        adapter, WorkflowInputs(config=cfg), workdir
    )
    assert extracted.is_success, extracted.issues
    extracted_types = {item.output_type for item in extracted.value}
    assert {
        "cuttag_macs3_narrow_idr_final_peak",
        "cuttag_seacr_consensus_final_peak",
    } <= extracted_types
    for item in extracted.value:
        desc = describe_browser_artifact(
            workflow_id=WORKFLOW_ID,
            output_type=item.output_type,
            name=Path(item.relative_path).name,
        )
        assert desc is not None
        assert desc.wash_u_type == "bed"
        assert desc.role == "interval"


# ---------------------------------------------------------------------------
# 2. per_group_fallback acceptance: narrow 2/3 routing across public surfaces
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "switch_mode",
    ["chipseq_idr_top", "reproducibility_idr_chipseq_narrow"],
)
def test_per_group_fallback_chipseq_narrow_preserves_idr_and_promotes_3rep_consensus(
    tmp_path,
    tmp_config,
    run_validator,
    run_snakemake,
    switch_mode,
):
    repro: dict = {
        "enabled": True,
        "mixed_replicate_policy": "per_group_fallback",
        "consensus": {"enabled": True},
    }
    extra_cfg: dict = {"reproducibility": repro}
    if switch_mode == "chipseq_idr_top":
        extra_cfg["chipseq_idr"] = True
    else:
        extra_cfg["chipseq_idr"] = False
        repro["idr"] = {"chipseq_narrow": True}

    workdir, config_path, samples_path, cfg, inline_samples = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name=f"fb_chip_{switch_mode}",
        groups=[
            # Both groups have an extra technical replicate. Only distinct
            # biological replicates determine the published result.
            ("GRP2", "chipseq", "narrow", (1, 1, 2), "PE"),
            ("GRP3", "chipseq", "narrow", (1, 2, 2, 3), "PE"),
        ],
        extra_config=extra_cfg,
    )
    outdir = Path(cfg["outdir"])

    assert run_validator(config_path).returncode == 0
    cli_val = _run_validate_cli(config_path)
    assert cli_val.returncode == 0, cli_val.stderr

    adapter = EncodeStyleWorkflowAdapter()
    assert adapter.validate(WorkflowInputs(config=cfg)).is_success
    inline_cfg = dict(cfg)
    inline_cfg.pop("samples", None)
    assert adapter.validate(
        WorkflowInputs(config=inline_cfg, samples=inline_samples)
    ).is_success
    sm = run_snakemake(config_path, extra_args=["--summary"])
    assert sm.returncode == 0, sm.stdout + sm.stderr
    dag_targets = _parse_summary(sm.stdout)

    # 2-rep ChIP narrow keeps conservative/optimal IDR and no consensus final.
    grp2_cons_peak = (
        outdir / "experiments/GRP2/06_reproducibility/consensus/"
        "GRP2.chipseq.macs3.narrow.consensus.narrowPeak"
    )
    grp2_cons_final = (
        outdir / "experiments/GRP2/06_reproducibility/final/"
        "GRP2.chipseq.macs3.narrow.replicate_validated.consensus.narrowPeak"
    )
    assert (
        str(outdir / "experiments/GRP2/06_idr/final/conservative.narrowPeak")
        in dag_targets
    )
    assert (
        str(outdir / "experiments/GRP2/06_idr/final/optimal.narrowPeak") in dag_targets
    )
    assert (
        str(outdir / "experiments/GRP2/06_idr/final/reproducibility_summary.tsv")
        in dag_targets
    )
    assert str(grp2_cons_peak) in dag_targets
    assert str(grp2_cons_final) not in dag_targets

    # 3-rep ChIP narrow runs 2-of-3 consensus across all 3 bioreps and promotes final.
    grp3_cons_peak = (
        outdir / "experiments/GRP3/06_reproducibility/consensus/"
        "GRP3.chipseq.macs3.narrow.consensus.narrowPeak"
    )
    grp3_cons_summary = (
        outdir / "experiments/GRP3/06_reproducibility/consensus/"
        "GRP3.chipseq.macs3.narrow.consensus.summary.tsv"
    )
    grp3_cons_final = (
        outdir / "experiments/GRP3/06_reproducibility/final/"
        "GRP3.chipseq.macs3.narrow.replicate_validated.consensus.narrowPeak"
    )
    for br in (1, 2, 3):
        assert (
            str(
                outdir / "experiments/GRP3/06_reproducibility/consensus/biorep_peaks/"
                f"GRP3.biorep{br}.chipseq.macs3.narrow_peaks.narrowPeak"
            )
            in dag_targets
        )
    assert str(grp3_cons_peak) in dag_targets
    assert str(grp3_cons_summary) in dag_targets
    assert str(grp3_cons_final) in dag_targets
    assert (
        str(outdir / "experiments/GRP3/06_idr/final/conservative.narrowPeak")
        not in dag_targets
    )
    assert (
        str(outdir / "experiments/GRP3/06_idr/final/optimal.narrowPeak")
        not in dag_targets
    )

    _assert_three_replicate_command(
        run_snakemake, config_path, grp3_cons_peak, minimum=2
    )

    rows, _, _ = build_manifest_rows(cfg, samples_path=str(samples_path))
    manifest_proc, cli_rows = _run_manifest_cli(
        config_path, workdir / "manifest_chip.tsv"
    )
    assert manifest_proc.returncode == 0, manifest_proc.stderr
    assert {(r["experiment_id"], r["output_type"], r["path"]) for r in rows} == {
        (r["experiment_id"], r["output_type"], r["path"]) for r in cli_rows
    }

    pairs = {(r["experiment_id"], r["output_type"]) for r in rows}
    assert ("GRP2", "idr_conservative") in pairs
    assert ("GRP2", "idr_optimal") in pairs
    assert ("GRP2", "idr_reproducibility_summary") in pairs
    assert ("GRP2", "chipseq_macs3_narrow_consensus_peak") in pairs
    assert ("GRP2", "chipseq_macs3_narrow_consensus_summary") in pairs
    assert ("GRP2", "chipseq_macs3_narrow_consensus_final_peak") not in pairs

    assert ("GRP3", "chipseq_macs3_narrow_consensus_peak") in pairs
    assert ("GRP3", "chipseq_macs3_narrow_consensus_summary") in pairs
    assert ("GRP3", "chipseq_macs3_narrow_consensus_final_peak") in pairs
    assert ("GRP3", "idr_conservative") not in pairs
    assert ("GRP3", "idr_optimal") not in pairs

    final_types = {
        "idr_conservative",
        "idr_optimal",
        "chipseq_macs3_narrow_consensus_final_peak",
    }
    _materialize_present_outputs(rows, output_types=final_types)
    extracted = _extract_workspace_artifacts(
        adapter, WorkflowInputs(config=cfg), workdir
    )
    assert extracted.is_success, extracted.issues
    extracted_pairs = {
        (item.metadata.get("experiment_id"), item.output_type)
        for item in extracted.value
    }
    assert {
        ("GRP2", "idr_conservative"),
        ("GRP2", "idr_optimal"),
        ("GRP3", "chipseq_macs3_narrow_consensus_final_peak"),
    } <= extracted_pairs
    for item in extracted.value:
        _assert_final_track(item)


@pytest.mark.parametrize(
    (
        "assay",
        "idr_key",
        "idr_final_type",
        "idr_summary_type",
        "cons_peak_type",
        "cons_summary_type",
        "cons_final_type",
        "seacr",
    ),
    [
        (
            "atac",
            "atac_narrow",
            "atac_macs3_narrow_idr_final_peak",
            "atac_macs3_narrow_idr_summary",
            "atac_macs3_narrow_consensus_peak",
            "atac_macs3_narrow_consensus_summary",
            "atac_macs3_narrow_consensus_final_peak",
            False,
        ),
        (
            "cuttag",
            "cuttag_narrow",
            "cuttag_macs3_narrow_idr_final_peak",
            "cuttag_macs3_narrow_idr_summary",
            "cuttag_macs3_narrow_consensus_peak",
            "cuttag_macs3_narrow_consensus_summary",
            "cuttag_macs3_narrow_consensus_final_peak",
            False,
        ),
        (
            "cuttag",
            "cuttag_narrow",
            "cuttag_macs3_narrow_idr_final_peak",
            "cuttag_macs3_narrow_idr_summary",
            "cuttag_macs3_narrow_consensus_peak",
            "cuttag_macs3_narrow_consensus_summary",
            "cuttag_macs3_narrow_consensus_final_peak",
            True,
        ),
    ],
    ids=["atac-narrow", "cuttag-narrow", "cuttag-narrow-seacr"],
)
def test_per_group_fallback_atac_and_cuttag_narrow_routing_and_seacr_coexistence(
    tmp_path,
    tmp_config,
    run_validator,
    run_snakemake,
    assay,
    idr_key,
    idr_final_type,
    idr_summary_type,
    cons_peak_type,
    cons_summary_type,
    cons_final_type,
    seacr,
):
    extra_cfg: dict = {
        "reproducibility": {
            "enabled": True,
            "mixed_replicate_policy": "per_group_fallback",
            "consensus": {"enabled": True, "min_replicates": 2},
            "idr": {idr_key: True},
        }
    }
    if seacr:
        extra_cfg["cuttag"] = {"seacr": {"enabled": True, "mode": "stringent"}}

    workdir, config_path, samples_path, cfg, inline_samples = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name=f"fb_{idr_key}_{'seacr' if seacr else 'noseacr'}",
        groups=[
            ("GRP2", assay, "narrow", (1, 2), "PE"),
            ("GRP3", assay, "narrow", (1, 2, 3), "PE"),
        ],
        extra_config=extra_cfg,
    )
    outdir = Path(cfg["outdir"])

    assert run_validator(config_path).returncode == 0
    cli_val = _run_validate_cli(config_path)
    assert cli_val.returncode == 0, cli_val.stderr

    adapter = EncodeStyleWorkflowAdapter()
    assert adapter.validate(WorkflowInputs(config=cfg)).is_success
    inline_cfg = dict(cfg)
    inline_cfg.pop("samples", None)
    assert adapter.validate(
        WorkflowInputs(config=inline_cfg, samples=inline_samples)
    ).is_success
    sm = run_snakemake(config_path, extra_args=["--summary"])
    assert sm.returncode == 0, sm.stdout + sm.stderr
    dag_targets = _parse_summary(sm.stdout)

    grp2_idr = str(
        outdir / f"experiments/GRP2/06_reproducibility/final/"
        f"GRP2.{assay}.macs3.narrow.replicate_validated.idr.narrowPeak"
    )
    grp2_cons = str(
        outdir / f"experiments/GRP2/06_reproducibility/consensus/"
        f"GRP2.{assay}.macs3.narrow.consensus.narrowPeak"
    )
    grp2_cons_final = str(
        outdir / f"experiments/GRP2/06_reproducibility/final/"
        f"GRP2.{assay}.macs3.narrow.replicate_validated.consensus.narrowPeak"
    )
    grp3_idr = str(
        outdir / f"experiments/GRP3/06_reproducibility/final/"
        f"GRP3.{assay}.macs3.narrow.replicate_validated.idr.narrowPeak"
    )
    grp3_cons = str(
        outdir / f"experiments/GRP3/06_reproducibility/consensus/"
        f"GRP3.{assay}.macs3.narrow.consensus.narrowPeak"
    )
    grp3_cons_final = str(
        outdir / f"experiments/GRP3/06_reproducibility/final/"
        f"GRP3.{assay}.macs3.narrow.replicate_validated.consensus.narrowPeak"
    )
    assert grp2_idr in dag_targets
    assert grp2_cons in dag_targets
    assert grp2_cons_final not in dag_targets

    for br in (1, 2, 3):
        assert (
            str(
                outdir / "experiments/GRP3/06_reproducibility/consensus/biorep_peaks/"
                f"GRP3.biorep{br}.{assay}.macs3.narrow_peaks.narrowPeak"
            )
            in dag_targets
        )
    assert grp3_cons in dag_targets
    assert grp3_cons_final in dag_targets
    assert grp3_idr not in dag_targets

    _assert_three_replicate_command(run_snakemake, config_path, grp3_cons, minimum=2)
    for wrong_final in (grp2_cons_final, grp3_idr):
        direct = run_snakemake(config_path, extra_args=[wrong_final])
        assert direct.returncode != 0, wrong_final

    if seacr:
        for exp in ("GRP2", "GRP3"):
            assert (
                str(
                    outdir / f"experiments/{exp}/06_reproducibility/final/"
                    f"{exp}.cuttag.seacr.stringent.replicate_validated.consensus.bed"
                )
                in dag_targets
            )
        _assert_three_replicate_command(
            run_snakemake,
            config_path,
            outdir / "experiments/GRP3/06_reproducibility/consensus/"
            "GRP3.cuttag.seacr.stringent.consensus.bed",
            minimum=2,
        )

    rows, _, _ = build_manifest_rows(cfg, samples_path=str(samples_path))
    manifest_proc, cli_rows = _run_manifest_cli(
        config_path, workdir / f"manifest_{idr_key}.tsv"
    )
    assert manifest_proc.returncode == 0, manifest_proc.stderr
    assert {(r["experiment_id"], r["output_type"], r["path"]) for r in rows} == {
        (r["experiment_id"], r["output_type"], r["path"]) for r in cli_rows
    }

    pairs = {(r["experiment_id"], r["output_type"]) for r in rows}
    assert ("GRP2", idr_final_type) in pairs
    assert ("GRP2", idr_summary_type) in pairs
    assert ("GRP2", cons_peak_type) in pairs
    assert ("GRP2", cons_summary_type) in pairs
    assert ("GRP2", cons_final_type) not in pairs

    assert ("GRP3", cons_peak_type) in pairs
    assert ("GRP3", cons_summary_type) in pairs
    assert ("GRP3", cons_final_type) in pairs
    assert ("GRP3", idr_final_type) not in pairs
    assert ("GRP3", idr_summary_type) not in pairs

    if seacr:
        assert ("GRP2", "cuttag_seacr_consensus_final_peak") in pairs
        assert ("GRP3", "cuttag_seacr_consensus_final_peak") in pairs

    expected_extract = {
        ("GRP2", idr_final_type),
        ("GRP3", cons_final_type),
    }
    if seacr:
        expected_extract |= {
            ("GRP2", "cuttag_seacr_consensus_final_peak"),
            ("GRP3", "cuttag_seacr_consensus_final_peak"),
        }
    final_types = {output_type for _, output_type in expected_extract}
    manifest_paths = {row["path"] for row in rows if row["output_type"] in final_types}
    assert manifest_paths == {
        path for path in dag_targets if ".replicate_validated." in Path(path).name
    }

    _materialize_present_outputs(rows, output_types=final_types)
    inline_inputs = WorkflowInputs(config=inline_cfg, samples=inline_samples)
    extracted = _extract_workspace_artifacts(adapter, inline_inputs, workdir)
    assert extracted.is_success, extracted.issues
    extracted_pairs = {
        (item.metadata.get("experiment_id"), item.output_type)
        for item in extracted.value
    }
    assert expected_extract <= extracted_pairs
    assert ("GRP2", cons_final_type) not in extracted_pairs
    assert ("GRP3", idr_final_type) not in extracted_pairs

    for item in extracted.value:
        _assert_final_track(item)


# ---------------------------------------------------------------------------
# 3. Global min_replicates=3 under per_group_fallback
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("pair_assay", "trio_assay", "seacr"),
    [
        ("chipseq", "chipseq", False),
        ("atac", "atac", False),
        ("cuttag", "cuttag", False),
        ("atac", "cuttag", True),
    ],
    ids=["chipseq", "atac", "cuttag", "atac-pair-cuttag-seacr-trio"],
)
def test_per_group_fallback_min3_skips_2rep_macs3_report_and_runs_3of3(
    tmp_path,
    tmp_config,
    run_validator,
    run_snakemake,
    pair_assay,
    trio_assay,
    seacr,
):
    workdir, config_path, samples_path, cfg, inline_samples = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name=f"min3_{pair_assay}_{trio_assay}",
        groups=[
            ("GRP2", pair_assay, "narrow", (1, 2), "PE"),
            ("GRP3", trio_assay, "narrow", (1, 2, 3), "PE"),
        ],
        extra_config={
            "cuttag": {"seacr": {"enabled": seacr, "mode": "stringent"}},
            "reproducibility": {
                "enabled": True,
                "mixed_replicate_policy": "per_group_fallback",
                "consensus": {"enabled": True, "min_replicates": 3},
                "idr": {f"{assay}_narrow": True for assay in (pair_assay, trio_assay)},
            },
        },
    )
    outdir = Path(cfg["outdir"])

    assert run_validator(config_path).returncode == 0
    assert _run_validate_cli(config_path).returncode == 0
    adapter = EncodeStyleWorkflowAdapter()
    assert adapter.validate(WorkflowInputs(config=cfg)).is_success
    inline_cfg = dict(cfg)
    inline_cfg.pop("samples", None)
    assert adapter.validate(
        WorkflowInputs(config=inline_cfg, samples=inline_samples)
    ).is_success

    sm = run_snakemake(config_path, extra_args=["--summary"])
    assert sm.returncode == 0, sm.stdout + sm.stderr
    dag_targets = _parse_summary(sm.stdout)

    # 2-rep IDR final remains; 2-rep auxiliary MACS3 consensus report is skipped.
    if pair_assay == "chipseq":
        assert (
            str(outdir / "experiments/GRP2/06_idr/final/conservative.narrowPeak")
            in dag_targets
        )
        assert (
            str(outdir / "experiments/GRP2/06_idr/final/optimal.narrowPeak")
            in dag_targets
        )
    else:
        assert (
            str(
                outdir / "experiments/GRP2/06_reproducibility/final/"
                f"GRP2.{pair_assay}.macs3.narrow.replicate_validated.idr.narrowPeak"
            )
            in dag_targets
        )
    assert (
        str(
            outdir / "experiments/GRP2/06_reproducibility/consensus/"
            f"GRP2.{pair_assay}.macs3.narrow.consensus.narrowPeak"
        )
        not in dag_targets
    )
    skipped_summary = (
        outdir / "experiments/GRP2/06_reproducibility/consensus/"
        f"GRP2.{pair_assay}.macs3.narrow.consensus.summary.tsv"
    )
    assert str(skipped_summary) not in dag_targets
    direct = run_snakemake(config_path, extra_args=[str(skipped_summary)])
    assert direct.returncode != 0, skipped_summary

    # 3-rep group runs 3-of-3 consensus and promotes primary final.
    grp3_final = (
        outdir / "experiments/GRP3/06_reproducibility/final/"
        f"GRP3.{trio_assay}.macs3.narrow.replicate_validated.consensus.narrowPeak"
    )
    assert str(grp3_final) in dag_targets
    _assert_three_replicate_command(run_snakemake, config_path, grp3_final, minimum=3)

    rows, _, _ = build_manifest_rows(cfg, samples_path=str(samples_path))
    manifest_proc, cli_rows = _run_manifest_cli(
        config_path, workdir / "manifest_min3.tsv"
    )
    assert manifest_proc.returncode == 0, manifest_proc.stderr
    assert {(r["experiment_id"], r["output_type"]) for r in rows} == {
        (r["experiment_id"], r["output_type"]) for r in cli_rows
    }

    pairs = {(r["experiment_id"], r["output_type"]) for r in rows}
    for suffix in ("peak", "summary", "final_peak"):
        assert ("GRP2", f"{pair_assay}_macs3_narrow_consensus_{suffix}") not in pairs
        assert ("GRP3", f"{trio_assay}_macs3_narrow_consensus_{suffix}") in pairs
    if seacr:
        assert ("GRP3", "cuttag_seacr_consensus_final_peak") in pairs
        seacr_final = (
            outdir / "experiments/GRP3/06_reproducibility/final/"
            "GRP3.cuttag.seacr.stringent.replicate_validated.consensus.bed"
        )
        assert str(seacr_final) in dag_targets
        _assert_three_replicate_command(
            run_snakemake, config_path, seacr_final, minimum=3
        )


@pytest.mark.parametrize("cuttag_idr", [True, False], ids=["idr-on", "idr-off"])
def test_per_group_fallback_min3_rejects_pe_cuttag_seacr_with_2_replicates(
    tmp_path, tmp_config, run_validator, run_snakemake, cuttag_idr
):
    workdir, config_path, samples_path, cfg, inline_samples = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name=f"min3_seacr_reject_{cuttag_idr}",
        groups=[
            ("GRP2", "cuttag", "narrow", (1, 2), "PE"),
            ("GRP3", "cuttag", "narrow", (1, 2, 3), "PE"),
            ("ATAC2", "atac", "narrow", (1, 2), "PE"),
        ],
        extra_config={
            "cuttag": {"seacr": {"enabled": True, "mode": "stringent"}},
            "reproducibility": {
                "enabled": True,
                "mixed_replicate_policy": "per_group_fallback",
                "consensus": {"enabled": True, "min_replicates": 3},
                "idr": {"cuttag_narrow": cuttag_idr, "atac_narrow": True},
            },
        },
    )

    _assert_all_public_entrypoints_reject(
        workdir=workdir,
        config_path=config_path,
        samples_path=samples_path,
        cfg=cfg,
        inline_samples=inline_samples,
        run_validator=run_validator,
        run_snakemake=run_snakemake,
    )


# ---------------------------------------------------------------------------
# 4. Bidirectional identity isolation & resumption for per_group_fallback
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "groups",
    [
        [("GRP2", "chipseq", "narrow", (1, 2), "PE")],
        [
            ("GRP2", "chipseq", "narrow", (1, 2), "PE"),
            ("GRP3", "chipseq", "narrow", (1, 2, 3), "PE"),
        ],
    ],
    ids=["pure2-same-samples", "mixed23"],
)
def test_per_group_fallback_identity_resumes_same_mode_and_isolates_from_strict(
    tmp_path, tmp_config, run_snakemake, snakefile, groups
):
    workdir, config_path, _, cfg, _ = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name="fb_identity",
        groups=groups,
        extra_config={
            "chipseq_idr": True,
            "reproducibility": {
                "enabled": True,
                "mixed_replicate_policy": "per_group_fallback",
                "consensus": {"enabled": True},
            },
        },
    )
    outdir = Path(cfg["outdir"])
    marker = outdir / ".helixweave/analysis_identity.json"
    eff = outdir / ".helixweave/effective_config.json"

    create = _run_snakemake_identity_targets(
        workdir, snakefile, config_path, marker, eff
    )
    assert create.returncode == 0, create.stderr
    marker_bytes = marker.read_bytes()
    eff_bytes = eff.read_bytes()

    # Same-mode resume: summary is ok / no update, re-exec does not rewrite bytes.
    sm_reuse = run_snakemake(
        config_path, extra_args=["--summary", str(marker), str(eff)]
    )
    assert sm_reuse.returncode == 0, sm_reuse.stderr
    reuse_summary = _parse_summary(sm_reuse.stdout)
    for rec in (marker, eff):
        assert reuse_summary[str(rec)]["status"] == "ok"
        assert reuse_summary[str(rec)]["plan"] == "no update"

    reexec = _run_snakemake_identity_targets(
        workdir, snakefile, config_path, marker, eff
    )
    assert reexec.returncode == 0, reexec.stderr
    assert marker.read_bytes() == marker_bytes
    assert eff.read_bytes() == eff_bytes

    # Negative control: changing threads=2 still marks identity/config params changed.
    cfg_t2 = dict(cfg)
    cfg_t2["threads"] = 2
    cfg_t2_path = workdir / "config_fb_t2.yaml"
    cfg_t2_path.write_text(yaml.safe_dump(cfg_t2, sort_keys=False), encoding="utf-8")
    sm_t2 = run_snakemake(cfg_t2_path, extra_args=["--summary", str(marker), str(eff)])
    assert sm_t2.returncode == 0, sm_t2.stderr
    t2_summary = _parse_summary(sm_t2.stdout)
    assert t2_summary[str(marker)]["status"] == "params changed"
    assert t2_summary[str(marker)]["plan"] == "update pending"

    # Bidirectional isolation: omitted and explicit strict must fail against fallback dir.
    for label, policy_dict in (
        ("omitted", {"enabled": True, "consensus": {"enabled": True}}),
        (
            "strict",
            {
                "enabled": True,
                "mixed_replicate_policy": "strict",
                "consensus": {"enabled": True},
            },
        ),
    ):
        cfg_other = dict(cfg)
        cfg_other["reproducibility"] = policy_dict
        other_path = workdir / f"config_{label}.yaml"
        other_path.write_text(
            yaml.safe_dump(cfg_other, sort_keys=False), encoding="utf-8"
        )
        refused = _run_snakemake_identity_targets(
            workdir, snakefile, other_path, marker, eff
        )
        assert refused.returncode != 0
        assert marker.read_bytes() == marker_bytes
        assert eff.read_bytes() == eff_bytes


# ---------------------------------------------------------------------------
# 5. Non-expansion & misconfiguration rejections under per_group_fallback
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("assay", "idr_key"),
    [
        ("chipseq", "chipseq_narrow"),
        ("atac", "atac_narrow"),
        ("cuttag", "cuttag_narrow"),
    ],
)
def test_per_group_fallback_rejects_pure_three_replicate_idr_without_silent_fallback(
    tmp_path, tmp_config, run_validator, run_snakemake, assay, idr_key
):
    workdir, config_path, samples_path, cfg, inline_samples = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name=f"pure3_reject_{idr_key}",
        groups=[("GRP3", assay, "narrow", (1, 2, 3), "PE")],
        extra_config={
            "reproducibility": {
                "enabled": True,
                "mixed_replicate_policy": "per_group_fallback",
                "consensus": {"enabled": True},
                "idr": {idr_key: True},
            }
        },
    )

    _assert_all_public_entrypoints_reject(
        workdir=workdir,
        config_path=config_path,
        samples_path=samples_path,
        cfg=cfg,
        inline_samples=inline_samples,
        run_validator=run_validator,
        run_snakemake=run_snakemake,
    )


@pytest.mark.parametrize(
    ("assay", "peak_mode", "idr_key", "other_bioreps"),
    [
        ("chipseq", "narrow", "chipseq_narrow", (1, 2, 3, 4)),
        ("atac", "narrow", "atac_narrow", (1, 2, 3, 4)),
        ("cuttag", "narrow", "cuttag_narrow", (1, 2, 3, 4)),
        ("chipseq", "broad", "chipseq_broad_experimental", (1, 2, 3)),
        ("cuttag", "broad", "cuttag_broad_experimental", (1, 2, 3)),
    ],
    ids=[
        "chipseq-narrow-4rep",
        "atac-narrow-4rep",
        "cuttag-narrow-4rep",
        "chipseq-broad-experimental-2+3",
        "cuttag-broad-experimental-2+3",
    ],
)
def test_per_group_fallback_does_not_expand_to_4plus_or_experimental_broad(
    tmp_path,
    tmp_config,
    run_validator,
    run_snakemake,
    assay,
    peak_mode,
    idr_key,
    other_bioreps,
):
    groups = [
        ("GRP2", assay, peak_mode, (1, 2), "PE"),
        ("GRP_OTHER", assay, peak_mode, other_bioreps, "PE"),
    ]
    idr_flags = {idr_key: True}
    if peak_mode == "broad":
        # Satisfy narrow IDR admission independently, so the broad group
        # actually reaches the unsupported replicate-count boundary.
        groups.append(("NARROW2", assay, "narrow", (1, 2), "PE"))
        idr_flags[f"{assay}_narrow"] = True
    workdir, config_path, samples_path, cfg, inline_samples = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name=f"no_expand_{idr_key}_{len(other_bioreps)}",
        groups=groups,
        extra_config={
            "reproducibility": {
                "enabled": True,
                "mixed_replicate_policy": "per_group_fallback",
                "consensus": {"enabled": True},
                "idr": idr_flags,
            }
        },
    )

    _assert_all_public_entrypoints_reject(
        workdir=workdir,
        config_path=config_path,
        samples_path=samples_path,
        cfg=cfg,
        inline_samples=inline_samples,
        run_validator=run_validator,
        run_snakemake=run_snakemake,
    )


@pytest.mark.parametrize(
    "repro_override",
    [
        {
            "enabled": False,
            "mixed_replicate_policy": "per_group_fallback",
            "consensus": {"enabled": True},
            "idr": {"atac_narrow": True},
        },
        {
            "enabled": True,
            "mixed_replicate_policy": "per_group_fallback",
            "consensus": {"enabled": True},
        },
        {
            "enabled": True,
            "mixed_replicate_policy": "per_group_fallback",
            "consensus": {"enabled": False},
            "idr": {"atac_narrow": True},
        },
        {
            "enabled": True,
            "mixed_replicate_policy": "invalid_policy",
            "consensus": {"enabled": True},
            "idr": {"atac_narrow": True},
        },
        {
            "enabled": True,
            "mixed_replicate_policy": "per_group_fallback",
            "consensus": {"enabled": True, "min_replicates": 4},
            "idr": {"atac_narrow": True},
        },
    ],
    ids=[
        "reproducibility-disabled",
        "no-narrow-idr",
        "consensus-disabled",
        "unknown-policy",
        "min-replicates-4",
    ],
)
def test_per_group_fallback_rejects_misconfiguration(
    tmp_path, tmp_config, run_validator, run_snakemake, repro_override
):
    workdir, config_path, samples_path, cfg, inline_samples = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name="misconfig",
        groups=[
            ("GRP2", "atac", "narrow", (1, 2), "PE"),
            ("GRP3", "atac", "narrow", (1, 2, 3), "PE"),
        ],
        extra_config={"reproducibility": repro_override},
    )

    _assert_all_public_entrypoints_reject(
        workdir=workdir,
        config_path=config_path,
        samples_path=samples_path,
        cfg=cfg,
        inline_samples=inline_samples,
        run_validator=run_validator,
        run_snakemake=run_snakemake,
    )


# ---------------------------------------------------------------------------
# 6. Lightweight real consensus execution for 3-rep fallback group
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("assay", "idr_key", "min_replicates"),
    [
        ("chipseq", "chipseq_narrow", 2),
        ("chipseq", "chipseq_narrow", 3),
        ("atac", "atac_narrow", 2),
        ("cuttag", "cuttag_narrow", 2),
    ],
)
def test_per_group_fallback_lightweight_real_consensus_execution(
    tmp_path,
    tmp_config,
    monkeypatch,
    snakemake_executable,
    snakefile,
    assay,
    idr_key,
    min_replicates,
):
    """Execute consensus_compute_narrow + consensus_final on tiny precomputed peaks."""
    bindir = tmp_path / "bin"
    bindir.mkdir()
    (bindir / "python3").symlink_to(sys.executable)
    monkeypatch.setenv("PATH", str(bindir) + os.pathsep + os.environ.get("PATH", ""))
    monkeypatch.setenv("PYTHONDONTWRITEBYTECODE", "1")
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / ".cache"))

    workdir, config_path, samples_path, cfg, _ = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name=f"exec_{assay}_min{min_replicates}",
        groups=[
            ("GRP2", assay, "narrow", (1, 2), "PE"),
            ("GRP3", assay, "narrow", (1, 2, 3), "PE"),
        ],
        extra_config={
            "trim": False,
            "multiqc": False,
            "qc": {
                "blacklist_filter": False,
                "frip": False,
                "library_complexity": False,
                "nrf_pbc": False,
                "signal_tracks": False,
                "summary": True,
                "cuttag_fragment_size": False,
            },
            "reproducibility": {
                "enabled": True,
                "mixed_replicate_policy": "per_group_fallback",
                "consensus": {"enabled": True, "min_replicates": min_replicates},
                "idr": {idr_key: True},
            },
        },
    )
    outdir = Path(cfg["outdir"])

    # Seed the real run-identity marker + effective_config before staging inputs.
    seed_run_identity(
        workdir, snakemake_executable, snakefile, config_path, outdir=str(outdir)
    )
    record = json.loads((outdir / ".helixweave/effective_config.json").read_bytes())
    routes = record["resolved_reproducibility_by_experiment"]
    assert {
        exp: (route["caller"], route["n_bioreps"], route["route"])
        for exp, route in routes.items()
    } == {"GRP2": ("macs3", 2, "idr"), "GRP3": ("macs3", 3, "consensus")}
    if min_replicates == 3:
        assert (
            routes["GRP2"]["auxiliary_consensus_skipped"]
            == "insufficient_biological_replicates"
        )
    else:
        assert "auxiliary_consensus_skipped" not in routes["GRP2"]
    assert "auxiliary_consensus_skipped" not in routes["GRP3"]

    # Stage tiny synthetic per-biorep narrowPeak inputs for GRP3 (3 bioreps):
    # - 3/3 peak at chr1:100-220 (all 3 bioreps, max signal 10)
    # - 2/3 peak at chr1:500-610 (requires bioreps 2 & 3, max signal 8)
    # - 1/3 singleton at chr1:900-1000 (biorep 3 only)
    biorep_lines = {
        1: "chr1\t100\t200\tp3_1\t100\t.\t5\t3\t2\t20\n",
        2: (
            "chr1\t120\t220\tp3_2\t100\t.\t10\t3\t2\t30\n"
            "chr1\t510\t610\tp2_2\t100\t.\t8\t3\t2\t35\n"
        ),
        3: (
            "chr1\t110\t210\tp3_3\t100\t.\t7\t3\t2\t15\n"
            "chr1\t500\t600\tp2_3\t100\t.\t4\t3\t2\t25\n"
            "chr1\t900\t1000\tp1_3\t100\t.\t6\t3\t2\t40\n"
        ),
    }
    staged_inputs = []
    for br, content in biorep_lines.items():
        peak_path = (
            outdir
            / "experiments/GRP3/06_reproducibility/consensus/biorep_peaks"
            / f"GRP3.biorep{br}.{assay}.macs3.narrow_peaks.narrowPeak"
        )
        peak_path.parent.mkdir(parents=True, exist_ok=True)
        peak_path.write_text(content, encoding="utf-8")
        staged_inputs.append(peak_path)

    cons_peak = (
        outdir
        / "experiments/GRP3/06_reproducibility/consensus"
        / f"GRP3.{assay}.macs3.narrow.consensus.narrowPeak"
    )
    cons_summary = (
        outdir
        / "experiments/GRP3/06_reproducibility/consensus"
        / f"GRP3.{assay}.macs3.narrow.consensus.summary.tsv"
    )
    cons_final = (
        outdir
        / "experiments/GRP3/06_reproducibility/final"
        / f"GRP3.{assay}.macs3.narrow.replicate_validated.consensus.narrowPeak"
    )

    exec_proc = subprocess.run(
        [
            snakemake_executable,
            "-s",
            snakefile,
            "--workflow-profile",
            "none",
            "--configfile",
            str(config_path),
            "--directory",
            str(workdir),
            "--cores",
            "1",
            "--drop-metadata",
            "--allowed-rules",
            "consensus_compute_narrow",
            "consensus_final",
            "--",
            str(cons_peak),
            str(cons_summary),
            str(cons_final),
        ],
        cwd=workdir,
        capture_output=True,
        text=True,
    )
    assert exec_proc.returncode == 0, exec_proc.stdout + exec_proc.stderr
    assert cons_final.read_bytes() == cons_peak.read_bytes()

    peak_lines = [
        line.split("\t")
        for line in cons_final.read_text(encoding="utf-8").splitlines()
        if line
    ]
    if min_replicates == 2:
        assert len(peak_lines) == 2
        assert peak_lines[0][:5] == ["chr1", "100", "220", "consensus_peak_1", "1000"]
        assert peak_lines[1][:5] == ["chr1", "500", "610", "consensus_peak_2", "667"]
        expected_support_dist = {"2": 1, "3": 1}
    else:
        assert len(peak_lines) == 1
        assert peak_lines[0][:5] == ["chr1", "100", "220", "consensus_peak_1", "1000"]
        expected_support_dist = {"3": 1}

    with cons_summary.open(encoding="utf-8", newline="") as handle:
        summary_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(summary_rows) == 1
    assert summary_rows[0]["experiment"] == "GRP3"
    assert summary_rows[0]["assay"] == assay
    assert summary_rows[0]["peak_mode"] == "narrow"
    assert summary_rows[0]["caller"] == "macs3"
    assert summary_rows[0]["n_bioreps"] == "3"
    assert summary_rows[0]["min_replicates"] == str(min_replicates)
    assert summary_rows[0]["biorep_labels"] == "1,2,3"
    assert json.loads(summary_rows[0]["support_distribution"]) == expected_support_dist
    assert json.loads(summary_rows[0]["source_peak_files"]) == [
        str(p) for p in staged_inputs
    ]
    assert summary_rows[0]["final_method"] == "consensus"
    assert Path(summary_rows[0]["final_output"]).resolve() == cons_final.resolve()

    rows, _, _ = build_manifest_rows(cfg, samples_path=str(samples_path))
    present_pairs = {
        (r["experiment_id"], r["output_type"]) for r in rows if r["status"] == "present"
    }
    expected_final_type = f"{assay}_macs3_narrow_consensus_final_peak"
    assert ("GRP3", expected_final_type) in present_pairs

    adapter = EncodeStyleWorkflowAdapter()
    extracted = _extract_workspace_artifacts(
        adapter, WorkflowInputs(config=cfg), workdir
    )
    assert extracted.is_success, extracted.issues
    extracted_final = [
        item
        for item in extracted.value
        if item.metadata.get("experiment_id") == "GRP3"
        and item.output_type == expected_final_type
    ]
    assert len(extracted_final) == 1
    assert workdir / extracted_final[0].relative_path == cons_final
    _assert_final_track(extracted_final[0])
