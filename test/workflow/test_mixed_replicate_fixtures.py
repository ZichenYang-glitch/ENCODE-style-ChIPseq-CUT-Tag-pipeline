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
from pathlib import Path

import pytest
import yaml

from conftest import seed_run_identity
from encode_pipeline.adapters.browser.catalog import describe_browser_artifact
from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.config.scientific_policy import SCIENTIFIC_POLICY_REVISION
from encode_pipeline.errors import ValidationError
from encode_pipeline.manifest.make import MANIFEST_COLUMNS, build_manifest_rows
from encode_pipeline.platform.adapters import WorkflowInputs


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
        for br in bioreps:
            fq2 = str(inputs["r2"]) if layout == "PE" else ""
            row = [
                f"{exp}_r{br}",
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
                "1",
            ]
            lines.append("\t".join(row))
    return "\n".join(lines) + "\n"


def _build_inline_samples(
    inputs: dict[str, Path],
    groups: list[tuple[str, str, str, tuple[int, ...], str]],
) -> list[dict[str, str]]:
    rows = []
    for exp, assay, peak_mode, bioreps, layout in groups:
        target = "H3K27me3" if peak_mode == "broad" else "CTCF"
        for br in bioreps:
            rows.append(
                {
                    "sample": f"{exp}_r{br}",
                    "layout": layout,
                    "fastq_1": str(inputs["r1"]),
                    "fastq_2": str(inputs["r2"]) if layout == "PE" else "",
                    "assay": assay,
                    "peak_mode": peak_mode,
                    "genome": "hs",
                    "bowtie2_index": str(inputs["bt2"]),
                    "target": target,
                    "experiment": exp,
                    "role": "treatment",
                    "biological_replicate": str(br),
                    "technical_replicate": "1",
                }
            )
    return rows


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
    inline_samples = _build_inline_samples(inputs, groups)
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
    (workdir / "config").mkdir(exist_ok=True)
    (workdir / "config" / "config.yaml").write_text(
        Path(config_path).read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    (workdir / "config" / "samples.tsv").write_text(
        Path(samples_path).read_text(encoding="utf-8"),
        encoding="utf-8",
    )
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


def _rule_shell_blocks(output: str) -> dict[str, list[str]]:
    blocks: dict[str, list[str]] = {}
    current: str | None = None
    buf: list[str] = []
    for line in output.splitlines():
        stripped = line.strip()
        if stripped.startswith(("rule ", "localrule ")) and stripped.endswith(":"):
            if current is not None:
                blocks.setdefault(current, []).append("\n".join(buf))
            current = stripped.split()[1].rstrip(":")
            buf = [line]
        elif current is not None:
            buf.append(line)
    if current is not None:
        blocks.setdefault(current, []).append("\n".join(buf))
    return blocks


def _scheduled_outputs_and_commands(run_snakemake, config_path):
    # Snakemake's summary exits before rendering any job shell commands.
    summary = run_snakemake(config_path, extra_args=["--summary"])
    assert summary.returncode == 0, summary.stdout + summary.stderr
    jobs = run_snakemake(config_path, extra_args=["--printshellcmds"])
    assert jobs.returncode == 0, jobs.stdout + jobs.stderr
    return _parse_summary(summary.stdout), _rule_shell_blocks(jobs.stdout + jobs.stderr)


def _materialize_present_outputs(
    manifest_rows: list[dict],
    *,
    output_types: set[str] | None = None,
) -> None:
    for row in manifest_rows:
        if output_types is not None and row["output_type"] not in output_types:
            continue
        path = Path(row["path"])
        if not row["path"] or path.suffix == "":
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("chr1\t10\t20\tp1\t1000\t.\t1.0\t-1\t-1\t5\n", encoding="utf-8")


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
            assert reader.fieldnames == MANIFEST_COLUMNS
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
    cross_mode = run_snakemake(
        cfg_fallback_path, extra_args=["--summary", str(marker), str(eff)]
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
    extracted = adapter.extract_artifacts(WorkflowInputs(config=cfg), workdir)
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
    "assay,idr_setting,minimum,seacr",
    [
        ("chipseq", "chipseq_idr", 2, False),
        ("chipseq", "chipseq_narrow", 2, False),
        ("atac", "atac_narrow", 2, False),
        ("cuttag", "cuttag_narrow", 2, False),
        ("cuttag", "cuttag_narrow", 2, True),
        ("chipseq", "chipseq_narrow", 3, False),
        ("atac", "atac_narrow", 3, False),
        ("cuttag", "cuttag_narrow", 3, False),
    ],
)
def test_mixed_pipeline_assigns_each_published_final_to_its_group(
    tmp_path,
    tmp_config,
    run_validator,
    run_snakemake,
    assay,
    idr_setting,
    minimum,
    seacr,
):
    repro = {
        "enabled": True,
        "mixed_replicate_policy": "per_group_fallback",
        "consensus": {"enabled": True, "min_replicates": minimum},
    }
    extra = {"reproducibility": repro}
    if idr_setting == "chipseq_idr":
        extra["chipseq_idr"] = True
    else:
        repro["idr"] = {idr_setting: True}
    if seacr:
        extra["cuttag"] = {"seacr": {"enabled": True, "mode": "stringent"}}
    work, config_path, samples_path, cfg, inline = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name=assay,
        groups=[
            ("GRP2", assay, "narrow", (1, 2), "PE"),
            ("GRP3", assay, "narrow", (1, 2, 3), "PE"),
        ],
        extra_config=extra,
    )
    validation = run_validator(config_path)
    assert validation.returncode == 0, validation.stderr
    cli = _run_validate_cli(config_path)
    assert cli.returncode == 0, cli.stderr
    assert "GRP2 (2 biological replicates, idr" in cli.stdout + cli.stderr
    assert "GRP3 (3 biological replicates, consensus" in cli.stdout + cli.stderr
    adapter = EncodeStyleWorkflowAdapter()
    assert adapter.validate(WorkflowInputs(config=cfg)).is_success
    inline_cfg = dict(cfg)
    inline_cfg.pop("samples")
    assert adapter.validate(
        WorkflowInputs(config=inline_cfg, samples=inline)
    ).is_success
    plan = adapter.plan_workspace(WorkflowInputs(config=cfg), work)
    assert plan.is_success, plan.issues
    persisted = yaml.safe_load(dict(plan.value.files)["config/config.yaml"])
    assert (
        persisted["reproducibility"]["mixed_replicate_policy"] == "per_group_fallback"
    )
    assert persisted["reproducibility"]["consensus"]["min_replicates"] == minimum

    dag, blocks = _scheduled_outputs_and_commands(run_snakemake, config_path)
    for rule in (
        ["consensus_compute_narrow", "consensus_compute_seacr"]
        if seacr
        else ["consensus_compute_narrow"]
    ):
        assert any(
            "GRP3" in job
            and "--bioreps 1 2 3" in job
            and f"--min-replicates {minimum}" in job
            for job in blocks[rule]
        )
    root = Path(cfg["outdir"]) / "experiments"
    for br in (1, 2, 3):
        peak = (
            root
            / "GRP3/06_reproducibility/consensus/biorep_peaks"
            / f"GRP3.biorep{br}.{assay}.macs3.narrow_peaks.narrowPeak"
        )
        assert str(peak) in dag

    rows, _, _ = build_manifest_rows(cfg, samples_path=str(samples_path))
    proc, cli_rows = _run_manifest_cli(config_path, work / "manifest.tsv")
    assert proc.returncode == 0, proc.stderr
    assert cli_rows == rows
    by_group = {(r["experiment_id"], r["output_type"]): r["path"] for r in rows}
    prefix = f"{assay}_macs3_narrow"
    final_types = {f"{prefix}_consensus_final_peak"}
    expected = {("GRP3", f"{prefix}_consensus_final_peak")}
    if assay == "chipseq":
        idr_types = {"idr_conservative", "idr_optimal"}
        summary_type = "idr_reproducibility_summary"
    else:
        idr_types = {f"{prefix}_idr_final_peak"}
        summary_type = f"{prefix}_idr_summary"
    final_types.update(idr_types)
    expected.update(("GRP2", t) for t in idr_types)
    assert ("GRP2", summary_type) in by_group
    assert ("GRP3", summary_type) not in by_group
    if seacr:
        final_types.add("cuttag_seacr_consensus_final_peak")
        expected.update(
            (group, "cuttag_seacr_consensus_final_peak") for group in ("GRP2", "GRP3")
        )
    assert {key for key in by_group if key[1] in final_types} == expected

    for group, scheduled in (("GRP2", minimum == 2), ("GRP3", True)):
        for kind, suffix in (("peak", "narrowPeak"), ("summary", "summary.tsv")):
            path = (
                root
                / group
                / "06_reproducibility/consensus"
                / f"{group}.{assay}.macs3.narrow.consensus.{suffix}"
            )
            assert (str(path) in dag) == scheduled
            assert ((group, f"{prefix}_consensus_{kind}") in by_group) == scheduled
    actual_finals = {
        p for p in dag if "/final/" in p and p.endswith((".narrowPeak", ".bed"))
    }
    assert actual_finals == {by_group[key] for key in expected}
    for group, kind in expected:
        path = Path(by_group[group, kind])
        assert path.is_relative_to(root / group)
        if kind.startswith("idr_"):
            assert path.name == kind.removeprefix("idr_") + ".narrowPeak"
        else:
            method = "idr" if kind.endswith("idr_final_peak") else "consensus"
            assert f".replicate_validated.{method}." in path.name
    # A directly requested final must respect the same experiment routing as rule all.
    trio_final = root / "GRP3/06_reproducibility/final"
    forbidden = [
        root
        / "GRP2/06_reproducibility/final"
        / f"GRP2.{assay}.macs3.narrow.replicate_validated.consensus.narrowPeak",
        root / "GRP3/06_idr/final/optimal.narrowPeak"
        if assay == "chipseq"
        else trio_final
        / f"GRP3.{assay}.macs3.narrow.replicate_validated.idr.narrowPeak",
    ]
    if minimum == 3:
        forbidden.append(
            root
            / "GRP2/06_reproducibility/consensus"
            / f"GRP2.{assay}.macs3.narrow.consensus.narrowPeak"
        )
    for target in forbidden:
        result = run_snakemake(config_path, extra_args=[str(target)])
        assert result.returncode != 0, f"incorrectly scheduled {target}"
        assert not target.exists()

    _materialize_present_outputs(rows, output_types=final_types)
    extracted = adapter.extract_artifacts(WorkflowInputs(config=cfg), work)
    assert extracted.is_success, extracted.issues
    final_artifacts = [
        item for item in extracted.value if item.output_type in final_types
    ]
    assert {
        (item.metadata["experiment_id"], item.output_type) for item in final_artifacts
    } == expected
    for item in final_artifacts:
        desc = describe_browser_artifact(
            workflow_id=WORKFLOW_ID,
            output_type=item.output_type,
            name=Path(item.relative_path).name,
        )
        assert desc is not None and desc.available
        assert (desc.wash_u_type, desc.role) == ("bed", "interval")
        assert desc.label and desc.detail and desc.score_meaning


def test_named_consensus_finals_match_default_dag_and_manifest(
    tmp_path, tmp_config, run_snakemake, snakemake_executable, snakefile
):
    work, config_path, samples_path, cfg, _ = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name="named_finals",
        groups=[
            ("CHIP2", "chipseq", "narrow", (1, 2), "PE"),
            ("ATAC3", "atac", "narrow", (1, 2, 3), "PE"),
            ("CUT3", "cuttag", "narrow", (1, 2, 3), "PE"),
            ("CUT1", "cuttag", "narrow", (1,), "PE"),
        ],
        extra_config={
            "chipseq_idr": True,
            "reproducibility": {
                "enabled": True,
                "mixed_replicate_policy": "per_group_fallback",
                "consensus": {"enabled": True, "min_replicates": 2},
                "idr": {"atac_narrow": False, "cuttag_narrow": False},
            },
            "cuttag": {"seacr": {"enabled": True, "mode": "stringent"}},
        },
    )
    dag = run_snakemake(config_path, extra_args=["--summary"])
    assert dag.returncode == 0, dag.stdout + dag.stderr
    scheduled = _parse_summary(dag.stdout)
    rows, _, _ = build_manifest_rows(cfg, samples_path=str(samples_path))
    manifest_paths = {row["path"] for row in rows}
    outdir = Path(cfg["outdir"])
    seed_run_identity(work, snakemake_executable, snakefile, config_path, str(outdir))

    peaks = {
        "bed": b"chr1\t100\t200\tpeak1\t667\t.\n",
        "broadPeak": b"chr1\t100\t200\tpeak1\t667\t.\t5\t-1\t-1\n",
        "narrowPeak": b"chr1\t100\t200\tpeak1\t667\t.\t5\t-1\t-1\t50\n",
    }
    ineligible = "not eligible for a MACS3 consensus final"
    for group, assay_caller_mode, suffix, rejection in [
        ("ATAC3", "atac.macs3.narrow", "narrowPeak", ineligible),
        ("CUT3", "cuttag.macs3.narrow", "narrowPeak", None),
        ("CUT1", "cuttag.macs3.narrow", "narrowPeak", ineligible),
        ("CHIP2", "chipseq.macs3.narrow", "narrowPeak", ineligible),
        ("CUT3", "atac.macs3.narrow", "narrowPeak", "assay wildcard mismatch"),
        ("CUT3", "cuttag.macs3.narrow", "broadPeak", "suffix wildcard mismatch"),
        ("CUT3", "cuttag.seacr.stringent", "bed", None),
        ("CUT3", "cuttag.seacr.relaxed", "bed", "SEACR consensus mode mismatch"),
        ("ATAC3", "cuttag.seacr.stringent", "bed", "not eligible for SEACR consensus"),
    ]:
        publishes = rejection is None
        base = outdir / "experiments" / group / "06_reproducibility"
        stem = f"{group}.{assay_caller_mode}"
        source = base / "consensus" / f"{stem}.consensus.{suffix}"
        final = base / "final" / f"{stem}.replicate_validated.consensus.{suffix}"
        assert (str(final) in scheduled) == publishes, final
        assert (str(final) in manifest_paths) == publishes, final
        # Stage each source so a named request exercises the publication rule itself.
        source.parent.mkdir(parents=True, exist_ok=True)
        peak = peaks[suffix]
        source.write_bytes(peak)
        result = subprocess.run(
            [
                snakemake_executable,
                "-s",
                snakefile,
                "--configfile",
                str(config_path),
                "--workflow-profile",
                "none",
                "--cores",
                "1",
                "--allowed-rules",
                "consensus_final",
                "consensus_final_seacr",
                "--",
                str(final),
            ],
            cwd=work,
            capture_output=True,
            text=True,
        )
        output = result.stdout + result.stderr
        assert source.read_bytes() == peak
        if rejection is None:
            assert result.returncode == 0, output
            assert final.read_bytes() == peak
        else:
            assert result.returncode != 0, final
            assert rejection in output
            assert not final.exists(), final


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
        ],
        extra_config={
            "cuttag": {"seacr": {"enabled": True, "mode": "stringent"}},
            "reproducibility": {
                "enabled": True,
                "mixed_replicate_policy": "per_group_fallback",
                "consensus": {"enabled": True, "min_replicates": 3},
                "idr": {"cuttag_narrow": cuttag_idr},
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

    # Launch against the same outdir from a fresh project config. Snakemake
    # merges CLI overrides onto the launch directory's config/config.yaml.
    launch = workdir / "independent-launch"
    (launch / "config").mkdir(parents=True)
    (launch / "config/config.yaml").write_text("{}\n")
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
        other_path = launch / f"config_{label}.yaml"
        other_path.write_text(
            yaml.safe_dump(cfg_other, sort_keys=False), encoding="utf-8"
        )
        refused = run_snakemake(
            other_path, extra_args=["--summary", str(marker), str(eff)]
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
    workdir, config_path, samples_path, cfg, inline_samples = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name=f"no_expand_{idr_key}_{len(other_bioreps)}",
        groups=[
            ("GRP2", assay, peak_mode, (1, 2), "PE"),
            ("GRP_OTHER", assay, peak_mode, other_bioreps, "PE"),
        ],
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
    "repro_override",
    [
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
    ids=["consensus-disabled", "unknown-policy", "min-replicates-4"],
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

    record = json.loads((outdir / ".helixweave/effective_config.json").read_text())
    routes = record["resolved_reproducibility_by_experiment"]
    assert routes["GRP3"] == {"route": "consensus", "caller": "macs3", "n_bioreps": 3}
    expected_pair = {"route": "idr", "caller": "macs3", "n_bioreps": 2}
    if min_replicates == 3:
        expected_pair["auxiliary_consensus_skipped"] = (
            "insufficient_biological_replicates"
        )
    assert routes["GRP2"] == expected_pair

    # The two-way hit occurs in reps 2 and 3. Overlapping singleton peaks
    # in rep 1 still count only once, even though their coordinates differ.
    biorep_lines = {
        1: (
            "chr1\t100\t200\tp3_1\t100\t.\t5\t3\t2\t20\n"
            "chr1\t900\t1000\tp1a\t100\t.\t6\t3\t2\t40\n"
            "chr1\t905\t1005\tp1b\t100\t.\t6\t3\t2\t40\n"
        ),
        2: (
            "chr1\t120\t220\tp3_2\t100\t.\t10\t3\t2\t30\n"
            "chr1\t510\t610\tp2_2\t100\t.\t8\t3\t2\t35\n"
        ),
        3: (
            "chr1\t110\t210\tp3_3\t100\t.\t7\t3\t2\t15\n"
            "chr1\t500\t600\tp2_3\t100\t.\t4\t3\t2\t25\n"
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

    rules = ["consensus_compute_narrow", "consensus_final"]
    targets = [cons_peak, cons_summary, cons_final]
    idr = outdir / "experiments/GRP2/06_idr"
    if assay == "chipseq":
        # Start the pair at precomputed, thresholded IDR peaks and run the
        # original summarizer in the same invocation as the three-rep final.
        for path, peaks in {
            "true_replicates/idr.thresholded.narrowPeak": biorep_lines[2],
            "pooled_pseudoreplicates/idr.thresholded.narrowPeak": biorep_lines[3],
            "self_pseudoreplicates/biorep1.idr.thresholded.narrowPeak": biorep_lines[2],
            "self_pseudoreplicates/biorep2.idr.thresholded.narrowPeak": biorep_lines[3],
        }.items():
            (idr / path).parent.mkdir(parents=True, exist_ok=True)
            (idr / path).write_text(peaks)
        rules.append("chipseq_idr_summary")
        targets.extend(
            [idr / "final/conservative.narrowPeak", idr / "final/optimal.narrowPeak"]
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
            *rules,
            "--",
            *map(str, targets),
        ],
        cwd=workdir,
        capture_output=True,
        text=True,
    )
    assert exec_proc.returncode == 0, exec_proc.stdout + exec_proc.stderr
    if assay == "chipseq":
        assert (idr / "final/conservative.narrowPeak").read_text() == biorep_lines[2]
        assert (idr / "final/optimal.narrowPeak").read_text() == biorep_lines[3]
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
    extracted = adapter.extract_artifacts(WorkflowInputs(config=cfg), workdir)
    assert extracted.is_success, extracted.issues
    extracted_final = [
        item
        for item in extracted.value
        if item.metadata.get("experiment_id") == "GRP3"
        and item.output_type == expected_final_type
    ]
    assert len(extracted_final) == 1
    desc = describe_browser_artifact(
        workflow_id=WORKFLOW_ID,
        output_type=extracted_final[0].output_type,
        name=Path(extracted_final[0].relative_path).name,
    )
    assert desc is not None
    assert desc.available is True
    assert desc.wash_u_type == "bed"
    assert desc.role == "interval"
