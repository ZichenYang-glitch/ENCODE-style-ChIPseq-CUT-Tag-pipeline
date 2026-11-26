"""Executable synthetic fixtures and verification for strict and candidate 2/3 routing.

Follows IRENE-51 / Draft PR #262 without changing production scientific policy
or defaults (SCIENTIFIC_POLICY_REVISION stays 7, default admission stays strict).
"""

from __future__ import annotations

import json
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

from encode_pipeline.adapters.browser.catalog import describe_browser_artifact
from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.config.scientific_policy import SCIENTIFIC_POLICY_REVISION
from encode_pipeline.errors import ValidationError
from encode_pipeline.manifest.make import build_manifest_rows
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.samples import replicates
from encode_pipeline.samples.replicates import validate_replicate_groups


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


def _build_fixture_workspace(
    tmp_path: Path,
    tmp_config,
    *,
    name: str,
    groups: list[tuple[str, str, str, tuple[int, ...], str]],
    extra_config: dict,
) -> tuple[Path, Path, Path, dict]:
    inputs = _write_synthetic_inputs(tmp_path / f"inp_{name}")
    samples_text = _render_samples_tsv(inputs, groups)
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
    return workdir, Path(config_path), Path(samples_path), loaded_cfg


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


# ---------------------------------------------------------------------------
# 1. Strict-policy 5-scenario comparison (baseline contract + admission parity)
# ---------------------------------------------------------------------------


def test_strict_policy_pure_two_chipseq_idr_and_non_rewrite_identity(
    tmp_path, tmp_config, run_validator, run_snakemake, snakefile
):
    assert SCIENTIFIC_POLICY_REVISION == 7
    workdir, config_path, samples_path, cfg = _build_fixture_workspace(
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
    create = subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "-s",
            snakefile,
            "--configfile",
            str(config_path),
            "--cores",
            "1",
            str(marker),
            str(eff),
        ],
        cwd=workdir,
        capture_output=True,
        text=True,
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

    reexec = subprocess.run(
        [
            sys.executable,
            "-m",
            "snakemake",
            "-s",
            snakefile,
            "--configfile",
            str(config_path),
            "--cores",
            "1",
            str(marker),
            str(eff),
        ],
        cwd=workdir,
        capture_output=True,
        text=True,
    )
    assert reexec.returncode == 0, reexec.stderr
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
    workdir, config_path, samples_path, cfg = _build_fixture_workspace(
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
    workdir, config_path, samples_path, cfg = _build_fixture_workspace(
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
    tmp_path, tmp_config, run_validator, run_snakemake, assay, peak_mode, idr_key
):
    workdir, config_path, samples_path, cfg = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name=f"reject_{idr_key}",
        groups=[
            ("GRP2", assay, peak_mode, (1, 2), "PE"),
            ("GRP3", assay, peak_mode, (1, 2, 3), "PE"),
        ],
        extra_config={
            "reproducibility": {
                "enabled": True,
                "consensus": {"enabled": True},
                "idr": {idr_key: True},
            }
        },
    )

    val = run_validator(config_path)
    assert val.returncode == 1
    assert "requires exactly 2" in val.stderr

    sm = run_snakemake(config_path, extra_args=["--summary"])
    assert sm.returncode == 1

    adapter = EncodeStyleWorkflowAdapter()
    path_res = adapter.validate(WorkflowInputs(config=cfg))
    assert not path_res.is_success
    assert path_res.issues[0].code == "ENCODE_SAMPLES_INVALID"

    inline_cfg = dict(cfg)
    inline_cfg.pop("samples", None)
    inline_res = adapter.validate(
        WorkflowInputs(
            config=inline_cfg,
            samples=[
                {
                    "sample": f"GRP{n}_{r}",
                    "layout": "PE",
                    "fastq_1": str(workdir / "r1.fq.gz"),
                    "fastq_2": str(workdir / "r2.fq.gz"),
                    "assay": assay,
                    "peak_mode": peak_mode,
                    "genome": "hs",
                    "bowtie2_index": cfg["bowtie2_index"],
                    "target": "H3K27me3" if peak_mode == "broad" else "CTCF",
                    "experiment": f"GRP{n}",
                    "role": "treatment",
                    "biological_replicate": str(r),
                    "technical_replicate": "1",
                }
                for n, reps in ((2, (1, 2)), (3, (1, 2, 3)))
                for r in reps
            ],
        )
    )
    assert not inline_res.is_success
    assert inline_res.issues[0].message == "Sample sheet is invalid."

    with pytest.raises(ValidationError, match="requires exactly 2"):
        build_manifest_rows(cfg, samples_path=str(samples_path))


def test_cuttag_independent_seacr_coexists_with_macs3_idr_d00_m00(
    tmp_path, tmp_config, run_validator, run_snakemake
):
    workdir, config_path, samples_path, cfg = _build_fixture_workspace(
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
        outdir
        / "experiments/CT2/06_reproducibility/final/CT2.cuttag.macs3.narrow.replicate_validated.idr.narrowPeak"
    )
    seacr_final = (
        outdir
        / "experiments/CT2/06_reproducibility/final/CT2.cuttag.seacr.stringent.replicate_validated.consensus.bed"
    )
    macs3_cons_final = (
        outdir
        / "experiments/CT2/06_reproducibility/final/CT2.cuttag.macs3.narrow.replicate_validated.consensus.narrowPeak"
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
            workflow_id="encode-style-chipseq-cuttag-atac-mnase",
            output_type=item.output_type,
            name=Path(item.relative_path).name,
        )
        assert desc is not None
        assert desc.wash_u_type == "bed"
        assert desc.role == "interval"


# ---------------------------------------------------------------------------
# 2. Candidate per-group routing verification (test-scoped acceptance evidence)
# ---------------------------------------------------------------------------


def _validate_candidate_groups(samples, use_control, *idr_args, **idr_kwargs):
    """Test only: check every group normally and apply strict IDR to the two-replicate group."""
    validate_replicate_groups(samples, use_control)
    bioreps = defaultdict(set)
    for sample in samples:
        if sample["role"] == "treatment":
            bioreps[sample["experiment"]].add(sample["biological_replicate"])
    assert sorted(map(len, bioreps.values())) == [2, 3]
    idr_samples = [
        row
        for row in samples
        if row["role"] == "control" or len(bioreps[row["experiment"]]) == 2
    ]
    validate_replicate_groups(idr_samples, use_control, *idr_args, **idr_kwargs)


@pytest.mark.parametrize(
    (
        "assay",
        "peak_mode",
        "idr_key",
        "idr_final_type",
        "consensus_final_type",
        "seacr",
    ),
    [
        (
            "cuttag",
            "narrow",
            "cuttag_narrow",
            "cuttag_macs3_narrow_idr_final_peak",
            "cuttag_macs3_narrow_consensus_final_peak",
            True,
        ),
        (
            "chipseq",
            "broad",
            "chipseq_broad_experimental",
            "chipseq_macs3_broad_idr_final_peak",
            "chipseq_macs3_broad_consensus_final_peak",
            False,
        ),
        (
            "cuttag",
            "broad",
            "cuttag_broad_experimental",
            "cuttag_macs3_broad_idr_final_peak",
            "cuttag_macs3_broad_consensus_final_peak",
            True,
        ),
    ],
    ids=["cuttag-narrow-seacr", "chipseq-broad", "cuttag-broad-seacr"],
)
def test_candidate_per_group_routing_dag_manifest_and_readonly_consumers(
    tmp_path,
    tmp_config,
    snakefile,
    run_validator,
    monkeypatch,
    assay,
    peak_mode,
    idr_key,
    idr_final_type,
    consensus_final_type,
    seacr,
):
    """Verify candidate per-group routing across DAG targets, consensus rules,
    manifest rows, artifact extraction, and WashU browser descriptions without
    altering production admission defaults.
    """
    extra_cfg = {
        "reproducibility": {
            "enabled": True,
            "consensus": {"enabled": True, "min_replicates": 2},
            "idr": {idr_key: True},
        }
    }
    if seacr:
        extra_cfg["cuttag"] = {"seacr": {"enabled": True, "mode": "stringent"}}

    workdir, config_path, samples_path, cfg = _build_fixture_workspace(
        tmp_path,
        tmp_config,
        name=f"cand_{idr_key}",
        groups=[
            ("GRP2", assay, peak_mode, (1, 2), "PE"),
            ("GRP3", assay, peak_mode, (1, 2, 3), "PE"),
        ],
        extra_config=extra_cfg,
    )
    outdir = cfg["outdir"]
    suffix = "narrowPeak" if peak_mode == "narrow" else "broadPeak"

    # Patch only the fixture's replicate admission. The real workflow derives its
    # experiment lists, expands the final targets and resolves their producing rules.
    monkeypatch.setattr(
        replicates, "validate_replicate_groups", _validate_candidate_groups
    )
    assert run_validator(config_path).returncode == 0
    sm = subprocess.run(
        [
            sys.executable,
            str(Path(__file__).resolve()),
            "--snakefile",
            snakefile,
            "--configfile",
            str(config_path),
            "--profile",
            "none",
            "--workflow-profile",
            "none",
        ],
        cwd=workdir,
        text=True,
        capture_output=True,
    )
    assert sm.returncode == 0, sm.stdout + sm.stderr
    dag_targets = _parse_summary(sm.stdout)

    grp2_idr = (
        f"{outdir}/experiments/GRP2/06_reproducibility/final/"
        f"GRP2.{assay}.macs3.{peak_mode}.replicate_validated.idr.{suffix}"
    )
    grp2_cons = (
        f"{outdir}/experiments/GRP2/06_reproducibility/final/"
        f"GRP2.{assay}.macs3.{peak_mode}.replicate_validated.consensus.{suffix}"
    )
    grp3_idr = (
        f"{outdir}/experiments/GRP3/06_reproducibility/final/"
        f"GRP3.{assay}.macs3.{peak_mode}.replicate_validated.idr.{suffix}"
    )
    grp3_cons = (
        f"{outdir}/experiments/GRP3/06_reproducibility/final/"
        f"GRP3.{assay}.macs3.{peak_mode}.replicate_validated.consensus.{suffix}"
    )
    assert grp2_idr in dag_targets
    assert grp2_cons not in dag_targets
    assert grp3_cons in dag_targets
    assert grp3_idr not in dag_targets
    if seacr:
        for exp in ("GRP2", "GRP3"):
            assert (
                f"{outdir}/experiments/{exp}/06_reproducibility/final/"
                f"{exp}.cuttag.seacr.stringent.replicate_validated.consensus.bed"
                in dag_targets
            )

    rows, _, _ = build_manifest_rows(cfg, samples_path=str(samples_path))
    pairs = {(r["experiment_id"], r["output_type"]) for r in rows}
    assert ("GRP2", idr_final_type) in pairs
    assert ("GRP2", consensus_final_type) not in pairs
    assert ("GRP3", consensus_final_type) in pairs
    assert ("GRP3", idr_final_type) not in pairs
    if seacr:
        assert ("GRP2", "cuttag_seacr_consensus_final_peak") in pairs
        assert ("GRP3", "cuttag_seacr_consensus_final_peak") in pairs

    expected_extract = {
        ("GRP2", idr_final_type),
        ("GRP3", consensus_final_type),
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

    adapter = EncodeStyleWorkflowAdapter()
    extracted = adapter.extract_artifacts(WorkflowInputs(config=cfg), workdir)
    assert extracted.is_success, extracted.issues
    extracted_pairs = {
        (item.metadata.get("experiment_id"), item.output_type)
        for item in extracted.value
    }
    assert expected_extract <= extracted_pairs
    assert ("GRP2", consensus_final_type) not in extracted_pairs
    assert ("GRP3", idr_final_type) not in extracted_pairs

    for item in extracted.value:
        desc = describe_browser_artifact(
            workflow_id="encode-style-chipseq-cuttag-atac-mnase",
            output_type=item.output_type,
            name=Path(item.relative_path).name,
        )
        assert desc is not None
        assert desc.wash_u_type == "bed"
        assert desc.role == "interval"


if __name__ == "__main__":
    # Snakemake runs in its own interpreter. Apply the same test-only admission
    # there while loading the unmodified production Snakefile from its checkout.
    from snakemake.cli import main

    with patch.object(
        replicates, "validate_replicate_groups", _validate_candidate_groups
    ):
        main(["--dry-run", "--summary", *sys.argv[1:]])
