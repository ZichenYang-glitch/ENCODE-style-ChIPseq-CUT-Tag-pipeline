"""Production-rule regression coverage for empty optional QC flags.

Snakemake's ``{params.name:q}`` drops an empty value entirely.  Before SCI-5
R1, ``analysis_policy`` and ``qc_summary`` always rendered their optional flags,
so the next token (``--output`` / ``--mnase-qc``) was consumed as the flag's
value and the production rules failed for default configurations.

These tests execute the real Snakemake rules with optional inputs enabled,
disabled, and not applicable to the assay.
"""

from __future__ import annotations

import csv
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from _tool_resolver import require_external_tools, resolve_tool
from conftest import seed_run_identity

from encode_pipeline.qc.contract import (
    ANALYSIS_POLICY_INPUT_STAGE_KEYS,
    FRAGMENT_QC_V2_COLUMNS,
    MITO_QC_V2_COLUMNS,
    MNASE_QC_V2_COLUMNS,
    POLICY_STAGE_COUNT_KEYS,
    SAMPLE_QC_V2_COLUMNS,
)

pytestmark = pytest.mark.real_execution

REPO = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO / "workflow" / "Snakefile"
SAMPLE = "S1"

SAMPLE_COLUMNS = (
    "sample",
    "fastq_1",
    "fastq_2",
    "layout",
    "assay",
    "target",
    "peak_mode",
    "genome",
    "bowtie2_index",
    "experiment",
    "biological_replicate",
)

PEAK_COLUMNS = ("sample", "peak_mode", "peaks", "blacklist_filtered_peaks")
FRIP_COLUMNS = ("sample", "total_reads", "reads_in_peaks", "frip", "bam", "peaks")
LIBRARY_COLUMNS = (
    "sample",
    "metrics_source",
    "unpaired_reads_examined",
    "read_pairs_examined",
    "secondary_or_supplementary_reads",
    "unmapped_reads",
    "unpaired_read_duplicates",
    "read_pair_duplicates",
    "read_pair_optical_duplicates",
    "percent_duplication",
    "estimated_library_size",
    "total_reads_examined",
    "duplicate_reads_estimate",
)
NRF_COLUMNS = (
    "sample",
    "total_fragments",
    "distinct_fragments",
    "one_read_fragments",
    "two_read_fragments",
    "nrf",
    "pbc1",
    "pbc2",
)


@pytest.fixture(scope="module")
def snakemake_executable() -> str:
    candidate = resolve_tool("snakemake", "SNAKEMAKE")
    resolved = shutil.which(candidate)
    if not resolved and os.path.isabs(candidate) and os.access(candidate, os.X_OK):
        resolved = candidate
    if not resolved:
        require_external_tools(["snakemake"], "QC optional-argv rule execution")
    return resolved


def _run(argv, cwd: Path):
    result = subprocess.run(
        [str(value) for value in argv],
        cwd=cwd,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"Command failed: {argv!r}\n"
        f"stdout:\n{result.stdout[-3000:]}\n"
        f"stderr:\n{result.stderr[-3000:]}"
    )
    return result


def _write_tsv(path: Path, columns, row) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerow({column: row.get(column, "NA") for column in columns})


def _write_project(
    tmp_path: Path,
    *,
    assay: str,
    cuttag_fragment_size: bool,
    cross_correlation: bool,
    outdir: str = "results",
    profile: str = "legacy",
) -> Path:
    target = {
        "chipseq": ("H3K27ac", "narrow"),
        "cuttag": ("H3K4me3", "narrow"),
        "mnase": ("Nucleosome", "nucleosome"),
    }[assay]
    samples = tmp_path / "config" / "samples.tsv"
    samples.parent.mkdir(parents=True)
    with samples.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=SAMPLE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(
            {
                "sample": SAMPLE,
                "fastq_1": "R1.fq",
                "fastq_2": "R2.fq",
                "layout": "PE",
                "assay": assay,
                "target": target[0],
                "peak_mode": target[1],
                "genome": "hs",
                "bowtie2_index": "idx",
                "experiment": "E1",
                "biological_replicate": "1",
            }
        )
    config = {
        "samples": "config/samples.tsv",
        "outdir": outdir,
        "threads": 1,
        "use_control": False,
        "multiqc": False,
        "replicate_analysis": False,
        "analysis_profile": profile,
        "config_schema_version": "2.0.0",
        "blacklist_mode": "off",
        "genome_resources": {"hs": {"effective_genome_size": "hs"}},
        "qc": {
            "frip": False,
            "signal_tracks": False,
            "summary": False,
            "cuttag_fragment_size": cuttag_fragment_size,
            "cross_correlation": cross_correlation,
        },
    }
    if profile == "updated":
        config["duplicate_backend"] = "picard"
        config["control_mode"] = "matched"
    config_path = tmp_path / "config" / "config.yaml"
    config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
    return config_path


def _seed_identity(
    tmp_path: Path,
    snakemake_executable: str,
    config_path: Path,
    outdir: str = "results",
):
    seed_run_identity(
        tmp_path, snakemake_executable, SNAKEFILE, config_path, outdir=outdir
    )


def _run_snakemake(
    tmp_path: Path,
    snakemake_executable: str,
    config_path: Path,
    rule: str,
    target: str,
):
    return _run(
        [
            snakemake_executable,
            "-s",
            SNAKEFILE,
            "--workflow-profile",
            "none",
            "--directory",
            tmp_path,
            "--configfile",
            config_path,
            "--cores",
            "1",
            "--printshellcmds",
            "--allowed-rules",
            rule,
            "--",
            target,
        ],
        cwd=tmp_path,
    )


def _write_stage_metrics(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "total_alignments": 6,
                "kept_alignments": 2,
                "dropped_alignments": 4,
                "dropped_templates": 2,
                "dropped_alignments_by_reason": {},
                "dropped_templates_by_reason": {},
            }
        ),
        encoding="utf-8",
    )


def _write_mito(path: Path) -> None:
    _write_tsv(
        path,
        MITO_QC_V2_COLUMNS,
        {
            "schema_version": "2",
            "sample": SAMPLE,
            "organelle_reads": "2",
            "total_mapped_reads": "8",
            "fraction": "0.250000",
            "mode": "exclude",
            "status": "ok",
            "input_stage": "pre_duplicate_mapq_bam",
            "matched_contigs": "chrM",
            "organelle_chroms_explicit": "false",
        },
    )


def _write_fragment_qc(path: Path) -> None:
    row = {column: "NA" for column in FRAGMENT_QC_V2_COLUMNS}
    row.update(
        {
            "schema_version": "2",
            "sample": SAMPLE,
            "layout": "PE",
            "fragment_count": "1",
            "fragment_mean": "150",
            "fragment_median": "150",
            "fragment_mode": "150",
            "fragment_min": "150",
            "fragment_max": "150",
            "fraction_lt_150": "0.000000",
            "fraction_150_300": "1.000000",
            "fraction_300_500": "0.000000",
            "fraction_ge_500": "0.000000",
            "fraction_lt_120": "0.000000",
            "status": "ok",
            "input_stage": "analysis_bam",
            "duplicate_analysis": "exclude",
            "blacklist_applied": "false",
            "templates_total": "4",
            "templates_valid": "3",
            "templates_invalid": "1",
            "invalid_tlen_zero": "0",
            "invalid_orphan": "1",
            "invalid_cross_contig": "0",
            "invalid_inconsistent": "0",
            "alignments_total": "8",
            "upstream_dropped_templates": "2",
            "upstream_dropped_alignments": "4",
            "upstream_alignments_total": "6",
            "upstream_alignments_kept": "2",
        }
    )
    _write_tsv(path, FRAGMENT_QC_V2_COLUMNS, row)


def _write_mnase_qc(path: Path) -> None:
    row = {column: "NA" for column in MNASE_QC_V2_COLUMNS}
    row.update(
        {
            "schema_version": "2",
            "sample": SAMPLE,
            "assay": "mnase",
            "peak_mode": "nucleosome",
            "sub_min": "80",
            "sub_max": "120",
            "mono_min": "140",
            "mono_max": "210",
            "di_min": "300",
            "di_max": "410",
            "dyad_min": "135",
            "dyad_max": "195",
            "analysis_status": "ok",
            "analysis_total_templates": "6",
            "analysis_valid_templates": "5",
            "analysis_invalid_templates": "1",
            "analysis_tlen_zero_templates": "0",
            "analysis_orphan_templates": "1",
            "analysis_cross_contig_templates": "0",
            "analysis_inconsistent_templates": "0",
            "analysis_classified_templates": "5",
            "analysis_unclassified_templates": "0",
            "analysis_fraction_classified": "1.000000",
            "window_overlap": "false",
            "sub_fragments": "1",
            "mono_fragments": "2",
            "di_fragments": "2",
            "upstream_dropped_templates": "0",
            "upstream_dropped_alignments": "0",
        }
    )
    _write_tsv(path, MNASE_QC_V2_COLUMNS, row)


def _write_policy(path: Path, assay: str, profile: str = "legacy") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    stage_counts = {
        "analysis_input_alignments": 6,
        "analysis_kept_alignments": 2,
        "analysis_valid_templates": 1,
        "analysis_invalid_templates": 0,
        "upstream_dropped_templates": 2,
        "upstream_dropped_alignments": 4,
        "blacklist_dropped_templates": 0,
        "blacklist_dropped_alignments": 0,
        "organelle_reads": 2,
        "organelle_total_mapped_reads": 8,
    }
    assert set(stage_counts) == set(POLICY_STAGE_COUNT_KEYS)
    policy = {
        "schema_version": 1,
        "scientific_policy_revision": 5,
        "sample": SAMPLE,
        "assay": assay,
        "profile": profile,
        "resolved_duplicate_analysis": "exclude",
        "resolved_max_insert": 700,
        "organelle_mode": "exclude",
        "organelle_status": "ok",
        "matched_organelle_contigs": ["chrM"],
        "blacklist_mode": "off",
        "blacklist_applied": False,
        "input_stages": {
            key: "not_applicable" for key in ANALYSIS_POLICY_INPUT_STAGE_KEYS
        },
        "stage_counts": stage_counts,
        "reference_summary": {"genome": "hs", "organelle_chroms_explicit": False},
        "qc_schema_versions": {
            "sample": "2",
            "mnase": "2",
            "fragment": "2",
            "mito": "2",
        },
    }
    path.write_text(
        json.dumps(policy, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _write_summary_inputs(
    tmp_path: Path, assay: str, *, outdir: str = "results", profile: str = "legacy"
) -> None:
    base = tmp_path / outdir / SAMPLE
    _write_tsv(
        base / "01_qc" / f"{SAMPLE}.peak_counts.tsv",
        PEAK_COLUMNS,
        {
            "sample": SAMPLE,
            "peak_mode": "narrow",
            "peaks": "3",
            "blacklist_filtered_peaks": "1",
        },
    )
    _write_tsv(
        base / "01_qc" / f"{SAMPLE}.frip.tsv",
        FRIP_COLUMNS,
        {
            "sample": SAMPLE,
            "total_reads": "10",
            "reads_in_peaks": "4",
            "frip": "0.4",
            "bam": f"{outdir}/{SAMPLE}/02_align/{SAMPLE}.analysis.bam",
            "peaks": f"{outdir}/{SAMPLE}/04_peaks/{SAMPLE}_peaks.narrowPeak",
        },
    )
    _write_tsv(
        base / "01_qc" / f"{SAMPLE}.library_complexity.tsv",
        LIBRARY_COLUMNS,
        {"sample": SAMPLE, "metrics_source": "picard"},
    )
    _write_tsv(
        base / "01_qc" / f"{SAMPLE}.nrf_pbc.tsv",
        NRF_COLUMNS,
        {
            "sample": SAMPLE,
            "total_fragments": "4",
            "distinct_fragments": "3",
            "one_read_fragments": "0",
            "two_read_fragments": "3",
            "nrf": "0.75",
            "pbc1": "0.75",
            "pbc2": "3.0",
        },
    )
    for name in ("final", "analysis"):
        bam = base / "02_align" / f"{SAMPLE}.{name}.bam"
        bam.parent.mkdir(parents=True, exist_ok=True)
        bam.write_bytes(b"fixture placeholder; qc_summary consumes it as a path\n")
    _write_mito(base / "05_qc" / f"{SAMPLE}.mito_fraction.tsv")
    _write_policy(base / "05_qc" / f"{SAMPLE}.analysis_policy.json", assay, profile)


def _write_cross_correlation(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "S1\t1000\t150\t0.9\t150\t0.8\t150\t0.1\t1.2\t1.1\t0\n",
        encoding="utf-8",
    )


@pytest.mark.parametrize(
    (
        "assay",
        "cuttag_fragment_size",
        "expected_valid",
        "expected_invalid",
        "fragment_stage",
        "mnase_stage",
    ),
    [
        ("chipseq", True, None, None, "not_applicable", "not_applicable"),
        ("cuttag", False, None, None, "not_applicable", "not_applicable"),
        ("cuttag", True, 3, 1, "analysis_bam", "not_applicable"),
        ("mnase", True, 5, 1, "not_applicable", "analysis_bam"),
    ],
)
def test_analysis_policy_rule_handles_optional_qc(
    tmp_path,
    snakemake_executable,
    assay,
    cuttag_fragment_size,
    expected_valid,
    expected_invalid,
    fragment_stage,
    mnase_stage,
):
    config_path = _write_project(
        tmp_path,
        assay=assay,
        cuttag_fragment_size=cuttag_fragment_size,
        cross_correlation=False,
    )
    _seed_identity(tmp_path, snakemake_executable, config_path)

    base = tmp_path / "results" / SAMPLE
    _write_stage_metrics(
        base / "05_qc" / f"{SAMPLE}.blacklist_dropped.stage_metrics.json"
    )
    _write_mito(base / "05_qc" / f"{SAMPLE}.mito_fraction.tsv")
    if assay == "cuttag" and cuttag_fragment_size:
        _write_fragment_qc(base / "01_qc" / f"{SAMPLE}.cuttag_fragment_size.tsv")
    if assay == "mnase":
        _write_mnase_qc(base / "01_qc" / f"{SAMPLE}.mnase_qc_summary.tsv")

    target = f"results/{SAMPLE}/05_qc/{SAMPLE}.analysis_policy.json"
    _run_snakemake(
        tmp_path, snakemake_executable, config_path, "analysis_policy", target
    )

    policy = json.loads((tmp_path / target).read_text(encoding="utf-8"))
    assert policy["sample"] == SAMPLE
    assert policy["assay"] == assay
    assert policy["input_stages"]["fragment_size"] == fragment_stage
    assert policy["input_stages"]["mnase_classification"] == mnase_stage
    if expected_valid is not None:
        assert policy["stage_counts"]["analysis_valid_templates"] == expected_valid
        assert policy["stage_counts"]["analysis_invalid_templates"] == expected_invalid


@pytest.mark.parametrize(
    (
        "assay",
        "cuttag_fragment_size",
        "cross_correlation",
        "fragment_valid",
        "cc_status",
    ),
    [
        ("chipseq", True, False, "NA", "not_enabled"),
        ("cuttag", True, False, "3", "not_enabled"),
        ("cuttag", False, True, "NA", "reported"),
        ("cuttag", True, True, "3", "reported"),
    ],
)
def test_qc_summary_rule_handles_optional_qc(
    tmp_path,
    snakemake_executable,
    assay,
    cuttag_fragment_size,
    cross_correlation,
    fragment_valid,
    cc_status,
):
    config_path = _write_project(
        tmp_path,
        assay=assay,
        cuttag_fragment_size=cuttag_fragment_size,
        cross_correlation=cross_correlation,
    )
    _seed_identity(tmp_path, snakemake_executable, config_path)
    _write_summary_inputs(tmp_path, assay)

    base = tmp_path / "results" / SAMPLE
    if assay == "cuttag" and cuttag_fragment_size:
        _write_fragment_qc(base / "01_qc" / f"{SAMPLE}.cuttag_fragment_size.tsv")
    if cross_correlation:
        _write_cross_correlation(
            base / "05_qc" / "cross_correlation" / f"{SAMPLE}.cc.qc"
        )

    target = f"results/{SAMPLE}/01_qc/{SAMPLE}.qc_summary.tsv"
    _run_snakemake(tmp_path, snakemake_executable, config_path, "qc_summary", target)

    with (tmp_path / target).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    assert reader.fieldnames == list(SAMPLE_QC_V2_COLUMNS)
    assert len(rows) == 1
    row = rows[0]
    assert row["schema_version"] == "2"
    assert row["sample"] == SAMPLE
    assert row["assay"] == assay
    assert row["raw_peak_count"] == "3"
    assert row["delivered_peak_count"] == "3"
    assert row["analysis_valid_fragments"] == fragment_valid
    assert row["cross_correlation_status"] == cc_status
    if cc_status == "reported":
        assert row["cross_correlation_nsc"] == "1.2"
        assert row["cross_correlation_rsc"] == "1.1"
    else:
        assert row["cross_correlation_nsc"] == "NA"
        assert row["cross_correlation_rsc"] == "NA"


POLICY_SENTINEL = "R1_POLICY_UNEXPECTED_SIDE_EFFECT"
SUMMARY_SENTINEL = "R1_SUMMARY_UNEXPECTED_SIDE_EFFECT"


def _outdir_for(path_kind: str, sentinel: str) -> str:
    if path_kind == "space":
        return "results with spaces"
    return f"results $(touch {sentinel}) 'single' \"double\""


@pytest.mark.parametrize("profile", ["legacy", "updated"])
@pytest.mark.parametrize("assay", ["cuttag", "mnase"])
@pytest.mark.parametrize("path_kind", ["space", "literal"])
def test_analysis_policy_rule_preserves_literal_paths(
    tmp_path, snakemake_executable, path_kind, assay, profile
):
    outdir = _outdir_for(path_kind, POLICY_SENTINEL)
    config_path = _write_project(
        tmp_path,
        assay=assay,
        cuttag_fragment_size=True,
        cross_correlation=False,
        outdir=outdir,
        profile=profile,
    )
    _seed_identity(tmp_path, snakemake_executable, config_path, outdir=outdir)

    base = tmp_path / outdir / SAMPLE
    _write_stage_metrics(
        base / "05_qc" / f"{SAMPLE}.blacklist_dropped.stage_metrics.json"
    )
    _write_mito(base / "05_qc" / f"{SAMPLE}.mito_fraction.tsv")
    if assay == "cuttag":
        _write_fragment_qc(base / "01_qc" / f"{SAMPLE}.cuttag_fragment_size.tsv")
    else:
        _write_mnase_qc(base / "01_qc" / f"{SAMPLE}.mnase_qc_summary.tsv")

    target = f"{outdir}/{SAMPLE}/05_qc/{SAMPLE}.analysis_policy.json"
    _run_snakemake(
        tmp_path, snakemake_executable, config_path, "analysis_policy", target
    )

    assert (tmp_path / outdir).is_dir()
    assert (tmp_path / target).is_file()
    assert not (tmp_path / POLICY_SENTINEL).exists(), (
        "analysis_policy expanded the optional QC path as shell code"
    )
    policy = json.loads((tmp_path / target).read_text(encoding="utf-8"))
    assert policy["profile"] == profile
    assert policy["input_stages"]["fragment_size"] == (
        "analysis_bam" if assay == "cuttag" else "not_applicable"
    )
    assert policy["input_stages"]["mnase_classification"] == (
        "analysis_bam" if assay == "mnase" else "not_applicable"
    )
    assert policy["stage_counts"]["analysis_valid_templates"] == (
        3 if assay == "cuttag" else 5
    )


@pytest.mark.parametrize("profile", ["legacy", "updated"])
@pytest.mark.parametrize("path_kind", ["space", "literal"])
@pytest.mark.parametrize(
    ("fragment_enabled", "cross_correlation_enabled"),
    [(True, False), (False, True), (True, True)],
)
def test_qc_summary_rule_preserves_literal_paths(
    tmp_path,
    snakemake_executable,
    path_kind,
    profile,
    fragment_enabled,
    cross_correlation_enabled,
):
    outdir = _outdir_for(path_kind, SUMMARY_SENTINEL)
    config_path = _write_project(
        tmp_path,
        assay="cuttag",
        cuttag_fragment_size=fragment_enabled,
        cross_correlation=cross_correlation_enabled,
        outdir=outdir,
        profile=profile,
    )
    _seed_identity(tmp_path, snakemake_executable, config_path, outdir=outdir)
    _write_summary_inputs(tmp_path, "cuttag", outdir=outdir, profile=profile)

    base = tmp_path / outdir / SAMPLE
    if fragment_enabled:
        _write_fragment_qc(base / "01_qc" / f"{SAMPLE}.cuttag_fragment_size.tsv")
    if cross_correlation_enabled:
        _write_cross_correlation(
            base / "05_qc" / "cross_correlation" / f"{SAMPLE}.cc.qc"
        )

    target = f"{outdir}/{SAMPLE}/01_qc/{SAMPLE}.qc_summary.tsv"
    _run_snakemake(tmp_path, snakemake_executable, config_path, "qc_summary", target)

    assert (tmp_path / outdir).is_dir()
    assert (tmp_path / target).is_file()
    assert not (tmp_path / SUMMARY_SENTINEL).exists(), (
        "qc_summary expanded the optional QC path as shell code"
    )
    with (tmp_path / target).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    assert reader.fieldnames == list(SAMPLE_QC_V2_COLUMNS)
    assert len(rows) == 1
    row = rows[0]
    assert row["assay"] == "cuttag"
    assert row["final_bam"] == f"{outdir}/{SAMPLE}/02_align/{SAMPLE}.final.bam"
    assert row["analysis_bam"] == f"{outdir}/{SAMPLE}/02_align/{SAMPLE}.analysis.bam"
    assert row["analysis_valid_fragments"] == ("3" if fragment_enabled else "NA")
    if cross_correlation_enabled:
        assert row["cross_correlation_status"] == "reported"
        assert row["cross_correlation_nsc"] == "1.2"
        assert row["cross_correlation_rsc"] == "1.1"
    else:
        assert row["cross_correlation_status"] == "not_enabled"
        assert row["cross_correlation_nsc"] == "NA"
        assert row["cross_correlation_rsc"] == "NA"
