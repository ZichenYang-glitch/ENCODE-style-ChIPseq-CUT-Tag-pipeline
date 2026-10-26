"""Production-chain regression for caller vs delivered peak counts.

SCI-5 R3: ``peak_counts`` used to overwrite ``peaks`` with the blacklist
filtered count in analysis mode.  The versioned sample summary then reported
``raw_peak_count`` equal to the delivered count.  These tests execute the real
``peak_counts`` rule followed by the real ``qc_summary`` rule so the component
writer and the assembler are covered together.

Each case seeds a fresh workdir with the run identity marker, stages a caller
narrowPeak (3 lines) and a blacklist-filtered narrowPeak (1 line), and checks
the actual TSV rows for analysis, off, and copies blacklist modes under both
analysis profiles.
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
    MITO_QC_V2_COLUMNS,
    POLICY_STAGE_COUNT_KEYS,
    SAMPLE_QC_V2_COLUMNS,
)

pytestmark = pytest.mark.real_execution

REPO = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO / "workflow" / "Snakefile"
SAMPLE = "S1"
RAW_OUTDIR = "results"

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

CALLER_PEAKS = (
    "chr1\t100\t200\tcaller_1\t50\t.\t1.0\t2.0\t3.0\t10\n"
    "chr1\t300\t400\tcaller_2\t50\t.\t1.0\t2.0\t3.0\t10\n"
    "chr1\t500\t600\tcaller_3\t50\t.\t1.0\t2.0\t3.0\t10\n"
)
FILTERED_PEAKS = "chr1\t300\t400\tcaller_2\t50\t.\t1.0\t2.0\t3.0\t10\n"


@pytest.fixture(scope="module")
def snakemake_executable() -> str:
    candidate = resolve_tool("snakemake", "SNAKEMAKE")
    resolved = shutil.which(candidate)
    if not resolved and os.path.isabs(candidate) and os.access(candidate, os.X_OK):
        resolved = candidate
    if not resolved:
        require_external_tools(["snakemake"], "peak_counts/qc_summary chain")
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
    profile: str,
    blacklist_mode: str,
) -> Path:
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
                "assay": "chipseq",
                "target": "H3K27ac",
                "peak_mode": "narrow",
                "genome": "hs",
                "bowtie2_index": "idx",
                "experiment": "E1",
                "biological_replicate": "1",
            }
        )

    blacklist = tmp_path / "config" / "genome_resources" / "synthetic_blacklist.bed"
    blacklist.parent.mkdir(parents=True)
    blacklist.write_text("chr1\t300\t400\n", encoding="utf-8")

    config = {
        "samples": "config/samples.tsv",
        "outdir": RAW_OUTDIR,
        "threads": 1,
        "use_control": False,
        "multiqc": False,
        "replicate_analysis": False,
        "analysis_profile": profile,
        "config_schema_version": "2.0.0",
        "blacklist_mode": blacklist_mode,
        "genome_resources": {
            "hs": {
                "effective_genome_size": "hs",
                "blacklist": str(blacklist),
            }
        },
        "qc": {
            "frip": False,
            "signal_tracks": False,
            "summary": False,
            "cuttag_fragment_size": True,
            "cross_correlation": False,
        },
    }
    if profile == "updated":
        config["duplicate_backend"] = "picard"
        config["control_mode"] = "matched"
    config_path = tmp_path / "config" / "config.yaml"
    config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
    return config_path


def _seed_identity(
    tmp_path: Path, snakemake_executable: str, config_path: Path
) -> None:
    seed_run_identity(
        tmp_path, snakemake_executable, SNAKEFILE, config_path, outdir=RAW_OUTDIR
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


def _output_root(tmp_path: Path, blacklist_mode: str) -> Path:
    if blacklist_mode == "analysis":
        return tmp_path / RAW_OUTDIR / "analysis_filtered"
    return tmp_path / RAW_OUTDIR


def _stage_peak_inputs(tmp_path: Path, blacklist_mode: str) -> Path:
    root = _output_root(tmp_path, blacklist_mode)
    caller = root / SAMPLE / "04_peaks" / SAMPLE / f"{SAMPLE}_peaks.narrowPeak"
    caller.parent.mkdir(parents=True, exist_ok=True)
    caller.write_text(CALLER_PEAKS, encoding="utf-8")

    filtered = (
        root
        / SAMPLE
        / "04_peaks"
        / f"{SAMPLE}_blacklist_filtered"
        / f"{SAMPLE}_peaks.blacklist_filtered.narrowPeak"
    )
    filtered.parent.mkdir(parents=True, exist_ok=True)
    filtered.write_text(FILTERED_PEAKS, encoding="utf-8")
    return root


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
            "mode": "report",
            "status": "ok",
            "input_stage": "pre_duplicate_mapq_bam",
            "matched_contigs": "chrM",
            "organelle_chroms_explicit": "false",
        },
    )


def _write_policy(path: Path, *, profile: str, blacklist_mode: str) -> None:
    stage_counts = {
        "analysis_input_alignments": 6,
        "analysis_kept_alignments": 4,
        "analysis_valid_templates": 3,
        "analysis_invalid_templates": 0,
        "upstream_dropped_templates": 2,
        "upstream_dropped_alignments": 4,
        "blacklist_dropped_templates": 1,
        "blacklist_dropped_alignments": 2,
        "organelle_reads": 2,
        "organelle_total_mapped_reads": 8,
    }
    assert set(stage_counts) == set(POLICY_STAGE_COUNT_KEYS)
    input_stages = {key: "not_applicable" for key in ANALYSIS_POLICY_INPUT_STAGE_KEYS}
    if blacklist_mode == "analysis":
        input_stages["blacklist"] = "final_bam_to_analysis_bam"
    policy = {
        "schema_version": 1,
        "scientific_policy_revision": 5,
        "sample": SAMPLE,
        "assay": "chipseq",
        "profile": profile,
        "resolved_duplicate_analysis": "exclude",
        "resolved_max_insert": 700,
        "organelle_mode": "report",
        "organelle_status": "ok",
        "matched_organelle_contigs": ["chrM"],
        "blacklist_mode": blacklist_mode,
        "blacklist_applied": blacklist_mode == "analysis",
        "input_stages": input_stages,
        "stage_counts": stage_counts,
        "reference_summary": {"genome": "hs", "organelle_chroms_explicit": False},
        "qc_schema_versions": {
            "sample": "2",
            "mnase": "2",
            "fragment": "2",
            "mito": "2",
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(policy, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _stage_qc_summary_inputs(
    tmp_path: Path, *, profile: str, blacklist_mode: str
) -> Path:
    root = _output_root(tmp_path, blacklist_mode)
    sample_root = root / SAMPLE
    raw_sample_root = tmp_path / RAW_OUTDIR / SAMPLE

    _write_tsv(
        sample_root / "01_qc" / f"{SAMPLE}.frip.tsv",
        FRIP_COLUMNS,
        {
            "sample": SAMPLE,
            "total_reads": "10",
            "reads_in_peaks": "4",
            "frip": "0.400000",
            "bam": f"{RAW_OUTDIR}/{SAMPLE}/02_align/{SAMPLE}.analysis.bam",
            "peaks": f"{SAMPLE}_peaks.narrowPeak",
        },
    )
    _write_tsv(
        raw_sample_root / "01_qc" / f"{SAMPLE}.library_complexity.tsv",
        LIBRARY_COLUMNS,
        {"sample": SAMPLE, "metrics_source": "picard"},
    )
    _write_tsv(
        raw_sample_root / "01_qc" / f"{SAMPLE}.nrf_pbc.tsv",
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
        bam = raw_sample_root / "02_align" / f"{SAMPLE}.{name}.bam"
        bam.parent.mkdir(parents=True, exist_ok=True)
        bam.write_bytes(b"fixture placeholder; qc_summary consumes it as a path\n")
    # In analysis mode OUTDIR is the analysis_filtered branch, so the analysis
    # BAM has to exist there even though final.bam stays under the raw outdir.
    analysis_bam = sample_root / "02_align" / f"{SAMPLE}.analysis.bam"
    analysis_bam.parent.mkdir(parents=True, exist_ok=True)
    analysis_bam.write_bytes(b"fixture placeholder; qc_summary consumes it as a path\n")
    _write_mito(sample_root / "05_qc" / f"{SAMPLE}.mito_fraction.tsv")
    _write_policy(
        sample_root / "05_qc" / f"{SAMPLE}.analysis_policy.json",
        profile=profile,
        blacklist_mode=blacklist_mode,
    )
    return root


def _read_tsv(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    return reader.fieldnames, rows


@pytest.mark.parametrize("profile", ["legacy", "updated"])
@pytest.mark.parametrize("blacklist_mode", ["analysis", "off", "copies"])
def test_peak_counts_chain_preserves_raw_and_delivered_counts(
    tmp_path, snakemake_executable, profile, blacklist_mode
):
    config_path = _write_project(
        tmp_path, profile=profile, blacklist_mode=blacklist_mode
    )
    _seed_identity(tmp_path, snakemake_executable, config_path)

    root = _stage_peak_inputs(tmp_path, blacklist_mode)
    sample_root = root / SAMPLE
    (sample_root / "01_qc").mkdir(parents=True, exist_ok=True)

    peak_target = (
        f"{root.relative_to(tmp_path)}/{SAMPLE}/01_qc/{SAMPLE}.peak_counts.tsv"
    )
    _run_snakemake(
        tmp_path, snakemake_executable, config_path, "peak_counts", peak_target
    )

    peak_header, peak_rows = _read_tsv(tmp_path / peak_target)
    assert peak_header == list(PEAK_COLUMNS)
    assert len(peak_rows) == 1
    peak_row = peak_rows[0]
    assert peak_row["sample"] == SAMPLE
    assert peak_row["peak_mode"] == "narrow"
    # The caller count is always the true caller count, never the filtered one.
    assert peak_row["peaks"] == "3"
    assert peak_row["blacklist_filtered_peaks"] == (
        "1" if blacklist_mode in ("analysis", "copies") else "NA"
    )

    _stage_qc_summary_inputs(tmp_path, profile=profile, blacklist_mode=blacklist_mode)

    summary_target = (
        f"{root.relative_to(tmp_path)}/{SAMPLE}/01_qc/{SAMPLE}.qc_summary.tsv"
    )
    _run_snakemake(
        tmp_path, snakemake_executable, config_path, "qc_summary", summary_target
    )

    summary_header, summary_rows = _read_tsv(tmp_path / summary_target)
    assert summary_header == list(SAMPLE_QC_V2_COLUMNS)
    assert len(summary_rows) == 1
    summary = summary_rows[0]

    expected_delivered = "1" if blacklist_mode == "analysis" else "3"
    assert summary["schema_version"] == "2"
    assert summary["sample"] == SAMPLE
    assert summary["assay"] == "chipseq"
    assert summary["blacklist_mode"] == blacklist_mode
    assert summary["blacklist_applied"] == (
        "true" if blacklist_mode == "analysis" else "false"
    )
    assert summary["raw_peak_count"] == "3"
    assert summary["delivered_peak_count"] == expected_delivered
    assert summary["peak_count"] == expected_delivered
    assert summary["blacklist_filtered_peak_count"] == (
        "1" if blacklist_mode in ("analysis", "copies") else "NA"
    )
