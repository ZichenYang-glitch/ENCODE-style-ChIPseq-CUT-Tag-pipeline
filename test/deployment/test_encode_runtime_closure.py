"""Executable contracts for the exported ENCODE scientific runtime bundle.

The real ``WorkflowBuildIdentityProvider.source_manifest()`` is the only input
to the exported tree.  Scripts run with ``python -S`` and without
``PYTHONPATH`` so an installed HelixWeave package cannot mask a missing runtime
module.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from encode_pipeline.qc.contract import (
    MITO_QC_V2_COLUMNS,
    SAMPLE_QC_V2_COLUMNS,
)
from encode_pipeline.qc.policy import parse_analysis_policy_bytes
from encode_pipeline.services.defaults import create_default_workflow_registry
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider


REPO_ROOT = Path(__file__).resolve().parents[2]
SHARED_QC_CLOSURE = {
    "src/encode_pipeline/__init__.py",
    "src/encode_pipeline/qc/__init__.py",
    "src/encode_pipeline/qc/contract.py",
    "src/encode_pipeline/qc/histone.py",
    "src/encode_pipeline/qc/policy.py",
}
HELP_SCRIPTS = (
    "calc_cuttag_fragment_size.py",
    "write_analysis_policy.py",
    "assemble_qc_summary.py",
    "mnase_qc_summary.py",
    "aggregate_qc_summary.py",
    "make_cross_correlation_tagalign.py",
    "cross_correlation_spp_window.py",
    "write_cross_correlation_qc.py",
    "summarize_cross_correlation.py",
)
PEAK_HEADER = ("sample", "peak_mode", "peaks", "blacklist_filtered_peaks")
FRIP_HEADER = ("sample", "total_reads", "reads_in_peaks", "frip", "bam", "peaks")
LIBRARY_HEADER = (
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
NRF_HEADER = (
    "sample",
    "total_fragments",
    "distinct_fragments",
    "one_read_fragments",
    "two_read_fragments",
    "nrf",
    "pbc1",
    "pbc2",
)


def _source_manifest() -> tuple[tuple[str, bytes], ...]:
    return WorkflowBuildIdentityProvider(
        create_default_workflow_registry(environ={}),
        project_root=REPO_ROOT,
    ).source_manifest()


def _export(bundle: Path, manifest, *, include_shared_closure: bool) -> Path:
    for logical_path, content in manifest:
        if not include_shared_closure and logical_path.startswith("src/"):
            continue
        destination = bundle.joinpath(*logical_path.split("/"))
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
    return bundle


def _clean_environment() -> dict[str, str]:
    environment = {
        key: value
        for key, value in os.environ.items()
        if key not in {"PYTHONPATH", "PYTHONHOME"}
    }
    environment["PYTHONNOUSERSITE"] = "1"
    return environment


def _run_script(bundle: Path, name: str, *arguments: str):
    return subprocess.run(
        [sys.executable, "-S", str(bundle / "scripts" / name), *arguments],
        capture_output=True,
        text=True,
        env=_clean_environment(),
        cwd=bundle,
        check=False,
    )


def _write_tsv(path: Path, header, rows) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return path


def _write_tsv_row(path: Path, header, values: dict[str, str]) -> Path:
    return _write_tsv(path, header, [[values.get(column, "NA") for column in header]])


@pytest.fixture(scope="module")
def manifest():
    return _source_manifest()


@pytest.fixture(scope="module")
def bundle(tmp_path_factory, manifest) -> Path:
    root = tmp_path_factory.mktemp("encode-runtime") / "bundle"
    return _export(root, manifest, include_shared_closure=True)


def test_real_manifest_carries_only_the_shared_qc_closure(manifest):
    logical_paths = {logical_path for logical_path, _content in manifest}
    shared_source_paths = {
        logical_path
        for logical_path in logical_paths
        if logical_path.startswith("src/")
    }

    assert "scripts/qc_common.py" in logical_paths
    assert SHARED_QC_CLOSURE <= shared_source_paths
    # No platform module outside the shared QC package enters the bundle.
    assert all(
        logical_path == "src/encode_pipeline/__init__.py"
        or logical_path.startswith("src/encode_pipeline/qc/")
        for logical_path in shared_source_paths
    )


def test_exported_bundle_has_the_layout_qc_common_resolves(bundle):
    assert (bundle / "scripts" / "qc_common.py").is_file()
    assert (bundle / "src" / "encode_pipeline" / "qc" / "contract.py").is_file()
    assert (bundle / "src" / "encode_pipeline" / "__init__.py").is_file()


@pytest.mark.parametrize("script", HELP_SCRIPTS)
def test_exported_scripts_import_without_installed_helixweave(bundle, script):
    completed = _run_script(bundle, script, "--help")

    assert completed.returncode == 0, completed.stderr


def test_bundle_without_the_shared_closure_fails_closed(tmp_path, manifest):
    exported = _export(
        tmp_path / "bundle",
        manifest,
        include_shared_closure=False,
    )

    completed = _run_script(exported, "calc_cuttag_fragment_size.py", "--help")

    assert completed.returncode != 0
    assert "No module named 'encode_pipeline'" in completed.stderr


def test_exported_bundle_writes_policy_and_qc_summary(tmp_path, bundle):
    inputs = tmp_path / "inputs"
    inputs.mkdir()
    stage_metrics = inputs / "stage_metrics.json"
    stage_metrics.write_text(
        json.dumps(
            {
                "total_alignments": 12,
                "kept_alignments": 10,
                "dropped_alignments": 2,
                "dropped_templates": 1,
                "dropped_alignments_by_reason": {"blacklist": 2},
                "dropped_templates_by_reason": {"blacklist": 1},
            }
        )
        + "\n",
        encoding="utf-8",
    )
    mito = _write_tsv_row(
        inputs / "mito.tsv",
        MITO_QC_V2_COLUMNS,
        {
            "schema_version": "2",
            "sample": "S1",
            "organelle_reads": "3",
            "total_mapped_reads": "100",
            "fraction": "0.030000",
            "mode": "report",
            "status": "ok",
            "input_stage": "pre_duplicate_mapq_bam",
            "matched_contigs": "chrM",
            "organelle_chroms_explicit": "false",
        },
    )
    policy = inputs / "analysis_policy.json"

    completed = _run_script(
        bundle,
        "write_analysis_policy.py",
        "--sample",
        "S1",
        "--assay",
        "chipseq",
        "--genome",
        "hs",
        "--profile",
        "legacy",
        "--policy-revision",
        "5",
        "--duplicate-analysis",
        "exclude",
        "--max-insert",
        "2000",
        "--organelle-mode",
        "report",
        "--blacklist-mode",
        "analysis",
        "--blacklist-applied",
        "true",
        "--stage-metrics",
        str(stage_metrics),
        "--mito",
        str(mito),
        "--output",
        str(policy),
    )

    assert completed.returncode == 0, completed.stderr
    record = parse_analysis_policy_bytes(policy.read_bytes())
    assert record["organelle_status"] == "ok"
    assert record["matched_organelle_contigs"] == ["chrM"]
    assert record["stage_counts"]["organelle_reads"] == 3
    assert record["stage_counts"]["organelle_total_mapped_reads"] == 100
    assert record["input_stages"]["fragment_size"] == "not_applicable"

    peak_counts = _write_tsv(
        inputs / "peak_counts.tsv",
        PEAK_HEADER,
        [["S1", "narrow", "15000", "12000"]],
    )
    frip = _write_tsv(
        inputs / "frip.tsv",
        FRIP_HEADER,
        [["S1", "10000000", "2500000", "0.250000", "bam", "peaks"]],
    )
    library = _write_tsv(
        inputs / "library_complexity.tsv",
        LIBRARY_HEADER,
        [
            [
                "S1",
                "picard",
                "100",
                "50",
                "0",
                "10",
                "5",
                "20",
                "1",
                "0.100000",
                "1000",
                "200",
                "30",
            ]
        ],
    )
    nrf_pbc = _write_tsv(
        inputs / "nrf_pbc.tsv",
        NRF_HEADER,
        [["S1", "100", "90", "5", "45", "0.900000", "0.950000", "0.900000"]],
    )
    summary = inputs / "qc_summary.tsv"

    completed = _run_script(
        bundle,
        "assemble_qc_summary.py",
        "--sample",
        "S1",
        "--assay",
        "chipseq",
        "--target",
        "CTCF",
        "--genome",
        "hs",
        "--layout",
        "PE",
        "--peak-mode",
        "narrow",
        "--use-control",
        "False",
        "--control-type",
        "none",
        "--final-bam",
        "results/S1/02_align/S1.final.bam",
        "--analysis-bam",
        "results/S1/02_align/S1.analysis.bam",
        "--peaks-file",
        "results/S1/04_peaks/S1/S1_peaks.narrowPeak",
        "--has-blacklist",
        "yes",
        "--blacklist",
        "/opt/genomes/blacklist.bed",
        "--bl-bam",
        "results/S1/02_align/S1.blacklist_filtered.bam",
        "--bl-peaks",
        "results/S1/04_peaks/S1.blacklist_filtered.narrowPeak",
        "--peak-counts",
        str(peak_counts),
        "--frip",
        str(frip),
        "--library-complexity",
        str(library),
        "--nrf-pbc",
        str(nrf_pbc),
        "--mito-tsv",
        str(mito),
        "--analysis-policy",
        str(policy),
        "--output",
        str(summary),
    )

    assert completed.returncode == 0, completed.stderr
    with summary.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert tuple(reader.fieldnames or ()) == SAMPLE_QC_V2_COLUMNS
        row = next(reader)
    assert row["schema_version"] == "2"
    assert row["raw_peak_count"] == "15000"
    assert row["delivered_peak_count"] == "12000"
    assert row["organelle_status"] == "ok"
    assert row["organelle_reads"] == "3"
    assert row["blacklist_applied"] == "true"
    assert row["cross_correlation_status"] == "not_enabled"
    assert row["scientific_policy_revision"] == "5"


def test_exported_cross_correlation_detail_and_project_consumer(bundle):
    meta = bundle / "meta.json"
    meta.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "status": "ready",
                "assay": "chipseq",
                "layout": "SE",
                "input_stage": "analysis_se",
                "read_length": 50,
                "target_type": "tf",
                "variable_read_length": False,
                "window_lower": -500,
                "window_upper": 60,
            }
        )
    )
    raw = bundle / "S1.cc.qc"
    raw.write_text("S1\t1000\t200\t0.9\t50\t0.5\t1500\t0.1\t1.5\t2\t2\n")
    detail = bundle / "detail.tsv"
    result = _run_script(
        bundle,
        "write_cross_correlation_qc.py",
        "--sample",
        "S1",
        "--assay",
        "chipseq",
        "--layout",
        "SE",
        "--cc-qc",
        str(raw),
        "--metadata",
        str(meta),
        "--output",
        str(detail),
    )
    assert result.returncode == 0, result.stderr
    with detail.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert (
        row["interpretation_status"] == "ok" and row["threshold_hint"] == "local_pass"
    )
    summary = bundle / "cross_correlation_summary.tsv"
    result = _run_script(
        bundle,
        "summarize_cross_correlation.py",
        "--input",
        str(detail),
        "--output",
        str(summary),
    )
    assert result.returncode == 0, result.stderr
    with summary.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["interpretation_status"] == "ok" and row["rsc"] == "2"
