"""Contracts for the general PE analysis-template component.

The script tests use a controlled fake ``samtools view`` so template
classification is exercised without an external tool; the real production
rules are covered by ``test/real_execution/test_analysis_template_counts.py``.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import stat
import subprocess
import sys

from encode_pipeline.qc.contract import (
    ANALYSIS_TEMPLATE_METRICS_COLUMNS,
    ANALYSIS_TEMPLATE_METRICS_V1_COLUMNS,
    FRAGMENT_QC_V2_COLUMNS,
    MITO_QC_V2_COLUMNS,
)
from encode_pipeline.qc.policy import parse_analysis_policy_bytes
from scripts import write_analysis_policy

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "calc_analysis_template_metrics.py"

SAM_HEADER = [
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:100000",
    "@SQ\tSN:chr2\tLN:100000",
]

VALID_PAIR = [
    "v1\t99\tchr1\t100\t60\t50M\t=\t250\t150\t*\t*",
    "v1\t147\tchr1\t250\t60\t50M\t=\t100\t-150\t*\t*",
]
NEGATIVE_TLEN_FIRST = [
    "n1\t147\tchr1\t300\t60\t50M\t=\t100\t-200\t*\t*",
    "n1\t99\tchr1\t100\t60\t50M\t=\t300\t200\t*\t*",
]
TLEN_ZERO = [
    "z1\t99\tchr1\t1500\t60\t50M\t=\t1600\t0\t*\t*",
    "z1\t147\tchr1\t1600\t60\t50M\t=\t1500\t0\t*\t*",
]
ORPHAN = ["o1\t99\tchr1\t1800\t60\t50M\t=\t1950\t150\t*\t*"]
CROSS_CONTIG = [
    "c1\t99\tchr1\t2100\t60\t50M\tchr2\t2300\t200\t*\t*",
    "c1\t147\tchr2\t2300\t60\t50M\tchr1\t2100\t-200\t*\t*",
]
INCONSISTENT = [
    "i1\t99\tchr1\t2500\t60\t50M\t=\t2700\t200\t*\t*",
    "i1\t147\tchr1\t2700\t60\t50M\t=\t2500\t-250\t*\t*",
]


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


def _read_row(path: Path) -> dict:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert list(reader.fieldnames or []) == list(ANALYSIS_TEMPLATE_METRICS_COLUMNS)
        rows = list(reader)
    assert len(rows) == 1
    return rows[0]


def _fake_samtools(tmp_path: Path, *, exit_code: int = 0) -> Path:
    path = tmp_path / "samtools"
    path.write_text(
        "#!/bin/sh\n"
        'if [ "$1" = "view" ]; then\n'
        f'  cat "$SAM_SOURCE"\n'
        f"  exit {exit_code}\n"
        "fi\n"
        "exit 2\n",
        encoding="utf-8",
    )
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


def _run_script(
    tmp_path: Path,
    *,
    records,
    layout: str = "PE",
    bam: Path | None = None,
    stage_metrics: Path | None = None,
    samtools_exit: int = 0,
    sample: str = "S1",
) -> dict:
    sam = tmp_path / "input.sam"
    sam.write_text("\n".join(SAM_HEADER + list(records)) + "\n", encoding="utf-8")
    samtools = _fake_samtools(tmp_path, exit_code=samtools_exit)
    output = tmp_path / "metrics.tsv"
    argv = [
        sys.executable,
        str(SCRIPT),
        "--sample",
        sample,
        "--layout",
        layout,
        "--bam",
        str(bam if bam is not None else tmp_path / "analysis.bam"),
        "--output",
        str(output),
    ]
    if stage_metrics is not None:
        argv[argv.index("--output") : argv.index("--output")] = [
            "--stage-metrics",
            str(stage_metrics),
        ]
    env = os.environ.copy()
    env["SAMTOOLS"] = str(samtools)
    env["SAM_SOURCE"] = str(sam)
    result = subprocess.run(argv, capture_output=True, text=True, env=env)
    assert result.returncode == 0, result.stderr
    return _read_row(output)


def _stage_metrics(path: Path) -> Path:
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
    return path


def test_valid_negative_tlen_is_kept_and_all_totals_conserve(tmp_path):
    bam = tmp_path / "analysis.bam"
    bam.write_bytes(b"placeholder; the fake samtools ignores it\n")
    stage = _stage_metrics(tmp_path / "stage.json")

    row = _run_script(
        tmp_path,
        records=(
            VALID_PAIR
            + NEGATIVE_TLEN_FIRST
            + TLEN_ZERO
            + ORPHAN
            + CROSS_CONTIG
            + INCONSISTENT
        ),
        bam=bam,
        stage_metrics=stage,
    )

    assert row["schema_version"] == "2"
    assert row["analysis_status"] == "ok"
    assert row["input_stage"] == "analysis_bam"
    assert row["analysis_total_templates"] == "6"
    assert row["analysis_valid_templates"] == "2"
    assert row["analysis_invalid_templates"] == "4"
    assert row["analysis_tlen_zero_templates"] == "1"
    assert row["analysis_orphan_templates"] == "1"
    assert row["analysis_cross_contig_templates"] == "1"
    assert row["analysis_inconsistent_templates"] == "1"
    assert row["alignments_total"] == "11"
    assert row["upstream_dropped_templates"] == "2"
    assert row["upstream_dropped_alignments"] == "4"
    assert int(row["analysis_valid_templates"]) + int(
        row["analysis_invalid_templates"]
    ) == int(row["analysis_total_templates"])
    assert sum(
        int(row[column])
        for column in (
            "analysis_tlen_zero_templates",
            "analysis_orphan_templates",
            "analysis_cross_contig_templates",
            "analysis_inconsistent_templates",
        )
    ) == int(row["analysis_invalid_templates"])


def test_zero_valid_templates_is_a_legitimate_conserved_zero(tmp_path):
    bam = tmp_path / "analysis.bam"
    bam.write_bytes(b"placeholder\n")

    row = _run_script(
        tmp_path,
        records=TLEN_ZERO + ORPHAN + CROSS_CONTIG + INCONSISTENT,
        bam=bam,
    )

    assert row["analysis_status"] == "no_valid_templates"
    assert row["analysis_total_templates"] == "4"
    assert row["analysis_valid_templates"] == "0"
    assert row["analysis_invalid_templates"] == "4"
    assert row["upstream_dropped_templates"] == "NA"
    assert row["upstream_dropped_alignments"] == "NA"


def test_single_end_is_explicitly_not_applicable_with_na_counts(tmp_path):
    stage = _stage_metrics(tmp_path / "stage.json")

    row = _run_script(
        tmp_path,
        records=[],
        layout="SE",
        bam=tmp_path / "missing.bam",
        stage_metrics=stage,
    )

    assert row["layout"] == "SE"
    assert row["analysis_status"] == "not_applicable"
    assert row["input_stage"] == "analysis_bam"
    for column in (
        "analysis_total_templates",
        "analysis_valid_templates",
        "analysis_invalid_templates",
        "analysis_tlen_zero_templates",
        "analysis_orphan_templates",
        "analysis_cross_contig_templates",
        "analysis_inconsistent_templates",
        "alignments_total",
    ):
        assert row[column] == "NA"
    assert row["upstream_dropped_templates"] == "2"
    assert row["upstream_dropped_alignments"] == "4"


def test_missing_bam_is_parse_failed_and_never_zero(tmp_path):
    row = _run_script(
        tmp_path,
        records=VALID_PAIR,
        bam=tmp_path / "missing.bam",
    )

    assert row["analysis_status"] == "parse_failed"
    for column in (
        "analysis_total_templates",
        "analysis_valid_templates",
        "analysis_invalid_templates",
        "analysis_tlen_zero_templates",
        "analysis_orphan_templates",
        "analysis_cross_contig_templates",
        "analysis_inconsistent_templates",
        "alignments_total",
    ):
        assert row[column] == "NA"


def test_unreadable_bam_is_parse_failed_and_never_zero(tmp_path):
    bam = tmp_path / "analysis.bam"
    bam.write_bytes(b"not a bam\n")

    row = _run_script(
        tmp_path,
        records=VALID_PAIR,
        bam=bam,
        samtools_exit=1,
    )

    assert row["analysis_status"] == "parse_failed"
    assert row["analysis_valid_templates"] == "NA"
    assert row["analysis_invalid_templates"] == "NA"
    assert row["analysis_total_templates"] == "NA"


def _policy_inputs(tmp_path: Path):
    stage = _stage_metrics(tmp_path / "stage.json")
    mito = tmp_path / "mito.tsv"
    _write_tsv(
        mito,
        MITO_QC_V2_COLUMNS,
        {
            "schema_version": "2",
            "sample": "S1",
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
    return stage, mito


def _run_writer(
    tmp_path: Path,
    monkeypatch,
    component: Path | None,
    *,
    max_insert: str = "2000",
    assay: str = "chipseq",
    fragment_qc: Path | None = None,
) -> dict:
    stage, mito = _policy_inputs(tmp_path)
    output = tmp_path / "policy.json"
    argv = [
        "write_analysis_policy.py",
        "--sample",
        "S1",
        "--assay",
        assay,
        "--genome",
        "hs",
        "--profile",
        "legacy",
        "--policy-revision",
        "5",
        "--duplicate-analysis",
        "exclude",
        "--max-insert",
        max_insert,
        "--organelle-mode",
        "report",
        "--blacklist-mode",
        "off",
        "--blacklist-applied",
        "false",
        "--stage-metrics",
        str(stage),
        "--mito",
        str(mito),
        "--output",
        str(output),
    ]
    if component is not None:
        argv[argv.index("--output") : argv.index("--output")] = [
            "--analysis-template-metrics",
            str(component),
        ]
    if fragment_qc is not None:
        argv[argv.index("--output") : argv.index("--output")] = [
            "--fragment-qc",
            str(fragment_qc),
        ]
    monkeypatch.setattr(sys, "argv", argv)
    assert write_analysis_policy.main() == 0
    return json.loads(output.read_text(encoding="utf-8"))


def _component_row(status="ok", valid="4", invalid="1", total="5"):
    return {
        "schema_version": "2",
        "sample": "S1",
        "layout": "PE",
        "analysis_status": status,
        "input_stage": "analysis_bam",
        "analysis_total_templates": total,
        "analysis_valid_templates": valid,
        "analysis_invalid_templates": invalid,
        "analysis_tlen_zero_templates": "0",
        "analysis_orphan_templates": "0",
        "analysis_cross_contig_templates": "0",
        "analysis_inconsistent_templates": "0",
        "alignments_total": "10",
        "upstream_dropped_templates": "2",
        "upstream_dropped_alignments": "4",
    }


def test_component_reader_accepts_versionless_v1_and_v2(tmp_path, monkeypatch):
    v1 = tmp_path / "v1.tsv"
    _write_tsv(
        v1,
        ANALYSIS_TEMPLATE_METRICS_V1_COLUMNS,
        {
            key: value
            for key, value in _component_row().items()
            if key != "schema_version"
        },
    )
    policy = _run_writer(tmp_path, monkeypatch, v1)
    assert policy["analysis_template_status"] == "ok"
    assert policy["stage_counts"]["analysis_valid_templates"] == 4
    assert policy["stage_counts"]["analysis_invalid_templates"] == 1

    v2 = tmp_path / "v2.tsv"
    _write_tsv(v2, ANALYSIS_TEMPLATE_METRICS_COLUMNS, _component_row())
    policy = _run_writer(tmp_path, monkeypatch, v2)
    assert policy["schema_version"] == 2
    assert policy["analysis_template_status"] == "ok"
    assert policy["input_stages"]["analysis_template_metrics"] == "analysis_bam"
    assert policy["stage_counts"]["analysis_valid_templates"] == 4


def test_v2_policy_round_trips_through_the_strict_parser(tmp_path, monkeypatch):
    component = tmp_path / "v2.tsv"
    _write_tsv(component, ANALYSIS_TEMPLATE_METRICS_COLUMNS, _component_row())

    policy = _run_writer(tmp_path, monkeypatch, component)
    parsed = parse_analysis_policy_bytes(json.dumps(policy).encode("utf-8"))

    assert parsed["schema_version"] == 2
    assert parsed["analysis_template_status"] == "ok"
    assert parsed["input_stages"]["analysis_template_metrics"] == "analysis_bam"
    assert parsed["stage_counts"]["analysis_valid_templates"] == 4


def test_unknown_component_version_is_parse_failed_not_counts(tmp_path, monkeypatch):
    component = tmp_path / "unknown.tsv"
    row = _component_row()
    row["schema_version"] = "999"
    _write_tsv(component, ANALYSIS_TEMPLATE_METRICS_COLUMNS, row)

    policy = _run_writer(tmp_path, monkeypatch, component)

    assert policy["analysis_template_status"] == "parse_failed"
    assert policy["input_stages"]["analysis_template_metrics"] == "analysis_bam"
    assert policy["stage_counts"]["analysis_valid_templates"] is None
    assert policy["stage_counts"]["analysis_invalid_templates"] is None


def test_single_end_component_is_policy_not_applicable(tmp_path, monkeypatch):
    component = tmp_path / "se.tsv"
    row = {column: "NA" for column in ANALYSIS_TEMPLATE_METRICS_COLUMNS}
    row.update(
        {
            "schema_version": "2",
            "sample": "S1",
            "layout": "SE",
            "analysis_status": "not_applicable",
            "input_stage": "analysis_bam",
        }
    )
    _write_tsv(component, ANALYSIS_TEMPLATE_METRICS_COLUMNS, row)

    policy = _run_writer(tmp_path, monkeypatch, component, max_insert="not_applicable")

    assert policy["analysis_template_status"] == "not_applicable"
    assert policy["input_stages"]["analysis_template_metrics"] == "not_applicable"
    assert policy["stage_counts"]["analysis_valid_templates"] is None
    assert policy["stage_counts"]["analysis_invalid_templates"] is None


def test_parse_failed_component_is_not_silently_zero(tmp_path, monkeypatch):
    component = tmp_path / "broken.tsv"
    row = {column: "NA" for column in ANALYSIS_TEMPLATE_METRICS_COLUMNS}
    row.update(
        {
            "schema_version": "2",
            "sample": "S1",
            "layout": "PE",
            "analysis_status": "parse_failed",
            "input_stage": "analysis_bam",
        }
    )
    _write_tsv(component, ANALYSIS_TEMPLATE_METRICS_COLUMNS, row)

    policy = _run_writer(tmp_path, monkeypatch, component)

    assert policy["analysis_template_status"] == "parse_failed"
    assert policy["stage_counts"]["analysis_valid_templates"] is None
    assert policy["stage_counts"]["analysis_invalid_templates"] is None


def test_absent_general_component_falls_back_to_legacy_fragment_source(
    tmp_path, monkeypatch
):
    fragment = tmp_path / "fragment.tsv"
    row = {column: "NA" for column in FRAGMENT_QC_V2_COLUMNS}
    row.update(
        {
            "schema_version": "2",
            "sample": "S1",
            "layout": "PE",
            "fragment_count": "3",
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
    _write_tsv(fragment, FRAGMENT_QC_V2_COLUMNS, row)

    policy = _run_writer(tmp_path, monkeypatch, None, fragment_qc=fragment)

    assert policy["analysis_template_status"] == "ok"
    assert policy["input_stages"]["analysis_template_metrics"] == "analysis_bam"
    assert policy["stage_counts"]["analysis_valid_templates"] == 3
    assert policy["stage_counts"]["analysis_invalid_templates"] == 1


def test_no_general_component_and_no_legacy_source_is_not_collected(
    tmp_path, monkeypatch
):
    policy = _run_writer(tmp_path, monkeypatch, None)

    assert policy["analysis_template_status"] == "not_collected"
    assert policy["input_stages"]["analysis_template_metrics"] == "not_applicable"
    assert policy["stage_counts"]["analysis_valid_templates"] is None
    assert policy["stage_counts"]["analysis_invalid_templates"] is None


def test_no_general_component_single_end_is_not_applicable(tmp_path, monkeypatch):
    policy = _run_writer(tmp_path, monkeypatch, None, max_insert="not_applicable")

    assert policy["analysis_template_status"] == "not_applicable"
    assert policy["input_stages"]["analysis_template_metrics"] == "not_applicable"
    assert policy["stage_counts"]["analysis_valid_templates"] is None
