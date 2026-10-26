"""Declared schema-version validation for the versioned QC readers."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import pytest

from encode_pipeline.qc.contract import (
    MITO_QC_V1_COLUMNS,
    MITO_QC_V2_COLUMNS,
    SAMPLE_QC_V1_COLUMNS,
    SAMPLE_QC_V2_COLUMNS,
)
from encode_pipeline.qc.policy import parse_analysis_policy_bytes
from scripts import aggregate_qc_summary, assemble_qc_summary, write_analysis_policy
from scripts.qc_common import read_versioned_tsv


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


def _write_tsv(path: Path, header, rows) -> Path:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return path


def _row(columns, **values) -> list[str]:
    return [str(values.get(column, "NA")) for column in columns]


def _mito_v2_row(declared: str) -> dict[str, str]:
    values = {column: "NA" for column in MITO_QC_V2_COLUMNS}
    values.update(
        {
            "schema_version": declared,
            "sample": "S1",
            "organelle_reads": "3",
            "total_mapped_reads": "100",
            "fraction": "0.030000",
            "mode": "report",
            "status": "ok",
            "input_stage": "pre_duplicate_mapq_bam",
            "matched_contigs": "chrM",
            "organelle_chroms_explicit": "false",
        }
    )
    return values


def _write_stage_metrics(path: Path) -> Path:
    path.write_text(
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
    return path


def _write_policy(path: Path) -> Path:
    record = {
        "schema_version": 1,
        "scientific_policy_revision": 5,
        "sample": "S1",
        "assay": "chipseq",
        "profile": "legacy",
        "resolved_duplicate_analysis": "exclude",
        "resolved_max_insert": 2000,
        "organelle_mode": "report",
        "organelle_status": "ok",
        "matched_organelle_contigs": ["chrM"],
        "blacklist_mode": "analysis",
        "blacklist_applied": True,
        "input_stages": {
            "fragment_size": "not_applicable",
            "mnase_classification": "not_applicable",
            "duplicate_metrics": "pre_duplicate_markdup_bam",
            "nrf_pbc": "pre_duplicate_mapq_bam",
            "organelle": "pre_duplicate_mapq_bam",
            "blacklist": "final_bam_to_analysis_bam",
            "peaks": "analysis_bam",
        },
        "stage_counts": {
            "analysis_input_alignments": 12,
            "analysis_kept_alignments": 10,
            "analysis_valid_templates": None,
            "analysis_invalid_templates": None,
            "upstream_dropped_templates": 1,
            "upstream_dropped_alignments": 2,
            "blacklist_dropped_templates": 1,
            "blacklist_dropped_alignments": 2,
            "organelle_reads": 3,
            "organelle_total_mapped_reads": 100,
        },
        "reference_summary": {
            "genome": "hs",
            "organelle_chroms_explicit": False,
        },
        "qc_schema_versions": {
            "sample": "2",
            "mnase": "2",
            "fragment": "2",
            "mito": "2",
        },
    }
    path.write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def _invoke(monkeypatch, program: str, main, args: list[str]) -> None:
    monkeypatch.setattr(sys, "argv", [program, *args])
    main()


def _component_files(tmp_path: Path) -> dict[str, Path]:
    paths = {
        "peak_counts": tmp_path / "peak_counts.tsv",
        "frip": tmp_path / "frip.tsv",
        "library_complexity": tmp_path / "library_complexity.tsv",
        "nrf_pbc": tmp_path / "nrf_pbc.tsv",
        "policy": tmp_path / "analysis_policy.json",
        "fragment": tmp_path / "fragment.tsv",
    }
    _write_tsv(paths["peak_counts"], PEAK_HEADER, [["S1", "narrow", "15000", "12000"]])
    _write_tsv(
        paths["frip"],
        FRIP_HEADER,
        [["S1", "10000000", "2500000", "0.250000", "bam", "peaks"]],
    )
    _write_tsv(
        paths["library_complexity"],
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
    _write_tsv(
        paths["nrf_pbc"],
        NRF_HEADER,
        [["S1", "100", "90", "5", "45", "0.900000", "0.950000", "0.900000"]],
    )
    _write_policy(paths["policy"])
    return paths


def _assemble_arguments(inputs: dict[str, Path], mito: Path, output: Path):
    return [
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
        str(inputs["peak_counts"]),
        "--frip",
        str(inputs["frip"]),
        "--library-complexity",
        str(inputs["library_complexity"]),
        "--nrf-pbc",
        str(inputs["nrf_pbc"]),
        "--mito-tsv",
        str(mito),
        "--analysis-policy",
        str(inputs["policy"]),
        "--output",
        str(output),
    ]


def test_read_versioned_tsv_keeps_v1_tables_readable(tmp_path):
    path = _write_tsv(
        tmp_path / "mito-v1.tsv",
        MITO_QC_V1_COLUMNS,
        [_row(MITO_QC_V1_COLUMNS, sample="S1", status="ok")],
    )

    version, row = read_versioned_tsv(str(path), MITO_QC_V1_COLUMNS, MITO_QC_V2_COLUMNS)

    assert version == "1"
    assert row["sample"] == "S1"
    assert "schema_version" not in row


def test_read_versioned_tsv_accepts_a_matching_v2_declaration(tmp_path):
    path = _write_tsv(
        tmp_path / "mito-v2.tsv",
        MITO_QC_V2_COLUMNS,
        [[_mito_v2_row("2")[column] for column in MITO_QC_V2_COLUMNS]],
    )

    version, row = read_versioned_tsv(str(path), MITO_QC_V1_COLUMNS, MITO_QC_V2_COLUMNS)

    assert version == "2"
    assert row["schema_version"] == "2"
    assert row["status"] == "ok"


@pytest.mark.parametrize("declared", ("999", "", "1"))
def test_read_versioned_tsv_rejects_an_unknown_v2_declaration(tmp_path, declared):
    path = _write_tsv(
        tmp_path / "mito-v2.tsv",
        MITO_QC_V2_COLUMNS,
        [[_mito_v2_row(declared)[column] for column in MITO_QC_V2_COLUMNS]],
    )

    with pytest.raises(ValueError, match="schema version is unsupported"):
        read_versioned_tsv(str(path), MITO_QC_V1_COLUMNS, MITO_QC_V2_COLUMNS)


def test_aggregate_keeps_v1_tables_readable(tmp_path, monkeypatch):
    source = _write_tsv(
        tmp_path / "sample-v1.tsv",
        SAMPLE_QC_V1_COLUMNS,
        [_row(SAMPLE_QC_V1_COLUMNS, sample="S1", assay="chipseq")],
    )
    output = tmp_path / "aggregate.tsv"

    _invoke(
        monkeypatch,
        "aggregate_qc_summary.py",
        aggregate_qc_summary.main,
        ["--output", str(output), str(source)],
    )

    with output.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert tuple(reader.fieldnames or ()) == SAMPLE_QC_V1_COLUMNS
        rows = list(reader)
    assert [row["sample"] for row in rows] == ["S1"]


def test_aggregate_rejects_an_unknown_v2_declaration(tmp_path, monkeypatch, capsys):
    values = {column: "NA" for column in SAMPLE_QC_V2_COLUMNS}
    values.update({"schema_version": "999", "sample": "S1"})
    source = _write_tsv(
        tmp_path / "sample-v2.tsv",
        SAMPLE_QC_V2_COLUMNS,
        [[values[column] for column in SAMPLE_QC_V2_COLUMNS]],
    )
    output = tmp_path / "aggregate.tsv"

    with pytest.raises(SystemExit) as raised:
        _invoke(
            monkeypatch,
            "aggregate_qc_summary.py",
            aggregate_qc_summary.main,
            ["--output", str(output), str(source)],
        )

    assert raised.value.code == 1
    assert "unsupported schema_version" in capsys.readouterr().err


def test_policy_writer_marks_an_unknown_v2_component_as_parse_failed(
    tmp_path, monkeypatch
):
    stage_metrics = _write_stage_metrics(tmp_path / "stage_metrics.json")
    mito = _write_tsv(
        tmp_path / "mito-v2.tsv",
        MITO_QC_V2_COLUMNS,
        [[_mito_v2_row("999")[column] for column in MITO_QC_V2_COLUMNS]],
    )
    output = tmp_path / "analysis_policy.json"

    _invoke(
        monkeypatch,
        "write_analysis_policy.py",
        write_analysis_policy.main,
        [
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
            str(output),
        ],
    )

    record = parse_analysis_policy_bytes(output.read_bytes())
    assert record["organelle_status"] == "parse_failed"
    assert record["matched_organelle_contigs"] == []
    assert record["stage_counts"]["organelle_reads"] is None
    assert record["stage_counts"]["organelle_total_mapped_reads"] is None


def test_assembly_exits_nonzero_for_an_unknown_v2_component_version(
    tmp_path, monkeypatch, capsys
):
    inputs = _component_files(tmp_path)
    mito = _write_tsv(
        tmp_path / "mito-v2.tsv",
        MITO_QC_V2_COLUMNS,
        [[_mito_v2_row("999")[column] for column in MITO_QC_V2_COLUMNS]],
    )
    output = tmp_path / "qc_summary.tsv"

    with pytest.raises(SystemExit) as raised:
        _invoke(
            monkeypatch,
            "assemble_qc_summary.py",
            assemble_qc_summary.main,
            _assemble_arguments(inputs, mito, output),
        )

    assert raised.value.code == 1
    assert "schema version is unsupported" in capsys.readouterr().err
    assert not output.exists()


def test_assembly_accepts_a_matching_v2_component_version(tmp_path, monkeypatch):
    inputs = _component_files(tmp_path)
    mito = _write_tsv(
        tmp_path / "mito-v2.tsv",
        MITO_QC_V2_COLUMNS,
        [[_mito_v2_row("2")[column] for column in MITO_QC_V2_COLUMNS]],
    )
    output = tmp_path / "qc_summary.tsv"

    _invoke(
        monkeypatch,
        "assemble_qc_summary.py",
        assemble_qc_summary.main,
        _assemble_arguments(inputs, mito, output),
    )

    with output.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert tuple(reader.fieldnames or ()) == SAMPLE_QC_V2_COLUMNS
        row = next(reader)
    assert row["schema_version"] == "2"
    assert row["organelle_status"] == "ok"
    assert row["organelle_reads"] == "3"
