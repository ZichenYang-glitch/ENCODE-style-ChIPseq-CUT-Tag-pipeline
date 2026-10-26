"""Producer contracts for versioned per-sample and aggregate QC summaries."""

from __future__ import annotations

import csv
import json
from pathlib import Path
import sys

import pytest

from encode_pipeline.adapters.encode_qc import (
    _QC_V2_HEADER as ADAPTER_QC_V2_HEADER,
)
from encode_pipeline.qc.contract import (
    MITO_QC_V2_COLUMNS,
    MNASE_QC_V2_COLUMNS,
    SAMPLE_QC_V1_COLUMNS,
    SAMPLE_QC_V2_COLUMNS,
)
from scripts import aggregate_qc_summary, assemble_qc_summary


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


def _write_tsv(path: Path, header, rows) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _policy(
    *,
    sample="S1",
    assay="chipseq",
    blacklist_mode="analysis",
    blacklist_applied=True,
    organelle_status="ok",
    duplicate="exclude",
    max_insert=2000,
) -> dict:
    return {
        "schema_version": 1,
        "scientific_policy_revision": 5,
        "sample": sample,
        "assay": assay,
        "profile": "legacy",
        "resolved_duplicate_analysis": duplicate,
        "resolved_max_insert": max_insert,
        "organelle_mode": "report",
        "organelle_status": organelle_status,
        "matched_organelle_contigs": ["chrM"] if organelle_status == "ok" else [],
        "blacklist_mode": blacklist_mode,
        "blacklist_applied": blacklist_applied,
        "input_stages": {
            "fragment_size": "analysis_bam",
            "mnase_classification": "not_applicable",
            "duplicate_metrics": "pre_duplicate_markdup_bam",
            "nrf_pbc": "pre_duplicate_mapq_bam",
            "organelle": "pre_duplicate_mapq_bam",
            "blacklist": (
                "final_bam_to_analysis_bam" if blacklist_applied else "not_applicable"
            ),
            "peaks": "analysis_bam",
        },
        "stage_counts": {
            "analysis_input_alignments": 120,
            "analysis_kept_alignments": 100,
            "analysis_valid_templates": 45,
            "analysis_invalid_templates": 5,
            "upstream_dropped_templates": 10,
            "upstream_dropped_alignments": 20,
            "blacklist_dropped_templates": 2,
            "blacklist_dropped_alignments": 4,
            "organelle_reads": 3 if organelle_status == "ok" else None,
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


def _write_policy(path: Path, **overrides) -> None:
    path.write_text(
        json.dumps(_policy(**overrides), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_mito(path: Path, *, status="ok") -> None:
    row = {
        "schema_version": "2",
        "sample": "S1",
        "organelle_reads": "3" if status == "ok" else "NA",
        "total_mapped_reads": "100",
        "fraction": "0.030000" if status == "ok" else "NA",
        "mode": "report",
        "status": status,
        "input_stage": "pre_duplicate_mapq_bam",
        "matched_contigs": "chrM" if status == "ok" else "NA",
        "organelle_chroms_explicit": "false",
    }
    _write_tsv(
        path, MITO_QC_V2_COLUMNS, [[row[column] for column in MITO_QC_V2_COLUMNS]]
    )


def _write_fragment_qc(path: Path) -> None:
    columns = (
        "schema_version",
        "sample",
        "layout",
        "fragment_count",
        "fragment_mean",
        "fragment_median",
        "fragment_mode",
        "fragment_min",
        "fragment_max",
        "fraction_lt_150",
        "fraction_150_300",
        "fraction_300_500",
        "fraction_ge_500",
        "fraction_lt_120",
        "status",
        "input_stage",
        "duplicate_analysis",
        "blacklist_applied",
        "templates_total",
        "templates_valid",
        "templates_invalid",
        "invalid_tlen_zero",
        "invalid_orphan",
        "invalid_cross_contig",
        "invalid_inconsistent",
        "alignments_total",
        "upstream_dropped_templates",
        "upstream_dropped_alignments",
        "upstream_alignments_total",
        "upstream_alignments_kept",
    )
    row = {
        "schema_version": "2",
        "sample": "S1",
        "layout": "PE",
        "fragment_count": "45",
        "status": "ok",
        "input_stage": "analysis_bam",
        "duplicate_analysis": "exclude",
        "blacklist_applied": "true",
        "templates_total": "50",
        "templates_valid": "45",
        "templates_invalid": "5",
        "invalid_tlen_zero": "1",
        "invalid_orphan": "2",
        "invalid_cross_contig": "1",
        "invalid_inconsistent": "1",
        "alignments_total": "100",
        "upstream_dropped_templates": "10",
        "upstream_dropped_alignments": "20",
        "upstream_alignments_total": "120",
        "upstream_alignments_kept": "100",
    }
    _write_tsv(path, columns, [[row.get(column, "NA") for column in columns]])


@pytest.fixture
def component_files(tmp_path):
    paths = {
        "peak_counts": tmp_path / "peak_counts.tsv",
        "frip": tmp_path / "frip.tsv",
        "library_complexity": tmp_path / "library_complexity.tsv",
        "nrf_pbc": tmp_path / "nrf_pbc.tsv",
        "mito": tmp_path / "mito.tsv",
        "policy": tmp_path / "analysis_policy.json",
        "fragment": tmp_path / "fragment.tsv",
    }
    _write_tsv(paths["peak_counts"], PEAK_HEADER, [["S1", "narrow", "15000", "12000"]])
    _write_tsv(
        paths["frip"],
        FRIP_HEADER,
        [["S1", "10000000", "2500000", "0.250000", "/bam", "/peaks"]],
    )
    _write_tsv(
        paths["library_complexity"],
        LIBRARY_HEADER,
        [
            [
                "S1",
                "picard_markduplicates",
                "1000",
                "4500000",
                "50000",
                "200000",
                "50000",
                "400000",
                "5000",
                "0.100000",
                "8500000",
                "9001000",
                "450000",
            ]
        ],
    )
    _write_tsv(
        paths["nrf_pbc"],
        NRF_HEADER,
        [
            [
                "S1",
                "9000000",
                "7500000",
                "6000000",
                "1500000",
                "0.833333",
                "0.800000",
                "4.000000",
            ]
        ],
    )
    _write_mito(paths["mito"])
    _write_policy(paths["policy"])
    _write_fragment_qc(paths["fragment"])
    return paths


def _invoke(monkeypatch, program: str, main, args: list[str]) -> None:
    monkeypatch.setattr(sys, "argv", [program, *args])
    main()


def _assemble(
    tmp_path: Path,
    monkeypatch,
    component_files,
    *,
    has_blacklist: str = "yes",
    input_overrides: dict[str, Path] | None = None,
    include_fragment: bool = True,
) -> Path:
    inputs = {**component_files, **(input_overrides or {})}
    output = tmp_path / "qc_summary.tsv"
    args = [
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
        has_blacklist,
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
        str(inputs["mito"]),
        "--analysis-policy",
        str(inputs["policy"]),
        "--fragment-qc",
        str(inputs["fragment"]) if include_fragment else "",
        "--output",
        str(output),
    ]
    _invoke(monkeypatch, "assemble_qc_summary.py", assemble_qc_summary.main, args)
    return output


def _read_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def test_complete_inputs_write_locked_v2_contract_with_v1_prefix(
    tmp_path,
    monkeypatch,
    component_files,
):
    output = _assemble(tmp_path, monkeypatch, component_files)
    header, rows = _read_rows(output)

    assert tuple(assemble_qc_summary._QC_SUMMARY_COLUMNS) == ADAPTER_QC_V2_HEADER
    assert header == list(SAMPLE_QC_V2_COLUMNS)
    assert tuple(header[1 : 1 + len(SAMPLE_QC_V1_COLUMNS)]) == SAMPLE_QC_V1_COLUMNS
    assert len(header) == 64
    assert len(rows) == 1
    row = rows[0]
    assert row["schema_version"] == "2"
    assert row["final_bam"] == "results/S1/02_align/S1.final.bam"
    assert row["analysis_bam"] == "results/S1/02_align/S1.analysis.bam"
    assert row["analysis_valid_fragments"] == "45"
    assert row["analysis_invalid_templates"] == "5"
    assert row["raw_peak_count"] == "15000"
    assert row["delivered_peak_count"] == "12000"
    assert row["organelle_reads"] == "3"
    assert row["organelle_total_mapped_reads"] == "100"
    assert row["organelle_fraction"] == "0.030000"
    assert row["duplicate_analysis"] == "exclude"
    assert row["resolved_max_insert"] == "2000"
    assert row["blacklist_applied"] == "true"
    assert row["cross_correlation_status"] == "not_enabled"
    assert row["scientific_policy_revision"] == "5"
    assert b"\r\n" not in output.read_bytes()


def test_no_blacklist_forces_all_blacklist_fields_to_na(
    tmp_path,
    monkeypatch,
    component_files,
):
    _write_policy(
        component_files["policy"],
        blacklist_mode="off",
        blacklist_applied=False,
    )
    output = _assemble(
        tmp_path,
        monkeypatch,
        component_files,
        has_blacklist="no",
    )
    _, rows = _read_rows(output)

    assert {
        column: rows[0][column]
        for column in (
            "blacklist",
            "blacklist_filtered_bam",
            "blacklist_filtered_peaks",
            "blacklist_filtered_peak_count",
        )
    } == {
        "blacklist": "NA",
        "blacklist_filtered_bam": "NA",
        "blacklist_filtered_peaks": "NA",
        "blacklist_filtered_peak_count": "NA",
    }
    assert rows[0]["delivered_peak_count"] == "15000"
    assert rows[0]["blacklist_applied"] == "false"


def test_fallback_library_complexity_preserves_na_values(
    tmp_path,
    monkeypatch,
    component_files,
):
    _write_tsv(
        component_files["library_complexity"],
        LIBRARY_HEADER,
        [["S1", "fallback", *(["NA"] * 11)]],
    )
    output = _assemble(tmp_path, monkeypatch, component_files)
    _, rows = _read_rows(output)

    assert rows[0]["metrics_source"] == "fallback"
    assert rows[0]["unpaired_reads_examined"] == "NA"
    assert rows[0]["estimated_library_size"] == "NA"


def test_zero_peaks_and_reads_preserve_zero_counts_and_na_frip(
    tmp_path,
    monkeypatch,
    component_files,
):
    _write_tsv(
        component_files["peak_counts"], PEAK_HEADER, [["S1", "narrow", "0", "NA"]]
    )
    _write_tsv(
        component_files["frip"],
        FRIP_HEADER,
        [["S1", "0", "0", "NA", "/bam", "/peaks"]],
    )
    output = _assemble(
        tmp_path,
        monkeypatch,
        component_files,
        has_blacklist="no",
    )
    _, rows = _read_rows(output)

    assert rows[0]["peak_count"] == "0"
    assert rows[0]["total_reads"] == "0"
    assert rows[0]["frip"] == "NA"


def test_mito_status_na_is_not_zero(
    tmp_path,
    monkeypatch,
    component_files,
):
    _write_mito(component_files["mito"], status="no_mapped_reads")
    _write_policy(component_files["policy"], organelle_status="no_mapped_reads")
    output = _assemble(tmp_path, monkeypatch, component_files)
    _, rows = _read_rows(output)

    assert rows[0]["organelle_status"] == "no_mapped_reads"
    assert rows[0]["organelle_reads"] == "NA"
    assert rows[0]["organelle_fraction"] == "NA"


def test_missing_component_file_fails_with_clear_error(
    tmp_path,
    monkeypatch,
    component_files,
    capsys,
):
    missing = tmp_path / "missing-frip.tsv"

    with pytest.raises(SystemExit) as raised:
        _assemble(
            tmp_path,
            monkeypatch,
            component_files,
            input_overrides={"frip": missing},
        )

    assert raised.value.code == 1
    assert "required input file not found" in capsys.readouterr().err


def test_aggregate_without_inputs_writes_header_only_lf_file(tmp_path, monkeypatch):
    output = tmp_path / "aggregate.tsv"
    _invoke(
        monkeypatch,
        "aggregate_qc_summary.py",
        aggregate_qc_summary.main,
        ["--output", str(output)],
    )

    header, rows = _read_rows(output)
    assert header == list(SAMPLE_QC_V2_COLUMNS)
    assert rows == []
    assert b"\r\n" not in output.read_bytes()


def test_aggregate_preserves_all_v1_inputs_as_v1(tmp_path, monkeypatch):
    source = tmp_path / "v1.tsv"
    output = tmp_path / "aggregate.tsv"
    values = {column: "NA" for column in SAMPLE_QC_V1_COLUMNS}
    values["sample"] = "S1"
    _write_tsv(
        source, SAMPLE_QC_V1_COLUMNS, [[values[c] for c in SAMPLE_QC_V1_COLUMNS]]
    )

    _invoke(
        monkeypatch,
        "aggregate_qc_summary.py",
        aggregate_qc_summary.main,
        ["--output", str(output), str(source)],
    )

    header, rows = _read_rows(output)
    assert header == list(SAMPLE_QC_V1_COLUMNS)
    assert rows[0]["sample"] == "S1"


def test_aggregate_mixed_versions_pads_v1_rows_to_v2(tmp_path, monkeypatch):
    v1 = tmp_path / "v1.tsv"
    v2 = tmp_path / "v2.tsv"
    output = tmp_path / "aggregate.tsv"
    v1_values = {column: "NA" for column in SAMPLE_QC_V1_COLUMNS}
    v1_values["sample"] = "S1"
    _write_tsv(v1, SAMPLE_QC_V1_COLUMNS, [[v1_values[c] for c in SAMPLE_QC_V1_COLUMNS]])
    v2_values = {column: "NA" for column in SAMPLE_QC_V2_COLUMNS}
    v2_values.update({"schema_version": "2", "sample": "S2"})
    _write_tsv(v2, SAMPLE_QC_V2_COLUMNS, [[v2_values[c] for c in SAMPLE_QC_V2_COLUMNS]])

    _invoke(
        monkeypatch,
        "aggregate_qc_summary.py",
        aggregate_qc_summary.main,
        ["--output", str(output), str(v1), str(v2)],
    )

    header, rows = _read_rows(output)
    assert header == list(SAMPLE_QC_V2_COLUMNS)
    assert [row["sample"] for row in rows] == ["S1", "S2"]
    assert rows[0]["schema_version"] == "NA"
    assert rows[1]["schema_version"] == "2"


def test_aggregate_rejects_mismatched_header(tmp_path, monkeypatch, capsys):
    bad = tmp_path / "bad.tsv"
    output = tmp_path / "aggregate.tsv"
    _write_tsv(bad, ["wrong", "header"], [["x", "y"]])

    with pytest.raises(SystemExit) as raised:
        _invoke(
            monkeypatch,
            "aggregate_qc_summary.py",
            aggregate_qc_summary.main,
            ["--output", str(output), str(bad)],
        )

    assert raised.value.code == 1
    assert "header mismatch" in capsys.readouterr().err.lower()


def test_aggregate_concatenates_inputs_in_order_without_duplicate_headers(
    tmp_path,
    monkeypatch,
):
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    output = tmp_path / "aggregate.tsv"
    for path, sample in ((first, "S1"), (second, "S2")):
        empty_values = {column: "NA" for column in ADAPTER_QC_V2_HEADER}
        empty_values.update({"schema_version": "2", "sample": sample})
        _write_tsv(
            path,
            ADAPTER_QC_V2_HEADER,
            [[empty_values[column] for column in ADAPTER_QC_V2_HEADER]],
        )

    _invoke(
        monkeypatch,
        "aggregate_qc_summary.py",
        aggregate_qc_summary.main,
        ["--output", str(output), str(first), str(second)],
    )

    header, rows = _read_rows(output)
    assert header == list(ADAPTER_QC_V2_HEADER)
    assert [row["sample"] for row in rows] == ["S1", "S2"]
    assert output.read_text(encoding="utf-8").count("sample\tassay\t") == 1


def test_mnase_v2_header_preserves_v1_prefix_and_adds_union_fields():
    from encode_pipeline.qc.contract import MNASE_QC_V1_COLUMNS

    assert tuple(MNASE_QC_V2_COLUMNS[:1]) == ("schema_version",)
    assert MNASE_QC_V2_COLUMNS[1 : 1 + len(MNASE_QC_V1_COLUMNS)] == MNASE_QC_V1_COLUMNS
    assert "analysis_unclassified_templates" in MNASE_QC_V2_COLUMNS
    assert "analysis_fraction_classified" in MNASE_QC_V2_COLUMNS
    assert "window_overlap" in MNASE_QC_V2_COLUMNS
