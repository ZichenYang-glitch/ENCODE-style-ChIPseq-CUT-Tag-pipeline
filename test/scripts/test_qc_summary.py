"""Producer contracts for per-sample and aggregate QC summary TSVs."""

from __future__ import annotations

import csv
from pathlib import Path
import sys

import pytest

from encode_pipeline.adapters.encode_qc import _QC_HEADER as ADAPTER_QC_HEADER
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


@pytest.fixture
def component_files(tmp_path):
    paths = {
        "peak_counts": tmp_path / "peak_counts.tsv",
        "frip": tmp_path / "frip.tsv",
        "library_complexity": tmp_path / "library_complexity.tsv",
        "nrf_pbc": tmp_path / "nrf_pbc.tsv",
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


def test_complete_inputs_write_locked_37_column_contract(
    tmp_path,
    monkeypatch,
    component_files,
):
    output = _assemble(tmp_path, monkeypatch, component_files)
    header, rows = _read_rows(output)

    assert tuple(assemble_qc_summary._QC_SUMMARY_COLUMNS) == ADAPTER_QC_HEADER
    assert (
        aggregate_qc_summary._QC_SUMMARY_COLUMNS
        == assemble_qc_summary._QC_SUMMARY_COLUMNS
    )
    assert header == list(ADAPTER_QC_HEADER)
    assert len(header) == 37
    assert len(rows) == 1
    assert {
        key: rows[0][key]
        for key in (
            "sample",
            "frip",
            "peak_count",
            "blacklist_filtered_peak_count",
            "nrf",
        )
    } == {
        "sample": "S1",
        "frip": "0.250000",
        "peak_count": "15000",
        "blacklist_filtered_peak_count": "12000",
        "nrf": "0.833333",
    }
    assert b"\r\n" not in output.read_bytes()


def test_no_blacklist_forces_all_blacklist_fields_to_na(
    tmp_path,
    monkeypatch,
    component_files,
):
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
    assert header == list(ADAPTER_QC_HEADER)
    assert rows == []
    assert b"\r\n" not in output.read_bytes()


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
    empty_values = {column: "NA" for column in ADAPTER_QC_HEADER}
    _write_tsv(
        first,
        ADAPTER_QC_HEADER,
        [[{**empty_values, "sample": "S1"}[column] for column in ADAPTER_QC_HEADER]],
    )
    _write_tsv(
        second,
        ADAPTER_QC_HEADER,
        [[{**empty_values, "sample": "S2"}[column] for column in ADAPTER_QC_HEADER]],
    )

    _invoke(
        monkeypatch,
        "aggregate_qc_summary.py",
        aggregate_qc_summary.main,
        ["--output", str(output), str(first), str(second)],
    )

    header, rows = _read_rows(output)
    assert header == list(ADAPTER_QC_HEADER)
    assert [row["sample"] for row in rows] == ["S1", "S2"]
    assert output.read_text(encoding="utf-8").count("sample\tassay\t") == 1
