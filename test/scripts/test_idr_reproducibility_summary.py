"""Behavioral tests for the unified IDR reproducibility summary producer."""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from scripts import idr_reproducibility_summary as idr_summary


SUMMARY_HEADER = (
    "experiment",
    "assay",
    "peak_mode",
    "caller",
    "bio_rep_a",
    "bio_rep_b",
    "true_peaks_Nt",
    "pooled_peaks_Np",
    "self1_peaks_N1",
    "self2_peaks_N2",
    "rescue_ratio",
    "self_consistency_ratio",
    "reproducibility_status",
    "final_method",
    "final_output",
)


def _peak_line(index: int, peak_mode: str) -> str:
    start = 100 + index * 100
    end = start + 50
    common = ["chr1", str(start), str(end), f"peak_{index + 1}", "1000", "."]
    if peak_mode == "broad":
        return "\t".join(
            common
            + [
                "5.0",
                "3.0",
                "2.0",
                "0.01",
                "0.005",
                str(start),
                str(end),
                "5.0",
                str(start + 10),
                str(end + 10),
                "4.5",
            ]
        )
    return "\t".join(common + ["5.0", "-1", "-1", "25"])


def _write_peaks(
    path: Path,
    count: int,
    peak_mode: str,
    *,
    prefix_lines: tuple[str, ...] = (),
) -> None:
    lines = [*prefix_lines, *(_peak_line(index, peak_mode) for index in range(count))]
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def _run_summary(
    tmp_path: Path,
    *,
    assay: str = "atac",
    peak_mode: str = "narrow",
    counts: tuple[int, int, int, int] = (10, 12, 9, 10),
    true_prefix_lines: tuple[str, ...] = (),
) -> tuple[list[str], dict[str, str], bytes, bytes]:
    suffix = "broadPeak" if peak_mode == "broad" else "narrowPeak"
    true_peaks = tmp_path / f"true.{suffix}"
    pooled_peaks = tmp_path / f"pooled.{suffix}"
    self1_peaks = tmp_path / f"self1.{suffix}"
    self2_peaks = tmp_path / f"self2.{suffix}"
    output_tsv = tmp_path / "summary.tsv"
    output_peak = tmp_path / f"final.{suffix}"

    _write_peaks(
        true_peaks,
        counts[0],
        peak_mode,
        prefix_lines=true_prefix_lines,
    )
    for path, count in zip(
        (pooled_peaks, self1_peaks, self2_peaks),
        counts[1:],
        strict=True,
    ):
        _write_peaks(path, count, peak_mode)

    final_output = (
        f"results/experiments/exp1/06_reproducibility/final/"
        f"exp1.{assay}.macs3.{peak_mode}.replicate_validated.idr.{suffix}"
    )
    idr_summary.main(
        [
            "--true-peaks",
            str(true_peaks),
            "--pooled-peaks",
            str(pooled_peaks),
            "--self1-peaks",
            str(self1_peaks),
            "--self2-peaks",
            str(self2_peaks),
            "--experiment",
            "exp1",
            "--assay",
            assay,
            "--caller",
            "macs3",
            "--peak-mode",
            peak_mode,
            "--bio-rep-a",
            "1",
            "--bio-rep-b",
            "2",
            "--final-method",
            "idr",
            "--final-output",
            final_output,
            "--output-tsv",
            str(output_tsv),
            "--output-peak",
            str(output_peak),
        ]
    )

    with output_tsv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        header = list(reader.fieldnames or [])
    assert len(rows) == 1
    return header, rows[0], true_peaks.read_bytes(), output_peak.read_bytes()


@pytest.mark.parametrize(
    ("assay", "peak_mode"),
    [
        pytest.param("atac", "narrow", id="atac-narrow"),
        pytest.param("cuttag", "narrow", id="cuttag-narrow"),
        pytest.param("chipseq", "broad", id="chipseq-broad"),
        pytest.param("cuttag", "broad", id="cuttag-broad"),
    ],
)
def test_summary_contract_for_every_supported_mode(tmp_path, assay, peak_mode):
    header, row, true_content, final_content = _run_summary(
        tmp_path,
        assay=assay,
        peak_mode=peak_mode,
    )

    suffix = "broadPeak" if peak_mode == "broad" else "narrowPeak"
    assert header == list(SUMMARY_HEADER)
    assert row == {
        "experiment": "exp1",
        "assay": assay,
        "peak_mode": peak_mode,
        "caller": "macs3",
        "bio_rep_a": "1",
        "bio_rep_b": "2",
        "true_peaks_Nt": "10",
        "pooled_peaks_Np": "12",
        "self1_peaks_N1": "9",
        "self2_peaks_N2": "10",
        "rescue_ratio": "1.200",
        "self_consistency_ratio": "1.111",
        "reproducibility_status": "pass",
        "final_method": "idr",
        "final_output": (
            "results/experiments/exp1/06_reproducibility/final/"
            f"exp1.{assay}.macs3.{peak_mode}.replicate_validated.idr.{suffix}"
        ),
    }
    assert final_content == true_content


@pytest.mark.parametrize(
    ("counts", "rescue_ratio", "self_ratio", "status"),
    [
        pytest.param((0, 3, 0, 0), "inf", "NA", "fail", id="zero-denominators"),
        pytest.param((10, 10, 10, 10), "1.000", "1.000", "pass", id="equal"),
        pytest.param((5, 20, 5, 5), "4.000", "1.000", "fail", id="ratio-fail"),
    ],
)
def test_summary_reports_ratio_edges(
    tmp_path,
    counts,
    rescue_ratio,
    self_ratio,
    status,
):
    _, row, _, _ = _run_summary(tmp_path, counts=counts)

    assert row["rescue_ratio"] == rescue_ratio
    assert row["self_consistency_ratio"] == self_ratio
    assert row["reproducibility_status"] == status


def test_broad_summary_ignores_headers_and_preserves_seventeen_column_peaks(tmp_path):
    prefix = ("# IDR thresholded broadPeak output", "track name=idr")
    _, row, true_content, final_content = _run_summary(
        tmp_path,
        assay="chipseq",
        peak_mode="broad",
        counts=(3, 4, 3, 3),
        true_prefix_lines=prefix,
    )

    assert row["true_peaks_Nt"] == "3"
    assert final_content == true_content
    data_lines = [
        line
        for line in final_content.decode().splitlines()
        if line and not line.startswith(("#", "track"))
    ]
    assert all(len(line.split("\t")) == 17 for line in data_lines)


@pytest.mark.parametrize(
    ("numerator", "denominator", "expected"),
    [(0, 0, "NA"), (3, 0, "inf"), (12, 10, "1.200")],
)
def test_compute_ratio_handles_zero_and_formats_values(
    numerator, denominator, expected
):
    assert idr_summary.compute_ratio(numerator, denominator) == expected
