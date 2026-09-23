"""Behavioral tests for the unified IDR reproducibility summary producer."""

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest

_SCRIPT = Path(__file__).resolve().parents[2] / "scripts/idr_reproducibility_summary.py"
_SPEC = importlib.util.spec_from_file_location("idr_reproducibility_summary", _SCRIPT)
idr_summary = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(idr_summary)


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


CHIPSEQ_SUMMARY_HEADER = (
    "experiment",
    "true_peaks(Nt)",
    "pooled_peaks(Np)",
    "self1_peaks(N1)",
    "self2_peaks(N2)",
    "rescue_ratio",
    "self_consistency_ratio",
    "reproducibility_status",
)


def _run_summary_cli(tmp_path, script, counts, *, assay="atac", peak_mode="narrow"):
    """Run the real CLI on unique peaks; retain argv, TSV and copied bytes."""
    repo = Path(__file__).resolve().parents[2]
    inputs = {}
    args = []
    for label, count in zip(("true", "pooled", "self1", "self2"), counts, strict=True):
        path = tmp_path / f"{label}.peaks"
        _write_peaks(path, count, peak_mode)
        inputs[label] = path
        args.extend([f"--{label}-peaks", str(path)])
    summary = tmp_path / "summary.tsv"
    args.extend(
        [
            "--experiment",
            "exp1",
            "--bio-rep-a",
            "1",
            "--bio-rep-b",
            "2",
            "--output-tsv",
            str(summary),
        ]
    )
    chipseq = script == "chipseq_idr_summary.py"
    copies = (
        {"conservative.peaks": "true", "optimal.peaks": "pooled"}
        if chipseq
        else {"final.peaks": "true"}
    )
    final_output = "results/final peaks.idr.narrowPeak"
    if chipseq:
        args.extend(
            [
                "--output-cons",
                str(tmp_path / "conservative.peaks"),
                "--output-opt",
                str(tmp_path / "optimal.peaks"),
            ]
        )
    else:
        args.extend(
            [
                "--caller",
                "macs3",
                "--final-method",
                "idr",
                "--final-output",
                final_output,
                "--output-peak",
                str(tmp_path / "final.peaks"),
            ]
        )
        # Exercise wrapper defaults, not explicit flags that bypass them.
        if script in ("idr_reproducibility_summary.py", "broad_idr_summary.py"):
            args.extend(["--assay", assay])
        if script == "idr_reproducibility_summary.py":
            args.extend(["--peak-mode", peak_mode])
    command = [sys.executable, "-B", "-S", str(repo / "scripts" / script), *args]
    result = subprocess.run(command, cwd=tmp_path, capture_output=True, text=True)
    (tmp_path / "cli.stdout").write_text(result.stdout, encoding="utf-8")
    (tmp_path / "cli.stderr").write_text(result.stderr, encoding="utf-8")
    (tmp_path / "command.json").write_text(
        json.dumps(
            {"argv": command, "cwd": str(tmp_path), "exit_code": result.returncode},
            indent=2,
        ),
        encoding="utf-8",
    )
    assert result.returncode == 0, result.stdout + result.stderr
    for output, source in copies.items():
        assert (tmp_path / output).read_bytes() == inputs[source].read_bytes()
    with summary.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    assert len(rows) == 1
    assert reader.fieldnames == list(
        CHIPSEQ_SUMMARY_HEADER if chipseq else SUMMARY_HEADER
    )
    expected = {"experiment": "exp1"}
    count_keys = CHIPSEQ_SUMMARY_HEADER[1:5] if chipseq else SUMMARY_HEADER[6:10]
    expected.update(zip(count_keys, map(str, counts), strict=True))
    if not chipseq:
        expected.update(
            assay=assay,
            peak_mode=peak_mode,
            caller="macs3",
            bio_rep_a="1",
            bio_rep_b="2",
            final_method="idr",
            final_output=final_output,
        )
    row = rows[0]
    ratios_and_status = {
        key: row.pop(key)
        for key in (
            "rescue_ratio",
            "self_consistency_ratio",
            "reproducibility_status",
        )
    }
    assert row == expected
    return ratios_and_status


@pytest.mark.parametrize(
    "script", ["chipseq_idr_summary.py", "idr_reproducibility_summary.py"]
)
@pytest.mark.parametrize(
    ("counts", "rescue", "self_ratio", "status"),
    [
        pytest.param((5000, 9998, 5, 5), "2.000", "1.000", "pass", id="rescue-below"),
        pytest.param((5, 5, 5000, 9998), "1.000", "2.000", "pass", id="self-below"),
        pytest.param(
            (5000, 9998, 5000, 9998), "2.000", "2.000", "pass", id="both-below"
        ),
        pytest.param((5000, 10000, 5, 5), "2.000", "1.000", "fail", id="rescue-exact"),
        pytest.param((5, 5, 5000, 10000), "1.000", "2.000", "fail", id="self-exact"),
        pytest.param((5000, 10001, 5, 5), "2.000", "1.000", "fail", id="rescue-above"),
        pytest.param((5, 5, 5000, 10001), "1.000", "2.000", "fail", id="self-above"),
        pytest.param((9998, 5000, 5, 5), "2.000", "1.000", "pass", id="rescue-swapped"),
        pytest.param((5, 5, 9998, 5000), "1.000", "2.000", "pass", id="self-swapped"),
        pytest.param(
            (3, 2, 3, 2), "1.500", "1.500", "pass", id="true-larger-than-pooled"
        ),
        pytest.param((3, 1, 3, 2), "3.000", "1.500", "fail", id="ordinary-fail"),
        pytest.param((0, 0, 0, 0), "NA", "NA", "fail", id="all-empty"),
        pytest.param((0, 0, 2, 2), "NA", "1.000", "fail", id="rescue-na"),
        pytest.param((2, 2, 0, 0), "1.000", "NA", "fail", id="self-na"),
        pytest.param((0, 3, 2, 2), "inf", "1.000", "fail", id="rescue-inf"),
        pytest.param((2, 2, 0, 3), "1.000", "inf", "fail", id="self-inf"),
        pytest.param((3, 0, 2, 2), "inf", "1.000", "fail", id="rescue-inf-swapped"),
        pytest.param((2, 2, 3, 0), "1.000", "inf", "fail", id="self-inf-swapped"),
    ],
)
def test_summary_cli_grades_raw_counts_and_preserves_contract(
    tmp_path, script, counts, rescue, self_ratio, status
):
    row = _run_summary_cli(tmp_path, script, counts)
    assert row == {
        "rescue_ratio": rescue,
        "self_consistency_ratio": self_ratio,
        "reproducibility_status": status,
    }


@pytest.mark.parametrize(
    ("script", "assay", "peak_mode"),
    [
        ("atac_idr_summary.py", "atac", "narrow"),
        ("cuttag_idr_summary.py", "cuttag", "narrow"),
        ("broad_idr_summary.py", "chipseq", "broad"),
    ],
)
def test_summary_wrapper_cli_delegates_grading_and_preserves_peaks(
    tmp_path, script, assay, peak_mode
):
    row = _run_summary_cli(
        tmp_path, script, (9998, 5000, 5, 5), assay=assay, peak_mode=peak_mode
    )
    assert row == {
        "rescue_ratio": "2.000",
        "self_consistency_ratio": "1.000",
        "reproducibility_status": "pass",
    }
