"""Tests for blacklist interval helpers and final peak filtering."""

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

_SCRIPT = Path(__file__).resolve().parents[2] / "scripts/blacklist_filter.py"
_SPEC = importlib.util.spec_from_file_location("blacklist_filter", _SCRIPT)
blacklist_filter = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(blacklist_filter)

REPO_ROOT = Path(__file__).resolve().parents[2]
NARROW_ROW = "chr1\t0\t100\tpeak1\t100\t.\t1\t2\t3\t50\n"


def test_blacklist_merges_and_checks_boundaries(tmp_path):
    path = tmp_path / "blacklist.bed"
    path.write_text(
        "chr1\t10\t20\nchr1\t15\t30\nchr2\t0\t5\ntrack name=test\n",
        encoding="utf-8",
    )
    blacklist = blacklist_filter.Blacklist(path)

    assert blacklist.overlaps("chr1", 0, 10) is False
    assert blacklist.overlaps("chr1", 0, 11) is True
    assert blacklist.overlaps("chr1", 30, 40) is False
    assert blacklist.overlaps("chr2", 4, 6) is True
    assert blacklist.overlaps("chr3", 0, 100) is False


def test_filter_peak_file_keeps_non_overlapping_rows(tmp_path):
    blacklist_path = tmp_path / "blacklist.bed"
    blacklist_path.write_text("chr1\t0\t10\n", encoding="utf-8")
    peaks = tmp_path / "peaks.narrowPeak"
    peaks.write_text(
        "chr1\t0\t100\tbad\t100\t.\t1\t2\t3\t50\n"
        "chr1\t50\t150\tgood\t100\t.\t1\t2\t3\t100\n",
        encoding="utf-8",
    )
    output = tmp_path / "filtered.narrowPeak"

    kept, dropped = blacklist_filter.Blacklist(blacklist_path).filter_peak_file(
        peaks, output
    )

    assert (kept, dropped) == (1, 1)
    assert "good" in output.read_text(encoding="utf-8")
    assert "bad" not in output.read_text(encoding="utf-8")


def _run_consensus(tmp_path, peak_text):
    peaks_a = tmp_path / "a.narrowPeak"
    peaks_b = tmp_path / "b.narrowPeak"
    peaks_a.write_text(peak_text, encoding="utf-8")
    peaks_b.write_text(peak_text, encoding="utf-8")
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr1\t0\t10\n", encoding="utf-8")
    output = tmp_path / "consensus.narrowPeak"
    summary = tmp_path / "summary.tsv"

    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts/compute_consensus.py"),
            "--peaks",
            str(peaks_a),
            str(peaks_b),
            "--bioreps",
            "1",
            "2",
            "--format",
            "narrowPeak",
            "--min-replicates",
            "2",
            "--reciprocal-overlap",
            "0.5",
            "--output",
            str(output),
            "--summary",
            str(summary),
            "--experiment",
            "exp1",
            "--assay",
            "chipseq",
            "--caller",
            "macs3",
            "--peak-mode",
            "narrow",
            "--blacklist",
            str(blacklist),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    with open(summary, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return output, rows[0]


def test_compute_consensus_filters_final_blacklist_peaks(tmp_path):
    output, summary = _run_consensus(tmp_path, NARROW_ROW)

    assert output.read_text(encoding="utf-8") == ""
    assert summary["consensus_peak_count"] == "0"
    distribution = json.loads(summary["support_distribution"])
    assert sum(distribution.values()) == 0


def test_compute_consensus_filters_components_and_support_distribution(tmp_path):
    peak_text = NARROW_ROW + "chr1\t500\t600\tpeak2\t100\t.\t1\t2\t3\t50\n"
    output, summary = _run_consensus(tmp_path, peak_text)

    rows = output.read_text(encoding="utf-8").strip().splitlines()
    assert len(rows) == 1
    assert "500" in rows[0]
    assert summary["consensus_peak_count"] == "1"
    distribution = json.loads(summary["support_distribution"])
    assert sum(distribution.values()) == 1


def test_chipseq_idr_summary_filters_final_peaks(tmp_path):
    peak = tmp_path / "idr.thresholded.narrowPeak"
    peak.write_text(NARROW_ROW, encoding="utf-8")
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr1\t0\t10\n", encoding="utf-8")
    conservative = tmp_path / "conservative.narrowPeak"
    optimal = tmp_path / "optimal.narrowPeak"

    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts/chipseq_idr_summary.py"),
            "--true-peaks",
            str(peak),
            "--pooled-peaks",
            str(peak),
            "--self1-peaks",
            str(peak),
            "--self2-peaks",
            str(peak),
            "--experiment",
            "exp1",
            "--bio-rep-a",
            "1",
            "--bio-rep-b",
            "2",
            "--output-tsv",
            str(tmp_path / "summary.tsv"),
            "--output-cons",
            str(conservative),
            "--output-opt",
            str(optimal),
            "--blacklist",
            str(blacklist),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert conservative.read_text(encoding="utf-8") == ""
    assert optimal.read_text(encoding="utf-8") == ""


def test_idr_reproducibility_summary_filters_final_peak(tmp_path):
    peak = tmp_path / "idr.thresholded.narrowPeak"
    peak.write_text(NARROW_ROW, encoding="utf-8")
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr1\t0\t10\n", encoding="utf-8")
    final_peak = tmp_path / "final.narrowPeak"

    result = subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts/idr_reproducibility_summary.py"),
            "--true-peaks",
            str(peak),
            "--pooled-peaks",
            str(peak),
            "--self1-peaks",
            str(peak),
            "--self2-peaks",
            str(peak),
            "--experiment",
            "exp1",
            "--assay",
            "atac",
            "--caller",
            "macs3",
            "--peak-mode",
            "narrow",
            "--bio-rep-a",
            "1",
            "--bio-rep-b",
            "2",
            "--final-method",
            "idr",
            "--final-output",
            str(final_peak),
            "--output-tsv",
            str(tmp_path / "summary.tsv"),
            "--output-peak",
            str(final_peak),
            "--blacklist",
            str(blacklist),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert final_peak.read_text(encoding="utf-8") == ""
