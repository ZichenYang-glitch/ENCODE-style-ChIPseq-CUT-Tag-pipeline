"""Real-tool contracts for CUT&Tag fragment-size reporting."""

import csv
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from _tool_resolver import require_external_tools, resolve_tool

TEST_ROOT = Path(__file__).resolve().parents[1]


REPO_ROOT = TEST_ROOT.parent
SCRIPT = REPO_ROOT / "scripts" / "calc_cuttag_fragment_size.py"
SAM_HEADER = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:10000"]

pytestmark = pytest.mark.real_execution


@pytest.fixture(scope="module")
def samtools_path():
    candidate = resolve_tool("samtools", "SAMTOOLS")
    resolved = shutil.which(candidate)
    if not resolved and os.path.isabs(candidate):
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            resolved = candidate
    if not resolved:
        require_external_tools(
            ["samtools"],
            "CUT&Tag fragment-size real execution",
        )
    return resolved


def _run(command, *, samtools_path=None):
    env = os.environ.copy()
    if samtools_path:
        env["SAMTOOLS"] = samtools_path
    result = subprocess.run(command, capture_output=True, text=True, env=env)
    assert result.returncode == 0, (
        f"Command failed: {command!r}\n"
        f"stdout:\n{result.stdout[-500:]}\n"
        f"stderr:\n{result.stderr[-500:]}"
    )
    return result


def _make_bam(tmp_path, sam_lines, samtools_path, bam_name="test.bam"):
    sam_path = tmp_path / f"{Path(bam_name).stem}.sam"
    bam_path = tmp_path / bam_name
    sam_path.write_text("\n".join(sam_lines) + "\n", encoding="utf-8")
    _run(
        [samtools_path, "view", "-bS", "-o", bam_path, sam_path],
        samtools_path=samtools_path,
    )
    return bam_path


def _run_report(sample, bam_path, layout, output_path, samtools_path=None):
    _run(
        [
            sys.executable,
            SCRIPT,
            "--sample",
            sample,
            "--bam",
            bam_path,
            "--layout",
            layout,
            "--output",
            output_path,
        ],
        samtools_path=samtools_path,
    )
    with output_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    return rows[0]


def _paired_reads(fragment_lengths):
    reads = []
    for index, fragment_length in enumerate(fragment_lengths):
        position = (index + 1) * 100
        reads.extend(
            [
                (
                    f"r{index}/1\t99\tchr1\t{position}\t60\t100M\t=\t"
                    f"{position + fragment_length}\t{fragment_length}\t*\t*"
                ),
                (
                    f"r{index}/2\t147\tchr1\t{position + fragment_length}\t"
                    f"60\t100M\t=\t{position}\t{-fragment_length}\t*\t*"
                ),
            ]
        )
    return reads


def _exclusive_fraction_sum(row):
    return sum(
        float(row[field])
        for field in (
            "fraction_lt_150",
            "fraction_150_300",
            "fraction_300_500",
            "fraction_ge_500",
        )
    )


def test_paired_end_report_summarizes_fragments_and_exclusive_bins(
    tmp_path, samtools_path
):
    bam_path = _make_bam(
        tmp_path,
        SAM_HEADER + _paired_reads([100, 200, 350, 450, 600]),
        samtools_path,
    )
    row = _run_report(
        "sample-1",
        bam_path,
        "PE",
        tmp_path / "report.tsv",
        samtools_path,
    )

    assert row["fragment_count"] == "5"
    assert row["fragment_median"] == "350.0"
    assert _exclusive_fraction_sum(row) == pytest.approx(1.0, abs=0.0002)
    assert row["fraction_lt_120"] == "0.2000"
    assert row["status"] == "ok"


def test_single_end_report_is_explicitly_not_supported_without_calling_samtools(
    tmp_path,
):
    row = _run_report(
        "sample-1",
        tmp_path / "unused.bam",
        "SE",
        tmp_path / "report.tsv",
    )

    assert row["status"] == "layout_not_supported"
    assert row["fragment_count"] == "NA"
    assert all(
        row[field] == "NA"
        for field in (
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
        )
    )


def test_paired_end_report_identifies_no_fragments(tmp_path, samtools_path):
    bam_path = _make_bam(
        tmp_path,
        SAM_HEADER + ["unmapped\t77\t*\t0\t0\t*\t*\t0\t0\t*\t*"],
        samtools_path,
    )
    row = _run_report(
        "sample-1",
        bam_path,
        "PE",
        tmp_path / "report.tsv",
        samtools_path,
    )

    assert row["status"] == "no_fragments"
    assert row["fragment_count"] == "0"


def test_secondary_alignments_and_read2_records_are_not_counted(
    tmp_path, samtools_path
):
    reads = [
        "r0/1\t99\tchr1\t100\t60\t100M\t=\t250\t150\t*\t*",
        "r0/2\t147\tchr1\t250\t60\t100M\t=\t100\t-150\t*\t*",
        "r0/1\t355\tchr1\t100\t0\t100M\t*\t0\t0\t*\t*",
        "r2/1\t99\tchr1\t500\t60\t100M\t=\t600\t100\t*\t*",
        "r2/2\t147\tchr1\t600\t60\t100M\t=\t500\t-100\t*\t*",
    ]
    bam_path = _make_bam(tmp_path, SAM_HEADER + reads, samtools_path)
    row = _run_report(
        "sample-1",
        bam_path,
        "PE",
        tmp_path / "report.tsv",
        samtools_path,
    )

    assert row["fragment_count"] == "2"


def test_one_proper_pair_contributes_one_fragment(tmp_path, samtools_path):
    bam_path = _make_bam(
        tmp_path,
        SAM_HEADER
        + [
            "r0/1\t99\tchr1\t100\t60\t100M\t=\t250\t150\t*\t*",
            "r0/2\t147\tchr1\t250\t60\t100M\t=\t100\t-150\t*\t*",
        ],
        samtools_path,
    )
    row = _run_report(
        "sample-1",
        bam_path,
        "PE",
        tmp_path / "report.tsv",
        samtools_path,
    )

    assert row["fragment_count"] == "1"


def test_fragment_bins_are_exclusive_and_diagnostic_subset_is_bounded(
    tmp_path, samtools_path
):
    bam_path = _make_bam(
        tmp_path,
        SAM_HEADER + _paired_reads([50, 120, 200, 250, 350, 450, 550, 650]),
        samtools_path,
    )
    row = _run_report(
        "sample-1",
        bam_path,
        "PE",
        tmp_path / "report.tsv",
        samtools_path,
    )

    assert _exclusive_fraction_sum(row) == pytest.approx(1.0, abs=0.0002)
    assert float(row["fraction_lt_120"]) <= float(row["fraction_lt_150"])
