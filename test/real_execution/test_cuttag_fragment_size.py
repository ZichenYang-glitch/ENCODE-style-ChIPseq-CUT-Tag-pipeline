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


def _run_report(
    sample,
    bam_path,
    layout,
    output_path,
    samtools_path=None,
    *,
    stage_metrics=None,
    duplicate_analysis="retain",
    blacklist_applied="false",
):
    command = [
        sys.executable,
        SCRIPT,
        "--sample",
        sample,
        "--bam",
        bam_path,
        "--layout",
        layout,
        "--duplicate-analysis",
        duplicate_analysis,
        "--blacklist-applied",
        blacklist_applied,
        "--output",
        output_path,
    ]
    if stage_metrics is not None:
        command[command.index("--output") : command.index("--output")] = [
            "--stage-metrics",
            str(stage_metrics),
        ]
    _run(command, samtools_path=samtools_path)
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
                    f"r{index}\t99\tchr1\t{position}\t60\t100M\t=\t"
                    f"{position + fragment_length}\t{fragment_length}\t*\t*"
                ),
                (
                    f"r{index}\t147\tchr1\t{position + fragment_length}\t"
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
        "r0\t99\tchr1\t100\t60\t100M\t=\t250\t150\t*\t*",
        "r0\t147\tchr1\t250\t60\t100M\t=\t100\t-150\t*\t*",
        "r0\t355\tchr1\t100\t0\t100M\t*\t0\t0\t*\t*",
        "r2\t99\tchr1\t500\t60\t100M\t=\t600\t100\t*\t*",
        "r2\t147\tchr1\t600\t60\t100M\t=\t500\t-100\t*\t*",
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
            "r0\t99\tchr1\t100\t60\t100M\t=\t250\t150\t*\t*",
            "r0\t147\tchr1\t250\t60\t100M\t=\t100\t-150\t*\t*",
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


def test_invalid_templates_never_enter_biological_length_bins(tmp_path, samtools_path):
    reads = [
        # valid 200 bp
        "valid\t99\tchr1\t101\t60\t50M\t=\t251\t200\t*\t*",
        "valid\t147\tchr1\t251\t60\t50M\t=\t101\t-200\t*\t*",
        # complete proper pair with zero TLEN
        "zero\t99\tchr1\t501\t60\t50M\t=\t501\t0\t*\t*",
        "zero\t147\tchr1\t501\t60\t50M\t=\t501\t0\t*\t*",
        # orphan
        "orphan\t99\tchr1\t701\t60\t50M\t=\t801\t100\t*\t*",
        # cross-contig mates
        "cross\t99\tchr1\t901\t60\t50M\tchr2\t1001\t100\t*\t*",
        "cross\t147\tchr2\t1001\t60\t50M\tchr1\t901\t-100\t*\t*",
        # inconsistent reciprocal coordinates
        "bad\t99\tchr1\t1101\t60\t50M\t=\t1301\t200\t*\t*",
        "bad\t147\tchr1\t1302\t60\t50M\t=\t1101\t-200\t*\t*",
    ]
    bam_path = _make_bam(
        tmp_path,
        SAM_HEADER + ["@SQ\tSN:chr2\tLN:10000"],
        samtools_path,
    )
    sam_path = tmp_path / "test.sam"
    sam_path.write_text(
        "\n".join(SAM_HEADER + ["@SQ\tSN:chr2\tLN:10000"] + reads) + "\n",
        encoding="utf-8",
    )
    _run(
        [samtools_path, "view", "-bS", "-o", bam_path, sam_path],
        samtools_path=samtools_path,
    )
    row = _run_report(
        "sample-1",
        bam_path,
        "PE",
        tmp_path / "report.tsv",
        samtools_path,
    )

    assert row["fragment_count"] == "1"
    assert row["fragment_min"] == "200"
    assert row["fragment_max"] == "200"
    assert row["templates_total"] == "5"
    assert row["templates_invalid"] == "4"
    assert row["invalid_tlen_zero"] == "1"
    assert row["invalid_orphan"] == "1"
    assert row["invalid_cross_contig"] == "1"
    assert row["invalid_inconsistent"] == "1"
    assert row["fraction_lt_150"] == "0.0000"


def test_retain_policy_counts_duplicate_flagged_template_once(tmp_path, samtools_path):
    reads = [
        "dup\t1123\tchr1\t101\t60\t50M\t=\t301\t200\t*\t*",
        "dup\t1171\tchr1\t301\t60\t50M\t=\t101\t-200\t*\t*",
    ]
    bam_path = _make_bam(tmp_path, SAM_HEADER + reads, samtools_path)
    row = _run_report(
        "sample-1",
        bam_path,
        "PE",
        tmp_path / "report.tsv",
        samtools_path,
        duplicate_analysis="retain",
    )

    assert row["duplicate_analysis"] == "retain"
    assert row["fragment_count"] == "1"
    assert row["templates_valid"] == "1"


def test_stage_metrics_report_upstream_drops_separately(tmp_path, samtools_path):
    import json

    bam_path = _make_bam(
        tmp_path,
        SAM_HEADER
        + [
            "valid\t99\tchr1\t101\t60\t50M\t=\t301\t200\t*\t*",
            "valid\t147\tchr1\t301\t60\t50M\t=\t101\t-200\t*\t*",
        ],
        samtools_path,
    )
    stage = tmp_path / "stage.json"
    stage.write_text(
        json.dumps(
            {
                "total_alignments": 10,
                "kept_alignments": 2,
                "dropped_alignments": 8,
                "dropped_templates": 4,
                "dropped_alignments_by_reason": {"invalid": 8},
                "dropped_templates_by_reason": {"invalid": 4},
            }
        ),
        encoding="utf-8",
    )
    row = _run_report(
        "sample-1",
        bam_path,
        "PE",
        tmp_path / "report.tsv",
        samtools_path,
        stage_metrics=stage,
    )

    assert row["alignments_total"] == "2"
    assert row["templates_total"] == "1"
    assert row["templates_invalid"] == "0"
    assert row["upstream_dropped_templates"] == "4"
    assert row["upstream_dropped_alignments"] == "8"
    assert row["upstream_alignments_total"] == "10"
    assert row["upstream_alignments_kept"] == "2"


def test_fragment_report_header_is_exact_v2(tmp_path, samtools_path):
    from encode_pipeline.qc.contract import FRAGMENT_QC_V2_COLUMNS

    bam_path = _make_bam(
        tmp_path,
        SAM_HEADER
        + [
            "valid\t99\tchr1\t101\t60\t50M\t=\t301\t200\t*\t*",
            "valid\t147\tchr1\t301\t60\t50M\t=\t101\t-200\t*\t*",
        ],
        samtools_path,
    )
    output = tmp_path / "report.tsv"
    _run_report("sample-1", bam_path, "PE", output, samtools_path)
    with output.open(newline="", encoding="utf-8") as handle:
        header = tuple(next(csv.reader(handle, delimiter="\t")))
    assert header == FRAGMENT_QC_V2_COLUMNS
