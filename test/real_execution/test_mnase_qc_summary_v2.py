"""Real-samtools contracts for versioned MNase QC fragment classification."""

from __future__ import annotations

import csv
import json
import os
import shutil
from pathlib import Path
import subprocess
import sys

import pytest

from _tool_resolver import require_external_tools, resolve_tool

TEST_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = TEST_ROOT.parent
SCRIPT = REPO_ROOT / "scripts" / "mnase_qc_summary.py"
SAM_HEADER = [
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:100000",
    "@SQ\tSN:chr2\tLN:100000",
]

pytestmark = pytest.mark.real_execution


@pytest.fixture(scope="module")
def samtools_path():
    candidate = resolve_tool("samtools", "SAMTOOLS")
    resolved = shutil.which(candidate)
    if not resolved and os.path.isabs(candidate):
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            resolved = candidate
    if not resolved:
        require_external_tools(["samtools"], "MNase QC real execution")
    return resolved


def _run(command, *, samtools_path):
    env = os.environ.copy()
    env["SAMTOOLS"] = samtools_path
    env["PATH"] = f"{Path(samtools_path).parent}{os.pathsep}{env.get('PATH', '')}"
    result = subprocess.run(command, capture_output=True, text=True, env=env)
    assert result.returncode == 0, (
        f"Command failed: {command!r}\n"
        f"stdout:\n{result.stdout[-500:]}\n"
        f"stderr:\n{result.stderr[-500:]}"
    )
    return result


def _make_bam(tmp_path, records, samtools_path):
    sam = tmp_path / "input.sam"
    bam = tmp_path / "input.bam"
    sam.write_text("\n".join(SAM_HEADER + records) + "\n", encoding="utf-8")
    _run(
        [samtools_path, "view", "-bS", "-o", bam, sam],
        samtools_path=samtools_path,
    )
    return bam


def _pair(qname, start, length):
    return [
        f"{qname}\t99\tchr1\t{start}\t60\t50M\t=\t{start + length}\t{length}\t*\t*",
        f"{qname}\t147\tchr1\t{start + length}\t60\t50M\t=\t{start}\t{-length}\t*\t*",
    ]


def _run_summary(
    tmp_path,
    samtools_path,
    *,
    bam,
    windows,
    stage_metrics=None,
    analysis_bam=None,
):
    output = tmp_path / "mnase.tsv"
    args = [
        sys.executable,
        SCRIPT,
        "--sample",
        "M1",
        "--assay",
        "mnase",
        "--peak-mode",
        "nucleosome",
        "--sub-min",
        str(windows[0][0]),
        "--sub-max",
        str(windows[0][1]),
        "--mono-min",
        str(windows[1][0]),
        "--mono-max",
        str(windows[1][1]),
        "--di-min",
        str(windows[2][0]),
        "--di-max",
        str(windows[2][1]),
        "--dyad-min",
        "130",
        "--dyad-max",
        "200",
        "--sub-bam",
        str(bam),
        "--mono-bam",
        str(bam),
        "--di-bam",
        str(bam),
        "--analysis-bam",
        str(analysis_bam if analysis_bam is not None else bam),
        "--output",
        str(output),
    ]
    if stage_metrics is not None:
        args[args.index("--output") : args.index("--output")] = [
            "--stage-metrics",
            str(stage_metrics),
        ]
    _run(args, samtools_path=samtools_path)
    with output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    return rows[0]


def test_windows_are_inclusive_and_unclassified_uses_the_union(tmp_path, samtools_path):
    records = []
    for index, length in enumerate((139, 140, 200, 201, 300, 400)):
        records.extend(_pair(f"p{index}", 1000 + index * 1000, length))
    bam = _make_bam(tmp_path, records, samtools_path)
    stage = tmp_path / "stage.json"
    stage.write_text(
        json.dumps(
            {
                "total_alignments": 12,
                "kept_alignments": 12,
                "dropped_alignments": 0,
                "dropped_templates": 0,
                "dropped_alignments_by_reason": {},
                "dropped_templates_by_reason": {},
            }
        ),
        encoding="utf-8",
    )
    row = _run_summary(
        tmp_path,
        samtools_path,
        bam=bam,
        windows=((1, 139), (140, 200), (201, 300)),
        stage_metrics=stage,
    )

    assert row["analysis_status"] == "ok"
    assert row["analysis_total_templates"] == "6"
    assert row["analysis_valid_templates"] == "6"
    assert row["analysis_invalid_templates"] == "0"
    assert row["sub_fragments"] == "1"
    assert row["mono_fragments"] == "2"
    assert row["di_fragments"] == "2"
    assert row["analysis_classified_templates"] == "5"
    assert row["analysis_unclassified_templates"] == "1"
    assert row["analysis_fraction_classified"] == "0.833333"
    assert row["window_overlap"] == "false"
    # Read counts remain read-record counts, not half of the fragment count.
    assert row["sub_reads"] == "12"
    assert row["mono_reads"] == "12"
    assert row["di_reads"] == "12"


def test_overlapping_windows_count_each_window_and_union(tmp_path, samtools_path):
    records = []
    for index, length in enumerate((139, 200, 300, 400, 450)):
        records.extend(_pair(f"p{index}", 1000 + index * 1000, length))
    bam = _make_bam(tmp_path, records, samtools_path)
    row = _run_summary(
        tmp_path,
        samtools_path,
        bam=bam,
        windows=((1, 200), (140, 300), (250, 400)),
    )

    assert row["window_overlap"] == "true"
    assert row["sub_fragments"] == "2"  # 139, 200
    assert row["mono_fragments"] == "2"  # 200, 300
    assert row["di_fragments"] == "2"  # 300, 400
    assert row["analysis_classified_templates"] == "4"  # 139,200,300,400
    assert row["analysis_unclassified_templates"] == "1"  # 450
    assert row["analysis_fraction_classified"] == "0.800000"


def test_invalid_templates_and_zero_tlen_are_counted_separately(
    tmp_path, samtools_path
):
    records = _pair("valid", 1000, 150)
    records += [
        "zero\t99\tchr1\t5000\t60\t50M\t=\t5000\t0\t*\t*",
        "zero\t147\tchr1\t5000\t60\t50M\t=\t5000\t0\t*\t*",
        "orphan\t99\tchr1\t6000\t60\t50M\t=\t6100\t100\t*\t*",
        "cross\t99\tchr1\t7000\t60\t50M\tchr2\t7100\t100\t*\t*",
        "cross\t147\tchr2\t7100\t60\t50M\tchr1\t7000\t-100\t*\t*",
        "bad\t99\tchr1\t8000\t60\t50M\t=\t8200\t200\t*\t*",
        "bad\t147\tchr1\t8201\t60\t50M\t=\t8000\t-200\t*\t*",
    ]
    bam = _make_bam(tmp_path, records, samtools_path)
    row = _run_summary(
        tmp_path,
        samtools_path,
        bam=bam,
        windows=((1, 139), (140, 200), (201, 300)),
    )

    assert row["analysis_total_templates"] == "5"
    assert row["analysis_valid_templates"] == "1"
    assert row["analysis_invalid_templates"] == "4"
    assert row["analysis_tlen_zero_templates"] == "1"
    assert row["analysis_orphan_templates"] == "1"
    assert row["analysis_cross_contig_templates"] == "1"
    assert row["analysis_inconsistent_templates"] == "1"
    assert row["mono_fragments"] == "1"
    assert row["analysis_classified_templates"] == "1"


def test_missing_analysis_bam_reports_parse_failed_not_zero(tmp_path, samtools_path):
    bam = _make_bam(tmp_path, _pair("valid", 1000, 150), samtools_path)
    row = _run_summary(
        tmp_path,
        samtools_path,
        bam=bam,
        windows=((1, 139), (140, 200), (201, 300)),
        analysis_bam=tmp_path / "missing.bam",
    )

    assert row["analysis_status"] == "parse_failed"
    assert row["analysis_valid_templates"] == "NA"
    assert row["analysis_invalid_templates"] == "NA"
    assert row["analysis_fraction_classified"] == "NA"
    assert row["mono_fragments"] == "NA"
