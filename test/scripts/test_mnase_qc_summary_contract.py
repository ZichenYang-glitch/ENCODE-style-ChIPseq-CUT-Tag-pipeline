"""Fast-tier contracts for the versioned MNase QC summary writer.

The real samtools path is covered by the real-execution tests.  These tests
use a controlled fake ``samtools`` so the v2 window classification, status,
stage counts, and read-count separation run in the deterministic fast tier.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import stat
import subprocess
import sys

from encode_pipeline.qc.contract import MNASE_QC_V2_COLUMNS

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "mnase_qc_summary.py"
HEADER = [
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:100000",
    "@SQ\tSN:chr2\tLN:100000",
]


def _pair(name: str, start: int, length: int) -> list[str]:
    return [
        f"{name}\t99\tchr1\t{start}\t60\t50M\t=\t{start + length}\t{length}\t*\t*",
        f"{name}\t147\tchr1\t{start + length}\t60\t50M\t=\t{start}\t-{length}\t*\t*",
    ]


INVALID = (
    _pair("z1", 1500, 0)
    + ["o1\t99\tchr1\t1800\t60\t50M\t=\t1950\t150\t*\t*"]
    + [
        "c1\t99\tchr1\t2100\t60\t50M\tchr2\t2300\t200\t*\t*",
        "c1\t147\tchr2\t2300\t60\t50M\tchr1\t2100\t-200\t*\t*",
    ]
    + _pair("i1", 2500, 200)
)
INVALID[-1] = "i1\t147\tchr1\t2700\t60\t50M\t=\t2500\t-250\t*\t*"


def _fake_samtools(tmp_path: Path, count: int = 12) -> Path:
    path = tmp_path / "samtools"
    path.write_text(
        "#!/bin/sh\n"
        'if [ "$1" = "view" ] && [ "$2" = "-c" ]; then\n'
        f"  echo {count}\n"
        "  exit 0\n"
        "fi\n"
        'if [ "$1" = "view" ]; then\n'
        '  cat "$SAM_SOURCE"\n'
        "  exit 0\n"
        "fi\n"
        "exit 2\n",
        encoding="utf-8",
    )
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


def _write_bam_placeholder(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"placeholder; fake samtools ignores it\n")
    return path


def _run_script(
    tmp_path: Path,
    *,
    records,
    windows=((1, 139), (140, 200), (201, 300)),
    analysis_bam: Path | None = None,
    sub_bam: Path | None = None,
    mono_bam: Path | None = None,
    di_bam: Path | None = None,
    stage_metrics: Path | None = None,
    count: int = 12,
) -> subprocess.CompletedProcess:
    sam = tmp_path / "input.sam"
    sam.write_text("\n".join(HEADER + list(records)) + "\n", encoding="utf-8")
    if analysis_bam is None:
        analysis_bam = _write_bam_placeholder(tmp_path / "analysis.bam")
    output = tmp_path / "mnase.tsv"
    argv = [
        sys.executable,
        str(SCRIPT),
        "--sample",
        "S1",
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
        "--analysis-bam",
        str(analysis_bam),
        "--output",
        str(output),
    ]
    for flag, path in (
        ("--sub-bam", sub_bam),
        ("--mono-bam", mono_bam),
        ("--di-bam", di_bam),
        ("--stage-metrics", stage_metrics),
    ):
        if path is not None:
            argv[argv.index("--output") : argv.index("--output")] = [
                flag,
                str(path),
            ]
    env = os.environ.copy()
    fake = _fake_samtools(tmp_path, count=count)
    env["SAMTOOLS"] = str(fake)
    env["SAM_SOURCE"] = str(sam)
    env["PATH"] = f"{tmp_path}{os.pathsep}{env.get('PATH', '')}"
    return subprocess.run(argv, capture_output=True, text=True, env=env)


def _read_row(path: Path) -> dict:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert list(reader.fieldnames or []) == list(MNASE_QC_V2_COLUMNS)
        rows = list(reader)
    assert len(rows) == 1
    return rows[0]


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


def test_non_overlapping_windows_classify_each_valid_template_once(tmp_path):
    stage = _stage_metrics(tmp_path / "stage.json")
    bam = _write_bam_placeholder(tmp_path / "sub.bam")
    records = _pair("a", 100, 100) + _pair("b", 300, 150)
    records += _pair("c", 600, 250) + _pair("d", 1000, 400)

    result = _run_script(
        tmp_path,
        records=records,
        sub_bam=bam,
        mono_bam=bam,
        di_bam=bam,
        stage_metrics=stage,
    )

    assert result.returncode == 0, result.stderr
    row = _read_row(tmp_path / "mnase.tsv")
    assert row["schema_version"] == "2"
    assert row["analysis_status"] == "ok"
    assert row["analysis_total_templates"] == "4"
    assert row["analysis_valid_templates"] == "4"
    assert row["analysis_invalid_templates"] == "0"
    assert row["sub_fragments"] == "1"
    assert row["mono_fragments"] == "1"
    assert row["di_fragments"] == "1"
    assert row["analysis_classified_templates"] == "3"
    assert row["analysis_unclassified_templates"] == "1"
    assert row["analysis_fraction_classified"] == "0.750000"
    assert row["window_overlap"] == "false"
    assert row["sub_reads"] == "12"
    assert row["mono_reads"] == "12"
    assert row["di_reads"] == "12"
    assert row["upstream_dropped_templates"] == "2"
    assert row["upstream_dropped_alignments"] == "4"


def test_overlapping_windows_use_the_union_for_classification(tmp_path):
    records = (
        _pair("a", 100, 139)
        + _pair("b", 300, 200)
        + _pair("c", 600, 300)
        + _pair("d", 1000, 400)
        + _pair("e", 1500, 450)
    )

    result = _run_script(
        tmp_path,
        records=records,
        windows=((1, 200), (140, 300), (250, 400)),
    )

    assert result.returncode == 0, result.stderr
    row = _read_row(tmp_path / "mnase.tsv")
    assert row["window_overlap"] == "true"
    assert row["sub_fragments"] == "2"
    assert row["mono_fragments"] == "2"
    assert row["di_fragments"] == "2"
    assert row["analysis_classified_templates"] == "4"
    assert row["analysis_unclassified_templates"] == "1"
    assert row["analysis_fraction_classified"] == "0.800000"


def test_zero_valid_templates_is_a_conserved_zero(tmp_path):
    result = _run_script(tmp_path, records=INVALID)

    assert result.returncode == 0, result.stderr
    row = _read_row(tmp_path / "mnase.tsv")
    assert row["analysis_status"] == "no_valid_templates"
    assert row["analysis_total_templates"] == "4"
    assert row["analysis_valid_templates"] == "0"
    assert row["analysis_invalid_templates"] == "4"
    assert row["analysis_classified_templates"] == "0"
    assert row["analysis_unclassified_templates"] == "0"
    assert row["analysis_fraction_classified"] == "NA"
    assert row["sub_fragments"] == "0"
    assert row["mono_fragments"] == "0"
    assert row["di_fragments"] == "0"


def test_missing_analysis_bam_is_parse_failed_with_na_counts(tmp_path):
    result = _run_script(
        tmp_path,
        records=_pair("a", 100, 150),
        analysis_bam=tmp_path / "missing.bam",
        sub_bam=tmp_path / "missing-sub.bam",
        mono_bam=tmp_path / "missing-mono.bam",
        di_bam=tmp_path / "missing-di.bam",
    )

    assert result.returncode == 0, result.stderr
    row = _read_row(tmp_path / "mnase.tsv")
    assert row["analysis_status"] == "parse_failed"
    for column in (
        "analysis_total_templates",
        "analysis_valid_templates",
        "analysis_invalid_templates",
        "analysis_classified_templates",
        "analysis_unclassified_templates",
        "analysis_fraction_classified",
        "sub_fragments",
        "mono_fragments",
        "di_fragments",
    ):
        assert row[column] == "NA"
    assert row["sub_reads"] == "NA"
    assert row["mono_reads"] == "NA"
    assert row["di_reads"] == "NA"
