"""Fast-tier contracts for the versioned CUT&Tag fragment-size writer.

The real execution path is covered by the samtools-backed real-execution
tests.  These tests use a controlled fake ``samtools view`` so the v2
counting, status, and stage-row contracts run in the deterministic fast tier.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import stat
import subprocess
import sys

from encode_pipeline.qc.contract import FRAGMENT_QC_V2_COLUMNS

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "calc_cuttag_fragment_size.py"
HEADER = [
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:100000",
    "@SQ\tSN:chr2\tLN:100000",
]
VALID = [
    "v1\t99\tchr1\t100\t60\t50M\t=\t250\t150\t*\t*",
    "v1\t147\tchr1\t250\t60\t50M\t=\t100\t-150\t*\t*",
    "v2\t99\tchr1\t400\t60\t50M\t=\t600\t200\t*\t*",
    "v2\t147\tchr1\t600\t60\t50M\t=\t400\t-200\t*\t*",
]
INVALID = [
    "z1\t99\tchr1\t1500\t60\t50M\t=\t1600\t0\t*\t*",
    "z1\t147\tchr1\t1600\t60\t50M\t=\t1500\t0\t*\t*",
    "o1\t99\tchr1\t1800\t60\t50M\t=\t1950\t150\t*\t*",
    "c1\t99\tchr1\t2100\t60\t50M\tchr2\t2300\t200\t*\t*",
    "c1\t147\tchr2\t2300\t60\t50M\tchr1\t2100\t-200\t*\t*",
    "i1\t99\tchr1\t2500\t60\t50M\t=\t2700\t200\t*\t*",
    "i1\t147\tchr1\t2700\t60\t50M\t=\t2500\t-250\t*\t*",
]


def _fake_samtools(tmp_path: Path) -> Path:
    path = tmp_path / "samtools"
    path.write_text(
        "#!/bin/sh\n"
        'if [ "$1" = "view" ]; then\n'
        '  cat "$SAM_SOURCE"\n'
        "  exit 0\n"
        "fi\n"
        "exit 2\n",
        encoding="utf-8",
    )
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return path


def _run_script(
    tmp_path: Path,
    *,
    records,
    layout: str = "PE",
    bam: Path | None = None,
    stage_metrics: Path | None = None,
) -> subprocess.CompletedProcess:
    sam = tmp_path / "input.sam"
    sam.write_text("\n".join(HEADER + list(records)) + "\n", encoding="utf-8")
    analysis_bam = bam if bam is not None else tmp_path / "analysis.bam"
    if bam is None:
        analysis_bam.write_bytes(b"placeholder; fake samtools ignores it\n")
    output = tmp_path / "fragment.tsv"
    argv = [
        sys.executable,
        str(SCRIPT),
        "--sample",
        "S1",
        "--bam",
        str(analysis_bam),
        "--layout",
        layout,
        "--duplicate-analysis",
        "exclude",
        "--blacklist-applied",
        "false",
        "--output",
        str(output),
    ]
    if stage_metrics is not None:
        argv[argv.index("--output") : argv.index("--output")] = [
            "--stage-metrics",
            str(stage_metrics),
        ]
    env = os.environ.copy()
    env["SAMTOOLS"] = str(_fake_samtools(tmp_path))
    env["SAM_SOURCE"] = str(sam)
    return subprocess.run(argv, capture_output=True, text=True, env=env)


def _read_row(path: Path) -> dict:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert list(reader.fieldnames or []) == list(FRAGMENT_QC_V2_COLUMNS)
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


def test_pe_counts_valid_and_invalid_templates_with_stage_metrics(tmp_path):
    stage = _stage_metrics(tmp_path / "stage.json")
    result = _run_script(
        tmp_path,
        records=VALID + INVALID,
        stage_metrics=stage,
    )

    assert result.returncode == 0, result.stderr
    row = _read_row(tmp_path / "fragment.tsv")
    assert row["schema_version"] == "2"
    assert row["sample"] == "S1"
    assert row["layout"] == "PE"
    assert row["status"] == "ok"
    assert row["fragment_count"] == "2"
    assert row["fragment_mean"] == "175.0"
    assert row["fragment_median"] == "175.0"
    assert row["fragment_mode"] == "150"
    assert row["fragment_min"] == "150"
    assert row["fragment_max"] == "200"
    assert row["fraction_150_300"] == "1.0000"
    assert row["fraction_lt_150"] == "0.0000"
    assert row["templates_total"] == "6"
    assert row["templates_valid"] == "2"
    assert row["templates_invalid"] == "4"
    assert row["invalid_tlen_zero"] == "1"
    assert row["invalid_orphan"] == "1"
    assert row["invalid_cross_contig"] == "1"
    assert row["invalid_inconsistent"] == "1"
    assert row["alignments_total"] == "11"
    assert row["upstream_dropped_templates"] == "2"
    assert row["upstream_dropped_alignments"] == "4"
    assert row["upstream_alignments_total"] == "6"
    assert row["upstream_alignments_kept"] == "2"


def test_pe_zero_valid_templates_is_not_a_biological_zero(tmp_path):
    result = _run_script(tmp_path, records=INVALID)

    assert result.returncode == 0, result.stderr
    row = _read_row(tmp_path / "fragment.tsv")
    assert row["status"] == "no_fragments"
    assert row["fragment_count"] == "0"
    for column in (
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
    ):
        assert row[column] == "NA"
    assert row["templates_total"] == "4"
    assert row["templates_valid"] == "0"
    assert row["templates_invalid"] == "4"
    assert row["upstream_dropped_templates"] == "NA"


def test_se_is_explicitly_not_supported_with_na_fragment_fields(tmp_path):
    result = _run_script(
        tmp_path,
        records=[],
        layout="SE",
        bam=tmp_path / "missing.bam",
    )

    assert result.returncode == 0, result.stderr
    row = _read_row(tmp_path / "fragment.tsv")
    assert row["layout"] == "SE"
    assert row["status"] == "layout_not_supported"
    assert row["fragment_count"] == "NA"
    assert row["fragment_mean"] == "NA"
    assert row["input_stage"] == "analysis_bam"
    assert row["duplicate_analysis"] == "exclude"
    assert row["blacklist_applied"] == "false"


def test_missing_pe_analysis_bam_fails_closed(tmp_path):
    result = _run_script(
        tmp_path,
        records=VALID,
        bam=tmp_path / "missing.bam",
    )

    assert result.returncode == 1
    assert "analysis BAM not found" in result.stderr
    assert not (tmp_path / "fragment.tsv").exists()
