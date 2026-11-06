"""Real-tool contracts for template-aware blacklist filtering."""

import os
import shutil
import subprocess
from pathlib import Path

import pytest
from _tool_resolver import require_external_tools, resolve_tool

pytestmark = pytest.mark.real_execution

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "filter_blacklist_templates.py"

BM_HEADER = [
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:10000",
]
SEQ = "A" * 50
QUAL = "I" * 50


def _samtools():
    candidate = resolve_tool("samtools", "SAMTOOLS")
    resolved = shutil.which(candidate)
    if not resolved and os.path.isabs(candidate) and os.access(candidate, os.X_OK):
        resolved = candidate
    if not resolved:
        require_external_tools(["samtools"], "blacklist template filter")
    return resolved


@pytest.fixture(scope="module")
def samtools_path():
    return _samtools()


def _run(command, *, env=None):
    result = subprocess.run(command, capture_output=True, text=True, env=env)
    assert result.returncode == 0, (
        f"Command failed: {command!r}\n"
        f"stdout:\n{result.stdout[-1000:]}\n"
        f"stderr:\n{result.stderr[-1000:]}"
    )
    return result


def _pair(name, start, mate_start, tlen):
    return [
        "\t".join(
            [
                name,
                "99",
                "chr1",
                str(start),
                "60",
                "50M",
                "=",
                str(mate_start),
                str(tlen),
                SEQ,
                QUAL,
            ]
        ),
        "\t".join(
            [
                name,
                "147",
                "chr1",
                str(mate_start),
                "60",
                "50M",
                "=",
                str(start),
                str(-tlen),
                SEQ,
                QUAL,
            ]
        ),
    ]


def _se(name, start):
    return "\t".join(
        [
            name,
            "0",
            "chr1",
            str(start),
            "60",
            "50M",
            "*",
            "0",
            "0",
            SEQ,
            QUAL,
        ]
    )


def _make_bam(tmp_path, samtools_path, records):
    sam = tmp_path / "input.sam"
    raw = tmp_path / "raw.bam"
    bam = tmp_path / "input.bam"
    sam.write_text("\n".join(BM_HEADER + records) + "\n", encoding="utf-8")
    _run([samtools_path, "view", "-b", "-o", str(raw), str(sam)])
    _run([samtools_path, "sort", "-o", str(bam), str(raw)])
    _run([samtools_path, "index", str(bam)])
    return bam


def _run_filter(tmp_path, samtools_path, bam, blacklist, layout):
    out = tmp_path / "filtered.bam"
    dropped = tmp_path / "dropped.tsv"
    script = (
        'python3 "$1" --bam "$2" --blacklist "$3" --layout "$4" '
        '--output "$5" --dropped-count "$6" --sample S1'
    )
    env = os.environ.copy()
    env["PATH"] = os.pathsep.join(
        [str(Path(samtools_path).parent), env.get("PATH", "")]
    )
    _run(
        [
            "bash",
            "-c",
            script,
            "bash",
            str(SCRIPT),
            str(bam),
            str(blacklist),
            layout,
            str(out),
            str(dropped),
        ],
        env=env,
    )
    return out, dropped


def test_pe_blacklist_removes_whole_fragment(tmp_path, samtools_path):
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text(
        "chr1\t0\t10\nchr1\t100\t110\nchr1\t300\t310\n",
        encoding="utf-8",
    )
    records = []
    # mate1 in blacklist -> whole template removed
    records += _pair("mate1_hit", 5, 200, 245)
    # mate2 in blacklist -> whole template removed
    records += _pair("mate2_hit", 250, 305, 195)
    # blacklist only in the insert gap -> whole template removed
    records += _pair("gap_hit", 50, 130, 100)
    # boundary: fragment ends before blacklist starts -> retained
    records += _pair("boundary_ok", 150, 200, 100)
    # normal pair -> retained
    records += _pair("normal", 400, 500, 150)
    # orphan/non-proper -> removed
    records.append(
        "\t".join(
            [
                "orphan",
                "73",
                "chr1",
                "600",
                "60",
                "50M",
                "*",
                "0",
                "0",
                SEQ,
                QUAL,
            ]
        )
    )
    bam = _make_bam(tmp_path, samtools_path, records)
    out, dropped = _run_filter(tmp_path, samtools_path, bam, blacklist, "PE")

    kept = _run([samtools_path, "view", str(out)]).stdout.splitlines()
    names = {line.split("\t")[0] for line in kept}
    assert names == {"boundary_ok", "normal"}
    assert len(kept) == 4

    fields = dropped.read_text(encoding="utf-8").strip().splitlines()
    assert fields[0].startswith("sample\tlayout\tdropped_templates")
    # Under R9, blacklist table only counts templates dropped by blacklist (3 hits; orphan is invalid)
    assert fields[1].split("\t")[:3] == ["S1", "PE", "3"]


def test_pe_incomplete_mate_is_removed(tmp_path, samtools_path):
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr1\t100\t110\n", encoding="utf-8")
    records = [
        _pair("incomplete", 400, 500, 150)[0],
        *_pair("normal", 700, 800, 150),
    ]
    bam = _make_bam(tmp_path, samtools_path, records)
    out, dropped = _run_filter(tmp_path, samtools_path, bam, blacklist, "PE")

    kept = _run([samtools_path, "view", str(out)]).stdout.splitlines()
    assert {line.split("\t")[0] for line in kept} == {"normal"}
    fields = dropped.read_text(encoding="utf-8").strip().splitlines()
    # Incomplete mate did not hit blacklist, so blacklist dropped_templates is 0
    assert fields[1].split("\t")[:3] == ["S1", "PE", "0"]


def test_se_blacklist_removes_overlapping_reads(tmp_path, samtools_path):
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr1\t100\t110\n", encoding="utf-8")
    records = [
        _se("se_hit", 105),
        _se("se_boundary", 150),
        _se("se_normal", 300),
    ]
    bam = _make_bam(tmp_path, samtools_path, records)
    out, dropped = _run_filter(tmp_path, samtools_path, bam, blacklist, "SE")

    kept = _run([samtools_path, "view", str(out)]).stdout.splitlines()
    names = {line.split("\t")[0] for line in kept}
    assert names == {"se_boundary", "se_normal"}
    fields = dropped.read_text(encoding="utf-8").strip().splitlines()
    assert fields[1].split("\t")[:3] == ["S1", "SE", "1"]
