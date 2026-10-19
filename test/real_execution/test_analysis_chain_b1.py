"""Real-tool contracts for the B1 updated analysis chain.

These tests pin the exact updated-profile commands, independent of Snakemake:
- PE mate-complete filter chain
- Picard mark-only duplicate handling
- analysis view duplicate exclusion
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest
from _tool_resolver import require_external_tools, resolve_tool

pytestmark = pytest.mark.real_execution

HEADER = [
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:10000",
    "@RG\tID:rg1\tSM:S1\tLB:lib1\tPL:ILLUMINA\tPU:unit1",
]
SEQ = "A" * 50
QUAL = "I" * 50


def _tool(name, env_var):
    candidate = resolve_tool(name, env_var)
    resolved = shutil.which(candidate)
    if not resolved and os.path.isabs(candidate):
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            resolved = candidate
    if not resolved:
        require_external_tools([name], "B1 analysis-chain real execution")
    return resolved


@pytest.fixture(scope="module")
def samtools_path():
    return _tool("samtools", "SAMTOOLS")


@pytest.fixture(scope="module")
def picard_path():
    return _tool("picard", "PICARD")


def _env(*tools):
    env = os.environ.copy()
    path_parts = [str(Path(tool).parent) for tool in tools]
    env["PATH"] = os.pathsep.join(path_parts + [env.get("PATH", "")])
    return env


def _run(command, *, env=None):
    result = subprocess.run(command, capture_output=True, text=True, env=env)
    assert result.returncode == 0, (
        f"Command failed: {command!r}\n"
        f"stdout:\n{result.stdout[-1000:]}\n"
        f"stderr:\n{result.stderr[-1000:]}"
    )
    return result


def _record(name, flag, pos, mapq, mate_pos=0, tlen=0, rnext=None):
    if rnext is None:
        rnext = "=" if mate_pos != 0 and tlen != 0 else "*"
    return "\t".join(
        [
            name,
            str(flag),
            "chr1",
            str(pos),
            str(mapq),
            "50M",
            rnext,
            str(mate_pos),
            str(tlen),
            SEQ,
            QUAL,
            "RG:Z:rg1",
        ]
    )


def _make_bam(tmp_path, lines, samtools_path):
    sam = tmp_path / "input.sam"
    raw = tmp_path / "input.raw.bam"
    bam = tmp_path / "input.bam"
    sam.write_text("\n".join(HEADER + lines) + "\n", encoding="utf-8")
    _run([samtools_path, "view", "-b", "-o", str(raw), str(sam)])
    _run([samtools_path, "sort", "-o", str(bam), str(raw)])
    _run([samtools_path, "index", str(bam)])
    return bam


def _updated_pe_filter(bam, tmp_path, samtools_path):
    filtered = tmp_path / "filtered.bam"
    names = tmp_path / "names.txt"
    out = tmp_path / "out.bam"
    _run(
        [
            samtools_path,
            "view",
            "-b",
            "-q",
            "30",
            "-F",
            "0xb0c",
            "-f",
            "2",
            "-o",
            str(filtered),
            str(bam),
        ]
    )
    records = _run([samtools_path, "view", str(filtered)]).stdout.splitlines()
    complete = {}
    for line in records:
        qname = line.split("\t")[0]
        complete[qname] = complete.get(qname, 0) + 1
    names.write_text(
        "\n".join(qname for qname, count in complete.items() if count == 2) + "\n",
        encoding="utf-8",
    )
    _run(
        [
            samtools_path,
            "view",
            "-b",
            "-N",
            str(names),
            "-o",
            str(out),
            str(filtered),
        ]
    )
    return out


def test_updated_pe_filter_keeps_only_complete_proper_pairs(tmp_path, samtools_path):
    lines = [
        # complete proper pair
        _record("proper", 99, 100, 60, mate_pos=200, tlen=150),
        _record("proper", 147, 200, 60, mate_pos=100, tlen=-150),
        # singleton / mate-unmapped
        _record("singleton", 73, 400, 60),
        # discordant
        _record("discordant", 97, 600, 60, mate_pos=900, tlen=350),
        _record("discordant", 145, 900, 60, mate_pos=600, tlen=-350),
        # QC-fail proper pair
        _record("qcfail", 611, 1100, 60, mate_pos=1200, tlen=150),
        _record("qcfail", 659, 1200, 60, mate_pos=1100, tlen=-150),
        # supplementary proper pair
        _record("supp", 2147, 1400, 60, mate_pos=1500, tlen=150),
        _record("supp", 2195, 1500, 60, mate_pos=1400, tlen=-150),
        # asymmetric MAPQ: only one mate survives MAPQ; QNAME completeness removes it
        _record("asym", 99, 1700, 60, mate_pos=1800, tlen=150),
        _record("asym", 147, 1800, 10, mate_pos=1700, tlen=-150),
    ]
    bam = _make_bam(tmp_path, lines, samtools_path)
    out = _updated_pe_filter(bam, tmp_path, samtools_path)
    records = _run([samtools_path, "view", str(out)]).stdout.strip().splitlines()
    assert len(records) == 2
    assert {line.split("\t")[0] for line in records} == {"proper"}
    assert {line.split("\t")[1] for line in records} == {"99", "147"}


def test_picard_mark_only_and_analysis_view(tmp_path, samtools_path, picard_path):
    lines = [
        _record("INST:1:FC:1:1101:100:200", 99, 100, 60, mate_pos=200, tlen=150),
        _record("INST:1:FC:1:1101:100:200", 147, 200, 60, mate_pos=100, tlen=-150),
        _record("INST:1:FC:1:1101:101:200", 99, 100, 60, mate_pos=200, tlen=150),
        _record("INST:1:FC:1:1101:101:200", 147, 200, 60, mate_pos=100, tlen=-150),
    ]
    bam = _make_bam(tmp_path, lines, samtools_path)
    marked = tmp_path / "marked.bam"
    metrics = tmp_path / "metrics.txt"
    env = _env(picard_path, samtools_path)
    _run(
        [
            picard_path,
            "MarkDuplicates",
            f"I={bam}",
            f"O={marked}",
            f"M={metrics}",
            "REMOVE_DUPLICATES=false",
        ],
        env=env,
    )
    flagged = _run(
        [samtools_path, "view", "-c", "-f", "0x400", str(marked)]
    ).stdout.strip()
    assert int(flagged) == 2
    unique = _run(
        [samtools_path, "view", "-c", "-F", "0x400", str(marked)]
    ).stdout.strip()
    assert int(unique) == 2
