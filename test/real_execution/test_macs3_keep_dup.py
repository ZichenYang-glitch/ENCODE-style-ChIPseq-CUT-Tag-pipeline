"""Prove MACS3 --keep-dup semantics for the SCI-2 BAM-layer duplicate contract.

The workflow decides duplicate policy at the BAM layer and renders
--keep-dup all on every MACS3 call. This test proves against the real
macs3 executable that:

- with --keep-dup all, two same-coordinate templates from two different
  libraries (RG/LB) both contribute to the pileup;
- with --keep-dup 1 (the previous implicit default behavior), one of them
  is collapsed away.

Usage:
    MACS3=/path/to/macs3 SAMTOOLS=/path/to/samtools \
        python3 -m pytest -m real_execution \
        test/real_execution/test_macs3_keep_dup.py -v
"""

import os
import shutil
import subprocess

import pytest
from _tool_resolver import require_external_tools, resolve_tool

pytestmark = pytest.mark.real_execution

CHROM_LEN = 1_000_000
GENOME_SIZE = "10000000"
READ_LEN = 50
FRAGMENT_LEN = 200
DUP_POS = 100_000
N_BACKGROUND = 40


def _resolve(name, env_var):
    candidate = resolve_tool(name, env_var)
    resolved = shutil.which(candidate)
    if not resolved and os.path.isabs(candidate):
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            resolved = candidate
    return resolved


@pytest.fixture(scope="module")
def macs3_path():
    resolved = _resolve("macs3", "MACS3")
    if not resolved:
        require_external_tools(["macs3"], "MACS3 keep-dup real execution")
    return resolved


@pytest.fixture(scope="module")
def samtools_path():
    resolved = _resolve("samtools", "SAMTOOLS")
    if not resolved:
        require_external_tools(["samtools"], "MACS3 keep-dup real execution")
    return resolved


def _write_dual_library_bam(tmp_path, samtools_path):
    """BAM with two libraries and one same-coordinate template pair each.

    fragDupA (LB1) and fragDupB (LB2) are proper PE pairs at identical
    coordinates — the same molecule position observed in two libraries.
    """
    lines = [
        "@HD\tVN:1.6\tSO:coordinate",
        f"@SQ\tSN:chr1\tLN:{CHROM_LEN}",
        "@RG\tID:rgA\tSM:S1\tLB:LB1",
        "@RG\tID:rgB\tSM:S1\tLB:LB2",
    ]
    seq = "A" * READ_LEN
    qual = "I" * READ_LEN

    def pair(name, rg, pos1):
        pos2 = pos1 + FRAGMENT_LEN - READ_LEN
        return [
            f"{name}\t99\tchr1\t{pos1}\t60\t{READ_LEN}M\t=\t{pos2}\t{FRAGMENT_LEN}\t{seq}\t{qual}\tRG:Z:{rg}",
            f"{name}\t147\tchr1\t{pos2}\t60\t{READ_LEN}M\t=\t{pos1}\t{-FRAGMENT_LEN}\t{seq}\t{qual}\tRG:Z:{rg}",
        ]

    lines += pair("fragDupA", "LB1", DUP_POS)
    lines += pair("fragDupB", "LB2", DUP_POS)
    for i in range(N_BACKGROUND):
        pos1 = 300_000 + i * 1000
        lines += pair(f"fragBg{i}", "LB1", pos1)

    sam_path = tmp_path / "dual.sam"
    sam_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    bam_path = tmp_path / "dual.sorted.bam"
    view = subprocess.run(
        [samtools_path, "view", "-b", "-o", str(bam_path), str(sam_path)],
        capture_output=True,
        text=True,
    )
    assert view.returncode == 0, f"samtools view failed: {view.stderr[-500:]}"
    return bam_path


def _run_macs3(macs3_path, bam_path, outdir, keep_dup, flag="--keep-dup"):
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        macs3_path,
        "callpeak",
        "-t",
        str(bam_path),
        "-f",
        "BAMPE",
        "-g",
        GENOME_SIZE,
        "-q",
        "0.05",
        "-B",
        "--nomodel",
        flag,
        keep_dup,
        "-n",
        "test",
        "--outdir",
        str(outdir),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, (
        f"macs3 callpeak failed: {cmd!r}\nstderr:\n{result.stderr[-1000:]}"
    )
    return outdir


def _pileup_at(outdir, position):
    path = outdir / "test_treat_pileup.bdg"
    assert path.exists(), f"missing {path}"
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        if len(fields) == 4 and int(fields[1]) <= position < int(fields[2]):
            return float(fields[3])
    return 0.0


def test_keep_dup_all_preserves_same_coordinate_templates_across_libraries(
    tmp_path, macs3_path, samtools_path
):
    bam = _write_dual_library_bam(tmp_path, samtools_path)
    center = DUP_POS + FRAGMENT_LEN // 2

    all_dir = _run_macs3(macs3_path, bam, tmp_path / "all", "all")
    one_dir = _run_macs3(macs3_path, bam, tmp_path / "one", "1")

    assert _pileup_at(all_dir, center) == 2.0
    assert _pileup_at(one_dir, center) == 1.0


def test_abbreviated_keep_d_is_the_same_flag_for_real_macs3(
    tmp_path, macs3_path, samtools_path
):
    """argparse prefix matching makes --keep-d identical to --keep-dup.

    The pinned macs3 collapses the cross-library duplicate pair under
    --keep-d 1 exactly like --keep-dup 1. This is why the validator must
    recognize and reject abbreviated spellings that conflict with the
    canonical --keep-dup all instead of only matching the full name.
    """
    bam = _write_dual_library_bam(tmp_path, samtools_path)
    center = DUP_POS + FRAGMENT_LEN // 2
    abbrev_dir = _run_macs3(macs3_path, bam, tmp_path / "abbrev", "1", flag="--keep-d")
    assert _pileup_at(abbrev_dir, center) == 1.0
