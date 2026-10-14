"""Pin MACS3 PE/SE flag semantics that the workflow's layout split relies on.

The ATAC/CUT&Tag MACS3 argument builders pass --nomodel --shift -100
--extsize 200 to SE (-f BAM) samples only. This test proves, against the
real macs3 executable, that:

- PE (-f BAMPE): adding --nomodel --shift -100 --extsize 200 changes
  nothing — _peaks.xls data rows and _treat_pileup.bdg are byte-identical,
  so removing them from the workflow is a zero-behavior-change cleanup.
- SE (-f BAM): adding --shift -100 (with fixed --extsize) changes the
  outputs, so the flags are meaningful exactly where the workflow keeps
  them.
"""

import os
import shutil
import subprocess

import pytest
from _tool_resolver import require_external_tools, resolve_tool

pytestmark = pytest.mark.real_execution

CHROM_LEN = 1_000_000
GENOME_SIZE = "10000000"
N_FRAGMENTS = 60
READ_LEN = 50
FRAGMENT_LEN = 200
CLUSTER_START = 100_000


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
        require_external_tools(["macs3"], "MACS3 PE/SE flag-semantics real execution")
    return resolved


@pytest.fixture(scope="module")
def samtools_path():
    resolved = _resolve("samtools", "SAMTOOLS")
    if not resolved:
        require_external_tools(
            ["samtools"], "MACS3 PE/SE flag-semantics real execution"
        )
    return resolved


def _write_pe_bam(tmp_path, samtools_path):
    """Build a coordinate-sorted BAM with N_FRAGMENTS proper PE pairs."""
    lines = [
        "@HD\tVN:1.6\tSO:coordinate",
        f"@SQ\tSN:chr1\tLN:{CHROM_LEN}",
    ]
    seq = "A" * READ_LEN
    qual = "I" * READ_LEN
    for i in range(N_FRAGMENTS):
        pos1 = CLUSTER_START + i * 10
        pos2 = pos1 + FRAGMENT_LEN - READ_LEN
        tlen = FRAGMENT_LEN
        lines.append(
            f"frag{i}\t99\tchr1\t{pos1}\t60\t{READ_LEN}M\t=\t{pos2}\t{tlen}\t{seq}\t{qual}"
        )
        lines.append(
            f"frag{i}\t147\tchr1\t{pos2}\t60\t{READ_LEN}M\t=\t{pos1}\t{-tlen}\t{seq}\t{qual}"
        )
    sam_path = tmp_path / "mini.sam"
    sam_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    bam_path = tmp_path / "mini.bam"
    view = subprocess.run(
        [samtools_path, "view", "-b", "-o", str(bam_path), str(sam_path)],
        capture_output=True,
        text=True,
    )
    assert view.returncode == 0, f"samtools view failed: {view.stderr[-500:]}"
    sorted_path = tmp_path / "mini.sorted.bam"
    sort = subprocess.run(
        [samtools_path, "sort", "-o", str(sorted_path), str(bam_path)],
        capture_output=True,
        text=True,
    )
    assert sort.returncode == 0, f"samtools sort failed: {sort.stderr[-500:]}"
    return sorted_path


def _run_macs3(macs3_path, bam_path, outdir, extra_flags):
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        macs3_path,
        "callpeak",
        "-t",
        str(bam_path),
        "-g",
        GENOME_SIZE,
        "-q",
        "0.05",
        "-B",
        "--nomodel",
        # Same -n for every run: _peaks.xls embeds the command line, outdir,
        # and peak names, so a differing name would pollute the comparison.
        "-n",
        "test",
        "--outdir",
        str(outdir),
        *extra_flags,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, (
        f"macs3 callpeak failed: {cmd!r}\nstderr:\n{result.stderr[-1000:]}"
    )
    return outdir


def _peaks_xls_data(outdir):
    """Data rows only: '#' comment lines embed argv/outdir and differ legitimately."""
    path = outdir / "test_peaks.xls"
    assert path.exists(), f"missing {path}"
    return b"".join(
        line
        for line in path.read_bytes().splitlines(keepends=True)
        if not line.startswith(b"#")
    )


def _treat_pileup(outdir):
    path = outdir / "test_treat_pileup.bdg"
    assert path.exists(), f"missing {path}"
    return path.read_bytes()


def test_pe_flags_are_byte_identical_noop(tmp_path, macs3_path, samtools_path):
    """-f BAMPE: --nomodel --shift -100 --extsize 200 must change nothing."""
    bam = _write_pe_bam(tmp_path, samtools_path)
    base_flags = ["-f", "BAMPE"]
    without = _run_macs3(macs3_path, bam, tmp_path / "pe_without", base_flags)
    with_flags = _run_macs3(
        macs3_path,
        bam,
        tmp_path / "pe_with",
        [*base_flags, "--shift", "-100", "--extsize", "200"],
    )
    assert _peaks_xls_data(without) == _peaks_xls_data(with_flags)
    assert _treat_pileup(without) == _treat_pileup(with_flags)


def test_se_shift_actually_changes_output(tmp_path, macs3_path, samtools_path):
    """-f BAM: adding --shift -100 at fixed --extsize must change the outputs."""
    bam = _write_pe_bam(tmp_path, samtools_path)
    base_flags = ["-f", "BAM", "--extsize", "200"]
    without = _run_macs3(macs3_path, bam, tmp_path / "se_without", base_flags)
    with_shift = _run_macs3(
        macs3_path,
        bam,
        tmp_path / "se_with",
        [*base_flags, "--shift", "-100"],
    )
    assert _treat_pileup(without) != _treat_pileup(with_shift)
