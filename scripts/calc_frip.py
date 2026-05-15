#!/usr/bin/env python3
"""Calculate FRiP (Fraction of Reads in Peaks) for a single sample.

Stage 3a: read-record based, not fragment-based.
Uses samtools and bedtools via subprocess pipes (no shell quoting).

Usage:
    python3 scripts/calc_frip.py \\
        --sample SAMPLE \\
        --bam sample.final.bam \\
        --peaks sample_peaks.narrowPeak \\
        --output sample.frip.tsv

Output is a tab-delimited TSV with columns:
    sample  total_reads  reads_in_peaks  frip  bam  peaks
"""

import argparse
import os
import subprocess
import sys


def _count_bam(bam_path: str) -> int:
    """Return total mapped read count via samtools view -c."""
    result = subprocess.run(
        ["samtools", "view", "-c", bam_path],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(
            f"ERROR: samtools view -c failed on {bam_path}: {result.stderr.strip()}",
            file=sys.stderr,
        )
        sys.exit(1)
    return int(result.stdout.strip())


def _count_reads_in_peaks(bam_path: str, peaks_path: str) -> int:
    """Return number of reads in a BAM that overlap a peak file.

    Uses a two-subprocess pipe to avoid shell quoting issues:
        bedtools intersect -u -abam BAM -b PEAKS | samtools view -c
    """
    bedtools = subprocess.Popen(
        ["bedtools", "intersect", "-u", "-abam", bam_path, "-b", peaks_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    samtools = subprocess.run(
        ["samtools", "view", "-c"],
        stdin=bedtools.stdout,
        capture_output=True,
        text=True,
    )
    bedtools.stdout.close()
    bedtools_rc = bedtools.wait()

    if bedtools_rc != 0:
        bedtools_err = bedtools.stderr.read().decode("utf-8", errors="replace")
        print(
            f"ERROR: bedtools intersect failed on {bam_path}: {bedtools_err.strip()}",
            file=sys.stderr,
        )
        sys.exit(1)
    if samtools.returncode != 0:
        print(
            f"ERROR: samtools view -c (pipe) failed: {samtools.stderr.strip()}",
            file=sys.stderr,
        )
        sys.exit(1)
    return int(samtools.stdout.strip())


def main():
    parser = argparse.ArgumentParser(
        description="Calculate FRiP (Fraction of Reads in Peaks)"
    )
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--bam", required=True, help="Path to BAM file")
    parser.add_argument("--peaks", required=True, help="Path to peak file")
    parser.add_argument("--output", required=True, help="Path to output TSV")
    args = parser.parse_args()

    if not os.path.isfile(args.bam):
        print(f"ERROR: BAM file not found: {args.bam}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(args.peaks):
        print(f"ERROR: Peak file not found: {args.peaks}", file=sys.stderr)
        sys.exit(1)

    total = _count_bam(args.bam)

    if total == 0:
        frip = "NA"
        reads_in_peaks = 0
    else:
        reads_in_peaks = _count_reads_in_peaks(args.bam, args.peaks)
        frip = f"{reads_in_peaks / total:.6f}"

    with open(args.output, "w") as fh:
        fh.write("sample\ttotal_reads\treads_in_peaks\tfrip\tbam\tpeaks\n")
        fh.write(
            f"{args.sample}\t{total}\t{reads_in_peaks}\t"
            f"{frip}\t{args.bam}\t{args.peaks}\n"
        )


if __name__ == "__main__":
    main()
