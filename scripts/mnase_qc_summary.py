#!/usr/bin/env python3
"""MNase-seq per-sample QC summary — stdlib-only.

Produces a single-row TSV with fragment stratification metadata,
read counts, and config-derived parameters.  Read counts are obtained
via ``samtools view -c``; missing/broken BAMs produce ``NA``.
"""

from __future__ import annotations

import argparse
import csv
import os
import subprocess


def _count_reads(bam_path: str) -> str:
    """Return read count (str) from ``samtools view -c``, or ``NA``."""
    if not bam_path or not os.path.isfile(bam_path):
        return "NA"
    try:
        cp = subprocess.run(
            ["samtools", "view", "-c", bam_path],
            capture_output=True,
            text=True,
            timeout=300,
        )
        if cp.returncode == 0:
            return cp.stdout.strip()
    except (OSError, subprocess.TimeoutExpired):
        pass
    return "NA"


def _file_status(path: str) -> str:
    """Return ``present`` or ``missing`` for a path."""
    if path and os.path.isfile(path):
        return "present"
    return "missing"


def main() -> None:
    p = argparse.ArgumentParser(
        description="MNase-seq per-sample QC summary (stdlib-only)."
    )
    p.add_argument("--sample", required=True)
    p.add_argument("--assay", required=True)
    p.add_argument("--peak-mode", required=True)
    p.add_argument("--sub-min", type=int, required=True)
    p.add_argument("--sub-max", type=int, required=True)
    p.add_argument("--mono-min", type=int, required=True)
    p.add_argument("--mono-max", type=int, required=True)
    p.add_argument("--di-min", type=int, required=True)
    p.add_argument("--di-max", type=int, required=True)
    p.add_argument("--dyad-min", type=int, required=True)
    p.add_argument("--dyad-max", type=int, required=True)
    p.add_argument("--sub-bam", default="")
    p.add_argument("--mono-bam", default="")
    p.add_argument("--di-bam", default="")
    p.add_argument("--dyad-bw", default="")
    p.add_argument("--mono-bw", default="")
    p.add_argument("--insert-size-metrics", default="")
    p.add_argument("--danpos3-enabled", default="false")
    p.add_argument("--inps-enabled", default="false")
    p.add_argument("--sem-enabled", default="false")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    header = [
        "sample",
        "assay",
        "peak_mode",
        "sub_min",
        "sub_max",
        "mono_min",
        "mono_max",
        "di_min",
        "di_max",
        "dyad_min",
        "dyad_max",
        "sub_bam",
        "mono_bam",
        "di_bam",
        "dyad_bigwig",
        "mono_bigwig",
        "insert_size_metrics",
        "caller_danpos3_enabled",
        "caller_inps_enabled",
        "caller_sem_enabled",
        "sub_reads",
        "mono_reads",
        "di_reads",
    ]

    row = {
        "sample": args.sample,
        "assay": args.assay,
        "peak_mode": args.peak_mode,
        "sub_min": str(args.sub_min),
        "sub_max": str(args.sub_max),
        "mono_min": str(args.mono_min),
        "mono_max": str(args.mono_max),
        "di_min": str(args.di_min),
        "di_max": str(args.di_max),
        "dyad_min": str(args.dyad_min),
        "dyad_max": str(args.dyad_max),
        "sub_bam": args.sub_bam,
        "mono_bam": args.mono_bam,
        "di_bam": args.di_bam,
        "dyad_bigwig": args.dyad_bw,
        "mono_bigwig": args.mono_bw,
        "insert_size_metrics": (
            args.insert_size_metrics
            if os.path.isfile(args.insert_size_metrics)
            else "NA"
        ),
        "caller_danpos3_enabled": args.danpos3_enabled,
        "caller_inps_enabled": args.inps_enabled,
        "caller_sem_enabled": args.sem_enabled,
        "sub_reads": _count_reads(args.sub_bam),
        "mono_reads": _count_reads(args.mono_bam),
        "di_reads": _count_reads(args.di_bam),
    }

    outdir = os.path.dirname(args.output)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerow(row)


if __name__ == "__main__":
    main()
