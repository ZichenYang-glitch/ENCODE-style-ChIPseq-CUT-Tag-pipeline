#!/usr/bin/env python3
"""MNase-seq per-sample QC summary — stdlib-only.

Produces one versioned row with fragment stratification metadata, read
counts, and config-derived parameters.  Read counts come from the
sub/mono/di BAMs via ``samtools view -c`` and remain separate from fragment
counts derived from the canonical analysis BAM.  Fragment windows may
overlap; unclassified templates are computed against the union of windows,
not as ``total - sum(window_counts)``.

Missing or malformed analysis input is reported as ``analysis_status``
``parse_failed`` with NA counts; it is never reported as a biological zero.
"""

from __future__ import annotations

import argparse
import os
import subprocess

from qc_common import (
    MNASE_QC_V2_COLUMNS,
    NA,
    format_fraction,
    load_stage_metrics,
    write_exact_tsv,
)
from qc_template import TemplateReader, classify_pe_template


def _count_reads(bam_path: str) -> str:
    """Return read count (str) from ``samtools view -c``, or ``NA``."""
    if not bam_path or not os.path.isfile(bam_path):
        return NA
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
    return NA


def _windows(args):
    return (
        ("sub", args.sub_min, args.sub_max),
        ("mono", args.mono_min, args.mono_max),
        ("di", args.di_min, args.di_max),
    )


def _windows_overlap(windows) -> bool:
    ordered = sorted(((low, high) for _, low, high in windows))
    for (_, high), (next_low, _) in zip(ordered, ordered[1:]):
        if next_low <= high:
            return True
    return False


def _classify_analysis_bam(path, windows):
    counts = {
        "analysis_total_templates": 0,
        "analysis_valid_templates": 0,
        "analysis_invalid_templates": 0,
        "analysis_tlen_zero_templates": 0,
        "analysis_orphan_templates": 0,
        "analysis_cross_contig_templates": 0,
        "analysis_inconsistent_templates": 0,
        "analysis_classified_templates": 0,
        "analysis_unclassified_templates": 0,
        "sub_fragments": 0,
        "mono_fragments": 0,
        "di_fragments": 0,
    }
    if not path or not os.path.isfile(path):
        return None, "parse_failed"
    try:
        with TemplateReader(path) as reader:
            for template in reader:
                counts["analysis_total_templates"] += 1
                classification = classify_pe_template(template)
                if not classification.valid:
                    counts["analysis_invalid_templates"] += 1
                    if classification.reason == "tlen_zero":
                        counts["analysis_tlen_zero_templates"] += 1
                    elif classification.reason == "orphan":
                        counts["analysis_orphan_templates"] += 1
                    elif classification.reason == "cross_contig":
                        counts["analysis_cross_contig_templates"] += 1
                    else:
                        counts["analysis_inconsistent_templates"] += 1
                    continue
                counts["analysis_valid_templates"] += 1
                length = classification.length
                in_any_window = False
                for name, low, high in windows:
                    if low <= length <= high:
                        counts[f"{name}_fragments"] += 1
                        in_any_window = True
                if in_any_window:
                    counts["analysis_classified_templates"] += 1
                else:
                    counts["analysis_unclassified_templates"] += 1
    except (OSError, ValueError, RuntimeError):
        return None, "parse_failed"
    if (
        counts["analysis_valid_templates"] + counts["analysis_invalid_templates"]
        != (counts["analysis_total_templates"])
    ):
        return None, "parse_failed"
    return counts, "ok" if counts["analysis_valid_templates"] else "no_valid_templates"


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
    p.add_argument("--analysis-bam", default="")
    p.add_argument("--stage-metrics", default="")
    p.add_argument("--insert-size-metrics", default="")
    p.add_argument("--danpos3-enabled", default="false")
    p.add_argument("--inps-enabled", default="false")
    p.add_argument("--sem-enabled", default="false")
    p.add_argument("--output", required=True)
    args = p.parse_args()

    windows = _windows(args)
    counts, analysis_status = _classify_analysis_bam(args.analysis_bam, windows)

    stage_metrics = None
    if args.stage_metrics:
        try:
            stage_metrics = load_stage_metrics(args.stage_metrics)
        except (OSError, ValueError):
            analysis_status = "parse_failed"

    def _stage(key):
        if stage_metrics is None:
            return NA
        return str(stage_metrics[key])

    def _count(key):
        if counts is None:
            return NA
        return str(counts[key])

    fraction_classified = (
        NA
        if counts is None
        else format_fraction(
            counts["analysis_classified_templates"],
            counts["analysis_valid_templates"],
        )
    )

    row = {
        "schema_version": "2",
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
            args.insert_size_metrics if os.path.isfile(args.insert_size_metrics) else NA
        ),
        "caller_danpos3_enabled": args.danpos3_enabled,
        "caller_inps_enabled": args.inps_enabled,
        "caller_sem_enabled": args.sem_enabled,
        "sub_reads": _count_reads(args.sub_bam),
        "mono_reads": _count_reads(args.mono_bam),
        "di_reads": _count_reads(args.di_bam),
        "analysis_bam": args.analysis_bam,
        "analysis_input_stage": "analysis_bam",
        "analysis_status": analysis_status,
        "analysis_total_templates": _count("analysis_total_templates"),
        "analysis_valid_templates": _count("analysis_valid_templates"),
        "analysis_invalid_templates": _count("analysis_invalid_templates"),
        "analysis_tlen_zero_templates": _count("analysis_tlen_zero_templates"),
        "analysis_orphan_templates": _count("analysis_orphan_templates"),
        "analysis_cross_contig_templates": _count("analysis_cross_contig_templates"),
        "analysis_inconsistent_templates": _count("analysis_inconsistent_templates"),
        "analysis_classified_templates": _count("analysis_classified_templates"),
        "analysis_unclassified_templates": _count("analysis_unclassified_templates"),
        "analysis_fraction_classified": fraction_classified,
        "window_overlap": "true" if _windows_overlap(windows) else "false",
        "sub_fragments": _count("sub_fragments"),
        "mono_fragments": _count("mono_fragments"),
        "di_fragments": _count("di_fragments"),
        "upstream_dropped_templates": _stage("dropped_templates"),
        "upstream_dropped_alignments": _stage("dropped_alignments"),
    }

    write_exact_tsv(args.output, MNASE_QC_V2_COLUMNS, row)


if __name__ == "__main__":
    main()
