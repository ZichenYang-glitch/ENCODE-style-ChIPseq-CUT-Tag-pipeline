#!/usr/bin/env python3
"""Stage 6b pooled experiment QC summary writer.

Replaces the inline shell printf logic in pooled_experiment_qc_summary.
Both the Snakemake rule (qc.smk) and the content-level test invoke this
script so that the generation logic is shared and testable.

Usage:
    python3 scripts/pooled_qc_summary.py \
        --experiment exp1 \
        --assay chipseq \
        --target H3K27me3 \
        --peak-mode narrow \
        --inferred-histone-class broad_like \
        --expected-peak-mode broad \
        --peak-mode-status mismatch \
        --n-bioreps 2 \
        --bio-rep-labels 1,2 \
        --pooled-bam /path/to/pooled.bam \
        --pooled-peaks-dir /path/to/pooled_peaks \
        --pooled-fe-bdg /path/to/FE.bdg \
        --pooled-ppois-bdg /path/to/ppois.bdg \
        --signal-tracks-status disabled \
        --output /path/to/output.tsv
"""

import argparse
import os
import sys


def count_peaks(peaks_dir, peak_mode, experiment):
    """Count lines in the pooled peak file. Exits on missing file."""
    suffix = "broadPeak" if peak_mode == "broad" else "narrowPeak"
    pk_file = os.path.join(
        peaks_dir, f"{experiment}_pooled_peaks.{suffix}")
    if not os.path.isfile(pk_file):
        print(f"ERROR: pooled peak file not found: {pk_file}", file=sys.stderr)
        sys.exit(1)
    count = 0
    with open(pk_file) as fh:
        for _ in fh:
            count += 1
    return count


def main():
    parser = argparse.ArgumentParser(
        description="Stage 6b pooled experiment QC summary writer")
    parser.add_argument("--experiment", required=True)
    parser.add_argument("--assay", required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--peak-mode", required=True)
    parser.add_argument("--inferred-histone-class", required=True)
    parser.add_argument("--expected-peak-mode", required=True)
    parser.add_argument("--peak-mode-status", required=True)
    parser.add_argument("--n-bioreps", type=int, required=True)
    parser.add_argument("--bio-rep-labels", required=True)
    parser.add_argument("--pooled-bam", required=True)
    parser.add_argument("--pooled-peaks-dir", required=True)
    parser.add_argument("--pooled-fe-bdg", default="NA")
    parser.add_argument("--pooled-ppois-bdg", default="NA")
    parser.add_argument("--signal-tracks-status", default="disabled")
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    pk_count = count_peaks(args.pooled_peaks_dir, args.peak_mode,
                           args.experiment)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    header = [
        "experiment", "assay", "target", "inferred_histone_class",
        "expected_peak_mode", "configured_peak_mode", "peak_mode_status",
        "biological_replicates", "biological_replicate_labels",
        "pooled_bam", "pooled_peaks", "pooled_peak_count",
        "pooled_FE_bdg", "pooled_ppois_bdg", "signal_tracks_status",
    ]
    data = [
        args.experiment,
        args.assay,
        args.target,
        args.inferred_histone_class,
        args.expected_peak_mode,
        args.peak_mode,
        args.peak_mode_status,
        str(args.n_bioreps),
        args.bio_rep_labels,
        args.pooled_bam,
        args.pooled_peaks_dir,
        str(pk_count),
        args.pooled_fe_bdg,
        args.pooled_ppois_bdg,
        args.signal_tracks_status,
    ]

    with open(args.output, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.write("\t".join(data) + "\n")

    print(f"Pooled QC summary written: {args.output}")


if __name__ == "__main__":
    main()
