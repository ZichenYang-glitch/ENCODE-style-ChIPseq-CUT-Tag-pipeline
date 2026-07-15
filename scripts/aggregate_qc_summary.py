#!/usr/bin/env python3
"""Aggregate per-sample QC summary TSVs into a project-level summary.

Replaces the shell head/tail concatenation in the stage3_qc_summary Snakemake rule.
Validates header consistency across all input files and concatenates data rows.
Uses the same 37-column header as assemble_qc_summary.py.

Usage:
    python3 scripts/aggregate_qc_summary.py \\
        --output results/multiqc/stage3_qc_summary.tsv \\
        results/S1/01_qc/S1.qc_summary.tsv \\
        results/S2/01_qc/S2.qc_summary.tsv
"""

import argparse
import csv
import os
import sys

_NA = "NA"

# Must match assemble_qc_summary.py exactly.
# test/scripts/test_qc_summary.py catches drift.
_QC_SUMMARY_COLUMNS = [
    "sample",
    "assay",
    "target",
    "genome",
    "layout",
    "peak_mode",
    "use_control",
    "control_type",
    "final_bam",
    "peaks",
    "blacklist",
    "blacklist_filtered_bam",
    "blacklist_filtered_peaks",
    "total_reads",
    "reads_in_peaks",
    "frip",
    "peak_count",
    "blacklist_filtered_peak_count",
    "metrics_source",
    "unpaired_reads_examined",
    "read_pairs_examined",
    "secondary_or_supplementary_reads",
    "unmapped_reads",
    "unpaired_read_duplicates",
    "read_pair_duplicates",
    "read_pair_optical_duplicates",
    "percent_duplication",
    "estimated_library_size",
    "total_reads_examined",
    "duplicate_reads_estimate",
    "total_fragments",
    "distinct_fragments",
    "one_read_fragments",
    "two_read_fragments",
    "nrf",
    "pbc1",
    "pbc2",
]

_QCOL = {c: i for i, c in enumerate(_QC_SUMMARY_COLUMNS)}


def _check_header_and_yield_rows(filepath):
    """Validate header and yield data rows from *filepath*.

    sys.exit(1) on header mismatch.
    """
    fh = open(filepath, newline="")
    reader = csv.reader(fh, delimiter="\t")
    try:
        header = next(reader)
    except StopIteration:
        fh.close()
        print(
            f"ERROR: empty or unreadable file: {filepath}",
            file=sys.stderr,
        )
        sys.exit(1)

    if header != _QC_SUMMARY_COLUMNS:
        fh.close()
        print(
            f"ERROR: header mismatch in {filepath}",
            file=sys.stderr,
        )
        print(f"  Expected {len(_QC_SUMMARY_COLUMNS)} columns", file=sys.stderr)
        print(f"  Got      {len(header)} columns", file=sys.stderr)
        for i, (exp, got) in enumerate(zip(_QC_SUMMARY_COLUMNS, header)):
            if exp != got:
                print(
                    f"  First diff at column {i}: expected {exp!r}, got {got!r}",
                    file=sys.stderr,
                )
                break
        sys.exit(1)

    for row in reader:
        if row:
            yield row
    fh.close()


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate per-sample QC summary TSVs into project-level summary"
    )
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("inputs", nargs="*", help="Per-sample qc_summary.tsv files")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    with open(args.output, "w", newline="") as out_fh:
        writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
        writer.writerow(_QC_SUMMARY_COLUMNS)

        if not args.inputs:
            # Zero inputs — header-only output (preserves legacy behavior)
            return

        for filepath in args.inputs:
            if not os.path.isfile(filepath):
                print(
                    f"ERROR: input file not found: {filepath}",
                    file=sys.stderr,
                )
                sys.exit(1)

            reader = _check_header_and_yield_rows(filepath)
            for row in reader:
                if row:  # skip empty lines
                    writer.writerow(row)


if __name__ == "__main__":
    main()
