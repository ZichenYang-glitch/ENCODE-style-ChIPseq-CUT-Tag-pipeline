#!/usr/bin/env python3
"""Assemble per-sample QC summary TSV from component QC metric files.

Replaces the shell printf + tail/cut logic in the qc_summary Snakemake rule.
Reads input TSVs by header (csv.DictReader), handles missing optional metrics
with "NA", and preserves the exact 37-column output contract.

Usage:
    python3 scripts/assemble_qc_summary.py \\
        --sample S1 --assay chipseq --target CTCF --genome hs --layout PE \\
        --peak-mode narrow --use-control False --control-type none \\
        --final-bam results/S1/02_align/S1.final.bam \\
        --peaks-file results/S1/04_peaks/S1/S1_peaks.narrowPeak \\
        --has-blacklist yes --blacklist /opt/genomes/bl.bed \\
        --bl-bam results/S1/02_align/S1.blacklist_filtered.bam \\
        --bl-peaks results/S1/04_peaks/S1_bl/S1_peaks.blacklist_filtered.narrowPeak \\
        --peak-counts results/S1/01_qc/S1.peak_counts.tsv \\
        --frip results/S1/01_qc/S1.frip.tsv \\
        --library-complexity results/S1/01_qc/S1.library_complexity.tsv \\
        --nrf-pbc results/S1/01_qc/S1.nrf_pbc.tsv \\
        --output results/S1/01_qc/S1.qc_summary.tsv
"""

import argparse
import csv
import os
import sys

_NA = "NA"

# Exact 37-column output order preserved from the legacy shell rule (lines 738-773).
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


def _read_tsv_columns(filepath, *column_names):
    """Read named columns from a TSV file via csv.DictReader.

    Returns a dict mapping each column name to its value, or _NA if the
    column is missing or empty.  Reads only the first data row.
    """
    result = {name: _NA for name in column_names}
    if not os.path.isfile(filepath):
        print(f"ERROR: required input file not found: {filepath}", file=sys.stderr)
        sys.exit(1)
    with open(filepath, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            return result
        for row in reader:
            for name in column_names:
                val = row.get(name, "").strip()
                if val:
                    result[name] = val
            break  # only first data row
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Assemble per-sample QC summary TSV from component metrics"
    )
    # Sample metadata
    parser.add_argument("--sample", required=True)
    parser.add_argument("--assay", required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--layout", required=True)
    parser.add_argument("--peak-mode", required=True)
    parser.add_argument("--use-control", required=True)
    parser.add_argument("--control-type", required=True)
    parser.add_argument("--final-bam", required=True)
    parser.add_argument("--peaks-file", required=True)
    parser.add_argument("--has-blacklist", required=True)
    parser.add_argument("--blacklist", default=_NA)
    parser.add_argument("--bl-bam", default=_NA)
    parser.add_argument("--bl-peaks", default=_NA)
    # Component TSV inputs
    parser.add_argument("--peak-counts", required=True)
    parser.add_argument("--frip", required=True)
    parser.add_argument("--library-complexity", required=True)
    parser.add_argument("--nrf-pbc", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()

    # ---- Read component TSVs by header ----

    peak_data = _read_tsv_columns(args.peak_counts, "peaks", "blacklist_filtered_peaks")
    frip_data = _read_tsv_columns(args.frip, "total_reads", "reads_in_peaks", "frip")
    lc_data = _read_tsv_columns(
        args.library_complexity,
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
    )
    nrf_data = _read_tsv_columns(
        args.nrf_pbc,
        "total_fragments",
        "distinct_fragments",
        "one_read_fragments",
        "two_read_fragments",
        "nrf",
        "pbc1",
        "pbc2",
    )

    # ---- Apply blacklist override (matches legacy shell lines 731-734) ----

    bl_bam = args.bl_bam
    bl_peaks = args.bl_peaks
    bl_path = args.blacklist
    bl_peak_count = peak_data["blacklist_filtered_peaks"]

    if args.has_blacklist == "no":
        bl_bam = _NA
        bl_peaks = _NA
        bl_path = _NA
        bl_peak_count = _NA

    # ---- Assemble into 37-column row ----

    row = {
        "sample": args.sample,
        "assay": args.assay,
        "target": args.target,
        "genome": args.genome,
        "layout": args.layout,
        "peak_mode": args.peak_mode,
        "use_control": args.use_control,
        "control_type": args.control_type,
        "final_bam": args.final_bam,
        "peaks": args.peaks_file,
        "blacklist": bl_path,
        "blacklist_filtered_bam": bl_bam,
        "blacklist_filtered_peaks": bl_peaks,
        "total_reads": frip_data["total_reads"],
        "reads_in_peaks": frip_data["reads_in_peaks"],
        "frip": frip_data["frip"],
        "peak_count": peak_data["peaks"],
        "blacklist_filtered_peak_count": bl_peak_count,
        # Library complexity
        "metrics_source": lc_data["metrics_source"],
        "unpaired_reads_examined": lc_data["unpaired_reads_examined"],
        "read_pairs_examined": lc_data["read_pairs_examined"],
        "secondary_or_supplementary_reads": lc_data["secondary_or_supplementary_reads"],
        "unmapped_reads": lc_data["unmapped_reads"],
        "unpaired_read_duplicates": lc_data["unpaired_read_duplicates"],
        "read_pair_duplicates": lc_data["read_pair_duplicates"],
        "read_pair_optical_duplicates": lc_data["read_pair_optical_duplicates"],
        "percent_duplication": lc_data["percent_duplication"],
        "estimated_library_size": lc_data["estimated_library_size"],
        "total_reads_examined": lc_data["total_reads_examined"],
        "duplicate_reads_estimate": lc_data["duplicate_reads_estimate"],
        # NRF/PBC
        "total_fragments": nrf_data["total_fragments"],
        "distinct_fragments": nrf_data["distinct_fragments"],
        "one_read_fragments": nrf_data["one_read_fragments"],
        "two_read_fragments": nrf_data["two_read_fragments"],
        "nrf": nrf_data["nrf"],
        "pbc1": nrf_data["pbc1"],
        "pbc2": nrf_data["pbc2"],
    }

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=_QC_SUMMARY_COLUMNS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
