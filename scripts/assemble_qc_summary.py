#!/usr/bin/env python3
"""Assemble the versioned per-sample QC summary TSV.

Reads component QC tables by exact header and writes the fixed v2 sample
contract.  Version 1 columns keep their original order after the leading
schema version; new analysis fields follow.  Values come from the same
analysis-policy record used by the workflow, so policy, stage counts and
user-visible metrics describe one resolved run.

When the fragment-specific component is absent, the policy's general
analysis-template counts populate ``analysis_valid_fragments`` and
``analysis_invalid_templates``.  Fragment-only reason counts and upstream
dropped counts remain NA in that case; they are available from the fragment
component when CUT&Tag fragment QC is enabled.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

from qc_common import (
    read_cross_correlation_detail,
    FRAGMENT_QC_V1_COLUMNS,
    FRAGMENT_QC_V2_COLUMNS,
    MITO_QC_V1_COLUMNS,
    MITO_QC_V2_COLUMNS,
    NA,
    SAMPLE_QC_V2_COLUMNS,
    read_versioned_tsv,
    write_exact_tsv,
)

from encode_pipeline.qc.policy import (
    AnalysisPolicyError,
    parse_analysis_policy_bytes,
)
from encode_pipeline.qc.contract import (
    CROSS_CORRELATION_EXECUTION_STATUSES,
    CROSS_CORRELATION_INTERPRETATION_STATUSES,
)

_QC_SUMMARY_COLUMNS = SAMPLE_QC_V2_COLUMNS


def _read_tsv_columns(filepath, *column_names):
    """Read named columns from a TSV file via csv.DictReader.

    Returns a dict mapping each column name to its value, or NA if the
    column is missing or empty.  Reads only the first data row.
    """
    result = {name: NA for name in column_names}
    if not os.path.isfile(filepath):
        print(f"ERROR: required input file not found: {filepath}", file=sys.stderr)
        sys.exit(1)
    with open(filepath, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            return result
        if len(set(reader.fieldnames)) != len(reader.fieldnames):
            print(f"ERROR: duplicate columns in {filepath}", file=sys.stderr)
            sys.exit(1)
        for row in reader:
            for name in column_names:
                val = row.get(name, "").strip()
                if val:
                    result[name] = val
            break  # only first data row
    return result


def _read_policy(path):
    try:
        with open(path, "rb") as handle:
            return parse_analysis_policy_bytes(handle.read())
    except (OSError, AnalysisPolicyError) as exc:
        print(f"ERROR: analysis policy is invalid: {exc}", file=sys.stderr)
        sys.exit(1)


def _read_optional_versioned(path, v1_columns, v2_columns):
    if not path:
        return None, None
    try:
        return read_versioned_tsv(path, v1_columns, v2_columns)
    except (OSError, ValueError) as exc:
        print(f"ERROR: QC source is invalid: {exc}", file=sys.stderr)
        sys.exit(1)


def _cross_correlation(path, enabled, detail_path=None):
    if not enabled:
        return NA, NA, "not_enabled"
    if detail_path:
        try:
            version, detail = read_cross_correlation_detail(detail_path)
        except ValueError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            sys.exit(1)
        if version == "1":
            return detail.get("nsc", NA), detail.get("rsc", NA), "reported"
        if detail is not None:
            execution = detail.get("execution_status", "")
            interpretation = detail.get("interpretation_status", "")
            if (
                execution not in CROSS_CORRELATION_EXECUTION_STATUSES
                or interpretation not in CROSS_CORRELATION_INTERPRETATION_STATUSES
            ):
                print(
                    "ERROR: cross-correlation detail status is invalid",
                    file=sys.stderr,
                )
                sys.exit(1)
            if execution == "reported":
                nsc = detail.get("nsc", NA) or NA
                rsc = detail.get("rsc", NA) or NA
                return nsc, rsc, interpretation
            return NA, NA, execution
    if not path or not os.path.isfile(path):
        return NA, NA, "missing"
    from parse_cross_correlation import parse_cc_qc_file

    parsed = parse_cc_qc_file(path)
    nsc = parsed.get("nsc")
    rsc = parsed.get("rsc")
    if nsc is None or rsc is None:
        return NA, NA, "parse_failed"
    return str(nsc), str(rsc), "reported"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--assay", required=True)
    parser.add_argument("--target", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--layout", required=True)
    parser.add_argument("--peak-mode", required=True)
    parser.add_argument("--use-control", required=True)
    parser.add_argument("--control-type", required=True)
    parser.add_argument("--final-bam", required=True)
    parser.add_argument("--analysis-bam", required=True)
    parser.add_argument("--peaks-file", required=True)
    parser.add_argument("--has-blacklist", required=True)
    parser.add_argument("--blacklist", default=NA)
    parser.add_argument("--bl-bam", default=NA)
    parser.add_argument("--bl-peaks", default=NA)
    parser.add_argument("--peak-counts", required=True)
    parser.add_argument("--frip", required=True)
    parser.add_argument("--library-complexity", required=True)
    parser.add_argument("--nrf-pbc", required=True)
    parser.add_argument("--mito-tsv", required=True)
    parser.add_argument("--analysis-policy", required=True)
    parser.add_argument("--fragment-qc", default="")
    parser.add_argument("--cross-correlation-qc", default="")
    parser.add_argument("--cross-correlation-detail", default="")
    parser.add_argument(
        "--cross-correlation-enabled", default="false", choices=("true", "false")
    )
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    policy = _read_policy(args.analysis_policy)
    if policy["sample"] != args.sample or policy["assay"] != args.assay:
        print("ERROR: analysis policy identity is inconsistent", file=sys.stderr)
        sys.exit(1)

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

    mito_version, mito = _read_optional_versioned(
        args.mito_tsv, MITO_QC_V1_COLUMNS, MITO_QC_V2_COLUMNS
    )
    fragment_version, fragment = _read_optional_versioned(
        args.fragment_qc, FRAGMENT_QC_V1_COLUMNS, FRAGMENT_QC_V2_COLUMNS
    )

    if mito is None:
        organelle = {
            "organelle_reads": NA,
            "total_mapped_reads": NA,
            "fraction": NA,
            "status": "missing",
            "matched_contigs": NA,
        }
    else:
        organelle = {
            "organelle_reads": mito["organelle_reads"],
            "total_mapped_reads": mito["total_mapped_reads"],
            "fraction": mito["fraction"],
            "status": mito["status"],
            "matched_contigs": (
                mito.get("matched_contigs", NA) if mito_version == "2" else NA
            ),
        }

    bl_bam = args.bl_bam
    bl_peaks = args.bl_peaks
    bl_path = args.blacklist
    bl_peak_count = peak_data["blacklist_filtered_peaks"]
    if args.has_blacklist == "no":
        bl_bam = NA
        bl_peaks = NA
        bl_path = NA
        bl_peak_count = NA

    blacklist_applied = policy["blacklist_applied"] and args.has_blacklist != "no"
    delivered_peak_count = (
        peak_data["blacklist_filtered_peaks"]
        if blacklist_applied
        else peak_data["peaks"]
    )
    if blacklist_applied and delivered_peak_count == NA:
        delivered_peak_count = NA

    fragment_fields = {
        "analysis_valid_fragments": NA,
        "analysis_invalid_templates": NA,
        "analysis_tlen_zero_templates": NA,
        "analysis_orphan_templates": NA,
        "analysis_cross_contig_templates": NA,
        "analysis_inconsistent_templates": NA,
        "analysis_upstream_dropped_templates": NA,
        "analysis_upstream_dropped_alignments": NA,
    }
    if fragment is not None and fragment_version == "2":
        fragment_fields = {
            "analysis_valid_fragments": fragment["templates_valid"],
            "analysis_invalid_templates": fragment["templates_invalid"],
            "analysis_tlen_zero_templates": fragment["invalid_tlen_zero"],
            "analysis_orphan_templates": fragment["invalid_orphan"],
            "analysis_cross_contig_templates": fragment["invalid_cross_contig"],
            "analysis_inconsistent_templates": fragment["invalid_inconsistent"],
            "analysis_upstream_dropped_templates": fragment[
                "upstream_dropped_templates"
            ],
            "analysis_upstream_dropped_alignments": fragment[
                "upstream_dropped_alignments"
            ],
        }
    elif fragment is not None:
        fragment_fields["analysis_valid_fragments"] = fragment["fragment_count"]
    else:
        # The fixed sample v2 header has no general analysis-template status
        # column.  When fragment QC is not applicable or not enabled, use the
        # policy's general analysis-template counts so ChIP-seq, ATAC-seq and
        # CUT&Tag with fragment QC off still report valid/invalid templates.
        # Reason-specific and upstream fragment fields remain NA because they
        # live in the fragment-specific component.
        general_valid = policy["stage_counts"].get("analysis_valid_templates")
        general_invalid = policy["stage_counts"].get("analysis_invalid_templates")
        if general_valid is not None:
            fragment_fields["analysis_valid_fragments"] = str(general_valid)
        if general_invalid is not None:
            fragment_fields["analysis_invalid_templates"] = str(general_invalid)

    nsc, rsc, cc_status = _cross_correlation(
        args.cross_correlation_qc,
        args.cross_correlation_enabled == "true",
        args.cross_correlation_detail,
    )

    row = {
        "schema_version": "2",
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
        "peak_count": delivered_peak_count,
        "blacklist_filtered_peak_count": bl_peak_count,
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
        "total_fragments": nrf_data["total_fragments"],
        "distinct_fragments": nrf_data["distinct_fragments"],
        "one_read_fragments": nrf_data["one_read_fragments"],
        "two_read_fragments": nrf_data["two_read_fragments"],
        "nrf": nrf_data["nrf"],
        "pbc1": nrf_data["pbc1"],
        "pbc2": nrf_data["pbc2"],
        "analysis_bam": args.analysis_bam,
        "analysis_input_stage": "analysis_bam",
        **fragment_fields,
        "raw_peak_count": peak_data["peaks"],
        "delivered_peak_count": delivered_peak_count,
        "duplicate_analysis": policy["resolved_duplicate_analysis"],
        "resolved_max_insert": str(policy["resolved_max_insert"]),
        "blacklist_mode": policy["blacklist_mode"],
        "blacklist_applied": "true" if blacklist_applied else "false",
        "organelle_mode": policy["organelle_mode"],
        "organelle_status": organelle["status"],
        "organelle_reads": organelle["organelle_reads"],
        "organelle_total_mapped_reads": organelle["total_mapped_reads"],
        "organelle_fraction": organelle["fraction"],
        "matched_organelle_contigs": organelle["matched_contigs"],
        "cross_correlation_nsc": nsc,
        "cross_correlation_rsc": rsc,
        "cross_correlation_status": cc_status,
        "scientific_policy_revision": str(policy["scientific_policy_revision"]),
    }

    write_exact_tsv(args.output, SAMPLE_QC_V2_COLUMNS, row)


if __name__ == "__main__":
    main()
