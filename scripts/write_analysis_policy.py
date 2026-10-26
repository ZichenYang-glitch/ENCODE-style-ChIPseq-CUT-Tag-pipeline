#!/usr/bin/env python3
"""Write the adapter-owned, versioned per-sample analysis policy record.

The record contains the values actually resolved for this run, the scientific
input stage behind each QC product, and stage counts.  It deliberately stores
roles and identifiers only: no absolute workspace paths, environment values,
or private configuration are copied into the artifact.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

from qc_common import (
    FRAGMENT_QC_V1_COLUMNS,
    FRAGMENT_QC_V2_COLUMNS,
    MITO_QC_V1_COLUMNS,
    MITO_QC_V2_COLUMNS,
    MNASE_QC_V1_COLUMNS,
    MNASE_QC_V2_COLUMNS,
    NA,
    load_stage_metrics,
    read_versioned_tsv,
)

SCHEMA_VERSION = 1
QC_SCHEMA_VERSIONS = {
    "sample": "2",
    "mnase": "2",
    "fragment": "2",
    "mito": "2",
}


def _optional_int(value):
    if value is None or value == NA or value == "":
        return None
    return int(value)


def _matched_contigs(value):
    if not value or value == NA:
        return []
    return [item for item in value.split(",") if item]


def _read_qc(path, v1, v2):
    if not path:
        return None, None
    try:
        version, row = read_versioned_tsv(path, v1, v2)
    except (OSError, ValueError):
        return "invalid", None
    return version, row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--assay", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--profile", required=True, choices=("legacy", "updated"))
    parser.add_argument("--policy-revision", type=int, required=True)
    parser.add_argument(
        "--duplicate-analysis", required=True, choices=("exclude", "retain")
    )
    parser.add_argument("--max-insert", required=True)
    parser.add_argument("--organelle-mode", required=True)
    parser.add_argument("--blacklist-mode", required=True)
    parser.add_argument("--blacklist-applied", required=True, choices=("true", "false"))
    parser.add_argument("--stage-metrics", required=True)
    parser.add_argument("--mito", required=True)
    parser.add_argument("--fragment-qc", default="")
    parser.add_argument("--mnase-qc", default="")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    stage_metrics = load_stage_metrics(args.stage_metrics)

    mito_version, mito = _read_qc(args.mito, MITO_QC_V1_COLUMNS, MITO_QC_V2_COLUMNS)
    if mito is None:
        organelle_status = (
            "parse_failed" if mito_version == "invalid" else "not_applicable"
        )
        matched = []
        organelle_reads = None
        organelle_total = None
    else:
        organelle_status = mito["status"]
        matched = (
            _matched_contigs(mito.get("matched_contigs")) if mito_version == "2" else []
        )
        organelle_reads = _optional_int(mito["organelle_reads"])
        organelle_total = _optional_int(mito["total_mapped_reads"])

    fragment_version, fragment = _read_qc(
        args.fragment_qc, FRAGMENT_QC_V1_COLUMNS, FRAGMENT_QC_V2_COLUMNS
    )
    mnase_version, mnase = _read_qc(
        args.mnase_qc, MNASE_QC_V1_COLUMNS, MNASE_QC_V2_COLUMNS
    )

    if fragment is not None and fragment_version == "2":
        valid_templates = _optional_int(fragment["templates_valid"])
        invalid_templates = _optional_int(fragment["templates_invalid"])
    elif mnase is not None and mnase_version == "2":
        valid_templates = _optional_int(mnase["analysis_valid_templates"])
        invalid_templates = _optional_int(mnase["analysis_invalid_templates"])
    else:
        valid_templates = None
        invalid_templates = None

    blacklist_dropped_templates = int(
        stage_metrics["dropped_templates_by_reason"].get("blacklist", 0)
    )
    blacklist_dropped_alignments = int(
        stage_metrics["dropped_alignments_by_reason"].get("blacklist", 0)
    )

    input_stages = {
        "fragment_size": "analysis_bam" if fragment is not None else "not_applicable",
        "mnase_classification": (
            "analysis_bam" if mnase is not None else "not_applicable"
        ),
        "duplicate_metrics": "pre_duplicate_markdup_bam",
        "nrf_pbc": "pre_duplicate_mapq_bam",
        "organelle": "pre_duplicate_mapq_bam",
        "blacklist": (
            "final_bam_to_analysis_bam"
            if args.blacklist_applied == "true"
            else "not_applicable"
        ),
        "peaks": (
            "analysis_bam"
            if args.assay in ("chipseq", "cuttag", "atac")
            else "not_applicable"
        ),
    }

    stage_counts = {
        "analysis_input_alignments": int(stage_metrics["total_alignments"]),
        "analysis_kept_alignments": int(stage_metrics["kept_alignments"]),
        "analysis_valid_templates": valid_templates,
        "analysis_invalid_templates": invalid_templates,
        "upstream_dropped_templates": int(stage_metrics["dropped_templates"]),
        "upstream_dropped_alignments": int(stage_metrics["dropped_alignments"]),
        "blacklist_dropped_templates": blacklist_dropped_templates,
        "blacklist_dropped_alignments": blacklist_dropped_alignments,
        "organelle_reads": organelle_reads,
        "organelle_total_mapped_reads": organelle_total,
    }

    record = {
        "schema_version": SCHEMA_VERSION,
        "scientific_policy_revision": args.policy_revision,
        "sample": args.sample,
        "assay": args.assay,
        "profile": args.profile,
        "resolved_duplicate_analysis": args.duplicate_analysis,
        "resolved_max_insert": (
            args.max_insert
            if args.max_insert == "not_applicable"
            else int(args.max_insert)
        ),
        "organelle_mode": args.organelle_mode,
        "organelle_status": organelle_status,
        "matched_organelle_contigs": matched,
        "blacklist_mode": args.blacklist_mode,
        "blacklist_applied": args.blacklist_applied == "true",
        "input_stages": input_stages,
        "stage_counts": stage_counts,
        "reference_summary": {
            "genome": args.genome,
            "organelle_chroms_explicit": bool(
                mito.get("organelle_chroms_explicit") == "true"
                if mito is not None and mito_version == "2"
                else False
            ),
        },
        "qc_schema_versions": QC_SCHEMA_VERSIONS,
    }

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", encoding="utf-8", newline="\n") as handle:
        json.dump(record, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
