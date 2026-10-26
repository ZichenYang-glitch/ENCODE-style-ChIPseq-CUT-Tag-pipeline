#!/usr/bin/env python3
"""General PE analysis-template counts from the canonical analysis BAM.

This component is the assay-independent source for the per-sample analysis
template counts recorded in the analysis policy.  It runs for every sample and
never depends on CUT&Tag-specific optional QC.

One complete valid PE template is counted once and contributes one length
using ``abs(TLEN)``; a legal negative TLEN is preserved.  TLEN=0, orphan mates,
cross-contig pairs and mutually inconsistent mates are counted in separate
invalid-reason columns and never enter a biological length.  The totals are
conserved: ``total = valid + invalid`` and
``invalid = tlen_zero + orphan + cross_contig + inconsistent``.
Alignments already dropped by the upstream analysis-view stage are reported
separately from invalid templates at the final analysis input, so neither is
double-counted.

Statuses:
    ok                  PE counts were computed and at least one template is valid.
    no_valid_templates  PE counts were computed, totals are conserved, and the
                        legitimate result is zero valid templates.
    not_applicable      SE input; template counts are NA by contract.
    parse_failed        the analysis BAM is missing/unreadable/malformed or the
                        computed counts do not conserve; all counts are NA and
                        never a fabricated zero.

Usage:
    python3 scripts/calc_analysis_template_metrics.py \
        --sample <sample_id> \
        --layout PE|SE \
        --bam <analysis.bam> \
        --stage-metrics <stage_metrics.json> \
        --output <output.tsv>
"""

from __future__ import annotations

import argparse
import os
import sys

from qc_common import (
    ANALYSIS_TEMPLATE_METRICS_COLUMNS,
    NA,
    load_stage_metrics,
    write_exact_tsv,
)
from qc_template import TemplateReader, classify_pe_template


def _upstream_counts(stage_metrics):
    if stage_metrics is None:
        return NA, NA
    return (
        str(stage_metrics["dropped_templates"]),
        str(stage_metrics["dropped_alignments"]),
    )


def _base_row(sample_id, layout, status, stage_metrics):
    upstream_templates, upstream_alignments = _upstream_counts(stage_metrics)
    row = {column: NA for column in ANALYSIS_TEMPLATE_METRICS_COLUMNS}
    row.update(
        {
            "schema_version": "2",
            "sample": sample_id,
            "layout": layout,
            "analysis_status": status,
            "input_stage": "analysis_bam",
            "upstream_dropped_templates": upstream_templates,
            "upstream_dropped_alignments": upstream_alignments,
        }
    )
    return row


def _classify_bam(bam):
    """Return (counts, status) for one analysis BAM; counts is None on failure."""
    counts = {
        "analysis_total_templates": 0,
        "analysis_valid_templates": 0,
        "analysis_invalid_templates": 0,
        "analysis_tlen_zero_templates": 0,
        "analysis_orphan_templates": 0,
        "analysis_cross_contig_templates": 0,
        "analysis_inconsistent_templates": 0,
    }
    if not bam or not os.path.isfile(bam):
        return None, "parse_failed"
    try:
        with TemplateReader(bam) as reader:
            for template in reader:
                counts["analysis_total_templates"] += 1
                classification = classify_pe_template(template)
                if classification.valid:
                    counts["analysis_valid_templates"] += 1
                    continue
                counts["analysis_invalid_templates"] += 1
                if classification.reason == "tlen_zero":
                    counts["analysis_tlen_zero_templates"] += 1
                elif classification.reason == "orphan":
                    counts["analysis_orphan_templates"] += 1
                elif classification.reason == "cross_contig":
                    counts["analysis_cross_contig_templates"] += 1
                else:
                    counts["analysis_inconsistent_templates"] += 1
            counts["alignments_total"] = reader.alignments_total
    except (OSError, ValueError, RuntimeError):
        return None, "parse_failed"

    invalid_reasons = sum(
        counts[key]
        for key in (
            "analysis_tlen_zero_templates",
            "analysis_orphan_templates",
            "analysis_cross_contig_templates",
            "analysis_inconsistent_templates",
        )
    )
    if (
        counts["analysis_valid_templates"] + counts["analysis_invalid_templates"]
        != counts["analysis_total_templates"]
        or invalid_reasons != counts["analysis_invalid_templates"]
    ):
        return None, "parse_failed"
    if counts["analysis_valid_templates"] == 0:
        return counts, "no_valid_templates"
    return counts, "ok"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--layout", required=True, choices=("PE", "SE"))
    parser.add_argument("--bam", required=True)
    parser.add_argument("--stage-metrics", default="")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    stage_metrics = None
    if args.stage_metrics:
        try:
            stage_metrics = load_stage_metrics(args.stage_metrics)
        except (OSError, ValueError) as exc:
            print(f"ERROR: stage metrics are invalid: {exc}", file=sys.stderr)
            sys.exit(1)

    if args.layout == "SE":
        print(
            "Analysis template metrics: SE layout is not applicable; "
            "template counts remain NA.",
            file=sys.stderr,
        )
        row = _base_row(args.sample, "SE", "not_applicable", stage_metrics)
    else:
        counts, status = _classify_bam(args.bam)
        row = _base_row(args.sample, "PE", status, stage_metrics)
        if counts is not None:
            row.update({key: str(value) for key, value in counts.items()})

    write_exact_tsv(args.output, ANALYSIS_TEMPLATE_METRICS_COLUMNS, row)


if __name__ == "__main__":
    main()
