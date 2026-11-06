#!/usr/bin/env python3
"""CUT&Tag fragment-size QC from the canonical analysis BAM.

For PE samples, one complete valid template contributes one length using
``abs(TLEN)``.  Legal negative TLEN is preserved.  Zero TLEN, orphan mates,
cross-contig pairs and mutually inconsistent mates are counted as invalid and
never enter a biological length bin.  The input BAM already reflects the
resolved duplicate/organelle/blacklist policy; this script never removes
0x400 records on its own.

For SE samples, fragment length is not defined and the v1 columns remain NA
with ``layout_not_supported``.

Usage:
    python3 scripts/calc_cuttag_fragment_size.py \
        --sample <sample_id> \
        --bam <analysis.bam> \
        --layout PE|SE \
        --duplicate-analysis exclude|retain \
        --blacklist-applied true|false \
        --stage-metrics <stage_metrics.json> \
        --output <output.tsv>
"""

from __future__ import annotations

import argparse
import os
import statistics
import sys

from qc_common import (
    FRAGMENT_QC_V2_COLUMNS,
    NA,
    load_stage_metrics,
    write_exact_tsv,
)
from qc_template import TemplateReader, classify_pe_template


def _compute_pe_stats(fragments):
    n = len(fragments)
    if n == 0:
        return {
            "fragment_count": "0",
            "fragment_mean": NA,
            "fragment_median": NA,
            "fragment_mode": NA,
            "fragment_min": NA,
            "fragment_max": NA,
            "fraction_lt_150": NA,
            "fraction_150_300": NA,
            "fraction_300_500": NA,
            "fraction_ge_500": NA,
            "fraction_lt_120": NA,
            "status": "no_fragments",
        }

    mean_val = statistics.mean(fragments)
    median_val = statistics.median(fragments)

    freq = {}
    for value in fragments:
        freq[value] = freq.get(value, 0) + 1
    max_count = max(freq.values())
    mode_val = min(value for value, count in freq.items() if count == max_count)

    lt_150 = sum(1 for value in fragments if value < 150)
    r150_300 = sum(1 for value in fragments if 150 <= value < 300)
    r300_500 = sum(1 for value in fragments if 300 <= value < 500)
    ge_500 = sum(1 for value in fragments if value >= 500)
    lt_120 = sum(1 for value in fragments if value < 120)

    return {
        "fragment_count": str(n),
        "fragment_mean": f"{mean_val:.1f}",
        "fragment_median": f"{median_val:.1f}",
        "fragment_mode": str(mode_val),
        "fragment_min": str(min(fragments)),
        "fragment_max": str(max(fragments)),
        "fraction_lt_150": f"{lt_150 / n:.4f}",
        "fraction_150_300": f"{r150_300 / n:.4f}",
        "fraction_300_500": f"{r300_500 / n:.4f}",
        "fraction_ge_500": f"{ge_500 / n:.4f}",
        "fraction_lt_120": f"{lt_120 / n:.4f}",
        "status": "ok",
    }


def _se_row(sample_id, duplicate_analysis, blacklist_applied):
    row = {column: NA for column in FRAGMENT_QC_V2_COLUMNS}
    row.update(
        {
            "schema_version": "2",
            "sample": sample_id,
            "layout": "SE",
            "fragment_count": NA,
            "status": "layout_not_supported",
            "input_stage": "analysis_bam",
            "duplicate_analysis": duplicate_analysis,
            "blacklist_applied": blacklist_applied,
        }
    )
    return row


def _counts_row(stage_metrics):
    if stage_metrics is None:
        return {
            "upstream_dropped_templates": NA,
            "upstream_dropped_alignments": NA,
            "upstream_alignments_total": NA,
            "upstream_alignments_kept": NA,
        }
    return {
        "upstream_dropped_templates": str(stage_metrics["dropped_templates"]),
        "upstream_dropped_alignments": str(stage_metrics["dropped_alignments"]),
        "upstream_alignments_total": str(stage_metrics["total_alignments"]),
        "upstream_alignments_kept": str(stage_metrics["kept_alignments"]),
    }


def main():
    parser = argparse.ArgumentParser(description="CUT&Tag fragment-size QC")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--layout", required=True, choices=("PE", "SE"))
    parser.add_argument(
        "--duplicate-analysis", default=NA, choices=("exclude", "retain", NA)
    )
    parser.add_argument(
        "--blacklist-applied", default=NA, choices=("true", "false", NA)
    )
    parser.add_argument("--stage-metrics", default="")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    stage_metrics = None
    if args.stage_metrics:
        stage_metrics = load_stage_metrics(args.stage_metrics)

    if args.layout == "SE":
        print(
            "CUT&Tag fragment-size QC: SE layout not supported; "
            "insert size cannot be reliably computed from SE alignments. "
            "Output will contain NA values.",
            file=sys.stderr,
        )
        row = _se_row(args.sample, args.duplicate_analysis, args.blacklist_applied)
    else:
        if not os.path.isfile(args.bam):
            print(f"ERROR: analysis BAM not found: {args.bam}", file=sys.stderr)
            sys.exit(1)

        fragments = []
        invalid = {
            "tlen_zero": 0,
            "orphan": 0,
            "cross_contig": 0,
            "inconsistent": 0,
        }
        templates_total = 0
        with TemplateReader(args.bam) as reader:
            for template in reader:
                templates_total += 1
                classification = classify_pe_template(template)
                if classification.valid:
                    fragments.append(classification.length)
                else:
                    invalid[classification.reason] += 1
            alignments_total = reader.alignments_total

        stats = _compute_pe_stats(fragments)
        row = {
            "schema_version": "2",
            "sample": args.sample,
            "layout": "PE",
            **stats,
            "input_stage": "analysis_bam",
            "duplicate_analysis": args.duplicate_analysis,
            "blacklist_applied": args.blacklist_applied,
            "templates_total": str(templates_total),
            "templates_valid": str(len(fragments)),
            "templates_invalid": str(templates_total - len(fragments)),
            "invalid_tlen_zero": str(invalid["tlen_zero"]),
            "invalid_orphan": str(invalid["orphan"]),
            "invalid_cross_contig": str(invalid["cross_contig"]),
            "invalid_inconsistent": str(invalid["inconsistent"]),
            "alignments_total": str(alignments_total),
            **_counts_row(stage_metrics),
        }
        invalid_total = sum(invalid.values())
        if templates_total != len(fragments) + invalid_total:
            raise AssertionError("template counts do not conserve")

    write_exact_tsv(args.output, FRAGMENT_QC_V2_COLUMNS, row)


if __name__ == "__main__":
    main()
