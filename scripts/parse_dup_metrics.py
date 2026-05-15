#!/usr/bin/env python3
"""Parse Picard MarkDuplicates metrics into a library-complexity TSV.

Stage 3b-1: duplication-derived QC metrics.
Handles Picard MarkDuplicates output; falls back to NA fields when metrics
are unavailable (e.g., samtools fallback text file).

Usage:
    python3 scripts/parse_dup_metrics.py \\
        --sample SAMPLE \\
        --metrics sample.dup_metrics.txt \\
        --output sample.library_complexity.tsv

Output columns:
    sample  metrics_source  unpaired_reads_examined  read_pairs_examined
    secondary_or_supplementary_reads  unmapped_reads  unpaired_read_duplicates
    read_pair_duplicates  read_pair_optical_duplicates  percent_duplication
    estimated_library_size  total_reads_examined  duplicate_reads_estimate
"""

import argparse
import os
import sys

# Expected Picard metric field names (case-insensitive match)
_PICARD_FIELDS = {
    "unpaired_reads_examined": (
        "unpaired_reads_examined",
        "unpaired_reads",
    ),
    "read_pairs_examined": (
        "read_pairs_examined",
        "read_pairs",
    ),
    "secondary_or_supplementary_reads": (
        "secondary_or_supplementary_rds",
        "secondary_or_supplementary_reads",
    ),
    "unmapped_reads": (
        "unmapped_reads",
    ),
    "unpaired_read_duplicates": (
        "unpaired_read_duplicates",
    ),
    "read_pair_duplicates": (
        "read_pair_duplicates",
    ),
    "read_pair_optical_duplicates": (
        "read_pair_optical_duplicates",
    ),
    "percent_duplication": (
        "percent_duplication",
    ),
    "estimated_library_size": (
        "estimated_library_size",
    ),
}

_NA = "NA"


def _empty_fields() -> dict[str, str]:
    """Return all output metric fields initialized to NA."""
    fields = {k: _NA for k in _PICARD_FIELDS}
    fields["total_reads_examined"] = _NA
    fields["duplicate_reads_estimate"] = _NA
    return fields


def _is_fallback_metrics(filepath: str) -> bool:
    """Return True if the metrics file is a fallback (non-Picard) message."""
    try:
        with open(filepath) as fh:
            head = fh.read(512).lower()
    except OSError:
        return True
    # Fallback patterns from duplicate_handling rule
    if "picard not found" in head:
        return True
    if "samtools markdup used" in head:
        return True
    if "duplicate metrics unavailable" in head:
        return True
    return False


def _parse_picard_metrics(filepath: str) -> dict[str, str]:
    """Parse Picard MarkDuplicates metrics into a dict of field → value.

    Picard format is tab-delimited with a header line and one data line,
    possibly preceded by comment lines starting with '#'.
    """
    result: dict[str, str] = {k: _NA for k in _PICARD_FIELDS}

    with open(filepath) as fh:
        header = None
        data = None
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("##") or line.startswith("#"):
                # Comment / METRICS CLASS line — skip
                continue
            if header is None:
                header = [h.strip().lower() for h in line.split("\t")]
                continue
            if data is None:
                data = [d.strip() for d in line.split("\t")]
                break

    if header is None or data is None:
        return result

    # Build a lookup of header index → values
    for i, col in enumerate(header):
        col_lower = col.lower()
        for key, aliases in _PICARD_FIELDS.items():
            if col_lower in aliases and i < len(data) and data[i]:
                result[key] = data[i]
                break

    return result


def _has_parsed_metric(parsed: dict[str, str]) -> bool:
    """Return True if at least one expected Picard metric was parsed."""
    return any(parsed.get(key, _NA) != _NA for key in _PICARD_FIELDS)


def _compute_derived(parsed: dict[str, str]) -> dict[str, str]:
    """Add derived fields: total_reads_examined, duplicate_reads_estimate."""
    result = dict(parsed)

    urp = parsed.get("unpaired_reads_examined", _NA)
    rp = parsed.get("read_pairs_examined", _NA)
    urd = parsed.get("unpaired_read_duplicates", _NA)
    rpd = parsed.get("read_pair_duplicates", _NA)

    # total_reads_examined = unpaired + 2 * read_pairs
    if urp != _NA and rp != _NA:
        try:
            result["total_reads_examined"] = str(
                int(urp) + 2 * int(rp)
            )
        except (ValueError, TypeError):
            result["total_reads_examined"] = _NA
    else:
        result["total_reads_examined"] = _NA

    # duplicate_reads_estimate = unpaired_dups + 2 * read_pair_dups
    if urd != _NA and rpd != _NA:
        try:
            result["duplicate_reads_estimate"] = str(
                int(urd) + 2 * int(rpd)
            )
        except (ValueError, TypeError):
            result["duplicate_reads_estimate"] = _NA
    else:
        result["duplicate_reads_estimate"] = _NA

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Parse Picard MarkDuplicates metrics for library complexity"
    )
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--metrics", required=True, help="Path to dup_metrics.txt")
    parser.add_argument("--output", required=True, help="Path to output TSV")
    args = parser.parse_args()

    if not os.path.isfile(args.metrics):
        print(f"ERROR: Metrics file not found: {args.metrics}", file=sys.stderr)
        sys.exit(1)

    if _is_fallback_metrics(args.metrics):
        fields = _empty_fields()
        source = "fallback"
    else:
        parsed = _parse_picard_metrics(args.metrics)
        if _has_parsed_metric(parsed):
            fields = _compute_derived(parsed)
            source = "picard_markduplicates"
        else:
            fields = _empty_fields()
            source = "unavailable"

    columns = [
        "sample",
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
    ]

    with open(args.output, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        row = [args.sample, source]
        for col in columns[2:]:
            row.append(fields.get(col, _NA))
        fh.write("\t".join(str(v) for v in row) + "\n")


if __name__ == "__main__":
    main()
