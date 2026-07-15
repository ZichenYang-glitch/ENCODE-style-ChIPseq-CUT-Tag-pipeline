#!/usr/bin/env python3
"""Convert a GTF annotation to a BED6 file of transcription start sites.

The script prefers transcript features. If no transcript records are present,
it falls back to gene features. Coordinates are converted from GTF's 1-based
closed intervals to BED's 0-based half-open intervals.
"""

import argparse
import sys


def parse_gtf_attributes(attr_text):
    """Return a dict of GTF attributes from column 9."""
    attrs = {}
    for item in attr_text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        parts = item.split(None, 1)
        if len(parts) != 2:
            continue
        key, value = parts
        attrs[key] = value.strip().strip('"')
    return attrs


def _bed_record(fields, index):
    """Return a BED6 TSS tuple from GTF fields, or None if invalid."""
    chrom, _, feature, start, end, _, strand, _, attrs_text = fields[:9]
    if strand not in ("+", "-"):
        return None
    try:
        start_i = int(start)
        end_i = int(end)
    except ValueError:
        return None
    if start_i < 1 or end_i < start_i:
        return None

    tss0 = start_i - 1 if strand == "+" else end_i - 1
    attrs = parse_gtf_attributes(attrs_text)
    name = attrs.get("transcript_id") or attrs.get("gene_id") or f"{feature}_{index}"
    return (chrom, tss0, tss0 + 1, name, "0", strand)


def collect_tss_records(gtf_path):
    """Collect transcript TSS records, falling back to gene records."""
    transcript_records = []
    gene_records = []
    with open(gtf_path) as fh:
        for index, line in enumerate(fh, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            if feature not in ("transcript", "gene"):
                continue
            record = _bed_record(fields, index)
            if record is None:
                continue
            if feature == "transcript":
                transcript_records.append(record)
            else:
                gene_records.append(record)

    records = transcript_records if transcript_records else gene_records
    return sorted(records, key=lambda r: (r[0], r[1], r[2], r[3], r[5]))


def write_bed(records, output_path):
    """Write BED6 records."""
    with open(output_path, "w") as out:
        for chrom, start, end, name, score, strand in records:
            out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Convert GTF transcript/gene records to a TSS BED6 file."
    )
    parser.add_argument("--gtf", required=True, help="Input GTF annotation")
    parser.add_argument("--output", required=True, help="Output BED6 path")
    args = parser.parse_args(argv)

    records = collect_tss_records(args.gtf)
    if not records:
        print(
            f"ERROR: no transcript or gene TSS records found in {args.gtf}",
            file=sys.stderr,
        )
        return 1
    write_bed(records, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
