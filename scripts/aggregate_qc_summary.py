#!/usr/bin/env python3
"""Aggregate per-sample QC summary TSVs into a project-level summary.

Version 1 inputs stay version 1 when every input is version 1.  If any input
is version 2, the project summary is written with the exact v2 header and v1
rows are padded with NA for fields that version did not define.  A v2 header
must declare ``schema_version=2`` in every row; unknown versions are rejected
rather than re-labelled.  Arbitrary extra columns are rejected.
"""

from __future__ import annotations

import argparse
import csv
import os
import sys

from qc_common import (
    NA,
    SAMPLE_QC_V1_COLUMNS,
    SAMPLE_QC_V2_COLUMNS,
)

_QC_SUMMARY_COLUMNS = SAMPLE_QC_V2_COLUMNS


def _read_table(filepath):
    if not os.path.isfile(filepath):
        print(f"ERROR: input file not found: {filepath}", file=sys.stderr)
        sys.exit(1)
    with open(filepath, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = tuple(next(reader))
        except StopIteration:
            print(f"ERROR: empty or unreadable file: {filepath}", file=sys.stderr)
            sys.exit(1)
        if len(set(header)) != len(header):
            print(f"ERROR: duplicate columns in {filepath}", file=sys.stderr)
            sys.exit(1)
        if header == SAMPLE_QC_V1_COLUMNS:
            version = "1"
        elif header == SAMPLE_QC_V2_COLUMNS:
            version = "2"
        else:
            print(f"ERROR: header mismatch in {filepath}", file=sys.stderr)
            print(
                f"  Expected {len(SAMPLE_QC_V1_COLUMNS)} v1 or "
                f"{len(SAMPLE_QC_V2_COLUMNS)} v2 columns",
                file=sys.stderr,
            )
            print(f"  Got      {len(header)} columns", file=sys.stderr)
            sys.exit(1)
        rows = []
        for row in reader:
            if not row:
                continue
            if len(row) != len(header):
                print(f"ERROR: row width mismatch in {filepath}", file=sys.stderr)
                sys.exit(1)
            parsed = dict(zip(header, row, strict=True))
            if version == "2" and parsed.get("schema_version") != "2":
                print(
                    f"ERROR: unsupported schema_version in {filepath}",
                    file=sys.stderr,
                )
                print("  A v2 header must declare schema_version=2", file=sys.stderr)
                sys.exit(1)
            rows.append(parsed)
    return version, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, help="Output TSV path")
    parser.add_argument("inputs", nargs="*", help="Per-sample qc_summary.tsv files")
    args = parser.parse_args()

    tables = []
    saw_v2 = False
    for filepath in args.inputs:
        version, rows = _read_table(filepath)
        saw_v2 = saw_v2 or version == "2"
        tables.append((version, rows))

    columns = SAMPLE_QC_V2_COLUMNS if (saw_v2 or not tables) else SAMPLE_QC_V1_COLUMNS
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for version, rows in tables:
            for row in rows:
                if version == "1" and columns is SAMPLE_QC_V2_COLUMNS:
                    row = {**{column: NA for column in columns}, **row}
                writer.writerow(row)


if __name__ == "__main__":
    main()
