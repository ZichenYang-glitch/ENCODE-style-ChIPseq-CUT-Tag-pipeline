#!/usr/bin/env python3
"""Summarize per-sample cross-correlation detail rows for MultiQC."""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import sys

_REPO_ROOT = Path(__file__).resolve().parents[1]
for _path in (_REPO_ROOT / "scripts", _REPO_ROOT / "src"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))


from qc_common import NA, read_cross_correlation_detail  # noqa: E402

SUMMARY_COLUMNS = (
    "sample",
    "cc_qc_file",
    "estimated_fragment_length",
    "candidate_count",
    "phantom_peak",
    "phantom_delta",
    "nsc",
    "rsc",
    "execution_status",
    "interpretation_status",
    "threshold_hint",
    "target_type",
    "read_length",
    "input_stage",
    "window_lower",
    "window_upper",
    "variable_read_length",
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", nargs="+", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows = []
    for path in args.input:
        _version, detail = read_cross_correlation_detail(path)
        row = {column: detail.get(column, NA) for column in SUMMARY_COLUMNS}
        row["sample"] = detail.get("sample", os.path.basename(path))
        rows.append(row)

    if not rows:
        print("ERROR: no cross-correlation detail inputs", file=sys.stderr)
        return 1
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=SUMMARY_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
