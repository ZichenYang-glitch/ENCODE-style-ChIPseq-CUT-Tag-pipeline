"""Shared helpers for versioned QC writers.

Scripts are executed directly by Snakemake, so the checkout ``src`` directory
is added explicitly before importing the contract module.
"""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import sys

_REPO_ROOT = Path(__file__).resolve().parents[1]
_SRC = _REPO_ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from encode_pipeline.qc.contract import (  # noqa: E402
    ANALYSIS_TEMPLATE_METRICS_COLUMNS,
    ANALYSIS_TEMPLATE_METRICS_V1_COLUMNS,
    CROSS_CORRELATION_QC_V1_COLUMNS,
    CROSS_CORRELATION_QC_V2_COLUMNS,
    FRAGMENT_QC_V1_COLUMNS,
    FRAGMENT_QC_V2_COLUMNS,
    MITO_QC_V1_COLUMNS,
    MITO_QC_V2_COLUMNS,
    MNASE_QC_V1_COLUMNS,
    MNASE_QC_V2_COLUMNS,
    SAMPLE_QC_V1_COLUMNS,
    SAMPLE_QC_V2_COLUMNS,
    validate_cross_correlation_detail,
)

NA = "NA"


def read_cross_correlation_metadata(path: str) -> dict:
    """Read the exact version of the input-length/window metadata."""
    with open(path, encoding="utf-8") as handle:
        value = json.load(handle)
    if (
        not isinstance(value, dict)
        or type(value.get("schema_version")) is not int
        or value["schema_version"] != 1
    ):
        raise ValueError("cross-correlation metadata version is unsupported")
    return value


def read_cross_correlation_detail(path: str):
    version, row = read_versioned_tsv(
        path, CROSS_CORRELATION_QC_V1_COLUMNS, CROSS_CORRELATION_QC_V2_COLUMNS
    )
    if version == "2":
        validate_cross_correlation_detail(row)
    return version, row


def read_versioned_tsv(path: str, v1_columns, v2_columns):
    """Return ``(schema_version, row)`` for one exact v1 or v2 table.

    The function rejects duplicate columns, unknown headers, extra columns,
    missing files, empty tables and multi-row files.  It never guesses a
    version from context: a v2 header must also declare ``schema_version=2``
    in its single row, so an unknown version cannot be re-labelled as current.
    """
    if not path or not os.path.isfile(path):
        raise ValueError(f"QC source is missing: {path}")
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = tuple(next(reader))
        except StopIteration as exc:
            raise ValueError(f"QC source is empty: {path}") from exc
        if len(set(header)) != len(header):
            raise ValueError(f"QC source has duplicate columns: {path}")
        if header == tuple(v1_columns):
            version = "1"
            expected = tuple(v1_columns)
        elif header == tuple(v2_columns):
            version = "2"
            expected = tuple(v2_columns)
        else:
            raise ValueError(f"QC source header is invalid: {path}")
        rows = [row for row in reader if row]
    if len(rows) != 1:
        raise ValueError(f"QC source must contain exactly one row: {path}")
    if len(rows[0]) != len(expected):
        raise ValueError(f"QC source row width is invalid: {path}")
    row = dict(zip(expected, rows[0], strict=True))
    if version == "2" and row.get("schema_version") != "2":
        raise ValueError(f"QC source schema version is unsupported: {path}")
    return version, row


def write_exact_tsv(path: str, columns, row):
    """Write one exact fixed-header row with LF endings."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerow({column: row.get(column, NA) for column in columns})


def format_fraction(numerator: int, denominator: int) -> str:
    """Return a six-decimal fraction, or NA when undefined."""
    if denominator <= 0:
        return NA
    return f"{numerator / denominator:.6f}"


def load_stage_metrics(path: str) -> dict:
    """Load and validate the template-filter stage metrics JSON."""
    if not path or not os.path.isfile(path):
        raise ValueError(f"stage metrics are missing: {path}")
    with open(path, encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError("stage metrics must be a JSON object")
    counts = {}
    for key in (
        "total_alignments",
        "kept_alignments",
        "dropped_alignments",
        "dropped_templates",
    ):
        value = data.get(key)
        if type(value) is not int or value < 0:
            raise ValueError(f"stage metric {key} is invalid")
        counts[key] = value
    for key in ("dropped_alignments_by_reason", "dropped_templates_by_reason"):
        value = data.get(key)
        if not isinstance(value, dict):
            raise ValueError(f"stage metric {key} is invalid")
        if any(type(item) is not int or item < 0 for item in value.values()):
            raise ValueError(f"stage metric {key} is invalid")
    if (
        counts["kept_alignments"] + counts["dropped_alignments"]
        != counts["total_alignments"]
    ):
        raise ValueError("stage metrics do not conserve alignments")
    return data


__all__ = [
    "ANALYSIS_TEMPLATE_METRICS_COLUMNS",
    "ANALYSIS_TEMPLATE_METRICS_V1_COLUMNS",
    "CROSS_CORRELATION_QC_V1_COLUMNS",
    "CROSS_CORRELATION_QC_V2_COLUMNS",
    "FRAGMENT_QC_V1_COLUMNS",
    "FRAGMENT_QC_V2_COLUMNS",
    "MITO_QC_V1_COLUMNS",
    "MITO_QC_V2_COLUMNS",
    "MNASE_QC_V1_COLUMNS",
    "MNASE_QC_V2_COLUMNS",
    "NA",
    "SAMPLE_QC_V1_COLUMNS",
    "SAMPLE_QC_V2_COLUMNS",
    "format_fraction",
    "load_stage_metrics",
    "read_versioned_tsv",
    "write_exact_tsv",
]
