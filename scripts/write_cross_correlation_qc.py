#!/usr/bin/env python3
"""Write the versioned per-sample cross-correlation detail table.

The raw ``.cc.qc`` file remains the phantompeakqualtools output.  This writer
combines it with the tagAlign/window metadata and emits the exact v2 detail
contract consumed by the sample QC summary, policy/reporting, and MultiQC.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys

_REPO_ROOT = Path(__file__).resolve().parents[1]
for _path in (_REPO_ROOT / "scripts", _REPO_ROOT / "src"):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from encode_pipeline.qc.contract import (  # noqa: E402
    CROSS_CORRELATION_QC_SCHEMA_VERSION,
    CROSS_CORRELATION_QC_V2_COLUMNS,
    CROSS_CORRELATION_THRESHOLD_HINTS,
    validate_cross_correlation_detail,
)

from parse_cross_correlation import (  # noqa: E402
    interpret_cross_correlation,
    parse_cc_qc_record,
)
from qc_common import read_cross_correlation_metadata, write_exact_tsv  # noqa: E402

NA = "NA"


def _fmt(value):
    if value is None:
        return NA
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def _read_metadata(path: str | None) -> dict:
    if not path:
        return {}
    try:
        value = read_cross_correlation_metadata(path)
    except (OSError, ValueError) as exc:
        return {"status": "parse_failed", "error": str(exc)}
    return value if isinstance(value, dict) else {"status": "parse_failed"}


def _execution_status(parsed: dict, metadata: dict, cc_qc: str, override: str | None):
    if override:
        return override
    if metadata.get("status") == "no_data":
        return "no_data"
    if metadata.get("status") == "parse_failed":
        return "parse_failed"
    if not cc_qc or not os.path.isfile(cc_qc):
        return "missing"
    if (
        parsed.get("nsc") is None
        and parsed.get("rsc") is None
        and parsed.get("candidate_status") != "ok"
    ):
        return "parse_failed"
    return "reported"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--assay", required=True)
    parser.add_argument("--layout", required=True)
    parser.add_argument("--cc-qc", default="")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--execution-status", default=None)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    metadata = _read_metadata(args.metadata)
    parsed = parse_cc_qc_record(args.cc_qc)
    execution = _execution_status(parsed, metadata, args.cc_qc, args.execution_status)
    interpreted = interpret_cross_correlation(
        parsed, metadata=metadata, execution_status=execution
    )
    if execution != "reported":
        for key in (
            "estimated_fragment_length",
            "phantom_peak",
            "nsc",
            "rsc",
            "quality_tag",
            "candidate_lengths",
            "candidate_correlations",
            "primary_correlation",
            "phantom_correlation",
            "phantom_delta",
        ):
            interpreted[key] = None
        interpreted["candidate_status"] = "missing"
        interpreted["threshold_hint"] = "not_evaluated"

    candidate_lengths = interpreted.get("candidate_lengths") or []
    candidate_correlations = interpreted.get("candidate_correlations") or []
    row = {
        "schema_version": CROSS_CORRELATION_QC_SCHEMA_VERSION,
        "sample": args.sample,
        "cc_qc_file": os.path.basename(args.cc_qc) if args.cc_qc else NA,
        "estimated_fragment_length": _fmt(interpreted.get("estimated_fragment_length")),
        "phantom_peak": _fmt(interpreted.get("phantom_peak")),
        "nsc": _fmt(interpreted.get("nsc")),
        "rsc": _fmt(interpreted.get("rsc")),
        "quality_tag": interpreted.get("quality_tag") or NA,
        "assay": args.assay,
        "layout": args.layout,
        "input_stage": metadata.get("input_stage")
        or ("analysis_r1" if args.layout == "PE" else "analysis_se"),
        "execution_status": execution,
        "interpretation_status": interpreted["interpretation_status"],
        "status_flags": interpreted["status_flags"],
        "threshold_hint": interpreted["threshold_hint"],
        "read_length": _fmt(metadata.get("read_length")),
        "read_length_method": metadata.get("read_length_method") or NA,
        "variable_read_length": "true"
        if metadata.get("variable_read_length")
        else "false",
        "target_type": metadata.get("target_type") or "unknown",
        "target_type_source": metadata.get("target_type_source") or NA,
        "target_inferred_histone_class": metadata.get("target_inferred_histone_class")
        or "unknown",
        "window_lower": _fmt(metadata.get("window_lower")),
        "window_upper": _fmt(metadata.get("window_upper")),
        "candidate_count": str(len(candidate_lengths)),
        "candidate_lengths": ",".join(_fmt(value) for value in candidate_lengths) or NA,
        "candidate_correlations": ",".join(
            _fmt(value) for value in candidate_correlations
        )
        or NA,
        "phantom_delta": _fmt(interpreted.get("phantom_delta")),
        "phantom_position_selected": "true"
        if interpreted.get("phantom_position_selected")
        else "false",
        "phantom_competes": "true" if interpreted.get("phantom_competes") else "false",
    }
    if row["threshold_hint"] not in CROSS_CORRELATION_THRESHOLD_HINTS:
        raise ValueError("cross-correlation threshold hint is invalid")
    validate_cross_correlation_detail(row)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    write_exact_tsv(args.output, CROSS_CORRELATION_QC_V2_COLUMNS, row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
