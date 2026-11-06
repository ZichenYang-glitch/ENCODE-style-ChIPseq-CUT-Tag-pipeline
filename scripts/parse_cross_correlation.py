#!/usr/bin/env python3
"""Parse phantompeakqualtools .cc.qc files into a project-level summary TSV.

Input:  one or more .cc.qc files (phantompeakqualtools / run_spp.R -out format)
Output: TSV with columns: sample, cc_qc_file, estimated_fragment_length,
        phantom_peak, nsc, rsc, quality_flag
"""

import argparse
import math
import os
import sys

# Quality thresholds (ENCODE-inspired)
NSC_THRESHOLD = 1.05
RSC_THRESHOLD = 0.8


def _quality_flag(nsc, rsc):
    """Return a quality flag string based on NSC and RSC thresholds.

    Returns one of: ok, low_nsc, low_rsc, low_nsc_low_rsc.
    If either value is None (unparseable), returns "parse_failed".
    """
    if nsc is None or rsc is None:
        return "parse_failed"
    low_nsc = nsc < NSC_THRESHOLD
    low_rsc = rsc < RSC_THRESHOLD
    if low_nsc and low_rsc:
        return "low_nsc_low_rsc"
    if low_nsc:
        return "low_nsc"
    if low_rsc:
        return "low_rsc"
    return "ok"


def _fmt(val):
    """Format a parsed value for TSV output.

    Returns the string representation, or "NA" for None.
    """
    if val is None:
        return "NA"
    return str(val)


def _try_float(raw):
    """Parse a string to a finite float, returning None otherwise."""
    if raw is None:
        return None
    try:
        value = float(raw)
    except (ValueError, TypeError):
        return None
    return value if math.isfinite(value) else None


def _primary_fragment_length(raw):
    """Return the first tool-ranked candidate, or None for an invalid list.

    Every token must be a finite number. Empty/NA/malformed tokens invalidate
    the field; never skip a bad first candidate to select a later one.
    """
    if raw is None:
        return None
    tokens = [token.strip() for token in str(raw).split(",")]
    if not tokens or any(_try_float(token) is None for token in tokens):
        return None
    return _try_float(tokens[0])


def parse_candidates(raw_lengths, raw_correlations):
    """Parse paired fragment-length candidates and their correlation scores.

    Returns ``(first_length, pairs, status)`` where ``pairs`` is a tuple of
    ``(length, correlation)`` tuples in the tool's ranked order.  The status is
    ``ok``, ``missing``, or ``invalid``.  Lengths and scores must have the same
    length and every token must be finite; ties are preserved in order.
    """
    if raw_lengths is None or raw_correlations is None:
        return None, (), "missing"
    lengths_raw = [token.strip() for token in str(raw_lengths).split(",")]
    scores_raw = [token.strip() for token in str(raw_correlations).split(",")]
    if not lengths_raw or not scores_raw or lengths_raw == [""] or scores_raw == [""]:
        return None, (), "missing"
    if len(lengths_raw) != len(scores_raw):
        return None, (), "invalid"
    lengths = [_try_float(token) for token in lengths_raw]
    scores = [_try_float(token) for token in scores_raw]
    if any(value is None for value in lengths + scores):
        return None, (), "invalid"
    pairs = tuple(zip(lengths, scores, strict=True))
    return pairs[0][0], pairs, "ok"


def parse_cc_qc_record(filepath):
    """Parse a phantompeakqualtools .cc.qc file.

    Handles ENCODE-style tab-separated formats produced by run_spp.R with the
    -out flag. Some versions emit a header followed by a data row:

        Filename  numReads  estFragLen  corr_estFragLen  PhantomPeak
        corr_phantomPeak  argmin_corr  min_corr  NSC  RSC  QualityTag

    Other common outputs are a single data row with the same column order and
    no header.

    Returns a dict with keys: estimated_fragment_length, phantom_peak,
    nsc, rsc, quality_tag.  Unparseable fields are None.
    """
    result = {
        "estimated_fragment_length": None,
        "phantom_peak": None,
        "nsc": None,
        "rsc": None,
        "quality_tag": None,
        "candidate_lengths": [],
        "candidate_correlations": [],
        "candidate_status": "missing",
        "primary_correlation": None,
        "phantom_correlation": None,
        "argmin_corr": None,
        "min_corr": None,
        "candidate_lengths_raw": None,
        "candidate_correlations_raw": None,
    }

    try:
        with open(filepath, "r") as fh:
            lines = [line.rstrip("\n").rstrip("\r") for line in fh if line.strip()]
    except OSError:
        return result

    if not lines:
        return result

    # Strategy: scan all lines. The data table is tab-separated and starts
    # with "Filename" (header) followed by the data row.  Comment lines
    # (starting with #) may contain human-readable key=value pairs from
    # older run_spp.R versions; try those as a fallback.

    # --- Header-based parsing ------------------------------------------------
    for i, line in enumerate(lines):
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        # The standard header
        if "Filename" in fields and "NSC" in fields and "RSC" in fields:
            # Map column names to indices
            col_map = {name.strip(): idx for idx, name in enumerate(fields)}
            # Data row is the next line (or we check subsequent lines)
            data_rows = [
                data_line
                for data_line in lines[i + 1 :]
                if not data_line.startswith("#")
                and len(data_line.split("\t")) >= len(fields)
            ]
            if len(data_rows) > 1:
                return result
            if not data_rows:
                continue
            data_line = data_rows[0]
            dfields = data_line.split("\t")

            def _get(name):
                idx = col_map.get(name)
                if idx is None or idx >= len(dfields):
                    return None
                return dfields[idx].strip()

            raw_lengths = _get("estFragLen")
            raw_scores = _get("corr_estFragLen")
            result["candidate_lengths_raw"] = raw_lengths
            result["candidate_correlations_raw"] = raw_scores
            first, pairs, candidate_status = parse_candidates(raw_lengths, raw_scores)
            result["estimated_fragment_length"] = first
            result["candidate_lengths"] = [pair[0] for pair in pairs]
            result["candidate_correlations"] = [pair[1] for pair in pairs]
            result["candidate_status"] = candidate_status
            result["primary_correlation"] = pairs[0][1] if pairs else None
            result["phantom_peak"] = _try_float(_get("PhantomPeak"))
            result["phantom_correlation"] = _try_float(_get("corr_phantomPeak"))
            result["argmin_corr"] = _try_float(_get("argmin_corr"))
            result["min_corr"] = _try_float(_get("min_corr"))
            result["nsc"] = _try_float(_get("NSC"))
            result["rsc"] = _try_float(_get("RSC"))
            result["quality_tag"] = _get("QualityTag") or None
            return result

    # --- Headerless ENCODE-style parsing -----------------------------------
    # Common run_spp.R -out output is one tab-separated data row:
    # Filename, numReads, estFragLen, corr_estFragLen, PhantomPeak,
    # corr_phantomPeak, argmin_corr, min_corr, NSC, RSC, QualityTag
    valid_rows = []
    for line in lines:
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 10:
            continue
        # Skip obvious headers and malformed rows.  numReads should be numeric.
        if fields[0].strip().lower() == "filename":
            continue
        if _try_float(fields[1].strip()) is None:
            continue
        valid_rows.append(fields)
    if len(valid_rows) > 1:
        return result
    if valid_rows:
        fields = valid_rows[0]
        raw_lengths = fields[2].strip()
        raw_scores = fields[3].strip()
        result["candidate_lengths_raw"] = raw_lengths
        result["candidate_correlations_raw"] = raw_scores
        first, pairs, candidate_status = parse_candidates(raw_lengths, raw_scores)
        result["estimated_fragment_length"] = first
        result["candidate_lengths"] = [pair[0] for pair in pairs]
        result["candidate_correlations"] = [pair[1] for pair in pairs]
        result["candidate_status"] = candidate_status
        result["primary_correlation"] = pairs[0][1] if pairs else None
        result["phantom_peak"] = _try_float(fields[4].strip())
        result["phantom_correlation"] = _try_float(fields[5].strip())
        result["argmin_corr"] = _try_float(fields[6].strip())
        result["min_corr"] = _try_float(fields[7].strip())
        result["nsc"] = _try_float(fields[8].strip())
        result["rsc"] = _try_float(fields[9].strip())
        result["quality_tag"] = fields[10].strip() if len(fields) > 10 else None
        return result

    # --- Comment-line fallback (older run_spp.R format) ---------------------
    # Example:
    # # cross-correlation peak = 175
    # # phantom peak = 225
    # # NSC = 1.08
    # # RSC = 1.45
    for line in lines:
        if not line.startswith("#"):
            continue
        stripped = line.lstrip("#").strip()
        if "=" in stripped:
            key, _, val = stripped.partition("=")
            key = key.strip().lower()
            val = val.strip()
            if "cross-correlation" in key and "peak" in key:
                result["estimated_fragment_length"] = _try_float(val)
            elif "phantom peak" in key:
                result["phantom_peak"] = _try_float(val)
            elif key == "nsc":
                result["nsc"] = _try_float(val)
            elif key == "rsc":
                result["rsc"] = _try_float(val)

    return result


def parse_cc_qc_file(filepath):
    """Compatibility wrapper returning the original five scalar fields."""
    record = parse_cc_qc_record(filepath)
    raw_lengths = record.get("candidate_lengths_raw")
    estimated = (
        _primary_fragment_length(raw_lengths)
        if raw_lengths is not None
        else record["estimated_fragment_length"]
    )
    return {
        "estimated_fragment_length": estimated,
        "phantom_peak": record["phantom_peak"],
        "nsc": record["nsc"],
        "rsc": record["rsc"],
        "quality_tag": record["quality_tag"],
    }


def sample_name_from_path(filepath):
    """Derive a sample name from a .cc.qc filename.

    Strips the .cc.qc suffix and returns the basename.
    Example: /path/to/s1.cc.qc → s1
    """
    basename = os.path.basename(filepath)
    if basename.endswith(".cc.qc"):
        return basename[: -len(".cc.qc")]
    # Fallback: strip the last two suffixes (e.g. .cc.qc from full path)
    name = basename
    while name.endswith(".qc") or name.endswith(".cc"):
        name = os.path.splitext(name)[0]
    return name


def main():
    parser = argparse.ArgumentParser(
        description="Parse phantompeakqualtools .cc.qc files into a project-level TSV"
    )
    parser.add_argument(
        "--input",
        nargs="+",
        required=True,
        help="One or more .cc.qc file paths",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV path",
    )
    args = parser.parse_args()

    rows = []
    for filepath in args.input:
        sample = sample_name_from_path(filepath)
        parsed = parse_cc_qc_file(filepath)
        flag = _quality_flag(parsed["nsc"], parsed["rsc"])
        rows.append(
            {
                "sample": sample,
                "cc_qc_file": os.path.basename(filepath),
                "estimated_fragment_length": _fmt(parsed["estimated_fragment_length"]),
                "phantom_peak": _fmt(parsed["phantom_peak"]),
                "nsc": _fmt(parsed["nsc"]),
                "rsc": _fmt(parsed["rsc"]),
                "quality_flag": flag,
            }
        )

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    with open(args.output, "w") as fh:
        header = [
            "sample",
            "cc_qc_file",
            "estimated_fragment_length",
            "phantom_peak",
            "nsc",
            "rsc",
            "quality_flag",
        ]
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(row[h] for h in header) + "\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())


# Local diagnostic thresholds only.  These are not ENCODE certification.
LOCAL_NSC_THRESHOLD = NSC_THRESHOLD
LOCAL_RSC_THRESHOLD = RSC_THRESHOLD
PHANTOM_TOLERANCE_MIN = 10
PHANTOM_TOLERANCE_FRACTION = 0.1


def interpret_cross_correlation(parsed, *, metadata, execution_status=None):
    """Classify one parsed .cc.qc record for the versioned detail contract.

    ``execution_status`` describes whether the tool ran and produced a usable
    record.  ``interpretation_status`` describes what the numbers mean for the
    declared input/target.  They are intentionally separate: a tool failure or
    malformed file must never be reported as a biological failure.
    """
    metadata = metadata or {}
    assay = str(metadata.get("assay") or "unknown")
    target_type = str(metadata.get("target_type") or "unknown")
    read_length = metadata.get("read_length")
    variable = bool(metadata.get("variable_read_length", False))
    candidate_status = parsed.get("candidate_status", "missing")
    nsc = parsed.get("nsc")
    rsc = parsed.get("rsc")
    first = parsed.get("estimated_fragment_length")
    phantom = parsed.get("phantom_peak")

    if execution_status is None:
        if nsc is None and rsc is None and candidate_status != "ok":
            execution_status = "parse_failed"
        else:
            execution_status = "reported"

    phantom_delta = None
    phantom_position_selected = False
    phantom_competes = False
    if first is not None and phantom is not None and read_length:
        phantom_delta = abs(first - phantom)
        tolerance = (
            max(
                PHANTOM_TOLERANCE_MIN,
                int(round(read_length * PHANTOM_TOLERANCE_FRACTION)),
            )
            if isinstance(read_length, (int, float)) and read_length > 0
            else PHANTOM_TOLERANCE_MIN
        )
        candidates = parsed.get("candidate_lengths", [])
        correlations = parsed.get("candidate_correlations", [])
        primary_corr = parsed.get("primary_correlation")
        phantom_corr = parsed.get("phantom_correlation")
        # run_spp.R prints the phantom score at lower precision than the
        # candidate scores. Compare at the tool's printed precision.
        tied_positions = [
            length
            for length, score in zip(candidates, correlations, strict=True)
            if primary_corr is not None
            and math.isclose(score, primary_corr, rel_tol=1e-6, abs_tol=1e-12)
        ]
        phantom_position_selected = abs(phantom - read_length) <= tolerance and any(
            abs(length - phantom) <= tolerance for length in tied_positions
        )
        phantom_competes = (
            primary_corr is not None
            and phantom_corr is not None
            and (
                primary_corr <= phantom_corr
                or math.isclose(primary_corr, phantom_corr, rel_tol=1e-6, abs_tol=1e-12)
            )
        )
    phantom_selected = phantom_position_selected and phantom_competes

    flags = []
    if execution_status != "reported":
        flags.append(execution_status)
        return {
            **parsed,
            "execution_status": execution_status,
            "interpretation_status": execution_status,
            "status_flags": ",".join(flags),
            "threshold_hint": "not_evaluated",
            "phantom_delta": phantom_delta,
            "phantom_position_selected": phantom_position_selected,
            "phantom_competes": phantom_competes,
            "phantom_selected": phantom_selected,
        }

    if variable:
        flags.append("variable_read_length")
    if phantom_selected:
        flags.append("phantom_selected")
    if assay != "chipseq":
        flags.append("non_chip_diagnostic")
    elif target_type not in ("tf", "histone"):
        flags.append("unknown_target")
    if candidate_status != "ok" or nsc is None or rsc is None:
        flags.append("missing_values")
    if nsc is not None and nsc < LOCAL_NSC_THRESHOLD:
        flags.append("low_nsc")
    if rsc is not None and rsc < LOCAL_RSC_THRESHOLD:
        flags.append("low_rsc")

    if "variable_read_length" in flags:
        interpretation = "variable_read_length"
    elif "phantom_selected" in flags:
        interpretation = "phantom_selected"
    elif "non_chip_diagnostic" in flags:
        interpretation = "non_chip_diagnostic"
    elif "unknown_target" in flags:
        interpretation = "unknown_target"
    elif "missing_values" in flags:
        interpretation = "missing_values"
    elif "low_nsc" in flags and "low_rsc" in flags:
        interpretation = "low_nsc_low_rsc"
    elif "low_nsc" in flags:
        interpretation = "low_nsc"
    elif "low_rsc" in flags:
        interpretation = "low_rsc"
    else:
        interpretation = "ok"

    definitive_chip = (
        assay == "chipseq"
        and target_type in ("tf", "histone")
        and not variable
        and not phantom_selected
        and candidate_status == "ok"
        and nsc is not None
        and rsc is not None
    )
    if definitive_chip:
        threshold_hint = (
            "local_pass"
            if nsc >= LOCAL_NSC_THRESHOLD and rsc >= LOCAL_RSC_THRESHOLD
            else "local_fail"
        )
    else:
        threshold_hint = "not_evaluated"

    return {
        **parsed,
        "execution_status": execution_status,
        "interpretation_status": interpretation,
        "status_flags": ",".join(flags) if flags else "ok",
        "threshold_hint": threshold_hint,
        "phantom_delta": phantom_delta,
        "phantom_position_selected": phantom_position_selected,
        "phantom_competes": phantom_competes,
        "phantom_selected": phantom_selected,
    }
