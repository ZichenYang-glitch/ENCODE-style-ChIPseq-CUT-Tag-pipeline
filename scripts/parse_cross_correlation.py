#!/usr/bin/env python3
"""Parse phantompeakqualtools .cc.qc files into a project-level summary TSV.

Input:  one or more .cc.qc files (phantompeakqualtools / run_spp.R -out format)
Output: TSV with columns: sample, cc_qc_file, estimated_fragment_length,
        phantom_peak, nsc, rsc, quality_flag
"""

import argparse
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
    """Parse a string to float, returning None on failure."""
    if raw is None:
        return None
    try:
        return float(raw)
    except (ValueError, TypeError):
        return None


def parse_cc_qc_file(filepath):
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
            for data_line in lines[i + 1:]:
                if data_line.startswith("#"):
                    continue
                dfields = data_line.split("\t")
                if len(dfields) < len(fields):
                    continue

                def _get(name):
                    idx = col_map.get(name)
                    if idx is None or idx >= len(dfields):
                        return None
                    return dfields[idx].strip()

                result["estimated_fragment_length"] = _try_float(_get("estFragLen"))
                result["phantom_peak"] = _try_float(_get("PhantomPeak"))
                result["nsc"] = _try_float(_get("NSC"))
                result["rsc"] = _try_float(_get("RSC"))
                result["quality_tag"] = _get("QualityTag") or None
                return result

    # --- Headerless ENCODE-style parsing -----------------------------------
    # Common run_spp.R -out output is one tab-separated data row:
    # Filename, numReads, estFragLen, corr_estFragLen, PhantomPeak,
    # corr_phantomPeak, argmin_corr, min_corr, NSC, RSC, QualityTag
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
        result["estimated_fragment_length"] = _try_float(fields[2].strip())
        result["phantom_peak"] = _try_float(fields[4].strip())
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


def sample_name_from_path(filepath):
    """Derive a sample name from a .cc.qc filename.

    Strips the .cc.qc suffix and returns the basename.
    Example: /path/to/s1.cc.qc → s1
    """
    basename = os.path.basename(filepath)
    if basename.endswith(".cc.qc"):
        return basename[:-len(".cc.qc")]
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
        "--input", nargs="+", required=True,
        help="One or more .cc.qc file paths",
    )
    parser.add_argument(
        "--output", required=True,
        help="Output TSV path",
    )
    args = parser.parse_args()

    rows = []
    for filepath in args.input:
        sample = sample_name_from_path(filepath)
        parsed = parse_cc_qc_file(filepath)
        flag = _quality_flag(parsed["nsc"], parsed["rsc"])
        rows.append({
            "sample": sample,
            "cc_qc_file": filepath,
            "estimated_fragment_length": _fmt(parsed["estimated_fragment_length"]),
            "phantom_peak": _fmt(parsed["phantom_peak"]),
            "nsc": _fmt(parsed["nsc"]),
            "rsc": _fmt(parsed["rsc"]),
            "quality_flag": flag,
        })

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)

    with open(args.output, "w") as fh:
        header = [
            "sample", "cc_qc_file", "estimated_fragment_length",
            "phantom_peak", "nsc", "rsc", "quality_flag",
        ]
        fh.write("\t".join(header) + "\n")
        for row in rows:
            fh.write("\t".join(row[h] for h in header) + "\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
