#!/usr/bin/env python3
"""Unit tests for scripts/parse_cross_correlation.py."""

import os
import sys
import tempfile
import shutil
import subprocess

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_REPO, "scripts"))
from parse_cross_correlation import (
    parse_cc_qc_file,
    sample_name_from_path,
    _quality_flag,
    _try_float,
)


def _write_file(path, lines):
    with open(path, "w") as fh:
        for line in lines:
            fh.write(line + "\n")


# --- _try_float ----------------------------------------------------------

def test_try_float_valid():
    assert _try_float("1.05") == 1.05
    assert _try_float("0.8") == 0.8
    assert _try_float("-1") == -1.0


def test_try_float_none():
    assert _try_float(None) is None


def test_try_float_invalid():
    assert _try_float("nope") is None
    assert _try_float("") is None


# --- _quality_flag -------------------------------------------------------

def test_quality_flag_ok():
    assert _quality_flag(1.1, 1.0) == "ok"


def test_quality_flag_low_nsc():
    assert _quality_flag(1.0, 1.0) == "low_nsc"


def test_quality_flag_low_rsc():
    assert _quality_flag(1.1, 0.5) == "low_rsc"


def test_quality_flag_low_both():
    assert _quality_flag(0.5, 0.5) == "low_nsc_low_rsc"


def test_quality_flag_parse_failed():
    assert _quality_flag(None, 1.0) == "parse_failed"
    assert _quality_flag(1.0, None) == "parse_failed"
    assert _quality_flag(None, None) == "parse_failed"


# --- sample_name_from_path -----------------------------------------------

def test_sample_name_standard():
    assert sample_name_from_path("/path/to/s1.cc.qc") == "s1"


def test_sample_name_basename_only():
    assert sample_name_from_path("s1.cc.qc") == "s1"


def test_sample_name_no_suffix_fallback():
    # Falls back to stripping known suffixes or returning basename
    name = sample_name_from_path("/path/to/s1.other")
    # Should return the basename unchanged
    assert name == "s1.other"


# --- parse_cc_qc_file ----------------------------------------------------

def test_parse_standard_header_format():
    """Tab-separated header + data format (standard run_spp.R -out)."""
    td = tempfile.mkdtemp(prefix="cc_parse_")
    try:
        fp = os.path.join(td, "s1.cc.qc")
        _write_file(fp, [
            "Filename\tnumReads\testFragLen\tcorr_estFragLen\tPhantomPeak\t"
            "corr_phantomPeak\targmin_corr\tmin_corr\tNSC\tRSC\tQualityTag",
            "s1.bam\t10000000\t175\t0.85\t225\t0.60\t7\t0.55\t1.12\t1.45\t-1",
        ])
        result = parse_cc_qc_file(fp)
        assert result["estimated_fragment_length"] == 175
        assert result["phantom_peak"] == 225
        assert result["nsc"] == 1.12
        assert result["rsc"] == 1.45
        assert result["quality_tag"] == "-1"
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_parse_headerless_encode_format():
    """Single-line ENCODE-style run_spp.R output without a header."""
    td = tempfile.mkdtemp(prefix="cc_parse_")
    try:
        fp = os.path.join(td, "s1.cc.qc")
        _write_file(fp, [
            "s1.tagAlign\t10000000\t175\t0.85\t225\t0.60\t7\t0.55\t1.12\t1.45\t-1",
        ])
        result = parse_cc_qc_file(fp)
        assert result["estimated_fragment_length"] == 175
        assert result["phantom_peak"] == 225
        assert result["nsc"] == 1.12
        assert result["rsc"] == 1.45
        assert result["quality_tag"] == "-1"
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_parse_comment_header_format():
    """Older comment-line format (# key = value)."""
    td = tempfile.mkdtemp(prefix="cc_parse_")
    try:
        fp = os.path.join(td, "s2.cc.qc")
        _write_file(fp, [
            "# Filename: s2.bam",
            "# cross-correlation peak = 200",
            "# phantom peak = 250",
            "# NSC = 2.30",
            "# RSC = 1.80",
        ])
        result = parse_cc_qc_file(fp)
        assert result["estimated_fragment_length"] == 200
        assert result["phantom_peak"] == 250
        assert result["nsc"] == 2.30
        assert result["rsc"] == 1.80
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_parse_malformed_returns_none():
    """Unparseable file → all fields None."""
    td = tempfile.mkdtemp(prefix="cc_parse_")
    try:
        fp = os.path.join(td, "bad.cc.qc")
        _write_file(fp, [
            "garbage line with no structure",
            "also bad",
        ])
        result = parse_cc_qc_file(fp)
        assert result["nsc"] is None
        assert result["rsc"] is None
        assert result["estimated_fragment_length"] is None
        assert result["phantom_peak"] is None
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_parse_empty_file():
    """Empty file → all None."""
    td = tempfile.mkdtemp(prefix="cc_parse_")
    try:
        fp = os.path.join(td, "empty.cc.qc")
        _write_file(fp, [])
        result = parse_cc_qc_file(fp)
        assert result["nsc"] is None
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_parse_nonexistent_file():
    """Missing file → all fields None (no exception)."""
    result = parse_cc_qc_file("/nonexistent/path.cc.qc")
    assert result["nsc"] is None
    assert result["rsc"] is None


def test_parse_header_with_columns_only():
    """Header-only file (no data rows) → all None."""
    td = tempfile.mkdtemp(prefix="cc_parse_")
    try:
        fp = os.path.join(td, "header_only.cc.qc")
        _write_file(fp, [
            "Filename\tnumReads\testFragLen\tNSC\tRSC\tPhantomPeak\tQualityTag",
        ])
        result = parse_cc_qc_file(fp)
        assert result["nsc"] is None
        assert result["rsc"] is None
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_parse_partial_data():
    """Some columns missing → missing fields None, present fields parsed."""
    td = tempfile.mkdtemp(prefix="cc_parse_")
    try:
        fp = os.path.join(td, "partial.cc.qc")
        _write_file(fp, [
            "Filename\tNSC\tRSC",
            "s1.bam\t1.08\t1.33",
        ])
        result = parse_cc_qc_file(fp)
        assert result["nsc"] == 1.08
        assert result["rsc"] == 1.33
        assert result["estimated_fragment_length"] is None
        assert result["phantom_peak"] is None
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_quality_flag_integration():
    """End-to-end: quality_flag derived correctly from parsed values."""
    td = tempfile.mkdtemp(prefix="cc_parse_")
    try:
        # Good sample
        fp_good = os.path.join(td, "good.cc.qc")
        _write_file(fp_good, [
            "Filename\tnumReads\testFragLen\tNSC\tRSC\tPhantomPeak\tQualityTag",
            "good.bam\t10000000\t175\t1.12\t1.45\t225\t-1",
        ])
        parsed = parse_cc_qc_file(fp_good)
        assert _quality_flag(parsed["nsc"], parsed["rsc"]) == "ok"

        # Low NSC
        fp_low = os.path.join(td, "low_nsc.cc.qc")
        _write_file(fp_low, [
            "Filename\tnumReads\testFragLen\tNSC\tRSC\tPhantomPeak\tQualityTag",
            "low.bam\t10000000\t175\t1.02\t1.20\t225\t-1",
        ])
        parsed = parse_cc_qc_file(fp_low)
        assert _quality_flag(parsed["nsc"], parsed["rsc"]) == "low_nsc"

        # Bad file → parse_failed
        fp_bad = os.path.join(td, "bad.cc.qc")
        _write_file(fp_bad, ["not parsable"])
        parsed = parse_cc_qc_file(fp_bad)
        assert _quality_flag(parsed["nsc"], parsed["rsc"]) == "parse_failed"
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_cli_writes_basename_for_cc_qc_file():
    """CLI summary writes basename, not machine-specific absolute paths."""
    td = tempfile.mkdtemp(prefix="cc_parse_")
    try:
        nested = os.path.join(td, "sample_dir")
        os.makedirs(nested)
        fp = os.path.join(nested, "s1.cc.qc")
        out = os.path.join(td, "summary.tsv")
        _write_file(fp, [
            "s1.tagAlign\t10000000\t175\t0.85\t225\t0.60\t7\t0.55\t1.12\t1.45\t-1",
        ])

        script = os.path.join(_REPO, "scripts", "parse_cross_correlation.py")
        result = subprocess.run(
            [sys.executable, script, "--input", fp, "--output", out],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0, result.stderr

        with open(out) as fh:
            rows = [line.rstrip("\n").split("\t") for line in fh]

        header = rows[0]
        data = rows[1]
        cc_idx = header.index("cc_qc_file")
        assert data[cc_idx] == "s1.cc.qc"
        assert td not in data[cc_idx]
        assert not os.path.isabs(data[cc_idx])
    finally:
        shutil.rmtree(td, ignore_errors=True)


# --- Runner ---------------------------------------------------------------

if __name__ == "__main__":
    import traceback

    tests = [
        # _try_float
        ("try_float_valid", test_try_float_valid),
        ("try_float_none", test_try_float_none),
        ("try_float_invalid", test_try_float_invalid),
        # _quality_flag
        ("quality_flag_ok", test_quality_flag_ok),
        ("quality_flag_low_nsc", test_quality_flag_low_nsc),
        ("quality_flag_low_rsc", test_quality_flag_low_rsc),
        ("quality_flag_low_both", test_quality_flag_low_both),
        ("quality_flag_parse_failed", test_quality_flag_parse_failed),
        # sample_name_from_path
        ("sample_name_standard", test_sample_name_standard),
        ("sample_name_basename_only", test_sample_name_basename_only),
        ("sample_name_no_suffix_fallback", test_sample_name_no_suffix_fallback),
        # parse_cc_qc_file
        ("parse_standard_header_format", test_parse_standard_header_format),
        ("parse_headerless_encode_format", test_parse_headerless_encode_format),
        ("parse_comment_header_format", test_parse_comment_header_format),
        ("parse_malformed_returns_none", test_parse_malformed_returns_none),
        ("parse_empty_file", test_parse_empty_file),
        ("parse_nonexistent_file", test_parse_nonexistent_file),
        ("parse_header_with_columns_only", test_parse_header_with_columns_only),
        ("parse_partial_data", test_parse_partial_data),
        ("quality_flag_integration", test_quality_flag_integration),
        ("cli_writes_basename_for_cc_qc_file", test_cli_writes_basename_for_cc_qc_file),
    ]

    passed = 0
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"PASS: {name}")
            passed += 1
        except AssertionError as exc:
            print(f"FAIL: {name} — {exc}")
            failed += 1
        except Exception:
            print(f"FAIL: {name} — unexpected error:")
            traceback.print_exc()
            failed += 1

    print(f"\n{passed} passed, {failed} failed, {passed + failed} total")
    sys.exit(0 if failed == 0 else 1)
