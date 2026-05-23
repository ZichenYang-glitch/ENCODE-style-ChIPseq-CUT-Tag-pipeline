#!/usr/bin/env python3
"""Stage 24 unit tests for assemble_qc_summary.py and aggregate_qc_summary.py.

Tests:
  1. Complete inputs → correct 37-column TSV
  2. NA blacklist fields when has_blacklist=no
  3. NA library complexity when Picard unavailable (all-NA input)
  4. Header-only aggregation when zero inputs
  5. Header mismatch detection in aggregation
  6. Missing input file → error exit
  7. Empty peak_counts (0 peaks) handled
  8. Byte-compare header with expected 37-column header
"""

import csv
import os
import shutil
import subprocess
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ASSEMBLE = os.path.join(REPO_ROOT, "scripts", "assemble_qc_summary.py")
AGGREGATE = os.path.join(REPO_ROOT, "scripts", "aggregate_qc_summary.py")

_NA = "NA"

# Expected 37-column header from assemble_qc_summary.py
_EXPECTED_HEADER = [
    "sample", "assay", "target", "genome", "layout", "peak_mode",
    "use_control", "control_type", "final_bam", "peaks",
    "blacklist", "blacklist_filtered_bam", "blacklist_filtered_peaks",
    "total_reads", "reads_in_peaks", "frip", "peak_count",
    "blacklist_filtered_peak_count",
    "metrics_source", "unpaired_reads_examined", "read_pairs_examined",
    "secondary_or_supplementary_reads", "unmapped_reads",
    "unpaired_read_duplicates", "read_pair_duplicates",
    "read_pair_optical_duplicates", "percent_duplication",
    "estimated_library_size", "total_reads_examined",
    "duplicate_reads_estimate",
    "total_fragments", "distinct_fragments", "one_read_fragments",
    "two_read_fragments", "nrf", "pbc1", "pbc2",
]

PASSED = 0
TOTAL = 0


def _record(name, passed):
    global PASSED, TOTAL
    TOTAL += 1
    if passed:
        PASSED += 1
        print("PASS: %s" % name)
    else:
        print("FAIL: %s" % name)


def _write_tsv(path, header, rows):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(header)
        for row in rows:
            w.writerow(row)


# ---------------------------------------------------------------------------
# Test 1: Complete inputs → correct 37-column TSV
# ---------------------------------------------------------------------------

def test_complete_inputs():
    workdir = tempfile.mkdtemp(prefix="s24_t1_", dir="/tmp")
    try:
        # Write component TSVs
        _write_tsv(
            os.path.join(workdir, "peak_counts.tsv"),
            ["sample", "peak_mode", "peaks", "blacklist_filtered_peaks"],
            [["S1", "narrow", "15000", "12000"]],
        )
        _write_tsv(
            os.path.join(workdir, "frip.tsv"),
            ["sample", "total_reads", "reads_in_peaks", "frip", "bam", "peaks"],
            [["S1", "10000000", "2500000", "0.250000", "/bam", "/peaks"]],
        )
        _write_tsv(
            os.path.join(workdir, "library_complexity.tsv"),
            ["sample", "metrics_source", "unpaired_reads_examined",
             "read_pairs_examined", "secondary_or_supplementary_reads",
             "unmapped_reads", "unpaired_read_duplicates",
             "read_pair_duplicates", "read_pair_optical_duplicates",
             "percent_duplication", "estimated_library_size",
             "total_reads_examined", "duplicate_reads_estimate"],
            [["S1", "picard_markduplicates", "1000", "4500000", "50000",
              "200000", "50000", "400000", "5000", "0.100000",
              "8500000", "9001000", "450000"]],
        )
        _write_tsv(
            os.path.join(workdir, "nrf_pbc.tsv"),
            ["sample", "total_fragments", "distinct_fragments",
             "one_read_fragments", "two_read_fragments",
             "nrf", "pbc1", "pbc2"],
            [["S1", "9000000", "7500000", "6000000", "1500000",
              "0.833333", "0.800000", "4.000000"]],
        )

        out = os.path.join(workdir, "qc_summary.tsv")
        result = subprocess.run(
            [sys.executable, ASSEMBLE,
             "--sample", "S1",
             "--assay", "chipseq",
             "--target", "CTCF",
             "--genome", "hs",
             "--layout", "PE",
             "--peak-mode", "narrow",
             "--use-control", "False",
             "--control-type", "none",
             "--final-bam", "results/S1/02_align/S1.final.bam",
             "--peaks-file", "results/S1/04_peaks/S1/S1_peaks.narrowPeak",
             "--has-blacklist", "yes",
             "--blacklist", "/opt/genomes/bl.bed",
             "--bl-bam", "results/S1/02_align/S1.blacklist_filtered.bam",
             "--bl-peaks", "results/S1/04_peaks/S1_bl/S1_peaks.bl.narrowPeak",
             "--peak-counts", os.path.join(workdir, "peak_counts.tsv"),
             "--frip", os.path.join(workdir, "frip.tsv"),
             "--library-complexity", os.path.join(workdir, "library_complexity.tsv"),
             "--nrf-pbc", os.path.join(workdir, "nrf_pbc.tsv"),
             "--output", out,
             ],
            capture_output=True, text=True,
        )

        passed = (result.returncode == 0)
        if not passed:
            print("  stderr:", result.stderr.strip()[-200:])
            _record("1-complete_inputs_37cols", False)
            return

        with open(out, newline="") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader)
            data = next(reader)

        passed = (header == _EXPECTED_HEADER
                  and len(data) == 37
                  and data[0] == "S1"
                  and data[15] == "0.250000"  # frip
                  and data[16] == "15000"       # peak_count
                  and data[17] == "12000"       # bl_peak_count
                  and data[34] == "0.833333"    # nrf (col 34)
                  )
        if not passed:
            print("  header_match=%s len=%d data[0]=%s" % (
                header == _EXPECTED_HEADER, len(data), data[0] if data else "NONE"))
        _record("1-complete_inputs_37cols", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 2: NA blacklist fields when has_blacklist=no
# ---------------------------------------------------------------------------

def test_na_blacklist():
    workdir = tempfile.mkdtemp(prefix="s24_t2_", dir="/tmp")
    try:
        _write_tsv(
            os.path.join(workdir, "peak_counts.tsv"),
            ["sample", "peak_mode", "peaks", "blacklist_filtered_peaks"],
            [["S1", "narrow", "100", "80"]],
        )
        _write_tsv(
            os.path.join(workdir, "frip.tsv"),
            ["sample", "total_reads", "reads_in_peaks", "frip", "bam", "peaks"],
            [["S1", "1000", "100", "0.100000", "/b", "/p"]],
        )
        _write_tsv(
            os.path.join(workdir, "library_complexity.tsv"),
            ["sample", "metrics_source", "unpaired_reads_examined",
             "read_pairs_examined", "secondary_or_supplementary_reads",
             "unmapped_reads", "unpaired_read_duplicates",
             "read_pair_duplicates", "read_pair_optical_duplicates",
             "percent_duplication", "estimated_library_size",
             "total_reads_examined", "duplicate_reads_estimate"],
            [["S1", "picard", "0", "500", "0", "0", "0", "0", "0",
              "0", "1000", "1000", "0"]],
        )
        _write_tsv(
            os.path.join(workdir, "nrf_pbc.tsv"),
            ["sample", "total_fragments", "distinct_fragments",
             "one_read_fragments", "two_read_fragments",
             "nrf", "pbc1", "pbc2"],
            [["S1", "800", "700", "600", "100", "0.875", "0.857", "6.0"]],
        )

        out = os.path.join(workdir, "qc_summary.tsv")
        result = subprocess.run(
            [sys.executable, ASSEMBLE,
             "--sample", "S1",
             "--assay", "chipseq", "--target", "CTCF", "--genome", "hs",
             "--layout", "PE", "--peak-mode", "narrow",
             "--use-control", "False", "--control-type", "none",
             "--final-bam", "/bam", "--peaks-file", "/peaks",
             "--has-blacklist", "no",
             "--peak-counts", os.path.join(workdir, "peak_counts.tsv"),
             "--frip", os.path.join(workdir, "frip.tsv"),
             "--library-complexity", os.path.join(workdir, "library_complexity.tsv"),
             "--nrf-pbc", os.path.join(workdir, "nrf_pbc.tsv"),
             "--output", out,
             ],
            capture_output=True, text=True,
        )

        if result.returncode != 0:
            _record("2-na_blacklist", False)
            return

        with open(out, newline="") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)  # header
            data = next(reader)

        # Columns 10,11,12 = blacklist, blacklist_filtered_bam, blacklist_filtered_peaks
        # Column 17 = blacklist_filtered_peak_count
        passed = (data[10] == _NA and data[11] == _NA and data[12] == _NA
                  and data[17] == _NA)
        if not passed:
            print("  bl=%s bl_bam=%s bl_peaks=%s bl_count=%s" % (
                data[10], data[11], data[12], data[17]))
        _record("2-na_blacklist", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 3: NA library complexity when Picard unavailable
# ---------------------------------------------------------------------------

def test_na_library_complexity():
    workdir = tempfile.mkdtemp(prefix="s24_t3_", dir="/tmp")
    try:
        _write_tsv(
            os.path.join(workdir, "peak_counts.tsv"),
            ["sample", "peak_mode", "peaks", "blacklist_filtered_peaks"],
            [["S1", "narrow", "50", _NA]],
        )
        _write_tsv(
            os.path.join(workdir, "frip.tsv"),
            ["sample", "total_reads", "reads_in_peaks", "frip", "bam", "peaks"],
            [["S1", "500", "50", "0.100000", "/b", "/p"]],
        )
        # All-NA library complexity (Picard unavailable)
        _write_tsv(
            os.path.join(workdir, "library_complexity.tsv"),
            ["sample", "metrics_source", "unpaired_reads_examined",
             "read_pairs_examined", "secondary_or_supplementary_reads",
             "unmapped_reads", "unpaired_read_duplicates",
             "read_pair_duplicates", "read_pair_optical_duplicates",
             "percent_duplication", "estimated_library_size",
             "total_reads_examined", "duplicate_reads_estimate"],
            [["S1", "fallback", _NA, _NA, _NA, _NA, _NA, _NA, _NA,
              _NA, _NA, _NA, _NA]],
        )
        _write_tsv(
            os.path.join(workdir, "nrf_pbc.tsv"),
            ["sample", "total_fragments", "distinct_fragments",
             "one_read_fragments", "two_read_fragments",
             "nrf", "pbc1", "pbc2"],
            [["S1", "400", "350", "300", "50", "0.875", "0.857", "6.0"]],
        )

        out = os.path.join(workdir, "qc_summary.tsv")
        result = subprocess.run(
            [sys.executable, ASSEMBLE,
             "--sample", "S1", "--assay", "chipseq", "--target", "CTCF",
             "--genome", "hs", "--layout", "PE", "--peak-mode", "narrow",
             "--use-control", "False", "--control-type", "none",
             "--final-bam", "/bam", "--peaks-file", "/peaks",
             "--has-blacklist", "yes",
             "--blacklist", "/bl.bed",
             "--peak-counts", os.path.join(workdir, "peak_counts.tsv"),
             "--frip", os.path.join(workdir, "frip.tsv"),
             "--library-complexity", os.path.join(workdir, "library_complexity.tsv"),
             "--nrf-pbc", os.path.join(workdir, "nrf_pbc.tsv"),
             "--output", out,
             ],
            capture_output=True, text=True,
        )

        if result.returncode != 0:
            _record("3-na_library_complexity", False)
            return

        with open(out, newline="") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)
            data = next(reader)

        # Column 18 = metrics_source (should be "fallback")
        # Column 19 = unpaired_reads_examined (should be NA)
        passed = (data[18] == "fallback" and data[19] == _NA)
        _record("3-na_library_complexity", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 4: Header-only aggregation when zero inputs
# ---------------------------------------------------------------------------

def test_aggregate_zero_inputs():
    workdir = tempfile.mkdtemp(prefix="s24_t4_", dir="/tmp")
    try:
        out = os.path.join(workdir, "stage3_qc_summary.tsv")
        result = subprocess.run(
            [sys.executable, AGGREGATE, "--output", out],
            capture_output=True, text=True,
        )

        if result.returncode != 0:
            _record("4-aggregate_zero_inputs", False)
            return

        with open(out, newline="") as f:
            reader = csv.reader(f, delimiter="\t")
            header = next(reader)
            rows = list(reader)

        passed = (header == _EXPECTED_HEADER and len(rows) == 0)
        _record("4-aggregate_zero_inputs", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 5: Header mismatch detection in aggregation
# ---------------------------------------------------------------------------

def test_aggregate_header_mismatch():
    workdir = tempfile.mkdtemp(prefix="s24_t5_", dir="/tmp")
    try:
        good = os.path.join(workdir, "good.tsv")
        bad = os.path.join(workdir, "bad.tsv")

        _write_tsv(good, _EXPECTED_HEADER, [["S1"] + [_NA] * 36])
        _write_tsv(bad, ["wrong", "header"], [["x", "y"]])

        out = os.path.join(workdir, "stage3.tsv")
        result = subprocess.run(
            [sys.executable, AGGREGATE, "--output", out, good, bad],
            capture_output=True, text=True,
        )

        passed = (result.returncode != 0 and "header mismatch" in result.stderr.lower())
        if not passed:
            print("  rc=%s stderr=%s" % (result.returncode, result.stderr[-200:]))
        _record("5-aggregate_header_mismatch", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 6: Missing input file → error exit
# ---------------------------------------------------------------------------

def test_missing_input_file():
    workdir = tempfile.mkdtemp(prefix="s24_t6_", dir="/tmp")
    try:
        # Create one valid file so the others don't fail first
        _write_tsv(
            os.path.join(workdir, "pc.tsv"),
            ["sample", "peak_mode", "peaks", "blacklist_filtered_peaks"],
            [["S1", "narrow", "10", _NA]],
        )

        out = os.path.join(workdir, "out.tsv")
        result = subprocess.run(
            [sys.executable, ASSEMBLE,
             "--sample", "S1", "--assay", "c", "--target", "T",
             "--genome", "hs", "--layout", "PE", "--peak-mode", "narrow",
             "--use-control", "False", "--control-type", "none",
             "--final-bam", "/b", "--peaks-file", "/p",
             "--has-blacklist", "no",
             "--peak-counts", os.path.join(workdir, "pc.tsv"),
             "--frip", os.path.join(workdir, "nonexistent.tsv"),
             "--library-complexity", os.path.join(workdir, "pc.tsv"),
             "--nrf-pbc", os.path.join(workdir, "pc.tsv"),
             "--output", out,
             ],
            capture_output=True, text=True,
        )

        passed = (result.returncode != 0 and "not found" in result.stderr.lower())
        if not passed:
            print("  rc=%s stderr=%s" % (result.returncode, result.stderr[-200:]))
        _record("6-missing_input_file", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 7: Empty peak_counts (0 peaks) handled
# ---------------------------------------------------------------------------

def test_zero_peaks():
    workdir = tempfile.mkdtemp(prefix="s24_t7_", dir="/tmp")
    try:
        _write_tsv(
            os.path.join(workdir, "peak_counts.tsv"),
            ["sample", "peak_mode", "peaks", "blacklist_filtered_peaks"],
            [["S1", "narrow", "0", _NA]],
        )
        _write_tsv(
            os.path.join(workdir, "frip.tsv"),
            ["sample", "total_reads", "reads_in_peaks", "frip", "bam", "peaks"],
            [["S1", "0", "0", _NA, "/b", "/p"]],
        )
        _write_tsv(
            os.path.join(workdir, "library_complexity.tsv"),
            ["sample", "metrics_source", "unpaired_reads_examined",
             "read_pairs_examined", "secondary_or_supplementary_reads",
             "unmapped_reads", "unpaired_read_duplicates",
             "read_pair_duplicates", "read_pair_optical_duplicates",
             "percent_duplication", "estimated_library_size",
             "total_reads_examined", "duplicate_reads_estimate"],
            [["S1", "fallback", _NA, _NA, _NA, _NA, _NA, _NA, _NA,
              _NA, _NA, _NA, _NA]],
        )
        _write_tsv(
            os.path.join(workdir, "nrf_pbc.tsv"),
            ["sample", "total_fragments", "distinct_fragments",
             "one_read_fragments", "two_read_fragments",
             "nrf", "pbc1", "pbc2"],
            [["S1", "0", "0", "0", "0", _NA, _NA, _NA]],
        )

        out = os.path.join(workdir, "qc_summary.tsv")
        result = subprocess.run(
            [sys.executable, ASSEMBLE,
             "--sample", "S1", "--assay", "c", "--target", "T",
             "--genome", "hs", "--layout", "PE", "--peak-mode", "narrow",
             "--use-control", "False", "--control-type", "none",
             "--final-bam", "/b", "--peaks-file", "/p",
             "--has-blacklist", "no",
             "--peak-counts", os.path.join(workdir, "peak_counts.tsv"),
             "--frip", os.path.join(workdir, "frip.tsv"),
             "--library-complexity", os.path.join(workdir, "library_complexity.tsv"),
             "--nrf-pbc", os.path.join(workdir, "nrf_pbc.tsv"),
             "--output", out,
             ],
            capture_output=True, text=True,
        )

        if result.returncode != 0:
            _record("7-zero_peaks", False)
            return

        with open(out, newline="") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader)
            data = next(reader)

        passed = (data[16] == "0"   # peak_count
                  and data[15] == _NA)  # frip (total_reads=0)
        _record("7-zero_peaks", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 8: Byte-compare header
# ---------------------------------------------------------------------------

def test_header_byte_compare():
    """Verify the header constant matches expected and output uses LF newlines."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "assemble_qc_summary", ASSEMBLE)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)

    script_header = mod._QC_SUMMARY_COLUMNS
    passed = (script_header == _EXPECTED_HEADER)
    if not passed:
        for i, (s, e) in enumerate(zip(script_header, _EXPECTED_HEADER)):
            if s != e:
                print("  First diff at col %d: script=%r expected=%r" % (i, s, e))
                break
    _record("8-header_byte_compare", passed)


# ---------------------------------------------------------------------------
# Test 9: No CRLF in output
# ---------------------------------------------------------------------------

def test_no_crlf():
    """Verify assembled output does not contain CRLF line terminators."""
    workdir = tempfile.mkdtemp(prefix="s24_t9_", dir="/tmp")
    try:
        _write_tsv(
            os.path.join(workdir, "peak_counts.tsv"),
            ["sample", "peak_mode", "peaks", "blacklist_filtered_peaks"],
            [["S1", "narrow", "10", _NA]],
        )
        _write_tsv(
            os.path.join(workdir, "frip.tsv"),
            ["sample", "total_reads", "reads_in_peaks", "frip", "bam", "peaks"],
            [["S1", "100", "10", "0.1", "/b", "/p"]],
        )
        lc_hdr = ["sample", "metrics_source", "unpaired_reads_examined",
                   "read_pairs_examined", "secondary_or_supplementary_reads",
                   "unmapped_reads", "unpaired_read_duplicates",
                   "read_pair_duplicates", "read_pair_optical_duplicates",
                   "percent_duplication", "estimated_library_size",
                   "total_reads_examined", "duplicate_reads_estimate"]
        _write_tsv(
            os.path.join(workdir, "library_complexity.tsv"),
            lc_hdr,
            [["S1", "fallback"] + [_NA] * 11],
        )
        _write_tsv(
            os.path.join(workdir, "nrf_pbc.tsv"),
            ["sample", "total_fragments", "distinct_fragments",
             "one_read_fragments", "two_read_fragments",
             "nrf", "pbc1", "pbc2"],
            [["S1", "10", "8", "6", "2", "0.8", "0.75", "3.0"]],
        )

        out = os.path.join(workdir, "qc_summary.tsv")
        result = subprocess.run(
            [sys.executable, ASSEMBLE,
             "--sample", "S1", "--assay", "c", "--target", "T",
             "--genome", "hs", "--layout", "PE", "--peak-mode", "narrow",
             "--use-control", "False", "--control-type", "none",
             "--final-bam", "results/S1/02_align/S1.final.bam",
             "--peaks-file", "results/S1/04_peaks/S1/S1_peaks.narrowPeak",
             "--has-blacklist", "no",
             "--peak-counts", os.path.join(workdir, "peak_counts.tsv"),
             "--frip", os.path.join(workdir, "frip.tsv"),
             "--library-complexity", os.path.join(workdir, "library_complexity.tsv"),
             "--nrf-pbc", os.path.join(workdir, "nrf_pbc.tsv"),
             "--output", out,
             ],
            capture_output=True, text=True,
        )

        if result.returncode != 0:
            _record("9-no_crlf", False)
            return

        with open(out, "rb") as fh:
            raw = fh.read()

        passed = b"\r\n" not in raw
        if not passed:
            print("  CRLF found in output")
        _record("9-no_crlf", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    print("Starting Stage 24 QC Summary Unit Tests\n")

    test_complete_inputs()
    test_na_blacklist()
    test_na_library_complexity()
    test_aggregate_zero_inputs()
    test_aggregate_header_mismatch()
    test_missing_input_file()
    test_zero_peaks()
    test_header_byte_compare()
    test_no_crlf()

    print("\nSummary: %d/%d tests passed." % (PASSED, TOTAL))

    pycache = os.path.join(REPO_ROOT, "scripts", "__pycache__")
    if os.path.isdir(pycache):
        shutil.rmtree(pycache)

    if PASSED < TOTAL:
        sys.exit(1)


if __name__ == "__main__":
    main()
