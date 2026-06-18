"""Stage 54 consensus engine unit tests.

31 synthetic tests covering core algorithm, graph topology, CLI validation,
output formats, edge cases, and summary TSV correctness.

All tests invoke scripts/compute_consensus.py via subprocess with
temporary synthetic peak files. No real data, no pytest, no DAG.
"""

import json
import os
import subprocess
import sys
import tempfile

SCRIPT = os.path.join(
    os.path.dirname(__file__), "..", "scripts", "compute_consensus.py"
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def narrowpeak_row(chrom, start, end, name=".", score=100, strand=".",
                   signal=5.0, pval="3.0", qval="2.0", peak=100):
    return "\t".join(str(x) for x in [
        chrom, start, end, name, score, strand, signal, pval, qval, peak,
    ])


def broadpeak_row(chrom, start, end, name=".", score=100, strand=".",
                  signal=5.0, pval="3.0", qval="2.0"):
    return "\t".join(str(x) for x in [
        chrom, start, end, name, score, strand, signal, pval, qval,
    ])


def bed_row(chrom, start, end, name=".", score=100, strand="."):
    return "\t".join(str(x) for x in [
        chrom, start, end, name, score, strand,
    ])


def bed3_row(chrom, start, end):
    return "\t".join(str(x) for x in [chrom, start, end])


def make_peak_file(content_lines, suffix=".narrowPeak"):
    """Write synthetic peak content to a temp file. Returns path."""
    fd, path = tempfile.mkstemp(suffix=suffix, prefix="consensus_test_")
    with os.fdopen(fd, "w") as fh:
        for line in content_lines:
            fh.write(line + "\n")
    return path


def run_consensus(peaks, bioreps, fmt="narrowPeak", min_rep=2, overlap=0.5,
                  experiment="test_exp", assay="chipseq", caller="macs3",
                  peak_mode="narrow", final_method="", final_output=""):
    """Run compute_consensus.py and return (returncode, stderr, peak_content, summary_content)."""
    out_fd, out_path = tempfile.mkstemp(suffix=".peaks", prefix="consensus_out_")
    os.close(out_fd)
    sum_fd, sum_path = tempfile.mkstemp(suffix=".tsv", prefix="consensus_sum_")
    os.close(sum_fd)

    cmd = [
        sys.executable, SCRIPT,
        "--peaks"] + peaks + [
        "--bioreps"] + bioreps + [
        "--format", fmt,
        "--min-replicates", str(min_rep),
        "--reciprocal-overlap", str(overlap),
        "--output", out_path,
        "--summary", sum_path,
        "--experiment", experiment,
        "--assay", assay,
        "--caller", caller,
        "--peak-mode", peak_mode,
        "--final-method", final_method,
        "--final-output", final_output,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    # Read output files
    peak_content = ""
    summary_content = ""
    if os.path.exists(out_path):
        with open(out_path) as fh:
            peak_content = fh.read()
    if os.path.exists(sum_path):
        with open(sum_path) as fh:
            summary_content = fh.read()

    return result.returncode, result.stderr, peak_content, summary_content


# ---------------------------------------------------------------------------
# Main test runner
# ---------------------------------------------------------------------------


def main():
    passed = 0
    total = 0

    def check(name, fn):
        nonlocal passed, total
        total += 1
        try:
            if fn():
                print("PASS: %s" % name)
                passed += 1
        except Exception as e:
            print("FAIL: %s\n   Unexpected error: %s" % (name, e))

    # ===================================================================
    # Core algorithm (T1–T8)
    # ===================================================================

    # T1: 2 bioreps, perfect overlap
    def t1():
        p1 = make_peak_file([
            narrowpeak_row("chr1", 100, 500, "peak1", signal=5.0, peak=200),
        ])
        p2 = make_peak_file([
            narrowpeak_row("chr1", 100, 500, "peak1", signal=5.0, peak=200),
        ])
        rc, stderr, peaks, summary = run_consensus(
            [p1, p2], ["1", "2"])
        if rc != 0:
            print("   exit code %d: %s" % (rc, stderr[:200]))
            return False
        lines = [l for l in peaks.strip().split("\n") if l.strip()]
        if len(lines) != 1:
            print("   expected 1 consensus peak, got %d" % len(lines))
            return False
        fields = lines[0].split("\t")
        if len(fields) != 10:
            print("   expected 10 narrowPeak columns, got %d" % len(fields))
            return False
        # score = 1000 * 2/2 = 1000
        if fields[4] != "1000":
            print("   expected score=1000, got %s" % fields[4])
            return False
        return True
    check("T1: 2 bioreps perfect overlap → 1 consensus, score=1000", t1)

    # T2: 2 bioreps, no overlap
    def t2():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 500)])
        p2 = make_peak_file([narrowpeak_row("chr1", 600, 1000)])
        rc, stderr, peaks, summary = run_consensus([p1, p2], ["1", "2"])
        if rc != 0:
            print("   exit code %d" % rc)
            return False
        lines = [l for l in peaks.strip().split("\n") if l.strip()]
        if len(lines) != 0:
            print("   expected 0 consensus peaks, got %d" % len(lines))
            return False
        # Summary should have consensus_peak_count=0
        if "consensus_peak_count" not in summary:
            print("   summary missing")
            return False
        return True
    check("T2: 2 bioreps no overlap → 0 consensus", t2)

    # T3: 3 bioreps, N-of-M=2
    def t3():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 500, signal=5.0, peak=200)])
        p2 = make_peak_file([narrowpeak_row("chr1", 100, 500, signal=5.0, peak=200)])
        p3 = make_peak_file([narrowpeak_row("chr1", 600, 1000, signal=5.0, peak=200)])
        rc, stderr, peaks, summary = run_consensus(
            [p1, p2, p3], ["1", "2", "3"])
        if rc != 0:
            print("   exit code %d" % rc)
            return False
        lines = [l for l in peaks.strip().split("\n") if l.strip()]
        if len(lines) != 1:
            print("   expected 1 consensus, got %d" % len(lines))
            return False
        fields = lines[0].split("\t")
        # score = int(1000 * 2/3 + 0.5) = 667
        if fields[4] != "667":
            print("   expected score=667, got %s" % fields[4])
            return False
        return True
    check("T3: 3 bioreps N-of-M=2 → 1 consensus, score=667", t3)

    # T4: 3 bioreps, all 3 overlap
    def t4():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 500, signal=5.0, peak=200)])
        p2 = make_peak_file([narrowpeak_row("chr1", 100, 500, signal=5.0, peak=200)])
        p3 = make_peak_file([narrowpeak_row("chr1", 100, 500, signal=5.0, peak=200)])
        rc, stderr, peaks, summary = run_consensus(
            [p1, p2, p3], ["1", "2", "3"])
        if rc != 0:
            print("   exit code %d" % rc)
            return False
        fields = peaks.strip().split("\t")
        if fields[4] != "1000":
            print("   expected score=1000 (3/3), got %s" % fields[4])
            return False
        return True
    check("T4: 3 bioreps all overlap → score=1000", t4)

    # T5: reciprocal_overlap fails (500bp peak, 50bp overlap = 10%)
    def t5():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 600)])  # 500bp
        p2 = make_peak_file([narrowpeak_row("chr1", 550, 1050)])  # 50bp overlap
        rc, stderr, peaks, summary = run_consensus([p1, p2], ["1", "2"])
        lines = [l for l in peaks.strip().split("\n") if l.strip()]
        if len(lines) != 0:
            print("   expected 0 consensus (overlap 10% < 50%), got %d" % len(lines))
            return False
        return True
    check("T5: reciprocal_overlap fails (10%)", t5)

    # T6: reciprocal_overlap passes (500bp peak, 400bp overlap = 80%)
    def t6():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 600)])  # 500bp
        p2 = make_peak_file([narrowpeak_row("chr1", 200, 700)])  # 400bp overlap
        rc, stderr, peaks, summary = run_consensus([p1, p2], ["1", "2"])
        lines = [l for l in peaks.strip().split("\n") if l.strip()]
        if len(lines) != 1:
            print("   expected 1 consensus, got %d" % len(lines))
            return False
        return True
    check("T6: reciprocal_overlap passes (80%)", t6)

    # T7: reciprocal_overlap at exact threshold 0.5
    def t7():
        # peak A: 100-300 (200bp), peak B: 200-400 (200bp)
        # overlap = 100bp. 100/200 = 0.5 for both
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 300)])
        p2 = make_peak_file([narrowpeak_row("chr1", 200, 400)])
        rc, stderr, peaks, summary = run_consensus(
            [p1, p2], ["1", "2"], overlap=0.5)
        lines = [l for l in peaks.strip().split("\n") if l.strip()]
        if len(lines) != 1:
            print("   expected 1 consensus at threshold, got %d" % len(lines))
            return False
        return True
    check("T7: reciprocal_overlap at exact threshold (0.5)", t7)

    # T8: Same biorep multiple peaks in one file → count once
    def t8():
        # biorep1 has 2 overlapping peaks, biorep2 has 1
        p1 = make_peak_file([
            narrowpeak_row("chr1", 100, 500, signal=5.0, peak=200),
            narrowpeak_row("chr1", 120, 480, signal=4.0, peak=150),
        ])
        p2 = make_peak_file([
            narrowpeak_row("chr1", 100, 500, signal=3.0, peak=200),
        ])
        rc, stderr, peaks, summary = run_consensus([p1, p2], ["1", "2"])
        if rc != 0:
            print("   exit code %d" % rc)
            return False
        lines = [l for l in peaks.strip().split("\n") if l.strip()]
        if len(lines) != 1:
            print("   expected 1 consensus, got %d" % len(lines))
            return False
        # Both biorep1 and biorep2 support → support=2, score=1000
        fields = lines[0].split("\t")
        if fields[4] != "1000":
            print("   expected score=1000 (support=2 distinct bioreps), got %s" % fields[4])
            return False
        return True
    check("T8: same biorep multiple peaks do not add extra support", t8)

    # ===================================================================
    # Graph topology (T9–T12)
    # ===================================================================

    # T9: Chained overlap A-B-C
    def t9():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 300)])  # A
        p2 = make_peak_file([narrowpeak_row("chr1", 200, 400)])  # B (overlaps A, C)
        p3 = make_peak_file([narrowpeak_row("chr1", 300, 500)])  # C (overlaps B, not A)
        rc, stderr, peaks, summary = run_consensus(
            [p1, p2, p3], ["1", "2", "3"])
        lines = [l for l in peaks.strip().split("\n") if l.strip()]
        if len(lines) != 1:
            print("   expected 1 chained component, got %d" % len(lines))
            return False
        fields = lines[0].split("\t")
        # Component: min_start=100, max_end=500
        if fields[1] != "100" or fields[2] != "500":
            print("   expected merged 100-500, got %s-%s" % (fields[1], fields[2]))
            return False
        if fields[4] != "1000":
            print("   expected score=1000 (3/3), got %s" % fields[4])
            return False
        return True
    check("T9: chained overlap A-B-C → 1 component, support=3", t9)

    # T10: Adjacent intervals, no overlap
    def t10():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        p2 = make_peak_file([narrowpeak_row("chr1", 200, 300)])
        rc, stderr, peaks, summary = run_consensus([p1, p2], ["1", "2"])
        lines = [l for l in peaks.strip().split("\n") if l.strip()]
        if len(lines) != 0:
            print("   expected 0 consensus, got %d" % len(lines))
            return False
        return True
    check("T10: adjacent intervals (100-200, 200-300) no overlap", t10)

    # T11: Cross-chromosome separation
    def t11():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        p2 = make_peak_file([narrowpeak_row("chr2", 100, 200)])
        rc, stderr, peaks, summary = run_consensus([p1, p2], ["1", "2"])
        lines = [l for l in peaks.strip().split("\n") if l.strip()]
        if len(lines) != 0:
            print("   expected 0 cross-chrom consensus, got %d" % len(lines))
            return False
        return True
    check("T11: cross-chromosome never clusters", t11)

    # T12: Unsorted input → deterministic sorted output
    def t12():
        p1 = make_peak_file([
            narrowpeak_row("chr2", 100, 300),
            narrowpeak_row("chr1", 500, 700),
        ])
        p2 = make_peak_file([
            narrowpeak_row("chr2", 100, 300),
            narrowpeak_row("chr1", 500, 700),
        ])
        rc, stderr, peaks, summary = run_consensus([p1, p2], ["1", "2"])
        lines = [l for l in peaks.strip().split("\n") if l.strip()]
        if len(lines) != 2:
            print("   expected 2 consensus peaks, got %d" % len(lines))
            return False
        # chr1 peak should come before chr2 peak
        if "chr1" not in lines[0] or "chr2" not in lines[1]:
            print("   expected sorted output (chr1 before chr2)")
            return False
        return True
    check("T12: unsorted input → sorted output", t12)

    # ===================================================================
    # CLI validation (T13–T17)
    # ===================================================================

    # T13: Mismatched peaks/bioreps count
    def t13():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        p2 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        rc, stderr, peaks, summary = run_consensus(
            [p1, p2], ["1"])  # 2 peaks, 1 biorep
        if rc == 0:
            print("   expected non-zero exit")
            return False
        if "biorep" not in stderr.lower() and "peak" not in stderr.lower():
            print("   expected error about count mismatch, got: %s" % stderr[:100])
            return False
        return True
    check("T13: mismatched peaks/bioreps → non-zero exit", t13)

    # T14: Duplicate biorep labels
    def t14():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        p2 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        rc, stderr, peaks, summary = run_consensus(
            [p1, p2], ["1", "1"])
        if rc == 0:
            print("   expected non-zero exit for duplicate labels")
            return False
        if "duplicate" not in stderr.lower():
            print("   expected 'duplicate' error, got: %s" % stderr[:100])
            return False
        return True
    check("T14: duplicate biorep labels → non-zero exit", t14)

    # T15: min-replicates=1 rejected
    def t15():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        p2 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        rc, stderr, peaks, summary = run_consensus(
            [p1, p2], ["1", "2"], min_rep=1)
        if rc == 0:
            print("   expected non-zero exit for min_replicates=1")
            return False
        if "min-replicates" not in stderr.lower():
            print("   expected error about min-replicates")
            return False
        return True
    check("T15: min-replicates=1 rejected", t15)

    # T16: min-replicates > n_bioreps rejected
    def t16():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        p2 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        rc, stderr, peaks, summary = run_consensus(
            [p1, p2], ["1", "2"], min_rep=3)
        if rc == 0:
            print("   expected non-zero exit")
            return False
        return True
    check("T16: min-replicates > n_bioreps rejected", t16)

    # T17: Only 1 peak file rejected
    def t17():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        rc, stderr, peaks, summary = run_consensus(
            [p1], ["1"])
        if rc == 0:
            print("   expected non-zero exit for single peak file")
            return False
        if "least two" not in stderr.lower() and "at least" not in stderr.lower():
            print("   expected 'at least two' error, got: %s" % stderr[:100])
            return False
        return True
    check("T17: only 1 peak file → non-zero exit", t17)

    # ===================================================================
    # Output format (T18–T21)
    # ===================================================================

    # T18: narrowPeak output shape
    def t18():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 500, signal=5.0, peak=200)])
        p2 = make_peak_file([narrowpeak_row("chr1", 100, 500, signal=5.0, peak=200)])
        rc, stderr, peaks, summary = run_consensus(
            [p1, p2], ["1", "2"], fmt="narrowPeak")
        if rc != 0:
            print("   exit code %d" % rc)
            return False
        lines = [l for l in peaks.strip().split("\n") if l.strip()]
        if len(lines) != 1:
            print("   expected 1 peak, got %d" % len(lines))
            return False
        fields = lines[0].split("\t")
        if len(fields) != 10:
            print("   expected 10 columns, got %d" % len(fields))
            return False
        # Summit should be valid (within [0, end-start))
        summit = int(fields[9])
        span = int(fields[2]) - int(fields[1])
        if not (0 <= summit < span):
            print("   summit %d out of [0, %d)" % (summit, span))
            return False
        return True
    check("T18: narrowPeak output 10 cols, valid summit", t18)

    # T19: broadPeak output shape
    def t19():
        p1 = make_peak_file([
            broadpeak_row("chr1", 100, 500, signal=5.0),
        ])
        p2 = make_peak_file([
            broadpeak_row("chr1", 100, 500, signal=5.0),
        ])
        rc, stderr, peaks, summary = run_consensus(
            [p1, p2], ["1", "2"], fmt="broadPeak")
        if rc != 0:
            print("   exit code %d" % rc)
            return False
        lines = [l for l in peaks.strip().split("\n") if l.strip()]
        if len(lines) != 1:
            print("   expected 1 peak, got %d" % len(lines))
            return False
        fields = lines[0].split("\t")
        if len(fields) != 9:
            print("   expected 9 broadPeak columns, got %d" % len(fields))
            return False
        return True
    check("T19: broadPeak output 9 cols", t19)

    # T20: BED output shape
    def t20():
        p1 = make_peak_file([bed_row("chr1", 100, 500)])
        p2 = make_peak_file([bed_row("chr1", 100, 500)])
        rc, stderr, peaks, summary = run_consensus(
            [p1, p2], ["1", "2"], fmt="bed")
        if rc != 0:
            print("   exit code %d" % rc)
            return False
        lines = [l for l in peaks.strip().split("\n") if l.strip()]
        if len(lines) != 1:
            print("   expected 1 BED peak, got %d" % len(lines))
            return False
        fields = lines[0].split("\t")
        if len(fields) != 6:
            print("   expected 6 BED6 columns, got %d" % len(fields))
            return False
        return True
    check("T20: BED output 6 cols (BED6)", t20)

    # T21: Score calculation
    def t21():
        # 2/2 = 1000
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        p2 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        rc, _, peaks, _ = run_consensus([p1, p2], ["1", "2"])
        if rc != 0 or peaks.strip().split("\t")[4] != "1000":
            print("   2/2 expected 1000")
            return False

        # 2/3 = int(1000*2/3 + 0.5) = 667
        p3 = make_peak_file([narrowpeak_row("chr1", 600, 700)])
        rc, _, peaks2, _ = run_consensus([p1, p2, p3], ["1", "2", "3"])
        if rc != 0 or peaks2.strip().split("\t")[4] != "667":
            print("   2/3 expected 667, got %s" % (
                peaks2.strip().split("\t")[4] if peaks2.strip() else "empty"))
            return False

        # 3/4 = int(1000*3/4 + 0.5) = 750
        p4 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        p5 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        p6 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        p7 = make_peak_file([narrowpeak_row("chr1", 600, 700)])
        rc, _, peaks3, _ = run_consensus(
            [p4, p5, p6, p7], ["1", "2", "3", "4"])
        if rc != 0 or peaks3.strip().split("\t")[4] != "750":
            print("   3/4 expected 750, got %s" % (
                peaks3.strip().split("\t")[4] if peaks3.strip() else "empty"))
            return False

        return True
    check("T21: score rounding 2/2=1000, 2/3=667, 3/4=750", t21)

    # ===================================================================
    # Edge cases (T22–T25)
    # ===================================================================

    # T22: Empty consensus output
    def t22():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        p2 = make_peak_file([narrowpeak_row("chr2", 100, 200)])
        rc, stderr, peaks, summary = run_consensus([p1, p2], ["1", "2"])
        if rc != 0:
            print("   exit code %d (should be 0 for empty consensus)" % rc)
            return False
        lines = [l for l in peaks.strip().split("\n") if l.strip()]
        if len(lines) != 0:
            print("   expected empty output")
            return False
        # Summary should have count=0 and distribution={}
        sum_lines = summary.strip().split("\n")
        if len(sum_lines) < 2:
            print("   expected header + data row in summary")
            return False
        data = sum_lines[1].split("\t")
        # column 7 is consensus_peak_count (0-indexed)
        if len(data) > 7 and data[7] != "0":
            print("   expected consensus_peak_count=0, got %s" % data[7])
            return False
        # column 8 is support_distribution
        if len(data) > 8 and data[8] != "{}":
            print("   expected support_distribution={}, got %s" % data[8])
            return False
        return True
    check("T22: empty consensus → output empty, summary count=0, dist={}", t22)

    # T23: Malformed line reports physical line number
    def t23():
        lines = [
            "# Header line 1",
            "",  # blank line 2
            narrowpeak_row("chr1", "bad_start", 200),  # line 3 — error
        ]
        p1 = make_peak_file(lines)
        p2 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        rc, stderr, peaks, summary = run_consensus([p1, p2], ["1", "2"])
        if rc == 0:
            print("   expected non-zero exit for malformed line")
            return False
        # Should report line 3 (physical), not line 1 (data)
        if "line 3" not in stderr:
            print("   expected 'line 3' in error, got: %s" % stderr[:200])
            return False
        return True
    check("T23: malformed line → non-zero exit with physical line number", t23)

    # T24: end <= start rejected
    def t24():
        p1 = make_peak_file([narrowpeak_row("chr1", 200, 100)])  # end < start
        p2 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        rc, stderr, peaks, summary = run_consensus([p1, p2], ["1", "2"])
        if rc == 0:
            print("   expected non-zero exit for end <= start")
            return False
        return True
    check("T24: end <= start → non-zero exit", t24)

    # T25: Summit fallback (midpoint used when source summit invalid)
    def t25():
        # Max signal is from p1 (5.0), but its summit is invalid.
        # Fallback to midpoint: (200-100)//2 = 50.
        for bad_peak in (9999, "not_an_int"):
            p1 = make_peak_file([
                narrowpeak_row("chr1", 100, 200, signal=5.0, peak=bad_peak)
            ])
            p2 = make_peak_file([
                narrowpeak_row("chr1", 100, 200, signal=3.0, peak=50)
            ])
            rc, stderr, peaks, summary = run_consensus([p1, p2], ["1", "2"])
            if rc != 0:
                print("   exit code %d" % rc)
                return False
            fields = peaks.strip().split("\t")
            summit = int(fields[9])
            span = int(fields[2]) - int(fields[1])
            if summit != 50:
                print("   expected midpoint summit=50, got %d" % summit)
                return False
            if not (0 <= summit < span):
                print("   summit %d out of [0, %d)" % (summit, span))
                return False
        return True
    check("T25: invalid summit → midpoint fallback", t25)

    # ===================================================================
    # Summary TSV (T26–T27)
    # ===================================================================

    # T26: Summary has all 13 columns
    def t26():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        p2 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        rc, stderr, peaks, summary = run_consensus([p1, p2], ["1", "2"])
        header = summary.split("\n")[0].split("\t")
        expected_cols = [
            "experiment", "assay", "peak_mode", "caller", "n_bioreps",
            "min_replicates", "reciprocal_overlap", "consensus_peak_count",
            "support_distribution", "biorep_labels", "source_peak_files",
            "final_method", "final_output",
        ]
        missing = [c for c in expected_cols if c not in header]
        if missing:
            print("   missing columns: %s" % missing)
            return False
        if len(header) != 13:
            print("   expected 13 columns, got %d" % len(header))
            return False
        return True
    check("T26: summary has all 13 columns", t26)

    # T27: Summary values correct
    def t27():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        p2 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        rc, stderr, peaks, summary = run_consensus(
            [p1, p2], ["A", "B"],
            experiment="myexp", assay="cuttag", caller="macs3",
            peak_mode="narrow", final_method="idr",
            final_output="../06_idr/final/conservative.narrowPeak")
        data = summary.split("\n")[1].split("\t")
        header = summary.split("\n")[0].split("\t")
        idx = {col: i for i, col in enumerate(header)}

        checks = [
            (data[idx["experiment"]] == "myexp", "experiment"),
            (data[idx["assay"]] == "cuttag", "assay"),
            (data[idx["caller"]] == "macs3", "caller"),
            (data[idx["peak_mode"]] == "narrow", "peak_mode"),
            (data[idx["n_bioreps"]] == "2", "n_bioreps"),
            (data[idx["min_replicates"]] == "2", "min_replicates"),
            (data[idx["consensus_peak_count"]] == "1", "consensus_peak_count"),
            (data[idx["final_method"]] == "idr", "final_method"),
            (data[idx["final_output"]] == "../06_idr/final/conservative.narrowPeak", "final_output"),
            (data[idx["biorep_labels"]] == "A,B", "biorep_labels"),
        ]
        for ok, name in checks:
            if not ok:
                print("   %s mismatch" % name)
                return False

        # support_distribution should be valid JSON
        try:
            dist = json.loads(data[idx["support_distribution"]])
            if dist != {"2": 1}:
                print("   expected support_distribution={'2': 1}, got %s" % dist)
                return False
        except json.JSONDecodeError:
            print("   support_distribution not valid JSON: %s" % data[idx["support_distribution"]])
            return False

        # source_peak_files should be valid JSON array
        try:
            sources = json.loads(data[idx["source_peak_files"]])
            if len(sources) != 2:
                print("   expected 2 source files, got %s" % sources)
                return False
        except json.JSONDecodeError:
            print("   source_peak_files not valid JSON")
            return False

        return True
    check("T27: summary values correct (experiment, assay, dist, labels, sources)", t27)

    # ===================================================================
    # Additional validation (T28–T31)
    # ===================================================================

    # T28: Invalid reciprocal_overlap rejected
    def t28():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        p2 = make_peak_file([narrowpeak_row("chr1", 100, 200)])

        # Test 0
        rc, stderr, _, _ = run_consensus([p1, p2], ["1", "2"], overlap=0)
        if rc == 0:
            print("   overlap=0 should be rejected")
            return False

        # Test 1.5
        rc, stderr, _, _ = run_consensus([p1, p2], ["1", "2"], overlap=1.5)
        if rc == 0:
            print("   overlap=1.5 should be rejected")
            return False

        return True
    check("T28: invalid reciprocal_overlap (0, 1.5) rejected", t28)

    # T29: Invalid format rejected
    def t29():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        p2 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        rc, stderr, _, _ = run_consensus([p1, p2], ["1", "2"], fmt="gff")
        if rc == 0:
            print("   invalid format should be rejected")
            return False
        return True
    check("T29: invalid format → non-zero exit", t29)

    # T30: BED3 input accepted, output BED6
    def t30():
        p1 = make_peak_file([bed3_row("chr1", 100, 500)])
        p2 = make_peak_file([bed3_row("chr1", 100, 500)])
        rc, stderr, peaks, summary = run_consensus(
            [p1, p2], ["1", "2"], fmt="bed")
        if rc != 0:
            print("   exit code %d: %s" % (rc, stderr[:200]))
            return False
        lines = [l for l in peaks.strip().split("\n") if l.strip()]
        if len(lines) != 1:
            print("   expected 1 BED peak from BED3 input, got %d" % len(lines))
            return False
        fields = lines[0].split("\t")
        if len(fields) != 6:
            print("   expected BED6 output from BED3 input, got %d cols" % len(fields))
            return False
        # Name should be synthesized (consensus_peak_1)
        if "consensus_peak" not in fields[3]:
            print("   expected synthesized name, got %s" % fields[3])
            return False
        return True
    check("T30: BED3 input → BED6 output with synthesized name", t30)

    # T31: Parent directories created
    def t31():
        p1 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        p2 = make_peak_file([narrowpeak_row("chr1", 100, 200)])
        tmpdir = tempfile.mkdtemp(prefix="consensus_nested_")
        nested_out = os.path.join(tmpdir, "a", "b", "consensus.narrowPeak")
        nested_sum = os.path.join(tmpdir, "a", "b", "summary.tsv")

        cmd = [
            sys.executable, SCRIPT,
            "--peaks", p1, p2,
            "--bioreps", "1", "2",
            "--format", "narrowPeak",
            "--output", nested_out,
            "--summary", nested_sum,
            "--experiment", "test",
            "--assay", "chipseq",
            "--caller", "macs3",
            "--peak-mode", "narrow",
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print("   exit code %d: %s" % (result.returncode, result.stderr[:200]))
            return False
        if not os.path.isfile(nested_out):
            print("   output file not created: %s" % nested_out)
            return False
        if not os.path.isfile(nested_sum):
            print("   summary file not created: %s" % nested_sum)
            return False
        return True
    check("T31: parent directories created for --output and --summary", t31)

    # ===================================================================
    # Summary
    # ===================================================================
    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
