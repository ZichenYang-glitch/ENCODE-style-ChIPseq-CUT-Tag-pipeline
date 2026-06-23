"""Stage 64 CUT&Tag narrow IDR summary script tests.

Verifies the 15-column TSV output, pass/fail logic, and final peak copy.
"""

import os
import shutil
import subprocess
import sys
import tempfile

SUMMARY_SCRIPT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "scripts", "cuttag_idr_summary.py",
)


def _write_peaks(path, count):
    """Write a narrowPeak file with *count* dummy entries."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        for i in range(count):
            f.write(
                "chr1\t%d\t%d\tpeak_%d\t1000\t.\t5.0\t-1\t-1\t255\n"
                % (i * 100, i * 100 + 50, i)
            )


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

    def run_summary(true_count, pooled_count, self1_count, self2_count):
        tmp = tempfile.mkdtemp(prefix="stage64_summary_")
        true_path = os.path.join(tmp, "true.narrowPeak")
        pooled_path = os.path.join(tmp, "pooled.narrowPeak")
        self1_path = os.path.join(tmp, "self1.narrowPeak")
        self2_path = os.path.join(tmp, "self2.narrowPeak")
        tsv_path = os.path.join(tmp, "summary.tsv")
        final_path = os.path.join(tmp, "final.narrowPeak")

        _write_peaks(true_path, true_count)
        _write_peaks(pooled_path, pooled_count)
        _write_peaks(self1_path, self1_count)
        _write_peaks(self2_path, self2_count)

        result = subprocess.run(
            [sys.executable, SUMMARY_SCRIPT,
             "--true-peaks", true_path,
             "--pooled-peaks", pooled_path,
             "--self1-peaks", self1_path,
             "--self2-peaks", self2_path,
             "--experiment", "cuttag_exp",
             "--assay", "cuttag",
             "--caller", "macs3",
             "--peak-mode", "narrow",
             "--bio-rep-a", "1",
             "--bio-rep-b", "2",
             "--final-method", "idr",
             "--final-output", "/results/final/cuttag.idr.narrowPeak",
             "--output-tsv", tsv_path,
             "--output-peak", final_path],
            capture_output=True, text=True,
        )

        if result.returncode != 0:
            print("   Script failed: %s" % result.stderr.strip()[:200])
            return None

        with open(tsv_path) as f:
            tsv_lines = f.read().strip().split("\n")

        with open(true_path) as f:
            true_content = f.read()

        with open(final_path) as f:
            final_content = f.read()

        return (result, tsv_lines, true_content, final_content)

    # S1: header + 1 data row, 15 columns
    def s1():
        r = run_summary(10, 12, 9, 10)
        if r is None:
            return False
        _, lines, _, _ = r
        return len(lines) == 2 and len(lines[0].split("\t")) == 15
    check("S1: 2 lines, 15 columns", s1)

    # S2: values match expected
    def s2():
        r = run_summary(10, 12, 9, 10)
        if r is None:
            return False
        _, lines, _, _ = r
        cols = lines[1].split("\t")
        expected = {
            0: "cuttag_exp", 1: "cuttag", 2: "narrow", 3: "macs3",
            4: "1", 5: "2",
            6: "10", 7: "12", 8: "9", 9: "10",
            10: "1.200", 11: "1.111", 12: "pass",
            13: "idr", 14: "/results/final/cuttag.idr.narrowPeak",
        }
        for idx, val in expected.items():
            if cols[idx] != val:
                print("   Col %d: expected %r, got %r" % (idx, val, cols[idx]))
                return False
        return True
    check("S2: values match", s2)

    # S3: final peak content identical to true-replicate
    def s3():
        r = run_summary(10, 12, 9, 10)
        if r is None:
            return False
        _, _, true_content, final_content = r
        return true_content == final_content
    check("S3: final peak is copy of true-rep", s3)

    # S4: 0/3/0/0 → rescue=inf, self=NA, status=fail
    def s4():
        r = run_summary(0, 3, 0, 0)
        if r is None:
            return False
        _, lines, _, _ = r
        cols = lines[1].split("\t")
        return cols[10] == "inf" and cols[11] == "NA" and cols[12] == "fail"
    check("S4: 0/3/0/0 → inf/NA/fail", s4)

    # S5: 10/10/10/10 → rescue=1.000, self=1.000, pass
    def s5():
        r = run_summary(10, 10, 10, 10)
        if r is None:
            return False
        _, lines, _, _ = r
        cols = lines[1].split("\t")
        return cols[10] == "1.000" and cols[11] == "1.000" and cols[12] == "pass"
    check("S5: 10/10/10/10 → 1.000/1.000/pass", s5)

    # S6: 5/20/5/5 → rescue=4.000, fail
    def s6():
        r = run_summary(5, 20, 5, 5)
        if r is None:
            return False
        _, lines, _, _ = r
        cols = lines[1].split("\t")
        return cols[10] == "4.000" and cols[12] == "fail"
    check("S6: 5/20/5/5 → rescue=4.000/fail", s6)

    # S7: assay column is cuttag
    def s7():
        r = run_summary(10, 12, 9, 10)
        if r is None:
            return False
        _, lines, _, _ = r
        return lines[1].split("\t")[1] == "cuttag"
    check("S7: assay=cuttag", s7)

    # S8: stdout mentions CUT&Tag, not ATAC
    def s8():
        r = run_summary(10, 12, 9, 10)
        if r is None:
            return False
        result, _, _, _ = r
        combined = result.stdout + result.stderr
        return "CUT&Tag" in combined
    check("S8: output mentions CUT&Tag", s8)

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
