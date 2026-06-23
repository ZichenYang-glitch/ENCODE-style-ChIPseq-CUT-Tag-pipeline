"""Stage 65 broad IDR summary script tests.

Verifies 15-column TSV, pass/fail logic, final peak copy, and
the 17-column broadPeak IDR output contract.
"""

import os
import subprocess
import sys
import tempfile

SUMMARY_SCRIPT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "scripts", "broad_idr_summary.py",
)

# Synthetic 17-column broadPeak-like IDR thresholded output
BROADPEAK_17COL_LINE = (
    "chr1\t100\t200\tpeak_1\t1000\t.\t5.0\t3.0\t2.0\t"
    "0.01\t0.005\t100\t200\t5.0\t150\t250\t4.5"
)


def _write_broadpeak(path, count, header_line=None):
    """Write a broadPeak IDR thresholded file with *count* entries.

    Uses the synthetic 17-column broadPeak format from IDR:
    chrom,start,end,name,score,strand,signalValue,pValue,qValue,
    localIDR,globalIDR,rep1_start,rep1_end,rep1_signalValue,
    rep2_start,rep2_end,rep2_signalValue
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        if header_line:
            f.write(header_line + "\n")
        for i in range(count):
            f.write(BROADPEAK_17COL_LINE + "\n")


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

    def run_summary(true_count, pooled_count, self1_count, self2_count,
                    assay="chipseq", header_line=None):
        tmp = tempfile.mkdtemp(prefix="stage65_summary_")
        true_path = os.path.join(tmp, "true.broadPeak")
        pooled_path = os.path.join(tmp, "pooled.broadPeak")
        self1_path = os.path.join(tmp, "self1.broadPeak")
        self2_path = os.path.join(tmp, "self2.broadPeak")
        tsv_path = os.path.join(tmp, "summary.tsv")
        final_path = os.path.join(tmp, "final.broadPeak")

        _write_broadpeak(true_path, true_count, header_line)
        _write_broadpeak(pooled_path, pooled_count, header_line)
        _write_broadpeak(self1_path, self1_count, header_line)
        _write_broadpeak(self2_path, self2_count, header_line)

        result = subprocess.run(
            [sys.executable, SUMMARY_SCRIPT,
             "--true-peaks", true_path,
             "--pooled-peaks", pooled_path,
             "--self1-peaks", self1_path,
             "--self2-peaks", self2_path,
             "--experiment", "broad_exp",
             "--assay", assay,
             "--caller", "macs3",
             "--peak-mode", "broad",
             "--bio-rep-a", "1",
             "--bio-rep-b", "2",
             "--final-method", "idr",
             "--final-output", "/results/final/broad.idr.broadPeak",
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
            1: "chipseq", 2: "broad", 3: "macs3",
            6: "10", 7: "12", 8: "9", 9: "10",
            10: "1.200", 11: "1.111", 12: "pass",
            13: "idr", 14: "/results/final/broad.idr.broadPeak",
        }
        for idx, val in expected.items():
            if cols[idx] != val:
                print("   Col %d: expected %r, got %r" % (idx, val, cols[idx]))
                return False
        return True
    check("S2: values match", s2)

    # S3: final peak is copy of true-rep (preserves all 17 columns)
    def s3():
        r = run_summary(3, 4, 3, 3)
        if r is None:
            return False
        _, _, true_content, final_content = r
        if true_content != final_content:
            print("   Final content differs from true-rep")
            return False
        # Verify 17 columns preserved
        lines = final_content.strip().split("\n")
        for line in lines:
            if len(line.split("\t")) != 17:
                print("   Expected 17 columns, got %d" % len(line.split("\t")))
                return False
        return True
    check("S3: final peak copies all 17 columns", s3)

    # S4: 0/3/0/0 → inf/NA/fail
    def s4():
        r = run_summary(0, 3, 0, 0)
        if r is None:
            return False
        _, lines, _, _ = r
        cols = lines[1].split("\t")
        return cols[10] == "inf" and cols[11] == "NA" and cols[12] == "fail"
    check("S4: 0/3/0/0 → inf/NA/fail", s4)

    # S5: 10/10/10/10 → 1.000/1.000/pass
    def s5():
        r = run_summary(10, 10, 10, 10)
        if r is None:
            return False
        _, lines, _, _ = r
        cols = lines[1].split("\t")
        return cols[10] == "1.000" and cols[11] == "1.000" and cols[12] == "pass"
    check("S5: 10/10/10/10 → 1.000/1.000/pass", s5)

    # S6: header/track lines ignored by count_peaks
    def s6():
        r = run_summary(3, 4, 3, 3, header_line="# IDR thresholded broadPeak output")
        if r is None:
            return False
        _, lines, _, _ = r
        cols = lines[1].split("\t")
        # Only 3 peaks counted despite header line
        return cols[6] == "3"
    check("S6: header lines ignored by count_peaks", s6)

    # S7: assay column is correct for cuttag
    def s7():
        r = run_summary(10, 12, 9, 10, assay="cuttag")
        if r is None:
            return False
        _, lines, _, _ = r
        return lines[1].split("\t")[1] == "cuttag"
    check("S7: assay=cuttag in TSV", s7)

    # S8: output mentions broad-peak IDR
    def s8():
        r = run_summary(10, 12, 9, 10)
        if r is None:
            return False
        result, _, _, _ = r
        return "broad-peak" in (result.stdout + result.stderr).lower()
    check("S8: output mentions broad-peak IDR", s8)

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
