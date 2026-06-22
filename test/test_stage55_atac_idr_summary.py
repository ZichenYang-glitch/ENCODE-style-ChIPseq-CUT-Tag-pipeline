"""Stage 55 ATAC IDR summary script tests."""

import os
import subprocess
import sys
import tempfile


SCRIPT = "scripts/atac_idr_summary.py"


def _write_peaks(path, count):
    with open(path, "w") as fh:
        for idx in range(count):
            start = 100 + idx * 100
            end = start + 50
            fh.write(
                "chr1\t%d\t%d\tpeak%d\t100\t.\t10\t-1\t-1\t25\n"
                % (start, end, idx + 1)
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
        except Exception as exc:
            print("FAIL: %s\n   Unexpected error: %s" % (name, exc))

    def run_summary(true_count=10, pooled_count=12, self1_count=9, self2_count=10):
        with tempfile.TemporaryDirectory(prefix="stage55_summary_") as tmp:
            true_peaks = os.path.join(tmp, "true.narrowPeak")
            pooled_peaks = os.path.join(tmp, "pooled.narrowPeak")
            self1_peaks = os.path.join(tmp, "self1.narrowPeak")
            self2_peaks = os.path.join(tmp, "self2.narrowPeak")
            final_peak = os.path.join(tmp, "final.narrowPeak")
            summary = os.path.join(tmp, "summary.tsv")

            _write_peaks(true_peaks, true_count)
            _write_peaks(pooled_peaks, pooled_count)
            _write_peaks(self1_peaks, self1_count)
            _write_peaks(self2_peaks, self2_count)

            result = subprocess.run(
                [
                    sys.executable,
                    SCRIPT,
                    "--true-peaks",
                    true_peaks,
                    "--pooled-peaks",
                    pooled_peaks,
                    "--self1-peaks",
                    self1_peaks,
                    "--self2-peaks",
                    self2_peaks,
                    "--experiment",
                    "atac_exp",
                    "--assay",
                    "atac",
                    "--caller",
                    "macs3",
                    "--peak-mode",
                    "narrow",
                    "--bio-rep-a",
                    "1",
                    "--bio-rep-b",
                    "2",
                    "--final-method",
                    "idr",
                    "--final-output",
                    final_peak,
                    "--output-tsv",
                    summary,
                    "--output-peak",
                    final_peak,
                ],
                capture_output=True,
                text=True,
            )
            if result.returncode != 0:
                return result, [], [], ""
            with open(summary) as fh:
                lines = [line.rstrip("\n").split("\t") for line in fh]
            with open(final_peak) as fh:
                final_content = fh.read()
            with open(true_peaks) as fh:
                true_content = fh.read()
            return result, lines, true_content, final_content

    def t1():
        result, lines, _, _ = run_summary()
        if result.returncode != 0:
            print("   script failed: %s" % result.stderr)
            return False
        if len(lines) != 2:
            print("   expected header + row, got %d lines" % len(lines))
            return False
        if len(lines[0]) != 15 or len(lines[1]) != 15:
            print("   expected 15 columns, got %d/%d" % (len(lines[0]), len(lines[1])))
            return False
        return True
    check("S1: writes 15-column summary TSV", t1)

    def t2():
        result, lines, _, _ = run_summary()
        row = dict(zip(lines[0], lines[1]))
        expected = {
            "experiment": "atac_exp",
            "assay": "atac",
            "peak_mode": "narrow",
            "caller": "macs3",
            "bio_rep_a": "1",
            "bio_rep_b": "2",
            "true_peaks_Nt": "10",
            "pooled_peaks_Np": "12",
            "self1_peaks_N1": "9",
            "self2_peaks_N2": "10",
            "rescue_ratio": "1.200",
            "self_consistency_ratio": "1.111",
            "reproducibility_status": "pass",
            "final_method": "idr",
        }
        mismatches = [
            (key, row.get(key), value)
            for key, value in expected.items()
            if row.get(key) != value
        ]
        if mismatches:
            print("   mismatches: %s" % mismatches)
            return False
        return result.returncode == 0
    check("S2: summary values and pass status", t2)

    def t3():
        result, _, true_content, final_content = run_summary()
        if result.returncode != 0:
            print("   script failed")
            return False
        if final_content != true_content:
            print("   final peak file differs from true-replicate peaks")
            return False
        return True
    check("S3: final peak copied from true-replicate IDR", t3)

    def t4():
        result, lines, _, _ = run_summary(true_count=0, pooled_count=3, self1_count=0, self2_count=0)
        row = dict(zip(lines[0], lines[1]))
        if row.get("rescue_ratio") != "inf":
            print("   expected rescue_ratio=inf, got %s" % row.get("rescue_ratio"))
            return False
        if row.get("self_consistency_ratio") != "NA":
            print("   expected self_consistency_ratio=NA, got %s" % row.get("self_consistency_ratio"))
            return False
        if row.get("reproducibility_status") != "fail":
            print("   expected fail status, got %s" % row.get("reproducibility_status"))
            return False
        return result.returncode == 0
    check("S4: zero denominators handled", t4)

    print("\n%d/%d tests passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
