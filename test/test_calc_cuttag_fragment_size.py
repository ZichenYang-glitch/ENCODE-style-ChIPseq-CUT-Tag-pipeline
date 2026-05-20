"""Unit tests for scripts/calc_cuttag_fragment_size.py."""

import os
import subprocess
import sys
import tempfile
from _tool_resolver import resolve_tool

SMT = resolve_tool("samtools", "SAMTOOLS")
SCRIPT = "scripts/calc_cuttag_fragment_size.py"


def _run(cmd):
    env = os.environ.copy()
    env.setdefault("SAMTOOLS", SMT)
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if r.returncode != 0:
        print(f"FAIL: {' '.join(cmd)}\n  {r.stderr.strip()[-300:]}")
        sys.exit(1)
    return r


def _make_bam(tmpdir, sam_lines, bam_name="test.bam"):
    sam = os.path.join(tmpdir, "test.sam")
    bam = os.path.join(tmpdir, bam_name)
    with open(sam, "w") as f:
        f.write("\n".join(sam_lines) + "\n")
    _run([SMT, "view", "-bS", "-o", bam, sam])
    return bam


def _run_script(sample, bam, layout, out):
    _run([sys.executable, SCRIPT,
          "--sample", sample, "--bam", bam,
          "--layout", layout, "--output", out])


def _read_tsv(path):
    lines = open(path).read().strip().split("\n")
    return dict(zip(lines[0].split("\t"), lines[1].split("\t")))


def main():
    print("Starting calc_cuttag_fragment_size tests\n")
    tests = 0
    passed = 0

    # --- Test 1: PE fixture, 5 properly paired pairs ---
    with tempfile.TemporaryDirectory(prefix="cuttag_test_") as td:
        hdr = ["@HD\tVN:1.6\tSO:coordinate", "@SQ\tSN:chr1\tLN:1000"]
        reads = []
        for i in range(5):
            pos = (i + 1) * 100
            tlen = [100, 200, 350, 450, 600][i]
            reads.append(
                f"r{i}/1\t99\tchr1\t{pos}\t60\t100M\t=\t{pos + tlen}\t{tlen}\t*\t*")
            reads.append(
                f"r{i}/2\t147\tchr1\t{pos + tlen}\t60\t100M\t=\t{pos}\t{-tlen}\t*\t*")
        bam = _make_bam(td, hdr + reads)
        out = os.path.join(td, "out.tsv")
        _run_script("s1", bam, "PE", out)
        d = _read_tsv(out)

        tests += 1
        if d["fragment_count"] == "5":
            print("PASS: PE fragment_count=5")
            passed += 1
        else:
            print("FAIL: fragment_count=%s" % d["fragment_count"])

        tests += 1
        # 120,200,350,450,600 -> mean=344.0, median=350.0
        if d["fragment_median"] == "350.0":
            print("PASS: PE median=350.0")
            passed += 1
        else:
            print("FAIL: median=%s" % d["fragment_median"])

        tests += 1
        # Bin check: lt_150=1/5=0.2000, 150_300=1/5=0.2000,
        # 300_500=2/5=0.4000, ge_500=1/5=0.2000 -> sum=1.0000
        s = float(d["fraction_lt_150"]) + float(d["fraction_150_300"]) + \
            float(d["fraction_300_500"]) + float(d["fraction_ge_500"])
        if abs(s - 1.0) <= 0.0002:
            print("PASS: bin sum=%.4f (tol OK)" % s)
            passed += 1
        else:
            print("FAIL: bin sum=%.4f" % s)

        tests += 1
        # lt_120 subset: only the 120bp fragment qualifies
        if d["fraction_lt_120"] == "0.2000":
            print("PASS: lt_120=0.2000 (diagnostic subset)")
            passed += 1
        else:
            print("FAIL: lt_120=%s" % d["fraction_lt_120"])

    # --- Test 2: SE -> layout_not_supported ---
    with tempfile.TemporaryDirectory(prefix="cuttag_test_") as td:
        bam = _make_bam(td, ["@HD\tVN:1.6", "@SQ\tSN:chr1\tLN:1000"])
        out = os.path.join(td, "out.tsv")
        _run_script("s1", bam, "SE", out)
        d = _read_tsv(out)

        tests += 1
        if d["status"] == "layout_not_supported" and d["fragment_count"] == "NA":
            print("PASS: SE -> layout_not_supported, all NA")
            passed += 1
        else:
            print("FAIL: SE status=%s count=%s" % (
                d["status"], d["fragment_count"]))

    # --- Test 3: zero fragments ---
    with tempfile.TemporaryDirectory(prefix="cuttag_test_") as td:
        reads = [
            "r0\t77\t*\t0\t0\t*\t*\t0\t0\t*\t*",   # unmapped
        ]
        bam = _make_bam(td, ["@HD\tVN:1.6", "@SQ\tSN:chr1\tLN:1000"] + reads)
        out = os.path.join(td, "out.tsv")
        _run_script("s1", bam, "PE", out)
        d = _read_tsv(out)

        tests += 1
        if d["status"] == "no_fragments" and d["fragment_count"] == "0":
            print("PASS: zero fragments -> no_fragments, count=0")
            passed += 1
        else:
            print("FAIL: status=%s count=%s" % (d["status"], d["fragment_count"]))

    # --- Test 4: mixed records (mapped + secondary) ---
    with tempfile.TemporaryDirectory(prefix="cuttag_test_") as td:
        # r0: proper pair, r1: secondary alignment, r2: proper pair
        reads = [
            "r0/1\t99\tchr1\t100\t60\t100M\t=\t250\t150\t*\t*",
            "r0/2\t147\tchr1\t250\t60\t100M\t=\t100\t-150\t*\t*",
            "r0/1\t355\tchr1\t100\t0\t100M\t*\t0\t0\t*\t*",    # secondary
            "r2/1\t99\tchr1\t500\t60\t100M\t=\t600\t100\t*\t*",
            "r2/2\t147\tchr1\t600\t60\t100M\t=\t500\t-100\t*\t*",
        ]
        bam = _make_bam(td, ["@HD\tVN:1.6", "@SQ\tSN:chr1\tLN:1000"] + reads)
        out = os.path.join(td, "out.tsv")
        _run_script("s1", bam, "PE", out)
        d = _read_tsv(out)

        tests += 1
        # Only r0/1 and r2/1 are counted (2 proper pairs)
        if d["fragment_count"] == "2":
            print("PASS: mixed records -> count=2 (secondary excluded)")
            passed += 1
        else:
            print("FAIL: mixed count=%s" % d["fragment_count"])

    # --- Test 5: read2 exclusion ---
    with tempfile.TemporaryDirectory(prefix="cuttag_test_") as td:
        reads = [
            "r0/1\t99\tchr1\t100\t60\t100M\t=\t250\t150\t*\t*",
            "r0/2\t147\tchr1\t250\t60\t100M\t=\t100\t-150\t*\t*",
        ]
        bam = _make_bam(td, ["@HD\tVN:1.6", "@SQ\tSN:chr1\tLN:1000"] + reads)
        out = os.path.join(td, "out.tsv")
        _run_script("s1", bam, "PE", out)
        d = _read_tsv(out)

        tests += 1
        if d["fragment_count"] == "1":
            print("PASS: read2 excluded, count=1")
            passed += 1
        else:
            print("FAIL: read2 exclusion count=%s" % d["fragment_count"])

    # --- Test 6: exclusive bin sum with tolerance ---
    with tempfile.TemporaryDirectory(prefix="cuttag_test_") as td:
        reads = []
        tlens = [50, 120, 200, 250, 350, 450, 550, 650]
        for i, t in enumerate(tlens):
            pos = (i + 1) * 200
            reads.append(
                f"r{i}/1\t99\tchr1\t{pos}\t60\t100M\t=\t{pos + t}\t{t}\t*\t*")
            reads.append(
                f"r{i}/2\t147\tchr1\t{pos + t}\t60\t100M\t=\t{pos}\t{-t}\t*\t*")
        bam = _make_bam(td, ["@HD\tVN:1.6", "@SQ\tSN:chr1\tLN:10000"] + reads)
        out = os.path.join(td, "out.tsv")
        _run_script("s1", bam, "PE", out)
        d = _read_tsv(out)

        s = float(d["fraction_lt_150"]) + float(d["fraction_150_300"]) + \
            float(d["fraction_300_500"]) + float(d["fraction_ge_500"])
        lt120 = float(d["fraction_lt_120"])
        lt150 = float(d["fraction_lt_150"])

        tests += 1
        if abs(s - 1.0) <= 0.0002 and lt120 <= lt150:
            print("PASS: bin sum=%.4f, lt_120(%.4f) ≤ lt_150(%.4f)" % (
                s, lt120, lt150))
            passed += 1
        else:
            print("FAIL: sum=%.4f lt120=%.4f lt150=%.4f" % (s, lt120, lt150))

    print("\nSummary: %d/%d tests passed." % (passed, tests))
    if passed < tests:
        sys.exit(1)


if __name__ == "__main__":
    main()
