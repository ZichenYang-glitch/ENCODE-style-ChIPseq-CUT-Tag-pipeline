"""Unit tests for scripts/split_pseudoreps.py.

Creates a tiny SAM fixture, converts to BAM, runs split_pseudoreps.py,
and verifies correctness properties.
"""

import os
import subprocess
import sys
import tempfile

SPLIT_SCRIPT = "scripts/split_pseudoreps.py"
SMT = os.environ.get("SAMTOOLS",
                      "/home/irenadler/miniconda3/envs/chipseq/bin/samtools")


def _run(cmd):
    """Run a command with SAMTOOLS env, exit on failure."""
    env = os.environ.copy()
    env.setdefault("SAMTOOLS", SMT)
    result = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if result.returncode != 0:
        print(f"FAIL: {' '.join(cmd)}")
        print(f"  {result.stderr.strip()[-300:]}")
        sys.exit(1)
    return result


def main():
    print("Starting split_pseudoreps.py unit tests\n")
    tests = 0
    passed = 0

    tmpdir = tempfile.mkdtemp(prefix="split_test_")

    sam_lines = [
        "@HD\tVN:1.6\tSO:coordinate",
        "@SQ\tSN:chr1\tLN:1000",
        "read1/1\t99\tchr1\t1\t60\t100M\t=\t100\t99\t*\t*",
        "read1/2\t147\tchr1\t100\t60\t100M\t=\t1\t-99\t*\t*",
        "read2\t0\tchr1\t500\t60\t100M\t*\t0\t0\t*\t*",
    ]

    sam_path = os.path.join(tmpdir, "test.sam")
    bam_path = os.path.join(tmpdir, "test.bam")
    pr1_bam = os.path.join(tmpdir, "test.pr1.bam")
    pr2_bam = os.path.join(tmpdir, "test.pr2.bam")

    with open(sam_path, "w") as f:
        f.write("\n".join(sam_lines) + "\n")

    # SAM -> BAM
    _run([SMT, "view", "-bS", "-o", bam_path, sam_path])

    # Test 1: split succeeds
    _run([sys.executable, SPLIT_SCRIPT,
          "--input", bam_path, "--out1", pr1_bam, "--out2", pr2_bam,
          "--seed", "42", "--threads", "1"])
    tests += 1
    if os.path.exists(pr1_bam) and os.path.exists(pr2_bam):
        print("PASS: split produces pr1.bam and pr2.bam")
        passed += 1
    else:
        print("FAIL: pr1.bam or pr2.bam missing")

    # Test 2: complementarity via samtools view -c
    n_input = int(_run([SMT, "view", "-c", bam_path]).stdout.strip())
    n_pr1 = int(_run([SMT, "view", "-c", pr1_bam]).stdout.strip())
    n_pr2 = int(_run([SMT, "view", "-c", pr2_bam]).stdout.strip())
    tests += 1
    if n_pr1 + n_pr2 == n_input:
        print(f"PASS: complementarity {n_pr1} + {n_pr2} == {n_input}")
        passed += 1
    else:
        print(f"FAIL: complementarity {n_pr1} + {n_pr2} != {n_input}")

    # Test 3: no read appears in both outputs
    def _read_names(b):
        names = set()
        for line in _run([SMT, "view", b]).stdout.strip().split("\n"):
            if line:
                names.add(line.split("\t")[0])
        return names

    pr1n = _read_names(pr1_bam)
    pr2n = _read_names(pr2_bam)
    overlap = pr1n & pr2n
    tests += 1
    if len(overlap) == 0:
        print("PASS: no read appears in both outputs")
        passed += 1
    else:
        print(f"FAIL: {len(overlap)} reads in both: {overlap}")

    # Test 4: PE mates stay together
    mate1_in_pr1 = "read1/1" in pr1n
    mate1_in_pr2 = "read1/1" in pr2n
    mate2_in_pr1 = "read1/2" in pr1n
    mate2_in_pr2 = "read1/2" in pr2n
    tests += 1
    if (mate1_in_pr1 and mate2_in_pr1) or (mate1_in_pr2 and mate2_in_pr2):
        print("PASS: PE mates stay together")
        passed += 1
    else:
        print(f"FAIL: PE mates separated: "
              f"r1/1 in pr1={mate1_in_pr1} pr2={mate1_in_pr2}, "
              f"r1/2 in pr1={mate2_in_pr1} pr2={mate2_in_pr2}")

    # Test 5: both .bai files produced
    tests += 1
    bai1 = pr1_bam + ".bai"
    bai2 = pr2_bam + ".bai"
    if os.path.isfile(bai1) and os.path.isfile(bai2):
        if os.path.getsize(bai1) > 0 and os.path.getsize(bai2) > 0:
            print("PASS: both .bai files produced")
            passed += 1
        else:
            print("FAIL: .bai files empty")
    else:
        print("FAIL: .bai files missing")

    # Test 6: same seed deterministic (compare samtools view output)
    pr1b = os.path.join(tmpdir, "test_b.pr1.bam")
    pr2b = os.path.join(tmpdir, "test_b.pr2.bam")
    _run([sys.executable, SPLIT_SCRIPT,
          "--input", bam_path, "--out1", pr1b, "--out2", pr2b,
          "--seed", "42", "--threads", "1"])

    v1 = _run([SMT, "view", pr1_bam]).stdout.strip()
    v2 = _run([SMT, "view", pr1b]).stdout.strip()
    tests += 1
    if v1 == v2:
        print("PASS: same seed deterministic")
        passed += 1
    else:
        print(f"FAIL: same seed non-deterministic "
              f"(run1={len(v1.split(chr(10)))} run2={len(v2.split(chr(10)))})")

    # Cleanup
    for f in os.listdir(tmpdir):
        os.unlink(os.path.join(tmpdir, f))
    os.rmdir(tmpdir)

    print(f"\nSummary: {passed}/{tests} tests passed.")
    if passed < tests:
        sys.exit(1)


if __name__ == "__main__":
    main()
