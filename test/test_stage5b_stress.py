"""Stage 5b stress tests — pseudorep IDR validation and dry-run DAG checks.

Reuses Stage 4b/5a harness: SNAKEMAKE env var with conda fallback, temp file
cleanup, scripts/__pycache__ removal.
"""

import subprocess
import os
import sys
import shutil

VALIDATOR = "scripts/validate_samples.py"
SNAKEFILE = "workflow/Snakefile"
SNAKEMAKE = os.environ.get(
    "SNAKEMAKE",
    "/home/irenadler/miniconda3/envs/chipseq/bin/snakemake",
)

BASE_CONFIG = """\
samples: "test_stage5b_samples.tsv"
use_control: false
threads: 1
stage4b: true
stage5: true
"""


def validate(name, config_yaml, samples_tsv, expect_fail=True, expected_error=""):
    with open("test_stage5b_config.yaml", "w") as f:
        f.write(config_yaml)
    with open("test_stage5b_samples.tsv", "w") as f:
        f.write(samples_tsv)
    result = subprocess.run(
        [sys.executable, VALIDATOR, "--config", "test_stage5b_config.yaml"],
        capture_output=True, text=True,
    )
    if expect_fail:
        if result.returncode == 0:
            print("FAIL: %s\n   Expected failure, got pass." % name)
            return False
        elif expected_error not in result.stderr and expected_error not in result.stdout:
            print("FAIL: %s\n   Expected error '%s' not found.\n   stderr: %s"
                  % (name, expected_error, result.stderr.strip()))
            return False
        else:
            print("PASS: %s" % name)
            return True
    else:
        if result.returncode != 0:
            print("FAIL: %s\n   Expected success, got failure.\n   stderr: %s"
                  % (name, result.stderr.strip()))
            return False
        else:
            print("PASS: %s" % name)
            return True


def dryrun(name, config_yaml, samples_tsv,
           expect_rules=None, forbid_rules=None):
    for fq in ["R1.fq", "R2.fq", "R3.fq", "R4.fq"]:
        if not os.path.exists(fq):
            with open(fq, "w") as f:
                f.write("")
    with open("test_stage5b_config.yaml", "w") as f:
        f.write(config_yaml)
    with open("test_stage5b_samples.tsv", "w") as f:
        f.write(samples_tsv)
    result = subprocess.run(
        [SNAKEMAKE, "-s", SNAKEFILE,
         "--configfile", "test_stage5b_config.yaml", "-n"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print("FAIL: %s\n   Dry-run failed."
              "\n   stdout: %s\n   stderr: %s"
              % (name, result.stdout.strip()[-500:],
                 result.stderr.strip()[-500:]))
        return False
    output = result.stdout + result.stderr
    if expect_rules:
        for rn in expect_rules:
            if rn not in output:
                print("FAIL: %s\n   Rule '%s' not found." % (name, rn))
                return False
    if forbid_rules:
        for rn in forbid_rules:
            if rn in output:
                print("FAIL: %s\n   Rule '%s' found but should be absent." % (name, rn))
                return False
    print("PASS: %s" % name)
    return True


def dryrun_printshell(name, config_yaml, samples_tsv,
                      expect_texts=None, forbid_texts=None):
    for fq in ["R1.fq", "R2.fq", "R3.fq", "R4.fq"]:
        if not os.path.exists(fq):
            with open(fq, "w") as f:
                f.write("")
    with open("test_stage5b_config.yaml", "w") as f:
        f.write(config_yaml)
    with open("test_stage5b_samples.tsv", "w") as f:
        f.write(samples_tsv)
    result = subprocess.run(
        [SNAKEMAKE, "-s", SNAKEFILE,
         "--configfile", "test_stage5b_config.yaml", "-n", "-p"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print("FAIL: %s\n   Dry-run -p failed."
              "\n   stdout: %s\n   stderr: %s"
              % (name, result.stdout.strip()[-500:],
                 result.stderr.strip()[-500:]))
        return False
    output = result.stdout + result.stderr
    if expect_texts:
        for t in expect_texts:
            if t not in output:
                print("FAIL: %s\n   Expected text '%s' not found." % (name, t))
                return False
    if forbid_texts:
        for t in forbid_texts:
            if t in output:
                print("FAIL: %s\n   Forbidden text '%s' found." % (name, t))
                return False
    print("PASS: %s" % name)
    return True


def cleanup():
    for f in ["test_stage5b_config.yaml", "test_stage5b_samples.tsv",
              "R1.fq", "R2.fq", "R3.fq", "R4.fq"]:
        if os.path.exists(f):
            os.remove(f)
    for d in ["__pycache__", "scripts/__pycache__"]:
        if os.path.isdir(d):
            shutil.rmtree(d)


def main():
    print("Starting Stage 5b IDR Stress Tests\n")
    tests = 0
    passed = 0

    HDR = ("sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\t"
           "peak_mode\tgenome\tbowtie2_index\texperiment\tcondition\t"
           "replicate\tbiological_replicate\ttechnical_replicate\trole\n")
    TWO_BIOREPS = HDR + (
        "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t1\t1\t1\ttreatment\n"
        "T2\tR2.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t2\t2\t1\ttreatment\n"
    )
    TWO_BIOREPS_ARB = HDR + (
        "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t2\t2\t1\ttreatment\n"
        "T2\tR2.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t4\t4\t1\ttreatment\n"
    )

    # ------------------------------------------------------------------
    # 1. idr.seed positive int validates
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "idr:\n  seed: 123\n"
    tests += 1
    if validate("seed positive int validates", cfg, TWO_BIOREPS, expect_fail=False):
        passed += 1

    # ------------------------------------------------------------------
    # 2. idr.seed negative rejected
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "idr:\n  seed: -5\n"
    tests += 1
    if validate("seed negative rejected", cfg, TWO_BIOREPS, expected_error="must be positive"):
        passed += 1

    # ------------------------------------------------------------------
    # 3. idr.seed zero rejected
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "idr:\n  seed: 0\n"
    tests += 1
    if validate("seed zero rejected", cfg, TWO_BIOREPS, expected_error="must be positive"):
        passed += 1

    # ------------------------------------------------------------------
    # 4. idr.seed non-integer string rejected
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "idr:\n  seed: high\n"
    tests += 1
    if validate("seed non-integer rejected", cfg, TWO_BIOREPS, expected_error="must be a positive integer"):
        passed += 1

    # ------------------------------------------------------------------
    # 5. idr.seed "42" (string) accepted
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "idr:\n  seed: '42'\n"
    tests += 1
    if validate('seed "42" string accepted', cfg, TWO_BIOREPS, expect_fail=False):
        passed += 1

    # ------------------------------------------------------------------
    # 6. split_pseudoreps in DAG
    # ------------------------------------------------------------------
    tests += 1
    if dryrun("split_pseudoreps in DAG", BASE_CONFIG, TWO_BIOREPS,
              expect_rules=["split_pseudoreps"]):
        passed += 1

    # ------------------------------------------------------------------
    # 7. Pseudorep paths use actual bio_rep labels (2 and 4)
    # ------------------------------------------------------------------
    tests += 1
    if dryrun("bio-rep 2 and 4 paths", BASE_CONFIG, TWO_BIOREPS_ARB,
              expect_rules=["split_pseudoreps", "macs3_idr_pseudorep"]):
        passed += 1

    tests += 1
    if dryrun_printshell("biorep2 and biorep4 in paths", BASE_CONFIG,
                         TWO_BIOREPS_ARB,
                         expect_texts=["biorep2", "biorep4"]):
        passed += 1

    # ------------------------------------------------------------------
    # 8. macs3_idr_pseudorep + idr_self_pseudoreps + idr_pooled_pseudoreps in DAG
    # ------------------------------------------------------------------
    tests += 1
    if dryrun("all 5b rules in DAG", BASE_CONFIG, TWO_BIOREPS,
              expect_rules=["macs3_idr_pseudorep", "idr_self_pseudoreps",
                            "idr_pooled_pseudoreps", "stage5b_summary"]):
        passed += 1

    # ------------------------------------------------------------------
    # 9. stage5b_summary in DAG
    # ------------------------------------------------------------------
    tests += 1
    if dryrun("stage5b_summary in DAG", BASE_CONFIG, TWO_BIOREPS,
              expect_rules=["stage5b_summary"]):
        passed += 1

    # ------------------------------------------------------------------
    # 10. split_pseudoreps.py --seed 42 visible
    # ------------------------------------------------------------------
    tests += 1
    if dryrun_printshell("split seed visible", BASE_CONFIG, TWO_BIOREPS,
                         expect_texts=["split_pseudoreps.py", "--seed 42"]):
        passed += 1

    # ------------------------------------------------------------------
    # 11. Self-IDR --idr-threshold visible
    # ------------------------------------------------------------------
    tests += 1
    if dryrun_printshell("self-IDR threshold visible", BASE_CONFIG, TWO_BIOREPS,
                         expect_texts=["--idr-threshold"]):
        passed += 1

    # ------------------------------------------------------------------
    # 12. Pooled-IDR --idr-threshold visible
    # ------------------------------------------------------------------
    tests += 1
    if dryrun_printshell("pooled-IDR threshold visible", BASE_CONFIG, TWO_BIOREPS,
                         expect_texts=["--idr-threshold"]):
        passed += 1

    # ------------------------------------------------------------------
    # 13. Unknown key in idr block rejected
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "idr:\n  seed: 42\n  extra: foo\n"
    tests += 1
    if validate("unknown idr key rejected", cfg, TWO_BIOREPS,
                expected_error="unknown key"):
        passed += 1

    print("\nSummary: %d/%d tests passed." % (passed, tests))
    cleanup()
    if passed < tests:
        sys.exit(1)


if __name__ == "__main__":
    main()
