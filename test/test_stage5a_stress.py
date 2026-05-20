"""Stage 5a stress tests — IDR config validation and dry-run DAG checks.

All tests use temp config/sample files. No real FASTQs or Bowtie2 indexes needed.
Reuses Stage 4b harness: SNAKEMAKE env var with conda fallback, temp file
cleanup, scripts/__pycache__ removal.
"""

import subprocess
import os
import sys
import shutil
from _tool_resolver import resolve_tool

VALIDATOR = "scripts/validate_samples.py"
SNAKEFILE = "workflow/Snakefile"
SNAKEMAKE = resolve_tool("snakemake", "SNAKEMAKE")

BASE_CONFIG = """\
samples: "test_stage5a_samples.tsv"
use_control: false
threads: 1
stage4b: true
"""


def validate(name, config_yaml, samples_tsv, expect_fail=True, expected_error=""):
    """Run validate_samples.py and check result."""
    with open("test_stage5a_config.yaml", "w") as f:
        f.write(config_yaml)
    with open("test_stage5a_samples.tsv", "w") as f:
        f.write(samples_tsv)

    result = subprocess.run(
        [sys.executable, VALIDATOR, "--config", "test_stage5a_config.yaml"],
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
    """Dry-run and check job listing for expected/forbidden rule names."""
    for fq in ["R1.fq", "R2.fq", "R3.fq", "R4.fq"]:
        if not os.path.exists(fq):
            with open(fq, "w") as f:
                f.write("")

    with open("test_stage5a_config.yaml", "w") as f:
        f.write(config_yaml)
    with open("test_stage5a_samples.tsv", "w") as f:
        f.write(samples_tsv)

    result = subprocess.run(
        [SNAKEMAKE, "-s", SNAKEFILE,
         "--configfile", "test_stage5a_config.yaml",
         "-n"],
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
                print("FAIL: %s\n   Rule '%s' not found in job listing."
                      % (name, rn))
                return False

    if forbid_rules:
        for rn in forbid_rules:
            if rn in output:
                print("FAIL: %s\n   Rule '%s' found but should be absent."
                      % (name, rn))
                return False

    print("PASS: %s" % name)
    return True


def dryrun_printshell(name, config_yaml, samples_tsv,
                      expect_texts=None, forbid_texts=None):
    """Dry-run with -p to inspect shell commands."""
    for fq in ["R1.fq", "R2.fq", "R3.fq", "R4.fq"]:
        if not os.path.exists(fq):
            with open(fq, "w") as f:
                f.write("")

    with open("test_stage5a_config.yaml", "w") as f:
        f.write(config_yaml)
    with open("test_stage5a_samples.tsv", "w") as f:
        f.write(samples_tsv)

    result = subprocess.run(
        [SNAKEMAKE, "-s", SNAKEFILE,
         "--configfile", "test_stage5a_config.yaml",
         "-n", "-p"],
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
    for f in ["test_stage5a_config.yaml", "test_stage5a_samples.tsv",
              "R1.fq", "R2.fq", "R3.fq", "R4.fq"]:
        if os.path.exists(f):
            os.remove(f)
    for d in ["__pycache__", "scripts/__pycache__"]:
        if os.path.isdir(d):
            shutil.rmtree(d)


def main():
    print("Starting Stage 5a IDR Stress Tests\n")
    tests = 0
    passed = 0

    # Common sample sheets
    SINGLE_HDR = ("sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\t"
                  "peak_mode\tgenome\tbowtie2_index\n")
    SINGLE = SINGLE_HDR + "S1\tR1.fq\t\tSE\tchipseq\tT\tnarrow\ths\tidx\n"

    HDR = ("sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\t"
           "peak_mode\tgenome\tbowtie2_index\texperiment\tcondition\t"
           "replicate\tbiological_replicate\ttechnical_replicate\trole\n")

    # 2 chipseq narrow bio-reps (labels 1 and 2)
    TWO_BIOREPS = HDR + (
        "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t1\t1\t1\ttreatment\n"
        "T2\tR2.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t2\t2\t1\ttreatment\n"
    )

    # 2 chipseq narrow bio-reps with arbitrary labels 2 and 4
    TWO_BIOREPS_ARB = HDR + (
        "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t2\t2\t1\ttreatment\n"
        "T2\tR2.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t4\t4\t1\ttreatment\n"
    )

    # ------------------------------------------------------------------
    # 1. No stage5 key → validates, same DAG
    # ------------------------------------------------------------------
    tests += 1
    if validate("no stage5 key validates", BASE_CONFIG, SINGLE, expect_fail=False):
        passed += 1

    tests += 1
    if dryrun("no stage5 key dry-runs", BASE_CONFIG, SINGLE):
        passed += 1

    # ------------------------------------------------------------------
    # 2. stage5: false → validates, same DAG
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "stage5: false\n"
    tests += 1
    if validate("stage5 false validates", cfg, SINGLE, expect_fail=False):
        passed += 1

    tests += 1
    if dryrun("stage5 false dry-runs", cfg, SINGLE):
        passed += 1

    # ------------------------------------------------------------------
    # 3. stage5: true + non-chipseq assay → rejected
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "stage5: true\n"
    samples = HDR + (
        "T1\tR1.fq\t\tSE\tcuttag\tH3K4me3\tnarrow\ths\tidx\texp1\tH3K4me3\t1\t1\t1\ttreatment\n"
        "T2\tR2.fq\t\tSE\tcuttag\tH3K4me3\tnarrow\ths\tidx\texp1\tH3K4me3\t2\t2\t1\ttreatment\n"
    )
    tests += 1
    if validate("non-chipseq rejected", cfg, samples, expected_error="chipseq only"):
        passed += 1

    # ------------------------------------------------------------------
    # 4. stage5: true + broad peak_mode → rejected
    # ------------------------------------------------------------------
    samples = HDR + (
        "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tbroad\ths\tidx\texp1\tH3K27ac\t1\t1\t1\ttreatment\n"
        "T2\tR2.fq\t\tSE\tchipseq\tH3K27ac\tbroad\ths\tidx\texp1\tH3K27ac\t2\t2\t1\ttreatment\n"
    )
    tests += 1
    if validate("broad peak rejected", cfg, samples, expected_error="narrowPeak only"):
        passed += 1

    # ------------------------------------------------------------------
    # 5. stage5: true + 1 bio-rep → rejected
    # ------------------------------------------------------------------
    samples = HDR + (
        "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t1\t1\t1\ttreatment\n"
    )
    tests += 1
    if validate("1 bio-rep rejected", cfg, samples, expected_error="requires exactly 2"):
        passed += 1

    # ------------------------------------------------------------------
    # 6. stage5: true + 3 bio-reps → rejected
    # ------------------------------------------------------------------
    samples = HDR + (
        "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t1\t1\t1\ttreatment\n"
        "T2\tR2.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t2\t2\t1\ttreatment\n"
        "T3\tR3.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t3\t3\t1\ttreatment\n"
    )
    tests += 1
    if validate("3 bio-reps rejected", cfg, samples, expected_error="requires exactly 2"):
        passed += 1

    # ------------------------------------------------------------------
    # 7. Invalid idr.threshold (0, 1.5, "high") → rejected
    # ------------------------------------------------------------------
    cfg_t = BASE_CONFIG + "stage5: true\nidr:\n  threshold: 0\n"
    tests += 1
    if validate("threshold 0 rejected", cfg_t, TWO_BIOREPS, expected_error="must be in (0, 1)"):
        passed += 1

    cfg_t = BASE_CONFIG + "stage5: true\nidr:\n  threshold: 1.5\n"
    tests += 1
    if validate("threshold 1.5 rejected", cfg_t, TWO_BIOREPS, expected_error="must be in (0, 1)"):
        passed += 1

    cfg_t = BASE_CONFIG + "stage5: true\nidr:\n  threshold: high\n"
    tests += 1
    if validate("threshold 'high' rejected", cfg_t, TWO_BIOREPS, expected_error="float in (0, 1)"):
        passed += 1

    # ------------------------------------------------------------------
    # 8. Invalid idr.rank → rejected
    # ------------------------------------------------------------------
    cfg_r = BASE_CONFIG + "stage5: true\nidr:\n  rank: score\n"
    tests += 1
    if validate("invalid rank rejected", cfg_r, TWO_BIOREPS, expected_error="must be 'p.value' or 'signal.value'"):
        passed += 1

    # ------------------------------------------------------------------
    # 9. stage5: true + stage4b: false → rejected
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG.replace("stage4b: true", "stage4b: false") + "stage5: true\n"
    tests += 1
    if validate("stage5 requires stage4b", cfg, TWO_BIOREPS, expected_error="stage5=true requires stage4b=true"):
        passed += 1

    # ------------------------------------------------------------------
    # 10. stage5: true + 2 chipseq narrow bio-reps → validates, IDR rules in DAG
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "stage5: true\n"
    tests += 1
    if validate("2 chipseq narrow validates", cfg, TWO_BIOREPS, expect_fail=False):
        passed += 1

    tests += 1
    if dryrun("IDR rules in DAG", cfg, TWO_BIOREPS,
              expect_rules=["macs3_idr_biorep", "idr_true_replicates"]):
        passed += 1

    # ------------------------------------------------------------------
    # 11. Arbitrary bio_rep labels (2 and 4) → correct paths in DAG
    # ------------------------------------------------------------------
    tests += 1
    if dryrun("arbitrary bio-rep labels", cfg, TWO_BIOREPS_ARB,
              expect_rules=["macs3_idr_biorep"]):
        passed += 1

    tests += 1
    if dryrun_printshell("bio-rep 2 and 4 paths", cfg, TWO_BIOREPS_ARB,
                         expect_texts=["biorep2_idr", "biorep4_idr"]):
        passed += 1

    # ------------------------------------------------------------------
    # 12. idr_macs3.pvalue → -p 0.1 in dry-run -p
    # ------------------------------------------------------------------
    tests += 1
    if dryrun_printshell("-p 0.1 in IDR MACS3", cfg, TWO_BIOREPS,
                         expect_texts=["-p 0.1"]):
        passed += 1

    # ------------------------------------------------------------------
    # 13. idr --idr-threshold in dry-run -p
    # ------------------------------------------------------------------
    tests += 1
    if dryrun_printshell("--idr-threshold visible", cfg, TWO_BIOREPS,
                         expect_texts=["--idr-threshold"]):
        passed += 1

    # ------------------------------------------------------------------
    print("\nSummary: %d/%d tests passed." % (passed, tests))
    cleanup()
    if passed < tests:
        sys.exit(1)


if __name__ == "__main__":
    main()
