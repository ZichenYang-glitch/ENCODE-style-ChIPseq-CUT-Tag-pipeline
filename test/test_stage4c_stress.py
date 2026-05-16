"""Stage 4c stress tests — tool parameterization validation and dry-run checks.

All tests use temp config/sample files. No real FASTQs or Bowtie2 indexes needed.
Reuses Stage 4b harness pattern: SNAKEMAKE env var with conda fallback,
temp file cleanup, scripts/__pycache__ removal.
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
samples: "test_stage4c_samples.tsv"
use_control: false
threads: 1
"""

# Minimal single-sample TSV (no replicate columns)
SINGLE_HDR = ("sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\t"
              "peak_mode\tgenome\tbowtie2_index\n")
SINGLE_SAMPLE = SINGLE_HDR + "S1\tR1.fq\t\tSE\tchipseq\tT\tnarrow\ths\tidx\n"


def validate(name, config_yaml, samples_tsv, expect_fail=True, expected_error=""):
    """Run validate_samples.py and check result."""
    with open("test_stage4c_config.yaml", "w") as f:
        f.write(config_yaml)
    with open("test_stage4c_samples.tsv", "w") as f:
        f.write(samples_tsv)

    result = subprocess.run(
        [sys.executable, VALIDATOR, "--config", "test_stage4c_config.yaml"],
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

    with open("test_stage4c_config.yaml", "w") as f:
        f.write(config_yaml)
    with open("test_stage4c_samples.tsv", "w") as f:
        f.write(samples_tsv)

    result = subprocess.run(
        [SNAKEMAKE, "-s", SNAKEFILE,
         "--configfile", "test_stage4c_config.yaml",
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
                print("FAIL: %s\n   Rule '%s' found in job listing but should be absent."
                      % (name, rn))
                return False

    print("PASS: %s" % name)
    return True


def dryrun_printshell(name, config_yaml, samples_tsv,
                      expect_texts=None, forbid_texts=None):
    """Dry-run with -p to inspect shell commands. Check for expected/forbidden text."""
    for fq in ["R1.fq", "R2.fq", "R3.fq", "R4.fq"]:
        if not os.path.exists(fq):
            with open(fq, "w") as f:
                f.write("")

    with open("test_stage4c_config.yaml", "w") as f:
        f.write(config_yaml)
    with open("test_stage4c_samples.tsv", "w") as f:
        f.write(samples_tsv)

    result = subprocess.run(
        [SNAKEMAKE, "-s", SNAKEFILE,
         "--configfile", "test_stage4c_config.yaml",
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
    for f in ["test_stage4c_config.yaml", "test_stage4c_samples.tsv",
              "R1.fq", "R2.fq", "R3.fq", "R4.fq"]:
        if os.path.exists(f):
            os.remove(f)
    for d in ["__pycache__", "scripts/__pycache__"]:
        if os.path.isdir(d):
            shutil.rmtree(d)


def main():
    print("Starting Stage 4c (Parameterization Foundation) Stress Tests\n")
    tests = 0
    passed = 0

    # ------------------------------------------------------------------
    # 1. No tool_parameters key → validates and dry-runs
    # ------------------------------------------------------------------
    tests += 1
    if validate("no tool_parameters key validates", BASE_CONFIG, SINGLE_SAMPLE,
                expect_fail=False):
        passed += 1

    tests += 1
    if dryrun("no tool_parameters key dry-runs", BASE_CONFIG, SINGLE_SAMPLE):
        passed += 1

    # ------------------------------------------------------------------
    # 2. Empty tool_parameters block
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "tool_parameters: {}\n"
    tests += 1
    if validate("empty tool_parameters validates", cfg, SINGLE_SAMPLE,
                expect_fail=False):
        passed += 1

    tests += 1
    if dryrun("empty tool_parameters dry-runs", cfg, SINGLE_SAMPLE):
        passed += 1

    # ------------------------------------------------------------------
    # 3. Unknown tool block rejected
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "tool_parameters:\n  bowtie:\n    extra_args: ''\n"
    samples = SINGLE_SAMPLE
    tests += 1
    if validate("unknown tool block rejected", cfg, samples,
                expected_error="unknown tool block"):
        passed += 1

    # ------------------------------------------------------------------
    # 4. Unknown key in known block rejected
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "tool_parameters:\n  trim_galore:\n    adapter: 'foo'\n"
    tests += 1
    if validate("unknown key rejected", cfg, samples,
                expected_error="unknown key"):
        passed += 1

    # ------------------------------------------------------------------
    # 5. Invalid type: quality not int
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "tool_parameters:\n  trim_galore:\n    quality: 'high'\n"
    tests += 1
    if validate("quality not int rejected", cfg, samples,
                expected_error="must be an integer"):
        passed += 1

    # ------------------------------------------------------------------
    # 6. Invalid type: qvalue == 0
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "tool_parameters:\n  macs3:\n    qvalue: 0\n"
    tests += 1
    if validate("qvalue zero rejected", cfg, samples,
                expected_error="must be positive"):
        passed += 1

    # ------------------------------------------------------------------
    # 7. Invalid enum: normalize_using RPGC rejected
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "tool_parameters:\n  bamcoverage:\n    normalize_using: RPGC\n"
    tests += 1
    if validate("RPGC normalize_using rejected", cfg, samples,
                expected_error="RPGC is not supported"):
        passed += 1

    # ------------------------------------------------------------------
    # 8. Invalid enum: mode unknown
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "tool_parameters:\n  bowtie2:\n    mode: extreme\n"
    tests += 1
    if validate("unknown bowtie2 mode rejected", cfg, samples,
                expected_error="must be one of"):
        passed += 1

    # ------------------------------------------------------------------
    # 9. String boolean accepted
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "tool_parameters:\n  bowtie2:\n    dovetail: 'true'\n"
    tests += 1
    if validate("string 'true' accepted as boolean", cfg, SINGLE_SAMPLE,
                expect_fail=False):
        passed += 1

    # ------------------------------------------------------------------
    # 10. Hex filter_flags accepted
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "tool_parameters:\n  samtools_filter:\n    filter_flags: '0xFF'\n"
    tests += 1
    if validate("hex filter_flags accepted", cfg, SINGLE_SAMPLE,
                expect_fail=False):
        passed += 1

    # ------------------------------------------------------------------
    # 11a. filter_flags "0" (decimal string) rejected
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "tool_parameters:\n  samtools_filter:\n    filter_flags: '0'\n"
    tests += 1
    if validate("filter_flags '0' rejected", cfg, SINGLE_SAMPLE,
                expected_error="must be positive"):
        passed += 1

    # ------------------------------------------------------------------
    # 11b. filter_flags "0x0" (hex string) rejected
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "tool_parameters:\n  samtools_filter:\n    filter_flags: '0x0'\n"
    tests += 1
    if validate("filter_flags '0x0' rejected", cfg, SINGLE_SAMPLE,
                expected_error="must be positive"):
        passed += 1

    # ------------------------------------------------------------------
    # 11c. extra_args appear in dry-run -p output
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "tool_parameters:\n  fastqc:\n    extra_args: '--nogroup'\n"
    tests += 1
    if dryrun_printshell("fastqc extra_args in -p output", cfg, SINGLE_SAMPLE,
                         expect_texts=["--nogroup"]):
        passed += 1

    # ------------------------------------------------------------------
    # 12. Default params visible in -p output
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG
    tests += 1
    if dryrun_printshell("default -q 0.01 visible", cfg, SINGLE_SAMPLE,
                         expect_texts=["-q 0.01"]):
        passed += 1

    # ------------------------------------------------------------------
    # 13. Opt-in params absent when not configured
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG
    tests += 1
    if dryrun_printshell("opt-in params absent", cfg, SINGLE_SAMPLE,
                         forbid_texts=["--quality", "--length", "--stringency",
                                       "--very-sensitive"]):
        passed += 1

    # ------------------------------------------------------------------
    # 14. Opt-in params present when configured
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "tool_parameters:\n  trim_galore:\n    quality: 25\n"
    tests += 1
    if dryrun_printshell("quality 25 present when configured", cfg, SINGLE_SAMPLE,
                         expect_texts=["--quality 25"]):
        passed += 1

    # ------------------------------------------------------------------
    # 15. MACS3 extra_args appended after policy args
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "tool_parameters:\n  macs3:\n    extra_args: '--keep-dup all'\n"
    tests += 1
    if dryrun_printshell("macs3 extra_args appended", cfg, SINGLE_SAMPLE,
                         expect_texts=["--keep-dup all"]):
        passed += 1

    # ------------------------------------------------------------------
    # 16. samtools_filter -F flag visible
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG
    tests += 1
    if dryrun_printshell("samtools_filter -F visible", cfg, SINGLE_SAMPLE,
                         expect_texts=["-F 0x904"]):
        passed += 1

    # ------------------------------------------------------------------
    print("\nSummary: %d/%d tests passed." % (passed, tests))
    cleanup()
    if passed < tests:
        sys.exit(1)


if __name__ == "__main__":
    main()
