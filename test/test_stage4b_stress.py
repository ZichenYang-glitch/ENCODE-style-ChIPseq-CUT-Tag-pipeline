"""Stage 4b stress tests — replicate-aware validation and dry-run DAG checks.

All tests use temp config/sample files and dry-run (-n) only.
No real FASTQs or Bowtie2 indexes are required.
"""

import subprocess
import os
import sys
import tempfile
from _tool_resolver import resolve_tool

VALIDATOR = "scripts/validate_samples.py"
SNAKEFILE = "workflow/Snakefile"
SNAKEMAKE = resolve_tool("snakemake", "SNAKEMAKE")

# Base config snippet (no controls, no genome_resources)
BASE_CONFIG = """\
samples: "test_stage4b_samples.tsv"
use_control: false
threads: 1
"""


def validate(name, config_yaml, samples_tsv, expect_fail=True, expected_error=""):
    """Run validate_samples.py and check result. Returns True on pass."""
    with open("test_stage4b_config.yaml", "w") as f:
        f.write(config_yaml)
    with open("test_stage4b_samples.tsv", "w") as f:
        f.write(samples_tsv)

    result = subprocess.run(
        [sys.executable, VALIDATOR, "--config", "test_stage4b_config.yaml"],
        capture_output=True, text=True,
    )

    if expect_fail:
        if result.returncode == 0:
            print("FAIL: %s\n   Expected failure, but got pass." % name)
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
            print("FAIL: %s\n   Expected success, but got failure.\n   stderr: %s"
                  % (name, result.stderr.strip()))
            return False
        else:
            print("PASS: %s" % name)
            return True


def dryrun(name, config_yaml, samples_tsv,
           expect_rules=None, forbid_rules=None):
    """Run snakemake -n and check job listing for rule presence/absence.

    expect_rules: list of rule names that MUST appear in job listing.
    forbid_rules: list of rule names that MUST NOT appear in job listing.
    Returns True on pass.
    """
    # Create dummy FASTQ files so Snakemake input checks pass during -n
    for fq in ["R1.fq", "R2.fq", "R3.fq", "R4.fq"]:
        if not os.path.exists(fq):
            with open(fq, "w") as f:
                f.write("")

    with open("test_stage4b_config.yaml", "w") as f:
        f.write(config_yaml)
    with open("test_stage4b_samples.tsv", "w") as f:
        f.write(samples_tsv)

    result = subprocess.run(
        [SNAKEMAKE, "-s", SNAKEFILE,
         "--configfile", "test_stage4b_config.yaml",
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


def cleanup():
    for f in ["test_stage4b_config.yaml", "test_stage4b_samples.tsv",
              "R1.fq", "R2.fq", "R3.fq", "R4.fq"]:
        if os.path.exists(f):
            os.remove(f)
    import shutil
    for d in ["__pycache__", "scripts/__pycache__"]:
        if os.path.isdir(d):
            shutil.rmtree(d)


def main():
    print("Starting Stage 4b Stress Tests\n")
    tests = 0
    passed = 0

    # ---- Header for samples (includes all columns) ----
    HDR = ("sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\t"
           "peak_mode\tgenome\tbowtie2_index\texperiment\tcondition\t"
           "replicate\tbiological_replicate\ttechnical_replicate\trole"
           "\tcontrol_sample\tcontrol_bam\n")

    SINGLE = HDR + (
        "S1\tR1.fq\t\tSE\tchipseq\tT\tnarrow\ths\tidx\t\t\t\t\t\t\t\t\n"
    )

    TWO_TREATMENT = HDR + (
        "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t1\t1\t1\ttreatment\t\t\n"
        "T2\tR2.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t2\t2\t1\ttreatment\t\t\n"
    )

    # Create a temp file for control_bam path tests
    tmp_ctrl_bam = tempfile.NamedTemporaryFile(suffix=".bam", delete=False)
    tmp_ctrl_bam.close()
    ctrl_bam_path = tmp_ctrl_bam.name

    # ------------------------------------------------------------------
    # 1. Baseline: minimal single-sample validates
    # ------------------------------------------------------------------
    tests += 1
    if validate("single-sample validates", BASE_CONFIG, SINGLE, expect_fail=False):
        passed += 1

    # ------------------------------------------------------------------
    # 2. Baseline: minimal single-sample dry-runs
    # ------------------------------------------------------------------
    tests += 1
    if dryrun("single-sample dry-runs", BASE_CONFIG, SINGLE):
        passed += 1

    # ------------------------------------------------------------------
    # 3. Two-replicate treatment validates
    # ------------------------------------------------------------------
    tests += 1
    if validate("two-replicate treatment validates", BASE_CONFIG, TWO_TREATMENT, expect_fail=False):
        passed += 1

    # ------------------------------------------------------------------
    # 4. Inconsistent assay rejected
    # ------------------------------------------------------------------
    samples = HDR + (
        "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t1\t1\t1\ttreatment\t\t\n"
        "T2\tR2.fq\t\tSE\tcuttag\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t2\t2\t1\ttreatment\t\t\n"
    )
    tests += 1
    if validate("inconsistent assay rejected", BASE_CONFIG, samples, expected_error="disagree on assay"):
        passed += 1

    # ------------------------------------------------------------------
    # 5. Inconsistent genome rejected
    # ------------------------------------------------------------------
    samples = HDR + (
        "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t1\t1\t1\ttreatment\t\t\n"
        "T2\tR2.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\tmm\tidx\texp1\tH3K27ac\t2\t2\t1\ttreatment\t\t\n"
    )
    tests += 1
    if validate("inconsistent genome rejected", BASE_CONFIG, samples, expected_error="disagree on genome"):
        passed += 1

    # ------------------------------------------------------------------
    # 6. Duplicate (bio_rep, tech_rep) rejected
    # ------------------------------------------------------------------
    samples = HDR + (
        "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t1\t1\t1\ttreatment\t\t\n"
        "T2\tR2.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t2\t1\t1\ttreatment\t\t\n"
    )
    tests += 1
    if validate("duplicate (bio_rep, tech_rep) rejected", BASE_CONFIG, samples, expected_error="duplicate"):
        passed += 1

    # ------------------------------------------------------------------
    # 7. Mixed control types rejected
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "use_control: true\n"
    samples = HDR + (
        "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t1\t1\t1\ttreatment\t\t%s\n"
        "T2\tR2.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t2\t2\t1\ttreatment\tCtrl1\t\n"
        "Ctrl1\tR3.fq\t\tSE\tchipseq\tInput\tnarrow\ths\tidx\texp1\tInput\t1\t1\t1\tcontrol\t\t\n"
    ) % ctrl_bam_path
    tests += 1
    if validate("mixed control types rejected", cfg, samples, expected_error="mixed control types"):
        passed += 1

    # ------------------------------------------------------------------
    # 8. Partial controls rejected
    # ------------------------------------------------------------------
    samples = HDR + (
        "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t1\t1\t1\ttreatment\tCtrl1\t\n"
        "T2\tR2.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t2\t2\t1\ttreatment\t\t\n"
        "Ctrl1\tR3.fq\t\tSE\tchipseq\tInput\tnarrow\ths\tidx\texp1\tInput\t1\t1\t1\tcontrol\t\t\n"
    )
    tests += 1
    if validate("partial controls rejected", cfg, samples, expected_error="partial controls"):
        passed += 1

    # ------------------------------------------------------------------
    # 9. Two-replicate dry-run shows pooled rules in jobs
    # ------------------------------------------------------------------
    tests += 1
    if dryrun("two-replicate shows pooled rules", BASE_CONFIG, TWO_TREATMENT,
              expect_rules=["pool_treatment_bam", "macs3_pooled_peaks"]):
        passed += 1

    # ------------------------------------------------------------------
    # 10. stage4b: false hides pooled rules from DAG
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "stage4b: false\n"
    tests += 1
    if dryrun("stage4b:false hides replicate rules from DAG", cfg, TWO_TREATMENT,
              forbid_rules=["pool_treatment_bam", "macs3_pooled_peaks"]):
        passed += 1

    # ------------------------------------------------------------------
    # 11. stage4b: true shows pooled rules in DAG
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "stage4b: true\n"
    tests += 1
    if dryrun("stage4b:true shows replicate rules in DAG", cfg, TWO_TREATMENT,
              expect_rules=["pool_treatment_bam", "macs3_pooled_peaks"]):
        passed += 1

    # ------------------------------------------------------------------
    # 12. Two-replicate + two-control validates
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "use_control: true\n"
    samples = HDR + (
        "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t1\t1\t1\ttreatment\tCtrl1\t\n"
        "T2\tR2.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t2\t2\t1\ttreatment\tCtrl2\t\n"
        "Ctrl1\tR3.fq\t\tSE\tchipseq\tInput\tnarrow\ths\tidx\texp1\tInput\t1\t1\t1\tcontrol\t\t\n"
        "Ctrl2\tR4.fq\t\tSE\tchipseq\tInput\tnarrow\ths\tidx\texp1\tInput\t2\t2\t1\tcontrol\t\t\n"
    )
    tests += 1
    if validate("two-rep + two-control validates", cfg, samples, expect_fail=False):
        passed += 1

    # ------------------------------------------------------------------
    # 13. Single-biorep with tech-reps: bio-rep BAM but no pooled peaks
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "stage4b: true\n"
    samples = HDR + (
        "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t1\t1\t1\ttreatment\t\t\n"
        "T2\tR2.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t1\t1\t2\ttreatment\t\t\n"
    )
    tests += 1
    if dryrun("single-biorep: merge_bam but no pooled peaks", cfg, samples,
              expect_rules=["merge_biorep_bam"],
              forbid_rules=["pool_treatment_bam", "macs3_pooled_peaks"]):
        passed += 1

    # ------------------------------------------------------------------
    # 14. control_sample from wrong experiment rejected
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "use_control: true\n"
    samples = HDR + (
        "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t1\t1\t1\ttreatment\tCtrl1\t\n"
        "T2\tR2.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\texp1\tH3K27ac\t2\t2\t1\ttreatment\tCtrl1\t\n"
        "Ctrl1\tR3.fq\t\tSE\tchipseq\tInput\tnarrow\ths\tidx\texp2\tInput\t1\t1\t1\tcontrol\t\t\n"
    )
    tests += 1
    if validate("control_sample wrong experiment", cfg, samples, expected_error="not"):
        passed += 1

    # ------------------------------------------------------------------
    # 15. stage4b invalid value rejected
    # ------------------------------------------------------------------
    cfg = BASE_CONFIG + "stage4b: maybe\n"
    tests += 1
    if validate("stage4b invalid value rejected", cfg, SINGLE, expected_error="stage4b must be true or false"):
        passed += 1

    # ------------------------------------------------------------------
    print("\nSummary: %d/%d tests passed." % (passed, tests))
    cleanup()
    os.unlink(ctrl_bam_path)
    if passed < tests:
        sys.exit(1)


if __name__ == "__main__":
    main()
