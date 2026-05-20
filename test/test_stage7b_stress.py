"""Stage 7b stress tests — SEACR sidecar scheduling and validation."""

import os
import shutil
import subprocess
import sys
from _tool_resolver import resolve_tool

SNAKEFILE = "workflow/Snakefile"
SNAKEMAKE = resolve_tool("snakemake", "SNAKEMAKE")

BASE_CONFIG = """\
samples: "test_stage7b_samples.tsv"
use_control: false
threads: 1
stage4b: true
stage5: false
"""

SEACR_ON = """\
cuttag:
  seacr:
    enabled: true
"""

SEACR_OFF = """\
cuttag:
  seacr:
    enabled: false
"""

HDR = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\t"
    "peak_mode\tgenome\tbowtie2_index\texperiment\tcondition\t"
    "replicate\tbiological_replicate\ttechnical_replicate\trole\n"
)

CHIPSEQ = HDR + (
    "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\t"
    "exp1\tH3K27ac\t1\t1\t1\ttreatment\n"
)

CUTTAG_PE = HDR + (
    "T1\tR1.fq\tR2.fq\tPE\tcuttag\tH3K4me3\tnarrow\ths\tidx\t"
    "exp1\tH3K4me3\t1\t1\t1\ttreatment\n"
)

CUTTAG_SE = HDR + (
    "T1\tR1.fq\t\tSE\tcuttag\tCTCF\tnarrow\ths\tidx\t"
    "exp1\tCTCF\t1\t1\t1\ttreatment\n"
)


def _write_inputs(config_yaml, samples_tsv):
    for fq in ["R1.fq", "R2.fq"]:
        if not os.path.exists(fq):
            with open(fq, "w") as f:
                f.write("")
    with open("test_stage7b_config.yaml", "w") as f:
        f.write(config_yaml)
    with open("test_stage7b_samples.tsv", "w") as f:
        f.write(samples_tsv)


def _snakemake(args, config_yaml, samples_tsv):
    _write_inputs(config_yaml, samples_tsv)
    return subprocess.run(
        [SNAKEMAKE, "-s", SNAKEFILE, "--configfile",
         "test_stage7b_config.yaml", *args],
        capture_output=True, text=True,
    )


def dryrun(name, cfg, tsv, expect=(), forbid=()):
    r = _snakemake(["-n"], cfg, tsv)
    o = r.stdout + r.stderr
    if r.returncode != 0:
        print("FAIL: %s\n   Dry-run failed.\n   %s" % (
            name, o.strip()[-400:]))
        return False
    for t in expect:
        if t not in o:
            print("FAIL: %s\n   Expected %r not found." % (name, t))
            return False
    for t in forbid:
        if t in o:
            print("FAIL: %s\n   Forbidden %r found." % (name, t))
            return False
    print("PASS: %s" % name)
    return True


def validate(name, cfg, tsv, expect_fail=True, expected_error=""):
    """Run validate_samples.py standalone."""
    _write_inputs(cfg, tsv)
    r = subprocess.run(
        [sys.executable, "scripts/validate_samples.py",
         "--config", "test_stage7b_config.yaml"],
        capture_output=True, text=True,
    )
    if expect_fail:
        if r.returncode == 0:
            print("FAIL: %s\n   Expected failure, got pass." % name)
            return False
        if expected_error not in r.stderr and expected_error not in r.stdout:
            print("FAIL: %s\n   Expected error %r not found.\n   %s" % (
                name, expected_error, r.stderr.strip()))
            return False
        print("PASS: %s" % name)
        return True
    else:
        if r.returncode != 0:
            print("FAIL: %s\n   Expected pass, got:\n   %s" % (
                name, r.stderr.strip()))
            return False
        print("PASS: %s" % name)
        return True


def cleanup():
    for f in ["test_stage7b_config.yaml", "test_stage7b_samples.tsv",
              "R1.fq", "R2.fq"]:
        if os.path.exists(f):
            os.remove(f)
    for d in ["__pycache__", "scripts/__pycache__"]:
        if os.path.isdir(d):
            shutil.rmtree(d)


def main():
    print("Starting Stage 7b Stress Tests\n")
    tests = 0
    passed = 0

    # 1. No cuttag block -> validates, no SEACR
    tests += 1
    if dryrun("no cuttag block", BASE_CONFIG, CUTTAG_PE,
              forbid=("seacr_bedgraph", "seacr_call")):
        passed += 1

    # 2. seacr.enabled: false -> no SEACR
    tests += 1
    if dryrun("seacr disabled", BASE_CONFIG + SEACR_OFF, CUTTAG_PE,
              forbid=("seacr_bedgraph", "seacr_call")):
        passed += 1

    # 3. seacr.enabled: true + CUT&Tag PE -> SEACR in DAG
    tests += 1
    if dryrun("seacr enabled PE", BASE_CONFIG + SEACR_ON, CUTTAG_PE,
              expect=("seacr_bedgraph", "seacr_call")):
        passed += 1

    # 4. seacr.enabled: true + ChIP-seq -> no SEACR
    tests += 1
    if dryrun("seacr ChIP-seq excluded", BASE_CONFIG + SEACR_ON, CHIPSEQ,
              forbid=("seacr_bedgraph", "seacr_call")):
        passed += 1

    # 5. seacr.enabled: true + CUT&Tag SE -> no SEACR
    tests += 1
    if dryrun("seacr SE excluded", BASE_CONFIG + SEACR_ON, CUTTAG_SE,
              forbid=("seacr_bedgraph", "seacr_call")):
        passed += 1

    # 6. Invalid mode rejected
    cfg = BASE_CONFIG + "cuttag:\n  seacr:\n    enabled: true\n    mode: extreme\n"
    tests += 1
    if validate("invalid mode", cfg, CUTTAG_PE,
                expected_error="must be 'stringent' or 'relaxed'"):
        passed += 1

    # 7. Invalid normalization rejected
    cfg = BASE_CONFIG + "cuttag:\n  seacr:\n    enabled: true\n    normalization: nope\n"
    tests += 1
    if validate("invalid norm", cfg, CUTTAG_PE,
                expected_error="must be 'non'"):
        passed += 1

    # 8. threshold:0 rejected
    cfg = BASE_CONFIG + SEACR_ON + "\n    threshold: 0\n"
    tests += 1
    if validate("threshold 0", cfg, CUTTAG_PE,
                expected_error="must be in (0, 1)"):
        passed += 1

    # 9. threshold:1 rejected
    cfg = BASE_CONFIG + SEACR_ON + "\n    threshold: 1\n"
    tests += 1
    if validate("threshold 1", cfg, CUTTAG_PE,
                expected_error="must be in (0, 1)"):
        passed += 1

    # 10. threshold: bad rejected
    cfg = BASE_CONFIG + SEACR_ON + "\n    threshold: bad\n"
    tests += 1
    if validate("threshold bad", cfg, CUTTAG_PE,
                expected_error="must be a float"):
        passed += 1

    # 11. --list-rules shows rules
    r = _snakemake(["--list-rules"], BASE_CONFIG + SEACR_ON, CUTTAG_PE)
    tests += 1
    if r.returncode == 0 and "seacr_bedgraph" in r.stdout + r.stderr:
        print("PASS: seacr rules in --list-rules")
        passed += 1
    else:
        print("FAIL: seacr rules not in --list-rules")

    # 12. Dry-run -p shows SEACR command details
    r = _snakemake(["-n", "-p"], BASE_CONFIG + SEACR_ON, CUTTAG_PE)
    out = r.stdout + r.stderr
    tests += 1
    shell_checks = ["SEACR_1.3.sh", "bedtools genomecov",
                    "04_peaks_seacr", "non stringent", "0.01"]
    ok = all(t in out for t in shell_checks)
    if r.returncode == 0 and ok:
        print("PASS: SEACR shell details visible in -p output")
        passed += 1
    else:
        missing = [t for t in shell_checks if t not in out]
        print("FAIL: SEACR shell check missing: %s" % missing)

    print("\nSummary: %d/%d tests passed." % (passed, tests))
    cleanup()
    if passed < tests:
        sys.exit(1)


if __name__ == "__main__":
    main()
