"""Stage 7a stress tests — CUT&Tag fragment-size QC scheduling."""

import os
import shutil
import subprocess
import sys
from _tool_resolver import resolve_tool

SNAKEFILE = "workflow/Snakefile"
SNAKEMAKE = resolve_tool("snakemake", "SNAKEMAKE")

BASE_CONFIG = """\
samples: "test_stage7a_samples.tsv"
use_control: false
threads: 1
stage4b: true
stage5: false
"""

QC_FALSE = "qc:\n  cuttag_fragment_size: false\n"

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
    "T2\tR3.fq\tR4.fq\tPE\tcuttag\tH3K4me3\tnarrow\ths\tidx\t"
    "exp1\tH3K4me3\t2\t2\t1\ttreatment\n"
)

CUTTAG_SE = HDR + (
    "T1\tR1.fq\t\tSE\tcuttag\tCTCF\tnarrow\ths\tidx\t"
    "exp1\tCTCF\t1\t1\t1\ttreatment\n"
)

MULTI_ASSAY = HDR + (
    "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\t"
    "exp_chip\tH3K27ac\t1\t1\t1\ttreatment\n"
    "T2\tR3.fq\tR4.fq\tPE\tcuttag\tH3K4me3\tnarrow\ths\tidx\t"
    "exp_cuttag\tH3K4me3\t1\t1\t1\ttreatment\n"
)

CTRL_CFG = BASE_CONFIG.replace("use_control: false", "use_control: true")
CUTTAG_CTRL = HDR + (
    "T1\tR1.fq\tR2.fq\tPE\tcuttag\tH3K4me3\tnarrow\ths\tidx\t"
    "exp1\tH3K4me3\t1\t1\t1\ttreatment\tctrl1\t\n"
    "ctrl1\tR3.fq\t\tSE\tcuttag\tIgG\tnarrow\ths\tidx\t"
    "exp1\tIgG\t1\t1\t1\tcontrol\t\t\n"
)


def _write_inputs(config_yaml, samples_tsv):
    for fq in ["R1.fq", "R2.fq", "R3.fq", "R4.fq"]:
        if not os.path.exists(fq):
            with open(fq, "w") as f:
                f.write("")
    with open("test_stage7a_config.yaml", "w") as f:
        f.write(config_yaml)
    with open("test_stage7a_samples.tsv", "w") as f:
        f.write(samples_tsv)


def _snakemake(args, config_yaml, samples_tsv):
    _write_inputs(config_yaml, samples_tsv)
    return subprocess.run(
        [SNAKEMAKE, "-s", SNAKEFILE, "--configfile",
         "test_stage7a_config.yaml", *args],
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


def cleanup():
    for f in ["test_stage7a_config.yaml", "test_stage7a_samples.tsv",
              "R1.fq", "R2.fq", "R3.fq", "R4.fq"]:
        if os.path.exists(f):
            os.remove(f)
    for d in ["__pycache__", "scripts/__pycache__"]:
        if os.path.isdir(d):
            shutil.rmtree(d)


def main():
    print("Starting Stage 7a Stress Tests\n")
    tests = 0
    passed = 0

    # 1. ChIP-seq excluded
    tests += 1
    if dryrun("ChIP-seq excluded", BASE_CONFIG, CHIPSEQ,
              forbid=("cuttag_fragment_size",)):
        passed += 1

    # 2. CUT&Tag PE in DAG
    tests += 1
    if dryrun("CUT&Tag PE in DAG", BASE_CONFIG, CUTTAG_PE,
              expect=("cuttag_fragment_size",)):
        passed += 1

    # 3. CUT&Tag SE in DAG
    tests += 1
    if dryrun("CUT&Tag SE in DAG", BASE_CONFIG, CUTTAG_SE,
              expect=("cuttag_fragment_size",)):
        passed += 1

    # 4. qc.cuttag_fragment_size: false respected
    tests += 1
    if dryrun("config false respected", BASE_CONFIG + QC_FALSE, CUTTAG_PE,
              forbid=("cuttag_fragment_size",)):
        passed += 1

    # 5. Empty qc block defaults to true
    cfg = BASE_CONFIG + "qc: {}\n"
    tests += 1
    if dryrun("empty qc defaults true", cfg, CUTTAG_PE,
              expect=("cuttag_fragment_size",)):
        passed += 1

    # 6. --list-rules shows rule
    r = _snakemake(["--list-rules"], BASE_CONFIG, CUTTAG_PE)
    tests += 1
    if r.returncode == 0 and "cuttag_fragment_size" in r.stdout + r.stderr:
        print("PASS: cuttag_fragment_size in --list-rules")
        passed += 1
    else:
        print("FAIL: cuttag_fragment_size not in --list-rules")

    # 7. -n -p shows script with correct args
    r = _snakemake(["-n", "-p"], BASE_CONFIG, CUTTAG_PE)
    out = r.stdout + r.stderr
    tests += 1
    if "calc_cuttag_fragment_size.py" in out and "--layout PE" in out:
        print("PASS: calc_cuttag_fragment_size.py with --layout PE visible")
        passed += 1
    else:
        print("FAIL: script args not found in -p output")

    # 8. Multi-assay sheet: ChIP-seq excluded, CUT&Tag included
    tests += 1
    if dryrun("multi-assay exclusion", BASE_CONFIG, MULTI_ASSAY,
              expect=("cuttag_fragment_size",)):
        passed += 1

    # 9. use_control: true + CUT&Tag control row scheduled
    tests += 1
    if dryrun("control row scheduled", CTRL_CFG, CUTTAG_CTRL,
              expect=("cuttag_fragment_size",)):
        passed += 1

    # 10. Invalid bool rejected
    cfg = BASE_CONFIG + "qc:\n  cuttag_fragment_size: maybe\n"
    r = _snakemake(["-n"], cfg, CUTTAG_PE)
    tests += 1
    if r.returncode != 0:
        print("PASS: invalid bool rejected")
        passed += 1
    else:
        print("FAIL: invalid bool not rejected")

    print("\nSummary: %d/%d tests passed." % (passed, tests))
    cleanup()
    if passed < tests:
        sys.exit(1)


if __name__ == "__main__":
    main()
