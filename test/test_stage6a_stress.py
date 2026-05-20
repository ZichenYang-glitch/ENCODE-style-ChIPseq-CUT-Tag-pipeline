"""Stage 6a stress tests — pooled experiment signal tracks.

Checks that multi-biorep experiments schedule pooled FE/ppois bedGraph signal
tracks when qc.signal_tracks is enabled, while single-sample and disabled-signal
configs keep the DAG unchanged.
"""

import os
import shutil
import subprocess
import sys
from _tool_resolver import resolve_tool

SNAKEFILE = "workflow/Snakefile"
SNAKEMAKE = resolve_tool("snakemake", "SNAKEMAKE")

BASE_CONFIG = """\
samples: "test_stage6a_samples.tsv"
use_control: false
threads: 1
stage4b: true
stage5: false
"""

QC_FALSE = """\
qc:
  signal_tracks: false
"""

HDR = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\t"
    "peak_mode\tgenome\tbowtie2_index\texperiment\tcondition\t"
    "replicate\tbiological_replicate\ttechnical_replicate\trole\n"
)

SINGLE_SAMPLE = HDR + (
    "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\t"
    "exp1\tH3K27ac\t1\t1\t1\ttreatment\n"
)

TWO_BIOREPS_NARROW = HDR + (
    "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\t"
    "exp1\tH3K27ac\t1\t1\t1\ttreatment\n"
    "T2\tR2.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\t"
    "exp1\tH3K27ac\t2\t2\t1\ttreatment\n"
)

TWO_BIOREPS_BROAD = HDR + (
    "T1\tR1.fq\t\tSE\tchipseq\tH3K27me3\tbroad\ths\tidx\t"
    "exp1\tH3K27me3\t1\t1\t1\ttreatment\n"
    "T2\tR2.fq\t\tSE\tchipseq\tH3K27me3\tbroad\ths\tidx\t"
    "exp1\tH3K27me3\t2\t2\t1\ttreatment\n"
)


def _write_inputs(config_yaml, samples_tsv):
    for fq in ["R1.fq", "R2.fq"]:
        if not os.path.exists(fq):
            with open(fq, "w") as f:
                f.write("")
    with open("test_stage6a_config.yaml", "w") as f:
        f.write(config_yaml)
    with open("test_stage6a_samples.tsv", "w") as f:
        f.write(samples_tsv)


def _snakemake(args, config_yaml, samples_tsv):
    _write_inputs(config_yaml, samples_tsv)
    return subprocess.run(
        [
            SNAKEMAKE,
            "-s",
            SNAKEFILE,
            "--configfile",
            "test_stage6a_config.yaml",
            *args,
        ],
        capture_output=True,
        text=True,
    )


def dryrun(name, config_yaml, samples_tsv, expect=(), forbid=(), printshell=False):
    args = ["-n", "-p"] if printshell else ["-n"]
    result = _snakemake(args, config_yaml, samples_tsv)
    output = result.stdout + result.stderr
    if result.returncode != 0:
        print(
            "FAIL: %s\n   Dry-run failed.\n   stdout: %s\n   stderr: %s"
            % (name, result.stdout.strip()[-500:], result.stderr.strip()[-500:])
        )
        return False
    for text in expect:
        if text not in output:
            print("FAIL: %s\n   Expected text %r not found." % (name, text))
            return False
    for text in forbid:
        if text in output:
            print("FAIL: %s\n   Forbidden text %r found." % (name, text))
            return False
    print("PASS: %s" % name)
    return True


def list_rules(name, expect=()):
    result = _snakemake(["--list-rules"], BASE_CONFIG, TWO_BIOREPS_NARROW)
    output = result.stdout + result.stderr
    if result.returncode != 0:
        print("FAIL: %s\n   --list-rules failed." % name)
        return False
    for text in expect:
        if text not in output:
            print("FAIL: %s\n   Rule %r not listed." % (name, text))
            return False
    print("PASS: %s" % name)
    return True


def cleanup():
    for f in [
        "test_stage6a_config.yaml",
        "test_stage6a_samples.tsv",
        "R1.fq",
        "R2.fq",
    ]:
        if os.path.exists(f):
            os.remove(f)
    for d in ["__pycache__", "scripts/__pycache__"]:
        if os.path.isdir(d):
            shutil.rmtree(d)


def main():
    print("Starting Stage 6a Pooled Signal Stress Tests\n")
    tests = 0
    passed = 0

    tests += 1
    if dryrun(
        "single sample has no pooled signal DAG",
        BASE_CONFIG,
        SINGLE_SAMPLE,
        forbid=("pooled_signal_track_fe", "pooled_signal_track_ppois"),
    ):
        passed += 1

    tests += 1
    if dryrun(
        "multi-biorep schedules pooled signal rules",
        BASE_CONFIG,
        TWO_BIOREPS_NARROW,
        expect=("pooled_signal_track_fe", "pooled_signal_track_ppois"),
    ):
        passed += 1

    tests += 1
    if dryrun(
        "pooled signal output paths visible",
        BASE_CONFIG,
        TWO_BIOREPS_NARROW,
        expect=("exp1.pooled.FE.bdg", "exp1.pooled.ppois.bdg"),
        printshell=True,
    ):
        passed += 1

    tests += 1
    if dryrun(
        "pooled MACS3 emits bedGraphs with -B",
        BASE_CONFIG,
        TWO_BIOREPS_NARROW,
        expect=("exp1_pooled", "-B"),
        printshell=True,
    ):
        passed += 1

    tests += 1
    if dryrun(
        "qc.signal_tracks false disables pooled signal rules",
        BASE_CONFIG + QC_FALSE,
        TWO_BIOREPS_NARROW,
        forbid=("pooled_signal_track_fe", "pooled_signal_track_ppois"),
    ):
        passed += 1

    tests += 1
    if dryrun(
        "broad pooled experiment schedules pooled signal rules",
        BASE_CONFIG,
        TWO_BIOREPS_BROAD,
        expect=("pooled_signal_track_fe", "pooled_signal_track_ppois"),
    ):
        passed += 1

    tests += 1
    if dryrun(
        "broad pooled peak mode still uses broadPeak",
        BASE_CONFIG,
        TWO_BIOREPS_BROAD,
        expect=("broadPeak", "exp1.pooled.FE.bdg"),
        printshell=True,
    ):
        passed += 1

    tests += 1
    if list_rules(
        "pooled signal rules listed",
        expect=("pooled_signal_track_fe", "pooled_signal_track_ppois"),
    ):
        passed += 1

    print("\nSummary: %d/%d tests passed." % (passed, tests))
    cleanup()
    if passed < tests:
        sys.exit(1)


if __name__ == "__main__":
    main()
