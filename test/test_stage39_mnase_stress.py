#!/usr/bin/env python3
"""Stage 39 MNase-seq stress tests: derived lists, target builder, DAG.

Validates:
- MNASE_SAMPLE_IDS / PEAK_SAMPLE_IDS derivation
- MNase targets scheduled in dry-run
- Peak-centric targets NOT scheduled for MNase samples
- Pooled MNase outputs for multi-biorep experiments
- No IDR for MNase experiments

Usage:
    SNAKEMAKE=/path/to/snakemake python3 test/test_stage39_mnase_stress.py
"""

import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
from _tool_resolver import resolve_tool

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAKEFILE = os.path.join(REPO_ROOT, "workflow", "Snakefile")
VALIDATOR = os.path.join(REPO_ROOT, "scripts", "validate_samples.py")
PROFILES_DIR = os.path.join(REPO_ROOT, "test", "profiles")
SNAKEMAKE = resolve_tool("snakemake", "SNAKEMAKE")


def _subprocess_env():
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env.setdefault("XDG_CACHE_HOME", os.path.join(tempfile.gettempdir(), "snakemake-cache"))
    return env


def _run(cmd, cwd, capture=True):
    kwargs = {"cwd": str(cwd), "env": _subprocess_env()}
    if capture:
        kwargs["capture_output"] = True
        kwargs["text"] = True
    return subprocess.run(cmd, **kwargs)


def _dry_run_mnase():
    """Set up the MNase profile and run a non-quiet dry-run. Returns (output, workdir)."""
    name = "mnase_pe_noctrl"
    profile_dir = os.path.join(PROFILES_DIR, name)
    workdir = tempfile.mkdtemp(prefix="stress39_", dir="/tmp")

    dest_samples = os.path.join(workdir, "samples.tsv")
    shutil.copy2(os.path.join(profile_dir, "samples.tsv"), dest_samples)

    config_src = os.path.join(profile_dir, "config.yaml")
    with open(config_src) as fh:
        content = fh.read()
    content = re.sub(
        r"^samples:.*$",
        f'samples: "{dest_samples}"',
        content,
        flags=re.MULTILINE,
    )
    dest_config = os.path.join(workdir, "config.yaml")
    with open(dest_config, "w") as fh:
        fh.write(content)

    for fn in ["R1.fq", "R2.fq"]:
        with open(os.path.join(workdir, fn), "w"):
            pass

    # Validate
    val = _run(
        [sys.executable, VALIDATOR, "--config", dest_config],
        cwd=workdir,
    )
    if val.returncode != 0:
        return None, workdir, "validation failed: " + val.stderr.strip()[-300:]

    # Non-quiet dry-run
    dry = _run(
        [SNAKEMAKE, "-s", SNAKEFILE, "--configfile", dest_config, "-n"],
        cwd=workdir,
    )
    if dry.returncode != 0:
        return None, workdir, "dry-run failed: " + (dry.stdout + dry.stderr).strip()[-300:]

    return (dry.stdout + dry.stderr), workdir, None


_PASS = 0
_FAIL = 0


def check(name, condition, detail=""):
    global _PASS, _FAIL
    if condition:
        print("PASS: %s" % name)
        _PASS += 1
    else:
        print("FAIL: %s%s" % (name, " — " + detail if detail else ""))
        _FAIL += 1


def main():
    global _PASS, _FAIL
    print("Starting Stage 39 MNase Stress Tests\n")

    output, workdir, err = _dry_run_mnase()
    if err:
        print("FATAL: cannot set up MNase profile: %s" % err)
        shutil.rmtree(workdir, ignore_errors=True)
        sys.exit(1)

    try:
        # --- MNase rules scheduled ---
        for rule in ["mnase_split_mono", "mnase_dyad_bigwig", "mnase_mono_bigwig"]:
            check(
                "MNase rule %s scheduled" % rule,
                rule in output,
            )

        # --- Pooled MNase rules scheduled (2 bioreps) ---
        for rule in ["mnase_pooled_mono", "mnase_pooled_dyad_bigwig", "mnase_pooled_mono_bigwig"]:
            check(
                "Pooled MNase rule %s scheduled" % rule,
                rule in output,
            )

        # --- Peak-centric rules NOT scheduled ---
        for rule in ["macs3_callpeak", "peak_counts", "frip",
                      "signal_track_fe", "signal_track_ppois",
                      "pooled_signal_track_fe", "pooled_signal_track_ppois",
                      "pooled_experiment_qc_summary"]:
            check(
                "Peak rule %s NOT scheduled" % rule,
                rule not in output,
            )

        # --- IDR NOT scheduled ---
        for rule in ["idr_true_replicates", "split_pseudoreps",
                      "idr_self_pseudoreps", "idr_pooled_pseudoreps"]:
            check(
                "IDR rule %s NOT scheduled" % rule,
                rule not in output,
            )

        # --- Preprocessing rules ARE scheduled ---
        for rule in ["fastqc", "trim_galore", "bowtie2_align",
                      "samtools_filter", "duplicate_handling",
                      "samtools_flagstat", "samtools_idxstats", "bamcoverage"]:
            check(
                "Preprocessing rule %s scheduled" % rule,
                rule in output,
            )

        # --- Replicate rules scheduled (stage4b) ---
        for rule in ["merge_biorep_bam", "pool_treatment_bam"]:
            check(
                "Replicate rule %s scheduled" % rule,
                rule in output,
            )

        # --- Direct pipeline.done target is assay-complete for MNase ---
        direct = _run(
            [
                SNAKEMAKE,
                "-s",
                SNAKEFILE,
                "--configfile",
                os.path.join(workdir, "config.yaml"),
                "-n",
                "results/M1/logs/M1.pipeline.done",
            ],
            cwd=workdir,
        )
        direct_output = direct.stdout + direct.stderr
        check(
            "Direct MNase pipeline.done dry-run succeeds",
            direct.returncode == 0,
            direct_output.strip()[-300:],
        )
        for rule in ["mnase_split_mono", "mnase_dyad_bigwig", "mnase_mono_bigwig"]:
            check(
                "Direct pipeline.done schedules %s" % rule,
                rule in direct_output,
            )

    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    # Cleanup
    pycache = os.path.join(REPO_ROOT, "scripts", "__pycache__")
    if os.path.isdir(pycache):
        shutil.rmtree(pycache)

    total = _PASS + _FAIL
    print("\nSummary: %d/%d Stage 39 stress tests passed." % (_PASS, total))
    if _FAIL:
        sys.exit(1)


if __name__ == "__main__":
    main()
