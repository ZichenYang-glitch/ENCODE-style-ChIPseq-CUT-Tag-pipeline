#!/usr/bin/env python3
"""Stage 8a smoke-test harness for ChIP-seq/CUT&Tag Snakemake profiles.

Run all committed test profiles through validate_samples.py and a
Snakemake dry-run.  Every temporary file is created under /tmp and
cleaned up afterward — nothing is written into the repository.

Usage:
    SNAKEMAKE=/path/to/snakemake python3 test/test_stage8_smoke_profiles.py
"""

import subprocess
import tempfile
import shutil
import os
import sys
import csv
import re


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAKEFILE = os.path.join(REPO_ROOT, "workflow", "Snakefile")
VALIDATOR = os.path.join(REPO_ROOT, "scripts", "validate_samples.py")
PROFILES_DIR = os.path.join(REPO_ROOT, "test", "profiles")

# SNAKEMAKE resolution: env var > known Conda path > bare "snakemake"
_SNAKEMAKE_ENV = os.environ.get("SNAKEMAKE", "")
_SNAKEMAKE_CONDA = "/home/irenadler/miniconda3/envs/chipseq/bin/snakemake"
if _SNAKEMAKE_ENV:
    SNAKEMAKE = _SNAKEMAKE_ENV
elif os.path.isfile(_SNAKEMAKE_CONDA):
    SNAKEMAKE = _SNAKEMAKE_CONDA
else:
    SNAKEMAKE = "snakemake"

# All seven profiles run by default.
PROFILES = [
    "chipseq_se_noctrl",
    "chipseq_pe_noctrl",
    "chipseq_pe_ctrlsample",
    "cuttag_pe_noctrl",
    "cuttag_pe_seacr",
    "chipseq_idr_dryrun",
    "chipseq_pe_external_ctrlbam",
]

# Profiles that get an additional non-quiet dry-run to verify specific
# rules are actually scheduled (not just defined in the Snakefile).
SCHEDULING_CHECKS = {
    "cuttag_pe_noctrl":      ["cuttag_fragment_size"],
    "cuttag_pe_seacr":       ["seacr_bedgraph", "seacr_call"],
    "chipseq_idr_dryrun":    [
        "macs3_idr_biorep", "idr_true_replicates",
        "split_pseudoreps", "idr_self_pseudoreps",
        "idr_pooled_pseudoreps", "stage5b_summary",
    ],
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _subprocess_env():
    """Return os.environ with PYTHONDONTWRITEBYTECODE=1."""
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    return env


def _discover_placeholders(samples_tsv_path):
    """Parse a samples.tsv and return paths that need empty placeholder files.

    Collects fastq_1, fastq_2 (when non-empty), and control_bam (when
    non-empty).  bowtie2_index is deliberately excluded — a fake single
    file is misleading for a multi-file index prefix.
    """
    paths = set()
    with open(samples_tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fq1 = (row.get("fastq_1") or "").strip()
            fq2 = (row.get("fastq_2") or "").strip()
            cb = (row.get("control_bam") or "").strip()
            if fq1:
                paths.add(fq1)
            if fq2:
                paths.add(fq2)
            if cb:
                paths.add(cb)
    return paths


def _rewrite_config(profile_config_path, workdir):
    """Read config.yaml and replace the ``samples:`` line with an
    absolute path pointing into *workdir*."""
    with open(profile_config_path) as fh:
        content = fh.read()
    abs_samples = os.path.join(workdir, "samples.tsv")
    content = re.sub(
        r"^samples:.*$",
        f'samples: "{abs_samples}"',
        content,
        flags=re.MULTILINE,
    )
    return content


def _run(cmd, cwd, capture=True):
    """Run a command with *cwd*, PYTHONDONTWRITEBYTECODE, and optional
    capture.  Returns a CompletedProcess."""
    kwargs = {"cwd": str(cwd), "env": _subprocess_env()}
    if capture:
        kwargs["capture_output"] = True
        kwargs["text"] = True
    return subprocess.run(cmd, **kwargs)


# ---------------------------------------------------------------------------
# One profile
# ---------------------------------------------------------------------------

def _test_profile(name):
    """Validate and dry-run a single profile.  Returns True on pass."""

    profile_dir = os.path.join(PROFILES_DIR, name)
    samples_tsv_src = os.path.join(profile_dir, "samples.tsv")
    config_yaml_src = os.path.join(profile_dir, "config.yaml")

    if not os.path.isdir(profile_dir):
        print("FAIL: %s — profile directory not found: %s" % (name, profile_dir))
        return False

    workdir = tempfile.mkdtemp(prefix="smoke_%s_" % name, dir="/tmp")

    try:
        # --- discover and create placeholder files in workdir ---
        placeholders = _discover_placeholders(samples_tsv_src)
        for rel_path in placeholders:
            placeholder_path = os.path.join(workdir, os.path.basename(rel_path))
            with open(placeholder_path, "w"):
                pass  # empty file

        # --- copy samples.tsv (paths stay relative, resolve against workdir) ---
        dest_samples = os.path.join(workdir, "samples.tsv")
        shutil.copy2(samples_tsv_src, dest_samples)

        # --- write config.yaml with absolute samples path ---
        rewritten_config = _rewrite_config(config_yaml_src, workdir)
        dest_config = os.path.join(workdir, "config.yaml")
        with open(dest_config, "w") as fh:
            fh.write(rewritten_config)

        # 1. validate (cwd=workdir — relative control_bam resolves here)
        try:
            val = _run(
                [sys.executable, VALIDATOR, "--config", dest_config],
                cwd=workdir,
            )
        except (OSError, FileNotFoundError, PermissionError) as exc:
            print("FAIL: %s — validation: %s" % (name, exc))
            return False
        if val.returncode != 0:
            tail = val.stderr.strip()[-500:]
            print("FAIL: %s — validation failed\n   %s" % (name, tail))
            return False

        # 2. quiet dry-run
        try:
            dry = _run(
                [SNAKEMAKE, "-s", SNAKEFILE, "--configfile", dest_config,
                 "-n", "--quiet"],
                cwd=workdir,
            )
        except (OSError, FileNotFoundError, PermissionError) as exc:
            print("FAIL: %s — quiet dry-run: %s" % (name, exc))
            return False
        if dry.returncode != 0:
            tail = (dry.stdout + dry.stderr).strip()[-500:]
            print("FAIL: %s — quiet dry-run failed\n   %s" % (name, tail))
            return False

        # 3. non-quiet dry-run for scheduling checks (selected profiles only)
        expected_rules = SCHEDULING_CHECKS.get(name, [])
        if expected_rules:
            try:
                sched = _run(
                    [SNAKEMAKE, "-s", SNAKEFILE, "--configfile", dest_config, "-n"],
                    cwd=workdir,
                )
            except (OSError, FileNotFoundError, PermissionError) as exc:
                print("FAIL: %s — scheduling dry-run: %s" % (name, exc))
                return False
            if sched.returncode != 0:
                tail = (sched.stdout + sched.stderr).strip()[-500:]
                print("FAIL: %s — scheduling dry-run failed\n   %s" % (name, tail))
                return False

            output = sched.stdout + sched.stderr
            for rule_name in expected_rules:
                if rule_name not in output:
                    print("FAIL: %s — rule %r not scheduled" % (name, rule_name))
                    return False

        print("PASS: %s" % name)
        return True

    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    print("Starting Stage 8a Smoke Profile Tests\n")
    passed = 0
    total = len(PROFILES)

    for name in PROFILES:
        if _test_profile(name):
            passed += 1

    print("\nSummary: %d/%d profiles passed." % (passed, total))

    # Remove scripts/__pycache__ if the import of validate_samples created it
    pycache = os.path.join(REPO_ROOT, "scripts", "__pycache__")
    if os.path.isdir(pycache):
        shutil.rmtree(pycache)

    if passed < total:
        sys.exit(1)


if __name__ == "__main__":
    main()
