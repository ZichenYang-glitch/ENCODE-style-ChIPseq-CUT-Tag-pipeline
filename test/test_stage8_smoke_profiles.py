#!/usr/bin/env python3
"""Pytest-based Stage 8 smoke tests for all standard Snakemake profiles.

This module replaces the standalone test_stage8_smoke_profiles.py script.
It reuses conftest.py fixtures and the shared prepare_profile_workdir helper
so that profile setup logic is not duplicated with test_dag_snapshots.py.
"""

import os
import shutil
import subprocess
import sys

import pytest

from _tool_resolver import resolve_tool
from conftest import SMOKE_PROFILES


SNAKEMAKE = resolve_tool("snakemake", "SNAKEMAKE")

# Profiles that get an additional non-quiet dry-run to verify specific
# rules are actually scheduled (not just defined in the Snakefile).
SCHEDULING_CHECKS = {
    "cuttag_pe_noctrl": ["cuttag_fragment_size"],
    "cuttag_pe_seacr": ["seacr_bedgraph", "seacr_call"],
    "chipseq_idr_dryrun": [
        "macs3_idr_biorep",
        "idr_true_replicates",
        "split_pseudoreps",
        "idr_self_pseudoreps",
        "idr_pooled_pseudoreps",
        "stage5b_summary",
    ],
}


def _subprocess_env():
    """Return os.environ with PYTHONDONTWRITEBYTECODE=1."""
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    return env


def _run(cmd, cwd, capture=True):
    """Run a command with cwd and PYTHONDONTWRITEBYTECODE."""
    kwargs = {"cwd": str(cwd), "env": _subprocess_env()}
    if capture:
        kwargs["capture_output"] = True
        kwargs["text"] = True
    return subprocess.run(cmd, **kwargs)


@pytest.mark.parametrize("profile", SMOKE_PROFILES)
def test_validate_samples(profile, profiles_dir, validator_script):
    """Run scripts/validate_samples.py against each profile config."""
    from conftest import prepare_profile_workdir

    profile_dir = os.path.join(profiles_dir, profile)
    workdir, dest_config = prepare_profile_workdir(profile_dir)
    try:
        result = _run(
            [sys.executable, validator_script, "--config", dest_config],
            cwd=workdir,
        )
        assert result.returncode == 0, (
            f"validate_samples.py failed for {profile}:\n{result.stderr}"
        )
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


@pytest.mark.parametrize("profile", SMOKE_PROFILES)
def test_quiet_dryrun(profile, profiles_dir, snakefile):
    """Run snakemake -n --quiet for each profile."""
    from conftest import prepare_profile_workdir

    profile_dir = os.path.join(profiles_dir, profile)
    workdir, dest_config = prepare_profile_workdir(profile_dir)
    try:
        result = _run(
            [
                SNAKEMAKE,
                "-s",
                snakefile,
                "--configfile",
                dest_config,
                "-n",
                "--quiet",
            ],
            cwd=workdir,
        )
        assert result.returncode == 0, (
            f"Quiet dry-run failed for {profile}:\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


@pytest.mark.parametrize("profile", SMOKE_PROFILES)
def test_scheduling_rules(profile, profiles_dir, snakefile):
    """Verify expected rules are scheduled for profiles with scheduling checks."""
    expected_rules = SCHEDULING_CHECKS.get(profile)
    if not expected_rules:
        pytest.skip(f"No scheduling checks configured for {profile}")

    from conftest import prepare_profile_workdir

    profile_dir = os.path.join(profiles_dir, profile)
    workdir, dest_config = prepare_profile_workdir(profile_dir)
    try:
        result = _run(
            [SNAKEMAKE, "-s", snakefile, "--configfile", dest_config, "-n"],
            cwd=workdir,
        )
        assert result.returncode == 0, (
            f"Scheduling dry-run failed for {profile}:\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
        output = result.stdout + result.stderr
        for rule_name in expected_rules:
            assert rule_name in output, (
                f"{profile}: expected rule {rule_name!r} not scheduled"
            )
    finally:
        shutil.rmtree(workdir, ignore_errors=True)
