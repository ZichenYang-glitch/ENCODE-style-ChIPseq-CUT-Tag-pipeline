#!/usr/bin/env python3
"""Stage 27c stress tests — CI workflow wiring verification."""

import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CI_WORKFLOW = os.path.join(REPO_ROOT, ".github", "workflows", "ci.yml")

PASSED = 0
TOTAL = 0


def _record(name, passed):
    global PASSED, TOTAL
    TOTAL += 1
    if passed:
        PASSED += 1
        print("PASS: %s" % name)
    else:
        print("FAIL: %s" % name)


# ---------------------------------------------------------------------------
# Test 1: CI workflow file exists
# ---------------------------------------------------------------------------

def test_ci_workflow_exists():
    passed = os.path.isfile(CI_WORKFLOW)
    _record("1-ci_workflow_exists", passed)


# ---------------------------------------------------------------------------
# Test 2: Fast-checks job includes all required test commands
# ---------------------------------------------------------------------------

def test_fast_checks_commands():
    with open(CI_WORKFLOW) as fh:
        content = fh.read()

    required = [
        "test/config/test_validation.py",
        "test_stage8_smoke_profiles.py",
        "test_no_hardcoded_paths.py",
        "test_stage22_bigwig_stress.py",
        "test_stage24_qc_summary_unit.py",
        "test/manifest/test_make_manifest.py",
        "test_stage27_public_validation_plan.py",
        "test_stage27b_metadata_ci_plan.py",
        "test_stage27c_ci_workflow.py",
    ]
    missing = [s for s in required if s not in content]
    passed = len(missing) == 0
    if not passed:
        print("  Missing tests:", missing)
    _record("2-fast_checks_commands", passed)


# ---------------------------------------------------------------------------
# Test 3: workflow_dispatch exists
# ---------------------------------------------------------------------------

def test_workflow_dispatch_exists():
    with open(CI_WORKFLOW) as fh:
        content = fh.read()

    passed = "workflow_dispatch" in content
    _record("3-workflow_dispatch_exists", passed)


# ---------------------------------------------------------------------------
# Test 4: No public data download commands in CI
# ---------------------------------------------------------------------------

def test_no_public_data_downloads():
    with open(CI_WORKFLOW) as fh:
        content = fh.read()

    forbidden = [
        "wget", "curl", "fastq-dump", "prefetch", "fasterq-dump",
        "gsutil", "aws s3", "ENA_", "SRR",
        "ENCSR000DYI", "ENCSR000AKB", "ENCSR254KDA", "GSE145187",
    ]
    found = [s for s in forbidden if s.lower() in content.lower()]
    passed = len(found) == 0
    if not passed:
        print("  Found forbidden patterns:", found)
    _record("4-no_public_data_downloads", passed)


# ---------------------------------------------------------------------------
# Test 5: Tiny real execution still in workflow_dispatch
# ---------------------------------------------------------------------------

def test_tiny_exec_in_workflow_dispatch():
    with open(CI_WORKFLOW) as fh:
        content = fh.read()

    passed = ("test_stage8b_tiny_execution.py" in content
              and "workflow_dispatch" in content)
    _record("5-tiny_exec_in_workflow_dispatch", passed)


# ---------------------------------------------------------------------------
# Test 6: Validate default config command present
# ---------------------------------------------------------------------------

def test_validate_config_command():
    with open(CI_WORKFLOW) as fh:
        content = fh.read()

    passed = "validate_samples.py --config config/config.yaml" in content
    _record("6-validate_config_command", passed)


# ---------------------------------------------------------------------------
# Test 7: No Snakemake rule changes in this stage
# ---------------------------------------------------------------------------
# (Not checking Snakemake rules here — just verifying this test is about CI.)

def test_ci_fast_uses_ci_fast_env():
    with open(CI_WORKFLOW) as fh:
        content = fh.read()

    passed = "ci-fast.yml" in content
    _record("7-ci_fast_uses_ci_fast_env", passed)


# ---------------------------------------------------------------------------
def main():
    print("Starting Stage 27c CI Workflow Tests\n")

    test_ci_workflow_exists()
    test_fast_checks_commands()
    test_workflow_dispatch_exists()
    test_no_public_data_downloads()
    test_tiny_exec_in_workflow_dispatch()
    test_validate_config_command()
    test_ci_fast_uses_ci_fast_env()

    print("\nSummary: %d/%d tests passed." % (PASSED, TOTAL))

    if PASSED < TOTAL:
        sys.exit(1)


if __name__ == "__main__":
    main()
