#!/usr/bin/env python3
"""Stage 27b stress tests — metadata verification helper and CI/CD plan."""

import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PREPARE_SCRIPT = os.path.join(REPO_ROOT, "scripts",
                              "prepare_public_validation_inputs.py")
CI_PLAN_DOC = os.path.join(REPO_ROOT, "docs", "release-checks",
                           "stage27b-metadata-ci-plan.md")
PLAN_DOC = os.path.join(REPO_ROOT, "docs", "release-checks",
                        "stage27-public-data-validation-plan.md")

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
# Test 1: CI plan doc exists
# ---------------------------------------------------------------------------

def test_ci_plan_doc_exists():
    passed = os.path.isfile(CI_PLAN_DOC)
    _record("1-ci_plan_doc_exists", passed)


# ---------------------------------------------------------------------------
# Test 2: Helper --check-metadata works without network
# ---------------------------------------------------------------------------

def test_check_metadata():
    result = subprocess.run(
        [sys.executable, PREPARE_SCRIPT, "--check-metadata"],
        capture_output=True, text=True,
    )
    output = result.stdout
    passed = (result.returncode == 0
              and "Metadata Checklist" in output
              and "ENCSR000DYI" in output
              and "ENCSR000AKB" in output
              and "GSE145187" in output
              and "No network requests" in output
              and "Acceptance Criteria" in output)
    if not passed:
        print("  rc=%s output_has_checklist=%s has_encsr=%s has_geo=%s no_net=%s"
              % (result.returncode,
                 "Metadata Checklist" in output,
                 "ENCSR000DYI" in output,
                 "GSE145187" in output,
                 "No network requests" in output))
    _record("2-check_metadata_no_network", passed)


# ---------------------------------------------------------------------------
# Test 3: CI plan mentions PR, workflow_dispatch, external/manual
# ---------------------------------------------------------------------------

def test_ci_plan_tiers():
    with open(CI_PLAN_DOC) as fh:
        content = fh.read()

    required = [
        "Tier 1", "Tier 2", "Tier 3",
        "workflow_dispatch",
        "External/Manual Public Data Run",
        "Fast PR checks",
        "pull_request",
    ]
    missing = [s for s in required if s not in content]
    passed = len(missing) == 0
    if not passed:
        print("  Missing:", missing)
    _record("3-ci_plan_tiers", passed)


# ---------------------------------------------------------------------------
# Test 4: No data files committed (same check as Stage 27a)
# ---------------------------------------------------------------------------

def test_no_data_files():
    forbidden = {".fq.gz", ".fastq.gz", ".bam", ".bai", ".bw", ".bdg"}
    found = []
    for root, dirs, files in os.walk(REPO_ROOT):
        dirs[:] = [d for d in dirs if not d.startswith(".")]
        for f in files:
            for ext in forbidden:
                if f.endswith(ext):
                    found.append(os.path.relpath(
                        os.path.join(root, f), REPO_ROOT))
    passed = len(found) == 0
    if not passed:
        print("  Found:", found)
    _record("4-no_data_files", passed)


# ---------------------------------------------------------------------------
# Test 5: All Stage 27a accessions still present in helper
# ---------------------------------------------------------------------------

def test_accessions_present():
    result = subprocess.run(
        [sys.executable, PREPARE_SCRIPT, "--json"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        _record("5-accessions_present", False)
        return

    import json
    data = json.loads(result.stdout)
    accessions = {d["accession"] for d in data}
    expected = {"ENCSR000DYI", "ENCSR000AKB", "ENCSR254KDA", "GSE145187"}
    passed = expected.issubset(accessions)
    if not passed:
        print("  Missing:", expected - accessions)
    _record("5-accessions_present", passed)


# ---------------------------------------------------------------------------
# Test 6: CI plan doc contains artifact policy
# ---------------------------------------------------------------------------

def test_ci_plan_artifact_policy():
    with open(CI_PLAN_DOC) as fh:
        content = fh.read()

    required = [
        "Artifact Policy",
        "NOT committed",
        "Committed to repo",
        "Why Public Data Is Not Run on Every PR",
    ]
    missing = [s for s in required if s not in content]
    passed = len(missing) == 0
    if not passed:
        print("  Missing:", missing)
    _record("6-ci_plan_artifact_policy", passed)


# ---------------------------------------------------------------------------
# Test 7: CI plan doc references both stage docs
# ---------------------------------------------------------------------------

def test_ci_plan_cross_refs():
    """Verify CI plan references the Stage 27a validation plan."""
    with open(CI_PLAN_DOC) as fh:
        content = fh.read()
    passed = "stage27-public-data-validation-plan" in content
    _record("7-ci_plan_cross_refs", passed)


# ---------------------------------------------------------------------------
def main():
    print("Starting Stage 27b Metadata CI Plan Tests\n")

    test_ci_plan_doc_exists()
    test_check_metadata()
    test_ci_plan_tiers()
    test_no_data_files()
    test_accessions_present()
    test_ci_plan_artifact_policy()
    test_ci_plan_cross_refs()

    print("\nSummary: %d/%d tests passed." % (PASSED, TOTAL))

    pycache = os.path.join(REPO_ROOT, "scripts", "__pycache__")
    if os.path.isdir(pycache):
        import shutil
        shutil.rmtree(pycache)

    if PASSED < TOTAL:
        sys.exit(1)


if __name__ == "__main__":
    main()
