#!/usr/bin/env python3
"""Stage 27a stress tests — public data validation plan and helper skeleton."""

import json
import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PREPARE_SCRIPT = os.path.join(REPO_ROOT, "scripts",
                              "prepare_public_validation_inputs.py")
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
# Test 1: Plan doc exists
# ---------------------------------------------------------------------------

def test_plan_doc_exists():
    passed = os.path.isfile(PLAN_DOC)
    _record("1-plan_doc_exists", passed)


# ---------------------------------------------------------------------------
# Test 2: No large data files committed
# ---------------------------------------------------------------------------

def test_no_data_files_committed():
    """Verify no FASTQ, BAM, BAI, BW, BDG, or HTML data files exist in the repo."""
    forbidden_exts = {".fq", ".fq.gz", ".fastq", ".fastq.gz",
                      ".bam", ".bai", ".bw", ".bdg",
                      ".html"}
    found = []
    skip_dirs = {".git", "frontend"}
    for root, dirs, files in os.walk(REPO_ROOT):
        # Skip .git, frontend source tree, and hidden dirs
        dirs[:] = [d for d in dirs if not d.startswith(".") and d not in skip_dirs]
        for f in files:
            for ext in forbidden_exts:
                if f.endswith(ext) or f.endswith(ext.replace(".fq", ".fastq")):
                    found.append(os.path.join(root, f))
    # Also check for common patterns
    for ext_check in [".fq.gz", ".bam", ".bai", ".bw", ".bdg"]:
        for root, dirs, files in os.walk(REPO_ROOT):
            dirs[:] = [d for d in dirs if not d.startswith(".") and d not in skip_dirs]
            for f in files:
                if f.endswith(ext_check):
                    found.append(os.path.join(root, f))

    found = list(set(found))
    passed = len(found) == 0
    if not passed:
        print("  Found:", found)
    _record("2-no_data_files_committed", passed)


# ---------------------------------------------------------------------------
# Test 3: Helper supports --dry-run
# ---------------------------------------------------------------------------

def test_helper_dry_run():
    result = subprocess.run(
        [sys.executable, PREPARE_SCRIPT, "--dry-run"],
        capture_output=True, text=True,
    )
    output = result.stdout + result.stderr
    passed = (result.returncode == 0
              and "[dry-run]" in output
              and "No downloads performed" in output)
    if not passed:
        print("  rc=%s output=%s" % (result.returncode, output[-200:]))
    _record("3-helper_dry_run", passed)


# ---------------------------------------------------------------------------
# Test 4: Helper emits deterministic TSV output
# ---------------------------------------------------------------------------

def test_helper_tsv_output():
    result1 = subprocess.run(
        [sys.executable, PREPARE_SCRIPT],
        capture_output=True, text=True,
    )
    result2 = subprocess.run(
        [sys.executable, PREPARE_SCRIPT],
        capture_output=True, text=True,
    )
    passed = (result1.returncode == 0
              and result2.returncode == 0
              and result1.stdout == result2.stdout
              and "tf_chip_cebpb" in result1.stdout
              and "ENCSR000DYI" in result1.stdout
              and "cuttag_h3k27me3" in result1.stdout)
    if not passed:
        print("  deterministic=%s has_tf=%s has_encsr=%s has_cuttag=%s"
              % (result1.stdout == result2.stdout,
                 "tf_chip_cebpb" in result1.stdout,
                 "ENCSR000DYI" in result1.stdout,
                 "cuttag_h3k27me3" in result1.stdout))
    _record("4-helper_tsv_output", passed)


# ---------------------------------------------------------------------------
# Test 5: Helper emits JSON output
# ---------------------------------------------------------------------------

def test_helper_json_output():
    result = subprocess.run(
        [sys.executable, PREPARE_SCRIPT, "--json"],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        _record("5-helper_json_output", False)
        return

    try:
        data = json.loads(result.stdout)
    except json.JSONDecodeError:
        _record("5-helper_json_output", False)
        return

    passed = (isinstance(data, list)
              and len(data) == 4
              and data[0]["queue"] == "tf_chip_cebpb"
              and data[0]["assay"] == "chipseq")
    _record("5-helper_json_output", passed)


# ---------------------------------------------------------------------------
# Test 6: Plan doc contains key sections
# ---------------------------------------------------------------------------

def test_plan_doc_sections():
    with open(PLAN_DOC) as fh:
        content = fh.read()

    required = [
        "Validation Goals",
        "Candidate Datasets",
        "ENCSR000DYI",
        "Genome / Build Assumptions",
        "Required Reference Resources",
        "Expected Config / Sample Sheet Shape",
        "Expected Outputs to Audit",
        "Known Runtime / Storage Constraints",
        "No Data Committed to Repo",
        "Execution Tiers",
    ]
    missing = [s for s in required if s not in content]
    passed = len(missing) == 0
    if not passed:
        print("  Missing sections:", missing)
    _record("6-plan_doc_sections", passed)


# ---------------------------------------------------------------------------
# Test 7: Accession examples present in plan
# ---------------------------------------------------------------------------

def test_plan_accessions():
    with open(PLAN_DOC) as fh:
        content = fh.read()

    accessions = ["ENCSR000DYI", "ENCSR000AKB", "ENCSR254KDA", "GSE145187"]
    missing = [a for a in accessions if a not in content]
    passed = len(missing) == 0
    if not passed:
        print("  Missing accessions:", missing)
    _record("7-plan_accessions", passed)


# ---------------------------------------------------------------------------
def main():
    print("Starting Stage 27a Public Validation Plan Tests\n")

    test_plan_doc_exists()
    test_no_data_files_committed()
    test_helper_dry_run()
    test_helper_tsv_output()
    test_helper_json_output()
    test_plan_doc_sections()
    test_plan_accessions()

    print("\nSummary: %d/%d tests passed." % (PASSED, TOTAL))

    pycache = os.path.join(REPO_ROOT, "scripts", "__pycache__")
    if os.path.isdir(pycache):
        import shutil
        shutil.rmtree(pycache)

    if PASSED < TOTAL:
        sys.exit(1)


if __name__ == "__main__":
    main()
