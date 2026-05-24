#!/usr/bin/env python3
"""Stage 32 stress tests — public data execution report scaffold."""

import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RELEASE_CHECKS = os.path.join(REPO_ROOT, "docs", "release-checks")
TEMPLATE = os.path.join(RELEASE_CHECKS, "public-data-execution-report-template.md")
STUBS_DIR = os.path.join(RELEASE_CHECKS, "public-data-runs")
HELPER = os.path.join(REPO_ROOT, "scripts", "prepare_public_validation_inputs.py")
PLAN27 = os.path.join(RELEASE_CHECKS, "stage27-public-data-validation-plan.md")

STUB_QUEUES = [
    "tf_chip_cebpb",
    "broad_histone_h3k27me3",
    "atac_keratinocyte",
    "cuttag_h3k27me3",
]

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
# Test 1: Template exists
# ---------------------------------------------------------------------------

def test_template_exists():
    passed = os.path.isfile(TEMPLATE)
    _record("1-template_exists", passed)


# ---------------------------------------------------------------------------
# Test 2: Four stubs exist
# ---------------------------------------------------------------------------

def test_four_stubs_exist():
    missing = []
    for q in STUB_QUEUES:
        path = os.path.join(STUBS_DIR, f"{q}.md")
        if not os.path.isfile(path):
            missing.append(q)
    passed = len(missing) == 0
    if not passed:
        print("  Missing:", missing)
    _record("2-four_stubs_exist", passed)


# ---------------------------------------------------------------------------
# Test 3: Each stub says "not executed yet"
# ---------------------------------------------------------------------------

def test_stubs_not_executed():
    missing = []
    for q in STUB_QUEUES:
        path = os.path.join(STUBS_DIR, f"{q}.md")
        with open(path) as fh:
            content = fh.read()
        if "not executed yet" not in content.lower():
            missing.append(q)
    passed = len(missing) == 0
    if not passed:
        print("  Missing 'not executed yet':", missing)
    _record("3-stubs_not_executed", passed)


# ---------------------------------------------------------------------------
# Test 4: Each stub links to template
# ---------------------------------------------------------------------------

def test_stubs_link_to_template():
    missing = []
    for q in STUB_QUEUES:
        path = os.path.join(STUBS_DIR, f"{q}.md")
        with open(path) as fh:
            content = fh.read()
        if "public-data-execution-report-template.md" not in content:
            missing.append(q)
    passed = len(missing) == 0
    if not passed:
        print("  Missing template link:", missing)
    _record("4-stubs_link_to_template", passed)


# ---------------------------------------------------------------------------
# Test 5: Stage 27 plan links to template and all stubs
# ---------------------------------------------------------------------------

def test_plan27_links():
    with open(PLAN27) as fh:
        content = fh.read()
    missing = []
    if "public-data-execution-report-template.md" not in content:
        missing.append("template")
    for q in STUB_QUEUES:
        if f"public-data-runs/{q}.md" not in content:
            missing.append(q)
    passed = len(missing) == 0
    if not passed:
        print("  Missing links:", missing)
    _record("5-plan27_links", passed)


# ---------------------------------------------------------------------------
# Test 6: Helper --report-stubs works
# ---------------------------------------------------------------------------

def test_helper_report_stubs():
    result = subprocess.run(
        [sys.executable, HELPER, "--report-stubs"],
        capture_output=True, text=True,
    )
    output = result.stdout
    passed = (result.returncode == 0
              and "Report Stubs" in output
              and all(q in output for q in STUB_QUEUES)
              and "template" in output.lower())
    if not passed:
        print("  rc=%s has_report=%s queues_present=%s" % (
            result.returncode,
            "Report Stubs" in output,
            all(q in output for q in STUB_QUEUES),
        ))
    _record("6-helper_report_stubs", passed)


# ---------------------------------------------------------------------------
# Test 7: No data files committed
# ---------------------------------------------------------------------------

def test_no_data_files():
    forbidden = {
        ".fq", ".fq.gz", ".fastq", ".fastq.gz",
        ".bam", ".bai",
        ".bw", ".bigWig",
        ".bdg", ".bedGraph",
        ".html",
        ".bt2", ".bt2l", ".rev.1.bt2", ".rev.2.bt2",
        ".rev.1.bt2l", ".rev.2.bt2l",
    }
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
    _record("7-no_data_files", passed)


# ---------------------------------------------------------------------------
def main():
    print("Starting Stage 32 Public Report Scaffold Tests\n")

    test_template_exists()
    test_four_stubs_exist()
    test_stubs_not_executed()
    test_stubs_link_to_template()
    test_plan27_links()
    test_helper_report_stubs()
    test_no_data_files()

    print("\nSummary: %d/%d tests passed." % (PASSED, TOTAL))

    pycache = os.path.join(REPO_ROOT, "scripts", "__pycache__")
    if os.path.isdir(pycache):
        import shutil
        shutil.rmtree(pycache)

    if PASSED < TOTAL:
        sys.exit(1)


if __name__ == "__main__":
    main()
