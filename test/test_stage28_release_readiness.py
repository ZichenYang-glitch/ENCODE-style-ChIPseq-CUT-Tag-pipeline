#!/usr/bin/env python3
"""Stage 28 release readiness tests — doc consistency and stale-phrase guard."""

import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

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


def _read(relpath):
    path = os.path.join(REPO_ROOT, relpath)
    if not os.path.isfile(path):
        return ""
    with open(path) as fh:
        return fh.read()


# ---------------------------------------------------------------------------
# Test 1: No stale "CI expansion designed but not wired"
# ---------------------------------------------------------------------------

def test_no_stale_ci_status():
    content = _read("docs/release-checks/stage27b-metadata-ci-plan.md")
    passed = "CI expansion designed but not wired" not in content
    _record("1-no_stale_ci_status", passed)


# ---------------------------------------------------------------------------
# Test 2: No stale "All stages 20-26 are complete"
# ---------------------------------------------------------------------------

def test_no_stale_stages_20_26():
    content = _read("ROADMAP_v0.2.md")
    passed = "All stages 20-26 are complete" not in content
    _record("2-no_stale_stages_20_26", passed)


# ---------------------------------------------------------------------------
# Test 3: No stale "manifest deferred"
# ---------------------------------------------------------------------------

def test_no_manifest_deferred():
    content = _read("docs/output-contract.md")
    passed = "manifest deferred" not in content.lower()
    _record("3-no_manifest_deferred", passed)


# ---------------------------------------------------------------------------
# Test 4: No stale "BigWig conversion not yet implemented"
# ---------------------------------------------------------------------------

def test_no_stale_bigwig_not_implemented():
    # Check active docs (not historical v0.1.0-beta CHANGELOG section)
    files_to_check = [
        "README.md",
        "KNOWN_ISSUES.md",
        "ROADMAP_v0.2.md",
        "docs/output-contract.md",
        "docs/qc-interpretation.md",
        "docs/configuration.md",
    ]
    failed = []
    for f in files_to_check:
        content = _read(f)
        if "BigWig conversion not yet implemented" in content:
            failed.append(f)
    passed = len(failed) == 0
    if not passed:
        print("  Stale in:", failed)
    _record("4-no_stale_bigwig_not_implemented", passed)


# ---------------------------------------------------------------------------
# Test 5: No stale "Stage 24a" in active docs
# ---------------------------------------------------------------------------

def test_no_stale_stage_24a():
    files_to_check = [
        "README.md",
        "KNOWN_ISSUES.md",
        "ROADMAP_v0.2.md",
        "docs/output-contract.md",
    ]
    failed = []
    for f in files_to_check:
        content = _read(f)
        if "Stage 24a" in content:
            failed.append(f)
    passed = len(failed) == 0
    if not passed:
        print("  Stale in:", failed)
    _record("5-no_stale_stage_24a", passed)


# ---------------------------------------------------------------------------
# Test 6: README contains key document links
# ---------------------------------------------------------------------------

def test_readme_contains_key_links():
    content = _read("README.md")
    required = [
        "ROADMAP_v0.2.md",
        "RELEASE_CHECKLIST.md",
        "docs/output-contract.md",
        "docs/assay-policy.md",
        "docs/idr-contract.md",
        "docs/release-checks/",
    ]
    missing = [s for s in required if s not in content]
    passed = len(missing) == 0
    if not passed:
        print("  Missing links:", missing)
    _record("6-readme_key_links", passed)


# ---------------------------------------------------------------------------
# Test 7: RELEASE_CHECKLIST contains Stage 22/24/25/27 test commands
# ---------------------------------------------------------------------------

def test_release_checklist_has_stage_tests():
    content = _read("RELEASE_CHECKLIST.md")
    required = [
        "test_stage22_bigwig_stress.py",
        "test_stage24_qc_summary_unit.py",
        "test_stage25_manifest_stress.py",
        "test_stage27_public_validation_plan.py",
        "test_stage27b_metadata_ci_plan.py",
        "test_stage27c_ci_workflow.py",
    ]
    missing = [s for s in required if s not in content]
    passed = len(missing) == 0
    if not passed:
        print("  Missing:", missing)
    _record("7-release_checklist_stage_tests", passed)


# ---------------------------------------------------------------------------
# Test 8: RELEASE_CHECKLIST does not contain expected-fail dry-run as must-pass
# ---------------------------------------------------------------------------

def test_release_checklist_no_fail_as_must_pass():
    content = _read("RELEASE_CHECKLIST.md")
    # Extract the first fenced bash code block after "## Automated Checks"
    auto_start = content.find("## Automated Checks")
    if auto_start < 0:
        _record("8-release_checklist_fail_note", False)
        return
    after_auto = content[auto_start:]
    fence_start = after_auto.find("```bash")
    if fence_start < 0:
        _record("8-release_checklist_fail_note", False)
        return
    fence_end = after_auto.find("```", fence_start + 7)
    if fence_end < 0:
        _record("8-release_checklist_fail_note", False)
        return
    code_block = after_auto[fence_start:fence_end]
    # The expected-fail dry-run must NOT be inside the automated checks code block
    dry_run_cmd = "snakemake -s workflow/Snakefile --configfile config/config.yaml -n --quiet"
    passed = dry_run_cmd not in code_block
    # And the explanatory note should still mention it's expected/non-blocking
    note_section = after_auto[fence_end:]
    passed = passed and ("expected to fail" in note_section.lower()
                        or "not a release-blocking" in note_section.lower())
    if not passed:
        print("  dry-run in code block: %s" % (dry_run_cmd in code_block))
        print("  note explains: %s" % ("expected to fail" in note_section.lower()))
    _record("8-release_checklist_fail_note", passed)


# ---------------------------------------------------------------------------
# Test 9: CHANGELOG [Unreleased] has no duplicate ### Added headings
# ---------------------------------------------------------------------------

def test_changelog_no_duplicate_headers():
    content = _read("CHANGELOG.md")
    errors = []

    # Both sections must exist
    if "## [Unreleased]" not in content:
        errors.append("[Unreleased] section missing")
    if "## [v0.2.0-rc1]" not in content:
        errors.append("[v0.2.0-rc1] section missing")

    def _extract_section(text, start_marker, end_marker=None):
        start = text.find(start_marker)
        if start < 0:
            return ""
        if end_marker:
            end = text.find(end_marker, start + len(start_marker))
            if end < 0:
                return text[start:]
            return text[start:end]
        return text[start:]

    def _count_duplicates(section_text):
        return section_text.count("### Added")

    # [Unreleased] section: between it and [v0.2.0-rc1]
    unreleased = _extract_section(content, "## [Unreleased]", "## [v0.2.0-rc1]")
    ur_added = _count_duplicates(unreleased)
    if ur_added > 1:
        errors.append("[Unreleased] has %d ### Added (max 1 allowed)" % ur_added)

    # [v0.2.0-rc1] section: between it and [v0.1.0-beta]
    rc1 = _extract_section(content, "## [v0.2.0-rc1]", "## [v0.1.0-beta]")
    rc1_added = _count_duplicates(rc1)
    if rc1_added != 1:
        errors.append("[v0.2.0-rc1] has %d ### Added (exactly 1 expected)" % rc1_added)

    passed = len(errors) == 0
    if not passed:
        for e in errors:
            print("  %s" % e)
    _record("9-changelog_no_dup_headers", passed)


# ---------------------------------------------------------------------------
# Test 10: ROADMAP contains 27a/27b/27c completed
# ---------------------------------------------------------------------------

def test_roadmap_has_27abc_completed():
    content = _read("ROADMAP_v0.2.md")
    required = ["27a", "27b", "27c", "Completed"]
    missing = [s for s in required if s not in content]
    passed = len(missing) == 0
    if not passed:
        print("  Missing:", missing)
    _record("10-roadmap_27abc_completed", passed)


# ---------------------------------------------------------------------------
# Test 11: RELEASE_CHECKLIST has public data execution as manual/external
# ---------------------------------------------------------------------------

def test_release_checklist_public_data_manual():
    content = _read("RELEASE_CHECKLIST.md")
    # Should either mention public data execution is manual or not list it as automated
    # No download commands like wget/curl in release checklist
    forbidden = ["wget", "curl", "fastq-dump", "prefetch"]
    found = [s for s in forbidden if s.lower() in content.lower()]
    passed = len(found) == 0
    if not passed:
        print("  Found download commands:", found)
    _record("11-release_checklist_no_downloads", passed)


# ---------------------------------------------------------------------------
def main():
    print("Starting Stage 28 Release Readiness Tests\n")

    test_no_stale_ci_status()
    test_no_stale_stages_20_26()
    test_no_manifest_deferred()
    test_no_stale_bigwig_not_implemented()
    test_no_stale_stage_24a()
    test_readme_contains_key_links()
    test_release_checklist_has_stage_tests()
    test_release_checklist_no_fail_as_must_pass()
    test_changelog_no_duplicate_headers()
    test_roadmap_has_27abc_completed()
    test_release_checklist_public_data_manual()

    print("\nSummary: %d/%d tests passed." % (PASSED, TOTAL))

    if PASSED < TOTAL:
        sys.exit(1)


if __name__ == "__main__":
    main()
