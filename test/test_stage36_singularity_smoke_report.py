#!/usr/bin/env python3
"""Stage 36 stress tests - Singularity runner smoke report."""

import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPORT = os.path.join(REPO_ROOT, "docs", "release-checks",
                      "stage36-singularity-runner-smoke.md")
README = os.path.join(REPO_ROOT, "containers", "README.md")
CHANGELOG = os.path.join(REPO_ROOT, "CHANGELOG.md")

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


def _read(path):
    with open(path) as fh:
        return fh.read()


def _tracked_files():
    try:
        result = subprocess.run(
            ["git", "ls-files"],
            cwd=REPO_ROOT,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        return [line.strip() for line in result.stdout.splitlines()
                if line.strip()]
    except Exception as exc:
        print("  Could not list git-tracked files:", exc)
        return []


# ---------------------------------------------------------------------------
# Test 1: Report exists
# ---------------------------------------------------------------------------

def test_report_exists():
    _record("1-report_exists", os.path.isfile(REPORT))


# ---------------------------------------------------------------------------
# Test 2: Report records SingularityCE version
# ---------------------------------------------------------------------------

def test_report_mentions_singularity_version():
    content = _read(REPORT)
    passed = ("singularity-ce version 4.1.1" in content
              and "singularity-container" in content)
    _record("2-report_singularity_version", passed)


# ---------------------------------------------------------------------------
# Test 3: Report records build command and definition file
# ---------------------------------------------------------------------------

def test_report_mentions_build_command():
    content = _read(REPORT)
    passed = ("sudo -E singularity build" in content
              and "containers/Apptainer.runner.def" in content)
    _record("3-report_build_command", passed)


# ---------------------------------------------------------------------------
# Test 4: Report records local SIF artifact name
# ---------------------------------------------------------------------------

def test_report_mentions_sif_name():
    content = _read(REPORT)
    passed = "chipseq-runner-stage36.sif" in content
    _record("4-report_sif_name", passed)


# ---------------------------------------------------------------------------
# Test 5: Report records pinned base image
# ---------------------------------------------------------------------------

def test_report_mentions_base_image():
    content = _read(REPORT)
    passed = "condaforge/miniforge3:24.11.3-0" in content
    _record("5-report_base_image", passed)


# ---------------------------------------------------------------------------
# Test 6: Report records Snakemake version
# ---------------------------------------------------------------------------

def test_report_mentions_snakemake_version():
    content = _read(REPORT)
    passed = "8.30.0" in content
    _record("6-report_snakemake_version", passed)


# ---------------------------------------------------------------------------
# Test 7: Report records PyYAML import smoke
# ---------------------------------------------------------------------------

def test_report_mentions_runner_ok():
    content = _read(REPORT)
    passed = "runner ok" in content.lower()
    _record("7-report_runner_ok", passed)


# ---------------------------------------------------------------------------
# Test 8: Report records bind-mounted dry-run result
# ---------------------------------------------------------------------------

def test_report_mentions_dry_run():
    content = _read(REPORT)
    passed = ("Bind-Mounted Dry-Run" in content
              and "total                         14" in content
              and "This was a dry-run" in content)
    _record("8-report_dry_run", passed)


# ---------------------------------------------------------------------------
# Test 9: Report documents HOME warning and correct guidance
# ---------------------------------------------------------------------------

def test_report_documents_home_warning():
    content = _read(REPORT)
    passed = ("Overriding HOME environment variable" in content
              and "--env HOME" in content
              and "--home" in content
              and "XDG_CACHE_HOME" in content)
    _record("9-report_home_warning", passed)


# ---------------------------------------------------------------------------
# Test 10: README links Stage 36 and no longer says Apptainer unexecuted
# ---------------------------------------------------------------------------

def test_readme_stage36_status():
    content = _read(README)
    passed = ("stage36-singularity-runner-smoke.md" in content
              and "Apptainer build is not executed yet" not in content)
    _record("10-readme_stage36_status", passed)


# ---------------------------------------------------------------------------
# Test 11: README Singularity guidance avoids --env HOME
# ---------------------------------------------------------------------------

def test_readme_singularity_home_guidance():
    content = _read(README)
    passed = ("--env XDG_CACHE_HOME=/conda_cache/xdg-cache" in content
              and "--env HOME=/conda_cache/home" not in content
              and "--home" in content)
    _record("11-readme_singularity_home_guidance", passed)


# ---------------------------------------------------------------------------
# Test 12: Changelog records Stage 36
# ---------------------------------------------------------------------------

def test_changelog_mentions_stage36():
    content = _read(CHANGELOG)
    passed = ("Singularity runner build smoke report (Stage 36)" in content
              and "singularity-ce version 4.1.1" in content)
    _record("12-changelog_stage36", passed)


# ---------------------------------------------------------------------------
# Test 13: No image artifacts are tracked by git
# ---------------------------------------------------------------------------

def test_no_tracked_image_artifacts():
    forbidden = (".sif", ".sqsh", ".img", ".oci",
                 ".tar", ".tar.gz", ".tar.zst", ".docker.tar")
    found = [path for path in _tracked_files() if path.endswith(forbidden)]
    passed = len(found) == 0
    if not passed:
        print("  Tracked image artifacts:", found)
    _record("13-no_tracked_image_artifacts", passed)


# ---------------------------------------------------------------------------
def main():
    print("Starting Stage 36 Singularity Runner Smoke Report Tests\n")

    test_report_exists()
    test_report_mentions_singularity_version()
    test_report_mentions_build_command()
    test_report_mentions_sif_name()
    test_report_mentions_base_image()
    test_report_mentions_snakemake_version()
    test_report_mentions_runner_ok()
    test_report_mentions_dry_run()
    test_report_documents_home_warning()
    test_readme_stage36_status()
    test_readme_singularity_home_guidance()
    test_changelog_mentions_stage36()
    test_no_tracked_image_artifacts()

    print("\nSummary: %d/%d tests passed." % (PASSED, TOTAL))

    pycache = os.path.join(REPO_ROOT, "scripts", "__pycache__")
    if os.path.isdir(pycache):
        import shutil
        shutil.rmtree(pycache, ignore_errors=True)

    if PASSED < TOTAL:
        sys.exit(1)


if __name__ == "__main__":
    main()
