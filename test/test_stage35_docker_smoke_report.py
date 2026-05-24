#!/usr/bin/env python3
"""Stage 35 stress tests — Docker runner smoke report."""

import os
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPORT = os.path.join(REPO_ROOT, "docs", "release-checks",
                      "stage35-docker-runner-smoke.md")
README = os.path.join(REPO_ROOT, "containers", "README.md")

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
    passed = os.path.isfile(REPORT)
    _record("1-report_exists", passed)


# ---------------------------------------------------------------------------
# Test 2: Report mentions hello-world verification
# ---------------------------------------------------------------------------

def test_report_mentions_hello_world():
    content = _read(REPORT)
    passed = "hello-world" in content.lower()
    _record("2-report_hello_world", passed)


# ---------------------------------------------------------------------------
# Test 3: Report mentions miniforge3:24.11.3-0
# ---------------------------------------------------------------------------

def test_report_mentions_base_image():
    content = _read(REPORT)
    passed = "condaforge/miniforge3:24.11.3-0" in content
    _record("3-report_base_image", passed)


# ---------------------------------------------------------------------------
# Test 4: Report mentions chipseq-runner:stage35-smoke
# ---------------------------------------------------------------------------

def test_report_mentions_image_tag():
    content = _read(REPORT)
    passed = "chipseq-runner:stage35-smoke" in content
    _record("4-report_image_tag", passed)


# ---------------------------------------------------------------------------
# Test 5: Report mentions snakemake 8.30.0
# ---------------------------------------------------------------------------

def test_report_mentions_snakemake_version():
    content = _read(REPORT)
    passed = "8.30.0" in content
    _record("5-report_snakemake_version", passed)


# ---------------------------------------------------------------------------
# Test 6: Report mentions runner ok
# ---------------------------------------------------------------------------

def test_report_mentions_runner_ok():
    content = _read(REPORT)
    passed = "runner ok" in content.lower()
    _record("6-report_runner_ok", passed)


# ---------------------------------------------------------------------------
# Test 7: Report documents PermissionError /.cache fix
# ---------------------------------------------------------------------------

def test_report_documents_cache_fix():
    content = _read(REPORT)
    passed = ("/.cache" in content
              and "PermissionError" in content
              and "HOME" in content
              and "XDG_CACHE_HOME" in content)
    if not passed:
        missing = []
        if "/.cache" not in content:
            missing.append("/.cache")
        if "PermissionError" not in content:
            missing.append("PermissionError")
        if "HOME" not in content:
            missing.append("HOME")
        if "XDG_CACHE_HOME" not in content:
            missing.append("XDG_CACHE_HOME")
        print("  Missing:", missing)
    _record("7-report_cache_fix", passed)


# ---------------------------------------------------------------------------
# Test 8: README Docker example includes HOME and XDG_CACHE_HOME
# ---------------------------------------------------------------------------

def test_readme_docker_example_has_cache_env():
    content = _read(README)
    passed = ("HOME=/conda_cache/home" in content
              and "XDG_CACHE_HOME=/conda_cache/xdg-cache" in content)
    if not passed:
        if "HOME=/conda_cache/home" not in content:
            print("  Missing HOME env var")
        if "XDG_CACHE_HOME=/conda_cache/xdg-cache" not in content:
            print("  Missing XDG_CACHE_HOME env var")
    _record("8-readme_cache_env", passed)


# ---------------------------------------------------------------------------
# Test 9: No image artifacts tracked
# ---------------------------------------------------------------------------

def test_no_image_artifacts():
    forbidden = (".sif", ".sqsh", ".img", ".oci",
                 ".tar", ".tar.gz", ".tar.zst", ".docker.tar")
    found = [path for path in _tracked_files() if path.endswith(forbidden)]
    passed = len(found) == 0
    if not passed:
        print("  Tracked image artifacts:", found)
    _record("9-no_image_artifacts", passed)


# ---------------------------------------------------------------------------
def main():
    print("Starting Stage 35 Docker Runner Smoke Report Tests\n")

    test_report_exists()
    test_report_mentions_hello_world()
    test_report_mentions_base_image()
    test_report_mentions_image_tag()
    test_report_mentions_snakemake_version()
    test_report_mentions_runner_ok()
    test_report_documents_cache_fix()
    test_readme_docker_example_has_cache_env()
    test_no_image_artifacts()

    print("\nSummary: %d/%d tests passed." % (PASSED, TOTAL))

    pycache = os.path.join(REPO_ROOT, "scripts", "__pycache__")
    if os.path.isdir(pycache):
        import shutil
        shutil.rmtree(pycache, ignore_errors=True)

    if PASSED < TOTAL:
        sys.exit(1)


if __name__ == "__main__":
    main()
