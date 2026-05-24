#!/usr/bin/env python3
"""Stage 37 stress tests — container UX docs and smoke script."""

import os
import stat
import subprocess
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
USAGE_DOC = os.path.join(REPO_ROOT, "docs", "container-usage.md")
SMOKE_SCRIPT = os.path.join(REPO_ROOT, "scripts", "smoke_container_runner.sh")
CONTAINER_README = os.path.join(REPO_ROOT, "containers", "README.md")
README = os.path.join(REPO_ROOT, "README.md")

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
            cwd=REPO_ROOT, check=True,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
        )
        return [line.strip() for line in result.stdout.splitlines() if line.strip()]
    except Exception as exc:
        print("  Could not list git-tracked files:", exc)
        return []


# ---------------------------------------------------------------------------
# Test 1: docs/container-usage.md exists
# ---------------------------------------------------------------------------

def test_usage_doc_exists():
    _record("1-usage_doc_exists", os.path.isfile(USAGE_DOC))


# ---------------------------------------------------------------------------
# Test 2: Doc covers Docker
# ---------------------------------------------------------------------------

def test_usage_doc_covers_docker():
    content = _read(USAGE_DOC)
    _record("2-usage_doc_docker", "docker run" in content.lower())


# ---------------------------------------------------------------------------
# Test 3: Doc covers Singularity/Apptainer
# ---------------------------------------------------------------------------

def test_usage_doc_covers_singularity():
    content = _read(USAGE_DOC)
    passed = "apptainer" in content.lower() and "singularity" in content.lower()
    _record("3-usage_doc_singularity", passed)


# ---------------------------------------------------------------------------
# Test 4: Doc mentions --conda-prefix
# ---------------------------------------------------------------------------

def test_usage_doc_mentions_conda_prefix():
    content = _read(USAGE_DOC)
    _record("4-usage_doc_conda_prefix", "--conda-prefix" in content)


# ---------------------------------------------------------------------------
# Test 5: Doc mentions bind mounts
# ---------------------------------------------------------------------------

def test_usage_doc_mentions_bind_mounts():
    content = _read(USAGE_DOC)
    _record("5-usage_doc_bind_mounts",
            "/workspace" in content and "/data" in content and "/conda_cache" in content)


# ---------------------------------------------------------------------------
# Test 6: Doc mentions HOME/XDG_CACHE_HOME
# ---------------------------------------------------------------------------

def test_usage_doc_mentions_cache_env():
    content = _read(USAGE_DOC)
    passed = "HOME" in content and "XDG_CACHE_HOME" in content
    _record("6-usage_doc_cache_env", passed)


# ---------------------------------------------------------------------------
# Test 7: Doc covers SingularityCE HOME warning / --home
# ---------------------------------------------------------------------------

def test_usage_doc_covers_home_warning():
    content = _read(USAGE_DOC)
    passed = "--home" in content and "HOME" in content
    _record("7-usage_doc_home_warning", passed)


# ---------------------------------------------------------------------------
# Test 8: Doc has troubleshooting section
# ---------------------------------------------------------------------------

def test_usage_doc_troubleshooting():
    content = _read(USAGE_DOC)
    passed = ("Troubleshooting" in content
              and "PermissionError" in content
              and "conda_cache" in content.lower())
    _record("8-usage_doc_troubleshooting", passed)


# ---------------------------------------------------------------------------
# Test 9: Doc mentions deferred publishing
# ---------------------------------------------------------------------------

def test_usage_doc_deferred_publishing():
    content = _read(USAGE_DOC)
    passed = ("deferred" in content.lower()
              or "Image Publishing" in content
              or "publishing is deferred" in content.lower())
    _record("9-usage_doc_deferred_publishing", passed)


# ---------------------------------------------------------------------------
# Test 10: Smoke script exists and is executable
# ---------------------------------------------------------------------------

def test_smoke_script_executable():
    passed = os.path.isfile(SMOKE_SCRIPT) and os.access(SMOKE_SCRIPT, os.X_OK)
    _record("10-smoke_script_executable", passed)


# ---------------------------------------------------------------------------
# Test 11: Smoke script has docker and singularity modes
# ---------------------------------------------------------------------------

def test_smoke_script_docker_singularity():
    content = _read(SMOKE_SCRIPT)
    passed = ("docker" in content and "singularity" in content
              and "MODE" in content)
    _record("11-smoke_script_modes", passed)


# ---------------------------------------------------------------------------
# Test 12: Smoke script uses /tmp workspace
# ---------------------------------------------------------------------------

def test_smoke_script_tmp_workspace():
    content = _read(SMOKE_SCRIPT)
    passed = "/tmp" in content and "mktemp" in content
    _record("12-smoke_script_tmp", passed)


# ---------------------------------------------------------------------------
# Test 13: Smoke script no download commands
# ---------------------------------------------------------------------------

def test_smoke_script_no_downloads():
    content = _read(SMOKE_SCRIPT)
    forbidden = ["wget ", "curl ", "fastq-dump", "prefetch"]
    passed = not any(cmd in content for cmd in forbidden)
    if not passed:
        found = [cmd for cmd in forbidden if cmd in content]
        print("  Found:", found)
    _record("13-smoke_script_no_downloads", passed)


# ---------------------------------------------------------------------------
# Test 14: Smoke script does not write .sif/.tar artifacts
# ---------------------------------------------------------------------------

def test_smoke_script_no_image_write():
    content = _read(SMOKE_SCRIPT)
    # Should not create .sif or .tar files
    passed = "docker build" not in content and "singularity build" not in content
    _record("14-smoke_script_no_image_write", passed)


# ---------------------------------------------------------------------------
# Test 15: No image artifacts tracked by git
# ---------------------------------------------------------------------------

def test_smoke_script_parses_sample_sheet():
    """Smoke script must create FASTQs matching the profile's samples.tsv."""
    content = _read(SMOKE_SCRIPT)
    # Must parse the sample sheet (awk / tail) not hardcode filenames
    passed = ("awk" in content and "samples.tsv" in content)
    if not passed:
        print("  Script should parse samples.tsv for FASTQ filenames")
    _record("16-smoke_script_parses_sample_sheet", passed)


def test_gitignore_unignores_smoke_script():
    """scripts/smoke_container_runner.sh must not be git-ignored."""
    result = subprocess.run(
        ["git", "check-ignore", "scripts/smoke_container_runner.sh"],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )
    # Exit code 0 means it IS ignored (bad); exit code 1 means NOT ignored (good)
    passed = result.returncode != 0
    if not passed:
        print("  Script is git-ignored! Output:", result.stdout.strip())
    _record("17-gitignore_unignores_smoke_script", passed)


def test_no_image_artifacts():
    forbidden = (".sif", ".sqsh", ".img", ".oci",
                 ".tar", ".tar.gz", ".tar.zst", ".docker.tar")
    found = [path for path in _tracked_files() if path.endswith(forbidden)]
    passed = len(found) == 0
    if not passed:
        print("  Tracked image artifacts:", found)
    _record("15-no_image_artifacts", passed)


# ---------------------------------------------------------------------------
def main():
    print("Starting Stage 37 Container UX Tests\n")

    test_usage_doc_exists()
    test_usage_doc_covers_docker()
    test_usage_doc_covers_singularity()
    test_usage_doc_mentions_conda_prefix()
    test_usage_doc_mentions_bind_mounts()
    test_usage_doc_mentions_cache_env()
    test_usage_doc_covers_home_warning()
    test_usage_doc_troubleshooting()
    test_usage_doc_deferred_publishing()
    test_smoke_script_executable()
    test_smoke_script_docker_singularity()
    test_smoke_script_tmp_workspace()
    test_smoke_script_no_downloads()
    test_smoke_script_no_image_write()
    test_smoke_script_parses_sample_sheet()
    test_gitignore_unignores_smoke_script()
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
