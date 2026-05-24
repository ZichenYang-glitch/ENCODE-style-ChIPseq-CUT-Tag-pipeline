#!/usr/bin/env python3
"""Stage 33 stress tests — containerization plan."""

import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RELEASE_CHECKS = os.path.join(REPO_ROOT, "docs", "release-checks")
SPEC = os.path.join(REPO_ROOT, "docs", "superpowers", "specs",
                    "2026-05-24-stage33-containerization-planning-design.md")
PLAN = os.path.join(REPO_ROOT, "docs", "superpowers", "plans",
                    "2026-05-24-stage33-containerization-planning.md")
REL_DOC = os.path.join(RELEASE_CHECKS, "stage33-containerization-plan.md")

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
# Test 1: Spec exists
# ---------------------------------------------------------------------------

def test_spec_exists():
    passed = os.path.isfile(SPEC)
    _record("1-spec_exists", passed)


# ---------------------------------------------------------------------------
# Test 2: Plan exists
# ---------------------------------------------------------------------------

def test_plan_exists():
    passed = os.path.isfile(PLAN)
    _record("2-plan_exists", passed)


# ---------------------------------------------------------------------------
# Test 3: Release checks doc exists
# ---------------------------------------------------------------------------

def test_rel_doc_exists():
    passed = os.path.isfile(REL_DOC)
    _record("3-rel_doc_exists", passed)


# ---------------------------------------------------------------------------
# Test 4: Plan mentions Docker
# ---------------------------------------------------------------------------

def test_plan_mentions_docker():
    with open(REL_DOC) as fh:
        content = fh.read()
    passed = "Docker" in content
    _record("4-plan_mentions_docker", passed)


# ---------------------------------------------------------------------------
# Test 5: Plan mentions Apptainer/Singularity
# ---------------------------------------------------------------------------

def test_plan_mentions_apptainer():
    with open(REL_DOC) as fh:
        content = fh.read()
    passed = "Apptainer" in content
    _record("5-plan_mentions_apptainer", passed)


# ---------------------------------------------------------------------------
# Test 6: Plan mentions runner image
# ---------------------------------------------------------------------------

def test_plan_mentions_runner_image():
    with open(REL_DOC) as fh:
        content = fh.read()
    passed = "runner" in content.lower() and "image" in content.lower()
    _record("6-plan_mentions_runner_image", passed)


# ---------------------------------------------------------------------------
# Test 7: Plan mentions full-tool image tradeoff
# ---------------------------------------------------------------------------

def test_plan_mentions_full_tool():
    with open(REL_DOC) as fh:
        content = fh.read()
    passed = "full-tool" in content.lower() or "full bioinformatics" in content.lower()
    _record("7-plan_mentions_full_tool", passed)


# ---------------------------------------------------------------------------
# Test 8: Plan mentions --use-conda
# ---------------------------------------------------------------------------

def test_plan_mentions_use_conda():
    with open(REL_DOC) as fh:
        content = fh.read()
    passed = "--use-conda" in content
    _record("8-plan_mentions_use_conda", passed)


# ---------------------------------------------------------------------------
# Test 9: Plan mentions CI scope / no builds
# ---------------------------------------------------------------------------

def test_plan_mentions_ci_scope():
    with open(REL_DOC) as fh:
        content = fh.read()
    passed = ("CI" in content
              and ("no builds" in content.lower()
                   or "no images" in content.lower()
                   or "not in this stage" in content.lower()
                   or "planning" in content.lower()))
    _record("9-plan_mentions_ci_scope", passed)


# ---------------------------------------------------------------------------
# Test 10: No public data downloads
# ---------------------------------------------------------------------------

def test_no_public_data_downloads():
    with open(REL_DOC) as fh:
        content = fh.read()
    forbidden = ["wget", "curl", "fastq-dump", "prefetch", "ENCSR", "GSE"]
    found = [f for f in forbidden if f.lower() in content.lower()]
    # ENCSR and GSE might appear in context — that's fine
    # Only fail if download commands are present
    download_cmds = ["wget", "curl", "fastq-dump", "prefetch"]
    found_dl = [f for f in download_cmds if f.lower() in content.lower()]
    passed = len(found_dl) == 0
    if not passed:
        print("  Found download commands:", found_dl)
    _record("10-no_public_data_downloads", passed)


# ---------------------------------------------------------------------------
# Test 11: No Docker image tarballs or large binary artifacts committed
# ---------------------------------------------------------------------------

_CONTAINER_ARTIFACT_SUFFIXES = {
    ".sif", ".sqsh", ".img", ".oci",
    ".tar", ".tar.gz", ".tar.zst", ".docker.tar",
}


def test_no_container_artifacts():
    found = []
    for root, dirs, files in os.walk(REPO_ROOT):
        dirs[:] = [d for d in dirs if not d.startswith(".")]
        for f in files:
            for ext in _CONTAINER_ARTIFACT_SUFFIXES:
                if f.endswith(ext):
                    found.append(os.path.relpath(
                        os.path.join(root, f), REPO_ROOT))
    passed = len(found) == 0
    if not passed:
        print("  Found:", found)
    _record("11-no_container_artifacts", passed)


def test_gitignore_covers_container_artifacts():
    gitignore_path = os.path.join(REPO_ROOT, ".gitignore")
    with open(gitignore_path) as fh:
        content = fh.read()
    # .gitignore patterns use *<suffix> to match any filename with that suffix
    patterns = {"*" + ext for ext in _CONTAINER_ARTIFACT_SUFFIXES}
    passed = all(p in content for p in patterns)
    if not passed:
        missing = [p for p in patterns if p not in content]
        print("  Missing .gitignore patterns:", missing)
    _record("12-gitignore_covers_container_artifacts", passed)


# ---------------------------------------------------------------------------
# Test 13: Plan mentions --conda-prefix
# ---------------------------------------------------------------------------

def test_plan_mentions_conda_prefix():
    with open(REL_DOC) as fh:
        content = fh.read()
    passed = "--conda-prefix" in content
    _record("13-plan_mentions_conda_prefix", passed)


# ---------------------------------------------------------------------------
# Test 14: Plan mentions writable cache / /conda_cache
# ---------------------------------------------------------------------------

def test_plan_mentions_writable_cache():
    with open(REL_DOC) as fh:
        content = fh.read()
    passed = "/conda_cache" in content and "writable" in content.lower()
    _record("14-plan_mentions_writable_cache", passed)


# ---------------------------------------------------------------------------
# Test 15: Plan does NOT contain :latest tag
# ---------------------------------------------------------------------------

def test_plan_no_latest_tags():
    with open(REL_DOC) as fh:
        content = fh.read()
    passed = ":latest" not in content
    if not passed:
        print("  Found ':latest' in plan doc")
    _record("15-plan_no_latest_tags", passed)


# ---------------------------------------------------------------------------
# Test 16: Plan mentions bind-mount / repo must be bind-mounted
# ---------------------------------------------------------------------------

def test_plan_mentions_bind_mount():
    with open(REL_DOC) as fh:
        content = fh.read()
    passed = ("bind-mount" in content.lower()
              or "bind mount" in content.lower()
              or ("workspace" in content.lower() and "bind" in content.lower()))
    _record("16-plan_mentions_bind_mount", passed)


# ---------------------------------------------------------------------------
# Test 17: .dockerignore exists
# ---------------------------------------------------------------------------

def test_dockerignore_exists():
    di_path = os.path.join(REPO_ROOT, ".dockerignore")
    passed = os.path.isfile(di_path)
    _record("17-dockerignore_exists", passed)


# ---------------------------------------------------------------------------
# Test 18: .dockerignore covers data and container artifacts
# ---------------------------------------------------------------------------

def test_dockerignore_covers_data_and_artifacts():
    di_path = os.path.join(REPO_ROOT, ".dockerignore")
    if not os.path.isfile(di_path):
        _record("18-dockerignore_covers_data_and_artifacts", False)
        return
    with open(di_path) as fh:
        content = fh.read()
    required = ["*.bam", "*.bw", "*.bdg", "*.fq.gz", "*.fastq.gz",
                ".sif", ".git/"]
    missing = [p for p in required if p not in content]
    passed = len(missing) == 0
    if not passed:
        print("  Missing .dockerignore patterns:", missing)
    _record("18-dockerignore_covers_data_and_artifacts", passed)


# ---------------------------------------------------------------------------
# Test 19: Apptainer example includes --pwd /workspace
# ---------------------------------------------------------------------------

def test_apptainer_example_pwd():
    with open(REL_DOC) as fh:
        content = fh.read()
    # Find the Apptainer usage code block and verify --pwd is present
    passed = "--pwd /workspace" in content
    _record("19-apptainer_example_pwd", passed)


# ---------------------------------------------------------------------------
def main():
    print("Starting Stage 33 Containerization Plan Tests\n")

    test_spec_exists()
    test_plan_exists()
    test_rel_doc_exists()
    test_plan_mentions_docker()
    test_plan_mentions_apptainer()
    test_plan_mentions_runner_image()
    test_plan_mentions_full_tool()
    test_plan_mentions_use_conda()
    test_plan_mentions_ci_scope()
    test_no_public_data_downloads()
    test_no_container_artifacts()
    test_gitignore_covers_container_artifacts()
    test_plan_mentions_conda_prefix()
    test_plan_mentions_writable_cache()
    test_plan_no_latest_tags()
    test_plan_mentions_bind_mount()
    test_dockerignore_exists()
    test_dockerignore_covers_data_and_artifacts()
    test_apptainer_example_pwd()

    print("\nSummary: %d/%d tests passed." % (PASSED, TOTAL))

    pycache = os.path.join(REPO_ROOT, "scripts", "__pycache__")
    if os.path.isdir(pycache):
        import shutil
        shutil.rmtree(pycache)

    if PASSED < TOTAL:
        sys.exit(1)


if __name__ == "__main__":
    main()
