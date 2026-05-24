#!/usr/bin/env python3
"""Stage 34 stress tests — runner container static files."""

import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONTAINERS_DIR = os.path.join(REPO_ROOT, "containers")
DOCKERFILE = os.path.join(CONTAINERS_DIR, "Dockerfile.runner")
APPTAINER_DEF = os.path.join(CONTAINERS_DIR, "Apptainer.runner.def")
README = os.path.join(CONTAINERS_DIR, "README.md")
DOCKERIGNORE = os.path.join(REPO_ROOT, ".dockerignore")

EXPECTED_FROM = "FROM condaforge/miniforge3:24.11.3-0"
EXPECTED_PATH = "/opt/conda/envs/chipseq-runner/bin"

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


# ---------------------------------------------------------------------------
# Test 1-3: Files exist
# ---------------------------------------------------------------------------

def test_dockerfile_exists():
    _record("1-dockerfile_exists", os.path.isfile(DOCKERFILE))


def test_apptainer_def_exists():
    _record("2-apptainer_def_exists", os.path.isfile(APPTAINER_DEF))


def test_readme_exists():
    _record("3-readme_exists", os.path.isfile(README))


# ---------------------------------------------------------------------------
# Test 4: Both files use identical FROM string
# ---------------------------------------------------------------------------

def test_identical_base_image():
    docker = _read(DOCKERFILE)
    apptainer = _read(APPTAINER_DEF)
    # Docker: "FROM condaforge/miniforge3:24.11.3-0"
    # Apptainer: "From: condaforge/miniforge3:24.11.3-0"
    passed = ("condaforge/miniforge3:24.11.3-0" in docker
              and "condaforge/miniforge3:24.11.3-0" in apptainer)
    _record("4-identical_base_image", passed)


# ---------------------------------------------------------------------------
# Test 5: No :latest or <pinned-version> in either file
# ---------------------------------------------------------------------------

def test_no_latest_or_placeholder():
    docker = _read(DOCKERFILE)
    apptainer = _read(APPTAINER_DEF)
    passed = (":latest" not in docker
              and ":latest" not in apptainer
              and "<pinned-version>" not in docker
              and "<pinned-version>" not in apptainer)
    if not passed:
        if ":latest" in docker or ":latest" in apptainer:
            print("  Found :latest")
        if "<pinned-version>" in docker or "<pinned-version>" in apptainer:
            print("  Found <pinned-version>")
    _record("5-no_latest_or_placeholder", passed)


# ---------------------------------------------------------------------------
# Test 6: Dockerfile does NOT contain COPY . . or COPY . or ADD .
# ---------------------------------------------------------------------------

def test_no_full_copy():
    docker = _read(DOCKERFILE)
    # Check each line independently to avoid false positives
    lines = docker.splitlines()
    full_copy = False
    for line in lines:
        stripped = line.strip()
        if stripped in ("COPY . .", "COPY .", "ADD .") or \
           stripped.startswith("COPY . ") or stripped.startswith("ADD . "):
            full_copy = True
            break
    passed = not full_copy
    if not passed:
        print("  Found full-context COPY/ADD pattern")
    _record("6-no_full_copy", passed)


# ---------------------------------------------------------------------------
# Test 7: Dockerfile copies workflow/envs/runner.yml
# ---------------------------------------------------------------------------

def test_dockerfile_copies_runner_yml():
    docker = _read(DOCKERFILE)
    passed = "workflow/envs/runner.yml" in docker
    _record("7-dockerfile_copies_runner_yml", passed)


# ---------------------------------------------------------------------------
# Test 8: Apptainer.def copies workflow/envs/runner.yml
# ---------------------------------------------------------------------------

def test_apptainer_copies_runner_yml():
    apptainer = _read(APPTAINER_DEF)
    passed = "workflow/envs/runner.yml" in apptainer
    _record("8-apptainer_copies_runner_yml", passed)


# ---------------------------------------------------------------------------
# Test 9: README mentions --conda-prefix
# ---------------------------------------------------------------------------

def test_readme_mentions_conda_prefix():
    readme = _read(README)
    passed = "--conda-prefix" in readme
    _record("9-readme_mentions_conda_prefix", passed)


# ---------------------------------------------------------------------------
# Test 10: README mentions Docker user mapping
# ---------------------------------------------------------------------------

def test_readme_mentions_docker_user():
    readme = _read(README)
    passed = "-u $(id -u):$(id -g)" in readme
    _record("10-readme_mentions_docker_user", passed)


# ---------------------------------------------------------------------------
# Test 11: README mentions --pwd /workspace (Apptainer)
# ---------------------------------------------------------------------------

def test_readme_mentions_apptainer_pwd():
    readme = _read(README)
    passed = "--pwd /workspace" in readme
    _record("11-readme_mentions_apptainer_pwd", passed)


# ---------------------------------------------------------------------------
# Test 12: README mentions bind-mount paths
# ---------------------------------------------------------------------------

def test_readme_mentions_bind_mounts():
    readme = _read(README)
    required = ["/workspace", "/data", "/reference", "/conda_cache"]
    missing = [p for p in required if p not in readme]
    passed = len(missing) == 0
    if not passed:
        print("  Missing bind-mount paths:", missing)
    _record("12-readme_bind_mounts", passed)


# ---------------------------------------------------------------------------
# Test 13: No wget/curl/fastq-dump in any container file
# ---------------------------------------------------------------------------

def test_no_download_commands():
    files = {
        "Dockerfile.runner": _read(DOCKERFILE),
        "Apptainer.runner.def": _read(APPTAINER_DEF),
        "README.md": _read(README),
    }
    forbidden = ["wget ", "curl ", "fastq-dump", "prefetch"]
    found = []
    for name, content in files.items():
        for cmd in forbidden:
            if cmd.lower() in content.lower():
                found.append(name)
    passed = len(found) == 0
    if not passed:
        print("  Download commands in:", found)
    _record("13-no_download_commands", passed)


# ---------------------------------------------------------------------------
# Test 14: .dockerignore excludes container artifacts and data files
# ---------------------------------------------------------------------------

def test_dockerignore_coverage():
    if not os.path.isfile(DOCKERIGNORE):
        _record("14-dockerignore_coverage", False)
        return
    content = _read(DOCKERIGNORE)
    required = [".git/", "*.bam", "*.sif"]
    missing = [p for p in required if p not in content]
    passed = len(missing) == 0
    if not passed:
        print("  Missing .dockerignore patterns:", missing)
    _record("14-dockerignore_coverage", passed)


# ---------------------------------------------------------------------------
# Test 15: Dockerfile has WORKDIR /workspace
# ---------------------------------------------------------------------------

def test_dockerfile_workdir():
    docker = _read(DOCKERFILE)
    passed = "WORKDIR /workspace" in docker
    _record("15-dockerfile_workdir", passed)


# ---------------------------------------------------------------------------
# Test 16: Dockerfile has ENTRYPOINT ["snakemake"]
# ---------------------------------------------------------------------------

def test_dockerfile_entrypoint():
    docker = _read(DOCKERFILE)
    passed = 'ENTRYPOINT ["snakemake"]' in docker
    _record("16-dockerfile_entrypoint", passed)


# ---------------------------------------------------------------------------
# Test 17: Dockerfile sets PATH to chipseq-runner env
# ---------------------------------------------------------------------------

def test_dockerfile_path():
    docker = _read(DOCKERFILE)
    passed = EXPECTED_PATH in docker
    _record("17-dockerfile_path", passed)


# ---------------------------------------------------------------------------
# Test 18: Apptainer %environment exports same PATH
# ---------------------------------------------------------------------------

def test_apptainer_env_path():
    apptainer = _read(APPTAINER_DEF)
    passed = EXPECTED_PATH in apptainer
    _record("18-apptainer_env_path", passed)


# ---------------------------------------------------------------------------
def main():
    print("Starting Stage 34 Runner Container Files Tests\n")

    test_dockerfile_exists()
    test_apptainer_def_exists()
    test_readme_exists()
    test_identical_base_image()
    test_no_latest_or_placeholder()
    test_no_full_copy()
    test_dockerfile_copies_runner_yml()
    test_apptainer_copies_runner_yml()
    test_readme_mentions_conda_prefix()
    test_readme_mentions_docker_user()
    test_readme_mentions_apptainer_pwd()
    test_readme_mentions_bind_mounts()
    test_no_download_commands()
    test_dockerignore_coverage()
    test_dockerfile_workdir()
    test_dockerfile_entrypoint()
    test_dockerfile_path()
    test_apptainer_env_path()

    print("\nSummary: %d/%d tests passed." % (PASSED, TOTAL))

    pycache = os.path.join(REPO_ROOT, "scripts", "__pycache__")
    if os.path.isdir(pycache):
        import shutil
        shutil.rmtree(pycache)

    if PASSED < TOTAL:
        sys.exit(1)


if __name__ == "__main__":
    main()
