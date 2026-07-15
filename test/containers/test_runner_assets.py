"""Durable static contracts for the supported workflow runner containers."""

from __future__ import annotations

from pathlib import Path
import re
import stat
import subprocess


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DOCKERFILE = PROJECT_ROOT / "containers" / "Dockerfile.runner"
APPTAINER_DEFINITION = PROJECT_ROOT / "containers" / "Apptainer.runner.def"
CONTAINER_README = PROJECT_ROOT / "containers" / "README.md"
USAGE_DOCUMENT = PROJECT_ROOT / "docs" / "container-usage.md"
SMOKE_SCRIPT = PROJECT_ROOT / "scripts" / "smoke_container_runner.sh"
RUNNER_PATH = "/opt/conda/envs/chipseq-runner/bin"
IMAGE_SUFFIXES = (
    ".sif",
    ".sqsh",
    ".img",
    ".oci",
    ".tar",
    ".tar.gz",
    ".tar.zst",
    ".docker.tar",
)


def _read(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def test_runner_definitions_share_a_pinned_base_and_runtime_path() -> None:
    dockerfile = _read(DOCKERFILE)
    apptainer = _read(APPTAINER_DEFINITION)
    docker_base = re.search(r"^FROM\s+(\S+)", dockerfile, re.MULTILINE)
    apptainer_base = re.search(
        r"^From:\s+(\S+)", apptainer, re.IGNORECASE | re.MULTILINE
    )

    assert docker_base is not None
    assert apptainer_base is not None
    assert docker_base.group(1) == apptainer_base.group(1)
    assert ":latest" not in docker_base.group(1)
    assert "<" not in docker_base.group(1)
    assert RUNNER_PATH in dockerfile
    assert RUNNER_PATH in apptainer


def test_runner_builds_use_explicit_locked_inputs_without_downloads() -> None:
    dockerfile = _read(DOCKERFILE)
    apptainer = _read(APPTAINER_DEFINITION)
    definitions = f"{dockerfile}\n{apptainer}"

    assert "workflow/envs/runner.yml" in dockerfile
    assert "conda env create -f /opt/pipeline/workflow/envs/runner.yml" in dockerfile
    assert "workflow/envs/runner.lock" in apptainer
    assert (
        "conda create -n chipseq-runner --file /opt/pipeline/workflow/envs/runner.lock"
        in apptainer
    )
    assert not re.search(r"^\s*(?:COPY|ADD)\s+\.\s", dockerfile, re.MULTILINE)
    assert not re.search(
        r"\b(?:curl|wget|prefetch|fastq-dump)\b", definitions, re.IGNORECASE
    )


def test_runner_entrypoints_preserve_the_workspace_contract() -> None:
    dockerfile = _read(DOCKERFILE)
    apptainer = _read(APPTAINER_DEFINITION)

    assert "WORKDIR /workspace" in dockerfile
    assert 'ENTRYPOINT ["snakemake"]' in dockerfile
    assert 'exec snakemake "$@"' in apptainer


def test_container_usage_documents_mount_and_cache_requirements() -> None:
    documentation = f"{_read(CONTAINER_README)}\n{_read(USAGE_DOCUMENT)}"

    for required in (
        "docker run",
        "Apptainer",
        "Singularity",
        "--conda-prefix",
        "/workspace",
        "/data",
        "/conda_cache",
        "HOME",
        "XDG_CACHE_HOME",
        "--home",
        "-u $(id -u):$(id -g)",
        "does not allow overriding `HOME` with `--env HOME=...`",
    ):
        assert required in documentation


def test_smoke_script_is_fail_closed_and_uses_an_ephemeral_workspace() -> None:
    script = _read(SMOKE_SCRIPT)

    assert SMOKE_SCRIPT.stat().st_mode & stat.S_IXUSR
    subprocess.run(["bash", "-n", SMOKE_SCRIPT], check=True)
    assert "set -euo pipefail" in script
    assert 'MODE="${1:-}"' in script
    assert 'MODE" != "docker"' in script
    assert 'MODE" != "singularity"' in script
    assert "mktemp -d /tmp/" in script
    assert "trap 'rm -rf \"$SMOKE_DIR\"' EXIT" in script
    assert "samples.tsv" in script
    assert "awk" in script
    assert not re.search(
        r"\b(?:curl|wget|prefetch|fastq-dump|docker build|singularity build)\b",
        script,
        re.IGNORECASE,
    )


def test_container_artifacts_stay_out_of_git_and_build_context() -> None:
    tracked = subprocess.run(
        ["git", "ls-files", "-z"],
        cwd=PROJECT_ROOT,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.split("\0")
    dockerignore = _read(PROJECT_ROOT / ".dockerignore").splitlines()
    gitignore = _read(PROJECT_ROOT / ".gitignore").splitlines()

    assert "scripts/smoke_container_runner.sh" in tracked
    assert not [path for path in tracked if path.endswith(IMAGE_SUFFIXES)]
    for pattern in (
        ".git/",
        ".snakemake/",
        "results/",
        "reference/",
        "*.fq.gz",
        "*.fastq.gz",
        "*.bam",
        "*.bw",
        "*.bdg",
        "*.sif",
        "*.tar",
    ):
        assert pattern in dockerignore
    for pattern in ("*.sif", "*.tar", "!scripts/smoke_container_runner.sh"):
        assert pattern in gitignore
