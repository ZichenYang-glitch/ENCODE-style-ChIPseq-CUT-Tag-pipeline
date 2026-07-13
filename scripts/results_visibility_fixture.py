"""Build the deterministic controlled project used by results visibility gates."""

from __future__ import annotations

import copy
from dataclasses import dataclass
from pathlib import Path
import shutil

import yaml

from encode_pipeline.adapters.encode_qc import _QC_HEADER


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_SENTINEL_NAME = ".encode-results-visibility-fixture"
FIXTURE_SENTINEL_CONTENT = "encode-results-visibility-fixture-v1\n"


@dataclass(frozen=True)
class ResultsVisibilityInputs:
    """Input documents and expected output shared by demo and browser gates."""

    samples_path: Path
    results_config: dict[str, object]
    cancel_config: dict[str, object]
    empty_config: dict[str, object]
    malformed_config: dict[str, object]
    expected_qc_summary: str


def prepare_results_visibility_fixture(
    project_root: Path,
    *,
    repository_root: Path = REPOSITORY_ROOT,
) -> ResultsVisibilityInputs:
    """Replace one sentinel-owned controlled project and return its inputs."""
    source_root = repository_root.expanduser().resolve()
    destination = _validate_destination(project_root, source_root)
    sentinel = destination / FIXTURE_SENTINEL_NAME
    if destination.exists():
        if (
            not sentinel.is_file()
            or sentinel.is_symlink()
            or sentinel.read_text(encoding="utf-8") != FIXTURE_SENTINEL_CONTENT
        ):
            raise ValueError("results visibility fixture directory is not owned")
        shutil.rmtree(destination)

    destination.mkdir(parents=True)
    sentinel.write_text(FIXTURE_SENTINEL_CONTENT, encoding="utf-8")
    try:
        _copy_project_contract(source_root, destination)
        expected_qc_summary = _build_qc_summary()
        _write_controlled_workflow(destination, expected_qc_summary)
        samples_path = (destination / "samples.tsv").resolve(strict=True)
        configs = tuple(
            _build_config(source_root, samples_path, threads)
            for threads in (1, 2, 3, 4)
        )
    except BaseException:
        shutil.rmtree(destination)
        raise

    return ResultsVisibilityInputs(
        samples_path=samples_path,
        results_config=configs[0],
        cancel_config=configs[1],
        empty_config=configs[2],
        malformed_config=configs[3],
        expected_qc_summary=expected_qc_summary,
    )


def _validate_destination(project_root: Path, repository_root: Path) -> Path:
    if not isinstance(project_root, Path):
        raise ValueError("fixture project root must be a Path")
    candidate = project_root.expanduser()
    if not candidate.is_absolute():
        candidate = Path.cwd() / candidate
    _reject_existing_symlink_components(candidate)
    destination = candidate.resolve(strict=False)
    home = Path.home().resolve()
    dangerous = {Path("/"), Path("/tmp"), home, repository_root}
    if destination in dangerous or destination.parent in dangerous:
        raise ValueError("fixture project root must be a dedicated child directory")
    if destination == repository_root or repository_root in destination.parents:
        allowed_parent = repository_root / ".local"
        if allowed_parent not in destination.parents:
            raise ValueError(
                "fixture project inside the repository must be under .local"
            )
    return destination


def _reject_existing_symlink_components(path: Path) -> None:
    current = path
    while current != current.parent:
        if current.is_symlink():
            raise ValueError("fixture project path must not contain symlinks")
        current = current.parent


def _copy_project_contract(source_root: Path, destination: Path) -> None:
    profile_root = source_root / "test" / "profiles" / "platform_worker_tiny"
    required_files = (
        source_root / "pyproject.toml",
        source_root / "docs" / "architecture" / "artifact-inventory.yaml",
        profile_root / "config.yaml",
        profile_root / "samples.tsv",
    )
    if not (source_root / "src" / "encode_pipeline").is_dir() or any(
        not path.is_file() or path.is_symlink() for path in required_files
    ):
        raise ValueError("repository does not contain the controlled fixture inputs")

    shutil.copy2(source_root / "pyproject.toml", destination / "pyproject.toml")
    shutil.copytree(
        source_root / "src" / "encode_pipeline",
        destination / "src" / "encode_pipeline",
    )
    inventory_target = destination / "docs" / "architecture"
    inventory_target.mkdir(parents=True)
    shutil.copy2(required_files[1], inventory_target / "artifact-inventory.yaml")
    profile_target = destination / "profiles" / "default"
    profile_target.mkdir(parents=True)
    profile_target.joinpath("config.yaml").write_text(
        "printshellcmds: true\ncores: 2\n",
        encoding="utf-8",
    )
    shutil.copy2(required_files[3], destination / "samples.tsv")


def _build_qc_summary() -> str:
    row = {column: "NA" for column in _QC_HEADER}
    row.update(
        {
            "sample": "C1",
            "assay": "chipseq",
            "total_reads": "1000",
            "frip": "0.125",
            "peak_count": "50",
            "percent_duplication": "0.2",
            "estimated_library_size": "9007199254740993",
            "nrf": "0.8",
            "pbc1": "0.75",
            "pbc2": "3.0",
        }
    )
    return (
        "\t".join(_QC_HEADER)
        + "\n"
        + "\t".join(row[column] for column in _QC_HEADER)
        + "\n"
    )


def _build_config(
    repository_root: Path,
    samples_path: Path,
    threads: int,
) -> dict[str, object]:
    profile_path = (
        repository_root / "test" / "profiles" / "platform_worker_tiny" / "config.yaml"
    )
    raw = yaml.safe_load(profile_path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("platform worker tiny config must be a mapping")
    config = copy.deepcopy(raw)
    config["samples"] = str(samples_path)
    config["outdir"] = "results"
    config["threads"] = threads
    return config


def _write_controlled_workflow(
    project_root: Path,
    expected_qc_summary: str,
) -> None:
    workflow_root = project_root / "workflow"
    scripts_root = project_root / "scripts"
    workflow_root.mkdir()
    scripts_root.mkdir()
    (workflow_root / "Snakefile").write_text(
        """
from pathlib import Path

HELPER = Path(workflow.basedir).parent / "scripts" / "results_visibility_task.py"
MODE_BY_THREADS = {1: "results", 2: "cancel", 3: "empty", 4: "malformed"}
MODE = MODE_BY_THREADS.get(int(config.get("threads", 1)))
if MODE is None:
    raise ValueError("controlled results fixture mode is invalid")

COMPLETE = "result/complete.txt"
MANIFEST = "results/multiqc/result_manifest.tsv"
QC_SUMMARY = "results/C1/01_qc/C1.qc_summary.tsv"
RESULT_OUTPUTS = [COMPLETE, MANIFEST]
if MODE in {"results", "malformed"}:
    RESULT_OUTPUTS.append(QC_SUMMARY)

rule all:
    input:
        RESULT_OUTPUTS

rule results_visibility_task:
    output:
        RESULT_OUTPUTS
    params:
        helper=str(HELPER),
        mode=MODE,
        complete=COMPLETE,
        manifest=MANIFEST,
        qc_summary=QC_SUMMARY
    shell:
        "python3 {params.helper:q} {params.mode:q} {params.complete:q} "
        "{params.manifest:q} {params.qc_summary:q}"
""".lstrip(),
        encoding="utf-8",
    )
    (scripts_root / "results_visibility_task.py").write_text(
        _task_script(expected_qc_summary),
        encoding="utf-8",
    )


def _task_script(expected_qc_summary: str) -> str:
    return f"""from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import time

EXPECTED_QC_SUMMARY = {expected_qc_summary!r}

mode = sys.argv[1]
complete = Path(sys.argv[2])
manifest = Path(sys.argv[3])
qc_summary = Path(sys.argv[4])
if mode not in {{"results", "cancel", "empty", "malformed"}}:
    raise SystemExit("controlled results fixture mode is invalid")
complete.parent.mkdir(parents=True, exist_ok=True)
manifest.parent.mkdir(parents=True, exist_ok=True)
qc_summary.parent.mkdir(parents=True, exist_ok=True)

if mode == "cancel":
    marker_root = Path(os.environ["TMPDIR"])
    child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(300)"])
    evidence = {{
        "shell_pid": os.getppid(),
        "helper_pid": os.getpid(),
        "child_pid": child.pid,
        "process_group": os.getpgrp(),
    }}
    temporary = marker_root / "browser-processes.json.tmp"
    temporary.write_text(json.dumps(evidence), encoding="utf-8")
    temporary.replace(marker_root / "browser-processes.json")
    print("browser-e2e-cancel-entered", flush=True)
    (marker_root / "browser-cancel-entered").write_text("entered\\n", encoding="utf-8")
    time.sleep(300)
    child.wait(timeout=5)
    complete.write_text("should-not-complete\\n", encoding="utf-8")
    manifest.write_text("should-not-complete\\n", encoding="utf-8")
    qc_summary.write_text("should-not-complete\\n", encoding="utf-8")
    (marker_root / "browser-helper-completed").write_text(
        "completed\\n", encoding="utf-8"
    )
    raise SystemExit(0)

print(f"browser-e2e-{{mode}}", flush=True)
complete.write_text("success\\n", encoding="utf-8")
has_qc = mode in {{"results", "malformed"}}
manifest_rows = [
    "output_type\\tstatus\\tpath",
]
if has_qc:
    manifest_rows.append(
        "qc_summary\\tpresent\\tresults/C1/01_qc/C1.qc_summary.tsv"
    )
manifest_rows.append(
    "result_manifest\\tpresent\\tresults/multiqc/result_manifest.tsv"
)
manifest.write_text("\\n".join(manifest_rows) + "\\n", encoding="utf-8")
if mode == "results":
    qc_summary.write_text(EXPECTED_QC_SUMMARY, encoding="utf-8")
elif mode == "malformed":
    qc_summary.write_text("sample\\tassay\\nC1\\tchipseq\\n", encoding="utf-8")
"""
