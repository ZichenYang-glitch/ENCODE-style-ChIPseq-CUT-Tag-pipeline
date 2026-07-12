#!/usr/bin/env python3
"""Prepare a controlled tiny workflow and exec the real local platform stack."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import sys
from uuid import uuid4

import yaml


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
PROFILE_ROOT = REPOSITORY_ROOT / "test" / "profiles" / "platform_worker_tiny"
OWNERSHIP_SENTINEL = ".encode-platform-playwright-owned"


def create_controlled_project(project_root: Path) -> None:
    project_root.mkdir(parents=True)
    shutil.copy2(REPOSITORY_ROOT / "pyproject.toml", project_root / "pyproject.toml")
    shutil.copytree(
        REPOSITORY_ROOT / "src" / "encode_pipeline",
        project_root / "src" / "encode_pipeline",
    )
    workflow_root = project_root / "workflow"
    profile_root = project_root / "profiles" / "default"
    scripts_root = project_root / "scripts"
    workflow_root.mkdir()
    profile_root.mkdir(parents=True)
    scripts_root.mkdir()
    (profile_root / "config.yaml").write_text(
        "printshellcmds: true\ncores: 2\n", encoding="utf-8"
    )
    (workflow_root / "Snakefile").write_text(
        """
from pathlib import Path

HELPER = Path(workflow.basedir).parent / "scripts" / "browser_task.py"
MODE = "cancel" if int(config.get("threads", 1)) == 2 else "success"

rule all:
    input:
        "result/complete.txt"

rule browser_task:
    output:
        "result/complete.txt"
    params:
        helper=str(HELPER),
        mode=MODE
    shell:
        "python3 {params.helper:q} {params.mode:q} {output:q}"
""".lstrip(),
        encoding="utf-8",
    )
    (scripts_root / "browser_task.py").write_text(
        """
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import time

mode = sys.argv[1]
output = Path(sys.argv[2])
output.parent.mkdir(parents=True, exist_ok=True)
if mode == "success":
    print("browser-e2e-success", flush=True)
    output.write_text("success\\n", encoding="utf-8")
    raise SystemExit(0)

marker_root = Path(os.environ["TMPDIR"])
child = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(300)"])
evidence = {
    "shell_pid": os.getppid(),
    "helper_pid": os.getpid(),
    "child_pid": child.pid,
    "process_group": os.getpgrp(),
}
temporary = marker_root / "browser-processes.json.tmp"
temporary.write_text(json.dumps(evidence), encoding="utf-8")
temporary.replace(marker_root / "browser-processes.json")
print("browser-e2e-cancel-entered", flush=True)
(marker_root / "browser-cancel-entered").write_text("entered\\n", encoding="utf-8")
time.sleep(300)
child.wait(timeout=5)
output.write_text("should-not-complete\\n", encoding="utf-8")
(marker_root / "browser-helper-completed").write_text("completed\\n", encoding="utf-8")
""".lstrip(),
        encoding="utf-8",
    )


def write_manifest(runtime_root: Path, project_root: Path, queue_name: str) -> None:
    base_config = yaml.safe_load(
        (PROFILE_ROOT / "config.yaml").read_text(encoding="utf-8")
    )
    samples_path = project_root / "samples.tsv"
    shutil.copy2(PROFILE_ROOT / "samples.tsv", samples_path)
    base_config["samples"] = str(samples_path)
    success_config = dict(base_config)
    success_config["threads"] = 1
    success_config["outdir"] = "browser-success"
    cancel_config = dict(base_config)
    cancel_config["threads"] = 2
    cancel_config["outdir"] = "browser-cancel"
    manifest = {
        "workflowId": "encode-style-chipseq-cuttag-atac-mnase",
        "samplesPath": str(samples_path),
        "successConfig": success_config,
        "cancelConfig": cancel_config,
        "runtimeRoot": str(runtime_root),
        "workspaceRoot": str(runtime_root / "workspaces"),
        "markerRoot": str(runtime_root / "tmp"),
        "queueName": queue_name,
    }
    (runtime_root / "runtime.json").write_text(
        json.dumps(manifest, sort_keys=True), encoding="utf-8"
    )


def prepare_owned_runtime_root(runtime_value: str, runtime_owner: str) -> Path:
    """Reset only the invocation-owned temporary runtime directory."""
    runtime_root = Path(runtime_value).resolve()
    if runtime_root in {Path("/"), Path("/tmp"), Path.home(), REPOSITORY_ROOT}:
        raise ValueError("E2E runtime root must be a dedicated temporary directory")
    sentinel = runtime_root / OWNERSHIP_SENTINEL
    if not sentinel.is_file() or sentinel.read_text(encoding="utf-8") != runtime_owner:
        raise ValueError("E2E runtime root is not owned by this Playwright invocation")
    shutil.rmtree(runtime_root)
    runtime_root.mkdir(parents=True)
    (runtime_root / OWNERSHIP_SENTINEL).write_text(runtime_owner, encoding="utf-8")
    return runtime_root


def main() -> None:
    runtime_value = os.environ.get("ENCODE_PIPELINE_E2E_ROOT")
    runtime_owner = os.environ.get("ENCODE_PIPELINE_E2E_OWNER")
    if not runtime_value or not runtime_owner:
        raise ValueError("Playwright runtime root and owner token are required")
    runtime_root = prepare_owned_runtime_root(runtime_value, runtime_owner)
    project_root = runtime_root / "project"
    create_controlled_project(project_root)
    queue_name = f"encode-pipeline-browser-{uuid4().hex}"
    write_manifest(runtime_root, project_root, queue_name)
    redis_url = os.environ.get(
        "ENCODE_PIPELINE_E2E_REDIS_URL", "redis://127.0.0.1:6380/0"
    )
    argv = [
        sys.executable,
        str(REPOSITORY_ROOT / "scripts" / "run_local_platform.py"),
        "--project-root",
        str(project_root),
        "--frontend-root",
        str(REPOSITORY_ROOT / "frontend"),
        "--runtime-root",
        str(runtime_root),
        "--redis-url",
        redis_url,
        "--queue-name",
        queue_name,
        "--api-port",
        "8010",
        "--frontend-port",
        "4173",
        "--readiness-timeout",
        "120",
    ]
    os.execvpe(sys.executable, argv, os.environ.copy())


if __name__ == "__main__":
    main()
