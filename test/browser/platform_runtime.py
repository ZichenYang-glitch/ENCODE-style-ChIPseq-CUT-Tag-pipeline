#!/usr/bin/env python3
"""Prepare a controlled tiny workflow and exec the real local platform stack."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import sys
from uuid import uuid4

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from scripts.results_visibility_fixture import (  # noqa: E402
    ResultsVisibilityInputs,
    prepare_results_visibility_fixture,
)


OWNERSHIP_SENTINEL = ".encode-platform-playwright-owned"


def write_manifest(
    runtime_root: Path,
    queue_name: str,
    inputs: ResultsVisibilityInputs,
) -> None:
    """Write the bounded browser-visible fixture contract."""
    manifest = {
        "workflowId": "encode-style-chipseq-cuttag-atac-mnase",
        "samplesPath": str(inputs.samples_path),
        "resultsConfig": inputs.results_config,
        "cancelConfig": inputs.cancel_config,
        "emptyConfig": inputs.empty_config,
        "malformedConfig": inputs.malformed_config,
        "expectedQcSummary": inputs.expected_qc_summary,
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
    inputs = prepare_results_visibility_fixture(project_root)
    queue_name = f"encode-pipeline-browser-{uuid4().hex}"
    write_manifest(runtime_root, queue_name, inputs)
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
