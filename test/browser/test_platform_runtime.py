from __future__ import annotations

import json
import os
from pathlib import Path
import signal
import subprocess
import sys

from scripts.run_local_platform import (
    RuntimeConfig,
    _session_process_groups,
    build_shared_environment,
    cleanup_service_sessions,
)
from scripts.results_visibility_fixture import prepare_results_visibility_fixture
from platform_runtime import write_manifest


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
OWNERSHIP_SENTINEL = ".encode-platform-playwright-owned"


def test_shared_environment_uses_one_absolute_lifecycle_configuration(tmp_path):
    project_root = tmp_path / "project"
    config = RuntimeConfig(
        project_root=project_root,
        frontend_root=tmp_path / "frontend",
        runtime_root=tmp_path / "runtime",
        redis_url="redis://127.0.0.1:6380/0",
        queue_name="browser-e2e",
        api_host="127.0.0.1",
        api_port=8010,
        frontend_host="127.0.0.1",
        frontend_port=4173,
        readiness_timeout=30,
    )

    environment = build_shared_environment(config, {"PATH": "/tools"})

    assert environment["ENCODE_PIPELINE_DATABASE_URL"] == (
        f"sqlite:///{tmp_path / 'runtime' / 'platform.db'}"
    )
    assert environment["ENCODE_PIPELINE_WORKSPACE_ROOT"] == str(
        tmp_path / "runtime" / "workspaces"
    )
    assert environment["ENCODE_PIPELINE_REDIS_URL"] == config.redis_url
    assert environment["ENCODE_PIPELINE_QUEUE_NAME"] == "browser-e2e"
    assert environment["PYTHONPATH"] == str(project_root / "src")
    assert environment["VITE_API_PROXY_TARGET"] == "http://127.0.0.1:8010"
    assert "127.0.0.1" in environment["NO_PROXY"]
    assert environment["NO_PROXY"] == environment["no_proxy"]


def test_session_process_groups_never_include_the_test_process_group(monkeypatch):
    own_group = 88
    monkeypatch.setattr("scripts.run_local_platform.Path.iterdir", lambda _path: ())
    monkeypatch.setattr("scripts.run_local_platform.os.getpgrp", lambda: own_group)

    groups = _session_process_groups(99)

    assert groups == ()
    assert own_group not in groups
    assert signal.SIGKILL.value > 0


def test_cleanup_service_sessions_signals_every_process_group(tmp_path, monkeypatch):
    session_id = 10001
    horse_group = 10002
    pid_file = tmp_path / "service-pids.json"
    pid_file.write_text(json.dumps({"worker": session_id}), encoding="utf-8")
    scans = iter(((horse_group, session_id), ()))
    monkeypatch.setattr(
        "scripts.run_local_platform._session_process_groups",
        lambda _session_id: next(scans, ()),
    )
    signals: list[tuple[int, signal.Signals]] = []
    monkeypatch.setattr(
        "scripts.run_local_platform.os.killpg",
        lambda group, signum: signals.append((group, signum)),
    )

    assert cleanup_service_sessions(pid_file) == 0
    assert (horse_group, signal.SIGTERM) in signals
    assert (session_id, signal.SIGTERM) in signals


def test_playwright_runtime_reset_rejects_a_mismatched_owner(tmp_path):
    runtime_root = tmp_path / "runtime"
    runtime_root.mkdir()
    unrelated = runtime_root / "unrelated.txt"
    unrelated.write_text("keep", encoding="utf-8")
    (runtime_root / OWNERSHIP_SENTINEL).write_text("owner-a", encoding="utf-8")
    environment = dict(os.environ)
    environment.update(
        {
            "ENCODE_PIPELINE_E2E_ROOT": str(runtime_root),
            "ENCODE_PIPELINE_E2E_OWNER": "owner-b",
        }
    )

    completed = subprocess.run(
        [sys.executable, str(REPOSITORY_ROOT / "test/browser/platform_runtime.py")],
        cwd=REPOSITORY_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )

    assert completed.returncode != 0
    assert "not owned" in completed.stderr
    assert unrelated.read_text(encoding="utf-8") == "keep"


def test_playwright_runtime_manifest_contains_only_controlled_fixture_inputs(
    tmp_path,
):
    runtime_root = tmp_path / "runtime"
    runtime_root.mkdir()
    inputs = prepare_results_visibility_fixture(runtime_root / "project")

    write_manifest(runtime_root, "browser-queue", inputs)

    raw = json.loads((runtime_root / "runtime.json").read_text(encoding="utf-8"))
    assert set(raw) == {
        "workflowId",
        "samplesPath",
        "resultsConfig",
        "cancelConfig",
        "emptyConfig",
        "malformedConfig",
        "expectedQcSummary",
        "runtimeRoot",
        "workspaceRoot",
        "markerRoot",
        "queueName",
    }
    assert raw["samplesPath"] == str(inputs.samples_path)
    assert raw["resultsConfig"] == inputs.results_config
    assert raw["cancelConfig"] == inputs.cancel_config
    assert raw["emptyConfig"] == inputs.empty_config
    assert raw["malformedConfig"] == inputs.malformed_config
    assert raw["expectedQcSummary"] == inputs.expected_qc_summary
    serialized = json.dumps(raw)
    assert "sqlite:" not in serialized
    assert "redis:" not in serialized
    assert "ENCODE_PIPELINE_" not in serialized
