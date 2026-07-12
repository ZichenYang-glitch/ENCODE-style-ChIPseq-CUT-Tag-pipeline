from __future__ import annotations

import signal

from scripts.run_local_platform import (
    RuntimeConfig,
    _session_process_groups,
    build_shared_environment,
)


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

    assert groups == (99,)
    assert own_group not in groups
    assert signal.SIGKILL.value > 0
