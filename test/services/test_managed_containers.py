"""Tests for the server-owned local Docker cleanup boundary."""

from __future__ import annotations

from contextlib import contextmanager
import json
from pathlib import Path
import stat
import sys
import time
from unittest.mock import patch

import pytest

from encode_pipeline.platform.managed_containers import (
    MANAGED_CONTAINER_SCOPE_LABEL,
    managed_container_endpoint_identity,
    managed_container_scope,
)
from encode_pipeline.services.managed_containers import ManagedContainerCleaner
import encode_pipeline.services.managed_containers as managed_containers_module


_CONTAINER_A = "a" * 64
_CONTAINER_B = "b" * 64


@contextmanager
def _fake_docker(
    tmp_path: Path,
    *,
    state: dict[str, object] | None = None,
    cleanup_timeout_seconds: float = 0.4,
    maximum_cleanup_timeout_seconds: float = 2.0,
    final_cleanup_timeout_seconds: float = 0.4,
):
    state_path = tmp_path / "docker-state.json"
    calls_path = tmp_path / "docker-calls.jsonl"
    configured = {
        "all": [],
        "running": [],
        "stubborn": False,
        "mode": "ok",
        "all_ps_count": 0,
        "late_on_all_ps": None,
        "late_id": None,
        **(state or {}),
    }
    state_path.write_text(json.dumps(configured), encoding="utf-8")
    executable = tmp_path / "docker"
    executable.write_text(
        f"""#!{sys.executable}
import json
import os
from pathlib import Path
import sys
import time

state_path = Path({str(state_path)!r})
calls_path = Path({str(calls_path)!r})
args = sys.argv[1:]
with calls_path.open("a", encoding="utf-8") as handle:
    handle.write(json.dumps({{"argv": args, "docker_host": os.environ.get("DOCKER_HOST")}}) + "\\n")
state = json.loads(state_path.read_text(encoding="utf-8"))
command = args[2]
if state["mode"] == "hang":
    time.sleep(60)
if state["mode"] == "fail":
    print("private socket diagnostic", file=sys.stderr)
    raise SystemExit(7)
if command == "ps":
    if state["mode"] == "malformed":
        print("not-a-container-id")
    else:
        if "--all" in args:
            state["all_ps_count"] += 1
            if state["all_ps_count"] == state["late_on_all_ps"]:
                late_id = state["late_id"]
                if late_id not in state["all"]:
                    state["all"].append(late_id)
                if late_id not in state["running"]:
                    state["running"].append(late_id)
            values = state["all"]
        else:
            values = state["running"]
        print("\\n".join(values))
elif command == "stop":
    if not state["stubborn"]:
        state["running"] = []
elif command == "kill":
    state["running"] = []
elif command == "rm":
    removed = set(args[4:])
    state["all"] = [value for value in state["all"] if value not in removed]
state_path.write_text(json.dumps(state), encoding="utf-8")
""",
        encoding="utf-8",
    )
    executable.chmod(0o755)
    socket_path = tmp_path / "docker.sock"
    socket_path.write_bytes(b"")
    real_lstat = managed_containers_module.os.lstat
    socket_stat = real_lstat(socket_path)
    socket_values = list(socket_stat)
    socket_values[0] = stat.S_IFSOCK | 0o660
    simulated_socket_stat = managed_containers_module.os.stat_result(socket_values)

    def lstat_with_simulated_socket(path):
        if Path(path) == socket_path:
            return simulated_socket_stat
        return real_lstat(path)

    with patch.object(
        managed_containers_module.os,
        "lstat",
        lstat_with_simulated_socket,
    ):
        cleaner = ManagedContainerCleaner(
            executable=executable,
            unix_socket=socket_path,
            cleanup_timeout_seconds=cleanup_timeout_seconds,
            maximum_cleanup_timeout_seconds=maximum_cleanup_timeout_seconds,
            final_cleanup_timeout_seconds=final_cleanup_timeout_seconds,
            stop_timeout_seconds=1,
        )
        yield cleaner, calls_path, state_path, executable, socket_path


def _calls(path: Path) -> list[dict[str, object]]:
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]


def test_managed_container_scope_is_deterministic_and_path_private(tmp_path):
    workspace = (tmp_path / "workspaces/run-1").resolve()

    first = managed_container_scope(workspace)

    assert first == managed_container_scope(workspace)
    assert len(first) == 64
    assert str(workspace) not in first
    assert managed_container_scope(workspace.parent / "run-2") != first
    assert MANAGED_CONTAINER_SCOPE_LABEL == "org.helixweave.workspace-scope"


def test_managed_container_endpoint_identity_is_deterministic_and_path_private(
    tmp_path,
):
    executable = (tmp_path / "bin/docker").resolve()
    socket = (tmp_path / "run/docker.sock").resolve()

    identity = managed_container_endpoint_identity(executable, socket)

    assert identity == managed_container_endpoint_identity(executable, socket)
    assert len(identity) == 64
    assert str(tmp_path) not in identity
    assert identity != managed_container_endpoint_identity(
        executable,
        (tmp_path / "other/docker.sock").resolve(),
    )


def test_cleaner_runs_bounded_stop_relist_kill_remove_sequence(tmp_path, monkeypatch):
    monkeypatch.setenv("DOCKER_HOST", "tcp://remote.example:2375")
    with _fake_docker(
        tmp_path,
        state={
            "all": [_CONTAINER_A, _CONTAINER_B],
            "running": [_CONTAINER_B],
            "stubborn": True,
        },
    ) as (cleaner, calls_path, state_path, _executable, socket_path):
        scope = managed_container_scope((tmp_path / "workspace").resolve())
        started = time.monotonic()

        result = cleaner.cleanup(scope)
        elapsed = time.monotonic() - started

        assert result.is_success
        calls = _calls(calls_path)
        commands = [call["argv"][2] for call in calls]
        assert commands[:6] == [
            "ps",
            "stop",
            "ps",
            "kill",
            "ps",
            "rm",
        ]
        assert commands[6:]
        assert set(commands[6:]) == {"ps"}
        assert elapsed >= cleaner.cleanup_timeout_seconds
        assert all(call["docker_host"] is None for call in calls)
        assert all(
            call["argv"][:2] == ["--host", f"unix://{socket_path}"] for call in calls
        )
        assert all(
            f"label={MANAGED_CONTAINER_SCOPE_LABEL}={scope}" in call["argv"]
            for call in calls
            if call["argv"][2] == "ps"
        )
        assert json.loads(state_path.read_text(encoding="utf-8"))["all"] == []


def test_cleaner_catches_container_registered_after_initial_empty_observation(
    tmp_path,
):
    with _fake_docker(
        tmp_path,
        state={
            "late_on_all_ps": 4,
            "late_id": _CONTAINER_A,
        },
    ) as (cleaner, calls_path, state_path, _executable, _socket_path):
        result = cleaner.cleanup("c" * 64)

        assert result.is_success
        calls = _calls(calls_path)
        commands = [call["argv"][2] for call in calls]
        assert commands[:4] == ["ps", "ps", "ps", "ps"]
        assert "stop" in commands[4:]
        assert "rm" in commands[4:]
        assert commands[-2:] == ["ps", "ps"]
        state = json.loads(state_path.read_text(encoding="utf-8"))
        assert state["all"] == []
        assert state["running"] == []


@pytest.mark.parametrize("mode", ("malformed", "fail"))
def test_cleaner_fails_closed_without_cli_diagnostics(tmp_path, mode):
    with _fake_docker(tmp_path, state={"mode": mode}) as (
        cleaner,
        _calls_path,
        _state_path,
        _executable,
        socket_path,
    ):
        result = cleaner.cleanup("c" * 64)

        assert result.is_failure
        rendered = str(result.issues[0].to_dict())
        assert result.issues[0].code == "MANAGED_CONTAINER_CLEANUP_FAILED"
        assert str(socket_path) not in rendered
        assert "private socket diagnostic" not in rendered


def test_cleaner_enforces_bounded_observation_and_final_phase_budgets(tmp_path):
    with _fake_docker(
        tmp_path,
        state={"mode": "hang"},
        cleanup_timeout_seconds=0.2,
        maximum_cleanup_timeout_seconds=0.2,
        final_cleanup_timeout_seconds=0.2,
    ) as (cleaner, _calls_path, _state_path, _executable, _socket_path):
        started = time.monotonic()

        result = cleaner.cleanup("d" * 64)

        elapsed = time.monotonic() - started
        assert result.is_failure
        assert result.issues[0].code == "MANAGED_CONTAINER_CLEANUP_FAILED"
        assert elapsed < 1.0


def test_cleaner_rejects_symlink_executable_and_non_socket(tmp_path):
    executable = tmp_path / "docker-real"
    executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    executable.chmod(0o755)
    symlink = tmp_path / "docker"
    symlink.symlink_to(executable)
    not_socket = tmp_path / "docker.sock"
    not_socket.write_text("not a socket", encoding="utf-8")

    with pytest.raises(ValueError, match="endpoints"):
        ManagedContainerCleaner(executable=symlink, unix_socket=not_socket)


def test_cleaner_detects_endpoint_replacement(tmp_path):
    with _fake_docker(tmp_path) as (
        cleaner,
        _calls_path,
        _state_path,
        executable,
        _socket_path,
    ):
        replacement = tmp_path / "replacement"
        replacement.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
        replacement.chmod(0o755)
        executable.unlink()
        replacement.rename(executable)

        verification = cleaner.verify_endpoint()
        result = cleaner.cleanup("e" * 64)

        assert verification.is_failure
        assert verification.issues[0].code == "MANAGED_CONTAINER_ENDPOINT_CHANGED"
        assert result.is_failure
        assert result.issues[0].code == "MANAGED_CONTAINER_CLEANUP_FAILED"
