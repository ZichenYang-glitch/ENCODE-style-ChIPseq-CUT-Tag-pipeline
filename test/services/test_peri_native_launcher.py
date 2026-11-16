"""Creation timeout ownership. These tests do not start Peri or systemd."""

from __future__ import annotations

from pathlib import Path
import queue
import subprocess
import threading
import time

import pytest

from encode_pipeline.services.peri_runtime import isolation
from encode_pipeline.services.peri_runtime.isolation import (
    _LaunchTicket,
    _wait_for_launch,
    native_launcher,
)
from encode_pipeline.services.peri_runtime.owner import NativeAuthoringBridge


_REAL_POPEN = subprocess.Popen


def _sleep() -> subprocess.Popen[bytes]:
    return _REAL_POPEN(
        ["/bin/sleep", "30"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        start_new_session=True,
    )


def _reap(process: subprocess.Popen[bytes] | None) -> None:
    if process is None or process.poll() is not None:
        return
    process.kill()
    process.wait(timeout=5)


def test_queued_timeout_does_not_start_the_process(monkeypatch) -> None:
    started = threading.Event()
    release = threading.Event()
    calls: list[list[str]] = []
    children: list[subprocess.Popen[bytes]] = []
    real_popen = subprocess.Popen

    def fake_popen(args, **kwargs):
        calls.append(list(args))
        if args == ["block"]:
            started.set()
            assert release.wait(5)
        child = _sleep()
        children.append(child)
        return child

    monkeypatch.setattr(isolation.subprocess, "Popen", fake_popen)
    launcher = native_launcher()
    holder: dict[str, object] = {}

    def run_first() -> None:
        try:
            holder["process"] = launcher.popen(["block"], timeout=5)
        except OSError as error:
            holder["error"] = error

    first = threading.Thread(target=run_first)
    second_error: dict[str, BaseException] = {}

    def run_second() -> None:
        try:
            launcher.popen(
                ["queued"],
                scope_unit="hw-peri-a1-queued",
                timeout=0.2,
            )
        except OSError as error:
            second_error["error"] = error

    second = threading.Thread(target=run_second)
    try:
        first.start()
        assert started.wait(5)
        second.start()
        second.join(2)
        assert second_error["error"].args == (
            "native launcher did not create the process",
        )
        assert calls == [["block"]]
        release.set()
        first.join(5)
        assert "error" not in holder
        assert calls == [["block"]]
    finally:
        release.set()
        first.join(5)
        monkeypatch.setattr(isolation.subprocess, "Popen", real_popen)
        for child in children:
            _reap(child)
        _reap(holder.get("process"))  # type: ignore[arg-type]


def test_inflight_timeout_reaps_the_process_and_scope(monkeypatch) -> None:
    entered = threading.Event()
    release = threading.Event()
    stopped: list[str | None] = []
    child: dict[str, subprocess.Popen[bytes]] = {}
    real_popen = subprocess.Popen

    def fake_popen(args, **kwargs):
        process = _sleep()
        child["process"] = process
        entered.set()
        assert release.wait(5)
        return process

    monkeypatch.setattr(isolation.subprocess, "Popen", fake_popen)
    monkeypatch.setattr(isolation, "stop_scope", lambda unit: stopped.append(unit))
    error: dict[str, BaseException] = {}

    def run() -> None:
        try:
            native_launcher().popen(
                ["/bin/sleep", "30"],
                scope_unit="hw-peri-a1-inflight",
                timeout=0.2,
            )
        except OSError as captured:
            error["error"] = captured

    worker = threading.Thread(target=run)
    try:
        worker.start()
        assert entered.wait(5)
        worker.join(2)
        assert error["error"].args == ("native launcher did not create the process",)
        release.set()
        deadline = time.monotonic() + 5
        while child["process"].poll() is None and time.monotonic() < deadline:
            time.sleep(0.02)
        assert child["process"].poll() is not None
        assert child["process"].stdin is not None and child["process"].stdin.closed
        assert child["process"].stdout is not None and child["process"].stdout.closed
        assert stopped == ["hw-peri-a1-inflight"]
    finally:
        release.set()
        worker.join(5)
        monkeypatch.setattr(isolation.subprocess, "Popen", real_popen)
        _reap(child.get("process"))


def test_delivered_process_survives_owner_shutdown(tmp_path: Path) -> None:
    process = native_launcher().popen(["/bin/sleep", "30"], timeout=5)
    bridge = NativeAuthoringBridge(binary=Path("/no/such/peri"), state_root=tmp_path)
    try:
        assert process.poll() is None
        bridge.shutdown()
        assert process.poll() is None
        assert bridge._sessions == {}
    finally:
        _reap(process)


def test_owner_shutdown_reaps_an_inflight_creation(tmp_path: Path, monkeypatch) -> None:
    entered = threading.Event()
    release = threading.Event()
    stopped: list[str | None] = []
    child: dict[str, subprocess.Popen[bytes]] = {}

    def fake_popen(args, **kwargs):
        process = _sleep()
        child["process"] = process
        entered.set()
        assert release.wait(5)
        return process

    monkeypatch.setattr(isolation.subprocess, "Popen", fake_popen)
    monkeypatch.setattr(isolation, "stop_scope", lambda unit: stopped.append(unit))
    error: dict[str, BaseException] = {}

    def run() -> None:
        try:
            native_launcher().popen(
                ["/bin/sleep", "30"],
                scope_unit="hw-peri-a1-shutdown",
                timeout=5,
            )
        except OSError as captured:
            error["error"] = captured

    worker = threading.Thread(target=run)
    bridge = NativeAuthoringBridge(binary=Path("/no/such/peri"), state_root=tmp_path)
    try:
        worker.start()
        assert entered.wait(5)
        bridge.shutdown()
        worker.join(2)
        assert error["error"].args == ("native launcher did not create the process",)
        release.set()
        deadline = time.monotonic() + 5
        while child["process"].poll() is None and time.monotonic() < deadline:
            time.sleep(0.02)
        assert child["process"].poll() is not None
        assert stopped == ["hw-peri-a1-shutdown"]
        follow = native_launcher().popen(["/bin/sleep", "30"], timeout=5)
        assert follow.poll() is None
        _reap(follow)
    finally:
        release.set()
        worker.join(5)
        _reap(child.get("process"))


def test_timeout_after_delivery_keeps_the_callers_process() -> None:
    process = _sleep()
    ticket = _LaunchTicket(lambda: process, "hw-peri-a1-delivered")
    assert ticket.claim_to_run()
    assert ticket.settle(process) is True

    original_get = ticket.box.get
    used = False

    def timeout_once(*args, **kwargs):
        nonlocal used
        if not used:
            used = True
            raise queue.Empty
        return original_get(*args, **kwargs)

    ticket.box.get = timeout_once  # type: ignore[method-assign]
    try:
        assert _wait_for_launch(ticket, 0.01) is process
        assert ticket.give_up() == "delivered"
        assert process.poll() is None
    finally:
        _reap(process)


def test_abandoned_running_ticket_does_not_return_the_process(monkeypatch) -> None:
    process = _sleep()
    ticket = _LaunchTicket(lambda: process, "hw-peri-a1-abandoned")
    stopped: list[str | None] = []
    assert ticket.claim_to_run()
    assert ticket.give_up() == "running"
    assert ticket.settle(process) is False
    monkeypatch.setattr(isolation, "stop_scope", lambda unit: stopped.append(unit))
    try:
        isolation._reap_unclaimed(process, ticket.scope_unit)
        assert process.poll() is not None
        assert stopped == ["hw-peri-a1-abandoned"]
    finally:
        _reap(process)


def test_abandoned_queued_ticket_never_starts() -> None:
    started = False

    def start() -> subprocess.Popen[bytes]:
        nonlocal started
        started = True
        return _sleep()

    ticket = _LaunchTicket(start, None)
    assert ticket.give_up() == "queued"
    assert ticket.claim_to_run() is False
    assert started is False


def test_launcher_thread_can_create_without_queueing_again() -> None:
    ticket = _LaunchTicket(
        lambda: native_launcher().popen(["/bin/sleep", "30"]),
        None,
    )
    native_launcher()._queue.put(ticket)
    process = ticket.box.get(timeout=5)
    try:
        assert isinstance(process, subprocess.Popen)
        assert process.poll() is None
    finally:
        _reap(process)


def test_spawn_records_scope_and_kills_when_qualification_fails(
    tmp_path: Path, monkeypatch
) -> None:
    process = _sleep()
    recorded: dict[str, object] = {}

    def fake_popen(args, **kwargs):
        recorded["unit"] = kwargs["scope_unit"]
        recorded["args"] = args
        return process

    monkeypatch.setattr(isolation._native_launcher, "popen", fake_popen)

    def unqualified(*_args, **_kwargs):
        raise OSError("native cgroup scope is not qualified")

    monkeypatch.setattr(isolation, "qualified_native_pid", unqualified)
    stopped: list[str] = []
    monkeypatch.setattr(isolation, "stop_scope", lambda unit: stopped.append(unit))
    monkeypatch.setenv("DBUS_SESSION_BUS_ADDRESS", "unix:path=/run/user/1000/bus")
    monkeypatch.setenv("XDG_RUNTIME_DIR", "/run/user/1000")
    binary = tmp_path / "peri"
    binary.write_bytes(b"peri")
    socket = tmp_path / "model.sock"
    socket.write_text("model", encoding="utf-8")
    try:
        with pytest.raises(OSError, match="not qualified"):
            isolation.spawn(binary, tmp_path / "state", socket, subprocess.DEVNULL)
        assert process.poll() is not None
        assert stopped == [recorded["unit"]]
        args = recorded["args"]
        assert isinstance(args, list)
        assert "--die-with-parent" in args
        assert "--property=MemoryMax=512M" in args
        assert "--property=MemorySwapMax=0" in args
        assert "--property=CPUQuota=100%" in args
        assert "--property=TasksMax=128" in args
    finally:
        _reap(process)


def test_spawn_returns_a_qualified_process(tmp_path: Path, monkeypatch) -> None:
    process = _sleep()

    def fake_popen(args, **kwargs):
        return process

    monkeypatch.setattr(isolation._native_launcher, "popen", fake_popen)
    monkeypatch.setattr(isolation, "qualified_native_pid", lambda *_a, **_k: 4321)
    monkeypatch.setenv("DBUS_SESSION_BUS_ADDRESS", "unix:path=/run/user/1000/bus")
    monkeypatch.setenv("XDG_RUNTIME_DIR", "/run/user/1000")
    binary = tmp_path / "peri"
    binary.write_bytes(b"peri")
    socket = tmp_path / "model.sock"
    socket.write_text("model", encoding="utf-8")
    try:
        created = isolation.spawn(
            binary, tmp_path / "state", socket, subprocess.DEVNULL
        )
        assert created is process
        assert created.peri_native_pid == 4321
        assert str(created.peri_scope_unit).startswith("hw-peri-a1-")
        assert isolation.isolation_paths_exist() in {True, False}
    finally:
        _reap(process)


def test_launcher_reports_its_thread_and_repeated_cancel() -> None:
    launcher = native_launcher()
    assert launcher.alive() is True
    assert isinstance(launcher.native_id(), int)
    ticket = _LaunchTicket(lambda: _sleep(), None)
    ticket.box.put(object())
    assert ticket.give_up() == "queued"
    assert ticket.give_up() == "abandoned"
    assert ticket.claim_to_run() is False


def test_delivered_result_missing_from_the_queue_is_an_error() -> None:
    ticket = _LaunchTicket(lambda: _sleep(), None)
    assert ticket.claim_to_run()
    ticket.phase = "settled"

    def empty(*_args, **_kwargs):
        raise queue.Empty

    ticket.box.get = empty  # type: ignore[method-assign]
    with pytest.raises(OSError, match="did not create"):
        _wait_for_launch(ticket, 0.01)


def test_launch_failure_returns_to_the_caller(monkeypatch) -> None:
    def fail_popen(*_args, **_kwargs):
        raise OSError("launch failed")

    monkeypatch.setattr(isolation.subprocess, "Popen", fail_popen)
    with pytest.raises(OSError, match="launch failed"):
        native_launcher().popen(["/bin/sleep", "30"], timeout=2)


def test_abandoned_non_process_is_discarded() -> None:
    started = threading.Event()
    release = threading.Event()

    def start():
        started.set()
        assert release.wait(5)
        return "not-a-process"

    ticket = _LaunchTicket(start, "hw-peri-a1-text")
    native_launcher()._queue.put(ticket)
    try:
        assert started.wait(5)
        assert ticket.give_up() == "running"
        release.set()
        deadline = time.monotonic() + 5
        while ticket.phase != "settled" and time.monotonic() < deadline:
            time.sleep(0.01)
        assert ticket.phase == "settled"
        assert ticket.abandoned is True
    finally:
        release.set()


def test_reap_tolerates_a_dead_process_and_closed_pipes(monkeypatch) -> None:
    class Streams:
        def close(self) -> None:
            raise OSError("close")

    class Process:
        stdin = None
        stdout = Streams()

        def poll(self) -> None:
            return None

        def kill(self) -> None:
            raise OSError("kill")

        def wait(self, timeout: float | None = None) -> None:
            raise subprocess.TimeoutExpired(cmd="sleep", timeout=timeout or 0)

    stopped: list[str | None] = []
    monkeypatch.setattr(isolation, "stop_scope", lambda unit: stopped.append(unit))
    isolation._reap_unclaimed(Process(), "hw-peri-a1-reap")  # type: ignore[arg-type]
    assert stopped == ["hw-peri-a1-reap"]


def test_spawn_continues_when_the_failed_child_does_not_exit(
    tmp_path: Path, monkeypatch
) -> None:
    process = _sleep()

    def fake_popen(args, **kwargs):
        return process

    def hang_wait(timeout: float | None = None) -> None:
        raise subprocess.TimeoutExpired(cmd="sleep", timeout=timeout or 0)

    real_wait = process.wait
    monkeypatch.setattr(isolation._native_launcher, "popen", fake_popen)

    def unqualified(*_args, **_kwargs):
        raise OSError("native cgroup scope is not qualified")

    monkeypatch.setattr(isolation, "qualified_native_pid", unqualified)
    monkeypatch.setattr(isolation, "stop_scope", lambda _unit: None)
    monkeypatch.setattr(process, "wait", hang_wait)
    monkeypatch.setenv("DBUS_SESSION_BUS_ADDRESS", "unix:path=/run/user/1000/bus")
    monkeypatch.setenv("XDG_RUNTIME_DIR", "/run/user/1000")
    binary = tmp_path / "peri"
    binary.write_bytes(b"peri")
    socket = tmp_path / "model.sock"
    socket.write_text("model", encoding="utf-8")
    try:
        with pytest.raises(OSError, match="not qualified"):
            isolation.spawn(binary, tmp_path / "state", socket, subprocess.DEVNULL)
    finally:
        process.wait = real_wait  # type: ignore[method-assign]
        _reap(process)


def test_spawn_requires_a_user_session(tmp_path: Path, monkeypatch) -> None:
    monkeypatch.delenv("DBUS_SESSION_BUS_ADDRESS", raising=False)
    monkeypatch.delenv("XDG_RUNTIME_DIR", raising=False)
    binary = tmp_path / "peri"
    binary.write_bytes(b"peri")
    socket = tmp_path / "model.sock"
    socket.write_text("model", encoding="utf-8")
    with pytest.raises(OSError, match="systemd user session"):
        isolation.spawn(binary, tmp_path / "state", socket, subprocess.DEVNULL)
