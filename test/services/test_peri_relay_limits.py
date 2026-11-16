"""Byte, connection, and deadline limits of the private relay."""

from __future__ import annotations

import socket
import threading

import pytest

from encode_pipeline.services.peri_runtime import relay as relay_module
from encode_pipeline.services.peri_runtime.relay import Relay, forward


def test_forward_transfers_until_disconnect_deadline_or_budget(tmp_path) -> None:
    upstream_path = tmp_path / "model.sock"
    listener = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    listener.bind(str(upstream_path))
    listener.listen(1)
    left, right = socket.socketpair()
    result: dict[str, str] = {}

    def run() -> None:
        result["value"] = forward(right, str(upstream_path), idle=1, total=2)

    worker = threading.Thread(target=run)
    worker.start()
    accepted, _address = listener.accept()
    try:
        left.sendall(b"ping")
        assert accepted.recv(4) == b"ping"
        accepted.sendall(b"pong")
        assert left.recv(4) == b"pong"
        left.close()
        worker.join(3)
        assert result["value"] == "disconnect"
    finally:
        accepted.close()
        listener.close()
        left.close()
        right.close()

    idle_left, idle_right = socket.socketpair()
    idle_listener = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    idle_path = tmp_path / "idle.sock"
    idle_listener.bind(str(idle_path))
    idle_listener.listen(1)
    idle_result: dict[str, str] = {}

    def run_idle() -> None:
        idle_result["value"] = forward(idle_right, str(idle_path), idle=0.05, total=1)

    idle_worker = threading.Thread(target=run_idle)
    idle_worker.start()
    idle_accepted, _address = idle_listener.accept()
    try:
        idle_worker.join(2)
        assert idle_result["value"] == "deadline"
    finally:
        idle_accepted.close()
        idle_listener.close()
        idle_left.close()
        idle_right.close()

    budget_left, budget_right = socket.socketpair()
    budget_listener = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    budget_path = tmp_path / "budget.sock"
    budget_listener.bind(str(budget_path))
    budget_listener.listen(1)
    budget_result: dict[str, str] = {}

    def run_budget() -> None:
        budget_result["value"] = forward(
            budget_right, str(budget_path), idle=2, total=5
        )

    budget_worker = threading.Thread(target=run_budget)
    budget_worker.start()
    budget_accepted, _address = budget_listener.accept()
    try:

        def drain() -> None:
            while True:
                try:
                    chunk = budget_accepted.recv(65536)
                except OSError:
                    return
                if not chunk:
                    return

        threading.Thread(target=drain, daemon=True).start()
        payload = b"x" * (relay_module.MAX_BYTES + 1)

        def send() -> None:
            try:
                budget_left.sendall(payload)
            except OSError:
                return

        sender = threading.Thread(target=send)
        sender.start()
        budget_worker.join(5)
        sender.join(5)
        assert budget_result["value"] == "byte_budget"
    finally:
        budget_accepted.close()
        budget_listener.close()
        budget_left.close()
        budget_right.close()


def test_relay_refuses_a_third_connection_and_closes(monkeypatch, tmp_path) -> None:
    held = threading.Event()
    release = threading.Event()

    def slow(client, *_args, **_kwargs) -> str:
        held.set()
        release.wait(2)
        client.close()
        return "held"

    monkeypatch.setattr(relay_module, "forward", slow)
    listener = socket.socket()
    listener.bind(("127.0.0.1", 0))
    listener.listen(4)
    port = listener.getsockname()[1]
    relay = Relay(listener, str(tmp_path / "missing.sock"))
    clients = []
    try:
        for _number in range(3):
            client = socket.socket()
            client.connect(("127.0.0.1", port))
            clients.append(client)
        assert held.wait(2)
        release.set()
    finally:
        release.set()
        for client in clients:
            client.close()
        relay.close()


def test_relay_main_stops_the_child(monkeypatch) -> None:
    monkeypatch.setattr(relay_module, "PORT", 0)

    class Child:
        def __init__(self) -> None:
            self.killed = False
            self.returncode = None

        def wait(self, timeout: float | None = None) -> int:
            del timeout
            return 0

        def poll(self) -> int | None:
            return None if not self.killed else 0

        def kill(self) -> None:
            self.killed = True

    child = Child()
    monkeypatch.setattr(relay_module.subprocess, "Popen", lambda *_a, **_k: child)
    assert relay_module.main() == 0
    assert child.killed is True


def test_relay_close_reports_a_stuck_worker(monkeypatch) -> None:
    listener = socket.socket()
    listener.bind(("127.0.0.1", 0))
    listener.listen(1)
    relay = Relay(listener, "/missing.sock")
    stuck = threading.Thread(target=lambda: threading.Event().wait(30), daemon=True)
    stuck.start()
    relay.threads.append(stuck)
    monkeypatch.setattr(relay_module, "TOTAL_SECONDS", 0)
    with pytest.raises(RuntimeError, match="relay close deadline"):
        relay.close()
