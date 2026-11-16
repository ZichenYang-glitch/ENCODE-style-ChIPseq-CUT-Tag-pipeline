"""Authorization and failure paths of one native session, without Peri."""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable
import json
import os
from pathlib import Path
import subprocess
import threading
import time
from types import SimpleNamespace

import pytest

from encode_pipeline.services.peri_runtime.identity import SERVER_NAME, WIRE_TOOLS
from encode_pipeline.services.peri_runtime.policy import process_policy
from encode_pipeline.services.peri_runtime.session import NativeRejected, NativeSession


def _sleep() -> subprocess.Popen[bytes]:
    return subprocess.Popen(
        ["/bin/sleep", "30"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        start_new_session=True,
    )


def _shell(process: subprocess.Popen[bytes] | None = None) -> NativeSession:
    session = NativeSession.__new__(NativeSession)
    session.gate = threading.RLock()
    session.write_lock = threading.Lock()
    session.pending_lock = threading.Lock()
    session.pending = {}
    session.request_id = 0
    session.state = "blocked"
    session.approve = False
    session.session_id = None
    session.server_id = f"{SERVER_NAME}-test"
    session.handler_calls = 0
    session.tool_arrivals = 0
    session.approvals = []
    session.denied_callbacks = []
    session.errors = []
    session.grants = Counter()
    session.used_requests = set()
    session.assistant_parts = []
    session.events = []
    session.event_bytes = 0
    session.active_turn = False
    session.policy = process_policy()
    session.qualification_egress = False
    session.last_tool_ok = False
    session.last_tool_payload = None
    session.tool_payloads = ()
    session.native_pid = None
    session.leader_pid = None
    session.scope_unit = "hw-peri-a1-session"
    session.scope_cgroup = None
    session.closed_resources = False
    session.reader = None
    session.stub = None
    session.log = None
    session.socket_dir = None
    session.authorization_check = None
    session.bootstrap_payload = b"\x00\x00\x00\x02{}"
    session.process = process
    session.tool_handler = lambda name, arguments, owner: {"ok": True, "name": name}
    return session


def _stop(process: subprocess.Popen[bytes] | None) -> None:
    if process is None or process.poll() is not None:
        return
    process.kill()
    process.wait(timeout=5)


def _facts(session: NativeSession, suppressed: int = 0) -> dict[str, object]:
    return {
        "version": 1,
        "prediction_enabled": False,
        "auto_compact_enabled": False,
        "tool_policy": session.policy,
        "prediction_started": 0,
        "prediction_suppressed": suppressed,
    }


def _catalog() -> dict[str, object]:
    return {
        "servers": [
            {
                "name": SERVER_NAME,
                "connectionStatus": "connected",
                "tools": [{"name": name} for name in WIRE_TOOLS],
            }
        ]
    }


def test_callbacks_reject_closed_mismatched_and_unknown_calls() -> None:
    session = _shell()
    session.state = "closed"
    assert session._callback("mcp/message", {})["error"]["code"] == -32602
    session.state = "ready"
    assert session._callback("mcp/message", "nope")["error"]["code"] == -32602
    session.session_id = "owned"
    rejected = session._callback("mcp/message", {"sessionId": "other"})
    assert rejected["error"]["message"] == "session rejected"
    unknown = session._callback("session/other", {})
    assert unknown["error"]["code"] == -32601
    discovered = session._mcp({"method": "server/discover", "params": {}})
    assert discovered["result"]["result"]["ttlMs"] == 0
    assert session._mcp({"method": "initialize", "params": {}})["result"]["result"]
    assert session._mcp({"method": "tools/list", "params": []})["result"]["result"]
    assert session._mcp({"method": "ping", "params": {}})["result"]["result"] == {}
    hidden = session._mcp({"method": "unknown", "params": {}})
    assert hidden["result"]["error"]["code"] == -32601


def test_permission_and_tool_call_follow_the_grant() -> None:
    session = _shell()
    session.state = "ready"
    session.approve = True
    session.active_turn = True
    name = f"mcp__{SERVER_NAME}__get_authoring_schema"
    denied = session._permission({"toolCall": {"title": "other", "rawInput": {}}})
    assert denied["result"]["outcome"]["outcome"] == "cancelled"
    allowed = session._permission(
        {
            "toolCall": {"title": name, "rawInput": {}},
            "options": [{"kind": "allow_once", "optionId": "allow"}],
        }
    )
    assert allowed["result"]["outcome"]["optionId"] == "allow"
    unauthorized = session._call_tool(
        {"requestId": "req-1"},
        {"name": "get_authoring_schema", "arguments": {"extra": True}},
    )
    assert unauthorized["result"]["error"]["code"] == -32602
    accepted = session._call_tool(
        {"requestId": "req-1"},
        {"name": "get_authoring_schema", "arguments": {}},
    )
    assert accepted["result"]["result"]["isError"] is False
    assert session.handler_calls == 1

    def fail_handler(*_args):
        raise RuntimeError("tool")

    session.tool_handler = fail_handler
    from encode_pipeline.services.peri_runtime.session import grant_key

    session.grants = Counter({grant_key(name, {}): 1})
    failed = session._call_tool(
        {"requestId": "req-2"},
        {"name": "get_authoring_schema", "arguments": {}},
    )
    assert failed["result"]["error"]["code"] == -32603


def test_frames_reads_and_rpc_failures(tmp_path: Path, monkeypatch) -> None:
    process = _sleep()
    session = _shell(process)
    try:
        with pytest.raises(NativeRejected, match="ACP_MALFORMED"):
            session.consume_line(b"not-json\n")
        with pytest.raises(NativeRejected, match="ACP_MALFORMED"):
            session.consume_line(b'{"jsonrpc":"1.0"}\n')
        with pytest.raises(NativeRejected, match="ACP_FRAME_TOO_LARGE"):
            session.consume_line(b"x" * (1024 * 1024 + 1))
        session.events = [{}] * 2000
        with pytest.raises(NativeRejected, match="ACP_EVENT_BUDGET"):
            session.consume_line(b'{"jsonrpc":"2.0","id":1}\n')
        session.events = []
        session.event_bytes = 0
        session.consume_line(
            json.dumps(
                {
                    "jsonrpc": "2.0",
                    "method": "session/update",
                    "params": {
                        "update": {
                            "sessionUpdate": "agent_message_chunk",
                            "content": {"type": "text", "text": "hello"},
                        }
                    },
                }
            ).encode()
        )
        session.consume_line(
            json.dumps(
                {
                    "jsonrpc": "2.0",
                    "method": "session/update",
                    "params": {
                        "update": {
                            "sessionUpdate": "agent_message",
                            "content": "there",
                        }
                    },
                }
            ).encode()
        )
        assert session.assistant_parts == ["hello", "there"]
        reply: list[dict[str, object]] = []
        session.pending[7] = type("Box", (), {"put": reply.append})()
        session.consume_line(b'{"jsonrpc":"2.0","id":7,"result":{"ok":true}}\n')
        assert reply[0]["result"] == {"ok": True}
        session.state = "closed"
        with pytest.raises(NativeRejected, match="CLOSED"):
            session.rpc("ping", {})
        session.state = "ready"

        def reject_remote(value) -> None:
            session.pending[value["id"]].put({"error": {"message": "no"}})

        session._send = reject_remote
        with pytest.raises(NativeRejected, match="ACP_REMOTE_REJECTION"):
            session.rpc("ping", {})
        session.state = "ready"
        session.errors.clear()

        def reject_shape(value) -> None:
            session.pending[value["id"]].put({"result": []})

        session._send = reject_shape
        with pytest.raises(NativeRejected, match="ACP_REMOTE_REJECTION"):
            session.rpc("ping", {})
        session.state = "ready"
        session.errors.clear()
        session._send = lambda value: None
        with pytest.raises(NativeRejected, match="RPC_TIMEOUT"):
            session.rpc("ping", {}, timeout=0.01)
        assert "RPC_TIMEOUT" in session.errors
        session.state = "ready"
        session.errors.clear()
        session._send = lambda value: (_ for _ in ()).throw(OSError("closed"))
        with pytest.raises(NativeRejected, match="RPC_TRANSPORT_FAILURE"):
            session.rpc("ping", {})
        read, write = os.pipe()
        os.set_blocking(write, False)
        session.process.stdin = os.fdopen(write, "wb", buffering=0)
        session.state = "ready"
        session._write_bytes(b"ok")
        assert os.read(read, 8) == b"ok"
        os.close(read)
        with pytest.raises(NativeRejected, match="ACP_WRITE_BUDGET"):
            session._write_bytes(b"x" * 262145)
        session.process.stdin = None
        with pytest.raises(NativeRejected, match="ACP_STDIN_CLOSED"):
            session._write_bytes(b"x")
        monkeypatch.setattr(
            "encode_pipeline.services.peri_runtime.session.select.select",
            lambda *_args, **_kwargs: ([], [], []),
        )
        session.process.stdin = os.fdopen(os.open(os.devnull, os.O_WRONLY), "wb")
        with pytest.raises(NativeRejected, match="ACP_WRITE_TIMEOUT"):
            session._write_bytes(b"x")
    finally:
        _stop(process)


def test_reader_eof_and_close_release_the_process(tmp_path: Path, monkeypatch) -> None:
    process = _sleep()
    session = _shell(process)
    log_path = tmp_path / "stderr.log"
    session.log = log_path.open("wb")
    session.socket_dir = tmp_path / "sock"
    session.socket_dir.mkdir()
    session.stub = SimpleNamespace(close=lambda: setattr(session, "stub_closed", True))
    stopped: list[str | None] = []
    monkeypatch.setattr(
        "encode_pipeline.services.peri_runtime.session.stop_scope",
        lambda unit: stopped.append(unit),
    )
    stdout = tmp_path / "out"
    stdout.write_bytes(b'{"jsonrpc":"2.0","method":"session/update","params":{}}\n')
    session.process.stdout = stdout.open("rb")
    session._read()
    assert "ACP_EOF" in session.errors
    session.state = "ready"
    session.errors.clear()
    session.process.stdout = type(
        "Broken",
        (),
        {
            "readline": lambda self, _limit: (_ for _ in ()).throw(OSError("broken")),
            "close": lambda self: None,
        },
    )()
    session._read()
    assert "ACP_TRANSPORT_FAILURE" in session.errors
    session.close()
    assert stopped == ["hw-peri-a1-session"]
    assert process.poll() is not None
    assert session.stub_closed is True
    session.close()
    assert session.closed_resources is True


def test_runtime_checks_and_initialization(monkeypatch) -> None:
    process = _sleep()
    session = _shell(process)
    monkeypatch.setattr(
        "encode_pipeline.services.peri_runtime.session.stop_scope",
        lambda _unit: None,
    )
    try:
        session.rpc = lambda method, params, timeout=15: (
            _facts(session) | {"version": 2}
        )
        with pytest.raises(NativeRejected, match="HOST_POLICY_MISMATCH"):
            session.verify_runtime()
        observed = {"ticks": 0}

        def clock() -> float:
            observed["ticks"] += 1
            return 0 if observed["ticks"] < 4 else 10

        monkeypatch.setattr(
            "encode_pipeline.services.peri_runtime.session.time.monotonic",
            clock,
        )
        monkeypatch.setattr(
            "encode_pipeline.services.peri_runtime.session.time.sleep",
            lambda _seconds: None,
        )
        session.state = "blocked"
        session.errors.clear()
        session.rpc = lambda method, params, timeout=15: _facts(session, suppressed=0)
        with pytest.raises(NativeRejected, match="TURN_END_POLICY_NOT_OBSERVED"):
            session.verify_runtime(minimum_suppressed=1)
        monkeypatch.setattr(
            "encode_pipeline.services.peri_runtime.session.time.monotonic",
            time.monotonic,
        )
        monkeypatch.setattr(
            "encode_pipeline.services.peri_runtime.session.time.sleep",
            time.sleep,
        )
        _stop(process)
        process = _sleep()
        session.process = process
        session.state = "blocked"
        session.errors.clear()
        session.session_id = "sid"
        session.rpc = lambda method, params, timeout=15: (
            _facts(session) if method == "peri/host_policy" else {"servers": []}
        )
        with pytest.raises(NativeRejected, match="LIVE_CATALOG_MISMATCH"):
            session.verify_runtime()
        _stop(process)
        process = _sleep()
        session.process = process
        session.state = "blocked"
        session.errors.clear()
        session.session_id = None
        calls: list[str] = []

        def rpc(method, params, timeout=15):
            calls.append(method)
            if method == "peri/host_policy":
                seen = calls.count("peri/host_policy")
                return _facts(session, 1 if seen >= 3 else 0)
            if method == "session/new":
                return {"sessionId": "sid-1"}
            if method == "session/set_config_option":
                return {"configOptions": [{"id": "mode", "currentValue": "default"}]}
            if method == "mcp/list":
                return _catalog()
            if method == "session/prompt":
                session.approvals.append({"allowed": False, "name": "qual"})
                return {}
            return {}

        session.rpc = rpc
        session.stub = SimpleNamespace(
            requests=[], failures=[], script=lambda *args, **kwargs: None
        )
        session.errors.clear()
        session.approvals.clear()
        session.handler_calls = 0
        session.tool_arrivals = 0
        session.state = "blocked"
        session.initialize(timeout=1)
        assert session.state == "ready"
        assert session.approve is True
        assert session.session_id == "sid-1"
        session.state = "blocked"
        session.session_id = None
        session.errors.clear()
        session.rpc = lambda method, params, timeout=15: (
            {} if method == "initialize" else _facts(session)
        )
        session.rpc = lambda method, params, timeout=15: {}
        with pytest.raises(NativeRejected, match="HOST_POLICY_MISMATCH"):
            session.initialize(timeout=1)
        process.wait(timeout=5)
        assert process.poll() is not None
    finally:
        _stop(process)


def test_initialize_rejects_bad_mode_activity_and_catalog(monkeypatch) -> None:
    monkeypatch.setattr(
        "encode_pipeline.services.peri_runtime.session.stop_scope",
        lambda _unit: None,
    )
    monkeypatch.setattr(
        "encode_pipeline.services.peri_runtime.session.time.sleep",
        lambda _seconds: None,
    )

    def run(script: Callable[[str, NativeSession], object]) -> None:
        process = _sleep()
        session = _shell(process)
        try:
            session.stub = SimpleNamespace(
                requests=[], failures=[], script=lambda *args, **kwargs: None
            )

            def rpc(method, params, timeout=15):
                return script(method, session)

            session.rpc = rpc
            with pytest.raises(NativeRejected):
                session.initialize(timeout=0)
        finally:
            _stop(process)

    run(lambda method, session: {} if method != "peri/host_policy" else _facts(session))
    run(
        lambda method, session: (
            {"sessionId": 4}
            if method == "session/new"
            else _facts(session)
            if method == "peri/host_policy"
            else {}
        )
    )
    run(
        lambda method, session: (
            {"sessionId": "sid"}
            if method == "session/new"
            else {"configOptions": [{"id": "mode", "currentValue": "other"}]}
            if method == "session/set_config_option"
            else _facts(session)
            if method == "peri/host_policy"
            else {}
        )
    )

    process = _sleep()
    session = _shell(process)
    try:
        clocks = iter([0, 0, 1])
        monkeypatch.setattr(
            "encode_pipeline.services.peri_runtime.session.time.monotonic",
            lambda: next(clocks, 2),
        )
        session.session_id = "sid"
        session.rpc = lambda method, params, timeout=1: {"servers": [{}]}
        with pytest.raises(NativeRejected, match="DISCOVERY_TIMEOUT"):
            session._wait_for_catalog(0)
        session.rpc = lambda method, params, timeout=1: {
            "servers": [
                {
                    "name": "other",
                    "connectionStatus": "connected",
                    "tools": [],
                }
            ]
        }
        with pytest.raises(NativeRejected, match="TOOL_CATALOG_MISMATCH"):
            session._wait_for_catalog(1)
        assert session._egress_allowed() is False
        session.state = "blocked"
        session.qualification_egress = True
        assert session._egress_allowed() is True
        session.state = "ready"
        session.qualification_egress = False
        session.authorization_check = lambda: True
        assert session._egress_allowed() is True
        session.authorization_check = lambda: (_ for _ in ()).throw(RuntimeError("x"))
        assert session._egress_allowed() is False
        session.session_id = None
        with pytest.raises(NativeRejected, match="BINDING_REJECTED"):
            session.prompt("hello")
        session.session_id = "sid"
        session.state = "ready"
        with pytest.raises(NativeRejected, match="PROMPT_SHAPE_REJECTED"):
            session.prompt("/quit")
        session.active_turn = True
        with pytest.raises(NativeRejected, match="TURN_BUSY"):
            session.prompt("hello")
        session.active_turn = False
        session.stub = None
        with pytest.raises(NativeRejected, match="STUB_CLOSED"):
            session.run_authoring_turn("answer")
        session.stub = SimpleNamespace(script=lambda *args, **kwargs: None)
        session.prompt = lambda text, timeout=45: "explained"
        session.last_tool_ok = False
        with pytest.raises(NativeRejected, match="TOOL_REJECTED"):
            session.run_authoring_turn("answer")
        calls = {"n": 0}

        def scripted_prompt(text, timeout=45):
            calls["n"] += 1
            session.last_tool_ok = True
            session.last_tool_payload = {"ok": True, "n": calls["n"]}
            return "done"

        session.prompt = scripted_prompt
        assert session.run_authoring_turn("answer") == "done"
        assert session.tool_payloads[1]["n"] == 2
    finally:
        _stop(process)
