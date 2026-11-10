"""Host-protocol coverage that does not require the pinned Peri binary."""

from __future__ import annotations

import json
import os
from pathlib import Path
import socket
import struct
import threading
import time

import pytest

from encode_pipeline.services.peri_runtime import isolation, owner as owner_module
from encode_pipeline.services.peri_runtime import session as session_module
from encode_pipeline.services.peri_runtime.identity import (
    CANONICAL_DRAFT_TOOL,
    CANONICAL_SCHEMA_TOOL,
    SCHEMA_TOOL,
    SERVER_NAME,
)
from encode_pipeline.services.peri_runtime.limits import (
    cgroup_is_qualified,
    cgroup_path_is_qualified,
    cgroup_root_for_pid,
    qualified_native_pid,
    stop_scope,
)
from encode_pipeline.services.peri_runtime.owner import NativeAuthoringBridge
from encode_pipeline.services.peri_runtime.policy import (
    SYNTHETIC_API_KEY,
    process_policy,
)
from encode_pipeline.services.peri_runtime.relay import MAX_BYTES, Relay, forward
from encode_pipeline.services.peri_runtime.socket_server import (
    OwnerServer,
    RemoteAuthoringBridge,
    main as owner_main,
)
from encode_pipeline.services.peri_runtime.session import NativeSession
from encode_pipeline.services.peri_runtime.stub import ModelStub
from encode_pipeline.services.authoring_bridge import (
    BridgeDisconnected,
    BridgeSessionRequest,
    BridgeStatus,
    BridgeTurnRequest,
)


class PipeProcess:
    def __init__(self) -> None:
        child_in_r, child_in_w = os.pipe()
        child_out_r, child_out_w = os.pipe()
        self.stdin = os.fdopen(child_in_w, "wb", buffering=0)
        self.stdout = os.fdopen(child_out_r, "rb", buffering=0)
        self._input = os.fdopen(child_in_r, "rb", buffering=0)
        self._output = os.fdopen(child_out_w, "wb", buffering=0)
        self.pid = 2**22
        self.args = ["bwrap", "peri"]
        self.peri_native_pid = os.getpid()
        self.peri_scope_unit = "hw-peri-a1-unittest"
        self.session: NativeSession | None = None
        self._alive = True

    def poll(self) -> int | None:
        return None if self._alive else 0

    def wait(self, timeout: float | None = None) -> int:
        del timeout
        self._alive = False
        return 0

    def kill(self) -> None:
        self._alive = False


def _read_exact(handle, count: int) -> bytes:
    chunks = []
    remaining = count
    while remaining:
        piece = handle.read(remaining)
        if not piece:
            break
        chunks.append(piece)
        remaining -= len(piece)
    return b"".join(chunks)


def _serve(process: PipeProcess) -> None:
    while process.session is None:
        time.sleep(0.005)
    session = process.session
    header = _read_exact(process._input, 4)
    if len(header) == 4:
        _read_exact(process._input, struct.unpack(">I", header)[0])
    policy_calls = 0
    prompts = 0
    next_id = 1000
    while True:
        line = process._input.readline()
        if not line:
            return
        try:
            message = json.loads(line)
        except json.JSONDecodeError:
            continue
        method = message.get("method")
        request_id = message.get("id")
        if method == "peri/host_policy":
            policy_calls += 1
            result = {
                "version": 1,
                "prediction_enabled": False,
                "auto_compact_enabled": False,
                "tool_policy": process_policy(),
                "prediction_started": 0,
                "prediction_suppressed": policy_calls,
            }
        elif method == "mcp/list":
            result = {
                "servers": [
                    {
                        "name": SERVER_NAME,
                        "connectionStatus": "connected",
                        "tools": [
                            {"name": "get_authoring_schema"},
                            {"name": "get_draft_summary"},
                        ],
                    }
                ]
            }
        elif method == "session/new":
            result = {"sessionId": "native-session"}
        elif method == "session/set_config_option":
            result = {"configOptions": [{"id": "mode", "currentValue": "default"}]}
        elif method == "session/prompt":
            prompts += 1
            if session.state == "ready" and session.stub is not None:
                session.stub.responses.clear()
            wire = SCHEMA_TOOL if prompts == 1 else "get_draft_summary"
            canonical = (
                CANONICAL_SCHEMA_TOOL if wire == SCHEMA_TOOL else CANONICAL_DRAFT_TOOL
            )
            _emit(
                process,
                {
                    "jsonrpc": "2.0",
                    "id": next_id,
                    "method": "session/request_permission",
                    "params": {
                        "sessionId": session.session_id,
                        "toolCall": {"title": canonical, "rawInput": {}},
                        "options": [
                            {"kind": "allow_once", "optionId": "allow_once"},
                            {"kind": "reject_once", "optionId": "reject_once"},
                        ],
                    },
                },
            )
            next_id += 1
            process._input.readline()
            if session.state == "ready":
                _emit(
                    process,
                    {
                        "jsonrpc": "2.0",
                        "id": next_id,
                        "method": "mcp/message",
                        "params": {
                            "serverId": session.server_id,
                            "sessionId": session.session_id,
                            "method": "tools/call",
                            "requestId": f"req-{prompts}",
                            "params": {"name": wire, "arguments": {}},
                        },
                    },
                )
                next_id += 1
                process._input.readline()
            _emit(
                process,
                {
                    "jsonrpc": "2.0",
                    "method": "session/update",
                    "params": {
                        "sessionId": session.session_id,
                        "update": {
                            "sessionUpdate": "agent_message_chunk",
                            "content": {"type": "text", "text": "answer"},
                        },
                    },
                },
            )
            result = {}
        elif method == "initialize":
            result = {}
        else:
            result = {}
        if request_id is not None and "method" in message:
            _emit(
                process,
                {"jsonrpc": "2.0", "id": request_id, "result": result},
            )


def _emit(process: PipeProcess, payload: dict) -> None:
    process._output.write(json.dumps(payload).encode() + b"\n")
    process._output.flush()


def test_fake_native_qualifies_reads_and_rejects_tool_failure(
    tmp_path: Path, monkeypatch
) -> None:
    created: list[PipeProcess] = []

    def fake_spawn(*_args, **_kwargs) -> PipeProcess:
        process = PipeProcess()
        created.append(process)
        threading.Thread(target=_serve, args=(process,), daemon=True).start()
        return process

    monkeypatch.setattr(session_module, "spawn", fake_spawn)
    monkeypatch.setattr(owner_module, "binary_matches", lambda _path: True)
    monkeypatch.setattr(owner_module, "isolation_paths_exist", lambda: True)

    def build(handler):
        real = NativeSession

        def wrapped(*args, **kwargs):
            session = real(*args, **kwargs)
            created[-1].session = session
            return session

        return wrapped

    monkeypatch.setattr(owner_module, "NativeSession", build(None))
    repository = type(
        "Repository",
        (),
        {
            "tool_read_state": staticmethod(lambda **_kwargs: ("ok", None, None)),
            "connection_still_authorized": staticmethod(lambda **_kwargs: True),
            "login_is_live": staticmethod(lambda *_args: True),
        },
    )()
    bridge = NativeAuthoringBridge(
        binary=tmp_path / "peri",
        state_root=tmp_path / "state",
        repository=repository,
    )

    def handler(**_kwargs):
        return {
            "contract_version": "hw-a1/1",
            "ok": True,
            "data": {"fields": []},
            "issues": [],
        }

    bridge.bind_tool_handler(handler)
    request = BridgeSessionRequest(
        sid="aas_" + "ab" * 16,
        generation=1,
        owner_epoch=bridge.owner_epoch(),
        capability_token="token",
        schema_digest="a" * 64,
        projection_digest="b" * 64,
        auth_session_digest="c" * 64,
        user_id="usr_" + "d" * 32,
    )
    try:
        admission = bridge.admit(request)
        assert admission.state == "ready"
        ready = created[-1].session
        assert ready is not None and ready.stub is not None
        ready.stub.responses.clear()
        try:
            answer = bridge.execute_turn(
                BridgeTurnRequest(
                    sid=request.sid,
                    turn_id="atr_" + "cd" * 16,
                    generation=1,
                    instance_id=admission.instance_id,
                    owner_epoch=request.owner_epoch,
                    capability_token="token",
                    question={
                        "intent": "suggest",
                        "field_tokens": [],
                        "preferences": [],
                    },
                    issue_codes=(),
                    schema_digest="a" * 64,
                    projection_digest="b" * 64,
                    base_revision=0,
                )
            )
        except BridgeDisconnected:
            failed = created[-1].session
            raise AssertionError(
                {
                    "errors": None if failed is None else list(failed.errors),
                    "denied": None if failed is None else list(failed.denied_callbacks),
                    "calls": None if failed is None else failed.handler_calls,
                    "ok": None if failed is None else failed.last_tool_ok,
                    "failures": None
                    if failed is None or failed.stub is None
                    else list(failed.stub.failures),
                }
            ) from None
        assert answer == "answer"
        session = created[-1].session
        assert session is not None
        assert session.handler_calls == 2
        stub = session.stub
        assert stub is not None
        _post_stub(stub, scripted=True)
        session.authorization_check = lambda: False
        _post_stub(stub, scripted=False)
        bridge._handler = lambda **_kwargs: {
            "contract_version": "hw-a1/1",
            "ok": False,
            "data": None,
            "issues": [],
        }
        try:
            bridge.execute_turn(
                BridgeTurnRequest(
                    sid=request.sid,
                    turn_id="atr_" + "ef" * 16,
                    generation=1,
                    instance_id=admission.instance_id,
                    owner_epoch=request.owner_epoch,
                    capability_token="token",
                    question={
                        "intent": "explain_issues",
                        "field_tokens": [],
                        "preferences": [],
                    },
                    issue_codes=("ENCODE_OPTIONS_INVALID",),
                    schema_digest="a" * 64,
                    projection_digest="b" * 64,
                    base_revision=0,
                )
            )
        except BridgeDisconnected:
            pass
        else:
            raise AssertionError("failed tool returned a proposal")
        assert bridge.sweep() >= 0
        bridge.drop_authorizations(("c" * 64,))
    finally:
        for process in created:
            process._output.close()
            process._alive = False
        bridge.shutdown()


def _post_stub(stub: ModelStub, *, scripted: bool) -> None:
    if scripted:
        stub.script("hosted answer")
    client = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    client.connect(str(stub.path))
    body = json.dumps(
        {
            "tools": [
                {"function": {"name": "SearchExtraTools"}},
                {"function": {"name": "ExecuteExtraTool"}},
            ]
        }
    ).encode()
    header = (
        "POST /v1/chat/completions HTTP/1.1\r\n"
        f"Authorization: Bearer {SYNTHETIC_API_KEY}\r\n"
        f"Content-Length: {len(body)}\r\n"
        "\r\n"
    ).encode()
    client.sendall(header + body)
    client.recv(4096)
    client.close()


def test_limits_and_relay_and_spawn_guard(tmp_path: Path, monkeypatch) -> None:
    root = tmp_path / "cgroup"
    root.mkdir()
    (root / "memory.max").write_text("536870912\n", encoding="utf-8")
    (root / "memory.swap.max").write_text("0\n", encoding="utf-8")
    (root / "pids.max").write_text("128\n", encoding="utf-8")
    (root / "cpu.max").write_text("100000 100000\n", encoding="utf-8")
    assert cgroup_path_is_qualified(root) is True
    (root / "memory.max").write_text("max\n", encoding="utf-8")
    assert cgroup_path_is_qualified(root) is False
    assert cgroup_root_for_pid(os.getpid()) is not None
    assert cgroup_is_qualified() is False
    try:
        qualified_native_pid("hw-peri-a1-missing", os.getpid(), timeout=0.05)
    except OSError as error:
        assert "not qualified" in str(error)
    stop_scope(None)
    stop_scope("hw-peri-a1-missing")

    upstream_path = tmp_path / "upstream.sock"
    upstream = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    upstream.bind(str(upstream_path))
    upstream.listen(1)
    left, right = socket.socketpair()
    outcome: list[str] = []

    def run() -> None:
        outcome.append(forward(left, str(upstream_path), idle=0.2, total=1))

    worker = threading.Thread(target=run)
    worker.start()
    accepted, _address = upstream.accept()
    right.sendall(b"ping")
    assert accepted.recv(4) == b"ping"
    right.close()
    worker.join(2)
    assert outcome == ["disconnect"]
    accepted.close()
    upstream.close()

    listener = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    listener.bind(str(tmp_path / "relay.sock"))
    listener.listen(1)
    relay = Relay(listener, str(upstream_path))
    relay.close()

    argv = isolation.command(
        tmp_path / "peri",
        tmp_path / "state",
        tmp_path / "model.sock",
    )
    assert argv[0].endswith("bwrap")
    assert "/runtime/peri" in argv

    monkeypatch.delenv("DBUS_SESSION_BUS_ADDRESS", raising=False)
    monkeypatch.delenv("XDG_RUNTIME_DIR", raising=False)
    try:
        isolation.spawn(
            tmp_path / "peri",
            tmp_path / "state",
            tmp_path / "model.sock",
            None,
        )
    except OSError as error:
        assert "systemd user session" in str(error)
    else:
        raise AssertionError("spawn proceeded without a user session")


def _serve_line(path: Path, body: bytes) -> None:
    ready = threading.Event()

    def run() -> None:
        server = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        server.bind(str(path))
        server.listen(1)
        ready.set()
        connection, _address = server.accept()
        connection.recv(65536)
        if body:
            connection.sendall(body)
        connection.close()
        server.close()

    threading.Thread(target=run, daemon=True).start()
    assert ready.wait(2)


def test_owner_entry_client_errors_and_relay_budget(
    tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.delenv("HELIXWEAVE_PERI_BINARY", raising=False)
    monkeypatch.delenv("HELIXWEAVE_PERI_OWNER_SOCKET", raising=False)
    with pytest.raises(SystemExit) as missing:
        owner_main()
    assert missing.value.code == 2

    monkeypatch.setenv("HELIXWEAVE_PERI_BINARY", str(tmp_path / "not-peri"))
    monkeypatch.setenv("HELIXWEAVE_PERI_OWNER_SOCKET", str(tmp_path / "owner.sock"))
    monkeypatch.setenv(
        "ENCODE_PIPELINE_DATABASE_URL",
        f"sqlite:///{tmp_path / 'owner.db'}",
    )
    (tmp_path / "not-peri").write_bytes(b"not-the-fixed-runtime")
    with pytest.raises(SystemExit) as rejected:
        owner_main()
    assert rejected.value.code == 1

    remote = RemoteAuthoringBridge(tmp_path / "missing.sock")
    assert remote.status().available is False
    with pytest.raises(BridgeDisconnected):
        remote.owner_epoch()

    refused = tmp_path / "refused.sock"
    _serve_line(refused, b"not-json\n")
    with pytest.raises(BridgeDisconnected):
        RemoteAuthoringBridge(refused).owner_epoch()
    unavailable = tmp_path / "unavailable.sock"
    _serve_line(unavailable, b'{"ok": true, "available": false}\n')
    assert RemoteAuthoringBridge(unavailable).status().available is False
    denied = tmp_path / "denied.sock"
    _serve_line(denied, b'{"ok": false}\n')
    with pytest.raises(BridgeDisconnected):
        RemoteAuthoringBridge(denied).admit(
            BridgeSessionRequest(
                sid="aas_" + "ab" * 16,
                generation=1,
                owner_epoch="epoch",
                capability_token="token",
                schema_digest="a" * 64,
                projection_digest="b" * 64,
            )
        )
    empty = tmp_path / "empty.sock"
    _serve_line(empty, b'{"ok": true, "answer": ""}\n')
    with pytest.raises(BridgeDisconnected):
        RemoteAuthoringBridge(empty).execute_turn(
            BridgeTurnRequest(
                sid="aas_" + "ab" * 16,
                turn_id="atr_" + "cd" * 16,
                generation=1,
                instance_id="instance01",
                owner_epoch="epoch",
                capability_token="token",
                question={},
                issue_codes=(),
                schema_digest="a" * 64,
                projection_digest="b" * 64,
                base_revision=0,
            )
        )
    closed = tmp_path / "closed.sock"
    _serve_line(closed, b"")
    RemoteAuthoringBridge(closed).close(
        BridgeSessionRequest(
            sid="aas_" + "ab" * 16,
            generation=1,
            owner_epoch="epoch",
            capability_token="token",
            schema_digest="a" * 64,
            projection_digest="b" * 64,
        )
    )
    counted = tmp_path / "counted.sock"
    _serve_line(counted, b'{"ok": true}\n')
    assert RemoteAuthoringBridge(counted).drop_authorizations(("c" * 64,)) == 0

    class HandleBridge:
        def status(self) -> BridgeStatus:
            return BridgeStatus(available=False, reason_code="BRIDGE_UNAVAILABLE")

        def owner_epoch(self) -> str:
            raise BridgeDisconnected

        def admit(self, request: BridgeSessionRequest):
            del request
            raise BridgeDisconnected

    server = OwnerServer(HandleBridge(), object(), tmp_path / "handle.sock")
    assert server.handle({"method": "status"})["available"] is False
    assert server.handle({"method": "open", "owner_epoch": "old"})["error"] == (
        "unavailable"
    )

    class ReadyBridge(HandleBridge):
        def owner_epoch(self) -> str:
            return "current"

    ready = OwnerServer(ReadyBridge(), object(), tmp_path / "ready.sock")
    assert ready.handle({"method": "drop_auth"})["error"] == "invalid_request"
    assert ready.handle({"method": "sweep"}) == {"ok": True, "closed": 0}
    assert ready.handle({"method": "unknown"})["error"] == "unknown_method"
    assert ready.handle({"method": "turn", "owner_epoch": "other"})["error"] == (
        "epoch_rejected"
    )
    assert ready.handle({"method": "open", "owner_epoch": "current"})["error"] == (
        "disconnected"
    )

    budget_path = tmp_path / "budget.sock"
    budget = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    budget.bind(str(budget_path))
    budget.listen(1)
    left, right = socket.socketpair()
    outcome: list[str] = []

    def overrun() -> None:
        outcome.append(forward(left, str(budget_path), idle=1, total=3))

    worker = threading.Thread(target=overrun)
    worker.start()
    accepted, _address = budget.accept()

    def drain() -> None:
        while True:
            try:
                if not accepted.recv(65536):
                    return
            except OSError:
                return

    threading.Thread(target=drain, daemon=True).start()
    try:
        right.sendall(b"x" * (MAX_BYTES + 1))
    except OSError:
        pass
    worker.join(4)
    assert outcome == ["byte_budget"]
    right.close()
    accepted.close()
    budget.close()

    idle_path = tmp_path / "idle.sock"
    idle_upstream = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    idle_upstream.bind(str(idle_path))
    idle_upstream.listen(1)
    quiet_left, quiet_right = socket.socketpair()
    assert forward(quiet_left, str(idle_path), idle=0.05, total=0.2) == "deadline"
    quiet_right.close()
    idle_upstream.close()
