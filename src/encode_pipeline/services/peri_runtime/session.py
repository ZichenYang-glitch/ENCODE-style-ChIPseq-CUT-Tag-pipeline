"""One isolated Peri ACP session for the restricted authoring bridge."""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable
import hashlib
import json
import os
from pathlib import Path
import queue
import select
import shutil
import signal
import struct
import subprocess
import tempfile
import threading
import time
import uuid

from encode_pipeline.services.peri_runtime.identity import (
    CANONICAL_TOOLS,
    SERVER_NAME,
    WIRE_TOOLS,
)
from encode_pipeline.services.peri_runtime.isolation import spawn
from encode_pipeline.services.peri_runtime.limits import cgroup_root_for_pid, stop_scope
from encode_pipeline.services.peri_runtime.policy import TOOLS, process_policy, settings
from encode_pipeline.services.peri_runtime.stub import ModelStub


DROP = object()
ToolHandler = Callable[[str, dict[str, object], "NativeSession"], dict[str, object]]


class NativeRejected(RuntimeError):
    """Stable bridge failure without a private exception payload."""


def grant_key(name: str, arguments: dict[str, object]) -> tuple[str, str]:
    encoded = json.dumps(arguments, sort_keys=True).encode()
    return name, hashlib.sha256(encoded).hexdigest()


class NativeSession:
    """Fixed-binary Peri process with the two A1 read tools."""

    def __init__(
        self,
        binary: Path,
        directory: Path,
        *,
        service_generation: int,
        tool_handler: ToolHandler,
    ) -> None:
        self.binary = binary
        self.directory = directory
        self.service_generation = service_generation
        self.tool_handler = tool_handler
        self.policy = process_policy()
        self.instance = uuid.uuid4().hex
        self.server_id = f"{SERVER_NAME}-{self.instance}"
        self.session_id: str | None = None
        self.state = "blocked"
        self.approve = False
        self.handler_calls = 0
        self.tool_arrivals = 0
        self.approvals: list[dict[str, object]] = []
        self.denied_callbacks: list[str] = []
        self.errors: list[str] = []
        self.grants: Counter[tuple[str, str]] = Counter()
        self.used_requests: set[str] = set()
        self.assistant_parts: list[str] = []
        self.events: list[dict[str, object]] = []
        self.event_bytes = 0
        self.active_turn = False
        self.pending: dict[int, queue.Queue[dict[str, object]]] = {}
        self.request_id = 0
        self.gate = threading.RLock()
        self.write_lock = threading.Lock()
        self.pending_lock = threading.Lock()
        self.process: subprocess.Popen[bytes] | None = None
        self.reader: threading.Thread | None = None
        self.stub: ModelStub | None = None
        self.log = None
        self.argv: list[str] | None = None
        self.closed_resources = False
        self.turn_context: object | None = None
        self.authorization_check: Callable[[], bool] | None = None
        self.qualification_egress = False
        self.last_tool_ok = False
        self.last_tool_payload: dict[str, object] | None = None
        self.tool_payloads: tuple[dict[str, object], ...] = ()
        self.auth_session_digest = ""
        self.user_id = ""
        self.bound_sid = ""
        self.native_pid: int | None = None
        self.leader_pid: int | None = None
        self.scope_unit: str | None = None
        self.scope_cgroup: str | None = None
        directory.mkdir(parents=True, exist_ok=False)
        self.socket_dir = Path(tempfile.mkdtemp(prefix="hw-peri-sock-"))
        try:
            self.log = (directory / "stderr.log").open("wb")
            self.stub = ModelStub(self.socket_dir / "model.sock")
            self.stub.on_catalog_failure = self.invalidate
            self.stub.authorize = self._egress_allowed
            self.process = spawn(
                binary,
                directory / "state",
                self.stub.path,
                self.log,
                policy=self.policy,
                prediction="off",
            )
            self.native_pid = getattr(self.process, "peri_native_pid", None)
            self.leader_pid = self.process.pid
            self.scope_unit = getattr(self.process, "peri_scope_unit", None)
            root = (
                None
                if self.native_pid is None
                else cgroup_root_for_pid(self.native_pid)
            )
            self.scope_cgroup = None if root is None else str(root)
            self.argv = list(self.process.args)
            payload = json.dumps(settings()).encode()
            os.set_blocking(self.process.stdin.fileno(), False)
            self.bootstrap_payload = struct.pack(">I", len(payload)) + payload
            self.reader = threading.Thread(target=self._read, daemon=True)
            self.reader.start()
        except BaseException:
            self.close()
            raise

    def invalidate(self, code: str) -> None:
        with self.gate:
            self.state = "closed"
            self.qualification_egress = False
            self.approve = False
            self.grants.clear()
            self.active_turn = False
            self.errors.append(code)
        if self.process is not None and self.process.poll() is None:
            try:
                os.killpg(self.process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass

    def _write_bytes(self, wire: bytes) -> None:
        if self.process is None or self.process.stdin is None:
            raise NativeRejected("ACP_STDIN_CLOSED")
        if len(wire) > 262144:
            raise NativeRejected("ACP_WRITE_BUDGET")
        deadline = time.monotonic() + 1
        with self.write_lock:
            while wire:
                remaining = deadline - time.monotonic()
                if (
                    remaining <= 0
                    or not select.select([], [self.process.stdin], [], remaining)[1]
                ):
                    raise NativeRejected("ACP_WRITE_TIMEOUT")
                try:
                    count = os.write(self.process.stdin.fileno(), wire)
                except BlockingIOError:
                    continue
                wire = wire[count:]

    def _send(self, value: dict[str, object]) -> None:
        self._write_bytes(json.dumps(value, separators=(",", ":")).encode() + b"\n")

    def consume_line(self, line: bytes) -> None:
        if len(line) > 1024 * 1024:
            raise NativeRejected("ACP_FRAME_TOO_LARGE")
        try:
            message = json.loads(line)
        except (ValueError, UnicodeError):
            raise NativeRejected("ACP_MALFORMED") from None
        if not isinstance(message, dict) or message.get("jsonrpc") != "2.0":
            raise NativeRejected("ACP_MALFORMED")
        if len(self.events) >= 2000 or self.event_bytes + len(line) > 8388608:
            raise NativeRejected("ACP_EVENT_BUDGET")
        self.event_bytes += len(line)
        self.events.append(message)
        self._record_assistant_text(message)
        if "method" in message and "id" in message:
            result = self._callback(str(message["method"]), message.get("params", {}))
            if result is not DROP:
                self._send({"jsonrpc": "2.0", "id": message["id"], **result})
        elif "id" in message:
            with self.pending_lock:
                target = self.pending.get(message["id"])
            if target is not None:
                target.put(message)

    def _record_assistant_text(self, message: dict[str, object]) -> None:
        if message.get("method") != "session/update":
            return
        params = message.get("params")
        if not isinstance(params, dict):
            return
        update = params.get("update")
        if not isinstance(update, dict):
            return
        if update.get("sessionUpdate") not in {"agent_message_chunk", "agent_message"}:
            return
        content = update.get("content")
        if isinstance(content, dict) and content.get("type") == "text":
            text = content.get("text")
            if isinstance(text, str):
                self.assistant_parts.append(text)
        elif isinstance(content, str):
            self.assistant_parts.append(content)

    def _read(self) -> None:
        try:
            assert self.process is not None and self.process.stdout is not None
            while line := self.process.stdout.readline(1024 * 1024 + 1):
                self.consume_line(line)
            if self.state != "closed":
                self.invalidate("ACP_EOF")
        except BaseException as error:
            self.invalidate(
                str(error)
                if isinstance(error, NativeRejected)
                else "ACP_TRANSPORT_FAILURE"
            )
        finally:
            with self.pending_lock:
                for target in self.pending.values():
                    target.put({"error": {"message": "transport closed"}})

    def rpc(
        self, method: str, params: dict[str, object], *, timeout: float = 15
    ) -> dict[str, object]:
        if self.state == "closed":
            raise NativeRejected("CLOSED")
        with self.pending_lock:
            self.request_id += 1
            request_id = self.request_id
            reply: queue.Queue[dict[str, object]] = queue.Queue()
            self.pending[request_id] = reply
        try:
            self._send(
                {"jsonrpc": "2.0", "id": request_id, "method": method, "params": params}
            )
            result = reply.get(timeout=timeout)
            if "error" in result:
                raise NativeRejected("ACP_REMOTE_REJECTION")
            payload = result.get("result", {})
            if not isinstance(payload, dict):
                raise NativeRejected("ACP_REMOTE_REJECTION")
            return payload
        except queue.Empty:
            self.invalidate("RPC_TIMEOUT")
            raise NativeRejected("RPC_TIMEOUT") from None
        except NativeRejected:
            self.invalidate("RPC_FAILURE")
            raise
        except (OSError, KeyError, TypeError):
            self.invalidate("RPC_TRANSPORT_FAILURE")
            raise NativeRejected("RPC_TRANSPORT_FAILURE") from None
        finally:
            with self.pending_lock:
                self.pending.pop(request_id, None)

    def _callback(self, method: str, params: object) -> dict[str, object] | object:
        with self.gate:
            if self.state == "closed" or not isinstance(params, dict):
                self.denied_callbacks.append("binding")
                return {"error": {"code": -32602, "message": "binding rejected"}}
            sid = params.get("sessionId", params.get("session_id"))
            if (
                sid is not None
                and self.session_id is not None
                and sid != self.session_id
            ):
                self.denied_callbacks.append("session")
                return {"error": {"code": -32602, "message": "session rejected"}}
            if method == "session/request_permission":
                return self._permission(params)
            if method != "mcp/message" or params.get("serverId") != self.server_id:
                self.denied_callbacks.append("method_or_server")
                return {"error": {"code": -32601, "message": "method not exposed"}}
            return self._mcp(params)

    def _permission(self, params: dict[str, object]) -> dict[str, object] | object:
        tool = params.get("toolCall")
        tool = tool if isinstance(tool, dict) else {}
        name = tool.get("title")
        arguments = tool.get("rawInput")
        if not isinstance(arguments, dict):
            arguments = {}
        eligible = (
            self.active_turn
            and self.state == "ready"
            and self.approve
            and name in CANONICAL_TOOLS
            and isinstance(arguments, dict)
        )
        self.approvals.append({"allowed": eligible, "name": name})
        expected = "allow_once" if eligible else "reject_once"
        options = params.get("options")
        option = None
        if isinstance(options, list):
            for item in options:
                if not isinstance(item, dict):
                    continue
                if item.get("kind") == expected or item.get("optionId") == expected:
                    option = item.get("optionId")
                    break
        if not isinstance(option, str):
            return {"result": {"outcome": {"outcome": "cancelled"}}}
        if eligible:
            self.grants[grant_key(str(name), arguments)] += 1
        return {"result": {"outcome": {"outcome": "selected", "optionId": option}}}

    def _mcp(self, params: dict[str, object]) -> dict[str, object] | object:
        inner = params.get("method")
        data = params.get("params", {})
        if not isinstance(data, dict):
            data = {}
        if inner == "server/discover":
            result = {
                "resultType": "complete",
                "supportedVersions": ["2026-07-28"],
                "capabilities": {"tools": {}},
                "ttlMs": 0,
                "cacheScope": "private",
                "_meta": {
                    "io.modelcontextprotocol/serverInfo": {
                        "name": SERVER_NAME,
                        "version": "1",
                    }
                },
            }
        elif inner == "initialize":
            result = {
                "protocolVersion": "2025-03-26",
                "capabilities": {"tools": {}},
                "serverInfo": {"name": SERVER_NAME, "version": "1"},
            }
        elif inner == "tools/list":
            result = {"tools": list(TOOLS)}
        elif inner == "ping":
            result = {}
        elif inner == "tools/call":
            return self._call_tool(params, data)
        else:
            self.denied_callbacks.append("inner_method")
            return {
                "result": {"error": {"code": -32601, "message": "method not exposed"}}
            }
        return {"result": {"result": result}}

    def _call_tool(
        self, params: dict[str, object], data: dict[str, object]
    ) -> dict[str, object]:
        self.tool_arrivals += 1
        name = data.get("name")
        arguments = data.get("arguments", {})
        if not isinstance(arguments, dict):
            arguments = {}
        request_id = params.get("requestId")
        canonical = f"mcp__{SERVER_NAME}__{name}"
        key = grant_key(canonical, arguments)
        authorized = (
            self.state == "ready"
            and self.active_turn
            and name in WIRE_TOOLS
            and arguments == {}
            and not (set(data) - {"name", "arguments", "_meta"})
            and isinstance(request_id, str)
            and request_id not in self.used_requests
            and self.grants[key] >= 1
        )
        if not authorized:
            self.denied_callbacks.append("tool_not_authorized")
            return {
                "result": {
                    "error": {"code": -32602, "message": "tool request rejected"}
                }
            }
        self.grants[key] -= 1
        self.used_requests.add(request_id)
        self.handler_calls += 1
        try:
            payload = self.tool_handler(str(name), arguments, self)
        except Exception:
            self.last_tool_ok = False
            self.last_tool_payload = None
            self.denied_callbacks.append("tool_handler_failed")
            return {
                "result": {
                    "error": {"code": -32603, "message": "tool request rejected"}
                }
            }
        self.last_tool_payload = payload if isinstance(payload, dict) else None
        self.last_tool_ok = isinstance(payload, dict) and payload.get("ok") is True
        return {
            "result": {
                "result": {
                    "content": [
                        {"type": "text", "text": json.dumps(payload, sort_keys=True)}
                    ],
                    "isError": False,
                }
            }
        }

    def verify_runtime(self, *, minimum_suppressed: int = 0) -> dict[str, object]:
        deadline = time.monotonic() + 2
        while True:
            facts = self.rpc("peri/host_policy", {})
            if (
                facts.get("version") != 1
                or facts.get("prediction_enabled") is not False
                or facts.get("auto_compact_enabled") is not False
                or facts.get("tool_policy") != self.policy
                or facts.get("prediction_started") != 0
            ):
                self.invalidate("HOST_POLICY_MISMATCH")
                raise NativeRejected("HOST_POLICY_MISMATCH")
            suppressed = facts.get("prediction_suppressed", -1)
            if isinstance(suppressed, int) and suppressed >= minimum_suppressed:
                break
            if time.monotonic() >= deadline:
                self.invalidate("TURN_END_POLICY_NOT_OBSERVED")
                raise NativeRejected("TURN_END_POLICY_NOT_OBSERVED")
            time.sleep(0.01)
        if self.session_id is None:
            return facts
        catalog = self.rpc("mcp/list", {"sessionId": self.session_id})
        servers = catalog.get("servers", [])
        if (
            not isinstance(servers, list)
            or len(servers) != 1
            or not isinstance(servers[0], dict)
            or servers[0].get("name") != SERVER_NAME
            or servers[0].get("connectionStatus") != "connected"
            or {
                item.get("name")
                for item in servers[0].get("tools", [])
                if isinstance(item, dict)
            }
            != WIRE_TOOLS
        ):
            self.invalidate("LIVE_CATALOG_MISMATCH")
            raise NativeRejected("LIVE_CATALOG_MISMATCH")
        return facts

    def initialize(self, *, timeout: float = 20) -> None:
        try:
            self._write_bytes(self.bootstrap_payload)
            self.rpc(
                "initialize",
                {
                    "protocolVersion": 1,
                    "clientCapabilities": {"_meta": {"peri.mcpControlV1": True}},
                },
                timeout=timeout,
            )
            self.verify_runtime()
            created = self.rpc(
                "session/new",
                {
                    "cwd": "/work",
                    "mcpServers": [
                        {"type": "acp", "name": SERVER_NAME, "serverId": self.server_id}
                    ],
                },
                timeout=timeout,
            )
            session_id = created.get("sessionId")
            if not isinstance(session_id, str):
                raise NativeRejected("SESSION_ID_REJECTED")
            self.session_id = session_id
            self.rpc(
                "session/set_mode",
                {"sessionId": self.session_id, "modeId": "default"},
                timeout=timeout,
            )
            modes = self.rpc(
                "session/set_config_option",
                {"sessionId": self.session_id, "configId": "mode", "value": "default"},
                timeout=timeout,
            )
            options = modes.get("configOptions")
            current = None
            if isinstance(options, list):
                for item in options:
                    if isinstance(item, dict) and item.get("id") == "mode":
                        current = item.get("currentValue")
            if current != "default":
                raise NativeRejected("MODE_VERIFICATION_FAILED")
            self._wait_for_catalog(timeout)
            assert self.stub is not None
            if self.stub.requests or self.handler_calls or self.tool_arrivals:
                raise NativeRejected("PRE_READY_ACTIVITY")
            before = self.verify_runtime()
            self.qualification_egress = True
            try:
                self.stub.script(
                    (
                        "SearchExtraTools",
                        {"query": "select:mcp__helixweave_a1__get_authoring_schema"},
                    ),
                    (
                        "ExecuteExtraTool",
                        {
                            "tool_name": "mcp__helixweave_a1__get_authoring_schema",
                            "params": {},
                        },
                    ),
                    "synthetic qualification complete",
                )
                self.rpc(
                    "session/prompt",
                    {
                        "sessionId": self.session_id,
                        "prompt": [
                            {"type": "text", "text": "A1 synthetic qualification only"}
                        ],
                    },
                    timeout=timeout,
                )
            finally:
                self.qualification_egress = False
            if (
                self.stub.failures
                or not self.approvals
                or any(item["allowed"] for item in self.approvals)
                or self.handler_calls
                or self.tool_arrivals
            ):
                raise NativeRejected("ACTUAL_DENIAL_NOT_VERIFIED")
            self.verify_runtime(
                minimum_suppressed=int(before.get("prediction_suppressed", 0)) + 1
            )
            with self.gate:
                if (
                    self.state != "blocked"
                    or self.process is None
                    or self.process.poll() is not None
                    or self.errors
                ):
                    raise NativeRejected("QUALIFICATION_CLOSED")
                self.state = "ready"
                self.approve = True
        except BaseException as error:
            self.invalidate("INITIALIZATION_FAILED")
            if isinstance(error, NativeRejected):
                raise
            raise NativeRejected("INITIALIZATION_TRANSPORT_FAILURE") from error

    def _wait_for_catalog(self, timeout: float) -> None:
        assert self.session_id is not None
        deadline = time.monotonic() + timeout
        while True:
            catalog = self.rpc(
                "mcp/list", {"sessionId": self.session_id}, timeout=timeout
            )
            servers = catalog.get("servers", [])
            if (
                isinstance(servers, list)
                and len(servers) == 1
                and isinstance(servers[0], dict)
            ):
                if servers[0].get("connectionStatus") == "connected":
                    names = {
                        item.get("name")
                        for item in servers[0].get("tools", [])
                        if isinstance(item, dict)
                    }
                    if servers[0].get("name") != SERVER_NAME or names != WIRE_TOOLS:
                        raise NativeRejected("TOOL_CATALOG_MISMATCH")
                    return
            if time.monotonic() >= deadline:
                raise NativeRejected("DISCOVERY_TIMEOUT")
            time.sleep(0.05)

    def _egress_allowed(self) -> bool:
        with self.gate:
            if self.state == "blocked" and self.qualification_egress:
                return True
            if self.state != "ready":
                return False
            check = self.authorization_check
        if check is None:
            return False
        try:
            return bool(check())
        except Exception:
            return False

    def run_authoring_turn(self, final_text: str) -> str:
        """Read both tools. A failed tool does not return the preset answer."""

        if self.stub is None:
            raise NativeRejected("STUB_CLOSED")
        self.last_tool_ok = False
        self._script_tool("get_authoring_schema", "schema read")
        self.prompt("Read the approved authoring schema.")
        schema_payload = self.last_tool_payload
        if not self.last_tool_ok or schema_payload is None:
            raise NativeRejected("TOOL_REJECTED")
        self.last_tool_ok = False
        self._script_tool("get_draft_summary", final_text)
        text = self.prompt("Read the approved draft projection.")
        draft_payload = self.last_tool_payload
        if not self.last_tool_ok or draft_payload is None:
            raise NativeRejected("TOOL_REJECTED")
        self.tool_payloads = (schema_payload, draft_payload)
        return text

    def _script_tool(self, wire_name: str, closing_text: str) -> None:
        assert self.stub is not None
        canonical = f"mcp__{SERVER_NAME}__{wire_name}"
        self.stub.script(
            ("SearchExtraTools", {"query": "select:" + canonical}),
            ("ExecuteExtraTool", {"tool_name": canonical, "params": {}}),
            closing_text,
        )

    def prompt(self, text: str, *, timeout: float = 45) -> str:
        if self.session_id is None or self.state != "ready":
            raise NativeRejected("BINDING_REJECTED")
        if text.lstrip().startswith("/") or len(text) > 4096:
            raise NativeRejected("PROMPT_SHAPE_REJECTED")
        with self.gate:
            if self.active_turn:
                raise NativeRejected("TURN_BUSY")
            self.active_turn = True
            self.grants.clear()
            self.assistant_parts.clear()
        try:
            before = self.verify_runtime()
            self.rpc(
                "session/prompt",
                {
                    "sessionId": self.session_id,
                    "prompt": [{"type": "text", "text": text}],
                },
                timeout=timeout,
            )
            self.verify_runtime(
                minimum_suppressed=int(before.get("prediction_suppressed", 0)) + 1
            )
            return "".join(self.assistant_parts).strip()
        finally:
            with self.gate:
                self.active_turn = False
                self.grants.clear()

    def close(self) -> None:
        if self.closed_resources:
            return
        stop_scope(self.scope_unit)
        self.invalidate("CLOSED_BY_HOST")
        if self.process is not None:
            try:
                self.process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                pass
            for stream in (self.process.stdin, self.process.stdout):
                if stream is None:
                    continue
                try:
                    stream.close()
                except (BrokenPipeError, OSError):
                    pass
        if self.reader is not None and self.reader is not threading.current_thread():
            self.reader.join(timeout=5)
        if self.stub is not None:
            self.stub.close()
        if self.log is not None:
            self.log.close()
        socket_dir = getattr(self, "socket_dir", None)
        if socket_dir is not None:
            shutil.rmtree(socket_dir, ignore_errors=True)
        self.closed_resources = True

    @property
    def pid(self) -> int | None:
        if self.native_pid is not None:
            return self.native_pid
        if self.process is None:
            return None
        return self.process.pid
