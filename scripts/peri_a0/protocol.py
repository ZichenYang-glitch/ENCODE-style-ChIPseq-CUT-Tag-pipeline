"""Narrow ACP/MCP experiment client; no platform integration or business consent."""

from collections import Counter
import hashlib
import json
import os
import queue
import signal
import struct
import subprocess
import threading
import time
import uuid

from configuration import TOOLS, settings
from sandbox import spawn
from stub import ModelStub
from observation import Observation

DROP = object()


class Rejected(RuntimeError):
    """Stable prototype failure code, without private exception payloads."""


def grant_key(name, arguments):
    return name, hashlib.sha256(
        json.dumps(arguments, sort_keys=True).encode()
    ).hexdigest()


class Client:
    def __init__(self, binary, directory, *, generation=1, tools=None):
        self.directory = directory
        directory.mkdir(parents=True, exist_ok=False)
        self.instance = uuid.uuid4().hex
        self.generation = generation
        self.server_id = "fixture-" + self.instance
        self.session = None
        self.state = "blocked"
        self.approve = False
        self.approval_failure = None
        self.callback_failure = None
        self.counter = self.lookups = self.tool_arrivals = 0
        self.approvals, self.denied_callbacks, self.events, self.errors = [], [], [], []
        self.declared_tools = TOOLS if tools is None else tools
        self.grants = Counter()
        self.used_requests = set()
        self.active_turn = False
        self.pending = {}
        self.request_id = 0
        self.gate = threading.RLock()
        self.write_lock = threading.Lock()
        self.pending_lock = threading.Lock()
        self.stub = self.log = None
        self.closed_resources = False
        self.process = self.reader = self.observer = None
        self.resources = None
        self.argv = None
        self.started_at = time.monotonic()
        self.ready_seconds = None
        self.closed_seconds = None
        self.cleanup_errors = []
        try:
            self.log = (directory / "stderr.log").open("wb")
            self.stub = ModelStub(directory / "model.sock")
            self.stub.on_catalog_failure = self.invalidate
            self.process = spawn(binary, directory / "state", self.stub.path, self.log)
            self.argv = self.process.args
            self.observer = Observation(self.process.pid)
            payload = json.dumps(settings()).encode()
            self.process.stdin.write(struct.pack(">I", len(payload)) + payload)
            self.process.stdin.flush()
            self.reader = threading.Thread(target=self._read, daemon=True)
            self.reader.start()
        except BaseException:
            self.close()
            raise

    def invalidate(self, code):
        with self.gate:
            if self.state != "closed" and hasattr(self, "started_at"):
                self.closed_seconds = time.monotonic() - self.started_at
            self.state = "closed"
            self.approve = False
            self.grants.clear()
            self.active_turn = False
            self.errors.append(code)
        if self.process is not None and self.process.poll() is None:
            try:
                os.killpg(self.process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass

    def _send(self, value):
        wire = json.dumps(value, separators=(",", ":")).encode() + b"\n"
        with self.write_lock:
            self.process.stdin.write(wire)
            self.process.stdin.flush()

    def consume_line(self, line):
        if len(line) > 1024 * 1024:
            raise Rejected("ACP_FRAME_TOO_LARGE")
        try:
            message = json.loads(line)
        except (ValueError, UnicodeError):
            raise Rejected("ACP_MALFORMED") from None
        if not isinstance(message, dict) or message.get("jsonrpc") != "2.0":
            raise Rejected("ACP_MALFORMED")
        if len(self.events) >= 2000:
            raise Rejected("ACP_EVENT_BUDGET")
        self.events.append(message)
        if "method" in message and "id" in message:
            result = self.callback(
                message["method"],
                message.get("params", {}),
                self.instance,
                self.generation,
            )
            if result is not DROP:
                self._send({"jsonrpc": "2.0", "id": message["id"], **result})
        elif "id" in message:
            with self.pending_lock:
                target = self.pending.get(message["id"])
            if target is not None:
                target.put(message)

    def _read(self):
        try:
            while line := self.process.stdout.readline(1024 * 1024 + 1):
                self.consume_line(line)
            if self.state != "closed":
                self.invalidate("ACP_EOF")
        except BaseException as error:
            self.invalidate(
                str(error) if isinstance(error, Rejected) else "ACP_TRANSPORT_FAILURE"
            )
        finally:
            with self.pending_lock:
                for target in self.pending.values():
                    target.put({"error": {"message": "transport closed"}})

    def rpc(self, method, params, *, timeout=15):
        if self.state == "closed":
            raise Rejected("CLOSED")
        with self.pending_lock:
            self.request_id += 1
            request_id = self.request_id
            reply = queue.Queue()
            self.pending[request_id] = reply
        try:
            self._send(
                {"jsonrpc": "2.0", "id": request_id, "method": method, "params": params}
            )
            result = reply.get(timeout=timeout)
            if "error" in result:
                raise Rejected("ACP_REMOTE_REJECTION")
            return result["result"]
        except queue.Empty:
            self.invalidate("RPC_TIMEOUT")
            raise Rejected("RPC_TIMEOUT") from None
        except Rejected:
            self.invalidate("RPC_FAILURE")
            raise
        except (OSError, KeyError, TypeError):
            self.invalidate("RPC_TRANSPORT_FAILURE")
            raise Rejected("RPC_TRANSPORT_FAILURE") from None
        finally:
            with self.pending_lock:
                self.pending.pop(request_id, None)

    def authorize_target(self, instance, generation, session):
        if (
            self.state != "ready"
            or instance != self.instance
            or generation != self.generation
            or session != self.session
            or session is None
        ):
            raise Rejected("BINDING_REJECTED")

    def invoke(self, method, params, *, instance=None, generation=None, timeout=45):
        with self.gate:
            self.authorize_target(
                instance or self.instance,
                self.generation if generation is None else generation,
                params.get("sessionId"),
            )
            if method != "session/prompt":
                raise Rejected("METHOD_NOT_EXPOSED")
            if set(params) != {"sessionId", "prompt"}:
                raise Rejected("UNKNOWN_PARAMETER")
            prompt = params["prompt"]
            if (
                not isinstance(prompt, list)
                or len(prompt) != 1
                or not isinstance(prompt[0], dict)
                or set(prompt[0]) != {"type", "text"}
                or prompt[0]["type"] != "text"
                or not isinstance(prompt[0]["text"], str)
                or len(prompt[0]["text"]) > 4096
            ):
                raise Rejected("PROMPT_SHAPE_REJECTED")
            if self.active_turn:
                raise Rejected("TURN_BUSY")
            self.active_turn = True
            self.grants.clear()
        try:
            return self.rpc(method, params, timeout=timeout)
        finally:
            with self.gate:
                self.active_turn = False
                self.grants.clear()

    def callback(self, method, params, instance, generation):
        with self.gate:
            return self._callback(method, params, instance, generation)

    def _callback(self, method, params, instance, generation):
        if (
            instance != self.instance
            or generation != self.generation
            or self.state == "closed"
            or not isinstance(params, dict)
        ):
            self.denied_callbacks.append("binding")
            return {"error": {"code": -32602, "message": "binding rejected"}}
        sid = params.get("sessionId", params.get("session_id"))
        if sid is not None and self.session is not None and sid != self.session:
            self.denied_callbacks.append("session")
            return {"error": {"code": -32602, "message": "session rejected"}}
        if method == "session/request_permission":
            tool = params.get("toolCall", {})
            name, arguments = tool.get("title"), tool.get("rawInput")
            eligible = (
                self.active_turn
                and sid == self.session
                and name in {"mcp__fixture__lookup", "mcp__fixture__increment"}
                and isinstance(arguments, dict)
            )
            allow = self.state == "ready" and self.approve and eligible
            self.approvals.append({"allowed": allow, "params": params})
            if self.approval_failure == "disconnect":
                self.invalidate("APPROVAL_DISCONNECTED")
                return DROP
            if self.approval_failure == "timeout":
                return DROP
            expected = "allow_once" if allow else "reject_once"
            option = next(
                (
                    item["optionId"]
                    for item in params.get("options", [])
                    if item.get("kind") == expected or item.get("optionId") == expected
                ),
                None,
            )
            if option is None:
                return {"result": {"outcome": {"outcome": "cancelled"}}}
            if allow:
                self.grants[grant_key(name, arguments)] += 1
            return {"result": {"outcome": {"outcome": "selected", "optionId": option}}}
        if method != "mcp/message" or params.get("serverId") != self.server_id:
            self.denied_callbacks.append("method_or_server")
            return {"error": {"code": -32601, "message": "method not exposed"}}
        inner, data = params.get("method"), params.get("params", {})
        if inner == "server/discover":
            result = {
                "resultType": "complete",
                "supportedVersions": ["2026-07-28"],
                "capabilities": {"tools": {}},
                "ttlMs": 0,
                "cacheScope": "private",
                "_meta": {
                    "io.modelcontextprotocol/serverInfo": {
                        "name": "a0-fixture",
                        "version": "1",
                    }
                },
            }
        elif inner == "initialize":
            result = {
                "protocolVersion": "2025-03-26",
                "capabilities": {"tools": {}},
                "serverInfo": {"name": "a0-fixture", "version": "1"},
            }
        elif inner == "tools/list":
            result = {"tools": self.declared_tools}
        elif inner == "ping":
            result = {}
        elif inner == "tools/call":
            self.tool_arrivals += 1
            if self.callback_failure:
                self.invalidate("TOOL_CALLBACK_DISCONNECTED")
                return DROP
            name, arguments = data.get("name"), data.get("arguments", {})
            valid = (
                name == "lookup"
                and arguments == {"object_id": "item-1"}
                or name == "increment"
                and arguments == {}
            )
            key = grant_key("mcp__fixture__" + str(name), arguments)
            request_id = params.get("requestId")
            if (
                self.state != "ready"
                or not self.active_turn
                or not valid
                or set(data) - {"name", "arguments", "_meta"}
                or not request_id
                or request_id in self.used_requests
                or self.grants[key] < 1
            ):
                self.denied_callbacks.append("tool_not_authorized")
                return {
                    "result": {
                        "error": {"code": -32602, "message": "tool request rejected"}
                    }
                }
            self.grants[key] -= 1
            self.used_requests.add(request_id)
            if name == "lookup":
                self.lookups += 1
                value = {"object_id": "item-1", "value": "synthetic"}
            else:
                self.counter += 1
                value = {"count": self.counter}
            result = {
                "content": [{"type": "text", "text": json.dumps(value)}],
                "isError": False,
            }
        else:
            self.denied_callbacks.append("inner_method")
            return {
                "result": {"error": {"code": -32601, "message": "method not exposed"}}
            }
        return {"result": {"result": result}}

    def initialize(self, *, mode="default", timeout=15):
        try:
            if mode != "default":
                raise Rejected("MODE_NOT_ALLOWED")
            self.rpc(
                "initialize",
                {
                    "protocolVersion": 1,
                    "clientCapabilities": {"_meta": {"peri.mcpControlV1": True}},
                },
                timeout=timeout,
            )
            result = self.rpc(
                "session/new",
                {
                    "cwd": "/work",
                    "mcpServers": [
                        {"type": "acp", "name": "fixture", "serverId": self.server_id}
                    ],
                },
                timeout=timeout,
            )
            self.session = result["sessionId"]
            self.rpc(
                "session/set_mode",
                {"sessionId": self.session, "modeId": mode},
                timeout=timeout,
            )
            modes = self.rpc(
                "session/set_config_option",
                {"sessionId": self.session, "configId": "mode", "value": mode},
                timeout=timeout,
            )
            option = next(
                item for item in modes["configOptions"] if item["id"] == "mode"
            )
            if option["currentValue"] != "default":
                raise Rejected("MODE_VERIFICATION_FAILED")
            deadline = time.monotonic() + timeout
            while True:
                catalog = self.rpc(
                    "mcp/list", {"sessionId": self.session}, timeout=timeout
                )
                servers = catalog.get("servers", [])
                if (
                    len(servers) == 1
                    and servers[0].get("connectionStatus") == "connected"
                ):
                    server = servers[0]
                    names = {item["name"] for item in server.get("tools", [])}
                    if (
                        self.declared_tools != TOOLS
                        or server.get("name") != "fixture"
                        or server.get("toolsCount") != len(TOOLS)
                        or names != {item["name"] for item in TOOLS}
                    ):
                        raise Rejected("TOOL_CATALOG_MISMATCH")
                    self.catalog = catalog
                    break
                if time.monotonic() >= deadline:
                    raise Rejected("DISCOVERY_TIMEOUT")
                time.sleep(0.05)
            if self.stub.requests or self.counter or self.lookups or self.tool_arrivals:
                raise Rejected("PRE_READY_ACTIVITY")
            # A fixed synthetic qualification prompt is internal test traffic,
            # never user input. Actual provider catalog and denial are checked
            # while blocked; a mismatch kills this instance before a response.
            self.stub.script(
                ("SearchExtraTools", {"query": "select:mcp__fixture__increment"}),
                (
                    "ExecuteExtraTool",
                    {"tool_name": "mcp__fixture__increment", "params": {}},
                ),
                "synthetic qualification complete",
            )
            self.rpc(
                "session/prompt",
                {
                    "sessionId": self.session,
                    "prompt": [
                        {"type": "text", "text": "A0 synthetic qualification only"}
                    ],
                },
                timeout=timeout,
            )
            if (
                self.stub.failures
                or not self.approvals
                or any(item["allowed"] for item in self.approvals)
                or self.counter
                or self.lookups
                or self.tool_arrivals
            ):
                raise Rejected("ACTUAL_DENIAL_NOT_VERIFIED")
            with self.gate:
                if (
                    self.state != "blocked"
                    or self.process.poll() is not None
                    or self.errors
                ):
                    raise Rejected("QUALIFICATION_CLOSED")
                self.state = "ready"
                self.ready_seconds = time.monotonic() - self.started_at
        except BaseException:
            self.invalidate("INITIALIZATION_FAILED")
            raise

    def prompt(self, text="A0 synthetic probe", *, timeout=45):
        return self.invoke(
            "session/prompt",
            {"sessionId": self.session, "prompt": [{"type": "text", "text": text}]},
            timeout=timeout,
        )

    def evidence(self):
        return {
            "instance": self.instance,
            "generation": self.generation,
            "state": self.state,
            "session": self.session,
            "argv": self.argv,
            "ready_seconds": self.ready_seconds,
            "closed_seconds": self.closed_seconds,
            "model_requests": self.stub.requests,
            "stub_failures": self.stub.failures,
            "counter": self.counter,
            "lookups": self.lookups,
            "tool_call_arrivals": self.tool_arrivals,
            "approvals": self.approvals,
            "denied_callbacks": self.denied_callbacks,
            "transport_errors": self.errors,
            "events": self.events,
            "cleanup_errors": self.cleanup_errors,
            "resources": self.resources,
            "exit_code": None if self.process is None else self.process.poll(),
        }

    def close(self):
        if self.closed_resources:
            return
        self.invalidate("CLOSED_BY_HOST")
        if self.process is not None:
            try:
                self.process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self.cleanup_errors.append("PROCESS_DID_NOT_EXIT")
            for stream in [self.process.stdin, self.process.stdout]:
                try:
                    stream.close()
                except (BrokenPipeError, OSError):
                    pass
        if self.reader is not None and self.reader is not threading.current_thread():
            self.reader.join(timeout=5)
            if self.reader.is_alive():
                self.cleanup_errors.append("ACP_READER_NOT_STOPPED")
        try:
            if self.stub is not None:
                self.stub.close()
        finally:
            if self.observer is not None:
                self.resources = self.observer.close()
                if self.resources["still_present"]:
                    self.cleanup_errors.append("DESCENDANT_STILL_PRESENT")
            if self.log is not None:
                self.log.close()
            self.closed_resources = True
