"""Native A0 scenarios. Each scenario has separate evidence and prerequisite errors."""

import json
import os
import signal
import time

from configuration import TOOLS
from fixtures import MARKERS, hostile_document, project_fixture
from protocol import Rejected


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def denied(call):
    try:
        call()
    except Rejected as error:
        return str(error)
    raise AssertionError("expected fail-closed rejection")


def tool_script(name, arguments):
    return [
        ("SearchExtraTools", {"query": "select:mcp__fixture__" + name}),
        (
            "ExecuteExtraTool",
            {"tool_name": "mcp__fixture__" + name, "params": arguments},
        ),
        "synthetic done",
    ]


def permission_and_tools(make):
    client = make("permissions")
    denied(lambda: client.prompt())
    require(not client.stub.requests, "blocked prompt reached model")
    client.initialize()
    client.stub.script(*tool_script("increment", {}))
    client.prompt()
    require(
        client.counter == 0 and client.tool_arrivals == 0,
        "Peri sent tools/call despite approval reject (broker masking is not a pass)",
    )
    require(
        any(not item["allowed"] for item in client.approvals),
        "Peri approval not observed",
    )
    client.approve = True
    client.stub.script(*tool_script("increment", {}))
    client.prompt()
    require(
        client.counter == 1 and client.tool_arrivals == 1,
        "approved tool not executed once",
    )
    client.stub.script(*tool_script("lookup", {"object_id": "item-1"}))
    client.prompt()
    require(client.lookups == 1, "allowed lookup failed")
    # Every attempt goes through the real model/tool path. These are not shell commands.
    for target, arguments in [
        ("Bash", {"command": "read synthetic secret"}),
        ("Read", {"file_path": "/private/a0/secret"}),
        ("SQL", {"query": "synthetic query"}),
        ("mcp__unregistered__write", {}),
        ("session/set_mode", {"modeId": "bypass"}),
        ("Agent", {"prompt": "synthetic escalation"}),
        ("Workflow", {}),
    ]:
        before = client.tool_arrivals
        client.stub.script(
            ("ExecuteExtraTool", {"tool_name": target, "params": arguments}),
            "synthetic refused",
        )
        client.prompt()
        require(client.tool_arrivals == before, "unregistered target reached fake MCP")
    for target in [
        "ReadContextArchive",
        "readcontextarchive",
        "DiscoverMCP",
        "mcp_read_resource",
        "meta.ReadContextArchive",
    ]:
        before = client.tool_arrivals
        client.stub.script((target, {}), "synthetic denied direct call")
        client.prompt()
        require(client.tool_arrivals == before, "forbidden direct tool reached broker")
        client.stub.script(
            ("ExecuteExtraTool", {"tool_name": target, "params": {}}),
            "synthetic denied wrapped call",
        )
        client.prompt()
        require(client.tool_arrivals == before, "forbidden wrapped tool reached broker")
    client.stub.script(
        (
            "SearchExtraTools",
            {"query": "ReadContextArchive DiscoverMCP mcp_read_resource"},
        ),
        "synthetic discovery checked",
    )
    client.prompt()
    for arguments in [{"extra": True}, {"object_id": "different-project"}]:
        client.stub.script(
            (
                "ExecuteExtraTool",
                {"tool_name": "mcp__fixture__increment", "params": arguments},
            ),
            "synthetic invalid",
        )
        client.prompt()
    require(
        client.counter == 1 and client.lookups == 1,
        "out-of-contract arguments executed",
    )
    for request in client.stub.requests:
        names = {
            entry["function"]["name"]
            for entry in request.get("body", {}).get("tools", [])
        }
        require(
            names == {"SearchExtraTools", "ExecuteExtraTool"},
            "unexpected actual direct tools: " + str(names),
        )
    wire = json.dumps([r.get("body") for r in client.stub.requests])
    require(
        "mcp__fixture__increment" in wire
        and ("inputSchema" in wire or "parameters" in wire),
        "deferred discovery not present in provider context",
    )
    require(not client.stub.failures, "unexpected provider traffic")
    return {
        "rejected_peri_tool_arrivals": 0,
        "allowed_increments": client.counter,
        "allowed_lookups": client.lookups,
    }


def independent_instances(make):
    first, second = make("instance-a"), make("instance-b")
    first.initialize()
    second.initialize()
    require(first.session != second.session, "same session identity")
    # Internal test control deliberately changes A only; not exposed through invoke().
    first.rpc("session/set_mode", {"sessionId": first.session, "modeId": "bypass"})
    # Do not set B back to Default: that would mask shared-mode leakage.
    second.stub.script(*tool_script("increment", {}))
    second.prompt()
    require(
        second.tool_arrivals == 0 and second.counter == 0,
        "A bypass changed B rejection",
    )
    first.rpc("session/set_mode", {"sessionId": first.session, "modeId": "default"})
    first.approve = True
    first.stub.script(*tool_script("increment", {}))
    first.prompt()
    second.stub.script(*tool_script("increment", {}))
    second.prompt()
    require(
        first.counter == 1 and second.counter == 0, "instance effects not independent"
    )
    for params, kwargs in [
        ({"sessionId": second.session, "prompt": []}, {}),
        ({"sessionId": first.session, "prompt": []}, {"instance": second.instance}),
        ({"sessionId": first.session, "prompt": []}, {"generation": 0}),
    ]:
        denied(
            lambda params=params, kwargs=kwargs: first.invoke(
                "session/prompt", params, **kwargs
            )
        )
    for method in [
        "session/load",
        "session/resume",
        "session/fork",
        "session/set_mode",
        "session/update_config",
        "peri/mcp/control",
    ]:
        denied(lambda method=method: first.invoke(method, {"sessionId": first.session}))
    callback = {
        "serverId": first.server_id,
        "requestId": "forged",
        "method": "tools/call",
        "params": {"name": "increment", "arguments": {}},
    }
    for instance, generation in [(second.instance, 1), (first.instance, 0)]:
        result = first.callback("mcp/message", callback, instance, generation)
        require("error" in result, "cross-instance or stale callback accepted")
    return {"a_counter": first.counter, "b_counter": second.counter}


def restart(make):
    old = make("old", generation=1)
    old.initialize()
    old_session = old.session
    os.killpg(old.process.pid, signal.SIGKILL)
    old.process.wait(timeout=5)
    deadline = time.monotonic() + 2
    while old.state != "closed" and time.monotonic() < deadline:
        time.sleep(0.01)
    require(old.state == "closed", "process exit did not close client")
    denied(lambda: old.prompt())
    new = make("new", generation=2)
    denied(lambda: new.prompt())
    require(not new.stub.requests, "new blocked instance made model request")
    new.initialize()
    denied(
        lambda: new.invoke(
            "session/prompt",
            {"sessionId": old_session, "prompt": []},
            instance=old.instance,
            generation=1,
        )
    )
    new.stub.script("fresh synthetic session")
    new.prompt()
    return {"old_closed": old.state, "new_generation": new.generation}


def initialization_failure(make, fault):
    tools = TOOLS
    if fault == "missing_tool":
        tools = TOOLS[:1]
    elif fault == "extra_tool":
        tools = TOOLS + [{"name": "extra", "inputSchema": {"type": "object"}}]
    elif fault == "replacement_tool":
        tools = [TOOLS[0], {"name": "Read", "inputSchema": {"type": "object"}}]
    client = make("instance", tools=tools)
    original = client._send
    injected = []
    if fault == "mode_error":

        def bad_mode(value):
            if value.get("method") == "session/set_mode":
                injected.append("mode_error")
                value["params"]["sessionId"] = "not-a-real-session"
            original(value)

        client._send = bad_mode
    elif fault == "mode_timeout":

        def lose_mode(value):
            if value.get("method") != "session/set_mode":
                original(value)
            else:
                injected.append("mode_timeout")

        client._send = lose_mode
    code = denied(
        lambda: client.initialize(
            mode="unknown" if fault == "unknown_mode" else "default",
            timeout=1 if fault == "mode_timeout" else 15,
        )
    )
    expected = {
        "unknown_mode": "MODE_NOT_ALLOWED",
        "mode_error": "ACP_REMOTE_REJECTION",
        "mode_timeout": "RPC_TIMEOUT",
        "missing_tool": "TOOL_CATALOG_MISMATCH",
        "extra_tool": "TOOL_CATALOG_MISMATCH",
        "replacement_tool": "TOOL_CATALOG_MISMATCH",
    }[fault]
    require(code == expected, "intended failure boundary not reached: " + code)
    if fault in {"mode_error", "mode_timeout"}:
        require(injected == [fault], "mode fault was never injected")
    require(client.state == "closed", "initialization failure remained usable")
    require(
        not client.stub.requests and client.tool_arrivals == 0,
        "work occurred on failed startup",
    )
    denied(lambda: client.prompt())
    return {
        "fault": fault,
        "stable_error": code,
        "injected": injected,
        "model_requests": 0,
        "tool_arrivals": 0,
    }


def interrupted_permission(make, fault):
    client = make("instance")
    client.initialize()
    client.approve = True
    client.approval_failure = fault
    baseline = len(client.approvals)
    client.stub.script(*tool_script("increment", {}))
    code = denied(lambda: client.prompt(timeout=2))
    require(len(client.approvals) > baseline, "new approval boundary not reached")
    expected = "APPROVAL_DISCONNECTED" if fault == "disconnect" else "RPC_TIMEOUT"
    require(expected in client.errors, "expected approval failure did not occur")
    require(
        client.counter == 0 and client.tool_arrivals == 0,
        "interrupted approval reached tool",
    )
    require(client.state == "closed", "interrupted permission did not close")
    return {"stable_error": code, "tool_arrivals": 0, "effects": 0}


def transport_failure(make, fault):
    client = make("instance")
    client.initialize()
    if fault == "malformed":
        # Fault injection at the same parser called by the real reader thread.
        try:
            client.consume_line(b"{invalid-json}\n")
        except Rejected as error:
            client.invalidate(str(error))
        code = denied(lambda: client.prompt())
    elif fault == "timeout":
        client.stub.script(None)
        code = denied(lambda: client.prompt(timeout=1))
    else:
        client.approve = True
        client.callback_failure = "disconnect"
        client.stub.script(*tool_script("increment", {}))
        code = denied(lambda: client.prompt(timeout=3))
        require(client.tool_arrivals == 1, "tool callback disconnect seam not reached")
    require(
        client.state == "closed" and client.counter == 0, "failure did not stop safely"
    )
    return {"fault": fault, "stable_error": code, "effects": 0}


def minimum_egress(make):
    client = make("minimum-egress")
    client.initialize()
    client.stub.script("safe synthetic explanation")
    projected = project_fixture(hostile_document())
    client.prompt("Explain this closed synthetic projection: " + json.dumps(projected))
    wire = json.dumps(client.stub.requests)
    for marker in MARKERS:
        require(marker not in wire, "forbidden synthetic field reached provider")
    require("INPUT_INVALID" in wire, "test did not reach actual provider")
    require(not client.stub.failures, "unexpected outbound HTTP")
    return {
        "forbidden_markers_absent": MARKERS,
        "scope": "closed structured fixture; no arbitrary free-text confidentiality guarantee",
    }
