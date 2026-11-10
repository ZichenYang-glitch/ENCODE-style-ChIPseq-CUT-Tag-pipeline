"""Additional checks for the explicitly authorized closeout gaps."""

import copy
import json
import os
from pathlib import Path
import resource

from configuration import TOOLS, process_policy
from cases import denied, require, tool_script


def invalid_policy(make, fault):
    policy = process_policy()
    prediction = "off"
    if fault == "unknown_identity":
        policy["tools"].append(dict(policy["tools"][0], name="unregistered_identity"))
    elif fault == "wrong_schema":
        policy["tools"][-1]["parameters"] = {"type": "integer"}
    elif fault == "unknown_field":
        policy["unexpected"] = True
    elif fault == "prediction_typo":
        prediction = "false"
    client = make("instance", policy=policy, prediction=prediction)
    code = denied(client.initialize)
    require(
        client.state == "closed"
        and not client.stub.requests
        and client.tool_arrivals == 0,
        "invalid startup produced work",
    )
    return {"fault": fault, "error": code, "provider_requests": 0, "tool_arrivals": 0}


def prediction_paths(make):
    client = make("instance")
    client.initialize()
    before = client.verify_runtime()
    client.stub.script("successful synthetic turn")
    client.prompt("normal synthetic turn")
    client.stub.script(("UnknownTool", {}), "synthetic tool error handled")
    client.prompt("synthetic unknown tool")
    facts = client.verify_runtime()
    require(
        facts["prediction_suppressed"] == before["prediction_suppressed"] + 2,
        "end paths not observed",
    )
    require(facts["prediction_started"] == 0, "auxiliary task started")
    require(
        all("prediction_directive" not in json.dumps(q) for q in client.stub.requests),
        "prediction requested model",
    )
    for text in ["/compact", " /compact", "/goal synthetic"]:
        denied(lambda text=text: client.prompt(text))
    for method in [
        "peri/host_policy",
        "session/set_config_option",
        "session/set_mode",
        "session/load",
        "session/resume",
        "session/fork",
    ]:
        denied(
            lambda method=method: client.invoke(method, {"sessionId": client.session})
        )
    require(not client.stub.failures, "unexpected model path")
    return {
        "effective_policy": facts,
        "extra_completed_turns": 2,
        "model_requests": len(client.stub.requests),
    }


def missing_broker(make):
    client = make("instance")
    client.initialize()
    client.approve = True
    client.approval_failure = "missing"
    before = len(client.approvals)
    client.stub.script(*tool_script("increment", {}))
    client.prompt()
    require(len(client.approvals) == before + 1, "approval error seam not reached")
    require(
        client.counter == 0 and client.tool_arrivals == 0,
        "unavailable approval executed",
    )
    return {
        "new_approval_requests": 1,
        "tools_call_arrivals": 0,
        "effects": 0,
        "boundary": "real ACP permission RPC returns method-not-found; not Rust None broker construction",
    }


def changed_catalog(make):
    client = make("instance")
    client.initialize()
    before = len(client.stub.requests)
    client.declared_tools = copy.deepcopy(TOOLS)
    client.declared_tools[-1]["inputSchema"] = {"type": "integer"}
    code = denied(client.prompt)
    require(
        code == "LIVE_CATALOG_MISMATCH", "runtime broker catalog mismatch not rejected"
    )
    require(
        len(client.stub.requests) == before and client.tool_arrivals == 0,
        "work after changed broker catalog",
    )
    return {
        "code": code,
        "provider_delta": 0,
        "boundary": "broker contract before prompt; native generation rejection separately tested in Rust",
    }


def native_resource_failure(make):
    client = make("instance")
    client.initialize()
    # Traverse only this instance's process tree; verify the exact mounted binary.
    pending = [client.process.pid]
    native = None
    while pending:
        pid = pending.pop()
        try:
            base = Path("/proc") / str(pid)
            pending.extend(
                map(int, (base / "task" / str(pid) / "children").read_text().split())
            )
            executable = os.readlink(base / "exe")
            if executable == "/runtime/peri" or executable.endswith("/release/peri"):
                native = pid
        except (FileNotFoundError, ProcessLookupError):
            continue
    require(native is not None, "owned native PID not located")
    before = len(client.stub.requests)
    old = resource.prlimit(native, resource.RLIMIT_NOFILE)
    resource.prlimit(native, resource.RLIMIT_NOFILE, (0, old[1]))
    client.stub.script("should not reach provider")
    code = denied(lambda: client.prompt(timeout=3))
    require(
        client.state == "closed" and client.counter == 0 and client.tool_arrivals == 0,
        "resource error left usable client/effect",
    )
    require(
        len(client.stub.requests) == before,
        "native opened provider connection after FD exhaustion",
    )
    return {
        "native_pid": native,
        "old_nofile": old,
        "injected_soft_limit": 0,
        "error": code,
        "provider_delta": 0,
        "boundary": "native RLIMIT_NOFILE exhaustion, not cgroup OOM",
    }
