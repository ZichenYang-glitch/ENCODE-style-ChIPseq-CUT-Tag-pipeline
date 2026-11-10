"""Run A0 against one pinned native binary; save synthetic evidence, fail closed."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import socket
import subprocess
import time

from configuration import CONTRACT, PERI_COMMIT
from resource_guard import inspect_limits
from protocol import Client, Rejected
from sandbox import ROOT, command


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def denied(call):
    try:
        call()
    except Rejected:
        return
    raise AssertionError("expected fail-closed rejection")


def os_probe(binary, directory):
    directory.mkdir()
    state = directory / "state"
    state.mkdir()
    canary = directory / "bridge-canary"
    canary.write_text("SYNTHETIC_BRIDGE_SECRET")
    other = directory / "other-session"
    other.mkdir()
    (other / "canary").write_text("SYNTHETIC_OTHER_SESSION")
    addresses = [directory / "model.sock", directory / "platform.sock"]
    sockets = []
    try:
        for address in addresses:
            listener = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
            listener.bind(str(address))
            listener.listen(1)
            sockets.append(listener)
        tcp = socket.socket()
        tcp.bind(("127.0.0.1", 0))
        tcp.listen(1)
        sockets.append(tcp)
        abstract = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        abstract_name = "a0-probe-" + str(os.getpid())
        abstract.bind("\0" + abstract_name)
        abstract.listen(1)
        sockets.append(abstract)
        (state / "probe.json").write_text(
            json.dumps(
                {
                    "forbidden_files": {
                        "bridge_file_denied": str(canary),
                        "other_session_denied": str(other / "canary"),
                    },
                    "tcp_port": tcp.getsockname()[1],
                    "socket_path": str(addresses[1]),
                    "abstract": abstract_name,
                    "host_pid": os.getpid(),
                }
            )
        )
        result = subprocess.run(
            command(binary, state, addresses[0], probe=ROOT / "os_probe.py"),
            env={"A0_BRIDGE_SECRET": "synthetic"},
            close_fds=True,
            capture_output=True,
            text=True,
            timeout=15,
        )
        require(
            result.returncode == 0,
            "OS isolation failed: " + result.stderr + result.stdout,
        )
        evidence = json.loads(result.stdout)
        require(all(evidence["checks"].values()), "OS isolation check failed")
        return evidence
    finally:
        for listener in sockets:
            listener.close()
        for address in addresses:
            address.unlink(missing_ok=True)


def main():
    from closeout_cases import (
        invalid_policy,
        prediction_paths,
        missing_broker,
        changed_catalog,
        native_resource_failure,
    )
    from cases import (
        permission_and_tools,
        independent_instances,
        restart,
        initialization_failure,
        interrupted_permission,
        transport_failure,
        minimum_egress,
    )

    parser = argparse.ArgumentParser()
    parser.add_argument("--build-manifest", type=Path, required=True)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--qualification-only", action="store_true")
    parser.add_argument("--targeted", action="store_true")
    args = parser.parse_args()
    limits = inspect_limits()
    build = json.loads(args.build_manifest.read_text())
    source = subprocess.check_output(
        ["git", "-C", str(args.source), "rev-parse", "HEAD"], text=True
    ).strip()
    require(
        source == PERI_COMMIT == build["source_commit"], "wrong Peri source identity"
    )
    require(
        build["success"] and bool(build["local_patch"]), "unqualified candidate build"
    )
    binary = Path(build["binary"])
    binary_digest = hashlib.file_digest(binary.open("rb"), "sha256").hexdigest()
    require(
        binary_digest == build["binary_sha256"], "binary does not match recorded build"
    )
    require(
        subprocess.check_output(
            ["git", "-C", str(args.source), "diff", "HEAD"], text=True
        )
        == build["local_patch"],
        "unrecorded source patch",
    )
    args.evidence.mkdir(parents=True, exist_ok=False)
    report = {
        "contract": CONTRACT,
        "aggregate_resource_limits": limits,
        "patch_sha256": build["patch_sha256"],
        "helixweave_baseline": "622778f64e59f935d822d92a91bade78807dabce",
        "source_commit": source,
        "binary_sha256": binary_digest,
        "build_manifest_sha256": hashlib.sha256(
            args.build_manifest.read_bytes()
        ).hexdigest(),
        "prototype_sha256": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(ROOT.glob("*.py"))
        },
        "started_at_unix": time.time(),
        "passed": False,
        "cases": {},
    }
    all_clients = []

    def save():
        (args.evidence / "report.json").write_text(json.dumps(report, indent=2) + "\n")

    def run_case(name, function):
        clients = []
        case_dir = args.evidence / name
        case_dir.mkdir()

        def make(label, **kwargs):
            client = Client(binary, case_dir / label, **kwargs)
            clients.append(client)
            all_clients.append(client)
            return client

        result = {"status": "failed"}
        try:
            result["detail"] = function(make)
            result["status"] = "passed"
        except BaseException as error:
            result["error"] = type(error).__name__ + ": " + str(error)
        finally:
            for client in clients:
                try:
                    client.close()
                except BaseException as error:
                    result["status"] = "failed"
                    result.setdefault("cleanup_failures", []).append(
                        type(error).__name__
                    )
                try:
                    evidence = client.evidence()
                    (client.directory / "evidence.json").write_text(
                        json.dumps(evidence, indent=2) + "\n"
                    )
                except BaseException as error:
                    result["status"] = "failed"
                    result.setdefault("evidence_failures", []).append(
                        type(error).__name__
                    )
                if client.cleanup_errors:
                    result["status"] = "failed"
                    result.setdefault("cleanup_failures", []).extend(
                        client.cleanup_errors
                    )
            report["cases"][name] = result
            save()
            print(name + ": " + result["status"], flush=True)
        return result["status"] == "passed"

    try:
        report["os_probe"] = os_probe(binary, args.evidence / "os-probe")
        for fault in [
            "unknown_mode",
            "mode_error",
            "mode_timeout",
            "missing_tool",
            "extra_tool",
            "replacement_tool",
        ]:
            run_case(
                "A_" + fault,
                lambda make, fault=fault: initialization_failure(make, fault),
            )
        for fault in [
            "unknown_identity",
            "wrong_schema",
            "unknown_field",
            "prediction_typo",
        ]:
            run_case(
                "A_policy_" + fault,
                lambda make, fault=fault: invalid_policy(make, fault),
            )
        qualified = run_case("A_B_permission_tools", permission_and_tools)
        dependent = [
            ("A_prediction_paths", prediction_paths),
            ("A_missing_broker", missing_broker),
            ("B_changed_catalog", changed_catalog),
            ("F_native_resource", native_resource_failure),
            ("C_independent_instances", independent_instances),
            ("D_restart", restart),
            ("E_minimum_egress", minimum_egress),
        ]
        dependent.extend(
            (
                "A_approval_" + fault,
                lambda make, fault=fault: interrupted_permission(make, fault),
            )
            for fault in ["disconnect", "timeout"]
        )
        dependent.extend(
            ("F_" + fault, lambda make, fault=fault: transport_failure(make, fault))
            for fault in ["malformed", "timeout", "callback_disconnect"]
        )
        for name, function in dependent:
            if args.qualification_only or (
                args.targeted and name != "A_prediction_paths"
            ):
                report["cases"][name] = {
                    "status": "not_measured",
                    "reason": "Qualification-only invocation; separate full run required",
                }
            elif qualified:
                run_case(name, function)
            else:
                report["cases"][name] = {
                    "status": "not_measured",
                    "reason": "Native qualification failed; no ready session permitted",
                }
        report["phase"] = "targeted" if args.targeted else "full"
        report["passed"] = all(
            c["status"] == "passed" or (args.targeted and c["status"] == "not_measured")
            for c in report["cases"].values()
        )
    except BaseException as error:
        report["error"] = type(error).__name__ + ": " + str(error)
    finally:
        for client in all_clients:
            try:
                client.close()
            except BaseException as error:
                report.setdefault("final_cleanup_failures", []).append(
                    type(error).__name__
                )
        report["finished_at_unix"] = time.time()
        report["all_processes_exited"] = all(
            c.process.poll() is not None for c in all_clients
        )
        report["passed"] = (
            report["passed"]
            and report["all_processes_exited"]
            and not report.get("final_cleanup_failures")
        )
        save()
        print(
            json.dumps(
                {
                    "passed": report["passed"],
                    "all_processes_exited": report["all_processes_exited"],
                },
                indent=2,
            )
        )
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
