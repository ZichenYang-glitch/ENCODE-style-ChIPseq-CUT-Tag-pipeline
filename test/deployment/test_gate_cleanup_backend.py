from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import runpy
import shutil
import signal
import stat

import pytest

from encode_pipeline.deployment.canonical import (
    canonical_identity,
    canonical_json_bytes,
)
from encode_pipeline.deployment.gate import (
    CLEANUP_PLAN_SCHEMA,
    CLEANUP_PLAN_IDENTITY_SCHEME,
    GATE_PROCESS_IDENTITY_SCHEME,
)


ROOT = Path(__file__).resolve().parents[2]
BACKEND = (
    ROOT
    / "src"
    / "encode_pipeline"
    / "deployment"
    / "templates"
    / "helixweave-gate-cleanup"
)
LAUNCHER = BACKEND.with_name("helixweave-gate-cleanup-launcher")
TASK = "task-0123456789abcdef0123456789abcdef"
OBSERVATION = f"sha256-{1:064x}"
RESOURCE_SCHEME = "helixweave-deployment-gate-resource-path-v1"
EXECUTOR_SCHEME = "helixweave-deployment-gate-cleanup-executor-path-v1"


def _mkdir(path: Path, mode: int = 0o700) -> Path:
    path.mkdir(mode=mode, parents=True, exist_ok=True)
    path.chmod(mode)
    return path


def _executable(path: Path, content: bytes) -> Path:
    _mkdir(path.parent)
    path.write_bytes(content)
    path.chmod(0o555)
    return path


def _file_identity(path: Path) -> str:
    return f"sha256-{hashlib.sha256(path.read_bytes()).hexdigest()}"


def _resource(kind: str, disposition: str, path: Path) -> dict[str, object]:
    observed = path.stat(follow_symlinks=False)
    return {
        "kind": kind,
        "disposition": disposition,
        "absolute_path": str(path),
        "path_identity": canonical_identity(
            {"kind": kind, "absolute_path": str(path)},
            scheme=RESOURCE_SCHEME,
        ),
        "device": observed.st_dev,
        "inode": observed.st_ino,
    }


def _proc_stat(pid: int, name: str, start_ticks: int) -> bytes:
    fields = [b"S", *(b"1" for _ in range(18)), str(start_ticks).encode("ascii")]
    return f"{pid} ({name}) ".encode("ascii") + b" ".join(fields) + b"\n"


def _process(
    *,
    proc_root: Path,
    task_root: Path,
    name: str,
    pid: int,
    start_ticks: int,
    socket_path: Path | None,
) -> dict[str, object]:
    process = _mkdir(proc_root / str(pid))
    fd_root = _mkdir(process / "fd")
    (process / "stat").write_bytes(_proc_stat(pid, name, start_ticks))
    executable = process / "exe"
    executable.write_bytes(f"{name}-executable\n".encode("ascii"))
    cmdline = f"/fixed/{name}".encode("ascii") + b"\0" + TASK.encode("ascii") + b"\0"
    (process / "cmdline").write_bytes(cmdline)
    pid_file = task_root / "pids" / f"{name}.pid"
    pid_file.write_text(f"{pid}\n", encoding="ascii")
    pid_file.chmod(0o600)
    executable_stat = executable.stat()
    socket_stat = socket_path.stat(follow_symlinks=False) if socket_path else None
    socket_kernel_inode = pid * 10 + 7 if socket_path else None
    if socket_path is not None:
        (fd_root / "3").symlink_to(f"socket:[{socket_kernel_inode}]")
        proc_net_unix = proc_root / "net/unix"
        if not proc_net_unix.exists():
            _mkdir(proc_net_unix.parent)
            proc_net_unix.write_text(
                "Num RefCount Protocol Flags Type St Inode Path\n",
                encoding="utf-8",
            )
        with proc_net_unix.open("a", encoding="utf-8") as stream:
            stream.write(
                "00000000: 00000002 00000000 00010000 "
                f"0001 01 {socket_kernel_inode} {socket_path}\n"
            )
    value: dict[str, object] = {
        "name": name,
        "task_identity": TASK,
        "pid": pid,
        "process_start_ticks": start_ticks,
        "executable_device": executable_stat.st_dev,
        "executable_inode": executable_stat.st_ino,
        "cmdline_identity": f"sha256-{hashlib.sha256(cmdline).hexdigest()}",
        "socket_device": socket_stat.st_dev if socket_stat else None,
        "socket_inode": socket_stat.st_ino if socket_stat else None,
        "socket_kernel_inode": socket_kernel_inode,
    }
    value["identity"] = canonical_identity(value, scheme=GATE_PROCESS_IDENTITY_SCHEME)
    return value


def _bind(path: Path) -> None:
    _mkdir(path.parent)
    path.write_bytes(b"test-only socket witness\n")
    path.chmod(0o600)
    return None


def _layout(tmp_path: Path) -> dict[str, object]:
    host = _mkdir(tmp_path / "host")
    gate = _mkdir(host / "var/lib/helixweave/operator/deployment-gates")
    task_root = _mkdir(gate / "tasks" / TASK)
    _mkdir(task_root / "pids")
    for relative in ("checkout", "venv", "redis", "runner/prepared"):
        directory = _mkdir(task_root / relative)
        (directory / "owned-evidence").write_bytes(b"delete\n")
    docker = _mkdir(task_root / "docker")
    docker_data = _mkdir(gate / "docker-data")
    history = _mkdir(gate / "evidence/history")
    (docker_data / "retained").write_bytes(b"docker\n")
    (history / "retained").write_bytes(b"history\n")

    redis_socket = _bind(task_root / "redis/redis.sock")
    docker_socket = _bind(docker / "docker.sock")
    proc_root = _mkdir(host / "proc")
    processes = [
        _process(
            proc_root=proc_root,
            task_root=task_root,
            name="dockerd",
            pid=1001,
            start_ticks=9001,
            socket_path=docker / "docker.sock",
        ),
        _process(
            proc_root=proc_root,
            task_root=task_root,
            name="redis",
            pid=1002,
            start_ticks=9002,
            socket_path=task_root / "redis/redis.sock",
        ),
        _process(
            proc_root=proc_root,
            task_root=task_root,
            name="runner",
            pid=1003,
            start_ticks=9003,
            socket_path=None,
        ),
    ]

    executor = _executable(
        host / "usr/libexec/helixweave-gate-cleanup", b"stable launcher\n"
    )
    closure = f"sha256-{2:064x}"
    selected_backend = _executable(
        host / "opt/helixweave/operator" / closure / "bin/helixweave-gate-cleanup",
        BACKEND.read_bytes(),
    )
    executor_stat = executor.stat()
    backend_stat = selected_backend.stat()
    executor_identity = _file_identity(executor)
    backend_identity = _file_identity(selected_backend)

    targets = [
        _resource("checkout", "delete", task_root / "checkout"),
        _resource("docker-data-root", "retain", docker_data),
        _resource("docker-images", "retain", docker_data),
        _resource("dockerd-socket", "delete", docker / "docker.sock"),
        _resource("historical-evidence", "retain", history),
        _resource("process-state", "delete", task_root / "pids"),
        _resource("redis-state", "delete", task_root / "redis"),
        _resource("runner-state", "delete", task_root / "runner"),
        _resource("venv", "delete", task_root / "venv"),
    ]
    plan: dict[str, object] = {
        "schema_version": CLEANUP_PLAN_SCHEMA,
        "task_identity": TASK,
        "observation_identity": OBSERVATION,
        "targets": targets,
        "processes": processes,
        "executor_path": str(executor),
        "executor": {
            "path_identity": canonical_identity(
                {"absolute_path": str(executor)}, scheme=EXECUTOR_SCHEME
            ),
            "file_identity": executor_identity,
            "device": executor_stat.st_dev,
            "inode": executor_stat.st_ino,
            "closure_identity": closure,
            "backend_identity": backend_identity,
            "backend_device": backend_stat.st_dev,
            "backend_inode": backend_stat.st_ino,
        },
    }
    plan["identity"] = canonical_identity(plan, scheme=CLEANUP_PLAN_IDENTITY_SCHEME)
    plan_path = task_root / "runner/prepared/cleanup-plan.json"
    plan_path.write_bytes(canonical_json_bytes(plan))
    plan_path.chmod(0o600)
    return {
        "host": host,
        "gate": gate,
        "task_root": task_root,
        "proc_root": proc_root,
        "executor": executor,
        "selected_backend": selected_backend,
        "closure": closure,
        "executor_identity": executor_identity,
        "backend_identity": backend_identity,
        "plan": plan,
        "plan_path": plan_path,
        "sockets": (redis_socket, docker_socket),
        "docker_data": docker_data,
        "history": history,
    }


def _module(layout: dict[str, object]) -> dict[str, object]:
    namespace = runpy.run_path(str(BACKEND))
    globals_ = namespace["main"].__globals__
    globals_.update(
        {
            "HOST_ROOT": layout["host"],
            "GATE_ROOT": layout["gate"],
            "EXECUTOR": layout["executor"],
            "OPERATOR_ROOT": Path(layout["selected_backend"]).parents[2],
            "SELECTED_BACKEND": layout["selected_backend"],
            "PROC_ROOT": layout["proc_root"],
            "ROOT_UID": os.getuid(),
            "ROOT_GID": os.getgid(),
            "GETEUID": lambda: 0,
            "PIDFD_OPEN": None,
            "PIDFD_SEND_SIGNAL": None,
            "POLL_SECONDS": 0.0,
            "SOCKET_MATCHER": stat.S_ISREG,
            "MOUNT_ID": lambda _descriptor: 1,
        }
    )
    return namespace


def _arguments(layout: dict[str, object]) -> tuple[str, ...]:
    plan = layout["plan"]
    return (
        "cleanup",
        "--task-id",
        TASK,
        "--plan-id",
        plan["identity"],
        "--executor-id",
        layout["executor_identity"],
        "--closure-id",
        layout["closure"],
        "--backend-id",
        layout["backend_identity"],
    )


def _rewrite_plan(layout: dict[str, object]) -> None:
    plan = layout["plan"]
    assert isinstance(plan, dict)
    plan.pop("identity", None)
    plan["identity"] = canonical_identity(
        plan,
        scheme=CLEANUP_PLAN_IDENTITY_SCHEME,
    )
    plan_path = Path(layout["plan_path"])
    plan_path.write_bytes(canonical_json_bytes(plan))
    plan_path.chmod(0o600)


def _close_sockets(layout: dict[str, object]) -> None:
    for value in layout["sockets"]:
        if value is not None:
            value.close()


def test_cleanup_stops_exact_processes_deletes_only_task_targets_and_retains_data(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    layout = _layout(tmp_path)
    namespace = _module(layout)
    calls: list[tuple[int, signal.Signals]] = []

    def stop(pid: int, selected_signal: signal.Signals) -> None:
        calls.append((pid, selected_signal))
        shutil.rmtree(Path(layout["proc_root"]) / str(pid))

    namespace["main"].__globals__["KILL"] = stop
    environment = dict(os.environ)
    try:
        assert namespace["main"](_arguments(layout)) == 0
    finally:
        os.environ.clear()
        os.environ.update(environment)
        _close_sockets(layout)

    receipt_text = capsys.readouterr().out
    receipt = json.loads(receipt_text)
    assert receipt_text == canonical_json_bytes(receipt).decode("utf-8")
    assert set(receipt) == {
        "deleted_kinds",
        "identity",
        "plan_identity",
        "retained_kinds",
        "schema_version",
        "status",
        "task_identity",
        "terminated_process_identities",
    }
    assert receipt["status"] == "cleaned"
    assert receipt["plan_identity"] == layout["plan"]["identity"]
    assert all(item[1] == signal.SIGTERM for item in calls)
    assert [item[0] for item in calls] == [1003, 1002, 1001]
    task_root = Path(layout["task_root"])
    assert not any(
        (task_root / relative).exists()
        for relative in ("checkout", "venv", "redis", "runner", "pids")
    )
    assert not (task_root / "docker/docker.sock").exists()
    assert (Path(layout["docker_data"]) / "retained").read_bytes() == b"docker\n"
    assert (Path(layout["history"]) / "retained").read_bytes() == b"history\n"
    assert str(tmp_path) not in receipt_text


def test_launcher_exec_failure_is_canonical_path_free_and_unavailable(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    values = (
        "cleanup",
        "--task-id",
        TASK,
        "--plan-id",
        OBSERVATION,
        "--executor-id",
        OBSERVATION,
        "--closure-id",
        OBSERVATION,
        "--backend-id",
        OBSERVATION,
    )
    monkeypatch.setattr(
        os,
        "execve",
        lambda *_args: (_ for _ in ()).throw(OSError("/private/backend token=secret")),
    )
    monkeypatch.setattr("sys.argv", ["helixweave-gate-cleanup", *values])

    with pytest.raises(SystemExit) as stopped:
        runpy.run_path(str(LAUNCHER))

    captured = capsys.readouterr()
    receipt = json.loads(captured.err)
    assert stopped.value.code == 69
    assert captured.out == ""
    assert captured.err == canonical_json_bytes(receipt).decode("utf-8")
    assert receipt == {
        "issue": {
            "code": "GATE_CLEANUP_BACKEND_UNAVAILABLE",
            "message": "Gate cleanup backend is unavailable.",
        },
        "schema_version": "helixweave-deployment-gate-cleanup-receipt-v1",
        "status": "error",
    }
    assert captured.err.count("\n") == 1
    assert "/opt/" not in captured.err
    assert "private" not in captured.err
    assert "secret" not in captured.err
    assert "OSError" not in captured.err


def test_process_identity_mismatch_emits_no_signal_and_deletes_nothing(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    layout = _layout(tmp_path)
    namespace = _module(layout)
    calls: list[tuple[int, signal.Signals]] = []
    namespace["main"].__globals__["KILL"] = lambda pid, sig: calls.append((pid, sig))
    (Path(layout["proc_root"]) / "1002/cmdline").write_bytes(b"changed\0")
    environment = dict(os.environ)
    try:
        assert namespace["main"](_arguments(layout)) == 65
    finally:
        os.environ.clear()
        os.environ.update(environment)
        _close_sockets(layout)

    error = json.loads(capsys.readouterr().err)
    assert error["issue"]["code"] == "GATE_CLEANUP_PROCESS_MISMATCH"
    assert calls == []
    assert (Path(layout["task_root"]) / "checkout/owned-evidence").exists()


def test_fallback_treats_identity_drift_as_alive_and_kills_the_same_instance(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    layout = _layout(tmp_path)
    namespace = _module(layout)
    globals_ = namespace["main"].__globals__
    calls: list[tuple[int, signal.Signals]] = []
    clock = [0.0]

    def selected_kill(pid: int, selected_signal: signal.Signals) -> None:
        calls.append((pid, selected_signal))
        process = Path(layout["proc_root"]) / str(pid)
        if selected_signal == signal.SIGTERM:
            (process / "cmdline").write_bytes(b"drifted-after-term\0")
            (process / "exe").write_bytes(b"drifted-after-term\n")
        else:
            shutil.rmtree(Path(layout["proc_root"]) / str(pid))

    globals_["KILL"] = selected_kill
    globals_["MONOTONIC"] = lambda: clock[0]
    globals_["SLEEP"] = lambda duration: clock.__setitem__(0, clock[0] + duration + 20)
    environment = dict(os.environ)
    try:
        assert namespace["main"](_arguments(layout)) == 0
    finally:
        os.environ.clear()
        os.environ.update(environment)
        _close_sockets(layout)

    assert calls == [
        (1003, signal.SIGTERM),
        (1002, signal.SIGTERM),
        (1001, signal.SIGTERM),
        (1003, signal.SIGKILL),
        (1002, signal.SIGKILL),
        (1001, signal.SIGKILL),
    ]
    assert json.loads(capsys.readouterr().out)["status"] == "cleaned"


@pytest.mark.parametrize(
    "mutation",
    (
        "duplicate-plan-key",
        "hardlinked-plan",
        "nan",
        "noncanonical",
        "oversized",
        "symlink-target",
        "wrong-mode",
    ),
)
def test_untrusted_plan_or_target_fails_before_signals(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    mutation: str,
) -> None:
    layout = _layout(tmp_path)
    namespace = _module(layout)
    calls: list[tuple[int, signal.Signals]] = []
    namespace["main"].__globals__["KILL"] = lambda pid, sig: calls.append((pid, sig))
    if mutation == "duplicate-plan-key":
        content = Path(layout["plan_path"]).read_bytes()
        Path(layout["plan_path"]).write_bytes(
            b'{"schema_version":"duplicate",' + content[1:]
        )
        Path(layout["plan_path"]).chmod(0o600)
    elif mutation == "hardlinked-plan":
        (tmp_path / "plan-alias").hardlink_to(layout["plan_path"])
    elif mutation == "nan":
        Path(layout["plan_path"]).write_bytes(b'{"value":NaN}\n')
    elif mutation == "noncanonical":
        content = Path(layout["plan_path"]).read_bytes()
        Path(layout["plan_path"]).write_bytes(b" " + content)
    elif mutation == "oversized":
        Path(layout["plan_path"]).write_bytes(b"x" * (64 * 1024 + 1))
    elif mutation == "wrong-mode":
        Path(layout["plan_path"]).chmod(0o644)
    else:
        checkout = Path(layout["task_root"]) / "checkout"
        shutil.rmtree(checkout)
        checkout.symlink_to(layout["docker_data"], target_is_directory=True)
    environment = dict(os.environ)
    try:
        assert namespace["main"](_arguments(layout)) == 65
    finally:
        os.environ.clear()
        os.environ.update(environment)
        _close_sockets(layout)

    error = json.loads(capsys.readouterr().err)
    assert error["issue"]["code"] in {
        "GATE_CLEANUP_PLAN_INVALID",
        "GATE_CLEANUP_TARGET_MISMATCH",
    }
    assert calls == []
    assert (Path(layout["docker_data"]) / "retained").exists()
    assert stat.S_ISDIR(Path(layout["history"]).stat().st_mode)


def test_target_alias_in_a_canonical_plan_fails_before_signals(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    layout = _layout(tmp_path)
    namespace = _module(layout)
    calls: list[tuple[int, signal.Signals]] = []
    namespace["main"].__globals__["KILL"] = lambda pid, sig: calls.append((pid, sig))
    plan = layout["plan"]
    assert isinstance(plan, dict)
    targets = plan["targets"]
    assert isinstance(targets, list)
    checkout = next(item for item in targets if item["kind"] == "checkout")
    retained = next(item for item in targets if item["kind"] == "docker-data-root")
    checkout["device"] = retained["device"]
    checkout["inode"] = retained["inode"]
    _rewrite_plan(layout)
    environment = dict(os.environ)
    try:
        assert namespace["main"](_arguments(layout)) == 65
    finally:
        os.environ.clear()
        os.environ.update(environment)
        _close_sockets(layout)

    assert json.loads(capsys.readouterr().err)["issue"]["code"] == (
        "GATE_CLEANUP_PLAN_INVALID"
    )
    assert calls == []
    assert (Path(layout["task_root"]) / "checkout/owned-evidence").exists()


def test_socket_owner_mismatch_emits_no_signal_and_deletes_nothing(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    layout = _layout(tmp_path)
    namespace = _module(layout)
    calls: list[tuple[int, signal.Signals]] = []
    namespace["main"].__globals__["KILL"] = lambda pid, sig: calls.append((pid, sig))
    proc_net_unix = Path(layout["proc_root"]) / "net/unix"
    proc_net_unix.write_text(
        proc_net_unix.read_text(encoding="utf-8").replace("10017", "99999"),
        encoding="utf-8",
    )
    environment = dict(os.environ)
    try:
        assert namespace["main"](_arguments(layout)) == 65
    finally:
        os.environ.clear()
        os.environ.update(environment)
        _close_sockets(layout)

    assert json.loads(capsys.readouterr().err)["issue"]["code"] == (
        "GATE_CLEANUP_PROCESS_MISMATCH"
    )
    assert calls == []
    assert (Path(layout["task_root"]) / "checkout/owned-evidence").exists()


def test_pre_signal_identity_drift_emits_no_signal_and_deletes_nothing(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    layout = _layout(tmp_path)
    namespace = _module(layout)
    globals_ = namespace["main"].__globals__
    original = globals_["_preflight_processes"]
    calls: list[tuple[int, signal.Signals]] = []

    def drift_after_preflight(records: object) -> object:
        handles = original(records)
        (Path(layout["proc_root"]) / "1002/cmdline").write_bytes(b"drifted\0")
        return handles

    globals_["_preflight_processes"] = drift_after_preflight
    globals_["KILL"] = lambda pid, sig: calls.append((pid, sig))
    environment = dict(os.environ)
    try:
        assert namespace["main"](_arguments(layout)) == 65
    finally:
        os.environ.clear()
        os.environ.update(environment)
        _close_sockets(layout)

    assert json.loads(capsys.readouterr().err)["issue"]["code"] == (
        "GATE_CLEANUP_PROCESS_MISMATCH"
    )
    assert calls == []
    assert (Path(layout["task_root"]) / "checkout/owned-evidence").exists()


@pytest.mark.parametrize(
    ("exit_on", "expected_status"),
    (("term", 0), ("kill", 0), ("never", 69)),
)
def test_pidfd_wait_controls_term_kill_and_timeout_without_host_signals(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    exit_on: str,
    expected_status: int,
) -> None:
    layout = _layout(tmp_path)
    namespace = _module(layout)
    globals_ = namespace["main"].__globals__
    calls: list[tuple[int, signal.Signals]] = []
    closed: list[int] = []
    exited: set[int] = set()
    clock = [0.0]

    def pidfd_open(pid: int, flags: int) -> int:
        assert flags == 0
        return pid + 10_000

    def pidfd_signal(
        descriptor: int,
        selected_signal: signal.Signals,
        siginfo: object,
        flags: int,
    ) -> None:
        assert siginfo is None
        assert flags == 0
        calls.append((descriptor, selected_signal))
        if exit_on == "term" and selected_signal == signal.SIGTERM:
            exited.add(descriptor)
        if exit_on == "kill" and selected_signal == signal.SIGKILL:
            exited.add(descriptor)

    globals_.update(
        {
            "PIDFD_OPEN": pidfd_open,
            "PIDFD_SEND_SIGNAL": pidfd_signal,
            "PIDFD_WAIT": lambda descriptor: descriptor in exited,
            "CLOSE_FD": closed.append,
            "KILL": lambda _pid, _signal: pytest.fail("fallback signal used"),
            "MONOTONIC": lambda: clock[0],
            "SLEEP": lambda duration: clock.__setitem__(
                0,
                clock[0] + duration + 20,
            ),
        }
    )
    environment = dict(os.environ)
    try:
        assert namespace["main"](_arguments(layout)) == expected_status
    finally:
        os.environ.clear()
        os.environ.update(environment)
        _close_sockets(layout)

    assert {descriptor for descriptor, _signal in calls} == {11001, 11002, 11003}
    assert all(selected == signal.SIGTERM for _, selected in calls[:3])
    if exit_on == "term":
        assert len(calls) == 3
    else:
        assert [selected for _, selected in calls[3:]] == [signal.SIGKILL] * 3
    assert sorted(closed) == [11001, 11002, 11003]
    checkout = Path(layout["task_root"]) / "checkout"
    if expected_status == 0:
        assert not checkout.exists()
        assert json.loads(capsys.readouterr().out)["status"] == "cleaned"
    else:
        assert (checkout / "owned-evidence").exists()
        assert json.loads(capsys.readouterr().err)["issue"]["code"] == (
            "GATE_CLEANUP_PROCESS_STOP_FAILED"
        )


def test_mount_id_mismatch_fails_before_signals_or_deletion(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    layout = _layout(tmp_path)
    namespace = _module(layout)
    globals_ = namespace["main"].__globals__
    calls: list[tuple[int, signal.Signals]] = []

    def mount_id(descriptor: int) -> int:
        target = os.readlink(f"/proc/self/fd/{descriptor}")
        return 2 if target.endswith("/checkout") else 1

    globals_["MOUNT_ID"] = mount_id
    globals_["KILL"] = lambda pid, sig: calls.append((pid, sig))
    environment = dict(os.environ)
    try:
        assert namespace["main"](_arguments(layout)) == 65
    finally:
        os.environ.clear()
        os.environ.update(environment)
        _close_sockets(layout)

    assert json.loads(capsys.readouterr().err)["issue"]["code"] == (
        "GATE_CLEANUP_TARGET_MISMATCH"
    )
    assert calls == []
    assert (Path(layout["task_root"]) / "checkout/owned-evidence").exists()


def test_mount_id_is_revalidated_after_stop_before_any_deletion(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    layout = _layout(tmp_path)
    namespace = _module(layout)
    globals_ = namespace["main"].__globals__
    original = globals_["_preflight_deletion"]
    calls: list[tuple[int, signal.Signals]] = []

    def switch_mount_after_preflight(plan: object) -> int:
        mount_id = original(plan)
        globals_["MOUNT_ID"] = lambda _descriptor: 2
        return mount_id

    def stop(pid: int, selected_signal: signal.Signals) -> None:
        calls.append((pid, selected_signal))
        shutil.rmtree(Path(layout["proc_root"]) / str(pid))

    globals_["_preflight_deletion"] = switch_mount_after_preflight
    globals_["KILL"] = stop
    environment = dict(os.environ)
    try:
        assert namespace["main"](_arguments(layout)) == 69
    finally:
        os.environ.clear()
        os.environ.update(environment)
        _close_sockets(layout)

    assert json.loads(capsys.readouterr().err)["issue"]["code"] == (
        "GATE_CLEANUP_DELETE_FAILED"
    )
    assert calls == [
        (1003, signal.SIGTERM),
        (1002, signal.SIGTERM),
        (1001, signal.SIGTERM),
    ]
    assert (Path(layout["task_root"]) / "checkout/owned-evidence").exists()


@pytest.mark.parametrize("limit", ("depth", "nodes"))
def test_tree_limits_fail_in_preflight_before_signals_or_deletion(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    limit: str,
) -> None:
    layout = _layout(tmp_path)
    namespace = _module(layout)
    globals_ = namespace["main"].__globals__
    calls: list[tuple[int, signal.Signals]] = []
    checkout = Path(layout["task_root"]) / "checkout"
    if limit == "depth":
        _mkdir(checkout / "nested")
        globals_["MAXIMUM_DELETE_DEPTH"] = 1
    else:
        globals_["MAXIMUM_DELETE_NODES"] = 1
    globals_["KILL"] = lambda pid, sig: calls.append((pid, sig))
    environment = dict(os.environ)
    try:
        assert namespace["main"](_arguments(layout)) == 65
    finally:
        os.environ.clear()
        os.environ.update(environment)
        _close_sockets(layout)

    assert json.loads(capsys.readouterr().err)["issue"]["code"] == (
        "GATE_CLEANUP_LIMIT_EXCEEDED"
    )
    assert calls == []
    assert (checkout / "owned-evidence").exists()
