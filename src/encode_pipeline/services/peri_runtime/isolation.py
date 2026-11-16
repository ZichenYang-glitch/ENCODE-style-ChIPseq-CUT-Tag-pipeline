"""Bubblewrap launch for one Peri process. There is no plain-subprocess fallback."""

from __future__ import annotations

from collections.abc import Callable
import json
import os
from pathlib import Path
import queue
import subprocess
import threading
from uuid import uuid4

from encode_pipeline.services.peri_runtime.limits import (
    qualified_native_pid,
    stop_scope,
)
from encode_pipeline.services.peri_runtime.policy import process_policy


Popen = subprocess.Popen
BWRAP = Path("/usr/bin/bwrap")
PRLIMIT = Path("/usr/bin/prlimit")
SYSTEMD_RUN = Path("/usr/bin/systemd-run")
PYTHON = Path("/usr/bin/python3.12")
PYTHON_LIB = Path("/usr/lib/python3.12")
LIB_DIR = Path("/lib/x86_64-linux-gnu")
LIB64 = Path("/lib64")
CA_BUNDLE = Path("/etc/ssl/certs/ca-certificates.crt")
RELAY_SCRIPT = Path(__file__).with_name("relay.py")


def isolation_paths_exist() -> bool:
    return all(
        path.exists()
        for path in (
            BWRAP,
            PRLIMIT,
            SYSTEMD_RUN,
            PYTHON,
            PYTHON_LIB,
            LIB_DIR,
            LIB64,
            CA_BUNDLE,
            RELAY_SCRIPT,
        )
    )


def command(
    binary: Path,
    state: Path,
    model_socket: Path,
    *,
    policy: dict[str, object] | None = None,
    prediction: str = "off",
) -> list[str]:
    args = [
        str(BWRAP),
        "--unshare-all",
        "--unshare-user",
        "--disable-userns",
        "--die-with-parent",
        "--new-session",
        "--cap-drop",
        "ALL",
        "--clearenv",
        "--proc",
        "/proc",
        "--dev",
        "/dev",
        "--tmpfs",
        "/tmp",
        "--dir",
        "/work",
        "--dir",
        "/runtime",
        "--ro-bind",
        str(binary.resolve()),
        "/runtime/peri",
        "--ro-bind",
        str(RELAY_SCRIPT.resolve()),
        "/runtime/relay.py",
        "--ro-bind",
        str(PYTHON),
        "/usr/bin/python3.12",
        "--ro-bind",
        str(PYTHON_LIB),
        "/usr/lib/python3.12",
        "--ro-bind",
        str(LIB_DIR),
        "/lib/x86_64-linux-gnu",
        "--ro-bind",
        str(LIB64),
        "/lib64",
        "--ro-bind",
        str(CA_BUNDLE),
        "/etc/ssl/certs/ca-certificates.crt",
        "--bind",
        str(state.resolve()),
        "/state",
        "--ro-bind",
        str(model_socket.resolve()),
        "/model.sock",
        "--setenv",
        "PATH",
        "/usr/bin",
        "--setenv",
        "HOME",
        "/state/home",
        "--setenv",
        "XDG_CONFIG_HOME",
        "/state/config",
        "--setenv",
        "XDG_CACHE_HOME",
        "/state/cache",
        "--setenv",
        "XDG_DATA_HOME",
        "/state/data",
        "--setenv",
        "TMPDIR",
        "/tmp",
        "--setenv",
        "LANG",
        "C.UTF-8",
        "--setenv",
        "TZ",
        "UTC",
        "--setenv",
        "PERI_MCP_BUILTIN",
        "off",
        "--setenv",
        "SSL_CERT_FILE",
        "/etc/ssl/certs/ca-certificates.crt",
        "--chdir",
        "/work",
        "--setenv",
        "PERI_PROCESS_TOOL_POLICY",
        json.dumps(process_policy() if policy is None else policy),
        "--setenv",
        "PERI_PREDICTION_POLICY",
        prediction,
        "--setenv",
        "TOKIO_WORKER_THREADS",
        "2",
        "--remount-ro",
        "/",
        "--",
        "/usr/bin/python3.12",
        "-I",
        "-S",
        "/runtime/relay.py",
    ]
    return args


class _LaunchTicket:
    """One creation. A process has one owner: the caller or the launcher."""

    def __init__(
        self,
        start: Callable[[], subprocess.Popen[bytes]],
        scope_unit: str | None,
    ) -> None:
        self.start = start
        self.scope_unit = scope_unit
        self._lock = threading.Lock()
        self.phase = "queued"
        self.abandoned = False
        self.box: queue.Queue[object] = queue.Queue(maxsize=1)

    def claim_to_run(self) -> bool:
        with self._lock:
            if self.abandoned or self.phase != "queued":
                return False
            self.phase = "running"
            return True

    def settle(self, outcome: object) -> bool:
        """Return whether the caller still owns a successful process."""

        with self._lock:
            self.phase = "settled"
            if self.abandoned:
                return False
            self.box.put(outcome)
            return True

    def give_up(self) -> str:
        with self._lock:
            if self.phase == "settled" and not self.abandoned:
                return "delivered"
            if self.abandoned:
                return "abandoned"
            self.abandoned = True
            phase = self.phase
            try:
                self.box.put_nowait(
                    OSError("native launcher did not create the process")
                )
            except queue.Full:
                pass
            return phase


def _wait_for_launch(ticket: _LaunchTicket, timeout: float) -> object:
    try:
        return ticket.box.get(timeout=timeout)
    except queue.Empty:
        if ticket.give_up() != "delivered":
            raise OSError("native launcher did not create the process") from None
        try:
            return ticket.box.get_nowait()
        except queue.Empty:
            raise OSError("native launcher did not create the process") from None


def _reap_unclaimed(process: subprocess.Popen[bytes], scope_unit: str | None) -> None:
    """Close a process the caller timed out or shut down before receiving."""

    stop_scope(scope_unit)
    if process.poll() is None:
        try:
            process.kill()
        except OSError:
            pass
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        pass
    for stream in (process.stdin, process.stdout):
        if stream is None:
            continue
        try:
            stream.close()
        except OSError:
            pass


class NativeLauncher:
    """Create native processes on one thread that outlives each request.

    bwrap --die-with-parent uses PR_SET_PDEATHSIG. Linux delivers that signal
    when the creating thread exits, even if the owner process is still alive.
    Request handlers are short-lived, so they must not be that parent.
    A caller that stops waiting does not leave the creation running: a queued
    request is dropped, and a process already being created is reaped here.
    """

    def __init__(self) -> None:
        self._queue: queue.Queue[_LaunchTicket | None] = queue.Queue()
        self._pending: set[_LaunchTicket] = set()
        self._pending_lock = threading.Lock()
        self._thread = threading.Thread(
            target=self._serve,
            name="peri-native-launcher",
            daemon=True,
        )
        self._thread.start()

    def _forget(self, ticket: _LaunchTicket) -> None:
        with self._pending_lock:
            self._pending.discard(ticket)

    def _serve(self) -> None:
        while True:
            ticket = self._queue.get()
            if ticket is None:
                return
            if not ticket.claim_to_run():
                self._forget(ticket)
                continue
            try:
                outcome = ticket.start()
            except BaseException as error:
                ticket.settle(error)
                self._forget(ticket)
                continue
            if ticket.settle(outcome):
                self._forget(ticket)
                continue
            if isinstance(outcome, Popen):
                _reap_unclaimed(outcome, ticket.scope_unit)
            self._forget(ticket)

    def popen(
        self,
        args: list[str],
        *,
        scope_unit: str | None = None,
        timeout: float = 30,
        **kwargs: object,
    ) -> subprocess.Popen[bytes]:
        if threading.get_ident() == self._thread.ident:
            return subprocess.Popen(args, **kwargs)  # type: ignore[arg-type]
        ticket = _LaunchTicket(
            lambda: subprocess.Popen(args, **kwargs),  # type: ignore[arg-type]
            scope_unit,
        )
        with self._pending_lock:
            self._pending.add(ticket)
        self._queue.put(ticket)
        outcome = _wait_for_launch(ticket, timeout)
        self._forget(ticket)
        if isinstance(outcome, BaseException):
            raise outcome
        return outcome  # type: ignore[return-value]

    def abandon_pending(self) -> None:
        """Drop queued creations and mark in-flight ones for launcher reap."""

        with self._pending_lock:
            pending = tuple(self._pending)
        for ticket in pending:
            ticket.give_up()

    def alive(self) -> bool:
        return self._thread.is_alive()

    def native_id(self) -> int | None:
        return self._thread.native_id


_native_launcher = NativeLauncher()


def native_launcher() -> NativeLauncher:
    return _native_launcher


def spawn(
    binary: Path,
    state: Path,
    model_socket: Path,
    stderr,
    *,
    policy: dict[str, object] | None = None,
    prediction: str = "off",
) -> subprocess.Popen[bytes]:
    for name in ("home", "config", "cache", "data"):
        (state / name).mkdir(parents=True, exist_ok=True)
    unit = f"hw-peri-a1-{uuid4().hex[:12]}"
    scope_env = {
        key: os.environ[key]
        for key in ("DBUS_SESSION_BUS_ADDRESS", "XDG_RUNTIME_DIR")
        if os.environ.get(key)
    }
    if len(scope_env) != 2:
        raise OSError("systemd user session is required for the native scope")
    process = _native_launcher.popen(
        [
            str(SYSTEMD_RUN),
            "--user",
            "--scope",
            "--collect",
            "--quiet",
            f"--unit={unit}",
            "--property=MemoryMax=512M",
            "--property=MemorySwapMax=0",
            "--property=CPUQuota=100%",
            "--property=TasksMax=128",
            "--",
            str(PRLIMIT),
            "--nofile=128:128",
            "--core=0:0",
            "--fsize=16777216:16777216",
            "--",
            *command(binary, state, model_socket, policy=policy, prediction=prediction),
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=stderr,
        close_fds=True,
        scope_unit=unit,
        env=scope_env,
        start_new_session=True,
    )
    try:
        native_pid = qualified_native_pid(unit, process.pid)
    except OSError:
        stop_scope(unit)
        if process.poll() is None:
            process.kill()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                pass
        raise
    process.peri_scope_unit = unit  # type: ignore[attr-defined]
    process.peri_native_pid = native_pid  # type: ignore[attr-defined]
    return process
