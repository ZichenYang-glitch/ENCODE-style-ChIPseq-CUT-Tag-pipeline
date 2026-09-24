"""POSIX process-group lifecycle for one private qualification attempt."""

from __future__ import annotations

import ctypes
from dataclasses import dataclass
import os
from pathlib import Path
import signal
import subprocess
import time
from collections.abc import Callable


@dataclass(frozen=True)
class Execution:
    returncode: int
    reason: str | None
    pid: int
    elapsed_seconds: float


def _reap_group(group: int) -> None:
    while True:
        try:
            if os.waitpid(-group, os.WNOHANG)[0] == 0:
                break
        except ChildProcessError:
            break


def _exists(group: int) -> bool:
    try:
        os.killpg(group, 0)
        return True
    except ProcessLookupError:
        return False


def execute(
    argv: list[str],
    cwd: Path,
    env: dict[str, str],
    private: Path,
    timeout: float,
    cancelled: Callable[[], bool],
) -> Execution:
    # Qualification is Linux-only. Subreaping is scoped to this call and only
    # this child's process group is waited/killed, never unrelated children.
    libc = ctypes.CDLL(None, use_errno=True)
    previous = ctypes.c_int()
    if libc.prctl(37, ctypes.byref(previous), 0, 0, 0) != 0:
        raise OSError("subreaper unavailable")
    if libc.prctl(36, 1, 0, 0, 0) != 0:
        raise OSError("subreaper unavailable")
    process = None
    start = time.monotonic()
    reason = None
    try:
        with (
            (private / "upstream.stdout").open("xb") as stdout,
            (private / "upstream.stderr").open("xb") as stderr,
        ):
            process = subprocess.Popen(
                argv,
                cwd=cwd,
                env=env,
                stdin=subprocess.DEVNULL,
                stdout=stdout,
                stderr=stderr,
                start_new_session=True,
            )
            while process.poll() is None:
                if cancelled():
                    reason = "cancelled"
                    break
                if time.monotonic() - start >= timeout:
                    reason = "timed_out"
                    break
                time.sleep(0.025)
            if reason is None:
                process.wait()
                # Finished parents must not leave scientific children running.
                for _ in range(10):
                    _reap_group(process.pid)
                    if not _exists(process.pid):
                        break
                    time.sleep(0.025)
                if _exists(process.pid):
                    reason = "child_survived_parent"
    finally:
        try:
            if process is not None:
                if _exists(process.pid):
                    try:
                        os.killpg(process.pid, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                    end = time.monotonic() + 1.0
                    while time.monotonic() < end:
                        # Popen must reap its direct child before group waitpid.
                        if process.poll() is not None:
                            _reap_group(process.pid)
                        if not _exists(process.pid):
                            break
                        time.sleep(0.025)
                    if _exists(process.pid):
                        try:
                            os.killpg(process.pid, signal.SIGKILL)
                        except ProcessLookupError:
                            pass
                process.wait()
                end = time.monotonic() + 2.0
                while time.monotonic() < end:
                    _reap_group(process.pid)
                    if not _exists(process.pid):
                        break
                    time.sleep(0.025)
                if _exists(process.pid):
                    raise RuntimeError("attempt_process_cleanup_incomplete")
        finally:
            if libc.prctl(36, previous.value, 0, 0, 0) != 0:
                raise OSError("subreaper_restore_failed")
    return Execution(process.returncode, reason, process.pid, time.monotonic() - start)
