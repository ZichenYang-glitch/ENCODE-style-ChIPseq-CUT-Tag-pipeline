"""Process-safe helpers for independent RQ worker integration tests."""

from __future__ import annotations

import os
from pathlib import Path
import signal
import subprocess
import sys


_WARM_SHUTDOWN_GRACE_SECONDS = 1.25
_FORCE_SHUTDOWN_GRACE_SECONDS = 5
_COMPLETED_PARENT_PIPE_GRACE_SECONDS = 0.1


def run_burst_worker(
    environment: dict[str, str],
    *,
    cwd: Path,
    timeout_seconds: float,
) -> subprocess.CompletedProcess[str]:
    """Run one burst worker with a bounded, RQ-aware shutdown path."""
    argv = [sys.executable, "-m", "encode_pipeline.workers.cli", "--burst"]
    process = subprocess.Popen(
        argv,
        cwd=cwd,
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    try:
        stdout, stderr = process.communicate(timeout=timeout_seconds)
    except subprocess.TimeoutExpired as exc:
        stdout, stderr = terminate_rq_worker(process)
        raise AssertionError(
            f"independent RQ worker timed out; stdout={stdout!r}; stderr={stderr!r}"
        ) from exc
    finally:
        if process.poll() is None:  # pragma: no cover - defensive test cleanup
            terminate_rq_worker(process)
    return subprocess.CompletedProcess(argv, process.returncode, stdout, stderr)


def terminate_rq_worker(process: subprocess.Popen[str]) -> tuple[str, str]:
    """Reap an RQ worker and its separately-grouped work horse.

    RQ moves its work horse into a new process group. The first SIGTERM asks
    the parent for a warm shutdown; a second signal after RQ's one-second
    guard requests a cold shutdown, which kills and waits for the horse group.
    Linux session process groups are discovered for a final defensive fallback.
    """
    session_id = process.pid
    if process.poll() is not None:
        try:
            return process.communicate(timeout=_COMPLETED_PARENT_PIPE_GRACE_SECONDS)
        except subprocess.TimeoutExpired:
            return _force_reap_worker_session(process, session_id)

    _signal_process(process, signal.SIGTERM)
    try:
        return process.communicate(timeout=_WARM_SHUTDOWN_GRACE_SECONDS)
    except subprocess.TimeoutExpired:
        _signal_process(process, signal.SIGTERM)

    try:
        return process.communicate(timeout=_FORCE_SHUTDOWN_GRACE_SECONDS)
    except subprocess.TimeoutExpired:
        return _force_reap_worker_session(process, session_id)


def _signal_process(process: subprocess.Popen[str], signum: signal.Signals) -> None:
    if process.poll() is not None:
        return
    try:
        process.send_signal(signum)
    except ProcessLookupError:
        pass


def _force_reap_worker_session(
    process: subprocess.Popen[str],
    session_id: int,
) -> tuple[str, str]:
    _kill_worker_session(session_id)
    if process.poll() is None:
        try:
            process.kill()
        except ProcessLookupError:
            pass
    try:
        return process.communicate(timeout=_FORCE_SHUTDOWN_GRACE_SECONDS)
    except subprocess.TimeoutExpired as exc:
        _kill_worker_session(session_id)
        raise AssertionError("RQ worker process session could not be reaped") from exc


def _kill_worker_session(session_id: int) -> None:
    for process_group in _session_process_groups(session_id):
        _kill_process_group(process_group)


def _session_process_groups(session_id: int) -> tuple[int, ...]:
    """Return Linux process groups still owned by the worker's session."""
    proc_root = Path("/proc")
    try:
        candidates = tuple(proc_root.iterdir())
    except OSError:
        return ()

    process_groups: set[int] = set()
    for candidate in candidates:
        if not candidate.name.isdigit():
            continue
        try:
            pid = int(candidate.name)
            if os.getsid(pid) != session_id:
                continue
            process_group = os.getpgid(pid)
        except (OSError, ValueError):
            continue
        if process_group > 0 and process_group != os.getpgrp():
            process_groups.add(process_group)
    return tuple(sorted(process_groups))


def _kill_process_group(process_group: int) -> None:
    if process_group <= 0 or process_group == os.getpgrp():
        return
    try:
        os.killpg(process_group, signal.SIGKILL)
    except ProcessLookupError:
        pass
