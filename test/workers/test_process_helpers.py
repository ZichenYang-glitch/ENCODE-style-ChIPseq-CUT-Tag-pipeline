"""Regression tests for independent worker process cleanup."""

from __future__ import annotations

import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import pytest

from .process_helpers import terminate_rq_worker


@pytest.mark.skipif(
    sys.platform != "linux",
    reason="orphaned RQ work-horse fallback uses Linux /proc session discovery",
)
def test_terminate_rq_worker_reaps_orphaned_child_process_group(tmp_path: Path):
    """A dead parent cannot leave a separately-grouped pipe holder alive."""
    marker = tmp_path / "child-pid"
    script = """
import os
from pathlib import Path
import sys
import time

child_pid = os.fork()
if child_pid == 0:
    os.setpgrp()
    Path(sys.argv[1]).write_text(str(os.getpid()), encoding="utf-8")
    time.sleep(60)
    os._exit(0)
os._exit(0)
"""
    process = subprocess.Popen(
        [sys.executable, "-c", script, str(marker)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    child_pid: int | None = None
    try:
        child_pid = _wait_for_child_pid(marker)
        assert process.wait(timeout=2) == 0

        started_at = time.monotonic()
        stdout, stderr = terminate_rq_worker(process)
        elapsed = time.monotonic() - started_at

        assert elapsed < 2
        assert stdout == ""
        assert stderr == ""
        _wait_for_process_group_exit(child_pid)
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=2)
        if child_pid is not None:
            try:
                os.killpg(child_pid, signal.SIGKILL)
            except ProcessLookupError:
                pass


def _wait_for_child_pid(marker: Path) -> int:
    deadline = time.monotonic() + 2
    while time.monotonic() < deadline:
        try:
            contents = marker.read_text(encoding="utf-8").strip()
            if contents:
                return int(contents)
        except (FileNotFoundError, ValueError):
            pass
        time.sleep(0.01)
    raise AssertionError("worker child pid marker was not completed")


def _wait_for_process_group_exit(process_group: int) -> None:
    """Assert that a process group disappears within a short reap window."""
    deadline = time.monotonic() + 2
    while time.monotonic() < deadline:
        try:
            os.killpg(process_group, 0)
        except ProcessLookupError:
            return
        time.sleep(0.01)
    raise AssertionError("worker child process group was not reaped")
