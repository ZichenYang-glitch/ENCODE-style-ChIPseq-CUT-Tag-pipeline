"""Bounded relaxation of the frozen-tree cleanup: an exit must not strand a run.

The F1 freeze proof treated every member as fatal while it could not be frozen.
A Nextflow submission burst makes a descendant exit inside that window routine
-- its stopped parent cannot reap it -- which used to abort ``kill_horse`` and
leave the run unacknowledged in ``running`` forever. These tests pin the
relaxation, its bounds, and the ordering evidence it depends on.

The failure branch itself is dynamic: a session-leader horse holding a real
unreaped zombie member, cleaned by the real ``kill_horse``.
"""

from __future__ import annotations

import logging
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import pytest

from encode_pipeline.services.process_runner import ProcessRunnerCleanupError
from encode_pipeline.workers import timeouts
from encode_pipeline.workers.timeouts import DurableWorker


def _stat(pid):
    try:
        value = Path(f"/proc/{pid}/stat").read_text()
    except FileNotFoundError:
        return None
    fields = value[value.rfind(")") + 2 :].split()
    return fields[0], int(fields[19])


def _wait_for(check, timeout=5.0):
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if check():
            return True
        time.sleep(0.005)
    return False


def _worker(pid):
    worker = object.__new__(DurableWorker)
    worker._horse_pid = pid
    worker.log = logging.getLogger("relaxed-owned-cleanup")
    return worker


def _finish(horse, member=None):
    if horse.poll() is None:
        try:
            os.killpg(horse.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    try:
        horse.wait(timeout=5)
    except subprocess.TimeoutExpired:  # pragma: no cover - fixture safety
        horse.kill()
        horse.wait(timeout=5)
    if member is not None:
        try:
            os.kill(member, signal.SIGKILL)
        except ProcessLookupError:
            pass
        try:
            os.waitpid(member, os.WNOHANG)
        except (ChildProcessError, OSError):
            pass


ZOMBIE_MEMBER_HORSE = (
    "import os,pathlib,sys,time\n"
    "pid=os.fork()\n"
    "if pid==0:\n"
    "    os._exit(0)\n"
    "pathlib.Path(sys.argv[1]).write_text(str(pid))\n"
    "time.sleep(30)\n"
)


def _start_zombie_horse(tmp_path):
    """A session-leader horse whose child exited and will never be reaped."""
    record = tmp_path / "zombie-member-pid"
    horse = subprocess.Popen(
        [sys.executable, "-I", "-B", "-c", ZOMBIE_MEMBER_HORSE, str(record)],
        start_new_session=True,
    )
    try:
        assert _wait_for(record.exists), "horse never recorded its member"
        assert horse.poll() is None
        member = int(record.read_text())
        assert _wait_for(lambda: (_stat(member) or (None,))[0] == "Z"), (
            "member never became a zombie"
        )
        identity = _stat(member)
        assert identity is not None
        assert os.getpgid(horse.pid) == horse.pid
        return horse, member, identity
    except BaseException:
        _finish(horse)
        raise


def test_exited_member_inside_freeze_window_still_confirms_cleanup(tmp_path, caplog):
    """The Gate-2 branch: a member exits, the run must still be acknowledged."""
    horse, member, identity = _start_zombie_horse(tmp_path)
    worker = _worker(horse.pid)
    caplog.set_level(logging.INFO, logger="relaxed-owned-cleanup")
    try:
        worker.kill_horse()

        assert worker._nested_cleanup_failed_pid is None
        assert worker._nested_cleanup_completed is not None
        report = worker._nested_cleanup_report
        assert report is not None
        assert report.horse_pid == horse.pid
        assert (memid := (member, identity[1])) in report.exited, report.exited
        # Registered identity set was complete for the root and every member.
        assert (horse.pid, report.horse_starttime) in report.registered
        assert memid in report.registered
        # The group signal still ran: RQ keeps its own horse wait status.
        assert horse.wait(timeout=5) == -signal.SIGKILL
        # Never silent: the relaxation is recorded in the worker log.
        messages = [record.getMessage() for record in caplog.records]
        assert any("retired" in message for message in messages), messages
    finally:
        _finish(horse, member)


def test_lost_root_stays_fatal_even_while_a_member_exits(tmp_path, monkeypatch):
    """Bounded: the retirement rule never covers the root."""
    horse, member, _identity = _start_zombie_horse(tmp_path)
    worker = _worker(horse.pid)
    original_signal = timeouts._signal_owned
    armed = False

    def kill_root_at_freeze(pid, starttime, sig):
        nonlocal armed
        original_signal(pid, starttime, sig)
        if pid == horse.pid and sig == signal.SIGSTOP and not armed:
            armed = True
            os.kill(horse.pid, signal.SIGKILL)
            assert _wait_for(lambda: (_stat(horse.pid) or (None,))[0] == "Z")

    monkeypatch.setattr(timeouts, "_signal_owned", kill_root_at_freeze)
    try:
        with pytest.raises(ProcessRunnerCleanupError):
            worker.kill_horse()
        assert armed
        assert worker._nested_cleanup_failed_pid == horse.pid
        assert getattr(worker, "_nested_cleanup_completed", None) is None
        assert getattr(worker, "_nested_cleanup_report", None) is None
    finally:
        _finish(horse, member)


def test_changed_member_identity_is_never_retired(tmp_path, monkeypatch):
    """Bounded: an identity that changed underneath us is still fatal."""
    record = tmp_path / "nested-member-pid"
    horse_source = (
        "import pathlib,subprocess,sys,time\n"
        "p=subprocess.Popen("
        "[sys.executable,'-I','-B','-c','import time;time.sleep(30)'],"
        "start_new_session=True)\n"
        f"pathlib.Path({str(record)!r}).write_text(str(p.pid))\n"
        "time.sleep(30)\n"
    )
    horse = subprocess.Popen(
        [sys.executable, "-I", "-B", "-c", horse_source], start_new_session=True
    )
    member = None
    original_read = timeouts._ProcessStat.read
    original_signal = timeouts._signal_owned
    armed = False
    try:
        assert _wait_for(record.exists)
        member = int(record.read_text())
        assert _wait_for(lambda: _stat(member) is not None)

        def read(pid):
            value = original_read(pid)
            if armed and str(pid) == str(member) and value is not None:
                return (value[0], value[1] + 1)
            return value

        def send(pid, starttime, sig):
            nonlocal armed
            original_signal(pid, starttime, sig)
            if pid == member and sig == signal.SIGSTOP:
                armed = True

        monkeypatch.setattr(timeouts._ProcessStat, "read", staticmethod(read))
        monkeypatch.setattr(timeouts, "_signal_owned", send)
        worker = _worker(horse.pid)
        with pytest.raises(ProcessRunnerCleanupError):
            worker.kill_horse()
        assert armed
        assert worker._nested_cleanup_failed_pid == horse.pid
        assert getattr(worker, "_nested_cleanup_report", None) is None
    finally:
        _finish(horse, member)


def test_member_visibility_separates_exit_from_identity_change(monkeypatch):
    """Only a truly absent or zombie entry may retire a registered member."""
    horse = subprocess.Popen(
        [sys.executable, "-I", "-B", "-c", "import time;time.sleep(30)"],
        start_new_session=True,
    )
    try:
        identity = _stat(horse.pid)
        assert identity is not None
        assert timeouts._member_visibility(horse.pid, identity[1]) == "live"
        assert timeouts._member_visibility(horse.pid, identity[1] + 1) == "changed"
        assert timeouts._member_visibility(999999, 1) == "gone"
        os.kill(horse.pid, signal.SIGSTOP)
        assert _wait_for(
            lambda: timeouts._member_visibility(horse.pid, identity[1]) == "frozen"
        )
        os.kill(horse.pid, signal.SIGCONT)
    finally:
        _finish(horse)


def test_uncovered_census_reports_foreign_session_child():
    """The detection hook for whatever the ownership proof cannot attribute."""
    child = subprocess.Popen(
        [sys.executable, "-I", "-B", "-c", "import time;time.sleep(30)"],
        start_new_session=True,
    )
    try:
        rows = [
            row
            for row in timeouts._uncovered_processes(horse_pid=1, registered=())
            if row["pid"] == child.pid
        ]
        assert rows, "a reparented foreign-session child must be reported"
        assert rows[0]["reason"] == "reparented live child outside the registered tree"
        assert rows[0]["pgrp"] == child.pid
    finally:
        child.terminate()
        child.wait(timeout=5)


def test_member_registration_precedes_the_group_signal(tmp_path, monkeypatch):
    """Every signalled identity is registered before the group kill."""
    record = tmp_path / "nested-member-pid"
    horse_source = (
        "import pathlib,subprocess,sys,time\n"
        "p=subprocess.Popen("
        "[sys.executable,'-I','-B','-c','import time;time.sleep(30)'],"
        "start_new_session=True)\n"
        f"pathlib.Path({str(record)!r}).write_text(str(p.pid))\n"
        "time.sleep(30)\n"
    )
    horse = subprocess.Popen(
        [sys.executable, "-I", "-B", "-c", horse_source], start_new_session=True
    )
    member = None
    original_descendants = timeouts._linux_descendants
    original_signal = timeouts._signal_owned
    original_killpg = os.killpg
    events = []
    try:
        assert _wait_for(record.exists)
        member = int(record.read_text())
        assert _wait_for(lambda: _stat(member) is not None)

        def descendants(pid):
            events.append(("descendants", pid))
            return original_descendants(pid)

        def send(pid, starttime, sig):
            original_signal(pid, starttime, sig)
            events.append(("signal", pid, int(sig)))

        def killpg(pgid, sig):
            events.append(("killpg", pgid, int(sig)))
            original_killpg(pgid, sig)

        monkeypatch.setattr(timeouts, "_linux_descendants", descendants)
        monkeypatch.setattr(timeouts, "_signal_owned", send)
        monkeypatch.setattr(timeouts.os, "killpg", killpg)

        worker = _worker(horse.pid)
        worker.kill_horse()
        report = worker._nested_cleanup_report
        assert report is not None

        group_signals = [i for i, event in enumerate(events) if event[0] == "killpg"]
        assert group_signals, events
        first_group_signal = group_signals[0]
        # No ownership discovery after the group signal.
        assert not any(
            event[0] == "descendants" for event in events[first_group_signal:]
        ), events
        # Every identity signalled with SIGKILL was already registered.
        registered_pids = {pid for pid, _starttime in report.registered}
        killed_pids = {
            event[1] for event in events if event[0] == "signal" and event[2] == 9
        }
        assert killed_pids
        assert killed_pids <= registered_pids, (killed_pids, registered_pids)
    finally:
        _finish(horse, member)
