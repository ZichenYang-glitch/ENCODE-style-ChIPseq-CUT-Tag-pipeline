"""Owned nested-session cleanup without scientific tools or Redis."""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

import pytest

from encode_pipeline.workers.timeouts import DurableWorker


def _stat(pid: int):
    try:
        value = Path(f"/proc/{pid}/stat").read_text()
    except FileNotFoundError:
        return None
    fields = value[value.rfind(")") + 2 :].split()
    return fields[0], int(fields[19])


def _worker(pid):
    worker = object.__new__(DurableWorker)
    worker._horse_pid = pid
    worker.log = logging.getLogger("nested-horse-test")
    return worker


def test_stop_reaps_nested_session_and_preserves_unrelated_child(tmp_path):
    """The worker must clean sessions which cannot receive its horse-group kill."""
    record = tmp_path / "owned.json"
    grandchild = "import time; time.sleep(30)"
    nested = (
        "import json,os,pathlib,subprocess,sys,time; "
        f"p=subprocess.Popen([sys.executable,'-I','-B','-c',{grandchild!r}]); "
        f"pathlib.Path({str(record)!r}).write_text(json.dumps([os.getpid(),p.pid])); "
        "time.sleep(30)"
    )
    horse_code = (
        "import subprocess,sys,time; "
        f"subprocess.Popen([sys.executable,'-I','-B','-c',{nested!r}],start_new_session=True); "
        "time.sleep(30)"
    )
    horse = subprocess.Popen(
        [sys.executable, "-I", "-B", "-c", horse_code], start_new_session=True
    )
    unrelated = subprocess.Popen([sys.executable, "-I", "-B", "-c", grandchild])
    owned = []
    try:
        deadline = time.monotonic() + 5
        while not record.exists():
            assert horse.poll() is None
            assert time.monotonic() < deadline
            time.sleep(0.01)
        owned = json.loads(record.read_text())
        identities = {pid: _stat(pid) for pid in owned}
        assert os.getpgid(owned[0]) == owned[0]
        assert owned[0] != horse.pid
        DurableWorker.kill_horse(_worker(horse.pid))
        # RQ still owns waiting/reaping its direct horse. The new cleanup must
        # not steal this wait status or rewrite the original stop signal.
        assert horse.wait(timeout=5) == -signal.SIGKILL
        assert unrelated.poll() is None
        assert all(_stat(pid) is None for pid in owned), identities
    finally:
        if owned:
            try:
                os.killpg(owned[0], signal.SIGKILL)
            except ProcessLookupError:
                pass
        if horse.poll() is None:
            os.killpg(horse.pid, signal.SIGKILL)
        horse.wait(timeout=5)
        unrelated.terminate()
        unrelated.wait(timeout=5)


def test_nested_stop_does_not_change_non_kill_signal(monkeypatch):
    worker = _worker(43210)
    seen = []
    monkeypatch.setattr(os, "getpgid", lambda _pid: 43210)
    monkeypatch.setattr(os, "killpg", lambda pid, sig: seen.append((pid, sig)))
    DurableWorker.kill_horse(worker, signal.SIGTERM)
    assert seen == [(43210, signal.SIGTERM)]


@pytest.mark.parametrize("cleanup_fails", [False, True])
def test_rq_wait_cannot_acknowledge_before_nested_cleanup(monkeypatch, cleanup_fails):
    """RQ waits in one thread while its pubsub thread kills the horse."""
    import threading

    from rq import Worker

    from encode_pipeline.services.process_runner import ProcessRunnerCleanupError
    import encode_pipeline.workers.timeouts as timeouts

    worker = _worker(43210)
    killed = threading.Event()
    release = threading.Event()
    waited = threading.Event()
    results = []
    errors = []
    monkeypatch.setattr(os, "getpgid", lambda _pid: 43210)
    monkeypatch.setattr(os, "killpg", lambda *_args: None)
    monkeypatch.setattr(timeouts._ProcessStat, "read", lambda _pid: (0, 100))

    def cleanup(_pid, _starttime, kill_horse):
        kill_horse()
        killed.set()
        assert release.wait(3)
        if cleanup_fails:
            raise ProcessRunnerCleanupError("injected cleanup refusal")

    monkeypatch.setattr(timeouts, "_kill_owned_horse_tree", cleanup)
    monkeypatch.setattr(Worker, "wait_for_horse", lambda _self: (43210, 9, None))

    def stop():
        try:
            worker.kill_horse()
        except ProcessRunnerCleanupError:
            errors.append("stop")

    def wait():
        try:
            results.append(worker.wait_for_horse())
        except ProcessRunnerCleanupError:
            errors.append("wait")
        finally:
            waited.set()

    stopping = threading.Thread(target=stop)
    waiting = threading.Thread(target=wait)
    stopping.start()
    try:
        assert killed.wait(3)
        waiting.start()
        assert not waited.wait(0.1), "RQ can acknowledge before cleanup finishes"
    finally:
        release.set()
        stopping.join(3)
        if waiting.ident is not None:
            waiting.join(3)
    assert not stopping.is_alive()
    assert not waiting.is_alive()
    if cleanup_fails:
        assert sorted(errors) == ["stop", "wait"]
        assert results == []
    else:
        assert errors == []
        assert results == [(43210, 9, None)]


def test_owned_signal_refuses_pid_reuse(monkeypatch):
    import encode_pipeline.workers.timeouts as timeouts

    sent = []
    monkeypatch.setattr(timeouts, "_linux_process_matches", lambda *_args: False)
    monkeypatch.setattr(os, "kill", lambda *args: sent.append(args))
    timeouts._signal_owned(43211, 999, signal.SIGKILL)
    assert sent == []


def _owned_pair(tmp_path):
    record = tmp_path / "nested-pid"
    child_code = "import time; time.sleep(30)"
    horse_code = (
        "import pathlib,subprocess,sys,time; "
        f"p=subprocess.Popen([sys.executable,'-I','-B','-c',{child_code!r}],start_new_session=True); "
        f"pathlib.Path({str(record)!r}).write_text(str(p.pid)); "
        "time.sleep(30)"
    )
    horse = subprocess.Popen(
        [sys.executable, "-I", "-B", "-c", horse_code], start_new_session=True
    )
    try:
        deadline = time.monotonic() + 5
        while not record.exists():
            assert horse.poll() is None
            assert time.monotonic() < deadline
            time.sleep(0.01)
        return horse, int(record.read_text())
    except BaseException:
        os.killpg(horse.pid, signal.SIGKILL)
        horse.wait(timeout=5)
        raise


def _finish_pair(horse, child):
    if _stat(child) is not None:
        try:
            os.killpg(child, signal.SIGKILL)
        except ProcessLookupError:
            pass
    if horse.poll() is None:
        os.killpg(horse.pid, signal.SIGKILL)
    horse.wait(timeout=5)


def test_freeze_failure_still_stops_owned_tree_and_refuses_confirmation(
    tmp_path, monkeypatch
):
    from encode_pipeline.services.process_runner import ProcessRunnerCleanupError
    import encode_pipeline.workers.timeouts as timeouts

    horse, child = _owned_pair(tmp_path)
    worker = _worker(horse.pid)

    def fail_state_read(*_args):
        raise PermissionError("injected /proc failure")

    monkeypatch.setattr(timeouts, "_stopped_or_gone", fail_state_read)
    try:
        with pytest.raises(ProcessRunnerCleanupError):
            worker.kill_horse()
        assert horse.wait(timeout=5) == -signal.SIGKILL
        assert _stat(child) is None
        assert worker._nested_cleanup_failed_pid == horse.pid
    finally:
        _finish_pair(horse, child)


def test_new_unrelated_child_during_subreaper_scope_is_not_adopted(
    tmp_path, monkeypatch
):
    import encode_pipeline.workers.timeouts as timeouts

    horse, child = _owned_pair(tmp_path)
    unrelated = []
    original_signal = timeouts._signal_owned

    def add_unrelated(pid, starttime, sig):
        original_signal(pid, starttime, sig)
        if pid == horse.pid and sig == signal.SIGSTOP:
            unrelated.append(
                subprocess.Popen(
                    [sys.executable, "-I", "-B", "-c", "import time; time.sleep(30)"]
                )
            )

    monkeypatch.setattr(timeouts, "_signal_owned", add_unrelated)
    try:
        _worker(horse.pid).kill_horse()
        assert horse.wait(timeout=5) == -signal.SIGKILL
        assert _stat(child) is None
        assert len(unrelated) == 1
        assert unrelated[0].poll() is None
    finally:
        _finish_pair(horse, child)
        for process in unrelated:
            process.terminate()
            process.wait(timeout=5)


def test_absent_horse_wait_preserves_rq_result(monkeypatch):
    from rq import Worker

    monkeypatch.setattr(Worker, "wait_for_horse", lambda _self: (None, None, None))
    assert _worker(0).wait_for_horse() == (None, None, None)


def test_stop_closes_session_spawn_window_after_initial_scan(tmp_path, monkeypatch):
    """A group created just after an empty scan must not escape the stop."""
    import encode_pipeline.workers.timeouts as timeouts

    trigger = tmp_path / "spawn"
    record = tmp_path / "spawned-pid"
    horse_code = (
        "import pathlib,subprocess,sys,time\n"
        f"trigger=pathlib.Path({str(trigger)!r})\n"
        "while not trigger.exists(): time.sleep(.001)\n"
        "p=subprocess.Popen([sys.executable,'-I','-B','-c','import time; time.sleep(30)'],start_new_session=True)\n"
        f"pathlib.Path({str(record)!r}).write_text(str(p.pid))\n"
        "time.sleep(30)\n"
    )
    horse = subprocess.Popen(
        [sys.executable, "-I", "-B", "-c", horse_code], start_new_session=True
    )
    scan = timeouts._linux_descendants
    armed = False

    def scan_then_enable_spawn(pid):
        nonlocal armed
        found = scan(pid)
        if pid == horse.pid and not armed:
            armed = True
            assert found == ()
            trigger.touch()
            # The old path permits the horse to start a new session after the
            # empty scan. A fixed path has already frozen it, so spawning cannot
            # advance and no arbitrary sleep is used to win a scheduling race.
            if _stat(horse.pid)[0] not in {"T", "t"}:
                deadline = time.monotonic() + 3
                while not record.exists():
                    assert time.monotonic() < deadline
                    time.sleep(0.005)
        return found

    monkeypatch.setattr(timeouts, "_linux_descendants", scan_then_enable_spawn)
    child = None
    try:
        _worker(horse.pid).kill_horse()
        assert horse.wait(timeout=5) == -signal.SIGKILL
        assert armed
        if record.exists():
            child = int(record.read_text())
            assert _stat(child) is None, "new session escaped after the empty scan"
    finally:
        if record.exists():
            child = int(record.read_text())
        if child is not None:
            _finish_pair(horse, child)
        elif horse.poll() is None:
            os.killpg(horse.pid, signal.SIGKILL)
            horse.wait(timeout=5)
