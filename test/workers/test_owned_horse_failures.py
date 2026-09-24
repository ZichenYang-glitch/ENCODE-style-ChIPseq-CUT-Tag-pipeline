"""Real RQ horse waits must not confirm cleanup from incomplete ownership proof."""

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
from rq import Worker

from encode_pipeline.services.process_runner import ProcessRunnerCleanupError
from encode_pipeline.workers import timeouts


def _state(pid):
    try:
        value = Path(f"/proc/{pid}/stat").read_text()
    except FileNotFoundError:
        return None
    fields = value.rsplit(")", 1)[1].split()
    return fields[0], int(fields[19])


def _wait_for(check):
    deadline = time.monotonic() + 5
    while not check():
        assert time.monotonic() < deadline
        time.sleep(0.005)


def _same(pid, identity):
    value = _state(pid)
    return value is not None and value[1] == identity[1]


def _exercise(tmp_path, monkeypatch, fault):
    """Only the indicated visibility/death point is injected; kill and wait4 run."""
    marker = tmp_path / "child-pid"
    horse_source = (
        "import pathlib,subprocess,sys,time\n"
        "p=subprocess.Popen([sys.executable,'-I','-B','-c','import time;time.sleep(30)'],start_new_session=True)\n"
        f"pathlib.Path({str(marker)!r}).write_text(str(p.pid))\n"
        "time.sleep(30)\n"
    )
    horse = subprocess.Popen(
        [sys.executable, "-I", "-B", "-c", horse_source], start_new_session=True
    )
    unrelated = subprocess.Popen(
        [sys.executable, "-I", "-B", "-c", "import time;time.sleep(30)"]
    )
    child = identity = None
    try:
        _wait_for(marker.exists)
        child = int(marker.read_text())
        identity = _state(child)
        assert identity is not None
        assert os.getpgid(child) == child
        worker = object.__new__(timeouts.DurableWorker)
        worker._horse_pid = horse.pid
        worker.log = logging.getLogger("owned-horse-failure")
        original_read = timeouts._ProcessStat.read
        original_signal = timeouts._signal_owned
        original_scandir = os.scandir
        original_text = Path.read_text
        real_rq_wait = Worker.wait_for_horse
        armed = False
        injected = 0
        rq_results = []

        def read(pid):
            nonlocal injected
            if fault == "unrelated-denied" and str(pid) == str(unrelated.pid):
                injected += 1
                raise PermissionError("unrelated proc entry is private")
            if (
                armed
                and str(pid) == str(child)
                and fault in {"denied", "malformed", "final-denied"}
            ):
                injected += 1
                if fault in {"denied", "final-denied"}:
                    raise PermissionError("controlled owned stat failure")
                return None
            return original_read(pid)

        def send(pid, starttime, sig):
            nonlocal armed, injected
            if (
                fault == "root-exit"
                and pid == horse.pid
                and sig == signal.SIGSTOP
                and not armed
            ):
                armed = True
                injected += 1
                os.kill(horse.pid, signal.SIGKILL)
                _wait_for(lambda: _state(horse.pid)[0] == "Z")
            original_signal(pid, starttime, sig)
            if pid == child and (
                (sig == signal.SIGSTOP and fault != "final-denied")
                or (sig == signal.SIGKILL and fault == "final-denied")
            ):
                armed = True

        def scan(path):
            nonlocal injected
            if fault == "scan-denied" and str(path) == f"/proc/{horse.pid}/task":
                injected += 1
                raise PermissionError("owned task enumeration denied")
            return original_scandir(path)

        def text(path, *args, **kwargs):
            nonlocal injected
            if fault in {"children-denied", "children-malformed"} and (
                str(path).startswith(f"/proc/{horse.pid}/task/")
                and path.name == "children"
            ):
                injected += 1
                if fault == "children-denied":
                    raise PermissionError("owned children list denied")
                return "not-a-pid"
            return original_text(path, *args, **kwargs)

        def wait(instance):
            value = real_rq_wait(instance)
            rq_results.append(value[:2])
            if value[0] == horse.pid:
                horse.returncode = os.waitstatus_to_exitcode(value[1])
            return value

        errors = {}
        with monkeypatch.context() as patch:
            patch.setattr(timeouts._ProcessStat, "read", staticmethod(read))
            patch.setattr(timeouts, "_signal_owned", send)
            patch.setattr(Worker, "wait_for_horse", wait)
            patch.setattr(os, "scandir", scan)
            patch.setattr(Path, "read_text", text)
            for name, operation in (
                ("kill", worker.kill_horse),
                ("wait", worker.wait_for_horse),
            ):
                try:
                    operation()
                except ProcessRunnerCleanupError:
                    errors[name] = "ProcessRunnerCleanupError"
        observation = {
            "fault": fault,
            "injected": injected,
            "errors": errors,
            "rq_results": rq_results,
            "horse_exit": horse.returncode,
            "child_identity": identity,
            "child_after": _state(child),
            "child_same_identity_alive": _same(child, identity),
            "failed_pid": getattr(worker, "_nested_cleanup_failed_pid", None),
            "unrelated_alive": unrelated.poll() is None,
        }
        (tmp_path / "observation.json").write_text(json.dumps(observation, indent=2))
        assert observation["unrelated_alive"]
        assert rq_results == [(horse.pid, signal.SIGKILL)]
        assert horse.returncode == -signal.SIGKILL
        if fault in {None, "unrelated-denied"}:
            assert errors == {}
            assert not observation["child_same_identity_alive"]
        else:
            assert injected > 0
            assert errors == {
                "kill": "ProcessRunnerCleanupError",
                "wait": "ProcessRunnerCleanupError",
            }
            assert observation["failed_pid"] == horse.pid
    finally:
        if child is not None and identity is not None and _same(child, identity):
            os.kill(child, signal.SIGKILL)
        if horse.poll() is None:
            horse.kill()
        horse.wait(timeout=5)
        if child is not None:
            try:
                os.waitpid(child, 0)
            except ChildProcessError:
                pass
        unrelated.terminate()
        unrelated.wait(timeout=5)
        if child is not None:
            _wait_for(lambda: _state(child) is None)
        (tmp_path / "cleanup.json").write_text(
            json.dumps(
                {
                    "horse_reaped": horse.poll() is not None,
                    "child_absent": child is None or _state(child) is None,
                    "unrelated_reaped": unrelated.poll() is not None,
                }
            )
        )


@pytest.mark.parametrize("fault", [None, "denied", "malformed", "root-exit"])
def test_rq_rejects_incomplete_owned_cleanup(tmp_path, monkeypatch, fault):
    _exercise(tmp_path, monkeypatch, fault)


@pytest.mark.parametrize(
    "fault",
    [
        "scan-denied",
        "children-denied",
        "children-malformed",
        "final-denied",
        "unrelated-denied",
    ],
)
def test_owned_discovery_and_final_visibility_boundaries(tmp_path, monkeypatch, fault):
    _exercise(tmp_path, monkeypatch, fault)


@pytest.mark.parametrize("visibility", ["missing", "replaced", "denied", "malformed"])
def test_signal_requires_confirmed_matching_identity(monkeypatch, visibility):
    sent = []

    def read(_pid):
        if visibility == "missing":
            raise FileNotFoundError
        if visibility == "replaced":
            return (1, 101)
        if visibility == "denied":
            raise PermissionError
        return None

    monkeypatch.setattr(timeouts._ProcessStat, "read", read)
    monkeypatch.setattr(os, "kill", lambda *args: sent.append(args))
    if visibility in {"denied", "malformed"}:
        with pytest.raises(ProcessRunnerCleanupError):
            timeouts._signal_owned(43210, 100, signal.SIGKILL)
    else:
        timeouts._signal_owned(43210, 100, signal.SIGKILL)
    assert sent == []


@pytest.mark.parametrize("visibility", ["missing", "denied", "malformed"])
def test_unknown_initial_root_never_signals_or_confirms(
    tmp_path, monkeypatch, visibility
):
    horse = subprocess.Popen(
        [sys.executable, "-I", "-B", "-c", "import time;time.sleep(30)"],
        start_new_session=True,
    )
    worker = object.__new__(timeouts.DurableWorker)
    worker._horse_pid = horse.pid
    worker.log = logging.getLogger("unknown-root-test")
    original = timeouts._ProcessStat.read
    try:

        def read(pid):
            if str(pid) == str(horse.pid):
                if visibility == "missing":
                    raise FileNotFoundError
                if visibility == "denied":
                    raise PermissionError
                return None
            return original(pid)

        with monkeypatch.context() as patch:
            patch.setattr(timeouts._ProcessStat, "read", read)
            with pytest.raises(ProcessRunnerCleanupError):
                worker.kill_horse()
        assert horse.poll() is None
        assert worker._nested_cleanup_failed_pid == horse.pid
        # Cleanup is the fixture's responsibility when identity is unknown.
        # The subsequent real RQ wait must still refuse acknowledgement.
        horse.kill()
        with pytest.raises(ProcessRunnerCleanupError):
            worker.wait_for_horse()
        horse.returncode = -signal.SIGKILL
    finally:
        if horse.poll() is None:
            horse.kill()
        horse.wait(timeout=5)
    (tmp_path / "cleanup.json").write_text(json.dumps({"horse_reaped": True}))


@pytest.mark.parametrize("lookup_error", [PermissionError, ProcessLookupError])
def test_unknown_group_lookup_refuses_later_rq_wait(
    tmp_path, monkeypatch, lookup_error
):
    """An unreadable or lost group is not evidence of detached-child cleanup."""
    import errno

    horse = subprocess.Popen(
        [sys.executable, "-I", "-B", "-c", "import time;time.sleep(30)"],
        start_new_session=True,
    )
    worker = object.__new__(timeouts.DurableWorker)
    worker._horse_pid = horse.pid
    worker.log = logging.getLogger("unknown-group-test")
    original = os.getpgid
    try:

        def getpgid(pid):
            if pid == horse.pid:
                code = errno.EPERM if lookup_error is PermissionError else errno.ESRCH
                raise lookup_error(code, "controlled group lookup failure")
            return original(pid)

        with monkeypatch.context() as patch:
            patch.setattr(os, "getpgid", getpgid)
            if lookup_error is PermissionError:
                with pytest.raises(PermissionError):
                    worker.kill_horse()
            else:
                assert worker.kill_horse() is None
        assert horse.poll() is None
        assert getattr(worker, "_nested_cleanup_failed_pid", None) == horse.pid
        horse.kill()
        with pytest.raises(ProcessRunnerCleanupError):
            worker.wait_for_horse()
        horse.returncode = -signal.SIGKILL
    finally:
        if horse.poll() is None:
            horse.kill()
        horse.wait(timeout=5)
    (tmp_path / "cleanup.json").write_text(json.dumps({"horse_reaped": True}))
