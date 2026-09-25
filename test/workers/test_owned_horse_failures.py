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

from .test_horse_stop_wait_race import acknowledgement_decision


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
                except BaseException as exc:  # noqa: BLE001 - no leak may pass
                    errors[name] = type(exc).__name__
        report = getattr(worker, "_nested_cleanup_report", None)
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
            "cleanup_completed": getattr(worker, "_nested_cleanup_completed", None),
            "report_confirmed": getattr(report, "confirmed", None),
            "unconfirmed_reason": getattr(report, "unconfirmed_reason", None),
            "unrelated_alive": unrelated.poll() is None,
        }
        (tmp_path / "observation.json").write_text(json.dumps(observation, indent=2))
        assert observation["unrelated_alive"]
        assert rq_results == [(horse.pid, signal.SIGKILL)]
        assert horse.returncode == -signal.SIGKILL
        if fault in {None, "unrelated-denied"}:
            assert errors == {}
            assert not observation["child_same_identity_alive"]
            assert observation["report_confirmed"] is True
            assert observation["cleanup_completed"] is not None
        else:
            assert injected > 0
            # The proof failure is reported, not raised: neither the kill nor
            # RQ's own wait lets it escape, so the worker keeps running.
            assert errors == {}
            assert observation["failed_pid"] == horse.pid
            assert observation["cleanup_completed"] is None
            assert observation["report_confirmed"] is False
            assert observation["unconfirmed_reason"]
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


@pytest.mark.parametrize(
    ("visibility", "reason"),
    [
        ("missing", "OWNED_TREE_ROOT_LOST_BEFORE_FREEZE"),
        ("denied", "OWNED_PROCESS_IDENTITY_UNAVAILABLE"),
        ("malformed", "OWNED_PROCESS_IDENTITY_INVALID"),
    ],
)
def test_unknown_initial_root_never_signals_or_confirms(
    tmp_path, monkeypatch, visibility, reason
):
    """An unidentified root is refused as data and is never signalled.

    ``kill_horse`` must not raise -- RQ calls it from its pubsub thread, where a
    raise ends the thread and, because the worker's own exit status is discarded,
    the worker exits ``0`` and takes a fail-fast platform down with it. The
    refusal is recorded instead, and nothing may be signalled for an identity the
    cleanup never established.
    """
    horse = subprocess.Popen(
        [sys.executable, "-I", "-B", "-c", "import time;time.sleep(30)"],
        start_new_session=True,
    )
    worker = object.__new__(timeouts.DurableWorker)
    worker._horse_pid = horse.pid
    # What RQ's stop command writes before it calls kill_horse.
    worker._nested_stopped_job_id = "unknown-root-job"
    worker.log = logging.getLogger("unknown-root-test")
    original = timeouts._ProcessStat.read
    original_kill = os.kill
    original_killpg = os.killpg
    signals = []
    group_signals = []
    try:

        def read(pid):
            if str(pid) == str(horse.pid):
                if visibility == "missing":
                    raise FileNotFoundError
                if visibility == "denied":
                    raise PermissionError
                return None
            return original(pid)

        def kill(pid, sig):
            signals.append((pid, int(sig)))
            return original_kill(pid, sig)

        def killpg(pgid, sig):
            group_signals.append((pgid, int(sig)))
            return original_killpg(pgid, sig)

        with monkeypatch.context() as patch:
            patch.setattr(timeouts._ProcessStat, "read", read)
            patch.setattr(os, "kill", kill)
            patch.setattr(os, "killpg", killpg)
            assert worker.kill_horse() is None
        assert signals == []
        assert group_signals == []
        assert horse.poll() is None
        assert worker._nested_cleanup_failed_pid == horse.pid
        assert getattr(worker, "_nested_cleanup_completed", None) is None
        report = worker._nested_cleanup_report
        assert report is not None
        assert report.confirmed is False
        assert report.unconfirmed_reason == reason
        # No identity was ever established, so nothing is vouched for.
        assert report.horse_starttime is None
        # Fail closed: a stop whose root identity is unknown never acknowledges.
        assert (
            acknowledgement_decision(worker, "unknown-root-job", horse.pid, None)
            is None
        )
        (tmp_path / "observation.json").write_text(
            json.dumps(
                {
                    "visibility": visibility,
                    "signals": signals,
                    "group_signals": group_signals,
                    "reason": reason,
                },
                indent=2,
            )
        )
        # Cleanup is the fixture's responsibility when identity is unknown. The
        # real RQ wait still reaps the horse -- it is the gate that refuses, so
        # the run is reported failed rather than stopped. Signal the raw pid
        # rather than through Popen, whose send_signal polls and would reap.
        os.kill(horse.pid, signal.SIGKILL)
        assert worker.wait_for_horse()[0] == horse.pid
        horse.returncode = -signal.SIGKILL
    finally:
        if horse.poll() is None:
            horse.kill()
        horse.wait(timeout=5)
    (tmp_path / "cleanup.json").write_text(json.dumps({"horse_reaped": True}))


@pytest.mark.parametrize("lookup_error", [PermissionError, ProcessLookupError])
def test_unknown_group_never_signals_a_guessed_parent_group(
    tmp_path, monkeypatch, lookup_error
):
    """A group lookup that does not prove the horse's own group is not trusted.

    An unreadable group (EPERM) falls back to the identity-guarded signal for
    this one process, and the frozen-tree proof still decides the outcome. A lost
    group (ESRCH) is refused outright: ``os.getpgid`` failed, so nothing about the
    detached children was proven, and "the process is gone" is not evidence that
    its session subtree was cleaned. In neither arm may the worker's own
    pre-``setpgrp`` parent group be signalled.
    """
    import errno

    horse = subprocess.Popen(
        [sys.executable, "-I", "-B", "-c", "import time;time.sleep(30)"],
        start_new_session=True,
    )
    worker = object.__new__(timeouts.DurableWorker)
    worker._horse_pid = horse.pid
    worker._nested_stopped_job_id = "unknown-group-job"
    worker.log = logging.getLogger("unknown-group-test")
    original_getpgid = os.getpgid
    original_kill = os.kill
    original_killpg = os.killpg
    signals = []
    group_signals = []
    identity = _state(horse.pid)
    assert identity is not None
    try:

        def getpgid(pid):
            if pid == horse.pid:
                code = errno.EPERM if lookup_error is PermissionError else errno.ESRCH
                raise lookup_error(code, "controlled group lookup failure")
            return original_getpgid(pid)

        def kill(pid, sig):
            signals.append((pid, int(sig)))
            return original_kill(pid, sig)

        def killpg(pgid, sig):
            group_signals.append((pgid, int(sig)))
            return original_killpg(pgid, sig)

        with monkeypatch.context() as patch:
            patch.setattr(os, "getpgid", getpgid)
            patch.setattr(os, "kill", kill)
            patch.setattr(os, "killpg", killpg)
            assert worker.kill_horse() is None
        # The horse's group was never proven, so it is never signalled as one.
        assert group_signals == []
        if lookup_error is PermissionError:
            # Bounded fallback: only this identity-guarded process, then the
            # frozen-tree proof decides whether the cleanup is confirmed.
            assert (horse.pid, int(signal.SIGKILL)) in signals, signals
            assert worker._nested_cleanup_failed_pid is None
            assert worker._nested_cleanup_completed == (
                "unknown-group-job",
                horse.pid,
                identity[1],
            )
            report = worker._nested_cleanup_report
            assert report is not None and report.confirmed is True
            assert (
                acknowledgement_decision(
                    worker, "unknown-group-job", horse.pid, identity[1]
                )
                == "unknown-group-job"
            )
        else:
            # Lost group: refused, and nothing was signalled on a guess.
            assert signals == []
            assert horse.poll() is None
            assert worker._nested_cleanup_failed_pid == horse.pid
            assert getattr(worker, "_nested_cleanup_completed", None) is None
            assert getattr(worker, "_nested_cleanup_report", None) is None
            assert (
                acknowledgement_decision(
                    worker, "unknown-group-job", horse.pid, identity[1]
                )
                is None
            )
        (tmp_path / "observation.json").write_text(
            json.dumps(
                {
                    "lookup_error": lookup_error.__name__,
                    "signals": signals,
                    "group_signals": group_signals,
                    "failed_pid": getattr(worker, "_nested_cleanup_failed_pid", None),
                    "completed": getattr(worker, "_nested_cleanup_completed", None),
                    "horse_state": _state(horse.pid),
                },
                indent=2,
            )
        )
        # Signal the raw pid instead of going through Popen: send_signal polls,
        # which would reap the fallback-killed zombie before RQ's own wait4 and
        # hide whether the real wait still reports this horse.
        os.kill(horse.pid, signal.SIGKILL)
        assert worker.wait_for_horse()[0] == horse.pid
        horse.returncode = -signal.SIGKILL
    finally:
        if horse.poll() is None:
            horse.kill()
        horse.wait(timeout=5)
    (tmp_path / "cleanup.json").write_text(json.dumps({"horse_reaped": True}))
