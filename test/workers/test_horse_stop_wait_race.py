"""Real RQ stop/wait cannot acknowledge an unfinished tree-cleanup decision."""

from __future__ import annotations

from contextlib import contextmanager
import ctypes
import json
import logging
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time

import pytest
from rq import Worker
from rq.command import handle_stop_job_command

from encode_pipeline.workers import timeouts


def _state(pid):
    try:
        fields = Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()
    except FileNotFoundError:
        return None
    return {"state": fields[0], "starttime": int(fields[19])}


def _wait_until(check):
    deadline = time.monotonic() + 5
    while not check():
        assert time.monotonic() < deadline
        time.sleep(0.005)


def acknowledgement_decision(worker, job_id, horse_pid, starttime):
    """Read RQ's stopped-job gate from inside RQ's monitor scope.

    ``monitor_work_horse`` reads ``_stopped_job_id`` to choose between RQ's
    "stopped by user" branch and its unexpected-termination branch, and
    ``handle_job_failure`` reads it again to choose STOPPED versus FAILED. Both
    reads happen inside the scope ``DurableWorker.monitor_work_horse`` installs,
    so this helper reproduces exactly those two reads. Outside that scope the
    property deliberately reports the raw marker: the pubsub thread and
    diagnostic readers must never block behind a paused kill.
    """
    token = timeouts._STOP_MONITOR.set((worker, job_id, horse_pid, starttime))
    try:
        return worker._stopped_job_id
    finally:
        timeouts._STOP_MONITOR.reset(token)


@contextmanager
def _fixture_subreaper():
    # Only test cleanup needs to reap this fixture's orphan by known PID. Do not
    # hold the product's thread lock or use adopted children for its discovery.
    libc = ctypes.CDLL(None, use_errno=True)
    before = ctypes.c_int()
    assert libc.prctl(37, ctypes.byref(before), 0, 0, 0) == 0
    assert libc.prctl(36, 1, 0, 0, 0) == 0
    try:
        yield
    finally:
        assert libc.prctl(36, before.value, 0, 0, 0) == 0


def _exercise(tmp_path, monkeypatch, *, barrier, root_exits, wait_after_stop=False):
    marker = tmp_path / "child-pid"
    entered = threading.Event()
    release = threading.Event()
    reaped = threading.Event()
    waited = threading.Event()
    observation = {"barrier": barrier, "root_exits": root_exits}
    errors = {}
    results = []
    child = identity = None
    stopping = waiting = None
    with _fixture_subreaper():
        source = (
            "import pathlib,subprocess,sys,time\n"
            "p=subprocess.Popen([sys.executable,'-I','-B','-c',"
            "'import time;time.sleep(30)'],start_new_session=True)\n"
            f"pathlib.Path({str(marker)!r}).write_text(str(p.pid))\n"
            "time.sleep(30)\n"
        )
        horse = subprocess.Popen(
            [sys.executable, "-I", "-B", "-c", source],
            start_new_session=barrier != "pre-setpgrp-retry",
        )
        unrelated = subprocess.Popen(
            [sys.executable, "-I", "-B", "-c", "import time;time.sleep(30)"]
        )
        try:
            _wait_until(marker.exists)
            root_identity = _state(horse.pid)
            assert root_identity is not None
            child = int(marker.read_text())
            identity = _state(child)
            assert identity is not None
            worker = object.__new__(timeouts.DurableWorker)
            worker._horse_pid = horse.pid
            worker._stopped_job_id = None
            # Isolate the command's job-id lookup; its real flag assignment,
            # kill, wait4, process groups and /proc ownership checks all run.
            worker.get_current_job_id = lambda: "concurrent-stop-job"
            worker.log = logging.getLogger("concurrent-stop-test")
            original_getpgid = os.getpgid
            original_kill = worker.kill_horse
            real_rq_wait = Worker.wait_for_horse
            lookups = 0

            def pause():
                entered.set()
                assert release.wait(5)

            def getpgid(pid):
                nonlocal lookups
                value = original_getpgid(pid)
                if pid == horse.pid and threading.current_thread() is stopping:
                    lookups += 1
                    if (barrier == "first-lookup" and lookups == 1) or (
                        barrier == "pre-setpgrp-retry" and lookups == 2
                    ):
                        pause()
                return value

            def before_kill(*args, **kwargs):
                pause()
                return original_kill(*args, **kwargs)

            def rq_wait(instance):
                # Observe the actual wait4 result, never substitute a result.
                value = real_rq_wait(instance)
                results.append(list(value[:2]))
                horse.returncode = os.waitstatus_to_exitcode(value[1])
                reaped.set()
                return value

            def stop():
                try:
                    handle_stop_job_command(worker, {"job_id": "concurrent-stop-job"})
                except BaseException as exc:
                    errors["stop"] = type(exc).__name__

            def wait():
                try:
                    worker.wait_for_horse()
                except BaseException as exc:
                    errors["wait"] = type(exc).__name__
                finally:
                    waited.set()

            with monkeypatch.context() as patch:
                patch.setattr(os, "getpgid", getpgid)
                patch.setattr(Worker, "wait_for_horse", rq_wait)
                if barrier == "before-kill-entry":
                    patch.setattr(worker, "kill_horse", before_kill)
                stopping = threading.Thread(target=stop, name="real-rq-stop")
                stopping.start()
                assert entered.wait(5)
                assert worker._stopped_job_id == "concurrent-stop-job"
                if root_exits:
                    # Controlled death after the real stop command begins. The
                    # detached session child is intentionally not killed here.
                    os.kill(horse.pid, signal.SIGKILL)
                else:
                    release.set()
                if wait_after_stop:
                    stopping.join(5)
                    assert not stopping.is_alive()
                    assert _state(horse.pid)["state"] == "Z"
                    observation["horse_zombie_before_wait"] = True
                waiting = threading.Thread(target=wait, name="real-rq-wait")
                waiting.start()
                if root_exits:
                    # First synchronize actual wait4 completion. The bounded
                    # event observation tests acknowledgement, not kill timing.
                    assert reaped.wait(5)
                    observation["wait_finished_before_release"] = waited.wait(0.2)
                    observation["wait_error_before_release"] = errors.get("wait")
                    observation["child_before_release"] = _state(child)
                release.set()
                stopping.join(5)
                waiting.join(5)
                assert not stopping.is_alive() and not waiting.is_alive()
            current = _state(child)
            observation.update(
                errors=errors,
                rq_results=results,
                horse_pid=horse.pid,
                root_starttime=root_identity["starttime"],
                child_before=identity,
                child_after=current,
                same_child_alive=current is not None
                and current["starttime"] == identity["starttime"],
                failed_pid=getattr(worker, "_nested_cleanup_failed_pid", None),
                unrelated_alive=unrelated.poll() is None,
            )
            (tmp_path / "observation.json").write_text(
                json.dumps(observation, indent=2)
            )
            assert results == [[horse.pid, signal.SIGKILL]]
            assert observation["unrelated_alive"]
            if root_exits:
                # The refusal is data, not an exception: neither the command
                # thread nor the RQ wait observed a raise, and RQ's own
                # acknowledgement gate returns nothing for this job.
                assert errors == {}, errors
                assert observation["failed_pid"] == horse.pid
                assert (
                    acknowledgement_decision(
                        worker,
                        "concurrent-stop-job",
                        horse.pid,
                        root_identity["starttime"],
                    )
                    is None
                )
                # Fail closed does not imply automatic recovery of this orphan.
                assert observation["same_child_alive"]
            else:
                assert errors == {}
                assert not observation["same_child_alive"]
        finally:
            release.set()
            for thread in (stopping, waiting):
                if thread is not None:
                    thread.join(5)
                    assert not thread.is_alive()
            if child is not None and identity is not None:
                current = _state(child)
                if (
                    current is not None
                    and current["starttime"] == identity["starttime"]
                ):
                    os.kill(child, signal.SIGKILL)
                try:
                    os.waitpid(child, 0)
                except ChildProcessError:
                    pass
            if horse.poll() is None:
                horse.kill()
            horse.wait(timeout=5)
            unrelated.terminate()
            unrelated.wait(timeout=5)
            cleanup = {
                "horse_absent": _state(horse.pid) is None,
                "child_absent": child is None or _state(child) is None,
                "unrelated_reaped": _state(unrelated.pid) is None,
            }
            (tmp_path / "cleanup.json").write_text(json.dumps(cleanup, indent=2))
            assert all(cleanup.values())


@pytest.mark.parametrize("root_exits", [False, True])
@pytest.mark.parametrize(
    "barrier", ["first-lookup", "pre-setpgrp-retry", "before-kill-entry"]
)
def test_real_rq_stop_wait_covers_early_preparation(
    tmp_path, monkeypatch, barrier, root_exits
):
    _exercise(tmp_path, monkeypatch, barrier=barrier, root_exits=root_exits)


def test_completed_real_stop_can_be_waited_after_horse_becomes_zombie(
    tmp_path, monkeypatch
):
    _exercise(
        tmp_path,
        monkeypatch,
        barrier="first-lookup",
        root_exits=False,
        wait_after_stop=True,
    )


@pytest.mark.parametrize("signum", [None, signal.SIGTERM])
def test_normal_wait_and_non_sigkill_preserve_real_rq_status(tmp_path, signum):
    source = "raise SystemExit(7)" if signum is None else "import time;time.sleep(30)"
    horse = subprocess.Popen(
        [sys.executable, "-I", "-B", "-c", source], start_new_session=True
    )
    worker = object.__new__(timeouts.DurableWorker)
    worker._horse_pid = horse.pid
    worker._stopped_job_id = None
    worker.log = logging.getLogger("normal-rq-wait")
    try:
        if signum is not None:
            worker.kill_horse(signum)
        result = worker.wait_for_horse()
        horse.returncode = os.waitstatus_to_exitcode(result[1])
        assert result[0] == horse.pid
        assert horse.returncode == (7 if signum is None else -signum)
        assert result[2] is not None
        assert getattr(worker, "_nested_cleanup_failed_pid", None) is None
        (tmp_path / "wait.json").write_text(
            json.dumps({"pid_status": list(result[:2]), "exit": horse.returncode})
        )
    finally:
        if horse.poll() is None:
            horse.kill()
        horse.wait(timeout=5)
        assert _state(horse.pid) is None
