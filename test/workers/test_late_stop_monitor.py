"""Original RQ monitor may consume a stop only after matching tree cleanup.

Real stop command, wait4, processes and monitor; Redis job lookup and terminal
callbacks are captured locally. The second-consumer arm enters original RQ
handle_job_failure and intercepts set_status before any Redis write.

A stop whose proof does not complete is refused as data, not raised: the monitor
must reach RQ's abnormal-termination branch, settle on FAILED, and leave the
worker able to run its next job. ``test_gate_rejects_a_cleanup_proof_that_is_not_this_monitor``
pins the matching rule itself (job, pid and starttime) directly.
"""

from contextlib import nullcontext
import json
import logging
import os
import signal
import subprocess
import sys
import threading
from types import SimpleNamespace

import pytest
from rq import Worker
from rq.command import handle_stop_job_command
from rq.job import JobStatus

from encode_pipeline.workers import timeouts
from .test_horse_stop_wait_race import (
    _fixture_subreaper,
    _state,
    _wait_until,
    acknowledgement_decision,
)


class _StatusWriteReached(Exception):
    pass


def _exercise(tmp_path, monkeypatch, boundary):
    marker = tmp_path / "child"
    reached = threading.Event()
    release = threading.Event()
    events = []
    errors = {}
    waits = []
    child = child_identity = None
    monitor_thread = stop_thread = None
    with _fixture_subreaper():
        source = (
            "import pathlib,subprocess,sys,time\n"
            "p=subprocess.Popen([sys.executable,'-I','-B','-c',"
            "'import time;time.sleep(30)'],start_new_session=True)\n"
            f"pathlib.Path({str(marker)!r}).write_text(str(p.pid))\n"
            "time.sleep(30)\n"
        )
        horse = subprocess.Popen(
            [sys.executable, "-I", "-B", "-c", source], start_new_session=True
        )
        sentinel = subprocess.Popen(
            [sys.executable, "-I", "-B", "-c", "import time;time.sleep(30)"]
        )
        worker = object.__new__(timeouts.DurableWorker)
        worker._horse_pid = horse.pid
        worker._stopped_job_id = None
        worker.name = "late-stop-monitor"
        worker.log = logging.getLogger(worker.name)
        worker.job_monitoring_interval = 60
        worker.death_penalty_class = lambda *a, **kw: nullcontext()
        worker.get_current_job_id = lambda: "late-stop-job"
        worker.set_current_job_working_time = lambda value: None
        worker.connection = SimpleNamespace(pipeline=lambda: nullcontext(object()))
        real_wait = Worker.wait_for_horse

        def rq_wait(instance):
            result = real_wait(instance)
            waits.append(list(result[:2]))
            horse.returncode = os.waitstatus_to_exitcode(result[1])
            return result

        def barrier():
            reached.set()
            assert release.wait(5)

        def get_status():
            if boundary == "after-wait":
                barrier()
            return JobStatus.STARTED

        def killed(*args):
            events.append("killed")
            if boundary == "before-failure-read":
                barrier()

        def set_status(status, **kwargs):
            events.append("status:" + status.value)
            # This is the boundary immediately before a real Redis mutation.
            raise _StatusWriteReached

        job = SimpleNamespace(
            id="late-stop-job",
            stopped_callback=True,
            get_status=get_status,
            execute_stopped_callback=lambda *a: events.append("stopped"),
            should_retry=False,
            set_status=set_status,
        )
        worker.handle_work_horse_killed = killed
        if boundary == "before-failure-read":
            worker.handle_job_failure = lambda job, queue, **kw: (
                Worker.handle_job_failure(
                    worker, job, queue, started_job_registry=object(), **kw
                )
            )
        else:
            worker.handle_job_failure = lambda *a, **kw: events.append("failure")

        def monitor():
            try:
                worker.monitor_work_horse(job, None)
            except BaseException as exc:
                errors["monitor"] = type(exc).__name__

        def stop():
            try:
                handle_stop_job_command(worker, {"job_id": job.id})
            except BaseException as exc:
                errors["stop"] = type(exc).__name__

        try:
            _wait_until(marker.exists)
            child = int(marker.read_text())
            child_identity = _state(child)
            root_identity = _state(horse.pid)
            assert child_identity is not None and root_identity is not None
            with monkeypatch.context() as patch:
                patch.setattr(Worker, "wait_for_horse", rq_wait)
                if boundary == "completed-stop":
                    stop_thread = threading.Thread(target=stop)
                    stop_thread.start()
                    stop_thread.join(5)
                    assert not stop_thread.is_alive()
                else:
                    os.kill(horse.pid, signal.SIGKILL)
                    _wait_until(lambda: _state(horse.pid)["state"] == "Z")
                monitor_thread = threading.Thread(target=monitor)
                monitor_thread.start()
                if boundary != "completed-stop":
                    assert reached.wait(5)
                    # Original monitor has actually returned from real wait4.
                    assert waits == [[horse.pid, signal.SIGKILL]]
                    stop_thread = threading.Thread(target=stop)
                    stop_thread.start()
                    stop_thread.join(5)
                    assert not stop_thread.is_alive()
                    assert worker._stopped_job_id == job.id
                    release.set()
                monitor_thread.join(5)
                assert not monitor_thread.is_alive()
            current = _state(child)
            alive = (
                current is not None
                and current["starttime"] == child_identity["starttime"]
            )
            report = getattr(worker, "_nested_cleanup_report", None)
            observed = {
                "boundary": boundary,
                "callbacks": events,
                "errors": errors,
                "wait4": waits,
                "same_child_alive": alive,
                "report_confirmed": getattr(report, "confirmed", None),
                "unconfirmed_reason": getattr(report, "unconfirmed_reason", None),
                "cleanup_proof_absent": report is None,
                "stopped_marker_cleared": (
                    worker.__dict__.get("_nested_stopped_job_id") is None
                ),
                "cleanup_markers_released": not {
                    "_nested_cleanup_failed_pid",
                    "_nested_cleanup_completed",
                    "_nested_cleanup_unconfirmed",
                }
                & set(worker.__dict__),
                "sentinel_alive": sentinel.poll() is None,
            }
            (tmp_path / "observation.json").write_text(json.dumps(observed, indent=2))
            assert observed["sentinel_alive"]
            assert waits == [[horse.pid, signal.SIGKILL]]
            assert observed["stopped_marker_cleared"]
            assert observed["cleanup_markers_released"]
            if boundary == "completed-stop":
                assert errors == {}
                assert events == ["stopped", "failure"]
                assert not alive
                assert observed["report_confirmed"] is True
            else:
                assert alive  # A refusal is not automatic orphan recovery.
                assert "stopped" not in events
                assert "status:stopped" not in events
                # RQ already reaped the horse before the stop arrived, so there
                # was nothing left to clean and no proof to produce. The stop is
                # still refused -- as data, never as a raise: this is the frame
                # that used to unwind work() and end the worker, leaving the run
                # ``running`` forever.
                assert observed["cleanup_proof_absent"] is True
                if boundary == "before-failure-read":
                    # Original RQ handle_job_failure ran and chose FAILED; the
                    # stub raises at the boundary immediately before the write.
                    assert events == ["killed", "status:failed"]
                    assert errors == {"monitor": "_StatusWriteReached"}
                else:
                    assert events == ["killed", "failure"]
                    assert errors == {}
        finally:
            release.set()
            for thread in (monitor_thread, stop_thread):
                if thread is not None:
                    thread.join(5)
                    assert not thread.is_alive()
            if child is not None and child_identity is not None:
                current = _state(child)
                if (
                    current is not None
                    and current["starttime"] == child_identity["starttime"]
                ):
                    os.kill(child, signal.SIGKILL)
                try:
                    os.waitpid(child, 0)
                except ChildProcessError:
                    pass
            if horse.poll() is None:
                horse.kill()
            horse.wait(timeout=5)
            sentinel.terminate()
            sentinel.wait(timeout=5)
            cleanup = {
                "horse_absent": _state(horse.pid) is None,
                "child_absent": child is None or _state(child) is None,
                "sentinel_absent": _state(sentinel.pid) is None,
            }
            (tmp_path / "cleanup.json").write_text(json.dumps(cleanup, indent=2))
            assert all(cleanup.values())


@pytest.mark.parametrize(
    "boundary", ["after-wait", "before-failure-read", "completed-stop"]
)
def test_original_rq_monitor_requires_proven_cleanup_for_late_stop(
    tmp_path, monkeypatch, boundary
):
    _exercise(tmp_path, monkeypatch, boundary)


@pytest.mark.parametrize("wrong_field", ["job", "pid", "starttime"])
def test_gate_rejects_a_cleanup_proof_that_is_not_this_monitor(wrong_field):
    """A confirmed cleanup belonging to anything else never acknowledges.

    The proof is keyed to the monitored job *and* to the root identity that
    monitor observed, so a proof for another job, pid or starttime must not be
    promoted into an acknowledgement. The matching proof is asserted in the
    same case so the refusal cannot be a blanket "never acknowledge".
    """
    root_starttime = 100
    worker = object.__new__(timeouts.DurableWorker)
    worker._horse_pid = 43210
    worker.log = logging.getLogger("late-stop-gate")
    worker._stopped_job_id = "late-stop-job"
    wrong = ["late-stop-job", 43210, root_starttime]
    wrong[{"job": 0, "pid": 1, "starttime": 2}[wrong_field]] = (
        "another-job" if wrong_field == "job" else -1
    )
    worker._nested_cleanup_completed = tuple(wrong)
    worker._nested_cleanup_failed_pid = None
    assert (
        acknowledgement_decision(worker, "late-stop-job", 43210, root_starttime) is None
    )
    worker._nested_cleanup_completed = ("late-stop-job", 43210, root_starttime)
    assert (
        acknowledgement_decision(worker, "late-stop-job", 43210, root_starttime)
        == "late-stop-job"
    )
