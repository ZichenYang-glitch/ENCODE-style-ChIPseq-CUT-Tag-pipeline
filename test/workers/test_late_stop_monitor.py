"""Original RQ monitor may consume a stop only after matching tree cleanup.

Real stop command, wait4, processes and monitor; Redis job lookup and terminal
callbacks are captured locally. The second-consumer arm enters original RQ
handle_job_failure and intercepts set_status before any Redis write.
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
from .test_horse_stop_wait_race import _fixture_subreaper, _state, _wait_until


class _StatusWriteReached(Exception):
    pass


def _exercise(tmp_path, monkeypatch, boundary, wrong_proof=None):
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
                    if wrong_proof is not None:
                        proof = [job.id, horse.pid, root_identity["starttime"]]
                        proof[{"job": 0, "pid": 1, "starttime": 2}[wrong_proof]] = (
                            "another-job" if wrong_proof == "job" else -1
                        )
                        worker._nested_cleanup_completed = tuple(proof)
                    release.set()
                monitor_thread.join(5)
                assert not monitor_thread.is_alive()
            current = _state(child)
            alive = (
                current is not None
                and current["starttime"] == child_identity["starttime"]
            )
            observed = {
                "boundary": boundary,
                "wrong_proof": wrong_proof,
                "callbacks": events,
                "errors": errors,
                "wait4": waits,
                "same_child_alive": alive,
                "completed": getattr(worker, "_nested_cleanup_completed", None),
                "sentinel_alive": sentinel.poll() is None,
            }
            (tmp_path / "observation.json").write_text(json.dumps(observed, indent=2))
            assert observed["sentinel_alive"]
            assert waits == [[horse.pid, signal.SIGKILL]]
            if boundary == "completed-stop":
                assert errors == {}
                assert events == ["stopped", "failure"]
                assert not alive
            else:
                assert alive  # A refusal is not automatic orphan recovery.
                assert "stopped" not in events
                assert "status:stopped" not in events
                assert errors.get("monitor") == "ProcessRunnerCleanupError"
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


@pytest.mark.parametrize("wrong_proof", ["job", "pid", "starttime"])
def test_monitor_rejects_unrelated_cleanup_proof(tmp_path, monkeypatch, wrong_proof):
    _exercise(tmp_path, monkeypatch, "after-wait", wrong_proof)
