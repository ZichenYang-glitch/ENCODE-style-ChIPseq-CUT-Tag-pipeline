"""The boundary that turns a failed ownership proof into data.

``DurableWorker.kill_horse`` can fail its ownership proof in many places (a
denied ``/proc`` read, a lost or changed member, an unstabilized tree). Every
one of those failures used to be *raised*, and a raise had two routes out of the
worker: RQ's pubsub thread (``handle_payload`` -> ``handle_command`` ->
``kill_horse``, where the thread's exception handler re-raises anything that is
not a Redis connection error) and RQ's own monitor/deadline path (where it
unwinds ``work()``'s bare ``except``). Either one ends the worker -- and because
``workers/cli.py`` discards ``work()``'s return value, the worker process still
exits ``0``. A fail-fast supervisor then reads a clean exit, tears the whole
platform down around it, and every later request sees a closed port.

These tests pin the boundary that stops the escape, the fail-closed behaviour it
must not weaken, and the fixed codes that make a refusal explainable.
"""

from __future__ import annotations

import ast
import errno
import json
from pathlib import Path
import signal

import pytest
from redis.exceptions import ConnectionError as RedisConnectionError

from encode_pipeline.services.process_runner import ProcessRunnerCleanupError
from encode_pipeline.workers import timeouts

from .test_horse_stop_wait_race import acknowledgement_decision


ROOT_STARTTIME = 500
HORSE_PID = 43210


def _worker():
    import logging

    worker = object.__new__(timeouts.DurableWorker)
    worker.name = "amplifier-boundary"
    worker.log = logging.getLogger("amplifier-boundary")
    worker._horse_pid = HORSE_PID
    worker.get_current_job_id = lambda: "boundary-job"
    return worker


def _stop_payload(job_id: str) -> dict:
    return {"data": json.dumps({"command": "stop-job", "job_id": job_id}).encode()}


def _cleanup_errors_in_module() -> set[str]:
    """Every message this module can pass to ``ProcessRunnerCleanupError``."""
    tree = ast.parse(Path(timeouts.__file__).read_text())
    messages = set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        function = node.func
        if not isinstance(function, ast.Name) or function.id != (
            "ProcessRunnerCleanupError"
        ):
            continue
        assert len(node.args) == 1 and not node.keywords, ast.dump(node)
        message = node.args[0]
        # A non-literal message would silently degrade to the fallback code.
        assert isinstance(message, ast.Constant) and isinstance(message.value, str), (
            ast.dump(node)
        )
        messages.add(message.value)
    return messages


def test_every_cleanup_failure_message_maps_to_a_fixed_code():
    """The code table must cover the module exactly, in both directions.

    ``_unconfirmed_reason_code`` relays no free-form message: an unmapped message
    degrades to the fallback code, which would erase the distinction between
    failure modes in the worker log. A new raise site must therefore add a code.
    """
    messages = _cleanup_errors_in_module()
    assert messages
    assert messages == set(timeouts._OWNED_CLEANUP_UNCONFIRMED_CODES)
    for message, code in timeouts._OWNED_CLEANUP_UNCONFIRMED_CODES.items():
        assert code and code.isupper()
        error = ProcessRunnerCleanupError(message)
        assert str(error) == message
        assert timeouts._unconfirmed_reason_code(error) == code
    denied = PermissionError(errno.EPERM, "controlled denial")
    assert timeouts._unconfirmed_reason_code(denied) == "OWNED_TREE_SIGNAL_DENIED"
    assert (
        timeouts._unconfirmed_reason_code(RuntimeError("unmapped"))
        == timeouts._OWNED_CLEANUP_UNCONFIRMED_FALLBACK
    )


@pytest.mark.parametrize("cleanup_fails", [False, True])
def test_stop_command_cleanup_failure_never_escapes_handle_payload(
    monkeypatch, cleanup_fails
):
    """A real ``stop-job`` command keeps its refusal inside the worker.

    The command thread must survive a failed proof: RQ's handler re-raises
    everything that is not a Redis connection error, and the worker's own exit
    status is discarded, so an escape reads as a clean shutdown to a fail-fast
    supervisor. The refusal is recorded as data and the stop stays unacknowledged.
    """
    worker = _worker()
    kills = []

    def cleanup(self, sig=signal.SIGKILL):
        kills.append((self.horse_pid, int(sig)))
        if cleanup_fails:
            self._nested_cleanup_failed_pid = self.horse_pid
            self._nested_cleanup_completed = None
            raise ProcessRunnerCleanupError(
                "Nested workflow cleanup could not be confirmed."
            )
        self._nested_cleanup_failed_pid = None
        self._nested_cleanup_completed = (
            self.__dict__.get("_nested_stopped_job_id"),
            self.horse_pid,
            ROOT_STARTTIME,
        )

    monkeypatch.setattr(timeouts.DurableWorker, "kill_horse", cleanup)
    # The command boundary itself must not raise for any of these outcomes.
    assert worker.handle_payload(_stop_payload("boundary-job")) is None
    assert kills == [(HORSE_PID, int(signal.SIGKILL))]
    report = getattr(worker, "_nested_cleanup_report", None)
    if cleanup_fails:
        assert report is not None
        assert report.confirmed is False
        assert report.unconfirmed_reason == "OWNED_TREE_CLEANUP_UNCONFIRMED"
        assert report.unconfirmed_detail.startswith(
            f"{ProcessRunnerCleanupError.__name__}: "
        )
        assert worker._nested_cleanup_failed_pid == HORSE_PID
        # Fail closed: the run is reported failed, never acknowledged as stopped.
        assert (
            acknowledgement_decision(worker, "boundary-job", HORSE_PID, ROOT_STARTTIME)
            is None
        )
    else:
        # Positive control: the same boundary acknowledges a confirmed cleanup, so
        # the refusal above cannot be a blanket "never acknowledge".
        assert report is None
        assert (
            acknowledgement_decision(worker, "boundary-job", HORSE_PID, ROOT_STARTTIME)
            == "boundary-job"
        )


def test_redis_transport_failures_still_reach_the_pubsub_retry_path(monkeypatch):
    """The guard must not swallow RQ's connection-retry contract.

    ``redis.exceptions`` errors are not ``OSError`` subclasses, so a transport
    failure raised while handling a command keeps propagating to RQ's pubsub
    exception handler, which is what reconnects the worker.
    """
    worker = _worker()

    def unavailable(self, sig=signal.SIGKILL):
        raise RedisConnectionError("controlled transport failure")

    monkeypatch.setattr(timeouts.DurableWorker, "kill_horse", unavailable)
    with pytest.raises(RedisConnectionError):
        worker.handle_payload(_stop_payload("boundary-job"))


def test_a_refused_stop_does_not_poison_the_next_job(monkeypatch):
    """Cleanup markers are keyed to one job and released when its monitor ends.

    A pid is reusable, so a failure marker left behind could refuse the
    acknowledgement of a later, unrelated job.
    """
    worker = _worker()

    def failing(self, sig=signal.SIGKILL):
        self._nested_cleanup_failed_pid = self.horse_pid
        self._nested_cleanup_completed = None
        raise ProcessRunnerCleanupError(
            "Nested workflow cleanup could not be confirmed."
        )

    monkeypatch.setattr(timeouts.DurableWorker, "kill_horse", failing)
    assert worker.handle_payload(_stop_payload("first-job")) is None
    assert (
        acknowledgement_decision(worker, "first-job", HORSE_PID, ROOT_STARTTIME) is None
    )

    # What ``monitor_work_horse`` does in its ``finally`` once the job is over.
    worker._release_consumed_stop("first-job")
    assert worker.__dict__.get("_nested_stopped_job_id") is None
    assert not {
        "_nested_cleanup_failed_pid",
        "_nested_cleanup_completed",
        "_nested_cleanup_unconfirmed",
    } & set(worker.__dict__)

    # The next job on the same, reused pid reaches an acknowledgement.
    worker._nested_stopped_job_id = "second-job"
    worker._nested_cleanup_completed = ("second-job", HORSE_PID, ROOT_STARTTIME)
    worker._nested_cleanup_failed_pid = None
    assert (
        acknowledgement_decision(worker, "second-job", HORSE_PID, ROOT_STARTTIME)
        == "second-job"
    )


def test_the_proof_function_itself_still_refuses_incomplete_cleanup():
    """The conversion happens at the boundary, not inside the proof.

    ``_kill_owned_horse_tree`` keeps raising, so its caller is the single place
    that decides how a refusal is delivered. If the proof instead returned an
    unconfirmed report, its ``registered`` completeness assumption would be
    silently reused by callers as if the tree had been proven.
    """
    kills = []
    with pytest.raises(ProcessRunnerCleanupError) as error:
        timeouts._kill_owned_horse_tree(
            1, ROOT_STARTTIME, lambda: kills.append(ROOT_STARTTIME)
        )
    assert str(error.value) == "Nested workflow root was lost."
    assert kills == []
