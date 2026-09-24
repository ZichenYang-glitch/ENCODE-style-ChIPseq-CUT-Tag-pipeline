"""Real local POSIX children for the private Hi-TrAC qualification lifecycle."""

import ctypes
import json
import os
from pathlib import Path
import signal
import sys
import time

import pytest

from encode_pipeline.adapters.hitrac_preprocess import calls, process


def subreaper(value=None):
    libc = ctypes.CDLL(None, use_errno=True)
    if value is not None:
        assert libc.prctl(36, value, 0, 0, 0) == 0
    current = ctypes.c_int()
    assert libc.prctl(37, ctypes.byref(current), 0, 0, 0) == 0
    return current.value


def identity(pid):
    try:
        stat = Path(f"/proc/{pid}/stat").read_text()
    except FileNotFoundError:
        return None
    fields = stat.rsplit(")", 1)[1].split()
    return {"pid": pid, "state": fields[0], "starttime": int(fields[19])}


MARK_PROCESS = """
import json, os, signal, time
from pathlib import Path

def mark(pid, role):
    fields = Path('/proc/%s/stat' % pid).read_text().rsplit(')', 1)[1].split()
    entry = {'pid': pid, 'starttime': int(fields[19]), 'role': role}
    fd = os.open(os.environ['HITRAC_TEST_PROCESSES'], os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
    try:
        os.write(fd, (json.dumps(entry) + '\\n').encode())
    finally:
        os.close(fd)
"""


def records(path):
    return (
        [json.loads(line) for line in path.read_text().splitlines()]
        if path.exists()
        else []
    )


def assert_gone(entries):
    for entry in entries:
        current = identity(entry["pid"])
        # PID reuse is allowed; a same-identity zombie is expressly not success.
        assert current is None or current["starttime"] != entry["starttime"], current


@pytest.fixture
def local_processes(tmp_path):
    """Emergency cleanup is limited to this test's recorded PID/start identities."""
    baseline = subreaper()
    marker = tmp_path / "processes.jsonl"
    private = tmp_path / "private"
    private.mkdir(mode=0o700)
    env = {
        "PATH": "/usr/bin:/bin",
        "LC_ALL": "C.UTF-8",
        "HITRAC_TEST_PROCESSES": str(marker),
        "PYTHONDONTWRITEBYTECODE": "1",
    }
    try:
        yield marker, private, env, baseline
    finally:
        for entry in records(marker):
            current = identity(entry["pid"])
            if current is not None and current["starttime"] == entry["starttime"]:
                try:
                    os.kill(entry["pid"], signal.SIGKILL)
                except ProcessLookupError:
                    pass
        deadline = time.monotonic() + 2
        remaining = records(marker)
        while remaining and time.monotonic() < deadline:
            waiting = []
            for entry in remaining:
                try:
                    os.waitpid(entry["pid"], os.WNOHANG)
                except ChildProcessError:
                    pass
                current = identity(entry["pid"])
                if current is not None and current["starttime"] == entry["starttime"]:
                    waiting.append(entry)
            remaining = waiting
            if remaining:
                time.sleep(0.01)
        subreaper(baseline)
        assert_gone(records(marker))


def launch(tmp_path, local_processes, body, *, timeout=3, cancelled=lambda: False):
    marker, private, env, _ = local_processes
    source = tmp_path / "child.py"
    source.write_text(MARK_PROCESS + "\nmark(os.getpid(), 'parent')\n" + body)
    execution = process.execute(
        [sys.executable, "-I", "-B", str(source)],
        tmp_path,
        env,
        private,
        timeout,
        cancelled,
    )
    return execution, records(marker)


@pytest.mark.parametrize("previous", [0, 1])
def test_normal_child_preserves_private_streams_and_subreaper(
    tmp_path, local_processes, previous
):
    subreaper(previous)
    execution, children = launch(
        tmp_path,
        local_processes,
        "import sys\nprint('private-scientific-output')\nprint('private-warning', file=sys.stderr)\n",
    )
    assert execution.returncode == 0
    assert execution.reason is None
    assert subreaper() == previous
    assert (
        local_processes[1] / "upstream.stdout"
    ).read_text() == "private-scientific-output\n"
    assert (local_processes[1] / "upstream.stderr").read_text() == "private-warning\n"
    assert_gone(children)


def test_nonzero_and_signal_are_not_rewritten_as_success(tmp_path, local_processes):
    execution, children = launch(
        tmp_path, local_processes, "os.kill(os.getpid(), signal.SIGTERM)\n"
    )
    assert execution.returncode == -signal.SIGTERM
    assert execution.reason is None
    assert subreaper() == local_processes[3]
    assert_gone(children)


def test_child_exit_status_is_retained(tmp_path, local_processes):
    execution, children = launch(tmp_path, local_processes, "raise SystemExit(73)\n")
    assert execution.returncode == 73
    assert execution.reason is None
    assert subreaper() == local_processes[3]
    assert_gone(children)


def test_successful_parent_and_grandchild_are_not_reported_as_leftovers(
    tmp_path, local_processes
):
    execution, children = launch(
        tmp_path,
        local_processes,
        "child = os.fork()\n"
        "if child == 0:\n"
        "    mark(os.getpid(), 'grandchild')\n"
        "    os._exit(0)\n"
        "os.waitpid(child, 0)\n",
    )
    assert execution.returncode == 0
    assert execution.reason is None
    assert {entry["role"] for entry in children} == {"parent", "grandchild"}
    assert subreaper() == local_processes[3]
    assert_gone(children)


TREE = """
child = os.fork()
if child == 0:
    signal.signal(signal.SIGTERM, signal.SIG_IGN)
    mark(os.getpid(), 'grandchild')
    while True:
        time.sleep(0.02)
"""


@pytest.mark.parametrize("reason", ["cancelled", "timed_out"])
def test_cancel_timeout_reap_ignoring_grandchild(tmp_path, local_processes, reason):
    marker = local_processes[0]
    execution, children = launch(
        tmp_path,
        local_processes,
        TREE + "while True:\n    time.sleep(0.02)\n",
        timeout=0.5 if reason == "timed_out" else 3,
        cancelled=(lambda: len(records(marker)) == 2)
        if reason == "cancelled"
        else lambda: False,
    )
    assert execution.reason == reason
    assert execution.returncode != 0
    assert {entry["role"] for entry in children} == {"parent", "grandchild"}
    assert subreaper() == local_processes[3]
    assert_gone(children)


def test_parent_exit_with_live_grandchild_is_rejected(tmp_path, local_processes):
    execution, children = launch(tmp_path, local_processes, TREE + "os._exit(0)\n")
    assert execution.returncode == 0
    assert execution.reason == "child_survived_parent"
    assert {entry["role"] for entry in children} == {"parent", "grandchild"}
    assert subreaper() == local_processes[3]
    assert_gone(children)


def test_wrapper_death_keeps_unfinished_receipt_and_reaps_real_child(
    tmp_path, local_processes
):
    marker, private, env, baseline = local_processes
    tool = tmp_path / "fake-samtools"
    tool.write_text(
        f"#!{sys.executable}\n"
        + MARK_PROCESS
        + "\nsignal.signal(signal.SIGTERM, signal.SIG_IGN)\n"
        + "mark(os.getpid(), 'real-tool')\nmark(os.getppid(), 'wrapper')\n"
        + "os.kill(os.getppid(), signal.SIGKILL)\n"
        + "while True:\n    time.sleep(0.02)\n"
    )
    tool.chmod(0o700)
    plan = [
        {"tool": "samtools", "stage": "view", "sample": "s000001", "argv": ["view"]}
    ]
    shim, config = calls.prepare_shims(
        tmp_path, {"samtools": tool}, plan, "test-runtime"
    )
    execution, children = launch(
        tmp_path,
        local_processes,
        f"import subprocess\nsubprocess.call([{str(shim / 'samtools')!r}, 'view'])\n",
    )
    assert execution.returncode == 0
    assert execution.reason == "child_survived_parent"
    with pytest.raises(calls.CallFailure) as error:
        calls.verify_calls(config)
    assert (error.value.code, error.value.sample, error.value.stage) == (
        "call_unfinished",
        "s000001",
        "view",
    )
    receipts = list((private / "calls").iterdir())
    assert len(receipts) == 1 and receipts[0].name.endswith(".start.json")
    assert {entry["role"] for entry in children} == {"parent", "wrapper", "real-tool"}
    assert subreaper() == baseline
    assert_gone(children)


def test_sigkill_from_tool_remains_signal_termination_at_wrapper(
    tmp_path, local_processes
):
    _, private, _, baseline = local_processes
    tool = tmp_path / "killed-samtools"
    tool.write_text(
        f"#!{sys.executable}\n"
        + MARK_PROCESS
        + "\nmark(os.getpid(), 'real-tool')\n"
        + "mark(os.getppid(), 'wrapper')\n"
        + "os.kill(os.getpid(), signal.SIGKILL)\n"
    )
    tool.chmod(0o700)
    plan = [
        {"tool": "samtools", "stage": "view", "sample": "s000001", "argv": ["view"]}
    ]
    shim, config = calls.prepare_shims(
        tmp_path, {"samtools": tool}, plan, "test-runtime"
    )
    result_file = private / "wrapper-returncode.json"
    execution, children = launch(
        tmp_path,
        local_processes,
        "import subprocess\n"
        + f"code = subprocess.call([{str(shim / 'samtools')!r}, 'view'])\n"
        + f"Path({str(result_file)!r}).write_text(json.dumps(code))\n",
    )
    assert execution.returncode == 0 and execution.reason is None
    with pytest.raises(calls.CallFailure) as error:
        calls.verify_calls(config)
    assert (error.value.code, error.value.sample, error.value.stage) == (
        "child_failed",
        "s000001",
        "view",
    )
    assert json.loads(result_file.read_text()) == -signal.SIGKILL
    assert subreaper() == baseline
    assert_gone(children)


def test_callback_error_still_cleans_attempt_and_restores_subreaper(
    tmp_path, local_processes
):
    def broken_callback():
        if len(records(local_processes[0])) == 2:
            raise ValueError("private-controller-error")
        return False

    with pytest.raises(ValueError, match="private-controller-error"):
        launch(
            tmp_path,
            local_processes,
            TREE + "while True:\n    time.sleep(0.02)\n",
            cancelled=broken_callback,
        )
    assert subreaper() == local_processes[3]
    assert_gone(records(local_processes[0]))


def test_grandchild_reaping_does_not_steal_direct_child_status(
    tmp_path, local_processes, monkeypatch
):
    real_reap = process._reap_group
    readiness = tmp_path / "signal-ready"

    def schedule_exit_between_poll_and_reap(group):
        # Real parent exit, controlled scheduling at the race boundary. The
        # original group-reaper must not reap the direct Popen child first.
        deadline = time.monotonic() + 0.5
        while time.monotonic() < deadline:
            current = identity(group)
            if current is None or current["state"] == "Z":
                break
            time.sleep(0.001)
        real_reap(group)

    monkeypatch.setattr(process, "_reap_group", schedule_exit_between_poll_and_reap)
    execution, children = launch(
        tmp_path,
        local_processes,
        "def terminate(sig, frame):\n"
        "    time.sleep(0.02)\n"
        "    os._exit(73)\n"
        "signal.signal(signal.SIGTERM, terminate)\n"
        + f"Path({str(readiness)!r}).write_text('ready')\n"
        + "while True:\n    time.sleep(0.02)\n",
        cancelled=readiness.exists,
    )
    assert execution.reason == "cancelled"
    assert execution.returncode == 73
    assert subreaper() == local_processes[3]
    assert_gone(children)


def test_cleanup_error_still_restores_subreaper(tmp_path, local_processes, monkeypatch):
    # Inject only failed cleanup verification; the recorded real child has exited.
    subreaper(0)
    monkeypatch.setattr(process, "_exists", lambda group: True)
    monkeypatch.setattr(process.os, "killpg", lambda group, sig: None)
    with pytest.raises(RuntimeError, match="attempt_process_cleanup_incomplete"):
        launch(tmp_path, local_processes, "pass\n")
    assert subreaper() == 0
    assert_gone(records(local_processes[0]))
