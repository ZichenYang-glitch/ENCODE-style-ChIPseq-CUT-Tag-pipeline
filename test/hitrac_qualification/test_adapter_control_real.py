"""Real H3 plans under runner timeout and RQ horse stop, without a Redis service.

Only the process-control fault (SIGSTOP after a scientific start receipt) is
injected. The admitted command, qualifier, tracPre2 and scientific tools remain
unchanged. These tests cover worker termination, not queue/lifecycle acceptance.
"""

import json
import logging
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
import traceback

import pytest

from encode_pipeline.adapters.hitrac_preprocess.admission import sha256_file
from encode_pipeline.adapters.hitrac_preprocess import entrypoint
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.defaults import create_default_process_runner
from encode_pipeline.workers.settings import WorkerSettings
from encode_pipeline.workers.timeouts import DurableWorker
import test_adapter_real as adapter_cases
from test_adapter_real import composed, record
from test_qualification_real import inject_tool, rows

pytestmark = pytest.mark.real_execution
coordinates = adapter_cases.coordinates


def _stat(pid):
    try:
        text = Path(f"/proc/{pid}/stat").read_text()
    except FileNotFoundError:
        return None
    fields = text[text.rfind(")") + 2 :].split()
    return {"pid": pid, "starttime": int(fields[19]), "group": int(fields[2])}


def _same(identity):
    current = _stat(identity["pid"])
    return current is not None and current["starttime"] == identity["starttime"]


def _tree(pid):
    found = {}
    pending = [pid]
    while pending:
        current = pending.pop()
        identity = _stat(current)
        if identity is None or current in found:
            continue
        found[current] = identity
        try:
            children = Path(f"/proc/{current}/task/{current}/children").read_text()
        except FileNotFoundError:
            continue
        pending.extend(map(int, children.split()))
    return found


def _scientific_start(workspace, horse, deadline):
    directory = workspace / "hitrac-attempt/private/calls"
    while time.monotonic() < deadline:
        for path in sorted(directory.glob("*.start.json")):
            try:
                receipt = json.loads(path.read_text())
            except json.JSONDecodeError:
                # Receipt creation is exclusive but its write may be in flight.
                continue
            if receipt["stage"] != "align" or not receipt["sample"]:
                continue
            wrapper = _stat(receipt["pid"])
            if wrapper is None:
                continue
            delegate_tree = _tree(receipt["pid"])
            if len(delegate_tree) < 2:
                # A receipt precedes subprocess.call; wait for the real delegate.
                continue
            scientific_group = wrapper["group"]
            assert scientific_group not in (horse, os.getpgrp())
            assert receipt["pid"] in _tree(horse)
            leader = _stat(scientific_group)
            assert leader is not None
            command = Path(f"/proc/{scientific_group}/cmdline").read_bytes()
            assert b"tracPre2.py" in command
            # Deterministic control fault: the real tool has actually started,
            # and cannot race to successful completion before the outer stop.
            os.killpg(scientific_group, signal.SIGSTOP)
            return {
                "receipt_path": str(path),
                "receipt": receipt,
                "delegate_tree_before_signal": list(delegate_tree.values()),
                "science_group": leader,
                "injected_signal": "SIGSTOP",
                "injected_monotonic": time.monotonic(),
                "observed_tree": list(_tree(horse).values()),
            }
        if _stat(horse) is None:
            raise AssertionError("horse exited before a scientific call started")
        time.sleep(0.005)
    raise AssertionError("scientific start receipt was not observed before deadline")


def _wait_horse(pid, deadline):
    while time.monotonic() < deadline:
        waited, status = os.waitpid(pid, os.WNOHANG)
        if waited:
            return os.waitstatus_to_exitcode(status)
        time.sleep(0.02)
    raise AssertionError("owned horse did not exit")


def _assert_stopped(identities):
    deadline = time.monotonic() + 3
    while time.monotonic() < deadline and any(map(_same, identities)):
        time.sleep(0.02)
    assert not [entry for entry in identities if _same(entry)]


@pytest.mark.parametrize("control", ["runner_timeout", "worker_kill_horse"])
def test_real_adapter_outer_control_cleans_nested_scientific_processes(
    coordinates, tmp_path, control
):
    adapter, runner, spec, workspace, inputs = composed(coordinates, tmp_path)
    if control == "runner_timeout":
        runner = create_default_process_runner(
            registry=WorkflowRegistry([adapter]),
            settings=WorkerSettings(
                database_url=f"sqlite:///{tmp_path / 'unused-control.db'}",
                redis_url="redis://127.0.0.1:1/0",
                queue_name="unused",
                workspace_root=tmp_path,
                job_timeout_seconds=30,
            ),
        )
    before = {
        path: sha256_file(Path(path))
        for row in inputs.samples
        for key, path in row.items()
        if key.startswith("fastq")
    }
    sentinel = subprocess.Popen(
        [sys.executable, "-I", "-B", "-c", "import time; time.sleep(120)"],
        cwd=tmp_path,
        start_new_session=True,
    )
    horse = os.fork()
    if horse == 0:
        try:
            os.setpgrp()
            result = runner.run(spec)
            record(tmp_path, result, spec)
        except BaseException:
            (tmp_path / "horse-error.txt").write_text(traceback.format_exc())
            os._exit(2)
        os._exit(0)

    captured = []
    reaped = False
    observation = {"control": control, "argv": spec.argv, "cwd": spec.cwd}
    try:
        deadline = time.monotonic() + 65
        observation.update(_scientific_start(workspace, horse, deadline))
        captured = observation["observed_tree"]
        assert len(captured) >= 4  # horse, qualifier, original script, shim
        assert os.getpgid(horse) == horse
        assert sentinel.poll() is None
        if control == "worker_kill_horse":
            # The real worker stop method needs only its owned horse identity;
            # no Redis connection or fabricated job success is involved.
            worker = object.__new__(DurableWorker)
            worker._horse_pid = horse
            worker.name = "hitrac-control-qualification"
            worker.log = logging.getLogger("test.hitrac-control")
            worker.kill_horse()
            observation["worker_stop_returned"] = True
        exit_code = _wait_horse(horse, deadline)
        reaped = True
        observation["horse_exit_code"] = exit_code
        assert exit_code == (0 if control == "runner_timeout" else -signal.SIGKILL)
        if control == "runner_timeout":
            result = json.loads((tmp_path / "runner-result.json").read_text())
            assert result["result"] is None
            assert [issue["code"] for issue in result["issues"]] == [
                "PROCESS_RUNNER_TIMEOUT"
            ]
        _assert_stopped(captured)
        observation["captured_processes_gone"] = True
        assert sentinel.poll() is None
        observation["unrelated_sentinel_alive"] = True
        attempt = workspace / "hitrac-attempt"
        assert not (attempt / "complete.json").exists()
        outcome = attempt / "private/outcome.json"
        if outcome.exists():
            observation["private_outcome"] = json.loads(outcome.read_text())
            assert observation["private_outcome"]["status"] != "complete"
        assert (attempt / "private/command.json").is_file()
        assert list((attempt / "private/calls").glob("*.start.json"))
        assert {path: sha256_file(Path(path)) for path in before} == before
        observation["original_fastq_sha256_unchanged"] = before
        assert adapter.execution_availability().execution == "not_configured"
    finally:
        (tmp_path / "control-observation.json").write_text(
            json.dumps(observation, indent=2) + "\n"
        )
        # Emergency cleanup is confined to identities created by this test.
        if not reaped:
            captured.extend(_tree(horse).values())
            for identity in reversed(captured):
                if _same(identity):
                    try:
                        os.kill(identity["pid"], signal.SIGKILL)
                    except ProcessLookupError:
                        pass
            try:
                os.waitpid(horse, 0)
            except ChildProcessError:
                pass
        sentinel.terminate()
        sentinel.wait(timeout=5)


def test_adapter_bound_entrypoint_rejects_original_exit73_fault(
    coordinates, tmp_path, monkeypatch, capsys
):
    """Real tools, narrow injected child fault, direct main (not subprocess CLI)."""
    _, _, spec, workspace, inputs = composed(coordinates, tmp_path)
    before = {
        path: sha256_file(Path(path))
        for row in inputs.samples
        for key, path in row.items()
        if key.startswith("fastq")
    }
    inject_tool(monkeypatch, coordinates, tmp_path, "view", partial=True)
    monkeypatch.chdir(workspace)
    monkeypatch.setattr(
        sys,
        "argv",
        ["hitrac-private-entry", *spec.argv[spec.argv.index("--request") :]],
    )
    handlers = {sig: signal.getsignal(sig) for sig in (signal.SIGTERM, signal.SIGINT)}
    try:
        code = entrypoint.main()
    finally:
        for sig, handler in handlers.items():
            signal.signal(sig, handler)
    captured = capsys.readouterr()
    assert code == 1
    response = json.loads(captured.out)
    assert response == {
        "status": "rejected",
        "reason_code": "child_failed",
        "sample": "s000001",
        "stage": "view",
    }
    attempt = workspace / "hitrac-attempt"
    assert not (attempt / "complete.json").exists()
    assert (attempt / "output/tracPre_summary.txt").is_file()
    assert len(rows(attempt / "output/s000001/s000001_all.bedpe.gz")) == 7
    assert len(rows(attempt / "output/s000001/s000001_unique.bedpe.gz")) == 4
    assert (tmp_path / "before-drop.bam").is_file()
    assert (
        json.loads((attempt / "private/execution.json").read_text())["returncode"] == 0
    )
    ends = [
        json.loads(path.read_text())
        for path in (attempt / "private/calls").glob("*.end.json")
    ]
    assert [item["returncode"] for item in ends].count(73) == 1
    assert {path: sha256_file(Path(path)) for path in before} == before
    (tmp_path / "entry-fault-observation.json").write_text(
        json.dumps(
            {
                "execution_boundary": "adapter plan and direct original entrypoint.main",
                "fault": "one admitted samtools delegate; real output then exit 73",
                "argv": spec.argv,
                "exit_code": code,
                "stdout": captured.out,
                "stderr": captured.err,
                "response": response,
                "original_fastq_sha256_unchanged": before,
            },
            indent=2,
        )
        + "\n"
    )
