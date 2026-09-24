"""Focused private-diagnostics contracts; no Docker, Redis or real worker.

These tests prove the two capture kinds exist, stay owner-only, stay outside the
uploadable evidence tree, and never replace the execution result they explain.
"""

from __future__ import annotations

import io
import json
import os
from pathlib import Path
import stat
from types import SimpleNamespace

import pytest

from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services import runs as runs_module

from . import platform_harness
from .cancellation_diagnostics import (
    WorkerStreamCapture,
    write_cancellation_snapshot,
)
from .test_runner_composition import _private_config_harness

PRIVATE_DIAGNOSTICS_ENV = "HELIXWEAVE_BULK_RNASEQ_PRIVATE_DIAGNOSTICS_DIR"


def _private_root(harness) -> Path:
    root = harness._resolve_private_diagnostics_root()
    assert root is not None
    return root


def _captures(harness, prefix: str) -> list[Path]:
    return sorted(
        path
        for path in _private_root(harness).iterdir()
        if path.name.startswith(f"{prefix}-")
    )


class _FakeWorkerProcess:
    returncode = None

    def __init__(self, pid: int = 4242) -> None:
        self.pid = pid

    def poll(self):
        return None


def _stub_worker_spawn(
    monkeypatch: pytest.MonkeyPatch, harness, *, error: bool = False
) -> dict[str, object]:
    monkeypatch.setattr(harness, "_worker_environment", lambda: {})
    recorded: dict[str, object] = {}

    def fake_popen(_argv, **kwargs):
        if error:
            raise OSError("worker could not be spawned")
        recorded.update(kwargs)
        return _FakeWorkerProcess()

    monkeypatch.setattr(platform_harness.subprocess, "Popen", fake_popen)
    return recorded


def test_worker_stream_capture_is_owner_only_and_sealed(tmp_path: Path) -> None:
    root = tmp_path / "private-diagnostics"
    root.mkdir(mode=0o700)
    capture = WorkerStreamCapture(root=root, ordinal=0)
    capture.stdout_handle.write("worker out\n")
    capture.stderr_handle.write("worker err\n")
    capture.close(returncode=7)
    capture.close(returncode=7)

    assert capture.stdout_path.read_text(encoding="utf-8") == "worker out\n"
    assert capture.stderr_path.read_text(encoding="utf-8") == "worker err\n"
    record = json.loads(
        (capture.directory / "worker-exit.json").read_text(encoding="utf-8")
    )
    assert record["returncode"] == 7
    assert stat.S_IMODE(capture.stdout_path.stat().st_mode) == 0o600
    assert stat.S_IMODE(capture.stderr_path.stat().st_mode) == 0o600
    assert stat.S_IMODE(capture.directory.stat().st_mode) == 0o700


def test_snapshot_is_timestamped_sealed_and_rejects_unsafe_labels(
    tmp_path: Path,
) -> None:
    root = tmp_path / "private-diagnostics"
    root.mkdir(mode=0o700)
    name = write_cancellation_snapshot(root=root, label="pre-close", document={"a": 1})
    assert name is not None
    assert name.startswith("pre-close-")
    destination = root / name
    payload = json.loads((destination / "snapshot.json").read_text(encoding="utf-8"))
    assert payload["label"] == "pre-close"
    assert payload["schema_version"] == "1.0.0"
    assert payload["captured_at_wall"].endswith("Z")
    assert isinstance(payload["captured_at_monotonic"], float)
    assert stat.S_IMODE(destination.stat().st_mode) == 0o700
    assert stat.S_IMODE((destination / "snapshot.json").stat().st_mode) == 0o600

    assert (
        write_cancellation_snapshot(root=root, label="../escape", document={}) is None
    )
    assert (
        write_cancellation_snapshot(
            root=tmp_path / "absent-root", label="pre-close", document={}
        )
        is None
    )
    assert [path.name for path in tmp_path.iterdir()] == ["private-diagnostics"]


def test_default_diagnostics_root_stays_outside_the_evidence_tree(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.delenv(PRIVATE_DIAGNOSTICS_ENV, raising=False)
    monkeypatch.delenv("RUNNER_TEMP", raising=False)
    harness = _private_config_harness(tmp_path, monkeypatch)

    root = _private_root(harness)
    evidence = harness.temporary_root / "evidence"
    assert root == harness.temporary_root / "private-diagnostics"
    assert "evidence" not in root.parts
    assert not root.is_relative_to(evidence)
    assert stat.S_IMODE(root.stat().st_mode) == 0o700


def test_start_worker_streams_output_into_the_private_root(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    harness = _private_config_harness(tmp_path, monkeypatch)
    recorded = _stub_worker_spawn(monkeypatch, harness)

    process = harness.start_worker()
    root = _private_root(harness)
    assert len(harness._worker_streams) == 1
    owner, capture = harness._worker_streams[0]
    assert owner is process
    assert recorded["stdout"] is capture.stdout_handle
    assert recorded["stderr"] is capture.stderr_handle
    assert recorded["stdout"] != platform_harness.subprocess.DEVNULL
    assert recorded["stderr"] != platform_harness.subprocess.DEVNULL
    assert recorded["stdin"] == platform_harness.subprocess.DEVNULL
    assert recorded["text"] is True
    for handle, path in (
        (capture.stdout_handle, capture.stdout_path),
        (capture.stderr_handle, capture.stderr_path),
    ):
        assert isinstance(handle, io.TextIOBase)
        assert path.parent.parent == root
        assert path.name in {"worker.stdout", "worker.stderr"}
        assert stat.S_IMODE(path.stat().st_mode) == 0o600

    harness._close_worker_streams(process, returncode=0)
    assert harness._worker_streams == []
    assert list(root.glob("worker-00-*/worker-exit.json"))


def test_start_worker_closes_streams_when_the_child_cannot_start(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    harness = _private_config_harness(tmp_path, monkeypatch)
    _stub_worker_spawn(monkeypatch, harness, error=True)

    with pytest.raises(OSError):
        harness.start_worker()
    root = _private_root(harness)
    assert harness._worker_processes == []
    assert harness._worker_streams == []
    assert list(root.glob("worker-00-*/worker-exit.json"))


def test_unusable_private_root_fails_closed_without_masking(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    harness = _private_config_harness(tmp_path, monkeypatch)
    monkeypatch.setenv(
        PRIVATE_DIAGNOSTICS_ENV,
        str(harness.temporary_root / "evidence" / "raw-diagnostics"),
    )
    assert harness._resolve_private_diagnostics_root() is None
    assert harness._private_diagnostics_error == "ValueError"

    recorded = _stub_worker_spawn(monkeypatch, harness)
    process = harness.start_worker()
    assert process.pid == 4242
    assert recorded["stdout"] == platform_harness.subprocess.DEVNULL
    assert harness._worker_streams == []
    assert harness._capture_diagnostics_snapshot("pre-close") is None


def test_rq_not_terminal_capture_records_durable_state(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    harness = _private_config_harness(tmp_path, monkeypatch)
    submitted = SimpleNamespace(run_id="run-capture", job_id="job-capture")
    harness._submitted.append(submitted)
    record = SimpleNamespace(
        run_id="run-capture",
        workflow_id="bulk-rnaseq",
        status=RunStatus.RUNNING,
        created_at=object(),
        updated_at=object(),
        started_at=object(),
        ended_at=None,
        current_stage="execution",
        cancellation_reason=None,
        error=None,
    )
    assignment = SimpleNamespace(
        run_id="run-capture",
        job_id="job-capture",
        backend="rq",
        queue_name="capture-queue",
        created_at=object(),
        managed_container_scope=None,
        dispatched_at=object(),
        claimed_at=object(),
        cancellation_requested_at=object(),
        cancellation_reason="User requested cancellation.",
        cancellation_acknowledged_at=None,
        requeue_requested_at=None,
        requeue_confirmed_at=None,
    )
    event = SimpleNamespace(
        sequence=1,
        event_type="cancellation_requested",
        timestamp=object(),
        status=None,
        stage=None,
        message="User requested cancellation.",
        context={},
        issue=None,
    )
    result_state = SimpleNamespace(
        artifact_revision=0,
        artifact_attempt_id=None,
        artifact_attempt_status=None,
        qc_revision=0,
        qc_attempt_id=None,
        qc_attempt_status=None,
    )

    class FakeRunService:
        def __init__(self, *_args, **_kwargs) -> None:
            return None

        @staticmethod
        def get_run(run_id):
            assert run_id == submitted.run_id
            return record

        @staticmethod
        def get_execution_assignment(run_id):
            assert run_id == submitted.run_id
            return assignment

        @staticmethod
        def list_events(run_id, *, limit):
            assert run_id == submitted.run_id
            assert limit == 1000
            return [event]

        @staticmethod
        def get_result_state(run_id):
            return result_state

        @staticmethod
        def list_artifacts(_run_id):
            return ()

        @staticmethod
        def list_qc_metrics(_run_id):
            return ()

    class FakeJob:
        @staticmethod
        def get_status(*, refresh: bool):
            assert refresh is True
            return platform_harness.JobStatus.STARTED

        @staticmethod
        def refresh() -> None:
            return None

    monkeypatch.setattr(runs_module, "RunService", FakeRunService)
    monkeypatch.setattr(
        platform_harness,
        "open_existing_run_persistence",
        lambda _url: SimpleNamespace(repository=object(), close=lambda: None),
    )
    monkeypatch.setattr(
        platform_harness,
        "_wait_for_rq_terminal_status",
        lambda _job: platform_harness.JobStatus.STARTED,
    )
    harness._run_queue = SimpleNamespace(
        _queue=SimpleNamespace(fetch_job=lambda _job_id: FakeJob())
    )

    with pytest.raises(
        AssertionError, match="accepted RQ job lacks a non-success terminal state"
    ):
        harness.collect_terminal(submitted)

    captures = _captures(harness, "rq-not-terminal")
    assert len(captures) == 1
    payload = json.loads((captures[0] / "snapshot.json").read_text(encoding="utf-8"))
    assert payload["label"] == "rq-not-terminal"
    assert payload["note"] == "rq_status=started"
    assert payload["job_ids"] == ["job-capture"]
    sections = payload["sections"]
    assert "section_error" not in sections["sqlite"]
    sqlite = sections["sqlite"]["run-capture"]
    assert sqlite["run"]["status"] == "running"
    assert sqlite["run"]["ended_at"] is None
    assert sqlite["run"]["current_stage"] == "execution"
    assert sqlite["assignment"]["cancellation_requested_at"] is not None
    assert sqlite["assignment"]["cancellation_acknowledged_at"] is None
    assert [entry["event_type"] for entry in sqlite["events"]] == [
        "cancellation_requested"
    ]
    assert sqlite["artifact_count"] == 0
    assert sections["rq_job"]["job-capture"]["status_refreshed"] == "started"
    assert sections["redis"] == {"unavailable": "harness is not open"}
    assert "workers" in sections["processes"]

    evidence = harness.temporary_root / "evidence"
    assert not _private_root(harness).is_relative_to(evidence)
    leaked = (
        []
        if not evidence.exists()
        else [path for path in evidence.rglob("*") if path.is_file()]
    )
    assert leaked == []


def test_close_captures_pre_close_state_before_cleanup(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    harness = _private_config_harness(tmp_path, monkeypatch)
    submitted = SimpleNamespace(run_id="run-close", job_id="job-close")
    harness._submitted.append(submitted)
    destructive: list[str] = []

    class FakeJob:
        @staticmethod
        def delete() -> None:
            destructive.append("job")

    harness._run_queue = SimpleNamespace(
        _queue=SimpleNamespace(
            fetch_job=lambda _job_id: FakeJob(),
            delete=lambda: destructive.append("queue"),
        )
    )
    harness._connection = SimpleNamespace(
        close=lambda: destructive.append("connection")
    )
    monkeypatch.setattr(
        platform_harness,
        "ManagedContainerCleaner",
        lambda **_kwargs: SimpleNamespace(
            cleanup=lambda _scope: SimpleNamespace(is_failure=False)
        ),
    )
    monkeypatch.setattr(
        platform_harness,
        "open_existing_run_persistence",
        lambda _url: SimpleNamespace(repository=object(), close=lambda: None),
    )

    class FakeRunService:
        def __init__(self, *_args, **_kwargs) -> None:
            return None

        @staticmethod
        def get_run(_run_id):
            return SimpleNamespace(
                run_id="run-close",
                status=RunStatus.RUNNING,
                ended_at=None,
                cancellation_reason=None,
                error=None,
            )

        @staticmethod
        def get_execution_assignment(_run_id):
            return None

        @staticmethod
        def list_events(_run_id, *, limit):
            assert limit == 1000
            return ()

        @staticmethod
        def get_result_state(_run_id):
            return None

        @staticmethod
        def list_artifacts(_run_id):
            return ()

        @staticmethod
        def list_qc_metrics(_run_id):
            return ()

    monkeypatch.setattr(runs_module, "RunService", FakeRunService)
    order: list[str] = []
    real_capture = harness._capture_diagnostics_snapshot

    def spy(label, note=None):
        order.append(f"capture:{label}")
        return real_capture(label, note=note)

    monkeypatch.setattr(harness, "_capture_diagnostics_snapshot", spy)

    # No reference profile config exists, so cleanup confirms; the capture must
    # still precede every destructive call.
    harness.close()

    assert order == ["capture:pre-close"]
    assert destructive == ["job", "queue", "connection"]
    captures = _captures(harness, "pre-close")
    assert len(captures) == 1
    payload = json.loads((captures[0] / "snapshot.json").read_text(encoding="utf-8"))
    assert payload["run_ids"] == ["run-close"]
    assert payload["worker_session_ids"] == []
    sections = payload["sections"]
    sqlite = sections["sqlite"]["run-close"]
    assert sqlite["run"]["status"] == "running"
    assert sqlite["run"]["ended_at"] is None
    assert sqlite["assignment"] is None
    assert sqlite["events"] == []
    assert sqlite["result_state"] is None
    assert sqlite["artifact_count"] == 0
    assert sections["rq_job"]["job-close"]["present"] is True
    assert sections["processes"]["harness_pid"] == os.getpid()
    # A section that cannot be read is recorded as degraded, never raised.
    assert sections["containers"] == {"endpoint_error": "AttributeError"}
    assert os.getuid() == captures[0].stat().st_uid
