import os
import subprocess
import sys
import textwrap
from dataclasses import FrozenInstanceError
from pathlib import Path

import pytest


SRC_ROOT = Path(__file__).resolve().parents[2] / "src"


def test_run_status_values_and_terminal_detection():
    from encode_pipeline.platform.runs import RunStatus

    assert RunStatus.CREATED.value == "created"
    assert not RunStatus.CREATED.is_terminal
    assert RunStatus.SUCCEEDED.is_terminal
    assert RunStatus.FAILED.is_terminal
    assert RunStatus.CANCELLED.is_terminal


def test_can_transition_accepts_pr99_path():
    from encode_pipeline.platform.runs import RunStatus, can_transition

    assert can_transition(RunStatus.CREATED, RunStatus.VALIDATING)
    assert can_transition(RunStatus.VALIDATING, RunStatus.PLANNED)
    assert can_transition(RunStatus.VALIDATING, RunStatus.FAILED)
    assert can_transition(RunStatus.PLANNED, RunStatus.QUEUED)
    assert can_transition(RunStatus.QUEUED, RunStatus.RUNNING)
    assert can_transition(RunStatus.QUEUED, RunStatus.FAILED)
    assert can_transition(RunStatus.RUNNING, RunStatus.SUCCEEDED)
    assert can_transition(RunStatus.RUNNING, RunStatus.FAILED)
    assert can_transition(RunStatus.RUNNING, RunStatus.CANCELLED)


def test_can_transition_rejects_invalid_transitions():
    from encode_pipeline.platform.runs import RunStatus, can_transition

    assert not can_transition(RunStatus.CREATED, RunStatus.RUNNING)
    assert not can_transition(RunStatus.PLANNED, RunStatus.FAILED)
    assert not can_transition(RunStatus.RUNNING, RunStatus.PLANNED)
    assert not can_transition(RunStatus.SUCCEEDED, RunStatus.FAILED)
    assert not can_transition(RunStatus.FAILED, RunStatus.CANCELLED)


def test_all_active_statuses_can_transition_to_cancelled():
    from encode_pipeline.platform.runs import RunStatus, can_transition

    active = {
        RunStatus.CREATED,
        RunStatus.VALIDATING,
        RunStatus.PLANNED,
        RunStatus.QUEUED,
        RunStatus.RUNNING,
    }
    for status in active:
        assert can_transition(status, RunStatus.CANCELLED)


def test_all_terminal_statuses_are_absorbing():
    from encode_pipeline.platform.runs import RunStatus, can_transition

    terminal = {
        RunStatus.SUCCEEDED,
        RunStatus.FAILED,
        RunStatus.CANCELLED,
    }
    for status in terminal:
        assert not can_transition(status, RunStatus.CREATED)
        assert not can_transition(status, RunStatus.RUNNING)
        assert not can_transition(status, RunStatus.CANCELLED)


def test_require_transition_raises_on_invalid():
    from encode_pipeline.platform.runs import RunStatus, require_transition

    with pytest.raises(ValueError, match="Illegal transition"):
        require_transition(RunStatus.SUCCEEDED, RunStatus.FAILED)


def test_run_record_stores_status_and_error():
    from datetime import datetime, timezone
    from encode_pipeline.platform.runs import RunRecord, RunStatus
    from encode_pipeline.platform.results import Issue

    now = datetime.now(timezone.utc)
    error = Issue(code="RUN_FAILED", message="Run failed")
    record = RunRecord(
        run_id="run-1",
        workflow_id="wf-1",
        inputs={"config": {"samples": "samples.tsv"}},
        status=RunStatus.FAILED,
        created_at=now,
        updated_at=now,
        started_at=None,
        ended_at=now,
        current_stage=None,
        cancellation_reason=None,
        error=error,
        tags={"env": "test"},
    )

    assert record.run_id == "run-1"
    assert record.status == RunStatus.FAILED
    assert record.error == error
    assert record.ended_at == now


def test_run_record_is_frozen():
    from datetime import datetime, timezone
    from encode_pipeline.platform.runs import RunRecord, RunStatus

    now = datetime.now(timezone.utc)
    record = RunRecord(
        run_id="run-1",
        workflow_id="wf-1",
        inputs={},
        status=RunStatus.CREATED,
        created_at=now,
        updated_at=now,
        started_at=None,
        ended_at=None,
        current_stage=None,
        cancellation_reason=None,
        error=None,
        tags={},
    )

    with pytest.raises(FrozenInstanceError):
        record.status = RunStatus.RUNNING


def test_run_record_to_dict_is_json_ready():
    from datetime import datetime, timezone
    from encode_pipeline.platform.runs import RunRecord, RunStatus

    now = datetime.now(timezone.utc)
    record = RunRecord(
        run_id="run-1",
        workflow_id="wf-1",
        inputs={"config": {}},
        status=RunStatus.CREATED,
        created_at=now,
        updated_at=now,
        started_at=None,
        ended_at=None,
        current_stage=None,
        cancellation_reason=None,
        error=None,
        tags={"env": "test"},
    )

    data = record.to_dict()
    assert data["run_id"] == "run-1"
    assert data["status"] == "created"
    assert data["created_at"] is now
    assert data["inputs"] == {"config": {}}
    assert data["tags"] == {"env": "test"}


def test_run_record_normalizes_string_status():
    from datetime import datetime, timezone
    from encode_pipeline.platform.runs import RunRecord, RunStatus

    now = datetime.now(timezone.utc)
    record = RunRecord(
        run_id="run-1",
        workflow_id="wf-1",
        inputs={},
        status="running",
        created_at=now,
        updated_at=now,
        started_at=None,
        ended_at=None,
        current_stage=None,
        cancellation_reason=None,
        error=None,
        tags={},
    )

    assert record.status is RunStatus.RUNNING
    assert record.to_dict()["status"] == "running"


def test_run_record_rejects_invalid_status():
    from datetime import datetime, timezone
    from encode_pipeline.platform.runs import RunRecord

    now = datetime.now(timezone.utc)
    with pytest.raises(ValueError, match="Invalid run status"):
        RunRecord(
            run_id="run-1",
            workflow_id="wf-1",
            inputs={},
            status="not-a-status",
            created_at=now,
            updated_at=now,
            started_at=None,
            ended_at=None,
            current_stage=None,
            cancellation_reason=None,
            error=None,
            tags={},
        )


def test_run_record_rejects_non_issue_error():
    from datetime import datetime, timezone
    from encode_pipeline.platform.runs import RunRecord

    now = datetime.now(timezone.utc)
    with pytest.raises(ValueError, match="RunRecord error must be an Issue or None"):
        RunRecord(
            run_id="run-1",
            workflow_id="wf-1",
            inputs={},
            status="created",
            created_at=now,
            updated_at=now,
            started_at=None,
            ended_at=None,
            current_stage=None,
            cancellation_reason=None,
            error={"code": "BAD"},
            tags={},
        )


def test_run_event_stores_sequence_and_issue():
    from datetime import datetime, timezone
    from encode_pipeline.platform.runs import RunEvent, RunStatus
    from encode_pipeline.platform.results import Issue

    now = datetime.now(timezone.utc)
    issue = Issue(code="STAGE_WARN", message="Stage warning", severity="warning")
    event = RunEvent(
        event_id="evt-1",
        run_id="run-1",
        sequence=3,
        event_type="issue_added",
        timestamp=now,
        status=RunStatus.RUNNING,
        stage="alignment",
        message="Warning during alignment.",
        context={"sample": "S1"},
        issue=issue,
    )

    assert event.sequence == 3
    assert event.event_type == "issue_added"
    assert event.issue == issue
    assert event.to_dict()["sequence"] == 3
    assert event.to_dict()["issue"]["code"] == "STAGE_WARN"


def test_run_event_context_is_defensively_copied():
    from datetime import datetime, timezone
    from encode_pipeline.platform.runs import RunEvent, RunStatus

    context = {"sample": "S1"}
    event = RunEvent(
        event_id="evt-1",
        run_id="run-1",
        sequence=1,
        event_type="status_changed",
        timestamp=datetime.now(timezone.utc),
        status=RunStatus.CREATED,
        stage=None,
        message="Run created.",
        context=context,
        issue=None,
    )

    context["sample"] = "S2"
    assert event.context == {"sample": "S1"}


def test_run_event_rejects_non_issue_issue():
    from datetime import datetime, timezone
    from encode_pipeline.platform.runs import RunEvent, RunStatus

    with pytest.raises(ValueError, match="RunEvent issue must be an Issue or None"):
        RunEvent(
            event_id="evt-1",
            run_id="run-1",
            sequence=1,
            event_type="status_changed",
            timestamp=datetime.now(timezone.utc),
            status=RunStatus.CREATED,
            stage=None,
            message="Run created.",
            context={},
            issue={"code": "BAD"},
        )


def test_run_log_chunk_stores_sequence_and_lines():
    from datetime import datetime, timezone
    from encode_pipeline.platform.runs import RunLogChunk

    now = datetime.now(timezone.utc)
    chunk = RunLogChunk(
        chunk_id="log-1",
        run_id="run-1",
        stream_name="stdout",
        sequence=7,
        timestamp=now,
        lines=["line 1", "line 2"],
    )

    assert chunk.stream_name == "stdout"
    assert chunk.sequence == 7
    assert chunk.lines == ("line 1", "line 2")
    assert chunk.to_dict()["lines"] == ["line 1", "line 2"]


def test_run_log_chunk_lines_are_defensively_copied():
    from datetime import datetime, timezone
    from encode_pipeline.platform.runs import RunLogChunk

    lines = ["line 1"]
    chunk = RunLogChunk(
        chunk_id="log-1",
        run_id="run-1",
        stream_name="stdout",
        sequence=1,
        timestamp=datetime.now(timezone.utc),
        lines=lines,
    )

    lines.append("line 2")
    assert list(chunk.lines) == ["line 1"]


def test_run_artifact_ref_stores_uri_and_metadata():
    from datetime import datetime, timezone
    from encode_pipeline.platform.runs import RunArtifactRef

    now = datetime.now(timezone.utc)
    ref = RunArtifactRef(
        artifact_id="art-1",
        run_id="run-1",
        artifact_type="file",
        name="peaks.narrowPeak",
        uri="run://runs/run-1/artifacts/peaks.narrowPeak",
        mime_type="text/plain",
        produced_at=now,
        metadata={"sample": "S1"},
    )

    assert ref.artifact_type == "file"
    assert ref.uri == "run://runs/run-1/artifacts/peaks.narrowPeak"
    assert ref.to_dict()["metadata"] == {"sample": "S1"}


def test_run_artifact_ref_metadata_is_defensively_copied():
    from datetime import datetime, timezone
    from encode_pipeline.platform.runs import RunArtifactRef

    metadata = {"sample": "S1"}
    ref = RunArtifactRef(
        artifact_id="art-1",
        run_id="run-1",
        artifact_type="file",
        name="peaks.narrowPeak",
        uri="run://runs/run-1/artifacts/peaks.narrowPeak",
        mime_type=None,
        produced_at=datetime.now(timezone.utc),
        metadata=metadata,
    )

    metadata["sample"] = "S2"
    assert ref.metadata == {"sample": "S1"}


def test_platform_package_exports_run_primitives():
    from encode_pipeline.platform import (
        RunArtifactRef,
        RunEvent,
        RunLogChunk,
        RunRecord,
        RunStatus,
        can_transition,
        require_transition,
    )

    assert [
        RunArtifactRef.__name__,
        RunEvent.__name__,
        RunLogChunk.__name__,
        RunRecord.__name__,
        RunStatus.CREATED.value,
        can_transition.__name__,
        require_transition.__name__,
    ] == [
        "RunArtifactRef",
        "RunEvent",
        "RunLogChunk",
        "RunRecord",
        "created",
        "can_transition",
        "require_transition",
    ]
    assert can_transition(RunStatus.CREATED, RunStatus.VALIDATING)
    assert require_transition(RunStatus.CREATED, RunStatus.VALIDATING) is None


def test_importing_platform_runs_does_not_import_services_or_adapters():
    code = """
        import sys
        import encode_pipeline.platform.runs
        forbidden = {
            "encode_pipeline.services",
            "encode_pipeline.adapters",
            "encode_pipeline.api",
            "encode_pipeline.frontend",
            "fastapi",
            "pydantic",
            "snakemake",
            "subprocess",
            "openai",
        }
        found = [name for name in sys.modules if any(name.startswith(p) for p in forbidden)]
        print(found)
    """
    env = dict(os.environ)
    env["PYTHONPATH"] = str(SRC_ROOT)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    proc = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(code)],
        capture_output=True,
        text=True,
        env=env,
    )

    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "[]"
