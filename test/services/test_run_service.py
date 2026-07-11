"""Tests for the in-memory run lifecycle service."""

from __future__ import annotations

import os
import subprocess
import sys
import textwrap
from datetime import datetime, timezone
from pathlib import Path

import pytest

from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagPreview,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
    WorkspacePlan,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.runs import RunArtifactRef, RunRecord, RunStatus
from encode_pipeline.services.runs import RunService


class FakeAdapter:
    """Minimal workflow adapter for service-layer tests."""

    def __init__(
        self,
        workflow_id: str = "fake",
        supports: tuple[str, ...] = ("validation",),
    ) -> None:
        self.metadata = WorkflowMetadata(
            workflow_id=workflow_id,
            name="Fake Workflow",
            version="1.0.0",
        )
        self.capabilities = WorkflowCapabilities(supports=supports)

    def schema(self) -> WorkflowSchema:
        return WorkflowSchema(config_schema={"type": "object"})

    def validate(self, inputs: WorkflowInputs) -> Result[object]:
        return Result.success({"validated": True})

    def preview_dag(self, inputs: WorkflowInputs) -> Result[DagPreview]:
        return Result.success(DagPreview())

    def plan_workspace(
        self,
        inputs: WorkflowInputs,
        workspace: str | Path,
    ) -> Result[WorkspacePlan]:
        return Result.success(WorkspacePlan(directories=[str(workspace)]))

    def build_command(self, plan: WorkspacePlan) -> Result[CommandSpec]:
        return Result.success(CommandSpec(argv=["run-workflow"]))


def test_service_rejects_non_workflow_registry_registry():
    with pytest.raises(ValueError, match="WorkflowRegistry"):
        RunService(registry=object())


def test_create_run_stores_record_and_emits_created_event():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    inputs = WorkflowInputs(config={"samples": "samples.tsv"})

    record = service.create_run("fake", inputs, tags={"env": "test"})

    assert isinstance(record, RunRecord)
    assert record.run_id == "run-1"
    assert record.workflow_id == "fake"
    assert record.status is RunStatus.CREATED
    assert record.inputs == inputs.to_dict()
    assert record.tags == {"env": "test"}
    assert record.started_at is None
    assert record.ended_at is None
    assert record.current_stage is None
    assert record.cancellation_reason is None
    assert record.error is None

    events = service.list_events("run-1")
    assert len(events) == 1
    assert events[0].sequence == 1
    assert events[0].event_type == "status_changed"
    assert events[0].status is RunStatus.CREATED
    assert events[0].message == "Run created."


def test_create_run_uses_canonical_workflow_id_from_adapter_metadata():
    registry = WorkflowRegistry(adapters=[FakeAdapter(workflow_id="canonical-fake")])
    service = RunService(registry=registry, id_factory=lambda: "run-1")

    record = service.create_run("  canonical-fake  ", WorkflowInputs(config={}))

    assert record.workflow_id == "canonical-fake"


def test_create_run_rejects_duplicate_generated_run_id():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))

    with pytest.raises(ValueError, match="Duplicate run_id"):
        service.create_run("fake", WorkflowInputs(config={}))


def test_create_run_unknown_workflow_raises_key_error():
    service = RunService(registry=WorkflowRegistry(), id_factory=lambda: "run-1")

    with pytest.raises(KeyError):
        service.create_run("missing", WorkflowInputs(config={}))


def test_get_run_returns_record():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    created = service.create_run("fake", WorkflowInputs(config={}))

    record = service.get_run("run-1")

    assert record is created


def test_get_run_unknown_raises_key_error():
    service = RunService(registry=WorkflowRegistry(), id_factory=lambda: "run-1")

    with pytest.raises(KeyError):
        service.get_run("missing")


def test_list_runs_returns_creation_order():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry)
    run_a = service.create_run("fake", WorkflowInputs(config={}), tags={"name": "a"})
    run_b = service.create_run("fake", WorkflowInputs(config={}), tags={"name": "b"})

    runs = service.list_runs()

    assert runs == (run_a, run_b)


def test_list_runs_empty_when_no_runs():
    service = RunService(registry=WorkflowRegistry())

    assert service.list_runs() == ()


def test_transition_run_follows_pr99_graph():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))

    record = service.transition_run("run-1", RunStatus.VALIDATING)
    assert record.status is RunStatus.VALIDATING

    record = service.transition_run("run-1", RunStatus.PLANNED)
    assert record.status is RunStatus.PLANNED

    record = service.transition_run("run-1", RunStatus.QUEUED)
    assert record.status is RunStatus.QUEUED

    record = service.transition_run("run-1", RunStatus.RUNNING)
    assert record.status is RunStatus.RUNNING
    assert record.started_at is not None

    record = service.transition_run("run-1", RunStatus.SUCCEEDED)
    assert record.status is RunStatus.SUCCEEDED
    assert record.ended_at is not None


def test_transition_run_rejects_illegal_transition():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))

    with pytest.raises(ValueError, match="Illegal transition"):
        service.transition_run("run-1", RunStatus.RUNNING)


def test_transition_run_updates_current_stage():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))

    record = service.transition_run(
        "run-1",
        RunStatus.VALIDATING,
        stage="validation",
    )

    assert record.current_stage == "validation"


def test_transition_run_failed_sets_error_issue():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))
    service.transition_run("run-1", RunStatus.VALIDATING)
    issue = Issue(code="VALIDATION_FAILED", message="Validation failed.")

    record = service.transition_run(
        "run-1",
        RunStatus.FAILED,
        issue=issue,
    )

    assert record.status is RunStatus.FAILED
    assert record.error is issue
    assert record.ended_at is not None


def test_transition_run_emits_one_status_changed_event():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))

    service.transition_run("run-1", RunStatus.VALIDATING)

    events = service.list_events("run-1")
    assert len(events) == 2
    assert events[1].event_type == "status_changed"
    assert events[1].status is RunStatus.VALIDATING


def test_cancel_run_transitions_active_run():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))

    record = service.cancel_run("run-1", reason="User requested cancellation.")

    assert record.status is RunStatus.CANCELLED
    assert record.cancellation_reason == "User requested cancellation."
    assert record.ended_at is not None


def test_cancel_run_is_idempotent_for_terminal_run():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))
    service.transition_run("run-1", RunStatus.VALIDATING)
    service.transition_run("run-1", RunStatus.FAILED)
    before_events = service.list_events("run-1")

    record = service.cancel_run("run-1", reason="Should be ignored.")

    assert record.status is RunStatus.FAILED
    assert record.cancellation_reason is None
    assert service.list_events("run-1") == before_events


def test_cancel_run_preserves_started_at_and_sets_ended_at_for_running_run():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))
    service.transition_run("run-1", RunStatus.VALIDATING)
    service.transition_run("run-1", RunStatus.PLANNED)
    service.transition_run("run-1", RunStatus.QUEUED)
    running = service.transition_run("run-1", RunStatus.RUNNING)

    record = service.cancel_run("run-1", reason="User requested cancellation.")

    assert record.status is RunStatus.CANCELLED
    assert record.started_at is running.started_at
    assert record.ended_at is not None
    assert record.cancellation_reason == "User requested cancellation."


def test_add_event_records_event_without_mutating_status():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))

    event = service.add_event(
        "run-1",
        "stage_started",
        "Alignment started.",
        stage="alignment",
    )

    assert event.event_type == "stage_started"
    assert event.sequence == 2
    assert event.stage == "alignment"
    record = service.get_run("run-1")
    assert record.status is RunStatus.CREATED
    assert record.current_stage is None


def test_list_events_returns_ordered_events():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))
    service.add_event("run-1", "stage_started", "Stage 1")
    service.add_event("run-1", "stage_completed", "Stage 2")

    events = service.list_events("run-1")

    assert [event.sequence for event in events] == [1, 2, 3]


def test_list_events_with_after_and_limit():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))
    service.add_event("run-1", "stage_started", "Stage 1")
    service.add_event("run-1", "stage_completed", "Stage 2")
    service.add_event("run-1", "issue_added", "Stage 3")

    events = service.list_events("run-1", after="evt-1", limit=1)

    assert len(events) == 1
    assert events[0].sequence == 2


def test_list_events_unknown_cursor_raises_key_error():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))

    with pytest.raises(KeyError):
        service.list_events("run-1", after="evt-missing")


def test_list_events_invalid_limit_raises_value_error():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))

    with pytest.raises(ValueError, match="limit"):
        service.list_events("run-1", limit=0)


def test_append_log_records_per_stream_sequence():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))

    chunk1 = service.append_log("run-1", "stdout", ["line 1"])
    chunk2 = service.append_log("run-1", "stdout", ["line 2"])
    chunk3 = service.append_log("run-1", "stderr", ["err 1"])

    assert chunk1.sequence == 1
    assert chunk2.sequence == 2
    assert chunk3.sequence == 1


def test_list_logs_returns_stream_chunks():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))
    service.append_log("run-1", "stdout", ["a"])
    service.append_log("run-1", "stderr", ["b"])

    stdout = service.list_logs("run-1", "stdout")
    stderr = service.list_logs("run-1", "stderr")

    assert len(stdout) == 1
    assert stdout[0].lines == ("a",)
    assert len(stderr) == 1
    assert stderr[0].lines == ("b",)


def test_list_logs_with_after_and_limit():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))
    service.append_log("run-1", "stdout", ["a"])
    service.append_log("run-1", "stdout", ["b"])
    service.append_log("run-1", "stdout", ["c"])

    chunks = service.list_logs("run-1", "stdout", after="log-1", limit=1)

    assert len(chunks) == 1
    assert chunks[0].lines == ("b",)


def test_list_logs_unknown_cursor_raises_key_error():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))

    with pytest.raises(KeyError):
        service.list_logs("run-1", "stdout", after="log-missing")


def test_list_logs_invalid_limit_raises_value_error():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))

    with pytest.raises(ValueError, match="limit"):
        service.list_logs("run-1", "stdout", limit=-1)


def test_record_artifact_stores_in_insertion_order():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))
    artifact1 = RunArtifactRef(
        artifact_id="art-1",
        run_id="run-1",
        artifact_type="file",
        name="peaks.narrowPeak",
        uri="run://runs/run-1/artifacts/peaks.narrowPeak",
        mime_type=None,
        produced_at=datetime.now(timezone.utc),
        metadata={"sample": "S1"},
    )
    artifact2 = RunArtifactRef(
        artifact_id="art-2",
        run_id="run-1",
        artifact_type="file",
        name="summary.json",
        uri="run://runs/run-1/artifacts/summary.json",
        mime_type="application/json",
        produced_at=datetime.now(timezone.utc),
        metadata={},
    )

    service.record_artifact("run-1", artifact1)
    service.record_artifact("run-1", artifact2)

    artifacts = service.list_artifacts("run-1")
    assert artifacts == (artifact1, artifact2)


def test_record_artifact_rejects_mismatched_run_id():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))
    artifact = RunArtifactRef(
        artifact_id="art-1",
        run_id="run-other",
        artifact_type="file",
        name="peaks.narrowPeak",
        uri="run://runs/run-other/artifacts/peaks.narrowPeak",
        mime_type=None,
        produced_at=datetime.now(timezone.utc),
        metadata={},
    )

    with pytest.raises(ValueError, match="run_id"):
        service.record_artifact("run-1", artifact)


def test_record_artifact_rejects_duplicate_artifact_id():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))
    artifact = RunArtifactRef(
        artifact_id="art-1",
        run_id="run-1",
        artifact_type="file",
        name="peaks.narrowPeak",
        uri="run://runs/run-1/artifacts/peaks.narrowPeak",
        mime_type=None,
        produced_at=datetime.now(timezone.utc),
        metadata={},
    )
    service.record_artifact("run-1", artifact)

    with pytest.raises(ValueError, match="artifact_id"):
        service.record_artifact("run-1", artifact)


def test_recover_interrupted_runs_marks_only_active_states_as_failed():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    run_ids = iter(
        [
            "run-created",
            "run-planned",
            "run-validating",
            "run-queued",
            "run-running",
            "run-succeeded",
        ]
    )
    service = RunService(registry=registry, id_factory=lambda: next(run_ids))
    created = service.create_run("fake", WorkflowInputs(config={}))
    planned = service.create_run("fake", WorkflowInputs(config={}))
    validating = service.create_run("fake", WorkflowInputs(config={}))
    queued = service.create_run("fake", WorkflowInputs(config={}))
    running = service.create_run("fake", WorkflowInputs(config={}))
    succeeded = service.create_run("fake", WorkflowInputs(config={}))

    service.transition_run(planned.run_id, RunStatus.VALIDATING)
    service.transition_run(planned.run_id, RunStatus.PLANNED)
    service.transition_run(validating.run_id, RunStatus.VALIDATING, stage="preflight")
    service.transition_run(queued.run_id, RunStatus.VALIDATING)
    service.transition_run(queued.run_id, RunStatus.PLANNED)
    service.transition_run(queued.run_id, RunStatus.QUEUED, stage="queue")
    service.transition_run(running.run_id, RunStatus.VALIDATING)
    service.transition_run(running.run_id, RunStatus.PLANNED)
    service.transition_run(running.run_id, RunStatus.QUEUED)
    service.transition_run(running.run_id, RunStatus.RUNNING, stage="run")
    service.transition_run(succeeded.run_id, RunStatus.VALIDATING)
    service.transition_run(succeeded.run_id, RunStatus.PLANNED)
    service.transition_run(succeeded.run_id, RunStatus.QUEUED)
    service.transition_run(succeeded.run_id, RunStatus.RUNNING)
    service.transition_run(succeeded.run_id, RunStatus.SUCCEEDED)

    recovered = service.recover_interrupted_runs()

    assert [record.run_id for record in recovered] == [
        validating.run_id,
        queued.run_id,
        running.run_id,
    ]
    assert service.get_run(created.run_id).status is RunStatus.CREATED
    assert service.get_run(planned.run_id).status is RunStatus.PLANNED
    assert service.get_run(succeeded.run_id).status is RunStatus.SUCCEEDED
    for record in recovered:
        assert record.status is RunStatus.FAILED
        assert record.ended_at is not None
        assert record.error == Issue(
            code="RUN_INTERRUPTED_BY_API_RESTART",
            message="Run was interrupted by an API restart.",
            source="run_service",
            path="run_id",
            hint="Review the run events and submit a new preflight if needed.",
        )
        event = service.list_events(record.run_id)[-1]
        assert event.event_type == "run_recovered_after_restart"
        assert event.status is RunStatus.FAILED
        assert event.context["reason_code"] == "API_RESTART_INTERRUPTED"
        assert event.issue == record.error


def test_recover_interrupted_runs_is_idempotent_after_terminal_transition():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))
    service.transition_run("run-1", RunStatus.VALIDATING)

    first = service.recover_interrupted_runs()
    events_after_first = service.list_events("run-1")
    second = service.recover_interrupted_runs()

    assert len(first) == 1
    assert second == ()
    assert service.list_events("run-1") == events_after_first


SRC_ROOT = Path(__file__).resolve().parents[2] / "src"


def test_services_package_exports_run_service():
    code = """
        from encode_pipeline.services import RunService
        print(RunService.__name__)
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
    assert proc.stdout.strip() == "RunService"


def test_importing_run_service_does_not_import_upper_layers():
    code = """
        import sys
        import encode_pipeline.services.runs

        forbidden = [
            "encode_pipeline.api",
            "encode_pipeline.frontend",
            "fastapi",
            "pydantic",
            "sqlalchemy",
            "alembic",
            "snakemake",
            "openai",
        ]
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
