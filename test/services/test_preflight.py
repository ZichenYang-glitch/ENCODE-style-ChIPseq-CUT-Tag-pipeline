"""Tests for the local preflight orchestrator."""

from __future__ import annotations

from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path

from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagPreview,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
    WorkspacePlan,
)
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Result
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.command_builder import CommandBuilder
from encode_pipeline.services.defaults import (
    create_default_workflow_build_identity_provider,
    create_default_workspace_planner,
)
from encode_pipeline.services.local_run_driver import LocalRunDriver
from encode_pipeline.services.materialization import WorkspaceMaterializer
from encode_pipeline.services.planning import ExecutionPlanner
from encode_pipeline.services.preflight import LocalPreflightService
from encode_pipeline.services.process_runner import ProcessResult, ProcessRunner
from encode_pipeline.services.runs import RunService


def _make_registry():
    return WorkflowRegistry(adapters=[EncodeStyleWorkflowAdapter()])


def _make_run_service(registry=None):
    return RunService(registry=registry or _make_registry())


def _make_command_builder(registry=None):
    return CommandBuilder(registry=registry or _make_registry())


def _make_valid_inputs(tmp_path: Path) -> WorkflowInputs:
    samples_tsv = tmp_path / "samples.tsv"
    samples_tsv.write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
        "S1\t/abs/S1_1.fq.gz\t/abs/S1_2.fq.gz\tPE\tchipseq\tH3K27ac\tnarrow\ths\t/abs/bt2/GRCh38\n",
        encoding="utf-8",
    )
    config = {
        "samples": str(samples_tsv),
        "threads": 1,
        "genome_resources": {"hs": {"effective_genome_size": "hs"}},
    }
    return WorkflowInputs(config=config, samples=str(samples_tsv), options={})


class _FakeProcessRunner(ProcessRunner):
    def __init__(self, *, exit_code: int = 0, stdout: str = "", stderr: str = ""):
        super().__init__(allowed_executables=("snakemake",))
        self._exit_code = exit_code
        self._stdout = stdout
        self._stderr = stderr

    def run(self, spec):
        return Result.success(
            ProcessResult(
                exit_code=self._exit_code,
                stdout=self._stdout,
                stderr=self._stderr,
            )
        )


class _ConfigurationPreflightAdapter:
    metadata = WorkflowMetadata(
        workflow_id="configuration-preflight",
        name="Configuration preflight",
        version="1.0.0",
        engines=("opaque-engine",),
    )
    capabilities = WorkflowCapabilities(
        supports=("validation", "workspace_plan", "command", "input_authoring")
    )

    def schema(self):
        return WorkflowSchema()

    def validate(self, inputs):
        return Result.success({})

    def preview_dag(self, inputs):
        return Result.success(DagPreview())

    def plan_workspace(self, inputs, workspace):
        return Result.success(WorkspacePlan(directories=("logs",)))

    def build_command(self, plan, workspace):
        return Result.success(
            CommandSpec(
                argv=("pinned-engine", "run"),
                cwd=str(workspace),
                preflight_argv=("pinned-engine", "config"),
                preflight_kind="configuration",
            )
        )

    def extract_artifacts(self, inputs, workspace):
        return Result.success(())


class _CallbackProcessRunner(ProcessRunner):
    def __init__(self, callback: Callable[[], None]):
        super().__init__(allowed_executables=("snakemake",))
        self._callback = callback

    def run(self, spec):
        self._callback()
        return Result.success(ProcessResult(exit_code=0, stdout="", stderr=""))


def _make_service(
    tmp_path: Path,
    *,
    runner: ProcessRunner | None = None,
    registry=None,
    run_service=None,
    build_identity_provider=None,
):
    registry = registry or _make_registry()
    run_service = run_service or _make_run_service(registry=registry)
    workspace_root = tmp_path / "workspaces"
    local_run_driver = LocalRunDriver(
        run_service=run_service,
        materializer=WorkspaceMaterializer(),
        command_builder=_make_command_builder(registry=registry),
        workspace_root=workspace_root,
        process_runner=runner or _FakeProcessRunner(),
    )
    build_identity_provider = (
        build_identity_provider
        or create_default_workflow_build_identity_provider(registry=registry)
    )
    return LocalPreflightService(
        run_service=run_service,
        execution_planner=ExecutionPlanner(run_service=run_service),
        workspace_planner=create_default_workspace_planner(registry=registry),
        local_run_driver=local_run_driver,
        build_identity_provider=build_identity_provider,
    )


def _build_identity(digest: str) -> WorkflowBuildIdentity:
    return WorkflowBuildIdentity(
        workflow_id="encode-style-chipseq-cuttag-atac-mnase",
        adapter_version="1.0.0",
        scheme="sha256-tree-v1",
        logical_entrypoint="workflow/Snakefile",
        digest=digest,
        captured_at=datetime.now(timezone.utc),
    )


class _SequencedBuildIdentityProvider:
    def __init__(self, *identities: WorkflowBuildIdentity) -> None:
        self._identities = iter(identities)

    def capture(self, _workflow_id: str):
        return Result.success(next(self._identities))


def test_preflight_transitions_run_to_planned(tmp_path):
    service = _make_service(tmp_path)
    inputs = _make_valid_inputs(tmp_path)
    record = service._run_service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", inputs
    )

    result = service.preflight(record.run_id)

    assert result.is_success is True
    assert result.value.status is RunStatus.PLANNED
    events = service._run_service.list_events(record.run_id)
    event_types = [e.event_type for e in events]
    assert "status_changed" in event_types
    assert "workspace_materialized" in event_types
    assert "command_built" in event_types
    assert "dry_run_completed" in event_types
    assert "preflight_completed" in event_types
    completed = [event for event in events if event.event_type == "preflight_completed"]
    assert len(completed) == 1
    assert completed[0].status is RunStatus.PLANNED
    assert completed[0].context["new_status"] == RunStatus.PLANNED.value
    assert completed[0].context["workflow_build_scheme"] == "sha256-tree-v1"
    assert len(completed[0].context["workflow_build_digest"]) == 64
    assert completed[0].message == "Local preflight completed; dry-run succeeded."
    assert completed[0].context["reason_code"] == "PREFLIGHT_COMPLETED"
    assert "preflight_kind" not in completed[0].context
    assert events[-1] == completed[0]
    identity = service._run_service.get_workflow_build_identity(record.run_id)
    assert identity is not None
    assert identity.digest == completed[0].context["workflow_build_digest"]


def test_configuration_preflight_uses_engine_neutral_completion_message(tmp_path):
    adapter = _ConfigurationPreflightAdapter()
    registry = WorkflowRegistry([adapter])
    service = _make_service(tmp_path, registry=registry)
    record = service._run_service.create_run(
        adapter.metadata.workflow_id,
        WorkflowInputs(config={}),
    )

    result = service.preflight(record.run_id)

    assert result.is_success
    completed = [
        event
        for event in service._run_service.list_events(record.run_id)
        if event.event_type == "preflight_completed"
    ]
    assert len(completed) == 2
    command_event = next(event for event in completed if event.status is None)
    lifecycle_event = next(
        event for event in completed if event.status is RunStatus.PLANNED
    )
    assert command_event.message == (
        "Workflow configuration preflight completed successfully."
    )
    assert lifecycle_event.message == (
        "Local workflow configuration preflight completed successfully."
    )
    assert lifecycle_event.context["preflight_kind"] == "configuration"
    assert lifecycle_event.context["reason_code"] == (
        "PREFLIGHT_CONFIGURATION_COMPLETED"
    )
    assert all("dry-run" not in event.message for event in completed)


def test_preflight_refuses_duplicate_trigger(tmp_path):
    service = _make_service(tmp_path)
    inputs = _make_valid_inputs(tmp_path)
    record = service._run_service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", inputs
    )
    service.preflight(record.run_id)

    result = service.preflight(record.run_id)

    assert result.is_failure is True
    assert result.issues[0].code == "PREFLIGHT_ALREADY_TRIGGERED"


def test_preflight_fails_when_dry_run_exits_nonzero(tmp_path):
    runner = _FakeProcessRunner(exit_code=1, stderr="dry-run error")
    service = _make_service(tmp_path, runner=runner)
    inputs = _make_valid_inputs(tmp_path)
    record = service._run_service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", inputs
    )

    result = service.preflight(record.run_id)

    assert result.is_failure is True
    updated = service._run_service.get_run(record.run_id)
    assert updated.status is RunStatus.FAILED
    assert updated.error.code == "PREFLIGHT_FAILED"
    assert updated.error.context == {"reason_code": "LOCAL_RUN_DRY_RUN_FAILED"}

    failure_events = [
        e
        for e in service._run_service.list_events(record.run_id)
        if e.event_type == "status_changed" and e.context.get("new_status") == "failed"
    ]
    assert failure_events
    assert failure_events[0].context["reason_code"] == "LOCAL_RUN_DRY_RUN_FAILED"

    logs = service._run_service.list_logs(record.run_id, "stderr")
    assert any("dry-run error" in line for chunk in logs for line in chunk.lines)


def test_preflight_fails_closed_when_workflow_source_changes_during_dry_run(
    tmp_path,
):
    provider = _SequencedBuildIdentityProvider(
        _build_identity("a" * 64),
        _build_identity("b" * 64),
    )
    service = _make_service(tmp_path, build_identity_provider=provider)
    record = service._run_service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase",
        _make_valid_inputs(tmp_path),
    )

    result = service.preflight(record.run_id)

    assert result.is_failure is True
    assert result.issues[0].code == "PREFLIGHT_WORKFLOW_BUILD_CHANGED"
    failed = service._run_service.get_run(record.run_id)
    assert failed.status is RunStatus.FAILED
    assert failed.error is not None
    assert failed.error.context == {"reason_code": "PREFLIGHT_WORKFLOW_BUILD_CHANGED"}
    assert service._run_service.get_workflow_build_identity(record.run_id) is None


def test_preflight_fails_when_workspace_collides(tmp_path):
    service = _make_service(tmp_path)
    inputs = _make_valid_inputs(tmp_path)
    record = service._run_service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", inputs
    )
    workspace_dir = service._local_run_driver.derive_workspace_dir(record.run_id)
    existing = workspace_dir / "config" / "config.yaml"
    existing.parent.mkdir(parents=True)
    existing.write_bytes(b"exists")

    result = service.preflight(record.run_id)

    assert result.is_failure is True
    assert service._run_service.get_run(record.run_id).status is RunStatus.FAILED


def test_run_preflight_expects_validating_status(tmp_path):
    service = _make_service(tmp_path)
    inputs = _make_valid_inputs(tmp_path)
    record = service._run_service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", inputs
    )

    result = service.run_preflight(record.run_id)

    assert result.is_failure is True
    assert result.issues[0].code == "PREFLIGHT_ALREADY_TRIGGERED"

    service._run_service.transition_run(
        record.run_id, RunStatus.VALIDATING, stage="preflight"
    )
    result = service.run_preflight(record.run_id)

    assert result.is_success is True
    assert result.value.status is RunStatus.PLANNED


def test_preflight_respects_cancellation_during_dry_run(tmp_path):
    registry = _make_registry()
    run_service = _make_run_service(registry=registry)
    cancelled = [False]

    def cancel_during_dry_run():
        cancelled[0] = True
        run_service.cancel_run(record.run_id, reason="User cancelled during preflight.")

    runner = _CallbackProcessRunner(callback=cancel_during_dry_run)
    service = _make_service(
        tmp_path, registry=registry, run_service=run_service, runner=runner
    )
    inputs = _make_valid_inputs(tmp_path)
    record = run_service.create_run("encode-style-chipseq-cuttag-atac-mnase", inputs)
    run_service.transition_run(record.run_id, RunStatus.VALIDATING, stage="preflight")

    result = service.run_preflight(record.run_id)

    assert cancelled[0] is True
    assert result.is_failure is True
    assert result.issues[0].code == "PREFLIGHT_CANCELLED"
    updated = run_service.get_run(record.run_id)
    assert updated.status is RunStatus.CANCELLED

    status_events = [
        e
        for e in run_service.list_events(record.run_id)
        if e.event_type == "status_changed"
    ]
    assert status_events[-1].context["new_status"] == "cancelled"
    assert not any(
        e.context["new_status"] in ("planned", "failed") for e in status_events
    )


def test_run_preflight_returns_cancelled_when_cancelled_before_worker_start(tmp_path):
    """Cancellation before the background worker starts must preserve CANCELLED
    and must not begin workspace materialization, command build, or dry-run.
    """
    service = _make_service(tmp_path)
    inputs = _make_valid_inputs(tmp_path)
    record = service._run_service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", inputs
    )
    service._run_service.transition_run(
        record.run_id, RunStatus.VALIDATING, stage="preflight"
    )
    service._run_service.cancel_run(record.run_id, reason="Cancelled before worker.")

    result = service.run_preflight(record.run_id)

    assert result.is_failure is True
    assert result.issues[0].code == "PREFLIGHT_CANCELLED"
    updated = service._run_service.get_run(record.run_id)
    assert updated.status is RunStatus.CANCELLED

    event_types = [
        e.event_type for e in service._run_service.list_events(record.run_id)
    ]
    assert "workspace_materialized" not in event_types
    assert "command_built" not in event_types
    assert "dry_run_completed" not in event_types
    assert "preflight_completed" not in event_types
