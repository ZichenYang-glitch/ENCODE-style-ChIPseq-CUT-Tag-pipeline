"""Tests for the fail-closed local run driver skeleton."""

from __future__ import annotations

import ast
from datetime import datetime, timezone
from pathlib import Path
import sys

import pytest

from encode_pipeline.platform.adapters import (
    CommandSpec,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
    WorkspacePlan,
)
from encode_pipeline.platform.planning import ExecutionPlan, PlanStatus
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.services.command_builder import CommandBuilder
from encode_pipeline.services.local_run_driver import LocalRunDriver
from encode_pipeline.services.materialization import WorkspaceMaterializer
from encode_pipeline.services.process_runner import ProcessResult, ProcessRunner
from encode_pipeline.services.runs import RunService


class FakeAdapter:
    metadata = WorkflowMetadata(
        workflow_id="fake",
        name="Fake Workflow",
        version="1.0.0",
    )
    capabilities = WorkflowCapabilities(supports=("validation",))

    def schema(self) -> WorkflowSchema:
        return WorkflowSchema()

    def validate(self, inputs: WorkflowInputs) -> Result[object]:
        return Result.success(None)

    def preview_dag(self, inputs: WorkflowInputs) -> Result[object]:
        return Result.success(None)

    def plan_workspace(
        self,
        inputs: WorkflowInputs,
        workspace: str,
    ) -> Result[WorkspacePlan]:
        return Result.success(WorkspacePlan())

    def build_command(self, plan: WorkspacePlan) -> Result[CommandSpec]:
        return Result.success(CommandSpec(argv=("run-workflow",)))

    def extract_artifacts(self, inputs, workspace):
        return Result.success(())


# ---------------------------------------------------------------------------
# Shared test fixtures
# ---------------------------------------------------------------------------


def _make_materializer():
    return WorkspaceMaterializer()


def _make_registry():
    from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter

    return WorkflowRegistry(adapters=[EncodeStyleWorkflowAdapter()])


def _make_command_builder():
    return CommandBuilder(registry=_make_registry())


def _make_run_service() -> RunService:
    return RunService(
        registry=WorkflowRegistry(adapters=[FakeAdapter()]),
        id_factory=lambda: "run-1",
    )


def _make_driver(
    workspace_root: Path | None = None,
    materializer=None,
    command_builder=None,
) -> LocalRunDriver:
    if workspace_root is None:
        workspace_root = Path("/tmp/test-workspaces")
    if materializer is None:
        materializer = _make_materializer()
    if command_builder is None:
        command_builder = _make_command_builder()
    return LocalRunDriver(
        run_service=_make_run_service(),
        materializer=materializer,
        command_builder=command_builder,
        workspace_root=workspace_root,
    )


def _make_run_service_with_encode_adapter():
    """RunService wired to EncodeStyleWorkflowAdapter for happy-path tests."""
    return RunService(registry=_make_registry(), id_factory=lambda: "run-1")


def _make_pending_plan_with_workspace(
    run_id: str = "run-1",
    workflow_id: str = "encode-style-chipseq-cuttag-atac-mnase",
    workspace_plan: WorkspacePlan | None = None,
) -> ExecutionPlan:
    if workspace_plan is None:
        workspace_plan = WorkspacePlan(
            directories=("logs",),
            files=(("config/config.yaml", b"use_control: false\n"),),
        )
    return ExecutionPlan(
        plan_id="plan-1",
        run_id=run_id,
        workflow_id=workflow_id,
        status=PlanStatus.PENDING,
        inputs_snapshot={},
        workspace_plan=workspace_plan,
        created_at=datetime.now(timezone.utc),
    )


def _make_plan(
    run_id: str = "run-1",
    workflow_id: str = "fake",
    status: PlanStatus = PlanStatus.UNSUPPORTED,
    command_spec: CommandSpec | None = None,
) -> ExecutionPlan:
    return ExecutionPlan(
        plan_id="plan-1",
        run_id=run_id,
        workflow_id=workflow_id,
        status=status,
        inputs_snapshot={},
        command_spec=command_spec,
        created_at=datetime.now(timezone.utc),
    )


# ---------------------------------------------------------------------------
# Constructor validation tests
# ---------------------------------------------------------------------------


def test_local_run_driver_requires_run_service():
    with pytest.raises(
        ValueError, match="LocalRunDriver requires a RunService instance"
    ):
        LocalRunDriver(
            run_service="not-a-run-service",
            materializer=_make_materializer(),
            command_builder=_make_command_builder(),
            workspace_root=Path("/tmp"),
        )


def test_local_run_driver_requires_materializer():
    with pytest.raises(
        ValueError, match="LocalRunDriver requires a WorkspaceMaterializer instance"
    ):
        LocalRunDriver(
            run_service=_make_run_service(),
            materializer="not-a-materializer",
            command_builder=_make_command_builder(),
            workspace_root=Path("/tmp"),
        )


def test_local_run_driver_requires_command_builder():
    with pytest.raises(
        ValueError, match="LocalRunDriver requires a CommandBuilder instance"
    ):
        LocalRunDriver(
            run_service=_make_run_service(),
            materializer=_make_materializer(),
            command_builder="not-a-builder",
            workspace_root=Path("/tmp"),
        )


def test_local_run_driver_requires_absolute_workspace_root():
    with pytest.raises(
        ValueError, match="LocalRunDriver requires an absolute workspace_root Path"
    ):
        LocalRunDriver(
            run_service=_make_run_service(),
            materializer=_make_materializer(),
            command_builder=_make_command_builder(),
            workspace_root=Path("relative/path"),
        )


def test_local_run_driver_constructor_does_not_create_directories(tmp_path):
    workspace_root = tmp_path / "nonexistent"
    LocalRunDriver(
        run_service=_make_run_service(),
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=workspace_root,
    )
    assert not workspace_root.exists()


# ---------------------------------------------------------------------------
# Constructor — process_runner parameter
# ---------------------------------------------------------------------------


class _FakeProcessRunner(ProcessRunner):
    """Test double that returns controlled stdout/stderr/exit_code."""

    def __init__(
        self,
        *,
        exit_code: int = 0,
        stdout: str = "",
        stderr: str = "",
        issues: list[Issue] | None = None,
    ):
        super().__init__(allowed_executables=("snakemake",))
        self._exit_code = exit_code
        self._stdout = stdout
        self._stderr = stderr
        self._issues = issues or []

    def run(self, spec, *, output_callback=None):
        if output_callback is not None and self._stdout:
            output_callback("stdout", tuple(self._stdout.splitlines()))
        if output_callback is not None and self._stderr:
            output_callback("stderr", tuple(self._stderr.splitlines()))
        return Result.success(
            ProcessResult(
                exit_code=self._exit_code,
                stdout=self._stdout,
                stderr=self._stderr,
            ),
            issues=self._issues,
        )


def test_constructor_defaults_process_runner_to_process_runner_instance():
    driver = LocalRunDriver(
        run_service=_make_run_service(),
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=Path("/tmp/test-workspaces"),
    )
    assert isinstance(driver._process_runner, ProcessRunner)
    assert driver._process_runner._allowed_executables == ("snakemake",)


def test_constructor_accepts_explicit_process_runner():
    explicit = ProcessRunner(allowed_executables=("snakemake",), timeout_seconds=10.0)
    driver = LocalRunDriver(
        run_service=_make_run_service(),
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=Path("/tmp/test-workspaces"),
        process_runner=explicit,
    )
    assert driver._process_runner is explicit
    assert driver._process_runner._timeout_seconds == 10.0


def test_constructor_rejects_non_process_runner():
    with pytest.raises(ValueError, match="ProcessRunner"):
        LocalRunDriver(
            run_service=_make_run_service(),
            materializer=_make_materializer(),
            command_builder=_make_command_builder(),
            workspace_root=Path("/tmp/test-workspaces"),
            process_runner="not-a-runner",
        )


# ---------------------------------------------------------------------------
# Dry-run flag conflict tests
# ---------------------------------------------------------------------------


def test_run_refuses_when_argv_contains_dash_n(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
    )
    plan = _make_plan(
        status=PlanStatus.PLANNED,
        workflow_id="encode-style-chipseq-cuttag-atac-mnase",
        command_spec=CommandSpec(argv=("snakemake", "-n")),
    )

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    assert result.issues[0].code == "LOCAL_RUN_DRY_RUN_FLAG_CONFLICT"
    assert result.issues[0].source == "local_run_driver"
    assert result.issues[0].path == "command_spec"


def test_run_refuses_when_argv_contains_dash_dash_dry_run(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
    )
    plan = _make_plan(
        status=PlanStatus.PLANNED,
        workflow_id="encode-style-chipseq-cuttag-atac-mnase",
        command_spec=CommandSpec(argv=("snakemake", "--dry-run")),
    )

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    assert result.issues[0].code == "LOCAL_RUN_DRY_RUN_FLAG_CONFLICT"


def test_flag_conflict_records_runner_refused_with_correct_reason_code(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
    )
    plan = _make_plan(
        status=PlanStatus.PLANNED,
        workflow_id="encode-style-chipseq-cuttag-atac-mnase",
        command_spec=CommandSpec(argv=("snakemake", "-n")),
    )

    driver.run("run-1", plan)
    events = service.list_events("run-1")
    refuse_event = [e for e in events if e.event_type == "runner_refused"][0]

    assert refuse_event.context["reason_code"] == "LOCAL_RUN_DRY_RUN_FLAG_CONFLICT"
    # No dry_run events
    event_types = [e.event_type for e in events]
    assert "dry_run_completed" not in event_types
    assert "dry_run_failed" not in event_types


def test_flag_conflict_does_not_call_process_runner(tmp_path):
    """ProcessRunner.run() must not be called when flag conflict is detected."""
    calls = []

    class SpyProcessRunner(ProcessRunner):
        def __init__(self):
            super().__init__(allowed_executables=("snakemake",))

        def run(self, spec):
            calls.append(("run", spec))
            return Result.success(ProcessResult(exit_code=0, stdout="", stderr=""))

    spy = SpyProcessRunner()
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=spy,
    )
    plan = _make_plan(
        status=PlanStatus.PLANNED,
        workflow_id="encode-style-chipseq-cuttag-atac-mnase",
        command_spec=CommandSpec(argv=("snakemake", "-n")),
    )

    driver.run("run-1", plan)

    assert len(calls) == 0


# ---------------------------------------------------------------------------
# Dry-run exit 0 tests
# ---------------------------------------------------------------------------


def test_dry_run_exit_zero_returns_planned_plan(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunner(),
    )
    plan = _make_pending_plan_with_workspace()

    result = driver.run("run-1", plan)

    assert result.is_success is True
    assert result.value.status is PlanStatus.PLANNED
    assert result.value.command_spec is not None


def test_dry_run_exit_zero_records_dry_run_completed_event(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunner(),
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)
    events = service.list_events("run-1")
    dry_run_events = [e for e in events if e.event_type == "dry_run_completed"]

    assert len(dry_run_events) == 1
    assert dry_run_events[0].context == {"exit_code": 0}
    assert dry_run_events[0].status is None


def test_dry_run_exit_zero_does_not_record_dry_run_failed(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunner(),
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)
    events = service.list_events("run-1")
    event_types = [e.event_type for e in events]

    assert "dry_run_failed" not in event_types


def test_dry_run_exit_zero_original_plan_not_mutated(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunner(),
    )
    plan = _make_pending_plan_with_workspace()

    # PENDING plans start without a command_spec
    assert plan.command_spec is None
    driver.run("run-1", plan)
    # The original plan object is unchanged (command_spec still None)
    assert plan.command_spec is None
    # A new planned_plan was returned from _prepare — that has the command_spec,
    # but this original plan object is untouched


# ---------------------------------------------------------------------------
# Dry-run nonzero exit tests
# ---------------------------------------------------------------------------


class _FakeProcessRunnerNonzero(ProcessRunner):
    """Test double returning nonzero exit with a PROCESS_RUNNER_NONZERO_EXIT warning."""

    def __init__(self, stdout: str = "", stderr: str = ""):
        super().__init__(allowed_executables=("snakemake",))
        self._stdout = stdout
        self._stderr = stderr

    def run(self, spec):
        return Result.success(
            ProcessResult(exit_code=3, stdout=self._stdout, stderr=self._stderr),
            issues=[
                Issue(
                    code="PROCESS_RUNNER_NONZERO_EXIT",
                    message="Subprocess exited with a non-zero code.",
                    severity="warning",
                    path="command_spec",
                    source="process_runner",
                    context={"exit_code": 3},
                )
            ],
        )


def test_dry_run_nonzero_returns_failure_with_wrapper(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunnerNonzero(
            stdout="nonzero stdout",
            stderr="nonzero stderr",
        ),
    )
    plan = _make_pending_plan_with_workspace()

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    assert result.issues[0].code == "LOCAL_RUN_DRY_RUN_FAILED"
    assert result.issues[0].source == "local_run_driver"
    assert result.issues[0].severity == "error"
    assert len(result.issues) == 2
    assert result.issues[1].code == "PROCESS_RUNNER_NONZERO_EXIT"

    stdout_chunks = service.list_logs("run-1", "stdout")
    stderr_chunks = service.list_logs("run-1", "stderr")
    assert stdout_chunks[0].lines == ("nonzero stdout",)
    assert stderr_chunks[0].lines == ("nonzero stderr",)


def test_dry_run_nonzero_records_dry_run_failed_event(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunnerNonzero(),
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)
    events = service.list_events("run-1")
    dry_run_events = [e for e in events if e.event_type == "dry_run_failed"]

    assert len(dry_run_events) == 1
    assert dry_run_events[0].context["reason_code"] == "LOCAL_RUN_DRY_RUN_FAILED"
    assert dry_run_events[0].context["exit_code"] == 3
    assert dry_run_events[0].context["issue_count"] == 1


def test_dry_run_nonzero_runner_refused_reason_is_dry_run_failed(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunnerNonzero(),
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)
    events = service.list_events("run-1")
    refuse_event = [e for e in events if e.event_type == "runner_refused"][0]

    assert refuse_event.context["reason_code"] == "LOCAL_RUN_DRY_RUN_FAILED"


def test_dry_run_nonzero_does_not_record_dry_run_completed(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunnerNonzero(),
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)
    events = service.list_events("run-1")
    event_types = [e.event_type for e in events]

    assert "dry_run_completed" not in event_types


# ---------------------------------------------------------------------------
# ProcessRunner failure tests
# ---------------------------------------------------------------------------


class _FakeProcessRunnerTimeout(ProcessRunner):
    """Test double simulating a ProcessRunner timeout."""

    def __init__(self):
        super().__init__(allowed_executables=("snakemake",))

    def run(self, spec):
        return Result.failure(
            [
                Issue(
                    code="PROCESS_RUNNER_TIMEOUT",
                    message="Subprocess timed out.",
                    severity="error",
                    path="command_spec",
                    source="process_runner",
                    context={"timeout_seconds": 300.0},
                )
            ]
        )


class _FakeProcessRunnerNotFound(ProcessRunner):
    """Test double simulating executable not found."""

    def __init__(self):
        super().__init__(allowed_executables=("snakemake",))

    def run(self, spec):
        return Result.failure(
            [
                Issue(
                    code="PROCESS_RUNNER_EXECUTABLE_NOT_FOUND",
                    message="Executable not found.",
                    severity="error",
                    path="command_spec.argv[0]",
                    source="process_runner",
                )
            ]
        )


def test_process_runner_timeout_returns_failure_with_wrapper(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunnerTimeout(),
    )
    plan = _make_pending_plan_with_workspace()

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    assert result.issues[0].code == "LOCAL_RUN_PROCESS_FAILED"
    assert result.issues[0].source == "local_run_driver"
    assert result.issues[0].severity == "error"
    assert len(result.issues) == 2
    assert result.issues[1].code == "PROCESS_RUNNER_TIMEOUT"


def test_process_runner_timeout_records_dry_run_failed_event(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunnerTimeout(),
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)
    events = service.list_events("run-1")
    dry_run_events = [e for e in events if e.event_type == "dry_run_failed"]

    assert len(dry_run_events) == 1
    assert dry_run_events[0].context["reason_code"] == "LOCAL_RUN_PROCESS_FAILED"
    assert dry_run_events[0].context["process_issue_code"] == "PROCESS_RUNNER_TIMEOUT"
    assert dry_run_events[0].context["issue_count"] == 1


def test_process_runner_timeout_runner_refused_reason_is_process_failed(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunnerTimeout(),
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)
    events = service.list_events("run-1")
    refuse_event = [e for e in events if e.event_type == "runner_refused"][0]

    assert refuse_event.context["reason_code"] == "LOCAL_RUN_PROCESS_FAILED"


def test_process_runner_not_found_returns_failure_with_wrapper(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunnerNotFound(),
    )
    plan = _make_pending_plan_with_workspace()

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    assert result.issues[0].code == "LOCAL_RUN_PROCESS_FAILED"
    assert result.issues[1].code == "PROCESS_RUNNER_EXECUTABLE_NOT_FOUND"


def test_process_runner_not_found_records_dry_run_failed_event(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunnerNotFound(),
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)
    events = service.list_events("run-1")
    dry_run_events = [e for e in events if e.event_type == "dry_run_failed"]

    assert (
        dry_run_events[0].context["process_issue_code"]
        == "PROCESS_RUNNER_EXECUTABLE_NOT_FOUND"
    )


# ---------------------------------------------------------------------------
# Existing tests updated for new constructor
# ---------------------------------------------------------------------------


def test_local_run_driver_unknown_run_raises_key_error():
    driver = _make_driver()
    plan = _make_plan(run_id="missing")

    with pytest.raises(KeyError):
        driver.run("missing", plan)


def test_local_run_driver_rejects_run_id_mismatch():
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=Path("/tmp/test-workspaces"),
    )
    plan = _make_plan(run_id="other-run")

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    assert len(result.issues) == 1
    issue = result.issues[0]
    assert issue.code == "LOCAL_RUN_PLAN_MISMATCH"
    assert issue.severity.value == "error"
    assert issue.source == "local_run_driver"
    assert issue.path == "plan"


def test_local_run_driver_rejects_workflow_id_mismatch():
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=Path("/tmp/test-workspaces"),
    )
    plan = _make_plan(workflow_id="other-workflow")

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    assert result.issues[0].code == "LOCAL_RUN_PLAN_MISMATCH"


def test_local_run_driver_rejects_unsupported_plan():
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=Path("/tmp/test-workspaces"),
    )
    plan = _make_plan(status=PlanStatus.UNSUPPORTED)

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "LOCAL_RUN_NOT_PLANNED"
    assert issue.severity.value == "error"
    assert issue.source == "local_run_driver"
    assert issue.path == "plan"


def test_local_run_driver_rejects_pending_plan_without_workspace():
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=Path("/tmp/test-workspaces"),
    )
    plan = _make_plan(status=PlanStatus.PENDING)

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    assert result.issues[0].code == "LOCAL_RUN_MISSING_WORKSPACE_PLAN"


def test_local_run_driver_rejects_planned_plan_missing_command_spec():
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=Path("/tmp/test-workspaces"),
    )
    plan = _make_plan(status=PlanStatus.PLANNED)

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "LOCAL_RUN_MISSING_COMMAND_SPEC"
    assert issue.severity.value == "error"
    assert issue.source == "local_run_driver"
    assert issue.path == "plan"


def test_local_run_driver_executes_planned_command():
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=Path("/tmp/test-workspaces"),
        process_runner=_FakeProcessRunner(stdout="started", stderr="warning"),
    )
    plan = _make_plan(
        status=PlanStatus.PLANNED,
        command_spec=CommandSpec(argv=("snakemake", "--cores", "4")),
    )

    result = driver.run("run-1", plan)

    assert result.is_success is True
    assert result.value is plan
    assert service.list_logs("run-1", "stdout")[0].lines == ("started",)
    assert service.list_logs("run-1", "stderr")[0].lines == ("warning",)


def test_local_run_driver_bounds_persisted_stream_output():
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    runner = ProcessRunner(
        allowed_executables=(sys.executable,),
        max_output_bytes=64,
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=Path("/tmp/test-workspaces"),
        process_runner=runner,
    )
    plan = _make_plan(
        status=PlanStatus.PLANNED,
        command_spec=CommandSpec(
            argv=(sys.executable, "-c", "print('x' * 10000)"),
        ),
    )

    result = driver.run("run-1", plan)

    assert result.is_success
    persisted = "".join(
        line for chunk in service.list_logs("run-1", "stdout") for line in chunk.lines
    )
    assert len(persisted.encode("utf-8")) <= 64
    assert any(
        issue.code == "PROCESS_RUNNER_OUTPUT_LIMIT_EXCEEDED" for issue in result.issues
    )
    truncation_events = [
        event
        for event in service.list_events("run-1", limit=100)
        if event.event_type == "execution_output_truncated"
    ]
    assert len(truncation_events) == 1
    assert truncation_events[0].stage == "execution"
    assert truncation_events[0].context == {
        "reason_code": "PROCESS_RUNNER_OUTPUT_LIMIT_EXCEEDED"
    }
    assert truncation_events[0].issue is not None
    assert truncation_events[0].issue.code == "PROCESS_RUNNER_OUTPUT_LIMIT_EXCEEDED"


def test_local_run_driver_persists_truncation_before_timeout_failure():
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    runner = ProcessRunner(
        allowed_executables=(sys.executable,),
        timeout_seconds=0.1,
        max_output_bytes=8,
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=Path("/tmp/test-workspaces"),
        process_runner=runner,
    )
    plan = _make_plan(
        status=PlanStatus.PLANNED,
        command_spec=CommandSpec(
            argv=(
                sys.executable,
                "-c",
                "import time; print('x' * 1000, flush=True); time.sleep(10)",
            ),
        ),
    )

    result = driver.run("run-1", plan)

    assert result.is_failure
    assert [issue.code for issue in result.issues] == [
        "LOCAL_RUN_PROCESS_FAILED",
        "PROCESS_RUNNER_TIMEOUT",
        "PROCESS_RUNNER_OUTPUT_LIMIT_EXCEEDED",
    ]
    truncation_events = [
        event
        for event in service.list_events("run-1", limit=100)
        if event.event_type == "execution_output_truncated"
    ]
    assert len(truncation_events) == 1
    assert truncation_events[0].issue is not None
    assert truncation_events[0].issue.code == "PROCESS_RUNNER_OUTPUT_LIMIT_EXCEEDED"


def test_local_run_driver_records_runner_refused_event():
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=Path("/tmp/test-workspaces"),
    )
    plan = _make_plan(status=PlanStatus.UNSUPPORTED)

    before = service.list_events("run-1")
    result = driver.run("run-1", plan)
    after = service.list_events("run-1")

    assert result.is_failure is True
    assert len(after) == len(before) + 1
    event = after[-1]
    assert event.event_type == "runner_refused"
    assert event.status is None
    assert event.context == {
        "reason_code": "LOCAL_RUN_NOT_PLANNED",
        "plan_status": "unsupported",
        "can_execute": False,
        "has_command_spec": False,
    }


def test_local_run_driver_refusal_event_has_no_secrets_or_paths():
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={"secret": "value"}))
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=Path("/tmp/test-workspaces"),
    )
    plan = _make_plan(
        status=PlanStatus.PLANNED,
        command_spec=CommandSpec(argv=("snakemake", "-n"), cwd="/tmp/results"),
    )

    driver.run("run-1", plan)
    event = service.list_events("run-1")[-1]

    context_text = str(event.context)
    assert "snakemake" not in context_text
    assert "/tmp/results" not in context_text
    assert "secret" not in context_text
    assert "value" not in context_text


def test_local_run_driver_does_not_transition_run_status():
    service = _make_run_service()
    record = service.create_run("fake", WorkflowInputs(config={}))
    original_status = record.status
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=Path("/tmp/test-workspaces"),
    )
    plan = _make_plan(status=PlanStatus.PLANNED)

    before_events = service.list_events("run-1")
    driver.run("run-1", plan)
    after = service.get_run("run-1")
    after_events = service.list_events("run-1")

    assert after.status == original_status
    new_events = after_events[len(before_events) :]
    assert not any(event.event_type == "status_changed" for event in new_events)


def test_local_run_driver_does_not_mutate_run_service_logs_or_artifacts():
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=Path("/tmp/test-workspaces"),
    )
    plan = _make_plan(status=PlanStatus.PLANNED)

    driver.run("run-1", plan)

    assert service.list_logs("run-1", "stdout") == ()
    assert service.list_artifacts("run-1") == ()


def test_local_run_driver_unknown_run_does_not_record_event():
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=Path("/tmp/test-workspaces"),
    )
    plan = _make_plan(run_id="missing")

    with pytest.raises(KeyError):
        driver.run("missing", plan)

    assert len(service.list_events("run-1")) == 1  # only the create event


# ---------------------------------------------------------------------------
# Happy path tests — PENDING + workspace_plan → materialize → build → refuse
# ---------------------------------------------------------------------------


def test_run_materializes_and_builds_command_then_returns_planned(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    workspace_root = tmp_path / "workspaces"
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=workspace_root,
        process_runner=_FakeProcessRunner(),
    )
    plan = _make_pending_plan_with_workspace()

    result = driver.run("run-1", plan)

    assert result.is_success is True
    assert result.value.status is PlanStatus.PLANNED
    # Verify filesystem state
    config_file = workspace_root / "run-1" / "config" / "config.yaml"
    assert config_file.is_file()
    assert config_file.read_text() == "use_control: false\n"
    # Verify events
    events = service.list_events("run-1")
    event_types = [e.event_type for e in events]
    assert "workspace_materialized" in event_types
    assert "command_built" in event_types
    assert "dry_run_completed" in event_types
    assert "runner_refused" not in event_types


# ---------------------------------------------------------------------------
# Event structure tests
# ---------------------------------------------------------------------------


def test_workspace_materialized_event_context(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
    )
    plan = _make_pending_plan_with_workspace(
        workspace_plan=WorkspacePlan(
            directories=("logs", "results"),
            files=(("config/config.yaml", b"k: v\n"), ("samples.tsv", b"a\tb\n")),
        ),
    )

    driver.run("run-1", plan)
    events = service.list_events("run-1")
    mat_event = [e for e in events if e.event_type == "workspace_materialized"][0]

    assert mat_event.status is None
    assert mat_event.context == {"directory_count": 2, "file_count": 2}
    # No paths in context
    context_text = str(mat_event.context)
    assert "workspaces" not in context_text
    assert "config" not in context_text


def test_command_built_event_context(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)
    events = service.list_events("run-1")
    cmd_event = [e for e in events if e.event_type == "command_built"][0]

    assert cmd_event.status is None
    assert cmd_event.context == {"has_command_spec": True, "plan_status": "planned"}


def test_no_status_changed_events(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
    )
    plan = _make_pending_plan_with_workspace()

    before = service.list_events("run-1")
    driver.run("run-1", plan)
    after = service.list_events("run-1")

    new_events = after[len(before) :]
    assert not any(e.event_type == "status_changed" for e in new_events)


# ---------------------------------------------------------------------------
# Materialization failure test
# ---------------------------------------------------------------------------


def test_run_materialization_failure_refuses_and_includes_underlying_issues(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
    )
    # A workspace_plan with an absolute file path triggers materialization failure
    plan = _make_pending_plan_with_workspace(
        workspace_plan=WorkspacePlan(
            files=(("/etc/passwd", b"bad"),),
        ),
    )

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    # Wrapper is first
    assert result.issues[0].code == "LOCAL_RUN_MATERIALIZATION_FAILED"
    assert result.issues[0].source == "local_run_driver"
    # Underlying issue(s) follow
    assert len(result.issues) > 1
    # reason_code is the wrapper, not underlying
    events = service.list_events("run-1")
    refuse_event = [e for e in events if e.event_type == "runner_refused"][0]
    assert refuse_event.context["reason_code"] == "LOCAL_RUN_MATERIALIZATION_FAILED"
    # No workspace_materialized event
    assert not any(e.event_type == "workspace_materialized" for e in events)
    # No command_built event
    assert not any(e.event_type == "command_built" for e in events)


# ---------------------------------------------------------------------------
# Command build failure test
# ---------------------------------------------------------------------------


def test_run_command_build_failure_refuses_and_includes_underlying_issues(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
    )
    # Missing config/config.yaml in workspace_plan causes command build failure
    plan = _make_pending_plan_with_workspace(
        workspace_plan=WorkspacePlan(
            directories=("logs",),
            files=(("samples.tsv", b"a\tb\n"),),
        ),
    )

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    assert result.issues[0].code == "LOCAL_RUN_COMMAND_BUILD_FAILED"
    assert result.issues[0].source == "local_run_driver"
    assert len(result.issues) > 1
    # reason_code is the wrapper
    events = service.list_events("run-1")
    refuse_event = [e for e in events if e.event_type == "runner_refused"][0]
    assert refuse_event.context["reason_code"] == "LOCAL_RUN_COMMAND_BUILD_FAILED"
    # workspace_materialized IS recorded (materialization succeeded)
    assert any(e.event_type == "workspace_materialized" for e in events)
    # command_built is NOT recorded
    assert not any(e.event_type == "command_built" for e in events)


# ---------------------------------------------------------------------------
# PENDING without workspace_plan test
# ---------------------------------------------------------------------------


def test_run_pending_without_workspace_plan_refuses_with_missing_workspace_plan(
    tmp_path,
):
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
    )
    plan = _make_plan(status=PlanStatus.PENDING)

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    assert result.issues[0].code == "LOCAL_RUN_MISSING_WORKSPACE_PLAN"
    assert result.issues[0].source == "local_run_driver"
    events = service.list_events("run-1")
    assert events[-1].event_type == "runner_refused"
    assert events[-1].context["reason_code"] == "LOCAL_RUN_MISSING_WORKSPACE_PLAN"


# ---------------------------------------------------------------------------
# Workspace directory derivation failure tests
# ---------------------------------------------------------------------------


def test_run_workspace_dir_derivation_with_nul_byte_is_unknown_run(tmp_path):
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    driver = _make_driver(workspace_root=tmp_path / "workspaces")
    plan = _make_plan(run_id="bad\0run")

    # NUL in run_id means get_run won't find it -> KeyError, no events
    with pytest.raises(KeyError):
        driver.run("bad\0run", plan)

    assert len(service.list_events("run-1")) == 1


def test_run_workspace_dir_derivation_with_traversal_refuses(tmp_path):
    bad_run_id = ".."
    service = RunService(
        registry=WorkflowRegistry(adapters=[FakeAdapter()]),
        id_factory=lambda: bad_run_id,
    )
    service.create_run("fake", WorkflowInputs(config={}))
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
    )
    plan = _make_pending_plan_with_workspace(run_id=bad_run_id, workflow_id="fake")

    result = driver.run(bad_run_id, plan)

    assert result.is_failure is True
    assert result.issues[0].code == "LOCAL_RUN_WORKSPACE_DIR_INVALID"
    assert result.issues[0].path == "run_id"
    assert result.issues[0].source == "local_run_driver"
    events = service.list_events(bad_run_id)
    assert events[-1].event_type == "runner_refused"
    assert events[-1].context["reason_code"] == "LOCAL_RUN_WORKSPACE_DIR_INVALID"


# ---------------------------------------------------------------------------
# Safety tests — no path/command leakage in events
# ---------------------------------------------------------------------------


def test_run_materialization_failure_event_has_no_paths(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "my-secret-workspaces",
    )
    plan = _make_pending_plan_with_workspace(
        workspace_plan=WorkspacePlan(
            files=(("/etc/passwd", b"bad"),),
        ),
    )

    driver.run("run-1", plan)
    events = service.list_events("run-1")

    context_text = str([e.context for e in events])
    assert "my-secret-workspaces" not in context_text
    assert "/etc/passwd" not in context_text


def test_run_prepare_success_event_has_no_command_leakage(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)
    events = service.list_events("run-1")

    context_text = str([e.context for e in events])
    assert "snakemake" not in context_text
    assert "--cores" not in context_text
    assert "Snakefile" not in context_text


# ---------------------------------------------------------------------------
# Safety tests
# ---------------------------------------------------------------------------


def test_dry_run_events_contain_no_stdout_or_stderr(tmp_path):
    class _FakeProcessRunnerWithOutput(ProcessRunner):
        def __init__(self):
            super().__init__(allowed_executables=("snakemake",))

        def run(self, spec):
            return Result.success(
                ProcessResult(
                    exit_code=0,
                    stdout="SENSITIVE_OUTPUT_abc123",
                    stderr="SENSITIVE_ERROR_xyz789",
                )
            )

    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunnerWithOutput(),
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)
    events = service.list_events("run-1")

    context_text = str([e.context for e in events])
    assert "SENSITIVE_OUTPUT_abc123" not in context_text
    assert "SENSITIVE_ERROR_xyz789" not in context_text


def test_dry_run_events_contain_no_paths(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    workspace_root = tmp_path / "my-secret-workspaces"
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=workspace_root,
        process_runner=_FakeProcessRunner(),
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)
    events = service.list_events("run-1")

    context_text = str([e.context for e in events])
    assert "my-secret-workspaces" not in context_text
    assert "Snakefile" not in context_text
    assert "config.yaml" not in context_text


def test_dry_run_events_contain_no_command_strings(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunner(),
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)
    events = service.list_events("run-1")

    context_text = str([e.context for e in events])
    assert "snakemake" not in context_text
    assert "--cores" not in context_text
    assert "-n" not in context_text
    assert "--dry-run" not in context_text


def test_dry_run_does_not_transition_run_status(tmp_path):
    service = _make_run_service_with_encode_adapter()
    record = service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    original_status = record.status
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunner(),
    )
    plan = _make_pending_plan_with_workspace()

    before_events = service.list_events("run-1")
    driver.run("run-1", plan)
    after = service.get_run("run-1")
    after_events = service.list_events("run-1")

    assert after.status == original_status
    new_events = after_events[len(before_events) :]
    assert not any(e.event_type == "status_changed" for e in new_events)


def test_dry_run_empty_output_does_not_append_logs(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunner(),
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)

    assert service.list_logs("run-1", "stdout") == ()
    assert service.list_logs("run-1", "stderr") == ()


def test_dry_run_success_appends_stdout(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunner(stdout="line 1\nline 2"),
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)

    stdout_chunks = service.list_logs("run-1", "stdout")
    assert len(stdout_chunks) == 1
    assert stdout_chunks[0].lines == ("line 1", "line 2")
    assert service.list_logs("run-1", "stderr") == ()


def test_dry_run_success_appends_stderr(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunner(stderr="error a\nerror b"),
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)

    stderr_chunks = service.list_logs("run-1", "stderr")
    assert len(stderr_chunks) == 1
    assert stderr_chunks[0].lines == ("error a", "error b")
    assert service.list_logs("run-1", "stdout") == ()


def test_dry_run_success_appends_both_streams(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunner(
            stdout="line 1",
            stderr="error 1",
        ),
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)

    stdout_chunks = service.list_logs("run-1", "stdout")
    stderr_chunks = service.list_logs("run-1", "stderr")
    assert len(stdout_chunks) == 1
    assert len(stderr_chunks) == 1
    assert stdout_chunks[0].lines == ("line 1",)
    assert stderr_chunks[0].lines == ("error 1",)


def test_dry_run_trailing_newline_does_not_add_empty_line(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunner(stdout="line 1\n"),
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)

    stdout_chunks = service.list_logs("run-1", "stdout")
    assert stdout_chunks[0].lines == ("line 1",)


def test_dry_run_process_failure_does_not_fabricate_logs(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunnerTimeout(),
    )
    plan = _make_pending_plan_with_workspace()

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    assert result.issues[0].code == "LOCAL_RUN_PROCESS_FAILED"
    assert service.list_logs("run-1", "stdout") == ()
    assert service.list_logs("run-1", "stderr") == ()


@pytest.mark.parametrize(
    "process_runner,event_type",
    [
        pytest.param(
            _FakeProcessRunner(stdout="out", stderr="err"),
            "dry_run_completed",
            id="exit-0",
        ),
        pytest.param(
            _FakeProcessRunnerNonzero(stdout="out", stderr="err"),
            "dry_run_failed",
            id="nonzero-exit",
        ),
    ],
)
def test_dry_run_appends_logs_before_event(tmp_path, process_runner, event_type):
    class SpyRunService(RunService):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            self.calls = []

        def append_log(self, run_id, stream_name, lines):
            self.calls.append(("append_log", stream_name, tuple(lines)))
            return super().append_log(run_id, stream_name, lines)

        def add_event(self, *, run_id, event_type, **kwargs):
            self.calls.append(("add_event", event_type))
            return super().add_event(run_id=run_id, event_type=event_type, **kwargs)

    service = SpyRunService(
        registry=_make_registry(),
        id_factory=lambda: "run-1",
    )
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=process_runner,
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)

    event_indices = {call: i for i, call in enumerate(service.calls)}
    assert (
        event_indices[("append_log", "stdout", ("out",))]
        < event_indices[("add_event", event_type)]
    )
    assert (
        event_indices[("append_log", "stderr", ("err",))]
        < event_indices[("add_event", event_type)]
    )


def test_dry_run_single_chunk_per_stream_per_invocation(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunner(stdout="a\n\nb", stderr="c"),
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)

    stdout_chunks = service.list_logs("run-1", "stdout")
    stderr_chunks = service.list_logs("run-1", "stderr")
    assert len(stdout_chunks) == 1
    assert len(stderr_chunks) == 1
    assert stdout_chunks[0].lines == ("a", "", "b")
    assert stderr_chunks[0].lines == ("c",)


def test_dry_run_flag_conflict_events_contain_no_command_leakage(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={})
    )
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
    )
    plan = _make_plan(
        status=PlanStatus.PLANNED,
        command_spec=CommandSpec(
            argv=("snakemake", "-n", "--configfile", "/secret/path.yaml")
        ),
    )

    driver.run("run-1", plan)
    events = service.list_events("run-1")

    context_text = str([e.context for e in events])
    assert "/secret/path.yaml" not in context_text
    assert "--configfile" not in context_text


# ---------------------------------------------------------------------------
# Import boundary test
# ---------------------------------------------------------------------------


def test_local_run_driver_import_boundary() -> None:
    source_path = (
        Path(__file__).resolve().parents[2]
        / "src/encode_pipeline/services/local_run_driver.py"
    )
    source = source_path.read_text(encoding="utf-8")
    tree = ast.parse(source)

    forbidden_modules = {
        "encode_pipeline.api",
        "encode_pipeline.frontend",
        "snakemake",
        "subprocess",
    }
    imported_modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_modules.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imported_modules.add(node.module)

    violations = [
        module
        for module in imported_modules
        for forbidden in forbidden_modules
        if module == forbidden or module.startswith(f"{forbidden}.")
    ]
    assert not violations, f"Forbidden imports found: {violations}"

    # Verify required imports are present
    assert "encode_pipeline.services.process_runner" in imported_modules
    assert "encode_pipeline.platform.adapters" in imported_modules
