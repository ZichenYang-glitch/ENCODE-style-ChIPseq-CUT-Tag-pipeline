"""Tests for the fail-closed local run driver skeleton."""

from __future__ import annotations

import ast
from datetime import datetime, timezone
from pathlib import Path

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
from encode_pipeline.platform.results import Result
from encode_pipeline.services.command_builder import CommandBuilder
from encode_pipeline.services.local_run_driver import LocalRunDriver
from encode_pipeline.services.materialization import WorkspaceMaterializer
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
    with pytest.raises(ValueError, match="LocalRunDriver requires a RunService instance"):
        LocalRunDriver(
            run_service="not-a-run-service",
            materializer=_make_materializer(),
            command_builder=_make_command_builder(),
            workspace_root=Path("/tmp"),
        )


def test_local_run_driver_requires_materializer():
    with pytest.raises(ValueError, match="LocalRunDriver requires a WorkspaceMaterializer instance"):
        LocalRunDriver(
            run_service=_make_run_service(),
            materializer="not-a-materializer",
            command_builder=_make_command_builder(),
            workspace_root=Path("/tmp"),
        )


def test_local_run_driver_requires_command_builder():
    with pytest.raises(ValueError, match="LocalRunDriver requires a CommandBuilder instance"):
        LocalRunDriver(
            run_service=_make_run_service(),
            materializer=_make_materializer(),
            command_builder="not-a-builder",
            workspace_root=Path("/tmp"),
        )


def test_local_run_driver_requires_absolute_workspace_root():
    with pytest.raises(ValueError, match="LocalRunDriver requires an absolute workspace_root Path"):
        LocalRunDriver(
            run_service=_make_run_service(),
            materializer=_make_materializer(),
            command_builder=_make_command_builder(),
            workspace_root=Path("relative/path"),
        )


def test_local_run_driver_constructor_does_not_create_directories(tmp_path):
    workspace_root = tmp_path / "nonexistent"
    driver = LocalRunDriver(
        run_service=_make_run_service(),
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=workspace_root,
    )
    assert not workspace_root.exists()


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


def test_local_run_driver_refuses_executable_plan_with_not_implemented():
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=Path("/tmp/test-workspaces"),
    )
    plan = _make_plan(
        status=PlanStatus.PLANNED,
        command_spec=CommandSpec(argv=("snakemake", "-n")),
    )

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "LOCAL_RUN_NOT_IMPLEMENTED"
    assert issue.severity.value == "error"
    assert issue.source == "local_run_driver"
    assert issue.path == "plan"


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
    new_events = after_events[len(before_events):]
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


def test_run_materializes_and_builds_command_then_refuses(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run("encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={}))
    workspace_root = tmp_path / "workspaces"
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=workspace_root,
    )
    plan = _make_pending_plan_with_workspace()

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    assert result.issues[0].code == "LOCAL_RUN_NOT_IMPLEMENTED"
    # Verify filesystem state
    config_file = workspace_root / "run-1" / "config" / "config.yaml"
    assert config_file.is_file()
    assert config_file.read_text() == "use_control: false\n"
    # Verify events
    events = service.list_events("run-1")
    event_types = [e.event_type for e in events]
    assert "workspace_materialized" in event_types
    assert "command_built" in event_types
    assert "runner_refused" in event_types
    # runner_refused must be last
    assert event_types[-1] == "runner_refused"


def test_run_runner_refused_context_after_successful_prepare(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run("encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={}))
    driver = LocalRunDriver(
        run_service=service,
        materializer=_make_materializer(),
        command_builder=_make_command_builder(),
        workspace_root=tmp_path / "workspaces",
    )
    plan = _make_pending_plan_with_workspace()

    driver.run("run-1", plan)
    events = service.list_events("run-1")
    refuse_event = [e for e in events if e.event_type == "runner_refused"][0]

    assert refuse_event.context["reason_code"] == "LOCAL_RUN_NOT_IMPLEMENTED"
    assert refuse_event.context["plan_status"] == "planned"
    assert refuse_event.context["can_execute"] is True
    assert refuse_event.context["has_command_spec"] is True


# ---------------------------------------------------------------------------
# Event structure tests
# ---------------------------------------------------------------------------


def test_workspace_materialized_event_context(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run("encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={}))
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
    service.create_run("encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={}))
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
    service.create_run("encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={}))
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

    new_events = after[len(before):]
    assert not any(e.event_type == "status_changed" for e in new_events)


# ---------------------------------------------------------------------------
# Materialization failure test
# ---------------------------------------------------------------------------


def test_run_materialization_failure_refuses_and_includes_underlying_issues(tmp_path):
    service = _make_run_service_with_encode_adapter()
    service.create_run("encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={}))
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
    service.create_run("encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={}))
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


def test_run_pending_without_workspace_plan_refuses_with_missing_workspace_plan(tmp_path):
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
    service.create_run("encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={}))
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
    service.create_run("encode-style-chipseq-cuttag-atac-mnase", WorkflowInputs(config={}))
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

    assert not any(
        module == forbidden or module.startswith(f"{forbidden}.")
        for module in imported_modules
        for forbidden in forbidden_modules
    )
