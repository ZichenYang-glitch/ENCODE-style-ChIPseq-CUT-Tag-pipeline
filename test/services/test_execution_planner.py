"""Tests for the ExecutionPlanner service boundary."""

import pytest

from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.runs import RunService


@pytest.fixture
def fake_adapter():
    from encode_pipeline.platform.adapters import (
        CommandSpec,
        DagPreview,
        WorkflowAdapter,
        WorkflowCapabilities,
        WorkflowMetadata,
        WorkflowSchema,
        WorkspacePlan,
    )
    from encode_pipeline.platform.results import Result

    class _FakeAdapter:
        metadata = WorkflowMetadata(
            workflow_id="stub",
            name="Stub Workflow",
            version="0.0.1",
        )
        capabilities = WorkflowCapabilities(supports=())

        def schema(self) -> WorkflowSchema:
            return WorkflowSchema()

        def validate(self, inputs: WorkflowInputs) -> Result[object]:
            return Result.success(None)

        def preview_dag(self, inputs: WorkflowInputs) -> Result[DagPreview]:
            raise AssertionError("ExecutionPlanner must not call preview_dag")

        def plan_workspace(
            self,
            inputs: WorkflowInputs,
            workspace: str,
        ) -> Result[WorkspacePlan]:
            raise AssertionError("ExecutionPlanner must not call plan_workspace")

        def build_command(self, plan: WorkspacePlan) -> Result[CommandSpec]:
            raise AssertionError("ExecutionPlanner must not call build_command")

    return _FakeAdapter()


@pytest.fixture
def run_service(fake_adapter):
    return RunService(registry=WorkflowRegistry(adapters=[fake_adapter]))


@pytest.fixture
def planner(run_service):
    from encode_pipeline.services.planning import ExecutionPlanner

    return ExecutionPlanner(run_service=run_service)


def test_plan_run_unknown_run_returns_failure():
    from encode_pipeline.services.planning import ExecutionPlanner

    run_service = RunService(registry=WorkflowRegistry(adapters=[]))
    planner = ExecutionPlanner(run_service=run_service)
    result = planner.plan_run("does-not-exist")

    assert result.is_failure is True
    assert result.value is None
    assert len(result.issues) == 1
    issue = result.issues[0]
    assert issue.code == "EXECUTION_RUN_NOT_FOUND"
    assert issue.message == "Run not found."
    assert issue.severity.value == "error"
    assert issue.path == "does-not-exist"
    assert issue.source == "execution_planner"


def test_plan_run_valid_run_returns_unsupported_plan(planner, run_service):
    inputs = WorkflowInputs(config={"genome": "hg38"}, samples=None, options={})
    record = run_service.create_run("stub", inputs)

    result = planner.plan_run(record.run_id)

    assert result.is_success is True
    plan = result.value
    assert plan is not None
    assert plan.run_id == record.run_id
    assert plan.workflow_id == record.workflow_id
    assert plan.status.value == "unsupported"
    assert plan.dag_preview is None
    assert plan.workspace_plan is None
    assert plan.command_spec is None
    assert plan.can_execute is False
    assert len(result.issues) == 1
    issue = result.issues[0]
    assert issue.code == "EXECUTION_PLANNING_UNSUPPORTED"
    assert issue.message == "Execution planning is not supported yet."
    assert issue.severity.value == "info"
    assert issue.path == "execution_plan"
    assert issue.source == "execution_planner"


def test_plan_run_does_not_change_run_status(planner, run_service):
    inputs = WorkflowInputs(config={}, samples=None, options={})
    record = run_service.create_run("stub", inputs)
    original_status = record.status
    original_updated_at = record.updated_at

    planner.plan_run(record.run_id)

    after = run_service.get_run(record.run_id)
    assert after.status == original_status
    assert after.updated_at == original_updated_at


def test_plan_run_does_not_add_events(planner, run_service):
    inputs = WorkflowInputs(config={}, samples=None, options={})
    record = run_service.create_run("stub", inputs)
    before_count = len(run_service.list_events(record.run_id))

    planner.plan_run(record.run_id)

    after_count = len(run_service.list_events(record.run_id))
    assert after_count == before_count


def test_plan_run_does_not_append_logs(planner, run_service):
    inputs = WorkflowInputs(config={}, samples=None, options={})
    record = run_service.create_run("stub", inputs)

    planner.plan_run(record.run_id)

    chunks = run_service.list_logs(record.run_id, "stdout")
    assert chunks == ()


def test_plan_run_does_not_record_artifacts(planner, run_service):
    inputs = WorkflowInputs(config={}, samples=None, options={})
    record = run_service.create_run("stub", inputs)

    planner.plan_run(record.run_id)

    artifacts = run_service.list_artifacts(record.run_id)
    assert artifacts == ()


def test_plan_run_does_not_call_adapter_planning_methods(planner, run_service):
    inputs = WorkflowInputs(config={}, samples=None, options={})
    record = run_service.create_run("stub", inputs)

    result = planner.plan_run(record.run_id)

    assert result.is_success is True
    plan = result.value
    assert plan.status.value == "unsupported"
    assert plan.dag_preview is None
    assert plan.workspace_plan is None
    assert plan.command_spec is None


def test_plan_run_defensively_copies_inputs(planner, run_service):
    inputs = WorkflowInputs(
        config={"samples": ["s1"]},
        samples=None,
        options={},
    )
    record = run_service.create_run("stub", inputs)

    result = planner.plan_run(record.run_id)
    plan = result.value

    original_inputs = run_service.get_run(record.run_id).inputs
    original_inputs["config"]["samples"].append("s2")

    assert plan.inputs_snapshot["config"] == {"samples": ["s1"]}


def test_default_factory_returns_planner(run_service):
    from encode_pipeline.services import ExecutionPlanner, create_default_execution_planner

    planner = create_default_execution_planner(run_service=run_service)
    assert isinstance(planner, ExecutionPlanner)


from pathlib import Path

from encode_pipeline.platform.adapters import WorkspacePlan
from encode_pipeline.platform.planning import ExecutionPlan, PlanStatus
from encode_pipeline.services.planning import ExecutionPlanner, WorkspacePlanner


def _make_execution_plan(run_service):
    inputs = WorkflowInputs(config={}, samples=None, options={})
    record = run_service.create_run("stub", inputs)
    planner = ExecutionPlanner(run_service=run_service)
    return planner.plan_run(record.run_id).value


def test_plan_workspace_returns_pending_plan_with_workspace_plan(run_service, tmp_path):
    input_plan = _make_execution_plan(run_service)
    workspace_planner = WorkspacePlanner()
    base_dir = tmp_path.resolve()

    result = workspace_planner.plan_workspace(input_plan, base_dir=base_dir)

    assert result.is_success is True
    plan = result.value
    assert plan is not None
    assert plan is not input_plan
    assert plan.status is PlanStatus.PENDING
    assert plan.workspace_plan is not None
    assert plan.workspace_plan.directories == ("logs", "results")
    assert plan.workspace_plan.files == ()
    assert plan.command_spec is None
    assert plan.can_execute is False
    assert plan.run_id == input_plan.run_id
    assert plan.workflow_id == input_plan.workflow_id
    assert plan.dag_preview == input_plan.dag_preview
    assert plan.inputs_snapshot == input_plan.inputs_snapshot


def test_plan_workspace_does_not_create_directories_or_files(run_service, tmp_path):
    input_plan = _make_execution_plan(run_service)
    workspace_planner = WorkspacePlanner()
    base_dir = tmp_path.resolve()

    workspace_planner.plan_workspace(input_plan, base_dir=base_dir)

    assert not (base_dir / "logs").exists()
    assert not (base_dir / "results").exists()


def test_plan_workspace_does_not_mutate_input_plan(run_service, tmp_path):
    input_plan = _make_execution_plan(run_service)
    workspace_planner = WorkspacePlanner()
    base_dir = tmp_path.resolve()

    workspace_planner.plan_workspace(input_plan, base_dir=base_dir)

    assert input_plan.status is PlanStatus.UNSUPPORTED
    assert input_plan.workspace_plan is None
    assert input_plan.command_spec is None


def test_plan_workspace_rejects_relative_base_dir(run_service, tmp_path):
    input_plan = _make_execution_plan(run_service)
    workspace_planner = WorkspacePlanner()

    result = workspace_planner.plan_workspace(input_plan, base_dir=Path("relative/path"))

    assert result.is_failure is True
    assert result.value is None
    assert len(result.issues) == 1
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_BASE_DIR_RELATIVE"
    assert issue.path == "base_dir"


def test_plan_workspace_does_not_mutate_run_service(run_service, tmp_path):
    input_plan = _make_execution_plan(run_service)
    before_record = run_service.get_run(input_plan.run_id)
    before_events = run_service.list_events(input_plan.run_id)
    before_logs = run_service.list_logs(input_plan.run_id, "stdout")
    before_artifacts = run_service.list_artifacts(input_plan.run_id)

    workspace_planner = WorkspacePlanner()
    workspace_planner.plan_workspace(input_plan, base_dir=tmp_path.resolve())

    after_record = run_service.get_run(input_plan.run_id)
    assert after_record.status == before_record.status
    assert after_record.updated_at == before_record.updated_at
    assert run_service.list_events(input_plan.run_id) == before_events
    assert run_service.list_logs(input_plan.run_id, "stdout") == before_logs
    assert run_service.list_artifacts(input_plan.run_id) == before_artifacts


def test_plan_workspace_does_not_call_adapter_planning_methods(run_service, tmp_path):
    # fake_adapter fixture already raises on planning methods.
    input_plan = _make_execution_plan(run_service)
    workspace_planner = WorkspacePlanner()

    result = workspace_planner.plan_workspace(input_plan, base_dir=tmp_path.resolve())

    assert result.is_success is True
    assert result.value.workspace_plan is not None
