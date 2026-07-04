"""Tests for the fail-closed local run driver skeleton."""

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
from encode_pipeline.platform.runs import RunStatus
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

    def plan_workspace(self, inputs: WorkflowInputs, workspace: str) -> Result[WorkspacePlan]:
        return Result.success(WorkspacePlan())

    def build_command(self, plan: WorkspacePlan) -> Result[CommandSpec]:
        return Result.success(CommandSpec(argv=("run-workflow",)))


def _make_run_service() -> RunService:
    return RunService(registry=WorkflowRegistry(adapters=[FakeAdapter()]), id_factory=lambda: "run-1")


def _make_driver(run_service: RunService):
    from encode_pipeline.services.local_run_driver import LocalRunDriver

    return LocalRunDriver(run_service=run_service)


def _create_run(status: PlanStatus, *, command_spec: CommandSpec | None = None) -> tuple[RunService, ExecutionPlan]:
    service = _make_run_service()
    record = service.create_run(
        workflow_id="fake",
        inputs=WorkflowInputs(config={}),
    )
    plan = ExecutionPlan(
        plan_id="plan-1",
        run_id=record.run_id,
        workflow_id=record.workflow_id,
        status=status,
        inputs_snapshot={},
        command_spec=command_spec,
    )
    return service, plan


def test_local_run_driver_requires_run_service():
    from encode_pipeline.services.local_run_driver import LocalRunDriver

    with pytest.raises(ValueError, match="LocalRunDriver requires a RunService instance"):
        LocalRunDriver(run_service="not-a-run-service")


def test_local_run_driver_refuses_mismatched_run_id():
    service, plan = _create_run(PlanStatus.PLANNED, command_spec=CommandSpec(argv=("run",)))
    driver = _make_driver(service)

    mismatched_plan = ExecutionPlan(
        plan_id=plan.plan_id,
        run_id="different-run",
        workflow_id=plan.workflow_id,
        status=plan.status,
        inputs_snapshot=plan.inputs_snapshot,
        command_spec=plan.command_spec,
    )

    result = driver.run(run_id="run-1", plan=mismatched_plan)

    assert result.is_failure
    issue = result.errors[0]
    assert issue.code == "LOCAL_RUN_PLAN_MISMATCH"
    assert issue.severity.value == "error"
    assert issue.source == "local_run_driver"
    assert issue.path == "plan"


def test_local_run_driver_refuses_mismatched_workflow_id():
    service, plan = _create_run(PlanStatus.PLANNED, command_spec=CommandSpec(argv=("run",)))
    driver = _make_driver(service)

    mismatched_plan = ExecutionPlan(
        plan_id=plan.plan_id,
        run_id=plan.run_id,
        workflow_id="different-workflow",
        status=plan.status,
        inputs_snapshot=plan.inputs_snapshot,
        command_spec=plan.command_spec,
    )

    result = driver.run(run_id="run-1", plan=mismatched_plan)

    assert result.is_failure
    assert result.errors[0].code == "LOCAL_RUN_PLAN_MISMATCH"


def test_local_run_driver_refuses_unplanned_status():
    service, plan = _create_run(PlanStatus.PENDING)
    driver = _make_driver(service)

    result = driver.run(run_id="run-1", plan=plan)

    assert result.is_failure
    issue = result.errors[0]
    assert issue.code == "LOCAL_RUN_NOT_PLANNED"
    assert issue.severity.value == "error"
    assert issue.source == "local_run_driver"
    assert issue.path == "plan"


def test_local_run_driver_refuses_planned_without_command_spec():
    service, plan = _create_run(PlanStatus.PLANNED, command_spec=None)
    driver = _make_driver(service)

    result = driver.run(run_id="run-1", plan=plan)

    assert result.is_failure
    issue = result.errors[0]
    assert issue.code == "LOCAL_RUN_MISSING_COMMAND_SPEC"
    assert issue.severity.value == "error"
    assert issue.source == "local_run_driver"
    assert issue.path == "plan"


def test_local_run_driver_refuses_even_when_executable():
    service, plan = _create_run(PlanStatus.PLANNED, command_spec=CommandSpec(argv=("run",)))
    driver = _make_driver(service)

    assert plan.can_execute

    result = driver.run(run_id="run-1", plan=plan)

    assert result.is_failure
    issue = result.errors[0]
    assert issue.code == "LOCAL_RUN_NOT_IMPLEMENTED"
    assert issue.severity.value == "error"
    assert issue.source == "local_run_driver"
    assert issue.path == "plan"


def test_local_run_driver_records_runner_refused_event_with_safe_context():
    service, plan = _create_run(PlanStatus.PLANNED, command_spec=CommandSpec(argv=("run",)))
    driver = _make_driver(service)

    result = driver.run(run_id="run-1", plan=plan)

    assert result.is_failure
    issue = result.errors[0]
    events = service.list_events(run_id="run-1")
    refused_events = [event for event in events if event.event_type == "runner_refused"]
    assert len(refused_events) == 1

    event = refused_events[0]
    assert event.status is None
    assert event.message == "Local execution refused."
    assert event.context == {
        "reason_code": issue.code,
        "plan_status": plan.status.value,
        "can_execute": plan.can_execute,
        "has_command_spec": plan.command_spec is not None,
    }
    # No raw execution surface leaked into the event context.
    assert "argv" not in event.context
    assert "cwd" not in event.context
    assert "env" not in event.context
    assert "config" not in event.context
    assert "samples" not in event.context


def test_local_run_driver_does_not_transition_run_status():
    service, plan = _create_run(PlanStatus.PLANNED, command_spec=CommandSpec(argv=("run",)))
    driver = _make_driver(service)

    record_before = service.get_run(run_id="run-1")
    driver.run(run_id="run-1", plan=plan)
    record_after = service.get_run(run_id="run-1")

    assert record_after.status is record_before.status is RunStatus.CREATED
    assert record_after.updated_at == record_before.updated_at

    events = service.list_events(run_id="run-1")
    status_events = [event for event in events if event.event_type == "status_changed"]
    assert len(status_events) == 1  # only the initial "created" event
