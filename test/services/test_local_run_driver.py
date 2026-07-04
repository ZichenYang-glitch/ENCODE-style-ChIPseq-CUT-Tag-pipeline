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
from encode_pipeline.services.local_run_driver import LocalRunDriver
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


def _make_run_service() -> RunService:
    return RunService(
        registry=WorkflowRegistry(adapters=[FakeAdapter()]),
        id_factory=lambda: "run-1",
    )


def _make_driver() -> LocalRunDriver:
    return LocalRunDriver(run_service=_make_run_service())


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


def test_local_run_driver_requires_run_service():
    with pytest.raises(ValueError, match="LocalRunDriver requires a RunService instance"):
        LocalRunDriver(run_service="not-a-run-service")


def test_local_run_driver_unknown_run_raises_key_error():
    driver = _make_driver()
    plan = _make_plan(run_id="missing")

    with pytest.raises(KeyError):
        driver.run("missing", plan)


def test_local_run_driver_rejects_run_id_mismatch():
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    driver = LocalRunDriver(run_service=service)
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
    driver = LocalRunDriver(run_service=service)
    plan = _make_plan(workflow_id="other-workflow")

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    assert result.issues[0].code == "LOCAL_RUN_PLAN_MISMATCH"


def test_local_run_driver_rejects_unsupported_plan():
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    driver = LocalRunDriver(run_service=service)
    plan = _make_plan(status=PlanStatus.UNSUPPORTED)

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "LOCAL_RUN_NOT_PLANNED"
    assert issue.severity.value == "error"
    assert issue.source == "local_run_driver"
    assert issue.path == "plan"


def test_local_run_driver_rejects_pending_plan():
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    driver = LocalRunDriver(run_service=service)
    plan = _make_plan(status=PlanStatus.PENDING)

    result = driver.run("run-1", plan)

    assert result.is_failure is True
    assert result.issues[0].code == "LOCAL_RUN_NOT_PLANNED"


def test_local_run_driver_rejects_planned_plan_missing_command_spec():
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    driver = LocalRunDriver(run_service=service)
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
    driver = LocalRunDriver(run_service=service)
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
    driver = LocalRunDriver(run_service=service)
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
    driver = LocalRunDriver(run_service=service)
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
    driver = LocalRunDriver(run_service=service)
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
    driver = LocalRunDriver(run_service=service)
    plan = _make_plan(status=PlanStatus.PLANNED)

    driver.run("run-1", plan)

    assert service.list_logs("run-1", "stdout") == ()
    assert service.list_artifacts("run-1") == ()


def test_local_run_driver_unknown_run_does_not_record_event():
    service = _make_run_service()
    service.create_run("fake", WorkflowInputs(config={}))
    driver = LocalRunDriver(run_service=service)
    plan = _make_plan(run_id="missing")

    with pytest.raises(KeyError):
        driver.run("missing", plan)

    assert len(service.list_events("run-1")) == 1  # only the create event


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
