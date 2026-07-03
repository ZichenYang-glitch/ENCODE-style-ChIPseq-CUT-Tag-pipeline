"""Tests for the synchronous stub execution driver."""

from __future__ import annotations

import ast
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
from encode_pipeline.platform.results import Result
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.stub_execution_driver import StubExecutionDriver


class FakeAdapter:
    """Minimal workflow adapter for stub driver tests."""

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
        workspace: str,
    ) -> Result[WorkspacePlan]:
        return Result.success(WorkspacePlan(directories=[workspace]))

    def build_command(self, plan: WorkspacePlan) -> Result[CommandSpec]:
        return Result.success(CommandSpec(argv=["run-workflow"]))


def test_advance_to_terminal_success_path():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))
    driver = StubExecutionDriver(service)

    record = driver.advance_to_terminal("run-1")

    assert record.status is RunStatus.SUCCEEDED
    assert record.current_stage == "run"
    assert record.started_at is not None
    assert record.ended_at is not None

    events = service.list_events("run-1")
    assert [event.status for event in events] == [
        RunStatus.CREATED,
        RunStatus.VALIDATING,
        RunStatus.PLANNED,
        RunStatus.QUEUED,
        RunStatus.RUNNING,
        RunStatus.SUCCEEDED,
    ]

    logs = service.list_logs("run-1", "stdout")
    assert len(logs) == 5
    assert logs[0].lines == ("[stub] validating inputs",)


def test_advance_to_terminal_failure_tag():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run(
        "fake",
        WorkflowInputs(config={}),
        tags={"stub_outcome": "failure"},
    )
    driver = StubExecutionDriver(service)

    record = driver.advance_to_terminal("run-1")

    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == "STUB_FAILURE"
    assert record.error.source == "stub_driver"
    assert record.ended_at is not None

    events = service.list_events("run-1")
    statuses = [event.status for event in events]
    assert statuses == [RunStatus.CREATED, RunStatus.VALIDATING, RunStatus.FAILED]

    logs = service.list_logs("run-1", "stdout")
    assert logs[-1].lines == ("[stub] marking failed state",)


def test_advance_to_terminal_is_idempotent_for_terminal_run():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))
    driver = StubExecutionDriver(service)

    driver.advance_to_terminal("run-1")
    before_events = service.list_events("run-1")
    before_logs = service.list_logs("run-1", "stdout")

    driver.advance_to_terminal("run-1")

    assert service.list_events("run-1") == before_events
    assert service.list_logs("run-1", "stdout") == before_logs


def test_advance_to_terminal_unknown_run_raises_key_error():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry)
    driver = StubExecutionDriver(service)

    with pytest.raises(KeyError):
        driver.advance_to_terminal("run-missing")


def test_stub_execution_driver_import_boundary() -> None:
    source_path = (
        Path(__file__).resolve().parents[2]
        / "src/encode_pipeline/services/stub_execution_driver.py"
    )
    source = source_path.read_text(encoding="utf-8")
    tree = ast.parse(source)

    forbidden_modules = {
        "encode_pipeline.api",
        "encode_pipeline.frontend",
        "snakemake",
        "subprocess",
        "openai",
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


def test_stub_output_contains_no_execution_like_wording():
    registry = WorkflowRegistry(adapters=[FakeAdapter()])
    service = RunService(registry=registry, id_factory=lambda: "run-1")
    service.create_run("fake", WorkflowInputs(config={}))
    driver = StubExecutionDriver(service)
    driver.advance_to_terminal("run-1")

    forbidden = [
        "execution",
        "execute",
        "command",
        "snakemake",
        "subprocess",
        "launch",
        "submit",
        "pipeline",
    ]
    text = ""
    for event in service.list_events("run-1"):
        text += event.message.lower()
    for chunk in service.list_logs("run-1", "stdout"):
        for line in chunk.lines:
            text += line.lower()

    for word in forbidden:
        assert word not in text, f"forbidden word {word!r} found in stub output"
