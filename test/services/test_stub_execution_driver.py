"""Tests for the synchronous stub execution driver."""

from __future__ import annotations

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
    assert record.current_stage is None
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
