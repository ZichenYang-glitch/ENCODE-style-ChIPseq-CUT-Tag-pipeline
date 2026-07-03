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
