"""Tests for the workflow info service."""

from __future__ import annotations

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
from encode_pipeline.services.workflow_info import WorkflowInfoService


class FakeAdapter:
    """Minimal workflow adapter for service-layer tests."""

    def __init__(
        self,
        workflow_id: str,
        name: str = "Fake Workflow",
        version: str = "1.0.0",
        supports: tuple[str, ...] = ("validation",),
    ) -> None:
        self.metadata = WorkflowMetadata(
            workflow_id=workflow_id,
            name=name,
            version=version,
        )
        self.capabilities = WorkflowCapabilities(supports=supports)

    def schema(self) -> WorkflowSchema:
        return WorkflowSchema(
            config_schema={"type": "object"},
            sample_schema={"type": "array"},
            option_schema={"debug": {"type": "boolean"}},
        )

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

    def extract_artifacts(self, inputs, workspace):
        return Result.success(())


def test_service_rejects_non_workflow_registry_registry():
    with pytest.raises(ValueError, match="WorkflowRegistry"):
        WorkflowInfoService(registry=object())


def test_list_workflows_preserves_registry_order():
    adapters = [
        FakeAdapter(workflow_id="alpha"),
        FakeAdapter(workflow_id="beta"),
        FakeAdapter(workflow_id="gamma"),
    ]
    service = WorkflowInfoService(registry=WorkflowRegistry(adapters=adapters))

    workflows = service.list_workflows()

    assert [metadata.workflow_id for metadata in workflows] == [
        "alpha",
        "beta",
        "gamma",
    ]
    assert all(isinstance(metadata, WorkflowMetadata) for metadata in workflows)


def test_list_workflows_returns_empty_list_when_registry_empty():
    service = WorkflowInfoService(registry=WorkflowRegistry())

    assert service.list_workflows() == []


def test_get_schema_success_returns_adapter_schema():
    adapter = FakeAdapter(workflow_id="fake")
    service = WorkflowInfoService(registry=WorkflowRegistry(adapters=[adapter]))

    result = service.get_schema("fake")

    assert result.is_success
    assert result.value == adapter.schema()


def test_get_schema_unknown_workflow_returns_workflow_not_found():
    service = WorkflowInfoService(registry=WorkflowRegistry())

    result = service.get_schema("missing")

    assert result.is_failure
    assert result.value is None
    assert result.errors == (
        Issue(
            code="WORKFLOW_NOT_FOUND",
            message="Workflow was not found.",
            source="workflow_info",
            path="workflow_id",
            context={"workflow_id": "missing"},
        ),
    )


def test_get_capabilities_success_returns_adapter_capabilities():
    adapter = FakeAdapter(workflow_id="fake", supports=("validation", "dag_preview"))
    service = WorkflowInfoService(registry=WorkflowRegistry(adapters=[adapter]))

    result = service.get_capabilities("fake")

    assert result.is_success
    assert result.value == adapter.capabilities


def test_get_capabilities_unknown_workflow_returns_workflow_not_found():
    service = WorkflowInfoService(registry=WorkflowRegistry())

    result = service.get_capabilities("missing")

    assert result.is_failure
    assert result.value is None
    assert result.errors == (
        Issue(
            code="WORKFLOW_NOT_FOUND",
            message="Workflow was not found.",
            source="workflow_info",
            path="workflow_id",
            context={"workflow_id": "missing"},
        ),
    )


def test_invalid_workflow_id_from_registry_get_propagates_value_error():
    service = WorkflowInfoService(registry=WorkflowRegistry())

    with pytest.raises(ValueError):
        service.get_schema("")

    with pytest.raises(ValueError):
        service.get_capabilities("")
