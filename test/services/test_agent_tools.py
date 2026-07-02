"""Tests for the read-only agent tool registry and tool builders."""

from __future__ import annotations

import asyncio

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
from encode_pipeline.services.agent_tools import (
    AgentTool,
    ReadOnlyToolRegistry,
    build_explain_issues_tool,
    build_get_workflow_schema_tool,
    build_list_workflows_tool,
    build_validate_workflow_input_tool,
)
from encode_pipeline.services.llm_client import MockLLMClient
from encode_pipeline.services.validation import ValidationService
from encode_pipeline.services.workflow_info import WorkflowInfoService


class FakeAdapter:
    """Minimal workflow adapter for service-layer tests."""

    def __init__(self, workflow_id: str = "fake") -> None:
        self.metadata = WorkflowMetadata(
            workflow_id=workflow_id,
            name="Fake Workflow",
            version="1.0.0",
        )
        self.capabilities = WorkflowCapabilities(supports=("validation",))
        self.last_inputs: WorkflowInputs | None = None

    def schema(self) -> WorkflowSchema:
        return WorkflowSchema(config_schema={"type": "object"})

    def validate(self, inputs: WorkflowInputs) -> Result[object]:
        self.last_inputs = inputs
        return Result.success({"validated": True})

    def preview_dag(self, inputs: WorkflowInputs) -> Result[DagPreview]:
        return Result.success(DagPreview())

    def plan_workspace(
        self,
        inputs: WorkflowInputs,
        workspace: str,
    ) -> Result[WorkspacePlan]:
        return Result.success(WorkspacePlan(directories=[str(workspace)]))

    def build_command(self, plan: WorkspacePlan) -> Result[CommandSpec]:
        return Result.success(CommandSpec(argv=["run-workflow"]))


def _dummy_handler(**kwargs) -> dict:
    return kwargs


def test_registry_accepts_read_only_tool():
    registry = ReadOnlyToolRegistry()
    tool = AgentTool(
        name="describe_workflows",
        description="Lists workflows.",
        parameters_schema={"type": "object"},
        handler=_dummy_handler,
    )

    registry.register(tool)

    assert registry.get("describe_workflows") is tool


def test_registry_rejects_read_only_false():
    registry = ReadOnlyToolRegistry()
    tool = AgentTool(
        name="safe_name",
        description="A tool.",
        parameters_schema={"type": "object"},
        handler=_dummy_handler,
        read_only=False,
    )

    with pytest.raises(ValueError, match="read_only=True"):
        registry.register(tool)


def test_registry_rejects_duplicate_names():
    registry = ReadOnlyToolRegistry()
    tool = AgentTool(
        name="unique_tool",
        description="First.",
        parameters_schema={"type": "object"},
        handler=_dummy_handler,
    )
    duplicate = AgentTool(
        name="unique_tool",
        description="Second.",
        parameters_schema={"type": "object"},
        handler=_dummy_handler,
    )
    registry.register(tool)

    with pytest.raises(ValueError, match="already registered"):
        registry.register(duplicate)


@pytest.mark.parametrize(
    "name",
    [
        "run_workflow",
        "submit_job",
        "execute_command",
        "delete_workflow",
        "update_config",
        "apply_patch",
        "write_file",
        "modify_state",
    ],
)
def test_registry_rejects_denied_prefixes_even_when_read_only(name: str):
    registry = ReadOnlyToolRegistry()
    tool = AgentTool(
        name=name,
        description="A read-only tool with a denied prefix.",
        parameters_schema={"type": "object"},
        handler=_dummy_handler,
    )

    with pytest.raises(ValueError, match="denied prefix"):
        registry.register(tool)


def test_list_tools_preserves_registration_order():
    registry = ReadOnlyToolRegistry()
    first = AgentTool(
        name="alpha_tool",
        description="First.",
        parameters_schema={"type": "object"},
        handler=_dummy_handler,
    )
    second = AgentTool(
        name="beta_tool",
        description="Second.",
        parameters_schema={"type": "object"},
        handler=_dummy_handler,
    )
    third = AgentTool(
        name="gamma_tool",
        description="Third.",
        parameters_schema={"type": "object"},
        handler=_dummy_handler,
    )
    registry.register(first)
    registry.register(second)
    registry.register(third)

    tools = registry.list_tools()

    assert tools == [first, second, third]


def test_get_unknown_tool_raises_key_error():
    registry = ReadOnlyToolRegistry()

    with pytest.raises(KeyError, match="unknown_tool"):
        registry.get("unknown_tool")


def test_build_list_workflows_tool_returns_agent_tool():
    service = WorkflowInfoService(registry=WorkflowRegistry(adapters=[FakeAdapter()]))

    tool = build_list_workflows_tool(service)

    assert isinstance(tool, AgentTool)
    assert tool.name == "list_workflows"
    assert tool.read_only is True
    assert callable(tool.handler)


def test_list_workflows_handler_returns_result_with_workflow_metadata():
    adapter = FakeAdapter(workflow_id="fake")
    service = WorkflowInfoService(registry=WorkflowRegistry(adapters=[adapter]))
    tool = build_list_workflows_tool(service)

    result = tool.handler()

    assert isinstance(result, Result)
    assert result.is_success
    assert result.value == [adapter.metadata]


def test_build_get_workflow_schema_tool_returns_agent_tool():
    service = WorkflowInfoService(registry=WorkflowRegistry(adapters=[FakeAdapter()]))

    tool = build_get_workflow_schema_tool(service)

    assert isinstance(tool, AgentTool)
    assert tool.name == "get_workflow_schema"
    assert tool.read_only is True
    assert callable(tool.handler)


def test_get_workflow_schema_handler_returns_result_from_service():
    adapter = FakeAdapter(workflow_id="fake")
    service = WorkflowInfoService(registry=WorkflowRegistry(adapters=[adapter]))
    tool = build_get_workflow_schema_tool(service)

    result = tool.handler(workflow_id="fake")

    assert isinstance(result, Result)
    assert result.is_success
    assert result.value == adapter.schema()


def test_get_workflow_schema_handler_returns_failure_for_unknown_workflow():
    service = WorkflowInfoService(registry=WorkflowRegistry())
    tool = build_get_workflow_schema_tool(service)

    result = tool.handler(workflow_id="missing")

    assert isinstance(result, Result)
    assert result.is_failure
    assert result.value is None
    assert any(issue.code == "WORKFLOW_NOT_FOUND" for issue in result.issues)


def test_build_validate_workflow_input_tool_returns_agent_tool():
    service = ValidationService(registry=WorkflowRegistry(adapters=[FakeAdapter()]))

    tool = build_validate_workflow_input_tool(service)

    assert isinstance(tool, AgentTool)
    assert tool.name == "validate_workflow_input"
    assert tool.read_only is True
    assert callable(tool.handler)


def test_validate_workflow_input_handler_returns_result_from_service():
    adapter = FakeAdapter(workflow_id="fake")
    service = ValidationService(registry=WorkflowRegistry(adapters=[adapter]))
    tool = build_validate_workflow_input_tool(service)

    result = tool.handler(workflow_id="fake", config={"mode": "test"})

    assert isinstance(result, Result)
    assert result.is_success
    assert result.value == {"validated": True}


def test_validate_workflow_input_handler_passes_config_samples_and_options():
    adapter = FakeAdapter(workflow_id="fake")
    service = ValidationService(registry=WorkflowRegistry(adapters=[adapter]))
    tool = build_validate_workflow_input_tool(service)

    tool.handler(
        workflow_id="fake",
        config={"mode": "test"},
        samples="/path/to/samples.tsv",
        options={"debug": True},
    )

    assert adapter.last_inputs is not None
    assert adapter.last_inputs.config == {"mode": "test"}
    assert adapter.last_inputs.samples == "/path/to/samples.tsv"
    assert adapter.last_inputs.options == {"debug": True}


def test_validate_workflow_input_handler_defaults_samples_and_options():
    adapter = FakeAdapter(workflow_id="fake")
    service = ValidationService(registry=WorkflowRegistry(adapters=[adapter]))
    tool = build_validate_workflow_input_tool(service)

    tool.handler(workflow_id="fake", config={"mode": "test"})

    assert adapter.last_inputs is not None
    assert adapter.last_inputs.config == {"mode": "test"}
    assert adapter.last_inputs.samples is None
    assert adapter.last_inputs.options == {}


def test_build_explain_issues_tool_returns_agent_tool():
    client = MockLLMClient()

    tool = build_explain_issues_tool(client)

    assert isinstance(tool, AgentTool)
    assert tool.name == "explain_issues"
    assert tool.read_only is True
    assert callable(tool.handler)


def _run_handler(tool: AgentTool, **kwargs):
    coro = tool.handler(**kwargs)
    return asyncio.run(coro)


def test_explain_issues_handler_returns_result_from_llm_client():
    client = MockLLMClient(response_text="This is the explanation.")
    tool = build_explain_issues_tool(client)

    result = _run_handler(tool, issues=[{"code": "X"}])

    assert isinstance(result, Result)
    assert result.is_success
    assert result.value == "This is the explanation."


def test_explain_issues_handler_uses_mock_llm_client_deterministically():
    client = MockLLMClient(response_text="deterministic")
    tool = build_explain_issues_tool(client)

    first = _run_handler(tool, issues=[])
    second = _run_handler(tool, issues=[])

    assert first.value == "deterministic"
    assert second.value == "deterministic"
    assert first.value == second.value
