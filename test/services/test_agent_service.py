"""Tests for the agent service."""

from __future__ import annotations

import asyncio

import pytest

from encode_pipeline.api.models import (
    AgentContext,
    AgentRequest,
    AgentResponse,
    AgentToolCall,
    IssueResponse,
)
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
from encode_pipeline.services.agent import AgentService
from encode_pipeline.services.agent_tools import (
    AgentTool,
    ReadOnlyToolRegistry,
    build_explain_issues_tool,
)
from encode_pipeline.services.llm_client import MockLLMClient
from encode_pipeline.services.validation import ValidationService
from encode_pipeline.services.workflow_info import WorkflowInfoService


class FakeAdapter:
    """Minimal workflow adapter for agent service tests."""

    def __init__(self, workflow_id: str = "fake") -> None:
        self.metadata = WorkflowMetadata(
            workflow_id=workflow_id,
            name="Fake Workflow",
            version="1.0.0",
        )
        self.capabilities = WorkflowCapabilities(supports=("validation",))

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
        return Result.success(WorkspacePlan(directories=[str(workspace)]))

    def build_command(self, plan: WorkspacePlan) -> Result[CommandSpec]:
        return Result.success(CommandSpec(argv=["run-workflow"]))


_DENIED_PREFIXES = (
    "run_",
    "submit_",
    "execute_",
    "delete_",
    "update_",
    "apply_",
    "write_",
    "modify_",
)


def _run_chat(
    service: AgentService,
    workflow_id: str,
    request: AgentRequest,
) -> AgentResponse:
    return asyncio.run(service.chat(workflow_id, request))


def test_chat_returns_agent_response_for_valid_workflow():
    adapter = FakeAdapter(workflow_id="fake")
    registry = WorkflowRegistry(adapters=[adapter])
    workflow_info = WorkflowInfoService(registry=registry)
    validation_service = ValidationService(registry=registry)
    llm_client = MockLLMClient(response_text="This is the LLM explanation.")
    service = AgentService(workflow_info, validation_service, llm_client)
    request = AgentRequest(session_id="sess-1", message="Explain the issues.")

    response = _run_chat(service, "fake", request)

    assert isinstance(response, AgentResponse)
    assert response.ok is True
    assert response.session_id == "sess-1"
    assert response.message == "This is the LLM explanation."
    assert response.suggestions == []
    assert response.issues == []
    assert len(response.tool_calls) == 1
    tool_call = response.tool_calls[0]
    assert isinstance(tool_call, AgentToolCall)
    assert tool_call.tool_name == "explain_issues"
    assert tool_call.read_only is True
    assert tool_call.input_summary == {"issues": []}
    assert tool_call.output_summary == "This is the LLM explanation."


def test_chat_unknown_workflow_returns_workflow_not_found():
    workflow_info = WorkflowInfoService(registry=WorkflowRegistry())
    validation_service = ValidationService(registry=WorkflowRegistry())
    llm_client = MockLLMClient(response_text="should not be used")
    service = AgentService(workflow_info, validation_service, llm_client)
    request = AgentRequest(session_id="sess-2", message="Explain the issues.")

    response = _run_chat(service, "missing", request)

    assert isinstance(response, AgentResponse)
    assert response.ok is False
    assert response.session_id == "sess-2"
    assert "not found" in response.message.lower()
    assert response.suggestions == []
    assert response.tool_calls == []
    assert len(response.issues) == 1
    issue = response.issues[0]
    assert issue.code == "WORKFLOW_NOT_FOUND"
    assert "missing" in (issue.context.get("workflow_id", "") or "")


def test_chat_tool_calls_are_read_only_and_not_mutating():
    adapter = FakeAdapter(workflow_id="fake")
    registry = WorkflowRegistry(adapters=[adapter])
    workflow_info = WorkflowInfoService(registry=registry)
    validation_service = ValidationService(registry=registry)
    llm_client = MockLLMClient(response_text="ok")
    service = AgentService(workflow_info, validation_service, llm_client)
    request = AgentRequest(message="Check workflow.")

    response = _run_chat(service, "fake", request)

    assert all(call.read_only is True for call in response.tool_calls)
    assert not any(
        any(call.tool_name.startswith(prefix) for prefix in _DENIED_PREFIXES)
        for call in response.tool_calls
    )


def test_chat_passes_context_issues_to_explain_tool():
    adapter = FakeAdapter(workflow_id="fake")
    registry = WorkflowRegistry(adapters=[adapter])
    workflow_info = WorkflowInfoService(registry=registry)
    validation_service = ValidationService(registry=registry)
    llm_client = MockLLMClient(response_text="ok")
    service = AgentService(workflow_info, validation_service, llm_client)
    issue = IssueResponse(code="FRIP_LOW", message="FRiP is low.")
    request = AgentRequest(
        session_id="sess-3",
        message="Why is FRiP low?",
        context=AgentContext(current_issues=[issue]),
    )

    response = _run_chat(service, "fake", request)

    assert response.ok is True
    assert len(response.tool_calls) == 1
    input_summary = response.tool_calls[0].input_summary
    assert "issues" in input_summary
    assert len(input_summary["issues"]) == 1
    assert input_summary["issues"][0]["code"] == "FRIP_LOW"


def test_chat_uses_mock_llm_client_deterministically():
    adapter = FakeAdapter(workflow_id="fake")
    registry = WorkflowRegistry(adapters=[adapter])
    workflow_info = WorkflowInfoService(registry=registry)
    validation_service = ValidationService(registry=registry)
    llm_client = MockLLMClient(response_text="deterministic")
    service = AgentService(workflow_info, validation_service, llm_client)
    request = AgentRequest(message="Hello.")

    first = _run_chat(service, "fake", request)
    second = _run_chat(service, "fake", request)

    assert first.message == "deterministic"
    assert second.message == "deterministic"
    assert first.tool_calls[0].output_summary == "deterministic"
    assert second.tool_calls[0].output_summary == "deterministic"
    assert first.model_dump() == second.model_dump()


def test_chat_accepts_optional_tool_registry():
    adapter = FakeAdapter(workflow_id="fake")
    registry = WorkflowRegistry(adapters=[adapter])
    workflow_info = WorkflowInfoService(registry=registry)
    validation_service = ValidationService(registry=registry)
    llm_client = MockLLMClient(response_text="custom registry")
    custom_registry = ReadOnlyToolRegistry()
    custom_registry.register(build_explain_issues_tool(llm_client))
    service = AgentService(
        workflow_info,
        validation_service,
        llm_client,
        tool_registry=custom_registry,
    )
    request = AgentRequest(message="Hello.")

    response = _run_chat(service, "fake", request)

    assert response.ok is True
    assert response.message == "custom registry"
    assert len(response.tool_calls) == 1
    assert response.tool_calls[0].tool_name == "explain_issues"
    assert response.tool_calls[0].output_summary == "custom registry"


def test_chat_does_not_mutate_request_context():
    adapter = FakeAdapter(workflow_id="fake")
    registry = WorkflowRegistry(adapters=[adapter])
    workflow_info = WorkflowInfoService(registry=registry)
    validation_service = ValidationService(registry=registry)
    llm_client = MockLLMClient(response_text="ok")
    service = AgentService(workflow_info, validation_service, llm_client)
    issue = IssueResponse(code="X", message="m")
    request = AgentRequest(
        message="Explain",
        context=AgentContext(current_issues=[issue]),
    )
    original_issues = list(request.context.current_issues)

    _run_chat(service, "fake", request)

    assert request.context.current_issues == original_issues
    assert request.message == "Explain"


def test_chat_response_message_matches_explain_issues_tool_output():
    adapter = FakeAdapter(workflow_id="fake")
    registry = WorkflowRegistry(adapters=[adapter])
    workflow_info = WorkflowInfoService(registry=registry)
    validation_service = ValidationService(registry=registry)
    llm_client = MockLLMClient(response_text="matched explanation")
    service = AgentService(workflow_info, validation_service, llm_client)
    request = AgentRequest(message="Explain.")

    response = _run_chat(service, "fake", request)

    assert response.message == "matched explanation"
    assert response.tool_calls[0].output_summary == response.message
