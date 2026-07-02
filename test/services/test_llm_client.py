"""Tests for the provider-neutral LLM client and agent API models."""

from __future__ import annotations

import asyncio

import pytest
from pydantic import ValidationError

from encode_pipeline.api.models import (
    AgentContext,
    AgentRequest,
    AgentResponse,
    AgentSuggestion,
    AgentToolCall,
    IssueResponse,
)
from encode_pipeline.services.llm_client import (
    LLMClient,
    LLMMessage,
    LLMResponse,
    MockLLMClient,
    ToolCall,
    ToolDefinition,
)


def _run(client: LLMClient, messages: list[LLMMessage], tools: list[ToolDefinition]) -> LLMResponse:
    return asyncio.run(client.complete(messages, tools))


def test_mock_llm_client_returns_configured_text_deterministically():
    client = MockLLMClient(response_text="hello")
    response = _run(client, [], [])
    assert response.content == "hello"
    assert response.tool_calls == ()


def test_mock_llm_client_returns_configured_tool_calls():
    calls = [ToolCall(name="explain_issues", arguments={"issue_codes": ["FRIP_LOW"]})]
    client = MockLLMClient(response_text="ok", tool_calls=calls)
    response = _run(client, [], [])
    assert response.content == "ok"
    assert len(response.tool_calls) == 1
    assert response.tool_calls[0].name == "explain_issues"
    assert response.tool_calls[0].arguments == {"issue_codes": ["FRIP_LOW"]}


def test_llm_response_defaults_to_empty_tool_calls():
    response = LLMResponse(content="just text")
    assert response.tool_calls == ()


def test_agent_suggestion_disclaimer_is_non_empty_and_mentions_model_generated():
    suggestion = AgentSuggestion(description="Try single-end mode.")
    assert suggestion.disclaimer
    assert "model-generated suggestion" in suggestion.disclaimer


def test_agent_tool_call_read_only_defaults_to_true():
    call = AgentToolCall(tool_name="list_workflows")
    assert call.read_only is True


def test_agent_tool_call_rejects_read_only_false():
    with pytest.raises(ValidationError):
        AgentToolCall(tool_name="list_workflows", read_only=False)


def test_agent_context_defaults_are_fresh_mutable_containers():
    ctx1 = AgentContext()
    ctx2 = AgentContext()
    issue = IssueResponse(code="X", message="m")
    ctx1.current_issues.append(issue)
    ctx1.current_config["key"] = "value"
    ctx1.current_schema["field"] = "value"
    assert ctx2.current_issues == []
    assert ctx2.current_config == {}
    assert ctx2.current_schema == {}
    assert ctx1.current_issues == [issue]
    assert ctx1.current_config == {"key": "value"}
    assert ctx1.current_schema == {"field": "value"}


def test_agent_request_round_trips_through_dump_and_validate():
    issue = IssueResponse(code="FRIP_LOW", message="FRiP is low.")
    request = AgentRequest(
        session_id="sess-1",
        message="Why is FRiP low?",
        context=AgentContext(
            current_issues=[issue],
            current_config={"assay": "chipseq"},
            current_schema={"type": "object"},
        ),
    )
    restored = AgentRequest.model_validate(request.model_dump())
    assert restored.session_id == "sess-1"
    assert restored.message == "Why is FRiP low?"
    assert restored.context is not None
    assert len(restored.context.current_issues) == 1
    assert restored.context.current_issues[0].code == "FRIP_LOW"
    assert restored.context.current_config == {"assay": "chipseq"}
    assert restored.context.current_schema == {"type": "object"}


def test_agent_response_round_trips_through_dump_and_validate():
    issue = IssueResponse(code="X", message="m")
    response = AgentResponse(
        ok=True,
        session_id="sess-2",
        message="Here is the answer.",
        suggestions=[
            AgentSuggestion(
                type="config_edit",
                description="Increase replicate count.",
                target_path="replicates",
            ),
        ],
        tool_calls=[
            AgentToolCall(
                tool_name="explain_issues",
                input_summary={"issue_codes": ["X"]},
                output_summary="explained",
            ),
        ],
        issues=[issue],
    )
    restored = AgentResponse.model_validate(response.model_dump())
    assert restored.ok is True
    assert restored.session_id == "sess-2"
    assert restored.message == "Here is the answer."
    assert len(restored.suggestions) == 1
    assert restored.suggestions[0].type == "config_edit"
    assert len(restored.tool_calls) == 1
    assert restored.tool_calls[0].tool_name == "explain_issues"
    assert len(restored.issues) == 1
    assert restored.issues[0].code == "X"


def test_agent_response_issues_can_hold_issue_response_instances():
    issue = IssueResponse(code="Y", message="n", severity="warning")
    response = AgentResponse(ok=True, message="m", issues=[issue])
    assert response.issues == [issue]
    assert isinstance(response.issues[0], IssueResponse)


def test_empty_message_fails_validation():
    with pytest.raises(ValidationError):
        AgentRequest(message="")


def test_missing_description_fails_validation():
    with pytest.raises(ValidationError):
        AgentSuggestion()


def test_empty_disclaimer_fails_validation():
    with pytest.raises(ValidationError):
        AgentSuggestion(description="Try single-end mode.", disclaimer="")
