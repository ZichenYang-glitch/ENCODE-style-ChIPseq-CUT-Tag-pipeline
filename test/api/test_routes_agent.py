"""Tests for the FastAPI workflow-scoped agent chat route."""

from __future__ import annotations

import ast
import importlib
import inspect
from collections.abc import Iterator
from pathlib import Path

import pytest

from encode_pipeline.api.main import create_app
from encode_pipeline.services.agent import AgentService
from encode_pipeline.services.defaults import create_default_agent_service
from encode_pipeline.services.llm_client import MockLLMClient
from encode_pipeline.services.validation import ValidationService
from encode_pipeline.services.workflow_info import WorkflowInfoService
from api_test_client import ApiTestClient

fastapi = pytest.importorskip("fastapi")

WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
DEFAULT_MOCK_RESPONSE = "This is a deterministic mock explanation from the workflow agent."


@pytest.fixture
def client() -> Iterator[ApiTestClient]:
    """Default app wired to the bundled ENCODE-style adapter and mock LLM."""
    app = create_app()
    with ApiTestClient(app) as tc:
        yield tc


@pytest.fixture
def deterministic_client() -> Iterator[ApiTestClient]:
    """App with a deterministic mock LLM client for reproducibility checks."""
    app = create_app()
    registry = app.state.registry
    workflow_info = WorkflowInfoService(registry=registry)
    validation_service = ValidationService(registry=registry)
    llm_client = MockLLMClient(response_text="Deterministic mock explanation.")
    app.state.agent_service = AgentService(workflow_info, validation_service, llm_client)
    with ApiTestClient(app) as tc:
        yield tc


def test_agent_chat_returns_agent_response_for_known_workflow(client: ApiTestClient) -> None:
    response = client.post(
        f"/api/v1/workflows/{WORKFLOW_ID}/agent/chat",
        json={"session_id": "sess-1", "message": "Explain the issues."},
    )
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert data["session_id"] == "sess-1"
    assert data["message"] == DEFAULT_MOCK_RESPONSE
    assert data["suggestions"] == []
    assert data["issues"] == []
    assert "tool_calls" in data


def test_agent_chat_message_is_deterministic(deterministic_client: ApiTestClient) -> None:
    payload = {"session_id": "sess-2", "message": "Explain the issues."}
    first = deterministic_client.post(
        f"/api/v1/workflows/{WORKFLOW_ID}/agent/chat",
        json=payload,
    )
    second = deterministic_client.post(
        f"/api/v1/workflows/{WORKFLOW_ID}/agent/chat",
        json=payload,
    )
    assert first.status_code == 200
    assert second.status_code == 200
    first_data = first.json()
    second_data = second.json()
    assert first_data["message"] == "Deterministic mock explanation."
    assert first_data["message"] == second_data["message"]
    assert first_data["tool_calls"] == second_data["tool_calls"]


def test_agent_chat_includes_read_only_explain_issues_tool_call(client: ApiTestClient) -> None:
    response = client.post(
        f"/api/v1/workflows/{WORKFLOW_ID}/agent/chat",
        json={"message": "Explain the issues."},
    )
    assert response.status_code == 200
    data = response.json()
    assert len(data["tool_calls"]) == 1
    tool_call = data["tool_calls"][0]
    assert tool_call["tool_name"] == "explain_issues"
    assert tool_call["read_only"] is True
    assert "input_summary" in tool_call
    assert "output_summary" in tool_call


def test_agent_chat_unknown_workflow_returns_ok_false(client: ApiTestClient) -> None:
    response = client.post(
        "/api/v1/workflows/missing-workflow/agent/chat",
        json={"message": "Explain the issues."},
    )
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is False
    assert data["session_id"] is None
    assert "missing-workflow" in data["message"]
    assert len(data["issues"]) == 1
    issue = data["issues"][0]
    assert issue["code"] == "WORKFLOW_NOT_FOUND"
    assert issue["context"]["workflow_id"] == "missing-workflow"
    assert data["tool_calls"] == []


def test_agent_chat_does_not_require_api_keys_or_network(client: ApiTestClient) -> None:
    """The default composition uses MockLLMClient, so no keys/network are needed."""
    response = client.post(
        f"/api/v1/workflows/{WORKFLOW_ID}/agent/chat",
        json={"message": "Explain the issues."},
    )
    assert response.status_code == 200
    data = response.json()
    assert "api key" not in data["message"].lower()
    assert "apikey" not in data["message"].lower()
    assert data["message"] == DEFAULT_MOCK_RESPONSE


def test_agent_chat_malformed_request_returns_400(client: ApiTestClient) -> None:
    response = client.post(
        f"/api/v1/workflows/{WORKFLOW_ID}/agent/chat",
        content="not-json",
        headers={"content-type": "application/json"},
    )
    assert response.status_code == 400
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "API_REQUEST_INVALID"


def test_agent_route_uses_app_state_services(client: ApiTestClient) -> None:
    """The endpoint should resolve the agent service from app.state, not recreate it."""
    app = client.app
    assert hasattr(app.state, "agent_service")
    assert isinstance(app.state.agent_service, AgentService)
    assert hasattr(app.state, "registry")
    assert hasattr(app.state, "validation_service")


def test_agent_route_registry_matches_validation_service_registry(client: ApiTestClient) -> None:
    """create_app should compose services from a single shared registry instance."""
    app = client.app
    agent_service = app.state.agent_service
    validation_service = app.state.validation_service
    registry = app.state.registry

    assert agent_service._workflow_info._registry is registry
    assert validation_service._registry is registry


def test_default_agent_service_uses_mock_llm_client() -> None:
    service = create_default_agent_service()
    assert isinstance(service, AgentService)
    assert isinstance(service._llm_client, MockLLMClient)


_FORBIDDEN_IMPORTS = {
    "encode_pipeline.adapters.encode",
    "encode_pipeline.config.validator",
    "encode_pipeline.samples",
    "snakemake",
    "subprocess",
    "openai",
    "anthropic",
    "google",
    "boto3",
    "azure",
    "langchain",
    "langgraph",
    "mastra",
    "transformers",
    "torch",
}

_FORBIDDEN_TOKENS = {
    "encode_pipeline.adapters.encode",
    "encode_pipeline.config.validator",
    "encode_pipeline.samples",
    "snakemake",
    "subprocess",
}


def _agent_route_source_path() -> Path:
    module = importlib.import_module("encode_pipeline.api.routes.agent")
    return Path(inspect.getfile(module))


def _collect_imported_dotted_names(source: str) -> set[str]:
    tree = ast.parse(source)
    imported: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                imported.add(alias.name)
        elif isinstance(node, ast.ImportFrom):
            if node.module is not None:
                prefix = node.module
                if node.level:
                    prefix = f"{'.' * node.level}{prefix}"
                imported.add(prefix)
    return imported


def test_agent_route_does_not_import_forbidden_modules() -> None:
    source = _agent_route_source_path().read_text(encoding="utf-8")
    imported = _collect_imported_dotted_names(source)
    violations = sorted(
        forbidden
        for forbidden in _FORBIDDEN_IMPORTS
        if any(name == forbidden or name.startswith(forbidden + ".") for name in imported)
    )
    assert not violations, f"agent route imports forbidden modules: {violations}"


def test_agent_route_does_not_reference_forbidden_tokens() -> None:
    source = _agent_route_source_path().read_text(encoding="utf-8")
    violations = sorted(token for token in _FORBIDDEN_TOKENS if token in source)
    assert not violations, f"agent route references forbidden tokens: {violations}"
