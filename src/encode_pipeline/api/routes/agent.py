"""Workflow-scoped agent chat route."""

from __future__ import annotations

from fastapi import APIRouter, Depends

from encode_pipeline.api.dependencies import get_agent_service
from encode_pipeline.api.models import AgentRequest, AgentResponse
from encode_pipeline.services.agent import AgentService


router = APIRouter(prefix="/workflows", tags=["agent"])


@router.post("/{workflow_id}/agent/chat", response_model=AgentResponse)
async def chat_with_workflow_agent(
    workflow_id: str,
    request_body: AgentRequest,
    agent_service: AgentService = Depends(get_agent_service),
) -> AgentResponse:
    """Chat with the workflow-scoped agent for a registered workflow.

    The endpoint is read-only and backed by a deterministic mock LLM client.
    It returns a structured envelope; an unknown workflow yields ``ok=false``
    with a ``WORKFLOW_NOT_FOUND`` issue rather than an HTTP 404.
    """
    return await agent_service.chat(workflow_id, request_body)
