"""FastAPI dependencies shared by API routes."""

from __future__ import annotations

from typing import TYPE_CHECKING

from fastapi import Request

from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.validation import ValidationService

if TYPE_CHECKING:
    from encode_pipeline.services.agent import AgentService
    from encode_pipeline.services.runs import RunService


def get_registry(request: Request) -> WorkflowRegistry:
    """Return the app registry."""
    return request.app.state.registry


def get_validation_service(request: Request) -> ValidationService:
    """Return the app validation service."""
    return request.app.state.validation_service


def get_agent_service(request: Request) -> "AgentService":
    """Return the app agent service."""
    return request.app.state.agent_service


def get_run_service(request: Request) -> "RunService":
    """Return the app run service."""
    return request.app.state.run_service
