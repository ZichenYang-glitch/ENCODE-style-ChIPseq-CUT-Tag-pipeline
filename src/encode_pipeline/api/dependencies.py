"""FastAPI dependencies shared by API routes."""

from __future__ import annotations

from fastapi import Request

from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.validation import ValidationService


def get_registry(request: Request) -> WorkflowRegistry:
    """Return the app registry."""
    return request.app.state.registry


def get_validation_service(request: Request) -> ValidationService:
    """Return the app validation service."""
    return request.app.state.validation_service
