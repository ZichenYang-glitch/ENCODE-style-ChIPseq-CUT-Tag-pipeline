"""Workflow validation MVP routes."""

from __future__ import annotations

from typing import Any

from fastapi import APIRouter, Depends, Request
from fastapi.responses import JSONResponse

from encode_pipeline.api.dependencies import get_registry, get_validation_service
from encode_pipeline.api.models import (
    IssueResponse,
    SchemaResponse,
    ValidationRequest,
    ValidationResponse,
    WorkflowCapabilityResponse,
    WorkflowListItem,
    WorkflowListResponse,
    WorkflowMetadataResponse,
)
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue
from encode_pipeline.services.validation import ValidationService


router = APIRouter(prefix="/workflows", tags=["workflows"])


VALIDATION_CAPABILITY = "validation"


def _workflow_not_found_issue(workflow_id: str) -> Issue:
    return Issue(
        code="WORKFLOW_NOT_FOUND",
        message="Workflow was not found.",
        severity="error",
        path="workflow_id",
        source="registry",
        context={"workflow_id": workflow_id},
    )


def _capability_unsupported_issue(workflow_id: str) -> Issue:
    return Issue(
        code="WORKFLOW_CAPABILITY_UNSUPPORTED",
        message="Workflow does not support validation.",
        severity="error",
        path="workflow.capabilities",
        source="registry",
        context={"workflow_id": workflow_id, "capability": VALIDATION_CAPABILITY},
    )


def _not_found_response(workflow_id: str, body: dict) -> JSONResponse:
    return JSONResponse(status_code=404, content=body)


def _conflict_response(workflow_id: str, body: dict) -> JSONResponse:
    return JSONResponse(status_code=409, content=body)


@router.get("/", response_model=WorkflowListResponse)
async def list_workflows(
    registry: WorkflowRegistry = Depends(get_registry),
) -> WorkflowListResponse:
    """List registered workflows with metadata and capabilities."""
    items: list[WorkflowListItem] = []
    for metadata in registry.list_metadata():
        adapter = registry.get(metadata.workflow_id)
        items.append(
            WorkflowListItem(
                metadata=WorkflowMetadataResponse(**metadata.to_dict()),
                capabilities=WorkflowCapabilityResponse(**adapter.capabilities.to_dict()),
            )
        )
    return WorkflowListResponse(ok=True, workflows=items, issues=[])


@router.get("/{workflow_id}/schema", response_model=SchemaResponse)
async def get_schema(
    workflow_id: str,
    registry: WorkflowRegistry = Depends(get_registry),
) -> SchemaResponse:
    """Return adapter-owned schema hints for one workflow."""
    try:
        adapter = registry.get(workflow_id)
    except KeyError:
        return _not_found_response(
            workflow_id,
            SchemaResponse(
                ok=False,
                workflow_id=workflow_id,
                schema_hints={},
                issues=[_workflow_not_found_issue(workflow_id).to_dict()],
            ).model_dump(),
        )

    return SchemaResponse(
        ok=True,
        workflow_id=workflow_id,
        schema_hints=dict(adapter.schema().to_dict()),
        issues=[],
    )


@router.post("/{workflow_id}/validate", response_model=ValidationResponse)
async def validate_workflow(
    workflow_id: str,
    request_body: ValidationRequest,
    validation_service: ValidationService = Depends(get_validation_service),
) -> ValidationResponse:
    """Validate submitted workflow inputs."""
    inputs = WorkflowInputs(
        config=request_body.config,
        samples=request_body.samples,
        options=request_body.options,
    )
    result = validation_service.validate(workflow_id, inputs)

    if result.is_failure and result.issues:
        first = result.issues[0]
        if first.code == "WORKFLOW_NOT_FOUND":
            return _not_found_response(
                workflow_id,
                ValidationResponse(
                    ok=False,
                    workflow_id=workflow_id,
                    value=None,
                    issues=[issue.to_dict() for issue in result.issues],
                ).model_dump(),
            )
        if first.code == "WORKFLOW_CAPABILITY_UNSUPPORTED":
            return _conflict_response(
                workflow_id,
                ValidationResponse(
                    ok=False,
                    workflow_id=workflow_id,
                    value=None,
                    issues=[issue.to_dict() for issue in result.issues],
                ).model_dump(),
            )

    return ValidationResponse(
        ok=result.is_success,
        workflow_id=workflow_id,
        value=result.value,
        issues=[issue.to_dict() for issue in result.issues],
    )
