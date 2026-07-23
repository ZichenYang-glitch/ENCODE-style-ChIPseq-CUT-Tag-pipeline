"""Workflow validation MVP routes."""

from __future__ import annotations

from fastapi import APIRouter, Depends
from fastapi.responses import JSONResponse

from encode_pipeline.api.dependencies import get_registry, get_validated_input_service
from encode_pipeline.api.models import (
    SchemaResponse,
    WorkflowDetailResponse,
    WorkflowSchemaResponse,
    ValidationRequest,
    ValidationResponse,
    ValidatedInputSnapshotResponse,
    WorkflowListItem,
    WorkflowListResponse,
)
from encode_pipeline.platform.adapters import VALIDATION_CAPABILITY, WorkflowInputs
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue
from encode_pipeline.services.validated_inputs import ValidatedInputService
from encode_pipeline.services.workflow_info import WorkflowInfoService


router = APIRouter(prefix="/workflows", tags=["workflows"])


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


def _descriptor_item(descriptor) -> WorkflowListItem:
    return WorkflowListItem.model_validate(descriptor.to_dict())


@router.get("/", response_model=WorkflowListResponse, operation_id="listWorkflows")
def list_workflows(
    registry: WorkflowRegistry = Depends(get_registry),
) -> WorkflowListResponse:
    """List registered workflows with metadata and capabilities."""
    items = [
        _descriptor_item(descriptor)
        for descriptor in WorkflowInfoService(registry).list_descriptors()
    ]
    return WorkflowListResponse(ok=True, workflows=items, issues=[])


@router.get(
    "/{workflow_id}",
    response_model=WorkflowDetailResponse,
    operation_id="getWorkflow",
    responses={
        404: {"model": WorkflowDetailResponse},
        503: {"model": WorkflowDetailResponse},
    },
)
def get_workflow(
    workflow_id: str,
    registry: WorkflowRegistry = Depends(get_registry),
) -> WorkflowDetailResponse | JSONResponse:
    """Return one safe workflow product descriptor."""
    result = WorkflowInfoService(registry).get_descriptor(workflow_id)
    if result.is_failure or result.value is None:
        issue = result.issues[0]
        status_code = 404 if issue.code == "WORKFLOW_NOT_FOUND" else 503
        return JSONResponse(
            status_code=status_code,
            content=WorkflowDetailResponse(
                ok=False,
                workflow_id=workflow_id,
                workflow=None,
                issues=[item.to_dict() for item in result.issues],
            ).model_dump(mode="json"),
        )
    return WorkflowDetailResponse(
        ok=True,
        workflow_id=workflow_id,
        workflow=_descriptor_item(result.value),
        issues=[],
    )


@router.get(
    "/{workflow_id}/schema",
    response_model=SchemaResponse,
    operation_id="getWorkflowSchema",
)
async def get_schema(
    workflow_id: str,
    registry: WorkflowRegistry = Depends(get_registry),
) -> SchemaResponse:
    """Return the adapter-owned authoring contract for one workflow."""
    try:
        adapter = registry.get(workflow_id)
    except KeyError:
        return _not_found_response(
            workflow_id,
            SchemaResponse(
                ok=False,
                workflow_id=workflow_id,
                schema=None,
                issues=[_workflow_not_found_issue(workflow_id).to_dict()],
            ).model_dump(by_alias=True),
        )

    return SchemaResponse(
        ok=True,
        workflow_id=workflow_id,
        schema=WorkflowSchemaResponse(**adapter.schema().to_dict()),
        issues=[],
    )


@router.post(
    "/{workflow_id}/validate",
    response_model=ValidationResponse,
    operation_id="validateWorkflow",
    responses={
        413: {
            "model": ValidationResponse,
            "description": "Request body too large.",
        },
        404: {"model": ValidationResponse},
        409: {"model": ValidationResponse},
        500: {"model": ValidationResponse},
        503: {"model": ValidationResponse},
    },
)
def validate_workflow(
    workflow_id: str,
    request_body: ValidationRequest,
    validation_service: ValidatedInputService = Depends(get_validated_input_service),
) -> ValidationResponse | JSONResponse:
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
                    snapshot=None,
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
                    snapshot=None,
                    issues=[issue.to_dict() for issue in result.issues],
                ).model_dump(),
            )
        status_by_code = {
            "VALIDATION_WORKFLOW_BUILD_CHANGED": 409,
            "VALIDATION_WORKFLOW_BUILD_UNAVAILABLE": 503,
            "VALIDATION_WORKFLOW_SCHEMA_UNAVAILABLE": 503,
            "VALIDATED_SNAPSHOT_PERSISTENCE_FAILED": 500,
        }
        status_code = status_by_code.get(first.code)
        if status_code is not None:
            return JSONResponse(
                status_code=status_code,
                content=ValidationResponse(
                    ok=False,
                    workflow_id=workflow_id,
                    value=None,
                    snapshot=None,
                    issues=[issue.to_dict() for issue in result.issues],
                ).model_dump(mode="json"),
            )

    return ValidationResponse(
        ok=result.is_success,
        workflow_id=workflow_id,
        value=None,
        snapshot=(
            ValidatedInputSnapshotResponse.from_snapshot(result.value)
            if result.value is not None
            else None
        ),
        issues=[issue.to_dict() for issue in result.issues],
    )
