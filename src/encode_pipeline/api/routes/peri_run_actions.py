"""Browser-only confirmation for two operations in a live Peri session."""

from __future__ import annotations

import json
from typing import Annotated, Literal

from fastapi import APIRouter, Depends
from fastapi.responses import JSONResponse
from pydantic import BaseModel, ConfigDict, Field

from encode_pipeline.api.dependencies import (
    AuthApiError,
    get_reference_profile_service,
    get_run_cancellation_service,
    get_run_service,
    get_validated_run_creation_service,
    require_authoring_actor,
    require_principal,
)
from encode_pipeline.api.models import RunCreateRequest, RunResponse
from encode_pipeline.api.routes.authoring import get_authoring_service
from encode_pipeline.api.routes.runs import cancel_run_response, create_run_response
from encode_pipeline.persistence.authoring import AuthoringSessionRecord
from encode_pipeline.platform.authentication import AuthenticatedPrincipal
from encode_pipeline.platform.authoring import (
    PERI_RUN_CONFIRMATION_INVALID,
    PERI_RUN_SNAPSHOT_CHANGED,
    AuthoringActor,
    AuthoringRejected,
    PeriRunAuthorization,
)
from encode_pipeline.services.authoring import AuthoringService
from encode_pipeline.services.peri_run_confirmation import (
    issue_confirmation,
    verified_confirmation_expiry,
)
from encode_pipeline.services.reference_profiles import ReferenceProfileService
from encode_pipeline.services.run_cancellation import RunCancellationService
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.validated_inputs import ValidatedRunCreationService


router = APIRouter(
    prefix="/agent/authoring-sessions/{sid}/run-actions", include_in_schema=False
)


class CreateTarget(RunCreateRequest):
    action: Literal["create"]

    def run_request(self) -> RunCreateRequest:
        return RunCreateRequest(snapshot_id=self.snapshot_id, tags=self.tags)


class CancelTarget(BaseModel):
    model_config = ConfigDict(extra="forbid")

    action: Literal["cancel"]
    run_id: str = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")


Target = Annotated[CreateTarget | CancelTarget, Field(discriminator="action")]


class PreviewRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    target: Target


class ConfirmRequest(PreviewRequest):
    confirmation: str = Field(pattern=r"^[0-9]{1,12}\.[0-9a-f]{64}$")


@router.post("/preview")
def preview_run_action(
    sid: str,
    body: PreviewRequest,
    actor: AuthoringActor = Depends(require_authoring_actor),
    authoring: AuthoringService = Depends(get_authoring_service),
    runs: RunService = Depends(get_run_service),
) -> JSONResponse:
    connection = authoring.require_run_action_connection(actor, sid)
    _require_target(connection, body.target, runs)
    target = body.target.model_dump()
    return JSONResponse(
        {
            "target": target,
            "workflow_id": connection.workflow_id,
            "confirmation": issue_confirmation(connection, target),
        },
        headers={"Cache-Control": "no-store"},
    )


@router.post("/confirm")
def confirm_run_action(
    sid: str,
    body: ConfirmRequest,
    actor: AuthoringActor = Depends(require_authoring_actor),
    principal: AuthenticatedPrincipal = Depends(require_principal),
    authoring: AuthoringService = Depends(get_authoring_service),
    runs: RunService = Depends(get_run_service),
    creation: ValidatedRunCreationService = Depends(get_validated_run_creation_service),
    cancellation: RunCancellationService = Depends(get_run_cancellation_service),
    references: ReferenceProfileService = Depends(get_reference_profile_service),
) -> JSONResponse:
    connection = authoring.require_run_action_connection(actor, sid)
    target = body.target
    expiry = verified_confirmation_expiry(
        connection, target.model_dump(), body.confirmation
    )
    if expiry is None:
        raise AuthoringRejected(PERI_RUN_CONFIRMATION_INVALID)
    _require_target(connection, target, runs)
    authorization = PeriRunAuthorization(
        actor=actor,
        sid=sid,
        generation=connection.generation,
        instance_id=connection.instance_id,
        owner_epoch=connection.bridge_owner_epoch,
        workflow_id=connection.workflow_id,
        action=target.action,
        target_id=target.snapshot_id
        if isinstance(target, CreateTarget)
        else target.run_id,
        expires_at=expiry,
    )
    if isinstance(target, CreateTarget):
        result = create_run_response(
            connection.workflow_id,
            target.run_request(),
            principal,
            creation,
            runs,
            references,
            peri_authorization=authorization,
        )
        return _summarize(result, default_status=201)
    result = cancel_run_response(
        target.run_id, principal, cancellation, peri_authorization=authorization
    )
    return _summarize(result, default_status=200)


def _require_target(
    connection: AuthoringSessionRecord, target: Target, runs: RunService
) -> None:
    if isinstance(target, CreateTarget):
        if (
            target.snapshot_id != connection.last_snapshot_id
            or connection.last_validated_revision != connection.revision
        ):
            raise AuthoringRejected(PERI_RUN_SNAPSHOT_CHANGED)
        return
    try:
        if runs.get_run(target.run_id).workflow_id == connection.workflow_id:
            return
    except KeyError:
        pass
    raise AuthApiError(404, "RUN_NOT_FOUND", "Run was not found in this workflow.")


def _summarize(
    result: RunResponse | JSONResponse, *, default_status: int
) -> JSONResponse:
    if isinstance(result, JSONResponse):
        status = result.status_code
        body = json.loads(result.body)
    else:
        status = default_status
        body = result.model_dump(mode="json")
    if status >= 400:
        return JSONResponse({"issues": body["issues"]}, status_code=status)
    return JSONResponse(
        {"run_id": body["run"]["run_id"], "stop_requested": status == 202},
        status_code=status,
    )
