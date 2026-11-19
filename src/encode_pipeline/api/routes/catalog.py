"""Slice A experiment catalog routes."""

from __future__ import annotations

from datetime import datetime
from typing import Annotated, Any

from fastapi import APIRouter, Depends, Header, Query, Request
from fastapi.responses import JSONResponse
from pydantic import WithJsonSchema
from starlette.concurrency import run_in_threadpool

from encode_pipeline.api.catalog_models import (
    CatalogAdoptProposalRequest,
    CatalogAdoptionReceiptListResponse,
    CatalogAdoptionReceiptResponse,
    CatalogArchiveExperimentResponse,
    CatalogArchiveRequest,
    CatalogCreateExperimentRequest,
    CatalogCreateTermRequest,
    CatalogErrorResponse,
    CatalogExperimentListResponse,
    CatalogExperimentResponse,
    CatalogExperimentRevisionListResponse,
    CatalogExperimentRunListResponse,
    CatalogPreviewProposalRequest,
    CatalogProposalPreviewResponse,
    CatalogProposalResponse,
    CatalogPublicationExperimentTraceResponse,
    CatalogRegisterProposalRequest,
    CatalogReviseExperimentRequest,
    CatalogReviseSampleRequest,
    CatalogReviseTermRequest,
    CatalogSampleListResponse,
    CatalogSampleRevisionListResponse,
    CatalogSampleRevisionResponse,
    CatalogSampleRunListResponse,
    CatalogTermListResponse,
    CatalogTermResponse,
)
from encode_pipeline.api.dependencies import (
    require_administrator,
    require_principal,
)
from encode_pipeline.api.models import ValidationResponse
from encode_pipeline.platform.authentication import AuthenticatedPrincipal
from encode_pipeline.services.experiment_catalog import ExperimentCatalogService


router = APIRouter(prefix="/catalog", tags=["catalog"])
# OpenAPI array query parameters repeat keys. anyOf [array, null] makes the
# generated client join values with a comma, which the filter rejects.
CatalogTermIds = Annotated[
    list[str] | None,
    Query(),
    WithJsonSchema({"type": "array", "items": {"type": "string"}}),
]
_ERRORS = {
    400: {"model": CatalogErrorResponse},
    403: {"model": CatalogErrorResponse},
    404: {"model": CatalogErrorResponse},
    409: {"model": CatalogErrorResponse},
    422: {"model": CatalogErrorResponse},
    500: {"model": ValidationResponse},
}
_M2_ERRORS = {
    **_ERRORS,
    410: {"model": CatalogErrorResponse},
    413: {"model": CatalogErrorResponse},
}


def get_catalog_service(request: Request) -> ExperimentCatalogService:
    return request.app.state.experiment_catalog_service


def _terms(
    species: list[str] | None,
    tissue: list[str] | None,
    cell_type: list[str] | None,
    cell_line: list[str] | None,
    treatment: list[str] | None,
    tag: list[str] | None,
    assay: list[str] | None = None,
    target: list[str] | None = None,
    batch: list[str] | None = None,
) -> dict[str, list[str] | None]:
    return {
        "species": species,
        "tissue": tissue,
        "cell_type": cell_type,
        "cell_line": cell_line,
        "treatment": treatment,
        "tag": tag,
        "assay": assay,
        "target": target,
        "batch": batch,
    }


def _statuses(
    species_status: str | None,
    tissue_status: str | None,
    cell_type_status: str | None,
    cell_line_status: str | None,
    treatment_status: str | None,
    tag_status: str | None,
    assay_status: str | None = None,
    target_status: str | None = None,
    batch_status: str | None = None,
) -> dict[str, str | None]:
    return {
        "species_status": species_status,
        "tissue_status": tissue_status,
        "cell_type_status": cell_type_status,
        "cell_line_status": cell_line_status,
        "treatment_status": treatment_status,
        "tag_status": tag_status,
        "assay_status": assay_status,
        "target_status": target_status,
        "batch_status": batch_status,
    }


def _field_statuses(values: dict[str, str | None]) -> dict[str, str | None]:
    return {key.removesuffix("_status"): value for key, value in values.items()}


def _page_query(
    *,
    cursor: str | None,
    limit: int | None,
    revision_view: str | None,
    include_archived: bool,
    term_type: str | None,
    created_from: datetime | None,
    created_before: datetime | None,
    catalog_state: str | None,
    field_terms: dict[str, list[str] | None],
    field_statuses: dict[str, str | None],
    project_id: str | None = None,
) -> dict[str, Any]:
    result = {
        "cursor": cursor,
        "limit": limit,
        "revision_view": revision_view,
        "include_archived": include_archived,
        "term_type": term_type,
        "created_from": created_from,
        "created_before": created_before,
        "catalog_state": catalog_state,
        "field_terms": field_terms,
        "field_statuses": _field_statuses(field_statuses),
    }
    if project_id is not None:
        result["project_id"] = project_id
    return result


async def _raw_body(request: Request) -> bytes:
    return await request.body()


def _idempotency_key(
    idempotency_key: Annotated[
        str,
        Header(
            alias="Idempotency-Key",
            min_length=1,
            max_length=128,
            pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$",
            description="Required key for one catalog write.",
        ),
    ],
) -> str:
    return idempotency_key


@router.get(
    "/projects/{project_id}/experiments",
    response_model=CatalogExperimentListResponse,
    operation_id="catalogListExperiments",
    responses=_ERRORS,
)
async def list_experiments(
    project_id: str,
    service: ExperimentCatalogService = Depends(get_catalog_service),
    cursor: str | None = Query(default=None),
    limit: int | None = Query(default=None),
    revision_view: str | None = Query(default=None),
    include_archived: bool = Query(default=False),
    created_from: datetime | None = Query(default=None),
    created_before: datetime | None = Query(default=None),
    species_term_id: CatalogTermIds = None,
    tissue_term_id: CatalogTermIds = None,
    cell_type_term_id: CatalogTermIds = None,
    cell_line_term_id: CatalogTermIds = None,
    treatment_term_id: CatalogTermIds = None,
    tag_term_id: CatalogTermIds = None,
    assay_term_id: CatalogTermIds = None,
    target_term_id: CatalogTermIds = None,
    batch_term_id: CatalogTermIds = None,
    species_status: str | None = Query(default=None),
    tissue_status: str | None = Query(default=None),
    cell_type_status: str | None = Query(default=None),
    cell_line_status: str | None = Query(default=None),
    treatment_status: str | None = Query(default=None),
    tag_status: str | None = Query(default=None),
    assay_status: str | None = Query(default=None),
    target_status: str | None = Query(default=None),
    batch_status: str | None = Query(default=None),
) -> JSONResponse:
    return JSONResponse(
        content=service.list_experiments(
            **_page_query(
                cursor=cursor,
                limit=limit,
                revision_view=revision_view,
                include_archived=include_archived,
                term_type=None,
                created_from=created_from,
                created_before=created_before,
                catalog_state=None,
                field_terms=_terms(
                    species_term_id,
                    tissue_term_id,
                    cell_type_term_id,
                    cell_line_term_id,
                    treatment_term_id,
                    tag_term_id,
                    assay_term_id,
                    target_term_id,
                    batch_term_id,
                ),
                field_statuses=_statuses(
                    species_status,
                    tissue_status,
                    cell_type_status,
                    cell_line_status,
                    treatment_status,
                    tag_status,
                    assay_status,
                    target_status,
                    batch_status,
                ),
                project_id=project_id,
            )
        )
    )


@router.get(
    "/projects/{project_id}/samples",
    response_model=CatalogSampleListResponse,
    operation_id="catalogListSamples",
    responses=_ERRORS,
)
async def list_samples(
    project_id: str,
    service: ExperimentCatalogService = Depends(get_catalog_service),
    cursor: str | None = Query(default=None),
    limit: int | None = Query(default=None),
    created_from: datetime | None = Query(default=None),
    created_before: datetime | None = Query(default=None),
    catalog_state: str | None = Query(default=None),
    species_term_id: CatalogTermIds = None,
    tissue_term_id: CatalogTermIds = None,
    cell_type_term_id: CatalogTermIds = None,
    cell_line_term_id: CatalogTermIds = None,
    treatment_term_id: CatalogTermIds = None,
    tag_term_id: CatalogTermIds = None,
    species_status: str | None = Query(default=None),
    tissue_status: str | None = Query(default=None),
    cell_type_status: str | None = Query(default=None),
    cell_line_status: str | None = Query(default=None),
    treatment_status: str | None = Query(default=None),
    tag_status: str | None = Query(default=None),
) -> JSONResponse:
    return JSONResponse(
        content=service.list_samples(
            **_page_query(
                cursor=cursor,
                limit=limit,
                revision_view=None,
                include_archived=False,
                term_type=None,
                created_from=created_from,
                created_before=created_before,
                catalog_state=catalog_state,
                field_terms=_terms(
                    species_term_id,
                    tissue_term_id,
                    cell_type_term_id,
                    cell_line_term_id,
                    treatment_term_id,
                    tag_term_id,
                ),
                field_statuses=_statuses(
                    species_status,
                    tissue_status,
                    cell_type_status,
                    cell_line_status,
                    treatment_status,
                    tag_status,
                ),
                project_id=project_id,
            )
        )
    )


@router.get(
    "/samples/{sample_id}/revisions",
    response_model=CatalogSampleRevisionListResponse,
    operation_id="catalogListSampleRevisions",
    responses=_ERRORS,
)
async def list_sample_revisions(
    sample_id: str,
    service: ExperimentCatalogService = Depends(get_catalog_service),
    cursor: str | None = Query(default=None),
    limit: int | None = Query(default=None),
) -> JSONResponse:
    return JSONResponse(
        content=service.list_sample_revisions(
            sample_id=sample_id,
            **_page_query(
                cursor=cursor,
                limit=limit,
                revision_view=None,
                include_archived=False,
                term_type=None,
                created_from=None,
                created_before=None,
                catalog_state=None,
                field_terms={},
                field_statuses={},
            ),
        )
    )


@router.get(
    "/samples/{sample_id}/revisions/{sample_revision_id}",
    response_model=CatalogSampleRevisionResponse,
    operation_id="catalogGetSampleRevision",
    responses=_ERRORS,
)
async def get_sample_revision(
    sample_id: str,
    sample_revision_id: str,
    _principal: AuthenticatedPrincipal = Depends(require_administrator),
    service: ExperimentCatalogService = Depends(get_catalog_service),
) -> JSONResponse:
    return JSONResponse(
        content=service.get_sample_revision(sample_id, sample_revision_id)
    )


@router.get(
    "/experiments/{experiment_id}/revisions",
    response_model=CatalogExperimentRevisionListResponse,
    operation_id="catalogListExperimentRevisions",
    responses=_ERRORS,
)
async def list_experiment_revisions(
    experiment_id: str,
    service: ExperimentCatalogService = Depends(get_catalog_service),
    cursor: str | None = Query(default=None),
    limit: int | None = Query(default=None),
    revision_view: str | None = Query(default=None),
) -> JSONResponse:
    return JSONResponse(
        content=service.list_experiment_revisions(
            experiment_id=experiment_id,
            **_page_query(
                cursor=cursor,
                limit=limit,
                revision_view=revision_view,
                include_archived=False,
                term_type=None,
                created_from=None,
                created_before=None,
                catalog_state=None,
                field_terms={},
                field_statuses={},
            ),
        )
    )


@router.get(
    "/experiments/{experiment_id}/revisions/{experiment_revision_id}",
    response_model=CatalogExperimentResponse,
    operation_id="catalogGetExperimentRevision",
    responses=_ERRORS,
)
async def get_experiment_revision(
    experiment_id: str,
    experiment_revision_id: str,
    _principal: AuthenticatedPrincipal = Depends(require_administrator),
    service: ExperimentCatalogService = Depends(get_catalog_service),
) -> JSONResponse:
    return JSONResponse(
        content=service.get_experiment_revision(experiment_id, experiment_revision_id)
    )


@router.get(
    "/terms",
    response_model=CatalogTermListResponse,
    operation_id="catalogListTerms",
    responses=_ERRORS,
)
async def list_terms(
    service: ExperimentCatalogService = Depends(get_catalog_service),
    cursor: str | None = Query(default=None),
    limit: int | None = Query(default=None),
    include_archived: bool = Query(default=False),
    term_type: str | None = Query(default=None),
    created_from: datetime | None = Query(default=None),
    created_before: datetime | None = Query(default=None),
) -> JSONResponse:
    return JSONResponse(
        content=service.list_terms(
            **_page_query(
                cursor=cursor,
                limit=limit,
                revision_view=None,
                include_archived=include_archived,
                term_type=term_type,
                created_from=created_from,
                created_before=created_before,
                catalog_state=None,
                field_terms={},
                field_statuses={},
            )
        )
    )


@router.get(
    "/terms/{term_id}",
    response_model=CatalogTermResponse,
    operation_id="catalogGetTerm",
    responses=_ERRORS,
)
async def get_term(
    term_id: str,
    _principal: AuthenticatedPrincipal = Depends(require_administrator),
    service: ExperimentCatalogService = Depends(get_catalog_service),
) -> JSONResponse:
    return JSONResponse(content=service.get_term(term_id))


@router.post(
    "/terms",
    response_model=CatalogTermResponse,
    operation_id="catalogCreateTerm",
    responses=_ERRORS,
)
async def create_term(
    request: Request,
    body: CatalogCreateTermRequest,
    idempotency_key: str = Depends(_idempotency_key),
    principal: AuthenticatedPrincipal = Depends(require_administrator),
    service: ExperimentCatalogService = Depends(get_catalog_service),
) -> JSONResponse:
    result = service.create_term(
        actor_user_id=principal.user_id,
        idempotency_key=idempotency_key,
        raw_body=await _raw_body(request),
        term_type=body.term_type,
        display_name=body.display_name,
        namespace=body.namespace,
        accession=body.accession,
        aliases=body.aliases,
    )
    return JSONResponse(content=result)


@router.post(
    "/terms/{term_id}/revisions",
    response_model=CatalogTermResponse,
    operation_id="catalogReviseTerm",
    responses=_ERRORS,
)
async def revise_term(
    term_id: str,
    request: Request,
    body: CatalogReviseTermRequest,
    idempotency_key: str = Depends(_idempotency_key),
    principal: AuthenticatedPrincipal = Depends(require_administrator),
    service: ExperimentCatalogService = Depends(get_catalog_service),
) -> JSONResponse:
    result = service.revise_term(
        actor_user_id=principal.user_id,
        idempotency_key=idempotency_key,
        raw_body=await _raw_body(request),
        term_id=term_id,
        expected_revision_number=body.expected_revision_number,
        change_kind=body.change_kind,
        display_name=body.display_name,
        namespace=body.namespace,
        accession=body.accession,
        aliases=body.aliases,
    )
    return JSONResponse(content=result)


@router.post(
    "/terms/{term_id}/archive",
    response_model=CatalogTermResponse,
    operation_id="catalogArchiveTerm",
    responses=_ERRORS,
)
async def archive_term(
    term_id: str,
    request: Request,
    body: CatalogArchiveRequest,
    idempotency_key: str = Depends(_idempotency_key),
    principal: AuthenticatedPrincipal = Depends(require_administrator),
    service: ExperimentCatalogService = Depends(get_catalog_service),
) -> JSONResponse:
    result = service.archive_term(
        actor_user_id=principal.user_id,
        idempotency_key=idempotency_key,
        raw_body=await _raw_body(request),
        term_id=term_id,
        expected_revision_number=body.expected_revision_number,
    )
    return JSONResponse(content=result)


@router.post(
    "/projects/{project_id}/experiments",
    response_model=CatalogExperimentResponse,
    operation_id="catalogCreateExperiment",
    responses=_ERRORS,
)
async def create_experiment(
    project_id: str,
    request: Request,
    body: CatalogCreateExperimentRequest,
    idempotency_key: str = Depends(_idempotency_key),
    principal: AuthenticatedPrincipal = Depends(require_administrator),
    service: ExperimentCatalogService = Depends(get_catalog_service),
) -> JSONResponse:
    result = service.create_experiment(
        actor_user_id=principal.user_id,
        idempotency_key=idempotency_key,
        raw_body=await _raw_body(request),
        project_id=project_id,
        stable_key=body.stable_key,
        display_name=body.display_name,
        assay=body.assay,
        target=body.target,
        batch=body.batch,
        samples=body.samples,
    )
    return JSONResponse(content=result)


@router.post(
    "/experiments/{experiment_id}/revisions",
    response_model=CatalogExperimentResponse,
    operation_id="catalogReviseExperiment",
    responses=_ERRORS,
)
async def revise_experiment(
    experiment_id: str,
    request: Request,
    body: CatalogReviseExperimentRequest,
    idempotency_key: str = Depends(_idempotency_key),
    principal: AuthenticatedPrincipal = Depends(require_administrator),
    service: ExperimentCatalogService = Depends(get_catalog_service),
) -> JSONResponse:
    result = service.revise_experiment(
        actor_user_id=principal.user_id,
        idempotency_key=idempotency_key,
        raw_body=await _raw_body(request),
        experiment_id=experiment_id,
        expected_revision_number=body.expected_revision_number,
        display_name=body.display_name,
        assay=body.assay,
        target=body.target,
        batch=body.batch,
        samples=body.samples,
    )
    return JSONResponse(content=result)


@router.post(
    "/experiments/{experiment_id}/archive",
    response_model=CatalogArchiveExperimentResponse,
    operation_id="catalogArchiveExperiment",
    responses=_ERRORS,
)
async def archive_experiment(
    experiment_id: str,
    request: Request,
    body: CatalogArchiveRequest,
    idempotency_key: str = Depends(_idempotency_key),
    principal: AuthenticatedPrincipal = Depends(require_administrator),
    service: ExperimentCatalogService = Depends(get_catalog_service),
) -> JSONResponse:
    result = service.archive_experiment(
        actor_user_id=principal.user_id,
        idempotency_key=idempotency_key,
        raw_body=await _raw_body(request),
        experiment_id=experiment_id,
        expected_revision_number=body.expected_revision_number,
    )
    return JSONResponse(content=result)


@router.post(
    "/samples/{sample_id}/catalog-revisions",
    response_model=CatalogSampleRevisionResponse,
    operation_id="catalogReviseSampleCatalog",
    responses=_ERRORS,
)
async def revise_sample_catalog(
    sample_id: str,
    request: Request,
    body: CatalogReviseSampleRequest,
    idempotency_key: str = Depends(_idempotency_key),
    principal: AuthenticatedPrincipal = Depends(require_administrator),
    service: ExperimentCatalogService = Depends(get_catalog_service),
) -> JSONResponse:
    result = service.revise_sample_catalog(
        actor_user_id=principal.user_id,
        idempotency_key=idempotency_key,
        raw_body=await _raw_body(request),
        sample_id=sample_id,
        expected_revision_number=body.expected_revision_number,
        display_name=body.display_name,
        catalog=body.catalog,
    )
    return JSONResponse(content=result)


@router.get(
    "/projects/{project_id}/sample-revisions/{sample_revision_id}/run-results",
    response_model=CatalogSampleRunListResponse,
    operation_id="catalogListSampleRevisionRunResults",
    responses=_ERRORS,
)
async def list_sample_run_results(
    project_id: str,
    sample_revision_id: str,
    service: ExperimentCatalogService = Depends(get_catalog_service),
    cursor: str | None = Query(default=None),
    limit: int | None = Query(default=None),
) -> JSONResponse:
    return JSONResponse(
        content=service.list_sample_run_results(
            project_id=project_id,
            sample_revision_id=sample_revision_id,
            **_page_query(
                cursor=cursor,
                limit=limit,
                revision_view=None,
                include_archived=False,
                term_type=None,
                created_from=None,
                created_before=None,
                catalog_state=None,
                field_terms={},
                field_statuses={},
            ),
        )
    )


@router.get(
    "/projects/{project_id}/experiment-revisions/{experiment_revision_id}/run-results",
    response_model=CatalogExperimentRunListResponse,
    operation_id="catalogListExperimentRevisionRunResults",
    responses=_ERRORS,
)
async def list_experiment_run_results(
    project_id: str,
    experiment_revision_id: str,
    service: ExperimentCatalogService = Depends(get_catalog_service),
    cursor: str | None = Query(default=None),
    limit: int | None = Query(default=None),
) -> JSONResponse:
    return JSONResponse(
        content=service.list_experiment_run_results(
            project_id=project_id,
            experiment_revision_id=experiment_revision_id,
            **_page_query(
                cursor=cursor,
                limit=limit,
                revision_view=None,
                include_archived=False,
                term_type=None,
                created_from=None,
                created_before=None,
                catalog_state=None,
                field_terms={},
                field_statuses={},
            ),
        )
    )


@router.post(
    "/projects/{project_id}/proposals",
    response_model=CatalogProposalResponse,
    operation_id="catalogRegisterProposal",
    responses=_M2_ERRORS,
)
async def register_proposal(
    project_id: str,
    request: Request,
    body: CatalogRegisterProposalRequest,
    idempotency_key: str = Depends(_idempotency_key),
    principal: AuthenticatedPrincipal = Depends(require_administrator),
    service: ExperimentCatalogService = Depends(get_catalog_service),
) -> JSONResponse:
    result = service.register_proposal(
        actor_user_id=principal.user_id,
        idempotency_key=idempotency_key,
        raw_body=await _raw_body(request),
        project_id=project_id,
        source=body.source,
        items=body.items,
        rule_version=body.rule_version,
        vocabulary_version=body.vocabulary_version,
        projection_version=body.projection_version,
    )
    return JSONResponse(content=result)


@router.get(
    "/proposals/{proposal_id}",
    response_model=CatalogProposalResponse,
    operation_id="catalogGetProposal",
    responses=_M2_ERRORS,
)
async def get_proposal(
    proposal_id: str,
    principal: AuthenticatedPrincipal = Depends(require_principal),
    service: ExperimentCatalogService = Depends(get_catalog_service),
) -> JSONResponse:
    return JSONResponse(
        content=service.get_proposal(
            proposal_id=proposal_id,
            include_candidates=principal.is_administrator,
        )
    )


@router.post(
    "/proposals/{proposal_id}/preview",
    response_model=CatalogProposalPreviewResponse,
    operation_id="catalogPreviewProposal",
    responses=_M2_ERRORS,
)
async def preview_proposal(
    proposal_id: str,
    request: Request,
    body: CatalogPreviewProposalRequest,
    principal: AuthenticatedPrincipal = Depends(require_administrator),
    service: ExperimentCatalogService = Depends(get_catalog_service),
) -> JSONResponse:
    result = service.preview_proposal(
        actor_user_id=principal.user_id,
        raw_body=await _raw_body(request),
        proposal_id=proposal_id,
        selected_item_ids=body.selected_item_ids,
    )
    return JSONResponse(content=result)


@router.post(
    "/proposals/{proposal_id}/adoptions",
    response_model=CatalogAdoptionReceiptResponse,
    operation_id="catalogAdoptProposal",
    responses=_M2_ERRORS,
)
async def adopt_proposal(
    proposal_id: str,
    request: Request,
    body: CatalogAdoptProposalRequest,
    idempotency_key: str = Depends(_idempotency_key),
    principal: AuthenticatedPrincipal = Depends(require_administrator),
    service: ExperimentCatalogService = Depends(get_catalog_service),
) -> JSONResponse:
    result = service.adopt_proposal(
        actor_user_id=principal.user_id,
        idempotency_key=idempotency_key,
        raw_body=await _raw_body(request),
        proposal_id=proposal_id,
        selected_item_ids=body.selected_item_ids,
        expected_revisions=body.expected_revisions,
        preview_digest=body.preview_digest,
    )
    return JSONResponse(content=result)


@router.get(
    "/projects/{project_id}/adoption-receipts",
    response_model=CatalogAdoptionReceiptListResponse,
    operation_id="catalogListAdoptionReceipts",
    responses=_M2_ERRORS,
)
async def list_adoption_receipts(
    project_id: str,
    principal: AuthenticatedPrincipal = Depends(require_administrator),
    service: ExperimentCatalogService = Depends(get_catalog_service),
    proposal_id: str | None = Query(default=None),
) -> JSONResponse:
    return JSONResponse(
        content=service.list_adoption_receipts(
            actor_user_id=principal.user_id,
            project_id=project_id,
            proposal_id=proposal_id,
        )
    )
@router.get(
    "/projects/{project_id}/artifact-publications/"
    "{run_id}/{artifact_id}/experiment-trace",
    response_model=CatalogPublicationExperimentTraceResponse,
    operation_id="catalogGetArtifactPublicationExperimentTrace",
    responses=_ERRORS,
)
async def get_publication_experiment_trace(
    project_id: str,
    run_id: str,
    artifact_id: str,
    generation: str = Query(
        min_length=76,
        max_length=76,
        pattern=r"^artifactgen-[0-9a-f]{64}$",
    ),
    service: ExperimentCatalogService = Depends(get_catalog_service),
) -> JSONResponse:
    """Return the experiment revision copied onto this publication's run.

    The copied revision id and payload digest are compared with that stored
    revision. A newer experiment revision does not replace them. A missing
    binding, or a binding that does not match the stored revision, returns
    the published row with an unresolved reason. The row stays
    ``run_input_scope``.
    """

    # This repository call opens SQLite. Keep it off the event loop.
    content = await run_in_threadpool(
        service.get_publication_experiment_trace,
        project_id=project_id,
        run_id=run_id,
        artifact_id=artifact_id,
        artifact_generation=generation,
    )
    return JSONResponse(content=content)
