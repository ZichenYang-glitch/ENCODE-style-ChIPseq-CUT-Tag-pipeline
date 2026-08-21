"""API v1 router aggregation."""

from __future__ import annotations

from fastapi import APIRouter, Depends

from encode_pipeline.api.dependencies import enforce_csrf, require_principal
from encode_pipeline.api.routes.agent import router as agent_router
from encode_pipeline.api.routes.artifacts import router as artifacts_router
from encode_pipeline.api.routes.artifact_publications import (
    router as artifact_publications_router,
)
from encode_pipeline.api.routes.auth import router as auth_router
from encode_pipeline.api.routes.preflight import router as preflight_router
from encode_pipeline.api.routes.qc_metrics import router as qc_metrics_router
from encode_pipeline.api.routes.runs import router as runs_router
from encode_pipeline.api.routes.workflows import router as workflows_router


api_v1_router = APIRouter()
api_v1_router.include_router(auth_router)
_member_dependencies = [Depends(require_principal), Depends(enforce_csrf)]
api_v1_router.include_router(workflows_router, dependencies=_member_dependencies)
api_v1_router.include_router(agent_router, dependencies=_member_dependencies)
api_v1_router.include_router(runs_router, dependencies=_member_dependencies)
api_v1_router.include_router(preflight_router, dependencies=_member_dependencies)
api_v1_router.include_router(artifacts_router, dependencies=_member_dependencies)
api_v1_router.include_router(
    artifact_publications_router,
    dependencies=_member_dependencies,
)
api_v1_router.include_router(qc_metrics_router, dependencies=_member_dependencies)
