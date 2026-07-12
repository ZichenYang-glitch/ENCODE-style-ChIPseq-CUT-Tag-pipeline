"""API v1 router aggregation."""

from __future__ import annotations

from fastapi import APIRouter

from encode_pipeline.api.routes.agent import router as agent_router
from encode_pipeline.api.routes.artifacts import router as artifacts_router
from encode_pipeline.api.routes.preflight import router as preflight_router
from encode_pipeline.api.routes.qc_metrics import router as qc_metrics_router
from encode_pipeline.api.routes.runs import router as runs_router
from encode_pipeline.api.routes.workflows import router as workflows_router


api_v1_router = APIRouter()
api_v1_router.include_router(workflows_router)
api_v1_router.include_router(agent_router)
api_v1_router.include_router(runs_router)
api_v1_router.include_router(preflight_router)
api_v1_router.include_router(artifacts_router)
api_v1_router.include_router(qc_metrics_router)
