"""API v1 router aggregation."""

from __future__ import annotations

from fastapi import APIRouter

from encode_pipeline.api.routes.agent import router as agent_router
from encode_pipeline.api.routes.workflows import router as workflows_router


api_v1_router = APIRouter()
api_v1_router.include_router(workflows_router)
api_v1_router.include_router(agent_router)
