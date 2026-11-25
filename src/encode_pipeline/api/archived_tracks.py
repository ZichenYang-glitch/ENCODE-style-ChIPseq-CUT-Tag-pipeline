"""Public models for experiment-scoped sealed track reads."""

from __future__ import annotations

from typing import Literal

from fastapi import Request
from fastapi.responses import JSONResponse
from pydantic import BaseModel, ConfigDict, Field

from encode_pipeline.api.catalog_models import (
    CatalogErrorResponse,
    CatalogIssueResponse,
)


ARCHIVED_TRACK_OPERATIONS = frozenset(
    {
        "catalogListArchivedTracks",
        "catalogGetArchivedTrack",
        "catalogReadArchivedTrackContent",
        "catalogHeadArchivedTrackContent",
        "catalogReadArchivedTrackIndex",
        "catalogHeadArchivedTrackIndex",
    }
)

PRIVATE_NO_STORE = "private, no-store"


class ArchivedTrackRunSource(BaseModel):
    """Run publication pinned by the seal. Paths and digests stay off this object."""

    model_config = ConfigDict(extra="forbid")

    kind: Literal["run"]
    run_id: str
    artifact_id: str
    generation: str
    revision: str
    index_artifact_id: str | None = None
    index_generation: str | None = None
    index_revision: str | None = None


class ArchivedTrackRegisteredSource(BaseModel):
    """Registered revision pinned by the seal. Paths and digests stay off this object."""

    model_config = ConfigDict(extra="forbid")

    kind: Literal["registered"]
    result_id: str
    result_revision_id: str
    index_result_id: str | None = None
    index_result_revision_id: str | None = None


class ArchivedTrackResponse(BaseModel):
    """One sealed track version. List items and the detail resource share this shape."""

    model_config = ConfigDict(extra="forbid")

    archive_version_id: str
    project_id: str
    experiment_id: str
    experiment_revision_id: str
    track_key: str
    version_number: int = Field(ge=1)
    display_name: str
    assembly: str
    track_type: Literal["bigwig", "bed", "longrange"]
    source: ArchivedTrackRunSource | ArchivedTrackRegisteredSource
    has_index: bool


class ArchivedTrackListResponse(BaseModel):
    """One page of sealed versions for an experiment."""

    model_config = ConfigDict(extra="forbid")

    items: list[ArchivedTrackResponse]
    next_offset: int | None = None


def archived_track_operation(request: Request) -> bool:
    """Return whether the matched route is an archived-track read."""

    route = request.scope.get("route")
    operation_id = getattr(route, "operation_id", None)
    if operation_id in ARCHIVED_TRACK_OPERATIONS:
        return True
    return "/archived-tracks" in request.url.path


def archived_track_json(status_code: int, code: str, message: str) -> JSONResponse:
    """Return a catalog failure that is not stored by shared caches."""

    body = CatalogErrorResponse(
        issues=[
            CatalogIssueResponse(
                code=code,
                message=message,
                severity="error",
                path=None,
                source="api",
                hint=None,
            )
        ]
    )
    return JSONResponse(
        status_code=status_code,
        content=body.model_dump(mode="json"),
        headers={"Cache-Control": PRIVATE_NO_STORE},
    )
