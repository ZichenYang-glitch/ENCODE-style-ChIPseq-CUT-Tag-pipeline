"""Read-only artifact-publication API contract and disclosure tests."""

from __future__ import annotations

import asyncio
from datetime import datetime, timezone
import inspect
from threading import get_ident

import httpx
from pydantic import ValidationError
import pytest

from api_test_client import ApiTestClient
from encode_pipeline.api.main import create_app
from encode_pipeline.api.models import ArtifactPublicationResponse
from encode_pipeline.api.routes.artifact_publications import (
    get_artifact_publication,
    list_artifact_publications,
)
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.artifact_publications import (
    ArtifactPublicationPage,
    ArtifactPublicationRunSampleBinding,
    ArtifactPublicationSummary,
    AssociatedRunSample,
)
from encode_pipeline.platform.data_registry import BindingMode, BindingProvenance
from encode_pipeline.platform.runs import RunArtifactRef, RunStatus
from encode_pipeline.services.artifact_publications import (
    ArtifactPublicationCursorInvalidError,
    ArtifactPublicationCursorNotFoundError,
    ArtifactPublicationDataInvalidError,
    ArtifactPublicationNotFoundError,
)


NOW = datetime(2026, 8, 7, 8, 0, tzinfo=timezone.utc)
PROJECT_ID = f"prj_{'1' * 32}"
SAMPLE_ID = f"smp_{'2' * 32}"
SAMPLE_REVISION_ID = f"smpr_{'3' * 32}"
GENERATION = f"artifactgen-{'a' * 64}"
OLD_GENERATION = f"artifactgen-{'b' * 64}"
REVISION = f"artifactrev-{'c' * 64}"


@pytest.fixture
def client():
    app = create_app()
    with ApiTestClient(app) as test_client:
        yield test_client


def _summary(*, current_generation: str = GENERATION) -> ArtifactPublicationSummary:
    return ArtifactPublicationSummary(
        run_id="run-a",
        project_id=PROJECT_ID,
        workflow_id="workflow-a",
        artifact_id="artifact-a",
        output_type="counts",
        artifact_generation=GENERATION,
        artifact_revision=REVISION,
        published_at=NOW,
        current_artifact_generation=current_generation,
        run_sample_binding=ArtifactPublicationRunSampleBinding(
            binding_mode=BindingMode.BOUND_V1,
            provenance=BindingProvenance.RESOLVED,
            associated_run_samples=(
                AssociatedRunSample(
                    sample_id=SAMPLE_ID,
                    sample_revision_id=SAMPLE_REVISION_ID,
                    revision_number=4,
                    ordinal=0,
                ),
            ),
        ),
    )


def _persisted_artifact(
    run_id: str,
    artifact_id: str,
    *,
    revision: str,
) -> RunArtifactRef:
    return RunArtifactRef(
        artifact_id=artifact_id,
        run_id=run_id,
        artifact_type="file",
        name=f"{artifact_id}.txt",
        uri=f"run://runs/{run_id}/artifacts/{artifact_id}",
        mime_type="text/plain",
        produced_at=NOW,
        revision=revision,
        metadata={
            "relative_path": f"results/private/{artifact_id}.txt",
            "output_type": "summary",
            "size_bytes": 7,
            "catalog_id": "summary",
            "scope": "project",
            "private_token": "SECRET_PERSISTED_METADATA",
        },
    )


def _create_succeeded_run(app) -> str:
    service = app.state.run_service
    run_id = service.create_run(
        "encode-style-chipseq-cuttag-atac-mnase",
        WorkflowInputs(config={"private_input": "/private/input.fastq.gz"}),
    ).run_id
    for status in (
        RunStatus.VALIDATING,
        RunStatus.PLANNED,
        RunStatus.QUEUED,
        RunStatus.RUNNING,
        RunStatus.SUCCEEDED,
    ):
        service.transition_run(run_id, status)
    return run_id


class _Service:
    private_path = "/private/publication/payload.bam"
    private_metadata = {"token": "SECRET_TOKEN"}

    def __init__(self, summary: ArtifactPublicationSummary | None = None) -> None:
        self.summary = summary
        self.list_kwargs = None

    def list_artifact_publications(self, **kwargs):
        self.list_kwargs = kwargs
        items = () if self.summary is None else (self.summary,)
        return ArtifactPublicationPage(
            artifact_publications=items,
            next_cursor=None,
        )

    def get_artifact_publication(self, **_kwargs):
        if self.summary is None:
            raise ArtifactPublicationNotFoundError
        return self.summary


def test_routes_are_sync_with_stable_operation_ids_and_empty_envelope(client) -> None:
    assert not inspect.iscoroutinefunction(list_artifact_publications)
    assert not inspect.iscoroutinefunction(get_artifact_publication)
    schema = client.get("/openapi.json").json()
    list_operation = schema["paths"]["/api/v1/artifact-publications"]["get"]
    detail_operation = schema["paths"][
        "/api/v1/artifact-publications/{run_id}/{artifact_id}"
    ]["get"]
    assert list_operation["operationId"] == "listArtifactPublications"
    assert detail_operation["operationId"] == "getArtifactPublication"

    client.app.state.artifact_publication_service = _Service()
    response = client.get("/api/v1/artifact-publications")

    assert response.status_code == 200
    assert response.json() == {
        "ok": True,
        "publications": [],
        "next_cursor": None,
        "issues": [],
    }


def test_list_projects_only_allowlisted_path_free_fields(client) -> None:
    service = _Service(_summary())
    client.app.state.artifact_publication_service = service

    response = client.get(
        "/api/v1/artifact-publications",
        params={
            "project_id": PROJECT_ID,
            "run_id": "run-a",
            "workflow_id": "workflow-a",
            "output_type": "counts",
            "associated_run_sample_revision_id": SAMPLE_REVISION_ID,
            "published_from": "2026-08-07T08:00:00Z",
            "published_before": "2026-08-07T09:00:00Z",
            "current_only": True,
        },
    )

    assert response.status_code == 200
    publication = response.json()["publications"][0]
    assert set(publication) == {
        "run_id",
        "project_id",
        "workflow_id",
        "artifact_id",
        "output_type",
        "artifact_generation",
        "artifact_revision",
        "published_at",
        "current_artifact_generation",
        "generation_status",
        "run_sample_binding",
    }
    assert publication["generation_status"] == "current"
    assert publication["run_sample_binding"] == {
        "binding_mode": "bound_v1",
        "provenance": "resolved",
        "associated_run_samples": [
            {
                "sample_id": SAMPLE_ID,
                "sample_revision_id": SAMPLE_REVISION_ID,
                "revision_number": 4,
                "ordinal": 0,
            }
        ],
    }
    assert service.list_kwargs["current_only"] is True
    assert "/private" not in response.text
    assert "SECRET_TOKEN" not in response.text
    for forbidden in (
        "uri",
        "path",
        "name",
        "mime",
        "metadata",
        "config",
        "environment",
        "display_name",
    ):
        assert forbidden not in publication


def test_detail_uses_exact_generation_and_derives_superseded_status(client) -> None:
    client.app.state.artifact_publication_service = _Service(
        _summary(current_generation=OLD_GENERATION)
    )

    response = client.get(
        "/api/v1/artifact-publications/run-a/artifact-a",
        params={"generation": GENERATION},
    )

    assert response.status_code == 200
    publication = response.json()["publication"]
    assert publication["run_id"] == "run-a"
    assert publication["artifact_generation"] == GENERATION
    assert publication["generation_status"] == "superseded"


@pytest.mark.parametrize(
    ("url", "params"),
    [
        ("/api/v1/artifact-publications", {"limit": 0}),
        ("/api/v1/artifact-publications", {"current_only": "maybe"}),
        ("/api/v1/artifact-publications", {"published_from": "not-a-date"}),
        ("/api/v1/artifact-publications", {"after": "not-a-cursor"}),
        ("/api/v1/artifact-publications/run-a/artifact-a", {}),
        (
            "/api/v1/artifact-publications/run-a/artifact-a",
            {"generation": "artifactgen-invalid"},
        ),
    ],
)
def test_invalid_requests_use_stable_400_envelopes_without_422(
    client,
    url,
    params,
) -> None:
    response = client.get(url, params=params)

    assert response.status_code == 400
    assert response.json()["ok"] is False
    assert response.json()["issues"][0]["code"] == "API_REQUEST_INVALID"
    assert "technical_message" not in response.text


@pytest.mark.parametrize(
    ("error", "code"),
    [
        (
            ArtifactPublicationCursorInvalidError(),
            "ARTIFACT_PUBLICATION_CURSOR_INVALID",
        ),
        (
            ArtifactPublicationCursorNotFoundError(),
            "ARTIFACT_PUBLICATION_CURSOR_NOT_FOUND",
        ),
        (
            ArtifactPublicationDataInvalidError(),
            "ARTIFACT_PUBLICATION_DATA_INVALID",
        ),
    ],
)
def test_list_maps_safe_publication_errors(client, error, code) -> None:
    class _FailingService(_Service):
        def list_artifact_publications(self, **_kwargs):
            raise error

    client.app.state.artifact_publication_service = _FailingService()
    response = client.get(
        "/api/v1/artifact-publications",
        params={"after": f"artifactpubcur_{'a' * 16}"},
    )

    assert response.status_code == (500 if code.endswith("DATA_INVALID") else 400)
    assert response.json()["publications"] == []
    assert response.json()["next_cursor"] is None
    assert response.json()["issues"][0]["code"] == code
    assert "technical_message" not in response.text


def test_missing_detail_returns_safe_404(client) -> None:
    client.app.state.artifact_publication_service = _Service()

    response = client.get(
        "/api/v1/artifact-publications/run-a/artifact-a",
        params={"generation": GENERATION},
    )

    assert response.status_code == 404
    assert response.json()["publication"] is None
    assert response.json()["issues"][0]["code"] == ("ARTIFACT_PUBLICATION_NOT_FOUND")


def test_repository_work_executes_off_request_event_loop(client) -> None:
    request_thread = get_ident()
    called_threads: list[int] = []

    class _ObservedService(_Service):
        def list_artifact_publications(self, **kwargs):
            called_threads.append(get_ident())
            return super().list_artifact_publications(**kwargs)

    client.app.state.artifact_publication_service = _ObservedService()
    response = client.get("/api/v1/artifact-publications")

    assert response.status_code == 200
    assert called_threads and called_threads[0] != request_thread


@pytest.mark.parametrize(
    "url",
    [
        "/api/v1/artifact-publications",
        f"/api/v1/artifact-publications/run-a/artifact-a?generation={GENERATION}",
    ],
)
def test_unexpected_failure_uses_route_specific_safe_500(client, url) -> None:
    private_text = "sqlite:////private/publications.db SECRET=value"

    class _FailingService(_Service):
        def list_artifact_publications(self, **_kwargs):
            raise RuntimeError(private_text)

        def get_artifact_publication(self, **_kwargs):
            raise RuntimeError(private_text)

    client.app.state.artifact_publication_service = _FailingService()

    async def request() -> httpx.Response:
        transport = httpx.ASGITransport(
            app=client.app,
            raise_app_exceptions=False,
        )
        async with httpx.AsyncClient(
            transport=transport,
            base_url="http://testserver",
        ) as async_client:
            return await async_client.get(url)

    response = asyncio.run(request())

    assert response.status_code == 500
    assert response.json()["ok"] is False
    assert response.json()["issues"][0]["code"] == "INTERNAL_SERVER_ERROR"
    assert private_text not in response.text
    assert "technical_message" not in response.text


def test_publication_response_models_forbid_unexpected_fields() -> None:
    payload = ArtifactPublicationResponse.from_summary(_summary()).model_dump()
    payload["uri"] = "file:///private/payload"

    with pytest.raises(ValidationError):
        ArtifactPublicationResponse.model_validate(payload)


def test_sqlite_publications_survive_reopen_and_cursor_boundary_must_still_match(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    first_app = create_app(database_url=database_url)
    run_id = _create_succeeded_run(first_app)
    original = (
        _persisted_artifact(run_id, "artifact-a", revision=REVISION),
        _persisted_artifact(run_id, "artifact-b", revision=REVISION),
    )
    first_app.state.run_service.replace_artifacts(run_id, original)
    old_generation = first_app.state.run_service.get_result_state(
        run_id
    ).artifact_generation
    assert old_generation is not None
    first_app.state.run_queue.close()
    first_app.state.persistence.close()

    reopened_app = create_app(database_url=database_url)
    with ApiTestClient(reopened_app) as reopened_client:
        first_page = reopened_client.get(
            "/api/v1/artifact-publications",
            params={"run_id": run_id, "limit": 1},
        )
        detail = reopened_client.get(
            f"/api/v1/artifact-publications/{run_id}/artifact-a",
            params={"generation": old_generation},
        )

        assert first_page.status_code == 200
        assert len(first_page.json()["publications"]) == 1
        cursor = first_page.json()["next_cursor"]
        assert cursor.startswith("artifactpubcur_")
        publication = first_page.json()["publications"][0]
        assert publication["project_id"] == f"prj_{'0' * 32}"
        assert publication["generation_status"] == "current"
        assert publication["run_sample_binding"] == {
            "binding_mode": "legacy_v1",
            "provenance": "unresolved",
            "associated_run_samples": [],
        }
        assert detail.status_code == 200
        assert detail.json()["publication"]["artifact_generation"] == (old_generation)
        assert "/private" not in first_page.text
        assert "SECRET_PERSISTED_METADATA" not in first_page.text

        replacement_revision = f"artifactrev-{'d' * 64}"
        replacement = tuple(
            _persisted_artifact(
                run_id,
                artifact.artifact_id,
                revision=replacement_revision,
            )
            for artifact in original
        )
        reopened_app.state.run_service.replace_artifacts(run_id, replacement)

        stale = reopened_client.get(
            "/api/v1/artifact-publications",
            params={"run_id": run_id, "limit": 1, "after": cursor},
        )
        history = reopened_client.get(
            "/api/v1/artifact-publications",
            params={"run_id": run_id, "current_only": False},
        )
        historical_detail = reopened_client.get(
            f"/api/v1/artifact-publications/{run_id}/artifact-a",
            params={"generation": old_generation},
        )

    assert stale.status_code == 400
    assert stale.json()["issues"][0]["code"] == (
        "ARTIFACT_PUBLICATION_CURSOR_NOT_FOUND"
    )
    assert len(history.json()["publications"]) == 4
    assert {item["generation_status"] for item in history.json()["publications"]} == {
        "current",
        "superseded",
    }
    assert historical_detail.status_code == 200
    assert historical_detail.json()["publication"]["generation_status"] == (
        "superseded"
    )
    reopened_app.state.run_queue.close()
    reopened_app.state.persistence.close()
