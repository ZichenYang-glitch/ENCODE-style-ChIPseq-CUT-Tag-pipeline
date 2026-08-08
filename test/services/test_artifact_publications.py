"""Artifact-publication query service and strict cursor tests."""

from __future__ import annotations

import base64
from datetime import datetime, timedelta, timezone
import json

import pytest

from encode_pipeline.platform.artifact_publications import (
    ARTIFACT_PUBLICATION_CURSOR_PREFIX,
    ArtifactPublicationFilters,
    ArtifactPublicationRunSampleBinding,
    ArtifactPublicationSummary,
    AssociatedRunSample,
    decode_artifact_publication_cursor,
    encode_artifact_publication_cursor,
)
from encode_pipeline.platform.data_registry import BindingMode, BindingProvenance
from encode_pipeline.services.artifact_publications import (
    ArtifactPublicationCursorInvalidError,
    ArtifactPublicationCursorNotFoundError,
    ArtifactPublicationDataInvalidError,
    ArtifactPublicationFilterInvalidError,
    ArtifactPublicationNotFoundError,
    ArtifactPublicationQueryService,
)


NOW = datetime(2026, 8, 7, 8, 0, tzinfo=timezone.utc)
PROJECT_ID = f"prj_{'1' * 32}"
OTHER_PROJECT_ID = f"prj_{'2' * 32}"
SAMPLE_ID = f"smp_{'3' * 32}"
SAMPLE_REVISION_ID = f"smpr_{'4' * 32}"
OTHER_SAMPLE_REVISION_ID = f"smpr_{'5' * 32}"
GENERATION_A = f"artifactgen-{'a' * 64}"
GENERATION_B = f"artifactgen-{'b' * 64}"
REVISION = f"artifactrev-{'c' * 64}"


def _binding() -> ArtifactPublicationRunSampleBinding:
    return ArtifactPublicationRunSampleBinding(
        binding_mode=BindingMode.BOUND_V1,
        provenance=BindingProvenance.RESOLVED,
        associated_run_samples=(
            AssociatedRunSample(
                sample_id=SAMPLE_ID,
                sample_revision_id=SAMPLE_REVISION_ID,
                revision_number=2,
                ordinal=0,
            ),
        ),
    )


def _summary(
    run_id: str,
    *,
    artifact_id: str = "artifact-a",
    project_id: str = PROJECT_ID,
    workflow_id: str = "workflow-a",
    output_type: str = "counts",
    generation: str = GENERATION_A,
    current_generation: str = GENERATION_A,
    published_at: datetime = NOW,
) -> ArtifactPublicationSummary:
    return ArtifactPublicationSummary(
        run_id=run_id,
        project_id=project_id,
        workflow_id=workflow_id,
        artifact_id=artifact_id,
        output_type=output_type,
        artifact_generation=generation,
        artifact_revision=REVISION,
        published_at=published_at,
        current_artifact_generation=current_generation,
        run_sample_binding=_binding(),
    )


def _summary_with_invalid_run_sample_binding() -> ArtifactPublicationSummary:
    summary = _summary("run-a")
    object.__setattr__(summary, "run_sample_binding", object())
    return summary


def _summary_with_invalid_associated_sample() -> ArtifactPublicationSummary:
    summary = _summary("run-a")
    object.__setattr__(
        summary.run_sample_binding,
        "associated_run_samples",
        (object(),),
    )
    return summary


def _matches(item: ArtifactPublicationSummary, filters: ArtifactPublicationFilters):
    return (
        (filters.project_id is None or filters.project_id == item.project_id)
        and (filters.run_id is None or filters.run_id == item.run_id)
        and (filters.workflow_id is None or filters.workflow_id == item.workflow_id)
        and (filters.output_type is None or filters.output_type == item.output_type)
        and (
            filters.associated_run_sample_revision_id is None
            or any(
                sample.sample_revision_id == filters.associated_run_sample_revision_id
                for sample in item.run_sample_binding.associated_run_samples
            )
        )
        and (
            filters.published_from is None
            or item.published_at >= filters.published_from
        )
        and (
            filters.published_before is None
            or item.published_at < filters.published_before
        )
        and (not filters.current_only or item.generation_status == "current")
    )


class _Repository:
    def __init__(self, *items: ArtifactPublicationSummary) -> None:
        self.items = list(items)

    def list_artifact_publications(self, *, filters, after, limit):
        eligible = sorted(
            (item for item in self.items if _matches(item, filters)),
            key=lambda item: item.cursor_position.key,
            reverse=True,
        )
        if after is not None:
            if not any(item.cursor_position.key == after.key for item in eligible):
                raise KeyError("boundary missing")
            eligible = [
                item for item in eligible if item.cursor_position.key < after.key
            ]
        return tuple(eligible[:limit])

    def get_artifact_publication(
        self,
        *,
        run_id,
        artifact_id,
        artifact_generation,
    ):
        for item in self.items:
            if (
                item.run_id == run_id
                and item.artifact_id == artifact_id
                and item.artifact_generation == artifact_generation
            ):
                return item
        raise KeyError("publication missing")


class _RaisingRepository(_Repository):
    def __init__(self, error: Exception) -> None:
        super().__init__()
        self.error = error

    def list_artifact_publications(self, *, filters, after, limit):
        raise self.error

    def get_artifact_publication(
        self,
        *,
        run_id,
        artifact_id,
        artifact_generation,
    ):
        raise self.error


def _rewrite_cursor(cursor: str, transform) -> str:
    body = cursor.removeprefix(ARTIFACT_PUBLICATION_CURSOR_PREFIX)
    raw = base64.urlsafe_b64decode(body + "=" * (-len(body) % 4))
    payload = json.loads(raw)
    transform(payload)
    rewritten = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    encoded = base64.urlsafe_b64encode(rewritten).decode().rstrip("=")
    return f"{ARTIFACT_PUBLICATION_CURSOR_PREFIX}{encoded}"


def test_service_pages_in_stable_order_and_binds_cursor_to_all_filters() -> None:
    repository = _Repository(
        _summary("run-a", published_at=NOW),
        _summary("run-c", published_at=NOW),
        _summary(
            "run-b",
            generation=GENERATION_A,
            current_generation=GENERATION_B,
            published_at=NOW + timedelta(minutes=1),
        ),
    )
    service = ArtifactPublicationQueryService(repository)

    first = service.list_artifact_publications(
        project_id=f" {PROJECT_ID} ",
        associated_run_sample_revision_id=SAMPLE_REVISION_ID,
        published_from=NOW - timedelta(seconds=1),
        published_before=NOW + timedelta(seconds=1),
        limit=1,
    )

    assert [item.run_id for item in first.artifact_publications] == ["run-c"]
    assert first.next_cursor is not None
    assert first.next_cursor.startswith(ARTIFACT_PUBLICATION_CURSOR_PREFIX)
    second = service.list_artifact_publications(
        project_id=PROJECT_ID,
        associated_run_sample_revision_id=SAMPLE_REVISION_ID,
        published_from=NOW - timedelta(seconds=1),
        published_before=NOW + timedelta(seconds=1),
        after=first.next_cursor,
        limit=1,
    )
    assert [item.run_id for item in second.artifact_publications] == ["run-a"]
    assert second.next_cursor is None

    with pytest.raises(ArtifactPublicationCursorInvalidError):
        service.list_artifact_publications(
            project_id=OTHER_PROJECT_ID,
            associated_run_sample_revision_id=SAMPLE_REVISION_ID,
            published_from=NOW - timedelta(seconds=1),
            published_before=NOW + timedelta(seconds=1),
            after=first.next_cursor,
        )


def test_filter_digest_normalizes_utc_and_binds_ordering() -> None:
    utc_filters = ArtifactPublicationFilters(published_from=NOW)
    offset_filters = ArtifactPublicationFilters(
        published_from=NOW.astimezone(timezone(timedelta(hours=8)))
    )
    cursor = encode_artifact_publication_cursor(
        _summary("run-a").cursor_position,
        filters=utc_filters,
    )

    assert utc_filters.sha256 == offset_filters.sha256
    assert (
        decode_artifact_publication_cursor(
            cursor,
            filters=offset_filters,
        ).position.run_id
        == "run-a"
    )


@pytest.mark.parametrize(
    "transform",
    [
        lambda payload: payload.__setitem__("extra", True),
        lambda payload: payload.__setitem__("v", True),
        lambda payload: payload["position"].__setitem__(
            "published_at", NOW.isoformat()
        ),
        lambda payload: payload["position"].__setitem__("run_id", 1),
        lambda payload: payload["position"].pop("artifact_id"),
    ],
)
def test_cursor_rejects_noncanonical_keys_types_and_timestamp(transform) -> None:
    filters = ArtifactPublicationFilters()
    cursor = encode_artifact_publication_cursor(
        _summary("run-a").cursor_position,
        filters=filters,
    )
    invalid = _rewrite_cursor(cursor, transform)

    with pytest.raises(ValueError):
        decode_artifact_publication_cursor(invalid, filters=filters)


def test_service_reports_deleted_or_filter_mismatched_boundary() -> None:
    repository = _Repository(_summary("run-a"), _summary("run-b"))
    service = ArtifactPublicationQueryService(repository)
    first = service.list_artifact_publications(limit=1)
    repository.items = [item for item in repository.items if item.run_id != "run-b"]

    with pytest.raises(ArtifactPublicationCursorNotFoundError):
        service.list_artifact_publications(after=first.next_cursor, limit=1)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"limit": 0},
        {"project_id": "/private/project"},
        {"published_from": NOW, "published_before": NOW},
        {"current_only": "true"},
    ],
)
def test_service_rejects_invalid_filters(kwargs) -> None:
    service = ArtifactPublicationQueryService(_Repository())

    with pytest.raises(ArtifactPublicationFilterInvalidError):
        service.list_artifact_publications(**kwargs)


def test_service_returns_exact_detail_and_reports_missing_identity() -> None:
    summary = _summary("run-a")
    service = ArtifactPublicationQueryService(_Repository(summary))

    assert (
        service.get_artifact_publication(
            run_id="run-a",
            artifact_id="artifact-a",
            artifact_generation=GENERATION_A,
        )
        == summary
    )
    with pytest.raises(ArtifactPublicationNotFoundError):
        service.get_artifact_publication(
            run_id="run-a",
            artifact_id="artifact-missing",
            artifact_generation=GENERATION_A,
        )


class _InvalidRepository(_Repository):
    def list_artifact_publications(self, *, filters, after, limit):
        return [_summary("run-a")]


class _StaticListRepository(_Repository):
    def __init__(self, result: tuple[object, ...]) -> None:
        super().__init__()
        self.result = result

    def list_artifact_publications(self, *, filters, after, limit):
        return self.result


@pytest.mark.parametrize(
    ("entrypoint", "repository_error"),
    [
        pytest.param(
            "list",
            KeyError("unexpected repository key"),
            id="list-key-error-without-cursor",
        ),
        pytest.param(
            "detail",
            TypeError("invalid repository row"),
            id="detail-type-error",
        ),
        pytest.param(
            "detail",
            ValueError("invalid repository row"),
            id="detail-value-error",
        ),
    ],
)
def test_service_fails_closed_for_repository_errors(
    entrypoint,
    repository_error: Exception,
) -> None:
    service = ArtifactPublicationQueryService(_RaisingRepository(repository_error))

    with pytest.raises(ArtifactPublicationDataInvalidError) as exc_info:
        if entrypoint == "list":
            service.list_artifact_publications()
        else:
            service.get_artifact_publication(
                run_id="run-a",
                artifact_id="artifact-a",
                artifact_generation=GENERATION_A,
            )

    assert exc_info.value.__cause__ is repository_error


def test_service_fails_closed_for_invalid_repository_collection() -> None:
    service = ArtifactPublicationQueryService(_InvalidRepository())

    with pytest.raises(ArtifactPublicationDataInvalidError):
        service.list_artifact_publications()


@pytest.mark.parametrize(
    "candidate_factory",
    [
        pytest.param(object, id="not-artifact-publication-summary"),
        pytest.param(
            _summary_with_invalid_run_sample_binding,
            id="invalid-run-sample-binding",
        ),
        pytest.param(
            _summary_with_invalid_associated_sample,
            id="invalid-associated-run-sample-member",
        ),
    ],
)
def test_service_fails_closed_for_invalid_repository_publications(
    candidate_factory,
) -> None:
    service = ArtifactPublicationQueryService(
        _StaticListRepository((candidate_factory(),))
    )

    with pytest.raises(ArtifactPublicationDataInvalidError) as exc_info:
        service.list_artifact_publications()

    assert isinstance(exc_info.value.__cause__, ValueError)


@pytest.mark.parametrize(
    ("kwargs", "item"),
    [
        pytest.param(
            {"project_id": OTHER_PROJECT_ID},
            _summary("run-a"),
            id="project",
        ),
        pytest.param(
            {"published_from": NOW + timedelta(seconds=1)},
            _summary("run-a"),
            id="published-from",
        ),
        pytest.param(
            {"published_before": NOW},
            _summary("run-a"),
            id="published-before",
        ),
        pytest.param(
            {},
            _summary(
                "run-a",
                generation=GENERATION_A,
                current_generation=GENERATION_B,
            ),
            id="current-only",
        ),
        pytest.param(
            {"associated_run_sample_revision_id": OTHER_SAMPLE_REVISION_ID},
            _summary("run-a"),
            id="associated-run-sample",
        ),
    ],
)
def test_service_revalidates_repository_rows_against_filters(
    kwargs,
    item: ArtifactPublicationSummary,
) -> None:
    service = ArtifactPublicationQueryService(_StaticListRepository((item,)))

    with pytest.raises(ArtifactPublicationDataInvalidError):
        service.list_artifact_publications(**kwargs)


@pytest.mark.parametrize(
    "items",
    [
        pytest.param(
            (_summary("run-a"), _summary("run-b")),
            id="not-descending",
        ),
        pytest.param(
            (
                _summary("run-c", published_at=NOW + timedelta(minutes=2)),
                _summary("run-b", published_at=NOW + timedelta(minutes=1)),
                _summary("run-a"),
            ),
            id="exceeds-limit-plus-one",
        ),
    ],
)
def test_service_fails_closed_for_repository_page_contract_violations(
    items: tuple[ArtifactPublicationSummary, ...],
) -> None:
    service = ArtifactPublicationQueryService(_StaticListRepository(items))

    with pytest.raises(ArtifactPublicationDataInvalidError):
        service.list_artifact_publications(limit=1)


def test_service_rejects_missing_detail_identity_before_repository_lookup() -> None:
    service = ArtifactPublicationQueryService(_Repository())

    with pytest.raises(ArtifactPublicationFilterInvalidError):
        service.get_artifact_publication(
            run_id=None,
            artifact_id="artifact-a",
            artifact_generation=GENERATION_A,
        )


def test_service_revalidates_exact_detail_repository_identity() -> None:
    requested = _summary("run-a")

    class _MisdirectedDetailRepository(_Repository):
        def get_artifact_publication(
            self,
            *,
            run_id,
            artifact_id,
            artifact_generation,
        ):
            return _summary("run-other")

    service = ArtifactPublicationQueryService(_MisdirectedDetailRepository(requested))

    with pytest.raises(ArtifactPublicationDataInvalidError):
        service.get_artifact_publication(
            run_id=requested.run_id,
            artifact_id=requested.artifact_id,
            artifact_generation=requested.artifact_generation,
        )
