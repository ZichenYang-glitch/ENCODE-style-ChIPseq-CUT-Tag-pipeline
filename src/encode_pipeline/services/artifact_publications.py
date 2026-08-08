"""Disclosure-safe artifact-publication query orchestration."""

from __future__ import annotations

from datetime import datetime

from encode_pipeline.platform.artifact_publications import (
    ArtifactPublicationFilters,
    ArtifactPublicationPage,
    ArtifactPublicationRepository,
    ArtifactPublicationRunSampleBinding,
    ArtifactPublicationSummary,
    AssociatedRunSample,
    decode_artifact_publication_cursor,
    encode_artifact_publication_cursor,
    normalize_artifact_publication_artifact_id,
    normalize_artifact_publication_run_id,
)
from encode_pipeline.platform.result_generations import validate_artifact_generation


class ArtifactPublicationFilterInvalidError(ValueError):
    """A publication request filter is malformed or unsupported."""


class ArtifactPublicationCursorInvalidError(ValueError):
    """A publication cursor is malformed or belongs to another query."""


class ArtifactPublicationCursorNotFoundError(KeyError):
    """A valid publication cursor boundary is no longer queryable."""


class ArtifactPublicationNotFoundError(KeyError):
    """An exact public artifact publication identity does not exist."""


class ArtifactPublicationDataInvalidError(ValueError):
    """Persisted publication data violates the public projection contract."""


class ArtifactPublicationQueryService:
    """Read-only service for list and exact publication queries."""

    def __init__(self, repository: ArtifactPublicationRepository) -> None:
        self._repository = repository

    def list_artifact_publications(
        self,
        *,
        project_id: str | None = None,
        run_id: str | None = None,
        workflow_id: str | None = None,
        output_type: str | None = None,
        associated_run_sample_revision_id: str | None = None,
        published_from: datetime | None = None,
        published_before: datetime | None = None,
        current_only: bool = True,
        after: str | None = None,
        limit: int = 50,
    ) -> ArtifactPublicationPage:
        """Return one filter-bound descending keyset page."""
        if (
            isinstance(limit, bool)
            or not isinstance(limit, int)
            or not 1 <= limit <= 100
        ):
            raise ArtifactPublicationFilterInvalidError(
                "artifact publication limit is invalid"
            )
        try:
            filters = ArtifactPublicationFilters(
                project_id=project_id,
                run_id=run_id,
                workflow_id=workflow_id,
                output_type=output_type,
                associated_run_sample_revision_id=(associated_run_sample_revision_id),
                published_from=published_from,
                published_before=published_before,
                current_only=current_only,
            )
        except (TypeError, ValueError) as exc:
            raise ArtifactPublicationFilterInvalidError(
                "artifact publication filters are invalid"
            ) from exc

        position = None
        if after is not None:
            try:
                position = decode_artifact_publication_cursor(
                    after,
                    filters=filters,
                ).position
            except (TypeError, ValueError) as exc:
                raise ArtifactPublicationCursorInvalidError(
                    "artifact publication cursor is invalid"
                ) from exc

        try:
            candidates = self._repository.list_artifact_publications(
                filters=filters,
                after=position,
                limit=limit + 1,
            )
            validated = _validated_publications(candidates, filters=filters)
        except KeyError as exc:
            if position is None:
                raise ArtifactPublicationDataInvalidError(
                    "persisted artifact publication data is invalid"
                ) from exc
            raise ArtifactPublicationCursorNotFoundError(
                "artifact publication cursor boundary does not exist"
            ) from exc
        except (TypeError, ValueError) as exc:
            raise ArtifactPublicationDataInvalidError(
                "persisted artifact publication data is invalid"
            ) from exc

        if len(validated) > limit + 1:
            raise ArtifactPublicationDataInvalidError(
                "publication repository exceeded the requested bound"
            )
        page_items = validated[:limit]
        next_cursor = None
        if len(validated) > limit:
            next_cursor = encode_artifact_publication_cursor(
                page_items[-1].cursor_position,
                filters=filters,
            )
        return ArtifactPublicationPage(
            artifact_publications=page_items,
            next_cursor=next_cursor,
        )

    def get_artifact_publication(
        self,
        *,
        run_id: str,
        artifact_id: str,
        artifact_generation: str,
    ) -> ArtifactPublicationSummary:
        """Return one exact public compound identity."""
        try:
            normalized_run_id = normalize_artifact_publication_run_id(run_id)
            normalized_artifact_id = normalize_artifact_publication_artifact_id(
                artifact_id
            )
            normalized_generation = validate_artifact_generation(artifact_generation)
            if normalized_run_id is None or normalized_artifact_id is None:
                raise ValueError("publication identity is required")
        except (TypeError, ValueError) as exc:
            raise ArtifactPublicationFilterInvalidError(
                "artifact publication identity is invalid"
            ) from exc

        try:
            candidate = self._repository.get_artifact_publication(
                run_id=normalized_run_id,
                artifact_id=normalized_artifact_id,
                artifact_generation=normalized_generation,
            )
        except KeyError as exc:
            raise ArtifactPublicationNotFoundError(
                "artifact publication does not exist"
            ) from exc
        except (TypeError, ValueError) as exc:
            raise ArtifactPublicationDataInvalidError(
                "persisted artifact publication data is invalid"
            ) from exc

        try:
            validated = _revalidate_publication(candidate)
            if (
                validated.run_id != normalized_run_id
                or validated.artifact_id != normalized_artifact_id
                or validated.artifact_generation != normalized_generation
            ):
                raise ValueError("publication repository returned another identity")
            return validated
        except (TypeError, ValueError) as exc:
            raise ArtifactPublicationDataInvalidError(
                "persisted artifact publication data is invalid"
            ) from exc


# The shorter name is retained for app composition and dependency injection.
ArtifactPublicationService = ArtifactPublicationQueryService


def _validated_publications(
    candidates: object,
    *,
    filters: ArtifactPublicationFilters,
) -> tuple[ArtifactPublicationSummary, ...]:
    if not isinstance(candidates, tuple):
        raise ValueError("publication repository must return a tuple")
    validated = tuple(_revalidate_publication(item) for item in candidates)
    previous_key = None
    for item in validated:
        key = item.cursor_position.key
        if previous_key is not None and previous_key <= key:
            raise ValueError("publication repository order is invalid")
        previous_key = key
        _validate_publication_matches_filters(item, filters)
    return validated


def _revalidate_publication(candidate: object) -> ArtifactPublicationSummary:
    if not isinstance(candidate, ArtifactPublicationSummary):
        raise ValueError("publication repository returned an invalid value")
    binding = candidate.run_sample_binding
    if not isinstance(binding, ArtifactPublicationRunSampleBinding):
        raise ValueError("publication Run sample binding is invalid")
    if any(
        not isinstance(item, AssociatedRunSample)
        for item in binding.associated_run_samples
    ):
        raise ValueError("publication associated Run samples are invalid")
    samples = tuple(
        AssociatedRunSample(
            sample_id=item.sample_id,
            sample_revision_id=item.sample_revision_id,
            revision_number=item.revision_number,
            ordinal=item.ordinal,
        )
        for item in binding.associated_run_samples
    )
    return ArtifactPublicationSummary(
        run_id=candidate.run_id,
        project_id=candidate.project_id,
        workflow_id=candidate.workflow_id,
        artifact_id=candidate.artifact_id,
        output_type=candidate.output_type,
        artifact_generation=candidate.artifact_generation,
        artifact_revision=candidate.artifact_revision,
        published_at=candidate.published_at,
        current_artifact_generation=candidate.current_artifact_generation,
        run_sample_binding=ArtifactPublicationRunSampleBinding(
            binding_mode=binding.binding_mode,
            provenance=binding.provenance,
            associated_run_samples=samples,
        ),
    )


def _validate_publication_matches_filters(
    item: ArtifactPublicationSummary,
    filters: ArtifactPublicationFilters,
) -> None:
    exact_matches = (
        (filters.project_id, item.project_id),
        (filters.run_id, item.run_id),
        (filters.workflow_id, item.workflow_id),
        (filters.output_type, item.output_type),
    )
    if any(
        expected is not None and expected != actual
        for expected, actual in exact_matches
    ):
        raise ValueError("publication repository row does not match filters")
    if (
        filters.published_from is not None
        and item.published_at < filters.published_from
    ):
        raise ValueError("publication repository row precedes published_from")
    if (
        filters.published_before is not None
        and item.published_at >= filters.published_before
    ):
        raise ValueError("publication repository row reaches published_before")
    if filters.current_only and item.generation_status != "current":
        raise ValueError("publication repository returned a superseded generation")
    if filters.associated_run_sample_revision_id is not None and all(
        sample.sample_revision_id != filters.associated_run_sample_revision_id
        for sample in item.run_sample_binding.associated_run_samples
    ):
        raise ValueError("publication repository row lacks associated Run sample")
