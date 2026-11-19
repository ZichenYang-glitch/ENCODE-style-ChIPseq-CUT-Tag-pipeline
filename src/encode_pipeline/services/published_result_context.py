"""Read-only context for one published result.

Publication visibility follows the publication detail. A publication in another
project never reaches the run, reference, QC, or event stores. A visible
legacy publication is explained as unbound here; the catalog project gate is
unchanged. Nothing in this reader creates a run, starts one, or triggers
preflight.
"""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Literal

from encode_pipeline.platform.data_registry import LEGACY_PROJECT_ID
from encode_pipeline.platform.experiment_catalog import (
    catalog_not_found,
    publication_experiment_trace,
)
from encode_pipeline.platform.reference_profiles import (
    ReferenceProfileRevisionSummary,
)
from encode_pipeline.platform.result_generations import (
    QcMetricCursor,
    encode_qc_metric_cursor,
)
from encode_pipeline.platform.runs import RunQcMetric, RunStatus
from encode_pipeline.services.artifact_publications import (
    ArtifactPublicationDataInvalidError,
    ArtifactPublicationFilterInvalidError,
    ArtifactPublicationNotFoundError,
    ArtifactPublicationQueryService,
)
from encode_pipeline.services.experiment_catalog import ExperimentCatalogService
from encode_pipeline.services.reference_profiles import ReferenceProfileService
from encode_pipeline.services.run_repositories import (
    ResultGenerationChangedError,
    _publication_trace_row,
)
from encode_pipeline.services.runs import RunService


_PUBLIC_CODE = re.compile(r"^[A-Z][A-Z0-9_]{0,127}$")
_PREFLIGHT_EVENT_TYPES = frozenset({"preflight_failed", "dry_run_failed"})
_QC_PAGE_LIMIT = 100
_EVENT_PAGE_LIMIT = 100
_EVENT_PAGE_CAP = 50
_MISSING_RESULT_STATE = "run result generation state is missing"


class PublishedResultEvidenceError(Exception):
    """Persisted result evidence could not be projected safely."""

    def __init__(self, code: str, message: str) -> None:
        super().__init__(message)
        self.code = code
        self.message = message


@dataclass(frozen=True)
class RecordedPreflightFailure:
    """Public projection of one already recorded preflight failure."""

    event_id: str
    event_type: str
    message: str
    stage: str | None = None
    code: str | None = None


@dataclass(frozen=True)
class PublishedResultContext:
    """Same-identity read of one published generation."""

    trace: dict
    run_status: str
    current_stage: str | None
    reference_evidence: Literal["absent", "present"]
    reference: ReferenceProfileRevisionSummary | None
    qc_evidence: Literal["absent", "present", "unverifiable"]
    qc_generation: str | None
    qc_metrics: tuple[RunQcMetric, ...]
    qc_next_cursor: str | None
    preflight_failure_evidence: Literal["absent", "present"]
    preflight_failures: tuple[RecordedPreflightFailure, ...]


class PublishedResultContextReader:
    """Compose catalog, run, reference, QC, and event reads."""

    def __init__(
        self,
        catalog: ExperimentCatalogService,
        runs: RunService,
        reference_profiles: ReferenceProfileService,
        publications: ArtifactPublicationQueryService,
    ) -> None:
        self._catalog = catalog
        self._runs = runs
        self._reference_profiles = reference_profiles
        self._publications = publications

    def read(
        self,
        *,
        project_id: str,
        run_id: str,
        artifact_id: str,
        artifact_generation: str,
    ) -> PublishedResultContext:
        publication = self._visible_publication(
            project_id=project_id,
            run_id=run_id,
            artifact_id=artifact_id,
            artifact_generation=artifact_generation,
        )
        if publication.project_id == LEGACY_PROJECT_ID:
            trace = publication_experiment_trace(
                project_id=project_id,
                run_id=run_id,
                artifact_id=artifact_id,
                artifact_generation=artifact_generation,
                result=_publication_trace_row(publication),
                binding_revision_id=None,
                binding_payload_digest=None,
                revision_project_id=None,
                revision_payload_digest=None,
            )
        else:
            trace = self._catalog.get_publication_experiment_trace(
                project_id=project_id,
                run_id=run_id,
                artifact_id=artifact_id,
                artifact_generation=artifact_generation,
            )
        try:
            record = self._runs.get_run(run_id)
        except KeyError as exc:
            raise PublishedResultEvidenceError(
                "RUN_RESULT_CONTEXT_INCONSISTENT",
                "Published result context could not be read.",
            ) from exc
        status = record.status
        run_status = status.value if isinstance(status, RunStatus) else str(status)
        reference_evidence, reference = self._reference(run_id)
        qc_evidence, generation, metrics, cursor = self._qc(
            run_id,
            artifact_generation,
        )
        failures = self._preflight_failures(run_id)
        return PublishedResultContext(
            trace=trace,
            run_status=run_status,
            current_stage=record.current_stage,
            reference_evidence=reference_evidence,
            reference=reference,
            qc_evidence=qc_evidence,
            qc_generation=generation,
            qc_metrics=metrics,
            qc_next_cursor=cursor,
            preflight_failure_evidence=("present" if failures else "absent"),
            preflight_failures=failures,
        )

    def _reference(
        self,
        run_id: str,
    ) -> tuple[Literal["absent", "present"], ReferenceProfileRevisionSummary | None]:
        binding = self._runs.get_run_reference_binding(run_id)
        if binding is None:
            return "absent", None
        try:
            summary = self._reference_profiles.get_revision_summary(binding.revision_id)
        except (KeyError, RuntimeError, ValueError) as exc:
            raise PublishedResultEvidenceError(
                "RUN_REFERENCE_EVIDENCE_INVALID",
                "Run reference evidence could not be read safely.",
            ) from exc
        if (
            summary.profile_id != binding.profile_id
            or summary.public_identity_sha256 != binding.revision_public_identity_sha256
        ):
            raise PublishedResultEvidenceError(
                "RUN_REFERENCE_EVIDENCE_INVALID",
                "Run reference evidence could not be read safely.",
            )
        return "present", summary

    def _visible_publication(
        self,
        *,
        project_id: str,
        run_id: str,
        artifact_id: str,
        artifact_generation: str,
    ):
        try:
            publication = self._publications.get_artifact_publication(
                run_id=run_id,
                artifact_id=artifact_id,
                artifact_generation=artifact_generation,
            )
        except ArtifactPublicationNotFoundError as exc:
            raise catalog_not_found() from exc
        except ArtifactPublicationFilterInvalidError as exc:
            raise catalog_not_found() from exc
        except ArtifactPublicationDataInvalidError as exc:
            raise PublishedResultEvidenceError(
                "RUN_RESULT_CONTEXT_INCONSISTENT",
                "Published result context could not be read.",
            ) from exc
        if publication.project_id != project_id:
            raise catalog_not_found()
        return publication

    def _qc(
        self,
        run_id: str,
        artifact_generation: str,
    ) -> tuple[
        Literal["absent", "present", "unverifiable"],
        str | None,
        tuple[RunQcMetric, ...],
        str | None,
    ]:
        try:
            state = self._runs.get_result_state(run_id)
        except ValueError as exc:
            if str(exc) == _MISSING_RESULT_STATE:
                return "absent", None, (), None
            raise PublishedResultEvidenceError(
                "RUN_QC_METRIC_DATA_INVALID",
                "Persisted QC metric data is unavailable.",
            ) from exc
        except KeyError as exc:
            raise PublishedResultEvidenceError(
                "RUN_QC_METRIC_DATA_INVALID",
                "Persisted QC metric data is unavailable.",
            ) from exc
        if state.artifact_generation != artifact_generation:
            return "unverifiable", None, (), None
        try:
            generation, metrics = self._runs.list_qc_metrics_page(
                run_id,
                expected_generation=state.qc_generation,
                limit=_QC_PAGE_LIMIT + 1,
            )
        except ResultGenerationChangedError:
            return "unverifiable", None, (), None
        except ValueError as exc:
            if str(exc) == _MISSING_RESULT_STATE:
                return "absent", None, (), None
            raise PublishedResultEvidenceError(
                "RUN_QC_METRIC_DATA_INVALID",
                "Persisted QC metric data is unavailable.",
            ) from exc
        except (AttributeError, KeyError, TypeError) as exc:
            raise PublishedResultEvidenceError(
                "RUN_QC_METRIC_DATA_INVALID",
                "Persisted QC metric data is unavailable.",
            ) from exc
        page = tuple(metrics)
        cursor = None
        if len(page) > _QC_PAGE_LIMIT:
            if generation is None:
                raise PublishedResultEvidenceError(
                    "RUN_QC_METRIC_DATA_INVALID",
                    "Persisted QC metric data is unavailable.",
                )
            page = page[:_QC_PAGE_LIMIT]
            cursor = encode_qc_metric_cursor(
                QcMetricCursor(
                    run_id=run_id,
                    qc_generation=generation,
                    after_metric_id=page[-1].metric_id,
                )
            )
        summary_bound = (
            generation is not None
            and generation == state.qc_generation
            and state.qc_artifact_generation == artifact_generation
        )
        if not summary_bound:
            if page:
                raise PublishedResultEvidenceError(
                    "RUN_QC_METRIC_DATA_INVALID",
                    "Persisted QC metric data is unavailable.",
                )
            return "absent", None, (), None
        return "present", generation, page, cursor

    def _preflight_failures(
        self,
        run_id: str,
    ) -> tuple[RecordedPreflightFailure, ...]:
        failures: list[RecordedPreflightFailure] = []
        after: str | None = None
        for _page in range(_EVENT_PAGE_CAP):
            try:
                events = self._runs.list_events(
                    run_id,
                    after=after,
                    limit=_EVENT_PAGE_LIMIT,
                )
            except KeyError as exc:
                raise PublishedResultEvidenceError(
                    "RUN_EVENT_READ_INCOMPLETE",
                    "Recorded preflight failures could not be read completely.",
                ) from exc
            for event in events:
                failure = _preflight_failure(event)
                if failure is not None:
                    failures.append(failure)
            if len(events) < _EVENT_PAGE_LIMIT:
                return tuple(failures)
            after = events[-1].event_id
        raise PublishedResultEvidenceError(
            "RUN_EVENT_READ_INCOMPLETE",
            "Recorded preflight failures could not be read completely.",
        )


def _preflight_failure(event: object) -> RecordedPreflightFailure | None:
    event_type = getattr(event, "event_type", None)
    issue = getattr(event, "issue", None)
    code = getattr(issue, "code", None) if issue is not None else None
    failed = event_type in _PREFLIGHT_EVENT_TYPES or code == "PREFLIGHT_FAILED"
    if not failed:
        return None
    public_code = (
        code if isinstance(code, str) and _PUBLIC_CODE.fullmatch(code) else None
    )
    stage = getattr(event, "stage", None)
    if not isinstance(stage, str) or "/" in stage or "\\" in stage or len(stage) > 64:
        stage = None
    return RecordedPreflightFailure(
        event_id=str(getattr(event, "event_id")),
        event_type=str(event_type),
        stage=stage,
        message=_public_message(getattr(event, "message", None)),
        code=public_code,
    )


def _public_message(message: object) -> str:
    if not isinstance(message, str):
        return "A preflight failure was recorded."
    text = message.strip()
    if (
        not text
        or len(text) > 240
        or "/" in text
        or "\\" in text
        or ".." in text
        or any(ord(character) < 32 for character in text)
    ):
        return "A preflight failure was recorded."
    return text
