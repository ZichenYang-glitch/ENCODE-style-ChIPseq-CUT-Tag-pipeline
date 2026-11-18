"""In-memory experiment catalog with the same write and query contract as SQL."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
import json
from threading import RLock
from typing import Any, Protocol

from encode_pipeline.platform.data_registry import (
    LEGACY_PROJECT_ID,
    SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
    ProjectKind,
    SampleRevision,
    build_sample_revision_payload_digest,
    canonical_sample_revision_payload,
)
from encode_pipeline.platform.experiment_catalog import (
    EXPERIMENT_FIELDS,
    SAMPLE_FIELDS,
    TERM_TYPES,
    CatalogCursorPosition,
    CatalogFilters,
    CatalogWrite,
    ExperimentPayload,
    ExperimentSamplePin,
    SampleCatalogBlock,
    TermLookup,
    TermPin,
    administrator_required,
    bind_catalog_target,
    block_matches,
    catalog_already_archived,
    catalog_cursor_not_found,
    catalog_experiment_archived,
    catalog_idempotency_conflict,
    ExperimentRevisionFreezeSource,
    build_experiment_revision_payload_digest,
    catalog_not_found,
    catalog_payload_invalid,
    catalog_revision_conflict,
    catalog_term_archived,
    collect_experiment_term_revision_ids,
    collect_term_revision_ids,
    decode_catalog_cursor,
    encode_catalog_cursor,
    interpret_experiment_fields,
    interpret_sample_catalog,
    interpret_stored_experiment_payload,
    render_sample_material,
    render_single,
    sample_catalog_object,
    slot_matches,
    validate_aliases,
    validate_display_name,
    validate_optional_external,
    validate_stable_key,
)
from encode_pipeline.services.data_registry_repositories import (
    DataRegistryConflictError,
    InMemoryDataRegistryRepository,
)


class ExperimentCatalogRepository(Protocol):
    """Catalog operations shared by the SQL and in-memory stores."""

    def create_term(self, write: CatalogWrite, **kwargs: Any) -> dict[str, Any]: ...

    def revise_term(self, write: CatalogWrite, **kwargs: Any) -> dict[str, Any]: ...

    def archive_term(self, write: CatalogWrite, **kwargs: Any) -> dict[str, Any]: ...

    def create_experiment(
        self, write: CatalogWrite, **kwargs: Any
    ) -> dict[str, Any]: ...

    def revise_experiment(
        self, write: CatalogWrite, **kwargs: Any
    ) -> dict[str, Any]: ...

    def archive_experiment(
        self, write: CatalogWrite, **kwargs: Any
    ) -> dict[str, Any]: ...

    def revise_sample_catalog(
        self, write: CatalogWrite, **kwargs: Any
    ) -> dict[str, Any]: ...

    def list_experiments(self, **kwargs: Any) -> dict[str, Any]: ...

    def list_samples(self, **kwargs: Any) -> dict[str, Any]: ...

    def list_terms(self, **kwargs: Any) -> dict[str, Any]: ...

    def get_term(self, term_id: str) -> dict[str, Any]: ...

    def list_sample_revisions(self, **kwargs: Any) -> dict[str, Any]: ...

    def get_sample_revision(
        self, sample_id: str, sample_revision_id: str
    ) -> dict[str, Any]: ...

    def list_experiment_revisions(self, **kwargs: Any) -> dict[str, Any]: ...

    def get_experiment_revision(
        self, experiment_id: str, experiment_revision_id: str
    ) -> dict[str, Any]: ...

    def list_sample_run_results(self, **kwargs: Any) -> dict[str, Any]: ...

    def list_experiment_run_results(self, **kwargs: Any) -> dict[str, Any]: ...


@dataclass
class _TermRevision:
    term_revision_id: str
    term_id: str
    revision_number: int
    display_name: str
    namespace: str | None
    accession: str | None
    aliases: tuple[str, ...]
    change_kind: str
    created_at: datetime


@dataclass
class _Term:
    term_id: str
    term_type: str
    created_at: datetime
    archived_at: datetime | None
    revisions: list[_TermRevision]


@dataclass
class _ExperimentRevision:
    experiment_revision_id: str
    experiment_id: str
    project_id: str
    revision_number: int
    canonical_payload: str
    payload: ExperimentPayload
    created_at: datetime


@dataclass
class _Experiment:
    experiment_id: str
    project_id: str
    stable_key: str
    created_at: datetime
    archived_at: datetime | None
    revisions: list[_ExperimentRevision]


@dataclass
class _Receipt:
    operation: str
    request_sha256: str
    target_key: str
    subject_id: str
    result_revision_id: str | None
    result_archived_at: datetime | None


class _MemoryTermLookup(TermLookup):
    def __init__(self, repository: "InMemoryExperimentCatalogRepository") -> None:
        self._repository = repository

    def __call__(self, term_revision_id: str) -> TermPin | None:
        revision = self._repository._term_revisions.get(term_revision_id)
        if revision is None:
            return None
        term = self._repository._terms[revision.term_id]
        return TermPin(
            term_id=term.term_id,
            term_revision_id=revision.term_revision_id,
            term_type=term.term_type,
            display_name=revision.display_name,
            archived=term.archived_at is not None,
        )


class InMemoryExperimentCatalogRepository:
    """Lock-protected catalog store used to parity-check the SQL repository."""

    def __init__(
        self,
        registry: InMemoryDataRegistryRepository,
        *,
        is_administrator: Callable[[str], bool] | None = None,
    ) -> None:
        if not isinstance(registry, InMemoryDataRegistryRepository):
            raise ValueError("registry must be an InMemoryDataRegistryRepository")
        self._registry = registry
        self._is_administrator = is_administrator or (lambda _user_id: True)
        self._terms: dict[str, _Term] = {}
        self._term_revisions: dict[str, _TermRevision] = {}
        self._experiments: dict[str, _Experiment] = {}
        self._experiment_revisions: dict[str, _ExperimentRevision] = {}
        self._receipts: dict[tuple[str, str], _Receipt] = {}
        self._lock = RLock()

    def create_term(
        self,
        write: CatalogWrite,
        *,
        term_id: str,
        term_revision_id: str,
        term_type: str,
        display_name: str,
        namespace: str | None,
        accession: str | None,
        aliases: tuple[str, ...],
    ) -> dict[str, Any]:
        def apply() -> dict[str, Any]:
            if term_type not in TERM_TYPES:
                raise catalog_payload_invalid()
            if term_id in self._terms or term_revision_id in self._term_revisions:
                raise catalog_revision_conflict()
            name = validate_display_name(display_name)
            revision = _TermRevision(
                term_revision_id=term_revision_id,
                term_id=term_id,
                revision_number=1,
                display_name=name,
                namespace=validate_optional_external(namespace),
                accession=validate_optional_external(accession),
                aliases=tuple(json.loads(self._alias_json(aliases))),
                change_kind="initial",
                created_at=write.now,
            )
            self._terms[term_id] = _Term(
                term_id=term_id,
                term_type=term_type,
                created_at=write.now,
                archived_at=None,
                revisions=[revision],
            )
            self._term_revisions[term_revision_id] = revision
            self._add_receipt(
                write,
                subject_id=term_id,
                result_revision_id=term_revision_id,
                result_archived_at=None,
            )
            return self._term_body(term_id, term_revision_id)

        write = bind_catalog_target(write, "vocabulary")
        return self._write(write, apply)

    def revise_term(
        self,
        write: CatalogWrite,
        *,
        term_id: str,
        term_revision_id: str,
        expected_revision_number: int,
        display_name: str,
        namespace: str | None,
        accession: str | None,
        aliases: tuple[str, ...],
    ) -> dict[str, Any]:
        def apply() -> dict[str, Any]:
            term = self._terms.get(term_id)
            if term is None:
                raise catalog_not_found()
            current = term.revisions[-1].revision_number
            _expect_revision(current, expected_revision_number)
            if term.archived_at is not None:
                raise catalog_term_archived()
            if term_revision_id in self._term_revisions:
                raise catalog_revision_conflict()
            revision = _TermRevision(
                term_revision_id=term_revision_id,
                term_id=term_id,
                revision_number=current + 1,
                display_name=validate_display_name(display_name),
                namespace=validate_optional_external(namespace),
                accession=validate_optional_external(accession),
                aliases=tuple(json.loads(self._alias_json(aliases))),
                change_kind="editorial",
                created_at=write.now,
            )
            term.revisions.append(revision)
            self._term_revisions[term_revision_id] = revision
            self._add_receipt(
                write,
                subject_id=term_id,
                result_revision_id=term_revision_id,
                result_archived_at=None,
            )
            return self._term_body(term_id, term_revision_id)

        write = bind_catalog_target(write, f"term:{term_id}")
        return self._write(write, apply)

    def archive_term(
        self,
        write: CatalogWrite,
        *,
        term_id: str,
        expected_revision_number: int,
    ) -> dict[str, Any]:
        def apply() -> dict[str, Any]:
            term = self._terms.get(term_id)
            if term is None:
                raise catalog_not_found()
            current = term.revisions[-1].revision_number
            _expect_revision(current, expected_revision_number)
            if term.archived_at is not None:
                raise catalog_already_archived()
            term.archived_at = write.now
            self._add_receipt(
                write,
                subject_id=term_id,
                result_revision_id=None,
                result_archived_at=write.now,
            )
            return self._term_body(term_id, term.revisions[-1].term_revision_id)

        write = bind_catalog_target(write, f"term:{term_id}")
        return self._write(write, apply)

    def create_experiment(
        self,
        write: CatalogWrite,
        *,
        project_id: str,
        experiment_id: str,
        experiment_revision_id: str,
        stable_key: str,
        display_name: str,
        assay: object,
        target: object,
        batch: object,
        samples: object,
    ) -> dict[str, Any]:
        def apply() -> dict[str, Any]:
            self._project(project_id)
            key = validate_stable_key(stable_key)
            if any(
                item.project_id == project_id and item.stable_key == key
                for item in self._experiments.values()
            ):
                raise catalog_revision_conflict()
            if (
                experiment_id in self._experiments
                or experiment_revision_id in self._experiment_revisions
            ):
                raise catalog_revision_conflict()
            payload = self._experiment_payload(
                project_id=project_id,
                display_name=validate_display_name(display_name),
                assay=assay,
                target=target,
                batch=batch,
                samples=samples,
            )
            canonical = payload.canonical_json()
            revision = _ExperimentRevision(
                experiment_revision_id=experiment_revision_id,
                experiment_id=experiment_id,
                project_id=project_id,
                revision_number=1,
                canonical_payload=canonical,
                payload=payload,
                created_at=write.now,
            )
            self._experiments[experiment_id] = _Experiment(
                experiment_id=experiment_id,
                project_id=project_id,
                stable_key=key,
                created_at=write.now,
                archived_at=None,
                revisions=[revision],
            )
            self._experiment_revisions[experiment_revision_id] = revision
            self._add_receipt(
                write,
                subject_id=experiment_id,
                result_revision_id=experiment_revision_id,
                result_archived_at=None,
            )
            return self._experiment_detail(
                experiment_revision_id,
                include_candidates=True,
            )

        write = bind_catalog_target(write, f"project:{project_id}")
        return self._write(write, apply)

    def revise_experiment(
        self,
        write: CatalogWrite,
        *,
        experiment_id: str,
        experiment_revision_id: str,
        expected_revision_number: int,
        display_name: str,
        assay: object,
        target: object,
        batch: object,
        samples: object,
    ) -> dict[str, Any]:
        def apply() -> dict[str, Any]:
            experiment = self._experiments.get(experiment_id)
            if experiment is None:
                raise catalog_not_found()
            self._project(experiment.project_id)
            current = experiment.revisions[-1].revision_number
            _expect_revision(current, expected_revision_number)
            if experiment.archived_at is not None:
                raise catalog_experiment_archived()
            if experiment_revision_id in self._experiment_revisions:
                raise catalog_revision_conflict()
            payload = self._experiment_payload(
                project_id=experiment.project_id,
                display_name=validate_display_name(display_name),
                assay=assay,
                target=target,
                batch=batch,
                samples=samples,
            )
            revision = _ExperimentRevision(
                experiment_revision_id=experiment_revision_id,
                experiment_id=experiment_id,
                project_id=experiment.project_id,
                revision_number=current + 1,
                canonical_payload=payload.canonical_json(),
                payload=payload,
                created_at=write.now,
            )
            experiment.revisions.append(revision)
            self._experiment_revisions[experiment_revision_id] = revision
            self._add_receipt(
                write,
                subject_id=experiment_id,
                result_revision_id=experiment_revision_id,
                result_archived_at=None,
            )
            return self._experiment_detail(
                experiment_revision_id,
                include_candidates=True,
            )

        write = bind_catalog_target(write, f"experiment:{experiment_id}")
        return self._write(write, apply)

    def archive_experiment(
        self,
        write: CatalogWrite,
        *,
        experiment_id: str,
        expected_revision_number: int,
    ) -> dict[str, Any]:
        def apply() -> dict[str, Any]:
            experiment = self._experiments.get(experiment_id)
            if experiment is None:
                raise catalog_not_found()
            self._project(experiment.project_id)
            current = experiment.revisions[-1].revision_number
            _expect_revision(current, expected_revision_number)
            if experiment.archived_at is not None:
                raise catalog_already_archived()
            experiment.archived_at = write.now
            self._add_receipt(
                write,
                subject_id=experiment_id,
                result_revision_id=None,
                result_archived_at=write.now,
            )
            return {
                "experiment_id": experiment.experiment_id,
                "archived_at": _stamp(write.now),
                "revision_number": current,
            }

        write = bind_catalog_target(write, f"experiment:{experiment_id}")
        return self._write(write, apply)

    def revise_sample_catalog(
        self,
        write: CatalogWrite,
        *,
        sample_id: str,
        sample_revision_id: str,
        expected_revision_number: int,
        display_name: str | None,
        catalog: object,
    ) -> dict[str, Any]:
        def apply() -> dict[str, Any]:
            try:
                sample = self._registry.get_sample(sample_id)
            except KeyError:
                raise catalog_not_found() from None
            self._project(sample.project_id)
            try:
                head_id, head_number = self._registry.get_sample_revision_head(
                    sample_id
                )
            except KeyError:
                raise catalog_not_found() from None
            _expect_revision(head_number, expected_revision_number)
            previous = self._registry.get_sample_revision(head_id)
            decoded = json.loads(previous.canonical_payload)
            attributes = {
                key: value
                for key, value in decoded["attributes"].items()
                if key != "catalog"
            }
            block = interpret_sample_catalog(
                catalog,
                _MemoryTermLookup(self),
                for_write=True,
            )
            if block is not None:
                attributes["catalog"] = sample_catalog_object(block)
            name = (
                validate_display_name(display_name)
                if display_name is not None
                else decoded["display_name"]
            )
            canonical = canonical_sample_revision_payload(
                display_name=name,
                attributes=attributes,
            )
            revision = SampleRevision(
                sample_revision_id=sample_revision_id,
                project_id=sample.project_id,
                sample_id=sample.sample_id,
                revision_number=previous.revision_number + 1,
                canonical_payload=canonical,
                payload_digest_scheme=SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
                payload_digest=build_sample_revision_payload_digest(canonical),
                created_at=write.now,
            )
            try:
                self._registry.append_sample_revision(
                    revision,
                    expected_previous_revision_number=previous.revision_number,
                )
            except DataRegistryConflictError:
                raise catalog_revision_conflict() from None
            self._registry.mark_catalog_revision_trusted(sample_revision_id)
            self._add_receipt(
                write,
                subject_id=sample_id,
                result_revision_id=sample_revision_id,
                result_archived_at=None,
            )
            return self._sample_revision_body(
                sample_revision_id,
                include_candidates=True,
            )

        write = bind_catalog_target(write, f"sample:{sample_id}")
        return self._write(write, apply)

    def list_experiments(
        self,
        *,
        filters: CatalogFilters,
        cursor: str | None,
        limit: int,
    ) -> dict[str, Any]:
        with self._lock:
            self._project(_required_project(filters))
            matched = [
                item
                for item in self._experiments.values()
                if item.project_id == filters.project_id
                and self._experiment_matches(item, filters)
            ]
            page, next_cursor = _page(
                matched,
                kind="experiments",
                filters=filters,
                cursor=cursor,
                limit=limit,
                created_at=lambda item: item.created_at,
                object_id=lambda item: item.experiment_id,
            )
            return {
                "ok": True,
                "project_id": filters.project_id,
                "revision_view": filters.revision_view,
                "experiments": [
                    self._experiment_detail(
                        item.revisions[-1].experiment_revision_id,
                        include_candidates=False,
                        revision_view=filters.revision_view,
                    )
                    for item in page
                ],
                "next_cursor": next_cursor,
            }

    def list_samples(
        self,
        *,
        filters: CatalogFilters,
        cursor: str | None,
        limit: int,
    ) -> dict[str, Any]:
        with self._lock:
            project_id = _required_project(filters)
            self._project(project_id)
            matched = [
                sample
                for sample in self._registry.list_samples(project_id)
                if self._sample_matches(sample, filters)
            ]
            page, next_cursor = _page(
                matched,
                kind="samples",
                filters=filters,
                cursor=cursor,
                limit=limit,
                created_at=lambda item: item.created_at,
                object_id=lambda item: item.sample_id,
            )
            return {
                "project_id": project_id,
                "samples": [self._sample_summary(sample.sample_id) for sample in page],
                "next_cursor": next_cursor,
            }

    def list_terms(
        self,
        *,
        filters: CatalogFilters,
        cursor: str | None,
        limit: int,
    ) -> dict[str, Any]:
        with self._lock:
            matched = [
                term
                for term in self._terms.values()
                if self._term_matches(term, filters)
            ]
            page, next_cursor = _page(
                matched,
                kind="terms",
                filters=filters,
                cursor=cursor,
                limit=limit,
                created_at=lambda item: item.created_at,
                object_id=lambda item: item.term_id,
            )
            terms = []
            for term in page:
                head = term.revisions[-1]
                terms.append(
                    {
                        "term_id": term.term_id,
                        "term_type": term.term_type,
                        "term_revision_id": head.term_revision_id,
                        "revision_number": head.revision_number,
                        "display_name": head.display_name,
                        "archived_at": _stamp(term.archived_at),
                    }
                )
            return {"terms": terms, "next_cursor": next_cursor}

    def get_term(self, term_id: str) -> dict[str, Any]:
        with self._lock:
            term = self._terms.get(term_id)
            if term is None:
                raise catalog_not_found()
            return self._term_body(term_id, term.revisions[-1].term_revision_id)

    def list_sample_revisions(
        self,
        *,
        sample_id: str,
        filters: CatalogFilters,
        cursor: str | None,
        limit: int,
    ) -> dict[str, Any]:
        with self._lock:
            self._project(self._sample_project(sample_id))
            revisions = self._registry.list_sample_revisions(sample_id)
            page, next_cursor = _page(
                list(revisions),
                kind="sample_revisions",
                filters=filters,
                cursor=cursor,
                limit=limit,
                created_at=lambda item: item.created_at,
                object_id=lambda item: item.sample_revision_id,
            )
            return {
                "sample_revisions": [
                    self._sample_revision_body(
                        item.sample_revision_id,
                        include_candidates=False,
                    )
                    for item in page
                ],
                "next_cursor": next_cursor,
            }

    def get_sample_revision(
        self, sample_id: str, sample_revision_id: str
    ) -> dict[str, Any]:
        with self._lock:
            try:
                revision = self._registry.get_sample_revision(sample_revision_id)
            except KeyError:
                raise catalog_not_found() from None
            if revision.sample_id != sample_id:
                raise catalog_not_found()
            self._project(revision.project_id)
            return self._sample_revision_body(
                sample_revision_id,
                include_candidates=True,
            )

    def list_experiment_revisions(
        self,
        *,
        experiment_id: str,
        filters: CatalogFilters,
        cursor: str | None,
        limit: int,
    ) -> dict[str, Any]:
        with self._lock:
            experiment = self._experiments.get(experiment_id)
            if experiment is None:
                raise catalog_not_found()
            self._project(experiment.project_id)
            page, next_cursor = _page(
                experiment.revisions,
                kind="experiment_revisions",
                filters=filters,
                cursor=cursor,
                limit=limit,
                created_at=lambda item: item.created_at,
                object_id=lambda item: item.experiment_revision_id,
            )
            return {
                "experiment_revisions": [
                    self._experiment_detail(
                        item.experiment_revision_id,
                        include_candidates=False,
                    )
                    for item in page
                ],
                "next_cursor": next_cursor,
            }

    def get_experiment_revision(
        self,
        experiment_id: str,
        experiment_revision_id: str,
    ) -> dict[str, Any]:
        with self._lock:
            revision = self._experiment_revisions.get(experiment_revision_id)
            if revision is None or revision.experiment_id != experiment_id:
                raise catalog_not_found()
            self._project(revision.project_id)
            return self._experiment_detail(
                experiment_revision_id,
                include_candidates=True,
            )

    def list_sample_run_results(
        self,
        *,
        project_id: str,
        sample_revision_id: str,
        filters: CatalogFilters,
        cursor: str | None,
        limit: int,
    ) -> dict[str, Any]:
        del limit
        position = decode_catalog_cursor(
            cursor,
            kind="sample_runs",
            filters=filters,
        )
        with self._lock:
            self._project(project_id)
            try:
                revision = self._registry.get_sample_revision(sample_revision_id)
            except KeyError:
                raise catalog_not_found() from None
            if revision.project_id != project_id:
                raise catalog_not_found()
            if position is not None:
                raise catalog_cursor_not_found()
            return {
                "ok": True,
                "project_id": project_id,
                "sample_revision_id": sample_revision_id,
                "results": [],
                "next_cursor": None,
            }

    def list_experiment_run_results(
        self,
        *,
        project_id: str,
        experiment_revision_id: str,
        filters: CatalogFilters,
        cursor: str | None,
        limit: int,
    ) -> dict[str, Any]:
        del limit
        position = decode_catalog_cursor(
            cursor,
            kind="experiment_runs",
            filters=filters,
        )
        with self._lock:
            self._project(project_id)
            revision = self._experiment_revisions.get(experiment_revision_id)
            if revision is None or revision.project_id != project_id:
                raise catalog_not_found()
            if position is not None:
                raise catalog_cursor_not_found()
            return {
                "ok": True,
                "project_id": project_id,
                "experiment_revision_id": experiment_revision_id,
                "trace_status": "unresolved",
                "reason": "experiment_snapshot_binding_absent",
                "results": [],
                "next_cursor": None,
            }

    def experiment_revision_freeze_source(
        self,
        experiment_revision_id: str,
    ) -> ExperimentRevisionFreezeSource | None:
        with self._lock:
            revision = self._experiment_revisions.get(experiment_revision_id)
            if revision is None:
                return None
            experiment = self._experiments.get(revision.experiment_id)
            if experiment is None or experiment.project_id != revision.project_id:
                return None
            try:
                project = self._registry.get_project(revision.project_id)
            except KeyError:
                project = None
            kind = None if project is None else project.kind
            kind_value = (
                None
                if kind is None
                else kind.value
                if isinstance(kind, ProjectKind)
                else str(kind)
            )
            project_active = (
                project is not None
                and project.archived_at is None
                and project.project_id != LEGACY_PROJECT_ID
                and kind_value == "user"
            )
            return ExperimentRevisionFreezeSource(
                experiment_revision_id=revision.experiment_revision_id,
                project_id=revision.project_id,
                project_active=project_active,
                experiment_archived=experiment.archived_at is not None,
                payload_digest=build_experiment_revision_payload_digest(
                    revision.canonical_payload
                ),
                sample_revision_ids=tuple(
                    pin.sample_revision_id for pin in revision.payload.samples
                ),
            )

    def rebuild_sample_projection(self, sample_revision_id: str) -> None:
        with self._lock:
            try:
                revision = self._registry.get_sample_revision(sample_revision_id)
            except KeyError:
                raise catalog_not_found() from None
            before = revision.canonical_payload
            self._block_for_revision(sample_revision_id)
            if revision.canonical_payload != before:
                raise catalog_payload_invalid()

    def rebuild_experiment_projection(self, experiment_revision_id: str) -> None:
        with self._lock:
            revision = self._experiment_revisions.get(experiment_revision_id)
            if revision is None:
                raise catalog_not_found()
            before = revision.canonical_payload
            payload = interpret_stored_experiment_payload(before)
            if payload.canonical_json() != before:
                raise catalog_payload_invalid()
            if revision.canonical_payload != before:
                raise catalog_payload_invalid()

    def _write(
        self, write: CatalogWrite, apply: Callable[[], dict[str, Any]]
    ) -> dict[str, Any]:
        with self._lock:
            if not self._is_administrator(write.actor_user_id):
                raise administrator_required()
            receipt = self._receipts.get((write.actor_user_id, write.idempotency_key))
            if not write.target_key:
                raise ValueError("catalog write target is required")
            if receipt is not None:
                if (
                    receipt.request_sha256 != write.request_sha256
                    or receipt.operation != write.operation
                    or receipt.target_key != write.target_key
                ):
                    raise catalog_idempotency_conflict()
                return self._replay(receipt)
            return apply()

    def _replay(self, receipt: _Receipt) -> dict[str, Any]:
        if receipt.operation in {
            "catalogCreateTerm",
            "catalogReviseTerm",
            "catalogArchiveTerm",
        }:
            term = self._terms[receipt.subject_id]
            revision_id = (
                receipt.result_revision_id or term.revisions[-1].term_revision_id
            )
            return self._term_body(receipt.subject_id, revision_id)
        if receipt.operation == "catalogArchiveExperiment":
            experiment = self._experiments[receipt.subject_id]
            return {
                "experiment_id": experiment.experiment_id,
                "archived_at": _stamp(receipt.result_archived_at),
                "revision_number": experiment.revisions[-1].revision_number,
            }
        if receipt.operation in {"catalogCreateExperiment", "catalogReviseExperiment"}:
            if receipt.result_revision_id is None:
                raise catalog_not_found()
            return self._experiment_detail(
                receipt.result_revision_id,
                include_candidates=True,
            )
        if receipt.operation == "catalogReviseSampleCatalog":
            if receipt.result_revision_id is None:
                raise catalog_not_found()
            return self._sample_revision_body(
                receipt.result_revision_id,
                include_candidates=True,
            )
        raise catalog_idempotency_conflict()

    def _add_receipt(
        self,
        write: CatalogWrite,
        *,
        subject_id: str,
        result_revision_id: str | None,
        result_archived_at: datetime | None,
    ) -> None:
        self._receipts[(write.actor_user_id, write.idempotency_key)] = _Receipt(
            operation=write.operation,
            request_sha256=write.request_sha256,
            target_key=write.target_key,
            subject_id=subject_id,
            result_revision_id=result_revision_id,
            result_archived_at=result_archived_at,
        )

    def _project(self, project_id: str) -> None:
        try:
            project = self._registry.get_project(project_id)
        except KeyError:
            raise catalog_not_found() from None
        kind = (
            project.kind.value
            if isinstance(project.kind, ProjectKind)
            else str(project.kind)
        )
        if (
            project.archived_at is not None
            or project.project_id == LEGACY_PROJECT_ID
            or kind != "user"
        ):
            raise catalog_not_found()

    def _sample_project(self, sample_id: str) -> str:
        try:
            return self._registry.get_sample(sample_id).project_id
        except KeyError:
            raise catalog_not_found() from None

    def _experiment_payload(
        self,
        *,
        project_id: str,
        display_name: str,
        assay: object,
        target: object,
        batch: object,
        samples: object,
    ) -> ExperimentPayload:
        assay_slot, target_slot, batch_slot, requests = interpret_experiment_fields(
            assay,
            target,
            batch,
            samples,
            _MemoryTermLookup(self),
            for_write=True,
        )
        pins: list[ExperimentSamplePin] = []
        for ordinal, (request, role) in enumerate(requests):
            sample_revision_id = str(request["sample_revision_id"])
            try:
                row = self._registry.get_sample_revision(sample_revision_id)
            except KeyError:
                raise catalog_payload_invalid() from None
            if row.project_id != project_id:
                raise catalog_payload_invalid()
            group = request.get("replicate_group")
            if group == "":
                group = None
            pins.append(
                ExperimentSamplePin(
                    sample_id=row.sample_id,
                    sample_revision_id=row.sample_revision_id,
                    payload_digest=row.payload_digest,
                    ordinal=ordinal,
                    replicate_group=group if isinstance(group, str) else None,
                    replicate_role=role,
                )
            )
        return ExperimentPayload(
            display_name=display_name,
            assay=assay_slot,
            target=target_slot,
            batch=batch_slot,
            samples=tuple(pins),
        )

    def _experiment_matches(
        self, experiment: _Experiment, filters: CatalogFilters
    ) -> bool:
        if not filters.include_archived and experiment.archived_at is not None:
            return False
        head = experiment.revisions[-1]
        created = _utc(head.created_at)
        if filters.created_from is not None and created < filters.created_from:
            return False
        if filters.created_before is not None and created >= filters.created_before:
            return False
        slots = {
            "assay": head.payload.assay,
            "target": head.payload.target,
            "batch": head.payload.batch,
        }
        for name in EXPERIMENT_FIELDS:
            field = filters.field(name)
            if field is not None and not slot_matches(slots[name], field):
                return False
        current = filters.revision_view == "current_sample_catalog"
        for name in SAMPLE_FIELDS:
            field = filters.field(name)
            if field is None:
                continue
            if not any(
                block_matches(
                    self._pin_block(pin, current=current),
                    field,
                )
                for pin in head.payload.samples
            ):
                return False
        return True

    def _sample_matches(self, sample: Any, filters: CatalogFilters) -> bool:
        created = _utc(sample.created_at)
        if filters.created_from is not None and created < filters.created_from:
            return False
        if filters.created_before is not None and created >= filters.created_before:
            return False
        head_id, _number = self._registry.get_sample_revision_head(sample.sample_id)
        block = self._block_for_revision(head_id)
        if filters.catalog_state == "absent" and block is not None:
            return False
        if filters.catalog_state == "present" and block is None:
            return False
        for field in filters.fields:
            if field.name not in SAMPLE_FIELDS:
                continue
            if not block_matches(block, field):
                return False
        return True

    def _term_matches(self, term: _Term, filters: CatalogFilters) -> bool:
        if not filters.include_archived and term.archived_at is not None:
            return False
        if filters.term_type is not None and term.term_type != filters.term_type:
            return False
        created = _utc(term.created_at)
        if filters.created_from is not None and created < filters.created_from:
            return False
        if filters.created_before is not None and created >= filters.created_before:
            return False
        return True

    def _pin_block(
        self, pin: ExperimentSamplePin, *, current: bool
    ) -> SampleCatalogBlock | None:
        revision_id = pin.sample_revision_id
        if current:
            head_id, _number = self._registry.get_sample_revision_head(pin.sample_id)
            revision_id = head_id
        return self._block_for_revision(revision_id)

    def _block_for_revision(self, sample_revision_id: str) -> SampleCatalogBlock | None:
        if not self._registry.catalog_revision_trusted(sample_revision_id):
            return None
        revision = self._registry.get_sample_revision(sample_revision_id)
        decoded = json.loads(revision.canonical_payload)
        attributes = decoded["attributes"]
        if "catalog" not in attributes:
            return None
        return interpret_sample_catalog(
            attributes.get("catalog"),
            _MemoryTermLookup(self),
            for_write=False,
        )

    def _names(self, revision_ids: set[str]) -> dict[str, str]:
        lookup = _MemoryTermLookup(self)
        names: dict[str, str] = {}
        for revision_id in revision_ids:
            pin = lookup(revision_id)
            if pin is None:
                raise catalog_payload_invalid()
            names[revision_id] = pin.display_name
        return names

    def _sample_summary(self, sample_id: str) -> dict[str, Any]:
        head_id, _number = self._registry.get_sample_revision_head(sample_id)
        return self._sample_revision_body(head_id, include_candidates=False)

    def _sample_revision_body(
        self,
        sample_revision_id: str,
        *,
        include_candidates: bool,
    ) -> dict[str, Any]:
        revision = self._registry.get_sample_revision(sample_revision_id)
        decoded = json.loads(revision.canonical_payload)
        block = self._block_for_revision(sample_revision_id)
        material = render_sample_material(
            block,
            self._names(collect_term_revision_ids(block)),
            include_candidates=include_candidates,
        )
        return {
            "sample_id": revision.sample_id,
            "sample_revision_id": revision.sample_revision_id,
            "revision_number": revision.revision_number,
            "display_name": decoded["display_name"],
            "created_at": _stamp(revision.created_at),
            **material,
        }

    def _experiment_detail(
        self,
        experiment_revision_id: str,
        *,
        include_candidates: bool,
        revision_view: str = "linked",
    ) -> dict[str, Any]:
        revision = self._experiment_revisions[experiment_revision_id]
        payload = revision.payload
        name_ids = collect_experiment_term_revision_ids(payload)
        current = revision_view == "current_sample_catalog"
        blocks = [self._pin_block(pin, current=current) for pin in payload.samples]
        for block in blocks:
            name_ids.update(collect_term_revision_ids(block))
        names = self._names(name_ids)
        samples = [
            self._render_pin(
                pin,
                block,
                names,
                revision_view=revision_view,
                include_candidates=include_candidates,
            )
            for pin, block in zip(payload.samples, blocks, strict=True)
        ]
        return {
            "experiment_id": revision.experiment_id,
            "experiment_revision_id": revision.experiment_revision_id,
            "revision_number": revision.revision_number,
            "display_name": payload.display_name,
            "assay": render_single(
                payload.assay,
                names,
                include_candidates=include_candidates,
            ),
            "target": render_single(
                payload.target,
                names,
                include_candidates=include_candidates,
            ),
            "batch": render_single(
                payload.batch,
                names,
                include_candidates=include_candidates,
            ),
            "samples": samples,
        }

    def _render_pin(
        self,
        pin: ExperimentSamplePin,
        block: SampleCatalogBlock | None,
        names: Mapping[str, str],
        *,
        revision_view: str,
        include_candidates: bool,
    ) -> dict[str, Any]:
        linked_id = pin.sample_revision_id
        current_id = linked_id
        if revision_view == "current_sample_catalog":
            current_id, _number = self._registry.get_sample_revision_head(pin.sample_id)
        interpreted_id = (
            current_id if revision_view == "current_sample_catalog" else linked_id
        )
        material = render_sample_material(block, names, include_candidates=False)
        if material.get("catalog_state") == "present":
            material.pop("catalog_state", None)
        body = {
            "sample_id": pin.sample_id,
            "sample_revision_id": interpreted_id,
            "ordinal": pin.ordinal,
            "replicate_group": pin.replicate_group,
            "replicate_role": render_single(
                pin.replicate_role,
                names,
                include_candidates=include_candidates,
            ),
        }
        body.update(material)
        if interpreted_id != linked_id:
            body["linked_sample_revision_id"] = linked_id
            body["current_sample_revision_id"] = current_id
        return body

    def _term_body(self, term_id: str, term_revision_id: str) -> dict[str, Any]:
        term = self._terms[term_id]
        revision = self._term_revisions[term_revision_id]
        return {
            "term_id": term.term_id,
            "term_type": term.term_type,
            "archived_at": _stamp(term.archived_at),
            "revision": {
                "term_revision_id": revision.term_revision_id,
                "revision_number": revision.revision_number,
                "display_name": revision.display_name,
                "namespace": revision.namespace,
                "accession": revision.accession,
                "aliases": list(revision.aliases),
                "change_kind": revision.change_kind,
                "created_at": _stamp(revision.created_at),
            },
        }

    @staticmethod
    def _alias_json(aliases: tuple[str, ...]) -> str:
        from encode_pipeline.platform.experiment_catalog import canonical_alias_json

        return canonical_alias_json(validate_aliases(list(aliases)))


def _expect_revision(current: int, expected: object) -> None:
    if isinstance(expected, bool) or not isinstance(expected, int) or expected < 1:
        raise catalog_payload_invalid()
    if current != expected:
        raise catalog_revision_conflict()


def _required_project(filters: CatalogFilters) -> str:
    if filters.project_id is None:
        raise catalog_not_found()
    return filters.project_id


def _page(
    items: Sequence[Any],
    *,
    kind: str,
    filters: CatalogFilters,
    cursor: str | None,
    limit: int,
    created_at: Callable[[Any], datetime],
    object_id: Callable[[Any], str],
) -> tuple[list[Any], str | None]:
    position = decode_catalog_cursor(cursor, kind=kind, filters=filters)
    ordered = sorted(items, key=lambda item: (_utc(created_at(item)), object_id(item)))
    if position is not None:
        hit = next(
            (item for item in ordered if object_id(item) == position.object_id),
            None,
        )
        if hit is None or _utc(created_at(hit)) != _utc(position.created_at):
            raise catalog_cursor_not_found()
        marker = (_utc(position.created_at), position.object_id)
        ordered = [
            item
            for item in ordered
            if (_utc(created_at(item)), object_id(item)) > marker
        ]
    page = list(ordered[:limit])
    if len(ordered) <= limit:
        return page, None
    last = page[-1]
    return page, encode_catalog_cursor(
        kind=kind,
        filters=filters,
        position=CatalogCursorPosition(
            created_at=_utc(created_at(last)),
            object_id=object_id(last),
        ),
    )


def _utc(value: datetime) -> datetime:
    if value.tzinfo is None or value.utcoffset() is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _stamp(value: datetime | None) -> str | None:
    if value is None:
        return None
    text_value = _utc(value).isoformat(timespec="microseconds")
    if text_value.endswith(".000000+00:00"):
        text_value = text_value.replace(".000000+00:00", "+00:00")
    return text_value.replace("+00:00", "Z")
