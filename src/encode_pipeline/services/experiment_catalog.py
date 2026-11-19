"""Orchestration for manual experiment-catalog maintenance."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from datetime import datetime, timezone
from hashlib import sha256
from typing import Any
from uuid import uuid4

from encode_pipeline.platform.experiment_catalog import (
    CatalogWrite,
    catalog_payload_invalid,
    parse_page_filters,
    parse_page_limit,
    validate_idempotency_key,
)


def _prefixed_id(prefix: str) -> str:
    return f"{prefix}{uuid4().hex}"


class ExperimentCatalogService:
    """Assign identities and delegate Slice A catalog commands to a repository."""

    def __init__(
        self,
        repository: Any,
        *,
        now_factory: Callable[[], datetime] | None = None,
        term_id_factory: Callable[[], str] | None = None,
        term_revision_id_factory: Callable[[], str] | None = None,
        experiment_id_factory: Callable[[], str] | None = None,
        experiment_revision_id_factory: Callable[[], str] | None = None,
        sample_revision_id_factory: Callable[[], str] | None = None,
    ) -> None:
        self._repository = repository
        self._now_factory = now_factory or (lambda: datetime.now(timezone.utc))
        self._term_id_factory = term_id_factory or (lambda: _prefixed_id("term_"))
        self._term_revision_id_factory = term_revision_id_factory or (
            lambda: _prefixed_id("trmv_")
        )
        self._experiment_id_factory = experiment_id_factory or (
            lambda: _prefixed_id("exp_")
        )
        self._experiment_revision_id_factory = experiment_revision_id_factory or (
            lambda: _prefixed_id("expr_")
        )
        self._sample_revision_id_factory = sample_revision_id_factory or (
            lambda: _prefixed_id("smpr_")
        )

    def create_term(
        self,
        *,
        actor_user_id: str,
        idempotency_key: str | None,
        raw_body: bytes,
        term_type: str,
        display_name: str,
        namespace: str | None,
        accession: str | None,
        aliases: Sequence[str],
    ) -> dict[str, Any]:
        return self._repository.create_term(
            self._write(actor_user_id, idempotency_key, raw_body, "catalogCreateTerm"),
            term_id=self._term_id_factory(),
            term_revision_id=self._term_revision_id_factory(),
            term_type=term_type,
            display_name=display_name,
            namespace=namespace,
            accession=accession,
            aliases=tuple(aliases),
        )

    def revise_term(
        self,
        *,
        actor_user_id: str,
        idempotency_key: str | None,
        raw_body: bytes,
        term_id: str,
        expected_revision_number: object,
        change_kind: str,
        display_name: str,
        namespace: str | None,
        accession: str | None,
        aliases: Sequence[str],
    ) -> dict[str, Any]:
        if change_kind != "editorial":
            raise catalog_payload_invalid()
        return self._repository.revise_term(
            self._write(actor_user_id, idempotency_key, raw_body, "catalogReviseTerm"),
            term_id=term_id,
            term_revision_id=self._term_revision_id_factory(),
            expected_revision_number=expected_revision_number,
            display_name=display_name,
            namespace=namespace,
            accession=accession,
            aliases=tuple(aliases),
        )

    def archive_term(
        self,
        *,
        actor_user_id: str,
        idempotency_key: str | None,
        raw_body: bytes,
        term_id: str,
        expected_revision_number: object,
    ) -> dict[str, Any]:
        return self._repository.archive_term(
            self._write(actor_user_id, idempotency_key, raw_body, "catalogArchiveTerm"),
            term_id=term_id,
            expected_revision_number=expected_revision_number,
        )

    def create_experiment(
        self,
        *,
        actor_user_id: str,
        idempotency_key: str | None,
        raw_body: bytes,
        project_id: str,
        stable_key: str,
        display_name: str,
        assay: object,
        target: object,
        batch: object,
        samples: object,
    ) -> dict[str, Any]:
        return self._repository.create_experiment(
            self._write(
                actor_user_id,
                idempotency_key,
                raw_body,
                "catalogCreateExperiment",
            ),
            project_id=project_id,
            experiment_id=self._experiment_id_factory(),
            experiment_revision_id=self._experiment_revision_id_factory(),
            stable_key=stable_key,
            display_name=display_name,
            assay=assay,
            target=target,
            batch=batch,
            samples=samples,
        )

    def revise_experiment(
        self,
        *,
        actor_user_id: str,
        idempotency_key: str | None,
        raw_body: bytes,
        experiment_id: str,
        expected_revision_number: object,
        display_name: str,
        assay: object,
        target: object,
        batch: object,
        samples: object,
    ) -> dict[str, Any]:
        return self._repository.revise_experiment(
            self._write(
                actor_user_id,
                idempotency_key,
                raw_body,
                "catalogReviseExperiment",
            ),
            experiment_id=experiment_id,
            experiment_revision_id=self._experiment_revision_id_factory(),
            expected_revision_number=expected_revision_number,
            display_name=display_name,
            assay=assay,
            target=target,
            batch=batch,
            samples=samples,
        )

    def archive_experiment(
        self,
        *,
        actor_user_id: str,
        idempotency_key: str | None,
        raw_body: bytes,
        experiment_id: str,
        expected_revision_number: object,
    ) -> dict[str, Any]:
        return self._repository.archive_experiment(
            self._write(
                actor_user_id,
                idempotency_key,
                raw_body,
                "catalogArchiveExperiment",
            ),
            experiment_id=experiment_id,
            expected_revision_number=expected_revision_number,
        )

    def revise_sample_catalog(
        self,
        *,
        actor_user_id: str,
        idempotency_key: str | None,
        raw_body: bytes,
        sample_id: str,
        expected_revision_number: object,
        display_name: str | None,
        catalog: object,
    ) -> dict[str, Any]:
        return self._repository.revise_sample_catalog(
            self._write(
                actor_user_id,
                idempotency_key,
                raw_body,
                "catalogReviseSampleCatalog",
            ),
            sample_id=sample_id,
            sample_revision_id=self._sample_revision_id_factory(),
            expected_revision_number=expected_revision_number,
            display_name=display_name,
            catalog=catalog,
        )

    def list_experiments(self, **query: Any) -> dict[str, Any]:
        return self._repository.list_experiments(**self._page("experiments", query))

    def list_samples(self, **query: Any) -> dict[str, Any]:
        return self._repository.list_samples(**self._page("samples", query))

    def list_terms(self, **query: Any) -> dict[str, Any]:
        return self._repository.list_terms(**self._page("terms", query))

    def get_term(self, term_id: str) -> dict[str, Any]:
        return self._repository.get_term(term_id)

    def list_sample_revisions(self, *, sample_id: str, **query: Any) -> dict[str, Any]:
        page = self._page("sample_revisions", query, owner_id=sample_id)
        return self._repository.list_sample_revisions(sample_id=sample_id, **page)

    def get_sample_revision(
        self, sample_id: str, sample_revision_id: str
    ) -> dict[str, Any]:
        return self._repository.get_sample_revision(sample_id, sample_revision_id)

    def list_experiment_revisions(
        self, *, experiment_id: str, **query: Any
    ) -> dict[str, Any]:
        page = self._page("experiment_revisions", query, owner_id=experiment_id)
        return self._repository.list_experiment_revisions(
            experiment_id=experiment_id,
            **page,
        )

    def get_experiment_revision(
        self, experiment_id: str, experiment_revision_id: str
    ) -> dict[str, Any]:
        return self._repository.get_experiment_revision(
            experiment_id,
            experiment_revision_id,
        )

    def list_sample_run_results(
        self, *, project_id: str, sample_revision_id: str, **query: Any
    ) -> dict[str, Any]:
        page = self._page("sample_runs", query, project_id=project_id)
        return self._repository.list_sample_run_results(
            project_id=project_id,
            sample_revision_id=sample_revision_id,
            **page,
        )

    def list_experiment_run_results(
        self, *, project_id: str, experiment_revision_id: str, **query: Any
    ) -> dict[str, Any]:
        page = self._page("experiment_runs", query, project_id=project_id)
        return self._repository.list_experiment_run_results(
            project_id=project_id,
            experiment_revision_id=experiment_revision_id,
            **page,
        )

    def _write(
        self,
        actor_user_id: str,
        idempotency_key: str | None,
        raw_body: bytes,
        operation: str,
    ) -> CatalogWrite:
        if not isinstance(raw_body, bytes):
            raise catalog_payload_invalid()
        return CatalogWrite(
            actor_user_id=actor_user_id,
            idempotency_key=validate_idempotency_key(idempotency_key),
            request_sha256=sha256(raw_body).hexdigest(),
            operation=operation,
            now=self._now_factory(),
        )

    def _page(
        self,
        kind: str,
        query: Mapping[str, Any],
        *,
        project_id: str | None = None,
        owner_id: str | None = None,
    ) -> dict[str, Any]:
        filters = parse_page_filters(
            kind=kind,
            project_id=project_id
            if project_id is not None
            else query.get("project_id"),
            revision_view=query.get("revision_view"),
            include_archived=bool(query.get("include_archived", False)),
            term_type=query.get("term_type"),
            created_from=query.get("created_from"),
            created_before=query.get("created_before"),
            catalog_state=query.get("catalog_state"),
            field_terms=query.get("field_terms") or {},
            field_statuses=query.get("field_statuses") or {},
            owner_id=owner_id,
        )
        return {
            "filters": filters,
            "cursor": query.get("cursor"),
            "limit": parse_page_limit(query.get("limit")),
        }
