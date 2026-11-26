"""Orchestration for manual experiment-catalog maintenance."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from datetime import datetime, timezone
from hashlib import sha256
import json
from typing import Any
from uuid import uuid4

from encode_pipeline.platform.catalog_proposal_adoption import (
    assert_request_within_limit,
)
from encode_pipeline.platform.experiment_catalog import (
    CatalogError,
    CatalogWrite,
    catalog_experiment_archived,
    catalog_not_found,
    catalog_payload_invalid,
    catalog_revision_conflict,
    catalog_term_archived,
    parse_page_filters,
    parse_page_limit,
    validate_idempotency_key,
)
from encode_pipeline.services.peri_runtime.catalog_classify import (
    CATALOG_FIELD,
    CATALOG_REASON,
    INSUFFICIENT_EVIDENCE,
    MAX_CLASSIFY_CANDIDATES,
    MODEL_PROPOSAL_SOURCE,
    BoundCandidate,
    CatalogTransportError,
    bind_target_candidates,
    clue_is_auditable,
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
        adoption_receipt_id_factory: Callable[[], str] | None = None,
        catalog_model: Any | None = None,
    ) -> None:
        self._repository = repository
        self._catalog_model = catalog_model
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
        self._adoption_receipt_id_factory = adoption_receipt_id_factory or (
            lambda: _prefixed_id("cadopt_")
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

    def register_proposal(
        self,
        *,
        actor_user_id: str,
        idempotency_key: str | None,
        raw_body: bytes,
        project_id: str,
        source: str,
        items: object,
        rule_version: str | None,
        vocabulary_version: str | None,
        projection_version: str | None,
    ) -> dict[str, Any]:
        assert_request_within_limit(raw_body)
        if source == MODEL_PROPOSAL_SOURCE:
            raise CatalogError(
                422,
                "catalog_payload_invalid",
                "Offline registration cannot claim a model source.",
                path="source",
            )
        return self._repository.register_proposal(
            self._write(
                actor_user_id,
                idempotency_key,
                raw_body,
                "catalogRegisterProposal",
            ),
            project_id=project_id,
            source=source,
            items=items,
            rule_version=rule_version,
            vocabulary_version=vocabulary_version,
            projection_version=projection_version,
        )

    def generate_classification_suggestion(
        self,
        *,
        actor_user_id: str,
        idempotency_key: str | None,
        raw_body: bytes,
        project_id: str,
        experiment_id: str,
        experiment_revision_id: str,
        field: str,
        reason_code: str,
        clue: str,
        vocabulary_term_revision_id: str | None,
    ) -> dict[str, Any]:
        """Ask the bound model for one target term, then store a pending proposal.

        Generation does not preview, select, or adopt, and it does not append
        an experiment revision.
        """

        assert_request_within_limit(raw_body)
        if self._catalog_model is None:
            raise CatalogError(
                409,
                "catalog_model_unbound",
                "The owner has no bound model for catalog classification.",
            )
        if not clue_is_auditable(clue):
            raise CatalogError(
                422,
                "catalog_payload_invalid",
                "Classification clue must name evidence, not an ordinal.",
                path="clue",
            )
        if field != "target" or reason_code != "TARGET_TERM_MATCH":
            raise CatalogError(
                422,
                "catalog_payload_invalid",
                "This slice only classifies experiment target.",
                path="field",
            )
        write = self._write(
            actor_user_id,
            idempotency_key,
            raw_body,
            "catalogGenerateClassificationSuggestion",
        )
        context = self._repository.classification_context(experiment_id)
        if context["project_id"] != project_id:
            raise catalog_not_found()
        if context["archived"]:
            raise catalog_experiment_archived()
        if context["experiment_revision_id"] != experiment_revision_id:
            raise catalog_revision_conflict()
        terms = self._classification_terms(context, vocabulary_term_revision_id)
        candidates = bind_target_candidates(
            terms, explicit_revision_id=vocabulary_term_revision_id
        )
        if vocabulary_term_revision_id is not None and len(candidates) != 1:
            raise CatalogError(
                422,
                "catalog_payload_invalid",
                "The selected target term has no auditable evidence token.",
                path="vocabulary_term_revision_id",
            )
        if vocabulary_term_revision_id is None and (
            context["target_terms_truncated"]
            or len(candidates) > MAX_CLASSIFY_CANDIDATES
            or not candidates
        ):
            raise CatalogError(
                422,
                "catalog_payload_invalid",
                "Select one current target term revision to classify.",
                path="vocabulary_term_revision_id",
            )
        if len(candidates) > MAX_CLASSIFY_CANDIDATES:
            candidates = candidates[:MAX_CLASSIFY_CANDIDATES]
        projection = {
            "field": CATALOG_FIELD,
            "reason_code": CATALOG_REASON,
            "clue": clue,
            "candidates": [
                {"token": item.token, "label": item.label} for item in candidates
            ],
        }
        digest = sha256(
            json.dumps(
                {
                    "project_id": project_id,
                    "experiment_id": experiment_id,
                    "experiment_revision_id": experiment_revision_id,
                    "clue": clue,
                    "vocabulary_term_revision_id": vocabulary_term_revision_id,
                    "candidates": [
                        {
                            "token": item.token,
                            "term_revision_id": item.term_revision_id,
                        }
                        for item in candidates
                    ],
                },
                sort_keys=True,
                separators=(",", ":"),
            ).encode()
        ).hexdigest()
        try:
            decision = self._catalog_model.classify(
                project_id=project_id,
                idempotency_key=write.idempotency_key,
                clue=clue,
                candidates=tuple((item.token, item.label) for item in candidates),
                projection_digest=digest,
            )
        except CatalogTransportError as exc:
            if exc.code == "catalog_idempotency_conflict":
                from encode_pipeline.platform.experiment_catalog import (
                    catalog_idempotency_conflict,
                )

                raise catalog_idempotency_conflict() from None
            raise CatalogError(409, exc.code, exc.message) from None
        confirmed = self._repository.classification_context(experiment_id)
        if (
            confirmed["project_id"] != project_id
            or confirmed["experiment_revision_id"] != experiment_revision_id
            or confirmed["revision_number"] != context["revision_number"]
        ):
            raise catalog_revision_conflict()
        body = _suggestion_envelope(
            context=context,
            clue=clue,
            projection=projection,
            decision=decision,
            proposal=None,
        )
        if decision == INSUFFICIENT_EVIDENCE:
            body["outcome"] = "insufficient_evidence"
            return body
        chosen = next((item for item in candidates if item.token == decision), None)
        if chosen is None:
            raise CatalogError(
                409,
                "catalog_model_malformed",
                "The model response was not a bound candidate or insufficient evidence.",
            )
        proposal = self._repository.register_proposal(
            write,
            project_id=project_id,
            source=MODEL_PROPOSAL_SOURCE,
            items=[_target_proposal_item(context, chosen, clue)],
            rule_version=context["rule_version"],
            vocabulary_version=context["vocabulary_version"],
            projection_version=context["projection_version"],
        )
        final = self._repository.classification_context(experiment_id)
        if (
            final["experiment_revision_id"] != experiment_revision_id
            or final["revision_number"] != context["revision_number"]
        ):
            raise catalog_revision_conflict()
        body["outcome"] = "proposal"
        body["decision"] = decision
        body["proposal"] = proposal
        return body

    def _classification_terms(
        self,
        context: Mapping[str, Any],
        vocabulary_term_revision_id: str | None,
    ) -> list[dict[str, Any]]:
        if vocabulary_term_revision_id is None:
            terms = context["target_terms"]
            if not isinstance(terms, list):
                raise catalog_payload_invalid()
            return terms
        found = self._repository.lookup_target_term(vocabulary_term_revision_id)
        if found is None:
            raise CatalogError(
                422,
                "catalog_payload_invalid",
                "The selected vocabulary revision is not a target term.",
                path="vocabulary_term_revision_id",
            )
        if found["archived"]:
            raise catalog_term_archived()
        if not found["requested_is_head"]:
            raise catalog_revision_conflict()
        return [found]

    def get_proposal(
        self,
        *,
        proposal_id: str,
        include_candidates: bool,
    ) -> dict[str, Any]:
        return self._repository.get_proposal(
            proposal_id=proposal_id,
            include_candidates=include_candidates,
            now=self._now_factory(),
        )

    def preview_proposal(
        self,
        *,
        actor_user_id: str,
        raw_body: bytes,
        proposal_id: str,
        selected_item_ids: object,
    ) -> dict[str, Any]:
        assert_request_within_limit(raw_body)
        return self._repository.preview_proposal(
            actor_user_id=actor_user_id,
            proposal_id=proposal_id,
            selected_item_ids=selected_item_ids,
            now=self._now_factory(),
        )

    def adopt_proposal(
        self,
        *,
        actor_user_id: str,
        idempotency_key: str | None,
        raw_body: bytes,
        proposal_id: str,
        selected_item_ids: object,
        expected_revisions: object,
        preview_digest: object,
    ) -> dict[str, Any]:
        assert_request_within_limit(raw_body)
        return self._repository.adopt_proposal(
            self._write(
                actor_user_id,
                idempotency_key,
                raw_body,
                "catalogAdoptProposal",
            ),
            proposal_id=proposal_id,
            selected_item_ids=selected_item_ids,
            expected_revisions=expected_revisions,
            preview_digest=preview_digest,
            experiment_revision_id_factory=self._experiment_revision_id_factory,
            sample_revision_id_factory=self._sample_revision_id_factory,
            receipt_id_factory=self._adoption_receipt_id_factory,
        )

    def list_adoption_receipts(
        self,
        *,
        actor_user_id: str,
        project_id: str,
        proposal_id: str | None,
    ) -> dict[str, Any]:
        return self._repository.list_adoption_receipts(
            actor_user_id=actor_user_id,
            project_id=project_id,
            proposal_id=proposal_id,
        )

    def set_active_curation_versions(
        self,
        *,
        rule_version: str,
        vocabulary_version: str,
        projection_version: str,
    ) -> None:
        self._repository.set_active_curation_versions(
            rule_version=rule_version,
            vocabulary_version=vocabulary_version,
            projection_version=projection_version,
        )

    def get_publication_experiment_trace(
        self,
        *,
        project_id: str,
        run_id: str,
        artifact_id: str,
        artifact_generation: str,
    ) -> dict[str, Any]:
        return self._repository.get_publication_experiment_trace(
            project_id=project_id,
            run_id=run_id,
            artifact_id=artifact_id,
            artifact_generation=artifact_generation,
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


def _suggestion_envelope(
    *,
    context: Mapping[str, Any],
    clue: str,
    projection: dict[str, Any],
    decision: str,
    proposal: dict[str, Any] | None,
) -> dict[str, Any]:
    return {
        "ok": True,
        "outcome": "insufficient_evidence",
        "field": "target",
        "reason_code": "TARGET_TERM_MATCH",
        "experiment_id": context["experiment_id"],
        "experiment_revision_id": context["experiment_revision_id"],
        "experiment_revision_number": context["revision_number"],
        "catalog_revision_changed": False,
        "projection": projection,
        "decision": decision,
        "proposal": proposal,
    }


def _target_proposal_item(
    context: Mapping[str, Any], chosen: BoundCandidate, clue: str
) -> dict[str, Any]:
    return {
        "item_id": "target",
        "target_kind": "experiment",
        "target_id": context["experiment_id"],
        "base_revision_id": context["experiment_revision_id"],
        "base_revision_number": context["revision_number"],
        "allowed_fields": ["target"],
        "before": {"target": context["target"]},
        "after": {
            "target": {
                "status": "known",
                "term_id": chosen.term_id,
                "term_revision_id": chosen.term_revision_id,
            }
        },
        "citations": {
            "rule_version": context["rule_version"],
            "vocabulary_version": context["vocabulary_version"],
            "reason_codes": ["TARGET_TERM_MATCH"],
            "source_refs": [
                {
                    "kind": "experiment_revision",
                    "ref_id": context["experiment_revision_id"],
                    "evidence_token": clue,
                },
                {
                    "kind": "term_revision",
                    "ref_id": chosen.term_revision_id,
                    "term_revision_id": chosen.term_revision_id,
                    "evidence_token": chosen.token,
                },
            ],
        },
    }
