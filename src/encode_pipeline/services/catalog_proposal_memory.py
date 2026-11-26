"""In-memory proposal adoption that follows the SQL repository's contract."""

from __future__ import annotations

from collections.abc import Callable
from copy import deepcopy
from datetime import datetime
import json
from typing import Any

from encode_pipeline.platform.catalog_proposal_adoption import (
    ActiveCurationVersions,
    SourceCatalog,
    TargetHead,
    adoption_receipt,
    applied_target,
    build_preview,
    build_proposal_view,
    build_registered_proposal,
    ensure_current_model_target,
    overlay_experiment_fields,
    overlay_sample_catalog,
    prepare_adoption,
    proposal_document,
    proposal_from_document,
    replayed_receipt,
)
from encode_pipeline.platform.data_registry import (
    SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
    SampleRevision,
    build_sample_revision_payload_digest,
    canonical_sample_revision_payload,
)
from encode_pipeline.platform.experiment_catalog import (
    CatalogError,
    CatalogWrite,
    catalog_not_found,
    catalog_payload_invalid,
    catalog_revision_conflict,
    interpret_sample_catalog,
    sample_catalog_object,
    validate_display_name,
)
from encode_pipeline.services.data_registry_repositories import (
    DataRegistryConflictError,
)
from encode_pipeline.services.experiment_catalog_repositories import (
    _ExperimentRevision,
    _MemoryTermLookup,
    _memory_target_term,
)


class _MemorySource(SourceCatalog):
    def __init__(self, repository) -> None:
        self._repository = repository

    def cited_experiment(self, experiment_id: str) -> tuple[str, bool] | None:
        experiment = self._repository._experiments.get(experiment_id)
        if experiment is None:
            return None
        return experiment.project_id, experiment.archived_at is not None

    def cited_experiment_revision(self, revision_id: str) -> tuple[str, bool] | None:
        revision = self._repository._experiment_revisions.get(revision_id)
        if revision is None:
            return None
        experiment = self._repository._experiments.get(revision.experiment_id)
        if experiment is None:
            return None
        return experiment.project_id, experiment.archived_at is not None

    def cited_sample(self, sample_id: str) -> str | None:
        try:
            return self._repository._registry.get_sample(sample_id).project_id
        except KeyError:
            return None

    def cited_sample_revision(self, revision_id: str) -> str | None:
        try:
            return self._repository._registry.get_sample_revision(
                revision_id
            ).project_id
        except KeyError:
            return None

    def cited_term_archived(self, revision_id: str) -> bool | None:
        revision = self._repository._term_revisions.get(revision_id)
        if revision is None:
            return None
        term = self._repository._terms.get(revision.term_id)
        if term is None:
            return None
        return term.archived_at is not None


def register_proposal(
    repository,
    write: CatalogWrite,
    *,
    project_id: str,
    proposal_id: str,
    source: str,
    items: object,
    rule_version: str | None,
    vocabulary_version: str | None,
    projection_version: str | None,
) -> dict[str, Any]:
    repository._project(project_id)
    proposal = build_registered_proposal(
        proposal_id=proposal_id,
        project_id=project_id,
        created_by_user_id=write.actor_user_id,
        source=source,
        items=items,
        created_at=write.now,
        active=repository._active_versions,
        rule_version=rule_version,
        vocabulary_version=vocabulary_version,
        projection_version=projection_version,
    )
    heads = _heads(repository, proposal)
    if proposal.source == "peri_bound_model":
        ensure_current_model_target(
            proposal,
            heads,
            repository._active_versions,
            lambda term_id: _memory_target_term(repository._terms[term_id]),
        )
    view = build_proposal_view(
        proposal,
        heads=heads,
        adopted_receipt=None,
        include_candidates=True,
        now=write.now,
        lookup=_MemoryTermLookup(repository),
    )
    repository._proposals[proposal.proposal_id] = {
        "document": proposal_document(proposal),
        "response": deepcopy(view),
        "project_id": project_id,
    }
    repository._add_receipt(
        write,
        subject_id=project_id,
        result_revision_id=None,
        result_archived_at=None,
    )
    return view


def get_proposal(
    repository,
    *,
    proposal_id: str,
    include_candidates: bool,
    now: datetime,
) -> dict[str, Any]:
    proposal = _load(repository, proposal_id)
    repository._project(proposal.project_id)
    return build_proposal_view(
        proposal,
        heads=_heads(repository, proposal),
        adopted_receipt=_adoption(repository, proposal.proposal_id),
        include_candidates=include_candidates,
        now=now,
        lookup=_MemoryTermLookup(repository),
    )


def preview_proposal(
    repository,
    *,
    proposal_id: str,
    selected_item_ids: object,
    now: datetime,
) -> dict[str, Any]:
    proposal = _load(repository, proposal_id)
    repository._project(proposal.project_id)
    if _adoption(repository, proposal.proposal_id) is not None:
        raise CatalogError(
            409,
            "catalog_proposal_invalidated",
            "This proposal has already been adopted.",
        )
    return build_preview(
        proposal,
        selected_item_ids=selected_item_ids,
        heads=_heads(repository, proposal),
        active=repository._active_versions,
        now=now,
        sources=_MemorySource(repository),
        lookup=_MemoryTermLookup(repository),
    )


def adopt_proposal(
    repository,
    write: CatalogWrite,
    *,
    proposal_id: str,
    selected_item_ids: object,
    expected_revisions: object,
    preview_digest: object,
    experiment_revision_id_factory: Callable[[], str],
    sample_revision_id_factory: Callable[[], str],
    receipt_id_factory: Callable[[], str],
) -> dict[str, Any]:
    proposal = _load(repository, proposal_id)
    repository._project(proposal.project_id)
    selected = prepare_adoption(
        proposal,
        selected_item_ids=selected_item_ids,
        expected_revisions=expected_revisions,
        preview_digest=preview_digest,
        heads=_heads(repository, proposal),
        active=repository._active_versions,
        now=write.now,
        adopted=_adoption(repository, proposal.proposal_id) is not None,
        sources=_MemorySource(repository),
        lookup=_MemoryTermLookup(repository),
    )
    if not isinstance(preview_digest, str):
        raise catalog_payload_invalid()
    applied = [
        _apply_item(
            repository,
            write,
            item,
            experiment_revision_id_factory=experiment_revision_id_factory,
            sample_revision_id_factory=sample_revision_id_factory,
        )
        for item in selected
    ]
    if proposal.proposal_id in repository._adoptions:
        raise CatalogError(
            409,
            "catalog_proposal_invalidated",
            "This proposal has already been adopted.",
        )
    body = adoption_receipt(
        proposal,
        selected=selected,
        preview_digest=preview_digest,
        actor_user_id=write.actor_user_id,
        idempotency_key=write.idempotency_key,
        receipt_id=receipt_id_factory(),
        applied_targets=applied,
        now=write.now,
    )
    repository._adoptions[proposal.proposal_id] = deepcopy(body)
    repository._add_receipt(
        write,
        subject_id=proposal.project_id,
        result_revision_id=None,
        result_archived_at=None,
    )
    return body


def list_adoption_receipts(
    repository,
    *,
    project_id: str,
    proposal_id: str | None,
) -> dict[str, Any]:
    repository._project(project_id)
    receipts = [
        deepcopy(body)
        for body in repository._adoptions.values()
        if body["project_id"] == project_id
        and (proposal_id is None or body["proposal_id"] == proposal_id)
    ]
    receipts.sort(key=lambda item: (item["adopted_at"], item["receipt_id"]))
    return {"ok": True, "project_id": project_id, "receipts": receipts}


def replay_register(repository, target_key: str) -> dict[str, Any]:
    stored = repository._proposals.get(target_key)
    if stored is None:
        from encode_pipeline.platform.experiment_catalog import (
            catalog_idempotency_conflict,
        )

        raise catalog_idempotency_conflict()
    return deepcopy(stored["response"])


def replay_adopt(repository, target_key: str) -> dict[str, Any]:
    stored = repository._adoptions.get(target_key)
    if stored is None:
        from encode_pipeline.platform.experiment_catalog import (
            catalog_idempotency_conflict,
        )

        raise catalog_idempotency_conflict()
    return replayed_receipt(deepcopy(stored))


def set_active_versions(
    repository,
    *,
    rule_version: str,
    vocabulary_version: str,
    projection_version: str,
) -> None:
    from encode_pipeline.platform.catalog_proposal_adoption import _require_version

    _require_version(rule_version)
    _require_version(vocabulary_version)
    _require_version(projection_version)
    repository._active_versions = ActiveCurationVersions(
        rule_version=rule_version,
        vocabulary_version=vocabulary_version,
        projection_version=projection_version,
    )


def _apply_item(
    repository,
    write: CatalogWrite,
    item,
    *,
    experiment_revision_id_factory: Callable[[], str],
    sample_revision_id_factory: Callable[[], str],
) -> dict[str, Any]:
    if item.target_kind == "experiment":
        return _apply_experiment(
            repository,
            write,
            item,
            experiment_revision_id_factory(),
        )
    if item.target_kind == "sample_catalog":
        return _apply_sample(
            repository,
            write,
            item,
            sample_revision_id_factory(),
        )
    raise catalog_payload_invalid()


def _apply_experiment(repository, write: CatalogWrite, item, revision_id: str):
    experiment = repository._experiments.get(item.target_id)
    if experiment is None:
        raise catalog_not_found()
    current = experiment.revisions[-1]
    if revision_id in repository._experiment_revisions:
        raise catalog_revision_conflict()
    fields = overlay_experiment_fields(current.canonical_payload, item)
    payload = repository._experiment_payload(
        project_id=experiment.project_id,
        display_name=validate_display_name(fields["display_name"]),
        assay=fields["assay"],
        target=fields["target"],
        batch=fields["batch"],
        samples=fields["samples"],
    )
    result_number = current.revision_number + 1
    revision = _ExperimentRevision(
        experiment_revision_id=revision_id,
        experiment_id=experiment.experiment_id,
        project_id=experiment.project_id,
        revision_number=result_number,
        canonical_payload=payload.canonical_json(),
        payload=payload,
        created_at=write.now,
    )
    experiment.revisions.append(revision)
    repository._experiment_revisions[revision_id] = revision
    return applied_target(
        item,
        base_revision_id=current.experiment_revision_id,
        base_revision_number=current.revision_number,
        result_revision_id=revision_id,
        result_revision_number=result_number,
    )


def _apply_sample(repository, write: CatalogWrite, item, revision_id: str):
    try:
        sample = repository._registry.get_sample(item.target_id)
        head_id, head_number = repository._registry.get_sample_revision_head(
            item.target_id
        )
        previous = repository._registry.get_sample_revision(head_id)
    except KeyError:
        raise catalog_not_found() from None
    decoded = json.loads(previous.canonical_payload)
    attributes = {
        key: value for key, value in decoded["attributes"].items() if key != "catalog"
    }
    display_name, catalog = overlay_sample_catalog(previous.canonical_payload, item)
    block = interpret_sample_catalog(
        catalog,
        _MemoryTermLookup(repository),
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
    result_number = head_number + 1
    revision = SampleRevision(
        sample_revision_id=revision_id,
        project_id=sample.project_id,
        sample_id=sample.sample_id,
        revision_number=result_number,
        canonical_payload=canonical,
        payload_digest_scheme=SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
        payload_digest=build_sample_revision_payload_digest(canonical),
        created_at=write.now,
    )
    try:
        repository._registry.append_sample_revision(
            revision,
            expected_previous_revision_number=head_number,
        )
    except DataRegistryConflictError:
        raise catalog_revision_conflict() from None
    repository._registry.mark_catalog_revision_trusted(revision_id)
    return applied_target(
        item,
        base_revision_id=head_id,
        base_revision_number=head_number,
        result_revision_id=revision_id,
        result_revision_number=result_number,
    )


def _heads(repository, proposal) -> dict[str, TargetHead]:
    return {
        item.target_id: _head(repository, item, proposal.project_id)
        for item in proposal.items
    }


def _head(repository, item, project_id: str) -> TargetHead:
    if item.target_kind == "experiment":
        experiment = repository._experiments.get(item.target_id)
        if experiment is None or experiment.project_id != project_id:
            raise catalog_not_found()
        current = experiment.revisions[-1]
        return TargetHead(
            target_kind="experiment",
            target_id=experiment.experiment_id,
            project_id=experiment.project_id,
            revision_number=current.revision_number,
            revision_id=current.experiment_revision_id,
            archived=experiment.archived_at is not None,
            canonical_payload=current.canonical_payload,
        )
    if item.target_kind == "sample_catalog":
        try:
            sample = repository._registry.get_sample(item.target_id)
            head_id, head_number = repository._registry.get_sample_revision_head(
                item.target_id
            )
            revision = repository._registry.get_sample_revision(head_id)
        except KeyError:
            raise catalog_not_found() from None
        if sample.project_id != project_id:
            raise catalog_not_found()
        return TargetHead(
            target_kind="sample_catalog",
            target_id=sample.sample_id,
            project_id=sample.project_id,
            revision_number=head_number,
            revision_id=head_id,
            archived=False,
            canonical_payload=revision.canonical_payload,
        )
    raise catalog_not_found()


def _load(repository, proposal_id: str):
    stored = repository._proposals.get(proposal_id)
    if stored is None:
        raise catalog_not_found()
    return proposal_from_document(stored["document"])


def _adoption(repository, proposal_id: str) -> dict[str, Any] | None:
    stored = repository._adoptions.get(proposal_id)
    if stored is None:
        return None
    return deepcopy(stored)
