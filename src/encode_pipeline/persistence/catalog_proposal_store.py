"""Session helpers for registering and adopting one offline catalog proposal."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from datetime import datetime
import json
from typing import Any

from sqlalchemy import select
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session

from encode_pipeline.persistence.experiment_catalog import (
    SessionTermLookup,
    _add_receipt,
    _experiment_payload,
    _insert_experiment_revision,
    _require_active_project,
    _sample_revision_row,
    _target_term_record,
    sync_sample_revision_state,
)
from encode_pipeline.persistence.models import (
    CatalogAdoptionReceiptRow,
    CatalogCurationVersionRow,
    CatalogWriteReceiptRow,
    CatalogProposalRow,
    CatalogTermRevisionRow,
    CatalogTermRow,
    ExperimentRevisionHeadRow,
    ExperimentRevisionRow,
    ExperimentRow,
    SampleRevisionHeadRow,
    SampleRevisionRow,
    SampleRow,
)
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
    catalog_idempotency_conflict,
    matches_committed_catalog_request,
    catalog_not_found,
    catalog_payload_invalid,
    catalog_revision_conflict,
    interpret_sample_catalog,
    sample_catalog_object,
    validate_display_name,
)


class AdoptionConflict(Exception):
    """A committed adoption receipt may already represent this attempt."""


_REVISION_INTEGRITY_TOKENS = (
    "experiment_revisions",
    "experiment_revision_heads",
    "experiment_revision_samples",
    "sample_revisions",
)
_RECEIPT_INTEGRITY_TOKENS = (
    "catalog_adoption_receipts",
    "catalog_write_receipts",
)


def read_active_versions(session: Session) -> ActiveCurationVersions:
    """Return the single active curation version row."""

    row = session.get(CatalogCurationVersionRow, 1)
    if row is None:
        raise CatalogError(
            409,
            "catalog_proposal_invalidated",
            "Active curation versions are not configured.",
        )
    return ActiveCurationVersions(
        rule_version=row.rule_version,
        vocabulary_version=row.vocabulary_version,
        projection_version=row.projection_version,
    )


def write_active_versions(session: Session, versions: ActiveCurationVersions) -> None:
    """Replace the active curation versions. This row is not an HTTP resource."""

    row = session.get(CatalogCurationVersionRow, 1)
    if row is None:
        session.add(
            CatalogCurationVersionRow(
                version_id=1,
                rule_version=versions.rule_version,
                vocabulary_version=versions.vocabulary_version,
                projection_version=versions.projection_version,
            )
        )
        return
    row.rule_version = versions.rule_version
    row.vocabulary_version = versions.vocabulary_version
    row.projection_version = versions.projection_version


def register_proposal(
    session: Session,
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
    """Insert one immutable proposal and its frozen registration view."""

    _require_active_project(session, project_id)
    active = read_active_versions(session)
    proposal = build_registered_proposal(
        proposal_id=proposal_id,
        project_id=project_id,
        created_by_user_id=write.actor_user_id,
        source=source,
        items=items,
        created_at=write.now,
        active=active,
        rule_version=rule_version,
        vocabulary_version=vocabulary_version,
        projection_version=projection_version,
    )
    lookup = SessionTermLookup(session)
    heads = _heads(session, proposal)
    model_source = proposal.source if proposal.source == "peri_bound_model" else None
    if model_source is not None:
        ensure_current_model_target(
            proposal,
            heads,
            active,
            lambda term_id: _target_term_record(session, term_id),
        )
    view = build_proposal_view(
        proposal,
        heads=heads,
        adopted_receipt=None,
        include_candidates=True,
        now=write.now,
        lookup=lookup,
    )
    session.add(
        CatalogProposalRow(
            proposal_id=proposal.proposal_id,
            project_id=proposal.project_id,
            created_by_user_id=proposal.created_by_user_id,
            source="agent" if model_source is not None else proposal.source,
            model_source=model_source,
            rule_version=proposal.rule_version,
            vocabulary_version=proposal.vocabulary_version,
            projection_version=proposal.projection_version,
            proposal_digest=proposal.proposal_digest,
            document_json=_dump(proposal_document(proposal)),
            response_json=_dump(view),
            register_actor_user_id=write.actor_user_id,
            register_idempotency_key=write.idempotency_key,
            created_at=proposal.created_at,
            expires_at=proposal.expires_at,
            retention_until=proposal.retention_until,
        )
    )
    session.flush()
    _add_receipt(
        session,
        write,
        subject_id=project_id,
        result_revision_id=None,
        result_archived_at=None,
    )
    return view


def get_proposal(
    session: Session,
    *,
    proposal_id: str,
    include_candidates: bool,
    now: datetime,
) -> dict[str, Any]:
    """Rebuild the live proposal view from current heads."""

    proposal = _load_proposal(session, proposal_id)
    _require_active_project(session, proposal.project_id)
    receipt = _adoption_body(session, proposal.proposal_id)
    return build_proposal_view(
        proposal,
        heads=_heads(session, proposal),
        adopted_receipt=receipt,
        include_candidates=include_candidates,
        now=now,
        lookup=SessionTermLookup(session),
    )


def preview_proposal(
    session: Session,
    *,
    proposal_id: str,
    selected_item_ids: object,
    now: datetime,
) -> dict[str, Any]:
    """Digest one subset without appending a catalog revision."""

    proposal = _load_proposal(session, proposal_id)
    _require_active_project(session, proposal.project_id)
    if _adoption_body(session, proposal.proposal_id) is not None:
        raise CatalogError(
            409,
            "catalog_proposal_invalidated",
            "This proposal has already been adopted.",
        )
    return build_preview(
        proposal,
        selected_item_ids=selected_item_ids,
        heads=_heads(session, proposal),
        active=read_active_versions(session),
        now=now,
        sources=_SqlSourceCatalog(session),
        lookup=SessionTermLookup(session),
    )


def adopt_proposal(
    session: Session,
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
    """Append the selected revisions and exactly one adoption receipt."""

    proposal = _load_proposal(session, proposal_id)
    _require_active_project(session, proposal.project_id)
    adopted = _adoption_body(session, proposal.proposal_id) is not None
    selected = prepare_adoption(
        proposal,
        selected_item_ids=selected_item_ids,
        expected_revisions=expected_revisions,
        preview_digest=preview_digest,
        heads=_heads(session, proposal),
        active=read_active_versions(session),
        now=write.now,
        adopted=adopted,
        sources=_SqlSourceCatalog(session),
        lookup=SessionTermLookup(session),
    )
    if not isinstance(preview_digest, str):
        raise catalog_payload_invalid()
    try:
        applied = [
            _apply_item(
                session,
                write,
                item,
                experiment_revision_id_factory=experiment_revision_id_factory,
                sample_revision_id_factory=sample_revision_id_factory,
            )
            for item in selected
        ]
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
        session.add(
            CatalogAdoptionReceiptRow(
                receipt_id=body["receipt_id"],
                proposal_id=proposal.proposal_id,
                project_id=proposal.project_id,
                actor_user_id=write.actor_user_id,
                idempotency_key=write.idempotency_key,
                response_json=_dump(body),
                created_at=write.now,
            )
        )
        session.flush()
        _add_receipt(
            session,
            write,
            subject_id=proposal.project_id,
            result_revision_id=None,
            result_archived_at=None,
        )
    except IntegrityError as exc:
        _raise_for_adoption_integrity(exc)
        raise
    return body


def list_adoption_receipts(
    session: Session,
    *,
    project_id: str,
    proposal_id: str | None,
) -> dict[str, Any]:
    """Return committed adoption receipts for one project."""

    _require_active_project(session, project_id)
    statement = select(CatalogAdoptionReceiptRow).where(
        CatalogAdoptionReceiptRow.project_id == project_id
    )
    if proposal_id is not None:
        statement = statement.where(
            CatalogAdoptionReceiptRow.proposal_id == proposal_id
        )
    rows = session.scalars(
        statement.order_by(
            CatalogAdoptionReceiptRow.created_at,
            CatalogAdoptionReceiptRow.receipt_id,
        )
    ).all()
    return {
        "ok": True,
        "project_id": project_id,
        "receipts": [_load(row.response_json) for row in rows],
    }


def replay_register(session: Session, actor_user_id: str, idempotency_key: str):
    """Return the registration view frozen when the proposal was inserted."""

    row = session.scalar(
        select(CatalogProposalRow).where(
            CatalogProposalRow.register_actor_user_id == actor_user_id,
            CatalogProposalRow.register_idempotency_key == idempotency_key,
        )
    )
    if row is None:
        raise catalog_idempotency_conflict()
    return _load(row.response_json)


def resolve_adoption_conflict(
    session: Session,
    write: CatalogWrite,
    proposal_id: str,
) -> dict[str, Any]:
    """Judge a rolled-back attempt against the committed catalog write receipt."""

    committed = session.get(
        CatalogWriteReceiptRow,
        (write.actor_user_id, write.idempotency_key),
    )
    if committed is not None:
        if not matches_committed_catalog_request(
            write,
            operation=committed.operation,
            target_key=committed.target_key,
            request_sha256=committed.request_sha256,
        ):
            raise catalog_idempotency_conflict()
        return replay_adopt(session, write.actor_user_id, write.idempotency_key)
    existing = _adoption_body(session, proposal_id)
    if existing is not None:
        raise CatalogError(
            409,
            "catalog_proposal_invalidated",
            "This proposal has already been adopted.",
        )
    raise catalog_revision_conflict()


def _reported_unique_columns(exc: IntegrityError) -> str:
    """Return the columns named by the database, without SQL or parameters."""

    orig = getattr(exc, "orig", None)
    args = getattr(orig, "args", None)
    if not args or not isinstance(args[0], str):
        return ""
    marker = "UNIQUE constraint failed:"
    if marker not in args[0]:
        return ""
    return args[0].split(marker, 1)[1].splitlines()[0]


def _raise_for_adoption_integrity(exc: IntegrityError) -> None:
    """Separate a head collision from the one committed adoption receipt."""

    reported = _reported_unique_columns(exc)
    if any(token in reported for token in _REVISION_INTEGRITY_TOKENS):
        raise catalog_revision_conflict() from exc
    if any(token in reported for token in _RECEIPT_INTEGRITY_TOKENS):
        raise AdoptionConflict() from exc
    raise exc


def replay_adopt(session: Session, actor_user_id: str, idempotency_key: str):
    """Return the original adoption body with the replay flag set."""

    row = session.scalar(
        select(CatalogAdoptionReceiptRow).where(
            CatalogAdoptionReceiptRow.actor_user_id == actor_user_id,
            CatalogAdoptionReceiptRow.idempotency_key == idempotency_key,
        )
    )
    if row is None:
        raise catalog_idempotency_conflict()
    return replayed_receipt(_load(row.response_json))


class _SqlSourceCatalog(SourceCatalog):
    """Prove citations against rows visible in the current transaction."""

    def __init__(self, session: Session) -> None:
        self._session = session

    def cited_experiment(self, experiment_id: str) -> tuple[str, bool] | None:
        row = self._session.get(ExperimentRow, experiment_id)
        if row is None:
            return None
        return row.project_id, row.archived_at is not None

    def cited_experiment_revision(self, revision_id: str) -> tuple[str, bool] | None:
        revision = self._session.get(ExperimentRevisionRow, revision_id)
        if revision is None:
            return None
        experiment = self._session.get(ExperimentRow, revision.experiment_id)
        if experiment is None:
            return None
        return experiment.project_id, experiment.archived_at is not None

    def cited_sample(self, sample_id: str) -> str | None:
        row = self._session.get(SampleRow, sample_id)
        if row is None:
            return None
        return row.project_id

    def cited_sample_revision(self, revision_id: str) -> str | None:
        row = self._session.get(SampleRevisionRow, revision_id)
        if row is None:
            return None
        return row.project_id

    def cited_term_archived(self, revision_id: str) -> bool | None:
        revision = self._session.get(CatalogTermRevisionRow, revision_id)
        if revision is None:
            return None
        term = self._session.get(CatalogTermRow, revision.term_id)
        if term is None:
            return None
        return term.archived_at is not None


def _apply_item(
    session: Session,
    write: CatalogWrite,
    item,
    *,
    experiment_revision_id_factory: Callable[[], str],
    sample_revision_id_factory: Callable[[], str],
) -> dict[str, Any]:
    if item.target_kind == "experiment":
        return _apply_experiment(
            session,
            write,
            item,
            experiment_revision_id_factory(),
        )
    if item.target_kind == "sample_catalog":
        return _apply_sample(
            session,
            write,
            item,
            sample_revision_id_factory(),
        )
    raise catalog_payload_invalid()


def _apply_experiment(session: Session, write: CatalogWrite, item, revision_id: str):
    experiment = session.get(ExperimentRow, item.target_id)
    head = session.get(ExperimentRevisionHeadRow, item.target_id)
    if experiment is None or head is None:
        raise catalog_not_found()
    previous = session.get(ExperimentRevisionRow, head.experiment_revision_id)
    if previous is None:
        raise catalog_not_found()
    fields = overlay_experiment_fields(previous.canonical_payload, item)
    payload = _experiment_payload(
        session,
        project_id=experiment.project_id,
        display_name=validate_display_name(fields["display_name"]),
        assay=fields["assay"],
        target=fields["target"],
        batch=fields["batch"],
        samples=fields["samples"],
    )
    base_number = head.revision_number
    base_revision_id = head.experiment_revision_id
    result_number = base_number + 1
    _insert_experiment_revision(
        session,
        experiment_id=experiment.experiment_id,
        project_id=experiment.project_id,
        experiment_revision_id=revision_id,
        revision_number=result_number,
        canonical=payload.canonical_json(),
        payload=payload,
        created_at=write.now,
        identity_created_at=experiment.created_at,
    )
    return applied_target(
        item,
        base_revision_id=base_revision_id,
        base_revision_number=base_number,
        result_revision_id=revision_id,
        result_revision_number=result_number,
    )


def _apply_sample(session: Session, write: CatalogWrite, item, revision_id: str):
    sample = session.get(SampleRow, item.target_id)
    head = session.get(SampleRevisionHeadRow, item.target_id)
    if sample is None or head is None:
        raise catalog_not_found()
    previous = session.get(SampleRevisionRow, head.sample_revision_id)
    if previous is None:
        raise catalog_not_found()
    decoded = json.loads(previous.canonical_payload)
    attributes = {
        key: value for key, value in decoded["attributes"].items() if key != "catalog"
    }
    display_name, catalog = overlay_sample_catalog(previous.canonical_payload, item)
    block = interpret_sample_catalog(
        catalog,
        SessionTermLookup(session),
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
    base_number = head.revision_number
    base_revision_id = head.sample_revision_id
    result_number = base_number + 1
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
    session.add(_sample_revision_row(revision))
    session.flush()
    sync_sample_revision_state(
        session,
        revision,
        allow_catalog=True,
        catalog_trusted=True,
    )
    return applied_target(
        item,
        base_revision_id=base_revision_id,
        base_revision_number=base_number,
        result_revision_id=revision_id,
        result_revision_number=result_number,
    )


def _heads(session: Session, proposal) -> dict[str, TargetHead]:
    return {
        item.target_id: _head(session, item, proposal.project_id)
        for item in proposal.items
    }


def _head(session: Session, item, project_id: str) -> TargetHead:
    if item.target_kind == "experiment":
        experiment = session.get(ExperimentRow, item.target_id)
        pointer = session.get(ExperimentRevisionHeadRow, item.target_id)
        if experiment is None or pointer is None or experiment.project_id != project_id:
            raise catalog_not_found()
        revision = session.get(ExperimentRevisionRow, pointer.experiment_revision_id)
        if revision is None:
            raise catalog_not_found()
        return TargetHead(
            target_kind="experiment",
            target_id=experiment.experiment_id,
            project_id=experiment.project_id,
            revision_number=pointer.revision_number,
            revision_id=pointer.experiment_revision_id,
            archived=experiment.archived_at is not None,
            canonical_payload=revision.canonical_payload,
        )
    if item.target_kind == "sample_catalog":
        sample = session.get(SampleRow, item.target_id)
        pointer = session.get(SampleRevisionHeadRow, item.target_id)
        if sample is None or pointer is None or sample.project_id != project_id:
            raise catalog_not_found()
        revision = session.get(SampleRevisionRow, pointer.sample_revision_id)
        if revision is None:
            raise catalog_not_found()
        return TargetHead(
            target_kind="sample_catalog",
            target_id=sample.sample_id,
            project_id=sample.project_id,
            revision_number=pointer.revision_number,
            revision_id=pointer.sample_revision_id,
            archived=False,
            canonical_payload=revision.canonical_payload,
        )
    raise catalog_not_found()


def _load_proposal(session: Session, proposal_id: str):
    row = session.get(CatalogProposalRow, proposal_id)
    if row is None:
        raise catalog_not_found()
    proposal = proposal_from_document(_load(row.document_json))
    if proposal.proposal_id != row.proposal_id:
        raise catalog_payload_invalid()
    return proposal


def _adoption_body(session: Session, proposal_id: str) -> dict[str, Any] | None:
    row = session.scalar(
        select(CatalogAdoptionReceiptRow).where(
            CatalogAdoptionReceiptRow.proposal_id == proposal_id
        )
    )
    if row is None:
        return None
    return _load(row.response_json)


def _dump(value: Mapping[str, Any]) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _load(value: str) -> dict[str, Any]:
    loaded = json.loads(value)
    if not isinstance(loaded, dict):
        raise catalog_payload_invalid()
    return loaded
