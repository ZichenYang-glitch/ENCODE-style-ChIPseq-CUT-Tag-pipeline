"""Persist and adopt an offline catalog proposal on the Slice A write path.

The preview digest, citation shape, and error codes come from the reviewed
IRENE-27 contract. This module does not open a database session.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
import json
from typing import Any, Protocol

from encode_pipeline.platform.catalog_curation_preview import (
    MAX_CHANGESET_ITEMS,
    MAX_PREVIEW_APPLY_BODY_BYTES,
    M2_CONTRACT_VERSION,
    M2_PROJECTION_VERSION,
    M2_RULE_VERSION,
    M2_VOCABULARY_VERSION,
    CatalogPreviewError,
    CatalogProposal,
    ProposalCitation,
    ProposalItem,
    SourceReference,
    TermRevisionRecord,
    _REPLICATE_GROUP_RE,
    _SAMPLE_PIN_KEYS,
    _check_id,
    _validate_and_render_slot,
    compute_preview_digest,
)
from encode_pipeline.platform.experiment_catalog import (
    CatalogError,
    TermLookup,
    catalog_experiment_archived,
    catalog_not_found,
    catalog_payload_invalid,
    catalog_revision_conflict,
    catalog_term_archived,
    interpret_experiment_fields,
    interpret_sample_catalog,
    validate_display_name,
)


@dataclass(frozen=True)
class ActiveCurationVersions:
    """Rule, vocabulary, and projection versions a proposal must still match."""

    rule_version: str = M2_RULE_VERSION
    vocabulary_version: str = M2_VOCABULARY_VERSION
    projection_version: str = M2_PROJECTION_VERSION


DEFAULT_CURATION_VERSIONS = ActiveCurationVersions()


@dataclass(frozen=True)
class TargetHead:
    """Current catalog head used to check one proposal target."""

    target_kind: str
    target_id: str
    project_id: str
    revision_number: int
    revision_id: str
    archived: bool
    canonical_payload: str


class SourceCatalog(Protocol):
    """Lookups that prove references still belong to the proposal's project."""

    def cited_experiment(self, experiment_id: str) -> tuple[str, bool] | None:
        """Return project id and archived flag, or None when the id is absent."""

    def cited_experiment_revision(self, revision_id: str) -> tuple[str, bool] | None:
        """Return the owning experiment's project id and archived flag."""

    def cited_sample(self, sample_id: str) -> str | None:
        """Return the sample's project id, or None when it is absent."""

    def cited_sample_revision(self, revision_id: str) -> str | None:
        """Return the owning sample's project id, or None when it is absent."""

    def cited_term_archived(self, revision_id: str) -> bool | None:
        """Return whether the term is archived, or None when the revision is absent."""


def stable_proposal_id(
    actor_user_id: str,
    idempotency_key: str,
    project_id: str,
) -> str:
    """Derive one proposal id so a retry binds the same Slice A receipt target."""

    digest = sha256(
        (f"catalog-proposal\n{actor_user_id}\n{idempotency_key}\n{project_id}").encode()
    ).hexdigest()
    return f"cprop_{digest[:32]}"


def assert_request_within_limit(raw_body: bytes) -> None:
    """Reject a catalog write whose raw body exceeds the reviewed 256 KiB cap."""

    if len(raw_body) > MAX_PREVIEW_APPLY_BODY_BYTES:
        raise CatalogError(
            413,
            "catalog_payload_invalid",
            "Request body exceeds the 256 KiB catalog changeset limit.",
        )


def format_utc(value: datetime) -> str:
    """Format a UTC timestamp the same way the reviewed preview contract does."""

    if value.tzinfo is None or value.utcoffset() is None:
        raise CatalogError(
            422,
            "catalog_payload_invalid",
            "The catalog payload is invalid.",
        )
    return (
        value.astimezone(timezone.utc)
        .isoformat(timespec="seconds")
        .replace("+00:00", "Z")
    )


def parse_utc(value: str) -> datetime:
    """Parse a contract timestamp that ends in Z."""

    if not isinstance(value, str) or not value.endswith("Z"):
        raise catalog_payload_invalid()
    parsed = datetime.fromisoformat(value[:-1] + "+00:00")
    if parsed.tzinfo is None:
        raise catalog_payload_invalid()
    return parsed.astimezone(timezone.utc)


def truncate_seconds(value: datetime) -> datetime:
    """Drop microseconds so stored stamps round-trip through format_utc."""

    if value.tzinfo is None or value.utcoffset() is None:
        raise catalog_payload_invalid()
    return value.astimezone(timezone.utc).replace(microsecond=0)


def build_registered_proposal(
    *,
    proposal_id: str,
    project_id: str,
    created_by_user_id: str,
    source: str,
    items: object,
    created_at: datetime,
    active: ActiveCurationVersions,
    rule_version: str | None,
    vocabulary_version: str | None,
    projection_version: str | None,
) -> CatalogProposal:
    """Validate one offline suggestion and bind it to the active versions by default."""

    if not isinstance(items, list) or not items:
        raise catalog_payload_invalid()
    _require_version(rule_version or active.rule_version)
    _require_version(vocabulary_version or active.vocabulary_version)
    _require_version(projection_version or active.projection_version)
    parsed = tuple(_item_from_object(item) for item in items)
    return _guard(
        lambda: CatalogProposal.create(
            proposal_id=proposal_id,
            project_id=project_id,
            created_by_user_id=created_by_user_id,
            source=source,
            items=parsed,
            created_at=truncate_seconds(created_at),
            rule_version=rule_version or active.rule_version,
            vocabulary_version=vocabulary_version or active.vocabulary_version,
            projection_version=projection_version or active.projection_version,
        )
    )


def proposal_document(proposal: CatalogProposal) -> dict[str, Any]:
    """Return the immutable JSON stored for one proposal."""

    return {
        "proposal_id": proposal.proposal_id,
        "project_id": proposal.project_id,
        "created_by_user_id": proposal.created_by_user_id,
        "source": proposal.source,
        "rule_version": proposal.rule_version,
        "vocabulary_version": proposal.vocabulary_version,
        "projection_version": proposal.projection_version,
        "proposal_digest": proposal.proposal_digest,
        "created_at": format_utc(proposal.created_at),
        "items": [_item_document(item) for item in proposal.items],
    }


def proposal_from_document(document: Mapping[str, Any]) -> CatalogProposal:
    """Rebuild a proposal and reject a document whose digest no longer matches."""

    items = tuple(_item_from_object(item) for item in document.get("items", ()))
    proposal = _guard(
        lambda: CatalogProposal.create(
            proposal_id=str(document.get("proposal_id")),
            project_id=str(document.get("project_id")),
            created_by_user_id=str(document.get("created_by_user_id")),
            source=str(document.get("source")),
            items=items,
            created_at=parse_utc(str(document.get("created_at"))),
            rule_version=str(document.get("rule_version")),
            vocabulary_version=str(document.get("vocabulary_version")),
            projection_version=str(document.get("projection_version")),
        )
    )
    if proposal.proposal_digest != document.get("proposal_digest"):
        raise catalog_payload_invalid()
    return proposal


def resolve_selected_items(
    proposal: CatalogProposal,
    selected_item_ids: object,
) -> tuple[ProposalItem, ...]:
    """Keep the proposal's item order and reject an empty or unknown selection."""

    if (
        not isinstance(selected_item_ids, Sequence)
        or isinstance(selected_item_ids, (str, bytes))
        or not selected_item_ids
    ):
        raise CatalogError(
            422,
            "catalog_payload_invalid",
            "At least one proposal item must be selected for preview or adoption.",
        )
    if len(selected_item_ids) > MAX_CHANGESET_ITEMS:
        raise CatalogError(
            422,
            "catalog_payload_invalid",
            f"Selection exceeds maximum of {MAX_CHANGESET_ITEMS} items.",
        )
    if any(not isinstance(item_id, str) for item_id in selected_item_ids):
        raise catalog_payload_invalid()
    if len(set(selected_item_ids)) != len(selected_item_ids):
        raise CatalogError(
            422,
            "catalog_payload_invalid",
            "Duplicate item_id values in selection.",
        )
    by_id = {item.item_id: item for item in proposal.items}
    unknown = [item_id for item_id in selected_item_ids if item_id not in by_id]
    if unknown:
        raise CatalogError(
            422,
            "catalog_payload_invalid",
            "Selection references an unknown proposal item_id.",
            path="selected_item_ids",
        )
    selected = set(selected_item_ids)
    return tuple(item for item in proposal.items if item.item_id in selected)


def proposal_state(
    proposal: CatalogProposal,
    *,
    adopted: bool,
    now: datetime,
) -> str:
    """Return pending, adopted, or expired. Adopted receipts outlive the window."""

    if adopted:
        return "adopted"
    if now >= proposal.expires_at:
        return "expired"
    return "pending"


def ensure_readable(proposal: CatalogProposal, *, adopted: bool, now: datetime) -> str:
    """Hide proposals past the seven-day retention window unless they were adopted."""

    state = proposal_state(proposal, adopted=adopted, now=now)
    if not adopted and now >= proposal.retention_until:
        raise catalog_not_found()
    return state


def ensure_adoptable(
    proposal: CatalogProposal, *, adopted: bool, now: datetime
) -> None:
    """Reject a second adoption and an expired window before any revision is written."""

    if adopted:
        raise CatalogError(
            409,
            "catalog_proposal_invalidated",
            "This proposal has already been adopted.",
            hint=(
                "Reopen the proposal or query adoption receipts "
                "to view the committed result."
            ),
        )
    if now >= proposal.retention_until:
        raise catalog_not_found()
    if now >= proposal.expires_at:
        raise CatalogError(
            410,
            "catalog_proposal_expired",
            "The proposal adoption window has expired.",
            hint=(
                "Read-only inspection remains available "
                "within the 7-day retention window."
            ),
        )


def ensure_active_versions(
    proposal: CatalogProposal,
    active: ActiveCurationVersions,
) -> None:
    """Reject a proposal whose rule, vocabulary, or citation version was retired."""

    if (
        proposal.rule_version != active.rule_version
        or proposal.vocabulary_version != active.vocabulary_version
        or proposal.projection_version != active.projection_version
    ):
        raise CatalogError(
            409,
            "catalog_proposal_invalidated",
            "Proposal rule, vocabulary, or projection version is no longer active.",
        )
    for item in proposal.items:
        citation = item.citation
        if (
            citation.rule_version != active.rule_version
            or citation.vocabulary_version != active.vocabulary_version
        ):
            raise CatalogError(
                409,
                "catalog_proposal_invalidated",
                "Item citation rule or vocabulary version is no longer active.",
                path=item.item_id,
            )


def verify_item_references(
    item: ProposalItem,
    project_id: str,
    catalog: SourceCatalog,
) -> None:
    """Validate cited sources and sample-revision pins for one selected item."""

    path = f"{item.item_id}.citations.source_refs"
    for ref in item.citation.source_refs:
        if ref.kind == "experiment_revision":
            _require_experiment_owner(
                catalog.cited_experiment_revision(ref.ref_id),
                project_id,
                path=path,
                message="Cited experiment revision is missing or outside this project.",
            )
        elif ref.kind == "sample_revision":
            _require_sample_owner(
                catalog.cited_sample_revision(ref.ref_id),
                project_id,
                path=path,
                message="Cited sample revision is missing or outside this project.",
            )
        elif ref.kind == "term_revision":
            archived = catalog.cited_term_archived(ref.ref_id)
            if archived is None:
                raise CatalogError(
                    409,
                    "catalog_proposal_invalidated",
                    "Cited term revision does not exist in the catalog.",
                    path=path,
                )
            if archived:
                raise CatalogError(
                    409,
                    "catalog_term_archived",
                    "Cited term revision is archived.",
                    path=path,
                )
        elif ref.kind == "catalog_projection_row":
            _verify_projection_ref(ref.ref_id, project_id, catalog, path=path)
        else:
            raise catalog_payload_invalid()

    if "samples" in item.allowed_fields:
        for side, fields in (("before", item.before), ("after", item.after)):
            for pin in fields["samples"]:
                _require_sample_owner(
                    catalog.cited_sample_revision(pin["sample_revision_id"]),
                    project_id,
                    path=f"{item.item_id}.{side}.samples.{pin['ordinal']}.sample_revision_id",
                    message="Sample-revision pin is missing or outside this project.",
                )


def render_item_fields(
    item: ProposalItem,
    lookup: TermLookup,
    *,
    include_candidates: bool,
    for_write: bool,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Render before/after slots, including display names for pinned terms."""

    records = _term_records(item, lookup)
    before: dict[str, Any] = {}
    after: dict[str, Any] = {}
    for field_name in item.allowed_fields:
        if field_name == "display_name":
            before[field_name] = item.before[field_name]
            after[field_name] = item.after[field_name]
            continue
        if field_name == "samples":
            before[field_name] = _render_sample_pins(
                item.before[field_name],
                records,
                include_candidates=include_candidates,
                for_write=False,
            )
            after[field_name] = _render_sample_pins(
                item.after[field_name],
                records,
                include_candidates=include_candidates,
                for_write=for_write,
            )
            continue
        before[field_name] = _render_slot(
            field_name,
            item.before[field_name],
            records,
            include_candidates=include_candidates,
            for_write=False,
        )
        after[field_name] = _render_slot(
            field_name,
            item.after[field_name],
            records,
            include_candidates=include_candidates,
            for_write=for_write,
        )
    return before, after


def project_fields_for_role(
    fields: Mapping[str, Any],
    *,
    include_candidates: bool,
) -> dict[str, Any]:
    """Drop ambiguous candidate pins from a member view."""

    if include_candidates:
        return dict(fields)
    projected: dict[str, Any] = {}
    for name, slot in fields.items():
        if isinstance(slot, Mapping) and "candidates" in slot:
            projected[name] = {
                key: value for key, value in slot.items() if key != "candidates"
            }
        else:
            projected[name] = slot
    return projected


def build_proposal_view(
    proposal: CatalogProposal,
    *,
    heads: Mapping[str, TargetHead],
    adopted_receipt: Mapping[str, Any] | None,
    include_candidates: bool,
    now: datetime,
    lookup: TermLookup,
) -> dict[str, Any]:
    """Return the role-projected proposal, including a prior receipt when adopted."""

    adopted = adopted_receipt is not None
    state = ensure_readable(proposal, adopted=adopted, now=now)
    items: list[dict[str, Any]] = []
    for item in proposal.items:
        head = heads.get(item.target_id)
        if head is None or head.project_id != proposal.project_id:
            raise catalog_not_found()
        before, after = render_item_fields(
            item,
            lookup,
            include_candidates=include_candidates,
            for_write=False,
        )
        items.append(
            {
                "item_id": item.item_id,
                "target_kind": item.target_kind,
                "target_id": item.target_id,
                "base_revision_id": item.base_revision_id,
                "base_revision_number": item.base_revision_number,
                "current_revision_id": head.revision_id,
                "current_revision_number": head.revision_number,
                "allowed_fields": list(item.allowed_fields),
                "before": project_fields_for_role(
                    before, include_candidates=include_candidates
                ),
                "after": project_fields_for_role(
                    after, include_candidates=include_candidates
                ),
                "citations": item.citation.to_object(
                    include_candidates=include_candidates,
                    candidate_term_revision_ids=item.candidate_term_revision_ids(),
                ),
            }
        )
    return {
        "ok": True,
        "contract_version": M2_CONTRACT_VERSION,
        "proposal_id": proposal.proposal_id,
        "project_id": proposal.project_id,
        "source": proposal.source,
        "proposal_state": state,
        "adoptable_by_caller": include_candidates and state == "pending",
        "rule_version": proposal.rule_version,
        "vocabulary_version": proposal.vocabulary_version,
        "projection_version": proposal.projection_version,
        "proposal_digest": proposal.proposal_digest,
        "created_at": format_utc(proposal.created_at),
        "expires_at": format_utc(proposal.expires_at),
        "retention_until": format_utc(proposal.retention_until),
        "adoption_receipt": None if adopted_receipt is None else dict(adopted_receipt),
        "items": items,
    }


def build_preview(
    proposal: CatalogProposal,
    *,
    selected_item_ids: Sequence[str],
    heads: Mapping[str, TargetHead],
    active: ActiveCurationVersions,
    now: datetime,
    sources: SourceCatalog,
    lookup: TermLookup,
) -> dict[str, Any]:
    """Build an immutable subset digest without writing a catalog revision."""

    ensure_adoptable(proposal, adopted=False, now=now)
    ensure_active_versions(proposal, active)
    selected = resolve_selected_items(proposal, selected_item_ids)
    canonical_ids = [item.item_id for item in selected]
    excluded = [
        item.item_id for item in proposal.items if item.item_id not in canonical_ids
    ]
    expected: dict[str, int] = {}
    preview_items: list[dict[str, Any]] = []
    conflicts: list[dict[str, Any]] = []
    for item in selected:
        head = _require_head(heads, item, proposal.project_id)
        if head.archived:
            raise catalog_experiment_archived()
        verify_item_references(item, proposal.project_id, sources)
        before, after = render_item_fields(
            item,
            lookup,
            include_candidates=True,
            for_write=True,
        )
        _assert_writable(item, head, lookup)
        expected[item.target_id] = item.base_revision_number
        conflict = _conflict(item, head)
        if conflict is not None:
            conflicts.append(conflict)
        preview_items.append(
            {
                "item_id": item.item_id,
                "target_kind": item.target_kind,
                "target_id": item.target_id,
                "expected_revision_number": item.base_revision_number,
                "base_revision_id": item.base_revision_id,
                "current_revision_number": head.revision_number,
                "current_revision_id": head.revision_id,
                "conflict_detected": conflict is not None,
                "allowed_fields": list(item.allowed_fields),
                "before": before,
                "after": after,
                "citations": item.citation.to_object(include_candidates=True),
            }
        )
    digest = compute_preview_digest(
        proposal_id=proposal.proposal_id,
        proposal_digest=proposal.proposal_digest,
        selected_item_ids=canonical_ids,
        expected_revisions=expected,
    )
    return {
        "ok": True,
        "contract_version": M2_CONTRACT_VERSION,
        "proposal_id": proposal.proposal_id,
        "project_id": proposal.project_id,
        "proposal_digest": proposal.proposal_digest,
        "preview_digest": digest,
        "selected_item_ids": canonical_ids,
        "excluded_item_ids": excluded,
        "expected_revisions": expected,
        "has_conflicts": bool(conflicts),
        "conflicts": conflicts,
        "expires_at": format_utc(proposal.expires_at),
        "items": preview_items,
    }


def prepare_adoption(
    proposal: CatalogProposal,
    *,
    selected_item_ids: object,
    expected_revisions: object,
    preview_digest: object,
    heads: Mapping[str, TargetHead],
    active: ActiveCurationVersions,
    now: datetime,
    adopted: bool,
    sources: SourceCatalog,
    lookup: TermLookup,
) -> tuple[ProposalItem, ...]:
    """Check the subset, digest, versions, sources, and heads before any write."""

    if not isinstance(preview_digest, str) or len(preview_digest) != 64:
        raise catalog_payload_invalid()
    if any(character not in "0123456789abcdef" for character in preview_digest):
        raise catalog_payload_invalid()
    ensure_adoptable(proposal, adopted=adopted, now=now)
    ensure_active_versions(proposal, active)
    selected = resolve_selected_items(proposal, selected_item_ids)
    canonical_ids = [item.item_id for item in selected]
    revisions = _expected_revision_map(expected_revisions)
    digest_ids: Sequence[str] = (
        selected_item_ids
        if isinstance(selected_item_ids, list)
        and list(selected_item_ids) != canonical_ids
        else canonical_ids
    )
    recomputed = compute_preview_digest(
        proposal_id=proposal.proposal_id,
        proposal_digest=proposal.proposal_digest,
        selected_item_ids=digest_ids,
        expected_revisions=revisions,
    )
    if recomputed != preview_digest:
        raise CatalogError(
            409,
            "catalog_preview_digest_mismatch",
            "Submitted selection or expected revisions do not match preview_digest.",
        )
    if set(revisions) != {item.target_id for item in selected}:
        raise CatalogError(
            422,
            "catalog_payload_invalid",
            "expected_revisions must match the selected proposal targets.",
        )
    conflicts: list[dict[str, Any]] = []
    for item in selected:
        caller_expected = revisions[item.target_id]
        if caller_expected != item.base_revision_number:
            raise CatalogError(
                409,
                "catalog_proposal_invalidated",
                "Caller expected revision does not match the server proposal baseline.",
                path=item.item_id,
            )
        head = _require_head(heads, item, proposal.project_id)
        if head.archived:
            raise catalog_experiment_archived()
        verify_item_references(item, proposal.project_id, sources)
        _assert_writable(item, head, lookup)
        conflict = _conflict(item, head)
        if conflict is not None:
            conflicts.append(conflict)
    if conflicts:
        raise catalog_revision_conflict(conflicts)
    return selected


def adoption_receipt(
    proposal: CatalogProposal,
    *,
    selected: Sequence[ProposalItem],
    preview_digest: str,
    actor_user_id: str,
    idempotency_key: str,
    receipt_id: str,
    applied_targets: Sequence[Mapping[str, Any]],
    now: datetime,
) -> dict[str, Any]:
    """Build the one adoption receipt stored for this proposal."""

    selected_ids = [item.item_id for item in selected]
    excluded = [
        item.item_id for item in proposal.items if item.item_id not in selected_ids
    ]
    return {
        "ok": True,
        "replayed": False,
        "contract_version": M2_CONTRACT_VERSION,
        "receipt_id": receipt_id,
        "proposal_id": proposal.proposal_id,
        "project_id": proposal.project_id,
        "proposal_digest": proposal.proposal_digest,
        "preview_digest": preview_digest,
        "actor_user_id": actor_user_id,
        "idempotency_key": idempotency_key,
        "selected_item_ids": selected_ids,
        "excluded_item_ids": excluded,
        "applied_targets": [dict(item) for item in applied_targets],
        "adopted_at": format_utc(truncate_seconds(now)),
    }


def replayed_receipt(receipt: Mapping[str, Any]) -> dict[str, Any]:
    """Return the original adoption body with the replay flag set."""

    body = dict(receipt)
    body["replayed"] = True
    return body


def overlay_experiment_fields(
    canonical_payload: str,
    item: ProposalItem,
) -> dict[str, Any]:
    """Copy the stored experiment and replace only the proposal's allowed fields."""

    try:
        decoded = json.loads(canonical_payload)
    except json.JSONDecodeError:
        raise catalog_payload_invalid() from None
    if not isinstance(decoded, dict):
        raise catalog_payload_invalid()
    display_name = decoded.get("display_name")
    assay = decoded.get("assay")
    target = decoded.get("target")
    batch = decoded.get("batch")
    samples = _stored_sample_requests(decoded.get("samples"))
    for field_name in item.allowed_fields:
        if field_name == "display_name":
            display_name = item.after[field_name]
        elif field_name == "assay":
            assay = item.after[field_name]
        elif field_name == "target":
            target = item.after[field_name]
        elif field_name == "batch":
            batch = item.after[field_name]
        elif field_name == "samples":
            # The write parser assigns ordinals from the validated proposal order.
            samples = _stored_sample_requests(item.after[field_name])
        else:
            raise catalog_payload_invalid()
    return {
        "display_name": display_name,
        "assay": assay,
        "target": target,
        "batch": batch,
        "samples": samples,
    }


def overlay_sample_catalog(
    canonical_payload: str,
    item: ProposalItem,
) -> tuple[str | None, dict[str, Any]]:
    """Copy non-catalog attributes implicitly by returning only the new catalog block."""

    try:
        decoded = json.loads(canonical_payload)
    except json.JSONDecodeError:
        raise catalog_payload_invalid() from None
    attributes = decoded.get("attributes") if isinstance(decoded, dict) else None
    if not isinstance(attributes, dict):
        raise catalog_payload_invalid()
    catalog = attributes.get("catalog")
    if not isinstance(catalog, dict):
        raise catalog_payload_invalid()
    merged = dict(catalog)
    display_name: str | None = None
    for field_name in item.allowed_fields:
        if field_name == "display_name":
            display_name = item.after[field_name]
        elif field_name in merged:
            merged[field_name] = item.after[field_name]
        else:
            raise catalog_payload_invalid()
    return display_name, merged


def applied_target(
    item: ProposalItem,
    *,
    base_revision_id: str,
    base_revision_number: int,
    result_revision_id: str,
    result_revision_number: int,
) -> dict[str, Any]:
    """Describe one revision appended for a selected proposal item."""

    return {
        "item_id": item.item_id,
        "target_kind": item.target_kind,
        "target_id": item.target_id,
        "base_revision_id": base_revision_id,
        "base_revision_number": base_revision_number,
        "result_revision_id": result_revision_id,
        "result_revision_number": result_revision_number,
        "changed_fields": list(item.allowed_fields),
    }


def _require_version(value: str) -> None:
    if (
        not isinstance(value, str)
        or not 1 <= len(value) <= 64
        or any(character.isspace() for character in value)
    ):
        raise catalog_payload_invalid()


def _item_document(item: ProposalItem) -> dict[str, Any]:
    body = item.canonical_object()
    refs: list[dict[str, Any]] = []
    for ref in item.citation.source_refs:
        encoded = ref.to_object()
        if ref.visibility != "shared":
            encoded["visibility"] = ref.visibility
        refs.append(encoded)
    body["citations"]["source_refs"] = refs
    return body


def _item_from_object(raw: object) -> ProposalItem:
    if not isinstance(raw, Mapping):
        raise catalog_payload_invalid()
    citation_raw = raw.get("citations")
    if not isinstance(citation_raw, Mapping):
        raise catalog_payload_invalid()
    refs_raw = citation_raw.get("source_refs")
    codes = citation_raw.get("reason_codes")
    if not isinstance(refs_raw, list) or not isinstance(codes, list):
        raise catalog_payload_invalid()
    allowed = raw.get("allowed_fields")
    if not isinstance(allowed, list):
        raise catalog_payload_invalid()
    item = _guard(
        lambda: ProposalItem(
            item_id=raw.get("item_id"),
            target_kind=raw.get("target_kind"),
            target_id=raw.get("target_id"),
            base_revision_id=raw.get("base_revision_id"),
            base_revision_number=raw.get("base_revision_number"),
            allowed_fields=tuple(allowed),
            before=raw.get("before"),
            after=raw.get("after"),
            citation=ProposalCitation(
                rule_version=citation_raw.get("rule_version"),
                vocabulary_version=citation_raw.get("vocabulary_version"),
                reason_codes=tuple(codes),
                source_refs=tuple(_source_ref(ref) for ref in refs_raw),
            ),
        )
    )
    if "samples" in item.allowed_fields:
        _validate_proposal_sample_pins(item.before["samples"])
        _validate_proposal_sample_pins(item.after["samples"])
    return item


def _validate_proposal_sample_pins(pins: list[object]) -> None:
    if not pins:
        raise catalog_payload_invalid()
    revision_ids: set[str] = set()
    for ordinal, pin in enumerate(pins):
        if not isinstance(pin, Mapping):
            raise catalog_payload_invalid()
        allowed = _SAMPLE_PIN_KEYS
        if "replicate_group" in pin:
            group = pin["replicate_group"]
            if not isinstance(group, str) or not _REPLICATE_GROUP_RE.fullmatch(group):
                raise catalog_payload_invalid()
            allowed = allowed | {"replicate_group"}
        if set(pin) != allowed:
            raise catalog_payload_invalid()
        position = pin["ordinal"]
        if (
            isinstance(position, bool)
            or not isinstance(position, int)
            or position != ordinal
        ):
            raise CatalogError(
                422,
                "catalog_payload_invalid",
                "Sample-revision pin ordinal must match its list position.",
                path="samples",
            )
        revision_id = _guard(
            lambda: _check_id(
                "sample_revision", pin["sample_revision_id"], path="samples"
            )
        )
        if revision_id in revision_ids:
            raise catalog_payload_invalid()
        revision_ids.add(revision_id)


def _source_ref(raw: object) -> SourceReference:
    if not isinstance(raw, Mapping):
        raise catalog_payload_invalid()
    return SourceReference(
        kind=raw.get("kind"),
        ref_id=raw.get("ref_id"),
        term_revision_id=raw.get("term_revision_id"),
        evidence_token=raw.get("evidence_token"),
        visibility=raw.get("visibility", "shared"),
    )


def _guard(callback: Callable[[], Any]) -> Any:
    try:
        return callback()
    except CatalogPreviewError as exc:
        raise _from_preview(exc) from exc


def _from_preview(exc: CatalogPreviewError) -> CatalogError:
    conflicts = [item.to_object() for item in exc.conflicts] or None
    code = exc.code
    if code == "catalog_term_archived":
        return catalog_term_archived()
    if code == "catalog_experiment_archived":
        return catalog_experiment_archived()
    return CatalogError(
        exc.status_code,
        code,
        exc.message,
        path=exc.path,
        hint=exc.hint,
        conflicts=conflicts,
    )


def _expected_revision_map(value: object) -> dict[str, int]:
    if not isinstance(value, Mapping) or not value:
        raise catalog_payload_invalid()
    parsed: dict[str, int] = {}
    for key, revision in value.items():
        if not isinstance(key, str):
            raise catalog_payload_invalid()
        if isinstance(revision, bool) or not isinstance(revision, int):
            raise CatalogError(
                422,
                "catalog_payload_invalid",
                "expected_revision_number must be an integer.",
            )
        parsed[key] = revision
    return parsed


def _require_head(
    heads: Mapping[str, TargetHead],
    item: ProposalItem,
    project_id: str,
) -> TargetHead:
    head = heads.get(item.target_id)
    if (
        head is None
        or head.project_id != project_id
        or head.target_kind != item.target_kind
    ):
        raise catalog_not_found()
    return head


def _conflict(item: ProposalItem, head: TargetHead) -> dict[str, Any] | None:
    if (
        head.revision_number == item.base_revision_number
        and head.revision_id == item.base_revision_id
    ):
        return None
    return {
        "target_kind": item.target_kind,
        "target_id": item.target_id,
        "expected_revision_number": item.base_revision_number,
        "current_revision_number": head.revision_number,
        "current_revision_id": head.revision_id,
        "expected_revision_id": item.base_revision_id,
    }


def _assert_writable(item: ProposalItem, head: TargetHead, lookup: TermLookup) -> None:
    if item.target_kind == "experiment":
        fields = overlay_experiment_fields(head.canonical_payload, item)
        validate_display_name(fields["display_name"])
        interpret_experiment_fields(
            fields["assay"],
            fields["target"],
            fields["batch"],
            fields["samples"],
            lookup,
            for_write=True,
        )
        return
    display_name, catalog = overlay_sample_catalog(head.canonical_payload, item)
    interpret_sample_catalog(catalog, lookup, for_write=True)
    if display_name is not None:
        validate_display_name(display_name)


def _stored_sample_requests(samples: object) -> list[dict[str, Any]]:
    if not isinstance(samples, list) or not samples:
        raise catalog_payload_invalid()
    requests: list[dict[str, Any]] = []
    for item in samples:
        if not isinstance(item, Mapping):
            raise catalog_payload_invalid()
        request: dict[str, Any] = {
            "sample_revision_id": item.get("sample_revision_id"),
            "replicate_role": item.get("replicate_role"),
        }
        if item.get("replicate_group"):
            request["replicate_group"] = item["replicate_group"]
        requests.append(request)
    return requests


def _term_records(
    item: ProposalItem,
    lookup: TermLookup,
) -> dict[str, TermRevisionRecord]:
    revision_ids: set[str] = set()
    for slot_map in (item.before, item.after):
        revision_ids.update(_revision_ids(slot_map))
    records: dict[str, TermRevisionRecord] = {}
    for revision_id in revision_ids:
        pin = lookup(revision_id)
        if pin is None:
            continue
        records[revision_id] = TermRevisionRecord(
            term_id=pin.term_id,
            term_revision_id=pin.term_revision_id,
            term_type=pin.term_type,
            display_name=pin.display_name,
            archived=pin.archived,
        )
    return records


def _revision_ids(value: object) -> set[str]:
    found: set[str] = set()
    if isinstance(value, Mapping):
        revision_id = value.get("term_revision_id")
        if isinstance(revision_id, str):
            found.add(revision_id)
        for child in value.values():
            found.update(_revision_ids(child))
    elif isinstance(value, list):
        for child in value:
            found.update(_revision_ids(child))
    return found


def _render_slot(
    field_name: str,
    slot: object,
    records: Mapping[str, TermRevisionRecord],
    *,
    include_candidates: bool,
    for_write: bool,
) -> dict[str, Any]:
    if not isinstance(slot, Mapping):
        raise catalog_payload_invalid()
    try:
        return _validate_and_render_slot(
            field_name,
            slot,
            records,
            include_candidates=include_candidates,
            for_write=for_write,
        )
    except CatalogPreviewError as exc:
        raise _from_preview(exc) from exc


def _render_sample_pins(
    pins: list[dict[str, Any]],
    records: Mapping[str, TermRevisionRecord],
    *,
    include_candidates: bool,
    for_write: bool,
) -> list[dict[str, Any]]:
    return [
        {
            **pin,
            "replicate_role": _render_slot(
                "replicate_role",
                pin["replicate_role"],
                records,
                include_candidates=include_candidates,
                for_write=for_write,
            ),
        }
        for pin in pins
    ]


def _require_experiment_owner(
    owner: tuple[str, bool] | None,
    project_id: str,
    *,
    path: str,
    message: str,
) -> None:
    if owner is None or owner[0] != project_id:
        raise CatalogError(
            409,
            "catalog_proposal_invalidated",
            message,
            path=path,
        )
    if owner[1]:
        raise CatalogError(
            409,
            "catalog_experiment_archived",
            "Cited experiment is archived.",
            path=path,
        )


def _require_sample_owner(
    owner: str | None,
    project_id: str,
    *,
    path: str,
    message: str,
) -> None:
    if owner is None or owner != project_id:
        raise CatalogError(
            409,
            "catalog_proposal_invalidated",
            message,
            path=path,
        )


def _verify_projection_ref(
    ref_id: str,
    project_id: str,
    catalog: SourceCatalog,
    *,
    path: str,
) -> None:
    experiment = catalog.cited_experiment(ref_id)
    if experiment is not None:
        _require_experiment_owner(
            experiment,
            project_id,
            path=path,
            message="Cited catalog projection row is missing or outside this project.",
        )
        return
    revision = catalog.cited_experiment_revision(ref_id)
    if revision is not None:
        _require_experiment_owner(
            revision,
            project_id,
            path=path,
            message="Cited catalog projection row is missing or outside this project.",
        )
        return
    sample = catalog.cited_sample(ref_id)
    if sample is None:
        sample = catalog.cited_sample_revision(ref_id)
    _require_sample_owner(
        sample,
        project_id,
        path=path,
        message="Cited catalog projection row does not exist in this project.",
    )
