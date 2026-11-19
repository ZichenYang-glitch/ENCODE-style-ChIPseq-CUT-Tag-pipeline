"""M2 offline catalog curation preview, subset selection, conflict, and idempotency contracts.

This module implements the standalone M2 offline preview domain contract and an
in-memory synthetic catalog fixture aligned with the M1 Slice A read/write
baseline. It intentionally defines no Alembic migrations, HTTP routes, or M1 B
snapshot/run bindings.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from hashlib import sha256
import json
import re
from typing import Any


M2_CONTRACT_VERSION = "hw-m2-catalog-preview/1"
M2_RULE_VERSION = "m2-catalog-rules/1"
M2_VOCABULARY_VERSION = "terms-snapshot/2026-11-18"
M2_PROJECTION_VERSION = "hw-m2-catalog/1"

MAX_CHANGESET_ITEMS = 100
MAX_PREVIEW_APPLY_BODY_BYTES = 256 * 1024
PROPOSAL_ADOPTION_WINDOW = timedelta(minutes=30)
EXPIRED_PROPOSAL_RETENTION = timedelta(days=7)

EXPERIMENT_ALLOWED_CURATION_FIELDS = frozenset(
    {"display_name", "assay", "target", "batch", "samples"}
)
SAMPLE_ALLOWED_CURATION_FIELDS = frozenset(
    {
        "display_name",
        "species",
        "tissue",
        "cell_type",
        "cell_line",
        "treatment",
        "tag",
    }
)
SAMPLE_SINGLE_FIELDS = frozenset({"species", "tissue"})
SAMPLE_MULTI_FIELDS = frozenset({"cell_type", "cell_line", "treatment", "tag"})
EXPERIMENT_SINGLE_FIELDS = frozenset({"assay", "target", "batch"})
# Term slots only. Experiment `samples` is a sample-revision pin list, not a term.
TERM_SLOT_FIELDS = (
    EXPERIMENT_SINGLE_FIELDS
    | SAMPLE_SINGLE_FIELDS
    | SAMPLE_MULTI_FIELDS
    | frozenset({"replicate_role"})
)
FIELD_STATUSES = frozenset(
    {"known", "unknown", "not_applicable", "not_detected", "ambiguous"}
)
NEGATIVE_STATUSES = frozenset({"unknown", "not_applicable", "not_detected"})
_SAMPLE_PIN_KEYS = frozenset({"ordinal", "sample_revision_id", "replicate_role"})
_REPLICATE_GROUP_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]{0,254}$")

ALLOWED_REASON_CODES = frozenset(
    {
        "BATCH_TERM_MATCH",
        "TARGET_TERM_MATCH",
        "CELL_TYPE_TERM_MATCH",
        "TAG_TERM_MATCH",
        "REPLICATE_ROLE_TERM_MATCH",
        "ALIAS_AMBIGUOUS_CANDIDATES",
        "MISSING_EVIDENCE_UNKNOWN",
    }
)

_ID_PATTERNS = {
    "project": re.compile(r"^prj_[0-9a-f]{32}$"),
    "experiment": re.compile(r"^exp_[0-9a-f]{32}$"),
    "experiment_revision": re.compile(r"^expr_[0-9a-f]{32}$"),
    "sample": re.compile(r"^smp_[0-9a-f]{32}$"),
    "sample_revision": re.compile(r"^smpr_[0-9a-f]{32}$"),
    "term": re.compile(r"^term_[0-9a-f]{32}$"),
    "term_revision": re.compile(r"^trmv_[0-9a-f]{32}$"),
    "proposal": re.compile(r"^cprop_[0-9a-f]{32}$"),
    "receipt": re.compile(r"^cadopt_[0-9a-f]{32}$"),
}
_IDEMPOTENCY_KEY_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_ITEM_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,63}$")
_EVIDENCE_TOKEN_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:-]{0,63}$")
_DIGEST_RE = re.compile(r"^[0-9a-f]{64}$")
_WINDOWS_DRIVE_RE = re.compile(r"^[A-Za-z]:[\\/]")


def _canonical_json(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _sha256_canonical(value: Any) -> str:
    return sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _format_utc(value: datetime) -> str:
    if value.tzinfo is None or value.utcoffset() is None:
        raise ValueError("timestamp must be timezone-aware")
    return (
        value.astimezone(timezone.utc)
        .isoformat(timespec="seconds")
        .replace("+00:00", "Z")
    )


def _check_id(kind: str, value: object, *, path: str | None = None) -> str:
    pattern = _ID_PATTERNS[kind]
    if not isinstance(value, str) or pattern.fullmatch(value) is None:
        raise catalog_preview_error(
            422,
            "catalog_payload_invalid",
            f"Invalid {kind} identifier.",
            path=path,
        )
    return value


def _assert_bounded_text(value: object, *, path: str | None = None) -> str:
    if not isinstance(value, str) or not value or value != value.strip():
        raise catalog_preview_error(
            422,
            "catalog_payload_invalid",
            "Bounded token must be a non-empty trimmed string.",
            path=path,
        )
    if len(value) > 255 or any(ord(ch) < 32 for ch in value):
        raise catalog_preview_error(
            422,
            "catalog_payload_invalid",
            "Bounded token exceeds length or control-character limits.",
            path=path,
        )
    # Reject host path leakage in proposal text or evidence tokens.
    if value.startswith(("/", "~", "\\")) or _WINDOWS_DRIVE_RE.match(value):
        raise catalog_preview_error(
            422,
            "catalog_payload_invalid",
            "Host filesystem paths are not allowed in catalog proposals.",
            path=path,
        )
    return value


def _assert_body_size_within_limit(value: Any) -> None:
    encoded = _canonical_json(value).encode("utf-8")
    if len(encoded) > MAX_PREVIEW_APPLY_BODY_BYTES:
        raise catalog_preview_error(
            413,
            "catalog_payload_invalid",
            "Request body exceeds the 256 KiB catalog changeset limit.",
        )


def _iter_term_slots(slot_map: Mapping[str, Any]) -> list[Mapping[str, Any]]:
    """Yield classification term slots, including replicate_role inside sample pins."""

    slots: list[Mapping[str, Any]] = []
    for field_name, slot_value in slot_map.items():
        if field_name == "display_name":
            continue
        if field_name == "samples" and isinstance(slot_value, list):
            for pin in slot_value:
                if not isinstance(pin, Mapping):
                    continue
                role = pin.get("replicate_role")
                if isinstance(role, Mapping):
                    slots.append(role)
            continue
        if isinstance(slot_value, Mapping):
            slots.append(slot_value)
    return slots


def _collect_ambiguous_candidate_revision_ids(
    *slot_maps: Mapping[str, Any],
) -> frozenset[str]:
    """Return every term_revision_id referenced inside any ambiguous slot's candidates."""

    candidate_ids: set[str] = set()
    for slot_map in slot_maps:
        for slot_value in _iter_term_slots(slot_map):
            if slot_value.get("status") != "ambiguous":
                continue
            raw_candidates = slot_value.get("candidates") or ()
            if isinstance(raw_candidates, (str, bytes)) or not isinstance(
                raw_candidates, Sequence
            ):
                continue
            for candidate in raw_candidates:
                if isinstance(candidate, Mapping):
                    rev_id = candidate.get("term_revision_id")
                    if isinstance(rev_id, str):
                        candidate_ids.add(rev_id)
    return frozenset(candidate_ids)


@dataclass(frozen=True)
class RevisionConflictDetail:
    """Bounded per-target optimistic concurrency mismatch."""

    target_kind: str
    target_id: str
    expected_revision_number: int
    current_revision_number: int
    current_revision_id: str
    expected_revision_id: str | None = None

    def to_object(self) -> dict[str, Any]:
        obj: dict[str, Any] = {
            "target_kind": self.target_kind,
            "target_id": self.target_id,
            "expected_revision_number": self.expected_revision_number,
            "current_revision_number": self.current_revision_number,
            "current_revision_id": self.current_revision_id,
        }
        if self.expected_revision_id is not None:
            obj["expected_revision_id"] = self.expected_revision_id
        return obj


class CatalogPreviewError(Exception):
    """Stable M2 preview/adoption failure matching M1 A issue envelope fields."""

    def __init__(
        self,
        status_code: int,
        code: str,
        message: str,
        *,
        path: str | None = None,
        hint: str | None = None,
        conflicts: Sequence[RevisionConflictDetail] = (),
    ) -> None:
        super().__init__(message)
        self.status_code = status_code
        self.code = code
        self.message = message
        self.path = path
        self.hint = hint
        self.conflicts = tuple(conflicts)

    def to_response(self) -> dict[str, Any]:
        body: dict[str, Any] = {
            "ok": False,
            "issues": [
                {
                    "code": self.code,
                    "message": self.message,
                    "severity": "error",
                    "path": self.path,
                    "source": "catalog_m2_preview",
                    "hint": self.hint,
                }
            ],
        }
        if self.conflicts:
            body["conflicts"] = [item.to_object() for item in self.conflicts]
        return body


def catalog_preview_error(
    status_code: int,
    code: str,
    message: str,
    *,
    path: str | None = None,
    hint: str | None = None,
    conflicts: Sequence[RevisionConflictDetail] = (),
) -> CatalogPreviewError:
    return CatalogPreviewError(
        status_code,
        code,
        message,
        path=path,
        hint=hint,
        conflicts=conflicts,
    )


@dataclass(frozen=True)
class CatalogActor:
    """Authenticated platform actor evaluated on each read or write transaction."""

    user_id: str
    role: str
    active: bool = True


@dataclass(frozen=True)
class TermRevisionRecord:
    """Immutable controlled-vocabulary revision aligned with M1 A."""

    term_id: str
    term_revision_id: str
    term_type: str
    display_name: str
    aliases: tuple[str, ...] = ()
    archived: bool = False


@dataclass(frozen=True)
class SourceReference:
    """Platform-verifiable provenance pointer for an offline suggestion."""

    kind: str
    ref_id: str
    term_revision_id: str | None = None
    evidence_token: str | None = None
    visibility: str = "shared"

    def __post_init__(self) -> None:
        if self.kind == "experiment_revision":
            _check_id("experiment_revision", self.ref_id, path="citations.source_refs")
        elif self.kind == "sample_revision":
            _check_id("sample_revision", self.ref_id, path="citations.source_refs")
        elif self.kind == "term_revision":
            _check_id("term_revision", self.ref_id, path="citations.source_refs")
        elif self.kind == "catalog_projection_row":
            _assert_bounded_text(self.ref_id, path="citations.source_refs")
        else:
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                "Unsupported source reference kind.",
                path="citations.source_refs",
            )
        if self.term_revision_id is not None:
            _check_id(
                "term_revision",
                self.term_revision_id,
                path="citations.source_refs",
            )
            if self.kind == "term_revision" and self.term_revision_id != self.ref_id:
                raise catalog_preview_error(
                    422,
                    "catalog_payload_invalid",
                    "term_revision source ref must match its term_revision_id.",
                    path="citations.source_refs",
                )
        if self.evidence_token is not None:
            _assert_bounded_text(self.evidence_token, path="citations.source_refs")
            if _EVIDENCE_TOKEN_RE.fullmatch(self.evidence_token) is None:
                raise catalog_preview_error(
                    422,
                    "catalog_payload_invalid",
                    "Invalid evidence token format.",
                    path="citations.source_refs",
                )
        if self.visibility not in {"shared", "admin_only"}:
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                "Unsupported source reference visibility.",
                path="citations.source_refs",
            )

    def is_admin_only(self, candidate_term_revision_ids: frozenset[str]) -> bool:
        """Return whether this ref identifies an ambiguous candidate hidden from members."""

        if self.visibility == "admin_only":
            return True
        if self.ref_id in candidate_term_revision_ids:
            return True
        return (
            self.term_revision_id is not None
            and self.term_revision_id in candidate_term_revision_ids
        )

    def to_object(self) -> dict[str, Any]:
        obj: dict[str, Any] = {"kind": self.kind, "ref_id": self.ref_id}
        if self.term_revision_id is not None:
            obj["term_revision_id"] = self.term_revision_id
        if self.evidence_token is not None:
            obj["evidence_token"] = self.evidence_token
        return obj


@dataclass(frozen=True)
class ProposalCitation:
    """Bounded suggestion provenance without free-text reasoning or host paths."""

    rule_version: str
    vocabulary_version: str
    reason_codes: tuple[str, ...]
    source_refs: tuple[SourceReference, ...]

    def __post_init__(self) -> None:
        _assert_bounded_text(self.rule_version, path="citations.rule_version")
        _assert_bounded_text(
            self.vocabulary_version, path="citations.vocabulary_version"
        )
        if not self.reason_codes:
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                "Proposal citation must include at least one reason code.",
                path="citations.reason_codes",
            )
        for code in self.reason_codes:
            if code not in ALLOWED_REASON_CODES:
                raise catalog_preview_error(
                    422,
                    "catalog_payload_invalid",
                    "Unknown proposal citation reason code.",
                    path="citations.reason_codes",
                )
        if not self.source_refs:
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                "Proposal citation must include at least one verifiable source ref.",
                path="citations.source_refs",
            )

    def to_object(
        self,
        *,
        include_candidates: bool = True,
        candidate_term_revision_ids: frozenset[str] = frozenset(),
    ) -> dict[str, Any]:
        refs = (
            self.source_refs
            if include_candidates
            else tuple(
                ref
                for ref in self.source_refs
                if not ref.is_admin_only(candidate_term_revision_ids)
            )
        )
        return {
            "rule_version": self.rule_version,
            "vocabulary_version": self.vocabulary_version,
            "reason_codes": list(self.reason_codes),
            "source_refs": [ref.to_object() for ref in refs],
        }


@dataclass(frozen=True)
class ProposalItem:
    """Single bounded object modification inside an offline M2 proposal."""

    item_id: str
    target_kind: str
    target_id: str
    base_revision_id: str
    base_revision_number: int
    allowed_fields: tuple[str, ...]
    before: dict[str, Any]
    after: dict[str, Any]
    citation: ProposalCitation

    def __post_init__(self) -> None:
        if _ITEM_ID_RE.fullmatch(self.item_id) is None:
            raise catalog_preview_error(
                422, "catalog_payload_invalid", "Invalid proposal item id."
            )
        if self.base_revision_number < 1:
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                "Base revision number must be at least 1.",
            )
        if self.target_kind == "experiment":
            _check_id("experiment", self.target_id, path="items.target_id")
            _check_id(
                "experiment_revision",
                self.base_revision_id,
                path="items.base_revision_id",
            )
            allowed_pool = EXPERIMENT_ALLOWED_CURATION_FIELDS
        elif self.target_kind == "sample_catalog":
            _check_id("sample", self.target_id, path="items.target_id")
            _check_id(
                "sample_revision",
                self.base_revision_id,
                path="items.base_revision_id",
            )
            allowed_pool = SAMPLE_ALLOWED_CURATION_FIELDS
        else:
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                "Unsupported curation target kind.",
                path="items.target_kind",
            )
        if not self.allowed_fields or len(set(self.allowed_fields)) != len(
            self.allowed_fields
        ):
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                "Allowed fields must be non-empty and unique.",
            )
        if not set(self.allowed_fields).issubset(allowed_pool):
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                "Proposal touches fields outside the M2 allowlist.",
            )
        if set(self.before) != set(self.allowed_fields) or set(self.after) != set(
            self.allowed_fields
        ):
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                "Before/after snapshots must match allowed_fields exactly.",
            )
        if self.target_kind == "experiment" and "samples" in self.allowed_fields:
            for side in (self.before["samples"], self.after["samples"]):
                if isinstance(side, Mapping) or not isinstance(side, list):
                    raise catalog_preview_error(
                        422,
                        "catalog_payload_invalid",
                        "Experiment samples are sample-revision associations, not a classification term.",
                        path="items.samples",
                    )

    def candidate_term_revision_ids(self) -> frozenset[str]:
        """Return candidate term_revision_ids from any ambiguous slot in this item."""

        return _collect_ambiguous_candidate_revision_ids(self.before, self.after)

    def canonical_object(self) -> dict[str, Any]:
        return {
            "item_id": self.item_id,
            "target_kind": self.target_kind,
            "target_id": self.target_id,
            "base_revision_id": self.base_revision_id,
            "base_revision_number": self.base_revision_number,
            "allowed_fields": list(self.allowed_fields),
            "before": self.before,
            "after": self.after,
            "citations": self.citation.to_object(include_candidates=True),
        }


@dataclass(frozen=True)
class CatalogProposal:
    """Immutable offline suggestion package separated from administrator adoption."""

    proposal_id: str
    project_id: str
    created_by_user_id: str
    source: str
    rule_version: str
    vocabulary_version: str
    projection_version: str
    items: tuple[ProposalItem, ...]
    proposal_digest: str
    created_at: datetime
    expires_at: datetime
    retention_until: datetime

    @classmethod
    def create(
        cls,
        *,
        proposal_id: str,
        project_id: str,
        created_by_user_id: str,
        source: str,
        items: Sequence[ProposalItem],
        created_at: datetime,
        rule_version: str = M2_RULE_VERSION,
        vocabulary_version: str = M2_VOCABULARY_VERSION,
        projection_version: str = M2_PROJECTION_VERSION,
    ) -> "CatalogProposal":
        _check_id("proposal", proposal_id, path="proposal_id")
        _check_id("project", project_id, path="project_id")
        if source not in {"offline_suggestion", "agent", "manual"}:
            raise catalog_preview_error(
                422, "catalog_payload_invalid", "Unsupported proposal source."
            )
        if not items or len(items) > MAX_CHANGESET_ITEMS:
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                f"Proposal must contain 1 to {MAX_CHANGESET_ITEMS} items.",
            )
        item_ids = [item.item_id for item in items]
        target_ids = [item.target_id for item in items]
        if len(set(item_ids)) != len(item_ids) or len(set(target_ids)) != len(
            target_ids
        ):
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                "Proposal item_id and target_id values must be unique per changeset.",
            )
        tuple_items = tuple(items)
        digest = compute_proposal_digest(
            project_id=project_id,
            rule_version=rule_version,
            vocabulary_version=vocabulary_version,
            projection_version=projection_version,
            items=tuple_items,
        )
        return cls(
            proposal_id=proposal_id,
            project_id=project_id,
            created_by_user_id=created_by_user_id,
            source=source,
            rule_version=rule_version,
            vocabulary_version=vocabulary_version,
            projection_version=projection_version,
            items=tuple_items,
            proposal_digest=digest,
            created_at=created_at,
            expires_at=created_at + PROPOSAL_ADOPTION_WINDOW,
            retention_until=created_at + EXPIRED_PROPOSAL_RETENTION,
        )


def compute_proposal_digest(
    *,
    project_id: str,
    rule_version: str,
    vocabulary_version: str,
    projection_version: str,
    items: Sequence[ProposalItem],
) -> str:
    """Compute the canonical SHA-256 digest of an immutable M2 proposal."""

    return _sha256_canonical(
        {
            "contract_version": M2_CONTRACT_VERSION,
            "project_id": project_id,
            "rule_version": rule_version,
            "vocabulary_version": vocabulary_version,
            "projection_version": projection_version,
            "items": [item.canonical_object() for item in items],
        }
    )


def compute_preview_digest(
    *,
    proposal_id: str,
    proposal_digest: str,
    selected_item_ids: Sequence[str],
    expected_revisions: Mapping[str, int],
) -> str:
    """Compute the immutable SHA-256 digest binding a selected subset and base revisions."""

    canonical_expected = {
        target_id: expected_revisions[target_id]
        for target_id in sorted(expected_revisions)
    }
    return _sha256_canonical(
        {
            "contract_version": M2_CONTRACT_VERSION,
            "proposal_id": proposal_id,
            "proposal_digest": proposal_digest,
            "selected_item_ids": list(selected_item_ids),
            "expected_revisions": canonical_expected,
        }
    )


def _validate_idempotency_key(value: str | None) -> str:
    if not isinstance(value, str) or _IDEMPOTENCY_KEY_RE.fullmatch(value) is None:
        raise catalog_preview_error(
            400,
            "API_REQUEST_INVALID",
            "Idempotency-Key header is required.",
        )
    return value


def _slot_is_multi(field_name: str) -> bool:
    return field_name in SAMPLE_MULTI_FIELDS


def _drop_rendered_display_names(value: Any) -> Any:
    """Keep verified term identity and drop preview-only display names."""

    if isinstance(value, list):
        return [_drop_rendered_display_names(item) for item in value]
    if isinstance(value, dict):
        return {
            key: _drop_rendered_display_names(item)
            for key, item in value.items()
            if key != "display_name"
        }
    return value


def _sample_state_for_revision(
    samples: Mapping[str, Any], revision_id: str
) -> Any | None:
    for state in samples.values():
        for revision in state.revisions:
            if revision["sample_revision_id"] == revision_id:
                return state
    return None


def _validate_and_render_sample_pins(
    pins: object,
    terms_by_revision: Mapping[str, TermRevisionRecord],
    samples: Mapping[str, Any],
    *,
    project_id: str,
    include_candidates: bool,
    for_write: bool,
) -> list[dict[str, Any]]:
    """Validate sample-revision pins. A term-shaped object is rejected here."""

    if isinstance(pins, Mapping):
        raise catalog_preview_error(
            422,
            "catalog_payload_invalid",
            "Experiment samples are sample-revision associations, not a classification term.",
            path="samples",
        )
    if not isinstance(pins, list) or not pins:
        raise catalog_preview_error(
            422,
            "catalog_payload_invalid",
            "Experiment samples must be a non-empty list of sample-revision pins.",
            path="samples",
        )
    rendered: list[dict[str, Any]] = []
    seen: set[str] = set()
    for index, pin in enumerate(pins):
        if not isinstance(pin, Mapping):
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                "Each sample association must be an object.",
                path="samples",
            )
        keys = set(pin)
        allowed = set(_SAMPLE_PIN_KEYS)
        if "replicate_group" in keys:
            allowed.add("replicate_group")
        if keys != allowed:
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                "Sample-revision pin contains unknown attributes.",
                path="samples",
            )
        ordinal = pin.get("ordinal")
        if (
            isinstance(ordinal, bool)
            or not isinstance(ordinal, int)
            or ordinal != index
        ):
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                "Sample-revision pin ordinal must match its list position.",
                path="samples",
            )
        revision_id = _check_id(
            "sample_revision", pin.get("sample_revision_id"), path="samples"
        )
        if revision_id in seen:
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                "Duplicate sample_revision_id in samples.",
                path="samples",
            )
        seen.add(revision_id)
        state = _sample_state_for_revision(samples, revision_id)
        if state is None or state.project_id != project_id:
            raise catalog_preview_error(
                409,
                "catalog_proposal_invalidated",
                "Sample-revision pin is missing or outside this project.",
                path="samples",
            )
        normalized: dict[str, Any] = {
            "ordinal": ordinal,
            "sample_revision_id": revision_id,
            "replicate_role": _validate_and_render_slot(
                "replicate_role",
                pin.get("replicate_role"),
                terms_by_revision,
                include_candidates=include_candidates,
                for_write=for_write,
            ),
        }
        if "replicate_group" in pin:
            group = pin.get("replicate_group")
            if (
                not isinstance(group, str)
                or _REPLICATE_GROUP_RE.fullmatch(group) is None
            ):
                raise catalog_preview_error(
                    422,
                    "catalog_payload_invalid",
                    "replicate_group must be a bounded token.",
                    path="samples",
                )
            normalized["replicate_group"] = group
        rendered.append(normalized)
    return rendered


def _validate_and_render_slot(
    field_name: str,
    slot_value: object,
    terms_by_revision: Mapping[str, TermRevisionRecord],
    *,
    include_candidates: bool,
    for_write: bool,
) -> dict[str, Any]:
    """Validate M1 A slot invariants and render bounded public/admin DTOs."""

    if field_name == "display_name":
        raise ValueError("display_name is scalar text, not a term slot")
    if field_name == "samples":
        raise catalog_preview_error(
            422,
            "catalog_payload_invalid",
            "Experiment samples are sample-revision associations, not a classification term.",
            path="samples",
        )
    if field_name not in TERM_SLOT_FIELDS:
        raise catalog_preview_error(
            422,
            "catalog_payload_invalid",
            f"{field_name} is not a classification term field.",
            path=field_name,
        )
    if not isinstance(slot_value, Mapping):
        raise catalog_preview_error(
            422,
            "catalog_payload_invalid",
            f"{field_name} must be a classification term slot.",
            path=field_name,
        )
    expected_type = "replicate_role" if field_name == "replicate_role" else field_name
    is_multi = _slot_is_multi(field_name)
    status = slot_value.get("status")
    if status not in FIELD_STATUSES:
        raise catalog_preview_error(
            422, "catalog_payload_invalid", f"Invalid status on {field_name}."
        )

    def resolve_pin(raw_pin: Mapping[str, Any]) -> dict[str, str]:
        term_id = _check_id("term", raw_pin.get("term_id"), path=field_name)
        term_rev_id = _check_id(
            "term_revision", raw_pin.get("term_revision_id"), path=field_name
        )
        record = terms_by_revision.get(term_rev_id)
        if record is None or record.term_id != term_id:
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                f"Unknown term revision reference on {field_name}.",
                path=field_name,
            )
        if record.term_type != expected_type:
            raise catalog_preview_error(
                409,
                "catalog_term_type_mismatch",
                f"Term type does not match {field_name}.",
                path=field_name,
            )
        if for_write and record.archived:
            raise catalog_preview_error(
                409,
                "catalog_term_archived",
                f"Archived term cannot be used in {field_name}.",
                path=field_name,
            )
        return {
            "term_id": record.term_id,
            "term_revision_id": record.term_revision_id,
            "display_name": record.display_name,
        }

    def resolve_nested_pin(raw_pin: object) -> dict[str, str]:
        if not isinstance(raw_pin, Mapping) or set(raw_pin) != {
            "term_id",
            "term_revision_id",
        }:
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                f"Term pin on {field_name} must contain only term_id and term_revision_id.",
                path=field_name,
            )
        return resolve_pin(raw_pin)

    if not is_multi:
        if status == "known":
            if set(slot_value) != {"status", "term_id", "term_revision_id"}:
                raise catalog_preview_error(
                    422,
                    "catalog_payload_invalid",
                    f"Known single slot {field_name} must specify one pinned term.",
                )
            resolved = resolve_pin(slot_value)
            return {"status": "known", **resolved}
        if status in NEGATIVE_STATUSES:
            if set(slot_value) != {"status"}:
                raise catalog_preview_error(
                    422,
                    "catalog_payload_invalid",
                    f"Negative slot {field_name} must not carry term pins.",
                )
            return {"status": status}
        if set(slot_value) not in ({"status"}, {"status", "candidates"}):
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                f"Invalid ambiguous single slot keys on {field_name}.",
            )
        rendered: dict[str, Any] = {"status": "ambiguous"}
        raw_candidates = slot_value.get("candidates", [])
        if not isinstance(raw_candidates, list):
            raise catalog_preview_error(
                422, "catalog_payload_invalid", "Ambiguous candidates must be a list."
            )
        resolved_candidates = [resolve_nested_pin(item) for item in raw_candidates]
        if include_candidates:
            rendered["candidates"] = resolved_candidates
        return rendered

    # Multi-value slot (`cell_type`, `cell_line`, `treatment`, `tag`)
    if status == "known":
        if set(slot_value) != {"status", "values"}:
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                f"Known multi slot {field_name} must only have status and values.",
            )
        raw_values = slot_value.get("values")
        if not isinstance(raw_values, list) or not raw_values:
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                f"Known multi slot {field_name} requires at least one pinned term.",
            )
        return {
            "status": "known",
            "values": [resolve_nested_pin(item) for item in raw_values],
        }
    if slot_value.get("values") != []:
        raise catalog_preview_error(
            422,
            "catalog_payload_invalid",
            f"Non-known multi slot {field_name} must have empty values.",
        )
    if status in NEGATIVE_STATUSES:
        if set(slot_value) != {"status", "values"}:
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                f"Negative multi slot {field_name} must not carry candidates.",
            )
        return {"status": status, "values": []}
    if set(slot_value) not in (
        {"status", "values"},
        {"status", "values", "candidates"},
    ):
        raise catalog_preview_error(
            422,
            "catalog_payload_invalid",
            f"Invalid ambiguous multi slot keys on {field_name}.",
        )
    raw_candidates = slot_value.get("candidates", [])
    if not isinstance(raw_candidates, list):
        raise catalog_preview_error(
            422, "catalog_payload_invalid", "Ambiguous candidates must be a list."
        )
    resolved_candidates = [resolve_nested_pin(item) for item in raw_candidates]
    rendered = {"status": "ambiguous", "values": []}
    if include_candidates:
        rendered["candidates"] = resolved_candidates
    return rendered


@dataclass(frozen=True)
class _IdempotencyReceipt:
    """Receipt keyed by (actor_user_id, idempotency_key) to match M1 A scope."""

    operation: str
    target_key: str
    request_sha256: str
    response: dict[str, Any]


@dataclass
class _ExperimentState:
    experiment_id: str
    project_id: str
    stable_key: str
    archived: bool
    revisions: list[dict[str, Any]]

    @property
    def head(self) -> dict[str, Any]:
        return self.revisions[-1]


@dataclass
class _SampleState:
    sample_id: str
    project_id: str
    stable_key: str
    revisions: list[dict[str, Any]]

    @property
    def head(self) -> dict[str, Any]:
        return self.revisions[-1]


class SyntheticCatalogPreviewStore:
    """In-memory M1 A-aligned catalog store supporting M2 preview, conflict, and idempotency."""

    def __init__(
        self,
        *,
        projects: Mapping[str, dict[str, Any]],
        actors: Mapping[str, CatalogActor],
        terms_by_revision: Mapping[str, TermRevisionRecord],
        experiments: Mapping[str, _ExperimentState],
        samples: Mapping[str, _SampleState],
        proposals: Mapping[str, CatalogProposal],
        active_rule_version: str = M2_RULE_VERSION,
        active_vocabulary_version: str = M2_VOCABULARY_VERSION,
        active_projection_version: str = M2_PROJECTION_VERSION,
        revision_id_factory: Callable[[str], str] | None = None,
    ) -> None:
        self._projects = dict(projects)
        self._actors = dict(actors)
        self._terms = dict(terms_by_revision)
        self._experiments = dict(experiments)
        self._samples = dict(samples)
        self._proposals = dict(proposals)
        self._active_rule_version = active_rule_version
        self._active_vocabulary_version = active_vocabulary_version
        self._active_projection_version = active_projection_version
        self._proposal_adoptions: dict[str, dict[str, Any]] = {}
        # M1 A key scope: (actor_user_id, idempotency_key), then check operation/target/body.
        self._idempotency_receipts: dict[tuple[str, str], _IdempotencyReceipt] = {}
        self._adoption_receipts: list[dict[str, Any]] = []
        self._counter = 0
        self._revision_id_factory = revision_id_factory or self._next_prefixed_id

    def _next_prefixed_id(self, prefix: str) -> str:
        self._counter += 1
        return f"{prefix}{self._counter:032x}"

    def set_actor_role(self, user_id: str, *, role: str, active: bool = True) -> None:
        """Update actor status/role to simulate admin revocation prior to commit."""

        existing = self._actors[user_id]
        self._actors[user_id] = CatalogActor(
            user_id=existing.user_id, role=role, active=active
        )

    def set_active_versions(
        self,
        *,
        rule_version: str | None = None,
        vocabulary_version: str | None = None,
        projection_version: str | None = None,
    ) -> None:
        """Advance active catalog rule/vocabulary versions to test proposal invalidation."""

        if rule_version is not None:
            self._active_rule_version = rule_version
        if vocabulary_version is not None:
            self._active_vocabulary_version = vocabulary_version
        if projection_version is not None:
            self._active_projection_version = projection_version

    def archive_term_revision(self, term_revision_id: str) -> None:
        """Mark a term revision's parent term as archived to test write rejection."""

        record = self._terms[term_revision_id]
        for rev_id, item in list(self._terms.items()):
            if item.term_id == record.term_id:
                self._terms[rev_id] = TermRevisionRecord(
                    term_id=item.term_id,
                    term_revision_id=item.term_revision_id,
                    term_type=item.term_type,
                    display_name=item.display_name,
                    aliases=item.aliases,
                    archived=True,
                )

    def get_experiment_head(self, experiment_id: str) -> dict[str, Any]:
        return deepcopy(self._experiments[experiment_id].head)

    def get_sample_head(self, sample_id: str) -> dict[str, Any]:
        return deepcopy(self._samples[sample_id].head)

    def _require_actor(self, user_id: str, *, require_admin: bool) -> CatalogActor:
        actor = self._actors.get(user_id)
        if actor is None or not actor.active:
            raise catalog_preview_error(
                401,
                "AUTHENTICATION_REQUIRED",
                "Active authenticated principal is required.",
            )
        if require_admin and actor.role != "administrator":
            raise catalog_preview_error(
                403,
                "ADMINISTRATOR_REQUIRED",
                "Administrator role is required to preview or adopt catalog changes.",
            )
        return actor

    def _require_active_project(self, project_id: str) -> None:
        project = self._projects.get(project_id)
        if project is None or not project.get("active") or project.get("legacy"):
            raise catalog_preview_error(
                404, "catalog_not_found", "The catalog project was not found."
            )

    def _require_proposal(
        self, proposal_id: str, *, now: datetime, for_adoption: bool
    ) -> CatalogProposal:
        proposal = self._proposals.get(proposal_id)
        if proposal is None:
            raise catalog_preview_error(
                404, "catalog_not_found", "The catalog proposal was not found."
            )
        self._require_active_project(proposal.project_id)
        adopted = self._proposal_adoptions.get(proposal_id)
        if adopted is not None:
            if for_adoption:
                raise catalog_preview_error(
                    409,
                    "catalog_proposal_invalidated",
                    "This proposal has already been adopted.",
                    hint="Reopen the proposal or query adoption receipts to view the committed result.",
                )
            # Adopted proposals and receipts remain readable beyond TTL as audit evidence.
            return proposal
        if now >= proposal.retention_until:
            raise catalog_preview_error(
                404, "catalog_not_found", "The catalog proposal was not found."
            )
        if for_adoption and now >= proposal.expires_at:
            raise catalog_preview_error(
                410,
                "catalog_proposal_expired",
                "The proposal adoption window has expired.",
                hint="Read-only inspection remains available within the 7-day retention window.",
            )
        return proposal

    def _verify_proposal_versions(self, proposal: CatalogProposal) -> None:
        if (
            proposal.rule_version != self._active_rule_version
            or proposal.vocabulary_version != self._active_vocabulary_version
            or proposal.projection_version != self._active_projection_version
        ):
            raise catalog_preview_error(
                409,
                "catalog_proposal_invalidated",
                "Proposal rule, vocabulary, or projection version is no longer active.",
            )
        for item in proposal.items:
            if (
                item.citation.rule_version != self._active_rule_version
                or item.citation.vocabulary_version != self._active_vocabulary_version
            ):
                raise catalog_preview_error(
                    409,
                    "catalog_proposal_invalidated",
                    "Item citation rule or vocabulary version is no longer active.",
                    path=item.item_id,
                )

    def _find_experiment_revision(
        self, revision_id: str
    ) -> tuple[_ExperimentState, dict[str, Any]] | None:
        for exp_state in self._experiments.values():
            for rev in exp_state.revisions:
                if rev["experiment_revision_id"] == revision_id:
                    return exp_state, rev
        return None

    def _find_sample_revision(
        self, revision_id: str
    ) -> tuple[_SampleState, dict[str, Any]] | None:
        for smp_state in self._samples.values():
            for rev in smp_state.revisions:
                if rev["sample_revision_id"] == revision_id:
                    return smp_state, rev
        return None

    def _verify_experiment_ref(
        self,
        exp_state: _ExperimentState | None,
        project_id: str,
        *,
        path: str,
        message_missing: str,
    ) -> None:
        if exp_state is None or exp_state.project_id != project_id:
            raise catalog_preview_error(
                409,
                "catalog_proposal_invalidated",
                message_missing,
                path=path,
            )
        if exp_state.archived:
            raise catalog_preview_error(
                409,
                "catalog_experiment_archived",
                "Cited experiment is archived.",
                path=path,
            )

    def _verify_sample_ref(
        self,
        smp_state: _SampleState | None,
        project_id: str,
        *,
        path: str,
        message_missing: str,
    ) -> None:
        if smp_state is None or smp_state.project_id != project_id:
            raise catalog_preview_error(
                409,
                "catalog_proposal_invalidated",
                message_missing,
                path=path,
            )

    def _verify_projection_row_ref(
        self, ref_id: str, project_id: str, *, path: str
    ) -> None:
        exp_state = self._experiments.get(ref_id)
        if exp_state is None:
            found_exp = self._find_experiment_revision(ref_id)
            exp_state = None if found_exp is None else found_exp[0]
        if exp_state is not None:
            self._verify_experiment_ref(
                exp_state,
                project_id,
                path=path,
                message_missing="Cited catalog projection row is missing or outside this project.",
            )
            return
        smp_state = self._samples.get(ref_id)
        if smp_state is None:
            found_smp = self._find_sample_revision(ref_id)
            smp_state = None if found_smp is None else found_smp[0]
        self._verify_sample_ref(
            smp_state,
            project_id,
            path=path,
            message_missing="Cited catalog projection row does not exist in this project.",
        )

    def _verify_source_references(self, item: ProposalItem, project_id: str) -> None:
        """Verify every cited source reference exists in this active project and is valid."""

        path = f"{item.item_id}.citations.source_refs"
        for ref in item.citation.source_refs:
            if ref.kind == "experiment_revision":
                found_exp = self._find_experiment_revision(ref.ref_id)
                self._verify_experiment_ref(
                    None if found_exp is None else found_exp[0],
                    project_id,
                    path=path,
                    message_missing="Cited experiment revision is missing or outside this project.",
                )
            elif ref.kind == "sample_revision":
                found_smp = self._find_sample_revision(ref.ref_id)
                self._verify_sample_ref(
                    None if found_smp is None else found_smp[0],
                    project_id,
                    path=path,
                    message_missing="Cited sample revision is missing or outside this project.",
                )
            elif ref.kind == "term_revision":
                term_rec = self._terms.get(ref.ref_id)
                if term_rec is None:
                    raise catalog_preview_error(
                        409,
                        "catalog_proposal_invalidated",
                        "Cited term revision does not exist in the catalog.",
                        path=path,
                    )
                if term_rec.archived:
                    raise catalog_preview_error(
                        409,
                        "catalog_term_archived",
                        "Cited term revision is archived.",
                        path=path,
                    )
            elif ref.kind == "catalog_projection_row":
                self._verify_projection_row_ref(ref.ref_id, project_id, path=path)

    def _lookup_idempotency(
        self,
        *,
        actor_user_id: str,
        idempotency_key: str,
        operation: str,
        target_key: str,
        request_sha256: str,
    ) -> dict[str, Any] | None:
        """M1 A idempotency lookup by (actor_user_id, idempotency_key)."""

        existing = self._idempotency_receipts.get((actor_user_id, idempotency_key))
        if existing is None:
            return None
        if (
            existing.operation != operation
            or existing.target_key != target_key
            or existing.request_sha256 != request_sha256
        ):
            raise catalog_preview_error(
                409,
                "catalog_idempotency_conflict",
                "The idempotency key was already used for a different operation, target, or request body.",
            )
        replayed = deepcopy(existing.response)
        replayed["replayed"] = True
        return replayed

    def _project_id_for_item(self, item: ProposalItem) -> str:
        if item.target_kind == "experiment":
            state = self._experiments.get(item.target_id)
            if state is None:
                raise catalog_preview_error(
                    404, "catalog_not_found", "Target experiment was not found."
                )
            return state.project_id
        sample_state = self._samples.get(item.target_id)
        if sample_state is None:
            raise catalog_preview_error(
                404, "catalog_not_found", "Target sample was not found."
            )
        return sample_state.project_id

    def _render_field(
        self,
        field_name: str,
        value: object,
        *,
        project_id: str,
        include_candidates: bool,
        for_write: bool,
    ) -> Any:
        if field_name == "display_name":
            return _assert_bounded_text(value, path="display_name")
        if field_name == "samples":
            return _validate_and_render_sample_pins(
                value,
                self._terms,
                self._samples,
                project_id=project_id,
                include_candidates=include_candidates,
                for_write=for_write,
            )
        return _validate_and_render_slot(
            field_name,
            value,
            self._terms,
            include_candidates=include_candidates,
            for_write=for_write,
        )

    def _stored_field(self, field_name: str, rendered: Any) -> Any:
        """Persist the verified value. Preview display names stay out of the revision."""

        if field_name == "display_name":
            return rendered
        return _drop_rendered_display_names(rendered)

    def _render_item_diff(
        self,
        item: ProposalItem,
        *,
        include_candidates: bool,
        for_write: bool,
    ) -> tuple[dict[str, Any], dict[str, Any]]:
        project_id = self._project_id_for_item(item)
        before_rendered: dict[str, Any] = {}
        after_rendered: dict[str, Any] = {}
        for field_name in item.allowed_fields:
            before_rendered[field_name] = self._render_field(
                field_name,
                item.before[field_name],
                project_id=project_id,
                include_candidates=include_candidates,
                for_write=False,
            )
            after_rendered[field_name] = self._render_field(
                field_name,
                item.after[field_name],
                project_id=project_id,
                include_candidates=include_candidates,
                for_write=for_write,
            )
        return before_rendered, after_rendered

    def _head_revision_info(
        self, item: ProposalItem, expected_project_id: str
    ) -> tuple[int, str, bool]:
        if item.target_kind == "experiment":
            state = self._experiments.get(item.target_id)
            if state is None or state.project_id != expected_project_id:
                raise catalog_preview_error(
                    404, "catalog_not_found", "Target experiment was not found."
                )
            head = state.head
            return (
                int(head["revision_number"]),
                str(head["experiment_revision_id"]),
                bool(state.archived),
            )
        sample_state = self._samples.get(item.target_id)
        if sample_state is None or sample_state.project_id != expected_project_id:
            raise catalog_preview_error(
                404, "catalog_not_found", "Target sample was not found."
            )
        head = sample_state.head
        return int(head["revision_number"]), str(head["sample_revision_id"]), False

    def get_proposal_view(
        self,
        *,
        actor_user_id: str,
        proposal_id: str,
        now: datetime,
    ) -> dict[str, Any]:
        """Return role-projected member or administrator view of an offline proposal."""

        actor = self._require_actor(actor_user_id, require_admin=False)
        proposal = self._require_proposal(proposal_id, now=now, for_adoption=False)
        is_admin = actor.role == "administrator"
        adopted_receipt = self._proposal_adoptions.get(proposal_id)
        if adopted_receipt is not None:
            proposal_state = "adopted"
        elif now >= proposal.expires_at:
            proposal_state = "expired"
        else:
            proposal_state = "pending"

        rendered_items: list[dict[str, Any]] = []
        for item in proposal.items:
            head_rev_num, head_rev_id, _ = self._head_revision_info(
                item, proposal.project_id
            )
            before_dto, after_dto = self._render_item_diff(
                item, include_candidates=is_admin, for_write=False
            )
            rendered_items.append(
                {
                    "item_id": item.item_id,
                    "target_kind": item.target_kind,
                    "target_id": item.target_id,
                    "base_revision_id": item.base_revision_id,
                    "base_revision_number": item.base_revision_number,
                    "current_revision_id": head_rev_id,
                    "current_revision_number": head_rev_num,
                    "allowed_fields": list(item.allowed_fields),
                    "before": before_dto,
                    "after": after_dto,
                    "citations": item.citation.to_object(
                        include_candidates=is_admin,
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
            "proposal_state": proposal_state,
            "adoptable_by_caller": is_admin and proposal_state == "pending",
            "rule_version": proposal.rule_version,
            "vocabulary_version": proposal.vocabulary_version,
            "projection_version": proposal.projection_version,
            "proposal_digest": proposal.proposal_digest,
            "created_at": _format_utc(proposal.created_at),
            "expires_at": _format_utc(proposal.expires_at),
            "retention_until": _format_utc(proposal.retention_until),
            "adoption_receipt": deepcopy(adopted_receipt),
            "items": rendered_items,
        }

    def manual_revise_experiment(
        self,
        *,
        actor_user_id: str,
        experiment_id: str,
        expected_revision_number: int,
        idempotency_key: str | None,
        updates: Mapping[str, Any],
        now: datetime,
    ) -> dict[str, Any]:
        """Simulate an M1 A manual experiment revision competing with an M2 preview."""

        actor = self._require_actor(actor_user_id, require_admin=True)
        validated_key = _validate_idempotency_key(idempotency_key)
        _check_id("experiment", experiment_id)
        request_sha = _sha256_canonical(
            {
                "experiment_id": experiment_id,
                "expected_revision_number": expected_revision_number,
                "updates": dict(updates),
            }
        )
        replayed = self._lookup_idempotency(
            actor_user_id=actor.user_id,
            idempotency_key=validated_key,
            operation="catalogReviseExperiment",
            target_key=experiment_id,
            request_sha256=request_sha,
        )
        if replayed is not None:
            return replayed

        state = self._experiments.get(experiment_id)
        if state is None:
            raise catalog_preview_error(
                404, "catalog_not_found", "The experiment was not found."
            )
        self._require_active_project(state.project_id)
        if state.archived:
            raise catalog_preview_error(
                409, "catalog_experiment_archived", "The experiment is archived."
            )
        head = state.head
        if int(head["revision_number"]) != expected_revision_number:
            raise catalog_preview_error(
                409,
                "catalog_revision_conflict",
                "The catalog revision changed before this write.",
                conflicts=[
                    RevisionConflictDetail(
                        target_kind="experiment",
                        target_id=experiment_id,
                        expected_revision_number=expected_revision_number,
                        current_revision_number=int(head["revision_number"]),
                        current_revision_id=str(head["experiment_revision_id"]),
                    )
                ],
            )
        if not set(updates).issubset(EXPERIMENT_ALLOWED_CURATION_FIELDS):
            raise catalog_preview_error(
                422, "catalog_payload_invalid", "Unsupported experiment field update."
            )
        normalized_updates = {
            key: self._stored_field(
                key,
                self._render_field(
                    key,
                    value,
                    project_id=state.project_id,
                    include_candidates=True,
                    for_write=True,
                ),
            )
            for key, value in updates.items()
        }
        new_payload = deepcopy(head["payload"])
        new_payload.update(normalized_updates)
        new_rev_num = int(head["revision_number"]) + 1
        new_rev_id = self._revision_id_factory("expr_")
        state.revisions.append(
            {
                "experiment_revision_id": new_rev_id,
                "revision_number": new_rev_num,
                "updated_at": _format_utc(now),
                "updated_by": actor.user_id,
                "payload": new_payload,
            }
        )
        result = {
            "ok": True,
            "replayed": False,
            "experiment_id": experiment_id,
            "experiment_revision_id": new_rev_id,
            "revision_number": new_rev_num,
        }
        self._idempotency_receipts[(actor.user_id, validated_key)] = (
            _IdempotencyReceipt(
                operation="catalogReviseExperiment",
                target_key=experiment_id,
                request_sha256=request_sha,
                response=deepcopy(result),
            )
        )
        return result

    def preview_proposal_subset(
        self,
        *,
        actor_user_id: str,
        proposal_id: str,
        selected_item_ids: Sequence[str],
        now: datetime,
    ) -> dict[str, Any]:
        """Generate an immutable subset preview and digest without mutating catalog state."""

        self._require_actor(actor_user_id, require_admin=True)
        _assert_body_size_within_limit(
            {"proposal_id": proposal_id, "selected_item_ids": list(selected_item_ids)}
        )
        proposal = self._require_proposal(proposal_id, now=now, for_adoption=True)
        self._verify_proposal_versions(proposal)
        selected_items = self._resolve_selected_items(proposal, selected_item_ids)
        canonical_selected_ids = [item.item_id for item in selected_items]
        excluded_ids = [
            item.item_id
            for item in proposal.items
            if item.item_id not in canonical_selected_ids
        ]

        expected_revisions: dict[str, int] = {}
        preview_items: list[dict[str, Any]] = []
        conflicts: list[dict[str, Any]] = []
        for item in selected_items:
            head_rev_num, head_rev_id, archived = self._head_revision_info(
                item, proposal.project_id
            )
            if archived:
                raise catalog_preview_error(
                    409,
                    "catalog_experiment_archived",
                    "Selected target experiment is archived.",
                    path=item.item_id,
                )
            self._verify_source_references(item, proposal.project_id)
            before_dto, after_dto = self._render_item_diff(
                item, include_candidates=True, for_write=True
            )
            expected_revisions[item.target_id] = item.base_revision_number
            has_conflict = (
                head_rev_num != item.base_revision_number
                or head_rev_id != item.base_revision_id
            )
            if has_conflict:
                conflicts.append(
                    RevisionConflictDetail(
                        target_kind=item.target_kind,
                        target_id=item.target_id,
                        expected_revision_number=item.base_revision_number,
                        current_revision_number=head_rev_num,
                        current_revision_id=head_rev_id,
                        expected_revision_id=item.base_revision_id,
                    ).to_object()
                )
            preview_items.append(
                {
                    "item_id": item.item_id,
                    "target_kind": item.target_kind,
                    "target_id": item.target_id,
                    "expected_revision_number": item.base_revision_number,
                    "base_revision_id": item.base_revision_id,
                    "current_revision_number": head_rev_num,
                    "current_revision_id": head_rev_id,
                    "conflict_detected": has_conflict,
                    "allowed_fields": list(item.allowed_fields),
                    "before": before_dto,
                    "after": after_dto,
                    "citations": item.citation.to_object(include_candidates=True),
                }
            )

        preview_digest = compute_preview_digest(
            proposal_id=proposal.proposal_id,
            proposal_digest=proposal.proposal_digest,
            selected_item_ids=canonical_selected_ids,
            expected_revisions=expected_revisions,
        )
        return {
            "ok": True,
            "contract_version": M2_CONTRACT_VERSION,
            "proposal_id": proposal.proposal_id,
            "project_id": proposal.project_id,
            "proposal_digest": proposal.proposal_digest,
            "preview_digest": preview_digest,
            "selected_item_ids": canonical_selected_ids,
            "excluded_item_ids": excluded_ids,
            "expected_revisions": expected_revisions,
            "has_conflicts": bool(conflicts),
            "conflicts": conflicts,
            "expires_at": _format_utc(proposal.expires_at),
            "items": preview_items,
        }

    def adopt_proposal_subset(
        self,
        *,
        actor_user_id: str,
        proposal_id: str,
        preview_digest: str,
        selected_item_ids: Sequence[str],
        expected_revisions: Mapping[str, int],
        idempotency_key: str | None,
        now: datetime,
    ) -> dict[str, Any]:
        """Atomically adopt a previewed subset with strict baseline, provenance, and idempotency checks."""

        validated_key = _validate_idempotency_key(idempotency_key)
        actor = self._require_actor(actor_user_id, require_admin=True)
        canonical_body = {
            "proposal_id": proposal_id,
            "preview_digest": preview_digest,
            "selected_item_ids": list(selected_item_ids),
            "expected_revisions": {
                key: expected_revisions[key] for key in sorted(expected_revisions)
            },
        }
        _assert_body_size_within_limit(canonical_body)
        if (
            not isinstance(preview_digest, str)
            or _DIGEST_RE.fullmatch(preview_digest) is None
        ):
            raise catalog_preview_error(
                422, "catalog_payload_invalid", "Invalid preview_digest format."
            )

        request_sha = _sha256_canonical(canonical_body)
        replayed = self._lookup_idempotency(
            actor_user_id=actor.user_id,
            idempotency_key=validated_key,
            operation="catalogAdoptProposal",
            target_key=proposal_id,
            request_sha256=request_sha,
        )
        if replayed is not None:
            return replayed

        proposal = self._require_proposal(proposal_id, now=now, for_adoption=True)
        self._verify_proposal_versions(proposal)
        selected_items = self._resolve_selected_items(proposal, selected_item_ids)
        canonical_selected_ids = [item.item_id for item in selected_items]
        if list(selected_item_ids) != canonical_selected_ids:
            recomputed = compute_preview_digest(
                proposal_id=proposal.proposal_id,
                proposal_digest=proposal.proposal_digest,
                selected_item_ids=selected_item_ids,
                expected_revisions=expected_revisions,
            )
        else:
            recomputed = compute_preview_digest(
                proposal_id=proposal.proposal_id,
                proposal_digest=proposal.proposal_digest,
                selected_item_ids=canonical_selected_ids,
                expected_revisions=expected_revisions,
            )
        if recomputed != preview_digest:
            raise catalog_preview_error(
                409,
                "catalog_preview_digest_mismatch",
                "Submitted selection or expected revisions do not match preview_digest.",
            )

        selected_target_ids = {item.target_id for item in selected_items}
        if set(expected_revisions) != selected_target_ids:
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                "expected_revisions must match the selected proposal targets.",
            )

        # Verify target ownership, source references, slot pins, AND both
        # base_revision_id and base_revision_number against server state before writing.
        conflicts: list[RevisionConflictDetail] = []
        for item in selected_items:
            caller_expected = expected_revisions[item.target_id]
            if isinstance(caller_expected, bool) or not isinstance(
                caller_expected, int
            ):
                raise catalog_preview_error(
                    422,
                    "catalog_payload_invalid",
                    "expected_revision_number must be an integer.",
                )
            if caller_expected != item.base_revision_number:
                raise catalog_preview_error(
                    409,
                    "catalog_proposal_invalidated",
                    "Caller expected revision does not match the server proposal baseline.",
                    path=item.item_id,
                )
            head_rev_num, head_rev_id, archived = self._head_revision_info(
                item, proposal.project_id
            )
            if archived:
                raise catalog_preview_error(
                    409,
                    "catalog_experiment_archived",
                    "Selected target experiment is archived.",
                    path=item.target_id,
                )
            self._verify_source_references(item, proposal.project_id)
            self._render_item_diff(item, include_candidates=True, for_write=True)
            if head_rev_num != item.base_revision_number or (
                head_rev_id != item.base_revision_id
            ):
                conflicts.append(
                    RevisionConflictDetail(
                        target_kind=item.target_kind,
                        target_id=item.target_id,
                        expected_revision_number=item.base_revision_number,
                        current_revision_number=head_rev_num,
                        current_revision_id=head_rev_id,
                        expected_revision_id=item.base_revision_id,
                    )
                )

        if conflicts:
            raise catalog_preview_error(
                409,
                "catalog_revision_conflict",
                "One or more selected catalog targets changed before adoption.",
                conflicts=conflicts,
            )

        # Atomically commit ONLY selected_items, proposal adoption state, and receipt.
        applied_targets: list[dict[str, Any]] = []
        for item in selected_items:
            applied_targets.append(
                self._apply_item(item, actor_user_id=actor.user_id, now=now)
            )

        excluded_ids = [
            item.item_id
            for item in proposal.items
            if item.item_id not in canonical_selected_ids
        ]
        receipt = {
            "ok": True,
            "replayed": False,
            "contract_version": M2_CONTRACT_VERSION,
            "receipt_id": self._revision_id_factory("cadopt_"),
            "proposal_id": proposal.proposal_id,
            "project_id": proposal.project_id,
            "proposal_digest": proposal.proposal_digest,
            "preview_digest": preview_digest,
            "actor_user_id": actor.user_id,
            "idempotency_key": validated_key,
            "selected_item_ids": canonical_selected_ids,
            "excluded_item_ids": excluded_ids,
            "applied_targets": applied_targets,
            "adopted_at": _format_utc(now),
        }
        self._proposal_adoptions[proposal.proposal_id] = deepcopy(receipt)
        self._idempotency_receipts[(actor.user_id, validated_key)] = (
            _IdempotencyReceipt(
                operation="catalogAdoptProposal",
                target_key=proposal.proposal_id,
                request_sha256=request_sha,
                response=deepcopy(receipt),
            )
        )
        self._adoption_receipts.append(deepcopy(receipt))
        return receipt

    def list_adoption_receipts(
        self,
        *,
        actor_user_id: str,
        project_id: str,
        proposal_id: str | None = None,
    ) -> dict[str, Any]:
        """Discover previously committed adoption receipts after dropped response or re-login."""

        self._require_actor(actor_user_id, require_admin=True)
        self._require_active_project(project_id)
        if proposal_id is not None:
            _check_id("proposal", proposal_id, path="proposal_id")
        matched = [
            deepcopy(receipt)
            for receipt in self._adoption_receipts
            if receipt["project_id"] == project_id
            and (proposal_id is None or receipt["proposal_id"] == proposal_id)
        ]
        return {
            "ok": True,
            "project_id": project_id,
            "receipts": matched,
        }

    def _resolve_selected_items(
        self, proposal: CatalogProposal, selected_item_ids: Sequence[str]
    ) -> tuple[ProposalItem, ...]:
        if (
            not isinstance(selected_item_ids, Sequence)
            or isinstance(selected_item_ids, (str, bytes))
            or not selected_item_ids
        ):
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                "At least one proposal item must be selected for preview or adoption.",
            )
        if len(selected_item_ids) > MAX_CHANGESET_ITEMS:
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                f"Selection exceeds maximum of {MAX_CHANGESET_ITEMS} items.",
            )
        if len(set(selected_item_ids)) != len(selected_item_ids):
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                "Duplicate item_id values in selection.",
            )
        by_id = {item.item_id: item for item in proposal.items}
        selected_set = set(selected_item_ids)
        unknown = [item_id for item_id in selected_item_ids if item_id not in by_id]
        if unknown:
            raise catalog_preview_error(
                422,
                "catalog_payload_invalid",
                "Selection references an unknown proposal item_id.",
                path="selected_item_ids",
            )
        # Preserve canonical proposal item ordering.
        return tuple(item for item in proposal.items if item.item_id in selected_set)

    def _apply_item(
        self, item: ProposalItem, *, actor_user_id: str, now: datetime
    ) -> dict[str, Any]:
        _, after_rendered = self._render_item_diff(
            item, include_candidates=True, for_write=True
        )
        if item.target_kind == "experiment":
            state = self._experiments[item.target_id]
            head = state.head
            new_payload = deepcopy(head["payload"])
            for field_name in item.allowed_fields:
                new_payload[field_name] = self._stored_field(
                    field_name, after_rendered[field_name]
                )
            new_rev_num = int(head["revision_number"]) + 1
            new_rev_id = self._revision_id_factory("expr_")
            state.revisions.append(
                {
                    "experiment_revision_id": new_rev_id,
                    "revision_number": new_rev_num,
                    "updated_at": _format_utc(now),
                    "updated_by": actor_user_id,
                    "payload": new_payload,
                }
            )
            return {
                "item_id": item.item_id,
                "target_kind": "experiment",
                "target_id": item.target_id,
                "base_revision_id": head["experiment_revision_id"],
                "base_revision_number": int(head["revision_number"]),
                "result_revision_id": new_rev_id,
                "result_revision_number": new_rev_num,
                "changed_fields": list(item.allowed_fields),
            }

        sample_state = self._samples[item.target_id]
        head = sample_state.head
        new_catalog = deepcopy(head["catalog"])
        new_display_name = head["display_name"]
        for field_name in item.allowed_fields:
            stored = self._stored_field(field_name, after_rendered[field_name])
            if field_name == "display_name":
                new_display_name = stored
            else:
                new_catalog[field_name] = stored
        new_rev_num = int(head["revision_number"]) + 1
        new_rev_id = self._revision_id_factory("smpr_")
        sample_state.revisions.append(
            {
                "sample_revision_id": new_rev_id,
                "revision_number": new_rev_num,
                "display_name": new_display_name,
                "updated_at": _format_utc(now),
                "updated_by": actor_user_id,
                "catalog": new_catalog,
                # Non-catalog attributes must be preserved verbatim across revisions.
                "extra_attributes": deepcopy(head["extra_attributes"]),
            }
        )
        return {
            "item_id": item.item_id,
            "target_kind": "sample_catalog",
            "target_id": item.target_id,
            "base_revision_id": head["sample_revision_id"],
            "base_revision_number": int(head["revision_number"]),
            "result_revision_id": new_rev_id,
            "result_revision_number": new_rev_num,
            "changed_fields": list(item.allowed_fields),
        }


SYNTHETIC_PROJECT_ID = "prj_11111111111111111111111111111111"
SYNTHETIC_OTHER_PROJECT_ID = "prj_22222222222222222222222222222222"
SYNTHETIC_PROPOSAL_ID = "cprop_0123456789abcdef0123456789abcdef"
SYNTHETIC_SECOND_PROPOSAL_ID = "cprop_22222222222222222222222222222222"
SYNTHETIC_EXP_1_ID = "exp_44444444444444444444444444444444"
SYNTHETIC_EXP_2_ID = "exp_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
SYNTHETIC_SMP_1_ID = "smp_22222222222222222222222222222222"
SYNTHETIC_SMP_CTRL_ID = "smp_abababababababababababababababab"
SYNTHETIC_BASE_TIME = datetime(2026, 11, 18, 10, 0, 0, tzinfo=timezone.utc)


def build_synthetic_m2_preview_fixture() -> SyntheticCatalogPreviewStore:
    """Build the canonical synthetic catalog table and offline proposal fixture."""

    terms = (
        TermRevisionRecord(
            term_id="term_eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee",
            term_revision_id="trmv_ffffffffffffffffffffffffffffffff",
            term_type="assay",
            display_name="CUT&Tag",
        ),
        TermRevisionRecord(
            term_id="term_12121212121212121212121212121212",
            term_revision_id="trmv_34343434343434343434343434343434",
            term_type="target",
            display_name="H3K27me3",
        ),
        TermRevisionRecord(
            term_id="term_23232323232323232323232323232323",
            term_revision_id="trmv_45454545454545454545454545454545",
            term_type="target",
            display_name="H3K4me3",
        ),
        TermRevisionRecord(
            term_id="term_34343434343434343434343434343434",
            term_revision_id="trmv_56565656565656565656565656565656",
            term_type="batch",
            display_name="batch-2026-11",
        ),
        TermRevisionRecord(
            term_id="term_45454545454545454545454545454545",
            term_revision_id="trmv_67676767676767676767676767676767",
            term_type="species",
            display_name="Mus musculus",
        ),
        TermRevisionRecord(
            term_id="term_56565656565656565656565656565656",
            term_revision_id="trmv_78787878787878787878787878787878",
            term_type="tissue",
            display_name="lymph node",
        ),
        TermRevisionRecord(
            term_id="term_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
            term_revision_id="trmv_bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
            term_type="cell_type",
            display_name="regulatory T cell",
            aliases=("Treg", "CD4+CD25+ T cell"),
        ),
        TermRevisionRecord(
            term_id="term_cccccccccccccccccccccccccccccccc",
            term_revision_id="trmv_dddddddddddddddddddddddddddddddd",
            term_type="cell_type",
            display_name="keratinocyte",
            aliases=("KC",),
        ),
        TermRevisionRecord(
            term_id="term_67676767676767676767676767676767",
            term_revision_id="trmv_89898989898989898989898989898989",
            term_type="replicate_role",
            display_name="treatment",
        ),
        TermRevisionRecord(
            term_id="term_78787878787878787878787878787878",
            term_revision_id="trmv_90909090909090909090909090909090",
            term_type="replicate_role",
            display_name="control",
        ),
        TermRevisionRecord(
            term_id="term_89898989898989898989898989898989",
            term_revision_id="trmv_10101010101010101010101010101010",
            term_type="tag",
            display_name="m2-curated",
        ),
    )
    terms_by_revision = {term.term_revision_id: term for term in terms}

    exp_1 = _ExperimentState(
        experiment_id=SYNTHETIC_EXP_1_ID,
        project_id=SYNTHETIC_PROJECT_ID,
        stable_key="treg-h3k27me3-cuttag",
        archived=False,
        revisions=[
            {
                "experiment_revision_id": "expr_44444444444444444444444444444441",
                "revision_number": 1,
                "payload": {
                    "schema_version": 1,
                    "display_name": "Treg H3K27me3 CUT&Tag",
                    "assay": {
                        "status": "known",
                        "term_id": "term_eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee",
                        "term_revision_id": "trmv_ffffffffffffffffffffffffffffffff",
                    },
                    "target": {
                        "status": "known",
                        "term_id": "term_12121212121212121212121212121212",
                        "term_revision_id": "trmv_34343434343434343434343434343434",
                    },
                    "batch": {"status": "unknown"},
                    "samples": [
                        {
                            "ordinal": 0,
                            "sample_revision_id": "smpr_33333333333333333333333333333331",
                            "replicate_group": "biorep-1",
                            "replicate_role": {
                                "status": "known",
                                "term_id": "term_67676767676767676767676767676767",
                                "term_revision_id": "trmv_89898989898989898989898989898989",
                            },
                        },
                        {
                            "ordinal": 1,
                            "sample_revision_id": "smpr_abababababababababababababababab",
                            "replicate_group": "biorep-1",
                            "replicate_role": {
                                "status": "known",
                                "term_id": "term_78787878787878787878787878787878",
                                "term_revision_id": "trmv_90909090909090909090909090909090",
                            },
                        },
                    ],
                },
            }
        ],
    )
    exp_2 = _ExperimentState(
        experiment_id=SYNTHETIC_EXP_2_ID,
        project_id=SYNTHETIC_PROJECT_ID,
        stable_key="treg-h3k4me3-cuttag",
        archived=False,
        revisions=[
            {
                "experiment_revision_id": "expr_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa1",
                "revision_number": 1,
                "payload": {
                    "schema_version": 1,
                    "display_name": "Treg H3K4me3 CUT&Tag",
                    "assay": {
                        "status": "known",
                        "term_id": "term_eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee",
                        "term_revision_id": "trmv_ffffffffffffffffffffffffffffffff",
                    },
                    "target": {"status": "unknown"},
                    "batch": {
                        "status": "known",
                        "term_id": "term_34343434343434343434343434343434",
                        "term_revision_id": "trmv_56565656565656565656565656565656",
                    },
                    "samples": [
                        {
                            "ordinal": 0,
                            "sample_revision_id": "smpr_33333333333333333333333333333331",
                            "replicate_group": "biorep-1",
                            "replicate_role": {
                                "status": "known",
                                "term_id": "term_67676767676767676767676767676767",
                                "term_revision_id": "trmv_89898989898989898989898989898989",
                            },
                        }
                    ],
                },
            }
        ],
    )

    sample_1 = _SampleState(
        sample_id=SYNTHETIC_SMP_1_ID,
        project_id=SYNTHETIC_PROJECT_ID,
        stable_key="treg-sample-1",
        revisions=[
            {
                "sample_revision_id": "smpr_33333333333333333333333333333331",
                "revision_number": 1,
                "display_name": "Treg sample 1",
                "catalog": {
                    "schema_version": 1,
                    "species": {
                        "status": "known",
                        "term_id": "term_45454545454545454545454545454545",
                        "term_revision_id": "trmv_67676767676767676767676767676767",
                    },
                    "tissue": {
                        "status": "known",
                        "term_id": "term_56565656565656565656565656565656",
                        "term_revision_id": "trmv_78787878787878787878787878787878",
                    },
                    "cell_type": {"status": "unknown", "values": []},
                    "cell_line": {"status": "not_applicable", "values": []},
                    "treatment": {"status": "not_detected", "values": []},
                    "tag": {"status": "unknown", "values": []},
                },
                "extra_attributes": {"legacy_note": "preserve-me-outside-catalog"},
            }
        ],
    )

    item_1 = ProposalItem(
        item_id="item_1",
        target_kind="experiment",
        target_id=SYNTHETIC_EXP_1_ID,
        base_revision_id="expr_44444444444444444444444444444441",
        base_revision_number=1,
        allowed_fields=("batch",),
        before={"batch": {"status": "unknown"}},
        after={
            "batch": {
                "status": "known",
                "term_id": "term_34343434343434343434343434343434",
                "term_revision_id": "trmv_56565656565656565656565656565656",
            }
        },
        citation=ProposalCitation(
            rule_version=M2_RULE_VERSION,
            vocabulary_version=M2_VOCABULARY_VERSION,
            reason_codes=("BATCH_TERM_MATCH",),
            source_refs=(
                SourceReference(
                    kind="experiment_revision",
                    ref_id="expr_44444444444444444444444444444441",
                    evidence_token="proj-exp-1-batch-slot",
                ),
                SourceReference(
                    kind="term_revision",
                    ref_id="trmv_56565656565656565656565656565656",
                    term_revision_id="trmv_56565656565656565656565656565656",
                    evidence_token="vocab-batch-2026-11",
                ),
            ),
        ),
    )
    item_2 = ProposalItem(
        item_id="item_2",
        target_kind="experiment",
        target_id=SYNTHETIC_EXP_2_ID,
        base_revision_id="expr_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa1",
        base_revision_number=1,
        allowed_fields=("target",),
        before={"target": {"status": "unknown"}},
        after={
            "target": {
                "status": "known",
                "term_id": "term_23232323232323232323232323232323",
                "term_revision_id": "trmv_45454545454545454545454545454545",
            }
        },
        citation=ProposalCitation(
            rule_version=M2_RULE_VERSION,
            vocabulary_version=M2_VOCABULARY_VERSION,
            reason_codes=("TARGET_TERM_MATCH",),
            source_refs=(
                SourceReference(
                    kind="experiment_revision",
                    ref_id="expr_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa1",
                    evidence_token="proj-exp-2-target-slot",
                ),
                SourceReference(
                    kind="term_revision",
                    ref_id="trmv_45454545454545454545454545454545",
                    term_revision_id="trmv_45454545454545454545454545454545",
                    evidence_token="vocab-h3k4me3",
                ),
            ),
        ),
    )
    item_3 = ProposalItem(
        item_id="item_3",
        target_kind="sample_catalog",
        target_id=SYNTHETIC_SMP_1_ID,
        base_revision_id="smpr_33333333333333333333333333333331",
        base_revision_number=1,
        allowed_fields=("cell_type", "tag"),
        before={
            "cell_type": {"status": "unknown", "values": []},
            "tag": {"status": "unknown", "values": []},
        },
        after={
            "cell_type": {
                "status": "ambiguous",
                "values": [],
                "candidates": [
                    {
                        "term_id": "term_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
                        "term_revision_id": "trmv_bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
                    },
                    {
                        "term_id": "term_cccccccccccccccccccccccccccccccc",
                        "term_revision_id": "trmv_dddddddddddddddddddddddddddddddd",
                    },
                ],
            },
            "tag": {
                "status": "known",
                "values": [
                    {
                        "term_id": "term_89898989898989898989898989898989",
                        "term_revision_id": "trmv_10101010101010101010101010101010",
                    }
                ],
            },
        },
        citation=ProposalCitation(
            rule_version=M2_RULE_VERSION,
            vocabulary_version=M2_VOCABULARY_VERSION,
            reason_codes=("ALIAS_AMBIGUOUS_CANDIDATES", "TAG_TERM_MATCH"),
            source_refs=(
                SourceReference(
                    kind="sample_revision",
                    ref_id="smpr_33333333333333333333333333333331",
                    evidence_token="proj-smp-1-catalog-slots",
                ),
                SourceReference(
                    kind="term_revision",
                    ref_id="trmv_bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
                    term_revision_id="trmv_bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
                    evidence_token="alias-candidate-treg",
                    visibility="admin_only",
                ),
                SourceReference(
                    kind="term_revision",
                    ref_id="trmv_dddddddddddddddddddddddddddddddd",
                    term_revision_id="trmv_dddddddddddddddddddddddddddddddd",
                    evidence_token="alias-candidate-kc",
                    visibility="admin_only",
                ),
                SourceReference(
                    kind="term_revision",
                    ref_id="trmv_10101010101010101010101010101010",
                    term_revision_id="trmv_10101010101010101010101010101010",
                    evidence_token="vocab-m2-curated",
                ),
            ),
        ),
    )

    proposal = CatalogProposal.create(
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        project_id=SYNTHETIC_PROJECT_ID,
        created_by_user_id="usr_member_01",
        source="offline_suggestion",
        items=(item_1, item_2, item_3),
        created_at=SYNTHETIC_BASE_TIME,
    )
    second_proposal = CatalogProposal.create(
        proposal_id=SYNTHETIC_SECOND_PROPOSAL_ID,
        project_id=SYNTHETIC_PROJECT_ID,
        created_by_user_id="usr_member_01",
        source="offline_suggestion",
        items=(item_2,),
        created_at=SYNTHETIC_BASE_TIME + timedelta(minutes=1),
    )

    return SyntheticCatalogPreviewStore(
        projects={
            SYNTHETIC_PROJECT_ID: {
                "active": True,
                "legacy": False,
                "display_name": "Treg 表观 (M2 合成目录夹具)",
            },
            SYNTHETIC_OTHER_PROJECT_ID: {
                "active": True,
                "legacy": False,
                "display_name": "Other Project (Cross-Project Isolation Check)",
            },
        },
        actors={
            "usr_member_01": CatalogActor(user_id="usr_member_01", role="member"),
            "usr_admin_a": CatalogActor(user_id="usr_admin_a", role="administrator"),
            "usr_admin_b": CatalogActor(user_id="usr_admin_b", role="administrator"),
        },
        terms_by_revision=terms_by_revision,
        experiments={SYNTHETIC_EXP_1_ID: exp_1, SYNTHETIC_EXP_2_ID: exp_2},
        samples={SYNTHETIC_SMP_1_ID: sample_1},
        proposals={
            proposal.proposal_id: proposal,
            second_proposal.proposal_id: second_proposal,
        },
    )


def export_m2_preview_contract_examples() -> dict[str, Any]:
    """Run the synthetic M2 fixture scenarios and return canonical contract examples."""

    t_preview = SYNTHETIC_BASE_TIME + timedelta(minutes=2)
    t_concurrent = SYNTHETIC_BASE_TIME + timedelta(minutes=3)
    t_adopt = SYNTHETIC_BASE_TIME + timedelta(minutes=4)
    t_retry = SYNTHETIC_BASE_TIME + timedelta(minutes=5)

    # Run 1: Role-projected views, subset adoption, single-outcome guard, replay, and M1 A idempotency scope.
    clean_store = build_synthetic_m2_preview_fixture()
    member_view = clean_store.get_proposal_view(
        actor_user_id="usr_member_01",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        now=t_preview,
    )
    admin_view = clean_store.get_proposal_view(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        now=t_preview,
    )
    admin_subset_preview = clean_store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        selected_item_ids=["item_1", "item_3"],
        now=t_preview,
    )
    disjoint_preview_admin_b = clean_store.preview_proposal_subset(
        actor_user_id="usr_admin_b",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        selected_item_ids=["item_2"],
        now=t_preview,
    )
    second_proposal_preview = clean_store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_SECOND_PROPOSAL_ID,
        selected_item_ids=["item_2"],
        now=t_preview,
    )
    adopt_receipt = clean_store.adopt_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        preview_digest=admin_subset_preview["preview_digest"],
        selected_item_ids=admin_subset_preview["selected_item_ids"],
        expected_revisions=admin_subset_preview["expected_revisions"],
        idempotency_key="idem-admin-a-adopt-1",
        now=t_adopt,
    )
    reopened_adopted_view = clean_store.get_proposal_view(
        actor_user_id="usr_admin_b",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        now=t_retry,
    )
    try:
        clean_store.adopt_proposal_subset(
            actor_user_id="usr_admin_b",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            preview_digest=disjoint_preview_admin_b["preview_digest"],
            selected_item_ids=disjoint_preview_admin_b["selected_item_ids"],
            expected_revisions=disjoint_preview_admin_b["expected_revisions"],
            idempotency_key="idem-admin-b-disjoint",
            now=t_retry,
        )
        raise AssertionError(
            "expected proposal already adopted error"
        )  # pragma: no cover
    except CatalogPreviewError as exc:
        disjoint_subset_after_adoption_conflict = exc.to_response()

    dropped_response_replay = clean_store.adopt_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        preview_digest=admin_subset_preview["preview_digest"],
        selected_item_ids=admin_subset_preview["selected_item_ids"],
        expected_revisions=admin_subset_preview["expected_revisions"],
        idempotency_key="idem-admin-a-adopt-1",
        now=t_retry,
    )
    try:
        clean_store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=SYNTHETIC_SECOND_PROPOSAL_ID,
            preview_digest=second_proposal_preview["preview_digest"],
            selected_item_ids=second_proposal_preview["selected_item_ids"],
            expected_revisions=second_proposal_preview["expected_revisions"],
            idempotency_key="idem-admin-a-adopt-1",
            now=t_retry,
        )
        raise AssertionError(
            "expected cross-proposal idempotency conflict"
        )  # pragma: no cover
    except CatalogPreviewError as exc:
        idempotency_conflict_response = exc.to_response()

    # Run 2: Two concurrent expected revisions producing a 409 atomic rollback.
    conflict_store = build_synthetic_m2_preview_fixture()
    preview_before_conflict = conflict_store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        selected_item_ids=["item_1", "item_3"],
        now=t_preview,
    )
    conflict_store.manual_revise_experiment(
        actor_user_id="usr_admin_b",
        experiment_id=SYNTHETIC_EXP_1_ID,
        expected_revision_number=1,
        idempotency_key="idem-admin-b-manual-1",
        updates={"display_name": "Treg H3K27me3 CUT&Tag (Manual Rev 2)"},
        now=t_concurrent,
    )
    try:
        conflict_store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            preview_digest=preview_before_conflict["preview_digest"],
            selected_item_ids=preview_before_conflict["selected_item_ids"],
            expected_revisions=preview_before_conflict["expected_revisions"],
            idempotency_key="idem-admin-a-conflict-1",
            now=t_adopt,
        )
        raise AssertionError("expected revision conflict")  # pragma: no cover
    except CatalogPreviewError as exc:
        concurrent_revision_conflict_response = exc.to_response()

    return {
        "contract_version": M2_CONTRACT_VERSION,
        "member_proposal_view": member_view,
        "admin_proposal_view": admin_view,
        "admin_subset_preview": admin_subset_preview,
        "concurrent_revision_conflict_response": concurrent_revision_conflict_response,
        "adopt_success_receipt": adopt_receipt,
        "reopened_adopted_proposal_view": reopened_adopted_view,
        "disjoint_subset_after_adoption_conflict": disjoint_subset_after_adoption_conflict,
        "adopt_dropped_response_replay": dropped_response_replay,
        "idempotency_conflict_response": idempotency_conflict_response,
    }
