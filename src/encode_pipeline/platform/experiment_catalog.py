"""Experiment catalog identities, payloads, cursors, and public projections.

Slice A freezes manual vocabulary and experiment maintenance. Slice B freezes an
experiment revision onto a validation snapshot only through
``freeze_experiment_snapshot_binding``; persistence copies that pair onto the run.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Protocol
from datetime import datetime, timezone
from hashlib import sha256
import base64
import binascii
import json
import re
from typing import Any


EXPERIMENT_REVISION_PAYLOAD_DIGEST_SCHEME = (
    "sha256-framed-experiment-revision-payload-v1"
)
CATALOG_CURSOR_PREFIX = "catcur_"
CATALOG_CURSOR_MAX_LENGTH = 1024
CATALOG_SCHEMA_VERSION = 1
DEFAULT_PAGE_SIZE = 20
MAX_PAGE_SIZE = 50

TERM_TYPES = (
    "species",
    "tissue",
    "cell_type",
    "cell_line",
    "treatment",
    "tag",
    "assay",
    "target",
    "batch",
    "replicate_role",
)
SAMPLE_SINGLE_FIELDS = ("species", "tissue")
SAMPLE_MULTI_FIELDS = ("cell_type", "cell_line", "treatment", "tag")
SAMPLE_FIELDS = SAMPLE_SINGLE_FIELDS + SAMPLE_MULTI_FIELDS
EXPERIMENT_FIELDS = ("assay", "target", "batch")
FIELD_STATUSES = (
    "known",
    "unknown",
    "not_applicable",
    "not_detected",
    "ambiguous",
)
NEGATIVE_STATUSES = ("unknown", "not_applicable", "not_detected")
REVISION_VIEWS = ("linked", "current_sample_catalog")

_EXPERIMENT_ID = re.compile(r"^exp_[0-9a-f]{32}$")
_EXPERIMENT_REVISION_ID = re.compile(r"^expr_[0-9a-f]{32}$")
_TERM_ID = re.compile(r"^term_[0-9a-f]{32}$")
_TERM_REVISION_ID = re.compile(r"^trmv_[0-9a-f]{32}$")
_SAMPLE_ID = re.compile(r"^smp_[0-9a-f]{32}$")
_SAMPLE_REVISION_ID = re.compile(r"^smpr_[0-9a-f]{32}$")
_PROJECT_ID = re.compile(r"^prj_[0-9a-f]{32}$")
_STABLE_KEY = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]{0,254}$")
_IDEMPOTENCY_KEY = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_CURSOR_BODY = re.compile(r"^[A-Za-z0-9_-]+$")
_CURSOR_KEYS = {"filter_sha256", "kind", "position", "v"}


class CatalogError(Exception):
    """Stable catalog failure whose message is safe to return to clients."""

    def __init__(
        self,
        status_code: int,
        code: str,
        message: str,
        *,
        path: str | None = None,
    ) -> None:
        super().__init__(message)
        self.status_code = status_code
        self.code = code
        self.message = message
        self.path = path


@dataclass(frozen=True)
class CatalogWrite:
    """One administrator write, including the raw-body idempotency receipt."""

    actor_user_id: str
    idempotency_key: str
    request_sha256: str
    operation: str
    now: datetime
    target_key: str = ""


def bind_catalog_target(write: CatalogWrite, target_key: str) -> CatalogWrite:
    """Bind one receipt to the project or resource the request names."""

    if (
        not isinstance(target_key, str)
        or not target_key
        or len(target_key) > 64
        or any(character.isspace() for character in target_key)
    ):
        raise ValueError("catalog target must identify one write scope")
    return replace(write, target_key=target_key)


class CatalogBlockRequiresCatalogCommand(ValueError):
    """Local sample edits must not author or drop a catalog block."""

    code = "catalog_block_requires_catalog_command"

    def __init__(self) -> None:
        super().__init__("Catalog block requires the catalog command.")


@dataclass(frozen=True)
class TermPin:
    """One existing term revision visible to catalog validation."""

    term_id: str
    term_revision_id: str
    term_type: str
    display_name: str
    archived: bool


@dataclass(frozen=True)
class PinnedTerm:
    term_id: str
    term_revision_id: str


@dataclass(frozen=True)
class SingleSlot:
    status: str
    term: PinnedTerm | None = None
    candidates: tuple[PinnedTerm, ...] = ()


@dataclass(frozen=True)
class MultiSlot:
    status: str
    values: tuple[PinnedTerm, ...] = ()
    candidates: tuple[PinnedTerm, ...] = ()


@dataclass(frozen=True)
class SampleCatalogBlock:
    species: SingleSlot
    tissue: SingleSlot
    cell_type: MultiSlot
    cell_line: MultiSlot
    treatment: MultiSlot
    tag: MultiSlot

    def slots(self) -> Mapping[str, SingleSlot | MultiSlot]:
        return {
            "species": self.species,
            "tissue": self.tissue,
            "cell_type": self.cell_type,
            "cell_line": self.cell_line,
            "treatment": self.treatment,
            "tag": self.tag,
        }


@dataclass(frozen=True)
class ExperimentSamplePin:
    sample_id: str
    sample_revision_id: str
    payload_digest: str
    ordinal: int
    replicate_group: str | None
    replicate_role: SingleSlot


@dataclass(frozen=True)
class ExperimentPayload:
    display_name: str
    assay: SingleSlot
    target: SingleSlot
    batch: SingleSlot
    samples: tuple[ExperimentSamplePin, ...]

    def canonical_json(self) -> str:
        return _canonical_json(self.to_object())

    def to_object(self) -> dict[str, Any]:
        samples: list[dict[str, Any]] = []
        for pin in self.samples:
            item: dict[str, Any] = {
                "ordinal": pin.ordinal,
                "payload_digest": pin.payload_digest,
                "replicate_role": _single_object(pin.replicate_role),
                "sample_revision_id": pin.sample_revision_id,
            }
            if pin.replicate_group is not None:
                item["replicate_group"] = pin.replicate_group
            samples.append(item)
        return {
            "assay": _single_object(self.assay),
            "batch": _single_object(self.batch),
            "display_name": self.display_name,
            "samples": samples,
            "schema_version": CATALOG_SCHEMA_VERSION,
            "target": _single_object(self.target),
        }


@dataclass(frozen=True)
class FieldFilter:
    name: str
    status: str | None
    term_ids: tuple[str, ...]


@dataclass(frozen=True)
class CatalogCursorPosition:
    created_at: datetime
    object_id: str


def catalog_not_found() -> CatalogError:
    return CatalogError(404, "catalog_not_found", "The catalog record was not found.")


def catalog_payload_invalid() -> CatalogError:
    return CatalogError(
        422,
        "catalog_payload_invalid",
        "The catalog payload is invalid.",
    )


def catalog_filter_invalid() -> CatalogError:
    return CatalogError(
        400,
        "catalog_filter_invalid",
        "The catalog filter is invalid.",
    )


def catalog_cursor_invalid() -> CatalogError:
    return CatalogError(
        400,
        "catalog_cursor_invalid",
        "The catalog cursor is invalid.",
    )


def catalog_cursor_not_found() -> CatalogError:
    return CatalogError(
        404,
        "catalog_cursor_not_found",
        "The catalog cursor no longer matches this query.",
    )


def catalog_revision_conflict() -> CatalogError:
    return CatalogError(
        409,
        "catalog_revision_conflict",
        "The catalog revision changed before this write.",
    )


def catalog_idempotency_conflict() -> CatalogError:
    return CatalogError(
        409,
        "catalog_idempotency_conflict",
        "The idempotency key was already used for a different request.",
    )


def catalog_already_archived() -> CatalogError:
    return CatalogError(
        409,
        "catalog_already_archived",
        "The catalog record is already archived.",
    )


def catalog_term_archived() -> CatalogError:
    return CatalogError(
        409,
        "catalog_term_archived",
        "The catalog term is archived.",
    )


def catalog_experiment_archived() -> CatalogError:
    return CatalogError(
        409,
        "catalog_experiment_archived",
        "The experiment is archived.",
    )


def catalog_sample_set_mismatch() -> CatalogError:
    return CatalogError(
        409,
        "catalog_sample_set_mismatch",
        "The sample revision set does not match the experiment revision.",
    )


@dataclass(frozen=True)
class ExperimentSnapshotBinding:
    """Experiment revision id and payload digest frozen onto one snapshot or run."""

    experiment_revision_id: str
    payload_digest: str

    def __post_init__(self) -> None:
        if _EXPERIMENT_REVISION_ID.fullmatch(self.experiment_revision_id) is None:
            raise ValueError("experiment_revision_id is invalid")
        if _DIGEST.fullmatch(self.payload_digest) is None:
            raise ValueError("payload_digest is invalid")


@dataclass(frozen=True)
class ExperimentRevisionFreezeSource:
    """Read-only view used to decide whether a snapshot may freeze one revision."""

    experiment_revision_id: str
    project_id: str
    project_active: bool
    experiment_archived: bool
    payload_digest: str
    sample_revision_ids: tuple[str, ...]
    association_matches_payload: bool


class ExperimentRevisionFreezeLookup(Protocol):
    """In-memory catalog view consulted inside the validation transaction."""

    def experiment_revision_freeze_source(
        self,
        experiment_revision_id: str,
    ) -> ExperimentRevisionFreezeSource | None: ...


def sample_revision_sets_equal(
    left: Sequence[str],
    right: Sequence[str],
) -> bool:
    """True when both sides are non-empty, duplicate-free, and equal ignoring order."""

    left_ids = tuple(left)
    right_ids = tuple(right)
    if not left_ids or not right_ids:
        return False
    if len(set(left_ids)) != len(left_ids) or len(set(right_ids)) != len(right_ids):
        return False
    return set(left_ids) == set(right_ids)


def experiment_sample_pin_identity(
    pin: ExperimentSamplePin,
) -> tuple[int, str, str, str | None, str, str | None, str | None]:
    """Ordinal, revision, digest, replicate group, and role stored on one pin."""

    term = pin.replicate_role.term
    return (
        pin.ordinal,
        pin.sample_revision_id,
        pin.payload_digest,
        pin.replicate_group,
        pin.replicate_role.status,
        None if term is None else term.term_id,
        None if term is None else term.term_revision_id,
    )


def association_pins_match(
    payload_pins: Sequence[
        tuple[int, str, str, str | None, str, str | None, str | None]
    ],
    mirror_pins: Sequence[
        tuple[int, str, str, str | None, str, str | None, str | None]
    ],
) -> bool:
    """True when every mirror row is one canonical sample pin.

    Each pin is ``(ordinal, sample_revision_id, payload_digest, replicate_group,
    role_status, role_term_id, role_term_revision_id)``. Order is ignored. A
    different ordinal, digest, replicate group, or role does not match, nor do
    a missing pin, an extra pin, or a duplicate ordinal or revision id.
    """

    payload = tuple(payload_pins)
    mirror = tuple(mirror_pins)
    if not payload or len(payload) != len(mirror):
        return False
    if len(set(payload)) != len(payload) or len(set(mirror)) != len(mirror):
        return False
    payload_ordinals = {item[0] for item in payload}
    mirror_ordinals = {item[0] for item in mirror}
    payload_ids = {item[1] for item in payload}
    mirror_ids = {item[1] for item in mirror}
    if (
        len(payload_ordinals) != len(payload)
        or len(mirror_ordinals) != len(mirror)
        or len(payload_ids) != len(payload)
        or len(mirror_ids) != len(mirror)
    ):
        return False
    return set(payload) == set(mirror)


def freeze_experiment_snapshot_binding(
    source: ExperimentRevisionFreezeSource | None,
    *,
    project_id: str,
    sample_revision_ids: Sequence[str],
) -> ExperimentSnapshotBinding:
    """Freeze one experiment revision or fail without a partial binding.

    Missing revisions, inactive projects, and cross-project ids are not found.
    An archived experiment is rejected before the sample-set comparison. The
    sample set is the canonical payload. The association mirror must match each
    pin's ordinal, revision id, digest, replicate group, and role; a drifted
    mirror is rejected before the request is compared. Equality ignores order
    and fails closed on an empty side, a duplicate, a missing pin, an extra pin,
    or the same sample at a different revision.
    """

    if (
        source is None
        or source.project_id != project_id
        or not source.project_active
        or _EXPERIMENT_REVISION_ID.fullmatch(source.experiment_revision_id) is None
        or _DIGEST.fullmatch(source.payload_digest) is None
    ):
        raise catalog_not_found()
    if source.experiment_archived:
        raise catalog_experiment_archived()
    if not source.association_matches_payload or not sample_revision_sets_equal(
        source.sample_revision_ids,
        sample_revision_ids,
    ):
        raise catalog_sample_set_mismatch()
    return ExperimentSnapshotBinding(
        experiment_revision_id=source.experiment_revision_id,
        payload_digest=source.payload_digest,
    )


def catalog_term_type_mismatch() -> CatalogError:
    return CatalogError(
        409,
        "catalog_term_type_mismatch",
        "The catalog term type does not match this field.",
    )


def administrator_required() -> CatalogError:
    return CatalogError(
        403,
        "ADMINISTRATOR_REQUIRED",
        "The operation requires an administrator.",
    )


def validate_idempotency_key(value: str | None) -> str:
    if not isinstance(value, str) or _IDEMPOTENCY_KEY.fullmatch(value) is None:
        raise CatalogError(
            400,
            "API_REQUEST_INVALID",
            "Idempotency-Key is required.",
        )
    return value


def parse_page_limit(value: int | None) -> int:
    if value is None:
        return DEFAULT_PAGE_SIZE
    if isinstance(value, bool) or not isinstance(value, int):
        raise catalog_filter_invalid()
    if value < 1 or value > MAX_PAGE_SIZE:
        raise catalog_filter_invalid()
    return value


def validate_display_name(value: object) -> str:
    if not isinstance(value, str):
        raise catalog_payload_invalid()
    normalized = value.strip()
    if (
        not normalized
        or normalized != value
        or len(normalized) > 255
        or any(ord(character) < 32 for character in normalized)
    ):
        raise catalog_payload_invalid()
    return normalized


def validate_stable_key(value: object) -> str:
    if not isinstance(value, str) or _STABLE_KEY.fullmatch(value) is None:
        raise catalog_payload_invalid()
    return value


def validate_optional_external(value: object) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise catalog_payload_invalid()
    if (
        not value
        or value != value.strip()
        or len(value) > 255
        or any(ord(character) < 32 for character in value)
    ):
        raise catalog_payload_invalid()
    return value


def validate_aliases(value: object) -> tuple[str, ...]:
    if value is None:
        return ()
    if not isinstance(value, list) or len(value) > 64:
        raise catalog_payload_invalid()
    aliases: list[str] = []
    for item in value:
        if not isinstance(item, str):
            raise catalog_payload_invalid()
        if (
            not item
            or item != item.strip()
            or len(item) > 255
            or any(ord(character) < 32 for character in item)
        ):
            raise catalog_payload_invalid()
        aliases.append(item)
    if len(set(aliases)) != len(aliases):
        raise catalog_payload_invalid()
    return tuple(aliases)


def canonical_alias_json(aliases: Sequence[str]) -> str:
    return _canonical_json(list(aliases))


def parse_page_filters(
    *,
    kind: str,
    project_id: str | None,
    revision_view: str | None,
    include_archived: bool,
    term_type: str | None,
    created_from: datetime | None,
    created_before: datetime | None,
    catalog_state: str | None,
    field_terms: Mapping[str, Sequence[str] | None],
    field_statuses: Mapping[str, str | None],
    owner_id: str | None = None,
) -> "CatalogFilters":
    if kind not in {
        "experiments",
        "samples",
        "terms",
        "sample_revisions",
        "experiment_revisions",
        "sample_runs",
        "experiment_runs",
    }:
        raise catalog_filter_invalid()
    if project_id is not None and _PROJECT_ID.fullmatch(project_id) is None:
        raise catalog_not_found()
    view = revision_view or "linked"
    if view not in REVISION_VIEWS:
        raise catalog_filter_invalid()
    if term_type is not None and term_type not in TERM_TYPES:
        raise catalog_filter_invalid()
    if catalog_state is not None and catalog_state not in {"absent", "present"}:
        raise catalog_filter_invalid()
    created_from = _optional_utc(created_from)
    created_before = _optional_utc(created_before)
    if (
        created_from is not None
        and created_before is not None
        and created_from >= created_before
    ):
        raise catalog_filter_invalid()
    filters: list[FieldFilter] = []
    names = tuple(dict.fromkeys((*field_terms.keys(), *field_statuses.keys())))
    for name in names:
        status = field_statuses.get(name)
        raw_terms = field_terms.get(name) or ()
        term_ids = tuple(dict.fromkeys(raw_terms))
        if status is not None and status not in FIELD_STATUSES:
            raise catalog_filter_invalid()
        for term_id in term_ids:
            if _TERM_ID.fullmatch(term_id) is None:
                raise catalog_filter_invalid()
        if status is not None and status != "known" and term_ids:
            raise catalog_filter_invalid()
        if status is None and not term_ids:
            continue
        filters.append(FieldFilter(name=name, status=status, term_ids=term_ids))
    return CatalogFilters(
        kind=kind,
        project_id=project_id,
        revision_view=view,
        include_archived=include_archived,
        term_type=term_type,
        created_from=created_from,
        created_before=created_before,
        catalog_state=catalog_state,
        owner_id=owner_id,
        fields=tuple(filters),
    )


@dataclass(frozen=True)
class CatalogFilters:
    kind: str
    project_id: str | None
    revision_view: str
    include_archived: bool
    term_type: str | None
    created_from: datetime | None
    created_before: datetime | None
    catalog_state: str | None
    owner_id: str | None
    fields: tuple[FieldFilter, ...]

    def field(self, name: str) -> FieldFilter | None:
        for item in self.fields:
            if item.name == name:
                return item
        return None

    @property
    def sha256(self) -> str:
        payload = {
            "catalog_state": self.catalog_state,
            "created_before": _optional_stamp(self.created_before),
            "created_from": _optional_stamp(self.created_from),
            "fields": [
                {
                    "name": item.name,
                    "status": item.status,
                    "term_ids": list(item.term_ids),
                }
                for item in self.fields
            ],
            "include_archived": self.include_archived,
            "kind": self.kind,
            "owner_id": self.owner_id,
            "project_id": self.project_id,
            "revision_view": self.revision_view,
            "term_type": self.term_type,
        }
        return sha256(_canonical_json(payload).encode("utf-8")).hexdigest()


def encode_catalog_cursor(
    *,
    kind: str,
    filters: CatalogFilters,
    position: CatalogCursorPosition,
) -> str:
    if filters.kind != kind:
        raise catalog_cursor_invalid()
    payload = {
        "filter_sha256": filters.sha256,
        "kind": kind,
        "position": {
            "created_at": _required_stamp(position.created_at),
            "object_id": position.object_id,
        },
        "v": 1,
    }
    raw = _canonical_json(payload).encode("utf-8")
    body = base64.urlsafe_b64encode(raw).decode("ascii").rstrip("=")
    encoded = f"{CATALOG_CURSOR_PREFIX}{body}"
    if len(encoded) > CATALOG_CURSOR_MAX_LENGTH or _CURSOR_BODY.fullmatch(body) is None:
        raise catalog_cursor_invalid()
    return encoded


def decode_catalog_cursor(
    value: str | None,
    *,
    kind: str,
    filters: CatalogFilters,
) -> CatalogCursorPosition | None:
    if value is None:
        return None
    if (
        not isinstance(value, str)
        or not value.startswith(CATALOG_CURSOR_PREFIX)
        or len(value) > CATALOG_CURSOR_MAX_LENGTH
    ):
        raise catalog_cursor_invalid()
    body = value[len(CATALOG_CURSOR_PREFIX) :]
    if _CURSOR_BODY.fullmatch(body) is None:
        raise catalog_cursor_invalid()
    padding = "=" * (-len(body) % 4)
    try:
        raw = base64.b64decode(body + padding, altchars=b"-_", validate=True)
        payload = json.loads(raw.decode("utf-8"))
    except (binascii.Error, UnicodeDecodeError, json.JSONDecodeError):
        raise catalog_cursor_invalid() from None
    if not isinstance(payload, dict) or set(payload) != _CURSOR_KEYS:
        raise catalog_cursor_invalid()
    if payload.get("v") != 1 or payload.get("kind") != kind:
        raise catalog_cursor_invalid()
    if payload.get("filter_sha256") != filters.sha256:
        raise catalog_cursor_invalid()
    position = payload.get("position")
    if not isinstance(position, dict) or set(position) != {"created_at", "object_id"}:
        raise catalog_cursor_invalid()
    created_at = _parse_stamp(position.get("created_at"))
    object_id = position.get("object_id")
    if created_at is None or not isinstance(object_id, str) or not object_id:
        raise catalog_cursor_invalid()
    decoded = CatalogCursorPosition(created_at=created_at, object_id=object_id)
    rebuilt = encode_catalog_cursor(kind=kind, filters=filters, position=decoded)
    if rebuilt != value:
        raise catalog_cursor_invalid()
    return decoded


class TermLookup:
    """Resolve a term revision. Missing revisions return None."""

    def __call__(self, term_revision_id: str) -> TermPin | None:
        raise NotImplementedError


def interpret_sample_catalog(
    value: object,
    lookup: TermLookup,
    *,
    for_write: bool,
) -> SampleCatalogBlock | None:
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise catalog_payload_invalid()
    if set(value) != {"schema_version", *SAMPLE_FIELDS}:
        raise catalog_payload_invalid()
    if value.get("schema_version") != CATALOG_SCHEMA_VERSION or isinstance(
        value.get("schema_version"), bool
    ):
        raise catalog_payload_invalid()
    return SampleCatalogBlock(
        species=_interpret_single(
            value["species"],
            expected_type="species",
            lookup=lookup,
            for_write=for_write,
        ),
        tissue=_interpret_single(
            value["tissue"],
            expected_type="tissue",
            lookup=lookup,
            for_write=for_write,
        ),
        cell_type=_interpret_multi(
            value["cell_type"],
            expected_type="cell_type",
            lookup=lookup,
            for_write=for_write,
        ),
        cell_line=_interpret_multi(
            value["cell_line"],
            expected_type="cell_line",
            lookup=lookup,
            for_write=for_write,
        ),
        treatment=_interpret_multi(
            value["treatment"],
            expected_type="treatment",
            lookup=lookup,
            for_write=for_write,
        ),
        tag=_interpret_multi(
            value["tag"],
            expected_type="tag",
            lookup=lookup,
            for_write=for_write,
        ),
    )


def sample_catalog_object(block: SampleCatalogBlock) -> dict[str, Any]:
    return {
        "schema_version": CATALOG_SCHEMA_VERSION,
        "species": _single_object(block.species),
        "tissue": _single_object(block.tissue),
        "cell_type": _multi_object(block.cell_type),
        "cell_line": _multi_object(block.cell_line),
        "treatment": _multi_object(block.treatment),
        "tag": _multi_object(block.tag),
    }


def canonical_catalog_json(block: SampleCatalogBlock) -> str:
    return _canonical_json(sample_catalog_object(block))


def interpret_experiment_fields(
    assay: object,
    target: object,
    batch: object,
    samples: object,
    lookup: TermLookup,
    *,
    for_write: bool,
) -> tuple[SingleSlot, SingleSlot, SingleSlot, tuple[tuple[object, SingleSlot], ...]]:
    """Validate experiment fields. Sample identity is resolved by the caller."""

    parsed_samples = _sample_requests(samples)
    roles: list[tuple[object, SingleSlot]] = []
    for request, role in parsed_samples:
        roles.append(
            (
                request,
                _interpret_single(
                    role,
                    expected_type="replicate_role",
                    lookup=lookup,
                    for_write=for_write,
                ),
            )
        )
    return (
        _interpret_single(
            assay,
            expected_type="assay",
            lookup=lookup,
            for_write=for_write,
        ),
        _interpret_single(
            target,
            expected_type="target",
            lookup=lookup,
            for_write=for_write,
        ),
        _interpret_single(
            batch,
            expected_type="batch",
            lookup=lookup,
            for_write=for_write,
        ),
        tuple(roles),
    )


def build_experiment_revision_payload_digest(canonical_payload: str) -> str:
    payload = interpret_stored_experiment_payload(canonical_payload)
    if payload.canonical_json() != canonical_payload:
        raise CatalogError(
            422,
            "catalog_payload_invalid",
            "The catalog payload is invalid.",
        )
    return _framed_sha256(
        EXPERIMENT_REVISION_PAYLOAD_DIGEST_SCHEME,
        canonical_payload,
    )


def interpret_stored_experiment_payload(canonical_payload: str) -> ExperimentPayload:
    try:
        decoded = json.loads(canonical_payload)
    except (json.JSONDecodeError, TypeError):
        raise catalog_payload_invalid() from None
    if not isinstance(decoded, dict):
        raise catalog_payload_invalid()
    if set(decoded) != {
        "assay",
        "batch",
        "display_name",
        "samples",
        "schema_version",
        "target",
    }:
        raise catalog_payload_invalid()
    if decoded["schema_version"] != CATALOG_SCHEMA_VERSION:
        raise catalog_payload_invalid()
    display_name = validate_display_name(decoded["display_name"])
    samples_raw = decoded["samples"]
    if not isinstance(samples_raw, list) or not samples_raw:
        raise catalog_payload_invalid()
    pins: list[ExperimentSamplePin] = []
    seen: set[str] = set()
    for ordinal, item in enumerate(samples_raw):
        if not isinstance(item, Mapping):
            raise catalog_payload_invalid()
        keys = set(item)
        if keys not in (
            {"ordinal", "payload_digest", "replicate_role", "sample_revision_id"},
            {
                "ordinal",
                "payload_digest",
                "replicate_group",
                "replicate_role",
                "sample_revision_id",
            },
        ):
            raise catalog_payload_invalid()
        if item["ordinal"] != ordinal:
            raise catalog_payload_invalid()
        sample_revision_id = item["sample_revision_id"]
        digest = item["payload_digest"]
        if (
            not isinstance(sample_revision_id, str)
            or _SAMPLE_REVISION_ID.fullmatch(sample_revision_id) is None
            or not isinstance(digest, str)
            or _DIGEST.fullmatch(digest) is None
        ):
            raise catalog_payload_invalid()
        if sample_revision_id in seen:
            raise catalog_payload_invalid()
        seen.add(sample_revision_id)
        group = item.get("replicate_group")
        if group == "":
            group = None
        elif group is not None and (
            not isinstance(group, str) or _STABLE_KEY.fullmatch(group) is None
        ):
            raise catalog_payload_invalid()
        role = _interpret_single(
            item["replicate_role"],
            expected_type="replicate_role",
            lookup=_RejectingLookup(),
            for_write=False,
        )
        pins.append(
            ExperimentSamplePin(
                sample_id="",
                sample_revision_id=sample_revision_id,
                payload_digest=digest,
                ordinal=ordinal,
                replicate_group=group if isinstance(group, str) else None,
                replicate_role=role,
            )
        )
    assay = _interpret_single(
        decoded["assay"],
        expected_type="assay",
        lookup=_RejectingLookup(),
        for_write=False,
    )
    target = _interpret_single(
        decoded["target"],
        expected_type="target",
        lookup=_RejectingLookup(),
        for_write=False,
    )
    batch = _interpret_single(
        decoded["batch"],
        expected_type="batch",
        lookup=_RejectingLookup(),
        for_write=False,
    )
    return ExperimentPayload(
        display_name=display_name,
        assay=assay,
        target=target,
        batch=batch,
        samples=tuple(pins),
    )


def render_single(
    slot: SingleSlot,
    names: Mapping[str, str],
    *,
    include_candidates: bool,
) -> dict[str, Any]:
    body: dict[str, Any] = {"status": slot.status}
    if slot.status == "known":
        if slot.term is None:
            raise catalog_payload_invalid()
        body["term_id"] = slot.term.term_id
        body["term_revision_id"] = slot.term.term_revision_id
        body["display_name"] = names[slot.term.term_revision_id]
    if include_candidates and slot.status == "ambiguous":
        body["candidates"] = [
            _candidate_object(item, names) for item in slot.candidates
        ]
    return body


def render_multi(
    slot: MultiSlot,
    names: Mapping[str, str],
    *,
    include_candidates: bool,
) -> dict[str, Any]:
    values: list[dict[str, str]] = []
    if slot.status == "known":
        values = [_candidate_object(item, names) for item in slot.values]
    body: dict[str, Any] = {"status": slot.status, "values": values}
    if include_candidates and slot.status == "ambiguous":
        body["candidates"] = [
            _candidate_object(item, names) for item in slot.candidates
        ]
    return body


def render_sample_material(
    block: SampleCatalogBlock | None,
    names: Mapping[str, str],
    *,
    include_candidates: bool,
) -> dict[str, Any]:
    if block is None:
        return {"catalog_state": "absent"}
    rendered = {
        name: (
            render_single(slot, names, include_candidates=include_candidates)
            if isinstance(slot, SingleSlot)
            else render_multi(slot, names, include_candidates=include_candidates)
        )
        for name, slot in block.slots().items()
    }
    rendered["catalog_state"] = "present"
    return rendered


def collect_term_revision_ids(block: SampleCatalogBlock | None) -> set[str]:
    found: set[str] = set()
    if block is None:
        return found
    for slot in block.slots().values():
        found.update(_slot_revision_ids(slot))
    return found


def collect_experiment_term_revision_ids(payload: ExperimentPayload) -> set[str]:
    found: set[str] = set()
    for slot in (payload.assay, payload.target, payload.batch):
        found.update(_slot_revision_ids(slot))
    for pin in payload.samples:
        found.update(_slot_revision_ids(pin.replicate_role))
    return found


class _RejectingLookup:
    """Parse stored pins without consulting the term dictionary."""

    def __call__(self, term_revision_id: str) -> TermPin | None:
        if _TERM_REVISION_ID.fullmatch(term_revision_id) is None:
            return None
        return TermPin(
            term_id="term_" + "0" * 32,
            term_revision_id=term_revision_id,
            term_type="replicate_role",
            display_name="stored",
            archived=False,
        )


def _sample_requests(
    samples: object,
) -> tuple[tuple[Mapping[str, Any], object], ...]:
    if not isinstance(samples, list) or not samples:
        raise catalog_payload_invalid()
    parsed: list[tuple[Mapping[str, Any], object]] = []
    seen: set[str] = set()
    for item in samples:
        if not isinstance(item, Mapping):
            raise catalog_payload_invalid()
        keys = set(item)
        allowed = {"sample_revision_id", "replicate_role"}
        if "replicate_group" in keys:
            allowed.add("replicate_group")
        if keys != allowed:
            raise catalog_payload_invalid()
        sample_revision_id = item.get("sample_revision_id")
        if (
            not isinstance(sample_revision_id, str)
            or _SAMPLE_REVISION_ID.fullmatch(sample_revision_id) is None
            or sample_revision_id in seen
        ):
            raise catalog_payload_invalid()
        seen.add(sample_revision_id)
        if "replicate_group" in item:
            group = item.get("replicate_group")
            if group == "":
                group = None
            elif not isinstance(group, str) or _STABLE_KEY.fullmatch(group) is None:
                raise catalog_payload_invalid()
        parsed.append((item, item.get("replicate_role")))
    return tuple(parsed)


def _interpret_single(
    value: object,
    *,
    expected_type: str,
    lookup: TermLookup,
    for_write: bool,
) -> SingleSlot:
    if not isinstance(value, Mapping):
        raise catalog_payload_invalid()
    status = value.get("status")
    if status not in FIELD_STATUSES:
        raise catalog_payload_invalid()
    keys = set(value)
    if status == "known":
        if keys != {"status", "term_id", "term_revision_id"}:
            raise catalog_payload_invalid()
        term = _require_pin(
            value.get("term_id"),
            value.get("term_revision_id"),
            expected_type=expected_type,
            lookup=lookup,
            for_write=for_write,
        )
        return SingleSlot(status="known", term=term)
    if status in NEGATIVE_STATUSES:
        if keys != {"status"}:
            raise catalog_payload_invalid()
        return SingleSlot(status=str(status))
    if keys not in ({"status", "candidates"}, {"status"}):
        raise catalog_payload_invalid()
    candidates = _candidates(
        value.get("candidates", []),
        expected_type=expected_type,
        lookup=lookup,
        for_write=for_write,
    )
    return SingleSlot(status="ambiguous", candidates=candidates)


def _interpret_multi(
    value: object,
    *,
    expected_type: str,
    lookup: TermLookup,
    for_write: bool,
) -> MultiSlot:
    if not isinstance(value, Mapping):
        raise catalog_payload_invalid()
    status = value.get("status")
    if status not in FIELD_STATUSES:
        raise catalog_payload_invalid()
    keys = set(value)
    if status == "known":
        if keys != {"status", "values"}:
            raise catalog_payload_invalid()
        values = _value_pins(
            value.get("values"),
            expected_type=expected_type,
            lookup=lookup,
            for_write=for_write,
        )
        if not values:
            raise catalog_payload_invalid()
        return MultiSlot(status="known", values=values)
    if "values" not in keys or value.get("values") not in ([], ()):
        raise catalog_payload_invalid()
    if status in NEGATIVE_STATUSES:
        if keys != {"status", "values"}:
            raise catalog_payload_invalid()
        return MultiSlot(status=str(status))
    if keys not in ({"status", "values", "candidates"}, {"status", "values"}):
        raise catalog_payload_invalid()
    candidates = _candidates(
        value.get("candidates", []),
        expected_type=expected_type,
        lookup=lookup,
        for_write=for_write,
    )
    return MultiSlot(status="ambiguous", candidates=candidates)


def _value_pins(
    value: object,
    *,
    expected_type: str,
    lookup: TermLookup,
    for_write: bool,
) -> tuple[PinnedTerm, ...]:
    if not isinstance(value, list):
        raise catalog_payload_invalid()
    pins: list[PinnedTerm] = []
    seen: set[str] = set()
    for item in value:
        if not isinstance(item, Mapping) or set(item) != {
            "term_id",
            "term_revision_id",
        }:
            raise catalog_payload_invalid()
        pin = _require_pin(
            item.get("term_id"),
            item.get("term_revision_id"),
            expected_type=expected_type,
            lookup=lookup,
            for_write=for_write,
        )
        if pin.term_id in seen:
            raise catalog_payload_invalid()
        seen.add(pin.term_id)
        pins.append(pin)
    return tuple(pins)


def _candidates(
    value: object,
    *,
    expected_type: str,
    lookup: TermLookup,
    for_write: bool,
) -> tuple[PinnedTerm, ...]:
    if value is None:
        return ()
    return _value_pins(
        value,
        expected_type=expected_type,
        lookup=lookup,
        for_write=for_write,
    )


def _require_pin(
    term_id: object,
    term_revision_id: object,
    *,
    expected_type: str,
    lookup: TermLookup,
    for_write: bool,
) -> PinnedTerm:
    if (
        not isinstance(term_id, str)
        or _TERM_ID.fullmatch(term_id) is None
        or not isinstance(term_revision_id, str)
        or _TERM_REVISION_ID.fullmatch(term_revision_id) is None
    ):
        raise catalog_payload_invalid()
    resolved = lookup(term_revision_id)
    if isinstance(lookup, _RejectingLookup):
        if resolved is None:
            raise catalog_payload_invalid()
        return PinnedTerm(term_id=term_id, term_revision_id=term_revision_id)
    if resolved is None or resolved.term_id != term_id:
        raise catalog_payload_invalid()
    if resolved.term_type != expected_type:
        raise catalog_term_type_mismatch()
    if for_write and resolved.archived:
        raise catalog_term_archived()
    return PinnedTerm(term_id=term_id, term_revision_id=term_revision_id)


def _single_object(slot: SingleSlot) -> dict[str, Any]:
    if slot.status == "known":
        if slot.term is None:
            raise catalog_payload_invalid()
        return {
            "status": "known",
            "term_id": slot.term.term_id,
            "term_revision_id": slot.term.term_revision_id,
        }
    if slot.status in NEGATIVE_STATUSES:
        return {"status": slot.status}
    body: dict[str, Any] = {"status": "ambiguous"}
    if slot.candidates:
        body["candidates"] = [
            {"term_id": item.term_id, "term_revision_id": item.term_revision_id}
            for item in slot.candidates
        ]
    return body


def _multi_object(slot: MultiSlot) -> dict[str, Any]:
    if slot.status == "known":
        return {
            "status": "known",
            "values": [
                {"term_id": item.term_id, "term_revision_id": item.term_revision_id}
                for item in slot.values
            ],
        }
    body: dict[str, Any] = {"status": slot.status, "values": []}
    if slot.status == "ambiguous" and slot.candidates:
        body["candidates"] = [
            {"term_id": item.term_id, "term_revision_id": item.term_revision_id}
            for item in slot.candidates
        ]
    return body


def _candidate_object(item: PinnedTerm, names: Mapping[str, str]) -> dict[str, str]:
    return {
        "term_id": item.term_id,
        "term_revision_id": item.term_revision_id,
        "display_name": names[item.term_revision_id],
    }


def _slot_revision_ids(slot: SingleSlot | MultiSlot) -> set[str]:
    found: set[str] = set()
    if isinstance(slot, SingleSlot) and slot.term is not None:
        found.add(slot.term.term_revision_id)
    if isinstance(slot, MultiSlot):
        found.update(item.term_revision_id for item in slot.values)
    found.update(item.term_revision_id for item in slot.candidates)
    return found


def _framed_sha256(*parts: str) -> str:
    digest = sha256()
    for part in parts:
        encoded = part.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return digest.hexdigest()


def _canonical_json(value: object) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def _required_stamp(value: datetime) -> str:
    return _require_utc(value).isoformat(timespec="microseconds").replace("+00:00", "Z")


def _optional_stamp(value: datetime | None) -> str | None:
    if value is None:
        return None
    return _required_stamp(value)


def _parse_stamp(value: object) -> datetime | None:
    if not isinstance(value, str) or not value.endswith("Z"):
        return None
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        return None
    if parsed.tzinfo is None:
        return None
    normalized = parsed.astimezone(timezone.utc)
    if _required_stamp(normalized) != value:
        return None
    return normalized


def _optional_utc(value: datetime | None) -> datetime | None:
    if value is None:
        return None
    return _require_utc(value)


def _require_utc(value: datetime) -> datetime:
    if (
        not isinstance(value, datetime)
        or value.tzinfo is None
        or value.utcoffset() is None
    ):
        raise catalog_filter_invalid()
    return value.astimezone(timezone.utc)


def slot_matches(slot: SingleSlot | MultiSlot, field: FieldFilter) -> bool:
    """Return whether one catalog slot satisfies one AND dimension."""

    if field.status is not None and slot.status != field.status:
        return False
    if not field.term_ids:
        return True
    if isinstance(slot, SingleSlot):
        return slot.term is not None and slot.term.term_id in field.term_ids
    return any(item.term_id in field.term_ids for item in slot.values)


def block_matches(block: SampleCatalogBlock | None, field: FieldFilter) -> bool:
    """Material filters never match an unclassified revision."""

    if block is None or field.name not in SAMPLE_FIELDS:
        return False
    return slot_matches(block.slots()[field.name], field)


def projection_rows_from_block(
    block: SampleCatalogBlock | None,
) -> tuple[dict[str, Any], tuple[dict[str, str], ...]]:
    """Return the scalar projection and known multi-value rows for one revision."""

    if block is None:
        return {"catalog_state": "absent"}, ()
    scalar: dict[str, Any] = {"catalog_state": "present"}
    values: list[dict[str, str]] = []
    for name in SAMPLE_SINGLE_FIELDS:
        slot = block.slots()[name]
        assert isinstance(slot, SingleSlot)
        scalar[f"{name}_status"] = slot.status
        scalar[f"{name}_term_id"] = None if slot.term is None else slot.term.term_id
        scalar[f"{name}_term_revision_id"] = (
            None if slot.term is None else slot.term.term_revision_id
        )
    for name in SAMPLE_MULTI_FIELDS:
        slot = block.slots()[name]
        assert isinstance(slot, MultiSlot)
        scalar[f"{name}_status"] = slot.status
        if slot.status == "known":
            for item in slot.values:
                values.append(
                    {
                        "field_name": name,
                        "term_id": item.term_id,
                        "term_revision_id": item.term_revision_id,
                    }
                )
    return scalar, tuple(values)


def block_from_projection(
    scalar: Mapping[str, Any] | None,
    values: Sequence[Mapping[str, str]],
) -> SampleCatalogBlock | None:
    if scalar is None or scalar.get("catalog_state") == "absent":
        return None
    if scalar.get("catalog_state") != "present":
        raise catalog_payload_invalid()

    def single(name: str) -> SingleSlot:
        status = scalar.get(f"{name}_status")
        if not isinstance(status, str):
            raise catalog_payload_invalid()
        term = None
        term_id = scalar.get(f"{name}_term_id")
        term_revision_id = scalar.get(f"{name}_term_revision_id")
        if status == "known":
            if not isinstance(term_id, str) or not isinstance(term_revision_id, str):
                raise catalog_payload_invalid()
            term = PinnedTerm(term_id=term_id, term_revision_id=term_revision_id)
        return SingleSlot(status=status, term=term)

    def multi(name: str) -> MultiSlot:
        status = scalar.get(f"{name}_status")
        if not isinstance(status, str):
            raise catalog_payload_invalid()
        pins = tuple(
            PinnedTerm(
                term_id=item["term_id"], term_revision_id=item["term_revision_id"]
            )
            for item in values
            if item["field_name"] == name
        )
        return MultiSlot(status=status, values=pins if status == "known" else ())

    return SampleCatalogBlock(
        species=single("species"),
        tissue=single("tissue"),
        cell_type=multi("cell_type"),
        cell_line=multi("cell_line"),
        treatment=multi("treatment"),
        tag=multi("tag"),
    )
