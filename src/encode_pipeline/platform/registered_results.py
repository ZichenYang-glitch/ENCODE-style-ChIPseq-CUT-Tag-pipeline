"""Independent identities for external and downstream registered results.

Registered results are not Runs. A revision records one observation of bytes
that already live in a bound storage pool, plus the source, method, and
reference the administrator actually supplied.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
import json
import re
from typing import Any


REVISION_DIGEST_SCHEME = "sha256-framed-registered-result-revision-v1"
SCHEMA_VERSION = 1

_PROJECT_ID = re.compile(r"^prj_[0-9a-f]{32}$")
_RESULT_ID = re.compile(r"^rreg_[0-9a-f]{32}$")
_REVISION_ID = re.compile(r"^rrev_[0-9a-f]{32}$")
_PROFILE_REVISION_ID = re.compile(r"^refpr_[0-9a-f]{32}$")
_STORAGE_POOL_ID = re.compile(r"^stgp_[0-9a-f]{32}$")
_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_STABLE_KEY = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]{0,254}$")
_FORMAT = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._+-]{0,31}$")
_ROLE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]{0,63}$")
_EVIDENCE = re.compile(r"^[^\x00-\x1f\x7f]{1,128}$")
_MAX_PARENTS = 32


class RegisteredResultPayloadError(ValueError):
    """The caller supplied metadata that cannot become a revision."""


@dataclass(frozen=True)
class ParentEdge:
    """One ordered parent pin stored with a downstream revision."""

    result_revision_id: str
    role: str
    revision_digest: str

    def __post_init__(self) -> None:
        _revision_id(self.result_revision_id)
        _role(self.role)
        _digest(self.revision_digest, "revision_digest")


@dataclass(frozen=True)
class ParentView:
    """Parent fields returned with a historical revision."""

    result_id: str
    result_revision_id: str
    revision_number: int
    revision_status: str
    current_result_revision_id: str
    display_name: str
    role: str
    revision_digest: str


@dataclass(frozen=True)
class RegisteredMetadata:
    """Normalized description stored in one revision payload."""

    display_name: str
    format: str
    origin: dict[str, Any]
    method: dict[str, Any]
    reference: dict[str, Any]

    def canonical_json(self) -> str:
        """Return the only byte sequence stored as this revision's payload."""

        return _canonical_json(
            {
                "schema_version": SCHEMA_VERSION,
                "display_name": self.display_name,
                "format": self.format,
                "origin": self.origin,
                "method": self.method,
                "reference": self.reference,
            }
        )

    @property
    def parent_refs(self) -> tuple[tuple[str, str], ...]:
        """Return ordered ``(result_revision_id, role)`` pins, if downstream."""

        if self.origin["kind"] != "downstream":
            return ()
        return tuple(
            (parent["result_revision_id"], parent["role"])
            for parent in self.origin["parents"]
        )


@dataclass(frozen=True)
class RevisionRecord:
    """One immutable revision plus the private coordinates needed to re-read it."""

    project_id: str
    project_archived_at: datetime | None
    result_id: str
    stable_key: str
    result_archived_at: datetime | None
    result_revision_id: str
    revision_number: int
    canonical_payload: str
    revision_digest: str
    created_at: datetime
    storage_pool_id: str
    relative_path: str
    size_bytes: int
    content_sha256: str
    head_result_revision_id: str
    parents: tuple[ParentView, ...]

    @property
    def revision_status(self) -> str:
        if self.result_revision_id == self.head_result_revision_id:
            return "current"
        return "superseded"


def validate_project_id(value: object) -> str:
    return _identifier(value, _PROJECT_ID, "project_id")


def validate_result_id(value: object) -> str:
    return _identifier(value, _RESULT_ID, "result_id")


def validate_result_revision_id(value: object) -> str:
    return _identifier(value, _REVISION_ID, "result_revision_id")


def validate_stable_key(value: object) -> str:
    if not isinstance(value, str) or _STABLE_KEY.fullmatch(value) is None:
        raise RegisteredResultPayloadError("stable_key is invalid")
    return value


def validate_idempotency_key(value: object) -> str:
    if (
        not isinstance(value, str)
        or not value
        or len(value) > 128
        or re.fullmatch(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$", value) is None
    ):
        raise RegisteredResultPayloadError("Idempotency-Key is invalid")
    return value


def parse_metadata(payload: Mapping[str, Any]) -> RegisteredMetadata:
    """Normalize one revision description. Unknown evidence stays unknown."""

    if not isinstance(payload, Mapping):
        raise RegisteredResultPayloadError("metadata must be an object")
    allowed = {
        "display_name",
        "format",
        "origin",
        "method",
        "reference",
        "schema_version",
    }
    if not set(payload) <= allowed or not {
        "display_name",
        "format",
        "origin",
        "method",
        "reference",
    } <= set(payload):
        raise RegisteredResultPayloadError("metadata is invalid")
    if "schema_version" in payload and payload["schema_version"] != SCHEMA_VERSION:
        raise RegisteredResultPayloadError("schema_version is unsupported")
    display_name = _display_name(payload.get("display_name"))
    result_format = _format_token(payload.get("format"))
    origin = _origin(payload.get("origin"))
    method = _method(payload.get("method"))
    reference = _reference(payload.get("reference"))
    return RegisteredMetadata(
        display_name=display_name,
        format=result_format,
        origin=origin,
        method=method,
        reference=reference,
    )


def metadata_from_canonical_json(canonical_payload: str) -> RegisteredMetadata:
    """Decode a stored payload and reject anything that is no longer canonical."""

    if not isinstance(canonical_payload, str) or not canonical_payload:
        raise RegisteredResultPayloadError("canonical payload is empty")
    try:
        decoded = json.loads(canonical_payload)
    except json.JSONDecodeError as exc:
        raise RegisteredResultPayloadError("canonical payload is not JSON") from exc
    if not isinstance(decoded, dict):
        raise RegisteredResultPayloadError("canonical payload must be an object")
    metadata = parse_metadata(decoded)
    if metadata.canonical_json() != canonical_payload:
        raise RegisteredResultPayloadError("canonical payload is not canonical")
    if decoded.get("schema_version") != SCHEMA_VERSION:
        raise RegisteredResultPayloadError("schema_version is unsupported")
    return metadata


def request_sha256(
    operation: str,
    project_id: str,
    result_id: str | None,
    body: Mapping[str, Any],
) -> str:
    """Hash the normalized write so key order does not change idempotency."""

    encoded = _canonical_json(
        {
            "body": body,
            "operation": operation,
            "project_id": validate_project_id(project_id),
            "result_id": None if result_id is None else validate_result_id(result_id),
        }
    ).encode("utf-8")
    return sha256(encoded).hexdigest()


def build_revision_digest(
    *,
    project_id: str,
    result_id: str,
    result_revision_id: str,
    revision_number: int,
    canonical_payload: str,
    storage_pool_id: str,
    relative_path: str,
    size_bytes: int,
    content_sha256: str,
    parents: tuple[ParentEdge, ...],
) -> str:
    """Frame every coordinate that makes this revision immutable."""

    if not isinstance(revision_number, int) or isinstance(revision_number, bool):
        raise RegisteredResultPayloadError("revision_number is invalid")
    if revision_number < 1:
        raise RegisteredResultPayloadError("revision_number is invalid")
    if not isinstance(size_bytes, int) or isinstance(size_bytes, bool):
        raise RegisteredResultPayloadError("size_bytes is invalid")
    if size_bytes < 0:
        raise RegisteredResultPayloadError("size_bytes is invalid")
    parent_document = _canonical_json(
        [
            {
                "ordinal": index,
                "result_revision_id": parent.result_revision_id,
                "revision_digest": parent.revision_digest,
                "role": parent.role,
            }
            for index, parent in enumerate(parents)
        ]
    )
    return _framed_sha256(
        REVISION_DIGEST_SCHEME,
        validate_project_id(project_id),
        validate_result_id(result_id),
        validate_result_revision_id(result_revision_id),
        str(revision_number),
        canonical_payload,
        _identifier(storage_pool_id, _STORAGE_POOL_ID, "storage_pool_id"),
        relative_path,
        str(size_bytes),
        _digest(content_sha256, "content_sha256"),
        parent_document,
    )


def public_revision(record: RevisionRecord) -> dict[str, Any]:
    """Return the bounded read DTO. Private coordinates stay out."""

    metadata = metadata_from_canonical_json(record.canonical_payload)
    if metadata.origin["kind"] == "downstream":
        origin: dict[str, Any] = {
            "kind": "downstream",
            "parents": [
                {
                    "result_id": parent.result_id,
                    "result_revision_id": parent.result_revision_id,
                    "revision_number": parent.revision_number,
                    "revision_status": parent.revision_status,
                    "current_result_revision_id": parent.current_result_revision_id,
                    "display_name": parent.display_name,
                    "role": parent.role,
                }
                for parent in record.parents
            ],
        }
    else:
        origin = metadata.origin
    created_at = record.created_at
    if created_at.tzinfo is None or created_at.utcoffset() is None:
        created_at = created_at.replace(tzinfo=timezone.utc)
    created_at = created_at.astimezone(timezone.utc)
    return {
        "project_id": record.project_id,
        "result_id": record.result_id,
        "result_revision_id": record.result_revision_id,
        "revision_number": record.revision_number,
        "revision_status": record.revision_status,
        "display_name": metadata.display_name,
        "format": metadata.format,
        "schema_version": SCHEMA_VERSION,
        "origin": origin,
        "method": metadata.method,
        "reference": metadata.reference,
        "bytes": {
            "size_bytes": record.size_bytes,
            "content_sha256": record.content_sha256,
        },
        "created_at": created_at.isoformat().replace("+00:00", "Z"),
    }


def _origin(value: object) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise RegisteredResultPayloadError("origin is invalid")
    fields = set(value)
    if fields not in ({"kind", "source"}, {"kind", "parents"}):
        raise RegisteredResultPayloadError("origin is invalid")
    kind = value.get("kind")
    if kind == "external":
        if set(value) != {"kind", "source"}:
            raise RegisteredResultPayloadError("origin is invalid")
        return {"kind": "external", "source": _source(value.get("source"))}
    if kind == "downstream":
        if set(value) != {"kind", "parents"}:
            raise RegisteredResultPayloadError("origin is invalid")
        return {"kind": "downstream", "parents": _parents(value.get("parents"))}
    raise RegisteredResultPayloadError("origin is invalid")


def _source(value: object) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise RegisteredResultPayloadError("source is invalid")
    status = value.get("status")
    if status == "unknown":
        if set(value) != {"status"}:
            raise RegisteredResultPayloadError("source is invalid")
        return {"status": "unknown"}
    if status == "known":
        if set(value) != {"status", "namespace", "accession", "version"}:
            raise RegisteredResultPayloadError("source is invalid")
        return {
            "status": "known",
            "namespace": _evidence(value.get("namespace"), "namespace"),
            "accession": _evidence(value.get("accession"), "accession"),
            "version": _evidence(value.get("version"), "version"),
        }
    raise RegisteredResultPayloadError("source is invalid")


def _parents(value: object) -> list[dict[str, str]]:
    if not isinstance(value, list) or not value or len(value) > _MAX_PARENTS:
        raise RegisteredResultPayloadError("parents are invalid")
    parents: list[dict[str, str]] = []
    seen: set[str] = set()
    for item in value:
        if not isinstance(item, Mapping) or set(item) != {
            "result_revision_id",
            "role",
        }:
            raise RegisteredResultPayloadError("parents are invalid")
        revision_id = _revision_id(item.get("result_revision_id"))
        if revision_id in seen:
            raise RegisteredResultPayloadError("parents are invalid")
        seen.add(revision_id)
        parents.append(
            {
                "result_revision_id": revision_id,
                "role": _role(item.get("role")),
            }
        )
    return parents


def _method(value: object) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise RegisteredResultPayloadError("method is invalid")
    status = value.get("status")
    if status == "unknown":
        if set(value) != {"status"}:
            raise RegisteredResultPayloadError("method is invalid")
        return {"status": "unknown"}
    if status == "known" and set(value) == {"status", "name", "version"}:
        return {
            "status": "known",
            "name": _evidence(value.get("name"), "name"),
            "version": _evidence(value.get("version"), "version"),
        }
    raise RegisteredResultPayloadError("method is invalid")


def _reference(value: object) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise RegisteredResultPayloadError("reference is invalid")
    status = value.get("status")
    if status == "unknown":
        if set(value) != {"status"}:
            raise RegisteredResultPayloadError("reference is invalid")
        return {"status": "unknown"}
    if status == "not_applicable":
        if set(value) != {"status"}:
            raise RegisteredResultPayloadError("reference is invalid")
        return {"status": "not_applicable"}
    if status != "known":
        raise RegisteredResultPayloadError("reference is invalid")
    kind = value.get("kind")
    if kind == "external" and set(value) == {
        "status",
        "kind",
        "namespace",
        "accession",
        "version",
    }:
        return {
            "status": "known",
            "kind": "external",
            "namespace": _evidence(value.get("namespace"), "namespace"),
            "accession": _evidence(value.get("accession"), "accession"),
            "version": _evidence(value.get("version"), "version"),
        }
    if kind == "reference_profile_revision" and set(value) == {
        "status",
        "kind",
        "reference_profile_revision_id",
        "public_identity_sha256",
    }:
        revision_id = value.get("reference_profile_revision_id")
        if (
            not isinstance(revision_id, str)
            or _PROFILE_REVISION_ID.fullmatch(revision_id) is None
        ):
            raise RegisteredResultPayloadError("reference is invalid")
        return {
            "status": "known",
            "kind": "reference_profile_revision",
            "reference_profile_revision_id": revision_id,
            "public_identity_sha256": _digest(
                value.get("public_identity_sha256"),
                "public_identity_sha256",
            ),
        }
    raise RegisteredResultPayloadError("reference is invalid")


def _display_name(value: object) -> str:
    if not isinstance(value, str):
        raise RegisteredResultPayloadError("display_name is invalid")
    normalized = value.strip()
    if (
        not normalized
        or len(normalized) > 255
        or any(ord(character) < 32 or ord(character) == 127 for character in normalized)
    ):
        raise RegisteredResultPayloadError("display_name is invalid")
    return normalized


def _format_token(value: object) -> str:
    if not isinstance(value, str) or _FORMAT.fullmatch(value) is None:
        raise RegisteredResultPayloadError("format is invalid")
    return value


def _role(value: object) -> str:
    if not isinstance(value, str) or _ROLE.fullmatch(value) is None:
        raise RegisteredResultPayloadError("role is invalid")
    return value


def _evidence(value: object, name: str) -> str:
    if (
        not isinstance(value, str)
        or value != value.strip()
        or _EVIDENCE.fullmatch(value) is None
        or "/" in value
        or "\\" in value
        or "://" in value
    ):
        raise RegisteredResultPayloadError(f"{name} is invalid")
    return value


def _revision_id(value: object) -> str:
    return _identifier(value, _REVISION_ID, "result_revision_id")


def _digest(value: object, name: str) -> str:
    if not isinstance(value, str) or _DIGEST.fullmatch(value) is None:
        raise RegisteredResultPayloadError(f"{name} is invalid")
    return value


def _identifier(value: object, pattern: re.Pattern[str], name: str) -> str:
    if not isinstance(value, str) or pattern.fullmatch(value) is None:
        raise RegisteredResultPayloadError(f"{name} is invalid")
    return value


def _canonical_json(value: object) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def _framed_sha256(*parts: str) -> str:
    digest = sha256()
    for part in parts:
        encoded = part.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return digest.hexdigest()
