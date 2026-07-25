"""Workflow-neutral Project, Sample, and immutable binding primitives."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from enum import Enum
from hashlib import sha256
import json
import math
import re
from typing import Any


LEGACY_PROJECT_ID = "prj_00000000000000000000000000000000"
LEGACY_PROJECT_DISPLAY_NAME = "Legacy Project"
SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME = "sha256-framed-sample-revision-payload-v1"
PROJECT_SAMPLE_BINDING_DIGEST_SCHEME = "sha256-framed-project-sample-binding-v1"

_PROJECT_ID = re.compile(r"^prj_[0-9a-f]{32}$")
_SAMPLE_ID = re.compile(r"^smp_[0-9a-f]{32}$")
_SAMPLE_REVISION_ID = re.compile(r"^smpr_[0-9a-f]{32}$")
_STABLE_KEY = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]{0,254}$")
_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_MAX_SAFE_JSON_INTEGER = 2**53 - 1


class ProjectKind(str, Enum):
    """Lifecycle category for an operator-owned or reserved Project."""

    USER = "user"
    SYSTEM = "system"


class BindingMode(str, Enum):
    """Compatibility mode for one immutable snapshot/run data binding."""

    LEGACY_V1 = "legacy_v1"
    BOUND_V1 = "bound_v1"


class BindingProvenance(str, Enum):
    """Whether exact SampleRevision provenance is present."""

    RESOLVED = "resolved"
    UNRESOLVED = "unresolved"


def _utc(value: datetime, name: str) -> datetime:
    if (
        not isinstance(value, datetime)
        or value.tzinfo is None
        or value.utcoffset() is None
    ):
        raise ValueError(f"{name} must be timezone-aware")
    return value.astimezone(timezone.utc)


def _non_empty_display_name(value: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError("display_name must be a non-empty string")
    normalized = value.strip()
    if len(normalized) > 255 or any(ord(character) < 32 for character in normalized):
        raise ValueError("display_name is invalid")
    return normalized


def _validate_id(value: str, pattern: re.Pattern[str], name: str) -> str:
    if not isinstance(value, str) or pattern.fullmatch(value) is None:
        raise ValueError(f"{name} must be an opaque {name.replace('_', ' ')}")
    return value


def _validate_digest(value: str, name: str) -> str:
    if not isinstance(value, str) or _DIGEST.fullmatch(value) is None:
        raise ValueError(f"{name} must be a lowercase SHA-256 digest")
    return value


def _normalize_json(value: object, *, depth: int = 0) -> Any:
    if depth > 64:
        raise ValueError("sample attributes must be JSON-safe")
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, int):
        if abs(value) > _MAX_SAFE_JSON_INTEGER:
            raise ValueError("sample attributes must use JSON-safe integers")
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("sample attributes must use JSON-safe numbers")
        return value
    if isinstance(value, (list, tuple)):
        return [_normalize_json(item, depth=depth + 1) for item in value]
    if isinstance(value, Mapping):
        normalized: dict[str, Any] = {}
        for key, item in value.items():
            if not isinstance(key, str):
                raise ValueError(
                    "sample attributes must use string JSON keys for JSON-safe data"
                )
            normalized[key] = _normalize_json(item, depth=depth + 1)
        return normalized
    raise ValueError("sample attributes must be JSON-safe")


def canonical_sample_revision_payload(
    *,
    display_name: str,
    attributes: Mapping[str, object],
) -> str:
    """Return canonical JSON for one complete immutable SampleRevision payload."""
    normalized_name = _non_empty_display_name(display_name)
    if not isinstance(attributes, Mapping):
        raise ValueError("sample attributes must be a JSON-safe object")
    payload = {
        "display_name": normalized_name,
        "attributes": _normalize_json(attributes),
    }
    return json.dumps(
        payload,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def _decode_sample_revision_payload(canonical_payload: str) -> dict[str, Any]:
    if not isinstance(canonical_payload, str) or not canonical_payload:
        raise ValueError("canonical_payload must be non-empty canonical JSON")
    try:
        decoded = json.loads(canonical_payload)
    except (json.JSONDecodeError, TypeError):
        raise ValueError("canonical_payload must contain canonical JSON") from None
    if not isinstance(decoded, dict) or set(decoded) != {
        "attributes",
        "display_name",
    }:
        raise ValueError("canonical_payload must contain display_name and attributes")
    attributes = decoded["attributes"]
    if not isinstance(attributes, dict):
        raise ValueError("canonical_payload attributes must be an object")
    try:
        rebuilt = canonical_sample_revision_payload(
            display_name=decoded["display_name"],
            attributes=attributes,
        )
    except (TypeError, ValueError):
        raise ValueError("canonical_payload is invalid") from None
    if rebuilt != canonical_payload:
        raise ValueError("canonical_payload is not canonical")
    return decoded


def _framed_sha256(*parts: str) -> str:
    digest = sha256()
    for part in parts:
        encoded = part.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return digest.hexdigest()


def build_sample_revision_payload_digest(canonical_payload: str) -> str:
    """Return the versioned length-framed digest for a canonical payload."""
    _decode_sample_revision_payload(canonical_payload)
    return _framed_sha256(
        SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
        canonical_payload,
    )


@dataclass(frozen=True)
class Project:
    """Stable lifecycle and query scope for data registry objects."""

    project_id: str
    display_name: str
    kind: ProjectKind | str
    created_at: datetime
    archived_at: datetime | None = None

    def __post_init__(self) -> None:
        project_id = _validate_id(self.project_id, _PROJECT_ID, "project_id")
        display_name = _non_empty_display_name(self.display_name)
        try:
            kind = (
                self.kind
                if isinstance(self.kind, ProjectKind)
                else ProjectKind(self.kind)
            )
        except (TypeError, ValueError):
            raise ValueError("kind must be user or system") from None
        created_at = _utc(self.created_at, "created_at")
        archived_at = (
            None if self.archived_at is None else _utc(self.archived_at, "archived_at")
        )
        if archived_at is not None and archived_at < created_at:
            raise ValueError("archived_at cannot precede created_at")

        if project_id == LEGACY_PROJECT_ID:
            if (
                display_name != LEGACY_PROJECT_DISPLAY_NAME
                or kind is not ProjectKind.SYSTEM
                or archived_at is not None
            ):
                raise ValueError("reserved Legacy Project is immutable")
        elif kind is ProjectKind.SYSTEM:
            raise ValueError("only the reserved Legacy Project may be system")

        object.__setattr__(self, "project_id", project_id)
        object.__setattr__(self, "display_name", display_name)
        object.__setattr__(self, "kind", kind)
        object.__setattr__(self, "created_at", created_at)
        object.__setattr__(self, "archived_at", archived_at)


def build_legacy_project(created_at: datetime) -> Project:
    """Return the deterministic reserved Legacy Project."""
    return Project(
        project_id=LEGACY_PROJECT_ID,
        display_name=LEGACY_PROJECT_DISPLAY_NAME,
        kind=ProjectKind.SYSTEM,
        created_at=created_at,
    )


@dataclass(frozen=True)
class Sample:
    """Stable logical sample identity, separate from adapter sample rows."""

    sample_id: str
    project_id: str
    stable_key: str
    created_at: datetime

    def __post_init__(self) -> None:
        sample_id = _validate_id(self.sample_id, _SAMPLE_ID, "sample_id")
        project_id = _validate_id(self.project_id, _PROJECT_ID, "project_id")
        if project_id == LEGACY_PROJECT_ID:
            raise ValueError("Legacy Project cannot own trusted Samples")
        if (
            not isinstance(self.stable_key, str)
            or _STABLE_KEY.fullmatch(self.stable_key) is None
        ):
            raise ValueError("stable_key must be a safe project-scoped key")
        object.__setattr__(self, "sample_id", sample_id)
        object.__setattr__(self, "project_id", project_id)
        object.__setattr__(self, "created_at", _utc(self.created_at, "created_at"))


@dataclass(frozen=True)
class SampleRevision:
    """One immutable append-only revision of workflow-neutral sample metadata."""

    sample_revision_id: str
    project_id: str
    sample_id: str
    revision_number: int
    canonical_payload: str
    payload_digest_scheme: str
    payload_digest: str
    created_at: datetime

    def __post_init__(self) -> None:
        _validate_id(
            self.sample_revision_id,
            _SAMPLE_REVISION_ID,
            "sample_revision_id",
        )
        project_id = _validate_id(self.project_id, _PROJECT_ID, "project_id")
        if project_id == LEGACY_PROJECT_ID:
            raise ValueError("Legacy Project cannot own trusted SampleRevisions")
        _validate_id(self.sample_id, _SAMPLE_ID, "sample_id")
        if (
            not isinstance(self.revision_number, int)
            or isinstance(self.revision_number, bool)
            or self.revision_number < 1
        ):
            raise ValueError("revision_number must be a positive integer")
        if self.payload_digest_scheme != SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME:
            raise ValueError("sample revision payload digest scheme is unsupported")
        _decode_sample_revision_payload(self.canonical_payload)
        payload_digest = _validate_digest(self.payload_digest, "payload_digest")
        if payload_digest != build_sample_revision_payload_digest(
            self.canonical_payload
        ):
            raise ValueError("payload digest does not match canonical payload")
        object.__setattr__(self, "created_at", _utc(self.created_at, "created_at"))

    def to_payload(self) -> dict[str, Any]:
        """Return a fresh decoded copy of the immutable payload."""
        return _decode_sample_revision_payload(self.canonical_payload)

    @property
    def display_name(self) -> str:
        """Return the revision-owned display name."""
        value = self.to_payload()["display_name"]
        assert isinstance(value, str)
        return value

    @property
    def attributes(self) -> dict[str, Any]:
        """Return a fresh copy of the revision-owned attributes."""
        value = self.to_payload()["attributes"]
        assert isinstance(value, dict)
        return value


@dataclass(frozen=True)
class ProjectSampleSelection:
    """Caller selection of exact ordered revisions for a non-Legacy Project."""

    project_id: str
    sample_revision_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        project_id = _validate_id(self.project_id, _PROJECT_ID, "project_id")
        if project_id == LEGACY_PROJECT_ID:
            raise ValueError("Legacy Project does not accept sample selections")
        revisions = tuple(self.sample_revision_ids)
        if not revisions:
            raise ValueError("selection must contain at least one sample revision")
        for revision_id in revisions:
            _validate_id(revision_id, _SAMPLE_REVISION_ID, "sample_revision_id")
        if len(set(revisions)) != len(revisions):
            raise ValueError("selection contains duplicate sample revisions")
        object.__setattr__(self, "project_id", project_id)
        object.__setattr__(self, "sample_revision_ids", revisions)


@dataclass(frozen=True)
class SampleRevisionBindingRef:
    """Digest-pinned reference to one exact SampleRevision."""

    sample_revision_id: str
    payload_digest: str

    def __post_init__(self) -> None:
        _validate_id(
            self.sample_revision_id,
            _SAMPLE_REVISION_ID,
            "sample_revision_id",
        )
        _validate_digest(self.payload_digest, "payload_digest")


def canonical_project_sample_binding(
    *,
    project_id: str,
    binding_mode: BindingMode,
    provenance: BindingProvenance,
    workflow_inputs_digest: str,
    sample_revisions: tuple[SampleRevisionBindingRef, ...],
) -> str:
    """Return canonical JSON for the complete Stage 2 binding evidence."""
    payload = {
        "binding_mode": binding_mode.value,
        "project_id": project_id,
        "provenance": provenance.value,
        "sample_revisions": [
            {
                "payload_digest": ref.payload_digest,
                "sample_revision_id": ref.sample_revision_id,
            }
            for ref in sample_revisions
        ],
        "workflow_inputs_digest": workflow_inputs_digest,
    }
    return json.dumps(
        payload,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    )


@dataclass(frozen=True)
class ProjectSampleBinding:
    """Immutable digest-pinned Project/Sample binding copied to a run."""

    project_id: str
    binding_mode: BindingMode | str
    provenance: BindingProvenance | str
    workflow_inputs_digest: str
    sample_revisions: tuple[SampleRevisionBindingRef, ...]
    digest_scheme: str
    digest: str

    def __post_init__(self) -> None:
        project_id = _validate_id(self.project_id, _PROJECT_ID, "project_id")
        try:
            binding_mode = (
                self.binding_mode
                if isinstance(self.binding_mode, BindingMode)
                else BindingMode(self.binding_mode)
            )
            provenance = (
                self.provenance
                if isinstance(self.provenance, BindingProvenance)
                else BindingProvenance(self.provenance)
            )
        except (TypeError, ValueError):
            raise ValueError("binding mode or provenance is unsupported") from None
        workflow_inputs_digest = _validate_digest(
            self.workflow_inputs_digest,
            "workflow_inputs_digest",
        )
        sample_revisions = tuple(self.sample_revisions)
        if any(
            not isinstance(ref, SampleRevisionBindingRef) for ref in sample_revisions
        ):
            raise ValueError("sample_revisions must contain binding references")
        sample_ids = tuple(ref.sample_revision_id for ref in sample_revisions)
        if len(set(sample_ids)) != len(sample_ids):
            raise ValueError("binding contains duplicate sample revisions")

        if binding_mode is BindingMode.LEGACY_V1:
            if (
                project_id != LEGACY_PROJECT_ID
                or provenance is not BindingProvenance.UNRESOLVED
                or sample_revisions
            ):
                raise ValueError("legacy binding evidence is inconsistent")
        elif (
            project_id == LEGACY_PROJECT_ID
            or provenance is not BindingProvenance.RESOLVED
            or not sample_revisions
        ):
            raise ValueError("bound binding evidence is inconsistent")

        if self.digest_scheme != PROJECT_SAMPLE_BINDING_DIGEST_SCHEME:
            raise ValueError("Project/Sample binding digest scheme is unsupported")
        digest = _validate_digest(self.digest, "digest")
        canonical_payload = canonical_project_sample_binding(
            project_id=project_id,
            binding_mode=binding_mode,
            provenance=provenance,
            workflow_inputs_digest=workflow_inputs_digest,
            sample_revisions=sample_revisions,
        )
        expected_digest = _framed_sha256(
            PROJECT_SAMPLE_BINDING_DIGEST_SCHEME,
            canonical_payload,
        )
        if digest != expected_digest:
            raise ValueError("binding digest does not match canonical evidence")

        object.__setattr__(self, "project_id", project_id)
        object.__setattr__(self, "binding_mode", binding_mode)
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(
            self,
            "workflow_inputs_digest",
            workflow_inputs_digest,
        )
        object.__setattr__(self, "sample_revisions", sample_revisions)

    @property
    def sample_revision_ids(self) -> tuple[str, ...]:
        """Return the exact ordered revision identities."""
        return tuple(ref.sample_revision_id for ref in self.sample_revisions)


def build_project_sample_binding(
    *,
    project_id: str,
    binding_mode: BindingMode,
    provenance: BindingProvenance,
    workflow_inputs_digest: str,
    sample_revisions: tuple[SampleRevisionBindingRef, ...],
) -> ProjectSampleBinding:
    """Build and verify a canonical Project/Sample binding."""
    canonical_payload = canonical_project_sample_binding(
        project_id=project_id,
        binding_mode=binding_mode,
        provenance=provenance,
        workflow_inputs_digest=workflow_inputs_digest,
        sample_revisions=tuple(sample_revisions),
    )
    return ProjectSampleBinding(
        project_id=project_id,
        binding_mode=binding_mode,
        provenance=provenance,
        workflow_inputs_digest=workflow_inputs_digest,
        sample_revisions=tuple(sample_revisions),
        digest_scheme=PROJECT_SAMPLE_BINDING_DIGEST_SCHEME,
        digest=_framed_sha256(
            PROJECT_SAMPLE_BINDING_DIGEST_SCHEME,
            canonical_payload,
        ),
    )


def build_legacy_project_sample_binding(
    workflow_inputs_digest: str,
) -> ProjectSampleBinding:
    """Build conservative unresolved binding evidence for a Legacy v1 payload."""
    return build_project_sample_binding(
        project_id=LEGACY_PROJECT_ID,
        binding_mode=BindingMode.LEGACY_V1,
        provenance=BindingProvenance.UNRESOLVED,
        workflow_inputs_digest=workflow_inputs_digest,
        sample_revisions=(),
    )
