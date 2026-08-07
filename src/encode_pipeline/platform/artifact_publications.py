"""Workflow-neutral, disclosure-safe artifact-publication query contracts."""

from __future__ import annotations

import base64
import binascii
from dataclasses import dataclass
from datetime import datetime, timezone
from hashlib import sha256
import json
import re
from typing import Protocol

from encode_pipeline.platform.data_registry import (
    LEGACY_PROJECT_ID,
    BindingMode,
    BindingProvenance,
)
from encode_pipeline.platform.result_generations import (
    validate_artifact_generation,
    validate_artifact_revision,
)


ARTIFACT_PUBLICATION_CURSOR_PREFIX = "artifactpubcur_"
ARTIFACT_PUBLICATION_CURSOR_MAX_LENGTH = 1024

_CURSOR_BODY_PATTERN = re.compile(r"^[A-Za-z0-9_-]+$")
_DIGEST_PATTERN = re.compile(r"^[0-9a-f]{64}$")
_PROJECT_ID_PATTERN = re.compile(r"^prj_[0-9a-f]{32}$")
_SAMPLE_ID_PATTERN = re.compile(r"^smp_[0-9a-f]{32}$")
_SAMPLE_REVISION_ID_PATTERN = re.compile(r"^smpr_[0-9a-f]{32}$")
_RUN_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_WORKFLOW_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,254}$")
_LOGICAL_ID_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,127}$")
_CURSOR_KEYS = frozenset({"v", "filter_sha256", "position"})
_POSITION_KEYS = frozenset(
    {"published_at", "run_id", "artifact_generation", "artifact_id"}
)
_FILTER_SORT = (
    ("published_at", "desc"),
    ("run_id", "desc"),
    ("artifact_generation", "desc"),
    ("artifact_id", "desc"),
)


def normalize_artifact_publication_project_id(value: str | None) -> str | None:
    """Normalize an optional opaque Project filter."""
    return _optional_token(value, _PROJECT_ID_PATTERN, "project_id")


def normalize_artifact_publication_run_id(value: str | None) -> str | None:
    """Normalize an optional opaque Run filter."""
    return _optional_token(value, _RUN_ID_PATTERN, "run_id")


def normalize_artifact_publication_workflow_id(value: str | None) -> str | None:
    """Normalize an optional workflow filter."""
    return _optional_token(value, _WORKFLOW_ID_PATTERN, "workflow_id")


def normalize_artifact_publication_output_type(value: str | None) -> str | None:
    """Normalize an optional exact output-type filter."""
    return _optional_token(value, _LOGICAL_ID_PATTERN, "output_type")


def normalize_artifact_publication_artifact_id(value: str | None) -> str | None:
    """Normalize an optional opaque artifact identity."""
    return _optional_token(value, _LOGICAL_ID_PATTERN, "artifact_id")


def normalize_associated_run_sample_revision_id(
    value: str | None,
) -> str | None:
    """Normalize an optional exact Run/SampleRevision association filter."""
    return _optional_token(
        value,
        _SAMPLE_REVISION_ID_PATTERN,
        "associated_run_sample_revision_id",
    )


@dataclass(frozen=True)
class AssociatedRunSample:
    """One ordered SampleRevision associated with the publication's Run."""

    sample_id: str
    sample_revision_id: str
    revision_number: int
    ordinal: int

    def __post_init__(self) -> None:
        _required_token(self.sample_id, _SAMPLE_ID_PATTERN, "sample_id")
        _required_token(
            self.sample_revision_id,
            _SAMPLE_REVISION_ID_PATTERN,
            "sample_revision_id",
        )
        if (
            isinstance(self.revision_number, bool)
            or not isinstance(self.revision_number, int)
            or self.revision_number < 1
        ):
            raise ValueError("revision_number must be a positive integer")
        if (
            isinstance(self.ordinal, bool)
            or not isinstance(self.ordinal, int)
            or self.ordinal < 0
        ):
            raise ValueError("ordinal must be a non-negative integer")


@dataclass(frozen=True)
class ArtifactPublicationRunSampleBinding:
    """Run-wide sample evidence; this is not per-artifact attribution."""

    binding_mode: BindingMode | str
    provenance: BindingProvenance | str
    associated_run_samples: tuple[AssociatedRunSample, ...]

    def __post_init__(self) -> None:
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
            raise ValueError(
                "run sample binding mode or provenance is invalid"
            ) from None
        samples = tuple(self.associated_run_samples)
        if any(not isinstance(item, AssociatedRunSample) for item in samples):
            raise ValueError("associated_run_samples must contain sample evidence")
        if any(item.ordinal != ordinal for ordinal, item in enumerate(samples)):
            raise ValueError("associated Run samples must use contiguous ordinals")
        revision_ids = tuple(item.sample_revision_id for item in samples)
        if len(set(revision_ids)) != len(revision_ids):
            raise ValueError("associated Run samples contain duplicate revisions")
        if binding_mode is BindingMode.LEGACY_V1:
            if provenance is not BindingProvenance.UNRESOLVED or samples:
                raise ValueError("legacy Run sample evidence is inconsistent")
        elif provenance is not BindingProvenance.RESOLVED or not samples:
            raise ValueError("bound Run sample evidence is inconsistent")
        object.__setattr__(self, "binding_mode", binding_mode)
        object.__setattr__(self, "provenance", provenance)
        object.__setattr__(self, "associated_run_samples", samples)


@dataclass(frozen=True)
class ArtifactPublicationFilters:
    """Normalized relational filters for public publication queries."""

    project_id: str | None = None
    run_id: str | None = None
    workflow_id: str | None = None
    output_type: str | None = None
    associated_run_sample_revision_id: str | None = None
    published_from: datetime | None = None
    published_before: datetime | None = None
    current_only: bool = True

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "project_id",
            normalize_artifact_publication_project_id(self.project_id),
        )
        object.__setattr__(
            self,
            "run_id",
            normalize_artifact_publication_run_id(self.run_id),
        )
        object.__setattr__(
            self,
            "workflow_id",
            normalize_artifact_publication_workflow_id(self.workflow_id),
        )
        object.__setattr__(
            self,
            "output_type",
            normalize_artifact_publication_output_type(self.output_type),
        )
        object.__setattr__(
            self,
            "associated_run_sample_revision_id",
            normalize_associated_run_sample_revision_id(
                self.associated_run_sample_revision_id
            ),
        )
        object.__setattr__(
            self,
            "published_from",
            _optional_utc(self.published_from, "published_from"),
        )
        object.__setattr__(
            self,
            "published_before",
            _optional_utc(self.published_before, "published_before"),
        )
        if not isinstance(self.current_only, bool):
            raise ValueError("current_only must be a boolean")
        if (
            self.published_from is not None
            and self.published_before is not None
            and self.published_before <= self.published_from
        ):
            raise ValueError("publication date range must be increasing")

    @property
    def sha256(self) -> str:
        """Return the canonical digest binding all filters and sort semantics."""
        payload = {
            "associated_run_sample_revision_id": (
                self.associated_run_sample_revision_id
            ),
            "current_only": self.current_only,
            "output_type": self.output_type,
            "project_id": self.project_id,
            "published_before": _optional_canonical_datetime(self.published_before),
            "published_from": _optional_canonical_datetime(self.published_from),
            "run_id": self.run_id,
            "sort": _FILTER_SORT,
            "v": 1,
            "workflow_id": self.workflow_id,
        }
        raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        return sha256(raw).hexdigest()


@dataclass(frozen=True)
class ArtifactPublicationCursorPosition:
    """Stable descending keyset position for one indexed publication."""

    published_at: datetime
    run_id: str
    artifact_generation: str
    artifact_id: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "published_at",
            _required_utc(self.published_at, "published_at"),
        )
        _required_token(self.run_id, _RUN_ID_PATTERN, "run_id")
        validate_artifact_generation(self.artifact_generation)
        _required_token(self.artifact_id, _LOGICAL_ID_PATTERN, "artifact_id")

    @property
    def key(self) -> tuple[datetime, str, str, str]:
        """Return the exact repository ordering key."""
        return (
            self.published_at,
            self.run_id,
            self.artifact_generation,
            self.artifact_id,
        )


@dataclass(frozen=True)
class ArtifactPublicationCursor:
    """Decoded keyset position bound to normalized query semantics."""

    filter_sha256: str
    position: ArtifactPublicationCursorPosition

    def __post_init__(self) -> None:
        if (
            not isinstance(self.filter_sha256, str)
            or _DIGEST_PATTERN.fullmatch(self.filter_sha256) is None
        ):
            raise ValueError("artifact publication cursor filter digest is invalid")
        if not isinstance(self.position, ArtifactPublicationCursorPosition):
            raise ValueError("artifact publication cursor position is invalid")


@dataclass(frozen=True)
class ArtifactPublicationSummary:
    """Allowlisted public projection of one indexed publication revision."""

    run_id: str
    project_id: str
    workflow_id: str
    artifact_id: str
    output_type: str
    artifact_generation: str
    artifact_revision: str
    published_at: datetime
    current_artifact_generation: str
    run_sample_binding: ArtifactPublicationRunSampleBinding

    def __post_init__(self) -> None:
        _required_token(self.run_id, _RUN_ID_PATTERN, "run_id")
        _required_token(self.project_id, _PROJECT_ID_PATTERN, "project_id")
        _required_token(self.workflow_id, _WORKFLOW_ID_PATTERN, "workflow_id")
        _required_token(self.artifact_id, _LOGICAL_ID_PATTERN, "artifact_id")
        _required_token(self.output_type, _LOGICAL_ID_PATTERN, "output_type")
        validate_artifact_generation(self.artifact_generation)
        validate_artifact_revision(self.artifact_revision)
        object.__setattr__(
            self,
            "published_at",
            _required_utc(self.published_at, "published_at"),
        )
        validate_artifact_generation(self.current_artifact_generation)
        if not isinstance(self.run_sample_binding, ArtifactPublicationRunSampleBinding):
            raise ValueError("run_sample_binding is invalid")
        is_legacy_project = self.project_id == LEGACY_PROJECT_ID
        is_legacy_binding = (
            self.run_sample_binding.binding_mode is BindingMode.LEGACY_V1
        )
        if is_legacy_project is not is_legacy_binding:
            raise ValueError("publication Project and Run sample binding differ")

    @property
    def generation_status(self) -> str:
        """Derive whether this publication belongs to the Run's current generation."""
        if self.artifact_generation == self.current_artifact_generation:
            return "current"
        return "superseded"

    @property
    def cursor_position(self) -> ArtifactPublicationCursorPosition:
        """Return the public compound identity's stable list position."""
        return ArtifactPublicationCursorPosition(
            published_at=self.published_at,
            run_id=self.run_id,
            artifact_generation=self.artifact_generation,
            artifact_id=self.artifact_id,
        )


@dataclass(frozen=True)
class ArtifactPublicationPage:
    """One bounded public publication page and an optional next cursor."""

    artifact_publications: tuple[ArtifactPublicationSummary, ...]
    next_cursor: str | None

    def __post_init__(self) -> None:
        if not isinstance(self.artifact_publications, tuple) or any(
            not isinstance(item, ArtifactPublicationSummary)
            for item in self.artifact_publications
        ):
            raise ValueError("publication page contains invalid summaries")
        if self.next_cursor is not None and not isinstance(self.next_cursor, str):
            raise ValueError("publication next_cursor must be a string or None")


class ArtifactPublicationRepository(Protocol):
    """Persistence contract consumed by the publication query service."""

    def list_artifact_publications(
        self,
        *,
        filters: ArtifactPublicationFilters,
        after: ArtifactPublicationCursorPosition | None,
        limit: int,
    ) -> tuple[ArtifactPublicationSummary, ...]:
        """Return descending rows, validating any boundary against filters.

        If ``after`` no longer identifies an indexed row satisfying ``filters``,
        implementations must raise ``KeyError`` rather than silently continuing.
        """
        ...

    def get_artifact_publication(
        self,
        *,
        run_id: str,
        artifact_id: str,
        artifact_generation: str,
    ) -> ArtifactPublicationSummary:
        """Return one exact public compound identity or raise ``KeyError``."""
        ...


def encode_artifact_publication_cursor(
    position: ArtifactPublicationCursorPosition,
    *,
    filters: ArtifactPublicationFilters,
) -> str:
    """Encode a canonical cursor containing only digest-bound public state."""
    if not isinstance(position, ArtifactPublicationCursorPosition) or not isinstance(
        filters, ArtifactPublicationFilters
    ):
        raise ValueError("artifact publication cursor state is invalid")
    payload = {
        "filter_sha256": filters.sha256,
        "position": {
            "artifact_generation": position.artifact_generation,
            "artifact_id": position.artifact_id,
            "published_at": _canonical_datetime(position.published_at),
            "run_id": position.run_id,
        },
        "v": 1,
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    body = base64.urlsafe_b64encode(raw).decode("ascii").rstrip("=")
    encoded = f"{ARTIFACT_PUBLICATION_CURSOR_PREFIX}{body}"
    if len(encoded) > ARTIFACT_PUBLICATION_CURSOR_MAX_LENGTH:
        raise ValueError("artifact publication cursor exceeds its public bound")
    return encoded


def decode_artifact_publication_cursor(
    value: str,
    *,
    filters: ArtifactPublicationFilters,
) -> ArtifactPublicationCursor:
    """Strictly decode, canonicalize, and bind a public publication cursor."""
    if not isinstance(filters, ArtifactPublicationFilters):
        raise ValueError("artifact publication cursor filters are invalid")
    if (
        not isinstance(value, str)
        or len(value) <= len(ARTIFACT_PUBLICATION_CURSOR_PREFIX)
        or len(value) > ARTIFACT_PUBLICATION_CURSOR_MAX_LENGTH
        or not value.startswith(ARTIFACT_PUBLICATION_CURSOR_PREFIX)
    ):
        raise ValueError("artifact publication cursor is invalid")
    body = value[len(ARTIFACT_PUBLICATION_CURSOR_PREFIX) :]
    if _CURSOR_BODY_PATTERN.fullmatch(body) is None:
        raise ValueError("artifact publication cursor is invalid")
    padding = "=" * (-len(body) % 4)
    try:
        raw = base64.b64decode(body + padding, altchars=b"-_", validate=True)
        payload = json.loads(raw.decode("utf-8"))
    except (binascii.Error, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("artifact publication cursor is invalid") from exc
    if not isinstance(payload, dict) or set(payload) != _CURSOR_KEYS:
        raise ValueError("artifact publication cursor payload is invalid")
    if payload["v"] != 1 or isinstance(payload["v"], bool):
        raise ValueError("artifact publication cursor version is invalid")
    filter_sha256 = payload["filter_sha256"]
    position_payload = payload["position"]
    if (
        not isinstance(filter_sha256, str)
        or _DIGEST_PATTERN.fullmatch(filter_sha256) is None
        or not isinstance(position_payload, dict)
        or set(position_payload) != _POSITION_KEYS
        or not isinstance(position_payload["published_at"], str)
    ):
        raise ValueError("artifact publication cursor payload is invalid")
    try:
        published_at = datetime.fromisoformat(
            position_payload["published_at"].replace("Z", "+00:00")
        )
        position = ArtifactPublicationCursorPosition(
            published_at=published_at,
            run_id=position_payload["run_id"],
            artifact_generation=position_payload["artifact_generation"],
            artifact_id=position_payload["artifact_id"],
        )
        cursor = ArtifactPublicationCursor(
            filter_sha256=filter_sha256,
            position=position,
        )
    except (TypeError, ValueError) as exc:
        raise ValueError("artifact publication cursor payload is invalid") from exc
    if encode_artifact_publication_cursor(position, filters=filters) != value:
        raise ValueError("artifact publication cursor is not canonical")
    if cursor.filter_sha256 != filters.sha256:
        raise ValueError("artifact publication cursor filters do not match")
    return cursor


def _required_token(
    value: str,
    pattern: re.Pattern[str],
    name: str,
) -> str:
    if not isinstance(value, str) or pattern.fullmatch(value) is None:
        raise ValueError(f"artifact publication {name} is invalid")
    return value


def _optional_token(
    value: str | None,
    pattern: re.Pattern[str],
    name: str,
) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise ValueError(f"artifact publication {name} filter is invalid")
    normalized = value.strip()
    _required_token(normalized, pattern, name)
    return normalized


def _required_utc(value: datetime, name: str) -> datetime:
    if (
        not isinstance(value, datetime)
        or value.tzinfo is None
        or value.utcoffset() is None
    ):
        raise ValueError(f"artifact publication {name} must be timezone-aware")
    return value.astimezone(timezone.utc)


def _optional_utc(value: datetime | None, name: str) -> datetime | None:
    return None if value is None else _required_utc(value, name)


def _canonical_datetime(value: datetime) -> str:
    return (
        _required_utc(value, "cursor timestamp")
        .isoformat(timespec="microseconds")
        .replace("+00:00", "Z")
    )


def _optional_canonical_datetime(value: datetime | None) -> str | None:
    return None if value is None else _canonical_datetime(value)
