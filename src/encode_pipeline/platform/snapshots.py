"""Immutable workflow-neutral evidence for successfully validated inputs."""

from __future__ import annotations

from dataclasses import dataclass, replace
from datetime import datetime, timezone
from hashlib import sha256
import json
import math
import re
from typing import Any

from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.runs import RunRecord


PAYLOAD_DIGEST_SCHEME = "sha256-framed-workflow-inputs-v1"
VALIDATION_EVIDENCE_OUTCOME = "adapter_validation_succeeded"
MAX_SAFE_JSON_INTEGER = 2**53 - 1
_SNAPSHOT_ID = re.compile(r"^vsnap_[0-9a-f]{32}$")
_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_ISSUE_CODE = re.compile(r"^[A-Z][A-Z0-9_]{0,127}$")


def _normalize_json(value: object, *, depth: int = 0) -> Any:
    if depth > 64:
        raise ValueError("workflow inputs must be JSON-safe")
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, int):
        if abs(value) > MAX_SAFE_JSON_INTEGER:
            raise ValueError("workflow inputs must use JSON-safe integers")
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("workflow inputs must use JSON-safe numbers")
        return value
    if isinstance(value, list):
        return [_normalize_json(item, depth=depth + 1) for item in value]
    if isinstance(value, dict):
        normalized: dict[str, Any] = {}
        for key, item in value.items():
            if not isinstance(key, str):
                raise ValueError("workflow inputs must use string JSON keys")
            normalized[key] = _normalize_json(item, depth=depth + 1)
        return normalized
    raise ValueError("workflow inputs must be JSON-safe")


def canonical_workflow_inputs_json(inputs: WorkflowInputs) -> str:
    """Return deterministic JSON for one complete Config/Samples/Options value."""
    if not isinstance(inputs, WorkflowInputs):
        raise ValueError("inputs must be WorkflowInputs")
    normalized = _normalize_json(inputs.to_dict())
    return json.dumps(
        normalized,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def build_workflow_inputs_digest(canonical_payload: str) -> str:
    """Return the versioned framed digest for canonical workflow input JSON."""
    if not isinstance(canonical_payload, str) or not canonical_payload:
        raise ValueError("canonical_payload must be a non-empty string")
    digest = sha256()
    for part in (PAYLOAD_DIGEST_SCHEME, canonical_payload):
        encoded = part.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return digest.hexdigest()


def workflow_inputs_from_canonical_json(canonical_payload: str) -> WorkflowInputs:
    """Validate and decode canonical snapshot payload without trusting storage."""
    if not isinstance(canonical_payload, str) or not canonical_payload:
        raise ValueError("canonical_payload must be a non-empty string")
    try:
        decoded = json.loads(canonical_payload)
    except (json.JSONDecodeError, TypeError):
        raise ValueError("canonical_payload must contain JSON") from None
    if not isinstance(decoded, dict) or set(decoded) != {
        "config",
        "samples",
        "options",
    }:
        raise ValueError("canonical_payload must contain complete workflow inputs")
    if not isinstance(decoded["config"], dict) or not isinstance(
        decoded["options"], dict
    ):
        raise ValueError("canonical_payload config and options must be objects")
    inputs = WorkflowInputs(
        config=_normalize_json(decoded["config"]),
        samples=_normalize_json(decoded["samples"]),
        options=_normalize_json(decoded["options"]),
    )
    if canonical_workflow_inputs_json(inputs) != canonical_payload:
        raise ValueError("canonical_payload is not in canonical form")
    return inputs


def _utc(value: datetime, name: str) -> datetime:
    if (
        not isinstance(value, datetime)
        or value.tzinfo is None
        or value.utcoffset() is None
    ):
        raise ValueError(f"{name} must be timezone-aware")
    return value.astimezone(timezone.utc)


@dataclass(frozen=True)
class ValidatedInputSnapshot:
    """Durable immutable authorization created by successful adapter validation."""

    snapshot_id: str
    workflow_id: str
    adapter_version: str
    schema_version: str
    schema_dialect: str
    workflow_build_identity: WorkflowBuildIdentity
    canonical_payload: str
    payload_digest_scheme: str
    payload_digest: str
    validation_outcome: str
    validation_issue_codes: tuple[str, ...]
    validated_at: datetime
    expires_at: datetime
    consumed_run_id: str | None = None
    consumed_at: datetime | None = None

    def __post_init__(self) -> None:
        if (
            not isinstance(self.snapshot_id, str)
            or _SNAPSHOT_ID.fullmatch(self.snapshot_id) is None
        ):
            raise ValueError("snapshot_id must be an opaque validated snapshot ID")
        for name in (
            "workflow_id",
            "adapter_version",
            "schema_version",
            "schema_dialect",
        ):
            value = getattr(self, name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{name} must be a non-empty string")
            object.__setattr__(self, name, value.strip())
        identity = self.workflow_build_identity
        if not isinstance(identity, WorkflowBuildIdentity):
            raise ValueError("workflow_build_identity must be a WorkflowBuildIdentity")
        if identity.workflow_id != self.workflow_id:
            raise ValueError("workflow build identity does not match workflow_id")
        if identity.adapter_version != self.adapter_version:
            raise ValueError("workflow build identity does not match adapter_version")

        if self.payload_digest_scheme != PAYLOAD_DIGEST_SCHEME:
            raise ValueError("payload digest scheme is unsupported")
        workflow_inputs_from_canonical_json(self.canonical_payload)
        if (
            not isinstance(self.payload_digest, str)
            or _DIGEST.fullmatch(self.payload_digest) is None
            or self.payload_digest
            != build_workflow_inputs_digest(self.canonical_payload)
        ):
            raise ValueError("payload digest does not match canonical payload")
        if self.validation_outcome != VALIDATION_EVIDENCE_OUTCOME:
            raise ValueError("validation evidence must record successful validation")
        issue_codes = tuple(self.validation_issue_codes)
        if len(issue_codes) > 100 or any(
            not isinstance(code, str) or _ISSUE_CODE.fullmatch(code) is None
            for code in issue_codes
        ):
            raise ValueError("validation issue codes are invalid")
        object.__setattr__(self, "validation_issue_codes", issue_codes)

        validated_at = _utc(self.validated_at, "validated_at")
        expires_at = _utc(self.expires_at, "expires_at")
        if expires_at <= validated_at:
            raise ValueError("expires_at must be later than validated_at")
        object.__setattr__(self, "validated_at", validated_at)
        object.__setattr__(self, "expires_at", expires_at)

        if (self.consumed_run_id is None) != (self.consumed_at is None):
            raise ValueError(
                "snapshot consumption ID and time must be recorded together"
            )
        if self.consumed_run_id is not None:
            if (
                not isinstance(self.consumed_run_id, str)
                or not self.consumed_run_id.strip()
            ):
                raise ValueError("consumed_run_id must be a non-empty string")
            consumed_at = _utc(self.consumed_at, "consumed_at")  # type: ignore[arg-type]
            if consumed_at < validated_at:
                raise ValueError("consumed_at cannot precede validation")
            object.__setattr__(self, "consumed_run_id", self.consumed_run_id.strip())
            object.__setattr__(self, "consumed_at", consumed_at)

    def to_workflow_inputs(self) -> WorkflowInputs:
        """Return a fresh decoded input value after verifying canonical storage."""
        return workflow_inputs_from_canonical_json(self.canonical_payload)

    def with_consumption(
        self,
        run_id: str,
        consumed_at: datetime,
    ) -> "ValidatedInputSnapshot":
        """Return this immutable snapshot with first-use evidence attached."""
        if self.consumed_run_id is not None:
            raise ValueError("validated snapshot is already consumed")
        return replace(self, consumed_run_id=run_id, consumed_at=consumed_at)


@dataclass(frozen=True)
class ValidatedSnapshotRunCreation:
    """Canonical result of atomically consuming a snapshot for one run."""

    record: RunRecord
    created: bool

    def __post_init__(self) -> None:
        if not isinstance(self.record, RunRecord):
            raise ValueError("record must be a RunRecord")
        if not isinstance(self.created, bool):
            raise ValueError("created must be a bool")
