"""Public-safe identities for durable artifact and QC result generations."""

from __future__ import annotations

import base64
import binascii
from dataclasses import dataclass
from hashlib import sha256
import json
import re
from secrets import token_hex
from typing import TYPE_CHECKING, Any, Iterable

if TYPE_CHECKING:
    from encode_pipeline.platform.runs import RunArtifactRef, RunQcMetric


ARTIFACT_REVISION_PREFIX = "artifactrev-"
ARTIFACT_PATH_IDENTITY_PREFIX = "artifactpath-"
ARTIFACT_PATH_IDENTITY_METADATA_KEY = "_platform_path_identity"
BOUNDED_ARTIFACT_CONTENT_REVISION_MAX_BYTES = 16 * 1024 * 1024
ARTIFACT_GENERATION_PREFIX = "artifactgen-"
QC_GENERATION_PREFIX = "qcgen-"
RESULT_ATTEMPT_PREFIX = "resultattempt-"
ARTIFACT_CURSOR_PREFIX = "artifactcur_"
QC_CURSOR_PREFIX = "qccur_"
RESULT_CURSOR_MAX_LENGTH = 1024

_HEX64 = re.compile(r"^[0-9a-f]{64}$")
_CURSOR_BODY = re.compile(r"^[A-Za-z0-9_-]+$")
_RUN_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_ARTIFACT_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,127}$")
_METRIC_ID = re.compile(r"^qcmetric-[0-9a-f]{64}$")
_ARTIFACT_CURSOR_KEYS = frozenset(
    {"v", "run_id", "artifact_generation", "after_artifact_id"}
)
_QC_CURSOR_KEYS = frozenset({"v", "run_id", "qc_generation", "after_metric_id"})


def new_result_attempt_id() -> str:
    """Return a caller-owned, disclosure-safe result attempt identity."""
    return f"{RESULT_ATTEMPT_PREFIX}{token_hex(32)}"


def validate_result_attempt_id(value: str) -> str:
    """Validate and return one result attempt identity."""
    return _validate_prefixed_digest(value, RESULT_ATTEMPT_PREFIX, "result attempt")


def validate_artifact_revision(value: str) -> str:
    """Validate and return one persisted artifact revision identity."""
    return _validate_prefixed_digest(
        value, ARTIFACT_REVISION_PREFIX, "artifact revision"
    )


def validate_artifact_path_identity(value: str) -> str:
    """Validate and return one disclosure-safe path descriptor identity."""
    return _validate_prefixed_digest(
        value, ARTIFACT_PATH_IDENTITY_PREFIX, "artifact path identity"
    )


def build_artifact_path_identity(
    *,
    parent_identities: Iterable[tuple[int, int, int]],
    file_identity: tuple[int, int, int, int, int, int, int, int, int],
) -> str:
    """Hash one descriptor-relative path chain without exposing host details."""
    parents = tuple(parent_identities)
    if not parents:
        raise ValueError("artifact path identity requires a parent chain")
    if any(not _valid_descriptor_tuple(value, 3) for value in parents):
        raise ValueError("artifact parent descriptor identity is invalid")
    if not _valid_descriptor_tuple(file_identity, 9):
        raise ValueError("artifact file descriptor identity is invalid")
    digest = _digest_json(
        {
            "file": file_identity,
            "parents": parents,
            "scheme": "descriptor-relative-path-sha256-v1",
        }
    )
    return f"{ARTIFACT_PATH_IDENTITY_PREFIX}{digest}"


def validate_artifact_generation(value: str) -> str:
    """Validate and return one public artifact generation token."""
    return _validate_prefixed_digest(
        value, ARTIFACT_GENERATION_PREFIX, "artifact generation"
    )


def validate_qc_generation(value: str) -> str:
    """Validate and return one public QC generation token."""
    return _validate_prefixed_digest(value, QC_GENERATION_PREFIX, "QC generation")


def build_artifact_content_revision(
    *,
    output_type: str,
    relative_path: str,
    content: bytes,
) -> str:
    """Bind one bounded QC source revision to its exact bytes."""
    if not isinstance(output_type, str) or not output_type:
        raise ValueError("artifact output_type must be a non-empty string")
    if not isinstance(relative_path, str) or not relative_path:
        raise ValueError("artifact relative_path must be a non-empty string")
    if not isinstance(content, bytes):
        raise ValueError("artifact content must be bytes")
    digest = _digest_json(
        {
            "content_sha256": sha256(content).hexdigest(),
            "output_type": output_type,
            "relative_path": relative_path,
            "scheme": "bounded-content-sha256-v1",
        }
    )
    return f"{ARTIFACT_REVISION_PREFIX}{digest}"


def build_artifact_descriptor_revision(
    *,
    output_type: str,
    relative_path: str,
    descriptor_identity: str,
) -> str:
    """Bind a large artifact to a disclosure-safe descriptor identity.

    The caller owns capture and validation of ``descriptor_identity``. This is a
    liveness/revision identity, not a full-content proof.
    """
    for name, value in (
        ("output_type", output_type),
        ("relative_path", relative_path),
        ("descriptor_identity", descriptor_identity),
    ):
        if not isinstance(value, str) or not value:
            raise ValueError(f"artifact {name} must be a non-empty string")
    digest = _digest_json(
        {
            "descriptor_identity": descriptor_identity,
            "output_type": output_type,
            "relative_path": relative_path,
            "scheme": "descriptor-revision-sha256-v1",
        }
    )
    return f"{ARTIFACT_REVISION_PREFIX}{digest}"


def artifact_manifest_digest(artifacts: Iterable[RunArtifactRef]) -> str:
    """Return the canonical digest of one complete artifact manifest."""
    values = tuple(artifacts)
    ordered = tuple(sorted(values, key=lambda artifact: artifact.artifact_id))
    if len({artifact.artifact_id for artifact in ordered}) != len(ordered):
        raise ValueError("artifact manifest contains a duplicate artifact_id")
    payload: list[dict[str, Any]] = []
    for artifact in ordered:
        validate_artifact_revision(artifact.revision)
        payload.append(
            {
                "artifact_id": artifact.artifact_id,
                "artifact_type": artifact.artifact_type,
                "metadata": artifact.to_dict()["metadata"],
                "mime_type": artifact.mime_type,
                "name": artifact.name,
                "produced_at": artifact.produced_at.isoformat(),
                "revision": artifact.revision,
                "run_id": artifact.run_id,
                "uri": artifact.uri,
            }
        )
    return _digest_json({"artifacts": payload, "scheme": "artifact-manifest-v1"})


def qc_metric_manifest_digest(metrics: Iterable[RunQcMetric]) -> str:
    """Return the canonical digest of one complete QC metric set."""
    values = tuple(metrics)
    ordered = tuple(sorted(values, key=lambda metric: metric.metric_id))
    if len({metric.metric_id for metric in ordered}) != len(ordered):
        raise ValueError("QC manifest contains a duplicate metric_id")
    payload = [
        {
            "assay": metric.assay,
            "display_name": metric.display_name,
            "experiment_id": metric.experiment_id,
            "metric_id": metric.metric_id,
            "metric_key": metric.metric_key,
            "produced_at": metric.produced_at.isoformat(),
            "qc_flag": metric.qc_flag,
            "run_id": metric.run_id,
            "sample_id": metric.sample_id,
            "scope": metric.scope,
            "source_artifact_id": metric.source_artifact_id,
            "unit": metric.unit,
            "value": _canonical_decimal(metric.value),
        }
        for metric in ordered
    ]
    return _digest_json({"metrics": payload, "scheme": "qc-metric-manifest-v1"})


def build_artifact_generation(
    *, run_id: str, revision: int, artifacts: Iterable[RunArtifactRef]
) -> str:
    """Return an ABA-safe artifact generation for a monotonic revision."""
    _validate_run_id(run_id)
    _validate_positive_revision(revision)
    digest = _digest_json(
        {
            "manifest_digest": artifact_manifest_digest(artifacts),
            "revision": revision,
            "run_id": run_id,
            "scheme": "artifact-generation-v1",
        }
    )
    return f"{ARTIFACT_GENERATION_PREFIX}{digest}"


def build_qc_generation(
    *,
    run_id: str,
    revision: int,
    artifact_generation: str,
    metrics: Iterable[RunQcMetric],
) -> str:
    """Return an ABA-safe QC generation bound to an artifact generation."""
    _validate_run_id(run_id)
    _validate_positive_revision(revision)
    validate_artifact_generation(artifact_generation)
    digest = _digest_json(
        {
            "artifact_generation": artifact_generation,
            "manifest_digest": qc_metric_manifest_digest(metrics),
            "revision": revision,
            "run_id": run_id,
            "scheme": "qc-generation-v1",
        }
    )
    return f"{QC_GENERATION_PREFIX}{digest}"


@dataclass(frozen=True)
class RunResultState:
    """Canonical persisted artifact/QC generation and attempt state."""

    run_id: str
    artifact_revision: int = 0
    artifact_generation: str | None = None
    artifact_manifest_digest: str | None = None
    artifact_attempt_id: str | None = None
    artifact_attempt_status: str | None = None
    artifact_outcome: str | None = None
    artifact_reason_code: str | None = None
    qc_revision: int = 0
    qc_generation: str | None = None
    qc_manifest_digest: str | None = None
    qc_attempt_id: str | None = None
    qc_attempt_status: str | None = None
    qc_attempt_artifact_generation: str | None = None
    qc_artifact_generation: str | None = None
    qc_outcome: str | None = None
    qc_reason_code: str | None = None

    def __post_init__(self) -> None:
        _validate_state_run_id(self.run_id)
        if any(
            value is not None
            for value in (
                self.artifact_generation,
                self.artifact_attempt_id,
                self.artifact_outcome,
                self.qc_generation,
                self.qc_attempt_id,
                self.qc_outcome,
            )
        ):
            _validate_run_id(self.run_id)
        _validate_state_revision(self.artifact_revision, "artifact")
        _validate_state_revision(self.qc_revision, "QC")
        _validate_bound_generation(
            revision=self.artifact_revision,
            generation=self.artifact_generation,
            digest=self.artifact_manifest_digest,
            generation_validator=validate_artifact_generation,
            name="artifact",
        )
        if self.qc_generation is not None:
            validate_qc_generation(self.qc_generation)
        if self.qc_manifest_digest is not None:
            _validate_hex_digest(self.qc_manifest_digest, "QC manifest")
        if (self.qc_generation is None) != (self.qc_manifest_digest is None):
            raise ValueError("QC generation and manifest digest must be paired")
        if self.qc_generation is None:
            if self.qc_artifact_generation is not None:
                raise ValueError("unbound QC state cannot name an artifact generation")
        else:
            if self.qc_revision <= 0 or self.qc_artifact_generation is None:
                raise ValueError("bound QC state requires a positive revision")
            validate_artifact_generation(self.qc_artifact_generation)
        _validate_attempt_pair(
            self.artifact_attempt_id,
            self.artifact_attempt_status,
            "artifact",
        )
        _validate_attempt_pair(self.qc_attempt_id, self.qc_attempt_status, "QC")
        if self.qc_attempt_id is None:
            if self.qc_attempt_artifact_generation is not None:
                raise ValueError("QC attempt generation requires an attempt")
        else:
            if self.qc_attempt_artifact_generation is None:
                raise ValueError("QC attempt requires an artifact generation")
            validate_artifact_generation(self.qc_attempt_artifact_generation)
        _validate_outcome(self.artifact_outcome, self.artifact_reason_code, "artifact")
        _validate_outcome(self.qc_outcome, self.qc_reason_code, "QC")


@dataclass(frozen=True)
class QcMetricCursor:
    """One QC keyset boundary bound to its run and generation."""

    run_id: str
    qc_generation: str
    after_metric_id: str

    def __post_init__(self) -> None:
        _validate_run_id(self.run_id)
        validate_qc_generation(self.qc_generation)
        if (
            not isinstance(self.after_metric_id, str)
            or _METRIC_ID.fullmatch(self.after_metric_id) is None
        ):
            raise ValueError("QC cursor metric ID is invalid")


@dataclass(frozen=True)
class ArtifactCursor:
    """One artifact keyset boundary bound to its run and generation."""

    run_id: str
    artifact_generation: str
    after_artifact_id: str

    def __post_init__(self) -> None:
        _validate_run_id(self.run_id)
        validate_artifact_generation(self.artifact_generation)
        if (
            not isinstance(self.after_artifact_id, str)
            or _ARTIFACT_ID.fullmatch(self.after_artifact_id) is None
        ):
            raise ValueError("artifact cursor artifact ID is invalid")


def encode_artifact_cursor(cursor: ArtifactCursor) -> str:
    """Encode a canonical public artifact cursor; it is not a credential."""
    if not isinstance(cursor, ArtifactCursor):
        raise ValueError("artifact cursor is invalid")
    payload = {
        "after_artifact_id": cursor.after_artifact_id,
        "artifact_generation": cursor.artifact_generation,
        "run_id": cursor.run_id,
        "v": 1,
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    body = base64.urlsafe_b64encode(raw).decode("ascii").rstrip("=")
    encoded = f"{ARTIFACT_CURSOR_PREFIX}{body}"
    if len(encoded) > RESULT_CURSOR_MAX_LENGTH:
        raise ValueError("artifact cursor exceeds its public bound")
    return encoded


def decode_artifact_cursor(
    value: str, *, run_id: str, artifact_generation: str
) -> ArtifactCursor:
    """Strictly decode and bind one artifact cursor to current result state."""
    _validate_run_id(run_id)
    validate_artifact_generation(artifact_generation)
    payload = _decode_cursor_payload(
        value,
        prefix=ARTIFACT_CURSOR_PREFIX,
        keys=_ARTIFACT_CURSOR_KEYS,
        name="artifact",
    )
    cursor = ArtifactCursor(
        run_id=payload["run_id"],
        artifact_generation=payload["artifact_generation"],
        after_artifact_id=payload["after_artifact_id"],
    )
    if encode_artifact_cursor(cursor) != value:
        raise ValueError("artifact cursor is not canonical")
    if cursor.run_id != run_id or cursor.artifact_generation != artifact_generation:
        raise ValueError("artifact cursor generation does not match")
    return cursor


def encode_qc_metric_cursor(cursor: QcMetricCursor) -> str:
    """Encode a canonical public QC cursor; it is not an auth credential."""
    if not isinstance(cursor, QcMetricCursor):
        raise ValueError("QC cursor is invalid")
    payload = {
        "after_metric_id": cursor.after_metric_id,
        "qc_generation": cursor.qc_generation,
        "run_id": cursor.run_id,
        "v": 1,
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    body = base64.urlsafe_b64encode(raw).decode("ascii").rstrip("=")
    encoded = f"{QC_CURSOR_PREFIX}{body}"
    if len(encoded) > RESULT_CURSOR_MAX_LENGTH:
        raise ValueError("QC cursor exceeds its public bound")
    return encoded


def decode_qc_metric_cursor(
    value: str, *, run_id: str, qc_generation: str
) -> QcMetricCursor:
    """Strictly decode and bind one QC cursor to current result state."""
    _validate_run_id(run_id)
    validate_qc_generation(qc_generation)
    if (
        not isinstance(value, str)
        or len(value) <= len(QC_CURSOR_PREFIX)
        or len(value) > RESULT_CURSOR_MAX_LENGTH
        or not value.startswith(QC_CURSOR_PREFIX)
    ):
        raise ValueError("QC cursor is invalid")
    body = value[len(QC_CURSOR_PREFIX) :]
    if _CURSOR_BODY.fullmatch(body) is None:
        raise ValueError("QC cursor is invalid")
    try:
        raw = base64.b64decode(
            body + "=" * (-len(body) % 4), altchars=b"-_", validate=True
        )
        payload = json.loads(raw.decode("utf-8"))
    except (binascii.Error, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("QC cursor is invalid") from exc
    if not isinstance(payload, dict) or set(payload) != _QC_CURSOR_KEYS:
        raise ValueError("QC cursor payload is invalid")
    if payload["v"] != 1 or isinstance(payload["v"], bool):
        raise ValueError("QC cursor version is invalid")
    cursor = QcMetricCursor(
        run_id=payload["run_id"],
        qc_generation=payload["qc_generation"],
        after_metric_id=payload["after_metric_id"],
    )
    if encode_qc_metric_cursor(cursor) != value:
        raise ValueError("QC cursor is not canonical")
    if cursor.run_id != run_id or cursor.qc_generation != qc_generation:
        raise ValueError("QC cursor generation does not match")
    return cursor


def _decode_cursor_payload(
    value: str,
    *,
    prefix: str,
    keys: frozenset[str],
    name: str,
) -> dict[str, Any]:
    if (
        not isinstance(value, str)
        or len(value) <= len(prefix)
        or len(value) > RESULT_CURSOR_MAX_LENGTH
        or not value.startswith(prefix)
    ):
        raise ValueError(f"{name} cursor is invalid")
    body = value[len(prefix) :]
    if _CURSOR_BODY.fullmatch(body) is None:
        raise ValueError(f"{name} cursor is invalid")
    try:
        raw = base64.b64decode(
            body + "=" * (-len(body) % 4), altchars=b"-_", validate=True
        )
        payload = json.loads(raw.decode("utf-8"))
    except (binascii.Error, UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError(f"{name} cursor is invalid") from exc
    if not isinstance(payload, dict) or set(payload) != keys:
        raise ValueError(f"{name} cursor payload is invalid")
    if payload["v"] != 1 or isinstance(payload["v"], bool):
        raise ValueError(f"{name} cursor version is invalid")
    return payload


def _validate_prefixed_digest(value: str, prefix: str, name: str) -> str:
    if (
        not isinstance(value, str)
        or not value.startswith(prefix)
        or _HEX64.fullmatch(value[len(prefix) :]) is None
    ):
        raise ValueError(f"{name} is invalid")
    return value


def _valid_descriptor_tuple(value: object, length: int) -> bool:
    return (
        isinstance(value, tuple)
        and len(value) == length
        and all(
            isinstance(item, int) and not isinstance(item, bool) and item >= 0
            for item in value
        )
    )


def _validate_run_id(value: str) -> None:
    if not isinstance(value, str) or _RUN_ID.fullmatch(value) is None:
        raise ValueError("result generation run_id is invalid")


def _validate_state_run_id(value: str) -> None:
    # A newly created run has an empty result state before any public result
    # identity exists. Preserve the platform's existing opaque run-id boundary
    # so workspace/path policy remains the authoritative rejection point; once
    # an attempt or generation is attached, the stricter public grammar applies.
    if not isinstance(value, str) or not value or len(value) > 128 or "\x00" in value:
        raise ValueError("result state run_id is invalid")


def _validate_positive_revision(value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise ValueError("result generation revision must be a positive integer")


def _validate_state_revision(value: int, name: str) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError(f"{name} revision must be a non-negative integer")


def _validate_hex_digest(value: str, name: str) -> None:
    if not isinstance(value, str) or _HEX64.fullmatch(value) is None:
        raise ValueError(f"{name} digest is invalid")


def _validate_bound_generation(
    *,
    revision: int,
    generation: str | None,
    digest: str | None,
    generation_validator,
    name: str,
) -> None:
    if (generation is None) != (digest is None):
        raise ValueError(f"{name} generation and manifest digest must be paired")
    if generation is None:
        if revision != 0:
            raise ValueError(f"unbound {name} state must have revision zero")
        return
    if revision <= 0:
        raise ValueError(f"bound {name} state requires a positive revision")
    generation_validator(generation)
    assert digest is not None
    _validate_hex_digest(digest, name)


def _validate_attempt_pair(
    attempt_id: str | None, status: str | None, name: str
) -> None:
    if (attempt_id is None) != (status is None):
        raise ValueError(f"{name} attempt ID and status must be paired")
    if attempt_id is not None:
        validate_result_attempt_id(attempt_id)
        if status not in {"pending", "succeeded", "failed"}:
            raise ValueError(f"{name} attempt status is invalid")


def _validate_outcome(outcome: str | None, reason: str | None, name: str) -> None:
    if outcome not in {None, "succeeded", "failed", "invalidated"}:
        raise ValueError(f"{name} outcome is invalid")
    if outcome == "failed":
        if not isinstance(reason, str) or not reason:
            raise ValueError(f"failed {name} outcome requires a reason")
    elif reason is not None:
        raise ValueError(f"non-failed {name} outcome cannot have a reason")


def _canonical_decimal(value: Any) -> str:
    if not hasattr(value, "is_finite") or not value.is_finite():
        raise ValueError("QC metric value must be a finite Decimal")
    text = format(value, "f")
    if "." in text:
        text = text.rstrip("0").rstrip(".")
    return "0" if text in {"", "-0"} else text


def _digest_json(payload: object) -> str:
    try:
        raw = json.dumps(
            payload,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")
    except (TypeError, ValueError) as exc:
        raise ValueError("result generation payload is not canonical JSON") from exc
    return sha256(raw).hexdigest()
