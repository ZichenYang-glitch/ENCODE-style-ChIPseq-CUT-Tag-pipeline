"""Workflow-neutral contracts for restricted authoring drafts and proposals.

The browser and bridge share these identities. HTTP status mapping stays in the
API. Tool envelopes use the same reason codes plus the tool-only codes documented
on ``TOOL_REASON_CODES``.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
import hashlib
import json
import re
from typing import Any

from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.snapshots import canonical_workflow_inputs_json


CONTRACT_VERSION = "hw-a1/1"
PROJECTION_VERSION = "hw-a1/1"
SERVER_NAME = "helixweave_a1"
SCHEMA_TOOL = "get_authoring_schema"
DRAFT_TOOL = "get_draft_summary"
CANONICAL_SCHEMA_TOOL = f"mcp__{SERVER_NAME}__{SCHEMA_TOOL}"
CANONICAL_DRAFT_TOOL = f"mcp__{SERVER_NAME}__{DRAFT_TOOL}"
META_SEARCH_TOOL = "SearchExtraTools"
META_EXECUTE_TOOL = "ExecuteExtraTool"

DRAFT_TTL = timedelta(days=7)
PROPOSAL_TTL = timedelta(minutes=30)
EXPIRED_RECORD_RETENTION = timedelta(days=7)

MAX_AUTHORING_REQUEST_BYTES = 256 * 1024
MAX_QUESTION_BYTES = 2 * 1024
MAX_PROJECTION_BYTES = 32 * 1024
MAX_TOOL_RESULT_BYTES = 16 * 1024
MAX_MODEL_ANSWER_BYTES = 8 * 1024
MAX_PROPOSAL_FIELDS = 16
DEFAULT_LIST_LIMIT = 20
MAX_LIST_LIMIT = 50

DRAFT_EXPIRED = "DRAFT_EXPIRED"
PROPOSAL_EXPIRED = "PROPOSAL_EXPIRED"
PROPOSAL_INVALIDATED = "PROPOSAL_INVALIDATED"
STALE_DRAFT = "STALE_DRAFT"
CONNECTION_OWNED_BY_OTHER_LOGIN = "CONNECTION_OWNED_BY_OTHER_LOGIN"
BRIDGE_UNAVAILABLE = "BRIDGE_UNAVAILABLE"
ADMISSION_FAILED = "ADMISSION_FAILED"
AUTHENTICATION_REQUIRED = "AUTHENTICATION_REQUIRED"
NOT_FOUND = "NOT_FOUND"
FORBIDDEN = "FORBIDDEN"
INVALID_REQUEST = "INVALID_REQUEST"
REQUEST_TOO_LARGE = "REQUEST_TOO_LARGE"
DECISION_CONFLICT = "DECISION_CONFLICT"
TURN_REQUEST_CONFLICT = "TURN_REQUEST_CONFLICT"
TURN_RESULT_CONFLICT = "TURN_RESULT_CONFLICT"
CONNECTION_NOT_READY = "CONNECTION_NOT_READY"
CONNECTION_LOST = "CONNECTION_LOST"
CAPABILITY_UNAVAILABLE = "CAPABILITY_UNAVAILABLE"
WORKFLOW_UNAVAILABLE = "WORKFLOW_UNAVAILABLE"
MODEL_OUTPUT_REJECTED = "MODEL_OUTPUT_REJECTED"
AUTH_SESSION_ENDED = "AUTH_SESSION_ENDED"
CLOSED_BY_OWNER = "CLOSED_BY_OWNER"
SAMPLE_SELECTION_REQUIRED = "SAMPLE_SELECTION_REQUIRED"
WORKFLOW_SAMPLE_PLAN_REQUIRED = "WORKFLOW_SAMPLE_PLAN_REQUIRED"
WORKFLOW_SAMPLE_PLAN_INVALID = "WORKFLOW_SAMPLE_PLAN_INVALID"
WORKFLOW_SAMPLE_INPUT_UNAVAILABLE = "WORKFLOW_SAMPLE_INPUT_UNAVAILABLE"
SCOPE_DENIED = "SCOPE_DENIED"
CAPABILITY_CHANGED = "CAPABILITY_CHANGED"

TOOL_REASON_CODES = frozenset(
    {SCOPE_DENIED, STALE_DRAFT, CAPABILITY_CHANGED, MODEL_OUTPUT_REJECTED}
)
SAFE_DETAIL_KEYS = frozenset(
    {"proposal_state", "draft_state", "draft_revision", "connection_state"}
)

_HTTP_STATUS = {
    DRAFT_EXPIRED: 410,
    PROPOSAL_EXPIRED: 410,
    PROPOSAL_INVALIDATED: 409,
    STALE_DRAFT: 409,
    CONNECTION_OWNED_BY_OTHER_LOGIN: 409,
    DECISION_CONFLICT: 409,
    TURN_REQUEST_CONFLICT: 409,
    TURN_RESULT_CONFLICT: 409,
    CONNECTION_NOT_READY: 409,
    CONNECTION_LOST: 409,
    BRIDGE_UNAVAILABLE: 503,
    ADMISSION_FAILED: 503,
    CAPABILITY_UNAVAILABLE: 503,
    CAPABILITY_CHANGED: 409,
    SCOPE_DENIED: 403,
    AUTHENTICATION_REQUIRED: 401,
    NOT_FOUND: 404,
    WORKFLOW_UNAVAILABLE: 404,
    FORBIDDEN: 403,
    INVALID_REQUEST: 422,
    REQUEST_TOO_LARGE: 413,
    MODEL_OUTPUT_REJECTED: 409,
}

_MESSAGES = {
    DRAFT_EXPIRED: "The authoring draft has expired.",
    PROPOSAL_EXPIRED: "The proposal adoption window has expired.",
    PROPOSAL_INVALIDATED: "The proposal can no longer be adopted.",
    STALE_DRAFT: "The draft revision does not match the server.",
    CONNECTION_OWNED_BY_OTHER_LOGIN: "Another login owns the active connection.",
    BRIDGE_UNAVAILABLE: "The authoring assistant is unavailable.",
    ADMISSION_FAILED: "The authoring assistant did not pass admission.",
    AUTHENTICATION_REQUIRED: "Authentication is required.",
    NOT_FOUND: "The authoring record was not found.",
    FORBIDDEN: "The authoring record is not available to this login.",
    INVALID_REQUEST: "The authoring request is invalid.",
    REQUEST_TOO_LARGE: "The authoring request exceeds the size limit.",
    DECISION_CONFLICT: "The proposal decision conflicts with the stored state.",
    TURN_REQUEST_CONFLICT: "The turn id was already used with different content.",
    TURN_RESULT_CONFLICT: "The turn already has a different result.",
    CONNECTION_NOT_READY: "The authoring connection is not ready.",
    CONNECTION_LOST: "The authoring connection closed before the turn finished.",
    CAPABILITY_UNAVAILABLE: "This workflow has no authoring assistant projection.",
    WORKFLOW_UNAVAILABLE: "The workflow is not available.",
    MODEL_OUTPUT_REJECTED: "The assistant response was rejected.",
    AUTH_SESSION_ENDED: "The login that owned the connection is no longer valid.",
    CLOSED_BY_OWNER: "The connection was closed.",
    SAMPLE_SELECTION_REQUIRED: (
        "Select at least one sample revision before validation."
    ),
    WORKFLOW_SAMPLE_PLAN_REQUIRED: (
        "Each selected sample revision needs an explicit workflow sample plan."
    ),
    WORKFLOW_SAMPLE_PLAN_INVALID: (
        "The workflow sample plan does not match the selected revisions."
    ),
    WORKFLOW_SAMPLE_INPUT_UNAVAILABLE: (
        "A selected input revision could not be verified."
    ),
    SCOPE_DENIED: "The assistant tool is outside the authorized scope.",
    CAPABILITY_CHANGED: "The authoring projection changed.",
}


class AuthoringRejected(Exception):
    """Stable authoring failure. ``detail`` is an allowlisted public hint."""

    def __init__(
        self,
        reason_code: str,
        *,
        detail: Mapping[str, object] | None = None,
    ) -> None:
        if reason_code not in _MESSAGES:
            raise ValueError("unknown authoring reason")
        safe: dict[str, object] = {}
        for key, value in dict(detail or {}).items():
            if key not in SAFE_DETAIL_KEYS:
                continue
            if isinstance(value, (str, int)) and not isinstance(value, bool):
                safe[key] = value
        self.reason_code = reason_code
        self.detail = safe
        super().__init__(reason_code)

    @property
    def status_code(self) -> int:
        return _HTTP_STATUS[self.reason_code]

    @property
    def message(self) -> str:
        return _MESSAGES[self.reason_code]


def authoring_message(reason_code: str) -> str:
    """Return the fixed public sentence for one reason code."""

    try:
        return _MESSAGES[reason_code]
    except KeyError:
        raise AuthoringRejected(INVALID_REQUEST) from None


@dataclass(frozen=True)
class AuthoringActor:
    """Current browser login. The digest never leaves the server."""

    user_id: str
    session_digest: str


@dataclass(frozen=True)
class FieldConstraint:
    """Approved scalar constraint. Schema defaults are intentionally absent."""

    value_type: str
    minimum: int | None = None
    maximum: int | None = None

    def public_dict(self) -> dict[str, object]:
        payload: dict[str, object] = {"type": self.value_type}
        if self.minimum is not None:
            payload["minimum"] = self.minimum
        if self.maximum is not None:
            payload["maximum"] = self.maximum
        return payload


@dataclass(frozen=True)
class AuthoringField:
    """One resource field the assistant may read and propose."""

    token: str
    label: str
    document: str
    key: str
    constraint: FieldConstraint

    def public_dict(self) -> dict[str, object]:
        return {
            "constraints": self.constraint.public_dict(),
            "label": self.label,
            "token": self.token,
            "type": self.constraint.value_type,
        }


@dataclass(frozen=True)
class IssueLabel:
    code: str
    label: str

    def public_dict(self) -> dict[str, str]:
        return {"code": self.code, "label": self.label}


def utc_now(clock: Any | None = None) -> datetime:
    """Return a whole-second UTC timestamp from an injected or system clock."""

    value = datetime.now(timezone.utc) if clock is None else clock()
    if (
        not isinstance(value, datetime)
        or value.tzinfo is None
        or value.utcoffset() is None
    ):
        raise ValueError("clock must return a timezone-aware datetime")
    return value.astimezone(timezone.utc).replace(microsecond=0)


def active_until(now: datetime, expires_at: datetime) -> bool:
    """Return whether ``now`` is strictly before ``expires_at``."""

    return utc_now(lambda: now) < utc_now(lambda: expires_at)


def canonical_json(value: object) -> str:
    """Serialize JSON with sorted keys and no insignificant whitespace."""

    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def parse_strict_json_object(payload: object, *, max_bytes: int) -> dict[str, Any]:
    """Accept one object. Reject duplicate keys, extra text, and oversize input."""

    if isinstance(payload, (bytes, bytearray)):
        raw_bytes = bytes(payload)
        text = raw_bytes.decode("utf-8")
    elif isinstance(payload, str):
        text = payload
        raw_bytes = text.encode("utf-8")
    else:
        raise AuthoringRejected(MODEL_OUTPUT_REJECTED)
    if len(raw_bytes) > max_bytes:
        raise AuthoringRejected(REQUEST_TOO_LARGE)
    try:
        decoded, end = json.JSONDecoder(
            object_pairs_hook=_reject_duplicate_keys
        ).raw_decode(text)
    except (json.JSONDecodeError, UnicodeError, ValueError):
        raise AuthoringRejected(MODEL_OUTPUT_REJECTED) from None
    if text[end:].strip():
        raise AuthoringRejected(MODEL_OUTPUT_REJECTED)
    if not isinstance(decoded, dict):
        raise AuthoringRejected(MODEL_OUTPUT_REJECTED)
    return decoded


def _reject_duplicate_keys(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    payload: dict[str, Any] = {}
    for key, value in pairs:
        if key in payload:
            raise ValueError("duplicate key")
        payload[key] = value
    return payload


_SAMPLE_REVISION_ID = re.compile(r"^smpr_[0-9a-f]{32}$")
_INPUT_FILE_REVISION_ID = re.compile(r"^inpfr_[0-9a-f]{32}$")
_SAMPLE_TOKEN = re.compile(r"^[A-Za-z0-9_.-]{1,64}$")
_WORKFLOW_SAMPLE_FIELDS = {
    "assay",
    "control_sample_revision_id",
    "layout",
    "peak_mode",
    "reads",
    "role",
    "sample_revision_id",
    "target",
}
_READ_FIELDS = {"input_file_revision_id", "read"}
_ASSAYS = {"chipseq", "cuttag", "atac", "mnase"}
_LAYOUTS = {"PE", "SE"}
_ROLES = {"treatment", "control"}
_PEAK_MODES = {"narrow", "broad", "nucleosome"}


def _normalize_workflow_samples(value: object) -> list[dict[str, Any]]:
    """Canonicalize explicit sample plans. Missing facts are rejected."""

    if value is None:
        return []
    if not isinstance(value, list) or len(value) > 256:
        raise ValueError("workflow_samples must be a bounded list")
    normalized: list[dict[str, Any]] = []
    seen_samples: set[str] = set()
    seen_reads: set[str] = set()
    for item in value:
        if not isinstance(item, Mapping) or set(item) != _WORKFLOW_SAMPLE_FIELDS:
            raise ValueError("workflow sample plan has an unexpected shape")
        sample_revision_id = item["sample_revision_id"]
        if (
            not isinstance(sample_revision_id, str)
            or _SAMPLE_REVISION_ID.fullmatch(sample_revision_id) is None
            or sample_revision_id in seen_samples
        ):
            raise ValueError("workflow sample revision is missing or repeated")
        seen_samples.add(sample_revision_id)
        assay = item["assay"]
        layout = item["layout"]
        role = item["role"]
        peak_mode = item["peak_mode"]
        target = item["target"]
        if (
            assay not in _ASSAYS
            or layout not in _LAYOUTS
            or role not in _ROLES
            or peak_mode not in _PEAK_MODES
            or not isinstance(target, str)
            or _SAMPLE_TOKEN.fullmatch(target) is None
        ):
            raise ValueError("workflow sample plan is missing an explicit field")
        if (assay == "atac" and peak_mode != "narrow") or (
            assay == "mnase" and (layout != "PE" or peak_mode != "nucleosome")
        ):
            raise ValueError("workflow sample assay and layout do not agree")
        if assay != "mnase" and peak_mode == "nucleosome":
            raise ValueError("workflow sample assay and peak mode do not agree")
        control = item["control_sample_revision_id"]
        if role == "control" and control is not None:
            raise ValueError("a control sample cannot name another control")
        if control is not None and (
            not isinstance(control, str)
            or _SAMPLE_REVISION_ID.fullmatch(control) is None
            or control == sample_revision_id
        ):
            raise ValueError("control sample revision is not explicit")
        reads = item["reads"]
        if not isinstance(reads, list):
            raise ValueError("workflow sample reads must be explicit")
        expected_reads = (1, 2) if layout == "PE" else (1,)
        observed: list[dict[str, Any]] = []
        numbers: list[int] = []
        for read in reads:
            if not isinstance(read, Mapping) or set(read) != _READ_FIELDS:
                raise ValueError("workflow sample read has an unexpected shape")
            number = read["read"]
            revision_id = read["input_file_revision_id"]
            if (
                number not in {1, 2}
                or not isinstance(revision_id, str)
                or _INPUT_FILE_REVISION_ID.fullmatch(revision_id) is None
                or revision_id in seen_reads
            ):
                raise ValueError("workflow sample read is missing or repeated")
            seen_reads.add(revision_id)
            numbers.append(number)
            observed.append({"input_file_revision_id": revision_id, "read": number})
        if tuple(sorted(numbers)) != expected_reads:
            raise ValueError("workflow sample reads do not match the layout")
        normalized.append(
            {
                "assay": assay,
                "control_sample_revision_id": control,
                "layout": layout,
                "peak_mode": peak_mode,
                "reads": sorted(observed, key=lambda read: read["read"]),
                "role": role,
                "sample_revision_id": sample_revision_id,
                "target": target,
            }
        )
    controls = {
        item["sample_revision_id"] for item in normalized if item["role"] == "control"
    }
    for item in normalized:
        control = item["control_sample_revision_id"]
        if control is not None and control not in controls:
            raise ValueError("control sample revision is not a control row")
    return normalized


def normalize_draft_content(content: Mapping[str, Any]) -> dict[str, Any]:
    """Return the canonical draft document used for equality and digests."""

    inputs = WorkflowInputs(
        config=content["config"],
        samples=content["samples"],
        options=content["options"],
    )
    decoded = json.loads(canonical_workflow_inputs_json(inputs))
    selections = []
    for item in content["input_selections"]:
        selections.append(
            {
                "input_file_revision_ids": list(item["input_file_revision_ids"]),
                "input_use_key": item["input_use_key"],
                "occurrence": item["occurrence"],
            }
        )
    reference = content["reference_profile_revision_id"]
    workflow_samples = _normalize_workflow_samples(content.get("workflow_samples", []))
    sample_revision_ids = list(content["sample_revision_ids"])
    if workflow_samples and sample_revision_ids != [
        item["sample_revision_id"] for item in workflow_samples
    ]:
        raise ValueError("workflow_samples must match sample_revision_ids")
    return {
        "config": decoded["config"],
        "input_selections": selections,
        "options": decoded["options"],
        "reference_profile_revision_id": reference,
        "sample_revision_ids": sample_revision_ids,
        "samples": decoded["samples"],
        "workflow_samples": workflow_samples,
    }


def empty_draft_content() -> dict[str, Any]:
    return normalize_draft_content(
        {
            "config": {},
            "input_selections": [],
            "options": {},
            "reference_profile_revision_id": None,
            "sample_revision_ids": [],
            "samples": None,
            "workflow_samples": [],
        }
    )


def draft_digest(content: Mapping[str, Any]) -> str:
    return sha256_text(canonical_json(normalize_draft_content(content)))


def schema_identity_digest(workflow_id: str, schema_payload: Mapping[str, Any]) -> str:
    """Digest the real schema document. This value is not sent to the model."""

    return sha256_text(
        canonical_json(
            {
                "schema": schema_payload,
                "workflow_id": workflow_id,
            }
        )
    )


def projection_identity_digest(
    workflow_id: str,
    fields: tuple[AuthoringField, ...],
    issues: tuple[IssueLabel, ...],
) -> str:
    return sha256_text(
        canonical_json(
            {
                "contract_version": PROJECTION_VERSION,
                "fields": [field.public_dict() for field in fields],
                "intents": ["explain_issues", "suggest"],
                "issue_codes": [issue.public_dict() for issue in issues],
                "workflow_id": workflow_id,
            }
        )
    )


def proposal_digest(
    *,
    sid: str,
    base_revision: int,
    schema_digest: str,
    projection_digest: str,
    changes: list[dict[str, Any]],
    reason_codes: list[str],
) -> str:
    ordered = sorted(changes, key=lambda item: str(item["field_token"]))
    return sha256_text(
        canonical_json(
            {
                "base_revision": base_revision,
                "changes": ordered,
                "projection_digest": projection_digest,
                "reason_codes": sorted(reason_codes),
                "schema_digest": schema_digest,
                "sid": sid,
            }
        )
    )


def adopt_request_digest(expected_revision: int, proposal_digest_value: str) -> str:
    return sha256_text(
        canonical_json(
            {
                "expected_revision": expected_revision,
                "proposal_digest": proposal_digest_value,
            }
        )
    )


def question_digest(
    question: Mapping[str, Any],
    *,
    expected_revision: int,
    schema_digest: str,
    projection_digest: str,
) -> str:
    return sha256_text(
        canonical_json(
            {
                "expected_revision": expected_revision,
                "projection_digest": projection_digest,
                "question": question,
                "schema_digest": schema_digest,
            }
        )
    )


def check_scalar(constraint: FieldConstraint, value: object) -> None:
    """Reject values outside the approved scalar constraint."""

    if constraint.value_type != "integer":
        raise AuthoringRejected(MODEL_OUTPUT_REJECTED)
    if isinstance(value, bool) or not isinstance(value, int):
        raise AuthoringRejected(MODEL_OUTPUT_REJECTED)
    if constraint.minimum is not None and value < constraint.minimum:
        raise AuthoringRejected(MODEL_OUTPUT_REJECTED)
    if constraint.maximum is not None and value > constraint.maximum:
        raise AuthoringRejected(MODEL_OUTPUT_REJECTED)


def apply_field_changes(
    content: Mapping[str, Any],
    changes: list[dict[str, Any]],
    fields: Mapping[str, AuthoringField],
) -> dict[str, Any]:
    """Return a new draft with approved scalar assignments applied."""

    if len(changes) > MAX_PROPOSAL_FIELDS:
        raise AuthoringRejected(MODEL_OUTPUT_REJECTED)
    seen: set[str] = set()
    updated = normalize_draft_content(content)
    for change in changes:
        token = change["field_token"]
        if not isinstance(token, str) or token in seen or token not in fields:
            raise AuthoringRejected(MODEL_OUTPUT_REJECTED)
        seen.add(token)
        field = fields[token]
        check_scalar(field.constraint, change["value"])
        document = updated[field.document]
        if not isinstance(document, dict):
            raise AuthoringRejected(MODEL_OUTPUT_REJECTED)
        document[field.key] = change["value"]
    return normalize_draft_content(updated)


def projected_field(
    content: Mapping[str, Any],
    field: AuthoringField,
) -> dict[str, object]:
    document = content.get(field.document)
    if not isinstance(document, dict) or field.key not in document:
        return {"present": False, "token": field.token, "value": None}
    return {"present": True, "token": field.token, "value": document[field.key]}


def sample_count(samples: object) -> int:
    """Count samples without copying their identifiers into a projection."""

    if samples is None:
        return 0
    if isinstance(samples, list):
        return len(samples)
    if isinstance(samples, str):
        return sum(1 for line in samples.splitlines() if line.strip())
    return 0


def tool_envelope(
    *,
    ok: bool,
    data: dict[str, Any] | None,
    issues: list[dict[str, str]] | None = None,
) -> dict[str, Any]:
    envelope = {
        "contract_version": CONTRACT_VERSION,
        "data": data,
        "issues": list(issues or []),
        "ok": ok,
    }
    encoded = canonical_json(envelope).encode("utf-8")
    if len(encoded) > MAX_TOOL_RESULT_BYTES:
        return {
            "contract_version": CONTRACT_VERSION,
            "data": None,
            "issues": [
                {
                    "code": CAPABILITY_CHANGED,
                    "retry_class": "fatal",
                    "severity": "error",
                }
            ],
            "ok": False,
        }
    return envelope


def tool_issue(code: str) -> dict[str, str]:
    retry = "stale" if code == STALE_DRAFT else "fatal"
    return {"code": code, "retry_class": retry, "severity": "error"}
