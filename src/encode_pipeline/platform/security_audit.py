"""Closed, disclosure-safe audit identities for security-sensitive actions."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from enum import Enum
from hashlib import sha256
import re
from secrets import token_hex

from encode_pipeline.platform.authentication import validate_user_id


_AUDIT_EVENT_ID = re.compile(r"^aevt_[0-9a-f]{32}$")
_RUN_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_ARTIFACT_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,127}$")
_REFERENCE_ID = re.compile(r"^refp_[0-9a-f]{32}$")
_STORAGE_ID = re.compile(r"^stgp_[0-9a-f]{32}$")
_RESOURCE_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_RESOURCE_DIGEST_DOMAIN = b"helixweave-security-audit-resource-v1\0"


class AuditAction(str, Enum):
    """The bounded action vocabulary required by the LAN authorization boundary."""

    LOGIN = "auth.login"
    LOGOUT = "auth.logout"
    RUN_CREATE = "run.create"
    RUN_START = "run.start"
    RUN_CANCEL = "run.cancel"
    ARTIFACT_DOWNLOAD = "artifact.download"
    REFERENCE_REGISTER = "reference.register"
    REFERENCE_ENABLE = "reference.enable"
    REFERENCE_DISABLE = "reference.disable"
    STORAGE_REGISTER = "storage.register"
    STORAGE_ARCHIVE = "storage.archive"
    ACCOUNT_CREATE = "account.create"
    ACCOUNT_ENABLE = "account.enable"
    ACCOUNT_DISABLE = "account.disable"
    ACCOUNT_PASSWORD_RESET = "account.password_reset"
    ACCOUNT_SESSIONS_REVOKE = "account.sessions_revoke"
    RUN_FAIL = "run.fail"
    RUN_REQUEUE = "run.requeue"


class AuditOutcome(str, Enum):
    """Publicly stable result categories without raw exception details."""

    SUCCEEDED = "succeeded"
    DENIED = "denied"
    FAILED = "failed"


class AuditActorKind(str, Enum):
    """Closed actor categories for browser users and local operator commands."""

    USER = "user"
    LOCAL_OPERATOR = "local_operator"
    UNAUTHENTICATED = "unauthenticated"


class AuditReasonCode(str, Enum):
    """Closed non-sensitive reasons; Stage B may deliberately extend this set."""

    INVALID_CREDENTIALS = "INVALID_CREDENTIALS"
    LOGIN_RATE_LIMITED = "LOGIN_RATE_LIMITED"
    AUTHENTICATION_REQUIRED = "AUTHENTICATION_REQUIRED"
    SESSION_INVALID = "SESSION_INVALID"
    CSRF_INVALID = "CSRF_INVALID"
    ADMINISTRATOR_REQUIRED = "ADMINISTRATOR_REQUIRED"
    ACCOUNT_DISABLED = "ACCOUNT_DISABLED"
    SETUP_REQUIRED = "SETUP_REQUIRED"
    OPERATION_CONFLICT = "OPERATION_CONFLICT"
    RESOURCE_NOT_FOUND = "RESOURCE_NOT_FOUND"
    INTERNAL_FAILURE = "INTERNAL_FAILURE"


class AuditResourceKind(str, Enum):
    """Resource categories whose identifiers are safe to retain."""

    ACCOUNT = "account"
    RUN = "run"
    ARTIFACT = "artifact"
    REFERENCE = "reference"
    STORAGE = "storage"


def new_audit_event_id() -> str:
    """Return a random path-safe identity for one durable audit event."""

    return f"aevt_{token_hex(16)}"


def _audit_event_id(value: object) -> str:
    if not isinstance(value, str) or _AUDIT_EVENT_ID.fullmatch(value) is None:
        raise ValueError("event_id must be an opaque audit event identifier")
    return value


def _action(value: AuditAction | str) -> AuditAction:
    try:
        return value if isinstance(value, AuditAction) else AuditAction(value)
    except (TypeError, ValueError):
        raise ValueError("action is not part of the audit vocabulary") from None


def _outcome(value: AuditOutcome | str) -> AuditOutcome:
    try:
        return value if isinstance(value, AuditOutcome) else AuditOutcome(value)
    except (TypeError, ValueError):
        raise ValueError("outcome is invalid") from None


def _actor_kind(value: AuditActorKind | str) -> AuditActorKind:
    try:
        return value if isinstance(value, AuditActorKind) else AuditActorKind(value)
    except (TypeError, ValueError):
        raise ValueError("actor kind is invalid") from None


def _resource_kind(value: AuditResourceKind | str) -> AuditResourceKind:
    try:
        return (
            value if isinstance(value, AuditResourceKind) else AuditResourceKind(value)
        )
    except (TypeError, ValueError):
        raise ValueError("resource kind is invalid") from None


def _reason_code(value: AuditReasonCode | str) -> AuditReasonCode:
    try:
        return value if isinstance(value, AuditReasonCode) else AuditReasonCode(value)
    except (TypeError, ValueError):
        raise ValueError("reason_code is not part of the audit vocabulary") from None


def _utc(value: object) -> datetime:
    if (
        not isinstance(value, datetime)
        or value.tzinfo is None
        or value.utcoffset() is None
    ):
        raise ValueError("occurred_at must be timezone-aware")
    return value.astimezone(timezone.utc)


@dataclass(frozen=True, init=False)
class AuditResource:
    """One domain-separated opaque audit target identity."""

    kind: AuditResourceKind | str
    resource_id: str

    def __init__(self, *_args: object, **_kwargs: object) -> None:
        raise TypeError("AuditResource can only be created by its bounded factory")

    @classmethod
    def _from_derived(
        cls,
        kind: AuditResourceKind,
        resource_id: str,
    ) -> AuditResource:
        expected_prefix = f"ares_{kind.value}_"
        if (
            not resource_id.startswith(expected_prefix)
            or _RESOURCE_DIGEST.fullmatch(resource_id[len(expected_prefix) :]) is None
        ):
            raise ValueError("resource_id must be an opaque audit resource identifier")
        instance = object.__new__(cls)
        object.__setattr__(instance, "kind", kind)
        object.__setattr__(instance, "resource_id", resource_id)
        return instance

    def to_dict(self) -> dict[str, str]:
        return {"kind": self.kind.value, "resource_id": self.resource_id}


def build_audit_resource(
    kind: AuditResourceKind | str,
    *stable_identities: object,
) -> AuditResource:
    """Derive a stable target without retaining the source resource identifier."""

    normalized_kind = _resource_kind(kind)
    normalized_identities = _stable_resource_identities(
        normalized_kind,
        stable_identities,
    )
    digest = sha256()
    digest.update(_RESOURCE_DIGEST_DOMAIN)
    digest.update(normalized_kind.value.encode("ascii"))
    for identity in normalized_identities:
        encoded = identity.encode("ascii")
        digest.update(len(encoded).to_bytes(2, "big"))
        digest.update(encoded)
    return AuditResource._from_derived(
        normalized_kind,
        f"ares_{normalized_kind.value}_{digest.hexdigest()}",
    )


def _stable_resource_identities(
    kind: AuditResourceKind,
    values: tuple[object, ...],
) -> tuple[str, ...]:
    if kind is AuditResourceKind.ACCOUNT:
        if len(values) != 1:
            raise ValueError("account audit resource requires one stable identity")
        return (validate_user_id(values[0]),)
    if kind is AuditResourceKind.RUN:
        patterns = (_RUN_ID,)
    elif kind is AuditResourceKind.ARTIFACT:
        patterns = (_RUN_ID, _ARTIFACT_ID)
    elif kind is AuditResourceKind.REFERENCE:
        patterns = (_REFERENCE_ID,)
    else:
        patterns = (_STORAGE_ID,)
    if len(values) != len(patterns):
        raise ValueError("audit resource has the wrong stable identity coordinates")
    normalized: list[str] = []
    for value, pattern in zip(values, patterns, strict=True):
        if not isinstance(value, str) or pattern.fullmatch(value) is None:
            raise ValueError("audit resource stable identity is invalid")
        normalized.append(value)
    return tuple(normalized)


_ACTION_ACTORS = {
    AuditAction.LOGIN: frozenset({AuditActorKind.USER, AuditActorKind.UNAUTHENTICATED}),
    AuditAction.LOGOUT: frozenset({AuditActorKind.USER}),
    AuditAction.RUN_CREATE: frozenset({AuditActorKind.USER}),
    AuditAction.RUN_START: frozenset({AuditActorKind.USER}),
    AuditAction.RUN_CANCEL: frozenset({AuditActorKind.USER}),
    AuditAction.ARTIFACT_DOWNLOAD: frozenset({AuditActorKind.USER}),
    AuditAction.REFERENCE_REGISTER: frozenset(
        {AuditActorKind.USER, AuditActorKind.LOCAL_OPERATOR}
    ),
    AuditAction.REFERENCE_ENABLE: frozenset(
        {AuditActorKind.USER, AuditActorKind.LOCAL_OPERATOR}
    ),
    AuditAction.REFERENCE_DISABLE: frozenset(
        {AuditActorKind.USER, AuditActorKind.LOCAL_OPERATOR}
    ),
    AuditAction.STORAGE_REGISTER: frozenset(
        {AuditActorKind.USER, AuditActorKind.LOCAL_OPERATOR}
    ),
    AuditAction.STORAGE_ARCHIVE: frozenset(
        {AuditActorKind.USER, AuditActorKind.LOCAL_OPERATOR}
    ),
    AuditAction.ACCOUNT_CREATE: frozenset(
        {AuditActorKind.USER, AuditActorKind.LOCAL_OPERATOR}
    ),
    AuditAction.ACCOUNT_ENABLE: frozenset(
        {AuditActorKind.USER, AuditActorKind.LOCAL_OPERATOR}
    ),
    AuditAction.ACCOUNT_DISABLE: frozenset(
        {AuditActorKind.USER, AuditActorKind.LOCAL_OPERATOR}
    ),
    AuditAction.ACCOUNT_PASSWORD_RESET: frozenset(
        {AuditActorKind.USER, AuditActorKind.LOCAL_OPERATOR}
    ),
    AuditAction.ACCOUNT_SESSIONS_REVOKE: frozenset(
        {AuditActorKind.USER, AuditActorKind.LOCAL_OPERATOR}
    ),
    AuditAction.RUN_FAIL: frozenset(
        {AuditActorKind.USER, AuditActorKind.LOCAL_OPERATOR}
    ),
    AuditAction.RUN_REQUEUE: frozenset(
        {AuditActorKind.USER, AuditActorKind.LOCAL_OPERATOR}
    ),
}

_ACTION_RESOURCE_KINDS = {
    AuditAction.RUN_CREATE: AuditResourceKind.RUN,
    AuditAction.RUN_START: AuditResourceKind.RUN,
    AuditAction.RUN_CANCEL: AuditResourceKind.RUN,
    AuditAction.ARTIFACT_DOWNLOAD: AuditResourceKind.ARTIFACT,
    AuditAction.REFERENCE_REGISTER: AuditResourceKind.REFERENCE,
    AuditAction.REFERENCE_ENABLE: AuditResourceKind.REFERENCE,
    AuditAction.REFERENCE_DISABLE: AuditResourceKind.REFERENCE,
    AuditAction.STORAGE_REGISTER: AuditResourceKind.STORAGE,
    AuditAction.STORAGE_ARCHIVE: AuditResourceKind.STORAGE,
    AuditAction.ACCOUNT_CREATE: AuditResourceKind.ACCOUNT,
    AuditAction.ACCOUNT_ENABLE: AuditResourceKind.ACCOUNT,
    AuditAction.ACCOUNT_DISABLE: AuditResourceKind.ACCOUNT,
    AuditAction.ACCOUNT_PASSWORD_RESET: AuditResourceKind.ACCOUNT,
    AuditAction.ACCOUNT_SESSIONS_REVOKE: AuditResourceKind.ACCOUNT,
    AuditAction.RUN_FAIL: AuditResourceKind.RUN,
    AuditAction.RUN_REQUEUE: AuditResourceKind.RUN,
}

_LOGIN_FAILURE_REASONS = frozenset(
    {
        AuditReasonCode.INVALID_CREDENTIALS,
        AuditReasonCode.LOGIN_RATE_LIMITED,
        AuditReasonCode.SETUP_REQUIRED,
        AuditReasonCode.INTERNAL_FAILURE,
    }
)


@dataclass(frozen=True)
class SecurityAuditEvent:
    """A closed audit record with no arbitrary payload or raw request fields."""

    event_id: str
    occurred_at: datetime
    action: AuditAction | str
    outcome: AuditOutcome | str
    actor_kind: AuditActorKind | str
    actor_user_id: str | None = None
    resource: AuditResource | None = None
    reason_code: AuditReasonCode | str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "event_id", _audit_event_id(self.event_id))
        object.__setattr__(self, "occurred_at", _utc(self.occurred_at))
        action = _action(self.action)
        object.__setattr__(self, "action", action)
        outcome = _outcome(self.outcome)
        actor_kind = _actor_kind(self.actor_kind)
        object.__setattr__(self, "outcome", outcome)
        object.__setattr__(self, "actor_kind", actor_kind)
        if self.actor_user_id is not None:
            object.__setattr__(
                self,
                "actor_user_id",
                validate_user_id(self.actor_user_id),
            )
        if (actor_kind is AuditActorKind.USER) != (self.actor_user_id is not None):
            raise ValueError("only a user actor may carry actor_user_id")
        if self.resource is not None and not isinstance(self.resource, AuditResource):
            raise ValueError("resource must be an AuditResource or None")
        expected_resource_kind = _ACTION_RESOURCE_KINDS.get(action)
        if expected_resource_kind is None:
            if self.resource is not None:
                raise ValueError("this audit action must not carry a resource")
        elif self.resource is None or self.resource.kind is not expected_resource_kind:
            raise ValueError("audit resource kind does not match the action")
        if actor_kind not in _ACTION_ACTORS[action]:
            raise ValueError("audit actor kind is not allowed for the action")
        if self.reason_code is not None:
            object.__setattr__(
                self,
                "reason_code",
                _reason_code(self.reason_code),
            )
        if (outcome is AuditOutcome.SUCCEEDED) != (self.reason_code is None):
            raise ValueError(
                "successful events omit reason_code; denied or failed events require it"
            )
        if action is AuditAction.LOGIN:
            if outcome is AuditOutcome.SUCCEEDED:
                if actor_kind is not AuditActorKind.USER:
                    raise ValueError("successful login requires the authenticated user")
            elif (
                actor_kind is not AuditActorKind.UNAUTHENTICATED
                or self.reason_code not in _LOGIN_FAILURE_REASONS
            ):
                raise ValueError(
                    "failed login must use a non-identifying actor and safe reason"
                )

    def to_dict(self) -> dict[str, object]:
        """Return the complete and intentionally narrow persisted representation."""

        return {
            "event_id": self.event_id,
            "occurred_at": self.occurred_at,
            "action": self.action.value,
            "outcome": self.outcome.value,
            "actor_kind": self.actor_kind.value,
            "actor_user_id": self.actor_user_id,
            "resource": self.resource.to_dict() if self.resource else None,
            "reason_code": self.reason_code.value if self.reason_code else None,
        }
