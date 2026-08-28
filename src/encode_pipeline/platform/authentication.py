"""Workflow-neutral local account and browser-session contracts."""

from __future__ import annotations

import base64
import binascii
from dataclasses import dataclass, field, replace
from datetime import datetime, timedelta, timezone
from enum import Enum
import re

from encode_pipeline.platform.notifications import normalize_notification_email


USERNAME_MIN_LENGTH = 3
USERNAME_MAX_LENGTH = 64
PASSWORD_HASH_MAX_LENGTH = 256

ARGON2ID_VERSION = 19
ARGON2ID_MIN_TIME_COST = 2
ARGON2ID_MAX_TIME_COST = 10
ARGON2ID_MIN_MEMORY_COST_KIB = 19 * 1024
ARGON2ID_MAX_MEMORY_COST_KIB = 256 * 1024
ARGON2ID_MIN_PARALLELISM = 1
ARGON2ID_MAX_PARALLELISM = 16
ARGON2ID_MIN_HASH_LENGTH = 32
ARGON2ID_MAX_HASH_LENGTH = 64
ARGON2ID_MIN_SALT_LENGTH = 16
ARGON2ID_MAX_SALT_LENGTH = 64

MIN_SESSION_LIFETIME = timedelta(minutes=5)
MAX_SESSION_LIFETIME = timedelta(days=7)

_USER_ID = re.compile(r"^usr_[0-9a-f]{32}$")
_USERNAME = re.compile(r"^[a-z][a-z0-9._-]*$")
_SHA256_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_ARGON2ID_HASH = re.compile(
    r"^\$argon2id\$v=19\$"
    r"m=([1-9][0-9]{0,6}),t=([1-9][0-9]{0,2}),p=([1-9][0-9]{0,2})\$"
    r"([A-Za-z0-9+/]{1,86})\$([A-Za-z0-9+/]{1,86})$"
)


class UserRole(str, Enum):
    """The complete role set for the trusted-laboratory LAN product."""

    ADMINISTRATOR = "administrator"
    MEMBER = "member"


class UserStatus(str, Enum):
    """Whether an account may establish or use authenticated sessions."""

    ENABLED = "enabled"
    DISABLED = "disabled"


class SessionRevocationReason(str, Enum):
    """Stable reasons for invalidating an opaque browser session."""

    LOGOUT = "logout"
    ALL_SESSIONS = "all_sessions"
    ACCOUNT_DISABLED = "account_disabled"
    PASSWORD_RESET = "password_reset"


def validate_user_id(value: object) -> str:
    """Return one path-safe opaque user identifier or reject it."""

    if not isinstance(value, str) or _USER_ID.fullmatch(value) is None:
        raise ValueError("user_id must be an opaque user identifier")
    return value


def normalize_username(value: object) -> str:
    """Return the canonical bounded ASCII login name.

    Usernames are deliberately narrower than human display names. The product does
    not currently have a display-name or email identity contract, and accepting a
    broad Unicode identifier here would create avoidable confusable-name behavior.
    """

    if not isinstance(value, str):
        raise ValueError("username must be a string")
    if len(value) > USERNAME_MAX_LENGTH + 2:
        raise ValueError("username is too long")
    if not value.isascii():
        raise ValueError("username must use ASCII characters")
    stripped = value.strip(" ")
    normalized = stripped.casefold()
    if not USERNAME_MIN_LENGTH <= len(normalized) <= USERNAME_MAX_LENGTH:
        raise ValueError(
            f"username must contain {USERNAME_MIN_LENGTH} to "
            f"{USERNAME_MAX_LENGTH} characters"
        )
    if _USERNAME.fullmatch(normalized) is None:
        raise ValueError("username contains unsupported characters")
    return normalized


def validate_sha256_digest(value: object, name: str) -> str:
    """Return a lowercase SHA-256 digest used in a durable security record."""

    if not isinstance(value, str) or _SHA256_DIGEST.fullmatch(value) is None:
        raise ValueError(f"{name} must be a lowercase SHA-256 digest")
    return value


def validate_password_hash(value: object) -> str:
    """Validate the complete bounded Argon2id PHC persistence contract."""

    if not isinstance(value, str) or len(value) > PASSWORD_HASH_MAX_LENGTH:
        raise ValueError("password_hash must be a bounded Argon2id PHC string")
    match = _ARGON2ID_HASH.fullmatch(value)
    if match is None:
        raise ValueError("password_hash must be a bounded Argon2id PHC string")
    memory_cost, time_cost, parallelism = (int(item) for item in match.groups()[:3])
    if not ARGON2ID_MIN_TIME_COST <= time_cost <= ARGON2ID_MAX_TIME_COST:
        raise ValueError("password_hash time cost is outside the accepted boundary")
    if not (
        ARGON2ID_MIN_MEMORY_COST_KIB <= memory_cost <= ARGON2ID_MAX_MEMORY_COST_KIB
    ):
        raise ValueError("password_hash memory cost is outside the accepted boundary")
    if not (ARGON2ID_MIN_PARALLELISM <= parallelism <= ARGON2ID_MAX_PARALLELISM):
        raise ValueError("password_hash parallelism is outside the accepted boundary")
    salt = _decode_phc_component(match.group(4), "salt")
    digest = _decode_phc_component(match.group(5), "digest")
    if not ARGON2ID_MIN_SALT_LENGTH <= len(salt) <= ARGON2ID_MAX_SALT_LENGTH:
        raise ValueError("password_hash salt length is outside the accepted boundary")
    if not ARGON2ID_MIN_HASH_LENGTH <= len(digest) <= ARGON2ID_MAX_HASH_LENGTH:
        raise ValueError("password_hash digest length is outside the accepted boundary")
    return value


def _decode_phc_component(value: str, name: str) -> bytes:
    padding = "=" * (-len(value) % 4)
    try:
        decoded = base64.b64decode(value + padding, validate=True)
    except (binascii.Error, ValueError):
        raise ValueError(f"password_hash {name} is not canonical base64") from None
    canonical = base64.b64encode(decoded).decode("ascii").rstrip("=")
    if canonical != value:
        raise ValueError(f"password_hash {name} is not canonical base64")
    return decoded


def _utc(value: object, name: str) -> datetime:
    if (
        not isinstance(value, datetime)
        or value.tzinfo is None
        or value.utcoffset() is None
    ):
        raise ValueError(f"{name} must be timezone-aware")
    return value.astimezone(timezone.utc)


def _role(value: UserRole | str) -> UserRole:
    try:
        return value if isinstance(value, UserRole) else UserRole(value)
    except (TypeError, ValueError):
        raise ValueError("role must be administrator or member") from None


def _status(value: UserStatus | str) -> UserStatus:
    try:
        return value if isinstance(value, UserStatus) else UserStatus(value)
    except (TypeError, ValueError):
        raise ValueError("status must be enabled or disabled") from None


def _revocation_reason(
    value: SessionRevocationReason | str | None,
) -> SessionRevocationReason | None:
    if value is None or isinstance(value, SessionRevocationReason):
        return value
    try:
        return SessionRevocationReason(value)
    except (TypeError, ValueError):
        raise ValueError("revocation_reason is invalid") from None


@dataclass(frozen=True)
class UserAccount:
    """Immutable canonical state for one local laboratory account."""

    user_id: str
    username: str
    role: UserRole | str
    status: UserStatus | str
    password_hash: str = field(repr=False)
    created_at: datetime
    updated_at: datetime
    password_changed_at: datetime
    notification_email: str | None = field(default=None, repr=False)
    terminal_email_enabled: bool = True

    def __post_init__(self) -> None:
        object.__setattr__(self, "user_id", validate_user_id(self.user_id))
        object.__setattr__(self, "username", normalize_username(self.username))
        object.__setattr__(self, "role", _role(self.role))
        object.__setattr__(self, "status", _status(self.status))
        object.__setattr__(
            self,
            "password_hash",
            validate_password_hash(self.password_hash),
        )
        created_at = _utc(self.created_at, "created_at")
        updated_at = _utc(self.updated_at, "updated_at")
        password_changed_at = _utc(
            self.password_changed_at,
            "password_changed_at",
        )
        if updated_at < created_at:
            raise ValueError("updated_at must not precede created_at")
        if not created_at <= password_changed_at <= updated_at:
            raise ValueError("password_changed_at must be within the account lifecycle")
        object.__setattr__(self, "created_at", created_at)
        object.__setattr__(self, "updated_at", updated_at)
        object.__setattr__(self, "password_changed_at", password_changed_at)
        object.__setattr__(
            self,
            "notification_email",
            (
                normalize_notification_email(self.notification_email)
                if self.notification_email is not None
                else None
            ),
        )
        if not isinstance(self.terminal_email_enabled, bool):
            raise ValueError("terminal_email_enabled must be boolean")

    @property
    def enabled(self) -> bool:
        return self.status is UserStatus.ENABLED

    def disable(self, changed_at: datetime) -> UserAccount:
        """Return an idempotently disabled account snapshot."""

        if self.status is UserStatus.DISABLED:
            return self
        return replace(
            self,
            status=UserStatus.DISABLED,
            updated_at=self._change_time(changed_at),
        )

    def enable(self, changed_at: datetime) -> UserAccount:
        """Return an idempotently enabled account snapshot."""

        if self.status is UserStatus.ENABLED:
            return self
        return replace(
            self,
            status=UserStatus.ENABLED,
            updated_at=self._change_time(changed_at),
        )

    def change_password(self, password_hash: str, changed_at: datetime) -> UserAccount:
        """Return the next password state; session revocation is repository-owned."""

        normalized_time = self._change_time(changed_at)
        return replace(
            self,
            password_hash=validate_password_hash(password_hash),
            password_changed_at=normalized_time,
            updated_at=normalized_time,
        )

    def change_notification_email(
        self,
        notification_email: object | None,
        changed_at: datetime,
    ) -> UserAccount:
        """Return a private non-identity contact update."""

        normalized = (
            normalize_notification_email(notification_email)
            if notification_email is not None
            else None
        )
        if normalized == self.notification_email:
            return self
        return replace(
            self,
            notification_email=normalized,
            updated_at=self._change_time(changed_at),
        )

    def change_terminal_email_enabled(
        self,
        enabled: object,
        changed_at: datetime,
    ) -> UserAccount:
        """Return the member-controlled terminal-email opt-out update."""

        if not isinstance(enabled, bool):
            raise ValueError("terminal_email_enabled must be boolean")
        if enabled is self.terminal_email_enabled:
            return self
        return replace(
            self,
            terminal_email_enabled=enabled,
            updated_at=self._change_time(changed_at),
        )

    def to_public_summary(self) -> dict[str, object]:
        """Return the exact path-safe account fields suitable for operator output."""

        return {
            "user_id": self.user_id,
            "username": self.username,
            "role": self.role.value,
            "status": self.status.value,
            "created_at": self.created_at,
            "updated_at": self.updated_at,
            "password_changed_at": self.password_changed_at,
        }

    def _change_time(self, value: datetime) -> datetime:
        changed_at = _utc(value, "changed_at")
        if changed_at < self.updated_at:
            raise ValueError("changed_at must not precede updated_at")
        return changed_at


@dataclass(frozen=True)
class TerminalEmailPreference:
    """Safe self-service projection that never exposes the private address."""

    terminal_email_enabled: bool
    address_configured: bool

    def __post_init__(self) -> None:
        if not isinstance(self.terminal_email_enabled, bool):
            raise ValueError("terminal_email_enabled must be boolean")
        if not isinstance(self.address_configured, bool):
            raise ValueError("address_configured must be boolean")

    @classmethod
    def from_account(cls, account: UserAccount) -> TerminalEmailPreference:
        if not isinstance(account, UserAccount):
            raise ValueError("account must be a UserAccount")
        return cls(
            terminal_email_enabled=account.terminal_email_enabled,
            address_configured=account.notification_email is not None,
        )


@dataclass(frozen=True)
class AuthenticatedPrincipal:
    """The only account identity services need for role checks and audit."""

    user_id: str
    username: str
    role: UserRole | str

    def __post_init__(self) -> None:
        object.__setattr__(self, "user_id", validate_user_id(self.user_id))
        object.__setattr__(self, "username", normalize_username(self.username))
        object.__setattr__(self, "role", _role(self.role))

    @classmethod
    def from_account(cls, account: UserAccount) -> AuthenticatedPrincipal:
        if not isinstance(account, UserAccount) or not account.enabled:
            raise ValueError("an authenticated principal requires an enabled account")
        return cls(
            user_id=account.user_id,
            username=account.username,
            role=account.role,
        )

    @property
    def is_administrator(self) -> bool:
        return self.role is UserRole.ADMINISTRATOR


@dataclass(frozen=True)
class SessionRecord:
    """Server-side state for one random opaque browser-session credential.

    Only digests of the browser-held session and CSRF values belong in this
    record. Raw values must stay at the request/response boundary.
    """

    session_digest: str = field(repr=False)
    csrf_digest: str = field(repr=False)
    user_id: str
    created_at: datetime
    expires_at: datetime
    revoked_at: datetime | None = None
    revocation_reason: SessionRevocationReason | str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "session_digest",
            validate_sha256_digest(self.session_digest, "session_digest"),
        )
        object.__setattr__(
            self,
            "csrf_digest",
            validate_sha256_digest(self.csrf_digest, "csrf_digest"),
        )
        object.__setattr__(self, "user_id", validate_user_id(self.user_id))
        created_at = _utc(self.created_at, "created_at")
        expires_at = _utc(self.expires_at, "expires_at")
        lifetime = expires_at - created_at
        if not MIN_SESSION_LIFETIME <= lifetime <= MAX_SESSION_LIFETIME:
            raise ValueError(
                "expires_at must be between five minutes and seven days after created_at"
            )
        if lifetime.microseconds != 0:
            raise ValueError("session lifetime must use whole seconds")
        revoked_at = (
            _utc(self.revoked_at, "revoked_at") if self.revoked_at is not None else None
        )
        reason = _revocation_reason(self.revocation_reason)
        if (revoked_at is None) != (reason is None):
            raise ValueError(
                "revoked_at and revocation_reason must either both be set or both be absent"
            )
        if revoked_at is not None and revoked_at < created_at:
            raise ValueError("revoked_at must not precede created_at")
        object.__setattr__(self, "created_at", created_at)
        object.__setattr__(self, "expires_at", expires_at)
        object.__setattr__(self, "revoked_at", revoked_at)
        object.__setattr__(self, "revocation_reason", reason)

    def active_at(self, checked_at: datetime) -> bool:
        """Return whether the server-side record is usable at ``checked_at``."""

        now = _utc(checked_at, "checked_at")
        return self.revoked_at is None and self.created_at <= now < self.expires_at

    def revoke(
        self,
        revoked_at: datetime,
        reason: SessionRevocationReason | str,
    ) -> SessionRecord:
        """Return an idempotently revoked session snapshot."""

        if self.revoked_at is not None:
            return self
        normalized_time = _utc(revoked_at, "revoked_at")
        if normalized_time < self.created_at:
            raise ValueError("revoked_at must not precede created_at")
        return replace(
            self,
            revoked_at=normalized_time,
            revocation_reason=_revocation_reason(reason),
        )


def authenticated_principal_for_session(
    session: SessionRecord,
    account: UserAccount | None,
    checked_at: datetime,
) -> AuthenticatedPrincipal | None:
    """Resolve an active session while enforcing immediate account disablement."""

    if not isinstance(session, SessionRecord):
        raise ValueError("session must be a SessionRecord")
    if (
        account is None
        or not isinstance(account, UserAccount)
        or account.user_id != session.user_id
        or not account.enabled
        or not session.active_at(checked_at)
    ):
        return None
    return AuthenticatedPrincipal.from_account(account)
