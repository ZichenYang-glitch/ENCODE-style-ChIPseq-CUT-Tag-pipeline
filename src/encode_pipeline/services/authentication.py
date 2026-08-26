"""Security primitives for local passwords and opaque browser sessions."""

from __future__ import annotations

from collections import OrderedDict, deque
from collections.abc import Callable
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from hashlib import sha256
import hmac
import math
import re
from secrets import token_hex, token_urlsafe
from threading import Lock
from time import monotonic

from argon2 import PasswordHasher, extract_parameters
from argon2.exceptions import InvalidHashError, VerificationError, VerifyMismatchError
from argon2.low_level import Type

from encode_pipeline.platform.authentication import (
    ARGON2ID_MAX_HASH_LENGTH,
    ARGON2ID_MAX_MEMORY_COST_KIB,
    ARGON2ID_MAX_PARALLELISM,
    ARGON2ID_MAX_SALT_LENGTH,
    ARGON2ID_MAX_TIME_COST,
    ARGON2ID_MIN_HASH_LENGTH,
    ARGON2ID_MIN_MEMORY_COST_KIB,
    ARGON2ID_MIN_PARALLELISM,
    ARGON2ID_MIN_SALT_LENGTH,
    ARGON2ID_MIN_TIME_COST,
    ARGON2ID_VERSION,
    MAX_SESSION_LIFETIME,
    MIN_SESSION_LIFETIME,
    AuthenticatedPrincipal,
    SessionRecord,
    UserAccount,
    normalize_username,
    validate_password_hash,
    validate_sha256_digest,
)


PASSWORD_HASH_SCHEME = "argon2id-v1"
PASSWORD_MIN_CHARACTERS = 15
PASSWORD_MAX_CHARACTERS = 128
PASSWORD_MAX_UTF8_BYTES = 512

SESSION_TOKEN_BYTES = 32
CSRF_TOKEN_BYTES = 32
SESSION_TOKEN_LENGTH = 43
CSRF_TOKEN_LENGTH = 43
DEFAULT_SESSION_LIFETIME = timedelta(hours=8)
MAX_LEGACY_PASSWORD_PROFILES = 2

DEVELOPMENT_SESSION_COOKIE_NAME = "helixweave_session"
DEVELOPMENT_CSRF_COOKIE_NAME = "helixweave_csrf"
SECURE_SESSION_COOKIE_NAME = "__Host-helixweave_session"
SECURE_CSRF_COOKIE_NAME = "__Host-helixweave_csrf"
CSRF_HEADER_NAME = "X-CSRF-Token"

INVALID_LOGIN_CODE = "INVALID_CREDENTIALS"
INVALID_LOGIN_MESSAGE = "The username or password is invalid."

_OPAQUE_TOKEN = re.compile(r"^[A-Za-z0-9_-]{43}$")
_HTTP_METHOD = re.compile(r"^[A-Za-z]{1,16}$")
_SESSION_DIGEST_DOMAIN = b"helixweave-session-token-v1\0"
_CSRF_DIGEST_DOMAIN = b"helixweave-csrf-token-v1\0"
_RATE_LIMIT_DOMAIN = b"helixweave-login-rate-limit-v1\0"
_SAFE_METHODS = frozenset({"GET", "HEAD", "OPTIONS"})


@dataclass(frozen=True)
class Argon2idParameters:
    """Versioned password-hash cost parameters for new and upgraded hashes."""

    time_cost: int = 3
    memory_cost_kib: int = 64 * 1024
    parallelism: int = 4
    hash_length: int = 32
    salt_length: int = 16

    def __post_init__(self) -> None:
        for name in (
            "time_cost",
            "memory_cost_kib",
            "parallelism",
            "hash_length",
            "salt_length",
        ):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool):
                raise ValueError(f"{name} must be an integer")
        if self.time_cost < ARGON2ID_MIN_TIME_COST:
            raise ValueError("time_cost must be at least 2")
        if self.time_cost > ARGON2ID_MAX_TIME_COST:
            raise ValueError("time_cost must not exceed 10")
        if self.memory_cost_kib < ARGON2ID_MIN_MEMORY_COST_KIB:
            raise ValueError("memory_cost_kib must be at least 19456")
        if self.memory_cost_kib > ARGON2ID_MAX_MEMORY_COST_KIB:
            raise ValueError("memory_cost_kib must not exceed 262144")
        if self.parallelism < ARGON2ID_MIN_PARALLELISM:
            raise ValueError("parallelism must be positive")
        if self.parallelism > ARGON2ID_MAX_PARALLELISM:
            raise ValueError("parallelism must not exceed 16")
        if self.hash_length < ARGON2ID_MIN_HASH_LENGTH:
            raise ValueError("hash_length must be at least 32")
        if self.hash_length > ARGON2ID_MAX_HASH_LENGTH:
            raise ValueError("hash_length must not exceed 64")
        if self.salt_length < ARGON2ID_MIN_SALT_LENGTH:
            raise ValueError("salt_length must be at least 16")
        if self.salt_length > ARGON2ID_MAX_SALT_LENGTH:
            raise ValueError("salt_length must not exceed 64")


@dataclass(frozen=True)
class PasswordVerification:
    """Internal verification result; an optional replacement hash is never repr'd."""

    verified: bool
    replacement_hash: str | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.verified, bool):
            raise ValueError("verified must be boolean")
        if not self.verified and self.replacement_hash is not None:
            raise ValueError("failed verification cannot carry a replacement hash")


class PasswordManager:
    """Bounded Argon2id hashing, verification, and parameter upgrades."""

    def __init__(
        self,
        parameters: Argon2idParameters | None = None,
        *,
        legacy_parameters: tuple[Argon2idParameters, ...] = (),
    ) -> None:
        if parameters is not None and not isinstance(parameters, Argon2idParameters):
            raise ValueError("parameters must be an Argon2id parameter profile")
        self.parameters = parameters or Argon2idParameters()
        if not isinstance(legacy_parameters, tuple) or any(
            not isinstance(profile, Argon2idParameters) for profile in legacy_parameters
        ):
            raise ValueError("legacy_parameters must be Argon2id parameter profiles")
        if len(legacy_parameters) > MAX_LEGACY_PASSWORD_PROFILES:
            raise ValueError("at most two legacy password profiles may be active")
        profiles = (self.parameters, *legacy_parameters)
        if len(set(profiles)) != len(profiles):
            raise ValueError("password parameter profiles must be unique")
        self._hashers = {profile: _password_hasher(profile) for profile in profiles}
        self._hasher = self._hashers[self.parameters]
        # Each approved profile gets a process-local dummy. Every bounded login
        # attempt verifies exactly one hash for every profile, substituting the
        # stored hash only in its matching slot. Account absence, disablement, and
        # a planned cost-rotation window therefore retain the same KDF work shape.
        self._dummy_hashes = {
            profile: hasher.hash(token_urlsafe(32))
            for profile, hasher in self._hashers.items()
        }

    def hash_password(self, password: object) -> str:
        """Hash a new password after enforcing the fixed resource boundary."""

        normalized = _password(password, require_minimum=True)
        return self._hasher.hash(normalized)

    def verify_password(
        self,
        password: object,
        encoded_hash: str | None,
        *,
        allow_rehash: bool = True,
    ) -> PasswordVerification:
        """Verify one candidate and optionally return an upgraded Argon2id hash."""

        if not isinstance(allow_rehash, bool):
            raise ValueError("allow_rehash must be boolean")
        try:
            normalized = _password(password, require_minimum=True)
        except ValueError:
            return PasswordVerification(verified=False)

        stored = _stored_password_hash(encoded_hash)
        stored_hash: str | None = None
        stored_profile: Argon2idParameters | None = None
        if stored is not None and stored[1] in self._hashers:
            stored_hash, stored_profile = stored

        stored_verified = False
        for profile, hasher in self._hashers.items():
            candidate_hash = (
                stored_hash
                if stored_profile == profile and stored_hash is not None
                else self._dummy_hashes[profile]
            )
            try:
                verified = bool(hasher.verify(candidate_hash, normalized))
            except VerifyMismatchError:
                verified = False
            except (InvalidHashError, VerificationError):
                verified = False
            if stored_profile == profile:
                stored_verified = verified

        if not stored_verified or stored_hash is None or stored_profile is None:
            return PasswordVerification(verified=False)
        replacement = (
            self._hasher.hash(normalized)
            if allow_rehash
            and (
                stored_profile != self.parameters
                or self._hasher.check_needs_rehash(stored_hash)
            )
            else None
        )
        return PasswordVerification(
            verified=True,
            replacement_hash=replacement,
        )


@dataclass(frozen=True)
class LoginAuthentication:
    """Internal result whose rejected form is identical for every account state."""

    principal: AuthenticatedPrincipal | None = None
    replacement_hash: str | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        if self.principal is not None and not isinstance(
            self.principal,
            AuthenticatedPrincipal,
        ):
            raise ValueError("principal must be an AuthenticatedPrincipal or None")
        if self.principal is None and self.replacement_hash is not None:
            raise ValueError("rejected login cannot carry a replacement hash")

    @property
    def authenticated(self) -> bool:
        return self.principal is not None


def authenticate_login_candidate(
    account: UserAccount | None,
    password: object,
    password_manager: PasswordManager,
) -> LoginAuthentication:
    """Apply constant failure semantics to absent, disabled, and mismatched users."""

    if account is not None and not isinstance(account, UserAccount):
        raise ValueError("account must be a UserAccount or None")
    if not isinstance(password_manager, PasswordManager):
        raise ValueError("password_manager must be a PasswordManager")
    verification = password_manager.verify_password(
        password,
        account.password_hash if account is not None else None,
        allow_rehash=bool(account is not None and account.enabled),
    )
    if account is None or not account.enabled or not verification.verified:
        return LoginAuthentication()
    return LoginAuthentication(
        principal=AuthenticatedPrincipal.from_account(account),
        replacement_hash=verification.replacement_hash,
    )


def new_user_id() -> str:
    """Return a random stable path-safe identity for a local user."""

    return f"usr_{token_hex(16)}"


@dataclass(frozen=True, repr=False, init=False)
class SessionSecrets:
    """Short-lived raw browser values that must never enter persistence or logs."""

    session_token: str = field(repr=False)
    csrf_token: str = field(repr=False)

    def __init__(self, *_args: object, **_kwargs: object) -> None:
        raise TypeError("SessionSecrets can only be created by the secure generator")

    @classmethod
    def _from_generated(cls, session_token: str, csrf_token: str) -> SessionSecrets:
        session_token = _opaque_token(session_token, "session token")
        csrf_token = _opaque_token(csrf_token, "CSRF token")
        if hmac.compare_digest(session_token, csrf_token):
            raise ValueError("session and CSRF tokens must be independent")
        instance = object.__new__(cls)
        object.__setattr__(instance, "session_token", session_token)
        object.__setattr__(instance, "csrf_token", csrf_token)
        return instance


def new_session_secrets() -> SessionSecrets:
    """Generate independent 256-bit session and synchronizer-CSRF values."""

    session_token = token_urlsafe(SESSION_TOKEN_BYTES)
    for _attempt in range(4):
        csrf_token = token_urlsafe(CSRF_TOKEN_BYTES)
        if not hmac.compare_digest(session_token, csrf_token):
            return SessionSecrets._from_generated(session_token, csrf_token)
    raise RuntimeError("failed to generate independent session values")


def digest_session_token(token: object) -> str:
    """Return the only representation of a session token allowed in SQLite."""

    normalized = _opaque_token(token, "session token")
    return _digest_token(_SESSION_DIGEST_DOMAIN, normalized)


def digest_csrf_token(token: object) -> str:
    """Return the server-side digest of a synchronizer-CSRF value."""

    normalized = _opaque_token(token, "CSRF token")
    return _digest_token(_CSRF_DIGEST_DOMAIN, normalized)


def csrf_token_matches(token: object, expected_digest: str) -> bool:
    """Compare a browser-provided CSRF token with the durable session digest."""

    expected = validate_sha256_digest(expected_digest, "expected_digest")
    valid = isinstance(token, str) and _OPAQUE_TOKEN.fullmatch(token) is not None
    candidate = token if valid else "A" * CSRF_TOKEN_LENGTH
    actual = _digest_token(_CSRF_DIGEST_DOMAIN, candidate)
    matches = hmac.compare_digest(actual, expected)
    return bool(valid and matches)


def csrf_request_is_valid(
    cookie_token: object,
    header_token: object,
    expected_digest: str,
) -> bool:
    """Require the readable cookie, request header, and server digest to agree."""

    cookie_valid = (
        isinstance(cookie_token, str)
        and _OPAQUE_TOKEN.fullmatch(cookie_token) is not None
    )
    header_valid = (
        isinstance(header_token, str)
        and _OPAQUE_TOKEN.fullmatch(header_token) is not None
    )
    cookie_candidate = cookie_token if cookie_valid else "A" * CSRF_TOKEN_LENGTH
    header_candidate = header_token if header_valid else "B" * CSRF_TOKEN_LENGTH
    browser_values_match = hmac.compare_digest(
        cookie_candidate,
        header_candidate,
    )
    digest_matches = csrf_token_matches(header_candidate, expected_digest)
    return bool(
        cookie_valid and header_valid and browser_values_match and digest_matches
    )


def new_session_record(
    *,
    user_id: str,
    secrets: SessionSecrets,
    created_at: datetime,
    lifetime: timedelta = DEFAULT_SESSION_LIFETIME,
) -> SessionRecord:
    """Create server-side state for newly generated browser credentials."""

    if not isinstance(secrets, SessionSecrets):
        raise ValueError("secrets must be SessionSecrets")
    normalized_lifetime = _session_lifetime(lifetime)
    normalized_created_at = _utc_datetime(created_at, "created_at")
    return SessionRecord(
        session_digest=digest_session_token(secrets.session_token),
        csrf_digest=digest_csrf_token(secrets.csrf_token),
        user_id=user_id,
        created_at=normalized_created_at,
        expires_at=normalized_created_at + normalized_lifetime,
    )


@dataclass(frozen=True)
class CookieDirective:
    """One secret-free cookie policy projection for the HTTP adapter."""

    name: str
    http_only: bool
    secure: bool
    same_site: str
    max_age_seconds: int
    path: str = "/"

    def __post_init__(self) -> None:
        allowed_names = {
            DEVELOPMENT_SESSION_COOKIE_NAME,
            DEVELOPMENT_CSRF_COOKIE_NAME,
            SECURE_SESSION_COOKIE_NAME,
            SECURE_CSRF_COOKIE_NAME,
        }
        if self.name not in allowed_names:
            raise ValueError("cookie name is not part of the auth contract")
        if not isinstance(self.http_only, bool) or not isinstance(self.secure, bool):
            raise ValueError("cookie flags must be boolean")
        session_cookie_names = {
            DEVELOPMENT_SESSION_COOKIE_NAME,
            SECURE_SESSION_COOKIE_NAME,
        }
        if self.http_only != (self.name in session_cookie_names):
            raise ValueError(
                "session cookies must be HttpOnly and CSRF cookies readable"
            )
        if self.same_site not in {"lax", "strict"}:
            raise ValueError("same_site must be lax or strict")
        if self.path != "/":
            raise ValueError("authentication cookies must be host-wide")
        if self.secure != self.name.startswith("__Host-"):
            raise ValueError("secure cookies must use the __Host- name variant")
        if (
            not isinstance(self.max_age_seconds, int)
            or isinstance(self.max_age_seconds, bool)
            or self.max_age_seconds <= 0
        ):
            raise ValueError("max_age_seconds must be positive")

    def to_response_kwargs(self) -> dict[str, object]:
        """Return Starlette-compatible flags without a cookie value or domain."""

        return {
            "key": self.name,
            "max_age": self.max_age_seconds,
            "path": self.path,
            "secure": self.secure,
            "httponly": self.http_only,
            "samesite": self.same_site,
        }


@dataclass(frozen=True)
class BrowserSessionCookiePolicy:
    """Explicit browser-cookie lifecycle with a readable CSRF companion cookie."""

    secure: bool
    same_site: str = "lax"
    lifetime: timedelta = DEFAULT_SESSION_LIFETIME

    def __post_init__(self) -> None:
        if not isinstance(self.secure, bool):
            raise ValueError("secure must be boolean")
        if self.same_site not in {"lax", "strict"}:
            raise ValueError("same_site must be lax or strict")
        _session_lifetime(self.lifetime)

    @property
    def session_cookie(self) -> CookieDirective:
        return CookieDirective(
            name=(
                SECURE_SESSION_COOKIE_NAME
                if self.secure
                else DEVELOPMENT_SESSION_COOKIE_NAME
            ),
            http_only=True,
            secure=self.secure,
            same_site=self.same_site,
            max_age_seconds=int(self.lifetime.total_seconds()),
        )

    @property
    def csrf_cookie(self) -> CookieDirective:
        # Browser code must echo this synchronizer value in CSRF_HEADER_NAME. It
        # has no authentication authority without the independent HttpOnly cookie.
        return CookieDirective(
            name=(
                SECURE_CSRF_COOKIE_NAME if self.secure else DEVELOPMENT_CSRF_COOKIE_NAME
            ),
            http_only=False,
            secure=self.secure,
            same_site=self.same_site,
            max_age_seconds=int(self.lifetime.total_seconds()),
        )


def request_requires_csrf(method: object) -> bool:
    """Return whether a cookie-authenticated HTTP method requires a CSRF token."""

    if not isinstance(method, str) or _HTTP_METHOD.fullmatch(method) is None:
        raise ValueError("HTTP method is invalid")
    return method.upper() not in _SAFE_METHODS


@dataclass(frozen=True)
class LoginRateLimitPolicy:
    """Bounded in-process limits for one API process on a single LAN host."""

    subject_attempts: int = 5
    client_attempts: int = 30
    window_seconds: float = 300.0
    max_subjects: int = 4096
    max_clients: int = 1024

    def __post_init__(self) -> None:
        for name in (
            "subject_attempts",
            "client_attempts",
            "max_subjects",
            "max_clients",
        ):
            value = getattr(self, name)
            if not isinstance(value, int) or isinstance(value, bool) or value <= 0:
                raise ValueError(f"{name} must be a positive integer")
        if self.subject_attempts > 10_000 or self.client_attempts > 10_000:
            raise ValueError("attempt limits must not exceed 10000")
        if self.max_subjects > 100_000 or self.max_clients > 100_000:
            raise ValueError("tracked identity limits must not exceed 100000")
        if (
            not isinstance(self.window_seconds, (int, float))
            or isinstance(self.window_seconds, bool)
            or not math.isfinite(self.window_seconds)
            or self.window_seconds <= 0
        ):
            raise ValueError("window_seconds must be positive and finite")
        if self.window_seconds > 86_400:
            raise ValueError("window_seconds must not exceed one day")


class BoundedLoginRateLimiter:
    """Thread-safe two-scope limiter that retains no raw client or username."""

    def __init__(
        self,
        policy: LoginRateLimitPolicy | None = None,
        *,
        clock: Callable[[], float] = monotonic,
    ) -> None:
        self.policy = policy or LoginRateLimitPolicy()
        self._clock = clock
        self._clients: OrderedDict[str, deque[float]] = OrderedDict()
        self._subjects: OrderedDict[str, deque[float]] = OrderedDict()
        self._lock = Lock()

    def allow_attempt(self, client_identity: object, username: object) -> bool:
        """Atomically consume one attempt if both client and subject scopes allow it."""

        client = _bounded_login_identity(client_identity)
        try:
            subject = normalize_username(username)
        except ValueError:
            subject = "<invalid>"
        client_key = _rate_limit_fingerprint("client", client)
        subject_key = _rate_limit_fingerprint("subject", subject)
        now = float(self._clock())
        if not math.isfinite(now):
            raise ValueError("rate limiter clock must be finite")
        with self._lock:
            self._purge_expired(self._clients, now)
            self._purge_expired(self._subjects, now)
            client_attempts = self._active_attempts(
                self._clients,
                client_key,
                now,
            )
            subject_attempts = self._active_attempts(
                self._subjects,
                subject_key,
                now,
            )
            if (
                len(client_attempts) >= self.policy.client_attempts
                or len(subject_attempts) >= self.policy.subject_attempts
                or (
                    client_key not in self._clients
                    and len(self._clients) >= self.policy.max_clients
                )
                or (
                    subject_key not in self._subjects
                    and len(self._subjects) >= self.policy.max_subjects
                )
            ):
                return False
            client_attempts.append(now)
            subject_attempts.append(now)
            self._store(
                self._clients,
                client_key,
                client_attempts,
            )
            self._store(
                self._subjects,
                subject_key,
                subject_attempts,
            )
            return True

    @property
    def tracked_key_count(self) -> int:
        """Return a safe size metric used by bounded-memory tests and diagnostics."""

        with self._lock:
            return len(self._clients) + len(self._subjects)

    def _active_attempts(
        self,
        entries: OrderedDict[str, deque[float]],
        key: str,
        now: float,
    ) -> deque[float]:
        attempts = entries.get(key, deque())
        cutoff = now - float(self.policy.window_seconds)
        while attempts and attempts[0] <= cutoff:
            attempts.popleft()
        return attempts

    def _purge_expired(
        self,
        entries: OrderedDict[str, deque[float]],
        now: float,
    ) -> None:
        cutoff = now - float(self.policy.window_seconds)
        while entries:
            _key, attempts = next(iter(entries.items()))
            if attempts and attempts[-1] > cutoff:
                break
            entries.popitem(last=False)

    @staticmethod
    def _store(
        entries: OrderedDict[str, deque[float]],
        key: str,
        attempts: deque[float],
    ) -> None:
        if key in entries:
            entries.move_to_end(key)
        entries[key] = attempts


def _password(value: object, *, require_minimum: bool) -> str:
    if not isinstance(value, str):
        raise ValueError("password must be a string")
    if require_minimum and len(value) < PASSWORD_MIN_CHARACTERS:
        raise ValueError("password is too short")
    if len(value) > PASSWORD_MAX_CHARACTERS:
        raise ValueError("password is too long")
    try:
        encoded = value.encode("utf-8")
    except UnicodeEncodeError:
        raise ValueError("password is not valid UTF-8") from None
    if len(encoded) > PASSWORD_MAX_UTF8_BYTES:
        raise ValueError("password is too long")
    return value


def _password_hasher(parameters: Argon2idParameters) -> PasswordHasher:
    return PasswordHasher(
        time_cost=parameters.time_cost,
        memory_cost=parameters.memory_cost_kib,
        parallelism=parameters.parallelism,
        hash_len=parameters.hash_length,
        salt_len=parameters.salt_length,
        type=Type.ID,
    )


def _stored_password_hash(
    value: object,
) -> tuple[str, Argon2idParameters] | None:
    try:
        encoded_hash = validate_password_hash(value)
        extracted = extract_parameters(encoded_hash)
        if extracted.type is not Type.ID or extracted.version != ARGON2ID_VERSION:
            return None
        profile = Argon2idParameters(
            time_cost=extracted.time_cost,
            memory_cost_kib=extracted.memory_cost,
            parallelism=extracted.parallelism,
            hash_length=extracted.hash_len,
            salt_length=extracted.salt_len,
        )
    except (InvalidHashError, ValueError):
        return None
    return encoded_hash, profile


def _opaque_token(value: object, name: str) -> str:
    if not isinstance(value, str) or _OPAQUE_TOKEN.fullmatch(value) is None:
        raise ValueError(f"{name} is invalid")
    return value


def _digest_token(domain: bytes, token: str) -> str:
    digest = sha256()
    digest.update(domain)
    digest.update(token.encode("ascii"))
    return digest.hexdigest()


def _session_lifetime(value: object) -> timedelta:
    if not isinstance(value, timedelta):
        raise ValueError("session lifetime must be a timedelta")
    if not MIN_SESSION_LIFETIME <= value <= MAX_SESSION_LIFETIME:
        raise ValueError("session lifetime must be between five minutes and seven days")
    if value.microseconds != 0:
        raise ValueError("session lifetime must use whole seconds")
    return value


def _utc_datetime(value: object, name: str) -> datetime:
    if (
        not isinstance(value, datetime)
        or value.tzinfo is None
        or value.utcoffset() is None
    ):
        raise ValueError(f"{name} must be timezone-aware")
    return value.astimezone(timezone.utc)


def _bounded_login_identity(value: object) -> str:
    if (
        not isinstance(value, str)
        or not value
        or len(value) > 255
        or any(ord(character) < 32 or ord(character) == 127 for character in value)
    ):
        return "<invalid>"
    try:
        value.encode("utf-8")
    except UnicodeEncodeError:
        return "<invalid>"
    return value


def _rate_limit_fingerprint(scope: str, *values: str) -> str:
    digest = sha256()
    digest.update(_RATE_LIMIT_DOMAIN)
    for value in (scope, *values):
        encoded = value.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return digest.hexdigest()
