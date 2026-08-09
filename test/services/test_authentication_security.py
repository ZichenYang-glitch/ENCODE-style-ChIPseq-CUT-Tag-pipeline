"""Focused service tests for local-password and browser-session security."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

import pytest
from argon2 import PasswordHasher, extract_parameters
from argon2.low_level import Type

import encode_pipeline.services.authentication as authentication_module
from encode_pipeline.platform.authentication import (
    UserAccount,
    UserRole,
    UserStatus,
    validate_user_id,
)
from encode_pipeline.services.authentication import (
    CSRF_HEADER_NAME,
    CSRF_TOKEN_LENGTH,
    DEFAULT_SESSION_LIFETIME,
    DEVELOPMENT_CSRF_COOKIE_NAME,
    DEVELOPMENT_SESSION_COOKIE_NAME,
    INVALID_LOGIN_CODE,
    INVALID_LOGIN_MESSAGE,
    PASSWORD_MAX_CHARACTERS,
    PASSWORD_MIN_CHARACTERS,
    SECURE_CSRF_COOKIE_NAME,
    SECURE_SESSION_COOKIE_NAME,
    SESSION_TOKEN_LENGTH,
    Argon2idParameters,
    BoundedLoginRateLimiter,
    BrowserSessionCookiePolicy,
    CookieDirective,
    LoginAuthentication,
    LoginRateLimitPolicy,
    PasswordManager,
    SessionSecrets,
    authenticate_login_candidate,
    csrf_request_is_valid,
    csrf_token_matches,
    digest_csrf_token,
    digest_session_token,
    new_session_record,
    new_session_secrets,
    new_user_id,
    request_requires_csrf,
)


USER_ID = "usr_11111111111111111111111111111111"
NOW = datetime(2026, 8, 9, 12, 0, tzinfo=timezone.utc)
PASSWORD = "correct horse battery staple"
WRONG_PASSWORD = "incorrect horse battery staple"
FAST_PARAMETERS = Argon2idParameters(
    time_cost=2,
    memory_cost_kib=19 * 1024,
    parallelism=1,
)


@pytest.fixture(scope="module")
def password_manager() -> PasswordManager:
    return PasswordManager(FAST_PARAMETERS)


@pytest.fixture(scope="module")
def password_hash(password_manager: PasswordManager) -> str:
    return password_manager.hash_password(PASSWORD)


def test_default_password_hash_is_versioned_argon2id() -> None:
    manager = PasswordManager()
    encoded = manager.hash_password(PASSWORD)
    parameters = extract_parameters(encoded)

    assert encoded.startswith("$argon2id$v=19$m=65536,t=3,p=4$")
    assert parameters.type is Type.ID
    assert parameters.version == 19
    assert parameters.time_cost == 3
    assert parameters.memory_cost == 64 * 1024
    assert parameters.parallelism == 4
    assert parameters.hash_len == 32
    assert parameters.salt_len == 16
    verification = manager.verify_password(PASSWORD, encoded)
    assert verification.verified is True
    assert verification.replacement_hash is None


@pytest.mark.parametrize(
    "overrides",
    (
        {"time_cost": 1},
        {"time_cost": True},
        {"memory_cost_kib": 19 * 1024 - 1},
        {"parallelism": 0},
        {"hash_length": 31},
        {"salt_length": 15},
        {"time_cost": 11},
        {"memory_cost_kib": 256 * 1024 + 1},
        {"parallelism": 17},
        {"hash_length": 65},
        {"salt_length": 65},
    ),
)
def test_argon2_parameters_cannot_fall_below_the_security_floor(
    overrides: dict[str, object],
) -> None:
    with pytest.raises(ValueError):
        Argon2idParameters(**overrides)  # type: ignore[arg-type]


def test_password_verification_fails_closed_for_wrong_or_malformed_hashes(
    password_manager: PasswordManager,
    password_hash: str,
) -> None:
    assert password_manager.verify_password(
        WRONG_PASSWORD,
        password_hash,
    ) == password_manager.verify_password(PASSWORD, None)

    for malformed_hash in (
        "not-a-password-hash",
        "$argon2id$broken",
        "$argon2id$" + "x" * 1024,
        "$argon2id$\x00secret",
    ):
        result = password_manager.verify_password(PASSWORD, malformed_hash)
        assert result.verified is False
        assert result.replacement_hash is None


def test_unbounded_or_unapproved_stored_phc_never_reaches_argon_verification(
    password_manager: PasswordManager,
    password_hash: str,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    unbounded = password_hash.replace("m=19456", "m=262145")
    unapproved_hasher = PasswordHasher(
        time_cost=3,
        memory_cost=19 * 1024,
        parallelism=1,
        hash_len=32,
        salt_len=16,
        type=Type.ID,
    )
    unapproved = unapproved_hasher.hash(PASSWORD)
    candidates: list[str] = []
    original_verify = PasswordHasher.verify

    def recording_verify(
        hasher: PasswordHasher,
        encoded_hash: str,
        password: str | bytes,
    ) -> bool:
        candidates.append(encoded_hash)
        return bool(original_verify(hasher, encoded_hash, password))

    monkeypatch.setattr(PasswordHasher, "verify", recording_verify)

    assert password_manager.verify_password(PASSWORD, unbounded).verified is False
    assert unbounded not in candidates
    candidates.clear()
    assert password_manager.verify_password(PASSWORD, unapproved).verified is False
    assert unapproved not in candidates
    assert candidates


def test_successful_verification_upgrades_outdated_parameters() -> None:
    old_hasher = PasswordHasher(
        time_cost=2,
        memory_cost=19 * 1024,
        parallelism=1,
        hash_len=32,
        salt_len=16,
        type=Type.ID,
    )
    old_hash = old_hasher.hash(PASSWORD)
    manager = PasswordManager(
        Argon2idParameters(
            time_cost=3,
            memory_cost_kib=19 * 1024,
            parallelism=1,
        ),
        legacy_parameters=(FAST_PARAMETERS,),
    )

    upgraded = manager.verify_password(PASSWORD, old_hash)

    assert upgraded.verified is True
    assert upgraded.replacement_hash is not None
    assert upgraded.replacement_hash != old_hash
    assert extract_parameters(upgraded.replacement_hash).time_cost == 3
    current = manager.verify_password(PASSWORD, upgraded.replacement_hash)
    assert current.verified is True
    assert current.replacement_hash is None


def test_every_approved_profile_has_the_same_verification_work_shape(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    current_profile = Argon2idParameters(
        time_cost=3,
        memory_cost_kib=19 * 1024,
        parallelism=1,
    )
    legacy_hasher = PasswordHasher(
        time_cost=FAST_PARAMETERS.time_cost,
        memory_cost=FAST_PARAMETERS.memory_cost_kib,
        parallelism=FAST_PARAMETERS.parallelism,
        hash_len=FAST_PARAMETERS.hash_length,
        salt_len=FAST_PARAMETERS.salt_length,
        type=Type.ID,
    )
    legacy_hash = legacy_hasher.hash(PASSWORD)
    manager = PasswordManager(
        current_profile,
        legacy_parameters=(FAST_PARAMETERS,),
    )
    calls: list[Argon2idParameters] = []
    original_verify = PasswordHasher.verify

    def recording_verify(
        hasher: PasswordHasher,
        encoded_hash: str,
        password: str | bytes,
    ) -> bool:
        parameters = extract_parameters(encoded_hash)
        calls.append(
            Argon2idParameters(
                time_cost=parameters.time_cost,
                memory_cost_kib=parameters.memory_cost,
                parallelism=parameters.parallelism,
                hash_length=parameters.hash_len,
                salt_length=parameters.salt_len,
            )
        )
        return bool(original_verify(hasher, encoded_hash, password))

    monkeypatch.setattr(PasswordHasher, "verify", recording_verify)
    expected_shape = [current_profile, FAST_PARAMETERS]

    for account, candidate in (
        (None, PASSWORD),
        (_account(password_hash=legacy_hash), WRONG_PASSWORD),
        (
            _account(password_hash=legacy_hash, status=UserStatus.DISABLED),
            PASSWORD,
        ),
        (
            _account(password_hash=legacy_hash, status=UserStatus.DISABLED),
            WRONG_PASSWORD,
        ),
    ):
        calls.clear()
        result = authenticate_login_candidate(account, candidate, manager)
        assert result == LoginAuthentication()
        assert calls == expected_shape

    monkeypatch.setattr(
        PasswordHasher,
        "hash",
        lambda *_args, **_kwargs: pytest.fail(
            "disabled-account verification must not perform a rehash"
        ),
    )
    disabled = _account(password_hash=legacy_hash, status=UserStatus.DISABLED)
    assert (
        authenticate_login_candidate(disabled, PASSWORD, manager)
        == LoginAuthentication()
    )


def test_password_profile_allowlist_is_small_and_unique() -> None:
    with pytest.raises(ValueError, match="parameters must"):
        PasswordManager("m=unbounded")  # type: ignore[arg-type]
    with pytest.raises(ValueError, match="unique"):
        PasswordManager(FAST_PARAMETERS, legacy_parameters=(FAST_PARAMETERS,))
    with pytest.raises(ValueError, match="at most two"):
        PasswordManager(
            Argon2idParameters(time_cost=5),
            legacy_parameters=(
                Argon2idParameters(time_cost=2),
                Argon2idParameters(time_cost=3),
                Argon2idParameters(time_cost=4),
            ),
        )


def test_password_input_has_fixed_character_and_utf8_boundaries(
    password_manager: PasswordManager,
) -> None:
    minimum = "a" * PASSWORD_MIN_CHARACTERS
    maximum_utf8 = "🔒" * PASSWORD_MAX_CHARACTERS

    assert password_manager.verify_password(
        minimum,
        password_manager.hash_password(minimum),
    ).verified
    assert password_manager.verify_password(
        maximum_utf8,
        password_manager.hash_password(maximum_utf8),
    ).verified

    for invalid in (
        None,
        "a" * (PASSWORD_MIN_CHARACTERS - 1),
        "a" * (PASSWORD_MAX_CHARACTERS + 1),
        "a" * (PASSWORD_MIN_CHARACTERS - 1) + "\ud800",
    ):
        with pytest.raises(ValueError):
            password_manager.hash_password(invalid)
        assert password_manager.verify_password(invalid, None).verified is False


def test_login_rejections_do_not_reveal_account_existence_or_status(
    password_manager: PasswordManager,
    password_hash: str,
) -> None:
    enabled = _account(password_hash=password_hash)
    disabled = _account(
        password_hash=password_hash,
        status=UserStatus.DISABLED,
    )

    absent = authenticate_login_candidate(None, PASSWORD, password_manager)
    mismatched = authenticate_login_candidate(
        enabled,
        WRONG_PASSWORD,
        password_manager,
    )
    disabled_result = authenticate_login_candidate(
        disabled,
        PASSWORD,
        password_manager,
    )

    assert absent == mismatched == disabled_result == LoginAuthentication()
    assert absent.authenticated is False
    assert INVALID_LOGIN_CODE == "INVALID_CREDENTIALS"
    assert INVALID_LOGIN_MESSAGE == "The username or password is invalid."
    assert enabled.username not in INVALID_LOGIN_MESSAGE
    assert disabled.status.value not in INVALID_LOGIN_MESSAGE

    accepted = authenticate_login_candidate(enabled, PASSWORD, password_manager)
    assert accepted.authenticated is True
    assert accepted.principal is not None
    assert accepted.principal.user_id == enabled.user_id
    assert accepted.principal.username == enabled.username


def test_new_user_id_is_random_opaque_and_path_safe() -> None:
    first = new_user_id()
    second = new_user_id()

    assert first != second
    assert validate_user_id(first) == first
    assert len(first) == 36
    assert "/" not in first and "\\" not in first


def test_new_session_values_are_independent_and_only_digests_are_persisted() -> None:
    first = new_session_secrets()
    second = new_session_secrets()
    first_record = new_session_record(
        user_id=USER_ID,
        secrets=first,
        created_at=NOW,
    )
    second_record = new_session_record(
        user_id=USER_ID,
        secrets=second,
        created_at=NOW,
    )

    assert len(first.session_token) == SESSION_TOKEN_LENGTH
    assert len(first.csrf_token) == CSRF_TOKEN_LENGTH
    assert first.session_token != first.csrf_token
    assert first.session_token != second.session_token
    assert first.csrf_token != second.csrf_token
    assert first_record.session_digest == digest_session_token(first.session_token)
    assert first_record.csrf_digest == digest_csrf_token(first.csrf_token)
    assert first_record.session_digest != second_record.session_digest
    assert first_record.csrf_digest != second_record.csrf_digest
    assert first_record.expires_at == NOW + DEFAULT_SESSION_LIFETIME
    assert first.session_token not in repr(first)
    assert first.csrf_token not in repr(first)
    assert first.session_token not in repr(first_record)
    assert first.csrf_token not in repr(first_record)
    assert not hasattr(first_record, "session_token")
    assert not hasattr(first_record, "csrf_token")

    with pytest.raises(ValueError, match="timezone-aware"):
        new_session_record(
            user_id=USER_ID,
            secrets=first,
            created_at=NOW.replace(tzinfo=None),
        )


def test_session_secrets_cannot_be_caller_selected_or_equal() -> None:
    with pytest.raises(TypeError, match="secure generator"):
        SessionSecrets(session_token="A" * 43, csrf_token="B" * 43)
    with pytest.raises(ValueError, match="independent"):
        SessionSecrets._from_generated("A" * 43, "A" * 43)


def test_csrf_secret_is_regenerated_if_the_rng_repeats(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    values = iter(("A" * 43, "A" * 43, "B" * 43))
    monkeypatch.setattr(
        authentication_module,
        "token_urlsafe",
        lambda _byte_count: next(values),
    )

    secrets = new_session_secrets()

    assert secrets.session_token == "A" * 43
    assert secrets.csrf_token == "B" * 43

    monkeypatch.setattr(
        authentication_module,
        "token_urlsafe",
        lambda _byte_count: "C" * 43,
    )
    with pytest.raises(RuntimeError, match="independent session values") as error:
        new_session_secrets()
    assert "C" * 43 not in str(error.value)


def test_session_and_csrf_digests_are_domain_separated_and_tamper_evident() -> None:
    secrets = new_session_secrets()
    tampered_csrf = _tamper(secrets.csrf_token)
    expected_csrf_digest = digest_csrf_token(secrets.csrf_token)

    assert digest_session_token(secrets.session_token) != digest_csrf_token(
        secrets.session_token
    )
    assert csrf_token_matches(
        secrets.csrf_token,
        expected_csrf_digest,
    )
    assert not csrf_token_matches(
        tampered_csrf,
        expected_csrf_digest,
    )
    assert not csrf_token_matches(None, expected_csrf_digest)
    assert csrf_request_is_valid(
        secrets.csrf_token,
        secrets.csrf_token,
        expected_csrf_digest,
    )
    assert not csrf_request_is_valid(
        secrets.csrf_token,
        tampered_csrf,
        expected_csrf_digest,
    )
    assert not csrf_request_is_valid(
        None,
        secrets.csrf_token,
        expected_csrf_digest,
    )
    with pytest.raises(ValueError, match="SHA-256"):
        csrf_token_matches(secrets.csrf_token, "raw-csrf-token")
    assert digest_session_token(_tamper(secrets.session_token)) != digest_session_token(
        secrets.session_token
    )


@pytest.mark.parametrize(
    "token",
    (None, "", "../session", "a" * 42, "a" * 44, "🔒" * 43),
)
def test_opaque_session_inputs_reject_malicious_or_wrong_length_values(
    token: object,
) -> None:
    with pytest.raises(ValueError):
        digest_session_token(token)


def test_browser_cookie_policy_keeps_auth_cookie_httponly_and_csrf_readable() -> None:
    policy = BrowserSessionCookiePolicy(
        secure=True,
        same_site="strict",
        lifetime=timedelta(hours=1),
    )

    assert SECURE_SESSION_COOKIE_NAME == "__Host-helixweave_session"
    assert SECURE_CSRF_COOKIE_NAME == "__Host-helixweave_csrf"
    assert CSRF_HEADER_NAME == "X-CSRF-Token"
    assert policy.session_cookie.http_only is True
    assert policy.csrf_cookie.http_only is False
    assert policy.session_cookie.secure is True
    assert policy.csrf_cookie.secure is True
    assert policy.session_cookie.max_age_seconds == 3600
    assert policy.csrf_cookie.max_age_seconds == 3600
    assert policy.session_cookie.to_response_kwargs() == {
        "key": SECURE_SESSION_COOKIE_NAME,
        "max_age": 3600,
        "path": "/",
        "secure": True,
        "httponly": True,
        "samesite": "strict",
    }
    assert "domain" not in policy.session_cookie.to_response_kwargs()
    assert "value" not in policy.session_cookie.to_response_kwargs()

    development = BrowserSessionCookiePolicy(secure=False)
    assert development.session_cookie.name == DEVELOPMENT_SESSION_COOKIE_NAME
    assert development.csrf_cookie.name == DEVELOPMENT_CSRF_COOKIE_NAME
    assert development.session_cookie.secure is False

    with pytest.raises(ValueError, match="__Host-"):
        CookieDirective(
            name=DEVELOPMENT_SESSION_COOKIE_NAME,
            http_only=True,
            secure=True,
            same_site="lax",
            max_age_seconds=60,
        )
    with pytest.raises(ValueError, match="HttpOnly"):
        CookieDirective(
            name=SECURE_SESSION_COOKIE_NAME,
            http_only=False,
            secure=True,
            same_site="lax",
            max_age_seconds=60,
        )


@pytest.mark.parametrize("method", ("GET", "head", "Options"))
def test_safe_cookie_authenticated_methods_do_not_require_csrf(method: str) -> None:
    assert request_requires_csrf(method) is False


@pytest.mark.parametrize("method", ("POST", "put", "PATCH", "delete"))
def test_side_effecting_cookie_authenticated_methods_require_csrf(method: str) -> None:
    assert request_requires_csrf(method) is True


@pytest.mark.parametrize("method", (None, "", "x" * 17, " POST", "P0ST"))
def test_csrf_method_parser_rejects_unbounded_or_non_string_input(
    method: object,
) -> None:
    with pytest.raises(ValueError, match="HTTP method"):
        request_requires_csrf(method)


@pytest.mark.parametrize(
    "lifetime",
    (
        timedelta(minutes=5) - timedelta(seconds=1),
        timedelta(days=7) + timedelta(seconds=1),
        timedelta(minutes=5, microseconds=1),
    ),
)
def test_cookie_lifetime_is_explicit_and_bounded(lifetime: timedelta) -> None:
    with pytest.raises(ValueError, match="session lifetime"):
        BrowserSessionCookiePolicy(secure=True, lifetime=lifetime)


def test_login_rate_limiter_enforces_subject_client_and_window_limits() -> None:
    clock = _Clock(100.0)
    limiter = BoundedLoginRateLimiter(
        LoginRateLimitPolicy(
            subject_attempts=2,
            client_attempts=3,
            window_seconds=10,
            max_subjects=10,
            max_clients=10,
        ),
        clock=clock,
    )

    assert limiter.allow_attempt("192.0.2.1", "Alice") is True
    assert limiter.allow_attempt("192.0.2.1", "alice") is True
    assert limiter.allow_attempt("192.0.2.1", "ALICE") is False
    assert limiter.allow_attempt("192.0.2.1", "bob") is True
    assert limiter.allow_attempt("192.0.2.1", "charlie") is False
    assert limiter.allow_attempt("192.0.2.2", "alice") is False
    assert limiter.allow_attempt("192.0.2.2", "bob") is True

    clock.advance(9)
    assert limiter.allow_attempt("192.0.2.1", "alice") is False
    clock.advance(1)
    assert limiter.allow_attempt("192.0.2.1", "alice") is True


def test_login_rate_limiter_retains_only_bounded_fingerprints() -> None:
    limiter = BoundedLoginRateLimiter(
        LoginRateLimitPolicy(
            subject_attempts=100,
            client_attempts=100,
            window_seconds=60,
            max_subjects=3,
            max_clients=2,
        ),
        clock=lambda: 1.0,
    )

    outcomes = [
        limiter.allow_attempt("client-0", f"member{index}") for index in range(20)
    ]
    assert outcomes[:3] == [True, True, True]
    assert not any(outcomes[3:])
    assert limiter.allow_attempt("client-1", "member0")
    assert not limiter.allow_attempt("client-2", "member0")

    assert limiter.tracked_key_count == 5
    retained_keys = (*limiter._clients, *limiter._subjects)
    assert all(
        len(key) == 64 and set(key) <= set("0123456789abcdef") for key in retained_keys
    )
    assert all("client" not in key and "member" not in key for key in retained_keys)


def test_rate_limit_capacity_fails_closed_without_evicting_a_live_subject() -> None:
    clock = _Clock(100.0)
    limiter = BoundedLoginRateLimiter(
        LoginRateLimitPolicy(
            subject_attempts=2,
            client_attempts=100,
            window_seconds=10,
            max_subjects=2,
            max_clients=10,
        ),
        clock=clock,
    )

    assert limiter.allow_attempt("client-1", "alice")
    assert limiter.allow_attempt("client-1", "alice")
    assert not limiter.allow_attempt("client-1", "alice")
    assert limiter.allow_attempt("client-2", "bob")
    assert not limiter.allow_attempt("client-3", "charlie")
    assert not limiter.allow_attempt("client-4", "alice")

    clock.advance(10)
    assert limiter.allow_attempt("client-3", "charlie")


def test_malicious_rate_limit_identities_share_bounded_failure_buckets() -> None:
    limiter = BoundedLoginRateLimiter(
        LoginRateLimitPolicy(
            subject_attempts=10,
            client_attempts=10,
            window_seconds=60,
            max_subjects=10,
            max_clients=10,
        ),
        clock=lambda: 1.0,
    )

    assert limiter.allow_attempt(None, None)
    assert limiter.allow_attempt("x" * 10_000, "y" * 10_000)
    assert limiter.allow_attempt("bad\x00client", "bad\x00username")
    assert limiter.tracked_key_count == 2


@pytest.mark.parametrize(
    "overrides",
    (
        {"subject_attempts": 0},
        {"client_attempts": True},
        {"window_seconds": 0},
        {"window_seconds": float("inf")},
        {"window_seconds": float("nan")},
        {"max_subjects": 0},
        {"max_clients": 0},
        {"subject_attempts": 10_001},
        {"client_attempts": 10_001},
        {"window_seconds": 86_401},
        {"max_subjects": 100_001},
        {"max_clients": 100_001},
    ),
)
def test_rate_limit_policy_rejects_invalid_resource_bounds(
    overrides: dict[str, object],
) -> None:
    with pytest.raises(ValueError):
        LoginRateLimitPolicy(**overrides)  # type: ignore[arg-type]


def test_rate_limiter_rejects_a_non_finite_clock() -> None:
    limiter = BoundedLoginRateLimiter(clock=lambda: float("nan"))

    with pytest.raises(ValueError, match="clock must be finite"):
        limiter.allow_attempt("client", "member")


def _account(
    *,
    password_hash: str,
    status: UserStatus = UserStatus.ENABLED,
) -> UserAccount:
    return UserAccount(
        user_id=USER_ID,
        username="member",
        role=UserRole.MEMBER,
        status=status,
        password_hash=password_hash,
        created_at=NOW,
        updated_at=NOW,
        password_changed_at=NOW,
    )


def _tamper(token: str) -> str:
    replacement = "A" if token[0] != "A" else "B"
    return replacement + token[1:]


class _Clock:
    def __init__(self, value: float) -> None:
        self.value = value

    def __call__(self) -> float:
        return self.value

    def advance(self, seconds: float) -> None:
        self.value += seconds
