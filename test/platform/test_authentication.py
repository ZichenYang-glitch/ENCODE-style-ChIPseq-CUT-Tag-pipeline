"""Contracts for local users, roles, and opaque server-side sessions."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

import pytest

from encode_pipeline.platform.authentication import (
    AuthenticatedPrincipal,
    SessionRecord,
    SessionRevocationReason,
    UserAccount,
    UserRole,
    UserStatus,
    authenticated_principal_for_session,
    normalize_username,
    validate_password_hash,
    validate_user_id,
)


USER_ID = "usr_11111111111111111111111111111111"
OTHER_USER_ID = "usr_22222222222222222222222222222222"
SESSION_DIGEST = "a" * 64
CSRF_DIGEST = "b" * 64
PASSWORD_HASH = (
    "$argon2id$v=19$m=65536,t=3,p=4$"
    "pGLQVK/BOQLC0oPSA8RQTg$TdiQWUP9gwBAI8iAXUT7oEtjOPqKjJhqyW0JS8ye/Ag"
)
REPLACEMENT_PASSWORD_HASH = (
    "$argon2id$v=19$m=65536,t=3,p=4$"
    "U9mvrmk+6rhLA+y7XE5MIg$0Tq7KEeNFcEhS2Hykdc8u6Gd/v2ZOzDs8BbxZ42sm/4"
)
NOW = datetime(2026, 8, 9, 12, 0, tzinfo=timezone.utc)


def test_role_status_and_revocation_vocabularies_are_closed() -> None:
    assert {role.value for role in UserRole} == {"administrator", "member"}
    assert {status.value for status in UserStatus} == {"enabled", "disabled"}
    assert {reason.value for reason in SessionRevocationReason} == {
        "logout",
        "all_sessions",
        "account_disabled",
        "password_reset",
    }


@pytest.mark.parametrize(
    ("value", "expected"),
    (
        ("Alice", "alice"),
        (" lab.member ", "lab.member"),
        ("A_B-c", "a_b-c"),
    ),
)
def test_username_normalization_is_bounded_ascii_and_case_insensitive(
    value: str,
    expected: str,
) -> None:
    assert normalize_username(value) == expected


@pytest.mark.parametrize(
    "value",
    (
        None,
        "",
        "ab",
        "a" * 65,
        "1member",
        "member name",
        "member/name",
        "m\x00ember",
        "\tmember\t",
        "\nmember",
        "mémber",
        "\N{NO-BREAK SPACE}member\N{NO-BREAK SPACE}",
        "ｍember",
    ),
)
def test_username_normalization_rejects_malicious_or_ambiguous_input(
    value: object,
) -> None:
    with pytest.raises(ValueError):
        normalize_username(value)


@pytest.mark.parametrize(
    "value",
    (
        "usr_" + "1" * 31,
        "usr_" + "G" * 32,
        "../usr_" + "1" * 32,
        "session-token",
        None,
    ),
)
def test_user_id_is_stable_opaque_and_path_safe(value: object) -> None:
    with pytest.raises(ValueError):
        validate_user_id(value)


def test_account_normalizes_role_status_and_lifecycle() -> None:
    offset = timezone(timedelta(hours=8))
    account = UserAccount(
        user_id=USER_ID,
        username=" Alice ",
        role="administrator",
        status="enabled",
        password_hash=PASSWORD_HASH,
        created_at=NOW.astimezone(offset),
        updated_at=NOW.astimezone(offset),
        password_changed_at=NOW.astimezone(offset),
    )

    assert account.username == "alice"
    assert account.role is UserRole.ADMINISTRATOR
    assert account.status is UserStatus.ENABLED
    assert account.created_at == NOW
    assert account.enabled is True
    assert account.to_public_summary() == {
        "user_id": USER_ID,
        "username": "alice",
        "role": "administrator",
        "status": "enabled",
        "created_at": NOW,
        "updated_at": NOW,
        "password_changed_at": NOW,
    }
    assert PASSWORD_HASH not in repr(account)
    assert "password_hash" not in account.to_public_summary()


@pytest.mark.parametrize(
    "encoded_hash",
    (
        "plaintext-password",
        "$argon2id$v=19$m=65536,t=3,p=4$encoded",
        PASSWORD_HASH.replace("argon2id", "argon2i"),
        PASSWORD_HASH.replace("v=19", "v=18"),
        PASSWORD_HASH.replace("m=65536", "m=262145"),
        PASSWORD_HASH.replace("t=3", "t=11"),
        PASSWORD_HASH.replace("p=4", "p=17"),
        "$argon2id$v=19$m=65536,t=3,p=4$QQ$QQ",
        PASSWORD_HASH + "=",
        "$argon2id$v=19$m=65536,t=3,p=4$" + "A" * 220,
    ),
)
def test_password_hash_contract_rejects_plaintext_malformed_and_unbounded_phc(
    encoded_hash: str,
) -> None:
    with pytest.raises(ValueError, match="password_hash"):
        validate_password_hash(encoded_hash)
    with pytest.raises(ValueError, match="password_hash"):
        _account(password_hash=encoded_hash)


@pytest.mark.parametrize("role", ("owner", "admin", "", None))
def test_account_role_is_a_closed_two_value_contract(role: object) -> None:
    with pytest.raises(ValueError, match="administrator or member"):
        _account(role=role)  # type: ignore[arg-type]


@pytest.mark.parametrize("status", ("locked", "deleted", "", None))
def test_account_status_is_a_closed_enabled_disabled_contract(status: object) -> None:
    with pytest.raises(ValueError, match="enabled or disabled"):
        _account(status=status)  # type: ignore[arg-type]


def test_account_lifecycle_rejects_naive_or_backwards_times() -> None:
    with pytest.raises(ValueError, match="timezone-aware"):
        _account(created_at=NOW.replace(tzinfo=None))
    with pytest.raises(ValueError, match="updated_at"):
        _account(updated_at=NOW - timedelta(seconds=1))
    with pytest.raises(ValueError, match="password_changed_at"):
        _account(password_changed_at=NOW + timedelta(seconds=1))


def test_account_enable_disable_and_password_change_are_immutable() -> None:
    account = _account()
    changed_at = NOW + timedelta(minutes=1)
    disabled = account.disable(changed_at)
    assert account.status is UserStatus.ENABLED
    assert disabled.status is UserStatus.DISABLED
    assert disabled.updated_at == changed_at
    assert disabled.disable(changed_at + timedelta(minutes=1)) is disabled

    enabled = disabled.enable(changed_at + timedelta(minutes=2))
    assert enabled.status is UserStatus.ENABLED
    replacement = enabled.change_password(
        REPLACEMENT_PASSWORD_HASH,
        changed_at + timedelta(minutes=3),
    )
    assert replacement.password_changed_at == changed_at + timedelta(minutes=3)
    assert replacement.updated_at == replacement.password_changed_at
    assert replacement.password_hash != account.password_hash


def test_principal_can_only_be_created_from_an_enabled_account() -> None:
    principal = AuthenticatedPrincipal.from_account(_account())
    assert principal.user_id == USER_ID
    assert principal.is_administrator is False

    with pytest.raises(ValueError, match="enabled account"):
        AuthenticatedPrincipal.from_account(_account(status=UserStatus.DISABLED))


def test_session_record_keeps_only_digests_and_has_explicit_expiry() -> None:
    session = _session()
    assert session.active_at(NOW - timedelta(microseconds=1)) is False
    assert session.active_at(NOW) is True
    assert session.active_at(NOW + timedelta(hours=12) - timedelta(microseconds=1))
    assert session.active_at(NOW + timedelta(hours=12)) is False
    assert SESSION_DIGEST not in repr(session)
    assert CSRF_DIGEST not in repr(session)
    assert not hasattr(session, "session_token")
    assert not hasattr(session, "csrf_token")


@pytest.mark.parametrize(
    ("field", "value"),
    (
        ("session_digest", "a" * 63),
        ("session_digest", "G" * 64),
        ("csrf_digest", "b" * 65),
        ("csrf_digest", "token"),
    ),
)
def test_session_record_rejects_non_digest_secret_storage(
    field: str,
    value: str,
) -> None:
    with pytest.raises(ValueError, match="SHA-256"):
        _session(**{field: value})


def test_session_revocation_is_idempotent_and_blocks_replay() -> None:
    session = _session()
    revoked_at = NOW + timedelta(minutes=5)
    revoked = session.revoke(revoked_at, SessionRevocationReason.LOGOUT)
    assert revoked.revoked_at == revoked_at
    assert revoked.revocation_reason is SessionRevocationReason.LOGOUT
    assert revoked.active_at(revoked_at) is False
    assert (
        revoked.revoke(
            revoked_at + timedelta(minutes=1),
            SessionRevocationReason.ALL_SESSIONS,
        )
        is revoked
    )
    assert (
        authenticated_principal_for_session(
            revoked,
            _account(),
            revoked_at,
        )
        is None
    )


def test_session_record_requires_consistent_revocation_fields() -> None:
    with pytest.raises(ValueError, match="both be set"):
        _session(revoked_at=NOW, revocation_reason=None)
    with pytest.raises(ValueError, match="both be set"):
        _session(
            revoked_at=None,
            revocation_reason=SessionRevocationReason.LOGOUT,
        )
    with pytest.raises(ValueError, match="expires_at"):
        _session(expires_at=NOW)


@pytest.mark.parametrize(
    "expires_at",
    (
        NOW + timedelta(minutes=5) - timedelta(seconds=1),
        NOW + timedelta(days=7, seconds=1),
        NOW + timedelta(minutes=5, microseconds=1),
    ),
)
def test_session_record_enforces_the_server_lifetime_boundary(
    expires_at: datetime,
) -> None:
    with pytest.raises(ValueError, match="five minutes|whole seconds"):
        _session(expires_at=expires_at)


def test_disabled_or_mismatched_account_invalidates_session_immediately() -> None:
    session = _session()
    assert authenticated_principal_for_session(session, _account(), NOW) is not None
    assert (
        authenticated_principal_for_session(
            session,
            _account(status=UserStatus.DISABLED),
            NOW,
        )
        is None
    )
    assert (
        authenticated_principal_for_session(
            session,
            _account(user_id=OTHER_USER_ID),
            NOW,
        )
        is None
    )
    assert authenticated_principal_for_session(session, None, NOW) is None


def _account(**overrides: object) -> UserAccount:
    values: dict[str, object] = {
        "user_id": USER_ID,
        "username": "member",
        "role": UserRole.MEMBER,
        "status": UserStatus.ENABLED,
        "password_hash": PASSWORD_HASH,
        "created_at": NOW,
        "updated_at": NOW,
        "password_changed_at": NOW,
    }
    values.update(overrides)
    return UserAccount(**values)  # type: ignore[arg-type]


def _session(**overrides: object) -> SessionRecord:
    values: dict[str, object] = {
        "session_digest": SESSION_DIGEST,
        "csrf_digest": CSRF_DIGEST,
        "user_id": USER_ID,
        "created_at": NOW,
        "expires_at": NOW + timedelta(hours=12),
    }
    values.update(overrides)
    return SessionRecord(**values)  # type: ignore[arg-type]
