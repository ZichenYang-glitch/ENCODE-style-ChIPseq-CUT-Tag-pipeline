"""SQLAlchemy LAN authentication repository behavior."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

import pytest
from sqlalchemy.exc import IntegrityError

from encode_pipeline.persistence import (
    SqlAlchemyAuthenticationRepository,
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.platform.authentication import (
    SessionRecord,
    SessionRevocationReason,
    UserAccount,
    UserRole,
    UserStatus,
)
from encode_pipeline.platform.security_audit import (
    AuditAction,
    AuditActorKind,
    AuditOutcome,
    AuditReasonCode,
    AuditResourceKind,
    SecurityAuditEvent,
    build_audit_resource,
    new_audit_event_id,
)
from encode_pipeline.services.authentication_repositories import (
    AuthenticationAccountConflictError,
)


USER_ID = "usr_" + "1" * 32
OTHER_USER_ID = "usr_" + "2" * 32
UNKNOWN_USER_ID = "usr_" + "9" * 32
SESSION_DIGEST = "a" * 64
OTHER_SESSION_DIGEST = "c" * 64
THIRD_SESSION_DIGEST = "d" * 64
CSRF_DIGEST = "b" * 64
PASSWORD_HASH = (
    "$argon2id$v=19$m=65536,t=3,p=4$"
    "pGLQVK/BOQLC0oPSA8RQTg$TdiQWUP9gwBAI8iAXUT7oEtjOPqKjJhqyW0JS8ye/Ag"
)
REPLACEMENT_PASSWORD_HASH = (
    "$argon2id$v=19$m=65536,t=3,p=4$"
    "U9mvrmk+6rhLA+y7XE5MIg$0Tq7KEeNFcEhS2Hykdc8u6Gd/v2ZOzDs8BbxZ42sm/4"
)
NOW = datetime(2026, 8, 18, 12, 0, tzinfo=timezone.utc)
LATER = datetime(2026, 8, 18, 13, 0, tzinfo=timezone.utc)
LATEST = datetime(2026, 8, 18, 14, 0, tzinfo=timezone.utc)


def _account(**overrides: object) -> UserAccount:
    fields: dict[str, object] = {
        "user_id": USER_ID,
        "username": "alice",
        "role": UserRole.MEMBER,
        "status": UserStatus.ENABLED,
        "password_hash": PASSWORD_HASH,
        "created_at": NOW,
        "updated_at": NOW,
        "password_changed_at": NOW,
    }
    fields.update(overrides)
    return UserAccount(**fields)  # type: ignore[arg-type]


def _session(**overrides: object) -> SessionRecord:
    fields: dict[str, object] = {
        "session_digest": SESSION_DIGEST,
        "csrf_digest": CSRF_DIGEST,
        "user_id": USER_ID,
        "created_at": NOW,
        "expires_at": NOW + timedelta(hours=8),
    }
    fields.update(overrides)
    return SessionRecord(**fields)  # type: ignore[arg-type]


def _audit_event(**overrides: object) -> SecurityAuditEvent:
    fields: dict[str, object] = {
        "event_id": new_audit_event_id(),
        "occurred_at": NOW,
        "action": AuditAction.ACCOUNT_CREATE,
        "outcome": AuditOutcome.SUCCEEDED,
        "actor_kind": AuditActorKind.LOCAL_OPERATOR,
        "resource": build_audit_resource(AuditResourceKind.ACCOUNT, USER_ID),
    }
    fields.update(overrides)
    return SecurityAuditEvent(**fields)  # type: ignore[arg-type]


@pytest.fixture
def repository(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'authentication.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    yield SqlAlchemyAuthenticationRepository(create_session_factory(engine))
    engine.dispose()


def test_account_round_trip_and_lookup_paths(repository) -> None:
    account = _account()
    assert repository.create_account(account) == account
    other = _account(
        user_id=OTHER_USER_ID,
        username="bob",
        role=UserRole.ADMINISTRATOR,
    )
    repository.create_account(other)

    assert repository.get_account_by_id(USER_ID) == account
    assert repository.get_account_by_username(" Alice ") == account
    assert repository.list_accounts() == (account, other)

    with pytest.raises(KeyError):
        repository.get_account_by_id(UNKNOWN_USER_ID)
    with pytest.raises(KeyError):
        repository.get_account_by_username("unknown")


def test_save_account_persists_state_and_requires_existing_account(
    repository,
) -> None:
    repository.create_account(_account())
    changed = _account().change_password(REPLACEMENT_PASSWORD_HASH, LATER)
    assert repository.save_account(changed) == 0
    assert repository.get_account_by_id(USER_ID) == changed

    disabled = changed.disable(LATEST)
    assert repository.save_account(disabled) == 0
    assert repository.get_account_by_id(USER_ID) == disabled
    assert repository.get_account_by_username("alice") == disabled

    with pytest.raises(KeyError):
        repository.save_account(_account(user_id=UNKNOWN_USER_ID))


def test_create_account_rejects_duplicate_identity_and_username(repository) -> None:
    repository.create_account(_account())
    with pytest.raises(AuthenticationAccountConflictError):
        repository.create_account(_account())
    with pytest.raises(AuthenticationAccountConflictError):
        repository.create_account(_account(user_id=OTHER_USER_ID, username=" ALICE "))


def test_has_enabled_administrator_tracks_account_state(repository) -> None:
    assert repository.has_enabled_administrator() is False
    repository.create_account(_account())
    assert repository.has_enabled_administrator() is False

    administrator = _account(
        user_id=OTHER_USER_ID,
        username="operator",
        role=UserRole.ADMINISTRATOR,
    )
    repository.create_account(administrator)
    assert repository.has_enabled_administrator() is True

    repository.save_account(administrator.disable(LATER))
    assert repository.has_enabled_administrator() is False

    repository.create_account(
        _account(
            user_id=UNKNOWN_USER_ID,
            username="second.op",
            role=UserRole.ADMINISTRATOR,
            status=UserStatus.DISABLED,
        )
    )
    assert repository.has_enabled_administrator() is False


def test_session_lifecycle_and_idempotent_revoke(repository) -> None:
    repository.create_account(_account())
    record = _session()
    repository.create_session(record)
    assert repository.get_session(SESSION_DIGEST) == record

    with pytest.raises(KeyError):
        repository.get_session(OTHER_SESSION_DIGEST)

    assert (
        repository.revoke_session(
            SESSION_DIGEST,
            LATER,
            SessionRevocationReason.LOGOUT,
        )
        is True
    )
    revoked = repository.get_session(SESSION_DIGEST)
    assert revoked.revoked_at == LATER
    assert revoked.revocation_reason is SessionRevocationReason.LOGOUT
    assert revoked.active_at(LATER) is False

    assert (
        repository.revoke_session(
            SESSION_DIGEST,
            LATEST,
            SessionRevocationReason.PASSWORD_RESET,
        )
        is False
    )
    assert repository.get_session(SESSION_DIGEST) == revoked

    with pytest.raises(KeyError):
        repository.revoke_session(
            OTHER_SESSION_DIGEST,
            LATER,
            SessionRevocationReason.LOGOUT,
        )


def test_create_session_requires_an_existing_account(repository) -> None:
    with pytest.raises(IntegrityError):
        repository.create_session(_session(user_id=UNKNOWN_USER_ID))
    with pytest.raises(KeyError):
        repository.get_session(SESSION_DIGEST)


def test_save_account_revokes_only_that_users_active_sessions(repository) -> None:
    repository.create_account(_account())
    repository.create_account(_account(user_id=OTHER_USER_ID, username="bob"))
    repository.create_session(_session())
    repository.create_session(_session(session_digest=OTHER_SESSION_DIGEST))
    repository.create_session(
        _session(
            session_digest=THIRD_SESSION_DIGEST,
            revoked_at=NOW,
            revocation_reason=SessionRevocationReason.LOGOUT,
        )
    )
    repository.create_session(
        _session(
            session_digest="e" * 64,
            user_id=OTHER_USER_ID,
        )
    )

    disabled = _account().disable(LATER)
    revoked_count = repository.save_account(
        disabled,
        revoke_sessions_reason=SessionRevocationReason.ACCOUNT_DISABLED,
        revoked_at=LATER,
    )
    assert revoked_count == 2
    for digest in (SESSION_DIGEST, OTHER_SESSION_DIGEST):
        revoked = repository.get_session(digest)
        assert revoked.revoked_at == LATER
        assert revoked.revocation_reason is SessionRevocationReason.ACCOUNT_DISABLED
    prerevoked = repository.get_session(THIRD_SESSION_DIGEST)
    assert prerevoked.revoked_at == NOW
    assert prerevoked.revocation_reason is SessionRevocationReason.LOGOUT
    other_user_session = repository.get_session("e" * 64)
    assert other_user_session.revoked_at is None
    assert other_user_session.revocation_reason is None
    assert repository.get_account_by_id(USER_ID) == disabled

    assert (
        repository.save_account(
            disabled,
            revoke_sessions_reason=SessionRevocationReason.ACCOUNT_DISABLED,
            revoked_at=LATER,
        )
        == 0
    )

    with pytest.raises(ValueError):
        repository.save_account(
            disabled,
            revoke_sessions_reason=SessionRevocationReason.ACCOUNT_DISABLED,
        )
    with pytest.raises(ValueError):
        repository.save_account(disabled, revoked_at=LATER)


def test_audit_rows_commit_or_roll_back_with_business_writes(repository) -> None:
    account = _account()
    created = _audit_event()
    repository.create_account(account, audit=created)
    repository.create_session(_session())
    assert repository.list_security_audit_events() == (created,)

    other = _account(user_id=OTHER_USER_ID, username="bob")
    with pytest.raises(AuthenticationAccountConflictError):
        repository.create_account(other, audit=created)
    with pytest.raises(KeyError):
        repository.get_account_by_id(OTHER_USER_ID)
    assert repository.list_security_audit_events() == (created,)

    disabled = account.disable(LATER)
    with pytest.raises(AuthenticationAccountConflictError):
        repository.save_account(disabled, audit=created)
    assert repository.get_account_by_id(USER_ID) == account
    assert repository.list_security_audit_events() == (created,)

    with pytest.raises(IntegrityError):
        repository.revoke_session(
            SESSION_DIGEST,
            LATER,
            SessionRevocationReason.LOGOUT,
            audit=created,
        )
    assert repository.get_session(SESSION_DIGEST).revoked_at is None
    assert repository.list_security_audit_events() == (created,)


def test_audit_event_round_trip_orders_recent_first_and_bounds_limit(
    repository,
) -> None:
    repository.create_account(_account())
    created = _audit_event()
    repository.create_account(
        _account(user_id=OTHER_USER_ID, username="bob"),
        audit=created,
    )
    repository.create_session(
        _session(),
        audit=_audit_event(
            occurred_at=LATER,
            action=AuditAction.LOGIN,
            actor_kind=AuditActorKind.USER,
            actor_user_id=USER_ID,
            resource=None,
        ),
    )
    revoked = _audit_event(
        occurred_at=LATEST,
        action=AuditAction.LOGOUT,
        actor_kind=AuditActorKind.USER,
        actor_user_id=USER_ID,
        resource=None,
    )
    repository.revoke_session(
        SESSION_DIGEST,
        LATEST,
        SessionRevocationReason.LOGOUT,
        audit=revoked,
    )

    events = repository.list_security_audit_events()
    assert [event.occurred_at for event in events] == [LATEST, LATER, NOW]
    assert events[0] == revoked
    assert events[2] == created
    assert [
        event.occurred_at for event in repository.list_security_audit_events(limit=2)
    ] == [LATEST, LATER]
    assert (
        repository.list_security_audit_events(limit=1000)
        == repository.list_security_audit_events(limit=3)
        == events
    )
    for invalid_limit in (0, -1, 1001, True):
        with pytest.raises(ValueError):
            repository.list_security_audit_events(limit=invalid_limit)


def test_failed_login_audit_round_trip_keeps_closed_shape(repository) -> None:
    repository.create_account(_account())
    denied = _audit_event(
        action=AuditAction.LOGIN,
        outcome=AuditOutcome.DENIED,
        actor_kind=AuditActorKind.UNAUTHENTICATED,
        resource=None,
        reason_code=AuditReasonCode.INVALID_CREDENTIALS,
    )
    repository.create_session(_session(), audit=denied)

    [event] = repository.list_security_audit_events()
    assert event == denied
    assert event.actor_user_id is None
    assert event.resource is None
    assert event.reason_code is AuditReasonCode.INVALID_CREDENTIALS
