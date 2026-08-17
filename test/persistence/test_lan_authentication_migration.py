"""LAN authentication migration integrity and reversibility tests."""

from __future__ import annotations

import pytest
from sqlalchemy import inspect, text
from sqlalchemy.exc import IntegrityError

from encode_pipeline.persistence import create_database_engine
from encode_pipeline.persistence.migration_admission import (
    verify_migration_execution_inventory,
)
from encode_pipeline.persistence.migrations import (
    downgrade_database,
    upgrade_database,
)


PRIOR_REVISION = "20260809_13"
AUTHENTICATION_REVISION = "20260818_14"
AUTHENTICATION_TABLES = {
    "auth_sessions",
    "security_audit_events",
    "user_accounts",
}
USER_ACCOUNT_COLUMNS = {
    "user_id",
    "username",
    "role",
    "status",
    "password_hash",
    "created_at",
    "updated_at",
    "password_changed_at",
}
AUTH_SESSION_COLUMNS = {
    "session_digest",
    "csrf_digest",
    "user_id",
    "created_at",
    "expires_at",
    "revoked_at",
    "revocation_reason",
}
SECURITY_AUDIT_COLUMNS = {
    "event_id",
    "occurred_at",
    "action",
    "outcome",
    "actor_kind",
    "actor_user_id",
    "resource_kind",
    "resource_id",
    "reason_code",
}
USER_ACCOUNT_CONSTRAINTS = {
    "ck_user_accounts_id",
    "ck_user_accounts_password_change_order",
    "ck_user_accounts_role",
    "ck_user_accounts_status",
    "ck_user_accounts_update_order",
    "ck_user_accounts_username",
}
AUTH_SESSION_CONSTRAINTS = {
    "ck_auth_sessions_csrf_digest",
    "ck_auth_sessions_expiry_order",
    "ck_auth_sessions_id",
    "ck_auth_sessions_revocation_order",
    "ck_auth_sessions_revocation_pair",
    "ck_auth_sessions_revocation_reason",
}
SECURITY_AUDIT_CONSTRAINTS = {
    "ck_security_audit_events_action",
    "ck_security_audit_events_actor_kind",
    "ck_security_audit_events_actor_pair",
    "ck_security_audit_events_id",
    "ck_security_audit_events_outcome",
    "ck_security_audit_events_outcome_reason_pair",
    "ck_security_audit_events_reason_code",
    "ck_security_audit_events_resource_kind",
    "ck_security_audit_events_resource_pair",
}
PASSWORD_HASH = (
    "$argon2id$v=19$m=65536,t=3,p=4$"
    "pGLQVK/BOQLC0oPSA8RQTg$TdiQWUP9gwBAI8iAXUT7oEtjOPqKjJhqyW0JS8ye/Ag"
)
VALID_ACCOUNT = {
    "user_id": "usr_" + "1" * 32,
    "username": "alice",
    "role": "member",
    "status": "enabled",
    "password_hash": PASSWORD_HASH,
    "created_at": "2026-08-18 10:00:00",
    "updated_at": "2026-08-18 10:00:00",
    "password_changed_at": "2026-08-18 10:00:00",
}
VALID_SESSION = {
    "session_digest": "a" * 64,
    "csrf_digest": "b" * 64,
    "user_id": VALID_ACCOUNT["user_id"],
    "created_at": "2026-08-18 10:00:00",
    "expires_at": "2026-08-18 18:00:00",
    "revoked_at": None,
    "revocation_reason": None,
}
VALID_AUDIT_EVENT = {
    "event_id": "aevt_" + "1" * 32,
    "occurred_at": "2026-08-18 10:00:00",
    "action": "auth.login",
    "outcome": "denied",
    "actor_kind": "unauthenticated",
    "actor_user_id": None,
    "resource_kind": None,
    "resource_id": None,
    "reason_code": "INVALID_CREDENTIALS",
}
INVALID_ACCOUNTS = (
    {"username": "ab"},
    {"username": "a" * 65},
    {"role": "owner"},
    {"status": "pending"},
    {"updated_at": "2026-08-18 09:00:00"},
    {"password_changed_at": "2026-08-18 11:00:00"},
)
INVALID_SESSIONS = (
    {"session_digest": "A" * 64},
    {"session_digest": "a" * 63},
    {"csrf_digest": "g" * 64},
    {"expires_at": "2026-08-18 10:00:00"},
    {"revoked_at": "2026-08-18 11:00:00"},
    {"revocation_reason": "logout"},
    {"revoked_at": "2026-08-18 11:00:00", "revocation_reason": "other"},
    {"revoked_at": "2026-08-18 09:00:00", "revocation_reason": "logout"},
)
INVALID_AUDIT_EVENTS = (
    {"event_id": "aevt_" + "1" * 31},
    {"event_id": "usr_" + "1" * 32},
    {"action": "auth.refresh"},
    {"outcome": "unknown"},
    {"actor_kind": "service"},
    {"actor_kind": "user"},
    {"actor_user_id": "usr_" + "2" * 32},
    {"resource_kind": "project", "resource_id": "ares_project_" + "3" * 64},
    {"resource_kind": "run"},
    {"resource_id": "ares_run_" + "3" * 64},
    {"reason_code": "SOME_REASON"},
    {"outcome": "succeeded"},
    {"outcome": "failed", "reason_code": None},
)


def test_rev14_is_the_sole_head_and_creates_empty_authentication_schema(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url, PRIOR_REVISION)
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        connection.execute(
            text(
                "INSERT INTO runs "
                "(run_id, workflow_id, inputs, status, created_at, updated_at, "
                "tags) VALUES ('queued-run', 'workflow-a', '{}', 'queued', "
                "'2026-08-18 09:00:00', '2026-08-18 09:00:00', '{}')"
            )
        )
    engine.dispose()

    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    inspector = inspect(engine)

    inventory = verify_migration_execution_inventory()
    assert inventory.heads == (AUTHENTICATION_REVISION,)
    assert AUTHENTICATION_TABLES <= set(inspector.get_table_names())
    assert {
        column["name"] for column in inspector.get_columns("user_accounts")
    } == USER_ACCOUNT_COLUMNS
    assert {
        column["name"] for column in inspector.get_columns("auth_sessions")
    } == AUTH_SESSION_COLUMNS
    assert {
        column["name"] for column in inspector.get_columns("security_audit_events")
    } == SECURITY_AUDIT_COLUMNS
    assert {
        constraint["name"]
        for constraint in inspector.get_check_constraints("user_accounts")
    } == USER_ACCOUNT_CONSTRAINTS
    assert {
        constraint["name"]
        for constraint in inspector.get_check_constraints("auth_sessions")
    } == AUTH_SESSION_CONSTRAINTS
    assert {
        constraint["name"]
        for constraint in inspector.get_check_constraints("security_audit_events")
    } == SECURITY_AUDIT_CONSTRAINTS
    assert [
        (constraint["name"], constraint["column_names"])
        for constraint in inspector.get_unique_constraints("user_accounts")
    ] == [("uq_user_accounts_username", ["username"])]
    assert [
        (index["name"], index["column_names"])
        for index in inspector.get_indexes("auth_sessions")
    ] == [("ix_auth_sessions_user", ["user_id"])]
    assert [
        (index["name"], index["column_names"])
        for index in inspector.get_indexes("security_audit_events")
    ] == [("ix_security_audit_events_occurred", ["occurred_at", "event_id"])]
    [session_foreign_key] = inspector.get_foreign_keys("auth_sessions")
    assert session_foreign_key["name"] == "fk_auth_sessions_user"
    assert session_foreign_key["constrained_columns"] == ["user_id"]
    assert session_foreign_key["referred_table"] == "user_accounts"
    assert session_foreign_key["referred_columns"] == ["user_id"]
    assert session_foreign_key["options"] == {"ondelete": "CASCADE"}
    nullable = {
        column["name"]: column["nullable"]
        for table in ("auth_sessions", "security_audit_events")
        for column in inspector.get_columns(table)
    }
    assert nullable == {
        "session_digest": False,
        "csrf_digest": False,
        "user_id": False,
        "created_at": False,
        "expires_at": False,
        "revoked_at": True,
        "revocation_reason": True,
        "event_id": False,
        "occurred_at": False,
        "action": False,
        "outcome": False,
        "actor_kind": False,
        "actor_user_id": True,
        "resource_kind": True,
        "resource_id": True,
        "reason_code": True,
    }

    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            AUTHENTICATION_REVISION
        )
        assert connection.scalar(text("SELECT count(*) FROM user_accounts")) == 0
        assert connection.scalar(text("SELECT count(*) FROM auth_sessions")) == 0
        assert (
            connection.scalar(text("SELECT count(*) FROM security_audit_events")) == 0
        )
        assert connection.execute(text("SELECT run_id, status FROM runs")).one() == (
            "queued-run",
            "queued",
        )
    engine.dispose()


def test_rev14_enforces_account_session_and_audit_constraints(tmp_path) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)

    with engine.begin() as connection:
        _insert_account(connection, VALID_ACCOUNT)
        _insert_session(connection, VALID_SESSION)
        _insert_audit_event(connection, VALID_AUDIT_EVENT)

    for overrides in INVALID_ACCOUNTS:
        account = dict(VALID_ACCOUNT) | overrides
        with pytest.raises(IntegrityError):
            with engine.begin() as connection:
                _insert_account(connection, account)

    duplicate_username = dict(VALID_ACCOUNT) | {"user_id": "usr_" + "9" * 32}
    with pytest.raises(IntegrityError):
        with engine.begin() as connection:
            _insert_account(connection, duplicate_username)

    for overrides in INVALID_SESSIONS:
        session = dict(VALID_SESSION) | {"session_digest": "d" * 64} | overrides
        with pytest.raises(IntegrityError):
            with engine.begin() as connection:
                _insert_session(connection, session)

    orphan_session = dict(VALID_SESSION) | {
        "session_digest": "e" * 64,
        "user_id": "usr_" + "8" * 32,
    }
    with pytest.raises(IntegrityError):
        with engine.begin() as connection:
            _insert_session(connection, orphan_session)

    for overrides in INVALID_AUDIT_EVENTS:
        event = dict(VALID_AUDIT_EVENT) | {"event_id": "aevt_" + "4" * 32} | overrides
        with pytest.raises(IntegrityError):
            with engine.begin() as connection:
                _insert_audit_event(connection, event)

    with engine.connect() as connection:
        assert connection.scalar(text("SELECT count(*) FROM user_accounts")) == 1
        assert connection.scalar(text("SELECT count(*) FROM auth_sessions")) == 1
        assert (
            connection.scalar(text("SELECT count(*) FROM security_audit_events")) == 1
        )
    engine.dispose()


def test_rev14_downgrade_drops_authentication_schema_and_upgrades_again(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        _insert_account(connection, VALID_ACCOUNT)
        _insert_session(connection, VALID_SESSION)
        _insert_audit_event(connection, VALID_AUDIT_EVENT)
    engine.dispose()

    downgrade_database(database_url, PRIOR_REVISION)
    engine = create_database_engine(database_url)
    assert not (AUTHENTICATION_TABLES & set(inspect(engine).get_table_names()))
    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            PRIOR_REVISION
        )
    engine.dispose()

    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    assert AUTHENTICATION_TABLES <= set(inspect(engine).get_table_names())
    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            AUTHENTICATION_REVISION
        )
        assert connection.scalar(text("SELECT count(*) FROM user_accounts")) == 0
        assert connection.scalar(text("SELECT count(*) FROM auth_sessions")) == 0
        assert (
            connection.scalar(text("SELECT count(*) FROM security_audit_events")) == 0
        )
    engine.dispose()


def _insert_account(connection, account: dict) -> None:
    connection.execute(
        text(
            "INSERT INTO user_accounts "
            "(user_id, username, role, status, password_hash, created_at, "
            "updated_at, password_changed_at) VALUES "
            "(:user_id, :username, :role, :status, :password_hash, :created_at, "
            ":updated_at, :password_changed_at)"
        ),
        account,
    )


def _insert_session(connection, session: dict) -> None:
    connection.execute(
        text(
            "INSERT INTO auth_sessions "
            "(session_digest, csrf_digest, user_id, created_at, expires_at, "
            "revoked_at, revocation_reason) VALUES "
            "(:session_digest, :csrf_digest, :user_id, :created_at, :expires_at, "
            ":revoked_at, :revocation_reason)"
        ),
        session,
    )


def _insert_audit_event(connection, event: dict) -> None:
    connection.execute(
        text(
            "INSERT INTO security_audit_events "
            "(event_id, occurred_at, action, outcome, actor_kind, actor_user_id, "
            "resource_kind, resource_id, reason_code) VALUES "
            "(:event_id, :occurred_at, :action, :outcome, :actor_kind, "
            ":actor_user_id, :resource_kind, :resource_id, :reason_code)"
        ),
        event,
    )
