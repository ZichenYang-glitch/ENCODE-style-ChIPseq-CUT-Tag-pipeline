"""Terminal-email notification migration integrity and reversibility tests."""

from __future__ import annotations

import pytest
from sqlalchemy import inspect, text
from sqlalchemy.exc import IntegrityError

from encode_pipeline.persistence import create_database_engine
from encode_pipeline.persistence.migration_admission import (
    verify_migration_execution_inventory,
)
from encode_pipeline.persistence.migrations import downgrade_database, upgrade_database


PRIOR_REVISION = "20260818_14"
NOTIFICATION_REVISION = "20260827_15"
USER_ID = "usr_" + "1" * 32
OTHER_USER_ID = "usr_" + "2" * 32
PASSWORD_HASH = (
    "$argon2id$v=19$m=65536,t=3,p=4$"
    "pGLQVK/BOQLC0oPSA8RQTg$TdiQWUP9gwBAI8iAXUT7oEtjOPqKjJhqyW0JS8ye/Ag"
)


def test_rev15_is_sole_head_and_preserves_legacy_rows_with_safe_defaults(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url, PRIOR_REVISION)
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        _insert_account(connection, USER_ID, "member")
        _insert_run(connection, "legacy-run")
    engine.dispose()

    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    inspector = inspect(engine)

    inventory = verify_migration_execution_inventory()
    assert inventory.heads == (NOTIFICATION_REVISION,)
    account_columns = {
        column["name"]: column for column in inspector.get_columns("user_accounts")
    }
    assert account_columns["notification_email"]["nullable"] is True
    assert account_columns["terminal_email_enabled"]["nullable"] is False
    run_columns = {column["name"]: column for column in inspector.get_columns("runs")}
    assert run_columns["requested_by_user_id"]["nullable"] is True
    assert any(
        index["name"] == "ix_runs_requested_by_user"
        and index["column_names"] == ["requested_by_user_id"]
        for index in inspector.get_indexes("runs")
    )
    requester_foreign_keys = [
        constraint
        for constraint in inspector.get_foreign_keys("runs")
        if constraint["constrained_columns"] == ["requested_by_user_id"]
    ]
    assert len(requester_foreign_keys) == 1
    assert requester_foreign_keys[0]["referred_table"] == "user_accounts"
    assert requester_foreign_keys[0]["referred_columns"] == ["user_id"]
    assert requester_foreign_keys[0]["options"] in (
        {},
        {"ondelete": "RESTRICT"},
    )

    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            NOTIFICATION_REVISION
        )
        account = connection.execute(
            text(
                "SELECT notification_email, terminal_email_enabled "
                "FROM user_accounts WHERE user_id = :user_id"
            ),
            {"user_id": USER_ID},
        ).one()
        assert account == (None, 1)
        assert (
            connection.scalar(
                text(
                    "SELECT requested_by_user_id FROM runs WHERE run_id = 'legacy-run'"
                )
            )
            is None
        )
        runs_ddl = connection.scalar(
            text("SELECT sql FROM sqlite_master WHERE type = 'table' AND name = 'runs'")
        )
        assert "ON DELETE RESTRICT" in runs_ddl
    engine.dispose()


def test_rev15_allows_duplicate_private_addresses_and_enforces_requester_fk(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        _insert_account(connection, USER_ID, "member-a")
        _insert_account(connection, OTHER_USER_ID, "member-b")
        connection.execute(
            text("UPDATE user_accounts SET notification_email = 'lab@example.org'")
        )
        _insert_run(connection, "owned-run", requested_by_user_id=USER_ID)

    with pytest.raises(IntegrityError):
        with engine.begin() as connection:
            _insert_run(
                connection,
                "orphan-run",
                requested_by_user_id="usr_" + "9" * 32,
            )
    engine.dispose()


def test_rev15_downgrade_preserves_legacy_account_and_run_data(tmp_path) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        _insert_account(connection, USER_ID, "member")
        connection.execute(
            text(
                "UPDATE user_accounts SET notification_email = 'member@example.org', "
                "terminal_email_enabled = 0 WHERE user_id = :user_id"
            ),
            {"user_id": USER_ID},
        )
        _insert_run(connection, "owned-run", requested_by_user_id=USER_ID)
    engine.dispose()

    downgrade_database(database_url, PRIOR_REVISION)
    engine = create_database_engine(database_url)
    inspector = inspect(engine)
    assert "notification_email" not in {
        column["name"] for column in inspector.get_columns("user_accounts")
    }
    assert "terminal_email_enabled" not in {
        column["name"] for column in inspector.get_columns("user_accounts")
    }
    assert "requested_by_user_id" not in {
        column["name"] for column in inspector.get_columns("runs")
    }
    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            PRIOR_REVISION
        )
        assert connection.execute(
            text("SELECT username, role FROM user_accounts WHERE user_id = :user_id"),
            {"user_id": USER_ID},
        ).one() == ("member", "member")
        assert (
            connection.scalar(
                text("SELECT status FROM runs WHERE run_id = 'owned-run'")
            )
            == "queued"
        )
    engine.dispose()


def _insert_account(connection, user_id: str, username: str) -> None:
    connection.execute(
        text(
            "INSERT INTO user_accounts "
            "(user_id, username, role, status, password_hash, created_at, "
            "updated_at, password_changed_at) VALUES "
            "(:user_id, :username, 'member', 'enabled', :password_hash, "
            "'2026-08-27 10:00:00', '2026-08-27 10:00:00', "
            "'2026-08-27 10:00:00')"
        ),
        {"user_id": user_id, "username": username, "password_hash": PASSWORD_HASH},
    )


def _insert_run(
    connection,
    run_id: str,
    *,
    requested_by_user_id: str | None = None,
) -> None:
    columns = "run_id, workflow_id, inputs, status, created_at, updated_at, tags"
    values = (
        ":run_id, 'workflow-a', '{}', 'queued', "
        "'2026-08-27 10:00:00', '2026-08-27 10:00:00', '{}'"
    )
    parameters = {"run_id": run_id}
    if requested_by_user_id is not None:
        columns += ", requested_by_user_id"
        values += ", :requested_by_user_id"
        parameters["requested_by_user_id"] = requested_by_user_id
    connection.execute(
        text(f"INSERT INTO runs ({columns}) VALUES ({values})"),
        parameters,
    )
