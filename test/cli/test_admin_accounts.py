"""Local operator account administration CLI contracts."""

from __future__ import annotations

import json

import pytest

from encode_pipeline.cli import admin
from encode_pipeline.persistence import (
    SqlAlchemyAuthenticationRepository,
    create_database_engine,
    create_session_factory,
)
from encode_pipeline.platform.security_audit import AuditAction, AuditActorKind
from encode_pipeline.services.authentication_service import AuthenticationService


def _database_url(tmp_path) -> str:
    return f"sqlite:///{tmp_path / 'platform.db'}"


def _reader(password: str):
    calls: list[str] = []

    def read(prompt: str) -> str:
        calls.append(prompt)
        return password

    read.calls = calls
    return read


def _auth_repository(database_url: str) -> SqlAlchemyAuthenticationRepository:
    engine = create_database_engine(database_url)
    return SqlAlchemyAuthenticationRepository(create_session_factory(engine))


def test_bootstrap_creates_the_first_administrator_with_safe_output(
    tmp_path, capsys
) -> None:
    database_url = _database_url(tmp_path)
    password = "correct horse battery staple"
    exit_code = admin.main(
        [
            "--database-url",
            database_url,
            "account",
            "bootstrap",
            "--username",
            "root-admin",
        ],
        password_reader=_reader(password),
    )
    assert exit_code == 0
    output = capsys.readouterr().out
    payload = json.loads(output)
    assert payload["username"] == "root-admin"
    assert payload["role"] == "administrator"
    assert payload["status"] == "enabled"
    assert password not in output
    assert "password_hash" not in output

    repository = _auth_repository(database_url)
    account = repository.get_account_by_username("root-admin")
    assert account.password_hash != password
    (event,) = repository.list_security_audit_events()
    assert event.action is AuditAction.ACCOUNT_CREATE
    assert event.actor_kind is AuditActorKind.LOCAL_OPERATOR

    login = AuthenticationService(repository=repository).login(
        "root-admin",
        password,
        client_identity="127.0.0.1",
    )
    assert login.principal.user_id == account.user_id


def test_bootstrap_rejects_a_second_administrator(tmp_path, capsys) -> None:
    database_url = _database_url(tmp_path)
    arguments = [
        "--database-url",
        database_url,
        "account",
        "bootstrap",
        "--username",
        "root-admin",
    ]
    assert admin.main(arguments, password_reader=_reader("x" * 15)) == 0
    capsys.readouterr()

    second = admin.main(
        [
            "--database-url",
            database_url,
            "account",
            "bootstrap",
            "--username",
            "other-admin",
        ],
        password_reader=_reader("y" * 15),
    )
    assert second == 1
    assert "password_hash" not in capsys.readouterr().out


def test_bootstrap_rejects_mismatched_confirmation(tmp_path, capsys) -> None:
    database_url = _database_url(tmp_path)
    responses = iter(("first password 1", "second password 2"))

    exit_code = admin.main(
        [
            "--database-url",
            database_url,
            "account",
            "bootstrap",
            "--username",
            "root-admin",
        ],
        password_reader=lambda _prompt: next(responses),
    )
    assert exit_code == 1
    repository = _auth_repository(database_url)
    assert repository.list_accounts() == ()


def test_reset_password_rotates_credentials_and_audits_the_operator(
    tmp_path, capsys
) -> None:
    database_url = _database_url(tmp_path)
    assert (
        admin.main(
            [
                "--database-url",
                database_url,
                "account",
                "bootstrap",
                "--username",
                "root-admin",
            ],
            password_reader=_reader("initial password phrase"),
        )
        == 0
    )
    service = AuthenticationService(repository=_auth_repository(database_url))
    session = service.login(
        "root-admin", "initial password phrase", client_identity="127.0.0.1"
    )
    capsys.readouterr()

    exit_code = admin.main(
        [
            "--database-url",
            database_url,
            "account",
            "reset-password",
            "--username",
            "root-admin",
        ],
        password_reader=_reader("a brand new password phrase"),
    )
    assert exit_code == 0
    output = capsys.readouterr().out
    assert "a brand new password phrase" not in output

    repository = _auth_repository(database_url)
    assert service.resolve_session(session.secrets.session_token) is None
    with pytest.raises(Exception):
        service.login(
            "root-admin", "initial password phrase", client_identity="127.0.0.1"
        )
    assert (
        service.login(
            "root-admin",
            "a brand new password phrase",
            client_identity="127.0.0.1",
        ).principal.username
        == "root-admin"
    )

    resets = [
        event
        for event in repository.list_security_audit_events()
        if event.action is AuditAction.ACCOUNT_PASSWORD_RESET
    ]
    assert len(resets) == 1
    assert resets[0].actor_kind is AuditActorKind.LOCAL_OPERATOR


def test_reset_password_rejects_an_unknown_account(tmp_path, capsys) -> None:
    database_url = _database_url(tmp_path)
    admin.main(
        [
            "--database-url",
            database_url,
            "account",
            "bootstrap",
            "--username",
            "root-admin",
        ],
        password_reader=_reader("initial password phrase"),
    )
    capsys.readouterr()

    exit_code = admin.main(
        [
            "--database-url",
            database_url,
            "account",
            "reset-password",
            "--username",
            "unknown-user",
        ],
        password_reader=_reader("a brand new password phrase"),
    )
    assert exit_code == 1


def test_account_list_emits_only_safe_summaries(tmp_path, capsys) -> None:
    database_url = _database_url(tmp_path)
    admin.main(
        [
            "--database-url",
            database_url,
            "account",
            "bootstrap",
            "--username",
            "root-admin",
        ],
        password_reader=_reader("initial password phrase"),
    )
    capsys.readouterr()

    assert admin.main(["--database-url", database_url, "account", "list"]) == 0
    output = capsys.readouterr().out
    payload = json.loads(output)
    assert [account["username"] for account in payload] == ["root-admin"]
    assert "password_hash" not in output
    assert "initial password phrase" not in output
