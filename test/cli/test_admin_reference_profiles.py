"""Administrator-only Reference Profile CLI contract tests."""

from __future__ import annotations

from contextlib import contextmanager
import hashlib
import json
from pathlib import Path
import sqlite3

import pytest

from encode_pipeline.cli import admin
from encode_pipeline.services.authentication_service import (
    AuthenticationActor,
)
from encode_pipeline.persistence.migrations import upgrade_database
from encode_pipeline.services.reference_profiles import ReferenceProfileAdminError


DATABASE_URL = "sqlite:////tmp/helixweave-reference-profile-test.db"


class RecordingReferenceProfiles:
    def __init__(self) -> None:
        self.calls: list[tuple[str, dict[str, object]]] = []

    def register(self, **kwargs):
        self.calls.append(("register", kwargs))
        return {
            "profile_id": "refp_" + "1" * 32,
            "revision_id": "refpr_" + "2" * 32,
            "public_identity_sha256": "3" * 64,
        }

    def verify(self, revision_id: str):
        self.calls.append(("verify", {"revision_id": revision_id}))
        return {"revision_id": revision_id, "enabled": False}

    def list(self):
        self.calls.append(("list", {}))
        return ()

    def enable(
        self,
        profile_id: str,
        *,
        revision_id: str | None = None,
        security_audit_actor=None,
    ):
        self.calls.append(
            (
                "enable",
                {
                    "profile_id": profile_id,
                    "revision_id": revision_id,
                    "security_audit_actor": security_audit_actor,
                },
            )
        )
        return {"profile_id": profile_id, "revision_id": revision_id, "enabled": True}

    def disable(self, profile_id: str, *, security_audit_actor=None):
        self.calls.append(
            (
                "disable",
                {
                    "profile_id": profile_id,
                    "security_audit_actor": security_audit_actor,
                },
            )
        )
        return {"profile_id": profile_id, "enabled": False}


def _factory(registry, observed):
    @contextmanager
    def factory(
        database_url: str,
        config_path: Path | None,
        allow_schema_upgrade: bool,
    ):
        observed.append((database_url, config_path, allow_schema_upgrade))
        yield registry

    return factory


def test_reference_profile_register_uses_explicit_private_config_and_safe_output(
    tmp_path,
    capsys,
) -> None:
    config_path = tmp_path / "operator-references.yaml"
    registry = RecordingReferenceProfiles()
    observed = []

    exit_code = admin.main(
        [
            "--database-url",
            DATABASE_URL,
            "--reference-profile-config",
            str(config_path),
            "reference-profile",
            "register",
            "--safe-key",
            "grch38",
            "--display-name",
            "GRCh38 primary",
            "--organism",
            "Homo sapiens",
            "--assembly",
            "GRCh38",
            "--config-key",
            "grch38-private",
        ],
        reference_profile_factory=_factory(registry, observed),
    )

    assert exit_code == 0
    assert observed == [(DATABASE_URL, config_path, True)]
    assert registry.calls == [
        (
            "register",
            {
                "safe_key": "grch38",
                "display_name": "GRCh38 primary",
                "organism": "Homo sapiens",
                "assembly": "GRCh38",
                "config_key": "grch38-private",
                "security_audit_actor": AuthenticationActor.local_operator(),
            },
        )
    ]
    output = capsys.readouterr().out
    assert str(config_path) not in output
    assert "grch38-private" not in output


def test_reference_profile_commands_have_exact_mutation_surface(
    tmp_path, capsys
) -> None:
    config_path = tmp_path / "operator-references.json"
    registry = RecordingReferenceProfiles()
    observed = []
    factory = _factory(registry, observed)
    profile_id = "refp_" + "1" * 32
    revision_id = "refpr_" + "2" * 32

    commands = (
        (["reference-profile", "list"], ("list", {})),
        (
            [
                "--reference-profile-config",
                str(config_path),
                "reference-profile",
                "verify",
                revision_id,
            ],
            ("verify", {"revision_id": revision_id}),
        ),
        (
            [
                "--reference-profile-config",
                str(config_path),
                "reference-profile",
                "enable",
                profile_id,
                "--revision-id",
                revision_id,
            ],
            (
                "enable",
                {
                    "profile_id": profile_id,
                    "revision_id": revision_id,
                    "security_audit_actor": AuthenticationActor.local_operator(),
                },
            ),
        ),
        (
            ["reference-profile", "disable", profile_id],
            (
                "disable",
                {
                    "profile_id": profile_id,
                    "security_audit_actor": AuthenticationActor.local_operator(),
                },
            ),
        ),
    )
    for command, expected in commands:
        assert (
            admin.main(
                ["--database-url", DATABASE_URL, *command],
                reference_profile_factory=factory,
            )
            == 0
        )
        assert registry.calls[-1] == expected
    assert [item[2] for item in observed] == [False, False, True, True]
    capsys.readouterr()


def test_reference_profile_failure_uses_stable_redacted_reason_code(
    tmp_path,
    capsys,
) -> None:
    secret = tmp_path / "operator" / "references.json"

    class Failing(RecordingReferenceProfiles):
        def verify(self, revision_id: str):
            raise ReferenceProfileAdminError("REFERENCE_PROFILE_CONFIG_INVALID")

    exit_code = admin.main(
        [
            "--database-url",
            DATABASE_URL,
            "--reference-profile-config",
            str(secret),
            "reference-profile",
            "verify",
            "refpr_" + "2" * 32,
        ],
        reference_profile_factory=_factory(Failing(), []),
    )

    assert exit_code == 1
    error = json.loads(capsys.readouterr().err)["error"]
    assert error["code"] == "REFERENCE_PROFILE_CONFIG_INVALID"
    assert str(secret) not in json.dumps(error)


def test_reference_profile_missing_private_config_is_stable_json(capsys) -> None:
    exit_code = admin.main(
        [
            "--database-url",
            DATABASE_URL,
            "reference-profile",
            "verify",
            "refpr_" + "2" * 32,
        ],
        reference_profile_factory=_factory(RecordingReferenceProfiles(), []),
    )

    assert exit_code == 1
    error = json.loads(capsys.readouterr().err)["error"]
    assert error == {
        "code": "REFERENCE_PROFILE_CONFIG_REQUIRED",
        "message": "Private Reference Profile configuration is required.",
    }


def test_reference_profile_unexpected_failure_does_not_disclose_path(
    tmp_path,
    capsys,
) -> None:
    secret = tmp_path / "operator" / "references.json"

    @contextmanager
    def failing_factory(
        _database_url: str,
        _config_path: Path | None,
        _allow_schema_upgrade: bool,
    ):
        raise OSError(f"could not read {secret}")
        yield RecordingReferenceProfiles()

    exit_code = admin.main(
        [
            "--database-url",
            DATABASE_URL,
            "--reference-profile-config",
            str(secret),
            "reference-profile",
            "verify",
            "refpr_" + "2" * 32,
        ],
        reference_profile_factory=failing_factory,
    )

    assert exit_code == 1
    error = json.loads(capsys.readouterr().err)["error"]
    assert error == {
        "code": "REFERENCE_PROFILE_COMMAND_FAILED",
        "message": "Reference Profile command failed.",
    }
    assert str(secret) not in json.dumps(error)


def test_reference_profile_arbitrary_factory_failure_is_redacted(
    tmp_path,
    capsys,
) -> None:
    secret = tmp_path / "operator" / "database.db"

    @contextmanager
    def failing_factory(
        _database_url: str,
        _config_path: Path | None,
        _allow_schema_upgrade: bool,
    ):
        raise Exception(f"driver leaked {secret}")
        yield RecordingReferenceProfiles()

    exit_code = admin.main(
        ["--database-url", DATABASE_URL, "reference-profile", "list"],
        reference_profile_factory=failing_factory,
    )

    captured = capsys.readouterr()
    assert exit_code == 1
    assert captured.out == ""
    assert json.loads(captured.err)["error"] == {
        "code": "REFERENCE_PROFILE_COMMAND_FAILED",
        "message": "Reference Profile command failed.",
    }
    assert str(secret) not in captured.err
    assert "Traceback" not in captured.err


def test_reference_profile_missing_database_parent_is_redacted_and_read_only(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    database_path = tmp_path / "missing" / "platform.db"

    exit_code = admin.main(
        [
            "--database-url",
            f"sqlite:///{database_path}",
            "reference-profile",
            "list",
        ]
    )

    captured = capsys.readouterr()
    assert exit_code == 1
    assert captured.out == ""
    assert json.loads(captured.err)["error"] == {
        "code": "REFERENCE_PROFILE_COMMAND_FAILED",
        "message": "Reference Profile command failed.",
    }
    assert str(database_path) not in captured.err
    assert "Traceback" not in captured.err
    assert not database_path.parent.exists()


@pytest.mark.parametrize("command", ["list", "verify"])
def test_reference_profile_read_commands_never_mutate_current_database(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    command: str,
) -> None:
    database_path = tmp_path / "platform.db"
    database_url = f"sqlite:///{database_path}"
    upgrade_database(database_url)
    before_digest = hashlib.sha256(database_path.read_bytes()).hexdigest()
    before_mtime_ns = database_path.stat().st_mtime_ns
    before_entries = tuple(sorted(path.name for path in tmp_path.iterdir()))
    arguments = ["--database-url", database_url]
    if command == "verify":
        arguments.extend(
            [
                "--reference-profile-config",
                str(tmp_path / "operator-references.json"),
            ]
        )
    arguments.extend(["reference-profile", command])
    if command == "verify":
        arguments.append("refpr_" + "2" * 32)

    exit_code = admin.main(arguments)

    captured = capsys.readouterr()
    if command == "list":
        assert exit_code == 0
        assert json.loads(captured.out) == []
    else:
        assert exit_code == 1
        assert json.loads(captured.err)["error"]["code"] == (
            "REFERENCE_PROFILE_UNAVAILABLE"
        )
    assert hashlib.sha256(database_path.read_bytes()).hexdigest() == before_digest
    assert database_path.stat().st_mtime_ns == before_mtime_ns
    assert tuple(sorted(path.name for path in tmp_path.iterdir())) == before_entries


def test_reference_profile_read_command_refuses_prior_schema_without_upgrade(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    database_path = tmp_path / "platform.db"
    database_url = f"sqlite:///{database_path}"
    upgrade_database(database_url, revision="20260726_10")
    before_digest = hashlib.sha256(database_path.read_bytes()).hexdigest()

    exit_code = admin.main(
        ["--database-url", database_url, "reference-profile", "list"]
    )

    error = json.loads(capsys.readouterr().err)["error"]
    assert exit_code == 1
    assert error == {
        "code": "REFERENCE_PROFILE_COMMAND_FAILED",
        "message": "Reference Profile command failed.",
    }
    assert hashlib.sha256(database_path.read_bytes()).hexdigest() == before_digest
    with sqlite3.connect(f"{database_path.as_uri()}?mode=ro", uri=True) as database:
        assert database.execute(
            "SELECT version_num FROM alembic_version"
        ).fetchone() == ("20260726_10",)
