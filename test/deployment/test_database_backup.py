from __future__ import annotations

import inspect
import os
from pathlib import Path
import sqlite3

import pytest

from encode_pipeline.deployment.canonical import (
    canonical_identity,
    canonical_json_bytes,
)
import encode_pipeline.deployment.database as database_module
from encode_pipeline.deployment.database import (
    DATABASE_WRITER_UNITS,
    WRITE_STOP_ISSUER,
    WRITE_STOP_WITNESS_IDENTITY_SCHEME,
    WRITE_STOP_WITNESS_SCHEMA,
    DatabaseSchemaAdmission,
    WriteStopWitness,
    backup_database,
    database_path_identity,
    inspect_database,
    write_stop_witness_path,
)
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.layout import DeploymentLayout


DEPLOYMENT_IDENTITY = f"sha256-{'a' * 64}"
STATE_IDENTITY = f"sha256-{'b' * 64}"
PROVIDER_IDENTITY = f"sha256-{'c' * 64}"
CONSUME_ONCE_IDENTITY = f"sha256-{'d' * 64}"
SERVICE_IDENTITIES = (f"sha256-{'e' * 64}", f"sha256-{'f' * 64}")
TASK_IDENTITY = f"task-{'1' * 32}"
CREATED_AT = "2026-08-09T12:00:00+00:00"
SCHEMA_HEAD = "native_resolver_head_42"


class AdmissionProvider:
    def __init__(self, heads: tuple[str, ...] = (SCHEMA_HEAD,)) -> None:
        self.heads = heads

    def admit(self, **values) -> DatabaseSchemaAdmission:
        return DatabaseSchemaAdmission.create(
            provider_identity=PROVIDER_IDENTITY,
            accepted_schema_heads=self.heads,
            **values,
        )


def _database(layout: DeploymentLayout, *, head: str = SCHEMA_HEAD) -> Path:
    layout.database.parent.mkdir(parents=True, mode=0o700)
    connection = sqlite3.connect(layout.database)
    try:
        connection.execute("CREATE TABLE alembic_version (version_num TEXT NOT NULL)")
        connection.execute(
            "INSERT INTO alembic_version (version_num) VALUES (?)", (head,)
        )
        connection.execute(
            "CREATE TABLE durable_state (id INTEGER PRIMARY KEY, value TEXT NOT NULL)"
        )
        connection.execute("INSERT INTO durable_state (value) VALUES ('preserved')")
        connection.commit()
    finally:
        connection.close()
    return layout.database


def _ready_layout(tmp_path: Path, *, head: str = SCHEMA_HEAD) -> DeploymentLayout:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    _database(layout, head=head)
    return layout


def _witness_document(
    layout: DeploymentLayout,
    *,
    database_device: int | None = None,
    database_inode: int | None = None,
) -> dict[str, object]:
    observed = layout.database.stat()
    value: dict[str, object] = {
        "schema_version": WRITE_STOP_WITNESS_SCHEMA,
        "issuer": WRITE_STOP_ISSUER,
        "consumption_state": "ready",
        "task_identity": TASK_IDENTITY,
        "deployment_identity": DEPLOYMENT_IDENTITY,
        "operator_state_identity": STATE_IDENTITY,
        "consume_once_identity": CONSUME_ONCE_IDENTITY,
        "database_path_identity": database_path_identity(layout.database),
        "database_device": observed.st_dev
        if database_device is None
        else database_device,
        "database_inode": observed.st_ino if database_inode is None else database_inode,
        "created_at": CREATED_AT,
        "writers": [
            {
                "unit": unit,
                "state": "stopped",
                "service_identity": SERVICE_IDENTITIES[index],
            }
            for index, unit in enumerate(DATABASE_WRITER_UNITS)
        ],
    }
    value["identity"] = canonical_identity(
        value, scheme=WRITE_STOP_WITNESS_IDENTITY_SCHEME
    )
    return value


def _write_witness(
    layout: DeploymentLayout,
    document: dict[str, object] | None = None,
) -> Path:
    path = write_stop_witness_path(layout, TASK_IDENTITY)
    path.parent.mkdir(parents=True, mode=0o700)
    path.parent.chmod(0o700)
    path.write_bytes(canonical_json_bytes(document or _witness_document(layout)))
    path.chmod(0o600)
    return path


def _allow_test_owner(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(database_module, "_ROOT_UID", os.getuid())
    monkeypatch.setattr(database_module, "_ROOT_GID", os.getgid())


def _backup(layout: DeploymentLayout, provider=None) -> None:
    backup_database(
        layout,
        task_identity=TASK_IDENTITY,
        expected_deployment_identity=DEPLOYMENT_IDENTITY,
        expected_operator_state_identity=STATE_IDENTITY,
        schema_provider=AdmissionProvider() if provider is None else provider,
    )


def test_inspection_uses_one_pinned_immutable_fd_and_dynamic_schema_head(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    layout = _ready_layout(tmp_path, head="arbitrary_native_head_9000")
    original_connect = sqlite3.connect
    observed_uris: list[str] = []

    def recording_connect(database, *args, **kwargs):
        observed_uris.append(str(database))
        return original_connect(database, *args, **kwargs)

    monkeypatch.setattr(database_module.sqlite3, "connect", recording_connect)
    inspection = inspect_database(layout.database)

    assert inspection.schema_heads == ("arbitrary_native_head_9000",)
    assert inspection.wal_state == "absent-or-empty"
    assert len(observed_uris) == 1
    assert observed_uris[0].startswith("file:/proc/self/fd/")
    assert "mode=ro&immutable=1" in observed_uris[0]


@pytest.mark.parametrize("suffix", ("-wal", "-shm"))
def test_nonempty_wal_or_shm_is_rejected_without_touching_it(
    tmp_path: Path, suffix: str
) -> None:
    layout = _ready_layout(tmp_path)
    sidecar = Path(f"{layout.database}{suffix}")
    sidecar.write_bytes(b"active writer evidence")

    with pytest.raises(DeploymentError) as captured:
        inspect_database(layout.database)

    assert captured.value.issue.code == "DATABASE_WAL_NOT_QUIESCENT"
    assert sidecar.read_bytes() == b"active writer evidence"


def test_backup_has_no_caller_minted_witness_or_schema_head_parameters() -> None:
    parameters = inspect.signature(backup_database).parameters

    assert "witness" not in parameters
    assert "accepted_schema_heads" not in parameters
    assert "platform_manifest" not in parameters
    assert "schema_provider" in parameters


def test_default_schema_provider_fails_closed_before_any_backup(
    tmp_path: Path,
) -> None:
    layout = _ready_layout(tmp_path)

    with pytest.raises(DeploymentError) as captured:
        backup_database(
            layout,
            task_identity=TASK_IDENTITY,
            expected_deployment_identity=DEPLOYMENT_IDENTITY,
            expected_operator_state_identity=STATE_IDENTITY,
        )

    assert captured.value.issue.code == "DATABASE_SCHEMA_ADMISSION_UNAVAILABLE"
    assert not layout.database_backups.exists()


def test_schema_provider_errors_are_redacted(tmp_path: Path) -> None:
    layout = _ready_layout(tmp_path)

    class LeakyProvider:
        def admit(self, **values):
            del values
            raise fail("PRIVATE_PROVIDER_ERROR", "/private/path token=secret")

    with pytest.raises(DeploymentError) as captured:
        _backup(layout, LeakyProvider())

    assert captured.value.issue.code == "DATABASE_SCHEMA_ADMISSION_FAILED"
    assert "private" not in str(captured.value)
    assert "secret" not in str(captured.value)


def test_valid_root_file_prerequisites_still_require_atomic_operator_consumption(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    layout = _ready_layout(tmp_path)
    _allow_test_owner(monkeypatch)
    _write_witness(layout)

    for _attempt in range(2):
        with pytest.raises(DeploymentError) as captured:
            _backup(layout)
        assert captured.value.issue.code == "DATABASE_BACKUP_OPERATOR_CONSUME_REQUIRED"

    assert not layout.database_backups.exists()
    assert write_stop_witness_path(layout, TASK_IDENTITY).exists()


def test_witness_requires_canonical_bytes_owner_mode_and_single_link(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    layout = _ready_layout(tmp_path)
    _allow_test_owner(monkeypatch)
    path = _write_witness(layout)

    path.chmod(0o644)
    with pytest.raises(DeploymentError) as wrong_mode:
        _backup(layout)
    assert wrong_mode.value.issue.code == "DATABASE_WRITE_STOP_WITNESS_UNTRUSTED"

    path.chmod(0o600)
    path.write_bytes(path.read_bytes() + b" ")
    with pytest.raises(DeploymentError) as noncanonical:
        _backup(layout)
    assert noncanonical.value.issue.code == "DATABASE_WRITE_STOP_WITNESS_INVALID"

    path.write_bytes(canonical_json_bytes(_witness_document(layout)))
    hardlink = path.with_name("second-link.json")
    os.link(path, hardlink)
    with pytest.raises(DeploymentError) as linked:
        _backup(layout)
    assert linked.value.issue.code == "DATABASE_WRITE_STOP_WITNESS_UNTRUSTED"


def test_untrusted_owner_and_symlink_witness_fail_closed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    layout = _ready_layout(tmp_path)
    path = _write_witness(layout)
    monkeypatch.setattr(database_module, "_ROOT_UID", os.getuid() + 1)
    monkeypatch.setattr(database_module, "_ROOT_GID", os.getgid())

    with pytest.raises(DeploymentError) as owner:
        _backup(layout)
    assert owner.value.issue.code == "DATABASE_WRITE_STOP_WITNESS_UNTRUSTED"

    monkeypatch.setattr(database_module, "_ROOT_UID", os.getuid())
    path.unlink()
    target = path.with_name("target.json")
    target.write_bytes(canonical_json_bytes(_witness_document(layout)))
    target.chmod(0o600)
    path.symlink_to(target)
    with pytest.raises(DeploymentError) as symlink:
        _backup(layout)
    assert symlink.value.issue.code == "DATABASE_WRITE_STOP_WITNESS_INVALID"


def test_witness_binds_database_device_inode_state_and_both_stopped_services(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    layout = _ready_layout(tmp_path)
    _allow_test_owner(monkeypatch)
    wrong_inode = _witness_document(
        layout, database_inode=layout.database.stat().st_ino + 1
    )
    _write_witness(layout, wrong_inode)

    with pytest.raises(DeploymentError) as mismatch:
        _backup(layout)
    assert mismatch.value.issue.code == "DATABASE_WRITE_STOP_WITNESS_MISMATCH"

    incomplete = _witness_document(layout)
    incomplete["writers"] = incomplete["writers"][:1]
    incomplete["identity"] = canonical_identity(
        {key: value for key, value in incomplete.items() if key != "identity"},
        scheme=WRITE_STOP_WITNESS_IDENTITY_SCHEME,
    )
    with pytest.raises(DeploymentError) as missing_writer:
        WriteStopWitness.from_dict(incomplete)
    assert missing_writer.value.issue.code == "DATABASE_WRITE_STOP_WITNESS_INVALID"


def test_schema_admission_is_identity_bound_to_operation_and_database(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    layout = _ready_layout(tmp_path)
    _allow_test_owner(monkeypatch)
    _write_witness(layout)

    class WrongDatabaseProvider(AdmissionProvider):
        def admit(self, **values) -> DatabaseSchemaAdmission:
            values["database_inode"] += 1
            return super().admit(**values)

    with pytest.raises(DeploymentError) as mismatch:
        _backup(layout, WrongDatabaseProvider())
    assert mismatch.value.issue.code == "DATABASE_SCHEMA_ADMISSION_MISMATCH"

    with pytest.raises(DeploymentError) as wrong_head:
        _backup(layout, AdmissionProvider(("different_native_head",)))
    assert wrong_head.value.issue.code == "DATABASE_SCHEMA_INCOMPATIBLE"
    assert not layout.database_backups.exists()


def test_database_path_symlink_and_writable_boundaries_are_rejected(
    tmp_path: Path,
) -> None:
    layout = _ready_layout(tmp_path)
    real = layout.database.with_name("real.db")
    layout.database.rename(real)
    layout.database.symlink_to(real)
    with pytest.raises(DeploymentError) as symlink:
        inspect_database(layout.database)
    assert symlink.value.issue.code == "DATABASE_PATH_INVALID"

    layout.database.unlink()
    real.rename(layout.database)
    layout.database.parent.chmod(0o777)
    with pytest.raises(DeploymentError) as parent:
        inspect_database(layout.database)
    assert parent.value.issue.code == "DATABASE_PATH_UNSAFE"

    layout.database.parent.chmod(0o700)
    layout.database.chmod(0o666)
    with pytest.raises(DeploymentError) as database:
        inspect_database(layout.database)
    assert database.value.issue.code == "DATABASE_INVALID"


@pytest.mark.parametrize("owner_index", (4, 5), ids=("uid", "gid"))
def test_database_owner_mismatch_fails_before_sqlite_or_backup_output(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    owner_index: int,
) -> None:
    layout = _ready_layout(tmp_path)
    expected_database = layout.database.stat()
    original_fstat = os.fstat
    sqlite_entered = False

    def mismatched_database_owner(descriptor: int) -> os.stat_result:
        observed = original_fstat(descriptor)
        if (
            observed.st_dev == expected_database.st_dev
            and observed.st_ino == expected_database.st_ino
        ):
            values = list(observed)
            values[owner_index] += 1
            return os.stat_result(values)
        return observed

    def forbidden_connect(*args, **kwargs):
        nonlocal sqlite_entered
        sqlite_entered = True
        raise AssertionError("SQLite must not be entered for an untrusted database")

    monkeypatch.setattr(database_module.os, "fstat", mismatched_database_owner)
    monkeypatch.setattr(database_module.sqlite3, "connect", forbidden_connect)

    with pytest.raises(DeploymentError) as captured:
        _backup(layout)

    assert captured.value.issue.code == "DATABASE_INVALID"
    assert sqlite_entered is False
    assert not layout.database_backups.exists()
