from __future__ import annotations

import inspect
import os
from pathlib import Path
import sqlite3
import stat

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
    DatabaseBackupReceipt,
    InvalidFreshDatabaseQuarantineReceipt,
    DatabaseQuarantineReceipt,
    DatabaseRestoreReceipt,
    WriteStopWitness,
    backup_database,
    consumed_write_stop_witness_path,
    database_content_identity,
    database_path_identity,
    fresh_database_candidate_path,
    inspect_database,
    publish_fresh_database,
    quarantine_invalid_fresh_database,
    quarantine_fresh_database,
    restore_database_backup,
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
FAILED_SCHEMA_HEAD = "failed_migration_head_43"


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
    return _database_at(layout.database, head=head, value="preserved")


def _database_at(path: Path, *, head: str, value: str) -> Path:
    path.parent.mkdir(parents=True, mode=0o700, exist_ok=True)
    connection = sqlite3.connect(path)
    try:
        connection.execute("CREATE TABLE alembic_version (version_num TEXT NOT NULL)")
        connection.execute(
            "INSERT INTO alembic_version (version_num) VALUES (?)", (head,)
        )
        connection.execute(
            "CREATE TABLE durable_state (id INTEGER PRIMARY KEY, value TEXT NOT NULL)"
        )
        connection.execute("INSERT INTO durable_state (value) VALUES (?)", (value,))
        connection.commit()
    finally:
        connection.close()
    return path


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


def _shared_live_boundary(layout: DeploymentLayout) -> None:
    layout.database.parent.chmod(0o2770)
    if layout.database.exists():
        layout.database.chmod(0o660)


def _fresh_candidate(
    layout: DeploymentLayout,
    *,
    head: str = SCHEMA_HEAD,
    value: str = "fresh",
) -> Path:
    path = fresh_database_candidate_path(layout, TASK_IDENTITY)
    _database_at(path, head=head, value=value)
    path.chmod(0o660)
    return path


def _fresh_layout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> DeploymentLayout:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    layout.database.parent.mkdir(parents=True, mode=0o2770)
    layout.database.parent.chmod(0o2770)
    _allow_test_owner(monkeypatch)
    return layout


def _database_value(path: Path) -> str:
    connection = sqlite3.connect(f"file:{path}?mode=ro", uri=True)
    try:
        return connection.execute("SELECT value FROM durable_state").fetchone()[0]
    finally:
        connection.close()


def _restore_fixture(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[DeploymentLayout, DatabaseBackupReceipt]:
    layout = _ready_layout(tmp_path)
    _shared_live_boundary(layout)
    _allow_test_owner(monkeypatch)
    _write_witness(layout)
    backup = _backup(layout)
    connection = sqlite3.connect(layout.database)
    try:
        connection.execute("UPDATE durable_state SET value = 'failed-current'")
        connection.execute(
            "UPDATE alembic_version SET version_num = ?", (FAILED_SCHEMA_HEAD,)
        )
        connection.commit()
    finally:
        connection.close()
    layout.database.chmod(0o660)
    return layout, backup


def _restore(
    layout: DeploymentLayout,
    backup: DatabaseBackupReceipt,
    *,
    fault=None,
) -> DatabaseRestoreReceipt:
    return restore_database_backup(
        layout,
        backup_identity=backup.backup_identity,
        expected_task_identity=TASK_IDENTITY,
        expected_deployment_identity=DEPLOYMENT_IDENTITY,
        expected_prior_state_identity=STATE_IDENTITY,
        expected_source_identity=database_content_identity(backup.source),
        expected_schema_heads=(SCHEMA_HEAD,),
        expected_database_uid=os.getuid(),
        expected_database_gid=os.getgid(),
        fault=fault,
    )


def _invalid_fresh_quarantine(
    layout: DeploymentLayout, *, fault=None
) -> InvalidFreshDatabaseQuarantineReceipt:
    return quarantine_invalid_fresh_database(
        layout,
        task_identity=TASK_IDENTITY,
        expected_deployment_identity=DEPLOYMENT_IDENTITY,
        expected_prior_state_identity=STATE_IDENTITY,
        expected_database_uid=os.getuid(),
        expected_database_gid=os.getgid(),
        fault=fault,
    )


def _backup(
    layout: DeploymentLayout, provider=None, *, fault=None
) -> DatabaseBackupReceipt:
    return backup_database(
        layout,
        task_identity=TASK_IDENTITY,
        expected_deployment_identity=DEPLOYMENT_IDENTITY,
        expected_operator_state_identity=STATE_IDENTITY,
        expected_database_uid=os.getuid(),
        expected_database_gid=os.getgid(),
        schema_provider=AdmissionProvider() if provider is None else provider,
        fault=fault,
    )


def _inspect_database(path: Path):
    return inspect_database(
        path,
        expected_owner_uid=os.getuid(),
        expected_owner_gid=os.getgid(),
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
    inspection = _inspect_database(layout.database)

    assert inspection.schema_heads == ("arbitrary_native_head_9000",)
    assert inspection.wal_state == "absent-or-empty"
    assert len(observed_uris) == 1
    assert observed_uris[0].startswith("file:/proc/self/fd/")
    assert "mode=ro&immutable=1" in observed_uris[0]


@pytest.mark.parametrize("suffix", ("-wal", "-shm", "-journal"))
def test_nonempty_sqlite_sidecar_is_rejected_without_touching_it(
    tmp_path: Path, suffix: str
) -> None:
    layout = _ready_layout(tmp_path)
    sidecar = Path(f"{layout.database}{suffix}")
    sidecar.write_bytes(b"active writer evidence")

    with pytest.raises(DeploymentError) as captured:
        _inspect_database(layout.database)

    assert captured.value.issue.code == "DATABASE_SIDECAR_NOT_QUIESCENT"
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
            expected_database_uid=os.getuid(),
            expected_database_gid=os.getgid(),
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


def test_valid_root_prerequisites_create_one_idempotent_immutable_backup(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    layout = _ready_layout(tmp_path)
    _allow_test_owner(monkeypatch)
    _write_witness(layout)

    first = _backup(layout)
    second = _backup(layout)

    assert first == second
    assert first.backup.sha256 == first.source.sha256
    assert not write_stop_witness_path(layout, TASK_IDENTITY).exists()
    assert consumed_write_stop_witness_path(layout, TASK_IDENTITY).exists()
    slot = layout.database_backups / first.backup_identity
    assert stat.S_IMODE(slot.stat().st_mode) == 0o500
    assert stat.S_IMODE((slot / "platform.db").stat().st_mode) == 0o400
    assert stat.S_IMODE((slot / "receipt.json").stat().st_mode) == 0o444
    connection = sqlite3.connect(f"file:{slot / 'platform.db'}?mode=ro", uri=True)
    try:
        assert connection.execute("SELECT value FROM durable_state").fetchone() == (
            "preserved",
        )
    finally:
        connection.close()


def test_backup_resumes_after_witness_consumption_without_replay(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    layout = _ready_layout(tmp_path)
    _allow_test_owner(monkeypatch)
    _write_witness(layout)

    def interrupt(phase: str) -> None:
        if phase == "witness-consumed":
            raise RuntimeError("private interruption")

    with pytest.raises(DeploymentError) as interrupted:
        _backup(layout, fault=interrupt)
    assert interrupted.value.issue.code == "DATABASE_BACKUP_FAILED"
    assert not write_stop_witness_path(layout, TASK_IDENTITY).exists()
    assert consumed_write_stop_witness_path(layout, TASK_IDENTITY).exists()

    receipt = _backup(layout)
    assert receipt.backup_identity.startswith("sha256-")
    assert len(tuple(layout.database_backups.glob("sha256-*"))) == 1


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
        _inspect_database(layout.database)
    assert symlink.value.issue.code == "DATABASE_PATH_INVALID"

    layout.database.unlink()
    real.rename(layout.database)
    layout.database.parent.chmod(0o777)
    with pytest.raises(DeploymentError) as parent:
        _inspect_database(layout.database)
    assert parent.value.issue.code == "DATABASE_PATH_UNSAFE"

    layout.database.parent.chmod(0o700)
    layout.database.chmod(0o666)
    with pytest.raises(DeploymentError) as database:
        _inspect_database(layout.database)
    assert database.value.issue.code == "DATABASE_INVALID"


def test_exact_setgid_shared_database_boundary_is_identity_bound(
    tmp_path: Path,
) -> None:
    layout = _ready_layout(tmp_path)
    layout.database.parent.chmod(0o2770)
    layout.database.chmod(0o660)
    sidecar = Path(f"{layout.database}-wal")
    sidecar.write_bytes(b"")
    sidecar.chmod(0o660)

    inspection = _inspect_database(layout.database)

    assert inspection.schema_heads == (SCHEMA_HEAD,)
    with pytest.raises(DeploymentError) as wrong_group:
        inspect_database(
            layout.database,
            expected_owner_uid=os.getuid(),
            expected_owner_gid=os.getgid() + 1,
        )
    assert wrong_group.value.issue.code == "DATABASE_PATH_UNSAFE"

    layout.database.parent.chmod(0o770)
    with pytest.raises(DeploymentError) as missing_setgid:
        _inspect_database(layout.database)
    assert missing_setgid.value.issue.code == "DATABASE_PATH_UNSAFE"


def test_shared_database_requires_exact_file_and_sidecar_modes(tmp_path: Path) -> None:
    layout = _ready_layout(tmp_path)
    layout.database.parent.chmod(0o2770)
    layout.database.chmod(0o640)
    with pytest.raises(DeploymentError) as database:
        _inspect_database(layout.database)
    assert database.value.issue.code == "DATABASE_INVALID"

    layout.database.chmod(0o660)
    sidecar = Path(f"{layout.database}-shm")
    sidecar.write_bytes(b"")
    sidecar.chmod(0o640)
    with pytest.raises(DeploymentError) as sidecar_mode:
        _inspect_database(layout.database)
    assert sidecar_mode.value.issue.code == "DATABASE_SIDECAR_NOT_QUIESCENT"


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


def test_recovery_primitives_expose_only_fixed_database_coordinates() -> None:
    publish_parameters = inspect.signature(publish_fresh_database).parameters
    restore_parameters = inspect.signature(restore_database_backup).parameters
    quarantine_parameters = inspect.signature(quarantine_fresh_database).parameters
    invalid_quarantine_parameters = inspect.signature(
        quarantine_invalid_fresh_database
    ).parameters

    assert "candidate_path" not in publish_parameters
    assert "database_path" not in publish_parameters
    assert "backup_path" not in restore_parameters
    assert "evidence_path" not in restore_parameters
    assert "source_path" not in quarantine_parameters
    assert not {
        "source_path",
        "candidate_path",
        "expected_candidate_identity",
        "expected_schema_heads",
    } & set(invalid_quarantine_parameters)

    layout = DeploymentLayout.isolated(Path("/tmp/fixed-database-coordinate"))
    with pytest.raises(DeploymentError) as captured:
        fresh_database_candidate_path(layout, "../../private")
    assert captured.value.issue.code == "DATABASE_FRESH_PUBLISH_REQUEST_INVALID"


def test_operator_parent_allows_dedicated_group_with_exact_traverse_mode(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    operator = tmp_path / "operator"
    operator.mkdir(mode=0o710)
    operator.chmod(0o710)
    monkeypatch.setattr(database_module, "_ROOT_UID", os.getuid())
    monkeypatch.setattr(database_module, "_ROOT_GID", os.getgid() + 1)

    database_module._require_or_create_operator_root_directory(operator, create=False)

    assert operator.stat().st_gid != database_module._ROOT_GID


@pytest.mark.parametrize("mode", (0o730, 0o711))
def test_operator_parent_rejects_write_or_world_traverse_modes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mode: int,
) -> None:
    operator = tmp_path / "operator"
    operator.mkdir(mode=mode)
    operator.chmod(mode)
    monkeypatch.setattr(database_module, "_ROOT_UID", os.getuid())

    with pytest.raises(DeploymentError) as captured:
        database_module._require_or_create_operator_root_directory(
            operator, create=False
        )

    assert captured.value.issue.code == "DATABASE_BACKUP_BOUNDARY_UNTRUSTED"


def test_fresh_candidate_is_atomically_published_and_retry_is_idempotent(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout = _fresh_layout(tmp_path, monkeypatch)
    candidate = _fresh_candidate(layout)
    candidate_stat = candidate.stat()
    candidate_inspection = _inspect_database(candidate)
    candidate_identity = database_content_identity(candidate_inspection)

    published = publish_fresh_database(
        layout,
        task_identity=TASK_IDENTITY,
        expected_candidate_identity=candidate_identity,
        target_schema_heads=(SCHEMA_HEAD,),
        expected_database_uid=os.getuid(),
        expected_database_gid=os.getgid(),
    )
    retried = publish_fresh_database(
        layout,
        task_identity=TASK_IDENTITY,
        expected_candidate_identity=candidate_identity,
        target_schema_heads=(SCHEMA_HEAD,),
        expected_database_uid=os.getuid(),
        expected_database_gid=os.getgid(),
    )

    assert published == retried
    assert database_content_identity(published) == candidate_identity
    assert not candidate.exists()
    assert layout.database.stat().st_ino == candidate_stat.st_ino
    assert stat.S_IMODE(layout.database.stat().st_mode) == 0o660
    assert _database_value(layout.database) == "fresh"


def test_fresh_publish_never_overwrites_a_racing_canonical_database(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout = _fresh_layout(tmp_path, monkeypatch)
    candidate = _fresh_candidate(layout)
    candidate_identity = database_content_identity(_inspect_database(candidate))
    original_rename = database_module._rename_noreplace_at

    def race_with_canonical(*args, **kwargs) -> None:
        _database_at(layout.database, head=SCHEMA_HEAD, value="race-winner")
        layout.database.chmod(0o660)
        original_rename(*args, **kwargs)

    monkeypatch.setattr(database_module, "_rename_noreplace_at", race_with_canonical)

    with pytest.raises(DeploymentError) as captured:
        publish_fresh_database(
            layout,
            task_identity=TASK_IDENTITY,
            expected_candidate_identity=candidate_identity,
            target_schema_heads=(SCHEMA_HEAD,),
            expected_database_uid=os.getuid(),
            expected_database_gid=os.getgid(),
        )

    assert captured.value.issue.code == "DATABASE_FRESH_PUBLISH_INVALID"
    assert candidate.exists()
    assert _database_value(layout.database) == "race-winner"


@pytest.mark.parametrize(
    "phase", ("fresh-candidate-synced", "fresh-database-published")
)
def test_fresh_publish_recovers_each_fsync_and_rename_window(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    phase: str,
) -> None:
    layout = _fresh_layout(tmp_path, monkeypatch)
    candidate = _fresh_candidate(layout)
    candidate_identity = database_content_identity(_inspect_database(candidate))

    def interrupt(observed: str) -> None:
        if observed == phase:
            raise RuntimeError("private interruption")

    with pytest.raises(DeploymentError) as captured:
        publish_fresh_database(
            layout,
            task_identity=TASK_IDENTITY,
            expected_candidate_identity=candidate_identity,
            target_schema_heads=(SCHEMA_HEAD,),
            expected_database_uid=os.getuid(),
            expected_database_gid=os.getgid(),
            fault=interrupt,
        )
    assert captured.value.issue.code == "DATABASE_FRESH_PUBLISH_FAILED"
    assert "private" not in str(captured.value)

    published = publish_fresh_database(
        layout,
        task_identity=TASK_IDENTITY,
        expected_candidate_identity=candidate_identity,
        target_schema_heads=(SCHEMA_HEAD,),
        expected_database_uid=os.getuid(),
        expected_database_gid=os.getgid(),
    )
    assert database_content_identity(published) == candidate_identity


@pytest.mark.parametrize(
    "fault", ("identity", "head", "sidecar", "hardlink", "symlink")
)
def test_fresh_publish_rejects_untrusted_candidate_and_sidecar_state(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fault: str,
) -> None:
    layout = _fresh_layout(tmp_path, monkeypatch)
    candidate = _fresh_candidate(layout)
    identity = database_content_identity(_inspect_database(candidate))
    heads = (SCHEMA_HEAD,)
    if fault == "identity":
        identity = f"sha256-{'9' * 64}"
    elif fault == "head":
        heads = (FAILED_SCHEMA_HEAD,)
    elif fault == "sidecar":
        sidecar = Path(f"{candidate}-wal")
        sidecar.write_bytes(b"writer")
        sidecar.chmod(0o660)
    elif fault == "hardlink":
        os.link(candidate, candidate.with_name("candidate-hardlink.db"))
    else:
        real = candidate.with_name("candidate-real.db")
        candidate.rename(real)
        candidate.symlink_to(real)

    with pytest.raises(DeploymentError):
        publish_fresh_database(
            layout,
            task_identity=TASK_IDENTITY,
            expected_candidate_identity=identity,
            target_schema_heads=heads,
            expected_database_uid=os.getuid(),
            expected_database_gid=os.getgid(),
        )

    assert not layout.database.exists()


def test_restore_revalidates_backup_preserves_failed_current_and_is_idempotent(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout, backup = _restore_fixture(tmp_path, monkeypatch)
    failed_bytes = layout.database.read_bytes()

    first = _restore(layout, backup)
    second = _restore(layout, backup)

    assert first == second
    assert isinstance(first, DatabaseRestoreReceipt)
    assert first.backup_identity == backup.backup_identity
    assert first.backup_receipt_identity == backup.identity
    assert _database_value(layout.database) == "preserved"
    assert _inspect_database(layout.database).schema_heads == (SCHEMA_HEAD,)
    evidence = (
        layout.data_root
        / "operator"
        / "database-recovery"
        / "restore"
        / TASK_IDENTITY
        / backup.backup_identity
    )
    receipt_bytes = (evidence / "receipt.json").read_bytes()
    assert (evidence / "platform.db").read_bytes() == failed_bytes
    assert stat.S_IMODE(evidence.stat().st_mode) == 0o500
    assert stat.S_IMODE((evidence / "platform.db").stat().st_mode) == 0o400
    assert stat.S_IMODE((evidence / "receipt.json").stat().st_mode) == 0o444
    assert canonical_json_bytes(first.to_dict()) == receipt_bytes
    assert str(tmp_path).encode() not in receipt_bytes


@pytest.mark.parametrize(
    "phase",
    (
        "restore-evidence-file-synced",
        "restore-evidence-ready",
        "restore-evidence-committed",
        "restore-copy-synced",
        "restore-ready-to-replace",
        "restore-live-replaced",
    ),
)
def test_restore_resumes_every_evidence_copy_fsync_and_replace_window(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    phase: str,
) -> None:
    layout, backup = _restore_fixture(tmp_path, monkeypatch)

    def interrupt(observed: str) -> None:
        if observed == phase:
            raise RuntimeError("/private/database token=secret")

    with pytest.raises(DeploymentError) as captured:
        _restore(layout, backup, fault=interrupt)
    assert captured.value.issue.code == "DATABASE_RESTORE_FAILED"
    assert "private" not in str(captured.value)
    assert "secret" not in str(captured.value)

    receipt = _restore(layout, backup)
    assert receipt.backup_identity == backup.backup_identity
    assert _database_value(layout.database) == "preserved"


@pytest.mark.parametrize("failure", ("copy", "rename", "directory-fsync"))
def test_restore_handles_injected_copy_rename_and_directory_fsync_failures(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    failure: str,
) -> None:
    layout, backup = _restore_fixture(tmp_path, monkeypatch)
    original_copy = database_module._copy_database_descriptor_at
    original_replace = database_module.os.replace
    original_fsync = database_module.os.fsync
    replaced = False
    fsync_failed = False

    if failure == "copy":

        def fail_restore_copy(source, directory_descriptor, name, **kwargs):
            if name.startswith(".restore-"):
                raise OSError("private copy failure")
            return original_copy(source, directory_descriptor, name, **kwargs)

        monkeypatch.setattr(
            database_module, "_copy_database_descriptor_at", fail_restore_copy
        )
    elif failure == "rename":

        def fail_replace(*args, **kwargs):
            raise OSError("private rename failure")

        monkeypatch.setattr(database_module.os, "replace", fail_replace)
    else:

        def observe_replace(*args, **kwargs):
            nonlocal replaced
            result = original_replace(*args, **kwargs)
            replaced = True
            return result

        def fail_post_replace_fsync(descriptor: int):
            nonlocal fsync_failed
            if replaced and not fsync_failed:
                fsync_failed = True
                raise OSError("private fsync failure")
            return original_fsync(descriptor)

        monkeypatch.setattr(database_module.os, "replace", observe_replace)
        monkeypatch.setattr(database_module.os, "fsync", fail_post_replace_fsync)

    with pytest.raises(DeploymentError) as captured:
        _restore(layout, backup)
    assert captured.value.issue.code == "DATABASE_RESTORE_FAILED"

    monkeypatch.setattr(database_module, "_copy_database_descriptor_at", original_copy)
    monkeypatch.setattr(database_module.os, "replace", original_replace)
    monkeypatch.setattr(database_module.os, "fsync", original_fsync)
    receipt = _restore(layout, backup)
    assert receipt.backup_identity == backup.backup_identity
    assert _database_value(layout.database) == "preserved"


@pytest.mark.parametrize(
    ("field", "value"),
    (
        ("expected_task_identity", f"task-{'2' * 32}"),
        ("expected_deployment_identity", f"sha256-{'2' * 64}"),
        ("expected_prior_state_identity", f"sha256-{'3' * 64}"),
        ("expected_source_identity", f"sha256-{'4' * 64}"),
        ("expected_schema_heads", ("wrong_head",)),
    ),
)
def test_restore_rejects_backup_binding_mismatches_without_touching_live(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    field: str,
    value,
) -> None:
    layout, backup = _restore_fixture(tmp_path, monkeypatch)
    before = layout.database.read_bytes()
    arguments = {
        "backup_identity": backup.backup_identity,
        "expected_task_identity": TASK_IDENTITY,
        "expected_deployment_identity": DEPLOYMENT_IDENTITY,
        "expected_prior_state_identity": STATE_IDENTITY,
        "expected_source_identity": database_content_identity(backup.source),
        "expected_schema_heads": (SCHEMA_HEAD,),
        "expected_database_uid": os.getuid(),
        "expected_database_gid": os.getgid(),
    }
    arguments[field] = value

    with pytest.raises(DeploymentError) as captured:
        restore_database_backup(layout, **arguments)

    assert captured.value.issue.code == "DATABASE_BACKUP_INVALID"
    assert layout.database.read_bytes() == before
    assert _database_value(layout.database) == "failed-current"


def test_restore_rejects_live_sidecars_and_tampered_backup_bytes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout, backup = _restore_fixture(tmp_path, monkeypatch)
    sidecar = Path(f"{layout.database}-wal")
    sidecar.write_bytes(b"writer")
    sidecar.chmod(0o660)

    with pytest.raises(DeploymentError) as sidecar_error:
        _restore(layout, backup)
    assert sidecar_error.value.issue.code == "DATABASE_SIDECAR_NOT_QUIESCENT"
    sidecar.unlink()

    slot = layout.database_backups / backup.backup_identity
    slot.chmod(0o700)
    backup_path = slot / "platform.db"
    backup_path.chmod(0o600)
    backup_path.write_bytes(backup_path.read_bytes() + b"tamper")
    backup_path.chmod(0o400)
    slot.chmod(0o500)
    with pytest.raises(DeploymentError) as tampered:
        _restore(layout, backup)
    assert tampered.value.issue.code in {
        "DATABASE_BACKUP_INVALID",
        "DATABASE_INVALID",
    }
    assert _database_value(layout.database) == "failed-current"


@pytest.mark.parametrize("fault", ("hardlink", "symlink"))
def test_restore_rejects_untrusted_backup_slot_inode_boundaries(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fault: str,
) -> None:
    layout, backup = _restore_fixture(tmp_path, monkeypatch)
    slot = layout.database_backups / backup.backup_identity
    slot.chmod(0o700)
    backup_path = slot / "platform.db"
    if fault == "hardlink":
        os.link(backup_path, slot / "second-link.db")
    else:
        real = slot / "real.db"
        backup_path.rename(real)
        backup_path.symlink_to(real)
    slot.chmod(0o500)

    with pytest.raises(DeploymentError) as captured:
        _restore(layout, backup)

    assert captured.value.issue.code in {
        "DATABASE_BACKUP_INVALID",
        "DATABASE_INVALID",
        "DATABASE_PATH_INVALID",
    }
    assert _database_value(layout.database) == "failed-current"


@pytest.mark.parametrize("published", (False, True), ids=("candidate", "canonical"))
def test_fresh_quarantine_preserves_evidence_and_never_touches_user_data(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    published: bool,
) -> None:
    layout = _fresh_layout(tmp_path, monkeypatch)
    candidate = _fresh_candidate(layout)
    candidate_identity = database_content_identity(_inspect_database(candidate))
    if published:
        publish_fresh_database(
            layout,
            task_identity=TASK_IDENTITY,
            expected_candidate_identity=candidate_identity,
            target_schema_heads=(SCHEMA_HEAD,),
            expected_database_uid=os.getuid(),
            expected_database_gid=os.getgid(),
        )
    sentinels = []
    for path in (
        layout.workspaces / "keep.txt",
        layout.artifacts / "keep.txt",
        layout.configuration_root / "reference-profiles.yaml",
    ):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"operator data")
        sentinels.append(path)

    first = quarantine_fresh_database(
        layout,
        task_identity=TASK_IDENTITY,
        expected_deployment_identity=DEPLOYMENT_IDENTITY,
        expected_prior_state_identity=STATE_IDENTITY,
        expected_candidate_identity=candidate_identity,
        expected_schema_heads=(SCHEMA_HEAD,),
        expected_database_uid=os.getuid(),
        expected_database_gid=os.getgid(),
    )
    second = quarantine_fresh_database(
        layout,
        task_identity=TASK_IDENTITY,
        expected_deployment_identity=DEPLOYMENT_IDENTITY,
        expected_prior_state_identity=STATE_IDENTITY,
        expected_candidate_identity=candidate_identity,
        expected_schema_heads=(SCHEMA_HEAD,),
        expected_database_uid=os.getuid(),
        expected_database_gid=os.getgid(),
    )

    assert first == second
    assert isinstance(first, DatabaseQuarantineReceipt)
    assert first.source_coordinate == ("canonical" if published else "candidate")
    assert not candidate.exists()
    assert not layout.database.exists()
    assert all(path.read_bytes() == b"operator data" for path in sentinels)
    evidence = (
        layout.data_root
        / "operator"
        / "database-recovery"
        / "fresh-quarantine"
        / TASK_IDENTITY
        / candidate_identity
    )
    assert stat.S_IMODE(evidence.stat().st_mode) == 0o500
    assert stat.S_IMODE((evidence / "platform.db").stat().st_mode) == 0o400
    receipt_bytes = (evidence / "receipt.json").read_bytes()
    assert canonical_json_bytes(first.to_dict()) == receipt_bytes
    assert str(tmp_path).encode() not in receipt_bytes


@pytest.mark.parametrize(
    "phase",
    (
        "quarantine-evidence-file-synced",
        "quarantine-evidence-ready",
        "quarantine-evidence-committed",
        "fresh-database-quarantined",
    ),
)
def test_fresh_quarantine_resumes_every_evidence_and_unlink_window(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    phase: str,
) -> None:
    layout = _fresh_layout(tmp_path, monkeypatch)
    candidate = _fresh_candidate(layout)
    candidate_identity = database_content_identity(_inspect_database(candidate))

    def interrupt(observed: str) -> None:
        if observed == phase:
            raise RuntimeError("private interruption")

    arguments = {
        "task_identity": TASK_IDENTITY,
        "expected_deployment_identity": DEPLOYMENT_IDENTITY,
        "expected_prior_state_identity": STATE_IDENTITY,
        "expected_candidate_identity": candidate_identity,
        "expected_schema_heads": (SCHEMA_HEAD,),
        "expected_database_uid": os.getuid(),
        "expected_database_gid": os.getgid(),
    }
    with pytest.raises(DeploymentError) as captured:
        quarantine_fresh_database(layout, **arguments, fault=interrupt)
    assert captured.value.issue.code == "DATABASE_QUARANTINE_FAILED"

    receipt = quarantine_fresh_database(layout, **arguments)
    assert receipt.candidate_identity == candidate_identity
    assert not candidate.exists()


@pytest.mark.parametrize(
    "fault", ("identity", "head", "sidecar", "hardlink", "symlink")
)
def test_fresh_quarantine_fails_closed_without_deleting_untrusted_source(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fault: str,
) -> None:
    layout = _fresh_layout(tmp_path, monkeypatch)
    candidate = _fresh_candidate(layout)
    identity = database_content_identity(_inspect_database(candidate))
    heads = (SCHEMA_HEAD,)
    if fault == "identity":
        identity = f"sha256-{'8' * 64}"
    elif fault == "head":
        heads = (FAILED_SCHEMA_HEAD,)
    elif fault == "sidecar":
        sidecar = Path(f"{candidate}-journal")
        sidecar.write_bytes(b"writer")
        sidecar.chmod(0o660)
    elif fault == "hardlink":
        os.link(candidate, candidate.with_name("candidate-hardlink.db"))
    else:
        real = candidate.with_name("candidate-real.db")
        candidate.rename(real)
        candidate.symlink_to(real)

    with pytest.raises(DeploymentError):
        quarantine_fresh_database(
            layout,
            task_identity=TASK_IDENTITY,
            expected_deployment_identity=DEPLOYMENT_IDENTITY,
            expected_prior_state_identity=STATE_IDENTITY,
            expected_candidate_identity=identity,
            expected_schema_heads=heads,
            expected_database_uid=os.getuid(),
            expected_database_gid=os.getgid(),
        )

    assert candidate.exists()


@pytest.mark.parametrize("candidate_kind", ("invalid", "headless"))
def test_invalid_fresh_quarantine_preserves_raw_candidate_and_hot_sidecars(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    candidate_kind: str,
) -> None:
    layout = _fresh_layout(tmp_path, monkeypatch)
    _database_at(layout.database, head=SCHEMA_HEAD, value="canonical")
    layout.database.chmod(0o660)
    canonical_bytes = layout.database.read_bytes()
    canonical_inode = layout.database.stat().st_ino
    candidate = fresh_database_candidate_path(layout, TASK_IDENTITY)
    if candidate_kind == "invalid":
        candidate.write_bytes(b"not a sqlite database")
    else:
        connection = sqlite3.connect(candidate)
        try:
            connection.execute("CREATE TABLE unfinished (value TEXT)")
            connection.commit()
        finally:
            connection.close()
    candidate.chmod(0o660)
    source_bytes = {"database": candidate.read_bytes()}
    for role, suffix, content in (
        ("journal", "-journal", b"hot journal"),
        ("shm", "-shm", b""),
        ("wal", "-wal", b"hot wal"),
    ):
        sidecar = Path(f"{candidate}{suffix}")
        sidecar.write_bytes(content)
        sidecar.chmod(0o660)
        source_bytes[role] = content

    first = _invalid_fresh_quarantine(layout)
    second = _invalid_fresh_quarantine(layout)

    assert first == second
    assert tuple(item.role for item in first.files) == (
        "database",
        "journal",
        "shm",
        "wal",
    )
    receipt_bytes = canonical_json_bytes(first.to_dict())
    assert str(tmp_path).encode() not in receipt_bytes
    assert b"schema_heads" not in receipt_bytes
    assert b"candidate_identity" not in receipt_bytes
    evidence = (
        layout.data_root
        / "operator"
        / "database-recovery"
        / "fresh-invalid"
        / TASK_IDENTITY
        / first.request_identity
    )
    assert stat.S_IMODE(evidence.stat().st_mode) == 0o500
    for item in first.files:
        preserved = evidence / database_module._INVALID_FRESH_FILE_NAMES[item.role]
        assert preserved.read_bytes() == source_bytes[item.role]
        assert stat.S_IMODE(preserved.stat().st_mode) == 0o400
    assert stat.S_IMODE((evidence / "receipt.json").stat().st_mode) == 0o444
    assert not candidate.exists()
    assert all(
        not Path(f"{candidate}{suffix}").exists()
        for suffix in ("-journal", "-shm", "-wal")
    )
    assert layout.database.read_bytes() == canonical_bytes
    assert layout.database.stat().st_ino == canonical_inode
    assert _database_value(layout.database) == "canonical"


@pytest.mark.parametrize(
    "phase",
    (
        "invalid-fresh-evidence-file-synced:database",
        "invalid-fresh-evidence-file-synced:wal",
        "invalid-fresh-evidence-receipt-synced",
        "invalid-fresh-evidence-directory-synced",
        "invalid-fresh-evidence-committed",
        "invalid-fresh-source-removed:database",
        "invalid-fresh-source-removed:wal",
        "invalid-fresh-source-directory-synced",
    ),
)
def test_invalid_fresh_quarantine_resumes_every_evidence_and_source_window(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    phase: str,
) -> None:
    layout = _fresh_layout(tmp_path, monkeypatch)
    candidate = fresh_database_candidate_path(layout, TASK_IDENTITY)
    candidate.write_bytes(b"partial failed migration")
    candidate.chmod(0o660)
    wal = Path(f"{candidate}-wal")
    wal.write_bytes(b"uncheckpointed bytes")
    wal.chmod(0o660)

    def interrupt(observed: str) -> None:
        if observed == phase:
            raise RuntimeError("private interruption")

    with pytest.raises(DeploymentError) as captured:
        _invalid_fresh_quarantine(layout, fault=interrupt)
    assert captured.value.issue.code == "DATABASE_INVALID_FRESH_QUARANTINE_FAILED"
    assert "private" not in str(captured.value)

    receipt = _invalid_fresh_quarantine(layout)
    assert tuple(item.role for item in receipt.files) == ("database", "wal")
    assert not candidate.exists()
    assert not wal.exists()


@pytest.mark.parametrize(
    "fault",
    (
        "main-symlink",
        "main-hardlink",
        "main-mode",
        "main-owner",
        "main-fifo",
        "sidecar-symlink",
        "sidecar-hardlink",
        "sidecar-mode",
        "sidecar-owner",
        "sidecar-fifo",
    ),
)
def test_invalid_fresh_quarantine_rejects_untrusted_main_and_sidecars(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fault: str,
) -> None:
    layout = _fresh_layout(tmp_path, monkeypatch)
    _database_at(layout.database, head=SCHEMA_HEAD, value="canonical")
    layout.database.chmod(0o660)
    canonical_bytes = layout.database.read_bytes()
    candidate = fresh_database_candidate_path(layout, TASK_IDENTITY)
    candidate.write_bytes(b"invalid candidate")
    candidate.chmod(0o660)
    selected = candidate
    if fault.startswith("sidecar-"):
        selected = Path(f"{candidate}-wal")
        selected.write_bytes(b"hot sidecar")
        selected.chmod(0o660)
    kind = fault.removeprefix("main-").removeprefix("sidecar-")
    if kind == "symlink":
        real = selected.with_name(f"{selected.name}.real")
        selected.rename(real)
        selected.symlink_to(real.name)
    elif kind == "hardlink":
        os.link(selected, selected.with_name(f"{selected.name}.link"))
    elif kind == "mode":
        selected.chmod(0o640)
    elif kind == "fifo":
        selected.unlink()
        os.mkfifo(selected, mode=0o660)
    else:
        expected = selected.stat()
        original_fstat = os.fstat

        def mismatched_owner(descriptor: int) -> os.stat_result:
            observed = original_fstat(descriptor)
            if (observed.st_dev, observed.st_ino) == (
                expected.st_dev,
                expected.st_ino,
            ):
                values = list(observed)
                values[4] += 1
                return os.stat_result(values)
            return observed

        monkeypatch.setattr(database_module.os, "fstat", mismatched_owner)

    with pytest.raises(DeploymentError) as captured:
        _invalid_fresh_quarantine(layout)

    assert captured.value.issue.code == (
        "DATABASE_INVALID_FRESH_QUARANTINE_SOURCE_INVALID"
    )
    assert layout.database.read_bytes() == canonical_bytes
    assert os.path.lexists(candidate)


def test_invalid_fresh_quarantine_requires_fixed_candidate_or_verified_receipt(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout = _fresh_layout(tmp_path, monkeypatch)
    _database_at(layout.database, head=SCHEMA_HEAD, value="canonical")
    layout.database.chmod(0o660)
    canonical_bytes = layout.database.read_bytes()

    with pytest.raises(DeploymentError) as captured:
        _invalid_fresh_quarantine(layout)

    assert captured.value.issue.code == (
        "DATABASE_INVALID_FRESH_QUARANTINE_SOURCE_MISSING"
    )
    assert layout.database.read_bytes() == canonical_bytes


def test_invalid_fresh_quarantine_rejects_tampered_retained_evidence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout = _fresh_layout(tmp_path, monkeypatch)
    candidate = fresh_database_candidate_path(layout, TASK_IDENTITY)
    candidate.write_bytes(b"failed candidate")
    candidate.chmod(0o660)
    receipt = _invalid_fresh_quarantine(layout)
    evidence = (
        layout.data_root
        / "operator"
        / "database-recovery"
        / "fresh-invalid"
        / TASK_IDENTITY
        / receipt.request_identity
    )
    retained = evidence / "database.raw"
    evidence.chmod(0o700)
    retained.chmod(0o600)
    retained.write_bytes(b"tampered")
    retained.chmod(0o400)
    evidence.chmod(0o500)

    with pytest.raises(DeploymentError) as captured:
        _invalid_fresh_quarantine(layout)

    assert captured.value.issue.code == (
        "DATABASE_INVALID_FRESH_QUARANTINE_EVIDENCE_INVALID"
    )
