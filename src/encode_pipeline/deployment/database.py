"""Read-only SQLite checks and the Phase-A backup authorization boundary.

The supported backup must atomically consume a root-operator write-stop
witness.  That consuming backend is intentionally deferred until PR #172 is
merged, so :func:`backup_database` validates every read-only prerequisite and
then fails closed without creating a backup.  This module has no migration,
restore, downgrade, deletion, or caller-minted authorization API.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import stat
from typing import Any, NoReturn, Protocol
from urllib.parse import quote

from encode_pipeline.deployment.canonical import (
    canonical_identity,
    canonical_json_bytes,
    without_key,
)
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.layout import DeploymentLayout


WRITE_STOP_WITNESS_SCHEMA = "helixweave-database-write-stop-witness-v2"
DATABASE_INSPECTION_SCHEMA = "helixweave-database-inspection-v1"
DATABASE_SCHEMA_ADMISSION_SCHEMA = "helixweave-database-schema-admission-v1"
WRITE_STOP_WITNESS_IDENTITY_SCHEME = (
    "helixweave-database-write-stop-witness-identity-v2"
)
DATABASE_PATH_IDENTITY_SCHEME = "helixweave-database-path-identity-v1"
DATABASE_SCHEMA_ADMISSION_IDENTITY_SCHEME = (
    "helixweave-database-schema-admission-identity-v1"
)
WRITE_STOP_ISSUER = "helixweave-root-operator-v1"
DATABASE_WRITER_UNITS = (
    "helixweave-api.service",
    "helixweave-worker.service",
)

_IDENTITY = re.compile(r"^sha256-[0-9a-f]{64}$")
_TASK_IDENTITY = re.compile(r"^task-[0-9a-f]{32}$")
_SCHEMA_HEAD = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_ROOT_UID = 0
_ROOT_GID = 0
_WITNESS_MODE = 0o600
_WITNESS_DIRECTORY_MODE = 0o700
_MAX_WITNESS_BYTES = 64 * 1024
_MAX_DATABASE_BYTES = 16 * 1024**4
_READ_CHUNK_BYTES = 1024 * 1024


def _document(value: object, *, code: str) -> dict[str, Any]:
    if not isinstance(value, Mapping) or any(not isinstance(key, str) for key in value):
        raise fail(code, "Database evidence is invalid.")
    return dict(value)


def _content_identity(value: object, *, code: str) -> str:
    if not isinstance(value, str) or _IDENTITY.fullmatch(value) is None:
        raise fail(code, "Database evidence is invalid.")
    return value


def _task_identity(value: object, *, code: str) -> str:
    if not isinstance(value, str) or _TASK_IDENTITY.fullmatch(value) is None:
        raise fail(code, "Database evidence is invalid.")
    return value


def _positive_integer(value: object, *, code: str) -> int:
    if (
        not isinstance(value, int)
        or isinstance(value, bool)
        or not 0 < value <= 2**63 - 1
    ):
        raise fail(code, "Database evidence is invalid.")
    return value


def _timestamp(value: object, *, code: str) -> str:
    if not isinstance(value, str) or len(value) > 64:
        raise fail(code, "Database evidence is invalid.")
    try:
        parsed = datetime.fromisoformat(value)
    except ValueError:
        raise fail(code, "Database evidence is invalid.") from None
    if parsed.tzinfo is None or parsed.utcoffset() is None:
        raise fail(code, "Database evidence is invalid.")
    return value


@dataclass(frozen=True, order=True)
class StoppedWriter:
    """Exact stopped-service evidence emitted by the privileged backend."""

    unit: str
    state: str
    service_identity: str

    @classmethod
    def from_dict(cls, raw: object) -> "StoppedWriter":
        code = "DATABASE_WRITE_STOP_WITNESS_INVALID"
        value = _document(raw, code=code)
        if (
            set(value) != {"unit", "state", "service_identity"}
            or value["unit"] not in DATABASE_WRITER_UNITS
            or value["state"] != "stopped"
        ):
            raise fail(code, "Database write-stop witness is invalid.")
        return cls(
            unit=value["unit"],
            state="stopped",
            service_identity=_content_identity(value["service_identity"], code=code),
        )

    def to_dict(self) -> dict[str, str]:
        return {
            "unit": self.unit,
            "state": self.state,
            "service_identity": self.service_identity,
        }


@dataclass(frozen=True)
class WriteStopWitness:
    """Parsed evidence; never accepted directly as backup authorization."""

    identity: str
    task_identity: str
    deployment_identity: str
    operator_state_identity: str
    consume_once_identity: str
    database_path_identity: str
    database_device: int
    database_inode: int
    created_at: str
    writers: tuple[StoppedWriter, ...]

    @classmethod
    def from_dict(cls, raw: object) -> "WriteStopWitness":
        code = "DATABASE_WRITE_STOP_WITNESS_INVALID"
        value = _document(raw, code=code)
        if set(value) != {
            "schema_version",
            "identity",
            "issuer",
            "consumption_state",
            "task_identity",
            "deployment_identity",
            "operator_state_identity",
            "consume_once_identity",
            "database_path_identity",
            "database_device",
            "database_inode",
            "created_at",
            "writers",
        }:
            raise fail(code, "Database write-stop witness is invalid.")
        raw_writers = value["writers"]
        if (
            value["schema_version"] != WRITE_STOP_WITNESS_SCHEMA
            or value["issuer"] != WRITE_STOP_ISSUER
            or value["consumption_state"] != "ready"
            or not isinstance(raw_writers, list)
        ):
            raise fail(code, "Database write-stop witness is invalid.")
        writers = tuple(StoppedWriter.from_dict(item) for item in raw_writers)
        if tuple(item.unit for item in writers) != DATABASE_WRITER_UNITS or len(
            {item.service_identity for item in writers}
        ) != len(writers):
            raise fail(code, "Database write-stop witness is invalid.")
        observed_identity = _content_identity(value["identity"], code=code)
        task_identity = _task_identity(value["task_identity"], code=code)
        deployment_identity = _content_identity(value["deployment_identity"], code=code)
        operator_state_identity = _content_identity(
            value["operator_state_identity"], code=code
        )
        consume_once_identity = _content_identity(
            value["consume_once_identity"], code=code
        )
        path_identity = _content_identity(value["database_path_identity"], code=code)
        database_device = _positive_integer(value["database_device"], code=code)
        database_inode = _positive_integer(value["database_inode"], code=code)
        created_at = _timestamp(value["created_at"], code=code)
        expected_identity = canonical_identity(
            without_key(value, "identity"),
            scheme=WRITE_STOP_WITNESS_IDENTITY_SCHEME,
        )
        if observed_identity != expected_identity:
            raise fail(code, "Database write-stop witness is invalid.")
        return cls(
            identity=observed_identity,
            task_identity=task_identity,
            deployment_identity=deployment_identity,
            operator_state_identity=operator_state_identity,
            consume_once_identity=consume_once_identity,
            database_path_identity=path_identity,
            database_device=database_device,
            database_inode=database_inode,
            created_at=created_at,
            writers=writers,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": WRITE_STOP_WITNESS_SCHEMA,
            "identity": self.identity,
            "issuer": WRITE_STOP_ISSUER,
            "consumption_state": "ready",
            "task_identity": self.task_identity,
            "deployment_identity": self.deployment_identity,
            "operator_state_identity": self.operator_state_identity,
            "consume_once_identity": self.consume_once_identity,
            "database_path_identity": self.database_path_identity,
            "database_device": self.database_device,
            "database_inode": self.database_inode,
            "created_at": self.created_at,
            "writers": [writer.to_dict() for writer in self.writers],
        }


@dataclass(frozen=True)
class DatabaseInspection:
    schema_heads: tuple[str, ...]
    size_bytes: int
    sha256: str
    device: int
    inode: int
    wal_state: str = "absent-or-empty"

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": DATABASE_INSPECTION_SCHEMA,
            "integrity": "ok",
            "schema_heads": list(self.schema_heads),
            "size_bytes": self.size_bytes,
            "sha256": self.sha256,
            "device": self.device,
            "inode": self.inode,
            "wal_state": self.wal_state,
        }


@dataclass(frozen=True)
class DatabaseSchemaAdmission:
    """Identity-bound output from the post-merge native schema resolver."""

    identity: str
    provider_identity: str
    task_identity: str
    deployment_identity: str
    operator_state_identity: str
    database_device: int
    database_inode: int
    accepted_schema_heads: tuple[str, ...]

    @classmethod
    def create(
        cls,
        *,
        provider_identity: str,
        task_identity: str,
        deployment_identity: str,
        operator_state_identity: str,
        database_device: int,
        database_inode: int,
        accepted_schema_heads: Sequence[str],
    ) -> "DatabaseSchemaAdmission":
        value: dict[str, object] = {
            "schema_version": DATABASE_SCHEMA_ADMISSION_SCHEMA,
            "provider_identity": provider_identity,
            "task_identity": task_identity,
            "deployment_identity": deployment_identity,
            "operator_state_identity": operator_state_identity,
            "database_device": database_device,
            "database_inode": database_inode,
            "accepted_schema_heads": list(accepted_schema_heads),
        }
        value["identity"] = canonical_identity(
            value, scheme=DATABASE_SCHEMA_ADMISSION_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "DatabaseSchemaAdmission":
        code = "DATABASE_SCHEMA_ADMISSION_INVALID"
        value = _document(raw, code=code)
        if (
            set(value)
            != {
                "schema_version",
                "identity",
                "provider_identity",
                "task_identity",
                "deployment_identity",
                "operator_state_identity",
                "database_device",
                "database_inode",
                "accepted_schema_heads",
            }
            or value["schema_version"] != DATABASE_SCHEMA_ADMISSION_SCHEMA
        ):
            raise fail(code, "Database schema admission is invalid.")
        observed = _content_identity(value["identity"], code=code)
        provider = _content_identity(value["provider_identity"], code=code)
        task = _task_identity(value["task_identity"], code=code)
        deployment = _content_identity(value["deployment_identity"], code=code)
        operator_state = _content_identity(value["operator_state_identity"], code=code)
        device = _positive_integer(value["database_device"], code=code)
        inode = _positive_integer(value["database_inode"], code=code)
        heads = _schema_heads(value["accepted_schema_heads"], code=code)
        if observed != canonical_identity(
            without_key(value, "identity"),
            scheme=DATABASE_SCHEMA_ADMISSION_IDENTITY_SCHEME,
        ):
            raise fail(code, "Database schema admission is invalid.")
        return cls(
            observed,
            provider,
            task,
            deployment,
            operator_state,
            device,
            inode,
            heads,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": DATABASE_SCHEMA_ADMISSION_SCHEMA,
            "identity": self.identity,
            "provider_identity": self.provider_identity,
            "task_identity": self.task_identity,
            "deployment_identity": self.deployment_identity,
            "operator_state_identity": self.operator_state_identity,
            "database_device": self.database_device,
            "database_inode": self.database_inode,
            "accepted_schema_heads": list(self.accepted_schema_heads),
        }


class SchemaAdmissionProvider(Protocol):
    def admit(
        self,
        *,
        task_identity: str,
        deployment_identity: str,
        operator_state_identity: str,
        database_device: int,
        database_inode: int,
    ) -> DatabaseSchemaAdmission: ...


class UnavailableSchemaAdmissionProvider:
    def admit(
        self,
        *,
        task_identity: str,
        deployment_identity: str,
        operator_state_identity: str,
        database_device: int,
        database_inode: int,
    ) -> DatabaseSchemaAdmission:
        del (
            task_identity,
            deployment_identity,
            operator_state_identity,
            database_device,
            database_inode,
        )
        raise fail(
            "DATABASE_SCHEMA_ADMISSION_UNAVAILABLE",
            "Trusted database schema admission is not installed.",
            recoverable=True,
        )


@dataclass
class _OpenedDatabase:
    path: Path
    directory_descriptor: int
    descriptor: int
    initial_stat: os.stat_result

    def close(self) -> None:
        os.close(self.descriptor)
        os.close(self.directory_descriptor)


def database_path_identity(path: Path) -> str:
    selected = _absolute_path(path, existing=True, code="DATABASE_PATH_INVALID")
    return canonical_identity(
        {"absolute_path": str(selected)}, scheme=DATABASE_PATH_IDENTITY_SCHEME
    )


def write_stop_witness_path(layout: DeploymentLayout, task_identity: str) -> Path:
    task = _task_identity(task_identity, code="DATABASE_BACKUP_REQUEST_INVALID")
    return layout.run_root / "operator" / "database-backups" / task / "write-stop.json"


def inspect_database(path: Path) -> DatabaseInspection:
    """Inspect one pinned database FD and reject active WAL/SHM sidecars."""
    opened = _open_database(path)
    try:
        _require_quiescent_sidecars(opened)
        inspection = _inspect_descriptor(opened)
        _require_quiescent_sidecars(opened)
        _require_unchanged(opened)
        return inspection
    except DeploymentError:
        raise
    except (OSError, sqlite3.DatabaseError):
        raise fail("DATABASE_INVALID", "Database is invalid.") from None
    finally:
        opened.close()


def backup_database(
    layout: DeploymentLayout,
    *,
    task_identity: str,
    expected_deployment_identity: str,
    expected_operator_state_identity: str,
    schema_provider: SchemaAdmissionProvider | None = None,
) -> NoReturn:
    """Validate backup prerequisites, then fail until atomic consumption exists.

    A content-valid witness can still be replayed.  Only the post-PR-172 root
    backend can atomically transition its ``consume_once_identity`` from ready
    to consumed while creating the backup, so Phase A never creates output.
    """
    if not isinstance(layout, DeploymentLayout):
        raise fail("DATABASE_BACKUP_REQUEST_INVALID", "Database backup is invalid.")
    task = _task_identity(task_identity, code="DATABASE_BACKUP_REQUEST_INVALID")
    deployment = _content_identity(
        expected_deployment_identity, code="DATABASE_BACKUP_REQUEST_INVALID"
    )
    operator_state = _content_identity(
        expected_operator_state_identity, code="DATABASE_BACKUP_REQUEST_INVALID"
    )
    opened = _open_database(layout.database)
    try:
        _require_quiescent_sidecars(opened)
        selected_provider = (
            UnavailableSchemaAdmissionProvider()
            if schema_provider is None
            else schema_provider
        )
        try:
            admission = DatabaseSchemaAdmission.from_dict(
                selected_provider.admit(
                    task_identity=task,
                    deployment_identity=deployment,
                    operator_state_identity=operator_state,
                    database_device=opened.initial_stat.st_dev,
                    database_inode=opened.initial_stat.st_ino,
                ).to_dict()
            )
        except DeploymentError as error:
            if error.issue.code == "DATABASE_SCHEMA_ADMISSION_UNAVAILABLE":
                raise
            raise fail(
                "DATABASE_SCHEMA_ADMISSION_FAILED",
                "Trusted database schema admission failed.",
            ) from None
        except Exception:
            raise fail(
                "DATABASE_SCHEMA_ADMISSION_FAILED",
                "Trusted database schema admission failed.",
            ) from None
        if (
            admission.task_identity != task
            or admission.deployment_identity != deployment
            or admission.operator_state_identity != operator_state
            or admission.database_device != opened.initial_stat.st_dev
            or admission.database_inode != opened.initial_stat.st_ino
        ):
            raise fail(
                "DATABASE_SCHEMA_ADMISSION_MISMATCH",
                "Database schema admission does not match this operation.",
            )
        witness = _read_write_stop_witness(
            write_stop_witness_path(layout, task),
            expected_task_identity=task,
            expected_deployment_identity=deployment,
            expected_operator_state_identity=operator_state,
            expected_database_path_identity=database_path_identity(layout.database),
            expected_database_stat=opened.initial_stat,
        )
        inspection = _inspect_descriptor(opened)
        if inspection.schema_heads != admission.accepted_schema_heads:
            raise fail(
                "DATABASE_SCHEMA_INCOMPATIBLE",
                "Database schema is not compatible with the verified platform contract.",
            )
        _require_quiescent_sidecars(opened)
        _require_unchanged(opened)
        if witness.consume_once_identity == witness.identity:
            raise fail(
                "DATABASE_WRITE_STOP_WITNESS_INVALID",
                "Database write-stop witness is invalid.",
            )
        raise fail(
            "DATABASE_BACKUP_OPERATOR_CONSUME_REQUIRED",
            "Database backup requires atomic root-operator witness consumption.",
            recoverable=True,
        )
    finally:
        opened.close()


def _read_write_stop_witness(
    path: Path,
    *,
    expected_task_identity: str,
    expected_deployment_identity: str,
    expected_operator_state_identity: str,
    expected_database_path_identity: str,
    expected_database_stat: os.stat_result,
) -> WriteStopWitness:
    selected = _absolute_path(
        path, existing=True, code="DATABASE_WRITE_STOP_WITNESS_INVALID"
    )
    parent_stat = selected.parent.lstat()
    if (
        not stat.S_ISDIR(parent_stat.st_mode)
        or stat.S_ISLNK(parent_stat.st_mode)
        or parent_stat.st_uid != _ROOT_UID
        or parent_stat.st_gid != _ROOT_GID
        or stat.S_IMODE(parent_stat.st_mode) != _WITNESS_DIRECTORY_MODE
    ):
        raise fail(
            "DATABASE_WRITE_STOP_WITNESS_UNTRUSTED",
            "Database write-stop witness is not trusted.",
        )
    directory_flags = (
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    file_flags = (
        os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    )
    directory_descriptor: int | None = None
    descriptor: int | None = None
    try:
        directory_descriptor = os.open(selected.parent, directory_flags)
        anchored_parent = os.fstat(directory_descriptor)
        if _inode(anchored_parent) != _inode(parent_stat):
            raise fail(
                "DATABASE_WRITE_STOP_WITNESS_UNTRUSTED",
                "Database write-stop witness is not trusted.",
            )
        descriptor = os.open(selected.name, file_flags, dir_fd=directory_descriptor)
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_uid != _ROOT_UID
            or before.st_gid != _ROOT_GID
            or stat.S_IMODE(before.st_mode) != _WITNESS_MODE
            or not 0 < before.st_size <= _MAX_WITNESS_BYTES
        ):
            raise fail(
                "DATABASE_WRITE_STOP_WITNESS_UNTRUSTED",
                "Database write-stop witness is not trusted.",
            )
        content = _read_bounded(descriptor, _MAX_WITNESS_BYTES)
        after = os.fstat(descriptor)
        path_after = os.stat(
            selected.name, dir_fd=directory_descriptor, follow_symlinks=False
        )
        if (
            _file_witness(before) != _file_witness(after)
            or _inode(after) != _inode(path_after)
            or len(content) != before.st_size
        ):
            raise fail(
                "DATABASE_WRITE_STOP_WITNESS_UNTRUSTED",
                "Database write-stop witness is not trusted.",
            )
    except DeploymentError:
        raise
    except OSError:
        raise fail(
            "DATABASE_WRITE_STOP_WITNESS_UNTRUSTED",
            "Database write-stop witness is not trusted.",
        ) from None
    finally:
        if descriptor is not None:
            os.close(descriptor)
        if directory_descriptor is not None:
            os.close(directory_descriptor)

    try:
        raw = json.loads(content, object_pairs_hook=_unique_object)
        witness = WriteStopWitness.from_dict(raw)
    except (DeploymentError, UnicodeDecodeError, ValueError, json.JSONDecodeError):
        raise fail(
            "DATABASE_WRITE_STOP_WITNESS_INVALID",
            "Database write-stop witness is invalid.",
        ) from None
    if canonical_json_bytes(witness.to_dict()) != content:
        raise fail(
            "DATABASE_WRITE_STOP_WITNESS_INVALID",
            "Database write-stop witness is invalid.",
        )
    if (
        witness.task_identity != expected_task_identity
        or witness.deployment_identity != expected_deployment_identity
        or witness.operator_state_identity != expected_operator_state_identity
        or witness.database_path_identity != expected_database_path_identity
        or witness.database_device != expected_database_stat.st_dev
        or witness.database_inode != expected_database_stat.st_ino
    ):
        raise fail(
            "DATABASE_WRITE_STOP_WITNESS_MISMATCH",
            "Database write-stop witness does not match this operation.",
        )
    return witness


def _open_database(path: Path) -> _OpenedDatabase:
    selected = _absolute_path(path, existing=True, code="DATABASE_PATH_INVALID")
    parent_stat = selected.parent.lstat()
    if (
        not stat.S_ISDIR(parent_stat.st_mode)
        or stat.S_ISLNK(parent_stat.st_mode)
        or parent_stat.st_mode & 0o022
    ):
        raise fail("DATABASE_PATH_UNSAFE", "Database path boundary is unsafe.")
    directory_flags = (
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    file_flags = (
        os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    )
    directory_descriptor: int | None = None
    descriptor: int | None = None
    try:
        directory_descriptor = os.open(selected.parent, directory_flags)
        if _inode(os.fstat(directory_descriptor)) != _inode(parent_stat):
            raise OSError
        descriptor = os.open(selected.name, file_flags, dir_fd=directory_descriptor)
        observed = os.fstat(descriptor)
        if (
            not stat.S_ISREG(observed.st_mode)
            or observed.st_nlink != 1
            or observed.st_uid != parent_stat.st_uid
            or observed.st_gid != parent_stat.st_gid
            or stat.S_IMODE(observed.st_mode) & 0o022
            or not 0 < observed.st_size <= _MAX_DATABASE_BYTES
        ):
            raise fail("DATABASE_INVALID", "Database is invalid.")
        return _OpenedDatabase(selected, directory_descriptor, descriptor, observed)
    except DeploymentError:
        if descriptor is not None:
            os.close(descriptor)
        if directory_descriptor is not None:
            os.close(directory_descriptor)
        raise
    except OSError:
        if descriptor is not None:
            os.close(descriptor)
        if directory_descriptor is not None:
            os.close(directory_descriptor)
        raise fail("DATABASE_UNAVAILABLE", "Database is unavailable.") from None


def _require_quiescent_sidecars(opened: _OpenedDatabase) -> None:
    for suffix in ("-wal", "-shm"):
        name = f"{opened.path.name}{suffix}"
        try:
            observed = os.stat(
                name,
                dir_fd=opened.directory_descriptor,
                follow_symlinks=False,
            )
        except FileNotFoundError:
            continue
        except OSError:
            raise fail(
                "DATABASE_WAL_NOT_QUIESCENT",
                "Database WAL/SHM state is not quiescent.",
            ) from None
        if not stat.S_ISREG(observed.st_mode) or observed.st_size != 0:
            raise fail(
                "DATABASE_WAL_NOT_QUIESCENT",
                "Database WAL/SHM state is not quiescent.",
            )


def _inspect_descriptor(opened: _OpenedDatabase) -> DatabaseInspection:
    before = os.fstat(opened.descriptor)
    uri_path = quote(f"/proc/self/fd/{opened.descriptor}", safe="/")
    connection = sqlite3.connect(
        f"file:{uri_path}?mode=ro&immutable=1", uri=True, timeout=1.0
    )
    try:
        connection.execute("PRAGMA query_only = ON")
        integrity = tuple(
            row[0] for row in connection.execute("PRAGMA integrity_check")
        )
        if integrity != ("ok",):
            raise fail("DATABASE_INTEGRITY_FAILED", "Database integrity failed.")
        heads = tuple(
            row[0]
            for row in connection.execute(
                "SELECT version_num FROM alembic_version ORDER BY version_num"
            )
        )
    except DeploymentError:
        raise
    except sqlite3.DatabaseError:
        raise fail("DATABASE_INVALID", "Database is invalid.") from None
    finally:
        connection.close()
    schema_heads = _schema_heads(heads, code="DATABASE_INVALID")
    sha256, size = _hash_descriptor(opened.descriptor)
    after = os.fstat(opened.descriptor)
    if _file_witness(before) != _file_witness(after):
        raise fail("DATABASE_CHANGED", "Database changed during inspection.")
    return DatabaseInspection(
        schema_heads=schema_heads,
        size_bytes=size,
        sha256=sha256,
        device=after.st_dev,
        inode=after.st_ino,
    )


def _require_unchanged(opened: _OpenedDatabase) -> None:
    after = os.fstat(opened.descriptor)
    try:
        path_after = os.stat(
            opened.path.name,
            dir_fd=opened.directory_descriptor,
            follow_symlinks=False,
        )
    except OSError:
        raise fail("DATABASE_CHANGED", "Database changed during inspection.") from None
    if _file_witness(opened.initial_stat) != _file_witness(after) or _inode(
        after
    ) != _inode(path_after):
        raise fail("DATABASE_CHANGED", "Database changed during inspection.")


def _hash_descriptor(descriptor: int) -> tuple[str, int]:
    os.lseek(descriptor, 0, os.SEEK_SET)
    digest = hashlib.sha256()
    size = 0
    while True:
        chunk = os.read(descriptor, _READ_CHUNK_BYTES)
        if not chunk:
            break
        size += len(chunk)
        if size > _MAX_DATABASE_BYTES:
            raise fail("DATABASE_INVALID", "Database is invalid.")
        digest.update(chunk)
    return digest.hexdigest(), size


def _read_bounded(descriptor: int, limit: int) -> bytes:
    chunks: list[bytes] = []
    remaining = limit + 1
    while remaining > 0:
        chunk = os.read(descriptor, min(_READ_CHUNK_BYTES, remaining))
        if not chunk:
            break
        chunks.append(chunk)
        remaining -= len(chunk)
    content = b"".join(chunks)
    if len(content) > limit:
        raise fail(
            "DATABASE_WRITE_STOP_WITNESS_INVALID",
            "Database write-stop witness is invalid.",
        )
    return content


def _schema_heads(value: object, *, code: str) -> tuple[str, ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise fail(code, "Database schema evidence is invalid.")
    heads = tuple(value)
    if (
        not heads
        or any(
            not isinstance(item, str) or _SCHEMA_HEAD.fullmatch(item) is None
            for item in heads
        )
        or tuple(sorted(set(heads))) != heads
    ):
        raise fail(code, "Database schema evidence is invalid.")
    return heads


def _absolute_path(path: Path, *, existing: bool, code: str) -> Path:
    if not isinstance(path, Path) or not path.is_absolute():
        raise fail(code, "Database path is invalid.")
    rendered = str(path)
    if len(rendered) > 4096 or any(
        character in rendered for character in ("\x00", "\n", "\r")
    ):
        raise fail(code, "Database path is invalid.")
    try:
        resolved = path.resolve(strict=existing)
    except OSError:
        raise fail(code, "Database path is invalid.") from None
    if resolved != path:
        raise fail(code, "Database path boundary is unsafe.")
    return path


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("duplicate key")
        result[key] = value
    return result


def _inode(value: os.stat_result) -> tuple[int, int]:
    return value.st_dev, value.st_ino


def _file_witness(value: os.stat_result) -> tuple[int, ...]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_uid,
        value.st_gid,
        value.st_nlink,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


__all__ = [
    "DATABASE_INSPECTION_SCHEMA",
    "DATABASE_SCHEMA_ADMISSION_SCHEMA",
    "DATABASE_WRITER_UNITS",
    "WRITE_STOP_ISSUER",
    "WRITE_STOP_WITNESS_IDENTITY_SCHEME",
    "WRITE_STOP_WITNESS_SCHEMA",
    "DatabaseInspection",
    "DatabaseSchemaAdmission",
    "SchemaAdmissionProvider",
    "StoppedWriter",
    "WriteStopWitness",
    "UnavailableSchemaAdmissionProvider",
    "backup_database",
    "database_path_identity",
    "inspect_database",
    "write_stop_witness_path",
]
