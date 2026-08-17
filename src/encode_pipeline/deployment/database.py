"""Pinned SQLite inspection and root-operator controlled offline backups.

Backups consume a root-owned write-stop witness exactly once and recovery
preserves displaced bytes in root-only immutable slots.  This module does not
run migrations or downgrades, and its cleanup is confined to fixed database
candidate coordinates rather than user workspace or artifact data.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager
import ctypes
from dataclasses import dataclass
from datetime import datetime, timezone
import errno
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import sqlite3
import stat
from typing import Any, Protocol
import uuid
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
DATABASE_BACKUP_RECEIPT_SCHEMA = "helixweave-database-backup-receipt-v1"
DATABASE_FILE_EVIDENCE_SCHEMA = "helixweave-database-file-evidence-v1"
DATABASE_RESTORE_RECEIPT_SCHEMA = "helixweave-database-restore-receipt-v1"
DATABASE_QUARANTINE_RECEIPT_SCHEMA = "helixweave-database-quarantine-receipt-v1"
DATABASE_RAW_FILE_EVIDENCE_SCHEMA = "helixweave-database-raw-file-evidence-v1"
DATABASE_INVALID_FRESH_QUARANTINE_RECEIPT_SCHEMA = (
    "helixweave-database-invalid-fresh-quarantine-receipt-v1"
)
WRITE_STOP_WITNESS_IDENTITY_SCHEME = (
    "helixweave-database-write-stop-witness-identity-v2"
)
DATABASE_PATH_IDENTITY_SCHEME = "helixweave-database-path-identity-v1"
DATABASE_SCHEMA_ADMISSION_IDENTITY_SCHEME = (
    "helixweave-database-schema-admission-identity-v1"
)
DATABASE_BACKUP_IDENTITY_SCHEME = "helixweave-database-backup-identity-v1"
DATABASE_BACKUP_RECEIPT_IDENTITY_SCHEME = (
    "helixweave-database-backup-receipt-identity-v1"
)
DATABASE_FILE_EVIDENCE_IDENTITY_SCHEME = "helixweave-database-file-evidence-identity-v1"
DATABASE_RESTORE_RECEIPT_IDENTITY_SCHEME = (
    "helixweave-database-restore-receipt-identity-v1"
)
DATABASE_QUARANTINE_RECEIPT_IDENTITY_SCHEME = (
    "helixweave-database-quarantine-receipt-identity-v1"
)
DATABASE_RAW_FILE_EVIDENCE_IDENTITY_SCHEME = (
    "helixweave-database-raw-file-evidence-identity-v1"
)
DATABASE_INVALID_FRESH_QUARANTINE_REQUEST_IDENTITY_SCHEME = (
    "helixweave-database-invalid-fresh-quarantine-request-identity-v1"
)
DATABASE_INVALID_FRESH_RAW_EVIDENCE_IDENTITY_SCHEME = (
    "helixweave-database-invalid-fresh-raw-evidence-identity-v1"
)
DATABASE_INVALID_FRESH_QUARANTINE_RECEIPT_IDENTITY_SCHEME = (
    "helixweave-database-invalid-fresh-quarantine-receipt-identity-v1"
)
DATABASE_CONTENT_IDENTITY_SCHEME = "helixweave-database-content-identity-v1"
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
_MAX_RECEIPT_BYTES = 256 * 1024
_MAX_DATABASE_BYTES = 16 * 1024**4
_READ_CHUNK_BYTES = 1024 * 1024
_RENAME_NOREPLACE = 1
_INVALID_FRESH_FILE_NAMES = {
    "database": "database.raw",
    "journal": "journal.raw",
    "shm": "shm.raw",
    "wal": "wal.raw",
}
DatabaseBackupFault = Callable[[str], None]
DatabaseRecoveryFault = Callable[[str], None]


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


def _nonnegative_integer(value: object, *, code: str) -> int:
    if (
        not isinstance(value, int)
        or isinstance(value, bool)
        or not 0 <= value <= 2**63 - 1
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

    @classmethod
    def from_dict(cls, raw: object) -> "DatabaseInspection":
        code = "DATABASE_BACKUP_RECEIPT_INVALID"
        value = _document(raw, code=code)
        if (
            set(value)
            != {
                "schema_version",
                "integrity",
                "schema_heads",
                "size_bytes",
                "sha256",
                "device",
                "inode",
                "wal_state",
            }
            or value["schema_version"] != DATABASE_INSPECTION_SCHEMA
            or value["integrity"] != "ok"
            or value["wal_state"] != "absent-or-empty"
            or not isinstance(value["sha256"], str)
            or re.fullmatch(r"[0-9a-f]{64}", value["sha256"]) is None
            or not isinstance(value["size_bytes"], int)
            or isinstance(value["size_bytes"], bool)
            or not 0 < value["size_bytes"] <= _MAX_DATABASE_BYTES
        ):
            raise fail(code, "Database backup receipt is invalid.")
        return cls(
            schema_heads=_schema_heads(value["schema_heads"], code=code),
            size_bytes=value["size_bytes"],
            sha256=value["sha256"],
            device=_positive_integer(value["device"], code=code),
            inode=_positive_integer(value["inode"], code=code),
            wal_state="absent-or-empty",
        )

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


@dataclass(frozen=True)
class DatabaseBackupReceipt:
    """Path-free identity proof for one durable offline database backup."""

    identity: str
    backup_identity: str
    task_identity: str
    deployment_identity: str
    operator_state_identity: str
    witness_identity: str
    consume_once_identity: str
    schema_admission_identity: str
    source: DatabaseInspection
    backup: DatabaseInspection

    @classmethod
    def create(
        cls,
        *,
        backup_identity: str,
        task_identity: str,
        deployment_identity: str,
        operator_state_identity: str,
        witness_identity: str,
        consume_once_identity: str,
        schema_admission_identity: str,
        source: DatabaseInspection,
        backup: DatabaseInspection,
    ) -> "DatabaseBackupReceipt":
        value: dict[str, object] = {
            "schema_version": DATABASE_BACKUP_RECEIPT_SCHEMA,
            "backup_identity": backup_identity,
            "task_identity": task_identity,
            "deployment_identity": deployment_identity,
            "operator_state_identity": operator_state_identity,
            "witness_identity": witness_identity,
            "consume_once_identity": consume_once_identity,
            "schema_admission_identity": schema_admission_identity,
            "source": source.to_dict(),
            "backup": backup.to_dict(),
        }
        value["identity"] = canonical_identity(
            value, scheme=DATABASE_BACKUP_RECEIPT_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "DatabaseBackupReceipt":
        code = "DATABASE_BACKUP_RECEIPT_INVALID"
        value = _document(raw, code=code)
        expected_keys = {
            "schema_version",
            "identity",
            "backup_identity",
            "task_identity",
            "deployment_identity",
            "operator_state_identity",
            "witness_identity",
            "consume_once_identity",
            "schema_admission_identity",
            "source",
            "backup",
        }
        if set(value) != expected_keys or value["schema_version"] != (
            DATABASE_BACKUP_RECEIPT_SCHEMA
        ):
            raise fail(code, "Database backup receipt is invalid.")
        identity = _content_identity(value["identity"], code=code)
        if identity != canonical_identity(
            without_key(value, "identity"),
            scheme=DATABASE_BACKUP_RECEIPT_IDENTITY_SCHEME,
        ):
            raise fail(code, "Database backup receipt is invalid.")
        return cls(
            identity=identity,
            backup_identity=_content_identity(value["backup_identity"], code=code),
            task_identity=_task_identity(value["task_identity"], code=code),
            deployment_identity=_content_identity(
                value["deployment_identity"], code=code
            ),
            operator_state_identity=_content_identity(
                value["operator_state_identity"], code=code
            ),
            witness_identity=_content_identity(value["witness_identity"], code=code),
            consume_once_identity=_content_identity(
                value["consume_once_identity"], code=code
            ),
            schema_admission_identity=_content_identity(
                value["schema_admission_identity"], code=code
            ),
            source=DatabaseInspection.from_dict(value["source"]),
            backup=DatabaseInspection.from_dict(value["backup"]),
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": DATABASE_BACKUP_RECEIPT_SCHEMA,
            "identity": self.identity,
            "backup_identity": self.backup_identity,
            "task_identity": self.task_identity,
            "deployment_identity": self.deployment_identity,
            "operator_state_identity": self.operator_state_identity,
            "witness_identity": self.witness_identity,
            "consume_once_identity": self.consume_once_identity,
            "schema_admission_identity": self.schema_admission_identity,
            "source": self.source.to_dict(),
            "backup": self.backup.to_dict(),
        }


@dataclass(frozen=True)
class DatabaseFileEvidence:
    """Path-free identity for trusted bytes, without claiming SQLite validity."""

    identity: str
    size_bytes: int
    sha256: str
    device: int
    inode: int

    @classmethod
    def create(
        cls,
        *,
        size_bytes: int,
        sha256: str,
        device: int,
        inode: int,
    ) -> "DatabaseFileEvidence":
        value: dict[str, object] = {
            "schema_version": DATABASE_FILE_EVIDENCE_SCHEMA,
            "size_bytes": size_bytes,
            "sha256": sha256,
            "device": device,
            "inode": inode,
        }
        value["identity"] = canonical_identity(
            value, scheme=DATABASE_FILE_EVIDENCE_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "DatabaseFileEvidence":
        code = "DATABASE_RECOVERY_RECEIPT_INVALID"
        value = _document(raw, code=code)
        if (
            set(value)
            != {
                "schema_version",
                "identity",
                "size_bytes",
                "sha256",
                "device",
                "inode",
            }
            or value["schema_version"] != DATABASE_FILE_EVIDENCE_SCHEMA
            or not isinstance(value["sha256"], str)
            or re.fullmatch(r"[0-9a-f]{64}", value["sha256"]) is None
            or not isinstance(value["size_bytes"], int)
            or isinstance(value["size_bytes"], bool)
            or not 0 < value["size_bytes"] <= _MAX_DATABASE_BYTES
        ):
            raise fail(code, "Database recovery receipt is invalid.")
        identity = _content_identity(value["identity"], code=code)
        if identity != canonical_identity(
            without_key(value, "identity"),
            scheme=DATABASE_FILE_EVIDENCE_IDENTITY_SCHEME,
        ):
            raise fail(code, "Database recovery receipt is invalid.")
        return cls(
            identity=identity,
            size_bytes=value["size_bytes"],
            sha256=value["sha256"],
            device=_positive_integer(value["device"], code=code),
            inode=_positive_integer(value["inode"], code=code),
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": DATABASE_FILE_EVIDENCE_SCHEMA,
            "identity": self.identity,
            "size_bytes": self.size_bytes,
            "sha256": self.sha256,
            "device": self.device,
            "inode": self.inode,
        }


@dataclass(frozen=True)
class DatabaseRawFileEvidence:
    """Path-free bytes and metadata for one possibly-invalid SQLite file."""

    identity: str
    size_bytes: int
    sha256: str
    device: int
    inode: int
    uid: int
    gid: int
    mode: int
    nlink: int
    mtime_ns: int
    ctime_ns: int

    @classmethod
    def create(
        cls,
        *,
        size_bytes: int,
        sha256: str,
        device: int,
        inode: int,
        uid: int,
        gid: int,
        mode: int,
        nlink: int,
        mtime_ns: int,
        ctime_ns: int,
    ) -> "DatabaseRawFileEvidence":
        value: dict[str, object] = {
            "schema_version": DATABASE_RAW_FILE_EVIDENCE_SCHEMA,
            "size_bytes": size_bytes,
            "sha256": sha256,
            "device": device,
            "inode": inode,
            "uid": uid,
            "gid": gid,
            "mode": mode,
            "nlink": nlink,
            "mtime_ns": mtime_ns,
            "ctime_ns": ctime_ns,
        }
        value["identity"] = canonical_identity(
            value, scheme=DATABASE_RAW_FILE_EVIDENCE_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "DatabaseRawFileEvidence":
        code = "DATABASE_INVALID_FRESH_QUARANTINE_RECEIPT_INVALID"
        value = _document(raw, code=code)
        if (
            set(value)
            != {
                "schema_version",
                "identity",
                "size_bytes",
                "sha256",
                "device",
                "inode",
                "uid",
                "gid",
                "mode",
                "nlink",
                "mtime_ns",
                "ctime_ns",
            }
            or value["schema_version"] != DATABASE_RAW_FILE_EVIDENCE_SCHEMA
            or not isinstance(value["sha256"], str)
            or re.fullmatch(r"[0-9a-f]{64}", value["sha256"]) is None
            or not isinstance(value["size_bytes"], int)
            or isinstance(value["size_bytes"], bool)
            or not 0 <= value["size_bytes"] <= _MAX_DATABASE_BYTES
            or not isinstance(value["mode"], int)
            or isinstance(value["mode"], bool)
            or value["mode"] not in {0o400, 0o660}
            or not isinstance(value["nlink"], int)
            or isinstance(value["nlink"], bool)
            or value["nlink"] != 1
        ):
            raise fail(code, "Invalid fresh database quarantine receipt is invalid.")
        identity = _content_identity(value["identity"], code=code)
        if identity != canonical_identity(
            without_key(value, "identity"),
            scheme=DATABASE_RAW_FILE_EVIDENCE_IDENTITY_SCHEME,
        ):
            raise fail(code, "Invalid fresh database quarantine receipt is invalid.")
        return cls(
            identity=identity,
            size_bytes=value["size_bytes"],
            sha256=value["sha256"],
            device=_positive_integer(value["device"], code=code),
            inode=_positive_integer(value["inode"], code=code),
            uid=_nonnegative_integer(value["uid"], code=code),
            gid=_nonnegative_integer(value["gid"], code=code),
            mode=value["mode"],
            nlink=value["nlink"],
            mtime_ns=_nonnegative_integer(value["mtime_ns"], code=code),
            ctime_ns=_nonnegative_integer(value["ctime_ns"], code=code),
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": DATABASE_RAW_FILE_EVIDENCE_SCHEMA,
            "identity": self.identity,
            "size_bytes": self.size_bytes,
            "sha256": self.sha256,
            "device": self.device,
            "inode": self.inode,
            "uid": self.uid,
            "gid": self.gid,
            "mode": self.mode,
            "nlink": self.nlink,
            "mtime_ns": self.mtime_ns,
            "ctime_ns": self.ctime_ns,
        }


@dataclass(frozen=True)
class InvalidFreshDatabaseFileEvidence:
    """One fixed candidate role and its immutable preserved copy."""

    role: str
    source: DatabaseRawFileEvidence
    evidence: DatabaseRawFileEvidence

    @classmethod
    def from_dict(cls, raw: object) -> "InvalidFreshDatabaseFileEvidence":
        code = "DATABASE_INVALID_FRESH_QUARANTINE_RECEIPT_INVALID"
        value = _document(raw, code=code)
        if (
            set(value) != {"role", "source", "evidence"}
            or not isinstance(value["role"], str)
            or value["role"] not in _INVALID_FRESH_FILE_NAMES
        ):
            raise fail(code, "Invalid fresh database quarantine receipt is invalid.")
        try:
            source = DatabaseRawFileEvidence.from_dict(value["source"])
            evidence = DatabaseRawFileEvidence.from_dict(value["evidence"])
        except DeploymentError:
            raise fail(
                code, "Invalid fresh database quarantine receipt is invalid."
            ) from None
        if (
            source.mode != 0o660
            or evidence.mode != 0o400
            or evidence.uid != _ROOT_UID
            or evidence.gid != _ROOT_GID
            or source.sha256 != evidence.sha256
            or source.size_bytes != evidence.size_bytes
        ):
            raise fail(code, "Invalid fresh database quarantine receipt is invalid.")
        return cls(role=value["role"], source=source, evidence=evidence)

    def to_dict(self) -> dict[str, object]:
        return {
            "role": self.role,
            "source": self.source.to_dict(),
            "evidence": self.evidence.to_dict(),
        }


@dataclass(frozen=True)
class InvalidFreshDatabaseQuarantineReceipt:
    """Root-owned proof that raw failed-candidate bytes were retained."""

    identity: str
    request_identity: str
    task_identity: str
    deployment_identity: str
    prior_state_identity: str
    raw_evidence_identity: str
    files: tuple[InvalidFreshDatabaseFileEvidence, ...]

    @classmethod
    def create(
        cls,
        *,
        request_identity: str,
        task_identity: str,
        deployment_identity: str,
        prior_state_identity: str,
        files: Sequence[InvalidFreshDatabaseFileEvidence],
    ) -> "InvalidFreshDatabaseQuarantineReceipt":
        ordered = tuple(sorted(files, key=lambda item: item.role))
        raw_evidence_identity = canonical_identity(
            {
                "files": [
                    {"role": item.role, "source": item.source.to_dict()}
                    for item in ordered
                ]
            },
            scheme=DATABASE_INVALID_FRESH_RAW_EVIDENCE_IDENTITY_SCHEME,
        )
        value: dict[str, object] = {
            "schema_version": DATABASE_INVALID_FRESH_QUARANTINE_RECEIPT_SCHEMA,
            "status": "raw-evidence-preserved",
            "request_identity": request_identity,
            "task_identity": task_identity,
            "deployment_identity": deployment_identity,
            "prior_state_identity": prior_state_identity,
            "raw_evidence_identity": raw_evidence_identity,
            "files": [item.to_dict() for item in ordered],
        }
        value["identity"] = canonical_identity(
            value,
            scheme=DATABASE_INVALID_FRESH_QUARANTINE_RECEIPT_IDENTITY_SCHEME,
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "InvalidFreshDatabaseQuarantineReceipt":
        code = "DATABASE_INVALID_FRESH_QUARANTINE_RECEIPT_INVALID"
        value = _document(raw, code=code)
        if (
            set(value)
            != {
                "schema_version",
                "identity",
                "status",
                "request_identity",
                "task_identity",
                "deployment_identity",
                "prior_state_identity",
                "raw_evidence_identity",
                "files",
            }
            or value["schema_version"]
            != DATABASE_INVALID_FRESH_QUARANTINE_RECEIPT_SCHEMA
            or value["status"] != "raw-evidence-preserved"
            or not isinstance(value["files"], list)
        ):
            raise fail(code, "Invalid fresh database quarantine receipt is invalid.")
        try:
            files = tuple(
                InvalidFreshDatabaseFileEvidence.from_dict(item)
                for item in value["files"]
            )
        except DeploymentError:
            raise fail(
                code, "Invalid fresh database quarantine receipt is invalid."
            ) from None
        roles = tuple(item.role for item in files)
        if not files or "database" not in roles or roles != tuple(sorted(set(roles))):
            raise fail(code, "Invalid fresh database quarantine receipt is invalid.")
        identity = _content_identity(value["identity"], code=code)
        request_identity = _content_identity(value["request_identity"], code=code)
        raw_evidence_identity = _content_identity(
            value["raw_evidence_identity"], code=code
        )
        if raw_evidence_identity != canonical_identity(
            {
                "files": [
                    {"role": item.role, "source": item.source.to_dict()}
                    for item in files
                ]
            },
            scheme=DATABASE_INVALID_FRESH_RAW_EVIDENCE_IDENTITY_SCHEME,
        ) or identity != canonical_identity(
            without_key(value, "identity"),
            scheme=DATABASE_INVALID_FRESH_QUARANTINE_RECEIPT_IDENTITY_SCHEME,
        ):
            raise fail(code, "Invalid fresh database quarantine receipt is invalid.")
        return cls(
            identity=identity,
            request_identity=request_identity,
            task_identity=_task_identity(value["task_identity"], code=code),
            deployment_identity=_content_identity(
                value["deployment_identity"], code=code
            ),
            prior_state_identity=_content_identity(
                value["prior_state_identity"], code=code
            ),
            raw_evidence_identity=raw_evidence_identity,
            files=files,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": DATABASE_INVALID_FRESH_QUARANTINE_RECEIPT_SCHEMA,
            "identity": self.identity,
            "status": "raw-evidence-preserved",
            "request_identity": self.request_identity,
            "task_identity": self.task_identity,
            "deployment_identity": self.deployment_identity,
            "prior_state_identity": self.prior_state_identity,
            "raw_evidence_identity": self.raw_evidence_identity,
            "files": [item.to_dict() for item in self.files],
        }


@dataclass(frozen=True)
class DatabaseRestoreReceipt:
    """Immutable proof that failed live bytes were preserved before restore."""

    identity: str
    backup_identity: str
    backup_receipt_identity: str
    task_identity: str
    deployment_identity: str
    prior_state_identity: str
    source_identity: str
    schema_heads: tuple[str, ...]
    failed_current: DatabaseFileEvidence
    evidence: DatabaseFileEvidence
    target: DatabaseInspection

    @classmethod
    def create(
        cls,
        *,
        backup_identity: str,
        backup_receipt_identity: str,
        task_identity: str,
        deployment_identity: str,
        prior_state_identity: str,
        source_identity: str,
        schema_heads: Sequence[str],
        failed_current: DatabaseFileEvidence,
        evidence: DatabaseFileEvidence,
        target: DatabaseInspection,
    ) -> "DatabaseRestoreReceipt":
        value: dict[str, object] = {
            "schema_version": DATABASE_RESTORE_RECEIPT_SCHEMA,
            "status": "recovery-evidence-preserved",
            "backup_identity": backup_identity,
            "backup_receipt_identity": backup_receipt_identity,
            "task_identity": task_identity,
            "deployment_identity": deployment_identity,
            "prior_state_identity": prior_state_identity,
            "source_identity": source_identity,
            "schema_heads": list(schema_heads),
            "failed_current": failed_current.to_dict(),
            "evidence": evidence.to_dict(),
            "target": target.to_dict(),
        }
        value["identity"] = canonical_identity(
            value, scheme=DATABASE_RESTORE_RECEIPT_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "DatabaseRestoreReceipt":
        code = "DATABASE_RESTORE_RECEIPT_INVALID"
        value = _document(raw, code=code)
        if (
            set(value)
            != {
                "schema_version",
                "identity",
                "status",
                "backup_identity",
                "backup_receipt_identity",
                "task_identity",
                "deployment_identity",
                "prior_state_identity",
                "source_identity",
                "schema_heads",
                "failed_current",
                "evidence",
                "target",
            }
            or value["schema_version"] != DATABASE_RESTORE_RECEIPT_SCHEMA
            or value["status"] != "recovery-evidence-preserved"
        ):
            raise fail(code, "Database restore receipt is invalid.")
        identity = _content_identity(value["identity"], code=code)
        if identity != canonical_identity(
            without_key(value, "identity"),
            scheme=DATABASE_RESTORE_RECEIPT_IDENTITY_SCHEME,
        ):
            raise fail(code, "Database restore receipt is invalid.")
        try:
            failed_current = DatabaseFileEvidence.from_dict(value["failed_current"])
            evidence = DatabaseFileEvidence.from_dict(value["evidence"])
            target = DatabaseInspection.from_dict(value["target"])
        except DeploymentError:
            raise fail(code, "Database restore receipt is invalid.") from None
        return cls(
            identity=identity,
            backup_identity=_content_identity(value["backup_identity"], code=code),
            backup_receipt_identity=_content_identity(
                value["backup_receipt_identity"], code=code
            ),
            task_identity=_task_identity(value["task_identity"], code=code),
            deployment_identity=_content_identity(
                value["deployment_identity"], code=code
            ),
            prior_state_identity=_content_identity(
                value["prior_state_identity"], code=code
            ),
            source_identity=_content_identity(value["source_identity"], code=code),
            schema_heads=_schema_heads(value["schema_heads"], code=code),
            failed_current=failed_current,
            evidence=evidence,
            target=target,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": DATABASE_RESTORE_RECEIPT_SCHEMA,
            "identity": self.identity,
            "status": "recovery-evidence-preserved",
            "backup_identity": self.backup_identity,
            "backup_receipt_identity": self.backup_receipt_identity,
            "task_identity": self.task_identity,
            "deployment_identity": self.deployment_identity,
            "prior_state_identity": self.prior_state_identity,
            "source_identity": self.source_identity,
            "schema_heads": list(self.schema_heads),
            "failed_current": self.failed_current.to_dict(),
            "evidence": self.evidence.to_dict(),
            "target": self.target.to_dict(),
        }


@dataclass(frozen=True)
class DatabaseQuarantineReceipt:
    """Immutable proof for a failed fresh candidate removed from live storage."""

    identity: str
    task_identity: str
    deployment_identity: str
    prior_state_identity: str
    candidate_identity: str
    schema_heads: tuple[str, ...]
    source_coordinate: str
    source: DatabaseFileEvidence
    evidence: DatabaseFileEvidence

    @classmethod
    def create(
        cls,
        *,
        task_identity: str,
        deployment_identity: str,
        prior_state_identity: str,
        candidate_identity: str,
        schema_heads: Sequence[str],
        source_coordinate: str,
        source: DatabaseFileEvidence,
        evidence: DatabaseFileEvidence,
    ) -> "DatabaseQuarantineReceipt":
        value: dict[str, object] = {
            "schema_version": DATABASE_QUARANTINE_RECEIPT_SCHEMA,
            "status": "quarantined",
            "task_identity": task_identity,
            "deployment_identity": deployment_identity,
            "prior_state_identity": prior_state_identity,
            "candidate_identity": candidate_identity,
            "schema_heads": list(schema_heads),
            "source_coordinate": source_coordinate,
            "source": source.to_dict(),
            "evidence": evidence.to_dict(),
        }
        value["identity"] = canonical_identity(
            value, scheme=DATABASE_QUARANTINE_RECEIPT_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "DatabaseQuarantineReceipt":
        code = "DATABASE_QUARANTINE_RECEIPT_INVALID"
        value = _document(raw, code=code)
        if (
            set(value)
            != {
                "schema_version",
                "identity",
                "status",
                "task_identity",
                "deployment_identity",
                "prior_state_identity",
                "candidate_identity",
                "schema_heads",
                "source_coordinate",
                "source",
                "evidence",
            }
            or value["schema_version"] != DATABASE_QUARANTINE_RECEIPT_SCHEMA
            or value["status"] != "quarantined"
            or value["source_coordinate"] not in {"candidate", "canonical"}
        ):
            raise fail(code, "Database quarantine receipt is invalid.")
        identity = _content_identity(value["identity"], code=code)
        if identity != canonical_identity(
            without_key(value, "identity"),
            scheme=DATABASE_QUARANTINE_RECEIPT_IDENTITY_SCHEME,
        ):
            raise fail(code, "Database quarantine receipt is invalid.")
        try:
            source = DatabaseFileEvidence.from_dict(value["source"])
            evidence = DatabaseFileEvidence.from_dict(value["evidence"])
        except DeploymentError:
            raise fail(code, "Database quarantine receipt is invalid.") from None
        return cls(
            identity=identity,
            task_identity=_task_identity(value["task_identity"], code=code),
            deployment_identity=_content_identity(
                value["deployment_identity"], code=code
            ),
            prior_state_identity=_content_identity(
                value["prior_state_identity"], code=code
            ),
            candidate_identity=_content_identity(
                value["candidate_identity"], code=code
            ),
            schema_heads=_schema_heads(value["schema_heads"], code=code),
            source_coordinate=value["source_coordinate"],
            source=source,
            evidence=evidence,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": DATABASE_QUARANTINE_RECEIPT_SCHEMA,
            "identity": self.identity,
            "status": "quarantined",
            "task_identity": self.task_identity,
            "deployment_identity": self.deployment_identity,
            "prior_state_identity": self.prior_state_identity,
            "candidate_identity": self.candidate_identity,
            "schema_heads": list(self.schema_heads),
            "source_coordinate": self.source_coordinate,
            "source": self.source.to_dict(),
            "evidence": self.evidence.to_dict(),
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
    shared_group_write: bool

    def close(self) -> None:
        os.close(self.descriptor)
        os.close(self.directory_descriptor)


@dataclass
class _OpenedRawDatabaseFile:
    role: str
    name: str
    descriptor: int
    initial_stat: os.stat_result

    def close(self) -> None:
        os.close(self.descriptor)


def database_path_identity(path: Path) -> str:
    selected = _absolute_path(path, existing=True, code="DATABASE_PATH_INVALID")
    return canonical_identity(
        {"absolute_path": str(selected)}, scheme=DATABASE_PATH_IDENTITY_SCHEME
    )


def write_stop_witness_path(layout: DeploymentLayout, task_identity: str) -> Path:
    task = _task_identity(task_identity, code="DATABASE_BACKUP_REQUEST_INVALID")
    return layout.run_root / "operator" / "database-backups" / task / "write-stop.json"


def consumed_write_stop_witness_path(
    layout: DeploymentLayout, task_identity: str
) -> Path:
    task = _task_identity(task_identity, code="DATABASE_BACKUP_REQUEST_INVALID")
    return (
        layout.run_root
        / "operator"
        / "database-backups"
        / task
        / "write-stop.consumed.json"
    )


def create_write_stop_witness(
    layout: DeploymentLayout,
    *,
    task_identity: str,
    deployment_identity: str,
    operator_state_identity: str,
    writers: Sequence[StoppedWriter],
    expected_database_uid: int,
    expected_database_gid: int,
) -> WriteStopWitness:
    """Persist one root-owned, consume-once witness after writers stop."""
    if not isinstance(layout, DeploymentLayout):
        raise fail(
            "DATABASE_WRITE_STOP_WITNESS_INVALID",
            "Database write-stop witness is invalid.",
        )
    task = _task_identity(task_identity, code="DATABASE_WRITE_STOP_WITNESS_INVALID")
    deployment = _content_identity(
        deployment_identity, code="DATABASE_WRITE_STOP_WITNESS_INVALID"
    )
    operator_state = _content_identity(
        operator_state_identity, code="DATABASE_WRITE_STOP_WITNESS_INVALID"
    )
    parsed_writers = tuple(
        StoppedWriter.from_dict(writer.to_dict()) for writer in writers
    )
    if tuple(writer.unit for writer in parsed_writers) != DATABASE_WRITER_UNITS:
        raise fail(
            "DATABASE_WRITE_STOP_WITNESS_INVALID",
            "Database write-stop witness is invalid.",
        )
    inspection = inspect_database(
        layout.database,
        expected_owner_uid=expected_database_uid,
        expected_owner_gid=expected_database_gid,
    )
    consume_once = canonical_identity(
        {"task_identity": task, "nonce": uuid.uuid4().hex},
        scheme="helixweave-database-consume-once-identity-v1",
    )
    value: dict[str, object] = {
        "schema_version": WRITE_STOP_WITNESS_SCHEMA,
        "issuer": WRITE_STOP_ISSUER,
        "consumption_state": "ready",
        "task_identity": task,
        "deployment_identity": deployment,
        "operator_state_identity": operator_state,
        "consume_once_identity": consume_once,
        "database_path_identity": database_path_identity(layout.database),
        "database_device": inspection.device,
        "database_inode": inspection.inode,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "writers": [writer.to_dict() for writer in parsed_writers],
    }
    value["identity"] = canonical_identity(
        value, scheme=WRITE_STOP_WITNESS_IDENTITY_SCHEME
    )
    witness = WriteStopWitness.from_dict(value)
    path = write_stop_witness_path(layout, task)
    root = layout.run_root / "operator"
    boundary = root / "database-backups"
    _require_or_create_root_directory(root, create=True, mode=0o711)
    for directory in (boundary, path.parent):
        _require_or_create_root_directory(directory, create=True, mode=0o700)
    _write_root_file(path, canonical_json_bytes(witness.to_dict()), mode=0o600)
    _fsync_directory_descriptor(path.parent)
    return witness


def inspect_database(
    path: Path,
    *,
    expected_owner_uid: int,
    expected_owner_gid: int,
) -> DatabaseInspection:
    """Inspect one pinned database FD and reject every SQLite sidecar."""
    opened = _open_database(
        path,
        expected_owner_uid=expected_owner_uid,
        expected_owner_gid=expected_owner_gid,
    )
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


def database_content_identity(inspection: DatabaseInspection) -> str:
    """Return the native DB identity used by schema observation and actions."""
    if not isinstance(inspection, DatabaseInspection):
        raise fail("DATABASE_INVALID", "Database is invalid.")
    return canonical_identity(
        {
            "device": inspection.device,
            "inode": inspection.inode,
            "heads": list(inspection.schema_heads),
        },
        scheme=DATABASE_CONTENT_IDENTITY_SCHEME,
    )


def fresh_database_candidate_path(layout: DeploymentLayout, task_identity: str) -> Path:
    """Return the only supported fresh-database candidate coordinate."""

    if not isinstance(layout, DeploymentLayout):
        raise fail(
            "DATABASE_FRESH_PUBLISH_REQUEST_INVALID",
            "Fresh database publication is invalid.",
        )
    task = _task_identity(task_identity, code="DATABASE_FRESH_PUBLISH_REQUEST_INVALID")
    return layout.database.parent / f".fresh-{task}.db"


def publish_fresh_database(
    layout: DeploymentLayout,
    *,
    task_identity: str,
    expected_candidate_identity: str,
    target_schema_heads: Sequence[str],
    expected_database_uid: int,
    expected_database_gid: int,
    fault: DatabaseRecoveryFault | None = None,
) -> DatabaseInspection:
    """Atomically publish the fixed fresh candidate, or verify an exact retry."""

    if not isinstance(layout, DeploymentLayout):
        raise fail(
            "DATABASE_FRESH_PUBLISH_REQUEST_INVALID",
            "Fresh database publication is invalid.",
        )
    task = _task_identity(task_identity, code="DATABASE_FRESH_PUBLISH_REQUEST_INVALID")
    expected_identity = _content_identity(
        expected_candidate_identity, code="DATABASE_FRESH_PUBLISH_REQUEST_INVALID"
    )
    heads = _schema_heads(
        target_schema_heads, code="DATABASE_FRESH_PUBLISH_REQUEST_INVALID"
    )
    _validate_owner_ids(
        expected_database_uid,
        expected_database_gid,
        code="DATABASE_FRESH_PUBLISH_REQUEST_INVALID",
    )
    candidate = fresh_database_candidate_path(layout, task)
    canonical = layout.database
    try:
        with _database_mutation_lock(layout):
            directory_descriptor = _open_live_database_directory(
                layout,
                expected_database_uid=expected_database_uid,
                expected_database_gid=expected_database_gid,
            )
            try:
                candidate_present = _entry_present_at(
                    directory_descriptor, candidate.name
                )
                canonical_present = _entry_present_at(
                    directory_descriptor, canonical.name
                )
                if canonical_present:
                    if candidate_present:
                        raise fail(
                            "DATABASE_FRESH_PUBLISH_INVALID",
                            "Fresh database publication is invalid.",
                        )
                    opened = _open_database(
                        canonical,
                        expected_owner_uid=expected_database_uid,
                        expected_owner_gid=expected_database_gid,
                    )
                    try:
                        _require_quiescent_sidecars(opened)
                        inspection = _inspect_descriptor(opened)
                        _verify_database_target(
                            inspection,
                            expected_identity=expected_identity,
                            expected_heads=heads,
                            code="DATABASE_FRESH_PUBLISH_INVALID",
                        )
                        _require_unchanged(opened)
                        os.fsync(opened.descriptor)
                        os.fsync(opened.directory_descriptor)
                        return inspection
                    finally:
                        opened.close()
                if not candidate_present:
                    raise fail(
                        "DATABASE_FRESH_PUBLISH_INVALID",
                        "Fresh database publication is invalid.",
                    )
                _require_quiescent_sidecar_names(
                    directory_descriptor,
                    canonical.name,
                    expected_owner_uid=expected_database_uid,
                    expected_owner_gid=expected_database_gid,
                )
                opened = _open_database(
                    candidate,
                    expected_owner_uid=expected_database_uid,
                    expected_owner_gid=expected_database_gid,
                )
                try:
                    _require_quiescent_sidecars(opened)
                    inspection = _inspect_descriptor(opened)
                    _verify_database_target(
                        inspection,
                        expected_identity=expected_identity,
                        expected_heads=heads,
                        code="DATABASE_FRESH_PUBLISH_INVALID",
                    )
                    _require_quiescent_sidecars(opened)
                    _require_unchanged(opened)
                    os.fsync(opened.descriptor)
                    if fault is not None:
                        fault("fresh-candidate-synced")
                    if _entry_present_at(directory_descriptor, canonical.name):
                        raise fail(
                            "DATABASE_FRESH_PUBLISH_INVALID",
                            "Fresh database publication is invalid.",
                        )
                    _rename_noreplace_at(
                        directory_descriptor,
                        candidate.name,
                        canonical.name,
                        conflict_code="DATABASE_FRESH_PUBLISH_INVALID",
                    )
                    os.fsync(directory_descriptor)
                    if fault is not None:
                        fault("fresh-database-published")
                finally:
                    opened.close()
            finally:
                os.close(directory_descriptor)
        published = inspect_database(
            canonical,
            expected_owner_uid=expected_database_uid,
            expected_owner_gid=expected_database_gid,
        )
        _verify_database_target(
            published,
            expected_identity=expected_identity,
            expected_heads=heads,
            code="DATABASE_FRESH_PUBLISH_INVALID",
        )
        return published
    except DeploymentError:
        raise
    except Exception:
        raise fail(
            "DATABASE_FRESH_PUBLISH_FAILED",
            "Fresh database could not be published.",
            recoverable=True,
        ) from None


def restore_database_backup(
    layout: DeploymentLayout,
    *,
    backup_identity: str,
    expected_task_identity: str,
    expected_deployment_identity: str,
    expected_prior_state_identity: str,
    expected_source_identity: str,
    expected_schema_heads: Sequence[str],
    expected_database_uid: int,
    expected_database_gid: int,
    fault: DatabaseRecoveryFault | None = None,
) -> DatabaseRestoreReceipt:
    """Preserve failed live bytes, then atomically restore one verified backup."""

    if not isinstance(layout, DeploymentLayout):
        raise fail("DATABASE_RESTORE_REQUEST_INVALID", "Database restore is invalid.")
    backup = _content_identity(backup_identity, code="DATABASE_RESTORE_REQUEST_INVALID")
    task = _task_identity(
        expected_task_identity, code="DATABASE_RESTORE_REQUEST_INVALID"
    )
    deployment = _content_identity(
        expected_deployment_identity, code="DATABASE_RESTORE_REQUEST_INVALID"
    )
    prior_state = _content_identity(
        expected_prior_state_identity, code="DATABASE_RESTORE_REQUEST_INVALID"
    )
    source_identity = _content_identity(
        expected_source_identity, code="DATABASE_RESTORE_REQUEST_INVALID"
    )
    heads = _schema_heads(
        expected_schema_heads, code="DATABASE_RESTORE_REQUEST_INVALID"
    )
    _validate_owner_ids(
        expected_database_uid,
        expected_database_gid,
        code="DATABASE_RESTORE_REQUEST_INVALID",
    )
    try:
        with _database_mutation_lock(layout):
            backup_opened, backup_receipt, target = _open_restore_backup(
                layout,
                backup_identity=backup,
                task=task,
                deployment=deployment,
                prior_state=prior_state,
                source_identity=source_identity,
                schema_heads=heads,
            )
            try:
                current = _open_database(
                    layout.database,
                    expected_owner_uid=expected_database_uid,
                    expected_owner_gid=expected_database_gid,
                )
                try:
                    _require_quiescent_sidecars(current)
                    current_evidence = _file_evidence(current)
                    evidence_receipt = _load_restore_evidence(
                        layout,
                        backup_identity=backup,
                        task=task,
                        deployment=deployment,
                        prior_state=prior_state,
                        source_identity=source_identity,
                        schema_heads=heads,
                        backup_receipt=backup_receipt,
                        target=target,
                    )
                    if _opened_database_matches_target(current, target, heads):
                        if evidence_receipt is None:
                            raise fail(
                                "DATABASE_RESTORE_EVIDENCE_MISSING",
                                "Database recovery evidence is missing.",
                            )
                        _require_quiescent_sidecars(current)
                        _require_unchanged(current)
                        os.fsync(current.descriptor)
                        os.fsync(current.directory_descriptor)
                        return evidence_receipt
                    if evidence_receipt is None:
                        evidence_receipt = _preserve_restore_evidence(
                            layout,
                            opened=current,
                            backup_identity=backup,
                            backup_receipt=backup_receipt,
                            task=task,
                            deployment=deployment,
                            prior_state=prior_state,
                            source_identity=source_identity,
                            schema_heads=heads,
                            target=target,
                            fault=fault,
                        )
                    elif current_evidence != evidence_receipt.failed_current:
                        raise fail(
                            "DATABASE_RESTORE_LIVE_MISMATCH",
                            "Live database does not match recovery evidence.",
                        )
                    temporary = layout.database.parent / f".restore-{task}.db"
                    if _entry_present_at(current.directory_descriptor, temporary.name):
                        temporary_inspection = inspect_database(
                            temporary,
                            expected_owner_uid=expected_database_uid,
                            expected_owner_gid=expected_database_gid,
                        )
                        _verify_backup_copy_target(
                            temporary_inspection,
                            target,
                            heads,
                            code="DATABASE_RESTORE_TEMP_INVALID",
                        )
                    else:
                        _copy_database_descriptor_at(
                            backup_opened.descriptor,
                            current.directory_descriptor,
                            temporary.name,
                            owner_uid=expected_database_uid,
                            owner_gid=expected_database_gid,
                            mode=0o660,
                            code="DATABASE_RESTORE_FAILED",
                        )
                        if fault is not None:
                            fault("restore-copy-synced")
                        temporary_inspection = inspect_database(
                            temporary,
                            expected_owner_uid=expected_database_uid,
                            expected_owner_gid=expected_database_gid,
                        )
                        _verify_backup_copy_target(
                            temporary_inspection,
                            target,
                            heads,
                            code="DATABASE_RESTORE_TEMP_INVALID",
                        )
                    _require_quiescent_sidecars(current)
                    _require_unchanged(current)
                    if fault is not None:
                        fault("restore-ready-to-replace")
                    os.replace(
                        temporary.name,
                        layout.database.name,
                        src_dir_fd=current.directory_descriptor,
                        dst_dir_fd=current.directory_descriptor,
                    )
                    os.fsync(current.directory_descriptor)
                    if fault is not None:
                        fault("restore-live-replaced")
                finally:
                    current.close()
            finally:
                backup_opened.close()
        restored = inspect_database(
            layout.database,
            expected_owner_uid=expected_database_uid,
            expected_owner_gid=expected_database_gid,
        )
        _verify_backup_copy_target(
            restored, target, heads, code="DATABASE_RESTORE_INVALID"
        )
        return evidence_receipt
    except DeploymentError:
        raise
    except Exception:
        raise fail(
            "DATABASE_RESTORE_FAILED",
            "Database restore could not be completed.",
            recoverable=True,
        ) from None


def quarantine_fresh_database(
    layout: DeploymentLayout,
    *,
    task_identity: str,
    expected_deployment_identity: str,
    expected_prior_state_identity: str,
    expected_candidate_identity: str,
    expected_schema_heads: Sequence[str],
    expected_database_uid: int,
    expected_database_gid: int,
    fault: DatabaseRecoveryFault | None = None,
) -> DatabaseQuarantineReceipt:
    """Preserve then remove only the fixed failed fresh candidate coordinate."""

    if not isinstance(layout, DeploymentLayout):
        raise fail(
            "DATABASE_QUARANTINE_REQUEST_INVALID",
            "Fresh database quarantine is invalid.",
        )
    task = _task_identity(task_identity, code="DATABASE_QUARANTINE_REQUEST_INVALID")
    deployment = _content_identity(
        expected_deployment_identity, code="DATABASE_QUARANTINE_REQUEST_INVALID"
    )
    prior_state = _content_identity(
        expected_prior_state_identity, code="DATABASE_QUARANTINE_REQUEST_INVALID"
    )
    candidate_identity = _content_identity(
        expected_candidate_identity, code="DATABASE_QUARANTINE_REQUEST_INVALID"
    )
    heads = _schema_heads(
        expected_schema_heads, code="DATABASE_QUARANTINE_REQUEST_INVALID"
    )
    _validate_owner_ids(
        expected_database_uid,
        expected_database_gid,
        code="DATABASE_QUARANTINE_REQUEST_INVALID",
    )
    candidate = fresh_database_candidate_path(layout, task)
    try:
        with _database_mutation_lock(layout):
            directory_descriptor = _open_live_database_directory(
                layout,
                expected_database_uid=expected_database_uid,
                expected_database_gid=expected_database_gid,
            )
            try:
                candidates = tuple(
                    (coordinate, path)
                    for coordinate, path in (
                        ("candidate", candidate),
                        ("canonical", layout.database),
                    )
                    if _entry_present_at(directory_descriptor, path.name)
                )
                receipt = _load_quarantine_evidence(
                    layout,
                    task=task,
                    deployment=deployment,
                    prior_state=prior_state,
                    candidate_identity=candidate_identity,
                    schema_heads=heads,
                )
                if not candidates:
                    if receipt is None:
                        raise fail(
                            "DATABASE_QUARANTINE_SOURCE_MISSING",
                            "Fresh database quarantine source is missing.",
                        )
                    os.fsync(directory_descriptor)
                    return receipt
                if len(candidates) != 1:
                    raise fail(
                        "DATABASE_QUARANTINE_SOURCE_INVALID",
                        "Fresh database quarantine source is invalid.",
                    )
                coordinate, source_path = candidates[0]
                opened = _open_database(
                    source_path,
                    expected_owner_uid=expected_database_uid,
                    expected_owner_gid=expected_database_gid,
                )
                try:
                    _require_quiescent_sidecars(opened)
                    inspection = _inspect_descriptor(opened)
                    _verify_database_target(
                        inspection,
                        expected_identity=candidate_identity,
                        expected_heads=heads,
                        code="DATABASE_QUARANTINE_SOURCE_INVALID",
                    )
                    source_evidence = _file_evidence(opened)
                    if receipt is None:
                        receipt = _preserve_quarantine_evidence(
                            layout,
                            opened=opened,
                            task=task,
                            deployment=deployment,
                            prior_state=prior_state,
                            candidate_identity=candidate_identity,
                            schema_heads=heads,
                            source_coordinate=coordinate,
                            fault=fault,
                        )
                    elif (
                        receipt.source_coordinate != coordinate
                        or receipt.source != source_evidence
                    ):
                        raise fail(
                            "DATABASE_QUARANTINE_SOURCE_MISMATCH",
                            "Fresh database does not match quarantine evidence.",
                        )
                    _require_quiescent_sidecars(opened)
                    _require_unchanged(opened)
                    os.unlink(source_path.name, dir_fd=opened.directory_descriptor)
                    os.fsync(opened.directory_descriptor)
                    if fault is not None:
                        fault("fresh-database-quarantined")
                finally:
                    opened.close()
            finally:
                os.close(directory_descriptor)
        return receipt
    except DeploymentError:
        raise
    except Exception:
        raise fail(
            "DATABASE_QUARANTINE_FAILED",
            "Fresh database could not be quarantined.",
            recoverable=True,
        ) from None


def quarantine_invalid_fresh_database(
    layout: DeploymentLayout,
    *,
    task_identity: str,
    expected_deployment_identity: str,
    expected_prior_state_identity: str,
    expected_database_uid: int,
    expected_database_gid: int,
    fault: DatabaseRecoveryFault | None = None,
) -> InvalidFreshDatabaseQuarantineReceipt:
    """Preserve raw failed-candidate bytes without asserting SQLite validity."""

    code = "DATABASE_INVALID_FRESH_QUARANTINE_REQUEST_INVALID"
    if not isinstance(layout, DeploymentLayout):
        raise fail(code, "Invalid fresh database quarantine is invalid.")
    task = _task_identity(task_identity, code=code)
    deployment = _content_identity(expected_deployment_identity, code=code)
    prior_state = _content_identity(expected_prior_state_identity, code=code)
    _validate_owner_ids(expected_database_uid, expected_database_gid, code=code)
    request_identity = canonical_identity(
        {
            "task_identity": task,
            "deployment_identity": deployment,
            "prior_state_identity": prior_state,
        },
        scheme=DATABASE_INVALID_FRESH_QUARANTINE_REQUEST_IDENTITY_SCHEME,
    )
    candidate = fresh_database_candidate_path(layout, task)
    try:
        with _database_mutation_lock(layout):
            directory_descriptor = _open_live_database_directory(
                layout,
                expected_database_uid=expected_database_uid,
                expected_database_gid=expected_database_gid,
            )
            try:
                receipt = _load_invalid_fresh_quarantine_evidence(
                    layout,
                    request_identity=request_identity,
                    task=task,
                    deployment=deployment,
                    prior_state=prior_state,
                    expected_database_uid=expected_database_uid,
                    expected_database_gid=expected_database_gid,
                )
                present = _invalid_fresh_roles_present(
                    directory_descriptor, candidate.name
                )
                if receipt is None:
                    if "database" not in present:
                        raise fail(
                            "DATABASE_INVALID_FRESH_QUARANTINE_SOURCE_MISSING",
                            "Invalid fresh database quarantine source is missing.",
                        )
                    opened = _open_invalid_fresh_sources(
                        directory_descriptor,
                        candidate.name,
                        expected_database_uid=expected_database_uid,
                        expected_database_gid=expected_database_gid,
                    )
                    try:
                        receipt = _preserve_invalid_fresh_quarantine_evidence(
                            layout,
                            opened=opened,
                            source_directory_descriptor=directory_descriptor,
                            request_identity=request_identity,
                            task=task,
                            deployment=deployment,
                            prior_state=prior_state,
                            expected_database_uid=expected_database_uid,
                            expected_database_gid=expected_database_gid,
                            fault=fault,
                        )
                    finally:
                        _close_raw_files(opened)

                remaining = _open_invalid_fresh_sources(
                    directory_descriptor,
                    candidate.name,
                    expected_database_uid=expected_database_uid,
                    expected_database_gid=expected_database_gid,
                )
                try:
                    _verify_invalid_fresh_sources_against_receipt(
                        directory_descriptor,
                        candidate.name,
                        remaining,
                        receipt,
                        expected_database_uid=expected_database_uid,
                        expected_database_gid=expected_database_gid,
                    )
                    for role in sorted(remaining):
                        source = remaining[role]
                        _require_raw_file_unchanged(directory_descriptor, source)
                        os.unlink(source.name, dir_fd=directory_descriptor)
                        if fault is not None:
                            fault(f"invalid-fresh-source-removed:{role}")
                    if _invalid_fresh_roles_present(
                        directory_descriptor, candidate.name
                    ):
                        raise fail(
                            "DATABASE_INVALID_FRESH_QUARANTINE_SOURCE_CHANGED",
                            "Invalid fresh database quarantine source changed.",
                        )
                    os.fsync(directory_descriptor)
                    if fault is not None:
                        fault("invalid-fresh-source-directory-synced")
                finally:
                    _close_raw_files(remaining)
                return receipt
            finally:
                os.close(directory_descriptor)
    except DeploymentError:
        raise
    except Exception:
        raise fail(
            "DATABASE_INVALID_FRESH_QUARANTINE_FAILED",
            "Invalid fresh database could not be quarantined.",
            recoverable=True,
        ) from None


def backup_database(
    layout: DeploymentLayout,
    *,
    task_identity: str,
    expected_deployment_identity: str,
    expected_operator_state_identity: str,
    expected_database_uid: int,
    expected_database_gid: int,
    schema_provider: SchemaAdmissionProvider | None = None,
    fault: DatabaseBackupFault | None = None,
) -> DatabaseBackupReceipt:
    """Consume one trusted write-stop witness and persist an offline backup.

    The operation is idempotent for a task identity.  A crash after witness
    consumption is resumed from the consumed witness; a committed backup is
    re-verified and returned without being overwritten.
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
    opened = _open_database(
        layout.database,
        expected_owner_uid=expected_database_uid,
        expected_owner_gid=expected_database_gid,
    )
    try:
        _require_quiescent_sidecars(opened)
        admission = _admit_schema(
            opened,
            task=task,
            deployment=deployment,
            operator_state=operator_state,
            schema_provider=schema_provider,
        )
        witness, witness_state = _load_witness_for_backup(
            layout,
            task=task,
            deployment=deployment,
            operator_state=operator_state,
            database_stat=opened.initial_stat,
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
        with _backup_operation_lock(layout, task):
            current_witness, current_state = _load_witness_for_backup(
                layout,
                task=task,
                deployment=deployment,
                operator_state=operator_state,
                database_stat=opened.initial_stat,
            )
            if current_witness != witness or current_state != witness_state:
                raise fail(
                    "DATABASE_WRITE_STOP_WITNESS_CHANGED",
                    "Database write-stop witness changed during backup admission.",
                    recoverable=True,
                )
            _require_quiescent_sidecars(opened)
            _require_unchanged(opened)
            backup_identity = canonical_identity(
                {
                    "task_identity": task,
                    "deployment_identity": deployment,
                    "operator_state_identity": operator_state,
                    "witness_identity": witness.identity,
                    "consume_once_identity": witness.consume_once_identity,
                    "schema_admission_identity": admission.identity,
                    "source": inspection.to_dict(),
                },
                scheme=DATABASE_BACKUP_IDENTITY_SCHEME,
            )
            committed = layout.database_backups / backup_identity
            if committed.exists() or committed.is_symlink():
                return _verify_committed_backup(
                    committed,
                    backup_identity=backup_identity,
                    task=task,
                    deployment=deployment,
                    operator_state=operator_state,
                    witness=witness,
                    admission=admission,
                    source=inspection,
                )
            if witness_state == "ready":
                _consume_witness(layout, task)
                if fault is not None:
                    fault("witness-consumed")
            return _create_backup_slot(
                layout,
                opened=opened,
                backup_identity=backup_identity,
                task=task,
                deployment=deployment,
                operator_state=operator_state,
                witness=witness,
                admission=admission,
                source=inspection,
                fault=fault,
            )
    except DeploymentError:
        raise
    except Exception:
        raise fail(
            "DATABASE_BACKUP_FAILED",
            "Database backup could not be completed.",
            recoverable=True,
        ) from None
    finally:
        opened.close()


def _admit_schema(
    opened: _OpenedDatabase,
    *,
    task: str,
    deployment: str,
    operator_state: str,
    schema_provider: SchemaAdmissionProvider | None,
) -> DatabaseSchemaAdmission:
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
    return admission


def _load_witness_for_backup(
    layout: DeploymentLayout,
    *,
    task: str,
    deployment: str,
    operator_state: str,
    database_stat: os.stat_result,
) -> tuple[WriteStopWitness, str]:
    ready = write_stop_witness_path(layout, task)
    consumed = consumed_write_stop_witness_path(layout, task)
    ready_present = _path_present_no_follow(ready)
    consumed_present = _path_present_no_follow(consumed)
    if ready_present == consumed_present:
        raise fail(
            "DATABASE_WRITE_STOP_WITNESS_INVALID",
            "Database write-stop witness is invalid.",
        )
    selected = ready if ready_present else consumed
    witness = _read_write_stop_witness(
        selected,
        expected_task_identity=task,
        expected_deployment_identity=deployment,
        expected_operator_state_identity=operator_state,
        expected_database_path_identity=database_path_identity(layout.database),
        expected_database_stat=database_stat,
    )
    return witness, "ready" if ready_present else "consumed"


def _consume_witness(layout: DeploymentLayout, task: str) -> None:
    ready = write_stop_witness_path(layout, task)
    consumed = consumed_write_stop_witness_path(layout, task)
    directory_flags = (
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    descriptor = -1
    try:
        descriptor = os.open(ready.parent, directory_flags)
        parent = os.fstat(descriptor)
        if (
            parent.st_uid != _ROOT_UID
            or parent.st_gid != _ROOT_GID
            or stat.S_IMODE(parent.st_mode) != _WITNESS_DIRECTORY_MODE
        ):
            raise OSError
        try:
            os.stat(consumed.name, dir_fd=descriptor, follow_symlinks=False)
        except FileNotFoundError:
            pass
        else:
            raise OSError
        os.rename(
            ready.name,
            consumed.name,
            src_dir_fd=descriptor,
            dst_dir_fd=descriptor,
        )
        os.fsync(descriptor)
    except OSError:
        raise fail(
            "DATABASE_WRITE_STOP_WITNESS_CONSUME_FAILED",
            "Database write-stop witness could not be consumed.",
            recoverable=True,
        ) from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _validate_owner_ids(uid: int, gid: int, *, code: str) -> None:
    if any(
        not isinstance(value, int) or isinstance(value, bool) or value < 0
        for value in (uid, gid)
    ):
        raise fail(code, "Database ownership is invalid.")


@contextmanager
def _database_mutation_lock(layout: DeploymentLayout) -> Iterator[None]:
    transactions = _ensure_operator_private_directory(
        layout, layout.operator_transactions, mode=0o700
    )
    path = transactions / "database-mutation.lock"
    flags = (
        os.O_RDWR
        | os.O_CREAT
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    try:
        descriptor = os.open(path, flags, 0o600)
    except OSError:
        raise fail(
            "DATABASE_MUTATION_BUSY",
            "Database mutation boundary is unavailable.",
            recoverable=True,
        ) from None
    try:
        observed = os.fstat(descriptor)
        if (
            not stat.S_ISREG(observed.st_mode)
            or observed.st_nlink != 1
            or observed.st_uid != _ROOT_UID
            or observed.st_gid != _ROOT_GID
            or stat.S_IMODE(observed.st_mode) != 0o600
        ):
            raise fail(
                "DATABASE_MUTATION_BOUNDARY_UNTRUSTED",
                "Database mutation boundary is not trusted.",
            )
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise fail(
                "DATABASE_MUTATION_BUSY",
                "Another database mutation is in progress.",
                recoverable=True,
            ) from None
        yield
    finally:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        finally:
            os.close(descriptor)


def _open_live_database_directory(
    layout: DeploymentLayout,
    *,
    expected_database_uid: int,
    expected_database_gid: int,
) -> int:
    selected = _absolute_path(
        layout.database.parent,
        existing=True,
        code="DATABASE_PATH_UNSAFE",
    )
    flags = (
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    descriptor = -1
    try:
        before = selected.lstat()
        descriptor = os.open(selected, flags)
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISDIR(opened.st_mode)
            or stat.S_ISLNK(before.st_mode)
            or _inode(before) != _inode(opened)
            or opened.st_uid != expected_database_uid
            or opened.st_gid != expected_database_gid
            or stat.S_IMODE(opened.st_mode) != 0o2770
        ):
            raise OSError
        return descriptor
    except OSError:
        if descriptor >= 0:
            os.close(descriptor)
        raise fail(
            "DATABASE_PATH_UNSAFE", "Database path boundary is unsafe."
        ) from None


def _entry_present_at(directory_descriptor: int, name: str) -> bool:
    try:
        os.stat(name, dir_fd=directory_descriptor, follow_symlinks=False)
    except FileNotFoundError:
        return False
    except OSError:
        raise fail(
            "DATABASE_PATH_UNSAFE", "Database path boundary is unsafe."
        ) from None
    return True


def _rename_noreplace_at(
    directory_descriptor: int,
    source_name: str,
    destination_name: str,
    *,
    conflict_code: str,
) -> None:
    """Atomically publish a live coordinate without overwriting a race winner."""

    try:
        libc = ctypes.CDLL(None, use_errno=True)
        renameat2 = libc.renameat2
        renameat2.argtypes = (
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint,
        )
        renameat2.restype = ctypes.c_int
        result = renameat2(
            directory_descriptor,
            os.fsencode(source_name),
            directory_descriptor,
            os.fsencode(destination_name),
            _RENAME_NOREPLACE,
        )
    except (AttributeError, OSError):
        raise fail(
            "DATABASE_FRESH_PUBLISH_FAILED",
            "Fresh database could not be published.",
            recoverable=True,
        ) from None
    if result == 0:
        return
    observed_errno = ctypes.get_errno()
    if observed_errno in {errno.EEXIST, errno.ENOTEMPTY}:
        raise fail(conflict_code, "Fresh database publication is invalid.")
    raise OSError(observed_errno, "renameat2")


def _require_quiescent_sidecar_names(
    directory_descriptor: int,
    database_name: str,
    *,
    expected_owner_uid: int,
    expected_owner_gid: int,
) -> None:
    for suffix in ("-wal", "-shm", "-journal"):
        try:
            observed = os.stat(
                f"{database_name}{suffix}",
                dir_fd=directory_descriptor,
                follow_symlinks=False,
            )
        except FileNotFoundError:
            continue
        except OSError:
            raise fail(
                "DATABASE_SIDECAR_NOT_QUIESCENT",
                "Database sidecar state is not quiescent.",
            ) from None
        if (
            not stat.S_ISREG(observed.st_mode)
            or observed.st_nlink != 1
            or observed.st_uid != expected_owner_uid
            or observed.st_gid != expected_owner_gid
            or stat.S_IMODE(observed.st_mode) != 0o660
            or observed.st_size != 0
        ):
            raise fail(
                "DATABASE_SIDECAR_NOT_QUIESCENT",
                "Database sidecar state is not quiescent.",
            )


def _verify_database_target(
    inspection: DatabaseInspection,
    *,
    expected_identity: str,
    expected_heads: tuple[str, ...],
    code: str,
) -> None:
    if (
        inspection.schema_heads != expected_heads
        or database_content_identity(inspection) != expected_identity
    ):
        raise fail(code, "Database identity or schema is invalid.")


def _verify_backup_copy_target(
    inspection: DatabaseInspection,
    target: DatabaseInspection,
    heads: tuple[str, ...],
    *,
    code: str,
) -> None:
    if (
        inspection.schema_heads != heads
        or inspection.schema_heads != target.schema_heads
        or inspection.sha256 != target.sha256
        or inspection.size_bytes != target.size_bytes
    ):
        raise fail(code, "Database restore target is invalid.")


def _file_evidence(opened: _OpenedDatabase) -> DatabaseFileEvidence:
    before = os.fstat(opened.descriptor)
    sha256, size = _hash_descriptor(opened.descriptor)
    after = os.fstat(opened.descriptor)
    if _file_witness(before) != _file_witness(after):
        raise fail("DATABASE_CHANGED", "Database changed during inspection.")
    return DatabaseFileEvidence.create(
        size_bytes=size,
        sha256=sha256,
        device=after.st_dev,
        inode=after.st_ino,
    )


def _opened_database_matches_target(
    opened: _OpenedDatabase,
    target: DatabaseInspection,
    heads: tuple[str, ...],
) -> bool:
    try:
        inspection = _inspect_descriptor(opened)
    except DeploymentError as error:
        if error.issue.code in {"DATABASE_INTEGRITY_FAILED", "DATABASE_INVALID"}:
            return False
        raise
    return (
        inspection.schema_heads == heads
        and inspection.sha256 == target.sha256
        and inspection.size_bytes == target.size_bytes
    )


def _open_restore_backup(
    layout: DeploymentLayout,
    *,
    backup_identity: str,
    task: str,
    deployment: str,
    prior_state: str,
    source_identity: str,
    schema_heads: tuple[str, ...],
) -> tuple[_OpenedDatabase, DatabaseBackupReceipt, DatabaseInspection]:
    operator = layout.data_root / "operator"
    _require_or_create_root_directory(layout.data_root, create=False, mode=None)
    _require_or_create_operator_root_directory(operator, create=False)
    _require_or_create_root_directory(layout.database_backups, create=False, mode=0o700)
    slot = layout.database_backups / backup_identity
    try:
        observed = slot.lstat()
    except OSError:
        raise fail("DATABASE_BACKUP_INVALID", "Database backup is invalid.") from None
    if (
        not stat.S_ISDIR(observed.st_mode)
        or stat.S_ISLNK(observed.st_mode)
        or observed.st_uid != _ROOT_UID
        or observed.st_gid != _ROOT_GID
        or stat.S_IMODE(observed.st_mode) != 0o500
    ):
        raise fail("DATABASE_BACKUP_INVALID", "Database backup is invalid.")
    opened = _open_database(
        slot / "platform.db",
        expected_owner_uid=_ROOT_UID,
        expected_owner_gid=_ROOT_GID,
    )
    try:
        if stat.S_IMODE(opened.initial_stat.st_mode) != 0o400 or set(
            os.listdir(opened.directory_descriptor)
        ) != {"platform.db", "receipt.json"}:
            raise fail("DATABASE_BACKUP_INVALID", "Database backup is invalid.")
        receipt = _read_backup_receipt_at(opened.directory_descriptor)
        _require_quiescent_sidecars(opened)
        inspection = _inspect_descriptor(opened)
        _require_unchanged(opened)
        expected_backup_identity = canonical_identity(
            {
                "task_identity": receipt.task_identity,
                "deployment_identity": receipt.deployment_identity,
                "operator_state_identity": receipt.operator_state_identity,
                "witness_identity": receipt.witness_identity,
                "consume_once_identity": receipt.consume_once_identity,
                "schema_admission_identity": receipt.schema_admission_identity,
                "source": receipt.source.to_dict(),
            },
            scheme=DATABASE_BACKUP_IDENTITY_SCHEME,
        )
        if (
            receipt.backup_identity != backup_identity
            or expected_backup_identity != backup_identity
            or receipt.task_identity != task
            or receipt.deployment_identity != deployment
            or receipt.operator_state_identity != prior_state
            or database_content_identity(receipt.source) != source_identity
            or receipt.source.schema_heads != schema_heads
            or receipt.backup != inspection
            or inspection.schema_heads != schema_heads
            or inspection.sha256 != receipt.source.sha256
            or inspection.size_bytes != receipt.source.size_bytes
        ):
            raise fail("DATABASE_BACKUP_INVALID", "Database backup is invalid.")
        return opened, receipt, inspection
    except Exception:
        opened.close()
        raise


def _read_backup_receipt_at(directory_descriptor: int) -> DatabaseBackupReceipt:
    return _read_canonical_receipt_at(
        directory_descriptor,
        "receipt.json",
        parser=DatabaseBackupReceipt.from_dict,
        code="DATABASE_BACKUP_INVALID",
        message="Database backup is invalid.",
    )


def _recovery_destination(
    layout: DeploymentLayout,
    *,
    operation: str,
    task: str,
    identity: str,
) -> Path:
    if operation not in {"restore", "fresh-quarantine", "fresh-invalid"}:
        raise fail(
            "DATABASE_RECOVERY_BOUNDARY_UNTRUSTED",
            "Database recovery boundary is not trusted.",
        )
    parent = _ensure_operator_private_directory(
        layout,
        layout.data_root / "operator" / "database-recovery" / operation / task,
        mode=0o700,
    )
    return parent / identity


def _preserve_restore_evidence(
    layout: DeploymentLayout,
    *,
    opened: _OpenedDatabase,
    backup_identity: str,
    backup_receipt: DatabaseBackupReceipt,
    task: str,
    deployment: str,
    prior_state: str,
    source_identity: str,
    schema_heads: tuple[str, ...],
    target: DatabaseInspection,
    fault: DatabaseRecoveryFault | None,
) -> DatabaseRestoreReceipt:
    destination = _recovery_destination(
        layout, operation="restore", task=task, identity=backup_identity
    )
    if _path_exists(destination):
        receipt = _load_restore_evidence(
            layout,
            backup_identity=backup_identity,
            task=task,
            deployment=deployment,
            prior_state=prior_state,
            source_identity=source_identity,
            schema_heads=schema_heads,
            backup_receipt=backup_receipt,
            target=target,
        )
        if receipt is None:
            raise fail(
                "DATABASE_RESTORE_EVIDENCE_INVALID",
                "Database recovery evidence is invalid.",
            )
        return receipt
    source = _file_evidence(opened)
    partial = destination.parent / f".partial-{backup_identity}-{uuid.uuid4().hex}"
    try:
        partial.mkdir(mode=0o700)
        os.chown(partial, _ROOT_UID, _ROOT_GID)
        directory_descriptor = _open_root_directory(partial, expected_mode=0o700)
        try:
            _copy_database_descriptor_at(
                opened.descriptor,
                directory_descriptor,
                "platform.db",
                owner_uid=_ROOT_UID,
                owner_gid=_ROOT_GID,
                mode=0o400,
                code="DATABASE_RESTORE_FAILED",
            )
            evidence_opened = _open_database(
                partial / "platform.db",
                expected_owner_uid=_ROOT_UID,
                expected_owner_gid=_ROOT_GID,
            )
            try:
                evidence = _file_evidence(evidence_opened)
            finally:
                evidence_opened.close()
            if (
                evidence.sha256 != source.sha256
                or evidence.size_bytes != source.size_bytes
            ):
                raise fail(
                    "DATABASE_RESTORE_EVIDENCE_INVALID",
                    "Database recovery evidence is invalid.",
                )
            if fault is not None:
                fault("restore-evidence-file-synced")
            receipt = DatabaseRestoreReceipt.create(
                backup_identity=backup_identity,
                backup_receipt_identity=backup_receipt.identity,
                task_identity=task,
                deployment_identity=deployment,
                prior_state_identity=prior_state,
                source_identity=source_identity,
                schema_heads=schema_heads,
                failed_current=source,
                evidence=evidence,
                target=target,
            )
            _write_root_file_at(
                directory_descriptor,
                "receipt.json",
                canonical_json_bytes(receipt.to_dict()),
                mode=0o444,
            )
            os.fchmod(directory_descriptor, 0o500)
            os.fsync(directory_descriptor)
            if fault is not None:
                fault("restore-evidence-ready")
        finally:
            os.close(directory_descriptor)
        _commit_recovery_slot(partial, destination)
        if fault is not None:
            fault("restore-evidence-committed")
        loaded = _load_restore_evidence(
            layout,
            backup_identity=backup_identity,
            task=task,
            deployment=deployment,
            prior_state=prior_state,
            source_identity=source_identity,
            schema_heads=schema_heads,
            backup_receipt=backup_receipt,
            target=target,
        )
        if loaded != receipt:
            raise fail(
                "DATABASE_RESTORE_EVIDENCE_INVALID",
                "Database recovery evidence is invalid.",
            )
        return receipt
    except DeploymentError:
        raise
    except OSError:
        raise fail(
            "DATABASE_RESTORE_FAILED",
            "Database restore could not be completed.",
            recoverable=True,
        ) from None


def _load_restore_evidence(
    layout: DeploymentLayout,
    *,
    backup_identity: str,
    task: str,
    deployment: str,
    prior_state: str,
    source_identity: str,
    schema_heads: tuple[str, ...],
    backup_receipt: DatabaseBackupReceipt,
    target: DatabaseInspection,
) -> DatabaseRestoreReceipt | None:
    destination = _recovery_destination(
        layout, operation="restore", task=task, identity=backup_identity
    )
    if not _path_exists(destination):
        return None
    opened = _open_recovery_evidence(destination)
    try:
        receipt = _read_canonical_receipt_at(
            opened.directory_descriptor,
            "receipt.json",
            parser=DatabaseRestoreReceipt.from_dict,
            code="DATABASE_RESTORE_EVIDENCE_INVALID",
            message="Database recovery evidence is invalid.",
        )
        evidence = _file_evidence(opened)
        if (
            receipt.backup_identity != backup_identity
            or receipt.backup_receipt_identity != backup_receipt.identity
            or receipt.task_identity != task
            or receipt.deployment_identity != deployment
            or receipt.prior_state_identity != prior_state
            or receipt.source_identity != source_identity
            or receipt.schema_heads != schema_heads
            or receipt.target != target
            or receipt.evidence != evidence
            or receipt.failed_current.sha256 != evidence.sha256
            or receipt.failed_current.size_bytes != evidence.size_bytes
        ):
            raise fail(
                "DATABASE_RESTORE_EVIDENCE_INVALID",
                "Database recovery evidence is invalid.",
            )
        return receipt
    finally:
        opened.close()


def _preserve_quarantine_evidence(
    layout: DeploymentLayout,
    *,
    opened: _OpenedDatabase,
    task: str,
    deployment: str,
    prior_state: str,
    candidate_identity: str,
    schema_heads: tuple[str, ...],
    source_coordinate: str,
    fault: DatabaseRecoveryFault | None,
) -> DatabaseQuarantineReceipt:
    destination = _recovery_destination(
        layout,
        operation="fresh-quarantine",
        task=task,
        identity=candidate_identity,
    )
    source = _file_evidence(opened)
    partial = destination.parent / f".partial-{candidate_identity}-{uuid.uuid4().hex}"
    try:
        partial.mkdir(mode=0o700)
        os.chown(partial, _ROOT_UID, _ROOT_GID)
        directory_descriptor = _open_root_directory(partial, expected_mode=0o700)
        try:
            _copy_database_descriptor_at(
                opened.descriptor,
                directory_descriptor,
                "platform.db",
                owner_uid=_ROOT_UID,
                owner_gid=_ROOT_GID,
                mode=0o400,
                code="DATABASE_QUARANTINE_FAILED",
            )
            evidence_opened = _open_database(
                partial / "platform.db",
                expected_owner_uid=_ROOT_UID,
                expected_owner_gid=_ROOT_GID,
            )
            try:
                evidence = _file_evidence(evidence_opened)
            finally:
                evidence_opened.close()
            if (
                evidence.sha256 != source.sha256
                or evidence.size_bytes != source.size_bytes
            ):
                raise fail(
                    "DATABASE_QUARANTINE_EVIDENCE_INVALID",
                    "Database quarantine evidence is invalid.",
                )
            if fault is not None:
                fault("quarantine-evidence-file-synced")
            receipt = DatabaseQuarantineReceipt.create(
                task_identity=task,
                deployment_identity=deployment,
                prior_state_identity=prior_state,
                candidate_identity=candidate_identity,
                schema_heads=schema_heads,
                source_coordinate=source_coordinate,
                source=source,
                evidence=evidence,
            )
            _write_root_file_at(
                directory_descriptor,
                "receipt.json",
                canonical_json_bytes(receipt.to_dict()),
                mode=0o444,
            )
            os.fchmod(directory_descriptor, 0o500)
            os.fsync(directory_descriptor)
            if fault is not None:
                fault("quarantine-evidence-ready")
        finally:
            os.close(directory_descriptor)
        _commit_recovery_slot(partial, destination)
        if fault is not None:
            fault("quarantine-evidence-committed")
        loaded = _load_quarantine_evidence(
            layout,
            task=task,
            deployment=deployment,
            prior_state=prior_state,
            candidate_identity=candidate_identity,
            schema_heads=schema_heads,
        )
        if loaded != receipt:
            raise fail(
                "DATABASE_QUARANTINE_EVIDENCE_INVALID",
                "Database quarantine evidence is invalid.",
            )
        return receipt
    except DeploymentError:
        raise
    except OSError:
        raise fail(
            "DATABASE_QUARANTINE_FAILED",
            "Fresh database could not be quarantined.",
            recoverable=True,
        ) from None


def _load_quarantine_evidence(
    layout: DeploymentLayout,
    *,
    task: str,
    deployment: str,
    prior_state: str,
    candidate_identity: str,
    schema_heads: tuple[str, ...],
) -> DatabaseQuarantineReceipt | None:
    destination = _recovery_destination(
        layout,
        operation="fresh-quarantine",
        task=task,
        identity=candidate_identity,
    )
    if not _path_exists(destination):
        return None
    opened = _open_recovery_evidence(destination)
    try:
        receipt = _read_canonical_receipt_at(
            opened.directory_descriptor,
            "receipt.json",
            parser=DatabaseQuarantineReceipt.from_dict,
            code="DATABASE_QUARANTINE_EVIDENCE_INVALID",
            message="Database quarantine evidence is invalid.",
        )
        evidence = _file_evidence(opened)
        if (
            receipt.task_identity != task
            or receipt.deployment_identity != deployment
            or receipt.prior_state_identity != prior_state
            or receipt.candidate_identity != candidate_identity
            or receipt.schema_heads != schema_heads
            or receipt.evidence != evidence
            or receipt.source.sha256 != evidence.sha256
            or receipt.source.size_bytes != evidence.size_bytes
        ):
            raise fail(
                "DATABASE_QUARANTINE_EVIDENCE_INVALID",
                "Database quarantine evidence is invalid.",
            )
        return receipt
    finally:
        opened.close()


def _invalid_fresh_source_names(candidate_name: str) -> dict[str, str]:
    return {
        "database": candidate_name,
        "journal": f"{candidate_name}-journal",
        "shm": f"{candidate_name}-shm",
        "wal": f"{candidate_name}-wal",
    }


def _invalid_fresh_roles_present(
    directory_descriptor: int, candidate_name: str
) -> tuple[str, ...]:
    return tuple(
        role
        for role, name in _invalid_fresh_source_names(candidate_name).items()
        if _entry_present_at(directory_descriptor, name)
    )


def _open_invalid_fresh_sources(
    directory_descriptor: int,
    candidate_name: str,
    *,
    expected_database_uid: int,
    expected_database_gid: int,
) -> dict[str, _OpenedRawDatabaseFile]:
    opened: dict[str, _OpenedRawDatabaseFile] = {}
    flags = (
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOFOLLOW", 0)
        | getattr(os, "O_NONBLOCK", 0)
    )
    try:
        for role, name in _invalid_fresh_source_names(candidate_name).items():
            if not _entry_present_at(directory_descriptor, name):
                continue
            descriptor = -1
            try:
                descriptor = os.open(name, flags, dir_fd=directory_descriptor)
                observed = os.fstat(descriptor)
                path_observed = os.stat(
                    name, dir_fd=directory_descriptor, follow_symlinks=False
                )
                if (
                    not stat.S_ISREG(observed.st_mode)
                    or observed.st_nlink != 1
                    or observed.st_uid != expected_database_uid
                    or observed.st_gid != expected_database_gid
                    or stat.S_IMODE(observed.st_mode) != 0o660
                    or not 0 <= observed.st_size <= _MAX_DATABASE_BYTES
                    or _inode(observed) != _inode(path_observed)
                ):
                    raise OSError
                opened[role] = _OpenedRawDatabaseFile(
                    role=role,
                    name=name,
                    descriptor=descriptor,
                    initial_stat=observed,
                )
                descriptor = -1
            finally:
                if descriptor >= 0:
                    os.close(descriptor)
        return opened
    except (OSError, DeploymentError):
        _close_raw_files(opened)
        raise fail(
            "DATABASE_INVALID_FRESH_QUARANTINE_SOURCE_INVALID",
            "Invalid fresh database quarantine source is invalid.",
        ) from None


def _close_raw_files(opened: Mapping[str, _OpenedRawDatabaseFile]) -> None:
    for item in opened.values():
        item.close()


def _raw_file_evidence(
    opened: _OpenedRawDatabaseFile,
    *,
    code: str,
    message: str,
) -> DatabaseRawFileEvidence:
    try:
        before = os.fstat(opened.descriptor)
        os.lseek(opened.descriptor, 0, os.SEEK_SET)
        digest = hashlib.sha256()
        size = 0
        while True:
            chunk = os.read(opened.descriptor, _READ_CHUNK_BYTES)
            if not chunk:
                break
            size += len(chunk)
            if size > _MAX_DATABASE_BYTES:
                raise OSError
            digest.update(chunk)
        after = os.fstat(opened.descriptor)
        if _file_witness(before) != _file_witness(after) or size != after.st_size:
            raise OSError
        return DatabaseRawFileEvidence.create(
            size_bytes=size,
            sha256=digest.hexdigest(),
            device=after.st_dev,
            inode=after.st_ino,
            uid=after.st_uid,
            gid=after.st_gid,
            mode=stat.S_IMODE(after.st_mode),
            nlink=after.st_nlink,
            mtime_ns=after.st_mtime_ns,
            ctime_ns=after.st_ctime_ns,
        )
    except (OSError, DeploymentError):
        raise fail(code, message) from None


def _require_raw_file_unchanged(
    directory_descriptor: int, opened: _OpenedRawDatabaseFile
) -> None:
    try:
        descriptor_observed = os.fstat(opened.descriptor)
        path_observed = os.stat(
            opened.name,
            dir_fd=directory_descriptor,
            follow_symlinks=False,
        )
        if _file_witness(opened.initial_stat) != _file_witness(
            descriptor_observed
        ) or _inode(descriptor_observed) != _inode(path_observed):
            raise OSError
    except OSError:
        raise fail(
            "DATABASE_INVALID_FRESH_QUARANTINE_SOURCE_CHANGED",
            "Invalid fresh database quarantine source changed.",
        ) from None


def _verify_invalid_fresh_sources_against_receipt(
    directory_descriptor: int,
    candidate_name: str,
    opened: Mapping[str, _OpenedRawDatabaseFile],
    receipt: InvalidFreshDatabaseQuarantineReceipt,
    *,
    expected_database_uid: int,
    expected_database_gid: int,
) -> None:
    expected = {item.role: item for item in receipt.files}
    if any(
        item.source.uid != expected_database_uid
        or item.source.gid != expected_database_gid
        or item.source.mode != 0o660
        or item.source.nlink != 1
        for item in receipt.files
    ) or not set(opened).issubset(expected):
        raise fail(
            "DATABASE_INVALID_FRESH_QUARANTINE_EVIDENCE_INVALID",
            "Invalid fresh database quarantine evidence is invalid.",
        )
    for role, source in opened.items():
        observed = _raw_file_evidence(
            source,
            code="DATABASE_INVALID_FRESH_QUARANTINE_SOURCE_CHANGED",
            message="Invalid fresh database quarantine source changed.",
        )
        if observed != expected[role].source:
            raise fail(
                "DATABASE_INVALID_FRESH_QUARANTINE_SOURCE_CHANGED",
                "Invalid fresh database quarantine source changed.",
            )
        _require_raw_file_unchanged(directory_descriptor, source)
    if set(_invalid_fresh_roles_present(directory_descriptor, candidate_name)) != set(
        opened
    ):
        raise fail(
            "DATABASE_INVALID_FRESH_QUARANTINE_SOURCE_CHANGED",
            "Invalid fresh database quarantine source changed.",
        )


def _preserve_invalid_fresh_quarantine_evidence(
    layout: DeploymentLayout,
    *,
    opened: Mapping[str, _OpenedRawDatabaseFile],
    source_directory_descriptor: int,
    request_identity: str,
    task: str,
    deployment: str,
    prior_state: str,
    expected_database_uid: int,
    expected_database_gid: int,
    fault: DatabaseRecoveryFault | None,
) -> InvalidFreshDatabaseQuarantineReceipt:
    destination = _recovery_destination(
        layout,
        operation="fresh-invalid",
        task=task,
        identity=request_identity,
    )
    partial = destination.parent / f".partial-{request_identity}-{uuid.uuid4().hex}"
    try:
        partial.mkdir(mode=0o700)
        os.chown(partial, _ROOT_UID, _ROOT_GID)
        directory_descriptor = _open_root_directory(partial, expected_mode=0o700)
        try:
            files: list[InvalidFreshDatabaseFileEvidence] = []
            for role in sorted(opened):
                source = opened[role]
                source_evidence = _raw_file_evidence(
                    source,
                    code="DATABASE_INVALID_FRESH_QUARANTINE_SOURCE_CHANGED",
                    message="Invalid fresh database quarantine source changed.",
                )
                evidence_name = _INVALID_FRESH_FILE_NAMES[role]
                _copy_database_descriptor_at(
                    source.descriptor,
                    directory_descriptor,
                    evidence_name,
                    owner_uid=_ROOT_UID,
                    owner_gid=_ROOT_GID,
                    mode=0o400,
                    code="DATABASE_INVALID_FRESH_QUARANTINE_FAILED",
                )
                evidence_opened = _open_raw_file_at(
                    directory_descriptor,
                    evidence_name,
                    role=role,
                    expected_uid=_ROOT_UID,
                    expected_gid=_ROOT_GID,
                    expected_mode=0o400,
                    code="DATABASE_INVALID_FRESH_QUARANTINE_EVIDENCE_INVALID",
                    message="Invalid fresh database quarantine evidence is invalid.",
                )
                try:
                    evidence = _raw_file_evidence(
                        evidence_opened,
                        code="DATABASE_INVALID_FRESH_QUARANTINE_EVIDENCE_INVALID",
                        message=(
                            "Invalid fresh database quarantine evidence is invalid."
                        ),
                    )
                finally:
                    evidence_opened.close()
                if (
                    source_evidence.sha256 != evidence.sha256
                    or source_evidence.size_bytes != evidence.size_bytes
                ):
                    raise fail(
                        "DATABASE_INVALID_FRESH_QUARANTINE_EVIDENCE_INVALID",
                        "Invalid fresh database quarantine evidence is invalid.",
                    )
                files.append(
                    InvalidFreshDatabaseFileEvidence(
                        role=role,
                        source=source_evidence,
                        evidence=evidence,
                    )
                )
                if fault is not None:
                    fault(f"invalid-fresh-evidence-file-synced:{role}")
            receipt = InvalidFreshDatabaseQuarantineReceipt.create(
                request_identity=request_identity,
                task_identity=task,
                deployment_identity=deployment,
                prior_state_identity=prior_state,
                files=files,
            )
            _write_root_file_at(
                directory_descriptor,
                "receipt.json",
                canonical_json_bytes(receipt.to_dict()),
                mode=0o444,
            )
            if fault is not None:
                fault("invalid-fresh-evidence-receipt-synced")
            os.fchmod(directory_descriptor, 0o500)
            os.fsync(directory_descriptor)
            if fault is not None:
                fault("invalid-fresh-evidence-directory-synced")
        finally:
            os.close(directory_descriptor)
        for source in opened.values():
            _require_raw_file_unchanged(source_directory_descriptor, source)
        _commit_recovery_slot(partial, destination)
        if fault is not None:
            fault("invalid-fresh-evidence-committed")
        loaded = _load_invalid_fresh_quarantine_evidence(
            layout,
            request_identity=request_identity,
            task=task,
            deployment=deployment,
            prior_state=prior_state,
            expected_database_uid=expected_database_uid,
            expected_database_gid=expected_database_gid,
        )
        if loaded != receipt:
            raise fail(
                "DATABASE_INVALID_FRESH_QUARANTINE_EVIDENCE_INVALID",
                "Invalid fresh database quarantine evidence is invalid.",
            )
        return receipt
    except DeploymentError:
        raise
    except OSError:
        raise fail(
            "DATABASE_INVALID_FRESH_QUARANTINE_FAILED",
            "Invalid fresh database could not be quarantined.",
            recoverable=True,
        ) from None


def _load_invalid_fresh_quarantine_evidence(
    layout: DeploymentLayout,
    *,
    request_identity: str,
    task: str,
    deployment: str,
    prior_state: str,
    expected_database_uid: int,
    expected_database_gid: int,
) -> InvalidFreshDatabaseQuarantineReceipt | None:
    destination = _recovery_destination(
        layout,
        operation="fresh-invalid",
        task=task,
        identity=request_identity,
    )
    if not _path_exists(destination):
        return None
    directory_descriptor = _open_root_directory(destination, expected_mode=0o500)
    try:
        receipt = _read_canonical_receipt_at(
            directory_descriptor,
            "receipt.json",
            parser=InvalidFreshDatabaseQuarantineReceipt.from_dict,
            code="DATABASE_INVALID_FRESH_QUARANTINE_EVIDENCE_INVALID",
            message="Invalid fresh database quarantine evidence is invalid.",
        )
        if (
            receipt.request_identity != request_identity
            or receipt.task_identity != task
            or receipt.deployment_identity != deployment
            or receipt.prior_state_identity != prior_state
            or any(
                item.source.uid != expected_database_uid
                or item.source.gid != expected_database_gid
                for item in receipt.files
            )
        ):
            raise fail(
                "DATABASE_INVALID_FRESH_QUARANTINE_EVIDENCE_INVALID",
                "Invalid fresh database quarantine evidence is invalid.",
            )
        expected_names = {"receipt.json"} | {
            _INVALID_FRESH_FILE_NAMES[item.role] for item in receipt.files
        }
        if set(os.listdir(directory_descriptor)) != expected_names:
            raise fail(
                "DATABASE_INVALID_FRESH_QUARANTINE_EVIDENCE_INVALID",
                "Invalid fresh database quarantine evidence is invalid.",
            )
        for item in receipt.files:
            evidence_opened = _open_raw_file_at(
                directory_descriptor,
                _INVALID_FRESH_FILE_NAMES[item.role],
                role=item.role,
                expected_uid=_ROOT_UID,
                expected_gid=_ROOT_GID,
                expected_mode=0o400,
                code="DATABASE_INVALID_FRESH_QUARANTINE_EVIDENCE_INVALID",
                message="Invalid fresh database quarantine evidence is invalid.",
            )
            try:
                evidence = _raw_file_evidence(
                    evidence_opened,
                    code="DATABASE_INVALID_FRESH_QUARANTINE_EVIDENCE_INVALID",
                    message="Invalid fresh database quarantine evidence is invalid.",
                )
            finally:
                evidence_opened.close()
            if evidence != item.evidence:
                raise fail(
                    "DATABASE_INVALID_FRESH_QUARANTINE_EVIDENCE_INVALID",
                    "Invalid fresh database quarantine evidence is invalid.",
                )
        return receipt
    finally:
        os.close(directory_descriptor)


def _open_raw_file_at(
    directory_descriptor: int,
    name: str,
    *,
    role: str,
    expected_uid: int,
    expected_gid: int,
    expected_mode: int,
    code: str,
    message: str,
) -> _OpenedRawDatabaseFile:
    descriptor = -1
    flags = (
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOFOLLOW", 0)
        | getattr(os, "O_NONBLOCK", 0)
    )
    try:
        descriptor = os.open(name, flags, dir_fd=directory_descriptor)
        observed = os.fstat(descriptor)
        path_observed = os.stat(
            name, dir_fd=directory_descriptor, follow_symlinks=False
        )
        if (
            not stat.S_ISREG(observed.st_mode)
            or observed.st_nlink != 1
            or observed.st_uid != expected_uid
            or observed.st_gid != expected_gid
            or stat.S_IMODE(observed.st_mode) != expected_mode
            or not 0 <= observed.st_size <= _MAX_DATABASE_BYTES
            or _inode(observed) != _inode(path_observed)
        ):
            raise OSError
        return _OpenedRawDatabaseFile(role, name, descriptor, observed)
    except OSError:
        if descriptor >= 0:
            os.close(descriptor)
        raise fail(code, message) from None


def _open_recovery_evidence(path: Path) -> _OpenedDatabase:
    try:
        observed = path.lstat()
    except OSError:
        raise fail(
            "DATABASE_RECOVERY_EVIDENCE_INVALID",
            "Database recovery evidence is invalid.",
        ) from None
    if (
        not stat.S_ISDIR(observed.st_mode)
        or stat.S_ISLNK(observed.st_mode)
        or observed.st_uid != _ROOT_UID
        or observed.st_gid != _ROOT_GID
        or stat.S_IMODE(observed.st_mode) != 0o500
    ):
        raise fail(
            "DATABASE_RECOVERY_EVIDENCE_INVALID",
            "Database recovery evidence is invalid.",
        )
    opened = _open_database(
        path / "platform.db",
        expected_owner_uid=_ROOT_UID,
        expected_owner_gid=_ROOT_GID,
    )
    if stat.S_IMODE(opened.initial_stat.st_mode) != 0o400 or set(
        os.listdir(opened.directory_descriptor)
    ) != {"platform.db", "receipt.json"}:
        opened.close()
        raise fail(
            "DATABASE_RECOVERY_EVIDENCE_INVALID",
            "Database recovery evidence is invalid.",
        )
    return opened


def _open_root_directory(path: Path, *, expected_mode: int) -> int:
    flags = (
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    descriptor = -1
    try:
        before = path.lstat()
        descriptor = os.open(path, flags)
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISDIR(opened.st_mode)
            or stat.S_ISLNK(before.st_mode)
            or _inode(before) != _inode(opened)
            or opened.st_uid != _ROOT_UID
            or opened.st_gid != _ROOT_GID
            or stat.S_IMODE(opened.st_mode) != expected_mode
        ):
            raise OSError
        return descriptor
    except OSError:
        if descriptor >= 0:
            os.close(descriptor)
        raise fail(
            "DATABASE_RECOVERY_BOUNDARY_UNTRUSTED",
            "Database recovery boundary is not trusted.",
        ) from None


def _commit_recovery_slot(partial: Path, destination: Path) -> None:
    parent_descriptor = _open_root_directory(destination.parent, expected_mode=0o700)
    try:
        if _entry_present_at(parent_descriptor, destination.name):
            raise fail(
                "DATABASE_RECOVERY_EVIDENCE_INVALID",
                "Database recovery evidence is invalid.",
            )
        os.rename(
            partial.name,
            destination.name,
            src_dir_fd=parent_descriptor,
            dst_dir_fd=parent_descriptor,
        )
        os.fsync(parent_descriptor)
    finally:
        os.close(parent_descriptor)


def _path_exists(path: Path) -> bool:
    try:
        path.lstat()
    except FileNotFoundError:
        return False
    except OSError:
        raise fail(
            "DATABASE_RECOVERY_EVIDENCE_INVALID",
            "Database recovery evidence is invalid.",
        ) from None
    return True


def _copy_database_descriptor_at(
    source: int,
    directory_descriptor: int,
    name: str,
    *,
    owner_uid: int,
    owner_gid: int,
    mode: int,
    code: str,
) -> None:
    flags = (
        os.O_WRONLY
        | os.O_CREAT
        | os.O_EXCL
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    descriptor = -1
    try:
        descriptor = os.open(name, flags, 0o600, dir_fd=directory_descriptor)
        os.fchown(descriptor, owner_uid, owner_gid)
        os.lseek(source, 0, os.SEEK_SET)
        while True:
            chunk = os.read(source, _READ_CHUNK_BYTES)
            if not chunk:
                break
            _write_all(descriptor, chunk)
        os.fchmod(descriptor, mode)
        os.fsync(descriptor)
    except OSError:
        raise fail(
            code, "Database bytes could not be copied.", recoverable=True
        ) from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _write_root_file_at(
    directory_descriptor: int,
    name: str,
    content: bytes,
    *,
    mode: int,
) -> None:
    flags = (
        os.O_WRONLY
        | os.O_CREAT
        | os.O_EXCL
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    descriptor = os.open(name, flags, 0o600, dir_fd=directory_descriptor)
    try:
        os.fchown(descriptor, _ROOT_UID, _ROOT_GID)
        _write_all(descriptor, content)
        os.fchmod(descriptor, mode)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _read_canonical_receipt_at(
    directory_descriptor: int,
    name: str,
    *,
    parser: Callable[[object], Any],
    code: str,
    message: str,
):
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = -1
    try:
        descriptor = os.open(name, flags, dir_fd=directory_descriptor)
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_uid != _ROOT_UID
            or before.st_gid != _ROOT_GID
            or stat.S_IMODE(before.st_mode) != 0o444
            or not 0 < before.st_size <= _MAX_RECEIPT_BYTES
        ):
            raise OSError
        content = _read_bounded(descriptor, _MAX_RECEIPT_BYTES)
        after = os.fstat(descriptor)
        path_after = os.stat(name, dir_fd=directory_descriptor, follow_symlinks=False)
        if (
            _file_witness(before) != _file_witness(after)
            or _inode(after) != _inode(path_after)
            or len(content) != before.st_size
        ):
            raise OSError
    except OSError:
        raise fail(code, message) from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    try:
        raw = json.loads(content, object_pairs_hook=_unique_object)
        receipt = parser(raw)
    except (DeploymentError, UnicodeDecodeError, ValueError, json.JSONDecodeError):
        raise fail(code, message) from None
    if canonical_json_bytes(receipt.to_dict()) != content:
        raise fail(code, message)
    return receipt


@contextmanager
def _backup_operation_lock(layout: DeploymentLayout, task: str) -> Iterator[None]:
    transactions = _ensure_operator_private_directory(
        layout, layout.operator_transactions, mode=0o700
    )
    path = transactions / f"database-backup-{task}.lock"
    flags = (
        os.O_RDWR
        | os.O_CREAT
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    try:
        descriptor = os.open(path, flags, 0o600)
    except OSError:
        raise fail(
            "DATABASE_BACKUP_BUSY",
            "Database backup operation is unavailable.",
            recoverable=True,
        ) from None
    try:
        observed = os.fstat(descriptor)
        if (
            not stat.S_ISREG(observed.st_mode)
            or observed.st_nlink != 1
            or observed.st_uid != _ROOT_UID
            or observed.st_gid != _ROOT_GID
            or stat.S_IMODE(observed.st_mode) != 0o600
        ):
            raise fail(
                "DATABASE_BACKUP_BOUNDARY_UNTRUSTED",
                "Database backup boundary is not trusted.",
            )
        try:
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise fail(
                "DATABASE_BACKUP_BUSY",
                "Another database backup operation is in progress.",
                recoverable=True,
            ) from None
        yield
    finally:
        try:
            fcntl.flock(descriptor, fcntl.LOCK_UN)
        finally:
            os.close(descriptor)


def _create_backup_slot(
    layout: DeploymentLayout,
    *,
    opened: _OpenedDatabase,
    backup_identity: str,
    task: str,
    deployment: str,
    operator_state: str,
    witness: WriteStopWitness,
    admission: DatabaseSchemaAdmission,
    source: DatabaseInspection,
    fault: DatabaseBackupFault | None,
) -> DatabaseBackupReceipt:
    root = _ensure_operator_private_directory(
        layout, layout.database_backups, mode=0o700
    )
    destination = root / backup_identity
    partial = root / f".partial-{backup_identity}-{uuid.uuid4().hex}"
    try:
        partial.mkdir(mode=0o700)
        os.chown(partial, _ROOT_UID, _ROOT_GID)
        backup_path = partial / "platform.db"
        _copy_database_descriptor(opened.descriptor, backup_path)
        if fault is not None:
            fault("backup-file-synced")
        backup = inspect_database(
            backup_path,
            expected_owner_uid=_ROOT_UID,
            expected_owner_gid=_ROOT_GID,
        )
        if (
            backup.sha256 != source.sha256
            or backup.size_bytes != source.size_bytes
            or backup.schema_heads != source.schema_heads
        ):
            raise fail("DATABASE_BACKUP_INVALID", "Database backup is invalid.")
        receipt = DatabaseBackupReceipt.create(
            backup_identity=backup_identity,
            task_identity=task,
            deployment_identity=deployment,
            operator_state_identity=operator_state,
            witness_identity=witness.identity,
            consume_once_identity=witness.consume_once_identity,
            schema_admission_identity=admission.identity,
            source=source,
            backup=backup,
        )
        _write_root_file(
            partial / "receipt.json",
            canonical_json_bytes(receipt.to_dict()),
            mode=0o444,
        )
        os.chmod(backup_path, 0o400)
        _fsync_path(backup_path)
        os.chmod(partial, 0o500)
        _fsync_directory_descriptor(partial)
        if fault is not None:
            fault("backup-slot-ready")
        os.replace(partial, destination)
        _fsync_directory_descriptor(root)
        if fault is not None:
            fault("backup-committed")
        return _verify_committed_backup(
            destination,
            backup_identity=backup_identity,
            task=task,
            deployment=deployment,
            operator_state=operator_state,
            witness=witness,
            admission=admission,
            source=source,
        )
    except DeploymentError:
        raise
    except OSError:
        raise fail(
            "DATABASE_BACKUP_FAILED",
            "Database backup could not be completed.",
            recoverable=True,
        ) from None


def _verify_committed_backup(
    root: Path,
    *,
    backup_identity: str,
    task: str,
    deployment: str,
    operator_state: str,
    witness: WriteStopWitness,
    admission: DatabaseSchemaAdmission,
    source: DatabaseInspection,
) -> DatabaseBackupReceipt:
    try:
        observed = root.lstat()
    except OSError:
        raise fail("DATABASE_BACKUP_INVALID", "Database backup is invalid.") from None
    if (
        not stat.S_ISDIR(observed.st_mode)
        or stat.S_ISLNK(observed.st_mode)
        or observed.st_uid != _ROOT_UID
        or observed.st_gid != _ROOT_GID
        or stat.S_IMODE(observed.st_mode) != 0o500
    ):
        raise fail("DATABASE_BACKUP_INVALID", "Database backup is invalid.")
    receipt = _read_backup_receipt(root / "receipt.json")
    backup = inspect_database(
        root / "platform.db",
        expected_owner_uid=_ROOT_UID,
        expected_owner_gid=_ROOT_GID,
    )
    if (
        receipt.backup_identity != backup_identity
        or receipt.task_identity != task
        or receipt.deployment_identity != deployment
        or receipt.operator_state_identity != operator_state
        or receipt.witness_identity != witness.identity
        or receipt.consume_once_identity != witness.consume_once_identity
        or receipt.schema_admission_identity != admission.identity
        or receipt.source != source
        or receipt.backup != backup
        or backup.sha256 != source.sha256
        or backup.size_bytes != source.size_bytes
        or backup.schema_heads != source.schema_heads
    ):
        raise fail("DATABASE_BACKUP_INVALID", "Database backup is invalid.")
    return receipt


def _read_backup_receipt(path: Path) -> DatabaseBackupReceipt:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_uid != _ROOT_UID
            or before.st_gid != _ROOT_GID
            or stat.S_IMODE(before.st_mode) != 0o444
            or not 0 < before.st_size <= _MAX_RECEIPT_BYTES
        ):
            raise OSError
        content = _read_bounded(descriptor, _MAX_RECEIPT_BYTES)
        after = os.fstat(descriptor)
    except OSError:
        raise fail("DATABASE_BACKUP_INVALID", "Database backup is invalid.") from None
    finally:
        if "descriptor" in locals():
            os.close(descriptor)
    if _file_witness(before) != _file_witness(after) or len(content) != before.st_size:
        raise fail("DATABASE_BACKUP_INVALID", "Database backup is invalid.")
    try:
        raw = json.loads(content, object_pairs_hook=_unique_object)
        receipt = DatabaseBackupReceipt.from_dict(raw)
    except (DeploymentError, UnicodeDecodeError, ValueError, json.JSONDecodeError):
        raise fail("DATABASE_BACKUP_INVALID", "Database backup is invalid.") from None
    if canonical_json_bytes(receipt.to_dict()) != content:
        raise fail("DATABASE_BACKUP_INVALID", "Database backup is invalid.")
    return receipt


def _copy_database_descriptor(source: int, destination: Path) -> None:
    flags = (
        os.O_WRONLY
        | os.O_CREAT
        | os.O_EXCL
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    descriptor = -1
    try:
        descriptor = os.open(destination, flags, 0o600)
        os.fchown(descriptor, _ROOT_UID, _ROOT_GID)
        os.lseek(source, 0, os.SEEK_SET)
        while True:
            chunk = os.read(source, _READ_CHUNK_BYTES)
            if not chunk:
                break
            _write_all(descriptor, chunk)
        os.fsync(descriptor)
    except OSError:
        raise fail(
            "DATABASE_BACKUP_FAILED",
            "Database backup could not be completed.",
            recoverable=True,
        ) from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _ensure_operator_private_directory(
    layout: DeploymentLayout, path: Path, *, mode: int
) -> Path:
    operator = layout.data_root / "operator"
    if path != operator and operator not in path.parents:
        raise fail(
            "DATABASE_BACKUP_BOUNDARY_UNTRUSTED",
            "Database backup boundary is not trusted.",
        )
    _require_or_create_root_directory(layout.data_root, create=False, mode=None)
    _require_or_create_operator_root_directory(operator, create=True)
    current = operator
    for part in path.relative_to(operator).parts:
        current = current / part
        _require_or_create_root_directory(current, create=True, mode=mode)
    return path


def _require_or_create_operator_root_directory(path: Path, *, create: bool) -> None:
    """Require the fixed root-owned, group-traversable operator boundary.

    The bootstrap owns the group assignment so administrators in its dedicated
    group can traverse to fixed helpers.  The group identity is deliberately
    not duplicated here: exact mode 0710 grants traversal but no mutation, and
    every database backup/recovery descendant remains root:root private.
    """

    try:
        observed = path.lstat()
    except FileNotFoundError:
        if not create:
            raise fail(
                "DATABASE_BACKUP_BOUNDARY_UNTRUSTED",
                "Database backup boundary is not trusted.",
            ) from None
        try:
            path.mkdir(mode=0o710)
            os.chown(path, _ROOT_UID, _ROOT_GID)
            observed = path.lstat()
        except OSError:
            raise fail(
                "DATABASE_BACKUP_BOUNDARY_UNTRUSTED",
                "Database backup boundary is not trusted.",
            ) from None
    except OSError:
        raise fail(
            "DATABASE_BACKUP_BOUNDARY_UNTRUSTED",
            "Database backup boundary is not trusted.",
        ) from None
    if (
        not stat.S_ISDIR(observed.st_mode)
        or stat.S_ISLNK(observed.st_mode)
        or observed.st_uid != _ROOT_UID
        or stat.S_IMODE(observed.st_mode) != 0o710
    ):
        raise fail(
            "DATABASE_BACKUP_BOUNDARY_UNTRUSTED",
            "Database backup boundary is not trusted.",
        )


def _require_or_create_root_directory(
    path: Path, *, create: bool, mode: int | None
) -> None:
    try:
        observed = path.lstat()
    except FileNotFoundError:
        if not create or mode is None:
            raise fail(
                "DATABASE_BACKUP_BOUNDARY_UNTRUSTED",
                "Database backup boundary is not trusted.",
            ) from None
        try:
            path.mkdir(mode=mode)
            os.chown(path, _ROOT_UID, _ROOT_GID)
            observed = path.lstat()
        except OSError:
            raise fail(
                "DATABASE_BACKUP_BOUNDARY_UNTRUSTED",
                "Database backup boundary is not trusted.",
            ) from None
    except OSError:
        raise fail(
            "DATABASE_BACKUP_BOUNDARY_UNTRUSTED",
            "Database backup boundary is not trusted.",
        ) from None
    if (
        not stat.S_ISDIR(observed.st_mode)
        or stat.S_ISLNK(observed.st_mode)
        or observed.st_uid != _ROOT_UID
        or observed.st_gid != _ROOT_GID
        or stat.S_IMODE(observed.st_mode) & 0o022
        or (mode is not None and stat.S_IMODE(observed.st_mode) != mode)
    ):
        raise fail(
            "DATABASE_BACKUP_BOUNDARY_UNTRUSTED",
            "Database backup boundary is not trusted.",
        )


def _write_root_file(path: Path, content: bytes, *, mode: int) -> None:
    flags = (
        os.O_WRONLY
        | os.O_CREAT
        | os.O_EXCL
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    descriptor = os.open(path, flags, 0o600)
    try:
        os.fchown(descriptor, _ROOT_UID, _ROOT_GID)
        _write_all(descriptor, content)
        os.fchmod(descriptor, mode)
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _write_all(descriptor: int, content: bytes) -> None:
    offset = 0
    while offset < len(content):
        written = os.write(descriptor, content[offset:])
        if written <= 0:
            raise OSError
        offset += written


def _fsync_path(path: Path) -> None:
    descriptor = os.open(
        path,
        os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _fsync_directory_descriptor(path: Path) -> None:
    descriptor = os.open(
        path,
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _path_present_no_follow(path: Path) -> bool:
    try:
        path.lstat()
    except FileNotFoundError:
        return False
    except OSError:
        raise fail(
            "DATABASE_WRITE_STOP_WITNESS_INVALID",
            "Database write-stop witness is invalid.",
        ) from None
    return True


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


def _open_database(
    path: Path,
    *,
    expected_owner_uid: int,
    expected_owner_gid: int,
) -> _OpenedDatabase:
    if any(
        not isinstance(value, int) or isinstance(value, bool) or value < 0
        for value in (expected_owner_uid, expected_owner_gid)
    ):
        raise fail("DATABASE_PATH_UNSAFE", "Database path boundary is unsafe.")
    selected = _absolute_path(path, existing=True, code="DATABASE_PATH_INVALID")
    parent_stat = selected.parent.lstat()
    parent_mode = stat.S_IMODE(parent_stat.st_mode)
    owner_only_parent = parent_mode & 0o022 == 0
    shared_parent = parent_mode == 0o2770
    if (
        not stat.S_ISDIR(parent_stat.st_mode)
        or stat.S_ISLNK(parent_stat.st_mode)
        or parent_stat.st_uid != expected_owner_uid
        or parent_stat.st_gid != expected_owner_gid
        or not (owner_only_parent or shared_parent)
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
        file_mode = stat.S_IMODE(observed.st_mode)
        owner_only_file = file_mode & 0o022 == 0
        shared_file = file_mode == 0o660
        shared_group_write = shared_parent and shared_file
        if (
            not stat.S_ISREG(observed.st_mode)
            or observed.st_nlink != 1
            or observed.st_uid != expected_owner_uid
            or observed.st_gid != expected_owner_gid
            or not (shared_group_write or (owner_only_parent and owner_only_file))
            or not 0 < observed.st_size <= _MAX_DATABASE_BYTES
        ):
            raise fail("DATABASE_INVALID", "Database is invalid.")
        return _OpenedDatabase(
            selected,
            directory_descriptor,
            descriptor,
            observed,
            shared_group_write,
        )
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
    for suffix in ("-wal", "-shm", "-journal"):
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
                "DATABASE_SIDECAR_NOT_QUIESCENT",
                "Database sidecar state is not quiescent.",
            ) from None
        mode = stat.S_IMODE(observed.st_mode)
        trusted_mode = mode == 0o660 if opened.shared_group_write else mode & 0o022 == 0
        if (
            not stat.S_ISREG(observed.st_mode)
            or observed.st_nlink != 1
            or observed.st_uid != opened.initial_stat.st_uid
            or observed.st_gid != opened.initial_stat.st_gid
            or not trusted_mode
            or observed.st_size != 0
        ):
            raise fail(
                "DATABASE_SIDECAR_NOT_QUIESCENT",
                "Database sidecar state is not quiescent.",
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
    "DATABASE_BACKUP_RECEIPT_SCHEMA",
    "DATABASE_INSPECTION_SCHEMA",
    "DATABASE_INVALID_FRESH_QUARANTINE_RECEIPT_SCHEMA",
    "DATABASE_QUARANTINE_RECEIPT_SCHEMA",
    "DATABASE_RAW_FILE_EVIDENCE_SCHEMA",
    "DATABASE_RESTORE_RECEIPT_SCHEMA",
    "DATABASE_SCHEMA_ADMISSION_SCHEMA",
    "DATABASE_WRITER_UNITS",
    "WRITE_STOP_ISSUER",
    "WRITE_STOP_WITNESS_IDENTITY_SCHEME",
    "WRITE_STOP_WITNESS_SCHEMA",
    "DatabaseBackupReceipt",
    "DatabaseFileEvidence",
    "DatabaseInspection",
    "DatabaseRawFileEvidence",
    "DatabaseQuarantineReceipt",
    "DatabaseRestoreReceipt",
    "DatabaseSchemaAdmission",
    "SchemaAdmissionProvider",
    "StoppedWriter",
    "InvalidFreshDatabaseFileEvidence",
    "InvalidFreshDatabaseQuarantineReceipt",
    "WriteStopWitness",
    "UnavailableSchemaAdmissionProvider",
    "backup_database",
    "create_write_stop_witness",
    "database_content_identity",
    "consumed_write_stop_witness_path",
    "database_path_identity",
    "fresh_database_candidate_path",
    "inspect_database",
    "publish_fresh_database",
    "quarantine_fresh_database",
    "quarantine_invalid_fresh_database",
    "restore_database_backup",
    "write_stop_witness_path",
]
