"""Rejection-path coverage for deployment evidence document parsers."""

from __future__ import annotations

import os

import pytest

from encode_pipeline.deployment import database as database_module
from encode_pipeline.deployment.canonical import canonical_identity
from encode_pipeline.deployment.database import (
    DATABASE_FILE_EVIDENCE_IDENTITY_SCHEME,
    DATABASE_QUARANTINE_RECEIPT_IDENTITY_SCHEME,
    DATABASE_RAW_FILE_EVIDENCE_IDENTITY_SCHEME,
    DATABASE_RESTORE_RECEIPT_IDENTITY_SCHEME,
    DATABASE_WRITER_UNITS,
    WRITE_STOP_WITNESS_IDENTITY_SCHEME,
    WRITE_STOP_ISSUER,
    WRITE_STOP_WITNESS_SCHEMA,
    DatabaseBackupReceipt,
    DatabaseFileEvidence,
    DatabaseInspection,
    DatabaseQuarantineReceipt,
    DatabaseRawFileEvidence,
    DatabaseRestoreReceipt,
    DatabaseSchemaAdmission,
    InvalidFreshDatabaseFileEvidence,
    InvalidFreshDatabaseQuarantineReceipt,
    StoppedWriter,
    WriteStopWitness,
)
from encode_pipeline.deployment.errors import DeploymentError

_ID = f"sha256-{'a' * 64}"
_ID_B = f"sha256-{'b' * 64}"
_ID_C = f"sha256-{'c' * 64}"
_ID_D = f"sha256-{'d' * 64}"
_TASK = f"task-{'a' * 32}"
_TS = "2026-01-02T03:04:05+00:00"
_INSPECTION = DatabaseInspection(
    schema_heads=("20260809_13",),
    size_bytes=128,
    sha256="a" * 64,
    device=10,
    inode=20,
)
_INSPECTION_B = DatabaseInspection(
    schema_heads=("20260809_13",),
    size_bytes=256,
    sha256="b" * 64,
    device=11,
    inode=21,
)


def _reject(parser, document: object, code: str) -> None:
    with pytest.raises(DeploymentError) as captured:
        parser(document)
    assert captured.value.issue.code == code


def _signed(document: dict[str, object], scheme: str) -> dict[str, object]:
    value = dict(document)
    value["identity"] = canonical_identity(
        {key: item for key, item in value.items() if key != "identity"},
        scheme=scheme,
    )
    return value


def _witness_document() -> dict[str, object]:
    value: dict[str, object] = {
        "schema_version": WRITE_STOP_WITNESS_SCHEMA,
        "issuer": WRITE_STOP_ISSUER,
        "consumption_state": "ready",
        "task_identity": _TASK,
        "deployment_identity": _ID,
        "operator_state_identity": _ID_B,
        "consume_once_identity": _ID_C,
        "database_path_identity": _ID_D,
        "database_device": 10,
        "database_inode": 20,
        "created_at": _TS,
        "writers": [
            {
                "unit": unit,
                "state": "stopped",
                "service_identity": f"sha256-{str(index) * 64}",
            }
            for index, unit in enumerate(DATABASE_WRITER_UNITS)
        ],
    }
    value["identity"] = canonical_identity(
        value, scheme=WRITE_STOP_WITNESS_IDENTITY_SCHEME
    )
    return value


def _file_evidence() -> DatabaseFileEvidence:
    document = _signed(
        {
            "schema_version": "helixweave-database-file-evidence-v1",
            "size_bytes": 128,
            "sha256": "a" * 64,
            "device": 10,
            "inode": 20,
        },
        DATABASE_FILE_EVIDENCE_IDENTITY_SCHEME,
    )
    return DatabaseFileEvidence.from_dict(document)


def _raw_file_evidence(*, mode: int = 0o660, uid: int = 0, gid: int = 0):
    document = _signed(
        {
            "schema_version": "helixweave-database-raw-file-evidence-v1",
            "size_bytes": 64,
            "sha256": "a" * 64,
            "device": 10,
            "inode": 20,
            "uid": uid,
            "gid": gid,
            "mode": mode,
            "nlink": 1,
            "mtime_ns": 1,
            "ctime_ns": 2,
        },
        DATABASE_RAW_FILE_EVIDENCE_IDENTITY_SCHEME,
    )
    return DatabaseRawFileEvidence.from_dict(document)


def test_stopped_writer_rejects_malformed_documents() -> None:
    code = "DATABASE_WRITE_STOP_WITNESS_INVALID"
    valid = StoppedWriter(
        unit=DATABASE_WRITER_UNITS[0], state="stopped", service_identity=_ID
    ).to_dict()
    _reject(StoppedWriter.from_dict, None, code)
    _reject(StoppedWriter.from_dict, ["unit"], code)
    _reject(StoppedWriter.from_dict, {1: "unit"}, code)
    _reject(StoppedWriter.from_dict, {}, code)
    _reject(StoppedWriter.from_dict, {**valid, "extra": 1}, code)
    _reject(StoppedWriter.from_dict, {**valid, "unit": "other.service"}, code)
    _reject(StoppedWriter.from_dict, {**valid, "state": "running"}, code)
    _reject(StoppedWriter.from_dict, {**valid, "service_identity": 1}, code)
    _reject(StoppedWriter.from_dict, {**valid, "service_identity": "sha256-zz"}, code)
    assert StoppedWriter.from_dict(valid).unit == DATABASE_WRITER_UNITS[0]


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: ["schema_version"],
        lambda d: {key: d[key] for key in d if key != "writers"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "issuer": "other"},
        lambda d: {**d, "consumption_state": "consumed"},
        lambda d: {**d, "writers": "stopped"},
        lambda d: {**d, "writers": [{}]},
        lambda d: {**d, "writers": list(reversed(d["writers"]))},
        lambda d: {
            **d,
            "writers": [{**writer, "service_identity": _ID} for writer in d["writers"]],
        },
        lambda d: {**d, "task_identity": "not-a-task"},
        lambda d: {**d, "task_identity": 1},
        lambda d: {**d, "database_device": True},
        lambda d: {**d, "database_device": 0},
        lambda d: {**d, "database_device": 2**63},
        lambda d: {**d, "database_inode": "20"},
        lambda d: {**d, "created_at": "not-a-date"},
        lambda d: {**d, "created_at": "2026-01-02T03:04:05"},
        lambda d: {**d, "created_at": "x" * 65},
        lambda d: {**d, "identity": _ID},
    ],
)
def test_write_stop_witness_rejects_malformed_documents(mutation) -> None:
    document = mutation(_witness_document())
    _reject(WriteStopWitness.from_dict, document, "DATABASE_WRITE_STOP_WITNESS_INVALID")


def test_write_stop_witness_accepts_valid_document() -> None:
    witness = WriteStopWitness.from_dict(_witness_document())
    assert witness.task_identity == _TASK
    assert len(witness.writers) == len(DATABASE_WRITER_UNITS)


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "schema_heads"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "integrity": "broken"},
        lambda d: {**d, "wal_state": "present"},
        lambda d: {**d, "sha256": 1},
        lambda d: {**d, "sha256": "z" * 64},
        lambda d: {**d, "size_bytes": True},
        lambda d: {**d, "size_bytes": 0},
        lambda d: {**d, "size_bytes": 2**63},
        lambda d: {**d, "device": True},
        lambda d: {**d, "device": 0},
        lambda d: {**d, "inode": 0},
        lambda d: {**d, "schema_heads": "20260809_13"},
        lambda d: {**d, "schema_heads": ["!!bad"]},
    ],
)
def test_database_inspection_rejects_malformed_documents(mutation) -> None:
    document = mutation(_INSPECTION.to_dict())
    _reject(DatabaseInspection.from_dict, document, "DATABASE_BACKUP_RECEIPT_INVALID")


def _schema_admission_document() -> dict[str, object]:
    return DatabaseSchemaAdmission.create(
        provider_identity=_ID,
        task_identity=_TASK,
        deployment_identity=_ID_B,
        operator_state_identity=_ID_C,
        database_device=10,
        database_inode=20,
        accepted_schema_heads=("20260809_13",),
    ).to_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "provider_identity"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "identity": _ID_D},
        lambda d: {**d, "provider_identity": "bad"},
        lambda d: {**d, "task_identity": "bad"},
        lambda d: {**d, "deployment_identity": 1},
        lambda d: {**d, "operator_state_identity": "bad"},
        lambda d: {**d, "database_device": True},
        lambda d: {**d, "database_inode": -1},
        lambda d: {**d, "accepted_schema_heads": ["!!bad"]},
        lambda d: {**d, "accepted_schema_heads": "20260809_13"},
    ],
)
def test_schema_admission_rejects_malformed_documents(mutation) -> None:
    document = mutation(_schema_admission_document())
    _reject(
        DatabaseSchemaAdmission.from_dict, document, "DATABASE_SCHEMA_ADMISSION_INVALID"
    )


def _backup_receipt_document() -> dict[str, object]:
    return DatabaseBackupReceipt.create(
        backup_identity=_ID,
        task_identity=_TASK,
        deployment_identity=_ID_B,
        operator_state_identity=_ID_C,
        witness_identity=_ID_D,
        consume_once_identity=_ID,
        schema_admission_identity=_ID_B,
        source=_INSPECTION,
        backup=_INSPECTION_B,
    ).to_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "source"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "identity": _ID_D},
        lambda d: {**d, "backup_identity": "bad"},
        lambda d: {**d, "task_identity": "bad"},
        lambda d: {**d, "deployment_identity": 1},
        lambda d: {**d, "operator_state_identity": "bad"},
        lambda d: {**d, "witness_identity": "bad"},
        lambda d: {**d, "consume_once_identity": "bad"},
        lambda d: {**d, "schema_admission_identity": "bad"},
        lambda d: {**d, "source": None},
        lambda d: {**d, "source": {**d["source"], "integrity": "broken"}},
        lambda d: {**d, "backup": {**d["backup"], "sha256": "z" * 64}},
    ],
)
def test_backup_receipt_rejects_malformed_documents(mutation) -> None:
    document = mutation(_backup_receipt_document())
    _reject(
        DatabaseBackupReceipt.from_dict, document, "DATABASE_BACKUP_RECEIPT_INVALID"
    )


def _file_evidence_document() -> dict[str, object]:
    return _file_evidence().to_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "sha256"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "sha256": 1},
        lambda d: {**d, "sha256": "z" * 64},
        lambda d: {**d, "size_bytes": True},
        lambda d: {**d, "size_bytes": 0},
        lambda d: {**d, "size_bytes": 2**63},
        lambda d: {**d, "identity": _ID_D},
        lambda d: {**d, "device": True},
        lambda d: {**d, "inode": 0},
    ],
)
def test_file_evidence_rejects_malformed_documents(mutation) -> None:
    document = mutation(_file_evidence_document())
    _reject(
        DatabaseFileEvidence.from_dict, document, "DATABASE_RECOVERY_RECEIPT_INVALID"
    )


def _raw_file_evidence_document() -> dict[str, object]:
    return _raw_file_evidence().to_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "mtime_ns"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "sha256": "z" * 64},
        lambda d: {**d, "size_bytes": True},
        lambda d: {**d, "size_bytes": -1},
        lambda d: {**d, "mode": True},
        lambda d: {**d, "mode": 0o644},
        lambda d: {**d, "nlink": True},
        lambda d: {**d, "nlink": 2},
        lambda d: {**d, "identity": _ID_D},
        lambda d: {**d, "device": 0},
        lambda d: {**d, "inode": True},
        lambda d: {**d, "uid": -1},
        lambda d: {**d, "gid": True},
        lambda d: {**d, "mtime_ns": -1},
        lambda d: {**d, "ctime_ns": "2"},
    ],
)
def test_raw_file_evidence_rejects_malformed_documents(mutation) -> None:
    document = mutation(_raw_file_evidence_document())
    _reject(
        DatabaseRawFileEvidence.from_dict,
        document,
        "DATABASE_INVALID_FRESH_QUARANTINE_RECEIPT_INVALID",
    )


def _invalid_fresh_file_document(monkeypatch: pytest.MonkeyPatch) -> dict[str, object]:
    monkeypatch.setattr(database_module, "_ROOT_UID", os.getuid())
    monkeypatch.setattr(database_module, "_ROOT_GID", os.getgid())
    source = _raw_file_evidence(mode=0o660)
    evidence = _raw_file_evidence(mode=0o400, uid=os.getuid(), gid=os.getgid())
    return InvalidFreshDatabaseFileEvidence(
        role="database", source=source, evidence=evidence
    ).to_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "role"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "role": 1},
        lambda d: {**d, "role": "other"},
        lambda d: {**d, "source": None},
        lambda d: {**d, "evidence": {**d["evidence"], "mode": 0o660}},
        lambda d: {**d, "evidence": {**d["evidence"], "uid": 999999}},
        lambda d: {**d, "evidence": {**d["evidence"], "gid": 999999}},
        lambda d: {**d, "source": {**d["source"], "sha256": "b" * 64}},
        lambda d: {**d, "source": {**d["source"], "size_bytes": 65}},
    ],
)
def test_invalid_fresh_file_evidence_rejects_malformed_documents(
    mutation, monkeypatch: pytest.MonkeyPatch
) -> None:
    document = mutation(_invalid_fresh_file_document(monkeypatch))
    _reject(
        InvalidFreshDatabaseFileEvidence.from_dict,
        document,
        "DATABASE_INVALID_FRESH_QUARANTINE_RECEIPT_INVALID",
    )


def _invalid_fresh_receipt_document(
    monkeypatch: pytest.MonkeyPatch,
) -> dict[str, object]:
    file_document = _invalid_fresh_file_document(monkeypatch)
    file_evidence = InvalidFreshDatabaseFileEvidence.from_dict(file_document)
    value: dict[str, object] = {
        "schema_version": "helixweave-database-invalid-fresh-quarantine-receipt-v1",
        "status": "raw-evidence-preserved",
        "request_identity": _ID,
        "task_identity": _TASK,
        "deployment_identity": _ID_B,
        "prior_state_identity": _ID_C,
        "files": [file_document],
    }
    value["raw_evidence_identity"] = canonical_identity(
        {
            "files": [
                {
                    "role": file_evidence.role,
                    "source": file_evidence.source.to_dict(),
                }
            ]
        },
        scheme="helixweave-database-invalid-fresh-raw-evidence-identity-v1",
    )
    return _signed(
        value, "helixweave-database-invalid-fresh-quarantine-receipt-identity-v1"
    )


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "files"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "status": "other"},
        lambda d: {**d, "files": "database"},
        lambda d: {**d, "files": [{"role": "other", "source": {}, "evidence": {}}]},
        lambda d: {**d, "files": []},
        lambda d: {**d, "identity": _ID_D},
        lambda d: {**d, "request_identity": "bad"},
        lambda d: {**d, "raw_evidence_identity": _ID_D},
        lambda d: {**d, "task_identity": "bad"},
        lambda d: {**d, "deployment_identity": 1},
        lambda d: {**d, "prior_state_identity": "bad"},
    ],
)
def test_invalid_fresh_receipt_rejects_malformed_documents(
    mutation, monkeypatch: pytest.MonkeyPatch
) -> None:
    document = mutation(_invalid_fresh_receipt_document(monkeypatch))
    _reject(
        InvalidFreshDatabaseQuarantineReceipt.from_dict,
        document,
        "DATABASE_INVALID_FRESH_QUARANTINE_RECEIPT_INVALID",
    )


def _restore_receipt_document() -> dict[str, object]:
    receipt = DatabaseRestoreReceipt(
        identity=_ID,
        backup_identity=_ID_B,
        backup_receipt_identity=_ID_C,
        task_identity=_TASK,
        deployment_identity=_ID,
        prior_state_identity=_ID_B,
        source_identity=_ID_C,
        schema_heads=("20260809_13",),
        failed_current=_file_evidence(),
        evidence=_file_evidence(),
        target=_INSPECTION,
    )
    return _signed(receipt.to_dict(), DATABASE_RESTORE_RECEIPT_IDENTITY_SCHEME)


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "target"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "status": "other"},
        lambda d: {**d, "identity": _ID_D},
        lambda d: {**d, "backup_identity": "bad"},
        lambda d: {**d, "backup_receipt_identity": 1},
        lambda d: {**d, "task_identity": "bad"},
        lambda d: {**d, "deployment_identity": 1},
        lambda d: {**d, "prior_state_identity": "bad"},
        lambda d: {**d, "source_identity": "bad"},
        lambda d: {**d, "schema_heads": ["!!bad"]},
        lambda d: {**d, "failed_current": None},
        lambda d: {**d, "evidence": {**d["evidence"], "size_bytes": 0}},
        lambda d: {**d, "target": {**d["target"], "wal_state": "present"}},
    ],
)
def test_restore_receipt_rejects_malformed_documents(mutation) -> None:
    document = mutation(_restore_receipt_document())
    _reject(
        DatabaseRestoreReceipt.from_dict, document, "DATABASE_RESTORE_RECEIPT_INVALID"
    )


def _quarantine_receipt_document() -> dict[str, object]:
    receipt = DatabaseQuarantineReceipt(
        identity=_ID,
        task_identity=_TASK,
        deployment_identity=_ID_B,
        prior_state_identity=_ID_C,
        candidate_identity=_ID,
        schema_heads=("20260809_13",),
        source_coordinate="backup",
        source=_file_evidence(),
        evidence=_file_evidence(),
    )
    return _signed(receipt.to_dict(), DATABASE_QUARANTINE_RECEIPT_IDENTITY_SCHEME)


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "source_coordinate"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "identity": _ID_D},
        lambda d: {**d, "task_identity": "bad"},
        lambda d: {**d, "deployment_identity": 1},
        lambda d: {**d, "prior_state_identity": "bad"},
        lambda d: {**d, "candidate_identity": "bad"},
        lambda d: {**d, "schema_heads": "20260809_13"},
        lambda d: {**d, "source": None},
        lambda d: {**d, "evidence": {**d["evidence"], "inode": 0}},
    ],
)
def test_quarantine_receipt_rejects_malformed_documents(mutation) -> None:
    document = mutation(_quarantine_receipt_document())
    _reject(
        DatabaseQuarantineReceipt.from_dict,
        document,
        "DATABASE_QUARANTINE_RECEIPT_INVALID",
    )
