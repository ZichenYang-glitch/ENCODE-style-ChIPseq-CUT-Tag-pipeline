"""Fixed unprivileged candidate-action protocol used by the root operator.

The root process owns request selection and state mutation.  Candidate code is
only ever run by the fixed systemd oneshot as the ``helixweave`` account and
returns a canonical, identity-bound attestation through a fixed coordinate.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import subprocess
import tempfile
from typing import Protocol

from encode_pipeline.deployment.canonical import (
    canonical_identity,
    canonical_json_bytes,
    without_key,
)
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import (
    BULK_RNASEQ_RUNTIME,
    COMPONENTS,
    ENCODE_RUNTIME,
)


ACTION_REQUEST_SCHEMA = "helixweave-operator-action-request-v1"
ACTION_RECEIPT_SCHEMA = "helixweave-operator-action-receipt-v2"
ACTION_REQUEST_IDENTITY_SCHEME = "helixweave-operator-action-request-identity-v1"
ACTION_RECEIPT_IDENTITY_SCHEME = "helixweave-operator-action-receipt-identity-v2"
DATABASE_PREPARE_REQUEST_SCHEMA = "helixweave-database-prepare-request-v2"
DATABASE_PREPARE_RECEIPT_SCHEMA = "helixweave-database-prepare-receipt-v1"
DATABASE_PREPARE_REQUEST_IDENTITY_SCHEME = (
    "helixweave-database-prepare-request-identity-v2"
)
DATABASE_PREPARE_RECEIPT_IDENTITY_SCHEME = (
    "helixweave-database-prepare-receipt-identity-v1"
)
ENCODE_RUNTIME_PREPARE_REQUEST_SCHEMA = "helixweave-encode-runtime-prepare-request-v1"
ENCODE_RUNTIME_PREPARE_RECEIPT_SCHEMA = "helixweave-encode-runtime-prepare-receipt-v1"
ENCODE_RUNTIME_INVENTORY_SCHEMA = "helixweave-encode-runtime-inventory-v1"
ENCODE_RUNTIME_PREPARE_REQUEST_IDENTITY_SCHEME = (
    "helixweave-encode-runtime-prepare-request-identity-v1"
)
ENCODE_RUNTIME_PREPARE_RECEIPT_IDENTITY_SCHEME = (
    "helixweave-encode-runtime-prepare-receipt-identity-v1"
)
ENCODE_RUNTIME_TREE_IDENTITY_SCHEME = "helixweave-encode-runtime-tree-identity-v1"
BULK_RUNTIME_PREPARE_REQUEST_SCHEMA = "helixweave-bulk-runtime-prepare-request-v2"
BULK_RUNTIME_PREPARE_RECEIPT_SCHEMA = "helixweave-bulk-runtime-prepare-receipt-v2"
BULK_RUNTIME_PREPARE_REQUEST_IDENTITY_SCHEME = (
    "helixweave-bulk-runtime-prepare-request-identity-v2"
)
BULK_RUNTIME_PREPARE_RECEIPT_IDENTITY_SCHEME = (
    "helixweave-bulk-runtime-prepare-receipt-identity-v2"
)
ACTION_UNIT = "helixweave-operator-action.service"
DATABASE_PREPARE_UNIT = "helixweave-db-prepare.service"
ENCODE_RUNTIME_PREPARE_UNIT = "helixweave-encode-runtime-prepare.service"
SYSTEMCTL = Path("/usr/bin/systemctl")
ACTION_PHASES = ("admit", "observe")
ACTION_OPERATIONS = ("activate", "rollback", "observe")
ACTION_STATUSES = ("admitted", "observed")
COMPATIBILITY_STATUSES = ("compatible", "incomplete", "incompatible")
DATABASE_MODES = ("fresh-candidate", "existing-live")
BULK_RUNTIME_PREPARE_OPERATIONS = ("activate", "rollback", "verify")
VERIFICATION_CHECKS = (
    "platform-native",
    "python-runtime",
    "frontend",
    "encode-runtime-native",
    "bulk-runtime-native",
    "database-schema",
    "configuration",
    "permissions",
    "references",
    "redis",
    "worker",
    "docker",
)
READINESS_STATUSES = ("ready", "not-ready", "not-applicable", "unavailable")
READINESS_REASONS = (
    "READY",
    "COMPONENT_NOT_ACTIVE",
    "CONTRACT_INVALID",
    "SCHEMA_INCOMPATIBLE",
    "CONFIGURATION_INVALID",
    "PERMISSION_INVALID",
    "REFERENCE_NOT_READY",
    "SERVICE_STOPPED",
    "SERVICE_IDENTITY_INVALID",
    "REDIS_UNAVAILABLE",
    "DOCKER_UNAVAILABLE",
    "FRONTEND_INVALID",
    "NOT_APPLICABLE",
)

_READINESS_REASON_BY_STATUS = {
    "ready": frozenset({"READY"}),
    "not-ready": frozenset(
        {
            "CONTRACT_INVALID",
            "SCHEMA_INCOMPATIBLE",
            "CONFIGURATION_INVALID",
            "PERMISSION_INVALID",
            "REFERENCE_NOT_READY",
            "SERVICE_STOPPED",
            "SERVICE_IDENTITY_INVALID",
            "REDIS_UNAVAILABLE",
            "DOCKER_UNAVAILABLE",
            "FRONTEND_INVALID",
        }
    ),
    "not-applicable": frozenset({"COMPONENT_NOT_ACTIVE", "NOT_APPLICABLE"}),
    "unavailable": frozenset(
        {
            "CONTRACT_INVALID",
            "SCHEMA_INCOMPATIBLE",
            "CONFIGURATION_INVALID",
            "PERMISSION_INVALID",
            "REFERENCE_NOT_READY",
            "SERVICE_STOPPED",
            "SERVICE_IDENTITY_INVALID",
            "REDIS_UNAVAILABLE",
            "DOCKER_UNAVAILABLE",
            "FRONTEND_INVALID",
        }
    ),
}

_IDENTITY = re.compile(r"^sha256-[0-9a-f]{64}$")
_TASK = re.compile(r"^task-[0-9a-f]{32}$")
_HEAD = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_MAX_DOCUMENT_BYTES = 1024 * 1024
_MAX_RUNTIME_INVENTORY_BYTES = 64 * 1024 * 1024
_MAX_RUNTIME_FILES = 500_000
_MAX_RUNTIME_FILE_BYTES = 8 * 1024 * 1024 * 1024
_MAX_RUNTIME_TOTAL_BYTES = 128 * 1024 * 1024 * 1024
_MAX_BULK_RUNTIME_IMAGES = 4_096
_ACTION_TIMEOUT_SECONDS = 300.0
_ENCODE_RUNTIME_UNIT_TIMEOUT_SECONDS = 1_860.0
_UNIT_QUIESCE_TIMEOUT_SECONDS = 60.0
_SAFE_ENVIRONMENT = {
    "LANG": "C.UTF-8",
    "LC_ALL": "C.UTF-8",
    "PATH": "/usr/sbin:/usr/bin:/sbin:/bin",
}


def _identity(value: object, *, code: str) -> str:
    if not isinstance(value, str) or _IDENTITY.fullmatch(value) is None:
        raise fail(code, "Operator action evidence is invalid.")
    return value


def _optional_identity(value: object, *, code: str) -> str | None:
    return None if value is None else _identity(value, code=code)


def _task(value: object, *, code: str) -> str:
    if not isinstance(value, str) or _TASK.fullmatch(value) is None:
        raise fail(code, "Operator action evidence is invalid.")
    return value


def _heads(value: object, *, code: str, allow_empty: bool = True) -> tuple[str, ...]:
    if not isinstance(value, list):
        raise fail(code, "Operator action evidence is invalid.")
    parsed = tuple(value)
    if (
        (not allow_empty and not parsed)
        or tuple(sorted(set(parsed))) != parsed
        or any(
            not isinstance(item, str) or _HEAD.fullmatch(item) is None
            for item in parsed
        )
    ):
        raise fail(code, "Operator action evidence is invalid.")
    return parsed


@dataclass(frozen=True)
class DeploymentActionRequest:
    identity: str
    phase: str
    operation: str
    component: str
    task_identity: str
    deployment_identity: str
    authority_platform_identity: str
    prior_state_identity: str
    candidate_state_identity: str
    candidate_active: dict[str, str | None]

    @classmethod
    def create(
        cls,
        *,
        phase: str,
        operation: str,
        component: str,
        task_identity: str,
        deployment_identity: str,
        authority_platform_identity: str,
        prior_state_identity: str,
        candidate_state_identity: str,
        candidate_active: dict[str, str | None],
    ) -> "DeploymentActionRequest":
        value: dict[str, object] = {
            "schema_version": ACTION_REQUEST_SCHEMA,
            "phase": phase,
            "operation": operation,
            "component": component,
            "task_identity": task_identity,
            "deployment_identity": deployment_identity,
            "authority_platform_identity": authority_platform_identity,
            "prior_state_identity": prior_state_identity,
            "candidate_state_identity": candidate_state_identity,
            "candidate_active": candidate_active,
        }
        value["identity"] = canonical_identity(
            value, scheme=ACTION_REQUEST_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "DeploymentActionRequest":
        code = "OPERATOR_ACTION_REQUEST_INVALID"
        if not isinstance(raw, dict) or set(raw) != {
            "schema_version",
            "identity",
            "phase",
            "operation",
            "component",
            "task_identity",
            "deployment_identity",
            "authority_platform_identity",
            "prior_state_identity",
            "candidate_state_identity",
            "candidate_active",
        }:
            raise fail(code, "Operator action request is invalid.")
        active = raw["candidate_active"]
        if (
            raw["schema_version"] != ACTION_REQUEST_SCHEMA
            or raw["phase"] not in ACTION_PHASES
            or raw["operation"] not in ACTION_OPERATIONS
            or raw["component"] not in COMPONENTS
            or (raw["phase"] == "observe") != (raw["operation"] == "observe")
            or not isinstance(active, dict)
            or set(active) != set(COMPONENTS)
            or any(
                value is not None
                and (not isinstance(value, str) or _IDENTITY.fullmatch(value) is None)
                for value in active.values()
            )
        ):
            raise fail(code, "Operator action request is invalid.")
        identity = _identity(raw["identity"], code=code)
        if identity != canonical_identity(
            without_key(raw, "identity"), scheme=ACTION_REQUEST_IDENTITY_SCHEME
        ):
            raise fail(code, "Operator action request is invalid.")
        return cls(
            identity=identity,
            phase=raw["phase"],
            operation=raw["operation"],
            component=raw["component"],
            task_identity=_task(raw["task_identity"], code=code),
            deployment_identity=_identity(raw["deployment_identity"], code=code),
            authority_platform_identity=_identity(
                raw["authority_platform_identity"], code=code
            ),
            prior_state_identity=_identity(raw["prior_state_identity"], code=code),
            candidate_state_identity=_identity(
                raw["candidate_state_identity"], code=code
            ),
            candidate_active={component: active[component] for component in COMPONENTS},
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": ACTION_REQUEST_SCHEMA,
            "identity": self.identity,
            "phase": self.phase,
            "operation": self.operation,
            "component": self.component,
            "task_identity": self.task_identity,
            "deployment_identity": self.deployment_identity,
            "authority_platform_identity": self.authority_platform_identity,
            "prior_state_identity": self.prior_state_identity,
            "candidate_state_identity": self.candidate_state_identity,
            "candidate_active": dict(self.candidate_active),
        }


@dataclass(frozen=True)
class ReadinessCheck:
    status: str
    reason_code: str
    identity: str | None = None

    @classmethod
    def from_dict(cls, raw: object) -> "ReadinessCheck":
        code = "OPERATOR_ACTION_RECEIPT_INVALID"
        if (
            not isinstance(raw, dict)
            or set(raw) != {"status", "reason_code", "identity"}
            or raw["status"] not in READINESS_STATUSES
            or raw["reason_code"] not in READINESS_REASONS
        ):
            raise fail(code, "Operator action receipt is invalid.")
        identity = _optional_identity(raw["identity"], code=code)
        if (raw["status"] == "ready") != (identity is not None) or raw[
            "reason_code"
        ] not in _READINESS_REASON_BY_STATUS[raw["status"]]:
            raise fail(code, "Operator action receipt is invalid.")
        return cls(raw["status"], raw["reason_code"], identity)

    def to_dict(self) -> dict[str, object]:
        return {
            "status": self.status,
            "reason_code": self.reason_code,
            "identity": self.identity,
        }


@dataclass(frozen=True)
class DeploymentActionReceipt:
    identity: str
    request_identity: str
    status: str
    compatibility: str
    database_before_identity: str | None
    database_after_identity: str | None
    accepted_schema_heads: tuple[str, ...]
    target_schema_heads: tuple[str, ...]
    migration_inventory_identity: str
    known_schema_revisions: tuple[str, ...]
    migration_required: bool
    rollback_supported: bool
    api_contract_sha256: str
    native_identities: dict[str, str | None]
    frontend_identity: str | None
    reference_compatibility_identity: str | None
    readiness: dict[str, ReadinessCheck]

    @classmethod
    def create(
        cls,
        *,
        request_identity: str,
        status: str,
        compatibility: str,
        database_before_identity: str | None = None,
        database_after_identity: str | None = None,
        accepted_schema_heads: Sequence[str] = (),
        target_schema_heads: Sequence[str] = (),
        migration_inventory_identity: str,
        known_schema_revisions: Sequence[str],
        migration_required: bool = False,
        rollback_supported: bool = False,
        api_contract_sha256: str,
        native_identities: dict[str, str | None],
        frontend_identity: str | None,
        reference_compatibility_identity: str | None,
        readiness: dict[str, ReadinessCheck],
    ) -> "DeploymentActionReceipt":
        value: dict[str, object] = {
            "schema_version": ACTION_RECEIPT_SCHEMA,
            "request_identity": request_identity,
            "status": status,
            "compatibility": compatibility,
            "database_before_identity": database_before_identity,
            "database_after_identity": database_after_identity,
            "accepted_schema_heads": list(accepted_schema_heads),
            "target_schema_heads": list(target_schema_heads),
            "migration_inventory_identity": migration_inventory_identity,
            "known_schema_revisions": list(known_schema_revisions),
            "migration_required": migration_required,
            "rollback_supported": rollback_supported,
            "api_contract_sha256": api_contract_sha256,
            "native_identities": native_identities,
            "frontend_identity": frontend_identity,
            "reference_compatibility_identity": reference_compatibility_identity,
            "readiness": {
                name: readiness[name].to_dict() for name in VERIFICATION_CHECKS
            },
        }
        value["identity"] = canonical_identity(
            value, scheme=ACTION_RECEIPT_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "DeploymentActionReceipt":
        code = "OPERATOR_ACTION_RECEIPT_INVALID"
        if not isinstance(raw, dict) or set(raw) != {
            "schema_version",
            "identity",
            "request_identity",
            "status",
            "compatibility",
            "database_before_identity",
            "database_after_identity",
            "accepted_schema_heads",
            "target_schema_heads",
            "migration_inventory_identity",
            "known_schema_revisions",
            "migration_required",
            "rollback_supported",
            "api_contract_sha256",
            "native_identities",
            "frontend_identity",
            "reference_compatibility_identity",
            "readiness",
        }:
            raise fail(code, "Operator action receipt is invalid.")
        readiness = raw["readiness"]
        native = raw["native_identities"]
        if (
            raw["schema_version"] != ACTION_RECEIPT_SCHEMA
            or raw["status"] not in ACTION_STATUSES
            or raw["compatibility"] not in COMPATIBILITY_STATUSES
            or not isinstance(raw["migration_required"], bool)
            or not isinstance(raw["rollback_supported"], bool)
            or not isinstance(readiness, dict)
            or set(readiness) != set(VERIFICATION_CHECKS)
            or not isinstance(native, dict)
            or set(native) != set(COMPONENTS)
            or any(
                value is not None
                and (not isinstance(value, str) or _IDENTITY.fullmatch(value) is None)
                for value in native.values()
            )
            or not isinstance(raw["api_contract_sha256"], str)
            or re.fullmatch(r"[0-9a-f]{64}", raw["api_contract_sha256"]) is None
        ):
            raise fail(code, "Operator action receipt is invalid.")
        identity = _identity(raw["identity"], code=code)
        if identity != canonical_identity(
            without_key(raw, "identity"), scheme=ACTION_RECEIPT_IDENTITY_SCHEME
        ):
            raise fail(code, "Operator action receipt is invalid.")
        accepted = _heads(raw["accepted_schema_heads"], code=code)
        target = _heads(raw["target_schema_heads"], code=code)
        known_revisions = _heads(raw["known_schema_revisions"], code=code)
        if (
            len(target) != 1
            or target[0] in known_revisions
            or any(head not in {*known_revisions, *target} for head in accepted)
            or raw["migration_required"] != (accepted != target)
            or raw["rollback_supported"] != (accepted == target)
        ):
            raise fail(code, "Operator action receipt is invalid.")
        return cls(
            identity=identity,
            request_identity=_identity(raw["request_identity"], code=code),
            status=raw["status"],
            compatibility=raw["compatibility"],
            database_before_identity=_optional_identity(
                raw["database_before_identity"], code=code
            ),
            database_after_identity=_optional_identity(
                raw["database_after_identity"], code=code
            ),
            accepted_schema_heads=accepted,
            target_schema_heads=target,
            migration_inventory_identity=_identity(
                raw["migration_inventory_identity"], code=code
            ),
            known_schema_revisions=known_revisions,
            migration_required=raw["migration_required"],
            rollback_supported=raw["rollback_supported"],
            api_contract_sha256=raw["api_contract_sha256"],
            native_identities={
                component: native[component] for component in COMPONENTS
            },
            frontend_identity=_optional_identity(raw["frontend_identity"], code=code),
            reference_compatibility_identity=_optional_identity(
                raw["reference_compatibility_identity"], code=code
            ),
            readiness={
                name: ReadinessCheck.from_dict(readiness[name])
                for name in VERIFICATION_CHECKS
            },
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": ACTION_RECEIPT_SCHEMA,
            "identity": self.identity,
            "request_identity": self.request_identity,
            "status": self.status,
            "compatibility": self.compatibility,
            "database_before_identity": self.database_before_identity,
            "database_after_identity": self.database_after_identity,
            "accepted_schema_heads": list(self.accepted_schema_heads),
            "target_schema_heads": list(self.target_schema_heads),
            "migration_inventory_identity": self.migration_inventory_identity,
            "known_schema_revisions": list(self.known_schema_revisions),
            "migration_required": self.migration_required,
            "rollback_supported": self.rollback_supported,
            "api_contract_sha256": self.api_contract_sha256,
            "native_identities": dict(self.native_identities),
            "frontend_identity": self.frontend_identity,
            "reference_compatibility_identity": self.reference_compatibility_identity,
            "readiness": {
                name: self.readiness[name].to_dict() for name in VERIFICATION_CHECKS
            },
        }


class DeploymentActionRunner(Protocol):
    def run(self, request: DeploymentActionRequest) -> DeploymentActionReceipt: ...


@dataclass(frozen=True)
class DatabasePrepareRequest:
    identity: str
    operation: str
    database_mode: str
    task_identity: str
    deployment_identity: str
    prior_state_identity: str
    candidate_state_identity: str
    action_receipt_identity: str
    backup_receipt_identity: str | None
    target_schema_heads: tuple[str, ...]

    @classmethod
    def create(
        cls,
        *,
        operation: str,
        database_mode: str,
        task_identity: str,
        deployment_identity: str,
        prior_state_identity: str,
        candidate_state_identity: str,
        action_receipt_identity: str,
        backup_receipt_identity: str | None,
        target_schema_heads: Sequence[str],
    ) -> "DatabasePrepareRequest":
        value: dict[str, object] = {
            "schema_version": DATABASE_PREPARE_REQUEST_SCHEMA,
            "operation": operation,
            "component": "platform",
            "database_mode": database_mode,
            "task_identity": task_identity,
            "deployment_identity": deployment_identity,
            "prior_state_identity": prior_state_identity,
            "candidate_state_identity": candidate_state_identity,
            "action_receipt_identity": action_receipt_identity,
            "backup_receipt_identity": backup_receipt_identity,
            "target_schema_heads": list(target_schema_heads),
        }
        value["identity"] = canonical_identity(
            value, scheme=DATABASE_PREPARE_REQUEST_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "DatabasePrepareRequest":
        code = "DATABASE_PREPARE_REQUEST_INVALID"
        if not isinstance(raw, dict) or set(raw) != {
            "schema_version",
            "identity",
            "operation",
            "component",
            "database_mode",
            "task_identity",
            "deployment_identity",
            "prior_state_identity",
            "candidate_state_identity",
            "action_receipt_identity",
            "backup_receipt_identity",
            "target_schema_heads",
        }:
            raise fail(code, "Database preparation request is invalid.")
        if (
            raw["schema_version"] != DATABASE_PREPARE_REQUEST_SCHEMA
            or raw["operation"] not in {"activate", "rollback"}
            or raw["component"] != "platform"
            or raw["database_mode"] not in DATABASE_MODES
            or (
                raw["database_mode"] == "fresh-candidate"
                and (
                    raw["operation"] != "activate"
                    or raw["backup_receipt_identity"] is not None
                )
            )
            or (
                raw["database_mode"] == "existing-live"
                and raw["backup_receipt_identity"] is None
            )
        ):
            raise fail(code, "Database preparation request is invalid.")
        identity = _identity(raw["identity"], code=code)
        if identity != canonical_identity(
            without_key(raw, "identity"),
            scheme=DATABASE_PREPARE_REQUEST_IDENTITY_SCHEME,
        ):
            raise fail(code, "Database preparation request is invalid.")
        target_schema_heads = _heads(
            raw["target_schema_heads"], code=code, allow_empty=False
        )
        if len(target_schema_heads) != 1:
            raise fail(code, "Database preparation request is invalid.")
        return cls(
            identity=identity,
            operation=raw["operation"],
            database_mode=raw["database_mode"],
            task_identity=_task(raw["task_identity"], code=code),
            deployment_identity=_identity(raw["deployment_identity"], code=code),
            prior_state_identity=_identity(raw["prior_state_identity"], code=code),
            candidate_state_identity=_identity(
                raw["candidate_state_identity"], code=code
            ),
            action_receipt_identity=_identity(
                raw["action_receipt_identity"], code=code
            ),
            backup_receipt_identity=_optional_identity(
                raw["backup_receipt_identity"], code=code
            ),
            target_schema_heads=target_schema_heads,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": DATABASE_PREPARE_REQUEST_SCHEMA,
            "identity": self.identity,
            "operation": self.operation,
            "component": "platform",
            "database_mode": self.database_mode,
            "task_identity": self.task_identity,
            "deployment_identity": self.deployment_identity,
            "prior_state_identity": self.prior_state_identity,
            "candidate_state_identity": self.candidate_state_identity,
            "action_receipt_identity": self.action_receipt_identity,
            "backup_receipt_identity": self.backup_receipt_identity,
            "target_schema_heads": list(self.target_schema_heads),
        }


@dataclass(frozen=True)
class DatabasePrepareReceipt:
    identity: str
    request_identity: str
    database_before_identity: str | None
    database_after_identity: str
    schema_heads: tuple[str, ...]

    @classmethod
    def create(
        cls,
        *,
        request_identity: str,
        database_before_identity: str | None,
        database_after_identity: str,
        schema_heads: Sequence[str],
    ) -> "DatabasePrepareReceipt":
        value: dict[str, object] = {
            "schema_version": DATABASE_PREPARE_RECEIPT_SCHEMA,
            "status": "prepared",
            "request_identity": request_identity,
            "database_before_identity": database_before_identity,
            "database_after_identity": database_after_identity,
            "schema_heads": list(schema_heads),
        }
        value["identity"] = canonical_identity(
            value, scheme=DATABASE_PREPARE_RECEIPT_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "DatabasePrepareReceipt":
        code = "DATABASE_PREPARE_RECEIPT_INVALID"
        if (
            not isinstance(raw, dict)
            or set(raw)
            != {
                "schema_version",
                "identity",
                "status",
                "request_identity",
                "database_before_identity",
                "database_after_identity",
                "schema_heads",
            }
            or (
                raw["schema_version"] != DATABASE_PREPARE_RECEIPT_SCHEMA
                or raw["status"] != "prepared"
            )
        ):
            raise fail(code, "Database preparation receipt is invalid.")
        identity = _identity(raw["identity"], code=code)
        if identity != canonical_identity(
            without_key(raw, "identity"),
            scheme=DATABASE_PREPARE_RECEIPT_IDENTITY_SCHEME,
        ):
            raise fail(code, "Database preparation receipt is invalid.")
        return cls(
            identity=identity,
            request_identity=_identity(raw["request_identity"], code=code),
            database_before_identity=_optional_identity(
                raw["database_before_identity"], code=code
            ),
            database_after_identity=_identity(
                raw["database_after_identity"], code=code
            ),
            schema_heads=_heads(raw["schema_heads"], code=code, allow_empty=False),
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": DATABASE_PREPARE_RECEIPT_SCHEMA,
            "identity": self.identity,
            "status": "prepared",
            "request_identity": self.request_identity,
            "database_before_identity": self.database_before_identity,
            "database_after_identity": self.database_after_identity,
            "schema_heads": list(self.schema_heads),
        }


class DatabasePreparer(Protocol):
    def prepare(self, request: DatabasePrepareRequest) -> DatabasePrepareReceipt: ...


def snakemake_environment_hash(
    conda_prefix: Path,
    environment_content: bytes,
) -> str:
    """Return Snakemake 8's final-prefix-bound full environment hash."""

    if (
        not isinstance(conda_prefix, Path)
        or not conda_prefix.is_absolute()
        or not isinstance(environment_content, bytes)
        or not environment_content
    ):
        raise fail(
            "ENCODE_RUNTIME_ENVIRONMENT_INVALID",
            "ENCODE runtime environment is invalid.",
            component=ENCODE_RUNTIME,
        )
    digest = hashlib.md5(usedforsecurity=False)
    digest.update(os.path.realpath(conda_prefix).encode("utf-8"))
    digest.update(environment_content)
    return digest.hexdigest()


def _safe_runtime_member(value: object) -> bool:
    if not isinstance(value, str) or value.startswith("/"):
        return False
    try:
        encoded = value.encode("utf-8", errors="strict")
        segments = value.split("/")
        encoded_segments = tuple(
            segment.encode("utf-8", errors="strict") for segment in segments
        )
    except UnicodeError:
        return False
    return (
        0 < len(encoded) <= 1024
        and all(character not in value for character in ("\x00", "\r", "\n"))
        and all(segment not in {"", ".", ".."} for segment in segments)
        and all(0 < len(segment) <= 255 for segment in encoded_segments)
    )


def _safe_runtime_symlink(path: str, target: str) -> bool:
    if (
        not target
        or target.startswith("/")
        or any(character in target for character in ("\x00", "\n", "\r"))
    ):
        return False
    parts = path.split("/")[:-1]
    for part in target.split("/"):
        if part in {"", "."}:
            continue
        if part == "..":
            if not parts:
                return False
            parts.pop()
            continue
        try:
            encoded = part.encode("utf-8", errors="strict")
        except UnicodeError:
            return False
        if not 0 < len(encoded) <= 255:
            return False
        parts.append(part)
    resolved = "/".join(parts)
    return _safe_runtime_member(resolved)


@dataclass(frozen=True, order=True)
class EncodeRuntimeEntry:
    path: str
    kind: str
    sha256: str
    size: int
    mode: int | None
    target: str | None

    @classmethod
    def from_dict(cls, raw: object) -> "EncodeRuntimeEntry":
        code = "ENCODE_RUNTIME_PREPARE_RECEIPT_INVALID"
        if (
            not isinstance(raw, dict)
            or set(raw) != {"path", "kind", "sha256", "size", "mode", "target"}
            or not isinstance(raw["path"], str)
            or not _safe_runtime_member(raw["path"])
            or raw["path"] == ".helixweave-runtime-inventory.json"
            or raw["kind"] not in {"file", "symlink"}
            or not isinstance(raw["sha256"], str)
            or re.fullmatch(r"[0-9a-f]{64}", raw["sha256"]) is None
            or not isinstance(raw["size"], int)
            or isinstance(raw["size"], bool)
            or not 0 <= raw["size"] <= _MAX_RUNTIME_FILE_BYTES
        ):
            raise fail(code, "ENCODE runtime preparation receipt is invalid.")
        if raw["kind"] == "file":
            if raw["mode"] not in {0o444, 0o555} or raw["target"] is not None:
                raise fail(code, "ENCODE runtime preparation receipt is invalid.")
        elif (
            raw["mode"] is not None
            or not isinstance(raw["target"], str)
            or not _safe_runtime_symlink(raw["path"], raw["target"])
            or raw["size"] != len(raw["target"].encode("utf-8"))
            or raw["sha256"]
            != hashlib.sha256(raw["target"].encode("utf-8")).hexdigest()
        ):
            raise fail(code, "ENCODE runtime preparation receipt is invalid.")
        return cls(
            raw["path"],
            raw["kind"],
            raw["sha256"],
            raw["size"],
            raw["mode"],
            raw["target"],
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "path": self.path,
            "kind": self.kind,
            "sha256": self.sha256,
            "size": self.size,
            "mode": self.mode,
            "target": self.target,
        }


@dataclass(frozen=True)
class EncodeRuntimeInventory:
    tree_identity: str
    entries: tuple[EncodeRuntimeEntry, ...]

    @classmethod
    def create(cls, entries: Sequence[EncodeRuntimeEntry]) -> "EncodeRuntimeInventory":
        ordered = tuple(entries)
        return cls.from_dict(
            {
                "schema_version": ENCODE_RUNTIME_INVENTORY_SCHEMA,
                "tree_identity": canonical_identity(
                    {"entries": [item.to_dict() for item in ordered]},
                    scheme=ENCODE_RUNTIME_TREE_IDENTITY_SCHEME,
                ),
                "entries": [item.to_dict() for item in ordered],
            }
        )

    @classmethod
    def from_dict(cls, raw: object) -> "EncodeRuntimeInventory":
        code = "ENCODE_RUNTIME_PREPARE_RECEIPT_INVALID"
        if not isinstance(raw, dict) or set(raw) != {
            "schema_version",
            "tree_identity",
            "entries",
        }:
            raise fail(code, "ENCODE runtime preparation receipt is invalid.")
        raw_entries = raw["entries"]
        if (
            raw["schema_version"] != ENCODE_RUNTIME_INVENTORY_SCHEMA
            or not isinstance(raw_entries, list)
            or not 0 < len(raw_entries) <= _MAX_RUNTIME_FILES
        ):
            raise fail(code, "ENCODE runtime preparation receipt is invalid.")
        entries = tuple(EncodeRuntimeEntry.from_dict(item) for item in raw_entries)
        if (
            tuple(sorted(entries, key=lambda item: item.path)) != entries
            or len({item.path for item in entries}) != len(entries)
            or sum(item.size for item in entries) > _MAX_RUNTIME_TOTAL_BYTES
        ):
            raise fail(code, "ENCODE runtime preparation receipt is invalid.")
        by_path = {item.path: item for item in entries}
        snakemake = by_path.get("runner/bin/snakemake")
        conda = by_path.get("runner/bin/conda")
        if (
            snakemake is None
            or snakemake.kind != "file"
            or snakemake.mode != 0o555
            or conda is None
            or conda.kind != "file"
            or conda.mode != 0o555
            or not any(item.path.startswith("conda-envs/") for item in entries)
        ):
            raise fail(code, "ENCODE runtime preparation receipt is invalid.")
        tree_identity = _identity(raw["tree_identity"], code=code)
        if tree_identity != canonical_identity(
            {"entries": [item.to_dict() for item in entries]},
            scheme=ENCODE_RUNTIME_TREE_IDENTITY_SCHEME,
        ):
            raise fail(code, "ENCODE runtime preparation receipt is invalid.")
        return cls(tree_identity, entries)

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": ENCODE_RUNTIME_INVENTORY_SCHEMA,
            "tree_identity": self.tree_identity,
            "entries": [item.to_dict() for item in self.entries],
        }

    def to_bytes(self) -> bytes:
        return canonical_json_bytes(self.to_dict())


@dataclass(frozen=True)
class EncodeRuntimePrepareRequest:
    identity: str
    task_identity: str
    deployment_identity: str
    authority_platform_identity: str
    prior_state_identity: str
    candidate_state_identity: str

    @classmethod
    def create(
        cls,
        *,
        task_identity: str,
        deployment_identity: str,
        authority_platform_identity: str,
        prior_state_identity: str,
        candidate_state_identity: str,
    ) -> "EncodeRuntimePrepareRequest":
        value: dict[str, object] = {
            "schema_version": ENCODE_RUNTIME_PREPARE_REQUEST_SCHEMA,
            "operation": "materialize",
            "component": "encode-runtime",
            "task_identity": task_identity,
            "deployment_identity": deployment_identity,
            "authority_platform_identity": authority_platform_identity,
            "prior_state_identity": prior_state_identity,
            "candidate_state_identity": candidate_state_identity,
        }
        value["identity"] = canonical_identity(
            value, scheme=ENCODE_RUNTIME_PREPARE_REQUEST_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "EncodeRuntimePrepareRequest":
        code = "ENCODE_RUNTIME_PREPARE_REQUEST_INVALID"
        if not isinstance(raw, dict) or set(raw) != {
            "schema_version",
            "identity",
            "operation",
            "component",
            "task_identity",
            "deployment_identity",
            "authority_platform_identity",
            "prior_state_identity",
            "candidate_state_identity",
        }:
            raise fail(code, "ENCODE runtime preparation request is invalid.")
        identity = _identity(raw["identity"], code=code)
        if (
            raw["schema_version"] != ENCODE_RUNTIME_PREPARE_REQUEST_SCHEMA
            or raw["operation"] != "materialize"
            or raw["component"] != "encode-runtime"
            or identity
            != canonical_identity(
                without_key(raw, "identity"),
                scheme=ENCODE_RUNTIME_PREPARE_REQUEST_IDENTITY_SCHEME,
            )
        ):
            raise fail(code, "ENCODE runtime preparation request is invalid.")
        return cls(
            identity=identity,
            task_identity=_task(raw["task_identity"], code=code),
            deployment_identity=_identity(raw["deployment_identity"], code=code),
            authority_platform_identity=_identity(
                raw["authority_platform_identity"], code=code
            ),
            prior_state_identity=_identity(raw["prior_state_identity"], code=code),
            candidate_state_identity=_identity(
                raw["candidate_state_identity"], code=code
            ),
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": ENCODE_RUNTIME_PREPARE_REQUEST_SCHEMA,
            "identity": self.identity,
            "operation": "materialize",
            "component": "encode-runtime",
            "task_identity": self.task_identity,
            "deployment_identity": self.deployment_identity,
            "authority_platform_identity": self.authority_platform_identity,
            "prior_state_identity": self.prior_state_identity,
            "candidate_state_identity": self.candidate_state_identity,
        }


@dataclass(frozen=True)
class EncodeRuntimePrepareReceipt:
    identity: str
    request_identity: str
    deployment_identity: str
    tree_identity: str
    inventory_sha256: str
    inventory_size: int
    entry_count: int

    @classmethod
    def create(
        cls,
        *,
        request_identity: str,
        deployment_identity: str,
        inventory: EncodeRuntimeInventory,
    ) -> "EncodeRuntimePrepareReceipt":
        if not isinstance(inventory, EncodeRuntimeInventory):
            raise fail(
                "ENCODE_RUNTIME_PREPARE_RECEIPT_INVALID",
                "ENCODE runtime preparation receipt is invalid.",
            )
        inventory_content = inventory.to_bytes()
        value: dict[str, object] = {
            "schema_version": ENCODE_RUNTIME_PREPARE_RECEIPT_SCHEMA,
            "status": "materialized",
            "request_identity": request_identity,
            "deployment_identity": deployment_identity,
            "tree_identity": inventory.tree_identity,
            "inventory_sha256": hashlib.sha256(inventory_content).hexdigest(),
            "inventory_size": len(inventory_content),
            "entry_count": len(inventory.entries),
        }
        value["identity"] = canonical_identity(
            value, scheme=ENCODE_RUNTIME_PREPARE_RECEIPT_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "EncodeRuntimePrepareReceipt":
        code = "ENCODE_RUNTIME_PREPARE_RECEIPT_INVALID"
        if not isinstance(raw, dict) or set(raw) != {
            "schema_version",
            "identity",
            "status",
            "request_identity",
            "deployment_identity",
            "tree_identity",
            "inventory_sha256",
            "inventory_size",
            "entry_count",
        }:
            raise fail(code, "ENCODE runtime preparation receipt is invalid.")
        if (
            not isinstance(raw["inventory_sha256"], str)
            or re.fullmatch(r"[0-9a-f]{64}", raw["inventory_sha256"]) is None
            or not isinstance(raw["inventory_size"], int)
            or isinstance(raw["inventory_size"], bool)
            or not 0 < raw["inventory_size"] <= _MAX_RUNTIME_INVENTORY_BYTES
            or not isinstance(raw["entry_count"], int)
            or isinstance(raw["entry_count"], bool)
            or not 0 < raw["entry_count"] <= _MAX_RUNTIME_FILES
        ):
            raise fail(code, "ENCODE runtime preparation receipt is invalid.")
        tree_identity = _identity(raw["tree_identity"], code=code)
        identity = _identity(raw["identity"], code=code)
        if (
            raw["schema_version"] != ENCODE_RUNTIME_PREPARE_RECEIPT_SCHEMA
            or raw["status"] != "materialized"
            or identity
            != canonical_identity(
                without_key(raw, "identity"),
                scheme=ENCODE_RUNTIME_PREPARE_RECEIPT_IDENTITY_SCHEME,
            )
        ):
            raise fail(code, "ENCODE runtime preparation receipt is invalid.")
        return cls(
            identity=identity,
            request_identity=_identity(raw["request_identity"], code=code),
            deployment_identity=_identity(raw["deployment_identity"], code=code),
            tree_identity=tree_identity,
            inventory_sha256=raw["inventory_sha256"],
            inventory_size=raw["inventory_size"],
            entry_count=raw["entry_count"],
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": ENCODE_RUNTIME_PREPARE_RECEIPT_SCHEMA,
            "identity": self.identity,
            "status": "materialized",
            "request_identity": self.request_identity,
            "deployment_identity": self.deployment_identity,
            "tree_identity": self.tree_identity,
            "inventory_sha256": self.inventory_sha256,
            "inventory_size": self.inventory_size,
            "entry_count": self.entry_count,
        }


class EncodeRuntimePreparer(Protocol):
    def prepare(
        self, request: EncodeRuntimePrepareRequest
    ) -> EncodeRuntimePrepareReceipt: ...


@dataclass(frozen=True)
class BulkRuntimePrepareRequest:
    """Root-selected request to make one admitted image set locally available."""

    identity: str
    operation: str
    task_identity: str
    candidate_bulk_identity: str
    authority_platform_identity: str
    prior_state_identity: str
    candidate_state_identity: str
    docker_service_identity: str
    docker_client_identity: str
    docker_endpoint_identity: str
    docker_daemon_uid: int
    docker_daemon_gid: int

    @classmethod
    def create(
        cls,
        *,
        operation: str,
        task_identity: str,
        candidate_bulk_identity: str,
        authority_platform_identity: str,
        prior_state_identity: str,
        candidate_state_identity: str,
        docker_service_identity: str,
        docker_client_identity: str,
        docker_endpoint_identity: str,
        docker_daemon_uid: int,
        docker_daemon_gid: int,
    ) -> "BulkRuntimePrepareRequest":
        value: dict[str, object] = {
            "schema_version": BULK_RUNTIME_PREPARE_REQUEST_SCHEMA,
            "operation": operation,
            "component": BULK_RNASEQ_RUNTIME,
            "task_identity": task_identity,
            "candidate_bulk_identity": candidate_bulk_identity,
            "authority_platform_identity": authority_platform_identity,
            "prior_state_identity": prior_state_identity,
            "candidate_state_identity": candidate_state_identity,
            "docker_service_identity": docker_service_identity,
            "docker_client_identity": docker_client_identity,
            "docker_endpoint_identity": docker_endpoint_identity,
            "docker_daemon_uid": docker_daemon_uid,
            "docker_daemon_gid": docker_daemon_gid,
        }
        value["identity"] = canonical_identity(
            value, scheme=BULK_RUNTIME_PREPARE_REQUEST_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "BulkRuntimePrepareRequest":
        code = "BULK_RUNTIME_PREPARE_REQUEST_INVALID"
        if not isinstance(raw, dict) or set(raw) != {
            "schema_version",
            "identity",
            "operation",
            "component",
            "task_identity",
            "candidate_bulk_identity",
            "authority_platform_identity",
            "prior_state_identity",
            "candidate_state_identity",
            "docker_service_identity",
            "docker_client_identity",
            "docker_endpoint_identity",
            "docker_daemon_uid",
            "docker_daemon_gid",
        }:
            raise fail(code, "Bulk runtime preparation request is invalid.")
        identity = _identity(raw["identity"], code=code)
        if (
            raw["schema_version"] != BULK_RUNTIME_PREPARE_REQUEST_SCHEMA
            or raw["operation"] not in BULK_RUNTIME_PREPARE_OPERATIONS
            or raw["component"] != BULK_RNASEQ_RUNTIME
            or any(
                not isinstance(value, int)
                or isinstance(value, bool)
                or not 0 <= value <= 2**32 - 1
                for value in (
                    raw["docker_daemon_uid"],
                    raw["docker_daemon_gid"],
                )
            )
            or identity
            != canonical_identity(
                without_key(raw, "identity"),
                scheme=BULK_RUNTIME_PREPARE_REQUEST_IDENTITY_SCHEME,
            )
        ):
            raise fail(code, "Bulk runtime preparation request is invalid.")
        return cls(
            identity=identity,
            operation=raw["operation"],
            task_identity=_task(raw["task_identity"], code=code),
            candidate_bulk_identity=_identity(
                raw["candidate_bulk_identity"], code=code
            ),
            authority_platform_identity=_identity(
                raw["authority_platform_identity"], code=code
            ),
            prior_state_identity=_identity(raw["prior_state_identity"], code=code),
            candidate_state_identity=_identity(
                raw["candidate_state_identity"], code=code
            ),
            docker_service_identity=_identity(
                raw["docker_service_identity"], code=code
            ),
            docker_client_identity=_identity(raw["docker_client_identity"], code=code),
            docker_endpoint_identity=_identity(
                raw["docker_endpoint_identity"], code=code
            ),
            docker_daemon_uid=raw["docker_daemon_uid"],
            docker_daemon_gid=raw["docker_daemon_gid"],
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": BULK_RUNTIME_PREPARE_REQUEST_SCHEMA,
            "identity": self.identity,
            "operation": self.operation,
            "component": BULK_RNASEQ_RUNTIME,
            "task_identity": self.task_identity,
            "candidate_bulk_identity": self.candidate_bulk_identity,
            "authority_platform_identity": self.authority_platform_identity,
            "prior_state_identity": self.prior_state_identity,
            "candidate_state_identity": self.candidate_state_identity,
            "docker_service_identity": self.docker_service_identity,
            "docker_client_identity": self.docker_client_identity,
            "docker_endpoint_identity": self.docker_endpoint_identity,
            "docker_daemon_uid": self.docker_daemon_uid,
            "docker_daemon_gid": self.docker_daemon_gid,
        }


@dataclass(frozen=True)
class BulkRuntimePrepareReceipt:
    """Path-free attestation for one fully inspected local image set."""

    identity: str
    request_identity: str
    candidate_bulk_identity: str
    runtime_identity: str
    image_set_identity: str
    image_count: int

    @classmethod
    def create(
        cls,
        *,
        request_identity: str,
        candidate_bulk_identity: str,
        runtime_identity: str,
        image_set_identity: str,
        image_count: int,
    ) -> "BulkRuntimePrepareReceipt":
        value: dict[str, object] = {
            "schema_version": BULK_RUNTIME_PREPARE_RECEIPT_SCHEMA,
            "status": "prepared",
            "request_identity": request_identity,
            "candidate_bulk_identity": candidate_bulk_identity,
            "runtime_identity": runtime_identity,
            "image_set_identity": image_set_identity,
            "image_count": image_count,
        }
        value["identity"] = canonical_identity(
            value, scheme=BULK_RUNTIME_PREPARE_RECEIPT_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "BulkRuntimePrepareReceipt":
        code = "BULK_RUNTIME_PREPARE_RECEIPT_INVALID"
        if not isinstance(raw, dict) or set(raw) != {
            "schema_version",
            "identity",
            "status",
            "request_identity",
            "candidate_bulk_identity",
            "runtime_identity",
            "image_set_identity",
            "image_count",
        }:
            raise fail(code, "Bulk runtime preparation receipt is invalid.")
        identity = _identity(raw["identity"], code=code)
        if (
            raw["schema_version"] != BULK_RUNTIME_PREPARE_RECEIPT_SCHEMA
            or raw["status"] != "prepared"
            or not isinstance(raw["image_count"], int)
            or isinstance(raw["image_count"], bool)
            or not 0 < raw["image_count"] <= _MAX_BULK_RUNTIME_IMAGES
            or identity
            != canonical_identity(
                without_key(raw, "identity"),
                scheme=BULK_RUNTIME_PREPARE_RECEIPT_IDENTITY_SCHEME,
            )
        ):
            raise fail(code, "Bulk runtime preparation receipt is invalid.")
        return cls(
            identity=identity,
            request_identity=_identity(raw["request_identity"], code=code),
            candidate_bulk_identity=_identity(
                raw["candidate_bulk_identity"], code=code
            ),
            runtime_identity=_identity(raw["runtime_identity"], code=code),
            image_set_identity=_identity(raw["image_set_identity"], code=code),
            image_count=raw["image_count"],
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": BULK_RUNTIME_PREPARE_RECEIPT_SCHEMA,
            "identity": self.identity,
            "status": "prepared",
            "request_identity": self.request_identity,
            "candidate_bulk_identity": self.candidate_bulk_identity,
            "runtime_identity": self.runtime_identity,
            "image_set_identity": self.image_set_identity,
            "image_count": self.image_count,
        }


ActionCommandRunner = Callable[[tuple[str, ...]], int]
UnitQuiescer = Callable[[], bool]


class SystemdDeploymentActionRunner:
    """Exchange canonical documents around one fixed systemd oneshot."""

    def __init__(
        self,
        layout: DeploymentLayout,
        *,
        service_uid: int,
        service_gid: int,
        root_uid: int = 0,
        root_gid: int = 0,
        command_runner: ActionCommandRunner | None = None,
    ) -> None:
        if any(
            not isinstance(value, int) or isinstance(value, bool) or value < 0
            for value in (service_uid, service_gid, root_uid, root_gid)
        ):
            raise fail(
                "OPERATOR_ACTION_BOUNDARY_INVALID", "Operator action is unavailable."
            )
        self.layout = layout
        self.service_uid = service_uid
        self.service_gid = service_gid
        self.root_uid = root_uid
        self.root_gid = root_gid
        self.command_runner = (
            self._run_systemctl if command_runner is None else command_runner
        )

    def run(self, request: DeploymentActionRequest) -> DeploymentActionReceipt:
        if not isinstance(request, DeploymentActionRequest):
            raise fail(
                "OPERATOR_ACTION_REQUEST_INVALID", "Operator action request is invalid."
            )
        root = self._require_root()
        self._replace_request(root, request)
        self._prepare_receipt(root)
        try:
            return_code = self.command_runner(
                (
                    str(SYSTEMCTL),
                    "--no-ask-password",
                    "start",
                    "--",
                    ACTION_UNIT,
                )
            )
        except Exception:
            raise fail(
                "OPERATOR_ACTION_FAILED",
                "Operator action failed.",
                recoverable=True,
            ) from None
        if return_code != 0:
            raise fail(
                "OPERATOR_ACTION_FAILED",
                "Operator action failed.",
                recoverable=True,
            )
        receipt = self._read_receipt(request)
        if receipt.request_identity != request.identity:
            raise fail(
                "OPERATOR_ACTION_RECEIPT_INVALID", "Operator action receipt is invalid."
            )
        return receipt

    def _require_root(self) -> Path:
        path = self.layout.operator_action_root
        try:
            observed = path.lstat()
        except OSError:
            raise fail(
                "OPERATOR_ACTION_BOUNDARY_INVALID", "Operator action is unavailable."
            ) from None
        if (
            not stat.S_ISDIR(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or observed.st_uid != self.root_uid
            or observed.st_gid != self.service_gid
            or stat.S_IMODE(observed.st_mode) != 0o750
        ):
            raise fail(
                "OPERATOR_ACTION_BOUNDARY_INVALID", "Operator action is unavailable."
            )
        return path

    def _replace_request(self, root: Path, request: DeploymentActionRequest) -> None:
        temporary = root / f".{request.task_identity}.request"
        content = canonical_json_bytes(request.to_dict())
        self._write_file(
            temporary,
            content,
            uid=self.root_uid,
            gid=self.service_gid,
            mode=0o640,
        )
        try:
            os.replace(temporary, self.layout.operator_action_request)
            _fsync_directory(root)
        except OSError:
            try:
                temporary.unlink()
            except OSError:
                pass
            raise fail(
                "OPERATOR_ACTION_BOUNDARY_INVALID", "Operator action is unavailable."
            ) from None

    def _prepare_receipt(self, root: Path) -> None:
        receipt = self.layout.operator_action_receipt
        try:
            receipt.unlink()
        except FileNotFoundError:
            pass
        except OSError:
            raise fail(
                "OPERATOR_ACTION_BOUNDARY_INVALID", "Operator action is unavailable."
            ) from None
        self._write_file(
            receipt,
            b"",
            uid=self.service_uid,
            gid=self.service_gid,
            mode=0o600,
        )
        _fsync_directory(root)

    def _read_receipt(
        self, request: DeploymentActionRequest
    ) -> DeploymentActionReceipt:
        path = self.layout.operator_action_receipt
        descriptor = -1
        try:
            descriptor = os.open(
                path,
                os.O_RDONLY
                | getattr(os, "O_CLOEXEC", 0)
                | getattr(os, "O_NOFOLLOW", 0),
            )
            before = os.fstat(descriptor)
            if (
                not stat.S_ISREG(before.st_mode)
                or before.st_nlink != 1
                or before.st_uid != self.service_uid
                or before.st_gid != self.service_gid
                or stat.S_IMODE(before.st_mode) != 0o600
                or not 0 < before.st_size <= _MAX_DOCUMENT_BYTES
            ):
                raise OSError
            content = _read_bounded(descriptor, _MAX_DOCUMENT_BYTES)
            after = os.fstat(descriptor)
            path_after = path.stat(follow_symlinks=False)
            if (
                _witness(before) != _witness(after)
                or (after.st_dev, after.st_ino)
                != (path_after.st_dev, path_after.st_ino)
                or len(content) != before.st_size
            ):
                raise OSError
        except OSError:
            raise fail(
                "OPERATOR_ACTION_RECEIPT_INVALID", "Operator action receipt is invalid."
            ) from None
        finally:
            if descriptor >= 0:
                os.close(descriptor)
        try:
            raw = json.loads(content, object_pairs_hook=_unique_object)
            receipt = DeploymentActionReceipt.from_dict(raw)
        except (DeploymentError, UnicodeDecodeError, ValueError, json.JSONDecodeError):
            raise fail(
                "OPERATOR_ACTION_RECEIPT_INVALID", "Operator action receipt is invalid."
            ) from None
        if canonical_json_bytes(receipt.to_dict()) != content:
            raise fail(
                "OPERATOR_ACTION_RECEIPT_INVALID", "Operator action receipt is invalid."
            )
        return receipt

    @staticmethod
    def _write_file(
        path: Path, content: bytes, *, uid: int, gid: int, mode: int
    ) -> None:
        descriptor = -1
        try:
            descriptor = os.open(
                path,
                os.O_WRONLY
                | os.O_CREAT
                | os.O_EXCL
                | getattr(os, "O_CLOEXEC", 0)
                | getattr(os, "O_NOFOLLOW", 0),
                0o600,
            )
            os.fchown(descriptor, uid, gid)
            _write_all(descriptor, content)
            os.fchmod(descriptor, mode)
            os.fsync(descriptor)
        except OSError:
            raise fail(
                "OPERATOR_ACTION_BOUNDARY_INVALID", "Operator action is unavailable."
            ) from None
        finally:
            if descriptor >= 0:
                os.close(descriptor)

    @staticmethod
    def _run_systemctl(argv: tuple[str, ...]) -> int:
        return _run_fixed_systemctl(argv, ACTION_UNIT)


class SystemdDatabasePreparer:
    """Run the only DB-writing candidate boundary after root backup."""

    def __init__(
        self,
        layout: DeploymentLayout,
        *,
        service_uid: int,
        service_gid: int,
        root_uid: int = 0,
        command_runner: ActionCommandRunner | None = None,
    ) -> None:
        if any(
            not isinstance(value, int) or isinstance(value, bool) or value < 0
            for value in (service_uid, service_gid, root_uid)
        ):
            raise fail(
                "DATABASE_PREPARE_BOUNDARY_INVALID",
                "Database preparation is unavailable.",
            )
        self.layout = layout
        self.service_uid = service_uid
        self.service_gid = service_gid
        self.root_uid = root_uid
        self.command_runner = (
            self._run_systemctl if command_runner is None else command_runner
        )

    def prepare(self, request: DatabasePrepareRequest) -> DatabasePrepareReceipt:
        if not isinstance(request, DatabasePrepareRequest):
            raise fail(
                "DATABASE_PREPARE_REQUEST_INVALID",
                "Database preparation request is invalid.",
            )
        root = self.layout.database_prepare_root
        try:
            observed = root.lstat()
        except OSError:
            raise fail(
                "DATABASE_PREPARE_BOUNDARY_INVALID",
                "Database preparation is unavailable.",
            ) from None
        if (
            not stat.S_ISDIR(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or observed.st_uid != self.root_uid
            or observed.st_gid != self.service_gid
            or stat.S_IMODE(observed.st_mode) != 0o710
        ):
            raise fail(
                "DATABASE_PREPARE_BOUNDARY_INVALID",
                "Database preparation is unavailable.",
            )
        temporary = root / f".{request.task_identity}.prepare"
        SystemdDeploymentActionRunner._write_file(
            temporary,
            canonical_json_bytes(request.to_dict()),
            uid=self.root_uid,
            gid=self.service_gid,
            mode=0o640,
        )
        try:
            os.replace(temporary, self.layout.database_prepare_request)
            try:
                self.layout.database_prepare_receipt.unlink()
            except FileNotFoundError:
                pass
            SystemdDeploymentActionRunner._write_file(
                self.layout.database_prepare_receipt,
                b"",
                uid=self.service_uid,
                gid=self.service_gid,
                mode=0o600,
            )
            _fsync_directory(root)
        except OSError:
            raise fail(
                "DATABASE_PREPARE_BOUNDARY_INVALID",
                "Database preparation is unavailable.",
            ) from None
        try:
            return_code = self.command_runner(
                (
                    str(SYSTEMCTL),
                    "--no-ask-password",
                    "start",
                    "--",
                    DATABASE_PREPARE_UNIT,
                )
            )
        except Exception:
            raise fail(
                "DATABASE_PREPARE_FAILED",
                "Database preparation failed.",
                recoverable=True,
            ) from None
        if return_code != 0:
            raise fail(
                "DATABASE_PREPARE_FAILED",
                "Database preparation failed.",
                recoverable=True,
            )
        raw = _read_service_document(
            self.layout.database_prepare_receipt,
            uid=self.service_uid,
            gid=self.service_gid,
        )
        try:
            receipt = DatabasePrepareReceipt.from_dict(raw)
        except DeploymentError:
            raise fail(
                "DATABASE_PREPARE_RECEIPT_INVALID",
                "Database preparation receipt is invalid.",
            ) from None
        if receipt.request_identity != request.identity:
            raise fail(
                "DATABASE_PREPARE_RECEIPT_INVALID",
                "Database preparation receipt is invalid.",
            )
        return receipt

    @staticmethod
    def _run_systemctl(argv: tuple[str, ...]) -> int:
        return _run_fixed_systemctl(argv, DATABASE_PREPARE_UNIT)


class SystemdEncodeRuntimePreparer:
    """Materialize one offline ENCODE runtime through a fixed unprivileged unit."""

    def __init__(
        self,
        layout: DeploymentLayout,
        *,
        service_uid: int,
        service_gid: int,
        root_uid: int = 0,
        root_gid: int = 0,
        command_runner: ActionCommandRunner | None = None,
        unit_quiescer: UnitQuiescer | None = None,
    ) -> None:
        if any(
            not isinstance(value, int) or isinstance(value, bool) or value < 0
            for value in (service_uid, service_gid, root_uid, root_gid)
        ):
            raise fail(
                "ENCODE_RUNTIME_PREPARE_BOUNDARY_INVALID",
                "ENCODE runtime preparation is unavailable.",
            )
        self.layout = layout
        self.service_uid = service_uid
        self.service_gid = service_gid
        self.root_uid = root_uid
        self.root_gid = root_gid
        self.command_runner = (
            self._run_systemctl if command_runner is None else command_runner
        )
        self.unit_quiescer = (
            self._quiesce_systemd_unit if unit_quiescer is None else unit_quiescer
        )
        if not callable(self.command_runner) or not callable(self.unit_quiescer):
            raise fail(
                "ENCODE_RUNTIME_PREPARE_BOUNDARY_INVALID",
                "ENCODE runtime preparation is unavailable.",
            )

    def prepare(
        self, request: EncodeRuntimePrepareRequest
    ) -> EncodeRuntimePrepareReceipt:
        if not isinstance(request, EncodeRuntimePrepareRequest):
            raise fail(
                "ENCODE_RUNTIME_PREPARE_REQUEST_INVALID",
                "ENCODE runtime preparation request is invalid.",
            )
        boundary = self._require_boundary()
        destination = self.layout.encode_runtime_active_root(
            request.deployment_identity
        )
        reused = self._reuse_materialization(destination, request)
        if reused is not None:
            return reused
        self._prepare_destination(destination)
        unit_dispatched = False
        try:
            self._replace_request(boundary, request)
            self._prepare_receipt(boundary)
            unit_dispatched = True
            return_code = self.command_runner(
                (
                    str(SYSTEMCTL),
                    "--no-ask-password",
                    "start",
                    "--",
                    ENCODE_RUNTIME_PREPARE_UNIT,
                )
            )
            if return_code != 0:
                raise fail(
                    "ENCODE_RUNTIME_PREPARE_FAILED",
                    "ENCODE runtime preparation failed.",
                    recoverable=True,
                )
            raw = _read_service_document(
                self.layout.encode_runtime_prepare_receipt,
                uid=self.service_uid,
                gid=self.service_gid,
                code="ENCODE_RUNTIME_PREPARE_RECEIPT_INVALID",
                message="ENCODE runtime preparation receipt is invalid.",
            )
            receipt = EncodeRuntimePrepareReceipt.from_dict(raw)
            if (
                receipt.request_identity != request.identity
                or receipt.deployment_identity != request.deployment_identity
            ):
                raise fail(
                    "ENCODE_RUNTIME_PREPARE_RECEIPT_INVALID",
                    "ENCODE runtime preparation receipt is invalid.",
                )
            inventory = self._read_inventory(
                destination,
                receipt,
                owner_uid=self.service_uid,
                owner_gid=self.service_gid,
                mode=0o600,
            )
            self._verify_and_freeze(destination, inventory)
            self._record_materialization(receipt)
        except Exception:
            if unit_dispatched and not self.unit_quiescer():
                raise fail(
                    "ENCODE_RUNTIME_RECOVERY_REQUIRED",
                    "ENCODE runtime materialization requires recovery.",
                    recoverable=True,
                ) from None
            self._quarantine_failed(destination, request)
            raise
        return receipt

    def _require_boundary(self) -> Path:
        boundaries = (
            (
                self.layout.encode_runtime_prepare_root,
                self.root_uid,
                self.service_gid,
                0o750,
            ),
            (
                self.layout.encode_runtime_materialized,
                self.root_uid,
                self.root_gid,
                0o555,
            ),
            (
                self.layout.data_root / "operator" / "runtime-materializations",
                self.root_uid,
                self.root_gid,
                0o700,
            ),
        )
        try:
            observed = tuple(
                (path, path.lstat(), uid, gid, mode)
                for path, uid, gid, mode in boundaries
            )
        except OSError:
            raise fail(
                "ENCODE_RUNTIME_PREPARE_BOUNDARY_INVALID",
                "ENCODE runtime preparation is unavailable.",
            ) from None
        for _path, witness, uid, gid, mode in observed:
            if (
                not stat.S_ISDIR(witness.st_mode)
                or stat.S_ISLNK(witness.st_mode)
                or witness.st_uid != uid
                or witness.st_gid != gid
                or stat.S_IMODE(witness.st_mode) != mode
            ):
                raise fail(
                    "ENCODE_RUNTIME_PREPARE_BOUNDARY_INVALID",
                    "ENCODE runtime preparation is unavailable.",
                )
        return self.layout.encode_runtime_prepare_root

    def _prepare_destination(self, destination: Path) -> None:
        try:
            destination.lstat()
        except FileNotFoundError:
            pass
        except OSError:
            raise fail(
                "ENCODE_RUNTIME_PREPARE_BOUNDARY_INVALID",
                "ENCODE runtime preparation is unavailable.",
            ) from None
        else:
            raise fail(
                "ENCODE_RUNTIME_ALREADY_MATERIALIZED",
                "ENCODE runtime is already materialized.",
                recoverable=True,
            )
        try:
            os.mkdir(destination, 0o700)
            os.chown(destination, self.service_uid, self.service_gid)
            os.chmod(destination, 0o700)
            _fsync_directory(self.layout.encode_runtime_materialized)
        except OSError:
            raise fail(
                "ENCODE_RUNTIME_PREPARE_BOUNDARY_INVALID",
                "ENCODE runtime preparation is unavailable.",
                recoverable=True,
            ) from None

    def _reuse_materialization(
        self,
        destination: Path,
        request: EncodeRuntimePrepareRequest,
    ) -> EncodeRuntimePrepareReceipt | None:
        try:
            destination.lstat()
        except FileNotFoundError:
            return None
        except OSError:
            raise fail(
                "ENCODE_RUNTIME_PREPARE_BOUNDARY_INVALID",
                "ENCODE runtime preparation is unavailable.",
            ) from None
        record = self.layout.encode_runtime_materialization_receipt(
            request.deployment_identity
        )
        descriptor = -1
        try:
            descriptor = os.open(
                record,
                os.O_RDONLY
                | getattr(os, "O_CLOEXEC", 0)
                | getattr(os, "O_NOFOLLOW", 0),
            )
            before = os.fstat(descriptor)
            if (
                not stat.S_ISREG(before.st_mode)
                or before.st_nlink != 1
                or before.st_uid != self.root_uid
                or before.st_gid != self.root_gid
                or stat.S_IMODE(before.st_mode) != 0o400
                or not 0 < before.st_size <= _MAX_DOCUMENT_BYTES
            ):
                raise OSError
            content = _read_bounded(descriptor, _MAX_DOCUMENT_BYTES)
            after = os.fstat(descriptor)
            path_after = record.stat(follow_symlinks=False)
            if (
                _witness(before) != _witness(after)
                or (after.st_dev, after.st_ino)
                != (path_after.st_dev, path_after.st_ino)
                or len(content) != before.st_size
            ):
                raise OSError
            raw = json.loads(content, object_pairs_hook=_unique_object)
            receipt = EncodeRuntimePrepareReceipt.from_dict(raw)
            if (
                canonical_json_bytes(receipt.to_dict()) != content
                or receipt.deployment_identity != request.deployment_identity
            ):
                raise OSError
            inventory = self._read_inventory(
                destination,
                receipt,
                owner_uid=self.root_uid,
                owner_gid=self.root_gid,
                mode=0o444,
            )
            self._verify_frozen_tree(destination, inventory)
        except (DeploymentError, OSError, UnicodeDecodeError, ValueError):
            raise fail(
                "ENCODE_RUNTIME_RECOVERY_REQUIRED",
                "ENCODE runtime materialization requires recovery.",
                recoverable=True,
            ) from None
        finally:
            if descriptor >= 0:
                os.close(descriptor)
        return EncodeRuntimePrepareReceipt.create(
            request_identity=request.identity,
            deployment_identity=request.deployment_identity,
            inventory=inventory,
        )

    def _record_materialization(self, receipt: EncodeRuntimePrepareReceipt) -> None:
        record = self.layout.encode_runtime_materialization_receipt(
            receipt.deployment_identity
        )
        SystemdDeploymentActionRunner._write_file(
            record,
            canonical_json_bytes(receipt.to_dict()),
            uid=self.root_uid,
            gid=self.root_gid,
            mode=0o400,
        )
        _fsync_directory(record.parent)

    def _read_inventory(
        self,
        root: Path,
        receipt: EncodeRuntimePrepareReceipt,
        *,
        owner_uid: int,
        owner_gid: int,
        mode: int,
    ) -> EncodeRuntimeInventory:
        try:
            inventory, content = _read_runtime_inventory_document(
                root,
                owner_uid=owner_uid,
                owner_gid=owner_gid,
                mode=mode,
            )
            if (
                len(content) != receipt.inventory_size
                or hashlib.sha256(content).hexdigest() != receipt.inventory_sha256
                or inventory.tree_identity != receipt.tree_identity
                or len(inventory.entries) != receipt.entry_count
            ):
                raise OSError
            return inventory
        except (DeploymentError, OSError, UnicodeDecodeError, ValueError):
            raise fail(
                "ENCODE_RUNTIME_PREPARE_RECEIPT_INVALID",
                "ENCODE runtime preparation receipt is invalid.",
            ) from None

    def _replace_request(
        self, boundary: Path, request: EncodeRuntimePrepareRequest
    ) -> None:
        temporary = boundary / f".{request.task_identity}.request"
        try:
            SystemdDeploymentActionRunner._write_file(
                temporary,
                canonical_json_bytes(request.to_dict()),
                uid=self.root_uid,
                gid=self.service_gid,
                mode=0o640,
            )
            os.replace(temporary, self.layout.encode_runtime_prepare_request)
            _fsync_directory(boundary)
        except (DeploymentError, OSError):
            try:
                temporary.unlink()
            except OSError:
                pass
            raise fail(
                "ENCODE_RUNTIME_PREPARE_BOUNDARY_INVALID",
                "ENCODE runtime preparation is unavailable.",
                recoverable=True,
            ) from None

    def _prepare_receipt(self, boundary: Path) -> None:
        receipt = self.layout.encode_runtime_prepare_receipt
        try:
            receipt.unlink()
        except FileNotFoundError:
            pass
        except OSError:
            raise fail(
                "ENCODE_RUNTIME_PREPARE_BOUNDARY_INVALID",
                "ENCODE runtime preparation is unavailable.",
            ) from None
        SystemdDeploymentActionRunner._write_file(
            receipt,
            b"",
            uid=self.service_uid,
            gid=self.service_gid,
            mode=0o600,
        )
        _fsync_directory(boundary)

    def _verify_and_freeze(self, root: Path, inventory: EncodeRuntimeInventory) -> None:
        expected = {item.path: item for item in inventory.entries}
        expected_directories = {
            parent.as_posix()
            for item in inventory.entries
            for parent in Path(item.path).parents
            if parent != Path(".")
        }
        observed_entries: set[str] = set()
        observed_directories: set[str] = set()
        try:
            for current, directories, filenames in os.walk(
                root, topdown=True, followlinks=False
            ):
                current_path = Path(current)
                relative_root = current_path.relative_to(root)
                if relative_root != Path("."):
                    observed_directories.add(relative_root.as_posix())
                for name in tuple(directories):
                    child = current_path / name
                    relative = child.relative_to(root).as_posix()
                    witness = child.lstat()
                    if stat.S_ISLNK(witness.st_mode):
                        directories.remove(name)
                        item = expected.get(relative)
                        if item is None or item.kind != "symlink":
                            raise OSError
                        self._verify_symlink(
                            root,
                            child,
                            item,
                            owner_uid=self.service_uid,
                            owner_gid=self.service_gid,
                        )
                        os.lchown(child, self.root_uid, self.root_gid)
                        observed_entries.add(relative)
                    elif (
                        not stat.S_ISDIR(witness.st_mode)
                        or witness.st_uid != self.service_uid
                        or witness.st_gid != self.service_gid
                        or stat.S_IMODE(witness.st_mode) & 0o077
                    ):
                        raise OSError
                for name in filenames:
                    child = current_path / name
                    relative = child.relative_to(root).as_posix()
                    if relative == ".helixweave-runtime-inventory.json":
                        continue
                    item = expected.get(relative)
                    if item is None:
                        raise OSError
                    witness = child.lstat()
                    if item.kind == "symlink":
                        if not stat.S_ISLNK(witness.st_mode):
                            raise OSError
                        self._verify_symlink(
                            root,
                            child,
                            item,
                            owner_uid=self.service_uid,
                            owner_gid=self.service_gid,
                        )
                        os.lchown(child, self.root_uid, self.root_gid)
                    else:
                        self._verify_and_freeze_file(child, item)
                    observed_entries.add(relative)
            if (
                observed_entries != set(expected)
                or observed_directories != expected_directories
            ):
                raise OSError
            self._verify_resolved_symlinks(root, inventory.entries)
            for relative in sorted(
                expected_directories,
                key=lambda value: (value.count("/"), value),
                reverse=True,
            ):
                directory = root / relative
                witness = directory.lstat()
                if (
                    not stat.S_ISDIR(witness.st_mode)
                    or stat.S_ISLNK(witness.st_mode)
                    or witness.st_uid != self.service_uid
                    or witness.st_gid != self.service_gid
                    or stat.S_IMODE(witness.st_mode) & 0o077
                ):
                    raise OSError
                os.chown(directory, self.root_uid, self.root_gid)
                os.chmod(directory, 0o555)
                _fsync_directory(directory)
            inventory_path = root / ".helixweave-runtime-inventory.json"
            os.chown(inventory_path, self.root_uid, self.root_gid)
            os.chmod(inventory_path, 0o444)
            os.chown(root, self.root_uid, self.root_gid)
            os.chmod(root, 0o555)
            _fsync_directory(root)
            _fsync_directory(self.layout.encode_runtime_materialized)
        except OSError:
            raise fail(
                "ENCODE_RUNTIME_PREPARE_RECEIPT_INVALID",
                "ENCODE runtime preparation receipt is invalid.",
                recoverable=True,
            ) from None

    def _verify_and_freeze_file(self, path: Path, item: EncodeRuntimeEntry) -> None:
        descriptor = -1
        try:
            descriptor = os.open(
                path,
                os.O_RDONLY
                | getattr(os, "O_CLOEXEC", 0)
                | getattr(os, "O_NOFOLLOW", 0),
            )
            before = os.fstat(descriptor)
            if (
                not stat.S_ISREG(before.st_mode)
                or before.st_nlink != 1
                or before.st_uid != self.service_uid
                or before.st_gid != self.service_gid
                or stat.S_IMODE(before.st_mode) & 0o077
                or before.st_size != item.size
            ):
                raise OSError
            digest = hashlib.sha256()
            remaining = item.size
            while remaining:
                chunk = os.read(descriptor, min(1024 * 1024, remaining))
                if not chunk:
                    raise OSError
                digest.update(chunk)
                remaining -= len(chunk)
            if os.read(descriptor, 1) or digest.hexdigest() != item.sha256:
                raise OSError
            after_read = os.fstat(descriptor)
            if _witness(before) != _witness(after_read):
                raise OSError
            os.fchown(descriptor, self.root_uid, self.root_gid)
            if item.mode is None:
                raise OSError
            os.fchmod(descriptor, item.mode)
            os.fsync(descriptor)
            path_after = path.stat(follow_symlinks=False)
            after = os.fstat(descriptor)
            if (after.st_dev, after.st_ino) != (
                path_after.st_dev,
                path_after.st_ino,
            ):
                raise OSError
        except OSError:
            raise
        finally:
            if descriptor >= 0:
                os.close(descriptor)

    def _verify_frozen_tree(
        self, root: Path, inventory: EncodeRuntimeInventory
    ) -> None:
        expected = {item.path: item for item in inventory.entries}
        expected_directories = {
            parent.as_posix()
            for item in inventory.entries
            for parent in Path(item.path).parents
            if parent != Path(".")
        }
        observed_entries: set[str] = set()
        observed_directories: set[str] = set()
        root_witness = root.lstat()
        if (
            not stat.S_ISDIR(root_witness.st_mode)
            or stat.S_ISLNK(root_witness.st_mode)
            or root_witness.st_uid != self.root_uid
            or root_witness.st_gid != self.root_gid
            or stat.S_IMODE(root_witness.st_mode) != 0o555
        ):
            raise OSError
        for current, directories, filenames in os.walk(
            root, topdown=True, followlinks=False
        ):
            current_path = Path(current)
            relative_root = current_path.relative_to(root)
            if relative_root != Path("."):
                observed_directories.add(relative_root.as_posix())
            for name in tuple(directories):
                child = current_path / name
                relative = child.relative_to(root).as_posix()
                witness = child.lstat()
                if stat.S_ISLNK(witness.st_mode):
                    directories.remove(name)
                    item = expected.get(relative)
                    if item is None or item.kind != "symlink":
                        raise OSError
                    self._verify_symlink(
                        root,
                        child,
                        item,
                        owner_uid=self.root_uid,
                        owner_gid=self.root_gid,
                    )
                    observed_entries.add(relative)
                elif (
                    not stat.S_ISDIR(witness.st_mode)
                    or witness.st_uid != self.root_uid
                    or witness.st_gid != self.root_gid
                    or stat.S_IMODE(witness.st_mode) != 0o555
                ):
                    raise OSError
            for name in filenames:
                path = current_path / name
                relative = path.relative_to(root).as_posix()
                if relative == ".helixweave-runtime-inventory.json":
                    continue
                item = expected.get(relative)
                if item is None:
                    raise OSError
                if item.kind == "symlink":
                    self._verify_symlink(
                        root,
                        path,
                        item,
                        owner_uid=self.root_uid,
                        owner_gid=self.root_gid,
                    )
                else:
                    self._verify_frozen_file(path, item)
                observed_entries.add(relative)
        if (
            observed_entries != set(expected)
            or observed_directories != expected_directories
        ):
            raise OSError
        self._verify_resolved_symlinks(root, inventory.entries)

    def _verify_frozen_file(self, path: Path, item: EncodeRuntimeEntry) -> None:
        descriptor = -1
        try:
            descriptor = os.open(
                path,
                os.O_RDONLY
                | getattr(os, "O_CLOEXEC", 0)
                | getattr(os, "O_NOFOLLOW", 0),
            )
            before = os.fstat(descriptor)
            if (
                not stat.S_ISREG(before.st_mode)
                or before.st_nlink != 1
                or before.st_uid != self.root_uid
                or before.st_gid != self.root_gid
                or item.mode is None
                or stat.S_IMODE(before.st_mode) != item.mode
                or before.st_size != item.size
            ):
                raise OSError
            digest = hashlib.sha256()
            remaining = item.size
            while remaining:
                chunk = os.read(descriptor, min(1024 * 1024, remaining))
                if not chunk:
                    raise OSError
                digest.update(chunk)
                remaining -= len(chunk)
            if os.read(descriptor, 1) or digest.hexdigest() != item.sha256:
                raise OSError
            after = os.fstat(descriptor)
            path_after = path.stat(follow_symlinks=False)
            if _witness(before) != _witness(after) or (after.st_dev, after.st_ino) != (
                path_after.st_dev,
                path_after.st_ino,
            ):
                raise OSError
        finally:
            if descriptor >= 0:
                os.close(descriptor)

    @staticmethod
    def _verify_symlink(
        root: Path,
        path: Path,
        item: EncodeRuntimeEntry,
        *,
        owner_uid: int,
        owner_gid: int,
    ) -> None:
        witness = path.lstat()
        target = os.readlink(path)
        if (
            item.kind != "symlink"
            or item.target is None
            or not stat.S_ISLNK(witness.st_mode)
            or witness.st_uid != owner_uid
            or witness.st_gid != owner_gid
            or target != item.target
            or not _safe_runtime_symlink(item.path, target)
            or hashlib.sha256(target.encode("utf-8")).hexdigest() != item.sha256
            or len(target.encode("utf-8")) != item.size
            or path.parent.resolve(strict=True).is_relative_to(
                root.resolve(strict=True)
            )
            is False
        ):
            raise OSError

    @staticmethod
    def _verify_resolved_symlinks(
        root: Path, entries: Sequence[EncodeRuntimeEntry]
    ) -> None:
        resolved_root = root.resolve(strict=True)
        for item in entries:
            if item.kind != "symlink":
                continue
            resolved = (root / item.path).resolve(strict=True)
            if resolved != resolved_root and not resolved.is_relative_to(resolved_root):
                raise OSError

    def _quarantine_failed(
        self,
        destination: Path,
        request: EncodeRuntimePrepareRequest,
    ) -> None:
        failed = self.layout.encode_runtime_failed(
            request.deployment_identity, request.task_identity
        )
        try:
            failed.lstat()
        except FileNotFoundError:
            pass
        except OSError:
            raise fail(
                "ENCODE_RUNTIME_RECOVERY_REQUIRED",
                "ENCODE runtime materialization requires recovery.",
                recoverable=True,
            ) from None
        else:
            raise fail(
                "ENCODE_RUNTIME_RECOVERY_REQUIRED",
                "ENCODE runtime materialization requires recovery.",
                recoverable=True,
            )
        try:
            os.rename(destination, failed)
            os.chown(failed, self.root_uid, self.root_gid)
            os.chmod(failed, 0o500)
            _fsync_directory(self.layout.encode_runtime_materialized)
        except OSError:
            raise fail(
                "ENCODE_RUNTIME_RECOVERY_REQUIRED",
                "ENCODE runtime materialization requires recovery.",
                recoverable=True,
            ) from None

    @staticmethod
    def _run_systemctl(argv: tuple[str, ...]) -> int:
        return _run_fixed_systemctl(argv, ENCODE_RUNTIME_PREPARE_UNIT)

    @staticmethod
    def _quiesce_systemd_unit() -> bool:
        return _quiesce_fixed_systemd_unit(ENCODE_RUNTIME_PREPARE_UNIT)


@dataclass(frozen=True)
class MaterializedEncodeRuntime:
    deployment_identity: str
    tree_identity: str
    entry_count: int


def verify_materialized_encode_runtime(
    layout: DeploymentLayout,
    deployment_identity: str,
    *,
    root_uid: int = 0,
    root_gid: int = 0,
) -> MaterializedEncodeRuntime:
    """Read-only verification of one frozen offline scientific runtime."""

    if not isinstance(layout, DeploymentLayout):
        raise fail(
            "ENCODE_RUNTIME_MATERIALIZATION_INVALID",
            "ENCODE runtime materialization is invalid.",
        )
    _identity(deployment_identity, code="ENCODE_RUNTIME_MATERIALIZATION_INVALID")
    root = layout.encode_runtime_active_root(deployment_identity)
    try:
        inventory, _content = _read_runtime_inventory_document(
            root,
            owner_uid=root_uid,
            owner_gid=root_gid,
            mode=0o444,
        )
        verifier = SystemdEncodeRuntimePreparer(
            layout,
            service_uid=root_uid,
            service_gid=root_gid,
            root_uid=root_uid,
            root_gid=root_gid,
            command_runner=lambda _argv: 1,
        )
        verifier._verify_frozen_tree(root, inventory)
        source = (
            layout.encode_runtimes
            / deployment_identity
            / "payload"
            / "contracts"
            / "encode-runtime"
        )
        index_content = _read_owned_bounded_file(
            source / "package-index.json",
            maximum=32 * 1024 * 1024,
            owner_uid=root_uid,
            owner_gid=root_gid,
            mode=0o444,
        )
        from encode_pipeline.deployment.native_contracts import (
            parse_encode_runtime_index,
        )

        index = parse_encode_runtime_index(index_content)
        entries = {item.path: item for item in inventory.entries}
        for executable in (
            "mamba-root/bin/activate",
            "runner/bin/conda",
            "runner/bin/snakemake",
            "runner/libexec/micromamba",
        ):
            item = entries.get(executable)
            if item is None or item.kind != "file" or item.mode != 0o555:
                raise OSError
        conda_prefix = root / "conda-envs"
        for environment in index.environments:
            environment_path = environment["environment_path"]
            if environment_path == "workflow/envs/runner.yml":
                continue
            environment_content = _read_owned_bounded_file(
                source / str(environment_path),
                maximum=_MAX_DOCUMENT_BYTES,
                owner_uid=root_uid,
                owner_gid=root_gid,
                mode=0o444,
            )
            full_hash = snakemake_environment_hash(
                conda_prefix,
                environment_content,
            )
            marker = entries.get(f"conda-envs/{full_hash}.env_setup_done")
            prefix = f"conda-envs/{full_hash}/"
            if (
                marker is None
                or marker.kind != "file"
                or marker.mode != 0o444
                or marker.size != 0
                or not any(path.startswith(prefix) for path in entries)
            ):
                raise OSError
    except Exception:
        raise fail(
            "ENCODE_RUNTIME_MATERIALIZATION_INVALID",
            "ENCODE runtime materialization is invalid.",
        ) from None
    return MaterializedEncodeRuntime(
        deployment_identity=deployment_identity,
        tree_identity=inventory.tree_identity,
        entry_count=len(inventory.entries),
    )


def _read_runtime_inventory_document(
    root: Path,
    *,
    owner_uid: int,
    owner_gid: int,
    mode: int,
) -> tuple[EncodeRuntimeInventory, bytes]:
    content = _read_owned_bounded_file(
        root / ".helixweave-runtime-inventory.json",
        maximum=_MAX_RUNTIME_INVENTORY_BYTES,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        mode=mode,
    )
    raw = json.loads(content, object_pairs_hook=_unique_object)
    inventory = EncodeRuntimeInventory.from_dict(raw)
    if inventory.to_bytes() != content:
        raise OSError
    return inventory, content


def _read_owned_bounded_file(
    path: Path,
    *,
    maximum: int,
    owner_uid: int,
    owner_gid: int,
    mode: int,
) -> bytes:
    descriptor = -1
    try:
        descriptor = os.open(
            path,
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
        )
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_uid != owner_uid
            or before.st_gid != owner_gid
            or stat.S_IMODE(before.st_mode) != mode
            or not 0 < before.st_size <= maximum
        ):
            raise OSError
        content = _read_bounded(descriptor, maximum)
        after = os.fstat(descriptor)
        path_after = path.stat(follow_symlinks=False)
        if (
            _witness(before) != _witness(after)
            or (after.st_dev, after.st_ino) != (path_after.st_dev, path_after.st_ino)
            or len(content) != before.st_size
        ):
            raise OSError
        return content
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _run_fixed_systemctl(argv: tuple[str, ...], unit: str) -> int:
    expected = (
        str(SYSTEMCTL),
        "--no-ask-password",
        "start",
        "--",
        unit,
    )
    if argv != expected:
        raise fail("OPERATOR_ACTION_COMMAND_INVALID", "Operator action is unavailable.")
    try:
        with tempfile.TemporaryFile() as stdout, tempfile.TemporaryFile() as stderr:
            completed = subprocess.run(
                argv,
                stdin=subprocess.DEVNULL,
                stdout=stdout,
                stderr=stderr,
                env=_SAFE_ENVIRONMENT,
                close_fds=True,
                timeout=(
                    _ENCODE_RUNTIME_UNIT_TIMEOUT_SECONDS
                    if unit == ENCODE_RUNTIME_PREPARE_UNIT
                    else _ACTION_TIMEOUT_SECONDS
                ),
                check=False,
            )
            if stdout.tell() or stderr.tell():
                return 1
            return completed.returncode
    except (OSError, subprocess.SubprocessError):
        return 1


def _quiesce_fixed_systemd_unit(unit: str) -> bool:
    if unit != ENCODE_RUNTIME_PREPARE_UNIT:
        return False
    stop = (
        str(SYSTEMCTL),
        "--no-ask-password",
        "stop",
        "--",
        unit,
    )
    show = (
        str(SYSTEMCTL),
        "--no-ask-password",
        "show",
        "--property=ActiveState",
        "--value",
        "--",
        unit,
    )
    try:
        with tempfile.TemporaryFile() as stdout, tempfile.TemporaryFile() as stderr:
            stopped = subprocess.run(
                stop,
                stdin=subprocess.DEVNULL,
                stdout=stdout,
                stderr=stderr,
                env=_SAFE_ENVIRONMENT,
                close_fds=True,
                timeout=_UNIT_QUIESCE_TIMEOUT_SECONDS,
                check=False,
            )
            if stopped.returncode != 0 or stdout.tell() or stderr.tell():
                return False
        with tempfile.TemporaryFile() as stdout, tempfile.TemporaryFile() as stderr:
            observed = subprocess.run(
                show,
                stdin=subprocess.DEVNULL,
                stdout=stdout,
                stderr=stderr,
                env=_SAFE_ENVIRONMENT,
                close_fds=True,
                timeout=_UNIT_QUIESCE_TIMEOUT_SECONDS,
                check=False,
            )
            if observed.returncode != 0 or stderr.tell() or stdout.tell() > 16:
                return False
            stdout.seek(0)
            return stdout.read() in {b"failed\n", b"inactive\n"}
    except (OSError, subprocess.SubprocessError):
        return False


def _read_service_document(
    path: Path,
    *,
    uid: int,
    gid: int,
    code: str = "DATABASE_PREPARE_RECEIPT_INVALID",
    message: str = "Database preparation receipt is invalid.",
) -> dict[str, object]:
    descriptor = -1
    try:
        descriptor = os.open(
            path,
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
        )
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_uid != uid
            or before.st_gid != gid
            or stat.S_IMODE(before.st_mode) != 0o600
            or not 0 < before.st_size <= _MAX_DOCUMENT_BYTES
        ):
            raise OSError
        content = _read_bounded(descriptor, _MAX_DOCUMENT_BYTES)
        after = os.fstat(descriptor)
        path_after = path.stat(follow_symlinks=False)
        if (
            _witness(before) != _witness(after)
            or (after.st_dev, after.st_ino) != (path_after.st_dev, path_after.st_ino)
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
    except (UnicodeDecodeError, ValueError, json.JSONDecodeError):
        raise fail(code, message) from None
    if not isinstance(raw, dict) or canonical_json_bytes(raw) != content:
        raise fail(code, message)
    return raw


def _write_all(descriptor: int, content: bytes) -> None:
    offset = 0
    while offset < len(content):
        written = os.write(descriptor, content[offset:])
        if written <= 0:
            raise OSError
        offset += written


def _read_bounded(descriptor: int, maximum: int) -> bytes:
    chunks: list[bytes] = []
    remaining = maximum + 1
    while remaining:
        chunk = os.read(descriptor, min(65536, remaining))
        if not chunk:
            break
        chunks.append(chunk)
        remaining -= len(chunk)
    content = b"".join(chunks)
    if len(content) > maximum:
        raise OSError
    return content


def _fsync_directory(path: Path) -> None:
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


def _witness(value: os.stat_result) -> tuple[int, ...]:
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


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError("duplicate key")
        value[key] = item
    return value


__all__ = [
    "ACTION_PHASES",
    "ACTION_RECEIPT_SCHEMA",
    "ACTION_REQUEST_SCHEMA",
    "ACTION_UNIT",
    "BULK_RUNTIME_PREPARE_OPERATIONS",
    "BULK_RUNTIME_PREPARE_RECEIPT_IDENTITY_SCHEME",
    "BULK_RUNTIME_PREPARE_RECEIPT_SCHEMA",
    "BULK_RUNTIME_PREPARE_REQUEST_IDENTITY_SCHEME",
    "BULK_RUNTIME_PREPARE_REQUEST_SCHEMA",
    "COMPATIBILITY_STATUSES",
    "DATABASE_MODES",
    "DATABASE_PREPARE_RECEIPT_SCHEMA",
    "DATABASE_PREPARE_REQUEST_SCHEMA",
    "DATABASE_PREPARE_UNIT",
    "ENCODE_RUNTIME_PREPARE_RECEIPT_SCHEMA",
    "ENCODE_RUNTIME_PREPARE_REQUEST_SCHEMA",
    "ENCODE_RUNTIME_PREPARE_UNIT",
    "ENCODE_RUNTIME_INVENTORY_SCHEMA",
    "ENCODE_RUNTIME_TREE_IDENTITY_SCHEME",
    "READINESS_REASONS",
    "READINESS_STATUSES",
    "VERIFICATION_CHECKS",
    "DatabasePrepareReceipt",
    "DatabasePrepareRequest",
    "DatabasePreparer",
    "BulkRuntimePrepareReceipt",
    "BulkRuntimePrepareRequest",
    "DeploymentActionReceipt",
    "DeploymentActionRequest",
    "DeploymentActionRunner",
    "EncodeRuntimeEntry",
    "EncodeRuntimeInventory",
    "EncodeRuntimePrepareReceipt",
    "EncodeRuntimePrepareRequest",
    "EncodeRuntimePreparer",
    "ReadinessCheck",
    "MaterializedEncodeRuntime",
    "SystemdDeploymentActionRunner",
    "SystemdDatabasePreparer",
    "SystemdEncodeRuntimePreparer",
    "snakemake_environment_hash",
    "verify_materialized_encode_runtime",
]
