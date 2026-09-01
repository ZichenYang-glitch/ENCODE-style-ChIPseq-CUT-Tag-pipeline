"""Fail-closed request boundary for the privileged deployment operator.

The sudo-facing process accepts only the small grammar defined here.  Host
mutation is deliberately delegated to an injected backend so request parsing,
path derivation, and public receipts can be tested without systemd or root.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
import fcntl
import grp
import hashlib
import json
import os
from pathlib import Path
import pwd
import re
import socket
import stat
import subprocess
import sys
import tempfile
from typing import Protocol

from encode_pipeline.deployment.bundle import BundleStore
from encode_pipeline.deployment.bulk_docker_boundary import (
    BulkDockerBoundary,
    observe_bulk_docker_boundary,
)
from encode_pipeline.deployment.canonical import (
    canonical_identity,
    canonical_json_bytes,
)
from encode_pipeline.deployment.database import (
    DatabaseSchemaAdmission,
    SchemaAdmissionProvider,
    StoppedWriter,
    backup_database,
    create_write_stop_witness,
    database_content_identity,
    fresh_database_candidate_path,
    inspect_database,
    publish_fresh_database,
    quarantine_invalid_fresh_database,
    quarantine_fresh_database,
    restore_database_backup,
)
from encode_pipeline.deployment.errors import DeploymentError, DeploymentIssue, fail
from encode_pipeline.deployment.filesystem import read_regular_file
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import (
    BULK_RNASEQ_RUNTIME,
    ENCODE_RUNTIME,
    PLATFORM,
    COMPONENTS,
    DeploymentState,
)
from encode_pipeline.deployment.operator_transaction import (
    OperatorJournalHandle,
    OperatorJournalStore,
    OperatorJournalSummary,
    OperatorTransaction,
)
from encode_pipeline.deployment.operator_action import (
    BulkRuntimePrepareReceipt,
    BulkRuntimePrepareRequest,
    DatabasePrepareRequest,
    DatabasePreparer,
    DeploymentActionReceipt,
    DeploymentActionRequest,
    DeploymentActionRunner,
    EncodeRuntimePrepareRequest,
    EncodeRuntimePreparer,
    ReadinessCheck,
    SystemdDatabasePreparer,
    SystemdDeploymentActionRunner,
    SystemdEncodeRuntimePreparer,
)
from encode_pipeline.deployment.operator_boundary import (
    StableBoundaryError,
    verify_stable_operator_boundary,
)
from encode_pipeline.deployment.state import (
    MAX_PLATFORM_ENV_BYTES,
    PLATFORM_ENV_FILENAME,
    PlatformEnvironment,
    StateStore,
    render_platform_environment,
)


OPERATOR_RECEIPT_SCHEMA = "helixweave-operator-receipt-v1"
SERVICE_IDENTITY_SCHEME = "helixweave-service-identity-v2"

STAGE = "stage"
ACTIVATE = "activate"
ROLLBACK = "rollback"
START = "start"
STATUS = "status"
STOP = "stop"
CLEANUP = "cleanup"
UNINSTALL = "uninstall"
OBSERVE = "observe"
VERIFY = "verify"

OPERATIONS = (
    STAGE,
    ACTIVATE,
    ROLLBACK,
    START,
    STATUS,
    STOP,
    CLEANUP,
    UNINSTALL,
    OBSERVE,
    VERIFY,
)
SERVICE_OPERATIONS = (START, STATUS, STOP, CLEANUP)
BOUND_SERVICE_OPERATIONS = (STOP, CLEANUP)

SERVICE_UNITS = (
    "helixweave-redis.service",
    "helixweave-docker-rootless.service",
    "helixweave-api.service",
    "helixweave-worker.service",
)
WRITER_UNITS = ("helixweave-api.service", "helixweave-worker.service")
DEPENDENCY_UNITS = (
    "helixweave-redis.service",
    "helixweave-docker-rootless.service",
)
SERVICE_STOP_ORDER = (
    "helixweave-api.service",
    "helixweave-worker.service",
    "helixweave-docker-rootless.service",
    "helixweave-redis.service",
)

SERVICE_SOCKET_NAMES: dict[str, tuple[str, ...]] = {
    "helixweave-api.service": ("api-http",),
    "helixweave-worker.service": (),
    "helixweave-redis.service": ("redis-queue",),
    "helixweave-docker-rootless.service": ("bulk-docker",),
}

SYSTEMCTL = Path("/usr/bin/systemctl")
SAFE_ENVIRONMENT = {
    "LANG": "C.UTF-8",
    "LC_ALL": "C.UTF-8",
    "PATH": "/usr/sbin:/usr/bin:/sbin:/bin",
}
_SYSTEMCTL_TIMEOUT_SECONDS = 15.0
_BULK_RUNTIME_SYSTEMD_TIMEOUT_SECONDS = 14_700.0
_MAX_COMMAND_OUTPUT = 64 * 1024
_MAX_OPERATOR_DOCUMENT_BYTES = 1024 * 1024
_MAX_REFERENCE_API_BYTES = 1024 * 1024
_REFERENCE_WORKFLOW_IDS = (
    "bulk-rnaseq",
    "encode-style-chipseq-cuttag-atac-mnase",
)
_BULK_RUNTIME_PREPARE_UNIT = "helixweave-bulk-runtime-prepare.service"
OPERATOR_OBSERVATION_SCHEMA = "helixweave-operator-observation-v3"
OPERATOR_OBSERVATION_IDENTITY_SCHEME = "helixweave-operator-observation-identity-v3"
OPERATOR_PERMISSIONS_IDENTITY_SCHEME = "helixweave-operator-permissions-v1"
OPERATOR_BOUNDARY_REASONS = ("READY", "BOUNDARY_INVALID")

_CONTENT_IDENTITY = re.compile(r"^sha256-[0-9a-f]{64}$")
_TASK_IDENTITY = re.compile(r"^task-[0-9a-f]{32}$")


def _valid_identity(value: object) -> bool:
    return isinstance(value, str) and _CONTENT_IDENTITY.fullmatch(value) is not None


def _valid_task_identity(value: object) -> bool:
    return isinstance(value, str) and _TASK_IDENTITY.fullmatch(value) is not None


def _positive_integer(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool) and value > 0


@dataclass(frozen=True, order=True)
class SocketWitness:
    """Stable socket evidence without publishing a private filesystem path."""

    name: str
    device: int
    inode: int
    kernel_inode: int

    def __post_init__(self) -> None:
        if (
            self.name
            not in {name for names in SERVICE_SOCKET_NAMES.values() for name in names}
            or not _positive_integer(self.device)
            or not _positive_integer(self.inode)
            or not _positive_integer(self.kernel_inode)
        ):
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_INVALID",
                "Service identity evidence is invalid.",
            )

    def to_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "device": self.device,
            "inode": self.inode,
            "kernel_inode": self.kernel_inode,
        }


@dataclass(frozen=True)
class ServiceIdentity:
    """Identity binding one unit to its task, process, executable, and sockets."""

    identity: str
    unit: str
    deployment_identity: str
    task_identity: str
    main_pid: int
    process_start_ticks: int
    executable_device: int
    executable_inode: int
    cmdline_identity: str
    boot_identity: str
    invocation_identity: str
    cgroup_identity: str
    sockets: tuple[SocketWitness, ...]

    @classmethod
    def create(
        cls,
        *,
        unit: str,
        deployment_identity: str,
        task_identity: str,
        main_pid: int,
        process_start_ticks: int,
        executable_device: int,
        executable_inode: int,
        cmdline_identity: str,
        boot_identity: str,
        invocation_identity: str,
        cgroup_identity: str,
        sockets: Sequence[SocketWitness],
    ) -> "ServiceIdentity":
        value: dict[str, object] = {
            "unit": unit,
            "deployment_identity": deployment_identity,
            "task_identity": task_identity,
            "main_pid": main_pid,
            "process_start_ticks": process_start_ticks,
            "executable_device": executable_device,
            "executable_inode": executable_inode,
            "cmdline_identity": cmdline_identity,
            "boot_identity": boot_identity,
            "invocation_identity": invocation_identity,
            "cgroup_identity": cgroup_identity,
            "sockets": [item.to_dict() for item in sockets],
        }
        value["identity"] = canonical_identity(
            value,
            scheme=SERVICE_IDENTITY_SCHEME,
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "ServiceIdentity":
        if not isinstance(raw, dict) or set(raw) != {
            "identity",
            "unit",
            "deployment_identity",
            "task_identity",
            "main_pid",
            "process_start_ticks",
            "executable_device",
            "executable_inode",
            "cmdline_identity",
            "boot_identity",
            "invocation_identity",
            "cgroup_identity",
            "sockets",
        }:
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_INVALID",
                "Service identity evidence is invalid.",
            )
        unit = raw["unit"]
        raw_sockets = raw["sockets"]
        if unit not in SERVICE_UNITS or not isinstance(raw_sockets, list):
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_INVALID",
                "Service identity evidence is invalid.",
            )
        try:
            sockets = tuple(
                SocketWitness(
                    name=item["name"],
                    device=item["device"],
                    inode=item["inode"],
                    kernel_inode=item["kernel_inode"],
                )
                for item in raw_sockets
                if isinstance(item, dict)
                and set(item) == {"name", "device", "inode", "kernel_inode"}
            )
        except (KeyError, TypeError):
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_INVALID",
                "Service identity evidence is invalid.",
            ) from None
        if (
            len(sockets) != len(raw_sockets)
            or tuple(sorted(sockets)) != sockets
            or len({item.name for item in sockets}) != len(sockets)
            or tuple(item.name for item in sockets) != SERVICE_SOCKET_NAMES[unit]
            or not _valid_identity(raw["deployment_identity"])
            or not _valid_task_identity(raw["task_identity"])
            or not _valid_identity(raw["cmdline_identity"])
            or not _valid_identity(raw["boot_identity"])
            or not _valid_identity(raw["invocation_identity"])
            or not _valid_identity(raw["cgroup_identity"])
            or not _positive_integer(raw["main_pid"])
            or not _positive_integer(raw["process_start_ticks"])
            or not _positive_integer(raw["executable_device"])
            or not _positive_integer(raw["executable_inode"])
        ):
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_INVALID",
                "Service identity evidence is invalid.",
            )
        without_identity = {
            key: value for key, value in raw.items() if key != "identity"
        }
        expected = canonical_identity(
            without_identity,
            scheme=SERVICE_IDENTITY_SCHEME,
        )
        if raw["identity"] != expected:
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_INVALID",
                "Service identity evidence is invalid.",
            )
        return cls(
            identity=expected,
            unit=unit,
            deployment_identity=raw["deployment_identity"],
            task_identity=raw["task_identity"],
            main_pid=raw["main_pid"],
            process_start_ticks=raw["process_start_ticks"],
            executable_device=raw["executable_device"],
            executable_inode=raw["executable_inode"],
            cmdline_identity=raw["cmdline_identity"],
            boot_identity=raw["boot_identity"],
            invocation_identity=raw["invocation_identity"],
            cgroup_identity=raw["cgroup_identity"],
            sockets=sockets,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "identity": self.identity,
            "unit": self.unit,
            "deployment_identity": self.deployment_identity,
            "task_identity": self.task_identity,
            "main_pid": self.main_pid,
            "process_start_ticks": self.process_start_ticks,
            "executable_device": self.executable_device,
            "executable_inode": self.executable_inode,
            "cmdline_identity": self.cmdline_identity,
            "boot_identity": self.boot_identity,
            "invocation_identity": self.invocation_identity,
            "cgroup_identity": self.cgroup_identity,
            "sockets": [item.to_dict() for item in self.sockets],
        }


@dataclass(frozen=True)
class OperatorRequest:
    operation: str
    task_identity: str
    deployment_identity: str
    component: str | None = None
    unit: str | None = None
    service_identity: str | None = None

    def __post_init__(self) -> None:
        component_operation = self.operation in {STAGE, ACTIVATE, ROLLBACK}
        service_operation = self.operation in SERVICE_OPERATIONS
        if (
            self.operation not in OPERATIONS
            or not _valid_task_identity(self.task_identity)
            or not _valid_identity(self.deployment_identity)
            or (component_operation and self.component not in COMPONENTS)
            or (not component_operation and self.component is not None)
            or (service_operation and self.unit not in SERVICE_UNITS)
            or (not service_operation and self.unit is not None)
            or (
                self.operation in BOUND_SERVICE_OPERATIONS
                and not _valid_identity(self.service_identity)
            )
            or (
                self.operation not in BOUND_SERVICE_OPERATIONS
                and self.service_identity is not None
            )
        ):
            raise fail("OPERATOR_REQUEST_INVALID", "Operator request is invalid.")


@dataclass(frozen=True)
class OperatorBoundaryObservation:
    status: str
    reason_code: str
    identity: str | None

    @classmethod
    def ready(cls, identity: str) -> "OperatorBoundaryObservation":
        return cls.from_dict(
            {"status": "ready", "reason_code": "READY", "identity": identity}
        )

    @classmethod
    def invalid(cls) -> "OperatorBoundaryObservation":
        return cls("not-ready", "BOUNDARY_INVALID", None)

    @classmethod
    def from_dict(cls, raw: object) -> "OperatorBoundaryObservation":
        if (
            not isinstance(raw, dict)
            or set(raw) != {"status", "reason_code", "identity"}
            or raw["status"] not in {"ready", "not-ready"}
            or raw["reason_code"] not in OPERATOR_BOUNDARY_REASONS
            or (raw["status"] == "ready") != _valid_identity(raw["identity"])
            or (raw["status"] == "ready") != (raw["reason_code"] == "READY")
        ):
            raise fail(
                "OPERATOR_OBSERVATION_INVALID", "Operator observation is invalid."
            )
        return cls(raw["status"], raw["reason_code"], raw["identity"])

    def to_dict(self) -> dict[str, object]:
        return {
            "status": self.status,
            "reason_code": self.reason_code,
            "identity": self.identity,
        }


class OperatorBoundaryObserver(Protocol):
    def observe(self) -> OperatorBoundaryObservation: ...


class HostOperatorBoundaryObserver:
    """Redact a full root-owned boundary verification to identity and reason."""

    def __init__(
        self,
        *,
        operator_root: Path = Path("/opt/helixweave/operator"),
        host_root: Path = Path("/"),
        root_uid: int = 0,
        root_gid: int = 0,
    ) -> None:
        self.operator_root = operator_root
        self.host_root = host_root
        self.root_uid = root_uid
        self.root_gid = root_gid

    def observe(self) -> OperatorBoundaryObservation:
        try:
            identity = verify_stable_operator_boundary(
                operator_root=self.operator_root,
                host_root=self.host_root,
                expected_uid=self.root_uid,
                expected_gid=self.root_gid,
            )
        except StableBoundaryError:
            return OperatorBoundaryObservation.invalid()
        return OperatorBoundaryObservation.ready(identity)


def _with_boundary_readiness(
    receipt: DeploymentActionReceipt,
    boundary: OperatorBoundaryObservation,
) -> DeploymentActionReceipt:
    """Bind a candidate receipt to the currently installed root boundary."""
    readiness = dict(receipt.readiness)
    candidate_permissions = readiness["permissions"]
    if boundary.status != "ready":
        readiness["permissions"] = ReadinessCheck("not-ready", "PERMISSION_INVALID")
    elif candidate_permissions.status != "ready":
        readiness["permissions"] = candidate_permissions
    else:
        readiness["permissions"] = ReadinessCheck(
            "ready",
            "READY",
            canonical_identity(
                {
                    "candidate_permissions_identity": candidate_permissions.identity,
                    "operator_boundary_identity": boundary.identity,
                },
                scheme=OPERATOR_PERMISSIONS_IDENTITY_SCHEME,
            ),
        )
    return DeploymentActionReceipt.create(
        request_identity=receipt.request_identity,
        status=receipt.status,
        compatibility=receipt.compatibility,
        database_before_identity=receipt.database_before_identity,
        database_after_identity=receipt.database_after_identity,
        accepted_schema_heads=receipt.accepted_schema_heads,
        target_schema_heads=receipt.target_schema_heads,
        migration_inventory_identity=receipt.migration_inventory_identity,
        known_schema_revisions=receipt.known_schema_revisions,
        migration_required=receipt.migration_required,
        rollback_supported=receipt.rollback_supported,
        api_contract_sha256=receipt.api_contract_sha256,
        native_identities=receipt.native_identities,
        frontend_identity=receipt.frontend_identity,
        reference_compatibility_identity=receipt.reference_compatibility_identity,
        readiness=readiness,
    )


@dataclass(frozen=True)
class OperatorObservation:
    identity: str
    state_identity: str
    active: dict[str, str | None]
    database_schema_identity: str | None
    database_schema_heads: tuple[str, ...]
    services: dict[str, str | None]
    operator_boundary: OperatorBoundaryObservation
    operator_pending_count: int
    operator_recovery_required_count: int

    @classmethod
    def create(
        cls,
        *,
        state_identity: str,
        active: dict[str, str | None],
        database_schema_identity: str | None,
        database_schema_heads: Sequence[str],
        services: dict[str, str | None],
        operator_boundary: OperatorBoundaryObservation | None = None,
        operator_pending_count: int = 0,
        operator_recovery_required_count: int = 0,
    ) -> "OperatorObservation":
        selected_boundary = (
            OperatorBoundaryObservation.invalid()
            if operator_boundary is None
            else operator_boundary
        )
        value: dict[str, object] = {
            "schema_version": OPERATOR_OBSERVATION_SCHEMA,
            "state_identity": state_identity,
            "active": active,
            "database_schema_identity": database_schema_identity,
            "database_schema_heads": list(database_schema_heads),
            "services": services,
            "operator_boundary": selected_boundary.to_dict(),
            "operator_pending_count": operator_pending_count,
            "operator_recovery_required_count": operator_recovery_required_count,
        }
        value["identity"] = canonical_identity(
            value, scheme=OPERATOR_OBSERVATION_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "OperatorObservation":
        if not isinstance(raw, dict) or set(raw) != {
            "schema_version",
            "identity",
            "state_identity",
            "active",
            "database_schema_identity",
            "database_schema_heads",
            "services",
            "operator_boundary",
            "operator_pending_count",
            "operator_recovery_required_count",
        }:
            raise fail(
                "OPERATOR_OBSERVATION_INVALID", "Operator observation is invalid."
            )
        active = raw["active"]
        services = raw["services"]
        heads = raw["database_schema_heads"]
        boundary = OperatorBoundaryObservation.from_dict(raw["operator_boundary"])
        if (
            raw["schema_version"] != OPERATOR_OBSERVATION_SCHEMA
            or not isinstance(active, dict)
            or set(active) != set(COMPONENTS)
            or any(
                value is not None and not _valid_identity(value)
                for value in active.values()
            )
            or not isinstance(services, dict)
            or set(services) != set(SERVICE_UNITS)
            or any(
                value is not None and not _valid_identity(value)
                for value in services.values()
            )
            or not isinstance(heads, list)
            or tuple(sorted(set(heads))) != tuple(heads)
            or any(
                not isinstance(head, str)
                or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", head) is None
                for head in heads
            )
            or not _valid_identity(raw["state_identity"])
            or (
                raw["database_schema_identity"] is not None
                and not _valid_identity(raw["database_schema_identity"])
            )
            or (raw["database_schema_identity"] is None and bool(heads))
            or not _valid_identity(raw["identity"])
            or not isinstance(raw["operator_pending_count"], int)
            or isinstance(raw["operator_pending_count"], bool)
            or raw["operator_pending_count"] not in {0, 1}
            or not isinstance(raw["operator_recovery_required_count"], int)
            or isinstance(raw["operator_recovery_required_count"], bool)
            or raw["operator_recovery_required_count"] not in {0, 1}
            or raw["operator_pending_count"] + raw["operator_recovery_required_count"]
            > 1
        ):
            raise fail(
                "OPERATOR_OBSERVATION_INVALID", "Operator observation is invalid."
            )
        expected = canonical_identity(
            {key: value for key, value in raw.items() if key != "identity"},
            scheme=OPERATOR_OBSERVATION_IDENTITY_SCHEME,
        )
        if raw["identity"] != expected:
            raise fail(
                "OPERATOR_OBSERVATION_INVALID", "Operator observation is invalid."
            )
        return cls(
            identity=expected,
            state_identity=raw["state_identity"],
            active={component: active[component] for component in COMPONENTS},
            database_schema_identity=raw["database_schema_identity"],
            database_schema_heads=tuple(heads),
            services={unit: services[unit] for unit in SERVICE_UNITS},
            operator_boundary=boundary,
            operator_pending_count=raw["operator_pending_count"],
            operator_recovery_required_count=raw["operator_recovery_required_count"],
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": OPERATOR_OBSERVATION_SCHEMA,
            "identity": self.identity,
            "state_identity": self.state_identity,
            "active": self.active,
            "database_schema_identity": self.database_schema_identity,
            "database_schema_heads": list(self.database_schema_heads),
            "services": self.services,
            "operator_boundary": self.operator_boundary.to_dict(),
            "operator_pending_count": self.operator_pending_count,
            "operator_recovery_required_count": self.operator_recovery_required_count,
        }

    def with_operator_boundary(
        self, boundary: OperatorBoundaryObservation
    ) -> "OperatorObservation":
        return OperatorObservation.create(
            state_identity=self.state_identity,
            active=self.active,
            database_schema_identity=self.database_schema_identity,
            database_schema_heads=self.database_schema_heads,
            services=self.services,
            operator_boundary=boundary,
            operator_pending_count=self.operator_pending_count,
            operator_recovery_required_count=self.operator_recovery_required_count,
        )

    def with_operator_journal(
        self,
        summary: OperatorJournalSummary,
    ) -> "OperatorObservation":
        return OperatorObservation.create(
            state_identity=self.state_identity,
            active=self.active,
            database_schema_identity=self.database_schema_identity,
            database_schema_heads=self.database_schema_heads,
            services=self.services,
            operator_boundary=self.operator_boundary,
            operator_pending_count=summary.pending_count,
            operator_recovery_required_count=summary.recovery_required_count,
        )


@dataclass(frozen=True)
class OperatorOutcome:
    """Whitelisted backend result; it cannot carry paths or exception strings."""

    state: str
    service: ServiceIdentity | None = None
    observation: OperatorObservation | None = None
    verification: DeploymentActionReceipt | None = None

    def __post_init__(self) -> None:
        if (
            self.state
            not in {
                "staged",
                "activated",
                "rolled-back",
                "running",
                "stopped",
                "clean",
                "uninstalled",
                "observed",
                "verified",
            }
            or (
                self.service is not None
                and not isinstance(self.service, ServiceIdentity)
            )
            or (
                self.observation is not None
                and not isinstance(self.observation, OperatorObservation)
            )
            or (
                self.verification is not None
                and not isinstance(self.verification, DeploymentActionReceipt)
            )
        ):
            raise fail("OPERATOR_RESULT_INVALID", "Operator result is invalid.")


class OperatorBackend(Protocol):
    """Privileged implementation seam; callers can never supply a command."""

    def execute(
        self,
        request: OperatorRequest,
        *,
        bundle_path: Path | None,
    ) -> OperatorOutcome: ...


class UnavailableOperatorBackend:
    """Phase-A default: fail closed until post-PR-172 host wiring is installed."""

    def execute(
        self,
        request: OperatorRequest,
        *,
        bundle_path: Path | None,
    ) -> OperatorOutcome:
        del request, bundle_path
        raise fail(
            "OPERATOR_BACKEND_UNAVAILABLE",
            "Operator backend is not installed.",
            recoverable=True,
        )


@dataclass(frozen=True)
class CommandResult:
    returncode: int
    stdout: bytes = b""


class CommandExecutor(Protocol):
    def run(self, argv: tuple[str, ...], *, timeout: float) -> CommandResult: ...


class BoundedCommandExecutor:
    """Execute fixed absolute argv with no stdin or caller environment."""

    def run(self, argv: tuple[str, ...], *, timeout: float) -> CommandResult:
        if (
            not argv
            or argv[0] != str(SYSTEMCTL)
            or not 0 < timeout <= _SYSTEMCTL_TIMEOUT_SECONDS
        ):
            raise fail("OPERATOR_COMMAND_INVALID", "Operator command is invalid.")
        try:
            with tempfile.TemporaryFile() as stdout, tempfile.TemporaryFile() as stderr:
                completed = subprocess.run(
                    argv,
                    stdin=subprocess.DEVNULL,
                    stdout=stdout,
                    stderr=stderr,
                    env=SAFE_ENVIRONMENT,
                    close_fds=True,
                    timeout=timeout,
                    check=False,
                )
                if stdout.tell() > _MAX_COMMAND_OUTPUT or stderr.tell() > (
                    _MAX_COMMAND_OUTPUT
                ):
                    raise fail("OPERATOR_COMMAND_FAILED", "Operator command failed.")
                stdout.seek(0)
                content = stdout.read(_MAX_COMMAND_OUTPUT + 1)
        except DeploymentError:
            raise
        except (OSError, subprocess.SubprocessError):
            raise fail(
                "OPERATOR_COMMAND_FAILED",
                "Operator command failed.",
                recoverable=True,
            ) from None
        return CommandResult(completed.returncode, content)


class BulkRuntimePreparer(Protocol):
    def observe_boundary(self) -> BulkDockerBoundary: ...

    def prepare(
        self, request: BulkRuntimePrepareRequest
    ) -> BulkRuntimePrepareReceipt: ...


class SystemdBulkRuntimePreparer:
    """Exchange one canonical bulk request through the fixed oneshot unit."""

    def __init__(
        self,
        layout: DeploymentLayout,
        *,
        service_uid: int,
        service_gid: int,
        root_uid: int = 0,
        root_gid: int = 0,
        command_runner: Callable[[tuple[str, ...]], int] | None = None,
        boundary_observer: Callable[[int, int], BulkDockerBoundary] | None = None,
    ) -> None:
        if any(
            not isinstance(value, int) or isinstance(value, bool) or value < 0
            for value in (service_uid, service_gid, root_uid, root_gid)
        ):
            raise fail(
                "BULK_RUNTIME_PREPARE_BOUNDARY_INVALID",
                "Bulk runtime preparation is unavailable.",
            )
        self.layout = layout
        self.service_uid = service_uid
        self.service_gid = service_gid
        self.root_uid = root_uid
        self.root_gid = root_gid
        self.command_runner = (
            self._run_systemctl if command_runner is None else command_runner
        )
        self.boundary_observer = (
            observe_bulk_docker_boundary
            if boundary_observer is None
            else boundary_observer
        )

    def observe_boundary(self) -> BulkDockerBoundary:
        try:
            boundary = self.boundary_observer(self.service_uid, self.service_gid)
        except DeploymentError:
            raise
        except Exception:
            raise fail(
                "BULK_RUNTIME_DOCKER_BOUNDARY_INVALID",
                "Bulk runtime Docker boundary is invalid.",
            ) from None
        if (
            not isinstance(boundary, BulkDockerBoundary)
            or boundary.daemon_uid != self.service_uid
            or boundary.daemon_gid != self.service_gid
            or not _valid_identity(boundary.client_identity)
            or not _valid_identity(boundary.endpoint_identity)
        ):
            raise fail(
                "BULK_RUNTIME_DOCKER_BOUNDARY_INVALID",
                "Bulk runtime Docker boundary is invalid.",
            )
        return boundary

    def prepare(self, request: BulkRuntimePrepareRequest) -> BulkRuntimePrepareReceipt:
        if not isinstance(request, BulkRuntimePrepareRequest):
            raise fail(
                "BULK_RUNTIME_PREPARE_REQUEST_INVALID",
                "Bulk runtime preparation request is invalid.",
            )
        root = self._require_root()
        lock = self._lock(root)
        try:
            self._replace_request(root, request)
            self._prepare_receipt(root)
            expected = (
                str(SYSTEMCTL),
                "--no-ask-password",
                "start",
                "--",
                _BULK_RUNTIME_PREPARE_UNIT,
            )
            try:
                return_code = self.command_runner(expected)
            except Exception:
                raise fail(
                    "BULK_RUNTIME_PREPARE_FAILED",
                    "Bulk runtime preparation failed.",
                    recoverable=True,
                ) from None
            if return_code != 0:
                raise fail(
                    "BULK_RUNTIME_PREPARE_FAILED",
                    "Bulk runtime preparation failed.",
                    recoverable=True,
                )
            receipt = self._read_receipt()
            if (
                receipt.request_identity != request.identity
                or receipt.candidate_bulk_identity != request.candidate_bulk_identity
            ):
                raise fail(
                    "BULK_RUNTIME_PREPARE_RECEIPT_INVALID",
                    "Bulk runtime preparation receipt is invalid.",
                )
            return receipt
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)
            os.close(lock)

    def _require_root(self) -> Path:
        root = self.layout.bulk_runtime_prepare_root
        try:
            observed = root.lstat()
        except OSError:
            raise fail(
                "BULK_RUNTIME_PREPARE_BOUNDARY_INVALID",
                "Bulk runtime preparation is unavailable.",
            ) from None
        if (
            not stat.S_ISDIR(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or observed.st_uid != self.root_uid
            or observed.st_gid != self.service_gid
            or stat.S_IMODE(observed.st_mode) != 0o750
        ):
            raise fail(
                "BULK_RUNTIME_PREPARE_BOUNDARY_INVALID",
                "Bulk runtime preparation is unavailable.",
            )
        return root

    def _lock(self, root: Path) -> int:
        descriptor = -1
        try:
            descriptor = os.open(
                root / "operation.lock",
                os.O_RDWR
                | os.O_CREAT
                | getattr(os, "O_CLOEXEC", 0)
                | getattr(os, "O_NOFOLLOW", 0),
                0o600,
            )
            os.fchown(descriptor, self.root_uid, self.root_gid)
            os.fchmod(descriptor, 0o600)
            observed = os.fstat(descriptor)
            if (
                not stat.S_ISREG(observed.st_mode)
                or observed.st_nlink != 1
                or observed.st_uid != self.root_uid
                or observed.st_gid != self.root_gid
                or stat.S_IMODE(observed.st_mode) != 0o600
            ):
                raise OSError
            fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            return descriptor
        except BlockingIOError:
            if descriptor >= 0:
                os.close(descriptor)
            raise fail(
                "BULK_RUNTIME_PREPARE_BUSY",
                "Bulk runtime preparation is already in progress.",
                recoverable=True,
            ) from None
        except OSError:
            if descriptor >= 0:
                os.close(descriptor)
            raise fail(
                "BULK_RUNTIME_PREPARE_BOUNDARY_INVALID",
                "Bulk runtime preparation is unavailable.",
            ) from None

    def _replace_request(self, root: Path, request: BulkRuntimePrepareRequest) -> None:
        temporary = root / f".{request.task_identity}.request"
        self._write_file(
            temporary,
            canonical_json_bytes(request.to_dict()),
            uid=self.root_uid,
            gid=self.service_gid,
            mode=0o640,
        )
        try:
            os.replace(temporary, self.layout.bulk_runtime_prepare_request)
            _fsync_directory(root)
        except OSError:
            try:
                temporary.unlink()
            except OSError:
                pass
            raise fail(
                "BULK_RUNTIME_PREPARE_BOUNDARY_INVALID",
                "Bulk runtime preparation is unavailable.",
            ) from None

    def _prepare_receipt(self, root: Path) -> None:
        receipt = self.layout.bulk_runtime_prepare_receipt
        try:
            observed = receipt.lstat()
        except FileNotFoundError:
            pass
        except OSError:
            raise fail(
                "BULK_RUNTIME_PREPARE_BOUNDARY_INVALID",
                "Bulk runtime preparation is unavailable.",
            ) from None
        else:
            if (
                not stat.S_ISREG(observed.st_mode)
                or stat.S_ISLNK(observed.st_mode)
                or observed.st_nlink != 1
                or observed.st_uid != self.service_uid
                or observed.st_gid != self.service_gid
                or stat.S_IMODE(observed.st_mode) != 0o600
            ):
                raise fail(
                    "BULK_RUNTIME_PREPARE_BOUNDARY_INVALID",
                    "Bulk runtime preparation is unavailable.",
                )
            try:
                receipt.unlink()
            except OSError:
                raise fail(
                    "BULK_RUNTIME_PREPARE_BOUNDARY_INVALID",
                    "Bulk runtime preparation is unavailable.",
                ) from None
        self._write_file(
            receipt,
            b"",
            uid=self.service_uid,
            gid=self.service_gid,
            mode=0o600,
        )
        _fsync_directory(root)

    def _read_receipt(self) -> BulkRuntimePrepareReceipt:
        path = self.layout.bulk_runtime_prepare_receipt
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
                or not 0 < before.st_size <= _MAX_OPERATOR_DOCUMENT_BYTES
            ):
                raise OSError
            content = _read_bounded_fd(descriptor, _MAX_OPERATOR_DOCUMENT_BYTES)
            after = os.fstat(descriptor)
            at_path = path.stat(follow_symlinks=False)
            if (
                len(content) != before.st_size
                or _file_witness(before) != _file_witness(after)
                or (after.st_dev, after.st_ino) != (at_path.st_dev, at_path.st_ino)
            ):
                raise OSError
            raw = json.loads(content, object_pairs_hook=_unique_object)
            receipt = BulkRuntimePrepareReceipt.from_dict(raw)
            if canonical_json_bytes(receipt.to_dict()) != content:
                raise ValueError
            return receipt
        except (
            DeploymentError,
            OSError,
            UnicodeError,
            ValueError,
            json.JSONDecodeError,
        ):
            raise fail(
                "BULK_RUNTIME_PREPARE_RECEIPT_INVALID",
                "Bulk runtime preparation receipt is invalid.",
            ) from None
        finally:
            if descriptor >= 0:
                os.close(descriptor)

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
                "BULK_RUNTIME_PREPARE_BOUNDARY_INVALID",
                "Bulk runtime preparation is unavailable.",
            ) from None
        finally:
            if descriptor >= 0:
                os.close(descriptor)

    @staticmethod
    def _run_systemctl(argv: tuple[str, ...]) -> int:
        expected = (
            str(SYSTEMCTL),
            "--no-ask-password",
            "start",
            "--",
            _BULK_RUNTIME_PREPARE_UNIT,
        )
        if argv != expected:
            raise fail(
                "BULK_RUNTIME_PREPARE_COMMAND_INVALID",
                "Bulk runtime preparation is unavailable.",
            )
        try:
            with tempfile.TemporaryFile() as stdout, tempfile.TemporaryFile() as stderr:
                completed = subprocess.run(
                    argv,
                    stdin=subprocess.DEVNULL,
                    stdout=stdout,
                    stderr=stderr,
                    env=SAFE_ENVIRONMENT,
                    close_fds=True,
                    timeout=_BULK_RUNTIME_SYSTEMD_TIMEOUT_SECONDS,
                    check=False,
                )
                if stdout.tell() or stderr.tell():
                    return 1
                return completed.returncode
        except (OSError, subprocess.SubprocessError):
            return 1


class FixedSystemctl:
    """The only production command surface exposed by the root backend."""

    _ACTIONS = frozenset({"start", "stop", "reset-failed"})
    _SHOW_PROPERTIES = (
        "ActiveState",
        "SubState",
        "MainPID",
        "InvocationID",
        "ControlGroup",
        "NeedDaemonReload",
    )

    def __init__(self, executor: CommandExecutor | None = None) -> None:
        self.executor = BoundedCommandExecutor() if executor is None else executor

    def control(self, action: str, unit: str) -> None:
        if action not in self._ACTIONS or unit not in SERVICE_UNITS:
            raise fail("OPERATOR_COMMAND_INVALID", "Operator command is invalid.")
        result = self.executor.run(
            (str(SYSTEMCTL), "--no-ask-password", action, "--", unit),
            timeout=_SYSTEMCTL_TIMEOUT_SECONDS,
        )
        if result.returncode != 0 or result.stdout:
            raise fail(
                "OPERATOR_SERVICE_CONTROL_FAILED",
                "Service action failed.",
                recoverable=True,
            )

    def show(self, unit: str) -> dict[str, str]:
        if unit not in SERVICE_UNITS:
            raise fail("OPERATOR_COMMAND_INVALID", "Operator command is invalid.")
        properties = ",".join(self._SHOW_PROPERTIES)
        result = self.executor.run(
            (
                str(SYSTEMCTL),
                "--no-pager",
                "--no-legend",
                "show",
                f"--property={properties}",
                "--",
                unit,
            ),
            timeout=_SYSTEMCTL_TIMEOUT_SECONDS,
        )
        if result.returncode != 0:
            raise fail(
                "OPERATOR_SERVICE_OBSERVE_FAILED",
                "Service status could not be observed.",
                recoverable=True,
            )
        try:
            text = result.stdout.decode("ascii")
        except UnicodeDecodeError:
            raise fail(
                "OPERATOR_SERVICE_OBSERVE_FAILED",
                "Service status could not be observed.",
            ) from None
        values: dict[str, str] = {}
        for line in text.splitlines():
            key, separator, value = line.partition("=")
            if (
                separator != "="
                or key not in self._SHOW_PROPERTIES
                or key in values
                or len(value) > 4096
                or any(character in value for character in ("\x00", "\r", "\n"))
            ):
                raise fail(
                    "OPERATOR_SERVICE_OBSERVE_FAILED",
                    "Service status could not be observed.",
                )
            values[key] = value
        if set(values) != set(self._SHOW_PROPERTIES):
            raise fail(
                "OPERATOR_SERVICE_OBSERVE_FAILED",
                "Service status could not be observed.",
            )
        if values["NeedDaemonReload"] not in {"yes", "no"}:
            raise fail(
                "OPERATOR_SERVICE_OBSERVE_FAILED",
                "Service status could not be observed.",
            )
        if values["NeedDaemonReload"] == "yes":
            raise fail(
                "OPERATOR_SYSTEMD_RELOAD_REQUIRED",
                "System service definitions require a daemon reload.",
                recoverable=True,
            )
        return values


class ServiceProbe(Protocol):
    def observe(
        self,
        *,
        unit: str,
        deployment_identity: str,
        task_identity: str,
    ) -> ServiceIdentity | None: ...


class LinuxServiceProbe:
    """Bind systemd, procfs, cgroup, executable, command line, and sockets."""

    _MAX_CGROUP_PROCS_BYTES = 64 * 1024
    _MAX_CGROUP_PROCS = 4096
    _MAX_PROC_FDS = 65_536
    _MAX_PROC_NET_UNIX_BYTES = 8 * 1024 * 1024
    _MAX_PROC_NET_UNIX_LINES = 131_072
    _SO_ACCEPTCON = 0x00010000

    _UNIX_SOCKETS = {
        "helixweave-redis.service": Path("/run/helixweave/redis/redis.sock"),
        "helixweave-docker-rootless.service": Path(
            "/run/helixweave/docker/docker.sock"
        ),
    }

    def __init__(
        self,
        systemctl: FixedSystemctl,
        *,
        proc_root: Path = Path("/proc"),
        cgroup_root: Path = Path("/sys/fs/cgroup"),
        unix_sockets: dict[str, Path] | None = None,
        filesystem_socket_stat: Callable[[Path], os.stat_result] | None = None,
    ) -> None:
        self.systemctl = systemctl
        self.proc_root = proc_root
        self.cgroup_root = cgroup_root
        self.unix_sockets = (
            dict(self._UNIX_SOCKETS) if unix_sockets is None else dict(unix_sockets)
        )
        self.filesystem_socket_stat = (
            (lambda path: path.stat(follow_symlinks=False))
            if filesystem_socket_stat is None
            else filesystem_socket_stat
        )
        if set(self.unix_sockets) != set(self._UNIX_SOCKETS) or any(
            not isinstance(path, Path) or not path.is_absolute()
            for path in self.unix_sockets.values()
        ):
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_INVALID",
                "Service identity evidence is invalid.",
            )

    def observe(
        self,
        *,
        unit: str,
        deployment_identity: str,
        task_identity: str,
    ) -> ServiceIdentity | None:
        values = self.systemctl.show(unit)
        active = values["ActiveState"]
        try:
            main_pid = int(values["MainPID"])
        except ValueError:
            raise fail(
                "OPERATOR_SERVICE_OBSERVE_FAILED",
                "Service status could not be observed.",
            ) from None
        if active in {"inactive", "failed", "deactivating"}:
            if main_pid != 0:
                raise fail(
                    "OPERATOR_SERVICE_OBSERVE_FAILED",
                    "Service status could not be observed.",
                )
            return None
        expected_cgroup = f"/system.slice/{unit}"
        if (
            active != "active"
            or values["SubState"] not in {"running", "listening"}
            or main_pid <= 0
            or values["ControlGroup"] != expected_cgroup
            or re.fullmatch(r"[0-9a-f]{32}", values["InvocationID"]) is None
        ):
            raise fail(
                "OPERATOR_SERVICE_OBSERVE_FAILED",
                "Service status could not be observed.",
            )
        try:
            process = self.proc_root / str(main_pid)
            raw_stat = (process / "stat").read_bytes()
            closing = raw_stat.rfind(b")")
            fields = raw_stat[closing + 2 :].split()
            start_ticks = int(fields[19])
            executable = (process / "exe").stat()
            cmdline = (process / "cmdline").read_bytes()
            boot = (self.proc_root / "sys/kernel/random/boot_id").read_bytes()
        except (OSError, ValueError, IndexError):
            raise fail(
                "OPERATOR_SERVICE_OBSERVE_FAILED",
                "Service status could not be observed.",
            ) from None
        if (
            closing < 1
            or not stat.S_ISREG(executable.st_mode)
            or executable.st_dev <= 0
            or executable.st_ino <= 0
            or not 0 < len(cmdline) <= _MAX_COMMAND_OUTPUT
            or not 0 < len(boot) <= 128
        ):
            raise fail(
                "OPERATOR_SERVICE_OBSERVE_FAILED",
                "Service status could not be observed.",
            )
        sockets = self._socket_witnesses(
            unit=unit,
            cgroup=expected_cgroup,
            main_pid=main_pid,
        )
        try:
            final_values = self.systemctl.show(unit)
            final_stat = (process / "stat").read_bytes()
            final_closing = final_stat.rfind(b")")
            final_fields = final_stat[final_closing + 2 :].split()
            final_start_ticks = int(final_fields[19])
            final_executable = (process / "exe").stat()
            final_cmdline = (process / "cmdline").read_bytes()
        except (OSError, ValueError, IndexError):
            raise fail(
                "OPERATOR_SERVICE_OBSERVE_FAILED",
                "Service status could not be observed.",
            ) from None
        if (
            final_closing < 1
            or final_values != values
            or final_start_ticks != start_ticks
            or _file_witness(final_executable) != _file_witness(executable)
            or final_cmdline != cmdline
        ):
            raise fail(
                "OPERATOR_SERVICE_OBSERVE_FAILED",
                "Service status could not be observed.",
            )
        return ServiceIdentity.create(
            unit=unit,
            deployment_identity=deployment_identity,
            task_identity=task_identity,
            main_pid=main_pid,
            process_start_ticks=start_ticks,
            executable_device=executable.st_dev,
            executable_inode=executable.st_ino,
            cmdline_identity=_bytes_identity(cmdline),
            boot_identity=_bytes_identity(boot.strip()),
            invocation_identity=_bytes_identity(values["InvocationID"].encode()),
            cgroup_identity=_bytes_identity(expected_cgroup.encode()),
            sockets=sockets,
        )

    def _socket_witnesses(
        self,
        *,
        unit: str,
        cgroup: str,
        main_pid: int,
    ) -> tuple[SocketWitness, ...]:
        names = SERVICE_SOCKET_NAMES[unit]
        if not names:
            return ()
        if unit == "helixweave-api.service":
            observed = self._api_socket(cgroup, main_pid=main_pid)
            return (
                SocketWitness(
                    "api-http",
                    observed.st_dev,
                    observed.st_ino,
                    observed.st_ino,
                ),
            )
        path = self.unix_sockets[unit]
        try:
            before = self.filesystem_socket_stat(path)
            pids_before = self._cgroup_pids(cgroup)
            if main_pid not in pids_before:
                raise OSError
            kernel_inode = self._listening_unix_socket_inode(path)
            if not self._pids_own_socket(pids_before, kernel_inode):
                raise OSError
            pids_after = self._cgroup_pids(cgroup)
            after = self.filesystem_socket_stat(path)
        except OSError:
            raise fail(
                "OPERATOR_SERVICE_OBSERVE_FAILED",
                "Service socket could not be observed.",
            ) from None
        if (
            not stat.S_ISSOCK(before.st_mode)
            or before.st_nlink != 1
            or _file_witness(before) != _file_witness(after)
            or pids_before != pids_after
        ):
            raise fail(
                "OPERATOR_SERVICE_OBSERVE_FAILED",
                "Service socket could not be observed.",
            )
        return (
            SocketWitness(
                names[0],
                after.st_dev,
                after.st_ino,
                kernel_inode,
            ),
        )

    def _api_socket(self, cgroup: str, *, main_pid: int) -> os.stat_result:
        inodes: set[int] = set()
        try:
            lines = (
                (self.proc_root / "net/tcp").read_text(encoding="ascii").splitlines()
            )
        except OSError:
            raise fail(
                "OPERATOR_SERVICE_OBSERVE_FAILED",
                "Service socket could not be observed.",
            ) from None
        for line in lines[1:]:
            fields = line.split()
            if len(fields) >= 10 and fields[1] == "0100007F:1F40" and fields[3] == "0A":
                try:
                    inodes.add(int(fields[9]))
                except ValueError:
                    pass
        if len(inodes) != 1:
            raise fail(
                "OPERATOR_SERVICE_OBSERVE_FAILED",
                "Service socket could not be observed.",
            )
        inode = next(iter(inodes))
        observed = self._cgroup_socket_stat(cgroup, inode, main_pid=main_pid)
        if observed is None:
            raise fail(
                "OPERATOR_SERVICE_OBSERVE_FAILED",
                "Service socket could not be observed.",
            )
        return observed

    def _listening_unix_socket_inode(self, socket_path: Path) -> int:
        try:
            content = _read_bounded_path(
                self.proc_root / "net/unix",
                self._MAX_PROC_NET_UNIX_BYTES,
            )
            lines = content.splitlines()
            if not lines or len(lines) > self._MAX_PROC_NET_UNIX_LINES:
                raise OSError
            rendered_path = str(socket_path)
            matches: list[int] = []
            for raw_line in lines[1:]:
                fields = raw_line.decode("utf-8").split(maxsplit=7)
                if len(fields) != 8 or fields[7] != rendered_path:
                    continue
                flags = int(fields[3], 16)
                kernel_inode = int(fields[6], 10)
                if (
                    fields[4] != "0001"
                    or fields[5] != "01"
                    or not flags & self._SO_ACCEPTCON
                    or not 0 < kernel_inode <= 2**63 - 1
                ):
                    raise OSError
                matches.append(kernel_inode)
            if len(matches) != 1:
                raise OSError
            return matches[0]
        except (OSError, UnicodeError, ValueError):
            raise OSError from None

    def _cgroup_pids(self, cgroup: str) -> frozenset[int]:
        content = _read_bounded_path(
            self.cgroup_root / cgroup.lstrip("/") / "cgroup.procs",
            self._MAX_CGROUP_PROCS_BYTES,
        )
        try:
            values = content.decode("ascii").split()
            pids = frozenset(int(value) for value in values)
        except (UnicodeDecodeError, ValueError):
            raise OSError from None
        if (
            not pids
            or len(values) != len(pids)
            or len(pids) > self._MAX_CGROUP_PROCS
            or any(pid <= 0 or pid > 2**31 - 1 for pid in pids)
        ):
            raise OSError
        return pids

    def _pids_own_socket(self, pids: frozenset[int], inode: int) -> bool:
        owned = False
        for pid in pids:
            directory = self.proc_root / str(pid) / "fd"
            try:
                count = 0
                with os.scandir(directory) as descriptors:
                    for descriptor in descriptors:
                        count += 1
                        if count > self._MAX_PROC_FDS or not descriptor.name.isdigit():
                            raise OSError
                        try:
                            target = os.readlink(descriptor.path)
                        except OSError:
                            raise OSError from None
                        if target == f"socket:[{inode}]":
                            owned = True
            except OSError:
                raise OSError from None
        return owned

    def _cgroup_socket_stat(
        self,
        cgroup: str,
        inode: int,
        *,
        main_pid: int,
    ) -> os.stat_result | None:
        try:
            before = self._cgroup_pids(cgroup)
            if main_pid not in before:
                return None
            observed: os.stat_result | None = None
            for pid in before:
                directory = self.proc_root / str(pid) / "fd"
                count = 0
                with os.scandir(directory) as descriptors:
                    for descriptor in descriptors:
                        count += 1
                        if count > self._MAX_PROC_FDS or not descriptor.name.isdigit():
                            raise OSError
                        if os.readlink(descriptor.path) != f"socket:[{inode}]":
                            continue
                        candidate = descriptor.stat(follow_symlinks=True)
                        if (
                            not stat.S_ISSOCK(candidate.st_mode)
                            or candidate.st_ino != inode
                        ):
                            raise OSError
                        observed = candidate
            after = self._cgroup_pids(cgroup)
            if before != after:
                return None
            return observed
        except OSError:
            return None


class ServiceController(Protocol):
    def start(self, request: OperatorRequest) -> ServiceIdentity: ...

    def recover_start(self, request: OperatorRequest) -> ServiceIdentity: ...

    def recover_observe(self, request: OperatorRequest) -> ServiceIdentity | None: ...

    def status(self, request: OperatorRequest) -> ServiceIdentity | None: ...

    def stop(self, request: OperatorRequest, *, cleanup: bool) -> None: ...

    def recover_stop(self, request: OperatorRequest, *, cleanup: bool) -> None: ...


class SystemdServiceController:
    """Persist start identity and re-observe it before every destructive action."""

    def __init__(
        self,
        layout: DeploymentLayout,
        *,
        systemctl: FixedSystemctl | None = None,
        probe: ServiceProbe | None = None,
        owner_uid: int = 0,
        owner_gid: int = 0,
    ) -> None:
        self.layout = layout
        self.systemctl = FixedSystemctl() if systemctl is None else systemctl
        self.probe = LinuxServiceProbe(self.systemctl) if probe is None else probe
        self.owner_uid = owner_uid
        self.owner_gid = owner_gid

    def start(self, request: OperatorRequest) -> ServiceIdentity:
        assert request.unit is not None
        existing = self.probe.observe(
            unit=request.unit,
            deployment_identity=request.deployment_identity,
            task_identity=request.task_identity,
        )
        if existing is not None:
            raise fail(
                "OPERATOR_SERVICE_ALREADY_RUNNING", "Service is already running."
            )
        self.systemctl.control("start", request.unit)
        service = self.probe.observe(
            unit=request.unit,
            deployment_identity=request.deployment_identity,
            task_identity=request.task_identity,
        )
        if service is None:
            raise fail(
                "OPERATOR_SERVICE_START_FAILED",
                "Service did not enter the running state.",
                recoverable=True,
            )
        self._persist_started_service(service)
        return service

    def recover_start(self, request: OperatorRequest) -> ServiceIdentity:
        """Adopt only this journal-bound start, or execute it idempotently."""
        service = self.recover_observe(request)
        if service is None:
            return self.start(request)
        return service

    def recover_observe(self, request: OperatorRequest) -> ServiceIdentity | None:
        """Adopt a running service only under an authoritative recovery journal."""
        assert request.unit is not None
        service = self.probe.observe(
            unit=request.unit,
            deployment_identity=request.deployment_identity,
            task_identity=request.task_identity,
        )
        if service is None:
            return None
        self._persist_started_service(service)
        return service

    def status(self, request: OperatorRequest) -> ServiceIdentity | None:
        assert request.unit is not None
        prior = self._read_identity(request.unit, required=False)
        if (
            prior is not None
            and prior.deployment_identity != request.deployment_identity
        ):
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_MISMATCH",
                "Service identity does not match this deployment.",
            )
        task = request.task_identity if prior is None else prior.task_identity
        service = self.probe.observe(
            unit=request.unit,
            deployment_identity=request.deployment_identity,
            task_identity=task,
        )
        if service is None:
            return None
        if prior is None:
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_UNAVAILABLE",
                "Running service has no trusted operator identity.",
            )
        if service.identity != prior.identity:
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_MISMATCH",
                "Service identity changed before the requested action.",
                recoverable=True,
            )
        return service

    def stop(self, request: OperatorRequest, *, cleanup: bool) -> None:
        assert request.unit is not None and request.service_identity is not None
        prior = self._read_identity(request.unit, required=True)
        if prior is None or prior.deployment_identity != request.deployment_identity:
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_MISMATCH",
                "Service identity does not match this deployment.",
            )
        current = self.probe.observe(
            unit=request.unit,
            deployment_identity=prior.deployment_identity,
            task_identity=prior.task_identity,
        )
        if current is None or current.identity != request.service_identity:
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_MISMATCH",
                "Service identity changed before the requested action.",
                recoverable=True,
            )
        self.systemctl.control("stop", request.unit)
        stopped = self.probe.observe(
            unit=request.unit,
            deployment_identity=prior.deployment_identity,
            task_identity=prior.task_identity,
        )
        if stopped is not None:
            raise fail(
                "OPERATOR_SERVICE_STOP_FAILED",
                "Service did not enter the stopped state.",
                recoverable=True,
            )
        if cleanup:
            self.systemctl.control("reset-failed", request.unit)

    def recover_stop(self, request: OperatorRequest, *, cleanup: bool) -> None:
        """Finish a journal-bound stop even when systemd already stopped it."""
        assert request.unit is not None and request.service_identity is not None
        prior = self._read_identity(request.unit, required=True)
        if prior is None or prior.deployment_identity != request.deployment_identity:
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_MISMATCH",
                "Service identity does not match this deployment.",
            )
        current = self.probe.observe(
            unit=request.unit,
            deployment_identity=prior.deployment_identity,
            task_identity=prior.task_identity,
        )
        if current is None:
            if cleanup:
                self.systemctl.control("reset-failed", request.unit)
            return
        if current.identity != request.service_identity:
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_MISMATCH",
                "Service identity changed before the requested action.",
                recoverable=True,
            )
        self.stop(request, cleanup=cleanup)

    def _persist_started_service(self, service: ServiceIdentity) -> None:
        try:
            self._write_identity(service)
        except DeploymentError:
            try:
                self.systemctl.control("stop", service.unit)
                stopped = self.probe.observe(
                    unit=service.unit,
                    deployment_identity=service.deployment_identity,
                    task_identity=service.task_identity,
                )
            except DeploymentError:
                raise fail(
                    "OPERATOR_SERVICE_RECOVERY_REQUIRED",
                    "Service start requires recovery.",
                    recoverable=True,
                ) from None
            if stopped is not None:
                raise fail(
                    "OPERATOR_SERVICE_RECOVERY_REQUIRED",
                    "Service start requires recovery.",
                    recoverable=True,
                )
            raise

    def _identity_path(self, unit: str) -> Path:
        return self.layout.service_identities / f"{unit}.json"

    def _write_identity(self, service: ServiceIdentity) -> None:
        directory = self._identity_directory()
        destination = self._identity_path(service.unit)
        temporary = directory / f".{service.unit}.{os.getpid()}.tmp"
        content = json.dumps(
            service.to_dict(), sort_keys=True, separators=(",", ":")
        ).encode()
        descriptor = -1
        try:
            descriptor = os.open(
                temporary,
                os.O_WRONLY
                | os.O_CREAT
                | os.O_EXCL
                | getattr(os, "O_CLOEXEC", 0)
                | getattr(os, "O_NOFOLLOW", 0),
                0o600,
            )
            os.fchown(descriptor, self.owner_uid, self.owner_gid)
            _write_all(descriptor, content)
            os.fsync(descriptor)
            os.close(descriptor)
            descriptor = -1
            os.replace(temporary, destination)
            _fsync_directory(directory)
        except OSError:
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_UNAVAILABLE",
                "Service identity could not be persisted.",
            ) from None
        finally:
            if descriptor >= 0:
                os.close(descriptor)
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass

    def _read_identity(self, unit: str, *, required: bool) -> ServiceIdentity | None:
        path = self._identity_path(unit)
        try:
            descriptor = os.open(
                path,
                os.O_RDONLY
                | getattr(os, "O_CLOEXEC", 0)
                | getattr(os, "O_NOFOLLOW", 0),
            )
        except FileNotFoundError:
            if not required:
                return None
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_UNAVAILABLE",
                "Service identity is unavailable.",
            ) from None
        except OSError:
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_UNAVAILABLE",
                "Service identity is unavailable.",
            ) from None
        try:
            observed = os.fstat(descriptor)
            if (
                not stat.S_ISREG(observed.st_mode)
                or observed.st_nlink != 1
                or observed.st_uid != self.owner_uid
                or observed.st_gid != self.owner_gid
                or stat.S_IMODE(observed.st_mode) != 0o600
                or not 0 < observed.st_size <= _MAX_COMMAND_OUTPUT
            ):
                raise OSError
            content = os.read(descriptor, _MAX_COMMAND_OUTPUT + 1)
            if len(content) != observed.st_size:
                raise OSError
        except OSError:
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_UNAVAILABLE",
                "Service identity is unavailable.",
            ) from None
        finally:
            os.close(descriptor)
        try:
            return ServiceIdentity.from_dict(json.loads(content))
        except (DeploymentError, UnicodeDecodeError, json.JSONDecodeError):
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_UNAVAILABLE",
                "Service identity is unavailable.",
            ) from None

    def _identity_directory(self) -> Path:
        path = self.layout.service_identities
        try:
            observed = path.lstat()
        except OSError:
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_UNAVAILABLE",
                "Service identity boundary is unavailable.",
            ) from None
        if (
            not stat.S_ISDIR(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or observed.st_uid != self.owner_uid
            or observed.st_gid != self.owner_gid
            or stat.S_IMODE(observed.st_mode) != 0o700
        ):
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_UNAVAILABLE",
                "Service identity boundary is unavailable.",
            )
        return path


class ObservationProvider(Protocol):
    def observe(self, request: OperatorRequest) -> OperatorObservation: ...


class FixedObservationProvider:
    """Read one descriptor-pinned state/schema/service snapshot without writes."""

    def __init__(
        self,
        layout: DeploymentLayout,
        service_controller: ServiceController,
        *,
        root_uid: int = 0,
        root_gid: int = 0,
        operator_group_gid: int,
        service_uid: int,
        service_gid: int,
    ) -> None:
        self.layout = layout
        self.service_controller = service_controller
        self.root_uid = root_uid
        self.root_gid = root_gid
        self.operator_group_gid = operator_group_gid
        self.service_uid = service_uid
        self.service_gid = service_gid

    def observe(self, request: OperatorRequest) -> OperatorObservation:
        if request.operation != OBSERVE:
            raise fail("OPERATOR_REQUEST_INVALID", "Operator request is invalid.")
        state = self._read_state(request.deployment_identity)
        try:
            inspection = inspect_database(
                self.layout.database,
                expected_owner_uid=self.service_uid,
                expected_owner_gid=self.service_gid,
            )
        except DeploymentError:
            inspection = None
        schema_identity = (
            None if inspection is None else database_content_identity(inspection)
        )
        active = {
            component: state.components[component].active for component in COMPONENTS
        }
        unit_deployments = {
            "helixweave-api.service": active["platform"],
            "helixweave-worker.service": active["platform"],
            "helixweave-redis.service": active["platform"],
            "helixweave-docker-rootless.service": active["bulk-rnaseq-runtime"],
        }
        services: dict[str, str | None] = {}
        for unit in SERVICE_UNITS:
            unit_deployment = unit_deployments[unit]
            if unit_deployment is None:
                services[unit] = None
                continue
            service = self.service_controller.status(
                OperatorRequest(
                    operation=STATUS,
                    task_identity=request.task_identity,
                    deployment_identity=unit_deployment,
                    unit=unit,
                )
            )
            services[unit] = None if service is None else service.identity
        return OperatorObservation.create(
            state_identity=state.identity,
            active=active,
            database_schema_identity=schema_identity,
            database_schema_heads=(
                () if inspection is None else inspection.schema_heads
            ),
            services=services,
        )

    def _read_state(self, expected_identity: str) -> DeploymentState:
        try:
            state = StateStore(
                self.layout,
                reader_gid=self.operator_group_gid,
                service_gid=self.service_gid,
            ).read(
                expected_owner_uid=self.root_uid,
                expected_owner_gid=self.root_gid,
            )
        except DeploymentError as error:
            raise fail(
                "OPERATOR_STATE_UNTRUSTED",
                "Operator state is not trusted.",
                recoverable=error.issue.recoverable,
            ) from None
        if state.identity != expected_identity:
            raise fail(
                "OPERATOR_STATE_IDENTITY_MISMATCH",
                "Operator state identity does not match the request.",
                recoverable=True,
            )
        return state


class DeploymentActionController(Protocol):
    def execute(
        self,
        request: OperatorRequest,
        *,
        journal: OperatorJournalHandle,
    ) -> OperatorOutcome: ...


class DeploymentVerificationController(Protocol):
    def verify(self, request: OperatorRequest) -> DeploymentActionReceipt: ...


class UnavailableDeploymentActionController:
    """Activation stays fail-closed until unprivileged DB/state wiring is injected."""

    def execute(
        self,
        request: OperatorRequest,
        *,
        journal: OperatorJournalHandle,
    ) -> OperatorOutcome:
        del request, journal
        raise fail(
            "OPERATOR_DEPLOYMENT_ACTION_UNAVAILABLE",
            "Deployment activation backend is not installed.",
            recoverable=True,
        )


class DeploymentConfigurationController(Protocol):
    def activate(
        self,
        *,
        state: DeploymentState,
        api_contract_sha256: str,
    ) -> PlatformEnvironment: ...


class BoundaryUninstaller(Protocol):
    def uninstall(self, *, before_control_removal: Callable[[], None]) -> None: ...

    def recover(self) -> None: ...


UNINSTALL_STABLE_BOUNDARY_FILES = (
    Path("/usr/libexec/helixweave-operator"),
    Path("/usr/libexec/helixweave-gate-cleanup"),
    Path("/usr/libexec/helixweave-db-prepare"),
    Path("/usr/libexec/helixweave-operator-action"),
    Path("/usr/libexec/helixweave-encode-runtime-prepare"),
    Path("/usr/libexec/helixweave-bulk-runtime-prepare"),
    Path("/usr/libexec/helixweave-active-service"),
    Path("/etc/sudoers.d/helixweave-operator"),
)
UNINSTALL_LINKED_BOUNDARY_TARGETS = {
    Path("/usr/lib/tmpfiles.d/helixweave.conf"): Path(
        "/opt/helixweave/operator/current/templates/helixweave.tmpfiles.conf"
    ),
    Path("/usr/lib/systemd/system/helixweave-api.service"): Path(
        "/opt/helixweave/operator/current/templates/helixweave-api.service.in"
    ),
    Path("/usr/lib/systemd/system/helixweave-worker.service"): Path(
        "/opt/helixweave/operator/current/templates/helixweave-worker.service.in"
    ),
    Path("/usr/lib/systemd/system/helixweave-redis.service"): Path(
        "/opt/helixweave/operator/current/templates/helixweave-redis.service"
    ),
    Path("/usr/lib/systemd/system/helixweave-docker-rootless.service"): Path(
        "/opt/helixweave/operator/current/templates/helixweave-docker-rootless.service"
    ),
    Path("/usr/lib/systemd/system/helixweave-db-prepare.service"): Path(
        "/opt/helixweave/operator/current/templates/helixweave-db-prepare.service.in"
    ),
    Path("/usr/lib/systemd/system/helixweave-operator-action.service"): Path(
        "/opt/helixweave/operator/current/templates/helixweave-operator-action.service"
    ),
    Path("/usr/lib/systemd/system/helixweave-encode-runtime-prepare.service"): Path(
        "/opt/helixweave/operator/current/templates/helixweave-encode-runtime-prepare.service"
    ),
    Path("/usr/lib/systemd/system/helixweave-bulk-runtime-prepare.service"): Path(
        "/opt/helixweave/operator/current/templates/helixweave-bulk-runtime-prepare.service"
    ),
    Path("/usr/lib/systemd/system/helixweave.target"): Path(
        "/opt/helixweave/operator/current/templates/helixweave.target"
    ),
    Path("/etc/helixweave/redis.conf"): Path(
        "/opt/helixweave/operator/current/templates/redis.conf"
    ),
}
UNINSTALL_BOUNDARY_FILES = (
    *UNINSTALL_STABLE_BOUNDARY_FILES,
    *UNINSTALL_LINKED_BOUNDARY_TARGETS,
)
_UNINSTALL_CONTROL_FILES = frozenset(
    {
        Path("/usr/libexec/helixweave-operator"),
        Path("/etc/sudoers.d/helixweave-operator"),
    }
)


class HostBoundaryUninstaller:
    """Remove only the fixed executable/service boundary, never user data."""

    def __init__(
        self,
        *,
        owner_uid: int = 0,
        owner_gid: int = 0,
        executor: CommandExecutor | None = None,
        root_prefix: Path = Path("/"),
    ) -> None:
        if (
            not isinstance(owner_uid, int)
            or isinstance(owner_uid, bool)
            or owner_uid < 0
            or not isinstance(owner_gid, int)
            or isinstance(owner_gid, bool)
            or owner_gid < 0
            or not isinstance(root_prefix, Path)
            or not root_prefix.is_absolute()
        ):
            raise fail("OPERATOR_UNINSTALL_INVALID", "Uninstall boundary is invalid.")
        self.owner_uid = owner_uid
        self.owner_gid = owner_gid
        self.executor = BoundedCommandExecutor() if executor is None else executor
        self.root_prefix = root_prefix
        self.paths = {
            path: root_prefix / str(path).lstrip("/")
            for path in UNINSTALL_BOUNDARY_FILES
        }

    def uninstall(self, *, before_control_removal: Callable[[], None]) -> None:
        if not callable(before_control_removal):
            raise fail("OPERATOR_UNINSTALL_INVALID", "Uninstall boundary is invalid.")
        staged: list[tuple[Path, Path]] = []
        committed = False
        try:
            candidates = self._preflight()
            for logical, path in candidates:
                if logical in _UNINSTALL_CONTROL_FILES:
                    continue
                pending = self._pending_path(path)
                if path.exists() or path.is_symlink():
                    os.replace(path, pending)
                    _fsync_directory(path.parent)
                staged.append((path, pending))
            result = self.executor.run(
                (str(SYSTEMCTL), "daemon-reload"),
                timeout=_SYSTEMCTL_TIMEOUT_SECONDS,
            )
            if result.returncode != 0 or result.stdout:
                raise OSError
            before_control_removal()
            committed = True
            for _path, pending in staged:
                pending.unlink()
                _fsync_directory(pending.parent)
            for logical in (
                Path("/etc/sudoers.d/helixweave-operator"),
                Path("/usr/libexec/helixweave-operator"),
            ):
                path = self.paths[logical]
                try:
                    path.unlink()
                except FileNotFoundError:
                    pass
                else:
                    _fsync_directory(path.parent)
        except OSError:
            if not committed:
                for path, pending in reversed(staged):
                    try:
                        if pending.exists() or pending.is_symlink():
                            os.replace(pending, path)
                            _fsync_directory(path.parent)
                    except OSError:
                        pass
            raise fail(
                "OPERATOR_UNINSTALL_FAILED",
                "Operator boundary could not be uninstalled.",
                recoverable=True,
            ) from None

    def recover(self) -> None:
        """Restore any pre-terminal uninstall staging while controls still exist."""
        restored: set[Path] = set()
        try:
            candidates = self._preflight()
            for logical, path in candidates:
                if logical in _UNINSTALL_CONTROL_FILES:
                    continue
                pending = self._pending_path(path)
                if pending.exists() or pending.is_symlink():
                    os.replace(pending, path)
                    restored.add(path.parent)
            for parent in sorted(restored, key=str):
                _fsync_directory(parent)
            result = self.executor.run(
                (str(SYSTEMCTL), "daemon-reload"),
                timeout=_SYSTEMCTL_TIMEOUT_SECONDS,
            )
            if result.returncode != 0 or result.stdout:
                raise OSError
        except OSError:
            raise fail(
                "OPERATOR_UNINSTALL_FAILED",
                "Operator boundary could not be uninstalled.",
                recoverable=True,
            ) from None

    def _preflight(self) -> tuple[tuple[Path, Path], ...]:
        candidates: list[tuple[Path, Path]] = []
        for logical in UNINSTALL_BOUNDARY_FILES:
            path = self.paths[logical]
            pending = self._pending_path(path)
            present = path.exists() or path.is_symlink()
            pending_present = pending.exists() or pending.is_symlink()
            if present and pending_present:
                raise OSError
            selected = path if present else pending if pending_present else None
            if selected is None:
                continue
            observed = selected.lstat()
            if observed.st_uid != self.owner_uid or observed.st_gid != self.owner_gid:
                raise OSError
            expected_target = UNINSTALL_LINKED_BOUNDARY_TARGETS.get(logical)
            if expected_target is None:
                if (
                    not stat.S_ISREG(observed.st_mode)
                    or stat.S_ISLNK(observed.st_mode)
                    or observed.st_nlink != 1
                    or stat.S_IMODE(observed.st_mode) & 0o022
                ):
                    raise OSError
            elif (
                not stat.S_ISLNK(observed.st_mode)
                or observed.st_nlink != 1
                or os.readlink(selected) != expected_target.as_posix()
            ):
                raise OSError
            candidates.append((logical, path))
        return tuple(candidates)

    @staticmethod
    def _pending_path(path: Path) -> Path:
        return path.with_name(f".{path.name}.helixweave-uninstall-pending")


class AtomicPlatformConfiguration:
    """Derive, but never separately activate, one generation environment."""

    def __init__(
        self,
        layout: DeploymentLayout,
        *,
        owner_uid: int,
        service_gid: int,
    ) -> None:
        self.layout = layout
        self.owner_uid = owner_uid
        self.service_gid = service_gid

    def activate(
        self,
        *,
        state: DeploymentState,
        api_contract_sha256: str,
    ) -> PlatformEnvironment:
        return render_platform_environment(
            self.layout,
            state,
            api_contract_sha256=api_contract_sha256,
        )


class _RootSchemaAdmissionProvider(SchemaAdmissionProvider):
    def __init__(
        self,
        receipt: DeploymentActionReceipt,
        accepted_schema_heads: tuple[str, ...],
    ) -> None:
        self.provider_identity = canonical_identity(
            {
                "action_receipt_identity": receipt.identity,
                "accepted_schema_heads": list(accepted_schema_heads),
            },
            scheme="helixweave-root-database-schema-admission-v1",
        )
        self.accepted_schema_heads = accepted_schema_heads

    def admit(
        self,
        *,
        task_identity: str,
        deployment_identity: str,
        operator_state_identity: str,
        database_device: int,
        database_inode: int,
    ) -> DatabaseSchemaAdmission:
        return DatabaseSchemaAdmission.create(
            provider_identity=self.provider_identity,
            task_identity=task_identity,
            deployment_identity=deployment_identity,
            operator_state_identity=operator_state_identity,
            database_device=database_device,
            database_inode=database_inode,
            accepted_schema_heads=self.accepted_schema_heads,
        )


def _candidate_schema_target(
    receipt: DeploymentActionReceipt,
) -> tuple[str, ...]:
    database = receipt.readiness["database-schema"]
    if (
        receipt.database_before_identity is not None
        or receipt.database_after_identity is not None
        or receipt.accepted_schema_heads
        or not receipt.target_schema_heads
        or not receipt.migration_required
        or receipt.rollback_supported
        or database.status != "unavailable"
        or database.reason_code != "SCHEMA_INCOMPATIBLE"
        or database.identity is not None
    ):
        raise fail(
            "OPERATOR_ACTION_RECEIPT_INVALID",
            "Operator action receipt is invalid.",
        )
    return receipt.target_schema_heads


def _plan_schema_transition(
    *,
    operation: str,
    component: str,
    database_exists: bool,
    current_schema_heads: tuple[str, ...],
    target_schema_heads: tuple[str, ...],
    prior_platform_schema_heads: tuple[str, ...] | None,
    legacy_adoption: bool = False,
    known_schema_revisions: tuple[str, ...] = (),
) -> bool:
    """Return whether a forward platform migration is required."""
    if len(target_schema_heads) != 1:
        raise fail(
            "DEPLOYMENT_SCHEMA_INCOMPATIBLE",
            "Database schema is not compatible with the deployment.",
        )
    if not database_exists:
        if (
            operation == ACTIVATE
            and component == PLATFORM
            and not current_schema_heads
            and prior_platform_schema_heads is None
        ):
            return True
        raise fail(
            "DEPLOYMENT_SCHEMA_INCOMPATIBLE",
            "Database schema is not compatible with the deployment.",
        )
    if legacy_adoption:
        if (
            operation != ACTIVATE
            or component != PLATFORM
            or prior_platform_schema_heads is not None
            or len(current_schema_heads) != 1
            or current_schema_heads[0]
            not in {*known_schema_revisions, *target_schema_heads}
        ):
            raise fail(
                "DEPLOYMENT_SCHEMA_INCOMPATIBLE",
                "Database schema is not compatible with the deployment.",
            )
        return current_schema_heads != target_schema_heads
    if component != PLATFORM:
        if current_schema_heads == target_schema_heads:
            return False
    elif operation == ROLLBACK:
        if (
            prior_platform_schema_heads is not None
            and current_schema_heads == prior_platform_schema_heads
            and current_schema_heads == target_schema_heads
        ):
            return False
    elif prior_platform_schema_heads is not None:
        if current_schema_heads != prior_platform_schema_heads:
            raise fail(
                "DEPLOYMENT_SCHEMA_INCOMPATIBLE",
                "Database schema is not compatible with the deployment.",
            )
        return current_schema_heads != target_schema_heads
    raise fail(
        "DEPLOYMENT_SCHEMA_INCOMPATIBLE",
        "Database schema is not compatible with the deployment.",
    )


def _terminal_recovery_record(
    record: OperatorTransaction,
    *,
    phase: str,
) -> OperatorTransaction:
    return OperatorTransaction.create(
        request_identity=record.request_identity,
        operation=record.operation,
        task_identity=record.task_identity,
        deployment_identity=record.deployment_identity,
        component=record.component,
        unit=record.unit,
        phase=phase,
        failure_phase=record.failure_phase,
        write_fence=record.write_fence,
        point_of_no_return=record.point_of_no_return,
        restart_units=record.restart_units,
        prior_running_units=record.prior_running_units,
        prior_active=record.prior_active,
        candidate_active=record.candidate_active,
        database_mode=record.database_mode,
        prior_state_identity=record.prior_state_identity,
        candidate_state_identity=record.candidate_state_identity,
        backup_slot_identity=record.backup_slot_identity,
        backup_receipt_identity=record.backup_receipt_identity,
        source_database_identity=record.source_database_identity,
        schema_before_identity=record.schema_before_identity,
        schema_after_identity=record.schema_after_identity,
        schema_before_heads=record.schema_before_heads,
        target_schema_heads=record.target_schema_heads,
        evidence=record.evidence,
    )


class HostDeploymentActionController:
    """Root state/config/service transaction around unprivileged candidate code."""

    def __init__(
        self,
        layout: DeploymentLayout,
        *,
        states: StateStore,
        services: ServiceController,
        action_runner: DeploymentActionRunner | None = None,
        database_preparer: DatabasePreparer | None = None,
        encode_runtime_preparer: EncodeRuntimePreparer | None = None,
        bulk_runtime_preparer: BulkRuntimePreparer | None = None,
        configuration: DeploymentConfigurationController | None = None,
        boundary_uninstaller: BoundaryUninstaller | None = None,
        root_uid: int = 0,
        root_gid: int = 0,
        service_uid: int | None = None,
        service_gid: int | None = None,
        candidate_uid: int | None = None,
        candidate_gid: int | None = None,
    ) -> None:
        self.layout = layout
        self.states = states
        self.services = services
        self.action_runner = action_runner
        self.database_preparer = database_preparer
        self.encode_runtime_preparer = encode_runtime_preparer
        self.bulk_runtime_preparer = bulk_runtime_preparer
        self.configuration = configuration
        self.boundary_uninstaller = boundary_uninstaller
        self.root_uid = root_uid
        self.root_gid = root_gid
        self.service_uid = service_uid
        self.service_gid = service_gid
        self.candidate_uid = candidate_uid
        self.candidate_gid = candidate_gid

    def execute(
        self,
        request: OperatorRequest,
        *,
        journal: OperatorJournalHandle,
    ) -> OperatorOutcome:
        if request.operation == UNINSTALL:
            return self._uninstall(request, journal=journal)
        if request.operation not in {ACTIVATE, ROLLBACK} or request.component is None:
            raise fail("OPERATOR_REQUEST_INVALID", "Operator request is invalid.")
        (
            action_runner,
            database_preparer,
            encode_runtime_preparer,
            bulk_runtime_preparer,
            configuration,
        ) = self._production_boundaries()
        with self.states.transaction(
            exclusive=True,
            expected_owner_uid=self.root_uid,
            expected_owner_gid=self.root_gid,
        ) as transaction:
            state = transaction.read()
            if transaction.pending_transactions():
                raise fail(
                    "DEPLOYMENT_RECOVERY_REQUIRED",
                    "An interrupted deployment operation requires recovery.",
                    recoverable=True,
                )
            slots = state.components[request.component]
            if request.operation == ACTIVATE:
                if slots.staged != request.deployment_identity:
                    raise fail(
                        "DEPLOYMENT_STAGED_IDENTITY_MISMATCH",
                        "Staged deployment identity does not match the request.",
                    )
                candidate = state.activate(request.component)
            else:
                if slots.previous != request.deployment_identity:
                    raise fail(
                        "DEPLOYMENT_PREVIOUS_IDENTITY_MISMATCH",
                        "Previous deployment identity does not match the request.",
                    )
                candidate = state.rollback(request.component)
            journal.advance(
                "candidate-selected",
                prior_active={
                    component: state.components[component].active
                    for component in COMPONENTS
                },
                candidate_active={
                    component: candidate.components[component].active
                    for component in COMPONENTS
                },
                evidence={
                    "prior_state_identity": state.identity,
                    "candidate_state_identity": candidate.identity,
                },
            )
            authority = candidate.components[PLATFORM].active
            if authority is None:
                raise fail(
                    "OPERATOR_ACTION_UNAVAILABLE",
                    "Operator action is unavailable.",
                    recoverable=True,
                )
            active = {
                component: candidate.components[component].active
                for component in COMPONENTS
            }
            action_request = DeploymentActionRequest.create(
                phase="admit",
                operation=request.operation,
                component=request.component,
                task_identity=request.task_identity,
                deployment_identity=request.deployment_identity,
                authority_platform_identity=authority,
                prior_state_identity=state.identity,
                candidate_state_identity=candidate.identity,
                candidate_active=active,
            )
            if request.operation == ACTIVATE and request.component == ENCODE_RUNTIME:
                runtime_request = EncodeRuntimePrepareRequest.create(
                    task_identity=request.task_identity,
                    deployment_identity=request.deployment_identity,
                    authority_platform_identity=authority,
                    prior_state_identity=state.identity,
                    candidate_state_identity=candidate.identity,
                )
                journal.advance("runtime-prepare-started")
                runtime_receipt = encode_runtime_preparer.prepare(runtime_request)
                if (
                    runtime_receipt.request_identity != runtime_request.identity
                    or runtime_receipt.deployment_identity
                    != request.deployment_identity
                ):
                    raise fail(
                        "ENCODE_RUNTIME_PREPARE_RECEIPT_INVALID",
                        "ENCODE runtime preparation receipt is invalid.",
                    )
                journal.advance(
                    "runtime-prepared",
                    evidence={
                        "runtime_prepare_receipt_identity": runtime_receipt.identity,
                        "runtime_tree_identity": runtime_receipt.tree_identity,
                    },
                )
            prior_platform_schema_heads = None
            if request.component == PLATFORM and slots.active is not None:
                prior_platform_schema_heads = self._observe_platform_schema_target(
                    action_runner,
                    state=state,
                    task_identity=request.task_identity,
                )
            admission = action_runner.run(action_request)
            target_schema_heads = _candidate_schema_target(admission)
            partial_fresh_assembly = not all(
                state.components[component].active is not None
                for component in COMPONENTS
            ) and not all(active[component] is not None for component in COMPONENTS)
            native_checks = {
                PLATFORM: "platform-native",
                ENCODE_RUNTIME: "encode-runtime-native",
                BULK_RNASEQ_RUNTIME: "bulk-runtime-native",
            }
            if (
                admission.status != "admitted"
                or (
                    admission.compatibility != "compatible"
                    and not (
                        partial_fresh_assembly
                        and admission.compatibility == "incomplete"
                    )
                )
                or admission.request_identity != action_request.identity
                or any(
                    (active[component] is None)
                    != (admission.native_identities[component] is None)
                    for component in COMPONENTS
                )
                or any(
                    active[component] is not None
                    and (
                        admission.readiness[native_checks[component]].status != "ready"
                        or admission.readiness[native_checks[component]].identity
                        != admission.native_identities[component]
                    )
                    for component in COMPONENTS
                )
            ):
                raise fail(
                    "DEPLOYMENT_COMPATIBILITY_FAILED",
                    "Deployment components are not mutually compatible.",
                )
            was_service_ready = all(
                state.components[component].active is not None
                for component in COMPONENTS
            )
            candidate_service_ready = all(
                active[component] is not None for component in COMPONENTS
            )
            start_initial = not was_service_ready and candidate_service_ready
            stopped = self._stop_affected_services(
                state=state,
                candidate=candidate,
                component=request.component,
                task_identity=request.task_identity,
                journal=journal,
                start_initial=start_initial,
            )
            database_exists = self._database_exists()
            database_before = None
            database_before_identity = None
            current_schema_heads: tuple[str, ...] = ()
            if database_exists:
                assert self.service_uid is not None
                assert self.service_gid is not None
                database_before = inspect_database(
                    self.layout.database,
                    expected_owner_uid=self.service_uid,
                    expected_owner_gid=self.service_gid,
                )
                database_before_identity = database_content_identity(database_before)
                current_schema_heads = database_before.schema_heads
            legacy_adoption = (
                request.operation == ACTIVATE
                and request.component == PLATFORM
                and slots.active is None
                and database_exists
            )
            migration_required = _plan_schema_transition(
                operation=request.operation,
                component=request.component,
                database_exists=database_exists,
                current_schema_heads=current_schema_heads,
                target_schema_heads=target_schema_heads,
                prior_platform_schema_heads=prior_platform_schema_heads,
                legacy_adoption=legacy_adoption,
                known_schema_revisions=admission.known_schema_revisions,
            )
            recovery_evidence = (
                {}
                if database_before_identity is None
                else {"source_database_identity": database_before_identity}
            )
            journal.advance(
                "writers-stopped",
                write_fence=request.component == PLATFORM,
                database_mode=(
                    "fresh-candidate"
                    if request.component == PLATFORM and not database_exists
                    else (
                        "existing-live"
                        if request.component == PLATFORM and database_exists
                        else "none"
                    )
                ),
                schema_before_heads=current_schema_heads,
                target_schema_heads=target_schema_heads,
                evidence=recovery_evidence,
            )
            if request.component == BULK_RNASEQ_RUNTIME:
                if start_initial:
                    self._start_initial_redis(
                        candidate=candidate,
                        task_identity=request.task_identity,
                    )
                self._prepare_bulk_runtime(
                    operation=request.operation,
                    task_identity=request.task_identity,
                    candidate=candidate,
                    authority_platform_identity=authority,
                    prior_state_identity=state.identity,
                    preparer=bulk_runtime_preparer,
                    journal=journal,
                )
            backup_receipt_identity = None
            backup_required = database_exists and (
                migration_required or legacy_adoption
            )
            if backup_required:
                writers = self._writer_witnesses(
                    stopped,
                    deployment_identity=slots.active,
                    task_identity=request.task_identity,
                )
                witness = create_write_stop_witness(
                    self.layout,
                    task_identity=request.task_identity,
                    deployment_identity=request.deployment_identity,
                    operator_state_identity=state.identity,
                    writers=writers,
                    expected_database_uid=self.service_uid,
                    expected_database_gid=self.service_gid,
                )
                journal.advance(
                    "witness-ready",
                    evidence={"witness_identity": witness.identity},
                )
                backup = backup_database(
                    self.layout,
                    task_identity=request.task_identity,
                    expected_deployment_identity=request.deployment_identity,
                    expected_operator_state_identity=state.identity,
                    expected_database_uid=self.service_uid,
                    expected_database_gid=self.service_gid,
                    schema_provider=_RootSchemaAdmissionProvider(
                        admission,
                        current_schema_heads,
                    ),
                )
                backup_receipt_identity = backup.identity
                journal.advance(
                    "backup-committed",
                    evidence={
                        "backup_slot_identity": backup.backup_identity,
                        "backup_receipt_identity": backup.identity,
                        "source_database_identity": database_before_identity,
                    },
                )
            if migration_required:
                prepare_request = DatabasePrepareRequest.create(
                    operation=request.operation,
                    database_mode=(
                        "existing-live" if database_exists else "fresh-candidate"
                    ),
                    task_identity=request.task_identity,
                    deployment_identity=request.deployment_identity,
                    prior_state_identity=state.identity,
                    candidate_state_identity=candidate.identity,
                    action_receipt_identity=admission.identity,
                    backup_receipt_identity=backup_receipt_identity,
                    target_schema_heads=target_schema_heads,
                )
                journal.advance("migration-started")
                prepared = database_preparer.prepare(prepare_request)
                prepared_path = (
                    self.layout.database
                    if database_exists
                    else fresh_database_candidate_path(
                        self.layout, request.task_identity
                    )
                )
                after = inspect_database(
                    prepared_path,
                    expected_owner_uid=self.service_uid,
                    expected_owner_gid=self.service_gid,
                )
                if (
                    prepared.request_identity != prepare_request.identity
                    or prepared.database_after_identity
                    != database_content_identity(after)
                    or prepared.database_before_identity != database_before_identity
                    or prepared.schema_heads != target_schema_heads
                    or after.schema_heads != target_schema_heads
                ):
                    raise fail(
                        "DATABASE_PREPARE_RECEIPT_INVALID",
                        "Database preparation receipt is invalid.",
                    )
                journal.advance(
                    "migration-verified",
                    evidence={
                        "schema_after_identity": prepared.database_after_identity,
                        "database_prepare_receipt_identity": prepared.identity,
                    },
                )
                if not database_exists:
                    published = publish_fresh_database(
                        self.layout,
                        task_identity=request.task_identity,
                        expected_candidate_identity=prepared.database_after_identity,
                        target_schema_heads=target_schema_heads,
                        expected_database_uid=self.service_uid,
                        expected_database_gid=self.service_gid,
                    )
                    if (
                        database_content_identity(published)
                        != prepared.database_after_identity
                        or published.schema_heads != target_schema_heads
                    ):
                        raise fail(
                            "DATABASE_FRESH_PUBLISH_INVALID",
                            "Fresh database publication is invalid.",
                        )
                    journal.advance("database-published")
            prepared_configuration = configuration.activate(
                state=candidate,
                api_contract_sha256=admission.api_contract_sha256,
            )
            transaction.commit(
                candidate,
                operation=f"{request.operation}-{request.component}",
                expected_current_identity=state.identity,
                platform_environment=prepared_configuration,
            )
            journal.advance(
                "state-committed",
                evidence={
                    "prior_state_identity": state.identity,
                    "candidate_state_identity": candidate.identity,
                    "configuration_identity": prepared_configuration.identity,
                    "action_receipt_identity": admission.identity,
                },
            )
        self._restart_services(
            stopped,
            candidate=candidate,
            task_identity=request.task_identity,
            journal=journal,
            start_initial=start_initial,
        )
        return OperatorOutcome(
            "activated" if request.operation == ACTIVATE else "rolled-back"
        )

    def recover(self, record: OperatorTransaction) -> OperatorTransaction:
        """Synchronously restore pre-PONR work or finish a post-PONR candidate."""
        if (
            not isinstance(record, OperatorTransaction)
            or record.phase != "recovery-required"
            or record.failure_phase is None
        ):
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        if record.operation == STAGE:
            return self._recover_stage(record)
        if record.operation in {START, STOP, CLEANUP, UNINSTALL}:
            return self._recover_direct_operation(record)
        if (
            record.operation not in {ACTIVATE, ROLLBACK}
            or record.component not in COMPONENTS
            or record.prior_state_identity is None
            or record.candidate_state_identity is None
            or record.prior_active is None
            or record.candidate_active is None
        ):
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        if record.point_of_no_return:
            self._resume_candidate(record)
            return _terminal_recovery_record(record, phase="complete")
        self._restore_prior(record)
        return _terminal_recovery_record(record, phase="aborted")

    def _recover_stage(self, record: OperatorTransaction) -> OperatorTransaction:
        if (
            record.component not in COMPONENTS
            or record.unit is not None
            or record.prior_state_identity is None
            or record.candidate_state_identity is None
            or record.prior_active is None
            or record.candidate_active != record.prior_active
        ):
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        terminal_phase = "complete"
        with self.states.transaction(
            exclusive=True,
            expected_owner_uid=self.root_uid,
            expected_owner_gid=self.root_gid,
        ) as transaction:
            current = transaction.read()
            if transaction.pending_transactions():
                try:
                    selected = transaction.recover(
                        prior_state_identity=record.prior_state_identity,
                        candidate_state_identity=record.candidate_state_identity,
                        desired="complete-candidate",
                    )
                except DeploymentError as error:
                    if error.issue.code != "DEPLOYMENT_RECOVERY_REQUIRED":
                        raise
                    selected = transaction.recover(
                        prior_state_identity=record.prior_state_identity,
                        candidate_state_identity=record.candidate_state_identity,
                        desired="restore-prior",
                    )
                    terminal_phase = "aborted"
            elif current.identity == record.candidate_state_identity:
                selected = current
            elif current.identity == record.prior_state_identity:
                selected = current
                terminal_phase = "aborted"
            else:
                raise fail(
                    "OPERATOR_RECOVERY_REQUIRED",
                    "Operator transaction requires recovery.",
                    recoverable=True,
                )
        if selected is None:
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        active = {
            component: selected.components[component].active for component in COMPONENTS
        }
        if active != record.prior_active or (
            terminal_phase == "complete"
            and selected.components[record.component].staged
            != record.deployment_identity
        ):
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        return _terminal_recovery_record(record, phase=terminal_phase)

    def _recover_direct_operation(
        self,
        record: OperatorTransaction,
    ) -> OperatorTransaction:
        if (
            record.component is not None
            or record.prior_state_identity is None
            or record.candidate_state_identity != record.prior_state_identity
            or record.prior_active is None
            or record.candidate_active != record.prior_active
            or (record.operation == UNINSTALL) != (record.unit is None)
            or (record.operation != UNINSTALL) != (record.unit is not None)
            or (
                record.operation == UNINSTALL
                and (
                    record.deployment_identity != record.prior_state_identity
                    or record.restart_units != record.prior_running_units
                )
            )
            or (
                record.operation != UNINSTALL and record.restart_units != (record.unit,)
            )
            or (record.operation == START and record.prior_running_units)
            or (
                record.operation in {STOP, CLEANUP}
                and record.prior_running_units != (record.unit,)
            )
        ):
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        with self.states.transaction(
            exclusive=False,
            expected_owner_uid=self.root_uid,
            expected_owner_gid=self.root_gid,
        ) as transaction:
            state = transaction.read()
            if (
                transaction.pending_transactions()
                or state.identity != record.prior_state_identity
                or {
                    component: state.components[component].active
                    for component in COMPONENTS
                }
                != record.prior_active
            ):
                raise fail(
                    "OPERATOR_RECOVERY_REQUIRED",
                    "Operator transaction requires recovery.",
                    recoverable=True,
                )
        if record.operation == UNINSTALL:
            self._recover_uninstall(record)
            return _terminal_recovery_record(record, phase="aborted")
        else:
            self._recover_service_operation(record)
        return _terminal_recovery_record(record, phase="complete")

    def _recover_service_operation(self, record: OperatorTransaction) -> None:
        assert record.unit is not None
        assert record.prior_active is not None
        component = (
            BULK_RNASEQ_RUNTIME
            if record.unit == "helixweave-docker-rootless.service"
            else PLATFORM
        )
        if record.prior_active[component] != record.deployment_identity:
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        request = OperatorRequest(
            operation=STATUS,
            task_identity=record.task_identity,
            deployment_identity=record.deployment_identity,
            unit=record.unit,
        )
        if record.operation == START:
            if isinstance(self.services, SystemdServiceController):
                self.services.recover_start(
                    OperatorRequest(
                        operation=START,
                        task_identity=record.task_identity,
                        deployment_identity=record.deployment_identity,
                        unit=record.unit,
                    )
                )
                return
            observed = self.services.status(request)
            if observed is None:
                self.services.start(
                    OperatorRequest(
                        operation=START,
                        task_identity=record.task_identity,
                        deployment_identity=record.deployment_identity,
                        unit=record.unit,
                    )
                )
            return
        expected_service_identity = record.evidence.get("service_identity")
        if expected_service_identity is None:
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        if isinstance(self.services, SystemdServiceController):
            self.services.recover_stop(
                OperatorRequest(
                    operation=STOP,
                    task_identity=record.task_identity,
                    deployment_identity=record.deployment_identity,
                    unit=record.unit,
                    service_identity=expected_service_identity,
                ),
                cleanup=record.operation == CLEANUP,
            )
            return
        observed = self.services.status(request)
        if observed is None:
            return
        if observed.identity != expected_service_identity:
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        self.services.stop(
            OperatorRequest(
                operation=STOP,
                task_identity=observed.task_identity,
                deployment_identity=record.deployment_identity,
                unit=record.unit,
                service_identity=observed.identity,
            ),
            cleanup=record.operation == CLEANUP,
        )

    def _recover_uninstall(self, record: OperatorTransaction) -> None:
        assert record.prior_active is not None
        if self.boundary_uninstaller is None:
            self.boundary_uninstaller = HostBoundaryUninstaller(
                owner_uid=self.root_uid,
                owner_gid=self.root_gid,
            )
        self.boundary_uninstaller.recover()
        for unit in record.prior_running_units:
            deployment_identity = self._service_deployment(record.prior_active, unit)
            if deployment_identity is None:
                raise fail(
                    "OPERATOR_RECOVERY_REQUIRED",
                    "Operator transaction requires recovery.",
                    recoverable=True,
                )
            service = self._status_for_recovery(
                unit=unit,
                deployment_identity=deployment_identity,
                task_identity=record.task_identity,
                tolerate_identity_mismatch=False,
            )
            if service is None:
                self.services.start(
                    OperatorRequest(
                        operation=START,
                        task_identity=record.task_identity,
                        deployment_identity=deployment_identity,
                        unit=unit,
                    )
                )

    def _restore_prior(self, record: OperatorTransaction) -> None:
        assert record.prior_active is not None
        assert record.candidate_active is not None
        self._stop_candidate_services(record)
        self._restore_database_before_ponr(record)
        self.states.recover_pending_transaction(
            prior_state_identity=record.prior_state_identity,
            candidate_state_identity=record.candidate_state_identity,
            desired="restore-prior",
            expected_owner_uid=self.root_uid,
            expected_owner_gid=self.root_gid,
        )
        for unit in record.prior_running_units:
            deployment = self._service_deployment(record.prior_active, unit)
            if deployment is None:
                raise fail(
                    "OPERATOR_RECOVERY_REQUIRED",
                    "Operator transaction requires recovery.",
                    recoverable=True,
                )
            observed = self._status_for_recovery(
                unit=unit,
                deployment_identity=deployment,
                task_identity=record.task_identity,
                tolerate_identity_mismatch=False,
            )
            if observed is None:
                self.services.start(
                    OperatorRequest(
                        operation=START,
                        task_identity=record.task_identity,
                        deployment_identity=deployment,
                        unit=unit,
                    )
                )

    def _resume_candidate(self, record: OperatorTransaction) -> None:
        assert record.candidate_active is not None
        self.states.recover_pending_transaction(
            prior_state_identity=record.prior_state_identity,
            candidate_state_identity=record.candidate_state_identity,
            desired="complete-candidate",
            expected_owner_uid=self.root_uid,
            expected_owner_gid=self.root_gid,
        )
        if not self._database_exists():
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        assert self.service_uid is not None
        assert self.service_gid is not None
        database = inspect_database(
            self.layout.database,
            expected_owner_uid=self.service_uid,
            expected_owner_gid=self.service_gid,
        )
        expected_database_identity = (
            record.schema_after_identity or record.source_database_identity
        )
        if (
            expected_database_identity is None
            or database_content_identity(database) != expected_database_identity
            or database.schema_heads != record.target_schema_heads
        ):
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        for unit in record.restart_units:
            deployment = self._service_deployment(record.candidate_active, unit)
            if deployment is None:
                raise fail(
                    "OPERATOR_RECOVERY_REQUIRED",
                    "Operator transaction requires recovery.",
                    recoverable=True,
                )
            observed = self._status_for_recovery(
                unit=unit,
                deployment_identity=deployment,
                task_identity=record.task_identity,
                tolerate_identity_mismatch=False,
            )
            if observed is None:
                self.services.start(
                    OperatorRequest(
                        operation=START,
                        task_identity=record.task_identity,
                        deployment_identity=deployment,
                        unit=unit,
                    )
                )

    def _stop_candidate_services(self, record: OperatorTransaction) -> None:
        assert record.candidate_active is not None
        assert record.prior_active is not None
        for unit in record.restart_units:
            observed = self._observe_recovery_service(record, unit)
            if observed is None:
                continue
            deployment, service = observed
            self.services.stop(
                OperatorRequest(
                    operation=STOP,
                    task_identity=service.task_identity,
                    deployment_identity=deployment,
                    unit=unit,
                    service_identity=service.identity,
                ),
                cleanup=False,
            )
        bulk_unit = "helixweave-docker-rootless.service"
        if record.component != BULK_RNASEQ_RUNTIME or bulk_unit in record.restart_units:
            return
        deployment = record.candidate_active[BULK_RNASEQ_RUNTIME]
        if deployment is None:
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        observed = self._status_for_recovery(
            unit=bulk_unit,
            deployment_identity=deployment,
            task_identity=record.task_identity,
            tolerate_identity_mismatch=False,
        )
        if observed is None:
            return
        bulk_service_identity = record.evidence.get("bulk_runtime_service_identity")
        if (
            bulk_service_identity is not None
            and observed.identity != bulk_service_identity
        ):
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        self.services.stop(
            OperatorRequest(
                operation=STOP,
                task_identity=observed.task_identity,
                deployment_identity=deployment,
                unit=bulk_unit,
                service_identity=observed.identity,
            ),
            cleanup=False,
        )

    def _observe_recovery_service(
        self,
        record: OperatorTransaction,
        unit: str,
    ) -> tuple[str, ServiceIdentity] | None:
        assert record.prior_active is not None
        assert record.candidate_active is not None
        deployments: list[str] = []
        for active in (record.candidate_active, record.prior_active):
            deployment = self._service_deployment(active, unit)
            if deployment is not None and deployment not in deployments:
                deployments.append(deployment)
        for deployment in deployments:
            observed = self._status_for_recovery(
                unit=unit,
                deployment_identity=deployment,
                task_identity=record.task_identity,
                tolerate_identity_mismatch=True,
            )
            if observed is not None:
                return deployment, observed
        return None

    def _restore_database_before_ponr(self, record: OperatorTransaction) -> None:
        if record.database_mode == "none":
            return
        assert self.service_uid is not None
        assert self.service_gid is not None
        if record.database_mode == "existing-live":
            if record.backup_slot_identity is None:
                if (
                    record.source_database_identity is None
                    or not record.schema_before_heads
                    or not self._database_exists()
                ):
                    raise fail(
                        "OPERATOR_RECOVERY_REQUIRED",
                        "Operator transaction requires recovery.",
                        recoverable=True,
                    )
                unchanged = inspect_database(
                    self.layout.database,
                    expected_owner_uid=self.service_uid,
                    expected_owner_gid=self.service_gid,
                )
                if (
                    database_content_identity(unchanged)
                    != record.source_database_identity
                    or unchanged.schema_heads != record.schema_before_heads
                ):
                    raise fail(
                        "OPERATOR_RECOVERY_REQUIRED",
                        "Operator transaction requires recovery.",
                        recoverable=True,
                    )
                return
            if (
                record.source_database_identity is None
                or not record.schema_before_heads
            ):
                raise fail(
                    "OPERATOR_RECOVERY_REQUIRED",
                    "Operator transaction requires recovery.",
                    recoverable=True,
                )
            restore_database_backup(
                self.layout,
                backup_identity=record.backup_slot_identity,
                expected_task_identity=record.task_identity,
                expected_deployment_identity=record.deployment_identity,
                expected_prior_state_identity=record.prior_state_identity,
                expected_source_identity=record.source_database_identity,
                expected_schema_heads=record.schema_before_heads,
                expected_database_uid=self.service_uid,
                expected_database_gid=self.service_gid,
            )
            return
        if record.database_mode != "fresh-candidate":
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        candidate = fresh_database_candidate_path(self.layout, record.task_identity)
        candidate_present = candidate.exists() or candidate.is_symlink()
        canonical_present = self._database_exists()
        if not candidate_present and not canonical_present:
            return
        if record.schema_after_identity is None:
            if canonical_present or not candidate_present:
                raise fail(
                    "OPERATOR_RECOVERY_REQUIRED",
                    "Operator transaction requires recovery.",
                    recoverable=True,
                )
            quarantine_invalid_fresh_database(
                self.layout,
                task_identity=record.task_identity,
                expected_deployment_identity=record.deployment_identity,
                expected_prior_state_identity=record.prior_state_identity,
                expected_database_uid=self.service_uid,
                expected_database_gid=self.service_gid,
            )
            return
        if not record.target_schema_heads:
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        quarantine_fresh_database(
            self.layout,
            task_identity=record.task_identity,
            expected_deployment_identity=record.deployment_identity,
            expected_prior_state_identity=record.prior_state_identity,
            expected_candidate_identity=record.schema_after_identity,
            expected_schema_heads=record.target_schema_heads,
            expected_database_uid=self.service_uid,
            expected_database_gid=self.service_gid,
        )

    @staticmethod
    def _service_deployment(
        active: Mapping[str, str | None],
        unit: str,
    ) -> str | None:
        return active[
            BULK_RNASEQ_RUNTIME
            if unit == "helixweave-docker-rootless.service"
            else PLATFORM
        ]

    def _status_for_recovery(
        self,
        *,
        unit: str,
        deployment_identity: str,
        task_identity: str,
        tolerate_identity_mismatch: bool,
    ) -> ServiceIdentity | None:
        try:
            return self.services.status(
                OperatorRequest(
                    operation=STATUS,
                    task_identity=task_identity,
                    deployment_identity=deployment_identity,
                    unit=unit,
                )
            )
        except DeploymentError as error:
            if isinstance(
                self.services, SystemdServiceController
            ) and error.issue.code in {
                "OPERATOR_SERVICE_IDENTITY_MISMATCH",
                "OPERATOR_SERVICE_IDENTITY_UNAVAILABLE",
            }:
                return self.services.recover_observe(
                    OperatorRequest(
                        operation=STATUS,
                        task_identity=task_identity,
                        deployment_identity=deployment_identity,
                        unit=unit,
                    )
                )
            if (
                tolerate_identity_mismatch
                and error.issue.code == "OPERATOR_SERVICE_IDENTITY_MISMATCH"
            ):
                return None
            raise

    def verify(self, request: OperatorRequest) -> DeploymentActionReceipt:
        if request.operation != VERIFY:
            raise fail("OPERATOR_REQUEST_INVALID", "Operator request is invalid.")
        (
            action_runner,
            _database_preparer,
            _encode_runtime_preparer,
            bulk_runtime_preparer,
            _configuration,
        ) = self._production_boundaries()
        with self.states.transaction(
            exclusive=False,
            expected_owner_uid=self.root_uid,
            expected_owner_gid=self.root_gid,
        ) as transaction:
            state = transaction.read()
            if state.identity != request.deployment_identity:
                raise fail(
                    "OPERATOR_STATE_IDENTITY_MISMATCH",
                    "Operator state identity does not match the request.",
                    recoverable=True,
                )
            active = {
                component: state.components[component].active
                for component in COMPONENTS
            }
            authority = active[PLATFORM]
            if authority is None:
                raise fail(
                    "OPERATOR_ACTION_UNAVAILABLE",
                    "Operator action is unavailable.",
                    recoverable=True,
                )
            action_request = DeploymentActionRequest.create(
                phase="observe",
                operation="observe",
                component=PLATFORM,
                task_identity=request.task_identity,
                deployment_identity=authority,
                authority_platform_identity=authority,
                prior_state_identity=state.identity,
                candidate_state_identity=state.identity,
                candidate_active=active,
            )
            receipt = action_runner.run(action_request)
            target_schema_heads = _candidate_schema_target(receipt)
            database_exists = self._database_exists()
            database = (
                inspect_database(
                    self.layout.database,
                    expected_owner_uid=self.service_uid,
                    expected_owner_gid=self.service_gid,
                )
                if database_exists
                else None
            )
            native_checks = {
                PLATFORM: "platform-native",
                ENCODE_RUNTIME: "encode-runtime-native",
                BULK_RNASEQ_RUNTIME: "bulk-runtime-native",
            }
            if (
                receipt.status != "observed"
                or receipt.request_identity != action_request.identity
                or receipt.compatibility != "compatible"
                or database is None
                or database.schema_heads != target_schema_heads
                or any(
                    (active[component] is None)
                    != (receipt.native_identities[component] is None)
                    for component in COMPONENTS
                )
                or any(
                    active[component] is not None
                    and (
                        receipt.readiness[native_checks[component]].status != "ready"
                        or receipt.readiness[native_checks[component]].identity
                        != receipt.native_identities[component]
                    )
                    for component in COMPONENTS
                )
            ):
                raise fail(
                    "OPERATOR_ACTION_RECEIPT_INVALID",
                    "Operator action receipt is invalid.",
                )
            database_identity = database_content_identity(database)
            readiness = dict(receipt.readiness)
            readiness["database-schema"] = ReadinessCheck(
                "ready", "READY", database_identity
            )
            readiness["configuration"] = self._configuration_readiness(
                state,
                api_contract_sha256=receipt.api_contract_sha256,
            )
            readiness["permissions"] = self._storage_readiness(active)
            readiness["redis"] = self._redis_readiness(
                deployment_identity=active[PLATFORM],
                task_identity=request.task_identity,
            )
            readiness["docker"] = self._bulk_runtime_readiness(
                state=state,
                task_identity=request.task_identity,
                preparer=bulk_runtime_preparer,
            )
            readiness["worker"] = self._service_readiness(
                unit="helixweave-worker.service",
                deployment_identity=active[PLATFORM],
                task_identity=request.task_identity,
            )
            readiness["references"] = _reference_api_readiness()
            return DeploymentActionReceipt.create(
                request_identity=receipt.request_identity,
                status=receipt.status,
                compatibility=receipt.compatibility,
                database_before_identity=None,
                database_after_identity=database_identity,
                accepted_schema_heads=database.schema_heads,
                target_schema_heads=target_schema_heads,
                migration_inventory_identity=receipt.migration_inventory_identity,
                known_schema_revisions=receipt.known_schema_revisions,
                migration_required=False,
                rollback_supported=True,
                api_contract_sha256=receipt.api_contract_sha256,
                native_identities=receipt.native_identities,
                frontend_identity=receipt.frontend_identity,
                reference_compatibility_identity=(
                    receipt.reference_compatibility_identity
                ),
                readiness=readiness,
            )

    def _configuration_readiness(
        self,
        state: DeploymentState,
        *,
        api_contract_sha256: str,
    ) -> ReadinessCheck:
        expected = render_platform_environment(
            self.layout,
            state,
            api_contract_sha256=api_contract_sha256,
        )
        path = self.layout.state_generations / state.identity / PLATFORM_ENV_FILENAME
        try:
            content, observed = read_regular_file(
                path,
                max_bytes=MAX_PLATFORM_ENV_BYTES,
                code="DEPLOYMENT_CONFIGURATION_INVALID",
            )
        except DeploymentError:
            return ReadinessCheck("not-ready", "CONFIGURATION_INVALID")
        if (
            content != expected.content
            or observed.st_uid != self.root_uid
            or observed.st_gid != self.service_gid
            or observed.st_nlink != 1
            or stat.S_IMODE(observed.st_mode) != 0o440
        ):
            return ReadinessCheck("not-ready", "CONFIGURATION_INVALID")
        return ReadinessCheck("ready", "READY", expected.identity)

    def _storage_readiness(
        self,
        active: Mapping[str, str | None],
    ) -> ReadinessCheck:
        evidence: list[dict[str, str]] = []
        try:
            store = BundleStore(self.layout)
            for component in COMPONENTS:
                identity = active[component]
                if identity is None:
                    continue
                manifest = store.verify_installed(
                    component,
                    identity,
                    expected_owner_uid=self.root_uid,
                    expected_owner_gid=self.root_gid,
                )
                if manifest.identity != identity:
                    raise ValueError
                evidence.append(
                    {
                        "component": component,
                        "deployment_identity": identity,
                    }
                )
        except (DeploymentError, OSError, ValueError):
            return ReadinessCheck("not-ready", "PERMISSION_INVALID")
        return ReadinessCheck(
            "ready",
            "READY",
            canonical_identity(
                evidence,
                scheme="helixweave-root-storage-readiness-v1",
            ),
        )

    def _service_readiness(
        self,
        *,
        unit: str,
        deployment_identity: str | None,
        task_identity: str,
    ) -> ReadinessCheck:
        if deployment_identity is None:
            return ReadinessCheck("not-applicable", "COMPONENT_NOT_ACTIVE")
        try:
            service = self.services.status(
                OperatorRequest(
                    operation=STATUS,
                    unit=unit,
                    deployment_identity=deployment_identity,
                    task_identity=task_identity,
                )
            )
        except DeploymentError:
            service = None
        if service is None:
            return ReadinessCheck("not-ready", "SERVICE_STOPPED")
        return ReadinessCheck("ready", "READY", service.identity)

    def _redis_readiness(
        self,
        *,
        deployment_identity: str | None,
        task_identity: str,
    ) -> ReadinessCheck:
        service = self._service_readiness(
            unit="helixweave-redis.service",
            deployment_identity=deployment_identity,
            task_identity=task_identity,
        )
        if service.status != "ready":
            return ReadinessCheck("unavailable", "REDIS_UNAVAILABLE")
        protocol = _unix_protocol_readiness(
            self.layout.run_root / "redis" / "redis.sock",
            request=b"*1\r\n$4\r\nPING\r\n",
            response_validator=lambda value: value == b"+PONG\r\n",
            unavailable_reason="REDIS_UNAVAILABLE",
            scheme="helixweave-redis-ping-observation-v1",
        )
        if protocol.status != "ready":
            return protocol
        return ReadinessCheck(
            "ready",
            "READY",
            canonical_identity(
                {
                    "service_identity": service.identity,
                    "protocol_identity": protocol.identity,
                },
                scheme="helixweave-root-service-readiness-v1",
            ),
        )

    def _bulk_runtime_readiness(
        self,
        *,
        state: DeploymentState,
        task_identity: str,
        preparer: BulkRuntimePreparer,
    ) -> ReadinessCheck:
        platform_identity = state.components[PLATFORM].active
        bulk_identity = state.components[BULK_RNASEQ_RUNTIME].active
        if bulk_identity is None:
            return ReadinessCheck("not-applicable", "COMPONENT_NOT_ACTIVE")
        if platform_identity is None:
            return ReadinessCheck("unavailable", "DOCKER_UNAVAILABLE")
        service = self._service_readiness(
            unit="helixweave-docker-rootless.service",
            deployment_identity=bulk_identity,
            task_identity=task_identity,
        )
        if service.status != "ready" or service.identity is None:
            return ReadinessCheck("unavailable", "DOCKER_UNAVAILABLE")
        try:
            boundary_before = preparer.observe_boundary()
            bulk_request = BulkRuntimePrepareRequest.create(
                operation="verify",
                task_identity=task_identity,
                candidate_bulk_identity=bulk_identity,
                authority_platform_identity=platform_identity,
                prior_state_identity=state.identity,
                candidate_state_identity=state.identity,
                docker_service_identity=service.identity,
                docker_client_identity=boundary_before.client_identity,
                docker_endpoint_identity=boundary_before.endpoint_identity,
                docker_daemon_uid=boundary_before.daemon_uid,
                docker_daemon_gid=boundary_before.daemon_gid,
            )
            receipt = preparer.prepare(bulk_request)
            observed = self.services.status(
                OperatorRequest(
                    operation=STATUS,
                    task_identity=task_identity,
                    deployment_identity=bulk_identity,
                    unit="helixweave-docker-rootless.service",
                )
            )
            boundary_after = preparer.observe_boundary()
            if (
                observed is None
                or observed.identity != service.identity
                or boundary_after != boundary_before
                or receipt.request_identity != bulk_request.identity
                or receipt.candidate_bulk_identity != bulk_identity
            ):
                raise ValueError
        except Exception:
            return ReadinessCheck("unavailable", "DOCKER_UNAVAILABLE")
        return ReadinessCheck(
            "ready",
            "READY",
            canonical_identity(
                {
                    "service_identity": service.identity,
                    "client_identity": boundary_before.client_identity,
                    "endpoint_identity": boundary_before.endpoint_identity,
                    "runtime_identity": receipt.runtime_identity,
                    "image_set_identity": receipt.image_set_identity,
                    "image_count": receipt.image_count,
                },
                scheme="helixweave-root-bulk-runtime-readiness-v1",
            ),
        )

    def _uninstall(
        self,
        request: OperatorRequest,
        *,
        journal: OperatorJournalHandle,
    ) -> OperatorOutcome:
        if request.operation != UNINSTALL or request.component is not None:
            raise fail("OPERATOR_REQUEST_INVALID", "Operator request is invalid.")
        with self.states.transaction(
            exclusive=True,
            expected_owner_uid=self.root_uid,
            expected_owner_gid=self.root_gid,
        ) as transaction:
            state = transaction.read()
            if state.identity != request.deployment_identity:
                raise fail(
                    "OPERATOR_STATE_IDENTITY_MISMATCH",
                    "Operator state identity does not match the request.",
                    recoverable=True,
                )
            if transaction.pending_transactions():
                raise fail(
                    "DEPLOYMENT_RECOVERY_REQUIRED",
                    "An interrupted deployment operation requires recovery.",
                    recoverable=True,
                )
            active_platform = state.components[PLATFORM].active
            active_bulk = state.components[BULK_RNASEQ_RUNTIME].active
            active = {
                component: state.components[component].active
                for component in COMPONENTS
            }
            running: dict[str, ServiceIdentity] = {}
            for unit in SERVICE_UNITS:
                deployment_identity = (
                    active_bulk
                    if unit == "helixweave-docker-rootless.service"
                    else active_platform
                )
                if deployment_identity is None:
                    continue
                service = self.services.status(
                    OperatorRequest(
                        operation=STATUS,
                        task_identity=request.task_identity,
                        deployment_identity=deployment_identity,
                        unit=unit,
                    )
                )
                if service is None:
                    continue
                running[unit] = service
            running_units = tuple(unit for unit in SERVICE_UNITS if unit in running)
            journal.advance(
                "service-stopping",
                restart_units=running_units,
                prior_running_units=running_units,
                prior_active=active,
                candidate_active=active,
                evidence={
                    "prior_state_identity": state.identity,
                    "candidate_state_identity": state.identity,
                },
            )
            for unit in running_units:
                service = running[unit]
                deployment_identity = (
                    active_bulk
                    if unit == "helixweave-docker-rootless.service"
                    else active_platform
                )
                assert deployment_identity is not None
                self.services.stop(
                    OperatorRequest(
                        operation=STOP,
                        task_identity=service.task_identity,
                        deployment_identity=deployment_identity,
                        unit=unit,
                        service_identity=service.identity,
                    ),
                    cleanup=False,
                )
            journal.advance("writers-stopped", write_fence=True)
        if self.boundary_uninstaller is None:
            self.boundary_uninstaller = HostBoundaryUninstaller(
                owner_uid=self.root_uid,
                owner_gid=self.root_gid,
            )
        self.boundary_uninstaller.uninstall(
            before_control_removal=journal.complete,
        )
        return OperatorOutcome("uninstalled")

    def _database_exists(self) -> bool:
        try:
            observed = self.layout.database.lstat()
        except FileNotFoundError:
            return False
        except OSError:
            raise fail(
                "DATABASE_INSPECTION_FAILED",
                "Database identity could not be inspected.",
                recoverable=True,
            ) from None
        if not stat.S_ISREG(observed.st_mode) or stat.S_ISLNK(observed.st_mode):
            raise fail(
                "DATABASE_INSPECTION_FAILED",
                "Database identity could not be inspected.",
            )
        return True

    def _observe_platform_schema_target(
        self,
        action_runner: DeploymentActionRunner,
        *,
        state: DeploymentState,
        task_identity: str,
    ) -> tuple[str, ...]:
        authority = state.components[PLATFORM].active
        if authority is None:
            raise fail(
                "DEPLOYMENT_SCHEMA_INCOMPATIBLE",
                "Database schema is not compatible with the deployment.",
            )
        active = {
            component: state.components[component].active for component in COMPONENTS
        }
        action_request = DeploymentActionRequest.create(
            phase="observe",
            operation="observe",
            component=PLATFORM,
            task_identity=task_identity,
            deployment_identity=authority,
            authority_platform_identity=authority,
            prior_state_identity=state.identity,
            candidate_state_identity=state.identity,
            candidate_active=active,
        )
        receipt = action_runner.run(action_request)
        if (
            receipt.status != "observed"
            or receipt.request_identity != action_request.identity
            or any(
                (active[component] is None)
                != (receipt.native_identities[component] is None)
                for component in COMPONENTS
            )
        ):
            raise fail(
                "OPERATOR_ACTION_RECEIPT_INVALID",
                "Operator action receipt is invalid.",
            )
        return _candidate_schema_target(receipt)

    def _production_boundaries(
        self,
    ) -> tuple[
        DeploymentActionRunner,
        DatabasePreparer,
        EncodeRuntimePreparer,
        BulkRuntimePreparer,
        DeploymentConfigurationController,
    ]:
        needs_service = (
            self.service_uid is None
            or self.service_gid is None
            or self.database_preparer is None
            or self.bulk_runtime_preparer is None
            or self.configuration is None
        )
        needs_candidate = (
            self.action_runner is None or self.encode_runtime_preparer is None
        )
        if needs_service:
            self.service_uid, self.service_gid = _resolve_account_identity(
                user="helixweave",
                group="helixweave",
                uid=self.service_uid,
                gid=self.service_gid,
            )
        if needs_candidate:
            self.candidate_uid, self.candidate_gid = _resolve_account_identity(
                user="helixweave-candidate",
                group="helixweave-candidate",
                uid=self.candidate_uid,
                gid=self.candidate_gid,
            )
        if (
            needs_service
            and needs_candidate
            and (
                self.service_uid == self.candidate_uid
                or self.service_gid == self.candidate_gid
            )
        ):
            raise fail(
                "OPERATOR_ACCOUNT_BOUNDARY_UNAVAILABLE",
                "Operator account boundary is unavailable.",
            )
        if self.action_runner is None:
            assert self.candidate_uid is not None
            assert self.candidate_gid is not None
            self.action_runner = SystemdDeploymentActionRunner(
                self.layout,
                service_uid=self.candidate_uid,
                service_gid=self.candidate_gid,
                root_uid=self.root_uid,
                root_gid=self.root_gid,
            )
        if self.database_preparer is None:
            assert self.service_uid is not None
            assert self.service_gid is not None
            self.database_preparer = SystemdDatabasePreparer(
                self.layout,
                service_uid=self.service_uid,
                service_gid=self.service_gid,
                root_uid=self.root_uid,
            )
        if self.encode_runtime_preparer is None:
            assert self.candidate_uid is not None
            assert self.candidate_gid is not None
            self.encode_runtime_preparer = SystemdEncodeRuntimePreparer(
                self.layout,
                service_uid=self.candidate_uid,
                service_gid=self.candidate_gid,
                root_uid=self.root_uid,
                root_gid=self.root_gid,
            )
        if self.bulk_runtime_preparer is None:
            assert self.service_uid is not None
            assert self.service_gid is not None
            self.bulk_runtime_preparer = SystemdBulkRuntimePreparer(
                self.layout,
                service_uid=self.service_uid,
                service_gid=self.service_gid,
                root_uid=self.root_uid,
                root_gid=self.root_gid,
            )
        if self.configuration is None:
            assert self.service_gid is not None
            self.configuration = AtomicPlatformConfiguration(
                self.layout,
                owner_uid=self.root_uid,
                service_gid=self.service_gid,
            )
        return (
            self.action_runner,
            self.database_preparer,
            self.encode_runtime_preparer,
            self.bulk_runtime_preparer,
            self.configuration,
        )

    def _start_initial_redis(
        self,
        *,
        candidate: DeploymentState,
        task_identity: str,
    ) -> ServiceIdentity:
        platform_identity = candidate.components[PLATFORM].active
        if platform_identity is None:
            raise fail(
                "OPERATOR_SERVICE_START_FAILED",
                "Service did not enter the running state.",
                recoverable=True,
            )
        unit = "helixweave-redis.service"
        status_request = OperatorRequest(
            operation=STATUS,
            task_identity=task_identity,
            deployment_identity=platform_identity,
            unit=unit,
        )
        observed = self.services.status(status_request)
        if observed is not None:
            return observed
        return self.services.start(
            OperatorRequest(
                operation=START,
                task_identity=task_identity,
                deployment_identity=platform_identity,
                unit=unit,
            )
        )

    def _prepare_bulk_runtime(
        self,
        *,
        operation: str,
        task_identity: str,
        candidate: DeploymentState,
        authority_platform_identity: str,
        prior_state_identity: str,
        preparer: BulkRuntimePreparer,
        journal: OperatorJournalHandle,
    ) -> BulkRuntimePrepareReceipt:
        candidate_bulk = candidate.components[BULK_RNASEQ_RUNTIME].active
        if candidate_bulk is None or operation not in {ACTIVATE, ROLLBACK}:
            raise fail(
                "BULK_RUNTIME_PREPARE_REQUEST_INVALID",
                "Bulk runtime preparation request is invalid.",
            )
        unit = "helixweave-docker-rootless.service"
        service = self.services.start(
            OperatorRequest(
                operation=START,
                task_identity=task_identity,
                deployment_identity=candidate_bulk,
                unit=unit,
            )
        )
        journal.advance(
            "bulk-runtime-prepare-started",
            evidence={"bulk_runtime_service_identity": service.identity},
        )
        boundary_before = preparer.observe_boundary()
        bulk_request = BulkRuntimePrepareRequest.create(
            operation=operation,
            task_identity=task_identity,
            candidate_bulk_identity=candidate_bulk,
            authority_platform_identity=authority_platform_identity,
            prior_state_identity=prior_state_identity,
            candidate_state_identity=candidate.identity,
            docker_service_identity=service.identity,
            docker_client_identity=boundary_before.client_identity,
            docker_endpoint_identity=boundary_before.endpoint_identity,
            docker_daemon_uid=boundary_before.daemon_uid,
            docker_daemon_gid=boundary_before.daemon_gid,
        )
        receipt = preparer.prepare(bulk_request)
        observed = self.services.status(
            OperatorRequest(
                operation=STATUS,
                task_identity=task_identity,
                deployment_identity=candidate_bulk,
                unit=unit,
            )
        )
        boundary_after = preparer.observe_boundary()
        if (
            observed is None
            or observed.identity != service.identity
            or boundary_after != boundary_before
            or receipt.request_identity != bulk_request.identity
            or receipt.candidate_bulk_identity != candidate_bulk
        ):
            raise fail(
                "BULK_RUNTIME_PREPARE_RECEIPT_INVALID",
                "Bulk runtime preparation receipt is invalid.",
            )
        journal.advance(
            "bulk-runtime-prepared",
            evidence={
                "bulk_runtime_prepare_receipt_identity": receipt.identity,
                "bulk_runtime_image_set_identity": receipt.image_set_identity,
            },
        )
        if unit not in journal.record.restart_units:
            self.services.stop(
                OperatorRequest(
                    operation=STOP,
                    task_identity=service.task_identity,
                    deployment_identity=candidate_bulk,
                    unit=unit,
                    service_identity=service.identity,
                ),
                cleanup=False,
            )
        return receipt

    def _stop_affected_services(
        self,
        *,
        state: DeploymentState,
        candidate: DeploymentState,
        component: str,
        task_identity: str,
        journal: OperatorJournalHandle,
        start_initial: bool,
    ) -> dict[str, ServiceIdentity | None]:
        units = (
            SERVICE_STOP_ORDER
            if start_initial
            else {
                PLATFORM: WRITER_UNITS,
                ENCODE_RUNTIME: WRITER_UNITS,
                BULK_RNASEQ_RUNTIME: (
                    "helixweave-api.service",
                    "helixweave-worker.service",
                    "helixweave-docker-rootless.service",
                ),
            }[component]
        )
        active_platform = state.components[PLATFORM].active
        active_bulk = state.components[BULK_RNASEQ_RUNTIME].active
        values: dict[str, ServiceIdentity | None] = {}
        for unit in units:
            deployment = (
                active_bulk
                if unit == "helixweave-docker-rootless.service"
                else active_platform
            )
            if deployment is None:
                values[unit] = None
                continue
            status_request = OperatorRequest(
                operation=STATUS,
                task_identity=task_identity,
                deployment_identity=deployment,
                unit=unit,
            )
            service = self.services.status(status_request)
            values[unit] = service
        restart_units = self._restart_unit_set(
            values,
            candidate=candidate,
            start_initial=start_initial,
        )
        prior_running_units = tuple(
            unit for unit in SERVICE_UNITS if values.get(unit) is not None
        )
        journal.advance(
            "service-stopping",
            restart_units=restart_units,
            prior_running_units=prior_running_units,
        )
        for unit, service in values.items():
            if service is None:
                continue
            deployment = (
                active_bulk
                if unit == "helixweave-docker-rootless.service"
                else active_platform
            )
            assert deployment is not None
            self.services.stop(
                OperatorRequest(
                    operation=STOP,
                    task_identity=service.task_identity,
                    deployment_identity=deployment,
                    unit=unit,
                    service_identity=service.identity,
                ),
                cleanup=False,
            )
        return values

    @staticmethod
    def _restart_unit_set(
        stopped: dict[str, ServiceIdentity | None],
        *,
        candidate: DeploymentState,
        start_initial: bool,
    ) -> tuple[str, ...]:
        return tuple(
            unit
            for unit in SERVICE_UNITS
            if unit in stopped
            and (stopped[unit] is not None or start_initial)
            and (
                candidate.components[BULK_RNASEQ_RUNTIME].active
                if unit == "helixweave-docker-rootless.service"
                else candidate.components[PLATFORM].active
            )
            is not None
        )

    @staticmethod
    def _writer_witnesses(
        stopped: dict[str, ServiceIdentity | None],
        *,
        deployment_identity: str | None,
        task_identity: str,
    ) -> tuple[StoppedWriter, ...]:
        values: list[StoppedWriter] = []
        for unit in ("helixweave-api.service", "helixweave-worker.service"):
            service = stopped.get(unit)
            identity = (
                service.identity
                if service is not None
                else canonical_identity(
                    {
                        "unit": unit,
                        "state": "already-stopped",
                        "deployment_identity": deployment_identity,
                        "task_identity": task_identity,
                    },
                    scheme="helixweave-stopped-service-witness-v1",
                )
            )
            values.append(StoppedWriter(unit, "stopped", identity))
        return tuple(values)

    def _restart_services(
        self,
        stopped: dict[str, ServiceIdentity | None],
        *,
        candidate: DeploymentState,
        task_identity: str,
        journal: OperatorJournalHandle,
        start_initial: bool,
    ) -> None:
        restart_units = self._restart_unit_set(
            stopped,
            candidate=candidate,
            start_initial=start_initial,
        )
        if journal.record.restart_units != restart_units:
            raise fail("OPERATOR_JOURNAL_INVALID", "Operator journal is invalid.")
        dependency_units = tuple(
            unit for unit in restart_units if unit in DEPENDENCY_UNITS
        )
        writer_units = tuple(unit for unit in restart_units if unit in WRITER_UNITS)
        last_service: ServiceIdentity | None = None
        if dependency_units:
            journal.advance("dependency-starting")
        for unit in dependency_units:
            deployment = (
                candidate.components[BULK_RNASEQ_RUNTIME].active
                if unit == "helixweave-docker-rootless.service"
                else candidate.components[PLATFORM].active
            )
            assert deployment is not None
            status_request = OperatorRequest(
                operation=STATUS,
                task_identity=task_identity,
                deployment_identity=deployment,
                unit=unit,
            )
            last_service = self.services.status(status_request)
            if last_service is None:
                last_service = self.services.start(
                    OperatorRequest(
                        operation=START,
                        task_identity=task_identity,
                        deployment_identity=deployment,
                        unit=unit,
                    )
                )
        if writer_units:
            journal.advance("service-starting", point_of_no_return=True)
        for unit in writer_units:
            deployment = (
                candidate.components[BULK_RNASEQ_RUNTIME].active
                if unit == "helixweave-docker-rootless.service"
                else candidate.components[PLATFORM].active
            )
            assert deployment is not None
            last_service = self.services.start(
                OperatorRequest(
                    operation=START,
                    task_identity=task_identity,
                    deployment_identity=deployment,
                    unit=unit,
                )
            )
        journal.advance(
            "service-observed",
            evidence=(
                {}
                if last_service is None
                else {"service_identity": last_service.identity}
            ),
        )


def _resolve_account_identity(
    *,
    user: str,
    group: str,
    uid: int | None,
    gid: int | None,
) -> tuple[int, int]:
    if (uid is None) != (gid is None):
        raise fail(
            "OPERATOR_ACCOUNT_BOUNDARY_UNAVAILABLE",
            "Operator account boundary is unavailable.",
        )
    if uid is None:
        try:
            account = pwd.getpwnam(user)
            group_entry = grp.getgrnam(group)
        except KeyError:
            raise fail(
                "OPERATOR_ACCOUNT_BOUNDARY_UNAVAILABLE",
                "Operator account boundary is unavailable.",
                recoverable=True,
            ) from None
        if (
            account.pw_name != user
            or group_entry.gr_name != group
            or account.pw_gid != group_entry.gr_gid
        ):
            raise fail(
                "OPERATOR_ACCOUNT_BOUNDARY_UNAVAILABLE",
                "Operator account boundary is unavailable.",
            )
        try:
            memberships = set(os.getgrouplist(user, group_entry.gr_gid))
        except (KeyError, OSError):
            raise fail(
                "OPERATOR_ACCOUNT_BOUNDARY_UNAVAILABLE",
                "Operator account boundary is unavailable.",
            ) from None
        if memberships != {group_entry.gr_gid}:
            raise fail(
                "OPERATOR_ACCOUNT_BOUNDARY_UNAVAILABLE",
                "Operator account boundary is unavailable.",
            )
        uid = account.pw_uid
        gid = group_entry.gr_gid
    if (
        not isinstance(uid, int)
        or isinstance(uid, bool)
        or uid <= 0
        or not isinstance(gid, int)
        or isinstance(gid, bool)
        or gid <= 0
    ):
        raise fail(
            "OPERATOR_ACCOUNT_BOUNDARY_UNAVAILABLE",
            "Operator account boundary is unavailable.",
        )
    return uid, gid


class HostOperatorBackend:
    """Fixed production backend; tests inject all host-facing controllers."""

    def __init__(
        self,
        *,
        layout: DeploymentLayout | None = None,
        bundle_store: BundleStore | None = None,
        service_controller: ServiceController | None = None,
        observation_provider: ObservationProvider | None = None,
        deployment_controller: DeploymentActionController | None = None,
        verification_controller: DeploymentVerificationController | None = None,
        journal_store: OperatorJournalStore | None = None,
        state_store: StateStore | None = None,
        root_uid: int = 0,
        root_gid: int = 0,
        operator_group_gid: int | None = None,
        service_uid: int | None = None,
        service_gid: int | None = None,
        api_uid: int | None = None,
        api_gid: int | None = None,
        candidate_uid: int | None = None,
        candidate_gid: int | None = None,
        boundary_observer: OperatorBoundaryObserver | None = None,
    ) -> None:
        self.layout = DeploymentLayout.supported() if layout is None else layout
        if operator_group_gid is None:
            try:
                operator_group_gid = grp.getgrnam("helixweave-operators").gr_gid
            except KeyError:
                raise fail(
                    "OPERATOR_ACCOUNT_BOUNDARY_UNAVAILABLE",
                    "Operator account boundary is unavailable.",
                    recoverable=True,
                ) from None
        if not isinstance(operator_group_gid, int) or operator_group_gid <= 0:
            raise fail(
                "OPERATOR_ACCOUNT_BOUNDARY_UNAVAILABLE",
                "Operator account boundary is unavailable.",
            )
        self.root_uid = root_uid
        self.root_gid = root_gid
        self.operator_group_gid = operator_group_gid
        self.service_uid, self.service_gid = _resolve_account_identity(
            user="helixweave",
            group="helixweave",
            uid=service_uid,
            gid=service_gid,
        )
        self.api_uid, self.api_gid = _resolve_account_identity(
            user="helixweave-api",
            group="helixweave-api",
            uid=api_uid,
            gid=api_gid,
        )
        self.candidate_uid, self.candidate_gid = _resolve_account_identity(
            user="helixweave-candidate",
            group="helixweave-candidate",
            uid=candidate_uid,
            gid=candidate_gid,
        )
        if (
            len({self.service_uid, self.api_uid, self.candidate_uid}) != 3
            or len({self.service_gid, self.api_gid, self.candidate_gid}) != 3
        ):
            raise fail(
                "OPERATOR_ACCOUNT_BOUNDARY_UNAVAILABLE",
                "Operator account boundary is unavailable.",
            )
        self.boundary_observer = (
            HostOperatorBoundaryObserver(root_uid=root_uid, root_gid=root_gid)
            if boundary_observer is None
            else boundary_observer
        )
        self.bundle_store = (
            BundleStore(self.layout) if bundle_store is None else bundle_store
        )
        self.service_controller = (
            SystemdServiceController(
                self.layout, owner_uid=root_uid, owner_gid=root_gid
            )
            if service_controller is None
            else service_controller
        )
        self.observation_provider = (
            FixedObservationProvider(
                self.layout,
                self.service_controller,
                root_uid=root_uid,
                root_gid=root_gid,
                operator_group_gid=operator_group_gid,
                service_uid=self.service_uid,
                service_gid=self.service_gid,
            )
            if observation_provider is None
            else observation_provider
        )
        self.state_store = (
            StateStore(
                self.layout,
                reader_gid=operator_group_gid,
                service_gid=self.service_gid,
            )
            if state_store is None
            else state_store
        )
        self.deployment_controller = (
            HostDeploymentActionController(
                self.layout,
                states=self.state_store,
                services=self.service_controller,
                root_uid=root_uid,
                root_gid=root_gid,
                service_uid=self.service_uid,
                service_gid=self.service_gid,
                candidate_uid=self.candidate_uid,
                candidate_gid=self.candidate_gid,
            )
            if deployment_controller is None
            else deployment_controller
        )
        if verification_controller is None:
            if not isinstance(
                self.deployment_controller, HostDeploymentActionController
            ):
                self.verification_controller = None
            else:
                self.verification_controller = self.deployment_controller
        else:
            self.verification_controller = verification_controller
        self.journal_store = (
            OperatorJournalStore(
                self.layout,
                owner_uid=root_uid,
                owner_gid=root_gid,
                recovery_controller=(
                    self.deployment_controller
                    if isinstance(
                        self.deployment_controller, HostDeploymentActionController
                    )
                    else None
                ),
            )
            if journal_store is None
            else journal_store
        )

    def execute(
        self,
        request: OperatorRequest,
        *,
        bundle_path: Path | None,
    ) -> OperatorOutcome:
        if request.operation == STATUS:
            service = self.service_controller.status(request)
            return OperatorOutcome(
                "stopped" if service is None else "running", service=service
            )
        if request.operation == OBSERVE:
            observation = self.observation_provider.observe(request)
            return OperatorOutcome(
                "observed",
                observation=observation.with_operator_boundary(
                    self.boundary_observer.observe()
                ).with_operator_journal(self.journal_store.summary()),
            )
        if request.operation == VERIFY:
            if self.verification_controller is None:
                raise fail(
                    "OPERATOR_ACTION_UNAVAILABLE",
                    "Operator action is unavailable.",
                    recoverable=True,
                )
            verification = self.verification_controller.verify(request)
            boundary = self.boundary_observer.observe()
            return OperatorOutcome(
                "verified",
                verification=_with_boundary_readiness(verification, boundary),
            )
        with self.journal_store.operation(
            operation=request.operation,
            task_identity=request.task_identity,
            deployment_identity=request.deployment_identity,
            component=request.component,
            unit=request.unit,
        ) as journal:
            outcome = self._execute_mutation(
                request,
                bundle_path=bundle_path,
                journal=journal,
            )
            if not journal.terminal:
                journal.complete()
            return outcome

    def _execute_mutation(
        self,
        request: OperatorRequest,
        *,
        bundle_path: Path | None,
        journal: OperatorJournalHandle,
    ) -> OperatorOutcome:
        if request.operation == STAGE:
            if request.component is None or bundle_path is None:
                raise fail("OPERATOR_REQUEST_INVALID", "Operator request is invalid.")
            expected = self.layout.ingress_bundle(
                request.component, request.deployment_identity
            )
            if bundle_path != expected:
                raise fail("OPERATOR_REQUEST_INVALID", "Operator request is invalid.")
            group_gid = self._require_ingress_directory(request.component)
            journal.advance("ingress-pinned")
            manifest = self.bundle_store.stage(
                bundle_path,
                expected_owner_gid=group_gid,
                expected_mode=0o440,
                expected_parent_uid=self.root_uid,
                expected_parent_gid=group_gid,
                expected_parent_mode=0o2770,
                expected_component=request.component,
                expected_identity=request.deployment_identity,
                installed_owner_uid=self.root_uid,
                installed_owner_gid=self.root_gid,
            )
            if (
                manifest.component != request.component
                or manifest.identity != request.deployment_identity
            ):
                raise fail(
                    "DEPLOYMENT_BUNDLE_IDENTITY_MISMATCH",
                    "Deployment bundle does not match the requested identity.",
                )
            journal.advance("payload-staged")
            with self.state_store.transaction(
                exclusive=True,
                expected_owner_uid=self.root_uid,
                expected_owner_gid=self.root_gid,
            ) as transaction:
                state = transaction.initialize()
                pending = transaction.pending_transactions()
                if pending:
                    raise fail(
                        "DEPLOYMENT_RECOVERY_REQUIRED",
                        "An interrupted deployment operation requires recovery.",
                        recoverable=True,
                    )
                slots = state.components[request.component]
                if slots.staged != request.deployment_identity:
                    staged = state.stage(
                        request.component,
                        request.deployment_identity,
                    )
                    active = {
                        component: state.components[component].active
                        for component in COMPONENTS
                    }
                    journal.advance(
                        "candidate-selected",
                        prior_active=active,
                        candidate_active=active,
                        evidence={
                            "prior_state_identity": state.identity,
                            "candidate_state_identity": staged.identity,
                        },
                    )
                    transaction.commit(
                        staged,
                        operation=f"stage-{request.component}",
                        expected_current_identity=state.identity,
                    )
                    journal.advance(
                        "state-committed",
                        evidence={
                            "prior_state_identity": state.identity,
                            "candidate_state_identity": staged.identity,
                        },
                    )
            return OperatorOutcome("staged")
        if request.operation == START:
            assert request.unit is not None
            state_identity, active = self._service_mutation_state(request)
            journal.advance(
                (
                    "service-starting"
                    if request.unit in WRITER_UNITS
                    else "dependency-starting"
                ),
                point_of_no_return=request.unit in WRITER_UNITS,
                restart_units=(request.unit,),
                prior_active=active,
                candidate_active=active,
                evidence={
                    "prior_state_identity": state_identity,
                    "candidate_state_identity": state_identity,
                },
            )
            service = self.service_controller.start(request)
            journal.advance(
                "service-observed", evidence={"service_identity": service.identity}
            )
            return OperatorOutcome("running", service=service)
        if request.operation in {STOP, CLEANUP}:
            assert request.unit is not None and request.service_identity is not None
            state_identity, active = self._service_mutation_state(request)
            self._validate_bound_service_mutation(request)
            journal.advance(
                "service-stopping",
                restart_units=(request.unit,),
                prior_running_units=(request.unit,),
                prior_active=active,
                candidate_active=active,
                evidence={
                    "prior_state_identity": state_identity,
                    "candidate_state_identity": state_identity,
                    "service_identity": request.service_identity,
                },
            )
            self.service_controller.stop(request, cleanup=request.operation == CLEANUP)
            journal.advance("writers-stopped")
            return OperatorOutcome(
                "clean" if request.operation == CLEANUP else "stopped"
            )
        if request.operation in {ACTIVATE, ROLLBACK, UNINSTALL}:
            return self.deployment_controller.execute(request, journal=journal)
        raise fail("OPERATOR_REQUEST_INVALID", "Operator request is invalid.")

    def _service_mutation_state(
        self,
        request: OperatorRequest,
    ) -> tuple[str, dict[str, str | None]]:
        assert request.unit is not None
        with self.state_store.transaction(
            exclusive=False,
            expected_owner_uid=self.root_uid,
            expected_owner_gid=self.root_gid,
        ) as transaction:
            state = transaction.read()
            if transaction.pending_transactions():
                raise fail(
                    "DEPLOYMENT_RECOVERY_REQUIRED",
                    "An interrupted deployment operation requires recovery.",
                    recoverable=True,
                )
            component = (
                BULK_RNASEQ_RUNTIME
                if request.unit == "helixweave-docker-rootless.service"
                else PLATFORM
            )
            if state.components[component].active != request.deployment_identity:
                raise fail(
                    "OPERATOR_STATE_IDENTITY_MISMATCH",
                    "Operator state identity does not match the request.",
                    recoverable=True,
                )
            active = {name: state.components[name].active for name in COMPONENTS}
            return state.identity, active

    def _validate_bound_service_mutation(self, request: OperatorRequest) -> None:
        assert request.unit is not None and request.service_identity is not None
        observed = self.service_controller.status(
            OperatorRequest(
                operation=STATUS,
                task_identity=request.task_identity,
                deployment_identity=request.deployment_identity,
                unit=request.unit,
            )
        )
        if observed is None or observed.identity != request.service_identity:
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_MISMATCH",
                "Service identity changed before the requested action.",
                recoverable=True,
            )

    def _require_ingress_directory(self, component: str) -> int:
        path = self.layout.ingress / component
        group_gid = self.operator_group_gid
        try:
            boundaries = (
                (self.layout.data_root / "operator", group_gid, 0o710),
                (self.layout.ingress, group_gid, 0o750),
                (path, group_gid, 0o2770),
            )
            observed_boundaries = tuple(
                (boundary.lstat(), expected_gid, expected_mode)
                for boundary, expected_gid, expected_mode in boundaries
            )
        except OSError:
            raise fail(
                "OPERATOR_INGRESS_UNAVAILABLE", "Operator ingress is unavailable."
            ) from None
        for observed, expected_gid, expected_mode in observed_boundaries:
            if (
                not stat.S_ISDIR(observed.st_mode)
                or stat.S_ISLNK(observed.st_mode)
                or observed.st_uid != self.root_uid
                or observed.st_gid != expected_gid
                or stat.S_IMODE(observed.st_mode) != expected_mode
            ):
                raise fail(
                    "OPERATOR_INGRESS_UNTRUSTED", "Operator ingress is not trusted."
                )
        return group_gid


def _unix_protocol_readiness(
    path: Path,
    *,
    request: bytes,
    response_validator: Callable[[bytes], bool],
    unavailable_reason: str,
    scheme: str,
) -> ReadinessCheck:
    try:
        observed = path.lstat()
        if not stat.S_ISSOCK(observed.st_mode) or observed.st_nlink != 1:
            raise OSError
        response = bytearray()
        with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as client:
            client.settimeout(2.0)
            client.connect(str(path))
            client.sendall(request)
            while len(response) <= 4096:
                chunk = client.recv(min(1024, 4097 - len(response)))
                if not chunk:
                    break
                response.extend(chunk)
                if response_validator(bytes(response)):
                    break
        final = path.lstat()
        if (
            len(response) > 4096
            or not response_validator(bytes(response))
            or _file_witness(observed) != _file_witness(final)
        ):
            raise OSError
        evidence = canonical_identity(
            {
                "device": observed.st_dev,
                "inode": observed.st_ino,
                "mode": stat.S_IMODE(observed.st_mode),
                "protocol_response_sha256": hashlib.sha256(response).hexdigest(),
            },
            scheme=scheme,
        )
        return ReadinessCheck("ready", "READY", evidence)
    except OSError:
        return ReadinessCheck("unavailable", unavailable_reason)


def _reference_api_readiness(
    *, fetcher: Callable[[str], object] | None = None
) -> ReadinessCheck:
    """Consume only public, path-free reference registrations from the live API."""

    selected_fetcher = _fixed_api_json if fetcher is None else fetcher
    try:
        workflows = selected_fetcher("/api/v1/workflows/")
        if (
            not isinstance(workflows, dict)
            or set(workflows) != {"ok", "workflows", "issues"}
            or workflows["ok"] is not True
            or workflows["issues"] != []
            or not isinstance(workflows["workflows"], list)
            or not 0 < len(workflows["workflows"]) <= 32
        ):
            raise ValueError
        workflow_ids: list[str] = []
        for item in workflows["workflows"]:
            if not isinstance(item, dict):
                raise ValueError
            metadata = item.get("metadata")
            if not isinstance(metadata, dict):
                raise ValueError
            workflow_id = metadata.get("workflow_id")
            if (
                not isinstance(workflow_id, str)
                or workflow_id not in _REFERENCE_WORKFLOW_IDS
                or workflow_id in workflow_ids
            ):
                raise ValueError
            workflow_ids.append(workflow_id)
        if tuple(sorted(workflow_ids)) != _REFERENCE_WORKFLOW_IDS:
            raise ValueError
        registrations: list[dict[str, object]] = []
        profile_id = re.compile(r"^refp_[0-9a-f]{32}$")
        revision_id = re.compile(r"^refpr_[0-9a-f]{32}$")
        digest = re.compile(r"^[0-9a-f]{64}$")
        for workflow_id in _REFERENCE_WORKFLOW_IDS:
            document = selected_fetcher(
                f"/api/v1/workflows/{workflow_id}/reference-profiles"
            )
            if (
                not isinstance(document, dict)
                or set(document) != {"ok", "workflow_id", "profiles", "issues"}
                or document["ok"] is not True
                or document["workflow_id"] != workflow_id
                or document["issues"] != []
                or not isinstance(document["profiles"], list)
                or len(document["profiles"]) > 256
            ):
                raise ValueError
            for profile in document["profiles"]:
                if not isinstance(profile, dict) or set(profile) != {
                    "profile_id",
                    "revision_id",
                    "revision_number",
                    "display_name",
                    "organism",
                    "assembly",
                    "identity_sha256",
                }:
                    raise ValueError
                strings = tuple(
                    profile[name] for name in ("display_name", "organism", "assembly")
                )
                if (
                    not isinstance(profile["profile_id"], str)
                    or profile_id.fullmatch(profile["profile_id"]) is None
                    or not isinstance(profile["revision_id"], str)
                    or revision_id.fullmatch(profile["revision_id"]) is None
                    or not isinstance(profile["identity_sha256"], str)
                    or digest.fullmatch(profile["identity_sha256"]) is None
                    or not isinstance(profile["revision_number"], int)
                    or isinstance(profile["revision_number"], bool)
                    or profile["revision_number"] <= 0
                    or any(
                        not isinstance(value, str)
                        or not 0 < len(value.encode("utf-8")) <= 1024
                        or any(ord(character) < 32 for character in value)
                        for value in strings
                    )
                ):
                    raise ValueError
                registrations.append(
                    {
                        "workflow_id": workflow_id,
                        "profile_id": profile["profile_id"],
                        "revision_id": profile["revision_id"],
                        "revision_number": profile["revision_number"],
                        "identity_sha256": profile["identity_sha256"],
                    }
                )
        if not registrations:
            raise ValueError
        registrations.sort(
            key=lambda item: (
                str(item["workflow_id"]),
                str(item["profile_id"]),
                int(item["revision_number"]),
            )
        )
        return ReadinessCheck(
            "ready",
            "READY",
            canonical_identity(
                registrations,
                scheme="helixweave-reference-readiness-observation-v1",
            ),
        )
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError):
        return ReadinessCheck("not-ready", "REFERENCE_NOT_READY")


def _fixed_api_json(path: str) -> object:
    if path not in {
        "/api/v1/workflows/",
        *{
            f"/api/v1/workflows/{workflow_id}/reference-profiles"
            for workflow_id in _REFERENCE_WORKFLOW_IDS
        },
    }:
        raise ValueError
    request = (
        f"GET {path} HTTP/1.0\r\nHost: 127.0.0.1\r\nConnection: close\r\n\r\n"
    ).encode("ascii")
    response = bytearray()
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as client:
        client.settimeout(2.0)
        client.connect(("127.0.0.1", 8000))
        client.sendall(request)
        while len(response) <= _MAX_REFERENCE_API_BYTES:
            chunk = client.recv(
                min(65536, _MAX_REFERENCE_API_BYTES + 1 - len(response))
            )
            if not chunk:
                break
            response.extend(chunk)
    if len(response) > _MAX_REFERENCE_API_BYTES:
        raise ValueError
    raw_head, separator, body = bytes(response).partition(b"\r\n\r\n")
    lines = raw_head.split(b"\r\n")
    if (
        separator != b"\r\n\r\n"
        or not lines
        or lines[0]
        not in {
            b"HTTP/1.0 200 OK",
            b"HTTP/1.1 200 OK",
        }
    ):
        raise ValueError
    headers: dict[bytes, bytes] = {}
    for line in lines[1:]:
        name, marker, value = line.partition(b":")
        name = name.strip().lower()
        value = value.strip()
        if (
            marker != b":"
            or not name
            or name in headers
            or len(name) > 128
            or len(value) > 4096
        ):
            raise ValueError
        headers[name] = value
    length = headers.get(b"content-length")
    if (
        length is None
        or not length.isdigit()
        or int(length) != len(body)
        or not 0 < len(body) <= _MAX_REFERENCE_API_BYTES
        or b"transfer-encoding" in headers
    ):
        raise ValueError
    return json.loads(
        body,
        object_pairs_hook=_unique_object,
        parse_constant=lambda _value: (_ for _ in ()).throw(ValueError()),
    )


def _bytes_identity(value: bytes) -> str:
    return f"sha256-{hashlib.sha256(value).hexdigest()}"


def _write_all(descriptor: int, content: bytes) -> None:
    offset = 0
    while offset < len(content):
        written = os.write(descriptor, content[offset:])
        if written <= 0:
            raise OSError
        offset += written


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


def _read_bounded_fd(descriptor: int, limit: int) -> bytes:
    chunks: list[bytes] = []
    remaining = limit + 1
    while remaining > 0:
        chunk = os.read(descriptor, min(65536, remaining))
        if not chunk:
            break
        chunks.append(chunk)
        remaining -= len(chunk)
    content = b"".join(chunks)
    if len(content) > limit:
        raise OSError
    return content


def _read_bounded_path(path: Path, limit: int) -> bytes:
    descriptor = os.open(
        path,
        os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        observed = os.fstat(descriptor)
        if not stat.S_ISREG(observed.st_mode):
            raise OSError
        return _read_bounded_fd(descriptor, limit)
    finally:
        os.close(descriptor)


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


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError("duplicate key")
        value[key] = item
    return value


def parse_request(argv: Sequence[str]) -> OperatorRequest:
    """Parse the exact positional grammar without echoing rejected input."""
    if not isinstance(argv, Sequence) or isinstance(argv, (str, bytes)):
        raise fail("OPERATOR_REQUEST_INVALID", "Operator request is invalid.")
    values = tuple(argv)
    if any(not isinstance(item, str) or len(item) > 256 for item in values):
        raise fail("OPERATOR_REQUEST_INVALID", "Operator request is invalid.")
    if not values:
        raise fail("OPERATOR_REQUEST_INVALID", "Operator request is invalid.")
    operation = values[0]
    if operation in {STAGE, ACTIVATE, ROLLBACK} and len(values) == 4:
        return OperatorRequest(
            operation=operation,
            component=values[1],
            deployment_identity=values[2],
            task_identity=values[3],
        )
    if operation in {START, STATUS} and len(values) == 4:
        return OperatorRequest(
            operation=operation,
            unit=values[1],
            deployment_identity=values[2],
            task_identity=values[3],
        )
    if operation in BOUND_SERVICE_OPERATIONS and len(values) == 5:
        return OperatorRequest(
            operation=operation,
            unit=values[1],
            deployment_identity=values[2],
            task_identity=values[3],
            service_identity=values[4],
        )
    if operation == UNINSTALL and len(values) == 3:
        return OperatorRequest(
            operation=operation,
            deployment_identity=values[1],
            task_identity=values[2],
        )
    if operation in {OBSERVE, VERIFY} and len(values) == 3:
        return OperatorRequest(
            operation=operation,
            deployment_identity=values[1],
            task_identity=values[2],
        )
    raise fail("OPERATOR_REQUEST_INVALID", "Operator request is invalid.")


def bundle_ingress_path(layout: DeploymentLayout, request: OperatorRequest) -> Path:
    """Derive the only admitted bundle path; no path is accepted from argv."""
    if request.operation != STAGE or request.component is None:
        raise fail("OPERATOR_REQUEST_INVALID", "Operator request is invalid.")
    return layout.ingress_bundle(request.component, request.deployment_identity)


def execute_request(
    argv: Sequence[str],
    *,
    backend: OperatorBackend,
    layout: DeploymentLayout | None = None,
) -> dict[str, object]:
    """Validate, execute through the injected backend, and bind the receipt."""
    request = parse_request(argv)
    selected_layout = DeploymentLayout.supported() if layout is None else layout
    bundle_path = (
        bundle_ingress_path(selected_layout, request)
        if request.operation == STAGE
        else None
    )
    try:
        outcome = backend.execute(request, bundle_path=bundle_path)
    except DeploymentError as error:
        raise fail(
            "OPERATOR_OPERATION_FAILED",
            "Operator action failed.",
            component=request.component,
            recoverable=error.issue.recoverable,
        ) from None
    except Exception:
        raise fail(
            "OPERATOR_OPERATION_FAILED",
            "Operator action failed.",
            component=request.component,
        ) from None

    if not isinstance(outcome, OperatorOutcome):
        raise fail("OPERATOR_RESULT_INVALID", "Operator result is invalid.")
    _verify_outcome(request, outcome)
    receipt: dict[str, object] = {
        "schema_version": OPERATOR_RECEIPT_SCHEMA,
        "operation": request.operation,
        "state": outcome.state,
        "task_identity": request.task_identity,
        "deployment_identity": request.deployment_identity,
    }
    if request.component is not None:
        receipt["component"] = request.component
    if request.unit is not None:
        receipt["unit"] = request.unit
    if request.service_identity is not None:
        receipt["prior_service_identity"] = request.service_identity
    if outcome.service is not None:
        receipt["service_identity"] = outcome.service.to_dict()
    if outcome.observation is not None:
        receipt["observation"] = outcome.observation.to_dict()
    if outcome.verification is not None:
        receipt["verification"] = outcome.verification.to_dict()
    return receipt


def _verify_outcome(request: OperatorRequest, outcome: OperatorOutcome) -> None:
    expected_states = {
        STAGE: {"staged"},
        ACTIVATE: {"activated"},
        ROLLBACK: {"rolled-back"},
        START: {"running"},
        STATUS: {"running", "stopped"},
        STOP: {"stopped"},
        CLEANUP: {"clean"},
        UNINSTALL: {"uninstalled"},
        OBSERVE: {"observed"},
        VERIFY: {"verified"},
    }
    if outcome.state not in expected_states[request.operation]:
        raise fail("OPERATOR_RESULT_INVALID", "Operator result is invalid.")
    service_required = request.operation == START or (
        request.operation == STATUS and outcome.state == "running"
    )
    if service_required != (outcome.service is not None):
        raise fail("OPERATOR_RESULT_INVALID", "Operator result is invalid.")
    observation_required = request.operation == OBSERVE
    if observation_required != (outcome.observation is not None):
        raise fail("OPERATOR_RESULT_INVALID", "Operator result is invalid.")
    if outcome.service is not None and outcome.observation is not None:
        raise fail("OPERATOR_RESULT_INVALID", "Operator result is invalid.")
    verification_required = request.operation == VERIFY
    if verification_required != (outcome.verification is not None):
        raise fail("OPERATOR_RESULT_INVALID", "Operator result is invalid.")
    if outcome.verification is not None and (
        outcome.service is not None or outcome.observation is not None
    ):
        raise fail("OPERATOR_RESULT_INVALID", "Operator result is invalid.")
    if (
        outcome.observation is not None
        and outcome.observation.state_identity != request.deployment_identity
    ):
        raise fail("OPERATOR_RESULT_INVALID", "Operator result is invalid.")
    if outcome.service is not None:
        if (
            outcome.service.unit != request.unit
            or outcome.service.deployment_identity != request.deployment_identity
            or (
                request.operation == START
                and outcome.service.task_identity != request.task_identity
            )
            or (
                request.service_identity is not None
                and outcome.service.identity != request.service_identity
            )
        ):
            raise fail("OPERATOR_RESULT_INVALID", "Operator result is invalid.")


def _error_receipt(issue: DeploymentIssue) -> dict[str, object]:
    return {
        "schema_version": OPERATOR_RECEIPT_SCHEMA,
        "status": "error",
        "issue": issue.to_dict(),
    }


def _write_json(stream, value: object) -> None:
    stream.write(json.dumps(value, sort_keys=True, separators=(",", ":")))
    stream.write("\n")
    stream.flush()


def main(
    argv: Sequence[str] | None = None,
    *,
    backend: OperatorBackend | None = None,
) -> int:
    """Run the sudo-facing helper without accepting caller environment values."""
    arguments = tuple(sys.argv[1:] if argv is None else argv)
    if os.geteuid() != 0:
        issue = DeploymentIssue(
            code="OPERATOR_PRIVILEGE_REQUIRED",
            message="Operator helper must run as root.",
        )
        _write_json(sys.stderr, _error_receipt(issue))
        return 77

    os.environ.clear()
    os.environ.update(
        {
            "LANG": "C.UTF-8",
            "LC_ALL": "C.UTF-8",
            "PATH": "/usr/sbin:/usr/bin:/sbin:/bin",
        }
    )
    try:
        receipt = execute_request(
            arguments,
            backend=HostOperatorBackend() if backend is None else backend,
        )
    except DeploymentError as error:
        _write_json(sys.stderr, _error_receipt(error.issue))
        if error.issue.code == "OPERATOR_REQUEST_INVALID":
            return 64
        if error.issue.code == "OPERATOR_OPERATION_FAILED":
            return 69 if error.issue.recoverable else 70
        return 70
    except Exception:
        issue = DeploymentIssue(
            code="OPERATOR_OPERATION_FAILED",
            message="Operator action failed.",
        )
        _write_json(sys.stderr, _error_receipt(issue))
        return 70
    _write_json(sys.stdout, receipt)
    return 0


__all__ = [
    "ACTIVATE",
    "BOUND_SERVICE_OPERATIONS",
    "CLEANUP",
    "BoundedCommandExecutor",
    "CommandExecutor",
    "CommandResult",
    "DeploymentActionController",
    "BoundaryUninstaller",
    "FixedSystemctl",
    "HostOperatorBackend",
    "HostBoundaryUninstaller",
    "HostDeploymentActionController",
    "HostOperatorBoundaryObserver",
    "LinuxServiceProbe",
    "OPERATIONS",
    "OPERATOR_RECEIPT_SCHEMA",
    "OBSERVE",
    "OperatorBackend",
    "OperatorBoundaryObservation",
    "OperatorBoundaryObserver",
    "OperatorObservation",
    "OperatorOutcome",
    "OperatorRequest",
    "ROLLBACK",
    "SERVICE_OPERATIONS",
    "SERVICE_SOCKET_NAMES",
    "SERVICE_UNITS",
    "ServiceController",
    "ServiceProbe",
    "START",
    "STATUS",
    "STOP",
    "STAGE",
    "ServiceIdentity",
    "SocketWitness",
    "UNINSTALL",
    "UNINSTALL_BOUNDARY_FILES",
    "UNINSTALL_LINKED_BOUNDARY_TARGETS",
    "VERIFY",
    "UnavailableOperatorBackend",
    "UnavailableDeploymentActionController",
    "SystemdServiceController",
    "bundle_ingress_path",
    "execute_request",
    "main",
    "parse_request",
]
