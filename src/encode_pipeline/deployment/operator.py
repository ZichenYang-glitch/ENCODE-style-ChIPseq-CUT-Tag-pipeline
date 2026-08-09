"""Fail-closed request boundary for the privileged deployment operator.

The sudo-facing process accepts only the small grammar defined here.  Host
mutation is deliberately delegated to an injected backend so request parsing,
path derivation, and public receipts can be tested without systemd or root.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
import json
import os
from pathlib import Path
import re
import sys
from typing import Protocol

from encode_pipeline.deployment.canonical import canonical_identity
from encode_pipeline.deployment.errors import DeploymentError, DeploymentIssue, fail
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import COMPONENTS


OPERATOR_RECEIPT_SCHEMA = "helixweave-operator-receipt-v1"
SERVICE_IDENTITY_SCHEME = "helixweave-service-identity-v1"

STAGE = "stage"
ACTIVATE = "activate"
ROLLBACK = "rollback"
START = "start"
STATUS = "status"
STOP = "stop"
CLEANUP = "cleanup"
UNINSTALL = "uninstall"

OPERATIONS = (
    STAGE,
    ACTIVATE,
    ROLLBACK,
    START,
    STATUS,
    STOP,
    CLEANUP,
    UNINSTALL,
)
SERVICE_OPERATIONS = (START, STATUS, STOP, CLEANUP)
BOUND_SERVICE_OPERATIONS = (STATUS, STOP, CLEANUP)

SERVICE_UNITS = (
    "helixweave-api.service",
    "helixweave-worker.service",
    "helixweave-redis.service",
    "helixweave-docker-rootless.service",
)

SERVICE_SOCKET_NAMES: dict[str, tuple[str, ...]] = {
    "helixweave-api.service": ("api-http",),
    "helixweave-worker.service": (),
    "helixweave-redis.service": ("redis-queue",),
    "helixweave-docker-rootless.service": ("bulk-docker",),
}

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

    def __post_init__(self) -> None:
        if (
            self.name
            not in {name for names in SERVICE_SOCKET_NAMES.values() for name in names}
            or not _positive_integer(self.device)
            or not _positive_integer(self.inode)
        ):
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_INVALID",
                "Service identity evidence is invalid.",
            )

    def to_dict(self) -> dict[str, object]:
        return {"name": self.name, "device": self.device, "inode": self.inode}


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
                )
                for item in raw_sockets
                if isinstance(item, dict) and set(item) == {"name", "device", "inode"}
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
class OperatorOutcome:
    """Whitelisted backend result; it cannot carry paths or exception strings."""

    state: str
    service: ServiceIdentity | None = None

    def __post_init__(self) -> None:
        if self.state not in {
            "staged",
            "activated",
            "rolled-back",
            "running",
            "stopped",
            "clean",
            "uninstalled",
        } or (
            self.service is not None and not isinstance(self.service, ServiceIdentity)
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
    if operation == START and len(values) == 4:
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
    raise fail("OPERATOR_REQUEST_INVALID", "Operator request is invalid.")


def bundle_ingress_path(layout: DeploymentLayout, request: OperatorRequest) -> Path:
    """Derive the only admitted bundle path; no path is accepted from argv."""
    if request.operation != STAGE or request.component is None:
        raise fail("OPERATOR_REQUEST_INVALID", "Operator request is invalid.")
    return (
        layout.ingress / request.component / request.deployment_identity / "bundle.tar"
    )


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
    }
    if outcome.state not in expected_states[request.operation]:
        raise fail("OPERATOR_RESULT_INVALID", "Operator result is invalid.")
    service_required = request.operation in {START, STATUS}
    if service_required != (outcome.service is not None):
        raise fail("OPERATOR_RESULT_INVALID", "Operator result is invalid.")
    if outcome.service is not None:
        if (
            outcome.service.unit != request.unit
            or outcome.service.deployment_identity != request.deployment_identity
            or outcome.service.task_identity != request.task_identity
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
            backend=UnavailableOperatorBackend() if backend is None else backend,
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
    "OPERATIONS",
    "OPERATOR_RECEIPT_SCHEMA",
    "OperatorBackend",
    "OperatorOutcome",
    "OperatorRequest",
    "ROLLBACK",
    "SERVICE_OPERATIONS",
    "SERVICE_SOCKET_NAMES",
    "SERVICE_UNITS",
    "START",
    "STATUS",
    "STOP",
    "STAGE",
    "ServiceIdentity",
    "SocketWitness",
    "UNINSTALL",
    "UnavailableOperatorBackend",
    "bundle_ingress_path",
    "execute_request",
    "main",
    "parse_request",
]
