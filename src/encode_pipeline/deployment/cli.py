"""Scriptable JSON command contract for supported deployment operations.

The production default is composed lazily so this contract module remains
independent from the fixed operator client.  Tests may inject a backend, but a
normal invocation always selects the supported host composition.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import json
from pathlib import Path
import re
import sys
from typing import Protocol, TextIO

from encode_pipeline.deployment.errors import (
    DeploymentError,
    DeploymentIssue,
    fail,
)
from encode_pipeline.deployment.doctor import (
    CHECKS as DOCTOR_CHECK_INVENTORY,
    DOCTOR_REPORT_SCHEMA,
    PUBLIC_REASON_CODES,
)
from encode_pipeline.deployment.models import COMPONENTS


CLI_RESULT_SCHEMA = "helixweave-deployment-cli-result-v1"

INSTALL = "install"
STATUS = "status"
DOCTOR = "doctor"
VERIFY = "verify"
UPGRADE = "upgrade"
ROLLBACK = "rollback"
COMMANDS = (INSTALL, STATUS, DOCTOR, VERIFY, UPGRADE, ROLLBACK)

EXIT_OK = 0
EXIT_DATA = 65
EXIT_UNAVAILABLE = 69
EXIT_OPERATION = 70
EXIT_PERMISSION = 77
EXIT_INCOMPATIBLE = 78
EXIT_USAGE = 64
_BACKEND_EXIT_CODES = {EXIT_OK, EXIT_DATA, EXIT_UNAVAILABLE, EXIT_INCOMPATIBLE}

OPERATION_RESULT_SCHEMA = "helixweave-deployment-operation-v1"
STATUS_RESULT_SCHEMA = "helixweave-deployment-status-v2"
VERIFY_RESULT_SCHEMA = "helixweave-deployment-verify-v1"
DOCTOR_RESULT_SCHEMA = DOCTOR_REPORT_SCHEMA

_IDENTITY = re.compile(r"^sha256-[0-9a-f]{64}$")
_PUBLIC_KEY = re.compile(r"^[a-z][a-z0-9_-]{0,63}$")
_PUBLIC_TEXT = re.compile(r"^[0-9A-Za-z][0-9A-Za-z ._:+-]{0,255}$")
_FORBIDDEN_KEYS = (
    "path",
    "secret",
    "password",
    "token",
    "environment",
    "exception",
    "traceback",
    "cmdline",
    "database_url",
    "reference_coordinate",
)
_MUTATION_STATES = {
    INSTALL: "staged",
    UPGRADE: "activated",
    ROLLBACK: "rolled-back",
}
_SERVICE_KEYS = ("api", "worker", "redis", "bulk_rnaseq_docker")
_DOCTOR_CHECKS = tuple(check_id for check_id, _category in DOCTOR_CHECK_INVENTORY)
_DOCTOR_CATEGORIES = dict(DOCTOR_CHECK_INVENTORY)
_PUBLIC_ISSUE_EXIT_CODES = {
    "DEPLOYMENT_COMMAND_INVALID": EXIT_USAGE,
    "DEPLOYMENT_BUNDLE_INVALID": EXIT_DATA,
    "DEPLOYMENT_MANIFEST_INVALID": EXIT_DATA,
    "DEPLOYMENT_RELEASE_INVALID": EXIT_DATA,
    "DEPLOYMENT_RESULT_INVALID": EXIT_OPERATION,
    "DEPLOYMENT_BUNDLE_UNAVAILABLE": EXIT_UNAVAILABLE,
    "DEPLOYMENT_RELEASE_UNAVAILABLE": EXIT_UNAVAILABLE,
    "DEPLOYMENT_STATE_UNAVAILABLE": EXIT_UNAVAILABLE,
    "DEPLOYMENT_OPERATOR_UNAVAILABLE": EXIT_UNAVAILABLE,
    "OPERATOR_OBSERVATION_UNAVAILABLE": EXIT_UNAVAILABLE,
    "DEPLOYMENT_INGRESS_UNAVAILABLE": EXIT_UNAVAILABLE,
    "DEPLOYMENT_INTEGRATION_DEFERRED": EXIT_UNAVAILABLE,
    "DEPLOYMENT_CONTRACT_ADMISSION_DEFERRED": EXIT_UNAVAILABLE,
    "DEPLOYMENT_SCHEMA_OBSERVATION_DEFERRED": EXIT_UNAVAILABLE,
    "DEPLOYMENT_BUSY": EXIT_UNAVAILABLE,
    "DEPLOYMENT_RECOVERY_REQUIRED": EXIT_UNAVAILABLE,
    "DEPLOYMENT_CAPACITY_INSUFFICIENT": EXIT_UNAVAILABLE,
    "DEPLOYMENT_COMPATIBILITY_FAILED": EXIT_INCOMPATIBLE,
    "DEPLOYMENT_SCHEMA_INCOMPATIBLE": EXIT_INCOMPATIBLE,
    "DEPLOYMENT_STAGED_IDENTITY_MISMATCH": EXIT_INCOMPATIBLE,
    "DEPLOYMENT_PREVIOUS_IDENTITY_MISMATCH": EXIT_INCOMPATIBLE,
    "DEPLOYMENT_BUNDLE_COMPONENT_MISMATCH": EXIT_INCOMPATIBLE,
    "FRONTEND_API_INCOMPATIBLE": EXIT_INCOMPATIBLE,
    "DEPLOYMENT_BUNDLE_SOURCE_INVALID": EXIT_DATA,
    "DEPLOYMENT_BUNDLE_SOURCE_CHANGED": EXIT_DATA,
    "DEPLOYMENT_BUNDLE_OUTPUT_INVALID": EXIT_DATA,
    "DEPLOYMENT_BUNDLE_OUTPUT_EXISTS": EXIT_DATA,
    "DEPLOYMENT_BUNDLE_BUILD_FAILED": EXIT_OPERATION,
    "DEPLOYMENT_COMPONENT_INVALID": EXIT_USAGE,
    "DEPLOYMENT_TASK_IDENTITY_INVALID": EXIT_OPERATION,
    "DEPLOYMENT_PRIVILEGE_REQUIRED": EXIT_PERMISSION,
    "DEPLOYMENT_PERMISSION_DENIED": EXIT_PERMISSION,
}


@dataclass(frozen=True)
class CommandRequest:
    command: str
    component: str | None = None
    bundle_path: Path | None = None
    deployment_identity: str | None = None

    def __post_init__(self) -> None:
        bundle_command = self.command in {INSTALL, UPGRADE}
        component_command = bundle_command or self.command == ROLLBACK
        if (
            self.command not in COMMANDS
            or component_command != (self.component is not None)
            or bundle_command != (self.bundle_path is not None)
            or (self.command == ROLLBACK) != (self.deployment_identity is not None)
            or (self.component is not None and self.component not in COMPONENTS)
            or (
                self.deployment_identity is not None
                and _IDENTITY.fullmatch(self.deployment_identity) is None
            )
        ):
            raise fail("DEPLOYMENT_COMMAND_INVALID", "Deployment command is invalid.")


@dataclass(frozen=True)
class PublicCommandResult:
    """Backend result whose values are checked before they reach stdout."""

    value: Mapping[str, object]
    exit_code: int = EXIT_OK

    def __post_init__(self) -> None:
        if self.exit_code not in _BACKEND_EXIT_CODES:
            raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
        _validate_result_document(self.value)


class CommandBackend(Protocol):
    def install(self, component: str, bundle_path: Path) -> PublicCommandResult: ...

    def status(self) -> PublicCommandResult: ...

    def doctor(self) -> PublicCommandResult: ...

    def verify(self) -> PublicCommandResult: ...

    def upgrade(self, component: str, bundle_path: Path) -> PublicCommandResult: ...

    def rollback(
        self, component: str, deployment_identity: str
    ) -> PublicCommandResult: ...


class UnavailableCommandBackend:
    """Explicit fail-closed backend retained for contract-level callers."""

    @staticmethod
    def _unavailable() -> PublicCommandResult:
        raise fail(
            "DEPLOYMENT_INTEGRATION_DEFERRED",
            "Deployment service is not ready.",
            recoverable=True,
        )

    def install(self, component: str, bundle_path: Path) -> PublicCommandResult:
        del component, bundle_path
        return self._unavailable()

    def status(self) -> PublicCommandResult:
        return self._unavailable()

    def doctor(self) -> PublicCommandResult:
        return self._unavailable()

    def verify(self) -> PublicCommandResult:
        return self._unavailable()

    def upgrade(self, component: str, bundle_path: Path) -> PublicCommandResult:
        del component, bundle_path
        return self._unavailable()

    def rollback(self, component: str, deployment_identity: str) -> PublicCommandResult:
        del component, deployment_identity
        return self._unavailable()


def parse_command(argv: Sequence[str]) -> CommandRequest:
    """Parse an exact flag order without echoing rejected arguments."""
    if not isinstance(argv, Sequence) or isinstance(argv, (str, bytes)):
        raise fail("DEPLOYMENT_COMMAND_INVALID", "Deployment command is invalid.")
    values = tuple(argv)
    if any(
        not isinstance(item, str)
        or not item
        or len(item) > 4096
        or any(ord(character) < 32 or ord(character) == 127 for character in item)
        for item in values
    ):
        raise fail("DEPLOYMENT_COMMAND_INVALID", "Deployment command is invalid.")
    if len(values) == 1 and values[0] in {STATUS, DOCTOR, VERIFY}:
        return CommandRequest(command=values[0])
    if (
        len(values) == 5
        and values[0] in {INSTALL, UPGRADE}
        and values[1] == "--component"
        and values[3] == "--bundle"
    ):
        return CommandRequest(
            command=values[0],
            component=values[2],
            bundle_path=_bundle_path(values[4]),
        )
    if (
        len(values) == 5
        and values[0] == ROLLBACK
        and values[1] == "--component"
        and values[3] == "--identity"
    ):
        return CommandRequest(
            command=ROLLBACK,
            component=values[2],
            deployment_identity=values[4],
        )
    raise fail("DEPLOYMENT_COMMAND_INVALID", "Deployment command is invalid.")


def execute_command(
    request: CommandRequest,
    *,
    backend: CommandBackend,
) -> PublicCommandResult:
    if request.command == INSTALL:
        assert request.component is not None and request.bundle_path is not None
        result = backend.install(request.component, request.bundle_path)
    elif request.command == STATUS:
        result = backend.status()
    elif request.command == DOCTOR:
        result = backend.doctor()
    elif request.command == VERIFY:
        result = backend.verify()
    elif request.command == UPGRADE:
        assert request.component is not None and request.bundle_path is not None
        result = backend.upgrade(request.component, request.bundle_path)
    else:
        assert (
            request.command == ROLLBACK
            and request.component is not None
            and request.deployment_identity is not None
        )
        result = backend.rollback(request.component, request.deployment_identity)
    if not isinstance(result, PublicCommandResult):
        raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
    _validate_command_result(request, result)
    return result


def success_envelope(
    request: CommandRequest,
    result: PublicCommandResult,
) -> dict[str, object]:
    return {
        "schema_version": CLI_RESULT_SCHEMA,
        "command": request.command,
        "status": "ok" if result.exit_code == EXIT_OK else "not-ready",
        "result": _public_copy(result.value),
    }


def _bundle_path(value: str) -> Path:
    path = Path(value)
    if (
        not path.is_absolute()
        or ".." in path.parts
        or "\x00" in value
        or "\n" in value
        or "\r" in value
    ):
        raise fail("DEPLOYMENT_COMMAND_INVALID", "Deployment command is invalid.")
    return path


def _public_copy(value: object, *, depth: int = 0) -> object:
    if depth > 16:
        raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
    if value is None or isinstance(value, bool):
        return value
    if isinstance(value, int) and not isinstance(value, bool):
        if -(2**63) <= value < 2**63:
            return value
        raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
    if isinstance(value, str):
        if _PUBLIC_TEXT.fullmatch(value) is None or "/" in value or "\\" in value:
            raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
        return value
    if isinstance(value, Mapping):
        copied: dict[str, object] = {}
        if len(value) > 256:
            raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
        for key, item in value.items():
            if (
                not isinstance(key, str)
                or _PUBLIC_KEY.fullmatch(key) is None
                or any(blocked in key for blocked in _FORBIDDEN_KEYS)
                or key in copied
            ):
                raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
            copied[key] = _public_copy(item, depth=depth + 1)
        return copied
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes)):
        if len(value) > 10_000:
            raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
        return [_public_copy(item, depth=depth + 1) for item in value]
    raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")


def _validate_result_document(value: Mapping[str, object]) -> None:
    if not isinstance(value, Mapping):
        raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
    schema = value.get("schema_version")
    if schema == OPERATION_RESULT_SCHEMA:
        _validate_operation_result(value)
    elif schema == STATUS_RESULT_SCHEMA:
        _validate_status_result(value)
    elif schema == DOCTOR_RESULT_SCHEMA:
        _validate_doctor_result(value)
    elif schema == VERIFY_RESULT_SCHEMA:
        _validate_verify_result(value)
    else:
        raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
    _public_copy(value)


def _exact_mapping(value: object, keys: set[str]) -> Mapping[str, object]:
    if (
        not isinstance(value, Mapping)
        or any(not isinstance(key, str) for key in value)
        or set(value) != keys
    ):
        raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
    return value


def _identity_or_none(value: object) -> bool:
    return value is None or (
        isinstance(value, str) and _IDENTITY.fullmatch(value) is not None
    )


def _validate_operation_result(value: Mapping[str, object]) -> None:
    document = _exact_mapping(
        value,
        {
            "schema_version",
            "operation",
            "state",
            "component",
            "deployment_identity",
            "state_identity",
        },
    )
    operation = document["operation"]
    if (
        operation not in _MUTATION_STATES
        or document["state"] != _MUTATION_STATES[operation]
        or document["component"] not in COMPONENTS
        or not _identity_or_none(document["deployment_identity"])
        or document["deployment_identity"] is None
        or not _identity_or_none(document["state_identity"])
        or document["state_identity"] is None
    ):
        raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")


def _validate_slot(value: object) -> None:
    if value is None:
        return
    slot = _exact_mapping(value, {"identity"})
    if not _identity_or_none(slot["identity"]) or slot["identity"] is None:
        raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")


def _validate_status_result(value: Mapping[str, object]) -> None:
    document = _exact_mapping(
        value,
        {
            "schema_version",
            "state_identity",
            "generation",
            "components",
            "interrupted",
            "partial_staging_count",
            "pending_transaction_count",
            "orphaned_deployment_count",
            "operator_pending_transaction_count",
            "operator_recovery_required_count",
            "database_schema_identity",
            "database_schema_reason_code",
            "services",
        },
    )
    if (
        not _identity_or_none(document["state_identity"])
        or document["state_identity"] is None
        or not isinstance(document["generation"], int)
        or isinstance(document["generation"], bool)
        or document["generation"] < 0
        or not isinstance(document["interrupted"], bool)
        or not _identity_or_none(document["database_schema_identity"])
        or document["database_schema_reason_code"]
        not in {"DATABASE_READY", "DATABASE_UNAVAILABLE"}
        or (
            (document["database_schema_identity"] is not None)
            != (document["database_schema_reason_code"] == "DATABASE_READY")
        )
    ):
        raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
    for key in (
        "partial_staging_count",
        "pending_transaction_count",
        "orphaned_deployment_count",
        "operator_pending_transaction_count",
        "operator_recovery_required_count",
    ):
        count = document[key]
        if not isinstance(count, int) or isinstance(count, bool) or count < 0:
            raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
    if (
        document["operator_pending_transaction_count"] not in {0, 1}
        or document["operator_recovery_required_count"] not in {0, 1}
        or document["operator_pending_transaction_count"]
        + document["operator_recovery_required_count"]
        > 1
    ):
        raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
    components = _exact_mapping(document["components"], set(COMPONENTS))
    for component in COMPONENTS:
        slots = _exact_mapping(components[component], {"active", "previous", "staged"})
        for slot in ("active", "previous", "staged"):
            _validate_slot(slots[slot])
    services = _exact_mapping(document["services"], set(_SERVICE_KEYS))
    for name in _SERVICE_KEYS:
        service = _exact_mapping(services[name], {"state", "identity"})
        if service["state"] not in {
            "running",
            "stopped",
            "unavailable",
        } or not _identity_or_none(service["identity"]):
            raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
        if (service["state"] == "running") != (service["identity"] is not None):
            raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")


def _validate_doctor_result(value: Mapping[str, object]) -> None:
    document = _exact_mapping(
        value,
        {"schema_version", "status", "ready", "checks"},
    )
    if (
        document["status"] not in {"healthy", "degraded", "unhealthy"}
        or not isinstance(document["ready"], bool)
        or (document["status"] == "healthy") != document["ready"]
        or not isinstance(document["checks"], Sequence)
        or isinstance(document["checks"], (str, bytes))
        or len(document["checks"]) != len(_DOCTOR_CHECKS)
    ):
        raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
    observed_ids: list[str] = []
    for item in document["checks"]:
        check = _exact_mapping(
            item,
            {"check_id", "category", "state", "reason_code"}
            | (
                {"evidence_identity"}
                if isinstance(item, Mapping) and "evidence_identity" in item
                else set()
            ),
        )
        if (
            check["check_id"] not in _DOCTOR_CHECKS
            or check["category"] != _DOCTOR_CATEGORIES[check["check_id"]]
            or check["state"] not in {"pass", "warning", "fail"}
            or not isinstance(check["reason_code"], str)
            or check["reason_code"] not in PUBLIC_REASON_CODES
            or (
                "evidence_identity" in check
                and not _identity_or_none(check["evidence_identity"])
            )
        ):
            raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
        observed_ids.append(check["check_id"])
    if tuple(observed_ids) != _DOCTOR_CHECKS:
        raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")


def _validate_verify_result(value: Mapping[str, object]) -> None:
    document = _exact_mapping(
        value,
        {
            "schema_version",
            "verified",
            "deployment",
            "frontend_identity",
            "database_schema_identity",
        },
    )
    if (
        not isinstance(document["verified"], bool)
        or not _identity_or_none(document["frontend_identity"])
        or not _identity_or_none(document["database_schema_identity"])
    ):
        raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
    _validate_status_result(document["deployment"])


def _validate_command_result(
    request: CommandRequest,
    result: PublicCommandResult,
) -> None:
    schema = result.value["schema_version"]
    expected = {
        INSTALL: OPERATION_RESULT_SCHEMA,
        UPGRADE: OPERATION_RESULT_SCHEMA,
        ROLLBACK: OPERATION_RESULT_SCHEMA,
        STATUS: STATUS_RESULT_SCHEMA,
        DOCTOR: DOCTOR_RESULT_SCHEMA,
        VERIFY: VERIFY_RESULT_SCHEMA,
    }[request.command]
    if schema != expected:
        raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
    if request.command in _MUTATION_STATES or request.command == STATUS:
        if result.exit_code != EXIT_OK:
            raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
    elif request.command == DOCTOR:
        expected_exit = EXIT_OK if result.value["ready"] else EXIT_UNAVAILABLE
        if result.exit_code != expected_exit:
            raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
    elif request.command == VERIFY:
        expected_exit = EXIT_OK if result.value["verified"] else EXIT_DATA
        if result.exit_code != expected_exit:
            raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")
    if request.command in _MUTATION_STATES:
        if (
            result.value["operation"] != request.command
            or result.value["component"] != request.component
            or (
                request.command == ROLLBACK
                and result.value["deployment_identity"] != request.deployment_identity
            )
        ):
            raise fail("DEPLOYMENT_RESULT_INVALID", "Deployment result is invalid.")


def _exit_code(issue: DeploymentIssue) -> int:
    return _PUBLIC_ISSUE_EXIT_CODES.get(issue.code, EXIT_OPERATION)


def _public_issue(issue: DeploymentIssue, exit_code: int) -> DeploymentIssue:
    code = issue.code if issue.code in _PUBLIC_ISSUE_EXIT_CODES else "DEPLOYMENT_FAILED"
    messages = {
        EXIT_USAGE: "Deployment command is invalid.",
        EXIT_DATA: "Deployment verification failed.",
        EXIT_UNAVAILABLE: "Deployment service is not ready.",
        EXIT_PERMISSION: "Deployment permission was denied.",
        EXIT_INCOMPATIBLE: "Deployment is not compatible.",
        EXIT_OPERATION: "Deployment operation failed.",
    }
    component = issue.component if issue.component in COMPONENTS else None
    return DeploymentIssue(
        code=code,
        message=messages[exit_code],
        component=component,
        recoverable=issue.recoverable,
    )


def _write_json(stream: TextIO, value: object) -> None:
    stream.write(
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        )
    )
    stream.write("\n")
    stream.flush()


def main(
    argv: Sequence[str] | None = None,
    *,
    backend: CommandBackend | None = None,
    stdout: TextIO | None = None,
    stderr: TextIO | None = None,
) -> int:
    """Execute one command and emit exactly one canonical JSON document."""
    arguments = tuple(sys.argv[1:] if argv is None else argv)
    output = sys.stdout if stdout is None else stdout
    errors = sys.stderr if stderr is None else stderr
    request: CommandRequest | None = None
    try:
        request = parse_command(arguments)
        selected_backend = backend
        if selected_backend is None:
            # The lazy import avoids a module cycle: the production backend
            # consumes the public result contract defined above.
            from encode_pipeline.deployment.backend import ProductionCommandBackend

            selected_backend = ProductionCommandBackend.supported()
        result = execute_command(
            request,
            backend=selected_backend,
        )
        envelope = success_envelope(request, result)
    except DeploymentError as error:
        exit_code = _exit_code(error.issue)
        public_issue = _public_issue(error.issue, exit_code)
        _write_json(
            errors,
            {
                "schema_version": CLI_RESULT_SCHEMA,
                "command": None if request is None else request.command,
                "status": "error",
                "issue": public_issue.to_dict(),
            },
        )
        return exit_code
    except Exception:
        _write_json(
            errors,
            {
                "schema_version": CLI_RESULT_SCHEMA,
                "command": None if request is None else request.command,
                "status": "error",
                "issue": DeploymentIssue(
                    code="DEPLOYMENT_FAILED",
                    message="Deployment operation failed.",
                ).to_dict(),
            },
        )
        return EXIT_OPERATION
    _write_json(output, envelope)
    return result.exit_code


__all__ = [
    "CLI_RESULT_SCHEMA",
    "COMMANDS",
    "DOCTOR",
    "DOCTOR_RESULT_SCHEMA",
    "EXIT_DATA",
    "EXIT_INCOMPATIBLE",
    "EXIT_OK",
    "EXIT_OPERATION",
    "EXIT_PERMISSION",
    "EXIT_UNAVAILABLE",
    "EXIT_USAGE",
    "INSTALL",
    "OPERATION_RESULT_SCHEMA",
    "ROLLBACK",
    "STATUS",
    "STATUS_RESULT_SCHEMA",
    "UPGRADE",
    "VERIFY",
    "VERIFY_RESULT_SCHEMA",
    "CommandBackend",
    "CommandRequest",
    "PublicCommandResult",
    "UnavailableCommandBackend",
    "execute_command",
    "main",
    "parse_command",
    "success_envelope",
]
