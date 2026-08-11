"""Supported deployment CLI composition and fixed sudo operator client."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
import grp
import hashlib
import json
import os
from pathlib import Path
import secrets
import stat
import subprocess

from encode_pipeline.deployment.bundle import BundleStore, MAX_BUNDLE_BYTES
from encode_pipeline.deployment.cli import (
    EXIT_DATA,
    EXIT_OK,
    EXIT_UNAVAILABLE,
    INSTALL,
    OPERATION_RESULT_SCHEMA,
    ROLLBACK,
    STATUS_RESULT_SCHEMA,
    UPGRADE,
    VERIFY_RESULT_SCHEMA,
    PublicCommandResult,
)
from encode_pipeline.deployment.doctor import (
    FAIL,
    PASS,
    WARNING,
    DeploymentDoctor,
    DeploymentSnapshot,
    DeploymentStateProbe,
    fixed_probe,
)
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.manager import (
    DeploymentManager,
    DeploymentOwnership,
    DeploymentStatus,
)
from encode_pipeline.deployment.models import (
    BULK_RNASEQ_RUNTIME,
    COMPONENTS,
    PLATFORM,
    BundleManifest,
)
from encode_pipeline.deployment.native_contracts import (
    PLATFORM_FRONTEND_CONTRACT,
)
from encode_pipeline.deployment.operator import (
    OPERATOR_RECEIPT_SCHEMA,
    SERVICE_UNITS,
    OperatorObservation,
    ServiceIdentity,
)
from encode_pipeline.deployment.operator_action import DeploymentActionReceipt


OPERATOR_EXECUTABLE = Path("/usr/libexec/helixweave-operator")
SUDO_EXECUTABLE = Path("/usr/bin/sudo")
OPERATOR_GROUP = "helixweave-operators"
SERVICE_GROUP = "helixweave"
_MAX_OPERATOR_RECEIPT_BYTES = 1024 * 1024
_OPERATOR_TIMEOUT_SECONDS = 2_700
_BULK_OPERATOR_TIMEOUT_SECONDS = 15_000
_TASK_PREFIX = "task-"

_SERVICE_NAMES = {
    "api": "helixweave-api.service",
    "worker": "helixweave-worker.service",
    "redis": "helixweave-redis.service",
    "bulk_rnaseq_docker": "helixweave-docker-rootless.service",
}


@dataclass(frozen=True)
class OperatorExecution:
    returncode: int
    stdout: bytes
    stderr: bytes


@dataclass(frozen=True)
class ServiceObservation:
    state: str
    identity: str | None = None

    def __post_init__(self) -> None:
        if self.state not in {"running", "stopped", "unavailable"} or (
            (self.state == "running") != (self.identity is not None)
        ):
            raise fail(
                "DEPLOYMENT_SERVICE_OBSERVATION_INVALID",
                "Deployment service observation is invalid.",
            )


class DeploymentOperatorClient:
    """Protocol-like base for the fixed operator calls used by the backend."""

    def stage(self, component: str, identity: str, task: str) -> None:
        raise NotImplementedError

    def activate(self, component: str, identity: str, task: str) -> None:
        raise NotImplementedError

    def rollback(self, component: str, identity: str, task: str) -> None:
        raise NotImplementedError

    def status_service(
        self,
        unit: str,
        deployment_identity: str,
        task: str,
    ) -> ServiceObservation:
        raise NotImplementedError

    def observe(self, state_identity: str, task: str) -> OperatorObservation:
        raise NotImplementedError

    def verify(self, state_identity: str, task: str) -> DeploymentActionReceipt:
        raise NotImplementedError


class SudoOperatorClient(DeploymentOperatorClient):
    """Invoke only the installed root-owned helper through ``sudo -n``."""

    def __init__(
        self,
        *,
        runner: Callable[[tuple[str, ...]], OperatorExecution] | None = None,
    ) -> None:
        self._runner = _run_operator if runner is None else runner

    def stage(self, component: str, identity: str, task: str) -> None:
        self._mutation("stage", component, identity, task, "staged")

    def activate(self, component: str, identity: str, task: str) -> None:
        self._mutation("activate", component, identity, task, "activated")

    def rollback(self, component: str, identity: str, task: str) -> None:
        self._mutation("rollback", component, identity, task, "rolled-back")

    def status_service(
        self,
        unit: str,
        deployment_identity: str,
        task: str,
    ) -> ServiceObservation:
        if unit not in SERVICE_UNITS:
            raise fail(
                "DEPLOYMENT_SERVICE_OBSERVATION_INVALID",
                "Deployment service observation is invalid.",
            )
        receipt = self._invoke(("status", unit, deployment_identity, task))
        common_keys = {
            "schema_version",
            "operation",
            "state",
            "task_identity",
            "deployment_identity",
            "unit",
        }
        state = receipt.get("state")
        expected_keys = common_keys | (
            {"service_identity"} if state == "running" else set()
        )
        if (
            set(receipt) != expected_keys
            or receipt.get("operation") != "status"
            or receipt.get("unit") != unit
            or receipt.get("deployment_identity") != deployment_identity
            or receipt.get("task_identity") != task
        ):
            raise _operator_failure()
        if state == "stopped":
            return ServiceObservation("stopped")
        raw_identity = receipt.get("service_identity")
        try:
            service = ServiceIdentity.from_dict(raw_identity)
        except DeploymentError:
            raise _operator_failure() from None
        if (
            state != "running"
            or service.unit != unit
            or service.deployment_identity != deployment_identity
        ):
            raise _operator_failure()
        return ServiceObservation("running", service.identity)

    def observe(self, state_identity: str, task: str) -> OperatorObservation:
        receipt = self._invoke(("observe", state_identity, task))
        if (
            set(receipt)
            != {
                "schema_version",
                "operation",
                "state",
                "task_identity",
                "deployment_identity",
                "observation",
            }
            or receipt.get("operation") != "observe"
            or receipt.get("state") != "observed"
            or receipt.get("task_identity") != task
            or receipt.get("deployment_identity") != state_identity
        ):
            raise _operator_failure()
        try:
            observation = OperatorObservation.from_dict(receipt.get("observation"))
        except DeploymentError:
            raise _operator_failure() from None
        if observation.state_identity != state_identity:
            raise _operator_failure()
        return observation

    def verify(self, state_identity: str, task: str) -> DeploymentActionReceipt:
        receipt = self._invoke(("verify", state_identity, task))
        if (
            set(receipt)
            != {
                "schema_version",
                "operation",
                "state",
                "task_identity",
                "deployment_identity",
                "verification",
            }
            or receipt.get("operation") != "verify"
            or receipt.get("state") != "verified"
            or receipt.get("task_identity") != task
            or receipt.get("deployment_identity") != state_identity
        ):
            raise _operator_failure()
        try:
            verification = DeploymentActionReceipt.from_dict(
                receipt.get("verification")
            )
        except DeploymentError:
            raise _operator_failure() from None
        if verification.status != "observed":
            raise _operator_failure()
        return verification

    def _mutation(
        self,
        operation: str,
        component: str,
        identity: str,
        task: str,
        state: str,
    ) -> None:
        receipt = self._invoke((operation, component, identity, task))
        if (
            set(receipt)
            != {
                "schema_version",
                "operation",
                "state",
                "task_identity",
                "deployment_identity",
                "component",
            }
            or receipt.get("operation") != operation
            or receipt.get("component") != component
            or receipt.get("deployment_identity") != identity
            or receipt.get("task_identity") != task
            or receipt.get("state") != state
            or "service_identity" in receipt
        ):
            raise _operator_failure()

    def _invoke(self, arguments: tuple[str, ...]) -> Mapping[str, object]:
        try:
            completed = self._runner(arguments)
        except Exception:
            raise _operator_failure(recoverable=True) from None
        if not isinstance(completed, OperatorExecution):
            raise _operator_failure()
        if (
            completed.returncode != 0
            or completed.stderr
            or not 0 < len(completed.stdout) <= _MAX_OPERATOR_RECEIPT_BYTES
        ):
            raise _operator_failure(recoverable=completed.returncode in {69, 77})
        try:
            document = json.loads(
                completed.stdout,
                object_pairs_hook=_unique_object,
                parse_constant=lambda _value: (_ for _ in ()).throw(ValueError()),
            )
        except (UnicodeError, ValueError, json.JSONDecodeError):
            raise _operator_failure() from None
        if (
            not isinstance(document, dict)
            or document.get("schema_version") != OPERATOR_RECEIPT_SCHEMA
            or completed.stdout
            != json.dumps(document, sort_keys=True, separators=(",", ":")).encode()
            + b"\n"
        ):
            raise _operator_failure()
        return document


class IngressPublisher:
    """Atomically publish a local bundle at the fixed operator ingress name."""

    def __init__(
        self,
        layout: DeploymentLayout,
        *,
        directory_uid: int,
        directory_gid: int,
        file_uid: int,
    ) -> None:
        self._layout = layout
        self._directory_uid = directory_uid
        self._directory_gid = directory_gid
        self._file_uid = file_uid

    @classmethod
    def supported(cls, layout: DeploymentLayout) -> "IngressPublisher":
        return cls(
            layout,
            directory_uid=0,
            directory_gid=_operator_group_gid(),
            file_uid=os.getuid(),
        )

    def publish(
        self,
        bundle_path: Path,
        manifest: BundleManifest,
        task: str,
    ) -> Path:
        directory = self._layout.ingress / manifest.component
        flags = (
            os.O_RDONLY
            | getattr(os, "O_CLOEXEC", 0)
            | getattr(os, "O_DIRECTORY", 0)
            | getattr(os, "O_NOFOLLOW", 0)
        )
        descriptor = -1
        temporary = f".{manifest.identity}.{task}.partial"
        target = f"{manifest.identity}.tar"
        try:
            before = directory.lstat()
            descriptor = os.open(directory, flags)
            opened = os.fstat(descriptor)
            if (
                not stat.S_ISDIR(opened.st_mode)
                or stat.S_ISLNK(before.st_mode)
                or stat.S_IMODE(opened.st_mode) != 0o2730
                or opened.st_uid != self._directory_uid
                or opened.st_gid != self._directory_gid
                or (opened.st_dev, opened.st_ino) != (before.st_dev, before.st_ino)
            ):
                raise OSError
            if self._existing_target(directory, target, manifest):
                return directory / target
            self._copy_to_temporary(
                bundle_path,
                descriptor,
                temporary,
            )
            copied = directory / temporary
            inspected = BundleStore(self._layout).inspect(
                copied,
                expected_owner_uid=self._file_uid,
                expected_owner_gid=self._directory_gid,
            )
            if inspected != manifest:
                raise OSError
            try:
                os.link(
                    temporary,
                    target,
                    src_dir_fd=descriptor,
                    dst_dir_fd=descriptor,
                    follow_symlinks=False,
                )
            except FileExistsError:
                if not self._existing_target(directory, target, manifest):
                    raise OSError
            os.unlink(temporary, dir_fd=descriptor)
            os.fsync(descriptor)
            return directory / target
        except DeploymentError:
            raise
        except OSError:
            raise fail(
                "DEPLOYMENT_INGRESS_UNAVAILABLE",
                "Deployment ingress is unavailable.",
                component=manifest.component,
                recoverable=True,
            ) from None
        finally:
            if descriptor >= 0:
                try:
                    os.unlink(temporary, dir_fd=descriptor)
                except FileNotFoundError:
                    pass
                os.close(descriptor)

    def _existing_target(
        self,
        directory: Path,
        target: str,
        manifest: BundleManifest,
    ) -> bool:
        path = directory / target
        try:
            observed = path.lstat()
        except FileNotFoundError:
            return False
        if (
            not stat.S_ISREG(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or observed.st_nlink != 1
            or stat.S_IMODE(observed.st_mode) != 0o440
            or observed.st_uid != self._file_uid
            or observed.st_gid != self._directory_gid
        ):
            raise OSError
        inspected = BundleStore(self._layout).inspect(
            path,
            expected_owner_uid=self._file_uid,
            expected_owner_gid=self._directory_gid,
        )
        return inspected == manifest

    def _copy_to_temporary(
        self,
        source: Path,
        directory_descriptor: int,
        temporary: str,
    ) -> None:
        source_flags = (
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
        )
        destination_flags = (
            os.O_WRONLY
            | os.O_CREAT
            | os.O_EXCL
            | getattr(os, "O_CLOEXEC", 0)
            | getattr(os, "O_NOFOLLOW", 0)
        )
        source_descriptor = os.open(source, source_flags)
        destination_descriptor = -1
        try:
            source_stat = os.fstat(source_descriptor)
            if (
                not stat.S_ISREG(source_stat.st_mode)
                or source_stat.st_nlink != 1
                or source_stat.st_mode & 0o022
                or not 0 < source_stat.st_size <= MAX_BUNDLE_BYTES
            ):
                raise OSError
            destination_descriptor = os.open(
                temporary,
                destination_flags,
                0o400,
                dir_fd=directory_descriptor,
            )
            total = 0
            while True:
                chunk = os.read(source_descriptor, 1024 * 1024)
                if not chunk:
                    break
                total += len(chunk)
                if total > MAX_BUNDLE_BYTES:
                    raise OSError
                written = 0
                while written < len(chunk):
                    written += os.write(destination_descriptor, chunk[written:])
            if total != source_stat.st_size:
                raise OSError
            os.fchown(destination_descriptor, self._file_uid, self._directory_gid)
            os.fchmod(destination_descriptor, 0o440)
            os.fsync(destination_descriptor)
            final = os.fstat(destination_descriptor)
            if (
                final.st_uid != self._file_uid
                or final.st_gid != self._directory_gid
                or final.st_nlink != 1
                or stat.S_IMODE(final.st_mode) != 0o440
            ):
                raise OSError
        finally:
            os.close(source_descriptor)
            if destination_descriptor >= 0:
                os.close(destination_descriptor)


class ProductionCommandBackend:
    """Compose read-only verification with the fixed privileged mutations."""

    def __init__(
        self,
        *,
        layout: DeploymentLayout,
        manager: DeploymentManager,
        operator: DeploymentOperatorClient,
        ingress: IngressPublisher,
        task_factory: Callable[[], str] | None = None,
    ) -> None:
        self.layout = layout
        self.manager = manager
        self.operator = operator
        self.ingress = ingress
        self._task_factory = _new_task if task_factory is None else task_factory

    @classmethod
    def supported(cls) -> "ProductionCommandBackend":
        layout = DeploymentLayout.supported()
        operator_group_gid = _operator_group_gid()
        manager = DeploymentManager(
            layout,
            ownership=DeploymentOwnership.root(),
            state_reader_gid=operator_group_gid,
            state_service_gid=_fixed_group_gid(SERVICE_GROUP),
            state_environment_content_verification=False,
        )
        return cls(
            layout=layout,
            manager=manager,
            operator=SudoOperatorClient(),
            ingress=IngressPublisher(
                layout,
                directory_uid=0,
                directory_gid=operator_group_gid,
                file_uid=os.getuid(),
            ),
        )

    def install(self, component: str, bundle_path: Path) -> PublicCommandResult:
        manifest = self._inspect_bundle(component, bundle_path)
        task = self._task()
        self.ingress.publish(bundle_path, manifest, task)
        self.operator.stage(component, manifest.identity, task)
        status = self.manager.status()
        if status.state.components[component].staged != manifest.identity:
            raise fail(
                "DEPLOYMENT_STATE_CHANGED",
                "Deployment state changed before the operation could complete.",
                component=component,
                recoverable=True,
            )
        return PublicCommandResult(
            _operation_result(INSTALL, component, manifest.identity, status)
        )

    def upgrade(self, component: str, bundle_path: Path) -> PublicCommandResult:
        manifest = self._inspect_bundle(component, bundle_path)
        task = self._task()
        stage_task = _operation_task(task, "stage")
        activate_task = _operation_task(task, "activate")
        self.ingress.publish(bundle_path, manifest, stage_task)
        self.operator.stage(component, manifest.identity, stage_task)
        self.operator.activate(component, manifest.identity, activate_task)
        status = self.manager.status()
        if status.state.components[component].active != manifest.identity:
            raise fail(
                "DEPLOYMENT_STATE_CHANGED",
                "Deployment state changed before the operation could complete.",
                component=component,
                recoverable=True,
            )
        return PublicCommandResult(
            _operation_result(UPGRADE, component, manifest.identity, status)
        )

    def rollback(
        self,
        component: str,
        deployment_identity: str,
    ) -> PublicCommandResult:
        before = self.manager.status()
        if before.state.components[component].previous != deployment_identity:
            raise fail(
                "DEPLOYMENT_PREVIOUS_IDENTITY_MISMATCH",
                "Previous deployment identity does not match the requested deployment.",
                component=component,
            )
        task = self._task()
        self.operator.rollback(component, deployment_identity, task)
        status = self.manager.status()
        if status.state.components[component].active != deployment_identity:
            raise fail(
                "DEPLOYMENT_STATE_CHANGED",
                "Deployment state changed before the operation could complete.",
                component=component,
                recoverable=True,
            )
        return PublicCommandResult(
            _operation_result(ROLLBACK, component, deployment_identity, status)
        )

    def status(self) -> PublicCommandResult:
        status = self.manager.status()
        observation = self._operator_observation(status)
        if observation is None:
            raise fail(
                "OPERATOR_OBSERVATION_UNAVAILABLE",
                "Operator observation is unavailable.",
                recoverable=True,
            )
        database_identity, database_reason = self._database_status(status, observation)
        services = self._service_statuses(status, observation)
        return PublicCommandResult(
            _status_result(
                status,
                database_identity,
                database_reason,
                services,
                observation,
            )
        )

    def verify(self) -> PublicCommandResult:
        status = self.manager.status()
        observation = self._operator_observation(status)
        if observation is None:
            raise fail(
                "OPERATOR_OBSERVATION_UNAVAILABLE",
                "Operator observation is unavailable.",
                recoverable=True,
            )
        services = self._service_statuses(status, observation)
        database_identity, database_reason = self._database_status(status, observation)
        frontend_identity = _frontend_identity(status)
        verification = self._operator_verification(
            status,
            observation,
            required=True,
        )
        try:
            verified = self.manager.verify_storage()
            database_identity, database_reason = self._database_status(
                verified, observation
            )
            if verification is not None:
                frontend_identity = verification.frontend_identity
                database_identity = verification.database_after_identity
                database_reason = (
                    "DATABASE_READY"
                    if database_identity is not None
                    else "DATABASE_UNAVAILABLE"
                )
            passed = verification is not None and _verification_ready(
                verification, verified
            )
        except DeploymentError:
            verified = status
            passed = False
        return PublicCommandResult(
            {
                "schema_version": VERIFY_RESULT_SCHEMA,
                "verified": passed,
                "deployment": _status_result(
                    verified,
                    database_identity,
                    database_reason,
                    services,
                    observation,
                ),
                "frontend_identity": frontend_identity,
                "database_schema_identity": database_identity,
            },
            EXIT_OK if passed else EXIT_DATA,
        )

    def doctor(self) -> PublicCommandResult:
        snapshot = DeploymentSnapshot(self.manager)
        try:
            status = snapshot.read()
        except DeploymentError:
            status = None
        observation = (
            self._operator_observation(status)
            if isinstance(status, DeploymentStatus)
            else None
        )
        services = (
            self._service_statuses(status, observation)
            if isinstance(status, DeploymentStatus)
            else _unavailable_services()
        )
        verification = (
            self._operator_verification(status, observation)
            if isinstance(status, DeploymentStatus) and observation is not None
            else None
        )
        probes = {
            "deployment-state": _deployment_state_probe(snapshot, observation),
            "configuration": _action_probe(
                verification,
                "configuration",
                ready_reason="CONFIGURATION_READY",
                failed_reason="CONFIGURATION_INVALID",
            ),
            "permissions": _action_probe(
                verification,
                "permissions",
                ready_reason="PERMISSIONS_READY",
                failed_reason="PERMISSIONS_INVALID",
            ),
            "database": _action_probe(
                verification,
                "database-schema",
                ready_reason="DATABASE_READY",
                failed_reason="DATABASE_SCHEMA_INCOMPATIBLE",
            ),
            "redis": _action_probe(
                verification,
                "redis",
                ready_reason="REDIS_READY",
                failed_reason="REDIS_UNAVAILABLE",
            ),
            "api": _service_probe(
                services["api"],
                ready_reason="API_READY",
                failed_reason="API_UNAVAILABLE",
            ),
            "worker": _service_probe(
                services["worker"],
                ready_reason="WORKER_READY",
                failed_reason="WORKER_UNAVAILABLE",
            ),
            "frontend": _combined_action_probe(
                verification,
                ("platform-native", "python-runtime", "frontend"),
                ready_reason="FRONTEND_READY",
                failed_reason="FRONTEND_ASSET_INTEGRITY_FAILED",
                evidence_identity=(
                    None if verification is None else verification.frontend_identity
                ),
            ),
            "encode-runtime": _action_probe(
                verification,
                "encode-runtime-native",
                ready_reason="RUNTIME_READY",
                failed_reason="RUNTIME_NOT_ACTIVE",
            ),
            "bulk-rnaseq-runtime": _combined_action_probe(
                verification,
                ("bulk-runtime-native", "docker"),
                ready_reason="RUNTIME_READY",
                failed_reason="RUNTIME_NOT_ACTIVE",
                evidence_identity=_native_identity(
                    verification,
                    BULK_RNASEQ_RUNTIME,
                ),
            ),
            "references": _reference_action_probe(verification),
        }
        report = DeploymentDoctor(probes).run()
        return PublicCommandResult(
            report.to_dict(),
            EXIT_OK if report.ready else EXIT_UNAVAILABLE,
        )

    def _inspect_bundle(self, component: str, path: Path) -> BundleManifest:
        if component not in COMPONENTS:
            raise fail(
                "DEPLOYMENT_COMPONENT_INVALID",
                "Deployment component is invalid.",
            )
        manifest = BundleStore(self.layout).inspect(path)
        if manifest.component != component:
            raise fail(
                "DEPLOYMENT_BUNDLE_COMPONENT_MISMATCH",
                "Deployment bundle component does not match the request.",
                component=component,
            )
        return manifest

    def _task(self) -> str:
        task = self._task_factory()
        if (
            not isinstance(task, str)
            or len(task) != 37
            or not task.startswith(_TASK_PREFIX)
            or any(character not in "0123456789abcdef" for character in task[5:])
        ):
            raise fail(
                "DEPLOYMENT_TASK_IDENTITY_INVALID",
                "Deployment task identity is invalid.",
            )
        return task

    def _database_status(
        self,
        status: DeploymentStatus,
        observation: OperatorObservation | None,
    ) -> tuple[str | None, str]:
        if observation is not None:
            return (
                observation.database_schema_identity,
                "DATABASE_READY"
                if observation.database_schema_identity is not None
                else "DATABASE_UNAVAILABLE",
            )
        return None, "DATABASE_UNAVAILABLE"

    def _operator_observation(
        self,
        status: DeploymentStatus,
    ) -> OperatorObservation | None:
        try:
            observation = self.operator.observe(status.state.identity, self._task())
        except Exception:
            return None
        active = {
            component: status.state.components[component].active
            for component in COMPONENTS
        }
        if (
            observation.state_identity != status.state.identity
            or observation.active != active
        ):
            return None
        return observation

    def _operator_verification(
        self,
        status: DeploymentStatus,
        observation: OperatorObservation | None,
        *,
        required: bool = False,
    ) -> DeploymentActionReceipt | None:
        try:
            verification = self.operator.verify(
                status.state.identity,
                self._task(),
            )
        except Exception:
            if required:
                raise fail(
                    "DEPLOYMENT_OPERATOR_UNAVAILABLE",
                    "Deployment operator is unavailable.",
                    recoverable=True,
                ) from None
            return None
        active = {
            component: status.state.components[component].active is not None
            for component in COMPONENTS
        }
        if any(
            (verification.native_identities[component] is not None) != active[component]
            for component in COMPONENTS
        ) or (
            observation is not None
            and verification.database_after_identity
            != observation.database_schema_identity
        ):
            return None
        return verification

    def _service_statuses(
        self,
        status: DeploymentStatus | None,
        observation: OperatorObservation | None = None,
    ) -> dict[str, ServiceObservation]:
        if status is None:
            return _unavailable_services()
        task = self._task()
        values: dict[str, ServiceObservation] = {}
        for name, unit in _SERVICE_NAMES.items():
            component = (
                BULK_RNASEQ_RUNTIME if name == "bulk_rnaseq_docker" else PLATFORM
            )
            deployment_identity = status.state.components[component].active
            if deployment_identity is None:
                values[name] = ServiceObservation("stopped")
                continue
            try:
                observed = self.operator.status_service(
                    unit,
                    deployment_identity,
                    task,
                )
                expected = None if observation is None else observation.services[unit]
                if observation is not None and (
                    (observed.state == "running" and observed.identity != expected)
                    or (observed.state == "stopped" and expected is not None)
                ):
                    values[name] = ServiceObservation("unavailable")
                else:
                    values[name] = observed
            except DeploymentError:
                values[name] = ServiceObservation("unavailable")
        return values


def _operation_result(
    operation: str,
    component: str,
    deployment_identity: str,
    status: DeploymentStatus,
) -> dict[str, object]:
    states = {
        INSTALL: "staged",
        UPGRADE: "activated",
        ROLLBACK: "rolled-back",
    }
    return {
        "schema_version": OPERATION_RESULT_SCHEMA,
        "operation": operation,
        "state": states[operation],
        "component": component,
        "deployment_identity": deployment_identity,
        "state_identity": status.state.identity,
    }


def _status_result(
    status: DeploymentStatus,
    database_identity: str | None,
    database_reason: str,
    services: Mapping[str, ServiceObservation],
    observation: OperatorObservation | None,
) -> dict[str, object]:
    pending_count = 0 if observation is None else observation.operator_pending_count
    recovery_count = (
        0 if observation is None else observation.operator_recovery_required_count
    )
    value = {
        "schema_version": STATUS_RESULT_SCHEMA,
        **status.to_dict(),
        "database_schema_identity": database_identity,
        "database_schema_reason_code": database_reason,
        "operator_pending_transaction_count": pending_count,
        "operator_recovery_required_count": recovery_count,
        "services": {
            name: {
                "state": services[name].state,
                "identity": services[name].identity,
            }
            for name in _SERVICE_NAMES
        },
    }
    value["interrupted"] = bool(value["interrupted"] or pending_count or recovery_count)
    return value


def _deployment_state_probe(
    snapshot: DeploymentSnapshot,
    observation: OperatorObservation | None,
):
    if observation is None:
        return fixed_probe(FAIL, "OPERATOR_UNAVAILABLE")
    if observation is not None and observation.operator_recovery_required_count:
        return fixed_probe(FAIL, "DEPLOYMENT_RECOVERY_REQUIRED")
    if observation is not None and observation.operator_pending_count:
        return fixed_probe(WARNING, "DEPLOYMENT_INTERRUPTED")
    return DeploymentStateProbe(snapshot)


def _frontend_identity(status: DeploymentStatus) -> str | None:
    manifest = status.manifests[PLATFORM]["active"]
    if manifest is None:
        return None
    identities = [
        item.identity
        for item in manifest.provides
        if item.contract == PLATFORM_FRONTEND_CONTRACT
    ]
    return identities[0] if len(identities) == 1 else None


def _service_probe(
    service: ServiceObservation,
    *,
    ready_reason: str,
    failed_reason: str,
):
    if service.state == "running":
        return fixed_probe(PASS, ready_reason, service.identity)
    return fixed_probe(FAIL, failed_reason)


def _action_probe(
    verification: DeploymentActionReceipt | None,
    check_name: str,
    *,
    ready_reason: str,
    failed_reason: str,
):
    if verification is None:
        return fixed_probe(FAIL, "OPERATOR_UNAVAILABLE")
    check = verification.readiness[check_name]
    if check.status == "ready":
        return fixed_probe(PASS, ready_reason, check.identity)
    return fixed_probe(FAIL, failed_reason)


def _combined_action_probe(
    verification: DeploymentActionReceipt | None,
    check_names: tuple[str, ...],
    *,
    ready_reason: str,
    failed_reason: str,
    evidence_identity: str | None,
):
    if verification is None:
        return fixed_probe(FAIL, "OPERATOR_UNAVAILABLE")
    checks = tuple(verification.readiness[name] for name in check_names)
    if all(check.status == "ready" for check in checks):
        evidence = evidence_identity or checks[0].identity
        return fixed_probe(PASS, ready_reason, evidence)
    return fixed_probe(FAIL, failed_reason)


def _native_identity(
    verification: DeploymentActionReceipt | None,
    component: str,
) -> str | None:
    if verification is None:
        return None
    return verification.native_identities[component]


def _reference_action_probe(
    verification: DeploymentActionReceipt | None,
):
    if verification is None:
        return fixed_probe(FAIL, "OPERATOR_UNAVAILABLE")
    check = verification.readiness["references"]
    if check.status == "ready":
        return fixed_probe(PASS, "REFERENCES_READY", check.identity)
    if check.status in {"not-ready", "not-applicable"}:
        return fixed_probe(WARNING, "REFERENCES_INCOMPLETE")
    return fixed_probe(FAIL, "REFERENCES_INVALID")


def _verification_ready(
    verification: DeploymentActionReceipt,
    status: DeploymentStatus,
) -> bool:
    if (
        verification.status != "observed"
        or verification.compatibility != "compatible"
        or verification.migration_required
        or status.state.components[PLATFORM].active is None
        or verification.database_after_identity is None
        or verification.accepted_schema_heads != verification.target_schema_heads
        or len(verification.target_schema_heads) != 1
        or verification.frontend_identity is None
        or verification.frontend_identity != _frontend_identity(status)
        or verification.reference_compatibility_identity is None
    ):
        return False
    active = {
        component: status.state.components[component].active is not None
        for component in COMPONENTS
    }
    required = {
        "platform-native": verification.native_identities[PLATFORM],
        "python-runtime": None,
        "frontend": verification.frontend_identity,
        "database-schema": verification.database_after_identity,
        "configuration": None,
        "permissions": None,
    }
    for component, check_name in (
        ("encode-runtime", "encode-runtime-native"),
        (BULK_RNASEQ_RUNTIME, "bulk-runtime-native"),
    ):
        check = verification.readiness[check_name]
        if active[component]:
            if (
                check.status != "ready"
                or check.identity != verification.native_identities[component]
            ):
                return False
        elif check.status != "not-applicable":
            return False
    for name, expected_identity in required.items():
        check = verification.readiness[name]
        if check.status != "ready" or (
            expected_identity is not None and check.identity != expected_identity
        ):
            return False
    return True


def _unavailable_services() -> dict[str, ServiceObservation]:
    return {name: ServiceObservation("unavailable") for name in _SERVICE_NAMES}


def _new_task() -> str:
    return f"{_TASK_PREFIX}{secrets.token_hex(16)}"


def _operation_task(parent: str, operation: str) -> str:
    if (
        not isinstance(parent, str)
        or len(parent) != 37
        or not parent.startswith(_TASK_PREFIX)
        or any(character not in "0123456789abcdef" for character in parent[5:])
        or operation not in {"stage", "activate"}
    ):
        raise fail(
            "DEPLOYMENT_TASK_IDENTITY_INVALID",
            "Deployment task identity is invalid.",
        )
    digest = hashlib.sha256(
        b"helixweave-deployment-operation-task-v1\0"
        + parent.encode("ascii")
        + b"\0"
        + operation.encode("ascii")
    ).hexdigest()
    return f"{_TASK_PREFIX}{digest[:32]}"


def _operator_group_gid() -> int:
    return _fixed_group_gid(OPERATOR_GROUP)


def _fixed_group_gid(name: str) -> int:
    try:
        group = grp.getgrnam(name)
    except KeyError:
        raise fail(
            "DEPLOYMENT_OPERATOR_UNAVAILABLE",
            "Deployment operator is unavailable.",
            recoverable=True,
        ) from None
    if not isinstance(group.gr_gid, int) or group.gr_gid <= 0:
        raise fail(
            "DEPLOYMENT_OPERATOR_UNAVAILABLE",
            "Deployment operator is unavailable.",
        )
    return group.gr_gid


def _run_operator(arguments: tuple[str, ...]) -> OperatorExecution:
    argv = (str(SUDO_EXECUTABLE), "-n", str(OPERATOR_EXECUTABLE), *arguments)
    timeout = (
        _BULK_OPERATOR_TIMEOUT_SECONDS
        if (
            arguments[:2]
            in {
                ("activate", BULK_RNASEQ_RUNTIME),
                ("rollback", BULK_RNASEQ_RUNTIME),
            }
            or arguments[:1] == ("verify",)
        )
        else _OPERATOR_TIMEOUT_SECONDS
    )
    try:
        completed = subprocess.run(
            argv,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=timeout,
            check=False,
            shell=False,
            cwd="/",
            env={
                "LANG": "C.UTF-8",
                "LC_ALL": "C.UTF-8",
                "PATH": "/usr/sbin:/usr/bin:/sbin:/bin",
            },
        )
    except (OSError, subprocess.TimeoutExpired):
        raise _operator_failure(recoverable=True) from None
    return OperatorExecution(
        completed.returncode,
        completed.stdout,
        completed.stderr,
    )


def _operator_failure(*, recoverable: bool = False) -> DeploymentError:
    return fail(
        "DEPLOYMENT_OPERATOR_UNAVAILABLE",
        "Deployment operator is unavailable.",
        recoverable=recoverable,
    )


def _unique_object(pairs: Sequence[tuple[str, object]]) -> dict[str, object]:
    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError
        value[key] = item
    return value


__all__ = [
    "DeploymentOperatorClient",
    "IngressPublisher",
    "OperatorExecution",
    "ProductionCommandBackend",
    "ServiceObservation",
    "SudoOperatorClient",
]
