"""Fail-closed contracts for the single supported deployment release Gate.

Requests are prepared only from a trusted, read-only observer.  The production
observer is unavailable during Phase A; consequently the public preparation
CLI cannot turn caller assertions into verified evidence or dispatch a Gate.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
import hashlib
import hmac
import json
import os
from pathlib import Path
import re
import stat
from typing import Any, Protocol

from encode_pipeline.deployment.canonical import (
    canonical_identity,
    canonical_json_bytes,
    without_key,
)
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.models import BULK_RNASEQ_RUNTIME, ENCODE_RUNTIME


CLEANUP_PLAN_SCHEMA = "helixweave-deployment-gate-cleanup-plan-v3"
GATE_OBSERVATION_SCHEMA = "helixweave-deployment-gate-observation-v2"
GATE_REQUEST_SCHEMA = "helixweave-deployment-gate-request-v3"
GATE_STAGE_VERIFICATION_SCHEMA = "helixweave-deployment-gate-stage-verification-v1"
GATE_STAGE_RECEIPT_SCHEMA = "helixweave-deployment-gate-stage-receipt-v2"
CLEANUP_PLAN_IDENTITY_SCHEME = "helixweave-deployment-gate-cleanup-plan-identity-v3"
GATE_OBSERVATION_IDENTITY_SCHEME = "helixweave-deployment-gate-observation-identity-v2"
GATE_REQUEST_IDENTITY_SCHEME = "helixweave-deployment-gate-request-identity-v3"
GATE_PROCESS_IDENTITY_SCHEME = "helixweave-deployment-gate-process-identity-v2"
GATE_STAGE_VERIFICATION_IDENTITY_SCHEME = (
    "helixweave-deployment-gate-stage-verification-identity-v1"
)
GATE_STAGE_RECEIPT_IDENTITY_SCHEME = (
    "helixweave-deployment-gate-stage-receipt-identity-v2"
)

SUPPORTED_GATE_ROOT = Path("/var/lib/helixweave/operator/deployment-gates")
CLEANUP_EXECUTOR = Path("/usr/libexec/helixweave-gate-cleanup")
OPERATOR_ROOT = Path("/opt/helixweave/operator")
FIXED_DISK_REQUIRED_BYTES = 20 * 1024**3
FIXED_LOAD_PER_CPU_MILLI = 2000

DELETE_RELATIVE_PATHS = {
    "runner-state": Path("runner"),
    "checkout": Path("checkout"),
    "venv": Path("venv"),
    "redis-state": Path("redis"),
    "dockerd-socket": Path("docker/docker.sock"),
    "process-state": Path("pids"),
}
DELETE_RESOURCE_KINDS = tuple(DELETE_RELATIVE_PATHS)
RETAIN_RESOURCE_KINDS = (
    "docker-data-root",
    "docker-images",
    "historical-evidence",
)
RESOURCE_KINDS = tuple(sorted((*DELETE_RESOURCE_KINDS, *RETAIN_RESOURCE_KINDS)))
REQUIRED_CACHE_KINDS = (
    "actions",
    "container-lock",
    "frontend-toolchain",
    "nextflow-toolchain",
    "python-toolchain",
)
RUNTIME_COMPONENTS = (BULK_RNASEQ_RUNTIME, ENCODE_RUNTIME)
PROCESS_NAMES = ("dockerd", "redis", "runner")
GATE_STAGES = (
    "identity",
    "checkout",
    "environment",
    "services",
    "runtime",
    "smoke",
    "evidence",
)

_IDENTITY = re.compile(r"^sha256-[0-9a-f]{64}$")
_TASK_IDENTITY = re.compile(r"^task-[0-9a-f]{32}$")
_HEAD_SHA = re.compile(r"^[0-9a-f]{40}$")
_SAFE_REASON = re.compile(r"^[A-Z][A-Z0-9_]{0,63}$")
_MAX_DISK_BYTES = 100 * 1024**4
_MAX_LOAD_MILLI = 1_000_000
_MAX_CONTROL_FILE_BYTES = 64 * 1024


def _document(raw: object, *, code: str) -> dict[str, Any]:
    if not isinstance(raw, Mapping) or any(not isinstance(key, str) for key in raw):
        raise fail(code, "Deployment Gate document is invalid.")
    return dict(raw)


def _identity(value: object, *, code: str) -> str:
    if not isinstance(value, str) or _IDENTITY.fullmatch(value) is None:
        raise fail(code, "Deployment Gate document is invalid.")
    return value


def _task_identity(value: object, *, code: str) -> str:
    if not isinstance(value, str) or _TASK_IDENTITY.fullmatch(value) is None:
        raise fail(code, "Deployment Gate document is invalid.")
    return value


def _positive_integer(value: object, *, code: str) -> int:
    if (
        not isinstance(value, int)
        or isinstance(value, bool)
        or not 0 < value <= 2**63 - 1
    ):
        raise fail(code, "Deployment Gate document is invalid.")
    return value


def _absolute_path(path: Path, *, code: str, existing: bool = False) -> Path:
    if (
        not isinstance(path, Path)
        or not path.is_absolute()
        or path == Path(path.anchor)
    ):
        raise fail(code, "Deployment Gate path is invalid.")
    rendered = str(path)
    if (
        "'" in rendered
        or len(rendered) > 4096
        or any(character in rendered for character in ("\x00", "\n", "\r"))
    ):
        raise fail(code, "Deployment Gate path is invalid.")
    try:
        if path.resolve(strict=existing) != path:
            raise fail(code, "Deployment Gate path boundary is unsafe.")
    except OSError:
        raise fail(code, "Deployment Gate path boundary is unsafe.") from None
    return path


def _strict_descendant(path: Path, parent: Path) -> bool:
    try:
        return path.relative_to(parent) != Path(".")
    except ValueError:
        return False


@dataclass(frozen=True)
class GatePolicy:
    """Fixed host policy; only ``supported`` is used by the public CLI."""

    gate_root: Path
    cleanup_executor: Path
    operator_root: Path
    operator_uid: int
    operator_gid: int
    root_uid: int
    root_gid: int
    disk_required_bytes: int = FIXED_DISK_REQUIRED_BYTES
    load_per_cpu_milli: int = FIXED_LOAD_PER_CPU_MILLI

    def __post_init__(self) -> None:
        _absolute_path(self.gate_root, code="GATE_POLICY_INVALID")
        _absolute_path(self.cleanup_executor, code="GATE_POLICY_INVALID")
        _absolute_path(self.operator_root, code="GATE_POLICY_INVALID")
        integers = (
            self.operator_uid,
            self.operator_gid,
            self.root_uid,
            self.root_gid,
            self.disk_required_bytes,
            self.load_per_cpu_milli,
        )
        if (
            any(
                not isinstance(value, int) or isinstance(value, bool) or value < 0
                for value in integers
            )
            or not 1024**3 <= self.disk_required_bytes <= _MAX_DISK_BYTES
            or not 1 <= self.load_per_cpu_milli <= 100_000
        ):
            raise fail("GATE_POLICY_INVALID", "Deployment Gate policy is invalid.")

    @classmethod
    def supported(cls) -> "GatePolicy":
        return cls(SUPPORTED_GATE_ROOT, CLEANUP_EXECUTOR, OPERATOR_ROOT, 0, 0, 0, 0)

    @classmethod
    def isolated(cls, root: Path, *, uid: int, gid: int) -> "GatePolicy":
        """Construct a test-only fixed layout; never selected by the CLI."""
        selected = _absolute_path(root, code="GATE_POLICY_INVALID")
        return cls(
            gate_root=selected / "deployment-gates",
            cleanup_executor=selected / "installed" / "helixweave-gate-cleanup",
            operator_root=selected / "installed" / "operator",
            operator_uid=uid,
            operator_gid=gid,
            root_uid=uid,
            root_gid=gid,
            disk_required_bytes=1024**3,
            load_per_cpu_milli=100_000,
        )

    def task_root(self, task_identity: str) -> Path:
        task = _task_identity(task_identity, code="GATE_REQUEST_INVALID")
        return self.gate_root / "tasks" / task

    def delete_paths(self, task_identity: str) -> dict[str, Path]:
        root = self.task_root(task_identity)
        return {
            kind: root / relative for kind, relative in DELETE_RELATIVE_PATHS.items()
        }

    def retain_paths(self) -> dict[str, Path]:
        return {
            "docker-data-root": self.gate_root / "docker-data",
            "docker-images": self.gate_root / "docker-data",
            "historical-evidence": self.gate_root / "evidence" / "history",
        }

    def output_directory(self, task_identity: str) -> Path:
        return self.task_root(task_identity) / "runner" / "prepared"


@dataclass(frozen=True, order=True)
class RuntimeIdentity:
    component: str
    identity: str

    @classmethod
    def from_dict(cls, raw: object, *, code: str) -> "RuntimeIdentity":
        value = _document(raw, code=code)
        if (
            set(value) != {"component", "identity"}
            or value["component"] not in RUNTIME_COMPONENTS
        ):
            raise fail(code, "Deployment Gate document is invalid.")
        return cls(value["component"], _identity(value["identity"], code=code))

    def to_dict(self) -> dict[str, str]:
        return {"component": self.component, "identity": self.identity}


@dataclass(frozen=True, order=True)
class ResourceEvidence:
    kind: str
    path_identity: str
    device: int
    inode: int

    @classmethod
    def create(
        cls, *, kind: str, path: Path, observed: os.stat_result
    ) -> "ResourceEvidence":
        if kind not in RESOURCE_KINDS:
            raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.")
        path_identity = canonical_identity(
            {"kind": kind, "absolute_path": str(path)},
            scheme="helixweave-deployment-gate-resource-path-v1",
        )
        return cls(
            kind,
            path_identity,
            _positive_integer(observed.st_dev, code="GATE_OBSERVATION_INVALID"),
            _positive_integer(observed.st_ino, code="GATE_OBSERVATION_INVALID"),
        )

    @classmethod
    def from_dict(cls, raw: object) -> "ResourceEvidence":
        code = "GATE_OBSERVATION_INVALID"
        value = _document(raw, code=code)
        if (
            set(value) != {"kind", "path_identity", "device", "inode"}
            or value["kind"] not in RESOURCE_KINDS
        ):
            raise fail(code, "Gate observation is invalid.")
        return cls(
            value["kind"],
            _identity(value["path_identity"], code=code),
            _positive_integer(value["device"], code=code),
            _positive_integer(value["inode"], code=code),
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "path_identity": self.path_identity,
            "device": self.device,
            "inode": self.inode,
        }


@dataclass(frozen=True, order=True)
class CacheEvidence:
    kind: str
    identity: str
    device: int
    inode: int

    @classmethod
    def from_dict(cls, raw: object) -> "CacheEvidence":
        code = "GATE_OBSERVATION_INVALID"
        value = _document(raw, code=code)
        if (
            set(value) != {"kind", "identity", "device", "inode"}
            or value["kind"] not in REQUIRED_CACHE_KINDS
        ):
            raise fail(code, "Gate observation is invalid.")
        return cls(
            value["kind"],
            _identity(value["identity"], code=code),
            _positive_integer(value["device"], code=code),
            _positive_integer(value["inode"], code=code),
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "identity": self.identity,
            "device": self.device,
            "inode": self.inode,
        }


@dataclass(frozen=True)
class GateProcessIdentity:
    identity: str
    name: str
    task_identity: str
    pid: int
    process_start_ticks: int
    executable_device: int
    executable_inode: int
    cmdline_identity: str
    socket_device: int | None
    socket_inode: int | None

    @classmethod
    def create(
        cls,
        *,
        name: str,
        task_identity: str,
        pid: int,
        process_start_ticks: int,
        executable_device: int,
        executable_inode: int,
        cmdline_identity: str,
        socket_device: int | None = None,
        socket_inode: int | None = None,
    ) -> "GateProcessIdentity":
        value: dict[str, object] = {
            "name": name,
            "task_identity": task_identity,
            "pid": pid,
            "process_start_ticks": process_start_ticks,
            "executable_device": executable_device,
            "executable_inode": executable_inode,
            "cmdline_identity": cmdline_identity,
            "socket_device": socket_device,
            "socket_inode": socket_inode,
        }
        value["identity"] = canonical_identity(
            value, scheme=GATE_PROCESS_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "GateProcessIdentity":
        code = "GATE_PROCESS_IDENTITY_INVALID"
        value = _document(raw, code=code)
        if (
            set(value)
            != {
                "identity",
                "name",
                "task_identity",
                "pid",
                "process_start_ticks",
                "executable_device",
                "executable_inode",
                "cmdline_identity",
                "socket_device",
                "socket_inode",
            }
            or value["name"] not in PROCESS_NAMES
        ):
            raise fail(code, "Deployment Gate process identity is invalid.")
        observed_identity = _identity(value["identity"], code=code)
        task = _task_identity(value["task_identity"], code=code)
        pid = _positive_integer(value["pid"], code=code)
        start = _positive_integer(value["process_start_ticks"], code=code)
        executable_device = _positive_integer(value["executable_device"], code=code)
        executable_inode = _positive_integer(value["executable_inode"], code=code)
        cmdline = _identity(value["cmdline_identity"], code=code)
        sockets = value["socket_device"], value["socket_inode"]
        if pid > 2**31 - 1 or (
            (value["name"] == "runner" and sockets != (None, None))
            or (
                value["name"] != "runner"
                and any(
                    not isinstance(item, int)
                    or isinstance(item, bool)
                    or not 0 < item <= 2**63 - 1
                    for item in sockets
                )
            )
        ):
            raise fail(code, "Deployment Gate process identity is invalid.")
        expected = canonical_identity(
            without_key(value, "identity"), scheme=GATE_PROCESS_IDENTITY_SCHEME
        )
        if observed_identity != expected:
            raise fail(code, "Deployment Gate process identity is invalid.")
        return cls(
            observed_identity,
            value["name"],
            task,
            pid,
            start,
            executable_device,
            executable_inode,
            cmdline,
            value["socket_device"],
            value["socket_inode"],
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "identity": self.identity,
            "name": self.name,
            "task_identity": self.task_identity,
            "pid": self.pid,
            "process_start_ticks": self.process_start_ticks,
            "executable_device": self.executable_device,
            "executable_inode": self.executable_inode,
            "cmdline_identity": self.cmdline_identity,
            "socket_device": self.socket_device,
            "socket_inode": self.socket_inode,
        }


@dataclass(frozen=True)
class CleanupExecutorEvidence:
    path_identity: str
    file_identity: str
    device: int
    inode: int
    closure_identity: str
    backend_identity: str
    backend_device: int
    backend_inode: int

    @classmethod
    def from_dict(cls, raw: object) -> "CleanupExecutorEvidence":
        code = "GATE_OBSERVATION_INVALID"
        value = _document(raw, code=code)
        if set(value) != {
            "path_identity",
            "file_identity",
            "device",
            "inode",
            "closure_identity",
            "backend_identity",
            "backend_device",
            "backend_inode",
        }:
            raise fail(code, "Gate observation is invalid.")
        return cls(
            _identity(value["path_identity"], code=code),
            _identity(value["file_identity"], code=code),
            _positive_integer(value["device"], code=code),
            _positive_integer(value["inode"], code=code),
            _identity(value["closure_identity"], code=code),
            _identity(value["backend_identity"], code=code),
            _positive_integer(value["backend_device"], code=code),
            _positive_integer(value["backend_inode"], code=code),
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "path_identity": self.path_identity,
            "file_identity": self.file_identity,
            "device": self.device,
            "inode": self.inode,
            "closure_identity": self.closure_identity,
            "backend_identity": self.backend_identity,
            "backend_device": self.backend_device,
            "backend_inode": self.backend_inode,
        }


@dataclass(frozen=True)
class GateObservation:
    identity: str
    task_identity: str
    head_sha: str
    release_identity: str
    runtime_identities: tuple[RuntimeIdentity, ...]
    environment_identity: str
    resources: tuple[ResourceEvidence, ...]
    processes: tuple[GateProcessIdentity, ...]
    caches: tuple[CacheEvidence, ...]
    cleanup_executor: CleanupExecutorEvidence
    disk_free_bytes: int
    disk_required_bytes: int
    load_milli: int
    load_limit_milli: int

    @classmethod
    def create(
        cls,
        *,
        task_identity: str,
        head_sha: str,
        release_identity: str,
        runtime_identities: Mapping[str, str],
        environment_identity: str,
        resources: Sequence[ResourceEvidence],
        processes: Sequence[GateProcessIdentity],
        caches: Sequence[CacheEvidence],
        cleanup_executor: CleanupExecutorEvidence,
        disk_free_bytes: int,
        disk_required_bytes: int,
        load_milli: int,
        load_limit_milli: int,
    ) -> "GateObservation":
        value: dict[str, object] = {
            "schema_version": GATE_OBSERVATION_SCHEMA,
            "task_identity": task_identity,
            "head_sha": head_sha,
            "release_identity": release_identity,
            "runtime_identities": [
                {"component": component, "identity": runtime_identities[component]}
                for component in RUNTIME_COMPONENTS
                if component in runtime_identities
            ],
            "environment_identity": environment_identity,
            "resources": [item.to_dict() for item in sorted(resources)],
            "processes": [
                item.to_dict() for item in sorted(processes, key=lambda item: item.name)
            ],
            "caches": [item.to_dict() for item in sorted(caches)],
            "cleanup_executor": cleanup_executor.to_dict(),
            "disk_free_bytes": disk_free_bytes,
            "disk_required_bytes": disk_required_bytes,
            "load_milli": load_milli,
            "load_limit_milli": load_limit_milli,
        }
        value["identity"] = canonical_identity(
            value, scheme=GATE_OBSERVATION_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "GateObservation":
        code = "GATE_OBSERVATION_INVALID"
        value = _document(raw, code=code)
        if (
            set(value)
            != {
                "schema_version",
                "identity",
                "task_identity",
                "head_sha",
                "release_identity",
                "runtime_identities",
                "environment_identity",
                "resources",
                "processes",
                "caches",
                "cleanup_executor",
                "disk_free_bytes",
                "disk_required_bytes",
                "load_milli",
                "load_limit_milli",
            }
            or value["schema_version"] != GATE_OBSERVATION_SCHEMA
        ):
            raise fail(code, "Gate observation is invalid.")
        sequence_keys = ("runtime_identities", "resources", "processes", "caches")
        if any(not isinstance(value[key], list) for key in sequence_keys):
            raise fail(code, "Gate observation is invalid.")
        task = _task_identity(value["task_identity"], code=code)
        if (
            not isinstance(value["head_sha"], str)
            or _HEAD_SHA.fullmatch(value["head_sha"]) is None
        ):
            raise fail(code, "Gate observation is invalid.")
        release = _identity(value["release_identity"], code=code)
        environment = _identity(value["environment_identity"], code=code)
        runtimes = tuple(
            RuntimeIdentity.from_dict(item, code=code)
            for item in value["runtime_identities"]
        )
        resources = tuple(
            ResourceEvidence.from_dict(item) for item in value["resources"]
        )
        processes = tuple(
            GateProcessIdentity.from_dict(item) for item in value["processes"]
        )
        caches = tuple(CacheEvidence.from_dict(item) for item in value["caches"])
        executor = CleanupExecutorEvidence.from_dict(value["cleanup_executor"])
        integers = tuple(
            value[key]
            for key in (
                "disk_free_bytes",
                "disk_required_bytes",
                "load_milli",
                "load_limit_milli",
            )
        )
        if (
            tuple(item.component for item in runtimes) != RUNTIME_COMPONENTS
            or tuple(item.kind for item in resources) != RESOURCE_KINDS
            or tuple(item.name for item in processes) != PROCESS_NAMES
            or any(item.task_identity != task for item in processes)
            or len({item.pid for item in processes}) != len(processes)
            or tuple(item.kind for item in caches) != REQUIRED_CACHE_KINDS
            or any(
                not isinstance(item, int) or isinstance(item, bool) for item in integers
            )
            or not 1024**3 <= value["disk_required_bytes"] <= _MAX_DISK_BYTES
            or not value["disk_required_bytes"]
            <= value["disk_free_bytes"]
            <= _MAX_DISK_BYTES
            or not 0 <= value["load_milli"] <= value["load_limit_milli"]
            or not 1 <= value["load_limit_milli"] <= _MAX_LOAD_MILLI
        ):
            raise fail(code, "Gate observation is invalid.")
        observed = _identity(value["identity"], code=code)
        if observed != canonical_identity(
            without_key(value, "identity"), scheme=GATE_OBSERVATION_IDENTITY_SCHEME
        ):
            raise fail(code, "Gate observation is invalid.")
        return cls(
            observed,
            task,
            value["head_sha"],
            release,
            runtimes,
            environment,
            resources,
            processes,
            caches,
            executor,
            *integers,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": GATE_OBSERVATION_SCHEMA,
            "identity": self.identity,
            "task_identity": self.task_identity,
            "head_sha": self.head_sha,
            "release_identity": self.release_identity,
            "runtime_identities": [item.to_dict() for item in self.runtime_identities],
            "environment_identity": self.environment_identity,
            "resources": [item.to_dict() for item in self.resources],
            "processes": [item.to_dict() for item in self.processes],
            "caches": [item.to_dict() for item in self.caches],
            "cleanup_executor": self.cleanup_executor.to_dict(),
            "disk_free_bytes": self.disk_free_bytes,
            "disk_required_bytes": self.disk_required_bytes,
            "load_milli": self.load_milli,
            "load_limit_milli": self.load_limit_milli,
        }


class GateObserver(Protocol):
    def observe(self, policy: GatePolicy, task_identity: str) -> GateObservation: ...


class UnavailableGateObserver:
    def observe(self, policy: GatePolicy, task_identity: str) -> GateObservation:
        del policy, task_identity
        raise fail(
            "GATE_OBSERVER_UNAVAILABLE",
            "Trusted deployment Gate observation is not installed.",
            recoverable=True,
        )


class FilesystemGateObserver:
    """Read fixed owned paths, proc identities, sockets, and cache objects."""

    def observe(self, policy: GatePolicy, task_identity: str) -> GateObservation:
        task = _task_identity(task_identity, code="GATE_OBSERVATION_INVALID")
        root = policy.task_root(task)
        _owned_directory(policy.gate_root, policy.operator_uid, policy.operator_gid)
        _owned_directory(root, policy.operator_uid, policy.operator_gid)
        identities = _read_canonical_json(
            root / "identities.json", policy.operator_uid, policy.operator_gid
        )
        if (
            set(identities)
            != {
                "schema_version",
                "head_sha",
                "release_identity",
                "runtime_identities",
            }
            or identities["schema_version"] != "helixweave-gate-identities-v1"
        ):
            raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.")
        environment = _read_canonical_json(
            root / "venv" / "install.json",
            policy.operator_uid,
            policy.operator_gid,
        )
        if (
            set(environment)
            != {"schema_version", "install_mode", "editable", "identity"}
            or environment["schema_version"] != "helixweave-gate-environment-v1"
            or environment["install_mode"] != "wheel"
            or environment["editable"] is not False
        ):
            raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.")

        resources = tuple(
            sorted(
                _observe_resources(policy, task),
                key=lambda item: item.kind,
            )
        )
        processes = tuple(
            sorted(
                (_observe_process(policy, task, name) for name in PROCESS_NAMES),
                key=lambda item: item.name,
            )
        )
        caches = tuple(_observe_cache(policy, kind) for kind in REQUIRED_CACHE_KINDS)
        executor = _observe_cleanup_executor(policy)
        try:
            statvfs = os.statvfs(policy.gate_root)
            disk_free = statvfs.f_bavail * statvfs.f_frsize
            load_milli = max(0, round(os.getloadavg()[0] * 1000))
        except OSError:
            raise fail(
                "GATE_OBSERVATION_INVALID", "Gate observation is invalid."
            ) from None
        cpu_count = os.cpu_count() or 1
        load_limit = min(_MAX_LOAD_MILLI, cpu_count * policy.load_per_cpu_milli)
        runtime_map = _runtime_map(identities["runtime_identities"])
        return GateObservation.create(
            task_identity=task,
            head_sha=identities["head_sha"],
            release_identity=identities["release_identity"],
            runtime_identities=runtime_map,
            environment_identity=environment["identity"],
            resources=resources,
            processes=processes,
            caches=caches,
            cleanup_executor=executor,
            disk_free_bytes=disk_free,
            disk_required_bytes=policy.disk_required_bytes,
            load_milli=load_milli,
            load_limit_milli=load_limit,
        )


@dataclass(frozen=True, order=True)
class CleanupTarget:
    kind: str
    disposition: str
    path: Path
    path_identity: str
    device: int
    inode: int

    @classmethod
    def create(
        cls,
        *,
        kind: str,
        disposition: str,
        path: Path,
        device: int,
        inode: int,
    ) -> "CleanupTarget":
        code = "GATE_CLEANUP_PLAN_INVALID"
        selected = _absolute_path(path, code=code)
        if (
            kind not in RESOURCE_KINDS
            or disposition not in ("delete", "retain")
            or (kind in DELETE_RESOURCE_KINDS) != (disposition == "delete")
        ):
            raise fail(code, "Deployment Gate cleanup plan is invalid.")
        return cls(
            kind,
            disposition,
            selected,
            canonical_identity(
                {"kind": kind, "absolute_path": str(selected)},
                scheme="helixweave-deployment-gate-resource-path-v1",
            ),
            _positive_integer(device, code=code),
            _positive_integer(inode, code=code),
        )

    @classmethod
    def from_dict(cls, raw: object) -> "CleanupTarget":
        code = "GATE_CLEANUP_PLAN_INVALID"
        value = _document(raw, code=code)
        if set(value) != {
            "kind",
            "disposition",
            "absolute_path",
            "path_identity",
            "device",
            "inode",
        } or not isinstance(value["absolute_path"], str):
            raise fail(code, "Deployment Gate cleanup plan is invalid.")
        created = cls.create(
            kind=value["kind"],
            disposition=value["disposition"],
            path=Path(value["absolute_path"]),
            device=value["device"],
            inode=value["inode"],
        )
        if value["path_identity"] != created.path_identity:
            raise fail(code, "Deployment Gate cleanup plan is invalid.")
        return created

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "disposition": self.disposition,
            "absolute_path": str(self.path),
            "path_identity": self.path_identity,
            "device": self.device,
            "inode": self.inode,
        }


@dataclass(frozen=True)
class CleanupPlan:
    identity: str
    task_identity: str
    observation_identity: str
    targets: tuple[CleanupTarget, ...]
    processes: tuple[GateProcessIdentity, ...]
    executor_path: Path
    executor: CleanupExecutorEvidence

    @classmethod
    def create(cls, policy: GatePolicy, observation: GateObservation) -> "CleanupPlan":
        task = observation.task_identity
        delete_paths = policy.delete_paths(task)
        retain_paths = policy.retain_paths()
        resource_evidence = {item.kind: item for item in observation.resources}
        targets = tuple(
            sorted(
                (
                    *(
                        CleanupTarget.create(
                            kind=kind,
                            disposition="delete",
                            path=delete_paths[kind],
                            device=resource_evidence[kind].device,
                            inode=resource_evidence[kind].inode,
                        )
                        for kind in DELETE_RESOURCE_KINDS
                    ),
                    *(
                        CleanupTarget.create(
                            kind=kind,
                            disposition="retain",
                            path=retain_paths[kind],
                            device=resource_evidence[kind].device,
                            inode=resource_evidence[kind].inode,
                        )
                        for kind in RETAIN_RESOURCE_KINDS
                    ),
                )
            )
        )
        cls._validate_targets(policy, task, targets)
        observed_paths = {
            item.kind: item.path_identity for item in observation.resources
        }
        if any(
            observed_paths.get(target.kind) != target.path_identity
            for target in targets
        ):
            raise fail(
                "GATE_CLEANUP_PLAN_INVALID",
                "Deployment Gate cleanup plan is invalid.",
            )
        expected_executor_path_identity = canonical_identity(
            {"absolute_path": str(policy.cleanup_executor)},
            scheme="helixweave-deployment-gate-cleanup-executor-path-v1",
        )
        if (
            observation.cleanup_executor.path_identity
            != expected_executor_path_identity
        ):
            raise fail(
                "GATE_CLEANUP_PLAN_INVALID",
                "Deployment Gate cleanup plan is invalid.",
            )
        value: dict[str, object] = {
            "schema_version": CLEANUP_PLAN_SCHEMA,
            "task_identity": task,
            "observation_identity": observation.identity,
            "targets": [item.to_dict() for item in targets],
            "processes": [item.to_dict() for item in observation.processes],
            "executor_path": str(policy.cleanup_executor),
            "executor": observation.cleanup_executor.to_dict(),
        }
        value["identity"] = canonical_identity(
            value, scheme=CLEANUP_PLAN_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "CleanupPlan":
        code = "GATE_CLEANUP_PLAN_INVALID"
        value = _document(raw, code=code)
        if (
            set(value)
            != {
                "schema_version",
                "identity",
                "task_identity",
                "observation_identity",
                "targets",
                "processes",
                "executor_path",
                "executor",
            }
            or value["schema_version"] != CLEANUP_PLAN_SCHEMA
        ):
            raise fail(code, "Deployment Gate cleanup plan is invalid.")
        if (
            not isinstance(value["targets"], list)
            or not isinstance(value["processes"], list)
            or not isinstance(value["executor_path"], str)
        ):
            raise fail(code, "Deployment Gate cleanup plan is invalid.")
        task = _task_identity(value["task_identity"], code=code)
        observation_identity = _identity(value["observation_identity"], code=code)
        targets = tuple(CleanupTarget.from_dict(item) for item in value["targets"])
        processes = tuple(
            GateProcessIdentity.from_dict(item) for item in value["processes"]
        )
        executor_path = _absolute_path(Path(value["executor_path"]), code=code)
        executor = CleanupExecutorEvidence.from_dict(value["executor"])
        expected_executor_path_identity = canonical_identity(
            {"absolute_path": str(executor_path)},
            scheme="helixweave-deployment-gate-cleanup-executor-path-v1",
        )
        if (
            tuple(item.kind for item in targets) != RESOURCE_KINDS
            or tuple(item.name for item in processes) != PROCESS_NAMES
            or any(item.task_identity != task for item in processes)
            or len({item.pid for item in processes}) != len(processes)
            or executor.path_identity != expected_executor_path_identity
        ):
            raise fail(code, "Deployment Gate cleanup plan is invalid.")
        cls._validate_overlap(targets)
        observed = _identity(value["identity"], code=code)
        if observed != canonical_identity(
            without_key(value, "identity"), scheme=CLEANUP_PLAN_IDENTITY_SCHEME
        ):
            raise fail(code, "Deployment Gate cleanup plan is invalid.")
        return cls(
            observed,
            task,
            observation_identity,
            targets,
            processes,
            executor_path,
            executor,
        )

    @staticmethod
    def _validate_targets(
        policy: GatePolicy,
        task_identity: str,
        targets: tuple[CleanupTarget, ...],
    ) -> None:
        expected = {
            **policy.delete_paths(task_identity),
            **policy.retain_paths(),
        }
        if any(target.path != expected[target.kind] for target in targets):
            raise fail(
                "GATE_CLEANUP_PLAN_INVALID",
                "Deployment Gate cleanup plan is invalid.",
            )
        CleanupPlan._validate_overlap(targets)

    @staticmethod
    def _validate_overlap(targets: tuple[CleanupTarget, ...]) -> None:
        deleted = [item.path for item in targets if item.disposition == "delete"]
        retained = [item.path for item in targets if item.disposition == "retain"]
        if any(
            deleted_path == retained_path
            or _strict_descendant(deleted_path, retained_path)
            or _strict_descendant(retained_path, deleted_path)
            for deleted_path in deleted
            for retained_path in retained
        ):
            raise fail(
                "GATE_CLEANUP_PLAN_INVALID",
                "Deployment Gate cleanup plan is invalid.",
            )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": CLEANUP_PLAN_SCHEMA,
            "identity": self.identity,
            "task_identity": self.task_identity,
            "observation_identity": self.observation_identity,
            "targets": [item.to_dict() for item in self.targets],
            "processes": [item.to_dict() for item in self.processes],
            "executor_path": str(self.executor_path),
            "executor": self.executor.to_dict(),
        }


def render_cleanup_script(plan: CleanupPlan, *, policy: GatePolicy) -> bytes:
    checked = CleanupPlan.from_dict(plan.to_dict())
    if checked.executor_path != policy.cleanup_executor:
        raise fail("GATE_CLEANUP_SCRIPT_INVALID", "Gate cleanup script is invalid.")
    content = (
        "#!/bin/sh\n"
        "set -eu\n"
        f"exec /usr/bin/sudo -n '{checked.executor_path}' cleanup "
        f"--task-id '{checked.task_identity}' "
        f"--plan-id '{checked.identity}' "
        f"--executor-id '{checked.executor.file_identity}' "
        f"--closure-id '{checked.executor.closure_identity}' "
        f"--backend-id '{checked.executor.backend_identity}'\n"
    )
    return content.encode("ascii")


def cleanup_script_identity(content: bytes) -> str:
    if not isinstance(content, bytes) or not 0 < len(content) <= 16 * 1024:
        raise fail("GATE_CLEANUP_SCRIPT_INVALID", "Gate cleanup script is invalid.")
    return f"sha256-{hashlib.sha256(content).hexdigest()}"


def verify_cleanup_script(
    plan: CleanupPlan, content: bytes, *, policy: GatePolicy
) -> str:
    if not hmac.compare_digest(render_cleanup_script(plan, policy=policy), content):
        raise fail("GATE_CLEANUP_SCRIPT_INVALID", "Gate cleanup script is invalid.")
    return cleanup_script_identity(content)


@dataclass(frozen=True)
class GateRequest:
    identity: str
    head_sha: str
    release_identity: str
    runtime_identities: tuple[RuntimeIdentity, ...]
    task_identity: str
    observation_identity: str
    environment_identity: str
    cleanup_plan_identity: str
    cleanup_script_identity: str

    @classmethod
    def create(
        cls,
        *,
        policy: GatePolicy,
        observation: GateObservation,
        cleanup_plan: CleanupPlan,
        cleanup_script: bytes,
    ) -> "GateRequest":
        observation = GateObservation.from_dict(observation.to_dict())
        plan = CleanupPlan.from_dict(cleanup_plan.to_dict())
        if (
            plan.task_identity != observation.task_identity
            or plan.observation_identity != observation.identity
        ):
            raise fail("GATE_REQUEST_INVALID", "Deployment Gate request is invalid.")
        script_identity = verify_cleanup_script(plan, cleanup_script, policy=policy)
        value: dict[str, object] = {
            "schema_version": GATE_REQUEST_SCHEMA,
            "head_sha": observation.head_sha,
            "release_identity": observation.release_identity,
            "runtime_identities": [
                item.to_dict() for item in observation.runtime_identities
            ],
            "task_identity": observation.task_identity,
            "observation_identity": observation.identity,
            "environment_identity": observation.environment_identity,
            "cleanup_plan_identity": plan.identity,
            "cleanup_script_identity": script_identity,
        }
        value["identity"] = canonical_identity(
            value, scheme=GATE_REQUEST_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "GateRequest":
        code = "GATE_REQUEST_INVALID"
        value = _document(raw, code=code)
        if (
            set(value)
            != {
                "schema_version",
                "identity",
                "head_sha",
                "release_identity",
                "runtime_identities",
                "task_identity",
                "observation_identity",
                "environment_identity",
                "cleanup_plan_identity",
                "cleanup_script_identity",
            }
            or value["schema_version"] != GATE_REQUEST_SCHEMA
        ):
            raise fail(code, "Deployment Gate request is invalid.")
        if (
            not isinstance(value["head_sha"], str)
            or _HEAD_SHA.fullmatch(value["head_sha"]) is None
            or not isinstance(value["runtime_identities"], list)
        ):
            raise fail(code, "Deployment Gate request is invalid.")
        runtimes = tuple(
            RuntimeIdentity.from_dict(item, code=code)
            for item in value["runtime_identities"]
        )
        if tuple(item.component for item in runtimes) != RUNTIME_COMPONENTS:
            raise fail(code, "Deployment Gate request is invalid.")
        observed = _identity(value["identity"], code=code)
        identities = tuple(
            _identity(value[key], code=code)
            for key in (
                "release_identity",
                "observation_identity",
                "environment_identity",
                "cleanup_plan_identity",
                "cleanup_script_identity",
            )
        )
        task = _task_identity(value["task_identity"], code=code)
        if observed != canonical_identity(
            without_key(value, "identity"), scheme=GATE_REQUEST_IDENTITY_SCHEME
        ):
            raise fail(code, "Deployment Gate request is invalid.")
        return cls(
            observed, value["head_sha"], identities[0], runtimes, task, *identities[1:]
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": GATE_REQUEST_SCHEMA,
            "identity": self.identity,
            "head_sha": self.head_sha,
            "release_identity": self.release_identity,
            "runtime_identities": [item.to_dict() for item in self.runtime_identities],
            "task_identity": self.task_identity,
            "observation_identity": self.observation_identity,
            "environment_identity": self.environment_identity,
            "cleanup_plan_identity": self.cleanup_plan_identity,
            "cleanup_script_identity": self.cleanup_script_identity,
        }


def prepare_gate_request(
    *,
    policy: GatePolicy,
    observer: GateObserver,
    task_identity: str,
    head_sha: str,
    release_identity: str,
    runtime_identities: Mapping[str, str],
) -> tuple[CleanupPlan, bytes, GateRequest]:
    """Observe fixed resources and prepare, but never dispatch, one request."""
    task = _task_identity(task_identity, code="GATE_REQUEST_INVALID")
    if not isinstance(head_sha, str) or _HEAD_SHA.fullmatch(head_sha) is None:
        raise fail("GATE_REQUEST_INVALID", "Deployment Gate request is invalid.")
    release = _identity(release_identity, code="GATE_REQUEST_INVALID")
    if set(runtime_identities) != set(RUNTIME_COMPONENTS):
        raise fail("GATE_REQUEST_INVALID", "Deployment Gate request is invalid.")
    expected_runtimes = {
        component: _identity(runtime_identities[component], code="GATE_REQUEST_INVALID")
        for component in RUNTIME_COMPONENTS
    }
    try:
        observation = GateObservation.from_dict(
            observer.observe(policy, task).to_dict()
        )
    except DeploymentError as error:
        if error.issue.code == "GATE_OBSERVER_UNAVAILABLE":
            raise
        raise fail(
            "GATE_OBSERVATION_FAILED", "Trusted Gate observation failed."
        ) from None
    except Exception:
        raise fail(
            "GATE_OBSERVATION_FAILED", "Trusted Gate observation failed."
        ) from None
    if (
        observation.task_identity != task
        or observation.head_sha != head_sha
        or observation.release_identity != release
        or {item.component: item.identity for item in observation.runtime_identities}
        != expected_runtimes
        or observation.disk_required_bytes != policy.disk_required_bytes
        or observation.load_limit_milli
        != min(
            _MAX_LOAD_MILLI,
            (os.cpu_count() or 1) * policy.load_per_cpu_milli,
        )
    ):
        raise fail(
            "GATE_OBSERVATION_MISMATCH", "Gate observation does not match request."
        )
    plan = CleanupPlan.create(policy, observation)
    script = render_cleanup_script(plan, policy=policy)
    request = GateRequest.create(
        policy=policy,
        observation=observation,
        cleanup_plan=plan,
        cleanup_script=script,
    )
    return plan, script, request


@dataclass(frozen=True)
class GateStageVerification:
    identity: str
    request_identity: str
    task_identity: str
    stage: str
    status: str
    reason_code: str
    verifier_identity: str
    evidence_identity: str

    @classmethod
    def create(
        cls,
        *,
        request_identity: str,
        task_identity: str,
        stage: str,
        status: str,
        reason_code: str,
        verifier_identity: str,
        evidence_identity: str,
    ) -> "GateStageVerification":
        value: dict[str, object] = {
            "schema_version": GATE_STAGE_VERIFICATION_SCHEMA,
            "request_identity": request_identity,
            "task_identity": task_identity,
            "stage": stage,
            "status": status,
            "reason_code": reason_code,
            "verifier_identity": verifier_identity,
            "evidence_identity": evidence_identity,
        }
        value["identity"] = canonical_identity(
            value, scheme=GATE_STAGE_VERIFICATION_IDENTITY_SCHEME
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "GateStageVerification":
        code = "GATE_STAGE_VERIFICATION_INVALID"
        value = _document(raw, code=code)
        if (
            set(value)
            != {
                "schema_version",
                "identity",
                "request_identity",
                "task_identity",
                "stage",
                "status",
                "reason_code",
                "verifier_identity",
                "evidence_identity",
            }
            or value["schema_version"] != GATE_STAGE_VERIFICATION_SCHEMA
        ):
            raise fail(code, "Gate stage verification is invalid.")
        if (
            value["stage"] not in (*GATE_STAGES, "cleanup")
            or value["status"] not in ("passed", "failed")
            or not isinstance(value["reason_code"], str)
            or _SAFE_REASON.fullmatch(value["reason_code"]) is None
        ):
            raise fail(code, "Gate stage verification is invalid.")
        observed = _identity(value["identity"], code=code)
        request = _identity(value["request_identity"], code=code)
        task = _task_identity(value["task_identity"], code=code)
        verifier = _identity(value["verifier_identity"], code=code)
        evidence = _identity(value["evidence_identity"], code=code)
        if observed != canonical_identity(
            without_key(value, "identity"),
            scheme=GATE_STAGE_VERIFICATION_IDENTITY_SCHEME,
        ):
            raise fail(code, "Gate stage verification is invalid.")
        return cls(
            observed,
            request,
            task,
            value["stage"],
            value["status"],
            value["reason_code"],
            verifier,
            evidence,
        )

    def to_dict(self) -> dict[str, str]:
        return {
            "schema_version": GATE_STAGE_VERIFICATION_SCHEMA,
            "identity": self.identity,
            "request_identity": self.request_identity,
            "task_identity": self.task_identity,
            "stage": self.stage,
            "status": self.status,
            "reason_code": self.reason_code,
            "verifier_identity": self.verifier_identity,
            "evidence_identity": self.evidence_identity,
        }


class GateStageVerifier(Protocol):
    def verify(self, request: GateRequest, stage: str) -> GateStageVerification: ...


@dataclass(frozen=True)
class GateStageReceipt:
    identity: str
    request_identity: str
    task_identity: str
    stage: str
    status: str
    reason_code: str
    verification_identity: str
    verifier_identity: str
    evidence_identity: str

    @classmethod
    def from_verification(
        cls, request: GateRequest, verification: GateStageVerification
    ) -> "GateStageReceipt":
        if (
            verification.request_identity != request.identity
            or verification.task_identity != request.task_identity
        ):
            raise fail("GATE_STAGE_INVALID", "Deployment Gate stage is invalid.")
        value: dict[str, object] = {
            "schema_version": GATE_STAGE_RECEIPT_SCHEMA,
            "request_identity": request.identity,
            "task_identity": request.task_identity,
            "stage": verification.stage,
            "status": verification.status,
            "reason_code": verification.reason_code,
            "verification_identity": verification.identity,
            "verifier_identity": verification.verifier_identity,
            "evidence_identity": verification.evidence_identity,
        }
        value["identity"] = canonical_identity(
            value, scheme=GATE_STAGE_RECEIPT_IDENTITY_SCHEME
        )
        return cls(
            value["identity"],
            request.identity,
            request.task_identity,
            verification.stage,
            verification.status,
            verification.reason_code,
            verification.identity,
            verification.verifier_identity,
            verification.evidence_identity,
        )

    def to_dict(self) -> dict[str, str]:
        return {
            "schema_version": GATE_STAGE_RECEIPT_SCHEMA,
            "identity": self.identity,
            "request_identity": self.request_identity,
            "task_identity": self.task_identity,
            "stage": self.stage,
            "status": self.status,
            "reason_code": self.reason_code,
            "verification_identity": self.verification_identity,
            "verifier_identity": self.verifier_identity,
            "evidence_identity": self.evidence_identity,
        }


@dataclass(frozen=True)
class DeploymentGateRun:
    request: GateRequest
    state: str
    next_stage: int
    dispatch_identity: str | None
    approval_identity: str | None
    failure_stage: str | None
    receipts: tuple[GateStageReceipt, ...]

    @classmethod
    def prepared(cls, request: GateRequest) -> "DeploymentGateRun":
        return cls(
            GateRequest.from_dict(request.to_dict()),
            "prepared",
            0,
            None,
            None,
            None,
            (),
        )

    def dispatch(self, dispatch_identity: str) -> "DeploymentGateRun":
        if self.state != "prepared" or self.dispatch_identity is not None:
            raise fail(
                "GATE_TRANSITION_INVALID", "Deployment Gate transition is invalid."
            )
        return replace(
            self,
            state="dispatched",
            dispatch_identity=_identity(
                dispatch_identity, code="GATE_TRANSITION_INVALID"
            ),
        )

    def approve(self, approval_identity: str) -> "DeploymentGateRun":
        if self.state != "dispatched" or self.approval_identity is not None:
            raise fail(
                "GATE_TRANSITION_INVALID", "Deployment Gate transition is invalid."
            )
        approval = _identity(approval_identity, code="GATE_TRANSITION_INVALID")
        if approval == self.dispatch_identity:
            raise fail(
                "GATE_TRANSITION_INVALID", "Deployment Gate transition is invalid."
            )
        return replace(self, state="approved", approval_identity=approval)

    def start(self) -> "DeploymentGateRun":
        if self.state != "approved":
            raise fail(
                "GATE_TRANSITION_INVALID", "Deployment Gate transition is invalid."
            )
        return replace(self, state="running")

    def verify_stage(
        self, stage: str, *, verifier: GateStageVerifier
    ) -> "DeploymentGateRun":
        if (
            self.state != "running"
            or self.next_stage >= len(GATE_STAGES)
            or stage != GATE_STAGES[self.next_stage]
        ):
            raise fail(
                "GATE_TRANSITION_INVALID", "Deployment Gate transition is invalid."
            )
        verification = _run_stage_verifier(verifier, self.request, stage)
        receipt = GateStageReceipt.from_verification(self.request, verification)
        if verification.status == "failed":
            return replace(
                self,
                state="failed-awaiting-cleanup",
                failure_stage=stage,
                receipts=(*self.receipts, receipt),
            )
        next_stage = self.next_stage + 1
        return replace(
            self,
            state=("awaiting-cleanup" if next_stage == len(GATE_STAGES) else "running"),
            next_stage=next_stage,
            receipts=(*self.receipts, receipt),
        )

    def verify_cleanup(self, *, verifier: GateStageVerifier) -> "DeploymentGateRun":
        if self.state not in {"awaiting-cleanup", "failed-awaiting-cleanup"}:
            raise fail(
                "GATE_TRANSITION_INVALID", "Deployment Gate transition is invalid."
            )
        verification = _run_stage_verifier(verifier, self.request, "cleanup")
        receipt = GateStageReceipt.from_verification(self.request, verification)
        if verification.status == "failed":
            state = "cleanup-failed"
        elif self.state == "failed-awaiting-cleanup":
            state = "failed"
        else:
            state = "succeeded"
        return replace(self, state=state, receipts=(*self.receipts, receipt))


def _run_stage_verifier(
    verifier: GateStageVerifier, request: GateRequest, stage: str
) -> GateStageVerification:
    try:
        verification = GateStageVerification.from_dict(
            verifier.verify(request, stage).to_dict()
        )
    except DeploymentError:
        raise fail(
            "GATE_STAGE_VERIFICATION_FAILED", "Gate stage verification failed."
        ) from None
    except Exception:
        raise fail(
            "GATE_STAGE_VERIFICATION_FAILED", "Gate stage verification failed."
        ) from None
    if (
        verification.request_identity != request.identity
        or verification.task_identity != request.task_identity
        or verification.stage != stage
    ):
        raise fail(
            "GATE_STAGE_VERIFICATION_INVALID", "Gate stage verification is invalid."
        )
    return verification


def _observe_resources(
    policy: GatePolicy, task_identity: str
) -> tuple[ResourceEvidence, ...]:
    paths = {**policy.delete_paths(task_identity), **policy.retain_paths()}
    observed: list[ResourceEvidence] = []
    for kind in RESOURCE_KINDS:
        path = paths[kind]
        if kind == "dockerd-socket":
            value = _owned_socket(path, policy.operator_uid, policy.operator_gid)
        else:
            value = _owned_directory(path, policy.operator_uid, policy.operator_gid)
        observed.append(ResourceEvidence.create(kind=kind, path=path, observed=value))
    return tuple(observed)


def _observe_process(
    policy: GatePolicy, task_identity: str, name: str
) -> GateProcessIdentity:
    root = policy.task_root(task_identity)
    pid_content, _pid_stat = _read_owned_file(
        root / "pids" / f"{name}.pid",
        policy.operator_uid,
        policy.operator_gid,
        mode=0o600,
        max_bytes=32,
    )
    try:
        rendered = pid_content.decode("ascii")
        if not rendered.endswith("\n") or not rendered[:-1].isdigit():
            raise ValueError
        pid = int(rendered[:-1])
    except (UnicodeDecodeError, ValueError):
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.") from None
    if not 0 < pid <= 2**31 - 1:
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.")
    proc_root = Path("/proc") / str(pid)
    try:
        start_before = _process_start_ticks(proc_root / "stat")
        executable = (proc_root / "exe").stat()
        cmdline = _read_proc_file(proc_root / "cmdline", 64 * 1024)
        start_after = _process_start_ticks(proc_root / "stat")
    except OSError:
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.") from None
    if start_before != start_after or task_identity.encode(
        "ascii"
    ) not in cmdline.split(b"\0"):
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.")
    socket_device: int | None = None
    socket_inode: int | None = None
    if name != "runner":
        socket_path = root / (
            "redis/redis.sock" if name == "redis" else "docker/docker.sock"
        )
        socket_stat = _owned_socket(
            socket_path, policy.operator_uid, policy.operator_gid
        )
        socket_device, socket_inode = socket_stat.st_dev, socket_stat.st_ino
    return GateProcessIdentity.create(
        name=name,
        task_identity=task_identity,
        pid=pid,
        process_start_ticks=start_before,
        executable_device=executable.st_dev,
        executable_inode=executable.st_ino,
        cmdline_identity=f"sha256-{hashlib.sha256(cmdline).hexdigest()}",
        socket_device=socket_device,
        socket_inode=socket_inode,
    )


def _observe_cache(policy: GatePolicy, kind: str) -> CacheEvidence:
    cache_root = policy.gate_root / "cache" / kind
    _owned_directory(cache_root, policy.operator_uid, policy.operator_gid)
    content, _observed = _read_owned_file(
        cache_root / "active.identity",
        policy.operator_uid,
        policy.operator_gid,
        mode=0o600,
        max_bytes=80,
    )
    try:
        rendered = content.decode("ascii")
    except UnicodeDecodeError:
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.") from None
    if not rendered.endswith("\n"):
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.")
    identity = _identity(rendered[:-1], code="GATE_OBSERVATION_INVALID")
    object_path = cache_root / "objects" / identity
    object_stat = _owned_directory(
        object_path, policy.operator_uid, policy.operator_gid
    )
    return CacheEvidence(kind, identity, object_stat.st_dev, object_stat.st_ino)


def _observe_cleanup_executor(policy: GatePolicy) -> CleanupExecutorEvidence:
    content_hash, observed = _hash_owned_file(
        policy.cleanup_executor,
        policy.root_uid,
        policy.root_gid,
        mode=0o555,
        max_bytes=16 * 1024**2,
        parent_modes=(0o700, 0o755),
    )
    path_identity = canonical_identity(
        {"absolute_path": str(policy.cleanup_executor)},
        scheme="helixweave-deployment-gate-cleanup-executor-path-v1",
    )
    closure_identity, backend_hash, backend = _observe_cleanup_backend(policy)
    return CleanupExecutorEvidence(
        path_identity,
        f"sha256-{content_hash}",
        observed.st_dev,
        observed.st_ino,
        closure_identity,
        f"sha256-{backend_hash}",
        backend.st_dev,
        backend.st_ino,
    )


def _observe_cleanup_backend(
    policy: GatePolicy,
) -> tuple[str, str, os.stat_result]:
    root = _owned_directory(
        policy.operator_root,
        policy.root_uid,
        policy.root_gid,
        modes=(0o700, 0o755),
    )
    flags = (
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    descriptor: int | None = None
    try:
        descriptor = os.open(policy.operator_root, flags)
        if _inode(os.fstat(descriptor)) != _inode(root):
            raise OSError
        before = os.stat("current", dir_fd=descriptor, follow_symlinks=False)
        target = os.readlink("current", dir_fd=descriptor)
        after = os.stat("current", dir_fd=descriptor, follow_symlinks=False)
    except OSError:
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.") from None
    finally:
        if descriptor is not None:
            os.close(descriptor)
    if (
        not stat.S_ISLNK(before.st_mode)
        or _file_witness(before) != _file_witness(after)
        or before.st_uid != policy.root_uid
        or before.st_gid != policy.root_gid
        or _IDENTITY.fullmatch(target) is None
    ):
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.")
    closure = policy.operator_root / target
    _owned_directory(closure, policy.root_uid, policy.root_gid, modes=(0o555,))
    manifest_content, _manifest_stat = _read_owned_file(
        closure / "closure.json",
        policy.root_uid,
        policy.root_gid,
        mode=0o444,
        max_bytes=1024 * 1024,
        parent_modes=(0o555,),
    )
    try:
        manifest = json.loads(manifest_content, object_pairs_hook=_unique_object)
    except (UnicodeDecodeError, ValueError, json.JSONDecodeError):
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.") from None
    if (
        not isinstance(manifest, dict)
        or set(manifest) != {"schema_version", "identity", "files"}
        or manifest["schema_version"] != "helixweave-operator-closure-v1"
        or manifest["identity"] != target
        or not isinstance(manifest["files"], list)
        or canonical_json_bytes(manifest) != manifest_content
    ):
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.")
    identity_document = without_key(manifest, "identity")
    calculated_identity = (
        "sha256-"
        + hashlib.sha256(
            json.dumps(
                identity_document,
                allow_nan=False,
                separators=(",", ":"),
                sort_keys=True,
            ).encode("utf-8")
        ).hexdigest()
    )
    if calculated_identity != target:
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.")
    backend_hash, backend = _hash_owned_file(
        closure / "bin" / "helixweave-gate-cleanup",
        policy.root_uid,
        policy.root_gid,
        mode=0o555,
        max_bytes=1024 * 1024,
        parent_modes=(0o555,),
    )
    backend_records = [
        item
        for item in manifest["files"]
        if isinstance(item, dict) and item.get("path") == "bin/helixweave-gate-cleanup"
    ]
    if backend_records != [
        {
            "mode": 0o555,
            "path": "bin/helixweave-gate-cleanup",
            "sha256": backend_hash,
            "size_bytes": backend.st_size,
        }
    ]:
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.")
    return target, backend_hash, backend


def _runtime_map(raw: object) -> dict[str, str]:
    code = "GATE_OBSERVATION_INVALID"
    if not isinstance(raw, list):
        raise fail(code, "Gate observation is invalid.")
    records = tuple(RuntimeIdentity.from_dict(item, code=code) for item in raw)
    if tuple(item.component for item in records) != RUNTIME_COMPONENTS:
        raise fail(code, "Gate observation is invalid.")
    return {item.component: item.identity for item in records}


def _owned_directory(
    path: Path,
    uid: int,
    gid: int,
    *,
    modes: tuple[int, ...] = (0o700,),
) -> os.stat_result:
    selected = _absolute_path(path, code="GATE_OBSERVATION_INVALID", existing=True)
    try:
        observed = selected.lstat()
        flags = (
            os.O_RDONLY
            | getattr(os, "O_CLOEXEC", 0)
            | getattr(os, "O_DIRECTORY", 0)
            | getattr(os, "O_NOFOLLOW", 0)
        )
        descriptor = os.open(selected, flags)
    except OSError:
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.") from None
    try:
        anchored = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if (
        not stat.S_ISDIR(observed.st_mode)
        or stat.S_ISLNK(observed.st_mode)
        or observed.st_uid != uid
        or observed.st_gid != gid
        or stat.S_IMODE(observed.st_mode) not in modes
        or _inode(observed) != _inode(anchored)
    ):
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.")
    return anchored


def _owned_socket(path: Path, uid: int, gid: int) -> os.stat_result:
    selected = _absolute_path(path, code="GATE_OBSERVATION_INVALID", existing=True)
    parent = _owned_directory(selected.parent, uid, gid)
    directory_flags = (
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    descriptor: int | None = None
    try:
        descriptor = os.open(selected.parent, directory_flags)
        if _inode(os.fstat(descriptor)) != _inode(parent):
            raise OSError
        observed = os.stat(selected.name, dir_fd=descriptor, follow_symlinks=False)
    except OSError:
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.") from None
    finally:
        if descriptor is not None:
            os.close(descriptor)
    if (
        not stat.S_ISSOCK(observed.st_mode)
        or observed.st_uid != uid
        or observed.st_gid != gid
        or observed.st_nlink != 1
    ):
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.")
    return observed


def _read_canonical_json(path: Path, uid: int, gid: int) -> dict[str, object]:
    content, _observed = _read_owned_file(
        path, uid, gid, mode=0o600, max_bytes=_MAX_CONTROL_FILE_BYTES
    )
    try:
        raw = json.loads(content, object_pairs_hook=_unique_object)
    except (UnicodeDecodeError, ValueError, json.JSONDecodeError):
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.") from None
    if not isinstance(raw, dict) or canonical_json_bytes(raw) != content:
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.")
    return raw


def _read_owned_file(
    path: Path,
    uid: int,
    gid: int,
    *,
    mode: int,
    max_bytes: int,
    parent_modes: tuple[int, ...] = (0o700,),
) -> tuple[bytes, os.stat_result]:
    selected = _absolute_path(path, code="GATE_OBSERVATION_INVALID", existing=True)
    parent = _owned_directory(selected.parent, uid, gid, modes=parent_modes)
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
        if _inode(os.fstat(directory_descriptor)) != _inode(parent):
            raise OSError
        descriptor = os.open(selected.name, file_flags, dir_fd=directory_descriptor)
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_uid != uid
            or before.st_gid != gid
            or before.st_nlink != 1
            or stat.S_IMODE(before.st_mode) != mode
            or not 0 < before.st_size <= max_bytes
        ):
            raise OSError
        content = _read_fd(descriptor, max_bytes)
        after = os.fstat(descriptor)
        path_after = os.stat(
            selected.name, dir_fd=directory_descriptor, follow_symlinks=False
        )
        if (
            _file_witness(before) != _file_witness(after)
            or _inode(after) != _inode(path_after)
            or len(content) != before.st_size
        ):
            raise OSError
        return content, after
    except OSError:
        raise fail("GATE_OBSERVATION_INVALID", "Gate observation is invalid.") from None
    finally:
        if descriptor is not None:
            os.close(descriptor)
        if directory_descriptor is not None:
            os.close(directory_descriptor)


def _hash_owned_file(
    path: Path,
    uid: int,
    gid: int,
    *,
    mode: int,
    max_bytes: int,
    parent_modes: tuple[int, ...] = (0o700,),
) -> tuple[str, os.stat_result]:
    content, observed = _read_owned_file(
        path,
        uid,
        gid,
        mode=mode,
        max_bytes=max_bytes,
        parent_modes=parent_modes,
    )
    return hashlib.sha256(content).hexdigest(), observed


def _read_proc_file(path: Path, max_bytes: int) -> bytes:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags)
    try:
        return _read_fd(descriptor, max_bytes)
    finally:
        os.close(descriptor)


def _read_fd(descriptor: int, max_bytes: int) -> bytes:
    chunks: list[bytes] = []
    remaining = max_bytes + 1
    while remaining > 0:
        chunk = os.read(descriptor, min(1024 * 1024, remaining))
        if not chunk:
            break
        chunks.append(chunk)
        remaining -= len(chunk)
    content = b"".join(chunks)
    if len(content) > max_bytes:
        raise OSError
    return content


def _process_start_ticks(path: Path) -> int:
    content = _read_proc_file(path, 4096)
    closing = content.rfind(b")")
    if closing < 0:
        raise OSError
    fields = content[closing + 2 :].split()
    if len(fields) < 20:
        raise OSError
    return int(fields[19])


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
    "CLEANUP_EXECUTOR",
    "CLEANUP_PLAN_SCHEMA",
    "DELETE_RESOURCE_KINDS",
    "FIXED_DISK_REQUIRED_BYTES",
    "FIXED_LOAD_PER_CPU_MILLI",
    "GATE_REQUEST_SCHEMA",
    "GATE_STAGE_RECEIPT_SCHEMA",
    "GATE_STAGES",
    "REQUIRED_CACHE_KINDS",
    "RETAIN_RESOURCE_KINDS",
    "SUPPORTED_GATE_ROOT",
    "CacheEvidence",
    "CleanupExecutorEvidence",
    "CleanupPlan",
    "CleanupTarget",
    "DeploymentGateRun",
    "FilesystemGateObserver",
    "GateObservation",
    "GateObserver",
    "GatePolicy",
    "GateProcessIdentity",
    "GateRequest",
    "GateStageReceipt",
    "GateStageVerification",
    "GateStageVerifier",
    "ResourceEvidence",
    "RuntimeIdentity",
    "UnavailableGateObserver",
    "cleanup_script_identity",
    "prepare_gate_request",
    "render_cleanup_script",
    "verify_cleanup_script",
]
