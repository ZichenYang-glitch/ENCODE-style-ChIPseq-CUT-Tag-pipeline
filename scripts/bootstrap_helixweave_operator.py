#!/usr/bin/env python3
"""Interactively bootstrap the fixed HelixWeave privileged boundary.

This program never invokes sudo and never reads a password.  An administrator
must explicitly run it from a local terminal through interactive sudo.  It
installs only a root-owned helper, its sudoers policy, and inert templates; it
does not start services or alter Docker configuration.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import dataclass
import fcntl
import grp
import hashlib
import json
import os
from pathlib import Path
import pwd
import re
import stat
import subprocess
import sys
from typing import Protocol
import uuid


BOOTSTRAP_RECEIPT_SCHEMA = "helixweave-operator-bootstrap-receipt-v1"
SAFE_ENVIRONMENT = {
    "LANG": "C.UTF-8",
    "LC_ALL": "C.UTF-8",
    "PATH": "/usr/sbin:/usr/bin:/sbin:/bin",
}
OPERATOR_GROUP = "helixweave-operators"
SERVICE_USER = "helixweave"
SERVICE_GROUP = "helixweave"
API_USER = "helixweave-api"
API_GROUP = "helixweave-api"
CANDIDATE_USER = "helixweave-candidate"
CANDIDATE_GROUP = "helixweave-candidate"
GROUPADD = Path("/usr/sbin/groupadd")
USERMOD = Path("/usr/sbin/usermod")
USERADD = Path("/usr/sbin/useradd")
SYSTEMCTL = Path("/usr/bin/systemctl")
SYSTEMD_TMPFILES = Path("/usr/bin/systemd-tmpfiles")
VISUDO = Path("/usr/sbin/visudo")
MAX_TEMPLATE_BYTES = 1024 * 1024
_USERNAME = re.compile(r"^[a-z_][a-z0-9_-]{0,30}\$?$")
_CONTENT_IDENTITY = re.compile(r"^sha256-[0-9a-f]{64}$")
_UNINSTALL_PENDING_SUFFIX = "helixweave-uninstall-pending"


@dataclass(frozen=True)
class InstallSpec:
    source_name: str
    destination: Path
    mode: int
    sudoers: bool = False
    package_source: bool = False
    snapshot: bool = False


@dataclass(frozen=True)
class SystemAccountSpec:
    user: str
    group: str
    home: str


SYSTEM_ACCOUNTS = (
    SystemAccountSpec(SERVICE_USER, SERVICE_GROUP, "/var/lib/helixweave/service-home"),
    SystemAccountSpec(API_USER, API_GROUP, "/nonexistent"),
    SystemAccountSpec(CANDIDATE_USER, CANDIDATE_GROUP, "/nonexistent"),
)


OPERATOR_ROOT = Path("/opt/helixweave/operator")
BOUNDARY_SPECS = (
    InstallSpec(
        "helixweave-operator-launcher",
        Path("/usr/libexec/helixweave-operator"),
        0o555,
        snapshot=True,
    ),
    InstallSpec(
        "helixweave-operator.sudoers",
        Path("/etc/sudoers.d/helixweave-operator"),
        0o440,
        sudoers=True,
        snapshot=True,
    ),
    InstallSpec(
        "helixweave-gate-cleanup-launcher",
        Path("/usr/libexec/helixweave-gate-cleanup"),
        0o555,
        snapshot=True,
    ),
    InstallSpec(
        "helixweave-db-prepare-launcher",
        Path("/usr/libexec/helixweave-db-prepare"),
        0o555,
        snapshot=True,
    ),
    InstallSpec(
        "helixweave-operator-action",
        Path("/usr/libexec/helixweave-operator-action"),
        0o555,
        snapshot=True,
    ),
    InstallSpec(
        "helixweave-encode-runtime-prepare",
        Path("/usr/libexec/helixweave-encode-runtime-prepare"),
        0o555,
        snapshot=True,
    ),
    InstallSpec(
        "helixweave-bulk-runtime-prepare",
        Path("/usr/libexec/helixweave-bulk-runtime-prepare"),
        0o555,
        snapshot=True,
    ),
    InstallSpec(
        "helixweave-active-service",
        Path("/usr/libexec/helixweave-active-service"),
        0o555,
        snapshot=True,
    ),
    InstallSpec(
        "helixweave.tmpfiles.conf",
        Path("/usr/lib/tmpfiles.d/helixweave.conf"),
        0o444,
    ),
    InstallSpec(
        "helixweave-api.service.in",
        Path("/usr/lib/systemd/system/helixweave-api.service"),
        0o444,
    ),
    InstallSpec(
        "helixweave-worker.service.in",
        Path("/usr/lib/systemd/system/helixweave-worker.service"),
        0o444,
    ),
    InstallSpec(
        "helixweave-redis.service",
        Path("/usr/lib/systemd/system/helixweave-redis.service"),
        0o444,
    ),
    InstallSpec(
        "helixweave-docker-rootless.service",
        Path("/usr/lib/systemd/system/helixweave-docker-rootless.service"),
        0o444,
    ),
    InstallSpec(
        "helixweave-db-prepare.service.in",
        Path("/usr/lib/systemd/system/helixweave-db-prepare.service"),
        0o444,
    ),
    InstallSpec(
        "helixweave-operator-action.service",
        Path("/usr/lib/systemd/system/helixweave-operator-action.service"),
        0o444,
    ),
    InstallSpec(
        "helixweave-encode-runtime-prepare.service",
        Path("/usr/lib/systemd/system/helixweave-encode-runtime-prepare.service"),
        0o444,
    ),
    InstallSpec(
        "helixweave-bulk-runtime-prepare.service",
        Path("/usr/lib/systemd/system/helixweave-bulk-runtime-prepare.service"),
        0o444,
    ),
    InstallSpec(
        "helixweave.target",
        Path("/usr/lib/systemd/system/helixweave.target"),
        0o444,
    ),
    InstallSpec(
        "redis.conf",
        Path("/etc/helixweave/redis.conf"),
        0o444,
    ),
)
TEMPLATE_NAMES = (
    "helixweave.tmpfiles.conf",
    "helixweave-api.service.in",
    "helixweave-worker.service.in",
    "helixweave-redis.service",
    "helixweave-docker-rootless.service",
    "helixweave-db-prepare.service.in",
    "helixweave-operator-action.service",
    "helixweave-encode-runtime-prepare.service",
    "helixweave-bulk-runtime-prepare.service",
    "helixweave.target",
    "redis.conf",
)
BOUNDARY_SNAPSHOT_SPECS = tuple(
    InstallSpec(
        spec.source_name,
        Path("boundary") / spec.source_name,
        spec.mode,
    )
    for spec in BOUNDARY_SPECS
    if spec.snapshot
)
SNAPSHOTTED_BOUNDARIES_BY_SOURCE = {
    spec.source_name: spec for spec in BOUNDARY_SPECS if spec.snapshot
}
STABLE_BOUNDARY_SPECS = tuple(spec for spec in BOUNDARY_SPECS if spec.snapshot)
LINKED_BOUNDARY_SPECS = tuple(spec for spec in BOUNDARY_SPECS if not spec.snapshot)
BOUNDARY_SNAPSHOT_PATHS = {
    **{
        spec.source_name: Path("boundary") / spec.source_name
        for spec in BOUNDARY_SPECS
        if spec.snapshot
    },
    **{name: Path("templates") / name for name in TEMPLATE_NAMES},
}
CLOSURE_SPECS = (
    InstallSpec("helixweave-operator", Path("bin/helixweave-operator"), 0o555),
    InstallSpec(
        "helixweave-gate-cleanup",
        Path("bin/helixweave-gate-cleanup"),
        0o555,
    ),
    InstallSpec(
        "helixweave-db-prepare",
        Path("bin/helixweave-db-prepare"),
        0o555,
    ),
    InstallSpec(
        "operator-package-init.py",
        Path("lib/encode_pipeline/__init__.py"),
        0o444,
    ),
    InstallSpec(
        "operator-deployment-init.py",
        Path("lib/encode_pipeline/deployment/__init__.py"),
        0o444,
    ),
    *tuple(
        InstallSpec(
            name,
            Path("lib/encode_pipeline/deployment") / name,
            0o444,
            package_source=True,
        )
        for name in (
            "bulk_docker_boundary.py",
            "bundle.py",
            "canonical.py",
            "database.py",
            "errors.py",
            "filesystem.py",
            "layout.py",
            "models.py",
            "operator.py",
            "operator_action.py",
            "operator_boundary.py",
            "operator_transaction.py",
            "state.py",
        )
    ),
    *BOUNDARY_SNAPSHOT_SPECS,
    *tuple(
        InstallSpec(name, Path("templates") / name, 0o444) for name in TEMPLATE_NAMES
    ),
)
# Kept as the complete fixed installation inventory for contract tests.
INSTALL_SPECS = (*BOUNDARY_SPECS, *CLOSURE_SPECS)


class BootstrapFailure(RuntimeError):
    """Private failure with a stable, path-free public code."""

    def __init__(self, code: str, message: str) -> None:
        self.code = code
        self.message = message
        super().__init__(code)


@dataclass(frozen=True)
class BootstrapResult:
    installed_files: int
    closure_identity: str


class BootstrapBackend(Protocol):
    def apply(self, *, operation: str, invoking_user: str) -> BootstrapResult: ...


CommandRunner = Callable[[Sequence[str]], bool]
SudoersValidator = Callable[[Path], bool]


class HostBootstrapBackend:
    """Root-only, atomic filesystem installer with fixed account commands."""

    def __init__(
        self,
        *,
        source_root: Path | None = None,
        package_source_root: Path | None = None,
        root_prefix: Path = Path("/"),
        owner_uid: int = 0,
        owner_gid: int = 0,
        command_runner: CommandRunner | None = None,
        sudoers_validator: SudoersValidator | None = None,
    ) -> None:
        self.source_root = _template_root() if source_root is None else source_root
        self.package_source_root = (
            _template_root().parent
            if package_source_root is None
            else package_source_root
        )
        self.root_prefix = root_prefix
        self.owner_uid = owner_uid
        self.owner_gid = owner_gid
        self.command_runner = (
            _run_fixed_command if command_runner is None else command_runner
        )
        self.sudoers_validator = (
            _validate_sudoers if sudoers_validator is None else sudoers_validator
        )
        if (
            not self.source_root.is_absolute()
            or not self.package_source_root.is_absolute()
            or not self.root_prefix.is_absolute()
        ):
            raise BootstrapFailure(
                "BOOTSTRAP_SOURCE_INVALID",
                "Operator bootstrap source is invalid.",
            )

    def apply(self, *, operation: str, invoking_user: str) -> BootstrapResult:
        if operation not in {"install", "update"}:
            raise BootstrapFailure(
                "BOOTSTRAP_REQUEST_INVALID", "Operator bootstrap request is invalid."
            )
        if not _USERNAME.fullmatch(invoking_user):
            raise BootstrapFailure(
                "BOOTSTRAP_CALLER_INVALID", "Interactive sudo caller is invalid."
            )
        if self.root_prefix == Path("/"):
            self._ensure_service_accounts()
            self._ensure_operator_group(invoking_user)
        source_content = {
            _source_key(spec): _read_template(self._source_path(spec))
            for spec in INSTALL_SPECS
        }
        with self._bootstrap_lock():
            repair_target = self._require_operation_state(operation)
            if repair_target is not None:
                self._repair_selected_boundaries(repair_target)
                if operation == "install":
                    if self.root_prefix == Path("/"):
                        self._activate_static_boundaries()
                    installed = len(CLOSURE_SPECS) + len(BOUNDARY_SPECS) + 1
                    return BootstrapResult(installed, repair_target)
            closure_identity = self._install_closure(source_content)
            if operation == "install":
                for spec in STABLE_BOUNDARY_SPECS:
                    validator = self.sudoers_validator if spec.sudoers else None
                    self._install_file(
                        self._destination(spec.destination),
                        source_content[_source_key(spec)],
                        mode=spec.mode,
                        validator=validator,
                    )
                for spec in LINKED_BOUNDARY_SPECS:
                    self._install_boundary_link(spec)
            self._require_stable_boundaries(closure_identity)
            self._require_linked_boundaries()
            self._switch_closure(closure_identity)
            if self.root_prefix == Path("/"):
                self._activate_static_boundaries()
        installed = len(CLOSURE_SPECS) + len(BOUNDARY_SPECS) + 1
        return BootstrapResult(installed, closure_identity)

    def _require_operation_state(self, operation: str) -> str | None:
        current = self._destination(OPERATOR_ROOT) / "current"
        installed = current.exists() or current.is_symlink()
        if operation == "install" and installed:
            target = self._selected_closure_identity()
            try:
                self._require_stable_boundaries(target)
                self._require_linked_boundaries()
            except BootstrapFailure as error:
                if error.code != "BOOTSTRAP_BOUNDARY_UPDATE_REQUIRED":
                    raise
                return target
            raise BootstrapFailure(
                "BOOTSTRAP_ALREADY_INSTALLED",
                "Operator boundary is already installed; use the update operation.",
            )
        if operation == "update" and not installed:
            raise BootstrapFailure(
                "BOOTSTRAP_NOT_INSTALLED",
                "Operator boundary is not installed.",
            )
        if operation == "update":
            target = self._selected_closure_identity()
            try:
                self._require_stable_boundaries(target)
                self._require_linked_boundaries()
            except BootstrapFailure as error:
                if error.code != "BOOTSTRAP_BOUNDARY_UPDATE_REQUIRED":
                    raise
                return target
        return None

    def _selected_closure_identity(self) -> str:
        current = self._destination(OPERATOR_ROOT) / "current"
        try:
            observed = current.lstat()
            target = os.readlink(current)
            closure = current.parent / target
            closure_stat = closure.lstat()
            boundary_stat = (closure / "boundary").lstat()
        except OSError:
            raise BootstrapFailure(
                "BOOTSTRAP_BOUNDARY_UPDATE_REQUIRED",
                "Stable operator boundary requires a reviewed bootstrap upgrade.",
            ) from None
        if (
            not stat.S_ISLNK(observed.st_mode)
            or observed.st_uid != self.owner_uid
            or observed.st_gid != self.owner_gid
            or _CONTENT_IDENTITY.fullmatch(target) is None
            or not stat.S_ISDIR(closure_stat.st_mode)
            or stat.S_ISLNK(closure_stat.st_mode)
            or closure_stat.st_uid != self.owner_uid
            or closure_stat.st_gid != self.owner_gid
            or stat.S_IMODE(closure_stat.st_mode) != 0o555
            or not stat.S_ISDIR(boundary_stat.st_mode)
            or stat.S_ISLNK(boundary_stat.st_mode)
            or boundary_stat.st_uid != self.owner_uid
            or boundary_stat.st_gid != self.owner_gid
            or stat.S_IMODE(boundary_stat.st_mode) != 0o555
        ):
            raise BootstrapFailure(
                "BOOTSTRAP_BOUNDARY_UPDATE_REQUIRED",
                "Stable operator boundary requires a reviewed bootstrap upgrade.",
            )
        return target

    def _repair_selected_boundaries(self, closure_identity: str) -> None:
        closure = self._destination(OPERATOR_ROOT) / closure_identity
        try:
            for spec in STABLE_BOUNDARY_SPECS:
                snapshot = closure / BOUNDARY_SNAPSHOT_PATHS[spec.source_name]
                content = _read_template(snapshot)
                self._verify_installed_file(snapshot, content, spec.mode)
                destination = self._destination(spec.destination)
                self._restore_pending_boundary(
                    destination,
                    expected_content=content,
                    expected_mode=spec.mode,
                )
                validator = self.sudoers_validator if spec.sudoers else None
                self._install_file(
                    destination,
                    content,
                    mode=spec.mode,
                    validator=validator,
                )
            for spec in LINKED_BOUNDARY_SPECS:
                self._restore_pending_boundary(
                    self._destination(spec.destination),
                    expected_target=self._boundary_link_target(spec).as_posix(),
                )
                self._install_boundary_link(spec)
            self._require_stable_boundaries(closure_identity)
            self._require_linked_boundaries()
        except (BootstrapFailure, KeyError):
            raise BootstrapFailure(
                "BOOTSTRAP_BOUNDARY_UPDATE_REQUIRED",
                "Stable operator boundary requires a reviewed bootstrap upgrade.",
            ) from None

    def _restore_pending_boundary(
        self,
        destination: Path,
        *,
        expected_content: bytes | None = None,
        expected_mode: int | None = None,
        expected_target: str | None = None,
    ) -> None:
        pending = destination.with_name(
            f".{destination.name}.{_UNINSTALL_PENDING_SUFFIX}"
        )
        try:
            pending_observed = pending.lstat()
        except FileNotFoundError:
            return
        except OSError:
            raise BootstrapFailure(
                "BOOTSTRAP_DESTINATION_UNSAFE",
                "Existing operator boundary is unsafe.",
            ) from None
        if destination.exists() or destination.is_symlink():
            raise BootstrapFailure(
                "BOOTSTRAP_DESTINATION_UNSAFE",
                "Existing operator boundary is unsafe.",
            )
        if (
            pending_observed.st_uid != self.owner_uid
            or pending_observed.st_gid != self.owner_gid
            or pending_observed.st_nlink != 1
        ):
            raise BootstrapFailure(
                "BOOTSTRAP_DESTINATION_UNSAFE",
                "Existing operator boundary is unsafe.",
            )
        if expected_target is not None:
            try:
                target = os.readlink(pending)
            except OSError:
                raise BootstrapFailure(
                    "BOOTSTRAP_DESTINATION_UNSAFE",
                    "Existing operator boundary is unsafe.",
                ) from None
            if (
                expected_content is not None
                or expected_mode is not None
                or not stat.S_ISLNK(pending_observed.st_mode)
                or target != expected_target
            ):
                raise BootstrapFailure(
                    "BOOTSTRAP_DESTINATION_UNSAFE",
                    "Existing operator boundary is unsafe.",
                )
        elif (
            expected_content is None
            or expected_mode is None
            or not stat.S_ISREG(pending_observed.st_mode)
            or stat.S_ISLNK(pending_observed.st_mode)
            or stat.S_IMODE(pending_observed.st_mode) != expected_mode
            or _read_template(pending) != expected_content
        ):
            raise BootstrapFailure(
                "BOOTSTRAP_DESTINATION_UNSAFE",
                "Existing operator boundary is unsafe.",
            )
        try:
            os.replace(pending, destination)
            _fsync_directory(destination.parent)
        except OSError:
            raise BootstrapFailure(
                "BOOTSTRAP_INSTALL_FAILED",
                "Operator boundary could not be installed.",
            ) from None

    def _require_stable_boundaries(
        self,
        closure_identity: str,
    ) -> None:
        if _CONTENT_IDENTITY.fullmatch(closure_identity) is None:
            raise BootstrapFailure(
                "BOOTSTRAP_RESULT_INVALID", "Operator bootstrap result is invalid."
            )
        closure = self._destination(OPERATOR_ROOT) / closure_identity
        try:
            for spec in STABLE_BOUNDARY_SPECS:
                snapshot = closure / BOUNDARY_SNAPSHOT_PATHS[spec.source_name]
                expected = _read_template(snapshot)
                self._verify_installed_file(
                    snapshot,
                    expected,
                    spec.mode,
                )
                self._verify_installed_file(
                    self._destination(spec.destination),
                    expected,
                    spec.mode,
                )
        except (BootstrapFailure, KeyError):
            raise BootstrapFailure(
                "BOOTSTRAP_BOUNDARY_UPDATE_REQUIRED",
                "Stable operator boundary requires a reviewed bootstrap upgrade.",
            ) from None

    def verify_selected_boundaries(self) -> str:
        """Verify installed boundaries against the atomically selected snapshot."""
        current = self._destination(OPERATOR_ROOT) / "current"
        try:
            observed = current.lstat()
            target = os.readlink(current)
        except OSError:
            raise BootstrapFailure(
                "BOOTSTRAP_BOUNDARY_UPDATE_REQUIRED",
                "Stable operator boundary requires a reviewed bootstrap upgrade.",
            ) from None
        if (
            not stat.S_ISLNK(observed.st_mode)
            or observed.st_uid != self.owner_uid
            or observed.st_gid != self.owner_gid
            or _CONTENT_IDENTITY.fullmatch(target) is None
        ):
            raise BootstrapFailure(
                "BOOTSTRAP_BOUNDARY_UPDATE_REQUIRED",
                "Stable operator boundary requires a reviewed bootstrap upgrade.",
            )
        self._require_stable_boundaries(target)
        self._require_linked_boundaries()
        return target

    def _require_linked_boundaries(self) -> None:
        try:
            for spec in LINKED_BOUNDARY_SPECS:
                path = self._destination(spec.destination)
                observed = path.lstat()
                target = os.readlink(path)
                if (
                    not stat.S_ISLNK(observed.st_mode)
                    or observed.st_uid != self.owner_uid
                    or observed.st_gid != self.owner_gid
                    or target != self._boundary_link_target(spec).as_posix()
                ):
                    raise OSError
        except OSError:
            raise BootstrapFailure(
                "BOOTSTRAP_BOUNDARY_UPDATE_REQUIRED",
                "Stable operator boundary requires a reviewed bootstrap upgrade.",
            ) from None

    @contextmanager
    def _bootstrap_lock(self):
        operator_root = self._destination(OPERATOR_ROOT)
        self._ensure_directory(operator_root)
        lock_path = operator_root / "bootstrap.lock"
        flags = (
            os.O_RDWR
            | os.O_CREAT
            | getattr(os, "O_CLOEXEC", 0)
            | getattr(os, "O_NOFOLLOW", 0)
        )
        try:
            descriptor = os.open(lock_path, flags, 0o600)
        except OSError:
            raise BootstrapFailure(
                "BOOTSTRAP_INSTALL_FAILED",
                "Operator boundary could not be installed.",
            ) from None
        try:
            observed = os.fstat(descriptor)
            if (
                not stat.S_ISREG(observed.st_mode)
                or observed.st_nlink != 1
                or observed.st_uid != self.owner_uid
                or observed.st_gid != self.owner_gid
                or stat.S_IMODE(observed.st_mode) != 0o600
            ):
                raise BootstrapFailure(
                    "BOOTSTRAP_DESTINATION_UNSAFE",
                    "Existing operator boundary is unsafe.",
                )
            try:
                fcntl.flock(descriptor, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise BootstrapFailure(
                    "BOOTSTRAP_BUSY",
                    "Another operator bootstrap is in progress.",
                ) from None
            yield
        finally:
            try:
                fcntl.flock(descriptor, fcntl.LOCK_UN)
            finally:
                os.close(descriptor)

    def _install_closure(self, source_content: Mapping[str, bytes]) -> str:
        records: list[dict[str, object]] = []
        for spec in CLOSURE_SPECS:
            content = source_content[_source_key(spec)]
            record: dict[str, object] = {
                "mode": spec.mode,
                "path": spec.destination.as_posix(),
                "sha256": hashlib.sha256(content).hexdigest(),
                "size_bytes": len(content),
            }
            boundary = SNAPSHOTTED_BOUNDARIES_BY_SOURCE.get(spec.source_name)
            if boundary is not None and spec.destination == (
                Path("boundary") / spec.source_name
            ):
                record["installed_path"] = boundary.destination.as_posix()
            records.append(record)
        records.sort(key=lambda item: str(item["path"]))
        identity_document = {
            "schema_version": "helixweave-operator-closure-v1",
            "files": records,
        }
        identity = (
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
        manifest = {
            **identity_document,
            "identity": identity,
        }
        manifest_bytes = (
            json.dumps(
                manifest,
                allow_nan=False,
                separators=(",", ":"),
                sort_keys=True,
            )
            + "\n"
        ).encode("utf-8")
        operator_root = self._destination(OPERATOR_ROOT)
        destination = operator_root / identity
        if destination.exists() or destination.is_symlink():
            self._verify_closure(destination, manifest_bytes, source_content)
            return identity

        partial = operator_root / f".partial-{identity}-{uuid.uuid4().hex}"
        try:
            partial.mkdir(mode=0o700)
            for spec in CLOSURE_SPECS:
                self._install_file(
                    partial / spec.destination,
                    source_content[_source_key(spec)],
                    mode=spec.mode,
                    validator=None,
                )
            self._install_file(
                partial / "closure.json",
                manifest_bytes,
                mode=0o444,
                validator=None,
            )
            directories = [item for item in partial.rglob("*") if item.is_dir()]
            for directory in sorted(
                directories,
                key=lambda item: len(item.parts),
                reverse=True,
            ):
                os.chmod(directory, 0o555)
                _fsync_directory(directory)
            os.chmod(partial, 0o555)
            _fsync_directory(partial)
            os.replace(partial, destination)
            _fsync_directory(operator_root)
        except BootstrapFailure:
            raise
        except OSError:
            raise BootstrapFailure(
                "BOOTSTRAP_INSTALL_FAILED",
                "Operator boundary could not be installed.",
            ) from None
        self._verify_closure(destination, manifest_bytes, source_content)
        return identity

    def _verify_closure(
        self,
        root: Path,
        manifest_bytes: bytes,
        source_content: Mapping[str, bytes],
    ) -> None:
        try:
            observed_root = root.lstat()
        except OSError:
            raise BootstrapFailure(
                "BOOTSTRAP_DESTINATION_UNSAFE",
                "Existing operator boundary is unsafe.",
            ) from None
        if (
            not stat.S_ISDIR(observed_root.st_mode)
            or stat.S_ISLNK(observed_root.st_mode)
            or observed_root.st_uid != self.owner_uid
            or observed_root.st_gid != self.owner_gid
            or stat.S_IMODE(observed_root.st_mode) != 0o555
        ):
            raise BootstrapFailure(
                "BOOTSTRAP_DESTINATION_UNSAFE",
                "Existing operator boundary is unsafe.",
            )
        expected = {"closure.json"}
        self._verify_installed_file(root / "closure.json", manifest_bytes, 0o444)
        for spec in CLOSURE_SPECS:
            expected.add(spec.destination.as_posix())
            self._verify_installed_file(
                root / spec.destination,
                source_content[_source_key(spec)],
                spec.mode,
            )
        observed_files: set[str] = set()
        for path in root.rglob("*"):
            item = path.lstat()
            if stat.S_ISLNK(item.st_mode):
                raise BootstrapFailure(
                    "BOOTSTRAP_DESTINATION_UNSAFE",
                    "Existing operator boundary is unsafe.",
                )
            if stat.S_ISDIR(item.st_mode):
                if (
                    item.st_uid != self.owner_uid
                    or item.st_gid != self.owner_gid
                    or stat.S_IMODE(item.st_mode) != 0o555
                ):
                    raise BootstrapFailure(
                        "BOOTSTRAP_DESTINATION_UNSAFE",
                        "Existing operator boundary is unsafe.",
                    )
                continue
            observed_files.add(path.relative_to(root).as_posix())
        if observed_files != expected:
            raise BootstrapFailure(
                "BOOTSTRAP_DESTINATION_UNSAFE",
                "Existing operator boundary is unsafe.",
            )

    def _verify_installed_file(
        self,
        path: Path,
        expected: bytes,
        mode: int,
    ) -> None:
        content = _read_template(path)
        observed = path.stat(follow_symlinks=False)
        if (
            content != expected
            or observed.st_uid != self.owner_uid
            or observed.st_gid != self.owner_gid
            or stat.S_IMODE(observed.st_mode) != mode
        ):
            raise BootstrapFailure(
                "BOOTSTRAP_DESTINATION_UNSAFE",
                "Existing operator boundary is unsafe.",
            )

    def _switch_closure(self, identity: str) -> None:
        if _CONTENT_IDENTITY.fullmatch(identity) is None:
            raise BootstrapFailure(
                "BOOTSTRAP_RESULT_INVALID", "Operator bootstrap result is invalid."
            )
        operator_root = self._destination(OPERATOR_ROOT)
        current = operator_root / "current"
        temporary = operator_root / f".current-{uuid.uuid4().hex}"
        try:
            if current.exists() or current.is_symlink():
                observed = current.lstat()
                target = os.readlink(current)
                if (
                    not stat.S_ISLNK(observed.st_mode)
                    or observed.st_uid != self.owner_uid
                    or observed.st_gid != self.owner_gid
                    or _CONTENT_IDENTITY.fullmatch(target) is None
                ):
                    raise BootstrapFailure(
                        "BOOTSTRAP_DESTINATION_UNSAFE",
                        "Existing operator boundary is unsafe.",
                    )
            os.symlink(identity, temporary)
            os.replace(temporary, current)
            _fsync_directory(operator_root)
        except BootstrapFailure:
            raise
        except OSError:
            raise BootstrapFailure(
                "BOOTSTRAP_INSTALL_FAILED",
                "Operator boundary could not be installed.",
            ) from None
        finally:
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass
            except OSError:
                raise BootstrapFailure(
                    "BOOTSTRAP_INSTALL_FAILED",
                    "Operator boundary could not be installed.",
                ) from None

    def _destination(self, absolute: Path) -> Path:
        if self.root_prefix == Path("/"):
            return absolute
        return self.root_prefix.joinpath(*absolute.parts[1:])

    def _source_path(self, spec: InstallSpec) -> Path:
        root = self.package_source_root if spec.package_source else self.source_root
        return root / spec.source_name

    @staticmethod
    def _boundary_link_target(spec: InstallSpec) -> Path:
        snapshot = BOUNDARY_SNAPSHOT_PATHS.get(spec.source_name)
        if spec not in LINKED_BOUNDARY_SPECS or snapshot is None:
            raise BootstrapFailure(
                "BOOTSTRAP_INSTALL_FAILED",
                "Operator boundary could not be installed.",
            )
        return OPERATOR_ROOT / "current" / snapshot

    def _install_boundary_link(self, spec: InstallSpec) -> None:
        destination = self._destination(spec.destination)
        target = self._boundary_link_target(spec).as_posix()
        self._ensure_directory(destination.parent)
        try:
            observed = destination.lstat()
        except FileNotFoundError:
            observed = None
        except OSError:
            raise BootstrapFailure(
                "BOOTSTRAP_INSTALL_FAILED",
                "Operator boundary could not be installed.",
            ) from None
        if observed is not None:
            try:
                installed_target = os.readlink(destination)
            except OSError:
                raise BootstrapFailure(
                    "BOOTSTRAP_DESTINATION_UNSAFE",
                    "Existing operator boundary is unsafe.",
                ) from None
            if (
                not stat.S_ISLNK(observed.st_mode)
                or observed.st_uid != self.owner_uid
                or observed.st_gid != self.owner_gid
                or installed_target != target
            ):
                raise BootstrapFailure(
                    "BOOTSTRAP_DESTINATION_UNSAFE",
                    "Existing operator boundary is unsafe.",
                )
            return
        temporary = destination.parent / f".{destination.name}.{uuid.uuid4().hex}.tmp"
        try:
            os.symlink(target, temporary)
            os.lchown(temporary, self.owner_uid, self.owner_gid)
            os.replace(temporary, destination)
            _fsync_directory(destination.parent)
        except OSError:
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass
            raise BootstrapFailure(
                "BOOTSTRAP_INSTALL_FAILED",
                "Operator boundary could not be installed.",
            ) from None

    def _ensure_operator_group(self, invoking_user: str) -> None:
        try:
            grp.getgrnam(OPERATOR_GROUP)
        except KeyError:
            if not self.command_runner((str(GROUPADD), "--system", OPERATOR_GROUP)):
                raise BootstrapFailure(
                    "BOOTSTRAP_ACCOUNT_FAILED",
                    "Operator account boundary could not be installed.",
                )
        if not self.command_runner(
            (
                str(USERMOD),
                "--append",
                "--groups",
                OPERATOR_GROUP,
                "--",
                invoking_user,
            )
        ):
            raise BootstrapFailure(
                "BOOTSTRAP_ACCOUNT_FAILED",
                "Operator account boundary could not be installed.",
            )

    def _ensure_service_accounts(self) -> None:
        accounts = tuple(self._ensure_system_account(spec) for spec in SYSTEM_ACCOUNTS)
        if len({account.pw_uid for account, _group in accounts}) != len(
            accounts
        ) or len({group.gr_gid for _account, group in accounts}) != len(accounts):
            raise BootstrapFailure(
                "BOOTSTRAP_ACCOUNT_FAILED",
                "Operator account boundary could not be installed.",
            )

    def _ensure_system_account(self, spec: SystemAccountSpec):
        try:
            group = grp.getgrnam(spec.group)
        except KeyError:
            if not self.command_runner((str(GROUPADD), "--system", spec.group)):
                raise BootstrapFailure(
                    "BOOTSTRAP_ACCOUNT_FAILED",
                    "Operator account boundary could not be installed.",
                )
            try:
                group = grp.getgrnam(spec.group)
            except KeyError:
                raise BootstrapFailure(
                    "BOOTSTRAP_ACCOUNT_FAILED",
                    "Operator account boundary could not be installed.",
                ) from None
        if (
            group.gr_name != spec.group
            or group.gr_gid <= 0
            or set(group.gr_mem) - {spec.user}
        ):
            raise BootstrapFailure(
                "BOOTSTRAP_ACCOUNT_FAILED",
                "Operator account boundary could not be installed.",
            )
        try:
            account = pwd.getpwnam(spec.user)
        except KeyError:
            if not self.command_runner(_useradd_command(spec)):
                raise BootstrapFailure(
                    "BOOTSTRAP_ACCOUNT_FAILED",
                    "Operator account boundary could not be installed.",
                )
            try:
                account = pwd.getpwnam(spec.user)
            except KeyError:
                raise BootstrapFailure(
                    "BOOTSTRAP_ACCOUNT_FAILED",
                    "Operator account boundary could not be installed.",
                ) from None
        if (
            account.pw_name != spec.user
            or account.pw_uid <= 0
            or account.pw_dir != spec.home
            or account.pw_shell != "/usr/sbin/nologin"
            or account.pw_gid != group.gr_gid
        ):
            raise BootstrapFailure(
                "BOOTSTRAP_ACCOUNT_FAILED",
                "Operator account boundary could not be installed.",
            )
        try:
            memberships = set(os.getgrouplist(spec.user, group.gr_gid))
        except (KeyError, OSError):
            raise BootstrapFailure(
                "BOOTSTRAP_ACCOUNT_FAILED",
                "Operator account boundary could not be installed.",
            ) from None
        if memberships != {group.gr_gid}:
            raise BootstrapFailure(
                "BOOTSTRAP_ACCOUNT_FAILED",
                "Operator account boundary could not be installed.",
            )
        return account, group

    def _activate_static_boundaries(self) -> None:
        commands = (
            (
                str(SYSTEMD_TMPFILES),
                "--create",
                "/usr/lib/tmpfiles.d/helixweave.conf",
            ),
            (str(SYSTEMCTL), "daemon-reload"),
        )
        if any(not self.command_runner(command) for command in commands):
            raise BootstrapFailure(
                "BOOTSTRAP_INSTALL_FAILED",
                "Operator boundary could not be installed.",
            )

    def _install_file(
        self,
        destination: Path,
        content: bytes,
        *,
        mode: int,
        validator: SudoersValidator | None,
    ) -> None:
        self._ensure_directory(destination.parent)
        try:
            observed = destination.lstat()
        except FileNotFoundError:
            observed = None
        except OSError:
            raise BootstrapFailure(
                "BOOTSTRAP_INSTALL_FAILED",
                "Operator boundary could not be installed.",
            ) from None
        if observed is not None and (
            not stat.S_ISREG(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or observed.st_nlink != 1
            or observed.st_uid != self.owner_uid
            or observed.st_gid != self.owner_gid
            or stat.S_IMODE(observed.st_mode) & 0o022
        ):
            raise BootstrapFailure(
                "BOOTSTRAP_DESTINATION_UNSAFE",
                "Existing operator boundary is unsafe.",
            )

        temporary = destination.parent / f".{destination.name}.{uuid.uuid4().hex}.tmp"
        flags = (
            os.O_WRONLY
            | os.O_CREAT
            | os.O_EXCL
            | getattr(os, "O_CLOEXEC", 0)
            | getattr(os, "O_NOFOLLOW", 0)
        )
        descriptor = -1
        try:
            descriptor = os.open(temporary, flags, 0o600)
            _write_all(descriptor, content)
            os.fchmod(descriptor, mode)
            os.fchown(descriptor, self.owner_uid, self.owner_gid)
            os.fsync(descriptor)
            os.close(descriptor)
            descriptor = -1
            if validator is not None and not validator(temporary):
                raise BootstrapFailure(
                    "BOOTSTRAP_SUDOERS_INVALID",
                    "Operator sudo policy is invalid.",
                )
            os.replace(temporary, destination)
            _fsync_directory(destination.parent)
        except BootstrapFailure:
            raise
        except OSError:
            raise BootstrapFailure(
                "BOOTSTRAP_INSTALL_FAILED",
                "Operator boundary could not be installed.",
            ) from None
        finally:
            if descriptor >= 0:
                os.close(descriptor)
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass
            except OSError:
                raise BootstrapFailure(
                    "BOOTSTRAP_INSTALL_FAILED",
                    "Operator boundary could not be installed.",
                ) from None

    def _ensure_directory(self, path: Path) -> None:
        current = self.root_prefix
        relative = path.relative_to(self.root_prefix)
        for part in relative.parts:
            current = current / part
            try:
                observed = current.lstat()
            except FileNotFoundError:
                try:
                    current.mkdir(mode=0o755)
                    os.chown(current, self.owner_uid, self.owner_gid)
                except OSError:
                    raise BootstrapFailure(
                        "BOOTSTRAP_INSTALL_FAILED",
                        "Operator boundary could not be installed.",
                    ) from None
                continue
            except OSError:
                raise BootstrapFailure(
                    "BOOTSTRAP_INSTALL_FAILED",
                    "Operator boundary could not be installed.",
                ) from None
            if (
                not stat.S_ISDIR(observed.st_mode)
                or stat.S_ISLNK(observed.st_mode)
                or observed.st_uid != self.owner_uid
                or observed.st_gid != self.owner_gid
                or stat.S_IMODE(observed.st_mode) & 0o022
            ):
                raise BootstrapFailure(
                    "BOOTSTRAP_DESTINATION_UNSAFE",
                    "Existing operator boundary is unsafe.",
                )


def _template_root() -> Path:
    checkout_root = Path(__file__).resolve().parents[1]
    candidate = checkout_root / "src" / "encode_pipeline" / "deployment" / "templates"
    if candidate.is_dir() and not candidate.is_symlink():
        return candidate
    try:
        from importlib.resources import files

        installed = files("encode_pipeline.deployment.templates")
        installed_path = Path(str(installed))
    except (ImportError, TypeError):
        raise BootstrapFailure(
            "BOOTSTRAP_SOURCE_INVALID", "Operator bootstrap source is invalid."
        ) from None
    if not installed_path.is_dir() or installed_path.is_symlink():
        raise BootstrapFailure(
            "BOOTSTRAP_SOURCE_INVALID", "Operator bootstrap source is invalid."
        )
    return installed_path


def _source_key(spec: InstallSpec) -> str:
    return f"package:{spec.source_name}" if spec.package_source else spec.source_name


def _read_template(path: Path) -> bytes:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError:
        raise BootstrapFailure(
            "BOOTSTRAP_SOURCE_INVALID", "Operator bootstrap source is invalid."
        ) from None
    try:
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or not 0 < before.st_size <= MAX_TEMPLATE_BYTES
            or stat.S_IMODE(before.st_mode) & 0o022
        ):
            raise BootstrapFailure(
                "BOOTSTRAP_SOURCE_INVALID", "Operator bootstrap source is invalid."
            )
        content = bytearray()
        while len(content) <= MAX_TEMPLATE_BYTES:
            chunk = os.read(
                descriptor, min(65536, MAX_TEMPLATE_BYTES + 1 - len(content))
            )
            if not chunk:
                break
            content.extend(chunk)
        after = os.fstat(descriptor)
        if (
            len(content) != before.st_size
            or len(content) > MAX_TEMPLATE_BYTES
            or _file_witness(before) != _file_witness(after)
        ):
            raise BootstrapFailure(
                "BOOTSTRAP_SOURCE_INVALID", "Operator bootstrap source is invalid."
            )
        return bytes(content)
    finally:
        os.close(descriptor)


def _file_witness(value: os.stat_result) -> tuple[int, ...]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_nlink,
        value.st_uid,
        value.st_gid,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def _write_all(descriptor: int, content: bytes) -> None:
    offset = 0
    while offset < len(content):
        written = os.write(descriptor, content[offset:])
        if written <= 0:
            raise OSError
        offset += written


def _fsync_directory(path: Path) -> None:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_DIRECTORY", 0)
    descriptor = os.open(path, flags)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _run_fixed_command(command: Sequence[str]) -> bool:
    allowed = {
        (str(GROUPADD), "--system", OPERATOR_GROUP),
        *((str(GROUPADD), "--system", spec.group) for spec in SYSTEM_ACCOUNTS),
        *(_useradd_command(spec) for spec in SYSTEM_ACCOUNTS),
        (str(SYSTEMD_TMPFILES), "--create", "/usr/lib/tmpfiles.d/helixweave.conf"),
        (str(SYSTEMCTL), "daemon-reload"),
    }
    values = tuple(command)
    usermod = (
        len(values) == 6
        and values[:4]
        == (
            str(USERMOD),
            "--append",
            "--groups",
            OPERATOR_GROUP,
        )
        and values[4] == "--"
    )
    if values not in allowed and not usermod:
        return False
    try:
        completed = subprocess.run(
            values,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            env=SAFE_ENVIRONMENT,
            close_fds=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return False
    return completed.returncode == 0


def _useradd_command(spec: SystemAccountSpec) -> tuple[str, ...]:
    return (
        str(USERADD),
        "--system",
        "--gid",
        spec.group,
        "--home-dir",
        spec.home,
        "--no-create-home",
        "--shell",
        "/usr/sbin/nologin",
        "--",
        spec.user,
    )


def _validate_sudoers(candidate: Path) -> bool:
    try:
        completed = subprocess.run(
            (str(VISUDO), "-c", "-f", str(candidate)),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            env=SAFE_ENVIRONMENT,
            close_fds=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return False
    return completed.returncode == 0


def _resolve_invoking_user(environ: Mapping[str, str]) -> str:
    username = environ.get("SUDO_USER", "")
    raw_uid = environ.get("SUDO_UID", "")
    if (
        not _USERNAME.fullmatch(username)
        or not raw_uid.isascii()
        or not raw_uid.isdigit()
    ):
        raise BootstrapFailure(
            "BOOTSTRAP_CALLER_INVALID", "Interactive sudo caller is invalid."
        )
    try:
        uid = int(raw_uid)
        account = pwd.getpwnam(username)
    except (KeyError, ValueError):
        raise BootstrapFailure(
            "BOOTSTRAP_CALLER_INVALID", "Interactive sudo caller is invalid."
        ) from None
    if uid <= 0 or account.pw_uid != uid or account.pw_name != username:
        raise BootstrapFailure(
            "BOOTSTRAP_CALLER_INVALID", "Interactive sudo caller is invalid."
        )
    return username


def run_bootstrap(
    argv: Sequence[str],
    *,
    backend: BootstrapBackend,
    effective_uid: int,
    stdin_is_tty: bool,
    stderr_is_tty: bool,
    environ: Mapping[str, str],
) -> dict[str, object]:
    if tuple(argv) not in {("install",), ("update",)}:
        raise BootstrapFailure(
            "BOOTSTRAP_REQUEST_INVALID", "Operator bootstrap request is invalid."
        )
    if effective_uid != 0:
        raise BootstrapFailure(
            "BOOTSTRAP_PRIVILEGE_REQUIRED",
            "Operator bootstrap must run through interactive sudo.",
        )
    if not stdin_is_tty or not stderr_is_tty:
        raise BootstrapFailure(
            "BOOTSTRAP_INTERACTIVE_REQUIRED",
            "Operator bootstrap requires a local interactive terminal.",
        )
    invoking_user = _resolve_invoking_user(environ)
    result = backend.apply(operation=argv[0], invoking_user=invoking_user)
    if (
        not isinstance(result, BootstrapResult)
        or not isinstance(result.installed_files, int)
        or isinstance(result.installed_files, bool)
        or result.installed_files < 1
        or _CONTENT_IDENTITY.fullmatch(result.closure_identity) is None
    ):
        raise BootstrapFailure(
            "BOOTSTRAP_RESULT_INVALID", "Operator bootstrap result is invalid."
        )
    return {
        "schema_version": BOOTSTRAP_RECEIPT_SCHEMA,
        "operation": argv[0],
        "status": "complete",
        "installed_files": result.installed_files,
        "closure_identity": result.closure_identity,
    }


def _write_receipt(stream, value: object) -> None:
    stream.write(json.dumps(value, sort_keys=True, separators=(",", ":")))
    stream.write("\n")
    stream.flush()


def main(argv: Sequence[str] | None = None) -> int:
    arguments = tuple(sys.argv[1:] if argv is None else argv)
    caller_environment = {
        "SUDO_USER": os.environ.get("SUDO_USER", ""),
        "SUDO_UID": os.environ.get("SUDO_UID", ""),
    }
    os.environ.clear()
    os.environ.update(SAFE_ENVIRONMENT)
    try:
        receipt = run_bootstrap(
            arguments,
            backend=HostBootstrapBackend(),
            effective_uid=os.geteuid(),
            stdin_is_tty=sys.stdin.isatty(),
            stderr_is_tty=sys.stderr.isatty(),
            environ=caller_environment,
        )
    except BootstrapFailure as error:
        _write_receipt(
            sys.stderr,
            {
                "schema_version": BOOTSTRAP_RECEIPT_SCHEMA,
                "status": "error",
                "issue": {"code": error.code, "message": error.message},
            },
        )
        if error.code == "BOOTSTRAP_REQUEST_INVALID":
            return 64
        if error.code in {
            "BOOTSTRAP_PRIVILEGE_REQUIRED",
            "BOOTSTRAP_INTERACTIVE_REQUIRED",
            "BOOTSTRAP_CALLER_INVALID",
        }:
            return 77
        return 70
    except Exception:
        _write_receipt(
            sys.stderr,
            {
                "schema_version": BOOTSTRAP_RECEIPT_SCHEMA,
                "status": "error",
                "issue": {
                    "code": "BOOTSTRAP_INSTALL_FAILED",
                    "message": "Operator boundary could not be installed.",
                },
            },
        )
        return 70
    _write_receipt(sys.stdout, receipt)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
