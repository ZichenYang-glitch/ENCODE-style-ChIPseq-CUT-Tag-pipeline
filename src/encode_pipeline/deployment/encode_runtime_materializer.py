"""Unprivileged, offline materialization of one admitted ENCODE runtime.

The root operator creates the exact content-addressed destination and then
starts this code as the dedicated service account.  This module never chooses
a source or destination from caller-controlled paths, never consults the
process environment, and deliberately leaves an interrupted destination in
place for operator recovery.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shlex
import stat
import subprocess
import time
from typing import Any, Protocol

from encode_pipeline.deployment.bundle import BundleStore
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.filesystem import (
    fsync_directory,
    read_regular_file,
)
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import ENCODE_RUNTIME
from encode_pipeline.deployment.native_contracts import (
    ENCODE_MICROMAMBA_PATH,
    ENCODE_RUNTIME_INDEX_PATH,
    ENCODE_RUNTIME_ROOT_PATH,
    EncodeRuntimePackageIndex,
    ProductionNativeContractResolver,
    parse_encode_runtime_index,
)
from encode_pipeline.deployment.operator_action import (
    EncodeRuntimeEntry,
    EncodeRuntimeInventory,
    EncodeRuntimePrepareReceipt,
    EncodeRuntimePrepareRequest,
    snakemake_environment_hash,
)


RUNTIME_INVENTORY_FILENAME = ".helixweave-runtime-inventory.json"
CONDA_COMPAT_RELATIVE_PATH = "runner/bin/conda"
MICROMAMBA_RELATIVE_PATH = "runner/libexec/micromamba"
SNAKEMAKE_ACTIVATE_RELATIVE_PATH = "mamba-root/bin/activate"
MATERIALIZATION_TIMEOUT_SECONDS = 1_700.0

_MAX_INDEX_BYTES = 32 * 1024 * 1024
_MAX_SOURCE_FILE_BYTES = 8 * 1024 * 1024 * 1024
_MAX_INVENTORY_BYTES = 64 * 1024 * 1024
_MAX_RUNTIME_FILES = 500_000
_MAX_RUNTIME_FILE_BYTES = 8 * 1024 * 1024 * 1024
_MAX_RUNTIME_TOTAL_BYTES = 128 * 1024 * 1024 * 1024
_MAX_RELATIVE_PATH_BYTES = 1024
_COPY_CHUNK_BYTES = 1024 * 1024


class MicromambaCommandRunner(Protocol):
    """The exact ``subprocess.run`` surface used by the materializer."""

    def __call__(
        self,
        argv: Sequence[str],
        *,
        stdin: int,
        stdout: int,
        stderr: int,
        cwd: Path,
        env: Mapping[str, str],
        close_fds: bool,
        timeout: float,
        check: bool,
        shell: bool,
        umask: int,
    ) -> Any: ...


@dataclass(frozen=True)
class _EnvironmentPlan:
    prefix: Path
    explicit_lock: Path
    explicit_content: bytes
    marker: Path | None


class OfflineEncodeRuntimeMaterializer:
    """Materialize one exact admitted bundle without privilege or network."""

    def __init__(
        self,
        layout: DeploymentLayout,
        *,
        service_uid: int,
        service_gid: int,
        installed_owner_uid: int = 0,
        installed_owner_gid: int = 0,
        timeout_seconds: float = MATERIALIZATION_TIMEOUT_SECONDS,
        command_runner: MicromambaCommandRunner = subprocess.run,
        monotonic_clock: Callable[[], float] = time.monotonic,
    ) -> None:
        identifiers = (
            service_uid,
            service_gid,
            installed_owner_uid,
            installed_owner_gid,
        )
        if (
            not isinstance(layout, DeploymentLayout)
            or any(
                not isinstance(value, int) or isinstance(value, bool) or value < 0
                for value in identifiers
            )
            or not isinstance(timeout_seconds, (int, float))
            or isinstance(timeout_seconds, bool)
            or not 1.0 <= float(timeout_seconds) <= MATERIALIZATION_TIMEOUT_SECONDS
            or not callable(command_runner)
            or not callable(monotonic_clock)
        ):
            raise fail(
                "ENCODE_RUNTIME_MATERIALIZER_INVALID",
                "ENCODE runtime materializer is invalid.",
                component=ENCODE_RUNTIME,
            )
        self._layout = layout
        self._service_uid = service_uid
        self._service_gid = service_gid
        self._installed_owner_uid = installed_owner_uid
        self._installed_owner_gid = installed_owner_gid
        self._timeout_seconds = float(timeout_seconds)
        self._command_runner = command_runner
        self._monotonic_clock = monotonic_clock

    def prepare(
        self, request: EncodeRuntimePrepareRequest
    ) -> EncodeRuntimePrepareReceipt:
        """Build the fixed identity directory and return bounded evidence."""

        if not isinstance(request, EncodeRuntimePrepareRequest):
            raise fail(
                "ENCODE_RUNTIME_PREPARE_REQUEST_INVALID",
                "ENCODE runtime preparation request is invalid.",
                component=ENCODE_RUNTIME,
            )
        destination = self._layout.encode_runtime_active_root(
            request.deployment_identity
        )
        self._require_empty_destination(destination)
        source_root, index = self._admit_source(request.deployment_identity)
        plans = self._build_plans(source_root, destination, index)
        micromamba = source_root.joinpath(*PurePosixPath(ENCODE_MICROMAMBA_PATH).parts)
        deadline = self._monotonic_clock() + self._timeout_seconds

        try:
            self._make_directory(destination / "mamba-root")
            self._make_directory(destination / "mamba-root" / "explicit-locks")
            self._make_directory(destination / "conda-envs")
            for plan in plans:
                self._write_new_file(plan.explicit_lock, plan.explicit_content, 0o600)
                self._run_create(micromamba, destination, plan, deadline=deadline)
                self._require_created_prefix(plan.prefix)
                if plan.marker is not None:
                    self._write_new_file(plan.marker, b"", 0o600)
            self._require_before_deadline(deadline)
            self._install_runtime_execution_seam(
                micromamba, destination, deadline=deadline
            )
            entries = self._admit_and_normalize_tree(destination, deadline=deadline)
            inventory = EncodeRuntimeInventory.create(entries)
            inventory_content = inventory.to_bytes()
            if not 0 < len(inventory_content) <= _MAX_INVENTORY_BYTES:
                raise fail(
                    "ENCODE_RUNTIME_OUTPUT_INVALID",
                    "ENCODE runtime output is invalid.",
                    component=ENCODE_RUNTIME,
                    recoverable=True,
                )
            self._write_new_file(
                destination / RUNTIME_INVENTORY_FILENAME,
                inventory_content,
                0o600,
            )
            self._require_before_deadline(deadline)
            fsync_directory(destination)
            return EncodeRuntimePrepareReceipt.create(
                request_identity=request.identity,
                deployment_identity=request.deployment_identity,
                inventory=inventory,
            )
        except DeploymentError:
            raise
        except (OSError, subprocess.SubprocessError, ValueError):
            raise fail(
                "ENCODE_RUNTIME_MATERIALIZATION_FAILED",
                "ENCODE runtime materialization failed.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            ) from None

    def _require_empty_destination(self, destination: Path) -> None:
        try:
            observed = destination.lstat()
            with os.scandir(destination) as children:
                empty = next(children, None) is None
            after = destination.lstat()
        except OSError:
            raise fail(
                "ENCODE_RUNTIME_DESTINATION_INVALID",
                "ENCODE runtime destination is invalid.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            ) from None
        if (
            not stat.S_ISDIR(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or observed.st_uid != self._service_uid
            or observed.st_gid != self._service_gid
            or stat.S_IMODE(observed.st_mode) != 0o700
            or not empty
            or _stat_witness(observed) != _stat_witness(after)
        ):
            raise fail(
                "ENCODE_RUNTIME_DESTINATION_INVALID",
                "ENCODE runtime destination is invalid.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            )

    def _admit_source(self, identity: str) -> tuple[Path, EncodeRuntimePackageIndex]:
        source_root = self._layout.component_store(ENCODE_RUNTIME) / identity
        try:
            manifest = BundleStore(self._layout).verify_installed(
                ENCODE_RUNTIME,
                identity,
                expected_owner_uid=self._installed_owner_uid,
                expected_owner_gid=self._installed_owner_gid,
            )
            facts = ProductionNativeContractResolver().resolve(source_root, manifest)
            if (
                facts.component != ENCODE_RUNTIME
                or facts.deployment_identity != identity
                or facts.contracts != manifest.provides
            ):
                raise ValueError
            content, _observed = read_regular_file(
                source_root.joinpath(*PurePosixPath(ENCODE_RUNTIME_INDEX_PATH).parts),
                max_bytes=_MAX_INDEX_BYTES,
                code="ENCODE_RUNTIME_SOURCE_INVALID",
            )
            index = parse_encode_runtime_index(content)
        except Exception:
            raise fail(
                "ENCODE_RUNTIME_SOURCE_INVALID",
                "ENCODE runtime source is invalid.",
                component=ENCODE_RUNTIME,
            ) from None
        return source_root, index

    def _build_plans(
        self,
        source_root: Path,
        destination: Path,
        index: EncodeRuntimePackageIndex,
    ) -> tuple[_EnvironmentPlan, ...]:
        project_root = source_root.joinpath(
            *PurePosixPath(ENCODE_RUNTIME_ROOT_PATH).parts
        )
        conda_prefix = destination / "conda-envs"
        lock_root = destination / "mamba-root" / "explicit-locks"
        packages = {item["url"]: item for item in index.packages}
        runner: _EnvironmentPlan | None = None
        rules: dict[str, _EnvironmentPlan] = {}
        try:
            for environment in index.environments:
                environment_path = str(environment["environment_path"])
                environment_content, _ = read_regular_file(
                    project_root.joinpath(*PurePosixPath(environment_path).parts),
                    max_bytes=_MAX_SOURCE_FILE_BYTES,
                    code="ENCODE_RUNTIME_SOURCE_INVALID",
                )
                explicit_content = _local_explicit_lock(
                    source_root, tuple(environment["packages"]), packages
                )
                if environment_path == "workflow/envs/runner.yml":
                    runner = _EnvironmentPlan(
                        prefix=destination / "runner",
                        explicit_lock=lock_root / "runner.lock",
                        explicit_content=explicit_content,
                        marker=None,
                    )
                    continue
                full_hash = snakemake_environment_hash(
                    conda_prefix, environment_content
                )
                plan = _EnvironmentPlan(
                    prefix=conda_prefix / full_hash,
                    explicit_lock=lock_root / f"{full_hash}.lock",
                    explicit_content=explicit_content,
                    marker=conda_prefix / f"{full_hash}.env_setup_done",
                )
                prior = rules.get(full_hash)
                if prior is not None and prior.explicit_content != explicit_content:
                    raise ValueError
                rules[full_hash] = plan
        except DeploymentError:
            raise
        except (KeyError, OSError, TypeError, ValueError):
            raise fail(
                "ENCODE_RUNTIME_SOURCE_INVALID",
                "ENCODE runtime source is invalid.",
                component=ENCODE_RUNTIME,
            ) from None
        if runner is None:
            raise fail(
                "ENCODE_RUNTIME_SOURCE_INVALID",
                "ENCODE runtime source is invalid.",
                component=ENCODE_RUNTIME,
            )
        return (runner, *(rules[key] for key in sorted(rules)))

    def _run_create(
        self,
        micromamba: Path,
        destination: Path,
        plan: _EnvironmentPlan,
        *,
        deadline: float,
    ) -> None:
        argv = (
            str(micromamba),
            "--no-rc",
            "--root-prefix",
            str(destination / "mamba-root"),
            "create",
            "--yes",
            "--offline",
            "--always-copy",
            "--prefix",
            str(plan.prefix),
            "--file",
            str(plan.explicit_lock),
        )
        try:
            remaining = deadline - self._monotonic_clock()
            if not 0 < remaining <= self._timeout_seconds:
                raise subprocess.TimeoutExpired(argv, timeout=self._timeout_seconds)
            completed = self._command_runner(
                argv,
                stdin=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                cwd=destination,
                env={},
                close_fds=True,
                timeout=remaining,
                check=False,
                shell=False,
                umask=0o077,
            )
        except (OSError, subprocess.SubprocessError):
            raise fail(
                "ENCODE_RUNTIME_MATERIALIZATION_FAILED",
                "ENCODE runtime materialization failed.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            ) from None
        return_code = getattr(completed, "returncode", None)
        if (
            not isinstance(return_code, int)
            or isinstance(return_code, bool)
            or return_code != 0
        ):
            raise fail(
                "ENCODE_RUNTIME_MATERIALIZATION_FAILED",
                "ENCODE runtime materialization failed.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            )

    def _require_created_prefix(self, prefix: Path) -> None:
        try:
            observed = prefix.lstat()
        except OSError:
            raise fail(
                "ENCODE_RUNTIME_OUTPUT_INVALID",
                "ENCODE runtime output is invalid.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            ) from None
        if (
            not stat.S_ISDIR(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or observed.st_uid != self._service_uid
            or observed.st_gid != self._service_gid
        ):
            raise fail(
                "ENCODE_RUNTIME_OUTPUT_INVALID",
                "ENCODE runtime output is invalid.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            )

    def _install_runtime_execution_seam(
        self, source: Path, root: Path, *, deadline: float
    ) -> None:
        """Install the immutable Snakemake 8.30 conda compatibility seam.

        Environment creation already happened through the indexed static
        micromamba above.  Runtime Snakemake must therefore be able to perform
        only its four fixed, read-only conda probes and activate an already
        materialized full-hash prefix.  Every other conda invocation fails.
        """

        self._make_directory(root / "runner" / "libexec")
        self._copy_indexed_micromamba(
            source,
            root.joinpath(*PurePosixPath(MICROMAMBA_RELATIVE_PATH).parts),
            deadline=deadline,
        )
        self._write_new_file(
            root.joinpath(*PurePosixPath(CONDA_COMPAT_RELATIVE_PATH).parts),
            _conda_compat_script(root),
            0o700,
        )
        self._make_directory(root / "mamba-root" / "bin")
        self._write_new_file(
            root.joinpath(*PurePosixPath(SNAKEMAKE_ACTIVATE_RELATIVE_PATH).parts),
            _snakemake_activate_script(root),
            0o700,
        )
        self._require_before_deadline(deadline)

    def _copy_indexed_micromamba(
        self, source: Path, destination: Path, *, deadline: float
    ) -> None:
        try:
            parent = destination.parent.lstat()
        except OSError:
            raise fail(
                "ENCODE_RUNTIME_OUTPUT_INVALID",
                "ENCODE runtime output is invalid.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            ) from None
        if not stat.S_ISDIR(parent.st_mode) or stat.S_ISLNK(parent.st_mode):
            raise fail(
                "ENCODE_RUNTIME_OUTPUT_INVALID",
                "ENCODE runtime output is invalid.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            )
        self._copy_new_regular(source, destination, mode=0o700, deadline=deadline)

    def _admit_and_normalize_tree(
        self, root: Path, *, deadline: float
    ) -> tuple[EncodeRuntimeEntry, ...]:
        entries: list[EncodeRuntimeEntry] = []
        total_size = 0
        node_count = 0
        pending = [root]
        while pending:
            current = pending.pop()
            try:
                with os.scandir(current) as iterator:
                    children = sorted(iterator, key=lambda item: item.name)
            except OSError:
                raise fail(
                    "ENCODE_RUNTIME_OUTPUT_INVALID",
                    "ENCODE runtime output is invalid.",
                    component=ENCODE_RUNTIME,
                    recoverable=True,
                ) from None
            for child in children:
                self._require_before_deadline(deadline)
                path = Path(child.path)
                relative = path.relative_to(root).as_posix()
                if relative == RUNTIME_INVENTORY_FILENAME:
                    raise fail(
                        "ENCODE_RUNTIME_OUTPUT_INVALID",
                        "ENCODE runtime output is invalid.",
                        component=ENCODE_RUNTIME,
                        recoverable=True,
                    )
                node_count += 1
                if (
                    node_count > _MAX_RUNTIME_FILES
                    or len(relative.encode("utf-8")) > _MAX_RELATIVE_PATH_BYTES
                ):
                    raise fail(
                        "ENCODE_RUNTIME_OUTPUT_INVALID",
                        "ENCODE runtime output is invalid.",
                        component=ENCODE_RUNTIME,
                        recoverable=True,
                    )
                observed = path.lstat()
                if (
                    observed.st_uid != self._service_uid
                    or observed.st_gid != self._service_gid
                ):
                    raise fail(
                        "ENCODE_RUNTIME_OUTPUT_INVALID",
                        "ENCODE runtime output is invalid.",
                        component=ENCODE_RUNTIME,
                        recoverable=True,
                    )
                if stat.S_ISDIR(observed.st_mode) and not stat.S_ISLNK(
                    observed.st_mode
                ):
                    os.chmod(path, 0o700, follow_symlinks=False)
                    pending.append(path)
                    continue
                if stat.S_ISREG(observed.st_mode):
                    if not 0 <= observed.st_size <= _MAX_RUNTIME_FILE_BYTES:
                        raise fail(
                            "ENCODE_RUNTIME_OUTPUT_INVALID",
                            "ENCODE runtime output is invalid.",
                            component=ENCODE_RUNTIME,
                            recoverable=True,
                        )
                    executable = bool(stat.S_IMODE(observed.st_mode) & 0o111)
                    if observed.st_nlink != 1:
                        self._break_hardlink(
                            path,
                            relative,
                            mode=0o700 if executable else 0o600,
                            deadline=deadline,
                        )
                        observed = path.lstat()
                    os.chmod(path, 0o700 if executable else 0o600)
                    digest, size = self._hash_file(path, deadline=deadline)
                    total_size += size
                    entries.append(
                        EncodeRuntimeEntry.from_dict(
                            {
                                "path": relative,
                                "kind": "file",
                                "sha256": digest,
                                "size": size,
                                "mode": 0o555 if executable else 0o444,
                                "target": None,
                            }
                        )
                    )
                    continue
                if stat.S_ISLNK(observed.st_mode):
                    try:
                        target = os.readlink(path)
                    except OSError:
                        raise fail(
                            "ENCODE_RUNTIME_OUTPUT_INVALID",
                            "ENCODE runtime output is invalid.",
                            component=ENCODE_RUNTIME,
                            recoverable=True,
                        ) from None
                    _require_safe_symlink(root, relative, path, target)
                    target_bytes = target.encode("utf-8")
                    total_size += len(target_bytes)
                    entries.append(
                        EncodeRuntimeEntry.from_dict(
                            {
                                "path": relative,
                                "kind": "symlink",
                                "sha256": hashlib.sha256(target_bytes).hexdigest(),
                                "size": len(target_bytes),
                                "mode": None,
                                "target": target,
                            }
                        )
                    )
                    continue
                raise fail(
                    "ENCODE_RUNTIME_OUTPUT_INVALID",
                    "ENCODE runtime output is invalid.",
                    component=ENCODE_RUNTIME,
                    recoverable=True,
                )
        if total_size > _MAX_RUNTIME_TOTAL_BYTES:
            raise fail(
                "ENCODE_RUNTIME_OUTPUT_INVALID",
                "ENCODE runtime output is invalid.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            )
        ordered = tuple(sorted(entries, key=lambda item: item.path))
        if not ordered:
            raise fail(
                "ENCODE_RUNTIME_OUTPUT_INVALID",
                "ENCODE runtime output is invalid.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            )
        return ordered

    def _break_hardlink(
        self, path: Path, relative: str, *, mode: int, deadline: float
    ) -> None:
        temporary = path.parent / (
            ".helixweave-hardlink-"
            + hashlib.sha256(relative.encode("utf-8")).hexdigest()
        )
        self._copy_new_regular(path, temporary, mode=mode, deadline=deadline)
        try:
            os.replace(temporary, path)
            fsync_directory(path.parent)
            observed = path.lstat()
        except OSError:
            raise fail(
                "ENCODE_RUNTIME_OUTPUT_INVALID",
                "ENCODE runtime output is invalid.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            ) from None
        if (
            observed.st_nlink != 1
            or observed.st_uid != self._service_uid
            or observed.st_gid != self._service_gid
        ):
            raise fail(
                "ENCODE_RUNTIME_OUTPUT_INVALID",
                "ENCODE runtime output is invalid.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            )

    def _make_directory(self, path: Path) -> None:
        try:
            path.mkdir(mode=0o700)
            observed = path.lstat()
        except OSError:
            raise fail(
                "ENCODE_RUNTIME_MATERIALIZATION_FAILED",
                "ENCODE runtime materialization failed.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            ) from None
        if (
            not stat.S_ISDIR(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or observed.st_uid != self._service_uid
            or observed.st_gid != self._service_gid
            or stat.S_IMODE(observed.st_mode) != 0o700
        ):
            raise fail(
                "ENCODE_RUNTIME_MATERIALIZATION_FAILED",
                "ENCODE runtime materialization failed.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            )

    @staticmethod
    def _write_new_file(path: Path, content: bytes, mode: int) -> None:
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
            _write_all(descriptor, content)
            os.fchmod(descriptor, mode)
            os.fsync(descriptor)
        except OSError:
            raise fail(
                "ENCODE_RUNTIME_MATERIALIZATION_FAILED",
                "ENCODE runtime materialization failed.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            ) from None
        finally:
            if descriptor >= 0:
                os.close(descriptor)
        fsync_directory(path.parent)

    def _copy_new_regular(
        self, source: Path, destination: Path, *, mode: int, deadline: float
    ) -> None:
        source_descriptor = -1
        destination_descriptor = -1
        try:
            source_descriptor = os.open(
                source,
                os.O_RDONLY
                | getattr(os, "O_CLOEXEC", 0)
                | getattr(os, "O_NOFOLLOW", 0),
            )
            before = os.fstat(source_descriptor)
            if (
                not stat.S_ISREG(before.st_mode)
                or before.st_nlink < 1
                or not 0 <= before.st_size <= _MAX_RUNTIME_FILE_BYTES
            ):
                raise OSError
            destination_descriptor = os.open(
                destination,
                os.O_WRONLY
                | os.O_CREAT
                | os.O_EXCL
                | getattr(os, "O_CLOEXEC", 0)
                | getattr(os, "O_NOFOLLOW", 0),
                0o600,
            )
            remaining = before.st_size
            while remaining:
                self._require_before_deadline(deadline)
                chunk = os.read(source_descriptor, min(_COPY_CHUNK_BYTES, remaining))
                if not chunk:
                    raise OSError
                _write_all(destination_descriptor, chunk)
                remaining -= len(chunk)
            if os.read(source_descriptor, 1):
                raise OSError
            after = os.fstat(source_descriptor)
            if _stat_witness(before) != _stat_witness(after):
                raise OSError
            os.fchmod(destination_descriptor, mode)
            os.fsync(destination_descriptor)
        except OSError:
            raise fail(
                "ENCODE_RUNTIME_OUTPUT_INVALID",
                "ENCODE runtime output is invalid.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            ) from None
        finally:
            if source_descriptor >= 0:
                os.close(source_descriptor)
            if destination_descriptor >= 0:
                os.close(destination_descriptor)
        fsync_directory(destination.parent)

    def _hash_file(self, path: Path, *, deadline: float) -> tuple[str, int]:
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
                or not 0 <= before.st_size <= _MAX_RUNTIME_FILE_BYTES
            ):
                raise OSError
            digest = hashlib.sha256()
            remaining = before.st_size
            while remaining:
                self._require_before_deadline(deadline)
                chunk = os.read(descriptor, min(_COPY_CHUNK_BYTES, remaining))
                if not chunk:
                    raise OSError
                digest.update(chunk)
                remaining -= len(chunk)
            if os.read(descriptor, 1):
                raise OSError
            after = os.fstat(descriptor)
            path_after = path.stat(follow_symlinks=False)
            if _stat_witness(before) != _stat_witness(after) or (
                after.st_dev,
                after.st_ino,
            ) != (path_after.st_dev, path_after.st_ino):
                raise OSError
            return digest.hexdigest(), before.st_size
        except OSError:
            raise fail(
                "ENCODE_RUNTIME_OUTPUT_INVALID",
                "ENCODE runtime output is invalid.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            ) from None
        finally:
            if descriptor >= 0:
                os.close(descriptor)

    def _require_before_deadline(self, deadline: float) -> None:
        if self._monotonic_clock() >= deadline:
            raise fail(
                "ENCODE_RUNTIME_MATERIALIZATION_FAILED",
                "ENCODE runtime materialization failed.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            )


def _local_explicit_lock(
    source_root: Path,
    urls: tuple[object, ...],
    packages: Mapping[object, Mapping[str, object]],
) -> bytes:
    lines = ["# platform: linux-64", "@EXPLICIT"]
    for url in urls:
        package = packages[url]
        archive_path = package["archive_path"]
        md5 = package["md5"]
        if not isinstance(archive_path, str) or not isinstance(md5, str):
            raise ValueError
        archive = source_root.joinpath(*PurePosixPath(archive_path).parts)
        lines.append(f"{archive.as_uri()}#{md5}")
    return ("\n".join(lines) + "\n").encode("utf-8")


def _conda_compat_script(root: Path) -> bytes:
    """Return the closed read-only command surface required by Snakemake 8.30."""

    root_prefix = root / "mamba-root"
    _require_script_coordinate(root_prefix)
    info = json.dumps(
        {"conda_prefix": str(root_prefix), "platform": "linux-64"},
        sort_keys=True,
        separators=(",", ":"),
    )
    config = json.dumps(
        {"get": {"channel_priority": "strict"}},
        sort_keys=True,
        separators=(",", ":"),
    )
    return (
        "#!/bin/sh\n"
        'if [ "$#" -eq 2 ] && [ "$1" = info ] && [ "$2" = --json ]; then\n'
        f"    printf '%s\\n' {shlex.quote(info)}\n"
        "    exit 0\n"
        "fi\n"
        'if [ "$#" -eq 1 ] && [ "$1" = --version ]; then\n'
        "    printf '%s\\n' 'conda 24.7.1'\n"
        "    exit 0\n"
        "fi\n"
        'if [ "$#" -eq 4 ] && [ "$1" = config ] '
        '&& [ "$2" = --get ] && [ "$3" = channel_priority ] '
        '&& [ "$4" = --json ]; then\n'
        f"    printf '%s\\n' {shlex.quote(config)}\n"
        "    exit 0\n"
        "fi\n"
        "exit 64\n"
    ).encode("utf-8")


def _snakemake_activate_script(root: Path) -> bytes:
    """Return a source-only activator for one pre-materialized rule prefix."""

    environment_root = root / "conda-envs"
    runner_bin = root / "runner" / "bin"
    _require_script_coordinate(environment_root)
    _require_script_coordinate(runner_bin)
    return (
        "#!/bin/sh\n"
        'if [ "$#" -ne 1 ]; then\n'
        "    return 64\n"
        "fi\n"
        f"_helixweave_environment_root={shlex.quote(str(environment_root))}\n"
        f"_helixweave_runner_bin={shlex.quote(str(runner_bin))}\n"
        "_helixweave_prefix=${1-}\n"
        'case "$_helixweave_prefix" in\n'
        '    "$_helixweave_environment_root"/*) ;;\n'
        "    *) return 64 ;;\n"
        "esac\n"
        '_helixweave_name=${_helixweave_prefix#"$_helixweave_environment_root"/}\n'
        'if [ "${#_helixweave_name}" -ne 32 ]; then\n'
        "    return 64\n"
        "fi\n"
        'case "$_helixweave_name" in\n'
        "    *[!0123456789abcdef]*) return 64 ;;\n"
        "esac\n"
        'if [ ! -d "$_helixweave_prefix/bin" ]; then\n'
        "    return 69\n"
        "fi\n"
        'CONDA_PREFIX="$_helixweave_prefix"\n'
        'CONDA_DEFAULT_ENV="$_helixweave_prefix"\n'
        "CONDA_SHLVL=1\n"
        'CONDA_EXE="$_helixweave_runner_bin/conda"\n'
        '_CONDA_EXE="$CONDA_EXE"\n'
        f"_CONDA_ROOT={shlex.quote(str(root / 'mamba-root'))}\n"
        'PATH="$_helixweave_prefix/bin:$_helixweave_runner_bin:'
        '/usr/sbin:/usr/bin:/sbin:/bin"\n'
        "export CONDA_PREFIX CONDA_DEFAULT_ENV CONDA_SHLVL CONDA_EXE "
        "_CONDA_EXE _CONDA_ROOT PATH\n"
        "unset _helixweave_environment_root _helixweave_runner_bin "
        "_helixweave_prefix _helixweave_name\n"
    ).encode("utf-8")


def _require_script_coordinate(path: Path) -> None:
    text = str(path)
    if not path.is_absolute() or any(
        character in text for character in ("\x00", "\n", "\r")
    ):
        raise fail(
            "ENCODE_RUNTIME_MATERIALIZATION_FAILED",
            "ENCODE runtime materialization failed.",
            component=ENCODE_RUNTIME,
            recoverable=True,
        )


def _require_safe_symlink(root: Path, relative: str, path: Path, target: str) -> None:
    if (
        not target
        or target.startswith("/")
        or len(target.encode("utf-8")) > _MAX_RELATIVE_PATH_BYTES
        or any(character in target for character in ("\x00", "\n", "\r"))
    ):
        raise fail(
            "ENCODE_RUNTIME_OUTPUT_INVALID",
            "ENCODE runtime output is invalid.",
            component=ENCODE_RUNTIME,
            recoverable=True,
        )
    resolved_parts = list(PurePosixPath(relative).parts[:-1])
    for part in PurePosixPath(target).parts:
        if part in {"", "."}:
            continue
        if part == "..":
            if not resolved_parts:
                raise fail(
                    "ENCODE_RUNTIME_OUTPUT_INVALID",
                    "ENCODE runtime output is invalid.",
                    component=ENCODE_RUNTIME,
                    recoverable=True,
                )
            resolved_parts.pop()
            continue
        resolved_parts.append(part)
    if not resolved_parts:
        raise fail(
            "ENCODE_RUNTIME_OUTPUT_INVALID",
            "ENCODE runtime output is invalid.",
            component=ENCODE_RUNTIME,
            recoverable=True,
        )
    try:
        resolved = path.resolve(strict=True)
        resolved.relative_to(root.resolve(strict=True))
    except (OSError, RuntimeError, ValueError):
        raise fail(
            "ENCODE_RUNTIME_OUTPUT_INVALID",
            "ENCODE runtime output is invalid.",
            component=ENCODE_RUNTIME,
            recoverable=True,
        ) from None


def _write_all(descriptor: int, content: bytes) -> None:
    offset = 0
    while offset < len(content):
        written = os.write(descriptor, content[offset:])
        if written <= 0:
            raise OSError
        offset += written


def _stat_witness(value: os.stat_result) -> tuple[int, ...]:
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
    "MATERIALIZATION_TIMEOUT_SECONDS",
    "RUNTIME_INVENTORY_FILENAME",
    "MicromambaCommandRunner",
    "OfflineEncodeRuntimeMaterializer",
    "snakemake_environment_hash",
]
