"""Offline platform-Python closure and fixed candidate service entry points.

The supported platform release uses the host's fixed CPython 3.12 executable,
but never its ``site-packages``.  Release construction expands an exact local,
hash-locked wheel set into an immutable application directory.  The target host
does not run pip and the root operator never imports candidate code.

This module deliberately keeps dependency versions in wheel ``METADATA`` and
file identities in the deployment bundle manifest.  The wheel lock binds only
filenames to source-byte hashes; it is not a second distribution catalogue.
"""

from __future__ import annotations

from collections import deque
from collections.abc import Callable, Sequence
import base64
import csv
from dataclasses import dataclass
from email.parser import BytesParser
from email.policy import compat32
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import platform
import pwd
import grp
import re
import shutil
import socket
import stat
import sys
import uuid
import zipfile

from encode_pipeline.deployment.canonical import (
    canonical_identity,
    canonical_json_bytes,
)
from encode_pipeline.deployment.bundle import (
    MAX_BUNDLE_BYTES,
    MAX_BUNDLE_FILES,
    MAX_MANIFEST_BYTES,
    MIN_FREE_SPACE_RESERVE,
)
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.operator_action import (
    BulkRuntimePrepareReceipt,
    BulkRuntimePrepareRequest,
    DatabasePrepareReceipt,
    DatabasePrepareRequest,
    DeploymentActionReceipt,
    DeploymentActionRequest,
    EncodeRuntimePrepareReceipt,
    EncodeRuntimePrepareRequest,
    ReadinessCheck,
    VERIFICATION_CHECKS,
)


SUPPORTED_PYTHON_EXECUTABLE = Path("/usr/bin/python3.12")
SUPPORTED_PYTHON_ABI = "cp312"
SUPPORTED_PYTHON_CACHE_TAG = "cpython-312"
SUPPORTED_PYTHON_IMPLEMENTATION_NAME = "cpython"
SUPPORTED_PYTHON_VERSION = (3, 12)
SUPPORTED_PYTHON_MINIMUM_VERSION = (3, 12, 3)
SUPPORTED_PYTHON_MARKER_FLOOR = "3.12.3"
SUPPORTED_PYTHON_MARKER_CEILING = "3.12.999999"
SUPPORTED_MACHINE = "x86_64"
SUPPORTED_SYSTEM = "linux"
MINIMUM_GLIBC = (2, 17)

PLATFORM_WHEEL_LOCK_SCHEMA = "helixweave-platform-wheel-lock-v1"
PLATFORM_WHEEL_LOCK_IDENTITY_SCHEME = "helixweave-platform-wheel-lock-identity-v1"
PLATFORM_RUNTIME_OBSERVATION_IDENTITY_SCHEME = (
    "helixweave-platform-python-observation-identity-v1"
)
PLATFORM_RUNTIME_ROOT = "payload/platform"
PLATFORM_RUNTIME_LOCK_PATH = "payload/contracts/platform/python-runtime-wheel-lock.json"
PLATFORM_RUNTIME_WHEELHOUSE_ROOT = f"{PLATFORM_RUNTIME_ROOT}/wheelhouse"
PLATFORM_RUNTIME_SITE_PACKAGES = f"{PLATFORM_RUNTIME_ROOT}/lib/python3.12/site-packages"
PLATFORM_RUNTIME_LAUNCHER_PATH = f"{PLATFORM_RUNTIME_ROOT}/bin/helixweave-service"
PLATFORM_RUNTIME_LAUNCHER_MEMBER = (
    "encode_pipeline/deployment/templates/helixweave-service"
)

PLATFORM_ACTION_REQUEST_PATH = Path("/run/helixweave/operator/action/request.json")
DATABASE_PREPARE_REQUEST_PATH = Path("/run/helixweave/database/prepare.json")
ENCODE_RUNTIME_PREPARE_REQUEST_PATH = Path(
    "/run/helixweave/operator/encode-runtime/request.json"
)
BULK_RUNTIME_PREPARE_REQUEST_PATH = Path(
    "/run/helixweave/operator/bulk-runtime/request.json"
)
SUPPORTED_DATABASE_PATH = Path("/var/lib/helixweave/database/live/platform.db")
SUPPORTED_REDIS_SOCKET = Path("/run/helixweave/redis/redis.sock")
SUPPORTED_DOCKER_SOCKET = Path("/run/helixweave/docker/docker.sock")
SUPPORTED_REFERENCE_CONFIG = Path("/etc/helixweave/reference-profiles.yaml")

_IDENTITY = re.compile(r"^sha256-[0-9a-f]{64}$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_WHEEL_FILENAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.+!-]{0,239}\.whl$")
_DISTRIBUTION = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,191}$")
_VERSION = re.compile(r"^[0-9A-Za-z][0-9A-Za-z._+!-]{0,191}$")
_MANYLINUX = re.compile(r"^manylinux_(\d+)_(\d+)_x86_64$")
_LEGACY_MANYLINUX = frozenset(
    {"manylinux1_x86_64", "manylinux2010_x86_64", "manylinux2014_x86_64"}
)

_MAX_LOCK_BYTES = 2 * 1024 * 1024
_MAX_WHEELS = 256
_MAX_WHEEL_BYTES = 1024 * 1024 * 1024
_MAX_WHEEL_FILES = 100_000
_MAX_WHEEL_MEMBER_BYTES = 512 * 1024 * 1024
_MAX_EXPANDED_BYTES = 4 * 1024 * 1024 * 1024
_MAX_REQUEST_BYTES = 64 * 1024
_READ_CHUNK = 1024 * 1024
_PLATFORM_ROOT_EXTRAS = frozenset({"api", "llm-openai"})


def _invalid_closure() -> DeploymentError:
    return fail(
        "PLATFORM_RUNTIME_CLOSURE_INVALID",
        "Platform Python runtime closure is invalid.",
        component="platform",
    )


def _invalid_lock() -> DeploymentError:
    return fail(
        "PLATFORM_RUNTIME_LOCK_INVALID",
        "Platform Python runtime lock is invalid.",
        component="platform",
    )


@dataclass(frozen=True, order=True)
class LockedWheel:
    filename: str
    sha256: str

    def to_dict(self) -> dict[str, str]:
        return {"filename": self.filename, "sha256": self.sha256}


@dataclass(frozen=True)
class PlatformWheelLock:
    wheels: tuple[LockedWheel, ...]

    @classmethod
    def from_bytes(cls, content: bytes) -> "PlatformWheelLock":
        try:
            value = json.loads(
                content,
                object_pairs_hook=_unique_object,
                parse_constant=_reject_json_constant,
            )
        except (UnicodeDecodeError, ValueError, json.JSONDecodeError):
            raise _invalid_lock() from None
        if (
            not isinstance(value, dict)
            or set(value) != {"schema_version", "python", "wheels"}
            or value["schema_version"] != PLATFORM_WHEEL_LOCK_SCHEMA
            or value["python"]
            != {
                "abi": SUPPORTED_PYTHON_ABI,
                "executable": str(SUPPORTED_PYTHON_EXECUTABLE),
                "implementation_name": SUPPORTED_PYTHON_IMPLEMENTATION_NAME,
                "machine": SUPPORTED_MACHINE,
                "minimum_python_full_version": SUPPORTED_PYTHON_MARKER_FLOOR,
                "python_version": ".".join(
                    str(item) for item in SUPPORTED_PYTHON_VERSION
                ),
                "system": SUPPORTED_SYSTEM,
            }
            or not isinstance(value["wheels"], list)
            or not 0 < len(value["wheels"]) <= _MAX_WHEELS
        ):
            raise _invalid_lock()
        wheels: list[LockedWheel] = []
        for raw in value["wheels"]:
            if (
                not isinstance(raw, dict)
                or set(raw) != {"filename", "sha256"}
                or not isinstance(raw["filename"], str)
                or _WHEEL_FILENAME.fullmatch(raw["filename"]) is None
                or PurePosixPath(raw["filename"]).name != raw["filename"]
                or not isinstance(raw["sha256"], str)
                or _SHA256.fullmatch(raw["sha256"]) is None
            ):
                raise _invalid_lock()
            wheels.append(LockedWheel(raw["filename"], raw["sha256"]))
        ordered = tuple(sorted(wheels))
        if tuple(wheels) != ordered or len({item.filename for item in ordered}) != len(
            ordered
        ):
            raise _invalid_lock()
        result = cls(ordered)
        if content != result.to_bytes():
            raise _invalid_lock()
        return result

    @property
    def identity(self) -> str:
        return canonical_identity(
            self.to_dict(), scheme=PLATFORM_WHEEL_LOCK_IDENTITY_SCHEME
        )

    @property
    def python_version(self) -> str:
        return ".".join(str(item) for item in SUPPORTED_PYTHON_VERSION)

    @property
    def minimum_python_full_version(self) -> str:
        return SUPPORTED_PYTHON_MARKER_FLOOR

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": PLATFORM_WHEEL_LOCK_SCHEMA,
            "python": {
                "abi": SUPPORTED_PYTHON_ABI,
                "executable": str(SUPPORTED_PYTHON_EXECUTABLE),
                "implementation_name": SUPPORTED_PYTHON_IMPLEMENTATION_NAME,
                "machine": SUPPORTED_MACHINE,
                "minimum_python_full_version": SUPPORTED_PYTHON_MARKER_FLOOR,
                "python_version": self.python_version,
                "system": SUPPORTED_SYSTEM,
            },
            "wheels": [item.to_dict() for item in self.wheels],
        }

    def to_bytes(self) -> bytes:
        return canonical_json_bytes(self.to_dict())


@dataclass(frozen=True, order=True)
class PlatformRuntimeFile:
    """One immutable source file ready for ``bundle_builder._path_payload``."""

    logical_path: str
    source: Path
    mode: int
    size_bytes: int
    sha256: str


@dataclass(frozen=True)
class PlatformRuntimeClosure:
    lock_identity: str
    files: tuple[PlatformRuntimeFile, ...]

    def source_for(self, logical_path: str) -> Path:
        matches = tuple(
            item.source for item in self.files if item.logical_path == logical_path
        )
        if len(matches) != 1:
            raise _invalid_closure()
        return matches[0]


@dataclass(frozen=True)
class _WheelDistribution:
    name: str
    version: str
    requires_python: str | None
    requirements: tuple[str, ...]


@dataclass(frozen=True)
class _WheelPlan:
    locked: LockedWheel
    source: Path
    distribution: _WheelDistribution
    members: tuple[tuple[str, str], ...]


def build_platform_wheel_lock(wheelhouse: Path) -> PlatformWheelLock:
    """Derive the reviewed canonical lock for one exact local wheelhouse.

    This operation is deliberately offline and non-mutating.  It validates the
    same wheel tags, RECORD bytes, dependency reachability, and collision rules
    as collection before returning filename/hash bindings for review.
    """

    _verify_build_python()
    source_root = _source_directory(wheelhouse)
    inventory = _wheelhouse_inventory(source_root)
    if not 0 < len(inventory) <= _MAX_WHEELS:
        raise _invalid_lock()
    locked = tuple(
        LockedWheel(
            name,
            hashlib.sha256(_read_regular_source(path, _MAX_WHEEL_BYTES)).hexdigest(),
        )
        for name, path in sorted(inventory.items())
    )
    lock = PlatformWheelLock(locked)
    plans = tuple(_plan_wheel(item, inventory[item.filename]) for item in lock.wheels)
    _verify_dependency_closure(plans, lock)
    _expected_installed_files(plans)
    for item in lock.wheels:
        if _hash_file(inventory[item.filename]) != item.sha256:
            raise _invalid_closure()
    return lock


def collect_platform_runtime_closure(
    wheelhouse: Path,
    lock_path: Path,
    destination: Path,
) -> PlatformRuntimeClosure:
    """Snapshot and expand one exact local wheel-only platform closure.

    ``destination`` must not exist.  Its parent is the only filesystem location
    mutated.  The returned logical paths already use deployment-bundle
    coordinates, so the shared bundle builder needs no path policy of its own.
    """

    _verify_build_python()
    wheelhouse = _source_directory(wheelhouse)
    lock_content = _read_regular_source(lock_path, _MAX_LOCK_BYTES)
    lock = PlatformWheelLock.from_bytes(lock_content)
    destination = _new_destination(destination)
    partial = destination.parent / f".{destination.name}.{uuid.uuid4().hex}.partial"
    try:
        expected = {item.filename: item for item in lock.wheels}
        observed = _wheelhouse_inventory(wheelhouse)
        if set(observed) != set(expected):
            raise _invalid_lock()
        source_plans = tuple(
            _plan_wheel(item, observed[item.filename]) for item in lock.wheels
        )
        _verify_dependency_closure(source_plans, lock)
        expected_site, expected_launcher = _expected_installed_files(source_plans)
        wheel_bytes = 0
        for item in lock.wheels:
            source = observed[item.filename]
            if _hash_file(source) != item.sha256:
                raise _invalid_lock()
            wheel_bytes += source.lstat().st_size
        _preflight_closure_capacity(
            destination.parent,
            file_count=2 + len(lock.wheels) + len(expected_site),
            payload_bytes=(
                len(lock_content)
                + wheel_bytes
                + sum(size for size, _sha256 in expected_site.values())
                + expected_launcher[0]
            ),
        )

        partial.mkdir(mode=0o700)
        lock_destination = partial.joinpath(
            *PurePosixPath(PLATFORM_RUNTIME_LOCK_PATH).parts
        )
        _write_bytes(lock_destination, lock_content, 0o444)
        wheel_root = partial.joinpath(
            *PurePosixPath(PLATFORM_RUNTIME_WHEELHOUSE_ROOT).parts
        )
        wheel_root.mkdir(parents=True, mode=0o700)
        snapshots: list[tuple[LockedWheel, Path]] = []
        for item in lock.wheels:
            target = wheel_root / item.filename
            _copy_hashed_source(observed[item.filename], target, item.sha256)
            snapshots.append((item, target))

        plans = tuple(_plan_wheel(item, source) for item, source in snapshots)
        _verify_dependency_closure(plans, lock)
        _extract_wheels(partial, plans)
        _install_candidate_launcher(partial, plans)
        _freeze_tree(partial)
        os.replace(partial, destination)
        _fsync_directory(destination.parent)
        return inspect_platform_runtime_closure(destination, lock.identity)
    except DeploymentError:
        _discard_partial(partial)
        raise
    except Exception:
        _discard_partial(partial)
        raise _invalid_closure() from None


def inspect_platform_runtime_closure(
    root: Path,
    expected_lock_identity: str | None = None,
) -> PlatformRuntimeClosure:
    """Re-hash a collected closure without importing or executing its wheels."""

    if not isinstance(root, Path) or not root.is_absolute():
        raise _invalid_closure()
    lock_path = root.joinpath(*PurePosixPath(PLATFORM_RUNTIME_LOCK_PATH).parts)
    lock = PlatformWheelLock.from_bytes(
        _read_installed_file(lock_path, _MAX_LOCK_BYTES)
    )
    if expected_lock_identity is not None and lock.identity != expected_lock_identity:
        raise _invalid_closure()
    expected_wheels = {item.filename: item.sha256 for item in lock.wheels}
    wheel_root = root.joinpath(*PurePosixPath(PLATFORM_RUNTIME_WHEELHOUSE_ROOT).parts)
    if _installed_wheel_hashes(wheel_root) != expected_wheels:
        raise _invalid_closure()

    snapshots = tuple((item, wheel_root / item.filename) for item in lock.wheels)
    plans = tuple(_plan_wheel(item, source) for item, source in snapshots)
    _verify_dependency_closure(plans, lock)
    expected_site, expected_launcher = _expected_installed_files(plans)

    # Only the closure-owned coordinates are enumerated.  A staged release also
    # contains manifest.json and other native contracts which belong to the
    # bundle store and must neither be accepted as closure bytes nor rejected as
    # unexpected closure members.
    _verify_fixed_directory(root)
    for logical in (
        "payload",
        "payload/contracts",
        "payload/contracts/platform",
        PLATFORM_RUNTIME_ROOT,
        f"{PLATFORM_RUNTIME_ROOT}/bin",
        f"{PLATFORM_RUNTIME_ROOT}/lib",
        f"{PLATFORM_RUNTIME_ROOT}/lib/python3.12",
        PLATFORM_RUNTIME_SITE_PACKAGES,
        PLATFORM_RUNTIME_WHEELHOUSE_ROOT,
    ):
        _verify_fixed_directory(root.joinpath(*PurePosixPath(logical).parts))

    site_root = root.joinpath(*PurePosixPath(PLATFORM_RUNTIME_SITE_PACKAGES).parts)
    observed_site = _installed_subtree_files(
        site_root,
        PLATFORM_RUNTIME_SITE_PACKAGES,
        expected_paths=set(expected_site),
    )
    if set(observed_site) != set(expected_site):
        raise _invalid_closure()
    for logical, (expected_size, expected_sha256) in expected_site.items():
        observed = observed_site[logical]
        if (
            observed.mode != 0o444
            or observed.size_bytes != expected_size
            or observed.sha256 != expected_sha256
        ):
            raise _invalid_closure()

    bin_root = root.joinpath(*PurePosixPath(f"{PLATFORM_RUNTIME_ROOT}/bin").parts)
    observed_bin = _installed_subtree_files(
        bin_root,
        f"{PLATFORM_RUNTIME_ROOT}/bin",
        expected_paths={PLATFORM_RUNTIME_LAUNCHER_PATH},
    )
    if set(observed_bin) != {PLATFORM_RUNTIME_LAUNCHER_PATH}:
        raise _invalid_closure()
    launcher = observed_bin[PLATFORM_RUNTIME_LAUNCHER_PATH]
    if (
        launcher.mode != 0o555
        or launcher.size_bytes != expected_launcher[0]
        or launcher.sha256 != expected_launcher[1]
    ):
        raise _invalid_closure()

    lock_record = _installed_runtime_file(root, PLATFORM_RUNTIME_LOCK_PATH, 0o444)
    wheel_records = tuple(
        _installed_runtime_file(
            root,
            f"{PLATFORM_RUNTIME_WHEELHOUSE_ROOT}/{item.filename}",
            0o444,
        )
        for item in lock.wheels
    )
    records = (
        lock_record,
        *wheel_records,
        *tuple(observed_site[path] for path in sorted(observed_site)),
        launcher,
    )
    return PlatformRuntimeClosure(lock.identity, tuple(sorted(records)))


def _plan_wheel(locked: LockedWheel, source: Path) -> _WheelPlan:
    try:
        from packaging.tags import parse_tag
        from packaging.utils import canonicalize_name, parse_wheel_filename
        from packaging.version import Version
    except ImportError:
        raise fail(
            "PLATFORM_RUNTIME_BUILD_TOOL_UNAVAILABLE",
            "Platform runtime build tooling is unavailable.",
            component="platform",
        ) from None

    try:
        parsed_name, parsed_version, _build, filename_tags = parse_wheel_filename(
            locked.filename
        )
        with zipfile.ZipFile(source) as archive:
            files = tuple(item for item in archive.infolist() if not item.is_dir())
            if not 0 < len(files) <= _MAX_WHEEL_FILES:
                raise ValueError
            names = tuple(item.filename for item in files)
            if len(set(names)) != len(names) or len(
                {name.casefold() for name in names}
            ) != len(names):
                raise ValueError
            for item in files:
                _wheel_member_path(item.filename)
                mode = item.external_attr >> 16
                kind = stat.S_IFMT(mode)
                if item.flag_bits & 0x1 or (kind and kind != stat.S_IFREG):
                    raise ValueError
                if not 0 <= item.file_size <= _MAX_WHEEL_MEMBER_BYTES:
                    raise ValueError
            if sum(item.file_size for item in files) > _MAX_EXPANDED_BYTES:
                raise ValueError
            dist_info_roots = {
                PurePosixPath(name).parts[0]
                for name in names
                if len(PurePosixPath(name).parts) >= 2
                and PurePosixPath(name).parts[0].endswith(".dist-info")
            }
            if len(dist_info_roots) != 1:
                raise ValueError
            dist_info = next(iter(dist_info_roots))
            metadata_name = f"{dist_info}/METADATA"
            wheel_name = f"{dist_info}/WHEEL"
            record_name = f"{dist_info}/RECORD"
            if not {metadata_name, wheel_name, record_name}.issubset(names):
                raise ValueError
            metadata_content = _verified_zip_bytes(archive, metadata_name)
            wheel_content = _verified_zip_bytes(archive, wheel_name)
            record_content = _verified_zip_bytes(archive, record_name)
            _verify_wheel_record(archive, files, record_name, record_content)
            metadata = BytesParser(policy=compat32).parsebytes(metadata_content)
            wheel_metadata = BytesParser(policy=compat32).parsebytes(wheel_content)
            name = metadata.get("Name")
            version = metadata.get("Version")
            if (
                not isinstance(name, str)
                or _DISTRIBUTION.fullmatch(name) is None
                or not isinstance(version, str)
                or _VERSION.fullmatch(version) is None
                or canonicalize_name(name) != canonicalize_name(str(parsed_name))
                or Version(version) != parsed_version
                or wheel_metadata.get("Wheel-Version") != "1.0"
                or wheel_metadata.get("Root-Is-Purelib") not in {"true", "false"}
            ):
                raise ValueError
            declared_tags = frozenset(
                tag
                for raw in wheel_metadata.get_all("Tag", [])
                for tag in parse_tag(raw)
            )
            if declared_tags != frozenset(filename_tags) or not any(
                _supported_wheel_tag(tag.interpreter, tag.abi, tag.platform)
                for tag in declared_tags
            ):
                raise ValueError
            members = tuple(
                sorted(
                    (item.filename, _installed_member_path(item.filename))
                    for item in files
                )
            )
            if any(path.endswith(".pth") for _member, path in members):
                raise ValueError
            return _WheelPlan(
                locked=locked,
                source=source,
                distribution=_WheelDistribution(
                    canonicalize_name(name),
                    version,
                    metadata.get("Requires-Python"),
                    tuple(metadata.get_all("Requires-Dist", [])),
                ),
                members=members,
            )
    except DeploymentError:
        raise
    except Exception:
        raise _invalid_closure() from None


def _verify_dependency_closure(
    plans: tuple[_WheelPlan, ...], lock: PlatformWheelLock
) -> None:
    try:
        from packaging.requirements import Requirement
        from packaging.specifiers import SpecifierSet
        from packaging.utils import canonicalize_name
        from packaging.version import Version
    except ImportError:
        raise fail(
            "PLATFORM_RUNTIME_BUILD_TOOL_UNAVAILABLE",
            "Platform runtime build tooling is unavailable.",
            component="platform",
        ) from None
    distributions = {item.distribution.name: item.distribution for item in plans}
    if len(distributions) != len(plans) or "helixweave" not in distributions:
        raise _invalid_closure()
    environment = {
        "implementation_name": SUPPORTED_PYTHON_IMPLEMENTATION_NAME,
        # Patch-sensitive marker variables are rejected below.  Supplying the
        # supported floor keeps packaging's evaluator total without making
        # the build host's micro version a release fact.
        "implementation_version": lock.minimum_python_full_version,
        "os_name": "posix",
        "platform_machine": SUPPORTED_MACHINE,
        "platform_python_implementation": "CPython",
        "platform_system": "Linux",
        "python_full_version": lock.minimum_python_full_version,
        "python_version": lock.python_version,
        "sys_platform": SUPPORTED_SYSTEM,
    }
    allowed_marker_variables = frozenset(
        {
            "extra",
            "implementation_name",
            "os_name",
            "platform_machine",
            "platform_python_implementation",
            "platform_system",
            "python_version",
            "sys_platform",
        }
    )
    selected_extras: dict[str, set[str]] = {"helixweave": set(_PLATFORM_ROOT_EXTRAS)}
    queue: deque[str] = deque(("helixweave",))
    reachable: set[str] = set()
    while queue:
        name = queue.popleft()
        distribution = distributions.get(name)
        if distribution is None:
            raise _invalid_closure()
        reachable.add(name)
        try:
            supports_python = distribution.requires_python is None or (
                _supports_all_python_312_patches(
                    SpecifierSet(distribution.requires_python),
                    Version,
                )
            )
        except Exception:
            raise _invalid_closure() from None
        if not supports_python:
            raise _invalid_closure()
        extras = selected_extras.get(name, set())
        for raw in distribution.requirements:
            try:
                requirement = Requirement(raw)
            except Exception:
                raise _invalid_closure() from None
            if requirement.url is not None:
                raise _invalid_closure()
            if requirement.marker is not None and not _marker_variables(
                requirement.marker
            ).issubset(allowed_marker_variables):
                raise _invalid_closure()
            marker_extras = ("", *sorted(extras))
            if requirement.marker is not None and not any(
                requirement.marker.evaluate({**environment, "extra": extra})
                for extra in marker_extras
            ):
                continue
            dependency = canonicalize_name(requirement.name)
            installed = distributions.get(dependency)
            if (
                installed is None
                or Version(installed.version) not in requirement.specifier
            ):
                raise _invalid_closure()
            requested = set(requirement.extras)
            previous = selected_extras.setdefault(dependency, set())
            changed = not requested.issubset(previous)
            previous.update(requested)
            if dependency not in reachable or changed:
                queue.append(dependency)
    if reachable != set(distributions):
        raise _invalid_closure()


def _supports_all_python_312_patches(specifiers: object, version_type: object) -> bool:
    """Require one patch-neutral Python 3.12 support interval.

    The portable ``cp312`` closure may be produced and consumed by different
    security-patched CPython 3.12 interpreters.  Exclusions and arbitrary
    equality operators would make dependency selection change between those
    hosts, so they are rejected even when the current builder happens to pass.
    """

    try:
        items = tuple(specifiers)
        if any(item.operator in {"!=", "==="} for item in items):
            return False
        floor = version_type(SUPPORTED_PYTHON_MARKER_FLOOR)
        ceiling = version_type(SUPPORTED_PYTHON_MARKER_CEILING)
        return floor in specifiers and ceiling in specifiers
    except Exception:
        return False


def _marker_variables(marker: object) -> frozenset[str]:
    raw = getattr(marker, "_markers", None)
    variables: set[str] = set()

    def visit(value: object) -> None:
        if isinstance(value, list):
            for item in value:
                visit(item)
            return
        if isinstance(value, str) and value in {"and", "or"}:
            return
        if not isinstance(value, tuple) or len(value) != 3:
            raise _invalid_closure()
        left, operator, right = value
        if type(operator).__name__ != "Op":
            raise _invalid_closure()
        for operand in (left, right):
            kind = type(operand).__name__
            rendered = getattr(operand, "value", None)
            if kind not in {"Value", "Variable"} or not isinstance(rendered, str):
                raise _invalid_closure()
            if kind == "Variable":
                variables.add(rendered)

    visit(raw)
    return frozenset(variables)


def _extract_wheels(root: Path, plans: tuple[_WheelPlan, ...]) -> None:
    site_root = root.joinpath(*PurePosixPath(PLATFORM_RUNTIME_SITE_PACKAGES).parts)
    site_root.mkdir(parents=True, mode=0o700)
    targets: dict[str, str] = {}
    casefolded: dict[str, str] = {}
    for plan in plans:
        for member, relative in plan.members:
            prior = targets.get(relative)
            folded = casefolded.get(relative.casefold())
            if prior is not None or (folded is not None and folded != relative):
                raise _invalid_closure()
            targets[relative] = member
            casefolded[relative.casefold()] = relative
    for plan in plans:
        with zipfile.ZipFile(plan.source) as archive:
            for member, relative in plan.members:
                target = site_root.joinpath(*PurePosixPath(relative).parts)
                _copy_zip_member(archive, member, target)


def _install_candidate_launcher(root: Path, plans: tuple[_WheelPlan, ...]) -> None:
    platform_wheels = tuple(
        plan for plan in plans if plan.distribution.name == "helixweave"
    )
    if len(platform_wheels) != 1:
        raise _invalid_closure()
    platform_wheel = platform_wheels[0]
    if PLATFORM_RUNTIME_LAUNCHER_MEMBER not in {
        member for member, _relative in platform_wheel.members
    }:
        raise _invalid_closure()
    with zipfile.ZipFile(platform_wheel.source) as archive:
        content = _verified_zip_bytes(archive, PLATFORM_RUNTIME_LAUNCHER_MEMBER)
    target = root.joinpath(*PurePosixPath(PLATFORM_RUNTIME_LAUNCHER_PATH).parts)
    _write_bytes(target, content, 0o555)


def _expected_installed_files(
    plans: tuple[_WheelPlan, ...],
) -> tuple[dict[str, tuple[int, str]], tuple[int, str]]:
    expected: dict[str, tuple[int, str]] = {}
    casefolded: set[str] = set()
    launcher: tuple[int, str] | None = None
    platform_count = 0
    for plan in plans:
        if plan.distribution.name == "helixweave":
            platform_count += 1
        with zipfile.ZipFile(plan.source) as archive:
            for member, relative in plan.members:
                logical = f"{PLATFORM_RUNTIME_SITE_PACKAGES}/{relative}"
                if logical in expected or logical.casefold() in casefolded:
                    raise _invalid_closure()
                content = _verified_zip_bytes(archive, member)
                expected[logical] = (len(content), hashlib.sha256(content).hexdigest())
                casefolded.add(logical.casefold())
                if (
                    plan.distribution.name == "helixweave"
                    and member == PLATFORM_RUNTIME_LAUNCHER_MEMBER
                ):
                    launcher = (len(content), hashlib.sha256(content).hexdigest())
    if platform_count != 1 or launcher is None or not expected:
        raise _invalid_closure()
    return expected, launcher


def _verify_wheel_record(
    archive: zipfile.ZipFile,
    members: tuple[zipfile.ZipInfo, ...],
    record_name: str,
    record_content: bytes,
) -> None:
    try:
        rows = tuple(csv.reader(io.StringIO(record_content.decode("utf-8"))))
    except (UnicodeDecodeError, csv.Error):
        raise _invalid_closure() from None
    records: dict[str, tuple[str, str]] = {}
    for row in rows:
        if len(row) != 3 or row[0] in records:
            raise _invalid_closure()
        _wheel_member_path(row[0])
        records[row[0]] = (row[1], row[2])
    names = {item.filename for item in members}
    if set(records) != names:
        raise _invalid_closure()
    for item in members:
        digest, size = records[item.filename]
        if item.filename == record_name and digest == size == "":
            continue
        if not size.isdecimal() or int(size) != item.file_size:
            raise _invalid_closure()
        prefix = "sha256="
        if not digest.startswith(prefix):
            raise _invalid_closure()
        try:
            expected = base64.urlsafe_b64decode(
                digest.removeprefix(prefix)
                + "=" * (-len(digest.removeprefix(prefix)) % 4)
            )
        except (ValueError, TypeError):
            raise _invalid_closure() from None
        observed = hashlib.sha256()
        with archive.open(item) as source:
            while chunk := source.read(_READ_CHUNK):
                observed.update(chunk)
        if observed.digest() != expected:
            raise _invalid_closure()


def _installed_member_path(member: str) -> str:
    parts = _wheel_member_path(member).parts
    if len(parts) >= 3 and parts[0].endswith(".data"):
        if parts[1] not in {"purelib", "platlib"}:
            raise _invalid_closure()
        parts = parts[2:]
    if not parts:
        raise _invalid_closure()
    rendered = PurePosixPath(*parts).as_posix()
    if rendered.endswith(".pth"):
        raise _invalid_closure()
    return rendered


def _wheel_member_path(value: str) -> PurePosixPath:
    if (
        not isinstance(value, str)
        or not value
        or len(value) > 1024
        or value.startswith("/")
        or value.endswith("/")
        or "//" in value
        or "\\" in value
        or "\x00" in value
    ):
        raise _invalid_closure()
    path = PurePosixPath(value)
    if any(part in {"", ".", ".."} or len(part) > 255 for part in path.parts):
        raise _invalid_closure()
    return path


def _supported_wheel_tag(interpreter: str, abi: str, platform_tag: str) -> bool:
    if abi == "none" and platform_tag == "any":
        return interpreter in {"py3", "py312", "cp312"}
    if abi not in {SUPPORTED_PYTHON_ABI, "abi3"}:
        return False
    if abi == SUPPORTED_PYTHON_ABI and interpreter != SUPPORTED_PYTHON_ABI:
        return False
    if abi == "abi3":
        matched = re.fullmatch(r"cp3(\d+)", interpreter)
        if matched is None or not 7 <= int(matched.group(1)) <= 12:
            return False
    if platform_tag in _LEGACY_MANYLINUX:
        return True
    matched_platform = _MANYLINUX.fullmatch(platform_tag)
    return (
        matched_platform is not None
        and (int(matched_platform.group(1)), int(matched_platform.group(2)))
        <= MINIMUM_GLIBC
    )


def verify_supported_python_runtime(
    *,
    executable: Path | None = None,
    version_info: Sequence[int] | None = None,
    implementation_name: str | None = None,
    implementation_version: str | None = None,
    cache_tag: str | None = None,
    machine: str | None = None,
    system: str | None = None,
    libc_version: tuple[int, int] | None = None,
) -> str:
    """Fail closed unless the process is the one supported host ABI."""

    selected_executable = Path(sys.executable) if executable is None else executable
    selected_version = tuple(
        sys.version_info[:3] if version_info is None else version_info
    )
    selected_implementation_name = (
        sys.implementation.name if implementation_name is None else implementation_name
    )
    selected_implementation_version = (
        _running_implementation_version()
        if implementation_version is None
        else implementation_version
    )
    selected_cache_tag = (
        sys.implementation.cache_tag if cache_tag is None else cache_tag
    )
    selected_machine = platform.machine() if machine is None else machine
    selected_system = sys.platform if system is None else system
    selected_libc = _glibc_version() if libc_version is None else libc_version
    rendered_version = ".".join(str(item) for item in selected_version)
    if (
        selected_executable != SUPPORTED_PYTHON_EXECUTABLE
        or len(selected_version) != 3
        or any(
            not isinstance(item, int) or isinstance(item, bool) or item < 0
            for item in selected_version
        )
        or selected_version[:2] != SUPPORTED_PYTHON_VERSION
        or selected_version < SUPPORTED_PYTHON_MINIMUM_VERSION
        or selected_implementation_name != SUPPORTED_PYTHON_IMPLEMENTATION_NAME
        or selected_implementation_version != rendered_version
        or selected_cache_tag != SUPPORTED_PYTHON_CACHE_TAG
        or selected_machine != SUPPORTED_MACHINE
        or selected_system != SUPPORTED_SYSTEM
        or selected_libc < MINIMUM_GLIBC
    ):
        raise fail(
            "PLATFORM_PYTHON_INCOMPATIBLE",
            "Platform Python runtime is incompatible.",
            component="platform",
        )
    return canonical_identity(
        {
            "executable": str(selected_executable),
            "version": list(selected_version),
            "implementation_name": selected_implementation_name,
            "implementation_version": selected_implementation_version,
            "cache_tag": selected_cache_tag,
            "machine": selected_machine,
            "system": selected_system,
            "libc": list(selected_libc),
        },
        scheme=PLATFORM_RUNTIME_OBSERVATION_IDENTITY_SCHEME,
    )


def candidate_service_main(
    argv: Sequence[str] | None = None,
    *,
    release_root: Path | None = None,
    api_runner: Callable[[], int] | None = None,
    worker_runner: Callable[[], int] | None = None,
    encode_runtime_preparer: Callable[[], int] | None = None,
    bulk_runtime_preparer: Callable[[], int] | None = None,
    action_probe: Callable[[DeploymentActionRequest, Path], DeploymentActionReceipt]
    | None = None,
    action_request_path: Path = PLATFORM_ACTION_REQUEST_PATH,
    database_request_path: Path = DATABASE_PREPARE_REQUEST_PATH,
) -> int:
    """Run one fixed candidate action; never accept a module or filesystem path."""

    arguments = tuple(sys.argv[1:] if argv is None else argv)
    if len(arguments) != 1 or arguments[0] not in {
        "api",
        "worker",
        "db-prepare",
        "operator-action",
        "encode-runtime-prepare",
        "bulk-runtime-prepare",
    }:
        _write_public_error(
            "PLATFORM_SERVICE_ACTION_INVALID", "Platform service action is invalid."
        )
        return 64
    selected_root = _candidate_release_root() if release_root is None else release_root
    try:
        if arguments[0] == "api":
            return (_run_api if api_runner is None else api_runner)()
        if arguments[0] == "worker":
            return (_run_worker if worker_runner is None else worker_runner)()
        if arguments[0] == "encode-runtime-prepare":
            if encode_runtime_preparer is not None:
                return encode_runtime_preparer()
            from encode_pipeline.deployment.encode_runtime_materializer import (
                OfflineEncodeRuntimeMaterializer,
            )
            from encode_pipeline.deployment.layout import DeploymentLayout

            request = _parse_root_request(
                _read_root_request(ENCODE_RUNTIME_PREPARE_REQUEST_PATH),
                EncodeRuntimePrepareRequest.from_dict,
            )
            if request.authority_platform_identity != selected_root.name:
                raise fail(
                    "ENCODE_RUNTIME_PREPARE_REQUEST_MISMATCH",
                    "ENCODE runtime preparation request does not match this release.",
                )
            receipt = OfflineEncodeRuntimeMaterializer(
                DeploymentLayout.supported(),
                service_uid=os.geteuid(),
                service_gid=os.getegid(),
            ).prepare(request)
            if (
                not isinstance(receipt, EncodeRuntimePrepareReceipt)
                or receipt.request_identity != request.identity
                or receipt.deployment_identity != request.deployment_identity
            ):
                raise fail(
                    "ENCODE_RUNTIME_PREPARE_RECEIPT_INVALID",
                    "ENCODE runtime preparation receipt is invalid.",
                )
            sys.stdout.buffer.write(canonical_json_bytes(receipt.to_dict()))
            sys.stdout.buffer.flush()
            return 0
        if arguments[0] == "bulk-runtime-prepare":
            if bulk_runtime_preparer is not None:
                return bulk_runtime_preparer()
            from encode_pipeline.deployment.bulk_runtime_materializer import (
                OfflineBulkRuntimeMaterializer,
            )
            from encode_pipeline.deployment.layout import DeploymentLayout

            request = _parse_root_request(
                _read_root_request(BULK_RUNTIME_PREPARE_REQUEST_PATH),
                BulkRuntimePrepareRequest.from_dict,
            )
            if request.authority_platform_identity != selected_root.name:
                raise fail(
                    "BULK_RUNTIME_PREPARE_REQUEST_MISMATCH",
                    "Bulk runtime preparation request does not match this release.",
                )
            receipt = OfflineBulkRuntimeMaterializer(
                DeploymentLayout.supported()
            ).prepare(request)
            if (
                not isinstance(receipt, BulkRuntimePrepareReceipt)
                or receipt.request_identity != request.identity
                or receipt.candidate_bulk_identity != request.candidate_bulk_identity
            ):
                raise fail(
                    "BULK_RUNTIME_PREPARE_RECEIPT_INVALID",
                    "Bulk runtime preparation receipt is invalid.",
                )
            sys.stdout.buffer.write(canonical_json_bytes(receipt.to_dict()))
            sys.stdout.buffer.flush()
            return 0
        if arguments[0] == "db-prepare":
            request = _parse_root_request(
                _read_root_request(database_request_path),
                DatabasePrepareRequest.from_dict,
            )
            if request.deployment_identity != selected_root.name:
                raise fail(
                    "DB_PREPARE_REQUEST_MISMATCH",
                    "Database preparation request does not match this release.",
                )
            receipt = prepare_candidate_database(request)
            sys.stdout.buffer.write(canonical_json_bytes(receipt.to_dict()))
            sys.stdout.buffer.flush()
            return 0
        request = _parse_root_request(
            _read_root_request(action_request_path),
            DeploymentActionRequest.from_dict,
        )
        if (
            request.authority_platform_identity != selected_root.name
            or request.candidate_active["platform"] != selected_root.name
        ):
            raise fail(
                "OPERATOR_ACTION_REQUEST_MISMATCH",
                "Operator action request does not match this release.",
            )
        receipt = (
            _default_platform_action_probe(request, selected_root)
            if action_probe is None
            else action_probe(request, selected_root)
        )
        if (
            not isinstance(receipt, DeploymentActionReceipt)
            or receipt.request_identity != request.identity
            or receipt.status
            != ("admitted" if request.phase == "admit" else "observed")
        ):
            raise fail(
                "OPERATOR_ACTION_RECEIPT_INVALID",
                "Operator action receipt is invalid.",
            )
        sys.stdout.buffer.write(canonical_json_bytes(receipt.to_dict()))
        sys.stdout.buffer.flush()
        return 0
    except DeploymentError as error:
        _write_public_error(
            error.issue.code, error.issue.message, error.issue.recoverable
        )
        return 69 if error.issue.recoverable else 65
    except Exception:
        _write_public_error(
            "PLATFORM_SERVICE_FAILED", "Platform service action failed."
        )
        return 70


def prepare_candidate_database(
    request: DatabasePrepareRequest,
    *,
    layout: DeploymentLayout | None = None,
    expected_database_uid: int | None = None,
    expected_database_gid: int | None = None,
) -> DatabasePrepareReceipt:
    """Migrate only the request-selected fixed database coordinate."""

    if not isinstance(request, DatabasePrepareRequest):
        raise fail("DB_PREPARE_REQUEST_INVALID", "Database request is invalid.")
    try:
        from encode_pipeline.deployment.database import (
            database_content_identity,
            fresh_database_candidate_path,
            inspect_database,
        )
        from encode_pipeline.persistence.migration_admission import (
            verify_migration_execution_inventory,
        )
        from encode_pipeline.persistence.migrations import upgrade_database

        selected_layout = DeploymentLayout.supported() if layout is None else layout
        if not isinstance(selected_layout, DeploymentLayout):
            raise ValueError
        inventory = verify_migration_execution_inventory()
        if len(inventory.heads) != 1 or request.target_schema_heads != inventory.heads:
            raise ValueError
        if expected_database_uid is None or expected_database_gid is None:
            expected_database_uid, expected_database_gid = (
                _supported_database_ownership()
            )
        canonical = selected_layout.database
        if request.database_mode == "fresh-candidate":
            if canonical.exists() or canonical.is_symlink():
                raise ValueError
            database_path = fresh_database_candidate_path(
                selected_layout, request.task_identity
            )
            before_identity = None
        else:
            database_path = canonical
            if not database_path.exists() or database_path.is_symlink():
                raise ValueError
            before_identity = database_content_identity(
                inspect_database(
                    database_path,
                    expected_owner_uid=expected_database_uid,
                    expected_owner_gid=expected_database_gid,
                )
            )
        if database_path.is_symlink():
            raise ValueError
        database_url = f"sqlite:///{database_path}"
        if database_path.exists() and request.database_mode == "fresh-candidate":
            resumed = inspect_database(
                database_path,
                expected_owner_uid=expected_database_uid,
                expected_owner_gid=expected_database_gid,
            )
            if resumed.schema_heads != inventory.heads:
                raise ValueError
            return DatabasePrepareReceipt.create(
                request_identity=request.identity,
                database_before_identity=None,
                database_after_identity=database_content_identity(resumed),
                schema_heads=inventory.heads,
            )
        upgrade_database(database_url)
        inspection = inspect_database(
            database_path,
            expected_owner_uid=expected_database_uid,
            expected_owner_gid=expected_database_gid,
        )
        if inspection.schema_heads != inventory.heads:
            raise fail(
                "DATABASE_SCHEMA_INCOMPATIBLE",
                "Database schema is incompatible with this release.",
            )
        return DatabasePrepareReceipt.create(
            request_identity=request.identity,
            database_before_identity=before_identity,
            database_after_identity=database_content_identity(inspection),
            schema_heads=inventory.heads,
        )
    except DeploymentError:
        raise
    except Exception:
        raise fail(
            "DB_PREPARE_FAILED",
            "Database preparation failed.",
            recoverable=True,
        ) from None


def _migration_inventory_evidence(
    inventory: object,
) -> tuple[str, tuple[str, ...], tuple[str, ...]]:
    """Bind the verified candidate graph without restating migration facts."""

    from encode_pipeline.persistence.migration_admission import (
        VerifiedMigrationExecutionInventory,
    )

    if (
        not isinstance(inventory, VerifiedMigrationExecutionInventory)
        or len(inventory.heads) != 1
        or not isinstance(inventory.contract_sha256, str)
        or _SHA256.fullmatch(inventory.contract_sha256) is None
    ):
        raise ValueError
    target = inventory.heads[0]
    revisions = {item.revision: item for item in inventory.revisions}
    if target not in revisions or len(revisions) != len(inventory.revisions):
        raise ValueError
    ancestors: set[str] = set()
    pending = list(revisions[target].down_revision)
    while pending:
        revision = pending.pop()
        if revision in ancestors:
            continue
        item = revisions.get(revision)
        if item is None or revision == target:
            raise ValueError
        ancestors.add(revision)
        pending.extend(item.down_revision)
    if ancestors | {target} != set(revisions):
        raise ValueError
    return (
        f"sha256-{inventory.contract_sha256}",
        (target,),
        tuple(sorted(ancestors)),
    )


def _default_platform_action_probe(
    request: DeploymentActionRequest,
    release_root: Path,
) -> DeploymentActionReceipt:
    from encode_pipeline.deployment.admission import (
        resolved_facts_compatibility,
        validate_resolved_facts,
    )
    from encode_pipeline.deployment.bundle import BundleStore
    from encode_pipeline.deployment.database import (
        database_content_identity,
        inspect_database,
    )
    from encode_pipeline.deployment.layout import DeploymentLayout
    from encode_pipeline.deployment.models import (
        BULK_RNASEQ_RUNTIME,
        COMPONENTS,
        ENCODE_RUNTIME,
        PLATFORM,
    )
    from encode_pipeline.deployment.native_contracts import (
        PLATFORM_FRONTEND_CONTRACT,
        PLATFORM_PYTHON_RUNTIME_CONTRACT,
        PLATFORM_REFERENCES_CONTRACT,
        ProductionNativeContractResolver,
    )

    if not isinstance(request, DeploymentActionRequest):
        raise fail(
            "OPERATOR_ACTION_REQUEST_INVALID", "Operator action request is invalid."
        )
    layout = DeploymentLayout.supported()
    expected_authority = (
        layout.component_store(PLATFORM) / request.authority_platform_identity
    )
    if release_root != expected_authority:
        raise fail(
            "OPERATOR_ACTION_REQUEST_MISMATCH",
            "Operator action request does not match this release.",
        )

    readiness = {
        check: ReadinessCheck("not-applicable", "NOT_APPLICABLE")
        for check in VERIFICATION_CHECKS
    }
    native_checks = {
        PLATFORM: "platform-native",
        ENCODE_RUNTIME: "encode-runtime-native",
        BULK_RNASEQ_RUNTIME: "bulk-runtime-native",
    }
    native_identities: dict[str, str | None] = {
        component: None for component in COMPONENTS
    }
    facts_by_component = {}
    present_resolution_failed = False
    store = BundleStore(layout)
    resolver = ProductionNativeContractResolver()
    for component in COMPONENTS:
        identity = request.candidate_active[component]
        check_id = native_checks[component]
        if identity is None:
            readiness[check_id] = ReadinessCheck(
                "not-applicable", "COMPONENT_NOT_ACTIVE"
            )
            continue
        try:
            root = layout.component_store(component) / identity
            manifest = store.verify_installed(
                component,
                identity,
                expected_owner_uid=0,
                expected_owner_gid=0,
            )
            facts = resolver.resolve(root, manifest)
            validate_resolved_facts(manifest, facts)
            evidence = _resolved_facts_identity(facts)
            if component == ENCODE_RUNTIME:
                from encode_pipeline.deployment.operator_action import (
                    verify_materialized_encode_runtime,
                )

                materialized = verify_materialized_encode_runtime(layout, identity)
                evidence = canonical_identity(
                    {
                        "native_contracts": evidence,
                        "materialized_tree": materialized.tree_identity,
                        "entry_count": materialized.entry_count,
                    },
                    scheme="helixweave-encode-runtime-readiness-v1",
                )
            facts_by_component[component] = facts
            native_identities[component] = evidence
            readiness[check_id] = ReadinessCheck("ready", "READY", evidence)
        except Exception:
            present_resolution_failed = True
            readiness[check_id] = ReadinessCheck("not-ready", "CONTRACT_INVALID")

    platform_facts = facts_by_component.get(PLATFORM)
    if platform_facts is None:
        raise fail(
            "DEPLOYMENT_CONTRACT_ADMISSION_FAILED",
            "Deployment contract admission failed.",
            component=PLATFORM,
        )
    compatibility = (
        "incompatible"
        if present_resolution_failed
        else resolved_facts_compatibility(
            tuple(
                facts_by_component[component]
                for component in COMPONENTS
                if component in facts_by_component
            )
        )
    )

    try:
        from encode_pipeline.persistence.migration_admission import (
            verify_migration_execution_inventory,
        )

        migration_inventory = verify_migration_execution_inventory()
        (
            migration_inventory_identity,
            target_heads,
            known_schema_revisions,
        ) = _migration_inventory_evidence(migration_inventory)
        if target_heads != platform_facts.database_heads:
            raise ValueError
    except Exception:
        raise fail(
            "DEPLOYMENT_CONTRACT_ADMISSION_FAILED",
            "Deployment contract admission failed.",
            component=PLATFORM,
        ) from None

    wheel_lock_identity = _fact_contract_identity(
        platform_facts, PLATFORM_PYTHON_RUNTIME_CONTRACT
    )
    try:
        python_observation = verify_supported_python_runtime()
        python_identity = canonical_identity(
            {
                "abi_observation": python_observation,
                "wheel_lock": wheel_lock_identity,
            },
            scheme="helixweave-platform-python-readiness-v1",
        )
        readiness["python-runtime"] = ReadinessCheck("ready", "READY", python_identity)
    except Exception:
        readiness["python-runtime"] = ReadinessCheck("not-ready", "CONTRACT_INVALID")

    api_contract_sha256 = "0" * 64
    frontend_identity = None
    try:
        from encode_pipeline.frontend_assets import load_packaged_frontend_assets

        assets = load_packaged_frontend_assets()
        expected_frontend = _fact_contract_identity(
            platform_facts, PLATFORM_FRONTEND_CONTRACT
        )
        if assets.manifest.identity != expected_frontend:
            raise ValueError
        frontend_identity = assets.manifest.identity
        api_contract_sha256 = assets.manifest.api_contract_sha256
        readiness["frontend"] = ReadinessCheck("ready", "READY", frontend_identity)
    except Exception:
        readiness["frontend"] = ReadinessCheck("not-ready", "FRONTEND_INVALID")

    accepted_heads: tuple[str, ...] = ()
    database_identity = None
    try:
        database_uid, database_gid = _supported_database_ownership()
        inspection = inspect_database(
            SUPPORTED_DATABASE_PATH,
            expected_owner_uid=database_uid,
            expected_owner_gid=database_gid,
        )
        accepted_heads = inspection.schema_heads
        database_identity = database_content_identity(inspection)
        readiness["database-schema"] = ReadinessCheck(
            "ready" if accepted_heads == target_heads else "not-ready",
            "READY" if accepted_heads == target_heads else "SCHEMA_INCOMPATIBLE",
            database_identity if accepted_heads == target_heads else None,
        )
    except Exception:
        database_identity = None
        readiness["database-schema"] = ReadinessCheck(
            "unavailable", "SCHEMA_INCOMPATIBLE"
        )

    readiness["redis"] = _redis_ping(SUPPORTED_REDIS_SOCKET)
    readiness["docker"] = _docker_ping(SUPPORTED_DOCKER_SOCKET)
    readiness["configuration"] = ReadinessCheck("not-applicable", "NOT_APPLICABLE")
    readiness["permissions"] = _permission_readiness(layout, request.candidate_active)
    readiness["worker"] = ReadinessCheck("not-applicable", "NOT_APPLICABLE")
    readiness["references"] = _reference_readiness(SUPPORTED_REFERENCE_CONFIG)

    reference_compatibility_identity = _fact_contract_identity(
        platform_facts, PLATFORM_REFERENCES_CONTRACT
    )
    if any(
        readiness[check].status != "ready"
        for check in ("platform-native", "python-runtime", "frontend")
    ):
        compatibility = "incompatible"
    return DeploymentActionReceipt.create(
        request_identity=request.identity,
        status="admitted" if request.phase == "admit" else "observed",
        compatibility=compatibility,
        database_before_identity=(
            database_identity if request.phase == "admit" else None
        ),
        database_after_identity=(
            database_identity if request.phase == "observe" else None
        ),
        accepted_schema_heads=accepted_heads,
        target_schema_heads=target_heads,
        migration_inventory_identity=migration_inventory_identity,
        known_schema_revisions=known_schema_revisions,
        migration_required=accepted_heads != target_heads,
        rollback_supported=accepted_heads == target_heads,
        api_contract_sha256=api_contract_sha256,
        native_identities=native_identities,
        frontend_identity=frontend_identity,
        reference_compatibility_identity=reference_compatibility_identity,
        readiness=readiness,
    )


def _supported_database_ownership() -> tuple[int, int]:
    try:
        account = pwd.getpwnam("helixweave")
        group = grp.getgrnam("helixweave")
    except KeyError:
        raise fail(
            "PLATFORM_SERVICE_ACCOUNT_UNAVAILABLE",
            "Platform service account is unavailable.",
        ) from None
    if (
        account.pw_name != "helixweave"
        or group.gr_name != "helixweave"
        or account.pw_uid <= 0
        or group.gr_gid <= 0
        or account.pw_gid != group.gr_gid
    ):
        raise fail(
            "PLATFORM_SERVICE_ACCOUNT_UNAVAILABLE",
            "Platform service account is unavailable.",
        )
    return account.pw_uid, group.gr_gid


def _resolved_facts_identity(facts) -> str:
    return canonical_identity(
        {
            "component": facts.component,
            "deployment_identity": facts.deployment_identity,
            "version": facts.version,
            "contracts": [item.to_dict() for item in facts.contracts],
            "requirements": [item.to_dict() for item in facts.requirements],
            "database_heads": list(facts.database_heads),
        },
        scheme="helixweave-resolved-native-facts-evidence-v1",
    )


def _fact_contract_identity(facts, contract: str) -> str:
    matches = tuple(
        item.identity for item in facts.contracts if item.contract == contract
    )
    if len(matches) != 1:
        raise ValueError
    return matches[0]


def _redis_ping(path: Path) -> ReadinessCheck:
    return _unix_protocol_probe(
        path,
        request=b"*1\r\n$4\r\nPING\r\n",
        response_validator=lambda value: value == b"+PONG\r\n",
        unavailable_reason="REDIS_UNAVAILABLE",
        scheme="helixweave-redis-ping-observation-v1",
    )


def _docker_ping(path: Path) -> ReadinessCheck:
    def valid(response: bytes) -> bool:
        head, separator, body = response.partition(b"\r\n\r\n")
        status = head.split(b"\r\n", 1)[0]
        return (
            separator == b"\r\n\r\n"
            and status in {b"HTTP/1.0 200 OK", b"HTTP/1.1 200 OK"}
            and body.strip() == b"OK"
        )

    return _unix_protocol_probe(
        path,
        request=b"GET /_ping HTTP/1.0\r\nHost: localhost\r\nConnection: close\r\n\r\n",
        response_validator=valid,
        unavailable_reason="DOCKER_UNAVAILABLE",
        scheme="helixweave-docker-ping-observation-v1",
    )


def _unix_protocol_probe(
    path: Path,
    *,
    request: bytes,
    response_validator: Callable[[bytes], bool],
    unavailable_reason: str,
    scheme: str,
) -> ReadinessCheck:
    try:
        observed = path.lstat()
        if not stat.S_ISSOCK(observed.st_mode):
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
        if len(response) > 4096 or not response_validator(bytes(response)):
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


def _configuration_readiness(path: Path) -> ReadinessCheck:
    try:
        observed = path.lstat()
        if (
            not stat.S_ISREG(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or observed.st_nlink != 1
            or observed.st_uid != 0
            or stat.S_IMODE(observed.st_mode) != 0o640
            or not 0 < observed.st_size <= 64 * 1024
        ):
            raise OSError
        evidence = canonical_identity(
            {
                "sha256": _hash_file(path),
                "mode": stat.S_IMODE(observed.st_mode),
            },
            scheme="helixweave-platform-configuration-observation-v1",
        )
        return ReadinessCheck("ready", "READY", evidence)
    except OSError:
        return ReadinessCheck("not-ready", "CONFIGURATION_INVALID")


def _permission_readiness(layout, active: dict[str, str | None]) -> ReadinessCheck:
    try:
        evidence = []
        for component, identity in active.items():
            if identity is None:
                continue
            observed = (layout.component_store(component) / identity).lstat()
            if (
                not stat.S_ISDIR(observed.st_mode)
                or stat.S_ISLNK(observed.st_mode)
                or observed.st_uid != 0
                or observed.st_gid != 0
                or stat.S_IMODE(observed.st_mode) != 0o555
            ):
                raise OSError
            evidence.append(
                {
                    "component": component,
                    "deployment_identity": identity,
                    "owner_uid": observed.st_uid,
                    "owner_gid": observed.st_gid,
                    "mode": stat.S_IMODE(observed.st_mode),
                }
            )
        identity = canonical_identity(
            evidence,
            scheme="helixweave-candidate-permission-observation-v1",
        )
        return ReadinessCheck("ready", "READY", identity)
    except OSError:
        return ReadinessCheck("not-ready", "PERMISSION_INVALID")


def _reference_readiness(path: Path) -> ReadinessCheck:
    try:
        from encode_pipeline.adapters import EncodeStyleWorkflowAdapter
        from encode_pipeline.adapters.bulk_rnaseq import BulkRnaSeqWorkflowAdapter
        from encode_pipeline.services.private_reference_profiles import (
            load_private_reference_profile_config,
        )

        config = load_private_reference_profile_config(path)
        adapters = {
            adapter.metadata.workflow_id: adapter
            for adapter in (
                EncodeStyleWorkflowAdapter(),
                BulkRnaSeqWorkflowAdapter(),
            )
        }
        identities = []
        for config_key in config.config_keys:
            for workflow_id in config.workflow_ids_for(config_key):
                adapter = adapters.get(workflow_id)
                if adapter is None:
                    raise ValueError
                result = adapter.verify_reference_profile_binding(
                    config.binding_for(config_key, workflow_id)
                )
                if result.is_failure or result.value is None:
                    raise ValueError
                identities.append(
                    {
                        "workflow_id": result.value.workflow_id,
                        "contract_version": result.value.contract_version,
                        "identity_scheme": result.value.identity_scheme,
                        "identity_sha256": result.value.identity_sha256,
                    }
                )
        evidence = canonical_identity(
            identities,
            scheme="helixweave-reference-readiness-observation-v1",
        )
        return ReadinessCheck("ready", "READY", evidence)
    except Exception:
        return ReadinessCheck("not-ready", "REFERENCE_NOT_READY")


def _run_api() -> int:
    import uvicorn

    uvicorn.run(
        "encode_pipeline.deployment.frontend:create_app",
        factory=True,
        host="127.0.0.1",
        port=8000,
    )
    return 0


def _run_worker() -> int:
    from encode_pipeline.workers.cli import main

    return main(())


def _candidate_release_root() -> Path:
    current = Path(__file__).absolute()
    parts = current.parts
    marker = ("payload", "platform", "lib", "python3.12", "site-packages")
    for index in range(len(parts) - len(marker)):
        if tuple(parts[index : index + len(marker)]) == marker:
            root = Path(*parts[:index])
            if _valid_identity(root.name):
                return root
    raise fail(
        "PLATFORM_RELEASE_IDENTITY_INVALID",
        "Platform release identity is invalid.",
    )


def _read_root_request(path: Path) -> bytes:
    if not isinstance(path, Path) or not path.is_absolute():
        raise fail(
            "PLATFORM_SERVICE_REQUEST_UNTRUSTED", "Service request is untrusted."
        )
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
            or before.st_uid != 0
            or before.st_gid != os.getegid()
            or stat.S_IMODE(before.st_mode) != 0o640
            or not 0 < before.st_size <= _MAX_REQUEST_BYTES
        ):
            raise OSError
        content = bytearray()
        while chunk := os.read(
            descriptor, min(_READ_CHUNK, _MAX_REQUEST_BYTES + 1 - len(content))
        ):
            content.extend(chunk)
            if len(content) > _MAX_REQUEST_BYTES:
                raise OSError
        after = os.fstat(descriptor)
        if len(content) != before.st_size or _witness(before) != _witness(after):
            raise OSError
        return bytes(content)
    except OSError:
        raise fail(
            "PLATFORM_SERVICE_REQUEST_UNTRUSTED", "Service request is untrusted."
        ) from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _parse_root_request(content: bytes, parser: Callable[[object], object]):
    try:
        value = json.loads(
            content,
            object_pairs_hook=_unique_object,
            parse_constant=_reject_json_constant,
        )
        request = parser(value)
        if canonical_json_bytes(request.to_dict()) != content:
            raise ValueError
        return request
    except DeploymentError:
        raise
    except (
        AttributeError,
        TypeError,
        UnicodeDecodeError,
        ValueError,
        json.JSONDecodeError,
    ):
        raise fail(
            "PLATFORM_SERVICE_REQUEST_UNTRUSTED", "Service request is untrusted."
        ) from None


def _write_public_error(code: str, message: str, recoverable: bool = False) -> None:
    value = {
        "schema_version": "helixweave-platform-service-error-v1",
        "status": "error",
        "issue": {"code": code, "message": message, "recoverable": recoverable},
    }
    sys.stderr.write(canonical_json_bytes(value).decode("utf-8"))
    sys.stderr.flush()


def _verify_build_python() -> None:
    if (
        tuple(sys.version_info[:2]) != SUPPORTED_PYTHON_VERSION
        or tuple(sys.version_info[:3]) < SUPPORTED_PYTHON_MINIMUM_VERSION
        or sys.implementation.name != SUPPORTED_PYTHON_IMPLEMENTATION_NAME
        or _running_implementation_version()
        != ".".join(str(item) for item in sys.version_info[:3])
        or sys.implementation.cache_tag != SUPPORTED_PYTHON_CACHE_TAG
        or platform.machine() != SUPPORTED_MACHINE
        or sys.platform != SUPPORTED_SYSTEM
        or _glibc_version() < MINIMUM_GLIBC
    ):
        raise fail(
            "PLATFORM_RUNTIME_BUILD_ABI_INCOMPATIBLE",
            "Platform runtime build ABI is incompatible.",
            component="platform",
        )


def _running_implementation_version() -> str:
    value = sys.implementation.version
    releaselevel = getattr(value, "releaselevel", None)
    serial = getattr(value, "serial", None)
    if releaselevel != "final" or serial != 0:
        return "unsupported"
    return f"{value.major}.{value.minor}.{value.micro}"


def _glibc_version() -> tuple[int, int]:
    try:
        rendered = os.confstr("CS_GNU_LIBC_VERSION")
        if not isinstance(rendered, str):
            raise ValueError
        matched = re.fullmatch(r"glibc (\d+)\.(\d+)(?:\.\d+)?", rendered)
        if matched is None:
            raise ValueError
        return int(matched.group(1)), int(matched.group(2))
    except (OSError, ValueError):
        return (0, 0)


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError
        value[key] = item
    return value


def _reject_json_constant(_value: str) -> object:
    raise ValueError


def _valid_identity(value: object) -> bool:
    return isinstance(value, str) and _IDENTITY.fullmatch(value) is not None


def _source_directory(path: Path) -> Path:
    if not isinstance(path, Path) or not path.is_absolute():
        raise _invalid_closure()
    try:
        observed = path.lstat()
    except OSError:
        raise _invalid_closure() from None
    if (
        not stat.S_ISDIR(observed.st_mode)
        or stat.S_ISLNK(observed.st_mode)
        or observed.st_mode & 0o022
    ):
        raise _invalid_closure()
    return path


def _new_destination(path: Path) -> Path:
    if (
        not isinstance(path, Path)
        or not path.is_absolute()
        or path.name in {"", ".", ".."}
        or path.exists()
        or path.is_symlink()
    ):
        raise _invalid_closure()
    parent = _source_directory(path.parent)
    return parent / path.name


def _preflight_closure_capacity(
    parent: Path, *, file_count: int, payload_bytes: int
) -> None:
    declared = payload_bytes + MAX_MANIFEST_BYTES
    if (
        isinstance(file_count, bool)
        or isinstance(payload_bytes, bool)
        or not 0 < file_count <= MAX_BUNDLE_FILES
        or payload_bytes < 0
        or declared > MAX_BUNDLE_BYTES
    ):
        raise _invalid_closure()
    try:
        available = shutil.disk_usage(parent).free
    except OSError:
        raise fail(
            "DEPLOYMENT_STORAGE_UNAVAILABLE",
            "Deployment storage is unavailable.",
            component="platform",
            recoverable=True,
        ) from None
    reserve = max(MIN_FREE_SPACE_RESERVE, declared // 20)
    if available < declared + reserve:
        raise fail(
            "DEPLOYMENT_CAPACITY_INSUFFICIENT",
            "Deployment storage capacity is insufficient.",
            component="platform",
            recoverable=True,
        )


def _read_regular_source(path: Path, maximum: int) -> bytes:
    if not isinstance(path, Path) or not path.is_absolute():
        raise _invalid_closure()
    descriptor = -1
    try:
        before_path = path.lstat()
        descriptor = os.open(
            path,
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
        )
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_mode & 0o022
            or not 0 < before.st_size <= maximum
            or _witness(before_path) != _witness(before)
        ):
            raise OSError
        content = bytearray()
        while chunk := os.read(
            descriptor, min(_READ_CHUNK, maximum + 1 - len(content))
        ):
            content.extend(chunk)
            if len(content) > maximum:
                raise OSError
        after = os.fstat(descriptor)
        if len(content) != before.st_size or _witness(before) != _witness(after):
            raise OSError
        return bytes(content)
    except OSError:
        raise _invalid_closure() from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _wheelhouse_inventory(root: Path) -> dict[str, Path]:
    try:
        entries = tuple(sorted(root.iterdir(), key=lambda item: item.name))
    except OSError:
        raise _invalid_closure() from None
    inventory: dict[str, Path] = {}
    for entry in entries:
        if _WHEEL_FILENAME.fullmatch(entry.name) is None:
            raise _invalid_lock()
        observed = entry.lstat()
        if (
            not stat.S_ISREG(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or observed.st_nlink != 1
            or observed.st_mode & 0o022
            or not 0 < observed.st_size <= _MAX_WHEEL_BYTES
        ):
            raise _invalid_closure()
        inventory[entry.name] = entry
    return inventory


def _copy_hashed_source(source: Path, destination: Path, expected_sha256: str) -> None:
    content = _read_regular_source(source, _MAX_WHEEL_BYTES)
    if hashlib.sha256(content).hexdigest() != expected_sha256:
        raise _invalid_lock()
    _write_bytes(destination, content, 0o444)


def _write_bytes(path: Path, content: bytes, mode: int) -> None:
    path.parent.mkdir(parents=True, mode=0o700, exist_ok=True)
    flags = (
        os.O_WRONLY
        | os.O_CREAT
        | os.O_EXCL
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    descriptor = -1
    try:
        descriptor = os.open(path, flags, 0o600)
        offset = 0
        while offset < len(content):
            written = os.write(descriptor, content[offset:])
            if written <= 0:
                raise OSError
            offset += written
        os.fchmod(descriptor, mode)
        os.fsync(descriptor)
    except OSError:
        raise _invalid_closure() from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _copy_zip_member(archive: zipfile.ZipFile, member: str, target: Path) -> None:
    target.parent.mkdir(parents=True, mode=0o700, exist_ok=True)
    flags = (
        os.O_WRONLY
        | os.O_CREAT
        | os.O_EXCL
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    descriptor = -1
    try:
        descriptor = os.open(target, flags, 0o600)
        total = 0
        with archive.open(member) as source:
            while chunk := source.read(_READ_CHUNK):
                total += len(chunk)
                if total > _MAX_WHEEL_MEMBER_BYTES:
                    raise OSError
                offset = 0
                while offset < len(chunk):
                    written = os.write(descriptor, chunk[offset:])
                    if written <= 0:
                        raise OSError
                    offset += written
        os.fchmod(descriptor, 0o444)
        os.fsync(descriptor)
    except OSError:
        raise _invalid_closure() from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _verified_zip_bytes(archive: zipfile.ZipFile, member: str) -> bytes:
    info = archive.getinfo(member)
    if not 0 <= info.file_size <= _MAX_WHEEL_MEMBER_BYTES:
        raise _invalid_closure()
    content = archive.read(info)
    if len(content) != info.file_size:
        raise _invalid_closure()
    return content


def _freeze_tree(root: Path) -> None:
    directories: list[Path] = []
    for directory, child_directories, filenames in os.walk(root, topdown=True):
        child_directories.sort()
        filenames.sort()
        current = Path(directory)
        directories.append(current)
        for name in filenames:
            path = current / name
            expected = (
                0o555
                if path.relative_to(root).as_posix() == PLATFORM_RUNTIME_LAUNCHER_PATH
                else 0o444
            )
            os.chmod(path, expected, follow_symlinks=False)
    for directory in reversed(directories):
        os.chmod(directory, 0o555, follow_symlinks=False)
        _fsync_directory(directory)


def _logical_path(relative: str) -> str:
    path = _wheel_member_path(relative)
    logical = path.as_posix()
    if not logical.startswith("payload/"):
        raise _invalid_closure()
    return logical


def _read_installed_file(path: Path, maximum: int) -> bytes:
    try:
        observed = path.lstat()
        if (
            not stat.S_ISREG(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or observed.st_nlink != 1
            or stat.S_IMODE(observed.st_mode) != 0o444
            or not 0 < observed.st_size <= maximum
        ):
            raise OSError
        content = path.read_bytes()
        if len(content) != observed.st_size:
            raise OSError
        return content
    except OSError:
        raise _invalid_closure() from None


def _verify_fixed_directory(path: Path) -> None:
    try:
        observed = path.lstat()
    except OSError:
        raise _invalid_closure() from None
    if (
        not stat.S_ISDIR(observed.st_mode)
        or stat.S_ISLNK(observed.st_mode)
        or stat.S_IMODE(observed.st_mode) != 0o555
    ):
        raise _invalid_closure()


def _installed_runtime_file(
    root: Path,
    logical_path: str,
    expected_mode: int,
) -> PlatformRuntimeFile:
    logical = _logical_path(logical_path)
    source = root.joinpath(*PurePosixPath(logical).parts)
    try:
        observed = source.lstat()
    except OSError:
        raise _invalid_closure() from None
    if (
        not stat.S_ISREG(observed.st_mode)
        or stat.S_ISLNK(observed.st_mode)
        or observed.st_nlink != 1
        or stat.S_IMODE(observed.st_mode) != expected_mode
        or not 0 <= observed.st_size <= _MAX_WHEEL_BYTES
    ):
        raise _invalid_closure()
    return PlatformRuntimeFile(
        logical_path=logical,
        source=source,
        mode=expected_mode,
        size_bytes=observed.st_size,
        sha256=_hash_file(source),
    )


def _installed_subtree_files(
    subtree: Path,
    logical_root: str,
    *,
    expected_paths: set[str],
) -> dict[str, PlatformRuntimeFile]:
    _verify_fixed_directory(subtree)
    result: dict[str, PlatformRuntimeFile] = {}
    casefolded: set[str] = set()
    observed_directories = {logical_root}
    try:
        iterator = os.walk(subtree, topdown=True, followlinks=False)
        for directory, directories, filenames in iterator:
            directories.sort()
            filenames.sort()
            current = Path(directory)
            _verify_fixed_directory(current)
            for name in directories:
                _verify_fixed_directory(current / name)
                relative = (current / name).relative_to(subtree).as_posix()
                observed_directories.add(_logical_path(f"{logical_root}/{relative}"))
            for name in filenames:
                source = current / name
                relative = source.relative_to(subtree).as_posix()
                logical = _logical_path(f"{logical_root}/{relative}")
                if logical in result or logical.casefold() in casefolded:
                    raise _invalid_closure()
                observed = source.lstat()
                if (
                    not stat.S_ISREG(observed.st_mode)
                    or stat.S_ISLNK(observed.st_mode)
                    or observed.st_nlink != 1
                    or stat.S_IMODE(observed.st_mode) not in {0o444, 0o555}
                    or not 0 <= observed.st_size <= _MAX_WHEEL_MEMBER_BYTES
                ):
                    raise _invalid_closure()
                result[logical] = PlatformRuntimeFile(
                    logical_path=logical,
                    source=source,
                    mode=stat.S_IMODE(observed.st_mode),
                    size_bytes=observed.st_size,
                    sha256=_hash_file(source),
                )
                casefolded.add(logical.casefold())
    except DeploymentError:
        raise
    except OSError:
        raise _invalid_closure() from None
    expected_directories = {logical_root}
    for logical in expected_paths:
        parent = PurePosixPath(logical).parent
        root_path = PurePosixPath(logical_root)
        while parent != root_path:
            expected_directories.add(parent.as_posix())
            parent = parent.parent
        expected_directories.add(logical_root)
    if observed_directories != expected_directories:
        raise _invalid_closure()
    return result


def _installed_wheel_hashes(root: Path) -> dict[str, str]:
    _verify_fixed_directory(root)
    try:
        entries = tuple(sorted(root.iterdir(), key=lambda item: item.name))
    except OSError:
        raise _invalid_closure() from None
    result: dict[str, str] = {}
    for path in entries:
        if _WHEEL_FILENAME.fullmatch(path.name) is None:
            raise _invalid_closure()
        observed = path.lstat()
        if (
            not stat.S_ISREG(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or observed.st_nlink != 1
            or stat.S_IMODE(observed.st_mode) != 0o444
            or not 0 < observed.st_size <= _MAX_WHEEL_BYTES
        ):
            raise _invalid_closure()
        result[path.name] = _hash_file(path)
    return result


def _hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as source:
            while chunk := source.read(_READ_CHUNK):
                digest.update(chunk)
    except OSError:
        raise _invalid_closure() from None
    return digest.hexdigest()


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


def _discard_partial(path: Path) -> None:
    if not path.exists() or path.is_symlink():
        return
    try:
        for directory, directories, filenames in os.walk(path, topdown=False):
            current = Path(directory)
            os.chmod(current, 0o700, follow_symlinks=False)
            for name in filenames:
                child = current / name
                os.chmod(child, 0o600, follow_symlinks=False)
                child.unlink()
            for name in directories:
                child = current / name
                os.chmod(child, 0o700, follow_symlinks=False)
                child.rmdir()
        path.rmdir()
    except OSError:
        # The caller still receives a stable failure; doctor can report an
        # interrupted build directory if a hostile local race prevented cleanup.
        pass


def _witness(value: os.stat_result) -> tuple[int, ...]:
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


__all__ = [
    "BULK_RUNTIME_PREPARE_REQUEST_PATH",
    "DATABASE_PREPARE_REQUEST_PATH",
    "PLATFORM_ACTION_REQUEST_PATH",
    "PLATFORM_RUNTIME_LAUNCHER_PATH",
    "PLATFORM_RUNTIME_LOCK_PATH",
    "PLATFORM_RUNTIME_SITE_PACKAGES",
    "PLATFORM_RUNTIME_WHEELHOUSE_ROOT",
    "SUPPORTED_PYTHON_EXECUTABLE",
    "SUPPORTED_PYTHON_IMPLEMENTATION_NAME",
    "SUPPORTED_PYTHON_MINIMUM_VERSION",
    "SUPPORTED_PYTHON_VERSION",
    "BulkRuntimePrepareReceipt",
    "BulkRuntimePrepareRequest",
    "DatabasePrepareReceipt",
    "DatabasePrepareRequest",
    "DeploymentActionReceipt",
    "DeploymentActionRequest",
    "LockedWheel",
    "PlatformRuntimeClosure",
    "PlatformRuntimeFile",
    "PlatformWheelLock",
    "ReadinessCheck",
    "VERIFICATION_CHECKS",
    "build_platform_wheel_lock",
    "candidate_service_main",
    "collect_platform_runtime_closure",
    "inspect_platform_runtime_closure",
    "prepare_candidate_database",
    "verify_supported_python_runtime",
]
