"""Read-only verification for the root-owned stable operator boundary."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat


OPERATOR_CLOSURE_SCHEMA = "helixweave-operator-closure-v1"
STABLE_BOUNDARY_IDENTITY_SCHEME = "helixweave-stable-operator-boundary-v1"
SUPPORTED_OPERATOR_ROOT = Path("/opt/helixweave/operator")
_IDENTITY = re.compile(r"sha256-[0-9a-f]{64}\Z")
_SHA256 = re.compile(r"[0-9a-f]{64}\Z")
_MAX_MANIFEST_BYTES = 1024 * 1024
_MAX_BOUNDARY_BYTES = 1024 * 1024
_EXPECTED_BOUNDARY_COUNT = 8


class StableBoundaryError(RuntimeError):
    """Private, path-free stable-boundary verification failure."""


@dataclass(frozen=True)
class _BoundaryRecord:
    snapshot: PurePosixPath
    installed: PurePosixPath
    mode: int
    sha256: str
    size_bytes: int


def verify_stable_operator_boundary(
    *,
    operator_root: Path = SUPPORTED_OPERATOR_ROOT,
    host_root: Path = Path("/"),
    expected_uid: int = 0,
    expected_gid: int = 0,
) -> str:
    """Return one aggregate identity after verifying every stable entry."""
    try:
        _validate_arguments(operator_root, host_root, expected_uid, expected_gid)
        _require_directory(
            operator_root,
            expected_uid=expected_uid,
            expected_gid=expected_gid,
            modes={0o755},
        )
        current = operator_root / "current"
        before = current.lstat()
        target = os.readlink(current)
        after = current.lstat()
        if (
            not stat.S_ISLNK(before.st_mode)
            or _stat_identity(before) != _stat_identity(after)
            or before.st_uid != expected_uid
            or before.st_gid != expected_gid
            or _IDENTITY.fullmatch(target) is None
        ):
            raise StableBoundaryError
        closure = operator_root / target
        _require_directory(
            closure,
            expected_uid=expected_uid,
            expected_gid=expected_gid,
            modes={0o555},
        )
        manifest_content = _read_owned_file(
            closure / "closure.json",
            expected_uid=expected_uid,
            expected_gid=expected_gid,
            expected_mode=0o444,
            expected_size=None,
            expected_sha256=None,
            maximum=_MAX_MANIFEST_BYTES,
        )
        manifest = json.loads(manifest_content, object_pairs_hook=_unique_object)
        records = _parse_manifest(manifest, target)
        if _canonical_bytes(manifest) != manifest_content:
            raise StableBoundaryError

        evidence: list[dict[str, object]] = []
        for record in records:
            snapshot_path = closure.joinpath(*record.snapshot.parts)
            installed_path = _host_path(host_root, record.installed)
            snapshot_content = _read_owned_file(
                snapshot_path,
                expected_uid=expected_uid,
                expected_gid=expected_gid,
                expected_mode=record.mode,
                expected_size=record.size_bytes,
                expected_sha256=record.sha256,
                maximum=_MAX_BOUNDARY_BYTES,
            )
            installed_content = _read_owned_file(
                installed_path,
                expected_uid=expected_uid,
                expected_gid=expected_gid,
                expected_mode=record.mode,
                expected_size=record.size_bytes,
                expected_sha256=record.sha256,
                maximum=_MAX_BOUNDARY_BYTES,
            )
            if installed_content != snapshot_content:
                raise StableBoundaryError
            evidence.append(
                {
                    "installed_path": record.installed.as_posix(),
                    "mode": record.mode,
                    "sha256": record.sha256,
                    "size_bytes": record.size_bytes,
                }
            )
        return _canonical_identity(
            {
                "closure_identity": target,
                "boundaries": evidence,
            },
            scheme=STABLE_BOUNDARY_IDENTITY_SCHEME,
        )
    except StableBoundaryError:
        raise
    except (OSError, ValueError, UnicodeError, json.JSONDecodeError):
        raise StableBoundaryError from None


def _parse_manifest(raw: object, expected_identity: str) -> tuple[_BoundaryRecord, ...]:
    if (
        not isinstance(raw, dict)
        or set(raw) != {"schema_version", "identity", "files"}
        or raw["schema_version"] != OPERATOR_CLOSURE_SCHEMA
        or raw["identity"] != expected_identity
        or not isinstance(raw["files"], list)
    ):
        raise StableBoundaryError
    identity_document = {
        "schema_version": raw["schema_version"],
        "files": raw["files"],
    }
    calculated = (
        f"sha256-{hashlib.sha256(_identity_bytes(identity_document)).hexdigest()}"
    )
    if calculated != expected_identity:
        raise StableBoundaryError

    all_paths: list[str] = []
    boundaries: list[_BoundaryRecord] = []
    installed_paths: set[str] = set()
    for item in raw["files"]:
        if not isinstance(item, dict):
            raise StableBoundaryError
        base_keys = {"mode", "path", "sha256", "size_bytes"}
        if frozenset(item) not in {
            frozenset(base_keys),
            frozenset(base_keys | {"installed_path"}),
        }:
            raise StableBoundaryError
        path = _relative_path(item["path"])
        mode = item["mode"]
        size_bytes = item["size_bytes"]
        digest = item["sha256"]
        if (
            not isinstance(mode, int)
            or isinstance(mode, bool)
            or not 0 <= mode <= 0o777
            or not isinstance(size_bytes, int)
            or isinstance(size_bytes, bool)
            or not 0 <= size_bytes <= 2**63 - 1
            or not isinstance(digest, str)
            or _SHA256.fullmatch(digest) is None
        ):
            raise StableBoundaryError
        all_paths.append(path.as_posix())
        installed_value = item.get("installed_path")
        if installed_value is None:
            continue
        installed = _installed_path(installed_value)
        if (
            len(path.parts) != 2
            or path.parts[0] != "boundary"
            or installed.as_posix() in installed_paths
        ):
            raise StableBoundaryError
        installed_paths.add(installed.as_posix())
        boundaries.append(_BoundaryRecord(path, installed, mode, digest, size_bytes))
    if (
        all_paths != sorted(all_paths)
        or len(set(all_paths)) != len(all_paths)
        or len(boundaries) != _EXPECTED_BOUNDARY_COUNT
    ):
        raise StableBoundaryError
    boundaries.sort(key=lambda item: item.installed.as_posix())
    return tuple(boundaries)


def _validate_arguments(
    operator_root: Path,
    host_root: Path,
    expected_uid: int,
    expected_gid: int,
) -> None:
    if (
        not isinstance(operator_root, Path)
        or not operator_root.is_absolute()
        or not isinstance(host_root, Path)
        or not host_root.is_absolute()
        or any(
            not isinstance(value, int) or isinstance(value, bool) or value < 0
            for value in (expected_uid, expected_gid)
        )
    ):
        raise StableBoundaryError


def _relative_path(value: object) -> PurePosixPath:
    if not isinstance(value, str) or not value or len(value.encode("utf-8")) > 1024:
        raise StableBoundaryError
    path = PurePosixPath(value)
    if (
        path.is_absolute()
        or path.as_posix() != value
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise StableBoundaryError
    return path


def _installed_path(value: object) -> PurePosixPath:
    if not isinstance(value, str) or not value or len(value.encode("utf-8")) > 1024:
        raise StableBoundaryError
    path = PurePosixPath(value)
    if (
        not path.is_absolute()
        or path.as_posix() != value
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise StableBoundaryError
    parent = path.parent.as_posix()
    if not (
        (parent == "/usr/libexec" and path.name.startswith("helixweave-"))
        or (parent == "/etc/sudoers.d" and path.name == "helixweave-operator")
    ):
        raise StableBoundaryError
    return path


def _host_path(host_root: Path, absolute: PurePosixPath) -> Path:
    return host_root.joinpath(*absolute.parts[1:])


def _require_directory(
    path: Path,
    *,
    expected_uid: int,
    expected_gid: int,
    modes: set[int],
) -> None:
    observed = path.lstat()
    if (
        not stat.S_ISDIR(observed.st_mode)
        or stat.S_ISLNK(observed.st_mode)
        or observed.st_uid != expected_uid
        or observed.st_gid != expected_gid
        or stat.S_IMODE(observed.st_mode) not in modes
    ):
        raise StableBoundaryError


def _read_owned_file(
    path: Path,
    *,
    expected_uid: int,
    expected_gid: int,
    expected_mode: int,
    expected_size: int | None,
    expected_sha256: str | None,
    maximum: int,
) -> bytes:
    _require_directory(
        path.parent,
        expected_uid=expected_uid,
        expected_gid=expected_gid,
        modes={0o555, 0o755, 0o750, 0o700},
    )
    descriptor = os.open(
        path,
        os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_uid != expected_uid
            or before.st_gid != expected_gid
            or stat.S_IMODE(before.st_mode) != expected_mode
            or not 0 < before.st_size <= maximum
            or (expected_size is not None and before.st_size != expected_size)
        ):
            raise StableBoundaryError
        content = _read_bounded(descriptor, maximum)
        after = os.fstat(descriptor)
        at_path = path.stat(follow_symlinks=False)
        if (
            len(content) != before.st_size
            or _stat_identity(before) != _stat_identity(after)
            or (after.st_dev, after.st_ino) != (at_path.st_dev, at_path.st_ino)
            or (
                expected_sha256 is not None
                and hashlib.sha256(content).hexdigest() != expected_sha256
            )
        ):
            raise StableBoundaryError
        return content
    finally:
        os.close(descriptor)


def _read_bounded(descriptor: int, maximum: int) -> bytes:
    chunks: list[bytes] = []
    remaining = maximum + 1
    while remaining:
        chunk = os.read(descriptor, min(64 * 1024, remaining))
        if not chunk:
            break
        chunks.append(chunk)
        remaining -= len(chunk)
    content = b"".join(chunks)
    if len(content) > maximum:
        raise StableBoundaryError
    return content


def _unique_object(pairs):
    value = {}
    for key, item in pairs:
        if key in value:
            raise ValueError
        value[key] = item
    return value


def _canonical_bytes(value: object) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("utf-8")


def _identity_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        allow_nan=False,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


def _canonical_identity(value: object, *, scheme: str) -> str:
    digest = hashlib.sha256()
    for item in (scheme.encode("ascii"), _canonical_bytes(value)):
        digest.update(len(item).to_bytes(8, byteorder="big", signed=False))
        digest.update(item)
    return f"sha256-{digest.hexdigest()}"


def _stat_identity(observed: os.stat_result) -> tuple[int, ...]:
    return (
        observed.st_dev,
        observed.st_ino,
        observed.st_mode,
        observed.st_nlink,
        observed.st_uid,
        observed.st_gid,
        observed.st_size,
        observed.st_mtime_ns,
        observed.st_ctime_ns,
    )


__all__ = [
    "STABLE_BOUNDARY_IDENTITY_SCHEME",
    "StableBoundaryError",
    "verify_stable_operator_boundary",
]
