"""Stdlib-only identity witness for the fixed supported Docker boundary."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import hashlib
import os
from pathlib import Path
import stat

from encode_pipeline.deployment.canonical import canonical_identity
from encode_pipeline.deployment.errors import fail


DOCKER_EXECUTABLE = Path("/usr/bin/docker")
DOCKER_SOCKET = Path("/run/helixweave/docker/docker.sock")
DOCKER_CLIENT_IDENTITY_SCHEME = "helixweave-bulk-docker-client-identity-v1"
DOCKER_ENDPOINT_IDENTITY_SCHEME = "helixweave-bulk-docker-endpoint-identity-v1"

_MAX_DOCKER_CLIENT_BYTES = 256 * 1024 * 1024
_MAX_PROC_NET_UNIX_BYTES = 4 * 1024 * 1024
_READ_CHUNK_BYTES = 1024 * 1024
_PROC_NET_UNIX = Path("/proc/net/unix")


@dataclass(frozen=True)
class BulkDockerBoundary:
    """Path-free identities for the fixed Docker client and daemon socket."""

    client_identity: str
    endpoint_identity: str
    daemon_uid: int
    daemon_gid: int


def observe_bulk_docker_boundary(
    daemon_uid: int,
    daemon_gid: int,
    *,
    _client_path: Path = DOCKER_EXECUTABLE,
    _socket_path: Path = DOCKER_SOCKET,
    _proc_unix_reader: Callable[[], bytes] | None = None,
    _socket_lstat: Callable[[Path], os.stat_result] | None = None,
    _client_owner_uid: int = 0,
    _client_owner_gid: int = 0,
) -> BulkDockerBoundary:
    """Observe the fixed client bytes and filesystem/kernel socket identity."""

    if (
        any(
            not isinstance(value, int)
            or isinstance(value, bool)
            or not 0 <= value <= 2**32 - 1
            for value in (
                daemon_uid,
                daemon_gid,
                _client_owner_uid,
                _client_owner_gid,
            )
        )
        or not isinstance(_client_path, Path)
        or not _client_path.is_absolute()
        or not isinstance(_socket_path, Path)
        or not _socket_path.is_absolute()
        or (_socket_lstat is not None and not callable(_socket_lstat))
    ):
        raise _boundary_invalid()
    descriptor = -1
    try:
        client_before = os.lstat(_client_path)
        descriptor = os.open(
            _client_path,
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
        )
        client_opened = os.fstat(descriptor)
        if (
            not stat.S_ISREG(client_before.st_mode)
            or client_before.st_nlink != 1
            or client_before.st_uid != _client_owner_uid
            or client_before.st_gid != _client_owner_gid
            or client_before.st_mode & 0o111 == 0
            or client_before.st_mode & 0o022
            or not 0 < client_before.st_size <= _MAX_DOCKER_CLIENT_BYTES
            or _stat_witness(client_before) != _stat_witness(client_opened)
        ):
            raise OSError
        client_sha256 = _hash_descriptor(
            descriptor, expected_size=client_opened.st_size
        )
        client_after = os.fstat(descriptor)
        if _stat_witness(client_after) != _stat_witness(client_opened):
            raise OSError
    except OSError:
        raise _boundary_invalid() from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)

    try:
        socket_observer = os.lstat if _socket_lstat is None else _socket_lstat
        socket_before = socket_observer(_socket_path)
        socket_mode = stat.S_IMODE(socket_before.st_mode)
        if (
            not stat.S_ISSOCK(socket_before.st_mode)
            or socket_before.st_nlink != 1
            or socket_before.st_uid != daemon_uid
            or socket_before.st_gid != daemon_gid
            or socket_mode not in {0o600, 0o660}
        ):
            raise OSError
        proc_content = (
            _read_proc_net_unix() if _proc_unix_reader is None else _proc_unix_reader()
        )
        kernel_inode = _socket_kernel_inode(proc_content, _socket_path)
        socket_after = socket_observer(_socket_path)
        if _stat_witness(socket_before) != _stat_witness(socket_after):
            raise OSError
    except (OSError, TypeError, ValueError):
        raise _boundary_invalid() from None

    client_identity = canonical_identity(
        {
            "kind": "regular",
            "dev": client_opened.st_dev,
            "ino": client_opened.st_ino,
            "mode": stat.S_IMODE(client_opened.st_mode),
            "nlink": client_opened.st_nlink,
            "uid": client_opened.st_uid,
            "gid": client_opened.st_gid,
            "size": client_opened.st_size,
            "ctime_ns": client_opened.st_ctime_ns,
            "sha256": client_sha256,
        },
        scheme=DOCKER_CLIENT_IDENTITY_SCHEME,
    )
    endpoint_identity = canonical_identity(
        {
            "kind": "socket",
            "dev": socket_after.st_dev,
            "ino": socket_after.st_ino,
            "mode": stat.S_IMODE(socket_after.st_mode),
            "nlink": socket_after.st_nlink,
            "uid": socket_after.st_uid,
            "gid": socket_after.st_gid,
            "ctime_ns": socket_after.st_ctime_ns,
            "kernel_inode": kernel_inode,
        },
        scheme=DOCKER_ENDPOINT_IDENTITY_SCHEME,
    )
    return BulkDockerBoundary(
        client_identity,
        endpoint_identity,
        daemon_uid,
        daemon_gid,
    )


def _hash_descriptor(descriptor: int, *, expected_size: int) -> str:
    os.lseek(descriptor, 0, os.SEEK_SET)
    digest = hashlib.sha256()
    total = 0
    while chunk := os.read(descriptor, _READ_CHUNK_BYTES):
        total += len(chunk)
        if total > expected_size:
            raise OSError
        digest.update(chunk)
    if total != expected_size:
        raise OSError
    os.lseek(descriptor, 0, os.SEEK_SET)
    return digest.hexdigest()


def _read_proc_net_unix() -> bytes:
    descriptor = -1
    try:
        descriptor = os.open(
            _PROC_NET_UNIX,
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
        )
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode) or before.st_uid != 0:
            raise OSError
        content = bytearray()
        while chunk := os.read(
            descriptor,
            min(
                _READ_CHUNK_BYTES,
                _MAX_PROC_NET_UNIX_BYTES + 1 - len(content),
            ),
        ):
            content.extend(chunk)
            if len(content) > _MAX_PROC_NET_UNIX_BYTES:
                raise OSError
        after = os.fstat(descriptor)
        if _stat_witness(before) != _stat_witness(after):
            raise OSError
        return bytes(content)
    except OSError:
        raise _boundary_invalid() from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _socket_kernel_inode(content: bytes, socket_path: Path) -> int:
    try:
        lines = content.decode("ascii", errors="strict").splitlines()
    except UnicodeDecodeError:
        raise ValueError from None
    if not lines or lines[0].split() != [
        "Num",
        "RefCount",
        "Protocol",
        "Flags",
        "Type",
        "St",
        "Inode",
        "Path",
    ]:
        raise ValueError
    matches: list[int] = []
    for line in lines[1:]:
        fields = line.split(maxsplit=7)
        if len(fields) != 8 or fields[7] != str(socket_path):
            continue
        if fields[3] != "00010000" or fields[4] != "0001" or fields[5] != "01":
            raise ValueError
        if not fields[6].isdigit() or int(fields[6]) <= 0:
            raise ValueError
        matches.append(int(fields[6]))
    if len(matches) != 1:
        raise ValueError
    return matches[0]


def _stat_witness(value: os.stat_result) -> tuple[int, ...]:
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


def _boundary_invalid():
    return fail(
        "BULK_RUNTIME_DOCKER_BOUNDARY_INVALID",
        "Bulk runtime Docker boundary is invalid.",
        component="bulk-rnaseq-runtime",
        recoverable=True,
    )


__all__ = [
    "BulkDockerBoundary",
    "DOCKER_CLIENT_IDENTITY_SCHEME",
    "DOCKER_ENDPOINT_IDENTITY_SCHEME",
    "DOCKER_EXECUTABLE",
    "DOCKER_SOCKET",
    "observe_bulk_docker_boundary",
]
