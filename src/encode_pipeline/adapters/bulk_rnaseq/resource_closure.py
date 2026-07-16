"""Fail-closed local resource closure verification for bulk RNA-seq.

This module deliberately owns filesystem inspection for adapter-private runtime
resources.  Authoring validation remains purely structural; callers invoke
these functions only while composing an execution workspace.
"""

from __future__ import annotations

from contextlib import suppress
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
from typing import Any

from encode_pipeline.platform.results import Issue, Result


RIBO_DATABASE_CLOSURE_SCHEME = "sha256-framed-ribo-database-v1"
SORTMERNA_INDEX_CLOSURE_SCHEME = "sha256-framed-sortmerna-index-v1"
SORTMERNA_INDEX_MANIFEST_FILENAME = "helixweave-sortmerna-index.json"
SORTMERNA_INDEX_MANIFEST_SCHEMA_VERSION = "1.0.0"
SORTMERNA_VERSION = "4.3.7"
SORTMERNA_NO_PREBUILT_INDEX_STRATEGY = (
    "nfcore-rnaseq-3.26.0-single-database-fastq-sortmerna-index-v1"
)
# Compatibility aliases for the PR150-era typed contract name. The file is a
# strict content manifest, not merely an unverified binding sidecar.
SORTMERNA_INDEX_BINDING_FILENAME = SORTMERNA_INDEX_MANIFEST_FILENAME
SORTMERNA_INDEX_BINDING_SCHEMA_VERSION = SORTMERNA_INDEX_MANIFEST_SCHEMA_VERSION

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_SAFE_VERSION = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._+-]{0,63}$")
_SAFE_EXECUTION_PATH_COMPONENT = re.compile(r"^[A-Za-z0-9._+@%=-]+$")
_RRNA_FASTA_NAME = re.compile(r"^[A-Za-z0-9._+@%=-]+\.(?:fa|fasta|fna)(?:\.gz)?$")
_URI_SCHEME = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]*:")
_READ_CHUNK_BYTES = 1024 * 1024


@dataclass(frozen=True)
class ResourceClosurePolicy:
    """Explicit resource bounds for local rRNA assets."""

    maximum_manifest_bytes: int = 1024 * 1024
    maximum_manifest_entries: int = 128
    maximum_relative_path_bytes: int = 1024
    maximum_input_file_bytes: int = 256 * 1024**3
    maximum_database_file_bytes: int = 32 * 1024**3
    maximum_database_total_bytes: int = 128 * 1024**3
    maximum_index_entries: int = 4096
    maximum_index_depth: int = 32
    maximum_index_file_bytes: int = 32 * 1024**3
    maximum_index_total_bytes: int = 256 * 1024**3
    maximum_binding_bytes: int = 64 * 1024

    def __post_init__(self) -> None:
        for name, value in vars(self).items():
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ValueError(f"{name} must be a positive integer")


DEFAULT_RESOURCE_CLOSURE_POLICY = ResourceClosurePolicy()


@dataclass(frozen=True)
class RiboDatabaseFile:
    """One regular database file observed relative to the manifest."""

    manifest_entry: str
    path: Path
    size_bytes: int
    sha256: str


@dataclass(frozen=True)
class VerifiedLocalFile:
    """Identity of one safely opened, stable local regular file."""

    path: Path
    size_bytes: int
    sha256: str


@dataclass(frozen=True)
class VerifiedLocalFileBytes:
    """Bounded bytes and identity from one stable no-follow file read."""

    path: Path
    size_bytes: int
    sha256: str
    content: bytes


@dataclass(frozen=True)
class NoPrebuiltSortMeRnaIndexRoute:
    """Deterministic composition evidence for an unqualified scientific route.

    One verified database file, the immutable upstream source, and an exact
    container permit one fixed index-building route to be composed.  This does
    not claim that a real scientific container acceptance gate has passed.
    """

    database_closure_sha256: str
    strategy: str = SORTMERNA_NO_PREBUILT_INDEX_STRATEGY
    sortmerna_version: str = SORTMERNA_VERSION
    deterministic_database_inputs: bool = True
    deterministic_composition_accepted: bool = False
    requires_runtime_index_build: bool = True
    scientific_runtime_accepted: bool = False


@dataclass(frozen=True)
class RiboDatabaseClosure:
    """Verified manifest content and every regular file it references."""

    manifest_sha256: str
    identity_sha256: str
    files: tuple[RiboDatabaseFile, ...]
    no_prebuilt_sortmerna_index: NoPrebuiltSortMeRnaIndexRoute


@dataclass(frozen=True)
class SortMeRnaIndexEntry:
    """One deterministic directory or regular-file entry in an index tree."""

    relative_path: str
    kind: str
    size_bytes: int | None = None
    sha256: str | None = None


@dataclass(frozen=True)
class SortMeRnaIndexClosure:
    """Verified SortMeRNA index tree bound to one database closure."""

    path: Path
    manifest_sha256: str
    identity_sha256: str
    database_closure_sha256: str
    sortmerna_version: str
    entries: tuple[SortMeRnaIndexEntry, ...]


@dataclass(frozen=True)
class _PathObservation:
    parts: tuple[str, ...]
    identities: tuple[tuple[int, ...], ...]
    is_directory: bool


@dataclass(frozen=True)
class _FileObservation:
    size_bytes: int
    sha256: str
    content: bytes
    path_observation: _PathObservation


class _ResourceFailure(Exception):
    """Internal control-flow exception whose details are never user-facing."""


class _DuplicateJsonKey(ValueError):
    pass


class _DescriptorRoot:
    """Read descendants through an already-open, symlink-safe directory."""

    def __init__(self, root_path: Path, descriptor: int) -> None:
        self.root_path = root_path
        self._descriptor = descriptor
        self._root_identity = _directory_identity(os.fstat(descriptor))
        self._closed = False

    @classmethod
    def open(cls, path: Path) -> "_DescriptorRoot":
        descriptor = _open_absolute_directory(path)
        return cls(path, descriptor)

    def __enter__(self) -> "_DescriptorRoot":
        return self

    def __exit__(self, *_args: object) -> None:
        self.close()

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        with suppress(OSError):
            os.close(self._descriptor)

    @property
    def descriptor(self) -> int:
        return self._descriptor

    def read_regular_file(
        self,
        parts: tuple[str, ...],
        *,
        maximum_bytes: int,
        capture: bool,
    ) -> _FileObservation:
        descriptors: list[int] = []
        try:
            current = self._descriptor
            identities: list[tuple[int, ...]] = [self._root_identity]
            for component in parts[:-1]:
                listed = os.stat(component, dir_fd=current, follow_symlinks=False)
                descriptor = os.open(component, _directory_flags(), dir_fd=current)
                descriptors.append(descriptor)
                opened = os.fstat(descriptor)
                if (
                    not stat.S_ISDIR(listed.st_mode)
                    or not stat.S_ISDIR(opened.st_mode)
                    or _directory_identity(listed) != _directory_identity(opened)
                ):
                    raise _ResourceFailure
                current = descriptor
                identities.append(_directory_identity(opened))

            listed = os.stat(parts[-1], dir_fd=current, follow_symlinks=False)
            if not stat.S_ISREG(listed.st_mode) or listed.st_nlink != 1:
                raise _ResourceFailure
            descriptor = os.open(parts[-1], _file_flags(), dir_fd=current)
            descriptors.append(descriptor)
            before = os.fstat(descriptor)
            if (
                not stat.S_ISREG(before.st_mode)
                or before.st_nlink != 1
                or _file_identity(listed) != _file_identity(before)
                or before.st_size > maximum_bytes
            ):
                raise _ResourceFailure

            digest = hashlib.sha256()
            chunks: list[bytes] = []
            remaining = before.st_size
            while remaining:
                chunk = os.read(descriptor, min(_READ_CHUNK_BYTES, remaining))
                if not chunk:
                    raise _ResourceFailure
                digest.update(chunk)
                if capture:
                    chunks.append(chunk)
                remaining -= len(chunk)
            if os.read(descriptor, 1):
                raise _ResourceFailure
            after = os.fstat(descriptor)
            if _file_identity(before) != _file_identity(after):
                raise _ResourceFailure
            identities.append(_file_identity(after))
            path_observation = _PathObservation(
                parts=parts,
                identities=tuple(identities),
                is_directory=False,
            )
            if self.reopen_identity(parts, is_directory=False) != path_observation:
                raise _ResourceFailure
            return _FileObservation(
                size_bytes=after.st_size,
                sha256=digest.hexdigest(),
                content=b"".join(chunks),
                path_observation=path_observation,
            )
        except _ResourceFailure:
            raise
        except (OSError, TypeError, ValueError) as error:
            raise _ResourceFailure from error
        finally:
            for descriptor in reversed(descriptors):
                with suppress(OSError):
                    os.close(descriptor)

    def reopen_identity(
        self,
        parts: tuple[str, ...],
        *,
        is_directory: bool,
    ) -> _PathObservation:
        descriptors: list[int] = []
        try:
            current = self._descriptor
            identities: list[tuple[int, ...]] = [
                _directory_identity(os.fstat(self._descriptor))
            ]
            for position, component in enumerate(parts):
                final = position == len(parts) - 1
                if final and not is_directory:
                    descriptor = os.open(component, _file_flags(), dir_fd=current)
                    info = os.fstat(descriptor)
                    if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
                        raise _ResourceFailure
                    identity = _file_identity(info)
                else:
                    descriptor = os.open(component, _directory_flags(), dir_fd=current)
                    info = os.fstat(descriptor)
                    if not stat.S_ISDIR(info.st_mode):
                        raise _ResourceFailure
                    identity = _directory_identity(info)
                descriptors.append(descriptor)
                current = descriptor
                identities.append(identity)
            return _PathObservation(parts, tuple(identities), is_directory)
        except _ResourceFailure:
            raise
        except (OSError, TypeError, ValueError) as error:
            raise _ResourceFailure from error
        finally:
            for descriptor in reversed(descriptors):
                with suppress(OSError):
                    os.close(descriptor)

    def verify_stable(self, observations: list[_PathObservation]) -> None:
        if _directory_identity(os.fstat(self._descriptor)) != self._root_identity:
            raise _ResourceFailure
        for observation in observations:
            reopened = self.reopen_identity(
                observation.parts,
                is_directory=observation.is_directory,
            )
            if reopened != observation:
                raise _ResourceFailure
        reopened_root = _open_absolute_directory(self.root_path)
        try:
            if _directory_identity(os.fstat(reopened_root)) != self._root_identity:
                raise _ResourceFailure
        finally:
            os.close(reopened_root)


@dataclass
class _IndexWalkState:
    policy: ResourceClosurePolicy
    entries: list[SortMeRnaIndexEntry]
    observations: list[_PathObservation]
    total_bytes: int = 0
    binding_content: bytes | None = None


def safe_regular_file_identity(
    path: str | Path,
    *,
    expected_sha256: str | None = None,
    maximum_bytes: int | None = None,
    policy: ResourceClosurePolicy = DEFAULT_RESOURCE_CLOSURE_POLICY,
) -> Result[VerifiedLocalFile]:
    """Stream-hash one local regular file through a symlink-safe descriptor.

    ``maximum_bytes`` lets a server raise or lower the explicit per-file bound
    for a configured deployment.  It is never inferred from submitted data.
    """

    try:
        _require_policy(policy)
        if expected_sha256 is not None:
            _require_sha256(expected_sha256)
        effective_maximum = (
            policy.maximum_input_file_bytes if maximum_bytes is None else maximum_bytes
        )
        if (
            isinstance(effective_maximum, bool)
            or not isinstance(effective_maximum, int)
            or effective_maximum <= 0
        ):
            raise _ResourceFailure
        file_path = _absolute_lexical_path(path)
        file_name = _safe_leaf_name(file_path.name)
        with _DescriptorRoot.open(file_path.parent) as reader:
            observed = reader.read_regular_file(
                (file_name,),
                maximum_bytes=effective_maximum,
                capture=False,
            )
            if expected_sha256 is not None and observed.sha256 != expected_sha256:
                raise _ResourceFailure
            reader.verify_stable([observed.path_observation])
        return Result.success(
            VerifiedLocalFile(
                path=file_path,
                size_bytes=observed.size_bytes,
                sha256=observed.sha256,
            )
        )
    except (_ResourceFailure, OSError, TypeError, ValueError, UnicodeError):
        return _failure("BULK_RNASEQ_LOCAL_FILE_INVALID")


def safe_regular_file_bytes(
    path: str | Path,
    *,
    maximum_bytes: int,
    expected_sha256: str | None = None,
    policy: ResourceClosurePolicy = DEFAULT_RESOURCE_CLOSURE_POLICY,
) -> Result[VerifiedLocalFileBytes]:
    """Read bounded local metadata bytes through the same no-follow boundary."""
    try:
        _require_policy(policy)
        if expected_sha256 is not None:
            _require_sha256(expected_sha256)
        if (
            isinstance(maximum_bytes, bool)
            or not isinstance(maximum_bytes, int)
            or maximum_bytes <= 0
        ):
            raise _ResourceFailure
        file_path = _absolute_lexical_path(path)
        file_name = _safe_leaf_name(file_path.name)
        with _DescriptorRoot.open(file_path.parent) as reader:
            observed = reader.read_regular_file(
                (file_name,),
                maximum_bytes=maximum_bytes,
                capture=True,
            )
            if expected_sha256 is not None and observed.sha256 != expected_sha256:
                raise _ResourceFailure
            reader.verify_stable([observed.path_observation])
        return Result.success(
            VerifiedLocalFileBytes(
                path=file_path,
                size_bytes=observed.size_bytes,
                sha256=observed.sha256,
                content=observed.content,
            )
        )
    except (_ResourceFailure, OSError, TypeError, ValueError, UnicodeError):
        return _failure("BULK_RNASEQ_LOCAL_FILE_INVALID")


def verify_ribo_database_manifest(
    path: str | Path,
    *,
    expected_manifest_sha256: str,
    policy: ResourceClosurePolicy = DEFAULT_RESOURCE_CLOSURE_POLICY,
) -> Result[RiboDatabaseClosure]:
    """Verify a bounded local rRNA manifest and its complete file closure."""

    try:
        _require_policy(policy)
        _require_sha256(expected_manifest_sha256)
        manifest_path = _absolute_lexical_path(path)
        manifest_name = _safe_leaf_name(manifest_path.name)
        with _DescriptorRoot.open(manifest_path.parent) as reader:
            manifest = reader.read_regular_file(
                (manifest_name,),
                maximum_bytes=policy.maximum_manifest_bytes,
                capture=True,
            )
            if manifest.sha256 != expected_manifest_sha256:
                raise _ResourceFailure
            entries = _parse_manifest(manifest.content, policy=policy)
            observations = [manifest.path_observation]
            files: list[RiboDatabaseFile] = []
            total_bytes = 0
            for entry, parts in entries:
                observed = reader.read_regular_file(
                    parts,
                    maximum_bytes=policy.maximum_database_file_bytes,
                    capture=False,
                )
                total_bytes += observed.size_bytes
                if total_bytes > policy.maximum_database_total_bytes:
                    raise _ResourceFailure
                observations.append(observed.path_observation)
                files.append(
                    RiboDatabaseFile(
                        manifest_entry=entry,
                        path=reader.root_path.joinpath(*parts),
                        size_bytes=observed.size_bytes,
                        sha256=observed.sha256,
                    )
                )
            reader.verify_stable(observations)

        identity = _database_closure_identity(manifest.sha256, files)
        route = NoPrebuiltSortMeRnaIndexRoute(
            database_closure_sha256=identity,
            deterministic_composition_accepted=len(files) == 1,
        )
        return Result.success(
            RiboDatabaseClosure(
                manifest_sha256=manifest.sha256,
                identity_sha256=identity,
                files=tuple(files),
                no_prebuilt_sortmerna_index=route,
            )
        )
    except (_ResourceFailure, OSError, TypeError, ValueError, UnicodeError):
        return _failure("BULK_RNASEQ_RIBO_DATABASE_INVALID")


def verify_sortmerna_index(
    path: str | Path,
    *,
    expected_index_sha256: str,
    database_closure: RiboDatabaseClosure,
    policy: ResourceClosurePolicy = DEFAULT_RESOURCE_CLOSURE_POLICY,
) -> Result[SortMeRnaIndexClosure]:
    """Verify an exact sidecar-manifested SortMeRNA index tree and binding."""

    try:
        _require_policy(policy)
        if not isinstance(database_closure, RiboDatabaseClosure):
            raise _ResourceFailure
        _require_sha256(expected_index_sha256)
        _require_sha256(database_closure.identity_sha256)
        index_path = _absolute_lexical_path(path)
        with _DescriptorRoot.open(index_path) as reader:
            manifest = reader.read_regular_file(
                (SORTMERNA_INDEX_MANIFEST_FILENAME,),
                maximum_bytes=policy.maximum_binding_bytes,
                capture=True,
            )
            if manifest.sha256 != expected_index_sha256:
                raise _ResourceFailure
            binding = _parse_index_binding(manifest.content, policy=policy)
            state = _IndexWalkState(
                policy=policy,
                entries=[],
                observations=[manifest.path_observation],
            )
            _walk_index_directory(reader, reader.descriptor, (), state)
            if (
                state.binding_content is None
                or state.binding_content != manifest.content
            ):
                raise _ResourceFailure
            reader.verify_stable(state.observations)

        if binding["database_closure_sha256"] != database_closure.identity_sha256:
            raise _ResourceFailure
        if binding["sortmerna_version"] != SORTMERNA_VERSION:
            raise _ResourceFailure
        observed_files = tuple(
            sorted(
                (
                    entry
                    for entry in state.entries
                    if entry.kind == "file"
                    and entry.relative_path != SORTMERNA_INDEX_MANIFEST_FILENAME
                ),
                key=lambda entry: entry.relative_path,
            )
        )
        observed_directories = {
            entry.relative_path for entry in state.entries if entry.kind == "directory"
        }
        declared_files = binding["files"]
        if len(observed_files) != len(declared_files):
            raise _ResourceFailure
        for observed, declared in zip(observed_files, declared_files, strict=True):
            if (
                observed.relative_path != declared["path"]
                or observed.size_bytes != declared["size_bytes"]
                or observed.sha256 != declared["sha256"]
            ):
                raise _ResourceFailure
        expected_directories = {
            PurePosixPath(item["path"]).parent.as_posix()
            for item in declared_files
            if PurePosixPath(item["path"]).parent.as_posix() != "."
        }
        for item in tuple(expected_directories):
            parent = PurePosixPath(item).parent
            while parent.as_posix() != ".":
                expected_directories.add(parent.as_posix())
                parent = parent.parent
        if observed_directories != expected_directories:
            raise _ResourceFailure
        identity = _index_closure_identity(manifest.sha256, observed_files)
        return Result.success(
            SortMeRnaIndexClosure(
                path=index_path,
                manifest_sha256=manifest.sha256,
                identity_sha256=identity,
                database_closure_sha256=database_closure.identity_sha256,
                sortmerna_version=binding["sortmerna_version"],
                entries=observed_files,
            )
        )
    except (_ResourceFailure, OSError, TypeError, ValueError, UnicodeError):
        return _failure("BULK_RNASEQ_SORTMERNA_INDEX_INVALID")


def _walk_index_directory(
    reader: _DescriptorRoot,
    descriptor: int,
    prefix: tuple[str, ...],
    state: _IndexWalkState,
) -> None:
    try:
        names = sorted(os.listdir(descriptor))
    except OSError as error:
        raise _ResourceFailure from error
    for raw_name in names:
        name = _safe_leaf_name(raw_name)
        parts = (*prefix, name)
        relative_path = PurePosixPath(*parts).as_posix()
        if (
            len(state.entries) >= state.policy.maximum_index_entries
            or len(parts) > state.policy.maximum_index_depth
            or len(relative_path.encode("utf-8"))
            > state.policy.maximum_relative_path_bytes
        ):
            raise _ResourceFailure
        try:
            listed = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
        except OSError as error:
            raise _ResourceFailure from error
        if stat.S_ISDIR(listed.st_mode):
            child = -1
            try:
                child = os.open(name, _directory_flags(), dir_fd=descriptor)
                before = os.fstat(child)
                if _directory_identity(listed) != _directory_identity(before):
                    raise _ResourceFailure
                state.entries.append(
                    SortMeRnaIndexEntry(relative_path=relative_path, kind="directory")
                )
                _walk_index_directory(reader, child, parts, state)
                after = os.fstat(child)
                if _directory_identity(before) != _directory_identity(after):
                    raise _ResourceFailure
                observation = reader.reopen_identity(parts, is_directory=True)
                if observation.identities[-1] != _directory_identity(after):
                    raise _ResourceFailure
                state.observations.append(observation)
            except _ResourceFailure:
                raise
            except OSError as error:
                raise _ResourceFailure from error
            finally:
                if child >= 0:
                    with suppress(OSError):
                        os.close(child)
            continue
        if not stat.S_ISREG(listed.st_mode) or listed.st_nlink != 1:
            raise _ResourceFailure

        maximum_bytes = state.policy.maximum_index_file_bytes
        capture = parts == (SORTMERNA_INDEX_BINDING_FILENAME,)
        if capture:
            maximum_bytes = min(maximum_bytes, state.policy.maximum_binding_bytes)
            if state.binding_content is not None:
                raise _ResourceFailure
        observed = reader.read_regular_file(
            parts,
            maximum_bytes=maximum_bytes,
            capture=capture,
        )
        if observed.path_observation.identities[-1] != _file_identity(listed):
            raise _ResourceFailure
        state.total_bytes += observed.size_bytes
        if state.total_bytes > state.policy.maximum_index_total_bytes:
            raise _ResourceFailure
        state.entries.append(
            SortMeRnaIndexEntry(
                relative_path=relative_path,
                kind="file",
                size_bytes=observed.size_bytes,
                sha256=observed.sha256,
            )
        )
        state.observations.append(observed.path_observation)
        if capture:
            state.binding_content = observed.content


def _parse_manifest(
    content: bytes,
    *,
    policy: ResourceClosurePolicy,
) -> tuple[tuple[str, tuple[str, ...]], ...]:
    try:
        text = content.decode("utf-8", errors="strict")
    except UnicodeDecodeError as error:
        raise _ResourceFailure from error
    if (
        not text
        or text.startswith("\ufeff")
        or any(separator in text for separator in ("\x85", "\u2028", "\u2029"))
    ):
        raise _ResourceFailure
    lines = text.splitlines()
    if not lines or len(lines) > policy.maximum_manifest_entries:
        raise _ResourceFailure
    entries: list[tuple[str, tuple[str, ...]]] = []
    seen: set[tuple[str, ...]] = set()
    staged_names: set[str] = set()
    logical_names: set[str] = set()
    for line in lines:
        if not line or line != line.strip():
            raise _ResourceFailure
        parts = _safe_relative_parts(
            line,
            maximum_bytes=policy.maximum_relative_path_bytes,
        )
        staged_name = parts[-1].removesuffix(".gz")
        logical_name = staged_name.rsplit(".", 1)[0]
        if (
            len(parts) > policy.maximum_index_depth
            or parts in seen
            or _RRNA_FASTA_NAME.fullmatch(parts[-1]) is None
            or staged_name in staged_names
            or logical_name in logical_names
        ):
            raise _ResourceFailure
        seen.add(parts)
        staged_names.add(staged_name)
        logical_names.add(logical_name)
        entries.append((line, parts))
    return tuple(entries)


def _parse_index_binding(
    content: bytes,
    *,
    policy: ResourceClosurePolicy,
) -> dict[str, Any]:
    try:
        value = json.loads(content.decode("utf-8"), object_pairs_hook=_unique_object)
    except (UnicodeDecodeError, json.JSONDecodeError, _DuplicateJsonKey) as error:
        raise _ResourceFailure from error
    if not isinstance(value, dict) or set(value) != {
        "schema_version",
        "database_closure_sha256",
        "sortmerna_version",
        "files",
    }:
        raise _ResourceFailure
    if value["schema_version"] != SORTMERNA_INDEX_MANIFEST_SCHEMA_VERSION:
        raise _ResourceFailure
    database_identity = value["database_closure_sha256"]
    sortmerna_version = value["sortmerna_version"]
    _require_sha256(database_identity)
    if not isinstance(sortmerna_version, str) or not _SAFE_VERSION.fullmatch(
        sortmerna_version
    ):
        raise _ResourceFailure
    raw_files = value["files"]
    if (
        not isinstance(raw_files, list)
        or not raw_files
        or len(raw_files) > policy.maximum_index_entries
    ):
        raise _ResourceFailure
    files: list[dict[str, Any]] = []
    paths: list[str] = []
    total_bytes = 0
    for item in raw_files:
        if not isinstance(item, dict) or set(item) != {
            "path",
            "size_bytes",
            "sha256",
        }:
            raise _ResourceFailure
        path = item["path"]
        size_bytes = item["size_bytes"]
        file_sha256 = item["sha256"]
        if not isinstance(path, str):
            raise _ResourceFailure
        parts = _safe_relative_parts(
            path,
            maximum_bytes=policy.maximum_relative_path_bytes,
        )
        if (
            parts == (SORTMERNA_INDEX_MANIFEST_FILENAME,)
            or isinstance(size_bytes, bool)
            or not isinstance(size_bytes, int)
            or size_bytes < 0
            or size_bytes > policy.maximum_index_file_bytes
        ):
            raise _ResourceFailure
        _require_sha256(file_sha256)
        total_bytes += size_bytes
        if total_bytes > policy.maximum_index_total_bytes:
            raise _ResourceFailure
        paths.append(path)
        files.append({"path": path, "size_bytes": size_bytes, "sha256": file_sha256})
    if paths != sorted(paths) or len(paths) != len(set(paths)):
        raise _ResourceFailure
    declared_paths = set(paths)
    declared_directories: set[str] = set()
    for path in paths:
        parent = PurePosixPath(path).parent
        while parent.as_posix() != ".":
            declared_directories.add(parent.as_posix())
            parent = parent.parent
    if declared_paths.intersection(declared_directories):
        raise _ResourceFailure
    return {
        "database_closure_sha256": database_identity,
        "sortmerna_version": sortmerna_version,
        "files": files,
    }


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            raise _DuplicateJsonKey
        value[key] = item
    return value


def _database_closure_identity(
    manifest_sha256: str,
    files: list[RiboDatabaseFile],
) -> str:
    digest = hashlib.sha256()
    _update_framed(digest, RIBO_DATABASE_CLOSURE_SCHEME.encode())
    _update_framed(digest, bytes.fromhex(manifest_sha256))
    for resource in files:
        _update_framed(digest, resource.manifest_entry.encode("utf-8"))
        _update_framed(digest, str(resource.size_bytes).encode("ascii"))
        _update_framed(digest, bytes.fromhex(resource.sha256))
    return digest.hexdigest()


def _index_closure_identity(
    manifest_sha256: str,
    entries: tuple[SortMeRnaIndexEntry, ...],
) -> str:
    digest = hashlib.sha256()
    _update_framed(digest, SORTMERNA_INDEX_CLOSURE_SCHEME.encode())
    _update_framed(digest, bytes.fromhex(manifest_sha256))
    for entry in entries:
        _update_framed(digest, entry.relative_path.encode("utf-8"))
        assert entry.kind == "file"
        assert entry.size_bytes is not None
        assert entry.sha256 is not None
        _update_framed(digest, str(entry.size_bytes).encode("ascii"))
        _update_framed(digest, bytes.fromhex(entry.sha256))
    return digest.hexdigest()


def _update_framed(digest: Any, value: bytes) -> None:
    digest.update(len(value).to_bytes(8, "big"))
    digest.update(value)


def _require_sha256(value: object) -> str:
    if not isinstance(value, str) or _SHA256.fullmatch(value) is None:
        raise _ResourceFailure
    return value


def _require_policy(value: object) -> ResourceClosurePolicy:
    if not isinstance(value, ResourceClosurePolicy):
        raise _ResourceFailure
    return value


def _absolute_lexical_path(value: str | Path) -> Path:
    if not isinstance(value, (str, Path)) or isinstance(value, bool):
        raise _ResourceFailure
    raw = os.fspath(value)
    if not raw or "\x00" in raw:
        raise _ResourceFailure
    path = Path(raw)
    if not path.is_absolute() or any(
        part in {"", ".", ".."}
        or _SAFE_EXECUTION_PATH_COMPONENT.fullmatch(part) is None
        for part in path.parts[1:]
    ):
        raise _ResourceFailure
    return path


def _safe_relative_parts(value: str, *, maximum_bytes: int) -> tuple[str, ...]:
    if (
        not value
        or value.startswith(("/", "//"))
        or "\\" in value
        or "\x00" in value
        or _URI_SCHEME.match(value)
        or len(value.encode("utf-8")) > maximum_bytes
        or any(ord(character) < 32 or ord(character) == 127 for character in value)
    ):
        raise _ResourceFailure
    raw_parts = value.split("/")
    if any(
        part in {"", ".", ".."}
        or _SAFE_EXECUTION_PATH_COMPONENT.fullmatch(part) is None
        for part in raw_parts
    ):
        raise _ResourceFailure
    parts = PurePosixPath(value).parts
    if tuple(raw_parts) != parts:
        raise _ResourceFailure
    return parts


def _safe_leaf_name(value: object) -> str:
    if (
        not isinstance(value, str)
        or not value
        or value in {".", ".."}
        or "/" in value
        or "\x00" in value
        or _SAFE_EXECUTION_PATH_COMPONENT.fullmatch(value) is None
        or any(ord(character) < 32 or ord(character) == 127 for character in value)
    ):
        raise _ResourceFailure
    return value


def _open_absolute_directory(path: Path) -> int:
    if not path.is_absolute() or any(
        part in {"", ".", ".."} for part in path.parts[1:]
    ):
        raise _ResourceFailure
    current = os.open("/", _directory_flags())
    try:
        for component in path.parts[1:]:
            descriptor = os.open(component, _directory_flags(), dir_fd=current)
            os.close(current)
            current = descriptor
        info = os.fstat(current)
        if not stat.S_ISDIR(info.st_mode):
            raise _ResourceFailure
        return current
    except Exception:
        with suppress(OSError):
            os.close(current)
        raise


def _directory_flags() -> int:
    if not hasattr(os, "O_DIRECTORY") or not hasattr(os, "O_NOFOLLOW"):
        raise _ResourceFailure
    return (
        os.O_RDONLY
        | os.O_DIRECTORY
        | os.O_NOFOLLOW
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOCTTY", 0)
    )


def _file_flags() -> int:
    if not hasattr(os, "O_NOFOLLOW"):
        raise _ResourceFailure
    return (
        os.O_RDONLY
        | os.O_NOFOLLOW
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NONBLOCK", 0)
        | getattr(os, "O_NOCTTY", 0)
    )


def _file_identity(info: os.stat_result) -> tuple[int, ...]:
    return (
        info.st_dev,
        info.st_ino,
        info.st_mode,
        info.st_nlink,
        info.st_size,
        info.st_mtime_ns,
        info.st_ctime_ns,
    )


def _directory_identity(info: os.stat_result) -> tuple[int, ...]:
    return (
        info.st_dev,
        info.st_ino,
        info.st_mode,
        info.st_nlink,
        info.st_size,
        info.st_mtime_ns,
        info.st_ctime_ns,
    )


def _failure(code: str) -> Result[Any]:
    if code == "BULK_RNASEQ_LOCAL_FILE_INVALID":
        message = "The local input resource could not be verified."
        path = "runtime.resource"
    elif code == "BULK_RNASEQ_SORTMERNA_INDEX_INVALID":
        message = "The SortMeRNA index resource could not be verified."
        path = "config.standard.ribosomal_rna_removal.sortmerna_index"
    else:
        message = "The ribosomal RNA database resource could not be verified."
        path = "config.standard.ribosomal_rna_removal.database_manifest"
    return Result.failure(
        [
            Issue(
                code=code,
                message=message,
                source="adapter",
                path=path,
                context={},
            )
        ]
    )
