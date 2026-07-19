"""Installed implementation identity for the bulk RNA-seq execution boundary."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from importlib import resources
import json
from pathlib import Path, PurePosixPath
import re
from typing import Any

from encode_pipeline.adapters.bulk_rnaseq.resource_closure import (
    safe_regular_file_bytes,
)
from encode_pipeline.platform.results import Issue, Result


EXECUTION_IMPLEMENTATION_MANIFEST_FILE = "execution-implementation-manifest-1.0.0.json"
EXECUTION_IMPLEMENTATION_SCHEMA_VERSION = "1.0.0"
EXECUTION_IMPLEMENTATION_SCHEME = "sha256-framed-execution-implementation-v1"

EXECUTION_MIGRATION_REVISION_PATHS = (
    "src/encode_pipeline/persistence/alembic/versions/20260711_01_run_persistence.py",
    "src/encode_pipeline/persistence/alembic/versions/20260711_02_run_execution_assignments.py",
    "src/encode_pipeline/persistence/alembic/versions/20260712_03_run_workflow_build_identities.py",
    "src/encode_pipeline/persistence/alembic/versions/20260712_04_run_cancellation_intent.py",
    "src/encode_pipeline/persistence/alembic/versions/20260712_05_run_qc_metrics.py",
    "src/encode_pipeline/persistence/alembic/versions/20260714_06_validated_input_snapshots.py",
    "src/encode_pipeline/persistence/alembic/versions/20260714_07_run_history_indexes.py",
    "src/encode_pipeline/persistence/alembic/versions/20260717_08_run_result_generations.py",
)

# This is an exact allowlist, not a recursive glob. A new execution dependency
# is added deliberately and requires regenerating the committed manifest.
EXECUTION_IMPLEMENTATION_PATHS = (
    "src/encode_pipeline/adapters/bulk_rnaseq/__init__.py",
    "src/encode_pipeline/adapters/bulk_rnaseq/adapter.py",
    "src/encode_pipeline/adapters/bulk_rnaseq/artifacts.py",
    "src/encode_pipeline/adapters/bulk_rnaseq/authoring.py",
    "src/encode_pipeline/adapters/bulk_rnaseq/execution.py",
    "src/encode_pipeline/adapters/bulk_rnaseq/execution_identity.py",
    "src/encode_pipeline/adapters/bulk_rnaseq/reference_closure.py",
    "src/encode_pipeline/adapters/bulk_rnaseq/resource_closure.py",
    "src/encode_pipeline/adapters/bulk_rnaseq/results_contract.py",
    "src/encode_pipeline/adapters/bulk_rnaseq/runtime_assets.py",
    "src/encode_pipeline/adapters/bulk_rnaseq/status_evidence.py",
    "src/encode_pipeline/adapters/bulk_rnaseq/qc.py",
    "src/encode_pipeline/adapters/bulk_rnaseq/qualification.py",
    "src/encode_pipeline/adapters/bulk_rnaseq/upstream.py",
    "src/encode_pipeline/adapters/bulk_rnaseq/validation.py",
    "src/encode_pipeline/platform/__init__.py",
    "src/encode_pipeline/platform/adapters.py",
    "src/encode_pipeline/platform/builds.py",
    "src/encode_pipeline/platform/execution.py",
    "src/encode_pipeline/platform/managed_containers.py",
    "src/encode_pipeline/platform/planning.py",
    "src/encode_pipeline/platform/registry.py",
    "src/encode_pipeline/platform/result_generations.py",
    "src/encode_pipeline/platform/results.py",
    "src/encode_pipeline/platform/runs.py",
    "src/encode_pipeline/contracts/nfcore_rnaseq/results-contract-3.26.0.json",
    "src/encode_pipeline/persistence/runtime.py",
    "src/encode_pipeline/persistence/database.py",
    "src/encode_pipeline/persistence/migrations.py",
    "src/encode_pipeline/persistence/models.py",
    "src/encode_pipeline/persistence/repositories.py",
    "src/encode_pipeline/persistence/alembic/env.py",
    *EXECUTION_MIGRATION_REVISION_PATHS,
    "src/encode_pipeline/services/artifact_extraction.py",
    "src/encode_pipeline/services/artifact_downloads.py",
    "src/encode_pipeline/services/command_builder.py",
    "src/encode_pipeline/services/defaults.py",
    "src/encode_pipeline/services/local_execution.py",
    "src/encode_pipeline/services/local_run_driver.py",
    "src/encode_pipeline/services/managed_containers.py",
    "src/encode_pipeline/services/materialization.py",
    "src/encode_pipeline/services/planning.py",
    "src/encode_pipeline/services/preflight.py",
    "src/encode_pipeline/services/process_runner.py",
    "src/encode_pipeline/services/qc_summary_indexing.py",
    "src/encode_pipeline/services/run_repositories.py",
    "src/encode_pipeline/services/runs.py",
    "src/encode_pipeline/services/workflow_builds.py",
    "src/encode_pipeline/workers/jobs.py",
    "src/encode_pipeline/workers/rq_queue.py",
    "src/encode_pipeline/workers/runtime.py",
    "src/encode_pipeline/workers/settings.py",
    "src/encode_pipeline/workers/timeouts.py",
)

_CONTRACT_PACKAGE = "encode_pipeline.contracts.nfcore_rnaseq"
_PACKAGE_PREFIX = PurePosixPath("src/encode_pipeline")
_MAXIMUM_MANIFEST_BYTES = 512 * 1024
_MAXIMUM_IMPLEMENTATION_FILE_BYTES = 8 * 1024 * 1024
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


@dataclass(frozen=True)
class ExecutionImplementationFile:
    """One exact installed source file bound by the implementation manifest."""

    path: str
    size_bytes: int
    sha256: str


@dataclass(frozen=True)
class VerifiedExecutionImplementation:
    """Verified identity of the complete controlled execution implementation."""

    manifest_sha256: str
    aggregate_sha256: str
    files: tuple[ExecutionImplementationFile, ...]


class _ImplementationFailure(Exception):
    pass


class _DuplicateJsonKey(ValueError):
    pass


def verify_execution_implementation(
    *,
    manifest_bytes: bytes | None = None,
    package_root: object | None = None,
) -> Result[VerifiedExecutionImplementation]:
    """Verify the committed manifest and every installed controlled file."""

    try:
        content = (
            _read_contract_manifest()
            if manifest_bytes is None
            else _bounded_bytes(manifest_bytes, _MAXIMUM_MANIFEST_BYTES)
        )
        manifest_sha256 = hashlib.sha256(content).hexdigest()
        manifest = _parse_manifest(content)
        root = (
            resources.files("encode_pipeline") if package_root is None else package_root
        )
        _verify_migration_revision_set(root)
        files: list[ExecutionImplementationFile] = []
        for item in manifest["files"]:
            observed = _read_package_file(root, item["path"])
            if (
                len(observed) != item["size_bytes"]
                or hashlib.sha256(observed).hexdigest() != item["sha256"]
            ):
                raise _ImplementationFailure
            files.append(
                ExecutionImplementationFile(
                    path=item["path"],
                    size_bytes=item["size_bytes"],
                    sha256=item["sha256"],
                )
            )
        return Result.success(
            VerifiedExecutionImplementation(
                manifest_sha256=manifest_sha256,
                aggregate_sha256=manifest["aggregate_sha256"],
                files=tuple(files),
            )
        )
    except (
        _ImplementationFailure,
        OSError,
        TypeError,
        ValueError,
        UnicodeError,
    ):
        return _failure()


def build_execution_implementation_manifest(project_root: Path) -> dict[str, Any]:
    """Build canonical manifest data from one checked-out repository root."""

    if not isinstance(project_root, Path) or not project_root.is_absolute():
        raise ValueError("project_root must be an absolute Path")
    try:
        _verify_migration_revision_set(project_root / "src/encode_pipeline")
    except _ImplementationFailure as error:
        raise ValueError("production migration revision set is invalid") from error
    files: list[dict[str, Any]] = []
    for logical_path in EXECUTION_IMPLEMENTATION_PATHS:
        result = safe_regular_file_bytes(
            project_root.joinpath(*PurePosixPath(logical_path).parts),
            maximum_bytes=_MAXIMUM_IMPLEMENTATION_FILE_BYTES,
        )
        if result.is_failure:
            raise ValueError("controlled implementation file is unavailable")
        files.append(
            {
                "path": logical_path,
                "size_bytes": result.value.size_bytes,
                "sha256": result.value.sha256,
            }
        )
    aggregate = execution_implementation_aggregate(files)
    return {
        "schema_version": EXECUTION_IMPLEMENTATION_SCHEMA_VERSION,
        "scheme": EXECUTION_IMPLEMENTATION_SCHEME,
        "file_count": len(files),
        "aggregate_sha256": aggregate,
        "files": files,
    }


def canonical_execution_manifest_bytes(manifest: object) -> bytes:
    """Serialize generated manifest data deterministically."""

    return (
        json.dumps(manifest, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        + "\n"
    ).encode()


def execution_implementation_aggregate(files: object) -> str:
    """Return the framed aggregate for already-normalized file records."""

    if not isinstance(files, list):
        raise ValueError("files must be a list")
    digest = hashlib.sha256()
    _frame(digest, EXECUTION_IMPLEMENTATION_SCHEME.encode())
    for item in files:
        if not isinstance(item, dict):
            raise ValueError("file entry must be an object")
        path = item.get("path")
        size_bytes = item.get("size_bytes")
        file_sha256 = item.get("sha256")
        if (
            not isinstance(path, str)
            or isinstance(size_bytes, bool)
            or not isinstance(size_bytes, int)
            or size_bytes <= 0
            or not isinstance(file_sha256, str)
            or _SHA256.fullmatch(file_sha256) is None
        ):
            raise ValueError("file entry is invalid")
        _frame(digest, path.encode())
        _frame(digest, str(size_bytes).encode())
        _frame(digest, bytes.fromhex(file_sha256))
    return digest.hexdigest()


def _read_contract_manifest() -> bytes:
    resource = resources.files(_CONTRACT_PACKAGE).joinpath(
        EXECUTION_IMPLEMENTATION_MANIFEST_FILE
    )
    return _read_traversable(resource, maximum_bytes=_MAXIMUM_MANIFEST_BYTES)


def _verify_migration_revision_set(root: object) -> None:
    resource = root
    for part in ("persistence", "alembic", "versions"):
        if isinstance(resource, Path):
            resource = resource / part
        else:
            joinpath = getattr(resource, "joinpath", None)
            if not callable(joinpath):
                raise _ImplementationFailure
            resource = joinpath(part)

    expected = {PurePosixPath(path).name for path in EXECUTION_MIGRATION_REVISION_PATHS}
    if isinstance(resource, Path):
        if resource.is_symlink() or not resource.is_dir():
            raise _ImplementationFailure
        entries = tuple(resource.iterdir())
    else:
        is_dir = getattr(resource, "is_dir", None)
        iterdir = getattr(resource, "iterdir", None)
        if not callable(is_dir) or not is_dir() or not callable(iterdir):
            raise _ImplementationFailure
        entries = tuple(iterdir())

    observed: set[str] = set()
    for entry in entries:
        name = getattr(entry, "name", None)
        if (
            not isinstance(name, str)
            or not name.endswith(".py")
            or name == "__init__.py"
        ):
            continue
        is_file = getattr(entry, "is_file", None)
        if not callable(is_file) or not is_file():
            raise _ImplementationFailure
        if isinstance(entry, Path) and entry.is_symlink():
            raise _ImplementationFailure
        observed.add(name)
    if observed != expected:
        raise _ImplementationFailure


def _read_package_file(root: object, logical_path: str) -> bytes:
    path = PurePosixPath(logical_path)
    try:
        relative = path.relative_to(_PACKAGE_PREFIX)
    except ValueError as error:
        raise _ImplementationFailure from error
    if isinstance(root, Path):
        result = safe_regular_file_bytes(
            root.joinpath(*relative.parts),
            maximum_bytes=_MAXIMUM_IMPLEMENTATION_FILE_BYTES,
        )
        if result.is_failure:
            raise _ImplementationFailure
        return result.value.content
    resource = root
    for part in relative.parts:
        joinpath = getattr(resource, "joinpath", None)
        if not callable(joinpath):
            raise _ImplementationFailure
        resource = joinpath(part)
    return _read_traversable(
        resource,
        maximum_bytes=_MAXIMUM_IMPLEMENTATION_FILE_BYTES,
    )


def _read_traversable(resource: object, *, maximum_bytes: int) -> bytes:
    if isinstance(resource, Path):
        result = safe_regular_file_bytes(resource, maximum_bytes=maximum_bytes)
        if result.is_failure:
            raise _ImplementationFailure
        return result.value.content
    is_file = getattr(resource, "is_file", None)
    read_bytes = getattr(resource, "read_bytes", None)
    if not callable(is_file) or not is_file() or not callable(read_bytes):
        raise _ImplementationFailure
    return _bounded_bytes(read_bytes(), maximum_bytes)


def _bounded_bytes(value: object, maximum_bytes: int) -> bytes:
    if not isinstance(value, bytes) or not value or len(value) > maximum_bytes:
        raise _ImplementationFailure
    return value


def _parse_manifest(content: bytes) -> dict[str, Any]:
    try:
        value = json.loads(content.decode("utf-8"), object_pairs_hook=_unique_object)
    except (UnicodeDecodeError, json.JSONDecodeError, _DuplicateJsonKey) as error:
        raise _ImplementationFailure from error
    if not isinstance(value, dict) or set(value) != {
        "schema_version",
        "scheme",
        "file_count",
        "aggregate_sha256",
        "files",
    }:
        raise _ImplementationFailure
    if (
        value["schema_version"] != EXECUTION_IMPLEMENTATION_SCHEMA_VERSION
        or value["scheme"] != EXECUTION_IMPLEMENTATION_SCHEME
        or isinstance(value["file_count"], bool)
        or value["file_count"] != len(EXECUTION_IMPLEMENTATION_PATHS)
        or not isinstance(value["files"], list)
        or len(value["files"]) != len(EXECUTION_IMPLEMENTATION_PATHS)
    ):
        raise _ImplementationFailure
    normalized: list[dict[str, Any]] = []
    for item in value["files"]:
        if not isinstance(item, dict) or set(item) != {
            "path",
            "size_bytes",
            "sha256",
        }:
            raise _ImplementationFailure
        path = item["path"]
        size_bytes = item["size_bytes"]
        file_sha256 = item["sha256"]
        if (
            not isinstance(path, str)
            or isinstance(size_bytes, bool)
            or not isinstance(size_bytes, int)
            or size_bytes <= 0
            or size_bytes > _MAXIMUM_IMPLEMENTATION_FILE_BYTES
            or not isinstance(file_sha256, str)
            or _SHA256.fullmatch(file_sha256) is None
        ):
            raise _ImplementationFailure
        normalized.append(
            {"path": path, "size_bytes": size_bytes, "sha256": file_sha256}
        )
    paths = tuple(item["path"] for item in normalized)
    if paths != EXECUTION_IMPLEMENTATION_PATHS:
        raise _ImplementationFailure
    aggregate = execution_implementation_aggregate(normalized)
    if value["aggregate_sha256"] != aggregate:
        raise _ImplementationFailure
    return {**value, "files": normalized}


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise _DuplicateJsonKey
        result[key] = value
    return result


def _frame(digest: Any, value: bytes) -> None:
    digest.update(len(value).to_bytes(8, "big"))
    digest.update(value)


def _failure() -> Result[VerifiedExecutionImplementation]:
    return Result.failure(
        [
            Issue(
                code="BULK_RNASEQ_EXECUTION_IMPLEMENTATION_INVALID",
                message="The bulk RNA-seq execution implementation could not be verified.",
                severity="error",
                path="runtime.implementation",
                source="adapter",
                context={},
            )
        ]
    )
