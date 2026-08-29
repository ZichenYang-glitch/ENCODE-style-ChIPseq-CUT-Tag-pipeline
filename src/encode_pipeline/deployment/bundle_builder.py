"""Deterministic offline bundle producer for the one supported topology.

Inputs are already-built wheel or scientific runtime bytes.  The producer does
not invent a release catalog: contract bindings are identities re-derived from
the same native bytes that the production resolver later verifies.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import io
from importlib import resources
import os
from pathlib import Path
import shutil
import stat
import tarfile
import tempfile
import uuid

from encode_pipeline.adapters import EncodeStyleWorkflowAdapter
from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import (
    RUNTIME_IDENTITY_FILE,
    RUNTIME_IDENTITY_SHA256,
    RuntimeAssetBinding,
    VerifiedRuntimeAssets,
    verify_packaged_runtime_asset_closure,
)
from encode_pipeline.deployment.canonical import canonical_json_bytes
from encode_pipeline.deployment.bundle import (
    MAX_BUNDLE_BYTES,
    MAX_BUNDLE_FILES,
    MAX_MANIFEST_BYTES,
    MIN_FREE_SPACE_RESERVE,
)
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.models import (
    BULK_RNASEQ_RUNTIME,
    ENCODE_RUNTIME,
    PLATFORM,
    BundleManifest,
    ContractDocument,
    FileRecord,
)
from encode_pipeline.deployment.native_contracts import (
    BULK_RUNTIME_CONTRACT,
    BULK_RUNTIME_IDENTITY_PATH,
    BULK_RUNTIME_ROOT_PATH,
    ENCODE_RUNTIME_CONTRACT,
    ENCODE_MICROMAMBA_PATH,
    ENCODE_PACKAGE_ARCHIVE_ROOT,
    ENCODE_RUNTIME_INDEX_PATH,
    ENCODE_RUNTIME_ROOT_PATH,
    ExplicitPackageCoordinate,
    PLATFORM_BULK_REFERENCES_PATH,
    PLATFORM_DISTRIBUTION_CONTRACT,
    PLATFORM_ENCODE_REFERENCES_PATH,
    PLATFORM_FRONTEND_CONTRACT,
    PLATFORM_FRONTEND_PATH,
    PLATFORM_METADATA_PATH,
    PLATFORM_MIGRATIONS_CONTRACT,
    PLATFORM_MIGRATIONS_PATH,
    PLATFORM_PYTHON_RUNTIME_CONTRACT,
    PLATFORM_REFERENCES_CONTRACT,
    PLATFORM_WHEEL_CONTRACT,
    PLATFORM_WHEEL_PATH,
    encode_environment_paths,
    encode_runtime_index_bytes,
    parse_encode_explicit_lock,
    parse_encode_runtime_index,
    verify_static_micromamba,
    verify_platform_wheel,
)
from encode_pipeline.deployment.platform_runtime import (
    PLATFORM_RUNTIME_LOCK_PATH,
    PLATFORM_RUNTIME_WHEELHOUSE_ROOT,
    PlatformRuntimeClosure,
    PlatformWheelLock,
    collect_platform_runtime_closure,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider


_MAX_PLATFORM_WHEEL_BYTES = 2 * 1024**3
_MAX_CONTRACT_BYTES = 2 * 1024 * 1024
_MAX_PAYLOAD_FILE_BYTES = 100 * 1024**3
_BULK_RUNTIME_DIRECTORIES = frozenset(
    {"source", "nextflow", "jdk", "plugins", "containers"}
)


@dataclass(frozen=True)
class _BytesPayload:
    logical_path: str
    content: bytes
    mode: int = 0o444


@dataclass(frozen=True)
class _PathPayload:
    logical_path: str
    source: Path
    mode: int
    size_bytes: int
    sha256: str
    witness: tuple[int, ...]


Payload = _BytesPayload | _PathPayload


def build_platform_bundle(
    wheel_path: Path,
    wheelhouse: Path,
    wheel_lock_path: Path,
    output_path: Path,
    *,
    scratch_root: Path,
) -> BundleManifest:
    """Build a canonical platform tar with its expanded, locked Python closure."""
    wheel = _read_source_bytes(wheel_path, maximum_bytes=_MAX_PLATFORM_WHEEL_BYTES)
    try:
        facts = verify_platform_wheel(wheel)
        lock_content = _read_source_bytes(
            wheel_lock_path,
            maximum_bytes=_MAX_CONTRACT_BYTES,
        )
        lock = PlatformWheelLock.from_bytes(lock_content)
    except Exception:
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=PLATFORM,
        ) from None
    wheel_sha256 = hashlib.sha256(wheel).hexdigest()
    matching_wheels = tuple(item for item in lock.wheels if item.sha256 == wheel_sha256)
    if len(matching_wheels) != 1:
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=PLATFORM,
        )
    scratch = _scratch_directory(scratch_root)
    try:
        with tempfile.TemporaryDirectory(
            prefix=".helixweave-platform-",
            dir=scratch,
            ignore_cleanup_errors=True,
        ) as temporary:
            closure_root = Path(temporary) / "closure"
            closure = collect_platform_runtime_closure(
                wheelhouse,
                wheel_lock_path,
                closure_root,
            )
            runtime_payload = _platform_runtime_payload(
                closure,
                closure_root=closure_root,
                lock=lock,
                lock_content=lock_content,
                platform_wheel_filename=matching_wheels[0].filename,
                platform_wheel_sha256=wheel_sha256,
            )
            payload: tuple[Payload, ...] = (
                _BytesPayload(PLATFORM_METADATA_PATH, facts.metadata),
                _BytesPayload(PLATFORM_WHEEL_PATH, wheel),
                _BytesPayload(PLATFORM_FRONTEND_PATH, facts.frontend_manifest),
                _BytesPayload(PLATFORM_MIGRATIONS_PATH, facts.migration_inventory),
                _BytesPayload(
                    PLATFORM_ENCODE_REFERENCES_PATH,
                    facts.encode_reference_source,
                ),
                _BytesPayload(
                    PLATFORM_BULK_REFERENCES_PATH,
                    facts.bulk_reference_source,
                ),
                *runtime_payload,
            )
            contracts = (
                ContractDocument(
                    PLATFORM_DISTRIBUTION_CONTRACT,
                    _raw_identity(facts.metadata),
                    PLATFORM_METADATA_PATH,
                ),
                ContractDocument(
                    PLATFORM_WHEEL_CONTRACT,
                    _raw_identity(wheel),
                    PLATFORM_WHEEL_PATH,
                ),
                ContractDocument(
                    PLATFORM_PYTHON_RUNTIME_CONTRACT,
                    closure.lock_identity,
                    PLATFORM_RUNTIME_LOCK_PATH,
                ),
                ContractDocument(
                    PLATFORM_FRONTEND_CONTRACT,
                    facts.frontend_identity,
                    PLATFORM_FRONTEND_PATH,
                ),
                ContractDocument(
                    PLATFORM_MIGRATIONS_CONTRACT,
                    _raw_identity(facts.migration_inventory),
                    PLATFORM_MIGRATIONS_PATH,
                ),
                ContractDocument(
                    PLATFORM_REFERENCES_CONTRACT,
                    facts.reference_identity,
                    PLATFORM_ENCODE_REFERENCES_PATH,
                ),
            )
            return _build_bundle(
                component=PLATFORM,
                contracts=contracts,
                payload=payload,
                output_path=output_path,
            )
    except DeploymentError:
        raise
    except Exception:
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=PLATFORM,
        ) from None


def build_encode_runtime_bundle(
    project_root: Path,
    micromamba_path: Path,
    archive_cache: Path,
    output_path: Path,
) -> BundleManifest:
    """Build the ENCODE source plus its fully offline native tool closure."""
    root = _source_directory(project_root)
    cache = _source_directory(archive_cache)
    try:
        adapter = EncodeStyleWorkflowAdapter(
            catalog_path=str(root / "docs" / "architecture" / "artifact-inventory.yaml")
        )
        registry = WorkflowRegistry(
            adapters=(adapter,),
            legacy_execution_fallbacks=(adapter,),
        )
        provider = WorkflowBuildIdentityProvider(registry, project_root=root)
        captured = provider.capture(adapter.metadata.workflow_id)
        source_manifest = provider.source_manifest()
    except Exception:
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=ENCODE_RUNTIME,
        ) from None
    if captured.is_failure:
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=ENCODE_RUNTIME,
        )
    workflow_identity = f"sha256-{captured.value.digest}"
    payload: list[Payload] = []
    for logical_path, content in source_manifest:
        source = root.joinpath(*Path(logical_path).parts)
        try:
            observed = source.lstat()
        except OSError:
            raise fail(
                "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
                "Deployment bundle source is invalid.",
                component=ENCODE_RUNTIME,
            ) from None
        if (
            not stat.S_ISREG(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or observed.st_nlink != 1
            or observed.st_mode & 0o022
        ):
            raise fail(
                "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
                "Deployment bundle source is invalid.",
                component=ENCODE_RUNTIME,
            )
        if _read_source_bytes(source, maximum_bytes=_MAX_PAYLOAD_FILE_BYTES) != content:
            raise fail(
                "DEPLOYMENT_BUNDLE_SOURCE_CHANGED",
                "Deployment bundle source changed during verification.",
                component=ENCODE_RUNTIME,
                recoverable=True,
            )
        payload.append(
            _BytesPayload(
                f"{ENCODE_RUNTIME_ROOT_PATH}/{logical_path}",
                content,
                0o555 if stat.S_IMODE(observed.st_mode) & 0o111 else 0o444,
            )
        )
    try:
        environment_paths = encode_environment_paths(source_manifest)
        source = dict(source_manifest)
        environment_packages: dict[str, tuple[ExplicitPackageCoordinate, ...]] = {}
        coordinates: dict[str, ExplicitPackageCoordinate] = {}
        basenames: dict[str, str] = {}
        for environment_path in environment_paths:
            lock_path = environment_path[:-4] + ".lock"
            packages = parse_encode_explicit_lock(source[lock_path])
            environment_packages[environment_path] = packages
            for package in packages:
                prior_url = basenames.get(package.filename)
                if prior_url is not None and prior_url != package.url:
                    raise ValueError
                basenames[package.filename] = package.url
                coordinates[package.url] = package
        cache_entries = tuple(sorted(cache.iterdir()))
        if {item.name for item in cache_entries} != set(basenames) or len(
            cache_entries
        ) != len(basenames):
            raise ValueError
        by_name = {item.name: item for item in cache_entries}
        package_documents: list[dict[str, object]] = []
        for url in sorted(coordinates):
            package = coordinates[url]
            archive_payload = _path_payload(
                by_name[package.filename],
                "unused",
                component=ENCODE_RUNTIME,
            )
            _verify_path_md5(archive_payload, package.md5)
            logical_path = (
                f"{ENCODE_PACKAGE_ARCHIVE_ROOT}/{archive_payload.sha256}/"
                f"{package.filename}"
            )
            archive_payload = _PathPayload(
                logical_path=logical_path,
                source=archive_payload.source,
                mode=0o444,
                size_bytes=archive_payload.size_bytes,
                sha256=archive_payload.sha256,
                witness=archive_payload.witness,
            )
            payload.append(archive_payload)
            package_documents.append(
                {
                    "url": package.url,
                    "platform": package.platform,
                    "filename": package.filename,
                    "md5": package.md5,
                    "archive_path": logical_path,
                    "size_bytes": archive_payload.size_bytes,
                    "sha256": archive_payload.sha256,
                }
            )
        micromamba = _path_payload(
            micromamba_path,
            ENCODE_MICROMAMBA_PATH,
            component=ENCODE_RUNTIME,
        )
        if micromamba.mode != 0o555:
            raise ValueError
        micromamba_content = _read_source_bytes(
            micromamba_path,
            maximum_bytes=_MAX_PLATFORM_WHEEL_BYTES,
        )
        if hashlib.sha256(micromamba_content).hexdigest() != micromamba.sha256:
            raise ValueError
        verify_static_micromamba(micromamba_content)
        payload.append(micromamba)
        environments = tuple(
            {
                "environment_path": environment_path,
                "environment_sha256": hashlib.sha256(
                    source[environment_path]
                ).hexdigest(),
                "lock_path": environment_path[:-4] + ".lock",
                "lock_sha256": hashlib.sha256(
                    source[environment_path[:-4] + ".lock"]
                ).hexdigest(),
                "packages": sorted(
                    package.url for package in environment_packages[environment_path]
                ),
            }
            for environment_path in environment_paths
        )
        index_content = encode_runtime_index_bytes(
            workflow_build_identity=workflow_identity,
            micromamba={
                "path": ENCODE_MICROMAMBA_PATH,
                "size_bytes": micromamba.size_bytes,
                "sha256": micromamba.sha256,
            },
            packages=tuple(package_documents),
            environments=environments,
        )
        index = parse_encode_runtime_index(index_content)
    except DeploymentError:
        raise
    except Exception:
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=ENCODE_RUNTIME,
        ) from None
    payload.append(_BytesPayload(ENCODE_RUNTIME_INDEX_PATH, index_content))
    contract = ContractDocument(
        ENCODE_RUNTIME_CONTRACT,
        index.identity,
        ENCODE_RUNTIME_INDEX_PATH,
    )
    return _build_bundle(
        component=ENCODE_RUNTIME,
        contracts=(contract,),
        payload=tuple(payload),
        output_path=output_path,
    )


def _platform_runtime_payload(
    closure: PlatformRuntimeClosure,
    *,
    closure_root: Path,
    lock: PlatformWheelLock,
    lock_content: bytes,
    platform_wheel_filename: str,
    platform_wheel_sha256: str,
) -> tuple[_PathPayload, ...]:
    if closure.lock_identity != lock.identity or not closure.files:
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=PLATFORM,
        )
    payload: list[_PathPayload] = []
    for item in closure.files:
        try:
            item.source.relative_to(closure_root)
        except (AttributeError, TypeError, ValueError):
            raise fail(
                "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
                "Deployment bundle source is invalid.",
                component=PLATFORM,
            ) from None
        captured = _path_payload(
            item.source,
            item.logical_path,
            component=PLATFORM,
        )
        if (
            captured.mode != item.mode
            or captured.size_bytes != item.size_bytes
            or captured.sha256 != item.sha256
        ):
            raise fail(
                "DEPLOYMENT_BUNDLE_SOURCE_CHANGED",
                "Deployment bundle source changed during verification.",
                component=PLATFORM,
                recoverable=True,
            )
        payload.append(captured)
    by_path = {item.logical_path: item for item in payload}
    if len(by_path) != len(payload):
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=PLATFORM,
        )
    lock_payload = by_path.get(PLATFORM_RUNTIME_LOCK_PATH)
    wheel_payload = by_path.get(
        f"{PLATFORM_RUNTIME_WHEELHOUSE_ROOT}/{platform_wheel_filename}"
    )
    if (
        lock_payload is None
        or wheel_payload is None
        or lock_payload.sha256 != hashlib.sha256(lock_content).hexdigest()
        or lock_payload.size_bytes != len(lock_content)
        or wheel_payload.sha256 != platform_wheel_sha256
    ):
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=PLATFORM,
        )
    return tuple(payload)


def build_bulk_rnaseq_runtime_bundle(
    runtime_root: Path,
    output_path: Path,
    *,
    docker_executable: Path = Path("/usr/bin/docker"),
    docker_socket: Path = Path("/run/helixweave/docker/docker.sock"),
    runtime_verifier=verify_packaged_runtime_asset_closure,
) -> BundleManifest:
    """Build the statically verified fixed nf-core/Nextflow runtime closure."""
    root = _source_directory(runtime_root)
    try:
        top_level = {path.name for path in root.iterdir()}
    except OSError:
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=BULK_RNASEQ_RUNTIME,
        ) from None
    if top_level != _BULK_RUNTIME_DIRECTORIES:
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=BULK_RNASEQ_RUNTIME,
        )
    try:
        binding = RuntimeAssetBinding(
            root=root,
            docker_executable=docker_executable,
            docker_socket=docker_socket,
        )
        result = runtime_verifier(binding)
    except Exception:
        result = None
    if (
        result is None
        or not hasattr(result, "is_success")
        or not result.is_success
        or not isinstance(result.value, VerifiedRuntimeAssets)
        or result.value.runtime_identity_sha256 != RUNTIME_IDENTITY_SHA256
    ):
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=BULK_RNASEQ_RUNTIME,
        )
    payload: list[Payload] = []
    for directory in sorted(_BULK_RUNTIME_DIRECTORIES):
        source_root = root / directory
        if not source_root.is_dir() or source_root.is_symlink():
            raise fail(
                "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
                "Deployment bundle source is invalid.",
                component=BULK_RNASEQ_RUNTIME,
            )
        for source in sorted(source_root.rglob("*")):
            if source.is_dir() and not source.is_symlink():
                continue
            relative = source.relative_to(root).as_posix()
            payload.append(
                _path_payload(
                    source,
                    f"{BULK_RUNTIME_ROOT_PATH}/{relative}",
                    component=BULK_RNASEQ_RUNTIME,
                )
            )
    try:
        runtime_identity = (
            resources.files("encode_pipeline.contracts.nfcore_rnaseq")
            .joinpath(RUNTIME_IDENTITY_FILE)
            .read_bytes()
        )
    except Exception:
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=BULK_RNASEQ_RUNTIME,
        ) from None
    if hashlib.sha256(runtime_identity).hexdigest() != RUNTIME_IDENTITY_SHA256:
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=BULK_RNASEQ_RUNTIME,
        )
    payload.append(_BytesPayload(BULK_RUNTIME_IDENTITY_PATH, runtime_identity))
    contract = ContractDocument(
        BULK_RUNTIME_CONTRACT,
        _raw_identity(runtime_identity),
        BULK_RUNTIME_IDENTITY_PATH,
    )
    return _build_bundle(
        component=BULK_RNASEQ_RUNTIME,
        contracts=(contract,),
        payload=tuple(payload),
        output_path=output_path,
    )


def _build_bundle(
    *,
    component: str,
    contracts: tuple[ContractDocument, ...],
    payload: tuple[Payload, ...],
    output_path: Path,
) -> BundleManifest:
    if not payload or len({item.logical_path for item in payload}) != len(payload):
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=component,
        )
    records = tuple(sorted(_record(item) for item in payload))
    manifest = BundleManifest.create(
        component=component,
        contracts=contracts,
        files=records,
    )
    _preflight_bundle_output(output_path, manifest)
    entries = {item.logical_path: item for item in payload}
    _write_canonical_tar(output_path, manifest, entries)
    return manifest


def _preflight_bundle_output(output_path: Path, manifest: BundleManifest) -> None:
    output = _output_path(output_path)
    manifest_content = canonical_json_bytes(manifest.to_dict())
    declared = sum(item.size_bytes for item in manifest.files) + MAX_MANIFEST_BYTES
    try:
        tar_bytes = _canonical_tar_size(manifest_content, manifest.files)
    except (OverflowError, ValueError):
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=manifest.component,
        ) from None
    if (
        len(manifest_content) > MAX_MANIFEST_BYTES
        or len(manifest.files) > MAX_BUNDLE_FILES
        or declared > MAX_BUNDLE_BYTES
        or tar_bytes > MAX_BUNDLE_BYTES
    ):
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=manifest.component,
        )
    try:
        available = shutil.disk_usage(output.parent).free
    except OSError:
        raise fail(
            "DEPLOYMENT_STORAGE_UNAVAILABLE",
            "Deployment storage is unavailable.",
            component=manifest.component,
            recoverable=True,
        ) from None
    reserve = max(MIN_FREE_SPACE_RESERVE, declared // 20)
    if available < max(declared, tar_bytes) + reserve:
        raise fail(
            "DEPLOYMENT_CAPACITY_INSUFFICIENT",
            "Deployment storage capacity is insufficient.",
            component=manifest.component,
            recoverable=True,
        )


def _canonical_tar_size(manifest_content: bytes, files: tuple[FileRecord, ...]) -> int:
    entries = (("manifest.json", len(manifest_content), 0o444),) + tuple(
        (item.path, item.size_bytes, item.mode) for item in files
    )
    offset = 0
    for name, size, mode in entries:
        header = _header(name, size, mode).tobuf(
            tarfile.GNU_FORMAT,
            tarfile.ENCODING,
            "surrogateescape",
        )
        offset += len(header) + (
            (size + tarfile.BLOCKSIZE - 1) // tarfile.BLOCKSIZE
        ) * (tarfile.BLOCKSIZE)
    return (
        (offset + 2 * tarfile.BLOCKSIZE + tarfile.RECORDSIZE - 1)
        // tarfile.RECORDSIZE
        * tarfile.RECORDSIZE
    )


def _write_canonical_tar(
    output_path: Path,
    manifest: BundleManifest,
    payload: dict[str, Payload],
) -> None:
    output = _output_path(output_path)
    temporary = f".{output.name}.{uuid.uuid4().hex}.partial"
    flags = (
        os.O_WRONLY
        | os.O_CREAT
        | os.O_EXCL
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    descriptor: int | None = None
    directory_descriptor: int | None = None
    linked = False
    published = False
    try:
        directory_descriptor = _open_output_directory(output.parent)
        descriptor = os.open(
            temporary,
            flags,
            0o400,
            dir_fd=directory_descriptor,
        )
        with os.fdopen(descriptor, "wb", closefd=False) as handle:
            with tarfile.open(
                fileobj=handle,
                mode="w",
                format=tarfile.GNU_FORMAT,
            ) as archive:
                manifest_content = canonical_json_bytes(manifest.to_dict())
                _add_bytes(archive, "manifest.json", manifest_content, 0o444)
                for record in manifest.files:
                    entry = payload[record.path]
                    if isinstance(entry, _BytesPayload):
                        _add_bytes(
                            archive,
                            record.path,
                            entry.content,
                            record.mode,
                        )
                    else:
                        _add_path(archive, entry, record)
            handle.flush()
            os.fchmod(descriptor, 0o440)
            os.fsync(descriptor)
        os.close(descriptor)
        descriptor = None
        try:
            os.link(
                temporary,
                output.name,
                src_dir_fd=directory_descriptor,
                dst_dir_fd=directory_descriptor,
                follow_symlinks=False,
            )
            linked = True
        except FileExistsError:
            raise fail(
                "DEPLOYMENT_BUNDLE_OUTPUT_EXISTS",
                "Deployment bundle output already exists.",
                component=manifest.component,
            ) from None
        os.unlink(temporary, dir_fd=directory_descriptor)
        os.fsync(directory_descriptor)
        published = True
    except DeploymentError:
        raise
    except (OSError, tarfile.TarError, ValueError):
        raise fail(
            "DEPLOYMENT_BUNDLE_BUILD_FAILED",
            "Deployment bundle could not be built.",
            component=manifest.component,
            recoverable=True,
        ) from None
    finally:
        if descriptor is not None:
            os.close(descriptor)
        if not published and directory_descriptor is not None:
            try:
                os.unlink(temporary, dir_fd=directory_descriptor)
            except OSError:
                pass
            if linked:
                try:
                    os.unlink(output.name, dir_fd=directory_descriptor)
                    os.fsync(directory_descriptor)
                except OSError:
                    pass
        if directory_descriptor is not None:
            os.close(directory_descriptor)


def _add_bytes(
    archive: tarfile.TarFile,
    name: str,
    content: bytes,
    mode: int,
) -> None:
    header = _header(name, len(content), mode)
    archive.addfile(header, io.BytesIO(content))


def _add_path(
    archive: tarfile.TarFile,
    entry: _PathPayload,
    record: FileRecord,
) -> None:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = -1
    try:
        descriptor = os.open(entry.source, flags)
        before = os.fstat(descriptor)
        if _file_witness(before) != entry.witness:
            raise OSError
        with os.fdopen(os.dup(descriptor), "rb") as handle:
            archive.addfile(
                _header(record.path, record.size_bytes, record.mode),
                handle,
            )
        after = os.fstat(descriptor)
        if _file_witness(after) != entry.witness:
            raise OSError
        if _hash_descriptor(descriptor) != record.sha256:
            raise OSError
    except OSError:
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_CHANGED",
            "Deployment bundle source changed during verification.",
            recoverable=True,
        ) from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _header(name: str, size: int, mode: int) -> tarfile.TarInfo:
    header = tarfile.TarInfo(name)
    header.type = tarfile.REGTYPE
    header.size = size
    header.mode = mode
    header.uid = 0
    header.gid = 0
    header.uname = ""
    header.gname = ""
    header.mtime = 0
    return header


def _record(entry: Payload) -> FileRecord:
    if isinstance(entry, _BytesPayload):
        return FileRecord(
            entry.logical_path,
            len(entry.content),
            hashlib.sha256(entry.content).hexdigest(),
            entry.mode,
        )
    return FileRecord(
        entry.logical_path,
        entry.size_bytes,
        entry.sha256,
        entry.mode,
    )


def _path_payload(source: Path, logical_path: str, *, component: str) -> _PathPayload:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = -1
    try:
        before = source.lstat()
        descriptor = os.open(source, flags)
        observed = os.fstat(descriptor)
        if (
            not stat.S_ISREG(observed.st_mode)
            or observed.st_nlink != 1
            or observed.st_mode & 0o022
            or not 0 <= observed.st_size <= _MAX_PAYLOAD_FILE_BYTES
            or _file_witness(before) != _file_witness(observed)
        ):
            raise OSError
        digest = _hash_descriptor(descriptor)
        after_descriptor = os.fstat(descriptor)
        after_path = source.lstat()
        if _file_witness(observed) != _file_witness(after_descriptor) or _file_witness(
            observed
        ) != _file_witness(after_path):
            raise OSError
    except OSError:
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=component,
        ) from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)
    return _PathPayload(
        logical_path=logical_path,
        source=source,
        mode=0o555 if stat.S_IMODE(observed.st_mode) & 0o111 else 0o444,
        size_bytes=observed.st_size,
        sha256=digest,
        witness=_file_witness(observed),
    )


def _read_source_bytes(path: Path, *, maximum_bytes: int) -> bytes:
    if not isinstance(path, Path) or not path.is_absolute():
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
        )
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = -1
    try:
        before = path.lstat()
        descriptor = os.open(path, flags)
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISREG(opened.st_mode)
            or opened.st_nlink != 1
            or opened.st_mode & 0o022
            or not 0 < opened.st_size <= maximum_bytes
            or _file_witness(before) != _file_witness(opened)
        ):
            raise OSError
        content = bytearray()
        while True:
            chunk = os.read(
                descriptor, min(1024 * 1024, maximum_bytes + 1 - len(content))
            )
            if not chunk:
                break
            content.extend(chunk)
            if len(content) > maximum_bytes:
                raise OSError
        after = os.fstat(descriptor)
        if len(content) != opened.st_size or _file_witness(opened) != _file_witness(
            after
        ):
            raise OSError
        return bytes(content)
    except OSError:
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
        ) from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _source_directory(path: Path) -> Path:
    if not isinstance(path, Path) or not path.is_absolute():
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
        )
    try:
        observed = path.lstat()
    except OSError:
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
        ) from None
    if (
        not stat.S_ISDIR(observed.st_mode)
        or stat.S_ISLNK(observed.st_mode)
        or observed.st_mode & 0o022
    ):
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
        )
    return path


def _scratch_directory(path: Path) -> Path:
    if not isinstance(path, Path) or not path.is_absolute():
        raise fail(
            "DEPLOYMENT_BUNDLE_OUTPUT_INVALID",
            "Deployment bundle output is invalid.",
        )
    try:
        observed = path.lstat()
    except OSError:
        raise fail(
            "DEPLOYMENT_BUNDLE_OUTPUT_INVALID",
            "Deployment bundle output is invalid.",
        ) from None
    if (
        not stat.S_ISDIR(observed.st_mode)
        or stat.S_ISLNK(observed.st_mode)
        or observed.st_uid != os.getuid()
        or stat.S_IMODE(observed.st_mode) & 0o077
    ):
        raise fail(
            "DEPLOYMENT_BUNDLE_OUTPUT_INVALID",
            "Deployment bundle output is invalid.",
        )
    return path


def _output_path(path: Path) -> Path:
    if (
        not isinstance(path, Path)
        or not path.is_absolute()
        or path.name in {"", ".", ".."}
    ):
        raise fail(
            "DEPLOYMENT_BUNDLE_OUTPUT_INVALID",
            "Deployment bundle output is invalid.",
        )
    try:
        parent = path.parent.lstat()
    except OSError:
        raise fail(
            "DEPLOYMENT_BUNDLE_OUTPUT_INVALID",
            "Deployment bundle output is invalid.",
        ) from None
    if (
        not stat.S_ISDIR(parent.st_mode)
        or stat.S_ISLNK(parent.st_mode)
        or parent.st_mode & 0o022
    ):
        raise fail(
            "DEPLOYMENT_BUNDLE_OUTPUT_INVALID",
            "Deployment bundle output is invalid.",
        )
    return path


def _open_output_directory(path: Path) -> int:
    flags = (
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    descriptor = -1
    try:
        before = path.lstat()
        descriptor = os.open(path, flags)
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISDIR(opened.st_mode)
            or stat.S_ISLNK(before.st_mode)
            or opened.st_mode & 0o022
            or (before.st_dev, before.st_ino) != (opened.st_dev, opened.st_ino)
        ):
            raise OSError
        return descriptor
    except OSError:
        if descriptor >= 0:
            os.close(descriptor)
        raise fail(
            "DEPLOYMENT_BUNDLE_OUTPUT_INVALID",
            "Deployment bundle output is invalid.",
        ) from None


def _hash_descriptor(descriptor: int) -> str:
    os.lseek(descriptor, 0, os.SEEK_SET)
    digest = hashlib.sha256()
    while True:
        chunk = os.read(descriptor, 1024 * 1024)
        if not chunk:
            break
        digest.update(chunk)
    return digest.hexdigest()


def _verify_path_md5(entry: _PathPayload, expected: str) -> None:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = -1
    try:
        descriptor = os.open(entry.source, flags)
        before = os.fstat(descriptor)
        if _file_witness(before) != entry.witness:
            raise OSError
        digest = hashlib.md5(usedforsecurity=False)
        while True:
            chunk = os.read(descriptor, 1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
        after = os.fstat(descriptor)
        if _file_witness(after) != entry.witness or digest.hexdigest() != expected:
            raise OSError
    except OSError:
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            "Deployment bundle source is invalid.",
            component=ENCODE_RUNTIME,
        ) from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _file_witness(observed: os.stat_result) -> tuple[int, ...]:
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


def _raw_identity(content: bytes) -> str:
    return f"sha256-{hashlib.sha256(content).hexdigest()}"


__all__ = [
    "build_bulk_rnaseq_runtime_bundle",
    "build_encode_runtime_bundle",
    "build_platform_bundle",
]
