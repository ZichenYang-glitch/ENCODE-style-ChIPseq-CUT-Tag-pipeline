"""Production resolution of deployment bindings to native product contracts.

The bundle manifest indexes bytes; it never becomes a second version catalog.
This module re-opens the indexed candidate bytes and delegates admission to the
existing migration, frontend, ENCODE build, and Bulk execution/runtime
contracts.  Resolved facts are process-local and are not serialized back into
the bundle.
"""

from __future__ import annotations

import ast
import base64
import csv
from dataclasses import dataclass
from email import policy
from email.parser import BytesParser
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import posixpath
import re
import sqlite3
import stat
import tempfile
from typing import Any, Callable
from urllib.parse import urlsplit
import zipfile

from encode_pipeline.adapters import EncodeStyleWorkflowAdapter
from encode_pipeline.adapters.bulk_rnaseq.execution_identity import (
    EXECUTION_IMPLEMENTATION_MANIFEST_FILE,
    verify_execution_implementation,
)
from encode_pipeline.adapters.bulk_rnaseq.qualification import (
    DEFAULT_EXECUTION_QUALIFICATION_FILE,
    load_default_execution_qualification,
)
from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import (
    RUNTIME_IDENTITY_FILE,
    RUNTIME_IDENTITY_SHA256,
    RuntimeAssetBinding,
    VerifiedRuntimeAssets,
    verify_runtime_asset_closure,
)
from encode_pipeline.deployment.admission import (
    DatabaseSchemaObservation,
    ResolvedContractFacts,
)
from encode_pipeline.deployment.canonical import (
    canonical_identity,
    canonical_json_bytes,
    without_key,
)
from encode_pipeline.deployment.database import DATABASE_CONTENT_IDENTITY_SCHEME
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.filesystem import read_regular_file
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import (
    BULK_RNASEQ_RUNTIME,
    ENCODE_RUNTIME,
    PLATFORM,
    BundleManifest,
    ContractIdentity,
    ContractRequirement,
    DeploymentState,
)
from encode_pipeline.deployment.platform_runtime import (
    PLATFORM_RUNTIME_LOCK_PATH,
    PlatformWheelLock,
    inspect_platform_runtime_closure,
)
from encode_pipeline.frontend_assets import parse_manifest_bytes
from encode_pipeline.persistence.migration_admission import (
    MIGRATION_EXECUTION_INVENTORY_FILE,
    MigrationInventoryTrustAnchor,
    VerifiedMigrationExecutionInventory,
    migration_inventory_trust_anchor_from_source,
    verify_migration_execution_inventory,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider


PLATFORM_DISTRIBUTION_CONTRACT = "helixweave.platform.distribution"
PLATFORM_WHEEL_CONTRACT = "helixweave.platform.python-wheel"
PLATFORM_PYTHON_RUNTIME_CONTRACT = "helixweave.platform.python-runtime"
PLATFORM_FRONTEND_CONTRACT = "helixweave.platform.frontend-assets"
PLATFORM_MIGRATIONS_CONTRACT = "helixweave.platform.database-migrations"
PLATFORM_REFERENCES_CONTRACT = "helixweave.platform.reference-compatibility"
ENCODE_RUNTIME_CONTRACT = "helixweave.runtime.encode"
BULK_RUNTIME_CONTRACT = "helixweave.runtime.bulk-rnaseq"

PLATFORM_METADATA_PATH = "payload/contracts/platform/distribution.METADATA"
PLATFORM_WHEEL_PATH = "payload/contracts/platform/helixweave.whl"
PLATFORM_FRONTEND_PATH = "payload/contracts/platform/frontend-asset-manifest.json"
PLATFORM_MIGRATIONS_PATH = (
    "payload/contracts/platform/migration-execution-inventory-1.0.0.json"
)
PLATFORM_ENCODE_REFERENCES_PATH = "payload/contracts/platform/encode_references.py"
PLATFORM_BULK_REFERENCES_PATH = "payload/contracts/platform/bulk_reference_profiles.py"
ENCODE_RUNTIME_ROOT_PATH = "payload/contracts/encode-runtime"
ENCODE_RUNTIME_ENTRYPOINT_PATH = f"{ENCODE_RUNTIME_ROOT_PATH}/workflow/Snakefile"
ENCODE_RUNTIME_INDEX_PATH = f"{ENCODE_RUNTIME_ROOT_PATH}/package-index.json"
ENCODE_RUNTIME_PAYLOAD_ROOT = "payload/runtime/encode"
ENCODE_MICROMAMBA_PATH = f"{ENCODE_RUNTIME_PAYLOAD_ROOT}/bin/micromamba"
ENCODE_PACKAGE_ARCHIVE_ROOT = f"{ENCODE_RUNTIME_PAYLOAD_ROOT}/packages/sha256"
BULK_RUNTIME_ROOT_PATH = "payload/runtime"
BULK_RUNTIME_IDENTITY_PATH = f"payload/contracts/bulk-rnaseq/{RUNTIME_IDENTITY_FILE}"

PLATFORM_WHEEL_MEMBER_FRONTEND = "encode_pipeline/frontend_assets/asset-manifest.json"
PLATFORM_WHEEL_MEMBER_MIGRATIONS = (
    f"encode_pipeline/persistence/{MIGRATION_EXECUTION_INVENTORY_FILE}"
)
PLATFORM_WHEEL_MEMBER_MIGRATION_ADMISSION = (
    "encode_pipeline/persistence/migration_admission.py"
)
PLATFORM_WHEEL_MEMBER_ENCODE_REFERENCES = (
    "encode_pipeline/adapters/encode_references.py"
)
PLATFORM_WHEEL_MEMBER_BULK_REFERENCES = (
    "encode_pipeline/adapters/bulk_rnaseq/reference_profiles.py"
)
PLATFORM_WHEEL_MEMBER_IMPLEMENTATION = (
    f"encode_pipeline/contracts/nfcore_rnaseq/{EXECUTION_IMPLEMENTATION_MANIFEST_FILE}"
)
PLATFORM_WHEEL_MEMBER_QUALIFICATION = (
    f"encode_pipeline/contracts/nfcore_rnaseq/{DEFAULT_EXECUTION_QUALIFICATION_FILE}"
)

REFERENCE_COMPATIBILITY_IDENTITY_SCHEME = (
    "helixweave-platform-reference-compatibility-v1"
)
ENCODE_RUNTIME_INDEX_SCHEMA = "helixweave-encode-runtime-package-index-v1"
ENCODE_RUNTIME_INDEX_IDENTITY_SCHEME = "helixweave-encode-runtime-closure-v1"

_MAX_CONTRACT_BYTES = 2 * 1024 * 1024
_MAX_WHEEL_BYTES = 2 * 1024**3
_MAX_WHEEL_FILES = 100_000
_MAX_WHEEL_MEMBER_BYTES = 512 * 1024 * 1024
_MAX_DATABASE_BYTES = 4 * 1024**4
_MAX_ENCODE_INDEX_BYTES = 32 * 1024 * 1024
_MAX_ENCODE_PACKAGES = 100_000
_IDENTITY = re.compile(r"^sha256-[0-9a-f]{64}$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_MD5 = re.compile(r"^[0-9a-f]{32}$")
# Conda package names may begin with an underscore (for example the
# ``_openmp_mutex`` and ``_python_abi3_support`` packages in the checked-in
# locks).  Keep dot and dash disallowed in the first position so a basename
# can never be confused with a relative path segment or an option.
_PACKAGE_FILENAME = re.compile(r"^[A-Za-z0-9_][A-Za-z0-9._+-]{0,254}$")
_SOURCE_ENV_PATH = re.compile(r"^workflow/envs/[A-Za-z0-9][A-Za-z0-9_.-]*\.yml$")
_CONDA_DIRECTIVE = re.compile(
    r"""(?m)^[ \t]*conda:[ \t]*(?:\n[ \t]+)?["']([^"']+\.yml)["'][ \t]*,?[ \t]*(?:\#.*)?$"""
)
_CONDA_LABEL = re.compile(r"(?m)^[ \t]*conda[ \t]*:")
_VERSION = re.compile(r"^[0-9A-Za-z][0-9A-Za-z._+-]{0,127}$")
_REVISION = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")


class _NativeContractFault(Exception):
    pass


@dataclass(frozen=True, order=True)
class ExplicitPackageCoordinate:
    """One exact conda archive coordinate from a linux-64 explicit lock."""

    url: str
    platform: str
    filename: str
    md5: str


@dataclass(frozen=True)
class EncodeRuntimePackageIndex:
    """Canonical offline package closure bound to one workflow build."""

    identity: str
    workflow_build_identity: str
    micromamba: dict[str, object]
    packages: tuple[dict[str, object], ...]
    environments: tuple[dict[str, object], ...]


def parse_encode_explicit_lock(content: bytes) -> tuple[ExplicitPackageCoordinate, ...]:
    """Parse the finite conda-lock explicit format accepted by deployment."""
    try:
        text = content.decode("utf-8")
    except UnicodeError:
        raise _NativeContractFault from None
    if "\r" in text or not text.endswith("\n"):
        raise _NativeContractFault
    lines = text.splitlines()
    if lines.count("@EXPLICIT") != 1 or "# platform: linux-64" not in lines:
        raise _NativeContractFault
    marker = lines.index("@EXPLICIT")
    if any(line and not line.startswith("#") for line in lines[:marker]):
        raise _NativeContractFault
    packages: list[ExplicitPackageCoordinate] = []
    for line in lines[marker + 1 :]:
        if not line or line.startswith("#"):
            raise _NativeContractFault
        parsed = urlsplit(line)
        parts = PurePosixPath(parsed.path).parts
        if (
            parsed.scheme != "https"
            or parsed.netloc != "conda.anaconda.org"
            or parsed.username is not None
            or parsed.password is not None
            or parsed.query
            or len(parts) < 4
            or parts[-2] not in {"linux-64", "noarch"}
            or _PACKAGE_FILENAME.fullmatch(parts[-1]) is None
            or not parts[-1].endswith((".conda", ".tar.bz2"))
            or _MD5.fullmatch(parsed.fragment) is None
        ):
            raise _NativeContractFault
        packages.append(
            ExplicitPackageCoordinate(
                url=line,
                platform=parts[-2],
                filename=parts[-1],
                md5=parsed.fragment,
            )
        )
    if not packages or len(packages) > _MAX_ENCODE_PACKAGES:
        raise _NativeContractFault
    urls = tuple(item.url for item in packages)
    if len(set(urls)) != len(urls):
        raise _NativeContractFault
    return tuple(packages)


def encode_environment_paths(
    source_manifest: tuple[tuple[str, bytes], ...],
) -> tuple[str, ...]:
    """Resolve static workflow conda directives plus the fixed runner env."""
    source = dict(source_manifest)
    if len(source) != len(source_manifest):
        raise _NativeContractFault
    environments = {"workflow/envs/runner.yml"}
    for logical_path, content in source_manifest:
        if not logical_path.startswith("workflow/rules/") or not logical_path.endswith(
            ".smk"
        ):
            continue
        try:
            text = content.decode("utf-8")
        except UnicodeError:
            raise _NativeContractFault from None
        matches = _CONDA_DIRECTIVE.findall(text)
        if len(matches) != len(_CONDA_LABEL.findall(text)):
            raise _NativeContractFault
        for relative in matches:
            resolved = posixpath.normpath(
                posixpath.join(posixpath.dirname(logical_path), relative)
            )
            if _SOURCE_ENV_PATH.fullmatch(resolved) is None:
                raise _NativeContractFault
            environments.add(resolved)
    for environment in environments:
        lock = environment[:-4] + ".lock"
        if environment not in source or lock not in source:
            raise _NativeContractFault
    return tuple(sorted(environments))


def verify_static_micromamba(content: bytes) -> None:
    """Require one static little-endian Linux x86_64 ELF executable."""
    if (
        len(content) < 120
        or content[:7] != b"\x7fELF\x02\x01\x01"
        or int.from_bytes(content[18:20], "little") != 62
    ):
        raise _NativeContractFault
    program_offset = int.from_bytes(content[32:40], "little")
    program_size = int.from_bytes(content[54:56], "little")
    program_count = int.from_bytes(content[56:58], "little")
    if (
        program_size != 56
        or not 1 <= program_count <= 128
        or program_offset < 64
        or program_offset + program_size * program_count > len(content)
    ):
        raise _NativeContractFault
    program_types = {
        int.from_bytes(
            content[
                program_offset + index * program_size : program_offset
                + index * program_size
                + 4
            ],
            "little",
        )
        for index in range(program_count)
    }
    if 1 not in program_types or 3 in program_types:
        raise _NativeContractFault


def encode_runtime_index_bytes(
    *,
    workflow_build_identity: str,
    micromamba: dict[str, object],
    packages: tuple[dict[str, object], ...],
    environments: tuple[dict[str, object], ...],
) -> bytes:
    """Create and self-validate the sole ENCODE native runtime contract."""
    value: dict[str, object] = {
        "schema_version": ENCODE_RUNTIME_INDEX_SCHEMA,
        "workflow_build_identity": workflow_build_identity,
        "micromamba": micromamba,
        "packages": list(packages),
        "environments": list(environments),
    }
    value["identity"] = canonical_identity(
        value,
        scheme=ENCODE_RUNTIME_INDEX_IDENTITY_SCHEME,
    )
    content = canonical_json_bytes(value)
    parse_encode_runtime_index(content)
    return content


def parse_encode_runtime_index(content: bytes) -> EncodeRuntimePackageIndex:
    """Strictly parse the canonical ENCODE offline closure index."""
    if not 0 < len(content) <= _MAX_ENCODE_INDEX_BYTES:
        raise _NativeContractFault
    raw = _strict_json(content)
    if (
        not isinstance(raw, dict)
        or set(raw)
        != {
            "schema_version",
            "identity",
            "workflow_build_identity",
            "micromamba",
            "packages",
            "environments",
        }
        or raw["schema_version"] != ENCODE_RUNTIME_INDEX_SCHEMA
        or canonical_json_bytes(raw) != content
        or not isinstance(raw["identity"], str)
        or raw["identity"]
        != canonical_identity(
            without_key(raw, "identity"),
            scheme=ENCODE_RUNTIME_INDEX_IDENTITY_SCHEME,
        )
        or not isinstance(raw["workflow_build_identity"], str)
        or _IDENTITY.fullmatch(raw["workflow_build_identity"]) is None
    ):
        raise _NativeContractFault
    micromamba = raw["micromamba"]
    if (
        not isinstance(micromamba, dict)
        or set(micromamba) != {"path", "size_bytes", "sha256"}
        or micromamba["path"] != ENCODE_MICROMAMBA_PATH
        or not isinstance(micromamba["size_bytes"], int)
        or not 0 < micromamba["size_bytes"] <= _MAX_WHEEL_BYTES
        or not isinstance(micromamba["sha256"], str)
        or _SHA256.fullmatch(micromamba["sha256"]) is None
    ):
        raise _NativeContractFault
    packages_raw = raw["packages"]
    if (
        not isinstance(packages_raw, list)
        or not 0 < len(packages_raw) <= _MAX_ENCODE_PACKAGES
    ):
        raise _NativeContractFault
    packages: list[dict[str, object]] = []
    coordinates: dict[str, ExplicitPackageCoordinate] = {}
    basenames: dict[str, str] = {}
    for item in packages_raw:
        if not isinstance(item, dict) or set(item) != {
            "url",
            "platform",
            "filename",
            "md5",
            "archive_path",
            "size_bytes",
            "sha256",
        }:
            raise _NativeContractFault
        coordinate = parse_encode_explicit_lock(
            b"# platform: linux-64\n@EXPLICIT\n"
            + str(item["url"]).encode("utf-8")
            + b"\n"
        )[0]
        expected_archive = (
            f"{ENCODE_PACKAGE_ARCHIVE_ROOT}/{item['sha256']}/{coordinate.filename}"
        )
        if (
            item["platform"] != coordinate.platform
            or item["filename"] != coordinate.filename
            or item["md5"] != coordinate.md5
            or item["archive_path"] != expected_archive
            or not isinstance(item["size_bytes"], int)
            or not 0 < item["size_bytes"] <= _MAX_WHEEL_BYTES
            or not isinstance(item["sha256"], str)
            or _SHA256.fullmatch(item["sha256"]) is None
            or coordinate.url in coordinates
            or (
                coordinate.filename in basenames
                and basenames[coordinate.filename] != coordinate.url
            )
        ):
            raise _NativeContractFault
        coordinates[coordinate.url] = coordinate
        basenames[coordinate.filename] = coordinate.url
        packages.append(item)
    if packages != sorted(packages, key=lambda item: item["url"]):
        raise _NativeContractFault
    environments_raw = raw["environments"]
    if not isinstance(environments_raw, list) or not environments_raw:
        raise _NativeContractFault
    environments: list[dict[str, object]] = []
    used_urls: set[str] = set()
    for item in environments_raw:
        if not isinstance(item, dict) or set(item) != {
            "environment_path",
            "environment_sha256",
            "lock_path",
            "lock_sha256",
            "packages",
        }:
            raise _NativeContractFault
        environment_path = item["environment_path"]
        expected_lock = (
            environment_path[:-4] + ".lock"
            if isinstance(environment_path, str)
            else None
        )
        urls = item["packages"]
        if (
            not isinstance(environment_path, str)
            or _SOURCE_ENV_PATH.fullmatch(environment_path) is None
            or item["lock_path"] != expected_lock
            or not isinstance(item["environment_sha256"], str)
            or _SHA256.fullmatch(item["environment_sha256"]) is None
            or not isinstance(item["lock_sha256"], str)
            or _SHA256.fullmatch(item["lock_sha256"]) is None
            or not isinstance(urls, list)
            or not urls
            or urls != sorted(set(urls))
            or any(url not in coordinates for url in urls)
        ):
            raise _NativeContractFault
        used_urls.update(urls)
        environments.append(item)
    if (
        environments != sorted(environments, key=lambda item: item["environment_path"])
        or len({item["environment_path"] for item in environments}) != len(environments)
        or "workflow/envs/runner.yml"
        not in {item["environment_path"] for item in environments}
        or used_urls != set(coordinates)
    ):
        raise _NativeContractFault
    return EncodeRuntimePackageIndex(
        identity=raw["identity"],
        workflow_build_identity=raw["workflow_build_identity"],
        micromamba=micromamba,
        packages=tuple(packages),
        environments=tuple(environments),
    )


@dataclass(frozen=True)
class VerifiedPlatformWheel:
    """Facts admitted from one exact candidate wheel."""

    version: str
    metadata: bytes
    frontend_manifest: bytes
    migration_inventory: bytes
    encode_reference_source: bytes
    bulk_reference_source: bytes
    frontend_identity: str
    migration: VerifiedMigrationExecutionInventory
    reference_identity: str
    bulk_runtime_identity: str


BulkRuntimeVerifier = Callable[[RuntimeAssetBinding], object]


class ProductionNativeContractResolver:
    """Resolve the finite component set from candidate-owned native bytes."""

    def __init__(
        self,
        *,
        docker_executable: Path = Path("/usr/bin/docker"),
        docker_socket: Path = Path("/run/helixweave/docker/docker.sock"),
        bulk_runtime_verifier: BulkRuntimeVerifier = verify_runtime_asset_closure,
    ) -> None:
        if (
            not isinstance(docker_executable, Path)
            or not docker_executable.is_absolute()
            or not isinstance(docker_socket, Path)
            or not docker_socket.is_absolute()
            or not callable(bulk_runtime_verifier)
        ):
            raise fail(
                "DEPLOYMENT_CONTRACT_RESOLVER_INVALID",
                "Deployment contract resolver is invalid.",
            )
        self._docker_executable = docker_executable
        self._docker_socket = docker_socket
        self._bulk_runtime_verifier = bulk_runtime_verifier

    def resolve(
        self,
        root: Path,
        manifest: BundleManifest,
    ) -> ResolvedContractFacts:
        if not isinstance(root, Path) or not root.is_absolute():
            raise fail(
                "DEPLOYMENT_CONTRACT_ADMISSION_FAILED",
                "Deployment contract admission failed.",
                component=getattr(manifest, "component", None),
            )
        try:
            if manifest.component == PLATFORM:
                return self._resolve_platform(root, manifest)
            if manifest.component == ENCODE_RUNTIME:
                return self._resolve_encode(root, manifest)
            if manifest.component == BULK_RNASEQ_RUNTIME:
                return self._resolve_bulk(root, manifest)
        except DeploymentError:
            raise
        except Exception:
            pass
        raise fail(
            "DEPLOYMENT_CONTRACT_ADMISSION_FAILED",
            "Deployment contract admission failed.",
            component=(
                manifest.component if isinstance(manifest, BundleManifest) else None
            ),
        )

    def _resolve_platform(
        self,
        root: Path,
        manifest: BundleManifest,
    ) -> ResolvedContractFacts:
        metadata = _contract_bytes(
            root,
            manifest,
            PLATFORM_DISTRIBUTION_CONTRACT,
            PLATFORM_METADATA_PATH,
            maximum_bytes=_MAX_CONTRACT_BYTES,
        )
        wheel = _contract_bytes(
            root,
            manifest,
            PLATFORM_WHEEL_CONTRACT,
            PLATFORM_WHEEL_PATH,
            maximum_bytes=_MAX_WHEEL_BYTES,
        )
        runtime_lock_content = _contract_bytes(
            root,
            manifest,
            PLATFORM_PYTHON_RUNTIME_CONTRACT,
            PLATFORM_RUNTIME_LOCK_PATH,
            maximum_bytes=_MAX_CONTRACT_BYTES,
        )
        frontend = _contract_bytes(
            root,
            manifest,
            PLATFORM_FRONTEND_CONTRACT,
            PLATFORM_FRONTEND_PATH,
            maximum_bytes=_MAX_CONTRACT_BYTES,
        )
        migrations = _contract_bytes(
            root,
            manifest,
            PLATFORM_MIGRATIONS_CONTRACT,
            PLATFORM_MIGRATIONS_PATH,
            maximum_bytes=_MAX_CONTRACT_BYTES,
        )
        encode_references = _contract_bytes(
            root,
            manifest,
            PLATFORM_REFERENCES_CONTRACT,
            PLATFORM_ENCODE_REFERENCES_PATH,
            maximum_bytes=_MAX_CONTRACT_BYTES,
        )
        bulk_references = _indexed_bytes(
            root,
            manifest,
            PLATFORM_BULK_REFERENCES_PATH,
            maximum_bytes=_MAX_CONTRACT_BYTES,
        )
        facts = verify_platform_wheel(wheel)
        runtime_lock = PlatformWheelLock.from_bytes(runtime_lock_content)
        wheel_sha256 = hashlib.sha256(wheel).hexdigest()
        if sum(item.sha256 == wheel_sha256 for item in runtime_lock.wheels) != 1:
            raise _NativeContractFault
        runtime_closure = inspect_platform_runtime_closure(
            root,
            expected_lock_identity=runtime_lock.identity,
        )
        _verify_platform_runtime_records(manifest, runtime_closure.files)
        if (
            metadata != facts.metadata
            or frontend != facts.frontend_manifest
            or migrations != facts.migration_inventory
            or encode_references != facts.encode_reference_source
            or bulk_references != facts.bulk_reference_source
        ):
            raise _NativeContractFault
        contracts = tuple(
            sorted(
                (
                    ContractIdentity(
                        PLATFORM_DISTRIBUTION_CONTRACT,
                        _raw_identity(metadata),
                    ),
                    ContractIdentity(PLATFORM_WHEEL_CONTRACT, _raw_identity(wheel)),
                    ContractIdentity(
                        PLATFORM_PYTHON_RUNTIME_CONTRACT,
                        runtime_closure.lock_identity,
                    ),
                    ContractIdentity(
                        PLATFORM_FRONTEND_CONTRACT,
                        facts.frontend_identity,
                    ),
                    ContractIdentity(
                        PLATFORM_MIGRATIONS_CONTRACT,
                        _raw_identity(migrations),
                    ),
                    ContractIdentity(
                        PLATFORM_REFERENCES_CONTRACT,
                        facts.reference_identity,
                    ),
                )
            )
        )
        return ResolvedContractFacts(
            component=PLATFORM,
            deployment_identity=manifest.identity,
            version=facts.version,
            contracts=contracts,
            requirements=(
                ContractRequirement(
                    BULK_RUNTIME_CONTRACT,
                    (facts.bulk_runtime_identity,),
                ),
            ),
            database_heads=facts.migration.heads,
        )

    def _resolve_encode(
        self,
        root: Path,
        manifest: BundleManifest,
    ) -> ResolvedContractFacts:
        index_content = _contract_bytes(
            root,
            manifest,
            ENCODE_RUNTIME_CONTRACT,
            ENCODE_RUNTIME_INDEX_PATH,
            maximum_bytes=_MAX_ENCODE_INDEX_BYTES,
        )
        index = parse_encode_runtime_index(index_content)
        project_root = root.joinpath(*PurePosixPath(ENCODE_RUNTIME_ROOT_PATH).parts)
        adapter = EncodeStyleWorkflowAdapter(
            catalog_path=str(
                project_root / "docs" / "architecture" / "artifact-inventory.yaml"
            )
        )
        registry = WorkflowRegistry(
            adapters=(adapter,),
            legacy_execution_fallbacks=(adapter,),
        )
        provider = WorkflowBuildIdentityProvider(registry, project_root=project_root)
        captured = provider.capture(adapter.metadata.workflow_id)
        if captured.is_failure:
            raise _NativeContractFault
        identity = captured.value
        source_manifest = provider.source_manifest()
        source_paths = {
            f"{ENCODE_RUNTIME_ROOT_PATH}/{logical_path}"
            for logical_path, _content in source_manifest
        }
        observed = f"sha256-{identity.digest}"
        if (
            _IDENTITY.fullmatch(observed) is None
            or index.workflow_build_identity != observed
        ):
            raise _NativeContractFault
        source = dict(source_manifest)
        environments = encode_environment_paths(source_manifest)
        if (
            tuple(item["environment_path"] for item in index.environments)
            != environments
        ):
            raise _NativeContractFault
        indexed_packages = {item["url"]: item for item in index.packages}
        for environment in index.environments:
            environment_path = environment["environment_path"]
            lock_path = environment["lock_path"]
            environment_bytes = source[environment_path]
            lock_bytes = source[lock_path]
            coordinates = parse_encode_explicit_lock(lock_bytes)
            if (
                environment["environment_sha256"]
                != hashlib.sha256(environment_bytes).hexdigest()
                or environment["lock_sha256"] != hashlib.sha256(lock_bytes).hexdigest()
                or environment["packages"] != sorted(item.url for item in coordinates)
            ):
                raise _NativeContractFault
        micromamba_path = index.micromamba["path"]
        micromamba = _indexed_bytes(
            root,
            manifest,
            micromamba_path,
            maximum_bytes=_MAX_WHEEL_BYTES,
        )
        micromamba_record = _single_record(manifest, micromamba_path)
        verify_static_micromamba(micromamba)
        if (
            micromamba_record.mode != 0o555
            or len(micromamba) != index.micromamba["size_bytes"]
            or hashlib.sha256(micromamba).hexdigest() != index.micromamba["sha256"]
        ):
            raise _NativeContractFault
        package_paths: set[str] = set()
        for url, package in indexed_packages.items():
            archive_path = package["archive_path"]
            archive = _indexed_bytes(
                root,
                manifest,
                archive_path,
                maximum_bytes=_MAX_WHEEL_BYTES,
            )
            if (
                _single_record(manifest, archive_path).mode != 0o444
                or len(archive) != package["size_bytes"]
                or hashlib.sha256(archive).hexdigest() != package["sha256"]
                or hashlib.md5(archive, usedforsecurity=False).hexdigest()
                != package["md5"]
                or parse_encode_explicit_lock(
                    b"# platform: linux-64\n@EXPLICIT\n" + url.encode("utf-8") + b"\n"
                )[0].filename
                != package["filename"]
            ):
                raise _NativeContractFault
            package_paths.add(archive_path)
        expected_paths = source_paths | {
            ENCODE_RUNTIME_INDEX_PATH,
            micromamba_path,
            *package_paths,
        }
        if {record.path for record in manifest.files} != expected_paths:
            raise _NativeContractFault
        return ResolvedContractFacts(
            component=ENCODE_RUNTIME,
            deployment_identity=manifest.identity,
            version=identity.adapter_version,
            contracts=(ContractIdentity(ENCODE_RUNTIME_CONTRACT, index.identity),),
        )

    def _resolve_bulk(
        self,
        root: Path,
        manifest: BundleManifest,
    ) -> ResolvedContractFacts:
        runtime_identity = _contract_bytes(
            root,
            manifest,
            BULK_RUNTIME_CONTRACT,
            BULK_RUNTIME_IDENTITY_PATH,
            maximum_bytes=_MAX_CONTRACT_BYTES,
        )
        if hashlib.sha256(runtime_identity).hexdigest() != RUNTIME_IDENTITY_SHA256:
            raise _NativeContractFault
        document = _strict_json(runtime_identity)
        source = document.get("source") if isinstance(document, dict) else None
        version = source.get("release") if isinstance(source, dict) else None
        if not isinstance(version, str) or _VERSION.fullmatch(version) is None:
            raise _NativeContractFault
        allowed_prefix = f"{BULK_RUNTIME_ROOT_PATH}/"
        if any(
            record.path != BULK_RUNTIME_IDENTITY_PATH
            and not record.path.startswith(allowed_prefix)
            for record in manifest.files
        ):
            raise _NativeContractFault
        binding = RuntimeAssetBinding(
            root=root.joinpath(*PurePosixPath(BULK_RUNTIME_ROOT_PATH).parts),
            docker_executable=self._docker_executable,
            docker_socket=self._docker_socket,
        )
        result = self._bulk_runtime_verifier(binding)
        if (
            not hasattr(result, "is_success")
            or not result.is_success
            or not isinstance(result.value, VerifiedRuntimeAssets)
            or result.value.runtime_identity_sha256 != RUNTIME_IDENTITY_SHA256
        ):
            raise _NativeContractFault
        return ResolvedContractFacts(
            component=BULK_RNASEQ_RUNTIME,
            deployment_identity=manifest.identity,
            version=version,
            contracts=(
                ContractIdentity(
                    BULK_RUNTIME_CONTRACT, _raw_identity(runtime_identity)
                ),
            ),
        )


class InventoryDatabaseSchemaObserver:
    """Observe one trusted SQLite path and bind heads to the native inventory."""

    def __init__(self, layout: DeploymentLayout) -> None:
        if not isinstance(layout, DeploymentLayout):
            raise fail(
                "DEPLOYMENT_SCHEMA_OBSERVER_INVALID",
                "Database schema observer is invalid.",
            )
        self._layout = layout

    def observe(self, state: DeploymentState) -> DatabaseSchemaObservation:
        if not isinstance(state, DeploymentState):
            raise fail(
                "DEPLOYMENT_SCHEMA_OBSERVATION_FAILED",
                "Database schema observation failed.",
            )
        inventory = verify_migration_execution_inventory()
        heads, database_identity = _observe_live_database(self._layout.database)
        return DatabaseSchemaObservation.create(
            provider_identity=f"sha256-{inventory.contract_sha256}",
            state_identity=state.identity,
            database_identity=database_identity,
            heads=heads,
        )


def verify_platform_wheel(wheel: bytes) -> VerifiedPlatformWheel:
    """Verify one wheel and re-use every existing native product admission."""
    if not isinstance(wheel, bytes) or not 0 < len(wheel) <= _MAX_WHEEL_BYTES:
        raise _NativeContractFault
    with zipfile.ZipFile(io.BytesIO(wheel), mode="r") as archive:
        members = _verify_wheel_record(archive)
        metadata_path = next(
            (
                name
                for name in members
                if name.endswith(".dist-info/METADATA") and name.count("/") == 1
            ),
            None,
        )
        if (
            metadata_path is None
            or sum(name.endswith(".dist-info/METADATA") for name in members) != 1
        ):
            raise _NativeContractFault
        metadata = _zip_bytes(archive, metadata_path, _MAX_CONTRACT_BYTES)
        version = _distribution_version(metadata)
        frontend = _zip_bytes(
            archive,
            PLATFORM_WHEEL_MEMBER_FRONTEND,
            _MAX_CONTRACT_BYTES,
        )
        migrations = _zip_bytes(
            archive,
            PLATFORM_WHEEL_MEMBER_MIGRATIONS,
            _MAX_CONTRACT_BYTES,
        )
        migration_admission = _zip_bytes(
            archive,
            PLATFORM_WHEEL_MEMBER_MIGRATION_ADMISSION,
            _MAX_CONTRACT_BYTES,
        )
        migration_trust_anchor = migration_inventory_trust_anchor_from_source(
            migration_admission
        )
        encode_references = _zip_bytes(
            archive,
            PLATFORM_WHEEL_MEMBER_ENCODE_REFERENCES,
            _MAX_CONTRACT_BYTES,
        )
        bulk_references = _zip_bytes(
            archive,
            PLATFORM_WHEEL_MEMBER_BULK_REFERENCES,
            _MAX_CONTRACT_BYTES,
        )
        implementation_manifest = _zip_bytes(
            archive,
            PLATFORM_WHEEL_MEMBER_IMPLEMENTATION,
            _MAX_CONTRACT_BYTES,
        )
        qualification_content = _zip_bytes(
            archive,
            PLATFORM_WHEEL_MEMBER_QUALIFICATION,
            _MAX_CONTRACT_BYTES,
        )
        package_root = zipfile.Path(archive, at="encode_pipeline/")
        implementation = verify_execution_implementation(
            manifest_bytes=implementation_manifest,
            package_root=package_root,
        )
        if implementation.is_failure:
            raise _NativeContractFault
        qualification = load_default_execution_qualification(
            implementation.value,
            content=qualification_content,
        )
        if qualification.is_failure:
            raise _NativeContractFault
        qualification_document = _strict_json(qualification_content)
        runtime = (
            qualification_document.get("runtime")
            if isinstance(qualification_document, dict)
            else None
        )
        runtime_sha256 = (
            runtime.get("identity_sha256") if isinstance(runtime, dict) else None
        )
        if (
            not isinstance(runtime_sha256, str)
            or re.fullmatch(r"[0-9a-f]{64}", runtime_sha256) is None
        ):
            raise _NativeContractFault
        migration = _verify_wheel_migrations(
            archive,
            migrations,
            trust_anchor=migration_trust_anchor,
        )
    frontend_manifest = parse_manifest_bytes(frontend)
    reference_identity = _reference_compatibility_identity(
        encode_references,
        bulk_references,
    )
    return VerifiedPlatformWheel(
        version=version,
        metadata=metadata,
        frontend_manifest=frontend,
        migration_inventory=migrations,
        encode_reference_source=encode_references,
        bulk_reference_source=bulk_references,
        frontend_identity=frontend_manifest.identity,
        migration=migration,
        reference_identity=reference_identity,
        bulk_runtime_identity=f"sha256-{runtime_sha256}",
    )


def _verify_wheel_migrations(
    archive: zipfile.ZipFile,
    inventory: bytes,
    *,
    trust_anchor: MigrationInventoryTrustAnchor,
) -> VerifiedMigrationExecutionInventory:
    prefix = "encode_pipeline/persistence/"
    with tempfile.TemporaryDirectory(prefix="helixweave-wheel-migrations-") as value:
        root = Path(value)
        for member in archive.infolist():
            if member.is_dir() or not member.filename.startswith(prefix):
                continue
            relative = member.filename.removeprefix(prefix)
            if relative == MIGRATION_EXECUTION_INVENTORY_FILE or relative.startswith(
                "alembic/"
            ):
                target = root.joinpath(*PurePosixPath(relative).parts)
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(
                    _zip_bytes(
                        archive,
                        member.filename,
                        _MAX_WHEEL_MEMBER_BYTES,
                        allow_empty=True,
                    )
                )
        return verify_migration_execution_inventory(
            persistence_root=root,
            inventory_bytes=inventory,
            trust_anchor=trust_anchor,
        )


def _verify_wheel_record(archive: zipfile.ZipFile) -> frozenset[str]:
    infos = archive.infolist()
    if not infos or len(infos) > _MAX_WHEEL_FILES:
        raise _NativeContractFault
    names: set[str] = set()
    regular: dict[str, zipfile.ZipInfo] = {}
    total = 0
    for info in infos:
        name = info.filename
        path = PurePosixPath(name.rstrip("/"))
        if (
            not name
            or name in names
            or path.is_absolute()
            or any(part in {"", ".", ".."} for part in path.parts)
            or "\\" in name
            or info.file_size < 0
            or info.file_size > _MAX_WHEEL_MEMBER_BYTES
        ):
            raise _NativeContractFault
        names.add(name)
        mode = info.external_attr >> 16
        if mode and stat.S_ISLNK(mode):
            raise _NativeContractFault
        if info.is_dir():
            continue
        total += info.file_size
        if total > _MAX_WHEEL_BYTES:
            raise _NativeContractFault
        regular[name] = info
    record_names = [
        name
        for name in regular
        if name.endswith(".dist-info/RECORD") and name.count("/") == 1
    ]
    if len(record_names) != 1:
        raise _NativeContractFault
    record_name = record_names[0]
    try:
        rows = tuple(
            csv.reader(
                io.StringIO(
                    _zip_bytes(archive, record_name, _MAX_CONTRACT_BYTES).decode(
                        "utf-8"
                    ),
                    newline="",
                )
            )
        )
    except (UnicodeError, csv.Error):
        raise _NativeContractFault from None
    declared: set[str] = set()
    for row in rows:
        if len(row) != 3:
            raise _NativeContractFault
        name, encoded_digest, encoded_size = row
        if name in declared or name not in regular:
            raise _NativeContractFault
        declared.add(name)
        if name == record_name:
            if encoded_digest or encoded_size:
                raise _NativeContractFault
            continue
        if not encoded_digest.startswith("sha256=") or not encoded_size.isdecimal():
            raise _NativeContractFault
        content = _zip_bytes(
            archive,
            name,
            _MAX_WHEEL_MEMBER_BYTES,
            allow_empty=True,
        )
        digest = base64.urlsafe_b64encode(hashlib.sha256(content).digest()).rstrip(b"=")
        try:
            declared_digest = encoded_digest.removeprefix("sha256=").encode("ascii")
        except UnicodeEncodeError:
            raise _NativeContractFault from None
        if declared_digest != digest or int(encoded_size) != len(content):
            raise _NativeContractFault
    if declared != set(regular):
        raise _NativeContractFault
    return frozenset(regular)


def _distribution_version(metadata: bytes) -> str:
    try:
        document = BytesParser(policy=policy.compat32).parsebytes(metadata)
    except Exception:
        raise _NativeContractFault from None
    names = document.get_all("Name", failobj=[])
    versions = document.get_all("Version", failobj=[])
    if (
        names != ["helixweave"]
        or len(versions) != 1
        or _VERSION.fullmatch(versions[0]) is None
    ):
        raise _NativeContractFault
    return versions[0]


def _reference_compatibility_identity(
    encode_source: bytes,
    bulk_source: bytes,
) -> str:
    values = {
        "encode": _python_string_constant(
            encode_source,
            "ENCODE_REFERENCE_BINDING_CONTRACT",
        ),
        "bulk_rnaseq": _python_string_constant(
            bulk_source,
            "BULK_RNASEQ_REFERENCE_BINDING_CONTRACT",
        ),
    }
    return canonical_identity(
        values,
        scheme=REFERENCE_COMPATIBILITY_IDENTITY_SCHEME,
    )


def _verify_platform_runtime_records(
    manifest: BundleManifest,
    files: object,
) -> None:
    try:
        captured = tuple(files)
        expected = {
            item.logical_path: (item.size_bytes, item.sha256, item.mode)
            for item in captured
        }
    except (AttributeError, TypeError):
        raise _NativeContractFault from None
    if not expected or len(expected) != len(captured):
        raise _NativeContractFault
    actual = {
        item.path: (item.size_bytes, item.sha256, item.mode)
        for item in manifest.files
        if item.path == PLATFORM_RUNTIME_LOCK_PATH
        or item.path.startswith("payload/platform/")
    }
    if actual != expected:
        raise _NativeContractFault


def _python_string_constant(source: bytes, name: str) -> str:
    try:
        tree = ast.parse(source.decode("utf-8"))
    except (SyntaxError, UnicodeError):
        raise _NativeContractFault from None
    values: list[str] = []
    for statement in tree.body:
        target = None
        value = None
        if isinstance(statement, ast.Assign) and len(statement.targets) == 1:
            target = statement.targets[0]
            value = statement.value
        elif isinstance(statement, ast.AnnAssign):
            target = statement.target
            value = statement.value
        if isinstance(target, ast.Name) and target.id == name:
            try:
                literal = ast.literal_eval(value)
            except (ValueError, TypeError):
                raise _NativeContractFault from None
            if not isinstance(literal, str) or _VERSION.fullmatch(literal) is None:
                raise _NativeContractFault
            values.append(literal)
    if len(values) != 1:
        raise _NativeContractFault
    return values[0]


def _observe_live_database(path: Path) -> tuple[tuple[str, ...], str]:
    try:
        parent = path.parent.lstat()
        before = path.lstat()
        if (
            not stat.S_ISDIR(parent.st_mode)
            or stat.S_ISLNK(parent.st_mode)
            or parent.st_mode & 0o022
            or not stat.S_ISREG(before.st_mode)
            or stat.S_ISLNK(before.st_mode)
            or before.st_nlink != 1
            or before.st_uid != parent.st_uid
            or before.st_gid != parent.st_gid
            or before.st_mode & 0o022
            or not 0 < before.st_size <= _MAX_DATABASE_BYTES
        ):
            raise _NativeContractFault
        connection = sqlite3.connect(
            f"{path.as_uri()}?mode=ro",
            uri=True,
            timeout=1.0,
        )
        try:
            connection.execute("PRAGMA query_only=ON")
            integrity = tuple(
                row[0] for row in connection.execute("PRAGMA integrity_check(1)")
            )
            heads = tuple(
                row[0]
                for row in connection.execute(
                    "SELECT version_num FROM alembic_version ORDER BY version_num"
                )
            )
        finally:
            connection.close()
        after = path.lstat()
        if (
            (before.st_dev, before.st_ino) != (after.st_dev, after.st_ino)
            or integrity != ("ok",)
            or not heads
            or tuple(sorted(set(heads))) != heads
            or any(
                not isinstance(head, str) or _REVISION.fullmatch(head) is None
                for head in heads
            )
        ):
            raise _NativeContractFault
    except (OSError, sqlite3.DatabaseError, _NativeContractFault):
        raise fail(
            "DEPLOYMENT_SCHEMA_OBSERVATION_FAILED",
            "Database schema observation failed.",
        ) from None
    database_identity = canonical_identity(
        {
            "device": before.st_dev,
            "inode": before.st_ino,
            "heads": list(heads),
        },
        scheme=DATABASE_CONTENT_IDENTITY_SCHEME,
    )
    return heads, database_identity


def _contract_bytes(
    root: Path,
    manifest: BundleManifest,
    contract: str,
    expected_path: str,
    *,
    maximum_bytes: int,
) -> bytes:
    bindings = [item for item in manifest.contracts if item.contract == contract]
    if len(bindings) != 1 or bindings[0].path != expected_path:
        raise _NativeContractFault
    return _indexed_bytes(
        root,
        manifest,
        expected_path,
        maximum_bytes=maximum_bytes,
    )


def _indexed_bytes(
    root: Path,
    manifest: BundleManifest,
    logical_path: str,
    *,
    maximum_bytes: int,
) -> bytes:
    records = [item for item in manifest.files if item.path == logical_path]
    if len(records) != 1 or records[0].size_bytes > maximum_bytes:
        raise _NativeContractFault
    content, observed = read_regular_file(
        root.joinpath(*PurePosixPath(logical_path).parts),
        max_bytes=maximum_bytes,
        code="DEPLOYMENT_CONTRACT_ADMISSION_FAILED",
    )
    record = records[0]
    if (
        len(content) != record.size_bytes
        or hashlib.sha256(content).hexdigest() != record.sha256
        or stat.S_IMODE(observed.st_mode) != record.mode
    ):
        raise _NativeContractFault
    return content


def _single_record(manifest: BundleManifest, logical_path: object):
    if not isinstance(logical_path, str):
        raise _NativeContractFault
    records = [item for item in manifest.files if item.path == logical_path]
    if len(records) != 1:
        raise _NativeContractFault
    return records[0]


def _zip_bytes(
    archive: zipfile.ZipFile,
    name: str,
    maximum_bytes: int,
    *,
    allow_empty: bool = False,
) -> bytes:
    try:
        info = archive.getinfo(name)
    except KeyError:
        raise _NativeContractFault from None
    if (
        info.is_dir()
        or not 0 <= info.file_size <= maximum_bytes
        or (not allow_empty and info.file_size == 0)
    ):
        raise _NativeContractFault
    try:
        content = archive.read(info)
    except (OSError, RuntimeError, zipfile.BadZipFile):
        raise _NativeContractFault from None
    if len(content) != info.file_size:
        raise _NativeContractFault
    return content


def _strict_json(content: bytes) -> Any:
    def unique(pairs: list[tuple[str, object]]) -> dict[str, object]:
        value: dict[str, object] = {}
        for key, item in pairs:
            if key in value:
                raise ValueError
            value[key] = item
        return value

    try:
        return json.loads(
            content.decode("utf-8"),
            object_pairs_hook=unique,
            parse_constant=lambda _value: (_ for _ in ()).throw(ValueError()),
        )
    except (UnicodeError, ValueError, json.JSONDecodeError):
        raise _NativeContractFault from None


def _raw_identity(content: bytes) -> str:
    identity = f"sha256-{hashlib.sha256(content).hexdigest()}"
    if _IDENTITY.fullmatch(identity) is None:
        raise _NativeContractFault
    return identity


__all__ = [
    "BULK_RUNTIME_CONTRACT",
    "BULK_RUNTIME_IDENTITY_PATH",
    "BULK_RUNTIME_ROOT_PATH",
    "ENCODE_RUNTIME_CONTRACT",
    "ENCODE_RUNTIME_ENTRYPOINT_PATH",
    "ENCODE_RUNTIME_INDEX_PATH",
    "ENCODE_RUNTIME_INDEX_SCHEMA",
    "ENCODE_RUNTIME_PAYLOAD_ROOT",
    "ENCODE_RUNTIME_ROOT_PATH",
    "ENCODE_MICROMAMBA_PATH",
    "ENCODE_PACKAGE_ARCHIVE_ROOT",
    "EncodeRuntimePackageIndex",
    "ExplicitPackageCoordinate",
    "InventoryDatabaseSchemaObserver",
    "PLATFORM_BULK_REFERENCES_PATH",
    "PLATFORM_DISTRIBUTION_CONTRACT",
    "PLATFORM_ENCODE_REFERENCES_PATH",
    "PLATFORM_FRONTEND_CONTRACT",
    "PLATFORM_FRONTEND_PATH",
    "PLATFORM_METADATA_PATH",
    "PLATFORM_MIGRATIONS_CONTRACT",
    "PLATFORM_MIGRATIONS_PATH",
    "PLATFORM_PYTHON_RUNTIME_CONTRACT",
    "PLATFORM_REFERENCES_CONTRACT",
    "PLATFORM_WHEEL_CONTRACT",
    "PLATFORM_WHEEL_PATH",
    "ProductionNativeContractResolver",
    "VerifiedPlatformWheel",
    "encode_environment_paths",
    "encode_runtime_index_bytes",
    "parse_encode_explicit_lock",
    "parse_encode_runtime_index",
    "verify_static_micromamba",
    "verify_platform_wheel",
]
