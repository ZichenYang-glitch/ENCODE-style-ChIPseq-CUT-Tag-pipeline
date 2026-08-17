"""Canonical identity and fail-closed loading for production frontend bytes."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import hashlib
from importlib import resources
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
from types import MappingProxyType
from typing import Any

from encode_pipeline.deployment.canonical import (
    canonical_identity,
    canonical_json_bytes,
    without_key,
)
from encode_pipeline.deployment.errors import DeploymentError, fail


FRONTEND_ASSET_MANIFEST_SCHEMA = "helixweave-frontend-assets-v1"
FRONTEND_ASSET_IDENTITY_SCHEME = "helixweave-frontend-assets-identity-v1"
MANIFEST_FILENAME = "asset-manifest.json"
STATIC_DIRECTORY = "static"
ENTRYPOINT = "index.html"

_MANIFEST_MAX_BYTES = 1024 * 1024
_METADATA_MAX_BYTES = 16 * 1024 * 1024
_FILE_MAX_BYTES = 32 * 1024 * 1024
_ASSET_SET_MAX_BYTES = 64 * 1024 * 1024
_ASSET_SET_MAX_FILES = 10_000
_READ_CHUNK_BYTES = 1024 * 1024

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_IDENTITY = re.compile(r"^sha256-[0-9a-f]{64}$")
_VERSION = re.compile(r"^[0-9A-Za-z][0-9A-Za-z._+-]{0,127}$")
_PATH_PART = re.compile(r"^[0-9A-Za-z][0-9A-Za-z._@+=~-]{0,255}$")


def _invalid_manifest(
    message: str = "Frontend asset manifest is invalid.",
) -> DeploymentError:
    return fail(
        "FRONTEND_ASSET_MANIFEST_INVALID",
        message,
        component="frontend",
    )


def _invalid_package() -> DeploymentError:
    return fail(
        "FRONTEND_ASSET_PACKAGE_INVALID",
        "Packaged frontend assets are invalid.",
        component="frontend",
    )


def _integrity_failed() -> DeploymentError:
    return fail(
        "FRONTEND_ASSET_INTEGRITY_FAILED",
        "Packaged frontend asset integrity verification failed.",
        component="frontend",
    )


def _json_pairs(pairs: list[tuple[str, object]]) -> dict[str, object]:
    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError("duplicate JSON key")
        value[key] = item
    return value


def _reject_json_constant(_value: str) -> object:
    raise ValueError("non-finite JSON number")


def _load_json_document(raw: bytes) -> object:
    return json.loads(
        raw.decode("utf-8"),
        object_pairs_hook=_json_pairs,
        parse_constant=_reject_json_constant,
    )


def canonical_json_sha256(raw: bytes) -> str:
    """Return the semantic SHA-256 of one strict JSON document."""
    if not isinstance(raw, bytes) or not 0 < len(raw) <= _METADATA_MAX_BYTES:
        raise _invalid_manifest()
    try:
        document = _load_json_document(raw)
        rendered = canonical_json_bytes(document)
    except (UnicodeDecodeError, ValueError, TypeError):
        raise _invalid_manifest() from None
    return hashlib.sha256(rendered).hexdigest()


def _metadata_sha256(raw: bytes) -> str:
    if not isinstance(raw, bytes) or not 0 < len(raw) <= _METADATA_MAX_BYTES:
        raise _invalid_manifest()
    return hashlib.sha256(raw).hexdigest()


def _object(raw: object) -> dict[str, Any]:
    if not isinstance(raw, Mapping) or any(not isinstance(key, str) for key in raw):
        raise _invalid_manifest()
    return dict(raw)


def _exact_keys(value: Mapping[str, object], expected: set[str]) -> None:
    if set(value) != expected:
        raise _invalid_manifest()


def _sha256(raw: object) -> str:
    if not isinstance(raw, str) or _SHA256.fullmatch(raw) is None:
        raise _invalid_manifest()
    return raw


def _identity(raw: object) -> str:
    if not isinstance(raw, str) or _IDENTITY.fullmatch(raw) is None:
        raise _invalid_manifest()
    return raw


def _version(raw: object) -> str:
    if not isinstance(raw, str) or _VERSION.fullmatch(raw) is None:
        raise _invalid_manifest()
    return raw


def _asset_path(raw: object) -> str:
    if not isinstance(raw, str) or not raw or len(raw) > 512:
        raise _invalid_manifest()
    if any(ord(character) < 32 or ord(character) == 127 for character in raw):
        raise _invalid_manifest()
    if "\\" in raw or raw.startswith("/") or raw.endswith("/") or "//" in raw:
        raise _invalid_manifest()
    path = PurePosixPath(raw)
    parts = path.parts
    if (
        not parts
        or path.is_absolute()
        or any(part in {"", ".", ".."} for part in parts)
        or any(_PATH_PART.fullmatch(part) is None for part in parts)
    ):
        raise _invalid_manifest()
    return raw


@dataclass(frozen=True, order=True)
class AssetRecord:
    """One immutable file in the frontend asset closure."""

    path: str
    size_bytes: int
    sha256: str

    @classmethod
    def from_dict(cls, raw: object) -> "AssetRecord":
        value = _object(raw)
        _exact_keys(value, {"path", "size_bytes", "sha256"})
        size = value["size_bytes"]
        if (
            not isinstance(size, int)
            or isinstance(size, bool)
            or not 0 <= size <= _FILE_MAX_BYTES
        ):
            raise _invalid_manifest()
        return cls(
            path=_asset_path(value["path"]),
            size_bytes=size,
            sha256=_sha256(value["sha256"]),
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "path": self.path,
            "size_bytes": self.size_bytes,
            "sha256": self.sha256,
        }


@dataclass(frozen=True)
class FrontendAssetManifest:
    """Canonical build identity and API compatibility for frontend assets."""

    identity: str
    frontend_version: str
    package_lock_sha256: str
    api_contract_sha256: str
    entrypoint: str
    files: tuple[AssetRecord, ...]

    @classmethod
    def create(
        cls,
        *,
        frontend_version: str,
        package_lock_sha256: str,
        api_contract_sha256: str,
        entrypoint: str,
        files: Sequence[AssetRecord],
    ) -> "FrontendAssetManifest":
        value: dict[str, object] = {
            "schema_version": FRONTEND_ASSET_MANIFEST_SCHEMA,
            "frontend_version": frontend_version,
            "package_lock_sha256": package_lock_sha256,
            "api_contract_sha256": api_contract_sha256,
            "entrypoint": entrypoint,
            "files": [item.to_dict() for item in sorted(files)],
        }
        value["identity"] = canonical_identity(
            value,
            scheme=FRONTEND_ASSET_IDENTITY_SCHEME,
        )
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "FrontendAssetManifest":
        value = _object(raw)
        _exact_keys(
            value,
            {
                "schema_version",
                "identity",
                "frontend_version",
                "package_lock_sha256",
                "api_contract_sha256",
                "entrypoint",
                "files",
            },
        )
        if value["schema_version"] != FRONTEND_ASSET_MANIFEST_SCHEMA:
            raise _invalid_manifest()
        raw_files = value["files"]
        if not isinstance(raw_files, Sequence) or isinstance(raw_files, (str, bytes)):
            raise _invalid_manifest()
        if not 0 < len(raw_files) <= _ASSET_SET_MAX_FILES:
            raise _invalid_manifest()
        files = tuple(AssetRecord.from_dict(item) for item in raw_files)
        if tuple(sorted(files)) != files or len({item.path for item in files}) != len(
            files
        ):
            raise _invalid_manifest()
        if sum(item.size_bytes for item in files) > _ASSET_SET_MAX_BYTES:
            raise _invalid_manifest()
        entrypoint = _asset_path(value["entrypoint"])
        if entrypoint != ENTRYPOINT or entrypoint not in {item.path for item in files}:
            raise _invalid_manifest()

        observed_identity = _identity(value["identity"])
        expected_identity = canonical_identity(
            without_key(value, "identity"),
            scheme=FRONTEND_ASSET_IDENTITY_SCHEME,
        )
        if observed_identity != expected_identity:
            raise fail(
                "FRONTEND_ASSET_IDENTITY_INVALID",
                "Frontend asset manifest identity is invalid.",
                component="frontend",
            )
        return cls(
            identity=observed_identity,
            frontend_version=_version(value["frontend_version"]),
            package_lock_sha256=_sha256(value["package_lock_sha256"]),
            api_contract_sha256=_sha256(value["api_contract_sha256"]),
            entrypoint=entrypoint,
            files=files,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": FRONTEND_ASSET_MANIFEST_SCHEMA,
            "identity": self.identity,
            "frontend_version": self.frontend_version,
            "package_lock_sha256": self.package_lock_sha256,
            "api_contract_sha256": self.api_contract_sha256,
            "entrypoint": self.entrypoint,
            "files": [item.to_dict() for item in self.files],
        }

    def to_bytes(self) -> bytes:
        return canonical_json_bytes(self.to_dict())


@dataclass(frozen=True)
class VerifiedFrontendAssets:
    """Manifest-bound in-memory bytes safe to serve without later path access."""

    manifest: FrontendAssetManifest
    content: Mapping[str, bytes]

    @classmethod
    def create(
        cls,
        manifest: FrontendAssetManifest,
        content: Mapping[str, bytes],
    ) -> "VerifiedFrontendAssets":
        if any(
            not isinstance(path, str) or not isinstance(value, bytes)
            for path, value in content.items()
        ):
            raise _integrity_failed()
        copied = {path: bytes(value) for path, value in content.items()}
        records = {item.path: item for item in manifest.files}
        if set(copied) != set(records):
            raise _integrity_failed()
        for path, value in copied.items():
            record = records[path]
            if (
                len(value) != record.size_bytes
                or hashlib.sha256(value).hexdigest() != record.sha256
            ):
                raise _integrity_failed()
        return cls(manifest=manifest, content=MappingProxyType(copied))


def parse_manifest_bytes(raw: bytes) -> FrontendAssetManifest:
    """Parse a byte-for-byte canonical frontend manifest."""
    if not isinstance(raw, bytes) or not 0 < len(raw) <= _MANIFEST_MAX_BYTES:
        raise _invalid_manifest()
    try:
        document = _load_json_document(raw)
        if canonical_json_bytes(document) != raw:
            raise _invalid_manifest()
        return FrontendAssetManifest.from_dict(document)
    except DeploymentError:
        raise
    except (UnicodeDecodeError, ValueError, TypeError):
        raise _invalid_manifest() from None


def _read_regular_file(path: Path, *, maximum_bytes: int) -> bytes:
    flags = os.O_RDONLY
    flags |= getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    descriptor: int | None = None
    try:
        before = os.lstat(path)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_size > maximum_bytes
        ):
            raise _invalid_package()
        descriptor = os.open(path, flags)
        opened = os.fstat(descriptor)
        if (
            not stat.S_ISREG(opened.st_mode)
            or opened.st_nlink != 1
            or opened.st_size > maximum_bytes
            or (opened.st_dev, opened.st_ino) != (before.st_dev, before.st_ino)
        ):
            raise _invalid_package()
        chunks: list[bytes] = []
        observed = 0
        while True:
            chunk = os.read(
                descriptor, min(_READ_CHUNK_BYTES, maximum_bytes + 1 - observed)
            )
            if not chunk:
                break
            observed += len(chunk)
            if observed > maximum_bytes:
                raise _invalid_package()
            chunks.append(chunk)
        after = os.lstat(path)
        if observed != opened.st_size or (after.st_dev, after.st_ino) != (
            opened.st_dev,
            opened.st_ino,
        ):
            raise _invalid_package()
        return b"".join(chunks)
    except DeploymentError:
        raise
    except OSError:
        raise _invalid_package() from None
    finally:
        if descriptor is not None:
            os.close(descriptor)


def _directory_files(root: Path) -> tuple[Path, ...]:
    def fail_on_walk_error(error: OSError) -> None:
        raise error

    try:
        root_stat = os.lstat(root)
        if not stat.S_ISDIR(root_stat.st_mode):
            raise _invalid_package()
        observed: list[Path] = []
        for directory, directory_names, file_names in os.walk(
            root,
            topdown=True,
            followlinks=False,
            onerror=fail_on_walk_error,
        ):
            directory_path = Path(directory)
            for name in sorted(directory_names):
                child = directory_path / name
                if not stat.S_ISDIR(os.lstat(child).st_mode):
                    raise _invalid_package()
            for name in sorted(file_names):
                child = directory_path / name
                if not stat.S_ISREG(os.lstat(child).st_mode):
                    raise _invalid_package()
                _asset_path(child.relative_to(root).as_posix())
                observed.append(child)
                if len(observed) > _ASSET_SET_MAX_FILES:
                    raise _invalid_package()
        return tuple(
            sorted(observed, key=lambda item: item.relative_to(root).as_posix())
        )
    except DeploymentError:
        raise
    except (OSError, ValueError):
        raise _invalid_package() from None


def _load_directory_content(root: Path) -> dict[str, bytes]:
    files = _directory_files(root)
    if not 0 < len(files) <= _ASSET_SET_MAX_FILES:
        raise _invalid_package()
    content: dict[str, bytes] = {}
    total = 0
    for path in files:
        relative = path.relative_to(root).as_posix()
        value = _read_regular_file(path, maximum_bytes=_FILE_MAX_BYTES)
        total += len(value)
        if total > _ASSET_SET_MAX_BYTES:
            raise _invalid_package()
        content[relative] = value
    return content


def build_frontend_assets(
    asset_root: Path,
    *,
    frontend_version: str,
    package_lock_bytes: bytes,
    openapi_bytes: bytes,
) -> VerifiedFrontendAssets:
    """Create canonical manifest-bound bytes from one completed Vite build."""
    content = _load_directory_content(asset_root)
    files = tuple(
        AssetRecord(
            path=path,
            size_bytes=len(value),
            sha256=hashlib.sha256(value).hexdigest(),
        )
        for path, value in sorted(content.items())
    )
    manifest = FrontendAssetManifest.create(
        frontend_version=_version(frontend_version),
        package_lock_sha256=_metadata_sha256(package_lock_bytes),
        api_contract_sha256=canonical_json_sha256(openapi_bytes),
        entrypoint=ENTRYPOINT,
        files=files,
    )
    return VerifiedFrontendAssets.create(manifest, content)


def verify_asset_directory(
    manifest: FrontendAssetManifest,
    asset_root: Path,
) -> VerifiedFrontendAssets:
    """Verify an exact on-disk asset closure and retain the verified bytes."""
    try:
        content = _load_directory_content(asset_root)
        return VerifiedFrontendAssets.create(manifest, content)
    except DeploymentError as error:
        if error.issue.code == "FRONTEND_ASSET_INTEGRITY_FAILED":
            raise
        raise _integrity_failed() from None


def load_packaged_frontend_assets(
    package: str = "encode_pipeline.frontend_assets",
) -> VerifiedFrontendAssets:
    """Load and verify the installed package's production frontend closure."""
    try:
        package_root = Path(os.fspath(resources.files(package)))
    except (ModuleNotFoundError, TypeError, ValueError):
        raise _invalid_package() from None
    manifest_raw = _read_regular_file(
        package_root / MANIFEST_FILENAME,
        maximum_bytes=_MANIFEST_MAX_BYTES,
    )
    manifest = parse_manifest_bytes(manifest_raw)
    return verify_asset_directory(manifest, package_root / STATIC_DIRECTORY)


def verify_frontend_api_contract(
    manifest: FrontendAssetManifest,
    schema: object,
) -> None:
    """Fail closed unless the active API exactly matches the built client contract."""
    try:
        digest = hashlib.sha256(canonical_json_bytes(schema)).hexdigest()
    except (TypeError, ValueError):
        raise fail(
            "FRONTEND_API_CONTRACT_INVALID",
            "Active API compatibility could not be verified.",
            component="frontend",
        ) from None
    if digest != manifest.api_contract_sha256:
        raise fail(
            "FRONTEND_API_INCOMPATIBLE",
            "The packaged frontend is incompatible with the active API.",
            component="frontend",
        )


__all__ = [
    "ENTRYPOINT",
    "FRONTEND_ASSET_IDENTITY_SCHEME",
    "FRONTEND_ASSET_MANIFEST_SCHEMA",
    "MANIFEST_FILENAME",
    "STATIC_DIRECTORY",
    "AssetRecord",
    "FrontendAssetManifest",
    "VerifiedFrontendAssets",
    "build_frontend_assets",
    "canonical_json_sha256",
    "load_packaged_frontend_assets",
    "parse_manifest_bytes",
    "verify_asset_directory",
    "verify_frontend_api_contract",
]
