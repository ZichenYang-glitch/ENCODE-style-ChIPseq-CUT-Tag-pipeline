"""Canonical offline bundle admission and immutable staging."""

from __future__ import annotations

from collections.abc import Callable
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import stat
import tarfile
import uuid

from encode_pipeline.deployment.canonical import canonical_json_bytes
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.filesystem import (
    create_directory,
    fsync_directory,
    hash_regular_file,
    read_regular_file,
)
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import BundleManifest, FileRecord


MANIFEST_FILENAME = "manifest.json"
MAX_MANIFEST_BYTES = 2 * 1024 * 1024
MAX_BUNDLE_BYTES = 200 * 1024**3
MAX_BUNDLE_FILES = 100_000
MIN_FREE_SPACE_RESERVE = 2 * 1024**3
FaultInjector = Callable[[str], None]
_CONTENT_IDENTITY = re.compile(r"^sha256-[0-9a-f]{64}$")


class BundleStore:
    """Admit one plain canonical tar without executing any bundled content."""

    def __init__(self, layout: DeploymentLayout) -> None:
        self.layout = layout

    def inspect(
        self,
        bundle_path: Path,
        *,
        expected_owner_uid: int | None = None,
        expected_owner_gid: int | None = None,
    ) -> BundleManifest:
        with self._open_bundle(
            bundle_path,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        ) as archive:
            manifest, _raw = self._read_manifest(archive)
            self._verify_member_inventory(archive, manifest)
            return manifest

    def stage(
        self,
        bundle_path: Path,
        *,
        expected_owner_uid: int | None = None,
        expected_owner_gid: int | None = None,
        installed_owner_uid: int | None = None,
        installed_owner_gid: int | None = None,
        fault: FaultInjector | None = None,
    ) -> BundleManifest:
        with self._open_bundle(
            bundle_path,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        ) as archive:
            manifest, raw_manifest = self._read_manifest(archive)
            members = self._verify_member_inventory(archive, manifest)
            store = self._prepare_store(
                manifest.component,
                expected_owner_uid=installed_owner_uid,
                expected_owner_gid=installed_owner_gid,
            )
            self._preflight_capacity(store, manifest)
            destination = store / manifest.identity
            if destination.exists() or destination.is_symlink():
                existing = self.verify_installed(
                    manifest.component,
                    manifest.identity,
                    expected_owner_uid=installed_owner_uid,
                    expected_owner_gid=installed_owner_gid,
                )
                if existing != manifest:
                    raise fail(
                        "DEPLOYMENT_IDENTITY_CONFLICT",
                        "Deployment identity conflicts with installed content.",
                        component=manifest.component,
                    )
                return existing

            partial = store / f".partial-{manifest.identity}-{uuid.uuid4().hex}"
            try:
                partial.mkdir(mode=0o700)
                self._write_file(partial / MANIFEST_FILENAME, raw_manifest, 0o444)
                for member, record in zip(members, manifest.files, strict=True):
                    target = partial.joinpath(*Path(record.path).parts)
                    self._create_relative_parents(partial, target.parent)
                    extracted = archive.extractfile(member)
                    if extracted is None:
                        raise fail(
                            "DEPLOYMENT_BUNDLE_INVALID", "Deployment bundle is invalid."
                        )
                    self._copy_member(extracted, target, record)
                self._freeze_directories(partial)
                if fault is not None:
                    fault("payload-verified")
                os.replace(partial, destination)
                fsync_directory(store)
            except DeploymentError:
                self._discard_partial(partial)
                raise
            except Exception:
                self._discard_partial(partial)
                raise fail(
                    "DEPLOYMENT_STAGE_FAILED",
                    "Deployment bundle could not be staged.",
                    component=manifest.component,
                    recoverable=True,
                ) from None
            if fault is not None:
                fault("release-committed")
            return self.verify_installed(
                manifest.component,
                manifest.identity,
                expected_owner_uid=installed_owner_uid,
                expected_owner_gid=installed_owner_gid,
            )

    def verify_installed(
        self,
        component: str,
        identity: str,
        *,
        expected_owner_uid: int | None = None,
        expected_owner_gid: int | None = None,
    ) -> BundleManifest:
        root, manifest = self._read_installed_manifest(
            component,
            identity,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )
        expected_paths = {MANIFEST_FILENAME}
        for record in manifest.files:
            expected_paths.add(record.path)
            path = root.joinpath(*Path(record.path).parts)
            digest, observed = hash_regular_file(
                path,
                expected_size=record.size_bytes,
                code="DEPLOYMENT_RELEASE_INVALID",
            )
            if (
                digest != record.sha256
                or stat.S_IMODE(observed.st_mode) != record.mode
                or (
                    expected_owner_uid is not None
                    and observed.st_uid != expected_owner_uid
                )
                or (
                    expected_owner_gid is not None
                    and observed.st_gid != expected_owner_gid
                )
            ):
                raise fail(
                    "DEPLOYMENT_RELEASE_INVALID",
                    "Deployment content is invalid.",
                    component=component,
                )
        observed_paths: set[str] = set()
        for path in root.rglob("*"):
            relative = path.relative_to(root).as_posix()
            observed = path.lstat()
            if stat.S_ISLNK(observed.st_mode):
                raise fail(
                    "DEPLOYMENT_RELEASE_INVALID",
                    "Deployment content is invalid.",
                    component=component,
                )
            if stat.S_ISDIR(observed.st_mode):
                if (
                    stat.S_IMODE(observed.st_mode) != 0o555
                    or (
                        expected_owner_uid is not None
                        and observed.st_uid != expected_owner_uid
                    )
                    or (
                        expected_owner_gid is not None
                        and observed.st_gid != expected_owner_gid
                    )
                ):
                    raise fail(
                        "DEPLOYMENT_RELEASE_INVALID",
                        "Deployment content is invalid.",
                        component=component,
                    )
                continue
            if not stat.S_ISREG(observed.st_mode):
                raise fail(
                    "DEPLOYMENT_RELEASE_INVALID",
                    "Deployment content is invalid.",
                    component=component,
                )
            observed_paths.add(relative)
        if observed_paths != expected_paths:
            raise fail(
                "DEPLOYMENT_RELEASE_INVALID",
                "Deployment content is invalid.",
                component=component,
            )
        return manifest

    def read_installed_manifest(
        self,
        component: str,
        identity: str,
        *,
        expected_owner_uid: int | None = None,
        expected_owner_gid: int | None = None,
    ) -> BundleManifest:
        """Read immutable slot metadata without hashing the payload closure."""
        _root, manifest = self._read_installed_manifest(
            component,
            identity,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )
        return manifest

    def _read_installed_manifest(
        self,
        component: str,
        identity: str,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> tuple[Path, BundleManifest]:
        if _CONTENT_IDENTITY.fullmatch(identity) is None:
            raise fail(
                "DEPLOYMENT_RELEASE_INVALID",
                "Deployment content is invalid.",
                component=component
                if component in {"platform", "encode-runtime", "bulk-rnaseq-runtime"}
                else None,
            )
        store = self._require_store(
            component,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )
        root = store / identity
        try:
            root_stat = root.lstat()
        except OSError:
            raise fail(
                "DEPLOYMENT_RELEASE_UNAVAILABLE",
                "Deployment content is unavailable.",
                component=component,
            ) from None
        if (
            not stat.S_ISDIR(root_stat.st_mode)
            or stat.S_ISLNK(root_stat.st_mode)
            or stat.S_IMODE(root_stat.st_mode) != 0o555
            or (
                expected_owner_uid is not None
                and root_stat.st_uid != expected_owner_uid
            )
            or (
                expected_owner_gid is not None
                and root_stat.st_gid != expected_owner_gid
            )
        ):
            raise fail(
                "DEPLOYMENT_RELEASE_INVALID",
                "Deployment content is invalid.",
                component=component,
            )
        raw_manifest, manifest_stat = read_regular_file(
            root / MANIFEST_FILENAME,
            max_bytes=MAX_MANIFEST_BYTES,
            code="DEPLOYMENT_RELEASE_INVALID",
        )
        if (
            stat.S_IMODE(manifest_stat.st_mode) != 0o444
            or (
                expected_owner_uid is not None
                and manifest_stat.st_uid != expected_owner_uid
            )
            or (
                expected_owner_gid is not None
                and manifest_stat.st_gid != expected_owner_gid
            )
        ):
            raise fail(
                "DEPLOYMENT_RELEASE_INVALID",
                "Deployment content is invalid.",
                component=component,
            )
        manifest = _parse_manifest(raw_manifest)
        if manifest.component != component or manifest.identity != identity:
            raise fail(
                "DEPLOYMENT_RELEASE_INVALID",
                "Deployment content is invalid.",
                component=component,
            )
        return root, manifest

    def _prepare_store(
        self,
        component: str,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> Path:
        store = self.layout.component_store(component)
        if (
            not self.layout.immutable_root.exists()
            and not self.layout.immutable_root.is_symlink()
        ):
            create_directory(self.layout.immutable_root)
        current = self.layout.immutable_root
        self._require_store_directory(
            current,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
            code="DEPLOYMENT_PATH_UNSAFE",
        )
        for part in store.relative_to(self.layout.immutable_root).parts:
            current = current / part
            if not current.exists() and not current.is_symlink():
                try:
                    current.mkdir(mode=0o755)
                except OSError:
                    raise fail(
                        "DEPLOYMENT_STORAGE_UNAVAILABLE",
                        "Deployment storage is unavailable.",
                    ) from None
            self._require_store_directory(
                current,
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
                code="DEPLOYMENT_PATH_UNSAFE",
            )
        return store

    def _require_store(
        self,
        component: str,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> Path:
        store = self.layout.component_store(component)
        current = self.layout.immutable_root
        self._require_store_directory(
            current,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
            code="DEPLOYMENT_RELEASE_INVALID",
        )
        for part in store.relative_to(self.layout.immutable_root).parts:
            current = current / part
            self._require_store_directory(
                current,
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
                code="DEPLOYMENT_RELEASE_INVALID",
            )
        return store

    @staticmethod
    def _require_store_directory(
        path: Path,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
        code: str,
    ) -> None:
        try:
            observed = path.lstat()
        except OSError:
            raise fail(
                "DEPLOYMENT_STORAGE_UNAVAILABLE",
                "Deployment storage is unavailable.",
            ) from None
        if (
            not stat.S_ISDIR(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or stat.S_IMODE(observed.st_mode) & 0o022
            or (
                expected_owner_uid is not None and observed.st_uid != expected_owner_uid
            )
            or (
                expected_owner_gid is not None and observed.st_gid != expected_owner_gid
            )
        ):
            raise fail(
                code,
                (
                    "Deployment content is invalid."
                    if code == "DEPLOYMENT_RELEASE_INVALID"
                    else "Deployment storage boundary is unsafe."
                ),
            )

    def partial_staging(self) -> tuple[tuple[str, str], ...]:
        values: list[tuple[str, str]] = []
        from encode_pipeline.deployment.models import COMPONENTS

        for component in COMPONENTS:
            store = self.layout.component_store(component)
            if not store.exists() or store.is_symlink():
                continue
            for path in sorted(store.iterdir()):
                if path.name.startswith(".partial-"):
                    values.append((component, path.name))
        return tuple(values)

    def installed_identities(self, component: str) -> tuple[str, ...]:
        store = self.layout.component_store(component)
        if not store.exists() or store.is_symlink():
            return ()
        identities: list[str] = []
        for path in sorted(store.iterdir()):
            if path.name.startswith("."):
                continue
            try:
                observed = path.lstat()
            except OSError:
                continue
            if stat.S_ISDIR(observed.st_mode) and not stat.S_ISLNK(observed.st_mode):
                identities.append(path.name)
        return tuple(identities)

    @contextmanager
    def _open_bundle(
        self,
        bundle_path: Path,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ):
        if not isinstance(bundle_path, Path) or not bundle_path.is_absolute():
            raise fail(
                "DEPLOYMENT_BUNDLE_PATH_INVALID", "Deployment bundle path is invalid."
            )
        flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
        try:
            descriptor = os.open(bundle_path, flags)
            observed = os.fstat(descriptor)
        except OSError:
            raise fail(
                "DEPLOYMENT_BUNDLE_UNAVAILABLE", "Deployment bundle is unavailable."
            ) from None
        if (
            not stat.S_ISREG(observed.st_mode)
            or observed.st_nlink != 1
            or not 0 < observed.st_size <= MAX_BUNDLE_BYTES
            or observed.st_size % tarfile.RECORDSIZE != 0
            or stat.S_IMODE(observed.st_mode) & 0o022
            or (
                expected_owner_uid is not None and observed.st_uid != expected_owner_uid
            )
            or (
                expected_owner_gid is not None and observed.st_gid != expected_owner_gid
            )
        ):
            os.close(descriptor)
            raise fail("DEPLOYMENT_BUNDLE_INVALID", "Deployment bundle is invalid.")
        before = (
            observed.st_dev,
            observed.st_ino,
            observed.st_size,
            observed.st_mtime_ns,
            observed.st_ctime_ns,
        )
        file_object = os.fdopen(descriptor, "rb", closefd=True)
        try:
            archive = tarfile.open(fileobj=file_object, mode="r:")
        except (OSError, tarfile.TarError):
            file_object.close()
            raise fail(
                "DEPLOYMENT_BUNDLE_INVALID", "Deployment bundle is invalid."
            ) from None
        try:
            yield archive
            after_stat = os.fstat(file_object.fileno())
            after = (
                after_stat.st_dev,
                after_stat.st_ino,
                after_stat.st_size,
                after_stat.st_mtime_ns,
                after_stat.st_ctime_ns,
            )
            if after != before:
                raise fail(
                    "DEPLOYMENT_BUNDLE_RACE",
                    "Deployment bundle changed during verification.",
                    recoverable=True,
                )
        finally:
            archive.close()
            file_object.close()

    @staticmethod
    def _read_manifest(archive: tarfile.TarFile) -> tuple[BundleManifest, bytes]:
        try:
            member = archive.next()
        except tarfile.TarError:
            raise fail(
                "DEPLOYMENT_BUNDLE_INVALID", "Deployment bundle is invalid."
            ) from None
        if (
            member is None
            or member.name != MANIFEST_FILENAME
            or member.type != tarfile.REGTYPE
            or not 0 < member.size <= MAX_MANIFEST_BYTES
            or member.mode != 0o444
            or member.uid != 0
            or member.gid != 0
            or member.uname != ""
            or member.gname != ""
            or member.mtime != 0
            or member.linkname != ""
            or member.pax_headers
            or member.devmajor != 0
            or member.devminor != 0
        ):
            raise fail("DEPLOYMENT_BUNDLE_INVALID", "Deployment bundle is invalid.")
        handle = archive.extractfile(member)
        if handle is None:
            raise fail("DEPLOYMENT_BUNDLE_INVALID", "Deployment bundle is invalid.")
        raw = handle.read(MAX_MANIFEST_BYTES + 1)
        return _parse_manifest(raw), raw

    @staticmethod
    def _verify_member_inventory(
        archive: tarfile.TarFile,
        manifest: BundleManifest,
    ) -> tuple[tarfile.TarInfo, ...]:
        members: list[tarfile.TarInfo] = []
        while True:
            try:
                member = archive.next()
            except tarfile.TarError:
                raise fail(
                    "DEPLOYMENT_BUNDLE_INVALID", "Deployment bundle is invalid."
                ) from None
            if member is None:
                break
            members.append(member)
            if len(members) > MAX_BUNDLE_FILES:
                raise fail("DEPLOYMENT_BUNDLE_INVALID", "Deployment bundle is invalid.")
        if len(members) != len(manifest.files):
            raise fail("DEPLOYMENT_BUNDLE_INVALID", "Deployment bundle is invalid.")
        for member, record in zip(members, manifest.files, strict=True):
            if (
                member.name != record.path
                or member.type != tarfile.REGTYPE
                or member.size != record.size_bytes
                or member.mode != record.mode
                or member.uid != 0
                or member.gid != 0
                or member.uname != ""
                or member.gname != ""
                or member.mtime != 0
                or member.linkname != ""
                or member.pax_headers
                or member.devmajor != 0
                or member.devminor != 0
            ):
                raise fail("DEPLOYMENT_BUNDLE_INVALID", "Deployment bundle is invalid.")
        expected_size = (
            (archive.offset + 1024 + tarfile.RECORDSIZE - 1)
            // tarfile.RECORDSIZE
            * tarfile.RECORDSIZE
        )
        file_object = archive.fileobj
        try:
            observed_size = os.fstat(file_object.fileno()).st_size
            remaining = expected_size - file_object.tell()
            if remaining < 0:
                raise ValueError
            trailing = file_object.read(remaining)
        except (OSError, ValueError):
            raise fail(
                "DEPLOYMENT_BUNDLE_INVALID", "Deployment bundle is invalid."
            ) from None
        if (
            observed_size != expected_size
            or len(trailing) != remaining
            or any(trailing)
        ):
            raise fail("DEPLOYMENT_BUNDLE_INVALID", "Deployment bundle is invalid.")
        return tuple(members)

    @staticmethod
    def _create_relative_parents(root: Path, parent: Path) -> None:
        relative = parent.relative_to(root)
        current = root
        for part in relative.parts:
            current = current / part
            try:
                observed = current.lstat()
            except FileNotFoundError:
                current.mkdir(mode=0o700)
                continue
            if not stat.S_ISDIR(observed.st_mode) or stat.S_ISLNK(observed.st_mode):
                raise fail("DEPLOYMENT_BUNDLE_INVALID", "Deployment bundle is invalid.")

    @staticmethod
    def _copy_member(source, target: Path, record: FileRecord) -> None:
        flags = (
            os.O_WRONLY
            | os.O_CREAT
            | os.O_EXCL
            | getattr(os, "O_CLOEXEC", 0)
            | getattr(os, "O_NOFOLLOW", 0)
        )
        digest = hashlib.sha256()
        observed_size = 0
        try:
            descriptor = os.open(target, flags, 0o600)
            try:
                while True:
                    chunk = source.read(1024 * 1024)
                    if not chunk:
                        break
                    observed_size += len(chunk)
                    if observed_size > record.size_bytes:
                        raise fail(
                            "DEPLOYMENT_BUNDLE_INVALID", "Deployment bundle is invalid."
                        )
                    digest.update(chunk)
                    written = 0
                    while written < len(chunk):
                        written += os.write(descriptor, chunk[written:])
                if (
                    observed_size != record.size_bytes
                    or digest.hexdigest() != record.sha256
                ):
                    raise fail(
                        "DEPLOYMENT_BUNDLE_INVALID", "Deployment bundle is invalid."
                    )
                os.fchmod(descriptor, record.mode)
                os.fsync(descriptor)
            finally:
                os.close(descriptor)
        except DeploymentError:
            raise
        except OSError:
            raise fail(
                "DEPLOYMENT_STAGE_FAILED", "Deployment bundle could not be staged."
            ) from None

    @staticmethod
    def _write_file(path: Path, content: bytes, mode: int) -> None:
        flags = (
            os.O_WRONLY
            | os.O_CREAT
            | os.O_EXCL
            | getattr(os, "O_CLOEXEC", 0)
            | getattr(os, "O_NOFOLLOW", 0)
        )
        try:
            descriptor = os.open(path, flags, 0o600)
            try:
                written = 0
                while written < len(content):
                    written += os.write(descriptor, content[written:])
                os.fchmod(descriptor, mode)
                os.fsync(descriptor)
            finally:
                os.close(descriptor)
        except OSError:
            raise fail(
                "DEPLOYMENT_STAGE_FAILED", "Deployment bundle could not be staged."
            ) from None

    @staticmethod
    def _freeze_directories(root: Path) -> None:
        directories = [path for path in root.rglob("*") if path.is_dir()]
        for path in sorted(directories, key=lambda item: len(item.parts), reverse=True):
            os.chmod(path, 0o555)
            fsync_directory(path)
        os.chmod(root, 0o555)
        fsync_directory(root)

    @staticmethod
    def _preflight_capacity(store: Path, manifest: BundleManifest) -> None:
        declared = sum(item.size_bytes for item in manifest.files) + MAX_MANIFEST_BYTES
        if declared > MAX_BUNDLE_BYTES:
            raise fail("DEPLOYMENT_BUNDLE_INVALID", "Deployment bundle is invalid.")
        try:
            available = shutil.disk_usage(store).free
        except OSError:
            raise fail(
                "DEPLOYMENT_STORAGE_UNAVAILABLE",
                "Deployment storage is unavailable.",
                recoverable=True,
            ) from None
        reserve = max(MIN_FREE_SPACE_RESERVE, declared // 20)
        if available < declared + reserve:
            raise fail(
                "DEPLOYMENT_CAPACITY_INSUFFICIENT",
                "Deployment storage capacity is insufficient.",
                recoverable=True,
            )

    @staticmethod
    def _discard_partial(path: Path) -> None:
        if (
            path.name.startswith(".partial-")
            and path.exists()
            and not path.is_symlink()
        ):
            shutil.rmtree(path)


def _parse_manifest(raw: bytes) -> BundleManifest:
    if len(raw) > MAX_MANIFEST_BYTES:
        raise fail("DEPLOYMENT_MANIFEST_INVALID", "Deployment document is invalid.")
    try:
        document = json.loads(raw, object_pairs_hook=_unique_object)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError):
        raise fail(
            "DEPLOYMENT_MANIFEST_INVALID", "Deployment document is invalid."
        ) from None
    manifest = BundleManifest.from_dict(document)
    if raw != canonical_json_bytes(manifest.to_dict()):
        raise fail("DEPLOYMENT_MANIFEST_INVALID", "Deployment document is invalid.")
    return manifest


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError("duplicate JSON key")
        value[key] = item
    return value


__all__ = ["BundleStore", "FaultInjector", "MANIFEST_FILENAME"]
