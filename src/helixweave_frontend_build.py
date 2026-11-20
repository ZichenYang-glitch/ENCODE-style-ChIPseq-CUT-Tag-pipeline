"""Stage the current frontend closure into each wheel build.

setuptools copies package data into ``build/lib`` and packs that tree.
``copy_file`` skips a destination whose modification time is not older, so a
later frontend build can leave the previous ``index.html`` and
``asset-manifest.json`` bytes in place. Removed chunk names have the same
problem when they remain beside the new files.

Each wheel build therefore copies the current static tree and manifest into
a new empty directory, then installs that copy into this build. Every
directory from ``build/lib`` through the frontend package must be a real
directory; a symlink anywhere on that path keeps the outside tree unchanged.
"""

from __future__ import annotations

import os
from pathlib import Path
import stat
import tempfile

from setuptools.command.build_py import build_py as _BuildPy
from setuptools.errors import BaseError


_PACKAGE = "encode_pipeline.frontend_assets"
_STATIC = "static"
_MANIFEST = "asset-manifest.json"
_STAGING_FAILED = "Cannot stage the current frontend assets for this wheel build."


class build_py(_BuildPy):
    """Package-data build that stages the current frontend closure."""

    def run(self) -> None:
        super().run()
        packages = set(self.packages or ())
        if (
            getattr(self, "editable_mode", False)
            or not self.build_lib
            or _PACKAGE not in packages
        ):
            return
        stage_frontend_assets(
            source_package=Path(self.get_package_dir(_PACKAGE)),
            build_lib=Path(self.build_lib),
        )


def stage_frontend_assets(*, source_package: Path, build_lib: Path) -> None:
    """Copy the current frontend assets and manifest into this build."""

    if _trees_overlap(source_package, build_lib):
        raise BaseError(_STAGING_FAILED)
    package_dir = _trusted_package_directory(build_lib)
    source_static = source_package / _STATIC
    source_manifest = source_package / _MANIFEST
    if not _real_directory(source_static) or not _regular_file(source_manifest):
        raise BaseError(_STAGING_FAILED)

    staging = Path(tempfile.mkdtemp(prefix="helixweave-frontend-"))
    try:
        staged_static = staging / _STATIC
        staged_manifest = staging / _MANIFEST
        _copy_regular_tree(source_static, staged_static)
        _copy_regular_file(source_manifest, staged_manifest)
        if _regular_files(source_static) != _regular_files(staged_static):
            raise BaseError(_STAGING_FAILED)
        _install_tree(package_dir, _STATIC, staged_static)
        _install_file(package_dir, _MANIFEST, staged_manifest)
    finally:
        _remove_owned_tree(staging)


def _trusted_package_directory(build_lib: Path) -> Path:
    current = build_lib
    if not _real_directory(current):
        raise BaseError(_STAGING_FAILED)
    for name in _PACKAGE.split("."):
        current = current / name
        if not _real_directory(current):
            raise BaseError(_STAGING_FAILED)
    return current


def _install_tree(parent: Path, name: str, staged: Path) -> None:
    if not _real_directory(parent):
        raise BaseError(_STAGING_FAILED)
    incoming: Path | None = parent / f".{name}.{os.getpid()}.incoming"
    if incoming.exists() or incoming.is_symlink():
        raise BaseError(_STAGING_FAILED)
    try:
        _copy_regular_tree(staged, incoming)
        destination = parent / name
        if destination.is_symlink():
            destination.unlink()
        elif _real_directory(destination):
            displaced = parent / f".{name}.{os.getpid()}.previous"
            os.rename(destination, displaced)
            try:
                os.rename(incoming, destination)
            except Exception:
                os.rename(displaced, destination)
                raise
            incoming = None
            _remove_owned_tree(displaced)
            return
        elif destination.exists():
            raise BaseError(_STAGING_FAILED)
        os.rename(incoming, destination)
        incoming = None
    finally:
        if incoming is not None and (incoming.exists() or incoming.is_symlink()):
            _remove_owned_tree(incoming)


def _install_file(parent: Path, name: str, staged: Path) -> None:
    if not _real_directory(parent):
        raise BaseError(_STAGING_FAILED)
    incoming = parent / f".{name}.{os.getpid()}.incoming"
    try:
        _copy_regular_file(staged, incoming)
        destination = parent / name
        if destination.is_symlink():
            destination.unlink()
        elif _real_directory(destination):
            raise BaseError(_STAGING_FAILED)
        os.replace(incoming, destination)
        incoming = None
    finally:
        if incoming is not None and (incoming.exists() or incoming.is_symlink()):
            if incoming.is_symlink() or _regular_file(incoming):
                incoming.unlink()


def _copy_regular_tree(source: Path, destination: Path) -> None:
    if not _real_directory(source) or destination.exists() or destination.is_symlink():
        raise BaseError(_STAGING_FAILED)
    destination.mkdir()
    for directory, dirnames, filenames in os.walk(source, followlinks=False):
        current = Path(directory)
        kept: list[str] = []
        for name in dirnames:
            child = current / name
            if not _real_directory(child):
                raise BaseError(_STAGING_FAILED)
            kept.append(name)
            (destination / child.relative_to(source)).mkdir()
        dirnames[:] = kept
        for name in filenames:
            child = current / name
            _copy_regular_file(child, destination / child.relative_to(source))


def _copy_regular_file(source: Path, destination: Path) -> None:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(source, flags)
    try:
        info = os.fstat(descriptor)
        if not stat.S_ISREG(info.st_mode):
            raise BaseError(_STAGING_FAILED)
        with destination.open("xb") as handle:
            while True:
                chunk = os.read(descriptor, 1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
    finally:
        os.close(descriptor)
    os.chmod(destination, 0o644)


def _remove_owned_tree(path: Path) -> None:
    try:
        mode = path.lstat().st_mode
    except OSError:
        return
    if stat.S_ISLNK(mode) or not stat.S_ISDIR(mode):
        path.unlink()
        return
    for directory, dirnames, filenames in os.walk(
        path, topdown=False, followlinks=False
    ):
        current = Path(directory)
        for name in filenames:
            child = current / name
            if _real_directory(child):
                continue
            child.unlink()
        for name in dirnames:
            child = current / name
            if child.is_symlink():
                child.unlink()
            elif _real_directory(child):
                child.rmdir()
    path.rmdir()


def _regular_files(root: Path) -> set[str]:
    files: set[str] = set()
    for directory, dirnames, filenames in os.walk(root, followlinks=False):
        current = Path(directory)
        kept: list[str] = []
        for name in dirnames:
            child = current / name
            if _real_directory(child):
                kept.append(name)
        dirnames[:] = kept
        for name in filenames:
            child = current / name
            if _regular_file(child):
                files.add(child.relative_to(root).as_posix())
    return files


def _regular_file(path: Path) -> bool:
    try:
        mode = path.lstat().st_mode
    except OSError:
        return False
    return stat.S_ISREG(mode) and not stat.S_ISLNK(mode)


def _real_directory(path: Path) -> bool:
    try:
        mode = path.lstat().st_mode
    except OSError:
        return False
    return stat.S_ISDIR(mode) and not stat.S_ISLNK(mode)


def _trees_overlap(source: Path, build: Path) -> bool:
    try:
        source_real = source.resolve()
        build_real = build.resolve()
    except OSError:
        return True
    return (
        source_real == build_real
        or _is_inside(source_real, build_real)
        or _is_inside(build_real, source_real)
    )


def _is_inside(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True
