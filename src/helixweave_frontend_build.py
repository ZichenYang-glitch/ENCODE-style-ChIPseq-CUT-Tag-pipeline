"""Keep one wheel's frontend static files aligned with the current source tree.

setuptools copies package data into ``build/lib`` and then packs that whole
tree. A later frontend build removes content-hashed chunks from source, but
the previous copies stay in ``build/lib`` and are shipped again. The install
check then reports ``FRONTEND_ASSET_INTEGRITY_FAILED``.

The command below removes only those extra files from the frontend static
directory of the build tree this command is already using. It does not remove
build directories, other packages, or a tree shared with the source checkout.
"""

from __future__ import annotations

import os
from pathlib import Path
import stat

from setuptools.command.build_py import build_py as _BuildPy


_PACKAGE = "encode_pipeline.frontend_assets"
_STATIC = "static"


class build_py(_BuildPy):
    """Package-data build that drops removed frontend static files."""

    def run(self) -> None:
        super().run()
        packages = set(self.packages or ())
        if (
            getattr(self, "editable_mode", False)
            or not self.build_lib
            or _PACKAGE not in packages
        ):
            return
        remove_stale_frontend_static_files(
            source_root=Path(self.get_package_dir(_PACKAGE)) / _STATIC,
            build_root=Path(self.build_lib, *_PACKAGE.split("."), _STATIC),
        )


def remove_stale_frontend_static_files(*, source_root: Path, build_root: Path) -> None:
    """Remove static files under ``build_root`` that ``source_root`` no longer has."""

    if not _real_directory(build_root) or not _real_directory(source_root):
        return
    if _trees_overlap(source_root, build_root):
        return
    source_files, source_directories = _snapshot(source_root)
    build_files, build_directories = _snapshot(build_root)
    for relative in build_files - source_files:
        _unlink_contained(build_root, relative)
    extras = sorted(
        build_directories - source_directories,
        key=lambda item: item.count("/"),
        reverse=True,
    )
    for relative in extras:
        directory = _contained(build_root, relative)
        if directory is None or not _real_directory(directory):
            continue
        if not any(directory.iterdir()):
            directory.rmdir()


def _real_directory(path: Path) -> bool:
    try:
        mode = path.lstat().st_mode
    except OSError:
        return False
    return stat.S_ISDIR(mode) and not stat.S_ISLNK(mode)


def _trees_overlap(source_root: Path, build_root: Path) -> bool:
    try:
        source = source_root.resolve()
        build = build_root.resolve()
    except OSError:
        return True
    return source == build or _is_inside(source, build) or _is_inside(build, source)


def _is_inside(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def _snapshot(root: Path) -> tuple[set[str], set[str]]:
    files: set[str] = set()
    directories: set[str] = set()
    for directory, dirnames, filenames in os.walk(root, followlinks=False):
        current = Path(directory)
        kept: list[str] = []
        for name in dirnames:
            child = current / name
            relative = child.relative_to(root).as_posix()
            mode = child.lstat().st_mode
            if stat.S_ISDIR(mode) and not stat.S_ISLNK(mode):
                directories.add(relative)
                kept.append(name)
            else:
                files.add(relative)
        dirnames[:] = kept
        for name in filenames:
            files.add((current / name).relative_to(root).as_posix())
    return files, directories


def _contained(root: Path, relative: str) -> Path | None:
    parts = tuple(relative.split("/"))
    if not parts or any(part in {"", ".", ".."} for part in parts):
        return None
    current = root
    try:
        for part in parts[:-1]:
            current = current / part
            mode = current.lstat().st_mode
            if not stat.S_ISDIR(mode) or stat.S_ISLNK(mode):
                return None
        return current / parts[-1]
    except OSError:
        return None


def _unlink_contained(root: Path, relative: str) -> None:
    path = _contained(root, relative)
    if path is None:
        return
    try:
        mode = path.lstat().st_mode
    except OSError:
        return
    if stat.S_ISDIR(mode) and not stat.S_ISLNK(mode):
        return
    path.unlink()
