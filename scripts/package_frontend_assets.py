#!/usr/bin/env python3
"""Package one completed Vite build without network or frontend tooling."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import stat
import sys
import tempfile
from typing import Any


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
sys.path.insert(0, str(SOURCE_ROOT))

from encode_pipeline.deployment.errors import DeploymentError  # noqa: E402
from encode_pipeline.frontend_assets import (  # noqa: E402
    build_frontend_assets,
    parse_manifest_bytes,
    verify_asset_directory,
)


DEFAULT_DIST_ROOT = REPOSITORY_ROOT / "frontend" / "dist"
DEFAULT_PACKAGE_JSON = REPOSITORY_ROOT / "frontend" / "package.json"
DEFAULT_PACKAGE_LOCK = REPOSITORY_ROOT / "frontend" / "package-lock.json"
DEFAULT_OPENAPI = REPOSITORY_ROOT / "frontend" / "openapi.json"
DEFAULT_OUTPUT_ROOT = SOURCE_ROOT / "encode_pipeline" / "frontend_assets"
_MAX_INPUT_BYTES = 16 * 1024 * 1024


def _read_input(path: Path) -> bytes:
    descriptor: int | None = None
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        before = os.lstat(path)
        if (
            not stat.S_ISREG(before.st_mode)
            or not 0 < before.st_size <= _MAX_INPUT_BYTES
        ):
            raise ValueError
        descriptor = os.open(path, flags)
        opened = os.fstat(descriptor)
        if not stat.S_ISREG(opened.st_mode) or (opened.st_dev, opened.st_ino) != (
            before.st_dev,
            before.st_ino,
        ):
            raise ValueError
        content = bytearray()
        while len(content) <= _MAX_INPUT_BYTES:
            chunk = os.read(
                descriptor, min(1024 * 1024, _MAX_INPUT_BYTES + 1 - len(content))
            )
            if not chunk:
                break
            content.extend(chunk)
        after = os.lstat(path)
        if (
            len(content) != opened.st_size
            or len(content) > _MAX_INPUT_BYTES
            or (after.st_dev, after.st_ino) != (opened.st_dev, opened.st_ino)
        ):
            raise ValueError
        return bytes(content)
    finally:
        if descriptor is not None:
            os.close(descriptor)


def _json_pairs(pairs: list[tuple[str, object]]) -> dict[str, object]:
    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError
        value[key] = item
    return value


def _frontend_version(package_json: bytes) -> str:
    value: Any = json.loads(
        package_json.decode("utf-8"),
        object_pairs_hook=_json_pairs,
        parse_constant=lambda _value: (_ for _ in ()).throw(ValueError()),
    )
    if (
        not isinstance(value, dict)
        or value.get("name") != "helixweave-frontend"
        or not isinstance(value.get("version"), str)
    ):
        raise ValueError
    return value["version"]


def _write_payload(root: Path, content: dict[str, bytes]) -> None:
    root.mkdir(parents=True)
    for relative, value in sorted(content.items()):
        destination = root.joinpath(*relative.split("/"))
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("xb") as handle:
            handle.write(value)
        destination.chmod(0o644)


def package_frontend_assets(
    *,
    dist_root: Path,
    package_json: Path,
    package_lock: Path,
    openapi: Path,
    output_root: Path,
) -> str:
    """Atomically replace only the package-owned generated frontend payload."""
    package_json_bytes = _read_input(package_json)
    package_lock_bytes = _read_input(package_lock)
    openapi_bytes = _read_input(openapi)
    verified = build_frontend_assets(
        dist_root,
        frontend_version=_frontend_version(package_json_bytes),
        package_lock_bytes=package_lock_bytes,
        openapi_bytes=openapi_bytes,
    )

    output_root.parent.mkdir(parents=True, exist_ok=True)
    if output_root.exists():
        output_mode = os.lstat(output_root).st_mode
        if not stat.S_ISDIR(output_mode):
            raise ValueError
    else:
        output_root.mkdir(mode=0o755)

    temporary_root = Path(
        tempfile.mkdtemp(prefix=".frontend-assets-", dir=output_root.parent)
    )
    staged_static = temporary_root / "static"
    previous_static = temporary_root / "previous-static"
    target_static = output_root / "static"
    manifest_path = output_root / "asset-manifest.json"
    temporary_manifest = output_root / f".asset-manifest.{os.getpid()}.tmp"
    installed_new_static = False
    moved_previous_static = False
    try:
        _write_payload(staged_static, dict(verified.content))
        manifest_bytes = verified.manifest.to_bytes()
        parsed = parse_manifest_bytes(manifest_bytes)
        verify_asset_directory(parsed, staged_static)

        with temporary_manifest.open("xb") as handle:
            handle.write(manifest_bytes)
        temporary_manifest.chmod(0o644)

        if os.path.lexists(target_static):
            os.replace(target_static, previous_static)
            moved_previous_static = True
        os.replace(staged_static, target_static)
        installed_new_static = True
        os.replace(temporary_manifest, manifest_path)
        if moved_previous_static:
            shutil.rmtree(previous_static, ignore_errors=True)
        return verified.manifest.identity
    except Exception:
        temporary_manifest.unlink(missing_ok=True)
        if installed_new_static and os.path.lexists(target_static):
            if target_static.is_dir() and not target_static.is_symlink():
                shutil.rmtree(target_static)
            else:
                target_static.unlink()
        if moved_previous_static and os.path.lexists(previous_static):
            os.replace(previous_static, target_static)
        raise
    finally:
        shutil.rmtree(temporary_root, ignore_errors=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Package a completed HelixWeave frontend build.",
    )
    parser.add_argument("--dist-root", type=Path, default=DEFAULT_DIST_ROOT)
    parser.add_argument("--package-json", type=Path, default=DEFAULT_PACKAGE_JSON)
    parser.add_argument("--package-lock", type=Path, default=DEFAULT_PACKAGE_LOCK)
    parser.add_argument("--openapi", type=Path, default=DEFAULT_OPENAPI)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    return parser


def main(argv: list[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    try:
        identity = package_frontend_assets(
            dist_root=arguments.dist_root,
            package_json=arguments.package_json,
            package_lock=arguments.package_lock,
            openapi=arguments.openapi,
            output_root=arguments.output_root,
        )
    except DeploymentError as error:
        print(
            f"frontend asset packaging failed [{error.issue.code}]: "
            f"{error.issue.message}",
            file=sys.stderr,
        )
        return 2
    except (OSError, UnicodeDecodeError, ValueError, json.JSONDecodeError):
        print(
            "frontend asset packaging failed [FRONTEND_ASSET_BUILD_INVALID]: "
            "Frontend build inputs are invalid.",
            file=sys.stderr,
        )
        return 2
    print(
        json.dumps(
            {
                "ok": True,
                "schema_version": "helixweave-frontend-package-receipt-v1",
                "frontend_identity": identity,
            },
            separators=(",", ":"),
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
