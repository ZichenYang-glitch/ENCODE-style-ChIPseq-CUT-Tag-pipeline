"""Repeated wheel builds ship the current frontend closure only."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

import pytest
from setuptools.errors import BaseError

from encode_pipeline.frontend_assets import build_frontend_assets
from helixweave_frontend_build import stage_frontend_assets


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
_PACKAGE_ROOT = "encode_pipeline/frontend_assets"
_LOCK = b'{"lockfileVersion":3}\n'
_OPENAPI = b"{}\n"
_GENERATION_A = {
    "index.html": b"<html>A</html>\n",
    "assets/chunk-A.js": b"console.log('A')\n",
    "assets/shared.js": b"console.log('shared')\n",
}
_GENERATION_B = {
    "index.html": b"<html>B</html>\n",
    "assets/chunk-B.js": b"console.log('B')\n",
    "assets/shared.js": b"console.log('shared')\n",
}
_FROZEN_MTIME = 1700000000
_PACKAGE_JSON = b'{"name":"helixweave-frontend","version":"0.4.1"}\n'


def _checkout(tmp_path: Path) -> Path:
    checkout = tmp_path / "checkout"
    checkout.mkdir()
    for name in ("pyproject.toml", "README.md", "LICENSE", "MANIFEST.in"):
        shutil.copy2(REPOSITORY_ROOT / name, checkout / name)
    shutil.copytree(
        REPOSITORY_ROOT / "src",
        checkout / "src",
        ignore=shutil.ignore_patterns("*.egg-info", "__pycache__", "*.pyc"),
    )
    return checkout


def _write_frontend(checkout: Path, files: dict[str, bytes]) -> str:
    dist = checkout / ".synthetic-dist"
    if dist.exists():
        shutil.rmtree(dist)
    for relative, content in files.items():
        destination = dist.joinpath(*relative.split("/"))
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
    verified = build_frontend_assets(
        dist,
        frontend_version="0.4.1",
        package_lock_bytes=_LOCK,
        openapi_bytes=_OPENAPI,
    )
    package = checkout / "src" / "encode_pipeline" / "frontend_assets"
    static = package / "static"
    if static.exists():
        shutil.rmtree(static)
    for relative, content in verified.content.items():
        destination = static.joinpath(*relative.split("/"))
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
    manifest_path = package / "asset-manifest.json"
    manifest_path.write_bytes(verified.manifest.to_bytes())
    shutil.rmtree(dist)
    return verified.manifest.identity


def _build_wheel(checkout: Path, outdir: Path, *, sdist: bool = False) -> Path:
    outdir.mkdir()
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "build",
            "--no-isolation",
            *(("--sdist", "--wheel") if sdist else ("--wheel",)),
            "--outdir",
            str(outdir),
            str(checkout),
        ],
        cwd=checkout,
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
    )
    assert completed.returncode == 0, completed.stderr
    wheels = tuple(outdir.glob("helixweave-*.whl"))
    assert len(wheels) == 1
    return wheels[0]


def _static_members(wheel: Path) -> dict[str, bytes]:
    prefix = f"{_PACKAGE_ROOT}/static/"
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        assert "helixweave_frontend_build.py" not in names
        assert not any(name.endswith("/helixweave_frontend_build.py") for name in names)
        return {
            name.removeprefix(prefix): archive.read(name)
            for name in names
            if name.startswith(prefix) and not name.endswith("/")
        }


def _manifest_bytes(wheel: Path) -> bytes:
    with zipfile.ZipFile(wheel) as archive:
        return archive.read(f"{_PACKAGE_ROOT}/asset-manifest.json")


def _load_packaged(tmp_path: Path, wheel: Path) -> dict[str, object]:
    installed = tmp_path / "installed"
    outside = tmp_path / "outside"
    outside.mkdir()
    with zipfile.ZipFile(wheel) as archive:
        archive.extractall(installed)
    code = """
import json
import sys

sys.path.insert(0, sys.argv[1])
from encode_pipeline.frontend_assets import load_packaged_frontend_assets

assets = load_packaged_frontend_assets()
print(json.dumps({
    "identity": assets.manifest.identity,
    "files": sorted(assets.content),
}))
"""
    completed = subprocess.run(
        [sys.executable, "-I", "-S", "-c", code, str(installed)],
        cwd=outside,
        env={
            "PATH": os.environ.get("PATH", ""),
            "PYTHONNOUSERSITE": "1",
            "PYTHONDONTWRITEBYTECODE": "1",
        },
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert completed.returncode == 0, completed.stderr
    receipt = json.loads(completed.stdout)
    assert isinstance(receipt, dict)
    return receipt


def test_repeated_wheel_build_replaces_previous_frontend_chunks(
    tmp_path: Path,
) -> None:
    checkout = _checkout(tmp_path)
    _write_frontend(checkout, _GENERATION_A)
    _build_wheel(checkout, tmp_path / "dist-a")

    identity_b = _write_frontend(checkout, _GENERATION_B)
    parallel_build = checkout / "parallel-build" / "lib" / "static" / "chunk-A.js"
    parallel_build.parent.mkdir(parents=True)
    parallel_build.write_bytes(b"parallel-build")
    user_directory = tmp_path / "user-build" / "static"
    user_directory.mkdir(parents=True)
    user_file = user_directory / "chunk-A.js"
    user_file.write_bytes(b"user-build")
    marker = checkout / "build" / "lib" / "encode_pipeline" / "services" / "marker.txt"
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_bytes(b"keep-other-package")

    wheel = _build_wheel(checkout, tmp_path / "dist-b")

    assert _static_members(wheel) == _GENERATION_B
    manifest = json.loads(_manifest_bytes(wheel).decode("utf-8"))
    assert manifest["identity"] == identity_b
    assert [item["path"] for item in manifest["files"]] == sorted(_GENERATION_B)
    receipt = _load_packaged(tmp_path, wheel)
    assert receipt == {"identity": identity_b, "files": sorted(_GENERATION_B)}
    stale = (
        checkout
        / "build"
        / "lib"
        / "encode_pipeline"
        / "frontend_assets"
        / "static"
        / "assets"
        / "chunk-A.js"
    )
    assert not stale.exists()
    assert parallel_build.read_bytes() == b"parallel-build"
    assert user_file.read_bytes() == b"user-build"
    assert marker.read_bytes() == b"keep-other-package"
    assert (checkout / "build").is_dir()


def test_clean_wheel_build_matches_frontend_manifest(tmp_path: Path) -> None:
    checkout = _checkout(tmp_path)
    identity = _write_frontend(checkout, _GENERATION_B)
    assert not (checkout / "build").exists()

    wheel = _build_wheel(checkout, tmp_path / "dist")

    assert _static_members(wheel) == _GENERATION_B
    receipt = _load_packaged(tmp_path, wheel)
    assert receipt == {"identity": identity, "files": sorted(_GENERATION_B)}
    source_static = checkout / "src" / "encode_pipeline" / "frontend_assets" / "static"
    built_static = (
        checkout / "build" / "lib" / "encode_pipeline" / "frontend_assets" / "static"
    )
    source_files = {
        path.relative_to(source_static).as_posix()
        for path in source_static.rglob("*")
        if path.is_file()
    }
    built_files = {
        path.relative_to(built_static).as_posix()
        for path in built_static.rglob("*")
        if path.is_file()
    }
    assert built_files == source_files


def _snapshot(root: Path) -> dict[str, bytes | None]:
    recorded: dict[str, bytes | None] = {}
    for directory, dirnames, filenames in os.walk(root, followlinks=False):
        current = Path(directory)
        relative_dir = current.relative_to(root).as_posix()
        if relative_dir != ".":
            recorded[relative_dir + "/"] = None
        for name in sorted(dirnames):
            child = current / name
            if child.is_symlink():
                recorded[child.relative_to(root).as_posix() + "/"] = b"symlink"
        for name in filenames:
            child = current / name
            relative = child.relative_to(root).as_posix()
            if child.is_symlink():
                recorded[relative] = b"symlink:" + os.readlink(child).encode()
            else:
                recorded[relative] = child.read_bytes()
    return recorded


def _external_tree(root: Path) -> Path:
    files = {
        "encode_pipeline/frontend_assets/asset-manifest.json": b'{"generation":"external"}\n',
        "encode_pipeline/frontend_assets/static/index.html": b"<html>external</html>\n",
        "encode_pipeline/frontend_assets/static/assets/chunk-A.js": b"console.log('external')\n",
        "notes/readme.txt": b"keep this tree\n",
        "notes/nested/data.bin": b"\x00\x01external",
    }
    for relative, content in files.items():
        destination = root.joinpath(*relative.split("/"))
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
    (root / "notes" / "empty").mkdir()
    return root


def _freeze_mtime(package: Path, timestamp: int) -> None:
    targets = [package / "asset-manifest.json", package / "static"]
    targets.extend((package / "static").rglob("*"))
    for path in targets:
        if path.is_symlink():
            continue
        os.utime(path, (timestamp, timestamp))


def _package_with_script(checkout: Path, files: dict[str, bytes]) -> str:
    work = checkout / ".synthetic-inputs"
    if work.exists():
        shutil.rmtree(work)
    dist = work / "dist"
    for relative, content in files.items():
        destination = dist.joinpath(*relative.split("/"))
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
    (work / "package.json").write_bytes(_PACKAGE_JSON)
    (work / "package-lock.json").write_bytes(_LOCK)
    (work / "openapi.json").write_bytes(_OPENAPI)
    script_dir = checkout / "scripts"
    script_dir.mkdir(exist_ok=True)
    script = script_dir / "package_frontend_assets.py"
    shutil.copy2(
        REPOSITORY_ROOT / "scripts" / "package_frontend_assets.py",
        script,
    )
    completed = subprocess.run(
        [
            sys.executable,
            str(script),
            "--dist-root",
            str(dist),
            "--package-json",
            str(work / "package.json"),
            "--package-lock",
            str(work / "package-lock.json"),
            "--openapi",
            str(work / "openapi.json"),
            "--output-root",
            str(checkout / "src" / "encode_pipeline" / "frontend_assets"),
        ],
        cwd=checkout,
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert completed.returncode == 0, completed.stderr
    manifest = json.loads(
        (
            checkout
            / "src"
            / "encode_pipeline"
            / "frontend_assets"
            / "asset-manifest.json"
        ).read_text(encoding="utf-8")
    )
    return str(manifest["identity"])


def _write_package(package: Path, files: dict[str, bytes], manifest: bytes) -> None:
    static = package / "static"
    static.mkdir(parents=True, exist_ok=True)
    for relative, content in files.items():
        destination = static.joinpath(*relative.split("/"))
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
    (package / "asset-manifest.json").write_bytes(manifest)


def test_same_mtime_rebuild_loads_the_replacement_frontend(tmp_path: Path) -> None:
    checkout = _checkout(tmp_path)
    identity_a = _package_with_script(checkout, _GENERATION_A)
    _freeze_mtime(
        checkout / "src" / "encode_pipeline" / "frontend_assets",
        _FROZEN_MTIME,
    )
    _build_wheel(checkout, tmp_path / "dist-a", sdist=True)

    identity_b = _package_with_script(checkout, _GENERATION_B)
    package = checkout / "src" / "encode_pipeline" / "frontend_assets"
    cached = checkout / "build" / "lib" / "encode_pipeline" / "frontend_assets"
    _freeze_mtime(package, _FROZEN_MTIME)
    _freeze_mtime(cached, _FROZEN_MTIME)
    assert (package / "static" / "index.html").read_bytes() == _GENERATION_B[
        "index.html"
    ]
    assert (cached / "static" / "index.html").read_bytes() == _GENERATION_A[
        "index.html"
    ]
    assert (package / "static" / "index.html").stat().st_mtime == _FROZEN_MTIME
    assert (cached / "static" / "index.html").stat().st_mtime == _FROZEN_MTIME
    assert (cached / "asset-manifest.json").read_bytes() != (
        package / "asset-manifest.json"
    ).read_bytes()

    external = _external_tree(tmp_path / "external")
    before = _snapshot(external)
    wheel = _build_wheel(checkout, tmp_path / "dist-b", sdist=True)

    assert _static_members(wheel) == _GENERATION_B
    manifest = json.loads(_manifest_bytes(wheel).decode("utf-8"))
    assert manifest["identity"] == identity_b
    assert identity_b != identity_a
    assert [item["path"] for item in manifest["files"]] == sorted(_GENERATION_B)
    receipt = _load_packaged(tmp_path, wheel)
    assert receipt == {"identity": identity_b, "files": sorted(_GENERATION_B)}
    assert _snapshot(external) == before


def test_frontend_staging_refreshes_files_with_equal_mtime(tmp_path: Path) -> None:
    source = tmp_path / "source"
    _write_package(source, _GENERATION_B, b"manifest-b")
    _freeze_mtime(source, _FROZEN_MTIME)
    build_lib = tmp_path / "build" / "lib"
    package = build_lib / "encode_pipeline" / "frontend_assets"
    _write_package(package, _GENERATION_A, b"manifest-a")
    _freeze_mtime(package, _FROZEN_MTIME)
    marker = build_lib / "encode_pipeline" / "services" / "marker.txt"
    marker.parent.mkdir(parents=True)
    marker.write_bytes(b"keep-other-package")
    external = _external_tree(tmp_path / "external")
    before = _snapshot(external)

    stage_frontend_assets(source_package=source, build_lib=build_lib)

    assert (package / "static" / "index.html").read_bytes() == _GENERATION_B[
        "index.html"
    ]
    assert (package / "static" / "assets" / "chunk-B.js").read_bytes() == (
        _GENERATION_B["assets/chunk-B.js"]
    )
    assert not (package / "static" / "assets" / "chunk-A.js").exists()
    assert (package / "asset-manifest.json").read_bytes() == b"manifest-b"
    assert marker.read_bytes() == b"keep-other-package"
    assert _snapshot(external) == before


def test_frontend_staging_leaves_symlinked_trees_unchanged(tmp_path: Path) -> None:
    source = tmp_path / "source"
    _write_package(source, _GENERATION_B, b"manifest-b")
    external = _external_tree(tmp_path / "external")
    before = _snapshot(external)

    linked_lib = tmp_path / "linked-lib"
    linked_lib.symlink_to(external, target_is_directory=True)
    with pytest.raises(BaseError):
        stage_frontend_assets(source_package=source, build_lib=linked_lib)
    assert _snapshot(external) == before

    build_lib = tmp_path / "build" / "lib"
    build_lib.mkdir(parents=True)
    (build_lib / "encode_pipeline").symlink_to(
        external / "encode_pipeline",
        target_is_directory=True,
    )
    with pytest.raises(BaseError):
        stage_frontend_assets(source_package=source, build_lib=build_lib)
    assert (build_lib / "encode_pipeline").is_symlink()
    assert _snapshot(external) == before

    build_lib = tmp_path / "real-build" / "lib"
    package = build_lib / "encode_pipeline" / "frontend_assets"
    package.mkdir(parents=True)
    (package / "asset-manifest.json").write_bytes(b"manifest-a")
    (package / "static").symlink_to(
        external / "encode_pipeline" / "frontend_assets" / "static",
        target_is_directory=True,
    )
    stage_frontend_assets(source_package=source, build_lib=build_lib)
    assert _snapshot(external) == before
    assert (package / "static").is_dir()
    assert not (package / "static").is_symlink()
    assert (package / "static" / "index.html").read_bytes() == _GENERATION_B[
        "index.html"
    ]
    assert (package / "asset-manifest.json").read_bytes() == b"manifest-b"

    overlapped = package
    (overlapped / "static" / "assets" / "chunk-A.js").write_bytes(b"do-not-delete")
    with pytest.raises(BaseError):
        stage_frontend_assets(source_package=overlapped, build_lib=build_lib)
    assert (
        overlapped / "static" / "assets" / "chunk-A.js"
    ).read_bytes() == b"do-not-delete"
