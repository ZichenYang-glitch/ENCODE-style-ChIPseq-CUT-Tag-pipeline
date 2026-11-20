"""Repeated wheel builds ship the current frontend closure only."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile

from encode_pipeline.frontend_assets import build_frontend_assets
from helixweave_frontend_build import remove_stale_frontend_static_files


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


def _build_wheel(checkout: Path, outdir: Path) -> Path:
    outdir.mkdir()
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "build",
            "--no-isolation",
            "--wheel",
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


def test_stale_frontend_pruner_leaves_other_directories_in_place(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source"
    build = tmp_path / "build-static"
    (source / "assets").mkdir(parents=True)
    (source / "index.html").write_bytes(b"B")
    (source / "assets" / "chunk-B.js").write_bytes(b"B")
    (build / "assets").mkdir(parents=True)
    (build / "index.html").write_bytes(b"B")
    (build / "assets" / "chunk-B.js").write_bytes(b"B")
    (build / "assets" / "chunk-A.js").write_bytes(b"A")
    (build / "vendor").mkdir()
    (build / "vendor" / "old.js").write_bytes(b"old")

    user = tmp_path / "user-files"
    user.mkdir()
    secret = user / "secret.txt"
    secret.write_bytes(b"keep")
    linked = build / "linked-user"
    linked.symlink_to(user, target_is_directory=True)
    parallel = tmp_path / "parallel-build" / "chunk-A.js"
    parallel.parent.mkdir()
    parallel.write_bytes(b"other-build")

    remove_stale_frontend_static_files(source_root=source, build_root=build)

    assert (build / "index.html").read_bytes() == b"B"
    assert (build / "assets" / "chunk-B.js").read_bytes() == b"B"
    assert not (build / "assets" / "chunk-A.js").exists()
    assert not (build / "vendor").exists()
    assert secret.read_bytes() == b"keep"
    assert user.is_dir()
    assert not linked.exists()
    assert parallel.read_bytes() == b"other-build"

    (build / "chunk-A.js").write_bytes(b"same-tree")
    remove_stale_frontend_static_files(source_root=build, build_root=build)
    assert (build / "chunk-A.js").read_bytes() == b"same-tree"

    outside = tmp_path / "outside-static"
    outside.mkdir()
    outside_file = outside / "chunk-A.js"
    outside_file.write_bytes(b"keep-outside")
    alias = tmp_path / "alias-static"
    alias.symlink_to(outside, target_is_directory=True)
    remove_stale_frontend_static_files(source_root=source, build_root=alias)
    assert outside_file.read_bytes() == b"keep-outside"
    assert alias.is_symlink()
