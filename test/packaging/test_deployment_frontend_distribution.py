"""Wheel/sdist delivery evidence for the precompiled production frontend."""

from __future__ import annotations

import base64
import csv
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import zipfile

import pytest

from encode_pipeline.frontend_assets import load_packaged_frontend_assets


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
PACKAGE_PREFIX = "encode_pipeline/frontend_assets"


def _build_distributions(tmp_path: Path) -> tuple[Path, Path]:
    source = tmp_path / "source"
    source.mkdir()
    for filename in ("pyproject.toml", "README.md", "LICENSE", "MANIFEST.in"):
        shutil.copy2(REPOSITORY_ROOT / filename, source / filename)
    shutil.copytree(
        REPOSITORY_ROOT / "src",
        source / "src",
        ignore=shutil.ignore_patterns("*.egg-info", "__pycache__", "*.pyc"),
    )
    destination = tmp_path / "dist"
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "build",
            "--no-isolation",
            "--sdist",
            "--wheel",
            "--outdir",
            str(destination),
            str(source),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    wheels = tuple(destination.glob("helixweave-*.whl"))
    sdists = tuple(destination.glob("helixweave-*.tar.gz"))
    assert len(wheels) == len(sdists) == 1
    return wheels[0], sdists[0]


def _record_rows(archive: zipfile.ZipFile) -> dict[str, tuple[str, str]]:
    record_name = next(
        name for name in archive.namelist() if name.endswith(".dist-info/RECORD")
    )
    return {
        row[0]: (row[1], row[2])
        for row in csv.reader(io.StringIO(archive.read(record_name).decode("utf-8")))
    }


def _record_digest(content: bytes) -> str:
    encoded = base64.urlsafe_b64encode(hashlib.sha256(content).digest()).rstrip(b"=")
    return f"sha256={encoded.decode('ascii')}"


@pytest.fixture(scope="module")
def frontend_distributions(
    tmp_path_factory: pytest.TempPathFactory,
) -> tuple[Path, Path]:
    return _build_distributions(tmp_path_factory.mktemp("frontend-distributions"))


def test_wheel_and_sdist_ship_the_exact_precompiled_frontend_closure(
    frontend_distributions: tuple[Path, Path],
) -> None:
    packaged = load_packaged_frontend_assets()
    wheel, sdist = frontend_distributions
    expected = {
        f"{PACKAGE_PREFIX}/asset-manifest.json": packaged.manifest.to_bytes(),
        **{
            f"{PACKAGE_PREFIX}/static/{path}": content
            for path, content in packaged.content.items()
        },
    }

    with zipfile.ZipFile(wheel) as archive:
        names = set(archive.namelist())
        rows = _record_rows(archive)
        for name, content in expected.items():
            assert name in names
            assert archive.read(name) == content
            assert rows[name] == (_record_digest(content), str(len(content)))
        frontend_members = {
            name for name in names if name.startswith(f"{PACKAGE_PREFIX}/")
        }
        expected_members = set(expected) | {
            f"{PACKAGE_PREFIX}/__init__.py",
            f"{PACKAGE_PREFIX}/manifest.py",
        }
        assert frontend_members == expected_members
        assert not any(name.endswith(".map") for name in frontend_members)
        assert not any("node_modules" in name for name in names)
        assert not any(name.startswith("frontend/") for name in names)

    with tarfile.open(sdist, "r:gz") as archive:
        members = {item.name: item for item in archive.getmembers() if item.isfile()}
        prefix = next(name.split("/", 1)[0] for name in members)
        for name, content in expected.items():
            member_name = f"{prefix}/src/{name}"
            assert member_name in members
            extracted = archive.extractfile(members[member_name])
            assert extracted is not None
            assert extracted.read() == content


def test_extracted_wheel_verifies_frontend_without_node_or_source_checkout(
    tmp_path: Path,
    frontend_distributions: tuple[Path, Path],
) -> None:
    wheel, _sdist = frontend_distributions
    installed = tmp_path / "installed"
    with zipfile.ZipFile(wheel) as archive:
        archive.extractall(installed)
    outside = tmp_path / "outside"
    outside.mkdir()
    code = """
import json
from encode_pipeline.frontend_assets import load_packaged_frontend_assets

assets = load_packaged_frontend_assets()
print(json.dumps({
    "identity": assets.manifest.identity,
    "version": assets.manifest.frontend_version,
    "files": sorted(assets.content),
}, sort_keys=True))
"""
    environment = {
        "PATH": "/node-and-npm-are-intentionally-unavailable",
        "PYTHONNOUSERSITE": "1",
        "PYTHONDONTWRITEBYTECODE": "1",
    }
    isolated_code = f"import sys; sys.path.insert(0, {str(installed)!r})\n" + code
    completed = subprocess.run(
        [sys.executable, "-I", "-S", "-c", isolated_code],
        cwd=outside,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    receipt = json.loads(completed.stdout)
    assert receipt["identity"].startswith("sha256-")
    assert receipt["version"] == "0.3.0"
    assert "index.html" in receipt["files"]
    assert all("node_modules" not in path for path in receipt["files"])
