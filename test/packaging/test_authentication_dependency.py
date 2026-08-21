"""Dependency and license contracts for the password-hashing adapter."""

from __future__ import annotations

from email.parser import BytesParser
from importlib import metadata as package_metadata
from pathlib import Path
import shutil
import subprocess
import sys
import tomllib
import zipfile

from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet
from packaging.version import Version


REPO_ROOT = Path(__file__).resolve().parents[2]
ARGON2_RANGE = SpecifierSet(">=25.1,<26")
CI_LOCK = REPO_ROOT / "workflow" / "envs" / "ci-fast.lock"
EXPECTED_RUNTIME_LICENSES = {
    "argon2-cffi": "MIT",
    "argon2-cffi-bindings": "MIT",
    "cffi": "MIT-0",
    "pycparser": "BSD-3-Clause",
}


def test_argon2_cffi_is_an_explicit_bounded_api_dependency() -> None:
    project = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))[
        "project"
    ]

    assert "argon2-cffi>=25.1,<26" in project["optional-dependencies"]["api"]
    assert "argon2-cffi>=25.1,<26" not in project["dependencies"]


def test_installed_argon2_chain_is_supported_and_permissively_licensed() -> None:
    distribution = package_metadata.distribution("argon2-cffi")
    metadata = distribution.metadata

    assert Version(distribution.version) in ARGON2_RANGE
    assert metadata["License-Expression"] == "MIT"
    python_range = SpecifierSet(metadata["Requires-Python"])
    for supported in ("3.10", "3.11", "3.12", "3.13"):
        assert Version(supported) in python_range
    for name, expected_license in EXPECTED_RUNTIME_LICENSES.items():
        assert package_metadata.metadata(name)["License-Expression"] == expected_license


def test_ci_lock_contains_the_complete_argon2_runtime_chain() -> None:
    lock = CI_LOCK.read_text(encoding="utf-8")

    assert "/argon2-cffi-25.1.0-" in lock
    assert "/argon2-cffi-bindings-25.1.0-" in lock
    assert "/cffi-2.1.1-" in lock
    assert "/pycparser-3.0-" in lock


def test_wheel_metadata_places_argon2_in_the_api_extra(tmp_path: Path) -> None:
    build_root = tmp_path / "source"
    build_root.mkdir()
    for filename in ("pyproject.toml", "README.md", "LICENSE", "MANIFEST.in"):
        shutil.copy2(REPO_ROOT / filename, build_root / filename)
    shutil.copytree(
        REPO_ROOT / "src",
        build_root / "src",
        ignore=shutil.ignore_patterns("*.egg-info", "__pycache__", "*.pyc"),
    )
    wheel_root = tmp_path / "wheel"
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "pip",
            "wheel",
            "--no-index",
            "--no-deps",
            "--no-build-isolation",
            "--wheel-dir",
            str(wheel_root),
            str(build_root),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    wheel = next(wheel_root.glob("helixweave-*.whl"))
    with zipfile.ZipFile(wheel) as archive:
        metadata_name = next(
            name for name in archive.namelist() if name.endswith(".dist-info/METADATA")
        )
        metadata = BytesParser().parsebytes(archive.read(metadata_name))

    requirements = [
        Requirement(value) for value in metadata.get_all("Requires-Dist") or ()
    ]
    argon2 = next(
        requirement for requirement in requirements if requirement.name == "argon2-cffi"
    )
    assert Version("25.1") in argon2.specifier
    assert Version("26") not in argon2.specifier
    assert argon2.marker is not None
    assert argon2.marker.evaluate({"extra": "api"}) is True
    assert argon2.marker.evaluate({"extra": ""}) is False
