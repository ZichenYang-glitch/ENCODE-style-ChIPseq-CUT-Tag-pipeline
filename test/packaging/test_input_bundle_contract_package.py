"""Distribution evidence for the pinned public Bundle contract."""

from __future__ import annotations

import hashlib
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile


REPO_ROOT = Path(__file__).resolve().parents[2]
SCHEMA_MEMBER = "encode_pipeline/contracts/omics_intake/intake-bundle-0.2.schema.json"
SCHEMA_SHA256 = "6dcda336c9f0ba763383ddd58bec280946f8970af8bd730eb758fab5e3a8dd71"


def test_wheel_ships_pinned_bundle_schema_and_runtime_validator(tmp_path: Path) -> None:
    build_root = tmp_path / "source"
    build_root.mkdir()
    shutil.copy2(REPO_ROOT / "pyproject.toml", build_root / "pyproject.toml")
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
    wheels = tuple(wheel_root.glob("encode_pipeline-*.whl"))
    assert len(wheels) == 1

    with zipfile.ZipFile(wheels[0]) as archive:
        schema = archive.read(SCHEMA_MEMBER)
        metadata_name = next(
            name for name in archive.namelist() if name.endswith(".dist-info/METADATA")
        )
        metadata = archive.read(metadata_name).decode("utf-8")

    assert hashlib.sha256(schema).hexdigest() == SCHEMA_SHA256
    assert "Requires-Dist: jsonschema" in metadata
