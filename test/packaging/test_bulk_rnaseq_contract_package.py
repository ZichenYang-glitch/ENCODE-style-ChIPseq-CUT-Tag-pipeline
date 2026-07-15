"""Distribution evidence for pinned nf-core/rnaseq contract resources."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import zipfile


REPO_ROOT = Path(__file__).resolve().parents[2]
RESOURCE_ROOT = "encode_pipeline/contracts/nfcore_rnaseq"
EXPECTED = {
    "nextflow-schema-3.26.0.json": (
        60_620,
        "8f2f84a25c0aec65a18234cf01acdd74f2385e8dfac8417e4bad23a70bfb4388",
    ),
    "samplesheet-schema-3.26.0.json": (
        3_218,
        "013b669b1a3d38709f548a2548350ecdabb3f2ba578f2b7de843105e2ce87a7d",
    ),
    "LICENSE-3.26.0.txt": (
        1_075,
        "49147db9f580fc819494d436920fd8db97589e67989d8a82d5d18a4deb8d49f9",
    ),
}


def test_wheel_ships_exact_pinned_nfcore_contracts(tmp_path: Path) -> None:
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
        for filename, (size, digest) in EXPECTED.items():
            content = archive.read(f"{RESOURCE_ROOT}/{filename}")
            assert len(content) == size
            assert hashlib.sha256(content).hexdigest() == digest
        provenance = json.loads(archive.read(f"{RESOURCE_ROOT}/provenance.json"))

    assert provenance["release"] == "3.26.0"
    assert provenance["commit"] == "e7ca46272c8f9d5ceee3f71759f4ba551d3217a4"
    assert {entry["sha256"] for entry in provenance["files"]} == {
        digest for _, digest in EXPECTED.values()
    }
