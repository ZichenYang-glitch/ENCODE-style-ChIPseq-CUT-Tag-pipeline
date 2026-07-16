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
UPSTREAM_EXPECTED = {
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
RUNTIME_EXPECTED = {
    "source-manifest-3.26.0.json": (
        177_796,
        "dc75d105ad26b381197268ef67c44da3107b694e788a3c693237924a86ead774",
    ),
    "container-inventory-3.26.0.json": (
        17_516,
        "3cae9f36c2b872958dd82b2233cc178aa53fe7ab6111f7bdb0346f8e2bbbe9cf",
    ),
    "container-process-audit-3.26.0.json": (
        30_536,
        "45a2cabd6ec74de359f8d53faa2ca9423ac97cdd89cfa04e2c79c005418fd1c2",
    ),
    "container-availability-lock-1.0.0.schema.json": (
        3_931,
        "1774c6dfbdefbb8d03ca99a76ad88e893a07bc913d2534855911353d210fa636",
    ),
    "runtime-identity-3.26.0.json": (
        5_340,
        "5196c09a97b15ca8b51497234c13a0a13521f3b1819db69fcf4891aafc3a59c1",
    ),
    "results-contract-3.26.0.json": (
        10_586,
        "2e2a64389b66a6fcfb59d83281562b8973bbf3846542a4f20ddaa3bda3e0fe9f",
    ),
}
EXPECTED = {**UPSTREAM_EXPECTED, **RUNTIME_EXPECTED}


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
        execution_manifest = json.loads(
            archive.read(
                f"{RESOURCE_ROOT}/execution-implementation-manifest-1.0.0.json"
            )
        )

    assert provenance["release"] == "3.26.0"
    assert provenance["commit"] == "e7ca46272c8f9d5ceee3f71759f4ba551d3217a4"
    assert {entry["sha256"] for entry in provenance["files"]} == {
        digest for _, digest in UPSTREAM_EXPECTED.values()
    }
    assert execution_manifest["schema_version"] == "1.0.0"
    assert execution_manifest["file_count"] == len(execution_manifest["files"])
    assert execution_manifest["files"]
