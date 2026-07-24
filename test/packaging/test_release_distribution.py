"""Release identity and installed-distribution contracts for HelixWeave."""

from __future__ import annotations

from email.parser import Parser
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tomllib
import zipfile

import yaml

from encode_pipeline import __version__
from encode_pipeline.api.main import create_app


REPO_ROOT = Path(__file__).resolve().parents[2]
RELEASE_VERSION = "0.3.0"
EXPECTED_CONSOLE_SCRIPTS = {
    "encode-dag": "encode_pipeline.cli.dag:main",
    "encode-manifest": "encode_pipeline.cli.manifest:main",
    "encode-validate": "encode_pipeline.cli.validate:main",
    "encode-worker": "encode_pipeline.workers.cli:main",
}


def _build_wheel(tmp_path: Path) -> Path:
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
    wheels = tuple(wheel_root.glob("encode_pipeline-*.whl"))
    assert len(wheels) == 1
    return wheels[0]


def test_release_identity_is_consistent_across_public_metadata(
    tmp_path: Path,
) -> None:
    pyproject = tomllib.loads(
        (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    )
    frontend = json.loads(
        (REPO_ROOT / "frontend/package.json").read_text(encoding="utf-8")
    )
    frontend_lock = json.loads(
        (REPO_ROOT / "frontend/package-lock.json").read_text(encoding="utf-8")
    )
    citation = yaml.safe_load((REPO_ROOT / "CITATION.cff").read_text(encoding="utf-8"))
    app = create_app(
        database_url=f"sqlite:///{tmp_path / 'release-identity.db'}",
        workspace_root=tmp_path / "workspaces",
        project_root=REPO_ROOT,
    )
    try:
        assert pyproject["project"]["name"] == "encode-pipeline"
        assert pyproject["project"]["version"] == RELEASE_VERSION
        assert pyproject["project"]["description"].startswith("HelixWeave")
        assert __version__ == RELEASE_VERSION
        assert app.title == "HelixWeave API"
        assert app.version == RELEASE_VERSION
        assert frontend["name"] == "helixweave-frontend"
        assert frontend["private"] is True
        assert frontend["version"] == RELEASE_VERSION
        assert frontend_lock["version"] == RELEASE_VERSION
        assert frontend_lock["packages"][""]["version"] == RELEASE_VERSION
        assert citation["title"] == "HelixWeave"
        assert citation["version"] == RELEASE_VERSION
    finally:
        app.state.run_queue.close()
        app.state.persistence.close()


def test_packaged_artifact_inventory_matches_documented_contract() -> None:
    documented = REPO_ROOT / "docs/architecture/artifact-inventory.yaml"
    packaged = REPO_ROOT / "src/encode_pipeline/artifacts/artifact-inventory.yaml"

    assert packaged.read_bytes() == documented.read_bytes()


def test_wheel_metadata_entrypoints_and_runtime_resources(tmp_path: Path) -> None:
    wheel = _build_wheel(tmp_path)

    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        metadata_name = next(
            name for name in names if name.endswith(".dist-info/METADATA")
        )
        entry_points_name = next(
            name for name in names if name.endswith(".dist-info/entry_points.txt")
        )
        metadata_text = archive.read(metadata_name).decode("utf-8")
        metadata = Parser().parsestr(metadata_text)
        entry_points = archive.read(entry_points_name).decode("utf-8")

        assert metadata["Name"] == "encode-pipeline"
        assert metadata["Version"] == RELEASE_VERSION
        assert metadata["Summary"].startswith("HelixWeave")
        assert metadata["Description-Content-Type"] == "text/markdown"
        assert "\n# HelixWeave\n" in metadata_text
        for name, target in EXPECTED_CONSOLE_SCRIPTS.items():
            assert f"{name} = {target}" in entry_points
        assert "encode_pipeline/artifacts/artifact-inventory.yaml" in names
        assert "encode_pipeline/persistence/alembic/script.py.mako" in names
        assert (
            sum(
                "/persistence/alembic/versions/" in name and name.endswith(".py")
                for name in names
            )
            == 9
        )
        assert not any(name.startswith("test/") for name in names)


def test_pep517_build_produces_bounded_wheel_and_sdist(tmp_path: Path) -> None:
    build_root = tmp_path / "source"
    build_root.mkdir()
    for filename in ("pyproject.toml", "README.md", "LICENSE", "MANIFEST.in"):
        shutil.copy2(REPO_ROOT / filename, build_root / filename)
    shutil.copytree(
        REPO_ROOT / "src",
        build_root / "src",
        ignore=shutil.ignore_patterns("*.egg-info", "__pycache__", "*.pyc"),
    )
    distribution_root = tmp_path / "dist"
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "build",
            "--no-isolation",
            "--sdist",
            "--wheel",
            "--outdir",
            str(distribution_root),
            str(build_root),
        ],
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert {path.name for path in distribution_root.iterdir()} == {
        "encode_pipeline-0.3.0-py3-none-any.whl",
        "encode_pipeline-0.3.0.tar.gz",
    }
    with tarfile.open(
        distribution_root / "encode_pipeline-0.3.0.tar.gz", "r:gz"
    ) as archive:
        names = archive.getnames()
    assert "encode_pipeline-0.3.0/README.md" in names
    assert (
        "encode_pipeline-0.3.0/src/encode_pipeline/artifacts/artifact-inventory.yaml"
    ) in names
    assert not any(name.startswith("encode_pipeline-0.3.0/test/") for name in names)


def test_extracted_wheel_supports_registry_and_openapi_outside_source_tree(
    tmp_path: Path,
) -> None:
    wheel = _build_wheel(tmp_path)
    installed = tmp_path / "installed"
    with zipfile.ZipFile(wheel) as archive:
        archive.extractall(installed)
    outside = tmp_path / "outside"
    outside.mkdir()
    runtime = tmp_path / "runtime"
    code = """
import json
from pathlib import Path
from encode_pipeline.api.main import create_app
from encode_pipeline.services.defaults import create_default_workflow_registry

registry = create_default_workflow_registry(environ={})
metadata = registry.list_metadata()
assert [item.workflow_id for item in metadata] == [
    "encode-style-chipseq-cuttag-atac-mnase",
    "bulk-rnaseq",
]
runtime = Path(__import__("os").environ["RELEASE_RUNTIME_ROOT"])
app = create_app(
    database_url=f"sqlite:///{runtime / 'platform.db'}",
    workspace_root=runtime / "workspaces",
)
try:
    schema = app.openapi()
    assert schema["info"]["title"] == "HelixWeave API"
    assert schema["info"]["version"] == "0.3.0"
    assert "/api/v1/workflows/" in schema["paths"]
finally:
    app.state.run_queue.close()
    app.state.persistence.close()
print(json.dumps([item.workflow_id for item in metadata]))
"""
    environment = {
        **os.environ,
        "PYTHONPATH": str(installed),
        "PYTHONNOUSERSITE": "1",
        "PYTHONDONTWRITEBYTECODE": "1",
        "RELEASE_RUNTIME_ROOT": str(runtime),
    }
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=outside,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert json.loads(completed.stdout) == [
        "encode-style-chipseq-cuttag-atac-mnase",
        "bulk-rnaseq",
    ]
