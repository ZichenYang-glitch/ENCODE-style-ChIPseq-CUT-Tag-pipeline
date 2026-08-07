"""Release identity and installed-distribution contracts for HelixWeave."""

from __future__ import annotations

from datetime import date
from email.parser import Parser
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tomllib
import venv
import zipfile

import pytest
import yaml

from encode_pipeline import __version__
from encode_pipeline.api.main import create_app


REPO_ROOT = Path(__file__).resolve().parents[2]
RELEASE_VERSION = "0.3.0"
RELEASE_DATE = date(2026, 7, 25)
EXPECTED_CONSOLE_SCRIPTS = {
    "encode-dag": "encode_pipeline.cli.dag:main",
    "encode-manifest": "encode_pipeline.cli.manifest:main",
    "encode-validate": "encode_pipeline.cli.validate:main",
    "encode-worker": "encode_pipeline.workers.cli:main",
    "helixweave": "encode_pipeline.cli.app:main",
}
CHECKOUT_BOOTSTRAP = REPO_ROOT / "scripts" / "checkout_bootstrap.py"


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
    wheels = tuple(wheel_root.glob("helixweave-*.whl"))
    assert len(wheels) == 1
    return wheels[0]


def _build_sdist(tmp_path: Path) -> Path:
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
            "--outdir",
            str(distribution_root),
            str(build_root),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    sdists = tuple(distribution_root.glob("helixweave-*.tar.gz"))
    assert len(sdists) == 1
    return sdists[0]


def _install_in_clean_room(tmp_path: Path, artifact: Path) -> Path:
    environment_root = tmp_path / "venv"
    venv.EnvBuilder(with_pip=False, system_site_packages=True).create(environment_root)
    python = environment_root / (
        "Scripts/python.exe" if os.name == "nt" else "bin/python"
    )
    environment = {
        key: value
        for key, value in os.environ.items()
        if key not in {"PYTHONPATH", "PYTHONHOME"}
    }
    installed = subprocess.run(
        [
            str(python),
            "-m",
            "pip",
            "install",
            "--ignore-installed",
            "--no-index",
            "--no-deps",
            "--no-build-isolation",
            str(artifact),
        ],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert installed.returncode == 0, installed.stderr

    configuration = environment_root / "pyvenv.cfg"
    content = configuration.read_text(encoding="utf-8")
    assert "include-system-site-packages = true" in content
    configuration.write_text(
        content.replace(
            "include-system-site-packages = true",
            "include-system-site-packages = false",
        ),
        encoding="utf-8",
    )
    return python


@pytest.mark.parametrize("artifact_kind", ("wheel", "sdist"))
def test_clean_room_artifact_provenance_mode(
    tmp_path: Path,
    artifact_kind: str,
) -> None:
    case_root = tmp_path / artifact_kind
    case_root.mkdir()
    artifact = (
        _build_wheel(case_root) if artifact_kind == "wheel" else _build_sdist(case_root)
    )
    python = _install_in_clean_room(case_root, artifact)

    provenance = subprocess.run(
        [
            str(python),
            "-I",
            "-S",
            str(CHECKOUT_BOOTSTRAP),
            "--repository-root",
            str(REPO_ROOT),
            "installed-artifact",
        ],
        cwd=case_root,
        capture_output=True,
        text=True,
        check=False,
    )
    imported = subprocess.run(
        [
            str(python),
            "-c",
            "import encode_pipeline; print(encode_pipeline.__version__)",
        ],
        cwd=case_root,
        capture_output=True,
        text=True,
        check=False,
    )

    assert provenance.returncode == 0, provenance.stderr
    assert provenance.stdout == ""
    assert provenance.stderr == ""
    assert imported.returncode == 0, imported.stderr
    assert imported.stdout.strip() == RELEASE_VERSION


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
    changelog = (REPO_ROOT / "CHANGELOG.md").read_text(encoding="utf-8")
    citation = yaml.safe_load((REPO_ROOT / "CITATION.cff").read_text(encoding="utf-8"))
    app = create_app(
        database_url=f"sqlite:///{tmp_path / 'release-identity.db'}",
        workspace_root=tmp_path / "workspaces",
        project_root=REPO_ROOT,
    )
    try:
        assert pyproject["project"]["name"] == "helixweave"
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
        assert citation["date-released"] == RELEASE_DATE
        assert "## [Unreleased]\n" in changelog
        assert f"## [{RELEASE_VERSION}] - {RELEASE_DATE.isoformat()}\n" in changelog
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

        assert metadata["Name"] == "helixweave"
        assert metadata["Version"] == RELEASE_VERSION
        assert metadata["Summary"].startswith("HelixWeave")
        assert metadata["Description-Content-Type"] == "text/markdown"
        assert "\n# HelixWeave\n" in metadata_text
        for name, target in EXPECTED_CONSOLE_SCRIPTS.items():
            assert f"{name} = {target}" in entry_points
        assert "encode_pipeline/__main__.py" in names
        assert "encode_pipeline/cli/app.py" in names
        assert "encode_pipeline/cli/local_platform.py" in names
        assert "encode_pipeline/cli/results_visibility_fixture.py" in names
        assert "encode_pipeline/artifacts/artifact-inventory.yaml" in names
        assert "encode_pipeline/persistence/alembic/script.py.mako" in names
        assert (
            "encode_pipeline/persistence/alembic/versions/"
            "20260726_09_project_sample_registry.py"
        ) in names
        assert (
            "encode_pipeline/persistence/alembic/versions/20260726_10_input_registry.py"
        ) in names
        assert (
            "encode_pipeline/persistence/alembic/versions/"
            "20260803_11_reference_profiles.py"
        ) in names
        assert (
            "encode_pipeline/persistence/alembic/versions/"
            "20260807_12_artifact_publications.py"
        ) in names
        assert (
            sum(
                "/persistence/alembic/versions/" in name and name.endswith(".py")
                for name in names
            )
            == 13
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
        "helixweave-0.3.0-py3-none-any.whl",
        "helixweave-0.3.0.tar.gz",
    }
    with tarfile.open(distribution_root / "helixweave-0.3.0.tar.gz", "r:gz") as archive:
        names = archive.getnames()
    assert "helixweave-0.3.0/README.md" in names
    assert (
        "helixweave-0.3.0/src/encode_pipeline/artifacts/artifact-inventory.yaml"
    ) in names
    assert not any(name.startswith("helixweave-0.3.0/test/") for name in names)


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
