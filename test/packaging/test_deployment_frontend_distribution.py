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
    (source / "scripts").mkdir()
    shutil.copy2(
        REPOSITORY_ROOT / "scripts" / "bootstrap_helixweave_operator.py",
        source / "scripts" / "bootstrap_helixweave_operator.py",
    )
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
        timeout=20,
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


def test_extracted_wheel_serves_frontend_and_api_without_node_or_checkout(
    tmp_path: Path,
    frontend_distributions: tuple[Path, Path],
) -> None:
    wheel, _sdist = frontend_distributions
    installed = tmp_path / "installed"
    with zipfile.ZipFile(wheel) as archive:
        archive.extractall(installed)
    outside = tmp_path / "outside"
    outside.mkdir()
    runtime = (
        tmp_path
        / "runtimes"
        / "encode"
        / ("sha256-" + "a" * 64)
        / "payload"
        / "contracts"
        / "encode-runtime"
    )
    shutil.copytree(REPOSITORY_ROOT / "workflow", runtime / "workflow")
    shutil.copytree(
        REPOSITORY_ROOT / "profiles" / "default",
        runtime / "profiles" / "default",
    )
    shutil.copytree(REPOSITORY_ROOT / "scripts", runtime / "scripts")
    (runtime / "docs" / "architecture").mkdir(parents=True)
    shutil.copy2(
        REPOSITORY_ROOT / "docs" / "architecture" / "artifact-inventory.yaml",
        runtime / "docs" / "architecture" / "artifact-inventory.yaml",
    )
    state_root = tmp_path / "state"
    state_root.mkdir()
    database = state_root / "platform.db"
    from encode_pipeline.persistence.runtime import open_run_persistence

    prepared = open_run_persistence(f"sqlite:///{database}")
    prepared.close()
    code = f"""
import asyncio
import json
import os
from pathlib import Path
import sys

import httpx

repository = Path({str(REPOSITORY_ROOT)!r}).resolve()
installed = Path({str(installed)!r}).resolve()
checkout_entries = {{repository, repository / "src"}}
sys.path[:] = [
    str(installed),
    *[
        item
        for item in sys.path
        if item and Path(item).resolve() not in checkout_entries
    ],
]

import encode_pipeline
import encode_pipeline.api.main as api_main
from encode_pipeline.deployment.frontend import create_app
from encode_pipeline.frontend_assets import load_packaged_frontend_assets
from encode_pipeline.persistence.authentication import (
    SqlAlchemyAuthenticationRepository,
)
from encode_pipeline.persistence.database import (
    create_database_engine,
    create_session_factory,
)
from encode_pipeline.services.authentication_service import (
    AccountAdministrationService,
)

assert Path(encode_pipeline.__file__).resolve().is_relative_to(installed)
assert Path(api_main.__file__).resolve().is_relative_to(installed)
assets = load_packaged_frontend_assets()
engine = create_database_engine(os.environ["ENCODE_PIPELINE_DATABASE_URL"])
try:
    AccountAdministrationService(
        repository=SqlAlchemyAuthenticationRepository(create_session_factory(engine))
    ).bootstrap_initial_administrator(
        "packaging-admin",
        "packaging admin password",
    )
finally:
    engine.dispose()
app = create_app()

async def verify():
    transport = httpx.ASGITransport(app=app)
    async with httpx.AsyncClient(
        transport=transport,
        base_url="http://helixweave.invalid",
    ) as client:
        index = await client.get("/", headers={{"accept": "text/html"}})
        javascript_path = next(
            path for path in sorted(assets.content) if path.endswith(".js")
        )
        script = await client.get("/" + javascript_path)
        login = await client.post(
            "/api/v1/auth/login",
            json={{
                "username": "packaging-admin",
                "password": "packaging admin password",
            }},
        )
        assert login.status_code == 200
        api = await client.get(
            "/api/v1/workflows/encode-style-chipseq-cuttag-atac-mnase/schema",
            headers={{"accept": "application/json"}},
        )
    assert index.status_code == script.status_code == api.status_code == 200
    assert index.content == assets.content["index.html"]
    assert script.content == assets.content[javascript_path]
    assert api.json()["workflow_id"] == "encode-style-chipseq-cuttag-atac-mnase"
    assert api.json()["schema"]["schema_version"]

try:
    asyncio.run(verify())
finally:
    app.state.run_queue.close()
    app.state.persistence.close()

print(json.dumps({{
    "api_status": 200,
    "frontend_identity": assets.manifest.identity,
}}, sort_keys=True))
"""
    environment = {
        "ENCODE_PIPELINE_DATABASE_URL": f"sqlite:///{database}",
        "ENCODE_PIPELINE_REDIS_URL": "redis://127.0.0.1:1/15",
        "ENCODE_PIPELINE_WORKSPACE_ROOT": str(state_root / "workspaces"),
        "HELIXWEAVE_ACTIVE_API_CONTRACT_SHA256": (
            load_packaged_frontend_assets().manifest.api_contract_sha256
        ),
        "HELIXWEAVE_ENCODE_RUNTIME_ROOT": str(runtime),
        "PATH": "/node-and-npm-are-intentionally-unavailable",
        "PYTHONNOUSERSITE": "1",
        "PYTHONDONTWRITEBYTECODE": "1",
    }
    completed = subprocess.run(
        [sys.executable, "-I", "-c", code],
        cwd=outside,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )

    assert completed.returncode == 0, completed.stderr
    receipt = json.loads(completed.stdout)
    assert receipt == {
        "api_status": 200,
        "frontend_identity": load_packaged_frontend_assets().manifest.identity,
    }
