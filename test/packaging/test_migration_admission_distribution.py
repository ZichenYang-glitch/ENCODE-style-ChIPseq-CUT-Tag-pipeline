"""Installed-artifact tests for migration execution admission."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import venv

import pytest

from encode_pipeline.persistence.migration_admission import (
    MIGRATION_EXECUTION_INVENTORY_FILE,
    build_migration_execution_inventory,
    canonical_migration_execution_inventory_bytes,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
CHECKOUT_BOOTSTRAP = REPOSITORY_ROOT / "scripts/checkout_bootstrap.py"
UNRELATED_CHILD = """\
from __future__ import annotations

from alembic import op
import sqlalchemy as sa

revision = "20990101_99"
down_revision = "20260803_11"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "review_catalog_entries",
        sa.Column("entry_id", sa.Text(), nullable=False),
        sa.PrimaryKeyConstraint("entry_id"),
    )


def downgrade() -> None:
    op.drop_table("review_catalog_entries")
"""


def _tampered_source_tree(
    tmp_path: Path,
    *,
    synchronize_inventory: bool,
) -> Path:
    source = tmp_path / "source"
    source.mkdir()
    for filename in (
        "pyproject.toml",
        "README.md",
        "LICENSE",
        "MANIFEST.in",
    ):
        shutil.copy2(REPOSITORY_ROOT / filename, source / filename)
    shutil.copytree(
        REPOSITORY_ROOT / "src",
        source / "src",
        ignore=shutil.ignore_patterns("*.egg-info", "__pycache__", "*.pyc"),
    )
    persistence_root = source / "src/encode_pipeline/persistence"
    child = persistence_root / "alembic/versions/20990101_99_review_catalog.py"
    child.write_text(UNRELATED_CHILD, encoding="utf-8")
    if synchronize_inventory:
        document = build_migration_execution_inventory(persistence_root)
        (persistence_root / MIGRATION_EXECUTION_INVENTORY_FILE).write_bytes(
            canonical_migration_execution_inventory_bytes(document)
        )
        # Deliberately do not synchronize the reviewed size/SHA constants in
        # migration_admission.py. The resulting distribution RECORD is internally
        # valid, but runtime admission must still reject it.
    return source


def _build_artifact(
    tmp_path: Path,
    *,
    artifact_kind: str,
    synchronize_inventory: bool,
) -> Path:
    source = _tampered_source_tree(
        tmp_path,
        synchronize_inventory=synchronize_inventory,
    )
    destination = tmp_path / "dist"
    destination.mkdir()
    if artifact_kind == "wheel":
        build_python = _python_with_module("setuptools")
        command = [
            str(build_python),
            "-m",
            "pip",
            "wheel",
            "--no-index",
            "--no-deps",
            "--no-build-isolation",
            "--wheel-dir",
            str(destination),
            str(source),
        ]
        pattern = "helixweave-*.whl"
    elif artifact_kind == "sdist":
        build_python = _python_with_module("build")
        command = [
            str(build_python),
            "-m",
            "build",
            "--no-isolation",
            "--sdist",
            "--outdir",
            str(destination),
            str(source),
        ]
        pattern = "helixweave-*.tar.gz"
    else:  # pragma: no cover - protects the parameter contract
        raise AssertionError("unknown artifact kind")
    completed = subprocess.run(
        command,
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    artifacts = tuple(destination.glob(pattern))
    assert len(artifacts) == 1
    return artifacts[0]


def _python_with_module(module: str) -> Path:
    candidates = [
        Path(sys.executable),
        Path(sys.base_prefix) / ("python.exe" if os.name == "nt" else "bin/python"),
    ]
    for candidate in dict.fromkeys(candidates):
        available = subprocess.run(
            [str(candidate), "-c", f"import {module}"],
            capture_output=True,
            text=True,
            check=False,
        )
        if available.returncode == 0:
            return candidate
    raise AssertionError(f"test environment cannot import {module}")


def _install_artifact(
    tmp_path: Path,
    artifact: Path,
) -> tuple[Path, Path]:
    environment_root = tmp_path / "venv"
    venv.EnvBuilder(with_pip=False, system_site_packages=True).create(environment_root)
    python = environment_root / (
        "Scripts/python.exe" if os.name == "nt" else "bin/python"
    )
    environment = {
        key: value
        for key, value in os.environ.items()
        if key not in {"PYTHONHOME", "PYTHONPATH"}
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
    located = subprocess.run(
        [
            str(python),
            "-c",
            "import sysconfig; print(sysconfig.get_path('purelib'))",
        ],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert located.returncode == 0, located.stderr
    site_packages = Path(located.stdout.strip())
    assert site_packages.is_dir()

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
    return python, site_packages


@pytest.mark.parametrize("artifact_kind", ("wheel", "sdist"))
@pytest.mark.parametrize(
    ("synchronize_inventory", "expected_reason"),
    (
        (False, "MIGRATION_REVISION_UNKNOWN"),
        (True, "MIGRATION_EXECUTION_INVENTORY_INVALID"),
    ),
)
def test_self_consistent_record_cannot_admit_migration_tamper(
    tmp_path: Path,
    artifact_kind: str,
    synchronize_inventory: bool,
    expected_reason: str,
) -> None:
    artifact = _build_artifact(
        tmp_path,
        artifact_kind=artifact_kind,
        synchronize_inventory=synchronize_inventory,
    )
    python, site_packages = _install_artifact(tmp_path, artifact)
    environment = {
        key: value
        for key, value in os.environ.items()
        if key not in {"PYTHONHOME", "PYTHONPATH"}
    }
    environment["HELIXWEAVE_TEST_SITE_PACKAGES"] = str(site_packages)
    provenance = subprocess.run(
        [
            str(python),
            "-I",
            "-S",
            str(CHECKOUT_BOOTSTRAP),
            "--repository-root",
            str(REPOSITORY_ROOT),
            "installed-artifact",
        ],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    probe = r"""
import importlib.util
import json
import os
from pathlib import Path
import sys

path = (
    Path(os.environ["HELIXWEAVE_TEST_SITE_PACKAGES"])
    / "encode_pipeline/persistence/migration_admission.py"
)
spec = importlib.util.spec_from_file_location(
    "_installed_migration_admission",
    path,
)
module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = module
spec.loader.exec_module(module)
try:
    module.verify_migration_execution_inventory()
except module.MigrationAdmissionError as error:
    print(
        json.dumps(
            {
                "status": "rejected",
                "reason_code": error.reason_code,
                "message": str(error),
            },
            sort_keys=True,
        )
    )
else:
    print(json.dumps({"status": "accepted"}, sort_keys=True))
"""
    admitted = subprocess.run(
        [str(python), "-I", "-S", "-c", probe],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert provenance.returncode == 0, provenance.stderr
    assert provenance.stdout == ""
    assert provenance.stderr == ""
    assert admitted.returncode == 0, admitted.stderr
    assert json.loads(admitted.stdout) == {
        "status": "rejected",
        "reason_code": expected_reason,
        "message": (f"migration execution admission failed [{expected_reason}]"),
    }
