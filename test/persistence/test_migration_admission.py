"""Fail-closed admission tests for source-owned Alembic execution inventory."""

from __future__ import annotations

from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys

import pytest

from encode_pipeline.persistence import migration_admission as admission
from encode_pipeline.persistence import migrations as migration_runner
from encode_pipeline.persistence.migration_admission import (
    MigrationAdmissionError,
    admitted_migration_script_location,
    build_migration_execution_inventory,
    canonical_migration_execution_inventory_bytes,
    verify_migration_execution_inventory,
)
from encode_pipeline.persistence.migrations import upgrade_database


PROJECT_ROOT = Path(__file__).resolve().parents[2]
UNKNOWN_CHILD = """\
from __future__ import annotations

import os
from pathlib import Path

from alembic import op
import sqlalchemy as sa

Path(os.environ["HELIXWEAVE_MIGRATION_IMPORT_MARKER"]).write_text(
    "unknown revision imported",
    encoding="utf-8",
)

revision = "20990101_99"
down_revision = "20260726_10"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.execute("UPDATE admission_sentinel SET value = 'mutated'")
    op.create_table(
        "unknown_revision_mutation",
        sa.Column("value", sa.Text(), nullable=False),
    )


def downgrade() -> None:
    op.drop_table("unknown_revision_mutation")
"""


def test_migration_runner_does_not_expose_a_raw_config_builder() -> None:
    assert not hasattr(migration_runner, "alembic_config")


def test_maximum_linear_revision_graph_validates_without_recursion() -> None:
    revisions = []
    previous: str | None = None
    for ordinal in range(admission._MAXIMUM_REVISIONS):
        revision = f"{20000000 + ordinal:08d}_01"
        revisions.append(
            {
                "branch_labels": [],
                "depends_on": [],
                "down_revision": [] if previous is None else [previous],
                "revision": revision,
            }
        )
        previous = revision

    assert admission._validate_graph_entries(revisions) == (
        ("20000000_01",),
        (previous,),
    )


def _database_state(path: Path) -> dict[str, object]:
    auxiliary = {}
    for suffix in ("", "-journal", "-shm", "-wal"):
        candidate = Path(f"{path}{suffix}")
        auxiliary[suffix or "main"] = (
            hashlib.sha256(candidate.read_bytes()).hexdigest()
            if candidate.is_file()
            else None
        )
    with sqlite3.connect(path) as connection:
        schema = connection.execute(
            "SELECT type, name, tbl_name, sql FROM sqlite_master ORDER BY type, name"
        ).fetchall()
        sentinel = connection.execute(
            "SELECT value FROM admission_sentinel ORDER BY rowid"
        ).fetchall()
        [revision] = connection.execute(
            "SELECT version_num FROM alembic_version"
        ).fetchone()
    return {
        "files": auxiliary,
        "schema": schema,
        "sentinel": sentinel,
        "revision": revision,
    }


def _copy_source_with_unknown_child(destination: Path) -> Path:
    package_root = destination / "encode_pipeline"
    shutil.copytree(
        PROJECT_ROOT / "src/encode_pipeline",
        package_root,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    child = package_root / "persistence/alembic/versions/20990101_99_unknown_child.py"
    child.write_text(UNKNOWN_CHILD, encoding="utf-8")
    return destination


def _startup_probe(
    *,
    source_root: Path,
    database_path: Path,
    marker: Path,
    workspace: Path,
    entrypoint: str,
) -> dict[str, object]:
    code = """
import json
import os
from pathlib import Path

database_url = os.environ["HELIXWEAVE_TEST_DATABASE_URL"]
entrypoint = os.environ["HELIXWEAVE_TEST_ENTRYPOINT"]
workspace = Path(os.environ["HELIXWEAVE_TEST_WORKSPACE"])
try:
    if entrypoint == "api":
        from encode_pipeline.api.main import create_app

        resource = create_app(database_url=database_url, workspace_root=workspace)
        resource.state.run_queue.close()
        resource.state.persistence.close()
    elif entrypoint == "worker":
        from encode_pipeline.workers.runtime import open_worker_runtime
        from encode_pipeline.workers.settings import load_worker_settings

        settings = load_worker_settings(
            {
                "ENCODE_PIPELINE_DATABASE_URL": database_url,
                "ENCODE_PIPELINE_WORKSPACE_ROOT": str(workspace),
                "ENCODE_PIPELINE_REDIS_URL": "redis://127.0.0.1:1/15",
            }
        )
        resource = open_worker_runtime(settings)
        resource.close()
    else:
        raise AssertionError("unknown test entrypoint")
except Exception as error:
    print(
        json.dumps(
            {
                "status": "rejected",
                "reason_code": getattr(error, "reason_code", None),
                "message": str(error),
            },
            sort_keys=True,
        )
    )
else:
    print(json.dumps({"status": "accepted"}, sort_keys=True))
"""
    environment = {
        **os.environ,
        "PYTHONPATH": str(source_root),
        "PYTHONDONTWRITEBYTECODE": "1",
        "HELIXWEAVE_MIGRATION_IMPORT_MARKER": str(marker),
        "HELIXWEAVE_TEST_DATABASE_URL": f"sqlite:///{database_path}",
        "HELIXWEAVE_TEST_ENTRYPOINT": entrypoint,
        "HELIXWEAVE_TEST_WORKSPACE": str(workspace),
    }
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=source_root.parent,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr
    return json.loads(completed.stdout.splitlines()[-1])


@pytest.mark.parametrize("entrypoint", ("api", "worker"))
def test_unknown_child_is_rejected_before_api_or_worker_database_mutation(
    tmp_path: Path,
    entrypoint: str,
) -> None:
    database_path = tmp_path / f"{entrypoint}.db"
    upgrade_database(f"sqlite:///{database_path}")
    with sqlite3.connect(database_path) as connection:
        connection.execute("CREATE TABLE admission_sentinel (value TEXT NOT NULL)")
        connection.execute("INSERT INTO admission_sentinel (value) VALUES ('original')")
    before = _database_state(database_path)
    source_root = _copy_source_with_unknown_child(tmp_path / f"{entrypoint}-source")
    marker = tmp_path / f"{entrypoint}-revision-imported"

    result = _startup_probe(
        source_root=source_root,
        database_path=database_path,
        marker=marker,
        workspace=tmp_path / f"{entrypoint}-workspaces",
        entrypoint=entrypoint,
    )
    after = _database_state(database_path)

    assert {
        "startup": result,
        "revision_imported": marker.exists(),
        "database_unchanged": after == before,
    } == {
        "startup": {
            "status": "rejected",
            "reason_code": "MIGRATION_REVISION_UNKNOWN",
            "message": "migration execution admission failed "
            "[MIGRATION_REVISION_UNKNOWN]",
        },
        "revision_imported": False,
        "database_unchanged": True,
    }


@pytest.mark.parametrize("entrypoint", ("api", "worker"))
def test_unknown_child_rejection_does_not_create_database_parent(
    tmp_path: Path,
    entrypoint: str,
) -> None:
    source_root = _copy_source_with_unknown_child(
        tmp_path / f"{entrypoint}-fresh-source"
    )
    database_path = tmp_path / f"{entrypoint}-missing/state/platform.db"
    marker = tmp_path / f"{entrypoint}-fresh-revision-imported"

    result = _startup_probe(
        source_root=source_root,
        database_path=database_path,
        marker=marker,
        workspace=tmp_path / f"{entrypoint}-fresh-workspaces",
        entrypoint=entrypoint,
    )

    assert result == {
        "status": "rejected",
        "reason_code": "MIGRATION_REVISION_UNKNOWN",
        "message": (
            "migration execution admission failed [MIGRATION_REVISION_UNKNOWN]"
        ),
    }
    assert marker.exists() is False
    assert database_path.parent.exists() is False
    assert database_path.exists() is False


def _copy_persistence(destination: Path) -> Path:
    persistence_root = destination / "persistence"
    shutil.copytree(
        PROJECT_ROOT / "src/encode_pipeline/persistence",
        persistence_root,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    return persistence_root


@pytest.mark.parametrize(
    ("mutation", "reason_code"),
    (
        ("unknown", "MIGRATION_REVISION_UNKNOWN"),
        ("missing", "MIGRATION_REVISION_MISSING"),
        ("damaged", "MIGRATION_REVISION_DIGEST_MISMATCH"),
        ("symlink", "MIGRATION_REVISION_UNSAFE_TYPE"),
        ("directory", "MIGRATION_REVISION_UNSAFE_TYPE"),
    ),
)
def test_revision_inventory_rejects_unadmitted_file_states(
    tmp_path: Path,
    mutation: str,
    reason_code: str,
) -> None:
    persistence_root = _copy_persistence(tmp_path / mutation)
    revisions = persistence_root / "alembic/versions"
    selected = revisions / "20260711_02_run_execution_assignments.py"
    if mutation == "unknown":
        (revisions / "20990101_99_unknown.py").write_text(
            "revision = '20990101_99'\n"
            "down_revision = '20260726_10'\n"
            "branch_labels = None\n"
            "depends_on = None\n",
            encoding="utf-8",
        )
    elif mutation == "missing":
        selected.unlink()
    elif mutation == "damaged":
        selected.write_bytes(selected.read_bytes() + b"\n# damaged\n")
    elif mutation == "symlink":
        target = tmp_path / "external-revision.py"
        target.write_bytes(selected.read_bytes())
        selected.unlink()
        selected.symlink_to(target)
    elif mutation == "directory":
        selected.unlink()
        selected.mkdir()
    else:  # pragma: no cover - protects the parameter contract
        raise AssertionError("unknown mutation")

    with pytest.raises(MigrationAdmissionError) as raised:
        verify_migration_execution_inventory(persistence_root=persistence_root)

    assert raised.value.reason_code == reason_code
    assert str(raised.value) == (
        f"migration execution admission failed [{reason_code}]"
    )
    assert str(persistence_root) not in str(raised.value)


def _trusted_inventory_bytes(
    monkeypatch: pytest.MonkeyPatch,
    document: dict[str, object],
) -> bytes:
    content = canonical_migration_execution_inventory_bytes(document)
    monkeypatch.setattr(
        admission,
        "MIGRATION_EXECUTION_INVENTORY_SIZE_BYTES",
        len(content),
    )
    monkeypatch.setattr(
        admission,
        "MIGRATION_EXECUTION_INVENTORY_SHA256",
        hashlib.sha256(content).hexdigest(),
    )
    return content


def _with_inventory_identity(document: dict[str, object]) -> dict[str, object]:
    payload = {
        key: value for key, value in document.items() if key != "inventory_sha256"
    }
    return {
        **payload,
        "inventory_sha256": hashlib.sha256(
            json.dumps(
                payload,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            ).encode("utf-8")
        ).hexdigest(),
    }


@pytest.mark.parametrize(
    "mutation",
    ("duplicate", "orphan", "cycle", "duplicate_branch_label"),
)
def test_inventory_rejects_duplicate_or_invalid_revision_graph(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
) -> None:
    persistence_root = _copy_persistence(tmp_path / mutation)
    document = json.loads(
        (persistence_root / admission.MIGRATION_EXECUTION_INVENTORY_FILE).read_text(
            encoding="utf-8"
        )
    )
    revisions = document["revisions"]
    if mutation == "duplicate":
        revisions[1]["revision"] = revisions[0]["revision"]
        expected = "MIGRATION_REVISION_DUPLICATE"
    elif mutation == "orphan":
        revisions[-1]["down_revision"] = ["missing_parent"]
        expected = "MIGRATION_REVISION_GRAPH_INVALID"
    elif mutation == "cycle":
        revisions[0]["down_revision"] = [revisions[-1]["revision"]]
        expected = "MIGRATION_REVISION_GRAPH_INVALID"
    elif mutation == "duplicate_branch_label":
        revisions[0]["branch_labels"] = ["shared"]
        revisions[1]["branch_labels"] = ["shared"]
        expected = "MIGRATION_REVISION_GRAPH_INVALID"
    else:  # pragma: no cover - protects the parameter contract
        raise AssertionError("unknown mutation")
    document = _with_inventory_identity(document)
    content = _trusted_inventory_bytes(monkeypatch, document)

    with pytest.raises(MigrationAdmissionError) as raised:
        verify_migration_execution_inventory(
            persistence_root=persistence_root,
            inventory_bytes=content,
        )

    assert raised.value.reason_code == expected


def test_inventory_contract_has_an_independent_source_anchor(
    tmp_path: Path,
) -> None:
    persistence_root = _copy_persistence(tmp_path / "contract-tamper")
    child = persistence_root / "alembic/versions/20990101_99_review_catalog.py"
    child.write_text(
        "revision = '20990101_99'\n"
        "down_revision = '20260726_10'\n"
        "branch_labels = None\n"
        "depends_on = None\n",
        encoding="utf-8",
    )
    updated = build_migration_execution_inventory(persistence_root)
    updated_bytes = canonical_migration_execution_inventory_bytes(updated)
    (persistence_root / admission.MIGRATION_EXECUTION_INVENTORY_FILE).write_bytes(
        updated_bytes
    )

    with pytest.raises(MigrationAdmissionError) as raised:
        verify_migration_execution_inventory(persistence_root=persistence_root)

    assert raised.value.reason_code == "MIGRATION_EXECUTION_INVENTORY_INVALID"
    assert hashlib.sha256(updated_bytes).hexdigest() != (
        admission.MIGRATION_EXECUTION_INVENTORY_SHA256
    )


def test_reviewed_unrelated_revision_executes_after_explicit_inventory_anchor_update(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    persistence_root = _copy_persistence(tmp_path / "reviewed-unrelated")
    child = persistence_root / "alembic/versions/20990101_99_review_catalog.py"
    child.write_text(
        "from alembic import op\n"
        "import sqlalchemy as sa\n"
        "revision = '20990101_99'\n"
        "down_revision = '20260726_10'\n"
        "branch_labels = None\n"
        "depends_on = None\n"
        "def upgrade():\n"
        "    op.create_table(\n"
        "        'review_catalog_entries',\n"
        "        sa.Column('entry_id', sa.Text(), nullable=False),\n"
        "        sa.PrimaryKeyConstraint('entry_id'),\n"
        "    )\n"
        "def downgrade():\n"
        "    op.drop_table('review_catalog_entries')\n",
        encoding="utf-8",
    )
    updated = canonical_migration_execution_inventory_bytes(
        build_migration_execution_inventory(persistence_root)
    )
    (persistence_root / admission.MIGRATION_EXECUTION_INVENTORY_FILE).write_bytes(
        updated
    )
    monkeypatch.setattr(
        admission,
        "MIGRATION_EXECUTION_INVENTORY_SIZE_BYTES",
        len(updated),
    )
    monkeypatch.setattr(
        admission,
        "MIGRATION_EXECUTION_INVENTORY_SHA256",
        hashlib.sha256(updated).hexdigest(),
    )

    @contextmanager
    def reviewed_snapshot():
        with admitted_migration_script_location(
            persistence_root=persistence_root
        ) as captured:
            yield captured

    monkeypatch.setattr(
        migration_runner,
        "admitted_migration_script_location",
        reviewed_snapshot,
    )
    database_path = tmp_path / "reviewed-unrelated.db"

    migration_runner.upgrade_database(f"sqlite:///{database_path}")

    with sqlite3.connect(database_path) as connection:
        [revision] = connection.execute(
            "SELECT version_num FROM alembic_version"
        ).fetchone()
        [table_count] = connection.execute(
            "SELECT count(*) FROM sqlite_master "
            "WHERE type='table' AND name='review_catalog_entries'"
        ).fetchone()
    assert revision == "20990101_99"
    assert table_count == 1


def test_snapshot_materialization_failure_is_stable_and_public_safe(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    private_path = tmp_path / "private-snapshot-detail"

    def fail_snapshot(_path: Path, _content: bytes) -> None:
        raise OSError(f"cannot write {private_path}")

    monkeypatch.setattr(admission, "_write_snapshot_file", fail_snapshot)

    with pytest.raises(MigrationAdmissionError) as raised:
        with admitted_migration_script_location():
            raise AssertionError("snapshot failure must occur before yielding")

    assert raised.value.reason_code == "MIGRATION_EXECUTION_SNAPSHOT_FAILED"
    assert str(raised.value) == (
        "migration execution admission failed [MIGRATION_EXECUTION_SNAPSHOT_FAILED]"
    )
    assert str(private_path) not in str(raised.value)


def test_percent_in_snapshot_parent_is_a_legal_private_path(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    private_temporary_root = tmp_path / "private%snapshot"
    private_temporary_root.mkdir()
    monkeypatch.setattr(
        admission.tempfile,
        "tempdir",
        str(private_temporary_root),
    )
    database_path = tmp_path / "percent-snapshot.db"

    upgrade_database(f"sqlite:///{database_path}")

    with sqlite3.connect(database_path) as connection:
        [revision] = connection.execute(
            "SELECT version_num FROM alembic_version"
        ).fetchone()
    assert revision == "20260726_10"


def test_runtime_metadata_rewrite_is_rejected_before_database_path_mutation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    persistence_root = _copy_persistence(tmp_path / "runtime-rewrite")
    selected = persistence_root / "alembic/versions/20260726_10_input_registry.py"
    selected.write_bytes(
        selected.read_bytes() + b"\nif True:\n" + b"    down_revision = None\n"
    )
    updated = canonical_migration_execution_inventory_bytes(
        build_migration_execution_inventory(persistence_root)
    )
    (persistence_root / admission.MIGRATION_EXECUTION_INVENTORY_FILE).write_bytes(
        updated
    )
    monkeypatch.setattr(
        admission,
        "MIGRATION_EXECUTION_INVENTORY_SIZE_BYTES",
        len(updated),
    )
    monkeypatch.setattr(
        admission,
        "MIGRATION_EXECUTION_INVENTORY_SHA256",
        hashlib.sha256(updated).hexdigest(),
    )

    @contextmanager
    def reviewed_snapshot():
        with admitted_migration_script_location(
            persistence_root=persistence_root
        ) as captured:
            yield captured

    monkeypatch.setattr(
        migration_runner,
        "admitted_migration_script_location",
        reviewed_snapshot,
    )
    database_path = tmp_path / "missing/private/platform.db"

    with pytest.raises(MigrationAdmissionError) as raised:
        migration_runner.upgrade_database(f"sqlite:///{database_path}")

    assert raised.value.reason_code == "MIGRATION_REVISION_GRAPH_INVALID"
    assert str(raised.value) == (
        "migration execution admission failed [MIGRATION_REVISION_GRAPH_INVALID]"
    )
    assert database_path.parent.exists() is False
    assert database_path.exists() is False


def test_validated_revision_modules_are_not_reloaded_for_execution(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    persistence_root = _copy_persistence(tmp_path / "single-load")
    marker = tmp_path / "revision-import-count"
    monkeypatch.setenv("HELIXWEAVE_MIGRATION_IMPORT_MARKER", str(marker))
    selected = persistence_root / "alembic/versions/20260726_10_input_registry.py"
    selected.write_bytes(
        selected.read_bytes()
        + b"\nimport os as _test_os\n"
        + b"from pathlib import Path as _TestPath\n"
        + b"_test_marker = _TestPath("
        + b"_test_os.environ['HELIXWEAVE_MIGRATION_IMPORT_MARKER'])\n"
        + b"_test_count = int(_test_marker.read_text() "
        + b"if _test_marker.exists() else '0') + 1\n"
        + b"_test_marker.write_text(str(_test_count))\n"
        + b"if _test_count > 1:\n"
        + b"    down_revision = None\n"
    )
    updated = canonical_migration_execution_inventory_bytes(
        build_migration_execution_inventory(persistence_root)
    )
    (persistence_root / admission.MIGRATION_EXECUTION_INVENTORY_FILE).write_bytes(
        updated
    )
    monkeypatch.setattr(
        admission,
        "MIGRATION_EXECUTION_INVENTORY_SIZE_BYTES",
        len(updated),
    )
    monkeypatch.setattr(
        admission,
        "MIGRATION_EXECUTION_INVENTORY_SHA256",
        hashlib.sha256(updated).hexdigest(),
    )

    @contextmanager
    def reviewed_snapshot():
        with admitted_migration_script_location(
            persistence_root=persistence_root
        ) as captured:
            yield captured

    monkeypatch.setattr(
        migration_runner,
        "admitted_migration_script_location",
        reviewed_snapshot,
    )
    database_path = tmp_path / "single-load.db"

    migration_runner.upgrade_database(f"sqlite:///{database_path}")

    assert marker.read_text(encoding="utf-8") == "1"
    with sqlite3.connect(database_path) as connection:
        [revision] = connection.execute(
            "SELECT version_num FROM alembic_version"
        ).fetchone()
    assert revision == "20260726_10"


def test_inventory_generator_is_byte_identical_and_does_not_import_revisions(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    persistence_root = _copy_persistence(tmp_path / "generator")
    marker = tmp_path / "generator-import-marker"
    monkeypatch.setenv("HELIXWEAVE_MIGRATION_IMPORT_MARKER", str(marker))
    child = persistence_root / "alembic/versions/20990101_99_review_catalog.py"
    child.write_text(UNKNOWN_CHILD, encoding="utf-8")

    first = canonical_migration_execution_inventory_bytes(
        build_migration_execution_inventory(persistence_root)
    )
    second = canonical_migration_execution_inventory_bytes(
        build_migration_execution_inventory(persistence_root)
    )

    assert first == second
    assert marker.exists() is False


def test_verified_snapshot_excludes_later_source_replacement(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    persistence_root = _copy_persistence(tmp_path / "snapshot-source")
    marker = tmp_path / "snapshot-import-marker"
    monkeypatch.setenv("HELIXWEAVE_MIGRATION_IMPORT_MARKER", str(marker))
    database_path = tmp_path / "snapshot.db"
    verified_heads: tuple[str, ...] = ()

    @contextmanager
    def admitted_then_replace_source():
        nonlocal verified_heads
        with admitted_migration_script_location(
            persistence_root=persistence_root
        ) as captured:
            verified_heads = captured[1].heads
            child = persistence_root / "alembic/versions/20990101_99_unknown_child.py"
            child.write_text(UNKNOWN_CHILD, encoding="utf-8")
            selected = (
                persistence_root
                / "alembic/versions/20260711_02_run_execution_assignments.py"
            )
            selected.write_bytes(
                selected.read_bytes()
                + b"\nraise RuntimeError('replaced after admission')\n"
            )
            yield captured

    monkeypatch.setattr(
        migration_runner,
        "admitted_migration_script_location",
        admitted_then_replace_source,
    )
    migration_runner.upgrade_database(f"sqlite:///{database_path}")

    with sqlite3.connect(database_path) as connection:
        [revision] = connection.execute(
            "SELECT version_num FROM alembic_version"
        ).fetchone()
        unknown_table = connection.execute(
            "SELECT count(*) FROM sqlite_master "
            "WHERE type='table' AND name='unknown_revision_mutation'"
        ).fetchone()[0]
    assert revision == verified_heads[0] == "20260726_10"
    assert unknown_table == 0
    assert marker.exists() is False
