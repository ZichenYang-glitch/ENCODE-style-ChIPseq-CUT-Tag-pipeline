"""Fail-closed admission tests for source-owned Alembic execution inventory."""

from __future__ import annotations

from contextlib import contextmanager
import errno
import hashlib
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
import sys
from types import SimpleNamespace

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


@pytest.mark.parametrize(
    "revisions",
    (
        [],
        [{"revision": "invalid/revision"}],
        [{"revision": "r1"}, {"revision": "r1"}],
        [
            {
                "revision": "r1",
                "down_revision": [],
                "branch_labels": [],
                "depends_on": [],
            },
            {
                "revision": "r2",
                "down_revision": ["r1"],
                "branch_labels": [],
                "depends_on": [],
            },
            {
                "revision": "r3",
                "down_revision": ["r1"],
                "branch_labels": [],
                "depends_on": [],
            },
        ],
    ),
    ids=("empty", "invalid-revision", "duplicate-revision", "multiple-heads"),
)
def test_graph_validator_rejects_empty_ambiguous_or_invalid_graphs(
    revisions: list[dict[str, object]],
) -> None:
    with pytest.raises(admission._AdmissionFailure) as raised:
        admission._validate_graph_entries(revisions)

    expected = (
        "MIGRATION_REVISION_DUPLICATE"
        if len(revisions) == 2
        else "MIGRATION_REVISION_GRAPH_INVALID"
    )
    assert raised.value.reason_code == expected


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


def _copy_source_with_fifo_inventory(destination: Path) -> Path:
    package_root = destination / "encode_pipeline"
    shutil.copytree(
        PROJECT_ROOT / "src/encode_pipeline",
        package_root,
        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"),
    )
    inventory = (
        package_root / "persistence" / admission.MIGRATION_EXECUTION_INVENTORY_FILE
    )
    inventory.unlink()
    os.mkfifo(inventory)
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
        timeout=3,
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


@pytest.mark.parametrize("entrypoint", ("api", "worker"))
def test_fifo_inventory_is_rejected_before_api_or_worker_database_creation(
    tmp_path: Path,
    entrypoint: str,
) -> None:
    source_root = _copy_source_with_fifo_inventory(
        tmp_path / f"{entrypoint}-fifo-source"
    )
    database_path = tmp_path / f"{entrypoint}-fifo-missing/state/platform.db"
    marker = tmp_path / f"{entrypoint}-fifo-revision-imported"

    result = _startup_probe(
        source_root=source_root,
        database_path=database_path,
        marker=marker,
        workspace=tmp_path / f"{entrypoint}-fifo-workspaces",
        entrypoint=entrypoint,
    )

    assert result == {
        "status": "rejected",
        "reason_code": "MIGRATION_EXECUTION_INVENTORY_INVALID",
        "message": (
            "migration execution admission failed "
            "[MIGRATION_EXECUTION_INVENTORY_INVALID]"
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


@pytest.mark.parametrize("mutation", ("missing", "symlink", "hardlink", "directory"))
def test_inventory_anchor_rejects_missing_or_nonregular_source(
    tmp_path: Path,
    mutation: str,
) -> None:
    persistence_root = _copy_persistence(tmp_path / f"inventory-{mutation}")
    inventory = persistence_root / admission.MIGRATION_EXECUTION_INVENTORY_FILE
    if mutation == "missing":
        inventory.unlink()
    elif mutation == "symlink":
        target = tmp_path / "external-inventory.json"
        target.write_bytes(inventory.read_bytes())
        inventory.unlink()
        inventory.symlink_to(target)
    elif mutation == "hardlink":
        os.link(inventory, tmp_path / "second-inventory-link.json")
    elif mutation == "directory":
        inventory.unlink()
        inventory.mkdir()
    else:  # pragma: no cover - protects the parameter contract
        raise AssertionError("unknown mutation")

    with pytest.raises(MigrationAdmissionError) as raised:
        verify_migration_execution_inventory(persistence_root=persistence_root)

    assert raised.value.reason_code == "MIGRATION_EXECUTION_INVENTORY_INVALID"
    assert str(raised.value) == (
        "migration execution admission failed [MIGRATION_EXECUTION_INVENTORY_INVALID]"
    )
    assert str(persistence_root) not in str(raised.value)


def _assert_descriptor_is_closed(descriptor: int) -> None:
    with pytest.raises(OSError) as raised:
        os.fstat(descriptor)
    assert raised.value.errno == errno.EBADF


def test_open_directory_at_closes_descriptor_when_fstat_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    parent_descriptor = os.open(tmp_path, admission._directory_flags())
    opened: list[int] = []
    real_open = os.open

    def tracked_open(*args, **kwargs):
        descriptor = real_open(*args, **kwargs)
        opened.append(descriptor)
        return descriptor

    try:
        with monkeypatch.context() as patch:
            patch.setattr(admission.os, "open", tracked_open)
            patch.setattr(
                admission.os,
                "fstat",
                lambda _descriptor: (_ for _ in ()).throw(OSError("fstat failed")),
            )
            with pytest.raises(admission._AdmissionFailure) as raised:
                admission._open_directory_at(
                    parent_descriptor,
                    ".",
                    reason_code="MIGRATION_EXECUTION_INVENTORY_INVALID",
                )
        assert raised.value.reason_code == "MIGRATION_EXECUTION_INVENTORY_INVALID"
        assert len(opened) == 1
        _assert_descriptor_is_closed(opened[0])
    finally:
        for descriptor in opened:
            try:
                os.close(descriptor)
            except OSError:
                pass
        os.close(parent_descriptor)


def test_open_absolute_directory_closes_unowned_descriptor_when_fstat_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    opened: list[int] = []
    real_open = os.open
    real_fstat = os.fstat

    def tracked_open(*args, **kwargs):
        descriptor = real_open(*args, **kwargs)
        opened.append(descriptor)
        return descriptor

    def fail_first_component(descriptor: int):
        if len(opened) >= 2 and descriptor == opened[-1]:
            raise OSError("fstat failed")
        return real_fstat(descriptor)

    with monkeypatch.context() as patch:
        patch.setattr(admission.os, "open", tracked_open)
        patch.setattr(admission.os, "fstat", fail_first_component)
        with pytest.raises(admission._AdmissionFailure) as raised:
            with admission._open_absolute_directory(tmp_path):
                raise AssertionError("fstat failure must prevent yielding")

    try:
        assert raised.value.reason_code == "MIGRATION_EXECUTION_INVENTORY_INVALID"
        assert len(opened) == 2
        for descriptor in opened:
            _assert_descriptor_is_closed(descriptor)
    finally:
        for descriptor in opened:
            try:
                os.close(descriptor)
            except OSError:
                pass


@pytest.mark.parametrize(
    ("source", "reason_code"),
    (
        (b"\xff", "MIGRATION_REVISION_GRAPH_INVALID"),
        (b"revision =", "MIGRATION_REVISION_GRAPH_INVALID"),
        (
            b"revision = 'r2'\ndown_revision = 'r1'\nbranch_labels = None\n",
            "MIGRATION_REVISION_GRAPH_INVALID",
        ),
        (
            b"revision = 'r2'\nrevision = 'r3'\n"
            b"down_revision = 'r1'\nbranch_labels = None\ndepends_on = None\n",
            "MIGRATION_REVISION_DUPLICATE",
        ),
        (
            b"revision = make_revision()\ndown_revision = 'r1'\n"
            b"branch_labels = None\ndepends_on = None\n",
            "MIGRATION_REVISION_GRAPH_INVALID",
        ),
        (
            b"revision = 'r3'\ndown_revision = ('r2', 'r1')\n"
            b"branch_labels = None\ndepends_on = None\n",
            "MIGRATION_REVISION_GRAPH_INVALID",
        ),
        (
            b"revision: str\ndown_revision = 'r1'\n"
            b"branch_labels = None\ndepends_on = None\n",
            "MIGRATION_REVISION_DUPLICATE",
        ),
        (
            b"revision = 'invalid/revision'\ndown_revision = None\n"
            b"branch_labels = None\ndepends_on = None\n",
            "MIGRATION_REVISION_GRAPH_INVALID",
        ),
        (
            b"revision = 'r2'\ndown_revision = 1\n"
            b"branch_labels = None\ndepends_on = None\n",
            "MIGRATION_REVISION_GRAPH_INVALID",
        ),
    ),
)
def test_revision_metadata_parser_rejects_ambiguous_or_dynamic_graph_values(
    source: bytes,
    reason_code: str,
) -> None:
    with pytest.raises(admission._AdmissionFailure) as raised:
        admission._parse_revision_metadata(source)

    assert raised.value.reason_code == reason_code


def test_revision_metadata_parser_accepts_literal_annotated_assignments() -> None:
    assert admission._parse_revision_metadata(
        b"revision: str = 'r2'\n"
        b"down_revision: str = 'r1'\n"
        b"branch_labels: tuple[str, ...] = ('reviewed',)\n"
        b"depends_on: list[str] = []\n"
    ) == {
        "revision": "r2",
        "down_revision": ("r1",),
        "branch_labels": ("reviewed",),
        "depends_on": (),
    }


@pytest.mark.parametrize("mutation", ("short_read", "appended", "identity"))
def test_regular_source_capture_rejects_descriptor_read_races(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
) -> None:
    source = tmp_path / "source.py"
    source.write_bytes(b"revision = 'r1'\n")
    directory_descriptor = os.open(tmp_path, admission._directory_flags())
    real_open = os.open
    real_read = os.read
    real_fstat = os.fstat
    opened: list[int] = []
    fstat_calls = 0

    def tracked_open(*args, **kwargs):
        descriptor = real_open(*args, **kwargs)
        opened.append(descriptor)
        return descriptor

    def raced_read(descriptor: int, maximum_bytes: int) -> bytes:
        if mutation == "short_read":
            return b""
        if mutation == "appended" and maximum_bytes == 1:
            return b"x"
        return real_read(descriptor, maximum_bytes)

    def raced_fstat(descriptor: int):
        nonlocal fstat_calls
        info = real_fstat(descriptor)
        fstat_calls += 1
        if mutation == "identity" and fstat_calls == 2:
            return SimpleNamespace(
                st_dev=info.st_dev,
                st_ino=info.st_ino + 1,
                st_mode=info.st_mode,
                st_nlink=info.st_nlink,
                st_size=info.st_size,
            )
        return info

    try:
        with monkeypatch.context() as patch:
            patch.setattr(admission.os, "open", tracked_open)
            patch.setattr(admission.os, "read", raced_read)
            patch.setattr(admission.os, "fstat", raced_fstat)
            with pytest.raises(admission._AdmissionFailure) as raised:
                admission._read_regular_at(
                    directory_descriptor,
                    source.name,
                    maximum_bytes=1024,
                    invalid_reason="MIGRATION_REVISION_UNSAFE_TYPE",
                )
        assert raised.value.reason_code == "MIGRATION_REVISION_UNSAFE_TYPE"
        assert len(opened) == 1
        _assert_descriptor_is_closed(opened[0])
    finally:
        for descriptor in opened:
            try:
                os.close(descriptor)
            except OSError:
                pass
        os.close(directory_descriptor)


def test_regular_source_capture_requires_nonblocking_nofollow_open(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "source.py"
    source.write_bytes(b"revision = 'r1'\n")
    directory_descriptor = os.open(tmp_path, admission._directory_flags())
    real_open = os.open
    observed_flags: list[int] = []

    def capture_flags(*args, **kwargs):
        observed_flags.append(args[1])
        return real_open(*args, **kwargs)

    try:
        monkeypatch.setattr(admission.os, "open", capture_flags)
        assert (
            admission._read_regular_at(
                directory_descriptor,
                source.name,
                maximum_bytes=1024,
                invalid_reason="MIGRATION_REVISION_UNSAFE_TYPE",
            )
            == source.read_bytes()
        )
    finally:
        os.close(directory_descriptor)

    assert len(observed_flags) == 1
    assert observed_flags[0] & os.O_NONBLOCK
    assert observed_flags[0] & os.O_NOFOLLOW


def test_regular_source_capture_fails_closed_without_nonblocking_support(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    source = tmp_path / "source.py"
    source.write_bytes(b"revision = 'r1'\n")
    directory_descriptor = os.open(tmp_path, admission._directory_flags())
    monkeypatch.delattr(admission.os, "O_NONBLOCK")

    try:
        with pytest.raises(admission._AdmissionFailure) as raised:
            admission._read_regular_at(
                directory_descriptor,
                source.name,
                maximum_bytes=1024,
                invalid_reason="MIGRATION_REVISION_UNSAFE_TYPE",
            )
    finally:
        os.close(directory_descriptor)

    assert raised.value.reason_code == "MIGRATION_REVISION_UNSAFE_TYPE"


@pytest.mark.parametrize(
    "operation",
    (
        "relative_root",
        "invalid_file_name",
        "invalid_directory_name",
        "missing_directory",
    ),
)
def test_descriptor_helpers_reject_unanchored_or_invalid_paths(
    tmp_path: Path,
    operation: str,
) -> None:
    directory_descriptor = os.open(tmp_path, admission._directory_flags())
    try:
        with pytest.raises(admission._AdmissionFailure) as raised:
            if operation == "relative_root":
                with admission._open_absolute_directory(Path("relative")):
                    raise AssertionError("relative root must not be opened")
            elif operation == "invalid_file_name":
                admission._read_regular_at(
                    directory_descriptor,
                    "../escape.py",
                    maximum_bytes=1024,
                    invalid_reason="MIGRATION_REVISION_UNSAFE_TYPE",
                )
            elif operation == "invalid_directory_name":
                admission._open_directory_at(
                    directory_descriptor,
                    "../escape",
                    reason_code="MIGRATION_EXECUTION_INVENTORY_INVALID",
                )
            elif operation == "missing_directory":
                admission._open_directory_at(
                    directory_descriptor,
                    "missing",
                    reason_code="MIGRATION_EXECUTION_INVENTORY_INVALID",
                )
            else:  # pragma: no cover - protects the parameter contract
                raise AssertionError("unknown operation")
    finally:
        os.close(directory_descriptor)

    expected = (
        "MIGRATION_REVISION_UNSAFE_TYPE"
        if operation == "invalid_file_name"
        else "MIGRATION_EXECUTION_INVENTORY_INVALID"
    )
    assert raised.value.reason_code == expected


def test_snapshot_writer_closes_descriptor_after_short_write(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    destination = tmp_path / "snapshot.py"
    opened: list[int] = []
    real_open = os.open

    def tracked_open(*args, **kwargs):
        descriptor = real_open(*args, **kwargs)
        opened.append(descriptor)
        return descriptor

    with monkeypatch.context() as patch:
        patch.setattr(admission.os, "open", tracked_open)
        patch.setattr(admission.os, "write", lambda _descriptor, _content: 0)
        with pytest.raises(OSError, match="snapshot write failed"):
            admission._write_snapshot_file(destination, b"content")

    try:
        assert len(opened) == 1
        _assert_descriptor_is_closed(opened[0])
    finally:
        for descriptor in opened:
            try:
                os.close(descriptor)
            except OSError:
                pass


@pytest.mark.parametrize(
    ("failure", "reason_code"),
    (
        (
            MigrationAdmissionError("MIGRATION_REVISION_UNKNOWN"),
            "MIGRATION_REVISION_UNKNOWN",
        ),
        (ValueError("private parser detail"), "MIGRATION_EXECUTION_INVENTORY_INVALID"),
    ),
)
def test_verifier_public_boundary_preserves_only_stable_reason_codes(
    monkeypatch: pytest.MonkeyPatch,
    failure: Exception,
    reason_code: str,
) -> None:
    @contextmanager
    def fail_open(_path: Path):
        raise failure
        yield 0

    monkeypatch.setattr(admission, "_open_absolute_directory", fail_open)

    with pytest.raises(MigrationAdmissionError) as raised:
        verify_migration_execution_inventory()

    assert raised.value.reason_code == reason_code
    assert str(raised.value) == (
        f"migration execution admission failed [{reason_code}]"
    )
    assert "private parser detail" not in str(raised.value)


@pytest.mark.parametrize(
    "failure",
    (
        admission._AdmissionFailure("MIGRATION_REVISION_GRAPH_INVALID"),
        OSError("private source path"),
    ),
)
def test_inventory_builder_translates_internal_source_failures(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    failure: Exception,
) -> None:
    def fail_root(_value: object) -> Path:
        raise failure

    monkeypatch.setattr(admission, "_require_absolute_path", fail_root)

    with pytest.raises(
        ValueError,
        match="^migration execution inventory source is invalid$",
    ) as raised:
        build_migration_execution_inventory(tmp_path)

    assert str(tmp_path) not in str(raised.value)
    assert "private source path" not in str(raised.value)


@pytest.mark.parametrize(
    ("mutation", "reason_code"),
    (
        ("unknown", "MIGRATION_REVISION_UNKNOWN"),
        ("missing", "MIGRATION_REVISION_MISSING"),
        ("damaged", "MIGRATION_REVISION_DIGEST_MISMATCH"),
        ("symlink", "MIGRATION_REVISION_UNSAFE_TYPE"),
        ("hardlink", "MIGRATION_REVISION_UNSAFE_TYPE"),
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
    elif mutation == "hardlink":
        os.link(selected, tmp_path / "second-revision-link.py")
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
    return _trust_raw_inventory_bytes(monkeypatch, content)


def _trust_raw_inventory_bytes(
    monkeypatch: pytest.MonkeyPatch,
    content: bytes,
) -> bytes:
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
    ("mutation", "reason_code"),
    (
        ("duplicate_json_key", "MIGRATION_EXECUTION_INVENTORY_INVALID"),
        ("invalid_utf8", "MIGRATION_EXECUTION_INVENTORY_INVALID"),
        ("noncanonical_json", "MIGRATION_EXECUTION_INVENTORY_INVALID"),
        ("top_level_fields", "MIGRATION_EXECUTION_INVENTORY_INVALID"),
        ("identity", "MIGRATION_EXECUTION_INVENTORY_INVALID"),
        ("revision_count", "MIGRATION_EXECUTION_INVENTORY_INVALID"),
        ("revision_order", "MIGRATION_EXECUTION_INVENTORY_INVALID"),
        ("revision_fields", "MIGRATION_EXECUTION_INVENTORY_INVALID"),
        ("revision_id", "MIGRATION_EXECUTION_INVENTORY_INVALID"),
        ("revision_list", "MIGRATION_REVISION_GRAPH_INVALID"),
        ("environment_shape", "MIGRATION_EXECUTION_INVENTORY_INVALID"),
        ("environment_fields", "MIGRATION_EXECUTION_INVENTORY_INVALID"),
        ("environment_size", "MIGRATION_EXECUTION_INVENTORY_INVALID"),
        ("declared_head", "MIGRATION_REVISION_GRAPH_INVALID"),
    ),
)
def test_trusted_inventory_rejects_ambiguous_or_noncanonical_documents(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
    reason_code: str,
) -> None:
    persistence_root = _copy_persistence(tmp_path / mutation)
    document = json.loads(
        (persistence_root / admission.MIGRATION_EXECUTION_INVENTORY_FILE).read_text(
            encoding="utf-8"
        )
    )
    if mutation == "duplicate_json_key":
        canonical = canonical_migration_execution_inventory_bytes(document)
        content = canonical.replace(
            b"{",
            b'{"schema_version":"1.0.0",',
            1,
        )
    elif mutation == "invalid_utf8":
        content = b"\xff"
    elif mutation == "noncanonical_json":
        content = (json.dumps(document, indent=2) + "\n").encode("utf-8")
    else:
        recompute_identity = True
        if mutation == "top_level_fields":
            del document["scheme"]
        elif mutation == "identity":
            document["inventory_sha256"] = "0" * 64
            recompute_identity = False
        elif mutation == "revision_count":
            document["revision_count"] += 1
        elif mutation == "revision_order":
            document["revisions"] = list(reversed(document["revisions"]))
        elif mutation == "revision_fields":
            del document["revisions"][0]["depends_on"]
        elif mutation == "revision_id":
            document["revisions"][0]["revision"] = "invalid/revision"
        elif mutation == "revision_list":
            document["revisions"][0]["branch_labels"] = ["z", "a"]
        elif mutation == "environment_shape":
            document["environment"] = []
        elif mutation == "environment_fields":
            document["environment"]["unexpected"] = True
        elif mutation == "environment_size":
            document["environment"]["size_bytes"] = 0
        elif mutation == "declared_head":
            document["heads"] = [document["revisions"][0]["revision"]]
        else:  # pragma: no cover - protects the parameter contract
            raise AssertionError("unknown mutation")
        content = canonical_migration_execution_inventory_bytes(
            _with_inventory_identity(document) if recompute_identity else document
        )
    content = _trust_raw_inventory_bytes(monkeypatch, content)

    with pytest.raises(MigrationAdmissionError) as raised:
        verify_migration_execution_inventory(
            persistence_root=persistence_root,
            inventory_bytes=content,
        )

    assert raised.value.reason_code == reason_code
    assert str(raised.value) == (
        f"migration execution admission failed [{reason_code}]"
    )
    assert str(persistence_root) not in str(raised.value)


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


def test_declared_metadata_must_match_the_captured_revision_source(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    persistence_root = _copy_persistence(tmp_path / "metadata-mismatch")
    document = json.loads(
        (persistence_root / admission.MIGRATION_EXECUTION_INVENTORY_FILE).read_text(
            encoding="utf-8"
        )
    )
    document["revisions"][0]["branch_labels"] = ["reviewed_branch"]
    content = _trusted_inventory_bytes(
        monkeypatch,
        _with_inventory_identity(document),
    )

    with pytest.raises(MigrationAdmissionError) as raised:
        verify_migration_execution_inventory(
            persistence_root=persistence_root,
            inventory_bytes=content,
        )

    assert raised.value.reason_code == "MIGRATION_REVISION_GRAPH_INVALID"


def test_revision_disappearing_after_inventory_listing_is_stably_missing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    persistence_root = _copy_persistence(tmp_path / "read-race")
    selected = "20260711_02_run_execution_assignments.py"
    real_read = admission._read_regular_at

    def disappear_after_listing(
        directory_descriptor: int,
        name: str,
        *,
        maximum_bytes: int,
        invalid_reason: str,
    ) -> bytes:
        if name == selected:
            raise FileNotFoundError(name)
        return real_read(
            directory_descriptor,
            name,
            maximum_bytes=maximum_bytes,
            invalid_reason=invalid_reason,
        )

    monkeypatch.setattr(admission, "_read_regular_at", disappear_after_listing)

    with pytest.raises(MigrationAdmissionError) as raised:
        verify_migration_execution_inventory(persistence_root=persistence_root)

    assert raised.value.reason_code == "MIGRATION_REVISION_MISSING"


def test_revision_stat_race_is_stably_rejected_as_unsafe(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    persistence_root = _copy_persistence(tmp_path / "stat-race")
    real_stat = os.stat
    calls = 0

    def fail_first_revision_stat(*args, **kwargs):
        nonlocal calls
        if kwargs.get("dir_fd") is not None:
            calls += 1
            if calls == 1:
                raise FileNotFoundError(args[0])
        return real_stat(*args, **kwargs)

    monkeypatch.setattr(admission.os, "stat", fail_first_revision_stat)

    with pytest.raises(MigrationAdmissionError) as raised:
        verify_migration_execution_inventory(persistence_root=persistence_root)

    assert raised.value.reason_code == "MIGRATION_REVISION_UNSAFE_TYPE"


@pytest.mark.parametrize(
    ("mutation", "reason_code"),
    (
        ("list_error", "MIGRATION_EXECUTION_INVENTORY_INVALID"),
        ("empty", "MIGRATION_REVISION_GRAPH_INVALID"),
        ("casefold_duplicate", "MIGRATION_REVISION_DUPLICATE"),
    ),
)
def test_revision_inventory_listing_rejects_unreadable_empty_or_ambiguous_sets(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
    reason_code: str,
) -> None:
    versions = tmp_path / "versions"
    versions.mkdir()
    if mutation == "casefold_duplicate":
        (versions / "Reviewed.py").write_text("revision = 'r1'\n", encoding="utf-8")
        (versions / "reviewed.py").write_text("revision = 'r2'\n", encoding="utf-8")
    descriptor = os.open(versions, admission._directory_flags())
    if mutation == "list_error":
        monkeypatch.setattr(
            admission.os,
            "listdir",
            lambda _descriptor: (_ for _ in ()).throw(OSError("list failed")),
        )

    try:
        with pytest.raises(admission._AdmissionFailure) as raised:
            admission._revision_source_names(descriptor)
    finally:
        os.close(descriptor)

    assert raised.value.reason_code == reason_code


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
    created_snapshots: list[Path] = []
    real_temporary_directory = admission.tempfile.TemporaryDirectory

    def track_temporary_directory(*args, **kwargs):
        temporary = real_temporary_directory(*args, **kwargs)
        created_snapshots.append(Path(temporary.name))
        return temporary

    def fail_snapshot(_path: Path, _content: bytes) -> None:
        raise OSError(f"cannot write {private_path}")

    monkeypatch.setattr(
        admission.tempfile,
        "TemporaryDirectory",
        track_temporary_directory,
    )
    monkeypatch.setattr(admission, "_write_snapshot_file", fail_snapshot)

    with pytest.raises(MigrationAdmissionError) as raised:
        with admitted_migration_script_location():
            raise AssertionError("snapshot failure must occur before yielding")

    assert raised.value.reason_code == "MIGRATION_EXECUTION_SNAPSHOT_FAILED"
    assert str(raised.value) == (
        "migration execution admission failed [MIGRATION_EXECUTION_SNAPSHOT_FAILED]"
    )
    assert str(private_path) not in str(raised.value)
    assert len(created_snapshots) == 1
    assert created_snapshots[0].exists() is False


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


@pytest.mark.parametrize(
    "runtime_mutation",
    (
        b"\nif True:\n    down_revision = None\n",
        b"\nraise RuntimeError('private runtime import detail')\n",
    ),
    ids=("metadata-rewrite", "import-failure"),
)
def test_runtime_revision_failure_is_rejected_before_database_path_mutation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    runtime_mutation: bytes,
) -> None:
    persistence_root = _copy_persistence(tmp_path / "runtime-rewrite")
    selected = persistence_root / "alembic/versions/20260726_10_input_registry.py"
    selected.write_bytes(selected.read_bytes() + runtime_mutation)
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
