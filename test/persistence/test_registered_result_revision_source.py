"""Execute the registered-result migration from its source file.

Admission runs a verified copy, so shard coverage would not attribute those
statements to this revision. Pointing Alembic at the source tree does.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from alembic.config import Config
from alembic.script import ScriptDirectory
from sqlalchemy import create_engine, text

from alembic import command

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
ALEMBIC = REPOSITORY_ROOT / "src/encode_pipeline/persistence/alembic"
MIGRATION = ALEMBIC / "versions" / "20261119_39_registered_results.py"


def _config(database_url: str) -> Config:
    config = Config()
    config.set_main_option("script_location", str(ALEMBIC))
    config.set_main_option("sqlalchemy.url", database_url.replace("%", "%%"))
    return config


def _names(database_url: str) -> set[str]:
    engine = create_engine(database_url)
    with engine.connect() as connection:
        names = set(
            connection.execute(
                text("SELECT name FROM sqlite_master WHERE type = 'table'")
            ).scalars()
        )
    engine.dispose()
    return names


def test_source_revision_upgrades_and_downgrades_around_snapshot_bindings(
    tmp_path: Path,
) -> None:
    source = MIGRATION.read_text(encoding="utf-8")
    assert 'revision = "20261119_39"' in source
    assert 'down_revision = "20261119_32"' in source
    database_url = f"sqlite:///{tmp_path / 'revision-39.db'}"
    command.upgrade(_config(database_url), "20261119_39")
    names = _names(database_url)
    assert "registered_results" in names
    assert "catalog_adoption_receipts" in names
    assert "snapshot_experiment_bindings" in names
    command.downgrade(_config(database_url), "20261119_32")
    names = _names(database_url)
    assert "registered_results" not in names
    assert "catalog_adoption_receipts" in names
    assert "snapshot_experiment_bindings" in names
    script = ScriptDirectory.from_config(_config(database_url))
    assert list(script.get_heads()) == ["20261119_39"]
    registered = script.get_revision("20261119_39")
    proposals = script.get_revision("20261119_32")
    assert registered is not None and proposals is not None
    assert registered.down_revision == "20261119_32"
    assert proposals.down_revision == "20261119_18"
    command.downgrade(_config(database_url), "20261119_18")
    names = _names(database_url)
    assert "catalog_proposals" not in names
    assert "catalog_adoption_receipts" not in names
    assert "snapshot_experiment_bindings" in names


def test_registered_result_downgrade_refuses_existing_download_audit(
    tmp_path: Path,
) -> None:
    from encode_pipeline.persistence.migrations import downgrade_database

    database_url = f"sqlite:///{tmp_path / 'audit-downgrade.db'}"
    command.upgrade(_config(database_url), "20261119_39")
    engine = create_engine(database_url)
    with engine.begin() as connection:
        for index in range(3):
            connection.execute(
                text(
                    """
                    INSERT INTO security_audit_events (
                        event_id, occurred_at, action, outcome, actor_kind,
                        actor_user_id, resource_kind, resource_id, reason_code
                    ) VALUES (
                        :event_id, :now, 'registered_result.download', 'succeeded',
                        'local_operator', NULL, 'registered_result', :resource, NULL
                    )
                    """
                ),
                {
                    "event_id": "aevt_" + f"{index:032x}",
                    "now": "2026-11-19 00:00:00",
                    "resource": "ares_registered_result_" + "ab" * 32,
                },
            )
    engine.dispose()
    before = _audit_snapshot(database_url)
    assert len(before["audit"]) == 3
    for downgrade in (
        lambda: command.downgrade(_config(database_url), "20261119_32"),
        lambda: downgrade_database(database_url, "20261119_32"),
    ):
        with pytest.raises(Exception) as caught:
            downgrade()
        assert "cannot downgrade" in _exception_text(caught.value)
        assert _audit_snapshot(database_url) == before
        assert "_alembic_tmp_security_audit_events" not in before["tables"]


def _audit_snapshot(database_url: str) -> dict[str, object]:
    engine = create_engine(database_url)
    with engine.connect() as connection:
        version = connection.scalar(text("SELECT version_num FROM alembic_version"))
        audit = connection.execute(
            text(
                "SELECT event_id, action, resource_kind FROM security_audit_events "
                "ORDER BY event_id"
            )
        ).fetchall()
        tables = set(
            connection.execute(
                text("SELECT name FROM sqlite_master WHERE type = 'table'")
            ).scalars()
        )
    engine.dispose()
    return {"version": version, "audit": audit, "tables": tables}


def _exception_text(error: BaseException) -> str:
    parts: list[str] = []
    current: BaseException | None = error
    seen = 0
    while current is not None and seen < 6:
        parts.append(str(current))
        current = current.__cause__ or current.__context__
        seen += 1
    return "\n".join(parts)
