"""Execute the registered-result migration from its source file.

Admission runs a verified copy, so shard coverage would not attribute those
statements to this revision. Pointing Alembic at the source tree does.
"""

from __future__ import annotations

from pathlib import Path

from alembic.config import Config
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
    assert 'down_revision = "20261119_18"' in source
    database_url = f"sqlite:///{tmp_path / 'revision-39.db'}"
    command.upgrade(_config(database_url), "20261119_39")
    names = _names(database_url)
    assert "registered_results" in names
    assert "snapshot_experiment_bindings" in names
    command.downgrade(_config(database_url), "20261119_18")
    names = _names(database_url)
    assert "registered_results" not in names
    assert "registered_result_revisions" not in names
    assert "snapshot_experiment_bindings" in names
