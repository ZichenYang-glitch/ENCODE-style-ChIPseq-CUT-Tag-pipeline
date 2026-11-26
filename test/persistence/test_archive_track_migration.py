"""The sealed-track schema upgrades a populated pre-archive database in place."""

from alembic import command
from alembic.config import Config
from pathlib import Path

from sqlalchemy import create_engine, inspect, text


def test_upgrade_from_registered_results_and_downgrade(tmp_path):
    url = f"sqlite:///{tmp_path / 'upgrade.db'}"
    cfg = Config()
    cfg.set_main_option(
        "script_location",
        str(
            Path(__file__).resolve().parents[2]
            / "src/encode_pipeline/persistence/alembic"
        ),
    )
    cfg.set_main_option("sqlalchemy.url", url)
    command.upgrade(cfg, "20261119_39")
    engine = create_engine(url)
    project = "prj_" + "2" * 32
    with engine.begin() as db:
        db.execute(
            text(
                "INSERT INTO projects(project_id, display_name, kind, created_at) VALUES (:id,'existing','user',CURRENT_TIMESTAMP)"
            ),
            {"id": project},
        )
    command.upgrade(cfg, "20261125_77")
    assert "archived_track_versions" in inspect(engine).get_table_names()
    with engine.connect() as db:
        assert (
            db.scalar(
                text("SELECT display_name FROM projects WHERE project_id = :id"),
                {"id": project},
            )
            == "existing"
        )
        assert db.scalar(text("SELECT count(*) FROM archived_track_versions")) == 0
    command.downgrade(cfg, "20261119_39")
    assert "archived_track_versions" not in inspect(engine).get_table_names()
    assert "registered_result_revisions" in inspect(engine).get_table_names()
    engine.dispose()
