"""Artifact-publication migration integrity and reversibility tests."""

from __future__ import annotations

from sqlalchemy import inspect, text
from sqlalchemy.exc import IntegrityError

from encode_pipeline.persistence import create_database_engine
from encode_pipeline.persistence.migrations import (
    downgrade_database,
    upgrade_database,
)


GENERATION = "artifactgen-" + "a" * 64
REVISION = "artifactrev-" + "b" * 64


def _insert_legacy_artifact(connection) -> None:
    connection.execute(
        text(
            "INSERT INTO runs "
            "(run_id, workflow_id, inputs, status, created_at, updated_at, tags) "
            "VALUES ('legacy-run', 'legacy-workflow', '{}', 'succeeded', "
            "'2026-08-07 00:00:00', '2026-08-07 00:00:00', '{}')"
        )
    )
    connection.execute(
        text(
            "INSERT INTO run_artifacts "
            "(artifact_id, run_id, artifact_type, name, uri, produced_at, "
            "revision, artifact_metadata) VALUES "
            "('artifact-a', 'legacy-run', 'file', 'artifact-a.txt', "
            "'run://legacy-run/artifact-a', '2026-08-07 00:00:00', "
            ":revision, '{}')"
        ),
        {"revision": REVISION},
    )


def test_rev11_to_rev12_adds_zero_backfill_append_only_publications(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'legacy.db'}"
    upgrade_database(database_url, "20260803_11")
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        _insert_legacy_artifact(connection)
    engine.dispose()

    upgrade_database(database_url, "20260807_12")
    engine = create_database_engine(database_url)
    inspector = inspect(engine)

    assert "artifact_publications" in inspector.get_table_names()
    assert inspector.get_pk_constraint("artifact_publications")[
        "constrained_columns"
    ] == ["run_id", "artifact_id", "artifact_generation"]
    assert inspector.get_foreign_keys("artifact_publications") == [
        {
            "name": "fk_artifact_publications_run",
            "constrained_columns": ["run_id"],
            "referred_schema": None,
            "referred_table": "runs",
            "referred_columns": ["run_id"],
            "options": {"ondelete": "RESTRICT"},
        }
    ]
    assert {
        item["name"] for item in inspector.get_indexes("artifact_publications")
    } == {
        "ix_artifact_publications_output_type_published",
        "ix_artifact_publications_published",
        "ix_artifact_publications_run_published",
    }

    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260807_12"
        )
        assert connection.scalar(text("SELECT count(*) FROM run_artifacts")) == 1
        assert (
            connection.scalar(text("SELECT count(*) FROM artifact_publications")) == 0
        )
        trigger_names = {
            row[0]
            for row in connection.execute(
                text(
                    "SELECT name FROM sqlite_master "
                    "WHERE type = 'trigger' AND tbl_name = 'artifact_publications'"
                )
            )
        }
        assert trigger_names == {
            "trg_artifact_publications_no_delete",
            "trg_artifact_publications_no_update",
        }
    engine.dispose()


def test_publication_fk_checks_and_append_only_triggers_fail_closed(tmp_path) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        connection.execute(
            text(
                "INSERT INTO runs "
                "(run_id, workflow_id, inputs, status, created_at, updated_at, tags) "
                "VALUES ('run-a', 'workflow-a', '{}', 'succeeded', "
                "'2026-08-07 00:00:00', '2026-08-07 00:00:00', '{}')"
            )
        )
        connection.execute(
            text(
                "INSERT INTO artifact_publications "
                "(run_id, artifact_id, artifact_generation, artifact_revision, "
                "output_type, published_at) VALUES "
                "('run-a', 'artifact-a', :generation, :revision, "
                "'summary', '2026-08-07 00:01:00')"
            ),
            {"generation": GENERATION, "revision": REVISION},
        )

    for statement in (
        "UPDATE artifact_publications SET output_type = 'other'",
        "DELETE FROM artifact_publications",
        "DELETE FROM runs WHERE run_id = 'run-a'",
    ):
        with engine.begin() as connection:
            try:
                connection.execute(text(statement))
            except IntegrityError:
                pass
            else:  # pragma: no cover - explicit database fail-closed guard
                raise AssertionError(f"statement unexpectedly succeeded: {statement}")

    invalid_rows = (
        ("artifact-a", "bad-generation", REVISION, "summary"),
        ("artifact-a", GENERATION, "bad-revision", "summary"),
        ("../artifact", GENERATION, REVISION, "summary"),
        ("artifact-a", GENERATION, REVISION, "bad/output"),
    )
    for artifact_id, generation, revision, output_type in invalid_rows:
        with engine.begin() as connection:
            try:
                connection.execute(
                    text(
                        "INSERT INTO artifact_publications "
                        "(run_id, artifact_id, artifact_generation, "
                        "artifact_revision, output_type, published_at) VALUES "
                        "('run-a', :artifact_id, :generation, :revision, "
                        ":output_type, '2026-08-07 00:02:00')"
                    ),
                    {
                        "artifact_id": artifact_id,
                        "generation": generation,
                        "revision": revision,
                        "output_type": output_type,
                    },
                )
            except IntegrityError:
                pass
            else:  # pragma: no cover - explicit database fail-closed guard
                raise AssertionError("invalid publication unexpectedly persisted")
    engine.dispose()


def test_rev12_downgrade_removes_only_publication_catalog(tmp_path) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        connection.execute(
            text(
                "INSERT INTO runs "
                "(run_id, workflow_id, inputs, status, created_at, updated_at, tags) "
                "VALUES ('run-a', 'workflow-a', '{}', 'succeeded', "
                "'2026-08-07 00:00:00', '2026-08-07 00:00:00', '{}')"
            )
        )
        connection.execute(
            text(
                "INSERT INTO artifact_publications "
                "(run_id, artifact_id, artifact_generation, artifact_revision, "
                "output_type, published_at) VALUES "
                "('run-a', 'artifact-a', :generation, :revision, "
                "'summary', '2026-08-07 00:01:00')"
            ),
            {"generation": GENERATION, "revision": REVISION},
        )
    engine.dispose()

    downgrade_database(database_url, "20260803_11")
    engine = create_database_engine(database_url)
    assert "artifact_publications" not in inspect(engine).get_table_names()
    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260803_11"
        )
        assert (
            connection.scalar(text("SELECT count(*) FROM runs WHERE run_id = 'run-a'"))
            == 1
        )
    engine.dispose()
