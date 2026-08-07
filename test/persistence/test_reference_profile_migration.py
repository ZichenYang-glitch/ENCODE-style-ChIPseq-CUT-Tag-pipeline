"""Reference Profile migration and legacy-null evidence tests."""

from __future__ import annotations

from sqlalchemy import inspect, text

from encode_pipeline.persistence import create_database_engine
from encode_pipeline.persistence.migrations import upgrade_database


REFERENCE_TABLES = {
    "reference_profiles",
    "reference_profile_revisions",
    "reference_profile_workflow_bindings",
    "snapshot_reference_bindings",
    "run_reference_bindings",
}


def test_rev10_to_rev11_adds_reference_catalog_without_legacy_backfill(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'legacy.db'}"
    upgrade_database(database_url, "20260726_10")
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        connection.execute(
            text(
                "INSERT INTO runs "
                "(run_id, workflow_id, inputs, status, created_at, updated_at, tags) "
                "VALUES ('legacy-run', 'legacy-workflow', '{}', 'succeeded', "
                "'2026-08-03 00:00:00', '2026-08-03 00:00:00', '{}')"
            )
        )
        connection.execute(
            text(
                "INSERT INTO validated_input_snapshots "
                "(snapshot_id, workflow_id, adapter_version, schema_version, "
                "schema_dialect, canonical_payload, payload_digest_scheme, "
                "payload_digest, validation_outcome, validation_issue_codes, "
                "validated_at, expires_at, build_adapter_version, build_scheme, "
                "build_logical_entrypoint, build_digest, build_captured_at, "
                "consumed_run_id, consumed_at) VALUES "
                "('vsnap_11111111111111111111111111111111', "
                "'legacy-workflow', '1.0.0', '1.0.0', "
                "'https://json-schema.org/draft/2020-12/schema', "
                ":payload, "
                "'sha256-framed-workflow-inputs-v1', :digest, "
                "'adapter_validation_succeeded', '[]', "
                "'2026-08-03 00:00:00', '2026-08-04 00:00:00', "
                "'1.0.0', 'sha256-tree-v1', 'legacy/workflow', :digest, "
                "'2026-08-03 00:00:00', NULL, NULL)"
            ),
            {
                "digest": "0" * 64,
                "payload": '{"config":{},"options":{},"samples":null}',
            },
        )
    engine.dispose()

    upgrade_database(database_url, "20260803_11")
    engine = create_database_engine(database_url)
    inspector = inspect(engine)

    assert REFERENCE_TABLES <= set(inspector.get_table_names())
    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260803_11"
        )
        assert connection.scalar(text("SELECT count(*) FROM reference_profiles")) == 0
        assert (
            connection.scalar(text("SELECT count(*) FROM run_reference_bindings")) == 0
        )
        assert (
            connection.scalar(text("SELECT count(*) FROM snapshot_reference_bindings"))
            == 0
        )
        assert (
            connection.scalar(
                text("SELECT count(*) FROM runs WHERE run_id = 'legacy-run'")
            )
            == 1
        )
        assert (
            connection.scalar(
                text(
                    "SELECT count(*) FROM validated_input_snapshots "
                    "WHERE snapshot_id = 'vsnap_11111111111111111111111111111111'"
                )
            )
            == 1
        )
    engine.dispose()
