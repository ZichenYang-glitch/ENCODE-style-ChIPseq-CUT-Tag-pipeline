"""Compatibility contract tests for Bulk RNA-seq execution persistence."""

from __future__ import annotations

import json
from pathlib import Path
import re

from alembic.config import Config
from alembic.script import ScriptDirectory
import pytest
from sqlalchemy import inspect, text

from encode_pipeline.adapters.bulk_rnaseq.persistence_projection import (
    PersistenceProjectionError,
    SCHEMA_PROJECTION_SCHEME,
    parse_schema_projection_spec,
    schema_capability_projection,
    schema_projection_sha256,
)
from encode_pipeline.persistence import create_database_engine
from encode_pipeline.persistence.migrations import upgrade_database
from encode_pipeline.persistence.models import Base


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = (
    PROJECT_ROOT
    / "src/encode_pipeline/contracts/nfcore_rnaseq"
    / "execution-persistence-contract-1.1.0.json"
)
EXPECTED_TOP_LEVEL_KEYS = {
    "schema_version",
    "contract_id",
    "contract_version",
    "minimum_supported_revision",
    "capabilities",
    "required_revisions",
    "schema_projection",
}
EXPECTED_CAPABILITIES = [
    "sqlite.run-aggregate/v1",
    "sqlite.validated-snapshot-atomic-consume/v1",
    "sqlite.durable-assignment-cancellation/v1",
    "sqlite.artifact-qc-generation/v1",
    "sqlite.project-sample-binding/v1",
    "sqlite.compatibility-input-binding/v1",
]
EXPECTED_REQUIRED_REVISIONS = [
    "20260711_01",
    "20260711_02",
    "20260712_03",
    "20260712_04",
    "20260712_05",
    "20260714_06",
    "20260717_08",
    "20260726_09",
    "20260726_10",
]
EXPECTED_REQUIRED_SCHEMA = {
    "projects": [
        "kind",
        "project_id",
    ],
    "run_artifacts": [
        "artifact_id",
        "artifact_metadata",
        "artifact_type",
        "mime_type",
        "name",
        "produced_at",
        "revision",
        "run_id",
        "uri",
    ],
    "run_events": [
        "context",
        "event_id",
        "event_type",
        "issue",
        "message",
        "run_id",
        "sequence",
        "stage",
        "status",
        "timestamp",
    ],
    "run_execution_assignments": [
        "backend",
        "cancellation_acknowledged_at",
        "cancellation_reason",
        "cancellation_requested_at",
        "claimed_at",
        "created_at",
        "dispatched_at",
        "job_id",
        "queue_name",
        "run_id",
    ],
    "run_input_bindings": [
        "adapter_contract_version",
        "binding_digest",
        "binding_digest_scheme",
        "binding_mode",
        "project_id",
        "project_sample_binding_digest",
        "run_id",
        "workflow_id",
        "workflow_inputs_digest",
    ],
    "run_logs": [
        "chunk_id",
        "lines",
        "run_id",
        "sequence",
        "stream_name",
        "timestamp",
    ],
    "run_project_bindings": [
        "binding_digest",
        "binding_digest_scheme",
        "binding_mode",
        "project_id",
        "provenance",
        "run_id",
        "workflow_inputs_digest",
    ],
    "run_qc_metrics": [
        "assay",
        "display_name",
        "experiment_id",
        "metric_id",
        "metric_key",
        "produced_at",
        "qc_flag",
        "run_id",
        "sample_id",
        "scope",
        "source_artifact_id",
        "unit",
        "value_text",
    ],
    "run_result_attempts": [
        "artifact_generation",
        "attempt_id",
        "result_kind",
        "run_id",
    ],
    "run_result_states": [
        "artifact_attempt_id",
        "artifact_attempt_status",
        "artifact_generation",
        "artifact_manifest_digest",
        "artifact_outcome",
        "artifact_reason_code",
        "artifact_revision",
        "qc_artifact_generation",
        "qc_attempt_artifact_generation",
        "qc_attempt_id",
        "qc_attempt_status",
        "qc_generation",
        "qc_manifest_digest",
        "qc_outcome",
        "qc_reason_code",
        "qc_revision",
        "run_id",
    ],
    "run_samples": [
        "ordinal",
        "payload_digest",
        "project_id",
        "run_id",
        "sample_revision_id",
    ],
    "run_workflow_build_identities": [
        "adapter_version",
        "captured_at",
        "digest",
        "logical_entrypoint",
        "run_id",
        "scheme",
        "workflow_id",
    ],
    "runs": [
        "cancellation_reason",
        "created_at",
        "current_stage",
        "ended_at",
        "error",
        "inputs",
        "run_id",
        "started_at",
        "status",
        "tags",
        "updated_at",
        "workflow_id",
    ],
    "sample_revisions": [
        "payload_digest",
        "project_id",
        "sample_id",
        "sample_revision_id",
    ],
    "samples": [
        "project_id",
        "sample_id",
    ],
    "snapshot_input_bindings": [
        "adapter_contract_version",
        "binding_digest",
        "binding_digest_scheme",
        "binding_mode",
        "project_id",
        "project_sample_binding_digest",
        "snapshot_id",
        "workflow_id",
        "workflow_inputs_digest",
    ],
    "snapshot_project_bindings": [
        "binding_digest",
        "binding_digest_scheme",
        "binding_mode",
        "project_id",
        "provenance",
        "snapshot_id",
        "workflow_inputs_digest",
    ],
    "snapshot_sample_revisions": [
        "ordinal",
        "payload_digest",
        "project_id",
        "sample_revision_id",
        "snapshot_id",
    ],
    "validated_input_snapshots": [
        "adapter_version",
        "build_adapter_version",
        "build_captured_at",
        "build_digest",
        "build_logical_entrypoint",
        "build_scheme",
        "canonical_payload",
        "consumed_at",
        "consumed_run_id",
        "expires_at",
        "payload_digest",
        "payload_digest_scheme",
        "schema_dialect",
        "schema_version",
        "snapshot_id",
        "validated_at",
        "validation_issue_codes",
        "validation_outcome",
        "workflow_id",
    ],
}
UNACTIVATED_MANAGED_INPUT_TABLES = {
    "input_file_revisions",
    "input_files",
    "project_storage_pool_bindings",
    "run_input_members",
    "run_input_uses",
    "snapshot_input_members",
    "snapshot_input_uses",
    "storage_pools",
}


def _load_contract() -> tuple[bytes, dict[str, object]]:
    raw = CONTRACT_PATH.read_bytes()
    return raw, json.loads(raw)


def _projection_spec(contract: dict[str, object]):
    return parse_schema_projection_spec(contract["schema_projection"])


def _required_schema(contract: dict[str, object]) -> dict[str, list[str]]:
    projection = contract["schema_projection"]
    assert isinstance(projection, dict)
    tables = projection["tables"]
    assert isinstance(tables, dict)
    required_schema: dict[str, list[str]] = {}
    for table_name, table in tables.items():
        assert isinstance(table_name, str)
        assert isinstance(table, dict)
        columns = table["columns"]
        assert isinstance(columns, list)
        assert all(isinstance(column, str) for column in columns)
        required_schema[table_name] = columns
    return required_schema


def _probe_spec(
    *,
    columns: list[str],
    primary_key: bool = False,
    foreign_keys: list[dict[str, object]] | None = None,
    unique_constraints: list[dict[str, object]] | None = None,
    indexes: list[str] | None = None,
    check_constraints: list[str] | None = None,
    triggers: list[str] | None = None,
    table_options: list[str] | None = None,
):
    tables: dict[str, object] = {}
    for foreign_key in foreign_keys or []:
        referred_table = foreign_key["referred_table"]
        referred_columns = foreign_key["referred_columns"]
        assert isinstance(referred_table, str)
        assert isinstance(referred_columns, list)
        tables[referred_table] = {
            "check_constraints": [],
            "columns": referred_columns,
            "foreign_keys": [],
            "indexes": [],
            "primary_key": False,
            "table_options": [],
            "triggers": [],
            "unique_constraints": [],
        }
    tables["probe"] = {
        "check_constraints": check_constraints or [],
        "columns": columns,
        "foreign_keys": foreign_keys or [],
        "indexes": indexes or [],
        "primary_key": primary_key,
        "table_options": table_options or [],
        "triggers": triggers or [],
        "unique_constraints": unique_constraints or [],
    }
    return parse_schema_projection_spec(
        {
            "scheme": SCHEMA_PROJECTION_SCHEME,
            "sha256": "0" * 64,
            "tables": dict(sorted(tables.items())),
        }
    )


def _canonical_json_bytes(value: object) -> bytes:
    return (
        json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        )
        + "\n"
    ).encode("utf-8")


def _string_leaves(value: object):
    if isinstance(value, str):
        yield value
    elif isinstance(value, list):
        for item in value:
            yield from _string_leaves(item)
    elif isinstance(value, dict):
        for key, item in value.items():
            yield key
            yield from _string_leaves(item)


def _missing_columns(
    required: dict[str, list[str]],
    available: dict[str, set[str]],
) -> dict[str, list[str]]:
    return {
        table_name: sorted(set(column_names) - available.get(table_name, set()))
        for table_name, column_names in required.items()
        if set(column_names) - available.get(table_name, set())
    }


def test_bulk_rnaseq_persistence_contract_is_exact_canonical_and_path_free():
    raw, contract = _load_contract()
    spec = _projection_spec(contract)

    assert set(contract) == EXPECTED_TOP_LEVEL_KEYS
    assert contract["schema_version"] == "1.1.0"
    assert contract["contract_id"] == "bulk-rnaseq-execution-persistence"
    assert contract["contract_version"] == "1.1.0"
    assert contract["minimum_supported_revision"] == "20260714_07"
    assert contract["capabilities"] == EXPECTED_CAPABILITIES
    assert contract["required_revisions"] == EXPECTED_REQUIRED_REVISIONS
    assert spec.scheme == SCHEMA_PROJECTION_SCHEME
    assert _required_schema(contract) == EXPECTED_REQUIRED_SCHEMA
    assert raw == _canonical_json_bytes(contract)
    assert "20260714_07" not in contract["required_revisions"]
    assert UNACTIVATED_MANAGED_INPUT_TABLES.isdisjoint(_required_schema(contract))
    for value in _string_leaves(contract):
        assert not value.startswith(("/", "\\", "./", "../", "~"))
        assert not re.match(r"^[A-Za-z]:[\\/]", value)
        assert "://" not in value


def test_declared_revisions_exist_below_the_sole_alembic_head():
    _, contract = _load_contract()
    config = Config()
    config.set_main_option(
        "script_location",
        str(
            (PROJECT_ROOT / "src/encode_pipeline/persistence/alembic").resolve()
        ).replace("%", "%%"),
    )
    scripts = ScriptDirectory.from_config(config)
    installed_revisions = {revision.revision for revision in scripts.walk_revisions()}
    heads = scripts.get_heads()

    assert len(heads) == 1
    sole_head_lineage = {
        revision.revision for revision in scripts.iterate_revisions(heads[0], "base")
    }
    required_revisions = set(contract["required_revisions"])
    minimum_supported_revision = contract["minimum_supported_revision"]

    assert required_revisions.issubset(installed_revisions)
    assert required_revisions.issubset(sole_head_lineage)
    assert minimum_supported_revision in installed_revisions
    assert minimum_supported_revision in sole_head_lineage


def test_required_schema_is_present_in_models_and_fresh_database(tmp_path: Path):
    _, contract = _load_contract()
    required_schema = _required_schema(contract)
    metadata_schema = {
        table_name: set(table.columns.keys())
        for table_name, table in Base.metadata.tables.items()
    }

    assert _missing_columns(required_schema, metadata_schema) == {}

    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    try:
        inspector = inspect(engine)
        database_schema = {
            table_name: {column["name"] for column in inspector.get_columns(table_name)}
            for table_name in inspector.get_table_names()
        }
    finally:
        engine.dispose()

    assert _missing_columns(required_schema, database_schema) == {}


def test_schema_capability_projection_binds_bulk_semantics_not_unrelated_tables(
    tmp_path: Path,
):
    _, contract = _load_contract()
    required_schema = _required_schema(contract)
    spec = _projection_spec(contract)
    migrated_url = f"sqlite:///{tmp_path / 'migrated.db'}"
    upgrade_database(migrated_url)
    migrated_engine = create_database_engine(migrated_url)
    modeled_engine = create_database_engine(f"sqlite:///{tmp_path / 'modeled.db'}")
    Base.metadata.create_all(modeled_engine)
    try:
        migrated_inspector = inspect(migrated_engine)
        modeled_inspector = inspect(modeled_engine)
        for table_name, required_columns in required_schema.items():
            migrated_columns = {
                column["name"]: column
                for column in migrated_inspector.get_columns(table_name)
            }
            modeled_columns = {
                column["name"]: column
                for column in modeled_inspector.get_columns(table_name)
            }
            for column_name in required_columns:
                assert {
                    key: migrated_columns[column_name][key]
                    for key in ("nullable", "primary_key")
                } == {
                    key: modeled_columns[column_name][key]
                    for key in ("nullable", "primary_key")
                }
                assert str(migrated_columns[column_name]["type"]) == str(
                    modeled_columns[column_name]["type"]
                )
        projection = schema_capability_projection(migrated_engine, spec)
        serialized_projection = json.dumps(projection, sort_keys=True)
        expected_digest = spec.sha256
        assert schema_projection_sha256(migrated_engine, spec) == expected_digest
        assert "CREATE TABLE" not in serialized_projection.upper()
        assert "create_table_sql" not in serialized_projection

        with migrated_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE TABLE reference_catalog_review_only "
                    "(reference_id TEXT PRIMARY KEY)"
                )
            )
            connection.execute(
                text(
                    "CREATE TRIGGER bulk_run_status_guard "
                    "BEFORE UPDATE OF status ON runs "
                    "BEGIN SELECT RAISE(ABORT, 'review-only'); END"
                )
            )
        assert schema_projection_sha256(migrated_engine, spec) == expected_digest
    finally:
        migrated_engine.dispose()
        modeled_engine.dispose()


@pytest.mark.parametrize(
    ("first_column", "second_column"),
    [
        ("status TEXT NOT NULL", "status INTEGER NOT NULL"),
        ("status TEXT NOT NULL", "status TEXT"),
        ("status TEXT DEFAULT 'pending'", "status TEXT DEFAULT 'ready'"),
        ("status TEXT COLLATE BINARY", "status TEXT COLLATE NOCASE"),
        ("status TEXT", "status TEXT PRIMARY KEY"),
    ],
)
def test_selected_column_semantics_change_projection(
    tmp_path: Path,
    first_column: str,
    second_column: str,
):
    first_engine = create_database_engine(f"sqlite:///{tmp_path / 'first.db'}")
    second_engine = create_database_engine(f"sqlite:///{tmp_path / 'second.db'}")
    try:
        with first_engine.begin() as connection:
            connection.execute(text(f"CREATE TABLE probe ({first_column})"))
        with second_engine.begin() as connection:
            connection.execute(text(f"CREATE TABLE probe ({second_column})"))
        spec = _probe_spec(columns=["status"])

        assert schema_projection_sha256(
            first_engine,
            spec,
        ) != schema_projection_sha256(
            second_engine,
            spec,
        )
    finally:
        first_engine.dispose()
        second_engine.dispose()


def test_selected_not_null_conflict_policy_changes_projection(tmp_path: Path):
    plain_engine = create_database_engine(f"sqlite:///{tmp_path / 'plain.db'}")
    abort_engine = create_database_engine(f"sqlite:///{tmp_path / 'abort.db'}")
    ignore_engine = create_database_engine(f"sqlite:///{tmp_path / 'ignore.db'}")
    spec = _probe_spec(columns=["status"])
    try:
        with plain_engine.begin() as connection:
            connection.execute(text("CREATE TABLE probe (status TEXT NOT NULL)"))
        with abort_engine.begin() as connection:
            connection.execute(
                text("CREATE TABLE probe (status TEXT NOT NULL ON CONFLICT ABORT)")
            )
        with ignore_engine.begin() as connection:
            connection.execute(
                text("CREATE TABLE probe (status TEXT NOT NULL ON CONFLICT IGNORE)")
            )

        assert schema_projection_sha256(
            plain_engine,
            spec,
        ) != schema_projection_sha256(ignore_engine, spec)
        assert schema_projection_sha256(
            plain_engine,
            spec,
        ) == schema_projection_sha256(abort_engine, spec)
    finally:
        plain_engine.dispose()
        abort_engine.dispose()
        ignore_engine.dispose()


def test_unselected_projection_probe_name_collision_does_not_change_projection(
    tmp_path: Path,
):
    plain_engine = create_database_engine(f"sqlite:///{tmp_path / 'plain.db'}")
    expanded_engine = create_database_engine(f"sqlite:///{tmp_path / 'expanded.db'}")
    spec = _probe_spec(columns=["status"])
    try:
        for engine in (plain_engine, expanded_engine):
            with engine.begin() as connection:
                connection.execute(text("CREATE TABLE probe (status TEXT)"))
        with expanded_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE INDEX helixweave_projection_0_probe_status ON probe(status)"
                )
            )

        assert schema_projection_sha256(
            plain_engine,
            spec,
        ) == schema_projection_sha256(expanded_engine, spec)
    finally:
        plain_engine.dispose()
        expanded_engine.dispose()


def test_projection_probe_avoids_case_insensitive_index_name_collision(
    tmp_path: Path,
):
    plain_engine = create_database_engine(f"sqlite:///{tmp_path / 'plain.db'}")
    expanded_engine = create_database_engine(f"sqlite:///{tmp_path / 'expanded.db'}")
    spec = _probe_spec(columns=["status"])
    try:
        for engine in (plain_engine, expanded_engine):
            with engine.begin() as connection:
                connection.execute(text("CREATE TABLE probe (status TEXT)"))
        with expanded_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE INDEX HELIXWEAVE_PROJECTION_0_PROBE_STATUS_0 "
                    "ON probe(status)"
                )
            )

        assert schema_projection_sha256(
            plain_engine,
            spec,
        ) == schema_projection_sha256(expanded_engine, spec)
    finally:
        plain_engine.dispose()
        expanded_engine.dispose()


def test_unselected_cross_namespace_duplicate_name_does_not_change_projection(
    tmp_path: Path,
):
    plain_engine = create_database_engine(f"sqlite:///{tmp_path / 'plain.db'}")
    expanded_engine = create_database_engine(f"sqlite:///{tmp_path / 'expanded.db'}")
    spec = _probe_spec(columns=["status"])
    try:
        for engine in (plain_engine, expanded_engine):
            with engine.begin() as connection:
                connection.execute(text("CREATE TABLE probe (status TEXT)"))
        with expanded_engine.begin() as connection:
            connection.execute(
                text("CREATE TRIGGER probe AFTER INSERT ON probe BEGIN SELECT 1; END")
            )

        assert schema_projection_sha256(
            plain_engine,
            spec,
        ) == schema_projection_sha256(expanded_engine, spec)
    finally:
        plain_engine.dispose()
        expanded_engine.dispose()


def test_unrelated_nullable_column_on_selected_table_does_not_change_projection(
    tmp_path: Path,
):
    _, contract = _load_contract()
    spec = _projection_spec(contract)
    engine = create_database_engine(f"sqlite:///{tmp_path / 'platform.db'}")
    upgrade_database(str(engine.url))
    try:
        before = schema_projection_sha256(engine, spec)
        with engine.begin() as connection:
            connection.execute(
                text("ALTER TABLE runs ADD COLUMN review_only_note TEXT")
            )

        assert schema_projection_sha256(engine, spec) == before
    finally:
        engine.dispose()


def test_unselected_components_on_selected_table_do_not_change_projection(
    tmp_path: Path,
):
    plain_engine = create_database_engine(f"sqlite:///{tmp_path / 'plain.db'}")
    expanded_engine = create_database_engine(f"sqlite:///{tmp_path / 'expanded.db'}")
    spec = _probe_spec(columns=["status"])
    try:
        with plain_engine.begin() as connection:
            connection.execute(text("CREATE TABLE probe (status TEXT)"))
        with expanded_engine.begin() as connection:
            connection.execute(
                text("CREATE TABLE parent (parent_id TEXT PRIMARY KEY);")
            )
            connection.execute(
                text(
                    "CREATE TABLE probe ("
                    "status TEXT, review_note TEXT, parent_id TEXT, "
                    "UNIQUE(review_note), CHECK(length(review_note) < 100), "
                    "FOREIGN KEY(parent_id) REFERENCES parent(parent_id)"
                    ")"
                )
            )
            connection.execute(
                text("CREATE INDEX ix_probe_review ON probe(review_note)")
            )
            connection.execute(
                text(
                    "CREATE TRIGGER trg_probe_review BEFORE UPDATE OF review_note "
                    "ON probe BEGIN SELECT RAISE(ABORT, 'review-only'); END"
                )
            )

        assert schema_projection_sha256(
            plain_engine,
            spec,
        ) == schema_projection_sha256(expanded_engine, spec)
    finally:
        plain_engine.dispose()
        expanded_engine.dispose()


def test_unselected_constraint_with_noncontract_name_does_not_change_projection(
    tmp_path: Path,
):
    plain_engine = create_database_engine(f"sqlite:///{tmp_path / 'plain.db'}")
    expanded_engine = create_database_engine(f"sqlite:///{tmp_path / 'expanded.db'}")
    spec = _probe_spec(columns=["status"])
    try:
        with plain_engine.begin() as connection:
            connection.execute(text("CREATE TABLE probe (status TEXT)"))
        with expanded_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE TABLE probe ("
                    "status TEXT, review_note TEXT, "
                    'CONSTRAINT "review-only" CHECK(review_note IS NULL)'
                    ")"
                )
            )

        assert schema_projection_sha256(
            plain_engine,
            spec,
        ) == schema_projection_sha256(expanded_engine, spec)
    finally:
        plain_engine.dispose()
        expanded_engine.dispose()


def test_unselected_shorthand_foreign_key_does_not_change_projection(
    tmp_path: Path,
):
    plain_engine = create_database_engine(f"sqlite:///{tmp_path / 'plain.db'}")
    expanded_engine = create_database_engine(f"sqlite:///{tmp_path / 'expanded.db'}")
    spec = _probe_spec(columns=["status"])
    try:
        with plain_engine.begin() as connection:
            connection.execute(text("CREATE TABLE probe (status TEXT)"))
        with expanded_engine.begin() as connection:
            connection.execute(text("CREATE TABLE parent (parent_id TEXT PRIMARY KEY)"))
            connection.execute(
                text(
                    "CREATE TABLE probe ("
                    "status TEXT, parent_id TEXT, "
                    "FOREIGN KEY(parent_id) REFERENCES parent"
                    ")"
                )
            )

        assert schema_projection_sha256(
            plain_engine,
            spec,
        ) == schema_projection_sha256(expanded_engine, spec)
    finally:
        plain_engine.dispose()
        expanded_engine.dispose()


def test_unselected_shorthand_foreign_key_is_ignored_beside_selected_fk(
    tmp_path: Path,
):
    plain_engine = create_database_engine(f"sqlite:///{tmp_path / 'plain.db'}")
    expanded_engine = create_database_engine(f"sqlite:///{tmp_path / 'expanded.db'}")
    selector = [
        {
            "columns": ["parent_id"],
            "referred_columns": ["parent_id"],
            "referred_table": "parent",
        }
    ]
    spec = _probe_spec(
        columns=["parent_id", "status"],
        foreign_keys=selector,
    )
    try:
        for engine in (plain_engine, expanded_engine):
            with engine.begin() as connection:
                connection.execute(
                    text("CREATE TABLE parent (parent_id TEXT PRIMARY KEY)")
                )
                connection.execute(
                    text("CREATE TABLE other (other_id TEXT PRIMARY KEY)")
                )
        with plain_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE TABLE probe ("
                    "status TEXT, parent_id TEXT, other_id TEXT, "
                    "FOREIGN KEY(parent_id) REFERENCES parent(parent_id)"
                    ")"
                )
            )
        with expanded_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE TABLE probe ("
                    "status TEXT, parent_id TEXT, other_id TEXT, "
                    "FOREIGN KEY(parent_id) REFERENCES parent(parent_id), "
                    "FOREIGN KEY(other_id) REFERENCES other"
                    ")"
                )
            )

        assert schema_projection_sha256(
            plain_engine,
            spec,
        ) == schema_projection_sha256(expanded_engine, spec)
    finally:
        plain_engine.dispose()
        expanded_engine.dispose()


def test_shorthand_fk_to_parent_pk_is_ignored_when_other_target_is_selected(
    tmp_path: Path,
):
    plain_engine = create_database_engine(f"sqlite:///{tmp_path / 'plain.db'}")
    expanded_engine = create_database_engine(f"sqlite:///{tmp_path / 'expanded.db'}")
    selector = [
        {
            "columns": ["parent_key"],
            "referred_columns": ["alternate_key"],
            "referred_table": "parent",
        }
    ]
    spec = _probe_spec(
        columns=["parent_key"],
        foreign_keys=selector,
    )
    try:
        for engine in (plain_engine, expanded_engine):
            with engine.begin() as connection:
                connection.execute(
                    text(
                        "CREATE TABLE parent ("
                        "parent_id TEXT PRIMARY KEY, "
                        "alternate_key TEXT UNIQUE"
                        ")"
                    )
                )
        with plain_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE TABLE probe ("
                    "parent_key TEXT, "
                    "FOREIGN KEY(parent_key) "
                    "REFERENCES parent(alternate_key)"
                    ")"
                )
            )
        with expanded_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE TABLE probe ("
                    "parent_key TEXT, "
                    "FOREIGN KEY(parent_key) "
                    "REFERENCES parent(alternate_key), "
                    "FOREIGN KEY(parent_key) REFERENCES parent"
                    ")"
                )
            )

        assert schema_projection_sha256(
            plain_engine,
            spec,
        ) == schema_projection_sha256(expanded_engine, spec)
    finally:
        plain_engine.dispose()
        expanded_engine.dispose()


def test_selected_shorthand_foreign_key_resolves_parent_primary_key(
    tmp_path: Path,
):
    explicit_engine = create_database_engine(f"sqlite:///{tmp_path / 'explicit.db'}")
    shorthand_engine = create_database_engine(f"sqlite:///{tmp_path / 'shorthand.db'}")
    selector = [
        {
            "columns": ["parent_id"],
            "referred_columns": ["parent_id"],
            "referred_table": "parent",
        }
    ]
    spec = _probe_spec(
        columns=["parent_id"],
        foreign_keys=selector,
    )
    try:
        for engine in (explicit_engine, shorthand_engine):
            with engine.begin() as connection:
                connection.execute(
                    text("CREATE TABLE parent (parent_id TEXT PRIMARY KEY)")
                )
        with explicit_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE TABLE probe ("
                    "parent_id TEXT, "
                    "FOREIGN KEY(parent_id) REFERENCES parent(parent_id)"
                    ")"
                )
            )
        with shorthand_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE TABLE probe ("
                    "parent_id TEXT, "
                    "FOREIGN KEY(parent_id) REFERENCES parent"
                    ")"
                )
            )

        assert schema_projection_sha256(
            explicit_engine,
            spec,
        ) == schema_projection_sha256(shorthand_engine, spec)
    finally:
        explicit_engine.dispose()
        shorthand_engine.dispose()


def test_selected_foreign_key_options_change_and_target_drift_rejects(tmp_path: Path):
    plain_engine = create_database_engine(f"sqlite:///{tmp_path / 'plain.db'}")
    cascade_engine = create_database_engine(f"sqlite:///{tmp_path / 'cascade.db'}")
    deferred_engine = create_database_engine(f"sqlite:///{tmp_path / 'deferred.db'}")
    wrong_engine = create_database_engine(f"sqlite:///{tmp_path / 'wrong.db'}")
    selector = [
        {
            "columns": ["parent_id"],
            "referred_columns": ["parent_id"],
            "referred_table": "parent",
        }
    ]
    spec = _probe_spec(columns=["parent_id"], foreign_keys=selector)
    try:
        for engine, foreign_key in (
            (
                plain_engine,
                "FOREIGN KEY(parent_id) REFERENCES parent(parent_id)",
            ),
            (
                cascade_engine,
                "FOREIGN KEY(parent_id) REFERENCES parent(parent_id) ON DELETE CASCADE",
            ),
            (
                wrong_engine,
                "FOREIGN KEY(parent_id) REFERENCES parent(other_id)",
            ),
            (
                deferred_engine,
                "FOREIGN KEY(parent_id) REFERENCES parent(parent_id) "
                "DEFERRABLE INITIALLY DEFERRED",
            ),
        ):
            with engine.begin() as connection:
                connection.execute(
                    text(
                        "CREATE TABLE parent "
                        "(parent_id TEXT PRIMARY KEY, other_id TEXT UNIQUE)"
                    )
                )
                connection.execute(
                    text(f"CREATE TABLE probe (parent_id TEXT NOT NULL, {foreign_key})")
                )

        assert schema_projection_sha256(
            plain_engine,
            spec,
        ) != schema_projection_sha256(cascade_engine, spec)
        assert schema_projection_sha256(
            plain_engine,
            spec,
        ) != schema_projection_sha256(deferred_engine, spec)
        with pytest.raises(PersistenceProjectionError):
            schema_projection_sha256(wrong_engine, spec)
    finally:
        plain_engine.dispose()
        cascade_engine.dispose()
        deferred_engine.dispose()
        wrong_engine.dispose()


def test_selected_primary_and_unique_conflict_semantics_are_bound(tmp_path: Path):
    plain_engine = create_database_engine(f"sqlite:///{tmp_path / 'plain.db'}")
    changed_engine = create_database_engine(f"sqlite:///{tmp_path / 'changed.db'}")
    primary_spec = _probe_spec(columns=["status"], primary_key=True)
    unique_spec = _probe_spec(
        columns=["status"],
        unique_constraints=[{"columns": ["status"], "name": "uq_probe_status"}],
    )
    try:
        with plain_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE TABLE probe ("
                    "status TEXT NOT NULL, "
                    "CONSTRAINT uq_probe_status UNIQUE(status)"
                    ")"
                )
            )
        with changed_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE TABLE probe ("
                    "status TEXT NOT NULL, "
                    "CONSTRAINT uq_probe_status UNIQUE(status) ON CONFLICT IGNORE"
                    ")"
                )
            )
        assert schema_projection_sha256(
            plain_engine,
            unique_spec,
        ) != schema_projection_sha256(changed_engine, unique_spec)
    finally:
        plain_engine.dispose()
        changed_engine.dispose()

    plain_engine = create_database_engine(f"sqlite:///{tmp_path / 'plain-pk.db'}")
    changed_engine = create_database_engine(f"sqlite:///{tmp_path / 'changed-pk.db'}")
    try:
        with plain_engine.begin() as connection:
            connection.execute(
                text("CREATE TABLE probe (status TEXT, PRIMARY KEY(status))")
            )
        with changed_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE TABLE probe ("
                    "status TEXT, PRIMARY KEY(status) ON CONFLICT IGNORE"
                    ")"
                )
            )
        assert schema_projection_sha256(
            plain_engine,
            primary_spec,
        ) != schema_projection_sha256(changed_engine, primary_spec)
    finally:
        plain_engine.dispose()
        changed_engine.dispose()


@pytest.mark.parametrize(
    "extra_unique",
    (
        "CONSTRAINT uq_review UNIQUE(status)",
        "CONSTRAINT uq_review UNIQUE(status COLLATE NOCASE)",
    ),
)
def test_unselected_unique_on_selected_columns_does_not_change_projection(
    tmp_path: Path,
    extra_unique: str,
):
    plain_engine = create_database_engine(f"sqlite:///{tmp_path / 'plain.db'}")
    expanded_engine = create_database_engine(f"sqlite:///{tmp_path / 'expanded.db'}")
    spec = _probe_spec(
        columns=["status"],
        unique_constraints=[{"columns": ["status"], "name": "uq_probe_status"}],
    )
    try:
        with plain_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE TABLE probe ("
                    "status TEXT, "
                    "CONSTRAINT uq_probe_status UNIQUE(status)"
                    ")"
                )
            )
        with expanded_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE TABLE probe ("
                    "status TEXT, "
                    "CONSTRAINT uq_probe_status UNIQUE(status), "
                    f"{extra_unique}"
                    ")"
                )
            )

        assert schema_projection_sha256(
            plain_engine,
            spec,
        ) == schema_projection_sha256(expanded_engine, spec)
    finally:
        plain_engine.dispose()
        expanded_engine.dispose()


def test_selected_unique_with_uninterpretable_term_fails_closed(tmp_path: Path):
    engine = create_database_engine(f"sqlite:///{tmp_path / 'probe.db'}")
    spec = _probe_spec(
        columns=["status"],
        unique_constraints=[{"columns": ["status"], "name": "uq_probe_status"}],
    )
    try:
        with engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE TABLE probe ("
                    "status TEXT, "
                    "CONSTRAINT uq_probe_status UNIQUE(status COLLATE NOCASE)"
                    ")"
                )
            )
        with pytest.raises(PersistenceProjectionError):
            schema_projection_sha256(engine, spec)
    finally:
        engine.dispose()


def test_selected_unique_index_check_trigger_and_table_option_are_bound(
    tmp_path: Path,
):
    first_engine = create_database_engine(f"sqlite:///{tmp_path / 'first.db'}")
    second_engine = create_database_engine(f"sqlite:///{tmp_path / 'second.db'}")
    spec = _probe_spec(
        columns=["status"],
        unique_constraints=[{"columns": ["status"], "name": "uq_probe_status"}],
        indexes=["ix_probe_status"],
        check_constraints=["ck_probe_status"],
        triggers=["trg_probe_status"],
        table_options=["strict"],
    )
    try:
        with first_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE TABLE probe ("
                    "status TEXT NOT NULL, "
                    "CONSTRAINT uq_probe_status UNIQUE(status), "
                    "CONSTRAINT ck_probe_status CHECK(length(status) > 0)"
                    ")"
                )
            )
            connection.execute(text("CREATE INDEX ix_probe_status ON probe(status)"))
            connection.execute(
                text(
                    "CREATE TRIGGER trg_probe_status BEFORE UPDATE ON probe "
                    "BEGIN SELECT RAISE(ABORT, 'first'); END"
                )
            )
        with second_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE TABLE probe ("
                    "status TEXT NOT NULL, "
                    "CONSTRAINT uq_probe_status UNIQUE(status), "
                    "CONSTRAINT ck_probe_status CHECK(length(status) > 1)"
                    ") STRICT"
                )
            )
            connection.execute(
                text("CREATE UNIQUE INDEX ix_probe_status ON probe(status)")
            )
            connection.execute(
                text(
                    "CREATE TRIGGER trg_probe_status BEFORE UPDATE ON probe "
                    "BEGIN SELECT RAISE(ABORT, 'second'); END"
                )
            )

        first_projection = schema_capability_projection(first_engine, spec)["tables"][
            "probe"
        ]
        second_projection = schema_capability_projection(second_engine, spec)["tables"][
            "probe"
        ]
        assert first_projection["indexes"] != second_projection["indexes"]
        assert (
            first_projection["check_constraints"]
            != (second_projection["check_constraints"])
        )
        assert first_projection["triggers"] != second_projection["triggers"]
        assert first_projection["table_options"] != second_projection["table_options"]
        assert schema_projection_sha256(
            first_engine,
            spec,
        ) != schema_projection_sha256(second_engine, spec)
    finally:
        first_engine.dispose()
        second_engine.dispose()


def test_missing_duplicate_ambiguous_or_unsupported_components_fail_closed(
    tmp_path: Path,
):
    engine = create_database_engine(f"sqlite:///{tmp_path / 'probe.db'}")
    try:
        with engine.begin() as connection:
            connection.execute(text("CREATE TABLE parent (parent_id TEXT PRIMARY KEY)"))
            connection.execute(
                text(
                    "CREATE TABLE probe ("
                    "status TEXT, parent_id TEXT, "
                    "FOREIGN KEY(parent_id) REFERENCES parent(parent_id), "
                    "FOREIGN KEY(parent_id) REFERENCES parent(parent_id)"
                    ")"
                )
            )

        with pytest.raises(PersistenceProjectionError):
            schema_projection_sha256(engine, _probe_spec(columns=["missing"]))
        with pytest.raises(PersistenceProjectionError):
            _probe_spec(columns=["status", "status"])
        with pytest.raises(PersistenceProjectionError):
            _probe_spec(columns=["status"], table_options=["autoincrement"])
        with pytest.raises(PersistenceProjectionError):
            _probe_spec(
                columns=["status"],
                unique_constraints=[
                    {"columns": ["review_note"], "name": "uq_probe_review"}
                ],
            )
        with pytest.raises(PersistenceProjectionError):
            schema_projection_sha256(
                engine,
                _probe_spec(
                    columns=["parent_id"],
                    foreign_keys=[
                        {
                            "columns": ["parent_id"],
                            "referred_columns": ["parent_id"],
                            "referred_table": "parent",
                        }
                    ],
                ),
            )
    finally:
        engine.dispose()


def test_selected_component_removal_fails_closed(tmp_path: Path):
    engine = create_database_engine(f"sqlite:///{tmp_path / 'probe.db'}")
    try:
        with engine.begin() as connection:
            connection.execute(text("CREATE TABLE probe (status TEXT NOT NULL)"))

        missing_unique = _probe_spec(
            columns=["status"],
            unique_constraints=[{"columns": ["status"], "name": "uq_probe_status"}],
        )
        missing_index = _probe_spec(
            columns=["status"],
            indexes=["ix_probe_status"],
        )
        missing_check = _probe_spec(
            columns=["status"],
            check_constraints=["ck_probe_status"],
        )
        missing_trigger = _probe_spec(
            columns=["status"],
            triggers=["trg_probe_status"],
        )

        for spec in (
            missing_unique,
            missing_index,
            missing_check,
            missing_trigger,
        ):
            with pytest.raises(PersistenceProjectionError):
                schema_projection_sha256(engine, spec)
    finally:
        engine.dispose()
