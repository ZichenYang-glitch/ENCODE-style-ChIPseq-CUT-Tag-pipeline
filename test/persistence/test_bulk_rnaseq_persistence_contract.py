"""Compatibility contract tests for Bulk RNA-seq execution persistence."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import re

from alembic.script import ScriptDirectory
from sqlalchemy import inspect, text

from encode_pipeline.persistence import create_database_engine
from encode_pipeline.persistence.migrations import alembic_config, upgrade_database
from encode_pipeline.persistence.models import Base


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONTRACT_PATH = (
    PROJECT_ROOT
    / "src/encode_pipeline/contracts/nfcore_rnaseq"
    / "execution-persistence-contract-1.0.0.json"
)
EXPECTED_TOP_LEVEL_KEYS = {
    "schema_version",
    "contract_id",
    "contract_version",
    "minimum_supported_revision",
    "capabilities",
    "required_revisions",
    "required_schema",
    "schema_projection_sha256",
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
EXPECTED_CONTRACT = {
    "schema_version": "1.0.0",
    "contract_id": "bulk-rnaseq-execution-persistence",
    "contract_version": "1.0.0",
    "minimum_supported_revision": "20260714_07",
    "capabilities": EXPECTED_CAPABILITIES,
    "required_revisions": EXPECTED_REQUIRED_REVISIONS,
    "required_schema": EXPECTED_REQUIRED_SCHEMA,
    "schema_projection_sha256": (
        "4d6dc42e39e4ea07ca597a4e651e88ca18734e45783f1382e97216fca64ac346"
    ),
}


def _load_contract() -> tuple[bytes, dict[str, object]]:
    raw = CONTRACT_PATH.read_bytes()
    return raw, json.loads(raw)


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


def _normalized_sql(value: object) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise AssertionError("reflected SQL must be text")
    return " ".join(value.split())


def _split_top_level_sql_list(value: str) -> list[str]:
    clauses: list[str] = []
    start = 0
    depth = 0
    quote: str | None = None
    index = 0
    while index < len(value):
        character = value[index]
        if quote is not None:
            if character == quote:
                if index + 1 < len(value) and value[index + 1] == quote:
                    index += 1
                else:
                    quote = None
        elif character in {"'", '"', "`"}:
            quote = character
        elif character == "[":
            quote = "]"
        elif character == "(":
            depth += 1
        elif character == ")":
            if depth == 0:
                raise AssertionError("unbalanced table DDL")
            depth -= 1
        elif character == "," and depth == 0:
            clauses.append(value[start:index].strip())
            start = index + 1
        index += 1
    clauses.append(value[start:].strip())
    if quote is not None or depth != 0 or any(not clause for clause in clauses):
        raise AssertionError("invalid table DDL")
    return clauses


def _canonical_create_table_sql(value: object) -> str:
    normalized = _normalized_sql(value)
    if normalized is None:
        raise AssertionError("table DDL is missing")
    opening = normalized.find("(")
    closing = normalized.rfind(")")
    if opening < 1 or closing <= opening:
        raise AssertionError("table DDL is invalid")
    clauses = _split_top_level_sql_list(normalized[opening + 1 : closing])
    table_constraint_prefixes = (
        "CONSTRAINT ",
        "PRIMARY KEY",
        "UNIQUE ",
        "CHECK ",
        "FOREIGN KEY",
    )
    columns: list[str] = []
    constraints: list[str] = []
    for clause in clauses:
        if clause.upper().startswith(table_constraint_prefixes):
            constraints.append(clause)
        else:
            columns.append(clause)
    canonical_body = ", ".join(
        (*columns, *sorted(constraints, key=lambda item: item.upper()))
    )
    return f"{normalized[: opening + 1]}{canonical_body}{normalized[closing:]}"


def _canonical_sort_key(value: object) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def _schema_capability_projection(
    engine,
    required_schema: dict[str, list[str]],
) -> dict[str, object]:
    inspector = inspect(engine)
    tables: dict[str, object] = {}
    with engine.connect() as connection:
        for table_name in sorted(required_schema):
            table_sql = connection.execute(
                text(
                    "SELECT sql FROM sqlite_master "
                    "WHERE type = 'table' AND name = :name"
                ),
                {"name": table_name},
            ).scalar_one()
            columns = [
                {
                    "name": column["name"],
                    "type": str(column["type"]),
                    "nullable": bool(column["nullable"]),
                    "default": _normalized_sql(column["default"]),
                    "primary_key_position": int(column["primary_key"]),
                }
                for column in inspector.get_columns(table_name)
            ]
            primary_key = inspector.get_pk_constraint(table_name)
            foreign_keys = [
                {
                    "columns": list(item["constrained_columns"]),
                    "referred_table": item["referred_table"],
                    "referred_columns": list(item["referred_columns"]),
                    "options": {
                        key: item.get("options", {}).get(key)
                        for key in (
                            "deferrable",
                            "initially",
                            "match",
                            "ondelete",
                            "onupdate",
                        )
                        if item.get("options", {}).get(key) is not None
                    },
                }
                for item in inspector.get_foreign_keys(table_name)
            ]
            unique_constraints = [
                list(item["column_names"])
                for item in inspector.get_unique_constraints(table_name)
            ]
            unique_indexes = [
                {
                    "columns": list(item["column_names"]),
                    "expressions": list(item.get("expressions") or ()),
                    "where": _normalized_sql(
                        str(item.get("dialect_options", {}).get("sqlite_where"))
                    )
                    if item.get("dialect_options", {}).get("sqlite_where") is not None
                    else None,
                }
                for item in inspector.get_indexes(table_name)
                if item.get("unique") is True
            ]
            check_constraints = [
                _normalized_sql(item["sqltext"])
                for item in inspector.get_check_constraints(table_name)
            ]
            triggers = [
                {
                    "name": row.name,
                    "sql": _normalized_sql(row.sql),
                }
                for row in connection.execute(
                    text(
                        "SELECT name, sql FROM sqlite_master "
                        "WHERE type = 'trigger' AND tbl_name = :name "
                        "ORDER BY name"
                    ),
                    {"name": table_name},
                )
            ]
            tables[table_name] = {
                "create_table_sql": _canonical_create_table_sql(table_sql),
                "columns": columns,
                "primary_key": list(primary_key["constrained_columns"]),
                "foreign_keys": sorted(foreign_keys, key=_canonical_sort_key),
                "unique_constraints": sorted(
                    unique_constraints,
                    key=_canonical_sort_key,
                ),
                "unique_indexes": sorted(unique_indexes, key=_canonical_sort_key),
                "check_constraints": sorted(check_constraints),
                "triggers": triggers,
                "table_options": {
                    "autoincrement": "AUTOINCREMENT" in table_sql.upper(),
                    "strict": re.search(r"\bSTRICT\b", table_sql.upper()) is not None,
                    "without_rowid": (
                        re.search(r"\bWITHOUT\s+ROWID\b", table_sql.upper()) is not None
                    ),
                },
            }
    return {
        "schema_version": "bulk-rnaseq-sqlite-capability-projection-v1",
        "tables": tables,
    }


def _schema_projection_sha256(
    engine,
    required_schema: dict[str, list[str]],
) -> str:
    return hashlib.sha256(
        _canonical_sort_key(_schema_capability_projection(engine, required_schema))
    ).hexdigest()


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

    assert set(contract) == EXPECTED_TOP_LEVEL_KEYS
    assert contract == EXPECTED_CONTRACT
    assert raw == _canonical_json_bytes(contract)
    assert "20260714_07" not in contract["required_revisions"]
    assert UNACTIVATED_MANAGED_INPUT_TABLES.isdisjoint(contract["required_schema"])
    for value in _string_leaves(contract):
        assert not value.startswith(("/", "\\", "./", "../", "~"))
        assert not re.match(r"^[A-Za-z]:[\\/]", value)
        assert "://" not in value


def test_declared_revisions_exist_below_the_sole_alembic_head():
    _, contract = _load_contract()
    scripts = ScriptDirectory.from_config(alembic_config("sqlite://"))
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
    required_schema = contract["required_schema"]
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
    required_schema = contract["required_schema"]
    migrated_url = f"sqlite:///{tmp_path / 'migrated.db'}"
    upgrade_database(migrated_url)
    migrated_engine = create_database_engine(migrated_url)
    modeled_engine = create_database_engine(f"sqlite:///{tmp_path / 'modeled.db'}")
    Base.metadata.create_all(modeled_engine)
    try:
        migrated_projection = _schema_capability_projection(
            migrated_engine,
            required_schema,
        )
        modeled_projection = _schema_capability_projection(
            modeled_engine,
            required_schema,
        )

        for table_name, required_columns in required_schema.items():
            migrated_columns = {
                column["name"]: column
                for column in migrated_projection["tables"][table_name]["columns"]
            }
            modeled_columns = {
                column["name"]: column
                for column in modeled_projection["tables"][table_name]["columns"]
            }
            for column_name in required_columns:
                assert {
                    key: migrated_columns[column_name][key]
                    for key in ("type", "nullable", "primary_key_position")
                } == {
                    key: modeled_columns[column_name][key]
                    for key in ("type", "nullable", "primary_key_position")
                }
        expected_digest = contract["schema_projection_sha256"]
        assert (
            _schema_projection_sha256(
                migrated_engine,
                required_schema,
            )
            == expected_digest
        )

        with migrated_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE TABLE reference_catalog_review_only "
                    "(reference_id TEXT PRIMARY KEY)"
                )
            )
        assert (
            _schema_projection_sha256(
                migrated_engine,
                required_schema,
            )
            == expected_digest
        )

        with migrated_engine.begin() as connection:
            connection.execute(
                text(
                    "CREATE TRIGGER bulk_run_status_guard "
                    "BEFORE UPDATE OF status ON runs "
                    "BEGIN SELECT RAISE(ABORT, 'review-only'); END"
                )
            )
        assert (
            _schema_projection_sha256(
                migrated_engine,
                required_schema,
            )
            != expected_digest
        )
    finally:
        migrated_engine.dispose()
        modeled_engine.dispose()


def test_schema_capability_projection_binds_column_ddl_semantics(tmp_path: Path):
    plain_engine = create_database_engine(f"sqlite:///{tmp_path / 'plain.db'}")
    collated_engine = create_database_engine(f"sqlite:///{tmp_path / 'collated.db'}")
    try:
        with plain_engine.begin() as connection:
            connection.execute(text("CREATE TABLE probe (status TEXT NOT NULL)"))
        with collated_engine.begin() as connection:
            connection.execute(
                text("CREATE TABLE probe (status TEXT COLLATE NOCASE NOT NULL)")
            )

        required_schema = {"probe": ["status"]}
        assert _schema_projection_sha256(
            plain_engine,
            required_schema,
        ) != _schema_projection_sha256(
            collated_engine,
            required_schema,
        )
    finally:
        plain_engine.dispose()
        collated_engine.dispose()


def test_table_ddl_projection_canonicalizes_table_constraint_order():
    first = (
        "CREATE TABLE probe ("
        "status TEXT NOT NULL, "
        "CONSTRAINT uq_probe_status UNIQUE (status), "
        "CONSTRAINT ck_probe_status CHECK (length(status) > 0)"
        ")"
    )
    second = (
        "CREATE TABLE probe ("
        "status TEXT NOT NULL, "
        "CONSTRAINT ck_probe_status CHECK (length(status) > 0), "
        "CONSTRAINT uq_probe_status UNIQUE (status)"
        ")"
    )

    assert _canonical_create_table_sql(first) == _canonical_create_table_sql(second)
