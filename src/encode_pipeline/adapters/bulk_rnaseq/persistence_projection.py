"""Contract-directed SQLite capability projection for Bulk RNA-seq."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import re
import sqlite3
from typing import Any

from sqlalchemy import inspect
from sqlalchemy.engine import Engine


SCHEMA_PROJECTION_SCHEME = "bulk-rnaseq-sqlite-capability-projection-v2"
SUPPORTED_TABLE_OPTIONS = frozenset({"strict", "without_rowid"})

_SCHEMA_NAME = re.compile(r"^[a-z][a-z0-9_]*$")
_COMPONENT_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


class PersistenceProjectionError(ValueError):
    """One requested schema component could not be projected safely."""


@dataclass(frozen=True)
class ForeignKeySelector:
    """Semantic selector for exactly one reflected foreign key."""

    columns: tuple[str, ...]
    referred_table: str
    referred_columns: tuple[str, ...]


@dataclass(frozen=True)
class UniqueConstraintSelector:
    """Selector for exactly one reflected unique constraint."""

    name: str | None
    columns: tuple[str, ...]


@dataclass(frozen=True)
class TableProjectionSpec:
    """Explicit schema components used by Bulk execution in one table."""

    columns: tuple[str, ...]
    primary_key: bool
    foreign_keys: tuple[ForeignKeySelector, ...]
    unique_constraints: tuple[UniqueConstraintSelector, ...]
    indexes: tuple[str, ...]
    check_constraints: tuple[str, ...]
    triggers: tuple[str, ...]
    table_options: tuple[str, ...]


@dataclass(frozen=True)
class SchemaProjectionSpec:
    """Canonical source-owned selectors and their expected projection digest."""

    scheme: str
    tables: tuple[tuple[str, TableProjectionSpec], ...]
    sha256: str


@dataclass(frozen=True)
class _SqlToken:
    kind: str
    value: str


@dataclass(frozen=True)
class _ConstraintClause:
    kind: str
    name: str | None
    tokens: tuple[_SqlToken, ...]


def parse_schema_projection_spec(value: object) -> SchemaProjectionSpec:
    """Parse one exact canonical selector document."""

    if not isinstance(value, dict) or set(value) != {"scheme", "sha256", "tables"}:
        raise PersistenceProjectionError("schema projection contract is invalid")
    scheme = value["scheme"]
    expected_sha256 = value["sha256"]
    tables = value["tables"]
    if (
        scheme != SCHEMA_PROJECTION_SCHEME
        or not isinstance(expected_sha256, str)
        or _SHA256.fullmatch(expected_sha256) is None
        or not isinstance(tables, dict)
        or not tables
        or tuple(tables) != tuple(sorted(tables))
    ):
        raise PersistenceProjectionError("schema projection contract is invalid")

    normalized_tables: list[tuple[str, TableProjectionSpec]] = []
    for table_name, table_value in tables.items():
        if (
            not isinstance(table_name, str)
            or _SCHEMA_NAME.fullmatch(table_name) is None
            or not isinstance(table_value, dict)
            or set(table_value)
            != {
                "check_constraints",
                "columns",
                "foreign_keys",
                "indexes",
                "primary_key",
                "table_options",
                "triggers",
                "unique_constraints",
            }
        ):
            raise PersistenceProjectionError("schema projection table is invalid")
        columns = _canonical_name_list(table_value["columns"], schema_names=True)
        if not columns:
            raise PersistenceProjectionError("schema projection columns are invalid")
        foreign_keys = _foreign_key_selectors(table_value["foreign_keys"])
        unique_constraints = _unique_constraint_selectors(
            table_value["unique_constraints"]
        )
        indexes = _canonical_name_list(table_value["indexes"])
        checks = _canonical_name_list(table_value["check_constraints"])
        triggers = _canonical_name_list(table_value["triggers"])
        options = _canonical_name_list(
            table_value["table_options"],
            schema_names=True,
        )
        if not set(options).issubset(SUPPORTED_TABLE_OPTIONS):
            raise PersistenceProjectionError("schema projection option is unsupported")
        primary_key = table_value["primary_key"]
        if not isinstance(primary_key, bool):
            raise PersistenceProjectionError("schema projection primary key is invalid")
        normalized_tables.append(
            (
                table_name,
                TableProjectionSpec(
                    columns=columns,
                    primary_key=primary_key,
                    foreign_keys=foreign_keys,
                    unique_constraints=unique_constraints,
                    indexes=indexes,
                    check_constraints=checks,
                    triggers=triggers,
                    table_options=options,
                ),
            )
        )
    normalized = SchemaProjectionSpec(
        scheme=scheme,
        tables=tuple(normalized_tables),
        sha256=expected_sha256,
    )
    table_specs = dict(normalized.tables)
    for _table_name, table_spec in normalized.tables:
        selected_columns = set(table_spec.columns)
        if any(
            not set(selector.columns).issubset(selected_columns)
            for selector in table_spec.foreign_keys
        ) or any(
            not set(selector.columns).issubset(selected_columns)
            for selector in table_spec.unique_constraints
        ):
            raise PersistenceProjectionError(
                "schema projection constraint column is not selected"
            )
        for selector in table_spec.foreign_keys:
            referred_spec = table_specs.get(selector.referred_table)
            if referred_spec is None or not set(selector.referred_columns).issubset(
                referred_spec.columns
            ):
                raise PersistenceProjectionError(
                    "schema projection foreign key target is not selected"
                )
    return normalized


def schema_projection_spec_document(spec: SchemaProjectionSpec) -> dict[str, Any]:
    """Return the canonical JSON-compatible representation of one spec."""

    if not isinstance(spec, SchemaProjectionSpec):
        raise PersistenceProjectionError("schema projection spec is invalid")
    tables: dict[str, Any] = {}
    for table_name, table in spec.tables:
        tables[table_name] = {
            "check_constraints": list(table.check_constraints),
            "columns": list(table.columns),
            "foreign_keys": [
                {
                    "columns": list(selector.columns),
                    "referred_columns": list(selector.referred_columns),
                    "referred_table": selector.referred_table,
                }
                for selector in table.foreign_keys
            ],
            "indexes": list(table.indexes),
            "primary_key": table.primary_key,
            "table_options": list(table.table_options),
            "triggers": list(table.triggers),
            "unique_constraints": [
                {
                    "columns": list(selector.columns),
                    "name": selector.name,
                }
                for selector in table.unique_constraints
            ],
        }
    return {
        "scheme": spec.scheme,
        "sha256": spec.sha256,
        "tables": tables,
    }


def with_schema_projection_sha256(
    spec: SchemaProjectionSpec,
    sha256: str,
) -> SchemaProjectionSpec:
    """Return the same selectors with one verified projection digest."""

    if not isinstance(sha256, str) or _SHA256.fullmatch(sha256) is None:
        raise PersistenceProjectionError("schema projection digest is invalid")
    return SchemaProjectionSpec(
        scheme=spec.scheme,
        tables=spec.tables,
        sha256=sha256,
    )


def schema_capability_projection(
    engine: Engine,
    spec: SchemaProjectionSpec,
) -> dict[str, Any]:
    """Project only the SQLite components explicitly selected by the contract."""

    if not isinstance(engine, Engine) or not isinstance(spec, SchemaProjectionSpec):
        raise PersistenceProjectionError("schema projection inputs are invalid")
    if engine.dialect.name != "sqlite":
        raise PersistenceProjectionError("schema projection dialect is unsupported")

    inspector = inspect(engine)
    observed_tables = inspector.get_table_names()
    if len(observed_tables) != len(set(observed_tables)):
        raise PersistenceProjectionError("schema projection table is ambiguous")
    selected_tables = {table_name for table_name, _table in spec.tables}
    if not selected_tables.issubset(observed_tables):
        raise PersistenceProjectionError("schema projection table is missing")

    clone = _sqlite_backup(engine)
    try:
        projected_tables = {
            table_name: _project_table(
                inspector,
                engine,
                clone,
                table_name,
                table_spec,
            )
            for table_name, table_spec in spec.tables
        }
    finally:
        clone.close()
    return {
        "scheme": spec.scheme,
        "tables": projected_tables,
    }


def schema_projection_sha256(engine: Engine, spec: SchemaProjectionSpec) -> str:
    """Return the canonical digest for one contract-directed projection."""

    projection = schema_capability_projection(engine, spec)
    return hashlib.sha256(_canonical_json_bytes(projection)).hexdigest()


def _project_table(
    inspector: Any,
    engine: Engine,
    clone: sqlite3.Connection,
    table_name: str,
    spec: TableProjectionSpec,
) -> dict[str, Any]:
    reflected_columns = inspector.get_columns(table_name)
    column_map = _unique_named_components(reflected_columns)
    with engine.connect() as connection:
        pragma_columns = (
            connection.exec_driver_sql(
                f"PRAGMA table_xinfo({_quote_identifier(table_name)})"
            )
            .mappings()
            .all()
        )
        pragma_map = _unique_named_components(pragma_columns)
        table_list_rows = [
            row
            for row in connection.exec_driver_sql("PRAGMA table_list").mappings()
            if row["schema"] == "main"
            and row["name"] == table_name
            and row["type"] == "table"
        ]
    if len(table_list_rows) != 1:
        raise PersistenceProjectionError("schema projection table is ambiguous")
    definition_terms = _table_definition_terms(engine, table_name)
    constraint_clauses = _constraint_clauses(definition_terms)

    projected_columns = []
    for column_name in spec.columns:
        reflected = column_map.get(column_name)
        pragma = pragma_map.get(column_name)
        if reflected is None or pragma is None:
            raise PersistenceProjectionError("schema projection column is missing")
        hidden = pragma.get("hidden")
        if isinstance(hidden, bool) or not isinstance(hidden, int) or hidden != 0:
            raise PersistenceProjectionError("schema projection column is unsupported")
        reflected_type = reflected.get("type")
        declared_type = pragma.get("type")
        nullable = reflected.get("nullable")
        primary_key_position = pragma.get("pk")
        if (
            reflected_type is None
            or not isinstance(declared_type, str)
            or not declared_type.strip()
            or not isinstance(nullable, bool)
            or isinstance(primary_key_position, bool)
            or not isinstance(primary_key_position, int)
            or primary_key_position < 0
        ):
            raise PersistenceProjectionError("schema projection column is invalid")
        projected_columns.append(
            {
                "collation": _column_collation(
                    clone,
                    table_name,
                    column_name,
                    len(projected_columns),
                ),
                "declared_not_null": bool(pragma.get("notnull")),
                "declared_type": _canonical_declared_type(declared_type),
                "default": _sql_scalar(pragma.get("dflt_value")),
                "name": column_name,
                "not_null_conflict_policy": _column_not_null_conflict_policy(
                    _selected_column_clause(definition_terms, column_name),
                    declared_not_null=bool(pragma.get("notnull")),
                ),
                "nullable": nullable,
                "primary_key_position": primary_key_position,
                "reflected_type": str(reflected_type).upper(),
            }
        )
    effective_collations = {
        column["name"]: column["collation"] for column in projected_columns
    }

    primary_key: dict[str, Any] | None = None
    if spec.primary_key:
        value = inspector.get_pk_constraint(table_name)
        primary_key_columns = _component_columns(value)
        if not set(primary_key_columns).issubset(spec.columns):
            raise PersistenceProjectionError(
                "schema projection primary key column is not selected"
            )
        primary_key = _selected_primary_key(
            engine,
            table_name,
            primary_key_columns,
            constraint_clauses,
            effective_collations,
        )
        if not primary_key:
            raise PersistenceProjectionError("schema projection primary key is missing")

    projected_foreign_keys = _select_foreign_keys(
        engine,
        table_name,
        inspector,
        spec.foreign_keys,
    )
    projected_unique_constraints = _select_unique_constraints(
        engine,
        table_name,
        inspector.get_unique_constraints(table_name),
        constraint_clauses,
        spec.unique_constraints,
        effective_collations,
    )
    projected_indexes = _select_indexes(
        engine,
        inspector.get_indexes(table_name),
        table_name,
        spec.indexes,
    )
    projected_checks = _select_checks(
        inspector.get_check_constraints(table_name),
        constraint_clauses,
        frozenset(column_map),
        frozenset(spec.columns),
        spec.check_constraints,
    )
    projected_triggers = _select_triggers(engine, table_name, spec.triggers)
    table_row = table_list_rows[0]
    available_options = {
        "strict": bool(table_row["strict"]),
        "without_rowid": bool(table_row["wr"]),
    }
    return {
        "check_constraints": projected_checks,
        "columns": projected_columns,
        "foreign_keys": projected_foreign_keys,
        "indexes": projected_indexes,
        "primary_key": primary_key,
        "table_options": {
            option: available_options[option] for option in spec.table_options
        },
        "triggers": projected_triggers,
        "unique_constraints": projected_unique_constraints,
    }


def _selected_primary_key(
    engine: Engine,
    table_name: str,
    columns: list[str],
    clauses: tuple[_ConstraintClause, ...],
    effective_collations: dict[str, str],
) -> dict[str, Any]:
    matches = [
        clause
        for clause in clauses
        if clause.kind == "primary_key"
        and _constraint_columns(clause, allow_conflict=True) == tuple(columns)
    ]
    if len(matches) != 1:
        raise PersistenceProjectionError(
            "schema projection primary key is missing or ambiguous"
        )
    return {
        "backing_index": _constraint_backing_index(
            engine,
            table_name,
            _constraint_index_columns(matches[0], effective_collations),
            origin="pk",
            required=False,
        ),
        "columns": columns,
        "definition_sha256": _clause_sha256(matches[0]),
    }


def _select_foreign_keys(
    engine: Engine,
    table_name: str,
    inspector: Any,
    selectors: tuple[ForeignKeySelector, ...],
) -> list[dict[str, Any]]:
    if not selectors:
        return []
    with engine.connect() as connection:
        rows = (
            connection.exec_driver_sql(
                f"PRAGMA foreign_key_list({_quote_identifier(table_name)})"
            )
            .mappings()
            .all()
        )
    grouped: dict[int, list[dict[str, Any]]] = {}
    for row in rows:
        identifier = row.get("id")
        if isinstance(identifier, bool) or not isinstance(identifier, int):
            raise PersistenceProjectionError("schema projection foreign key is invalid")
        grouped.setdefault(identifier, []).append(dict(row))

    normalized = []
    for identifier in sorted(grouped):
        sequence_values = [part.get("seq") for part in grouped[identifier]]
        if any(
            isinstance(value, bool) or not isinstance(value, int)
            for value in sequence_values
        ):
            raise PersistenceProjectionError("schema projection foreign key is invalid")
        parts = sorted(grouped[identifier], key=lambda row: row["seq"])
        if [part["seq"] for part in parts] != list(range(len(parts))):
            raise PersistenceProjectionError("schema projection foreign key is invalid")
        referred_tables = {part.get("table") for part in parts}
        columns = [part.get("from") for part in parts]
        if len(referred_tables) != 1 or any(
            not isinstance(column, str) or not column for column in columns
        ):
            continue
        referred_table = next(iter(referred_tables))
        if not any(
            tuple(columns) == selector.columns
            and referred_table == selector.referred_table
            for selector in selectors
        ):
            continue
        on_updates = {part.get("on_update") for part in parts}
        on_deletes = {part.get("on_delete") for part in parts}
        matches = {part.get("match") for part in parts}
        referred_columns = _resolved_foreign_key_target(
            inspector,
            referred_table,
            [part.get("to") for part in parts],
            local_column_count=len(columns),
        )
        if not any(
            tuple(columns) == selector.columns
            and referred_table == selector.referred_table
            and referred_columns == selector.referred_columns
            for selector in selectors
        ):
            continue
        if len(on_updates) != 1 or len(on_deletes) != 1 or len(matches) != 1:
            raise PersistenceProjectionError("schema projection foreign key is invalid")
        on_update = next(iter(on_updates))
        on_delete = next(iter(on_deletes))
        match = next(iter(matches))
        if (
            not isinstance(referred_table, str)
            or not referred_table
            or not isinstance(on_update, str)
            or not on_update
            or not isinstance(on_delete, str)
            or not on_delete
            or not isinstance(match, str)
            or not match
        ):
            raise PersistenceProjectionError("schema projection foreign key is invalid")
        normalized.append(
            {
                "columns": columns,
                "options": {
                    "match": match.upper(),
                    "on_delete": on_delete.upper(),
                    "on_update": on_update.upper(),
                },
                "referred_columns": list(referred_columns),
                "referred_table": referred_table,
            }
        )
    pragma_matches = []
    for selector in selectors:
        matches = [
            item
            for item in normalized
            if tuple(item["columns"]) == selector.columns
            and item["referred_table"] == selector.referred_table
            and tuple(item["referred_columns"]) == selector.referred_columns
        ]
        if len(matches) != 1:
            raise PersistenceProjectionError(
                "schema projection foreign key is missing or ambiguous"
            )
        pragma_matches.append(matches[0])

    reflected = inspector.get_foreign_keys(table_name)
    reflected_foreign_keys = []
    for item in reflected:
        raw_columns = item.get("constrained_columns")
        raw_referred_table = item.get("referred_table")
        if not isinstance(raw_columns, (list, tuple)) or not any(
            tuple(raw_columns) == selector.columns
            and raw_referred_table == selector.referred_table
            for selector in selectors
        ):
            continue
        columns = _component_columns(item)
        referred_table = raw_referred_table
        referred_columns = _resolved_foreign_key_target(
            inspector,
            referred_table,
            item.get("referred_columns"),
            local_column_count=len(columns),
        )
        if not any(
            tuple(columns) == selector.columns
            and referred_table == selector.referred_table
            and referred_columns == selector.referred_columns
            for selector in selectors
        ):
            continue
        options = item.get("options") or {}
        if (
            not isinstance(referred_table, str)
            or not referred_table
            or not isinstance(options, dict)
            or any(
                key
                not in {
                    "deferrable",
                    "initially",
                    "match",
                    "ondelete",
                    "onupdate",
                }
                and value is not None
                for key, value in options.items()
            )
        ):
            raise PersistenceProjectionError("schema projection foreign key is invalid")
        deferrable = options.get("deferrable")
        initially = options.get("initially")
        if deferrable is not None and not isinstance(deferrable, bool):
            raise PersistenceProjectionError("schema projection foreign key is invalid")
        if initially is not None and (
            not isinstance(initially, str)
            or initially.upper() not in {"DEFERRED", "IMMEDIATE"}
        ):
            raise PersistenceProjectionError("schema projection foreign key is invalid")
        reflected_foreign_keys.append(
            {
                "columns": columns,
                "deferrable": deferrable,
                "initially": initially.upper() if isinstance(initially, str) else None,
                "options": {
                    key: str(options[key]).upper()
                    for key in ("match", "ondelete", "onupdate")
                    if options.get(key) is not None
                },
                "referred_columns": list(referred_columns),
                "referred_table": referred_table,
            }
        )
    selected = []
    for selector, pragma_match in zip(selectors, pragma_matches, strict=True):
        reflected_matches = [
            item
            for item in reflected_foreign_keys
            if tuple(item["columns"]) == selector.columns
            and item["referred_table"] == selector.referred_table
            and tuple(item["referred_columns"]) == selector.referred_columns
        ]
        if len(reflected_matches) != 1:
            raise PersistenceProjectionError(
                "schema projection foreign key is missing or ambiguous"
            )
        reflected_match = reflected_matches[0]
        reflected_options = reflected_match["options"]
        for reflected_key, pragma_key in (
            ("match", "match"),
            ("ondelete", "on_delete"),
            ("onupdate", "on_update"),
        ):
            reflected_value = reflected_options.get(reflected_key)
            if (
                reflected_value is not None
                and reflected_value != pragma_match["options"][pragma_key]
            ):
                raise PersistenceProjectionError(
                    "schema projection foreign key reflection is inconsistent"
                )
        selected.append(
            {
                **pragma_match,
                "options": {
                    **pragma_match["options"],
                    "deferrable": reflected_match["deferrable"],
                    "initially": reflected_match["initially"],
                },
            }
        )
    return selected


def _resolved_foreign_key_target(
    inspector: Any,
    referred_table: object,
    values: object,
    *,
    local_column_count: int,
) -> tuple[str, ...]:
    if (
        not isinstance(referred_table, str)
        or not referred_table
        or not isinstance(values, (list, tuple))
        or len(values) != local_column_count
    ):
        raise PersistenceProjectionError("schema projection foreign key is invalid")
    if all(isinstance(value, str) and value for value in values):
        return tuple(values)
    if not values or not all(value is None for value in values):
        raise PersistenceProjectionError("schema projection foreign key is invalid")
    try:
        primary_key = inspector.get_pk_constraint(referred_table)
        columns = tuple(_component_columns(primary_key))
    except Exception:
        raise PersistenceProjectionError(
            "schema projection foreign key target is invalid"
        ) from None
    if len(columns) != local_column_count:
        raise PersistenceProjectionError(
            "schema projection foreign key target is invalid"
        )
    return columns


def _select_unique_constraints(
    engine: Engine,
    table_name: str,
    observed: list[dict[str, Any]],
    clauses: tuple[_ConstraintClause, ...],
    selectors: tuple[UniqueConstraintSelector, ...],
    effective_collations: dict[str, str],
) -> list[dict[str, Any]]:
    normalized = [
        {
            "columns": _component_columns(item),
            "name": _optional_component_name(item.get("name")),
        }
        for item in observed
    ]
    selected = []
    for selector in selectors:
        matches = [
            item
            for item in normalized
            if item["name"] == selector.name
            and tuple(item["columns"]) == selector.columns
        ]
        if len(matches) != 1:
            raise PersistenceProjectionError(
                "schema projection unique constraint is missing or ambiguous"
            )
        clause_matches = [
            clause
            for clause in clauses
            if clause.kind == "unique"
            and clause.name == selector.name
            and _constraint_columns(clause, allow_conflict=True) == selector.columns
        ]
        if len(clause_matches) != 1:
            raise PersistenceProjectionError(
                "schema projection unique constraint is missing or ambiguous"
            )
        selected.append(
            {
                **matches[0],
                "backing_index": _constraint_backing_index(
                    engine,
                    table_name,
                    _constraint_index_columns(
                        clause_matches[0],
                        effective_collations,
                    ),
                    origin="u",
                    required=True,
                ),
                "definition_sha256": _clause_sha256(clause_matches[0]),
            }
        )
    return selected


def _select_indexes(
    engine: Engine,
    observed: list[dict[str, Any]],
    table_name: str,
    selectors: tuple[str, ...],
) -> list[dict[str, Any]]:
    index_map = _selected_named_components(observed, selectors)
    selected = []
    with engine.connect() as connection:
        for name in selectors:
            item = index_map.get(name)
            if item is None:
                raise PersistenceProjectionError("schema projection index is missing")
            columns = _component_columns(item)
            expressions = item.get("expressions") or ()
            if expressions and tuple(expressions) != tuple(columns):
                raise PersistenceProjectionError(
                    "schema projection expression index is unsupported"
                )
            pragma_rows = (
                connection.exec_driver_sql(
                    f"PRAGMA index_xinfo({_quote_identifier(name)})"
                )
                .mappings()
                .all()
            )
            key_rows = [
                row
                for row in pragma_rows
                if row.get("key") == 1 and row.get("cid") is not None
            ]
            if len(key_rows) != len(columns):
                raise PersistenceProjectionError("schema projection index is invalid")
            index_columns = []
            for expected_name, row in zip(columns, key_rows, strict=True):
                if row.get("name") != expected_name:
                    raise PersistenceProjectionError(
                        "schema projection index is ambiguous"
                    )
                collation = row.get("coll")
                if not isinstance(collation, str) or not collation:
                    raise PersistenceProjectionError(
                        "schema projection index collation is invalid"
                    )
                index_columns.append(
                    {
                        "collation": collation.upper(),
                        "descending": bool(row.get("desc")),
                        "name": expected_name,
                    }
                )
            dialect_options = item.get("dialect_options") or {}
            if not isinstance(dialect_options, dict):
                raise PersistenceProjectionError("schema projection index is invalid")
            where = dialect_options.get("sqlite_where")
            selected.append(
                {
                    "columns": index_columns,
                    "name": name,
                    "partial_where": _sql_scalar(str(where))
                    if where is not None
                    else None,
                    "table": table_name,
                    "unique": bool(item.get("unique")),
                }
            )
    return selected


def _select_checks(
    observed: list[dict[str, Any]],
    clauses: tuple[_ConstraintClause, ...],
    observed_columns: frozenset[str],
    selected_columns: frozenset[str],
    selectors: tuple[str, ...],
) -> list[dict[str, str]]:
    by_name = _selected_named_components(observed, selectors)
    selected = []
    for name in selectors:
        item = by_name.get(name)
        if item is None:
            raise PersistenceProjectionError(
                "schema projection check constraint is missing"
            )
        sqltext = _sql_scalar(item.get("sqltext"))
        if sqltext is None:
            raise PersistenceProjectionError(
                "schema projection check constraint is invalid"
            )
        clause_matches = [
            clause
            for clause in clauses
            if clause.kind == "check" and clause.name == name
        ]
        if len(clause_matches) != 1:
            raise PersistenceProjectionError(
                "schema projection check constraint is missing or ambiguous"
            )
        participants = _constraint_participant_columns(
            clause_matches[0],
            observed_columns,
        )
        if not participants.issubset(selected_columns):
            raise PersistenceProjectionError(
                "schema projection check constraint column is not selected"
            )
        selected.append(
            {
                "definition_sha256": _clause_sha256(clause_matches[0]),
                "name": name,
                "sqltext": sqltext,
            }
        )
    return selected


def _constraint_backing_index(
    engine: Engine,
    table_name: str,
    expected_columns: tuple[dict[str, Any], ...],
    *,
    origin: str,
    required: bool,
) -> dict[str, Any] | None:
    candidates = []
    with engine.connect() as connection:
        index_rows = (
            connection.exec_driver_sql(
                f"PRAGMA index_list({_quote_identifier(table_name)})"
            )
            .mappings()
            .all()
        )
        for index_row in index_rows:
            if (
                index_row.get("unique") != 1
                or index_row.get("origin") != origin
                or index_row.get("partial") != 0
            ):
                continue
            index_name = index_row.get("name")
            if not isinstance(index_name, str) or not index_name:
                raise PersistenceProjectionError(
                    "schema projection constraint index is invalid"
                )
            key_rows = [
                row
                for row in connection.exec_driver_sql(
                    f"PRAGMA index_xinfo({_quote_identifier(index_name)})"
                ).mappings()
                if row.get("key") == 1
            ]
            projected_columns = []
            for row in key_rows:
                name = row.get("name")
                collation = row.get("coll")
                descending = row.get("desc")
                column_id = row.get("cid")
                if (
                    not isinstance(name, str)
                    or not name
                    or not isinstance(collation, str)
                    or not collation
                    or descending not in (0, 1)
                    or isinstance(column_id, bool)
                    or not isinstance(column_id, int)
                    or column_id < 0
                ):
                    raise PersistenceProjectionError(
                        "schema projection constraint index is unsupported"
                    )
                projected_columns.append(
                    {
                        "collation": collation.upper(),
                        "descending": bool(descending),
                        "name": name,
                    }
                )
            if tuple(projected_columns) == expected_columns:
                candidates.append(
                    {
                        "columns": projected_columns,
                        "origin": origin,
                    }
                )
    if required and not candidates:
        raise PersistenceProjectionError(
            "schema projection constraint index is missing"
        )
    if not candidates:
        return None
    return {
        "columns": [dict(column) for column in expected_columns],
        "origin": origin,
    }


def _table_definition_terms(
    engine: Engine,
    table_name: str,
) -> tuple[tuple[_SqlToken, ...], ...]:
    with engine.connect() as connection:
        rows = (
            connection.exec_driver_sql(
                "SELECT sql FROM sqlite_schema WHERE type = 'table' AND name = ?",
                (table_name,),
            )
            .mappings()
            .all()
        )
    if len(rows) != 1 or not isinstance(rows[0].get("sql"), str):
        raise PersistenceProjectionError("schema projection table SQL is unavailable")
    tokens = _tokenize_sql(rows[0]["sql"])
    if (
        len(tokens) < 5
        or not _is_keyword(tokens[0], "CREATE")
        or not _is_keyword(tokens[1], "TABLE")
        or _sql_identifier(tokens[2]) != table_name
        or tokens[3] != _SqlToken("symbol", "(")
    ):
        raise PersistenceProjectionError("schema projection table SQL is unsupported")
    closing = _matching_parenthesis(tokens, 3)
    suffix = tokens[closing + 1 :]
    if suffix and not _supported_table_suffix(suffix):
        raise PersistenceProjectionError("schema projection table SQL is unsupported")
    return _split_sql_terms(tokens[4:closing])


def _constraint_clauses(
    definition_terms: tuple[tuple[_SqlToken, ...], ...],
) -> tuple[_ConstraintClause, ...]:
    clauses = []
    for clause_tokens in definition_terms:
        clause = _classify_constraint_clause(clause_tokens)
        if clause is not None:
            clauses.append(clause)
    return tuple(clauses)


def _selected_column_clause(
    definition_terms: tuple[tuple[_SqlToken, ...], ...],
    column_name: str,
) -> tuple[_SqlToken, ...]:
    matches = []
    for clause_tokens in definition_terms:
        if _classify_constraint_clause(clause_tokens) is not None:
            continue
        try:
            observed_name = _sql_identifier(clause_tokens[0])
        except (IndexError, PersistenceProjectionError):
            continue
        if observed_name == column_name:
            matches.append(clause_tokens)
    if len(matches) != 1:
        raise PersistenceProjectionError(
            "schema projection column clause is missing or ambiguous"
        )
    return matches[0]


def _column_not_null_conflict_policy(
    clause_tokens: tuple[_SqlToken, ...],
    *,
    declared_not_null: bool,
) -> str | None:
    positions = []
    depth = 0
    for index, token in enumerate(clause_tokens):
        if token == _SqlToken("symbol", "("):
            depth += 1
        elif token == _SqlToken("symbol", ")"):
            depth -= 1
            if depth < 0:
                raise PersistenceProjectionError(
                    "schema projection column clause is invalid"
                )
        elif (
            depth == 0
            and _is_keyword(token, "NOT")
            and index + 1 < len(clause_tokens)
            and _is_keyword(clause_tokens[index + 1], "NULL")
        ):
            positions.append(index)
    if depth != 0 or len(positions) > 1:
        raise PersistenceProjectionError("schema projection column clause is ambiguous")
    if not positions:
        return None
    if not declared_not_null:
        raise PersistenceProjectionError(
            "schema projection column nullability is inconsistent"
        )
    cursor = positions[0] + 2
    if cursor >= len(clause_tokens) or not _is_keyword(
        clause_tokens[cursor],
        "ON",
    ):
        return "ABORT"
    if cursor + 2 >= len(clause_tokens) or not _is_keyword(
        clause_tokens[cursor + 1], "CONFLICT"
    ):
        raise PersistenceProjectionError(
            "schema projection column conflict policy is unsupported"
        )
    policy = clause_tokens[cursor + 2]
    for supported in ("ROLLBACK", "ABORT", "FAIL", "IGNORE", "REPLACE"):
        if _is_keyword(policy, supported):
            return supported
    raise PersistenceProjectionError(
        "schema projection column conflict policy is unsupported"
    )


def _classify_constraint_clause(
    tokens: tuple[_SqlToken, ...],
) -> _ConstraintClause | None:
    if not tokens:
        raise PersistenceProjectionError("schema projection table clause is invalid")
    offset = 0
    name = None
    if _is_keyword(tokens[0], "CONSTRAINT"):
        if len(tokens) < 3:
            raise PersistenceProjectionError(
                "schema projection constraint clause is invalid"
            )
        try:
            name = _sql_identifier(tokens[1])
        except PersistenceProjectionError:
            return None
        offset = 2
    remaining = tokens[offset:]
    if (
        len(remaining) >= 2
        and _is_keyword(remaining[0], "PRIMARY")
        and _is_keyword(remaining[1], "KEY")
    ):
        kind = "primary_key"
    elif remaining and _is_keyword(remaining[0], "UNIQUE"):
        kind = "unique"
    elif (
        len(remaining) >= 2
        and _is_keyword(remaining[0], "FOREIGN")
        and _is_keyword(remaining[1], "KEY")
    ):
        kind = "foreign_key"
    elif remaining and _is_keyword(remaining[0], "CHECK"):
        kind = "check"
    else:
        return None
    return _ConstraintClause(kind=kind, name=name, tokens=tokens)


def _constraint_columns(
    clause: _ConstraintClause,
    *,
    allow_conflict: bool,
) -> tuple[str, ...]:
    tokens = clause.tokens
    offset = 2 if tokens and _is_keyword(tokens[0], "CONSTRAINT") else 0
    if clause.kind == "primary_key":
        offset += 2
    elif clause.kind == "unique":
        offset += 1
    else:
        raise PersistenceProjectionError(
            "schema projection constraint clause is unsupported"
        )
    if offset >= len(tokens) or tokens[offset] != _SqlToken("symbol", "("):
        raise PersistenceProjectionError(
            "schema projection constraint clause is invalid"
        )
    closing = _matching_parenthesis(tokens, offset)
    columns = []
    for term in _split_sql_terms(tokens[offset + 1 : closing]):
        if not term:
            raise PersistenceProjectionError(
                "schema projection constraint column is invalid"
            )
        column = _sql_identifier(term[0])
        cursor = 1
        if cursor + 1 < len(term) and _is_keyword(term[cursor], "COLLATE"):
            _sql_identifier(term[cursor + 1])
            cursor += 2
        if cursor < len(term) and (
            _is_keyword(term[cursor], "ASC") or _is_keyword(term[cursor], "DESC")
        ):
            cursor += 1
        if cursor != len(term):
            raise PersistenceProjectionError(
                "schema projection constraint column is unsupported"
            )
        columns.append(column)
    tail = tokens[closing + 1 :]
    if tail:
        if (
            not allow_conflict
            or len(tail) != 3
            or not _is_keyword(tail[0], "ON")
            or not _is_keyword(tail[1], "CONFLICT")
            or not any(
                _is_keyword(tail[2], policy)
                for policy in ("ROLLBACK", "ABORT", "FAIL", "IGNORE", "REPLACE")
            )
        ):
            raise PersistenceProjectionError(
                "schema projection constraint clause is unsupported"
            )
    if not columns or len(columns) != len(set(columns)):
        raise PersistenceProjectionError(
            "schema projection constraint column is invalid"
        )
    return tuple(columns)


def _constraint_index_columns(
    clause: _ConstraintClause,
    effective_collations: dict[str, str],
) -> tuple[dict[str, Any], ...]:
    columns = _constraint_columns(clause, allow_conflict=True)
    tokens = clause.tokens
    offset = 2 if tokens and _is_keyword(tokens[0], "CONSTRAINT") else 0
    offset += 2 if clause.kind == "primary_key" else 1
    closing = _matching_parenthesis(tokens, offset)
    terms = _split_sql_terms(tokens[offset + 1 : closing])
    projected = []
    for column, term in zip(columns, terms, strict=True):
        collation = effective_collations.get(column)
        if not isinstance(collation, str) or not collation:
            raise PersistenceProjectionError(
                "schema projection constraint column is not selected"
            )
        cursor = 1
        if cursor + 1 < len(term) and _is_keyword(term[cursor], "COLLATE"):
            collation = _sql_identifier(term[cursor + 1]).upper()
            cursor += 2
        descending = False
        if cursor < len(term) and (
            _is_keyword(term[cursor], "ASC") or _is_keyword(term[cursor], "DESC")
        ):
            descending = _is_keyword(term[cursor], "DESC")
            cursor += 1
        if cursor != len(term):
            raise PersistenceProjectionError(
                "schema projection constraint column is unsupported"
            )
        projected.append(
            {
                "collation": collation,
                "descending": descending,
                "name": column,
            }
        )
    return tuple(projected)


def _constraint_participant_columns(
    clause: _ConstraintClause,
    observed_columns: frozenset[str],
) -> frozenset[str]:
    return frozenset(
        identifier
        for token in clause.tokens
        if token.kind in {"word", "identifier"}
        for identifier in (_sql_identifier(token),)
        if identifier in observed_columns
    )


def _clause_sha256(clause: _ConstraintClause) -> str:
    return hashlib.sha256(
        _canonical_json_bytes(
            [
                {
                    "kind": token.kind,
                    "value": token.value,
                }
                for token in clause.tokens
            ]
        )
    ).hexdigest()


def _tokenize_sql(value: str) -> tuple[_SqlToken, ...]:
    tokens: list[_SqlToken] = []
    cursor = 0
    while cursor < len(value):
        character = value[cursor]
        if character.isspace():
            cursor += 1
            continue
        if value.startswith("--", cursor):
            newline = value.find("\n", cursor + 2)
            cursor = len(value) if newline < 0 else newline + 1
            continue
        if value.startswith("/*", cursor):
            closing = value.find("*/", cursor + 2)
            if closing < 0:
                raise PersistenceProjectionError(
                    "schema projection SQL comment is invalid"
                )
            cursor = closing + 2
            continue
        if character in {"'", '"', "`", "["}:
            token, cursor = _quoted_sql_token(value, cursor)
            tokens.append(token)
            continue
        if character.isalpha() or character == "_":
            end = cursor + 1
            while end < len(value) and (
                value[end].isalnum() or value[end] in {"_", "$"}
            ):
                end += 1
            tokens.append(_SqlToken("word", value[cursor:end].upper()))
            cursor = end
            continue
        if character.isdigit() or (
            character == "." and cursor + 1 < len(value) and value[cursor + 1].isdigit()
        ):
            end = cursor + 1
            while end < len(value) and (
                value[end].isalnum() or value[end] in {".", "_", "+", "-"}
            ):
                end += 1
            tokens.append(_SqlToken("number", value[cursor:end].upper()))
            cursor = end
            continue
        two_character = value[cursor : cursor + 2]
        if two_character in {"<=", ">=", "!=", "<>", "==", "||", "<<", ">>"}:
            tokens.append(_SqlToken("symbol", two_character))
            cursor += 2
            continue
        if character in "(),.+-*/%<>=~&|":
            tokens.append(_SqlToken("symbol", character))
            cursor += 1
            continue
        raise PersistenceProjectionError("schema projection SQL token is unsupported")
    if not tokens:
        raise PersistenceProjectionError("schema projection SQL is invalid")
    return tuple(tokens)


def _quoted_sql_token(value: str, start: int) -> tuple[_SqlToken, int]:
    opening = value[start]
    closing = "]" if opening == "[" else opening
    kind = "string" if opening == "'" else "identifier"
    characters = []
    cursor = start + 1
    while cursor < len(value):
        character = value[cursor]
        if character == closing:
            if cursor + 1 < len(value) and value[cursor + 1] == closing:
                characters.append(closing)
                cursor += 2
                continue
            normalized = "".join(characters)
            if kind == "identifier":
                normalized = normalized.casefold()
            return _SqlToken(kind, normalized), cursor + 1
        if character == "\x00":
            break
        characters.append(character)
        cursor += 1
    raise PersistenceProjectionError("schema projection quoted SQL is invalid")


def _matching_parenthesis(tokens: tuple[_SqlToken, ...], opening: int) -> int:
    if opening >= len(tokens) or tokens[opening] != _SqlToken("symbol", "("):
        raise PersistenceProjectionError("schema projection SQL nesting is invalid")
    depth = 0
    for index in range(opening, len(tokens)):
        token = tokens[index]
        if token == _SqlToken("symbol", "("):
            depth += 1
        elif token == _SqlToken("symbol", ")"):
            depth -= 1
            if depth == 0:
                return index
            if depth < 0:
                break
    raise PersistenceProjectionError("schema projection SQL nesting is invalid")


def _split_sql_terms(
    tokens: tuple[_SqlToken, ...],
) -> tuple[tuple[_SqlToken, ...], ...]:
    terms = []
    start = 0
    depth = 0
    for index, token in enumerate(tokens):
        if token == _SqlToken("symbol", "("):
            depth += 1
        elif token == _SqlToken("symbol", ")"):
            depth -= 1
            if depth < 0:
                raise PersistenceProjectionError(
                    "schema projection SQL nesting is invalid"
                )
        elif token == _SqlToken("symbol", ",") and depth == 0:
            terms.append(tokens[start:index])
            start = index + 1
    if depth != 0:
        raise PersistenceProjectionError("schema projection SQL nesting is invalid")
    terms.append(tokens[start:])
    if any(not term for term in terms):
        raise PersistenceProjectionError("schema projection SQL list is invalid")
    return tuple(terms)


def _sql_identifier(token: _SqlToken) -> str:
    if token.kind == "word":
        identifier = token.value.casefold()
    elif token.kind == "identifier":
        identifier = token.value
    else:
        raise PersistenceProjectionError("schema projection identifier is invalid")
    if _SCHEMA_NAME.fullmatch(identifier) is None:
        raise PersistenceProjectionError("schema projection identifier is unsupported")
    return identifier


def _is_keyword(token: _SqlToken, keyword: str) -> bool:
    return token.kind == "word" and token.value == keyword


def _supported_table_suffix(tokens: tuple[_SqlToken, ...]) -> bool:
    normalized = tuple((token.kind, token.value) for token in tokens)
    return normalized in {
        (("word", "STRICT"),),
        (("word", "WITHOUT"), ("word", "ROWID")),
        (
            ("word", "WITHOUT"),
            ("word", "ROWID"),
            ("symbol", ","),
            ("word", "STRICT"),
        ),
        (
            ("word", "STRICT"),
            ("symbol", ","),
            ("word", "WITHOUT"),
            ("word", "ROWID"),
        ),
    }


def _select_triggers(
    engine: Engine,
    table_name: str,
    selectors: tuple[str, ...],
) -> list[dict[str, str]]:
    selected = []
    with engine.connect() as connection:
        for name in selectors:
            rows = (
                connection.exec_driver_sql(
                    "SELECT name, sql FROM sqlite_schema "
                    "WHERE type = 'trigger' AND tbl_name = ? AND name = ?",
                    (table_name, name),
                )
                .mappings()
                .all()
            )
            if len(rows) != 1:
                raise PersistenceProjectionError(
                    "schema projection trigger is missing or ambiguous"
                )
            sql = _sql_scalar(rows[0].get("sql"))
            if sql is None:
                raise PersistenceProjectionError("schema projection trigger is invalid")
            selected.append({"name": name, "sql": sql})
    return selected


def _sqlite_backup(engine: Engine) -> sqlite3.Connection:
    clone = sqlite3.connect(":memory:")
    connection = engine.raw_connection()
    try:
        raw = getattr(connection, "driver_connection", None)
        if not isinstance(raw, sqlite3.Connection):
            raise PersistenceProjectionError(
                "schema projection connection is unsupported"
            )
        raw.backup(clone)
    except (PersistenceProjectionError, sqlite3.Error):
        clone.close()
        raise
    finally:
        connection.close()
    return clone


def _column_collation(
    clone: sqlite3.Connection,
    table_name: str,
    column_name: str,
    ordinal: int,
) -> str:
    try:
        object_rows = clone.execute(
            "SELECT name FROM sqlite_schema WHERE name IS NOT NULL"
        ).fetchall()
        if any(
            len(row) != 1 or not isinstance(row[0], str) or not row[0]
            for row in object_rows
        ):
            raise PersistenceProjectionError(
                "schema projection object inventory is invalid"
            )
        names = {row[0].casefold() for row in object_rows}
        prefix = f"helixweave_projection_{ordinal}_{table_name}_{column_name}"
        index_name = next(
            (
                f"{prefix}_{candidate}"
                for candidate in range(len(names) + 1)
                if f"{prefix}_{candidate}".casefold() not in names
            ),
            None,
        )
        if index_name is None:
            raise PersistenceProjectionError(
                "schema projection probe name is unavailable"
            )
        clone.execute(
            f"CREATE INDEX {_quote_identifier(index_name)} "
            f"ON {_quote_identifier(table_name)}({_quote_identifier(column_name)})"
        )
        rows = clone.execute(
            f"PRAGMA index_xinfo({_quote_identifier(index_name)})"
        ).fetchall()
    except sqlite3.Error as error:
        raise PersistenceProjectionError(
            "schema projection collation could not be interpreted"
        ) from error
    matches = [
        row for row in rows if len(row) == 6 and row[5] == 1 and row[2] == column_name
    ]
    if len(matches) != 1 or not isinstance(matches[0][4], str) or not matches[0][4]:
        raise PersistenceProjectionError("schema projection collation is ambiguous")
    return matches[0][4].upper()


def _foreign_key_selectors(value: object) -> tuple[ForeignKeySelector, ...]:
    if not isinstance(value, list):
        raise PersistenceProjectionError("schema projection foreign keys are invalid")
    selectors = []
    documents = []
    for item in value:
        if not isinstance(item, dict) or set(item) != {
            "columns",
            "referred_columns",
            "referred_table",
        }:
            raise PersistenceProjectionError(
                "schema projection foreign key selector is invalid"
            )
        columns = _ordered_schema_names(item["columns"])
        referred_columns = _ordered_schema_names(item["referred_columns"])
        referred_table = item["referred_table"]
        if (
            not columns
            or len(columns) != len(referred_columns)
            or not isinstance(referred_table, str)
            or _SCHEMA_NAME.fullmatch(referred_table) is None
        ):
            raise PersistenceProjectionError(
                "schema projection foreign key selector is invalid"
            )
        selector = ForeignKeySelector(
            columns=columns,
            referred_table=referred_table,
            referred_columns=referred_columns,
        )
        selectors.append(selector)
        documents.append(
            {
                "columns": list(columns),
                "referred_columns": list(referred_columns),
                "referred_table": referred_table,
            }
        )
    _require_canonical_documents(documents)
    return tuple(selectors)


def _unique_constraint_selectors(
    value: object,
) -> tuple[UniqueConstraintSelector, ...]:
    if not isinstance(value, list):
        raise PersistenceProjectionError(
            "schema projection unique constraints are invalid"
        )
    selectors = []
    documents = []
    for item in value:
        if not isinstance(item, dict) or set(item) != {"columns", "name"}:
            raise PersistenceProjectionError(
                "schema projection unique constraint selector is invalid"
            )
        columns = _ordered_schema_names(item["columns"])
        name = item["name"]
        if not columns or (
            name is not None
            and (not isinstance(name, str) or _COMPONENT_NAME.fullmatch(name) is None)
        ):
            raise PersistenceProjectionError(
                "schema projection unique constraint selector is invalid"
            )
        selector = UniqueConstraintSelector(name=name, columns=columns)
        selectors.append(selector)
        documents.append({"columns": list(columns), "name": name})
    _require_canonical_documents(documents)
    return tuple(selectors)


def _canonical_name_list(
    value: object,
    *,
    schema_names: bool = False,
) -> tuple[str, ...]:
    if not isinstance(value, list):
        raise PersistenceProjectionError("schema projection names are invalid")
    pattern = _SCHEMA_NAME if schema_names else _COMPONENT_NAME
    if any(
        not isinstance(item, str) or pattern.fullmatch(item) is None for item in value
    ):
        raise PersistenceProjectionError("schema projection name is invalid")
    if value != sorted(set(value)):
        raise PersistenceProjectionError(
            "schema projection names are duplicate or non-canonical"
        )
    return tuple(value)


def _ordered_schema_names(value: object) -> tuple[str, ...]:
    if (
        not isinstance(value, list)
        or not value
        or any(
            not isinstance(item, str) or _SCHEMA_NAME.fullmatch(item) is None
            for item in value
        )
        or len(value) != len(set(value))
    ):
        raise PersistenceProjectionError("schema projection columns are invalid")
    return tuple(value)


def _require_canonical_documents(documents: list[dict[str, Any]]) -> None:
    encoded = [_canonical_json_bytes(item) for item in documents]
    if encoded != sorted(set(encoded)):
        raise PersistenceProjectionError(
            "schema projection selectors are duplicate or non-canonical"
        )


def _unique_named_components(
    values: Any,
) -> dict[str, Any]:
    if not isinstance(values, (list, tuple)):
        values = list(values)
    result: dict[str, Any] = {}
    for item in values:
        if not isinstance(item, dict):
            item = dict(item)
        name = item.get("name")
        if not isinstance(name, str) or not name or name in result:
            raise PersistenceProjectionError(
                "schema projection component is unnamed or ambiguous"
            )
        result[name] = item
    return result


def _selected_named_components(
    values: Any,
    selected_names: tuple[str, ...],
) -> dict[str, Any]:
    selected = set(selected_names)
    matches: dict[str, list[dict[str, Any]]] = {name: [] for name in selected_names}
    for item in values:
        if not isinstance(item, dict):
            item = dict(item)
        name = item.get("name")
        if isinstance(name, str) and name in selected:
            matches[name].append(item)
    if any(len(items) != 1 for items in matches.values()):
        raise PersistenceProjectionError(
            "schema projection component is missing or ambiguous"
        )
    return {name: matches[name][0] for name in selected_names}


def _component_columns(
    value: dict[str, Any],
    *,
    key: str = "constrained_columns",
) -> list[str]:
    if key == "constrained_columns" and key not in value:
        key = "column_names"
    columns = value.get(key)
    if (
        not isinstance(columns, (list, tuple))
        or not columns
        or any(not isinstance(column, str) or not column for column in columns)
    ):
        raise PersistenceProjectionError("schema projection component is invalid")
    return list(columns)


def _optional_component_name(value: object) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str) or not value:
        raise PersistenceProjectionError("schema projection component name is invalid")
    return value


def _canonical_declared_type(value: str) -> str:
    normalized = " ".join(value.split()).upper()
    if not normalized:
        raise PersistenceProjectionError("schema projection type is invalid")
    return normalized


def _sql_scalar(value: object) -> str | None:
    if value is None:
        return None
    if not isinstance(value, str):
        raise PersistenceProjectionError("schema projection SQL value is invalid")
    normalized = value.strip()
    if not normalized or "\x00" in normalized:
        raise PersistenceProjectionError("schema projection SQL value is invalid")
    return normalized


def _quote_identifier(value: str) -> str:
    if not isinstance(value, str) or _COMPONENT_NAME.fullmatch(value) is None:
        raise PersistenceProjectionError("schema projection identifier is invalid")
    return f'"{value}"'


def _canonical_json_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
