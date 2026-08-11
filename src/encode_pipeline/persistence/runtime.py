"""Runtime composition for durable local workflow-platform state."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path
import sqlite3
from urllib.parse import quote

from sqlalchemy import Engine
from sqlalchemy.engine import make_url

from encode_pipeline.persistence.database import (
    create_database_engine,
    create_session_factory,
)
from encode_pipeline.persistence.data_registry import (
    SqlAlchemyDataRegistryRepository,
)
from encode_pipeline.persistence.migrations import upgrade_database
from encode_pipeline.persistence.migration_admission import (
    MigrationAdmissionError,
    VerifiedMigrationExecutionInventory,
    verify_migration_execution_inventory,
)
from encode_pipeline.persistence.input_registry import (
    SqlAlchemyInputRegistryRepository,
)
from encode_pipeline.persistence.repositories import SqlAlchemyRunRepository
from encode_pipeline.persistence.reference_profiles import (
    SqlAlchemyReferenceProfileRepository,
)


DATABASE_URL_ENV = "ENCODE_PIPELINE_DATABASE_URL"


class DatabaseSchemaNotReadyError(RuntimeError):
    """Stable rejection used by supported services before any repository opens."""

    def __init__(self, reason_code: str) -> None:
        self.reason_code = reason_code
        super().__init__(f"database schema admission failed [{reason_code}]")

    def __str__(self) -> str:
        return f"database schema admission failed [{self.reason_code}]"


@dataclass(frozen=True)
class RunPersistence:
    """Owned SQLite resources injected into one API process."""

    database_url: str
    engine: Engine
    repository: SqlAlchemyRunRepository
    data_registry_repository: SqlAlchemyDataRegistryRepository
    input_registry_repository: SqlAlchemyInputRegistryRepository
    reference_profile_repository: SqlAlchemyReferenceProfileRepository

    def close(self) -> None:
        """Release pooled database connections during API shutdown."""
        self.engine.dispose()


def resolve_database_url(database_url: str | None = None) -> str:
    """Return an explicit or environment-configured local platform database URL."""
    configured: str | None
    if database_url is not None:
        configured = database_url
    else:
        configured = os.getenv(DATABASE_URL_ENV)
        if configured is None or not configured.strip():
            configured = f"sqlite:///{Path.home() / '.encode-pipeline' / 'platform.db'}"
    if not isinstance(configured, str) or not configured.strip():
        raise ValueError("database_url must be a non-empty string")
    resolved_url = configured.strip()
    try:
        parsed_url = make_url(resolved_url)
    except Exception as exc:
        raise ValueError("database_url must be a valid SQLAlchemy URL") from exc
    if parsed_url.get_backend_name() != "sqlite" or parsed_url.database in {
        None,
        "",
        ":memory:",
    }:
        raise ValueError("database_url must point to a file-backed SQLite database")
    if not Path(parsed_url.database).is_absolute():
        raise ValueError("database_url must point to an absolute SQLite database path")
    return resolved_url


def open_run_persistence(database_url: str | None = None) -> RunPersistence:
    """Migrate and open the durable repository for one API process."""
    resolved_url = resolve_database_url(database_url)
    upgrade_database(resolved_url)
    return _open_persistence(resolved_url)


def open_existing_run_persistence(
    database_url: str | None = None,
) -> RunPersistence:
    """Open an already prepared database without executing Alembic.

    Supported API and worker services use this entry point.  Only the operator
    database-prepare transaction may migrate, after writers are stopped and a
    verified backup exists.
    """
    resolved_url = resolve_database_url(database_url)
    inventory = _verify_existing_schema(resolved_url)
    persistence = _open_persistence(resolved_url)
    try:
        observed_heads = _read_engine_heads(persistence.engine)
        if observed_heads != inventory.heads:
            raise DatabaseSchemaNotReadyError("DATABASE_SCHEMA_NOT_CURRENT")
    except Exception:
        persistence.close()
        raise
    return persistence


def _open_persistence(resolved_url: str) -> RunPersistence:
    engine = create_database_engine(resolved_url)
    session_factory = create_session_factory(engine)
    repository = SqlAlchemyRunRepository(session_factory)
    data_registry_repository = SqlAlchemyDataRegistryRepository(session_factory)
    input_registry_repository = SqlAlchemyInputRegistryRepository(session_factory)
    reference_profile_repository = SqlAlchemyReferenceProfileRepository(session_factory)
    return RunPersistence(
        database_url=resolved_url,
        engine=engine,
        repository=repository,
        data_registry_repository=data_registry_repository,
        input_registry_repository=input_registry_repository,
        reference_profile_repository=reference_profile_repository,
    )


def _verify_existing_schema(
    database_url: str,
) -> VerifiedMigrationExecutionInventory:
    """Compare a read-only SQLite observation with the admitted sole head."""
    try:
        inventory = verify_migration_execution_inventory()
        database = make_url(database_url).database
        if database in {None, "", ":memory:"}:
            raise DatabaseSchemaNotReadyError("DATABASE_SCHEMA_UNAVAILABLE")
        path = Path(database)
        if path.is_symlink() or not path.is_file():
            raise DatabaseSchemaNotReadyError("DATABASE_SCHEMA_UNAVAILABLE")
        connection = sqlite3.connect(
            f"file:{quote(path.as_posix(), safe='/')}?mode=ro",
            uri=True,
            timeout=5,
        )
        try:
            connection.execute("PRAGMA query_only=ON")
            rows = connection.execute(
                "SELECT version_num FROM alembic_version ORDER BY version_num"
            ).fetchall()
        finally:
            connection.close()
        heads = tuple(row[0] for row in rows)
        if (
            len(inventory.heads) != 1
            or heads != inventory.heads
            or any(not isinstance(head, str) for head in heads)
        ):
            raise DatabaseSchemaNotReadyError("DATABASE_SCHEMA_NOT_CURRENT")
        return inventory
    except DatabaseSchemaNotReadyError:
        raise
    except MigrationAdmissionError:
        raise DatabaseSchemaNotReadyError("DATABASE_SCHEMA_ADMISSION_FAILED") from None
    except (OSError, sqlite3.Error, TypeError, ValueError):
        raise DatabaseSchemaNotReadyError("DATABASE_SCHEMA_UNAVAILABLE") from None


def _read_engine_heads(engine: Engine) -> tuple[str, ...]:
    try:
        with engine.connect() as connection:
            rows = connection.exec_driver_sql(
                "SELECT version_num FROM alembic_version ORDER BY version_num"
            ).all()
        heads = tuple(row[0] for row in rows)
        if any(not isinstance(head, str) for head in heads):
            raise ValueError("invalid database head")
        return heads
    except DatabaseSchemaNotReadyError:
        raise
    except Exception:
        raise DatabaseSchemaNotReadyError("DATABASE_SCHEMA_UNAVAILABLE") from None
