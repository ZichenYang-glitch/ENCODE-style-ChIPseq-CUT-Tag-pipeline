"""SQL persistence adapters for the workflow platform."""

from encode_pipeline.persistence.database import (
    create_database_engine,
    create_session_factory,
)
from encode_pipeline.persistence.data_registry import (
    SqlAlchemyDataRegistryRepository,
)
from encode_pipeline.persistence.input_registry import (
    SqlAlchemyInputRegistryRepository,
)
from encode_pipeline.persistence.migrations import upgrade_database
from encode_pipeline.persistence.repositories import SqlAlchemyRunRepository
from encode_pipeline.persistence.reference_profiles import (
    SqlAlchemyReferenceProfileRepository,
)
from encode_pipeline.persistence.runtime import (
    DATABASE_URL_ENV,
    DatabaseSchemaNotReadyError,
    RunPersistence,
    open_existing_run_persistence,
    open_run_persistence,
    resolve_database_url,
)

__all__ = [
    "DATABASE_URL_ENV",
    "DatabaseSchemaNotReadyError",
    "RunPersistence",
    "SqlAlchemyDataRegistryRepository",
    "SqlAlchemyInputRegistryRepository",
    "SqlAlchemyReferenceProfileRepository",
    "SqlAlchemyRunRepository",
    "create_database_engine",
    "create_session_factory",
    "open_run_persistence",
    "open_existing_run_persistence",
    "resolve_database_url",
    "upgrade_database",
]
