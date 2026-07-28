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
from encode_pipeline.persistence.runtime import (
    DATABASE_URL_ENV,
    RunPersistence,
    open_run_persistence,
    resolve_database_url,
)

__all__ = [
    "DATABASE_URL_ENV",
    "RunPersistence",
    "SqlAlchemyDataRegistryRepository",
    "SqlAlchemyInputRegistryRepository",
    "SqlAlchemyRunRepository",
    "create_database_engine",
    "create_session_factory",
    "open_run_persistence",
    "resolve_database_url",
    "upgrade_database",
]
