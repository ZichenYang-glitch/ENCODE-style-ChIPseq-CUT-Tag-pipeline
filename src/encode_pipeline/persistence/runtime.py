"""Runtime composition for durable local workflow-platform state."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from sqlalchemy import Engine
from sqlalchemy.engine import make_url

from encode_pipeline.persistence.database import (
    create_database_engine,
    create_session_factory,
)
from encode_pipeline.persistence.migrations import upgrade_database
from encode_pipeline.persistence.repositories import SqlAlchemyRunRepository


DATABASE_URL_ENV = "ENCODE_PIPELINE_DATABASE_URL"


@dataclass(frozen=True)
class RunPersistence:
    """Owned SQLite resources injected into one API process."""

    database_url: str
    engine: Engine
    repository: SqlAlchemyRunRepository

    def close(self) -> None:
        """Release pooled database connections during API shutdown."""
        self.engine.dispose()


def resolve_database_url(database_url: str | None = None) -> str:
    """Return an explicit or environment-configured local platform database URL."""
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
    return resolved_url


def open_run_persistence(database_url: str | None = None) -> RunPersistence:
    """Migrate and open the durable repository for one API process."""
    resolved_url = resolve_database_url(database_url)
    upgrade_database(resolved_url)
    engine = create_database_engine(resolved_url)
    repository = SqlAlchemyRunRepository(create_session_factory(engine))
    return RunPersistence(
        database_url=resolved_url,
        engine=engine,
        repository=repository,
    )
