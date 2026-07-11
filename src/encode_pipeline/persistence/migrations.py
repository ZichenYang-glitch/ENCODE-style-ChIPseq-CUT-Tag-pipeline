"""Programmatic Alembic entry points for platform database migrations."""

from __future__ import annotations

from pathlib import Path

from alembic import command
from alembic.config import Config


def alembic_config(database_url: str) -> Config:
    """Build an Alembic config without relying on a process working directory."""
    if not isinstance(database_url, str) or not database_url.strip():
        raise ValueError("database_url must be a non-empty string")
    config = Config()
    config.set_main_option(
        "script_location",
        str(Path(__file__).with_name("alembic")),
    )
    config.set_main_option("sqlalchemy.url", database_url.strip().replace("%", "%%"))
    return config


def upgrade_database(database_url: str, revision: str = "head") -> None:
    """Upgrade a platform database to the requested revision."""
    command.upgrade(alembic_config(database_url), revision)


def downgrade_database(database_url: str, revision: str) -> None:
    """Downgrade a platform database to the requested revision."""
    command.downgrade(alembic_config(database_url), revision)
