"""Programmatic Alembic entry points for platform database migrations."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from encode_pipeline.persistence.migration_admission import (
    MigrationAdmissionError,
    VerifiedMigrationExecutionInventory,
    admitted_migration_script_location,
)


def _admitted_alembic_config(
    database_url: str,
    *,
    script_location: Path,
) -> Any:
    """Build an internal config for an already-admitted script snapshot."""
    if not isinstance(database_url, str) or not database_url.strip():
        raise ValueError("database_url must be a non-empty string")
    if not isinstance(script_location, Path) or not script_location.is_absolute():
        raise ValueError("script_location must be an absolute pathlib.Path")
    from alembic.config import Config

    config = Config()
    config.set_main_option(
        "script_location",
        str(script_location).replace("%", "%%"),
    )
    config.set_main_option("sqlalchemy.url", database_url.strip().replace("%", "%%"))
    return config


def upgrade_database(database_url: str, revision: str = "head") -> None:
    """Upgrade a platform database to the requested revision."""
    with admitted_migration_script_location() as (script_location, inventory):
        config = _admitted_alembic_config(
            database_url,
            script_location=script_location,
        )
        scripts = _validated_script_directory(config, inventory)
        _ensure_database_parent(database_url)
        _run_validated_upgrade(config, scripts, revision)


def downgrade_database(database_url: str, revision: str) -> None:
    """Downgrade a platform database to the requested revision."""
    with admitted_migration_script_location() as (script_location, inventory):
        config = _admitted_alembic_config(
            database_url,
            script_location=script_location,
        )
        scripts = _validated_script_directory(config, inventory)
        _ensure_database_parent(database_url)
        _run_validated_downgrade(config, scripts, revision)


def _validated_script_directory(
    config: Any,
    inventory: VerifiedMigrationExecutionInventory,
) -> Any:
    """Return the exact Alembic graph validated before database path effects."""
    try:
        from alembic.script import ScriptDirectory

        scripts = ScriptDirectory.from_config(config)
        observed_heads = tuple(sorted(scripts.get_heads()))
        observed_bases = tuple(sorted(scripts.get_bases()))
        observed_revisions: dict[
            str,
            tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]],
        ] = {}
        for script in scripts.walk_revisions(base="base", head="heads"):
            module = script.module
            revision = module.revision
            if not isinstance(revision, str) or revision in observed_revisions:
                raise ValueError("ambiguous runtime revision")
            observed_revisions[revision] = (
                _runtime_revision_tuple(module.down_revision),
                _runtime_revision_tuple(module.branch_labels),
                _runtime_revision_tuple(module.depends_on),
            )
        expected_revisions = {
            item.revision: (
                item.down_revision,
                item.branch_labels,
                item.depends_on,
            )
            for item in inventory.revisions
        }
        if (
            observed_heads != inventory.heads
            or observed_bases != inventory.bases
            or observed_revisions != expected_revisions
        ):
            raise ValueError("runtime revision graph differs from inventory")
        return scripts
    except MigrationAdmissionError:
        raise
    except Exception:
        raise MigrationAdmissionError("MIGRATION_REVISION_GRAPH_INVALID") from None


def _run_validated_upgrade(config: Any, scripts: Any, revision: str) -> None:
    """Run an upgrade without constructing a second ScriptDirectory."""
    from alembic.runtime.environment import EnvironmentContext

    def upgrade_revisions(current_revision: str, _context: Any) -> Any:
        return scripts._upgrade_revs(revision, current_revision)

    with EnvironmentContext(
        config,
        scripts,
        fn=upgrade_revisions,
        destination_rev=revision,
    ):
        scripts.run_env()


def _run_validated_downgrade(config: Any, scripts: Any, revision: str) -> None:
    """Run a downgrade without constructing a second ScriptDirectory."""
    from alembic.runtime.environment import EnvironmentContext

    def downgrade_revisions(current_revision: str, _context: Any) -> Any:
        return scripts._downgrade_revs(revision, current_revision)

    with EnvironmentContext(
        config,
        scripts,
        fn=downgrade_revisions,
        destination_rev=revision,
    ):
        scripts.run_env()


def _runtime_revision_tuple(value: object) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        normalized = (value,)
    elif isinstance(value, (tuple, list)):
        normalized = tuple(value)
    else:
        raise ValueError("unsupported runtime revision metadata")
    if (
        any(not isinstance(item, str) or not item for item in normalized)
        or tuple(sorted(set(normalized))) != normalized
    ):
        raise ValueError("invalid runtime revision metadata")
    return normalized


def _ensure_database_parent(database_url: str) -> None:
    """Create a file-backed SQLite parent before Alembic opens its own engine."""
    from encode_pipeline.persistence.database import create_database_engine

    engine = create_database_engine(database_url)
    engine.dispose()
