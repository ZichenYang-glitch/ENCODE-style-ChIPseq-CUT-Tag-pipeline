"""Administrator-only local catalog and Reference Profile mutations."""

from __future__ import annotations

import argparse
from collections.abc import Callable, Iterable, Iterator, Mapping, Sequence
from contextlib import AbstractContextManager, contextmanager
import csv
from dataclasses import asdict, dataclass, is_dataclass
from datetime import datetime, timezone
from enum import Enum
import json
import os
from pathlib import Path
import sys
from typing import TYPE_CHECKING, Protocol, cast

if TYPE_CHECKING:
    from encode_pipeline.services.authentication_service import (
        AccountAdministrationService,
    )


@dataclass(frozen=True)
class SampleImportRow:
    """One validated, workflow-neutral row from an administrator TSV intake."""

    sample_key: str
    display_name: str
    attributes: dict[str, object]


class RegistryAdmin(Protocol):
    """Narrow command surface kept independent from SQL persistence details."""

    def create_project(self, *, display_name: str) -> object: ...

    def list_projects(self, *, include_archived: bool) -> Sequence[object]: ...

    def archive_project(self, *, project_id: str) -> object: ...

    def import_samples(
        self,
        *,
        project_id: str,
        rows: tuple[SampleImportRow, ...],
    ) -> Sequence[object]: ...

    def list_samples(self, *, project_id: str) -> Sequence[object]: ...

    def revise_sample(
        self,
        *,
        sample_id: str,
        display_name: str,
        attributes: dict[str, object],
    ) -> object: ...


class InputRegistryAdmin(Protocol):
    """Narrow local-only surface for StoragePool and InputFile mutations."""

    def register_storage_pool(
        self,
        *,
        display_name: str,
        config_key: str,
    ) -> object: ...

    def bind_project_storage_pool(
        self,
        *,
        project_id: str,
        storage_pool_id: str,
    ) -> object: ...

    def register_input_file(
        self,
        *,
        project_id: str,
        stable_key: str,
        pool_relative_path: str,
    ) -> object: ...


class ReferenceProfileAdmin(Protocol):
    """Narrow local-only surface for operator-prepared reference mutations."""

    def register(
        self,
        *,
        safe_key: str,
        display_name: str,
        organism: str,
        assembly: str,
        config_key: str,
    ) -> object: ...

    def verify(self, revision_id: str) -> object: ...

    def list(self) -> Sequence[object]: ...

    def enable(
        self,
        profile_id: str,
        *,
        revision_id: str | None = None,
    ) -> object: ...

    def disable(self, profile_id: str) -> object: ...


class RunRecoveryAdmin(Protocol):
    """Narrow administrator-only surface for explicit run recovery."""

    def diagnose(self, run_id: str) -> object: ...

    def summarize(self) -> object: ...

    def fail_run(
        self,
        run_id: str,
        *,
        expected_status: object,
        expected_assignment: object,
    ) -> object: ...

    def requeue_run(
        self,
        run_id: str,
        *,
        expected_status: object,
        expected_assignment: object,
    ) -> object: ...


class _ReferenceProfileCliError(RuntimeError):
    """Stable, path-free error for CLI-only admission failures."""

    def __init__(self, reason_code: str, message: str) -> None:
        self.reason_code = reason_code
        super().__init__(message)


class _RunRecoveryCliError(RuntimeError):
    """Fixed, disclosure-safe administrator recovery failure."""

    def __init__(self, reason_code: str) -> None:
        self.reason_code = reason_code
        super().__init__(_RUN_RECOVERY_ERROR_MESSAGES[reason_code])


class _RunRecoveryDatabaseUnavailable(RuntimeError):
    """Internal marker for a read-only database admission failure."""

    code = "RUN_RECOVERY_DATA_INVALID"


class _DataRegistryServiceLike(Protocol):
    def create_project(self, display_name: str) -> object: ...

    def list_projects(self, *, include_archived: bool) -> Sequence[object]: ...

    def archive_project(self, project_id: str) -> object: ...

    def import_samples(
        self,
        project_id: str,
        *,
        rows: Iterable[object],
    ) -> Sequence[object]: ...

    def list_samples(self, project_id: str) -> Sequence[object]: ...

    def revise_sample(
        self,
        sample_id: str,
        *,
        display_name: str,
        attributes: Mapping[str, object],
    ) -> object: ...


RegistryFactory = Callable[[str], AbstractContextManager[RegistryAdmin]]
InputRegistryFactory = Callable[
    [str, Path | None],
    AbstractContextManager[InputRegistryAdmin],
]
ReferenceProfileFactory = Callable[
    [str, Path | None, bool],
    AbstractContextManager[ReferenceProfileAdmin],
]
RunRecoveryFactory = Callable[..., AbstractContextManager[RunRecoveryAdmin]]


_RUN_STATUS_VALUES = (
    "created",
    "validating",
    "planned",
    "queued",
    "running",
    "succeeded",
    "failed",
    "cancelled",
)
_RUN_RECOVERY_ERROR_MESSAGES = {
    "RUN_RECOVERY_NOT_FOUND": "Run was not found.",
    "RUN_RECOVERY_CONFLICT": "Run recovery preconditions did not match.",
    "RUN_RECOVERY_NOT_SAFE": "Requested run recovery action is not safe.",
    "RUN_RECOVERY_QUEUE_UNAVAILABLE": "Run recovery queue is unavailable.",
    "RUN_RECOVERY_CLEANUP_FAILED": "Run recovery cleanup failed.",
    "RUN_RECOVERY_DATA_INVALID": "Run recovery data is invalid.",
    "RUN_RECOVERY_INTERNAL_ERROR": "Run recovery command failed.",
}


class _ServiceRegistryAdmin:
    """Translate stable CLI commands to the workflow-neutral registry service."""

    def __init__(
        self,
        service: object,
        *,
        sample_import_spec_factory: Callable[..., object],
    ) -> None:
        self._service = cast(_DataRegistryServiceLike, service)
        self._sample_import_spec_factory = sample_import_spec_factory

    def create_project(self, *, display_name: str) -> object:
        return self._service.create_project(display_name=display_name)

    def list_projects(self, *, include_archived: bool) -> Sequence[object]:
        return self._service.list_projects(include_archived=include_archived)

    def archive_project(self, *, project_id: str) -> object:
        return self._service.archive_project(project_id)

    def import_samples(
        self,
        *,
        project_id: str,
        rows: tuple[SampleImportRow, ...],
    ) -> Sequence[object]:
        specs = tuple(
            self._sample_import_spec_factory(
                stable_key=row.sample_key,
                display_name=row.display_name,
                attributes=row.attributes,
            )
            for row in rows
        )
        return self._service.import_samples(
            project_id,
            rows=specs,
        )

    def list_samples(self, *, project_id: str) -> Sequence[object]:
        return self._service.list_samples(project_id)

    def revise_sample(
        self,
        *,
        sample_id: str,
        display_name: str,
        attributes: dict[str, object],
    ) -> object:
        return self._service.revise_sample(
            sample_id,
            display_name=display_name,
            attributes=attributes,
        )


@contextmanager
def _open_registry(database_url: str) -> Iterator[RegistryAdmin]:
    # These imports stay local so help and the existing local-platform entry point
    # do not initialize SQLAlchemy or depend on administrator-only composition.
    from encode_pipeline.persistence.data_registry import (
        SqlAlchemyDataRegistryRepository,
    )
    from encode_pipeline.persistence.database import (
        create_database_engine,
        create_session_factory,
    )
    from encode_pipeline.persistence.migrations import upgrade_database
    from encode_pipeline.persistence.runtime import resolve_database_url
    from encode_pipeline.services.data_registry import (
        DataRegistryService,
        SampleImportSpec,
    )

    resolved_url = resolve_database_url(database_url)
    upgrade_database(resolved_url)
    engine = create_database_engine(resolved_url)
    try:
        repository = SqlAlchemyDataRegistryRepository(create_session_factory(engine))
        yield _ServiceRegistryAdmin(
            DataRegistryService(repository=repository),
            sample_import_spec_factory=SampleImportSpec,
        )
    finally:
        engine.dispose()


@contextmanager
def _open_account_administration(
    database_url: str,
) -> Iterator["AccountAdministrationService"]:
    """Compose the local account administration service for operator commands."""
    from encode_pipeline.persistence.authentication import (
        SqlAlchemyAuthenticationRepository,
    )
    from encode_pipeline.persistence.database import (
        create_database_engine,
        create_session_factory,
    )
    from encode_pipeline.persistence.migrations import upgrade_database
    from encode_pipeline.persistence.runtime import resolve_database_url
    from encode_pipeline.services.authentication_service import (
        AccountAdministrationService,
    )

    resolved_url = resolve_database_url(database_url)
    upgrade_database(resolved_url)
    engine = create_database_engine(resolved_url)
    try:
        repository = SqlAlchemyAuthenticationRepository(create_session_factory(engine))
        yield AccountAdministrationService(repository=repository)
    finally:
        engine.dispose()


def _read_interactive_password(
    reader: Callable[[str], str],
    *,
    prompt: str = "Password: ",
    confirmation_prompt: str = "Confirm password: ",
) -> str:
    """Read a password interactively without echo, argv, or environment leaks."""
    first = reader(prompt)
    if not isinstance(first, str) or not first:
        raise ValueError("password must not be empty")
    second = reader(confirmation_prompt)
    if first != second:
        raise ValueError("passwords do not match")
    return first


@contextmanager
def _open_input_registry(
    database_url: str,
    storage_pool_config: Path | None,
) -> Iterator[InputRegistryAdmin]:
    """Compose the private input registry lazily for administrator commands."""
    # The concrete composition is completed with the Stage 3 registry service.
    # Keeping the import local preserves the lightweight help boundary.
    from encode_pipeline.services.input_registry import open_input_registry_admin

    with open_input_registry_admin(
        database_url=database_url,
        storage_pool_config=storage_pool_config,
    ) as registry:
        yield cast(InputRegistryAdmin, registry)


@contextmanager
def _open_reference_profiles(
    database_url: str,
    reference_profile_config: Path | None,
    allow_schema_upgrade: bool,
) -> Iterator[ReferenceProfileAdmin]:
    """Compose the private Reference Profile service only for admin commands."""
    from encode_pipeline.persistence.database import (
        create_database_engine,
        create_session_factory,
    )
    from encode_pipeline.persistence.migrations import upgrade_database
    from encode_pipeline.persistence.reference_profiles import (
        SqlAlchemyReferenceProfileRepository,
    )
    from encode_pipeline.persistence.migration_admission import (
        verify_migration_execution_inventory,
    )
    from encode_pipeline.persistence.runtime import resolve_database_url
    from encode_pipeline.services.defaults import create_default_workflow_registry
    from encode_pipeline.services.private_reference_profiles import (
        PrivateReferenceProfileConfigError,
        load_private_reference_profile_config,
    )
    from encode_pipeline.services.reference_profiles import ReferenceProfileService

    resolved_url = resolve_database_url(database_url)
    if allow_schema_upgrade:
        upgrade_database(resolved_url)
        engine = create_database_engine(resolved_url)
    else:
        import sqlite3

        from sqlalchemy import create_engine
        from sqlalchemy.engine import make_url

        inventory = verify_migration_execution_inventory()
        if len(inventory.heads) != 1:
            raise RuntimeError("Reference Profile database schema is unavailable.")
        database_path = Path(cast(str, make_url(resolved_url).database))

        def open_read_only_database():
            connection = sqlite3.connect(
                f"{database_path.as_uri()}?mode=ro",
                uri=True,
                check_same_thread=False,
                timeout=30,
            )
            connection.execute("PRAGMA query_only=ON")
            connection.execute("PRAGMA foreign_keys=ON")
            return connection

        engine = create_engine(
            "sqlite://",
            creator=open_read_only_database,
            future=True,
        )
        try:
            with engine.connect() as connection:
                observed_heads = tuple(
                    sorted(
                        connection.exec_driver_sql(
                            "SELECT version_num FROM alembic_version"
                        ).scalars()
                    )
                )
            if observed_heads != inventory.heads:
                raise RuntimeError("Reference Profile database schema is unavailable.")
        except Exception:
            engine.dispose()
            raise
    try:
        repository = SqlAlchemyReferenceProfileRepository(
            create_session_factory(engine)
        )
        registry = create_default_workflow_registry()

        def private_config_provider():
            if reference_profile_config is None:
                raise PrivateReferenceProfileConfigError()
            return load_private_reference_profile_config(reference_profile_config)

        yield ReferenceProfileService(
            repository=repository,
            private_config_provider=private_config_provider,
            adapter_provider=registry.get,
        )
    finally:
        engine.dispose()


@contextmanager
def _open_run_recovery(
    database_url: str,
    *,
    read_only: bool,
    redis_url: str | None = None,
    queue_name: str | None = None,
) -> Iterator[RunRecoveryAdmin]:
    """Compose recovery lazily without letting diagnosis migrate or create state."""
    import sqlite3

    from sqlalchemy import create_engine
    from sqlalchemy.engine import make_url

    from encode_pipeline.persistence.database import (
        create_database_engine,
        create_session_factory,
    )
    from encode_pipeline.persistence.migration_admission import (
        verify_migration_execution_inventory,
    )
    from encode_pipeline.persistence.migrations import upgrade_database
    from encode_pipeline.persistence.repositories import SqlAlchemyRunRepository
    from encode_pipeline.persistence.runtime import (
        DATABASE_URL_ENV,
        resolve_database_url,
    )
    from encode_pipeline.services.managed_containers import ManagedContainerCleaner
    from encode_pipeline.services.run_recovery import RunRecoveryService
    from encode_pipeline.workers.rq_queue import RqRunQueue
    from encode_pipeline.workers.settings import (
        QUEUE_NAME_ENV,
        REDIS_URL_ENV,
        load_worker_settings,
    )

    resolved_url = resolve_database_url(database_url)
    if read_only:
        try:
            inventory = verify_migration_execution_inventory()
            if len(inventory.heads) != 1:
                raise ValueError("migration inventory has multiple heads")
            database_path = Path(cast(str, make_url(resolved_url).database))
            if not database_path.is_file():
                raise FileNotFoundError("database is unavailable")
        except Exception:
            raise _RunRecoveryDatabaseUnavailable(
                "Run recovery database is unavailable."
            ) from None

        def open_read_only_database():
            connection = sqlite3.connect(
                f"{database_path.as_uri()}?mode=ro",
                uri=True,
                check_same_thread=False,
                timeout=30,
            )
            connection.execute("PRAGMA query_only=ON")
            connection.execute("PRAGMA foreign_keys=ON")
            return connection

        engine = create_engine(
            "sqlite://",
            creator=open_read_only_database,
            future=True,
        )
        try:
            with engine.connect() as connection:
                observed_heads = tuple(
                    sorted(
                        connection.exec_driver_sql(
                            "SELECT version_num FROM alembic_version"
                        ).scalars()
                    )
                )
        except Exception:
            engine.dispose()
            raise _RunRecoveryDatabaseUnavailable(
                "Run recovery database is unavailable."
            ) from None
        if observed_heads != inventory.heads:
            engine.dispose()
            raise _RunRecoveryDatabaseUnavailable(
                "Run recovery database is unavailable."
            )
    else:
        upgrade_database(resolved_url)
        engine = create_database_engine(resolved_url)

    settings_environment = dict(os.environ)
    settings_environment[DATABASE_URL_ENV] = resolved_url
    if redis_url is not None:
        settings_environment[REDIS_URL_ENV] = redis_url
    if queue_name is not None:
        settings_environment[QUEUE_NAME_ENV] = queue_name
    queue = None
    try:
        repository = SqlAlchemyRunRepository(create_session_factory(engine))
        settings = load_worker_settings(settings_environment)
        queue = RqRunQueue(settings)
        cleanup = None
        cleanup_endpoint_identity = None
        cleaner = None
        if settings.managed_docker_executable is not None:
            try:
                cleaner = ManagedContainerCleaner(
                    executable=settings.managed_docker_executable,
                    unix_socket=settings.managed_docker_socket,
                )
            except (OSError, ValueError):
                cleaner = None
        if cleaner is not None:
            cleanup_endpoint_identity = cleaner.endpoint_identity

            def cleanup(scope: str) -> bool:
                try:
                    return cleaner.cleanup(scope).is_success
                except (OSError, RuntimeError, ValueError):
                    return False

        yield RunRecoveryService(
            repository,
            queue,
            cleanup=cleanup,
            cleanup_endpoint_identity=cleanup_endpoint_identity,
        )
    finally:
        if queue is not None:
            try:
                queue.close()
            except Exception:
                pass
        engine.dispose()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="helixweave admin",
        description=(
            "Manage the local HelixWeave data catalog and run recovery. Catalog "
            "mutations target the explicitly selected SQLite database; recovery "
            "coordinates its configured local queue and cleanup boundary."
        ),
    )
    parser.add_argument(
        "--database-url",
        required=True,
        help="Absolute file-backed SQLite URL selected by the local administrator.",
    )
    parser.add_argument(
        "--storage-pool-config",
        type=Path,
        help=(
            "Private operator StoragePool mapping. Required only for commands "
            "that inspect local file bytes."
        ),
    )
    parser.add_argument(
        "--reference-profile-config",
        type=Path,
        help=(
            "Private operator Reference Profile bindings. Required only for "
            "register, verify, and enable."
        ),
    )
    resource_parsers = parser.add_subparsers(dest="resource", required=True)

    project = resource_parsers.add_parser("project", help="Manage Projects.")
    project_commands = project.add_subparsers(
        dest="project_command",
        required=True,
    )
    project_create = project_commands.add_parser("create", help="Create a Project.")
    project_create.add_argument("--display-name", required=True)
    project_list = project_commands.add_parser("list", help="List Projects.")
    project_list.add_argument("--include-archived", action="store_true")
    project_archive = project_commands.add_parser(
        "archive",
        help="Archive a Project without deleting history.",
    )
    project_archive.add_argument("project_id")
    project_bind_pool = project_commands.add_parser(
        "bind-storage-pool",
        help="Bind one Project to one approved StoragePool.",
    )
    project_bind_pool.add_argument("--project-id", required=True)
    project_bind_pool.add_argument("--storage-pool-id", required=True)

    sample = resource_parsers.add_parser("sample", help="Manage Samples.")
    sample_commands = sample.add_subparsers(dest="sample_command", required=True)
    sample_import = sample_commands.add_parser(
        "import",
        help="Import Samples from a reviewed local TSV.",
    )
    sample_import.add_argument("--project-id", required=True)
    sample_import.add_argument("--tsv", required=True, type=Path)
    sample_list = sample_commands.add_parser("list", help="List Project Samples.")
    sample_list.add_argument("--project-id", required=True)
    sample_revise = sample_commands.add_parser(
        "revise",
        help="Append an immutable Sample revision.",
    )
    sample_revise.add_argument("sample_id")
    sample_revise.add_argument("--display-name", required=True)
    sample_revise.add_argument("--attributes-json", default="{}")

    storage_pool = resource_parsers.add_parser(
        "storage-pool",
        help="Manage approved local StoragePools.",
    )
    storage_pool_commands = storage_pool.add_subparsers(
        dest="storage_pool_command",
        required=True,
    )
    storage_pool_register = storage_pool_commands.add_parser(
        "register",
        help="Register an approved private configuration key.",
    )
    storage_pool_register.add_argument("--display-name", required=True)
    storage_pool_register.add_argument("--config-key", required=True)

    input_file = resource_parsers.add_parser(
        "input-file",
        help="Register immutable revisions of local regular files.",
    )
    input_file_commands = input_file.add_subparsers(
        dest="input_file_command",
        required=True,
    )
    input_file_register = input_file_commands.add_parser(
        "register",
        help="Register or revise one pool-relative regular file.",
    )
    input_file_register.add_argument("--project-id", required=True)
    input_file_register.add_argument("--stable-key", required=True)
    input_file_register.add_argument("--pool-relative-path", required=True)

    reference_profile = resource_parsers.add_parser(
        "reference-profile",
        help="Manage operator-prepared immutable Reference Profiles.",
    )
    reference_commands = reference_profile.add_subparsers(
        dest="reference_profile_command",
        required=True,
    )
    reference_register = reference_commands.add_parser(
        "register",
        help="Create a stable profile or append a verified revision.",
    )
    reference_register.add_argument("--safe-key", required=True)
    reference_register.add_argument("--display-name", required=True)
    reference_register.add_argument("--organism", required=True)
    reference_register.add_argument("--assembly", required=True)
    reference_register.add_argument("--config-key", required=True)
    reference_verify = reference_commands.add_parser(
        "verify",
        help="Verify one revision without mutating the database.",
    )
    reference_verify.add_argument("revision_id")
    reference_commands.add_parser("list", help="List path-free profile metadata.")
    reference_enable = reference_commands.add_parser(
        "enable",
        help="Verify and enable one exact revision.",
    )
    reference_enable.add_argument("profile_id")
    reference_enable.add_argument("--revision-id")
    reference_disable = reference_commands.add_parser(
        "disable",
        help="Disable new use without deleting history.",
    )
    reference_disable.add_argument("profile_id")

    account = resource_parsers.add_parser(
        "account",
        help="Bootstrap and recover local user accounts.",
    )
    account_commands = account.add_subparsers(dest="account_command", required=True)
    account_bootstrap = account_commands.add_parser(
        "bootstrap",
        help="Create the unique first administrator (interactive password).",
    )
    account_bootstrap.add_argument("--username", required=True)
    account_reset = account_commands.add_parser(
        "reset-password",
        help="Reset one account password interactively and revoke its sessions.",
    )
    account_reset.add_argument("--username", required=True)
    account_commands.add_parser(
        "list",
        help="List safe account summaries without secret material.",
    )

    run = resource_parsers.add_parser(
        "run",
        help="Diagnose and explicitly recover durable Runs.",
    )
    run_commands = run.add_subparsers(dest="run_command", required=True)
    run_diagnose = run_commands.add_parser(
        "diagnose",
        help="Inspect one run and its exact queue identity without mutation.",
    )
    run_diagnose.add_argument("run_id")
    run_diagnose.add_argument(
        "--queue-name",
        help="Select the configured local queue without exposing its endpoint.",
    )
    for command, help_text in (
        ("fail", "Fail one recoverable run after exact precondition checks."),
        ("requeue", "Requeue one recoverable run exactly once."),
    ):
        mutation = run_commands.add_parser(command, help=help_text)
        mutation.add_argument("run_id")
        mutation.add_argument(
            "--expected-status",
            required=True,
            choices=_RUN_STATUS_VALUES,
        )
        mutation.add_argument("--job-id", required=True)
        mutation.add_argument("--backend", required=True)
        mutation.add_argument("--queue-name", required=True)
    return parser


def _parse_json_object(raw_value: str, *, context: str) -> dict[str, object]:
    try:
        parsed = json.loads(raw_value)
    except (json.JSONDecodeError, TypeError) as exc:
        raise ValueError(f"{context} must contain a valid JSON object") from exc
    if not isinstance(parsed, dict):
        raise ValueError(f"{context} must contain a JSON object")
    return parsed


def _read_sample_rows(path: Path) -> tuple[SampleImportRow, ...]:
    try:
        handle = path.open("r", encoding="utf-8-sig", newline="")
    except OSError as exc:
        raise ValueError("sample TSV could not be read") from exc
    with handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = tuple(reader.fieldnames or ())
        required = {"sample_key", "display_name"}
        allowed = required | {"attributes_json"}
        if any(not fieldname.strip() for fieldname in fieldnames):
            raise ValueError("sample TSV has an empty column name")
        duplicates = sorted(
            {fieldname for fieldname in fieldnames if fieldnames.count(fieldname) > 1}
        )
        if duplicates:
            raise ValueError(
                "sample TSV has duplicate column(s): " + ", ".join(duplicates)
            )
        unknown = sorted(set(fieldnames) - allowed)
        if unknown:
            raise ValueError("sample TSV has unknown column(s): " + ", ".join(unknown))
        missing = sorted(required - set(fieldnames))
        if missing:
            raise ValueError(
                "sample TSV is missing required column(s): " + ", ".join(missing)
            )
        rows: list[SampleImportRow] = []
        seen_keys: set[str] = set()
        for row_number, row in enumerate(reader, start=2):
            if any(key is None for key in row):
                raise ValueError(
                    f"sample TSV row {row_number} has unexpected trailing field(s)"
                )
            sample_key = (row.get("sample_key") or "").strip()
            display_name = (row.get("display_name") or "").strip()
            if not sample_key:
                raise ValueError(f"sample TSV row {row_number} has an empty sample_key")
            if not display_name:
                raise ValueError(
                    f"sample TSV row {row_number} has an empty display_name"
                )
            if sample_key in seen_keys:
                raise ValueError(
                    f"sample TSV row {row_number} repeats sample_key {sample_key!r}"
                )
            seen_keys.add(sample_key)
            raw_attributes = (row.get("attributes_json") or "").strip() or "{}"
            attributes = _parse_json_object(
                raw_attributes,
                context=f"sample TSV row {row_number} attributes_json",
            )
            rows.append(
                SampleImportRow(
                    sample_key=sample_key,
                    display_name=display_name,
                    attributes=attributes,
                )
            )
    if not rows:
        raise ValueError("sample TSV must contain at least one data row")
    return tuple(rows)


def _jsonable(value: object) -> object:
    if is_dataclass(value) and not isinstance(value, type):
        return _jsonable(asdict(value))
    if isinstance(value, Enum):
        return _jsonable(value.value)
    if isinstance(value, datetime):
        normalized = value.astimezone(timezone.utc)
        return normalized.isoformat(timespec="microseconds").replace("+00:00", "Z")
    if isinstance(value, Mapping):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise TypeError(f"administrator result contains unsupported {type(value).__name__}")


def _write_result(result: object) -> None:
    json.dump(
        _jsonable(result),
        sys.stdout,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    sys.stdout.write("\n")


def _write_command_error(error: Exception) -> None:
    code = getattr(error, "reason_code", "registry-command-failed")
    if not isinstance(code, str) or not code or len(code) > 128:
        code = "registry-command-failed"
    json.dump(
        {
            "error": {
                "code": code,
                "message": str(error),
            }
        },
        sys.stderr,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    sys.stderr.write("\n")


def _value(value: object) -> object:
    """Return one explicit scalar enum value for a fixed CLI projection."""
    if isinstance(value, Enum):
        return value.value
    return value


def _state_value(value: object) -> object:
    state = getattr(value, "state", value)
    return _value(state)


def _run_recovery_diagnostic_payload(diagnostic: object) -> dict[str, object]:
    """Project only the reviewed public recovery fields."""
    assignment = getattr(diagnostic, "assignment")
    assignment_payload: dict[str, object] | None = None
    if assignment is not None:
        assignment_payload = {
            "job_id": getattr(assignment, "job_id"),
            "backend": getattr(assignment, "backend"),
            "queue_name": getattr(assignment, "queue_name"),
            "dispatched": getattr(assignment, "dispatched_at") is not None,
            "claimed": getattr(assignment, "claimed_at") is not None,
            "cancellation_requested": (
                getattr(assignment, "cancellation_requested_at") is not None
            ),
            "cancellation_acknowledged": (
                getattr(assignment, "cancellation_acknowledged_at") is not None
            ),
            "requeue_requested": (
                getattr(assignment, "requeue_requested_at", None) is not None
            ),
            "requeue_confirmed": (
                getattr(assignment, "requeue_confirmed_at", None) is not None
            ),
        }
    return {
        "run_id": getattr(diagnostic, "run_id"),
        "workflow_id": getattr(diagnostic, "workflow_id"),
        "status": _value(getattr(diagnostic, "status")),
        "diagnosis_code": getattr(diagnostic, "diagnosis_code"),
        "assignment": assignment_payload,
        "queue": {
            "state": _state_value(getattr(diagnostic, "queue_evidence")),
        },
        "result_indexing": _state_value(getattr(diagnostic, "result_indexing")),
        "cleanup": _state_value(getattr(diagnostic, "cleanup")),
        "gaps": sorted(_value(item) for item in getattr(diagnostic, "gaps")),
        "allowed_actions": sorted(
            _value(item) for item in getattr(diagnostic, "allowed_actions")
        ),
    }


def _run_recovery_action_payload(
    result: object,
    *,
    action: str,
    previous_status: object,
) -> dict[str, object]:
    payload: dict[str, object] = {
        "action": action,
        "run_id": getattr(result, "run_id"),
        "status": _value(getattr(result, "status")),
        "reason_code": getattr(result, "reason_code"),
        "changed": getattr(result, "changed"),
    }
    if action == "fail":
        payload["previous_status"] = _value(
            getattr(result, "previous_status", previous_status)
        )
        return payload
    assignment = getattr(result, "assignment")
    payload.update(
        {
            "job_id": getattr(assignment, "job_id"),
            "requeue_requested": (
                getattr(assignment, "requeue_requested_at", None) is not None
            ),
            "requeue_confirmed": (
                getattr(assignment, "requeue_confirmed_at", None) is not None
            ),
        }
    )
    return payload


def _recovery_preconditions_match(
    diagnostic: object,
    *,
    expected_status: object,
    job_id: str,
    backend: str,
    queue_name: str,
) -> bool:
    assignment = getattr(diagnostic, "assignment", None)
    return bool(
        _value(getattr(diagnostic, "status", None)) == _value(expected_status)
        and assignment is not None
        and getattr(assignment, "job_id", None) == job_id
        and getattr(assignment, "backend", None) == backend
        and getattr(assignment, "queue_name", None) == queue_name
    )


def _write_run_recovery_error(error: Exception) -> None:
    candidate = getattr(error, "code", None)
    if candidate not in _RUN_RECOVERY_ERROR_MESSAGES:
        candidate = getattr(error, "reason_code", None)
    code = (
        candidate
        if candidate in _RUN_RECOVERY_ERROR_MESSAGES
        else "RUN_RECOVERY_INTERNAL_ERROR"
    )
    _write_command_error(_RunRecoveryCliError(cast(str, code)))


def main(
    argv: Sequence[str] | None = None,
    *,
    registry_factory: RegistryFactory = _open_registry,
    input_registry_factory: InputRegistryFactory = _open_input_registry,
    reference_profile_factory: ReferenceProfileFactory = _open_reference_profiles,
    run_recovery_factory: RunRecoveryFactory | None = None,
    account_administration_factory: (
        "Callable[[str], AbstractContextManager[object]] | None"
    ) = None,
    password_reader: Callable[[str], str] | None = None,
) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    rows: tuple[SampleImportRow, ...] | None = None
    attributes: dict[str, object] | None = None
    try:
        if args.resource == "sample" and args.sample_command == "import":
            rows = _read_sample_rows(args.tsv)
        elif args.resource == "sample" and args.sample_command == "revise":
            attributes = _parse_json_object(
                args.attributes_json,
                context="--attributes-json",
            )
    except ValueError as exc:
        parser.error(str(exc))

    try:
        if args.resource == "account":
            if account_administration_factory is None:
                account_administration_factory = _open_account_administration
            reader = password_reader
            if reader is None:
                import getpass

                reader = getpass.getpass
            with account_administration_factory(args.database_url) as administration:
                if args.account_command == "bootstrap":
                    password = _read_interactive_password(reader)
                    account = administration.bootstrap_initial_administrator(
                        args.username,
                        password,
                    )
                    result = account.to_public_summary()
                elif args.account_command == "reset-password":
                    password = _read_interactive_password(reader)
                    account = administration.reset_password_for_username(
                        args.username,
                        password,
                    )
                    result = account.to_public_summary()
                else:
                    result = [
                        account.to_public_summary()
                        for account in administration.list_accounts()
                    ]
        elif args.resource == "run":
            if run_recovery_factory is None:
                run_recovery_factory = _open_run_recovery
            read_only = args.run_command == "diagnose"
            recovery_factory_options: dict[str, object] = {"read_only": read_only}
            if args.queue_name is not None:
                recovery_factory_options["queue_name"] = args.queue_name
            with run_recovery_factory(
                args.database_url,
                **recovery_factory_options,
            ) as recovery:
                diagnostic = recovery.diagnose(args.run_id)
                if args.run_command == "diagnose":
                    result = _run_recovery_diagnostic_payload(diagnostic)
                else:
                    from encode_pipeline.platform.runs import RunStatus

                    expected_status = RunStatus(args.expected_status)
                    if not _recovery_preconditions_match(
                        diagnostic,
                        expected_status=expected_status,
                        job_id=args.job_id,
                        backend=args.backend,
                        queue_name=args.queue_name,
                    ):
                        raise _RunRecoveryCliError("RUN_RECOVERY_CONFLICT")
                    assignment = getattr(diagnostic, "assignment")
                    if args.run_command == "fail":
                        action_result = recovery.fail_run(
                            args.run_id,
                            expected_status=expected_status,
                            expected_assignment=assignment,
                        )
                    else:
                        action_result = recovery.requeue_run(
                            args.run_id,
                            expected_status=expected_status,
                            expected_assignment=assignment,
                        )
                    result = _run_recovery_action_payload(
                        action_result,
                        action=args.run_command,
                        previous_status=expected_status,
                    )
        elif args.resource == "reference-profile":
            if (
                args.reference_profile_command in {"register", "verify", "enable"}
                and args.reference_profile_config is None
            ):
                raise _ReferenceProfileCliError(
                    "REFERENCE_PROFILE_CONFIG_REQUIRED",
                    "Private Reference Profile configuration is required.",
                )
            with reference_profile_factory(
                args.database_url,
                args.reference_profile_config,
                args.reference_profile_command in {"register", "enable", "disable"},
            ) as profiles:
                if args.reference_profile_command == "register":
                    result = profiles.register(
                        safe_key=args.safe_key,
                        display_name=args.display_name,
                        organism=args.organism,
                        assembly=args.assembly,
                        config_key=args.config_key,
                    )
                elif args.reference_profile_command == "verify":
                    result = profiles.verify(args.revision_id)
                elif args.reference_profile_command == "list":
                    result = profiles.list()
                elif args.reference_profile_command == "enable":
                    result = profiles.enable(
                        args.profile_id,
                        revision_id=args.revision_id,
                    )
                else:
                    result = profiles.disable(args.profile_id)
        elif args.resource in {"storage-pool", "input-file"} or (
            args.resource == "project" and args.project_command == "bind-storage-pool"
        ):
            if (
                args.resource == "storage-pool" or args.resource == "input-file"
            ) and args.storage_pool_config is None:
                parser.error(
                    "--storage-pool-config is required for local storage commands"
                )
            with input_registry_factory(
                args.database_url,
                args.storage_pool_config,
            ) as input_registry:
                if args.resource == "storage-pool":
                    result = input_registry.register_storage_pool(
                        display_name=args.display_name,
                        config_key=args.config_key,
                    )
                elif args.resource == "project":
                    result = input_registry.bind_project_storage_pool(
                        project_id=args.project_id,
                        storage_pool_id=args.storage_pool_id,
                    )
                else:
                    result = input_registry.register_input_file(
                        project_id=args.project_id,
                        stable_key=args.stable_key,
                        pool_relative_path=args.pool_relative_path,
                    )
        else:
            with registry_factory(args.database_url) as registry:
                if args.resource == "project":
                    if args.project_command == "create":
                        result = registry.create_project(display_name=args.display_name)
                    elif args.project_command == "list":
                        result = registry.list_projects(
                            include_archived=args.include_archived
                        )
                    else:
                        result = registry.archive_project(project_id=args.project_id)
                elif args.sample_command == "import":
                    assert rows is not None
                    result = registry.import_samples(
                        project_id=args.project_id,
                        rows=rows,
                    )
                elif args.sample_command == "list":
                    result = registry.list_samples(project_id=args.project_id)
                else:
                    assert attributes is not None
                    result = registry.revise_sample(
                        sample_id=args.sample_id,
                        display_name=args.display_name,
                        attributes=attributes,
                    )
    except Exception as exc:
        if args.resource == "run":
            _write_run_recovery_error(exc)
            return 1
        if args.resource == "reference-profile":
            if not isinstance(getattr(exc, "reason_code", None), str):
                exc = _ReferenceProfileCliError(
                    "REFERENCE_PROFILE_COMMAND_FAILED",
                    "Reference Profile command failed.",
                )
        elif not isinstance(exc, (LookupError, OSError, RuntimeError, ValueError)):
            raise
        _write_command_error(exc)
        return 1
    _write_result(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
