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
from pathlib import Path
import sys
from typing import Protocol, cast


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


class _ReferenceProfileCliError(RuntimeError):
    """Stable, path-free error for CLI-only admission failures."""

    def __init__(self, reason_code: str, message: str) -> None:
        self.reason_code = reason_code
        super().__init__(message)


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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="helixweave admin",
        description=(
            "Manage the local HelixWeave data catalog. These commands mutate only "
            "the explicitly selected SQLite database."
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


def main(
    argv: Sequence[str] | None = None,
    *,
    registry_factory: RegistryFactory = _open_registry,
    input_registry_factory: InputRegistryFactory = _open_input_registry,
    reference_profile_factory: ReferenceProfileFactory = _open_reference_profiles,
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
        if args.resource == "reference-profile":
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
