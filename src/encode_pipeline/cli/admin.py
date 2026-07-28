"""Administrator-only local commands for Project and Sample catalog mutations."""

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
    json.dump(
        {
            "error": {
                "code": "registry-command-failed",
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
    except (LookupError, OSError, RuntimeError, ValueError) as exc:
        _write_command_error(exc)
        return 1
    _write_result(result)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
