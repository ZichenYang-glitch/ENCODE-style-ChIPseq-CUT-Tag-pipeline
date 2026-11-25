"""Operator commands for sealing and restoring exact historical tracks."""

from __future__ import annotations

import argparse
from collections.abc import Callable, Mapping, Sequence
import json
import os
from pathlib import Path
import sys

from encode_pipeline.persistence.runtime import DATABASE_URL_ENV
from encode_pipeline.services.archive_sealing import (
    ArchiveSealError,
    require_operator_receipt,
)
from encode_pipeline.workers.settings import (
    ARCHIVE_ROOT_ENV,
    STORAGE_POOL_CONFIG_ENV,
    WORKSPACE_ROOT_ENV,
    WorkerSettings,
    load_worker_settings,
)


_EXAMPLES = """
examples:
  helixweave archive seal registered \\
    --project-id prj_... --experiment-id exp_... \\
    --experiment-revision-id expr_... --track-key signal \\
    --display-name signal.bw --assembly hg38 --track-type bigwig \\
    --result-id rreg_... --result-revision-id rrev_...

  helixweave archive seal run \\
    --project-id prj_... --experiment-id exp_... \\
    --experiment-revision-id expr_... --track-key signal \\
    --display-name signal.bw --assembly hg38 --track-type bigwig \\
    --run-id RUN --artifact-id signal --generation artifactgen-... \\
    --revision artifactrev-...

Bytes are read from that revision or from the published run artifact.
When the registered path has been retired, its recorded digest selects and
verifies the corresponding object in the configured archive. For a relocated
run publication, supply its operator directory with --relocation-root.
Longrange tracks require both the data revision and its tabix index revision.
A filename, a relative path, and the latest revision are not arguments.
ENCODE_PIPELINE_ARCHIVE_ROOT, ENCODE_PIPELINE_DATABASE_URL, and
ENCODE_PIPELINE_STORAGE_POOL_CONFIG supply the operator locations.
""".strip()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="helixweave archive",
        description="Seal a track version or restore its verified archive copy.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_EXAMPLES,
    )
    commands = parser.add_subparsers(dest="command", required=True)
    scope = _scope_parser()
    restore = commands.add_parser(
        "restore",
        parents=[scope],
        help="Verify all members of a restored archive version and recover its original link.",
    )
    restore.add_argument("--archive-version-id", required=True)
    seal = commands.add_parser(
        "seal",
        help="Copy one exact source and publish its sealed version.",
        description=(
            "Copy one exact source and publish its sealed version. "
            "This command does not accept a filename."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_EXAMPLES,
    )
    kinds = seal.add_subparsers(dest="kind", required=True)
    common = _common_parser(scope)
    registered = kinds.add_parser(
        "registered",
        parents=[common],
        help="Seal one registered result revision and an optional index revision.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_EXAMPLES,
    )
    registered.add_argument("--result-id", required=True)
    registered.add_argument("--result-revision-id", required=True)
    registered.add_argument("--index-result-id")
    registered.add_argument("--index-result-revision-id")
    run = kinds.add_parser(
        "run",
        parents=[common],
        help="Seal one published run artifact and an optional index artifact.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_EXAMPLES,
    )
    run.add_argument("--run-id", required=True)
    run.add_argument("--artifact-id", required=True)
    run.add_argument("--generation", required=True)
    run.add_argument("--revision", required=True)
    run.add_argument("--index-artifact-id")
    run.add_argument("--index-generation")
    run.add_argument("--index-revision")
    run.add_argument(
        "--relocation-root",
        type=Path,
        help=(
            "Read a switched copy of this exact publication from the relocation "
            "directory, including after the run publishes a newer generation."
        ),
    )
    return parser


def main(
    argv: Sequence[str] | None = None,
    *,
    sealer_factory: Callable[[WorkerSettings], object] | None = None,
) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)
    persistence = None
    try:
        settings = settings_from_args(args)
        if sealer_factory is None:
            sealer, persistence = compose_sealer(
                settings,
                relocation_root=getattr(args, "relocation_root", None),
            )
        else:
            sealer = sealer_factory(settings)
        if args.command == "restore":
            receipt = sealer.restore_version(  # type: ignore[attr-defined]
                project_id=args.project_id,
                experiment_id=args.experiment_id,
                archive_version_id=args.archive_version_id,
            )
        elif args.kind == "registered":
            receipt = sealer.seal_registered(  # type: ignore[attr-defined]
                project_id=args.project_id,
                experiment_id=args.experiment_id,
                experiment_revision_id=args.experiment_revision_id,
                track_key=args.track_key,
                display_name=args.display_name,
                assembly=args.assembly,
                track_type=args.track_type,
                result_id=args.result_id,
                result_revision_id=args.result_revision_id,
                index_result_id=args.index_result_id,
                index_result_revision_id=args.index_result_revision_id,
            )
        else:
            receipt = sealer.seal_run(  # type: ignore[attr-defined]
                project_id=args.project_id,
                experiment_id=args.experiment_id,
                experiment_revision_id=args.experiment_revision_id,
                track_key=args.track_key,
                display_name=args.display_name,
                assembly=args.assembly,
                track_type=args.track_type,
                run_id=args.run_id,
                artifact_id=args.artifact_id,
                generation=args.generation,
                revision=args.revision,
                index_artifact_id=args.index_artifact_id,
                index_generation=args.index_generation,
                index_revision=args.index_revision,
            )
        _write_json(require_operator_receipt(receipt))
        return 0
    except ArchiveSealError as exc:
        _write_error(exc)
        return 1
    finally:
        if persistence is not None:
            persistence.close()


def settings_from_args(args: argparse.Namespace) -> WorkerSettings:
    """Load WorkerSettings, with explicit flags overriding the environment."""

    environ = dict(os.environ)
    if args.database_url:
        environ[DATABASE_URL_ENV] = args.database_url
    if not str(environ.get(DATABASE_URL_ENV, "")).strip():
        raise ArchiveSealError(
            "archive_invalid",
            "Set --database-url or ENCODE_PIPELINE_DATABASE_URL.",
        )
    if args.archive_root is not None:
        environ[ARCHIVE_ROOT_ENV] = str(args.archive_root)
    if getattr(args, "storage_pool_config", None) is not None:
        environ[STORAGE_POOL_CONFIG_ENV] = str(args.storage_pool_config)
    if getattr(args, "workspace_root", None) is not None:
        environ[WORKSPACE_ROOT_ENV] = str(args.workspace_root)
    try:
        return load_worker_settings(environ)
    except ValueError as exc:
        message = str(exc)
        if not message or any(character in message for character in ("/", "\\")):
            message = "Worker settings are not valid."
        raise ArchiveSealError("archive_invalid", message) from None


def compose_sealer(
    settings: WorkerSettings,
    *,
    relocation_root: Path | None = None,
) -> tuple[object, object]:
    """Open the archive store from WorkerSettings.archive_root."""

    if settings.archive_root is None:
        raise ArchiveSealError(
            "archive_invalid",
            "ENCODE_PIPELINE_ARCHIVE_ROOT is not configured.",
        )
    if settings.archive_root.is_symlink() or not settings.archive_root.is_dir():
        raise ArchiveSealError(
            "archive_offline",
            "The archive root is unavailable.",
        )
    if relocation_root is not None and (
        not isinstance(relocation_root, Path)
        or not relocation_root.is_absolute()
        or any(part in {"", ".", ".."} for part in relocation_root.parts[1:])
    ):
        raise ArchiveSealError(
            "archive_invalid",
            "The relocation root is not a canonical absolute path.",
        )
    try:
        from encode_pipeline.persistence.archived_tracks import ArchiveRepository
        from encode_pipeline.services.archive_objects import ArchiveObjects
    except ImportError:
        raise ArchiveSealError(
            "archive_offline",
            "The archive store is unavailable.",
        ) from None
    try:
        from encode_pipeline.persistence.database import create_session_factory
        from encode_pipeline.persistence.registered_results import (
            SqlAlchemyRegisteredResultRepository,
        )
        from encode_pipeline.persistence.runtime import open_existing_run_persistence
        from encode_pipeline.services.artifact_downloads import ArtifactDownloadService
        from encode_pipeline.services.defaults import create_default_workflow_registry
        from encode_pipeline.services.runs import RunService

        persistence = open_existing_run_persistence(settings.database_url)
    except ArchiveSealError:
        raise
    except Exception:
        raise ArchiveSealError(
            "archive_offline",
            "The platform database is unavailable.",
        ) from None
    try:
        session_factory = create_session_factory(persistence.engine)
        sealer = _service(
            settings,
            objects=ArchiveObjects(settings.archive_root),
            repository=ArchiveRepository(session_factory),
            session_factory=session_factory,
            persistence=persistence,
            downloads=ArtifactDownloadService(
                run_service=RunService(
                    create_default_workflow_registry(),
                    repository=persistence.repository,
                ),
                workspace_root=settings.workspace_root,
            ),
            relocation_root=relocation_root,
            registered_revisions=SqlAlchemyRegisteredResultRepository(session_factory),
        )
    except ArchiveSealError:
        persistence.close()
        raise
    except Exception:
        persistence.close()
        raise ArchiveSealError(
            "archive_offline",
            "The archive store is unavailable.",
        ) from None
    return sealer, persistence


def _service(
    settings: WorkerSettings,
    *,
    objects: object,
    repository: object,
    session_factory: object,
    persistence: object,
    downloads: object,
    relocation_root: Path | None,
    registered_revisions: object,
) -> object:
    from encode_pipeline.services.archive_sealing import (
        ArchiveSealingService,
        CatalogExperimentGate,
        InputRegistryPoolLookup,
        ReferenceAssemblyLookup,
        SwitchedRunCopies,
    )

    return ArchiveSealingService(
        archive_root=settings.archive_root,
        objects=objects,
        repository=repository,
        experiments=CatalogExperimentGate(
            persistence.catalog_repository,  # type: ignore[attr-defined]
            session_factory,
        ),
        downloads=downloads,
        switched_copies=(
            SwitchedRunCopies(settings.workspace_root, relocation_root)
            if relocation_root is not None
            else None
        ),
        registered_revisions=registered_revisions,
        pool_lookup=InputRegistryPoolLookup(
            persistence.input_registry_repository  # type: ignore[attr-defined]
        ),
        pool_config_path=settings.storage_pool_config,
        assemblies=ReferenceAssemblyLookup(session_factory),
    )


def _scope_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "--database-url",
        help="File-backed SQLite URL. Overrides ENCODE_PIPELINE_DATABASE_URL.",
    )
    parser.add_argument(
        "--archive-root",
        type=Path,
        help="Absolute archive directory. Overrides ENCODE_PIPELINE_ARCHIVE_ROOT.",
    )
    parser.add_argument("--project-id", required=True)
    parser.add_argument("--experiment-id", required=True)
    return parser


def _common_parser(scope: argparse.ArgumentParser) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False, parents=[scope])
    parser.add_argument(
        "--storage-pool-config",
        type=Path,
        help=(
            "Private storage-pool map. Overrides "
            "ENCODE_PIPELINE_STORAGE_POOL_CONFIG. Required for registered bytes."
        ),
    )
    parser.add_argument(
        "--workspace-root",
        type=Path,
        help=(
            "Run workspace directory. Overrides ENCODE_PIPELINE_WORKSPACE_ROOT. "
            "Required when sealing a run artifact."
        ),
    )
    parser.add_argument("--experiment-revision-id", required=True)
    parser.add_argument("--track-key", required=True)
    parser.add_argument("--display-name", required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument(
        "--track-type",
        required=True,
        choices=("bigwig", "bed", "longrange"),
    )
    return parser


def _write_json(payload: Mapping[str, object]) -> None:
    json.dump(
        payload,
        sys.stdout,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    sys.stdout.write("\n")


def _write_error(error: ArchiveSealError) -> None:
    message = error.message if "/" not in error.message else "The seal failed."
    json.dump(
        {"error": {"code": error.code, "message": message}},
        sys.stderr,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    )
    sys.stderr.write("\n")
