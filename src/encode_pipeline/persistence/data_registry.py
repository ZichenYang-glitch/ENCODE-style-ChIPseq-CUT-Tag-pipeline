"""SQLAlchemy implementation of the Project and Sample registry boundary."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timezone

from sqlalchemy import select
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session, sessionmaker

from encode_pipeline.persistence.models import (
    ProjectRow,
    SampleRevisionRow,
    SampleRow,
)
from encode_pipeline.platform.data_registry import (
    LEGACY_PROJECT_ID,
    BindingMode,
    BindingProvenance,
    Project,
    ProjectKind,
    ProjectSampleBinding,
    ProjectSampleSelection,
    Sample,
    SampleRevision,
    SampleRevisionBindingRef,
    build_project_sample_binding,
)
from encode_pipeline.services.data_registry_repositories import (
    DataRegistryConflictError,
)


class SqlAlchemyDataRegistryRepository:
    """Persist registry objects without exposing ORM rows to callers."""

    def __init__(self, session_factory: sessionmaker[Session]) -> None:
        if not isinstance(session_factory, sessionmaker):
            raise ValueError("session_factory must be a SQLAlchemy sessionmaker")
        self._session_factory = session_factory

    def create_project(self, project: Project) -> Project:
        if not isinstance(project, Project):
            raise ValueError("project must be a Project")
        if project.project_id == LEGACY_PROJECT_ID:
            raise DataRegistryConflictError("reserved Legacy Project already exists")
        try:
            with self._session_factory.begin() as session:
                _begin_write(session)
                session.add(_project_row(project))
                session.flush()
        except IntegrityError as exc:
            raise DataRegistryConflictError(
                f"Duplicate project_id: {project.project_id!r}"
            ) from exc
        return project

    def get_project(self, project_id: str) -> Project:
        with self._session_factory() as session:
            row = session.get(ProjectRow, project_id)
            if row is None:
                raise KeyError(f"Unknown project_id: {project_id!r}")
            return _project_from_row(row)

    def list_projects(
        self,
        *,
        include_archived: bool = False,
    ) -> tuple[Project, ...]:
        with self._session_factory() as session:
            statement = select(ProjectRow)
            if not include_archived:
                statement = statement.where(ProjectRow.archived_at.is_(None))
            rows = session.scalars(
                statement.order_by(ProjectRow.created_at, ProjectRow.project_id)
            ).all()
            return tuple(_project_from_row(row) for row in rows)

    def archive_project(
        self,
        project_id: str,
        *,
        archived_at: datetime,
    ) -> Project:
        with self._session_factory.begin() as session:
            _begin_write(session)
            row = session.get(ProjectRow, project_id)
            if row is None:
                raise KeyError(f"Unknown project_id: {project_id!r}")
            current = _project_from_row(row)
            if project_id == LEGACY_PROJECT_ID:
                raise ValueError("reserved Legacy Project cannot be archived")
            if current.archived_at is not None:
                raise DataRegistryConflictError("Project is already archived")
            archived = replace(current, archived_at=archived_at)
            row.archived_at = archived.archived_at
            session.flush()
            return archived

    def create_sample(
        self,
        sample: Sample,
        initial_revision: SampleRevision,
    ) -> tuple[Sample, SampleRevision]:
        return self.create_samples(((sample, initial_revision),))[0]

    def create_samples(
        self,
        entries: tuple[tuple[Sample, SampleRevision], ...],
    ) -> tuple[tuple[Sample, SampleRevision], ...]:
        entries = tuple(entries)
        try:
            with self._session_factory.begin() as session:
                _begin_write(session)
                self._validate_new_sample_entries(session, entries)
                for sample, _initial_revision in entries:
                    session.add(_sample_row(sample))
                session.flush()
                for _sample, initial_revision in entries:
                    session.add(_sample_revision_row(initial_revision))
                session.flush()
        except IntegrityError as exc:
            raise DataRegistryConflictError(
                "Sample stable_key or registry identity already exists"
            ) from exc
        return entries

    def get_sample(self, sample_id: str) -> Sample:
        with self._session_factory() as session:
            row = session.get(SampleRow, sample_id)
            if row is None:
                raise KeyError(f"Unknown sample_id: {sample_id!r}")
            return _sample_from_row(row)

    def list_samples(self, project_id: str) -> tuple[Sample, ...]:
        with self._session_factory() as session:
            if session.get(ProjectRow, project_id) is None:
                raise KeyError(f"Unknown project_id: {project_id!r}")
            rows = session.scalars(
                select(SampleRow)
                .where(SampleRow.project_id == project_id)
                .order_by(SampleRow.created_at, SampleRow.sample_id)
            ).all()
            return tuple(_sample_from_row(row) for row in rows)

    def get_sample_revision(self, sample_revision_id: str) -> SampleRevision:
        with self._session_factory() as session:
            row = session.get(SampleRevisionRow, sample_revision_id)
            if row is None:
                raise KeyError(f"Unknown sample_revision_id: {sample_revision_id!r}")
            return _sample_revision_from_row(row)

    def list_sample_revisions(
        self,
        sample_id: str,
    ) -> tuple[SampleRevision, ...]:
        with self._session_factory() as session:
            if session.get(SampleRow, sample_id) is None:
                raise KeyError(f"Unknown sample_id: {sample_id!r}")
            rows = session.scalars(
                select(SampleRevisionRow)
                .where(SampleRevisionRow.sample_id == sample_id)
                .order_by(SampleRevisionRow.revision_number)
            ).all()
            return tuple(_sample_revision_from_row(row) for row in rows)

    def append_sample_revision(
        self,
        revision: SampleRevision,
        *,
        expected_previous_revision_number: int,
    ) -> SampleRevision:
        if not isinstance(revision, SampleRevision):
            raise ValueError("revision must be a SampleRevision")
        if (
            not isinstance(expected_previous_revision_number, int)
            or isinstance(expected_previous_revision_number, bool)
            or expected_previous_revision_number < 1
        ):
            raise ValueError("expected_previous_revision_number must be positive")
        try:
            with self._session_factory.begin() as session:
                _begin_write(session)
                sample_row = session.get(SampleRow, revision.sample_id)
                if sample_row is None:
                    raise KeyError(f"Unknown sample_id: {revision.sample_id!r}")
                sample = _sample_from_row(sample_row)
                project_row = session.get(ProjectRow, sample.project_id)
                if project_row is None:
                    raise ValueError("Sample references a missing Project")
                _require_active_user_project(_project_from_row(project_row))
                if (
                    revision.project_id != sample.project_id
                    or revision.sample_id != sample.sample_id
                ):
                    raise ValueError("revision must belong to the existing Sample")
                if (
                    session.get(
                        SampleRevisionRow,
                        revision.sample_revision_id,
                    )
                    is not None
                ):
                    raise DataRegistryConflictError("Duplicate sample_revision_id")
                previous_row = session.scalar(
                    select(SampleRevisionRow)
                    .where(SampleRevisionRow.sample_id == sample.sample_id)
                    .order_by(SampleRevisionRow.revision_number.desc())
                    .limit(1)
                )
                if previous_row is None:
                    raise ValueError("Sample has no initial revision")
                previous = _sample_revision_from_row(previous_row)
                if previous.revision_number != expected_previous_revision_number:
                    raise DataRegistryConflictError(
                        "Sample revision changed concurrently"
                    )
                if revision.revision_number != previous.revision_number + 1:
                    raise ValueError("Sample revision number must append exactly once")
                if revision.created_at < previous.created_at:
                    raise ValueError("new Sample revision cannot precede history")
                session.add(_sample_revision_row(revision))
                session.flush()
        except IntegrityError as exc:
            raise DataRegistryConflictError(
                "Sample revision changed concurrently"
            ) from exc
        return revision

    def resolve_project_sample_selection(
        self,
        selection: ProjectSampleSelection,
        *,
        workflow_inputs_digest: str,
    ) -> ProjectSampleBinding:
        if not isinstance(selection, ProjectSampleSelection):
            raise ValueError("selection must be a ProjectSampleSelection")
        with self._session_factory() as session:
            _begin_consistent_read(session)
            project_row = session.get(ProjectRow, selection.project_id)
            if project_row is None:
                raise KeyError(f"Unknown project_id: {selection.project_id!r}")
            _require_active_user_project(_project_from_row(project_row))
            rows = session.scalars(
                select(SampleRevisionRow).where(
                    SampleRevisionRow.sample_revision_id.in_(
                        selection.sample_revision_ids
                    )
                )
            ).all()
            by_id = {row.sample_revision_id: row for row in rows}
            revisions: list[SampleRevision] = []
            for revision_id in selection.sample_revision_ids:
                row = by_id.get(revision_id)
                if row is None:
                    raise KeyError(f"Unknown sample_revision_id: {revision_id!r}")
                revision = _sample_revision_from_row(row)
                if revision.project_id != selection.project_id:
                    raise ValueError(
                        "all SampleRevisions must belong to the same Project"
                    )
                revisions.append(revision)
            return build_project_sample_binding(
                project_id=selection.project_id,
                binding_mode=BindingMode.BOUND_V1,
                provenance=BindingProvenance.RESOLVED,
                workflow_inputs_digest=workflow_inputs_digest,
                sample_revisions=tuple(
                    SampleRevisionBindingRef(
                        sample_revision_id=revision.sample_revision_id,
                        payload_digest=revision.payload_digest,
                    )
                    for revision in revisions
                ),
            )

    @staticmethod
    def _validate_new_sample_entries(
        session: Session,
        entries: tuple[tuple[Sample, SampleRevision], ...],
    ) -> None:
        pending_sample_ids: set[str] = set()
        pending_revision_ids: set[str] = set()
        pending_stable_identities: set[tuple[str, str]] = set()
        projects: dict[str, Project] = {}
        for sample, initial_revision in entries:
            if not isinstance(sample, Sample):
                raise ValueError("sample must be a Sample")
            if not isinstance(initial_revision, SampleRevision):
                raise ValueError("initial_revision must be a SampleRevision")
            project = projects.get(sample.project_id)
            if project is None:
                project_row = session.get(ProjectRow, sample.project_id)
                if project_row is None:
                    raise KeyError(f"Unknown project_id: {sample.project_id!r}")
                project = _project_from_row(project_row)
                _require_active_user_project(project)
                projects[sample.project_id] = project
            if (
                session.get(SampleRow, sample.sample_id) is not None
                or sample.sample_id in pending_sample_ids
            ):
                raise DataRegistryConflictError(
                    f"Duplicate sample_id: {sample.sample_id!r}"
                )
            stable_identity = (sample.project_id, sample.stable_key)
            if (
                stable_identity in pending_stable_identities
                or session.scalar(
                    select(SampleRow.sample_id).where(
                        SampleRow.project_id == sample.project_id,
                        SampleRow.stable_key == sample.stable_key,
                    )
                )
                is not None
            ):
                raise DataRegistryConflictError(
                    "Sample stable_key already exists in this Project"
                )
            if (
                session.get(
                    SampleRevisionRow,
                    initial_revision.sample_revision_id,
                )
                is not None
                or initial_revision.sample_revision_id in pending_revision_ids
            ):
                raise DataRegistryConflictError("Duplicate sample_revision_id")
            if (
                initial_revision.project_id != sample.project_id
                or initial_revision.sample_id != sample.sample_id
                or initial_revision.revision_number != 1
            ):
                raise ValueError(
                    "initial revision must be revision 1 of the new Sample"
                )
            if initial_revision.created_at < sample.created_at:
                raise ValueError("initial revision cannot precede Sample creation")
            pending_sample_ids.add(sample.sample_id)
            pending_revision_ids.add(initial_revision.sample_revision_id)
            pending_stable_identities.add(stable_identity)


def _require_active_user_project(project: Project) -> None:
    if project.project_id == LEGACY_PROJECT_ID:
        raise ValueError("Legacy Project cannot receive trusted registry data")
    if project.archived_at is not None:
        raise ValueError("Project is archived")


def _begin_write(session: Session) -> None:
    """Serialize SQLite read-then-write registry operations."""
    if session.get_bind().dialect.name == "sqlite":
        session.connection().exec_driver_sql("BEGIN IMMEDIATE")


def _begin_consistent_read(session: Session) -> None:
    """Pin one SQLite snapshot across a selection resolution."""
    if session.get_bind().dialect.name == "sqlite":
        session.connection().exec_driver_sql("BEGIN")


def _project_row(project: Project) -> ProjectRow:
    return ProjectRow(
        project_id=project.project_id,
        display_name=project.display_name,
        kind=(
            project.kind.value
            if isinstance(project.kind, ProjectKind)
            else project.kind
        ),
        created_at=project.created_at,
        archived_at=project.archived_at,
    )


def _project_from_row(row: ProjectRow) -> Project:
    return Project(
        project_id=row.project_id,
        display_name=row.display_name,
        kind=row.kind,
        created_at=_as_utc(row.created_at),
        archived_at=_optional_utc(row.archived_at),
    )


def _sample_row(sample: Sample) -> SampleRow:
    return SampleRow(
        sample_id=sample.sample_id,
        project_id=sample.project_id,
        stable_key=sample.stable_key,
        created_at=sample.created_at,
    )


def _sample_from_row(row: SampleRow) -> Sample:
    return Sample(
        sample_id=row.sample_id,
        project_id=row.project_id,
        stable_key=row.stable_key,
        created_at=_as_utc(row.created_at),
    )


def _sample_revision_row(revision: SampleRevision) -> SampleRevisionRow:
    return SampleRevisionRow(
        sample_revision_id=revision.sample_revision_id,
        project_id=revision.project_id,
        sample_id=revision.sample_id,
        revision_number=revision.revision_number,
        canonical_payload=revision.canonical_payload,
        payload_digest_scheme=revision.payload_digest_scheme,
        payload_digest=revision.payload_digest,
        created_at=revision.created_at,
    )


def _sample_revision_from_row(row: SampleRevisionRow) -> SampleRevision:
    return SampleRevision(
        sample_revision_id=row.sample_revision_id,
        project_id=row.project_id,
        sample_id=row.sample_id,
        revision_number=row.revision_number,
        canonical_payload=row.canonical_payload,
        payload_digest_scheme=row.payload_digest_scheme,
        payload_digest=row.payload_digest,
        created_at=_as_utc(row.created_at),
    )


def _as_utc(value: datetime) -> datetime:
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _optional_utc(value: datetime | None) -> datetime | None:
    return _as_utc(value) if value is not None else None
