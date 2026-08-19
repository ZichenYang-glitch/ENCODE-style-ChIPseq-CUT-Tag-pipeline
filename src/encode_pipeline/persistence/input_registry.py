"""SQLAlchemy implementation of the StoragePool and InputFile registry."""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import replace
from datetime import datetime, timezone

from sqlalchemy import select
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session, sessionmaker

from encode_pipeline.persistence.authentication import (
    insert_security_audit_event,
)
from encode_pipeline.platform.security_audit import SecurityAuditEvent
from encode_pipeline.persistence.models import (
    InputFileRevisionRow,
    InputFileRow,
    ProjectRow,
    ProjectStoragePoolBindingRow,
    StoragePoolRow,
)
from encode_pipeline.platform.data_registry import (
    LEGACY_PROJECT_ID,
    Project,
    ProjectKind,
)
from encode_pipeline.platform.input_registry import (
    InputFile,
    InputFileRevision,
    InputFileRevisionBindingRef,
    InputProvenanceMode,
    InputUseBindingEnvelope,
    InputUseBindingPlan,
    ProjectStoragePoolBinding,
    StoragePool,
    build_input_use_binding,
    build_input_use_binding_envelope,
)
from encode_pipeline.services.input_registry_repositories import (
    InputRegistryConflictError,
)


class SqlAlchemyInputRegistryRepository:
    """Persist logical pool and regular-file metadata without exposing rows."""

    def __init__(self, session_factory: sessionmaker[Session]) -> None:
        if not isinstance(session_factory, sessionmaker):
            raise ValueError("session_factory must be a SQLAlchemy sessionmaker")
        self._session_factory = session_factory

    def create_storage_pool(
        self,
        storage_pool: StoragePool,
        *,
        security_audit: SecurityAuditEvent | None = None,
    ) -> StoragePool:
        if not isinstance(storage_pool, StoragePool):
            raise ValueError("storage_pool must be a StoragePool")
        try:
            with self._session_factory.begin() as session:
                _begin_write(session)
                session.add(_storage_pool_row(storage_pool))
                session.flush()
                if security_audit is not None:
                    insert_security_audit_event(session, security_audit)
        except IntegrityError as exc:
            raise InputRegistryConflictError(
                "StoragePool ID or config_key already exists"
            ) from exc
        return storage_pool

    def get_storage_pool(self, storage_pool_id: str) -> StoragePool:
        with self._session_factory() as session:
            row = session.get(StoragePoolRow, storage_pool_id)
            if row is None:
                raise KeyError(f"Unknown storage_pool_id: {storage_pool_id!r}")
            return _storage_pool_from_row(row)

    def list_storage_pools(
        self,
        *,
        include_archived: bool = False,
    ) -> tuple[StoragePool, ...]:
        with self._session_factory() as session:
            statement = select(StoragePoolRow)
            if not include_archived:
                statement = statement.where(StoragePoolRow.archived_at.is_(None))
            rows = session.scalars(
                statement.order_by(
                    StoragePoolRow.created_at,
                    StoragePoolRow.storage_pool_id,
                )
            ).all()
            return tuple(_storage_pool_from_row(row) for row in rows)

    def archive_storage_pool(
        self,
        storage_pool_id: str,
        *,
        archived_at: datetime,
    ) -> StoragePool:
        with self._session_factory.begin() as session:
            _begin_write(session)
            row = session.get(StoragePoolRow, storage_pool_id)
            if row is None:
                raise KeyError(f"Unknown storage_pool_id: {storage_pool_id!r}")
            current = _storage_pool_from_row(row)
            if current.archived_at is not None:
                raise InputRegistryConflictError("StoragePool is already archived")
            archived = replace(current, archived_at=archived_at)
            row.archived_at = archived.archived_at
            session.flush()
            return archived

    def create_project_storage_pool_binding(
        self,
        binding: ProjectStoragePoolBinding,
    ) -> ProjectStoragePoolBinding:
        if not isinstance(binding, ProjectStoragePoolBinding):
            raise ValueError("binding must be a ProjectStoragePoolBinding")
        try:
            with self._session_factory.begin() as session:
                _begin_write(session)
                project_row = session.get(ProjectRow, binding.project_id)
                if project_row is None:
                    raise KeyError(f"Unknown project_id: {binding.project_id!r}")
                project = _project_from_row(project_row)
                _require_active_user_project(project)
                pool_row = session.get(
                    StoragePoolRow,
                    binding.storage_pool_id,
                )
                if pool_row is None:
                    raise KeyError(
                        f"Unknown storage_pool_id: {binding.storage_pool_id!r}"
                    )
                pool = _storage_pool_from_row(pool_row)
                _require_active_storage_pool(pool)
                if (
                    binding.bound_at < project.created_at
                    or binding.bound_at < pool.created_at
                ):
                    raise ValueError(
                        "Project StoragePool binding cannot precede its scope"
                    )
                if (
                    session.get(
                        ProjectStoragePoolBindingRow,
                        binding.project_id,
                    )
                    is not None
                ):
                    raise InputRegistryConflictError(
                        "Project is already bound to a StoragePool"
                    )
                session.add(_project_storage_pool_binding_row(binding))
                session.flush()
        except IntegrityError as exc:
            raise InputRegistryConflictError(
                "Project is already bound to a StoragePool"
            ) from exc
        return binding

    def get_project_storage_pool_binding(
        self,
        project_id: str,
    ) -> ProjectStoragePoolBinding:
        with self._session_factory() as session:
            row = session.get(ProjectStoragePoolBindingRow, project_id)
            if row is None:
                raise KeyError(f"Unknown Project StoragePool binding: {project_id!r}")
            return _project_storage_pool_binding_from_row(row)

    def get_active_project_storage_pool_binding(
        self,
        project_id: str,
    ) -> ProjectStoragePoolBinding:
        """Return one binding only while its Project and pool are active."""
        with self._session_factory() as session:
            _begin_consistent_read(session)
            project_row = session.get(ProjectRow, project_id)
            if project_row is None:
                raise KeyError(f"Unknown project_id: {project_id!r}")
            _require_active_user_project(_project_from_row(project_row))
            binding_row = session.get(ProjectStoragePoolBindingRow, project_id)
            if binding_row is None:
                raise KeyError(f"Unknown Project StoragePool binding: {project_id!r}")
            binding = _project_storage_pool_binding_from_row(binding_row)
            pool_row = session.get(StoragePoolRow, binding.storage_pool_id)
            if pool_row is None:
                raise ValueError("Project binding references a missing StoragePool")
            _require_active_storage_pool(_storage_pool_from_row(pool_row))
            return binding

    def create_input_file(
        self,
        input_file: InputFile,
        initial_revision: InputFileRevision,
    ) -> tuple[InputFile, InputFileRevision]:
        if not isinstance(input_file, InputFile):
            raise ValueError("input_file must be an InputFile")
        if not isinstance(initial_revision, InputFileRevision):
            raise ValueError("initial_revision must be an InputFileRevision")
        _validate_initial_revision(input_file, initial_revision)
        try:
            with self._session_factory.begin() as session:
                _begin_write(session)
                _require_active_registration_scope(session, input_file)
                session.add(_input_file_row(input_file))
                session.flush()
                session.add(_input_file_revision_row(initial_revision))
                session.flush()
        except IntegrityError as exc:
            raise InputRegistryConflictError(
                "InputFile stable_key or registry identity already exists"
            ) from exc
        return input_file, initial_revision

    def get_input_file(self, input_file_id: str) -> InputFile:
        with self._session_factory() as session:
            row = session.get(InputFileRow, input_file_id)
            if row is None:
                raise KeyError(f"Unknown input_file_id: {input_file_id!r}")
            return _input_file_from_row(row)

    def get_input_file_by_stable_key(
        self,
        project_id: str,
        stable_key: str,
    ) -> InputFile:
        with self._session_factory() as session:
            row = session.scalar(
                select(InputFileRow).where(
                    InputFileRow.project_id == project_id,
                    InputFileRow.stable_key == stable_key,
                )
            )
            if row is None:
                raise KeyError(
                    "Unknown InputFile stable_key for Project: "
                    f"{project_id!r}, {stable_key!r}"
                )
            return _input_file_from_row(row)

    def list_input_files(
        self,
        project_id: str,
        *,
        include_archived: bool = False,
    ) -> tuple[InputFile, ...]:
        with self._session_factory() as session:
            if session.get(ProjectRow, project_id) is None:
                raise KeyError(f"Unknown project_id: {project_id!r}")
            statement = select(InputFileRow).where(
                InputFileRow.project_id == project_id
            )
            if not include_archived:
                statement = statement.where(InputFileRow.archived_at.is_(None))
            rows = session.scalars(
                statement.order_by(
                    InputFileRow.created_at,
                    InputFileRow.input_file_id,
                )
            ).all()
            return tuple(_input_file_from_row(row) for row in rows)

    def archive_input_file(
        self,
        input_file_id: str,
        *,
        archived_at: datetime,
    ) -> InputFile:
        with self._session_factory.begin() as session:
            _begin_write(session)
            row = session.get(InputFileRow, input_file_id)
            if row is None:
                raise KeyError(f"Unknown input_file_id: {input_file_id!r}")
            current = _input_file_from_row(row)
            project_row = session.get(ProjectRow, current.project_id)
            if project_row is None:
                raise ValueError("InputFile references a missing Project")
            _require_active_user_project(_project_from_row(project_row))
            if current.archived_at is not None:
                raise InputRegistryConflictError("InputFile is already archived")
            archived = replace(current, archived_at=archived_at)
            row.archived_at = archived.archived_at
            session.flush()
            return archived

    def get_input_file_revision(
        self,
        input_file_revision_id: str,
    ) -> InputFileRevision:
        with self._session_factory() as session:
            row = session.get(
                InputFileRevisionRow,
                input_file_revision_id,
            )
            if row is None:
                raise KeyError(
                    f"Unknown input_file_revision_id: {input_file_revision_id!r}"
                )
            return _input_file_revision_from_row(row)

    def list_input_file_revisions(
        self,
        input_file_id: str,
    ) -> tuple[InputFileRevision, ...]:
        with self._session_factory() as session:
            if session.get(InputFileRow, input_file_id) is None:
                raise KeyError(f"Unknown input_file_id: {input_file_id!r}")
            rows = session.scalars(
                select(InputFileRevisionRow)
                .where(InputFileRevisionRow.input_file_id == input_file_id)
                .order_by(InputFileRevisionRow.revision_number)
            ).all()
            return tuple(_input_file_revision_from_row(row) for row in rows)

    def append_input_file_revision(
        self,
        revision: InputFileRevision,
        *,
        expected_previous_revision_number: int,
    ) -> InputFileRevision:
        if not isinstance(revision, InputFileRevision):
            raise ValueError("revision must be an InputFileRevision")
        if (
            not isinstance(expected_previous_revision_number, int)
            or isinstance(expected_previous_revision_number, bool)
            or expected_previous_revision_number < 1
        ):
            raise ValueError("expected_previous_revision_number must be positive")
        try:
            with self._session_factory.begin() as session:
                _begin_write(session)
                input_file_row = session.get(
                    InputFileRow,
                    revision.input_file_id,
                )
                if input_file_row is None:
                    raise KeyError(f"Unknown input_file_id: {revision.input_file_id!r}")
                input_file = _input_file_from_row(input_file_row)
                _require_active_registration_scope(session, input_file)
                if input_file.archived_at is not None:
                    raise ValueError("InputFile is archived")
                if (
                    revision.project_id != input_file.project_id
                    or revision.storage_pool_id != input_file.storage_pool_id
                ):
                    raise ValueError(
                        "revision must belong to the existing InputFile scope"
                    )
                if (
                    session.get(
                        InputFileRevisionRow,
                        revision.input_file_revision_id,
                    )
                    is not None
                ):
                    raise InputRegistryConflictError("Duplicate input_file_revision_id")
                previous_row = session.scalar(
                    select(InputFileRevisionRow)
                    .where(
                        InputFileRevisionRow.input_file_id == input_file.input_file_id
                    )
                    .order_by(InputFileRevisionRow.revision_number.desc())
                    .limit(1)
                )
                if previous_row is None:
                    raise ValueError("InputFile has no initial revision")
                previous = _input_file_revision_from_row(previous_row)
                if previous.revision_number != expected_previous_revision_number:
                    raise InputRegistryConflictError(
                        "InputFile revision changed concurrently"
                    )
                if revision.revision_number != previous.revision_number + 1:
                    raise ValueError(
                        "InputFile revision number must append exactly once"
                    )
                if revision.created_at < previous.created_at:
                    raise ValueError("new InputFile revision cannot precede history")
                session.add(_input_file_revision_row(revision))
                session.flush()
        except IntegrityError as exc:
            raise InputRegistryConflictError(
                "InputFile revision changed concurrently"
            ) from exc
        return revision

    def resolve_input_use_binding_plan(
        self,
        plan: InputUseBindingPlan,
        *,
        project_sample_binding_digest: str,
        workflow_inputs_digest: str,
    ) -> InputUseBindingEnvelope:
        """Resolve one plan in a consistent read transaction."""
        with self._session_factory() as session:
            _begin_consistent_read(session)
            return resolve_input_use_binding_plan_in_session(
                session,
                plan,
                project_sample_binding_digest=project_sample_binding_digest,
                workflow_inputs_digest=workflow_inputs_digest,
            )

    @contextmanager
    def locked_input_use_binding_plan(
        self,
        plan: InputUseBindingPlan,
        *,
        project_sample_binding_digest: str,
        workflow_inputs_digest: str,
    ) -> Iterator[InputUseBindingEnvelope]:
        """Serialize resolution until the caller stores snapshot evidence."""
        with self._session_factory.begin() as session:
            _begin_write(session)
            yield resolve_input_use_binding_plan_in_session(
                session,
                plan,
                project_sample_binding_digest=project_sample_binding_digest,
                workflow_inputs_digest=workflow_inputs_digest,
            )


def resolve_input_use_binding_plan_in_session(
    session: Session,
    plan: InputUseBindingPlan,
    *,
    project_sample_binding_digest: str,
    workflow_inputs_digest: str,
) -> InputUseBindingEnvelope:
    """Resolve input evidence inside an existing caller-owned transaction."""
    if not isinstance(session, Session):
        raise ValueError("session must be a SQLAlchemy Session")
    if not isinstance(plan, InputUseBindingPlan):
        raise ValueError("plan must be an InputUseBindingPlan")
    project_row = session.get(ProjectRow, plan.project_id)
    if project_row is None:
        raise KeyError(f"Unknown project_id: {plan.project_id!r}")
    _require_active_user_project(_project_from_row(project_row))
    binding: ProjectStoragePoolBinding | None = None
    if any(
        planned_use.provenance_mode is InputProvenanceMode.MANAGED_REVISION_V1
        for planned_use in plan.input_uses
    ):
        binding_row = session.get(
            ProjectStoragePoolBindingRow,
            plan.project_id,
        )
        if binding_row is None:
            raise KeyError(f"Unknown Project StoragePool binding: {plan.project_id!r}")
        binding = _project_storage_pool_binding_from_row(binding_row)
        pool_row = session.get(StoragePoolRow, binding.storage_pool_id)
        if pool_row is None:
            raise ValueError("Project binding references a missing StoragePool")
        _require_active_storage_pool(_storage_pool_from_row(pool_row))

    resolved_uses = []
    for planned_use in plan.input_uses:
        if planned_use.provenance_mode is InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1:
            resolved_uses.append(build_input_use_binding(planned_use, members=()))
            continue
        assert binding is not None
        if planned_use.closure_contract_version != "regular_file_v1":
            raise ValueError("managed closure contract is unsupported by this registry")
        members = []
        for revision_id in planned_use.input_file_revision_ids:
            revision_row = session.get(InputFileRevisionRow, revision_id)
            if revision_row is None:
                raise KeyError(f"Unknown input_file_revision_id: {revision_id!r}")
            revision = _input_file_revision_from_row(revision_row)
            input_file_row = session.get(
                InputFileRow,
                revision.input_file_id,
            )
            if input_file_row is None:
                raise ValueError("InputFile revision references a missing InputFile")
            input_file = _input_file_from_row(input_file_row)
            _require_active_binding_member(
                plan_project_id=plan.project_id,
                binding=binding,
                input_file=input_file,
                revision=revision,
            )
            members.append(_input_file_revision_binding_ref(revision))
        resolved_uses.append(
            build_input_use_binding(
                planned_use,
                members=tuple(members),
            )
        )
    return build_input_use_binding_envelope(
        project_id=plan.project_id,
        project_sample_binding_digest=project_sample_binding_digest,
        workflow_id=plan.workflow_id,
        adapter_contract_version=plan.adapter_contract_version,
        workflow_inputs_digest=workflow_inputs_digest,
        input_uses=tuple(resolved_uses),
    )


def _require_active_registration_scope(
    session: Session,
    input_file: InputFile,
) -> None:
    project_row = session.get(ProjectRow, input_file.project_id)
    if project_row is None:
        raise KeyError(f"Unknown project_id: {input_file.project_id!r}")
    project = _project_from_row(project_row)
    _require_active_user_project(project)
    pool_row = session.get(StoragePoolRow, input_file.storage_pool_id)
    if pool_row is None:
        raise KeyError(f"Unknown storage_pool_id: {input_file.storage_pool_id!r}")
    pool = _storage_pool_from_row(pool_row)
    _require_active_storage_pool(pool)
    binding_row = session.get(
        ProjectStoragePoolBindingRow,
        input_file.project_id,
    )
    if binding_row is None:
        raise KeyError(
            f"Unknown Project StoragePool binding: {input_file.project_id!r}"
        )
    if binding_row.storage_pool_id != input_file.storage_pool_id:
        raise ValueError("InputFile StoragePool does not match the Project binding")
    binding = _project_storage_pool_binding_from_row(binding_row)
    if (
        input_file.created_at < project.created_at
        or input_file.created_at < pool.created_at
        or input_file.created_at < binding.bound_at
    ):
        raise ValueError("InputFile cannot precede its authorization scope")


def _require_active_binding_member(
    *,
    plan_project_id: str,
    binding: ProjectStoragePoolBinding,
    input_file: InputFile,
    revision: InputFileRevision,
) -> None:
    if input_file.archived_at is not None:
        raise ValueError("InputFile is archived")
    if (
        input_file.project_id != plan_project_id
        or input_file.storage_pool_id != binding.storage_pool_id
        or revision.project_id != plan_project_id
        or revision.storage_pool_id != binding.storage_pool_id
        or revision.input_file_id != input_file.input_file_id
    ):
        raise ValueError(
            "InputFile revision is not authorized by the Project pool binding"
        )


def _input_file_revision_binding_ref(
    revision: InputFileRevision,
) -> InputFileRevisionBindingRef:
    """Map one regular file to the fixed generic member coordinate."""
    return InputFileRevisionBindingRef(
        logical_member_key="file",
        input_file_id=revision.input_file_id,
        input_file_revision_id=revision.input_file_revision_id,
        revision_digest=revision.digest,
        size_bytes=revision.size_bytes,
        content_sha256=revision.content_sha256,
    )


def _validate_initial_revision(
    input_file: InputFile,
    revision: InputFileRevision,
) -> None:
    if (
        revision.input_file_id != input_file.input_file_id
        or revision.project_id != input_file.project_id
        or revision.storage_pool_id != input_file.storage_pool_id
        or revision.revision_number != 1
    ):
        raise ValueError("initial revision must be revision 1 of the new InputFile")
    if revision.created_at < input_file.created_at:
        raise ValueError("initial revision cannot precede InputFile creation")


def _require_active_user_project(project: Project) -> None:
    if project.project_id == LEGACY_PROJECT_ID:
        raise ValueError("Legacy Project cannot receive trusted InputFiles")
    if project.archived_at is not None:
        raise ValueError("Project is archived")


def _require_active_storage_pool(storage_pool: StoragePool) -> None:
    if storage_pool.archived_at is not None:
        raise ValueError("StoragePool is archived")


def _begin_write(session: Session) -> None:
    """Serialize SQLite read-then-write registry operations."""
    if session.get_bind().dialect.name == "sqlite":
        session.connection().exec_driver_sql("BEGIN IMMEDIATE")


def _begin_consistent_read(session: Session) -> None:
    """Pin one SQLite snapshot across an input binding resolution."""
    if session.get_bind().dialect.name == "sqlite":
        session.connection().exec_driver_sql("BEGIN")


def _storage_pool_row(storage_pool: StoragePool) -> StoragePoolRow:
    return StoragePoolRow(
        storage_pool_id=storage_pool.storage_pool_id,
        config_key=storage_pool.config_key,
        display_name=storage_pool.display_name,
        created_at=storage_pool.created_at,
        archived_at=storage_pool.archived_at,
    )


def _storage_pool_from_row(row: StoragePoolRow) -> StoragePool:
    return StoragePool(
        storage_pool_id=row.storage_pool_id,
        config_key=row.config_key,
        display_name=row.display_name,
        created_at=_as_utc(row.created_at),
        archived_at=_optional_utc(row.archived_at),
    )


def _project_storage_pool_binding_row(
    binding: ProjectStoragePoolBinding,
) -> ProjectStoragePoolBindingRow:
    return ProjectStoragePoolBindingRow(
        project_id=binding.project_id,
        storage_pool_id=binding.storage_pool_id,
        bound_at=binding.bound_at,
    )


def _project_storage_pool_binding_from_row(
    row: ProjectStoragePoolBindingRow,
) -> ProjectStoragePoolBinding:
    return ProjectStoragePoolBinding(
        project_id=row.project_id,
        storage_pool_id=row.storage_pool_id,
        bound_at=_as_utc(row.bound_at),
    )


def _input_file_row(input_file: InputFile) -> InputFileRow:
    return InputFileRow(
        input_file_id=input_file.input_file_id,
        project_id=input_file.project_id,
        storage_pool_id=input_file.storage_pool_id,
        stable_key=input_file.stable_key,
        created_at=input_file.created_at,
        archived_at=input_file.archived_at,
    )


def _input_file_from_row(row: InputFileRow) -> InputFile:
    return InputFile(
        input_file_id=row.input_file_id,
        project_id=row.project_id,
        storage_pool_id=row.storage_pool_id,
        stable_key=row.stable_key,
        created_at=_as_utc(row.created_at),
        archived_at=_optional_utc(row.archived_at),
    )


def _input_file_revision_row(
    revision: InputFileRevision,
) -> InputFileRevisionRow:
    return InputFileRevisionRow(
        input_file_revision_id=revision.input_file_revision_id,
        input_file_id=revision.input_file_id,
        project_id=revision.project_id,
        storage_pool_id=revision.storage_pool_id,
        revision_number=revision.revision_number,
        relative_path=revision.relative_path,
        size_bytes=revision.size_bytes,
        content_sha256=revision.content_sha256,
        digest_scheme=revision.digest_scheme,
        digest=revision.digest,
        created_at=revision.created_at,
    )


def _input_file_revision_from_row(
    row: InputFileRevisionRow,
) -> InputFileRevision:
    return InputFileRevision(
        input_file_revision_id=row.input_file_revision_id,
        input_file_id=row.input_file_id,
        project_id=row.project_id,
        storage_pool_id=row.storage_pool_id,
        revision_number=row.revision_number,
        relative_path=row.relative_path,
        size_bytes=row.size_bytes,
        content_sha256=row.content_sha256,
        digest_scheme=row.digest_scheme,
        digest=row.digest,
        created_at=_as_utc(row.created_at),
    )


def _project_from_row(row: ProjectRow) -> Project:
    return Project(
        project_id=row.project_id,
        display_name=row.display_name,
        kind=ProjectKind(row.kind),
        created_at=_as_utc(row.created_at),
        archived_at=_optional_utc(row.archived_at),
    )


def _as_utc(value: datetime) -> datetime:
    if value.tzinfo is None or value.utcoffset() is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _optional_utc(value: datetime | None) -> datetime | None:
    return None if value is None else _as_utc(value)
