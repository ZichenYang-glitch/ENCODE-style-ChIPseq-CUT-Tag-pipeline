"""SQLAlchemy persistence for stable, append-only Reference Profiles."""

from __future__ import annotations

from datetime import datetime, timezone

from sqlalchemy import select
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session, sessionmaker

from encode_pipeline.persistence.models import (
    ReferenceProfileRevisionRow,
    ReferenceProfileRow,
    ReferenceProfileWorkflowBindingRow,
)
from encode_pipeline.platform.reference_profiles import (
    ReferenceProfile,
    ReferenceProfileRevision,
    ReferenceProfileWorkflowBinding,
)
from encode_pipeline.services.reference_profile_repositories import (
    ReferenceProfileConflictError,
)


class SqlAlchemyReferenceProfileRepository:
    """Persist private catalog coordinates without exposing ORM rows."""

    def __init__(self, session_factory: sessionmaker[Session]) -> None:
        if not isinstance(session_factory, sessionmaker):
            raise ValueError("session_factory must be a SQLAlchemy sessionmaker")
        self._session_factory = session_factory

    def create_profile(
        self,
        profile: ReferenceProfile,
        initial_revision: ReferenceProfileRevision,
    ) -> tuple[ReferenceProfile, ReferenceProfileRevision]:
        _validate_initial(profile, initial_revision)
        try:
            with self._session_factory.begin() as session:
                _begin_write(session)
                session.add(_profile_row(profile))
                session.flush()
                session.add(_revision_row(initial_revision))
                session.flush()
                session.add_all(_binding_rows(initial_revision))
                session.flush()
        except IntegrityError as exc:
            raise ReferenceProfileConflictError(
                "reference profile identity already exists"
            ) from exc
        return profile, initial_revision

    def get_profile(self, profile_id: str) -> ReferenceProfile:
        with self._session_factory() as session:
            row = session.get(ReferenceProfileRow, profile_id)
            if row is None:
                raise KeyError("unknown reference profile")
            return _profile_from_row(row)

    def get_profile_by_safe_key(self, safe_key: str) -> ReferenceProfile:
        with self._session_factory() as session:
            row = session.scalar(
                select(ReferenceProfileRow).where(
                    ReferenceProfileRow.safe_key == safe_key
                )
            )
            if row is None:
                raise KeyError("unknown reference profile")
            return _profile_from_row(row)

    def get_profile_for_revision(self, revision_id: str) -> ReferenceProfile:
        with self._session_factory() as session:
            revision = session.get(ReferenceProfileRevisionRow, revision_id)
            if revision is None:
                raise KeyError("unknown reference profile revision")
            profile = session.get(ReferenceProfileRow, revision.profile_id)
            if profile is None:
                raise ValueError("reference revision has no profile")
            return _profile_from_row(profile)

    def list_profiles(self) -> tuple[ReferenceProfile, ...]:
        with self._session_factory() as session:
            rows = session.scalars(
                select(ReferenceProfileRow).order_by(
                    ReferenceProfileRow.created_at,
                    ReferenceProfileRow.profile_id,
                )
            ).all()
            return tuple(_profile_from_row(row) for row in rows)

    def get_revision(self, revision_id: str) -> ReferenceProfileRevision:
        with self._session_factory() as session:
            row = session.get(ReferenceProfileRevisionRow, revision_id)
            if row is None:
                raise KeyError("unknown reference profile revision")
            return _revision_from_row(session, row)

    def list_revisions(
        self,
        profile_id: str,
    ) -> tuple[ReferenceProfileRevision, ...]:
        with self._session_factory() as session:
            if session.get(ReferenceProfileRow, profile_id) is None:
                raise KeyError("unknown reference profile")
            rows = session.scalars(
                select(ReferenceProfileRevisionRow)
                .where(ReferenceProfileRevisionRow.profile_id == profile_id)
                .order_by(ReferenceProfileRevisionRow.revision_number)
            ).all()
            return tuple(_revision_from_row(session, row) for row in rows)

    def append_revision(
        self,
        revision: ReferenceProfileRevision,
        *,
        expected_previous_revision_number: int,
    ) -> ReferenceProfileRevision:
        if not isinstance(revision, ReferenceProfileRevision):
            raise ValueError("revision must be a ReferenceProfileRevision")
        if (
            not isinstance(expected_previous_revision_number, int)
            or isinstance(expected_previous_revision_number, bool)
            or expected_previous_revision_number < 1
        ):
            raise ValueError("expected_previous_revision_number must be positive")
        try:
            with self._session_factory.begin() as session:
                _begin_write(session)
                profile = session.get(ReferenceProfileRow, revision.profile_id)
                if profile is None:
                    raise KeyError("unknown reference profile")
                current = session.scalar(
                    select(ReferenceProfileRevisionRow)
                    .where(
                        ReferenceProfileRevisionRow.profile_id == revision.profile_id
                    )
                    .order_by(ReferenceProfileRevisionRow.revision_number.desc())
                    .limit(1)
                )
                if current is None:
                    raise ValueError("reference profile has no revision")
                if current.revision_number != expected_previous_revision_number:
                    raise ReferenceProfileConflictError(
                        "reference revision changed concurrently"
                    )
                if revision.revision_number != current.revision_number + 1:
                    raise ValueError("revision number must be append-only")
                if _as_utc(revision.created_at) < _as_utc(current.created_at):
                    raise ValueError("revision cannot precede previous revision")
                session.add(_revision_row(revision))
                session.flush()
                session.add_all(_binding_rows(revision))
                session.flush()
        except IntegrityError as exc:
            raise ReferenceProfileConflictError(
                "reference revision identity already exists"
            ) from exc
        return revision

    def set_enabled_revision(
        self,
        profile_id: str,
        revision_id: str | None,
    ) -> ReferenceProfile:
        try:
            with self._session_factory.begin() as session:
                _begin_write(session)
                profile = session.get(ReferenceProfileRow, profile_id)
                if profile is None:
                    raise KeyError("unknown reference profile")
                if revision_id is not None:
                    revision = session.get(ReferenceProfileRevisionRow, revision_id)
                    if revision is None or revision.profile_id != profile_id:
                        raise ValueError("enabled revision does not belong to profile")
                profile.enabled_revision_id = revision_id
                session.flush()
                result = _profile_from_row(profile)
        except IntegrityError as exc:
            raise ReferenceProfileConflictError(
                "enabled reference revision conflicted"
            ) from exc
        return result

    def list_enabled_for_workflow(
        self,
        workflow_id: str,
    ) -> tuple[tuple[ReferenceProfile, ReferenceProfileRevision], ...]:
        with self._session_factory() as session:
            rows = session.execute(
                select(ReferenceProfileRow, ReferenceProfileRevisionRow)
                .join(
                    ReferenceProfileRevisionRow,
                    ReferenceProfileRevisionRow.revision_id
                    == ReferenceProfileRow.enabled_revision_id,
                )
                .join(
                    ReferenceProfileWorkflowBindingRow,
                    (
                        ReferenceProfileWorkflowBindingRow.revision_id
                        == ReferenceProfileRevisionRow.revision_id
                    )
                    & (
                        ReferenceProfileWorkflowBindingRow.profile_id
                        == ReferenceProfileRevisionRow.profile_id
                    ),
                )
                .where(ReferenceProfileWorkflowBindingRow.workflow_id == workflow_id)
                .order_by(
                    ReferenceProfileRow.created_at,
                    ReferenceProfileRow.profile_id,
                )
            ).all()
            return tuple(
                (_profile_from_row(profile), _revision_from_row(session, revision))
                for profile, revision in rows
            )


def _validate_initial(
    profile: ReferenceProfile,
    revision: ReferenceProfileRevision,
) -> None:
    if not isinstance(profile, ReferenceProfile):
        raise ValueError("profile must be a ReferenceProfile")
    if not isinstance(revision, ReferenceProfileRevision):
        raise ValueError("initial_revision must be a ReferenceProfileRevision")
    if profile.enabled_revision_id is not None:
        raise ValueError("a new ReferenceProfile must begin disabled")
    if (
        revision.profile_id != profile.profile_id
        or revision.revision_number != 1
        or revision.created_at < profile.created_at
    ):
        raise ValueError("initial revision does not match ReferenceProfile")


def _profile_row(profile: ReferenceProfile) -> ReferenceProfileRow:
    return ReferenceProfileRow(
        profile_id=profile.profile_id,
        safe_key=profile.safe_key,
        created_at=profile.created_at,
        enabled_revision_id=profile.enabled_revision_id,
    )


def _profile_from_row(row: ReferenceProfileRow) -> ReferenceProfile:
    return ReferenceProfile(
        profile_id=row.profile_id,
        safe_key=row.safe_key,
        created_at=_as_utc(row.created_at),
        enabled_revision_id=row.enabled_revision_id,
    )


def _revision_row(revision: ReferenceProfileRevision) -> ReferenceProfileRevisionRow:
    return ReferenceProfileRevisionRow(
        revision_id=revision.revision_id,
        profile_id=revision.profile_id,
        revision_number=revision.revision_number,
        display_name=revision.display_name,
        organism=revision.organism,
        assembly=revision.assembly,
        config_key=revision.config_key,
        public_identity_scheme=revision.public_identity_scheme,
        public_identity_sha256=revision.public_identity_sha256,
        created_at=revision.created_at,
    )


def _binding_rows(
    revision: ReferenceProfileRevision,
) -> tuple[ReferenceProfileWorkflowBindingRow, ...]:
    return tuple(
        ReferenceProfileWorkflowBindingRow(
            profile_id=revision.profile_id,
            revision_id=revision.revision_id,
            workflow_id=binding.workflow_id,
            contract_version=binding.contract_version,
            identity_scheme=binding.identity_scheme,
            identity_sha256=binding.identity_sha256,
        )
        for binding in revision.workflow_bindings
    )


def _revision_from_row(
    session: Session,
    row: ReferenceProfileRevisionRow,
) -> ReferenceProfileRevision:
    bindings = session.scalars(
        select(ReferenceProfileWorkflowBindingRow)
        .where(
            ReferenceProfileWorkflowBindingRow.profile_id == row.profile_id,
            ReferenceProfileWorkflowBindingRow.revision_id == row.revision_id,
        )
        .order_by(ReferenceProfileWorkflowBindingRow.workflow_id)
    ).all()
    return ReferenceProfileRevision(
        revision_id=row.revision_id,
        profile_id=row.profile_id,
        revision_number=row.revision_number,
        display_name=row.display_name,
        organism=row.organism,
        assembly=row.assembly,
        config_key=row.config_key,
        workflow_bindings=tuple(
            ReferenceProfileWorkflowBinding(
                workflow_id=binding.workflow_id,
                contract_version=binding.contract_version,
                identity_scheme=binding.identity_scheme,
                identity_sha256=binding.identity_sha256,
            )
            for binding in bindings
        ),
        public_identity_scheme=row.public_identity_scheme,
        public_identity_sha256=row.public_identity_sha256,
        created_at=_as_utc(row.created_at),
    )


def _as_utc(value: datetime) -> datetime:
    if value.tzinfo is None or value.utcoffset() is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _begin_write(session: Session) -> None:
    if session.get_bind().dialect.name == "sqlite":
        session.connection().exec_driver_sql("BEGIN IMMEDIATE")
