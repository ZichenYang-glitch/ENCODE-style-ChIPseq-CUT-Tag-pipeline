"""SQLite publication and scoped lookup of immutable experiment track versions."""

from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
from uuid import uuid4

from sqlalchemy import func, select
from sqlalchemy.orm import defer, Session, sessionmaker

from encode_pipeline.persistence.models import (
    ArchivedTrackVersionRow as Row,
    ArtifactPublicationRow,
    ExperimentRevisionRow,
    ExperimentRow,
    ProjectRow,
    RegisteredResultRevisionRow,
    RegisteredResultRow,
    RunExperimentBindingRow,
    RunProjectBindingRow,
)
from encode_pipeline.platform.archived_tracks import (
    ArchiveStoreError,
    ArchiveVersion,
    SealedMember,
)


def _scope(session: Session, project_id: str, experiment_id: str) -> ExperimentRow:
    experiment = session.get(ExperimentRow, experiment_id)
    project = session.get(ProjectRow, project_id)
    if (
        project is None
        or project.archived_at is not None
        or experiment is None
        or experiment.project_id != project_id
    ):
        raise ArchiveStoreError("archive_version_not_found")
    return experiment


def _scoped_row(
    session: Session, project_id: str, experiment_id: str, archive_version_id: str
) -> Row:
    _scope(session, project_id, experiment_id)
    row = session.get(Row, archive_version_id)
    if (
        row is None
        or row.project_id != project_id
        or row.experiment_id != experiment_id
    ):
        raise ArchiveStoreError("archive_version_not_found")
    return row


def _version(row: Row, *, include_members: bool = True) -> ArchiveVersion:
    members = {}
    if include_members:
        try:
            members = {
                role: SealedMember(
                    props["size_bytes"],
                    props["sha256"],
                    tuple(props["chunks"]),
                    tuple(props["identity"]),
                )
                for role, props in row.members.items()
            }
            if set(members) != (
                {"data", "index"} if row.track_type == "longrange" else {"data"}
            ):
                raise ArchiveStoreError("archive_integrity")
        except (KeyError, TypeError, ValueError, AttributeError):
            raise ArchiveStoreError("archive_integrity") from None
    return ArchiveVersion(
        row.project_id,
        row.experiment_id,
        row.experiment_revision_id,
        row.archive_version_id,
        row.track_key,
        row.version_number,
        row.display_name,
        row.assembly,
        row.track_type,
        dict(row.source),
        members,
    )


def _rebind_members(
    row: Row, previous: ArchiveVersion, members: dict[str, SealedMember]
) -> ArchiveVersion:
    if previous.members.keys() != members.keys() or any(
        not member.same_content(members[role])
        for role, member in previous.members.items()
    ):
        raise ArchiveStoreError("archive_version_conflict")
    if previous.members == members:
        return previous
    row.members = {role: asdict(member) for role, member in members.items()}
    return _version(row)


def _source_matches(
    session: Session,
    project: str,
    exp_revision: ExperimentRevisionRow,
    source: dict[str, str],
    members: dict[str, SealedMember],
) -> None:
    kind = source.get("kind")
    if kind not in {"registered", "run"}:
        raise ArchiveStoreError("archive_invalid")
    fields = (
        ("result_id", "result_revision_id")
        if kind == "registered"
        else (
            "artifact_id",
            "generation",
            "revision",
        )
    )
    allowed = {"kind", *fields}
    if kind == "run":
        allowed.add("run_id")
    if "index" in members:
        allowed.update("index_" + field for field in fields)
    if set(source) != allowed or any(not isinstance(v, str) for v in source.values()):
        raise ArchiveStoreError("archive_invalid")
    if "index" in members and (
        members["index"].sha256 == members["data"].sha256
        or all(source[field] == source["index_" + field] for field in fields)
    ):
        raise ArchiveStoreError("archive_invalid")
    if kind == "run":
        run_project = session.get(RunProjectBindingRow, source["run_id"])
        run_experiment = session.get(RunExperimentBindingRow, source["run_id"])
        if (
            run_project is None
            or run_project.project_id != project
            or run_experiment is None
            or run_experiment.experiment_revision_id
            != exp_revision.experiment_revision_id
            or run_experiment.payload_digest != exp_revision.payload_digest
        ):
            raise ArchiveStoreError("archive_version_not_found")
    for role, member in members.items():
        prefix = "" if role == "data" else "index_"
        if kind == "registered":
            result = session.get(RegisteredResultRow, source[prefix + "result_id"])
            rev = session.get(
                RegisteredResultRevisionRow, source[prefix + "result_revision_id"]
            )
            if (
                result is None
                or result.project_id != project
                or rev is None
                or rev.result_id != result.result_id
            ):
                raise ArchiveStoreError("archive_version_not_found")
            if (
                member.size_bytes != rev.size_bytes
                or member.sha256 != rev.content_sha256
            ):
                raise ArchiveStoreError("archive_integrity")
        else:
            publication = session.get(
                ArtifactPublicationRow,
                {
                    "run_id": source["run_id"],
                    "artifact_id": source[prefix + "artifact_id"],
                    "artifact_generation": source[prefix + "generation"],
                },
            )
            if (
                publication is None
                or publication.artifact_revision != source[prefix + "revision"]
            ):
                raise ArchiveStoreError("archive_version_not_found")


class ArchiveRepository:
    def __init__(self, session_factory: sessionmaker[Session]):
        self._sessions = session_factory

    def list_versions(
        self,
        project_id: str,
        experiment_id: str,
        *,
        limit: int = 50,
        offset: int = 0,
    ) -> tuple[list[ArchiveVersion], int | None]:
        if (
            type(limit) is not int
            or type(offset) is not int
            or not 1 <= limit <= 100
            or offset < 0
        ):
            raise ArchiveStoreError("archive_invalid")
        with self._sessions() as session:
            _scope(session, project_id, experiment_id)
            statement = (
                select(Row)
                .options(defer(Row.members))
                .where(Row.project_id == project_id, Row.experiment_id == experiment_id)
                .order_by(Row.track_key, Row.version_number.desc())
                .limit(limit + 1)
                .offset(offset)
            )
            rows = session.scalars(statement).all()
            return (
                [_version(row, include_members=False) for row in rows[:limit]],
                offset + limit if len(rows) > limit else None,
            )

    def get_version(
        self, project_id: str, experiment_id: str, archive_version_id: str
    ) -> ArchiveVersion:
        with self._sessions() as session:
            return _version(
                _scoped_row(session, project_id, experiment_id, archive_version_id)
            )

    def rebind_restored(
        self, original: ArchiveVersion, verified: dict[str, SealedMember]
    ) -> ArchiveVersion:
        """Commit the verified identities together, only on the requested row."""
        with self._sessions() as session, session.begin():
            session.connection().exec_driver_sql("BEGIN IMMEDIATE")
            row = _scoped_row(
                session,
                original.project_id,
                original.experiment_id,
                original.archive_version_id,
            )
            current = _version(row)
            if current.members != original.members and current.members != verified:
                raise ArchiveStoreError("archive_version_conflict")
            return _rebind_members(row, current, verified)

    def publish(
        self,
        *,
        project_id: str,
        experiment_id: str,
        experiment_revision_id: str,
        track_key: str,
        display_name: str,
        assembly: str,
        track_type: str,
        source: dict[str, str],
        members: dict[str, SealedMember],
    ) -> ArchiveVersion:
        expected = {"data", "index"} if track_type == "longrange" else {"data"}
        if (
            not isinstance(members, dict)
            or not isinstance(source, dict)
            or not all(
                isinstance(key, str) and isinstance(value, str)
                for key, value in source.items()
            )
            or set(members) != expected
            or track_type not in {"bigwig", "bed", "longrange"}
            or assembly not in {"hg19", "hg38", "mm9", "mm10", "mm39"}
            or any(not isinstance(member, SealedMember) for member in members.values())
            or any(
                not isinstance(name, str)
                or not name.strip()
                or len(name) > 255
                or any(c in "/\\" or not c.isprintable() for c in name)
                for name in (track_key, display_name)
            )
        ):
            raise ArchiveStoreError("archive_invalid")
        identity = sha256(
            json.dumps(source, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        with self._sessions() as session, session.begin():
            # SQLite owns numbering across all operator processes, including retries.
            session.connection().exec_driver_sql("BEGIN IMMEDIATE")
            experiment = _scope(session, project_id, experiment_id)
            revision = session.get(ExperimentRevisionRow, experiment_revision_id)
            if (
                revision is None
                or revision.experiment_id != experiment_id
                or revision.project_id != project_id
                or experiment.archived_at is not None
            ):
                raise ArchiveStoreError("archive_version_not_found")
            _source_matches(session, project_id, revision, source, members)
            series = (Row.experiment_id == experiment_id, Row.track_key == track_key)
            existing = session.scalar(
                select(Row).where(*series, Row.source_identity == identity)
            )
            if existing is not None:
                previous = _version(existing)
                if (
                    previous.experiment_revision_id != experiment_revision_id
                    or previous.assembly != assembly
                    or previous.track_type != track_type
                    or previous.display_name != display_name
                ):
                    raise ArchiveStoreError("archive_version_conflict")
                return _rebind_members(existing, previous, members)
            number = (
                session.scalar(select(func.max(Row.version_number)).where(*series)) or 0
            )
            row = Row(
                archive_version_id="aver_" + uuid4().hex,
                project_id=project_id,
                experiment_id=experiment_id,
                experiment_revision_id=experiment_revision_id,
                track_key=track_key,
                version_number=number + 1,
                display_name=display_name,
                assembly=assembly,
                track_type=track_type,
                source_identity=identity,
                source=source,
                members={key: asdict(value) for key, value in members.items()},
            )
            session.add(row)
            session.flush()
            return _version(row)
