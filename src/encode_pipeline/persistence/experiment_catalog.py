"""SQLite experiment catalog repository for Slice A."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from datetime import datetime, timezone
import json
from typing import Any

from sqlalchemy import and_, or_, select, text
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session, aliased, sessionmaker

from encode_pipeline.persistence.models import (
    ArtifactPublicationRow,
    CatalogTermRevisionRow,
    CatalogTermRow,
    CatalogWriteReceiptRow,
    ExperimentCatalogProjectionRow,
    ExperimentRevisionHeadRow,
    ExperimentRevisionRow,
    ExperimentRevisionSampleRow,
    ExperimentRow,
    ProjectRow,
    RunExperimentBindingRow,
    SnapshotExperimentBindingRow,
    RunProjectBindingRow,
    RunResultStateRow,
    RunRow,
    RunSampleRow,
    SampleCatalogProjectionRow,
    SampleCatalogValueRow,
    SampleRevisionCatalogTrustRow,
    SampleRevisionHeadRow,
    SampleRevisionRow,
    SampleRow,
    UserAccountRow,
)
from encode_pipeline.platform.data_registry import (
    LEGACY_PROJECT_ID,
    SampleRevision,
)
from encode_pipeline.platform.experiment_catalog import (
    EXPERIMENT_FIELDS,
    EXPERIMENT_REVISION_PAYLOAD_DIGEST_SCHEME,
    SAMPLE_FIELDS,
    SAMPLE_MULTI_FIELDS,
    SAMPLE_SINGLE_FIELDS,
    CatalogCursorPosition,
    CatalogError,
    CatalogFilters,
    CatalogWrite,
    ExperimentPayload,
    ExperimentSamplePin,
    PinnedTerm,
    SampleCatalogBlock,
    SingleSlot,
    TermLookup,
    TermPin,
    administrator_required,
    bind_catalog_target,
    matches_committed_catalog_request,
    block_from_projection,
    build_experiment_revision_payload_digest,
    canonical_alias_json,
    catalog_already_archived,
    catalog_cursor_not_found,
    catalog_experiment_archived,
    catalog_not_found,
    catalog_payload_invalid,
    publication_experiment_trace,
    catalog_revision_conflict,
    catalog_idempotency_conflict,
    collect_experiment_term_revision_ids,
    collect_term_revision_ids,
    decode_catalog_cursor,
    encode_catalog_cursor,
    interpret_experiment_fields,
    interpret_sample_catalog,
    interpret_stored_experiment_payload,
    _single_object,
    projection_rows_from_block,
    render_sample_material,
    render_single,
    sample_catalog_object,
    validate_aliases,
    validate_display_name,
    validate_optional_external,
    validate_stable_key,
)
from encode_pipeline.platform.data_registry import (
    canonical_sample_revision_payload,
    build_sample_revision_payload_digest,
    SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
)


class SessionTermLookup(TermLookup):
    """Resolve term revisions inside the caller's transaction."""

    def __init__(self, session: Session) -> None:
        self._session = session
        self._cache: dict[str, TermPin | None] = {}

    def __call__(self, term_revision_id: str) -> TermPin | None:
        if term_revision_id in self._cache:
            return self._cache[term_revision_id]
        row = self._session.get(CatalogTermRevisionRow, term_revision_id)
        if row is None:
            self._cache[term_revision_id] = None
            return None
        term = self._session.get(CatalogTermRow, row.term_id)
        if term is None:
            self._cache[term_revision_id] = None
            return None
        pin = TermPin(
            term_id=term.term_id,
            term_revision_id=row.term_revision_id,
            term_type=term.term_type,
            display_name=row.display_name,
            archived=term.archived_at is not None,
        )
        self._cache[term_revision_id] = pin
        return pin


def sync_sample_revision_state(
    session: Session,
    revision: SampleRevision,
    *,
    allow_catalog: bool,
    catalog_trusted: bool | None = None,
) -> None:
    """Point the sample head at this revision and refresh its projection."""

    attributes = json.loads(revision.canonical_payload)["attributes"]
    if not isinstance(attributes, dict):
        raise catalog_payload_invalid()
    if "catalog" in attributes and not allow_catalog:
        from encode_pipeline.platform.experiment_catalog import (
            CatalogBlockRequiresCatalogCommand,
        )

        raise CatalogBlockRequiresCatalogCommand()
    head = session.get(SampleRevisionHeadRow, revision.sample_id)
    if head is None:
        session.add(
            SampleRevisionHeadRow(
                sample_id=revision.sample_id,
                project_id=revision.project_id,
                sample_revision_id=revision.sample_revision_id,
                revision_number=revision.revision_number,
            )
        )
    else:
        head.sample_revision_id = revision.sample_revision_id
        head.revision_number = revision.revision_number
        head.project_id = revision.project_id
    if catalog_trusted is None:
        catalog_trusted = _previous_catalog_trusted(session, revision)
    _record_catalog_trust(session, revision.sample_revision_id, catalog_trusted)
    replace_sample_projection(
        session,
        revision,
        catalog_trusted=catalog_trusted,
    )
    session.flush()


def replace_sample_projection(
    session: Session,
    revision: SampleRevision,
    *,
    catalog_trusted: bool,
) -> None:
    """Project a sample revision without promoting an untrusted catalog key."""

    payload = revision.canonical_payload
    attributes = json.loads(payload)["attributes"]
    block = None
    if catalog_trusted and isinstance(attributes, dict) and "catalog" in attributes:
        block = interpret_sample_catalog(
            attributes["catalog"],
            SessionTermLookup(session),
            for_write=False,
        )
    _write_sample_projection(
        session,
        revision,
        block,
        catalog_trusted=catalog_trusted,
    )
    if session.get(SampleRevisionRow, revision.sample_revision_id) is not None:
        stored = session.get(SampleRevisionRow, revision.sample_revision_id)
        if stored is not None and stored.canonical_payload != payload:
            raise catalog_payload_invalid()


def _previous_catalog_trusted(session: Session, revision: SampleRevision) -> bool:
    """A copied block stays classified only when the previous write was trusted."""

    if revision.revision_number <= 1:
        return False
    previous_id = session.scalar(
        select(SampleRevisionRow.sample_revision_id).where(
            SampleRevisionRow.sample_id == revision.sample_id,
            SampleRevisionRow.revision_number == revision.revision_number - 1,
        )
    )
    if previous_id is None:
        return False
    return _catalog_trust(session, previous_id)


def _catalog_trust(session: Session, sample_revision_id: str) -> bool:
    """Read the trust bit stored on the revision, not on its projection."""

    row = session.get(SampleRevisionCatalogTrustRow, sample_revision_id)
    return row is not None and bool(row.catalog_trusted)


def _record_catalog_trust(
    session: Session,
    sample_revision_id: str,
    trusted: bool,
) -> None:
    """Keep the classification source after the projection tables are removed."""

    session.add(
        SampleRevisionCatalogTrustRow(
            sample_revision_id=sample_revision_id,
            catalog_trusted=1 if trusted else 0,
        )
    )


class SqlAlchemyExperimentCatalogRepository:
    """Transactional catalog store backed by the platform SQLite database."""

    def __init__(self, session_factory: sessionmaker[Session]) -> None:
        if not isinstance(session_factory, sessionmaker):
            raise ValueError("session_factory must be a SQLAlchemy sessionmaker")
        self._session_factory = session_factory

    def create_term(
        self,
        write: CatalogWrite,
        *,
        term_id: str,
        term_revision_id: str,
        term_type: str,
        display_name: str,
        namespace: str | None,
        accession: str | None,
        aliases: tuple[str, ...],
    ) -> dict[str, Any]:
        def apply(session: Session) -> dict[str, Any]:
            if term_type not in {
                "species",
                "tissue",
                "cell_type",
                "cell_line",
                "treatment",
                "tag",
                "assay",
                "target",
                "batch",
                "replicate_role",
            }:
                raise catalog_payload_invalid()
            name = validate_display_name(display_name)
            session.add(
                CatalogTermRow(
                    term_id=term_id,
                    term_type=term_type,
                    created_at=write.now,
                    archived_at=None,
                )
            )
            session.flush()
            session.add(
                _term_revision_row(
                    term_id=term_id,
                    term_revision_id=term_revision_id,
                    revision_number=1,
                    display_name=name,
                    namespace=validate_optional_external(namespace),
                    accession=validate_optional_external(accession),
                    aliases=aliases,
                    change_kind="initial",
                    created_at=write.now,
                )
            )
            session.flush()
            _add_receipt(
                session,
                write,
                subject_id=term_id,
                result_revision_id=term_revision_id,
                result_archived_at=None,
            )
            return _term_body(session, term_id, term_revision_id)

        write = bind_catalog_target(write, "vocabulary")
        return self._write(write, apply)

    def revise_term(
        self,
        write: CatalogWrite,
        *,
        term_id: str,
        term_revision_id: str,
        expected_revision_number: int,
        display_name: str,
        namespace: str | None,
        accession: str | None,
        aliases: tuple[str, ...],
    ) -> dict[str, Any]:
        def apply(session: Session) -> dict[str, Any]:
            term = session.get(CatalogTermRow, term_id)
            if term is None:
                raise catalog_not_found()
            current = _term_head_number(session, term_id)
            _expect_revision(current, expected_revision_number)
            if term.archived_at is not None:
                raise catalog_term_archived_error()
            name = validate_display_name(display_name)
            session.add(
                _term_revision_row(
                    term_id=term_id,
                    term_revision_id=term_revision_id,
                    revision_number=current + 1,
                    display_name=name,
                    namespace=validate_optional_external(namespace),
                    accession=validate_optional_external(accession),
                    aliases=aliases,
                    change_kind="editorial",
                    created_at=write.now,
                )
            )
            session.flush()
            _add_receipt(
                session,
                write,
                subject_id=term_id,
                result_revision_id=term_revision_id,
                result_archived_at=None,
            )
            return _term_body(session, term_id, term_revision_id)

        write = bind_catalog_target(write, f"term:{term_id}")
        return self._write(write, apply)

    def archive_term(
        self,
        write: CatalogWrite,
        *,
        term_id: str,
        expected_revision_number: int,
    ) -> dict[str, Any]:
        def apply(session: Session) -> dict[str, Any]:
            term = session.get(CatalogTermRow, term_id)
            if term is None:
                raise catalog_not_found()
            current = _term_head_number(session, term_id)
            _expect_revision(current, expected_revision_number)
            if term.archived_at is not None:
                raise catalog_already_archived()
            term.archived_at = write.now
            session.flush()
            _add_receipt(
                session,
                write,
                subject_id=term_id,
                result_revision_id=None,
                result_archived_at=write.now,
            )
            return _term_body(session, term_id, _term_head_id(session, term_id))

        write = bind_catalog_target(write, f"term:{term_id}")
        return self._write(write, apply)

    def create_experiment(
        self,
        write: CatalogWrite,
        *,
        project_id: str,
        experiment_id: str,
        experiment_revision_id: str,
        stable_key: str,
        display_name: str,
        assay: object,
        target: object,
        batch: object,
        samples: object,
    ) -> dict[str, Any]:
        def apply(session: Session) -> dict[str, Any]:
            _require_active_project(session, project_id)
            key = validate_stable_key(stable_key)
            name = validate_display_name(display_name)
            payload = _experiment_payload(
                session,
                project_id=project_id,
                display_name=name,
                assay=assay,
                target=target,
                batch=batch,
                samples=samples,
            )
            canonical = payload.canonical_json()
            session.add(
                ExperimentRow(
                    experiment_id=experiment_id,
                    project_id=project_id,
                    stable_key=key,
                    created_at=write.now,
                    archived_at=None,
                )
            )
            session.flush()
            _insert_experiment_revision(
                session,
                experiment_id=experiment_id,
                project_id=project_id,
                experiment_revision_id=experiment_revision_id,
                revision_number=1,
                canonical=canonical,
                payload=payload,
                created_at=write.now,
                identity_created_at=write.now,
            )
            _add_receipt(
                session,
                write,
                subject_id=experiment_id,
                result_revision_id=experiment_revision_id,
                result_archived_at=None,
            )
            return self._experiment_detail(
                session,
                experiment_revision_id,
                include_candidates=True,
            )

        write = bind_catalog_target(write, f"project:{project_id}")
        try:
            return self._write(write, apply)
        except IntegrityError as exc:
            raise catalog_revision_conflict() from exc

    def revise_experiment(
        self,
        write: CatalogWrite,
        *,
        experiment_id: str,
        experiment_revision_id: str,
        expected_revision_number: int,
        display_name: str,
        assay: object,
        target: object,
        batch: object,
        samples: object,
    ) -> dict[str, Any]:
        def apply(session: Session) -> dict[str, Any]:
            experiment = session.get(ExperimentRow, experiment_id)
            if experiment is None:
                raise catalog_not_found()
            _require_active_project(session, experiment.project_id)
            head = session.get(ExperimentRevisionHeadRow, experiment_id)
            if head is None:
                raise catalog_not_found()
            _expect_revision(head.revision_number, expected_revision_number)
            if experiment.archived_at is not None:
                raise catalog_experiment_archived()
            payload = _experiment_payload(
                session,
                project_id=experiment.project_id,
                display_name=validate_display_name(display_name),
                assay=assay,
                target=target,
                batch=batch,
                samples=samples,
            )
            canonical = payload.canonical_json()
            _insert_experiment_revision(
                session,
                experiment_id=experiment_id,
                project_id=experiment.project_id,
                experiment_revision_id=experiment_revision_id,
                revision_number=head.revision_number + 1,
                canonical=canonical,
                payload=payload,
                created_at=write.now,
                identity_created_at=experiment.created_at,
            )
            _add_receipt(
                session,
                write,
                subject_id=experiment_id,
                result_revision_id=experiment_revision_id,
                result_archived_at=None,
            )
            return self._experiment_detail(
                session,
                experiment_revision_id,
                include_candidates=True,
            )

        write = bind_catalog_target(write, f"experiment:{experiment_id}")
        return self._write(write, apply)

    def archive_experiment(
        self,
        write: CatalogWrite,
        *,
        experiment_id: str,
        expected_revision_number: int,
    ) -> dict[str, Any]:
        def apply(session: Session) -> dict[str, Any]:
            experiment = session.get(ExperimentRow, experiment_id)
            if experiment is None:
                raise catalog_not_found()
            _require_active_project(session, experiment.project_id)
            head = session.get(ExperimentRevisionHeadRow, experiment_id)
            if head is None:
                raise catalog_not_found()
            _expect_revision(head.revision_number, expected_revision_number)
            if experiment.archived_at is not None:
                raise catalog_already_archived()
            experiment.archived_at = write.now
            session.flush()
            _add_receipt(
                session,
                write,
                subject_id=experiment_id,
                result_revision_id=None,
                result_archived_at=write.now,
            )
            return {
                "experiment_id": experiment.experiment_id,
                "archived_at": _stamp(write.now),
                "revision_number": head.revision_number,
            }

        write = bind_catalog_target(write, f"experiment:{experiment_id}")
        return self._write(write, apply)

    def revise_sample_catalog(
        self,
        write: CatalogWrite,
        *,
        sample_id: str,
        sample_revision_id: str,
        expected_revision_number: int,
        display_name: str | None,
        catalog: object,
    ) -> dict[str, Any]:
        def apply(session: Session) -> dict[str, Any]:
            sample = session.get(SampleRow, sample_id)
            if sample is None:
                raise catalog_not_found()
            _require_active_project(session, sample.project_id)
            head = session.get(SampleRevisionHeadRow, sample_id)
            if head is None:
                raise catalog_not_found()
            _expect_revision(head.revision_number, expected_revision_number)
            previous = session.get(SampleRevisionRow, head.sample_revision_id)
            if previous is None:
                raise catalog_not_found()
            decoded = json.loads(previous.canonical_payload)
            attributes = {
                key: value
                for key, value in decoded["attributes"].items()
                if key != "catalog"
            }
            block = interpret_sample_catalog(
                catalog,
                SessionTermLookup(session),
                for_write=True,
            )
            if block is not None:
                attributes["catalog"] = sample_catalog_object(block)
            name = (
                validate_display_name(display_name)
                if display_name is not None
                else decoded["display_name"]
            )
            canonical = canonical_sample_revision_payload(
                display_name=name,
                attributes=attributes,
            )
            revision = SampleRevision(
                sample_revision_id=sample_revision_id,
                project_id=sample.project_id,
                sample_id=sample.sample_id,
                revision_number=previous.revision_number + 1,
                canonical_payload=canonical,
                payload_digest_scheme=SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
                payload_digest=build_sample_revision_payload_digest(canonical),
                created_at=write.now,
            )
            session.add(_sample_revision_row(revision))
            session.flush()
            sync_sample_revision_state(
                session,
                revision,
                allow_catalog=True,
                catalog_trusted=True,
            )
            _add_receipt(
                session,
                write,
                subject_id=sample_id,
                result_revision_id=sample_revision_id,
                result_archived_at=None,
            )
            return self._sample_revision_body(
                session,
                sample_revision_id,
                include_candidates=True,
            )

        write = bind_catalog_target(write, f"sample:{sample_id}")
        return self._write(write, apply)

    def list_experiments(
        self,
        *,
        filters: CatalogFilters,
        cursor: str | None,
        limit: int,
    ) -> dict[str, Any]:
        position = decode_catalog_cursor(
            cursor,
            kind="experiments",
            filters=filters,
        )
        with self._session_factory() as session:
            _require_active_project(session, _required_project(filters))
            statement = _experiment_statement(filters)
            if position is not None:
                _require_cursor_row(
                    session,
                    statement,
                    ExperimentRow.experiment_id,
                    ExperimentRow.created_at,
                    position,
                )
                statement = statement.where(
                    _after(
                        ExperimentRow.created_at,
                        ExperimentRow.experiment_id,
                        position,
                    )
                )
            rows = session.execute(
                statement.order_by(
                    ExperimentRow.created_at,
                    ExperimentRow.experiment_id,
                ).limit(limit + 1)
            ).all()
            selected = rows[:limit]
            experiments = [
                self._experiment_summary(
                    session,
                    row.experiment_revision_id,
                    revision_view=filters.revision_view,
                    include_candidates=False,
                )
                for row in selected
            ]
            return {
                "ok": True,
                "project_id": filters.project_id,
                "revision_view": filters.revision_view,
                "experiments": experiments,
                "next_cursor": _next_cursor(
                    rows,
                    limit,
                    kind="experiments",
                    filters=filters,
                    created_at=lambda row: row.created_at,
                    object_id=lambda row: row.experiment_id,
                ),
            }

    def list_samples(
        self,
        *,
        filters: CatalogFilters,
        cursor: str | None,
        limit: int,
    ) -> dict[str, Any]:
        position = decode_catalog_cursor(cursor, kind="samples", filters=filters)
        with self._session_factory() as session:
            _require_active_project(session, _required_project(filters))
            statement = _sample_statement(filters)
            if position is not None:
                _require_cursor_row(
                    session,
                    statement,
                    SampleRow.sample_id,
                    SampleRow.created_at,
                    position,
                )
                statement = statement.where(
                    _after(SampleRow.created_at, SampleRow.sample_id, position)
                )
            rows = session.execute(
                statement.order_by(SampleRow.created_at, SampleRow.sample_id).limit(
                    limit + 1
                )
            ).all()
            samples = [
                self._sample_summary(session, row.sample_revision_id)
                for row in rows[:limit]
            ]
            return {
                "project_id": filters.project_id,
                "samples": samples,
                "next_cursor": _next_cursor(
                    rows,
                    limit,
                    kind="samples",
                    filters=filters,
                    created_at=lambda row: row.created_at,
                    object_id=lambda row: row.sample_id,
                ),
            }

    def list_terms(
        self,
        *,
        filters: CatalogFilters,
        cursor: str | None,
        limit: int,
    ) -> dict[str, Any]:
        position = decode_catalog_cursor(cursor, kind="terms", filters=filters)
        with self._session_factory() as session:
            statement = _term_statement(filters)
            if position is not None:
                _require_cursor_row(
                    session,
                    statement,
                    CatalogTermRow.term_id,
                    CatalogTermRow.created_at,
                    position,
                )
                statement = statement.where(
                    _after(
                        CatalogTermRow.created_at,
                        CatalogTermRow.term_id,
                        position,
                    )
                )
            rows = session.execute(
                statement.order_by(
                    CatalogTermRow.created_at,
                    CatalogTermRow.term_id,
                ).limit(limit + 1)
            ).all()
            terms = []
            for row in rows[:limit]:
                revision_id = _term_head_id(session, row.term_id)
                body = _term_body(session, row.term_id, revision_id)
                revision = body["revision"]
                terms.append(
                    {
                        "term_id": body["term_id"],
                        "term_type": body["term_type"],
                        "term_revision_id": revision["term_revision_id"],
                        "revision_number": revision["revision_number"],
                        "display_name": revision["display_name"],
                        "archived_at": body["archived_at"],
                    }
                )
            return {
                "terms": terms,
                "next_cursor": _next_cursor(
                    rows,
                    limit,
                    kind="terms",
                    filters=filters,
                    created_at=lambda row: row.created_at,
                    object_id=lambda row: row.term_id,
                ),
            }

    def get_term(self, term_id: str) -> dict[str, Any]:
        with self._session_factory() as session:
            if session.get(CatalogTermRow, term_id) is None:
                raise catalog_not_found()
            return _term_body(session, term_id, _term_head_id(session, term_id))

    def list_sample_revisions(
        self,
        *,
        sample_id: str,
        filters: CatalogFilters,
        cursor: str | None,
        limit: int,
    ) -> dict[str, Any]:
        return self._list_owned_revisions(
            kind="sample_revisions",
            owner_id=sample_id,
            filters=filters,
            cursor=cursor,
            limit=limit,
            model=SampleRevisionRow,
            owner_column=SampleRevisionRow.sample_id,
            id_column=SampleRevisionRow.sample_revision_id,
            created_column=SampleRevisionRow.created_at,
            render=lambda session, revision_id: self._sample_revision_body(
                session,
                revision_id,
                include_candidates=False,
            ),
        )

    def get_sample_revision(
        self, sample_id: str, sample_revision_id: str
    ) -> dict[str, Any]:
        with self._session_factory() as session:
            row = session.get(SampleRevisionRow, sample_revision_id)
            if row is None or row.sample_id != sample_id:
                raise catalog_not_found()
            _require_active_project(session, row.project_id)
            return self._sample_revision_body(
                session,
                sample_revision_id,
                include_candidates=True,
            )

    def list_experiment_revisions(
        self,
        *,
        experiment_id: str,
        filters: CatalogFilters,
        cursor: str | None,
        limit: int,
    ) -> dict[str, Any]:
        return self._list_owned_revisions(
            kind="experiment_revisions",
            owner_id=experiment_id,
            filters=filters,
            cursor=cursor,
            limit=limit,
            model=ExperimentRevisionRow,
            owner_column=ExperimentRevisionRow.experiment_id,
            id_column=ExperimentRevisionRow.experiment_revision_id,
            created_column=ExperimentRevisionRow.created_at,
            render=lambda session, revision_id: self._experiment_summary(
                session,
                revision_id,
                revision_view=filters.revision_view,
                include_candidates=False,
            ),
        )

    def get_experiment_revision(
        self,
        experiment_id: str,
        experiment_revision_id: str,
    ) -> dict[str, Any]:
        with self._session_factory() as session:
            row = session.get(ExperimentRevisionRow, experiment_revision_id)
            if row is None or row.experiment_id != experiment_id:
                raise catalog_not_found()
            _require_active_project(session, row.project_id)
            return self._experiment_detail(
                session,
                experiment_revision_id,
                include_candidates=True,
            )

    def classification_context(
        self, experiment_id: str, write: CatalogWrite
    ) -> dict[str, Any]:
        """Current experiment head and the active curation versions."""

        from encode_pipeline.persistence.catalog_proposal_store import (
            read_active_versions,
        )

        with self._session_factory() as session:
            _checked_write_receipt(session, write)
            experiment = session.get(ExperimentRow, experiment_id)
            if experiment is None:
                raise catalog_not_found()
            _require_active_project(session, experiment.project_id)
            head = session.get(ExperimentRevisionHeadRow, experiment_id)
            if head is None:
                raise catalog_not_found()
            revision = session.get(ExperimentRevisionRow, head.experiment_revision_id)
            if revision is None:
                raise catalog_not_found()
            payload = interpret_stored_experiment_payload(revision.canonical_payload)
            active = read_active_versions(session)
            return {
                "project_id": experiment.project_id,
                "experiment_id": experiment.experiment_id,
                "experiment_revision_id": revision.experiment_revision_id,
                "revision_number": revision.revision_number,
                "archived": experiment.archived_at is not None,
                "target": _single_object(payload.target),
                "rule_version": active.rule_version,
                "vocabulary_version": active.vocabulary_version,
                "projection_version": active.projection_version,
            }

    def project_target_terms(self, project_id: str, limit: int) -> list[dict[str, Any]]:
        """Only offer a vocabulary revision recorded on an active project head."""

        projection = ExperimentCatalogProjectionRow
        head = ExperimentRevisionHeadRow
        used_revisions = (
            select(projection.target_term_revision_id)
            .join(
                head, head.experiment_revision_id == projection.experiment_revision_id
            )
            .join(ExperimentRow, ExperimentRow.experiment_id == head.experiment_id)
            .where(
                head.project_id == project_id,
                ExperimentRow.archived_at.is_(None),
                projection.target_status == "known",
            )
        )
        latest_revision = (
            select(CatalogTermRevisionRow.term_revision_id)
            .where(CatalogTermRevisionRow.term_id == CatalogTermRow.term_id)
            .order_by(CatalogTermRevisionRow.revision_number.desc())
            .limit(1)
            .correlate(CatalogTermRow)
            .scalar_subquery()
        )
        with self._session_factory() as session:
            # Read labels from the same snapshot that proved the revision was used.
            rows = session.execute(
                select(
                    CatalogTermRow.term_id,
                    CatalogTermRevisionRow.term_revision_id,
                    CatalogTermRevisionRow.display_name,
                    CatalogTermRevisionRow.aliases_json,
                )
                .join(
                    CatalogTermRevisionRow,
                    CatalogTermRevisionRow.term_id == CatalogTermRow.term_id,
                )
                .where(
                    CatalogTermRow.term_type == "target",
                    CatalogTermRow.archived_at.is_(None),
                    CatalogTermRevisionRow.term_revision_id.in_(used_revisions),
                    CatalogTermRevisionRow.term_revision_id == latest_revision,
                )
                .order_by(CatalogTermRow.created_at, CatalogTermRow.term_id)
                .limit(limit)
            )
            return [
                {
                    "term_id": term_id,
                    "term_revision_id": revision_id,
                    "display_name": name,
                    "aliases": json.loads(aliases),
                }
                for term_id, revision_id, name, aliases in rows
            ]

    def lookup_target_term(self, term_revision_id: str) -> dict[str, Any] | None:
        with self._session_factory() as session:
            revision = session.get(CatalogTermRevisionRow, term_revision_id)
            if revision is None:
                return None
            term = session.get(CatalogTermRow, revision.term_id)
            if term is None or term.term_type != "target":
                return None
            record = _target_term_record(session, term.term_id)
            record["requested_is_head"] = record["term_revision_id"] == term_revision_id
            return record

    def list_sample_run_results(
        self,
        *,
        project_id: str,
        sample_revision_id: str,
        filters: CatalogFilters,
        cursor: str | None,
        limit: int,
    ) -> dict[str, Any]:
        position = decode_catalog_cursor(
            cursor,
            kind="sample_runs",
            filters=filters,
        )
        with self._session_factory() as session:
            _require_active_project(session, project_id)
            revision = session.get(SampleRevisionRow, sample_revision_id)
            if revision is None or revision.project_id != project_id:
                raise catalog_not_found()
            results = _sample_trace_results(session, project_id, sample_revision_id)
            start = 0
            if position is not None:
                start = _trace_index(results, position)
            page = results[start : start + limit]
            next_cursor = None
            if start + limit < len(results):
                last = page[-1]
                next_cursor = encode_catalog_cursor(
                    kind="sample_runs",
                    filters=filters,
                    position=CatalogCursorPosition(
                        created_at=_parse_trace_stamp(last),
                        object_id=_trace_object_id(last),
                    ),
                )
            return {
                "ok": True,
                "project_id": project_id,
                "sample_revision_id": sample_revision_id,
                "results": page,
                "next_cursor": next_cursor,
            }

    def list_experiment_run_results(
        self,
        *,
        project_id: str,
        experiment_revision_id: str,
        filters: CatalogFilters,
        cursor: str | None,
        limit: int,
    ) -> dict[str, Any]:
        position = decode_catalog_cursor(
            cursor,
            kind="experiment_runs",
            filters=filters,
        )
        with self._session_factory() as session:
            _require_active_project(session, project_id)
            row = session.get(ExperimentRevisionRow, experiment_revision_id)
            if row is None or row.project_id != project_id:
                raise catalog_not_found()
            snapshot_bound = session.scalar(
                select(SnapshotExperimentBindingRow.snapshot_id)
                .where(
                    SnapshotExperimentBindingRow.experiment_revision_id
                    == experiment_revision_id,
                    SnapshotExperimentBindingRow.payload_digest == row.payload_digest,
                )
                .limit(1)
            )
            trusted_run_ids = list(
                session.scalars(
                    select(RunExperimentBindingRow.run_id).where(
                        RunExperimentBindingRow.experiment_revision_id
                        == experiment_revision_id,
                        RunExperimentBindingRow.payload_digest == row.payload_digest,
                    )
                ).all()
            )
            if snapshot_bound is None and not trusted_run_ids:
                if position is not None:
                    raise catalog_cursor_not_found()
                return {
                    "ok": True,
                    "project_id": project_id,
                    "experiment_revision_id": experiment_revision_id,
                    "trace_status": "unresolved",
                    "reason": "experiment_snapshot_binding_absent",
                    "results": [],
                    "next_cursor": None,
                }
            results = _published_trace_rows(session, project_id, trusted_run_ids)
            start = 0
            if position is not None:
                start = _trace_index(results, position)
            page = results[start : start + limit]
            next_cursor = None
            if start + limit < len(results):
                last = page[-1]
                next_cursor = encode_catalog_cursor(
                    kind="experiment_runs",
                    filters=filters,
                    position=CatalogCursorPosition(
                        created_at=_parse_trace_stamp(last),
                        object_id=_trace_object_id(last),
                    ),
                )
            return {
                "ok": True,
                "project_id": project_id,
                "experiment_revision_id": experiment_revision_id,
                "trace_status": "experiment_snapshot",
                "results": page,
                "next_cursor": next_cursor,
            }

    def get_publication_experiment_trace(
        self,
        *,
        project_id: str,
        run_id: str,
        artifact_id: str,
        artifact_generation: str,
    ) -> dict[str, Any]:
        with self._session_factory() as session:
            _require_active_project(session, project_id)
            matches = [
                row
                for row in _published_trace_rows(session, project_id, [run_id])
                if row["artifact_id"] == artifact_id
                and row["artifact_generation"] == artifact_generation
            ]
            if len(matches) != 1:
                raise catalog_not_found()
            binding = session.get(RunExperimentBindingRow, run_id)
            revision = None
            if binding is not None:
                revision = session.get(
                    ExperimentRevisionRow,
                    binding.experiment_revision_id,
                )
            return publication_experiment_trace(
                project_id=project_id,
                run_id=run_id,
                artifact_id=artifact_id,
                artifact_generation=artifact_generation,
                result=matches[0],
                binding_revision_id=(
                    None if binding is None else binding.experiment_revision_id
                ),
                binding_payload_digest=(
                    None if binding is None else binding.payload_digest
                ),
                revision_project_id=(None if revision is None else revision.project_id),
                revision_payload_digest=(
                    None if revision is None else revision.payload_digest
                ),
            )

    def rebuild_sample_projection(self, sample_revision_id: str) -> None:
        with self._session_factory.begin() as session:
            _begin_write(session)
            row = session.get(SampleRevisionRow, sample_revision_id)
            if row is None:
                raise catalog_not_found()
            before = row.canonical_payload
            revision = _sample_revision_from_row(row)
            trusted = _catalog_trust(session, sample_revision_id)
            replace_sample_projection(
                session,
                revision,
                catalog_trusted=trusted,
            )
            session.flush()
            if row.canonical_payload != before:
                raise catalog_payload_invalid()

    def rebuild_experiment_projection(self, experiment_revision_id: str) -> None:
        with self._session_factory.begin() as session:
            _begin_write(session)
            row = session.get(ExperimentRevisionRow, experiment_revision_id)
            if row is None:
                raise catalog_not_found()
            before = row.canonical_payload
            payload = interpret_stored_experiment_payload(before)
            if payload.canonical_json() != before:
                raise catalog_payload_invalid()
            pins = _pins_from_stored_payload(session, row.project_id, payload)
            session.query(ExperimentRevisionSampleRow).filter(
                ExperimentRevisionSampleRow.experiment_revision_id
                == experiment_revision_id
            ).delete(synchronize_session=False)
            session.query(ExperimentCatalogProjectionRow).filter(
                ExperimentCatalogProjectionRow.experiment_revision_id
                == experiment_revision_id
            ).delete(synchronize_session=False)
            session.flush()
            _write_experiment_children(
                session,
                experiment_id=row.experiment_id,
                project_id=row.project_id,
                experiment_revision_id=experiment_revision_id,
                payload=ExperimentPayload(
                    display_name=payload.display_name,
                    assay=payload.assay,
                    target=payload.target,
                    batch=payload.batch,
                    samples=pins,
                ),
            )
            session.flush()
            if row.canonical_payload != before:
                raise catalog_payload_invalid()

    def explain_experiment_query(self, filters: CatalogFilters) -> str:
        with self._session_factory() as session:
            return _explain(session, _experiment_statement(filters))

    def explain_sample_trace(self, project_id: str, sample_revision_id: str) -> str:
        with self._session_factory() as session:
            lookup = select(RunSampleRow.run_id).where(
                RunSampleRow.sample_revision_id == sample_revision_id,
                RunSampleRow.project_id == project_id,
            )
            run_ids = [row for row in session.scalars(lookup).all()] or ["run-plan"]
            publications = _publication_statement(project_id, run_ids)
            return _explain(session, lookup) + "\n" + _explain(session, publications)

    def register_proposal(
        self,
        write: CatalogWrite,
        *,
        project_id: str,
        source: str,
        items: object,
        rule_version: str | None,
        vocabulary_version: str | None,
        projection_version: str | None,
    ) -> dict[str, Any]:
        from encode_pipeline.persistence.catalog_proposal_store import (
            register_proposal,
        )
        from encode_pipeline.platform.catalog_proposal_adoption import (
            stable_proposal_id,
        )

        proposal_id = stable_proposal_id(
            write.actor_user_id,
            write.idempotency_key,
            project_id,
        )
        write = bind_catalog_target(write, proposal_id)

        def apply(session: Session) -> dict[str, Any]:
            return register_proposal(
                session,
                write,
                project_id=project_id,
                proposal_id=proposal_id,
                source=source,
                items=items,
                rule_version=rule_version,
                vocabulary_version=vocabulary_version,
                projection_version=projection_version,
            )

        return self._write(write, apply)

    def get_proposal(
        self,
        *,
        proposal_id: str,
        include_candidates: bool,
        now: datetime,
    ) -> dict[str, Any]:
        from encode_pipeline.persistence.catalog_proposal_store import get_proposal

        with self._session_factory() as session:
            return get_proposal(
                session,
                proposal_id=proposal_id,
                include_candidates=include_candidates,
                now=now,
            )

    def preview_proposal(
        self,
        *,
        actor_user_id: str,
        proposal_id: str,
        selected_item_ids: object,
        now: datetime,
    ) -> dict[str, Any]:
        from encode_pipeline.persistence.catalog_proposal_store import (
            preview_proposal,
        )

        with self._session_factory() as session:
            _require_administrator(session, actor_user_id)
            return preview_proposal(
                session,
                proposal_id=proposal_id,
                selected_item_ids=selected_item_ids,
                now=now,
            )

    def adopt_proposal(
        self,
        write: CatalogWrite,
        *,
        proposal_id: str,
        selected_item_ids: object,
        expected_revisions: object,
        preview_digest: object,
        experiment_revision_id_factory,
        sample_revision_id_factory,
        receipt_id_factory,
    ) -> dict[str, Any]:
        from encode_pipeline.persistence.catalog_proposal_store import (
            AdoptionConflict,
            adopt_proposal,
            resolve_adoption_conflict,
        )

        write = bind_catalog_target(write, proposal_id)

        def apply(session: Session) -> dict[str, Any]:
            return adopt_proposal(
                session,
                write,
                proposal_id=proposal_id,
                selected_item_ids=selected_item_ids,
                expected_revisions=expected_revisions,
                preview_digest=preview_digest,
                experiment_revision_id_factory=experiment_revision_id_factory,
                sample_revision_id_factory=sample_revision_id_factory,
                receipt_id_factory=receipt_id_factory,
            )

        try:
            return self._write(write, apply)
        except AdoptionConflict:
            with self._session_factory() as session:
                return resolve_adoption_conflict(session, write, proposal_id)

    def list_adoption_receipts(
        self,
        *,
        actor_user_id: str,
        project_id: str,
        proposal_id: str | None,
    ) -> dict[str, Any]:
        from encode_pipeline.persistence.catalog_proposal_store import (
            list_adoption_receipts,
        )

        with self._session_factory() as session:
            _require_administrator(session, actor_user_id)
            return list_adoption_receipts(
                session,
                project_id=project_id,
                proposal_id=proposal_id,
            )

    def set_active_curation_versions(
        self,
        *,
        rule_version: str,
        vocabulary_version: str,
        projection_version: str,
    ) -> None:
        from encode_pipeline.persistence.catalog_proposal_store import (
            write_active_versions,
        )
        from encode_pipeline.platform.catalog_proposal_adoption import (
            ActiveCurationVersions,
            _require_version,
        )

        _require_version(rule_version)
        _require_version(vocabulary_version)
        _require_version(projection_version)
        versions = ActiveCurationVersions(
            rule_version=rule_version,
            vocabulary_version=vocabulary_version,
            projection_version=projection_version,
        )
        with self._session_factory.begin() as session:
            write_active_versions(session, versions)

    def _write(self, write: CatalogWrite, apply) -> dict[str, Any]:
        try:
            with self._session_factory.begin() as session:
                _begin_write(session)
                receipt = _checked_write_receipt(session, write)
                if receipt is not None:
                    return _replay(self, session, receipt)
                return apply(session)
        except IntegrityError as exc:
            raise catalog_revision_conflict() from exc

    def _list_owned_revisions(
        self,
        *,
        kind: str,
        owner_id: str,
        filters: CatalogFilters,
        cursor: str | None,
        limit: int,
        model: type,
        owner_column,
        id_column,
        created_column,
        render,
    ) -> dict[str, Any]:
        del model
        position = decode_catalog_cursor(cursor, kind=kind, filters=filters)
        with self._session_factory() as session:
            if (
                filters.revision_view == "current_sample_catalog"
                and session.get_bind().dialect.name == "sqlite"
            ):
                # Resolve term names and sample heads from one SQLite read snapshot.
                session.connection().exec_driver_sql("BEGIN")
            owner = _owner_project(session, kind, owner_id)
            _require_active_project(session, owner)
            statement = select(id_column, created_column).where(
                owner_column == owner_id
            )
            if position is not None:
                _require_cursor_row(
                    session,
                    statement,
                    id_column,
                    created_column,
                    position,
                )
                statement = statement.where(_after(created_column, id_column, position))
            rows = session.execute(
                statement.order_by(created_column, id_column).limit(limit + 1)
            ).all()
            revisions = [render(session, row[0]) for row in rows[:limit]]
            key = (
                "sample_revisions"
                if kind == "sample_revisions"
                else "experiment_revisions"
            )
            return {
                key: revisions,
                "next_cursor": _next_cursor(
                    rows,
                    limit,
                    kind=kind,
                    filters=filters,
                    created_at=lambda row: row[1],
                    object_id=lambda row: row[0],
                ),
            }

    def _sample_summary(
        self, session: Session, sample_revision_id: str
    ) -> dict[str, Any]:
        body = self._sample_revision_body(
            session,
            sample_revision_id,
            include_candidates=False,
        )
        return body

    def _sample_revision_body(
        self,
        session: Session,
        sample_revision_id: str,
        *,
        include_candidates: bool,
    ) -> dict[str, Any]:
        row = session.get(SampleRevisionRow, sample_revision_id)
        if row is None:
            raise catalog_not_found()
        decoded = json.loads(row.canonical_payload)
        projection = session.get(SampleCatalogProjectionRow, sample_revision_id)
        trusted = projection is not None and bool(projection.catalog_trusted)
        block = None
        if trusted and "catalog" in decoded["attributes"]:
            block = interpret_sample_catalog(
                decoded["attributes"]["catalog"],
                SessionTermLookup(session),
                for_write=False,
            )
        names = _names_for(session, collect_term_revision_ids(block))
        material = render_sample_material(
            block,
            names,
            include_candidates=include_candidates,
        )
        return {
            "sample_id": row.sample_id,
            "sample_revision_id": row.sample_revision_id,
            "revision_number": row.revision_number,
            "display_name": decoded["display_name"],
            "created_at": _stamp(row.created_at),
            **material,
        }

    def _experiment_summary(
        self,
        session: Session,
        experiment_revision_id: str,
        *,
        revision_view: str,
        include_candidates: bool,
    ) -> dict[str, Any]:
        return self._experiment_detail(
            session,
            experiment_revision_id,
            include_candidates=include_candidates,
            revision_view=revision_view,
        )

    def _experiment_detail(
        self,
        session: Session,
        experiment_revision_id: str,
        *,
        include_candidates: bool,
        revision_view: str = "linked",
    ) -> dict[str, Any]:
        row = session.get(ExperimentRevisionRow, experiment_revision_id)
        if row is None:
            raise catalog_not_found()
        payload = interpret_stored_experiment_payload(row.canonical_payload)
        pins = session.scalars(
            select(ExperimentRevisionSampleRow)
            .where(
                ExperimentRevisionSampleRow.experiment_revision_id
                == experiment_revision_id
            )
            .order_by(ExperimentRevisionSampleRow.ordinal)
        ).all()
        names = _names_for(
            session, _experiment_name_ids(session, payload, pins, revision_view)
        )
        samples = [
            _render_experiment_sample(
                session,
                pin,
                names,
                revision_view=revision_view,
                include_candidates=include_candidates,
            )
            for pin in pins
        ]
        return {
            "experiment_id": row.experiment_id,
            "experiment_revision_id": row.experiment_revision_id,
            "revision_number": row.revision_number,
            "display_name": payload.display_name,
            "assay": render_single(
                payload.assay,
                names,
                include_candidates=include_candidates,
            ),
            "target": render_single(
                payload.target,
                names,
                include_candidates=include_candidates,
            ),
            "batch": render_single(
                payload.batch,
                names,
                include_candidates=include_candidates,
            ),
            "samples": samples,
        }


def _render_experiment_sample(
    session: Session,
    pin: ExperimentRevisionSampleRow,
    names: Mapping[str, str],
    *,
    revision_view: str,
    include_candidates: bool,
) -> dict[str, Any]:
    linked_id = pin.sample_revision_id
    current_id = linked_id
    if revision_view == "current_sample_catalog":
        head = session.get(SampleRevisionHeadRow, pin.sample_id)
        if head is not None:
            current_id = head.sample_revision_id
    interpreted_id = (
        current_id if revision_view == "current_sample_catalog" else linked_id
    )
    block = _block_for_revision(session, interpreted_id)
    material = render_sample_material(
        block,
        names,
        include_candidates=False,
    )
    if material.get("catalog_state") == "present":
        material.pop("catalog_state", None)
    role = SingleSlot(
        status=pin.replicate_role_status,
        term=(
            None
            if pin.replicate_role_term_id is None
            or pin.replicate_role_term_revision_id is None
            else PinnedTerm(
                term_id=pin.replicate_role_term_id,
                term_revision_id=pin.replicate_role_term_revision_id,
            )
        ),
    )
    if include_candidates and role.status == "ambiguous":
        stored = interpret_stored_experiment_payload(
            session.get(
                ExperimentRevisionRow, pin.experiment_revision_id
            ).canonical_payload
        )
        for item in stored.samples:
            if item.ordinal == pin.ordinal:
                role = item.replicate_role
                break
    body = {
        "sample_id": pin.sample_id,
        "sample_revision_id": interpreted_id,
        "ordinal": pin.ordinal,
        "replicate_group": pin.replicate_group,
        "replicate_role": render_single(
            role,
            names,
            include_candidates=include_candidates,
        ),
    }
    body.update(material)
    if interpreted_id != linked_id or (
        revision_view == "current_sample_catalog" and current_id != linked_id
    ):
        body["linked_sample_revision_id"] = linked_id
        body["current_sample_revision_id"] = current_id
    return body


def _experiment_payload(
    session: Session,
    *,
    project_id: str,
    display_name: str,
    assay: object,
    target: object,
    batch: object,
    samples: object,
) -> ExperimentPayload:
    lookup = SessionTermLookup(session)
    assay_slot, target_slot, batch_slot, requests = interpret_experiment_fields(
        assay,
        target,
        batch,
        samples,
        lookup,
        for_write=True,
    )
    pins: list[ExperimentSamplePin] = []
    for ordinal, (request, role) in enumerate(requests):
        sample_revision_id = str(request["sample_revision_id"])
        row = session.get(SampleRevisionRow, sample_revision_id)
        if row is None or row.project_id != project_id:
            raise catalog_payload_invalid()
        group = request.get("replicate_group")
        if group == "":
            group = None
        pins.append(
            ExperimentSamplePin(
                sample_id=row.sample_id,
                sample_revision_id=row.sample_revision_id,
                payload_digest=row.payload_digest,
                ordinal=ordinal,
                replicate_group=group if isinstance(group, str) else None,
                replicate_role=role,
            )
        )
    return ExperimentPayload(
        display_name=display_name,
        assay=assay_slot,
        target=target_slot,
        batch=batch_slot,
        samples=tuple(pins),
    )


def _insert_experiment_revision(
    session: Session,
    *,
    experiment_id: str,
    project_id: str,
    experiment_revision_id: str,
    revision_number: int,
    canonical: str,
    payload: ExperimentPayload,
    created_at: datetime,
    identity_created_at: datetime,
) -> None:
    digest = build_experiment_revision_payload_digest(canonical)
    session.add(
        ExperimentRevisionRow(
            experiment_revision_id=experiment_revision_id,
            experiment_id=experiment_id,
            project_id=project_id,
            revision_number=revision_number,
            canonical_payload=canonical,
            payload_digest_scheme=EXPERIMENT_REVISION_PAYLOAD_DIGEST_SCHEME,
            payload_digest=digest,
            created_at=created_at,
        )
    )
    session.flush()
    head = session.get(ExperimentRevisionHeadRow, experiment_id)
    if head is None:
        session.add(
            ExperimentRevisionHeadRow(
                experiment_id=experiment_id,
                project_id=project_id,
                experiment_revision_id=experiment_revision_id,
                revision_number=revision_number,
                created_at=identity_created_at,
                revision_created_at=created_at,
            )
        )
    else:
        head.experiment_revision_id = experiment_revision_id
        head.revision_number = revision_number
        head.revision_created_at = created_at
    _write_experiment_children(
        session,
        experiment_id=experiment_id,
        project_id=project_id,
        experiment_revision_id=experiment_revision_id,
        payload=payload,
    )
    session.flush()


def _write_experiment_children(
    session: Session,
    *,
    experiment_id: str,
    project_id: str,
    experiment_revision_id: str,
    payload: ExperimentPayload,
) -> None:
    for pin in payload.samples:
        session.add(
            ExperimentRevisionSampleRow(
                experiment_revision_id=experiment_revision_id,
                ordinal=pin.ordinal,
                project_id=project_id,
                sample_id=pin.sample_id,
                sample_revision_id=pin.sample_revision_id,
                payload_digest=pin.payload_digest,
                replicate_group=pin.replicate_group,
                replicate_role_status=pin.replicate_role.status,
                replicate_role_term_id=(
                    None
                    if pin.replicate_role.term is None
                    else pin.replicate_role.term.term_id
                ),
                replicate_role_term_revision_id=(
                    None
                    if pin.replicate_role.term is None
                    else pin.replicate_role.term.term_revision_id
                ),
            )
        )
    session.add(
        ExperimentCatalogProjectionRow(
            experiment_revision_id=experiment_revision_id,
            experiment_id=experiment_id,
            project_id=project_id,
            assay_status=payload.assay.status,
            assay_term_id=_term_id(payload.assay),
            assay_term_revision_id=_term_revision_id(payload.assay),
            target_status=payload.target.status,
            target_term_id=_term_id(payload.target),
            target_term_revision_id=_term_revision_id(payload.target),
            batch_status=payload.batch.status,
            batch_term_id=_term_id(payload.batch),
            batch_term_revision_id=_term_revision_id(payload.batch),
        )
    )


def _write_sample_projection(
    session: Session,
    revision: SampleRevision,
    block: SampleCatalogBlock | None,
    *,
    catalog_trusted: bool,
) -> None:
    session.query(SampleCatalogValueRow).filter(
        SampleCatalogValueRow.sample_revision_id == revision.sample_revision_id
    ).delete(synchronize_session=False)
    session.query(SampleCatalogProjectionRow).filter(
        SampleCatalogProjectionRow.sample_revision_id == revision.sample_revision_id
    ).delete(synchronize_session=False)
    session.flush()
    scalar, values = projection_rows_from_block(block)
    session.add(
        SampleCatalogProjectionRow(
            sample_revision_id=revision.sample_revision_id,
            project_id=revision.project_id,
            sample_id=revision.sample_id,
            catalog_state=scalar["catalog_state"],
            catalog_trusted=1 if catalog_trusted else 0,
            species_status=scalar.get("species_status"),
            species_term_id=scalar.get("species_term_id"),
            species_term_revision_id=scalar.get("species_term_revision_id"),
            tissue_status=scalar.get("tissue_status"),
            tissue_term_id=scalar.get("tissue_term_id"),
            tissue_term_revision_id=scalar.get("tissue_term_revision_id"),
            cell_type_status=scalar.get("cell_type_status"),
            cell_line_status=scalar.get("cell_line_status"),
            treatment_status=scalar.get("treatment_status"),
            tag_status=scalar.get("tag_status"),
        )
    )
    for value in values:
        session.add(
            SampleCatalogValueRow(
                sample_revision_id=revision.sample_revision_id,
                field_name=value["field_name"],
                term_id=value["term_id"],
                term_revision_id=value["term_revision_id"],
            )
        )
    session.flush()


def _install_sqlite_index_hint() -> None:
    """SQLite's compiler drops FROM hints, so INDEXED BY never reaches EXPLAIN."""
    from sqlalchemy.dialects.sqlite.base import SQLiteCompiler

    current = SQLiteCompiler.get_from_hint_text
    if getattr(current, "_catalog_passthrough", False):
        return

    def get_from_hint_text(self, table, hint):
        return hint

    get_from_hint_text._catalog_passthrough = True  # type: ignore[attr-defined]
    SQLiteCompiler.get_from_hint_text = get_from_hint_text  # type: ignore[method-assign]


def _experiment_statement(filters: CatalogFilters):
    _install_sqlite_index_hint()
    selected = (
        ExperimentRow.experiment_id,
        ExperimentRow.created_at,
        ExperimentRevisionHeadRow.experiment_revision_id,
        ExperimentRevisionHeadRow.revision_created_at,
    )
    head_on = ExperimentRevisionHeadRow.experiment_id == ExperimentRow.experiment_id
    projection_on = (
        ExperimentCatalogProjectionRow.experiment_revision_id
        == ExperimentRevisionHeadRow.experiment_revision_id
    )
    assay = filters.field("assay")
    force_assay_index = assay is not None and bool(assay.term_ids)
    if force_assay_index:
        statement = (
            select(*selected)
            .select_from(ExperimentCatalogProjectionRow)
            .join(ExperimentRevisionHeadRow, projection_on)
            .join(ExperimentRow, head_on)
            .where(ExperimentCatalogProjectionRow.project_id == filters.project_id)
        )
    else:
        statement = (
            select(*selected)
            .join(ExperimentRevisionHeadRow, head_on)
            .join(ExperimentCatalogProjectionRow, projection_on)
        )
    statement = statement.where(ExperimentRow.project_id == filters.project_id)
    if not filters.include_archived:
        statement = statement.where(ExperimentRow.archived_at.is_(None))
    if filters.created_from is not None:
        statement = statement.where(
            ExperimentRevisionHeadRow.revision_created_at >= filters.created_from
        )
    if filters.created_before is not None:
        statement = statement.where(
            ExperimentRevisionHeadRow.revision_created_at < filters.created_before
        )
    for name in EXPERIMENT_FIELDS:
        field = filters.field(name)
        if field is None:
            continue
        status_column = getattr(ExperimentCatalogProjectionRow, f"{name}_status")
        term_column = getattr(ExperimentCatalogProjectionRow, f"{name}_term_id")
        if field.status is not None:
            statement = statement.where(status_column == field.status)
        if field.term_ids:
            statement = statement.where(term_column.in_(field.term_ids))
    assay = filters.field("assay")
    if assay is not None and assay.term_ids:
        statement = statement.with_hint(
            ExperimentCatalogProjectionRow,
            "INDEXED BY ix_experiment_catalog_projection_assay",
            dialect_name="sqlite",
        )
    current = filters.revision_view == "current_sample_catalog"
    for name in SAMPLE_FIELDS:
        field = filters.field(name)
        if field is None:
            continue
        statement = statement.where(
            _sample_dimension_exists(name, field, current=current)
        )
    return statement


def _sample_dimension_exists(name: str, field, *, current: bool):
    sample_head = aliased(SampleRevisionHeadRow)
    revision_id = (
        sample_head.sample_revision_id
        if current
        else ExperimentRevisionSampleRow.sample_revision_id
    )
    if (
        name in SAMPLE_MULTI_FIELDS
        and field.term_ids
        and field.status in {None, "known"}
    ):
        statement = (
            select(SampleCatalogValueRow.sample_revision_id)
            .select_from(ExperimentRevisionSampleRow)
            .where(
                ExperimentRevisionSampleRow.experiment_revision_id
                == ExperimentRevisionHeadRow.experiment_revision_id
            )
        )
        if current:
            statement = statement.join(
                sample_head,
                sample_head.sample_id == ExperimentRevisionSampleRow.sample_id,
            )
        statement = statement.join(
            SampleCatalogValueRow,
            and_(
                SampleCatalogValueRow.sample_revision_id == revision_id,
                SampleCatalogValueRow.field_name == name,
                SampleCatalogValueRow.term_id.in_(field.term_ids),
            ),
        )
        return statement.exists()
    statement = (
        select(SampleCatalogProjectionRow.sample_revision_id)
        .select_from(ExperimentRevisionSampleRow)
        .where(
            ExperimentRevisionSampleRow.experiment_revision_id
            == ExperimentRevisionHeadRow.experiment_revision_id
        )
    )
    if current:
        statement = statement.join(
            sample_head,
            sample_head.sample_id == ExperimentRevisionSampleRow.sample_id,
        )
    statement = statement.join(
        SampleCatalogProjectionRow,
        SampleCatalogProjectionRow.sample_revision_id == revision_id,
    )
    if field.status is not None:
        statement = statement.where(
            getattr(SampleCatalogProjectionRow, f"{name}_status") == field.status
        )
    if field.term_ids and name in SAMPLE_SINGLE_FIELDS:
        statement = statement.where(
            getattr(SampleCatalogProjectionRow, f"{name}_term_id").in_(field.term_ids)
        )
    return statement.exists()


def _sample_statement(filters: CatalogFilters):
    statement = (
        select(
            SampleRow.sample_id,
            SampleRow.created_at,
            SampleRevisionHeadRow.sample_revision_id,
        )
        .join(
            SampleRevisionHeadRow,
            SampleRevisionHeadRow.sample_id == SampleRow.sample_id,
        )
        .outerjoin(
            SampleCatalogProjectionRow,
            SampleCatalogProjectionRow.sample_revision_id
            == SampleRevisionHeadRow.sample_revision_id,
        )
        .where(SampleRow.project_id == filters.project_id)
    )
    if filters.catalog_state == "absent":
        statement = statement.where(
            or_(
                SampleCatalogProjectionRow.sample_revision_id.is_(None),
                SampleCatalogProjectionRow.catalog_state == "absent",
            )
        )
    elif filters.catalog_state == "present":
        statement = statement.where(
            SampleCatalogProjectionRow.catalog_state == "present"
        )
    if filters.created_from is not None:
        statement = statement.where(SampleRow.created_at >= filters.created_from)
    if filters.created_before is not None:
        statement = statement.where(SampleRow.created_at < filters.created_before)
    for name in SAMPLE_FIELDS:
        field = filters.field(name)
        if field is None:
            continue
        if (
            name in SAMPLE_MULTI_FIELDS
            and field.term_ids
            and field.status in {None, "known"}
        ):
            statement = statement.where(
                select(SampleCatalogValueRow.sample_revision_id)
                .where(
                    SampleCatalogValueRow.sample_revision_id
                    == SampleRevisionHeadRow.sample_revision_id,
                    SampleCatalogValueRow.field_name == name,
                    SampleCatalogValueRow.term_id.in_(field.term_ids),
                )
                .exists()
            )
            continue
        if field.status is not None:
            statement = statement.where(
                getattr(SampleCatalogProjectionRow, f"{name}_status") == field.status
            )
        if field.term_ids and name in SAMPLE_SINGLE_FIELDS:
            statement = statement.where(
                getattr(SampleCatalogProjectionRow, f"{name}_term_id").in_(
                    field.term_ids
                )
            )
    return statement


def _term_statement(filters: CatalogFilters):
    statement = select(
        CatalogTermRow.term_id,
        CatalogTermRow.created_at,
        CatalogTermRow.term_type,
    )
    if not filters.include_archived:
        statement = statement.where(CatalogTermRow.archived_at.is_(None))
    if filters.term_type is not None:
        statement = statement.where(CatalogTermRow.term_type == filters.term_type)
    if filters.created_from is not None:
        statement = statement.where(CatalogTermRow.created_at >= filters.created_from)
    if filters.created_before is not None:
        statement = statement.where(CatalogTermRow.created_at < filters.created_before)
    return statement


def _publication_statement(project_id: str, run_ids: Sequence[str]):
    return (
        select(ArtifactPublicationRow.run_id)
        .join(RunRow, RunRow.run_id == ArtifactPublicationRow.run_id)
        .join(
            RunProjectBindingRow,
            RunProjectBindingRow.run_id == ArtifactPublicationRow.run_id,
        )
        .where(
            ArtifactPublicationRow.run_id.in_(tuple(run_ids)),
            RunProjectBindingRow.project_id == project_id,
        )
    )


def _sample_trace_results(
    session: Session,
    project_id: str,
    sample_revision_id: str,
) -> list[dict[str, Any]]:
    run_ids = list(
        session.scalars(
            select(RunSampleRow.run_id).where(
                RunSampleRow.sample_revision_id == sample_revision_id,
                RunSampleRow.project_id == project_id,
            )
        ).all()
    )
    if not run_ids:
        return []
    published = _published_trace_rows(session, project_id, run_ids)
    published_run_ids = {item["run_id"] for item in published}
    unpublished = []
    bindings = session.scalars(
        select(RunProjectBindingRow).where(
            RunProjectBindingRow.run_id.in_(run_ids),
            RunProjectBindingRow.project_id == project_id,
        )
    ).all()
    for binding in bindings:
        if binding.run_id in published_run_ids:
            continue
        unpublished.append(
            {
                "run_id": binding.run_id,
                "provenance": binding.provenance,
                "publication_evidence": "absent",
                "_sort": (1, binding.run_id),
            }
        )
    unpublished.sort(key=lambda item: item["run_id"], reverse=True)
    results = []
    for item in [*published, *unpublished]:
        item.pop("_sort", None)
        results.append(item)
    return results


def _published_trace_rows(
    session: Session,
    project_id: str,
    run_ids: Sequence[str],
) -> list[dict[str, Any]]:
    if not run_ids:
        return []
    publication_rows = session.execute(
        select(
            ArtifactPublicationRow,
            RunRow.workflow_id,
            RunProjectBindingRow.binding_mode,
            RunProjectBindingRow.provenance,
            RunResultStateRow.artifact_generation,
        )
        .join(RunRow, RunRow.run_id == ArtifactPublicationRow.run_id)
        .join(
            RunProjectBindingRow,
            RunProjectBindingRow.run_id == ArtifactPublicationRow.run_id,
        )
        .join(
            RunResultStateRow,
            RunResultStateRow.run_id == ArtifactPublicationRow.run_id,
        )
        .where(
            ArtifactPublicationRow.run_id.in_(run_ids),
            RunProjectBindingRow.project_id == project_id,
        )
    ).all()
    samples_by_run = _associated_samples(session, run_ids)
    published = []
    for publication, workflow_id, binding_mode, provenance, current in publication_rows:
        current_generation = current or publication.artifact_generation
        published.append(
            {
                "run_id": publication.run_id,
                "provenance": provenance,
                "publication_evidence": "present",
                "artifact_association_kind": "run_input_scope",
                "workflow_id": workflow_id,
                "artifact_id": publication.artifact_id,
                "output_type": publication.output_type,
                "artifact_generation": publication.artifact_generation,
                "artifact_revision": publication.artifact_revision,
                "published_at": _stamp(publication.published_at),
                "current_artifact_generation": current_generation,
                "generation_status": (
                    "current"
                    if publication.artifact_generation == current_generation
                    else "superseded"
                ),
                "run_sample_binding": {
                    "binding_mode": binding_mode,
                    "provenance": provenance,
                    "associated_run_samples": samples_by_run.get(
                        publication.run_id, []
                    ),
                },
                "_sort": (
                    0,
                    _as_utc(publication.published_at),
                    publication.run_id,
                    publication.artifact_generation,
                    publication.artifact_id,
                ),
            }
        )
    published.sort(
        key=lambda item: (
            item["_sort"][1],
            item["_sort"][2],
            item["_sort"][3],
            item["_sort"][4],
        ),
        reverse=True,
    )
    for item in published:
        item.pop("_sort", None)
    return published


def _associated_samples(
    session: Session,
    run_ids: Sequence[str],
) -> dict[str, list[dict[str, Any]]]:
    rows = session.execute(
        select(RunSampleRow, SampleRevisionRow)
        .join(
            SampleRevisionRow,
            and_(
                SampleRevisionRow.project_id == RunSampleRow.project_id,
                SampleRevisionRow.sample_revision_id == RunSampleRow.sample_revision_id,
                SampleRevisionRow.payload_digest == RunSampleRow.payload_digest,
            ),
        )
        .where(RunSampleRow.run_id.in_(tuple(run_ids)))
        .order_by(RunSampleRow.run_id, RunSampleRow.ordinal)
    ).all()
    grouped: dict[str, list[dict[str, Any]]] = {}
    for run_sample, revision in rows:
        grouped.setdefault(run_sample.run_id, []).append(
            {
                "sample_id": revision.sample_id,
                "sample_revision_id": revision.sample_revision_id,
                "revision_number": revision.revision_number,
                "ordinal": run_sample.ordinal,
            }
        )
    return grouped


def _trace_object_id(item: Mapping[str, Any]) -> str:
    if item["publication_evidence"] == "absent":
        return "unpublished\n" + item["run_id"]
    return "\n".join(
        (
            "publication",
            item["run_id"],
            item["artifact_generation"],
            item["artifact_id"],
        )
    )


def _parse_trace_stamp(item: Mapping[str, Any]) -> datetime:
    if item["publication_evidence"] == "absent":
        return datetime(1970, 1, 1, tzinfo=timezone.utc)
    value = item["published_at"]
    parsed = datetime.fromisoformat(str(value).replace("Z", "+00:00"))
    return _as_utc(parsed)


def _trace_index(
    results: Sequence[Mapping[str, Any]], position: CatalogCursorPosition
) -> int:
    for index, item in enumerate(results):
        if _trace_object_id(item) != position.object_id:
            continue
        if _as_utc(_parse_trace_stamp(item)) != _as_utc(position.created_at):
            continue
        return index + 1
    raise catalog_cursor_not_found()


def _block_for_revision(
    session: Session,
    sample_revision_id: str,
) -> SampleCatalogBlock | None:
    scalar = session.get(SampleCatalogProjectionRow, sample_revision_id)
    if scalar is None or scalar.catalog_state == "absent":
        return None
    values = session.execute(
        select(
            SampleCatalogValueRow.field_name,
            SampleCatalogValueRow.term_id,
            SampleCatalogValueRow.term_revision_id,
        ).where(SampleCatalogValueRow.sample_revision_id == sample_revision_id)
    ).all()
    return block_from_projection(
        {
            "catalog_state": scalar.catalog_state,
            "species_status": scalar.species_status,
            "species_term_id": scalar.species_term_id,
            "species_term_revision_id": scalar.species_term_revision_id,
            "tissue_status": scalar.tissue_status,
            "tissue_term_id": scalar.tissue_term_id,
            "tissue_term_revision_id": scalar.tissue_term_revision_id,
            "cell_type_status": scalar.cell_type_status,
            "cell_line_status": scalar.cell_line_status,
            "treatment_status": scalar.treatment_status,
            "tag_status": scalar.tag_status,
        },
        [
            {
                "field_name": field_name,
                "term_id": term_id,
                "term_revision_id": term_revision_id,
            }
            for field_name, term_id, term_revision_id in values
        ],
    )


def _pins_from_stored_payload(
    session: Session,
    project_id: str,
    payload: ExperimentPayload,
) -> tuple[ExperimentSamplePin, ...]:
    pins: list[ExperimentSamplePin] = []
    for pin in payload.samples:
        row = session.get(SampleRevisionRow, pin.sample_revision_id)
        if (
            row is None
            or row.project_id != project_id
            or row.payload_digest != pin.payload_digest
        ):
            raise catalog_payload_invalid()
        pins.append(
            ExperimentSamplePin(
                sample_id=row.sample_id,
                sample_revision_id=row.sample_revision_id,
                payload_digest=row.payload_digest,
                ordinal=pin.ordinal,
                replicate_group=pin.replicate_group,
                replicate_role=pin.replicate_role,
            )
        )
    return tuple(pins)


def _experiment_name_ids(session, payload, pins, revision_view: str) -> set[str]:
    found = collect_experiment_term_revision_ids(payload)
    for pin in pins:
        interpreted = pin.sample_revision_id
        if revision_view == "current_sample_catalog":
            head = session.get(SampleRevisionHeadRow, pin.sample_id)
            if head is not None:
                interpreted = head.sample_revision_id
        block = _block_for_revision(session, interpreted)
        found.update(collect_term_revision_ids(block))
        if pin.replicate_role_term_revision_id:
            found.add(pin.replicate_role_term_revision_id)
    return found


def _names_for(session: Session, revision_ids: set[str]) -> dict[str, str]:
    if not revision_ids:
        return {}
    rows = session.scalars(
        select(CatalogTermRevisionRow).where(
            CatalogTermRevisionRow.term_revision_id.in_(tuple(revision_ids))
        )
    ).all()
    return {row.term_revision_id: row.display_name for row in rows}


def _term_revision_row(
    *,
    term_id: str,
    term_revision_id: str,
    revision_number: int,
    display_name: str,
    namespace: str | None,
    accession: str | None,
    aliases: tuple[str, ...],
    change_kind: str,
    created_at: datetime,
) -> CatalogTermRevisionRow:
    return CatalogTermRevisionRow(
        term_revision_id=term_revision_id,
        term_id=term_id,
        revision_number=revision_number,
        display_name=display_name,
        namespace=namespace,
        accession=accession,
        aliases_json=canonical_alias_json(validate_aliases(list(aliases))),
        change_kind=change_kind,
        created_at=created_at,
    )


def _term_body(session: Session, term_id: str, term_revision_id: str) -> dict[str, Any]:
    term = session.get(CatalogTermRow, term_id)
    revision = session.get(CatalogTermRevisionRow, term_revision_id)
    if term is None or revision is None:
        raise catalog_not_found()
    return {
        "term_id": term.term_id,
        "term_type": term.term_type,
        "archived_at": _stamp(term.archived_at),
        "revision": {
            "term_revision_id": revision.term_revision_id,
            "revision_number": revision.revision_number,
            "display_name": revision.display_name,
            "namespace": revision.namespace,
            "accession": revision.accession,
            "aliases": json.loads(revision.aliases_json),
            "change_kind": revision.change_kind,
            "created_at": _stamp(revision.created_at),
        },
    }


def _term_head_number(session: Session, term_id: str) -> int:
    value = session.scalar(
        select(CatalogTermRevisionRow.revision_number)
        .where(CatalogTermRevisionRow.term_id == term_id)
        .order_by(CatalogTermRevisionRow.revision_number.desc())
        .limit(1)
    )
    if value is None:
        raise catalog_not_found()
    return int(value)


def _term_head_id(session: Session, term_id: str) -> str:
    value = session.scalar(
        select(CatalogTermRevisionRow.term_revision_id)
        .where(CatalogTermRevisionRow.term_id == term_id)
        .order_by(CatalogTermRevisionRow.revision_number.desc())
        .limit(1)
    )
    if value is None:
        raise catalog_not_found()
    return str(value)


def _checked_write_receipt(
    session: Session, write: CatalogWrite
) -> CatalogWriteReceiptRow | None:
    _require_administrator(session, write.actor_user_id)
    if not write.target_key:
        raise ValueError("catalog write target is required")
    receipt = session.get(
        CatalogWriteReceiptRow, (write.actor_user_id, write.idempotency_key)
    )
    if receipt is not None and not matches_committed_catalog_request(
        write,
        operation=receipt.operation,
        target_key=receipt.target_key,
        request_sha256=receipt.request_sha256,
    ):
        raise catalog_idempotency_conflict()
    return receipt


def _replay(
    repository, session: Session, receipt: CatalogWriteReceiptRow
) -> dict[str, Any]:
    operation = receipt.operation
    if operation in {"catalogCreateTerm", "catalogReviseTerm", "catalogArchiveTerm"}:
        revision_id = receipt.result_revision_id or _term_head_id(
            session,
            receipt.subject_id,
        )
        return _term_body(session, receipt.subject_id, revision_id)
    if operation == "catalogArchiveExperiment":
        return {
            "experiment_id": receipt.subject_id,
            "archived_at": _stamp(receipt.result_archived_at),
            "revision_number": _experiment_head_number(session, receipt.subject_id),
        }
    if operation in {"catalogCreateExperiment", "catalogReviseExperiment"}:
        if receipt.result_revision_id is None:
            raise catalog_not_found()
        return repository._experiment_detail(
            session,
            receipt.result_revision_id,
            include_candidates=True,
        )
    if operation == "catalogReviseSampleCatalog":
        if receipt.result_revision_id is None:
            raise catalog_not_found()
        return repository._sample_revision_body(
            session,
            receipt.result_revision_id,
            include_candidates=True,
        )
    if operation in {
        "catalogRegisterProposal",
        "catalogGenerateClassificationSuggestion",
    }:
        from encode_pipeline.persistence.catalog_proposal_store import replay_register

        return replay_register(
            session,
            receipt.actor_user_id,
            receipt.idempotency_key,
        )
    if operation == "catalogAdoptProposal":
        from encode_pipeline.persistence.catalog_proposal_store import replay_adopt

        return replay_adopt(
            session,
            receipt.actor_user_id,
            receipt.idempotency_key,
        )
    raise catalog_idempotency_conflict()


def _target_term_record(session: Session, term_id: str) -> dict[str, Any]:
    revision_id = _term_head_id(session, term_id)
    body = _term_body(session, term_id, revision_id)
    revision = body["revision"]
    return {
        "term_id": body["term_id"],
        "term_revision_id": revision["term_revision_id"],
        "display_name": revision["display_name"],
        "aliases": list(revision["aliases"]),
        "archived": body["archived_at"] is not None,
        "requested_is_head": True,
    }


def _experiment_head_number(session: Session, experiment_id: str) -> int:
    head = session.get(ExperimentRevisionHeadRow, experiment_id)
    if head is None:
        raise catalog_not_found()
    return head.revision_number


def _add_receipt(
    session: Session,
    write: CatalogWrite,
    *,
    subject_id: str,
    result_revision_id: str | None,
    result_archived_at: datetime | None,
) -> None:
    session.add(
        CatalogWriteReceiptRow(
            actor_user_id=write.actor_user_id,
            idempotency_key=write.idempotency_key,
            operation=write.operation,
            request_sha256=write.request_sha256,
            target_key=write.target_key,
            subject_id=subject_id,
            result_revision_id=result_revision_id,
            result_archived_at=result_archived_at,
            created_at=write.now,
        )
    )
    session.flush()


def _require_administrator(session: Session, user_id: str) -> None:
    account = session.get(UserAccountRow, user_id)
    if (
        account is None
        or account.status != "enabled"
        or account.role != "administrator"
    ):
        raise administrator_required()


def _require_active_project(session: Session, project_id: str) -> ProjectRow:
    row = session.get(ProjectRow, project_id)
    if (
        row is None
        or row.archived_at is not None
        or row.project_id == LEGACY_PROJECT_ID
        or row.kind != "user"
    ):
        raise catalog_not_found()
    return row


def _owner_project(session: Session, kind: str, owner_id: str) -> str:
    if kind == "sample_revisions":
        row = session.get(SampleRow, owner_id)
        if row is None:
            raise catalog_not_found()
        return row.project_id
    row = session.get(ExperimentRow, owner_id)
    if row is None:
        raise catalog_not_found()
    return row.project_id


def _expect_revision(current: int, expected: object) -> None:
    if isinstance(expected, bool) or not isinstance(expected, int) or expected < 1:
        raise catalog_payload_invalid()
    if current != expected:
        raise catalog_revision_conflict()


def catalog_term_archived_error() -> CatalogError:
    from encode_pipeline.platform.experiment_catalog import catalog_term_archived

    return catalog_term_archived()


def _required_project(filters: CatalogFilters) -> str:
    if filters.project_id is None:
        raise catalog_not_found()
    return filters.project_id


def _require_cursor_row(
    session, statement, id_column, created_column, position
) -> None:
    row = session.execute(statement.where(id_column == position.object_id)).first()
    if row is None:
        raise catalog_cursor_not_found()
    created_at = getattr(row, created_column.key, None)
    if created_at is None:
        created_at = row[1]
    if _as_utc(created_at) != _as_utc(position.created_at):
        raise catalog_cursor_not_found()


def _after(created_column, id_column, position: CatalogCursorPosition):
    return or_(
        created_column > position.created_at,
        and_(
            created_column == position.created_at,
            id_column > position.object_id,
        ),
    )


def _next_cursor(
    rows, limit: int, *, kind: str, filters: CatalogFilters, created_at, object_id
):
    if len(rows) <= limit:
        return None
    last = rows[limit - 1]
    return encode_catalog_cursor(
        kind=kind,
        filters=filters,
        position=CatalogCursorPosition(
            created_at=_as_utc(created_at(last)),
            object_id=object_id(last),
        ),
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


def _term_id(slot: SingleSlot) -> str | None:
    return None if slot.term is None else slot.term.term_id


def _term_revision_id(slot: SingleSlot) -> str | None:
    return None if slot.term is None else slot.term.term_revision_id


def _explain(session: Session, statement) -> str:
    compiled = statement.compile(
        session.get_bind(),
        compile_kwargs={"literal_binds": True},
    )
    rows = session.execute(text("EXPLAIN QUERY PLAN " + str(compiled))).all()
    details = []
    for row in rows:
        details.append(str(row[3]) if len(row) > 3 else str(tuple(row)))
    return "\n".join(details)


def _begin_write(session: Session) -> None:
    if session.get_bind().dialect.name == "sqlite":
        session.connection().exec_driver_sql("BEGIN IMMEDIATE")


def _as_utc(value: datetime) -> datetime:
    if value.tzinfo is None or value.utcoffset() is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _stamp(value: datetime | None) -> str | None:
    if value is None:
        return None
    text_value = _as_utc(value).isoformat(timespec="microseconds")
    if text_value.endswith(".000000+00:00"):
        text_value = text_value.replace(".000000+00:00", "+00:00")
    return text_value.replace("+00:00", "Z")
