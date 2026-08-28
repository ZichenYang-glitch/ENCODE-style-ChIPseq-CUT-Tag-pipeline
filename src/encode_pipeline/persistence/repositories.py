"""SQLAlchemy implementation of the workflow run repository boundary."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import replace
from datetime import datetime, timezone
from threading import RLock

from sqlalchemy import and_, case, delete, func, or_, select, update
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session, sessionmaker

from encode_pipeline.persistence.models import (
    ArtifactPublicationRow,
    ProjectRow,
    ReferenceProfileRow,
    RunReferenceBindingRow,
    RunInputBindingRow,
    RunInputMemberRow,
    RunInputUseRow,
    RunArtifactRow,
    RunEventRow,
    RunExecutionAssignmentRow,
    RunLogRow,
    RunQcMetricRow,
    RunResultAttemptRow,
    RunResultStateRow,
    RunProjectBindingRow,
    RunRow,
    RunSampleRow,
    RunWorkflowBuildIdentityRow,
    SampleRevisionRow,
    SnapshotProjectBindingRow,
    SnapshotReferenceBindingRow,
    SnapshotInputBindingRow,
    SnapshotInputMemberRow,
    SnapshotInputUseRow,
    SnapshotSampleRevisionRow,
    ValidatedInputSnapshotRow,
)
from encode_pipeline.persistence.authentication import (
    insert_security_audit_event,
)
from encode_pipeline.persistence.input_registry import (
    resolve_input_use_binding_plan_in_session,
)
from encode_pipeline.platform.artifact_publications import (
    ArtifactPublicationCursorPosition,
    ArtifactPublicationFilters,
    ArtifactPublicationRunSampleBinding,
    ArtifactPublicationSummary,
    AssociatedRunSample,
)
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.data_registry import (
    BindingMode,
    BindingProvenance,
    ProjectKind,
    ProjectSampleBinding,
    ProjectSampleSelection,
    SampleRevision,
    SampleRevisionBindingRef,
    build_legacy_project_sample_binding,
    build_project_sample_binding,
)
from encode_pipeline.platform.execution import (
    RunExecutionAssignment,
    RunExecutionCancellationRequest,
    RunExecutionClaim,
    RunExecutionOwnership,
    RunExecutionRequeuePreparation,
    RunExecutionStopAcknowledgement,
)
from encode_pipeline.platform.input_registry import (
    InputBindingContractMode,
    InputFileRevisionBindingRef,
    InputProvenanceMode,
    InputUseBinding,
    InputUseBindingEnvelope,
    InputUseBindingPlan,
    build_compatibility_input_binding,
)
from encode_pipeline.platform.results import Issue
from encode_pipeline.platform.reference_profiles import ReferenceProfileRevisionBinding
from encode_pipeline.platform.result_generations import (
    RunResultState,
    artifact_manifest_digest,
    build_artifact_generation,
    build_qc_generation,
    qc_metric_manifest_digest,
    validate_artifact_generation,
    validate_qc_generation,
    validate_result_attempt_id,
)
from encode_pipeline.platform.run_history import RunHistoryCursor, RunSummary
from encode_pipeline.platform.security_audit import SecurityAuditEvent
from encode_pipeline.platform.runs import (
    RunArtifactRef,
    RunEvent,
    RunLogChunk,
    RunQcMetric,
    RunRecord,
    RunStatus,
)
from encode_pipeline.platform.snapshots import (
    ValidatedInputSnapshot,
    ValidatedSnapshotRunCreation,
)
from encode_pipeline.services.run_repositories import (
    ConcurrentRunUpdateError,
    InputBindingSelectionError,
    ProjectSampleSelectionError,
    ReferenceBindingSelectionError,
    ResultGenerationChangedError,
    RunEventDraft,
    ValidatedSnapshotExpiredError,
    ValidatedSnapshotReplayConflictError,
    _artifact_publication_output_type,
    _assignment_with_bound_cleanup_identity,
    _event_with_context,
    _legacy_run_inputs_digest,
    _require_current_attempt,
    _sorted_artifacts,
    _state_after_artifact_change,
    _validate_qc_metric_fields,
    _validated_expected_artifacts,
    _validated_qc_replacement,
    _validate_snapshot_consumption_candidate,
    _validate_snapshot_linked_run,
    canonical_decimal_text,
    decimal_from_canonical_text,
)


class SqlAlchemyRunRepository:
    """Persist run aggregates without exposing ORM rows to service callers."""

    def __init__(self, session_factory: sessionmaker[Session]) -> None:
        if not isinstance(session_factory, sessionmaker):
            raise ValueError("session_factory must be a SQLAlchemy sessionmaker")
        self._session_factory = session_factory
        self._lock = RLock()

    def contains_run(self, run_id: str) -> bool:
        with self._session_factory() as session:
            return (
                session.scalar(select(RunRow.id).where(RunRow.run_id == run_id))
                is not None
            )

    def create_run(
        self,
        record: RunRecord,
        event: RunEventDraft,
        *,
        requested_by_user_id: str | None = None,
        security_audit: SecurityAuditEvent | None = None,
    ) -> RunEvent:
        _validate_security_audit(security_audit)
        with self._lock:
            try:
                with self._session_factory.begin() as session:
                    _begin_write(session)
                    session.add(
                        _run_row(
                            record,
                            requested_by_user_id=requested_by_user_id,
                        )
                    )
                    session.flush()
                    session.add(RunResultStateRow(run_id=record.run_id))
                    binding = build_legacy_project_sample_binding(
                        _legacy_run_inputs_digest(record)
                    )
                    _add_run_data_binding(
                        session,
                        record.run_id,
                        binding,
                        created_at=record.created_at,
                    )
                    session.flush()
                    _add_run_input_binding(
                        session,
                        record.run_id,
                        build_compatibility_input_binding(
                            project_id=binding.project_id,
                            project_sample_binding_digest=binding.digest,
                            workflow_id=record.workflow_id,
                            adapter_contract_version=None,
                            workflow_inputs_digest=binding.workflow_inputs_digest,
                        ),
                        created_at=record.created_at,
                    )
                    session.flush()
                    inserted = self._insert_event(session, record.run_id, event)
                    if security_audit is not None:
                        insert_security_audit_event(session, security_audit)
                    return inserted
            except IntegrityError as exc:
                raise ValueError(f"Duplicate run_id: {record.run_id!r}") from exc

    def create_validated_input_snapshot(
        self,
        snapshot: ValidatedInputSnapshot,
        *,
        project_sample_selection: ProjectSampleSelection | None = None,
        input_use_binding_plan: InputUseBindingPlan | None = None,
        reference_binding: ReferenceProfileRevisionBinding | None = None,
    ) -> ValidatedInputSnapshot:
        validated = _validated_input_snapshot_from_row(
            _validated_input_snapshot_row(snapshot)
        )
        if input_use_binding_plan is not None:
            if project_sample_selection is None:
                raise InputBindingSelectionError(
                    "managed input selection requires an exact Project binding"
                )
            if (
                input_use_binding_plan.project_id != project_sample_selection.project_id
                or input_use_binding_plan.workflow_id != validated.workflow_id
            ):
                raise InputBindingSelectionError(
                    "input-use plan does not match the validation scope"
                )
        if reference_binding is not None and (
            not isinstance(reference_binding, ReferenceProfileRevisionBinding)
            or reference_binding.workflow_id != validated.workflow_id
        ):
            raise ReferenceBindingSelectionError(
                "Reference Profile selection does not match the validation scope"
            )
        try:
            with self._lock, self._session_factory.begin() as session:
                _begin_write(session)
                session.add(_validated_input_snapshot_row(validated))
                session.flush()
                if project_sample_selection is None:
                    binding = build_legacy_project_sample_binding(
                        validated.payload_digest
                    )
                else:
                    try:
                        binding = _resolve_project_sample_selection(
                            session,
                            project_sample_selection,
                            workflow_inputs_digest=validated.payload_digest,
                        )
                    except (KeyError, ValueError) as exc:
                        raise ProjectSampleSelectionError(
                            "Project/SampleRevision selection is not eligible"
                        ) from exc
                _add_snapshot_data_binding(
                    session,
                    validated.snapshot_id,
                    binding,
                    created_at=validated.validated_at,
                )
                session.flush()
                if input_use_binding_plan is None:
                    input_binding = build_compatibility_input_binding(
                        project_id=binding.project_id,
                        project_sample_binding_digest=binding.digest,
                        workflow_id=validated.workflow_id,
                        adapter_contract_version=None,
                        workflow_inputs_digest=validated.payload_digest,
                    )
                else:
                    try:
                        input_binding = resolve_input_use_binding_plan_in_session(
                            session,
                            input_use_binding_plan,
                            project_sample_binding_digest=binding.digest,
                            workflow_inputs_digest=validated.payload_digest,
                        )
                    except (KeyError, ValueError) as exc:
                        raise InputBindingSelectionError(
                            "input-use revision selection is not eligible"
                        ) from exc
                _add_snapshot_input_binding(
                    session,
                    validated.snapshot_id,
                    input_binding,
                    created_at=validated.validated_at,
                )
                session.flush()
                if reference_binding is not None:
                    profile = session.get(
                        ReferenceProfileRow,
                        reference_binding.profile_id,
                    )
                    if (
                        profile is None
                        or profile.enabled_revision_id != reference_binding.revision_id
                    ):
                        raise ReferenceBindingSelectionError(
                            "Reference Profile selection is not eligible"
                        )
                    session.add(
                        _snapshot_reference_binding_row(
                            validated.snapshot_id,
                            reference_binding,
                        )
                    )
                    session.flush()
        except (InputBindingSelectionError, ReferenceBindingSelectionError):
            raise
        except IntegrityError as exc:
            raise ValueError(
                f"Duplicate validated snapshot ID: {validated.snapshot_id!r}"
            ) from exc
        return validated

    def get_validated_input_snapshot(
        self,
        snapshot_id: str,
    ) -> ValidatedInputSnapshot:
        with self._session_factory() as session:
            row = session.get(ValidatedInputSnapshotRow, snapshot_id)
            if row is None:
                raise KeyError(snapshot_id)
            return _validated_input_snapshot_from_row(row)

    def get_validated_input_binding(
        self,
        snapshot_id: str,
    ) -> ProjectSampleBinding:
        with self._session_factory() as session:
            snapshot_row = session.get(ValidatedInputSnapshotRow, snapshot_id)
            if snapshot_row is None:
                raise KeyError(snapshot_id)
            return _snapshot_data_binding_from_rows(
                session,
                snapshot_id,
                expected_workflow_inputs_digest=snapshot_row.payload_digest,
            )

    def get_validated_input_use_binding(
        self,
        snapshot_id: str,
    ) -> InputUseBindingEnvelope:
        with self._session_factory() as session:
            snapshot_row = session.get(ValidatedInputSnapshotRow, snapshot_id)
            if snapshot_row is None:
                raise KeyError(snapshot_id)
            data_binding = _snapshot_data_binding_from_rows(
                session,
                snapshot_id,
                expected_workflow_inputs_digest=snapshot_row.payload_digest,
            )
            return _snapshot_input_binding_from_rows(
                session,
                snapshot_id,
                expected_data_binding=data_binding,
                expected_workflow_id=snapshot_row.workflow_id,
                expected_workflow_inputs_digest=snapshot_row.payload_digest,
            )

    def get_validated_reference_binding(
        self,
        snapshot_id: str,
    ) -> ReferenceProfileRevisionBinding | None:
        with self._session_factory() as session:
            if session.get(ValidatedInputSnapshotRow, snapshot_id) is None:
                raise KeyError(snapshot_id)
            row = session.get(SnapshotReferenceBindingRow, snapshot_id)
            return None if row is None else _reference_binding_from_row(row)

    def consume_validated_input_snapshot(
        self,
        snapshot_id: str,
        *,
        workflow_id: str,
        expected_build_identity: WorkflowBuildIdentity,
        record: RunRecord,
        consumed_at: datetime,
        event: RunEventDraft,
        requested_by_user_id: str | None = None,
        security_audit: SecurityAuditEvent | None = None,
    ) -> ValidatedSnapshotRunCreation:
        _validate_security_audit(security_audit)
        try:
            with self._lock, self._session_factory.begin() as session:
                _begin_write(session)
                row = session.get(ValidatedInputSnapshotRow, snapshot_id)
                if row is None:
                    raise KeyError(snapshot_id)
                snapshot = _validated_input_snapshot_from_row(row)
                snapshot_binding = _snapshot_data_binding_from_rows(
                    session,
                    snapshot_id,
                    expected_workflow_inputs_digest=snapshot.payload_digest,
                )
                snapshot_input_binding = _snapshot_input_binding_from_rows(
                    session,
                    snapshot_id,
                    expected_data_binding=snapshot_binding,
                    expected_workflow_id=snapshot.workflow_id,
                    expected_workflow_inputs_digest=snapshot.payload_digest,
                )
                snapshot_reference_row = session.get(
                    SnapshotReferenceBindingRow,
                    snapshot_id,
                )
                snapshot_reference_binding = (
                    None
                    if snapshot_reference_row is None
                    else _reference_binding_from_row(snapshot_reference_row)
                )
                if snapshot.workflow_id != workflow_id:
                    raise KeyError(snapshot_id)
                _validate_snapshot_consumption_candidate(
                    snapshot,
                    workflow_id=workflow_id,
                    expected_build_identity=expected_build_identity,
                    record=record,
                    consumed_at=consumed_at,
                )
                if snapshot.consumed_run_id is not None:
                    run_row = session.scalar(
                        select(RunRow).where(RunRow.run_id == snapshot.consumed_run_id)
                    )
                    if run_row is None:
                        raise ValueError("validated snapshot references a missing run")
                    current = _record_from_row(run_row)
                    _validate_snapshot_linked_run(snapshot, current)
                    if dict(current.tags) != dict(record.tags):
                        raise ValidatedSnapshotReplayConflictError(
                            "validated snapshot replay metadata differs"
                        )
                    try:
                        run_binding = _run_data_binding_from_rows(
                            session,
                            current.run_id,
                        )
                        run_input_binding = _run_input_binding_from_rows(
                            session,
                            current.run_id,
                            expected_data_binding=run_binding,
                            expected_workflow_id=current.workflow_id,
                        )
                        run_reference_row = session.get(
                            RunReferenceBindingRow,
                            current.run_id,
                        )
                        run_reference_binding = (
                            None
                            if run_reference_row is None
                            else _reference_binding_from_row(run_reference_row)
                        )
                    except (KeyError, ValueError) as exc:
                        raise ValidatedSnapshotReplayConflictError(
                            "validated snapshot replay binding evidence is invalid"
                        ) from exc
                    if (
                        run_binding != snapshot_binding
                        or run_input_binding != snapshot_input_binding
                        or run_reference_binding != snapshot_reference_binding
                    ):
                        raise ValidatedSnapshotReplayConflictError(
                            "validated snapshot replay binding evidence differs"
                        )
                    return ValidatedSnapshotRunCreation(
                        record=current,
                        created=False,
                    )
                if consumed_at >= snapshot.expires_at:
                    raise ValidatedSnapshotExpiredError(
                        "validated snapshot expired before first use"
                    )
                if snapshot_reference_binding is not None:
                    profile = session.get(
                        ReferenceProfileRow,
                        snapshot_reference_binding.profile_id,
                    )
                    if (
                        profile is None
                        or profile.enabled_revision_id
                        != snapshot_reference_binding.revision_id
                    ):
                        raise ReferenceBindingSelectionError(
                            "Reference Profile selection is not eligible"
                        )

                session.add(
                    _run_row(
                        record,
                        requested_by_user_id=requested_by_user_id,
                    )
                )
                session.flush()
                session.add(RunResultStateRow(run_id=record.run_id))
                _add_run_data_binding(
                    session,
                    record.run_id,
                    snapshot_binding,
                    created_at=consumed_at,
                )
                session.flush()
                _add_run_input_binding(
                    session,
                    record.run_id,
                    snapshot_input_binding,
                    created_at=consumed_at,
                )
                session.flush()
                if snapshot_reference_binding is not None:
                    session.add(
                        _run_reference_binding_row(
                            record.run_id,
                            snapshot_reference_binding,
                        )
                    )
                    session.flush()
                created_event = self._insert_event(
                    session,
                    record.run_id,
                    event,
                )
                if created_event.sequence != 1:
                    raise ValueError("validated run initial event is invalid")
                row.consumed_run_id = record.run_id
                row.consumed_at = consumed_at
                session.flush()
                if security_audit is not None:
                    insert_security_audit_event(session, security_audit)
                return ValidatedSnapshotRunCreation(record=record, created=True)
        except IntegrityError as exc:
            raise ValueError("validated snapshot could not create a run") from exc

    def get_run(self, run_id: str) -> RunRecord:
        with self._session_factory() as session:
            return _record_from_row(self._require_run(session, run_id))

    def get_run_requester_user_id(self, run_id: str) -> str | None:
        """Return the immutable create-time requester without exposing it on a run."""

        with self._session_factory() as session:
            row = session.scalar(select(RunRow).where(RunRow.run_id == run_id))
            if row is None:
                raise KeyError(run_id)
            return row.requested_by_user_id

    def get_run_data_binding(self, run_id: str) -> ProjectSampleBinding:
        with self._session_factory() as session:
            if session.scalar(select(RunRow.id).where(RunRow.run_id == run_id)) is None:
                raise KeyError(run_id)
            return _run_data_binding_from_rows(session, run_id)

    def get_run_input_use_binding(self, run_id: str) -> InputUseBindingEnvelope:
        with self._session_factory() as session:
            run_row = session.scalar(select(RunRow).where(RunRow.run_id == run_id))
            if run_row is None:
                raise KeyError(run_id)
            data_binding = _run_data_binding_from_rows(session, run_id)
            return _run_input_binding_from_rows(
                session,
                run_id,
                expected_data_binding=data_binding,
                expected_workflow_id=run_row.workflow_id,
            )

    def get_run_reference_binding(
        self,
        run_id: str,
    ) -> ReferenceProfileRevisionBinding | None:
        with self._session_factory() as session:
            if session.scalar(select(RunRow.id).where(RunRow.run_id == run_id)) is None:
                raise KeyError(run_id)
            row = session.get(RunReferenceBindingRow, run_id)
            return None if row is None else _reference_binding_from_row(row)

    def list_runs(self) -> tuple[RunRecord, ...]:
        with self._session_factory() as session:
            rows = session.scalars(select(RunRow).order_by(RunRow.id)).all()
            return tuple(_record_from_row(row) for row in rows)

    def list_run_summaries(
        self,
        *,
        after: RunHistoryCursor | None = None,
        limit: int = 50,
        workflow_id: str | None = None,
        status: RunStatus | None = None,
    ) -> tuple[RunSummary, ...]:
        """Return a bounded summary-only keyset query from SQLite."""
        if isinstance(limit, bool) or not isinstance(limit, int) or limit <= 0:
            raise ValueError("run summary limit must be a positive integer")
        if workflow_id is not None and not isinstance(workflow_id, str):
            raise ValueError("run summary workflow_id filter must be a string")
        if status is not None and not isinstance(status, RunStatus):
            raise ValueError("run summary status filter must be a RunStatus")
        if after is not None and not isinstance(after, RunHistoryCursor):
            raise ValueError("run summary cursor is invalid")
        if after is not None and (
            after.workflow_id != workflow_id or after.status is not status
        ):
            raise ValueError("run summary cursor filters do not match")

        with self._session_factory() as session:
            if after is not None:
                boundary_row = session.execute(
                    _run_summary_select().where(RunRow.run_id == after.run_id)
                ).one_or_none()
                if boundary_row is None:
                    raise KeyError(after.run_id)
                boundary = _summary_from_selected_row(boundary_row)
                if boundary.created_at != after.created_at or not _summary_matches(
                    boundary,
                    workflow_id=workflow_id,
                    status=status,
                ):
                    raise KeyError(after.run_id)

            statement = _run_summary_select()
            if workflow_id is not None:
                statement = statement.where(RunRow.workflow_id == workflow_id)
            if status is not None:
                statement = statement.where(RunRow.status == status.value)
            if after is not None:
                statement = statement.where(
                    or_(
                        RunRow.created_at < after.created_at,
                        and_(
                            RunRow.created_at == after.created_at,
                            RunRow.run_id < after.run_id,
                        ),
                    )
                )
            rows = session.execute(
                statement.order_by(
                    RunRow.created_at.desc(),
                    RunRow.run_id.desc(),
                ).limit(limit)
            ).all()
            return tuple(_summary_from_selected_row(row) for row in rows)

    def list_artifact_publications(
        self,
        *,
        filters: ArtifactPublicationFilters,
        after: ArtifactPublicationCursorPosition | None,
        limit: int,
    ) -> tuple[ArtifactPublicationSummary, ...]:
        if not isinstance(filters, ArtifactPublicationFilters):
            raise ValueError("artifact publication filters are invalid")
        if after is not None and not isinstance(
            after, ArtifactPublicationCursorPosition
        ):
            raise ValueError("artifact publication cursor position is invalid")
        if isinstance(limit, bool) or not isinstance(limit, int) or limit <= 0:
            raise ValueError("artifact publication limit must be a positive integer")

        with self._session_factory() as session:
            _begin_consistent_read(session)
            statement = _artifact_publication_select().where(
                *_artifact_publication_filter_conditions(filters)
            )
            if after is not None:
                boundary = session.execute(
                    statement.where(
                        ArtifactPublicationRow.run_id == after.run_id,
                        ArtifactPublicationRow.artifact_id == after.artifact_id,
                        ArtifactPublicationRow.artifact_generation
                        == after.artifact_generation,
                        ArtifactPublicationRow.published_at == after.published_at,
                    )
                ).one_or_none()
                if boundary is None:
                    raise KeyError(after.key)
                statement = statement.where(
                    _artifact_publication_after_condition(after)
                )
            rows = session.execute(
                statement.order_by(
                    ArtifactPublicationRow.published_at.desc(),
                    ArtifactPublicationRow.run_id.desc(),
                    ArtifactPublicationRow.artifact_generation.desc(),
                    ArtifactPublicationRow.artifact_id.desc(),
                ).limit(limit)
            ).all()
            return _artifact_publication_summaries_from_rows(session, rows)

    def get_artifact_publication(
        self,
        *,
        run_id: str,
        artifact_id: str,
        artifact_generation: str,
    ) -> ArtifactPublicationSummary:
        validate_artifact_generation(artifact_generation)
        with self._session_factory() as session:
            _begin_consistent_read(session)
            row = session.execute(
                _artifact_publication_select().where(
                    ArtifactPublicationRow.run_id == run_id,
                    ArtifactPublicationRow.artifact_id == artifact_id,
                    ArtifactPublicationRow.artifact_generation == artifact_generation,
                )
            ).one_or_none()
            if row is None:
                raise KeyError((run_id, artifact_id, artifact_generation))
            return _artifact_publication_summaries_from_rows(session, (row,))[0]

    def update_run(
        self,
        record: RunRecord,
        *,
        expected_status: RunStatus,
        event: RunEventDraft,
        security_audit: SecurityAuditEvent | None = None,
    ) -> RunEvent:
        _validate_security_audit(security_audit)
        with self._lock, self._session_factory.begin() as session:
            _begin_write(session)
            result = session.execute(
                update(RunRow)
                .where(
                    RunRow.run_id == record.run_id,
                    RunRow.status == expected_status.value,
                )
                .values(**_run_values(record))
            )
            if result.rowcount != 1:
                if (
                    session.scalar(
                        select(RunRow.id).where(RunRow.run_id == record.run_id)
                    )
                    is None
                ):
                    raise KeyError(record.run_id)
                raise ConcurrentRunUpdateError(
                    f"Run {record.run_id!r} changed while it was being updated."
                )
            inserted = self._insert_event(session, record.run_id, event)
            if security_audit is not None:
                insert_security_audit_event(session, security_audit)
            return inserted

    def complete_preflight(
        self,
        record: RunRecord,
        identity: WorkflowBuildIdentity,
        *,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent:
        """Atomically persist a build identity, PLANNED run, and event."""
        if record.status is not RunStatus.PLANNED:
            raise ValueError("completed preflight record must be planned")
        if identity.workflow_id != record.workflow_id:
            raise ValueError("workflow build identity does not match the run")

        try:
            with self._lock, self._session_factory.begin() as session:
                _begin_write(session)
                if session.get(RunWorkflowBuildIdentityRow, record.run_id) is not None:
                    raise ValueError("run already has a workflow build identity")
                result = session.execute(
                    update(RunRow)
                    .where(
                        RunRow.run_id == record.run_id,
                        RunRow.status == expected_status.value,
                    )
                    .values(**_run_values(record))
                )
                if result.rowcount != 1:
                    if (
                        session.scalar(
                            select(RunRow.id).where(RunRow.run_id == record.run_id)
                        )
                        is None
                    ):
                        raise KeyError(record.run_id)
                    raise ConcurrentRunUpdateError(
                        f"Run {record.run_id!r} changed while preflight completed."
                    )
                session.add(_workflow_build_identity_row(record.run_id, identity))
                session.flush()
                return self._insert_event(session, record.run_id, event)
        except IntegrityError as exc:
            raise ValueError("workflow build identity could not be persisted") from exc

    def get_workflow_build_identity(
        self,
        run_id: str,
    ) -> WorkflowBuildIdentity | None:
        with self._session_factory() as session:
            row = session.get(RunWorkflowBuildIdentityRow, run_id)
            return _workflow_build_identity_from_row(row) if row is not None else None

    def add_event(self, run_id: str, event: RunEventDraft) -> RunEvent:
        with self._lock, self._session_factory.begin() as session:
            _begin_write(session)
            self._require_run(session, run_id)
            return self._insert_event(session, run_id, event)

    def fail_interrupted_run_if_unowned(
        self,
        record: RunRecord,
        *,
        expected_status: RunStatus,
        required_ownership: RunExecutionOwnership | None,
        event: RunEventDraft,
    ) -> bool:
        with self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, record.run_id)
            if RunStatus(current.status) is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {record.run_id!r} changed while it was being recovered."
                )
            assignment = session.get(RunExecutionAssignmentRow, record.run_id)
            if _assignment_row_has_ownership(assignment, required_ownership):
                return False
            result = session.execute(
                update(RunRow)
                .where(
                    RunRow.run_id == record.run_id,
                    RunRow.status == expected_status.value,
                )
                .values(**_run_values(record))
            )
            if result.rowcount != 1:
                raise ConcurrentRunUpdateError(
                    f"Run {record.run_id!r} changed while it was being recovered."
                )
            self._insert_event(session, record.run_id, event)
            return True

    def list_events(
        self,
        run_id: str,
        *,
        after: str | None = None,
        limit: int = 50,
    ) -> tuple[RunEvent, ...]:
        with self._session_factory() as session:
            self._require_run(session, run_id)
            statement = select(RunEventRow).where(RunEventRow.run_id == run_id)
            if after is not None:
                cursor_sequence = session.scalar(
                    select(RunEventRow.sequence).where(
                        RunEventRow.run_id == run_id,
                        RunEventRow.event_id == after,
                    )
                )
                if cursor_sequence is None:
                    raise KeyError(after)
                statement = statement.where(RunEventRow.sequence > cursor_sequence)
            rows = session.scalars(
                statement.order_by(RunEventRow.sequence).limit(limit)
            ).all()
            return tuple(_event_from_row(row) for row in rows)

    def append_log(
        self,
        run_id: str,
        stream_name: str,
        lines: Iterable[str],
    ) -> RunLogChunk:
        normalized_lines = tuple(lines)
        with self._lock, self._session_factory.begin() as session:
            _begin_write(session)
            self._require_run(session, run_id)
            sequence = self._next_log_sequence(session, run_id, stream_name)
            timestamp = datetime.now(timezone.utc)
            row = RunLogRow(
                chunk_id=f"log-{sequence}",
                run_id=run_id,
                stream_name=stream_name,
                sequence=sequence,
                timestamp=timestamp,
                lines=list(normalized_lines),
            )
            session.add(row)
            session.flush()
            return RunLogChunk(
                chunk_id=row.chunk_id,
                run_id=run_id,
                stream_name=stream_name,
                sequence=sequence,
                timestamp=timestamp,
                lines=normalized_lines,
            )

    def list_logs(
        self,
        run_id: str,
        stream_name: str = "stdout",
        *,
        after: str | None = None,
        limit: int = 50,
    ) -> tuple[RunLogChunk, ...]:
        with self._session_factory() as session:
            self._require_run(session, run_id)
            statement = select(RunLogRow).where(
                RunLogRow.run_id == run_id,
                RunLogRow.stream_name == stream_name,
            )
            if after is not None:
                cursor_sequence = session.scalar(
                    select(RunLogRow.sequence).where(
                        RunLogRow.run_id == run_id,
                        RunLogRow.stream_name == stream_name,
                        RunLogRow.chunk_id == after,
                    )
                )
                if cursor_sequence is None:
                    raise KeyError(after)
                statement = statement.where(RunLogRow.sequence > cursor_sequence)
            rows = session.scalars(
                statement.order_by(RunLogRow.sequence).limit(limit)
            ).all()
            return tuple(_log_from_row(row) for row in rows)

    def get_result_state(self, run_id: str) -> RunResultState:
        with self._session_factory() as session:
            self._require_run(session, run_id)
            return _result_state_from_row(self._require_result_state(session, run_id))

    def begin_artifact_result_attempt(
        self,
        run_id: str,
        *,
        attempt_id: str,
        expected_status: RunStatus,
    ) -> RunResultState:
        validate_result_attempt_id(attempt_id)
        with self._lock, self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            if RunStatus(current.status) is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer eligible for artifact indexing."
                )
            row = self._require_result_state(session, run_id)
            state = _result_state_from_row(row)
            previous = session.get(RunResultAttemptRow, attempt_id)
            if previous is not None:
                if (
                    previous.run_id == run_id
                    and previous.result_kind == "artifact"
                    and previous.artifact_generation is None
                    and state.artifact_attempt_id == attempt_id
                ):
                    return state
                raise ConcurrentRunUpdateError(
                    "artifact result attempt was already superseded"
                )
            session.add(
                RunResultAttemptRow(
                    attempt_id=attempt_id,
                    run_id=run_id,
                    result_kind="artifact",
                    artifact_generation=None,
                )
            )
            session.flush()
            state = replace(
                state,
                artifact_attempt_id=attempt_id,
                artifact_attempt_status="pending",
            )
            _apply_result_state(row, state)
            session.flush()
            return state

    def begin_qc_result_attempt(
        self,
        run_id: str,
        *,
        attempt_id: str,
        expected_artifact_generation: str,
        expected_artifacts: tuple[RunArtifactRef, ...],
        expected_status: RunStatus,
    ) -> RunResultState:
        validate_result_attempt_id(attempt_id)
        validate_artifact_generation(expected_artifact_generation)
        with self._lock, self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            if RunStatus(current.status) is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer eligible for QC indexing."
                )
            row = self._require_result_state(session, run_id)
            state = _result_state_from_row(row)
            if state.artifact_generation != expected_artifact_generation:
                raise ResultGenerationChangedError(
                    f"Run {run_id!r} artifact generation changed before QC indexing."
                )
            current_artifacts = _sorted_artifacts(
                _artifact_from_row(artifact_row)
                for artifact_row in session.scalars(
                    select(RunArtifactRow).where(RunArtifactRow.run_id == run_id)
                ).all()
            )
            if current_artifacts != _validated_expected_artifacts(
                run_id,
                expected_artifacts,
            ):
                raise ResultGenerationChangedError(
                    f"Run {run_id!r} artifact manifest changed before QC indexing."
                )
            previous = session.get(RunResultAttemptRow, attempt_id)
            if previous is not None:
                if (
                    previous.run_id == run_id
                    and previous.result_kind == "qc"
                    and previous.artifact_generation == expected_artifact_generation
                    and state.qc_attempt_id == attempt_id
                    and state.qc_attempt_artifact_generation
                    == expected_artifact_generation
                ):
                    return state
                raise ConcurrentRunUpdateError(
                    "QC result attempt was already superseded"
                )
            session.add(
                RunResultAttemptRow(
                    attempt_id=attempt_id,
                    run_id=run_id,
                    result_kind="qc",
                    artifact_generation=expected_artifact_generation,
                )
            )
            session.flush()
            state = replace(
                state,
                qc_attempt_id=attempt_id,
                qc_attempt_status="pending",
                qc_attempt_artifact_generation=expected_artifact_generation,
            )
            _apply_result_state(row, state)
            session.flush()
            return state

    def replace_artifacts(
        self,
        run_id: str,
        artifacts: tuple[RunArtifactRef, ...],
        *,
        attempt_id: str,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent | None:
        validate_result_attempt_id(attempt_id)
        if expected_status is not RunStatus.SUCCEEDED:
            raise ConcurrentRunUpdateError(
                "artifact publication requires a succeeded run"
            )
        try:
            with self._lock, self._session_factory.begin() as session:
                _begin_write(session)
                current = self._require_run(session, run_id)
                if RunStatus(current.status) is not expected_status:
                    raise ConcurrentRunUpdateError(
                        f"Run {run_id!r} is no longer eligible for artifact replacement."
                    )
                sorted_replacement = _validated_expected_artifacts(run_id, artifacts)
                state_row = self._require_result_state(session, run_id)
                state = _result_state_from_row(state_row)
                _require_current_attempt(
                    state.artifact_attempt_id, attempt_id, "artifact"
                )
                manifest_digest = artifact_manifest_digest(sorted_replacement)
                if state.artifact_attempt_status == "succeeded":
                    if state.artifact_manifest_digest == manifest_digest:
                        return None
                    raise ConcurrentRunUpdateError("artifact attempt already committed")
                if state.artifact_attempt_status != "pending":
                    raise ConcurrentRunUpdateError(
                        "artifact attempt is no longer current"
                    )
                changed = state.artifact_manifest_digest != manifest_digest
                should_emit = changed or state.artifact_outcome != "succeeded"
                had_qc_rows = (
                    session.scalar(
                        select(func.count())
                        .select_from(RunQcMetricRow)
                        .where(RunQcMetricRow.run_id == run_id)
                    )
                    or 0
                ) > 0
                had_qc_state = (
                    had_qc_rows
                    or state.qc_generation is not None
                    or (state.qc_outcome is not None)
                )
                if changed:
                    revision = state.artifact_revision + 1
                    generation = build_artifact_generation(
                        run_id=run_id,
                        revision=revision,
                        artifacts=sorted_replacement,
                    )
                    session.execute(
                        delete(RunArtifactRow).where(RunArtifactRow.run_id == run_id)
                    )
                    session.add_all(
                        [_artifact_row(artifact) for artifact in sorted_replacement]
                    )
                    session.execute(
                        delete(RunQcMetricRow).where(RunQcMetricRow.run_id == run_id)
                    )
                    state = _state_after_artifact_change(
                        state,
                        artifacts=sorted_replacement,
                        artifact_revision=revision,
                        artifact_generation=generation,
                        attempt_id=attempt_id,
                        attempt_status="succeeded",
                        outcome="succeeded",
                    )
                    _apply_result_state(state_row, state)
                    session.flush()
                    if had_qc_state:
                        self._insert_event(
                            session,
                            run_id,
                            _event_with_context(
                                RunEventDraft(
                                    event_type="qc_metrics_invalidated",
                                    message="Workflow QC metrics invalidated.",
                                    status=RunStatus.SUCCEEDED,
                                    stage="qc_summary_indexing",
                                ),
                                artifact_generation=generation,
                                qc_generation=state.qc_generation,
                            ),
                        )
                else:
                    state = replace(
                        state,
                        artifact_attempt_status="succeeded",
                        artifact_outcome="succeeded",
                        artifact_reason_code=None,
                    )
                    _apply_result_state(state_row, state)
                    session.flush()
                if not should_emit:
                    return None
                result_event = self._insert_event(
                    session,
                    run_id,
                    _event_with_context(
                        event,
                        attempt_id=attempt_id,
                        artifact_generation=state.artifact_generation,
                        artifact_count=len(sorted_replacement),
                    ),
                )
                if changed:
                    if state.artifact_generation is None:
                        raise ValueError("artifact publication generation is missing")
                    session.add_all(
                        [
                            _artifact_publication_row(
                                artifact,
                                artifact_generation=state.artifact_generation,
                                published_at=result_event.timestamp,
                            )
                            for artifact in sorted_replacement
                        ]
                    )
                    session.flush()
                return result_event
        except IntegrityError as exc:
            raise ValueError("artifact replacement could not be persisted") from exc

    def record_artifact_failure(
        self,
        run_id: str,
        *,
        attempt_id: str,
        reason_code: str,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent | None:
        validate_result_attempt_id(attempt_id)
        with self._lock, self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            if RunStatus(current.status) is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer eligible for artifact failure."
                )
            row = self._require_result_state(session, run_id)
            state = _result_state_from_row(row)
            _require_current_attempt(state.artifact_attempt_id, attempt_id, "artifact")
            if state.artifact_attempt_status == "succeeded":
                return None
            if state.artifact_attempt_status == "failed":
                if state.artifact_reason_code == reason_code:
                    return None
                raise ConcurrentRunUpdateError("artifact attempt already failed")
            state = replace(
                state,
                artifact_attempt_status="failed",
                artifact_outcome="failed",
                artifact_reason_code=reason_code,
            )
            _apply_result_state(row, state)
            session.flush()
            return self._insert_event(
                session,
                run_id,
                _event_with_context(
                    event, attempt_id=attempt_id, reason_code=reason_code
                ),
            )

    def list_artifacts(
        self,
        run_id: str,
        *,
        after: str | None = None,
        limit: int | None = None,
    ) -> tuple[RunArtifactRef, ...]:
        _, artifacts = self.list_artifacts_page(
            run_id,
            expected_generation=None,
            after=after,
            limit=limit,
        )
        return artifacts

    def list_artifacts_page(
        self,
        run_id: str,
        *,
        expected_generation: str | None,
        after: str | None = None,
        limit: int | None = None,
    ) -> tuple[str | None, tuple[RunArtifactRef, ...]]:
        if expected_generation is not None:
            validate_artifact_generation(expected_generation)
        with self._session_factory() as session:
            _begin_consistent_read(session)
            self._require_run(session, run_id)
            state = _result_state_from_row(self._require_result_state(session, run_id))
            if (
                expected_generation is not None
                and state.artifact_generation != expected_generation
            ):
                raise ResultGenerationChangedError("artifact generation changed")
            has_rows = session.scalar(
                select(func.count())
                .select_from(RunArtifactRow)
                .where(RunArtifactRow.run_id == run_id)
            )
            if has_rows and state.artifact_generation is None:
                raise ValueError("artifact generation is unbound")
            if after is not None:
                cursor = session.scalar(
                    select(RunArtifactRow.artifact_id).where(
                        RunArtifactRow.run_id == run_id,
                        RunArtifactRow.artifact_id == after,
                    )
                )
                if cursor is None:
                    raise KeyError((run_id, after))
            statement = (
                select(RunArtifactRow)
                .where(
                    RunArtifactRow.run_id == run_id,
                    *(
                        (RunArtifactRow.artifact_id > after,)
                        if after is not None
                        else ()
                    ),
                )
                .order_by(RunArtifactRow.artifact_id)
            )
            if limit is not None:
                statement = statement.limit(limit)
            rows = session.scalars(statement).all()
            return state.artifact_generation, tuple(
                _artifact_from_row(row) for row in rows
            )

    def get_artifact(self, run_id: str, artifact_id: str) -> RunArtifactRef:
        _, artifact = self.get_artifact_at_generation(
            run_id,
            artifact_id,
            expected_generation=None,
        )
        return artifact

    def get_artifact_at_generation(
        self,
        run_id: str,
        artifact_id: str,
        *,
        expected_generation: str | None,
    ) -> tuple[str, RunArtifactRef]:
        if expected_generation is not None:
            validate_artifact_generation(expected_generation)
        with self._session_factory() as session:
            _begin_consistent_read(session)
            self._require_run(session, run_id)
            state = _result_state_from_row(self._require_result_state(session, run_id))
            if (
                expected_generation is not None
                and state.artifact_generation != expected_generation
            ):
                raise ResultGenerationChangedError("artifact generation changed")
            row = session.scalar(
                select(RunArtifactRow).where(
                    RunArtifactRow.run_id == run_id,
                    RunArtifactRow.artifact_id == artifact_id,
                )
            )
            if row is None:
                raise KeyError((run_id, artifact_id))
            if state.artifact_generation is None:
                raise ValueError("artifact generation is unbound")
            return state.artifact_generation, _artifact_from_row(row)

    def replace_qc_metrics(
        self,
        run_id: str,
        metrics: tuple[RunQcMetric, ...],
        *,
        attempt_id: str,
        expected_artifact_generation: str,
        expected_artifacts: tuple[RunArtifactRef, ...],
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent | None:
        validate_result_attempt_id(attempt_id)
        validate_artifact_generation(expected_artifact_generation)
        try:
            with self._lock, self._session_factory.begin() as session:
                _begin_write(session)
                current = self._require_run(session, run_id)
                if RunStatus(current.status) is not expected_status:
                    raise ConcurrentRunUpdateError(
                        f"Run {run_id!r} is no longer eligible for QC replacement."
                    )
                state_row = self._require_result_state(session, run_id)
                state = _result_state_from_row(state_row)
                _require_current_attempt(state.qc_attempt_id, attempt_id, "QC")
                if (
                    state.artifact_generation != expected_artifact_generation
                    or state.qc_attempt_artifact_generation
                    != expected_artifact_generation
                ):
                    raise ResultGenerationChangedError(
                        f"Run {run_id!r} artifact generation changed during QC indexing."
                    )
                current_artifacts = _sorted_artifacts(
                    _artifact_from_row(row)
                    for row in session.scalars(
                        select(RunArtifactRow).where(RunArtifactRow.run_id == run_id)
                    ).all()
                )
                if current_artifacts != _validated_expected_artifacts(
                    run_id,
                    expected_artifacts,
                ):
                    raise ResultGenerationChangedError(
                        f"Run {run_id!r} artifact generation changed during QC indexing."
                    )
                replacement = _validated_qc_replacement(
                    run_id,
                    metrics,
                    {artifact.artifact_id for artifact in current_artifacts},
                )
                manifest_digest = qc_metric_manifest_digest(replacement)
                if state.qc_attempt_status == "succeeded":
                    if (
                        state.qc_manifest_digest == manifest_digest
                        and state.qc_artifact_generation == expected_artifact_generation
                    ):
                        return None
                    raise ConcurrentRunUpdateError("QC attempt already committed")
                if state.qc_attempt_status != "pending":
                    raise ConcurrentRunUpdateError("QC attempt is no longer current")
                changed = (
                    state.qc_manifest_digest != manifest_digest
                    or state.qc_artifact_generation != expected_artifact_generation
                    or state.qc_outcome != "succeeded"
                )
                if not changed:
                    state = replace(state, qc_attempt_status="succeeded")
                    _apply_result_state(state_row, state)
                    session.flush()
                    return None
                revision = state.qc_revision + 1
                generation = build_qc_generation(
                    run_id=run_id,
                    revision=revision,
                    artifact_generation=expected_artifact_generation,
                    metrics=replacement,
                )
                session.execute(
                    delete(RunQcMetricRow).where(RunQcMetricRow.run_id == run_id)
                )
                session.add_all([_qc_metric_row(metric) for metric in replacement])
                state = replace(
                    state,
                    qc_revision=revision,
                    qc_generation=generation,
                    qc_manifest_digest=manifest_digest,
                    qc_attempt_status="succeeded",
                    qc_artifact_generation=expected_artifact_generation,
                    qc_outcome="succeeded",
                    qc_reason_code=None,
                )
                _apply_result_state(state_row, state)
                session.flush()
                return self._insert_event(
                    session,
                    run_id,
                    _event_with_context(
                        event,
                        attempt_id=attempt_id,
                        artifact_generation=expected_artifact_generation,
                        qc_generation=generation,
                        metric_count=len(replacement),
                    ),
                )
        except IntegrityError as exc:
            raise ValueError("QC replacement could not be persisted") from exc

    def list_qc_metrics(
        self,
        run_id: str,
        *,
        after: str | None = None,
        limit: int | None = None,
    ) -> tuple[RunQcMetric, ...]:
        _, metrics = self.list_qc_metrics_page(
            run_id,
            expected_generation=None,
            after=after,
            limit=limit,
        )
        return metrics

    def list_qc_metrics_page(
        self,
        run_id: str,
        *,
        expected_generation: str | None,
        after: str | None = None,
        limit: int | None = None,
    ) -> tuple[str | None, tuple[RunQcMetric, ...]]:
        if expected_generation is not None:
            validate_qc_generation(expected_generation)
        with self._session_factory() as session:
            _begin_consistent_read(session)
            self._require_run(session, run_id)
            state = _result_state_from_row(self._require_result_state(session, run_id))
            if (
                expected_generation is not None
                and state.qc_generation != expected_generation
            ):
                raise ResultGenerationChangedError("QC generation changed")
            has_rows = session.scalar(
                select(func.count())
                .select_from(RunQcMetricRow)
                .where(RunQcMetricRow.run_id == run_id)
            )
            if has_rows and (
                state.qc_generation is None
                or state.qc_artifact_generation != state.artifact_generation
            ):
                raise ValueError("QC generation is unbound")
            if after is not None:
                cursor_row = session.scalar(
                    select(RunQcMetricRow).where(
                        RunQcMetricRow.run_id == run_id,
                        RunQcMetricRow.metric_id == after,
                    )
                )
                if cursor_row is None:
                    raise KeyError((run_id, after))
                _qc_metric_from_row(cursor_row)
            statement = (
                select(RunQcMetricRow)
                .where(
                    RunQcMetricRow.run_id == run_id,
                    *((RunQcMetricRow.metric_id > after,) if after is not None else ()),
                )
                .order_by(RunQcMetricRow.metric_id)
            )
            if limit is not None:
                statement = statement.limit(limit)
            rows = session.scalars(statement).all()
            return state.qc_generation, tuple(_qc_metric_from_row(row) for row in rows)

    def list_qc_metrics_for_terminal_notification(
        self,
        run_id: str,
        *,
        expected_generation: str,
        limit: int,
    ) -> tuple[int, tuple[RunQcMetric, ...]]:
        """Read one bounded, stable projection of the current successful QC set."""

        validate_qc_generation(expected_generation)
        if isinstance(limit, bool) or not isinstance(limit, int) or limit <= 0:
            raise ValueError("terminal notification QC limit must be positive")
        with self._session_factory() as session:
            _begin_consistent_read(session)
            self._require_run(session, run_id)
            state = _result_state_from_row(self._require_result_state(session, run_id))
            if (
                state.qc_generation != expected_generation
                or state.qc_outcome != "succeeded"
            ):
                raise ResultGenerationChangedError("QC generation changed")
            total_count = session.scalar(
                select(func.count())
                .select_from(RunQcMetricRow)
                .where(RunQcMetricRow.run_id == run_id)
            )
            flag_order = case(
                (RunQcMetricRow.qc_flag == "fail", 0),
                (RunQcMetricRow.qc_flag == "warning", 1),
                (RunQcMetricRow.qc_flag == "pass", 2),
                else_=3,
            )
            scope_order = case(
                (RunQcMetricRow.scope == "run", 0),
                (RunQcMetricRow.scope == "experiment", 1),
                (RunQcMetricRow.scope == "sample", 2),
                else_=3,
            )
            rows = session.scalars(
                select(RunQcMetricRow)
                .where(RunQcMetricRow.run_id == run_id)
                .order_by(
                    flag_order,
                    scope_order,
                    func.coalesce(
                        RunQcMetricRow.experiment_id,
                        RunQcMetricRow.sample_id,
                        "",
                    ),
                    func.lower(RunQcMetricRow.display_name),
                    RunQcMetricRow.metric_key,
                    RunQcMetricRow.metric_id,
                )
                .limit(limit)
            ).all()
            return int(total_count or 0), tuple(
                _qc_metric_from_row(row) for row in rows
            )

    def record_qc_metrics_failure(
        self,
        run_id: str,
        *,
        attempt_id: str,
        expected_artifact_generation: str,
        reason_code: str,
        expected_status: RunStatus,
        event: RunEventDraft,
    ) -> RunEvent | None:
        validate_result_attempt_id(attempt_id)
        validate_artifact_generation(expected_artifact_generation)
        with self._lock, self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            if RunStatus(current.status) is not expected_status:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer eligible for QC failure."
                )
            row = self._require_result_state(session, run_id)
            state = _result_state_from_row(row)
            _require_current_attempt(state.qc_attempt_id, attempt_id, "QC")
            if (
                state.artifact_generation != expected_artifact_generation
                or state.qc_attempt_artifact_generation != expected_artifact_generation
            ):
                raise ResultGenerationChangedError(
                    "QC attempt artifact generation changed"
                )
            if state.qc_attempt_status == "succeeded":
                return None
            if state.qc_attempt_status == "failed":
                if state.qc_reason_code == reason_code:
                    return None
                raise ConcurrentRunUpdateError("QC attempt already failed")
            empty: tuple[RunQcMetric, ...] = ()
            revision = state.qc_revision + 1
            generation = build_qc_generation(
                run_id=run_id,
                revision=revision,
                artifact_generation=expected_artifact_generation,
                metrics=empty,
            )
            session.execute(
                delete(RunQcMetricRow).where(RunQcMetricRow.run_id == run_id)
            )
            state = replace(
                state,
                qc_revision=revision,
                qc_generation=generation,
                qc_manifest_digest=qc_metric_manifest_digest(empty),
                qc_attempt_status="failed",
                qc_artifact_generation=expected_artifact_generation,
                qc_outcome="failed",
                qc_reason_code=reason_code,
            )
            _apply_result_state(row, state)
            session.flush()
            return self._insert_event(
                session,
                run_id,
                _event_with_context(
                    event,
                    attempt_id=attempt_id,
                    artifact_generation=expected_artifact_generation,
                    qc_generation=generation,
                    reason_code=reason_code,
                ),
            )

    def ensure_execution_assignment(
        self,
        assignment: RunExecutionAssignment,
        *,
        expected_status: RunStatus,
    ) -> RunExecutionAssignment:
        try:
            with self._session_factory.begin() as session:
                _begin_write(session)
                current = self._require_run(session, assignment.run_id)
                if RunStatus(current.status) is not expected_status:
                    raise ConcurrentRunUpdateError(
                        f"Run {assignment.run_id!r} is no longer assignable."
                    )
                reference_binding = session.get(
                    RunReferenceBindingRow,
                    assignment.run_id,
                )
                if reference_binding is not None:
                    profile = session.get(
                        ReferenceProfileRow,
                        reference_binding.profile_id,
                    )
                    if (
                        profile is None
                        or profile.enabled_revision_id != reference_binding.revision_id
                    ):
                        raise ReferenceBindingSelectionError(
                            "Reference Profile selection is not eligible"
                        )
                row = RunExecutionAssignmentRow(
                    run_id=assignment.run_id,
                    job_id=assignment.job_id,
                    backend=assignment.backend,
                    queue_name=assignment.queue_name,
                    created_at=assignment.created_at,
                    dispatched_at=assignment.dispatched_at,
                    claimed_at=assignment.claimed_at,
                    cancellation_requested_at=(assignment.cancellation_requested_at),
                    cancellation_reason=assignment.cancellation_reason,
                    cancellation_acknowledged_at=(
                        assignment.cancellation_acknowledged_at
                    ),
                    requeue_requested_at=assignment.requeue_requested_at,
                    requeue_confirmed_at=assignment.requeue_confirmed_at,
                    managed_container_scope=assignment.managed_container_scope,
                    managed_container_endpoint_identity=(
                        assignment.managed_container_endpoint_identity
                    ),
                )
                session.add(row)
                session.flush()
                return _execution_assignment_from_row(row)
        except IntegrityError as exc:
            with self._session_factory() as session:
                existing = session.get(
                    RunExecutionAssignmentRow,
                    assignment.run_id,
                )
                if existing is not None:
                    return _execution_assignment_from_row(existing)
                if (
                    session.scalar(
                        select(RunRow.run_id).where(RunRow.run_id == assignment.run_id)
                    )
                    is None
                ):
                    raise KeyError(assignment.run_id) from exc
                assigned_run_id = session.scalar(
                    select(RunExecutionAssignmentRow.run_id).where(
                        RunExecutionAssignmentRow.job_id == assignment.job_id
                    )
                )
                if assigned_run_id is not None:
                    raise ValueError(
                        f"Execution job_id {assignment.job_id!r} is already "
                        f"assigned to run {assigned_run_id!r}."
                    ) from exc
            raise

    def get_execution_assignment(
        self,
        run_id: str,
    ) -> RunExecutionAssignment | None:
        with self._session_factory() as session:
            row = session.get(RunExecutionAssignmentRow, run_id)
            return _execution_assignment_from_row(row) if row is not None else None

    def bind_execution_cleanup_identity(
        self,
        run_id: str,
        *,
        expected_status: RunStatus,
        expected_assignment: RunExecutionAssignment,
        managed_container_scope: str | None,
        managed_container_endpoint_identity: str | None,
    ) -> RunExecutionAssignment:
        with self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            row = session.get(RunExecutionAssignmentRow, run_id)
            if row is None:
                raise KeyError(run_id)
            assignment = _execution_assignment_from_row(row)
            _require_recovery_row_expectation(
                run_id=run_id,
                current=current,
                assignment=assignment,
                expected_status=expected_status,
                expected_assignment=expected_assignment,
            )
            updated_assignment = _assignment_with_bound_cleanup_identity(
                assignment,
                expected_status=expected_status,
                managed_container_scope=managed_container_scope,
                managed_container_endpoint_identity=(
                    managed_container_endpoint_identity
                ),
            )
            row.managed_container_scope = updated_assignment.managed_container_scope
            row.managed_container_endpoint_identity = (
                updated_assignment.managed_container_endpoint_identity
            )
            session.flush()
            return _execution_assignment_from_row(row)

    def mark_execution_dispatched(
        self,
        run_id: str,
        *,
        job_id: str,
        dispatched_at: datetime,
        allowed_statuses: frozenset[RunStatus],
    ) -> RunExecutionAssignment:
        with self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            row = session.get(RunExecutionAssignmentRow, run_id)
            if row is None:
                raise KeyError(run_id)
            if row.job_id != job_id:
                raise ValueError("job_id does not match the execution assignment")
            if row.dispatched_at is not None:
                return _execution_assignment_from_row(row)
            if RunStatus(current.status) not in allowed_statuses:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer dispatchable."
                )
            row.dispatched_at = dispatched_at
            session.flush()
            return _execution_assignment_from_row(row)

    def queue_dispatched_run(
        self,
        record: RunRecord,
        *,
        expected_status: RunStatus,
        job_id: str,
        backend: str,
        queue_name: str,
        event: RunEventDraft,
        security_audit: SecurityAuditEvent | None = None,
    ) -> bool:
        """Queue a dispatched planned run and append exactly one event."""
        _validate_security_audit(security_audit)
        if expected_status is not RunStatus.PLANNED:
            raise ValueError("expected_status must be planned")
        if record.status is not RunStatus.QUEUED:
            raise ValueError("queued record must have queued status")

        with self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, record.run_id)
            assignment = session.get(RunExecutionAssignmentRow, record.run_id)
            if assignment is None:
                raise KeyError(record.run_id)
            _require_assignment_row_identity(
                assignment,
                job_id=job_id,
                backend=backend,
                queue_name=queue_name,
            )
            if assignment.dispatched_at is None:
                raise ValueError("execution assignment has not been dispatched")
            if RunStatus(current.status) is not expected_status:
                return False

            result = session.execute(
                update(RunRow)
                .where(
                    RunRow.run_id == record.run_id,
                    RunRow.status == expected_status.value,
                )
                .values(**_run_values(record))
            )
            if result.rowcount != 1:
                raise ConcurrentRunUpdateError(
                    f"Run {record.run_id!r} changed while it was being queued."
                )
            self._insert_event(session, record.run_id, event)
            if security_audit is not None:
                insert_security_audit_event(session, security_audit)
            return True

    def claim_execution_assignment(
        self,
        run_id: str,
        *,
        job_id: str,
        backend: str,
        queue_name: str,
        claimed_at: datetime,
        allowed_statuses: frozenset[RunStatus],
        event: RunEventDraft,
    ) -> RunExecutionClaim:
        with self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            row = session.get(RunExecutionAssignmentRow, run_id)
            if row is None:
                raise KeyError(run_id)
            _require_assignment_row_identity(
                row,
                job_id=job_id,
                backend=backend,
                queue_name=queue_name,
            )
            assignment = _execution_assignment_from_row(row)
            if assignment.claimed_at is not None:
                return RunExecutionClaim(assignment=assignment, acquired=False)

            current_status = RunStatus(current.status)
            if current_status not in allowed_statuses:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} is no longer claimable."
                )
            row.dispatched_at = row.dispatched_at or claimed_at
            row.claimed_at = claimed_at
            session.flush()
            self._insert_event(
                session,
                run_id,
                replace(event, status=current_status),
            )
            return RunExecutionClaim(
                assignment=_execution_assignment_from_row(row),
                acquired=True,
            )

    def request_execution_cancellation(
        self,
        run_id: str,
        *,
        job_id: str,
        backend: str,
        queue_name: str,
        requested_at: datetime,
        reason: str,
        event: RunEventDraft,
        security_audit: SecurityAuditEvent | None = None,
    ) -> RunExecutionCancellationRequest:
        _validate_security_audit(security_audit)
        with self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            row = session.get(RunExecutionAssignmentRow, run_id)
            if row is None:
                raise KeyError(run_id)
            _require_assignment_row_identity(
                row,
                job_id=job_id,
                backend=backend,
                queue_name=queue_name,
            )
            assignment = _execution_assignment_from_row(row)
            if assignment.claimed_at is None:
                raise ValueError("execution assignment has not been claimed")
            record = _record_from_row(current)
            if record.status is not RunStatus.RUNNING:
                if record.status.is_terminal:
                    return RunExecutionCancellationRequest(
                        assignment=assignment,
                        record=record,
                        created=False,
                    )
                raise ConcurrentRunUpdateError(f"Run {run_id!r} is no longer running.")
            if assignment.cancellation_requested_at is not None:
                return RunExecutionCancellationRequest(
                    assignment=assignment,
                    record=record,
                    created=False,
                )

            row.cancellation_requested_at = requested_at
            row.cancellation_reason = reason
            self._insert_event(session, run_id, event)
            if security_audit is not None:
                insert_security_audit_event(session, security_audit)
            session.flush()
            return RunExecutionCancellationRequest(
                assignment=_execution_assignment_from_row(row),
                record=record,
                created=True,
            )

    def acknowledge_execution_stop(
        self,
        run_id: str,
        *,
        job_id: str,
        backend: str,
        queue_name: str,
        acknowledged_at: datetime,
        cancellation_event: RunEventDraft,
        unexpected_stop_event: RunEventDraft,
    ) -> RunExecutionStopAcknowledgement:
        with self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            row = session.get(RunExecutionAssignmentRow, run_id)
            if row is None:
                raise KeyError(run_id)
            _require_assignment_row_identity(
                row,
                job_id=job_id,
                backend=backend,
                queue_name=queue_name,
            )
            assignment = _execution_assignment_from_row(row)
            record = _record_from_row(current)
            if record.status.is_terminal:
                return RunExecutionStopAcknowledgement(
                    assignment=assignment,
                    record=record,
                    transitioned=False,
                )

            if assignment.cancellation_requested_at is not None:
                if assignment.claimed_at is None:
                    raise ValueError("execution assignment has not been claimed")
                if record.status is not RunStatus.RUNNING:
                    raise ConcurrentRunUpdateError(
                        f"Run {run_id!r} is no longer running."
                    )
                reason = assignment.cancellation_reason
                assert reason is not None
                row.cancellation_acknowledged_at = acknowledged_at
                updated = replace(
                    record,
                    status=RunStatus.CANCELLED,
                    updated_at=acknowledged_at,
                    ended_at=acknowledged_at,
                    cancellation_reason=reason,
                    error=None,
                )
                context = dict(cancellation_event.context)
                context.update(
                    {
                        "previous_status": RunStatus.RUNNING.value,
                        "new_status": RunStatus.CANCELLED.value,
                        "cancellation_reason": reason,
                    }
                )
                event = replace(cancellation_event, context=context)
            else:
                if record.status not in {
                    RunStatus.PLANNED,
                    RunStatus.QUEUED,
                    RunStatus.RUNNING,
                }:
                    raise ConcurrentRunUpdateError(
                        f"Run {run_id!r} is not worker-owned."
                    )
                updated = replace(
                    record,
                    status=RunStatus.FAILED,
                    updated_at=acknowledged_at,
                    ended_at=acknowledged_at,
                    cancellation_reason=None,
                    error=unexpected_stop_event.issue,
                )
                context = dict(unexpected_stop_event.context)
                context.update(
                    {
                        "previous_status": record.status.value,
                        "new_status": RunStatus.FAILED.value,
                    }
                )
                event = replace(unexpected_stop_event, context=context)

            result = session.execute(
                update(RunRow)
                .where(
                    RunRow.run_id == run_id,
                    RunRow.status == record.status.value,
                )
                .values(**_run_values(updated))
            )
            if result.rowcount != 1:
                raise ConcurrentRunUpdateError(
                    f"Run {run_id!r} changed while stop was acknowledged."
                )
            self._insert_event(session, run_id, event)
            session.flush()
            return RunExecutionStopAcknowledgement(
                assignment=_execution_assignment_from_row(row),
                record=updated,
                transitioned=True,
            )

    def prepare_execution_requeue(
        self,
        run_id: str,
        *,
        expected_status: RunStatus,
        expected_assignment: RunExecutionAssignment,
        requested_at: datetime,
        event: RunEventDraft,
        security_audit: SecurityAuditEvent | None = None,
    ) -> RunExecutionRequeuePreparation:
        _validate_security_audit(security_audit)
        with self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            row = session.get(RunExecutionAssignmentRow, run_id)
            if row is None:
                raise KeyError(run_id)
            assignment = _execution_assignment_from_row(row)
            _require_recovery_row_expectation(
                run_id=run_id,
                current=current,
                assignment=assignment,
                expected_status=expected_status,
                expected_assignment=expected_assignment,
            )
            if assignment.requeue_requested_at is not None:
                return RunExecutionRequeuePreparation(
                    assignment=assignment,
                    created=False,
                )
            _require_sql_requeue_candidate(
                assignment,
                expected_status=expected_status,
            )
            row.requeue_requested_at = requested_at
            self._insert_event(session, run_id, event)
            if security_audit is not None:
                insert_security_audit_event(session, security_audit)
            session.flush()
            return RunExecutionRequeuePreparation(
                assignment=_execution_assignment_from_row(row),
                created=True,
            )

    def confirm_execution_requeue(
        self,
        run_id: str,
        *,
        expected_status: RunStatus,
        expected_assignment: RunExecutionAssignment,
        confirmed_at: datetime,
        event: RunEventDraft,
        security_audit: SecurityAuditEvent | None = None,
    ) -> RunExecutionAssignment:
        _validate_security_audit(security_audit)
        with self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, run_id)
            row = session.get(RunExecutionAssignmentRow, run_id)
            if row is None:
                raise KeyError(run_id)
            assignment = _execution_assignment_from_row(row)
            _require_requeue_confirmation_row_expectation(
                run_id=run_id,
                current=current,
                assignment=assignment,
                expected_status=expected_status,
                expected_assignment=expected_assignment,
            )
            if assignment.requeue_requested_at is None:
                raise ValueError("execution requeue has not been requested")
            if assignment.requeue_confirmed_at is not None:
                return assignment
            if confirmed_at < assignment.requeue_requested_at:
                raise ValueError("requeue confirmation cannot precede request")
            row.requeue_confirmed_at = confirmed_at
            self._insert_event(session, run_id, event)
            if security_audit is not None:
                insert_security_audit_event(session, security_audit)
            session.flush()
            return _execution_assignment_from_row(row)

    def fail_run_by_recovery(
        self,
        record: RunRecord,
        *,
        expected_status: RunStatus,
        expected_assignment: RunExecutionAssignment,
        event: RunEventDraft,
        security_audit: SecurityAuditEvent | None = None,
    ) -> bool:
        _validate_security_audit(security_audit)
        with self._session_factory.begin() as session:
            _begin_write(session)
            current = self._require_run(session, record.run_id)
            row = session.get(RunExecutionAssignmentRow, record.run_id)
            if row is None:
                raise KeyError(record.run_id)
            assignment = _execution_assignment_from_row(row)
            _require_recovery_row_expectation(
                run_id=record.run_id,
                current=current,
                assignment=assignment,
                expected_status=expected_status,
                expected_assignment=expected_assignment,
            )
            _require_sql_recovery_failure_record(
                _record_from_row(current),
                record,
                expected_status=expected_status,
            )
            result = session.execute(
                update(RunRow)
                .where(
                    RunRow.run_id == record.run_id,
                    RunRow.status == expected_status.value,
                )
                .values(**_run_values(record))
            )
            if result.rowcount != 1:
                raise ConcurrentRunUpdateError(
                    f"Run {record.run_id!r} changed while recovery was committed."
                )
            self._insert_event(session, record.run_id, event)
            if security_audit is not None:
                insert_security_audit_event(session, security_audit)
            return True

    @staticmethod
    def _require_run(session: Session, run_id: str) -> RunRow:
        row = session.scalar(select(RunRow).where(RunRow.run_id == run_id))
        if row is None:
            raise KeyError(run_id)
        return row

    @staticmethod
    def _require_result_state(session: Session, run_id: str) -> RunResultStateRow:
        row = session.get(RunResultStateRow, run_id)
        if row is None:
            raise ValueError("run result generation state is missing")
        return row

    @staticmethod
    def _next_log_sequence(
        session: Session,
        run_id: str,
        stream_name: str,
    ) -> int:
        latest = session.scalar(
            select(func.max(RunLogRow.sequence)).where(
                RunLogRow.run_id == run_id,
                RunLogRow.stream_name == stream_name,
            )
        )
        return (latest or 0) + 1

    @staticmethod
    def _insert_event(
        session: Session,
        run_id: str,
        draft: RunEventDraft,
    ) -> RunEvent:
        latest = session.scalar(
            select(func.max(RunEventRow.sequence)).where(RunEventRow.run_id == run_id)
        )
        sequence = (latest or 0) + 1
        timestamp = datetime.now(timezone.utc)
        row = RunEventRow(
            event_id=f"evt-{sequence}",
            run_id=run_id,
            sequence=sequence,
            event_type=draft.event_type,
            timestamp=timestamp,
            status=draft.status.value if draft.status is not None else None,
            stage=draft.stage,
            message=draft.message,
            context=dict(draft.context),
            issue=_issue_to_json(draft.issue),
        )
        session.add(row)
        session.flush()
        return RunEvent(
            event_id=row.event_id,
            run_id=run_id,
            sequence=sequence,
            event_type=draft.event_type,
            timestamp=timestamp,
            status=draft.status,
            stage=draft.stage,
            message=draft.message,
            context=draft.context,
            issue=draft.issue,
        )


def _run_row(
    record: RunRecord,
    *,
    requested_by_user_id: str | None = None,
) -> RunRow:
    return RunRow(
        run_id=record.run_id,
        requested_by_user_id=requested_by_user_id,
        **_run_values(record),
    )


def _workflow_build_identity_row(
    run_id: str,
    identity: WorkflowBuildIdentity,
) -> RunWorkflowBuildIdentityRow:
    return RunWorkflowBuildIdentityRow(
        run_id=run_id,
        workflow_id=identity.workflow_id,
        adapter_version=identity.adapter_version,
        scheme=identity.scheme,
        logical_entrypoint=identity.logical_entrypoint,
        digest=identity.digest,
        captured_at=identity.captured_at,
    )


def _validated_input_snapshot_row(
    snapshot: ValidatedInputSnapshot,
) -> ValidatedInputSnapshotRow:
    if not isinstance(snapshot, ValidatedInputSnapshot):
        raise ValueError("snapshot must be a ValidatedInputSnapshot")
    identity = snapshot.workflow_build_identity
    return ValidatedInputSnapshotRow(
        snapshot_id=snapshot.snapshot_id,
        workflow_id=snapshot.workflow_id,
        adapter_version=snapshot.adapter_version,
        schema_version=snapshot.schema_version,
        schema_dialect=snapshot.schema_dialect,
        canonical_payload=snapshot.canonical_payload,
        payload_digest_scheme=snapshot.payload_digest_scheme,
        payload_digest=snapshot.payload_digest,
        validation_outcome=snapshot.validation_outcome,
        validation_issue_codes=list(snapshot.validation_issue_codes),
        validated_at=snapshot.validated_at,
        expires_at=snapshot.expires_at,
        build_adapter_version=identity.adapter_version,
        build_scheme=identity.scheme,
        build_logical_entrypoint=identity.logical_entrypoint,
        build_digest=identity.digest,
        build_captured_at=identity.captured_at,
        consumed_run_id=snapshot.consumed_run_id,
        consumed_at=snapshot.consumed_at,
    )


def _snapshot_reference_binding_row(
    snapshot_id: str,
    binding: ReferenceProfileRevisionBinding,
) -> SnapshotReferenceBindingRow:
    if not isinstance(binding, ReferenceProfileRevisionBinding):
        raise ValueError("binding must be exact Reference Profile evidence")
    return SnapshotReferenceBindingRow(
        snapshot_id=snapshot_id,
        **_reference_binding_values(binding),
    )


def _run_reference_binding_row(
    run_id: str,
    binding: ReferenceProfileRevisionBinding,
) -> RunReferenceBindingRow:
    if not isinstance(binding, ReferenceProfileRevisionBinding):
        raise ValueError("binding must be exact Reference Profile evidence")
    return RunReferenceBindingRow(
        run_id=run_id,
        **_reference_binding_values(binding),
    )


def _reference_binding_values(
    binding: ReferenceProfileRevisionBinding,
) -> dict[str, object]:
    return {
        "profile_id": binding.profile_id,
        "revision_id": binding.revision_id,
        "workflow_id": binding.workflow_id,
        "revision_public_identity_scheme": (binding.revision_public_identity_scheme),
        "revision_public_identity_sha256": (binding.revision_public_identity_sha256),
        "adapter_contract_version": binding.adapter_contract_version,
        "adapter_identity_scheme": binding.adapter_identity_scheme,
        "adapter_identity_sha256": binding.adapter_identity_sha256,
        "binding_digest_scheme": binding.binding_digest_scheme,
        "binding_digest": binding.binding_digest,
        "bound_at": binding.bound_at,
    }


def _resolve_project_sample_selection(
    session: Session,
    selection: ProjectSampleSelection,
    *,
    workflow_inputs_digest: str,
) -> ProjectSampleBinding:
    if not isinstance(selection, ProjectSampleSelection):
        raise ValueError("selection must be a ProjectSampleSelection")
    project_row = session.get(ProjectRow, selection.project_id)
    if project_row is None:
        raise KeyError(selection.project_id)
    if ProjectKind(project_row.kind) is not ProjectKind.USER:
        raise ValueError("Legacy Project cannot receive trusted registry data")
    if project_row.archived_at is not None:
        raise ValueError("Project is archived")

    revision_rows = session.scalars(
        select(SampleRevisionRow).where(
            SampleRevisionRow.sample_revision_id.in_(selection.sample_revision_ids)
        )
    ).all()
    by_id = {row.sample_revision_id: row for row in revision_rows}
    refs: list[SampleRevisionBindingRef] = []
    for revision_id in selection.sample_revision_ids:
        revision_row = by_id.get(revision_id)
        if revision_row is None:
            raise KeyError(revision_id)
        revision = _sample_revision_from_storage_row(revision_row)
        if revision.project_id != selection.project_id:
            raise ValueError("all SampleRevisions must belong to the same Project")
        refs.append(
            SampleRevisionBindingRef(
                sample_revision_id=revision.sample_revision_id,
                payload_digest=revision.payload_digest,
            )
        )
    return build_project_sample_binding(
        project_id=selection.project_id,
        binding_mode=BindingMode.BOUND_V1,
        provenance=BindingProvenance.RESOLVED,
        workflow_inputs_digest=workflow_inputs_digest,
        sample_revisions=tuple(refs),
    )


def _add_snapshot_data_binding(
    session: Session,
    snapshot_id: str,
    binding: ProjectSampleBinding,
    *,
    created_at: datetime,
) -> None:
    if not isinstance(binding, ProjectSampleBinding):
        raise ValueError("binding must be a ProjectSampleBinding")
    session.add(
        SnapshotProjectBindingRow(
            snapshot_id=snapshot_id,
            project_id=binding.project_id,
            binding_mode=BindingMode(binding.binding_mode).value,
            provenance=BindingProvenance(binding.provenance).value,
            workflow_inputs_digest=binding.workflow_inputs_digest,
            binding_digest_scheme=binding.digest_scheme,
            binding_digest=binding.digest,
            created_at=created_at,
        )
    )
    for ordinal, ref in enumerate(binding.sample_revisions):
        session.add(
            SnapshotSampleRevisionRow(
                snapshot_id=snapshot_id,
                project_id=binding.project_id,
                sample_revision_id=ref.sample_revision_id,
                payload_digest=ref.payload_digest,
                ordinal=ordinal,
            )
        )


def _add_run_data_binding(
    session: Session,
    run_id: str,
    binding: ProjectSampleBinding,
    *,
    created_at: datetime,
) -> None:
    if not isinstance(binding, ProjectSampleBinding):
        raise ValueError("binding must be a ProjectSampleBinding")
    session.add(
        RunProjectBindingRow(
            run_id=run_id,
            project_id=binding.project_id,
            binding_mode=BindingMode(binding.binding_mode).value,
            provenance=BindingProvenance(binding.provenance).value,
            workflow_inputs_digest=binding.workflow_inputs_digest,
            binding_digest_scheme=binding.digest_scheme,
            binding_digest=binding.digest,
            created_at=created_at,
        )
    )
    for ordinal, ref in enumerate(binding.sample_revisions):
        session.add(
            RunSampleRow(
                run_id=run_id,
                project_id=binding.project_id,
                sample_revision_id=ref.sample_revision_id,
                payload_digest=ref.payload_digest,
                ordinal=ordinal,
            )
        )


def _add_snapshot_input_binding(
    session: Session,
    snapshot_id: str,
    binding: InputUseBindingEnvelope,
    *,
    created_at: datetime,
) -> None:
    if not isinstance(binding, InputUseBindingEnvelope):
        raise ValueError("binding must be an InputUseBindingEnvelope")
    session.add(
        SnapshotInputBindingRow(
            snapshot_id=snapshot_id,
            project_id=binding.project_id,
            workflow_id=binding.workflow_id,
            adapter_contract_version=binding.adapter_contract_version,
            binding_mode=InputBindingContractMode(binding.contract_mode).value,
            workflow_inputs_digest=binding.workflow_inputs_digest,
            project_sample_binding_digest=binding.project_sample_binding_digest,
            binding_digest_scheme=binding.digest_scheme,
            binding_digest=binding.digest,
            created_at=created_at,
        )
    )
    session.flush()
    for use_ordinal, input_use in enumerate(binding.input_uses):
        session.add(
            SnapshotInputUseRow(
                snapshot_id=snapshot_id,
                project_id=binding.project_id,
                ordinal=use_ordinal,
                input_use_key=input_use.key,
                occurrence=input_use.occurrence,
                capability_version=input_use.capability_version,
                closure_contract_version=input_use.closure_contract_version,
                provenance_mode=InputProvenanceMode(input_use.provenance_mode).value,
                closure_digest_scheme=input_use.closure_digest_scheme,
                closure_digest=input_use.closure_digest,
            )
        )
        session.flush()
        for member_ordinal, member in enumerate(input_use.members):
            session.add(
                SnapshotInputMemberRow(
                    snapshot_id=snapshot_id,
                    project_id=binding.project_id,
                    use_ordinal=use_ordinal,
                    member_ordinal=member_ordinal,
                    logical_member_key=member.logical_member_key,
                    input_file_id=member.input_file_id,
                    input_file_revision_id=member.input_file_revision_id,
                    revision_digest=member.revision_digest,
                    size_bytes=member.size_bytes,
                    content_sha256=member.content_sha256,
                )
            )


def _add_run_input_binding(
    session: Session,
    run_id: str,
    binding: InputUseBindingEnvelope,
    *,
    created_at: datetime,
) -> None:
    if not isinstance(binding, InputUseBindingEnvelope):
        raise ValueError("binding must be an InputUseBindingEnvelope")
    session.add(
        RunInputBindingRow(
            run_id=run_id,
            project_id=binding.project_id,
            workflow_id=binding.workflow_id,
            adapter_contract_version=binding.adapter_contract_version,
            binding_mode=InputBindingContractMode(binding.contract_mode).value,
            workflow_inputs_digest=binding.workflow_inputs_digest,
            project_sample_binding_digest=binding.project_sample_binding_digest,
            binding_digest_scheme=binding.digest_scheme,
            binding_digest=binding.digest,
            created_at=created_at,
        )
    )
    session.flush()
    for use_ordinal, input_use in enumerate(binding.input_uses):
        session.add(
            RunInputUseRow(
                run_id=run_id,
                project_id=binding.project_id,
                ordinal=use_ordinal,
                input_use_key=input_use.key,
                occurrence=input_use.occurrence,
                capability_version=input_use.capability_version,
                closure_contract_version=input_use.closure_contract_version,
                provenance_mode=InputProvenanceMode(input_use.provenance_mode).value,
                closure_digest_scheme=input_use.closure_digest_scheme,
                closure_digest=input_use.closure_digest,
            )
        )
        session.flush()
        for member_ordinal, member in enumerate(input_use.members):
            session.add(
                RunInputMemberRow(
                    run_id=run_id,
                    project_id=binding.project_id,
                    use_ordinal=use_ordinal,
                    member_ordinal=member_ordinal,
                    logical_member_key=member.logical_member_key,
                    input_file_id=member.input_file_id,
                    input_file_revision_id=member.input_file_revision_id,
                    revision_digest=member.revision_digest,
                    size_bytes=member.size_bytes,
                    content_sha256=member.content_sha256,
                )
            )


def _snapshot_data_binding_from_rows(
    session: Session,
    snapshot_id: str,
    *,
    expected_workflow_inputs_digest: str,
) -> ProjectSampleBinding:
    row = session.get(SnapshotProjectBindingRow, snapshot_id)
    if row is None:
        raise ValueError("validated snapshot binding evidence is missing")
    if row.workflow_inputs_digest != expected_workflow_inputs_digest:
        raise ValueError("validated snapshot binding input digest differs")
    associations = session.scalars(
        select(SnapshotSampleRevisionRow)
        .where(SnapshotSampleRevisionRow.snapshot_id == snapshot_id)
        .order_by(SnapshotSampleRevisionRow.ordinal)
    ).all()
    refs = _sample_binding_refs_from_rows(
        session,
        associations,
        expected_project_id=row.project_id,
    )
    return _project_sample_binding_from_values(
        project_id=row.project_id,
        binding_mode=row.binding_mode,
        provenance=row.provenance,
        workflow_inputs_digest=row.workflow_inputs_digest,
        digest_scheme=row.binding_digest_scheme,
        digest=row.binding_digest,
        sample_revisions=refs,
    )


def _run_data_binding_from_rows(
    session: Session,
    run_id: str,
) -> ProjectSampleBinding:
    row = session.get(RunProjectBindingRow, run_id)
    if row is None:
        raise ValueError("run binding evidence is missing")
    associations = session.scalars(
        select(RunSampleRow)
        .where(RunSampleRow.run_id == run_id)
        .order_by(RunSampleRow.ordinal)
    ).all()
    refs = _sample_binding_refs_from_rows(
        session,
        associations,
        expected_project_id=row.project_id,
    )
    return _project_sample_binding_from_values(
        project_id=row.project_id,
        binding_mode=row.binding_mode,
        provenance=row.provenance,
        workflow_inputs_digest=row.workflow_inputs_digest,
        digest_scheme=row.binding_digest_scheme,
        digest=row.binding_digest,
        sample_revisions=refs,
    )


def _snapshot_input_binding_from_rows(
    session: Session,
    snapshot_id: str,
    *,
    expected_data_binding: ProjectSampleBinding,
    expected_workflow_id: str,
    expected_workflow_inputs_digest: str,
) -> InputUseBindingEnvelope:
    row = session.get(SnapshotInputBindingRow, snapshot_id)
    if row is None:
        raise ValueError("validated snapshot input binding evidence is missing")
    use_rows = session.scalars(
        select(SnapshotInputUseRow)
        .where(SnapshotInputUseRow.snapshot_id == snapshot_id)
        .order_by(SnapshotInputUseRow.ordinal)
    ).all()
    member_rows = session.scalars(
        select(SnapshotInputMemberRow)
        .where(SnapshotInputMemberRow.snapshot_id == snapshot_id)
        .order_by(
            SnapshotInputMemberRow.use_ordinal,
            SnapshotInputMemberRow.member_ordinal,
        )
    ).all()
    return _input_binding_from_rows(
        project_id=row.project_id,
        workflow_id=row.workflow_id,
        adapter_contract_version=row.adapter_contract_version,
        binding_mode=row.binding_mode,
        workflow_inputs_digest=row.workflow_inputs_digest,
        project_sample_binding_digest=row.project_sample_binding_digest,
        binding_digest_scheme=row.binding_digest_scheme,
        binding_digest=row.binding_digest,
        use_rows=use_rows,
        member_rows=member_rows,
        expected_data_binding=expected_data_binding,
        expected_workflow_id=expected_workflow_id,
        expected_workflow_inputs_digest=expected_workflow_inputs_digest,
    )


def _run_input_binding_from_rows(
    session: Session,
    run_id: str,
    *,
    expected_data_binding: ProjectSampleBinding,
    expected_workflow_id: str,
) -> InputUseBindingEnvelope:
    row = session.get(RunInputBindingRow, run_id)
    if row is None:
        raise ValueError("run input binding evidence is missing")
    use_rows = session.scalars(
        select(RunInputUseRow)
        .where(RunInputUseRow.run_id == run_id)
        .order_by(RunInputUseRow.ordinal)
    ).all()
    member_rows = session.scalars(
        select(RunInputMemberRow)
        .where(RunInputMemberRow.run_id == run_id)
        .order_by(
            RunInputMemberRow.use_ordinal,
            RunInputMemberRow.member_ordinal,
        )
    ).all()
    return _input_binding_from_rows(
        project_id=row.project_id,
        workflow_id=row.workflow_id,
        adapter_contract_version=row.adapter_contract_version,
        binding_mode=row.binding_mode,
        workflow_inputs_digest=row.workflow_inputs_digest,
        project_sample_binding_digest=row.project_sample_binding_digest,
        binding_digest_scheme=row.binding_digest_scheme,
        binding_digest=row.binding_digest,
        use_rows=use_rows,
        member_rows=member_rows,
        expected_data_binding=expected_data_binding,
        expected_workflow_id=expected_workflow_id,
        expected_workflow_inputs_digest=expected_data_binding.workflow_inputs_digest,
    )


def _input_binding_from_rows(
    *,
    project_id: str,
    workflow_id: str,
    adapter_contract_version: str | None,
    binding_mode: str,
    workflow_inputs_digest: str,
    project_sample_binding_digest: str,
    binding_digest_scheme: str,
    binding_digest: str,
    use_rows: Iterable[SnapshotInputUseRow | RunInputUseRow],
    member_rows: Iterable[SnapshotInputMemberRow | RunInputMemberRow],
    expected_data_binding: ProjectSampleBinding,
    expected_workflow_id: str,
    expected_workflow_inputs_digest: str,
) -> InputUseBindingEnvelope:
    if (
        project_id != expected_data_binding.project_id
        or project_sample_binding_digest != expected_data_binding.digest
        or workflow_id != expected_workflow_id
        or workflow_inputs_digest != expected_workflow_inputs_digest
    ):
        raise ValueError("input binding scope differs from Project/Sample evidence")

    normalized_uses = tuple(use_rows)
    normalized_members = tuple(member_rows)
    members_by_use: dict[int, list[InputFileRevisionBindingRef]] = {}
    next_member_ordinal: dict[int, int] = {}
    for member_row in normalized_members:
        expected_member_ordinal = next_member_ordinal.get(member_row.use_ordinal, 0)
        if member_row.member_ordinal != expected_member_ordinal:
            raise ValueError("input binding member order is invalid")
        next_member_ordinal[member_row.use_ordinal] = expected_member_ordinal + 1
        members_by_use.setdefault(member_row.use_ordinal, []).append(
            InputFileRevisionBindingRef(
                logical_member_key=member_row.logical_member_key,
                input_file_id=member_row.input_file_id,
                input_file_revision_id=member_row.input_file_revision_id,
                revision_digest=member_row.revision_digest,
                size_bytes=member_row.size_bytes,
                content_sha256=member_row.content_sha256,
            )
        )

    input_uses: list[InputUseBinding] = []
    for expected_ordinal, use_row in enumerate(normalized_uses):
        if use_row.ordinal != expected_ordinal:
            raise ValueError("input binding use order is invalid")
        if use_row.project_id != project_id:
            raise ValueError("input binding use Project differs")
        input_uses.append(
            InputUseBinding(
                key=use_row.input_use_key,
                occurrence=use_row.occurrence,
                capability_version=use_row.capability_version,
                closure_contract_version=use_row.closure_contract_version,
                provenance_mode=use_row.provenance_mode,
                members=tuple(members_by_use.pop(expected_ordinal, ())),
                closure_digest_scheme=use_row.closure_digest_scheme,
                closure_digest=use_row.closure_digest,
            )
        )
    if members_by_use:
        raise ValueError("input binding contains members for an unknown use")
    return InputUseBindingEnvelope(
        project_id=project_id,
        project_sample_binding_digest=project_sample_binding_digest,
        workflow_id=workflow_id,
        adapter_contract_version=adapter_contract_version,
        workflow_inputs_digest=workflow_inputs_digest,
        contract_mode=binding_mode,
        input_uses=tuple(input_uses),
        digest_scheme=binding_digest_scheme,
        digest=binding_digest,
    )


def _sample_binding_refs_from_rows(
    session: Session,
    rows: Iterable[SnapshotSampleRevisionRow | RunSampleRow],
    *,
    expected_project_id: str,
) -> tuple[SampleRevisionBindingRef, ...]:
    refs: list[SampleRevisionBindingRef] = []
    for expected_ordinal, row in enumerate(rows):
        if row.ordinal != expected_ordinal:
            raise ValueError("binding sample revision order is invalid")
        if row.project_id != expected_project_id:
            raise ValueError("binding sample revision Project differs")
        revision_row = session.get(SampleRevisionRow, row.sample_revision_id)
        if revision_row is None:
            raise ValueError("binding SampleRevision evidence is missing")
        revision = _sample_revision_from_storage_row(revision_row)
        if (
            revision.project_id != expected_project_id
            or revision.payload_digest != row.payload_digest
        ):
            raise ValueError("binding SampleRevision evidence differs")
        refs.append(
            SampleRevisionBindingRef(
                sample_revision_id=row.sample_revision_id,
                payload_digest=row.payload_digest,
            )
        )
    return tuple(refs)


def _sample_revision_from_storage_row(
    row: SampleRevisionRow,
) -> SampleRevision:
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


def _project_sample_binding_from_values(
    *,
    project_id: str,
    binding_mode: str,
    provenance: str,
    workflow_inputs_digest: str,
    digest_scheme: str,
    digest: str,
    sample_revisions: tuple[SampleRevisionBindingRef, ...],
) -> ProjectSampleBinding:
    return ProjectSampleBinding(
        project_id=project_id,
        binding_mode=binding_mode,
        provenance=provenance,
        workflow_inputs_digest=workflow_inputs_digest,
        sample_revisions=sample_revisions,
        digest_scheme=digest_scheme,
        digest=digest,
    )


def _validate_security_audit(event: SecurityAuditEvent | None) -> None:
    if event is not None and not isinstance(event, SecurityAuditEvent):
        raise ValueError("security_audit must be a SecurityAuditEvent or None")


def _begin_write(session: Session) -> None:
    """Serialize SQLite writers before read-then-increment sequence allocation."""
    bind = session.get_bind()
    if bind.dialect.name == "sqlite":
        session.connection().exec_driver_sql("BEGIN IMMEDIATE")


def _begin_consistent_read(session: Session) -> None:
    """Establish one SQLite snapshot before a generation-bound multi-read."""
    bind = session.get_bind()
    if bind.dialect.name == "sqlite":
        # Python sqlite3 legacy transaction control does not emit BEGIN for
        # SELECT. Without an explicit read transaction, a writer can commit
        # between the generation and row queries and produce a mixed page.
        session.connection().exec_driver_sql("BEGIN")


def _run_values(record: RunRecord) -> dict[str, object]:
    return {
        "workflow_id": record.workflow_id,
        "inputs": dict(record.inputs),
        "status": record.status.value,
        "created_at": record.created_at,
        "updated_at": record.updated_at,
        "started_at": record.started_at,
        "ended_at": record.ended_at,
        "current_stage": record.current_stage,
        "cancellation_reason": record.cancellation_reason,
        "error": _issue_to_json(record.error),
        "tags": dict(record.tags),
    }


def _record_from_row(row: RunRow) -> RunRecord:
    return RunRecord(
        run_id=row.run_id,
        workflow_id=row.workflow_id,
        inputs=row.inputs,
        status=RunStatus(row.status),
        created_at=_as_utc(row.created_at),
        updated_at=_as_utc(row.updated_at),
        started_at=_optional_utc(row.started_at),
        ended_at=_optional_utc(row.ended_at),
        current_stage=row.current_stage,
        cancellation_reason=row.cancellation_reason,
        error=_issue_from_json(row.error),
        tags=row.tags,
    )


def _run_summary_select():
    return select(
        RunRow.run_id,
        RunRow.workflow_id,
        RunRow.status,
        RunRow.created_at,
        RunRow.updated_at,
        RunRow.started_at,
        RunRow.ended_at,
        RunRow.current_stage,
    )


def _artifact_publication_select():
    return (
        select(
            ArtifactPublicationRow,
            RunRow.workflow_id,
            RunProjectBindingRow.project_id,
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
    )


def _artifact_publication_filter_conditions(
    filters: ArtifactPublicationFilters,
) -> tuple[object, ...]:
    conditions: list[object] = []
    if filters.project_id is not None:
        conditions.append(RunProjectBindingRow.project_id == filters.project_id)
    if filters.run_id is not None:
        conditions.append(ArtifactPublicationRow.run_id == filters.run_id)
    if filters.workflow_id is not None:
        conditions.append(RunRow.workflow_id == filters.workflow_id)
    if filters.output_type is not None:
        conditions.append(ArtifactPublicationRow.output_type == filters.output_type)
    if filters.associated_run_sample_revision_id is not None:
        conditions.append(
            select(1)
            .select_from(RunSampleRow)
            .where(
                RunSampleRow.run_id == ArtifactPublicationRow.run_id,
                RunSampleRow.sample_revision_id
                == filters.associated_run_sample_revision_id,
            )
            .exists()
        )
    if filters.published_from is not None:
        conditions.append(ArtifactPublicationRow.published_at >= filters.published_from)
    if filters.published_before is not None:
        conditions.append(
            ArtifactPublicationRow.published_at < filters.published_before
        )
    if filters.current_only:
        conditions.append(
            ArtifactPublicationRow.artifact_generation
            == RunResultStateRow.artifact_generation
        )
    return tuple(conditions)


def _artifact_publication_after_condition(
    after: ArtifactPublicationCursorPosition,
):
    return or_(
        ArtifactPublicationRow.published_at < after.published_at,
        and_(
            ArtifactPublicationRow.published_at == after.published_at,
            ArtifactPublicationRow.run_id < after.run_id,
        ),
        and_(
            ArtifactPublicationRow.published_at == after.published_at,
            ArtifactPublicationRow.run_id == after.run_id,
            ArtifactPublicationRow.artifact_generation < after.artifact_generation,
        ),
        and_(
            ArtifactPublicationRow.published_at == after.published_at,
            ArtifactPublicationRow.run_id == after.run_id,
            ArtifactPublicationRow.artifact_generation == after.artifact_generation,
            ArtifactPublicationRow.artifact_id < after.artifact_id,
        ),
    )


def _artifact_publication_summaries_from_rows(
    session: Session,
    rows: Iterable[object],
) -> tuple[ArtifactPublicationSummary, ...]:
    selected_rows = tuple(rows)
    if not selected_rows:
        return ()
    run_ids: set[str] = set()
    for row in selected_rows:
        try:
            publication = row[0]  # type: ignore[index]
        except (IndexError, TypeError) as exc:
            raise ValueError("persisted artifact publication row is invalid") from exc
        if not isinstance(publication, ArtifactPublicationRow):
            raise ValueError("persisted artifact publication row is invalid")
        run_ids.add(publication.run_id)

    sample_rows = session.execute(
        select(RunSampleRow, SampleRevisionRow)
        .outerjoin(
            SampleRevisionRow,
            and_(
                SampleRevisionRow.project_id == RunSampleRow.project_id,
                SampleRevisionRow.sample_revision_id == RunSampleRow.sample_revision_id,
                SampleRevisionRow.payload_digest == RunSampleRow.payload_digest,
            ),
        )
        .where(RunSampleRow.run_id.in_(run_ids))
        .order_by(RunSampleRow.run_id, RunSampleRow.ordinal)
    ).all()
    samples_by_run: dict[str, list[AssociatedRunSample]] = {}
    for run_sample, revision in sample_rows:
        if revision is None:
            raise ValueError("artifact publication SampleRevision evidence is missing")
        samples_by_run.setdefault(run_sample.run_id, []).append(
            AssociatedRunSample(
                sample_id=revision.sample_id,
                sample_revision_id=revision.sample_revision_id,
                revision_number=revision.revision_number,
                ordinal=run_sample.ordinal,
            )
        )

    summaries: list[ArtifactPublicationSummary] = []
    bindings_by_run: dict[
        str, tuple[str, str, str, ArtifactPublicationRunSampleBinding]
    ] = {}
    for row in selected_rows:
        try:
            (
                publication,
                workflow_id,
                project_id,
                binding_mode,
                provenance,
                current_artifact_generation,
            ) = row  # type: ignore[misc]
        except (TypeError, ValueError) as exc:
            raise ValueError("persisted artifact publication row is invalid") from exc
        cached = bindings_by_run.get(publication.run_id)
        if cached is None:
            run_sample_binding = ArtifactPublicationRunSampleBinding(
                binding_mode=binding_mode,
                provenance=provenance,
                associated_run_samples=tuple(
                    samples_by_run.get(publication.run_id, ())
                ),
            )
            bindings_by_run[publication.run_id] = (
                workflow_id,
                project_id,
                current_artifact_generation,
                run_sample_binding,
            )
        else:
            (
                cached_workflow_id,
                cached_project_id,
                cached_generation,
                run_sample_binding,
            ) = cached
            if (
                cached_workflow_id != workflow_id
                or cached_project_id != project_id
                or cached_generation != current_artifact_generation
            ):
                raise ValueError("artifact publication Run evidence is inconsistent")
        summaries.append(
            ArtifactPublicationSummary(
                run_id=publication.run_id,
                project_id=project_id,
                workflow_id=workflow_id,
                artifact_id=publication.artifact_id,
                output_type=publication.output_type,
                artifact_generation=publication.artifact_generation,
                artifact_revision=publication.artifact_revision,
                published_at=_as_utc(publication.published_at),
                current_artifact_generation=current_artifact_generation,
                run_sample_binding=run_sample_binding,
            )
        )
    return tuple(summaries)


def _summary_from_selected_row(row: object) -> RunSummary:
    mapping = getattr(row, "_mapping", None)
    if mapping is None:
        raise ValueError("persisted run summary row is invalid")
    return RunSummary(
        run_id=mapping["run_id"],
        workflow_id=mapping["workflow_id"],
        status=mapping["status"],
        created_at=_as_utc(mapping["created_at"]),
        updated_at=_as_utc(mapping["updated_at"]),
        started_at=_optional_utc(mapping["started_at"]),
        ended_at=_optional_utc(mapping["ended_at"]),
        current_stage=mapping["current_stage"],
    )


def _summary_matches(
    summary: RunSummary,
    *,
    workflow_id: str | None,
    status: RunStatus | None,
) -> bool:
    return (workflow_id is None or summary.workflow_id == workflow_id) and (
        status is None or summary.status is status
    )


def _validated_input_snapshot_from_row(
    row: ValidatedInputSnapshotRow,
) -> ValidatedInputSnapshot:
    issue_codes = row.validation_issue_codes
    if not isinstance(issue_codes, list):
        raise ValueError("validated snapshot issue evidence is invalid")
    identity = WorkflowBuildIdentity(
        workflow_id=row.workflow_id,
        adapter_version=row.build_adapter_version,
        scheme=row.build_scheme,
        logical_entrypoint=row.build_logical_entrypoint,
        digest=row.build_digest,
        captured_at=_as_utc(row.build_captured_at),
    )
    return ValidatedInputSnapshot(
        snapshot_id=row.snapshot_id,
        workflow_id=row.workflow_id,
        adapter_version=row.adapter_version,
        schema_version=row.schema_version,
        schema_dialect=row.schema_dialect,
        workflow_build_identity=identity,
        canonical_payload=row.canonical_payload,
        payload_digest_scheme=row.payload_digest_scheme,
        payload_digest=row.payload_digest,
        validation_outcome=row.validation_outcome,
        validation_issue_codes=tuple(issue_codes),
        validated_at=_as_utc(row.validated_at),
        expires_at=_as_utc(row.expires_at),
        consumed_run_id=row.consumed_run_id,
        consumed_at=_optional_utc(row.consumed_at),
    )


def _reference_binding_from_row(
    row: SnapshotReferenceBindingRow | RunReferenceBindingRow,
) -> ReferenceProfileRevisionBinding:
    return ReferenceProfileRevisionBinding(
        profile_id=row.profile_id,
        revision_id=row.revision_id,
        workflow_id=row.workflow_id,
        revision_public_identity_scheme=row.revision_public_identity_scheme,
        revision_public_identity_sha256=row.revision_public_identity_sha256,
        adapter_contract_version=row.adapter_contract_version,
        adapter_identity_scheme=row.adapter_identity_scheme,
        adapter_identity_sha256=row.adapter_identity_sha256,
        binding_digest_scheme=row.binding_digest_scheme,
        binding_digest=row.binding_digest,
        bound_at=_as_utc(row.bound_at),
    )


def _event_from_row(row: RunEventRow) -> RunEvent:
    return RunEvent(
        event_id=row.event_id,
        run_id=row.run_id,
        sequence=row.sequence,
        event_type=row.event_type,
        timestamp=_as_utc(row.timestamp),
        status=RunStatus(row.status) if row.status is not None else None,
        stage=row.stage,
        message=row.message,
        context=row.context,
        issue=_issue_from_json(row.issue),
    )


def _log_from_row(row: RunLogRow) -> RunLogChunk:
    return RunLogChunk(
        chunk_id=row.chunk_id,
        run_id=row.run_id,
        stream_name=row.stream_name,
        sequence=row.sequence,
        timestamp=_as_utc(row.timestamp),
        lines=tuple(row.lines),
    )


def _artifact_from_row(row: RunArtifactRow) -> RunArtifactRef:
    if row.revision is None:
        raise ValueError("persisted artifact generation is unbound")
    return RunArtifactRef(
        artifact_id=row.artifact_id,
        run_id=row.run_id,
        artifact_type=row.artifact_type,
        name=row.name,
        uri=row.uri,
        mime_type=row.mime_type,
        produced_at=_as_utc(row.produced_at),
        revision=row.revision,
        metadata=row.artifact_metadata,
    )


def _artifact_row(artifact: RunArtifactRef) -> RunArtifactRow:
    return RunArtifactRow(
        artifact_id=artifact.artifact_id,
        run_id=artifact.run_id,
        artifact_type=artifact.artifact_type,
        name=artifact.name,
        uri=artifact.uri,
        mime_type=artifact.mime_type,
        produced_at=artifact.produced_at,
        revision=artifact.revision,
        artifact_metadata=artifact.to_dict()["metadata"],
    )


def _artifact_publication_row(
    artifact: RunArtifactRef,
    *,
    artifact_generation: str,
    published_at: datetime,
) -> ArtifactPublicationRow:
    ArtifactPublicationCursorPosition(
        published_at=published_at,
        run_id=artifact.run_id,
        artifact_generation=artifact_generation,
        artifact_id=artifact.artifact_id,
    )
    return ArtifactPublicationRow(
        run_id=artifact.run_id,
        artifact_id=artifact.artifact_id,
        artifact_generation=artifact_generation,
        artifact_revision=artifact.revision,
        output_type=_artifact_publication_output_type(artifact),
        published_at=published_at,
    )


def _result_state_from_row(row: RunResultStateRow) -> RunResultState:
    return RunResultState(
        run_id=row.run_id,
        artifact_revision=row.artifact_revision,
        artifact_generation=row.artifact_generation,
        artifact_manifest_digest=row.artifact_manifest_digest,
        artifact_attempt_id=row.artifact_attempt_id,
        artifact_attempt_status=row.artifact_attempt_status,
        artifact_outcome=row.artifact_outcome,
        artifact_reason_code=row.artifact_reason_code,
        qc_revision=row.qc_revision,
        qc_generation=row.qc_generation,
        qc_manifest_digest=row.qc_manifest_digest,
        qc_attempt_id=row.qc_attempt_id,
        qc_attempt_status=row.qc_attempt_status,
        qc_attempt_artifact_generation=row.qc_attempt_artifact_generation,
        qc_artifact_generation=row.qc_artifact_generation,
        qc_outcome=row.qc_outcome,
        qc_reason_code=row.qc_reason_code,
    )


def _apply_result_state(row: RunResultStateRow, state: RunResultState) -> None:
    if row.run_id != state.run_id:
        raise ValueError("result state run_id does not match")
    row.artifact_revision = state.artifact_revision
    row.artifact_generation = state.artifact_generation
    row.artifact_manifest_digest = state.artifact_manifest_digest
    row.artifact_attempt_id = state.artifact_attempt_id
    row.artifact_attempt_status = state.artifact_attempt_status
    row.artifact_outcome = state.artifact_outcome
    row.artifact_reason_code = state.artifact_reason_code
    row.qc_revision = state.qc_revision
    row.qc_generation = state.qc_generation
    row.qc_manifest_digest = state.qc_manifest_digest
    row.qc_attempt_id = state.qc_attempt_id
    row.qc_attempt_status = state.qc_attempt_status
    row.qc_attempt_artifact_generation = state.qc_attempt_artifact_generation
    row.qc_artifact_generation = state.qc_artifact_generation
    row.qc_outcome = state.qc_outcome
    row.qc_reason_code = state.qc_reason_code


def _qc_metric_row(metric: RunQcMetric) -> RunQcMetricRow:
    return RunQcMetricRow(
        metric_id=metric.metric_id,
        run_id=metric.run_id,
        metric_key=metric.metric_key,
        display_name=metric.display_name,
        value_text=canonical_decimal_text(metric.value),
        unit=metric.unit,
        scope=metric.scope,
        sample_id=metric.sample_id,
        experiment_id=metric.experiment_id,
        assay=metric.assay,
        qc_flag=metric.qc_flag,
        source_artifact_id=metric.source_artifact_id,
        produced_at=metric.produced_at,
    )


def _qc_metric_from_row(row: RunQcMetricRow) -> RunQcMetric:
    metric = RunQcMetric(
        metric_id=row.metric_id,
        run_id=row.run_id,
        metric_key=row.metric_key,
        display_name=row.display_name,
        value=decimal_from_canonical_text(row.value_text),
        unit=row.unit,
        scope=row.scope,
        sample_id=row.sample_id,
        experiment_id=row.experiment_id,
        assay=row.assay,
        qc_flag=row.qc_flag,
        source_artifact_id=row.source_artifact_id,
        produced_at=_as_utc(row.produced_at),
    )
    _validate_qc_metric_fields(metric)
    return metric


def _execution_assignment_from_row(
    row: RunExecutionAssignmentRow,
) -> RunExecutionAssignment:
    return RunExecutionAssignment(
        run_id=row.run_id,
        job_id=row.job_id,
        backend=row.backend,
        queue_name=row.queue_name,
        created_at=_as_utc(row.created_at),
        dispatched_at=_optional_utc(row.dispatched_at),
        claimed_at=_optional_utc(row.claimed_at),
        cancellation_requested_at=_optional_utc(row.cancellation_requested_at),
        cancellation_reason=row.cancellation_reason,
        cancellation_acknowledged_at=_optional_utc(row.cancellation_acknowledged_at),
        requeue_requested_at=_optional_utc(row.requeue_requested_at),
        requeue_confirmed_at=_optional_utc(row.requeue_confirmed_at),
        managed_container_scope=row.managed_container_scope,
        managed_container_endpoint_identity=row.managed_container_endpoint_identity,
    )


def _workflow_build_identity_from_row(
    row: RunWorkflowBuildIdentityRow,
) -> WorkflowBuildIdentity:
    return WorkflowBuildIdentity(
        workflow_id=row.workflow_id,
        adapter_version=row.adapter_version,
        scheme=row.scheme,
        logical_entrypoint=row.logical_entrypoint,
        digest=row.digest,
        captured_at=_as_utc(row.captured_at),
    )


def _require_assignment_row_identity(
    row: RunExecutionAssignmentRow,
    *,
    job_id: str,
    backend: str,
    queue_name: str,
) -> None:
    if row.job_id != job_id or row.backend != backend or row.queue_name != queue_name:
        raise ValueError("execution assignment identity does not match")


def _require_recovery_row_expectation(
    *,
    run_id: str,
    current: RunRow,
    assignment: RunExecutionAssignment,
    expected_status: RunStatus,
    expected_assignment: RunExecutionAssignment,
) -> None:
    if expected_assignment.run_id != run_id:
        raise ValueError("expected execution assignment does not match run_id")
    if (
        RunStatus(current.status) is not expected_status
        or assignment != expected_assignment
    ):
        raise ConcurrentRunUpdateError(
            f"Run {run_id!r} changed while recovery was being committed."
        )


def _require_requeue_confirmation_row_expectation(
    *,
    run_id: str,
    current: RunRow,
    assignment: RunExecutionAssignment,
    expected_status: RunStatus,
    expected_assignment: RunExecutionAssignment,
) -> None:
    if (
        expected_assignment.run_id != run_id
        or expected_assignment.requeue_requested_at is None
        or not _recovery_status_is_same_or_advanced(
            expected_status,
            RunStatus(current.status),
        )
        or not _assignment_is_confirmation_advance(
            expected_assignment,
            assignment,
        )
    ):
        raise ConcurrentRunUpdateError(
            f"Run {run_id!r} changed while recovery was being confirmed."
        )


def _recovery_status_is_same_or_advanced(
    expected: RunStatus,
    current: RunStatus,
) -> bool:
    if current is expected:
        return True
    if expected is RunStatus.QUEUED:
        return current in {
            RunStatus.RUNNING,
            RunStatus.SUCCEEDED,
            RunStatus.FAILED,
            RunStatus.CANCELLED,
        }
    if expected is RunStatus.RUNNING:
        return current.is_terminal
    return False


def _assignment_is_confirmation_advance(
    expected: RunExecutionAssignment,
    current: RunExecutionAssignment,
) -> bool:
    if (
        current.run_id != expected.run_id
        or current.job_id != expected.job_id
        or current.backend != expected.backend
        or current.queue_name != expected.queue_name
        or current.created_at != expected.created_at
        or current.dispatched_at != expected.dispatched_at
        or current.requeue_requested_at != expected.requeue_requested_at
    ):
        return False
    if expected.claimed_at is not None and current.claimed_at != expected.claimed_at:
        return False
    if expected.managed_container_scope is not None and (
        current.managed_container_scope != expected.managed_container_scope
        or current.managed_container_endpoint_identity
        != expected.managed_container_endpoint_identity
    ):
        return False
    if expected.cancellation_requested_at is not None and (
        current.cancellation_requested_at != expected.cancellation_requested_at
        or current.cancellation_reason != expected.cancellation_reason
    ):
        return False
    if (
        expected.cancellation_acknowledged_at is not None
        and current.cancellation_acknowledged_at
        != expected.cancellation_acknowledged_at
    ):
        return False
    if (
        expected.requeue_confirmed_at is not None
        and current.requeue_confirmed_at != expected.requeue_confirmed_at
    ):
        return False
    return True


def _require_sql_recovery_failure_record(
    current: RunRecord,
    failed: RunRecord,
    *,
    expected_status: RunStatus,
) -> None:
    if expected_status not in {RunStatus.QUEUED, RunStatus.RUNNING}:
        raise ValueError("only active runs can be failed by recovery")
    if (
        failed.status is not RunStatus.FAILED
        or failed.run_id != current.run_id
        or failed.workflow_id != current.workflow_id
        or failed.inputs != current.inputs
        or failed.created_at != current.created_at
        or failed.started_at != current.started_at
        or failed.current_stage != current.current_stage
        or failed.cancellation_reason != current.cancellation_reason
        or failed.tags != current.tags
        or failed.error is None
        or failed.ended_at is None
    ):
        raise ValueError("recovery failure record changes immutable run evidence")


def _require_sql_requeue_candidate(
    assignment: RunExecutionAssignment,
    *,
    expected_status: RunStatus,
) -> None:
    if expected_status is not RunStatus.QUEUED:
        raise ValueError("only queued runs can be requeued")
    if assignment.dispatched_at is None:
        raise ValueError("execution assignment has not been dispatched")
    if assignment.claimed_at is not None:
        raise ValueError("claimed execution assignments cannot be requeued")
    if assignment.cancellation_requested_at is not None:
        raise ValueError("cancelling execution assignments cannot be requeued")


def _assignment_row_has_ownership(
    row: RunExecutionAssignmentRow | None,
    required_ownership: RunExecutionOwnership | None,
) -> bool:
    if required_ownership is None:
        return False
    if row is None or row.backend != "rq":
        return False
    if required_ownership is RunExecutionOwnership.DISPATCHED:
        return row.dispatched_at is not None
    return row.claimed_at is not None


def _issue_to_json(issue: Issue | None) -> dict[str, object] | None:
    return issue.to_dict() if issue is not None else None


def _issue_from_json(payload: dict[str, object] | None) -> Issue | None:
    return Issue(**payload) if payload is not None else None


def _as_utc(value: datetime) -> datetime:
    if value.tzinfo is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _optional_utc(value: datetime | None) -> datetime | None:
    return _as_utc(value) if value is not None else None
