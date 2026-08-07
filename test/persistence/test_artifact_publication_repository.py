"""SQL/InMemory parity for atomic artifact publication and provenance queries."""

from __future__ import annotations

from dataclasses import dataclass, replace
from datetime import datetime, timedelta, timezone
from hashlib import sha256

import pytest
from sqlalchemy import event

import encode_pipeline.persistence.repositories as repository_module
from encode_pipeline.persistence import (
    SqlAlchemyDataRegistryRepository,
    SqlAlchemyRunRepository,
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.artifact_publications import (
    ArtifactPublicationCursorPosition,
    ArtifactPublicationFilters,
)
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.data_registry import (
    LEGACY_PROJECT_ID,
    SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
    BindingMode,
    BindingProvenance,
    Project,
    ProjectKind,
    ProjectSampleSelection,
    Sample,
    SampleRevision,
    build_sample_revision_payload_digest,
    canonical_sample_revision_payload,
)
from encode_pipeline.platform.runs import RunArtifactRef, RunRecord, RunStatus
from encode_pipeline.platform.snapshots import (
    PAYLOAD_DIGEST_SCHEME,
    VALIDATION_EVIDENCE_OUTCOME,
    ValidatedInputSnapshot,
    build_workflow_inputs_digest,
    canonical_workflow_inputs_json,
)
from encode_pipeline.services.data_registry_repositories import (
    DataRegistryRepository,
    InMemoryDataRegistryRepository,
)
from encode_pipeline.services.run_repositories import (
    ConcurrentRunUpdateError,
    InMemoryRunRepository,
    RunEventDraft,
    RunRepository,
)


NOW = datetime(2026, 8, 7, 8, 0, tzinfo=timezone.utc)
PROJECT_ID = "prj_11111111111111111111111111111111"
SAMPLE_A_ID = "smp_22222222222222222222222222222222"
SAMPLE_B_ID = "smp_33333333333333333333333333333333"
REVISION_A_ID = "smpr_44444444444444444444444444444444"
REVISION_B_ID = "smpr_55555555555555555555555555555555"
UNBOUND_REVISION_ID = "smpr_66666666666666666666666666666666"


@dataclass
class _RepositoryCase:
    repository: RunRepository
    data_registry: DataRegistryRepository
    engine: object | None = None


@pytest.fixture(params=("memory", "sql"), ids=("memory", "sql"))
def repository_case(request, tmp_path):
    if request.param == "memory":
        data_registry = InMemoryDataRegistryRepository(legacy_created_at=NOW)
        yield _RepositoryCase(
            repository=InMemoryRunRepository(
                data_registry_repository=data_registry,
            ),
            data_registry=data_registry,
        )
        return

    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    session_factory = create_session_factory(engine)
    yield _RepositoryCase(
        repository=SqlAlchemyRunRepository(session_factory),
        data_registry=SqlAlchemyDataRegistryRepository(session_factory),
        engine=engine,
    )
    engine.dispose()


def _attempt(label: str) -> str:
    return "resultattempt-" + sha256(label.encode()).hexdigest()


def _record(
    run_id: str,
    *,
    workflow_id: str = "workflow-a",
    status: RunStatus = RunStatus.SUCCEEDED,
) -> RunRecord:
    return RunRecord(
        run_id=run_id,
        workflow_id=workflow_id,
        inputs={"config": {}, "samples": None, "options": {}},
        status=status,
        created_at=NOW,
        updated_at=NOW,
        started_at=NOW if status is RunStatus.SUCCEEDED else None,
        ended_at=NOW if status.is_terminal else None,
        current_stage=None,
        cancellation_reason=None,
        error=None,
        tags={},
    )


def _event(
    event_type: str = "status_changed",
    status: RunStatus = RunStatus.CREATED,
) -> RunEventDraft:
    return RunEventDraft(
        event_type=event_type,
        message=event_type.replace("_", " "),
        status=status,
    )


def _artifact(
    run_id: str,
    artifact_id: str,
    *,
    seed: str,
    artifact_type: str = "file",
    output_type: object = "summary",
) -> RunArtifactRef:
    metadata = (
        {}
        if output_type is _MISSING
        else {
            "output_type": output_type,
            "sample_revision_id": UNBOUND_REVISION_ID,
        }
    )
    return RunArtifactRef(
        artifact_id=artifact_id,
        run_id=run_id,
        artifact_type=artifact_type,
        name=f"{artifact_id}.txt",
        uri=f"run://runs/{run_id}/artifacts/{artifact_id}",
        mime_type="text/plain",
        produced_at=NOW,
        revision="artifactrev-" + sha256(seed.encode()).hexdigest(),
        metadata=metadata,
    )


_MISSING = object()


def _publish(
    repository: RunRepository,
    run_id: str,
    artifacts: tuple[RunArtifactRef, ...],
    *,
    label: str,
    expected_status: RunStatus = RunStatus.SUCCEEDED,
):
    attempt_id = _attempt(label)
    repository.begin_artifact_result_attempt(
        run_id,
        attempt_id=attempt_id,
        expected_status=expected_status,
    )
    result = repository.replace_artifacts(
        run_id,
        artifacts,
        attempt_id=attempt_id,
        expected_status=expected_status,
        event=_event("artifacts_indexed", RunStatus.SUCCEEDED),
    )
    return attempt_id, result


def test_publication_is_atomic_append_only_and_tracks_current_generation(
    repository_case,
) -> None:
    repository = repository_case.repository
    repository.create_run(_record("run-a"), _event())
    first_a = _artifact("run-a", "artifact-a", seed="a1")
    first_b = _artifact(
        "run-a",
        "artifact-b",
        seed="b1",
        artifact_type="report",
        output_type=_MISSING,
    )

    first_attempt, first_event = _publish(
        repository,
        "run-a",
        (first_b, first_a),
        label="first",
    )
    assert first_event is not None
    first_generation = repository.get_result_state("run-a").artifact_generation
    assert first_generation is not None
    first_page = repository.list_artifact_publications(
        filters=ArtifactPublicationFilters(),
        after=None,
        limit=10,
    )
    assert {item.artifact_id for item in first_page} == {
        "artifact-a",
        "artifact-b",
    }
    assert {item.published_at for item in first_page} == {first_event.timestamp}
    assert {item.output_type for item in first_page} == {"summary", "report"}
    assert all(item.project_id == LEGACY_PROJECT_ID for item in first_page)
    assert all(
        item.run_sample_binding.binding_mode is BindingMode.LEGACY_V1
        and item.run_sample_binding.provenance is BindingProvenance.UNRESOLVED
        and item.run_sample_binding.associated_run_samples == ()
        for item in first_page
    )

    assert (
        repository.replace_artifacts(
            "run-a",
            (first_b, first_a),
            attempt_id=first_attempt,
            expected_status=RunStatus.SUCCEEDED,
            event=_event("artifacts_indexed", RunStatus.SUCCEEDED),
        )
        is None
    )
    assert (
        len(
            repository.list_artifact_publications(
                filters=ArtifactPublicationFilters(current_only=False),
                after=None,
                limit=10,
            )
        )
        == 2
    )

    second_a = _artifact("run-a", "artifact-a", seed="a2")
    _, second_event = _publish(
        repository,
        "run-a",
        (second_a,),
        label="second",
    )
    assert second_event is not None
    second_generation = repository.get_result_state("run-a").artifact_generation
    assert second_generation not in {None, first_generation}
    assert [
        item.artifact_revision
        for item in repository.list_artifact_publications(
            filters=ArtifactPublicationFilters(),
            after=None,
            limit=10,
        )
    ] == [second_a.revision]
    historical = repository.list_artifact_publications(
        filters=ArtifactPublicationFilters(current_only=False),
        after=None,
        limit=10,
    )
    assert len(historical) == 3
    assert {
        (item.artifact_generation, item.artifact_id, item.artifact_revision)
        for item in historical
    } == {
        (first_generation, "artifact-a", first_a.revision),
        (first_generation, "artifact-b", first_b.revision),
        (second_generation, "artifact-a", second_a.revision),
    }

    _, empty_event = _publish(repository, "run-a", (), label="empty")
    assert empty_event is not None
    empty_generation = repository.get_result_state("run-a").artifact_generation
    assert empty_generation not in {None, first_generation, second_generation}
    assert (
        repository.list_artifact_publications(
            filters=ArtifactPublicationFilters(),
            after=None,
            limit=10,
        )
        == ()
    )
    still_historical = repository.get_artifact_publication(
        run_id="run-a",
        artifact_id="artifact-a",
        artifact_generation=first_generation,
    )
    assert still_historical.artifact_revision == first_a.revision
    assert still_historical.current_artifact_generation == empty_generation
    assert still_historical.generation_status == "superseded"


def test_failed_partial_stale_and_non_succeeded_replacements_do_not_publish(
    repository_case,
) -> None:
    repository = repository_case.repository
    repository.create_run(_record("run-a"), _event())
    original = _artifact("run-a", "artifact-original", seed="original")
    _publish(repository, "run-a", (original,), label="original")

    bad_attempt = _attempt("bad-partial")
    repository.begin_artifact_result_attempt(
        "run-a",
        attempt_id=bad_attempt,
        expected_status=RunStatus.SUCCEEDED,
    )
    before_state = repository.get_result_state("run-a")
    before_events = repository.list_events("run-a")
    before_publications = repository.list_artifact_publications(
        filters=ArtifactPublicationFilters(current_only=False),
        after=None,
        limit=10,
    )
    with pytest.raises(ValueError, match="output_type"):
        repository.replace_artifacts(
            "run-a",
            (
                _artifact("run-a", "artifact-good", seed="good"),
                _artifact(
                    "run-a",
                    "artifact-bad",
                    seed="bad",
                    output_type="",
                ),
            ),
            attempt_id=bad_attempt,
            expected_status=RunStatus.SUCCEEDED,
            event=_event("artifacts_indexed", RunStatus.SUCCEEDED),
        )
    assert repository.list_artifacts("run-a") == (original,)
    assert repository.get_result_state("run-a") == before_state
    assert repository.list_events("run-a") == before_events
    assert (
        repository.list_artifact_publications(
            filters=ArtifactPublicationFilters(current_only=False),
            after=None,
            limit=10,
        )
        == before_publications
    )

    current_attempt = _attempt("current")
    repository.begin_artifact_result_attempt(
        "run-a",
        attempt_id=current_attempt,
        expected_status=RunStatus.SUCCEEDED,
    )
    with pytest.raises(ConcurrentRunUpdateError, match="no longer current"):
        repository.replace_artifacts(
            "run-a",
            (_artifact("run-a", "artifact-stale", seed="stale"),),
            attempt_id=bad_attempt,
            expected_status=RunStatus.SUCCEEDED,
            event=_event("artifacts_indexed", RunStatus.SUCCEEDED),
        )
    assert (
        repository.list_artifact_publications(
            filters=ArtifactPublicationFilters(current_only=False),
            after=None,
            limit=10,
        )
        == before_publications
    )

    for status in (RunStatus.FAILED, RunStatus.CANCELLED):
        run_id = f"run-{status.value}"
        repository.create_run(
            _record(run_id, status=status),
            _event("status_changed", status),
        )
        non_succeeded_attempt = _attempt(status.value)
        repository.begin_artifact_result_attempt(
            run_id,
            attempt_id=non_succeeded_attempt,
            expected_status=status,
        )
        with pytest.raises(ConcurrentRunUpdateError, match="succeeded run"):
            repository.replace_artifacts(
                run_id,
                (_artifact(run_id, f"artifact-{status.value}", seed=status.value),),
                attempt_id=non_succeeded_attempt,
                expected_status=status,
                event=_event("artifacts_indexed", status),
            )
        assert (
            repository.list_artifact_publications(
                filters=ArtifactPublicationFilters(run_id=run_id),
                after=None,
                limit=10,
            )
            == ()
        )


def test_publication_filters_cursor_and_run_sample_binding_use_relations(
    repository_case,
) -> None:
    repository = repository_case.repository
    _create_bound_succeeded_run(repository_case, "run-bound")
    artifact = _artifact("run-bound", "artifact-a", seed="bound")
    _, indexed_event = _publish(
        repository,
        "run-bound",
        (artifact,),
        label="bound",
    )
    assert indexed_event is not None
    [summary] = repository.list_artifact_publications(
        filters=ArtifactPublicationFilters(
            project_id=PROJECT_ID,
            workflow_id="workflow-bound",
            output_type="summary",
            associated_run_sample_revision_id=REVISION_B_ID,
            published_from=indexed_event.timestamp,
            published_before=indexed_event.timestamp + timedelta(microseconds=1),
        ),
        after=None,
        limit=10,
    )
    assert summary.run_id == "run-bound"
    assert [
        (item.sample_id, item.sample_revision_id, item.revision_number, item.ordinal)
        for item in summary.run_sample_binding.associated_run_samples
    ] == [
        (SAMPLE_B_ID, REVISION_B_ID, 1, 0),
        (SAMPLE_A_ID, REVISION_A_ID, 1, 1),
    ]
    assert (
        repository.list_artifact_publications(
            filters=ArtifactPublicationFilters(
                associated_run_sample_revision_id=UNBOUND_REVISION_ID,
            ),
            after=None,
            limit=10,
        )
        == ()
    )

    with pytest.raises(KeyError):
        repository.list_artifact_publications(
            filters=ArtifactPublicationFilters(output_type="other"),
            after=summary.cursor_position,
            limit=10,
        )
    with pytest.raises(KeyError):
        repository.list_artifact_publications(
            filters=ArtifactPublicationFilters(),
            after=ArtifactPublicationCursorPosition(
                published_at=summary.published_at + timedelta(microseconds=1),
                run_id=summary.run_id,
                artifact_generation=summary.artifact_generation,
                artifact_id=summary.artifact_id,
            ),
            limit=10,
        )


def test_sql_publication_insert_failure_rolls_back_generation_event_and_artifacts(
    tmp_path,
    monkeypatch,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    repository = SqlAlchemyRunRepository(create_session_factory(engine))
    repository.create_run(_record("run-a"), _event())
    original = _artifact("run-a", "artifact-original", seed="original")
    _publish(repository, "run-a", (original,), label="original")
    attempt_id = _attempt("database-failure")
    repository.begin_artifact_result_attempt(
        "run-a",
        attempt_id=attempt_id,
        expected_status=RunStatus.SUCCEEDED,
    )
    before_state = repository.get_result_state("run-a")
    before_events = repository.list_events("run-a")
    before_publications = repository.list_artifact_publications(
        filters=ArtifactPublicationFilters(current_only=False),
        after=None,
        limit=10,
    )
    original_builder = repository_module._artifact_publication_row

    def invalid_publication(*args, **kwargs):
        row = original_builder(*args, **kwargs)
        row.output_type = "bad/output"
        return row

    monkeypatch.setattr(
        repository_module,
        "_artifact_publication_row",
        invalid_publication,
    )
    with pytest.raises(ValueError, match="could not be persisted"):
        repository.replace_artifacts(
            "run-a",
            (_artifact("run-a", "artifact-new", seed="new"),),
            attempt_id=attempt_id,
            expected_status=RunStatus.SUCCEEDED,
            event=_event("artifacts_indexed", RunStatus.SUCCEEDED),
        )
    assert repository.get_result_state("run-a") == before_state
    assert repository.list_artifacts("run-a") == (original,)
    assert repository.list_events("run-a") == before_events
    assert (
        repository.list_artifact_publications(
            filters=ArtifactPublicationFilters(current_only=False),
            after=None,
            limit=10,
        )
        == before_publications
    )
    engine.dispose()


def test_sql_publication_page_uses_bounded_queries_without_n_plus_one(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    repository = SqlAlchemyRunRepository(create_session_factory(engine))
    repository.create_run(_record("run-a"), _event())
    _publish(
        repository,
        "run-a",
        tuple(
            _artifact("run-a", f"artifact-{index}", seed=str(index))
            for index in range(5)
        ),
        label="page",
    )
    statements: list[str] = []

    def capture(_connection, _cursor, statement, _parameters, _context, _many):
        if statement.lstrip().upper().startswith("SELECT"):
            statements.append(statement)

    event.listen(engine, "before_cursor_execute", capture)
    try:
        page = repository.list_artifact_publications(
            filters=ArtifactPublicationFilters(),
            after=None,
            limit=10,
        )
    finally:
        event.remove(engine, "before_cursor_execute", capture)
    assert len(page) == 5
    assert len(statements) == 2
    engine.dispose()


def _create_bound_succeeded_run(
    case: _RepositoryCase,
    run_id: str,
) -> None:
    data_registry = case.data_registry
    data_registry.create_project(
        Project(
            project_id=PROJECT_ID,
            display_name="Pilot",
            kind=ProjectKind.USER,
            created_at=NOW,
        )
    )
    for sample_id, revision_id, stable_key in (
        (SAMPLE_A_ID, REVISION_A_ID, "sample-a"),
        (SAMPLE_B_ID, REVISION_B_ID, "sample-b"),
    ):
        payload = canonical_sample_revision_payload(
            display_name=stable_key,
            attributes={},
        )
        data_registry.create_sample(
            Sample(
                sample_id=sample_id,
                project_id=PROJECT_ID,
                stable_key=stable_key,
                created_at=NOW,
            ),
            SampleRevision(
                sample_revision_id=revision_id,
                project_id=PROJECT_ID,
                sample_id=sample_id,
                revision_number=1,
                canonical_payload=payload,
                payload_digest_scheme=SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
                payload_digest=build_sample_revision_payload_digest(payload),
                created_at=NOW,
            ),
        )

    workflow_id = "workflow-bound"
    build_identity = WorkflowBuildIdentity(
        workflow_id=workflow_id,
        adapter_version="1.0.0",
        scheme="sha256-tree-v1",
        logical_entrypoint="workflow/Snakefile",
        digest="a" * 64,
        captured_at=NOW,
    )
    canonical_payload = canonical_workflow_inputs_json(
        WorkflowInputs(config={}, samples=None, options={})
    )
    snapshot = ValidatedInputSnapshot(
        snapshot_id="vsnap_77777777777777777777777777777777",
        workflow_id=workflow_id,
        adapter_version="1.0.0",
        schema_version="1.0.0",
        schema_dialect="https://json-schema.org/draft/2020-12/schema",
        workflow_build_identity=build_identity,
        canonical_payload=canonical_payload,
        payload_digest_scheme=PAYLOAD_DIGEST_SCHEME,
        payload_digest=build_workflow_inputs_digest(canonical_payload),
        validation_outcome=VALIDATION_EVIDENCE_OUTCOME,
        validation_issue_codes=(),
        validated_at=NOW,
        expires_at=NOW + timedelta(hours=1),
    )
    case.repository.create_validated_input_snapshot(
        snapshot,
        project_sample_selection=ProjectSampleSelection(
            project_id=PROJECT_ID,
            sample_revision_ids=(REVISION_B_ID, REVISION_A_ID),
        ),
    )
    created_record = RunRecord(
        run_id=run_id,
        workflow_id=workflow_id,
        inputs=snapshot.to_workflow_inputs().to_dict(),
        status=RunStatus.CREATED,
        created_at=NOW + timedelta(minutes=1),
        updated_at=NOW + timedelta(minutes=1),
        started_at=None,
        ended_at=None,
        current_stage=None,
        cancellation_reason=None,
        error=None,
        tags={},
    )
    case.repository.consume_validated_input_snapshot(
        snapshot.snapshot_id,
        workflow_id=workflow_id,
        expected_build_identity=build_identity,
        record=created_record,
        consumed_at=NOW + timedelta(minutes=1),
        event=_event(),
    )
    case.repository.update_run(
        replace(
            created_record,
            status=RunStatus.SUCCEEDED,
            updated_at=NOW + timedelta(minutes=2),
            started_at=NOW + timedelta(minutes=1),
            ended_at=NOW + timedelta(minutes=2),
        ),
        expected_status=RunStatus.CREATED,
        event=_event("status_changed", RunStatus.SUCCEEDED),
    )
