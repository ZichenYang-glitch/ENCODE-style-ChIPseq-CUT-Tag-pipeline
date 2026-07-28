"""SQLite integration for immutable snapshot/run Project/Sample bindings."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from itertools import count

import pytest
from sqlalchemy import text
from sqlalchemy.exc import IntegrityError

from encode_pipeline.persistence import open_run_persistence
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.data_registry import (
    LEGACY_PROJECT_ID,
    BindingMode,
    ProjectSampleSelection,
)
from encode_pipeline.platform.input_registry import (
    AdapterInputUseContract,
    InputFileRevisionSelection,
    InputProvenanceMode,
    InputUseDeclaration,
    build_input_use_binding_plan,
)
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.platform.snapshots import (
    PAYLOAD_DIGEST_SCHEME,
    VALIDATION_EVIDENCE_OUTCOME,
    ValidatedInputSnapshot,
    build_workflow_inputs_digest,
    canonical_workflow_inputs_json,
)
from encode_pipeline.services.data_registry import DataRegistryService
from encode_pipeline.services.data_registry_repositories import (
    InMemoryDataRegistryRepository,
)
from encode_pipeline.services.input_file_access import FileObservation
from encode_pipeline.services.input_registry import InputRegistryService
from encode_pipeline.services.input_registry_repositories import (
    InMemoryInputRegistryRepository,
)
from encode_pipeline.services.run_repositories import (
    InMemoryRunRepository,
    ProjectSampleSelectionError,
    RunEventDraft,
    ValidatedSnapshotReplayConflictError,
)


NOW = datetime(2026, 7, 26, 9, 0, tzinfo=timezone.utc)
WORKFLOW_ID = "workflow+" + "v" * 247


def _ids(prefix: str):
    serial = count(1)
    return lambda: f"{prefix}_{next(serial):032x}"


def _identity() -> WorkflowBuildIdentity:
    return WorkflowBuildIdentity(
        workflow_id=WORKFLOW_ID,
        adapter_version="1.0.0",
        scheme="sha256-tree-v1",
        logical_entrypoint="workflow/Snakefile",
        digest="a" * 64,
        captured_at=NOW,
    )


def _snapshot(serial: int) -> ValidatedInputSnapshot:
    payload = canonical_workflow_inputs_json(
        WorkflowInputs(
            config={"threads": 1},
            samples=[{"sample": "adapter-private"}],
        )
    )
    return ValidatedInputSnapshot(
        snapshot_id=f"vsnap_{serial:032x}",
        workflow_id=WORKFLOW_ID,
        adapter_version="1.0.0",
        schema_version="1.0.0",
        schema_dialect="https://json-schema.org/draft/2020-12/schema",
        workflow_build_identity=_identity(),
        canonical_payload=payload,
        payload_digest_scheme=PAYLOAD_DIGEST_SCHEME,
        payload_digest=build_workflow_inputs_digest(payload),
        validation_outcome=VALIDATION_EVIDENCE_OUTCOME,
        validation_issue_codes=(),
        validated_at=NOW,
        expires_at=NOW + timedelta(minutes=30),
    )


def _record(run_id: str) -> RunRecord:
    return RunRecord(
        run_id=run_id,
        workflow_id=WORKFLOW_ID,
        inputs=_snapshot(1).to_workflow_inputs().to_dict(),
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


def _event() -> RunEventDraft:
    return RunEventDraft(
        event_type="status_changed",
        message="Run created.",
        status=RunStatus.CREATED,
        context={"previous_status": None, "new_status": "created"},
    )


def _registry(persistence) -> DataRegistryService:
    return DataRegistryService(
        repository=persistence.data_registry_repository,
        project_id_factory=_ids("prj"),
        sample_id_factory=_ids("smp"),
        sample_revision_id_factory=_ids("smpr"),
        now_factory=lambda: NOW,
    )


def _input_registry(persistence) -> InputRegistryService:
    return InputRegistryService(
        repository=persistence.input_registry_repository,
        storage_pool_id_factory=_ids("stgp"),
        input_file_id_factory=_ids("inpf"),
        input_file_revision_id_factory=_ids("inpfr"),
        now_factory=lambda: NOW,
    )


def _consume(repository, snapshot_id: str, run_id: str):
    return repository.consume_validated_input_snapshot(
        snapshot_id,
        workflow_id=WORKFLOW_ID,
        expected_build_identity=_identity(),
        record=_record(run_id),
        consumed_at=NOW + timedelta(minutes=1),
        event=_event(),
    )


def test_sql_bound_snapshot_survives_restart_and_copies_exact_order(tmp_path) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    first = open_run_persistence(database_url)
    registry = _registry(first)
    project = registry.create_project("Pilot")
    sample_a = registry.import_sample(
        project.project_id,
        stable_key="sample-a",
        display_name="Sample A",
        attributes={"cohort": "pilot"},
    )
    sample_b = registry.import_sample(
        project.project_id,
        stable_key="sample-b",
        display_name="Sample B",
        attributes={"cohort": "pilot"},
    )
    selection = ProjectSampleSelection(
        project_id=project.project_id,
        sample_revision_ids=(
            sample_b.revision.sample_revision_id,
            sample_a.revision.sample_revision_id,
        ),
    )
    snapshot = _snapshot(1)
    first.repository.create_validated_input_snapshot(
        snapshot,
        project_sample_selection=selection,
    )
    expected = first.repository.get_validated_input_binding(snapshot.snapshot_id)
    first.close()

    second = open_run_persistence(database_url)
    assert second.repository.get_validated_input_binding(snapshot.snapshot_id) == (
        expected
    )
    _registry(second).archive_project(project.project_id)
    created = _consume(second.repository, snapshot.snapshot_id, "run-bound")

    assert created.created is True
    assert expected.binding_mode is BindingMode.BOUND_V1
    assert expected.sample_revision_ids == selection.sample_revision_ids
    assert second.repository.get_run_data_binding("run-bound") == expected
    replay = _consume(second.repository, snapshot.snapshot_id, "run-replay")
    assert replay.created is False
    assert replay.record.run_id == "run-bound"
    second.close()


def test_sql_managed_input_binding_survives_restart_and_input_archive(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    first = open_run_persistence(database_url)
    registry = _registry(first)
    project = registry.create_project("Pilot")
    sample = registry.import_sample(
        project.project_id,
        stable_key="sample-a",
        display_name="Sample A",
        attributes={},
    )
    input_registry = _input_registry(first)
    pool = input_registry.create_storage_pool(
        display_name="Approved ingress",
        config_key="ingress-primary",
    )
    input_registry.bind_project_storage_pool(
        project.project_id,
        pool.storage_pool_id,
    )
    registered = input_registry.register_input_file(
        project.project_id,
        stable_key="reads-r1",
        observation=FileObservation(
            relative_path="reads/r1.fastq.gz",
            size_bytes=5,
            content_sha256="a" * 64,
            path_fingerprint=((1, 2), (3, 4, 5, 1, 5, 6, 7)),
        ),
    )
    plan = build_input_use_binding_plan(
        project_id=project.project_id,
        workflow_id=WORKFLOW_ID,
        contract=AdapterInputUseContract(
            adapter_contract_version="workflow-a-inputs-v1",
            declarations=(
                InputUseDeclaration(
                    key="primary-reads",
                    occurrence=0,
                    capability_version="regular-file-v1",
                    closure_contract_version="regular_file_v1",
                    allowed_provenance_modes=(InputProvenanceMode.MANAGED_REVISION_V1,),
                ),
            ),
            allows_mixed=False,
        ),
        selections=(
            InputFileRevisionSelection(
                input_use_key="primary-reads",
                occurrence=0,
                input_file_revision_ids=(registered.revision.input_file_revision_id,),
            ),
        ),
    )
    snapshot = _snapshot(8)
    first.repository.create_validated_input_snapshot(
        snapshot,
        project_sample_selection=ProjectSampleSelection(
            project_id=project.project_id,
            sample_revision_ids=(sample.revision.sample_revision_id,),
        ),
        input_use_binding_plan=plan,
    )
    frozen = first.repository.get_validated_input_use_binding(snapshot.snapshot_id)
    first.close()

    second = open_run_persistence(database_url)
    assert (
        second.repository.get_validated_input_use_binding(snapshot.snapshot_id)
        == frozen
    )
    _input_registry(second).archive_input_file(registered.input_file.input_file_id)
    created = _consume(second.repository, snapshot.snapshot_id, "run-managed")

    assert created.created is True
    assert frozen.fully_managed is True
    assert frozen.adapter_contract_version == "workflow-a-inputs-v1"
    assert frozen.input_uses[0].members[0].input_file_revision_id == (
        registered.revision.input_file_revision_id
    )
    assert second.repository.get_run_input_use_binding("run-managed") == frozen
    second.close()


def test_sql_selection_failure_rolls_back_snapshot_and_binding(tmp_path) -> None:
    persistence = open_run_persistence(f"sqlite:///{tmp_path / 'platform.db'}")
    registry = _registry(persistence)
    owning_project = registry.create_project("Owner")
    selected = registry.import_sample(
        owning_project.project_id,
        stable_key="sample-a",
        display_name="Sample A",
        attributes={},
    )
    other_project = registry.create_project("Other")
    snapshot = _snapshot(2)

    with pytest.raises(ProjectSampleSelectionError):
        persistence.repository.create_validated_input_snapshot(
            snapshot,
            project_sample_selection=ProjectSampleSelection(
                project_id=other_project.project_id,
                sample_revision_ids=(selected.revision.sample_revision_id,),
            ),
        )

    with persistence.engine.connect() as connection:
        assert (
            connection.scalar(
                text(
                    "SELECT count(*) FROM validated_input_snapshots "
                    "WHERE snapshot_id=:snapshot_id"
                ),
                {"snapshot_id": snapshot.snapshot_id},
            )
            == 0
        )
        assert (
            connection.scalar(
                text(
                    "SELECT count(*) FROM snapshot_project_bindings "
                    "WHERE snapshot_id=:snapshot_id"
                ),
                {"snapshot_id": snapshot.snapshot_id},
            )
            == 0
        )
    persistence.close()


def test_sql_direct_run_and_multiple_runs_have_total_queryable_bindings(
    tmp_path,
) -> None:
    persistence = open_run_persistence(f"sqlite:///{tmp_path / 'platform.db'}")
    persistence.repository.create_run(_record("run-direct"), _event())
    direct = persistence.repository.get_run_data_binding("run-direct")
    assert direct.project_id == LEGACY_PROJECT_ID
    assert direct.binding_mode is BindingMode.LEGACY_V1
    direct_inputs = persistence.repository.get_run_input_use_binding("run-direct")
    assert direct_inputs.project_sample_binding_digest == direct.digest
    assert direct_inputs.contract_mode.value == "compatibility_unresolved_v1"

    registry = _registry(persistence)
    project = registry.create_project("Pilot")
    imported = registry.import_sample(
        project.project_id,
        stable_key="sample-a",
        display_name="Sample A",
        attributes={},
    )
    selection = ProjectSampleSelection(
        project_id=project.project_id,
        sample_revision_ids=(imported.revision.sample_revision_id,),
    )
    for serial, run_id in ((3, "run-a"), (4, "run-b")):
        snapshot = _snapshot(serial)
        persistence.repository.create_validated_input_snapshot(
            snapshot,
            project_sample_selection=selection,
        )
        _consume(persistence.repository, snapshot.snapshot_id, run_id)

    with persistence.engine.connect() as connection:
        assert connection.scalar(text("SELECT count(*) FROM run_project_bindings")) == 3
        assert connection.scalar(text("SELECT count(*) FROM run_input_bindings")) == 3
        assert (
            connection.scalar(
                text(
                    "SELECT count(*) FROM run_samples "
                    "WHERE sample_revision_id=:revision_id"
                ),
                {"revision_id": imported.revision.sample_revision_id},
            )
            == 2
        )
    persistence.close()


def test_sql_direct_run_delete_cascades_immutable_input_binding(tmp_path) -> None:
    persistence = open_run_persistence(f"sqlite:///{tmp_path / 'platform.db'}")
    persistence.repository.create_run(_record("run-direct"), _event())

    with pytest.raises(IntegrityError, match="Run input binding is immutable"):
        with persistence.engine.begin() as connection:
            connection.execute(
                text("DELETE FROM run_input_bindings WHERE run_id='run-direct'")
            )

    with persistence.engine.begin() as connection:
        connection.execute(text("DELETE FROM runs WHERE run_id='run-direct'"))
    with persistence.engine.connect() as connection:
        assert (
            connection.scalar(
                text(
                    "SELECT count(*) FROM run_input_bindings WHERE run_id='run-direct'"
                )
            )
            == 0
        )
    persistence.close()


def test_sql_snapshot_replay_fails_closed_on_binding_corruption(tmp_path) -> None:
    persistence = open_run_persistence(f"sqlite:///{tmp_path / 'platform.db'}")
    snapshot = _snapshot(5)
    persistence.repository.create_validated_input_snapshot(snapshot)
    _consume(persistence.repository, snapshot.snapshot_id, "run-legacy")
    with persistence.engine.begin() as connection:
        connection.execute(
            text("DROP TRIGGER IF EXISTS trg_run_input_bindings_no_update")
        )
        connection.execute(
            text(
                "UPDATE run_input_bindings SET binding_digest=:digest "
                "WHERE run_id='run-legacy'"
            ),
            {"digest": "f" * 64},
        )

    with pytest.raises(ValidatedSnapshotReplayConflictError, match="binding"):
        _consume(persistence.repository, snapshot.snapshot_id, "run-replay")

    assert [run.run_id for run in persistence.repository.list_runs()] == ["run-legacy"]
    persistence.close()


@pytest.mark.parametrize("token_length", (129, 255))
def test_contract_token_lengths_have_in_memory_and_sqlite_parity(
    tmp_path,
    token_length: int,
) -> None:
    token = "a" * token_length
    contract = AdapterInputUseContract(
        adapter_contract_version=token,
        declarations=(
            InputUseDeclaration(
                key="reference-index",
                occurrence=0,
                capability_version=token,
                closure_contract_version=token,
                allowed_provenance_modes=(
                    InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1,
                ),
            ),
        ),
        allows_mixed=False,
    )

    in_memory_data = InMemoryDataRegistryRepository()
    in_memory_registry = DataRegistryService(
        repository=in_memory_data,
        project_id_factory=_ids("prj"),
        sample_id_factory=_ids("smp"),
        sample_revision_id_factory=_ids("smpr"),
        now_factory=lambda: NOW,
    )
    in_memory_project = in_memory_registry.create_project("In-memory")
    in_memory_sample = in_memory_registry.import_sample(
        in_memory_project.project_id,
        stable_key="sample-a",
        display_name="Sample A",
        attributes={},
    )
    in_memory_plan = build_input_use_binding_plan(
        project_id=in_memory_project.project_id,
        workflow_id=WORKFLOW_ID,
        contract=contract,
        selections=(),
    )
    in_memory_inputs = InMemoryInputRegistryRepository(
        project_repository=in_memory_data,
    )
    in_memory_runs = InMemoryRunRepository(
        data_registry_repository=in_memory_data,
        input_registry_repository=in_memory_inputs,
    )
    in_memory_snapshot = _snapshot(1000 + token_length)
    in_memory_runs.create_validated_input_snapshot(
        in_memory_snapshot,
        project_sample_selection=ProjectSampleSelection(
            project_id=in_memory_project.project_id,
            sample_revision_ids=(in_memory_sample.revision.sample_revision_id,),
        ),
        input_use_binding_plan=in_memory_plan,
    )
    in_memory_binding = in_memory_runs.get_validated_input_use_binding(
        in_memory_snapshot.snapshot_id
    )
    _consume(
        in_memory_runs,
        in_memory_snapshot.snapshot_id,
        f"run-in-memory-{token_length}",
    )
    in_memory_run_binding = in_memory_runs.get_run_input_use_binding(
        f"run-in-memory-{token_length}"
    )

    sqlite = open_run_persistence(
        f"sqlite:///{tmp_path / f'platform-{token_length}.db'}"
    )
    sqlite_project = _registry(sqlite).create_project("SQLite")
    sqlite_sample = _registry(sqlite).import_sample(
        sqlite_project.project_id,
        stable_key="sample-a",
        display_name="Sample A",
        attributes={},
    )
    sqlite_plan = build_input_use_binding_plan(
        project_id=sqlite_project.project_id,
        workflow_id=WORKFLOW_ID,
        contract=contract,
        selections=(),
    )
    sqlite_snapshot = _snapshot(2000 + token_length)
    sqlite.repository.create_validated_input_snapshot(
        sqlite_snapshot,
        project_sample_selection=ProjectSampleSelection(
            project_id=sqlite_project.project_id,
            sample_revision_ids=(sqlite_sample.revision.sample_revision_id,),
        ),
        input_use_binding_plan=sqlite_plan,
    )
    sqlite_binding = sqlite.repository.get_validated_input_use_binding(
        sqlite_snapshot.snapshot_id
    )
    _consume(
        sqlite.repository,
        sqlite_snapshot.snapshot_id,
        f"run-sqlite-{token_length}",
    )
    sqlite_run_binding = sqlite.repository.get_run_input_use_binding(
        f"run-sqlite-{token_length}"
    )
    sqlite.close()

    for binding in (
        in_memory_binding,
        in_memory_run_binding,
        sqlite_binding,
        sqlite_run_binding,
    ):
        assert binding.adapter_contract_version == token
        assert binding.input_uses[0].capability_version == token
        assert binding.input_uses[0].closure_contract_version == token


@pytest.mark.parametrize(
    "field",
    (
        "adapter_contract_version",
        "capability_version",
        "closure_contract_version",
    ),
)
def test_contract_tokens_over_255_are_rejected_by_domain(field: str) -> None:
    too_long = "a" * 256
    if field == "adapter_contract_version":
        with pytest.raises(ValueError, match="adapter_contract_version"):
            AdapterInputUseContract(
                adapter_contract_version=too_long,
                declarations=(),
                allows_mixed=False,
            )
        return

    with pytest.raises(ValueError, match=field):
        InputUseDeclaration(
            key="reference-index",
            occurrence=0,
            capability_version=(too_long if field == "capability_version" else "v1"),
            closure_contract_version=(
                too_long if field == "closure_contract_version" else "v1"
            ),
            allowed_provenance_modes=(InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1,),
        )


def test_sql_selection_rejects_tampered_sample_revision_payload(tmp_path) -> None:
    persistence = open_run_persistence(f"sqlite:///{tmp_path / 'platform.db'}")
    registry = _registry(persistence)
    project = registry.create_project("Pilot")
    imported = registry.import_sample(
        project.project_id,
        stable_key="sample-a",
        display_name="Sample A",
        attributes={},
    )
    with persistence.engine.begin() as connection:
        connection.execute(
            text("DROP TRIGGER IF EXISTS trg_sample_revisions_no_update")
        )
        connection.execute(
            text(
                "UPDATE sample_revisions SET canonical_payload=:payload "
                "WHERE sample_revision_id=:revision_id"
            ),
            {
                "payload": '{"attributes":{},"display_name":"Tampered"}',
                "revision_id": imported.revision.sample_revision_id,
            },
        )

    with pytest.raises(ProjectSampleSelectionError):
        persistence.repository.create_validated_input_snapshot(
            _snapshot(6),
            project_sample_selection=ProjectSampleSelection(
                project_id=project.project_id,
                sample_revision_ids=(imported.revision.sample_revision_id,),
            ),
        )
    persistence.close()


def test_sql_bound_snapshot_read_rejects_tampered_revision_payload(tmp_path) -> None:
    persistence = open_run_persistence(f"sqlite:///{tmp_path / 'platform.db'}")
    registry = _registry(persistence)
    project = registry.create_project("Pilot")
    imported = registry.import_sample(
        project.project_id,
        stable_key="sample-a",
        display_name="Sample A",
        attributes={},
    )
    snapshot = _snapshot(7)
    persistence.repository.create_validated_input_snapshot(
        snapshot,
        project_sample_selection=ProjectSampleSelection(
            project_id=project.project_id,
            sample_revision_ids=(imported.revision.sample_revision_id,),
        ),
    )
    with persistence.engine.begin() as connection:
        connection.execute(
            text("DROP TRIGGER IF EXISTS trg_sample_revisions_no_update")
        )
        connection.execute(
            text(
                "UPDATE sample_revisions SET canonical_payload=:payload "
                "WHERE sample_revision_id=:revision_id"
            ),
            {
                "payload": '{"attributes":{},"display_name":"Tampered"}',
                "revision_id": imported.revision.sample_revision_id,
            },
        )

    with pytest.raises(ValueError, match="payload digest"):
        persistence.repository.get_validated_input_binding(snapshot.snapshot_id)
    with pytest.raises(ValueError, match="payload digest"):
        _consume(persistence.repository, snapshot.snapshot_id, "run-tampered")
    assert persistence.repository.contains_run("run-tampered") is False
    persistence.close()
