"""SQLite integration for immutable snapshot/run Project/Sample bindings."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from itertools import count

import pytest
from sqlalchemy import text

from encode_pipeline.persistence import open_run_persistence
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.data_registry import (
    LEGACY_PROJECT_ID,
    BindingMode,
    ProjectSampleSelection,
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
from encode_pipeline.services.run_repositories import (
    ProjectSampleSelectionError,
    RunEventDraft,
    ValidatedSnapshotReplayConflictError,
)


NOW = datetime(2026, 7, 26, 9, 0, tzinfo=timezone.utc)


def _ids(prefix: str):
    serial = count(1)
    return lambda: f"{prefix}_{next(serial):032x}"


def _identity() -> WorkflowBuildIdentity:
    return WorkflowBuildIdentity(
        workflow_id="workflow-a",
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
        workflow_id="workflow-a",
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
        workflow_id="workflow-a",
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


def _consume(repository, snapshot_id: str, run_id: str):
    return repository.consume_validated_input_snapshot(
        snapshot_id,
        workflow_id="workflow-a",
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


def test_sql_snapshot_replay_fails_closed_on_binding_corruption(tmp_path) -> None:
    persistence = open_run_persistence(f"sqlite:///{tmp_path / 'platform.db'}")
    snapshot = _snapshot(5)
    persistence.repository.create_validated_input_snapshot(snapshot)
    _consume(persistence.repository, snapshot.snapshot_id, "run-legacy")
    with persistence.engine.begin() as connection:
        connection.execute(
            text(
                "UPDATE run_project_bindings SET binding_digest=:digest "
                "WHERE run_id='run-legacy'"
            ),
            {"digest": "f" * 64},
        )

    with pytest.raises(ValidatedSnapshotReplayConflictError, match="binding"):
        _consume(persistence.repository, snapshot.snapshot_id, "run-replay")

    assert [run.run_id for run in persistence.repository.list_runs()] == ["run-legacy"]
    persistence.close()


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
