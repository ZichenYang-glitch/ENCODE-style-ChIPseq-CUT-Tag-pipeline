"""Snapshot-to-run Project/Sample binding integration contracts."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

import pytest

from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.data_registry import (
    LEGACY_PROJECT_ID,
    BindingMode,
    BindingProvenance,
    Project,
    ProjectKind,
    ProjectSampleSelection,
    Sample,
    SampleRevision,
    SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
    build_sample_revision_payload_digest,
    canonical_sample_revision_payload,
)
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.platform.snapshots import (
    PAYLOAD_DIGEST_SCHEME,
    VALIDATION_EVIDENCE_OUTCOME,
    ValidatedInputSnapshot,
    build_workflow_inputs_digest,
    canonical_workflow_inputs_json,
)
from encode_pipeline.services.run_repositories import (
    InMemoryRunRepository,
    ProjectSampleSelectionError,
    RunEventDraft,
    ValidatedSnapshotReplayConflictError,
)
from encode_pipeline.services.data_registry_repositories import (
    InMemoryDataRegistryRepository,
)


NOW = datetime(2026, 7, 26, 8, 0, tzinfo=timezone.utc)
PROJECT_ID = "prj_11111111111111111111111111111111"
OTHER_PROJECT_ID = "prj_22222222222222222222222222222222"
SAMPLE_A_ID = "smp_33333333333333333333333333333333"
SAMPLE_B_ID = "smp_44444444444444444444444444444444"
REVISION_A_ID = "smpr_55555555555555555555555555555555"
REVISION_B_ID = "smpr_66666666666666666666666666666666"


def _identity() -> WorkflowBuildIdentity:
    return WorkflowBuildIdentity(
        workflow_id="workflow-a",
        adapter_version="1.0.0",
        scheme="sha256-tree-v1",
        logical_entrypoint="workflow/Snakefile",
        digest="a" * 64,
        captured_at=NOW,
    )


def _snapshot() -> ValidatedInputSnapshot:
    canonical_payload = canonical_workflow_inputs_json(
        WorkflowInputs(
            config={"threads": 1},
            samples=[{"sample": "adapter-private-row"}],
        )
    )
    return ValidatedInputSnapshot(
        snapshot_id="vsnap_0123456789abcdef0123456789abcdef",
        workflow_id="workflow-a",
        adapter_version="1.0.0",
        schema_version="1.0.0",
        schema_dialect="https://json-schema.org/draft/2020-12/schema",
        workflow_build_identity=_identity(),
        canonical_payload=canonical_payload,
        payload_digest_scheme=PAYLOAD_DIGEST_SCHEME,
        payload_digest=build_workflow_inputs_digest(canonical_payload),
        validation_outcome=VALIDATION_EVIDENCE_OUTCOME,
        validation_issue_codes=(),
        validated_at=NOW,
        expires_at=NOW + timedelta(minutes=30),
    )


def _record(run_id: str) -> RunRecord:
    return RunRecord(
        run_id=run_id,
        workflow_id="workflow-a",
        inputs=_snapshot().to_workflow_inputs().to_dict(),
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


def _revision(
    *,
    sample_id: str,
    revision_id: str,
    display_name: str,
) -> SampleRevision:
    payload = canonical_sample_revision_payload(
        display_name=display_name,
        attributes={"cohort": "pilot"},
    )
    return SampleRevision(
        sample_revision_id=revision_id,
        project_id=PROJECT_ID,
        sample_id=sample_id,
        revision_number=1,
        canonical_payload=payload,
        payload_digest_scheme=SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
        payload_digest=build_sample_revision_payload_digest(payload),
        created_at=NOW,
    )


def _data_registry() -> InMemoryDataRegistryRepository:
    repository = InMemoryDataRegistryRepository(legacy_created_at=NOW)
    repository.create_project(
        Project(
            project_id=PROJECT_ID,
            display_name="Pilot",
            kind=ProjectKind.USER,
            created_at=NOW,
        )
    )
    repository.create_project(
        Project(
            project_id=OTHER_PROJECT_ID,
            display_name="Other",
            kind=ProjectKind.USER,
            created_at=NOW,
        )
    )
    for sample_id, stable_key, revision_id, display_name in (
        (SAMPLE_A_ID, "sample-a", REVISION_A_ID, "Sample A"),
        (SAMPLE_B_ID, "sample-b", REVISION_B_ID, "Sample B"),
    ):
        repository.create_sample(
            Sample(
                sample_id=sample_id,
                project_id=PROJECT_ID,
                stable_key=stable_key,
                created_at=NOW,
            ),
            _revision(
                sample_id=sample_id,
                revision_id=revision_id,
                display_name=display_name,
            ),
        )
    return repository


def test_unselected_snapshot_gets_explicit_unresolved_legacy_binding() -> None:
    repository = InMemoryRunRepository()
    snapshot = _snapshot()

    repository.create_validated_input_snapshot(snapshot)

    binding = repository.get_validated_input_binding(snapshot.snapshot_id)
    assert binding.project_id == LEGACY_PROJECT_ID
    assert binding.binding_mode is BindingMode.LEGACY_V1
    assert binding.provenance is BindingProvenance.UNRESOLVED
    assert binding.workflow_inputs_digest == snapshot.payload_digest
    assert binding.sample_revisions == ()
    assert binding.input_revisions == ()


def test_snapshot_consumption_copies_exact_legacy_binding_to_run() -> None:
    repository = InMemoryRunRepository()
    snapshot = _snapshot()
    repository.create_validated_input_snapshot(snapshot)

    creation = repository.consume_validated_input_snapshot(
        snapshot.snapshot_id,
        workflow_id="workflow-a",
        expected_build_identity=_identity(),
        record=_record("run-1"),
        consumed_at=NOW + timedelta(minutes=1),
        event=_event(),
    )

    assert creation.created is True
    assert repository.get_run_data_binding("run-1") == (
        repository.get_validated_input_binding(snapshot.snapshot_id)
    )


def test_direct_compatibility_run_gets_explicit_legacy_binding() -> None:
    repository = InMemoryRunRepository()
    record = _record("run-direct")

    repository.create_run(record, _event())

    binding = repository.get_run_data_binding(record.run_id)
    assert binding.project_id == LEGACY_PROJECT_ID
    assert binding.binding_mode is BindingMode.LEGACY_V1
    assert binding.provenance is BindingProvenance.UNRESOLVED
    assert binding.sample_revisions == ()


def test_bound_snapshot_freezes_order_and_consumption_copies_exact_evidence() -> None:
    data_registry = _data_registry()
    repository = InMemoryRunRepository(data_registry_repository=data_registry)
    snapshot = _snapshot()
    selection = ProjectSampleSelection(
        project_id=PROJECT_ID,
        sample_revision_ids=(REVISION_B_ID, REVISION_A_ID),
    )

    repository.create_validated_input_snapshot(
        snapshot,
        project_sample_selection=selection,
    )
    frozen = repository.get_validated_input_binding(snapshot.snapshot_id)
    data_registry.archive_project(
        PROJECT_ID,
        archived_at=NOW + timedelta(seconds=30),
    )
    creation = repository.consume_validated_input_snapshot(
        snapshot.snapshot_id,
        workflow_id="workflow-a",
        expected_build_identity=_identity(),
        record=_record("run-bound"),
        consumed_at=NOW + timedelta(minutes=1),
        event=_event(),
    )

    assert creation.created is True
    assert frozen.project_id == PROJECT_ID
    assert frozen.binding_mode is BindingMode.BOUND_V1
    assert frozen.provenance is BindingProvenance.RESOLVED
    assert frozen.sample_revision_ids == (REVISION_B_ID, REVISION_A_ID)
    assert repository.get_run_data_binding("run-bound") == frozen


def test_selection_is_resolved_before_snapshot_exists_and_fails_atomically() -> None:
    data_registry = _data_registry()
    repository = InMemoryRunRepository(data_registry_repository=data_registry)
    snapshot = _snapshot()
    cross_project = ProjectSampleSelection(
        project_id=OTHER_PROJECT_ID,
        sample_revision_ids=(REVISION_A_ID,),
    )

    with pytest.raises(ProjectSampleSelectionError, match="not eligible"):
        repository.create_validated_input_snapshot(
            snapshot,
            project_sample_selection=cross_project,
        )

    with pytest.raises(KeyError):
        repository.get_validated_input_snapshot(snapshot.snapshot_id)
    with pytest.raises(KeyError):
        repository.get_validated_input_binding(snapshot.snapshot_id)


def test_snapshot_replay_requires_identical_run_binding_evidence() -> None:
    data_registry = _data_registry()
    repository = InMemoryRunRepository(data_registry_repository=data_registry)
    snapshot = _snapshot()
    repository.create_validated_input_snapshot(
        snapshot,
        project_sample_selection=ProjectSampleSelection(
            project_id=PROJECT_ID,
            sample_revision_ids=(REVISION_A_ID,),
        ),
    )
    repository.consume_validated_input_snapshot(
        snapshot.snapshot_id,
        workflow_id="workflow-a",
        expected_build_identity=_identity(),
        record=_record("run-bound"),
        consumed_at=NOW + timedelta(minutes=1),
        event=_event(),
    )
    stored = repository._run_data_bindings["run-bound"]  # type: ignore[attr-defined]
    object.__setattr__(stored, "digest", "f" * 64)

    with pytest.raises(ValidatedSnapshotReplayConflictError, match="binding"):
        repository.consume_validated_input_snapshot(
            snapshot.snapshot_id,
            workflow_id="workflow-a",
            expected_build_identity=_identity(),
            record=_record("run-replay"),
            consumed_at=NOW + timedelta(minutes=1),
            event=_event(),
        )
