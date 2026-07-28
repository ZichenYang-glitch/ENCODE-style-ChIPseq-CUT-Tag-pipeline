"""Snapshot-to-run Project/Sample binding integration contracts."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from threading import Event

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
from encode_pipeline.services.run_repositories import (
    InMemoryRunRepository,
    ProjectSampleSelectionError,
    RunEventDraft,
    ValidatedSnapshotReplayConflictError,
)
from encode_pipeline.services.data_registry_repositories import (
    InMemoryDataRegistryRepository,
)
from encode_pipeline.services.input_file_access import FileObservation
from encode_pipeline.services.input_registry import InputRegistryService
from encode_pipeline.services.input_registry_repositories import (
    InMemoryInputRegistryRepository,
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


def test_every_snapshot_gets_explicit_compatibility_input_envelope() -> None:
    repository = InMemoryRunRepository()
    snapshot = _snapshot()

    repository.create_validated_input_snapshot(snapshot)

    project_samples = repository.get_validated_input_binding(snapshot.snapshot_id)
    inputs = repository.get_validated_input_use_binding(snapshot.snapshot_id)
    assert inputs.project_id == project_samples.project_id
    assert inputs.project_sample_binding_digest == project_samples.digest
    assert inputs.workflow_id == snapshot.workflow_id
    assert inputs.workflow_inputs_digest == snapshot.payload_digest
    assert inputs.contract_mode.value == "compatibility_unresolved_v1"
    assert inputs.input_uses == ()
    assert inputs.fully_managed is False


def test_declared_managed_input_revision_is_frozen_then_copied_after_archive() -> None:
    data_registry = _data_registry()
    input_repository = InMemoryInputRegistryRepository(project_repository=data_registry)
    input_service = InputRegistryService(
        repository=input_repository,
        storage_pool_id_factory=lambda: "stgp_" + "7" * 32,
        input_file_id_factory=lambda: "inpf_" + "8" * 32,
        input_file_revision_id_factory=lambda: "inpfr_" + "9" * 32,
        now_factory=lambda: NOW,
    )
    pool = input_service.create_storage_pool(
        display_name="Approved ingress",
        config_key="ingress-primary",
    )
    input_service.bind_project_storage_pool(PROJECT_ID, pool.storage_pool_id)
    registered = input_service.register_input_file(
        PROJECT_ID,
        stable_key="reads-r1",
        observation=FileObservation(
            relative_path="reads/r1.fastq.gz",
            size_bytes=5,
            content_sha256="a" * 64,
            path_fingerprint=((1, 2), (3, 4, 5, 1, 5, 6, 7)),
        ),
    )
    contract = AdapterInputUseContract(
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
    )
    plan = build_input_use_binding_plan(
        project_id=PROJECT_ID,
        workflow_id="workflow-a",
        contract=contract,
        selections=(
            InputFileRevisionSelection(
                input_use_key="primary-reads",
                occurrence=0,
                input_file_revision_ids=(registered.revision.input_file_revision_id,),
            ),
        ),
    )
    repository = InMemoryRunRepository(
        data_registry_repository=data_registry,
        input_registry_repository=input_repository,
    )
    snapshot = _snapshot()

    repository.create_validated_input_snapshot(
        snapshot,
        project_sample_selection=ProjectSampleSelection(
            project_id=PROJECT_ID,
            sample_revision_ids=(REVISION_A_ID,),
        ),
        input_use_binding_plan=plan,
    )
    frozen = repository.get_validated_input_use_binding(snapshot.snapshot_id)
    input_service.archive_input_file(registered.input_file.input_file_id)
    created = repository.consume_validated_input_snapshot(
        snapshot.snapshot_id,
        workflow_id="workflow-a",
        expected_build_identity=_identity(),
        record=_record("run-managed-input"),
        consumed_at=NOW + timedelta(minutes=1),
        event=_event(),
    )

    assert created.created is True
    assert frozen.fully_managed is True
    assert frozen.input_uses[0].members[0].input_file_revision_id == (
        registered.revision.input_file_revision_id
    )
    assert repository.get_run_input_use_binding("run-managed-input") == frozen


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
    input_binding = repository.get_run_input_use_binding(record.run_id)
    assert input_binding.project_id == LEGACY_PROJECT_ID
    assert input_binding.project_sample_binding_digest == binding.digest
    assert input_binding.contract_mode.value == "compatibility_unresolved_v1"
    assert input_binding.input_uses == ()


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


def test_in_memory_snapshot_selection_is_atomic_with_project_archive() -> None:
    resolved = Event()
    release_resolution = Event()
    archive_started = Event()
    archive_finished = Event()

    class PausingRegistry(InMemoryDataRegistryRepository):
        def resolve_project_sample_selection(
            self,
            selection,
            *,
            workflow_inputs_digest,
        ):
            binding = super().resolve_project_sample_selection(
                selection,
                workflow_inputs_digest=workflow_inputs_digest,
            )
            resolved.set()
            if not release_resolution.wait(timeout=5):
                raise RuntimeError("test did not release selection resolution")
            return binding

    data_registry = PausingRegistry(legacy_created_at=NOW)
    seed = _data_registry()
    for project in seed.list_projects(include_archived=True):
        if project.project_id != LEGACY_PROJECT_ID:
            data_registry.create_project(project)
    for sample_id in (SAMPLE_A_ID, SAMPLE_B_ID):
        sample = seed.get_sample(sample_id)
        data_registry.create_sample(
            sample,
            seed.list_sample_revisions(sample_id)[0],
        )
    repository = InMemoryRunRepository(data_registry_repository=data_registry)
    snapshot = _snapshot()
    selection = ProjectSampleSelection(
        project_id=PROJECT_ID,
        sample_revision_ids=(REVISION_A_ID,),
    )

    def archive_project() -> None:
        archive_started.set()
        data_registry.archive_project(
            PROJECT_ID,
            archived_at=NOW + timedelta(minutes=1),
        )
        archive_finished.set()

    with ThreadPoolExecutor(max_workers=2) as pool:
        snapshot_future = pool.submit(
            repository.create_validated_input_snapshot,
            snapshot,
            project_sample_selection=selection,
        )
        assert resolved.wait(timeout=2)
        archive_future = pool.submit(archive_project)
        assert archive_started.wait(timeout=2)
        archive_was_blocked = not archive_finished.wait(timeout=0.25)
        release_resolution.set()
        snapshot_future.result(timeout=2)
        archive_future.result(timeout=2)

    assert archive_was_blocked
    assert (
        repository.get_validated_input_binding(snapshot.snapshot_id).project_id
        == PROJECT_ID
    )
    assert data_registry.get_project(PROJECT_ID).archived_at is not None


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
