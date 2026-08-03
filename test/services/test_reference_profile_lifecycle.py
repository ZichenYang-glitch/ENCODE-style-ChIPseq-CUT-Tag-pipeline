"""Reference Profile selection, snapshot, and run lifecycle integration tests."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path
from types import SimpleNamespace

import pytest

from encode_pipeline.persistence import open_run_persistence
from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagPreview,
    WorkflowAvailability,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
    WorkspacePlan,
)
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.reference_profiles import (
    AdapterReferenceBindingIdentity,
    BoundWorkflowReference,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Result
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.private_reference_profiles import (
    PrivateReferenceProfileConfig,
)
from encode_pipeline.services.reference_profile_repositories import (
    InMemoryReferenceProfileRepository,
)
from encode_pipeline.services.reference_profiles import ReferenceProfileService
from encode_pipeline.services.reference_profile_runtime import (
    ReferenceProfileBindingService,
    ReferenceProfileRuntimeResolver,
)
from encode_pipeline.services.run_repositories import InMemoryRunRepository
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.run_submission import RunSubmissionService
from encode_pipeline.services.run_submission import RunExecutionUnavailableError
from encode_pipeline.services.validated_inputs import (
    ValidatedInputService,
    ValidatedRunCreationService,
    ValidatedSnapshotExecutionUnavailableError,
)
from encode_pipeline.services.validation import ValidationService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider
from encode_pipeline.workers.jobs import _capture_current_workflow_build


NOW = datetime(2026, 8, 3, 10, 0, tzinfo=timezone.utc)
WORKFLOW_ID = "reference-test-workflow"
PROFILE_ID = "refp_" + "1" * 32
REVISION_IDS = (
    "refpr_" + "2" * 32,
    "refpr_" + "3" * 32,
    "refpr_" + "4" * 32,
)
PRIVATE_REFERENCE_PATH = "/operator/private/references/GRCh38.fa"


class _ReferenceCapableAdapter:
    """Small adapter that makes private resolution observable to the tests."""

    metadata = WorkflowMetadata(
        workflow_id=WORKFLOW_ID,
        name="Reference test workflow",
        version="1.0.0",
    )
    capabilities = WorkflowCapabilities(
        supports=("validation", "workspace_plan", "command")
    )

    def __init__(
        self,
        *,
        build_digest: str = "f" * 64,
        validated_inputs: list[WorkflowInputs] | None = None,
    ) -> None:
        self.build_digest = build_digest
        self.validated_inputs = validated_inputs if validated_inputs is not None else []

    def schema(self) -> WorkflowSchema:
        return WorkflowSchema(schema_version="1.0.0")

    def validate(self, inputs: WorkflowInputs) -> Result[object]:
        self.validated_inputs.append(inputs)
        if "resolved_reference_path" not in inputs.config:
            raise AssertionError("validation did not receive the private binding")
        return Result.success({"accepted": True})

    def preview_dag(self, inputs: WorkflowInputs) -> Result[DagPreview]:
        return Result.success(DagPreview())

    def plan_workspace(
        self,
        inputs: WorkflowInputs,
        workspace: str | Path,
    ) -> Result[WorkspacePlan]:
        return Result.success(WorkspacePlan(directories=[str(workspace)]))

    def build_command(
        self,
        plan: WorkspacePlan,
        workspace: str | Path,
    ) -> Result[CommandSpec]:
        return Result.success(
            CommandSpec(argv=["reference-test"]),
        )

    def extract_artifacts(self, inputs, workspace):
        return Result.success(())

    def execution_availability(self) -> WorkflowAvailability:
        return WorkflowAvailability(
            execution="available",
            reason_code="WORKFLOW_EXECUTION_READY",
        )

    def capture_build_identity(self) -> Result[WorkflowBuildIdentity]:
        return Result.success(
            WorkflowBuildIdentity(
                workflow_id=WORKFLOW_ID,
                adapter_version="1.0.0",
                scheme="sha256-tree-v1",
                logical_entrypoint="tests/reference-workflow",
                digest=self.build_digest,
                captured_at=NOW,
            )
        )

    def verify_reference_profile_binding(self, payload):
        return Result.success(
            AdapterReferenceBindingIdentity(
                workflow_id=WORKFLOW_ID,
                contract_version="reference-test-binding-v1",
                identity_sha256=str(payload["asset_sha256"]),
            )
        )

    def bind_reference_profile(self, inputs: WorkflowInputs, payload):
        identity = self.verify_reference_profile_binding(payload).value
        assert identity is not None
        bound_adapter = type(self)(
            build_digest=identity.identity_sha256,
            validated_inputs=self.validated_inputs,
        )
        bound_inputs = WorkflowInputs(
            config={
                **inputs.config,
                "resolved_reference_path": str(payload["private_path"]),
            },
            samples=inputs.samples,
            options=inputs.options,
        )
        return Result.success(
            BoundWorkflowReference(
                inputs=bound_inputs,
                adapter=bound_adapter,
                identity=identity,
            )
        )


class _LifecycleHarness:
    def __init__(self, tmp_path: Path, backend: str) -> None:
        self.adapter = _ReferenceCapableAdapter()
        self.config_holder = [
            PrivateReferenceProfileConfig(
                {
                    "revision-one": {
                        WORKFLOW_ID: {
                            "private_path": PRIVATE_REFERENCE_PATH,
                            "asset_sha256": "a" * 64,
                        }
                    },
                    "revision-two": {
                        WORKFLOW_ID: {
                            "private_path": PRIVATE_REFERENCE_PATH,
                            "asset_sha256": "b" * 64,
                        }
                    },
                }
            )
        ]
        self.persistence = None
        if backend == "memory":
            self.profile_repository = InMemoryReferenceProfileRepository()
            self.run_repository = InMemoryRunRepository(
                reference_profile_repository=self.profile_repository
            )
        elif backend == "sqlite":
            self.persistence = open_run_persistence(
                f"sqlite:///{tmp_path / 'reference-lifecycle.db'}"
            )
            self.profile_repository = self.persistence.reference_profile_repository
            self.run_repository = self.persistence.repository
        else:  # pragma: no cover - test helper misuse
            raise AssertionError(backend)

        revision_ids = iter(REVISION_IDS)
        self.profiles = ReferenceProfileService(
            repository=self.profile_repository,
            private_config_provider=lambda: self.config_holder[0],
            adapter_provider=lambda workflow_id: (
                self.adapter if workflow_id == WORKFLOW_ID else None
            ),
            profile_id_factory=lambda: PROFILE_ID,
            revision_id_factory=lambda: next(revision_ids),
            now_factory=lambda: NOW,
        )
        self.bindings = ReferenceProfileBindingService(
            repository=self.profile_repository,
            private_config_provider=lambda: self.config_holder[0],
            adapter_provider=lambda workflow_id: (
                self.adapter if workflow_id == WORKFLOW_ID else None
            ),
            now_factory=lambda: NOW,
        )
        self.registry = WorkflowRegistry([self.adapter])
        self.resolver = ReferenceProfileRuntimeResolver(
            self.run_repository,
            self.registry,
            self.bindings,
        )
        self.builds = WorkflowBuildIdentityProvider(
            self.registry,
            project_root=tmp_path.resolve(),
        )
        snapshot_ids = iter(
            (
                "vsnap_" + "5" * 32,
                "vsnap_" + "6" * 32,
                "vsnap_" + "7" * 32,
            )
        )
        self.validation = ValidatedInputService(
            registry=self.registry,
            validation_service=ValidationService(self.registry),
            build_identity_provider=self.builds,
            repository=self.run_repository,
            reference_profile_binding_service=self.bindings,
            reference_profile_catalog=self.profiles,
            snapshot_id_factory=lambda: next(snapshot_ids),
            clock=lambda: NOW,
            snapshot_ttl=timedelta(hours=1),
        )
        run_ids = iter(("run-reference-1", "run-reference-2"))
        self.runs = RunService(
            self.registry,
            id_factory=lambda: next(run_ids),
            repository=self.run_repository,
        )
        self.creation = ValidatedRunCreationService(
            run_service=self.runs,
            build_identity_provider=self.builds,
            reference_profile_binding_service=self.bindings,
            clock=lambda: NOW + timedelta(minutes=1),
        )

    def close(self) -> None:
        if self.persistence is not None:
            self.persistence.close()

    def register(self, config_key: str, display_name: str):
        return self.profiles.register(
            safe_key="human-reference",
            display_name=display_name,
            organism="Homo sapiens",
            assembly="GRCh38",
            config_key=config_key,
        )


@pytest.mark.parametrize("backend", ("memory", "sqlite"))
def test_exact_revision_is_frozen_path_free_and_atomically_copied_to_run(
    tmp_path: Path,
    backend: str,
) -> None:
    harness = _LifecycleHarness(tmp_path, backend)
    try:
        revision = harness.register("revision-one", "GRCh38 primary")
        harness.profiles.enable(PROFILE_ID, revision_id=revision.revision_id)
        submitted = WorkflowInputs(
            config={"user_setting": "visible"},
            samples=[{"sample": "S1"}],
            options={},
        )

        validated = harness.validation.validate(
            WORKFLOW_ID,
            submitted,
            reference_profile_revision_id=revision.revision_id,
        )

        assert validated.is_success
        assert validated.value is not None
        snapshot = validated.value
        snapshot_binding = harness.run_repository.get_validated_reference_binding(
            snapshot.snapshot_id
        )
        assert snapshot_binding is not None
        assert snapshot_binding.revision_id == revision.revision_id
        assert snapshot_binding.revision_public_identity_sha256 == (
            revision.public_identity_sha256
        )
        assert snapshot.workflow_build_identity.digest == "a" * 64
        assert (
            harness.adapter.validated_inputs[-1].config["resolved_reference_path"]
            == PRIVATE_REFERENCE_PATH
        )
        assert PRIVATE_REFERENCE_PATH not in snapshot.canonical_payload
        assert PRIVATE_REFERENCE_PATH not in str(
            snapshot.to_workflow_inputs().to_dict()
        )

        created = harness.creation.create_run(
            WORKFLOW_ID,
            snapshot.snapshot_id,
            tags={"purpose": "reference-test"},
        )

        assert created.created is True
        run_binding = harness.runs.get_run_reference_binding(created.record.run_id)
        assert run_binding == snapshot_binding
        assert PRIVATE_REFERENCE_PATH not in str(created.record.inputs)

        harness.profiles.disable(PROFILE_ID)
        replay = harness.creation.create_run(
            WORKFLOW_ID,
            snapshot.snapshot_id,
            tags={"purpose": "reference-test"},
        )
        assert replay.created is False
        assert replay.record == created.record
        assert harness.runs.get_run_reference_binding(replay.record.run_id) == (
            snapshot_binding
        )
    finally:
        harness.close()


@pytest.mark.parametrize("backend", ("memory", "sqlite"))
def test_start_rechecks_bound_build_identity_not_registered_base(
    tmp_path: Path,
    backend: str,
) -> None:
    harness = _LifecycleHarness(tmp_path, backend)

    class _Queue:
        backend = "rq"
        queue_name = "reference-test"

        def enqueue_execution(self, assignment: RunExecutionAssignment) -> str:
            return assignment.job_id

    try:
        revision = harness.register("revision-one", "GRCh38 primary")
        harness.profiles.enable(PROFILE_ID, revision_id=revision.revision_id)
        snapshot = harness.validation.validate(
            WORKFLOW_ID,
            WorkflowInputs(config={"user_setting": "visible"}),
            reference_profile_revision_id=revision.revision_id,
        ).value
        assert snapshot is not None
        created = harness.creation.create_run(WORKFLOW_ID, snapshot.snapshot_id)
        harness.runs.transition_run(created.record.run_id, RunStatus.VALIDATING)
        harness.runs.complete_preflight(
            created.record.run_id,
            snapshot.workflow_build_identity,
        )
        submission = RunSubmissionService(
            harness.runs,
            _Queue(),
            build_identity_provider=harness.builds,
            reference_profile_resolver=harness.resolver,
        )

        queued = submission.start_run(created.record.run_id)

        assert queued.status is RunStatus.QUEUED
        assert (
            harness.runs.get_workflow_build_identity(created.record.run_id).digest
            == "a" * 64
        )
    finally:
        harness.close()


def test_start_rejects_reference_evidence_without_runtime_resolver(
    tmp_path: Path,
) -> None:
    harness = _LifecycleHarness(tmp_path, "memory")

    class _Queue:
        backend = "rq"
        queue_name = "reference-test"

        def enqueue_execution(self, assignment: RunExecutionAssignment) -> str:
            return assignment.job_id

    try:
        revision = harness.register("revision-one", "GRCh38 primary")
        harness.profiles.enable(PROFILE_ID, revision_id=revision.revision_id)
        snapshot = harness.validation.validate(
            WORKFLOW_ID,
            WorkflowInputs(config={"user_setting": "visible"}),
            reference_profile_revision_id=revision.revision_id,
        ).value
        assert snapshot is not None
        created = harness.creation.create_run(WORKFLOW_ID, snapshot.snapshot_id)
        harness.runs.transition_run(created.record.run_id, RunStatus.VALIDATING)
        harness.runs.complete_preflight(
            created.record.run_id,
            snapshot.workflow_build_identity,
        )
        submission = RunSubmissionService(
            harness.runs,
            _Queue(),
            build_identity_provider=harness.builds,
        )

        with pytest.raises(RunExecutionUnavailableError):
            submission.start_run(created.record.run_id)

        assert harness.runs.get_run(created.record.run_id).status is RunStatus.PLANNED
        assert harness.runs.get_execution_assignment(created.record.run_id) is None
    finally:
        harness.close()


def test_worker_rechecks_bound_build_after_profile_is_disabled(tmp_path: Path) -> None:
    harness = _LifecycleHarness(tmp_path, "memory")
    try:
        revision = harness.register("revision-one", "GRCh38 primary")
        harness.profiles.enable(PROFILE_ID, revision_id=revision.revision_id)
        snapshot = harness.validation.validate(
            WORKFLOW_ID,
            WorkflowInputs(config={"user_setting": "visible"}),
            reference_profile_revision_id=revision.revision_id,
        ).value
        assert snapshot is not None
        record = harness.creation.create_run(WORKFLOW_ID, snapshot.snapshot_id).record
        harness.profiles.disable(PROFILE_ID)

        current = _capture_current_workflow_build(
            SimpleNamespace(
                registry=harness.registry,
                reference_profile_resolver=harness.resolver,
                build_identity_provider=harness.builds,
            ),
            record,
        )

        assert current.is_success
        assert current.value.digest == "a" * 64
    finally:
        harness.close()


def test_missing_disabled_stale_and_identity_drift_fail_closed(tmp_path: Path) -> None:
    harness = _LifecycleHarness(tmp_path, "memory")
    try:
        first = harness.register("revision-one", "GRCh38 r1")
        inputs = WorkflowInputs(config={"user_setting": "visible"})

        missing = harness.validation.validate(WORKFLOW_ID, inputs)
        assert [issue.code for issue in missing.issues] == [
            "REFERENCE_PROFILE_REQUIRED"
        ]

        disabled = harness.validation.validate(
            WORKFLOW_ID,
            inputs,
            reference_profile_revision_id=first.revision_id,
        )
        assert [issue.code for issue in disabled.issues] == [
            "REFERENCE_PROFILE_DISABLED"
        ]

        harness.profiles.enable(PROFILE_ID, revision_id=first.revision_id)
        second = harness.register("revision-two", "GRCh38 r2")
        harness.profiles.enable(PROFILE_ID, revision_id=second.revision_id)
        stale = harness.validation.validate(
            WORKFLOW_ID,
            inputs,
            reference_profile_revision_id=first.revision_id,
        )
        assert [issue.code for issue in stale.issues] == ["REFERENCE_PROFILE_STALE"]

        harness.config_holder[0] = PrivateReferenceProfileConfig(
            {
                "revision-two": {
                    WORKFLOW_ID: {
                        "private_path": PRIVATE_REFERENCE_PATH,
                        "asset_sha256": "c" * 64,
                    }
                }
            }
        )
        drifted = harness.validation.validate(
            WORKFLOW_ID,
            inputs,
            reference_profile_revision_id=second.revision_id,
        )
        assert [issue.code for issue in drifted.issues] == [
            "REFERENCE_PROFILE_IDENTITY_MISMATCH"
        ]
        assert harness.runs.list_runs() == ()
    finally:
        harness.close()


def test_selection_change_changes_frozen_binding_without_polluting_payload(
    tmp_path: Path,
) -> None:
    harness = _LifecycleHarness(tmp_path, "memory")
    try:
        inputs = WorkflowInputs(config={"user_setting": "same"})
        first = harness.register("revision-one", "GRCh38 r1")
        harness.profiles.enable(PROFILE_ID, revision_id=first.revision_id)
        first_result = harness.validation.validate(
            WORKFLOW_ID,
            inputs,
            reference_profile_revision_id=first.revision_id,
        )
        assert first_result.value is not None

        second = harness.register("revision-two", "GRCh38 r2")
        harness.profiles.enable(PROFILE_ID, revision_id=second.revision_id)
        second_result = harness.validation.validate(
            WORKFLOW_ID,
            inputs,
            reference_profile_revision_id=second.revision_id,
        )
        assert second_result.value is not None

        first_snapshot = first_result.value
        second_snapshot = second_result.value
        first_binding = harness.run_repository.get_validated_reference_binding(
            first_snapshot.snapshot_id
        )
        second_binding = harness.run_repository.get_validated_reference_binding(
            second_snapshot.snapshot_id
        )
        assert first_snapshot.snapshot_id != second_snapshot.snapshot_id
        assert first_snapshot.canonical_payload == second_snapshot.canonical_payload
        assert first_snapshot.payload_digest == second_snapshot.payload_digest
        assert first_binding is not None and second_binding is not None
        assert first_binding.revision_id != second_binding.revision_id
        assert first_binding.binding_digest != second_binding.binding_digest
        assert PRIVATE_REFERENCE_PATH not in first_snapshot.canonical_payload
        assert PRIVATE_REFERENCE_PATH not in second_snapshot.canonical_payload
    finally:
        harness.close()


def test_disable_blocks_unconsumed_snapshot_but_preserves_snapshot_evidence(
    tmp_path: Path,
) -> None:
    harness = _LifecycleHarness(tmp_path, "memory")
    try:
        revision = harness.register("revision-one", "GRCh38 r1")
        harness.profiles.enable(PROFILE_ID, revision_id=revision.revision_id)
        result = harness.validation.validate(
            WORKFLOW_ID,
            WorkflowInputs(config={"user_setting": "visible"}),
            reference_profile_revision_id=revision.revision_id,
        )
        assert result.value is not None
        frozen = harness.run_repository.get_validated_reference_binding(
            result.value.snapshot_id
        )

        harness.profiles.disable(PROFILE_ID)

        with pytest.raises(ValidatedSnapshotExecutionUnavailableError):
            harness.creation.create_run(WORKFLOW_ID, result.value.snapshot_id)
        assert harness.runs.list_runs() == ()
        assert (
            harness.run_repository.get_validated_reference_binding(
                result.value.snapshot_id
            )
            == frozen
        )
    finally:
        harness.close()
