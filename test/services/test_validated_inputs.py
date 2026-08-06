from __future__ import annotations

from datetime import datetime, timedelta, timezone

import pytest

from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagPreview,
    WorkflowCapabilities,
    WorkflowAvailability,
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
    WorkspacePlan,
)
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.data_registry import ProjectSampleSelection
from encode_pipeline.platform.input_registry import InputFileRevisionSelection
from encode_pipeline.platform.input_registry import (
    AdapterInputUseContract,
    InputProvenanceMode,
    InputUseBindingPlan,
    InputUseDeclaration,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.services.run_repositories import InMemoryRunRepository
from encode_pipeline.services.run_repositories import InputBindingSelectionError
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.validated_inputs import (
    ValidatedInputService,
    ValidatedRunCreationService,
    ValidatedSnapshotExecutionUnavailableError,
    ValidatedSnapshotNotFoundError,
    ValidatedSnapshotStaleError,
)
from encode_pipeline.services.validation import ValidationService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider


NOW = datetime(2026, 7, 14, 11, 0, tzinfo=timezone.utc)


class FakeAdapter:
    def __init__(
        self,
        result: Result[object] | None = None,
        *,
        trace: list[str] | None = None,
    ) -> None:
        self.metadata = WorkflowMetadata(
            workflow_id="workflow-a",
            name="Workflow A",
            version="1.0.0",
        )
        self.capabilities = WorkflowCapabilities(
            supports=(
                "validation",
                "workspace_plan",
                "command",
            )
        )
        self.result = result or Result.success({"accepted": True})
        self.calls = 0
        self.trace = trace
        self.execution_status = "available"

    def schema(self) -> WorkflowSchema:
        if self.trace is not None:
            self.trace.append("schema")
        return WorkflowSchema(config_schema={"type": "object"})

    def validate(self, inputs: WorkflowInputs) -> Result[object]:
        if self.trace is not None:
            self.trace.append("validate")
        self.calls += 1
        return self.result

    def preview_dag(self, inputs: WorkflowInputs):
        return Result.success(DagPreview())

    def plan_workspace(self, inputs: WorkflowInputs, workspace):
        return Result.success(WorkspacePlan(directories=[str(workspace)]))

    def build_command(self, plan: WorkspacePlan, workspace):
        return Result.success(CommandSpec(argv=["run"]))

    def extract_artifacts(self, inputs, workspace):
        return Result.success(())

    def execution_availability(self) -> WorkflowAvailability:
        reason = {
            "available": "WORKFLOW_EXECUTION_READY",
            "not_configured": "WORKFLOW_EXECUTION_NOT_CONFIGURED",
            "unavailable": "WORKFLOW_EXECUTION_UNAVAILABLE",
        }[self.execution_status]
        return WorkflowAvailability(
            execution=self.execution_status,
            reason_code=reason,
        )

    def capture_build_identity(self):
        return Result.success(_identity())


class FakeBuildProvider:
    def __init__(self, results, *, trace: list[str] | None = None) -> None:
        self.results = list(results)
        self.calls = 0
        self.trace = trace

    def capture_executable(self, workflow_id: str):
        if self.trace is not None:
            self.trace.append("capture")
        self.calls += 1
        result = self.results[min(self.calls - 1, len(self.results) - 1)]
        return result


class DeclaringAdapter(FakeAdapter):
    def __init__(
        self,
        provenance_mode: InputProvenanceMode,
    ) -> None:
        super().__init__()
        self.provenance_mode = provenance_mode
        self.declaration_calls = 0

    def declare_input_uses(
        self,
        inputs: WorkflowInputs,
        validated: object,
    ) -> Result[AdapterInputUseContract]:
        self.declaration_calls += 1
        assert validated == {"accepted": True}
        return Result.success(
            AdapterInputUseContract(
                adapter_contract_version="workflow-a-inputs-v1",
                declarations=(
                    InputUseDeclaration(
                        key="primary_reads",
                        occurrence=0,
                        capability_version="regular-file-v1",
                        closure_contract_version="regular_file_v1",
                        allowed_provenance_modes=(self.provenance_mode,),
                    ),
                ),
                allows_mixed=False,
            )
        )


def _identity(digest: str = "a" * 64) -> WorkflowBuildIdentity:
    return WorkflowBuildIdentity(
        workflow_id="workflow-a",
        adapter_version="1.0.0",
        scheme="sha256-tree-v1",
        logical_entrypoint="workflow/Snakefile",
        digest=digest,
        captured_at=NOW,
    )


def _services(adapter=None, build_results=None, *, trace: list[str] | None = None):
    adapter = adapter or FakeAdapter()
    registry = WorkflowRegistry([adapter])
    repository = InMemoryRunRepository()
    provider = FakeBuildProvider(
        build_results or [Result.success(_identity()), Result.success(_identity())],
        trace=trace,
    )
    validation = ValidatedInputService(
        registry=registry,
        validation_service=ValidationService(registry),
        build_identity_provider=provider,
        repository=repository,
        snapshot_id_factory=lambda: "vsnap_0123456789abcdef0123456789abcdef",
        clock=lambda: NOW,
        snapshot_ttl=timedelta(minutes=30),
    )
    run_ids = iter(("run-1", "run-2", "run-3"))
    run_service = RunService(
        registry,
        id_factory=lambda: next(run_ids),
        repository=repository,
    )
    creation = ValidatedRunCreationService(
        run_service=run_service,
        build_identity_provider=provider,
        clock=lambda: NOW + timedelta(minutes=1),
    )
    return validation, creation, repository, provider


def test_run_creation_rejects_build_provider_from_different_registry() -> None:
    run_registry = WorkflowRegistry([FakeAdapter()])
    identity_registry = WorkflowRegistry([FakeAdapter()])
    run_service = RunService(run_registry)
    provider = WorkflowBuildIdentityProvider(identity_registry)

    with pytest.raises(ValueError, match="registry"):
        ValidatedRunCreationService(
            run_service=run_service,
            build_identity_provider=provider,
        )


def test_successful_validation_persists_snapshot_with_warning_evidence() -> None:
    adapter = FakeAdapter(
        Result.success(
            {"accepted": True},
            issues=[
                Issue(
                    code="INPUT_WARNING",
                    message="Input uses a supported warning case.",
                    severity="warning",
                )
            ],
        )
    )
    service, _, repository, provider = _services(adapter=adapter)
    inputs = WorkflowInputs(
        config={"threads": 1},
        samples=[{"sample": "S1"}],
        options={},
    )

    result = service.validate("workflow-a", inputs)

    assert result.is_success
    assert result.value is not None
    assert result.value.validation_issue_codes == ("INPUT_WARNING",)
    assert result.value.expires_at == NOW + timedelta(minutes=30)
    assert repository.get_validated_input_snapshot(result.value.snapshot_id) == (
        result.value
    )
    assert provider.calls == 2


def test_successful_validation_passes_project_sample_selection_to_repository() -> None:
    service, _, repository, _ = _services()
    selection = ProjectSampleSelection(
        project_id="prj_11111111111111111111111111111111",
        sample_revision_ids=("smpr_22222222222222222222222222222222",),
    )
    observed: list[ProjectSampleSelection | None] = []
    original = repository.create_validated_input_snapshot

    def record_selection(
        snapshot,
        *,
        project_sample_selection=None,
        reference_binding=None,
    ):
        assert reference_binding is None
        observed.append(project_sample_selection)
        return original(snapshot)

    repository.create_validated_input_snapshot = record_selection  # type: ignore[method-assign]

    result = service.validate(
        "workflow-a",
        WorkflowInputs(config={}),
        project_sample_selection=selection,
    )

    assert result.is_success
    assert observed == [selection]


def test_input_revision_selection_requires_exact_adapter_input_use_contract() -> None:
    service, _, repository, _ = _services()

    result = service.validate(
        "workflow-a",
        WorkflowInputs(config={}),
        input_file_revision_selections=(
            InputFileRevisionSelection(
                input_use_key="primary_reads",
                occurrence=0,
                input_file_revision_ids=("inpfr_33333333333333333333333333333333",),
            ),
        ),
    )

    assert result.is_failure
    assert result.issues[0].code == "INPUT_USE_CAPABILITY_UNAVAILABLE"
    with pytest.raises(KeyError):
        repository.get_validated_input_snapshot(
            "vsnap_0123456789abcdef0123456789abcdef"
        )


def test_adapter_declared_transitional_use_is_planned_without_path_inference() -> None:
    class TransitionalAdapter(FakeAdapter):
        def declare_input_uses(
            self,
            inputs: WorkflowInputs,
            validated: object,
        ) -> Result[AdapterInputUseContract]:
            assert validated == {"accepted": True}
            return Result.success(
                AdapterInputUseContract(
                    adapter_contract_version="workflow-a-inputs-v1",
                    declarations=(
                        InputUseDeclaration(
                            key="private-execution-closure",
                            occurrence=0,
                            capability_version="trusted-local-v1",
                            closure_contract_version="private_closure_v1",
                            allowed_provenance_modes=(
                                InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1,
                            ),
                        ),
                    ),
                    allows_mixed=False,
                )
            )

    service, _, repository, _ = _services(adapter=TransitionalAdapter())
    observed: list[InputUseBindingPlan | None] = []
    original = repository.create_validated_input_snapshot

    def record_plan(
        snapshot,
        *,
        project_sample_selection=None,
        input_use_binding_plan=None,
        reference_binding=None,
    ):
        assert reference_binding is None
        observed.append(input_use_binding_plan)
        return original(snapshot)

    repository.create_validated_input_snapshot = record_plan  # type: ignore[method-assign]

    result = service.validate(
        "workflow-a",
        WorkflowInputs(
            config={"private_path": "/operator/adapter-owned"},
        ),
        project_sample_selection=ProjectSampleSelection(
            project_id="prj_11111111111111111111111111111111",
            sample_revision_ids=("smpr_22222222222222222222222222222222",),
        ),
    )

    assert result.is_success
    assert len(observed) == 1
    assert observed[0] is not None
    assert observed[0].input_uses[0].provenance_mode is (
        InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1
    )
    assert observed[0].input_uses[0].input_file_revision_ids == ()


@pytest.mark.parametrize(
    "declaration_behavior",
    ("exception", "malformed", "failure", "wrong_value"),
)
def test_input_use_declaration_boundary_is_stable_and_redacted(
    declaration_behavior: str,
) -> None:
    private_detail = "/operator/private/input-use-contract"

    class BrokenDeclaringAdapter(FakeAdapter):
        def declare_input_uses(
            self,
            inputs: WorkflowInputs,
            validated: object,
        ) -> object:
            if declaration_behavior == "exception":
                raise RuntimeError(private_detail)
            if declaration_behavior == "malformed":
                return object()
            if declaration_behavior == "failure":
                return Result.failure(
                    [
                        Issue(
                            code="PRIVATE_ADAPTER_FAILURE",
                            message=private_detail,
                            technical_message=private_detail,
                        )
                    ]
                )
            return Result.success({"private_path": private_detail})

    service, _, repository, provider = _services(adapter=BrokenDeclaringAdapter())

    result = service.validate(
        "workflow-a",
        WorkflowInputs(config={"legacy_path": "/trusted/local/input"}),
        project_sample_selection=ProjectSampleSelection(
            project_id="prj_11111111111111111111111111111111",
            sample_revision_ids=("smpr_22222222222222222222222222222222",),
        ),
    )

    assert result.is_failure
    assert [issue.code for issue in result.issues] == [
        "INPUT_USE_CAPABILITY_UNAVAILABLE"
    ]
    assert private_detail not in str(result.issues[0].to_dict())
    assert provider.calls == 1
    with pytest.raises(KeyError):
        repository.get_validated_input_snapshot(
            "vsnap_0123456789abcdef0123456789abcdef"
        )


def test_managed_selection_fails_closed_without_scientific_validator() -> None:
    class DeclaringOnlyAdapter(FakeAdapter):
        def declare_input_uses(
            self,
            inputs: WorkflowInputs,
            validated: object,
        ) -> Result[AdapterInputUseContract]:
            return Result.success(
                AdapterInputUseContract(
                    adapter_contract_version="workflow-a-inputs-v1",
                    declarations=(
                        InputUseDeclaration(
                            key="primary_reads",
                            occurrence=0,
                            capability_version="regular-file-v1",
                            closure_contract_version="regular_file_v1",
                            allowed_provenance_modes=(
                                InputProvenanceMode.MANAGED_REVISION_V1,
                            ),
                        ),
                    ),
                    allows_mixed=False,
                )
            )

    service, _, repository, _ = _services(adapter=DeclaringOnlyAdapter())

    result = service.validate(
        "workflow-a",
        WorkflowInputs(config={"legacy_path": "/trusted/local/reads.fastq.gz"}),
        project_sample_selection=ProjectSampleSelection(
            project_id="prj_11111111111111111111111111111111",
            sample_revision_ids=("smpr_22222222222222222222222222222222",),
        ),
        input_file_revision_selections=(
            InputFileRevisionSelection(
                input_use_key="primary_reads",
                occurrence=0,
                input_file_revision_ids=("inpfr_33333333333333333333333333333333",),
            ),
        ),
    )

    assert result.is_failure
    assert result.issues[0].code == "INPUT_USE_CAPABILITY_UNAVAILABLE"
    with pytest.raises(KeyError):
        repository.get_validated_input_snapshot(
            "vsnap_0123456789abcdef0123456789abcdef"
        )


def test_managed_selection_is_not_ignored_when_execution_is_unavailable() -> None:
    class UnavailableDeclaringAdapter(FakeAdapter):
        def __init__(self) -> None:
            super().__init__()
            self.execution_status = "not_configured"

        def declare_input_uses(
            self,
            inputs: WorkflowInputs,
            validated: object,
        ) -> Result[AdapterInputUseContract]:
            raise AssertionError(
                "unqualified managed selection must fail before declaration"
            )

    adapter = UnavailableDeclaringAdapter()
    service, _, repository, provider = _services(adapter=adapter)

    result = service.validate(
        "workflow-a",
        WorkflowInputs(config={"legacy_path": "/trusted/local/reads.fastq.gz"}),
        project_sample_selection=ProjectSampleSelection(
            project_id="prj_11111111111111111111111111111111",
            sample_revision_ids=("smpr_22222222222222222222222222222222",),
        ),
        input_file_revision_selections=(
            InputFileRevisionSelection(
                input_use_key="primary_reads",
                occurrence=0,
                input_file_revision_ids=("inpfr_33333333333333333333333333333333",),
            ),
        ),
    )

    assert result.is_failure
    assert result.issues[0].code == "INPUT_USE_CAPABILITY_UNAVAILABLE"
    assert adapter.calls == 0
    assert provider.calls == 0
    with pytest.raises(KeyError):
        repository.get_validated_input_snapshot(
            "vsnap_0123456789abcdef0123456789abcdef"
        )


def test_managed_selection_requires_exact_project_before_validation() -> None:
    adapter = DeclaringAdapter(InputProvenanceMode.MANAGED_REVISION_V1)
    service, _, repository, provider = _services(adapter=adapter)

    result = service.validate(
        "workflow-a",
        WorkflowInputs(config={}),
        input_file_revision_selections=(
            InputFileRevisionSelection(
                input_use_key="primary_reads",
                occurrence=0,
                input_file_revision_ids=("inpfr_33333333333333333333333333333333",),
            ),
        ),
    )

    assert result.is_failure
    assert result.issues[0].code == "INPUT_BINDING_SELECTION_INVALID"
    assert result.issues[0].path == "project_id"
    assert adapter.calls == 0
    assert adapter.declaration_calls == 0
    assert provider.calls == 0
    with pytest.raises(KeyError):
        repository.get_validated_input_snapshot(
            "vsnap_0123456789abcdef0123456789abcdef"
        )


def test_selection_must_match_the_declared_input_use_contract() -> None:
    adapter = DeclaringAdapter(InputProvenanceMode.MANAGED_REVISION_V1)
    service, _, repository, provider = _services(adapter=adapter)

    result = service.validate(
        "workflow-a",
        WorkflowInputs(config={}),
        project_sample_selection=ProjectSampleSelection(
            project_id="prj_11111111111111111111111111111111",
            sample_revision_ids=("smpr_22222222222222222222222222222222",),
        ),
        input_file_revision_selections=(
            InputFileRevisionSelection(
                input_use_key="undeclared_reads",
                occurrence=0,
                input_file_revision_ids=("inpfr_33333333333333333333333333333333",),
            ),
        ),
    )

    assert result.is_failure
    assert result.issues[0].code == "INPUT_BINDING_SELECTION_INVALID"
    assert result.issues[0].path == "input_selections"
    assert adapter.calls == 1
    assert adapter.declaration_calls == 1
    assert provider.calls == 1
    with pytest.raises(KeyError):
        repository.get_validated_input_snapshot(
            "vsnap_0123456789abcdef0123456789abcdef"
        )


def test_second_build_capture_failure_leaves_no_validated_snapshot() -> None:
    service, _, repository, provider = _services(
        build_results=[
            Result.success(_identity()),
            Result.failure(
                [
                    Issue(
                        code="PRIVATE_BUILD_FAILURE",
                        message="/operator/private/workflow",
                    )
                ]
            ),
        ]
    )

    result = service.validate("workflow-a", WorkflowInputs(config={}))

    assert result.is_failure
    assert result.issues[0].code == "VALIDATION_WORKFLOW_BUILD_UNAVAILABLE"
    assert "/operator/private/workflow" not in str(result.issues[0].to_dict())
    assert provider.calls == 2
    with pytest.raises(KeyError):
        repository.get_validated_input_snapshot(
            "vsnap_0123456789abcdef0123456789abcdef"
        )


def test_repository_input_binding_rejection_is_stable_and_redacted() -> None:
    adapter = DeclaringAdapter(InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1)
    service, _, repository, provider = _services(adapter=adapter)

    def reject_input_binding(
        _snapshot,
        *,
        project_sample_selection=None,
        input_use_binding_plan=None,
        reference_binding=None,
    ):
        assert reference_binding is None
        assert project_sample_selection is not None
        assert input_use_binding_plan is not None
        raise InputBindingSelectionError("private revision /operator/input")

    repository.create_validated_input_snapshot = reject_input_binding  # type: ignore[method-assign]

    result = service.validate(
        "workflow-a",
        WorkflowInputs(config={}),
        project_sample_selection=ProjectSampleSelection(
            project_id="prj_11111111111111111111111111111111",
            sample_revision_ids=("smpr_22222222222222222222222222222222",),
        ),
    )

    assert result.is_failure
    assert result.issues[0].code == "INPUT_BINDING_SELECTION_INVALID"
    assert "/operator/input" not in str(result.issues[0].to_dict())
    assert adapter.declaration_calls == 1
    assert provider.calls == 2


@pytest.mark.parametrize("include_legacy_path", (False, True))
def test_managed_selection_is_unavailable_until_execution_handoff_is_qualified(
    include_legacy_path: bool,
) -> None:
    class ManagedAdapter(FakeAdapter):
        managed_validation_calls = 0

        def declare_input_uses(
            self,
            inputs: WorkflowInputs,
            validated: object,
        ) -> Result[AdapterInputUseContract]:
            return Result.success(
                AdapterInputUseContract(
                    adapter_contract_version="workflow-a-inputs-v1",
                    declarations=(
                        InputUseDeclaration(
                            key="primary_reads",
                            occurrence=0,
                            capability_version="regular-file-v1",
                            closure_contract_version="regular_file_v1",
                            allowed_provenance_modes=(
                                InputProvenanceMode.MANAGED_REVISION_V1,
                            ),
                        ),
                    ),
                    allows_mixed=False,
                )
            )

        def validate_managed_input_uses(
            self,
            inputs: WorkflowInputs,
            binding: object,
        ) -> Result[object]:
            self.managed_validation_calls += 1
            return Result.success({"scientifically_valid": True})

    adapter = ManagedAdapter()
    service, _, repository, _ = _services(adapter=adapter)
    inputs = WorkflowInputs(
        config=(
            {"legacy_path": "/trusted/local/reads.fastq.gz"}
            if include_legacy_path
            else {"assay": "test"}
        )
    )

    result = service.validate(
        "workflow-a",
        inputs,
        project_sample_selection=ProjectSampleSelection(
            project_id="prj_11111111111111111111111111111111",
            sample_revision_ids=("smpr_22222222222222222222222222222222",),
        ),
        input_file_revision_selections=(
            InputFileRevisionSelection(
                input_use_key="primary_reads",
                occurrence=0,
                input_file_revision_ids=("inpfr_33333333333333333333333333333333",),
            ),
        ),
    )

    assert result.is_failure
    assert result.issues[0].code == "INPUT_USE_CAPABILITY_UNAVAILABLE"
    assert adapter.managed_validation_calls == 0
    with pytest.raises(KeyError):
        repository.get_validated_input_snapshot(
            "vsnap_0123456789abcdef0123456789abcdef"
        )


def test_encode_validation_requires_exact_reference_before_snapshot(
    tmp_path,
) -> None:
    adapter = EncodeStyleWorkflowAdapter()
    registry = WorkflowRegistry(
        [adapter],
        legacy_execution_fallbacks=(adapter,),
    )
    repository = InMemoryRunRepository()
    identity = WorkflowBuildIdentity(
        workflow_id=adapter.metadata.workflow_id,
        adapter_version=adapter.metadata.version,
        scheme="sha256-tree-v1",
        logical_entrypoint="workflow/Snakefile",
        digest="c" * 64,
        captured_at=NOW,
    )
    provider = FakeBuildProvider([Result.success(identity), Result.success(identity)])
    service = ValidatedInputService(
        registry=registry,
        validation_service=ValidationService(registry),
        build_identity_provider=provider,
        repository=repository,
        snapshot_id_factory=lambda: "vsnap_abcdef0123456789abcdef0123456789",
        clock=lambda: NOW,
    )
    submitted_config = {
        "replicate_analysis": {"enabled": False},
        "chipseq_idr": {"enabled": False},
    }
    inputs = WorkflowInputs(
        config=submitted_config,
        samples=[
            {
                "sample": "S1",
                "fastq_1": str((tmp_path / "S1.R1.fastq.gz").resolve()),
                "fastq_2": str((tmp_path / "S1.R2.fastq.gz").resolve()),
                "layout": "PE",
                "assay": "chipseq",
                "target": "CTCF",
                "peak_mode": "narrow",
            }
        ],
    )

    result = service.validate(adapter.metadata.workflow_id, inputs)

    assert result.is_failure
    assert [issue.code for issue in result.issues] == ["REFERENCE_PROFILE_REQUIRED"]
    with pytest.raises(KeyError):
        repository.get_validated_input_snapshot(
            "vsnap_abcdef0123456789abcdef0123456789"
        )


def test_schema_contract_is_read_inside_the_stable_build_capture_window() -> None:
    trace: list[str] = []
    adapter = FakeAdapter(trace=trace)
    service, _, _, _ = _services(adapter=adapter, trace=trace)

    result = service.validate("workflow-a", WorkflowInputs(config={}))

    assert result.is_success
    assert trace == ["capture", "schema", "validate", "capture"]


def test_authoring_validation_succeeds_without_execution_and_persists_no_snapshot():
    adapter = FakeAdapter()
    adapter.execution_status = "not_configured"
    service, _, repository, provider = _services(adapter=adapter)

    result = service.validate("workflow-a", WorkflowInputs(config={}))

    assert result.is_success
    assert result.value is None
    assert adapter.calls == 1
    assert provider.calls == 0
    assert repository._validated_input_snapshots == {}


def test_unconsumed_snapshot_cannot_create_run_after_execution_becomes_unavailable():
    adapter = FakeAdapter()
    validation, creation, repository, _ = _services(adapter=adapter)
    validated = validation.validate("workflow-a", WorkflowInputs(config={}))
    assert validated.is_success and validated.value is not None
    adapter.execution_status = "unavailable"

    with pytest.raises(ValidatedSnapshotExecutionUnavailableError):
        creation.create_run("workflow-a", validated.value.snapshot_id)

    assert repository._runs == {}


def test_schema_contract_failure_is_sanitized_before_adapter_validation() -> None:
    private_detail = "private-schema-secret"

    class BrokenSchemaAdapter(FakeAdapter):
        def schema(self) -> WorkflowSchema:
            raise RuntimeError(private_detail)

    adapter = BrokenSchemaAdapter()
    service, _, repository, provider = _services(adapter=adapter)

    result = service.validate("workflow-a", WorkflowInputs(config={}))

    assert result.is_failure
    assert result.errors[0].code == "VALIDATION_WORKFLOW_SCHEMA_UNAVAILABLE"
    assert private_detail not in result.errors[0].message
    assert adapter.calls == 0
    assert provider.calls == 1
    assert repository._validated_input_snapshots == {}


def test_unexpected_snapshot_storage_failure_is_sanitized() -> None:
    private_detail = "private-database-secret"
    service, _, repository, provider = _services()

    def fail_storage(_snapshot) -> None:
        raise Exception(private_detail)

    repository.create_validated_input_snapshot = fail_storage  # type: ignore[method-assign]

    result = service.validate("workflow-a", WorkflowInputs(config={}))

    assert result.is_failure
    assert result.errors[0].code == "VALIDATED_SNAPSHOT_PERSISTENCE_FAILED"
    assert private_detail not in result.errors[0].message
    assert provider.calls == 2


def test_failed_adapter_validation_creates_no_snapshot() -> None:
    adapter = FakeAdapter(
        Result.failure([Issue(code="INPUT_INVALID", message="Input is invalid.")])
    )
    service, _, repository, provider = _services(adapter=adapter)

    result = service.validate("workflow-a", WorkflowInputs(config={}))

    assert result.is_failure
    assert result.errors[0].code == "INPUT_INVALID"
    assert repository._validated_input_snapshots == {}
    assert provider.calls == 1


def test_build_change_during_validation_fails_closed_without_snapshot() -> None:
    service, _, repository, _ = _services(
        build_results=[
            Result.success(_identity("a" * 64)),
            Result.success(_identity("b" * 64)),
        ]
    )

    result = service.validate("workflow-a", WorkflowInputs(config={}))

    assert result.is_failure
    assert result.errors[0].code == "VALIDATION_WORKFLOW_BUILD_CHANGED"
    assert repository._validated_input_snapshots == {}


def test_build_capture_failure_is_sanitized_and_skips_adapter() -> None:
    private_detail = "private-build-source-secret"
    adapter = FakeAdapter()
    service, _, repository, _ = _services(
        adapter=adapter,
        build_results=[
            Result.failure(
                [
                    Issue(
                        code="WORKFLOW_BUILD_SOURCE_UNAVAILABLE",
                        message=private_detail,
                        technical_message=private_detail,
                    )
                ]
            )
        ],
    )

    result = service.validate("workflow-a", WorkflowInputs(config={}))

    assert result.is_failure
    assert result.errors[0].code == "VALIDATION_WORKFLOW_BUILD_UNAVAILABLE"
    assert private_detail not in result.errors[0].message
    assert adapter.calls == 0
    assert repository._validated_input_snapshots == {}


def test_snapshot_only_creation_is_idempotent_and_reuses_canonical_inputs() -> None:
    validation, creation, repository, _ = _services()
    validated = validation.validate(
        "workflow-a",
        WorkflowInputs(
            config={"threads": 1},
            samples=[{"sample": "S1"}],
            options={},
        ),
    ).value
    assert validated is not None

    first = creation.create_run(
        "workflow-a",
        validated.snapshot_id,
        tags={"owner": "lab"},
    )
    replay = creation.create_run(
        "workflow-a",
        validated.snapshot_id,
        tags={"owner": "lab"},
    )

    assert first.created is True
    assert replay.created is False
    assert replay.record == first.record
    assert first.record.inputs == validated.to_workflow_inputs().to_dict()
    assert len(repository.list_runs()) == 1


def test_cross_workflow_snapshot_is_reported_as_not_found() -> None:
    validation, creation, _, _ = _services()
    snapshot = validation.validate("workflow-a", WorkflowInputs(config={})).value
    assert snapshot is not None

    with pytest.raises(ValidatedSnapshotNotFoundError):
        creation.create_run("workflow-b", snapshot.snapshot_id)


def test_current_build_drift_rejects_unconsumed_snapshot_without_run() -> None:
    validation, creation, repository, provider = _services()
    snapshot = validation.validate("workflow-a", WorkflowInputs(config={})).value
    assert snapshot is not None
    provider.results.append(Result.success(_identity("b" * 64)))

    with pytest.raises(ValidatedSnapshotStaleError):
        creation.create_run("workflow-a", snapshot.snapshot_id)

    assert repository.list_runs() == ()
