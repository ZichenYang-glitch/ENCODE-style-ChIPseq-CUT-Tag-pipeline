from __future__ import annotations

from datetime import datetime, timedelta, timezone

import pytest

from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagPreview,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
    WorkspacePlan,
)
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.services.run_repositories import InMemoryRunRepository
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.validated_inputs import (
    ValidatedInputService,
    ValidatedRunCreationService,
    ValidatedSnapshotNotFoundError,
    ValidatedSnapshotStaleError,
)
from encode_pipeline.services.validation import ValidationService


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
        self.capabilities = WorkflowCapabilities(supports=("validation",))
        self.result = result or Result.success({"accepted": True})
        self.calls = 0
        self.trace = trace

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

    def build_command(self, plan: WorkspacePlan):
        return Result.success(CommandSpec(argv=["run"]))

    def extract_artifacts(self, inputs, workspace):
        return Result.success(())


class FakeBuildProvider:
    def __init__(self, results, *, trace: list[str] | None = None) -> None:
        self.results = list(results)
        self.calls = 0
        self.trace = trace

    def capture(self, workflow_id: str):
        if self.trace is not None:
            self.trace.append("capture")
        self.calls += 1
        result = self.results[min(self.calls - 1, len(self.results) - 1)]
        return result


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


def test_schema_contract_is_read_inside_the_stable_build_capture_window() -> None:
    trace: list[str] = []
    adapter = FakeAdapter(trace=trace)
    service, _, _, _ = _services(adapter=adapter, trace=trace)

    result = service.validate("workflow-a", WorkflowInputs(config={}))

    assert result.is_success
    assert trace == ["capture", "schema", "validate", "capture"]


def test_schema_contract_failure_is_sanitized_before_adapter_validation() -> None:
    class BrokenSchemaAdapter(FakeAdapter):
        def schema(self) -> WorkflowSchema:
            raise RuntimeError("private schema failure /home/user/workflow")

    adapter = BrokenSchemaAdapter()
    service, _, repository, provider = _services(adapter=adapter)

    result = service.validate("workflow-a", WorkflowInputs(config={}))

    assert result.is_failure
    assert result.errors[0].code == "VALIDATION_WORKFLOW_SCHEMA_UNAVAILABLE"
    assert "/home/" not in result.errors[0].message
    assert adapter.calls == 0
    assert provider.calls == 1
    assert repository._validated_input_snapshots == {}


def test_unexpected_snapshot_storage_failure_is_sanitized() -> None:
    service, _, repository, provider = _services()

    def fail_storage(_snapshot) -> None:
        raise Exception("private database failure /home/user/platform.sqlite3")

    repository.create_validated_input_snapshot = fail_storage  # type: ignore[method-assign]

    result = service.validate("workflow-a", WorkflowInputs(config={}))

    assert result.is_failure
    assert result.errors[0].code == "VALIDATED_SNAPSHOT_PERSISTENCE_FAILED"
    assert "/home/" not in result.errors[0].message
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
    adapter = FakeAdapter()
    service, _, repository, _ = _services(
        adapter=adapter,
        build_results=[
            Result.failure(
                [
                    Issue(
                        code="WORKFLOW_BUILD_SOURCE_UNAVAILABLE",
                        message="private /home/user/source",
                        technical_message="private /home/user/source",
                    )
                ]
            )
        ],
    )

    result = service.validate("workflow-a", WorkflowInputs(config={}))

    assert result.is_failure
    assert result.errors[0].code == "VALIDATION_WORKFLOW_BUILD_UNAVAILABLE"
    assert "/home/" not in result.errors[0].message
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
