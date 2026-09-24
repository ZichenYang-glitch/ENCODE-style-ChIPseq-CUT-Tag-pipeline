"""H3 adapter consumers with real input/reference bytes and a runtime-loader stub.

These are contract tests: the narrow fixed-tool admission loader is substituted;
no scientific executable, Docker service or Redis instance is required.
"""

from __future__ import annotations

from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import sys

import pytest

from encode_pipeline.adapters.hitrac_preprocess import execution, qualification
from encode_pipeline.adapters.hitrac_preprocess.adapter import HiTracPreprocessAdapter
from encode_pipeline.adapters.hitrac_preprocess.admission import (
    Sample,
    INDEX_SUFFIXES,
    RuntimeBinding,
)
from encode_pipeline.adapters.hitrac_preprocess.deployment import (
    load_default_hitrac_adapter,
)
from encode_pipeline.platform.adapters import (
    CommandSpec,
    QcSummaryExtractingAdapter,
    ReferenceProfileBindingAdapter,
    WorkflowAdapter,
    WorkflowAvailabilityProvidingAdapter,
    WorkflowAvailability,
    WorkflowBuildIdentityProvidingAdapter,
    WorkflowInputs,
    WorkspacePlan,
)
from encode_pipeline.platform.planning import ExecutionPlan, PlanStatus
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.command_builder import CommandBuilder
from encode_pipeline.services.defaults import (
    create_default_process_runner,
    create_default_workflow_registry,
)
from encode_pipeline.services.materialization import WorkspaceMaterializer
from encode_pipeline.services.planning import WorkspacePlanner
from encode_pipeline.services.private_reference_profiles import (
    PrivateReferenceProfileConfig,
)
from encode_pipeline.services.reference_profile_runtime import (
    ReferenceProfileBindingService,
)
from encode_pipeline.services.reference_profiles import ReferenceProfileService
from encode_pipeline.services.run_repositories import InMemoryRunRepository
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.run_submission import (
    RunSubmissionService,
    RunExecutionUnavailableError,
)
from encode_pipeline.services.validated_inputs import (
    ValidatedInputService,
    ValidatedRunCreationService,
    ValidatedSnapshotStaleError,
)
from encode_pipeline.services.validation import ValidationService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider
from encode_pipeline.services.workflow_info import WorkflowInfoService
from encode_pipeline.workers.settings import WorkerSettings

WORKFLOW = "hitrac-preprocess"


def _digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_fastq(path, mate, *, sequence="ACGT"):
    path.write_bytes(
        gzip.compress(f"@read/{mate}\n{sequence}\n+\nIIII\n".encode(), mtime=0)
    )


@pytest.fixture
def bound_case(tmp_path, monkeypatch):
    runtime_path = tmp_path / "operator-runtime.json"
    runtime_path.write_text('{"test_fixture":"fixed-tool-admission-stub"}\n')
    runtime_sha = _digest(runtime_path)
    # Only this expensive tool-tree admission dependency is replaced. Runtime
    # binding file, local interpreter, H3 implementation and reference/input
    # hashing remain the production implementations.
    monkeypatch.setattr(
        execution,
        "load_runtime_binding",
        lambda path: RuntimeBinding(
            prefix=tmp_path,
            python=Path(sys.executable),
            script=tmp_path / "tracPre2.py",
            tools={},
            lock_sha256="c" * 64,
            binding_sha256=runtime_sha,
        ),
    )
    runtime = execution.admit_runtime(runtime_path, runtime_sha)
    fasta = tmp_path / "reference.fa"
    fasta.write_text(">chrA\nACGTACGT\n")
    prefix = tmp_path / "reference"
    files = {"fasta": _digest(fasta)}
    for suffix in INDEX_SUFFIXES:
        index = Path(f"{prefix}{suffix}")
        index.write_text(f"synthetic-not-a-scientific-index:{suffix}\n")
        files[suffix] = _digest(index)
    reference = tmp_path / "reference-binding.json"
    reference.write_text(
        json.dumps(
            {
                "schema_version": "hitrac-reference-binding-v1",
                "fasta": str(fasta),
                "prefix": str(prefix),
                "files": files,
                "contigs": {"chrA": 8},
            }
        )
    )
    payload = {
        "schema_version": "hitrac-reference-profile-v1",
        "binding": str(reference),
        "sha256": _digest(reference),
    }
    rows = []
    for sample in ("zeta", "alpha"):
        paths = [tmp_path / f"{sample}-R{i}.fastq.gz" for i in (1, 2)]
        for mate, path in enumerate(paths, 1):
            _write_fastq(path, mate)
        rows.append(
            {"sample_id": sample, "fastq_1": str(paths[0]), "fastq_2": str(paths[1])}
        )
    inputs = WorkflowInputs(config={}, samples=rows)
    adapter = HiTracPreprocessAdapter(runtime=runtime)
    bound = adapter.bind_reference_profile(inputs, payload)
    assert bound.is_success, bound.issues
    return {
        "adapter": adapter,
        "bound": bound.value.adapter,
        "runtime": runtime,
        "runtime_path": runtime_path,
        "runtime_sha": runtime_sha,
        "payload": payload,
        "inputs": inputs,
        "root": tmp_path,
        "reference": reference,
        "fasta": fasta,
    }


def _registry(adapter):
    return WorkflowRegistry(adapters=(adapter,))


def _catalog(case):
    registry = _registry(case["adapter"])
    private = PrivateReferenceProfileConfig({"tiny": {WORKFLOW: case["payload"]}})
    catalog = ReferenceProfileService(
        private_config_provider=lambda: private,
        adapter_provider=registry.get,
    )
    reference = catalog.register(
        safe_key="tiny",
        display_name="Tiny synthetic contract",
        organism="synthetic",
        assembly="tiny",
        config_key="tiny",
    )
    catalog.enable(reference.profile_id, revision_id=reference.revision_id)
    binder = ReferenceProfileBindingService(
        repository=catalog.repository,
        private_config_provider=lambda: private,
        adapter_provider=registry.get,
    )
    return registry, reference, binder


def test_default_registry_results_stay_unconfigured_without_runtime():
    registry = create_default_workflow_registry(environ={})
    adapter = registry.get(WORKFLOW)
    assert isinstance(adapter, WorkflowAdapter)
    assert isinstance(adapter, ReferenceProfileBindingAdapter)
    assert isinstance(adapter, WorkflowAvailabilityProvidingAdapter)
    assert isinstance(adapter, WorkflowBuildIdentityProvidingAdapter)
    # Authoring-only composition: execution-owned capabilities, and the optional
    # QC protocol they gate, arrive only with an admitted runtime.
    assert not isinstance(adapter, QcSummaryExtractingAdapter)
    assert adapter.capabilities.supports == ("validation", "input_authoring")
    descriptor = WorkflowInfoService(registry).get_descriptor(WORKFLOW).value
    assert descriptor.availability.execution == "not_configured"
    assert descriptor.availability.reason_code == "WORKFLOW_EXECUTION_NOT_CONFIGURED"
    assert descriptor.capabilities.supports == ("validation", "input_authoring")
    assert (
        descriptor.upstream_identity.revision
        == "de6cc732fa00b408551b9f4272933640c08447f1"
    )
    assert adapter.extract_artifacts(
        WorkflowInputs(config={}), Path("/unused")
    ).is_failure


def test_results_runtime_can_open_while_h3_base_and_invalid_runtime_stay_closed(
    bound_case,
):
    case = bound_case
    loaded = load_default_hitrac_adapter(
        {
            "HELIXWEAVE_HITRAC_RUNTIME_BINDING": str(case["runtime_path"]),
            "HELIXWEAVE_HITRAC_RUNTIME_SHA256": case["runtime_sha"],
        }
    )
    assert loaded._runtime is not None
    assert isinstance(loaded, QcSummaryExtractingAdapter)
    assert loaded.requires_atomic_result_publication() is True
    ready = WorkflowInfoService(_registry(loaded)).get_descriptor(WORKFLOW).value
    assert ready.availability.execution == "available"
    assert {"artifact_extract", "qc_summary_extract"}.issubset(
        ready.capabilities.supports
    )
    # Runtime admission cannot manufacture a per-run reference/input binding.
    assert loaded.capture_build_identity().is_failure
    for adapter in (case["adapter"], case["bound"]):
        descriptor = (
            WorkflowInfoService(_registry(adapter)).get_descriptor(WORKFLOW).value
        )
        assert descriptor.availability.execution == "not_configured"
        assert descriptor.capabilities.supports == ("validation", "input_authoring")
        assert (
            WorkflowBuildIdentityProvider(_registry(adapter))
            .capture_executable(WORKFLOW)
            .is_failure
        )
    assert case["bound"].capture_build_identity().is_success
    case["runtime_path"].write_text("changed-runtime")
    assert loaded.execution_availability().execution == "unavailable"
    assert (
        load_default_hitrac_adapter(
            {
                "HELIXWEAVE_HITRAC_RUNTIME_BINDING": str(case["runtime_path"]),
                "HELIXWEAVE_HITRAC_RUNTIME_SHA256": case["runtime_sha"],
            }
        )
        .execution_availability()
        .execution
        == "not_configured"
    )


def test_validation_with_real_reference_selection_cannot_issue_public_snapshot(
    bound_case,
):
    case = bound_case
    registry, reference, binder = _catalog(case)
    repository = InMemoryRunRepository()
    service = ValidatedInputService(
        registry=registry,
        validation_service=ValidationService(registry),
        build_identity_provider=WorkflowBuildIdentityProvider(registry),
        repository=repository,
        reference_profile_binding_service=binder,
    )
    no_reference = service.validate(WORKFLOW, case["inputs"])
    assert no_reference.is_failure
    assert no_reference.issues[0].code == "REFERENCE_PROFILE_REQUIRED"
    checked = service.validate(
        WORKFLOW,
        case["inputs"],
        reference_profile_revision_id=reference.revision_id,
    )
    assert checked.is_success
    assert checked.value is None


def test_private_plan_uses_actual_platform_consumers_and_keeps_attempt_new(bound_case):
    case = bound_case
    adapter = case["bound"]
    registry = _registry(adapter)
    workspace = case["root"] / "workspace"
    before = {
        Path(row[key]): _digest(Path(row[key]))
        for row in case["inputs"].samples
        for key in ("fastq_1", "fastq_2")
    }
    plan = ExecutionPlan(
        plan_id="plan",
        run_id="run",
        workflow_id=WORKFLOW,
        status=PlanStatus.PENDING,
        inputs_snapshot=case["inputs"].to_dict(),
    )
    planned = WorkspacePlanner(registry).plan_workspace(plan, workspace)
    assert planned.is_success, planned.issues
    assert planned.value.workspace_plan.directories == ()
    request = json.loads(dict(planned.value.workspace_plan.files)[execution.REQUEST])
    assert list(request["samples"]) == ["s000001", "s000002"]
    assert [request["samples"][s]["display_id"] for s in request["samples"]] == [
        "zeta",
        "alpha",
    ]
    assert request["inputs"]["options"] == {"threads": 2, "mapq": 10}
    assert request["reference_sha256"] == case["payload"]["sha256"]
    assert request["runtime_sha256"] == case["runtime_sha"]
    materialized = WorkspaceMaterializer().materialize(
        planned.value.workspace_plan, workspace
    )
    assert materialized.is_success
    assert not (workspace / execution.ATTEMPT).exists()
    command = CommandBuilder(registry).build_command(planned.value, workspace)
    assert command.is_success, command.issues
    spec = command.value.command_spec
    assert spec.argv[:4] == (str(Path(sys.executable).absolute()), "-I", "-S", "-B")
    assert spec.argv[4] == str(execution.ROOT / "scripts/run_hitrac_preprocess.py")
    assert spec.argv[5:] == (
        "--request",
        str(workspace / execution.REQUEST),
        "--sha256",
        _digest(workspace / execution.REQUEST),
    )
    assert spec.cwd == str(workspace)
    assert spec.env == {"PYTHONDONTWRITEBYTECODE": "1"}
    assert str(case["root"]) not in json.dumps(spec.to_dict())
    assert spec.managed_container_scope is None
    assert {path: _digest(path) for path in before} == before


@pytest.mark.parametrize(
    "kind", ["input", "reference", "runtime", "implementation", "interpreter"]
)
def test_live_identity_drift_invalidates_frozen_binding_and_plan(
    bound_case, monkeypatch, kind
):
    case = bound_case
    adapter = case["bound"]
    workspace = case["root"] / "workspace"
    plan = adapter.plan_workspace(case["inputs"], workspace)
    original = adapter.capture_build_identity()
    assert original.is_success
    if kind == "input":
        _write_fastq(Path(case["inputs"].samples[0]["fastq_1"]), 1, sequence="TGCA")
    elif kind == "reference":
        case["fasta"].write_text(">chrA\nTGCATGCA\n")
    elif kind == "runtime":
        case["runtime_path"].write_text("changed")
    elif kind == "implementation":
        original_source = execution.source_identity
        monkeypatch.setattr(
            execution,
            "source_identity",
            lambda: {**original_source(), "sha256": "0" * 64},
        )
    else:
        runtime = replace(case["runtime"], python_sha256="0" * 64)
        adapter = HiTracPreprocessAdapter(
            runtime=runtime, binding=replace(adapter._binding, runtime=runtime)
        )
    assert adapter.capture_build_identity().is_failure
    assert adapter.plan_workspace(case["inputs"], workspace).is_failure
    assert adapter.build_command(plan.value, workspace).is_failure


@pytest.mark.parametrize(
    "change", ["options", "order", "extra_file", "request", "workspace", "old_attempt"]
)
def test_command_rejects_wrong_inputs_tampered_plan_or_attempt(bound_case, change):
    case = bound_case
    adapter = case["bound"]
    workspace = case["root"] / "workspace"
    original = adapter.plan_workspace(case["inputs"], workspace)
    assert original.is_success
    plan = original.value
    if change == "options":
        different = WorkflowInputs(
            config={}, samples=case["inputs"].samples, options={"mapq": 11}
        )
        assert adapter.plan_workspace(different, workspace).is_failure
        return
    if change == "order":
        different = WorkflowInputs(
            config={}, samples=list(reversed(case["inputs"].samples))
        )
        assert adapter.plan_workspace(different, workspace).is_failure
        return
    if change == "extra_file":
        plan = WorkspacePlan(files=(*plan.files, ("injected", b"payload")))
    elif change == "request":
        request = json.loads(plan.files[0][1])
        request["inputs"]["options"]["mapq"] = 11
        plan = WorkspacePlan(files=((execution.REQUEST, json.dumps(request).encode()),))
    elif change == "workspace":
        workspace = case["root"] / "another"
    else:
        attempt = workspace / execution.ATTEMPT
        attempt.mkdir(parents=True)
        retained = attempt / "prior-evidence"
        retained.write_bytes(b"retain me")
    assert adapter.build_command(plan, workspace).is_failure
    if change == "old_attempt":
        assert retained.read_bytes() == b"retain me"


@pytest.mark.parametrize("change", ["hash", "schema", "extra", "index"])
def test_bad_reference_selection_never_becomes_execution_authority(bound_case, change):
    case = bound_case
    payload = dict(case["payload"])
    if change == "hash":
        payload["sha256"] = "0" * 64
    elif change == "schema":
        payload["schema_version"] = "unknown"
    elif change == "extra":
        payload["PRIVATE_KEY"] = "/private/path"
    else:
        case["reference"].write_text(
            case["reference"].read_text().replace('".2.bt2":', '".2.bt2l":')
        )
        payload["sha256"] = _digest(case["reference"])
    result = case["adapter"].bind_reference_profile(case["inputs"], payload)
    assert result.is_failure
    assert result.issues[0].code == "HITRAC_REFERENCE_BINDING_INVALID"
    assert str(case["root"]) not in json.dumps(result.to_dict())
    assert "PRIVATE_KEY" not in json.dumps(result.to_dict())


def test_default_runner_authorizes_admitted_interpreter_without_docker(bound_case):
    case = bound_case
    registry = _registry(case["adapter"])
    settings = WorkerSettings(
        database_url=f"sqlite:///{case['root'] / 'unused.db'}",
        redis_url="unix:///not-connected.sock",
        queue_name="test",
        workspace_root=case["root"] / "workspaces",
    )
    runner = create_default_process_runner(registry=registry, settings=settings)
    assert str(Path(sys.executable).absolute()) in runner._allowed_executables
    assert settings.managed_docker_executable is None
    result = runner.run(CommandSpec(argv=(str(case["root"] / "foreign-python"),)))
    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_EXECUTABLE_NOT_ALLOWED"
    case["runtime_path"].write_text("drifted")
    closed = create_default_process_runner(registry=registry, settings=settings)
    assert str(Path(sys.executable).absolute()) not in closed._allowed_executables


class _QualificationOnlyAvailableAdapter(HiTracPreprocessAdapter):
    """Test-only injection to exercise snapshot contracts before H4 publication.

    This class is not production composition and is never in the default registry.
    All binding, validation, planning and identity behavior is the actual adapter.
    """

    def execution_availability(self):
        return WorkflowAvailability()


@pytest.fixture
def qualification_snapshot(bound_case):
    case = bound_case
    private_case = dict(case)
    private_case["adapter"] = _QualificationOnlyAvailableAdapter(
        runtime=case["runtime"]
    )
    registry, reference, binder = _catalog(private_case)
    # Real repositories/services, retaining catalog eligibility at persistence.
    repository = InMemoryRunRepository(reference_profile_repository=binder._repository)
    builds = WorkflowBuildIdentityProvider(registry)
    service = ValidatedInputService(
        registry=registry,
        validation_service=ValidationService(registry),
        build_identity_provider=builds,
        repository=repository,
        reference_profile_binding_service=binder,
    )
    validated = service.validate(
        WORKFLOW,
        case["inputs"],
        reference_profile_revision_id=reference.revision_id,
    )
    assert validated.is_success, validated.issues
    assert validated.value is not None
    creation = ValidatedRunCreationService(
        run_service=RunService(registry=registry, repository=repository),
        build_identity_provider=builds,
        reference_profile_binding_service=binder,
    )
    return case, validated.value, repository, binder, creation


def test_private_snapshot_preserves_original_payload_and_exact_reference_on_replay(
    qualification_snapshot,
):
    case, snapshot, repository, binder, creation = qualification_snapshot
    assert snapshot.to_workflow_inputs().to_dict() == case["inputs"].to_dict()
    assert (
        snapshot.workflow_build_identity.scheme == "sha256-hitrac-execution-binding-v1"
    )
    evidence = repository.get_validated_reference_binding(snapshot.snapshot_id)
    first_bound = binder.resolve_evidence(
        evidence, snapshot.to_workflow_inputs(), require_enabled=True
    )
    second_bound = binder.resolve_evidence(
        evidence, snapshot.to_workflow_inputs(), require_enabled=True
    )
    assert first_bound.is_success and second_bound.is_success
    first_identity = (
        first_bound.value.bound_reference.adapter.capture_build_identity().value
    )
    second_identity = (
        second_bound.value.bound_reference.adapter.capture_build_identity().value
    )
    assert first_identity.matches(second_identity)
    assert snapshot.workflow_build_identity.matches(first_identity)
    first = creation.create_run(WORKFLOW, snapshot.snapshot_id)
    second = creation.create_run(WORKFLOW, snapshot.snapshot_id)
    assert first.created is True
    assert second.created is False
    assert first.record.run_id == second.record.run_id
    assert len(repository.list_runs()) == 1
    assert repository.get_run_reference_binding(first.record.run_id) == evidence
    assert len(repository.list_events(first.record.run_id)) == 1
    assert case["adapter"].execution_availability().execution == "not_configured"


def test_private_snapshot_refuses_new_valid_input_bytes_after_validation(
    qualification_snapshot,
):
    case, snapshot, repository, _, creation = qualification_snapshot
    _write_fastq(Path(case["inputs"].samples[0]["fastq_1"]), 1, sequence="TGCA")
    with pytest.raises(ValidatedSnapshotStaleError):
        creation.create_run(WORKFLOW, snapshot.snapshot_id)
    assert repository.list_runs() == ()
    assert (
        repository.get_validated_input_snapshot(snapshot.snapshot_id).consumed_run_id
        is None
    )


@pytest.mark.parametrize(
    ("kind", "reason"),
    [
        ("inputs", "input_binding_changed"),
        ("implementation", "implementation_binding_changed"),
    ],
)
def test_staging_binding_mismatch_refuses_science_before_first_tool(
    bound_case, monkeypatch, kind, reason
):
    case = bound_case
    runtime = execution.load_runtime_binding(case["runtime_path"])
    monkeypatch.setattr(qualification, "load_runtime_binding", lambda path: runtime)

    def forbidden_science(*args, **kwargs):
        pytest.fail("science must not begin with mismatched staged binding")

    monkeypatch.setattr(qualification, "execute", forbidden_science)
    expected = {
        "runtime_binding_sha256": runtime.binding_sha256,
        "runtime_lock_sha256": runtime.lock_sha256,
        "reference_binding_sha256": case["payload"]["sha256"],
        "samples": json.loads(case["bound"]._binding.samples_json),
    }
    expected_implementation = qualification.implementation_identity()["sha256"]
    if kind == "inputs":
        expected["samples"]["s000001"]["sha256"]["r1"] = "0" * 64
    else:
        expected_implementation = "0" * 64
    attempt = case["root"] / "staging-mismatch"
    result = qualification.qualify(
        runtime_binding=case["runtime_path"],
        reference_binding=case["reference"],
        reference_sha256=case["payload"]["sha256"],
        samples=[
            Sample(row["sample_id"], Path(row["fastq_1"]), Path(row["fastq_2"]))
            for row in case["inputs"].samples
        ],
        attempt=attempt,
        expected_input_identity=expected,
        expected_implementation_sha256=expected_implementation,
    )
    assert result == {"status": "rejected", "reason_code": reason}
    assert (attempt / "input/s000001_R1.fastq.gz").is_file()
    assert not (attempt / "output").exists()
    assert not (attempt / "complete.json").exists()
    assert json.loads((attempt / "private/outcome.json").read_text()) == result


def test_public_start_refuses_even_a_planned_run_with_trusted_build(bound_case):
    case = bound_case
    adapter = case["bound"]
    registry = _registry(adapter)
    repository = InMemoryRunRepository()
    runs = RunService(registry=registry, repository=repository)
    record = runs.create_run(WORKFLOW, case["inputs"])
    runs.transition_run(record.run_id, RunStatus.VALIDATING)
    identity = adapter.capture_build_identity()
    assert identity.is_success
    planned = runs.complete_preflight(record.run_id, identity.value)
    assert planned.status is RunStatus.PLANNED
    before_events = repository.list_events(record.run_id)

    class Queue:
        backend = "rq"
        queue_name = "not-connected"

        def __init__(self):
            self.calls = []

        def enqueue_execution(self, assignment):
            self.calls.append(assignment)
            raise AssertionError("public H3 execution must not reach the queue")

    queue = Queue()
    submission = RunSubmissionService(
        runs,
        queue,
        build_identity_provider=WorkflowBuildIdentityProvider(registry),
    )
    with pytest.raises(RunExecutionUnavailableError):
        submission.start_run(record.run_id)
    assert queue.calls == []
    assert runs.get_run(record.run_id) == planned
    assert repository.list_events(record.run_id) == before_events
    assert repository.get_workflow_build_identity(record.run_id) == identity.value
