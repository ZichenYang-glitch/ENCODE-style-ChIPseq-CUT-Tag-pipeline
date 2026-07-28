"""Fail-closed verification for exact managed input revisions."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from hashlib import sha256
import json
from pathlib import Path
from typing import Callable

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
from encode_pipeline.platform.input_registry import (
    InputFile,
    InputFileRevision,
    InputFileRevisionBindingRef,
    InputProvenanceMode,
    InputUseBindingEnvelope,
    PlannedInputUse,
    ProjectStoragePoolBinding,
    StoragePool,
    build_compatibility_input_binding,
    build_input_file_revision,
    build_input_use_binding,
    build_input_use_binding_envelope,
)
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.platform.snapshots import (
    build_workflow_inputs_digest,
    canonical_workflow_inputs_json,
)
from encode_pipeline.services.managed_input_verification import (
    ManagedInputVerificationService,
)


NOW = datetime(2026, 7, 26, 12, 0, tzinfo=timezone.utc)
PROJECT_ID = "prj_11111111111111111111111111111111"
POOL_ID = "stgp_22222222222222222222222222222222"
INPUT_FILE_ID = "inpf_33333333333333333333333333333333"
REVISION_ID = "inpfr_44444444444444444444444444444444"
WORKFLOW_ID = "managed-test"
CONFIG_KEY = "private-ingress"
RELATIVE_PATH = "reads/example.fastq.gz"
PROJECT_SAMPLE_DIGEST = "5" * 64


class _RunRepository:
    def __init__(
        self,
        run: RunRecord,
        binding: InputUseBindingEnvelope,
    ) -> None:
        self.run = run
        self.binding = binding
        self.run_reads = 0
        self.binding_reads = 0

    def get_run(self, run_id: str) -> RunRecord:
        self.run_reads += 1
        assert run_id == self.run.run_id
        return self.run

    def get_run_input_use_binding(
        self,
        run_id: str,
    ) -> InputUseBindingEnvelope:
        self.binding_reads += 1
        assert run_id == self.run.run_id
        return self.binding


class _InputRegistryRepository:
    def __init__(
        self,
        *,
        pool: StoragePool,
        project_binding: ProjectStoragePoolBinding,
        input_file: InputFile,
        revision: InputFileRevision,
    ) -> None:
        self.pool = pool
        self.project_binding = project_binding
        self.input_file = input_file
        self.revision = revision
        self.calls: list[tuple[str, str]] = []

    def get_project_storage_pool_binding(
        self,
        project_id: str,
    ) -> ProjectStoragePoolBinding:
        self.calls.append(("project_binding", project_id))
        return self.project_binding

    def get_storage_pool(self, storage_pool_id: str) -> StoragePool:
        self.calls.append(("pool", storage_pool_id))
        return self.pool

    def get_input_file(self, input_file_id: str) -> InputFile:
        self.calls.append(("input_file", input_file_id))
        return self.input_file

    def get_input_file_revision(
        self,
        input_file_revision_id: str,
    ) -> InputFileRevision:
        self.calls.append(("revision", input_file_revision_id))
        return self.revision


class _ExplodingInputRegistry:
    def __getattr__(self, name: str) -> object:
        raise AssertionError(f"input registry was read through {name}")


class _WorkflowAdapter:
    metadata = WorkflowMetadata(
        workflow_id=WORKFLOW_ID,
        name="Managed Test",
        version="1.0.0",
    )
    capabilities = WorkflowCapabilities(supports=())

    def schema(self) -> WorkflowSchema:
        return WorkflowSchema(config_schema={"type": "object"})

    def validate(self, inputs: WorkflowInputs) -> Result[object]:
        return Result.success({"validated": True})

    def preview_dag(self, inputs: WorkflowInputs) -> Result[DagPreview]:
        return Result.success(DagPreview())

    def plan_workspace(
        self,
        inputs: WorkflowInputs,
        workspace: str | Path,
    ) -> Result[WorkspacePlan]:
        return Result.success(WorkspacePlan())

    def build_command(
        self,
        plan: WorkspacePlan,
        workspace: str | Path,
    ) -> Result[CommandSpec]:
        return Result.success(CommandSpec(argv=("true",)))

    def extract_artifacts(
        self,
        inputs: WorkflowInputs,
        workspace: str | Path,
    ) -> Result[tuple]:
        return Result.success(())


class _ManagedAdapter(_WorkflowAdapter):
    def __init__(
        self,
        result_factory: Callable[[], object] | None = None,
    ) -> None:
        self._result_factory = result_factory or (
            lambda: Result.success({"scientifically_valid": True})
        )
        self.calls: list[tuple[WorkflowInputs, InputUseBindingEnvelope]] = []

    def validate_managed_input_uses(
        self,
        inputs: WorkflowInputs,
        binding: InputUseBindingEnvelope,
    ) -> object:
        self.calls.append((inputs, binding))
        return self._result_factory()


def _inputs() -> WorkflowInputs:
    return WorkflowInputs(
        config={"assay": "test"},
        samples=[{"sample_id": "s1"}],
        options={"threads": 1},
    )


def _inputs_digest(inputs: WorkflowInputs) -> str:
    return build_workflow_inputs_digest(canonical_workflow_inputs_json(inputs))


def _entities(
    contents: bytes,
    *,
    archived: bool = False,
) -> tuple[
    StoragePool,
    ProjectStoragePoolBinding,
    InputFile,
    InputFileRevision,
]:
    archived_at = NOW + timedelta(minutes=1) if archived else None
    pool = StoragePool(
        storage_pool_id=POOL_ID,
        config_key=CONFIG_KEY,
        display_name="Private ingress",
        created_at=NOW,
        archived_at=archived_at,
    )
    project_binding = ProjectStoragePoolBinding(
        project_id=PROJECT_ID,
        storage_pool_id=POOL_ID,
        bound_at=NOW,
    )
    input_file = InputFile(
        input_file_id=INPUT_FILE_ID,
        project_id=PROJECT_ID,
        storage_pool_id=POOL_ID,
        stable_key="example",
        created_at=NOW,
        archived_at=archived_at,
    )
    revision = build_input_file_revision(
        input_file_revision_id=REVISION_ID,
        input_file_id=INPUT_FILE_ID,
        project_id=PROJECT_ID,
        storage_pool_id=POOL_ID,
        revision_number=1,
        relative_path=RELATIVE_PATH,
        size_bytes=len(contents),
        content_sha256=sha256(contents).hexdigest(),
        created_at=NOW,
    )
    return pool, project_binding, input_file, revision


def _managed_use(revision: InputFileRevision):
    member = InputFileRevisionBindingRef(
        logical_member_key="file",
        input_file_id=revision.input_file_id,
        input_file_revision_id=revision.input_file_revision_id,
        revision_digest=revision.digest,
        size_bytes=revision.size_bytes,
        content_sha256=revision.content_sha256,
    )
    return build_input_use_binding(
        PlannedInputUse(
            key="reads",
            occurrence=0,
            capability_version="reads_v1",
            closure_contract_version="regular_file_v1",
            provenance_mode=InputProvenanceMode.MANAGED_REVISION_V1,
            input_file_revision_ids=(revision.input_file_revision_id,),
        ),
        members=(member,),
    )


def _transitional_use():
    return build_input_use_binding(
        PlannedInputUse(
            key="reference_index",
            occurrence=0,
            capability_version="reference_v1",
            closure_contract_version="directory_index_v1",
            provenance_mode=InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1,
            input_file_revision_ids=(),
        ),
        members=(),
    )


def _declared_binding(
    inputs: WorkflowInputs,
    revision: InputFileRevision,
    *,
    mixed: bool = False,
) -> InputUseBindingEnvelope:
    uses = (_managed_use(revision),)
    if mixed:
        uses += (_transitional_use(),)
    return build_input_use_binding_envelope(
        project_id=PROJECT_ID,
        project_sample_binding_digest=PROJECT_SAMPLE_DIGEST,
        workflow_id=WORKFLOW_ID,
        adapter_contract_version="managed-test-v1",
        workflow_inputs_digest=_inputs_digest(inputs),
        input_uses=uses,
    )


def _all_transitional_binding(
    inputs: WorkflowInputs,
) -> InputUseBindingEnvelope:
    return build_input_use_binding_envelope(
        project_id=PROJECT_ID,
        project_sample_binding_digest=PROJECT_SAMPLE_DIGEST,
        workflow_id=WORKFLOW_ID,
        adapter_contract_version="managed-test-v1",
        workflow_inputs_digest=_inputs_digest(inputs),
        input_uses=(_transitional_use(),),
    )


def _run(inputs: WorkflowInputs) -> RunRecord:
    return RunRecord(
        run_id="run-managed-test",
        workflow_id=WORKFLOW_ID,
        inputs=inputs.to_dict(),
        status=RunStatus.QUEUED,
        created_at=NOW,
        updated_at=NOW,
        started_at=None,
        ended_at=None,
        current_stage=None,
        cancellation_reason=None,
        error=None,
    )


def _write_config(tmp_path: Path, root: Path) -> Path:
    path = tmp_path / "storage-pools.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": "storage-pool-roots-v1",
                "storage_pool_roots": {CONFIG_KEY: str(root)},
            }
        ),
        encoding="utf-8",
    )
    return path


def _managed_fixture(
    tmp_path: Path,
    *,
    mixed: bool = False,
    archived: bool = False,
) -> tuple[
    ManagedInputVerificationService,
    _ManagedAdapter,
    _InputRegistryRepository,
    InputUseBindingEnvelope,
    WorkflowInputs,
    Path,
]:
    root = tmp_path / "pool"
    input_path = root / RELATIVE_PATH
    input_path.parent.mkdir(parents=True)
    contents = b"exact immutable bytes\n"
    input_path.write_bytes(contents)
    pool, project_binding, input_file, revision = _entities(
        contents,
        archived=archived,
    )
    inputs = _inputs()
    binding = _declared_binding(inputs, revision, mixed=mixed)
    run_repository = _RunRepository(_run(inputs), binding)
    input_repository = _InputRegistryRepository(
        pool=pool,
        project_binding=project_binding,
        input_file=input_file,
        revision=revision,
    )
    adapter = _ManagedAdapter()
    service = ManagedInputVerificationService(
        run_repository=run_repository,
        input_registry_repository=input_repository,
        storage_pool_config_path=_write_config(tmp_path, root),
    )
    return service, adapter, input_repository, binding, inputs, input_path


def _assert_redacted_failure(result: Result[None], *secrets: str) -> None:
    assert result.is_failure
    assert result.value is None
    assert len(result.issues) == 1
    issue = result.issues[0]
    assert issue.code == "MANAGED_INPUT_VERIFICATION_FAILED"
    assert issue.source == "input_registry"
    assert issue.path == "input_binding"
    rendered = str(issue.to_dict())
    for secret in secrets:
        assert secret not in rendered


def test_compatibility_run_is_noop_without_catalog_config_or_adapter_lookup(
    tmp_path: Path,
) -> None:
    inputs = _inputs()
    binding = build_compatibility_input_binding(
        project_id=PROJECT_ID,
        project_sample_binding_digest=PROJECT_SAMPLE_DIGEST,
        workflow_id=WORKFLOW_ID,
        adapter_contract_version=None,
        workflow_inputs_digest=_inputs_digest(inputs),
    )
    runs = _RunRepository(_run(inputs), binding)
    service = ManagedInputVerificationService(
        run_repository=runs,
        input_registry_repository=_ExplodingInputRegistry(),
        storage_pool_config_path=tmp_path / "must-not-be-read.json",
    )

    result = service.verify_run("run-managed-test")

    assert result == Result.success(None)
    assert runs.run_reads == 1
    assert runs.binding_reads == 1


def test_all_transitional_binding_is_noop_without_catalog_or_validator() -> None:
    inputs = _inputs()
    service = ManagedInputVerificationService(
        run_repository=object(),
        input_registry_repository=_ExplodingInputRegistry(),
        storage_pool_config_path=None,
    )

    result = service.verify_binding(
        inputs,
        _all_transitional_binding(inputs),
        object(),
    )

    assert result == Result.success(None)


def test_managed_run_fails_closed_before_unqualified_execution_handoff(
    tmp_path: Path,
) -> None:
    service, adapter, repository, binding, inputs, _input_path = _managed_fixture(
        tmp_path,
        mixed=True,
        archived=True,
    )

    result = service.verify_run("run-managed-test")

    assert result.is_failure
    assert result.value is None
    assert [issue.code for issue in result.issues] == [
        "MANAGED_INPUT_EXECUTION_UNAVAILABLE"
    ]
    assert result.issues[0].source == "input_registry"
    assert result.issues[0].path == "input_binding"
    assert repository.calls == []
    assert adapter.calls == []
    assert binding.input_uses[0].provenance_mode is (
        InputProvenanceMode.MANAGED_REVISION_V1
    )
    assert len(binding.input_uses) == 2
    assert inputs == _inputs()


def test_component_verification_rechecks_archived_bytes_and_mixed_envelope(
    tmp_path: Path,
) -> None:
    service, adapter, repository, binding, inputs, _input_path = _managed_fixture(
        tmp_path,
        mixed=True,
        archived=True,
    )

    result = service.verify_binding(inputs, binding, adapter)

    assert result == Result.success(None)
    assert repository.calls == [
        ("project_binding", PROJECT_ID),
        ("revision", REVISION_ID),
        ("input_file", INPUT_FILE_ID),
        ("pool", POOL_ID),
    ]
    assert adapter.calls == [(inputs, binding)]
    assert len(adapter.calls[0][1].input_uses) == 2


@pytest.mark.parametrize("replacement", [b"changed bytes\n", b""])
def test_changed_managed_bytes_fail_closed_and_skip_adapter(
    tmp_path: Path,
    replacement: bytes,
) -> None:
    service, adapter, _repository, binding, inputs, input_path = _managed_fixture(
        tmp_path
    )
    input_path.write_bytes(replacement)

    result = service.verify_binding(inputs, binding, adapter)

    _assert_redacted_failure(result, str(input_path), RELATIVE_PATH, CONFIG_KEY)
    assert adapter.calls == []


def test_symlink_replacement_fails_closed_and_skip_adapter(tmp_path: Path) -> None:
    service, adapter, _repository, binding, inputs, input_path = _managed_fixture(
        tmp_path
    )
    outside = tmp_path / "outside-secret.fastq.gz"
    outside.write_bytes(b"exact immutable bytes\n")
    input_path.unlink()
    input_path.symlink_to(outside)

    result = service.verify_binding(inputs, binding, adapter)

    _assert_redacted_failure(
        result,
        str(outside),
        str(input_path),
        RELATIVE_PATH,
        CONFIG_KEY,
    )
    assert adapter.calls == []


def test_missing_private_config_fails_closed_without_leaking_evidence(
    tmp_path: Path,
) -> None:
    service, adapter, repository, binding, inputs, _input_path = _managed_fixture(
        tmp_path
    )
    service = ManagedInputVerificationService(
        run_repository=object(),
        input_registry_repository=repository,
        storage_pool_config_path=None,
    )

    result = service.verify_binding(inputs, binding, adapter)

    _assert_redacted_failure(result, RELATIVE_PATH, CONFIG_KEY, str(tmp_path))
    assert adapter.calls == []


def test_absent_managed_validator_fails_after_platform_byte_verification(
    tmp_path: Path,
) -> None:
    service, _adapter, repository, binding, inputs, _input_path = _managed_fixture(
        tmp_path
    )
    service = ManagedInputVerificationService(
        run_repository=object(),
        input_registry_repository=repository,
        storage_pool_config_path=tmp_path / "storage-pools.json",
    )

    result = service.verify_binding(inputs, binding, _WorkflowAdapter())

    _assert_redacted_failure(result, RELATIVE_PATH, CONFIG_KEY)


@pytest.mark.parametrize(
    "adapter_result",
    [
        object(),
        Result(value=None),
        Result.failure(
            [
                Issue(
                    code="ADAPTER_SECRET",
                    message="secret adapter path /operator/private/input",
                    technical_message="secret exception text",
                )
            ]
        ),
    ],
)
def test_adapter_failure_or_malformed_result_is_redacted(
    tmp_path: Path,
    adapter_result: object,
) -> None:
    service, _adapter, repository, binding, inputs, _input_path = _managed_fixture(
        tmp_path
    )
    adapter = _ManagedAdapter(lambda: adapter_result)
    service = ManagedInputVerificationService(
        run_repository=object(),
        input_registry_repository=repository,
        storage_pool_config_path=tmp_path / "storage-pools.json",
    )

    result = service.verify_binding(inputs, binding, adapter)

    _assert_redacted_failure(
        result,
        "ADAPTER_SECRET",
        "/operator/private/input",
        "secret exception text",
        RELATIVE_PATH,
        CONFIG_KEY,
    )


def test_adapter_exception_is_redacted(tmp_path: Path) -> None:
    service, _adapter, repository, binding, inputs, _input_path = _managed_fixture(
        tmp_path
    )

    def explode() -> object:
        raise RuntimeError("secret adapter exception /operator/private/input")

    adapter = _ManagedAdapter(explode)
    service = ManagedInputVerificationService(
        run_repository=object(),
        input_registry_repository=repository,
        storage_pool_config_path=tmp_path / "storage-pools.json",
    )

    result = service.verify_binding(inputs, binding, adapter)

    _assert_redacted_failure(
        result,
        "secret adapter exception",
        "/operator/private/input",
        RELATIVE_PATH,
        CONFIG_KEY,
    )


def test_catalog_evidence_mismatch_fails_closed_before_file_or_adapter(
    tmp_path: Path,
) -> None:
    service, adapter, repository, binding, inputs, _input_path = _managed_fixture(
        tmp_path
    )
    repository.revision = build_input_file_revision(
        input_file_revision_id="inpfr_99999999999999999999999999999999",
        input_file_id=repository.revision.input_file_id,
        project_id=repository.revision.project_id,
        storage_pool_id=repository.revision.storage_pool_id,
        revision_number=repository.revision.revision_number,
        relative_path=repository.revision.relative_path,
        size_bytes=repository.revision.size_bytes,
        content_sha256=repository.revision.content_sha256,
        created_at=repository.revision.created_at,
    )

    result = service.verify_binding(inputs, binding, adapter)

    _assert_redacted_failure(result, RELATIVE_PATH, CONFIG_KEY)
    assert adapter.calls == []
