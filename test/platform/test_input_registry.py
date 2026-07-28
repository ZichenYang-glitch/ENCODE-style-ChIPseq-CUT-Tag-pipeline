"""Tests for workflow-neutral storage and input-binding primitives."""

from __future__ import annotations

from dataclasses import FrozenInstanceError, replace
from datetime import datetime, timedelta, timezone

import pytest

from encode_pipeline.platform.input_registry import (
    INPUT_CLOSURE_DIGEST_SCHEME,
    INPUT_FILE_REVISION_DIGEST_SCHEME,
    INPUT_USE_BINDING_DIGEST_SCHEME,
    AdapterInputUseContract,
    InputBindingContractMode,
    InputFile,
    InputFileRevisionBindingRef,
    InputFileRevisionSelection,
    InputProvenanceMode,
    InputUseBinding,
    InputUseBindingEnvelope,
    InputUseBindingPlan,
    InputUseDeclaration,
    PlannedInputUse,
    ProjectStoragePoolBinding,
    StoragePool,
    build_input_closure_digest,
    build_compatibility_input_binding,
    build_input_file_revision,
    build_input_use_binding,
    build_input_use_binding_envelope,
    build_input_use_binding_plan,
    validate_input_file_relative_path,
)


NOW = datetime(2026, 7, 26, 8, 0, tzinfo=timezone.utc)
PROJECT_ID = "prj_" + "1" * 32
POOL_ID = "stgp_" + "2" * 32
INPUT_FILE_ID = "inpf_" + "3" * 32
REVISION_ID = "inpfr_" + "4" * 32
PROJECT_SAMPLE_DIGEST = "5" * 64
WORKFLOW_INPUTS_DIGEST = "6" * 64
CONTENT_SHA256 = "7" * 64
MAX_COLLECTION_ITEMS = 256


def _revision_ids(count: int) -> tuple[str, ...]:
    return tuple(f"inpfr_{index:032x}" for index in range(count))


def _declarations(count: int) -> tuple[InputUseDeclaration, ...]:
    return tuple(
        InputUseDeclaration(
            key="reads",
            occurrence=index,
            capability_version="regular-file-v1",
            closure_contract_version="regular_file_v1",
            allowed_provenance_modes=(InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1,),
        )
        for index in range(count)
    )


def _planned_uses(count: int) -> tuple[PlannedInputUse, ...]:
    return tuple(
        PlannedInputUse(
            key="reads",
            occurrence=index,
            capability_version="regular-file-v1",
            closure_contract_version="regular_file_v1",
            provenance_mode=InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1,
            input_file_revision_ids=(),
        )
        for index in range(count)
    )


def _selections(count: int) -> tuple[InputFileRevisionSelection, ...]:
    return tuple(
        InputFileRevisionSelection(
            input_use_key="reads",
            occurrence=index,
            input_file_revision_ids=(REVISION_ID,),
        )
        for index in range(count)
    )


def _members(count: int) -> tuple[InputFileRevisionBindingRef, ...]:
    return tuple(
        InputFileRevisionBindingRef(
            logical_member_key=f"member-{index}",
            input_file_id=f"inpf_{index:032x}",
            input_file_revision_id=f"inpfr_{index:032x}",
            revision_digest="8" * 64,
            size_bytes=index,
            content_sha256=CONTENT_SHA256,
        )
        for index in range(count)
    )


def _transitional_bindings(count: int) -> tuple[InputUseBinding, ...]:
    return tuple(
        InputUseBinding(
            key="reads",
            occurrence=index,
            capability_version="regular-file-v1",
            closure_contract_version="regular_file_v1",
            provenance_mode=InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1,
            members=(),
            closure_digest_scheme=None,
            closure_digest=None,
        )
        for index in range(count)
    )


def test_storage_records_use_opaque_ids_and_normalize_utc_lifecycle() -> None:
    offset = timezone(timedelta(hours=8))
    pool = StoragePool(
        storage_pool_id=POOL_ID,
        config_key="primary-inputs",
        display_name=" Primary inputs ",
        created_at=datetime(2026, 7, 26, 16, 0, tzinfo=offset),
        archived_at=datetime(2026, 7, 27, 16, 0, tzinfo=offset),
    )
    binding = ProjectStoragePoolBinding(
        project_id=PROJECT_ID,
        storage_pool_id=POOL_ID,
        bound_at=datetime(2026, 7, 26, 16, 0, tzinfo=offset),
    )
    input_file = InputFile(
        input_file_id=INPUT_FILE_ID,
        project_id=PROJECT_ID,
        storage_pool_id=POOL_ID,
        stable_key="reads-R1",
        created_at=datetime(2026, 7, 26, 16, 0, tzinfo=offset),
    )

    assert pool.display_name == "Primary inputs"
    assert pool.created_at == binding.bound_at == input_file.created_at == NOW
    assert pool.created_at.tzinfo is timezone.utc
    with pytest.raises(FrozenInstanceError):
        input_file.stable_key = "replacement"  # type: ignore[misc]


@pytest.mark.parametrize(
    "relative_path",
    (
        "",
        "/absolute.fastq.gz",
        r"windows\path.fastq.gz",
        ".",
        "..",
        "reads/../escape.fastq.gz",
        "reads/./alias.fastq.gz",
        "reads//empty.fastq.gz",
        "reads/control\x00.fastq.gz",
    ),
)
def test_input_file_relative_path_is_nested_posix_and_traversal_safe(
    relative_path: str,
) -> None:
    with pytest.raises(ValueError, match="relative POSIX"):
        validate_input_file_relative_path(relative_path)


def test_input_file_revision_digest_is_immutable_and_path_sensitive() -> None:
    revision = build_input_file_revision(
        input_file_revision_id=REVISION_ID,
        input_file_id=INPUT_FILE_ID,
        project_id=PROJECT_ID,
        storage_pool_id=POOL_ID,
        revision_number=1,
        relative_path="cohort/reads-R1.fastq.gz",
        size_bytes=123,
        content_sha256=CONTENT_SHA256,
        created_at=NOW,
    )
    moved = build_input_file_revision(
        input_file_revision_id="inpfr_" + "8" * 32,
        input_file_id=INPUT_FILE_ID,
        project_id=PROJECT_ID,
        storage_pool_id=POOL_ID,
        revision_number=2,
        relative_path="archive/reads-R1.fastq.gz",
        size_bytes=123,
        content_sha256=CONTENT_SHA256,
        created_at=NOW,
    )

    assert revision.digest_scheme == INPUT_FILE_REVISION_DIGEST_SCHEME
    assert len(revision.digest) == 64
    assert revision.digest != moved.digest
    with pytest.raises(ValueError, match="digest"):
        replace(revision, digest="0" * 64)


def _contract(*, allows_mixed: bool = True) -> AdapterInputUseContract:
    return AdapterInputUseContract(
        adapter_contract_version="encode-v1",
        declarations=(
            InputUseDeclaration(
                key="reads",
                occurrence=0,
                capability_version="regular-file-v1",
                closure_contract_version="regular_file_v1",
                allowed_provenance_modes=(InputProvenanceMode.MANAGED_REVISION_V1,),
            ),
            InputUseDeclaration(
                key="genome-index",
                occurrence=0,
                capability_version="bowtie2-prefix-v1",
                closure_contract_version="bowtie2_prefix_v1",
                allowed_provenance_modes=(
                    InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1,
                ),
            ),
        ),
        allows_mixed=allows_mixed,
    )


def test_contract_builds_ordered_explicit_mixed_plan_without_path_inference() -> None:
    plan = build_input_use_binding_plan(
        project_id=PROJECT_ID,
        workflow_id="encode-chipseq",
        contract=_contract(),
        selections=(
            InputFileRevisionSelection(
                input_use_key="reads",
                occurrence=0,
                input_file_revision_ids=(REVISION_ID,),
            ),
        ),
    )

    assert [item.key for item in plan.input_uses] == ["reads", "genome-index"]
    assert plan.input_uses[0].provenance_mode is (
        InputProvenanceMode.MANAGED_REVISION_V1
    )
    assert plan.input_uses[0].input_file_revision_ids == (REVISION_ID,)
    assert plan.input_uses[1].provenance_mode is (
        InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1
    )
    assert plan.input_uses[1].input_file_revision_ids == ()
    assert plan.fully_managed is False


def test_workflow_id_matches_existing_workflow_metadata_contract() -> None:
    workflow_id = "workflow+" + "v" * 247
    plan = build_input_use_binding_plan(
        project_id=PROJECT_ID,
        workflow_id=f" {workflow_id} ",
        contract=_contract(),
        selections=(
            InputFileRevisionSelection(
                input_use_key="reads",
                occurrence=0,
                input_file_revision_ids=(REVISION_ID,),
            ),
        ),
    )
    compatibility = build_compatibility_input_binding(
        project_id=PROJECT_ID,
        project_sample_binding_digest=PROJECT_SAMPLE_DIGEST,
        workflow_id=f" {workflow_id} ",
        adapter_contract_version=None,
        workflow_inputs_digest=WORKFLOW_INPUTS_DIGEST,
    )

    assert len(workflow_id) == 256
    assert plan.workflow_id == workflow_id
    assert compatibility.workflow_id == workflow_id


def test_contract_rejects_duplicate_uses_extra_selections_and_forbidden_mix() -> None:
    declaration = _contract().declarations[0]
    with pytest.raises(ValueError, match="duplicate"):
        AdapterInputUseContract(
            adapter_contract_version="v1",
            declarations=(declaration, declaration),
            allows_mixed=False,
        )
    with pytest.raises(ValueError, match="undeclared"):
        build_input_use_binding_plan(
            project_id=PROJECT_ID,
            workflow_id="encode-chipseq",
            contract=_contract(),
            selections=(
                InputFileRevisionSelection(
                    input_use_key="other",
                    occurrence=0,
                    input_file_revision_ids=(REVISION_ID,),
                ),
            ),
        )
    with pytest.raises(ValueError, match="mixed"):
        build_input_use_binding_plan(
            project_id=PROJECT_ID,
            workflow_id="encode-chipseq",
            contract=_contract(allows_mixed=False),
            selections=(
                InputFileRevisionSelection(
                    input_use_key="reads",
                    occurrence=0,
                    input_file_revision_ids=(REVISION_ID,),
                ),
            ),
        )


def test_regular_file_managed_use_requires_exactly_one_opaque_revision() -> None:
    with pytest.raises(ValueError, match="exactly one"):
        build_input_use_binding_plan(
            project_id=PROJECT_ID,
            workflow_id="encode-chipseq",
            contract=_contract(),
            selections=(
                InputFileRevisionSelection(
                    input_use_key="reads",
                    occurrence=0,
                    input_file_revision_ids=(
                        REVISION_ID,
                        "inpfr_" + "9" * 32,
                    ),
                ),
            ),
        )


def test_domain_collections_accept_exactly_256_items() -> None:
    revision_ids = _revision_ids(MAX_COLLECTION_ITEMS)
    declarations = _declarations(MAX_COLLECTION_ITEMS)
    planned_uses = _planned_uses(MAX_COLLECTION_ITEMS)
    members = _members(MAX_COLLECTION_ITEMS)
    bindings = _transitional_bindings(MAX_COLLECTION_ITEMS)

    selection = InputFileRevisionSelection(
        input_use_key="bundle",
        occurrence=0,
        input_file_revision_ids=revision_ids,
    )
    contract = AdapterInputUseContract(
        adapter_contract_version="adapter-v1",
        declarations=declarations,
        allows_mixed=True,
    )
    planned = PlannedInputUse(
        key="bundle",
        occurrence=0,
        capability_version="bundle-v1",
        closure_contract_version="bundle_v1",
        provenance_mode=InputProvenanceMode.MANAGED_REVISION_V1,
        input_file_revision_ids=revision_ids,
    )
    plan = InputUseBindingPlan(
        project_id=PROJECT_ID,
        workflow_id="test-workflow",
        adapter_contract_version="adapter-v1",
        input_uses=planned_uses,
    )
    closure_digest = build_input_closure_digest(
        closure_contract_version="bundle_v1",
        members=members,
    )
    binding = InputUseBinding(
        key="bundle",
        occurrence=0,
        capability_version="bundle-v1",
        closure_contract_version="bundle_v1",
        provenance_mode=InputProvenanceMode.MANAGED_REVISION_V1,
        members=members,
        closure_digest_scheme=INPUT_CLOSURE_DIGEST_SCHEME,
        closure_digest=closure_digest,
    )
    envelope = build_input_use_binding_envelope(
        project_id=PROJECT_ID,
        project_sample_binding_digest=PROJECT_SAMPLE_DIGEST,
        workflow_id="test-workflow",
        adapter_contract_version="adapter-v1",
        workflow_inputs_digest=WORKFLOW_INPUTS_DIGEST,
        input_uses=bindings,
    )

    assert len(selection.input_file_revision_ids) == MAX_COLLECTION_ITEMS
    assert len(contract.declarations) == MAX_COLLECTION_ITEMS
    assert len(planned.input_file_revision_ids) == MAX_COLLECTION_ITEMS
    assert len(plan.input_uses) == MAX_COLLECTION_ITEMS
    assert len(binding.members) == MAX_COLLECTION_ITEMS
    assert len(envelope.input_uses) == MAX_COLLECTION_ITEMS


@pytest.mark.parametrize(
    "factory",
    (
        pytest.param(
            lambda: InputFileRevisionSelection(
                input_use_key="bundle",
                occurrence=0,
                input_file_revision_ids=_revision_ids(MAX_COLLECTION_ITEMS + 1),
            ),
            id="selection-revisions",
        ),
        pytest.param(
            lambda: AdapterInputUseContract(
                adapter_contract_version="adapter-v1",
                declarations=_declarations(MAX_COLLECTION_ITEMS + 1),
                allows_mixed=True,
            ),
            id="adapter-declarations",
        ),
        pytest.param(
            lambda: PlannedInputUse(
                key="bundle",
                occurrence=0,
                capability_version="bundle-v1",
                closure_contract_version="bundle_v1",
                provenance_mode=InputProvenanceMode.MANAGED_REVISION_V1,
                input_file_revision_ids=_revision_ids(MAX_COLLECTION_ITEMS + 1),
            ),
            id="planned-revisions",
        ),
        pytest.param(
            lambda: InputUseBindingPlan(
                project_id=PROJECT_ID,
                workflow_id="test-workflow",
                adapter_contract_version="adapter-v1",
                input_uses=_planned_uses(MAX_COLLECTION_ITEMS + 1),
            ),
            id="plan-input-uses",
        ),
        pytest.param(
            lambda: build_input_use_binding_plan(
                project_id=PROJECT_ID,
                workflow_id="test-workflow",
                contract=AdapterInputUseContract(
                    adapter_contract_version="adapter-v1",
                    declarations=(),
                    allows_mixed=True,
                ),
                selections=_selections(MAX_COLLECTION_ITEMS + 1),
            ),
            id="plan-selections",
        ),
        pytest.param(
            lambda: build_input_closure_digest(
                closure_contract_version="bundle_v1",
                members=_members(MAX_COLLECTION_ITEMS + 1),
            ),
            id="closure-digest-members",
        ),
        pytest.param(
            lambda: InputUseBinding(
                key="bundle",
                occurrence=0,
                capability_version="bundle-v1",
                closure_contract_version="bundle_v1",
                provenance_mode=InputProvenanceMode.MANAGED_REVISION_V1,
                members=_members(MAX_COLLECTION_ITEMS + 1),
                closure_digest_scheme=INPUT_CLOSURE_DIGEST_SCHEME,
                closure_digest="0" * 64,
            ),
            id="binding-members",
        ),
        pytest.param(
            lambda: InputUseBindingEnvelope(
                project_id=PROJECT_ID,
                project_sample_binding_digest=PROJECT_SAMPLE_DIGEST,
                workflow_id="test-workflow",
                adapter_contract_version="adapter-v1",
                workflow_inputs_digest=WORKFLOW_INPUTS_DIGEST,
                contract_mode=InputBindingContractMode.DECLARED_INPUT_USES_V1,
                input_uses=_transitional_bindings(MAX_COLLECTION_ITEMS + 1),
                digest_scheme=INPUT_USE_BINDING_DIGEST_SCHEME,
                digest="0" * 64,
            ),
            id="envelope-input-uses",
        ),
        pytest.param(
            lambda: build_input_use_binding_envelope(
                project_id=PROJECT_ID,
                project_sample_binding_digest=PROJECT_SAMPLE_DIGEST,
                workflow_id="test-workflow",
                adapter_contract_version="adapter-v1",
                workflow_inputs_digest=WORKFLOW_INPUTS_DIGEST,
                input_uses=_transitional_bindings(MAX_COLLECTION_ITEMS + 1),
            ),
            id="envelope-builder-input-uses",
        ),
    ),
)
def test_domain_collections_reject_257_items(factory) -> None:
    with pytest.raises(ValueError, match="at most 256"):
        factory()


def test_declared_envelope_digest_binds_project_sample_and_ordered_input_uses() -> None:
    plan = build_input_use_binding_plan(
        project_id=PROJECT_ID,
        workflow_id="encode-chipseq",
        contract=_contract(),
        selections=(
            InputFileRevisionSelection(
                input_use_key="reads",
                occurrence=0,
                input_file_revision_ids=(REVISION_ID,),
            ),
        ),
    )
    member = InputFileRevisionBindingRef(
        logical_member_key="file",
        input_file_id=INPUT_FILE_ID,
        input_file_revision_id=REVISION_ID,
        revision_digest="8" * 64,
        size_bytes=123,
        content_sha256=CONTENT_SHA256,
    )
    resolved = (
        build_input_use_binding(plan.input_uses[0], members=(member,)),
        build_input_use_binding(plan.input_uses[1], members=()),
    )

    envelope = build_input_use_binding_envelope(
        project_id=PROJECT_ID,
        project_sample_binding_digest=PROJECT_SAMPLE_DIGEST,
        workflow_id=plan.workflow_id,
        adapter_contract_version=plan.adapter_contract_version,
        workflow_inputs_digest=WORKFLOW_INPUTS_DIGEST,
        input_uses=resolved,
    )
    changed_sample_binding = build_input_use_binding_envelope(
        project_id=PROJECT_ID,
        project_sample_binding_digest="9" * 64,
        workflow_id=plan.workflow_id,
        adapter_contract_version=plan.adapter_contract_version,
        workflow_inputs_digest=WORKFLOW_INPUTS_DIGEST,
        input_uses=resolved,
    )

    assert envelope.contract_mode is InputBindingContractMode.DECLARED_INPUT_USES_V1
    assert envelope.digest_scheme == INPUT_USE_BINDING_DIGEST_SCHEME
    assert envelope.fully_managed is False
    assert envelope.digest != changed_sample_binding.digest
    with pytest.raises(ValueError, match="digest"):
        replace(envelope, digest="0" * 64)


def test_compatibility_envelope_is_explicitly_unresolved_and_has_no_uses() -> None:
    envelope = build_compatibility_input_binding(
        project_id=PROJECT_ID,
        project_sample_binding_digest=PROJECT_SAMPLE_DIGEST,
        workflow_id="encode-chipseq",
        adapter_contract_version=None,
        workflow_inputs_digest=WORKFLOW_INPUTS_DIGEST,
    )

    assert envelope.contract_mode is (
        InputBindingContractMode.COMPATIBILITY_UNRESOLVED_V1
    )
    assert envelope.input_uses == ()
    assert envelope.fully_managed is False
    assert len(envelope.digest) == 64


def test_transitional_binding_cannot_smuggle_members_or_closure_evidence() -> None:
    with pytest.raises(ValueError, match="transitional"):
        InputUseBinding(
            key="genome-index",
            occurrence=0,
            capability_version="bowtie2-prefix-v1",
            closure_contract_version="bowtie2_prefix_v1",
            provenance_mode=InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1,
            members=(
                InputFileRevisionBindingRef(
                    logical_member_key="index",
                    input_file_id=INPUT_FILE_ID,
                    input_file_revision_id=REVISION_ID,
                    revision_digest="8" * 64,
                    size_bytes=123,
                    content_sha256=CONTENT_SHA256,
                ),
            ),
            closure_digest_scheme=None,
            closure_digest=None,
        )
