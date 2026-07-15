"""Read-only Omics Intake Bundle contract consumption tests."""

from __future__ import annotations

from copy import deepcopy
import os
import sys

import pytest
import yaml

import encode_pipeline.services.input_bundle_imports as input_bundle_module
from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagPreview,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
    WorkspacePlan,
)
from encode_pipeline.platform.input_bundles import InputBundleMapping
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Result
from encode_pipeline.services.input_bundle_imports import InputBundleImportService
from input_bundle_support import (
    bowtie2_shard_paths,
    configure_bowtie2_shards,
    create_runnable_bundle,
    read_bundle,
    refresh_artifact,
    tree_digest,
    write_bundle,
)


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


class ContractOnlyAdapter:
    """Minimal adapter proving that the generic service owns no ENCODE fields."""

    metadata = WorkflowMetadata(
        workflow_id="contract-only",
        name="Contract-only workflow",
        version="1.0.0",
    )
    capabilities = WorkflowCapabilities(supports=("validation", "input_bundle_import"))

    def schema(self) -> WorkflowSchema:
        return WorkflowSchema(config_schema={"type": "object"})

    def validate(self, inputs: WorkflowInputs) -> Result[object]:
        return Result.success({"validated": inputs.to_dict()})

    def import_input_bundle(self, bundle) -> Result[InputBundleMapping]:
        relative_path = bundle.files[0].relative_path
        return Result.success(
            InputBundleMapping(
                inputs=WorkflowInputs(
                    config={"external_workflow": bundle.workflow_name},
                    samples=[{"sample": "imported"}],
                ),
                required_project_files=(relative_path,),
            )
        )

    def preview_dag(self, inputs: WorkflowInputs) -> Result[DagPreview]:
        return Result.success(DagPreview())

    def plan_workspace(
        self,
        inputs: WorkflowInputs,
        workspace: str,
    ) -> Result[WorkspacePlan]:
        return Result.success(WorkspacePlan(directories=(workspace,)))

    def build_command(self, plan: WorkspacePlan) -> Result[CommandSpec]:
        return Result.success(CommandSpec(argv=("run",)))

    def extract_artifacts(
        self, inputs: WorkflowInputs, workspace: str
    ) -> Result[tuple]:
        return Result.success(())


def _service() -> InputBundleImportService:
    return InputBundleImportService(WorkflowRegistry((EncodeStyleWorkflowAdapter(),)))


def _error_code(result) -> str:
    assert result.is_failure
    assert len(result.errors) == 1
    return result.errors[0].code


def _bowtie2_source_paths(result) -> tuple[str, ...]:
    assert result.is_success, result.issues
    return tuple(
        item.relative_path
        for item in result.value.source_files
        if item.relative_path.startswith("resources/mm10/index/genome.")
    )


def test_formal_v020_release_provenance_is_pinned() -> None:
    assert input_bundle_module.OMICS_INTAKE_BUNDLE_RELEASE_TAG == "v0.2.0"
    assert (
        input_bundle_module.OMICS_INTAKE_BUNDLE_RELEASE_TAG_OBJECT
        == "140a454d1313b19b322a825a1feebbb1494297c7"
    )
    assert (
        input_bundle_module.OMICS_INTAKE_BUNDLE_RELEASE_COMMIT
        == "32680c12465f543214ed7e0173c639e0d40c7113"
    )
    assert (
        input_bundle_module.OMICS_INTAKE_BUNDLE_RELEASE_TREE
        == "48aba2f48fa88fc37dab19c10f0ce70f2641add2"
    )
    assert input_bundle_module.OMICS_INTAKE_BUNDLE_CONTRACT_VERSION == "0.2"
    assert (
        input_bundle_module.OMICS_INTAKE_BUNDLE_SCHEMA_SHA256
        == "6dcda336c9f0ba763383ddd58bec280946f8970af8bd730eb758fab5e3a8dd71"
    )


def test_runnable_bundle_maps_to_fresh_validated_inputs_without_source_writes(
    tmp_path,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    before = tree_digest(bundle_path.parent)

    first = _service().inspect(bundle_path, WORKFLOW_ID)
    second = _service().inspect(bundle_path, WORKFLOW_ID)
    normalized = _service().inspect(bundle_path, f"  {WORKFLOW_ID}  ")

    assert first.is_success, first.issues
    assert second.is_success, second.issues
    assert normalized.is_success, normalized.issues
    assert first.value is not second.value
    assert first.value.inputs is not second.value.inputs
    assert first.value.inputs.to_dict() == second.value.inputs.to_dict()
    assert first.value.identity.contract_version == "0.2"
    assert first.value.identity.producer_name == "omics-intake"
    assert first.value.identity.producer_version == "0.2.0"
    assert first.value.workflow_id == WORKFLOW_ID
    assert normalized.value.workflow_id == WORKFLOW_ID
    assert first.value.inputs.options == {"strict_inputs": True}
    assert "samples" not in first.value.inputs.config
    assert "outdir" not in first.value.inputs.config
    row = first.value.inputs.samples[0]
    assert row["fastq_1"] == str(
        bundle_path.parent / "downloads/ena/SRR000001.fastq.gz"
    )
    assert row["bowtie2_index"] == str(
        bundle_path.parent / "resources/mm10/index/genome"
    )
    assert len(first.value.source_files) == 7
    assert sum(item.contract_bound for item in first.value.source_files) == 1
    assert tree_digest(bundle_path.parent) == before


@pytest.mark.parametrize(
    ("standard_count", "large_count", "selected_suffix"),
    (
        pytest.param(6, 0, ".bt2", id="standard"),
        pytest.param(0, 6, ".bt2l", id="large"),
        pytest.param(6, 6, ".bt2", id="standard-wins-when-both-complete"),
        pytest.param(3, 6, ".bt2l", id="large-wins-over-partial-standard"),
    ),
)
def test_service_safely_selects_complete_bowtie2_file_set(
    tmp_path,
    standard_count,
    large_count,
    selected_suffix,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    configure_bowtie2_shards(
        bundle_path,
        standard_count=standard_count,
        large_count=large_count,
    )

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    selected = _bowtie2_source_paths(result)
    assert len(selected) == 6
    assert all(path.endswith(selected_suffix) for path in selected)
    if selected_suffix == ".bt2":
        assert not any(path.endswith(".bt2l") for path in selected)


def test_two_incomplete_bowtie2_file_sets_fail_stably(tmp_path) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    configure_bowtie2_shards(
        bundle_path,
        standard_count=5,
        large_count=5,
    )

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == "INPUT_BUNDLE_SOURCE_UNSAFE"


@pytest.mark.parametrize("source_type", ("symlink", "fifo", "directory"))
def test_unsafe_standard_shard_never_falls_back_to_complete_large_set(
    tmp_path,
    source_type,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    configure_bowtie2_shards(
        bundle_path,
        standard_count=6,
        large_count=6,
    )
    target = bowtie2_shard_paths(bundle_path.parent)[0]
    target.unlink()
    if source_type == "symlink":
        outside = tmp_path / "outside.bt2"
        outside.write_bytes(b"tiny-index\n")
        target.symlink_to(outside)
    if source_type == "fifo":
        os.mkfifo(target)
    if source_type == "directory":
        target.mkdir()

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == "INPUT_BUNDLE_SOURCE_UNSAFE"


def test_missing_bowtie2_shard_without_complete_alternative_fails_closed(
    tmp_path,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    bowtie2_shard_paths(bundle_path.parent)[0].unlink()

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == "INPUT_BUNDLE_SOURCE_UNSAFE"


def test_bowtie2_candidate_rejects_symlinked_intermediate_directory(
    tmp_path,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    index_directory = bundle_path.parent / "resources/mm10/index"
    outside = tmp_path / "outside-index"
    index_directory.rename(outside)
    index_directory.symlink_to(outside, target_is_directory=True)

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == "INPUT_BUNDLE_SOURCE_UNSAFE"


def test_bowtie2_candidate_rejects_parent_traversal_before_resolution(
    tmp_path,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    samples_path = bundle_path.parent / "samples.encode.tsv"
    samples_path.write_text(
        samples_path.read_text(encoding="utf-8").replace(
            "resources/mm10/index/genome",
            "../outside/genome",
        ),
        encoding="utf-8",
    )
    refresh_artifact(bundle_path, "sample_sheet")

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == "ENCODE_INPUT_BUNDLE_PATH_UNSAFE"


def test_replacement_between_candidate_listing_and_open_fails_closed(
    tmp_path,
    monkeypatch,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    target = bowtie2_shard_paths(bundle_path.parent)[0]
    parked = tmp_path / "parked.bt2"
    body = target.read_bytes()
    real_open = input_bundle_module.os.open
    swapped = False

    def racing_open(path, flags, mode=0o777, *, dir_fd=None):
        nonlocal swapped
        if path == target.name and dir_fd is not None and not swapped:
            swapped = True
            target.rename(parked)
            target.write_bytes(body)
        return real_open(path, flags, mode, dir_fd=dir_fd)

    monkeypatch.setattr(input_bundle_module.os, "open", racing_open)

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    assert swapped is True
    assert _error_code(result) == "INPUT_BUNDLE_SOURCE_UNSAFE"


def test_validation_window_same_content_replacement_fails_closed(
    tmp_path,
    monkeypatch,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    target = bowtie2_shard_paths(bundle_path.parent)[0]
    parked = tmp_path / "validated-old.bt2"
    adapter = EncodeStyleWorkflowAdapter()
    real_validate = adapter.validate

    def replacing_validate(inputs):
        validation = real_validate(inputs)
        body = target.read_bytes()
        target.rename(parked)
        target.write_bytes(body)
        assert target.stat().st_ino != parked.stat().st_ino
        return validation

    monkeypatch.setattr(adapter, "validate", replacing_validate)
    service = InputBundleImportService(WorkflowRegistry((adapter,)))

    result = service.inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == "INPUT_BUNDLE_SOURCE_UNSAFE"


def test_validation_window_candidate_selection_change_fails_closed(
    tmp_path,
    monkeypatch,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    configure_bowtie2_shards(
        bundle_path,
        standard_count=5,
        large_count=6,
    )
    missing_standard = bowtie2_shard_paths(bundle_path.parent)[-1]
    adapter = EncodeStyleWorkflowAdapter()
    real_validate = adapter.validate

    def completing_validate(inputs):
        missing_standard.write_bytes(b"tiny-index\n")
        return real_validate(inputs)

    monkeypatch.setattr(adapter, "validate", completing_validate)
    service = InputBundleImportService(WorkflowRegistry((adapter,)))

    result = service.inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == "INPUT_BUNDLE_SOURCE_UNSAFE"


def test_workflow_neutral_service_dispatches_mapping_to_registered_adapter(
    tmp_path,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    adapter = ContractOnlyAdapter()
    service = InputBundleImportService(WorkflowRegistry((adapter,)))

    result = service.inspect(bundle_path, adapter.metadata.workflow_id)

    assert result.is_success, result.issues
    assert result.value.workflow_id == "contract-only"
    assert result.value.inputs.config == {"external_workflow": "encode-epigenomics"}
    assert result.value.inputs.samples == [{"sample": "imported"}]
    assert len(result.value.source_files) == 1
    assert result.value.source_files[0].contract_bound is True


@pytest.mark.parametrize(
    ("mutation", "expected_code"),
    [
        ("unknown-version", "INPUT_BUNDLE_CONTRACT_UNSUPPORTED"),
        ("extra-field", "INPUT_BUNDLE_DOCUMENT_INVALID"),
        ("not-runnable", "INPUT_BUNDLE_HANDOFF_NOT_READY"),
        ("stale-validation", "INPUT_BUNDLE_DOCUMENT_INVALID"),
        ("malformed-timestamp", "INPUT_BUNDLE_DOCUMENT_INVALID"),
        ("error-summary", "INPUT_BUNDLE_HANDOFF_NOT_READY"),
        ("needs-review-summary", "INPUT_BUNDLE_HANDOFF_NOT_READY"),
        ("transfer-summary", "INPUT_BUNDLE_HANDOFF_NOT_READY"),
        ("canonical-identity", "INPUT_BUNDLE_INTEGRITY_INVALID"),
        ("render-contract", "ENCODE_INPUT_BUNDLE_RENDER_UNSUPPORTED"),
    ],
)
def test_bundle_contract_and_adapter_mismatches_fail_closed(
    tmp_path,
    mutation,
    expected_code,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    payload = read_bundle(bundle_path)
    if mutation == "unknown-version":
        payload["bundle_schema_version"] = "0.3"
    elif mutation == "extra-field":
        payload["private"] = str(tmp_path / "secret")
    elif mutation == "not-runnable":
        payload["readiness"]["execution_status"] = "not_ready"
    elif mutation == "stale-validation":
        payload["validations"][1]["freshness"] = "stale"
    elif mutation == "malformed-timestamp":
        payload["validations"][0]["completed_at"] = "2026-07-15 00:00:00"
    elif mutation in {"error-summary", "needs-review-summary", "transfer-summary"}:
        severity = {
            "error-summary": "error",
            "needs-review-summary": "needs_review",
            "transfer-summary": "warning",
        }[mutation]
        payload["issues"] = [
            {
                "severity": severity,
                "code": "EXTERNAL_FILE_BINDING_REQUIRED",
                "owner_kind": "transfer"
                if mutation == "transfer-summary"
                else "data_file",
                "count": 1,
            }
        ]
    elif mutation == "canonical-identity":
        payload["canonical_project"]["identity"] = "sha256:" + "9" * 64
    else:
        payload["workflow"]["render_contract"] = "encode-render-v4"
    write_bundle(bundle_path, payload)

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == expected_code
    assert str(tmp_path) not in result.errors[0].message
    assert result.errors[0].context == {}


def test_duplicate_planned_destination_is_rejected_beyond_json_schema(tmp_path) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    payload = read_bundle(bundle_path)
    second_local = deepcopy(payload["files"][0])
    second_local["file_id"] = "insdc:SRR000002:1"
    second_planned = deepcopy(payload["files"][1])
    second_planned["file_id"] = "insdc:SRR000002:1"
    payload["files"].extend((second_local, second_planned))
    payload["files"].sort(
        key=lambda item: (
            item["file_id"],
            item["kind"],
            item["path"] or "",
            item["scope"],
        )
    )
    write_bundle(bundle_path, payload)

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == "INPUT_BUNDLE_DOCUMENT_INVALID"


def test_encode_adapter_rejects_unreferenced_bundle_fastq_binding(tmp_path) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    root = bundle_path.parent
    original = root / "downloads/ena/SRR000001.fastq.gz"
    extra_relative = "downloads/ena/SRR000002.fastq.gz"
    (root / extra_relative).write_bytes(original.read_bytes())
    payload = read_bundle(bundle_path)
    for record in tuple(payload["files"]):
        extra = deepcopy(record)
        extra["file_id"] = "insdc:SRR000002:1"
        extra["path"] = extra_relative
        payload["files"].append(extra)
    payload["files"].sort(
        key=lambda item: (
            item["file_id"],
            item["kind"],
            item["path"] or "",
            item["scope"],
        )
    )
    write_bundle(bundle_path, payload)

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == "ENCODE_INPUT_BUNDLE_REFERENCE_UNVERIFIED"


def test_encode_adapter_rejects_duplicate_fastq_consumption(tmp_path) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    samples_path = bundle_path.parent / "samples.encode.tsv"
    lines = samples_path.read_text(encoding="utf-8").splitlines()
    samples_path.write_text("\n".join((*lines, lines[-1])) + "\n", encoding="utf-8")
    refresh_artifact(bundle_path, "sample_sheet")

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == "ENCODE_INPUT_BUNDLE_REFERENCE_UNVERIFIED"


def test_changed_artifact_bytes_are_rejected_before_adapter_mapping(tmp_path) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    (bundle_path.parent / "samples.encode.tsv").write_text(
        "changed private content\n", encoding="utf-8"
    )

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == "INPUT_BUNDLE_INTEGRITY_INVALID"


def test_absolute_path_inside_digest_valid_tsv_is_rejected_without_echo(
    tmp_path,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    samples_path = bundle_path.parent / "samples.encode.tsv"
    text = samples_path.read_text(encoding="utf-8")
    private_path = str(tmp_path / "private" / "read.fastq.gz")
    samples_path.write_text(
        text.replace("downloads/ena/SRR000001.fastq.gz", private_path),
        encoding="utf-8",
    )
    refresh_artifact(bundle_path, "sample_sheet")

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == "ENCODE_INPUT_BUNDLE_PATH_UNSAFE"
    assert private_path not in result.errors[0].message
    assert result.errors[0].context == {}


@pytest.mark.parametrize(
    ("mutation", "expected_code"),
    [
        ("missing-binding", "ENCODE_INPUT_BUNDLE_REFERENCE_UNVERIFIED"),
        ("binding-role", "ENCODE_INPUT_BUNDLE_REFERENCE_UNVERIFIED"),
        ("sample-genome", "ENCODE_INPUT_BUNDLE_INVALID"),
        ("sample-header", "ENCODE_INPUT_BUNDLE_INVALID"),
        ("config-samples", "ENCODE_INPUT_BUNDLE_INVALID"),
        ("config-outdir", "ENCODE_INPUT_BUNDLE_INVALID"),
        ("config-not-mapping", "ENCODE_INPUT_BUNDLE_INVALID"),
        ("config-duplicate-key", "ENCODE_INPUT_BUNDLE_INVALID"),
        ("config-alias", "ENCODE_INPUT_BUNDLE_INVALID"),
        ("config-too-deep", "ENCODE_INPUT_BUNDLE_INVALID"),
        ("missing-index", "INPUT_BUNDLE_SOURCE_UNSAFE"),
    ],
)
def test_encode_render_contract_mapping_rejects_incomplete_or_mismatched_inputs(
    tmp_path,
    mutation,
    expected_code,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    root = bundle_path.parent
    if mutation in {"missing-binding", "sample-genome", "sample-header"}:
        samples_path = root / "samples.encode.tsv"
        text = samples_path.read_text(encoding="utf-8")
        if mutation == "missing-binding":
            text = text.replace(
                "downloads/ena/SRR000001.fastq.gz",
                "downloads/ena/UNBOUND.fastq.gz",
            )
        elif mutation == "sample-genome":
            text = text.replace("\tmm10\t", "\thg38\t")
        else:
            text = text.replace("fastq_1", "read_one", 1)
        samples_path.write_text(text, encoding="utf-8")
        refresh_artifact(bundle_path, "sample_sheet")
    elif mutation == "binding-role":
        payload = read_bundle(bundle_path)
        for record in payload["files"]:
            record["role"] = "control"
        write_bundle(bundle_path, payload)
    elif mutation.startswith("config-"):
        config_path = root / "config.encode.yaml"
        if mutation == "config-not-mapping":
            config_path.write_text("- invalid\n", encoding="utf-8")
        elif mutation == "config-duplicate-key":
            config_path.write_text(
                config_path.read_text(encoding="utf-8") + "threads: 8\n",
                encoding="utf-8",
            )
        elif mutation == "config-alias":
            config_path.write_text(
                config_path.read_text(encoding="utf-8")
                + "recursive: &recursive [*recursive]\n",
                encoding="utf-8",
            )
        elif mutation == "config-too-deep":
            config_path.write_text(
                config_path.read_text(encoding="utf-8")
                + "nested: "
                + "[" * 70
                + "0"
                + "]" * 70
                + "\n",
                encoding="utf-8",
            )
        else:
            config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
            if mutation == "config-samples":
                config["samples"] = "samples.tsv"
            else:
                config["outdir"] = "private-results"
            config_path.write_text(yaml.safe_dump(config), encoding="utf-8")
        refresh_artifact(bundle_path, "workflow_config")
    else:
        (root / "resources/mm10/index/genome.1.bt2").unlink()

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == expected_code
    assert str(tmp_path) not in result.errors[0].message
    assert result.errors[0].context == {}


def test_symlinked_artifact_and_bundle_are_rejected(tmp_path) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    samples_path = bundle_path.parent / "samples.encode.tsv"
    outside = tmp_path / "outside.tsv"
    outside.write_bytes(samples_path.read_bytes())
    samples_path.unlink()
    samples_path.symlink_to(outside)

    artifact_result = _service().inspect(bundle_path, WORKFLOW_ID)
    assert _error_code(artifact_result) == "INPUT_BUNDLE_SOURCE_UNSAFE"

    linked_root = tmp_path / "linked"
    linked_root.mkdir()
    linked_bundle = linked_root / "intake-bundle.json"
    linked_bundle.symlink_to(bundle_path)
    bundle_result = _service().inspect(linked_bundle, WORKFLOW_ID)
    assert _error_code(bundle_result) == "INPUT_BUNDLE_SOURCE_UNSAFE"


def test_missing_workflow_and_undeclared_capability_fail_without_side_effects(
    tmp_path,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    before = tree_digest(bundle_path.parent)

    missing = _service().inspect(bundle_path, "missing")

    assert _error_code(missing) == "INPUT_BUNDLE_WORKFLOW_NOT_FOUND"
    assert tree_digest(bundle_path.parent) == before

    adapter = ContractOnlyAdapter()
    adapter.capabilities = WorkflowCapabilities(supports=("validation",))
    with pytest.raises(ValueError, match="input_bundle_import"):
        WorkflowRegistry((adapter,))


@pytest.mark.parametrize(
    ("body", "expected_code"),
    [
        (
            b'{"bundle_schema_version":"0.2","value":NaN}',
            "INPUT_BUNDLE_DOCUMENT_INVALID",
        ),
        (b"\xff", "INPUT_BUNDLE_DOCUMENT_INVALID"),
        (
            b'{"bundle_schema_version":"0.2","bundle_schema_version":"0.2"}',
            "INPUT_BUNDLE_DOCUMENT_INVALID",
        ),
    ],
)
def test_non_canonical_json_encodings_are_rejected(
    tmp_path, body, expected_code
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    bundle_path.write_bytes(body)

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == expected_code


def test_oversized_bundle_and_fifo_source_file_are_rejected_without_blocking(
    tmp_path,
) -> None:
    oversized_path = create_runnable_bundle(tmp_path / "oversized")
    with oversized_path.open("wb") as handle:
        handle.truncate(16 * 1024 * 1024 + 1)
    assert _error_code(_service().inspect(oversized_path, WORKFLOW_ID)) == (
        "INPUT_BUNDLE_TOO_LARGE"
    )

    fifo_path = create_runnable_bundle(tmp_path / "fifo")
    fastq = fifo_path.parent / "downloads/ena/SRR000001.fastq.gz"
    fastq.unlink()
    os.mkfifo(fastq)

    assert _error_code(_service().inspect(fifo_path, WORKFLOW_ID)) == (
        "INPUT_BUNDLE_SOURCE_UNSAFE"
    )


def test_oversized_required_file_is_rejected_before_content_read(tmp_path) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    fastq = bundle_path.parent / "downloads/ena/SRR000001.fastq.gz"
    with fastq.open("wb") as handle:
        handle.truncate(256 * 1024**3 + 1)

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == "INPUT_BUNDLE_TOO_LARGE"


def test_bundle_inspection_does_not_load_omics_intake_python_modules(tmp_path) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    before = {
        name
        for name in sys.modules
        if name == "omics_intake" or name.startswith("omics_intake.")
    }

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    after = {
        name
        for name in sys.modules
        if name == "omics_intake" or name.startswith("omics_intake.")
    }
    assert result.is_success, result.issues
    assert after == before
