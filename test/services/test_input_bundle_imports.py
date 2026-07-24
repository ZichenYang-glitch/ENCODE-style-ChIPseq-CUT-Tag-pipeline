"""Read-only Omics Intake Bundle contract consumption tests."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path
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
from encode_pipeline.platform.results import Issue, Result
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


@pytest.mark.parametrize(
    ("failure_boundary", "expected_code"),
    (
        pytest.param(
            "capability-drift",
            "INPUT_BUNDLE_CAPABILITY_UNSUPPORTED",
            id="capability-drift",
        ),
        pytest.param(
            "contract-loader",
            "INPUT_BUNDLE_DOCUMENT_INVALID",
            id="unexpected-contract-loader-failure",
        ),
        pytest.param(
            "import-raises",
            "INPUT_BUNDLE_ADAPTER_FAILED",
            id="importer-raises",
        ),
        pytest.param(
            "import-non-result",
            "INPUT_BUNDLE_ADAPTER_FAILED",
            id="importer-non-result",
        ),
        pytest.param(
            "import-wrong-value",
            "INPUT_BUNDLE_ADAPTER_FAILED",
            id="importer-wrong-value",
        ),
        pytest.param(
            "observe-raises",
            "INPUT_BUNDLE_SOURCE_UNSAFE",
            id="observation-raises",
        ),
        pytest.param(
            "validate-raises",
            "INPUT_BUNDLE_ADAPTER_FAILED",
            id="validator-raises",
        ),
        pytest.param(
            "validate-non-result",
            "INPUT_BUNDLE_ADAPTER_FAILED",
            id="validator-non-result",
        ),
        pytest.param(
            "validate-failure",
            "CONTROLLED_VALIDATION_FAILURE",
            id="validator-failure",
        ),
    ),
)
def test_import_service_fail_closed_boundary_matrix(
    tmp_path,
    monkeypatch,
    failure_boundary,
    expected_code,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    adapter = ContractOnlyAdapter()
    registry = WorkflowRegistry((adapter,))
    service = InputBundleImportService(registry)

    if failure_boundary == "capability-drift":
        adapter.capabilities = WorkflowCapabilities(supports=("validation",))
    elif failure_boundary == "contract-loader":
        monkeypatch.setattr(
            input_bundle_module,
            "_load_contract_view",
            lambda _body: (_ for _ in ()).throw(RuntimeError("private parser failure")),
        )
    elif failure_boundary == "import-raises":
        monkeypatch.setattr(
            adapter,
            "import_input_bundle",
            lambda _bundle: (_ for _ in ()).throw(
                RuntimeError("private adapter failure")
            ),
        )
    elif failure_boundary == "import-non-result":
        monkeypatch.setattr(
            adapter,
            "import_input_bundle",
            lambda _bundle: object(),
        )
    elif failure_boundary == "import-wrong-value":
        monkeypatch.setattr(
            adapter,
            "import_input_bundle",
            lambda _bundle: Result.success(object()),
        )
    elif failure_boundary == "observe-raises":
        monkeypatch.setattr(
            input_bundle_module,
            "_observe_mapping_files",
            lambda *_args: (_ for _ in ()).throw(
                RuntimeError("private observation failure")
            ),
        )
    elif failure_boundary == "validate-raises":
        monkeypatch.setattr(
            adapter,
            "validate",
            lambda _inputs: (_ for _ in ()).throw(
                RuntimeError("private validator failure")
            ),
        )
    elif failure_boundary == "validate-non-result":
        monkeypatch.setattr(adapter, "validate", lambda _inputs: object())
    else:
        monkeypatch.setattr(
            adapter,
            "validate",
            lambda _inputs: Result.failure(
                [
                    Issue(
                        code="CONTROLLED_VALIDATION_FAILURE",
                        message="Input validation failed.",
                        source="test",
                    )
                ]
            ),
        )

    result = service.inspect(bundle_path, adapter.metadata.workflow_id)

    assert _error_code(result) == expected_code
    assert str(tmp_path) not in result.errors[0].message
    assert result.errors[0].context == {}


@pytest.mark.parametrize(
    "unsafe_source",
    (
        pytest.param("not-a-path", id="non-path"),
        pytest.param("wrong-name", id="wrong-filename"),
        pytest.param("directory-open-error", id="directory-open-error"),
    ),
)
def test_import_service_rejects_unsafe_source_coordinates(
    tmp_path,
    monkeypatch,
    unsafe_source,
) -> None:
    if unsafe_source == "not-a-path":
        bundle_path = "intake-bundle.json"
    elif unsafe_source == "wrong-name":
        bundle_path = tmp_path / "bundle.json"
    else:
        bundle_path = tmp_path / "intake-bundle.json"
        monkeypatch.setattr(
            input_bundle_module,
            "_open_absolute_directory",
            lambda _path: (_ for _ in ()).throw(OSError("private source open failure")),
        )

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == "INPUT_BUNDLE_SOURCE_UNSAFE"
    assert result.errors[0].context == {}


@pytest.mark.parametrize(
    ("mutation", "expected_code"),
    (
        pytest.param(
            "artifact-order",
            "INPUT_BUNDLE_DOCUMENT_INVALID",
            id="artifact-order",
        ),
        pytest.param(
            "file-count-limit",
            "INPUT_BUNDLE_TOO_LARGE",
            id="file-count-limit",
        ),
        pytest.param(
            "unbound-file",
            "INPUT_BUNDLE_HANDOFF_NOT_READY",
            id="unbound-file",
        ),
        pytest.param(
            "unverified-local",
            "INPUT_BUNDLE_HANDOFF_NOT_READY",
            id="unverified-local",
        ),
        pytest.param(
            "no-local-files",
            "INPUT_BUNDLE_HANDOFF_NOT_READY",
            id="no-local-files",
        ),
        pytest.param(
            "planned-path-mismatch",
            "INPUT_BUNDLE_INTEGRITY_INVALID",
            id="planned-path-mismatch",
        ),
        pytest.param(
            "validation-order",
            "INPUT_BUNDLE_DOCUMENT_INVALID",
            id="validation-order",
        ),
        pytest.param(
            "validation-not-ready",
            "INPUT_BUNDLE_HANDOFF_NOT_READY",
            id="validation-not-ready",
        ),
        pytest.param(
            "issue-count-limit",
            "INPUT_BUNDLE_TOO_LARGE",
            id="issue-count-limit",
        ),
        pytest.param(
            "issue-order",
            "INPUT_BUNDLE_DOCUMENT_INVALID",
            id="issue-order",
        ),
    ),
)
def test_bundle_semantics_rejects_noncanonical_handoff_states(
    tmp_path,
    monkeypatch,
    mutation,
    expected_code,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    payload = read_bundle(bundle_path)
    if mutation == "artifact-order":
        payload["artifacts"].reverse()
    elif mutation == "file-count-limit":
        monkeypatch.setattr(input_bundle_module, "_MAX_FILE_REFERENCES", 1)
    elif mutation == "unbound-file":
        payload["files"][0]["kind"] = "unbound"
        payload["files"][0]["scope"] = "external"
        payload["files"][0]["path"] = None
        payload["files"].sort(
            key=lambda item: (
                item["file_id"],
                item["kind"],
                item["path"] or "",
                item["scope"],
            )
        )
    elif mutation == "unverified-local":
        payload["files"][0]["checksum"]["status"] = "declared"
    elif mutation == "no-local-files":
        payload["files"] = [
            record for record in payload["files"] if record["kind"] == "planned"
        ]
    elif mutation == "planned-path-mismatch":
        payload["files"][1]["path"] = "downloads/ena/OTHER.fastq.gz"
    elif mutation == "validation-order":
        payload["validations"].reverse()
    elif mutation == "validation-not-ready":
        payload["validations"][0]["reason_code"] = "CONTROLLED_FAILURE"
    elif mutation == "issue-count-limit":
        monkeypatch.setattr(input_bundle_module, "_MAX_ISSUE_SUMMARIES", 0)
        payload["issues"] = [
            {
                "severity": "warning",
                "code": "CONTROLLED_WARNING",
                "owner_kind": "data_file",
                "count": 1,
            }
        ]
    else:
        payload["issues"] = [
            {
                "severity": "warning",
                "code": code,
                "owner_kind": "data_file",
                "count": 1,
            }
            for code in ("Z_WARNING", "A_WARNING")
        ]

    with pytest.raises(input_bundle_module._InputBundleFailure) as raised:
        input_bundle_module._validate_semantics(payload)

    assert raised.value.code == expected_code


@pytest.mark.parametrize(
    "body",
    (
        pytest.param(b"[]", id="top-level-array"),
        pytest.param(str(2**63).encode("ascii"), id="integer-out-of-range"),
        pytest.param(
            b"[" * 65 + b"0" + b"]" * 65,
            id="nesting-out-of-range",
        ),
    ),
)
def test_contract_loader_rejects_noncanonical_json_shapes(
    tmp_path,
    body,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    bundle_path.write_bytes(body)

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == "INPUT_BUNDLE_DOCUMENT_INVALID"


def test_contract_loader_reports_pinned_schema_unavailability(
    tmp_path,
    monkeypatch,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    monkeypatch.setattr(
        input_bundle_module,
        "_bundle_validator",
        lambda: (_ for _ in ()).throw(RuntimeError("schema unavailable")),
    )

    result = _service().inspect(bundle_path, WORKFLOW_ID)

    assert _error_code(result) == "INPUT_BUNDLE_CONTRACT_UNAVAILABLE"


@pytest.mark.parametrize(
    "schema_failure",
    (
        pytest.param("digest", id="schema-digest"),
        pytest.param("coordinate", id="schema-coordinate"),
    ),
)
def test_pinned_schema_identity_must_match_bundled_contract(
    monkeypatch,
    schema_failure,
) -> None:
    if schema_failure == "digest":
        schema_bytes = b"{}"
        expected_message = "digest mismatch"
    else:
        schema_bytes = json.dumps(
            {
                "$id": "urn:wrong-contract",
                "$schema": "https://json-schema.org/draft/2020-12/schema",
            }
        ).encode("utf-8")
        monkeypatch.setattr(
            input_bundle_module,
            "OMICS_INTAKE_BUNDLE_SCHEMA_SHA256",
            hashlib.sha256(schema_bytes).hexdigest(),
        )
        expected_message = "coordinate mismatch"

    class SchemaResource:
        def joinpath(self, _resource):
            return self

        def read_bytes(self):
            return schema_bytes

    monkeypatch.setattr(
        input_bundle_module,
        "files",
        lambda _package: SchemaResource(),
    )
    input_bundle_module._bundle_validator.cache_clear()

    with pytest.raises(ValueError, match=expected_message):
        input_bundle_module._bundle_validator()

    input_bundle_module._bundle_validator.cache_clear()


def test_invalid_calendar_timestamp_is_not_rfc3339() -> None:
    assert not input_bundle_module._is_rfc3339_datetime("2026-02-30T00:00:00Z")


def test_verified_bundle_rejects_nonmapping_workflow_identity(tmp_path) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    payload = read_bundle(bundle_path)
    contract = input_bundle_module._validate_semantics(payload)
    malformed = replace(
        contract,
        payload={**payload, "workflow": "private-invalid-workflow"},
    )

    with input_bundle_module._BundleProjectReader.open(bundle_path) as reader:
        with pytest.raises(input_bundle_module._InputBundleFailure) as raised:
            input_bundle_module._load_verified_bundle(
                reader,
                malformed,
                bundle_sha256="a" * 64,
            )

    assert raised.value.code == "INPUT_BUNDLE_HANDOFF_NOT_READY"


def test_reader_lifecycle_and_root_binding_fail_closed(
    tmp_path,
    monkeypatch,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    reader = input_bundle_module._BundleProjectReader.open(bundle_path)
    reader.close()
    reader.close()

    reader = input_bundle_module._BundleProjectReader.open(bundle_path)
    monkeypatch.setattr(
        input_bundle_module,
        "_open_absolute_directory",
        lambda _path: (_ for _ in ()).throw(OSError("root disappeared")),
    )
    with pytest.raises(input_bundle_module._InputBundleFailure) as raised:
        reader.verify_root_binding()
    reader.close()

    assert raised.value.code == "INPUT_BUNDLE_SOURCE_UNSAFE"


def test_reader_rejects_rebound_project_root(tmp_path, monkeypatch) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    other_root = tmp_path / "other-root"
    other_root.mkdir()
    other_descriptor = input_bundle_module._open_absolute_directory(other_root)
    reader = input_bundle_module._BundleProjectReader.open(bundle_path)
    monkeypatch.setattr(
        input_bundle_module,
        "_open_absolute_directory",
        lambda _path: os.dup(other_descriptor),
    )
    try:
        with pytest.raises(input_bundle_module._InputBundleFailure) as raised:
            reader.verify_root_binding()
    finally:
        reader.close()
        os.close(other_descriptor)

    assert raised.value.code == "INPUT_BUNDLE_SOURCE_UNSAFE"


def test_reader_rejects_invalid_paths_digests_and_reopen_identity(
    tmp_path,
    monkeypatch,
) -> None:
    bundle_path = create_runnable_bundle(tmp_path / "bundle")
    relative_path = "downloads/ena/SRR000001.fastq.gz"
    size_bytes = (bundle_path.parent / relative_path).stat().st_size

    with input_bundle_module._BundleProjectReader.open(bundle_path) as reader:
        with pytest.raises(input_bundle_module._InputBundleFailure) as invalid_path:
            reader.observe_file(
                "../private.fastq.gz",
                expected_size=None,
                expected_sha256=None,
                maximum_bytes=1024,
            )
        assert invalid_path.value.code == "INPUT_BUNDLE_DOCUMENT_INVALID"

        with pytest.raises(input_bundle_module._InputBundleFailure) as wrong_digest:
            reader.observe_file(
                relative_path,
                expected_size=size_bytes,
                expected_sha256="0" * 64,
                maximum_bytes=1024,
            )
        assert wrong_digest.value.code == "INPUT_BUNDLE_INTEGRITY_INVALID"

        monkeypatch.setattr(reader, "_reopen_identity", lambda _path: ())
        with pytest.raises(input_bundle_module._InputBundleFailure) as rebound:
            reader.observe_file(
                relative_path,
                expected_size=size_bytes,
                expected_sha256=None,
                maximum_bytes=1024,
            )
        assert rebound.value.code == "INPUT_BUNDLE_SOURCE_UNSAFE"


def test_safe_open_primitives_reject_unsupported_or_non_directory_paths(
    tmp_path,
    monkeypatch,
) -> None:
    with pytest.raises(ValueError, match="absolute and normalized"):
        input_bundle_module._open_absolute_directory(
            Path(f"{tmp_path}/../{tmp_path.name}")
        )

    with monkeypatch.context() as scoped:
        scoped.setattr(
            input_bundle_module,
            "_directory_flags",
            lambda: os.O_RDONLY | os.O_NOFOLLOW,
        )
        with pytest.raises(OSError, match="not a directory"):
            input_bundle_module._open_absolute_directory(Path("/dev/null"))

    with monkeypatch.context() as scoped:
        scoped.delattr(os, "O_DIRECTORY")
        with pytest.raises(OSError, match="directory flags"):
            input_bundle_module._directory_flags()

    with monkeypatch.context() as scoped:
        scoped.delattr(os, "O_NOFOLLOW")
        with pytest.raises(OSError, match="file flags"):
            input_bundle_module._file_flags()
