"""Tests for fail-closed default bulk RNA-seq operator composition."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
import json

from encode_pipeline.adapters.bulk_rnaseq import BulkRnaSeqResultsWorkflowAdapter
from encode_pipeline.adapters.bulk_rnaseq.deployment import (
    MANAGED_DOCKER_EXECUTABLE_ENV,
    MANAGED_DOCKER_SOCKET_ENV,
    RUNTIME_ROOT_ENV,
    TRANSCRIPTOME_BINDING_MANIFEST_ENV,
    local_execution_configuration,
    load_default_bulk_rnaseq_adapter,
)
from encode_pipeline.platform.adapters import WorkflowAvailability
from encode_pipeline.services.workflow_info import WorkflowInfoService
from encode_pipeline.platform.registry import WorkflowRegistry


def _manifest(tmp_path):
    path = (tmp_path / "transcriptome-binding.json").resolve()
    path.write_text(
        json.dumps(
            {
                "schema_version": "1.0.0",
                "reference_id": "tiny",
                "fasta_sha256": "a" * 64,
                "gtf_sha256": "b" * 64,
                "transcript_fasta": str((tmp_path / "transcripts.fa").resolve()),
                "transcript_fasta_sha256": "c" * 64,
            }
        ),
        encoding="utf-8",
    )
    return path


def _environment(tmp_path):
    return {
        RUNTIME_ROOT_ENV: str((tmp_path / "runtime").resolve()),
        TRANSCRIPTOME_BINDING_MANIFEST_ENV: str(_manifest(tmp_path)),
        MANAGED_DOCKER_EXECUTABLE_ENV: "/usr/bin/docker",
        MANAGED_DOCKER_SOCKET_ENV: str((tmp_path / "docker.sock").resolve()),
    }


def _enable_exact_head_qualification(monkeypatch):
    monkeypatch.setattr(
        "encode_pipeline.adapters.bulk_rnaseq.deployment."
        "_DEFAULT_EXECUTION_EXACT_HEAD_QUALIFIED",
        True,
    )


def test_absent_coordinates_keep_authoring_available_and_execution_not_configured():
    adapter = load_default_bulk_rnaseq_adapter({})

    assert adapter.capabilities.supports == ("validation", "input_authoring")
    assert adapter.execution_availability().to_dict() == {
        "authoring": "available",
        "execution": "not_configured",
        "reason_code": "WORKFLOW_EXECUTION_NOT_CONFIGURED",
    }


def test_partial_coordinates_fail_closed_without_exposing_coordinate_values(tmp_path):
    private_root = str((tmp_path / "private-runtime").resolve())

    adapter = load_default_bulk_rnaseq_adapter({RUNTIME_ROOT_ENV: private_root})

    availability = adapter.execution_availability()
    assert availability.execution == "unavailable"
    assert private_root not in repr(adapter)
    assert private_root not in repr(availability)


def test_complete_coordinates_remain_unavailable_until_exact_head_qualification(
    tmp_path,
    monkeypatch,
):
    manifest_reads = 0

    def record_manifest_read(_path):
        nonlocal manifest_reads
        manifest_reads += 1
        raise AssertionError("pending qualification must fail before private reads")

    monkeypatch.setattr(
        "encode_pipeline.adapters.bulk_rnaseq.deployment._load_transcriptome_binding",
        record_manifest_read,
    )
    environment = _environment(tmp_path)
    private_root = environment[RUNTIME_ROOT_ENV]

    adapter = load_default_bulk_rnaseq_adapter(environment)

    assert manifest_reads == 0
    assert not isinstance(adapter, BulkRnaSeqResultsWorkflowAdapter)
    assert adapter.capabilities.supports == ("validation", "input_authoring")
    assert adapter.execution_availability().to_dict() == {
        "authoring": "available",
        "execution": "unavailable",
        "reason_code": "WORKFLOW_EXECUTION_UNAVAILABLE",
    }
    assert local_execution_configuration(adapter) is None
    assert private_root not in repr(adapter)


def test_pending_qualification_checks_coordinate_keys_without_reading_values():
    class KeysOnlyCoordinates(Mapping[str, str]):
        def __iter__(self) -> Iterator[str]:
            return iter(
                (
                    RUNTIME_ROOT_ENV,
                    TRANSCRIPTOME_BINDING_MANIFEST_ENV,
                    MANAGED_DOCKER_EXECUTABLE_ENV,
                    MANAGED_DOCKER_SOCKET_ENV,
                )
            )

        def __len__(self) -> int:
            return 4

        def __getitem__(self, _key: str) -> str:
            raise AssertionError("pending qualification must not read private values")

    adapter = load_default_bulk_rnaseq_adapter(KeysOnlyCoordinates())

    assert adapter.execution_availability().execution == "unavailable"
    assert local_execution_configuration(adapter) is None


def test_malformed_binding_manifest_fails_closed(tmp_path, monkeypatch):
    _enable_exact_head_qualification(monkeypatch)
    environment = _environment(tmp_path)
    manifest = tmp_path / "transcriptome-binding.json"
    manifest.write_text('{"schema_version":"1.0.0","private_path":"/secret"}')

    adapter = load_default_bulk_rnaseq_adapter(environment)

    assert adapter.execution_availability().execution == "unavailable"
    assert adapter.capabilities.supports == ("validation", "input_authoring")


def test_complete_coordinates_declare_execution_only_after_admission(
    tmp_path,
    monkeypatch,
):
    _enable_exact_head_qualification(monkeypatch)
    monkeypatch.setattr(
        BulkRnaSeqResultsWorkflowAdapter,
        "execution_availability",
        lambda _self: WorkflowAvailability(),
    )

    adapter = load_default_bulk_rnaseq_adapter(_environment(tmp_path))

    assert isinstance(adapter, BulkRnaSeqResultsWorkflowAdapter)
    assert adapter.capabilities.supports == (
        "validation",
        "input_authoring",
        "workspace_plan",
        "command",
        "artifact_extract",
        "qc_summary_extract",
    )
    assert adapter.execution_binding is not None
    assert adapter.execution_binding.assets.root == (tmp_path / "runtime").resolve()


def test_failed_admission_is_rechecked_without_declaring_public_execution_capabilities(
    tmp_path,
    monkeypatch,
):
    _enable_exact_head_qualification(monkeypatch)
    ready = False

    def current_availability(_self):
        if ready:
            return WorkflowAvailability()
        return WorkflowAvailability(
            execution="unavailable",
            reason_code="WORKFLOW_EXECUTION_UNAVAILABLE",
        )

    monkeypatch.setattr(
        BulkRnaSeqResultsWorkflowAdapter,
        "execution_availability",
        current_availability,
    )

    adapter = load_default_bulk_rnaseq_adapter(_environment(tmp_path))

    assert isinstance(adapter, BulkRnaSeqResultsWorkflowAdapter)
    assert adapter.execution_availability().execution == "unavailable"
    descriptor = WorkflowInfoService(WorkflowRegistry((adapter,))).get_descriptor(
        "bulk-rnaseq"
    )
    assert descriptor.is_success
    assert descriptor.value.capabilities.supports == (
        "validation",
        "input_authoring",
    )

    ready = True
    recovered = WorkflowInfoService(WorkflowRegistry((adapter,))).get_descriptor(
        "bulk-rnaseq"
    )
    assert recovered.is_success
    assert recovered.value.availability.execution == "available"
    assert recovered.value.capabilities.supports == adapter.capabilities.supports


def test_unexpected_admission_error_fails_closed(tmp_path, monkeypatch):
    _enable_exact_head_qualification(monkeypatch)

    def raise_private_error(_path):
        raise RuntimeError("private runtime coordinate must not escape")

    monkeypatch.setattr(
        "encode_pipeline.adapters.bulk_rnaseq.deployment._load_transcriptome_binding",
        raise_private_error,
    )

    adapter = load_default_bulk_rnaseq_adapter(_environment(tmp_path))

    assert not isinstance(adapter, BulkRnaSeqResultsWorkflowAdapter)
    assert adapter.execution_availability().to_dict() == {
        "authoring": "available",
        "execution": "unavailable",
        "reason_code": "WORKFLOW_EXECUTION_UNAVAILABLE",
    }
