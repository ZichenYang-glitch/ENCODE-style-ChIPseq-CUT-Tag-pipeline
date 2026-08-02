"""Tests for fail-closed default bulk RNA-seq operator composition."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import replace
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
from encode_pipeline.adapters.bulk_rnaseq.execution_identity import (
    verify_execution_implementation,
)
from encode_pipeline.adapters.bulk_rnaseq.qualification import (
    load_default_execution_qualification,
)
from encode_pipeline.platform.adapters import WorkflowAvailability
from encode_pipeline.platform.results import Result
from encode_pipeline.services.defaults import create_default_workflow_registry
from encode_pipeline.services.workflow_info import WorkflowInfoService
from encode_pipeline.platform.registry import WorkflowRegistry


UNTRUSTED_QUALIFICATION_ENV = "HELIXWEAVE_BULK_RNASEQ_EXECUTION_QUALIFIED"


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


def test_missing_source_candidate_cannot_be_overridden_by_environment(
    tmp_path,
    monkeypatch,
):
    manifest_reads = 0

    def record_manifest_read(_path):
        nonlocal manifest_reads
        manifest_reads += 1
        raise AssertionError("missing qualification must fail before private reads")

    monkeypatch.setattr(
        "encode_pipeline.adapters.bulk_rnaseq.deployment._load_transcriptome_binding",
        record_manifest_read,
    )

    def missing_candidate():
        raise OSError("source candidate is missing")

    monkeypatch.setattr(
        "encode_pipeline.adapters.bulk_rnaseq.qualification."
        "_read_qualification_resource",
        missing_candidate,
    )
    environment = _environment(tmp_path)
    environment[UNTRUSTED_QUALIFICATION_ENV] = "true"
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


def test_stale_source_candidate_checks_keys_without_reading_private_values(
    monkeypatch,
):
    implementation = verify_execution_implementation()
    assert implementation.is_success
    stale = replace(implementation.value, aggregate_sha256="0" * 64)
    monkeypatch.setattr(
        "encode_pipeline.adapters.bulk_rnaseq.deployment."
        "verify_execution_implementation",
        lambda: Result.success(stale),
    )

    class KeysOnlyCoordinates(Mapping[str, str]):
        def __iter__(self) -> Iterator[str]:
            return iter(
                (
                    RUNTIME_ROOT_ENV,
                    TRANSCRIPTOME_BINDING_MANIFEST_ENV,
                    MANAGED_DOCKER_EXECUTABLE_ENV,
                    MANAGED_DOCKER_SOCKET_ENV,
                    UNTRUSTED_QUALIFICATION_ENV,
                )
            )

        def __len__(self) -> int:
            return 5

        def __getitem__(self, _key: str) -> str:
            raise AssertionError("stale qualification must not read private values")

    adapter = load_default_bulk_rnaseq_adapter(KeysOnlyCoordinates())

    assert adapter.execution_availability().execution == "unavailable"
    assert local_execution_configuration(adapter) is None


def test_source_enabled_complete_coordinates_project_ready_default_registry(
    tmp_path,
    monkeypatch,
):
    monkeypatch.setattr(
        BulkRnaSeqResultsWorkflowAdapter,
        "execution_availability",
        lambda _self: WorkflowAvailability(),
    )

    registry = create_default_workflow_registry(environ=_environment(tmp_path))
    descriptor = WorkflowInfoService(registry).get_descriptor("bulk-rnaseq")

    assert descriptor.is_success
    assert descriptor.value is not None
    assert descriptor.value.availability.to_dict() == {
        "authoring": "available",
        "execution": "available",
        "reason_code": "WORKFLOW_EXECUTION_READY",
    }
    assert {
        "workspace_plan",
        "command",
        "artifact_extract",
        "qc_summary_extract",
    } <= set(descriptor.value.capabilities.supports)


def test_malformed_binding_manifest_fails_closed(tmp_path):
    environment = _environment(tmp_path)
    manifest = tmp_path / "transcriptome-binding.json"
    manifest.write_text('{"schema_version":"1.0.0","private_path":"/secret"}')

    adapter = load_default_bulk_rnaseq_adapter(environment)

    assert adapter.execution_availability().execution == "unavailable"
    assert adapter.capabilities.supports == ("validation", "input_authoring")


def test_exact_source_candidate_composes_adapter_that_can_report_ready(
    tmp_path,
    monkeypatch,
):
    implementation = verify_execution_implementation()
    assert implementation.is_success
    qualification = load_default_execution_qualification(implementation.value)
    assert qualification.is_success
    assert qualification.value.implementation.matches(implementation.value)

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
    assert (
        adapter.execution_binding.implementation_qualification
        == qualification.value.implementation
    )


def test_failed_admission_is_rechecked_without_declaring_public_execution_capabilities(
    tmp_path,
    monkeypatch,
):
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
