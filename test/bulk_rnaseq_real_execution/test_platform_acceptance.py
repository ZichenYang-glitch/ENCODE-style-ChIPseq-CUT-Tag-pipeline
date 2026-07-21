"""Explicit full STAR+Salmon acceptance through SQLite, Redis, RQ, and worker."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from .platform_harness import PlatformAcceptanceHarness
from .support import (
    AcceptanceEvidence,
    AcceptanceFixture,
    load_acceptance_fixture,
    require_gate_settings,
    write_acceptance_evidence,
)


pytestmark = pytest.mark.bulk_rnaseq_real_execution
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def test_controlled_tiny_star_salmon_sortmerna_platform_acceptance(
    tmp_path: Path,
) -> None:
    """Prove platform execution mechanics, not biological validity."""
    settings = require_gate_settings()
    fixture = load_acceptance_fixture(settings.fixture_manifest)

    with PlatformAcceptanceHarness(
        gate_settings=settings,
        repository_root=REPOSITORY_ROOT,
        temporary_root=(tmp_path / "platform").resolve(),
    ) as harness:
        first = harness.execute(fixture)
        second = harness.execute(fixture)

    _assert_required_surface(first, fixture)
    _assert_required_surface(second, fixture)
    for evidence in (first, second):
        _assert_full_trace_contract(
            harness.workspace_root / evidence.values.run_id / "reports/trace.txt"
        )
    first_values = first.values
    second_values = second.values

    assert first_values.tested_head == second_values.tested_head
    assert first_values.workflow_build_digest == second_values.workflow_build_digest
    assert first_values.fixture_acceptance_manifest_sha256 == (
        second_values.fixture_acceptance_manifest_sha256
    )
    assert first_values.fixture_source_manifest_sha256 == (
        second_values.fixture_source_manifest_sha256
    )
    assert first_values.fixture_source_identity_sha256 == (
        second_values.fixture_source_identity_sha256
    )
    assert first_values.fixture_index_provenance_manifest_sha256 == (
        second_values.fixture_index_provenance_manifest_sha256
    )
    assert first_values.fixture_index_provenance_identity_sha256 == (
        second_values.fixture_index_provenance_identity_sha256
    )
    assert first_values.validated_snapshot_id != second_values.validated_snapshot_id
    assert first_values.validated_payload_digest == (
        second_values.validated_payload_digest
    )
    assert first_values.cache_identity_sha256 == second_values.cache_identity_sha256
    assert first_values.input_closure_sha256 == second_values.input_closure_sha256
    assert first_values.ribo_database_closure_sha256 == (
        second_values.ribo_database_closure_sha256
    )
    assert first_values.execution_implementation_manifest_sha256 == (
        second_values.execution_implementation_manifest_sha256
    )
    assert first_values.execution_implementation_aggregate_sha256 == (
        second_values.execution_implementation_aggregate_sha256
    )
    assert first_values.container_process_audit_sha256 == (
        second_values.container_process_audit_sha256
    )
    assert first_values.artifact_output_types == second_values.artifact_output_types
    assert first_values.qc_metric_keys == second_values.qc_metric_keys
    assert first_values.qc_sample_ids == second_values.qc_sample_ids
    assert first_values.artifact_sample_output_types == (
        second_values.artifact_sample_output_types
    )
    assert first_values.qc_sample_metric_keys == second_values.qc_sample_metric_keys
    assert first_values.qc_sample_metric_values == (
        second_values.qc_sample_metric_values
    )

    assert first_values.run_id != second_values.run_id
    assert first_values.job_id != second_values.job_id
    assert first_values.workspace_identity_sha256 != (
        second_values.workspace_identity_sha256
    )
    assert first_values.artifact_attempt_id != second_values.artifact_attempt_id
    assert first_values.artifact_generation != second_values.artifact_generation
    assert first_values.qc_attempt_id != second_values.qc_attempt_id
    assert first_values.qc_generation != second_values.qc_generation
    assert first_values.artifact_attempt_id != first_values.qc_attempt_id
    assert second_values.artifact_attempt_id != second_values.qc_attempt_id

    evidence_root = (tmp_path / "evidence").resolve()
    for index, evidence in enumerate((first, second), start=1):
        destination = evidence_root / f"platform-run-{index}.json"
        write_acceptance_evidence(evidence, destination)
        rendered = destination.read_text(encoding="utf-8")
        payload = json.loads(rendered)
        assert AcceptanceEvidence.from_dict(payload) == evidence
        assert "/" not in rendered
        for private_value in (
            settings.redis_url,
            str(settings.runtime_root),
            str(settings.fixture_manifest),
            str(settings.docker_executable),
            str(settings.docker_socket),
            str(harness.workspace_root),
        ):
            assert private_value not in rendered


def _assert_required_surface(
    evidence: AcceptanceEvidence,
    fixture: AcceptanceFixture,
) -> None:
    values = evidence.values
    assert values.lifecycle_status == "succeeded"
    assert values.artifact_revision == 1
    assert values.qc_revision == 1
    assert values.cleanup_confirmed is True
    assert set(fixture.required_artifact_output_types) <= set(
        values.artifact_output_types
    )
    assert set(fixture.required_qc_metric_keys) <= set(values.qc_metric_keys)
    assert set(fixture.required_sample_ids) <= set(values.qc_sample_ids)
    assert set(fixture.required_artifact_sample_output_types) <= set(
        values.artifact_sample_output_types
    )
    assert set(fixture.required_qc_sample_metric_keys) <= set(
        values.qc_sample_metric_keys
    )
    assert set(fixture.required_qc_sample_metric_values) <= set(
        values.qc_sample_metric_values
    )


def _assert_full_trace_contract(trace_path: Path) -> None:
    content = trace_path.read_bytes()
    if not content or len(content) > 16 * 1024 * 1024:
        raise AssertionError("full acceptance trace is missing or outside its bound")
    lines = content.decode("utf-8", errors="strict").splitlines()
    header = lines[0].split("\t")
    try:
        name_index = header.index("name")
        status_index = header.index("status")
    except ValueError:
        raise AssertionError("full acceptance trace fields are incomplete") from None
    rows = [line.split("\t") for line in lines[1:] if line]
    if not rows or any(len(row) != len(header) for row in rows):
        raise AssertionError("full acceptance trace rows are malformed")
    if any(row[status_index] != "COMPLETED" for row in rows):
        raise AssertionError("full acceptance trace contains an incomplete process")
    processes = tuple(row[name_index].split(" (", 1)[0] for row in rows)
    for required in ("STAR_ALIGN", "SALMON_QUANT", "FASTQC", "SORTMERNA"):
        assert any(name.endswith(f":{required}") for name in processes)
    assert not any(
        name.endswith(
            (
                ":STAR_GENOMEGENERATE",
                ":SALMON_INDEX",
                ":SORTMERNA_INDEX",
                ":MAKE_TRANSCRIPTS_FASTA",
                ":RSEM_PREPAREREFERENCE",
            )
        )
        for name in processes
    )
