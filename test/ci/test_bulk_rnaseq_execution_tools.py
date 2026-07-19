"""Behavioral coverage for the committed bulk RNA-seq contract tools."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

import scripts.audit_bulk_rnaseq_container_processes as container_audit
import scripts.generate_bulk_rnaseq_execution_manifest as manifest_generator
from encode_pipeline.adapters.bulk_rnaseq.upstream import (
    PLATFORM_OWNED_NATIVE_PARAMETERS,
)


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _canonical_sha256(value: object) -> str:
    return _sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    )


@pytest.fixture
def container_audit_case(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> tuple[Path, Path, Path]:
    source_root = tmp_path / "source"
    files = {
        "workflow.nf": (
            b"include { INCLUDED as PARABRICKS_STARGENOMEGENERATE } "
            b"from './modules/included'\n"
        ),
        "modules/included/main.nf": (
            b'process INCLUDED {\n  container "acme/included:1"\n'
            b"  input:\n  path reads\n}\n"
        ),
        "modules/excluded/main.nf": (
            b'process EXCLUDED {\n  container "quay.io/acme/excluded:2"\n'
            b"  input:\n  path reads\n}\n"
        ),
        "conf/modules/prepare_genome.config": (
            b"process {\n withName: '.*PARABRICKS_STARGENOMEGENERATE' {\n"
            b"  container = 'quay.io/acme/default:1'\n }\n}\n"
        ),
        "conf/arm.config": (
            b"process {\n withName: 'INCLUDED' {\n"
            b"  container = 'quay.io/acme/arm:1'\n }\n}\n"
        ),
        "modules/local/star_genomeparams_upgrade/tests/nextflow.legacy_index.config": (
            b"process {\n withName: 'INCLUDED' {\n"
            b"  container = 'quay.io/acme/test-a:1'\n }\n}\n"
        ),
        "modules/nf-core/parabricks/rnafq2bam/tests/nextflow.config": (
            b"process {\n withName: 'INCLUDED' {\n"
            b"  container = 'quay.io/acme/test-b:1'\n }\n}\n"
        ),
        "subworkflows/local/align_star/tests/nextflow.upgraded_legacy.config": (
            b"process {\n withName: 'INCLUDED' {\n"
            b"  container = 'quay.io/acme/test-c:1'\n }\n}\n"
        ),
    }
    for relative, content in files.items():
        path = source_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)

    manifest_entries = [
        {
            "path": relative,
            "size_bytes": len(content),
            "sha256": _sha256(content),
            "executable": False,
        }
        for relative, content in sorted(files.items())
    ]
    source_manifest = tmp_path / "source-manifest.json"
    source_manifest.write_text(
        json.dumps({"files": manifest_entries}),
        encoding="utf-8",
    )
    included_path = "modules/included/main.nf"
    inventory = tmp_path / "container-inventory.json"
    inventory.write_text(
        json.dumps(
            {
                "entries": [
                    {
                        "process": "INCLUDED",
                        "source_file": included_path,
                        "source_file_sha256": _sha256(files[included_path]),
                        "image_coordinate": "quay.io/acme/included:1",
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    committed_audit = tmp_path / "container-process-audit.json"
    replacements = {
        "SOURCE_MANIFEST": source_manifest,
        "CONTAINER_INVENTORY": inventory,
        "COMMITTED_AUDIT": committed_audit,
        "SOURCE_TREE_SHA256": _canonical_sha256(manifest_entries),
        "PROCESS_UNIVERSE_COUNT": 2,
        "INCLUDED_PROCESS_COUNT": 1,
        "EXCLUDED_PROCESS_COUNT": 1,
        "CONFIG_CONTAINER_ASSIGNMENT_COUNT": 5,
        "DEFAULT_CONFIG_CONTAINER_ASSIGNMENT_COUNT": 1,
        "EXCLUSIONS": {"EXCLUDED": "unsupported_test_route"},
    }
    for name, value in replacements.items():
        monkeypatch.setattr(container_audit, name, value)
    return source_root, inventory, committed_audit


def test_container_audit_normalizes_and_rejects_ambiguous_image_coordinates():
    assert (
        container_audit._normalise_image_coordinate("biocontainers/fastqc:1.0")
        == "quay.io/biocontainers/fastqc:1.0"
    )
    assert (
        container_audit._normalise_image_coordinate("quay.io/example/tool:2.0")
        == "quay.io/example/tool:2.0"
    )

    for invalid in (
        "quay.io/example/tool@sha256:abc",
        "quay.io/example/tool:latest value",
        "quay.io/example/tool",
        "localhost:5000",
    ):
        with pytest.raises(container_audit.AuditError):
            container_audit._normalise_image_coordinate(invalid)


def test_rsem_transcript_preparation_exclusion_requires_platform_owned_input():
    assert "transcript_fasta" in PLATFORM_OWNED_NATIVE_PARAMETERS
    assert container_audit.EXCLUSIONS["RSEM_PREPAREREFERENCE"] == (
        "platform_owned_transcript_fasta"
    )
    assert container_audit.EXCLUSIONS["RSEM_CALCULATEEXPRESSION"] == (
        "unsupported_rsem_route"
    )


def test_container_audit_extracts_fixed_and_conditional_docker_coordinates():
    fixed = """
process FASTQC {
    container "quay.io/biocontainers/fastqc:0.12.1"
    input:
    path reads
}
"""
    conditional = """
process TOOL {
    container "${ true ? 'biocontainers/tool:1.0' : 'quay.io/example/tool:2.0' }"
    input:
    path reads
}
"""

    assert container_audit._declared_image_coordinates(fixed) == (
        "quay.io/biocontainers/fastqc:0.12.1",
    )
    assert container_audit._declared_image_coordinates(conditional) == (
        "quay.io/biocontainers/tool:1.0",
        "quay.io/example/tool:2.0",
    )
    with pytest.raises(container_audit.AuditError, match="no container directive"):
        container_audit._declared_image_coordinates("process TOOL { input: }")
    with pytest.raises(container_audit.AuditError, match="no Docker image coordinate"):
        container_audit._declared_image_coordinates(
            'process TOOL {\n container "${ params.container }"\n input:\n path reads\n}'
        )


def test_container_audit_loads_only_json_objects_and_hashes_canonically(tmp_path):
    object_path = tmp_path / "object.json"
    object_path.write_text('{"b": 2, "a": 1}\n', encoding="utf-8")
    list_path = tmp_path / "list.json"
    list_path.write_text("[]\n", encoding="utf-8")

    assert container_audit._load_object(object_path) == {"a": 1, "b": 2}
    with pytest.raises(container_audit.AuditError, match="JSON object"):
        container_audit._load_object(list_path)
    assert container_audit._canonical_sha256({"b": 2, "a": 1}) == (
        container_audit._canonical_sha256({"a": 1, "b": 2})
    )


def test_container_audit_builds_the_complete_deterministic_universe(
    container_audit_case,
):
    source_root, _, _ = container_audit_case

    first = container_audit.build_audit(source_root)

    assert first == container_audit.build_audit(source_root)
    assert first["process_universe_count"] == 2
    assert [entry["disposition"] for entry in first["processes"]] == [
        "excluded",
        "included",
    ]
    assert first["processes"][0]["exclusion_reason"] == "unsupported_test_route"
    assert first["processes"][1]["image_coordinates"] == ["quay.io/acme/included:1"]
    [override] = first["default_config_container_assignments"]
    assert override["platform_override"] == "deny"
    assert override["matched_aliases"][0]["alias"] == ("PARABRICKS_STARGENOMEGENERATE")


def test_container_audit_fails_closed_for_source_or_inventory_drift(
    container_audit_case,
    capsys,
):
    source_root, inventory, _ = container_audit_case
    included_module = source_root / "modules/included/main.nf"
    original_module = included_module.read_bytes()
    included_module.write_text("tampered\n", encoding="utf-8")

    assert container_audit.main([str(source_root)]) == 1
    assert capsys.readouterr().err == (
        "audit failed: source file identity differs from the pinned manifest\n"
    )

    included_module.write_bytes(original_module)
    inventory_payload = json.loads(inventory.read_bytes())
    inventory_payload["entries"][0]["image_coordinate"] = "quay.io/acme/wrong:9"
    inventory.write_text(json.dumps(inventory_payload), encoding="utf-8")
    with pytest.raises(
        container_audit.AuditError, match="inventory differs from source"
    ):
        container_audit.build_audit(source_root)


def test_container_audit_check_accepts_exact_bytes_and_rejects_drift(
    container_audit_case,
    capsys,
):
    source_root, _, committed_audit = container_audit_case
    content = container_audit._canonical_bytes(container_audit.build_audit(source_root))
    value = json.loads(content)
    committed_audit.write_bytes(content)

    assert container_audit.main([str(source_root), "--check"]) == 0
    accepted = capsys.readouterr()
    assert accepted.err == ""
    assert accepted.out.splitlines() == [
        "processes=2",
        "config_container_assignments=5",
        "config_container_assignments_sha256="
        + value["config_container_assignments_sha256"],
        "default_config_container_assignments=1",
        "default_config_container_assignments_sha256="
        + value["default_config_container_assignments_sha256"],
        "audit_sha256=" + _sha256(content),
    ]

    committed_audit.write_text("{}\n", encoding="utf-8")
    assert container_audit.main([str(source_root), "--check"]) == 1
    rejected = capsys.readouterr()
    assert rejected.out == ""
    assert rejected.err == "committed container process audit differs\n"


def test_execution_manifest_generator_writes_exact_canonical_bytes(
    tmp_path,
    monkeypatch,
    capsys,
):
    destination = (
        tmp_path
        / "src/encode_pipeline/contracts/nfcore_rnaseq"
        / "execution-manifest.json"
    )
    destination.parent.mkdir(parents=True)
    manifest = {"file_count": 2, "aggregate_sha256": "c" * 64}
    content = b'{"aggregate_sha256":"value","file_count":2}\n'

    monkeypatch.setattr(manifest_generator, "PROJECT_ROOT", tmp_path)
    monkeypatch.setattr(
        manifest_generator,
        "EXECUTION_IMPLEMENTATION_MANIFEST_FILE",
        destination.name,
    )
    monkeypatch.setattr(
        manifest_generator,
        "build_execution_implementation_manifest",
        lambda project_root: manifest if project_root == tmp_path else None,
    )
    monkeypatch.setattr(
        manifest_generator,
        "canonical_execution_manifest_bytes",
        lambda value: content if value is manifest else b"",
    )

    assert manifest_generator.main() == 0

    assert destination.read_bytes() == content
    output = capsys.readouterr().out.splitlines()
    assert output == [
        "files=2",
        f"aggregate_sha256={'c' * 64}",
        f"manifest_sha256={hashlib.sha256(content).hexdigest()}",
    ]
