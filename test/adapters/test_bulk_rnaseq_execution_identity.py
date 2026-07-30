"""Installed implementation identity tests for the bulk RNA-seq runtime."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil

import pytest

import encode_pipeline.adapters.bulk_rnaseq.execution as execution_module
from encode_pipeline.adapters.bulk_rnaseq import (
    BulkRnaSeqExecutionBinding,
    BulkRnaSeqTranscriptomeBinding,
    BulkRnaSeqWorkflowAdapter,
    RuntimeAssetBinding,
)
from encode_pipeline.adapters.bulk_rnaseq.execution_identity import (
    EXECUTION_IMPLEMENTATION_MANIFEST_FILE,
    EXECUTION_IMPLEMENTATION_PATHS,
    EXECUTION_PERSISTENCE_CONTRACT_PATH,
    ExecutionImplementationQualification,
    build_execution_implementation_manifest,
    canonical_execution_manifest_bytes,
    qualify_execution_implementation,
    verify_execution_implementation,
)
from encode_pipeline.adapters.bulk_rnaseq.qualification import (
    DEFAULT_EXECUTION_QUALIFICATION_FILE,
    BulkRnaSeqExecutionMode,
    canonical_qualification_bytes,
    load_default_execution_qualification,
    qualification_document,
)
from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import (
    VerifiedRuntimeAssets,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
MANIFEST_PATH = (
    PROJECT_ROOT
    / "src/encode_pipeline/contracts/nfcore_rnaseq"
    / EXECUTION_IMPLEMENTATION_MANIFEST_FILE
)
QUALIFICATION_PATH = (
    PROJECT_ROOT
    / "src/encode_pipeline/contracts/nfcore_rnaseq"
    / DEFAULT_EXECUTION_QUALIFICATION_FILE
)

PRODUCTION_PERSISTENCE_PATHS = frozenset(
    {
        EXECUTION_PERSISTENCE_CONTRACT_PATH,
        "src/encode_pipeline/persistence/__init__.py",
        "src/encode_pipeline/persistence/runtime.py",
        "src/encode_pipeline/persistence/database.py",
        "src/encode_pipeline/persistence/data_registry.py",
        "src/encode_pipeline/persistence/migrations.py",
        "src/encode_pipeline/persistence/models.py",
        "src/encode_pipeline/persistence/repositories.py",
        "src/encode_pipeline/persistence/alembic/env.py",
        "src/encode_pipeline/persistence/alembic/versions/20260711_01_run_persistence.py",
        "src/encode_pipeline/persistence/alembic/versions/20260711_02_run_execution_assignments.py",
        "src/encode_pipeline/persistence/alembic/versions/20260712_03_run_workflow_build_identities.py",
        "src/encode_pipeline/persistence/alembic/versions/20260712_04_run_cancellation_intent.py",
        "src/encode_pipeline/persistence/alembic/versions/20260712_05_run_qc_metrics.py",
        "src/encode_pipeline/persistence/alembic/versions/20260714_06_validated_input_snapshots.py",
        "src/encode_pipeline/persistence/alembic/versions/20260717_08_run_result_generations.py",
        "src/encode_pipeline/persistence/alembic/versions/20260726_09_project_sample_registry.py",
        "src/encode_pipeline/platform/data_registry.py",
        "src/encode_pipeline/services/data_registry_repositories.py",
    }
)
PRODUCTION_RESULT_EXTRACTION_PATHS = frozenset(
    {
        "src/encode_pipeline/platform/result_generations.py",
        "src/encode_pipeline/services/artifact_extraction.py",
        "src/encode_pipeline/services/run_repositories.py",
        "src/encode_pipeline/services/runs.py",
    }
)
PRODUCT_ONLY_RESULT_SURFACE_PATHS = frozenset(
    {
        "src/encode_pipeline/api/routes/agent.py",
        "src/encode_pipeline/api/routes/artifacts.py",
        "src/encode_pipeline/api/routes/qc_metrics.py",
        "src/encode_pipeline/services/artifact_downloads.py",
    }
)


def _assets() -> VerifiedRuntimeAssets:
    root = Path("/runtime")
    return VerifiedRuntimeAssets(
        root=root,
        source_tree=root / "source/rnaseq",
        nextflow_executable=root / "nextflow/nextflow-25.04.3-dist",
        jdk_archive=root / "jdk/corretto.tar.gz",
        jdk_tree=root / "jdk/corretto",
        java_executable=root / "jdk/corretto/bin/java",
        plugin_root=root / "plugins",
        plugin_archive=root / "plugins/nf-schema-2.5.1.zip",
        plugin_meta=root / "plugins/nf-schema-2.5.1-meta.json",
        plugin_tree=root / "plugins/nf-schema-2.5.1",
        container_lock=root / "containers/availability-lock.json",
        containers=(),
        source_tree_sha256="1" * 64,
        runtime_identity_sha256="2" * 64,
        nextflow_sha256="3" * 64,
        jdk_archive_sha256="8" * 64,
        jdk_tree_sha256="9" * 64,
        java_executable_sha256="a" * 64,
        plugin_archive_sha256="4" * 64,
        plugin_tree_sha256="5" * 64,
        container_inventory_sha256="6" * 64,
        container_lock_sha256="7" * 64,
    )


def _copy_controlled_implementation(destination: Path) -> Path:
    package_root = destination / "src/encode_pipeline"
    for logical_path in EXECUTION_IMPLEMENTATION_PATHS:
        source = PROJECT_ROOT.joinpath(*PurePosixPath(logical_path).parts)
        target = destination.joinpath(*PurePosixPath(logical_path).parts)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    return package_root


def test_committed_execution_implementation_manifest_verifies_installed_files():
    result = verify_execution_implementation()

    assert result.is_success
    assert result.value is not None
    assert len(result.value.files) == len(EXECUTION_IMPLEMENTATION_PATHS)
    assert tuple(item.path for item in result.value.files) == (
        EXECUTION_IMPLEMENTATION_PATHS
    )
    assert (
        result.value.manifest_sha256
        == hashlib.sha256(MANIFEST_PATH.read_bytes()).hexdigest()
    )
    assert len(result.value.aggregate_sha256) == 64


def test_product_only_result_surfaces_are_not_scientific_implementation_files():
    assert PRODUCT_ONLY_RESULT_SURFACE_PATHS.isdisjoint(EXECUTION_IMPLEMENTATION_PATHS)


def test_committed_execution_manifest_is_exact_canonical_build():
    expected = canonical_execution_manifest_bytes(
        build_execution_implementation_manifest(PROJECT_ROOT)
    )

    assert MANIFEST_PATH.read_bytes() == expected


def test_committed_source_qualification_is_canonical_path_free_and_exact():
    implementation = verify_execution_implementation()
    assert implementation.is_success
    document = json.loads(QUALIFICATION_PATH.read_bytes())

    assert QUALIFICATION_PATH.read_bytes() == canonical_qualification_bytes(document)
    assert QUALIFICATION_PATH.name == "default-execution-qualification-1.1.0.json"
    assert document["schema_version"] == "1.1.0"
    result = load_default_execution_qualification(implementation.value)

    assert result.is_success
    assert result.value is not None
    assert result.value.implementation.matches(implementation.value)
    assert (
        result.value.implementation.persistence_contract_version
        == implementation.value.persistence_contract.contract_version
    )
    assert (
        result.value.implementation.persistence_contract_sha256
        == implementation.value.persistence_contract.sha256
    )
    serialized = json.dumps(document)
    persistence = document["execution_implementation"]["persistence_contract"]
    assert persistence["path"] == EXECUTION_PERSISTENCE_CONTRACT_PATH
    assert persistence["version"] == "1.1.0"
    assert set(persistence["schema_projection"]) == {"scheme", "sha256", "tables"}
    assert "required_schema" not in persistence
    assert str(PROJECT_ROOT) not in serialized
    assert "/tmp/" not in serialized
    assert "HELIXWEAVE_BULK_RNASEQ_RUNTIME_ROOT" not in serialized


def test_persistence_contract_change_changes_identity_and_stales_qualification(
    tmp_path: Path,
):
    original = verify_execution_implementation()
    assert original.is_success

    project = tmp_path / "changed-persistence-contract"
    package_root = _copy_controlled_implementation(project)
    contract_path = project.joinpath(
        *PurePosixPath(EXECUTION_PERSISTENCE_CONTRACT_PATH).parts
    )
    contract = json.loads(contract_path.read_bytes())
    contract["contract_version"] = "1.1.1"
    contract_path.write_bytes(
        json.dumps(
            contract,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        + b"\n"
    )

    changed_manifest = build_execution_implementation_manifest(project)
    changed = verify_execution_implementation(
        manifest_bytes=canonical_execution_manifest_bytes(changed_manifest),
        package_root=package_root,
    )

    assert changed.is_success
    assert changed.value.aggregate_sha256 != original.value.aggregate_sha256
    assert changed.value.persistence_contract.contract_version == "1.1.1"
    stale = load_default_execution_qualification(
        changed.value,
        content=QUALIFICATION_PATH.read_bytes(),
    )
    assert stale.is_failure
    assert stale.errors[0].code == "BULK_RNASEQ_EXECUTION_QUALIFICATION_INVALID"


def test_source_qualification_rejects_stale_malformed_or_extended_documents():
    implementation = verify_execution_implementation()
    assert implementation.is_success
    exact = qualification_document(implementation.value)
    stale = deepcopy(exact)
    stale["execution_implementation"]["aggregate_sha256"] = "0" * 64
    extended = deepcopy(exact)
    extended["private_runtime_path"] = str(PROJECT_ROOT)
    corrupted = deepcopy(exact)
    corrupted["record_sha256"] = "2" * 64
    malformed_documents = (
        b'{"schema_version":"1.1.0","schema_version":"1.1.0"}',
        canonical_qualification_bytes(stale),
        canonical_qualification_bytes(extended),
        canonical_qualification_bytes(corrupted),
    )

    for content in malformed_documents:
        result = load_default_execution_qualification(
            implementation.value,
            content=content,
        )
        assert result.is_failure
        assert result.errors[0].code == "BULK_RNASEQ_EXECUTION_QUALIFICATION_INVALID"
        public = json.dumps(result.to_dict())
        assert str(PROJECT_ROOT) not in public
        assert "0" * 64 not in public
        assert "2" * 64 not in public
        assert result.errors[0].technical_message is None
        assert result.errors[0].context == {}


def test_execution_implementation_qualification_requires_every_coordinate():
    implementation = verify_execution_implementation()
    assert implementation.is_success
    exact = ExecutionImplementationQualification.from_verified(implementation.value)

    admitted = qualify_execution_implementation(implementation.value, exact)

    assert admitted.is_success
    assert admitted.value is implementation.value
    mismatches = (
        ExecutionImplementationQualification(
            manifest_sha256="0" * 64,
            aggregate_sha256=exact.aggregate_sha256,
            file_count=exact.file_count,
            persistence_contract_version=exact.persistence_contract_version,
            persistence_contract_sha256=exact.persistence_contract_sha256,
        ),
        ExecutionImplementationQualification(
            manifest_sha256=exact.manifest_sha256,
            aggregate_sha256="0" * 64,
            file_count=exact.file_count,
            persistence_contract_version=exact.persistence_contract_version,
            persistence_contract_sha256=exact.persistence_contract_sha256,
        ),
        ExecutionImplementationQualification(
            manifest_sha256=exact.manifest_sha256,
            aggregate_sha256=exact.aggregate_sha256,
            file_count=exact.file_count + 1,
            persistence_contract_version=exact.persistence_contract_version,
            persistence_contract_sha256=exact.persistence_contract_sha256,
        ),
        ExecutionImplementationQualification(
            manifest_sha256=exact.manifest_sha256,
            aggregate_sha256=exact.aggregate_sha256,
            file_count=exact.file_count,
            persistence_contract_version="9.9.9",
            persistence_contract_sha256=exact.persistence_contract_sha256,
        ),
        ExecutionImplementationQualification(
            manifest_sha256=exact.manifest_sha256,
            aggregate_sha256=exact.aggregate_sha256,
            file_count=exact.file_count,
            persistence_contract_version=exact.persistence_contract_version,
            persistence_contract_sha256="0" * 64,
        ),
    )
    for mismatch in mismatches:
        result = qualify_execution_implementation(implementation.value, mismatch)
        assert result.is_failure
        assert result.errors[0].code == "BULK_RNASEQ_EXECUTION_IMPLEMENTATION_INVALID"
        assert "0" * 64 not in json.dumps(result.to_dict())

    assert qualify_execution_implementation(implementation.value, object()).is_failure


def test_manifest_or_installed_code_mismatch_fails_closed_without_leak(
    tmp_path: Path,
):
    original_bytes = MANIFEST_PATH.read_bytes()
    manifest = json.loads(original_bytes)
    tampered = deepcopy(manifest)
    tampered["files"][0]["sha256"] = "0" * 64
    tampered_result = verify_execution_implementation(
        manifest_bytes=canonical_execution_manifest_bytes(tampered),
    )

    package_root = _copy_controlled_implementation(tmp_path / "project")
    adapter = package_root / "adapters/bulk_rnaseq/adapter.py"
    adapter.write_bytes(adapter.read_bytes() + b"\n# changed after manifest capture\n")
    code_result = verify_execution_implementation(
        manifest_bytes=original_bytes,
        package_root=package_root,
    )

    for result in (tampered_result, code_result):
        assert result.is_failure
        assert result.errors[0].code == ("BULK_RNASEQ_EXECUTION_IMPLEMENTATION_INVALID")
        public = json.dumps(result.to_dict())
        assert str(tmp_path) not in public
        assert "0" * 64 not in public
        assert result.errors[0].technical_message is None
        assert result.errors[0].context == {}


def test_adapter_version_mode_manifest_and_code_each_change_build_digest(
    tmp_path: Path,
):
    original = verify_execution_implementation()
    assert original.is_success
    project = tmp_path / "changed-project"
    package_root = _copy_controlled_implementation(project)
    adapter = package_root / "adapters/bulk_rnaseq/adapter.py"
    adapter.write_bytes(adapter.read_bytes() + b"\n# intentional identity change\n")
    changed_manifest = build_execution_implementation_manifest(project)
    changed = verify_execution_implementation(
        manifest_bytes=canonical_execution_manifest_bytes(changed_manifest),
        package_root=package_root,
    )
    assert changed.is_success
    assert changed.value.aggregate_sha256 != original.value.aggregate_sha256

    original_digest = execution_module._runtime_build_digest(
        _assets(),
        adapter_version="1.0.0",
        adapter_variant="runtime-v1",
        execution_mode=BulkRnaSeqExecutionMode.STANDARD,
        implementation=original.value,
        transcriptome_identity_sha256="f" * 64,
    )
    changed_version_digest = execution_module._runtime_build_digest(
        _assets(),
        adapter_version="1.0.1",
        adapter_variant="runtime-v1",
        execution_mode=BulkRnaSeqExecutionMode.STANDARD,
        implementation=original.value,
        transcriptome_identity_sha256="f" * 64,
    )
    changed_code_digest = execution_module._runtime_build_digest(
        _assets(),
        adapter_version="1.0.0",
        adapter_variant="runtime-v1",
        execution_mode=BulkRnaSeqExecutionMode.STANDARD,
        implementation=changed.value,
        transcriptome_identity_sha256="f" * 64,
    )
    changed_variant_digest = execution_module._runtime_build_digest(
        _assets(),
        adapter_version="1.0.0",
        adapter_variant="results-v1",
        execution_mode=BulkRnaSeqExecutionMode.STANDARD,
        implementation=original.value,
        transcriptome_identity_sha256="f" * 64,
    )
    changed_mode_digest = execution_module._runtime_build_digest(
        _assets(),
        adapter_version="1.0.0",
        adapter_variant="runtime-v1",
        execution_mode=BulkRnaSeqExecutionMode.RAPID_QUANT,
        implementation=original.value,
        transcriptome_identity_sha256="f" * 64,
    )

    assert (
        len(
            {
                original_digest,
                changed_version_digest,
                changed_code_digest,
                changed_variant_digest,
                changed_mode_digest,
            }
        )
        == 5
    )


def test_production_sqlite_replacement_is_bound_and_changes_results_build_identity(
    tmp_path: Path,
):
    assert PRODUCTION_PERSISTENCE_PATHS.issubset(EXECUTION_IMPLEMENTATION_PATHS)
    original_bytes = MANIFEST_PATH.read_bytes()
    original = verify_execution_implementation(manifest_bytes=original_bytes)
    assert original.is_success

    project = tmp_path / "changed-production-persistence"
    package_root = _copy_controlled_implementation(project)
    repository = package_root / "persistence/repositories.py"
    repository.write_bytes(
        repository.read_bytes() + b"\n# intentional QC invalidation identity change\n"
    )

    stale = verify_execution_implementation(
        manifest_bytes=original_bytes,
        package_root=package_root,
    )
    assert stale.is_failure
    assert stale.errors[0].code == "BULK_RNASEQ_EXECUTION_IMPLEMENTATION_INVALID"

    changed_manifest = build_execution_implementation_manifest(project)
    changed = verify_execution_implementation(
        manifest_bytes=canonical_execution_manifest_bytes(changed_manifest),
        package_root=package_root,
    )
    assert changed.is_success
    assert changed.value.aggregate_sha256 != original.value.aggregate_sha256
    original_digest = execution_module._runtime_build_digest(
        _assets(),
        adapter_version="1.0.0",
        adapter_variant="results-v1",
        execution_mode=BulkRnaSeqExecutionMode.STANDARD,
        implementation=original.value,
        transcriptome_identity_sha256="f" * 64,
    )
    changed_digest = execution_module._runtime_build_digest(
        _assets(),
        adapter_version="1.0.0",
        adapter_variant="results-v1",
        execution_mode=BulkRnaSeqExecutionMode.STANDARD,
        implementation=changed.value,
        transcriptome_identity_sha256="f" * 64,
    )
    assert changed_digest != original_digest


def test_production_artifact_extraction_is_bound_and_stale_manifest_fails_closed(
    tmp_path: Path,
):
    assert PRODUCTION_RESULT_EXTRACTION_PATHS.issubset(EXECUTION_IMPLEMENTATION_PATHS)
    original_bytes = MANIFEST_PATH.read_bytes()
    original = verify_execution_implementation(manifest_bytes=original_bytes)
    assert original.is_success

    project = tmp_path / "changed-artifact-extraction"
    package_root = _copy_controlled_implementation(project)
    extraction_service = package_root / "services/artifact_extraction.py"
    extraction_service.write_bytes(
        extraction_service.read_bytes()
        + b"\n# intentional artifact extraction identity change\n"
    )

    stale = verify_execution_implementation(
        manifest_bytes=original_bytes,
        package_root=package_root,
    )
    assert stale.is_failure
    assert stale.errors[0].code == "BULK_RNASEQ_EXECUTION_IMPLEMENTATION_INVALID"

    changed_manifest = build_execution_implementation_manifest(project)
    changed = verify_execution_implementation(
        manifest_bytes=canonical_execution_manifest_bytes(changed_manifest),
        package_root=package_root,
    )
    assert changed.is_success
    assert changed.value.aggregate_sha256 != original.value.aggregate_sha256


def test_unrelated_migration_revision_does_not_change_or_reject_identity(
    tmp_path: Path,
):
    original_bytes = MANIFEST_PATH.read_bytes()
    original = verify_execution_implementation(manifest_bytes=original_bytes)
    assert original.is_success

    project = tmp_path / "unexpected-migration"
    package_root = _copy_controlled_implementation(project)
    extra = package_root / "persistence/alembic/versions/20990101_99_unlisted.py"
    extra.parent.mkdir(parents=True, exist_ok=True)
    extra.write_text(
        "from alembic import op\n"
        "import sqlalchemy as sa\n\n"
        "revision = '20990101_99'\n"
        "down_revision = '20260726_10'\n"
        "branch_labels = None\n"
        "depends_on = None\n\n"
        "def upgrade():\n"
        "    op.add_column(\n"
        "        'runs',\n"
        "        sa.Column('review_only_note', sa.Text(), nullable=True),\n"
        "    )\n\n"
        "def downgrade():\n"
        "    op.drop_column('runs', 'review_only_note')\n",
        encoding="utf-8",
    )

    installed = verify_execution_implementation(
        manifest_bytes=original_bytes,
        package_root=package_root,
    )
    regenerated = canonical_execution_manifest_bytes(
        build_execution_implementation_manifest(project)
    )

    assert installed.is_success
    assert installed.value.aggregate_sha256 == original.value.aggregate_sha256
    assert regenerated == original_bytes


def test_mock_agent_route_change_does_not_change_scientific_aggregate(
    tmp_path: Path,
):
    original_bytes = MANIFEST_PATH.read_bytes()
    original = verify_execution_implementation(manifest_bytes=original_bytes)
    assert original.is_success

    project = tmp_path / "product-only-route"
    package_root = _copy_controlled_implementation(project)
    route = package_root / "api/routes/agent.py"
    route.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(PROJECT_ROOT / "src/encode_pipeline/api/routes/agent.py", route)
    route.write_bytes(route.read_bytes() + b"\n# product-only mock route change\n")

    changed_manifest = build_execution_implementation_manifest(project)
    changed = verify_execution_implementation(
        manifest_bytes=canonical_execution_manifest_bytes(changed_manifest),
        package_root=package_root,
    )

    assert changed.is_success
    assert changed.value.aggregate_sha256 == original.value.aggregate_sha256
    assert canonical_execution_manifest_bytes(changed_manifest) == original_bytes


def test_missing_listed_production_migration_revision_fails_closed(
    tmp_path: Path,
):
    project = tmp_path / "missing-migration"
    package_root = _copy_controlled_implementation(project)
    missing = (
        package_root
        / "persistence/alembic/versions/20260726_09_project_sample_registry.py"
    )
    missing.unlink()

    installed = verify_execution_implementation(
        manifest_bytes=MANIFEST_PATH.read_bytes(),
        package_root=package_root,
    )

    assert installed.is_failure
    assert installed.errors[0].code == "BULK_RNASEQ_EXECUTION_IMPLEMENTATION_INVALID"
    with pytest.raises(ValueError, match="controlled implementation"):
        build_execution_implementation_manifest(project)


def test_capture_fails_before_assets_when_implementation_is_invalid(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    current = verify_execution_implementation()
    assert current.is_success
    qualification = ExecutionImplementationQualification.from_verified(current.value)
    invalid = verify_execution_implementation(manifest_bytes=b"{}")
    assert invalid.is_failure
    monkeypatch.setattr(
        execution_module,
        "verify_execution_implementation",
        lambda: invalid,
    )

    def unexpected_asset_verification(_binding):
        raise AssertionError("assets must not be read after implementation mismatch")

    monkeypatch.setattr(
        execution_module,
        "_acquire_runtime_assets",
        unexpected_asset_verification,
    )
    binding = BulkRnaSeqExecutionBinding(
        assets=RuntimeAssetBinding(root=(tmp_path / "runtime").resolve()),
        transcriptome=BulkRnaSeqTranscriptomeBinding(
            reference_id="tiny-ref",
            fasta_sha256="a" * 64,
            gtf_sha256="b" * 64,
            transcript_fasta=(tmp_path / "transcripts.fa").resolve(),
            transcript_fasta_sha256="c" * 64,
        ),
        implementation_qualification=qualification,
    )

    result = BulkRnaSeqWorkflowAdapter(execution=binding).capture_build_identity()

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_EXECUTION_IMPLEMENTATION_INVALID"
    assert str(tmp_path) not in json.dumps(result.to_dict())
