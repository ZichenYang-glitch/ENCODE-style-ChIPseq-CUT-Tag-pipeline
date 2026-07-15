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
    BulkRnaSeqWorkflowAdapter,
    RuntimeAssetBinding,
)
from encode_pipeline.adapters.bulk_rnaseq.execution_identity import (
    EXECUTION_IMPLEMENTATION_MANIFEST_FILE,
    EXECUTION_IMPLEMENTATION_PATHS,
    build_execution_implementation_manifest,
    canonical_execution_manifest_bytes,
    verify_execution_implementation,
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


def _assets() -> VerifiedRuntimeAssets:
    root = Path("/runtime")
    return VerifiedRuntimeAssets(
        root=root,
        source_tree=root / "source/rnaseq",
        nextflow_executable=root / "nextflow/nextflow-25.04.3-dist",
        plugin_root=root / "plugins",
        plugin_archive=root / "plugins/nf-schema-2.5.1.zip",
        plugin_meta=root / "plugins/nf-schema-2.5.1-meta.json",
        plugin_tree=root / "plugins/nf-schema-2.5.1",
        container_lock=root / "containers/availability-lock.json",
        containers=(),
        source_tree_sha256="1" * 64,
        runtime_identity_sha256="2" * 64,
        nextflow_sha256="3" * 64,
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


def test_adapter_version_manifest_and_code_identity_each_change_build_digest(
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
        implementation=original.value,
    )
    changed_version_digest = execution_module._runtime_build_digest(
        _assets(),
        adapter_version="1.0.1",
        implementation=original.value,
    )
    changed_code_digest = execution_module._runtime_build_digest(
        _assets(),
        adapter_version="1.0.0",
        implementation=changed.value,
    )

    assert len({original_digest, changed_version_digest, changed_code_digest}) == 3


def test_capture_fails_before_assets_when_implementation_is_invalid(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
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
        "verify_runtime_assets",
        unexpected_asset_verification,
    )
    binding = BulkRnaSeqExecutionBinding(
        assets=RuntimeAssetBinding(root=(tmp_path / "runtime").resolve())
    )

    result = BulkRnaSeqWorkflowAdapter(execution=binding).capture_build_identity()

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_EXECUTION_IMPLEMENTATION_INVALID"
    assert str(tmp_path) not in json.dumps(result.to_dict())
