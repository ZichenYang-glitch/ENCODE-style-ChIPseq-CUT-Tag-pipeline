"""Execution-preparation boundaries for supported rRNA removal routes."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

import encode_pipeline.adapters.bulk_rnaseq.execution as execution_module
from encode_pipeline.adapters.bulk_rnaseq import (
    BulkRnaSeqExecutionBinding,
    BulkRnaSeqWorkflowAdapter,
    RuntimeAssetBinding,
)
from encode_pipeline.adapters.bulk_rnaseq.resource_closure import (
    SORTMERNA_INDEX_BINDING_FILENAME,
    SORTMERNA_INDEX_BINDING_SCHEMA_VERSION,
    SORTMERNA_NO_PREBUILT_INDEX_STRATEGY,
    verify_ribo_database_manifest,
)
from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import (
    VerifiedContainerAsset,
    VerifiedRuntimeAssets,
)
from encode_pipeline.platform.adapters import WorkflowInputs, WorkspacePlan
from encode_pipeline.platform.results import Result


def _sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _write(path: Path, value: bytes) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(value)
    return _sha256(value)


def _ribo_database(
    root: Path,
    *,
    first_sequence: bytes = b"ACGT",
    multiple: bool = True,
) -> tuple[dict[str, str], tuple[Path, ...], str]:
    files = [root / "database/z-ribosomal.fa"]
    if multiple:
        files.append(root / "database/a-ribosomal.fa")
    _write(files[0], b">z-ribosomal\n" + first_sequence + b"\n")
    if multiple:
        _write(files[1], b">a-ribosomal\nTGCA\n")
    manifest = root / "database-manifest.txt"
    content = b"".join(f"database/{path.name}\n".encode() for path in files)
    manifest_sha256 = _write(manifest, content)
    result = verify_ribo_database_manifest(
        manifest,
        expected_manifest_sha256=manifest_sha256,
    )
    assert result.is_success
    return (
        {"path": str(manifest), "identity_sha256": manifest_sha256},
        tuple(files),
        result.value.identity_sha256,
    )


def _sortmerna_index(root: Path, *, database_closure_sha256: str) -> dict[str, str]:
    payload = root / "index/ref.idx"
    payload_sha256 = _write(payload, b"tiny-sortmerna-index\n")
    binding = {
        "schema_version": SORTMERNA_INDEX_BINDING_SCHEMA_VERSION,
        "database_closure_sha256": database_closure_sha256,
        "sortmerna_version": "4.3.7",
        "files": [
            {
                "path": "index/ref.idx",
                "size_bytes": payload.stat().st_size,
                "sha256": payload_sha256,
            }
        ],
    }
    manifest = root / SORTMERNA_INDEX_BINDING_FILENAME
    manifest_sha256 = _write(
        manifest,
        (json.dumps(binding, sort_keys=True, separators=(",", ":")) + "\n").encode(),
    )
    return {"path": str(root), "identity_sha256": manifest_sha256}


def _ribo_config(
    database_manifest: dict[str, str],
    *,
    tool: str,
    save_filtered_reads: bool,
    sortmerna_index: dict[str, str] | None = None,
) -> dict[str, object]:
    value: dict[str, object] = {
        "enabled": True,
        "tool": tool,
        "save_filtered_reads": save_filtered_reads,
        "database_manifest": database_manifest,
    }
    if sortmerna_index is not None:
        value["sortmerna_index"] = sortmerna_index
    return value


def _inputs(
    root: Path,
    *,
    layout: str,
    ribosomal_rna_removal: dict[str, object],
) -> WorkflowInputs:
    fasta = root / "inputs/reference.fa"
    gtf = root / "inputs/reference.gtf"
    reference = {
        "reference_id": "tiny-reference",
        "fasta": str(fasta),
        "fasta_sha256": _write(fasta, b">chr1\nACGT\n"),
        "gtf": str(gtf),
        "gtf_sha256": _write(
            gtf,
            b'chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n',
        ),
    }
    fastq_1 = root / "inputs/sample.R1.fastq.gz"
    _write(fastq_1, b"tiny-r1")
    sample = {
        "sample": "S1",
        "library": "lib1",
        "lane": "L001",
        "layout": layout,
        "fastq_1": str(fastq_1),
        "strandedness": "auto",
        "platform": "ILLUMINA",
    }
    if layout == "PE":
        fastq_2 = root / "inputs/sample.R2.fastq.gz"
        _write(fastq_2, b"tiny-r2")
        sample["fastq_2"] = str(fastq_2)
    return WorkflowInputs(
        config={
            "standard": {
                "reference": reference,
                "ribosomal_rna_removal": ribosomal_rna_removal,
            }
        },
        samples=[sample],
        options={},
    )


@pytest.fixture
def execution_binding(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    root = (tmp_path / "runtime").resolve()
    binding = BulkRnaSeqExecutionBinding(assets=RuntimeAssetBinding(root=root))
    verified = VerifiedRuntimeAssets(
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
        containers=(
            VerifiedContainerAsset(
                process="SORTMERNA",
                image="sha256:" + "8" * 64,
                oci_digest="sha256:" + "8" * 64,
                local_asset=root / "containers/assets/sortmerna.tar",
                local_sha256="9" * 64,
                size_bytes=1,
                distribution_manifest=(
                    root / "containers/assets/sortmerna.manifest.json"
                ),
                distribution_manifest_sha256="8" * 64,
                config_digest="sha256:" + "8" * 64,
                runtime_image="sha256:" + "8" * 64,
                rootfs_diff_ids=("sha256:" + "a" * 64,),
            ),
        ),
        source_tree_sha256="1" * 64,
        runtime_identity_sha256="2" * 64,
        nextflow_sha256="3" * 64,
        jdk_archive_sha256="a" * 64,
        jdk_tree_sha256="b" * 64,
        java_executable_sha256="c" * 64,
        plugin_archive_sha256="4" * 64,
        plugin_tree_sha256="5" * 64,
        container_inventory_sha256="6" * 64,
        container_lock_sha256="7" * 64,
    )
    monkeypatch.setattr(
        execution_module,
        "_acquire_runtime_assets",
        lambda _binding: Result.success(verified),
    )
    return binding


def _plan(
    inputs: WorkflowInputs,
    workspace: Path,
    *,
    binding: BulkRnaSeqExecutionBinding,
):
    return BulkRnaSeqWorkflowAdapter(execution=binding).plan_workspace(
        inputs,
        workspace,
    )


def _planned_bytes(plan: WorkspacePlan, relative_path: str) -> bytes:
    return dict(plan.files)[relative_path]


def test_bowtie2_se_can_save_filtered_reads_with_canonical_deterministic_plan(
    tmp_path: Path,
    execution_binding: BulkRnaSeqExecutionBinding,
):
    manifest, database_files, _ = _ribo_database(tmp_path / "rrna")
    inputs = _inputs(
        tmp_path,
        layout="SE",
        ribosomal_rna_removal=_ribo_config(
            manifest,
            tool="bowtie2",
            save_filtered_reads=True,
        ),
    )
    workspace = (tmp_path / "workspace").resolve()

    first = _plan(inputs, workspace, binding=execution_binding)
    second = _plan(inputs, workspace, binding=execution_binding)

    assert first.is_success and second.is_success
    assert first.value == second.value
    params_bytes = _planned_bytes(first.value, "config/params.json")
    assert params_bytes == _planned_bytes(second.value, "config/params.json")
    params = json.loads(params_bytes)
    assert params["remove_ribo_rna"] is True
    assert params["ribo_removal_tool"] == "bowtie2"
    assert params["save_non_ribo_reads"] is True
    assert params["ribo_database_manifest"] == str(
        workspace / "config/ribo-database-manifest.txt"
    )
    assert "sortmerna_index" not in params
    assert (
        params_bytes
        == (json.dumps(params, sort_keys=True, separators=(",", ":")) + "\n").encode()
    )
    assert (
        _planned_bytes(first.value, "config/ribo-database-manifest.txt")
        == ("\n".join(str(path) for path in database_files) + "\n").encode()
    )


def test_bowtie2_pe_without_saved_reads_preserves_upstream_filtering_boundary(
    tmp_path: Path,
    execution_binding: BulkRnaSeqExecutionBinding,
):
    manifest, _, _ = _ribo_database(tmp_path / "rrna")
    inputs = _inputs(
        tmp_path,
        layout="PE",
        ribosomal_rna_removal=_ribo_config(
            manifest,
            tool="bowtie2",
            save_filtered_reads=False,
        ),
    )

    result = _plan(
        inputs,
        (tmp_path / "workspace").resolve(),
        binding=execution_binding,
    )

    assert result.is_success
    params = json.loads(_planned_bytes(result.value, "config/params.json"))
    assert params["remove_ribo_rna"] is True
    assert params["ribo_removal_tool"] == "bowtie2"
    assert params["save_non_ribo_reads"] is False


def test_bowtie2_pe_saved_reads_fail_closed_before_execution_preparation(
    tmp_path: Path,
    execution_binding: BulkRnaSeqExecutionBinding,
):
    manifest, _, _ = _ribo_database(tmp_path / "rrna")
    inputs = _inputs(
        tmp_path,
        layout="PE",
        ribosomal_rna_removal=_ribo_config(
            manifest,
            tool="bowtie2",
            save_filtered_reads=True,
        ),
    )

    result = _plan(
        inputs,
        (tmp_path / "workspace").resolve(),
        binding=execution_binding,
    )

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_RRNA_CONFLICT"
    assert result.errors[0].path == (
        "config.standard.ribosomal_rna_removal.save_filtered_reads"
    )


@pytest.mark.parametrize("tool", ["sortmerna", "bowtie2"])
@pytest.mark.parametrize(
    "hostile_name",
    ["$(touch-PWN).fa", "`touch-PWN`.fa", "rrna;touch-PWN.fa", "rrna db.fa"],
)
def test_rrna_routes_reject_shell_active_database_filenames(
    tmp_path: Path,
    execution_binding: BulkRnaSeqExecutionBinding,
    tool: str,
    hostile_name: str,
):
    root = tmp_path / "rrna"
    database = root / "database" / hostile_name
    _write(database, b">rrna\nACGT\n")
    manifest_path = root / "database-manifest.txt"
    manifest_sha256 = _write(
        manifest_path,
        f"database/{hostile_name}\n".encode(),
    )
    inputs = _inputs(
        tmp_path,
        layout="SE",
        ribosomal_rna_removal=_ribo_config(
            {"path": str(manifest_path), "identity_sha256": manifest_sha256},
            tool=tool,
            save_filtered_reads=False,
        ),
    )

    result = _plan(
        inputs,
        (tmp_path / "workspace").resolve(),
        binding=execution_binding,
    )

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_RIBO_DATABASE_INVALID"
    assert hostile_name not in json.dumps(result.to_dict())


def test_sortmerna_accepts_index_bound_to_exact_database_closure(
    tmp_path: Path,
    execution_binding: BulkRnaSeqExecutionBinding,
):
    manifest, _, database_closure_sha256 = _ribo_database(tmp_path / "rrna")
    index = _sortmerna_index(
        tmp_path / "sortmerna-index",
        database_closure_sha256=database_closure_sha256,
    )
    inputs = _inputs(
        tmp_path,
        layout="PE",
        ribosomal_rna_removal=_ribo_config(
            manifest,
            tool="sortmerna",
            save_filtered_reads=False,
            sortmerna_index=index,
        ),
    )

    result = _plan(
        inputs,
        (tmp_path / "workspace").resolve(),
        binding=execution_binding,
    )

    assert result.is_success
    params = json.loads(_planned_bytes(result.value, "config/params.json"))
    assert params["ribo_removal_tool"] == "sortmerna"
    assert params["sortmerna_index"] == index["path"]
    assert params["ribo_database_manifest"].endswith(
        "/config/ribo-database-manifest.txt"
    )
    config = _planned_bytes(result.value, "config/platform.nextflow.config").decode()
    assert ".*:PREPARE_GENOME:SORTMERNA_INDEX" not in config
    assert ".*:FASTQ_REMOVE_RRNA:SORTMERNA_INDEX" not in config
    execution_identity = json.loads(
        _planned_bytes(result.value, "config/execution-identity.json")
    )
    cache_identity = json.loads(
        _planned_bytes(result.value, "engine/cache-identity.json")
    )
    assert execution_identity["sortmerna_index_build_strategy"] is None
    assert cache_identity["sortmerna_index_build_strategy"] is None


def test_sortmerna_single_database_without_prebuilt_index_has_deterministic_route(
    tmp_path: Path,
    execution_binding: BulkRnaSeqExecutionBinding,
):
    manifest, _, database_closure_sha256 = _ribo_database(
        tmp_path / "rrna",
        multiple=False,
    )
    inputs = _inputs(
        tmp_path,
        layout="SE",
        ribosomal_rna_removal=_ribo_config(
            manifest,
            tool="sortmerna",
            save_filtered_reads=False,
        ),
    )

    workspace = (tmp_path / "workspace").resolve()
    first = _plan(
        inputs,
        workspace,
        binding=execution_binding,
    )
    second = _plan(inputs, workspace, binding=execution_binding)

    assert first.is_success and second.is_success
    assert first.value == second.value
    params = json.loads(_planned_bytes(first.value, "config/params.json"))
    assert "sortmerna_index" not in params
    config = _planned_bytes(first.value, "config/platform.nextflow.config").decode()
    assert config.count(".*:PREPARE_GENOME:SORTMERNA_INDEX") == 1
    assert config.count(".*:FASTQ_REMOVE_RRNA:SORTMERNA_INDEX") == 1
    assert "ext.when = false" in config
    assert "container = 'sha256:" + "8" * 64 + "'" in config
    assert "ext.args = '--index 1'" in config
    assert "cpus = 1" in config
    assert "maxRetries = 0" in config
    assert "errorStrategy = 'terminate'" in config
    execution_identity = json.loads(
        _planned_bytes(first.value, "config/execution-identity.json")
    )
    cache_identity = json.loads(
        _planned_bytes(first.value, "engine/cache-identity.json")
    )
    assert (
        execution_identity["sortmerna_index_build_strategy"]
        == SORTMERNA_NO_PREBUILT_INDEX_STRATEGY
    )
    assert (
        cache_identity["sortmerna_index_build_strategy"]
        == SORTMERNA_NO_PREBUILT_INDEX_STRATEGY
    )
    assert (
        execution_identity["ribo_database_closure_sha256"]
        == database_closure_sha256
        == cache_identity["ribo_database_closure_sha256"]
    )


def test_sortmerna_multiple_databases_without_prebuilt_index_fail_closed(
    tmp_path: Path,
    execution_binding: BulkRnaSeqExecutionBinding,
):
    manifest, _, _ = _ribo_database(tmp_path / "rrna", multiple=True)
    inputs = _inputs(
        tmp_path,
        layout="SE",
        ribosomal_rna_removal=_ribo_config(
            manifest,
            tool="sortmerna",
            save_filtered_reads=False,
        ),
    )

    result = _plan(
        inputs,
        (tmp_path / "workspace").resolve(),
        binding=execution_binding,
    )

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_SORTMERNA_INDEX_BUILD_UNQUALIFIED"


def test_sortmerna_rejects_index_bound_to_different_database_closure(
    tmp_path: Path,
    execution_binding: BulkRnaSeqExecutionBinding,
):
    selected_manifest, _, _ = _ribo_database(tmp_path / "selected-rrna")
    _, _, other_database_closure = _ribo_database(
        tmp_path / "other-rrna",
        first_sequence=b"CCCC",
    )
    index = _sortmerna_index(
        tmp_path / "sortmerna-index",
        database_closure_sha256=other_database_closure,
    )
    inputs = _inputs(
        tmp_path,
        layout="SE",
        ribosomal_rna_removal=_ribo_config(
            selected_manifest,
            tool="sortmerna",
            save_filtered_reads=False,
            sortmerna_index=index,
        ),
    )

    result = _plan(
        inputs,
        (tmp_path / "workspace").resolve(),
        binding=execution_binding,
    )

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_SORTMERNA_INDEX_INVALID"
    assert str(tmp_path) not in json.dumps(result.to_dict())
