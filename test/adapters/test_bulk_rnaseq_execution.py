"""Workspace and command tests for the offline bulk RNA-seq runtime binding."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

import encode_pipeline.adapters.bulk_rnaseq.execution as execution_module
from encode_pipeline.adapters.bulk_rnaseq import (
    BulkRnaSeqExecutionBinding,
    BulkRnaSeqResultsWorkflowAdapter,
    BulkRnaSeqWorkflowAdapter,
    RuntimeAssetBinding,
)
from encode_pipeline.adapters.bulk_rnaseq.reference_closure import (
    REFERENCE_INDEX_MANIFEST,
)
from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import (
    RuntimeAssetAdmission,
    VerifiedContainerAsset,
    VerifiedRuntimeAssets,
)
from encode_pipeline.platform.adapters import (
    QcSourceArtifact,
    QcSourceDocument,
    WorkflowInputs,
    WorkspacePlan,
)
from encode_pipeline.platform.managed_containers import (
    managed_container_endpoint_identity,
    managed_container_scope,
)
from encode_pipeline.platform.results import Result
from encode_pipeline.services.process_runner import _subprocess_environment
from encode_pipeline.testing.adapter_conformance import (
    AdapterConformanceCase,
    verify_adapter_conformance,
)


def _sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _write(path: Path, value: bytes) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(value)
    return _sha256(value)


def _inputs(
    tmp_path: Path,
    *,
    layouts: tuple[str, ...] = ("PE",),
    standard_updates: dict[str, object] | None = None,
) -> WorkflowInputs:
    fasta = tmp_path / "inputs/reference.fa"
    gtf = tmp_path / "inputs/reference.gtf"
    reference = {
        "reference_id": "tiny-ref",
        "fasta": str(fasta),
        "fasta_sha256": _write(fasta, b">chr1\nACGT\n"),
        "gtf": str(gtf),
        "gtf_sha256": _write(gtf, b'chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n'),
    }
    standard = {"reference": reference, **(standard_updates or {})}
    samples = []
    for index, layout in enumerate(layouts, start=1):
        read1_name = f"S{index}_1" if layout == "PE" else f"S{index}"
        fastq_1 = tmp_path / f"inputs/{read1_name}.fastq.gz"
        row = {
            "sample": f"S{index}",
            "library": "lib1",
            "lane": "L001",
            "layout": layout,
            "fastq_1": str(fastq_1),
            "strandedness": "auto",
            "platform": "ILLUMINA",
        }
        _write(fastq_1, f"R1-{index}".encode())
        if layout == "PE":
            fastq_2 = tmp_path / f"inputs/S{index}_2.fastq.gz"
            _write(fastq_2, f"R2-{index}".encode())
            row["fastq_2"] = str(fastq_2)
        samples.append(row)
    return WorkflowInputs(config={"standard": standard}, samples=samples, options={})


@pytest.fixture
def composed_runtime(tmp_path: Path, monkeypatch):
    root = tmp_path / "runtime"
    binding = BulkRnaSeqExecutionBinding(assets=RuntimeAssetBinding(root=root))
    containers = (
        VerifiedContainerAsset(
            process="FASTQC",
            image="quay.io/biocontainers/fastqc:0.12.1@sha256:" + "1" * 64,
            oci_digest="sha256:" + "1" * 64,
            local_asset=root / "containers/assets/fastqc.tar",
            local_sha256="2" * 64,
            size_bytes=1,
            distribution_manifest=root / "containers/assets/fastqc.manifest.json",
            distribution_manifest_sha256="1" * 64,
            config_digest="sha256:" + "a" * 64,
            runtime_image="sha256:" + "a" * 64,
            rootfs_diff_ids=("sha256:" + "b" * 64,),
        ),
        VerifiedContainerAsset(
            process="STAR_ALIGN",
            image="quay.io/biocontainers/star:2.7.11b@sha256:" + "3" * 64,
            oci_digest="sha256:" + "3" * 64,
            local_asset=root / "containers/assets/star.tar",
            local_sha256="4" * 64,
            size_bytes=1,
            distribution_manifest=root / "containers/assets/star.manifest.json",
            distribution_manifest_sha256="3" * 64,
            config_digest="sha256:" + "c" * 64,
            runtime_image="sha256:" + "c" * 64,
            rootfs_diff_ids=("sha256:" + "d" * 64,),
        ),
        VerifiedContainerAsset(
            process="STAR_GENOMEGENERATE",
            image=(
                "community.wave.seqera.io/library/htslib_samtools_star_gawk:"
                + "ae438e9a604351a4@sha256:"
                + "e" * 64
            ),
            oci_digest="sha256:" + "e" * 64,
            local_asset=root / "containers/assets/star-generate.tar",
            local_sha256="f" * 64,
            size_bytes=1,
            distribution_manifest=(
                root / "containers/assets/star-generate.manifest.json"
            ),
            distribution_manifest_sha256="e" * 64,
            config_digest="sha256:" + "1" * 64,
            runtime_image="sha256:" + "1" * 64,
            rootfs_diff_ids=("sha256:" + "2" * 64,),
        ),
    )
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
        containers=containers,
        source_tree_sha256="5" * 64,
        runtime_identity_sha256="8" * 64,
        nextflow_sha256="9" * 64,
        jdk_archive_sha256="c" * 64,
        jdk_tree_sha256="d" * 64,
        java_executable_sha256="e" * 64,
        plugin_archive_sha256="a" * 64,
        plugin_tree_sha256="6" * 64,
        container_inventory_sha256="b" * 64,
        container_lock_sha256="7" * 64,
    )
    monkeypatch.setattr(
        execution_module,
        "_acquire_runtime_assets",
        lambda _binding: Result.success(verified),
    )
    return binding, verified


def _file_bytes(plan, path: str) -> bytes:
    return dict(plan.files)[path]


def _assert_samplesheet_and_params(
    plan: WorkspacePlan,
    workspace: Path,
    expected_rows: list[str],
) -> None:
    assert _file_bytes(plan, "config/samplesheet.csv").decode().splitlines() == [
        "sample,fastq_1,fastq_2,strandedness,seq_platform",
        *expected_rows,
    ]
    params = json.loads(_file_bytes(plan, "config/params.json"))
    assert params["input"] == str(workspace / "config/samplesheet.csv")
    assert params["outdir"] == str(workspace / "results")
    assert params["custom_config_base"] == ""
    assert params["igenomes_ignore"] is True
    assert params["monochrome_logs"] is True
    assert "seq_platform" not in params


def test_runtime_capabilities_are_truthful_and_default_remains_contract_only(
    composed_runtime,
):
    binding, _ = composed_runtime

    assert BulkRnaSeqWorkflowAdapter().capabilities.supports == (
        "validation",
        "input_authoring",
    )
    assert BulkRnaSeqWorkflowAdapter(execution=binding).capabilities.supports == (
        "validation",
        "input_authoring",
        "workspace_plan",
        "command",
    )


def test_execution_binding_owns_one_process_local_runtime_admission(
    tmp_path: Path,
) -> None:
    assets = RuntimeAssetBinding(root=(tmp_path / "runtime").resolve())
    first = BulkRnaSeqExecutionBinding(assets=assets)
    worker = BulkRnaSeqExecutionBinding(assets=assets)

    assert isinstance(first.runtime_admission, RuntimeAssetAdmission)
    assert first.runtime_admission.binding is assets
    assert worker.runtime_admission is not first.runtime_admission


def test_build_plan_and_command_share_the_binding_admission(
    tmp_path: Path,
    composed_runtime,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    binding, verified = composed_runtime
    observed: list[RuntimeAssetAdmission] = []

    def acquire(execution_binding: BulkRnaSeqExecutionBinding):
        observed.append(execution_binding.runtime_admission)
        return Result.success(verified)

    monkeypatch.setattr(execution_module, "_acquire_runtime_assets", acquire)
    adapter = BulkRnaSeqWorkflowAdapter(execution=binding)
    workspace = (tmp_path / "workspace").resolve()

    assert adapter.capture_build_identity().is_success
    plan = adapter.plan_workspace(_inputs(tmp_path), workspace)
    assert plan.is_success
    assert adapter.build_command(plan.value, workspace).is_success

    assert len(observed) == 3
    assert all(item is binding.runtime_admission for item in observed)


def test_multiqc_identity_conflict_fails_before_runtime_admission(
    tmp_path: Path,
    composed_runtime,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    binding, _ = composed_runtime
    inputs = _inputs(tmp_path, layouts=("PE", "SE"))
    inputs.samples[0]["sample"] = "A"
    inputs.samples[1]["sample"] = "A_1"

    def acquire(_binding: BulkRnaSeqExecutionBinding):
        pytest.fail("runtime admission must not run for invalid sample identities")

    monkeypatch.setattr(execution_module, "_acquire_runtime_assets", acquire)

    result = BulkRnaSeqWorkflowAdapter(execution=binding).plan_workspace(
        inputs,
        (tmp_path / "workspace").resolve(),
    )

    assert result.is_failure
    assert [issue.code for issue in result.issues] == [
        "BULK_RNASEQ_MULTIQC_SAMPLE_IDENTITY_CONFLICT"
    ]


def test_composed_adapter_passes_reusable_conformance(tmp_path: Path, composed_runtime):
    binding, _ = composed_runtime
    valid = _inputs(tmp_path)
    invalid = _inputs(tmp_path / "invalid")
    invalid.config["standard"]["reference"]["fasta_sha256"] = "bad"

    verify_adapter_conformance(
        AdapterConformanceCase(
            adapter=BulkRnaSeqWorkflowAdapter(execution=binding),
            valid_inputs=valid,
            invalid_inputs=invalid,
            planning_workspace=(tmp_path / "workspace").resolve(),
            artifact_workspace=(tmp_path / "artifacts").resolve(),
        )
    )


def test_results_adapter_passes_artifact_and_qc_conformance(
    tmp_path: Path,
    composed_runtime,
):
    binding, _ = composed_runtime
    standard_updates = {
        "trimming": {"enabled": False, "tool": "trimgalore"},
        "qc": {"enabled": False},
        "outputs": {"bigwig": False},
    }
    valid = _inputs(tmp_path, standard_updates=standard_updates)
    invalid = _inputs(tmp_path / "invalid", standard_updates=standard_updates)
    invalid.config["standard"]["reference"]["fasta_sha256"] = "bad"
    artifact_workspace = (tmp_path / "artifacts").resolve()
    artifacts = (
        "pipeline_info/nf_core_rnaseq_software_mqc_versions.yml",
        "star_salmon/S1.sorted.bam",
        "star_salmon/S1.sorted.bam.bai",
        "star_salmon/log/S1.Log.final.out",
        "star_salmon/log/S1.SJ.out.tab",
        "star_salmon/S1/quant.sf",
        "star_salmon/S1/quant.genes.sf",
        "star_salmon/S1/aux_info/meta_info.json",
        "star_salmon/salmon.merged.gene_counts.tsv",
        "star_salmon/salmon.merged.gene_counts_length_scaled.tsv",
        "star_salmon/salmon.merged.gene_counts_scaled.tsv",
        "star_salmon/salmon.merged.gene_lengths.tsv",
        "star_salmon/salmon.merged.gene_tpm.tsv",
        "star_salmon/salmon.merged.transcript_counts.tsv",
        "star_salmon/salmon.merged.transcript_lengths.tsv",
        "star_salmon/salmon.merged.transcript_tpm.tsv",
        "star_salmon/salmon.merged.tx2gene.tsv",
        "star_salmon/salmon.merged.tx2gene_augmented.tsv",
    )
    for relative_path in artifacts:
        _write(artifact_workspace / "results" / relative_path, b"fixture")
    star_content = b"""\
Number of input reads | 1000
Uniquely mapped reads number | 800
Uniquely mapped reads % | 80.00%
Number of reads mapped to multiple loci | 100
% of reads mapped to multiple loci | 10.00%
% of reads mapped to too many loci | 2.00%
% of reads unmapped: too many mismatches | 1.00%
% of reads unmapped: too short | 6.00%
% of reads unmapped: other | 1.00%
"""
    salmon_content = json.dumps(
        {
            "salmon_version": "1.10.3",
            "mapping_type": "alignment",
            "num_libraries": 1,
            "num_processed": 1000,
            "num_mapped": 800,
            "percent_mapped": 80,
        },
        separators=(",", ":"),
    ).encode()
    source_metadata = {"scope": "sample", "sample_id": "S1", "assay": "bulk-rnaseq"}
    qc_sources = (
        QcSourceDocument(
            source=QcSourceArtifact(
                artifact_id="artifact-star",
                output_type="bulk_rnaseq.star.log_final",
                relative_path="results/star_salmon/log/S1.Log.final.out",
                metadata=source_metadata,
            ),
            content=star_content,
        ),
        QcSourceDocument(
            source=QcSourceArtifact(
                artifact_id="artifact-salmon",
                output_type="bulk_rnaseq.salmon.meta_info",
                relative_path="results/star_salmon/S1/aux_info/meta_info.json",
                metadata=source_metadata,
            ),
            content=salmon_content,
        ),
    )

    verify_adapter_conformance(
        AdapterConformanceCase(
            adapter=BulkRnaSeqResultsWorkflowAdapter(execution=binding),
            valid_inputs=valid,
            invalid_inputs=invalid,
            planning_workspace=(tmp_path / "workspace").resolve(),
            artifact_workspace=artifact_workspace,
            qc_sources=qc_sources,
        )
    )


def test_results_composition_is_bound_into_build_and_workspace_identity(
    tmp_path: Path,
    composed_runtime,
):
    binding, _ = composed_runtime
    runtime_adapter = BulkRnaSeqWorkflowAdapter(execution=binding)
    results_adapter = BulkRnaSeqResultsWorkflowAdapter(execution=binding)

    runtime_build = runtime_adapter.capture_build_identity()
    results_build = results_adapter.capture_build_identity()

    assert runtime_build.is_success
    assert results_build.is_success
    assert runtime_build.value.digest != results_build.value.digest

    workspace = (tmp_path / "workspace").resolve()
    results_plan = results_adapter.plan_workspace(_inputs(tmp_path), workspace)
    assert results_plan.is_success
    assert results_adapter.build_command(results_plan.value, workspace).is_success

    mismatched = runtime_adapter.build_command(results_plan.value, workspace)
    assert mismatched.is_failure
    assert mismatched.errors[0].code == "BULK_RNASEQ_COMMAND_IDENTITY_MISMATCH"


@pytest.mark.parametrize("layout", ["SE", "PE"])
def test_workspace_serializes_single_layout_with_per_row_platform(
    layout: str, tmp_path: Path, composed_runtime
):
    binding, _ = composed_runtime
    inputs = _inputs(tmp_path, layouts=(layout,))
    adapter = BulkRnaSeqWorkflowAdapter(execution=binding)
    workspace = (tmp_path / "workspace").resolve()

    result = adapter.plan_workspace(inputs, workspace)

    assert result.is_success
    fastq_1 = (
        f"{tmp_path}/inputs/S1_1.fastq.gz"
        if layout == "PE"
        else f"{tmp_path}/inputs/S1.fastq.gz"
    )
    fastq_2 = f"{tmp_path}/inputs/S1_2.fastq.gz" if layout == "PE" else ""
    _assert_samplesheet_and_params(
        result.value,
        workspace,
        [
            f"S1,{fastq_1},{fastq_2},auto,ILLUMINA",
        ],
    )


def test_workspace_serializes_repeated_lanes_with_per_row_platform(
    tmp_path: Path, composed_runtime
):
    binding, _ = composed_runtime
    inputs = _inputs(tmp_path)
    second_lane = dict(inputs.samples[0])
    second_lane["lane"] = "L002"
    second_lane["fastq_1"] = str(tmp_path / "inputs/S1_1.lib1.L002.fastq.gz")
    second_lane["fastq_2"] = str(tmp_path / "inputs/S1_2.lib1.L002.fastq.gz")
    _write(Path(second_lane["fastq_1"]), b"R1-lane-2")
    _write(Path(second_lane["fastq_2"]), b"R2-lane-2")
    inputs.samples.append(second_lane)
    workspace = (tmp_path / "workspace").resolve()

    result = BulkRnaSeqWorkflowAdapter(execution=binding).plan_workspace(
        inputs, workspace
    )

    assert result.is_success
    _assert_samplesheet_and_params(
        result.value,
        workspace,
        [
            (
                f"S1,{tmp_path}/inputs/S1_1.fastq.gz,"
                f"{tmp_path}/inputs/S1_2.fastq.gz,auto,ILLUMINA"
            ),
            (
                f"S1,{tmp_path}/inputs/S1_1.lib1.L002.fastq.gz,"
                f"{tmp_path}/inputs/S1_2.lib1.L002.fastq.gz,auto,ILLUMINA"
            ),
        ],
    )


def test_workspace_serializes_mixed_layout_deterministically(
    tmp_path: Path, composed_runtime
):
    binding, _ = composed_runtime
    inputs = _inputs(tmp_path, layouts=("SE", "PE"))
    adapter = BulkRnaSeqWorkflowAdapter(execution=binding)
    workspace = (tmp_path / "workspace").resolve()

    first = adapter.plan_workspace(inputs, workspace)
    second = adapter.plan_workspace(inputs, workspace)

    assert first.is_success
    assert second.is_success
    assert first.value == second.value
    _assert_samplesheet_and_params(
        first.value,
        workspace,
        [
            f"S1,{tmp_path}/inputs/S1.fastq.gz,,auto,ILLUMINA",
            (
                f"S2,{tmp_path}/inputs/S2_1.fastq.gz,"
                f"{tmp_path}/inputs/S2_2.fastq.gz,auto,ILLUMINA"
            ),
        ],
    )


def test_workspace_maps_umi_without_caller_runtime_tokens(
    tmp_path: Path, composed_runtime
):
    binding, _ = composed_runtime
    inputs = _inputs(
        tmp_path,
        standard_updates={
            "umi": {
                "enabled": True,
                "mode": "read_sequence",
                "deduplication_tool": "umitools",
                "extraction_method": "string",
                "barcode_pattern": "NNNNNN",
            }
        },
    )

    result = BulkRnaSeqWorkflowAdapter(execution=binding).plan_workspace(
        inputs, (tmp_path / "workspace").resolve()
    )

    assert result.is_success
    params = json.loads(_file_bytes(result.value, "config/params.json"))
    assert params["with_umi"] is True
    assert params["umitools_bc_pattern"] == "NNNNNN"
    assert "nextflow_args" not in params
    assert "resume" not in params


def test_workspace_rejects_prebuilt_index_incompatible_with_run_shaping_params(
    tmp_path: Path,
    composed_runtime,
):
    binding, verified = composed_runtime
    inputs = _inputs(tmp_path)
    reference = inputs.config["standard"]["reference"]
    index_root = tmp_path / "inputs/star-index"
    index_file = index_root / "index/Genome"
    index_sha256 = _write(index_file, b"tiny-star-index")
    producer_image = next(
        item.image
        for item in verified.containers
        if item.process == "STAR_GENOMEGENERATE"
    )
    manifest = {
        "schema_version": "1.0.0",
        "index_kind": "star",
        "producer": {
            "process": "STAR_GENOMEGENERATE",
            "tool": "star",
            "tool_version": "2.7.11b",
            "build_contract": "nfcore-rnaseq-3.26.0-star-genomegenerate-default-v1",
            "container_image": producer_image,
        },
        "build_parameters": {
            "extra_args": [],
            "genome_sa_index_nbases": 9,
            "genome_sa_index_nbases_strategy": "nfcore-rnaseq-3.26.0-auto",
            "sjdb_gtf_feature": "exon",
            "sjdb_overhang": 100,
            "skip_gtf_filter": False,
            "skip_gtf_transcript_filter": False,
        },
        "reference": {
            "fasta_sha256": reference["fasta_sha256"],
            "gtf_sha256": reference["gtf_sha256"],
        },
        "files": [
            {
                "path": "index/Genome",
                "size_bytes": index_file.stat().st_size,
                "sha256": index_sha256,
            }
        ],
    }
    manifest_bytes = json.dumps(
        manifest,
        sort_keys=True,
        separators=(",", ":"),
    ).encode()
    manifest_sha256 = _write(
        index_root / REFERENCE_INDEX_MANIFEST,
        manifest_bytes,
    )
    reference["star_index"] = {
        "path": str(index_root),
        "identity_sha256": manifest_sha256,
    }
    inputs.config["advanced"] = {"skip_gtf_filter": True}

    result = BulkRnaSeqWorkflowAdapter(execution=binding).plan_workspace(
        inputs,
        (tmp_path / "workspace").resolve(),
    )

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_REFERENCE_INVALID"


def test_command_owns_nextflow_paths_profile_reports_and_no_pull(
    tmp_path: Path, composed_runtime
):
    binding, verified = composed_runtime
    adapter = BulkRnaSeqWorkflowAdapter(execution=binding)
    workspace = (tmp_path / "workspace").resolve()
    plan = adapter.plan_workspace(_inputs(tmp_path), workspace).value

    result = adapter.build_command(plan, workspace)

    assert result.is_success
    command = result.value
    assert command.argv[0] == str(verified.nextflow_executable)
    assert command.preflight_argv[0] == command.argv[0]
    assert command.argv.count("-C") == 1
    assert command.preflight_argv.count("-C") == 1
    assert "-c" not in command.argv
    assert "-c" not in command.preflight_argv
    expected_config = str(workspace / "config/platform.nextflow.config")
    assert command.argv[command.argv.index("-C") + 1] == expected_config
    assert (
        command.preflight_argv[command.preflight_argv.index("-C") + 1]
        == expected_config
    )
    assert command.argv.index("-C") < command.argv.index("run")
    assert command.preflight_argv.index("-C") < command.preflight_argv.index("config")
    assert command.preflight_kind == "configuration"
    assert command.preflight_managed_logs == (
        ("nextflow_preflight", str(workspace / "logs/nextflow-preflight.log")),
    )
    assert command.execution_managed_logs == (
        ("nextflow", str(workspace / "logs/nextflow.log")),
    )
    assert command.cwd == str(workspace / "engine/launch")
    assert command.env["NXF_OFFLINE"] == "true"
    assert command.env["NXF_HOME"] == str(workspace / "engine/nxf-home")
    assert command.env["NXF_PLUGINS_DIR"] == str(verified.plugin_tree.parent)
    assert command.env["JAVA_HOME"] == str(verified.jdk_tree)
    assert command.env["NXF_JAVA_HOME"] == str(verified.jdk_tree)
    assert command.env["JAVA_CMD"] == str(verified.java_executable)
    assert command.env["PATH"] == f"{verified.jdk_tree / 'bin'}:/usr/bin:/bin"
    assert command.env["LD_LIBRARY_PATH"] == ""
    assert command.env["CONDA_PREFIX"] == ""
    assert command.env["MAMBA_ROOT_PREFIX"] == ""
    assert command.argv[command.argv.index("run") + 1] == str(verified.source_tree)
    assert command.argv[command.argv.index("-profile") + 1] == "docker"
    assert "-resume" not in command.argv
    assert not any("nf-core/rnaseq" == token for token in command.argv)
    config = _file_bytes(plan, "config/platform.nextflow.config").decode()
    effective_lines = [
        line for line in config.splitlines() if line and not line.startswith("//")
    ]
    source_include = f"includeConfig '{verified.source_tree}/nextflow.config'"
    assert effective_lines[0] == source_include
    assert config.count("includeConfig ") == 1
    assert config.index(source_include) < config.index("process.executor = 'local'")
    assert "--pull=never --network=none" in config
    assert "process.stageInMode = 'copy'" in config
    assert "wave.enabled = false" in config
    assert "tower.enabled = false" in config
    assert "withName: 'FASTQC'" in config
    assert verified.containers[0].runtime_image in config
    assert "docker.registry = ''" in config
    assert "docker.remove = true" in config
    assert "withLabel: '!helixweave_verified_container_v1'" in config
    assert f"container = 'sha256:{'0' * 64}'" in config
    parabricks_deny = (
        "    withName: '.*PARABRICKS_STARGENOMEGENERATE' {\n"
        f"        container = 'sha256:{'0' * 64}'\n"
        "    }"
    )
    assert config.count(parabricks_deny) == 1
    assert config.index("withName: 'STAR_GENOMEGENERATE'") < config.index(
        parabricks_deny
    )
    assert "htslib_samtools_star_gawk:4de2f983041d42e6" not in config
    assert verified.containers[0].image not in config
    execution_identity = json.loads(_file_bytes(plan, "config/execution-identity.json"))
    cache_identity = json.loads(_file_bytes(plan, "engine/cache-identity.json"))
    captured = adapter.capture_build_identity()
    assert captured.is_success
    assert execution_identity["build_identity_sha256"] == captured.value.digest
    assert cache_identity["workflow_build_sha256"] == captured.value.digest
    assert execution_identity["adapter_version"] == captured.value.adapter_version
    assert (
        cache_identity["execution_implementation_manifest_sha256"]
        == execution_identity["execution_implementation_manifest_sha256"]
    )
    assert (
        cache_identity["execution_implementation_aggregate_sha256"]
        == execution_identity["execution_implementation_aggregate_sha256"]
    )
    workspace_identity = execution_identity["workspace_identity_sha256"]
    assert len(workspace_identity) == 64
    assert f"--label=org.helixweave.workspace-scope={workspace_identity}" in config


def test_command_environment_overrides_host_jvm_and_conda_poison(
    tmp_path: Path,
    composed_runtime,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    binding, verified = composed_runtime
    for name, value in {
        "JAVA_HOME": "/host/java",
        "JAVA_TOOL_OPTIONS": "-javaagent:/host/poison.jar",
        "JDK_JAVA_OPTIONS": "-Dhost.poison=true",
        "LD_LIBRARY_PATH": "/host/conda/lib",
        "PATH": "/host/conda/bin:/usr/bin",
        "CONDA_PREFIX": "/host/conda",
        "CONDA_DEFAULT_ENV": "host",
        "CONDA_EXE": "/host/conda/bin/conda",
        "MAMBA_ROOT_PREFIX": "/host/mamba",
        "_JAVA_OPTIONS": "-Xbootclasspath/a:/host/poison.jar",
    }.items():
        monkeypatch.setenv(name, value)
    workspace = (tmp_path / "workspace").resolve()
    adapter = BulkRnaSeqWorkflowAdapter(execution=binding)
    plan = adapter.plan_workspace(_inputs(tmp_path), workspace).value
    command = adapter.build_command(plan, workspace).value

    environment = _subprocess_environment(command.env)

    assert environment.is_success
    assert environment.value["JAVA_HOME"] == str(verified.jdk_tree)
    assert environment.value["NXF_JAVA_HOME"] == str(verified.jdk_tree)
    assert environment.value["JAVA_CMD"] == str(verified.java_executable)
    assert environment.value["PATH"] == (f"{verified.jdk_tree / 'bin'}:/usr/bin:/bin")
    assert environment.value["LD_LIBRARY_PATH"] == ""
    assert environment.value["CONDA_PREFIX"] == ""
    assert environment.value["CONDA_DEFAULT_ENV"] == ""
    assert environment.value["CONDA_EXE"] == ""
    assert environment.value["MAMBA_ROOT_PREFIX"] == ""
    assert "JAVA_TOOL_OPTIONS" not in environment.value
    assert "JDK_JAVA_OPTIONS" not in environment.value
    assert "_JAVA_OPTIONS" not in environment.value


def test_same_inputs_in_distinct_workspaces_get_unique_server_run_labels(
    tmp_path: Path,
    composed_runtime,
):
    binding, _ = composed_runtime
    adapter = BulkRnaSeqWorkflowAdapter(execution=binding)
    inputs = _inputs(tmp_path)

    first = adapter.plan_workspace(inputs, (tmp_path / "workspace-a").resolve())
    second = adapter.plan_workspace(inputs, (tmp_path / "workspace-b").resolve())

    assert first.is_success and second.is_success
    first_identity = json.loads(
        _file_bytes(first.value, "config/execution-identity.json")
    )
    second_identity = json.loads(
        _file_bytes(second.value, "config/execution-identity.json")
    )
    assert (
        first_identity["input_identity_sha256"]
        == second_identity["input_identity_sha256"]
    )
    assert (
        first_identity["build_identity_sha256"]
        == second_identity["build_identity_sha256"]
    )
    assert (
        first_identity["workspace_identity_sha256"]
        != second_identity["workspace_identity_sha256"]
    )
    first_config = _file_bytes(first.value, "config/platform.nextflow.config").decode()
    second_config = _file_bytes(
        second.value, "config/platform.nextflow.config"
    ).decode()
    assert (
        "--label=org.helixweave.workspace-scope="
        f"{first_identity['workspace_identity_sha256']}" in first_config
    )
    assert (
        "--label=org.helixweave.workspace-scope="
        f"{second_identity['workspace_identity_sha256']}" in second_config
    )


def test_command_redacts_workspace_runtime_and_submitted_paths(
    tmp_path: Path, composed_runtime
):
    binding, verified = composed_runtime
    adapter = BulkRnaSeqWorkflowAdapter(execution=binding)
    workspace = (tmp_path / "workspace").resolve()
    inputs = _inputs(tmp_path)
    plan = adapter.plan_workspace(inputs, workspace).value

    command = adapter.build_command(plan, workspace).value
    serialized = command.to_dict()

    assert str(workspace) in command.redaction_values
    assert str(verified.root) in command.redaction_values
    assert str(verified.jdk_archive) in command.redaction_values
    assert str(verified.jdk_tree) in command.redaction_values
    assert str(verified.java_executable) in command.redaction_values
    assert inputs.samples[0]["fastq_1"] in command.redaction_values
    assert command.managed_container_scope == managed_container_scope(workspace)
    assert (
        command.managed_container_endpoint_identity
        == managed_container_endpoint_identity(
            binding.assets.docker_executable,
            binding.assets.docker_socket,
        )
    )
    assert serialized["redaction_count"] == len(command.redaction_values)
    assert serialized["managed_container_cleanup"] is True
    assert str(workspace) not in json.dumps(serialized)
    assert str(verified.jdk_tree) not in json.dumps(serialized)
    assert inputs.samples[0]["fastq_1"] not in json.dumps(serialized)


def test_command_rejects_tampered_container_audit_identity(
    tmp_path: Path, composed_runtime
):
    binding, _verified = composed_runtime
    adapter = BulkRnaSeqWorkflowAdapter(execution=binding)
    workspace = (tmp_path / "workspace").resolve()
    plan = adapter.plan_workspace(_inputs(tmp_path), workspace).value
    identity = json.loads(_file_bytes(plan, "config/execution-identity.json"))
    identity["container_process_audit_sha256"] = "0" * 64
    replacement = (
        json.dumps(identity, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()
    tampered = WorkspacePlan(
        directories=plan.directories,
        files=tuple(
            (path, replacement if path == "config/execution-identity.json" else value)
            for path, value in plan.files
        ),
    )

    result = adapter.build_command(tampered, workspace)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_COMMAND_IDENTITY_MISMATCH"


def test_command_rejects_tampered_workspace_contract_bytes(
    tmp_path: Path, composed_runtime
):
    binding, _verified = composed_runtime
    adapter = BulkRnaSeqWorkflowAdapter(execution=binding)
    workspace = (tmp_path / "workspace").resolve()
    plan = adapter.plan_workspace(_inputs(tmp_path), workspace).value
    tampered = WorkspacePlan(
        directories=plan.directories,
        files=tuple(
            (path, b"{}\n" if path == "config/params.json" else value)
            for path, value in plan.files
        ),
    )

    result = adapter.build_command(tampered, workspace)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_COMMAND_IDENTITY_MISMATCH"


def test_build_identity_changes_with_container_lock_or_process_audit(
    composed_runtime, monkeypatch
):
    binding, verified = composed_runtime
    adapter = BulkRnaSeqWorkflowAdapter(execution=binding)
    first = adapter.capture_build_identity()
    changed = VerifiedRuntimeAssets(
        **{
            **verified.__dict__,
            "container_lock_sha256": "8" * 64,
        }
    )
    monkeypatch.setattr(
        execution_module,
        "_acquire_runtime_assets",
        lambda _binding: Result.success(changed),
    )

    second = adapter.capture_build_identity()

    assert first.is_success and second.is_success
    assert first.value.scheme == "sha256-bulk-rnaseq-runtime-v1"
    assert first.value.logical_entrypoint == "main.nf"
    assert first.value.digest != second.value.digest

    changed_audit = VerifiedRuntimeAssets(
        **{
            **verified.__dict__,
            "container_process_audit_sha256": "e" * 64,
        }
    )
    monkeypatch.setattr(
        execution_module,
        "_acquire_runtime_assets",
        lambda _binding: Result.success(changed_audit),
    )
    third = adapter.capture_build_identity()

    assert third.is_success
    assert len({first.value.digest, second.value.digest, third.value.digest}) == 3


def test_jdk_identity_changes_build_and_cache_identity(
    tmp_path: Path,
    composed_runtime,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    binding, verified = composed_runtime
    adapter = BulkRnaSeqWorkflowAdapter(execution=binding)
    workspace = (tmp_path / "workspace").resolve()
    inputs = _inputs(tmp_path)
    first = adapter.plan_workspace(inputs, workspace)
    changed = VerifiedRuntimeAssets(
        **{
            **verified.__dict__,
            "jdk_tree_sha256": "0" * 64,
        }
    )
    monkeypatch.setattr(
        execution_module,
        "_acquire_runtime_assets",
        lambda _binding: Result.success(changed),
    )

    second = adapter.plan_workspace(inputs, workspace)

    assert first.is_success and second.is_success
    first_cache = json.loads(_file_bytes(first.value, "engine/cache-identity.json"))
    second_cache = json.loads(_file_bytes(second.value, "engine/cache-identity.json"))
    assert first_cache["workflow_build_sha256"] != second_cache["workflow_build_sha256"]
    assert first_cache["identity_sha256"] != second_cache["identity_sha256"]


def test_resume_is_server_owned_and_unavailable_without_attempt_lifecycle(
    tmp_path: Path,
):
    assets = RuntimeAssetBinding(root=(tmp_path / "runtime").resolve())

    with pytest.raises(ValueError, match="attempt/session lifecycle"):
        BulkRnaSeqExecutionBinding(assets=assets, resume_enabled=True)


def test_workspace_rejects_missing_fastq_and_symlinked_reference(
    tmp_path: Path, composed_runtime
):
    binding, _ = composed_runtime
    adapter = BulkRnaSeqWorkflowAdapter(execution=binding)
    inputs = _inputs(tmp_path)
    Path(inputs.samples[0]["fastq_1"]).unlink()

    missing = adapter.plan_workspace(inputs, (tmp_path / "workspace").resolve())

    assert missing.is_failure
    assert missing.errors[0].code == "BULK_RNASEQ_FASTQ_INVALID"

    inputs = _inputs(tmp_path / "linked")
    fasta = Path(inputs.config["standard"]["reference"]["fasta"])
    real_fasta = fasta.with_name("real.fa")
    fasta.rename(real_fasta)
    fasta.symlink_to(real_fasta)
    linked = adapter.plan_workspace(inputs, (tmp_path / "workspace2").resolve())
    assert linked.is_failure
    assert linked.errors[0].code == "BULK_RNASEQ_REFERENCE_INVALID"
