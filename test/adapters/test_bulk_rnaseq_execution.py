"""Workspace and command tests for the offline bulk RNA-seq runtime binding."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

import encode_pipeline.adapters.bulk_rnaseq.execution as execution_module
from encode_pipeline.adapters.bulk_rnaseq import (
    BulkRnaSeqExecutionBinding,
    BulkRnaSeqRapidQuantQualificationAdapter,
    BulkRnaSeqResultsWorkflowAdapter,
    BulkRnaSeqTranscriptomeBinding,
    BulkRnaSeqWorkflowAdapter,
    RuntimeAssetBinding,
)
from encode_pipeline.adapters.bulk_rnaseq.execution_identity import (
    VerifiedExecutionImplementation,
)
from encode_pipeline.adapters.bulk_rnaseq.qualification import (
    RAPID_QUANT_MODE_OWNED_PARAMETERS,
    RAPID_QUANT_PROFILE_OWNED_PARAMETERS,
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
    MANAGED_CONTAINER_SCOPE_LABEL,
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


def _transcriptome_binding(tmp_path: Path) -> BulkRnaSeqTranscriptomeBinding:
    transcript_fasta = tmp_path / "inputs/transcripts.fa"
    transcript_fasta_sha256 = _write(transcript_fasta, b">tx1\nACGT\n")
    return BulkRnaSeqTranscriptomeBinding(
        reference_id="tiny-ref",
        fasta_sha256=_sha256(b">chr1\nACGT\n"),
        gtf_sha256=_sha256(b'chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n'),
        transcript_fasta=transcript_fasta,
        transcript_fasta_sha256=transcript_fasta_sha256,
    )


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
    binding = BulkRnaSeqExecutionBinding(
        assets=RuntimeAssetBinding(root=root),
        transcriptome=_transcriptome_binding(tmp_path),
    )
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


@pytest.fixture
def admitted_execution_implementation(monkeypatch: pytest.MonkeyPatch):
    implementation = VerifiedExecutionImplementation(
        manifest_sha256="3" * 64,
        aggregate_sha256="4" * 64,
        files=(),
    )
    monkeypatch.setattr(
        execution_module,
        "verify_execution_implementation",
        lambda: Result.success(implementation),
    )
    return implementation


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
    assert params["pipelines_testdata_base_path"] == str(workspace / "config")
    assert params["igenomes_base"] == str(workspace / "config")
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
    transcriptome = _transcriptome_binding(tmp_path)
    first = BulkRnaSeqExecutionBinding(
        assets=assets,
        transcriptome=transcriptome,
    )
    worker = BulkRnaSeqExecutionBinding(
        assets=assets,
        transcriptome=transcriptome,
    )

    assert isinstance(first.runtime_admission, RuntimeAssetAdmission)
    assert first.runtime_admission.binding is assets
    assert worker.runtime_admission is not first.runtime_admission


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("reference_id", "bad reference"),
        ("fasta_sha256", "bad"),
        ("gtf_sha256", "bad"),
        ("transcript_fasta", Path("relative.fa")),
        ("transcript_fasta_sha256", "bad"),
    ],
)
def test_transcriptome_binding_is_closed_and_canonical(
    tmp_path: Path,
    field: str,
    value: object,
) -> None:
    document = {
        "reference_id": "tiny-ref",
        "fasta_sha256": "a" * 64,
        "gtf_sha256": "b" * 64,
        "transcript_fasta": (tmp_path / "transcripts.fa").resolve(),
        "transcript_fasta_sha256": "c" * 64,
    }
    document[field] = value

    with pytest.raises(ValueError):
        BulkRnaSeqTranscriptomeBinding(**document)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("reference_id", "another-ref"),
        ("fasta_sha256", "a" * 64),
        ("gtf_sha256", "b" * 64),
    ],
)
def test_workspace_rejects_each_transcriptome_reference_binding_mismatch(
    tmp_path: Path,
    composed_runtime,
    admitted_execution_implementation,
    field: str,
    value: str,
) -> None:
    binding, _verified = composed_runtime
    transcriptome_document = {
        name: getattr(binding.transcriptome, name)
        for name in BulkRnaSeqTranscriptomeBinding.__dataclass_fields__
    }
    transcriptome_document[field] = value
    mismatched_binding = BulkRnaSeqExecutionBinding(
        assets=binding.assets,
        transcriptome=BulkRnaSeqTranscriptomeBinding(**transcriptome_document),
    )
    adapter = BulkRnaSeqWorkflowAdapter(execution=mismatched_binding)
    inputs = _inputs(tmp_path)
    workspace = (tmp_path / "workspace").resolve()

    result = adapter.plan_workspace(inputs, workspace)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_TRANSCRIPTOME_REFERENCE_MISMATCH"


def test_workspace_rejects_transcriptome_content_drift(
    tmp_path: Path,
    composed_runtime,
    admitted_execution_implementation,
) -> None:
    binding, _verified = composed_runtime
    adapter = BulkRnaSeqWorkflowAdapter(execution=binding)
    inputs = _inputs(tmp_path)
    workspace = (tmp_path / "workspace").resolve()
    binding.transcriptome.transcript_fasta.write_bytes(b">tx1\nTGCA\n")

    result = adapter.plan_workspace(inputs, workspace)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_TRANSCRIPTOME_INVALID"


@pytest.mark.parametrize("path_kind", ["missing", "symlink"])
def test_workspace_rejects_missing_or_symlinked_transcriptome(
    tmp_path: Path,
    composed_runtime,
    admitted_execution_implementation,
    path_kind: str,
) -> None:
    binding, _verified = composed_runtime
    transcript = (tmp_path / f"deployment/{path_kind}/transcripts.fa").resolve()
    expected_sha256 = _sha256(b">tx1\nACGT\n")
    if path_kind == "symlink":
        target = (tmp_path / "deployment/target.fa").resolve()
        _write(target, b">tx1\nACGT\n")
        transcript.parent.mkdir(parents=True)
        transcript.symlink_to(target)
    transcriptome = BulkRnaSeqTranscriptomeBinding(
        reference_id=binding.transcriptome.reference_id,
        fasta_sha256=binding.transcriptome.fasta_sha256,
        gtf_sha256=binding.transcriptome.gtf_sha256,
        transcript_fasta=transcript,
        transcript_fasta_sha256=expected_sha256,
    )
    adapter = BulkRnaSeqWorkflowAdapter(
        execution=BulkRnaSeqExecutionBinding(
            assets=binding.assets,
            transcriptome=transcriptome,
        )
    )

    result = adapter.plan_workspace(
        _inputs(tmp_path),
        (tmp_path / "workspace").resolve(),
    )

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_TRANSCRIPTOME_INVALID"


def test_transcriptome_change_updates_input_and_cache_identity(
    tmp_path: Path,
    composed_runtime,
    admitted_execution_implementation,
) -> None:
    binding, _verified = composed_runtime
    inputs = _inputs(tmp_path)
    workspace = (tmp_path / "workspace").resolve()
    first = BulkRnaSeqWorkflowAdapter(execution=binding).plan_workspace(
        inputs,
        workspace,
    )
    changed_content = b">tx1\nTGCA\n"
    changed_sha256 = _write(
        binding.transcriptome.transcript_fasta,
        changed_content,
    )
    changed_transcriptome = BulkRnaSeqTranscriptomeBinding(
        reference_id=binding.transcriptome.reference_id,
        fasta_sha256=binding.transcriptome.fasta_sha256,
        gtf_sha256=binding.transcriptome.gtf_sha256,
        transcript_fasta=binding.transcriptome.transcript_fasta,
        transcript_fasta_sha256=changed_sha256,
    )
    changed_adapter = BulkRnaSeqWorkflowAdapter(
        execution=BulkRnaSeqExecutionBinding(
            assets=binding.assets,
            transcriptome=changed_transcriptome,
        )
    )

    second = changed_adapter.plan_workspace(inputs, workspace)

    assert first.is_success and second.is_success
    first_execution = json.loads(
        _file_bytes(first.value, "config/execution-identity.json")
    )
    second_execution = json.loads(
        _file_bytes(second.value, "config/execution-identity.json")
    )
    first_cache = json.loads(_file_bytes(first.value, "engine/cache-identity.json"))
    second_cache = json.loads(_file_bytes(second.value, "engine/cache-identity.json"))
    assert (
        first_execution["input_identity_sha256"]
        != second_execution["input_identity_sha256"]
    )
    assert first_cache["identity_sha256"] != second_cache["identity_sha256"]


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
Started job on | Jul 17 00:00:00
Started mapping on | Jul 17 00:00:01
Finished on | Jul 17 00:00:02
Mapping speed, Million of reads per hour | 1800.00
Number of input reads | 1000
Average input read length | 100
UNIQUE READS:
Uniquely mapped reads number | 800
Uniquely mapped reads % | 80.00%
Average mapped length | 99.00
Number of splices: Total | 100
Number of splices: Annotated (sjdb) | 90
Number of splices: GT/AG | 80
Number of splices: GC/AG | 10
Number of splices: AT/AC | 5
Number of splices: Non-canonical | 5
Mismatch rate per base, % | 0.10%
Deletion rate per base | 0.01%
Deletion average length | 1.00
Insertion rate per base | 0.01%
Insertion average length | 1.00
MULTI-MAPPING READS:
Number of reads mapped to multiple loci | 100
% of reads mapped to multiple loci | 10.00%
Number of reads mapped to too many loci | 20
% of reads mapped to too many loci | 2.00%
UNMAPPED READS:
Number of reads unmapped: too many mismatches | 10
% of reads unmapped: too many mismatches | 1.00%
Number of reads unmapped: too short | 60
% of reads unmapped: too short | 6.00%
Number of reads unmapped: other | 10
% of reads unmapped: other | 1.00%
CHIMERIC READS:
Number of chimeric reads | 0
% of chimeric reads | 0.00%
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


def test_rapid_quant_qualification_requires_server_runtime_and_stays_runtime_only(
    tmp_path: Path,
    composed_runtime,
    admitted_execution_implementation,
):
    binding, _ = composed_runtime

    with pytest.raises(ValueError, match="BulkRnaSeqExecutionBinding"):
        BulkRnaSeqRapidQuantQualificationAdapter(  # type: ignore[arg-type]
            execution=None
        )

    adapter = BulkRnaSeqRapidQuantQualificationAdapter(execution=binding)

    assert adapter.capabilities.supports == (
        "validation",
        "input_authoring",
        "workspace_plan",
        "command",
    )
    assert (
        adapter.extract_artifacts(_inputs(tmp_path), (tmp_path / "artifacts").resolve())
        .errors[0]
        .code
        == "BULK_RNASEQ_CAPABILITY_UNSUPPORTED"
    )


def test_rapid_quant_profile_is_exact_and_cannot_be_overridden_by_params_file(
    tmp_path: Path,
    composed_runtime,
    admitted_execution_implementation,
):
    binding, _ = composed_runtime
    inputs = _inputs(tmp_path)
    workspace = (tmp_path / "workspace").resolve()
    standard = BulkRnaSeqWorkflowAdapter(execution=binding)
    rapid = BulkRnaSeqRapidQuantQualificationAdapter(execution=binding)

    standard_plan = standard.plan_workspace(inputs, workspace)
    rapid_plan = rapid.plan_workspace(inputs, workspace)

    assert standard_plan.is_success and rapid_plan.is_success
    standard_params = json.loads(_file_bytes(standard_plan.value, "config/params.json"))
    rapid_params = json.loads(_file_bytes(rapid_plan.value, "config/params.json"))
    assert RAPID_QUANT_PROFILE_OWNED_PARAMETERS == frozenset(
        {
            "pseudo_aligner",
            "skip_alignment",
            "skip_bigwig",
            "skip_biotype_qc",
            "skip_deseq2_qc",
            "skip_dupradar",
            "skip_markduplicates",
            "skip_multiqc",
            "skip_preseq",
            "skip_qualimap",
            "skip_quantification_merge",
            "skip_rseqc",
            "skip_stringtie",
        }
    )
    assert set(
        standard_params
    ) & RAPID_QUANT_PROFILE_OWNED_PARAMETERS == RAPID_QUANT_PROFILE_OWNED_PARAMETERS - {
        "pseudo_aligner"
    }
    assert set(rapid_params).isdisjoint(RAPID_QUANT_PROFILE_OWNED_PARAMETERS)
    assert RAPID_QUANT_MODE_OWNED_PARAMETERS == (
        RAPID_QUANT_PROFILE_OWNED_PARAMETERS | {"skip_pseudo_alignment"}
    )
    assert "skip_pseudo_alignment" not in rapid_params
    assert standard_params["transcript_fasta"] == rapid_params["transcript_fasta"]
    assert standard_params["transcript_fasta"].endswith("/inputs/transcripts.fa")
    assert rapid_params == {
        name: value
        for name, value in standard_params.items()
        if name not in RAPID_QUANT_MODE_OWNED_PARAMETERS
    }
    assert "config_profile_name" not in rapid_params
    assert "config_profile_description" not in rapid_params

    standard_identity = json.loads(
        _file_bytes(standard_plan.value, "config/execution-identity.json")
    )
    rapid_identity = json.loads(
        _file_bytes(rapid_plan.value, "config/execution-identity.json")
    )
    standard_cache = json.loads(
        _file_bytes(standard_plan.value, "engine/cache-identity.json")
    )
    rapid_cache = json.loads(
        _file_bytes(rapid_plan.value, "engine/cache-identity.json")
    )
    assert standard_identity["execution_mode"] == "standard-v1"
    assert rapid_identity["execution_mode"] == "rapid-quant-v1"
    assert standard_cache["execution_mode"] == "standard-v1"
    assert rapid_cache["execution_mode"] == "rapid-quant-v1"
    assert standard_identity["schema_version"] == "1.1.0"
    assert rapid_identity["schema_version"] == "1.1.0"
    assert (
        standard_identity["build_identity_sha256"]
        != rapid_identity["build_identity_sha256"]
    )
    assert (
        standard_identity["workspace_contract_sha256"]
        != rapid_identity["workspace_contract_sha256"]
    )
    assert standard_cache["identity_sha256"] != rapid_cache["identity_sha256"]

    standard_command = standard.build_command(standard_plan.value, workspace)
    command = rapid.build_command(rapid_plan.value, workspace)

    assert standard_command.is_success and command.is_success
    assert "-profile" not in standard_command.value.argv
    assert "-profile" not in standard_command.value.preflight_argv
    assert command.value.argv[command.value.argv.index("-profile") + 1] == (
        "rapid_quant"
    )
    assert (
        command.value.preflight_argv[command.value.preflight_argv.index("-profile") + 1]
        == "rapid_quant"
    )
    assert command.value.argv.count("-C") == 1
    assert command.value.preflight_argv.count("-C") == 1
    assert "-resume" not in command.value.argv


def test_rapid_quant_mode_is_bound_fail_closed_across_builders(
    tmp_path: Path,
    composed_runtime,
    admitted_execution_implementation,
):
    binding, _ = composed_runtime
    inputs = _inputs(tmp_path)
    workspace = (tmp_path / "workspace").resolve()
    standard = BulkRnaSeqWorkflowAdapter(execution=binding)
    rapid = BulkRnaSeqRapidQuantQualificationAdapter(execution=binding)
    standard_plan = standard.plan_workspace(inputs, workspace).value
    rapid_plan = rapid.plan_workspace(inputs, workspace).value

    assert standard.build_command(rapid_plan, workspace).errors[0].code == (
        "BULK_RNASEQ_COMMAND_IDENTITY_MISMATCH"
    )
    assert rapid.build_command(standard_plan, workspace).errors[0].code == (
        "BULK_RNASEQ_COMMAND_IDENTITY_MISMATCH"
    )

    identity = json.loads(_file_bytes(rapid_plan, "config/execution-identity.json"))
    identity["execution_mode"] = "standard-v1"
    replacement = (
        json.dumps(identity, sort_keys=True, separators=(",", ":")) + "\n"
    ).encode()
    tampered = WorkspacePlan(
        directories=rapid_plan.directories,
        files=tuple(
            (
                path,
                replacement if path == "config/execution-identity.json" else value,
            )
            for path, value in rapid_plan.files
        ),
    )

    assert rapid.build_command(tampered, workspace).errors[0].code == (
        "BULK_RNASEQ_COMMAND_IDENTITY_MISMATCH"
    )


def test_rapid_quant_uses_the_complete_runtime_admission_boundary(
    tmp_path: Path,
    composed_runtime,
    admitted_execution_implementation,
    monkeypatch: pytest.MonkeyPatch,
):
    binding, verified = composed_runtime
    observed: list[RuntimeAssetAdmission] = []

    def acquire(execution_binding: BulkRnaSeqExecutionBinding):
        observed.append(execution_binding.runtime_admission)
        return Result.success(verified)

    monkeypatch.setattr(execution_module, "_acquire_runtime_assets", acquire)
    adapter = BulkRnaSeqRapidQuantQualificationAdapter(execution=binding)
    workspace = (tmp_path / "workspace").resolve()

    assert adapter.capture_build_identity().is_success
    plan = adapter.plan_workspace(_inputs(tmp_path), workspace)
    assert plan.is_success
    assert adapter.build_command(plan.value, workspace).is_success
    assert observed == [binding.runtime_admission] * 3


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
        "schema_version": "1.1.0",
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
    expected_isolation_prefix = (
        str(verified.network_isolation_executable),
        "--user",
        "--map-current-user",
        "--net",
        "--",
    )
    assert command.argv[:5] == expected_isolation_prefix
    assert command.argv[5] == str(verified.nextflow_executable)
    assert command.preflight_argv[0] == command.argv[0]
    assert command.preflight_argv[:6] == command.argv[:6]
    assert "--fork" not in command.argv
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
    assert "-profile" not in command.argv
    assert "-profile" not in command.preflight_argv
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
    assert config.count("docker.runOptions = ") == 1
    assert (
        "docker.runOptions = '--pull=never --network=none "
        f"--user={binding.container_uid}:{binding.container_gid} "
        f"--label={MANAGED_CONTAINER_SCOPE_LABEL}="
    ) in config
    assert "process.stageInMode = 'copy'" in config
    assert (
        "    resourceLimits = [\n"
        "        cpus: 4,\n"
        "        memory: '24.GB',\n"
        "        time: '2.h'\n"
        "    ]"
    ) in config
    assert "wave.enabled = false" in config
    assert "shifter.enabled = false" in config
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


def test_network_isolation_identity_changes_build_and_cache_identity(
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
            "network_isolation_executable_sha256": "0" * 64,
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
        BulkRnaSeqExecutionBinding(
            assets=assets,
            transcriptome=_transcriptome_binding(tmp_path),
            resume_enabled=True,
        )


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
