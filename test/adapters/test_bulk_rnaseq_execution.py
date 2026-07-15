"""Workspace and command tests for the offline bulk RNA-seq runtime binding."""

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
from encode_pipeline.adapters.bulk_rnaseq.reference_closure import (
    REFERENCE_INDEX_MANIFEST,
)
from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import (
    VerifiedContainerAsset,
    VerifiedRuntimeAssets,
)
from encode_pipeline.platform.adapters import WorkflowInputs, WorkspacePlan
from encode_pipeline.platform.managed_containers import (
    managed_container_endpoint_identity,
    managed_container_scope,
)
from encode_pipeline.platform.results import Result
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
        fastq_1 = tmp_path / f"inputs/S{index}.R1.fastq.gz"
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
            fastq_2 = tmp_path / f"inputs/S{index}.R2.fastq.gz"
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
        plugin_root=root / "plugins",
        plugin_archive=root / "plugins/nf-schema-2.5.1.zip",
        plugin_meta=root / "plugins/nf-schema-2.5.1-meta.json",
        plugin_tree=root / "plugins/nf-schema-2.5.1",
        container_lock=root / "containers/availability-lock.json",
        containers=containers,
        source_tree_sha256="5" * 64,
        runtime_identity_sha256="8" * 64,
        nextflow_sha256="9" * 64,
        plugin_archive_sha256="a" * 64,
        plugin_tree_sha256="6" * 64,
        container_inventory_sha256="b" * 64,
        container_lock_sha256="7" * 64,
    )
    monkeypatch.setattr(
        execution_module,
        "verify_runtime_assets",
        lambda _binding: Result.success(verified),
    )
    return binding, verified


def _file_bytes(plan, path: str) -> bytes:
    return dict(plan.files)[path]


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


def test_workspace_serializes_se_pe_and_mixed_layout_deterministically(
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
    assert _file_bytes(first.value, "config/samplesheet.csv").decode().splitlines() == [
        "sample,fastq_1,fastq_2,strandedness",
        f"S1,{tmp_path}/inputs/S1.R1.fastq.gz,,auto",
        f"S2,{tmp_path}/inputs/S2.R1.fastq.gz,{tmp_path}/inputs/S2.R2.fastq.gz,auto",
    ]
    params = json.loads(_file_bytes(first.value, "config/params.json"))
    assert params["input"] == str(workspace / "config/samplesheet.csv")
    assert params["outdir"] == str(workspace / "results")
    assert params["custom_config_base"] == ""
    assert params["igenomes_ignore"] is True
    assert params["monochrome_logs"] is True


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
    assert command.argv[command.argv.index("run") + 1] == str(verified.source_tree)
    assert command.argv[command.argv.index("-profile") + 1] == "docker"
    assert "-resume" not in command.argv
    assert not any("nf-core/rnaseq" == token for token in command.argv)
    config = _file_bytes(plan, "config/platform.nextflow.config").decode()
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
        "verify_runtime_assets",
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
        "verify_runtime_assets",
        lambda _binding: Result.success(changed_audit),
    )
    third = adapter.capture_build_identity()

    assert third.is_success
    assert len({first.value.digest, second.value.digest, third.value.digest}) == 3


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
