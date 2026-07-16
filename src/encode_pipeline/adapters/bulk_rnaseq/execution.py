"""Offline workspace and command contract for pinned nf-core/rnaseq."""

from __future__ import annotations

from collections.abc import Mapping
import csv
from dataclasses import dataclass, field
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import re
from typing import Any

from encode_pipeline.adapters.bulk_rnaseq.execution_identity import (
    VerifiedExecutionImplementation,
    verify_execution_implementation,
)
from encode_pipeline.adapters.bulk_rnaseq.reference_closure import (
    ReferenceClosure,
    verify_reference_closure,
)
from encode_pipeline.adapters.bulk_rnaseq.resource_closure import (
    DEFAULT_RESOURCE_CLOSURE_POLICY,
    SORTMERNA_NO_PREBUILT_INDEX_STRATEGY,
    ResourceClosurePolicy,
    RiboDatabaseClosure,
    SortMeRnaIndexClosure,
    safe_regular_file_identity,
    verify_ribo_database_manifest,
    verify_sortmerna_index,
)
from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import (
    CONTAINER_DEFAULT_CONFIG_DENIED_SELECTORS,
    CONTAINER_RESERVED_DEFAULT_DENY_LABEL,
    DOCKER_REQUIRED_RUN_OPTIONS,
    NEXTFLOW_VERSION,
    NF_SCHEMA_VERSION,
    SOURCE_MANIFEST_SHA256,
    RuntimeAssetBinding,
    RuntimeAssetAdmission,
    RuntimeAssetDoctorReport,
    VerifiedRuntimeAssets,
)
from encode_pipeline.adapters.bulk_rnaseq.upstream import (
    NFCORE_RNASEQ_COMMIT,
    UPSTREAM_PARAMETER_SCHEMA_SHA256,
    UPSTREAM_SAMPLESHEET_SCHEMA_SHA256,
)
from encode_pipeline.adapters.bulk_rnaseq.validation import (
    validate_bulk_rnaseq_inputs,
)
from encode_pipeline.platform.adapters import CommandSpec, WorkflowInputs, WorkspacePlan
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.managed_containers import (
    MANAGED_CONTAINER_SCOPE_LABEL,
    managed_container_endpoint_identity,
    managed_container_scope,
)
from encode_pipeline.platform.results import Issue, Result


BUILD_IDENTITY_SCHEME = "sha256-bulk-rnaseq-runtime-v1"
LOGICAL_ENTRYPOINT = "main.nf"
WORKSPACE_SCHEMA_VERSION = "1.0.0"
_UNAVAILABLE_CONTAINER_IMAGE = f"sha256:{'0' * 64}"
_WORKSPACE_DIRECTORIES = (
    "config",
    "engine/cache",
    "engine/home",
    "engine/launch",
    "engine/nxf-home",
    "engine/tmp",
    "engine/work",
    "logs",
    "reports",
    "results",
)

_SAFE_WORKSPACE = re.compile(r"^/[A-Za-z0-9._/+@%=-]+$")
_SAFE_PROCESS = re.compile(r"^[A-Z][A-Z0-9_]*(?::[A-Z][A-Z0-9_]*)*$")
_SAFE_AUDITED_SELECTOR = re.compile(r"^\.\*[A-Z][A-Z0-9_]*$")


@dataclass(frozen=True)
class BulkRnaSeqExecutionBinding:
    """Deployment-owned execution settings; no field comes from WorkflowInputs."""

    assets: RuntimeAssetBinding
    container_uid: int = field(default_factory=os.getuid)
    container_gid: int = field(default_factory=os.getgid)
    resume_enabled: bool = False
    resource_policy: ResourceClosurePolicy = DEFAULT_RESOURCE_CLOSURE_POLICY
    runtime_admission: RuntimeAssetAdmission = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.assets, RuntimeAssetBinding):
            raise ValueError("assets must be a RuntimeAssetBinding")
        for name in ("container_uid", "container_gid"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"{name} must be a non-negative integer")
        if self.resume_enabled is not False:
            raise ValueError(
                "resume requires a durable attempt/session lifecycle and is not enabled"
            )
        if not isinstance(self.resource_policy, ResourceClosurePolicy):
            raise ValueError("resource_policy must be a ResourceClosurePolicy")
        object.__setattr__(
            self,
            "runtime_admission",
            RuntimeAssetAdmission(self.assets),
        )


@dataclass(frozen=True)
class VerifiedInputClosure:
    """Private resource evidence used to construct one deterministic workspace."""

    reference: ReferenceClosure
    fastq_files: tuple[tuple[str, int, str], ...]
    ribo_database: RiboDatabaseClosure | None
    sortmerna_index: SortMeRnaIndexClosure | None
    sortmerna_index_build_strategy: str | None
    identity_sha256: str


def plan_bulk_rnaseq_workspace(
    inputs: WorkflowInputs,
    workspace: str | Path,
    *,
    binding: BulkRnaSeqExecutionBinding,
    adapter_version: str,
    adapter_variant: str,
) -> Result[WorkspacePlan]:
    """Verify resources and return deterministic, server-owned workspace bytes."""
    validation = validate_bulk_rnaseq_inputs(inputs)
    if validation.is_failure:
        return Result.failure(validation.issues)
    try:
        workspace_path = _safe_workspace_path(workspace)
    except ValueError:
        return _failure("BULK_RNASEQ_WORKSPACE_INVALID", "workspace")

    implementation_result = verify_execution_implementation()
    if implementation_result.is_failure:
        return Result.failure(implementation_result.issues)
    assets_result = _acquire_runtime_assets(binding)
    if assets_result.is_failure:
        return _failure("BULK_RNASEQ_RUNTIME_UNAVAILABLE", "runtime")
    normalized = validation.value
    input_result = _verify_input_closure(
        inputs,
        binding=binding,
        assets=assets_result.value,
        normalized=normalized,
    )
    if input_result.is_failure:
        return input_result

    verified_inputs = input_result.value
    implementation = implementation_result.value
    build_identity_sha256 = _runtime_build_digest(
        assets_result.value,
        adapter_version=adapter_version,
        adapter_variant=adapter_variant,
        implementation=implementation,
    )
    workspace_identity_sha256 = managed_container_scope(workspace_path)
    params = _runtime_params(
        normalized,
        inputs=inputs,
        closure=verified_inputs,
        workspace=workspace_path,
    )
    samplesheet = _samplesheet_bytes(normalized["samples"])
    platform_config = _platform_config_bytes(
        assets_result.value,
        workspace=workspace_path,
        binding=binding,
        workspace_identity=workspace_identity_sha256,
        build_sortmerna_index=(
            verified_inputs.sortmerna_index_build_strategy is not None
        ),
    )
    cache_identity = _cache_identity_document(
        normalized=normalized,
        inputs=verified_inputs,
        adapter_version=adapter_version,
        implementation=implementation,
        workflow_build_sha256=build_identity_sha256,
    )
    files: list[tuple[str, bytes]] = [
        ("config/samplesheet.csv", samplesheet),
        ("config/params.json", _json_bytes(params)),
        ("config/platform.nextflow.config", platform_config),
        ("engine/cache-identity.json", _json_bytes(cache_identity)),
    ]
    if verified_inputs.ribo_database is not None:
        files.append(
            (
                "config/ribo-database-manifest.txt",
                _ribo_manifest_bytes(verified_inputs.ribo_database),
            )
        )
    workspace_contract_sha256 = _workspace_contract_digest(
        directories=_WORKSPACE_DIRECTORIES,
        files=files,
    )
    execution_identity = {
        "schema_version": WORKSPACE_SCHEMA_VERSION,
        "workflow_id": "bulk-rnaseq",
        "adapter_version": adapter_version,
        "build_identity_sha256": build_identity_sha256,
        "input_identity_sha256": verified_inputs.identity_sha256,
        "ribo_database_closure_sha256": (
            verified_inputs.ribo_database.identity_sha256
            if verified_inputs.ribo_database is not None
            else None
        ),
        "sortmerna_index_build_strategy": (
            verified_inputs.sortmerna_index_build_strategy
        ),
        "cache_identity_sha256": cache_identity["identity_sha256"],
        "execution_implementation_manifest_sha256": (implementation.manifest_sha256),
        "execution_implementation_aggregate_sha256": (implementation.aggregate_sha256),
        "container_process_audit_sha256": (
            assets_result.value.container_process_audit_sha256
        ),
        "workspace_identity_sha256": workspace_identity_sha256,
        "workspace_contract_sha256": workspace_contract_sha256,
        "resume_enabled": False,
    }
    files.append(("config/execution-identity.json", _json_bytes(execution_identity)))

    return Result.success(
        WorkspacePlan(
            directories=_WORKSPACE_DIRECTORIES,
            files=tuple(sorted(files)),
        )
    )


def build_bulk_rnaseq_command(
    plan: WorkspacePlan,
    workspace: str | Path,
    *,
    binding: BulkRnaSeqExecutionBinding,
    adapter_version: str,
    adapter_variant: str,
) -> Result[CommandSpec]:
    """Build fixed shell-free Nextflow argv from one adapter-owned plan."""
    if not isinstance(plan, WorkspacePlan):
        return _failure("BULK_RNASEQ_COMMAND_INVALID", "workspace_plan")
    try:
        workspace_path = _safe_workspace_path(workspace)
        params = _planned_json(plan, "config/params.json")
        identity = _planned_json(plan, "config/execution-identity.json")
        cache_identity = _planned_json(plan, "engine/cache-identity.json")
        submitted_resource_paths = _planned_resource_paths(plan)
        workspace_contract_sha256 = _workspace_contract_digest(
            directories=plan.directories,
            files=tuple(
                (path, value)
                for path, value in plan.files
                if path != "config/execution-identity.json"
            ),
        )
    except (TypeError, ValueError, json.JSONDecodeError):
        return _failure("BULK_RNASEQ_COMMAND_INVALID", "workspace_plan")
    implementation_result = verify_execution_implementation()
    if implementation_result.is_failure:
        return Result.failure(implementation_result.issues)
    assets_result = _acquire_runtime_assets(binding)
    if assets_result.is_failure:
        return _failure("BULK_RNASEQ_RUNTIME_UNAVAILABLE", "runtime")
    assets = assets_result.value
    implementation = implementation_result.value
    expected_build = _runtime_build_digest(
        assets,
        adapter_version=adapter_version,
        adapter_variant=adapter_variant,
        implementation=implementation,
    )
    if (
        not _execution_identity_shape_is_valid(identity)
        or identity.get("schema_version") != WORKSPACE_SCHEMA_VERSION
        or identity.get("workflow_id") != "bulk-rnaseq"
        or identity.get("adapter_version") != adapter_version
        or identity.get("build_identity_sha256") != expected_build
        or identity.get("execution_implementation_manifest_sha256")
        != implementation.manifest_sha256
        or identity.get("execution_implementation_aggregate_sha256")
        != implementation.aggregate_sha256
        or identity.get("container_process_audit_sha256")
        != assets.container_process_audit_sha256
        or identity.get("workspace_identity_sha256")
        != managed_container_scope(workspace_path)
        or identity.get("workspace_contract_sha256") != workspace_contract_sha256
        or identity.get("cache_identity_sha256")
        != cache_identity.get("identity_sha256")
        or identity.get("input_identity_sha256")
        != cache_identity.get("input_closure_sha256")
        or identity.get("ribo_database_closure_sha256")
        != cache_identity.get("ribo_database_closure_sha256")
        or identity.get("sortmerna_index_build_strategy")
        != cache_identity.get("sortmerna_index_build_strategy")
        or not _cache_identity_is_valid(
            cache_identity,
            adapter_version=adapter_version,
            implementation=implementation,
            expected_build=expected_build,
        )
    ):
        return _failure("BULK_RNASEQ_COMMAND_IDENTITY_MISMATCH", "workspace_plan")
    if identity.get("resume_enabled") is not False:
        return _failure("BULK_RNASEQ_RESUME_UNAVAILABLE", "workspace_plan")

    executable = str(assets.nextflow_executable)
    config_path = workspace_path / "config/platform.nextflow.config"
    params_path = workspace_path / "config/params.json"
    source_path = assets.source_tree
    launch_path = workspace_path / "engine/launch"
    log_path = workspace_path / "logs/nextflow.log"
    preflight_log = workspace_path / "logs/nextflow-preflight.log"
    input_identity = identity.get("input_identity_sha256")
    if not isinstance(input_identity, str) or len(input_identity) != 64:
        return _failure("BULK_RNASEQ_COMMAND_INVALID", "workspace_plan")

    argv = (
        executable,
        "-log",
        str(log_path),
        "-C",
        str(config_path),
        "run",
        str(source_path),
        "-profile",
        "docker",
        "-offline",
        "-params-file",
        str(params_path),
        "-work-dir",
        str(workspace_path / "engine/work"),
        "-name",
        f"helixweave-{input_identity[:16]}",
    )
    preflight_argv = (
        executable,
        "-log",
        str(preflight_log),
        "-C",
        str(config_path),
        "config",
        str(source_path),
        "-profile",
        "docker",
    )
    env = {
        "CLASSPATH": "",
        "CONDA_DEFAULT_ENV": "",
        "CONDA_EXE": "",
        "CONDA_PREFIX": "",
        "CONDA_PYTHON_EXE": "",
        "CONDA_SHLVL": "0",
        "HOME": str(workspace_path / "engine/home"),
        "JAVA_CMD": str(assets.java_executable),
        "JAVA_HOME": str(assets.jdk_tree),
        "LD_LIBRARY_PATH": "",
        "MAMBA_EXE": "",
        "MAMBA_ROOT_PREFIX": "",
        "NXF_ANSI_LOG": "false",
        "NXF_CACHE_DIR": str(workspace_path / "engine/cache"),
        "NXF_DISABLE_CHECK_LATEST": "true",
        "NXF_HOME": str(workspace_path / "engine/nxf-home"),
        "NXF_JAVA_HOME": str(assets.jdk_tree),
        "NXF_OFFLINE": "true",
        "NXF_OPTS": "",
        "NXF_PLUGINS_DIR": str(assets.plugin_root),
        "NXF_TEMP": str(workspace_path / "engine/tmp"),
        "PATH": f"{assets.jdk_tree / 'bin'}:/usr/bin:/bin",
        "_CE_CONDA": "",
        "_CE_M": "",
        "_CONDA_EXE": "",
        "_CONDA_ROOT": "",
    }
    redactions = {
        str(workspace_path),
        str(assets.root),
        str(assets.source_tree),
        str(assets.nextflow_executable),
        str(assets.jdk_archive),
        str(assets.jdk_tree),
        str(assets.java_executable),
        str(assets.plugin_tree),
        str(binding.assets.docker_executable),
        str(binding.assets.docker_socket),
        f"unix://{binding.assets.docker_socket}",
    }
    redactions.update(_absolute_json_strings(params))
    redactions.update(submitted_resource_paths)
    return Result.success(
        CommandSpec(
            argv=argv,
            cwd=str(launch_path),
            env=env,
            preflight_argv=preflight_argv,
            preflight_kind="configuration",
            preflight_managed_logs=(("nextflow_preflight", str(preflight_log)),),
            execution_managed_logs=(("nextflow", str(log_path)),),
            redaction_values=tuple(sorted(redactions)),
            managed_container_scope=managed_container_scope(workspace_path),
            managed_container_endpoint_identity=managed_container_endpoint_identity(
                binding.assets.docker_executable,
                binding.assets.docker_socket,
            ),
        )
    )


def capture_bulk_rnaseq_build_identity(
    *,
    binding: BulkRnaSeqExecutionBinding,
    adapter_version: str,
    adapter_variant: str,
) -> Result[WorkflowBuildIdentity]:
    """Capture source, engine, plugin, container and adapter contract identity."""
    implementation_result = verify_execution_implementation()
    if implementation_result.is_failure:
        return Result.failure(implementation_result.issues)
    assets_result = _acquire_runtime_assets(binding)
    if assets_result.is_failure:
        return _failure("BULK_RNASEQ_RUNTIME_UNAVAILABLE", "runtime")
    try:
        identity = WorkflowBuildIdentity(
            workflow_id="bulk-rnaseq",
            adapter_version=adapter_version,
            scheme=BUILD_IDENTITY_SCHEME,
            logical_entrypoint=LOGICAL_ENTRYPOINT,
            digest=_runtime_build_digest(
                assets_result.value,
                adapter_version=adapter_version,
                adapter_variant=adapter_variant,
                implementation=implementation_result.value,
            ),
            captured_at=datetime.now(timezone.utc),
        )
    except ValueError:
        return _failure("BULK_RNASEQ_RUNTIME_UNAVAILABLE", "runtime")
    return Result.success(identity)


def doctor_bulk_rnaseq_runtime(binding: BulkRnaSeqExecutionBinding):
    """Return the redacted runtime-asset doctor report."""
    implementation = verify_execution_implementation()
    assets = binding.runtime_admission.doctor()
    return RuntimeAssetDoctorReport(
        ready=implementation.is_success and assets.ready,
        issues=(*implementation.issues, *assets.issues),
    )


def _acquire_runtime_assets(
    binding: BulkRnaSeqExecutionBinding,
) -> Result[VerifiedRuntimeAssets]:
    return binding.runtime_admission.acquire()


def _verify_input_closure(
    inputs: WorkflowInputs,
    *,
    binding: BulkRnaSeqExecutionBinding,
    assets: VerifiedRuntimeAssets,
    normalized: Mapping[str, Any],
) -> Result[VerifiedInputClosure]:
    standard = inputs.config["standard"]
    reference_result = verify_reference_closure(
        standard["reference"],
        producer_images={item.process: item.image for item in assets.containers},
        index_build_parameters=normalized["nfcore_params"],
        policy=binding.resource_policy,
    )
    if reference_result.is_failure:
        return Result.failure(reference_result.issues)

    fastqs: list[tuple[str, int, str]] = []
    assert isinstance(inputs.samples, list)
    for row in inputs.samples:
        for field_name in ("fastq_1", "fastq_2"):
            path = row.get(field_name)
            if not path:
                continue
            result = safe_regular_file_identity(
                path,
                policy=binding.resource_policy,
            )
            if result.is_failure:
                return _failure("BULK_RNASEQ_FASTQ_INVALID", "samples")
            fastqs.append(
                (str(result.value.path), result.value.size_bytes, result.value.sha256)
            )

    ribo_database = None
    sortmerna_index = None
    sortmerna_index_build_strategy = None
    ribo = standard.get("ribosomal_rna_removal", {"enabled": False})
    if ribo["enabled"]:
        manifest = ribo["database_manifest"]
        result = verify_ribo_database_manifest(
            manifest["path"],
            expected_manifest_sha256=manifest["identity_sha256"],
            policy=binding.resource_policy,
        )
        if result.is_failure:
            return Result.failure(result.issues)
        ribo_database = result.value
        if ribo["tool"] == "sortmerna":
            index = ribo.get("sortmerna_index")
            if index is None:
                route = ribo_database.no_prebuilt_sortmerna_index
                if (
                    len(ribo_database.files) != 1
                    or not route.deterministic_database_inputs
                    or not route.deterministic_composition_accepted
                    or route.strategy != SORTMERNA_NO_PREBUILT_INDEX_STRATEGY
                ):
                    return _failure(
                        "BULK_RNASEQ_SORTMERNA_INDEX_BUILD_UNQUALIFIED",
                        "config.standard.ribosomal_rna_removal.sortmerna_index",
                    )
                sortmerna_index_build_strategy = route.strategy
            else:
                index_result = verify_sortmerna_index(
                    index["path"],
                    expected_index_sha256=index["identity_sha256"],
                    database_closure=ribo_database,
                    policy=binding.resource_policy,
                )
                if index_result.is_failure:
                    return Result.failure(index_result.issues)
                sortmerna_index = index_result.value

    digest = hashlib.sha256()
    _frame(digest, b"bulk-rnaseq-input-closure-v1")
    _frame(digest, bytes.fromhex(reference_result.value.identity_sha256))
    for path, size, file_sha256 in sorted(fastqs):
        _frame(digest, path.encode())
        _frame(digest, str(size).encode())
        _frame(digest, bytes.fromhex(file_sha256))
    if ribo_database is not None:
        _frame(digest, bytes.fromhex(ribo_database.identity_sha256))
    if sortmerna_index is not None:
        _frame(digest, bytes.fromhex(sortmerna_index.identity_sha256))
    if sortmerna_index_build_strategy is not None:
        _frame(digest, b"sortmerna-index-build-strategy")
        _frame(digest, sortmerna_index_build_strategy.encode())
    return Result.success(
        VerifiedInputClosure(
            reference=reference_result.value,
            fastq_files=tuple(sorted(fastqs)),
            ribo_database=ribo_database,
            sortmerna_index=sortmerna_index,
            sortmerna_index_build_strategy=sortmerna_index_build_strategy,
            identity_sha256=digest.hexdigest(),
        )
    )


def _runtime_params(
    normalized: Mapping[str, Any],
    *,
    inputs: WorkflowInputs,
    closure: VerifiedInputClosure,
    workspace: Path,
) -> dict[str, Any]:
    params = dict(normalized["nfcore_params"])
    params.update(
        {
            "custom_config_base": "",
            "input": str(workspace / "config/samplesheet.csv"),
            "monochrome_logs": True,
            "outdir": str(workspace / "results"),
            "trace_report_suffix": "helixweave",
            "validate_params": True,
        }
    )
    ribo = inputs.config["standard"].get("ribosomal_rna_removal", {"enabled": False})
    if ribo["enabled"]:
        params["ribo_database_manifest"] = str(
            workspace / "config/ribo-database-manifest.txt"
        )
        if closure.sortmerna_index is not None:
            params["sortmerna_index"] = str(closure.sortmerna_index.path)
    return {name: params[name] for name in sorted(params)}


def _samplesheet_bytes(samples: list[Mapping[str, str]]) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(("sample", "fastq_1", "fastq_2", "strandedness", "seq_platform"))
    for row in sorted(
        samples,
        key=lambda item: (item["sample"], item["library"], item["lane"]),
    ):
        writer.writerow(
            (
                row["sample"],
                row["fastq_1"],
                row["fastq_2"],
                row["strandedness"],
                row["seq_platform"],
            )
        )
    return output.getvalue().encode()


def _platform_config_bytes(
    assets: VerifiedRuntimeAssets,
    *,
    workspace: Path,
    binding: BulkRnaSeqExecutionBinding,
    workspace_identity: str,
    build_sortmerna_index: bool,
) -> bytes:
    run_options = (
        *DOCKER_REQUIRED_RUN_OPTIONS,
        f"--user={binding.container_uid}:{binding.container_gid}",
        f"--label={MANAGED_CONTAINER_SCOPE_LABEL}={workspace_identity}",
    )
    lines = [
        "// Generated by HelixWeave; user configuration is never included.",
        (
            "includeConfig "
            f"{_groovy_string(str(assets.source_tree / 'nextflow.config'))}"
        ),
        "process.executor = 'local'",
        "docker.enabled = true",
        "docker.registry = ''",
        "docker.remove = true",
        f"docker.runOptions = {_groovy_string(' '.join(run_options))}",
        "apptainer.enabled = false",
        "charliecloud.enabled = false",
        "conda.enabled = false",
        "podman.enabled = false",
        "singularity.enabled = false",
        "spack.enabled = false",
        "wave.enabled = false",
        "fusion.enabled = false",
        "tower.enabled = false",
        "timeline.enabled = true",
        'timeline.file = "${launchDir}/../../reports/timeline.html"',
        "report.enabled = true",
        'report.file = "${launchDir}/../../reports/report.html"',
        "trace.enabled = true",
        'trace.file = "${launchDir}/../../reports/trace.txt"',
        "dag.enabled = true",
        'dag.file = "${launchDir}/../../reports/dag.html"',
        "process.stageInMode = 'copy'",
        "process {",
        (f"    withLabel: '!{CONTAINER_RESERVED_DEFAULT_DENY_LABEL}' {{"),
        f"        container = '{_UNAVAILABLE_CONTAINER_IMAGE}'",
        "    }",
    ]
    for container in assets.containers:
        if _SAFE_PROCESS.fullmatch(container.process) is None:
            raise ValueError("verified container process name is invalid")
        lines.extend(
            (
                f"    withName: {_groovy_string(container.process)} {{",
                f"        container = {_groovy_string(container.runtime_image)}",
                "    }",
            )
        )
    # Nextflow gives an alias-specific selector higher priority than an
    # original-process selector (and than the negative-label default deny).
    # Repeat every immutable, audited default-config selector after the
    # verified allowlist so an unsupported alias cannot restore its upstream
    # network image coordinate.
    for selector in CONTAINER_DEFAULT_CONFIG_DENIED_SELECTORS:
        if _SAFE_AUDITED_SELECTOR.fullmatch(selector) is None:
            raise ValueError("audited container selector is invalid")
        lines.extend(
            (
                f"    withName: {_groovy_string(selector)} {{",
                f"        container = '{_UNAVAILABLE_CONTAINER_IMAGE}'",
                "    }",
            )
        )
    if build_sortmerna_index:
        sortmerna_container = _container_execution_coordinate(
            assets,
            process="SORTMERNA",
        )
        lines.extend(
            (
                "    withName: '.*:PREPARE_GENOME:SORTMERNA_INDEX' {",
                "        ext.when = false",
                "    }",
                "    withName: '.*:FASTQ_REMOVE_RRNA:SORTMERNA_INDEX' {",
                f"        container = {_groovy_string(sortmerna_container)}",
                "        ext.args = '--index 1'",
                "        cpus = 1",
                "        maxRetries = 0",
                "        errorStrategy = 'terminate'",
                "    }",
            )
        )
    lines.extend(("}", ""))
    return "\n".join(lines).encode()


def _cache_identity_document(
    *,
    normalized: Mapping[str, Any],
    inputs: VerifiedInputClosure,
    adapter_version: str,
    implementation: VerifiedExecutionImplementation,
    workflow_build_sha256: str,
) -> dict[str, Any]:
    coordinates = {
        "schema_version": WORKSPACE_SCHEMA_VERSION,
        "adapter_version": adapter_version,
        "workflow_build_sha256": workflow_build_sha256,
        "execution_implementation_manifest_sha256": (implementation.manifest_sha256),
        "execution_implementation_aggregate_sha256": (implementation.aggregate_sha256),
        "input_closure_sha256": inputs.identity_sha256,
        "ribo_database_closure_sha256": (
            inputs.ribo_database.identity_sha256
            if inputs.ribo_database is not None
            else None
        ),
        "sortmerna_index_build_strategy": inputs.sortmerna_index_build_strategy,
        "normalized_inputs_sha256": _canonical_digest(normalized),
        "nextflow_version": NEXTFLOW_VERSION,
        "resume_scope": "single-run-workspace",
    }
    return {**coordinates, "identity_sha256": _canonical_digest(coordinates)}


def _runtime_build_digest(
    assets: VerifiedRuntimeAssets,
    *,
    adapter_version: str,
    adapter_variant: str,
    implementation: VerifiedExecutionImplementation,
) -> str:
    if not isinstance(adapter_version, str) or not adapter_version.strip():
        raise ValueError("adapter_version must be non-empty")
    if (
        not isinstance(adapter_variant, str)
        or re.fullmatch(r"[a-z][a-z0-9-]{0,63}", adapter_variant) is None
    ):
        raise ValueError("adapter_variant is invalid")
    digest = hashlib.sha256()
    for value in (
        BUILD_IDENTITY_SCHEME.encode(),
        b"adapter-version",
        adapter_version.strip().encode(),
        b"adapter-variant",
        adapter_variant.encode(),
        b"execution-implementation-manifest",
        bytes.fromhex(implementation.manifest_sha256),
        b"execution-implementation-aggregate",
        bytes.fromhex(implementation.aggregate_sha256),
        b"sortmerna-no-prebuilt-index-strategy",
        SORTMERNA_NO_PREBUILT_INDEX_STRATEGY.encode(),
    ):
        _frame(digest, value)
    for value in (
        NFCORE_RNASEQ_COMMIT,
        assets.source_tree_sha256,
        SOURCE_MANIFEST_SHA256,
        NEXTFLOW_VERSION,
        assets.nextflow_sha256,
        assets.jdk_archive_sha256,
        assets.jdk_tree_sha256,
        assets.java_executable_sha256,
        NF_SCHEMA_VERSION,
        assets.plugin_archive_sha256,
        assets.plugin_tree_sha256,
        assets.container_inventory_sha256,
        assets.container_process_audit_sha256,
        assets.container_lock_sha256,
        assets.runtime_identity_sha256,
        UPSTREAM_PARAMETER_SCHEMA_SHA256,
        UPSTREAM_SAMPLESHEET_SCHEMA_SHA256,
    ):
        _frame(digest, value.encode())
    return digest.hexdigest()


def _container_execution_coordinate(
    assets: VerifiedRuntimeAssets,
    *,
    process: str,
) -> str:
    """Return one verified runtime coordinate for an adapter-owned selector."""
    matches = [item for item in assets.containers if item.process == process]
    if len(matches) != 1:
        raise ValueError("verified process container is missing or duplicated")
    coordinate = matches[0].runtime_image
    if not isinstance(coordinate, str) or not coordinate:
        raise ValueError("verified process container coordinate is invalid")
    return coordinate


def _ribo_manifest_bytes(closure: RiboDatabaseClosure) -> bytes:
    return ("\n".join(str(item.path) for item in closure.files) + "\n").encode()


def _workspace_contract_digest(
    *,
    directories: tuple[str, ...],
    files: list[tuple[str, bytes]] | tuple[tuple[str, bytes], ...],
) -> str:
    expected_files = {
        "config/samplesheet.csv",
        "config/params.json",
        "config/platform.nextflow.config",
        "engine/cache-identity.json",
    }
    paths = tuple(path for path, _value in files)
    if (
        directories != _WORKSPACE_DIRECTORIES
        or len(paths) != len(set(paths))
        or set(paths)
        not in (
            expected_files,
            expected_files | {"config/ribo-database-manifest.txt"},
        )
    ):
        raise ValueError("workspace contract file set is invalid")
    digest = hashlib.sha256()
    _frame(digest, b"sha256-framed-bulk-rnaseq-workspace-contract-v1")
    for directory in directories:
        _frame(digest, directory.encode())
    for path, value in sorted(files):
        if not isinstance(value, bytes):
            raise ValueError("workspace contract file bytes are invalid")
        _frame(digest, path.encode())
        _frame(digest, value)
    return digest.hexdigest()


def _cache_identity_is_valid(
    value: Mapping[str, Any],
    *,
    adapter_version: str,
    implementation: VerifiedExecutionImplementation,
    expected_build: str,
) -> bool:
    expected_keys = {
        "schema_version",
        "adapter_version",
        "workflow_build_sha256",
        "execution_implementation_manifest_sha256",
        "execution_implementation_aggregate_sha256",
        "input_closure_sha256",
        "ribo_database_closure_sha256",
        "sortmerna_index_build_strategy",
        "normalized_inputs_sha256",
        "nextflow_version",
        "resume_scope",
        "identity_sha256",
    }
    if set(value) != expected_keys:
        return False
    coordinates = {
        name: item for name, item in value.items() if name != "identity_sha256"
    }
    return (
        value.get("schema_version") == WORKSPACE_SCHEMA_VERSION
        and value.get("adapter_version") == adapter_version
        and value.get("workflow_build_sha256") == expected_build
        and value.get("execution_implementation_manifest_sha256")
        == implementation.manifest_sha256
        and value.get("execution_implementation_aggregate_sha256")
        == implementation.aggregate_sha256
        and isinstance(value.get("input_closure_sha256"), str)
        and re.fullmatch(r"[0-9a-f]{64}", value["input_closure_sha256"]) is not None
        and isinstance(value.get("normalized_inputs_sha256"), str)
        and re.fullmatch(r"[0-9a-f]{64}", value["normalized_inputs_sha256"]) is not None
        and (
            value.get("ribo_database_closure_sha256") is None
            or isinstance(value.get("ribo_database_closure_sha256"), str)
            and re.fullmatch(r"[0-9a-f]{64}", value["ribo_database_closure_sha256"])
            is not None
        )
        and value.get("sortmerna_index_build_strategy")
        in {None, SORTMERNA_NO_PREBUILT_INDEX_STRATEGY}
        and value.get("nextflow_version") == NEXTFLOW_VERSION
        and value.get("resume_scope") == "single-run-workspace"
        and value.get("identity_sha256") == _canonical_digest(coordinates)
    )


def _execution_identity_shape_is_valid(value: Mapping[str, Any]) -> bool:
    expected_keys = {
        "schema_version",
        "workflow_id",
        "adapter_version",
        "build_identity_sha256",
        "input_identity_sha256",
        "ribo_database_closure_sha256",
        "sortmerna_index_build_strategy",
        "cache_identity_sha256",
        "execution_implementation_manifest_sha256",
        "execution_implementation_aggregate_sha256",
        "container_process_audit_sha256",
        "workspace_identity_sha256",
        "workspace_contract_sha256",
        "resume_enabled",
    }
    if set(value) != expected_keys:
        return False
    digest_names = (
        "build_identity_sha256",
        "input_identity_sha256",
        "cache_identity_sha256",
        "execution_implementation_manifest_sha256",
        "execution_implementation_aggregate_sha256",
        "container_process_audit_sha256",
        "workspace_identity_sha256",
        "workspace_contract_sha256",
    )
    return all(
        isinstance(value.get(name), str)
        and re.fullmatch(r"[0-9a-f]{64}", value[name]) is not None
        for name in digest_names
    )


def _planned_json(plan: WorkspacePlan, relative_path: str) -> dict[str, Any]:
    value = json.loads(_planned_bytes(plan, relative_path))
    if not isinstance(value, dict):
        raise ValueError("planned JSON must be an object")
    return value


def _planned_bytes(plan: WorkspacePlan, relative_path: str) -> bytes:
    matches = [content for path, content in plan.files if path == relative_path]
    if len(matches) != 1 or not isinstance(matches[0], bytes):
        raise ValueError("planned file missing, duplicated, or invalid")
    return matches[0]


def _planned_resource_paths(plan: WorkspacePlan) -> set[str]:
    """Recover submitted paths from adapter-owned bytes for output redaction."""
    rows = list(
        csv.reader(
            io.StringIO(_planned_bytes(plan, "config/samplesheet.csv").decode("utf-8")),
            strict=True,
        )
    )
    if not rows or rows[0] != [
        "sample",
        "fastq_1",
        "fastq_2",
        "strandedness",
        "seq_platform",
    ]:
        raise ValueError("planned samplesheet header is invalid")
    paths: set[str] = set()
    for row in rows[1:]:
        if len(row) != 5 or row[4] != "ILLUMINA":
            raise ValueError("planned samplesheet row is invalid")
        for value in row[1:3]:
            if value:
                if not value.startswith("/"):
                    raise ValueError("planned input path is invalid")
                paths.add(value)

    ribo_matches = [
        content
        for path, content in plan.files
        if path == "config/ribo-database-manifest.txt"
    ]
    if len(ribo_matches) > 1:
        raise ValueError("planned rRNA manifest is duplicated")
    if ribo_matches:
        for value in ribo_matches[0].decode("utf-8").splitlines():
            if not value.startswith("/"):
                raise ValueError("planned rRNA database path is invalid")
            paths.add(value)
    return paths


def _absolute_json_strings(value: object) -> set[str]:
    found: set[str] = set()
    if isinstance(value, str) and value.startswith("/"):
        found.add(value)
    elif isinstance(value, list):
        for item in value:
            found.update(_absolute_json_strings(item))
    elif isinstance(value, Mapping):
        for item in value.values():
            found.update(_absolute_json_strings(item))
    return found


def _safe_workspace_path(value: str | Path) -> Path:
    if not isinstance(value, (str, Path)) or isinstance(value, bool):
        raise ValueError("workspace must be a path")
    raw = os.fspath(value)
    path = Path(raw)
    if (
        not path.is_absolute()
        or _SAFE_WORKSPACE.fullmatch(raw) is None
        or any(part in {"", ".", ".."} for part in path.parts[1:])
    ):
        raise ValueError("workspace path is unsafe")
    return path


def _groovy_string(value: str) -> str:
    if any(character in value for character in ("'", "\\", "\n", "\r", "$")):
        raise ValueError("unsafe fixed configuration value")
    return f"'{value}'"


def _json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        + "\n"
    ).encode()


def _canonical_digest(value: object) -> str:
    return hashlib.sha256(_json_bytes(value)).hexdigest()


def _frame(digest: Any, value: bytes) -> None:
    digest.update(len(value).to_bytes(8, "big"))
    digest.update(value)


def _failure(code: str, path: str) -> Result[Any]:
    return Result.failure(
        [
            Issue(
                code=code,
                message="Bulk RNA-seq execution preparation failed.",
                severity="error",
                path=path,
                source="adapter",
            )
        ]
    )
