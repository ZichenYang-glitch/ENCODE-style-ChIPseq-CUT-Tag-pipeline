"""Source-owned qualification record for the Bulk RNA-seq execution boundary."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass
from enum import Enum
import hashlib
from importlib import resources
import json
from pathlib import Path
import re
from typing import Any

from encode_pipeline.adapters.bulk_rnaseq.authoring import SCHEMA_VERSION
from encode_pipeline.adapters.bulk_rnaseq.execution_identity import (
    EXECUTION_IMPLEMENTATION_MANIFEST_FILE,
    EXECUTION_IMPLEMENTATION_SCHEMA_VERSION,
    EXECUTION_IMPLEMENTATION_SCHEME,
    EXECUTION_PERSISTENCE_CONTRACT_PATH,
    ExecutionImplementationQualification,
    VerifiedExecutionImplementation,
)
from encode_pipeline.adapters.bulk_rnaseq.reference_closure import (
    REFERENCE_CLOSURE_SCHEME,
    REFERENCE_INDEX_MANIFEST,
    REFERENCE_INDEX_SCHEMA_VERSION,
)
from encode_pipeline.adapters.bulk_rnaseq.resource_closure import (
    DEFAULT_RESOURCE_CLOSURE_POLICY,
    RIBO_DATABASE_CLOSURE_SCHEME,
    SORTMERNA_INDEX_CLOSURE_SCHEME,
    SORTMERNA_INDEX_MANIFEST_SCHEMA_VERSION,
    SORTMERNA_NO_PREBUILT_INDEX_STRATEGY,
    SORTMERNA_VERSION,
    safe_regular_file_bytes,
)
from encode_pipeline.adapters.bulk_rnaseq.results_contract import (
    RESULTS_CONTRACT_FILE,
    RESULTS_CONTRACT_SHA256,
    RESULTS_CONTRACT_SIZE,
)
from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import (
    CONTAINER_CONFIG_ASSIGNMENTS_SHA256,
    CONTAINER_CONFIG_ASSIGNMENT_COUNT,
    CONTAINER_DEFAULT_CONFIG_ASSIGNMENTS_SHA256,
    CONTAINER_DEFAULT_CONFIG_ASSIGNMENT_COUNT,
    CONTAINER_DEFAULT_CONFIG_DENIED_SELECTORS,
    CONTAINER_EXCLUDED_PROCESS_COUNT,
    CONTAINER_INVENTORY_ENTRIES_SHA256,
    CONTAINER_INVENTORY_SHA256,
    CONTAINER_LOCK_SCHEMA_SHA256,
    CONTAINER_PROCESS_AUDIT_PROCESSES_SHA256,
    CONTAINER_PROCESS_AUDIT_SHA256,
    CONTAINER_PROCESS_COUNT,
    CONTAINER_PROCESS_UNIVERSE_COUNT,
    CONTAINER_RESERVED_DEFAULT_DENY_LABEL,
    DOCKER_REQUIRED_RUN_OPTIONS,
    JAVA_EXECUTABLE_SHA256,
    JAVA_VERSION_OUTPUT_SHA256,
    JDK_ARCHITECTURE,
    JDK_ARCHIVE_SHA256,
    JDK_RUNTIME_VERSION,
    JDK_TREE_FILE_COUNT,
    JDK_TREE_SHA256,
    JDK_VENDOR,
    JDK_VERSION,
    NETWORK_ISOLATION_EXECUTABLE_SHA256,
    NETWORK_ISOLATION_REQUIRED_ARGS,
    NETWORK_ISOLATION_TOOL,
    NETWORK_ISOLATION_VERSION,
    NETWORK_ISOLATION_VERSION_OUTPUT_SHA256,
    NEXTFLOW_SHA256,
    NEXTFLOW_VERSION,
    NF_SCHEMA_ARCHIVE_SHA256,
    NF_SCHEMA_VERSION,
    RUNTIME_IDENTITY_SHA256,
    SOURCE_FILE_COUNT,
    SOURCE_MANIFEST_SHA256,
    SOURCE_TREE_SHA256,
)
from encode_pipeline.adapters.bulk_rnaseq.upstream import (
    NFCORE_RNASEQ_COMMIT,
    NFCORE_RNASEQ_RELEASE,
    UPSTREAM_PARAMETER_SCHEMA_SHA256,
    UPSTREAM_SAMPLESHEET_SCHEMA_SHA256,
)
from encode_pipeline.platform.results import Issue, Result


DEFAULT_EXECUTION_QUALIFICATION_FILE = "default-execution-qualification-1.0.0.json"
DEFAULT_EXECUTION_QUALIFICATION_SCHEMA_VERSION = "1.0.0"
DEFAULT_EXECUTION_QUALIFICATION_POLICY_ID = (
    "protected-bulk-rnaseq-product-qualification-v1"
)
DEFAULT_EXECUTION_QUALIFICATION_POLICY_SCHEMA_VERSION = "1.0.0"
PROTECTED_GATE_EVIDENCE_SCHEMA_VERSION = "1.1.0"
PROTECTED_GATE_FIXTURE_SCHEMA_VERSION = "1.1.0"
TRANSCRIPTOME_BINDING_SCHEMA_VERSION = "1.0.0"
BUILD_IDENTITY_SCHEME = "sha256-bulk-rnaseq-runtime-v1"
WORKSPACE_SCHEMA_VERSION = "1.1.0"

_CONTRACT_PACKAGE = "encode_pipeline.contracts.nfcore_rnaseq"
_MAXIMUM_QUALIFICATION_BYTES = 64 * 1024
_SHA256 = re.compile(r"^[0-9a-f]{64}$")


class BulkRnaSeqExecutionMode(Enum):
    """Closed deployment policy; workflow inputs cannot select these modes."""

    STANDARD = "standard-v1"
    RAPID_QUANT = "rapid-quant-v1"


# nf-core/rnaseq 3.26.0's immutable ``conf/rapid_quant.config`` owns these
# values. Omitting every key from params.json is necessary because params-file
# values would otherwise override the selected profile.
RAPID_QUANT_PROFILE_OWNED_PARAMETERS = frozenset(
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

# The pinned profile relies on the base-config default ``false`` for this
# upstream branch guard. The authoring snapshot deliberately sets it to true
# for standard runs, so qualification mode must omit it alongside the values
# written directly by rapid_quant.config or Salmon would never be scheduled.
RAPID_QUANT_MODE_OWNED_PARAMETERS = RAPID_QUANT_PROFILE_OWNED_PARAMETERS | {
    "skip_pseudo_alignment"
}

_NEXTFLOW_PROFILE_BY_MODE = {
    BulkRnaSeqExecutionMode.STANDARD: None,
    BulkRnaSeqExecutionMode.RAPID_QUANT: "rapid_quant",
}


@dataclass(frozen=True)
class BulkRnaSeqExecutionQualification:
    """Exact, source-owned execution candidate used by live admission."""

    implementation: ExecutionImplementationQualification
    record_sha256: str

    def __post_init__(self) -> None:
        if (
            not isinstance(self.implementation, ExecutionImplementationQualification)
            or not isinstance(self.record_sha256, str)
            or _SHA256.fullmatch(self.record_sha256) is None
        ):
            raise ValueError("bulk RNA-seq execution qualification is invalid")


def execution_mode_identity(mode: BulkRnaSeqExecutionMode) -> str:
    """Return the stable identity coordinate for one server-owned mode."""
    if not isinstance(mode, BulkRnaSeqExecutionMode):
        raise ValueError("bulk RNA-seq execution mode is invalid")
    return mode.value


def nextflow_profile_for_mode(mode: BulkRnaSeqExecutionMode) -> str | None:
    """Return the exact scientific profile, if this mode requires one."""
    try:
        return _NEXTFLOW_PROFILE_BY_MODE[mode]
    except (KeyError, TypeError) as error:
        raise ValueError("bulk RNA-seq execution mode is invalid") from error


def runtime_params_for_mode(
    params: Mapping[str, Any],
    *,
    mode: BulkRnaSeqExecutionMode,
) -> dict[str, Any]:
    """Remove only immutable-profile-owned params for qualification mode."""
    execution_mode_identity(mode)
    excluded = (
        RAPID_QUANT_MODE_OWNED_PARAMETERS
        if mode is BulkRnaSeqExecutionMode.RAPID_QUANT
        else frozenset()
    )
    return {name: value for name, value in params.items() if name not in excluded}


def load_default_execution_qualification(
    implementation: VerifiedExecutionImplementation,
    *,
    content: bytes | None = None,
) -> Result[BulkRnaSeqExecutionQualification]:
    """Load and exactly match the fixed source candidate without private input."""
    try:
        if not isinstance(implementation, VerifiedExecutionImplementation):
            raise ValueError("verified implementation is required")
        raw = _read_qualification_resource() if content is None else _bounded(content)
        document = _strict_json(raw)
        expected = qualification_document(implementation)
        if raw != canonical_qualification_bytes(document) or _canonical_json(
            document
        ) != _canonical_json(expected):
            raise ValueError("execution qualification is stale or non-canonical")
        return Result.success(
            BulkRnaSeqExecutionQualification(
                implementation=ExecutionImplementationQualification.from_verified(
                    implementation
                ),
                record_sha256=expected["record_sha256"],
            )
        )
    except (
        OSError,
        TypeError,
        UnicodeError,
        ValueError,
        json.JSONDecodeError,
    ):
        return _qualification_failure()


def canonical_qualification_bytes(document: object) -> bytes:
    """Serialize one source candidate deterministically."""
    return _canonical_json(document) + b"\n"


def qualification_document(
    implementation: VerifiedExecutionImplementation,
) -> dict[str, Any]:
    """Build the exact path-free candidate document for an explicit update."""
    if not isinstance(implementation, VerifiedExecutionImplementation):
        raise ValueError("verified implementation is required")
    payload = _qualification_payload(implementation)
    return {
        **payload,
        "record_sha256": hashlib.sha256(_canonical_json(payload)).hexdigest(),
    }


def _qualification_payload(
    implementation: VerifiedExecutionImplementation,
) -> dict[str, Any]:
    contract = implementation.persistence_contract
    required_schema = {
        table: list(columns) for table, columns in contract.required_schema
    }
    return {
        "schema_version": DEFAULT_EXECUTION_QUALIFICATION_SCHEMA_VERSION,
        "candidate_status": "source-enabled",
        "qualification_policy": {
            "id": DEFAULT_EXECUTION_QUALIFICATION_POLICY_ID,
            "schema_version": (DEFAULT_EXECUTION_QUALIFICATION_POLICY_SCHEMA_VERSION),
            "protected_gate_evidence_schema_version": (
                PROTECTED_GATE_EVIDENCE_SCHEMA_VERSION
            ),
            "protected_gate_fixture_schema_version": (
                PROTECTED_GATE_FIXTURE_SCHEMA_VERSION
            ),
            "exact_head_required": True,
            "execution_modes": sorted(mode.value for mode in BulkRnaSeqExecutionMode),
            "nextflow_profiles": {
                mode.value: _NEXTFLOW_PROFILE_BY_MODE[mode]
                for mode in sorted(BulkRnaSeqExecutionMode, key=lambda item: item.value)
            },
            "rapid_quant_mode_owned_parameters": sorted(
                RAPID_QUANT_MODE_OWNED_PARAMETERS
            ),
        },
        "execution_implementation": {
            "manifest_file": EXECUTION_IMPLEMENTATION_MANIFEST_FILE,
            "schema_version": EXECUTION_IMPLEMENTATION_SCHEMA_VERSION,
            "scheme": EXECUTION_IMPLEMENTATION_SCHEME,
            "manifest_sha256": implementation.manifest_sha256,
            "aggregate_sha256": implementation.aggregate_sha256,
            "file_count": len(implementation.files),
            "persistence_contract": {
                "path": EXECUTION_PERSISTENCE_CONTRACT_PATH,
                "id": contract.contract_id,
                "version": contract.contract_version,
                "sha256": contract.sha256,
                "minimum_supported_revision": contract.minimum_supported_revision,
                "capabilities": list(contract.capabilities),
                "required_revisions": list(contract.required_revisions),
                "required_schema": required_schema,
                "schema_projection_sha256": contract.schema_projection_sha256,
            },
        },
        "adapter_command_contract": {
            "workflow_id": "bulk-rnaseq",
            "authoring_schema_version": SCHEMA_VERSION,
            "build_identity_scheme": BUILD_IDENTITY_SCHEME,
            "workspace_schema_version": WORKSPACE_SCHEMA_VERSION,
            "logical_entrypoint": "main.nf",
            "offline_required": True,
            "runtime_pull_allowed": False,
        },
        "runtime": {
            "identity_sha256": RUNTIME_IDENTITY_SHA256,
            "source": {
                "project": "nf-core/rnaseq",
                "release": NFCORE_RNASEQ_RELEASE,
                "commit": NFCORE_RNASEQ_COMMIT,
                "manifest_sha256": SOURCE_MANIFEST_SHA256,
                "tree_sha256": SOURCE_TREE_SHA256,
                "file_count": SOURCE_FILE_COUNT,
            },
            "nextflow": {
                "version": NEXTFLOW_VERSION,
                "sha256": NEXTFLOW_SHA256,
            },
            "jdk": {
                "vendor": JDK_VENDOR,
                "version": JDK_VERSION,
                "runtime_version": JDK_RUNTIME_VERSION,
                "architecture": JDK_ARCHITECTURE,
                "archive_sha256": JDK_ARCHIVE_SHA256,
                "tree_sha256": JDK_TREE_SHA256,
                "tree_file_count": JDK_TREE_FILE_COUNT,
                "java_sha256": JAVA_EXECUTABLE_SHA256,
                "java_version_output_sha256": JAVA_VERSION_OUTPUT_SHA256,
            },
            "plugin": {
                "id": "nf-schema",
                "version": NF_SCHEMA_VERSION,
                "archive_sha256": NF_SCHEMA_ARCHIVE_SHA256,
            },
            "containers": {
                "lock_schema_sha256": CONTAINER_LOCK_SCHEMA_SHA256,
                "inventory_sha256": CONTAINER_INVENTORY_SHA256,
                "inventory_entries_sha256": CONTAINER_INVENTORY_ENTRIES_SHA256,
                "process_count": CONTAINER_PROCESS_COUNT,
                "process_audit_sha256": CONTAINER_PROCESS_AUDIT_SHA256,
                "process_audit_processes_sha256": (
                    CONTAINER_PROCESS_AUDIT_PROCESSES_SHA256
                ),
                "process_universe_count": CONTAINER_PROCESS_UNIVERSE_COUNT,
                "excluded_process_count": CONTAINER_EXCLUDED_PROCESS_COUNT,
                "config_assignment_count": CONTAINER_CONFIG_ASSIGNMENT_COUNT,
                "config_assignments_sha256": CONTAINER_CONFIG_ASSIGNMENTS_SHA256,
                "default_config_assignment_count": (
                    CONTAINER_DEFAULT_CONFIG_ASSIGNMENT_COUNT
                ),
                "default_config_assignments_sha256": (
                    CONTAINER_DEFAULT_CONFIG_ASSIGNMENTS_SHA256
                ),
                "default_config_denied_selectors": list(
                    CONTAINER_DEFAULT_CONFIG_DENIED_SELECTORS
                ),
                "reserved_default_deny_label": CONTAINER_RESERVED_DEFAULT_DENY_LABEL,
                "required_run_options": list(DOCKER_REQUIRED_RUN_OPTIONS),
                "runtime_pull_allowed": False,
            },
            "network_isolation": {
                "tool": NETWORK_ISOLATION_TOOL,
                "version": NETWORK_ISOLATION_VERSION,
                "executable_sha256": NETWORK_ISOLATION_EXECUTABLE_SHA256,
                "version_output_sha256": NETWORK_ISOLATION_VERSION_OUTPUT_SHA256,
                "required_args": list(NETWORK_ISOLATION_REQUIRED_ARGS),
            },
        },
        "reference_identity_protocol": {
            "transcriptome_binding_schema_version": (
                TRANSCRIPTOME_BINDING_SCHEMA_VERSION
            ),
            "reference_closure_scheme": REFERENCE_CLOSURE_SCHEME,
            "reference_index_manifest": REFERENCE_INDEX_MANIFEST,
            "reference_index_schema_version": REFERENCE_INDEX_SCHEMA_VERSION,
            "ribo_database_closure_scheme": RIBO_DATABASE_CLOSURE_SCHEME,
            "sortmerna_index_closure_scheme": SORTMERNA_INDEX_CLOSURE_SCHEME,
            "sortmerna_index_manifest_schema_version": (
                SORTMERNA_INDEX_MANIFEST_SCHEMA_VERSION
            ),
            "sortmerna_version": SORTMERNA_VERSION,
            "sortmerna_no_prebuilt_index_strategy": (
                SORTMERNA_NO_PREBUILT_INDEX_STRATEGY
            ),
            "resource_closure_policy": asdict(DEFAULT_RESOURCE_CLOSURE_POLICY),
            "upstream_parameter_schema_sha256": UPSTREAM_PARAMETER_SCHEMA_SHA256,
            "upstream_samplesheet_schema_sha256": (UPSTREAM_SAMPLESHEET_SCHEMA_SHA256),
        },
        "result_contract": {
            "file": RESULTS_CONTRACT_FILE,
            "size_bytes": RESULTS_CONTRACT_SIZE,
            "sha256": RESULTS_CONTRACT_SHA256,
        },
    }


def _read_qualification_resource() -> bytes:
    resource = resources.files(_CONTRACT_PACKAGE).joinpath(
        DEFAULT_EXECUTION_QUALIFICATION_FILE
    )
    if isinstance(resource, Path):
        result = safe_regular_file_bytes(
            resource,
            maximum_bytes=_MAXIMUM_QUALIFICATION_BYTES,
        )
        if result.is_failure:
            raise ValueError("qualification resource is unavailable")
        return result.value.content
    is_file = getattr(resource, "is_file", None)
    read_bytes = getattr(resource, "read_bytes", None)
    if not callable(is_file) or not is_file() or not callable(read_bytes):
        raise ValueError("qualification resource is unavailable")
    return _bounded(read_bytes())


def _bounded(value: object) -> bytes:
    if (
        not isinstance(value, bytes)
        or not value
        or len(value) > _MAXIMUM_QUALIFICATION_BYTES
    ):
        raise ValueError("qualification resource bounds are invalid")
    return value


def _strict_json(content: bytes) -> object:
    def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, item in pairs:
            if key in value:
                raise ValueError("duplicate qualification field")
            value[key] = item
        return value

    def reject_constant(_value: str) -> None:
        raise ValueError("non-finite qualification number")

    return json.loads(
        _bounded(content).decode("utf-8"),
        object_pairs_hook=unique_object,
        parse_constant=reject_constant,
    )


def _canonical_json(value: object) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


def _qualification_failure() -> Result[BulkRnaSeqExecutionQualification]:
    return Result.failure(
        (
            Issue(
                code="BULK_RNASEQ_EXECUTION_QUALIFICATION_INVALID",
                message="The bulk RNA-seq execution qualification is invalid.",
                severity="error",
                path="runtime.qualification",
                source="adapter",
                context={},
            ),
        )
    )
