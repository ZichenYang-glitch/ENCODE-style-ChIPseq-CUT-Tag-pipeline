"""Private composition and evidence helpers for the real bulk RNA-seq gate.

Nothing in this module is imported by the product registry or API.  The gate
uses the production adapter, services, SQLite repository, RQ job, and worker
runtime through the deployment-owned composition seams already exposed by the
platform.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import asdict, dataclass, field
from hashlib import sha256
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import subprocess
from typing import Any, Callable

from encode_pipeline.adapters.bulk_rnaseq import (
    BulkRnaSeqExecutionBinding,
    BulkRnaSeqResultsWorkflowAdapter,
    BulkRnaSeqTranscriptomeBinding,
    RuntimeAssetBinding,
)
from encode_pipeline.adapters.bulk_rnaseq.execution import WORKSPACE_SCHEMA_VERSION
from encode_pipeline.adapters.bulk_rnaseq.reference_closure import (
    verify_reference_index,
)
from encode_pipeline.adapters.bulk_rnaseq.resource_closure import (
    verify_ribo_database_manifest,
    verify_sortmerna_index,
)
from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import NEXTFLOW_VERSION
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.managed_containers import (
    MANAGED_CONTAINER_SCOPE_LABEL,
    managed_container_scope,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.result_generations import (
    artifact_manifest_digest,
    build_artifact_generation,
    build_qc_generation,
    qc_metric_manifest_digest,
    validate_artifact_generation,
    validate_qc_generation,
    validate_result_attempt_id,
)
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.artifact_downloads import ArtifactDownloadService
from encode_pipeline.services.managed_containers import ManagedContainerCleaner
from encode_pipeline.services.process_runner import ProcessRunner
from encode_pipeline.services.run_repositories import (
    canonical_decimal_text,
    decimal_from_canonical_text,
)
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider


REQUIRE_REAL_EXECUTION_ENV = "HELIXWEAVE_REQUIRE_BULK_RNASEQ_REAL_EXECUTION"
RUNTIME_ROOT_ENV = "HELIXWEAVE_BULK_RNASEQ_RUNTIME_ROOT"
FIXTURE_MANIFEST_ENV = "HELIXWEAVE_BULK_RNASEQ_FIXTURE_MANIFEST"
TEST_REDIS_URL_ENV = "ENCODE_PIPELINE_TEST_REDIS_URL"
MANAGED_DOCKER_EXECUTABLE_ENV = "ENCODE_PIPELINE_MANAGED_DOCKER_EXECUTABLE"
MANAGED_DOCKER_SOCKET_ENV = "ENCODE_PIPELINE_MANAGED_DOCKER_SOCKET"

FIXTURE_MANIFEST_SCHEMA_VERSION = "1.1.0"
ACCEPTANCE_EVIDENCE_SCHEMA_VERSION = "1.1.0"
SOURCE_MANIFEST_FILENAME = "helixweave-bulk-rnaseq-tiny-source.json"
INDEX_PROVENANCE_FILENAME = "helixweave-bulk-rnaseq-index-provenance.json"
SOURCE_MANIFEST_SCHEMA_VERSION = "1.0.0"
INDEX_PROVENANCE_SCHEMA_VERSION = "1.0.0"
FIXTURE_ID = "helixweave-bulk-rnaseq-tiny-v1"
SOURCE_IDENTITY_SCHEME = "sha256-framed-bulk-rnaseq-tiny-source-v1"
PROVENANCE_IDENTITY_SCHEME = "sha256-framed-bulk-rnaseq-index-provenance-v1"
_MAX_FIXTURE_MANIFEST_BYTES = 1024 * 1024
_MAX_FIXTURE_SOURCE_FILE_BYTES = 64 * 1024 * 1024
_MAX_FIXTURE_SOURCE_TOTAL_BYTES = 256 * 1024 * 1024
_MAX_EVIDENCE_ARTIFACT_BYTES = 512 * 1024 * 1024
_MAX_EVIDENCE_TOTAL_BYTES = 2 * 1024 * 1024 * 1024
_MAX_IDENTITY_DOCUMENT_BYTES = 64 * 1024
_HEX64 = re.compile(r"^[0-9a-f]{64}$")
_GIT_COMMIT_ID = re.compile(r"^(?:[0-9a-f]{40}|[0-9a-f]{64})$")
_PUBLIC_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,255}$")
_OUTPUT_TYPE = re.compile(r"^[a-z][a-z0-9_.-]{0,255}$")
_METRIC_KEY = re.compile(r"^[a-z][a-z0-9_.-]{0,255}$")
_CONTAINER_ID = re.compile(r"^[0-9a-f]{64}$")
_IMMUTABLE_IMAGE = re.compile(r"^[^@\s]+@sha256:[0-9a-f]{64}$")
_WINDOWS_ABSOLUTE = re.compile(r"^[A-Za-z]:[\\/]")
_EMBEDDED_ABSOLUTE = re.compile(r"(?:^|[=,\s])/(?!/)")
_EXPECTED_RAW_FASTQC_VALUES = (
    ("PE1", "fastqc.raw.read1.total_sequences", "768"),
    ("PE1", "fastqc.raw.read2.total_sequences", "768"),
    ("SE1", "fastqc.raw.single.total_sequences", "384"),
)
_EXPECTED_INDEX_PRODUCER_IMAGES = {
    "star": (
        "community.wave.seqera.io/library/htslib_samtools_star_gawk:"
        "ae438e9a604351a4@sha256:"
        "4a468118dbd7491a69bf9813c68233afa8558d1f3380fd8cab03e0e3d3135190"
    ),
    "salmon": (
        "quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2@sha256:"
        "f83ebb158845ee8138d793347f83b92c75e83c58dd8f4600c6fea2a2453ef08e"
    ),
    "sortmerna": (
        "community.wave.seqera.io/library/sortmerna:4.3.7--b730cad73fc42b8e"
        "@sha256:"
        "3c873f2a4c007c17b3b30aedab6b1d0d0670a9c629033686508c2d80d780a4af"
    ),
}


class AcceptanceEvidenceStaleError(ValueError):
    """The supplied evidence no longer matches its canonical closure."""


@dataclass(frozen=True)
class GateSettings:
    """Private operator coordinates required by the explicit real gate."""

    runtime_root: Path
    fixture_manifest: Path
    redis_url: str = field(repr=False)
    docker_executable: Path = field(repr=False)
    docker_socket: Path = field(repr=False)

    def __post_init__(self) -> None:
        for name in (
            "runtime_root",
            "fixture_manifest",
            "docker_executable",
            "docker_socket",
        ):
            value = getattr(self, name)
            if (
                not isinstance(value, Path)
                or not value.is_absolute()
                or str(value) != str(Path(value))
                or any(part == ".." for part in value.parts)
                or any(character in str(value) for character in ("\x00", "\n", "\r"))
            ):
                raise ValueError(f"{name} must be a canonical absolute Path")
        if self.docker_executable.name != "docker":
            raise ValueError("docker_executable must name the Docker CLI")
        if not isinstance(self.redis_url, str) or not self.redis_url.strip():
            raise ValueError("redis_url must be non-empty")


@dataclass(frozen=True)
class ResultsComposition:
    """One explicit result-capable registry and its exact runtime binding."""

    binding: BulkRnaSeqExecutionBinding
    registry: WorkflowRegistry
    build_identity_provider: WorkflowBuildIdentityProvider


@dataclass(frozen=True)
class AcceptanceFixture:
    """Controlled full STAR+Salmon fixture and expected public result surface."""

    workflow_inputs: WorkflowInputs
    transcriptome: BulkRnaSeqTranscriptomeBinding
    acceptance_manifest_sha256: str
    source_manifest_sha256: str
    source_identity_sha256: str
    index_provenance_manifest_sha256: str
    index_provenance_identity_sha256: str
    required_artifact_output_types: tuple[str, ...]
    required_qc_metric_keys: tuple[str, ...]
    required_sample_ids: tuple[str, ...]
    required_artifact_sample_output_types: tuple[tuple[str, str], ...]
    required_qc_sample_metric_keys: tuple[tuple[str, str], ...]
    required_qc_sample_metric_values: tuple[tuple[str, str, str], ...]


@dataclass(frozen=True)
class AcceptanceEvidenceValues:
    """Path-free coordinates that make prior acceptance evidence stale."""

    tested_head: str
    workflow_build_digest: str
    fixture_acceptance_manifest_sha256: str
    fixture_source_manifest_sha256: str
    fixture_source_identity_sha256: str
    fixture_index_provenance_manifest_sha256: str
    fixture_index_provenance_identity_sha256: str
    validated_snapshot_id: str
    validated_payload_digest: str
    cache_identity_sha256: str
    input_closure_sha256: str
    ribo_database_closure_sha256: str
    workspace_identity_sha256: str
    workspace_contract_sha256: str
    execution_implementation_manifest_sha256: str
    execution_implementation_aggregate_sha256: str
    container_process_audit_sha256: str
    run_id: str
    job_id: str
    lifecycle_status: str
    artifact_attempt_id: str
    artifact_revision: int
    artifact_generation: str
    artifact_manifest_digest: str
    artifact_content_sha256: tuple[tuple[str, str], ...]
    qc_attempt_id: str
    qc_revision: int
    qc_generation: str
    qc_manifest_digest: str
    qc_artifact_generation: str
    artifact_output_types: tuple[str, ...]
    qc_metric_keys: tuple[str, ...]
    qc_sample_ids: tuple[str, ...]
    artifact_sample_output_types: tuple[tuple[str, str], ...]
    qc_sample_metric_keys: tuple[tuple[str, str], ...]
    qc_sample_metric_values: tuple[tuple[str, str, str], ...]
    cleanup_confirmed: bool

    def __post_init__(self) -> None:
        for name in (
            "workflow_build_digest",
            "fixture_acceptance_manifest_sha256",
            "fixture_source_manifest_sha256",
            "fixture_source_identity_sha256",
            "fixture_index_provenance_manifest_sha256",
            "fixture_index_provenance_identity_sha256",
            "validated_payload_digest",
            "cache_identity_sha256",
            "input_closure_sha256",
            "ribo_database_closure_sha256",
            "workspace_identity_sha256",
            "workspace_contract_sha256",
            "execution_implementation_manifest_sha256",
            "execution_implementation_aggregate_sha256",
            "container_process_audit_sha256",
            "artifact_manifest_digest",
            "qc_manifest_digest",
        ):
            _require_digest(getattr(self, name), name)
        if _GIT_COMMIT_ID.fullmatch(self.tested_head) is None:
            raise ValueError("tested_head must be a Git commit identity")
        for name in ("validated_snapshot_id", "run_id", "job_id"):
            if (
                not isinstance(getattr(self, name), str)
                or _PUBLIC_ID.fullmatch(getattr(self, name)) is None
            ):
                raise ValueError(f"{name} is invalid")
        if self.lifecycle_status != "succeeded":
            raise ValueError("acceptance evidence requires a succeeded lifecycle")
        validate_result_attempt_id(self.artifact_attempt_id)
        validate_result_attempt_id(self.qc_attempt_id)
        if self.artifact_attempt_id == self.qc_attempt_id:
            raise ValueError("artifact and QC attempt identities must be distinct")
        validate_artifact_generation(self.artifact_generation)
        validate_artifact_generation(self.qc_artifact_generation)
        validate_qc_generation(self.qc_generation)
        if self.qc_artifact_generation != self.artifact_generation:
            raise ValueError("QC evidence must bind the accepted artifact generation")
        for name in ("artifact_revision", "qc_revision"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise ValueError(f"{name} must be positive")
        content = tuple(self.artifact_content_sha256)
        if (
            not content
            or content != tuple(sorted(content))
            or len({artifact_id for artifact_id, _digest in content}) != len(content)
        ):
            raise ValueError("artifact content evidence must be non-empty and sorted")
        for artifact_id, digest in content:
            if _PUBLIC_ID.fullmatch(artifact_id) is None:
                raise ValueError("artifact content identity is invalid")
            _require_digest(digest, "artifact content")
        object.__setattr__(self, "artifact_content_sha256", content)
        artifact_types = _closed_sorted_tokens(
            self.artifact_output_types,
            pattern=_OUTPUT_TYPE,
            name="artifact output types",
        )
        metric_keys = _closed_sorted_tokens(
            self.qc_metric_keys,
            pattern=_METRIC_KEY,
            name="QC metric keys",
        )
        sample_ids = _closed_sorted_tokens(
            self.qc_sample_ids,
            pattern=_PUBLIC_ID,
            name="QC sample IDs",
        )
        artifact_coordinates = _closed_sorted_pairs(
            self.artifact_sample_output_types,
            first_pattern=_PUBLIC_ID,
            second_pattern=_OUTPUT_TYPE,
            name="artifact sample/output coordinates",
        )
        metric_coordinates = _closed_sorted_pairs(
            self.qc_sample_metric_keys,
            first_pattern=_PUBLIC_ID,
            second_pattern=_METRIC_KEY,
            name="QC sample/metric coordinates",
        )
        metric_values = _closed_sorted_metric_values(
            self.qc_sample_metric_values,
            name="QC sample/metric values",
        )
        object.__setattr__(self, "artifact_output_types", artifact_types)
        object.__setattr__(self, "qc_metric_keys", metric_keys)
        object.__setattr__(self, "qc_sample_ids", sample_ids)
        object.__setattr__(
            self,
            "artifact_sample_output_types",
            artifact_coordinates,
        )
        object.__setattr__(self, "qc_sample_metric_keys", metric_coordinates)
        object.__setattr__(self, "qc_sample_metric_values", metric_values)
        if self.cleanup_confirmed is not True:
            raise ValueError("acceptance evidence requires confirmed cleanup")


@dataclass(frozen=True)
class AcceptanceEvidence:
    """Canonical, self-digesting, path-free evidence from one accepted run."""

    values: AcceptanceEvidenceValues
    evidence_sha256: str

    def __post_init__(self) -> None:
        if not isinstance(self.values, AcceptanceEvidenceValues):
            raise ValueError("values must be AcceptanceEvidenceValues")
        _require_digest(self.evidence_sha256, "acceptance evidence")
        expected = _evidence_digest(self.values)
        if self.evidence_sha256 != expected:
            raise AcceptanceEvidenceStaleError(
                "acceptance evidence does not match its canonical closure"
            )

    @classmethod
    def create(cls, values: AcceptanceEvidenceValues) -> AcceptanceEvidence:
        """Create evidence bound to every supplied acceptance coordinate."""
        return cls(values=values, evidence_sha256=_evidence_digest(values))

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> AcceptanceEvidence:
        """Validate one serialized evidence document and its closure digest."""
        if not isinstance(payload, Mapping):
            raise ValueError("acceptance evidence must be an object")
        expected_keys = {
            "schema_version",
            "evidence_sha256",
            *AcceptanceEvidenceValues.__dataclass_fields__,
        }
        if set(payload) != expected_keys:
            raise ValueError("acceptance evidence fields are invalid")
        if payload["schema_version"] != ACCEPTANCE_EVIDENCE_SCHEMA_VERSION:
            raise ValueError("acceptance evidence schema is unsupported")
        normalized = dict(payload)
        normalized.pop("schema_version")
        evidence_digest = normalized.pop("evidence_sha256")
        for name in (
            "artifact_content_sha256",
            "artifact_output_types",
            "qc_metric_keys",
            "qc_sample_ids",
            "artifact_sample_output_types",
            "qc_sample_metric_keys",
            "qc_sample_metric_values",
        ):
            value = normalized[name]
            if not isinstance(value, list):
                raise ValueError(f"{name} must be an array")
            if name in {
                "artifact_content_sha256",
                "artifact_sample_output_types",
                "qc_sample_metric_keys",
                "qc_sample_metric_values",
            }:
                expected_length = 3 if name == "qc_sample_metric_values" else 2
                if any(
                    not isinstance(item, list) or len(item) != expected_length
                    for item in value
                ):
                    raise ValueError(f"{name} entries are invalid")
                normalized[name] = tuple(tuple(item) for item in value)
            else:
                normalized[name] = tuple(value)
        try:
            values = AcceptanceEvidenceValues(**normalized)
            return cls(values=values, evidence_sha256=evidence_digest)
        except AcceptanceEvidenceStaleError:
            raise
        except (TypeError, ValueError):
            raise ValueError("acceptance evidence values are invalid") from None

    def to_dict(self) -> dict[str, object]:
        """Return one canonical JSON-compatible document without private paths."""
        payload = _values_payload(self.values)
        return {
            "schema_version": ACCEPTANCE_EVIDENCE_SCHEMA_VERSION,
            **payload,
            "evidence_sha256": self.evidence_sha256,
        }

    def assert_matches(self, current: AcceptanceEvidence) -> None:
        """Reject reuse when any current acceptance coordinate has changed."""
        if not isinstance(current, AcceptanceEvidence):
            raise ValueError("current must be AcceptanceEvidence")
        if self.evidence_sha256 != current.evidence_sha256:
            raise AcceptanceEvidenceStaleError("acceptance evidence is stale")


def require_gate_settings(
    environ: Mapping[str, str] | None = None,
    *,
    _socket_probe: Callable[[Path], bool] | None = None,
) -> GateSettings:
    """Load an explicitly enabled gate or fail; this function never skips."""
    source = os.environ if environ is None else environ
    required = (
        REQUIRE_REAL_EXECUTION_ENV,
        RUNTIME_ROOT_ENV,
        FIXTURE_MANIFEST_ENV,
        TEST_REDIS_URL_ENV,
        MANAGED_DOCKER_EXECUTABLE_ENV,
        MANAGED_DOCKER_SOCKET_ENV,
    )
    for name in required:
        value = source.get(name)
        if value is None or not value.strip():
            raise AssertionError(f"{name} is required by the explicit real gate")
    if source[REQUIRE_REAL_EXECUTION_ENV] != "1":
        raise AssertionError(f"{REQUIRE_REAL_EXECUTION_ENV} must be exactly 1")
    settings = GateSettings(
        runtime_root=Path(source[RUNTIME_ROOT_ENV]),
        fixture_manifest=Path(source[FIXTURE_MANIFEST_ENV]),
        redis_url=source[TEST_REDIS_URL_ENV],
        docker_executable=Path(source[MANAGED_DOCKER_EXECUTABLE_ENV]),
        docker_socket=Path(source[MANAGED_DOCKER_SOCKET_ENV]),
    )
    if not settings.runtime_root.is_dir() or settings.runtime_root.is_symlink():
        raise AssertionError(f"{RUNTIME_ROOT_ENV} must name the staged runtime root")
    if (
        not settings.fixture_manifest.is_file()
        or settings.fixture_manifest.is_symlink()
    ):
        raise AssertionError(
            f"{FIXTURE_MANIFEST_ENV} must name the generated fixture manifest"
        )
    if (
        not settings.docker_executable.is_file()
        or settings.docker_executable.is_symlink()
        or not os.access(settings.docker_executable, os.X_OK)
    ):
        raise AssertionError(
            f"{MANAGED_DOCKER_EXECUTABLE_ENV} must name the fixed Docker CLI"
        )
    socket_probe = _is_unix_socket if _socket_probe is None else _socket_probe
    if not socket_probe(settings.docker_socket):
        raise AssertionError(
            f"{MANAGED_DOCKER_SOCKET_ENV} must name the local Docker socket"
        )
    return settings


def _is_unix_socket(path: Path) -> bool:
    try:
        return stat.S_ISSOCK(path.stat().st_mode)
    except OSError:
        return False


def build_results_composition(
    settings: GateSettings,
    *,
    project_root: Path | None = None,
) -> ResultsComposition:
    """Build the explicit results registry without changing platform defaults."""
    if not isinstance(settings, GateSettings):
        raise ValueError("settings must be GateSettings")
    fixture = load_acceptance_fixture(settings.fixture_manifest)
    binding = BulkRnaSeqExecutionBinding(
        assets=RuntimeAssetBinding(
            root=settings.runtime_root,
            docker_executable=settings.docker_executable,
            docker_socket=settings.docker_socket,
        ),
        transcriptome=fixture.transcriptome,
    )
    registry = WorkflowRegistry((BulkRnaSeqResultsWorkflowAdapter(execution=binding),))
    provider = WorkflowBuildIdentityProvider(registry, project_root=project_root)
    return ResultsComposition(
        binding=binding,
        registry=registry,
        build_identity_provider=provider,
    )


def build_acceptance_process_runner(
    *,
    settings: GateSettings,
    binding: BulkRnaSeqExecutionBinding,
    timeout_seconds: float,
    passthrough_exceptions: tuple[type[BaseException], ...] = (),
) -> ProcessRunner:
    """Bind the private gate runner to its admitted isolation and Docker assets."""
    if not isinstance(settings, GateSettings):
        raise ValueError("settings must be GateSettings")
    try:
        assets = binding.assets
        network_isolation_executable = assets.network_isolation_executable
    except AttributeError:
        raise ValueError("execution binding lacks network isolation assets") from None
    if (
        assets.docker_executable != settings.docker_executable
        or assets.docker_socket != settings.docker_socket
    ):
        raise ValueError("execution binding Docker endpoint differs from gate settings")
    if (
        not isinstance(network_isolation_executable, Path)
        or not network_isolation_executable.is_absolute()
        or network_isolation_executable.name != "unshare"
    ):
        raise ValueError("network isolation executable must be exact unshare")
    cleaner = ManagedContainerCleaner(
        executable=settings.docker_executable,
        unix_socket=settings.docker_socket,
    )
    return ProcessRunner(
        allowed_executables=(str(network_isolation_executable),),
        timeout_seconds=timeout_seconds,
        passthrough_exceptions=passthrough_exceptions,
        managed_container_cleaner=cleaner,
    )


def managed_container_ids(
    cleaner: ManagedContainerCleaner,
    scope: str,
    *,
    all_containers: bool,
) -> tuple[str, ...]:
    """Read one managed scope without mutating its containers."""
    if not isinstance(scope, str) or _HEX64.fullmatch(scope) is None:
        raise ValueError("managed container scope is invalid")
    if cleaner.verify_endpoint().is_failure:
        raise AssertionError("managed Docker endpoint changed during acceptance")
    argv = [
        str(cleaner.executable),
        "--host",
        cleaner.local_docker_host,
        "ps",
    ]
    if all_containers:
        argv.append("--all")
    argv.extend(
        (
            "--quiet",
            "--no-trunc",
            "--filter",
            f"label={MANAGED_CONTAINER_SCOPE_LABEL}={scope}",
        )
    )
    try:
        completed = subprocess.run(
            tuple(argv),
            shell=False,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=5,
            check=False,
            env={"PATH": "/usr/bin:/bin"},
        )
    except (OSError, subprocess.SubprocessError):
        raise AssertionError("managed container activity cannot be audited") from None
    if completed.returncode != 0:
        raise AssertionError("managed container activity cannot be audited")
    identities = tuple(
        line.strip() for line in completed.stdout.splitlines() if line.strip()
    )
    if (
        len(identities) > 128
        or len(identities) != len(set(identities))
        or any(_CONTAINER_ID.fullmatch(identity) is None for identity in identities)
    ):
        raise AssertionError("managed container activity evidence is invalid")
    return identities


def assert_no_managed_containers(
    cleaner: ManagedContainerCleaner,
    scope: str,
) -> None:
    """Fail before hygiene cleanup if execution left a scoped container."""
    if managed_container_ids(cleaner, scope, all_containers=True):
        raise AssertionError("execution left a managed container residual")


def collect_success_evidence(
    *,
    run_service: RunService,
    run_id: str,
    expected_job_id: str,
    validated_snapshot_id: str,
    fixture: AcceptanceFixture,
    workspace_root: Path,
    repository_root: Path,
    cleaner: ManagedContainerCleaner,
) -> AcceptanceEvidence:
    """Reopen and bind one complete successful platform result projection."""
    if not isinstance(run_service, RunService):
        raise ValueError("run_service must be RunService")
    if not isinstance(fixture, AcceptanceFixture):
        raise ValueError("fixture must be AcceptanceFixture")
    if (
        not isinstance(expected_job_id, str)
        or _PUBLIC_ID.fullmatch(expected_job_id) is None
    ):
        raise ValueError("expected job identity is invalid")
    if (
        not isinstance(validated_snapshot_id, str)
        or _PUBLIC_ID.fullmatch(validated_snapshot_id) is None
    ):
        raise ValueError("validated snapshot identity is invalid")
    if not isinstance(workspace_root, Path) or not workspace_root.is_absolute():
        raise ValueError("workspace_root must be absolute")
    if not isinstance(repository_root, Path) or not repository_root.is_absolute():
        raise ValueError("repository_root must be absolute")
    if not isinstance(cleaner, ManagedContainerCleaner):
        raise ValueError("cleaner must be ManagedContainerCleaner")

    record = run_service.get_run(run_id)
    if record.status is not RunStatus.SUCCEEDED:
        raise AssertionError("accepted run lifecycle is not succeeded")
    assignment = run_service.get_execution_assignment(run_id)
    if (
        assignment is None
        or assignment.job_id != expected_job_id
        or assignment.dispatched_at is None
        or assignment.claimed_at is None
    ):
        raise AssertionError("accepted run lacks durable RQ ownership evidence")
    build_identity = run_service.get_workflow_build_identity(run_id)
    if build_identity is None:
        raise AssertionError("accepted run lacks workflow build identity")
    snapshot = run_service.get_validated_input_snapshot(validated_snapshot_id)
    if (
        snapshot.consumed_run_id != run_id
        or snapshot.consumed_at is None
        or snapshot.workflow_build_identity.digest != build_identity.digest
        or snapshot.to_workflow_inputs().to_dict() != record.inputs
    ):
        raise AssertionError("accepted run does not match its consumed snapshot")
    events = run_service.list_events(run_id, limit=1000)
    status_history = tuple(
        event.status
        for event in events
        if event.event_type == "status_changed" and event.status is not None
    )
    if status_history[-3:] != (
        RunStatus.QUEUED,
        RunStatus.RUNNING,
        RunStatus.SUCCEEDED,
    ):
        raise AssertionError("accepted run lacks the canonical lifecycle history")
    event_types = tuple(event.event_type for event in events)
    for event_type in (
        "worker_dependencies_rebuilt",
        "artifacts_indexed",
        "qc_metrics_indexed",
    ):
        if event_types.count(event_type) != 1:
            raise AssertionError(
                f"accepted run must persist exactly one {event_type} event"
            )
    if any(
        event_type in event_types
        for event_type in (
            "artifact_extraction_failed",
            "execution_cleanup_failed",
            "qc_metrics_indexing_failed",
        )
    ):
        raise AssertionError("accepted run contains a result or cleanup failure")
    state = run_service.get_result_state(run_id)
    if (
        state.artifact_attempt_status != "succeeded"
        or state.artifact_outcome != "succeeded"
        or state.qc_attempt_status != "succeeded"
        or state.qc_outcome != "succeeded"
        or state.artifact_attempt_id is None
        or state.qc_attempt_id is None
        or state.artifact_generation is None
        or state.artifact_manifest_digest is None
        or state.qc_generation is None
        or state.qc_manifest_digest is None
        or state.qc_artifact_generation is None
        or state.qc_attempt_artifact_generation != state.artifact_generation
        or state.qc_artifact_generation != state.artifact_generation
    ):
        raise AssertionError("accepted run lacks complete artifact/QC attempt evidence")

    artifacts = run_service.list_artifacts(run_id)
    metrics = run_service.list_qc_metrics(run_id)
    if not artifacts or not metrics:
        raise AssertionError("accepted run must persist non-empty artifacts and QC")
    if artifact_manifest_digest(artifacts) != state.artifact_manifest_digest:
        raise AssertionError("artifact manifest digest does not match persisted state")
    if qc_metric_manifest_digest(metrics) != state.qc_manifest_digest:
        raise AssertionError("QC manifest digest does not match persisted state")
    if (
        build_artifact_generation(
            run_id=run_id,
            revision=state.artifact_revision,
            artifacts=artifacts,
        )
        != state.artifact_generation
    ):
        raise AssertionError("artifact generation does not match its revision closure")
    if (
        build_qc_generation(
            run_id=run_id,
            revision=state.qc_revision,
            artifact_generation=state.artifact_generation,
            metrics=metrics,
        )
        != state.qc_generation
    ):
        raise AssertionError("QC generation does not match its revision closure")
    artifact_ids = {artifact.artifact_id for artifact in artifacts}
    if any(metric.source_artifact_id not in artifact_ids for metric in metrics):
        raise AssertionError("QC metric source is outside the artifact generation")

    workspace = workspace_root / run_id
    cache_identity = _read_identity_document(
        workspace / "engine" / "cache-identity.json"
    )
    execution_identity = _read_identity_document(
        workspace / "config" / "execution-identity.json"
    )
    workspace_identity = _verified_workspace_identity(
        cache_identity=cache_identity,
        execution_identity=execution_identity,
        workspace=workspace,
        workflow_build_digest=build_identity.digest,
    )
    assert_no_managed_containers(cleaner, workspace_identity)
    content_hashes = _hash_indexed_artifacts(
        run_service=run_service,
        workspace_root=workspace_root,
        run_id=run_id,
        artifact_generation=state.artifact_generation,
        artifacts=artifacts,
    )
    raw_output_types = tuple(
        artifact.metadata.get("output_type") for artifact in artifacts
    )
    if any(not isinstance(value, str) for value in raw_output_types):
        raise AssertionError("persisted artifact output type is invalid")
    output_types = tuple(sorted(set(raw_output_types)))
    metric_keys = tuple(sorted({metric.metric_key for metric in metrics}))
    qc_sample_ids = tuple(
        sorted({metric.sample_id for metric in metrics if metric.sample_id is not None})
    )
    artifact_sample_output_types: list[tuple[str, str]] = []
    for artifact in artifacts:
        output_type = artifact.metadata.get("output_type")
        sample_id = artifact.metadata.get("sample_id")
        if sample_id is None:
            continue
        if not isinstance(sample_id, str) or not isinstance(output_type, str):
            raise AssertionError("persisted artifact sample coordinate is invalid")
        artifact_sample_output_types.append((sample_id, output_type))
    qc_sample_metric_keys = tuple(
        sorted(
            (metric.sample_id, metric.metric_key)
            for metric in metrics
            if metric.sample_id is not None
        )
    )
    qc_sample_metric_values = tuple(
        sorted(
            (
                metric.sample_id,
                metric.metric_key,
                canonical_decimal_text(metric.value),
            )
            for metric in metrics
            if metric.sample_id is not None
        )
    )

    values = AcceptanceEvidenceValues(
        tested_head=read_tested_head(repository_root),
        workflow_build_digest=build_identity.digest,
        fixture_acceptance_manifest_sha256=fixture.acceptance_manifest_sha256,
        fixture_source_manifest_sha256=fixture.source_manifest_sha256,
        fixture_source_identity_sha256=fixture.source_identity_sha256,
        fixture_index_provenance_manifest_sha256=(
            fixture.index_provenance_manifest_sha256
        ),
        fixture_index_provenance_identity_sha256=(
            fixture.index_provenance_identity_sha256
        ),
        validated_snapshot_id=snapshot.snapshot_id,
        validated_payload_digest=snapshot.payload_digest,
        cache_identity_sha256=_identity_digest(cache_identity, "identity_sha256"),
        input_closure_sha256=_identity_digest(
            cache_identity,
            "input_closure_sha256",
        ),
        ribo_database_closure_sha256=_identity_digest(
            cache_identity,
            "ribo_database_closure_sha256",
        ),
        workspace_identity_sha256=workspace_identity,
        workspace_contract_sha256=_identity_digest(
            execution_identity,
            "workspace_contract_sha256",
        ),
        execution_implementation_manifest_sha256=_identity_digest(
            execution_identity,
            "execution_implementation_manifest_sha256",
        ),
        execution_implementation_aggregate_sha256=_identity_digest(
            execution_identity,
            "execution_implementation_aggregate_sha256",
        ),
        container_process_audit_sha256=_identity_digest(
            execution_identity,
            "container_process_audit_sha256",
        ),
        run_id=run_id,
        job_id=assignment.job_id,
        lifecycle_status=record.status.value,
        artifact_attempt_id=state.artifact_attempt_id,
        artifact_revision=state.artifact_revision,
        artifact_generation=state.artifact_generation,
        artifact_manifest_digest=state.artifact_manifest_digest,
        artifact_content_sha256=content_hashes,
        qc_attempt_id=state.qc_attempt_id,
        qc_revision=state.qc_revision,
        qc_generation=state.qc_generation,
        qc_manifest_digest=state.qc_manifest_digest,
        qc_artifact_generation=state.qc_artifact_generation,
        artifact_output_types=output_types,
        qc_metric_keys=metric_keys,
        qc_sample_ids=qc_sample_ids,
        artifact_sample_output_types=tuple(sorted(artifact_sample_output_types)),
        qc_sample_metric_keys=qc_sample_metric_keys,
        qc_sample_metric_values=qc_sample_metric_values,
        cleanup_confirmed=True,
    )
    return AcceptanceEvidence.create(values)


def read_tested_head(repository_root: Path) -> str:
    """Return the exact clean checked-out commit without disclosing its path."""
    try:
        head = subprocess.run(
            ("git", "rev-parse", "--verify", "HEAD^{commit}"),
            cwd=repository_root,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=10,
            check=False,
        )
        status = subprocess.run(
            ("git", "status", "--porcelain=v1", "--untracked-files=all"),
            cwd=repository_root,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            text=True,
            timeout=10,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        raise AssertionError("tested Git HEAD could not be resolved") from None
    value = head.stdout.strip()
    if head.returncode != 0 or _GIT_COMMIT_ID.fullmatch(value) is None:
        raise AssertionError("tested Git HEAD could not be resolved")
    if status.returncode != 0 or status.stdout:
        raise AssertionError("real acceptance requires a clean exact Git HEAD")
    return value


def write_acceptance_evidence(evidence: AcceptanceEvidence, path: Path) -> None:
    """Write one canonical path-free evidence document for CI archival."""
    if not isinstance(evidence, AcceptanceEvidence):
        raise ValueError("evidence must be AcceptanceEvidence")
    if not isinstance(path, Path) or not path.is_absolute():
        raise ValueError("evidence output path must be absolute")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            evidence.to_dict(),
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def load_acceptance_fixture(path: Path) -> AcceptanceFixture:
    """Load a bounded manifest for the full STAR+Salmon/SortMeRNA gate."""
    if not isinstance(path, Path) or not path.is_absolute():
        raise ValueError("fixture manifest path must be absolute")
    raw = _read_bounded_regular_file(
        path,
        maximum_bytes=_MAX_FIXTURE_MANIFEST_BYTES,
        error_message="fixture manifest is unavailable",
    )
    try:
        document = json.loads(raw, object_pairs_hook=_unique_object)
    except (UnicodeError, ValueError, json.JSONDecodeError):
        raise AssertionError("fixture manifest is invalid JSON") from None
    expected_keys = {
        "schema_version",
        "source_manifest_sha256",
        "source_identity_sha256",
        "index_provenance_manifest_sha256",
        "index_provenance_identity_sha256",
        "workflow_inputs",
        "transcriptome_binding",
        "required_artifact_output_types",
        "required_qc_metric_keys",
        "required_sample_ids",
        "required_artifact_sample_output_types",
        "required_qc_sample_metric_keys",
        "required_qc_sample_metric_values",
    }
    if not isinstance(document, dict) or set(document) != expected_keys:
        raise AssertionError("fixture manifest fields are invalid")
    if document["schema_version"] != FIXTURE_MANIFEST_SCHEMA_VERSION:
        raise AssertionError("fixture manifest schema is unsupported")
    inputs_document = document["workflow_inputs"]
    if not isinstance(inputs_document, dict) or set(inputs_document) != {
        "config",
        "samples",
        "options",
    }:
        raise AssertionError("fixture workflow inputs are incomplete")
    if (
        not isinstance(inputs_document["config"], dict)
        or not isinstance(inputs_document["samples"], list)
        or not isinstance(inputs_document["options"], dict)
    ):
        raise AssertionError("fixture workflow input types are invalid")
    _require_full_fixture_routes(inputs_document)
    (
        source_manifest_sha256,
        source_identity_sha256,
        provenance_manifest_sha256,
        provenance_identity_sha256,
    ) = _verify_fixture_identity_closure(
        manifest_path=path,
        acceptance=document,
        workflow_inputs=inputs_document,
    )
    transcriptome_document = document["transcriptome_binding"]
    if not isinstance(transcriptome_document, dict) or set(transcriptome_document) != {
        "reference_id",
        "fasta_sha256",
        "gtf_sha256",
        "transcript_fasta",
        "transcript_fasta_sha256",
    }:
        raise AssertionError("fixture transcriptome binding is incomplete")
    try:
        transcriptome = BulkRnaSeqTranscriptomeBinding(
            reference_id=transcriptome_document["reference_id"],
            fasta_sha256=transcriptome_document["fasta_sha256"],
            gtf_sha256=transcriptome_document["gtf_sha256"],
            transcript_fasta=Path(transcriptome_document["transcript_fasta"]),
            transcript_fasta_sha256=transcriptome_document["transcript_fasta_sha256"],
        )
    except (TypeError, ValueError):
        raise AssertionError("fixture transcriptome binding is invalid") from None
    artifact_types = _manifest_tokens(
        document["required_artifact_output_types"],
        pattern=_OUTPUT_TYPE,
        name="required artifact output types",
    )
    metric_keys = _manifest_tokens(
        document["required_qc_metric_keys"],
        pattern=_METRIC_KEY,
        name="required QC metric keys",
    )
    sample_ids = _manifest_tokens(
        document["required_sample_ids"],
        pattern=_PUBLIC_ID,
        name="required sample IDs",
    )
    artifact_coordinates = _manifest_pairs(
        document["required_artifact_sample_output_types"],
        first_pattern=_PUBLIC_ID,
        second_pattern=_OUTPUT_TYPE,
        name="required artifact sample/output coordinates",
    )
    metric_coordinates = _manifest_pairs(
        document["required_qc_sample_metric_keys"],
        first_pattern=_PUBLIC_ID,
        second_pattern=_METRIC_KEY,
        name="required QC sample/metric coordinates",
    )
    metric_values = _manifest_metric_values(
        document["required_qc_sample_metric_values"],
        name="required QC sample/metric values",
    )
    if metric_values != _EXPECTED_RAW_FASTQC_VALUES:
        raise AssertionError("fixture raw FastQC oracle differs from its fixed source")
    observed_samples = {
        row["sample"]
        for row in inputs_document["samples"]
        if isinstance(row, dict) and isinstance(row.get("sample"), str)
    }
    if set(sample_ids) != observed_samples:
        raise AssertionError("required sample IDs must cover the fixture exactly")
    if (
        {sample for sample, _output_type in artifact_coordinates} != observed_samples
        or {sample for sample, _metric_key in metric_coordinates} != observed_samples
        or {output_type for _sample, output_type in artifact_coordinates}
        - set(artifact_types)
        or {metric_key for _sample, metric_key in metric_coordinates} - set(metric_keys)
        or {(sample, metric_key) for sample, metric_key, _value in metric_values}
        - set(metric_coordinates)
    ):
        raise AssertionError("fixture sample evidence matrix is inconsistent")
    _require_sample_evidence_routes(
        inputs_document,
        artifact_coordinates=artifact_coordinates,
        metric_coordinates=metric_coordinates,
    )
    for prefix, values in (
        ("bulk_rnaseq.star.", artifact_types),
        ("bulk_rnaseq.salmon.", artifact_types),
        ("bulk_rnaseq.fastqc.", artifact_types),
        ("bulk_rnaseq.rrna.sortmerna.", artifact_types),
        ("star.", metric_keys),
        ("salmon.", metric_keys),
        ("fastqc.", metric_keys),
        ("trimming.", metric_keys),
    ):
        if not any(value.startswith(prefix) for value in values):
            raise AssertionError(f"fixture evidence must include {prefix} outputs")
    return AcceptanceFixture(
        workflow_inputs=WorkflowInputs(
            config=inputs_document["config"],
            samples=inputs_document["samples"],
            options=inputs_document["options"],
        ),
        transcriptome=transcriptome,
        acceptance_manifest_sha256=sha256(raw).hexdigest(),
        source_manifest_sha256=source_manifest_sha256,
        source_identity_sha256=source_identity_sha256,
        index_provenance_manifest_sha256=provenance_manifest_sha256,
        index_provenance_identity_sha256=provenance_identity_sha256,
        required_artifact_output_types=artifact_types,
        required_qc_metric_keys=metric_keys,
        required_sample_ids=sample_ids,
        required_artifact_sample_output_types=artifact_coordinates,
        required_qc_sample_metric_keys=metric_coordinates,
        required_qc_sample_metric_values=metric_values,
    )


def _verify_fixture_identity_closure(
    *,
    manifest_path: Path,
    acceptance: Mapping[str, object],
    workflow_inputs: Mapping[str, object],
) -> tuple[str, str, str, str]:
    """Bind the acceptance oracle to its path-free source/index provenance."""
    try:
        source_manifest_sha256 = _require_digest(
            acceptance["source_manifest_sha256"],
            "fixture source manifest",
        )
        source_identity_sha256 = _require_digest(
            acceptance["source_identity_sha256"],
            "fixture source identity",
        )
        provenance_manifest_sha256 = _require_digest(
            acceptance["index_provenance_manifest_sha256"],
            "fixture index provenance manifest",
        )
        provenance_identity_sha256 = _require_digest(
            acceptance["index_provenance_identity_sha256"],
            "fixture index provenance identity",
        )
    except (KeyError, ValueError):
        raise AssertionError("fixture identity coordinates are invalid") from None

    fixture_root = manifest_path.parent
    if fixture_root.resolve(strict=False) != fixture_root:
        raise AssertionError("fixture root must be canonical")
    source_raw = _read_bounded_regular_file(
        fixture_root / SOURCE_MANIFEST_FILENAME,
        maximum_bytes=_MAX_FIXTURE_MANIFEST_BYTES,
        error_message="fixture source manifest is unavailable",
    )
    provenance_raw = _read_bounded_regular_file(
        fixture_root / INDEX_PROVENANCE_FILENAME,
        maximum_bytes=_MAX_FIXTURE_MANIFEST_BYTES,
        error_message="fixture index provenance is unavailable",
    )
    if sha256(source_raw).hexdigest() != source_manifest_sha256:
        raise AssertionError("fixture source manifest digest changed")
    if sha256(provenance_raw).hexdigest() != provenance_manifest_sha256:
        raise AssertionError("fixture index provenance digest changed")
    source = _fixture_json_object(source_raw, "fixture source manifest")
    provenance = _fixture_json_object(provenance_raw, "fixture index provenance")
    _verify_source_manifest(
        source,
        fixture_root=fixture_root,
        acceptance=acceptance,
        workflow_inputs=workflow_inputs,
        expected_identity=source_identity_sha256,
    )
    _verify_index_provenance(
        provenance,
        fixture_root=fixture_root,
        acceptance=acceptance,
        workflow_inputs=workflow_inputs,
        expected_source_manifest_sha256=source_manifest_sha256,
        expected_source_identity_sha256=source_identity_sha256,
        expected_identity=provenance_identity_sha256,
    )
    return (
        source_manifest_sha256,
        source_identity_sha256,
        provenance_manifest_sha256,
        provenance_identity_sha256,
    )


def _fixture_json_object(raw: bytes, name: str) -> dict[str, object]:
    try:
        value = json.loads(raw, object_pairs_hook=_unique_object)
    except (UnicodeError, ValueError, json.JSONDecodeError):
        raise AssertionError(f"{name} is invalid JSON") from None
    if not isinstance(value, dict):
        raise AssertionError(f"{name} must be an object")
    return value


def _verify_source_manifest(
    source: Mapping[str, object],
    *,
    fixture_root: Path,
    acceptance: Mapping[str, object],
    workflow_inputs: Mapping[str, object],
    expected_identity: str,
) -> None:
    expected_keys = {
        "schema_version",
        "fixture_id",
        "biological_validity",
        "limitations",
        "generator",
        "reference",
        "samples",
        "files",
        "index_build_contracts",
        "source_identity_sha256",
    }
    if (
        set(source) != expected_keys
        or source.get("schema_version") != SOURCE_MANIFEST_SCHEMA_VERSION
        or source.get("fixture_id") != FIXTURE_ID
        or source.get("biological_validity") is not False
        or source.get("source_identity_sha256") != expected_identity
    ):
        raise AssertionError("fixture source identity is invalid")
    identity_payload = dict(source)
    identity_payload.pop("source_identity_sha256")
    if _framed_json_identity(SOURCE_IDENTITY_SCHEME, identity_payload) != (
        expected_identity
    ):
        raise AssertionError("fixture source identity closure changed")

    source_samples = source.get("samples")
    acceptance_samples = workflow_inputs.get("samples")
    if not isinstance(source_samples, list) or not isinstance(acceptance_samples, list):
        raise AssertionError("fixture source samples are invalid")
    expanded_samples: list[dict[str, object]] = []
    for value in source_samples:
        if not isinstance(value, dict):
            raise AssertionError("fixture source samples are invalid")
        row = dict(value)
        for key in ("fastq_1", "fastq_2"):
            if key not in row:
                continue
            relative = _fixture_relative_path(row[key], "fixture FASTQ")
            row[key] = str(fixture_root.joinpath(*relative.parts))
        expanded_samples.append(row)
    if expanded_samples != acceptance_samples:
        raise AssertionError("fixture source samples differ from accepted inputs")

    files = source.get("files")
    if not isinstance(files, list) or not files:
        raise AssertionError("fixture source file closure is invalid")
    observed_paths: list[str] = []
    total_bytes = 0
    file_digests: dict[str, str] = {}
    for value in files:
        if not isinstance(value, dict) or set(value) != {
            "path",
            "size_bytes",
            "sha256",
        }:
            raise AssertionError("fixture source file closure is invalid")
        relative = _fixture_relative_path(value["path"], "fixture source file")
        relative_text = relative.as_posix()
        expected_size = value["size_bytes"]
        expected_digest = value["sha256"]
        if (
            isinstance(expected_size, bool)
            or not isinstance(expected_size, int)
            or expected_size <= 0
            or expected_size > _MAX_FIXTURE_SOURCE_FILE_BYTES
        ):
            raise AssertionError("fixture source file size is invalid")
        try:
            _require_digest(expected_digest, "fixture source file")
        except ValueError:
            raise AssertionError("fixture source file digest is invalid") from None
        raw = _read_bounded_regular_file(
            fixture_root.joinpath(*relative.parts),
            maximum_bytes=_MAX_FIXTURE_SOURCE_FILE_BYTES,
            error_message="fixture source file is unavailable",
        )
        total_bytes += len(raw)
        if (
            len(raw) != expected_size
            or sha256(raw).hexdigest() != expected_digest
            or total_bytes > _MAX_FIXTURE_SOURCE_TOTAL_BYTES
        ):
            raise AssertionError("fixture source file closure changed")
        observed_paths.append(relative_text)
        file_digests[relative_text] = expected_digest
    if observed_paths != sorted(observed_paths) or len(observed_paths) != len(
        set(observed_paths)
    ):
        raise AssertionError("fixture source file closure is not canonical")
    _verify_source_bindings(
        fixture_root=fixture_root,
        acceptance=acceptance,
        workflow_inputs=workflow_inputs,
        acceptance_samples=acceptance_samples,
        file_digests=file_digests,
    )


def _verify_source_bindings(
    *,
    fixture_root: Path,
    acceptance: Mapping[str, object],
    workflow_inputs: Mapping[str, object],
    acceptance_samples: list[object],
    file_digests: Mapping[str, str],
) -> None:
    for row in acceptance_samples:
        if not isinstance(row, dict):
            raise AssertionError("fixture accepted samples are invalid")
        for key in ("fastq_1", "fastq_2"):
            if key not in row:
                continue
            relative = _fixture_descendant(
                row[key],
                fixture_root=fixture_root,
                name="accepted FASTQ",
            )
            if relative.as_posix() not in file_digests:
                raise AssertionError("accepted FASTQ is outside the source closure")
    config = workflow_inputs.get("config")
    try:
        standard = config["standard"]  # type: ignore[index]
        reference = standard["reference"]
        rrna = standard["ribosomal_rna_removal"]
        bindings = (
            (reference["fasta"], reference["fasta_sha256"]),
            (reference["gtf"], reference["gtf_sha256"]),
            (
                rrna["database_manifest"]["path"],
                rrna["database_manifest"]["identity_sha256"],
            ),
        )
        transcriptome = acceptance["transcriptome_binding"]
        bindings = (
            *bindings,
            (
                transcriptome["transcript_fasta"],
                transcriptome["transcript_fasta_sha256"],
            ),
        )
    except (KeyError, TypeError):
        raise AssertionError("fixture source bindings are incomplete") from None
    for raw_path, expected_digest in bindings:
        relative = _fixture_descendant(
            raw_path,
            fixture_root=fixture_root,
            name="accepted source binding",
        )
        if file_digests.get(relative.as_posix()) != expected_digest:
            raise AssertionError("accepted source binding differs from its source")


def _verify_index_provenance(
    provenance: Mapping[str, object],
    *,
    fixture_root: Path,
    acceptance: Mapping[str, object],
    workflow_inputs: Mapping[str, object],
    expected_source_manifest_sha256: str,
    expected_source_identity_sha256: str,
    expected_identity: str,
) -> None:
    if (
        set(provenance)
        != {
            "schema_version",
            "fixture_id",
            "biological_validity",
            "source_manifest_sha256",
            "source_identity_sha256",
            "indexes",
            "provenance_identity_sha256",
        }
        or provenance.get("schema_version") != INDEX_PROVENANCE_SCHEMA_VERSION
        or provenance.get("fixture_id") != FIXTURE_ID
        or provenance.get("biological_validity") is not False
        or provenance.get("source_manifest_sha256") != expected_source_manifest_sha256
        or provenance.get("source_identity_sha256") != expected_source_identity_sha256
        or provenance.get("provenance_identity_sha256") != expected_identity
    ):
        raise AssertionError("fixture index provenance identity is invalid")
    identity_payload = dict(provenance)
    identity_payload.pop("provenance_identity_sha256")
    if _framed_json_identity(PROVENANCE_IDENTITY_SCHEME, identity_payload) != (
        expected_identity
    ):
        raise AssertionError("fixture index provenance closure changed")

    indexes = provenance.get("indexes")
    expected_routes = {
        "star": ("STAR_GENOMEGENERATE", "2.7.11b", "star_index"),
        "salmon": ("SALMON_INDEX", "1.10.3", "salmon_index"),
        "sortmerna": ("SORTMERNA", "4.3.7", "sortmerna_index"),
    }
    if not isinstance(indexes, dict) or set(indexes) != set(expected_routes):
        raise AssertionError("fixture index provenance routes are invalid")
    config = workflow_inputs.get("config")
    try:
        standard = config["standard"]  # type: ignore[index]
        reference = standard["reference"]
        rrna = standard["ribosomal_rna_removal"]
        accepted_indexes = {
            "star": reference["star_index"],
            "salmon": reference["salmon_index"],
            "sortmerna": rrna["sortmerna_index"],
        }
    except (KeyError, TypeError):
        raise AssertionError("fixture accepted index bindings are incomplete") from None
    expected_index_keys = {
        "relative_index_root",
        "producer_process",
        "tool",
        "tool_version",
        "container_image",
        "immutable_config_sha256",
        "command_argv",
        "command_argv_sha256",
        "sidecar_manifest_sha256",
        "index_closure_sha256",
    }
    for kind, (process, version, binding_name) in expected_routes.items():
        value = indexes[kind]
        accepted = accepted_indexes[kind]
        if not isinstance(value, dict) or set(value) != expected_index_keys:
            raise AssertionError("fixture index provenance entry is invalid")
        argv = value.get("command_argv")
        try:
            _require_digest(
                value.get("immutable_config_sha256"),
                "fixture index config",
            )
            command_digest = _require_digest(
                value.get("command_argv_sha256"),
                "fixture index command",
            )
            sidecar_digest = _require_digest(
                value.get("sidecar_manifest_sha256"),
                "fixture index sidecar",
            )
            _require_digest(
                value.get("index_closure_sha256"),
                "fixture index closure",
            )
        except ValueError:
            raise AssertionError("fixture index provenance digest is invalid") from None
        if (
            value.get("producer_process") != process
            or value.get("tool") != kind
            or value.get("tool_version") != version
            or not isinstance(value.get("container_image"), str)
            or _IMMUTABLE_IMAGE.fullmatch(value["container_image"]) is None
            or value.get("container_image") != _EXPECTED_INDEX_PRODUCER_IMAGES[kind]
            or not isinstance(argv, list)
            or not argv
            or sha256(_canonical_json_bytes(argv)).hexdigest() != command_digest
            or not isinstance(accepted, dict)
            or accepted.get("identity_sha256") != sidecar_digest
        ):
            raise AssertionError("fixture index provenance entry changed")
        _require_path_free_provenance_argv(argv, fixture_root=fixture_root)
        relative = _fixture_relative_path(
            value.get("relative_index_root"),
            "fixture index root",
        )
        accepted_relative = _fixture_descendant(
            accepted.get("path"),
            fixture_root=fixture_root,
            name="accepted index root",
        )
        if relative != accepted_relative:
            raise AssertionError("fixture index provenance binding is invalid")
    _verify_actual_index_closures(
        indexes=indexes,
        acceptance=acceptance,
        workflow_inputs=workflow_inputs,
    )


def _verify_actual_index_closures(
    *,
    indexes: Mapping[str, object],
    acceptance: Mapping[str, object],
    workflow_inputs: Mapping[str, object],
) -> None:
    try:
        config = workflow_inputs["config"]
        standard = config["standard"]  # type: ignore[index]
        reference = standard["reference"]
        rrna = standard["ribosomal_rna_removal"]
        transcriptome = acceptance["transcriptome_binding"]
        annotation_style = reference.get("annotation_style", "ensembl")
        fasta_sha256 = reference["fasta_sha256"]
        gtf_sha256 = reference["gtf_sha256"]
        transcript_fasta_sha256 = transcriptome["transcript_fasta_sha256"]
    except (AttributeError, KeyError, TypeError):
        raise AssertionError("fixture index closure inputs are incomplete") from None
    for kind in ("star", "salmon"):
        binding = reference[f"{kind}_index"]
        provenance = indexes[kind]
        if not isinstance(binding, dict) or not isinstance(provenance, dict):
            raise AssertionError("fixture reference index binding is invalid")
        result = verify_reference_index(
            binding.get("path"),
            kind=kind,
            expected_manifest_sha256=binding.get("identity_sha256"),
            fasta_sha256=fasta_sha256,
            gtf_sha256=gtf_sha256,
            transcript_fasta_sha256=(
                transcript_fasta_sha256 if kind == "salmon" else None
            ),
            annotation_style=annotation_style,
            expected_container_image=provenance["container_image"],
        )
        if (
            result.is_failure
            or result.value is None
            or result.value.identity_sha256 != provenance.get("index_closure_sha256")
        ):
            raise AssertionError("fixture reference index closure changed")

    database = rrna["database_manifest"]
    sortmerna_binding = rrna["sortmerna_index"]
    sortmerna_provenance = indexes["sortmerna"]
    if (
        not isinstance(database, dict)
        or not isinstance(sortmerna_binding, dict)
        or not isinstance(sortmerna_provenance, dict)
    ):
        raise AssertionError("fixture SortMeRNA index binding is invalid")
    database_result = verify_ribo_database_manifest(
        database.get("path"),
        expected_manifest_sha256=database.get("identity_sha256"),
    )
    if database_result.is_failure or database_result.value is None:
        raise AssertionError("fixture rRNA database closure changed")
    sortmerna_result = verify_sortmerna_index(
        sortmerna_binding.get("path"),
        expected_index_sha256=sortmerna_binding.get("identity_sha256"),
        database_closure=database_result.value,
    )
    if (
        sortmerna_result.is_failure
        or sortmerna_result.value is None
        or sortmerna_result.value.identity_sha256
        != sortmerna_provenance.get("index_closure_sha256")
    ):
        raise AssertionError("fixture SortMeRNA index closure changed")


def _fixture_relative_path(value: object, name: str) -> PurePosixPath:
    if not isinstance(value, str):
        raise AssertionError(f"{name} path is invalid")
    relative = PurePosixPath(value)
    if (
        relative.is_absolute()
        or not relative.parts
        or any(part in {"", ".", ".."} for part in relative.parts)
    ):
        raise AssertionError(f"{name} path is invalid")
    return relative


def _require_path_free_provenance_argv(
    value: object,
    *,
    fixture_root: Path,
) -> None:
    if not isinstance(value, list) or not value or len(value) > 512:
        raise AssertionError("fixture index provenance argv is invalid")
    total_bytes = 0
    for token in value:
        if not isinstance(token, str) or not token or len(token) > 4096:
            raise AssertionError("fixture index provenance argv is invalid")
        total_bytes += len(token.encode("utf-8"))
        if (
            total_bytes > 64 * 1024
            or str(fixture_root) in token
            or any(character in token for character in ("\x00", "\n", "\r"))
            or _WINDOWS_ABSOLUTE.match(token)
            or _EMBEDDED_ABSOLUTE.search(token)
            or token.startswith("file:")
        ):
            raise AssertionError("fixture index provenance argv is not path-free")


def _fixture_descendant(
    value: object,
    *,
    fixture_root: Path,
    name: str,
) -> PurePosixPath:
    if not isinstance(value, str):
        raise AssertionError(f"{name} path is invalid")
    path = Path(value)
    if not path.is_absolute() or path.resolve(strict=False) != path:
        raise AssertionError(f"{name} path is invalid")
    try:
        relative = path.relative_to(fixture_root)
    except ValueError:
        raise AssertionError(f"{name} is outside the fixture root") from None
    return _fixture_relative_path(relative.as_posix(), name)


def _canonical_json_bytes(value: object) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")


def _framed_json_identity(scheme: str, payload: Mapping[str, object]) -> str:
    digest = sha256()
    for value in (scheme.encode("ascii"), _canonical_json_bytes(payload)):
        digest.update(len(value).to_bytes(8, "big"))
        digest.update(value)
    return digest.hexdigest()


def _require_full_fixture_routes(inputs: Mapping[str, object]) -> None:
    samples = inputs["samples"]
    assert isinstance(samples, list)
    if not samples or any(not isinstance(row, dict) for row in samples):
        raise AssertionError("fixture samples must be non-empty objects")
    layouts = {row.get("layout") for row in samples}
    if layouts != {"SE", "PE"}:
        raise AssertionError("fixture must contain representative SE and PE samples")
    lanes: dict[str, set[str]] = {}
    for row in samples:
        sample = row.get("sample")
        lane = row.get("lane")
        if not isinstance(sample, str) or not isinstance(lane, str):
            raise AssertionError("fixture sample and lane identities are required")
        lanes.setdefault(sample, set()).add(lane)
    if not any(len(values) >= 2 for values in lanes.values()):
        raise AssertionError("fixture must contain a repeated lane sample")

    config = inputs["config"]
    assert isinstance(config, dict)
    advanced = config.get("advanced")
    if (
        not isinstance(advanced, dict)
        or set(advanced) != {"min_mapped_reads", "min_trimmed_reads"}
        or type(advanced.get("min_mapped_reads")) is not int
        or type(advanced.get("min_trimmed_reads")) is not int
        or advanced["min_mapped_reads"] != 0
        or advanced["min_trimmed_reads"] != 1
    ):
        raise AssertionError("fixture must bind the tiny execution thresholds")
    standard = config.get("standard")
    if not isinstance(standard, dict):
        raise AssertionError("fixture standard configuration is required")
    analysis = standard.get("analysis")
    if not isinstance(analysis, dict) or analysis != {
        "alignment": "star",
        "quantification": "salmon",
    }:
        raise AssertionError("fixture must use the canonical STAR+Salmon route")
    trimming = standard.get("trimming")
    if not isinstance(trimming, dict) or trimming.get("enabled") is not True:
        raise AssertionError("fixture must enable the standard trimming route")
    reference = standard.get("reference")
    if not isinstance(reference, dict) or not all(
        isinstance(reference.get(name), dict) for name in ("star_index", "salmon_index")
    ):
        raise AssertionError("fixture must bind STAR and Salmon indexes")
    rrna = standard.get("ribosomal_rna_removal")
    if (
        not isinstance(rrna, dict)
        or rrna.get("enabled") is not True
        or rrna.get("tool") != "sortmerna"
        or not isinstance(rrna.get("database_manifest"), dict)
        or not isinstance(rrna.get("sortmerna_index"), dict)
    ):
        raise AssertionError("fixture must bind a real SortMeRNA database and index")


def _manifest_tokens(
    value: object, *, pattern: re.Pattern[str], name: str
) -> tuple[str, ...]:
    if not isinstance(value, list):
        raise AssertionError(f"{name} must be an array")
    if any(not isinstance(item, str) for item in value):
        raise AssertionError(f"{name} are invalid")
    if (
        not value
        or len(value) != len(set(value))
        or any(pattern.fullmatch(item) is None for item in value)
    ):
        raise AssertionError(f"{name} are invalid") from None
    return tuple(sorted(value))


def _manifest_pairs(
    value: object,
    *,
    first_pattern: re.Pattern[str],
    second_pattern: re.Pattern[str],
    name: str,
) -> tuple[tuple[str, str], ...]:
    if not isinstance(value, list) or any(
        not isinstance(item, list) or len(item) != 2 for item in value
    ):
        raise AssertionError(f"{name} must be an array of pairs")
    try:
        return _closed_sorted_pairs(
            tuple(tuple(item) for item in value),
            first_pattern=first_pattern,
            second_pattern=second_pattern,
            name=name,
        )
    except ValueError:
        raise AssertionError(f"{name} are invalid") from None


def _manifest_metric_values(
    value: object,
    *,
    name: str,
) -> tuple[tuple[str, str, str], ...]:
    if not isinstance(value, list) or any(
        not isinstance(item, list) or len(item) != 3 for item in value
    ):
        raise AssertionError(f"{name} must be an array of triples")
    try:
        return _closed_sorted_metric_values(
            tuple(tuple(item) for item in value),
            name=name,
        )
    except ValueError:
        raise AssertionError(f"{name} are invalid") from None


def _require_sample_evidence_routes(
    inputs: Mapping[str, object],
    *,
    artifact_coordinates: tuple[tuple[str, str], ...],
    metric_coordinates: tuple[tuple[str, str], ...],
) -> None:
    samples = inputs["samples"]
    assert isinstance(samples, list)
    layouts = {
        row["sample"]: row["layout"]
        for row in samples
        if isinstance(row, dict)
        and isinstance(row.get("sample"), str)
        and isinstance(row.get("layout"), str)
    }
    artifact_set = set(artifact_coordinates)
    metric_set = set(metric_coordinates)
    common_artifacts = {
        "bulk_rnaseq.salmon.meta_info",
        "bulk_rnaseq.salmon.quant_gene",
        "bulk_rnaseq.star.bam",
        "bulk_rnaseq.star.log_final",
    }
    common_metrics = {
        "salmon.mapping_fraction",
        "salmon.processed_fragments",
        "star.input_templates",
        "star.uniquely_mapped_template_fraction",
    }
    for sample, layout in layouts.items():
        roles = ("single",) if layout == "SE" else ("read1", "read2")
        required_artifacts = common_artifacts | {
            *(f"bulk_rnaseq.fastqc.raw.{role}.zip" for role in roles),
            *(f"bulk_rnaseq.rrna.sortmerna.filtered.{role}" for role in roles),
        }
        required_metrics = common_metrics | {
            *(f"fastqc.raw.{role}.total_sequences" for role in roles),
            *(f"trimming.{role}.retained_reads" for role in roles),
        }
        if any(
            (sample, output_type) not in artifact_set
            for output_type in required_artifacts
        ):
            raise AssertionError(
                "fixture artifact sample matrix omits a required route"
            )
        if any(
            (sample, metric_key) not in metric_set for metric_key in required_metrics
        ):
            raise AssertionError("fixture QC sample matrix omits a required route")


def _closed_sorted_tokens(
    value: object,
    *,
    pattern: re.Pattern[str],
    name: str,
) -> tuple[str, ...]:
    if isinstance(value, str):
        raise ValueError(f"{name} must be a sequence")
    try:
        values = tuple(value)
    except TypeError:
        raise ValueError(f"{name} must be a sequence") from None
    if any(not isinstance(item, str) for item in values):
        raise ValueError(f"{name} are invalid")
    if (
        not values
        or values != tuple(sorted(values))
        or len(values) != len(set(values))
        or any(pattern.fullmatch(item) is None for item in values)
    ):
        raise ValueError(f"{name} are invalid")
    return values


def _closed_sorted_pairs(
    value: object,
    *,
    first_pattern: re.Pattern[str],
    second_pattern: re.Pattern[str],
    name: str,
) -> tuple[tuple[str, str], ...]:
    if isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be a sequence")
    try:
        pairs = tuple(tuple(item) for item in value)  # type: ignore[arg-type]
    except TypeError:
        raise ValueError(f"{name} must be a sequence") from None
    if (
        not pairs
        or pairs != tuple(sorted(pairs))
        or len(pairs) != len(set(pairs))
        or any(
            len(pair) != 2
            or not isinstance(pair[0], str)
            or not isinstance(pair[1], str)
            or first_pattern.fullmatch(pair[0]) is None
            or second_pattern.fullmatch(pair[1]) is None
            for pair in pairs
        )
    ):
        raise ValueError(f"{name} are invalid")
    return pairs


def _closed_sorted_metric_values(
    value: object,
    *,
    name: str,
) -> tuple[tuple[str, str, str], ...]:
    if isinstance(value, (str, bytes)):
        raise ValueError(f"{name} must be a sequence")
    try:
        values = tuple(tuple(item) for item in value)  # type: ignore[arg-type]
    except TypeError:
        raise ValueError(f"{name} must be a sequence") from None
    if not values or values != tuple(sorted(values)) or len(values) != len(set(values)):
        raise ValueError(f"{name} are invalid")
    for item in values:
        if (
            len(item) != 3
            or not isinstance(item[0], str)
            or _PUBLIC_ID.fullmatch(item[0]) is None
            or not isinstance(item[1], str)
            or _METRIC_KEY.fullmatch(item[1]) is None
            or not isinstance(item[2], str)
        ):
            raise ValueError(f"{name} are invalid")
        try:
            decimal_from_canonical_text(item[2])
        except ValueError:
            raise ValueError(f"{name} are invalid") from None
    if len({(item[0], item[1]) for item in values}) != len(values):
        raise ValueError(f"{name} are invalid")
    return values


def _unique_object(values: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in values:
        if key in result:
            raise ValueError("duplicate JSON key")
        result[key] = value
    return result


def _require_digest(value: object, name: str) -> str:
    if not isinstance(value, str) or _HEX64.fullmatch(value) is None:
        raise ValueError(f"{name} digest is invalid")
    return value


def _values_payload(values: AcceptanceEvidenceValues) -> dict[str, object]:
    payload = asdict(values)
    payload["artifact_content_sha256"] = [
        list(item) for item in values.artifact_content_sha256
    ]
    payload["artifact_output_types"] = list(values.artifact_output_types)
    payload["qc_metric_keys"] = list(values.qc_metric_keys)
    payload["qc_sample_ids"] = list(values.qc_sample_ids)
    payload["artifact_sample_output_types"] = [
        list(item) for item in values.artifact_sample_output_types
    ]
    payload["qc_sample_metric_keys"] = [
        list(item) for item in values.qc_sample_metric_keys
    ]
    payload["qc_sample_metric_values"] = [
        list(item) for item in values.qc_sample_metric_values
    ]
    return payload


def _evidence_digest(values: AcceptanceEvidenceValues) -> str:
    raw = json.dumps(
        {
            "schema_version": ACCEPTANCE_EVIDENCE_SCHEMA_VERSION,
            **_values_payload(values),
        },
        allow_nan=False,
        ensure_ascii=True,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return sha256(raw).hexdigest()


def _read_identity_document(path: Path) -> dict[str, object]:
    raw = _read_bounded_regular_file(
        path,
        maximum_bytes=_MAX_IDENTITY_DOCUMENT_BYTES,
        error_message="workspace identity document is unavailable",
    )
    try:
        value = json.loads(raw, object_pairs_hook=_unique_object)
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError):
        raise AssertionError("workspace identity document is invalid") from None
    if not isinstance(value, dict):
        raise AssertionError("workspace identity document is invalid")
    return value


def _read_bounded_regular_file(
    path: Path,
    *,
    maximum_bytes: int,
    error_message: str,
) -> bytes:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = -1
    try:
        descriptor = os.open(path, flags)
        info = os.fstat(descriptor)
        if (
            not stat.S_ISREG(info.st_mode)
            or info.st_size <= 0
            or info.st_size > maximum_bytes
        ):
            raise OSError
        chunks: list[bytes] = []
        remaining = maximum_bytes + 1
        while remaining:
            chunk = os.read(descriptor, min(64 * 1024, remaining))
            if not chunk:
                break
            chunks.append(chunk)
            remaining -= len(chunk)
        raw = b"".join(chunks)
        if len(raw) != info.st_size or len(raw) > maximum_bytes:
            raise OSError
        return raw
    except (OSError, TypeError, ValueError):
        raise AssertionError(error_message) from None
    finally:
        if descriptor >= 0:
            try:
                os.close(descriptor)
            except OSError:
                pass


def _verified_workspace_identity(
    *,
    cache_identity: Mapping[str, object],
    execution_identity: Mapping[str, object],
    workspace: Path,
    workflow_build_digest: str,
) -> str:
    cache_keys = {
        "schema_version",
        "adapter_version",
        "execution_mode",
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
    execution_keys = {
        "schema_version",
        "workflow_id",
        "adapter_version",
        "execution_mode",
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
    if set(cache_identity) != cache_keys or set(execution_identity) != execution_keys:
        raise AssertionError("workspace identity fields are invalid")
    cache_coordinates = {
        name: value
        for name, value in cache_identity.items()
        if name != "identity_sha256"
    }
    cache_digest = _identity_digest(cache_identity, "identity_sha256")
    digest_fields = (
        "workflow_build_sha256",
        "execution_implementation_manifest_sha256",
        "execution_implementation_aggregate_sha256",
        "input_closure_sha256",
        "normalized_inputs_sha256",
    )
    for name in digest_fields:
        _identity_digest(cache_identity, name)
    ribo_identity = cache_identity.get("ribo_database_closure_sha256")
    if ribo_identity is not None:
        try:
            _require_digest(ribo_identity, "rRNA database closure")
        except ValueError:
            raise AssertionError("workspace cache identity is invalid") from None
    if (
        cache_identity.get("schema_version") != WORKSPACE_SCHEMA_VERSION
        or not isinstance(cache_identity.get("adapter_version"), str)
        or not cache_identity["adapter_version"]
        or not isinstance(cache_identity.get("execution_mode"), str)
        or not cache_identity["execution_mode"]
        or cache_identity.get("workflow_build_sha256") != workflow_build_digest
        or cache_identity.get("nextflow_version") != NEXTFLOW_VERSION
        or cache_identity.get("resume_scope") != "single-run-workspace"
        or cache_identity.get("sortmerna_index_build_strategy") not in {None}
        or cache_digest != sha256(_workspace_json_bytes(cache_coordinates)).hexdigest()
    ):
        raise AssertionError("workspace cache identity closure is invalid")

    workspace_identity = _identity_digest(
        execution_identity,
        "workspace_identity_sha256",
    )
    expected_workspace_identity = managed_container_scope(workspace)
    for name in (
        "build_identity_sha256",
        "input_identity_sha256",
        "cache_identity_sha256",
        "execution_implementation_manifest_sha256",
        "execution_implementation_aggregate_sha256",
        "container_process_audit_sha256",
        "workspace_contract_sha256",
    ):
        _identity_digest(execution_identity, name)
    if (
        execution_identity.get("schema_version") != WORKSPACE_SCHEMA_VERSION
        or execution_identity.get("workflow_id") != "bulk-rnaseq"
        or execution_identity.get("adapter_version")
        != cache_identity.get("adapter_version")
        or execution_identity.get("execution_mode")
        != cache_identity.get("execution_mode")
        or execution_identity.get("build_identity_sha256") != workflow_build_digest
        or execution_identity.get("input_identity_sha256")
        != cache_identity.get("input_closure_sha256")
        or execution_identity.get("ribo_database_closure_sha256") != ribo_identity
        or execution_identity.get("sortmerna_index_build_strategy")
        != cache_identity.get("sortmerna_index_build_strategy")
        or execution_identity.get("cache_identity_sha256") != cache_digest
        or execution_identity.get("execution_implementation_manifest_sha256")
        != cache_identity.get("execution_implementation_manifest_sha256")
        or execution_identity.get("execution_implementation_aggregate_sha256")
        != cache_identity.get("execution_implementation_aggregate_sha256")
        or execution_identity.get("resume_enabled") is not False
        or workspace_identity != expected_workspace_identity
    ):
        raise AssertionError("workspace execution identity closure is inconsistent")
    return workspace_identity


def _workspace_json_bytes(value: object) -> bytes:
    return (
        json.dumps(
            value,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode()


def _identity_digest(document: Mapping[str, object], key: str) -> str:
    try:
        return _require_digest(document[key], key)
    except (KeyError, ValueError):
        raise AssertionError("workspace identity digest is invalid") from None


def _hash_indexed_artifacts(
    *,
    run_service: RunService,
    workspace_root: Path,
    run_id: str,
    artifact_generation: str,
    artifacts: tuple[Any, ...],
) -> tuple[tuple[str, str], ...]:
    downloader = ArtifactDownloadService(
        run_service=run_service,
        workspace_root=workspace_root,
    )
    total_bytes = 0
    hashes: list[tuple[str, str]] = []
    for artifact in sorted(artifacts, key=lambda value: value.artifact_id):
        prepared = downloader.prepare(
            run_id,
            artifact.artifact_id,
            expected_generation=artifact_generation,
            expected_revision=artifact.revision,
        )
        if prepared.is_failure or prepared.value is None:
            raise AssertionError(
                "indexed artifact cannot be reopened at its generation"
            )
        plan = prepared.value
        if plan.size_bytes > _MAX_EVIDENCE_ARTIFACT_BYTES:
            plan.close()
            raise AssertionError("tiny acceptance artifact exceeds the evidence bound")
        total_bytes += plan.size_bytes
        if total_bytes > _MAX_EVIDENCE_TOTAL_BYTES:
            plan.close()
            raise AssertionError("tiny acceptance artifacts exceed the evidence bound")
        digest = sha256()
        try:
            for chunk in plan.iter_bytes():
                digest.update(chunk)
        finally:
            plan.close()
        hashes.append((artifact.artifact_id, digest.hexdigest()))
    return tuple(hashes)
