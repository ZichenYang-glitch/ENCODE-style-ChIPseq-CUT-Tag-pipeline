"""Safe post-success indexing of adapter-owned numeric QC summaries."""

from __future__ import annotations

from collections.abc import Mapping
from decimal import Decimal
import os
from pathlib import Path, PurePosixPath
import re
import stat

from encode_pipeline.platform.adapters import (
    QC_SUMMARY_EXTRACT_CAPABILITY,
    ExtractedQcMetricCandidate,
    QcSourceArtifact,
    QcSourceDocument,
    QcSummaryExtractingAdapter,
    WorkflowInputs,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.runs import (
    RunArtifactRef,
    RunQcMetric,
    RunStatus,
    build_qc_metric_id,
    validate_qc_identifier_token,
)
from encode_pipeline.services.run_repositories import (
    ConcurrentRunUpdateError,
    canonical_decimal_text,
)
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider


_MAX_SOURCE_BYTES = 1024 * 1024
_SAFE_ARTIFACT_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,127}$")
_SAFE_RUN_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_SAFE_OUTPUT_TYPE = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,127}$")
_SAFE_METRIC_KEY = re.compile(r"^[a-z][a-z0-9_]*(?:\.[a-z][a-z0-9_]*)*$")
_ALLOWED_UNITS = frozenset({"count", "fraction", "ratio"})
_ALLOWED_SCOPES = frozenset({"run", "sample", "experiment"})
_ALLOWED_FLAGS = frozenset({"pass", "warning", "fail"})
_SOURCE_METADATA_KEYS = ("scope", "sample_id", "experiment_id", "assay")


class QcSummaryIndexingService:
    """Validate and persist a complete run-scoped numeric QC index."""

    def __init__(
        self,
        *,
        run_service: RunService,
        registry: WorkflowRegistry,
        build_identity_provider: WorkflowBuildIdentityProvider,
        workspace_root: Path,
    ) -> None:
        if not isinstance(run_service, RunService):
            raise ValueError("run_service must be a RunService")
        if not isinstance(registry, WorkflowRegistry):
            raise ValueError("registry must be a WorkflowRegistry")
        if not isinstance(build_identity_provider, WorkflowBuildIdentityProvider):
            raise ValueError("build_identity_provider is invalid")
        if not isinstance(workspace_root, Path) or not workspace_root.is_absolute():
            raise ValueError("workspace_root must be an absolute Path")
        self._run_service = run_service
        self._registry = registry
        self._build_identity_provider = build_identity_provider
        self._workspace_root = workspace_root

    def index(
        self,
        run_id: str,
        artifacts: tuple[RunArtifactRef, ...],
    ) -> Result[tuple[RunQcMetric, ...]]:
        """Index one SUCCEEDED run without changing its terminal lifecycle."""
        try:
            record = self._run_service.get_run(run_id)
        except KeyError:
            return self._result_failure("QC_INDEXING_RUN_NOT_FOUND")
        if record.status is not RunStatus.SUCCEEDED:
            return self._result_failure("QC_INDEXING_RUN_NOT_SUCCEEDED")

        try:
            adapter = self._registry.get(record.workflow_id)
            if (
                QC_SUMMARY_EXTRACT_CAPABILITY not in adapter.capabilities.supports
                or not isinstance(adapter, QcSummaryExtractingAdapter)
            ):
                return self._fail(run_id, "QC_INDEXING_UNSUPPORTED")
            if not self._build_matches(record.run_id, record.workflow_id):
                return self._fail(run_id, "QC_INDEXING_BUILD_MISMATCH")
            expected_artifacts = self._validated_artifact_generation(run_id, artifacts)
            if expected_artifacts != tuple(
                sorted(
                    self._run_service.list_artifacts(run_id),
                    key=lambda item: item.artifact_id,
                )
            ):
                return self._fail(
                    run_id,
                    "QC_INDEXING_ARTIFACT_GENERATION_MISMATCH",
                )
            source_types = self._validated_source_types(
                adapter.qc_source_output_types()
            )
            workspace = self._workspace_for_run(run_id)
            try:
                documents = self._source_documents(
                    expected_artifacts,
                    source_types,
                    workspace,
                )
            except (OSError, TypeError, ValueError):
                return self._fail(
                    run_id,
                    "QC_INDEXING_SOURCE_VALIDATION_FAILED",
                )
            inputs = self._reconstruct_inputs(record.inputs)
            candidate_result = adapter.extract_qc_metrics(inputs, documents)
            if candidate_result.is_failure:
                return self._fail(run_id, "QC_INDEXING_ADAPTER_FAILED")
            metrics = self._build_metrics(
                run_id,
                record.ended_at,
                candidate_result.value,
                {document.source.artifact_id for document in documents},
            )
            self._run_service.replace_qc_metrics(
                run_id,
                metrics,
                expected_artifacts=expected_artifacts,
            )
            return Result.success(metrics)
        except ConcurrentRunUpdateError:
            return self._fail(run_id, "QC_INDEXING_STATE_CHANGED")
        except (KeyError, OSError, TypeError, ValueError):
            return self._fail(run_id, "QC_INDEXING_VALIDATION_FAILED")
        except Exception:
            return self._fail(run_id, "QC_INDEXING_UNEXPECTED_ERROR")

    def record_unexpected_failure(self, run_id: str) -> None:
        """Record a best-effort public-safe worker defensive failure."""
        self._record_failure_event(run_id, "QC_INDEXING_UNEXPECTED_ERROR")

    def record_artifact_source_failure(self, run_id: str) -> None:
        """Record that a complete artifact generation was unavailable."""
        self._record_failure_event(run_id, "QC_INDEXING_ARTIFACTS_UNAVAILABLE")

    def _build_matches(self, run_id: str, workflow_id: str) -> bool:
        persisted = self._run_service.get_workflow_build_identity(run_id)
        if persisted is None:
            return False
        current = self._build_identity_provider.capture(workflow_id)
        return current.is_success and persisted.matches(current.value)

    @staticmethod
    def _reconstruct_inputs(snapshot: Mapping[str, object]) -> WorkflowInputs:
        config = snapshot.get("config")
        samples = snapshot.get("samples")
        options = snapshot.get("options", {})
        if not isinstance(config, Mapping) or not isinstance(options, Mapping):
            raise ValueError("run inputs are malformed")
        return WorkflowInputs(config=config, samples=samples, options=options)

    @staticmethod
    def _validated_artifact_generation(
        run_id: str,
        artifacts: tuple[RunArtifactRef, ...],
    ) -> tuple[RunArtifactRef, ...]:
        if not isinstance(artifacts, tuple):
            raise ValueError("artifacts must be a tuple")
        seen: set[str] = set()
        for artifact in artifacts:
            if not isinstance(artifact, RunArtifactRef) or artifact.run_id != run_id:
                raise ValueError("artifact does not match the run")
            if artifact.artifact_id in seen:
                raise ValueError("artifact generation contains duplicate ids")
            seen.add(artifact.artifact_id)
        return tuple(sorted(artifacts, key=lambda item: item.artifact_id))

    @staticmethod
    def _validated_source_types(values) -> frozenset[str]:
        if not isinstance(values, tuple):
            raise ValueError("QC source output types must be a tuple")
        seen: set[str] = set()
        for value in values:
            if (
                not isinstance(value, str)
                or _SAFE_OUTPUT_TYPE.fullmatch(value) is None
                or value in seen
            ):
                raise ValueError("QC source output type is invalid")
            seen.add(value)
        return frozenset(seen)

    def _workspace_for_run(self, run_id: str) -> Path:
        if (
            not isinstance(run_id, str)
            or _SAFE_RUN_ID.fullmatch(run_id) is None
            or run_id in {".", ".."}
        ):
            raise ValueError("run_id is not workspace-safe")
        return self._workspace_root / run_id

    def _source_documents(
        self,
        artifacts: tuple[RunArtifactRef, ...],
        source_types: frozenset[str],
        workspace: Path,
    ) -> tuple[QcSourceDocument, ...]:
        documents: list[QcSourceDocument] = []
        seen_paths: set[str] = set()
        for artifact in artifacts:
            output_type = artifact.metadata.get("output_type")
            if output_type not in source_types:
                continue
            if artifact.artifact_type != "file":
                raise ValueError("QC source artifact is not a file")
            if (
                not isinstance(output_type, str)
                or _SAFE_OUTPUT_TYPE.fullmatch(output_type) is None
                or _SAFE_ARTIFACT_ID.fullmatch(artifact.artifact_id) is None
            ):
                raise ValueError("QC source identity is invalid")
            relative_path = self._validated_relative_path(
                artifact.metadata.get("relative_path")
            )
            if relative_path in seen_paths:
                raise ValueError("QC source path is duplicated")
            seen_paths.add(relative_path)
            expected_size = self._validated_source_size(
                artifact.metadata.get("size_bytes")
            )
            metadata = self._validated_source_metadata(artifact.metadata)
            content = self._read_bounded_regular_file(
                workspace,
                relative_path,
                expected_size,
            )
            documents.append(
                QcSourceDocument(
                    source=QcSourceArtifact(
                        artifact_id=artifact.artifact_id,
                        output_type=output_type,
                        relative_path=relative_path,
                        metadata=metadata,
                    ),
                    content=content,
                )
            )
        return tuple(documents)

    @staticmethod
    def _validated_relative_path(value: object) -> str:
        if (
            not isinstance(value, str)
            or not value
            or len(value) > 2048
            or "\x00" in value
            or "\\" in value
            or value.startswith(("/", "~"))
            or (len(value) > 1 and value[1] == ":")
        ):
            raise ValueError("QC source path is invalid")
        path = PurePosixPath(value)
        if (
            path.as_posix() != value
            or path.is_absolute()
            or not path.parts
            or path.parts[0] != "results"
            or any(part in {"", ".", ".."} for part in value.split("/"))
        ):
            raise ValueError("QC source path is invalid")
        return value

    @staticmethod
    def _validated_source_metadata(metadata: Mapping[str, object]) -> dict[str, str]:
        if not isinstance(metadata, Mapping):
            raise ValueError("QC source metadata is invalid")
        result: dict[str, str] = {}
        for key in _SOURCE_METADATA_KEYS:
            value = metadata.get(key)
            if value is None:
                continue
            try:
                validated = validate_qc_identifier_token(value)
            except ValueError:
                raise ValueError("QC source metadata is invalid")
            result[key] = validated
        return result

    @staticmethod
    def _validated_source_size(value: object) -> int:
        if (
            not isinstance(value, int)
            or isinstance(value, bool)
            or value < 0
            or value > _MAX_SOURCE_BYTES
        ):
            raise ValueError("QC source size is invalid")
        return value

    @staticmethod
    def _read_bounded_regular_file(
        workspace: Path,
        relative_path: str,
        expected_size: int,
    ) -> bytes:
        required_flags = ("O_DIRECTORY", "O_NOFOLLOW", "O_CLOEXEC", "O_NONBLOCK")
        if any(not hasattr(os, name) for name in required_flags):
            raise OSError("safe descriptor flags are unavailable")
        directory_flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
        final_flags = os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC
        components = (*workspace.parts[1:], *PurePosixPath(relative_path).parts)
        if not components:
            raise ValueError("QC source path is empty")
        descriptors: list[int] = []
        try:
            current = os.open("/", directory_flags)
            descriptors.append(current)
            for index, component in enumerate(components):
                final = index == len(components) - 1
                descriptor = os.open(
                    component,
                    final_flags if final else directory_flags,
                    dir_fd=current,
                )
                descriptors.append(descriptor)
                current = descriptor
                info = os.fstat(descriptor)
                if not final and not stat.S_ISDIR(info.st_mode):
                    raise ValueError("QC source parent is not a directory")
            before = os.fstat(current)
            if (
                not stat.S_ISREG(before.st_mode)
                or before.st_size != expected_size
                or before.st_size > _MAX_SOURCE_BYTES
            ):
                raise ValueError("QC source is not a bounded regular file")
            chunks: list[bytes] = []
            remaining = expected_size + 1
            while remaining:
                chunk = os.read(current, min(65536, remaining))
                if not chunk:
                    break
                chunks.append(chunk)
                remaining -= len(chunk)
            content = b"".join(chunks)
            after = os.fstat(current)
            identity_before = (
                before.st_dev,
                before.st_ino,
                before.st_mode,
                before.st_size,
                before.st_mtime_ns,
            )
            identity_after = (
                after.st_dev,
                after.st_ino,
                after.st_mode,
                after.st_size,
                after.st_mtime_ns,
            )
            if (
                len(content) > _MAX_SOURCE_BYTES
                or len(content) != expected_size
                or len(content) != before.st_size
                or identity_before != identity_after
            ):
                raise ValueError("QC source changed while it was read")
            return content
        finally:
            for descriptor in reversed(descriptors):
                try:
                    os.close(descriptor)
                except OSError:
                    pass

    def _build_metrics(
        self,
        run_id: str,
        ended_at,
        candidates,
        source_artifact_ids: set[str],
    ) -> tuple[RunQcMetric, ...]:
        if ended_at is None or not isinstance(candidates, tuple):
            raise ValueError("QC candidate collection is invalid")
        metrics: list[RunQcMetric] = []
        seen_ids: set[str] = set()
        for candidate in candidates:
            if not isinstance(candidate, ExtractedQcMetricCandidate):
                raise ValueError("QC candidate type is invalid")
            self._validate_candidate(candidate, source_artifact_ids)
            metric_id = build_qc_metric_id(
                candidate.metric_key,
                candidate.scope,
                candidate.sample_id,
                candidate.experiment_id,
            )
            if metric_id in seen_ids:
                raise ValueError("QC semantic metric is duplicated")
            seen_ids.add(metric_id)
            metrics.append(
                RunQcMetric(
                    metric_id=metric_id,
                    run_id=run_id,
                    metric_key=candidate.metric_key,
                    display_name=candidate.display_name,
                    value=candidate.value,
                    unit=candidate.unit,
                    scope=candidate.scope,
                    sample_id=candidate.sample_id,
                    experiment_id=candidate.experiment_id,
                    assay=candidate.assay,
                    qc_flag=candidate.qc_flag,
                    source_artifact_id=candidate.source_artifact_id,
                    produced_at=ended_at,
                )
            )
        return tuple(sorted(metrics, key=lambda metric: metric.metric_id))

    @staticmethod
    def _validate_candidate(
        candidate: ExtractedQcMetricCandidate,
        source_artifact_ids: set[str],
    ) -> None:
        if (
            not isinstance(candidate.metric_key, str)
            or len(candidate.metric_key) > 128
            or _SAFE_METRIC_KEY.fullmatch(candidate.metric_key) is None
        ):
            raise ValueError("QC metric key is invalid")
        if (
            not isinstance(candidate.display_name, str)
            or not 1 <= len(candidate.display_name) <= 255
            or not candidate.display_name.isprintable()
            or "/" in candidate.display_name
            or "\\" in candidate.display_name
        ):
            raise ValueError("QC metric display name is invalid")
        if not isinstance(candidate.value, Decimal):
            raise ValueError("QC metric value is invalid")
        canonical_decimal_text(candidate.value)
        if candidate.unit not in _ALLOWED_UNITS:
            raise ValueError("QC metric unit is invalid")
        if candidate.scope not in _ALLOWED_SCOPES:
            raise ValueError("QC metric scope is invalid")
        for value in (candidate.sample_id, candidate.experiment_id, candidate.assay):
            if value is not None:
                try:
                    validate_qc_identifier_token(value)
                except ValueError:
                    raise ValueError("QC metric identifier is invalid") from None
        if candidate.scope == "run" and (
            candidate.sample_id is not None or candidate.experiment_id is not None
        ):
            raise ValueError("run QC scope cannot have sample or experiment IDs")
        if candidate.scope == "sample" and candidate.sample_id is None:
            raise ValueError("sample QC scope requires sample_id")
        if candidate.scope == "experiment" and (
            candidate.sample_id is not None or candidate.experiment_id is None
        ):
            raise ValueError("experiment QC scope requires only experiment_id")
        if candidate.qc_flag is not None and candidate.qc_flag not in _ALLOWED_FLAGS:
            raise ValueError("QC metric flag is invalid")
        if candidate.source_artifact_id not in source_artifact_ids:
            raise ValueError("QC metric source artifact is invalid")

    def _fail(
        self,
        run_id: str,
        reason_code: str,
    ) -> Result[tuple[RunQcMetric, ...]]:
        self._record_failure_event(run_id, reason_code)
        return self._result_failure(reason_code)

    def _record_failure_event(self, run_id: str, reason_code: str) -> None:
        try:
            if self._run_service.get_run(run_id).status is RunStatus.SUCCEEDED:
                self._run_service.record_qc_metrics_failure(
                    run_id,
                    reason_code=reason_code,
                )
        except Exception:
            return

    @staticmethod
    def _result_failure(
        reason_code: str,
    ) -> Result[tuple[RunQcMetric, ...]]:
        return Result.failure(
            [
                Issue(
                    code="RUN_QC_INDEXING_FAILED",
                    message="Workflow QC metrics could not be indexed.",
                    severity="error",
                    path="qc_metrics",
                    source="qc_summary_indexing_service",
                    context={"reason_code": reason_code},
                )
            ]
        )
