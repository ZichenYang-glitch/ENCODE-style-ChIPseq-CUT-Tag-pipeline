"""Safe post-success extraction of adapter-owned artifact candidates."""

from __future__ import annotations

from collections.abc import Mapping
from hashlib import sha256
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat
from typing import cast
from urllib.parse import quote

from encode_pipeline.platform.adapters import (
    ARTIFACT_EXTRACT_CAPABILITY,
    ExtractedArtifactCandidate,
    MAX_SAMPLE_ROWS,
    QC_SUMMARY_EXTRACT_CAPABILITY,
    QcSummaryExtractingAdapter,
    SamplePayload,
    WorkflowAdapter,
    WorkflowInputs,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.result_generations import (
    build_artifact_content_revision,
    build_artifact_descriptor_revision,
)
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.runs import RunArtifactRef, RunStatus
from encode_pipeline.services.run_repositories import ConcurrentRunUpdateError
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider


_RESERVED_METADATA = frozenset({"relative_path", "output_type", "size_bytes"})
_SAFE_METADATA_KEY = re.compile(r"^[A-Za-z][A-Za-z0-9_]{0,63}$")
_SAFE_OUTPUT_TYPE = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,127}$")
_SAFE_MIME_TYPE = re.compile(
    r"^[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]{0,126}/"
    r"[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]{0,126}$"
)
_MAX_ARTIFACT_CANDIDATES = MAX_SAMPLE_ROWS * 128 + 128
_MAX_ARTIFACT_PATH_COMPONENTS = 32
_MAX_ARTIFACT_BYTES = 1024**4
_MAX_QC_REVISION_SOURCE_BYTES = 16 * 1024 * 1024
_MAX_TOTAL_QC_REVISION_BYTES = 256 * 1024 * 1024


class ArtifactExtractionService:
    """Validate and persist a complete run-scoped artifact index."""

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
        if (
            not isinstance(workspace_root, Path)
            or not workspace_root.is_absolute()
            or any(part in {"", ".", ".."} for part in workspace_root.parts[1:])
        ):
            raise ValueError("workspace_root must be an absolute Path")
        self._run_service = run_service
        self._registry = registry
        self._build_identity_provider = build_identity_provider
        self._workspace_root = workspace_root

    def begin_attempt(self, run_id: str) -> str:
        """Begin and return one caller-owned artifact indexing attempt."""
        state = self._run_service.begin_artifact_result_attempt(run_id)
        if state.artifact_attempt_id is None:
            raise ValueError("artifact attempt identity is unavailable")
        return state.artifact_attempt_id

    def extract(
        self,
        run_id: str,
        *,
        attempt_id: str | None = None,
    ) -> Result[tuple[RunArtifactRef, ...]]:
        """Index one canonical SUCCEEDED run without changing its lifecycle."""
        try:
            record = self._run_service.get_run(run_id)
        except KeyError:
            return self._result_failure("ARTIFACT_EXTRACTION_RUN_NOT_FOUND")
        if record.status is not RunStatus.SUCCEEDED:
            return self._result_failure("ARTIFACT_EXTRACTION_RUN_NOT_SUCCEEDED")

        try:
            if attempt_id is None:
                attempt_id = self.begin_attempt(run_id)
            else:
                self._run_service.begin_artifact_result_attempt(
                    run_id,
                    attempt_id=attempt_id,
                )
            adapter = self._registry.get(record.workflow_id)
            if ARTIFACT_EXTRACT_CAPABILITY not in adapter.capabilities.supports:
                return self._fail(
                    run_id,
                    "ARTIFACT_EXTRACTION_UNSUPPORTED",
                    attempt_id,
                )
            if not self._build_matches(record.run_id, record.workflow_id):
                return self._fail(
                    run_id,
                    "ARTIFACT_EXTRACTION_BUILD_MISMATCH",
                    attempt_id,
                )
            qc_source_types = self._qc_source_types(adapter)
            inputs = self._reconstruct_inputs(record.inputs)
            workspace = self._workspace_for_run(run_id)
            candidates_result = adapter.extract_artifacts(inputs, workspace)
            if candidates_result.is_failure:
                return self._fail(
                    run_id,
                    "ARTIFACT_EXTRACTION_ADAPTER_FAILED",
                    attempt_id,
                )
            references = self._build_references(
                record.run_id,
                record.ended_at,
                workspace,
                candidates_result.value,
                qc_source_types=qc_source_types,
            )
            self._run_service.replace_artifacts(
                run_id,
                references,
                attempt_id=attempt_id,
            )
            return Result.success(references)
        except ConcurrentRunUpdateError:
            return self._fail(
                run_id,
                "ARTIFACT_EXTRACTION_STATE_CHANGED",
                attempt_id,
            )
        except (KeyError, OSError, TypeError, ValueError):
            return self._fail(
                run_id,
                "ARTIFACT_EXTRACTION_VALIDATION_FAILED",
                attempt_id,
            )
        except Exception:
            return self._fail(
                run_id,
                "ARTIFACT_EXTRACTION_UNEXPECTED_ERROR",
                attempt_id,
            )

    def record_unexpected_failure(self, run_id: str, *, attempt_id: str) -> None:
        """Best-effort safe event for a worker-side defensive catch."""
        self._record_failure_event(
            run_id,
            "ARTIFACT_EXTRACTION_UNEXPECTED_ERROR",
            attempt_id,
        )

    @staticmethod
    def _qc_source_types(adapter: WorkflowAdapter) -> frozenset[str]:
        if QC_SUMMARY_EXTRACT_CAPABILITY not in adapter.capabilities.supports:
            return frozenset()
        if not isinstance(adapter, QcSummaryExtractingAdapter):
            raise ValueError("QC adapter contract is incomplete")
        values = adapter.qc_source_output_types()
        if not isinstance(values, tuple):
            raise ValueError("QC source output types must be a tuple")
        result: set[str] = set()
        for value in values:
            if (
                not isinstance(value, str)
                or _SAFE_OUTPUT_TYPE.fullmatch(value) is None
                or value in result
            ):
                raise ValueError("QC source output type is invalid")
            result.add(value)
        return frozenset(result)

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
        return WorkflowInputs(
            config=config,
            samples=cast(SamplePayload, samples),
            options=options,
        )

    def _workspace_for_run(self, run_id: str) -> Path:
        if (
            not isinstance(run_id, str)
            or not run_id
            or "\x00" in run_id
            or "/" in run_id
            or "\\" in run_id
            or run_id in {".", ".."}
        ):
            raise ValueError("run_id is not workspace-safe")
        self._require_directory_without_symlink_components(self._workspace_root)
        workspace = self._workspace_root / run_id
        self._require_directory_without_symlink_components(workspace)
        return workspace

    @staticmethod
    def _require_directory_without_symlink_components(path: Path) -> None:
        current = Path(path.anchor)
        component_stat = os.lstat(current)
        for part in path.parts[1:]:
            current /= part
            component_stat = os.lstat(current)
            if stat.S_ISLNK(component_stat.st_mode):
                raise ValueError("workspace path contains a symlink")
        if not stat.S_ISDIR(component_stat.st_mode):
            raise ValueError("workspace path is not a directory")

    def _build_references(
        self,
        run_id: str,
        ended_at,
        workspace: Path,
        candidates,
        *,
        qc_source_types: frozenset[str],
    ) -> tuple[RunArtifactRef, ...]:
        if ended_at is None:
            raise ValueError("succeeded run has no ended_at")
        if (
            not isinstance(candidates, tuple)
            or len(candidates) > _MAX_ARTIFACT_CANDIDATES
        ):
            raise ValueError("adapter candidates must be a tuple")
        workspace_stat = os.lstat(workspace)
        if stat.S_ISLNK(workspace_stat.st_mode) or not stat.S_ISDIR(
            workspace_stat.st_mode
        ):
            raise ValueError("run workspace is invalid")
        references: list[RunArtifactRef] = []
        seen_paths: set[str] = set()
        seen_ids: set[str] = set()
        total_qc_source_bytes = 0
        for candidate in candidates:
            if not isinstance(candidate, ExtractedArtifactCandidate):
                raise ValueError("adapter candidate type is invalid")
            relative_path = self._validate_relative_path(candidate.relative_path)
            if relative_path in seen_paths:
                raise ValueError("duplicate artifact path")
            seen_paths.add(relative_path)
            content_bound = candidate.output_type in qc_source_types
            target, size_bytes, revision = self._safe_regular_file(
                workspace,
                relative_path,
                output_type=candidate.output_type,
                content_bound=content_bound,
            )
            if content_bound:
                total_qc_source_bytes += size_bytes
                if total_qc_source_bytes > _MAX_TOTAL_QC_REVISION_BYTES:
                    raise ValueError("QC source files exceed the total byte limit")
            artifact_id = self._artifact_id(candidate.output_type, relative_path)
            if artifact_id in seen_ids:
                raise ValueError("duplicate artifact identity")
            seen_ids.add(artifact_id)
            metadata = self._validated_metadata(
                candidate.metadata,
                workspace,
            )
            metadata.update(
                {
                    "relative_path": relative_path,
                    "output_type": candidate.output_type,
                    "size_bytes": size_bytes,
                }
            )
            mime_type = candidate.mime_type
            if mime_type is not None and _SAFE_MIME_TYPE.fullmatch(mime_type) is None:
                raise ValueError("artifact MIME type is invalid")
            name = target.name
            if not name or len(name) > 255:
                raise ValueError("artifact name is invalid")
            references.append(
                RunArtifactRef(
                    artifact_id=artifact_id,
                    run_id=run_id,
                    artifact_type="file",
                    name=name,
                    uri=(
                        f"run://runs/{quote(run_id, safe='')}/artifacts/{artifact_id}"
                    ),
                    mime_type=mime_type,
                    produced_at=ended_at,
                    revision=revision,
                    metadata=metadata,
                )
            )
        return tuple(
            sorted(
                references,
                key=lambda item: (
                    str(item.metadata["output_type"]),
                    str(item.metadata["relative_path"]),
                ),
            )
        )

    @staticmethod
    def _validate_relative_path(value: str) -> str:
        if (
            not isinstance(value, str)
            or not value
            or len(value) > 2048
            or "\x00" in value
            or "\\" in value
            or value.startswith(("/", "~"))
            or (len(value) > 1 and value[1] == ":")
        ):
            raise ValueError("artifact path is invalid")
        path = PurePosixPath(value)
        if (
            path.as_posix() != value
            or path.is_absolute()
            or not path.parts
            or path.parts[0] != "results"
            or len(path.parts) > _MAX_ARTIFACT_PATH_COMPONENTS
            or any(part in {"", ".", ".."} for part in value.split("/"))
        ):
            raise ValueError("artifact path is invalid")
        return value

    @staticmethod
    def _safe_regular_file(
        workspace: Path,
        relative_path: str,
        *,
        output_type: str,
        content_bound: bool,
    ) -> tuple[Path, int, str]:
        required_flags = (
            "O_DIRECTORY",
            "O_NOFOLLOW",
            "O_CLOEXEC",
            "O_NONBLOCK",
            "O_PATH",
        )
        if any(not hasattr(os, name) for name in required_flags):
            raise OSError("safe descriptor flags are unavailable")
        directory_flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
        final_flags = (
            os.O_RDONLY | os.O_NONBLOCK | os.O_NOFOLLOW | os.O_CLOEXEC
            if content_bound
            else os.O_PATH | os.O_NOFOLLOW | os.O_CLOEXEC
        )
        components = (*workspace.parts[1:], *PurePosixPath(relative_path).parts)
        if not components:
            raise ValueError("artifact path is empty")

        def open_chain() -> tuple[list[int], tuple[os.stat_result, ...]]:
            descriptors: list[int] = []
            infos: list[os.stat_result] = []
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
                        raise ValueError("artifact parent is not a directory")
                    infos.append(info)
                return descriptors, tuple(infos)
            except Exception:
                for descriptor in reversed(descriptors):
                    try:
                        os.close(descriptor)
                    except OSError:
                        pass
                raise

        first_descriptors: list[int] = []
        reopened_descriptors: list[int] = []
        try:
            first_descriptors, first_infos = open_chain()
            file_info = first_infos[-1]
            if (
                not stat.S_ISREG(file_info.st_mode)
                or file_info.st_size < 0
                or file_info.st_size > _MAX_ARTIFACT_BYTES
            ):
                raise ValueError("artifact is not a bounded regular file")
            if content_bound and file_info.st_size > _MAX_QC_REVISION_SOURCE_BYTES:
                raise ValueError("QC source artifact exceeds its content bound")
            content: bytes | None = None
            if content_bound:
                descriptor = first_descriptors[-1]
                chunks: list[bytes] = []
                remaining = file_info.st_size + 1
                while remaining:
                    chunk = os.read(descriptor, min(65536, remaining))
                    if not chunk:
                        break
                    chunks.append(chunk)
                    remaining -= len(chunk)
                content = b"".join(chunks)
                after_read = os.fstat(descriptor)
                if (
                    len(content) != file_info.st_size
                    or ArtifactExtractionService._descriptor_identity(after_read)
                    != ArtifactExtractionService._descriptor_identity(file_info)
                ):
                    raise ValueError("artifact changed while its revision was read")
            reopened_descriptors, reopened_infos = open_chain()
            if tuple(
                ArtifactExtractionService._descriptor_identity(info)
                for info in first_infos
            ) != tuple(
                ArtifactExtractionService._descriptor_identity(info)
                for info in reopened_infos
            ):
                raise ValueError("artifact path changed while it was inspected")
            if content is not None:
                revision = build_artifact_content_revision(
                    output_type=output_type,
                    relative_path=relative_path,
                    content=content,
                )
            else:
                descriptor_identity = ":".join(
                    str(value)
                    for value in ArtifactExtractionService._descriptor_identity(
                        file_info
                    )
                )
                revision = build_artifact_descriptor_revision(
                    output_type=output_type,
                    relative_path=relative_path,
                    descriptor_identity=descriptor_identity,
                )
            return workspace / relative_path, file_info.st_size, revision
        finally:
            for descriptor in reversed((*first_descriptors, *reopened_descriptors)):
                try:
                    os.close(descriptor)
                except OSError:
                    pass

    @staticmethod
    def _descriptor_identity(info: os.stat_result) -> tuple[int, ...]:
        return (
            info.st_dev,
            info.st_ino,
            info.st_mode,
            info.st_nlink,
            info.st_uid,
            info.st_gid,
            info.st_size,
            info.st_mtime_ns,
            info.st_ctime_ns,
        )

    @staticmethod
    def _artifact_id(output_type: str, relative_path: str) -> str:
        if (
            not isinstance(output_type, str)
            or _SAFE_OUTPUT_TYPE.fullmatch(output_type) is None
        ):
            raise ValueError("artifact output type is invalid")
        digest = sha256()
        for value in (output_type.encode(), relative_path.encode()):
            digest.update(len(value).to_bytes(8, "big"))
            digest.update(value)
        return f"artifact-{digest.hexdigest()}"

    @staticmethod
    def _validated_metadata(
        metadata: Mapping[str, object],
        workspace: Path,
    ) -> dict[str, object]:
        if not isinstance(metadata, Mapping):
            raise ValueError("artifact metadata is invalid")
        result: dict[str, object] = {}
        workspace_value = str(workspace)
        for key, value in metadata.items():
            if (
                not isinstance(key, str)
                or key in _RESERVED_METADATA
                or _SAFE_METADATA_KEY.fullmatch(key) is None
            ):
                raise ValueError("artifact metadata key is invalid")
            if value is None or isinstance(value, bool) or isinstance(value, int):
                result[key] = value
            elif isinstance(value, float) and math.isfinite(value):
                result[key] = value
            elif isinstance(value, str):
                if (
                    len(value) > 512
                    or "\x00" in value
                    or value.startswith(("/", "~"))
                    or (len(value) > 1 and value[1] == ":")
                    or workspace_value in value
                ):
                    raise ValueError("artifact metadata value is invalid")
                result[key] = value
            else:
                raise ValueError("artifact metadata value is invalid")
        return result

    def _fail(
        self,
        run_id: str,
        reason_code: str,
        attempt_id: str | None,
    ):
        if attempt_id is not None:
            self._record_failure_event(run_id, reason_code, attempt_id)
        return self._result_failure(reason_code)

    def _record_failure_event(
        self,
        run_id: str,
        reason_code: str,
        attempt_id: str,
    ) -> None:
        try:
            current = self._run_service.get_run(run_id)
            if current.status is not RunStatus.SUCCEEDED:
                return
            self._run_service.record_artifact_failure(
                run_id,
                attempt_id=attempt_id,
                reason_code=reason_code,
            )
        except Exception:
            return

    @staticmethod
    def _result_failure(reason_code: str):
        return Result.failure(
            [
                Issue(
                    code="RUN_ARTIFACT_EXTRACTION_FAILED",
                    message="Workflow artifacts could not be indexed.",
                    severity="error",
                    path="artifacts",
                    source="artifact_extraction_service",
                    context={"reason_code": reason_code},
                )
            ]
        )
