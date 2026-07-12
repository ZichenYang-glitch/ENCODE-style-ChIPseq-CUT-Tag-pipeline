"""Safe post-success extraction of adapter-owned artifact candidates."""

from __future__ import annotations

from collections.abc import Mapping
from hashlib import sha256
import math
import os
from pathlib import Path, PurePosixPath
import re
import stat
from urllib.parse import quote

from encode_pipeline.platform.adapters import (
    ExtractedArtifactCandidate,
    WorkflowInputs,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.runs import RunArtifactRef, RunStatus
from encode_pipeline.services.run_repositories import ConcurrentRunUpdateError
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider


_RESERVED_METADATA = frozenset({"relative_path", "output_type", "size_bytes"})
_SAFE_METADATA_KEY = re.compile(r"^[A-Za-z][A-Za-z0-9_]{0,63}$")


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
        if not isinstance(workspace_root, Path) or not workspace_root.is_absolute():
            raise ValueError("workspace_root must be an absolute Path")
        self._run_service = run_service
        self._registry = registry
        self._build_identity_provider = build_identity_provider
        self._workspace_root = workspace_root

    def extract(self, run_id: str) -> Result[tuple[RunArtifactRef, ...]]:
        """Index one canonical SUCCEEDED run without changing its lifecycle."""
        try:
            record = self._run_service.get_run(run_id)
        except KeyError:
            return self._result_failure("ARTIFACT_EXTRACTION_RUN_NOT_FOUND")
        if record.status is not RunStatus.SUCCEEDED:
            return self._result_failure("ARTIFACT_EXTRACTION_RUN_NOT_SUCCEEDED")

        try:
            adapter = self._registry.get(record.workflow_id)
            if "artifact_extract" not in adapter.capabilities.supports:
                return self._fail(run_id, "ARTIFACT_EXTRACTION_UNSUPPORTED")
            if not self._build_matches(record.run_id, record.workflow_id):
                return self._fail(run_id, "ARTIFACT_EXTRACTION_BUILD_MISMATCH")
            inputs = self._reconstruct_inputs(record.inputs)
            workspace = self._workspace_for_run(run_id)
            candidates_result = adapter.extract_artifacts(inputs, workspace)
            if candidates_result.is_failure:
                return self._fail(run_id, "ARTIFACT_EXTRACTION_ADAPTER_FAILED")
            references = self._build_references(
                record.run_id,
                record.ended_at,
                workspace,
                candidates_result.value,
            )
            self._run_service.replace_artifacts(run_id, references)
            return Result.success(references)
        except ConcurrentRunUpdateError:
            return self._fail(run_id, "ARTIFACT_EXTRACTION_STATE_CHANGED")
        except (KeyError, OSError, TypeError, ValueError):
            return self._fail(run_id, "ARTIFACT_EXTRACTION_VALIDATION_FAILED")
        except Exception:
            return self._fail(run_id, "ARTIFACT_EXTRACTION_UNEXPECTED_ERROR")

    def record_unexpected_failure(self, run_id: str) -> None:
        """Best-effort safe event for a worker-side defensive catch."""
        self._record_failure_event(run_id, "ARTIFACT_EXTRACTION_UNEXPECTED_ERROR")

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
        return self._workspace_root / run_id

    def _build_references(
        self,
        run_id: str,
        ended_at,
        workspace: Path,
        candidates,
    ) -> tuple[RunArtifactRef, ...]:
        if ended_at is None:
            raise ValueError("succeeded run has no ended_at")
        if not isinstance(candidates, tuple):
            raise ValueError("adapter candidates must be a tuple")
        workspace_stat = os.lstat(workspace)
        if stat.S_ISLNK(workspace_stat.st_mode) or not stat.S_ISDIR(
            workspace_stat.st_mode
        ):
            raise ValueError("run workspace is invalid")
        workspace_resolved = workspace.resolve(strict=True)
        references: list[RunArtifactRef] = []
        seen_paths: set[str] = set()
        seen_ids: set[str] = set()
        for candidate in candidates:
            if not isinstance(candidate, ExtractedArtifactCandidate):
                raise ValueError("adapter candidate type is invalid")
            relative_path = self._validate_relative_path(candidate.relative_path)
            if relative_path in seen_paths:
                raise ValueError("duplicate artifact path")
            seen_paths.add(relative_path)
            target, size_bytes = self._safe_regular_file(
                workspace,
                workspace_resolved,
                relative_path,
            )
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
            if mime_type is not None and (
                not mime_type or len(mime_type) > 255 or "\x00" in mime_type
            ):
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
            or any(part in {"", ".", ".."} for part in value.split("/"))
        ):
            raise ValueError("artifact path is invalid")
        return value

    @staticmethod
    def _safe_regular_file(
        workspace: Path,
        workspace_resolved: Path,
        relative_path: str,
    ) -> tuple[Path, int]:
        target = workspace
        parts = PurePosixPath(relative_path).parts
        for index, part in enumerate(parts):
            target = target / part
            file_stat = os.lstat(target)
            if stat.S_ISLNK(file_stat.st_mode):
                raise ValueError("artifact path contains a symlink")
            if index < len(parts) - 1 and not stat.S_ISDIR(file_stat.st_mode):
                raise ValueError("artifact parent is not a directory")
        if not stat.S_ISREG(file_stat.st_mode):
            raise ValueError("artifact is not a regular file")
        resolved = target.resolve(strict=True)
        try:
            resolved.relative_to(workspace_resolved)
        except ValueError as exc:
            raise ValueError("artifact escaped the run workspace") from exc
        return target, file_stat.st_size

    @staticmethod
    def _artifact_id(output_type: str, relative_path: str) -> str:
        if (
            not isinstance(output_type, str)
            or not output_type
            or len(output_type) > 128
            or "\x00" in output_type
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

    def _fail(self, run_id: str, reason_code: str):
        self._record_failure_event(run_id, reason_code)
        return self._result_failure(reason_code)

    def _record_failure_event(self, run_id: str, reason_code: str) -> None:
        try:
            current = self._run_service.get_run(run_id)
            if current.status is not RunStatus.SUCCEEDED:
                return
            self._run_service.add_event(
                run_id,
                "artifact_extraction_failed",
                "Workflow artifacts could not be indexed.",
                status=RunStatus.SUCCEEDED,
                stage="artifact_extraction",
                context={"reason_code": reason_code},
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
