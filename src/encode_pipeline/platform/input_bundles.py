"""Workflow-neutral values for inspected external input bundles."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path, PurePosixPath
import re

from encode_pipeline.platform.adapters import WorkflowInputs


_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_FINGERPRINT = re.compile(r"^sha256:[0-9a-f]{64}$")
_SAFE_RELATIVE_PATH = re.compile(
    r"^[A-Za-z0-9_][A-Za-z0-9._-]*(/[A-Za-z0-9_][A-Za-z0-9._-]*)*$"
)


@dataclass(frozen=True)
class InputBundleIdentity:
    """Public contract coordinate and content identity for one inspected Bundle."""

    contract_id: str
    contract_version: str
    schema_sha256: str
    bundle_sha256: str
    producer_name: str
    producer_version: str
    canonical_project_identity: str

    def __post_init__(self) -> None:
        for name in (
            "contract_id",
            "contract_version",
            "producer_name",
            "producer_version",
        ):
            value = getattr(self, name)
            if not isinstance(value, str) or not value or value != value.strip():
                raise ValueError(f"{name} must be a non-empty normalized string")
        for name in ("schema_sha256", "bundle_sha256"):
            if _SHA256.fullmatch(getattr(self, name)) is None:
                raise ValueError(f"{name} must be a lowercase SHA-256 digest")
        if _FINGERPRINT.fullmatch(self.canonical_project_identity) is None:
            raise ValueError("canonical_project_identity must be a sha256 fingerprint")

    def to_dict(self) -> dict[str, str]:
        """Return a JSON-safe projection without local filesystem information."""
        return {
            "contract_id": self.contract_id,
            "contract_version": self.contract_version,
            "schema_sha256": self.schema_sha256,
            "bundle_sha256": self.bundle_sha256,
            "producer_name": self.producer_name,
            "producer_version": self.producer_version,
            "canonical_project_identity": self.canonical_project_identity,
        }


@dataclass(frozen=True)
class InputBundleArtifact:
    """Digest-verified validator-facing artifact carried by a Bundle."""

    role: str
    relative_path: str
    media_type: str
    size_bytes: int
    sha256: str
    content: bytes

    def __post_init__(self) -> None:
        for name in ("role", "media_type"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value:
                raise ValueError(f"{name} must be a non-empty string")
        validate_input_bundle_relative_path(self.relative_path)
        if (
            isinstance(self.size_bytes, bool)
            or not isinstance(self.size_bytes, int)
            or self.size_bytes < 0
        ):
            raise ValueError("size_bytes must be a non-negative integer")
        if _SHA256.fullmatch(self.sha256) is None:
            raise ValueError("sha256 must be a lowercase SHA-256 digest")
        if not isinstance(self.content, bytes):
            raise ValueError("content must be bytes")
        if len(self.content) != self.size_bytes:
            raise ValueError("artifact content size does not match its contract")
        if hashlib.sha256(self.content).hexdigest() != self.sha256:
            raise ValueError("artifact content digest does not match its contract")


@dataclass(frozen=True)
class InputBundleFile:
    """One producer-verified project-local data-file binding."""

    file_id: str
    relative_path: str
    file_format: str
    role: str
    read_number: int | None
    size_bytes: int
    sha256: str

    def __post_init__(self) -> None:
        for name in ("file_id", "file_format", "role"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value:
                raise ValueError(f"{name} must be a non-empty string")
        validate_input_bundle_relative_path(self.relative_path)
        if self.read_number not in (None, 1, 2):
            raise ValueError("read_number must be 1, 2, or None")
        if (
            isinstance(self.size_bytes, bool)
            or not isinstance(self.size_bytes, int)
            or self.size_bytes < 0
        ):
            raise ValueError("size_bytes must be a non-negative integer")
        if _SHA256.fullmatch(self.sha256) is None:
            raise ValueError("sha256 must be a lowercase SHA-256 digest")


@dataclass(frozen=True)
class WorkflowInputBundle:
    """Verified public Bundle content offered to a trusted workflow adapter."""

    identity: InputBundleIdentity
    workflow_name: str
    genome: str
    render_contract: str
    project_root: Path
    artifacts: tuple[InputBundleArtifact, ...]
    files: tuple[InputBundleFile, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.identity, InputBundleIdentity):
            raise ValueError("identity must be an InputBundleIdentity")
        for name in ("workflow_name", "genome", "render_contract"):
            value = getattr(self, name)
            if not isinstance(value, str) or not value:
                raise ValueError(f"{name} must be a non-empty string")
        if (
            not isinstance(self.project_root, Path)
            or not self.project_root.is_absolute()
        ):
            raise ValueError("project_root must be an absolute Path")
        artifacts = tuple(self.artifacts)
        files = tuple(self.files)
        if not all(isinstance(item, InputBundleArtifact) for item in artifacts):
            raise ValueError("artifacts must contain InputBundleArtifact values")
        if not all(isinstance(item, InputBundleFile) for item in files):
            raise ValueError("files must contain InputBundleFile values")
        artifact_keys = [(item.role, item.relative_path) for item in artifacts]
        if artifact_keys != sorted(artifact_keys) or len(artifact_keys) != len(
            set(artifact_keys)
        ):
            raise ValueError("artifacts must be sorted and unique")
        file_keys = [(item.relative_path, item.file_id) for item in files]
        if file_keys != sorted(file_keys) or len(file_keys) != len(set(file_keys)):
            raise ValueError("files must be sorted and unique")
        object.__setattr__(self, "artifacts", artifacts)
        object.__setattr__(self, "files", files)

    def artifact(self, role: str) -> InputBundleArtifact:
        """Return the unique verified artifact for a public role."""
        matches = tuple(item for item in self.artifacts if item.role == role)
        if len(matches) != 1:
            raise KeyError(role)
        return matches[0]

    def file_for_path(self, relative_path: str) -> InputBundleFile | None:
        """Return the producer-verified local binding for a relative path."""
        validate_input_bundle_relative_path(relative_path)
        matches = tuple(
            item for item in self.files if item.relative_path == relative_path
        )
        if len(matches) > 1:
            raise ValueError("bundle file path is ambiguous")
        return matches[0] if matches else None

    def project_path(self, relative_path: str) -> Path:
        """Lexically map one public project-relative path under the inspected root."""
        validate_input_bundle_relative_path(relative_path)
        return self.project_root.joinpath(*PurePosixPath(relative_path).parts)


@dataclass(frozen=True)
class InputBundleFileSetAlternatives:
    """Ordered complete-file-set alternatives declared by one adapter.

    Each inner tuple is an atomic candidate: every path must be safely
    available. The outer tuple order is the adapter's deterministic
    preference order. This value carries paths only; it never carries a
    reader, descriptor, callback, or filesystem observation.
    """

    alternatives: tuple[tuple[str, ...], ...]

    def __post_init__(self) -> None:
        if not isinstance(self.alternatives, (tuple, list)) or any(
            not isinstance(candidate, (tuple, list)) for candidate in self.alternatives
        ):
            raise ValueError("alternatives must contain complete path tuples")
        alternatives = tuple(tuple(candidate) for candidate in self.alternatives)
        if not alternatives:
            raise ValueError("alternatives must contain at least one file set")
        for candidate in alternatives:
            if not candidate:
                raise ValueError("each alternative file set must be non-empty")
            for path in candidate:
                validate_input_bundle_relative_path(path)
            if candidate != tuple(sorted(set(candidate))):
                raise ValueError(
                    "alternative file sets must contain sorted, unique paths"
                )
        if len(alternatives) != len(set(alternatives)):
            raise ValueError("alternative file sets must be unique")
        object.__setattr__(self, "alternatives", alternatives)


@dataclass(frozen=True)
class InputBundleMapping:
    """Adapter-owned mapping plus every project file needed for validation."""

    inputs: WorkflowInputs
    required_project_files: tuple[str, ...]
    required_project_file_sets: tuple[InputBundleFileSetAlternatives, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.inputs, WorkflowInputs):
            raise ValueError("inputs must be WorkflowInputs")
        paths = tuple(self.required_project_files)
        for path in paths:
            validate_input_bundle_relative_path(path)
        file_sets = tuple(self.required_project_file_sets)
        if not all(
            isinstance(item, InputBundleFileSetAlternatives) for item in file_sets
        ):
            raise ValueError(
                "required_project_file_sets must contain "
                "InputBundleFileSetAlternatives values"
            )
        if not paths and not file_sets:
            raise ValueError("an input Bundle mapping must require project files")
        if paths != tuple(sorted(set(paths))):
            raise ValueError("required_project_files must be sorted and unique")
        if file_sets != tuple(
            sorted(file_sets, key=lambda item: item.alternatives)
        ) or len(file_sets) != len(set(file_sets)):
            raise ValueError("required_project_file_sets must be sorted and unique")
        declared_paths = list(paths)
        declared_paths.extend(
            path
            for file_set in file_sets
            for candidate in file_set.alternatives
            for path in candidate
        )
        if len(declared_paths) != len(set(declared_paths)):
            raise ValueError(
                "required paths and alternative file-set paths must be disjoint"
            )
        object.__setattr__(self, "required_project_files", paths)
        object.__setattr__(self, "required_project_file_sets", file_sets)


@dataclass(frozen=True)
class InputBundleFileObservation:
    """Consumer-side digest observation made during one read-only inspection."""

    relative_path: str
    size_bytes: int
    sha256: str
    contract_bound: bool

    def __post_init__(self) -> None:
        validate_input_bundle_relative_path(self.relative_path)
        if (
            isinstance(self.size_bytes, bool)
            or not isinstance(self.size_bytes, int)
            or self.size_bytes < 0
        ):
            raise ValueError("size_bytes must be a non-negative integer")
        if _SHA256.fullmatch(self.sha256) is None:
            raise ValueError("sha256 must be a lowercase SHA-256 digest")
        if not isinstance(self.contract_bound, bool):
            raise ValueError("contract_bound must be a boolean")


@dataclass(frozen=True)
class ImportedWorkflowInputs:
    """Ephemeral read-only import result; no snapshot or run is implied."""

    identity: InputBundleIdentity
    workflow_id: str
    inputs: WorkflowInputs
    source_files: tuple[InputBundleFileObservation, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.identity, InputBundleIdentity):
            raise ValueError("identity must be an InputBundleIdentity")
        if not isinstance(self.workflow_id, str) or not self.workflow_id:
            raise ValueError("workflow_id must be a non-empty string")
        if not isinstance(self.inputs, WorkflowInputs):
            raise ValueError("inputs must be WorkflowInputs")
        source_files = tuple(self.source_files)
        if not all(
            isinstance(item, InputBundleFileObservation) for item in source_files
        ):
            raise ValueError(
                "source_files must contain InputBundleFileObservation values"
            )
        keys = [item.relative_path for item in source_files]
        if keys != sorted(keys) or len(keys) != len(set(keys)):
            raise ValueError("source_files must be sorted and unique")
        object.__setattr__(self, "source_files", source_files)


def validate_input_bundle_relative_path(value: str) -> str:
    """Validate the closed project-relative POSIX path syntax in Bundle 0.2."""
    if (
        not isinstance(value, str)
        or not value
        or len(value) > 512
        or _SAFE_RELATIVE_PATH.fullmatch(value) is None
        or "\\" in value
        or any(ord(character) < 32 or ord(character) == 127 for character in value)
    ):
        raise ValueError("input bundle path is invalid")
    path = PurePosixPath(value)
    if (
        path.is_absolute()
        or path.as_posix() != value
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise ValueError("input bundle path is invalid")
    return value
