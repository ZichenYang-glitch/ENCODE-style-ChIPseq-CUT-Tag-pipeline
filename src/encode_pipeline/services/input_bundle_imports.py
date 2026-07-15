"""Read-only inspection of versioned Omics Intake Bundle contracts."""

from __future__ import annotations

from contextlib import suppress
from dataclasses import dataclass
from datetime import datetime
from functools import lru_cache
import hashlib
from importlib.resources import files
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
from typing import Any, cast

from jsonschema import Draft202012Validator, FormatChecker

from encode_pipeline.platform.adapters import (
    INPUT_BUNDLE_IMPORT_CAPABILITY,
    InputBundleImportingAdapter,
    WorkflowInputs,
)
from encode_pipeline.platform.input_bundles import (
    ImportedWorkflowInputs,
    InputBundleArtifact,
    InputBundleFile,
    InputBundleFileObservation,
    InputBundleIdentity,
    InputBundleMapping,
    WorkflowInputBundle,
    validate_input_bundle_relative_path,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result


OMICS_INTAKE_BUNDLE_FILENAME = "intake-bundle.json"
OMICS_INTAKE_BUNDLE_CONTRACT_VERSION = "0.2"
OMICS_INTAKE_BUNDLE_CONTRACT_ID = "urn:omics-intake:schema:intake-bundle:0.2"
OMICS_INTAKE_BUNDLE_SCHEMA_SHA256 = (
    "6dcda336c9f0ba763383ddd58bec280946f8970af8bd730eb758fab5e3a8dd71"
)
# The release label and annotated tag object are audit provenance only. Runtime
# acceptance remains offline and pinned to the contract version and schema digest.
OMICS_INTAKE_BUNDLE_RELEASE_TAG = "v0.2.0"
OMICS_INTAKE_BUNDLE_RELEASE_TAG_OBJECT = "140a454d1313b19b322a825a1feebbb1494297c7"
OMICS_INTAKE_BUNDLE_RELEASE_COMMIT = "32680c12465f543214ed7e0173c639e0d40c7113"
OMICS_INTAKE_BUNDLE_RELEASE_TREE = "48aba2f48fa88fc37dab19c10f0ce70f2641add2"

_SCHEMA_PACKAGE = "encode_pipeline.contracts.omics_intake"
_SCHEMA_RESOURCE = "intake-bundle-0.2.schema.json"
_MAX_BUNDLE_BYTES = 16 * 1024 * 1024
_MAX_ARTIFACT_BYTES = 2 * 1024 * 1024
_MAX_FILE_REFERENCES = 2_000
_MAX_ISSUE_SUMMARIES = 10_000
_MAX_JSON_DEPTH = 64
_MAX_INTEGER = 2**63 - 1
_MAX_REQUIRED_FILE_BYTES = 256 * 1024**3
_MAX_REQUIRED_FILES_BYTES = 1024**4
_READ_CHUNK_BYTES = 1024 * 1024
_RFC3339_DATETIME = re.compile(
    r"^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}"
    r"(?:\.[0-9]+)?(?:Z|[+-][0-9]{2}:[0-9]{2})$"
)


_PUBLIC_FAILURE_MESSAGES = {
    "INPUT_BUNDLE_SOURCE_UNSAFE": (
        "The input Bundle source is missing, changed, or unsafe."
    ),
    "INPUT_BUNDLE_TOO_LARGE": "The input Bundle exceeds a supported limit.",
    "INPUT_BUNDLE_DOCUMENT_INVALID": (
        "The input Bundle does not satisfy its closed document contract."
    ),
    "INPUT_BUNDLE_CONTRACT_UNSUPPORTED": (
        "The input Bundle contract version is not supported."
    ),
    "INPUT_BUNDLE_CONTRACT_UNAVAILABLE": (
        "The pinned input Bundle contract is unavailable."
    ),
    "INPUT_BUNDLE_INTEGRITY_INVALID": (
        "The input Bundle integrity evidence does not match its current files."
    ),
    "INPUT_BUNDLE_HANDOFF_NOT_READY": (
        "The input Bundle is not ready for fail-closed local consumption."
    ),
    "INPUT_BUNDLE_WORKFLOW_NOT_FOUND": "The target workflow was not found.",
    "INPUT_BUNDLE_CAPABILITY_UNSUPPORTED": (
        "The target workflow does not support input Bundle import."
    ),
    "INPUT_BUNDLE_ADAPTER_FAILED": (
        "The target workflow could not map the input Bundle safely."
    ),
}


class _InputBundleFailure(RuntimeError):
    def __init__(self, code: str) -> None:
        super().__init__(code)
        self.code = code


class _InputBundlePathMissing(_InputBundleFailure):
    """A candidate path was absent before any leaf object was observed."""

    def __init__(self) -> None:
        super().__init__("INPUT_BUNDLE_SOURCE_UNSAFE")


@dataclass(frozen=True)
class _BundleContractView:
    payload: dict[str, Any]
    canonical_record: dict[str, Any]
    artifact_records: tuple[dict[str, Any], ...]
    local_file_records: tuple[dict[str, Any], ...]


@dataclass(frozen=True)
class _ReaderFileObservation:
    size_bytes: int
    sha256: str
    path_identity: tuple[tuple[int, ...], ...]


@dataclass(frozen=True)
class _ObservedMappingFile:
    relative_path: str
    size_bytes: int
    sha256: str
    contract_bound: bool
    path_identity: tuple[tuple[int, ...], ...]

    def public(self) -> InputBundleFileObservation:
        return InputBundleFileObservation(
            relative_path=self.relative_path,
            size_bytes=self.size_bytes,
            sha256=self.sha256,
            contract_bound=self.contract_bound,
        )


@dataclass(frozen=True)
class _MissingMappingFile:
    relative_path: str


@dataclass(frozen=True)
class _MappingFilesObservation:
    source_files: tuple[InputBundleFileObservation, ...]
    private_path_evidence: tuple[_ObservedMappingFile | _MissingMappingFile, ...]
    selected_file_sets: tuple[tuple[str, ...], ...]


class InputBundleImportService:
    """Inspect and map one local Bundle without writing source or platform state."""

    def __init__(self, registry: WorkflowRegistry) -> None:
        if not isinstance(registry, WorkflowRegistry):
            raise ValueError("registry must be a WorkflowRegistry")
        self._registry = registry

    def inspect(
        self,
        bundle_path: Path,
        workflow_id: str,
    ) -> Result[ImportedWorkflowInputs]:
        """Return ephemeral validated inputs and source observations only."""
        try:
            adapter = self._registry.get(workflow_id)
        except (KeyError, ValueError):
            return _failure("INPUT_BUNDLE_WORKFLOW_NOT_FOUND")
        if (
            INPUT_BUNDLE_IMPORT_CAPABILITY not in adapter.capabilities.supports
            or not isinstance(adapter, InputBundleImportingAdapter)
        ):
            return _failure("INPUT_BUNDLE_CAPABILITY_UNSUPPORTED")

        try:
            reader = _BundleProjectReader.open(bundle_path)
        except _InputBundleFailure as error:
            return _failure(error.code)

        with reader:
            try:
                bundle_bytes = reader.read_bytes(
                    OMICS_INTAKE_BUNDLE_FILENAME,
                    maximum_bytes=_MAX_BUNDLE_BYTES,
                )
                bundle_sha256 = hashlib.sha256(bundle_bytes).hexdigest()
                contract = _load_contract_view(bundle_bytes)
                bundle = _load_verified_bundle(
                    reader,
                    contract,
                    bundle_sha256=bundle_sha256,
                )
            except _InputBundleFailure as error:
                return _failure(error.code)
            except Exception:
                return _failure("INPUT_BUNDLE_DOCUMENT_INVALID")

            importer = cast(InputBundleImportingAdapter, adapter)
            try:
                mapped = importer.import_input_bundle(bundle)
            except Exception:
                return _failure("INPUT_BUNDLE_ADAPTER_FAILED")
            if not isinstance(mapped, Result):
                return _failure("INPUT_BUNDLE_ADAPTER_FAILED")
            if mapped.is_failure:
                return Result.failure(mapped.issues)
            if not isinstance(mapped.value, InputBundleMapping):
                return _failure("INPUT_BUNDLE_ADAPTER_FAILED")

            try:
                observations = _observe_mapping_files(
                    reader,
                    bundle,
                    mapped.value,
                )
            except _InputBundleFailure as error:
                return _failure(error.code)
            except Exception:
                return _failure("INPUT_BUNDLE_SOURCE_UNSAFE")

            try:
                validation = adapter.validate(mapped.value.inputs)
            except Exception:
                return _failure("INPUT_BUNDLE_ADAPTER_FAILED")
            if not isinstance(validation, Result):
                return _failure("INPUT_BUNDLE_ADAPTER_FAILED")
            if validation.is_failure:
                return Result.failure(validation.issues)

            try:
                confirmed_observations = _observe_mapping_files(
                    reader,
                    bundle,
                    mapped.value,
                )
                if confirmed_observations != observations:
                    raise _InputBundleFailure("INPUT_BUNDLE_SOURCE_UNSAFE")
                reader.verify_root_binding()
            except _InputBundleFailure as error:
                return _failure(error.code)

            return Result.success(
                ImportedWorkflowInputs(
                    identity=bundle.identity,
                    workflow_id=adapter.metadata.workflow_id,
                    inputs=WorkflowInputs(**mapped.value.inputs.to_dict()),
                    source_files=observations.source_files,
                ),
                issues=(*mapped.issues, *validation.issues),
            )


def _load_contract_view(bundle_bytes: bytes) -> _BundleContractView:
    payload = _strict_json_loads(bundle_bytes)
    if not isinstance(payload, dict):
        raise _InputBundleFailure("INPUT_BUNDLE_DOCUMENT_INVALID")
    version = payload.get("bundle_schema_version")
    if version != OMICS_INTAKE_BUNDLE_CONTRACT_VERSION:
        raise _InputBundleFailure("INPUT_BUNDLE_CONTRACT_UNSUPPORTED")

    try:
        validator = _bundle_validator()
    except Exception as error:
        raise _InputBundleFailure("INPUT_BUNDLE_CONTRACT_UNAVAILABLE") from error
    if next(validator.iter_errors(payload), None) is not None:
        raise _InputBundleFailure("INPUT_BUNDLE_DOCUMENT_INVALID")
    return _validate_semantics(payload)


@lru_cache(maxsize=1)
def _bundle_validator() -> Draft202012Validator:
    schema_bytes = files(_SCHEMA_PACKAGE).joinpath(_SCHEMA_RESOURCE).read_bytes()
    if hashlib.sha256(schema_bytes).hexdigest() != OMICS_INTAKE_BUNDLE_SCHEMA_SHA256:
        raise ValueError("pinned Bundle schema digest mismatch")
    schema = json.loads(schema_bytes)
    if (
        not isinstance(schema, dict)
        or schema.get("$id") != OMICS_INTAKE_BUNDLE_CONTRACT_ID
        or schema.get("$schema") != "https://json-schema.org/draft/2020-12/schema"
    ):
        raise ValueError("pinned Bundle schema coordinate mismatch")
    Draft202012Validator.check_schema(schema)
    format_checker = FormatChecker()
    format_checker.checks("date-time")(_is_rfc3339_datetime)
    return Draft202012Validator(schema, format_checker=format_checker)


def _is_rfc3339_datetime(value: object) -> bool:
    if not isinstance(value, str) or _RFC3339_DATETIME.fullmatch(value) is None:
        return False
    normalized = value[:-1] + "+00:00" if value.endswith("Z") else value
    try:
        parsed = datetime.fromisoformat(normalized)
    except ValueError:
        return False
    return parsed.tzinfo is not None and parsed.utcoffset() is not None


def _strict_json_loads(document: bytes) -> object:
    def object_pairs(pairs: list[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("duplicate JSON key")
            result[key] = value
        return result

    def parse_integer(value: str) -> int:
        parsed = int(value)
        if abs(parsed) > _MAX_INTEGER:
            raise ValueError("JSON integer exceeds the supported range")
        return parsed

    def reject_number(_value: str) -> float:
        raise ValueError("non-integer JSON numbers are not supported")

    try:
        value = json.loads(
            document.decode("utf-8"),
            object_pairs_hook=object_pairs,
            parse_int=parse_integer,
            parse_float=reject_number,
            parse_constant=reject_number,
        )
        if _json_depth(value) > _MAX_JSON_DEPTH:
            raise ValueError("JSON nesting exceeds the supported depth")
        return value
    except (UnicodeError, ValueError, json.JSONDecodeError, RecursionError) as error:
        raise _InputBundleFailure("INPUT_BUNDLE_DOCUMENT_INVALID") from error


def _json_depth(value: object) -> int:
    if isinstance(value, dict):
        return 1 + max((_json_depth(item) for item in value.values()), default=0)
    if isinstance(value, list):
        return 1 + max((_json_depth(item) for item in value), default=0)
    return 0


def _validate_semantics(payload: dict[str, Any]) -> _BundleContractView:
    canonical = payload["canonical_project"]
    canonical_record = canonical["record"]
    if canonical["identity"] != f"sha256:{canonical_record['sha256']}":
        raise _InputBundleFailure("INPUT_BUNDLE_INTEGRITY_INVALID")

    artifacts = tuple(payload["artifacts"])
    artifact_keys = [
        (item["role"], item["path"], item["media_type"]) for item in artifacts
    ]
    if artifact_keys != [
        ("sample_sheet", "samples.encode.tsv", "text/tab-separated-values"),
        ("workflow_config", "config.encode.yaml", "text/yaml"),
    ]:
        raise _InputBundleFailure("INPUT_BUNDLE_DOCUMENT_INVALID")

    file_records = tuple(payload["files"])
    if len(file_records) > _MAX_FILE_REFERENCES:
        raise _InputBundleFailure("INPUT_BUNDLE_TOO_LARGE")
    file_keys = [
        (item["file_id"], item["kind"], item["path"] or "", item["scope"])
        for item in file_records
    ]
    file_kinds = [(item["file_id"], item["kind"]) for item in file_records]
    planned_paths = [item["path"] for item in file_records if item["kind"] == "planned"]
    local_paths = [item["path"] for item in file_records if item["kind"] == "local"]
    if (
        file_keys != sorted(file_keys)
        or len(file_keys) != len(set(file_keys))
        or len(file_kinds) != len(set(file_kinds))
        or len(planned_paths) != len(set(planned_paths))
        or len(local_paths) != len(set(local_paths))
    ):
        raise _InputBundleFailure("INPUT_BUNDLE_DOCUMENT_INVALID")

    readiness = payload["readiness"]
    if (
        readiness["export_status"] != "ready"
        or readiness["execution_status"] != "runnable"
        or readiness["validator_status"] != "passed"
        or readiness["strict_validator_status"] != "passed"
        or readiness["unresolved_requirements"] != []
    ):
        raise _InputBundleFailure("INPUT_BUNDLE_HANDOFF_NOT_READY")

    local_by_id: dict[str, dict[str, Any]] = {}
    planned_by_id: dict[str, dict[str, Any]] = {}
    for record in file_records:
        if record["scope"] != "project" or record["kind"] == "unbound":
            raise _InputBundleFailure("INPUT_BUNDLE_HANDOFF_NOT_READY")
        if record["kind"] == "planned":
            planned_by_id[record["file_id"]] = record
            continue
        checksum = record["checksum"]
        if (
            record["kind"] != "local"
            or not isinstance(record["path"], str)
            or not isinstance(record["size_bytes"], int)
            or checksum
            != {
                "status": "verified",
                "algorithm": "sha256",
                "digest": checksum.get("digest"),
            }
        ):
            raise _InputBundleFailure("INPUT_BUNDLE_HANDOFF_NOT_READY")
        local_by_id[record["file_id"]] = record

    if not local_by_id or set(local_by_id) != {
        item["file_id"] for item in file_records
    }:
        raise _InputBundleFailure("INPUT_BUNDLE_HANDOFF_NOT_READY")
    for file_id, planned in planned_by_id.items():
        if planned["path"] != local_by_id[file_id]["path"]:
            raise _InputBundleFailure("INPUT_BUNDLE_INTEGRITY_INVALID")

    validations = tuple(payload["validations"])
    if [item["mode"] for item in validations] != ["non_strict", "strict"]:
        raise _InputBundleFailure("INPUT_BUNDLE_DOCUMENT_INVALID")
    if any(
        item["freshness"] != "passed"
        or not item["recorded"]
        or not item["invoked"]
        or item["outcome"] != "passed"
        or item["reason_code"] is not None
        or item["return_code"] not in (None, 0)
        or item["completed_at"] is None
        or item["identity_state"] != "recorded"
        or item["identity_fingerprint"] is None
        or item["artifact_fingerprint"] is None
        or (item["mode"] == "strict" and item["strict_inputs_fingerprint"] is None)
        for item in validations
    ):
        raise _InputBundleFailure("INPUT_BUNDLE_HANDOFF_NOT_READY")

    issues = tuple(payload["issues"])
    if len(issues) > _MAX_ISSUE_SUMMARIES:
        raise _InputBundleFailure("INPUT_BUNDLE_TOO_LARGE")
    issue_keys = [
        (item["severity"], item["code"], item["owner_kind"] or "") for item in issues
    ]
    if issue_keys != sorted(issue_keys) or len(issue_keys) != len(set(issue_keys)):
        raise _InputBundleFailure("INPUT_BUNDLE_DOCUMENT_INVALID")
    if any(
        item["severity"] in {"error", "needs_review"}
        or item["owner_kind"] == "transfer"
        for item in issues
    ):
        raise _InputBundleFailure("INPUT_BUNDLE_HANDOFF_NOT_READY")

    return _BundleContractView(
        payload=payload,
        canonical_record=canonical_record,
        artifact_records=artifacts,
        local_file_records=tuple(
            sorted(
                local_by_id.values(), key=lambda item: (item["path"], item["file_id"])
            )
        ),
    )


def _load_verified_bundle(
    reader: "_BundleProjectReader",
    contract: _BundleContractView,
    *,
    bundle_sha256: str,
) -> WorkflowInputBundle:
    canonical = contract.canonical_record
    reader.read_bytes(
        canonical["path"],
        maximum_bytes=_MAX_BUNDLE_BYTES,
        expected_size=canonical["size_bytes"],
        expected_sha256=canonical["sha256"],
    )

    artifacts = []
    for record in contract.artifact_records:
        content = reader.read_bytes(
            record["path"],
            maximum_bytes=_MAX_ARTIFACT_BYTES,
            expected_size=record["size_bytes"],
            expected_sha256=record["sha256"],
        )
        artifacts.append(
            InputBundleArtifact(
                role=record["role"],
                relative_path=record["path"],
                media_type=record["media_type"],
                size_bytes=record["size_bytes"],
                sha256=record["sha256"],
                content=content,
            )
        )

    local_files = tuple(
        InputBundleFile(
            file_id=record["file_id"],
            relative_path=record["path"],
            file_format=record["file_format"],
            role=record["role"],
            read_number=record["read_number"],
            size_bytes=record["size_bytes"],
            sha256=record["checksum"]["digest"],
        )
        for record in contract.local_file_records
    )
    payload = contract.payload
    workflow = payload["workflow"]
    if not isinstance(workflow, dict):
        raise _InputBundleFailure("INPUT_BUNDLE_HANDOFF_NOT_READY")
    identity = InputBundleIdentity(
        contract_id=OMICS_INTAKE_BUNDLE_CONTRACT_ID,
        contract_version=OMICS_INTAKE_BUNDLE_CONTRACT_VERSION,
        schema_sha256=OMICS_INTAKE_BUNDLE_SCHEMA_SHA256,
        bundle_sha256=bundle_sha256,
        producer_name=payload["producer"]["name"],
        producer_version=payload["producer"]["version"],
        canonical_project_identity=payload["canonical_project"]["identity"],
    )
    return WorkflowInputBundle(
        identity=identity,
        workflow_name=workflow["name"],
        genome=workflow["genome"],
        render_contract=workflow["render_contract"],
        project_root=reader.project_root,
        artifacts=tuple(artifacts),
        files=local_files,
    )


def _observe_mapping_files(
    reader: "_BundleProjectReader",
    bundle: WorkflowInputBundle,
    mapping: InputBundleMapping,
) -> _MappingFilesObservation:
    private_evidence: list[_ObservedMappingFile | _MissingMappingFile] = []
    selected_files: list[_ObservedMappingFile] = []
    selected_file_sets: list[tuple[str, ...]] = []
    observed_bytes = 0

    def observe(relative_path: str) -> _ObservedMappingFile:
        nonlocal observed_bytes
        binding = bundle.file_for_path(relative_path)
        remaining_bytes = _MAX_REQUIRED_FILES_BYTES - observed_bytes
        observation = reader.observe_file(
            relative_path,
            expected_size=None if binding is None else binding.size_bytes,
            expected_sha256=None if binding is None else binding.sha256,
            maximum_bytes=min(_MAX_REQUIRED_FILE_BYTES, remaining_bytes),
        )
        observed_bytes += observation.size_bytes
        return _ObservedMappingFile(
            relative_path=relative_path,
            size_bytes=observation.size_bytes,
            sha256=observation.sha256,
            contract_bound=binding is not None,
            path_identity=observation.path_identity,
        )

    for relative_path in mapping.required_project_files:
        observed = observe(relative_path)
        private_evidence.append(observed)
        selected_files.append(observed)

    for file_set in mapping.required_project_file_sets:
        selected: tuple[_ObservedMappingFile, ...] | None = None
        for candidate in file_set.alternatives:
            candidate_evidence: list[_ObservedMappingFile | _MissingMappingFile] = []
            candidate_observations: list[_ObservedMappingFile] = []
            complete = True
            for relative_path in candidate:
                try:
                    observed = observe(relative_path)
                except _InputBundlePathMissing:
                    complete = False
                    candidate_evidence.append(
                        _MissingMappingFile(relative_path=relative_path)
                    )
                    continue
                candidate_evidence.append(observed)
                candidate_observations.append(observed)
            private_evidence.extend(candidate_evidence)
            if complete:
                selected = tuple(candidate_observations)
                selected_file_sets.append(candidate)
                break
        if selected is None:
            raise _InputBundleFailure("INPUT_BUNDLE_SOURCE_UNSAFE")
        selected_files.extend(selected)

    ordered_private = tuple(
        sorted(private_evidence, key=lambda item: item.relative_path)
    )
    ordered_selected = tuple(
        sorted(selected_files, key=lambda item: item.relative_path)
    )
    return _MappingFilesObservation(
        source_files=tuple(item.public() for item in ordered_selected),
        private_path_evidence=ordered_private,
        selected_file_sets=tuple(selected_file_sets),
    )


class _BundleProjectReader:
    def __init__(self, project_root: Path, root_descriptor: int) -> None:
        self.project_root = project_root
        self._root_descriptor = root_descriptor
        self._root_identity = _directory_identity(os.fstat(root_descriptor))
        self._closed = False

    @classmethod
    def open(cls, bundle_path: Path) -> "_BundleProjectReader":
        if not isinstance(bundle_path, Path):
            raise _InputBundleFailure("INPUT_BUNDLE_SOURCE_UNSAFE")
        absolute = (
            bundle_path if bundle_path.is_absolute() else Path.cwd() / bundle_path
        )
        if absolute.name != OMICS_INTAKE_BUNDLE_FILENAME:
            raise _InputBundleFailure("INPUT_BUNDLE_SOURCE_UNSAFE")
        project_root = absolute.parent
        try:
            descriptor = _open_absolute_directory(project_root)
        except (OSError, ValueError) as error:
            raise _InputBundleFailure("INPUT_BUNDLE_SOURCE_UNSAFE") from error
        return cls(project_root, descriptor)

    def __enter__(self) -> "_BundleProjectReader":
        return self

    def __exit__(self, *_args: object) -> None:
        self.close()

    def close(self) -> None:
        if self._closed:
            return
        self._closed = True
        with suppress(OSError):
            os.close(self._root_descriptor)

    def read_bytes(
        self,
        relative_path: str,
        *,
        maximum_bytes: int,
        expected_size: int | None = None,
        expected_sha256: str | None = None,
    ) -> bytes:
        observation, content = self._read_file(
            relative_path,
            maximum_bytes=maximum_bytes,
            expected_size=expected_size,
            expected_sha256=expected_sha256,
            capture=True,
        )
        assert len(content) == observation.size_bytes
        assert hashlib.sha256(content).hexdigest() == observation.sha256
        return content

    def observe_file(
        self,
        relative_path: str,
        *,
        expected_size: int | None,
        expected_sha256: str | None,
        maximum_bytes: int,
    ) -> _ReaderFileObservation:
        observation, content = self._read_file(
            relative_path,
            maximum_bytes=maximum_bytes,
            expected_size=expected_size,
            expected_sha256=expected_sha256,
            capture=False,
        )
        assert content == b""
        return observation

    def verify_root_binding(self) -> None:
        try:
            descriptor = _open_absolute_directory(self.project_root)
        except (OSError, ValueError) as error:
            raise _InputBundleFailure("INPUT_BUNDLE_SOURCE_UNSAFE") from error
        try:
            if _directory_identity(os.fstat(descriptor)) != self._root_identity:
                raise _InputBundleFailure("INPUT_BUNDLE_SOURCE_UNSAFE")
        finally:
            os.close(descriptor)

    def _read_file(
        self,
        relative_path: str,
        *,
        maximum_bytes: int,
        expected_size: int | None,
        expected_sha256: str | None,
        capture: bool,
    ) -> tuple[_ReaderFileObservation, bytes]:
        try:
            validate_input_bundle_relative_path(relative_path)
        except ValueError as error:
            raise _InputBundleFailure("INPUT_BUNDLE_DOCUMENT_INVALID") from error
        descriptors: list[int] = []
        leaf_was_listed = False
        try:
            current = self._root_descriptor
            path_identity: list[tuple[int, ...]] = [self._root_identity]
            parts = PurePosixPath(relative_path).parts
            for component in parts[:-1]:
                descriptor = os.open(component, _directory_flags(), dir_fd=current)
                descriptors.append(descriptor)
                current = descriptor
                path_identity.append(_directory_identity(os.fstat(descriptor)))
            listed = os.stat(parts[-1], dir_fd=current, follow_symlinks=False)
            leaf_was_listed = True
            descriptor = os.open(parts[-1], _file_flags(), dir_fd=current)
            descriptors.append(descriptor)
            before = os.fstat(descriptor)
            if (
                not stat.S_ISREG(listed.st_mode)
                or not stat.S_ISREG(before.st_mode)
                or listed.st_nlink != 1
                or before.st_nlink != 1
                or _file_identity(listed) != _file_identity(before)
            ):
                raise _InputBundleFailure("INPUT_BUNDLE_SOURCE_UNSAFE")
            size_bytes = before.st_size
            if size_bytes > maximum_bytes:
                raise _InputBundleFailure("INPUT_BUNDLE_TOO_LARGE")
            if expected_size is not None and size_bytes != expected_size:
                raise _InputBundleFailure("INPUT_BUNDLE_INTEGRITY_INVALID")

            digest = hashlib.sha256()
            chunks: list[bytes] = []
            remaining = size_bytes
            while remaining:
                chunk = os.read(descriptor, min(_READ_CHUNK_BYTES, remaining))
                if not chunk:
                    raise _InputBundleFailure("INPUT_BUNDLE_SOURCE_UNSAFE")
                digest.update(chunk)
                if capture:
                    chunks.append(chunk)
                remaining -= len(chunk)
            if os.read(descriptor, 1):
                raise _InputBundleFailure("INPUT_BUNDLE_SOURCE_UNSAFE")
            after = os.fstat(descriptor)
            observed = digest.hexdigest()
            if _file_identity(before) != _file_identity(after):
                raise _InputBundleFailure("INPUT_BUNDLE_SOURCE_UNSAFE")
            if expected_sha256 is not None and observed != expected_sha256:
                raise _InputBundleFailure("INPUT_BUNDLE_INTEGRITY_INVALID")
            observed_path_identity = tuple((*path_identity, _file_identity(after)))
            reopened = self._reopen_identity(relative_path)
            if reopened != observed_path_identity:
                raise _InputBundleFailure("INPUT_BUNDLE_SOURCE_UNSAFE")
            return (
                _ReaderFileObservation(
                    size_bytes=size_bytes,
                    sha256=observed,
                    path_identity=observed_path_identity,
                ),
                b"".join(chunks),
            )
        except _InputBundleFailure:
            raise
        except FileNotFoundError as error:
            if not leaf_was_listed:
                raise _InputBundlePathMissing() from error
            raise _InputBundleFailure("INPUT_BUNDLE_SOURCE_UNSAFE") from error
        except (OSError, TypeError, ValueError) as error:
            raise _InputBundleFailure("INPUT_BUNDLE_SOURCE_UNSAFE") from error
        finally:
            for descriptor in reversed(descriptors):
                with suppress(OSError):
                    os.close(descriptor)

    def _reopen_identity(self, relative_path: str) -> tuple[tuple[int, ...], ...]:
        descriptors: list[int] = []
        try:
            current = self._root_descriptor
            path_identity: list[tuple[int, ...]] = [self._root_identity]
            parts = PurePosixPath(relative_path).parts
            for component in parts[:-1]:
                descriptor = os.open(component, _directory_flags(), dir_fd=current)
                descriptors.append(descriptor)
                current = descriptor
                path_identity.append(_directory_identity(os.fstat(descriptor)))
            descriptor = os.open(parts[-1], _file_flags(), dir_fd=current)
            descriptors.append(descriptor)
            info = os.fstat(descriptor)
            if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
                raise OSError("unsafe source type")
            path_identity.append(_file_identity(info))
            return tuple(path_identity)
        finally:
            for descriptor in reversed(descriptors):
                with suppress(OSError):
                    os.close(descriptor)


def _open_absolute_directory(path: Path) -> int:
    if not path.is_absolute() or any(
        part in {"", ".", ".."} for part in path.parts[1:]
    ):
        raise ValueError("directory path must be absolute and normalized")
    current = os.open("/", _directory_flags())
    try:
        for component in path.parts[1:]:
            descriptor = os.open(component, _directory_flags(), dir_fd=current)
            os.close(current)
            current = descriptor
        info = os.fstat(current)
        if not stat.S_ISDIR(info.st_mode):
            raise OSError("source root is not a directory")
        return current
    except Exception:
        with suppress(OSError):
            os.close(current)
        raise


def _directory_flags() -> int:
    required = ("O_DIRECTORY", "O_NOFOLLOW")
    if any(not hasattr(os, name) for name in required):
        raise OSError("safe directory flags are unavailable")
    return os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | getattr(os, "O_CLOEXEC", 0)


def _file_flags() -> int:
    if not hasattr(os, "O_NOFOLLOW"):
        raise OSError("safe file flags are unavailable")
    return (
        os.O_RDONLY
        | os.O_NOFOLLOW
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NONBLOCK", 0)
        | getattr(os, "O_NOCTTY", 0)
    )


def _file_identity(info: os.stat_result) -> tuple[int, ...]:
    return (
        info.st_dev,
        info.st_ino,
        info.st_mode,
        info.st_nlink,
        info.st_size,
        info.st_mtime_ns,
        info.st_ctime_ns,
    )


def _directory_identity(info: os.stat_result) -> tuple[int, int]:
    return (info.st_dev, info.st_ino)


def _failure(code: str) -> Result[ImportedWorkflowInputs]:
    message = _PUBLIC_FAILURE_MESSAGES.get(
        code,
        _PUBLIC_FAILURE_MESSAGES["INPUT_BUNDLE_DOCUMENT_INVALID"],
    )
    return Result.failure(
        [
            Issue(
                code=code,
                message=message,
                source="input_bundle_import",
                path="bundle",
                context={},
            )
        ]
    )
