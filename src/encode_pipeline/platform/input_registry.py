"""Workflow-neutral local input registry and immutable binding primitives."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from enum import Enum
from hashlib import sha256
import json
from pathlib import PurePosixPath
import re


INPUT_FILE_REVISION_DIGEST_SCHEME = "sha256-framed-input-file-revision-v1"
INPUT_CLOSURE_DIGEST_SCHEME = "sha256-framed-input-closure-v1"
INPUT_USE_BINDING_DIGEST_SCHEME = "sha256-framed-input-use-binding-envelope-v1"

_PROJECT_ID = re.compile(r"^prj_[0-9a-f]{32}$")
_STORAGE_POOL_ID = re.compile(r"^stgp_[0-9a-f]{32}$")
_INPUT_FILE_ID = re.compile(r"^inpf_[0-9a-f]{32}$")
_INPUT_FILE_REVISION_ID = re.compile(r"^inpfr_[0-9a-f]{32}$")
_SAFE_KEY = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]{0,254}$")
_CONTRACT_TOKEN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:+/-]{0,254}$")
_DIGEST = re.compile(r"^[0-9a-f]{64}$")
_MAX_RELATIVE_PATH_BYTES = 4096
_MAX_COMPONENT_BYTES = 255
_MAX_INPUT_COLLECTION_ITEMS = 256


class InputProvenanceMode(str, Enum):
    """Per-input-use provenance classification frozen into a snapshot."""

    MANAGED_REVISION_V1 = "managed_revision_v1"
    TRANSITIONAL_UNMANAGED_V1 = "transitional_unmanaged_v1"


class InputBindingContractMode(str, Enum):
    """Top-level compatibility mode for input provenance evidence."""

    COMPATIBILITY_UNRESOLVED_V1 = "compatibility_unresolved_v1"
    DECLARED_INPUT_USES_V1 = "declared_input_uses_v1"


def _utc(value: datetime, name: str) -> datetime:
    if (
        not isinstance(value, datetime)
        or value.tzinfo is None
        or value.utcoffset() is None
    ):
        raise ValueError(f"{name} must be timezone-aware")
    return value.astimezone(timezone.utc)


def _validate_id(value: str, pattern: re.Pattern[str], name: str) -> str:
    if not isinstance(value, str) or pattern.fullmatch(value) is None:
        raise ValueError(f"{name} must be an opaque {name.replace('_', ' ')}")
    return value


def _validate_digest(value: str, name: str) -> str:
    if not isinstance(value, str) or _DIGEST.fullmatch(value) is None:
        raise ValueError(f"{name} must be a lowercase SHA-256 digest")
    return value


def _safe_key(value: str, name: str) -> str:
    if not isinstance(value, str) or _SAFE_KEY.fullmatch(value) is None:
        raise ValueError(f"{name} must be a safe opaque key")
    return value


def _workflow_id(value: str) -> str:
    if not isinstance(value, str):
        raise ValueError("workflow_id must be a string")
    normalized = value.strip()
    if not normalized:
        raise ValueError("workflow_id must be non-empty")
    return normalized


def validate_input_file_stable_key(stable_key: str) -> str:
    """Validate one stable logical InputFile key without observing file bytes."""
    return _safe_key(stable_key, "stable_key")


def _contract_token(value: str, name: str) -> str:
    if not isinstance(value, str) or _CONTRACT_TOKEN.fullmatch(value) is None:
        raise ValueError(f"{name} must be a safe versioned contract token")
    return value


def _display_name(value: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError("display_name must be a non-empty string")
    normalized = value.strip()
    if len(normalized) > 255 or any(
        ord(character) < 32 or ord(character) == 127 for character in normalized
    ):
        raise ValueError("display_name is invalid")
    return normalized


def _positive_integer(value: int, name: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _non_negative_integer(value: int, name: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValueError(f"{name} must be a non-negative integer")
    return value


def _framed_sha256(*parts: str) -> str:
    digest = sha256()
    for part in parts:
        encoded = part.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return digest.hexdigest()


def validate_input_file_relative_path(relative_path: str) -> str:
    """Validate a nested POSIX path without resolving it against a filesystem."""
    if (
        not isinstance(relative_path, str)
        or not relative_path
        or relative_path.startswith("/")
        or "\\" in relative_path
        or any(
            ord(character) < 32 or ord(character) == 127 for character in relative_path
        )
        or len(relative_path.encode("utf-8")) > _MAX_RELATIVE_PATH_BYTES
    ):
        raise ValueError("relative_path must be a safe relative POSIX path")
    components = relative_path.split("/")
    if (
        any(component in {"", ".", ".."} for component in components)
        or any(
            len(component.encode("utf-8")) > _MAX_COMPONENT_BYTES
            for component in components
        )
        or PurePosixPath(relative_path).is_absolute()
    ):
        raise ValueError("relative_path must be a safe relative POSIX path")
    return relative_path


@dataclass(frozen=True)
class StoragePool:
    """Opaque database identity for one operator-configured storage root."""

    storage_pool_id: str
    config_key: str
    display_name: str
    created_at: datetime
    archived_at: datetime | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "storage_pool_id",
            _validate_id(self.storage_pool_id, _STORAGE_POOL_ID, "storage_pool_id"),
        )
        object.__setattr__(self, "config_key", _safe_key(self.config_key, "config_key"))
        object.__setattr__(self, "display_name", _display_name(self.display_name))
        created_at = _utc(self.created_at, "created_at")
        archived_at = (
            None if self.archived_at is None else _utc(self.archived_at, "archived_at")
        )
        if archived_at is not None and archived_at < created_at:
            raise ValueError("archived_at cannot precede created_at")
        object.__setattr__(self, "created_at", created_at)
        object.__setattr__(self, "archived_at", archived_at)


@dataclass(frozen=True)
class ProjectStoragePoolBinding:
    """Current approved pool binding for one Project."""

    project_id: str
    storage_pool_id: str
    bound_at: datetime

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "project_id",
            _validate_id(self.project_id, _PROJECT_ID, "project_id"),
        )
        object.__setattr__(
            self,
            "storage_pool_id",
            _validate_id(self.storage_pool_id, _STORAGE_POOL_ID, "storage_pool_id"),
        )
        object.__setattr__(self, "bound_at", _utc(self.bound_at, "bound_at"))


@dataclass(frozen=True)
class InputFile:
    """Stable logical identity for one Project-owned regular input file."""

    input_file_id: str
    project_id: str
    storage_pool_id: str
    stable_key: str
    created_at: datetime
    archived_at: datetime | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "input_file_id",
            _validate_id(self.input_file_id, _INPUT_FILE_ID, "input_file_id"),
        )
        object.__setattr__(
            self,
            "project_id",
            _validate_id(self.project_id, _PROJECT_ID, "project_id"),
        )
        object.__setattr__(
            self,
            "storage_pool_id",
            _validate_id(self.storage_pool_id, _STORAGE_POOL_ID, "storage_pool_id"),
        )
        object.__setattr__(
            self,
            "stable_key",
            validate_input_file_stable_key(self.stable_key),
        )
        created_at = _utc(self.created_at, "created_at")
        archived_at = (
            None if self.archived_at is None else _utc(self.archived_at, "archived_at")
        )
        if archived_at is not None and archived_at < created_at:
            raise ValueError("archived_at cannot precede created_at")
        object.__setattr__(self, "created_at", created_at)
        object.__setattr__(self, "archived_at", archived_at)


@dataclass(frozen=True)
class InputFileRevision:
    """One append-only observed regular-file revision."""

    input_file_revision_id: str
    input_file_id: str
    project_id: str
    storage_pool_id: str
    revision_number: int
    relative_path: str
    size_bytes: int
    content_sha256: str
    digest_scheme: str
    digest: str
    created_at: datetime

    def __post_init__(self) -> None:
        _validate_id(
            self.input_file_revision_id,
            _INPUT_FILE_REVISION_ID,
            "input_file_revision_id",
        )
        _validate_id(self.input_file_id, _INPUT_FILE_ID, "input_file_id")
        _validate_id(self.project_id, _PROJECT_ID, "project_id")
        _validate_id(self.storage_pool_id, _STORAGE_POOL_ID, "storage_pool_id")
        _positive_integer(self.revision_number, "revision_number")
        validate_input_file_relative_path(self.relative_path)
        _non_negative_integer(self.size_bytes, "size_bytes")
        _validate_digest(self.content_sha256, "content_sha256")
        if self.digest_scheme != INPUT_FILE_REVISION_DIGEST_SCHEME:
            raise ValueError("input file revision digest scheme is unsupported")
        digest = _validate_digest(self.digest, "digest")
        expected = _input_file_revision_digest(
            input_file_revision_id=self.input_file_revision_id,
            input_file_id=self.input_file_id,
            project_id=self.project_id,
            storage_pool_id=self.storage_pool_id,
            revision_number=self.revision_number,
            relative_path=self.relative_path,
            size_bytes=self.size_bytes,
            content_sha256=self.content_sha256,
        )
        if digest != expected:
            raise ValueError("input file revision digest does not match evidence")
        object.__setattr__(self, "created_at", _utc(self.created_at, "created_at"))


def _input_file_revision_digest(
    *,
    input_file_revision_id: str,
    input_file_id: str,
    project_id: str,
    storage_pool_id: str,
    revision_number: int,
    relative_path: str,
    size_bytes: int,
    content_sha256: str,
) -> str:
    return _framed_sha256(
        INPUT_FILE_REVISION_DIGEST_SCHEME,
        input_file_revision_id,
        input_file_id,
        project_id,
        storage_pool_id,
        str(revision_number),
        relative_path,
        str(size_bytes),
        content_sha256,
    )


def build_input_file_revision(
    *,
    input_file_revision_id: str,
    input_file_id: str,
    project_id: str,
    storage_pool_id: str,
    revision_number: int,
    relative_path: str,
    size_bytes: int,
    content_sha256: str,
    created_at: datetime,
) -> InputFileRevision:
    """Build an immutable revision and its deterministic evidence digest."""
    digest = _input_file_revision_digest(
        input_file_revision_id=input_file_revision_id,
        input_file_id=input_file_id,
        project_id=project_id,
        storage_pool_id=storage_pool_id,
        revision_number=revision_number,
        relative_path=relative_path,
        size_bytes=size_bytes,
        content_sha256=content_sha256,
    )
    return InputFileRevision(
        input_file_revision_id=input_file_revision_id,
        input_file_id=input_file_id,
        project_id=project_id,
        storage_pool_id=storage_pool_id,
        revision_number=revision_number,
        relative_path=relative_path,
        size_bytes=size_bytes,
        content_sha256=content_sha256,
        digest_scheme=INPUT_FILE_REVISION_DIGEST_SCHEME,
        digest=digest,
        created_at=created_at,
    )


@dataclass(frozen=True)
class InputFileRevisionSelection:
    """Safe public selection of opaque revisions for one declared input use."""

    input_use_key: str
    occurrence: int
    input_file_revision_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "input_use_key",
            _safe_key(self.input_use_key, "input_use_key"),
        )
        object.__setattr__(
            self,
            "occurrence",
            _non_negative_integer(self.occurrence, "occurrence"),
        )
        revisions = tuple(self.input_file_revision_ids)
        if not revisions:
            raise ValueError("selection must contain at least one input file revision")
        if len(revisions) > _MAX_INPUT_COLLECTION_ITEMS:
            raise ValueError("selection must contain at most 256 input file revisions")
        for revision_id in revisions:
            _validate_id(
                revision_id,
                _INPUT_FILE_REVISION_ID,
                "input_file_revision_id",
            )
        if len(set(revisions)) != len(revisions):
            raise ValueError("selection contains duplicate input file revisions")
        object.__setattr__(self, "input_file_revision_ids", revisions)


@dataclass(frozen=True)
class InputUseDeclaration:
    """Adapter-owned declaration for one exact execution input use."""

    key: str
    occurrence: int
    capability_version: str
    closure_contract_version: str
    allowed_provenance_modes: tuple[InputProvenanceMode | str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "key", _safe_key(self.key, "input use key"))
        object.__setattr__(
            self,
            "occurrence",
            _non_negative_integer(self.occurrence, "occurrence"),
        )
        object.__setattr__(
            self,
            "capability_version",
            _contract_token(self.capability_version, "capability_version"),
        )
        object.__setattr__(
            self,
            "closure_contract_version",
            _contract_token(
                self.closure_contract_version,
                "closure_contract_version",
            ),
        )
        try:
            modes = tuple(
                mode
                if isinstance(mode, InputProvenanceMode)
                else InputProvenanceMode(mode)
                for mode in self.allowed_provenance_modes
            )
        except (TypeError, ValueError):
            raise ValueError("allowed provenance mode is unsupported") from None
        if not modes:
            raise ValueError("at least one provenance mode must be allowed")
        if len(set(modes)) != len(modes):
            raise ValueError("allowed provenance modes contain a duplicate")
        object.__setattr__(self, "allowed_provenance_modes", modes)


@dataclass(frozen=True)
class AdapterInputUseContract:
    """Exact adapter contract for its ordered execution-required input uses."""

    adapter_contract_version: str
    declarations: tuple[InputUseDeclaration, ...]
    allows_mixed: bool

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "adapter_contract_version",
            _contract_token(
                self.adapter_contract_version,
                "adapter_contract_version",
            ),
        )
        declarations = tuple(self.declarations)
        if len(declarations) > _MAX_INPUT_COLLECTION_ITEMS:
            raise ValueError("declarations must contain at most 256 input uses")
        if any(not isinstance(item, InputUseDeclaration) for item in declarations):
            raise ValueError("declarations must contain InputUseDeclaration values")
        coordinates = tuple((item.key, item.occurrence) for item in declarations)
        if len(set(coordinates)) != len(coordinates):
            raise ValueError("input use declarations contain a duplicate coordinate")
        if not isinstance(self.allows_mixed, bool):
            raise ValueError("allows_mixed must be a boolean")
        object.__setattr__(self, "declarations", declarations)


@dataclass(frozen=True)
class PlannedInputUse:
    """Server-side plan for resolving one adapter-declared input use."""

    key: str
    occurrence: int
    capability_version: str
    closure_contract_version: str
    provenance_mode: InputProvenanceMode | str
    input_file_revision_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "key", _safe_key(self.key, "input use key"))
        object.__setattr__(
            self,
            "occurrence",
            _non_negative_integer(self.occurrence, "occurrence"),
        )
        object.__setattr__(
            self,
            "capability_version",
            _contract_token(self.capability_version, "capability_version"),
        )
        closure = _contract_token(
            self.closure_contract_version,
            "closure_contract_version",
        )
        object.__setattr__(self, "closure_contract_version", closure)
        try:
            mode = (
                self.provenance_mode
                if isinstance(self.provenance_mode, InputProvenanceMode)
                else InputProvenanceMode(self.provenance_mode)
            )
        except (TypeError, ValueError):
            raise ValueError("input provenance mode is unsupported") from None
        revisions = tuple(self.input_file_revision_ids)
        if len(revisions) > _MAX_INPUT_COLLECTION_ITEMS:
            raise ValueError(
                "planned use must contain at most 256 input file revisions"
            )
        for revision_id in revisions:
            _validate_id(
                revision_id,
                _INPUT_FILE_REVISION_ID,
                "input_file_revision_id",
            )
        if len(set(revisions)) != len(revisions):
            raise ValueError("planned use contains duplicate input file revisions")
        if mode is InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1 and revisions:
            raise ValueError("transitional input use cannot select managed revisions")
        if (
            mode is InputProvenanceMode.MANAGED_REVISION_V1
            and closure == "regular_file_v1"
            and len(revisions) != 1
        ):
            raise ValueError(
                "managed regular_file_v1 use requires exactly one revision"
            )
        object.__setattr__(self, "provenance_mode", mode)
        object.__setattr__(self, "input_file_revision_ids", revisions)


@dataclass(frozen=True)
class InputUseBindingPlan:
    """Validated ordered plan awaiting repository authorization and resolution."""

    project_id: str
    workflow_id: str
    adapter_contract_version: str
    input_uses: tuple[PlannedInputUse, ...]

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "project_id",
            _validate_id(self.project_id, _PROJECT_ID, "project_id"),
        )
        object.__setattr__(
            self,
            "workflow_id",
            _workflow_id(self.workflow_id),
        )
        object.__setattr__(
            self,
            "adapter_contract_version",
            _contract_token(
                self.adapter_contract_version,
                "adapter_contract_version",
            ),
        )
        uses = tuple(self.input_uses)
        if len(uses) > _MAX_INPUT_COLLECTION_ITEMS:
            raise ValueError("input_uses must contain at most 256 input uses")
        if any(not isinstance(item, PlannedInputUse) for item in uses):
            raise ValueError("input_uses must contain PlannedInputUse values")
        coordinates = tuple((item.key, item.occurrence) for item in uses)
        if len(set(coordinates)) != len(coordinates):
            raise ValueError("planned input uses contain a duplicate coordinate")
        object.__setattr__(self, "input_uses", uses)

    @property
    def fully_managed(self) -> bool:
        """Whether every declared execution input use has managed provenance."""
        return bool(self.input_uses) and all(
            item.provenance_mode is InputProvenanceMode.MANAGED_REVISION_V1
            for item in self.input_uses
        )


def build_input_use_binding_plan(
    *,
    project_id: str,
    workflow_id: str,
    contract: AdapterInputUseContract,
    selections: tuple[InputFileRevisionSelection, ...],
) -> InputUseBindingPlan:
    """Match opaque public selections to an exact ordered adapter contract."""
    if not isinstance(contract, AdapterInputUseContract):
        raise ValueError("contract must be an AdapterInputUseContract")
    selected = tuple(selections)
    if len(selected) > _MAX_INPUT_COLLECTION_ITEMS:
        raise ValueError("selections must contain at most 256 input uses")
    if any(not isinstance(item, InputFileRevisionSelection) for item in selected):
        raise ValueError("selections must contain InputFileRevisionSelection values")
    coordinates = tuple((item.input_use_key, item.occurrence) for item in selected)
    if len(set(coordinates)) != len(coordinates):
        raise ValueError("selections contain a duplicate input use coordinate")
    selected_by_coordinate = {
        (item.input_use_key, item.occurrence): item for item in selected
    }
    declared_coordinates = {
        (item.key, item.occurrence) for item in contract.declarations
    }
    if set(selected_by_coordinate) - declared_coordinates:
        raise ValueError("selection references an undeclared input use")

    planned: list[PlannedInputUse] = []
    for declaration in contract.declarations:
        selection = selected_by_coordinate.get(
            (declaration.key, declaration.occurrence)
        )
        mode = (
            InputProvenanceMode.MANAGED_REVISION_V1
            if selection is not None
            else InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1
        )
        if mode not in declaration.allowed_provenance_modes:
            raise ValueError("input use does not allow the requested provenance mode")
        planned.append(
            PlannedInputUse(
                key=declaration.key,
                occurrence=declaration.occurrence,
                capability_version=declaration.capability_version,
                closure_contract_version=declaration.closure_contract_version,
                provenance_mode=mode,
                input_file_revision_ids=(
                    () if selection is None else selection.input_file_revision_ids
                ),
            )
        )
    modes = {item.provenance_mode for item in planned}
    if len(modes) > 1 and not contract.allows_mixed:
        raise ValueError("adapter contract does not allow mixed provenance")
    return InputUseBindingPlan(
        project_id=project_id,
        workflow_id=workflow_id,
        adapter_contract_version=contract.adapter_contract_version,
        input_uses=tuple(planned),
    )


@dataclass(frozen=True)
class InputFileRevisionBindingRef:
    """Path-free immutable member evidence for a managed input closure."""

    logical_member_key: str
    input_file_id: str
    input_file_revision_id: str
    revision_digest: str
    size_bytes: int
    content_sha256: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "logical_member_key",
            _safe_key(self.logical_member_key, "logical_member_key"),
        )
        _validate_id(self.input_file_id, _INPUT_FILE_ID, "input_file_id")
        _validate_id(
            self.input_file_revision_id,
            _INPUT_FILE_REVISION_ID,
            "input_file_revision_id",
        )
        _validate_digest(self.revision_digest, "revision_digest")
        _non_negative_integer(self.size_bytes, "size_bytes")
        _validate_digest(self.content_sha256, "content_sha256")


def build_input_closure_digest(
    *,
    closure_contract_version: str,
    members: tuple[InputFileRevisionBindingRef, ...],
) -> str:
    """Build generic ordered closure evidence without interpreting its science."""
    closure_contract_version = _contract_token(
        closure_contract_version,
        "closure_contract_version",
    )
    normalized = tuple(members)
    if len(normalized) > _MAX_INPUT_COLLECTION_ITEMS:
        raise ValueError("members must contain at most 256 input file revisions")
    if any(not isinstance(item, InputFileRevisionBindingRef) for item in normalized):
        raise ValueError("members must contain InputFileRevisionBindingRef values")
    canonical = json.dumps(
        [
            {
                "content_sha256": member.content_sha256,
                "input_file_id": member.input_file_id,
                "input_file_revision_id": member.input_file_revision_id,
                "logical_member_key": member.logical_member_key,
                "revision_digest": member.revision_digest,
                "size_bytes": member.size_bytes,
            }
            for member in normalized
        ],
        sort_keys=True,
        separators=(",", ":"),
    )
    return _framed_sha256(
        INPUT_CLOSURE_DIGEST_SCHEME,
        closure_contract_version,
        canonical,
    )


@dataclass(frozen=True)
class InputUseBinding:
    """Resolved path-free evidence for one exact adapter input use."""

    key: str
    occurrence: int
    capability_version: str
    closure_contract_version: str
    provenance_mode: InputProvenanceMode | str
    members: tuple[InputFileRevisionBindingRef, ...]
    closure_digest_scheme: str | None
    closure_digest: str | None

    def __post_init__(self) -> None:
        object.__setattr__(self, "key", _safe_key(self.key, "input use key"))
        object.__setattr__(
            self,
            "occurrence",
            _non_negative_integer(self.occurrence, "occurrence"),
        )
        object.__setattr__(
            self,
            "capability_version",
            _contract_token(self.capability_version, "capability_version"),
        )
        closure_contract_version = _contract_token(
            self.closure_contract_version,
            "closure_contract_version",
        )
        object.__setattr__(
            self,
            "closure_contract_version",
            closure_contract_version,
        )
        try:
            mode = (
                self.provenance_mode
                if isinstance(self.provenance_mode, InputProvenanceMode)
                else InputProvenanceMode(self.provenance_mode)
            )
        except (TypeError, ValueError):
            raise ValueError("input provenance mode is unsupported") from None
        members = tuple(self.members)
        if len(members) > _MAX_INPUT_COLLECTION_ITEMS:
            raise ValueError("members must contain at most 256 input file revisions")
        if any(not isinstance(item, InputFileRevisionBindingRef) for item in members):
            raise ValueError("members must contain InputFileRevisionBindingRef values")
        coordinates = tuple(item.logical_member_key for item in members)
        revision_ids = tuple(item.input_file_revision_id for item in members)
        if len(set(coordinates)) != len(coordinates):
            raise ValueError("input closure contains a duplicate logical member")
        if len(set(revision_ids)) != len(revision_ids):
            raise ValueError("input closure contains a duplicate revision")

        if mode is InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1:
            if (
                members
                or self.closure_digest_scheme is not None
                or self.closure_digest is not None
            ):
                raise ValueError(
                    "transitional input use cannot contain members or closure evidence"
                )
        else:
            if closure_contract_version == "regular_file_v1" and len(members) != 1:
                raise ValueError(
                    "managed regular_file_v1 use requires exactly one member"
                )
            if self.closure_digest_scheme != INPUT_CLOSURE_DIGEST_SCHEME:
                raise ValueError("input closure digest scheme is unsupported")
            digest = _validate_digest(self.closure_digest, "closure_digest")  # type: ignore[arg-type]
            expected = build_input_closure_digest(
                closure_contract_version=closure_contract_version,
                members=members,
            )
            if digest != expected:
                raise ValueError("closure digest does not match ordered members")
        object.__setattr__(self, "provenance_mode", mode)
        object.__setattr__(self, "members", members)


def build_input_use_binding(
    planned_use: PlannedInputUse,
    *,
    members: tuple[InputFileRevisionBindingRef, ...],
) -> InputUseBinding:
    """Resolve one plan to path-free evidence after repository authorization."""
    if not isinstance(planned_use, PlannedInputUse):
        raise ValueError("planned_use must be a PlannedInputUse")
    members = tuple(members)
    if len(members) > _MAX_INPUT_COLLECTION_ITEMS:
        raise ValueError("members must contain at most 256 input file revisions")
    if tuple(member.input_file_revision_id for member in members) != (
        planned_use.input_file_revision_ids
    ):
        raise ValueError("resolved members do not match the planned revisions")
    if planned_use.provenance_mode is InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1:
        return InputUseBinding(
            key=planned_use.key,
            occurrence=planned_use.occurrence,
            capability_version=planned_use.capability_version,
            closure_contract_version=planned_use.closure_contract_version,
            provenance_mode=planned_use.provenance_mode,
            members=(),
            closure_digest_scheme=None,
            closure_digest=None,
        )
    return InputUseBinding(
        key=planned_use.key,
        occurrence=planned_use.occurrence,
        capability_version=planned_use.capability_version,
        closure_contract_version=planned_use.closure_contract_version,
        provenance_mode=planned_use.provenance_mode,
        members=members,
        closure_digest_scheme=INPUT_CLOSURE_DIGEST_SCHEME,
        closure_digest=build_input_closure_digest(
            closure_contract_version=planned_use.closure_contract_version,
            members=members,
        ),
    )


def _canonical_input_binding(
    *,
    project_id: str,
    project_sample_binding_digest: str,
    workflow_id: str,
    adapter_contract_version: str | None,
    workflow_inputs_digest: str,
    contract_mode: InputBindingContractMode,
    input_uses: tuple[InputUseBinding, ...],
) -> str:
    return json.dumps(
        {
            "adapter_contract_version": adapter_contract_version,
            "contract_mode": contract_mode.value,
            "input_uses": [
                {
                    "capability_version": item.capability_version,
                    "closure_contract_version": item.closure_contract_version,
                    "closure_digest": item.closure_digest,
                    "closure_digest_scheme": item.closure_digest_scheme,
                    "key": item.key,
                    "members": [
                        {
                            "content_sha256": member.content_sha256,
                            "input_file_id": member.input_file_id,
                            "input_file_revision_id": member.input_file_revision_id,
                            "logical_member_key": member.logical_member_key,
                            "revision_digest": member.revision_digest,
                            "size_bytes": member.size_bytes,
                        }
                        for member in item.members
                    ],
                    "occurrence": item.occurrence,
                    "provenance_mode": InputProvenanceMode(item.provenance_mode).value,
                }
                for item in input_uses
            ],
            "project_id": project_id,
            "project_sample_binding_digest": project_sample_binding_digest,
            "workflow_id": workflow_id,
            "workflow_inputs_digest": workflow_inputs_digest,
        },
        sort_keys=True,
        separators=(",", ":"),
    )


@dataclass(frozen=True)
class InputUseBindingEnvelope:
    """Immutable digest-pinned input provenance copied from snapshot to run."""

    project_id: str
    project_sample_binding_digest: str
    workflow_id: str
    adapter_contract_version: str | None
    workflow_inputs_digest: str
    contract_mode: InputBindingContractMode | str
    input_uses: tuple[InputUseBinding, ...]
    digest_scheme: str
    digest: str

    def __post_init__(self) -> None:
        project_id = _validate_id(self.project_id, _PROJECT_ID, "project_id")
        project_sample_binding_digest = _validate_digest(
            self.project_sample_binding_digest,
            "project_sample_binding_digest",
        )
        workflow_id = _workflow_id(self.workflow_id)
        workflow_inputs_digest = _validate_digest(
            self.workflow_inputs_digest,
            "workflow_inputs_digest",
        )
        try:
            mode = (
                self.contract_mode
                if isinstance(self.contract_mode, InputBindingContractMode)
                else InputBindingContractMode(self.contract_mode)
            )
        except (TypeError, ValueError):
            raise ValueError("input binding contract mode is unsupported") from None
        uses = tuple(self.input_uses)
        if len(uses) > _MAX_INPUT_COLLECTION_ITEMS:
            raise ValueError("input_uses must contain at most 256 input uses")
        if any(not isinstance(item, InputUseBinding) for item in uses):
            raise ValueError("input_uses must contain InputUseBinding values")
        coordinates = tuple((item.key, item.occurrence) for item in uses)
        if len(set(coordinates)) != len(coordinates):
            raise ValueError("input binding contains a duplicate input use")

        adapter_contract_version = self.adapter_contract_version
        if mode is InputBindingContractMode.COMPATIBILITY_UNRESOLVED_V1:
            if uses:
                raise ValueError("compatibility input binding cannot contain uses")
            if adapter_contract_version is not None:
                adapter_contract_version = _contract_token(
                    adapter_contract_version,
                    "adapter_contract_version",
                )
        else:
            if adapter_contract_version is None:
                raise ValueError("declared input binding requires adapter contract")
            adapter_contract_version = _contract_token(
                adapter_contract_version,
                "adapter_contract_version",
            )
        if self.digest_scheme != INPUT_USE_BINDING_DIGEST_SCHEME:
            raise ValueError("input binding digest scheme is unsupported")
        digest = _validate_digest(self.digest, "digest")
        canonical = _canonical_input_binding(
            project_id=project_id,
            project_sample_binding_digest=project_sample_binding_digest,
            workflow_id=workflow_id,
            adapter_contract_version=adapter_contract_version,
            workflow_inputs_digest=workflow_inputs_digest,
            contract_mode=mode,
            input_uses=uses,
        )
        expected = _framed_sha256(INPUT_USE_BINDING_DIGEST_SCHEME, canonical)
        if digest != expected:
            raise ValueError("input binding digest does not match canonical evidence")
        object.__setattr__(self, "project_id", project_id)
        object.__setattr__(
            self,
            "project_sample_binding_digest",
            project_sample_binding_digest,
        )
        object.__setattr__(self, "workflow_id", workflow_id)
        object.__setattr__(
            self,
            "adapter_contract_version",
            adapter_contract_version,
        )
        object.__setattr__(
            self,
            "workflow_inputs_digest",
            workflow_inputs_digest,
        )
        object.__setattr__(self, "contract_mode", mode)
        object.__setattr__(self, "input_uses", uses)

    @property
    def fully_managed(self) -> bool:
        """Whether this declared envelope has no transitional required use."""
        return (
            self.contract_mode is InputBindingContractMode.DECLARED_INPUT_USES_V1
            and bool(self.input_uses)
            and all(
                item.provenance_mode is InputProvenanceMode.MANAGED_REVISION_V1
                for item in self.input_uses
            )
        )


def _build_input_binding(
    *,
    project_id: str,
    project_sample_binding_digest: str,
    workflow_id: str,
    adapter_contract_version: str | None,
    workflow_inputs_digest: str,
    contract_mode: InputBindingContractMode,
    input_uses: tuple[InputUseBinding, ...],
) -> InputUseBindingEnvelope:
    workflow_id = _workflow_id(workflow_id)
    input_uses = tuple(input_uses)
    if len(input_uses) > _MAX_INPUT_COLLECTION_ITEMS:
        raise ValueError("input_uses must contain at most 256 input uses")
    canonical = _canonical_input_binding(
        project_id=project_id,
        project_sample_binding_digest=project_sample_binding_digest,
        workflow_id=workflow_id,
        adapter_contract_version=adapter_contract_version,
        workflow_inputs_digest=workflow_inputs_digest,
        contract_mode=contract_mode,
        input_uses=input_uses,
    )
    return InputUseBindingEnvelope(
        project_id=project_id,
        project_sample_binding_digest=project_sample_binding_digest,
        workflow_id=workflow_id,
        adapter_contract_version=adapter_contract_version,
        workflow_inputs_digest=workflow_inputs_digest,
        contract_mode=contract_mode,
        input_uses=input_uses,
        digest_scheme=INPUT_USE_BINDING_DIGEST_SCHEME,
        digest=_framed_sha256(INPUT_USE_BINDING_DIGEST_SCHEME, canonical),
    )


def build_input_use_binding_envelope(
    *,
    project_id: str,
    project_sample_binding_digest: str,
    workflow_id: str,
    adapter_contract_version: str,
    workflow_inputs_digest: str,
    input_uses: tuple[InputUseBinding, ...],
) -> InputUseBindingEnvelope:
    """Build declared per-input-use provenance evidence."""
    return _build_input_binding(
        project_id=project_id,
        project_sample_binding_digest=project_sample_binding_digest,
        workflow_id=workflow_id,
        adapter_contract_version=adapter_contract_version,
        workflow_inputs_digest=workflow_inputs_digest,
        contract_mode=InputBindingContractMode.DECLARED_INPUT_USES_V1,
        input_uses=tuple(input_uses),
    )


def build_compatibility_input_binding(
    *,
    project_id: str,
    project_sample_binding_digest: str,
    workflow_id: str,
    adapter_contract_version: str | None,
    workflow_inputs_digest: str,
) -> InputUseBindingEnvelope:
    """Build explicit unresolved evidence without synthesizing input uses."""
    return _build_input_binding(
        project_id=project_id,
        project_sample_binding_digest=project_sample_binding_digest,
        workflow_id=workflow_id,
        adapter_contract_version=adapter_contract_version,
        workflow_inputs_digest=workflow_inputs_digest,
        contract_mode=InputBindingContractMode.COMPATIBILITY_UNRESOLVED_V1,
        input_uses=(),
    )
