"""Workflow-neutral Reference Profile identities and immutable evidence."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from hashlib import sha256
import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from encode_pipeline.platform.adapters import WorkflowInputs


ADAPTER_REFERENCE_BINDING_IDENTITY_SCHEME = "sha256-framed-adapter-reference-binding-v1"
REFERENCE_PROFILE_REVISION_IDENTITY_SCHEME = (
    "sha256-framed-reference-profile-revision-v1"
)
REFERENCE_PROFILE_BINDING_DIGEST_SCHEME = "sha256-framed-reference-profile-binding-v1"

_PROFILE_ID = re.compile(r"^refp_[0-9a-f]{32}$")
_REVISION_ID = re.compile(r"^refpr_[0-9a-f]{32}$")
_SAFE_KEY = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]{0,254}$")
_CONTRACT = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:+/-]{0,254}$")
_DIGEST = re.compile(r"^[0-9a-f]{64}$")


def _identifier(value: str, pattern: re.Pattern[str], name: str) -> str:
    if not isinstance(value, str) or pattern.fullmatch(value) is None:
        raise ValueError(f"{name} must be an opaque {name.replace('_', ' ')}")
    return value


def _safe_key(value: str, name: str) -> str:
    if not isinstance(value, str) or _SAFE_KEY.fullmatch(value) is None:
        raise ValueError(f"{name} must be a safe_key")
    return value


def _contract(value: str, name: str) -> str:
    if not isinstance(value, str) or _CONTRACT.fullmatch(value) is None:
        raise ValueError(f"{name} must be a safe versioned contract token")
    return value


def _digest(value: str, name: str) -> str:
    if not isinstance(value, str) or _DIGEST.fullmatch(value) is None:
        raise ValueError(f"{name} must be a lowercase SHA-256 digest")
    return value


def _text(value: str, name: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{name} must be a string")
    normalized = value.strip()
    if (
        not normalized
        or len(normalized) > 255
        or any(ord(character) < 32 or ord(character) == 127 for character in normalized)
    ):
        raise ValueError(f"{name} is invalid")
    return normalized


def _workflow_id(value: str) -> str:
    return _contract(value.strip() if isinstance(value, str) else value, "workflow_id")


def _utc(value: datetime, name: str) -> datetime:
    if (
        not isinstance(value, datetime)
        or value.tzinfo is None
        or value.utcoffset() is None
    ):
        raise ValueError(f"{name} must be timezone-aware")
    return value.astimezone(timezone.utc)


def _positive(value: int, name: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _framed_sha256(*parts: str) -> str:
    digest = sha256()
    for part in parts:
        encoded = part.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return digest.hexdigest()


@dataclass(frozen=True)
class AdapterReferenceBindingIdentity:
    """Adapter-owned identity for one verified private reference binding."""

    workflow_id: str
    contract_version: str
    identity_sha256: str
    identity_scheme: str = ADAPTER_REFERENCE_BINDING_IDENTITY_SCHEME

    def __post_init__(self) -> None:
        object.__setattr__(self, "workflow_id", _workflow_id(self.workflow_id))
        object.__setattr__(
            self,
            "contract_version",
            _contract(self.contract_version, "contract_version"),
        )
        if self.identity_scheme != ADAPTER_REFERENCE_BINDING_IDENTITY_SCHEME:
            raise ValueError("identity_scheme is unsupported")
        object.__setattr__(
            self,
            "identity_sha256",
            _digest(self.identity_sha256, "identity_sha256"),
        )


@dataclass(frozen=True)
class BoundWorkflowReference:
    """Private, process-local adapter binding paired with its public identity."""

    inputs: WorkflowInputs = field(repr=False)
    adapter: object = field(repr=False)
    identity: AdapterReferenceBindingIdentity

    def __post_init__(self) -> None:
        if self.adapter is None:
            raise ValueError("adapter must not be None")
        if not isinstance(self.identity, AdapterReferenceBindingIdentity):
            raise ValueError("identity must be an AdapterReferenceBindingIdentity")


@dataclass(frozen=True)
class ReferenceProfileWorkflowBinding:
    """Path-free stored identity for one revision/workflow compatibility edge."""

    workflow_id: str
    contract_version: str
    identity_sha256: str
    identity_scheme: str = ADAPTER_REFERENCE_BINDING_IDENTITY_SCHEME

    def __post_init__(self) -> None:
        identity = AdapterReferenceBindingIdentity(
            workflow_id=self.workflow_id,
            contract_version=self.contract_version,
            identity_scheme=self.identity_scheme,
            identity_sha256=self.identity_sha256,
        )
        object.__setattr__(self, "workflow_id", identity.workflow_id)
        object.__setattr__(self, "contract_version", identity.contract_version)
        object.__setattr__(self, "identity_scheme", identity.identity_scheme)
        object.__setattr__(self, "identity_sha256", identity.identity_sha256)

    @classmethod
    def from_adapter_identity(
        cls,
        identity: AdapterReferenceBindingIdentity,
    ) -> ReferenceProfileWorkflowBinding:
        if not isinstance(identity, AdapterReferenceBindingIdentity):
            raise ValueError("identity must be an AdapterReferenceBindingIdentity")
        return cls(
            workflow_id=identity.workflow_id,
            contract_version=identity.contract_version,
            identity_scheme=identity.identity_scheme,
            identity_sha256=identity.identity_sha256,
        )

    def to_adapter_identity(self) -> AdapterReferenceBindingIdentity:
        return AdapterReferenceBindingIdentity(
            workflow_id=self.workflow_id,
            contract_version=self.contract_version,
            identity_scheme=self.identity_scheme,
            identity_sha256=self.identity_sha256,
        )


@dataclass(frozen=True)
class ReferenceProfile:
    """Stable logical profile with an optional currently enabled revision."""

    profile_id: str
    safe_key: str
    created_at: datetime
    enabled_revision_id: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "profile_id",
            _identifier(self.profile_id, _PROFILE_ID, "profile_id"),
        )
        object.__setattr__(self, "safe_key", _safe_key(self.safe_key, "safe_key"))
        object.__setattr__(self, "created_at", _utc(self.created_at, "created_at"))
        if self.enabled_revision_id is not None:
            object.__setattr__(
                self,
                "enabled_revision_id",
                _identifier(
                    self.enabled_revision_id,
                    _REVISION_ID,
                    "enabled_revision_id",
                ),
            )


@dataclass(frozen=True)
class ReferenceProfileRevisionSummary:
    """Safe public projection for API, frontend, CLI, and run detail."""

    profile_id: str
    revision_id: str
    safe_key: str
    revision_number: int
    display_name: str
    organism: str
    assembly: str
    public_identity_scheme: str
    public_identity_sha256: str
    compatible_workflow_ids: tuple[str, ...]
    enabled: bool
    created_at: datetime

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "profile_id", _identifier(self.profile_id, _PROFILE_ID, "profile_id")
        )
        object.__setattr__(
            self,
            "revision_id",
            _identifier(self.revision_id, _REVISION_ID, "revision_id"),
        )
        object.__setattr__(self, "safe_key", _safe_key(self.safe_key, "safe_key"))
        object.__setattr__(
            self, "revision_number", _positive(self.revision_number, "revision_number")
        )
        for name in ("display_name", "organism", "assembly"):
            object.__setattr__(self, name, _text(getattr(self, name), name))
        if self.public_identity_scheme != REFERENCE_PROFILE_REVISION_IDENTITY_SCHEME:
            raise ValueError("public_identity_scheme is unsupported")
        object.__setattr__(
            self,
            "public_identity_sha256",
            _digest(self.public_identity_sha256, "public_identity_sha256"),
        )
        workflows = tuple(_workflow_id(item) for item in self.compatible_workflow_ids)
        if workflows != tuple(sorted(set(workflows))):
            raise ValueError("compatible_workflow_ids must be sorted and unique")
        object.__setattr__(self, "compatible_workflow_ids", workflows)
        if not isinstance(self.enabled, bool):
            raise ValueError("enabled must be boolean")
        object.__setattr__(self, "created_at", _utc(self.created_at, "created_at"))

    def to_public_dict(self) -> dict[str, object]:
        return {
            "profile_id": self.profile_id,
            "revision_id": self.revision_id,
            "safe_key": self.safe_key,
            "revision_number": self.revision_number,
            "display_name": self.display_name,
            "organism": self.organism,
            "assembly": self.assembly,
            "public_identity_scheme": self.public_identity_scheme,
            "public_identity_sha256": self.public_identity_sha256,
            "compatible_workflow_ids": list(self.compatible_workflow_ids),
            "enabled": self.enabled,
            "created_at": self.created_at,
        }


@dataclass(frozen=True)
class ReferenceProfileRevision:
    """Append-only metadata and adapter identities for one operator revision."""

    revision_id: str
    profile_id: str
    revision_number: int
    display_name: str
    organism: str
    assembly: str
    config_key: str = field(repr=False)
    workflow_bindings: tuple[ReferenceProfileWorkflowBinding, ...]
    public_identity_scheme: str
    public_identity_sha256: str
    created_at: datetime

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "revision_id",
            _identifier(self.revision_id, _REVISION_ID, "revision_id"),
        )
        object.__setattr__(
            self, "profile_id", _identifier(self.profile_id, _PROFILE_ID, "profile_id")
        )
        object.__setattr__(
            self, "revision_number", _positive(self.revision_number, "revision_number")
        )
        for name in ("display_name", "organism", "assembly"):
            object.__setattr__(self, name, _text(getattr(self, name), name))
        object.__setattr__(self, "config_key", _safe_key(self.config_key, "config_key"))
        bindings = tuple(self.workflow_bindings)
        if not bindings or any(
            not isinstance(binding, ReferenceProfileWorkflowBinding)
            for binding in bindings
        ):
            raise ValueError("workflow_bindings must contain at least one binding")
        workflows = tuple(binding.workflow_id for binding in bindings)
        if len(workflows) != len(set(workflows)):
            raise ValueError("workflow_bindings contain a duplicate workflow")
        bindings = tuple(sorted(bindings, key=lambda item: item.workflow_id))
        object.__setattr__(self, "workflow_bindings", bindings)
        if self.public_identity_scheme != REFERENCE_PROFILE_REVISION_IDENTITY_SCHEME:
            raise ValueError("public_identity_scheme is unsupported")
        expected = build_reference_profile_revision_identity(
            profile_id=self.profile_id,
            revision_number=self.revision_number,
            display_name=self.display_name,
            organism=self.organism,
            assembly=self.assembly,
            workflow_bindings=bindings,
        )
        if self.public_identity_sha256 != expected:
            raise ValueError("public_identity_sha256 does not match revision")
        object.__setattr__(self, "created_at", _utc(self.created_at, "created_at"))

    @property
    def workflow_ids(self) -> tuple[str, ...]:
        return tuple(binding.workflow_id for binding in self.workflow_bindings)

    def binding_for(self, workflow_id: str) -> ReferenceProfileWorkflowBinding:
        normalized = _workflow_id(workflow_id)
        for binding in self.workflow_bindings:
            if binding.workflow_id == normalized:
                return binding
        raise KeyError(normalized)

    def public_summary(
        self,
        *,
        safe_key: str,
        enabled: bool,
    ) -> ReferenceProfileRevisionSummary:
        return ReferenceProfileRevisionSummary(
            profile_id=self.profile_id,
            revision_id=self.revision_id,
            safe_key=safe_key,
            revision_number=self.revision_number,
            display_name=self.display_name,
            organism=self.organism,
            assembly=self.assembly,
            public_identity_scheme=self.public_identity_scheme,
            public_identity_sha256=self.public_identity_sha256,
            compatible_workflow_ids=self.workflow_ids,
            enabled=enabled,
            created_at=self.created_at,
        )


def build_reference_profile_revision_identity(
    *,
    profile_id: str,
    revision_number: int,
    display_name: str,
    organism: str,
    assembly: str,
    workflow_bindings: tuple[ReferenceProfileWorkflowBinding, ...],
) -> str:
    profile_id = _identifier(profile_id, _PROFILE_ID, "profile_id")
    revision_number = _positive(revision_number, "revision_number")
    display_name = _text(display_name, "display_name")
    organism = _text(organism, "organism")
    assembly = _text(assembly, "assembly")
    bindings = tuple(sorted(workflow_bindings, key=lambda item: item.workflow_id))
    parts = [
        REFERENCE_PROFILE_REVISION_IDENTITY_SCHEME,
        profile_id,
        str(revision_number),
        display_name,
        organism,
        assembly,
        str(len(bindings)),
    ]
    for binding in bindings:
        if not isinstance(binding, ReferenceProfileWorkflowBinding):
            raise ValueError("workflow_bindings are invalid")
        parts.extend(
            (
                binding.workflow_id,
                binding.contract_version,
                binding.identity_scheme,
                binding.identity_sha256,
            )
        )
    return _framed_sha256(*parts)


def build_reference_profile_revision(
    *,
    revision_id: str,
    profile_id: str,
    revision_number: int,
    display_name: str,
    organism: str,
    assembly: str,
    config_key: str,
    workflow_bindings: tuple[ReferenceProfileWorkflowBinding, ...],
    created_at: datetime,
) -> ReferenceProfileRevision:
    identity = build_reference_profile_revision_identity(
        profile_id=profile_id,
        revision_number=revision_number,
        display_name=display_name,
        organism=organism,
        assembly=assembly,
        workflow_bindings=workflow_bindings,
    )
    return ReferenceProfileRevision(
        revision_id=revision_id,
        profile_id=profile_id,
        revision_number=revision_number,
        display_name=display_name,
        organism=organism,
        assembly=assembly,
        config_key=config_key,
        workflow_bindings=workflow_bindings,
        public_identity_scheme=REFERENCE_PROFILE_REVISION_IDENTITY_SCHEME,
        public_identity_sha256=identity,
        created_at=created_at,
    )


@dataclass(frozen=True)
class ReferenceProfileRevisionBinding:
    """Path-free exact revision evidence frozen into a snapshot and run."""

    profile_id: str
    revision_id: str
    workflow_id: str
    revision_public_identity_scheme: str
    revision_public_identity_sha256: str
    adapter_contract_version: str
    adapter_identity_scheme: str
    adapter_identity_sha256: str
    binding_digest_scheme: str
    binding_digest: str
    bound_at: datetime

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "profile_id", _identifier(self.profile_id, _PROFILE_ID, "profile_id")
        )
        object.__setattr__(
            self,
            "revision_id",
            _identifier(self.revision_id, _REVISION_ID, "revision_id"),
        )
        object.__setattr__(self, "workflow_id", _workflow_id(self.workflow_id))
        if self.revision_public_identity_scheme != (
            REFERENCE_PROFILE_REVISION_IDENTITY_SCHEME
        ):
            raise ValueError("revision_public_identity_scheme is unsupported")
        object.__setattr__(
            self,
            "revision_public_identity_sha256",
            _digest(
                self.revision_public_identity_sha256,
                "revision_public_identity_sha256",
            ),
        )
        object.__setattr__(
            self,
            "adapter_contract_version",
            _contract(self.adapter_contract_version, "adapter_contract_version"),
        )
        if self.adapter_identity_scheme != ADAPTER_REFERENCE_BINDING_IDENTITY_SCHEME:
            raise ValueError("adapter_identity_scheme is unsupported")
        object.__setattr__(
            self,
            "adapter_identity_sha256",
            _digest(self.adapter_identity_sha256, "adapter_identity_sha256"),
        )
        if self.binding_digest_scheme != REFERENCE_PROFILE_BINDING_DIGEST_SCHEME:
            raise ValueError("binding_digest_scheme is unsupported")
        expected = _framed_sha256(
            REFERENCE_PROFILE_BINDING_DIGEST_SCHEME,
            self.profile_id,
            self.revision_id,
            self.workflow_id,
            self.revision_public_identity_scheme,
            self.revision_public_identity_sha256,
            self.adapter_contract_version,
            self.adapter_identity_scheme,
            self.adapter_identity_sha256,
        )
        if self.binding_digest != expected:
            raise ValueError("binding_digest does not match reference evidence")
        object.__setattr__(self, "bound_at", _utc(self.bound_at, "bound_at"))


def build_reference_profile_revision_binding(
    *,
    profile_id: str,
    revision_id: str,
    workflow_id: str,
    revision_public_identity_sha256: str,
    adapter_identity: AdapterReferenceBindingIdentity,
    bound_at: datetime,
) -> ReferenceProfileRevisionBinding:
    if not isinstance(adapter_identity, AdapterReferenceBindingIdentity):
        raise ValueError("adapter_identity must be an AdapterReferenceBindingIdentity")
    if _workflow_id(workflow_id) != adapter_identity.workflow_id:
        raise ValueError("adapter identity workflow does not match evidence")
    digest = _framed_sha256(
        REFERENCE_PROFILE_BINDING_DIGEST_SCHEME,
        profile_id,
        revision_id,
        workflow_id,
        REFERENCE_PROFILE_REVISION_IDENTITY_SCHEME,
        revision_public_identity_sha256,
        adapter_identity.contract_version,
        adapter_identity.identity_scheme,
        adapter_identity.identity_sha256,
    )
    return ReferenceProfileRevisionBinding(
        profile_id=profile_id,
        revision_id=revision_id,
        workflow_id=workflow_id,
        revision_public_identity_scheme=REFERENCE_PROFILE_REVISION_IDENTITY_SCHEME,
        revision_public_identity_sha256=revision_public_identity_sha256,
        adapter_contract_version=adapter_identity.contract_version,
        adapter_identity_scheme=adapter_identity.identity_scheme,
        adapter_identity_sha256=adapter_identity.identity_sha256,
        binding_digest_scheme=REFERENCE_PROFILE_BINDING_DIGEST_SCHEME,
        binding_digest=digest,
        bound_at=bound_at,
    )


@dataclass(frozen=True)
class ResolvedReferenceProfile:
    """Exact public evidence paired with private process-local adapter inputs."""

    summary: ReferenceProfileRevisionSummary
    evidence: ReferenceProfileRevisionBinding
    bound_reference: BoundWorkflowReference = field(repr=False)

    def __post_init__(self) -> None:
        if not isinstance(self.summary, ReferenceProfileRevisionSummary):
            raise ValueError("summary must be a ReferenceProfileRevisionSummary")
        if not isinstance(self.evidence, ReferenceProfileRevisionBinding):
            raise ValueError("evidence must be a ReferenceProfileRevisionBinding")
        if not isinstance(self.bound_reference, BoundWorkflowReference):
            raise ValueError("bound_reference must be a BoundWorkflowReference")
        if (
            self.summary.profile_id != self.evidence.profile_id
            or self.summary.revision_id != self.evidence.revision_id
            or self.bound_reference.identity.workflow_id != self.evidence.workflow_id
            or self.bound_reference.identity.identity_sha256
            != self.evidence.adapter_identity_sha256
        ):
            raise ValueError("resolved reference identities do not match")
