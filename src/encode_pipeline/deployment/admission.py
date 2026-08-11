"""Finite native-contract resolution seam for deployment bundle facts.

The deployment manifest only binds native contract identities to indexed files.
It is not allowed to restate distribution versions, migration heads, runtime
locks, or reference compatibility.  A component-specific resolver must derive
those facts from their existing authoritative documents before verification or
activation.  The deferred resolver remains an explicit metadata-reader fallback;
bundle production and the fixed operator candidate action inject the production
candidate-byte resolver at their respective trust boundaries.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Literal, Protocol

from encode_pipeline.deployment.canonical import canonical_identity
from encode_pipeline.deployment.errors import fail
from encode_pipeline.deployment.models import (
    PLATFORM,
    BundleManifest,
    ContractIdentity,
    ContractRequirement,
    DeploymentState,
)


_IDENTITY = re.compile(r"^sha256-[0-9a-f]{64}$")
_VERSION = re.compile(r"^[a-z0-9][a-z0-9._+-]{0,127}$")
_HEAD = re.compile(r"^[a-z0-9][a-z0-9._+-]{0,127}$")
DATABASE_SCHEMA_OBSERVATION_IDENTITY_SCHEME = (
    "helixweave-database-schema-observation-identity-v1"
)
CompatibilityStatus = Literal["compatible", "incomplete", "incompatible"]


@dataclass(frozen=True)
class ResolvedContractFacts:
    """Facts re-derived from one component's authoritative documents."""

    component: str
    deployment_identity: str
    version: str
    contracts: tuple[ContractIdentity, ...]
    requirements: tuple[ContractRequirement, ...] = ()
    database_heads: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if (
            _IDENTITY.fullmatch(self.deployment_identity) is None
            or _VERSION.fullmatch(self.version) is None
            or tuple(sorted(set(self.contracts))) != self.contracts
            or tuple(sorted(set(self.requirements))) != self.requirements
            or len({item.contract for item in self.requirements})
            != len(self.requirements)
            or tuple(sorted(set(self.database_heads))) != self.database_heads
            or any(_HEAD.fullmatch(item) is None for item in self.database_heads)
        ):
            raise fail(
                "DEPLOYMENT_CONTRACT_FACTS_INVALID",
                "Deployment contract facts are invalid.",
            )


class NativeContractResolver(Protocol):
    """Resolve only the repository's known native contract types."""

    def resolve(
        self,
        root: Path,
        manifest: BundleManifest,
    ) -> ResolvedContractFacts: ...


class DeferredNativeContractResolver:
    """Explicit fail-closed resolver for incomplete/bootstrap composition."""

    def resolve(
        self,
        root: Path,
        manifest: BundleManifest,
    ) -> ResolvedContractFacts:
        del root, manifest
        raise fail(
            "DEPLOYMENT_CONTRACT_ADMISSION_DEFERRED",
            "Deployment contract admission is not connected.",
            recoverable=True,
        )


@dataclass(frozen=True)
class DatabaseSchemaObservation:
    """Observed database state bound to the state generation being changed."""

    identity: str
    provider_identity: str
    state_identity: str
    database_identity: str
    heads: tuple[str, ...]

    @classmethod
    def create(
        cls,
        *,
        provider_identity: str,
        state_identity: str,
        database_identity: str,
        heads: tuple[str, ...],
    ) -> "DatabaseSchemaObservation":
        value = {
            "provider_identity": provider_identity,
            "state_identity": state_identity,
            "database_identity": database_identity,
            "heads": list(heads),
        }
        return cls(
            identity=canonical_identity(
                value,
                scheme=DATABASE_SCHEMA_OBSERVATION_IDENTITY_SCHEME,
            ),
            provider_identity=provider_identity,
            state_identity=state_identity,
            database_identity=database_identity,
            heads=heads,
        )

    def __post_init__(self) -> None:
        value = {
            "provider_identity": self.provider_identity,
            "state_identity": self.state_identity,
            "database_identity": self.database_identity,
            "heads": list(self.heads),
        }
        if (
            _IDENTITY.fullmatch(self.identity) is None
            or _IDENTITY.fullmatch(self.provider_identity) is None
            or _IDENTITY.fullmatch(self.state_identity) is None
            or _IDENTITY.fullmatch(self.database_identity) is None
            or tuple(sorted(set(self.heads))) != self.heads
            or any(_HEAD.fullmatch(item) is None for item in self.heads)
            or self.identity
            != canonical_identity(
                value,
                scheme=DATABASE_SCHEMA_OBSERVATION_IDENTITY_SCHEME,
            )
        ):
            raise fail(
                "DEPLOYMENT_SCHEMA_OBSERVATION_INVALID",
                "Database schema observation is invalid.",
            )


class DatabaseSchemaObserver(Protocol):
    def observe(self, state: DeploymentState) -> DatabaseSchemaObservation: ...


class DeferredDatabaseSchemaObserver:
    def observe(self, state: DeploymentState) -> DatabaseSchemaObservation:
        del state
        raise fail(
            "DEPLOYMENT_SCHEMA_OBSERVATION_DEFERRED",
            "Database schema observation is not connected.",
            recoverable=True,
        )


def validate_resolved_facts(
    manifest: BundleManifest,
    facts: ResolvedContractFacts,
) -> None:
    """Prove the resolver result is for exactly the indexed deployment."""
    if (
        not isinstance(facts, ResolvedContractFacts)
        or facts.component != manifest.component
        or facts.deployment_identity != manifest.identity
        or facts.contracts != manifest.provides
        or (manifest.component == PLATFORM) != bool(facts.database_heads)
    ):
        raise fail(
            "DEPLOYMENT_CONTRACT_ADMISSION_FAILED",
            "Deployment contract admission failed.",
            component=manifest.component,
        )


def validate_schema_observation(
    state: DeploymentState,
    observation: DatabaseSchemaObservation,
) -> None:
    if (
        not isinstance(observation, DatabaseSchemaObservation)
        or observation.state_identity != state.identity
    ):
        raise fail(
            "DEPLOYMENT_SCHEMA_OBSERVATION_FAILED",
            "Database schema observation does not match deployment state.",
        )


def resolved_facts_compatibility(
    facts: tuple[ResolvedContractFacts, ...],
) -> CompatibilityStatus:
    """Classify only provider/requirement facts derived from native contracts.

    A requirement without a provider is the sole incomplete case.  Once a
    provider is present, duplicate or unaccepted identities are incompatible;
    they can never be treated as partial assembly.
    """

    provided: dict[str, str] = {}
    for component in facts:
        for item in component.contracts:
            if item.contract in provided and provided[item.contract] != item.identity:
                return "incompatible"
            provided[item.contract] = item.identity
    missing_provider = False
    for component in facts:
        for requirement in component.requirements:
            identity = provided.get(requirement.contract)
            if identity is None:
                missing_provider = True
            elif identity not in requirement.accepted_identities:
                return "incompatible"
    return "incomplete" if missing_provider else "compatible"


def compatible_resolved_facts(facts: tuple[ResolvedContractFacts, ...]) -> bool:
    """Compatibility-only wrapper retained for existing non-receipt callers."""

    return resolved_facts_compatibility(facts) == "compatible"


__all__ = [
    "DATABASE_SCHEMA_OBSERVATION_IDENTITY_SCHEME",
    "DatabaseSchemaObservation",
    "DatabaseSchemaObserver",
    "DeferredDatabaseSchemaObserver",
    "DeferredNativeContractResolver",
    "CompatibilityStatus",
    "NativeContractResolver",
    "ResolvedContractFacts",
    "compatible_resolved_facts",
    "resolved_facts_compatibility",
    "validate_resolved_facts",
    "validate_schema_observation",
]
