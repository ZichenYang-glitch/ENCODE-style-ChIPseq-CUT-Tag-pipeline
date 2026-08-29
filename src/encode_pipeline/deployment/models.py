"""Versioned deployment manifest and mutable slot-state contracts.

Release manifests index bytes and point at existing content-addressed scientific,
migration, and frontend contracts. They intentionally do not restate facts owned
by those contracts. Deployment state owns only active, previous, and staged slot
selection.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timezone
import re
from typing import Any

from encode_pipeline.deployment.canonical import canonical_identity, without_key
from encode_pipeline.deployment.errors import fail


BUNDLE_SCHEMA_VERSION = "helixweave-deployment-bundle-v1"
STATE_SCHEMA_VERSION = "helixweave-deployment-state-v1"
BUNDLE_IDENTITY_SCHEME = "helixweave-deployment-bundle-identity-v1"
STATE_IDENTITY_SCHEME = "helixweave-deployment-state-identity-v1"

PLATFORM = "platform"
ENCODE_RUNTIME = "encode-runtime"
BULK_RNASEQ_RUNTIME = "bulk-rnaseq-runtime"
COMPONENTS = (PLATFORM, ENCODE_RUNTIME, BULK_RNASEQ_RUNTIME)
REQUIRED_PROVIDERS = {
    PLATFORM: frozenset(
        {
            "helixweave.platform.distribution",
            "helixweave.platform.python-wheel",
            "helixweave.platform.python-runtime",
            "helixweave.platform.frontend-assets",
            "helixweave.platform.database-migrations",
            "helixweave.platform.reference-compatibility",
        }
    ),
    ENCODE_RUNTIME: frozenset({"helixweave.runtime.encode"}),
    BULK_RNASEQ_RUNTIME: frozenset({"helixweave.runtime.bulk-rnaseq"}),
}

_IDENTITY = re.compile(r"^sha256-[0-9a-f]{64}$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_CONTRACT = re.compile(r"^[a-z0-9][a-z0-9._:/+-]{0,191}$")
_RELATIVE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._/+@=$-]{0,511}$")


def _object(value: object, *, code: str) -> dict[str, Any]:
    if not isinstance(value, Mapping) or any(not isinstance(key, str) for key in value):
        raise fail(code, "Deployment document is invalid.")
    return dict(value)


def _exact_keys(
    value: Mapping[str, object],
    required: set[str],
    *,
    optional: set[str] = frozenset(),
    code: str,
) -> None:
    keys = set(value)
    if not required.issubset(keys) or not keys.issubset(required | optional):
        raise fail(code, "Deployment document is invalid.")


def _string(value: object, pattern: re.Pattern[str], *, code: str) -> str:
    if not isinstance(value, str) or pattern.fullmatch(value) is None:
        raise fail(code, "Deployment document is invalid.")
    return value


def _identity(value: object, *, code: str) -> str:
    return _string(value, _IDENTITY, code=code)


def _relative(value: object, *, code: str) -> str:
    rendered = _string(value, _RELATIVE, code=code)
    parts = rendered.split("/")
    if (
        rendered.startswith("/")
        or rendered.endswith("/")
        or "//" in rendered
        or any(part in {"", ".", ".."} for part in parts)
        or "\\" in rendered
    ):
        raise fail(code, "Deployment document is invalid.")
    return rendered


@dataclass(frozen=True, order=True)
class FileRecord:
    path: str
    size_bytes: int
    sha256: str
    mode: int

    @classmethod
    def from_dict(cls, raw: object) -> "FileRecord":
        code = "DEPLOYMENT_MANIFEST_INVALID"
        value = _object(raw, code=code)
        _exact_keys(
            value,
            {"path", "size_bytes", "sha256", "mode"},
            code=code,
        )
        path = _relative(value["path"], code=code)
        if not path.startswith("payload/"):
            raise fail(code, "Deployment document is invalid.")
        size = value["size_bytes"]
        mode = value["mode"]
        if (
            not isinstance(size, int)
            or isinstance(size, bool)
            or not 0 <= size <= 100 * 1024**3
            or not isinstance(mode, int)
            or isinstance(mode, bool)
            or mode not in {0o444, 0o555}
        ):
            raise fail(code, "Deployment document is invalid.")
        return cls(
            path=path,
            size_bytes=size,
            sha256=_string(value["sha256"], _SHA256, code=code),
            mode=mode,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "path": self.path,
            "size_bytes": self.size_bytes,
            "sha256": self.sha256,
            "mode": self.mode,
        }


@dataclass(frozen=True, order=True)
class ContractIdentity:
    contract: str
    identity: str

    @classmethod
    def from_dict(cls, raw: object) -> "ContractIdentity":
        code = "DEPLOYMENT_MANIFEST_INVALID"
        value = _object(raw, code=code)
        _exact_keys(value, {"contract", "identity"}, code=code)
        return cls(
            contract=_string(value["contract"], _CONTRACT, code=code),
            identity=_identity(value["identity"], code=code),
        )

    def to_dict(self) -> dict[str, str]:
        return {"contract": self.contract, "identity": self.identity}


@dataclass(frozen=True, order=True)
class ContractDocument:
    """One native contract identity bound to one indexed document path."""

    contract: str
    identity: str
    path: str

    @classmethod
    def from_dict(cls, raw: object) -> "ContractDocument":
        code = "DEPLOYMENT_MANIFEST_INVALID"
        value = _object(raw, code=code)
        _exact_keys(value, {"contract", "identity", "path"}, code=code)
        path = _relative(value["path"], code=code)
        if not path.startswith("payload/contracts/"):
            raise fail(code, "Deployment document is invalid.")
        return cls(
            contract=_string(value["contract"], _CONTRACT, code=code),
            identity=_identity(value["identity"], code=code),
            path=path,
        )

    def to_dict(self) -> dict[str, str]:
        return {
            "contract": self.contract,
            "identity": self.identity,
            "path": self.path,
        }


@dataclass(frozen=True, order=True)
class ContractRequirement:
    contract: str
    accepted_identities: tuple[str, ...]

    def __post_init__(self) -> None:
        if (
            not isinstance(self.contract, str)
            or _CONTRACT.fullmatch(self.contract) is None
            or not self.accepted_identities
            or tuple(sorted(set(self.accepted_identities))) != self.accepted_identities
            or any(
                _IDENTITY.fullmatch(item) is None for item in self.accepted_identities
            )
        ):
            raise fail(
                "DEPLOYMENT_CONTRACT_FACTS_INVALID",
                "Deployment contract facts are invalid.",
            )

    @classmethod
    def from_dict(cls, raw: object) -> "ContractRequirement":
        code = "DEPLOYMENT_MANIFEST_INVALID"
        value = _object(raw, code=code)
        _exact_keys(value, {"contract", "accepted_identities"}, code=code)
        raw_identities = value["accepted_identities"]
        if not isinstance(raw_identities, Sequence) or isinstance(
            raw_identities, (str, bytes)
        ):
            raise fail(code, "Deployment document is invalid.")
        identities = tuple(_identity(item, code=code) for item in raw_identities)
        if not identities or tuple(sorted(set(identities))) != identities:
            raise fail(code, "Deployment document is invalid.")
        return cls(
            contract=_string(value["contract"], _CONTRACT, code=code),
            accepted_identities=identities,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "contract": self.contract,
            "accepted_identities": list(self.accepted_identities),
        }


@dataclass(frozen=True)
class BundleManifest:
    identity: str
    component: str
    contracts: tuple[ContractDocument, ...]
    files: tuple[FileRecord, ...]

    @classmethod
    def from_dict(cls, raw: object) -> "BundleManifest":
        code = "DEPLOYMENT_MANIFEST_INVALID"
        value = _object(raw, code=code)
        _exact_keys(
            value,
            {
                "schema_version",
                "identity",
                "component",
                "contracts",
                "files",
            },
            code=code,
        )
        if value["schema_version"] != BUNDLE_SCHEMA_VERSION:
            raise fail(code, "Deployment document is invalid.")
        component = value["component"]
        if component not in COMPONENTS:
            raise fail(code, "Deployment document is invalid.")

        contracts = _ordered_records(value["contracts"], ContractDocument, code=code)
        files = _ordered_records(value["files"], FileRecord, code=code)
        file_paths = {item.path for item in files}
        casefolded_paths = {item.casefold() for item in file_paths}
        contract_paths = {item.path for item in contracts}
        if (
            not files
            or len(file_paths) != len(files)
            or not contracts
            or len(contract_paths) != len(contracts)
            or not contract_paths.issubset(file_paths)
            or len(casefolded_paths) != len(file_paths)
        ):
            raise fail(code, "Deployment document is invalid.")

        observed_identity = _identity(value["identity"], code=code)
        expected_identity = canonical_identity(
            without_key(value, "identity"), scheme=BUNDLE_IDENTITY_SCHEME
        )
        if observed_identity != expected_identity:
            raise fail(code, "Deployment document identity is invalid.")

        manifest = cls(
            identity=observed_identity,
            component=component,
            contracts=contracts,
            files=files,
        )
        manifest._verify_contract_uniqueness()
        provided_contracts = {item.contract for item in manifest.provides}
        if provided_contracts != REQUIRED_PROVIDERS[component]:
            raise fail(code, "Deployment document is invalid.")
        return manifest

    @classmethod
    def create(
        cls,
        *,
        component: str,
        contracts: Sequence[ContractDocument],
        files: Sequence[FileRecord],
    ) -> "BundleManifest":
        value: dict[str, object] = {
            "schema_version": BUNDLE_SCHEMA_VERSION,
            "component": component,
            "contracts": [item.to_dict() for item in sorted(contracts)],
            "files": [item.to_dict() for item in sorted(files)],
        }
        value["identity"] = canonical_identity(value, scheme=BUNDLE_IDENTITY_SCHEME)
        return cls.from_dict(value)

    def _verify_contract_uniqueness(self) -> None:
        provided = [item.contract for item in self.provides]
        if len(provided) != len(set(provided)):
            raise fail("DEPLOYMENT_MANIFEST_INVALID", "Deployment document is invalid.")

    @property
    def provides(self) -> tuple[ContractIdentity, ...]:
        return tuple(
            ContractIdentity(item.contract, item.identity) for item in self.contracts
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": BUNDLE_SCHEMA_VERSION,
            "identity": self.identity,
            "component": self.component,
            "contracts": [item.to_dict() for item in self.contracts],
            "files": [item.to_dict() for item in self.files],
        }


def _ordered_records(value: object, record_type, *, code: str) -> tuple[Any, ...]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise fail(code, "Deployment document is invalid.")
    items = tuple(record_type.from_dict(item) for item in value)
    if tuple(sorted(set(items))) != items:
        raise fail(code, "Deployment document is invalid.")
    return items


@dataclass(frozen=True)
class ComponentSlots:
    active: str | None = None
    previous: str | None = None
    staged: str | None = None

    @classmethod
    def from_dict(cls, raw: object) -> "ComponentSlots":
        code = "DEPLOYMENT_STATE_INVALID"
        value = _object(raw, code=code)
        _exact_keys(value, {"active", "previous", "staged"}, code=code)
        identities: list[str | None] = []
        for key in ("active", "previous", "staged"):
            item = value[key]
            identities.append(None if item is None else _identity(item, code=code))
        nonempty = [item for item in identities if item is not None]
        if len(nonempty) != len(set(nonempty)):
            raise fail(code, "Deployment state is invalid.")
        return cls(*identities)

    def to_dict(self) -> dict[str, str | None]:
        return {
            "active": self.active,
            "previous": self.previous,
            "staged": self.staged,
        }


@dataclass(frozen=True)
class DeploymentState:
    identity: str
    generation: int
    created_at: str
    components: Mapping[str, ComponentSlots]

    @classmethod
    def initial(cls) -> "DeploymentState":
        return cls.create(
            generation=0,
            created_at=datetime.now(timezone.utc).isoformat(),
            components={component: ComponentSlots() for component in COMPONENTS},
        )

    @classmethod
    def create(
        cls,
        *,
        generation: int,
        created_at: str,
        components: Mapping[str, ComponentSlots],
    ) -> "DeploymentState":
        value: dict[str, object] = {
            "schema_version": STATE_SCHEMA_VERSION,
            "generation": generation,
            "created_at": created_at,
            "components": {
                component: components[component].to_dict() for component in COMPONENTS
            },
        }
        value["identity"] = canonical_identity(value, scheme=STATE_IDENTITY_SCHEME)
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "DeploymentState":
        code = "DEPLOYMENT_STATE_INVALID"
        value = _object(raw, code=code)
        _exact_keys(
            value,
            {"schema_version", "identity", "generation", "created_at", "components"},
            code=code,
        )
        if value["schema_version"] != STATE_SCHEMA_VERSION:
            raise fail(code, "Deployment state is invalid.")
        generation = value["generation"]
        if (
            not isinstance(generation, int)
            or isinstance(generation, bool)
            or generation < 0
        ):
            raise fail(code, "Deployment state is invalid.")
        created_at = value["created_at"]
        if not isinstance(created_at, str) or len(created_at) > 64:
            raise fail(code, "Deployment state is invalid.")
        try:
            timestamp = datetime.fromisoformat(created_at)
        except ValueError:
            raise fail(code, "Deployment state is invalid.") from None
        if timestamp.tzinfo is None or timestamp.utcoffset() is None:
            raise fail(code, "Deployment state is invalid.")

        raw_components = _object(value["components"], code=code)
        if set(raw_components) != set(COMPONENTS):
            raise fail(code, "Deployment state is invalid.")
        components = {
            component: ComponentSlots.from_dict(raw_components[component])
            for component in COMPONENTS
        }
        observed_identity = _identity(value["identity"], code=code)
        if observed_identity != canonical_identity(
            without_key(value, "identity"), scheme=STATE_IDENTITY_SCHEME
        ):
            raise fail(code, "Deployment state identity is invalid.")
        return cls(
            identity=observed_identity,
            generation=generation,
            created_at=created_at,
            components=components,
        )

    def stage(self, component: str, identity: str) -> "DeploymentState":
        _require_component(component)
        _identity(identity, code="DEPLOYMENT_STATE_INVALID")
        slots = self.components[component]
        if identity in {slots.active, slots.previous}:
            raise fail(
                "DEPLOYMENT_STAGE_CONFLICT",
                "Deployment identity is already active or previous.",
                component=component,
            )
        return self._replace(
            component, ComponentSlots(slots.active, slots.previous, identity)
        )

    def activate(self, component: str) -> "DeploymentState":
        _require_component(component)
        slots = self.components[component]
        if slots.staged is None:
            raise fail(
                "DEPLOYMENT_STAGE_MISSING",
                "No verified staged deployment is available.",
                component=component,
            )
        return self._replace(
            component,
            ComponentSlots(
                active=slots.staged,
                previous=slots.active,
                staged=None,
            ),
        )

    def rollback(self, component: str) -> "DeploymentState":
        _require_component(component)
        slots = self.components[component]
        if slots.previous is None:
            raise fail(
                "DEPLOYMENT_PREVIOUS_MISSING",
                "No compatible previous deployment is available.",
                component=component,
            )
        return self._replace(
            component,
            ComponentSlots(
                active=slots.previous,
                previous=slots.active,
                staged=slots.staged,
            ),
        )

    def _replace(self, component: str, slots: ComponentSlots) -> "DeploymentState":
        components = dict(self.components)
        components[component] = slots
        return self.create(
            generation=self.generation + 1,
            created_at=datetime.now(timezone.utc).isoformat(),
            components=components,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": STATE_SCHEMA_VERSION,
            "identity": self.identity,
            "generation": self.generation,
            "created_at": self.created_at,
            "components": {
                component: self.components[component].to_dict()
                for component in COMPONENTS
            },
        }


def _require_component(component: str) -> None:
    if component not in COMPONENTS:
        raise fail("DEPLOYMENT_COMPONENT_INVALID", "Deployment component is invalid.")


__all__ = [
    "BULK_RNASEQ_RUNTIME",
    "BUNDLE_SCHEMA_VERSION",
    "BundleManifest",
    "COMPONENTS",
    "ComponentSlots",
    "ContractDocument",
    "ContractIdentity",
    "ContractRequirement",
    "DeploymentState",
    "ENCODE_RUNTIME",
    "FileRecord",
    "PLATFORM",
    "STATE_SCHEMA_VERSION",
]
