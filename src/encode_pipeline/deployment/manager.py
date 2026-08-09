"""Workflow-neutral staging, compatibility, activation, and rollback service."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from pathlib import Path

from encode_pipeline.deployment.admission import (
    DatabaseSchemaObservation,
    DatabaseSchemaObserver,
    DeferredDatabaseSchemaObserver,
    DeferredNativeContractResolver,
    NativeContractResolver,
    ResolvedContractFacts,
    compatible_resolved_facts,
    validate_resolved_facts,
    validate_schema_observation,
)
from encode_pipeline.deployment.bundle import BundleStore
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import (
    COMPONENTS,
    BundleManifest,
    ComponentSlots,
    DeploymentState,
)
from encode_pipeline.deployment.state import StateStore


FaultInjector = Callable[[str], None]


@dataclass(frozen=True)
class DeploymentStatus:
    state: DeploymentState
    manifests: Mapping[str, Mapping[str, BundleManifest | None]]
    partial_staging: tuple[tuple[str, str], ...]
    pending_transactions: tuple[str, ...]
    orphaned_deployments: tuple[tuple[str, str], ...]

    @property
    def interrupted(self) -> bool:
        return bool(
            self.partial_staging
            or self.pending_transactions
            or self.orphaned_deployments
        )

    def to_dict(self) -> dict[str, object]:
        components: dict[str, object] = {}
        for component in COMPONENTS:
            slots = self.state.components[component]
            component_manifests = self.manifests[component]
            values: dict[str, object] = {}
            for slot in ("active", "previous", "staged"):
                identity = getattr(slots, slot)
                manifest = component_manifests[slot]
                values[slot] = (
                    None
                    if identity is None or manifest is None
                    else {"identity": identity}
                )
            components[component] = values
        return {
            "state_identity": self.state.identity,
            "generation": self.state.generation,
            "components": components,
            "interrupted": self.interrupted,
            "partial_staging_count": len(self.partial_staging),
            "pending_transaction_count": len(self.pending_transactions),
            "orphaned_deployment_count": len(self.orphaned_deployments),
        }


@dataclass(frozen=True)
class DeploymentOwnership:
    """Mandatory installed-content ownership for one manager instance."""

    uid: int
    gid: int

    def __post_init__(self) -> None:
        if (
            not isinstance(self.uid, int)
            or isinstance(self.uid, bool)
            or self.uid < 0
            or not isinstance(self.gid, int)
            or isinstance(self.gid, bool)
            or self.gid < 0
        ):
            raise fail(
                "DEPLOYMENT_OWNERSHIP_INVALID",
                "Deployment ownership policy is invalid.",
            )

    @classmethod
    def root(cls) -> "DeploymentOwnership":
        return cls(0, 0)


class DeploymentManager:
    """Coordinate immutable stores; service and migration control stay external."""

    def __init__(
        self,
        layout: DeploymentLayout,
        *,
        ownership: DeploymentOwnership,
        contract_resolver: NativeContractResolver | None = None,
        schema_observer: DatabaseSchemaObserver | None = None,
    ) -> None:
        self.layout = layout
        self.ownership = ownership
        self.bundles = BundleStore(layout)
        self.states = StateStore(layout)
        self.contract_resolver = (
            DeferredNativeContractResolver()
            if contract_resolver is None
            else contract_resolver
        )
        self.schema_observer = (
            DeferredDatabaseSchemaObserver()
            if schema_observer is None
            else schema_observer
        )

    def stage(
        self,
        bundle_path: Path,
        *,
        expected_owner_uid: int | None = None,
        expected_owner_gid: int | None = None,
        fault: FaultInjector | None = None,
    ) -> BundleManifest:
        with self.states.transaction(
            exclusive=True,
            expected_owner_uid=self.ownership.uid,
            expected_owner_gid=self.ownership.gid,
        ) as transaction:
            state = transaction.initialize()
            self._require_clean_state(transaction.pending_transactions())
            manifest = self.bundles.stage(
                bundle_path,
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
                installed_owner_uid=self.ownership.uid,
                installed_owner_gid=self.ownership.gid,
                fault=_prefixed_fault(fault, "bundle"),
            )
            slots = state.components[manifest.component]
            if slots.staged == manifest.identity:
                return manifest
            staged = state.stage(manifest.component, manifest.identity)
            transaction.commit(
                staged,
                operation=f"stage-{manifest.component}",
                expected_current_identity=state.identity,
                fault=_prefixed_fault(fault, "state"),
            )
            return manifest

    def activate(
        self,
        component: str,
        *,
        expected_staged_identity: str,
        fault: FaultInjector | None = None,
    ) -> DeploymentState:
        if component not in COMPONENTS:
            raise fail(
                "DEPLOYMENT_COMPONENT_INVALID",
                "Deployment component is invalid.",
            )
        with self.states.transaction(
            exclusive=True,
            expected_owner_uid=self.ownership.uid,
            expected_owner_gid=self.ownership.gid,
        ) as transaction:
            state = transaction.read()
            self._require_clean_state(transaction.pending_transactions())
            if state.components[component].staged != expected_staged_identity:
                raise fail(
                    "DEPLOYMENT_STAGED_IDENTITY_MISMATCH",
                    "Staged deployment identity does not match the requested deployment.",
                    component=component,
                )
            candidate = state.activate(component)
            self._verify_compatibility(
                candidate,
                base_state=state,
            )
            transaction.commit(
                candidate,
                operation=f"activate-{component}",
                expected_current_identity=state.identity,
                fault=fault,
            )
            return candidate

    def rollback(
        self,
        component: str,
        *,
        expected_previous_identity: str,
        fault: FaultInjector | None = None,
    ) -> DeploymentState:
        if component not in COMPONENTS:
            raise fail(
                "DEPLOYMENT_COMPONENT_INVALID",
                "Deployment component is invalid.",
            )
        with self.states.transaction(
            exclusive=True,
            expected_owner_uid=self.ownership.uid,
            expected_owner_gid=self.ownership.gid,
        ) as transaction:
            state = transaction.read()
            self._require_clean_state(transaction.pending_transactions())
            if state.components[component].previous != expected_previous_identity:
                raise fail(
                    "DEPLOYMENT_PREVIOUS_IDENTITY_MISMATCH",
                    "Previous deployment identity does not match the requested deployment.",
                    component=component,
                )
            candidate = state.rollback(component)
            self._verify_compatibility(
                candidate,
                base_state=state,
            )
            transaction.commit(
                candidate,
                operation=f"rollback-{component}",
                expected_current_identity=state.identity,
                fault=fault,
            )
            return candidate

    def status(self) -> DeploymentStatus:
        """Read one bounded metadata snapshot without hashing release payloads."""
        return self._snapshot(full_verification=False)

    def verify(self) -> DeploymentStatus:
        return self._snapshot(full_verification=True)

    def _snapshot(
        self,
        *,
        full_verification: bool,
    ) -> DeploymentStatus:
        with self.states.transaction(
            exclusive=False,
            expected_owner_uid=self.ownership.uid,
            expected_owner_gid=self.ownership.gid,
        ) as transaction:
            state = transaction.read()
            manifests: dict[str, dict[str, BundleManifest | None]] = {}
            for component in COMPONENTS:
                slots = state.components[component]
                manifests[component] = self._slot_manifests(
                    component,
                    slots,
                    full_verification=full_verification,
                )
            if full_verification:
                self._verify_compatibility(
                    state,
                    base_state=state,
                    verified_manifests=manifests,
                )
            return DeploymentStatus(
                state=state,
                manifests=manifests,
                partial_staging=self.bundles.partial_staging(),
                pending_transactions=transaction.pending_transactions(),
                orphaned_deployments=self._orphaned_deployments(
                    transaction.referenced_identities()
                ),
            )

    def _orphaned_deployments(
        self,
        historical_identities: Mapping[str, frozenset[str]],
    ) -> tuple[tuple[str, str], ...]:
        values: list[tuple[str, str]] = []
        for component in COMPONENTS:
            for identity in self.bundles.installed_identities(component):
                if identity not in historical_identities[component]:
                    values.append((component, identity))
        return tuple(values)

    def _verify_compatibility(
        self,
        state: DeploymentState,
        *,
        base_state: DeploymentState,
        verified_manifests: Mapping[str, Mapping[str, BundleManifest | None]]
        | None = None,
    ) -> None:
        resolved: list[ResolvedContractFacts] = []
        for component in COMPONENTS:
            identity = state.components[component].active
            if identity is None:
                continue
            if verified_manifests is None:
                manifest = self.bundles.verify_installed(
                    component,
                    identity,
                    expected_owner_uid=self.ownership.uid,
                    expected_owner_gid=self.ownership.gid,
                )
            else:
                manifest = verified_manifests[component]["active"]
                if manifest is None or manifest.identity != identity:
                    raise fail(
                        "DEPLOYMENT_RELEASE_INVALID",
                        "Deployment content is invalid.",
                        component=component,
                    )
            resolved.append(self._resolve_contracts(manifest))
        if not compatible_resolved_facts(tuple(resolved)):
            raise fail(
                "DEPLOYMENT_COMPATIBILITY_FAILED",
                "Deployment components are not mutually compatible.",
            )
        if not any(facts.database_heads for facts in resolved):
            return
        heads = self._observe_database_schema(base_state).heads
        for facts in resolved:
            if facts.database_heads and (
                len(heads) != 1 or heads[0] not in facts.database_heads
            ):
                raise fail(
                    "DEPLOYMENT_SCHEMA_INCOMPATIBLE",
                    "Database schema is not compatible with the requested deployment.",
                    component=facts.component,
                )

    def _slot_manifests(
        self,
        component: str,
        slots: ComponentSlots,
        *,
        full_verification: bool,
    ) -> dict[str, BundleManifest | None]:
        values: dict[str, BundleManifest | None] = {}
        for slot in ("active", "previous", "staged"):
            identity = getattr(slots, slot)
            if identity is None:
                values[slot] = None
                continue
            reader = (
                self.bundles.verify_installed
                if full_verification
                else self.bundles.read_installed_manifest
            )
            values[slot] = reader(
                component,
                identity,
                expected_owner_uid=self.ownership.uid,
                expected_owner_gid=self.ownership.gid,
            )
            if full_verification:
                assert values[slot] is not None
                self._resolve_contracts(values[slot])
        return values

    def _resolve_contracts(
        self,
        manifest: BundleManifest,
    ) -> ResolvedContractFacts:
        try:
            facts = self.contract_resolver.resolve(
                self.layout.component_store(manifest.component) / manifest.identity,
                manifest,
            )
        except Exception as error:
            if (
                isinstance(error, DeploymentError)
                and error.issue.code == "DEPLOYMENT_CONTRACT_ADMISSION_DEFERRED"
            ):
                raise
            raise fail(
                "DEPLOYMENT_CONTRACT_ADMISSION_FAILED",
                "Deployment contract admission failed.",
                component=manifest.component,
            ) from None
        validate_resolved_facts(manifest, facts)
        return facts

    def admit_manifest(self, manifest: BundleManifest) -> ResolvedContractFacts:
        """Resolve one snapshot manifest without hashing unrelated payloads."""
        return self._resolve_contracts(manifest)

    def _observe_database_schema(
        self,
        state: DeploymentState,
    ) -> DatabaseSchemaObservation:
        try:
            observation = self.schema_observer.observe(state)
        except Exception as error:
            if (
                isinstance(error, DeploymentError)
                and error.issue.code == "DEPLOYMENT_SCHEMA_OBSERVATION_DEFERRED"
            ):
                raise
            raise fail(
                "DEPLOYMENT_SCHEMA_OBSERVATION_FAILED",
                "Database schema observation failed.",
            ) from None
        validate_schema_observation(state, observation)
        return observation

    @staticmethod
    def _require_clean_state(pending_transactions: tuple[str, ...]) -> None:
        if pending_transactions:
            raise fail(
                "DEPLOYMENT_RECOVERY_REQUIRED",
                "An interrupted deployment operation requires recovery.",
                recoverable=True,
            )


def _prefixed_fault(
    fault: FaultInjector | None,
    prefix: str,
) -> FaultInjector | None:
    if fault is None:
        return None

    def inject(point: str) -> None:
        fault(f"{prefix}:{point}")

    return inject


__all__ = ["DeploymentManager", "DeploymentOwnership", "DeploymentStatus"]
