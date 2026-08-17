"""Fixed supported host layout and test-only isolated layout construction."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from encode_pipeline.deployment.errors import fail
from encode_pipeline.deployment.models import (
    BULK_RNASEQ_RUNTIME,
    ENCODE_RUNTIME,
    PLATFORM,
)


@dataclass(frozen=True)
class DeploymentLayout:
    immutable_root: Path
    configuration_root: Path
    data_root: Path
    log_root: Path
    run_root: Path

    def __post_init__(self) -> None:
        roots = (
            self.immutable_root,
            self.configuration_root,
            self.data_root,
            self.log_root,
            self.run_root,
        )
        if any(not isinstance(root, Path) or not root.is_absolute() for root in roots):
            raise fail("DEPLOYMENT_LAYOUT_INVALID", "Deployment layout is invalid.")
        rendered = [str(root) for root in roots]
        if len(rendered) != len(set(rendered)) or any(
            any(part in {"", ".", ".."} for part in root.parts[1:])
            or "\x00" in str(root)
            or "\n" in str(root)
            or "\r" in str(root)
            for root in roots
        ):
            raise fail("DEPLOYMENT_LAYOUT_INVALID", "Deployment layout is invalid.")

    @classmethod
    def supported(cls) -> "DeploymentLayout":
        return cls(
            immutable_root=Path("/opt/helixweave"),
            configuration_root=Path("/etc/helixweave"),
            data_root=Path("/var/lib/helixweave"),
            log_root=Path("/var/log/helixweave"),
            run_root=Path("/run/helixweave"),
        )

    @classmethod
    def isolated(cls, root: Path) -> "DeploymentLayout":
        if not isinstance(root, Path) or not root.is_absolute():
            raise fail("DEPLOYMENT_LAYOUT_INVALID", "Deployment layout is invalid.")
        return cls(
            immutable_root=root / "opt",
            configuration_root=root / "etc",
            data_root=root / "var" / "lib",
            log_root=root / "var" / "log",
            run_root=root / "run",
        )

    @property
    def platform_releases(self) -> Path:
        return self.immutable_root / "releases" / "platform"

    @property
    def encode_runtimes(self) -> Path:
        return self.immutable_root / "runtimes" / "encode"

    @property
    def encode_runtime_materialized(self) -> Path:
        return self.data_root / "scientific" / "encode"

    @property
    def bulk_rnaseq_runtimes(self) -> Path:
        return self.immutable_root / "runtimes" / "bulk-rnaseq"

    def component_store(self, component: str) -> Path:
        stores = {
            PLATFORM: self.platform_releases,
            ENCODE_RUNTIME: self.encode_runtimes,
            BULK_RNASEQ_RUNTIME: self.bulk_rnaseq_runtimes,
        }
        try:
            return stores[component]
        except KeyError:
            raise fail(
                "DEPLOYMENT_COMPONENT_INVALID", "Deployment component is invalid."
            ) from None

    @property
    def state_generations(self) -> Path:
        return self.data_root / "deployment" / "generations"

    @property
    def current_state_link(self) -> Path:
        return self.data_root / "deployment" / "current"

    @property
    def state_lock(self) -> Path:
        return self.data_root / "operator" / "state" / "operation.lock"

    @property
    def state_transactions(self) -> Path:
        return self.data_root / "operator" / "state" / "transactions"

    @property
    def operator_transactions(self) -> Path:
        return self.data_root / "operator" / "transactions"

    @property
    def operator_transaction_active(self) -> Path:
        return self.operator_transactions / "active.json"

    @property
    def operator_transaction_history(self) -> Path:
        return self.operator_transactions / "history"

    @property
    def operator_transaction_lock(self) -> Path:
        return self.operator_transactions / "operator.lock"

    @property
    def operator_action_root(self) -> Path:
        return self.run_root / "operator" / "action"

    @property
    def operator_action_request(self) -> Path:
        return self.operator_action_root / "request.json"

    @property
    def operator_action_receipt(self) -> Path:
        return self.operator_action_root / "receipt.json"

    @property
    def database_prepare_root(self) -> Path:
        return self.run_root / "database"

    @property
    def database_prepare_request(self) -> Path:
        return self.database_prepare_root / "prepare.json"

    @property
    def database_prepare_receipt(self) -> Path:
        return self.database_prepare_root / "prepare-receipt.json"

    @property
    def encode_runtime_prepare_root(self) -> Path:
        return self.run_root / "operator" / "encode-runtime"

    @property
    def encode_runtime_prepare_request(self) -> Path:
        return self.encode_runtime_prepare_root / "request.json"

    @property
    def encode_runtime_prepare_receipt(self) -> Path:
        return self.encode_runtime_prepare_root / "receipt.json"

    @property
    def bulk_runtime_prepare_root(self) -> Path:
        return self.run_root / "operator" / "bulk-runtime"

    @property
    def bulk_runtime_prepare_request(self) -> Path:
        return self.bulk_runtime_prepare_root / "request.json"

    @property
    def bulk_runtime_prepare_receipt(self) -> Path:
        return self.bulk_runtime_prepare_root / "receipt.json"

    def encode_runtime_failed(self, identity: str, task_identity: str) -> Path:
        return self.encode_runtime_materialized / f".failed.{identity}.{task_identity}"

    def encode_runtime_active_root(self, identity: str) -> Path:
        return self.encode_runtime_materialized / identity

    def encode_runtime_materialization_receipt(self, identity: str) -> Path:
        return (
            self.data_root
            / "operator"
            / "runtime-materializations"
            / f"{identity}.json"
        )

    @property
    def transactions(self) -> Path:
        """Compatibility name for the root-only operator journal directory."""
        return self.operator_transactions

    @property
    def service_identities(self) -> Path:
        return self.run_root / "services"

    @property
    def ingress(self) -> Path:
        return self.data_root / "operator" / "ingress"

    def ingress_bundle(self, component: str, identity: str) -> Path:
        """Return the fixed, flat ingress coordinate for one admitted bundle."""
        return self.ingress / component / f"{identity}.tar"

    @property
    def database(self) -> Path:
        return self.data_root / "database" / "live" / "platform.db"

    @property
    def database_backups(self) -> Path:
        return self.data_root / "operator" / "database-backups"

    @property
    def workspaces(self) -> Path:
        return self.data_root / "workspaces"

    @property
    def artifacts(self) -> Path:
        return self.data_root / "artifacts"

    @property
    def reference_profile_config(self) -> Path:
        return self.configuration_root / "reference-profiles.yaml"


__all__ = ["DeploymentLayout"]
