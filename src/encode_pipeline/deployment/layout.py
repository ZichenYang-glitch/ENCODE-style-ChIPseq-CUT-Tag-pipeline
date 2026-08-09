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
        return self.data_root / "operator" / "state" / "generations"

    @property
    def current_state_link(self) -> Path:
        return self.data_root / "operator" / "state" / "current"

    @property
    def state_lock(self) -> Path:
        return self.data_root / "operator" / "state" / "operation.lock"

    @property
    def transactions(self) -> Path:
        return self.data_root / "operator" / "transactions"

    @property
    def ingress(self) -> Path:
        return self.data_root / "operator" / "ingress"

    @property
    def database(self) -> Path:
        return self.data_root / "database" / "platform.db"

    @property
    def database_backups(self) -> Path:
        return self.data_root / "database" / "backups"

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
