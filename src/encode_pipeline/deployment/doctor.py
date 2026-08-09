"""Stable, path-free diagnostics for the supported single-host deployment.

The doctor coordinates read-only probes.  Concrete Redis, worker, database, and
reference probes are connected only after the PR #172 runtime composition is
available; this module fixes their public contract without opening or mutating
those systems itself.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
import re
from typing import Protocol

from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.manager import DeploymentManager
from encode_pipeline.deployment.manager import DeploymentStatus
from encode_pipeline.deployment.models import COMPONENTS
from encode_pipeline.frontend_assets import load_packaged_frontend_assets


DOCTOR_REPORT_SCHEMA = "helixweave-deployment-doctor-v1"

PASS = "pass"
WARNING = "warning"
FAIL = "fail"
CHECK_STATES = (PASS, WARNING, FAIL)

CHECKS = (
    ("deployment-state", "deployment"),
    ("configuration", "host"),
    ("permissions", "host"),
    ("database", "database"),
    ("redis", "queue"),
    ("worker", "queue"),
    ("frontend", "platform"),
    ("encode-runtime", "runtime"),
    ("bulk-rnaseq-runtime", "runtime"),
    ("references", "reference"),
)

_CHECK_IDS = tuple(check_id for check_id, _category in CHECKS)
_IDENTITY = re.compile(r"^sha256-[0-9a-f]{64}$")
PUBLIC_REASON_CODES = frozenset(
    {
        "READY",
        "REFERENCES_INCOMPLETE",
        "DEPLOYMENT_INTERRUPTED",
        "DEPLOYMENT_STATE_READY",
        "DEPLOYMENT_STATE_INVALID",
        "DEPLOYMENT_STATE_UNAVAILABLE",
        "DEPLOYMENT_BUSY",
        "RUNTIME_NOT_ACTIVE",
        "RUNTIME_READY",
        "DEPLOYMENT_CONTRACT_ADMISSION_DEFERRED",
        "DEPLOYMENT_CONTRACT_ADMISSION_FAILED",
        "FRONTEND_READY",
        "FRONTEND_ASSET_MANIFEST_INVALID",
        "FRONTEND_ASSET_PACKAGE_INVALID",
        "FRONTEND_ASSET_INTEGRITY_FAILED",
        "DOCTOR_CHECK_FAILED",
        "DOCTOR_RESULT_INVALID",
    }
)


@dataclass(frozen=True)
class ProbeResult:
    """One validated probe observation without free-form public text."""

    state: str
    reason_code: str
    evidence_identity: str | None = None

    def __post_init__(self) -> None:
        if (
            self.state not in CHECK_STATES
            or self.reason_code not in PUBLIC_REASON_CODES
            or (
                self.evidence_identity is not None
                and _IDENTITY.fullmatch(self.evidence_identity) is None
            )
        ):
            raise fail("DOCTOR_RESULT_INVALID", "Doctor result is invalid.")


class DoctorProbe(Protocol):
    """A bounded read-only check returning only public evidence."""

    def __call__(self) -> ProbeResult: ...


@dataclass(frozen=True)
class DoctorCheck:
    check_id: str
    category: str
    state: str
    reason_code: str
    evidence_identity: str | None = None

    def to_dict(self) -> dict[str, object]:
        value: dict[str, object] = {
            "check_id": self.check_id,
            "category": self.category,
            "state": self.state,
            "reason_code": self.reason_code,
        }
        if self.evidence_identity is not None:
            value["evidence_identity"] = self.evidence_identity
        return value


@dataclass(frozen=True)
class DoctorReport:
    status: str
    checks: tuple[DoctorCheck, ...]

    @classmethod
    def create(cls, checks: tuple[DoctorCheck, ...]) -> "DoctorReport":
        if tuple(item.check_id for item in checks) != _CHECK_IDS:
            raise fail("DOCTOR_RESULT_INVALID", "Doctor result is invalid.")
        states = {item.state for item in checks}
        if FAIL in states:
            status = "unhealthy"
        elif WARNING in states:
            status = "degraded"
        else:
            status = "healthy"
        return cls(status=status, checks=checks)

    @property
    def ready(self) -> bool:
        return self.status == "healthy"

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": DOCTOR_REPORT_SCHEMA,
            "status": self.status,
            "ready": self.ready,
            "checks": [item.to_dict() for item in self.checks],
        }


class DeploymentDoctor:
    """Run the exact supported check inventory in a stable order."""

    def __init__(self, probes: Mapping[str, DoctorProbe]) -> None:
        if set(probes) != set(_CHECK_IDS):
            raise fail("DOCTOR_PROBES_INVALID", "Doctor probes are invalid.")
        self._probes = dict(probes)

    def run(self) -> DoctorReport:
        checks: list[DoctorCheck] = []
        for check_id, category in CHECKS:
            result = self._run_probe(self._probes[check_id])
            checks.append(
                DoctorCheck(
                    check_id=check_id,
                    category=category,
                    state=result.state,
                    reason_code=result.reason_code,
                    evidence_identity=result.evidence_identity,
                )
            )
        return DoctorReport.create(tuple(checks))

    @staticmethod
    def _run_probe(probe: DoctorProbe) -> ProbeResult:
        try:
            result = probe()
        except DeploymentError as error:
            code = error.issue.code
            if code not in PUBLIC_REASON_CODES:
                code = "DOCTOR_CHECK_FAILED"
            return ProbeResult(FAIL, code)
        except Exception:
            return ProbeResult(FAIL, "DOCTOR_CHECK_FAILED")
        if not isinstance(result, ProbeResult):
            return ProbeResult(FAIL, "DOCTOR_RESULT_INVALID")
        return result


class DeploymentStateProbe:
    """Verify immutable deployment state and report interrupted transactions."""

    def __init__(
        self,
        snapshot: "DeploymentSnapshot",
    ) -> None:
        self._snapshot = snapshot

    def __call__(self) -> ProbeResult:
        status = self._snapshot.read()
        if status.interrupted:
            return ProbeResult(
                WARNING,
                "DEPLOYMENT_INTERRUPTED",
                status.state.identity,
            )
        return ProbeResult(
            PASS,
            "DEPLOYMENT_STATE_READY",
            status.state.identity,
        )


class RuntimeProbe:
    """Verify one active immutable scientific runtime without exposing its root."""

    def __init__(
        self,
        snapshot: "DeploymentSnapshot",
        component: str,
    ) -> None:
        if component not in COMPONENTS or component == "platform":
            raise fail("DOCTOR_PROBES_INVALID", "Doctor probes are invalid.")
        self._snapshot = snapshot
        self._component = component

    def __call__(self) -> ProbeResult:
        status = self._snapshot.read()
        active = status.manifests[self._component]["active"]
        if active is None:
            return ProbeResult(FAIL, "RUNTIME_NOT_ACTIVE")
        self._snapshot.admit(self._component)
        return ProbeResult(PASS, "RUNTIME_READY", active.identity)


class DeploymentSnapshot:
    """One lightweight state/manifest metadata read shared by doctor probes.

    A snapshot is intentionally single-report scoped.  Full payload ownership,
    mode, and hashes remain the separate ``verify`` operation.
    """

    def __init__(
        self,
        manager: DeploymentManager,
    ) -> None:
        self._manager = manager
        self._value: DeploymentStatus | None = None
        self._admitted: set[str] = set()

    def read(self) -> DeploymentStatus:
        if self._value is None:
            self._value = self._manager.status()
        return self._value

    def admit(self, component: str) -> None:
        if component in self._admitted:
            return
        manifest = self.read().manifests[component]["active"]
        if manifest is None:
            raise fail("RUNTIME_NOT_ACTIVE", "Scientific runtime is not active.")
        self._manager.admit_manifest(manifest)
        self._admitted.add(component)


def frontend_probe() -> ProbeResult:
    """Verify the package-owned production frontend closure in memory."""
    assets = load_packaged_frontend_assets()
    return ProbeResult(
        PASS,
        "FRONTEND_READY",
        assets.manifest.identity,
    )


def fixed_probe(
    state: str,
    reason_code: str,
    evidence_identity: str | None = None,
) -> Callable[[], ProbeResult]:
    """Build a constant probe for adapters and deterministic contract tests."""
    result = ProbeResult(state, reason_code, evidence_identity)
    return lambda: result


__all__ = [
    "CHECKS",
    "CHECK_STATES",
    "DOCTOR_REPORT_SCHEMA",
    "FAIL",
    "PASS",
    "PUBLIC_REASON_CODES",
    "WARNING",
    "DeploymentDoctor",
    "DeploymentSnapshot",
    "DeploymentStateProbe",
    "DoctorCheck",
    "DoctorProbe",
    "DoctorReport",
    "ProbeResult",
    "RuntimeProbe",
    "fixed_probe",
    "frontend_probe",
]
