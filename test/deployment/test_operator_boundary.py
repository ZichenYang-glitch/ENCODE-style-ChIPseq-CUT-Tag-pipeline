from __future__ import annotations

from dataclasses import dataclass
import errno
import importlib.util
import hashlib
import json
import os
from pathlib import Path
import runpy
import shutil
import sqlite3
import stat
import subprocess
import sys
from types import SimpleNamespace

import pytest

import encode_pipeline.deployment.database as database_module
import encode_pipeline.deployment.operator as operator_module
from encode_pipeline.deployment.bundle import BundleStore
from encode_pipeline.deployment.bulk_docker_boundary import BulkDockerBoundary
from encode_pipeline.deployment.canonical import canonical_json_bytes
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.database import (
    database_content_identity,
    fresh_database_candidate_path,
    inspect_database,
)
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import DeploymentState
from encode_pipeline.deployment.operator import (
    CommandResult,
    FixedObservationProvider,
    FixedSystemctl,
    HostBoundaryUninstaller,
    HostDeploymentActionController,
    HostOperatorBoundaryObserver,
    HostOperatorBackend,
    LinuxServiceProbe,
    OperatorOutcome,
    OperatorBoundaryObservation,
    OperatorObservation,
    OperatorRequest,
    SERVICE_UNITS,
    SystemdServiceController,
    SystemdBulkRuntimePreparer,
    UNINSTALL_BOUNDARY_FILES,
    UNINSTALL_LINKED_BOUNDARY_TARGETS,
    ServiceIdentity,
    SocketWitness,
    bundle_ingress_path,
    execute_request,
    parse_request,
    _candidate_schema_target,
    _plan_schema_transition,
    _reference_api_readiness,
    _resolve_account_identity,
)
from encode_pipeline.deployment.operator_boundary import (
    StableBoundaryError,
    verify_stable_operator_boundary,
)
from encode_pipeline.deployment.operator_action import (
    VERIFICATION_CHECKS,
    BulkRuntimePrepareReceipt,
    BulkRuntimePrepareRequest,
    DatabasePrepareReceipt,
    DatabasePrepareRequest,
    DeploymentActionReceipt,
    DeploymentActionRequest,
    EncodeRuntimeEntry,
    EncodeRuntimeInventory,
    EncodeRuntimePrepareReceipt,
    EncodeRuntimePrepareRequest,
    ReadinessCheck,
    SystemdEncodeRuntimePreparer,
    snakemake_environment_hash,
    verify_materialized_encode_runtime,
)
from encode_pipeline.deployment.operator_transaction import (
    OperatorJournalStore,
    OperatorTransaction,
)
from encode_pipeline.deployment.state import StateStore, render_platform_environment
from .support import manifest_for, write_bundle


ROOT = Path(__file__).resolve().parents[2]
TEMPLATES = ROOT / "src" / "encode_pipeline" / "deployment" / "templates"
BOOTSTRAP_SCRIPT = ROOT / "scripts" / "bootstrap_helixweave_operator.py"
IDENTITY = f"sha256-{'a' * 64}"
SERVICE_IDENTITY = f"sha256-{'b' * 64}"
THIRD_IDENTITY = f"sha256-{'c' * 64}"
OLD_PLATFORM_IDENTITY = f"sha256-{'d' * 64}"
OLD_ENCODE_IDENTITY = f"sha256-{'e' * 64}"
OLD_BULK_IDENTITY = f"sha256-{'f' * 64}"
TASK_IDENTITY = f"task-{'c' * 32}"


def _encode_runtime_inventory() -> EncodeRuntimeInventory:
    full_hash = "d" * 32
    link_target = "snakemake"
    entries = (
        EncodeRuntimeEntry(
            f"conda-envs/{full_hash}.env_setup_done",
            "file",
            hashlib.sha256(b"").hexdigest(),
            0,
            0o444,
            None,
        ),
        EncodeRuntimeEntry(
            f"conda-envs/{full_hash}/share/Lorem ipsum.txt",
            "file",
            hashlib.sha256(b"").hexdigest(),
            0,
            0o444,
            None,
        ),
        EncodeRuntimeEntry(
            "runner/bin/conda",
            "file",
            hashlib.sha256(b"").hexdigest(),
            0,
            0o555,
            None,
        ),
        EncodeRuntimeEntry(
            "runner/bin/snakemake",
            "file",
            hashlib.sha256(b"snakemake").hexdigest(),
            len(b"snakemake"),
            0o555,
            None,
        ),
        EncodeRuntimeEntry(
            "runner/bin/snakemake-link",
            "symlink",
            hashlib.sha256(link_target.encode()).hexdigest(),
            len(link_target),
            None,
            link_target,
        ),
    )
    return EncodeRuntimeInventory.create(entries)


def _load_bootstrap_module():
    spec = importlib.util.spec_from_file_location(
        "test_bootstrap_helixweave_operator",
        BOOTSTRAP_SCRIPT,
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


bootstrap = _load_bootstrap_module()


@dataclass
class RecordingBackend:
    outcome: OperatorOutcome
    request: OperatorRequest | None = None
    bundle_path: Path | None = None

    def execute(
        self,
        request: OperatorRequest,
        *,
        bundle_path: Path | None,
    ) -> OperatorOutcome:
        self.request = request
        self.bundle_path = bundle_path
        return self.outcome


def _worker_service_identity() -> ServiceIdentity:
    return ServiceIdentity.create(
        unit="helixweave-worker.service",
        deployment_identity=IDENTITY,
        task_identity=TASK_IDENTITY,
        main_pid=1234,
        process_start_ticks=5678,
        executable_device=42,
        executable_inode=84,
        cmdline_identity=f"sha256-{'d' * 64}",
        boot_identity=f"sha256-{'e' * 64}",
        invocation_identity=f"sha256-{'f' * 64}",
        cgroup_identity=f"sha256-{'1' * 64}",
        sockets=(),
    )


@pytest.mark.parametrize(
    ("argv", "operation", "component", "unit", "service_identity"),
    [
        (
            ("stage", "platform", IDENTITY, TASK_IDENTITY),
            "stage",
            "platform",
            None,
            None,
        ),
        (
            ("activate", "encode-runtime", IDENTITY, TASK_IDENTITY),
            "activate",
            "encode-runtime",
            None,
            None,
        ),
        (
            ("rollback", "bulk-rnaseq-runtime", IDENTITY, TASK_IDENTITY),
            "rollback",
            "bulk-rnaseq-runtime",
            None,
            None,
        ),
        (
            ("start", "helixweave-worker.service", IDENTITY, TASK_IDENTITY),
            "start",
            None,
            "helixweave-worker.service",
            None,
        ),
        (
            (
                "status",
                "helixweave-worker.service",
                IDENTITY,
                TASK_IDENTITY,
            ),
            "status",
            None,
            "helixweave-worker.service",
            None,
        ),
        (("uninstall", IDENTITY, TASK_IDENTITY), "uninstall", None, None, None),
    ],
)
def test_request_grammar_is_fixed(
    argv: tuple[str, ...],
    operation: str,
    component: str | None,
    unit: str | None,
    service_identity: str | None,
) -> None:
    request = parse_request(argv)

    assert request.operation == operation
    assert request.component == component
    assert request.unit == unit
    assert request.service_identity == service_identity
    assert request.deployment_identity == IDENTITY
    assert request.task_identity == TASK_IDENTITY


@pytest.mark.parametrize(
    "argv",
    [
        (),
        ("stage", "platform", "/tmp/private-bundle.tar", TASK_IDENTITY),
        ("stage", "platform", IDENTITY, TASK_IDENTITY, "--extra"),
        ("stage", "unknown", IDENTITY, TASK_IDENTITY),
        ("start", "/bin/sh", IDENTITY, TASK_IDENTITY),
        ("start", "helixweave-api.service", IDENTITY, "task-bad;id"),
        (
            "status",
            "helixweave-api.service",
            IDENTITY,
            TASK_IDENTITY,
            "--environment=LD_PRELOAD=/tmp/attack.so",
        ),
        ("exec", "/bin/sh", "-c", "id"),
        ("stage\nstart", "platform", IDENTITY, TASK_IDENTITY),
    ],
)
def test_request_grammar_rejects_paths_shells_environment_and_extra_args(
    argv: tuple[str, ...],
) -> None:
    with pytest.raises(DeploymentError) as caught:
        parse_request(argv)

    assert caught.value.issue.code == "OPERATOR_REQUEST_INVALID"
    assert caught.value.issue.message == "Operator request is invalid."
    assert "/tmp" not in str(caught.value)


def test_stage_derives_the_bundle_from_fixed_ingress_and_identity(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path)
    backend = RecordingBackend(OperatorOutcome("staged"))

    receipt = execute_request(
        ("stage", "platform", IDENTITY, TASK_IDENTITY),
        backend=backend,
        layout=layout,
    )

    assert backend.bundle_path == layout.ingress / "platform" / f"{IDENTITY}.tar"
    assert receipt == {
        "schema_version": "helixweave-operator-receipt-v1",
        "operation": "stage",
        "state": "staged",
        "task_identity": TASK_IDENTITY,
        "deployment_identity": IDENTITY,
        "component": "platform",
    }


def test_non_stage_requests_cannot_derive_an_ingress_path() -> None:
    request = parse_request(("activate", "platform", IDENTITY, TASK_IDENTITY))

    with pytest.raises(DeploymentError) as caught:
        bundle_ingress_path(DeploymentLayout.supported(), request)

    assert caught.value.issue.code == "OPERATOR_REQUEST_INVALID"


@dataclass
class FakeServiceController:
    service: ServiceIdentity | None = None
    stopped: tuple[str, bool] | None = None

    def start(self, request: OperatorRequest) -> ServiceIdentity:
        assert self.service is not None
        return self.service

    def status(self, request: OperatorRequest) -> ServiceIdentity | None:
        return self.service

    def stop(self, request: OperatorRequest, *, cleanup: bool) -> None:
        assert request.unit is not None
        self.stopped = (request.unit, cleanup)


def _host_backend(
    layout: DeploymentLayout,
    *,
    service_controller: FakeServiceController | None = None,
    observation_provider=None,
    verification_controller=None,
    boundary_observer=None,
    state_store: StateStore | None = None,
) -> HostOperatorBackend:
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    return HostOperatorBackend(
        layout=layout,
        service_controller=service_controller or FakeServiceController(),
        observation_provider=observation_provider,
        verification_controller=verification_controller,
        boundary_observer=boundary_observer,
        state_store=state_store,
        journal_store=OperatorJournalStore(
            layout, owner_uid=owner_uid, owner_gid=owner_gid
        ),
        root_uid=owner_uid,
        root_gid=owner_gid,
        operator_group_gid=owner_gid,
        service_uid=owner_uid + 10,
        service_gid=owner_gid + 10,
        api_uid=owner_uid + 20,
        api_gid=owner_gid + 20,
        candidate_uid=owner_uid + 30,
        candidate_gid=owner_gid + 30,
    )


def _operator_ingress(layout: DeploymentLayout, component: str) -> Path:
    path = layout.ingress / component
    path.mkdir(parents=True)
    (layout.data_root / "operator").chmod(0o710)
    layout.ingress.chmod(0o750)
    path.chmod(0o2770)
    return path


def test_host_backend_binds_distinct_service_api_and_candidate_identities(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    backend = _host_backend(layout)

    assert len({backend.service_uid, backend.api_uid, backend.candidate_uid}) == 3
    assert len({backend.service_gid, backend.api_gid, backend.candidate_gid}) == 3
    controller = backend.deployment_controller
    assert isinstance(controller, HostDeploymentActionController)
    assert (controller.service_uid, controller.service_gid) == (
        backend.service_uid,
        backend.service_gid,
    )
    assert (controller.candidate_uid, controller.candidate_gid) == (
        backend.candidate_uid,
        backend.candidate_gid,
    )


def test_host_backend_rejects_reused_runtime_account_identity(tmp_path: Path) -> None:
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    with pytest.raises(DeploymentError) as caught:
        HostOperatorBackend(
            layout=DeploymentLayout.isolated(tmp_path / "host"),
            service_controller=FakeServiceController(),
            root_uid=owner_uid,
            root_gid=owner_gid,
            operator_group_gid=owner_gid,
            service_uid=owner_uid + 10,
            service_gid=owner_gid + 10,
            api_uid=owner_uid + 10,
            api_gid=owner_gid + 20,
            candidate_uid=owner_uid + 30,
            candidate_gid=owner_gid + 30,
        )

    assert caught.value.issue.code == "OPERATOR_ACCOUNT_BOUNDARY_UNAVAILABLE"


def test_production_account_resolution_rejects_supplementary_group_drift(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    account = SimpleNamespace(pw_name="helixweave-candidate", pw_uid=4101, pw_gid=3101)
    group = SimpleNamespace(gr_name="helixweave-candidate", gr_gid=3101)
    monkeypatch.setattr(
        _resolve_account_identity.__globals__["pwd"],
        "getpwnam",
        lambda _name: account,
    )
    monkeypatch.setattr(
        _resolve_account_identity.__globals__["grp"],
        "getgrnam",
        lambda _name: group,
    )
    monkeypatch.setattr(
        _resolve_account_identity.__globals__["os"],
        "getgrouplist",
        lambda _user, primary: [primary, 9999],
    )

    with pytest.raises(DeploymentError) as caught:
        _resolve_account_identity(
            user="helixweave-candidate",
            group="helixweave-candidate",
            uid=None,
            gid=None,
        )

    assert caught.value.issue.code == "OPERATOR_ACCOUNT_BOUNDARY_UNAVAILABLE"


def test_host_backend_stages_only_the_requested_flat_manifest_identity(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manifest, payload = manifest_for("encode-runtime")
    _operator_ingress(layout, manifest.component)
    bundle = write_bundle(
        layout.ingress_bundle(manifest.component, manifest.identity),
        manifest,
        payload,
    )
    bundle.chmod(0o440)
    backend = _host_backend(layout)

    receipt = execute_request(
        ("stage", manifest.component, manifest.identity, TASK_IDENTITY),
        backend=backend,
        layout=layout,
    )

    assert receipt["state"] == "staged"
    assert (
        layout.component_store(manifest.component) / manifest.identity / "manifest.json"
    ).is_file()
    journal = layout.operator_transaction_history / f"{TASK_IDENTITY}.json"
    assert json.loads(journal.read_text())["phase"] == "complete"
    assert not layout.operator_transaction_active.exists()


def test_host_backend_rejects_legacy_ingress_mode_before_staging(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manifest, payload = manifest_for("encode-runtime")
    ingress = _operator_ingress(layout, manifest.component)
    ingress.chmod(0o2730)
    bundle = write_bundle(
        layout.ingress_bundle(manifest.component, manifest.identity),
        manifest,
        payload,
    )
    bundle.chmod(0o440)

    with pytest.raises(DeploymentError) as captured:
        _host_backend(layout).execute(
            parse_request(
                ("stage", manifest.component, manifest.identity, TASK_IDENTITY)
            ),
            bundle_path=bundle,
        )

    assert captured.value.issue.code == "OPERATOR_INGRESS_UNTRUSTED"
    assert not (layout.component_store(manifest.component) / manifest.identity).exists()


def test_host_backend_fsyncs_service_start_point_of_no_return_before_start(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    states = StateStore(layout, reader_gid=owner_gid, service_gid=owner_gid)
    with states.transaction(
        exclusive=True,
        expected_owner_uid=owner_uid,
        expected_owner_gid=owner_gid,
    ) as transaction:
        initial = transaction.initialize()
        staged = initial.stage("platform", IDENTITY)
        transaction.commit(
            staged,
            operation="stage-platform",
            expected_current_identity=initial.identity,
        )
        active = staged.activate("platform")
        transaction.commit(
            active,
            operation="activate-platform",
            expected_current_identity=staged.identity,
            platform_environment=render_platform_environment(
                layout,
                active,
                api_contract_sha256="a" * 64,
            ),
        )
    witnessed: list[dict[str, object]] = []

    class WitnessingServiceController(FakeServiceController):
        def start(self, request: OperatorRequest) -> ServiceIdentity:
            witnessed.append(json.loads(layout.operator_transaction_active.read_text()))
            return _worker_service_identity()

    backend = _host_backend(
        layout,
        service_controller=WitnessingServiceController(),
        state_store=states,
    )

    backend.execute(
        parse_request(("start", "helixweave-worker.service", IDENTITY, TASK_IDENTITY)),
        bundle_path=None,
    )

    assert witnessed[0]["phase"] == "service-starting"
    assert witnessed[0]["point_of_no_return"] is True
    assert witnessed[0]["restart_units"] == ["helixweave-worker.service"]
    assert not layout.operator_transaction_active.exists()


def test_host_backend_rejects_start_for_non_active_deployment_identity(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    states, _state = _state_with_active_components(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
    )

    class NeverStart(FakeServiceController):
        def start(self, request: OperatorRequest) -> ServiceIdentity:
            raise AssertionError(f"unexpected start: {request.unit}")

    backend = _host_backend(
        layout,
        service_controller=NeverStart(),
        state_store=states,
    )
    with pytest.raises(DeploymentError) as captured:
        backend.execute(
            parse_request(
                ("start", "helixweave-worker.service", IDENTITY, TASK_IDENTITY)
            ),
            bundle_path=None,
        )

    assert captured.value.issue.code == "OPERATOR_STATE_IDENTITY_MISMATCH"
    archived = json.loads(
        (layout.operator_transaction_history / f"{TASK_IDENTITY}.json").read_text()
    )
    assert archived["phase"] == "aborted"


@pytest.mark.parametrize("operation", ("stop", "cleanup"))
def test_bound_service_identity_is_validated_before_the_unsafe_journal_phase(
    tmp_path: Path,
    operation: str,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / operation)
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    states, _state = _state_with_active_components(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
    )
    unit = "helixweave-worker.service"
    service = _TrackingServices._identity(
        OperatorRequest(
            operation="start",
            task_identity=TASK_IDENTITY,
            deployment_identity=OLD_PLATFORM_IDENTITY,
            unit=unit,
        )
    )
    services = FakeServiceController(service=service)
    backend = _host_backend(
        layout,
        service_controller=services,
        state_store=states,
    )

    with pytest.raises(DeploymentError) as captured:
        backend.execute(
            OperatorRequest(
                operation=operation,
                task_identity=TASK_IDENTITY,
                deployment_identity=OLD_PLATFORM_IDENTITY,
                unit=unit,
                service_identity=THIRD_IDENTITY,
            ),
            bundle_path=None,
        )

    assert captured.value.issue.code == "OPERATOR_SERVICE_IDENTITY_MISMATCH"
    assert services.stopped is None
    assert not layout.operator_transaction_active.exists()
    first = json.loads(
        (layout.operator_transaction_history / f"{TASK_IDENTITY}.json").read_text()
    )
    assert first["phase"] == "aborted"
    assert first["failure_phase"] is None

    retry_task = "task-" + "9" * 32
    outcome = backend.execute(
        OperatorRequest(
            operation=operation,
            task_identity=retry_task,
            deployment_identity=OLD_PLATFORM_IDENTITY,
            unit=unit,
            service_identity=service.identity,
        ),
        bundle_path=None,
    )

    assert outcome.state == ("clean" if operation == "cleanup" else "stopped")
    assert services.stopped == (unit, operation == "cleanup")
    assert (
        json.loads(
            (layout.operator_transaction_history / f"{retry_task}.json").read_text()
        )["phase"]
        == "complete"
    )


def test_host_backend_rejects_wrong_manifest_before_staging_any_release(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manifest, payload = manifest_for("encode-runtime")
    _operator_ingress(layout, manifest.component)
    requested = IDENTITY
    assert requested != manifest.identity
    bundle = write_bundle(
        layout.ingress_bundle(manifest.component, requested), manifest, payload
    )
    bundle.chmod(0o440)

    with pytest.raises(DeploymentError) as captured:
        _host_backend(layout).execute(
            parse_request(("stage", manifest.component, requested, TASK_IDENTITY)),
            bundle_path=bundle,
        )

    assert captured.value.issue.code == "DEPLOYMENT_BUNDLE_IDENTITY_MISMATCH"
    assert not (layout.component_store(manifest.component) / manifest.identity).exists()
    journal = json.loads(
        (layout.operator_transaction_history / f"{TASK_IDENTITY}.json").read_text()
    )
    assert journal["phase"] == "aborted"
    assert not layout.operator_transaction_active.exists()


def test_host_backend_status_is_read_only_and_accepts_stopped(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    backend = _host_backend(layout, service_controller=FakeServiceController())

    outcome = backend.execute(
        parse_request(("status", "helixweave-worker.service", IDENTITY, TASK_IDENTITY)),
        bundle_path=None,
    )

    assert outcome == OperatorOutcome("stopped")
    assert not layout.transactions.exists()


@dataclass
class RecordingCommandExecutor:
    result: CommandResult
    calls: list[tuple[tuple[str, ...], float]]

    def run(self, argv: tuple[str, ...], *, timeout: float) -> CommandResult:
        self.calls.append((argv, timeout))
        return self.result


def test_systemctl_surface_builds_only_fixed_absolute_argv() -> None:
    executor = RecordingCommandExecutor(CommandResult(0), [])
    systemctl = FixedSystemctl(executor)

    systemctl.control("stop", "helixweave-worker.service")

    argv, timeout = executor.calls[0]
    assert argv == (
        "/usr/bin/systemctl",
        "--no-ask-password",
        "stop",
        "--",
        "helixweave-worker.service",
    )
    assert timeout == 15.0
    with pytest.raises(DeploymentError) as injected:
        systemctl.control("stop; /bin/sh", "helixweave-worker.service")
    assert injected.value.issue.code == "OPERATOR_COMMAND_INVALID"


def test_systemctl_observation_fails_closed_until_daemon_reload_completes() -> None:
    output = (
        b"ActiveState=inactive\n"
        b"SubState=dead\n"
        b"MainPID=0\n"
        b"InvocationID=\n"
        b"ControlGroup=\n"
        b"NeedDaemonReload=yes\n"
    )
    executor = RecordingCommandExecutor(CommandResult(0, output), [])

    with pytest.raises(DeploymentError) as captured:
        FixedSystemctl(executor).show("helixweave-worker.service")

    assert captured.value.issue.code == "OPERATOR_SYSTEMD_RELOAD_REQUIRED"
    assert "--property=" in executor.calls[0][0][4]
    assert "NeedDaemonReload" in executor.calls[0][0][4]


def test_service_start_stops_synchronously_when_identity_persistence_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    service = _worker_service_identity()

    class SequencedProbe:
        def __init__(self) -> None:
            self.values = iter((None, service, None))

        def observe(self, **_kwargs):
            return next(self.values)

    class RecordingSystemctl:
        def __init__(self) -> None:
            self.calls: list[tuple[str, str]] = []

        def control(self, action: str, unit: str) -> None:
            self.calls.append((action, unit))

    systemctl = RecordingSystemctl()
    controller = SystemdServiceController(
        DeploymentLayout.isolated(tmp_path / "host"),
        systemctl=systemctl,
        probe=SequencedProbe(),
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
    )

    def fail_write(_service: ServiceIdentity) -> None:
        raise fail(
            "OPERATOR_SERVICE_IDENTITY_UNAVAILABLE",
            "Service identity could not be persisted.",
        )

    monkeypatch.setattr(controller, "_write_identity", fail_write)
    with pytest.raises(DeploymentError) as captured:
        controller.start(
            OperatorRequest(
                operation="start",
                unit="helixweave-worker.service",
                deployment_identity=IDENTITY,
                task_identity=TASK_IDENTITY,
            )
        )

    assert captured.value.issue.code == "OPERATOR_SERVICE_IDENTITY_UNAVAILABLE"
    assert systemctl.calls == [
        ("start", "helixweave-worker.service"),
        ("stop", "helixweave-worker.service"),
    ]


def test_cleanup_recovery_resets_failed_unit_after_stop_already_completed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    service = _worker_service_identity()

    class StoppedProbe:
        def observe(self, **_kwargs):
            return None

    class RecordingSystemctl:
        def __init__(self) -> None:
            self.calls: list[tuple[str, str]] = []

        def control(self, action: str, unit: str) -> None:
            self.calls.append((action, unit))

    systemctl = RecordingSystemctl()
    controller = SystemdServiceController(
        DeploymentLayout.isolated(tmp_path / "host"),
        systemctl=systemctl,
        probe=StoppedProbe(),
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
    )
    monkeypatch.setattr(
        controller,
        "_read_identity",
        lambda _unit, *, required: service if required else None,
    )

    controller.recover_stop(
        OperatorRequest(
            operation="stop",
            unit="helixweave-worker.service",
            deployment_identity=service.deployment_identity,
            task_identity=TASK_IDENTITY,
            service_identity=service.identity,
        ),
        cleanup=True,
    )

    assert systemctl.calls == [("reset-failed", "helixweave-worker.service")]


def test_service_status_ignores_stale_identity_only_after_confirming_stopped(
    tmp_path: Path,
) -> None:
    prior = _worker_service_identity()

    class Probe:
        running = False

        def observe(
            self,
            *,
            unit: str,
            deployment_identity: str,
            task_identity: str,
        ) -> ServiceIdentity | None:
            if not self.running:
                return None
            return ServiceIdentity.create(
                unit=unit,
                deployment_identity=deployment_identity,
                task_identity=task_identity,
                main_pid=1234,
                process_start_ticks=5678,
                executable_device=42,
                executable_inode=84,
                cmdline_identity=OLD_PLATFORM_IDENTITY,
                boot_identity=OLD_ENCODE_IDENTITY,
                invocation_identity=OLD_BULK_IDENTITY,
                cgroup_identity=THIRD_IDENTITY,
                sockets=(),
            )

    probe = Probe()
    layout = DeploymentLayout.isolated(tmp_path / "host")
    layout.service_identities.mkdir(parents=True, mode=0o700)
    layout.service_identities.chmod(0o700)
    controller = SystemdServiceController(
        layout,
        systemctl=SimpleNamespace(),
        probe=probe,
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
    )
    controller._write_identity(prior)
    request = OperatorRequest(
        operation="status",
        unit=prior.unit,
        deployment_identity=OLD_PLATFORM_IDENTITY,
        task_identity=f"task-{'d' * 32}",
    )

    assert controller.status(request) is None

    probe.running = True
    with pytest.raises(DeploymentError) as captured:
        controller.status(request)
    assert captured.value.issue.code == "OPERATOR_SERVICE_IDENTITY_MISMATCH"


def test_uninstall_removes_only_fixed_boundary_and_preserves_data(
    tmp_path: Path,
) -> None:
    root = tmp_path / "host"
    for fixed in UNINSTALL_BOUNDARY_FILES:
        path = root / str(fixed).lstrip("/")
        path.parent.mkdir(parents=True, exist_ok=True)
        linked_target = UNINSTALL_LINKED_BOUNDARY_TARGETS.get(fixed)
        if linked_target is None:
            path.write_bytes(b"boundary")
            path.chmod(0o444)
        else:
            path.symlink_to(linked_target)
    secret = root / "etc/helixweave/secrets.env"
    secret.write_text("SECRET=preserve\n")
    database = root / "var/lib/helixweave/database/live/platform.db"
    database.parent.mkdir(parents=True)
    database.write_bytes(b"preserve")
    executor = RecordingCommandExecutor(CommandResult(0), [])
    uninstaller = HostBoundaryUninstaller(
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
        executor=executor,
        root_prefix=root,
    )

    callback_observations: list[tuple[bool, bool]] = []

    def before_control_removal() -> None:
        callback_observations.append(
            (
                (root / "usr/libexec/helixweave-operator").is_file(),
                (root / "etc/sudoers.d/helixweave-operator").is_file(),
            )
        )

    uninstaller.uninstall(before_control_removal=before_control_removal)

    assert callback_observations == [(True, True)]
    assert all(
        not (root / str(path).lstrip("/")).exists()
        and not (root / str(path).lstrip("/")).is_symlink()
        for path in UNINSTALL_BOUNDARY_FILES
    )
    assert secret.read_text() == "SECRET=preserve\n"
    assert database.read_bytes() == b"preserve"
    assert executor.calls == [(("/usr/bin/systemctl", "daemon-reload"), 15.0)]


def test_uninstall_link_contract_matches_bootstrap_boundary_inventory() -> None:
    assert UNINSTALL_LINKED_BOUNDARY_TARGETS == {
        spec.destination: bootstrap.HostBootstrapBackend._boundary_link_target(spec)
        for spec in bootstrap.LINKED_BOUNDARY_SPECS
    }


def test_uninstall_preflights_every_boundary_before_removing_anything(
    tmp_path: Path,
) -> None:
    root = tmp_path / "host"
    for fixed in UNINSTALL_BOUNDARY_FILES:
        path = root / str(fixed).lstrip("/")
        path.parent.mkdir(parents=True, exist_ok=True)
        linked_target = UNINSTALL_LINKED_BOUNDARY_TARGETS.get(fixed)
        if linked_target is None:
            path.write_bytes(b"boundary")
            path.chmod(0o444)
        else:
            path.symlink_to(linked_target)
    invalid = root / "etc/helixweave/redis.conf"
    invalid.unlink()
    invalid.symlink_to("/tmp/untrusted")
    executor = RecordingCommandExecutor(CommandResult(0), [])

    with pytest.raises(DeploymentError) as captured:
        HostBoundaryUninstaller(
            owner_uid=os.getuid(),
            owner_gid=os.getgid(),
            executor=executor,
            root_prefix=root,
        ).uninstall(before_control_removal=lambda: None)

    assert captured.value.issue.code == "OPERATOR_UNINSTALL_FAILED"
    assert executor.calls == []
    for fixed in UNINSTALL_BOUNDARY_FILES:
        path = root / str(fixed).lstrip("/")
        assert path.exists() or path.is_symlink()


def test_uninstall_restores_staged_boundaries_when_daemon_reload_fails(
    tmp_path: Path,
) -> None:
    root = tmp_path / "host"
    for fixed in UNINSTALL_BOUNDARY_FILES:
        path = root / str(fixed).lstrip("/")
        path.parent.mkdir(parents=True, exist_ok=True)
        linked_target = UNINSTALL_LINKED_BOUNDARY_TARGETS.get(fixed)
        if linked_target is None:
            path.write_bytes(b"boundary")
            path.chmod(0o444)
        else:
            path.symlink_to(linked_target)
    executor = RecordingCommandExecutor(CommandResult(1), [])

    with pytest.raises(DeploymentError) as captured:
        HostBoundaryUninstaller(
            owner_uid=os.getuid(),
            owner_gid=os.getgid(),
            executor=executor,
            root_prefix=root,
        ).uninstall(before_control_removal=lambda: None)

    assert captured.value.issue.code == "OPERATOR_UNINSTALL_FAILED"
    for fixed in UNINSTALL_BOUNDARY_FILES:
        path = root / str(fixed).lstrip("/")
        assert path.exists() or path.is_symlink()
        assert not path.with_name(f".{path.name}.helixweave-uninstall-pending").exists()


@pytest.mark.parametrize("repair_operation", ("install", "update"))
def test_bootstrap_repairs_a_post_commit_partial_uninstall_without_journal_deadlock(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    repair_operation: str,
) -> None:
    root = tmp_path / repair_operation
    root.mkdir()
    bootstrap_backend = bootstrap.HostBootstrapBackend(
        source_root=TEMPLATES,
        root_prefix=root,
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
        command_runner=lambda command: pytest.fail(f"unexpected command: {command}"),
        sudoers_validator=lambda _path: True,
    )
    installed = bootstrap_backend.apply(
        operation="install",
        invoking_user="labadmin",
    )
    layout = DeploymentLayout.isolated(root)
    (layout.data_root / "operator").mkdir(parents=True, mode=0o700)
    journals = OperatorJournalStore(
        layout,
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
    )
    uninstaller = HostBoundaryUninstaller(
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
        executor=RecordingCommandExecutor(CommandResult(0), []),
        root_prefix=root,
    )
    sudoers_parent = root / "etc/sudoers.d"
    original_fsync = operator_module._fsync_directory
    injected = False

    def fail_after_sudoers_unlink(path: Path) -> None:
        nonlocal injected
        if path == sudoers_parent and not injected:
            injected = True
            raise OSError("injected post-commit fsync failure")
        original_fsync(path)

    monkeypatch.setattr(operator_module, "_fsync_directory", fail_after_sudoers_unlink)
    active = {
        "platform": OLD_PLATFORM_IDENTITY,
        "encode-runtime": OLD_ENCODE_IDENTITY,
        "bulk-rnaseq-runtime": OLD_BULK_IDENTITY,
    }
    state_identity = IDENTITY

    with pytest.raises(DeploymentError) as captured:
        with journals.operation(
            operation="uninstall",
            task_identity=TASK_IDENTITY,
            deployment_identity=state_identity,
            component=None,
            unit=None,
        ) as journal:
            journal.advance(
                "writers-stopped",
                write_fence=True,
                prior_active=active,
                candidate_active=active,
                evidence={
                    "prior_state_identity": state_identity,
                    "candidate_state_identity": state_identity,
                },
            )
            uninstaller.uninstall(before_control_removal=journal.complete)

    assert captured.value.issue.code == "OPERATOR_UNINSTALL_FAILED"
    assert injected
    assert not layout.operator_transaction_active.exists()
    assert (
        json.loads(
            (layout.operator_transaction_history / f"{TASK_IDENTITY}.json").read_text()
        )["phase"]
        == "complete"
    )
    assert not (root / "etc/sudoers.d/helixweave-operator").exists()
    assert (root / "usr/libexec/helixweave-operator").is_file()

    repaired = bootstrap_backend.apply(
        operation=repair_operation,
        invoking_user="labadmin",
    )

    assert repaired.closure_identity == installed.closure_identity
    assert bootstrap_backend.verify_selected_boundaries() == installed.closure_identity
    assert (root / "etc/sudoers.d/helixweave-operator").is_file()
    assert (root / "usr/libexec/helixweave-operator").is_file()


@pytest.mark.parametrize(
    "error",
    [
        RuntimeError("/private/reference secret=value"),
        fail("BACKEND_PRIVATE", "/private/reference secret=value"),
    ],
)
def test_backend_errors_are_replaced_by_one_redacted_failure(error: Exception) -> None:
    class FailingBackend:
        def execute(self, request, *, bundle_path):
            del request, bundle_path
            raise error

    with pytest.raises(DeploymentError) as caught:
        execute_request(
            ("stage", "platform", IDENTITY, TASK_IDENTITY),
            backend=FailingBackend(),
        )

    assert caught.value.issue.code == "OPERATOR_OPERATION_FAILED"
    assert caught.value.issue.message == "Operator action failed."
    assert "private" not in str(caught.value)
    assert "secret" not in str(caught.value)


def test_start_receipt_binds_task_process_executable_cmdline_and_socket_contract() -> (
    None
):
    service = _worker_service_identity()
    backend = RecordingBackend(OperatorOutcome("running", service=service))

    receipt = execute_request(
        ("start", "helixweave-worker.service", IDENTITY, TASK_IDENTITY),
        backend=backend,
    )

    assert receipt["service_identity"] == service.to_dict()
    assert service.identity.startswith("sha256-")
    assert service.main_pid == 1234
    assert service.process_start_ticks == 5678
    assert service.executable_device == 42
    assert service.executable_inode == 84
    assert service.cmdline_identity == f"sha256-{'d' * 64}"
    assert service.invocation_identity == f"sha256-{'f' * 64}"
    assert service.cgroup_identity == f"sha256-{'1' * 64}"


def test_status_observes_the_current_service_identity_without_a_prior_identity() -> (
    None
):
    service = _worker_service_identity()
    backend = RecordingBackend(OperatorOutcome("running", service=service))

    receipt = execute_request(
        (
            "status",
            "helixweave-worker.service",
            IDENTITY,
            TASK_IDENTITY,
        ),
        backend=backend,
    )

    assert receipt["state"] == "running"
    assert receipt["service_identity"] == service.to_dict()


def test_status_can_report_stopped_without_a_service_identity() -> None:
    receipt = execute_request(
        ("status", "helixweave-worker.service", IDENTITY, TASK_IDENTITY),
        backend=RecordingBackend(OperatorOutcome("stopped")),
    )

    assert receipt["state"] == "stopped"
    assert "service_identity" not in receipt


@pytest.mark.parametrize(
    ("operation", "state"), (("stop", "stopped"), ("cleanup", "clean"))
)
def test_bound_service_receipt_retains_the_prior_identity(
    operation: str,
    state: str,
) -> None:
    receipt = execute_request(
        (
            operation,
            "helixweave-worker.service",
            IDENTITY,
            TASK_IDENTITY,
            SERVICE_IDENTITY,
        ),
        backend=RecordingBackend(OperatorOutcome(state)),
    )

    assert receipt["prior_service_identity"] == SERVICE_IDENTITY


def test_service_socket_contract_is_unit_specific() -> None:
    api_socket = SocketWitness(name="api-http", device=1, inode=2, kernel_inode=3)

    with pytest.raises(DeploymentError) as caught:
        ServiceIdentity.create(
            unit="helixweave-worker.service",
            deployment_identity=IDENTITY,
            task_identity=TASK_IDENTITY,
            main_pid=1,
            process_start_ticks=1,
            executable_device=1,
            executable_inode=1,
            cmdline_identity=IDENTITY,
            boot_identity=IDENTITY,
            invocation_identity=IDENTITY,
            cgroup_identity=IDENTITY,
            sockets=(api_socket,),
        )

    assert caught.value.issue.code == "OPERATOR_SERVICE_IDENTITY_INVALID"


def _unix_socket_probe(
    tmp_path: Path,
    *,
    kernel_inode: int = 4567,
    owned_inode: int | None = None,
    flags: str = "00010000",
) -> tuple[LinuxServiceProbe, Path]:
    proc_root = tmp_path / "proc"
    cgroup_root = tmp_path / "cgroup"
    fd_root = proc_root / "123/fd"
    fd_root.mkdir(parents=True)
    group = cgroup_root / "system.slice/helixweave-redis.service"
    group.mkdir(parents=True)
    (group / "cgroup.procs").write_text("123\n", encoding="ascii")
    (fd_root / "3").symlink_to(f"socket:[{owned_inode or kernel_inode}]")
    socket_path = tmp_path / "run/redis.sock"
    socket_path.parent.mkdir(parents=True)
    filesystem_witness = SimpleNamespace(
        st_mode=stat.S_IFSOCK | 0o660,
        st_nlink=1,
        st_dev=41,
        st_ino=84,
        st_uid=os.getuid(),
        st_gid=os.getgid(),
        st_size=0,
        st_mtime_ns=1,
        st_ctime_ns=2,
    )
    (proc_root / "net").mkdir()
    (proc_root / "net/unix").write_text(
        "Num RefCount Protocol Flags Type St Inode Path\n"
        f"00000000: 00000002 00000000 {flags} 0001 01 "
        f"{kernel_inode} {socket_path}\n",
        encoding="utf-8",
    )
    probe = LinuxServiceProbe(
        SimpleNamespace(),
        proc_root=proc_root,
        cgroup_root=cgroup_root,
        unix_sockets={
            "helixweave-redis.service": socket_path,
            "helixweave-docker-rootless.service": tmp_path / "run/docker.sock",
        },
        filesystem_socket_stat=lambda _path: filesystem_witness,
    )
    return probe, socket_path


def test_linux_service_probe_binds_filesystem_and_kernel_socket_inodes(
    tmp_path: Path,
) -> None:
    probe, _socket_path = _unix_socket_probe(tmp_path)
    witnesses = probe._socket_witnesses(
        unit="helixweave-redis.service",
        cgroup="/system.slice/helixweave-redis.service",
        main_pid=123,
    )

    assert witnesses == (
        SocketWitness(
            "redis-queue",
            41,
            84,
            4567,
        ),
    )


def test_linux_service_probe_rejects_kernel_socket_owned_by_another_cgroup(
    tmp_path: Path,
) -> None:
    probe, _socket_path = _unix_socket_probe(tmp_path, owned_inode=9999)
    with pytest.raises(DeploymentError) as caught:
        probe._socket_witnesses(
            unit="helixweave-redis.service",
            cgroup="/system.slice/helixweave-redis.service",
            main_pid=123,
        )

    assert caught.value.issue.code == "OPERATOR_SERVICE_OBSERVE_FAILED"


def test_linux_service_probe_rejects_stale_non_listening_proc_socket(
    tmp_path: Path,
) -> None:
    probe, _socket_path = _unix_socket_probe(tmp_path, flags="00000000")
    with pytest.raises(DeploymentError) as caught:
        probe._socket_witnesses(
            unit="helixweave-redis.service",
            cgroup="/system.slice/helixweave-redis.service",
            main_pid=123,
        )

    assert caught.value.issue.code == "OPERATOR_SERVICE_OBSERVE_FAILED"


def test_linux_service_probe_rejects_filesystem_socket_swap(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    probe, _socket_path = _unix_socket_probe(tmp_path)
    observations = iter(
        (
            probe.filesystem_socket_stat(_socket_path),
            SimpleNamespace(
                st_mode=stat.S_IFSOCK | 0o660,
                st_nlink=1,
                st_dev=41,
                st_ino=85,
                st_uid=os.getuid(),
                st_gid=os.getgid(),
                st_size=0,
                st_mtime_ns=3,
                st_ctime_ns=4,
            ),
        )
    )
    monkeypatch.setattr(
        probe, "filesystem_socket_stat", lambda _path: next(observations)
    )
    with pytest.raises(DeploymentError) as caught:
        probe._socket_witnesses(
            unit="helixweave-redis.service",
            cgroup="/system.slice/helixweave-redis.service",
            main_pid=123,
        )

    assert caught.value.issue.code == "OPERATOR_SERVICE_OBSERVE_FAILED"


@pytest.mark.parametrize("bounded_source", ("fds", "proc-net-unix"))
def test_linux_service_probe_applies_socket_scan_bounds_incrementally(
    tmp_path: Path,
    bounded_source: str,
) -> None:
    probe, _socket_path = _unix_socket_probe(tmp_path, owned_inode=9999)
    if bounded_source == "fds":
        probe._MAX_PROC_FDS = 1
        (probe.proc_root / "123/fd/4").symlink_to("socket:[8888]")
    else:
        probe._MAX_PROC_NET_UNIX_BYTES = 32
    with pytest.raises(DeploymentError) as caught:
        probe._socket_witnesses(
            unit="helixweave-redis.service",
            cgroup="/system.slice/helixweave-redis.service",
            main_pid=123,
        )

    assert caught.value.issue.code == "OPERATOR_SERVICE_OBSERVE_FAILED"


def test_only_supported_units_are_exposed() -> None:
    assert SERVICE_UNITS == (
        "helixweave-redis.service",
        "helixweave-docker-rootless.service",
        "helixweave-api.service",
        "helixweave-worker.service",
    )
    assert "docker.service" not in SERVICE_UNITS


def test_encode_runtime_inventory_is_bounded_canonical_and_path_free() -> None:
    inventory = _encode_runtime_inventory()
    request = EncodeRuntimePrepareRequest.create(
        task_identity=TASK_IDENTITY,
        deployment_identity=IDENTITY,
        authority_platform_identity=SERVICE_IDENTITY,
        prior_state_identity=IDENTITY,
        candidate_state_identity=SERVICE_IDENTITY,
    )
    receipt = EncodeRuntimePrepareReceipt.create(
        request_identity=request.identity,
        deployment_identity=IDENTITY,
        inventory=inventory,
    )

    assert EncodeRuntimeInventory.from_dict(inventory.to_dict()) == inventory
    assert EncodeRuntimePrepareReceipt.from_dict(receipt.to_dict()) == receipt
    assert receipt.entry_count == len(inventory.entries)
    assert receipt.inventory_size == len(inventory.to_bytes())
    assert "entries" not in receipt.to_dict()
    assert b"/var/" not in inventory.to_bytes()
    assert any(" " in item.path for item in inventory.entries)


def test_encode_runtime_inventory_rejects_a_symlink_escape() -> None:
    raw = _encode_runtime_inventory().to_dict()
    escaped = dict(raw["entries"][-1])
    target = "../../../../etc/passwd"
    escaped.update(
        target=target,
        size=len(target),
        sha256=hashlib.sha256(target.encode()).hexdigest(),
    )
    raw["entries"][-1] = escaped

    with pytest.raises(DeploymentError) as caught:
        EncodeRuntimeInventory.from_dict(raw)

    assert caught.value.issue.code == "ENCODE_RUNTIME_PREPARE_RECEIPT_INVALID"


def test_encode_runtime_root_is_frozen_in_place_with_safe_relative_symlinks(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    inventory = _encode_runtime_inventory()
    root = layout.encode_runtime_active_root(IDENTITY)
    root.mkdir(parents=True)
    for entry in inventory.entries:
        path = root / entry.path
        path.parent.mkdir(parents=True, exist_ok=True)
        if entry.kind == "symlink":
            path.symlink_to(entry.target)
        else:
            path.write_bytes(b"snakemake" if entry.path.endswith("/snakemake") else b"")
            path.chmod(0o600)
    inventory_path = root / ".helixweave-runtime-inventory.json"
    inventory_path.write_bytes(inventory.to_bytes())
    inventory_path.chmod(0o600)
    for path in sorted(root.rglob("*"), key=lambda item: len(item.parts), reverse=True):
        if path.is_dir() and not path.is_symlink():
            path.chmod(0o700)
    root.chmod(0o700)
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    preparer = SystemdEncodeRuntimePreparer(
        layout,
        service_uid=owner_uid,
        service_gid=owner_gid,
        root_uid=owner_uid,
        root_gid=owner_gid,
        command_runner=lambda _argv: pytest.fail("systemd must not run"),
    )

    preparer._verify_and_freeze(root, inventory)
    preparer._verify_frozen_tree(root, inventory)

    assert stat.S_IMODE(root.stat().st_mode) == 0o555
    assert stat.S_IMODE((root / "runner/bin/snakemake").stat().st_mode) == 0o555
    assert stat.S_IMODE(inventory_path.stat().st_mode) == 0o444
    assert os.readlink(root / "runner/bin/snakemake-link") == "snakemake"


def test_materialized_encode_verifier_uses_final_prefix_hash_and_real_markers(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    environment_content = b"name: tool\ndependencies:\n  - samtools=1.20\n"
    root = layout.encode_runtime_active_root(IDENTITY)
    full_hash = snakemake_environment_hash(
        root / "conda-envs",
        environment_content,
    )
    entries = EncodeRuntimeInventory.create(
        (
            EncodeRuntimeEntry(
                f"conda-envs/{full_hash}.env_setup_done",
                "file",
                hashlib.sha256(b"").hexdigest(),
                0,
                0o444,
                None,
            ),
            EncodeRuntimeEntry(
                f"conda-envs/{full_hash}/bin/samtools",
                "file",
                hashlib.sha256(b"samtools").hexdigest(),
                len(b"samtools"),
                0o555,
                None,
            ),
            EncodeRuntimeEntry(
                "mamba-root/bin/activate",
                "file",
                hashlib.sha256(b"activate").hexdigest(),
                len(b"activate"),
                0o555,
                None,
            ),
            EncodeRuntimeEntry(
                "runner/bin/conda",
                "file",
                hashlib.sha256(b"conda").hexdigest(),
                len(b"conda"),
                0o555,
                None,
            ),
            EncodeRuntimeEntry(
                "runner/bin/snakemake",
                "file",
                hashlib.sha256(b"snakemake").hexdigest(),
                len(b"snakemake"),
                0o555,
                None,
            ),
            EncodeRuntimeEntry(
                "runner/libexec/micromamba",
                "file",
                hashlib.sha256(b"micromamba").hexdigest(),
                len(b"micromamba"),
                0o555,
                None,
            ),
        )
    )
    root.mkdir(parents=True)
    content_by_path = {
        f"conda-envs/{full_hash}.env_setup_done": b"",
        f"conda-envs/{full_hash}/bin/samtools": b"samtools",
        "mamba-root/bin/activate": b"activate",
        "runner/bin/conda": b"conda",
        "runner/bin/snakemake": b"snakemake",
        "runner/libexec/micromamba": b"micromamba",
    }
    for entry in entries.entries:
        path = root / entry.path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content_by_path[entry.path])
        path.chmod(0o700 if entry.mode == 0o555 else 0o600)
    inventory_path = root / ".helixweave-runtime-inventory.json"
    inventory_path.write_bytes(entries.to_bytes())
    inventory_path.chmod(0o600)
    for path in sorted(root.rglob("*"), key=lambda item: len(item.parts), reverse=True):
        if path.is_dir() and not path.is_symlink():
            path.chmod(0o700)
    root.chmod(0o700)
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    preparer = SystemdEncodeRuntimePreparer(
        layout,
        service_uid=owner_uid,
        service_gid=owner_gid,
        root_uid=owner_uid,
        root_gid=owner_gid,
        command_runner=lambda _argv: pytest.fail("systemd must not run"),
    )
    preparer._verify_and_freeze(root, entries)

    source = layout.encode_runtimes / IDENTITY / "payload/contracts/encode-runtime"
    environment_path = source / "workflow/envs/tool.yml"
    environment_path.parent.mkdir(parents=True)
    environment_path.write_bytes(environment_content)
    environment_path.chmod(0o444)
    index_path = source / "package-index.json"
    index_path.write_bytes(b"{}")
    index_path.chmod(0o444)
    import encode_pipeline.deployment.native_contracts as native_contracts

    monkeypatch.setattr(
        native_contracts,
        "parse_encode_runtime_index",
        lambda _content: SimpleNamespace(
            environments=(
                {
                    "environment_path": "workflow/envs/runner.yml",
                },
                {
                    "environment_path": "workflow/envs/tool.yml",
                },
            )
        ),
    )

    observed = verify_materialized_encode_runtime(
        layout,
        IDENTITY,
        root_uid=owner_uid,
        root_gid=owner_gid,
    )

    assert observed.deployment_identity == IDENTITY
    assert observed.tree_identity == entries.tree_identity


def test_encode_runtime_failed_build_is_quarantined_without_deletion(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    destination = layout.encode_runtime_active_root(IDENTITY)
    destination.mkdir(parents=True)
    evidence = destination / "evidence.txt"
    evidence.write_text("preserve")
    request = EncodeRuntimePrepareRequest.create(
        task_identity=TASK_IDENTITY,
        deployment_identity=IDENTITY,
        authority_platform_identity=SERVICE_IDENTITY,
        prior_state_identity=IDENTITY,
        candidate_state_identity=SERVICE_IDENTITY,
    )
    preparer = SystemdEncodeRuntimePreparer(
        layout,
        service_uid=os.getuid(),
        service_gid=os.getgid(),
        root_uid=os.getuid(),
        root_gid=os.getgid(),
        command_runner=lambda _argv: pytest.fail("systemd must not run"),
    )

    preparer._quarantine_failed(destination, request)

    failed = layout.encode_runtime_failed(IDENTITY, TASK_IDENTITY)
    assert not destination.exists()
    assert (failed / "evidence.txt").read_text() == "preserve"
    assert stat.S_IMODE(failed.stat().st_mode) == 0o500


def test_systemd_bulk_preparer_uses_fixed_canonical_exchange_and_unit(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    root = layout.bulk_runtime_prepare_root
    root.mkdir(parents=True)
    root.chmod(0o750)
    request = BulkRuntimePrepareRequest.create(
        operation="verify",
        task_identity=TASK_IDENTITY,
        candidate_bulk_identity=THIRD_IDENTITY,
        authority_platform_identity=IDENTITY,
        prior_state_identity=SERVICE_IDENTITY,
        candidate_state_identity=SERVICE_IDENTITY,
        docker_service_identity=IDENTITY,
        docker_client_identity=SERVICE_IDENTITY,
        docker_endpoint_identity=THIRD_IDENTITY,
        docker_daemon_uid=owner_uid,
        docker_daemon_gid=owner_gid,
    )
    expected_receipt = BulkRuntimePrepareReceipt.create(
        request_identity=request.identity,
        candidate_bulk_identity=THIRD_IDENTITY,
        runtime_identity=IDENTITY,
        image_set_identity=SERVICE_IDENTITY,
        image_count=2,
    )
    commands: list[tuple[str, ...]] = []

    def run(argv: tuple[str, ...]) -> int:
        commands.append(argv)
        layout.bulk_runtime_prepare_receipt.write_bytes(
            (
                json.dumps(
                    expected_receipt.to_dict(), sort_keys=True, separators=(",", ":")
                )
                + "\n"
            ).encode()
        )
        layout.bulk_runtime_prepare_receipt.chmod(0o600)
        return 0

    preparer = SystemdBulkRuntimePreparer(
        layout,
        service_uid=owner_uid,
        service_gid=owner_gid,
        root_uid=owner_uid,
        root_gid=owner_gid,
        command_runner=run,
        boundary_observer=lambda _uid, _gid: BulkDockerBoundary(
            SERVICE_IDENTITY,
            THIRD_IDENTITY,
            owner_uid,
            owner_gid,
        ),
    )

    observed = preparer.prepare(request)

    assert observed == expected_receipt
    assert commands == [
        (
            "/usr/bin/systemctl",
            "--no-ask-password",
            "start",
            "--",
            "helixweave-bulk-runtime-prepare.service",
        )
    ]
    assert (
        layout.bulk_runtime_prepare_request.read_bytes()
        == (
            json.dumps(request.to_dict(), sort_keys=True, separators=(",", ":")) + "\n"
        ).encode()
    )
    assert stat.S_IMODE(layout.bulk_runtime_prepare_request.stat().st_mode) == 0o640
    assert stat.S_IMODE(layout.bulk_runtime_prepare_receipt.stat().st_mode) == 0o600
    assert stat.S_IMODE((root / "operation.lock").stat().st_mode) == 0o600


@pytest.mark.parametrize(
    ("operation", "database_mode", "backup_receipt_identity"),
    (
        ("activate", "fresh-candidate", None),
        ("activate", "existing-live", OLD_PLATFORM_IDENTITY),
        ("rollback", "existing-live", OLD_PLATFORM_IDENTITY),
    ),
)
def test_database_prepare_dispatcher_accepts_canonical_request(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    operation: str,
    database_mode: str,
    backup_receipt_identity: str | None,
) -> None:
    request = DatabasePrepareRequest.create(
        operation=operation,
        database_mode=database_mode,
        task_identity=TASK_IDENTITY,
        deployment_identity=IDENTITY,
        prior_state_identity=SERVICE_IDENTITY,
        candidate_state_identity=THIRD_IDENTITY,
        action_receipt_identity=IDENTITY,
        backup_receipt_identity=backup_receipt_identity,
        target_schema_heads=("schema-v1",),
    )
    request_path = tmp_path / "prepare.json"
    request_path.write_bytes(canonical_json_bytes(request.to_dict()))
    request_path.chmod(0o640)
    namespace = runpy.run_path(str(TEMPLATES / "helixweave-db-prepare"))
    read_request = namespace["_read_request"]
    read_request.__globals__["REQUEST"] = request_path
    original_fstat = os.fstat

    def root_owned_fstat(descriptor: int):
        observed = original_fstat(descriptor)
        return SimpleNamespace(
            st_mode=observed.st_mode,
            st_nlink=observed.st_nlink,
            st_uid=0,
            st_gid=os.getegid(),
            st_size=observed.st_size,
            st_dev=observed.st_dev,
            st_ino=observed.st_ino,
            st_mtime_ns=observed.st_mtime_ns,
            st_ctime_ns=observed.st_ctime_ns,
        )

    with monkeypatch.context() as scoped:
        scoped.setattr(read_request.__globals__["os"], "fstat", root_owned_fstat)
        observed = read_request()

    assert observed == request.to_dict()


@pytest.mark.parametrize("encoding", ("missing-newline", "extra-newline", "pretty"))
def test_database_prepare_dispatcher_rejects_noncanonical_request(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    encoding: str,
) -> None:
    request = DatabasePrepareRequest.create(
        operation="rollback",
        database_mode="existing-live",
        task_identity=TASK_IDENTITY,
        deployment_identity=IDENTITY,
        prior_state_identity=SERVICE_IDENTITY,
        candidate_state_identity=THIRD_IDENTITY,
        action_receipt_identity=IDENTITY,
        backup_receipt_identity=OLD_PLATFORM_IDENTITY,
        target_schema_heads=("schema-v1",),
    )
    if encoding == "missing-newline":
        content = canonical_json_bytes(request.to_dict())[:-1]
    elif encoding == "extra-newline":
        content = canonical_json_bytes(request.to_dict()) + b"\n"
    else:
        content = (
            json.dumps(request.to_dict(), indent=2, sort_keys=True) + "\n"
        ).encode()
    request_path = tmp_path / "prepare.json"
    request_path.write_bytes(content)
    request_path.chmod(0o640)
    namespace = runpy.run_path(str(TEMPLATES / "helixweave-db-prepare"))
    read_request = namespace["_read_request"]
    read_request.__globals__["REQUEST"] = request_path
    original_fstat = os.fstat

    def root_owned_fstat(descriptor: int):
        observed = original_fstat(descriptor)
        return SimpleNamespace(
            st_mode=observed.st_mode,
            st_nlink=observed.st_nlink,
            st_uid=0,
            st_gid=os.getegid(),
            st_size=observed.st_size,
            st_dev=observed.st_dev,
            st_ino=observed.st_ino,
            st_mtime_ns=observed.st_mtime_ns,
            st_ctime_ns=observed.st_ctime_ns,
        )

    with monkeypatch.context() as scoped:
        scoped.setattr(read_request.__globals__["os"], "fstat", root_owned_fstat)
        with pytest.raises(ValueError):
            read_request()


@pytest.mark.parametrize(
    ("encoding", "expected_exit"),
    (
        ("canonical", 0),
        ("missing-newline", 65),
        ("extra-newline", 65),
        ("pretty", 65),
    ),
)
def test_database_prepare_dispatcher_validates_canonical_child_receipt(
    monkeypatch: pytest.MonkeyPatch,
    encoding: str,
    expected_exit: int,
) -> None:
    request = DatabasePrepareRequest.create(
        operation="activate",
        database_mode="existing-live",
        task_identity=TASK_IDENTITY,
        deployment_identity=IDENTITY,
        prior_state_identity=SERVICE_IDENTITY,
        candidate_state_identity=THIRD_IDENTITY,
        action_receipt_identity=IDENTITY,
        backup_receipt_identity=OLD_PLATFORM_IDENTITY,
        target_schema_heads=("schema-v1",),
    )
    receipt = DatabasePrepareReceipt.create(
        request_identity=request.identity,
        database_before_identity=SERVICE_IDENTITY,
        database_after_identity=THIRD_IDENTITY,
        schema_heads=("schema-v1",),
    )
    if encoding == "canonical":
        content = canonical_json_bytes(receipt.to_dict())
    elif encoding == "missing-newline":
        content = canonical_json_bytes(receipt.to_dict())[:-1]
    elif encoding == "extra-newline":
        content = canonical_json_bytes(receipt.to_dict()) + b"\n"
    else:
        content = (
            json.dumps(receipt.to_dict(), indent=2, sort_keys=True) + "\n"
        ).encode()
    namespace = runpy.run_path(str(TEMPLATES / "helixweave-db-prepare"))
    main = namespace["main"]
    written: list[bytes] = []

    def run(_argv, *, stdout, **_kwargs):
        stdout.write(content)
        return SimpleNamespace(returncode=0)

    with monkeypatch.context() as scoped:
        scoped.setattr(main.__globals__["os"], "environ", {})
        scoped.setattr(main.__globals__["subprocess"], "run", run)
        scoped.setitem(main.__globals__, "_read_request", request.to_dict)
        scoped.setitem(main.__globals__, "_launcher", lambda _identity: Path("/fixed"))
        scoped.setitem(main.__globals__, "_write_receipt", written.append)
        observed_exit = main([])

    assert observed_exit == expected_exit
    assert written == ([content] if expected_exit == 0 else [])


def test_database_prepare_dispatcher_canonical_renderer_matches_authority() -> None:
    namespace = runpy.run_path(str(TEMPLATES / "helixweave-db-prepare"))
    render = namespace["_canonical_json_bytes"]
    value = {"unicode": "核", "nested": {"b": 2, "a": 1}}

    assert render(value) == canonical_json_bytes(value)
    with pytest.raises(ValueError):
        render({"invalid": float("nan")})


def test_database_prepare_dispatcher_loads_in_isolated_stdlib_python() -> None:
    program = (
        "import json,runpy,sys\n"
        "namespace=runpy.run_path(sys.argv[1])\n"
        "value={'unicode':'核','nested':{'b':2,'a':1}}\n"
        "expected=(json.dumps(value,allow_nan=False,ensure_ascii=False,"
        "separators=(',',':'),sort_keys=True)+'\\n').encode('utf-8')\n"
        "raise SystemExit(0 if namespace['_canonical_json_bytes'](value)==expected else 1)\n"
    )

    completed = subprocess.run(
        (
            "/usr/bin/python3",
            "-I",
            "-S",
            "-c",
            program,
            str(TEMPLATES / "helixweave-db-prepare"),
        ),
        stdin=subprocess.DEVNULL,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        cwd="/",
        env={"LANG": "C.UTF-8", "LC_ALL": "C.UTF-8", "PATH": "/usr/bin"},
        timeout=10,
        check=False,
    )

    assert completed.returncode == 0
    assert completed.stdout == b""
    assert completed.stderr == b""


_CANONICAL_LAUNCHERS = (
    "helixweave-operator-action",
    "helixweave-encode-runtime-prepare",
)


def _canonical_launcher_documents(
    launcher: str,
) -> tuple[dict[str, object], dict[str, object]]:
    if launcher == "helixweave-operator-action":
        request = DeploymentActionRequest.create(
            phase="admit",
            operation="activate",
            component="platform",
            task_identity=TASK_IDENTITY,
            deployment_identity=IDENTITY,
            authority_platform_identity=IDENTITY,
            prior_state_identity=SERVICE_IDENTITY,
            candidate_state_identity=THIRD_IDENTITY,
            candidate_active={
                "platform": IDENTITY,
                "encode-runtime": SERVICE_IDENTITY,
                "bulk-rnaseq-runtime": THIRD_IDENTITY,
            },
        )
        receipt = DeploymentActionReceipt.create(
            request_identity=request.identity,
            status="admitted",
            compatibility="compatible",
            database_before_identity=IDENTITY,
            accepted_schema_heads=("head",),
            target_schema_heads=("head",),
            migration_inventory_identity=SERVICE_IDENTITY,
            known_schema_revisions=("ancestor",),
            migration_required=False,
            rollback_supported=True,
            api_contract_sha256="a" * 64,
            native_identities={
                "platform": IDENTITY,
                "encode-runtime": SERVICE_IDENTITY,
                "bulk-rnaseq-runtime": THIRD_IDENTITY,
            },
            frontend_identity=IDENTITY,
            reference_compatibility_identity=SERVICE_IDENTITY,
            readiness={
                check: ReadinessCheck("ready", "READY", IDENTITY)
                for check in VERIFICATION_CHECKS
            },
        )
        return request.to_dict(), receipt.to_dict()
    if launcher == "helixweave-encode-runtime-prepare":
        request = EncodeRuntimePrepareRequest.create(
            task_identity=TASK_IDENTITY,
            deployment_identity=IDENTITY,
            authority_platform_identity=SERVICE_IDENTITY,
            prior_state_identity=IDENTITY,
            candidate_state_identity=SERVICE_IDENTITY,
        )
        receipt = EncodeRuntimePrepareReceipt.create(
            request_identity=request.identity,
            deployment_identity=IDENTITY,
            inventory=_encode_runtime_inventory(),
        )
        return request.to_dict(), receipt.to_dict()
    raise AssertionError(launcher)


def _document_encoding(value: dict[str, object], encoding: str) -> bytes:
    canonical = canonical_json_bytes(value)
    if encoding == "canonical":
        return canonical
    if encoding == "missing-newline":
        return canonical[:-1]
    if encoding == "extra-newline":
        return canonical + b"\n"
    if encoding == "pretty":
        return (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()
    raise AssertionError(encoding)


@pytest.mark.parametrize("launcher", _CANONICAL_LAUNCHERS)
@pytest.mark.parametrize(
    "encoding", ("canonical", "missing-newline", "extra-newline", "pretty")
)
def test_candidate_launcher_validates_canonical_request_bytes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    launcher: str,
    encoding: str,
) -> None:
    request, _receipt = _canonical_launcher_documents(launcher)
    path = tmp_path / "request.json"
    path.write_bytes(_document_encoding(request, encoding))
    path.chmod(0o640)
    namespace = runpy.run_path(str(TEMPLATES / launcher))
    read_request = namespace["_request"]
    read_request.__globals__["REQUEST"] = path
    original_fstat = os.fstat

    def root_owned_fstat(descriptor: int):
        observed = original_fstat(descriptor)
        return SimpleNamespace(
            st_mode=observed.st_mode,
            st_nlink=observed.st_nlink,
            st_uid=0,
            st_gid=os.getegid(),
            st_size=observed.st_size,
            st_dev=observed.st_dev,
            st_ino=observed.st_ino,
            st_mtime_ns=observed.st_mtime_ns,
            st_ctime_ns=observed.st_ctime_ns,
        )

    with monkeypatch.context() as scoped:
        scoped.setattr(read_request.__globals__["os"], "fstat", root_owned_fstat)
        if encoding == "canonical":
            assert read_request() == request
        else:
            with pytest.raises(ValueError):
                read_request()


@pytest.mark.parametrize("launcher", _CANONICAL_LAUNCHERS)
@pytest.mark.parametrize(
    ("encoding", "expected_exit"),
    (
        ("canonical", 0),
        ("missing-newline", 65),
        ("extra-newline", 65),
        ("pretty", 65),
    ),
)
def test_candidate_launcher_validates_canonical_child_receipt_bytes(
    monkeypatch: pytest.MonkeyPatch,
    launcher: str,
    encoding: str,
    expected_exit: int,
) -> None:
    request, receipt = _canonical_launcher_documents(launcher)
    content = _document_encoding(receipt, encoding)
    namespace = runpy.run_path(str(TEMPLATES / launcher))
    main = namespace["main"]
    written: list[bytes] = []

    def run(_argv, *, stdout, **_kwargs):
        stdout.write(content)
        return SimpleNamespace(returncode=0)

    with monkeypatch.context() as scoped:
        scoped.setattr(main.__globals__["os"], "environ", {})
        scoped.setattr(main.__globals__["subprocess"], "run", run)
        scoped.setitem(main.__globals__, "_request", lambda: request)
        scoped.setitem(main.__globals__, "_launcher", lambda _identity: Path("/fixed"))
        scoped.setitem(main.__globals__, "_write_receipt", written.append)
        observed_exit = main([])

    assert observed_exit == expected_exit
    assert written == ([content] if expected_exit == 0 else [])


@pytest.mark.parametrize("launcher", _CANONICAL_LAUNCHERS)
def test_candidate_launcher_canonical_renderer_matches_authority(
    launcher: str,
) -> None:
    namespace = runpy.run_path(str(TEMPLATES / launcher))
    render = namespace["_canonical_json_bytes"]
    value = {"unicode": "核", "nested": {"b": 2, "a": 1}}

    assert render(value) == canonical_json_bytes(value)
    with pytest.raises(ValueError):
        render({"invalid": float("nan")})


@dataclass
class _FreshActionRunner:
    target_schema: str = "schema-v1"
    prior_schema: str | None = None
    known_schema_revisions: tuple[str, ...] = ()
    compatibility: str | None = None

    def run(self, request):
        target_schema = (
            self.prior_schema
            if request.phase == "observe" and self.prior_schema is not None
            else self.target_schema
        )
        known_schema_revisions = tuple(
            head for head in self.known_schema_revisions if head != target_schema
        )
        native = {
            component: request.candidate_active[component]
            for component in (
                "platform",
                "encode-runtime",
                "bulk-rnaseq-runtime",
            )
        }
        readiness = {
            name: ReadinessCheck("ready", "READY", IDENTITY)
            for name in VERIFICATION_CHECKS
        }
        for component, check in (
            ("platform", "platform-native"),
            ("encode-runtime", "encode-runtime-native"),
            ("bulk-rnaseq-runtime", "bulk-runtime-native"),
        ):
            evidence = native[component]
            readiness[check] = (
                ReadinessCheck("not-applicable", "COMPONENT_NOT_ACTIVE")
                if evidence is None
                else ReadinessCheck("ready", "READY", evidence)
            )
        readiness["database-schema"] = ReadinessCheck(
            "unavailable", "SCHEMA_INCOMPATIBLE"
        )
        return DeploymentActionReceipt.create(
            request_identity=request.identity,
            status=("observed" if request.phase == "observe" else "admitted"),
            compatibility=(
                self.compatibility
                or (
                    "compatible"
                    if all(
                        value is not None for value in request.candidate_active.values()
                    )
                    else "incomplete"
                )
            ),
            accepted_schema_heads=(),
            target_schema_heads=(target_schema,),
            migration_inventory_identity=IDENTITY,
            known_schema_revisions=known_schema_revisions,
            migration_required=True,
            rollback_supported=False,
            api_contract_sha256="a" * 64,
            native_identities=native,
            frontend_identity=native["platform"],
            reference_compatibility_identity=IDENTITY,
            readiness=readiness,
        )


class _FreshRuntimePreparer:
    def prepare(self, request):
        return EncodeRuntimePrepareReceipt.create(
            request_identity=request.identity,
            deployment_identity=request.deployment_identity,
            inventory=_encode_runtime_inventory(),
        )


@dataclass
class _FreshBulkRuntimePreparer:
    requests: list[object] | None = None

    def observe_boundary(self) -> BulkDockerBoundary:
        return BulkDockerBoundary(IDENTITY, SERVICE_IDENTITY, os.getuid(), os.getgid())

    def prepare(self, request):
        if self.requests is not None:
            self.requests.append(request)
        return BulkRuntimePrepareReceipt.create(
            request_identity=request.identity,
            candidate_bulk_identity=request.candidate_bulk_identity,
            runtime_identity=THIRD_IDENTITY,
            image_set_identity=SERVICE_IDENTITY,
            image_count=1,
        )


@dataclass
class _FreshDatabasePreparer:
    layout: DeploymentLayout
    owner_uid: int
    owner_gid: int
    requests: list[DatabasePrepareRequest] | None = None
    fail_after_write: bool = False

    def prepare(self, request) -> DatabasePrepareReceipt:
        database = (
            fresh_database_candidate_path(self.layout, request.task_identity)
            if request.database_mode == "fresh-candidate"
            else self.layout.database
        )
        database.parent.mkdir(parents=True, exist_ok=True)
        database.parent.chmod(0o2770)
        before_identity = None
        if request.database_mode == "existing-live":
            before_identity = database_content_identity(
                inspect_database(
                    database,
                    expected_owner_uid=self.owner_uid,
                    expected_owner_gid=self.owner_gid,
                )
            )
        with sqlite3.connect(database) as connection:
            if request.database_mode == "fresh-candidate":
                connection.execute(
                    "CREATE TABLE alembic_version (version_num VARCHAR(32) NOT NULL)"
                )
            else:
                connection.execute("DELETE FROM alembic_version")
            connection.execute(
                "INSERT INTO alembic_version (version_num) VALUES (?)",
                (request.target_schema_heads[0],),
            )
        database.chmod(0o660)
        if self.requests is not None:
            self.requests.append(request)
        if self.fail_after_write:
            raise fail(
                "TEST_DATABASE_PREPARE_FAILED",
                "Database preparation failed.",
                recoverable=True,
            )
        inspection = inspect_database(
            database,
            expected_owner_uid=self.owner_uid,
            expected_owner_gid=self.owner_gid,
        )
        return DatabasePrepareReceipt.create(
            request_identity=request.identity,
            database_before_identity=before_identity,
            database_after_identity=database_content_identity(inspection),
            schema_heads=inspection.schema_heads,
        )


@dataclass
class _FreshConfiguration:
    states: list[str]
    layout: DeploymentLayout
    failure_code: str | None = None

    def activate(self, *, state, api_contract_sha256):
        assert api_contract_sha256 == "a" * 64
        if self.failure_code is not None:
            raise fail(
                self.failure_code,
                "Test configuration failed.",
                recoverable=True,
            )
        self.states.append(state.identity)
        return render_platform_environment(
            self.layout,
            state,
            api_contract_sha256=api_contract_sha256,
        )


@dataclass
class _FreshServices:
    started: list[str]
    running: dict[str, object] | None = None

    def __post_init__(self) -> None:
        if self.running is None:
            self.running = {}

    def status(self, request):
        assert self.running is not None
        return self.running.get(request.unit)

    def stop(self, request, *, cleanup):
        assert self.running is not None
        assert not cleanup
        service = self.running.get(request.unit)
        if service is None or service.identity != request.service_identity:
            raise AssertionError("fresh service identity mismatch")
        del self.running[request.unit]

    def start(self, request):
        assert self.running is not None
        self.started.append(request.unit)
        identity = _TrackingServices._identity(request)
        self.running[request.unit] = identity
        return identity


@dataclass
class _TrackingServices:
    running: dict[str, object]
    calls: list[tuple[str, str, str]]
    fail_start_once: str | None = None
    start_failed: bool = False

    @staticmethod
    def _identity(request: OperatorRequest):
        content = (
            f"{request.unit}:{request.deployment_identity}:{request.task_identity}"
        ).encode()
        return SimpleNamespace(
            identity=f"sha256-{hashlib.sha256(content).hexdigest()}",
            unit=request.unit,
            deployment_identity=request.deployment_identity,
            task_identity=request.task_identity,
        )

    def status(self, request: OperatorRequest):
        assert request.unit is not None
        self.calls.append(("status", request.unit, request.deployment_identity))
        service = self.running.get(request.unit)
        if service is None:
            return None
        if service.deployment_identity != request.deployment_identity:
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_MISMATCH",
                "Service identity does not match this deployment.",
            )
        return service

    def stop(self, request: OperatorRequest, *, cleanup: bool):
        assert request.unit is not None and not cleanup
        self.calls.append(("stop", request.unit, request.deployment_identity))
        service = self.running.get(request.unit)
        if (
            service is None
            or service.deployment_identity != request.deployment_identity
            or service.identity != request.service_identity
        ):
            raise fail(
                "OPERATOR_SERVICE_IDENTITY_MISMATCH",
                "Service identity does not match this deployment.",
            )
        del self.running[request.unit]

    def start(self, request: OperatorRequest):
        assert request.unit is not None
        self.calls.append(("start", request.unit, request.deployment_identity))
        if self.fail_start_once == request.unit and not self.start_failed:
            self.start_failed = True
            raise fail(
                "TEST_SERVICE_START_FAILED",
                "Test service start failed.",
                recoverable=True,
            )
        if request.unit in self.running:
            raise fail("OPERATOR_SERVICE_ALREADY_RUNNING", "Service is running.")
        service = self._identity(request)
        self.running[request.unit] = service
        return service


def _supported_state_store(
    layout: DeploymentLayout,
    owner_uid: int,
    owner_gid: int,
) -> StateStore:
    return StateStore(
        layout,
        reader_gid=owner_gid,
        service_gid=owner_gid,
    )


def _state_with_active_components(
    layout: DeploymentLayout,
    *,
    owner_uid: int,
    owner_gid: int,
    staged_component: str | None = None,
    staged_identity: str | None = None,
) -> tuple[StateStore, DeploymentState]:
    states = _supported_state_store(layout, owner_uid, owner_gid)
    with states.transaction(
        exclusive=True,
        expected_owner_uid=owner_uid,
        expected_owner_gid=owner_gid,
    ) as transaction:
        state = transaction.initialize()
        for component, identity in (
            ("platform", OLD_PLATFORM_IDENTITY),
            ("encode-runtime", OLD_ENCODE_IDENTITY),
            ("bulk-rnaseq-runtime", OLD_BULK_IDENTITY),
        ):
            staged = state.stage(component, identity)
            transaction.commit(
                staged,
                operation=f"stage-{component}",
                expected_current_identity=state.identity,
            )
            state = staged
            active = state.activate(component)
            transaction.commit(
                active,
                operation=f"activate-{component}",
                expected_current_identity=state.identity,
                platform_environment=render_platform_environment(
                    layout,
                    active,
                    api_contract_sha256="a" * 64,
                ),
            )
            state = active
        if staged_component is not None:
            assert staged_identity is not None
            staged = state.stage(staged_component, staged_identity)
            transaction.commit(
                staged,
                operation=f"stage-{staged_component}",
                expected_current_identity=state.identity,
            )
            state = staged
    return states, state


def test_stage_state_commit_handoff_crash_is_reconciled_as_complete(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    states = _supported_state_store(layout, owner_uid, owner_gid)
    with states.transaction(
        exclusive=True,
        expected_owner_uid=owner_uid,
        expected_owner_gid=owner_gid,
    ) as transaction:
        prior = transaction.initialize()
        candidate = prior.stage("platform", IDENTITY)
        transaction.commit(
            candidate,
            operation="stage-platform",
            expected_current_identity=prior.identity,
        )
    active = {
        component: prior.components[component].active
        for component in ("platform", "encode-runtime", "bulk-rnaseq-runtime")
    }
    seed = OperatorJournalStore(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
    )
    seed._directories(create=True)
    interrupted = OperatorTransaction.create(
        request_identity=THIRD_IDENTITY,
        operation="stage",
        task_identity=TASK_IDENTITY,
        deployment_identity=IDENTITY,
        component="platform",
        unit=None,
        phase="candidate-selected",
        prior_active=active,
        candidate_active=active,
        prior_state_identity=prior.identity,
        candidate_state_identity=candidate.identity,
    )
    seed._write_new(layout.operator_transaction_active, interrupted)
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=FakeServiceController(),
        root_uid=owner_uid,
        root_gid=owner_gid,
    )
    journals = OperatorJournalStore(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        recovery_controller=controller,
    )
    next_task = "task-" + "8" * 32

    with journals.operation(
        operation="stage",
        task_identity=next_task,
        deployment_identity=THIRD_IDENTITY,
        component="encode-runtime",
        unit=None,
    ) as journal:
        journal.complete()

    recovered = json.loads(
        (layout.operator_transaction_history / f"{TASK_IDENTITY}.json").read_text()
    )
    assert recovered["phase"] == "complete"
    assert recovered["failure_phase"] == "candidate-selected"
    assert states.read().identity == candidate.identity
    assert states.read().components["platform"].staged == IDENTITY
    assert not layout.operator_transaction_active.exists()


def test_direct_start_recovery_idempotently_finishes_requested_service(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    states, state = _state_with_active_components(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
    )
    services = _TrackingServices({}, [])
    active = {
        component: state.components[component].active
        for component in ("platform", "encode-runtime", "bulk-rnaseq-runtime")
    }
    record = OperatorTransaction.create(
        request_identity=IDENTITY,
        operation="start",
        task_identity=TASK_IDENTITY,
        deployment_identity=OLD_PLATFORM_IDENTITY,
        component=None,
        unit="helixweave-worker.service",
        phase="recovery-required",
        failure_phase="service-starting",
        point_of_no_return=True,
        restart_units=("helixweave-worker.service",),
        prior_active=active,
        candidate_active=active,
        prior_state_identity=state.identity,
        candidate_state_identity=state.identity,
    )
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=services,
        root_uid=owner_uid,
        root_gid=owner_gid,
    )

    recovered = controller.recover(record)

    assert recovered.phase == "complete"
    assert services.running["helixweave-worker.service"].deployment_identity == (
        OLD_PLATFORM_IDENTITY
    )
    assert [call[:2] for call in services.calls] == [
        ("status", "helixweave-worker.service"),
        ("start", "helixweave-worker.service"),
    ]


def test_direct_start_failure_after_side_effect_is_synchronously_recovered(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    states, state = _state_with_active_components(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
    )

    class FailAfterStart(_TrackingServices):
        def start(self, request: OperatorRequest):
            super().start(request)
            raise fail(
                "TEST_AFTER_START_FAILED",
                "Test start failed.",
                recoverable=True,
            )

    services = FailAfterStart({}, [])
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=services,
        root_uid=owner_uid,
        root_gid=owner_gid,
    )
    journals = OperatorJournalStore(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        recovery_controller=controller,
    )
    active = {
        component: state.components[component].active
        for component in ("platform", "encode-runtime", "bulk-rnaseq-runtime")
    }

    with pytest.raises(DeploymentError) as captured:
        with journals.operation(
            operation="start",
            task_identity=TASK_IDENTITY,
            deployment_identity=OLD_PLATFORM_IDENTITY,
            component=None,
            unit="helixweave-worker.service",
        ) as journal:
            journal.advance(
                "service-starting",
                point_of_no_return=True,
                restart_units=("helixweave-worker.service",),
                prior_active=active,
                candidate_active=active,
                evidence={
                    "prior_state_identity": state.identity,
                    "candidate_state_identity": state.identity,
                },
            )
            services.start(
                OperatorRequest(
                    operation="start",
                    task_identity=TASK_IDENTITY,
                    deployment_identity=OLD_PLATFORM_IDENTITY,
                    unit="helixweave-worker.service",
                )
            )

    assert captured.value.issue.code == "TEST_AFTER_START_FAILED"
    assert services.running["helixweave-worker.service"].deployment_identity == (
        OLD_PLATFORM_IDENTITY
    )
    archived = json.loads(
        (layout.operator_transaction_history / f"{TASK_IDENTITY}.json").read_text()
    )
    assert archived["phase"] == "complete"
    assert archived["failure_phase"] == "service-starting"


def test_direct_stop_recovery_idempotently_finishes_requested_stop(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    states, state = _state_with_active_components(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
    )
    services = _TrackingServices({}, [])
    request = OperatorRequest(
        operation="start",
        task_identity=TASK_IDENTITY,
        deployment_identity=OLD_PLATFORM_IDENTITY,
        unit="helixweave-worker.service",
    )
    service = services._identity(request)
    services.running["helixweave-worker.service"] = service
    active = {
        component: state.components[component].active
        for component in ("platform", "encode-runtime", "bulk-rnaseq-runtime")
    }
    record = OperatorTransaction.create(
        request_identity=IDENTITY,
        operation="stop",
        task_identity=TASK_IDENTITY,
        deployment_identity=OLD_PLATFORM_IDENTITY,
        component=None,
        unit="helixweave-worker.service",
        phase="recovery-required",
        failure_phase="service-stopping",
        restart_units=("helixweave-worker.service",),
        prior_running_units=("helixweave-worker.service",),
        prior_active=active,
        candidate_active=active,
        prior_state_identity=state.identity,
        candidate_state_identity=state.identity,
        evidence={"service_identity": service.identity},
    )
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=services,
        root_uid=owner_uid,
        root_gid=owner_gid,
    )

    recovered = controller.recover(record)

    assert recovered.phase == "complete"
    assert "helixweave-worker.service" not in services.running
    assert [call[:2] for call in services.calls] == [
        ("status", "helixweave-worker.service"),
        ("stop", "helixweave-worker.service"),
    ]


def test_uninstall_recovery_restores_the_boundary_after_stopping_services(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    states, state = _state_with_active_components(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
    )
    services = _TrackingServices({}, [])
    active = {
        component: state.components[component].active
        for component in ("platform", "encode-runtime", "bulk-rnaseq-runtime")
    }

    @dataclass
    class Boundary:
        uninstall_calls: int = 0
        recover_calls: int = 0

        def uninstall(self, *, before_control_removal) -> None:
            self.uninstall_calls += 1
            before_control_removal()

        def recover(self) -> None:
            self.recover_calls += 1

    boundary = Boundary()
    record = OperatorTransaction.create(
        request_identity=IDENTITY,
        operation="uninstall",
        task_identity=TASK_IDENTITY,
        deployment_identity=state.identity,
        component=None,
        unit=None,
        phase="recovery-required",
        failure_phase="writers-stopped",
        write_fence=True,
        prior_active=active,
        candidate_active=active,
        prior_state_identity=state.identity,
        candidate_state_identity=state.identity,
    )
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=services,
        boundary_uninstaller=boundary,
        root_uid=owner_uid,
        root_gid=owner_gid,
    )

    recovered = controller.recover(record)

    assert recovered.phase == "aborted"
    assert boundary.uninstall_calls == 0
    assert boundary.recover_calls == 1


def test_uninstall_failure_is_synchronously_restored_and_journal_aborted(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    states, state = _state_with_active_components(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
    )
    services = _TrackingServices({}, [])

    @dataclass
    class FailOnceBoundary:
        uninstall_calls: int = 0
        recover_calls: int = 0

        def uninstall(self, *, before_control_removal) -> None:
            del before_control_removal
            self.uninstall_calls += 1
            raise fail(
                "TEST_UNINSTALL_FAILED",
                "Test uninstall failed.",
                recoverable=True,
            )

        def recover(self) -> None:
            self.recover_calls += 1

    boundary = FailOnceBoundary()
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=services,
        boundary_uninstaller=boundary,
        root_uid=owner_uid,
        root_gid=owner_gid,
    )
    journals = OperatorJournalStore(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        recovery_controller=controller,
    )
    request = OperatorRequest(
        operation="uninstall",
        task_identity=TASK_IDENTITY,
        deployment_identity=state.identity,
    )

    with pytest.raises(DeploymentError) as captured:
        with journals.operation(
            operation="uninstall",
            task_identity=TASK_IDENTITY,
            deployment_identity=state.identity,
            component=None,
            unit=None,
        ) as journal:
            controller.execute(request, journal=journal)

    assert captured.value.issue.code == "TEST_UNINSTALL_FAILED"
    assert boundary.uninstall_calls == 1
    assert boundary.recover_calls == 1
    assert not layout.operator_transaction_active.exists()
    archived = json.loads(
        (layout.operator_transaction_history / f"{TASK_IDENTITY}.json").read_text()
    )
    assert archived["phase"] == "aborted"
    assert archived["failure_phase"] == "writers-stopped"


def _write_operator_test_database(
    layout: DeploymentLayout,
    *,
    heads: tuple[str, ...],
    value: str,
) -> None:
    layout.run_root.mkdir(parents=True, mode=0o755, exist_ok=True)
    layout.run_root.chmod(0o755)
    layout.database.parent.mkdir(parents=True, exist_ok=True)
    layout.database.parent.chmod(0o2770)
    with sqlite3.connect(layout.database) as connection:
        connection.execute(
            "CREATE TABLE alembic_version (version_num VARCHAR(128) NOT NULL)"
        )
        connection.executemany(
            "INSERT INTO alembic_version (version_num) VALUES (?)",
            ((head,) for head in heads),
        )
        connection.execute("CREATE TABLE durable_state (value VARCHAR(128) NOT NULL)")
        connection.execute("INSERT INTO durable_state (value) VALUES (?)", (value,))
    layout.database.chmod(0o660)


def _database_value(layout: DeploymentLayout) -> str:
    with sqlite3.connect(f"file:{layout.database}?mode=ro", uri=True) as connection:
        row = connection.execute("SELECT value FROM durable_state").fetchone()
    assert row is not None
    return str(row[0])


def test_observation_preserves_path_free_journal_summary_when_database_is_absent(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    state = StateStore(
        layout,
        reader_gid=owner_gid,
        service_gid=owner_gid,
    ).initialize(
        expected_owner_uid=owner_uid,
        expected_owner_gid=owner_gid,
    )
    provider = FixedObservationProvider(
        layout,
        _FreshServices([]),
        root_uid=owner_uid,
        root_gid=owner_gid,
        operator_group_gid=owner_gid,
        service_uid=owner_uid,
        service_gid=owner_gid,
    )
    journals = OperatorJournalStore(layout, owner_uid=owner_uid, owner_gid=owner_gid)

    with journals.operation(
        operation="stage",
        task_identity=TASK_IDENTITY,
        deployment_identity=IDENTITY,
        component="platform",
        unit=None,
    ) as journal:
        observation = provider.observe(
            OperatorRequest(
                operation="observe",
                deployment_identity=state.identity,
                task_identity=TASK_IDENTITY,
            )
        ).with_operator_journal(journals.summary())
        journal.complete()

    assert observation.database_schema_identity is None
    assert observation.database_schema_heads == ()
    assert observation.operator_pending_count == 1
    assert observation.operator_recovery_required_count == 0


def test_fresh_three_component_assembly_commits_partial_state_without_starting(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    monkeypatch.setattr(database_module, "_ROOT_UID", owner_uid)
    monkeypatch.setattr(database_module, "_ROOT_GID", owner_gid)
    states = StateStore(
        layout,
        reader_gid=owner_gid,
        service_gid=owner_gid,
    )
    with states.transaction(
        exclusive=True,
        expected_owner_uid=owner_uid,
        expected_owner_gid=owner_gid,
    ) as transaction:
        state = transaction.initialize()
        for component, identity in (
            ("platform", IDENTITY),
            ("encode-runtime", SERVICE_IDENTITY),
            ("bulk-rnaseq-runtime", THIRD_IDENTITY),
        ):
            staged = state.stage(component, identity)
            transaction.commit(
                staged,
                operation=f"stage-{component}",
                expected_current_identity=state.identity,
            )
            state = staged
    services = _FreshServices([])
    configuration = _FreshConfiguration([], layout)
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=services,
        action_runner=_FreshActionRunner(),
        database_preparer=_FreshDatabasePreparer(layout, owner_uid, owner_gid),
        encode_runtime_preparer=_FreshRuntimePreparer(),
        bulk_runtime_preparer=_FreshBulkRuntimePreparer(),
        configuration=configuration,
        root_uid=owner_uid,
        root_gid=owner_gid,
        service_uid=owner_uid,
        service_gid=owner_gid,
    )
    journals = OperatorJournalStore(layout, owner_uid=owner_uid, owner_gid=owner_gid)

    for index, (component, identity) in enumerate(
        (
            ("platform", IDENTITY),
            ("encode-runtime", SERVICE_IDENTITY),
            ("bulk-rnaseq-runtime", THIRD_IDENTITY),
        ),
        start=1,
    ):
        task = f"task-{index:032x}"
        request = OperatorRequest(
            operation="activate",
            component=component,
            deployment_identity=identity,
            task_identity=task,
        )
        with journals.operation(
            operation="activate",
            task_identity=task,
            deployment_identity=identity,
            component=component,
            unit=None,
        ) as journal:
            outcome = controller.execute(request, journal=journal)
            journal.advance("complete")
        assert outcome.state == "activated"
        if index < 3:
            assert services.started == []

    assert services.started == list(SERVICE_UNITS)
    final = states.read(
        expected_owner_uid=owner_uid,
        expected_owner_gid=owner_gid,
    )
    assert all(
        final.components[component].active is not None for component in final.components
    )


@pytest.mark.parametrize(
    ("current_head", "migration_count"),
    (("schema-v0", 1), ("schema-v1", 0)),
)
def test_legacy_database_is_backed_up_before_adoption_or_forward_migration(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    current_head: str,
    migration_count: int,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    monkeypatch.setattr(database_module, "_ROOT_UID", owner_uid)
    monkeypatch.setattr(database_module, "_ROOT_GID", owner_gid)
    states = _supported_state_store(layout, owner_uid, owner_gid)
    with states.transaction(
        exclusive=True,
        expected_owner_uid=owner_uid,
        expected_owner_gid=owner_gid,
    ) as transaction:
        prior = transaction.initialize()
        staged = prior.stage("platform", IDENTITY)
        transaction.commit(
            staged,
            operation="stage-platform",
            expected_current_identity=prior.identity,
        )
    _write_operator_test_database(
        layout,
        heads=(current_head,),
        value="legacy-preserved",
    )
    prepare_requests: list[DatabasePrepareRequest] = []
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=_TrackingServices({}, []),
        action_runner=_FreshActionRunner(
            target_schema="schema-v1",
            prior_schema="schema-v0",
            known_schema_revisions=("schema-v0",),
        ),
        database_preparer=_FreshDatabasePreparer(
            layout,
            owner_uid,
            owner_gid,
            requests=prepare_requests,
        ),
        encode_runtime_preparer=_FreshRuntimePreparer(),
        configuration=_FreshConfiguration([], layout),
        root_uid=owner_uid,
        root_gid=owner_gid,
        service_uid=owner_uid,
        service_gid=owner_gid,
    )
    journals = OperatorJournalStore(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        recovery_controller=controller,
    )
    request = OperatorRequest(
        operation="activate",
        component="platform",
        deployment_identity=IDENTITY,
        task_identity=TASK_IDENTITY,
    )

    with journals.operation(
        operation="activate",
        task_identity=TASK_IDENTITY,
        deployment_identity=IDENTITY,
        component="platform",
        unit=None,
    ) as journal:
        outcome = controller.execute(request, journal=journal)
        journal.complete()

    assert outcome.state == "activated"
    assert len(prepare_requests) == migration_count
    assert inspect_database(
        layout.database,
        expected_owner_uid=owner_uid,
        expected_owner_gid=owner_gid,
    ).schema_heads == ("schema-v1",)
    assert _database_value(layout) == "legacy-preserved"
    assert any(layout.database_backups.rglob("receipt.json"))


@pytest.mark.parametrize(
    "heads",
    (("unknown-head",), ("schema-v0", "schema-v1")),
)
def test_legacy_database_rejects_unknown_or_multiple_heads_without_mutation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    heads: tuple[str, ...],
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    monkeypatch.setattr(database_module, "_ROOT_UID", owner_uid)
    monkeypatch.setattr(database_module, "_ROOT_GID", owner_gid)
    states = _supported_state_store(layout, owner_uid, owner_gid)
    with states.transaction(
        exclusive=True,
        expected_owner_uid=owner_uid,
        expected_owner_gid=owner_gid,
    ) as transaction:
        prior = transaction.initialize()
        staged = prior.stage("platform", IDENTITY)
        transaction.commit(
            staged,
            operation="stage-platform",
            expected_current_identity=prior.identity,
        )
    _write_operator_test_database(layout, heads=heads, value="untouched")
    before = database_content_identity(
        inspect_database(
            layout.database,
            expected_owner_uid=owner_uid,
            expected_owner_gid=owner_gid,
        )
    )
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=_TrackingServices({}, []),
        action_runner=_FreshActionRunner(
            target_schema="schema-v1",
            prior_schema="schema-v0",
            known_schema_revisions=("schema-v0",),
        ),
        database_preparer=_FreshDatabasePreparer(layout, owner_uid, owner_gid),
        encode_runtime_preparer=_FreshRuntimePreparer(),
        configuration=_FreshConfiguration([], layout),
        root_uid=owner_uid,
        root_gid=owner_gid,
        service_uid=owner_uid,
        service_gid=owner_gid,
    )
    journals = OperatorJournalStore(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        recovery_controller=controller,
    )

    with pytest.raises(DeploymentError) as caught:
        with journals.operation(
            operation="activate",
            task_identity=TASK_IDENTITY,
            deployment_identity=IDENTITY,
            component="platform",
            unit=None,
        ) as journal:
            controller.execute(
                OperatorRequest(
                    operation="activate",
                    component="platform",
                    deployment_identity=IDENTITY,
                    task_identity=TASK_IDENTITY,
                ),
                journal=journal,
            )

    assert caught.value.issue.code == "DEPLOYMENT_SCHEMA_INCOMPATIBLE"
    assert (
        database_content_identity(
            inspect_database(
                layout.database,
                expected_owner_uid=owner_uid,
                expected_owner_gid=owner_gid,
            )
        )
        == before
    )
    assert (
        states.read(
            expected_owner_uid=owner_uid,
            expected_owner_gid=owner_gid,
        ).identity
        == staged.identity
    )
    history = json.loads(
        (layout.operator_transaction_history / f"{TASK_IDENTITY}.json").read_text()
    )
    assert history["phase"] == "aborted"


def test_pre_ponr_failure_synchronously_restores_database_state_and_writers(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    monkeypatch.setattr(database_module, "_ROOT_UID", owner_uid)
    monkeypatch.setattr(database_module, "_ROOT_GID", owner_gid)
    states, prior = _state_with_active_components(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        staged_component="platform",
        staged_identity=IDENTITY,
    )
    _write_operator_test_database(
        layout,
        heads=("schema-v0",),
        value="prior-data",
    )
    services = _TrackingServices({}, [])
    for unit in ("helixweave-api.service", "helixweave-worker.service"):
        service_request = OperatorRequest(
            operation="start",
            unit=unit,
            deployment_identity=OLD_PLATFORM_IDENTITY,
            task_identity=f"task-{'9' * 32}",
        )
        services.running[unit] = services._identity(service_request)
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=services,
        action_runner=_FreshActionRunner(
            target_schema="schema-v1",
            prior_schema="schema-v0",
            known_schema_revisions=("schema-v0",),
        ),
        database_preparer=_FreshDatabasePreparer(layout, owner_uid, owner_gid),
        encode_runtime_preparer=_FreshRuntimePreparer(),
        configuration=_FreshConfiguration(
            [], layout, failure_code="TEST_CONFIGURATION_FAILED"
        ),
        root_uid=owner_uid,
        root_gid=owner_gid,
        service_uid=owner_uid,
        service_gid=owner_gid,
    )
    journals = OperatorJournalStore(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        recovery_controller=controller,
    )

    with pytest.raises(DeploymentError) as caught:
        with journals.operation(
            operation="activate",
            task_identity=TASK_IDENTITY,
            deployment_identity=IDENTITY,
            component="platform",
            unit=None,
        ) as journal:
            controller.execute(
                OperatorRequest(
                    operation="activate",
                    component="platform",
                    deployment_identity=IDENTITY,
                    task_identity=TASK_IDENTITY,
                ),
                journal=journal,
            )

    assert caught.value.issue.code == "TEST_CONFIGURATION_FAILED"
    assert (
        states.read(
            expected_owner_uid=owner_uid,
            expected_owner_gid=owner_gid,
        ).identity
        == prior.identity
    )
    database = inspect_database(
        layout.database,
        expected_owner_uid=owner_uid,
        expected_owner_gid=owner_gid,
    )
    assert database.schema_heads == ("schema-v0",)
    assert _database_value(layout) == "prior-data"
    assert {
        unit: service.deployment_identity for unit, service in services.running.items()
    } == {
        "helixweave-api.service": OLD_PLATFORM_IDENTITY,
        "helixweave-worker.service": OLD_PLATFORM_IDENTITY,
    }
    assert [
        (action, unit)
        for action, unit, _deployment in services.calls
        if action in {"stop", "start"}
    ] == [
        ("stop", "helixweave-api.service"),
        ("stop", "helixweave-worker.service"),
        ("start", "helixweave-api.service"),
        ("start", "helixweave-worker.service"),
    ]
    history = json.loads(
        (layout.operator_transaction_history / f"{TASK_IDENTITY}.json").read_text()
    )
    assert history["phase"] == "aborted"
    assert history["failure_phase"] == "migration-verified"
    assert not layout.operator_transaction_active.exists()


def test_post_ponr_failure_finishes_candidate_without_restoring_prior(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    monkeypatch.setattr(database_module, "_ROOT_UID", owner_uid)
    monkeypatch.setattr(database_module, "_ROOT_GID", owner_gid)
    states, prior = _state_with_active_components(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        staged_component="platform",
        staged_identity=IDENTITY,
    )
    _write_operator_test_database(
        layout,
        heads=("schema-v1",),
        value="candidate-data",
    )
    services = _TrackingServices({}, [], fail_start_once="helixweave-api.service")
    for unit in ("helixweave-api.service", "helixweave-worker.service"):
        service_request = OperatorRequest(
            operation="start",
            unit=unit,
            deployment_identity=OLD_PLATFORM_IDENTITY,
            task_identity=f"task-{'8' * 32}",
        )
        services.running[unit] = services._identity(service_request)
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=services,
        action_runner=_FreshActionRunner(target_schema="schema-v1"),
        database_preparer=_FreshDatabasePreparer(layout, owner_uid, owner_gid),
        encode_runtime_preparer=_FreshRuntimePreparer(),
        configuration=_FreshConfiguration([], layout),
        root_uid=owner_uid,
        root_gid=owner_gid,
        service_uid=owner_uid,
        service_gid=owner_gid,
    )
    journals = OperatorJournalStore(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        recovery_controller=controller,
    )

    with pytest.raises(DeploymentError) as caught:
        with journals.operation(
            operation="activate",
            task_identity=TASK_IDENTITY,
            deployment_identity=IDENTITY,
            component="platform",
            unit=None,
        ) as journal:
            controller.execute(
                OperatorRequest(
                    operation="activate",
                    component="platform",
                    deployment_identity=IDENTITY,
                    task_identity=TASK_IDENTITY,
                ),
                journal=journal,
            )

    assert caught.value.issue.code == "TEST_SERVICE_START_FAILED"
    recovered = states.read(
        expected_owner_uid=owner_uid,
        expected_owner_gid=owner_gid,
    )
    assert recovered.components["platform"].active == IDENTITY
    assert recovered.components["platform"].previous == OLD_PLATFORM_IDENTITY
    assert _database_value(layout) == "candidate-data"
    assert {
        unit: service.deployment_identity for unit, service in services.running.items()
    } == {
        "helixweave-api.service": IDENTITY,
        "helixweave-worker.service": IDENTITY,
    }
    assert [
        (action, unit)
        for action, unit, _deployment in services.calls
        if action == "start"
    ] == [
        ("start", "helixweave-api.service"),
        ("start", "helixweave-api.service"),
        ("start", "helixweave-worker.service"),
    ]
    history = json.loads(
        (layout.operator_transaction_history / f"{TASK_IDENTITY}.json").read_text()
    )
    assert history["phase"] == "complete"
    assert history["point_of_no_return"] is True
    assert history["failure_phase"] == "service-starting"


def test_failed_fresh_prepare_is_quarantined_before_returning_to_caller(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    monkeypatch.setattr(database_module, "_ROOT_UID", owner_uid)
    monkeypatch.setattr(database_module, "_ROOT_GID", owner_gid)
    states = _supported_state_store(layout, owner_uid, owner_gid)
    with states.transaction(
        exclusive=True,
        expected_owner_uid=owner_uid,
        expected_owner_gid=owner_gid,
    ) as transaction:
        initial = transaction.initialize()
        prior = initial.stage("platform", IDENTITY)
        transaction.commit(
            prior,
            operation="stage-platform",
            expected_current_identity=initial.identity,
        )
    services = _TrackingServices({}, [])
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=services,
        action_runner=_FreshActionRunner(),
        database_preparer=_FreshDatabasePreparer(
            layout,
            owner_uid,
            owner_gid,
            fail_after_write=True,
        ),
        encode_runtime_preparer=_FreshRuntimePreparer(),
        configuration=_FreshConfiguration([], layout),
        root_uid=owner_uid,
        root_gid=owner_gid,
        service_uid=owner_uid,
        service_gid=owner_gid,
    )
    journals = OperatorJournalStore(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        recovery_controller=controller,
    )

    with pytest.raises(DeploymentError) as caught:
        with journals.operation(
            operation="activate",
            task_identity=TASK_IDENTITY,
            deployment_identity=IDENTITY,
            component="platform",
            unit=None,
        ) as journal:
            controller.execute(
                OperatorRequest(
                    operation="activate",
                    component="platform",
                    deployment_identity=IDENTITY,
                    task_identity=TASK_IDENTITY,
                ),
                journal=journal,
            )

    assert caught.value.issue.code == "TEST_DATABASE_PREPARE_FAILED"
    assert not layout.database.exists()
    assert not fresh_database_candidate_path(layout, TASK_IDENTITY).exists()
    assert (
        states.read(
            expected_owner_uid=owner_uid,
            expected_owner_gid=owner_gid,
        ).identity
        == prior.identity
    )
    assert any(
        (layout.data_root / "operator" / "database-recovery").rglob("receipt.json")
    )
    history = json.loads(
        (layout.operator_transaction_history / f"{TASK_IDENTITY}.json").read_text()
    )
    assert history["phase"] == "aborted"
    assert history["failure_phase"] == "migration-started"


def test_encode_upgrade_restarts_only_api_and_worker_in_dependency_order(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    monkeypatch.setattr(database_module, "_ROOT_UID", owner_uid)
    monkeypatch.setattr(database_module, "_ROOT_GID", owner_gid)
    states, _prior = _state_with_active_components(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        staged_component="encode-runtime",
        staged_identity=IDENTITY,
    )
    _write_operator_test_database(
        layout,
        heads=("schema-v1",),
        value="preserved",
    )
    services = _TrackingServices({}, [])
    deployments = {
        "helixweave-redis.service": OLD_PLATFORM_IDENTITY,
        "helixweave-docker-rootless.service": OLD_BULK_IDENTITY,
        "helixweave-api.service": OLD_PLATFORM_IDENTITY,
        "helixweave-worker.service": OLD_PLATFORM_IDENTITY,
    }
    for unit, deployment in deployments.items():
        service_request = OperatorRequest(
            operation="start",
            unit=unit,
            deployment_identity=deployment,
            task_identity=f"task-{'7' * 32}",
        )
        services.running[unit] = services._identity(service_request)
    daemon_identities = {
        unit: services.running[unit].identity
        for unit in (
            "helixweave-redis.service",
            "helixweave-docker-rootless.service",
        )
    }
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=services,
        action_runner=_FreshActionRunner(target_schema="schema-v1"),
        database_preparer=_FreshDatabasePreparer(layout, owner_uid, owner_gid),
        encode_runtime_preparer=_FreshRuntimePreparer(),
        configuration=_FreshConfiguration([], layout),
        root_uid=owner_uid,
        root_gid=owner_gid,
        service_uid=owner_uid,
        service_gid=owner_gid,
    )
    journals = OperatorJournalStore(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        recovery_controller=controller,
    )

    with journals.operation(
        operation="activate",
        task_identity=TASK_IDENTITY,
        deployment_identity=IDENTITY,
        component="encode-runtime",
        unit=None,
    ) as journal:
        outcome = controller.execute(
            OperatorRequest(
                operation="activate",
                component="encode-runtime",
                deployment_identity=IDENTITY,
                task_identity=TASK_IDENTITY,
            ),
            journal=journal,
        )
        journal.complete()

    assert outcome.state == "activated"
    assert [
        (action, unit)
        for action, unit, _deployment in services.calls
        if action in {"stop", "start"}
    ] == [
        ("stop", "helixweave-api.service"),
        ("stop", "helixweave-worker.service"),
        ("start", "helixweave-api.service"),
        ("start", "helixweave-worker.service"),
    ]
    assert not any(
        unit in {"helixweave-redis.service", "helixweave-docker-rootless.service"}
        for _action, unit, _deployment in services.calls
    )
    assert {
        unit: services.running[unit].identity for unit in daemon_identities
    } == daemon_identities
    final = states.read(
        expected_owner_uid=owner_uid,
        expected_owner_gid=owner_gid,
    )
    assert final.components["encode-runtime"].active == IDENTITY


def test_bulk_upgrade_materializes_after_stop_and_restarts_docker_before_writers(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    monkeypatch.setattr(database_module, "_ROOT_UID", owner_uid)
    monkeypatch.setattr(database_module, "_ROOT_GID", owner_gid)
    states, _prior = _state_with_active_components(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        staged_component="bulk-rnaseq-runtime",
        staged_identity=THIRD_IDENTITY,
    )
    _write_operator_test_database(
        layout,
        heads=("schema-v1",),
        value="preserved",
    )
    services = _TrackingServices({}, [])
    for unit, deployment in {
        "helixweave-redis.service": OLD_PLATFORM_IDENTITY,
        "helixweave-docker-rootless.service": OLD_BULK_IDENTITY,
        "helixweave-api.service": OLD_PLATFORM_IDENTITY,
        "helixweave-worker.service": OLD_PLATFORM_IDENTITY,
    }.items():
        service_request = OperatorRequest(
            operation="start",
            unit=unit,
            deployment_identity=deployment,
            task_identity=f"task-{'6' * 32}",
        )
        services.running[unit] = services._identity(service_request)
    bulk_requests: list[object] = []
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=services,
        action_runner=_FreshActionRunner(target_schema="schema-v1"),
        database_preparer=_FreshDatabasePreparer(layout, owner_uid, owner_gid),
        encode_runtime_preparer=_FreshRuntimePreparer(),
        bulk_runtime_preparer=_FreshBulkRuntimePreparer(bulk_requests),
        configuration=_FreshConfiguration([], layout),
        root_uid=owner_uid,
        root_gid=owner_gid,
        service_uid=owner_uid,
        service_gid=owner_gid,
    )
    journals = OperatorJournalStore(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        recovery_controller=controller,
    )

    with journals.operation(
        operation="activate",
        task_identity=TASK_IDENTITY,
        deployment_identity=THIRD_IDENTITY,
        component="bulk-rnaseq-runtime",
        unit=None,
    ) as journal:
        outcome = controller.execute(
            OperatorRequest(
                operation="activate",
                component="bulk-rnaseq-runtime",
                deployment_identity=THIRD_IDENTITY,
                task_identity=TASK_IDENTITY,
            ),
            journal=journal,
        )
        journal.complete()

    assert outcome.state == "activated"
    assert len(bulk_requests) == 1
    assert bulk_requests[0].operation == "activate"
    assert bulk_requests[0].candidate_bulk_identity == THIRD_IDENTITY
    assert [
        (action, unit)
        for action, unit, _deployment in services.calls
        if action in {"stop", "start"}
    ] == [
        ("stop", "helixweave-api.service"),
        ("stop", "helixweave-worker.service"),
        ("stop", "helixweave-docker-rootless.service"),
        ("start", "helixweave-docker-rootless.service"),
        ("start", "helixweave-api.service"),
        ("start", "helixweave-worker.service"),
    ]
    assert not any(unit == "helixweave-redis.service" for _, unit, _ in services.calls)
    assert (
        states.read(
            expected_owner_uid=owner_uid,
            expected_owner_gid=owner_gid,
        )
        .components["bulk-rnaseq-runtime"]
        .active
        == THIRD_IDENTITY
    )


def test_bulk_start_side_effect_is_synchronously_stopped_if_journal_update_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    monkeypatch.setattr(database_module, "_ROOT_UID", owner_uid)
    monkeypatch.setattr(database_module, "_ROOT_GID", owner_gid)
    states, prior = _state_with_active_components(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        staged_component="bulk-rnaseq-runtime",
        staged_identity=THIRD_IDENTITY,
    )
    _write_operator_test_database(
        layout,
        heads=("schema-v1",),
        value="preserved",
    )

    class StartThenFailServices(_TrackingServices):
        def start(self, request):
            service = super().start(request)
            if request.unit == "helixweave-docker-rootless.service":
                raise fail(
                    "TEST_AFTER_DOCKER_START_FAILED",
                    "Test Docker start failed.",
                    recoverable=True,
                )
            return service

    services = StartThenFailServices({}, [])
    for unit in ("helixweave-api.service", "helixweave-worker.service"):
        service_request = OperatorRequest(
            operation="start",
            unit=unit,
            deployment_identity=OLD_PLATFORM_IDENTITY,
            task_identity=f"task-{'4' * 32}",
        )
        services.running[unit] = services._identity(service_request)
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=services,
        action_runner=_FreshActionRunner(target_schema="schema-v1"),
        database_preparer=_FreshDatabasePreparer(layout, owner_uid, owner_gid),
        encode_runtime_preparer=_FreshRuntimePreparer(),
        bulk_runtime_preparer=_FreshBulkRuntimePreparer(),
        configuration=_FreshConfiguration([], layout),
        root_uid=owner_uid,
        root_gid=owner_gid,
        service_uid=owner_uid,
        service_gid=owner_gid,
    )
    journals = OperatorJournalStore(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        recovery_controller=controller,
    )

    with pytest.raises(DeploymentError) as captured:
        with journals.operation(
            operation="activate",
            task_identity=TASK_IDENTITY,
            deployment_identity=THIRD_IDENTITY,
            component="bulk-rnaseq-runtime",
            unit=None,
        ) as journal:
            controller.execute(
                OperatorRequest(
                    operation="activate",
                    component="bulk-rnaseq-runtime",
                    deployment_identity=THIRD_IDENTITY,
                    task_identity=TASK_IDENTITY,
                ),
                journal=journal,
            )

    assert captured.value.issue.code == "TEST_AFTER_DOCKER_START_FAILED"
    assert "helixweave-docker-rootless.service" not in services.running
    assert {
        unit: service.deployment_identity for unit, service in services.running.items()
    } == {
        "helixweave-api.service": OLD_PLATFORM_IDENTITY,
        "helixweave-worker.service": OLD_PLATFORM_IDENTITY,
    }
    assert (
        states.read(
            expected_owner_uid=owner_uid,
            expected_owner_gid=owner_gid,
        ).identity
        == prior.identity
    )
    history = json.loads(
        (layout.operator_transaction_history / f"{TASK_IDENTITY}.json").read_text()
    )
    assert history["phase"] == "aborted"
    assert history["failure_phase"] == "writers-stopped"


def test_incompatible_partial_assembly_is_rejected_without_side_effects(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    states = _supported_state_store(layout, owner_uid, owner_gid)
    with states.transaction(
        exclusive=True,
        expected_owner_uid=owner_uid,
        expected_owner_gid=owner_gid,
    ) as transaction:
        initial = transaction.initialize()
        prior = initial.stage("platform", IDENTITY)
        transaction.commit(
            prior,
            operation="stage-platform",
            expected_current_identity=initial.identity,
        )
    services = _TrackingServices({}, [])
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=services,
        action_runner=_FreshActionRunner(compatibility="incompatible"),
        database_preparer=_FreshDatabasePreparer(layout, owner_uid, owner_gid),
        encode_runtime_preparer=_FreshRuntimePreparer(),
        configuration=_FreshConfiguration([], layout),
        root_uid=owner_uid,
        root_gid=owner_gid,
        service_uid=owner_uid,
        service_gid=owner_gid,
    )
    journals = OperatorJournalStore(layout, owner_uid=owner_uid, owner_gid=owner_gid)

    with pytest.raises(DeploymentError) as caught:
        with journals.operation(
            operation="activate",
            task_identity=TASK_IDENTITY,
            deployment_identity=IDENTITY,
            component="platform",
            unit=None,
        ) as journal:
            controller.execute(
                OperatorRequest(
                    operation="activate",
                    component="platform",
                    deployment_identity=IDENTITY,
                    task_identity=TASK_IDENTITY,
                ),
                journal=journal,
            )

    assert caught.value.issue.code == "DEPLOYMENT_COMPATIBILITY_FAILED"
    assert services.calls == []
    assert not layout.database.exists()
    assert (
        states.read(
            expected_owner_uid=owner_uid,
            expected_owner_gid=owner_gid,
        ).identity
        == prior.identity
    )


def test_root_verify_replaces_candidate_unavailable_database_and_configuration(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    monkeypatch.setattr(database_module, "_ROOT_UID", owner_uid)
    monkeypatch.setattr(database_module, "_ROOT_GID", owner_gid)
    states, state = _state_with_active_components(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
    )
    _write_operator_test_database(
        layout,
        heads=("schema-v1",),
        value="verified",
    )
    monkeypatch.setattr(
        BundleStore,
        "verify_installed",
        lambda _store, _component, identity, **_expected: SimpleNamespace(
            identity=identity
        ),
    )
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=_TrackingServices({}, []),
        action_runner=_FreshActionRunner(target_schema="schema-v1"),
        database_preparer=_FreshDatabasePreparer(layout, owner_uid, owner_gid),
        encode_runtime_preparer=_FreshRuntimePreparer(),
        configuration=_FreshConfiguration([], layout),
        root_uid=owner_uid,
        root_gid=owner_gid,
        service_uid=owner_uid,
        service_gid=owner_gid,
    )

    receipt = controller.verify(
        OperatorRequest(
            operation="verify",
            deployment_identity=state.identity,
            task_identity=TASK_IDENTITY,
        )
    )

    assert receipt.status == "observed"
    assert receipt.compatibility == "compatible"
    assert receipt.database_after_identity == database_content_identity(
        inspect_database(
            layout.database,
            expected_owner_uid=owner_uid,
            expected_owner_gid=owner_gid,
        )
    )
    assert (
        receipt.accepted_schema_heads == receipt.target_schema_heads == ("schema-v1",)
    )
    assert receipt.migration_required is False
    assert receipt.rollback_supported is True
    assert receipt.readiness["database-schema"].status == "ready"
    assert receipt.readiness["configuration"].status == "ready"
    assert receipt.readiness["permissions"].status == "ready"
    assert receipt.readiness["redis"].reason_code == "REDIS_UNAVAILABLE"
    assert receipt.readiness["docker"].reason_code == "DOCKER_UNAVAILABLE"


def test_root_verify_holds_state_lock_through_bulk_runtime_observation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    monkeypatch.setattr(database_module, "_ROOT_UID", owner_uid)
    monkeypatch.setattr(database_module, "_ROOT_GID", owner_gid)
    states, state = _state_with_active_components(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
    )
    _write_operator_test_database(
        layout,
        heads=("schema-v1",),
        value="verified",
    )
    monkeypatch.setattr(
        BundleStore,
        "verify_installed",
        lambda _store, _component, identity, **_expected: SimpleNamespace(
            identity=identity
        ),
    )
    lock_failures: list[str] = []

    def assert_shared_lock_is_held() -> None:
        with pytest.raises(DeploymentError) as captured:
            with states.transaction(
                exclusive=True,
                expected_owner_uid=owner_uid,
                expected_owner_gid=owner_gid,
            ):
                pytest.fail("exclusive state lock unexpectedly acquired")
        lock_failures.append(captured.value.issue.code)

    class LockAssertingActionRunner:
        def run(self, request):
            assert_shared_lock_is_held()
            return _FreshActionRunner(target_schema="schema-v1").run(request)

    class LockAssertingBulkRuntimePreparer(_FreshBulkRuntimePreparer):
        def prepare(self, request):
            assert_shared_lock_is_held()
            return super().prepare(request)

    services = _TrackingServices({}, [])
    docker_request = OperatorRequest(
        operation="start",
        unit="helixweave-docker-rootless.service",
        deployment_identity=OLD_BULK_IDENTITY,
        task_identity=TASK_IDENTITY,
    )
    services.running[docker_request.unit] = services._identity(docker_request)
    bulk_requests: list[object] = []
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=services,
        action_runner=LockAssertingActionRunner(),
        database_preparer=_FreshDatabasePreparer(layout, owner_uid, owner_gid),
        encode_runtime_preparer=_FreshRuntimePreparer(),
        bulk_runtime_preparer=LockAssertingBulkRuntimePreparer(bulk_requests),
        configuration=_FreshConfiguration([], layout),
        root_uid=owner_uid,
        root_gid=owner_gid,
        service_uid=owner_uid,
        service_gid=owner_gid,
    )

    receipt = controller.verify(
        OperatorRequest(
            operation="verify",
            deployment_identity=state.identity,
            task_identity=TASK_IDENTITY,
        )
    )

    assert lock_failures == ["DEPLOYMENT_BUSY", "DEPLOYMENT_BUSY"]
    assert len(bulk_requests) == 1
    assert bulk_requests[0].operation == "verify"
    assert bulk_requests[0].candidate_bulk_identity == OLD_BULK_IDENTITY
    assert receipt.readiness["docker"].status == "ready"
    assert not any(action in {"start", "stop"} for action, _, _ in services.calls)


def test_root_reference_probe_consumes_only_bounded_public_api_documents() -> None:
    workflow_ids = (
        "bulk-rnaseq",
        "encode-style-chipseq-cuttag-atac-mnase",
    )
    documents: dict[str, object] = {
        "/api/v1/workflows/": {
            "ok": True,
            "workflows": [
                {"metadata": {"workflow_id": workflow_id}}
                for workflow_id in workflow_ids
            ],
            "issues": [],
        }
    }
    for index, workflow_id in enumerate(workflow_ids, start=1):
        documents[f"/api/v1/workflows/{workflow_id}/reference-profiles"] = {
            "ok": True,
            "workflow_id": workflow_id,
            "profiles": [
                {
                    "profile_id": f"refp_{index:032x}",
                    "revision_id": f"refpr_{index:032x}",
                    "revision_number": index,
                    "display_name": f"public reference {index}",
                    "organism": "Mus musculus",
                    "assembly": "mm10",
                    "identity_sha256": f"{index:064x}",
                }
            ],
            "issues": [],
        }

    observed = _reference_api_readiness(fetcher=documents.__getitem__)

    assert observed.status == "ready"
    assert observed.reason_code == "READY"
    assert observed.identity is not None
    assert "/" not in observed.identity
    documents["/api/v1/workflows/bulk-rnaseq/reference-profiles"] = {
        "private_path": "/references/private"
    }
    assert _reference_api_readiness(fetcher=documents.__getitem__) == ReadinessCheck(
        "not-ready", "REFERENCE_NOT_READY"
    )


def test_root_schema_plan_allows_only_fresh_platform_creation() -> None:
    assert _plan_schema_transition(
        operation="activate",
        component="platform",
        database_exists=False,
        current_schema_heads=(),
        target_schema_heads=("schema-v1",),
        prior_platform_schema_heads=None,
    )
    with pytest.raises(DeploymentError) as caught:
        _plan_schema_transition(
            operation="activate",
            component="encode-runtime",
            database_exists=False,
            current_schema_heads=(),
            target_schema_heads=("schema-v1",),
            prior_platform_schema_heads=None,
        )
    assert caught.value.issue.code == "DEPLOYMENT_SCHEMA_INCOMPATIBLE"


def test_root_rejects_candidate_database_observation_as_an_authority() -> None:
    readiness = {
        name: ReadinessCheck("ready", "READY", IDENTITY) for name in VERIFICATION_CHECKS
    }
    receipt = DeploymentActionReceipt.create(
        request_identity=IDENTITY,
        status="admitted",
        compatibility="compatible",
        database_before_identity=IDENTITY,
        accepted_schema_heads=("schema-v1",),
        target_schema_heads=("schema-v1",),
        migration_inventory_identity=IDENTITY,
        known_schema_revisions=(),
        migration_required=False,
        rollback_supported=True,
        api_contract_sha256="a" * 64,
        native_identities={
            component: IDENTITY
            for component in (
                "platform",
                "encode-runtime",
                "bulk-rnaseq-runtime",
            )
        },
        frontend_identity=IDENTITY,
        reference_compatibility_identity=IDENTITY,
        readiness=readiness,
    )

    with pytest.raises(DeploymentError) as caught:
        _candidate_schema_target(receipt)

    assert caught.value.issue.code == "OPERATOR_ACTION_RECEIPT_INVALID"


def test_root_schema_plan_never_migrates_for_runtime_only_activation() -> None:
    assert not _plan_schema_transition(
        operation="activate",
        component="bulk-rnaseq-runtime",
        database_exists=True,
        current_schema_heads=("schema-v1",),
        target_schema_heads=("schema-v1",),
        prior_platform_schema_heads=None,
    )
    with pytest.raises(DeploymentError) as caught:
        _plan_schema_transition(
            operation="activate",
            component="bulk-rnaseq-runtime",
            database_exists=True,
            current_schema_heads=("unknown-head",),
            target_schema_heads=("schema-v1",),
            prior_platform_schema_heads=None,
        )
    assert caught.value.issue.code == "DEPLOYMENT_SCHEMA_INCOMPATIBLE"


def test_root_schema_plan_migrates_only_from_the_active_platform_head() -> None:
    assert _plan_schema_transition(
        operation="activate",
        component="platform",
        database_exists=True,
        current_schema_heads=("schema-v1",),
        target_schema_heads=("schema-v2",),
        prior_platform_schema_heads=("schema-v1",),
    )
    with pytest.raises(DeploymentError) as caught:
        _plan_schema_transition(
            operation="activate",
            component="platform",
            database_exists=True,
            current_schema_heads=("unknown-head",),
            target_schema_heads=("schema-v2",),
            prior_platform_schema_heads=("schema-v1",),
        )
    assert caught.value.issue.code == "DEPLOYMENT_SCHEMA_INCOMPATIBLE"


def test_root_schema_plan_refuses_schema_changing_rollback() -> None:
    with pytest.raises(DeploymentError) as caught:
        _plan_schema_transition(
            operation="rollback",
            component="platform",
            database_exists=True,
            current_schema_heads=("schema-v2",),
            target_schema_heads=("schema-v1",),
            prior_platform_schema_heads=None,
        )
    assert caught.value.issue.code == "DEPLOYMENT_SCHEMA_INCOMPATIBLE"


def test_templates_encode_one_bounded_hybrid_topology() -> None:
    api = (TEMPLATES / "helixweave-api.service.in").read_text()
    worker = (TEMPLATES / "helixweave-worker.service.in").read_text()
    redis = (TEMPLATES / "helixweave-redis.service").read_text()
    docker = (TEMPLATES / "helixweave-docker-rootless.service").read_text()
    target = (TEMPLATES / "helixweave.target").read_text()

    for unit in (api, worker, redis, docker):
        assert "ExecStart=/bin/sh" not in unit
        assert "ExecStart=/bin/bash" not in unit
        assert "docker compose" not in unit.lower()
    assert "ExecStart=/usr/libexec/helixweave-active-service api" in api
    assert "ExecStart=/usr/libexec/helixweave-active-service worker" in worker
    assert "User=helixweave-api\nGroup=helixweave\n" in api
    assert "SupplementaryGroups=helixweave-api" in api
    assert "User=helixweave\nGroup=helixweave\n" in worker
    assert "UMask=0007" in api
    assert "UMask=0007" in worker
    assert "static frontend" in api
    assert "npm" not in api
    assert "vite" not in api.lower()
    assert "EnvironmentFile=-/etc/helixweave/secrets.env" in api
    assert "EnvironmentFile=/etc/helixweave/platform.env" not in api
    assert "EnvironmentFile=/etc/helixweave/platform.env" not in worker
    assert "InaccessiblePaths=/var/lib/helixweave/operator " in api
    assert (
        "InaccessiblePaths=/etc/helixweave/secrets.env /var/lib/helixweave/operator "
        in worker
    )
    active_service = (TEMPLATES / "helixweave-active-service").read_text()
    assert (
        'STATE_CURRENT = Path("/var/lib/helixweave/deployment/current")'
        in active_service
    )
    assert "InaccessiblePaths=/etc/helixweave/secrets.env" in worker
    assert "/run/helixweave/docker" in api
    assert "-/var/run/docker.sock" in api
    assert "port 0" in (TEMPLATES / "redis.conf").read_text()
    assert "dockerd-rootless.sh" in docker
    assert "Type=simple" in docker
    assert "Delegate=yes" in docker
    assert "unix:///run/helixweave/docker/docker.sock" in docker
    assert "--data-root=/var/lib/helixweave/docker-rootless" in docker
    docker_exec_start = next(
        line for line in docker.splitlines() if line.startswith("ExecStart=")
    )
    assert docker_exec_start == (
        "ExecStart=/usr/bin/dockerd-rootless.sh "
        "--host=unix:///run/helixweave/docker/docker.sock "
        "--data-root=/var/lib/helixweave/docker-rootless "
        "--exec-root=/run/helixweave/docker/exec "
        "--pidfile=/run/helixweave/docker/docker.pid "
        "--bridge=none --iptables=false --ip-forward=false --ip-masq=false "
        "--group=0"
    )
    assert docker.count("--group=0") == 1
    assert "/var/run/docker.sock" not in docker
    assert (
        "Requires=helixweave-redis.service helixweave-docker-rootless.service" in target
    )
    assert not any(
        path.name == "helixweave-frontend.service" for path in TEMPLATES.iterdir()
    )


def test_active_service_accepts_only_the_generation_bound_environment() -> None:
    namespace = runpy.run_path(str(TEMPLATES / "helixweave-active-service"))
    state = DeploymentState.initial().stage("platform", IDENTITY).activate("platform")
    active = {
        component: state.components[component].active for component in state.components
    }
    environment = render_platform_environment(
        DeploymentLayout.supported(),
        state,
        api_contract_sha256="a" * 64,
    )

    assert (
        namespace["_parse_environment"](
            environment.content,
            state.identity,
            active,
        )["HELIXWEAVE_DEPLOYMENT_IDENTITY"]
        == state.identity
    )
    with pytest.raises(ValueError):
        namespace["_parse_environment"](
            environment.content.replace(
                b"ENCODE_PIPELINE_QUEUE_NAME=encode-runs\n",
                b"UNTRUSTED=value\n",
            ),
            state.identity,
            active,
        )
    with pytest.raises(ValueError):
        namespace["_parse_environment"](
            environment.content.replace(IDENTITY.encode(), SERVICE_IDENTITY.encode()),
            state.identity,
            active,
        )


@pytest.mark.parametrize(
    "error",
    (
        OSError(errno.ENOENT, "missing", "/private/release/bin/service"),
        OSError(errno.EACCES, "denied", "/private/release/bin/service"),
        OSError(errno.ENOMEM, "token=secret"),
    ),
)
def test_active_service_redacts_execve_errors(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    error: OSError,
) -> None:
    namespace = runpy.run_path(str(TEMPLATES / "helixweave-active-service"))
    main = namespace["main"]
    binary = Path("/private/release/bin/helixweave-service")

    def fail_execve(*_args) -> None:
        raise error

    with monkeypatch.context() as scoped:
        scoped.setitem(
            main.__globals__,
            "_active_generation",
            lambda: (IDENTITY, {"PRIVATE_COORDINATE": "/private/value"}),
        )
        scoped.setitem(main.__globals__, "_launcher", lambda _identity: binary)
        scoped.setattr(main.__globals__["os"], "environ", {"PRIVATE_TOKEN": "secret"})
        scoped.setattr(main.__globals__["os"], "execve", fail_execve)
        observed = main(("api",))

    captured = capsys.readouterr()
    assert observed == 69
    assert captured.out == ""
    assert captured.err == ""


def test_active_service_retains_abnormal_execve_return_contract(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    namespace = runpy.run_path(str(TEMPLATES / "helixweave-active-service"))
    main = namespace["main"]
    binary = Path("/opt/helixweave/releases/platform/identity/bin/service")
    calls: list[tuple[object, ...]] = []

    with monkeypatch.context() as scoped:
        scoped.setitem(
            main.__globals__,
            "_active_generation",
            lambda: (IDENTITY, {"HELIXWEAVE_DEPLOYMENT_IDENTITY": IDENTITY}),
        )
        scoped.setitem(main.__globals__, "_launcher", lambda _identity: binary)
        scoped.setattr(main.__globals__["os"], "environ", {})
        scoped.setattr(
            main.__globals__["os"],
            "execve",
            lambda *arguments: calls.append(arguments),
        )

        observed = main(("worker",))

    assert observed == 70
    assert calls == [
        (
            binary,
            (str(binary), "worker"),
            {
                "HELIXWEAVE_DEPLOYMENT_IDENTITY": IDENTITY,
                "PYTHONNOUSERSITE": "1",
            },
        )
    ]


def test_sudoers_and_helper_templates_do_not_expose_a_shell_or_environment() -> None:
    policy = (TEMPLATES / "helixweave-operator.sudoers").read_text()
    helper = (TEMPLATES / "helixweave-operator").read_text()
    gate_cleanup = (TEMPLATES / "helixweave-gate-cleanup").read_text()
    gate_launcher = (TEMPLATES / "helixweave-gate-cleanup-launcher").read_text()
    launcher = (TEMPLATES / "helixweave-operator-launcher").read_text()

    assert policy.count("/usr/libexec/helixweave-operator") == 10
    assert "NOPASSWD: HELIXWEAVE_OPERATOR" in policy
    assert "!setenv" in policy
    assert "/bin/sh" not in policy
    assert "/bin/bash" not in policy
    assert "python" not in policy.lower()
    assert policy.count("/usr/libexec/helixweave-gate-cleanup") == 1
    assert (
        "cleanup --task-id * --plan-id * --executor-id * --closure-id * "
        "--backend-id *" in policy
    )
    assert helper.startswith("#!/usr/bin/python3 -I\n")
    assert "subprocess" not in helper
    assert "sys.path.insert(0, str(library))" in helper
    assert "from encode_pipeline.deployment.operator import main" in helper
    assert "os.execve(" in launcher
    assert "/opt/helixweave/operator/current/bin/helixweave-operator" in launcher
    assert "encode_pipeline" not in launcher
    assert gate_cleanup.startswith("#!/usr/bin/python3 -I\n")
    assert "subprocess" not in gate_cleanup
    assert "encode_pipeline" not in gate_cleanup
    assert "os.execve(" in gate_launcher
    assert 'BACKEND_ROOT = "/opt/helixweave/operator"' in gate_launcher
    assert 'sys.argv[8] != "--closure-id"' in gate_launcher
    assert "/current/" not in gate_launcher


def test_gate_cleanup_helper_accepts_only_the_fixed_identity_grammar() -> None:
    namespace = runpy.run_path(str(TEMPLATES / "helixweave-gate-cleanup"))
    valid = namespace["_valid"]

    assert valid(
        (
            "cleanup",
            "--task-id",
            TASK_IDENTITY,
            "--plan-id",
            IDENTITY,
            "--executor-id",
            SERVICE_IDENTITY,
            "--closure-id",
            IDENTITY,
            "--backend-id",
            SERVICE_IDENTITY,
        )
    )
    assert not valid(
        (
            "cleanup",
            "--task-id",
            TASK_IDENTITY,
            "--plan-id",
            IDENTITY,
            "--executor-id",
            SERVICE_IDENTITY,
            "--closure-id",
            IDENTITY,
            "--backend-id",
            SERVICE_IDENTITY,
            "--command=/bin/sh",
        )
    )


def test_gate_cleanup_launcher_execs_the_plan_bound_closure_not_current(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[object, ...]] = []
    values = (
        "cleanup",
        "--task-id",
        TASK_IDENTITY,
        "--plan-id",
        IDENTITY,
        "--executor-id",
        SERVICE_IDENTITY,
        "--closure-id",
        IDENTITY,
        "--backend-id",
        SERVICE_IDENTITY,
    )
    monkeypatch.setattr(sys, "argv", ["helixweave-gate-cleanup", *values])
    monkeypatch.setattr(os, "execve", lambda *arguments: calls.append(arguments))

    runpy.run_path(str(TEMPLATES / "helixweave-gate-cleanup-launcher"))

    assert len(calls) == 1
    _python, argv, environment = calls[0]
    assert argv[5] == (
        f"/opt/helixweave/operator/{IDENTITY}/bin/helixweave-gate-cleanup"
    )
    assert "/current/" not in argv[5]
    assert tuple(argv[6:]) == values
    assert set(environment) == {"LANG", "LC_ALL", "PATH"}


def test_deployed_operator_helper_loads_only_the_root_owned_closure(
    tmp_path: Path,
) -> None:
    root = tmp_path / "host"
    root.mkdir()
    backend = bootstrap.HostBootstrapBackend(
        source_root=TEMPLATES,
        root_prefix=root,
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
        command_runner=lambda command: pytest.fail(f"unexpected command: {command}"),
        sudoers_validator=lambda path: True,
    )
    result = backend.apply(operation="install", invoking_user="labadmin")
    script = (
        root
        / "opt/helixweave/operator"
        / result.closure_identity
        / "bin/helixweave-operator"
    )
    completed = subprocess.run(
        (sys.executable, "-I", "-S", "-E", "-B", str(script), "status"),
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 77
    receipt = json.loads(completed.stderr)
    assert receipt["issue"]["code"] == "OPERATOR_PRIVILEGE_REQUIRED"


def test_tmpfiles_separates_immutable_data_and_external_references() -> None:
    content = (TEMPLATES / "helixweave.tmpfiles.conf").read_text()

    assert "/var/lib/helixweave 0755 root root" in content
    assert "/run/helixweave 0755 root root" in content
    assert "/opt/helixweave/releases/platform 0555 root root" in content
    assert "/opt/helixweave/runtimes/bulk-rnaseq 0555 root root" in content
    assert "/etc/helixweave/secrets.env 0640 root helixweave-api" in content
    assert "z /etc/helixweave/secrets.env 0640 root helixweave-api" in content
    assert "/var/lib/helixweave/database/live 2770 helixweave helixweave" in content
    assert (
        "z /var/lib/helixweave/database/live/platform.db 0660 "
        "helixweave helixweave" in content
    )
    assert "/var/lib/helixweave/operator/database-backups 0700 root root" in content
    assert "/var/lib/helixweave/deployment 0755 root root" in content
    assert "/var/lib/helixweave/deployment/generations 0755 root root" in content
    assert "/var/lib/helixweave/operator 0710 root helixweave-operators" in content
    assert (
        "/var/lib/helixweave/operator/state 0750 root helixweave-operators" in content
    )
    assert "/operator/transactions/history 0700 root root" in content
    assert "/var/lib/helixweave/workspaces 2770 helixweave helixweave" in content
    assert "/var/lib/helixweave/artifacts 2770 helixweave helixweave" in content
    assert "/operator/action 0750 root helixweave-candidate" in content
    assert "/operator/encode-runtime 0750 root helixweave-candidate" in content
    for component in ("platform", "encode-runtime", "bulk-rnaseq-runtime"):
        assert (
            f"/operator/ingress/{component} 2770 root helixweave-operators" in content
        )
    assert " 2730 root helixweave-operators" not in content
    assert "reference" not in content.lower()

    for unit_name in ("helixweave-api.service.in", "helixweave-worker.service.in"):
        unit = (TEMPLATES / unit_name).read_text()
        assert "ReadWritePaths=/var/lib/helixweave /run/helixweave" not in unit
        assert "/var/lib/helixweave/operator" in unit


def test_candidate_units_have_a_distinct_identity_and_no_secret_or_runtime_access() -> (
    None
):
    action = (TEMPLATES / "helixweave-operator-action.service").read_text()
    encode = (TEMPLATES / "helixweave-encode-runtime-prepare.service").read_text()
    database = (TEMPLATES / "helixweave-db-prepare.service.in").read_text()

    for unit in (action, encode):
        assert "User=helixweave-candidate\nGroup=helixweave-candidate\n" in unit
        assert "/etc/helixweave" in unit
        assert "/run/helixweave/redis" in unit
        assert "/run/helixweave/docker" in unit
        assert "-/var/run/docker.sock" in unit
        assert "/var/lib/helixweave/workspaces" in unit
        assert "/var/lib/helixweave/artifacts" in unit
    assert "SupplementaryGroups=" not in action
    assert "SupplementaryGroups=" not in encode
    assert "InaccessiblePaths=/etc/helixweave " in action
    assert "/var/lib/helixweave/database" in action
    assert "User=helixweave\nGroup=helixweave\n" in database
    assert "UMask=0007" in database
    assert "InaccessiblePaths=/etc/helixweave/secrets.env" in database


@dataclass
class RecordingBootstrapBackend:
    installed: int = 11
    closure_identity: str = IDENTITY
    operation: str | None = None
    invoking_user: str | None = None

    def apply(self, *, operation: str, invoking_user: str):
        self.operation = operation
        self.invoking_user = invoking_user
        return bootstrap.BootstrapResult(self.installed, self.closure_identity)


def test_bootstrap_requires_root_local_tty_and_verified_sudo_identity(
    monkeypatch,
) -> None:
    backend = RecordingBootstrapBackend()
    monkeypatch.setattr(
        bootstrap.pwd,
        "getpwnam",
        lambda name: SimpleNamespace(pw_uid=1001, pw_name=name),
    )

    receipt = bootstrap.run_bootstrap(
        ("install",),
        backend=backend,
        effective_uid=0,
        stdin_is_tty=True,
        stderr_is_tty=True,
        environ={"SUDO_USER": "labadmin", "SUDO_UID": "1001"},
    )

    assert receipt == {
        "schema_version": "helixweave-operator-bootstrap-receipt-v1",
        "operation": "install",
        "status": "complete",
        "installed_files": 11,
        "closure_identity": IDENTITY,
    }
    assert backend.operation == "install"
    assert backend.invoking_user == "labadmin"


@pytest.mark.parametrize(
    ("argv", "uid", "stdin_tty", "stderr_tty", "environment", "code"),
    [
        (("install", "/tmp/root"), 0, True, True, {}, "BOOTSTRAP_REQUEST_INVALID"),
        (("install",), 1001, True, True, {}, "BOOTSTRAP_PRIVILEGE_REQUIRED"),
        (("install",), 0, False, True, {}, "BOOTSTRAP_INTERACTIVE_REQUIRED"),
        (("update",), 0, True, False, {}, "BOOTSTRAP_INTERACTIVE_REQUIRED"),
        (
            ("install",),
            0,
            True,
            True,
            {"SUDO_USER": "root", "SUDO_UID": "0"},
            "BOOTSTRAP_CALLER_INVALID",
        ),
    ],
)
def test_bootstrap_rejects_noninteractive_or_caller_controlled_scope(
    argv,
    uid,
    stdin_tty,
    stderr_tty,
    environment,
    code,
) -> None:
    backend = RecordingBootstrapBackend()

    with pytest.raises(bootstrap.BootstrapFailure) as caught:
        bootstrap.run_bootstrap(
            argv,
            backend=backend,
            effective_uid=uid,
            stdin_is_tty=stdin_tty,
            stderr_is_tty=stderr_tty,
            environ=environment,
        )

    assert caught.value.code == code
    assert backend.operation is None


def test_bootstrap_installs_only_fixed_files_with_root_equivalent_modes(
    tmp_path: Path,
) -> None:
    validator_calls: list[Path] = []

    def validate(path: Path) -> bool:
        validator_calls.append(path)
        return path.read_text().startswith("User_Alias HELIXWEAVE_OPERATORS")

    backend = bootstrap.HostBootstrapBackend(
        source_root=TEMPLATES,
        root_prefix=tmp_path,
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
        command_runner=lambda command: pytest.fail(f"unexpected command: {command}"),
        sudoers_validator=validate,
    )

    installed = backend.apply(operation="install", invoking_user="labadmin")

    assert installed.installed_files == (
        len(bootstrap.CLOSURE_SPECS) + len(bootstrap.BOUNDARY_SPECS) + 1
    )
    helper = tmp_path / "usr" / "libexec" / "helixweave-operator"
    gate_cleanup = tmp_path / "usr" / "libexec" / "helixweave-gate-cleanup"
    sudoers = tmp_path / "etc" / "sudoers.d" / "helixweave-operator"
    operator_root = tmp_path / "opt" / "helixweave" / "operator"
    current = operator_root / "current"
    closure = operator_root / os.readlink(current)
    assert os.readlink(current) == installed.closure_identity
    assert stat.S_IMODE(helper.stat().st_mode) == 0o555
    assert stat.S_IMODE(gate_cleanup.stat().st_mode) == 0o555
    assert (
        gate_cleanup.read_text()
        == (TEMPLATES / "helixweave-gate-cleanup-launcher").read_text()
    )
    assert stat.S_IMODE(sudoers.stat().st_mode) == 0o440
    assert stat.S_IMODE(closure.stat().st_mode) == 0o555
    assert (
        stat.S_IMODE(
            (closure / "templates" / "helixweave-api.service.in").stat().st_mode
        )
        == 0o444
    )
    assert (closure / "bin" / "helixweave-operator").read_text() == (
        TEMPLATES / "helixweave-operator"
    ).read_text()
    assert (closure / "bin" / "helixweave-gate-cleanup").read_text() == (
        TEMPLATES / "helixweave-gate-cleanup"
    ).read_text()
    assert (closure / "boundary" / "helixweave-gate-cleanup-launcher").read_text() == (
        TEMPLATES / "helixweave-gate-cleanup-launcher"
    ).read_text()
    for spec in bootstrap.STABLE_BOUNDARY_SPECS:
        installed_path = tmp_path.joinpath(*spec.destination.parts[1:])
        snapshot = closure / "boundary" / spec.source_name
        assert installed_path.read_bytes() == snapshot.read_bytes()
    for spec in bootstrap.LINKED_BOUNDARY_SPECS:
        installed_path = tmp_path.joinpath(*spec.destination.parts[1:])
        assert installed_path.is_symlink()
        assert (
            os.readlink(installed_path)
            == (
                bootstrap.OPERATOR_ROOT
                / "current"
                / bootstrap.BOUNDARY_SNAPSHOT_PATHS[spec.source_name]
            ).as_posix()
        )
    assert (closure / "closure.json").is_file()
    assert backend.verify_selected_boundaries() == installed.closure_identity
    assert len(validator_calls) == 1
    assert not list(tmp_path.rglob("*.tmp"))
    assert not list(tmp_path.rglob(".partial-*"))


@dataclass(frozen=True)
class _StaticObservationProvider:
    observation: OperatorObservation

    def observe(self, _request: OperatorRequest) -> OperatorObservation:
        return self.observation


@dataclass(frozen=True)
class _StaticVerificationController:
    receipt: DeploymentActionReceipt

    def verify(self, _request: OperatorRequest) -> DeploymentActionReceipt:
        return self.receipt


def _static_verification_receipt() -> DeploymentActionReceipt:
    readiness = {
        name: ReadinessCheck("ready", "READY", IDENTITY) for name in VERIFICATION_CHECKS
    }
    return DeploymentActionReceipt.create(
        request_identity=IDENTITY,
        status="observed",
        compatibility="compatible",
        database_before_identity=IDENTITY,
        database_after_identity=IDENTITY,
        accepted_schema_heads=("schema-v1",),
        target_schema_heads=("schema-v1",),
        migration_inventory_identity=IDENTITY,
        known_schema_revisions=(),
        migration_required=False,
        rollback_supported=True,
        api_contract_sha256="a" * 64,
        native_identities={
            "platform": IDENTITY,
            "encode-runtime": SERVICE_IDENTITY,
            "bulk-rnaseq-runtime": THIRD_IDENTITY,
        },
        frontend_identity=IDENTITY,
        reference_compatibility_identity=IDENTITY,
        readiness=readiness,
    )


def test_host_observe_and_verify_recheck_all_selected_stable_boundary_bytes(
    tmp_path: Path,
) -> None:
    root = tmp_path / "host"
    root.mkdir()
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    bootstrap_backend = bootstrap.HostBootstrapBackend(
        source_root=TEMPLATES,
        root_prefix=root,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        command_runner=lambda command: pytest.fail(f"unexpected command: {command}"),
        sudoers_validator=lambda _path: True,
    )
    installed = bootstrap_backend.apply(operation="install", invoking_user="labadmin")
    observer = HostOperatorBoundaryObserver(
        operator_root=root / "opt/helixweave/operator",
        host_root=root,
        root_uid=owner_uid,
        root_gid=owner_gid,
    )
    healthy = observer.observe()
    assert healthy.status == "ready"
    assert healthy.identity == verify_stable_operator_boundary(
        operator_root=root / "opt/helixweave/operator",
        host_root=root,
        expected_uid=owner_uid,
        expected_gid=owner_gid,
    )
    assert installed.closure_identity == os.readlink(
        root / "opt/helixweave/operator/current"
    )

    observation = OperatorObservation.create(
        state_identity=IDENTITY,
        active={
            "platform": IDENTITY,
            "encode-runtime": SERVICE_IDENTITY,
            "bulk-rnaseq-runtime": THIRD_IDENTITY,
        },
        database_schema_identity=IDENTITY,
        database_schema_heads=("schema-v1",),
        services={unit: None for unit in SERVICE_UNITS},
    )
    backend = _host_backend(
        DeploymentLayout.isolated(tmp_path / "deployment"),
        observation_provider=_StaticObservationProvider(observation),
        verification_controller=_StaticVerificationController(
            _static_verification_receipt()
        ),
        boundary_observer=observer,
    )
    observed = backend.execute(
        parse_request(("observe", IDENTITY, TASK_IDENTITY)), bundle_path=None
    )
    verified = backend.execute(
        parse_request(("verify", IDENTITY, TASK_IDENTITY)), bundle_path=None
    )
    assert observed.observation is not None
    assert observed.observation.operator_boundary == healthy
    assert verified.verification is not None
    permission_evidence = verified.verification.readiness["permissions"]
    assert permission_evidence.status == "ready"
    assert permission_evidence.reason_code == "READY"
    assert permission_evidence.identity not in {None, healthy.identity, IDENTITY}

    for spec in bootstrap.STABLE_BOUNDARY_SPECS:
        installed_path = root.joinpath(*spec.destination.parts[1:])
        original = installed_path.read_bytes()
        installed_path.chmod(0o600)
        installed_path.write_bytes(original + b"\n# tampered\n")
        installed_path.chmod(spec.mode)
        try:
            with pytest.raises(StableBoundaryError):
                verify_stable_operator_boundary(
                    operator_root=root / "opt/helixweave/operator",
                    host_root=root,
                    expected_uid=owner_uid,
                    expected_gid=owner_gid,
                )
            observed = backend.execute(
                parse_request(("observe", IDENTITY, TASK_IDENTITY)), bundle_path=None
            )
            verified = backend.execute(
                parse_request(("verify", IDENTITY, TASK_IDENTITY)), bundle_path=None
            )
            assert observed.observation is not None
            assert observed.observation.operator_boundary == (
                OperatorBoundaryObservation.invalid()
            )
            assert verified.verification is not None
            assert verified.verification.readiness["permissions"] == ReadinessCheck(
                "not-ready", "PERMISSION_INVALID"
            )
            rendered = json.dumps(
                {
                    "observation": observed.observation.to_dict(),
                    "verification": verified.verification.to_dict(),
                }
            )
            assert "/usr/" not in rendered
            assert "/opt/" not in rendered
        finally:
            installed_path.chmod(0o600)
            installed_path.write_bytes(original)
            installed_path.chmod(spec.mode)


def test_bootstrap_refuses_symlinked_destination_ancestors(tmp_path: Path) -> None:
    outside = tmp_path / "outside"
    outside.mkdir()
    (tmp_path / "usr").mkdir()
    (tmp_path / "usr" / "libexec").symlink_to(outside, target_is_directory=True)
    backend = bootstrap.HostBootstrapBackend(
        source_root=TEMPLATES,
        root_prefix=tmp_path,
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
        sudoers_validator=lambda path: True,
    )

    with pytest.raises(bootstrap.BootstrapFailure) as caught:
        backend.apply(operation="install", invoking_user="labadmin")

    assert caught.value.code == "BOOTSTRAP_DESTINATION_UNSAFE"
    assert list(outside.iterdir()) == []


def test_bootstrap_source_contains_no_password_or_sudo_invocation() -> None:
    source = BOOTSTRAP_SCRIPT.read_text()

    assert "getpass" not in source
    assert "input(" not in source
    assert 'Path("/usr/bin/sudo")' not in source
    assert 'subprocess.run(("sudo"' not in source
    assert "shell=True" not in source
    assert "stdin=subprocess.DEVNULL" in source


def test_bootstrap_creates_three_fixed_distinct_runtime_accounts(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    groups: dict[str, SimpleNamespace] = {}
    accounts: dict[str, SimpleNamespace] = {}
    commands: list[tuple[str, ...]] = []

    def get_group(name: str):
        if name not in groups:
            raise KeyError(name)
        return groups[name]

    def get_account(name: str):
        if name not in accounts:
            raise KeyError(name)
        return accounts[name]

    def run(command) -> bool:
        values = tuple(command)
        commands.append(values)
        if values[:2] == (str(bootstrap.GROUPADD), "--system"):
            name = values[2]
            groups[name] = SimpleNamespace(
                gr_name=name,
                gr_gid=3000 + len(groups),
                gr_mem=(),
            )
            return True
        if values[:2] == (str(bootstrap.USERADD), "--system"):
            name = values[-1]
            group = groups[values[3]]
            accounts[name] = SimpleNamespace(
                pw_name=name,
                pw_uid=4000 + len(accounts),
                pw_gid=group.gr_gid,
                pw_dir=values[5],
                pw_shell="/usr/sbin/nologin",
            )
            return True
        return False

    monkeypatch.setattr(bootstrap.grp, "getgrnam", get_group)
    monkeypatch.setattr(bootstrap.pwd, "getpwnam", get_account)
    backend = bootstrap.HostBootstrapBackend(
        source_root=TEMPLATES,
        root_prefix=tmp_path,
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
        command_runner=run,
        sudoers_validator=lambda _path: True,
    )

    backend._ensure_service_accounts()

    assert tuple(groups) == ("helixweave", "helixweave-api", "helixweave-candidate")
    assert tuple(accounts) == tuple(groups)
    assert len({account.pw_uid for account in accounts.values()}) == 3
    assert len({account.pw_gid for account in accounts.values()}) == 3
    assert all(
        "--no-create-home" in command for command in commands if "useradd" in command[0]
    )
    assert all("--groups" not in command for command in commands)


@pytest.mark.parametrize("unexpected_member", ("labadmin", "helixweave-candidate"))
def test_bootstrap_rejects_data_group_members_outside_fixed_service_boundary(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    unexpected_member: str,
) -> None:
    groups = {
        spec.group: SimpleNamespace(
            gr_name=spec.group,
            gr_gid=3100 + index,
            gr_mem=(unexpected_member,) if spec.group == "helixweave" else (),
        )
        for index, spec in enumerate(bootstrap.SYSTEM_ACCOUNTS)
    }
    accounts = {
        spec.user: SimpleNamespace(
            pw_name=spec.user,
            pw_uid=4100 + index,
            pw_gid=groups[spec.group].gr_gid,
            pw_dir=spec.home,
            pw_shell="/usr/sbin/nologin",
        )
        for index, spec in enumerate(bootstrap.SYSTEM_ACCOUNTS)
    }
    monkeypatch.setattr(bootstrap.grp, "getgrnam", groups.__getitem__)
    monkeypatch.setattr(bootstrap.pwd, "getpwnam", accounts.__getitem__)
    backend = bootstrap.HostBootstrapBackend(
        source_root=TEMPLATES,
        root_prefix=tmp_path,
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
        command_runner=lambda command: pytest.fail(f"unexpected command: {command}"),
        sudoers_validator=lambda _path: True,
    )

    with pytest.raises(bootstrap.BootstrapFailure) as caught:
        backend._ensure_service_accounts()

    assert caught.value.code == "BOOTSTRAP_ACCOUNT_FAILED"


def test_bootstrap_rejects_reused_fixed_runtime_uid(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    groups = {
        spec.group: SimpleNamespace(
            gr_name=spec.group,
            gr_gid=3200 + index,
            gr_mem=(),
        )
        for index, spec in enumerate(bootstrap.SYSTEM_ACCOUNTS)
    }
    accounts = {
        spec.user: SimpleNamespace(
            pw_name=spec.user,
            pw_uid=4200 if index < 2 else 4202,
            pw_gid=groups[spec.group].gr_gid,
            pw_dir=spec.home,
            pw_shell="/usr/sbin/nologin",
        )
        for index, spec in enumerate(bootstrap.SYSTEM_ACCOUNTS)
    }
    monkeypatch.setattr(bootstrap.grp, "getgrnam", groups.__getitem__)
    monkeypatch.setattr(bootstrap.pwd, "getpwnam", accounts.__getitem__)
    backend = bootstrap.HostBootstrapBackend(
        source_root=TEMPLATES,
        root_prefix=tmp_path,
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
        command_runner=lambda command: pytest.fail(f"unexpected command: {command}"),
        sudoers_validator=lambda _path: True,
    )

    with pytest.raises(bootstrap.BootstrapFailure) as caught:
        backend._ensure_service_accounts()

    assert caught.value.code == "BOOTSTRAP_ACCOUNT_FAILED"


def test_bootstrap_rejects_candidate_supplementary_group_membership(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    spec = next(
        item
        for item in bootstrap.SYSTEM_ACCOUNTS
        if item.user == "helixweave-candidate"
    )
    group = SimpleNamespace(gr_name=spec.group, gr_gid=3301, gr_mem=())
    account = SimpleNamespace(
        pw_name=spec.user,
        pw_uid=4301,
        pw_gid=group.gr_gid,
        pw_dir=spec.home,
        pw_shell="/usr/sbin/nologin",
    )
    monkeypatch.setattr(bootstrap.grp, "getgrnam", lambda _name: group)
    monkeypatch.setattr(bootstrap.pwd, "getpwnam", lambda _name: account)
    monkeypatch.setattr(
        bootstrap.os,
        "getgrouplist",
        lambda _user, primary: [primary, 9999],
    )
    backend = bootstrap.HostBootstrapBackend(
        source_root=TEMPLATES,
        root_prefix=tmp_path,
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
        command_runner=lambda command: pytest.fail(f"unexpected command: {command}"),
        sudoers_validator=lambda _path: True,
    )

    with pytest.raises(bootstrap.BootstrapFailure) as caught:
        backend._ensure_system_account(spec)

    assert caught.value.code == "BOOTSTRAP_ACCOUNT_FAILED"


@pytest.mark.parametrize(
    "helper_name",
    (
        "helixweave-active-service",
        "helixweave-db-prepare-launcher",
        "helixweave-operator-action",
        "helixweave-encode-runtime-prepare",
    ),
)
def test_bootstrap_closure_identity_includes_every_execution_critical_helper(
    tmp_path: Path,
    helper_name: str,
) -> None:
    original_source = tmp_path / "original-source"
    changed_source = tmp_path / "changed-source"
    shutil.copytree(TEMPLATES, original_source)
    shutil.copytree(TEMPLATES, changed_source)
    target = changed_source / helper_name
    target.write_bytes(target.read_bytes() + b"\n# changed stable boundary\n")

    identities = []
    for name, source in (("original", original_source), ("changed", changed_source)):
        root = tmp_path / name
        root.mkdir()
        backend = bootstrap.HostBootstrapBackend(
            source_root=source,
            root_prefix=root,
            owner_uid=os.getuid(),
            owner_gid=os.getgid(),
            command_runner=lambda command: pytest.fail(
                f"unexpected command: {command}"
            ),
            sudoers_validator=lambda _path: True,
        )
        identities.append(
            backend.apply(
                operation="install", invoking_user="labadmin"
            ).closure_identity
        )

    assert identities[0] != identities[1]


def test_bootstrap_update_switches_only_after_a_complete_new_closure(
    tmp_path: Path,
    monkeypatch,
) -> None:
    source = tmp_path / "source"
    shutil.copytree(TEMPLATES, source)
    root = tmp_path / "host"
    root.mkdir()
    backend = bootstrap.HostBootstrapBackend(
        source_root=source,
        root_prefix=root,
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
        command_runner=lambda command: pytest.fail(f"unexpected command: {command}"),
        sudoers_validator=lambda path: True,
    )
    backend.apply(operation="install", invoking_user="labadmin")
    current = root / "opt" / "helixweave" / "operator" / "current"
    first = os.readlink(current)
    stable_bytes = {
        spec.destination: root.joinpath(*spec.destination.parts[1:]).read_bytes()
        for spec in bootstrap.STABLE_BOUNDARY_SPECS
    }
    linked_targets = {
        spec.destination: os.readlink(root.joinpath(*spec.destination.parts[1:]))
        for spec in bootstrap.LINKED_BOUNDARY_SPECS
    }
    target = source / "helixweave.target"
    target.write_text(target.read_text() + "\n# revised\n")

    original_switch = backend._switch_closure

    def interrupted(identity: str) -> None:
        assert identity != first
        raise bootstrap.BootstrapFailure(
            "BOOTSTRAP_INSTALL_FAILED",
            "Operator boundary could not be installed.",
        )

    monkeypatch.setattr(backend, "_switch_closure", interrupted)
    with pytest.raises(bootstrap.BootstrapFailure):
        backend.apply(operation="update", invoking_user="labadmin")
    assert os.readlink(current) == first
    assert {
        spec.destination: root.joinpath(*spec.destination.parts[1:]).read_bytes()
        for spec in bootstrap.STABLE_BOUNDARY_SPECS
    } == stable_bytes
    assert {
        spec.destination: os.readlink(root.joinpath(*spec.destination.parts[1:]))
        for spec in bootstrap.LINKED_BOUNDARY_SPECS
    } == linked_targets

    monkeypatch.setattr(backend, "_switch_closure", original_switch)
    backend.apply(operation="update", invoking_user="labadmin")
    second = os.readlink(current)
    assert second != first
    assert (current.parent / second / "closure.json").is_file()


def test_bootstrap_update_refuses_to_mix_a_new_stable_boundary_with_old_bytes(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source"
    shutil.copytree(TEMPLATES, source)
    root = tmp_path / "host"
    root.mkdir()
    backend = bootstrap.HostBootstrapBackend(
        source_root=source,
        root_prefix=root,
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
        command_runner=lambda command: pytest.fail(f"unexpected command: {command}"),
        sudoers_validator=lambda path: True,
    )
    backend.apply(operation="install", invoking_user="labadmin")
    current = root / "opt" / "helixweave" / "operator" / "current"
    original = os.readlink(current)
    launcher = source / "helixweave-operator-launcher"
    launcher.write_text(launcher.read_text() + "\n# incompatible boundary\n")

    installed = root / "usr" / "libexec" / "helixweave-operator"
    original_bytes = installed.read_bytes()
    linked_targets = {
        spec.destination: os.readlink(root.joinpath(*spec.destination.parts[1:]))
        for spec in bootstrap.LINKED_BOUNDARY_SPECS
    }

    with pytest.raises(bootstrap.BootstrapFailure) as caught:
        backend.apply(operation="update", invoking_user="labadmin")

    assert caught.value.code == "BOOTSTRAP_BOUNDARY_UPDATE_REQUIRED"
    assert os.readlink(current) == original
    assert installed.read_bytes() == original_bytes
    assert {
        spec.destination: os.readlink(root.joinpath(*spec.destination.parts[1:]))
        for spec in bootstrap.LINKED_BOUNDARY_SPECS
    } == linked_targets
    with pytest.raises(bootstrap.BootstrapFailure) as reinstall:
        backend.apply(operation="install", invoking_user="labadmin")
    assert reinstall.value.code == "BOOTSTRAP_ALREADY_INSTALLED"


def _recovery_service_identity(
    *,
    deployment_identity: str,
    task_identity: str = TASK_IDENTITY,
    unit: str = "helixweave-worker.service",
) -> ServiceIdentity:
    return ServiceIdentity.create(
        unit=unit,
        deployment_identity=deployment_identity,
        task_identity=task_identity,
        main_pid=4321,
        process_start_ticks=8765,
        executable_device=31,
        executable_inode=41,
        cmdline_identity=IDENTITY,
        boot_identity=SERVICE_IDENTITY,
        invocation_identity=THIRD_IDENTITY,
        cgroup_identity=OLD_PLATFORM_IDENTITY,
        sockets=(),
    )


@dataclass
class _RecoveryProbe:
    service: ServiceIdentity
    running: bool
    calls: list[tuple[str, str, str]]

    def observe(
        self,
        *,
        unit: str,
        deployment_identity: str,
        task_identity: str,
    ) -> ServiceIdentity | None:
        self.calls.append((unit, deployment_identity, task_identity))
        if not self.running:
            return None
        assert unit == self.service.unit
        assert deployment_identity == self.service.deployment_identity
        assert task_identity == self.service.task_identity
        return self.service


@dataclass
class _RecoverySystemctl:
    probe: _RecoveryProbe
    calls: list[tuple[str, str]]

    def control(self, action: str, unit: str) -> None:
        self.calls.append((action, unit))
        if action == "start":
            self.probe.running = True
        elif action == "stop":
            self.probe.running = False
        else:
            assert action == "reset-failed"


def _systemd_recovery_controller(
    layout: DeploymentLayout,
    *,
    service: ServiceIdentity,
    running: bool,
) -> tuple[SystemdServiceController, _RecoveryProbe, _RecoverySystemctl]:
    layout.service_identities.mkdir(parents=True, mode=0o700)
    layout.service_identities.chmod(0o700)
    probe = _RecoveryProbe(service, running, [])
    systemctl = _RecoverySystemctl(probe, [])
    controller = SystemdServiceController(
        layout,
        systemctl=systemctl,
        probe=probe,
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
    )
    return controller, probe, systemctl


@pytest.mark.parametrize("already_running", (False, True))
def test_direct_start_recovery_uses_systemd_identity_boundary_idempotently(
    tmp_path: Path,
    already_running: bool,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / str(already_running))
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    states, state = _state_with_active_components(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
    )
    service = _recovery_service_identity(deployment_identity=OLD_PLATFORM_IDENTITY)
    services, probe, systemctl = _systemd_recovery_controller(
        layout,
        service=service,
        running=already_running,
    )
    active = {
        component: state.components[component].active
        for component in ("platform", "encode-runtime", "bulk-rnaseq-runtime")
    }
    record = OperatorTransaction.create(
        request_identity=IDENTITY,
        operation="start",
        task_identity=TASK_IDENTITY,
        deployment_identity=OLD_PLATFORM_IDENTITY,
        component=None,
        unit=service.unit,
        phase="recovery-required",
        failure_phase="service-starting",
        point_of_no_return=True,
        restart_units=(service.unit,),
        prior_active=active,
        candidate_active=active,
        prior_state_identity=state.identity,
        candidate_state_identity=state.identity,
    )
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=services,
        root_uid=owner_uid,
        root_gid=owner_gid,
    )

    recovered = controller.recover(record)

    assert recovered.phase == "complete"
    assert probe.running
    assert systemctl.calls == ([] if already_running else [("start", service.unit)])
    identity_path = layout.service_identities / f"{service.unit}.json"
    assert ServiceIdentity.from_dict(json.loads(identity_path.read_text())) == service
    assert stat.S_IMODE(identity_path.stat().st_mode) == 0o600


def test_direct_cleanup_recovery_finishes_a_live_systemd_stop(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "cleanup")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    states, state = _state_with_active_components(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
    )
    service = _recovery_service_identity(deployment_identity=OLD_PLATFORM_IDENTITY)
    services, probe, systemctl = _systemd_recovery_controller(
        layout,
        service=service,
        running=True,
    )
    assert (
        services.recover_observe(
            OperatorRequest(
                operation="status",
                task_identity=TASK_IDENTITY,
                deployment_identity=OLD_PLATFORM_IDENTITY,
                unit=service.unit,
            )
        )
        == service
    )
    probe.calls.clear()
    systemctl.calls.clear()
    active = {
        component: state.components[component].active
        for component in ("platform", "encode-runtime", "bulk-rnaseq-runtime")
    }
    record = OperatorTransaction.create(
        request_identity=IDENTITY,
        operation="cleanup",
        task_identity=TASK_IDENTITY,
        deployment_identity=OLD_PLATFORM_IDENTITY,
        component=None,
        unit=service.unit,
        phase="recovery-required",
        failure_phase="service-stopping",
        restart_units=(service.unit,),
        prior_running_units=(service.unit,),
        prior_active=active,
        candidate_active=active,
        prior_state_identity=state.identity,
        candidate_state_identity=state.identity,
        evidence={"service_identity": service.identity},
    )
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=services,
        root_uid=owner_uid,
        root_gid=owner_gid,
    )

    recovered = controller.recover(record)

    assert recovered.phase == "complete"
    assert not probe.running
    assert systemctl.calls == [
        ("stop", service.unit),
        ("reset-failed", service.unit),
    ]
    assert len(probe.calls) == 3


def test_pre_ponr_recovery_adopts_and_stops_unpersisted_candidate_service(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "candidate")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    states, prior = _state_with_active_components(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
    )
    candidate = prior.stage("platform", IDENTITY).activate("platform")
    service = _recovery_service_identity(deployment_identity=IDENTITY)
    services, probe, systemctl = _systemd_recovery_controller(
        layout,
        service=service,
        running=True,
    )
    prior_active = {
        component: prior.components[component].active
        for component in ("platform", "encode-runtime", "bulk-rnaseq-runtime")
    }
    candidate_active = {
        component: candidate.components[component].active
        for component in ("platform", "encode-runtime", "bulk-rnaseq-runtime")
    }
    record = OperatorTransaction.create(
        request_identity=SERVICE_IDENTITY,
        operation="activate",
        task_identity=TASK_IDENTITY,
        deployment_identity=IDENTITY,
        component="platform",
        unit=None,
        phase="recovery-required",
        failure_phase="writers-stopped",
        write_fence=True,
        restart_units=(service.unit,),
        prior_active=prior_active,
        candidate_active=candidate_active,
        prior_state_identity=prior.identity,
        candidate_state_identity=candidate.identity,
    )
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=services,
        root_uid=owner_uid,
        root_gid=owner_gid,
    )

    recovered = controller.recover(record)

    assert recovered.phase == "aborted"
    assert states.read() == prior
    assert not probe.running
    assert systemctl.calls == [("stop", service.unit)]
    identity_path = layout.service_identities / f"{service.unit}.json"
    assert ServiceIdentity.from_dict(json.loads(identity_path.read_text())) == service


def test_stage_recovery_falls_back_when_candidate_generation_was_not_written(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "stage")
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    states = _supported_state_store(layout, owner_uid, owner_gid)
    prior = states.initialize(
        expected_owner_uid=owner_uid,
        expected_owner_gid=owner_gid,
    )
    candidate = prior.stage("platform", IDENTITY)

    def interrupt(point: str) -> None:
        if point == "transaction-prepared":
            raise RuntimeError("candidate generation was not written")

    with pytest.raises(RuntimeError, match="candidate generation was not written"):
        states.commit(
            candidate,
            operation="stage-platform",
            expected_current_identity=prior.identity,
            expected_owner_uid=owner_uid,
            expected_owner_gid=owner_gid,
            fault=interrupt,
        )
    prior_active = {
        component: prior.components[component].active
        for component in ("platform", "encode-runtime", "bulk-rnaseq-runtime")
    }
    record = OperatorTransaction.create(
        request_identity=THIRD_IDENTITY,
        operation="stage",
        task_identity=TASK_IDENTITY,
        deployment_identity=IDENTITY,
        component="platform",
        unit=None,
        phase="recovery-required",
        failure_phase="candidate-selected",
        prior_active=prior_active,
        candidate_active=prior_active,
        prior_state_identity=prior.identity,
        candidate_state_identity=candidate.identity,
    )
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=FakeServiceController(),
        root_uid=owner_uid,
        root_gid=owner_gid,
    )

    recovered = controller.recover(record)

    assert recovered.phase == "aborted"
    assert states.read() == prior
    assert states.pending_transactions() == ()
    assert len(tuple(layout.state_transactions.glob("*.recovered.json"))) == 1
    assert not (layout.state_generations / candidate.identity).exists()


@pytest.mark.parametrize("daemon_reload_status", (0, 1))
def test_uninstall_recovery_restores_real_pending_boundary_before_reload(
    tmp_path: Path,
    daemon_reload_status: int,
) -> None:
    root = tmp_path / str(daemon_reload_status)
    layout = DeploymentLayout.isolated(root)
    owner_uid = os.getuid()
    owner_gid = os.getgid()
    states, state = _state_with_active_components(
        layout,
        owner_uid=owner_uid,
        owner_gid=owner_gid,
    )
    for logical in UNINSTALL_BOUNDARY_FILES:
        path = root / str(logical).lstrip("/")
        path.parent.mkdir(parents=True, exist_ok=True)
        linked_target = UNINSTALL_LINKED_BOUNDARY_TARGETS.get(logical)
        if linked_target is None:
            path.write_bytes(f"boundary:{logical}\n".encode())
            path.chmod(0o444)
        else:
            path.symlink_to(linked_target)
        if logical not in {
            Path("/usr/libexec/helixweave-operator"),
            Path("/etc/sudoers.d/helixweave-operator"),
        }:
            os.replace(
                path,
                path.with_name(f".{path.name}.helixweave-uninstall-pending"),
            )
    executor = RecordingCommandExecutor(CommandResult(daemon_reload_status), [])
    uninstaller = HostBoundaryUninstaller(
        owner_uid=owner_uid,
        owner_gid=owner_gid,
        executor=executor,
        root_prefix=root,
    )
    active = {
        component: state.components[component].active
        for component in ("platform", "encode-runtime", "bulk-rnaseq-runtime")
    }
    record = OperatorTransaction.create(
        request_identity=IDENTITY,
        operation="uninstall",
        task_identity=TASK_IDENTITY,
        deployment_identity=state.identity,
        component=None,
        unit=None,
        phase="recovery-required",
        failure_phase="writers-stopped",
        write_fence=True,
        prior_active=active,
        candidate_active=active,
        prior_state_identity=state.identity,
        candidate_state_identity=state.identity,
    )
    controller = HostDeploymentActionController(
        layout,
        states=states,
        services=_TrackingServices({}, []),
        boundary_uninstaller=uninstaller,
        root_uid=owner_uid,
        root_gid=owner_gid,
    )

    if daemon_reload_status == 0:
        assert controller.recover(record).phase == "aborted"
    else:
        with pytest.raises(DeploymentError) as captured:
            controller.recover(record)
        assert captured.value.issue.code == "OPERATOR_UNINSTALL_FAILED"

    assert executor.calls == [(("/usr/bin/systemctl", "daemon-reload"), 15.0)]
    for logical in UNINSTALL_BOUNDARY_FILES:
        path = root / str(logical).lstrip("/")
        assert path.exists() or path.is_symlink()
        assert not path.with_name(f".{path.name}.helixweave-uninstall-pending").exists()
