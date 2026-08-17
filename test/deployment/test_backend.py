from __future__ import annotations

from dataclasses import dataclass
import json
import os
from pathlib import Path
import stat
import subprocess
from types import SimpleNamespace

import pytest

import encode_pipeline.deployment.backend as backend_module
from encode_pipeline.deployment.admission import (
    DeferredDatabaseSchemaObserver,
    DeferredNativeContractResolver,
)
from encode_pipeline.deployment.backend import (
    IngressPublisher,
    OperatorExecution,
    ProductionCommandBackend,
    ServiceObservation,
    SudoOperatorClient,
)
from encode_pipeline.deployment.cli import OPERATION_RESULT_SCHEMA, STATUS_RESULT_SCHEMA
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import PLATFORM
from encode_pipeline.deployment.operator import (
    OPERATOR_RECEIPT_SCHEMA,
    SERVICE_UNITS,
    OperatorObservation,
    ServiceIdentity,
)
from encode_pipeline.deployment.operator_action import (
    VERIFICATION_CHECKS,
    DeploymentActionReceipt,
    ReadinessCheck,
)
from .support import manager_for, manifest_for, write_bundle


IDENTITY = f"sha256-{'a' * 64}"
OTHER_IDENTITY = f"sha256-{'b' * 64}"
TASK = f"task-{'1' * 32}"
SERVICE_TASK = f"task-{'2' * 32}"


def _receipt(**values: object) -> bytes:
    return (
        json.dumps(values, sort_keys=True, separators=(",", ":")).encode("utf-8")
        + b"\n"
    )


def _service(*, deployment_identity: str = IDENTITY) -> ServiceIdentity:
    return ServiceIdentity.create(
        unit="helixweave-worker.service",
        deployment_identity=deployment_identity,
        task_identity=SERVICE_TASK,
        main_pid=1234,
        process_start_ticks=42,
        executable_device=8,
        executable_inode=9,
        cmdline_identity=OTHER_IDENTITY,
        boot_identity=IDENTITY,
        invocation_identity=OTHER_IDENTITY,
        cgroup_identity=IDENTITY,
        sockets=(),
    )


def _verification(
    *,
    status: str = "observed",
    native_platform_identity: str | None = IDENTITY,
    database_identity: str = OTHER_IDENTITY,
    frontend_identity: str = IDENTITY,
    readiness_overrides: dict[str, ReadinessCheck] | None = None,
) -> DeploymentActionReceipt:
    readiness = {
        name: ReadinessCheck("ready", "READY", IDENTITY) for name in VERIFICATION_CHECKS
    }
    for name in (
        "encode-runtime-native",
        "bulk-runtime-native",
        "worker",
        "docker",
    ):
        readiness[name] = ReadinessCheck(
            "not-applicable",
            "NOT_APPLICABLE",
        )
    readiness["database-schema"] = ReadinessCheck(
        "ready",
        "READY",
        database_identity,
    )
    readiness["frontend"] = ReadinessCheck(
        "ready",
        "READY",
        frontend_identity,
    )
    readiness.update(readiness_overrides or {})
    return DeploymentActionReceipt.create(
        request_identity=IDENTITY,
        status=status,
        compatibility="compatible",
        database_after_identity=database_identity,
        accepted_schema_heads=("schema-v1",),
        target_schema_heads=("schema-v1",),
        migration_inventory_identity=IDENTITY,
        known_schema_revisions=(),
        rollback_supported=True,
        api_contract_sha256="a" * 64,
        native_identities={
            "platform": native_platform_identity,
            "encode-runtime": None,
            "bulk-rnaseq-runtime": None,
        },
        frontend_identity=frontend_identity,
        reference_compatibility_identity=IDENTITY,
        readiness=readiness,
    )


def test_operator_status_is_unbound_but_receipt_and_service_are_fully_bound() -> None:
    service = _service()
    calls: list[tuple[str, ...]] = []

    def run(arguments: tuple[str, ...]) -> OperatorExecution:
        calls.append(arguments)
        return OperatorExecution(
            0,
            _receipt(
                schema_version=OPERATOR_RECEIPT_SCHEMA,
                operation="status",
                state="running",
                task_identity=TASK,
                deployment_identity=IDENTITY,
                unit="helixweave-worker.service",
                service_identity=service.to_dict(),
            ),
            b"",
        )

    observed = SudoOperatorClient(runner=run).status_service(
        "helixweave-worker.service",
        IDENTITY,
        TASK,
    )

    assert calls == [("status", "helixweave-worker.service", IDENTITY, TASK)]
    assert observed == ServiceObservation("running", service.identity)
    assert service.task_identity != TASK


def test_operator_client_rejects_extra_or_mismatched_receipt_fields() -> None:
    content = _receipt(
        schema_version=OPERATOR_RECEIPT_SCHEMA,
        operation="stage",
        state="staged",
        task_identity=TASK,
        deployment_identity=IDENTITY,
        component=PLATFORM,
        unexpected="value",
    )
    client = SudoOperatorClient(
        runner=lambda _arguments: OperatorExecution(0, content, b"")
    )

    with pytest.raises(DeploymentError) as caught:
        client.stage(PLATFORM, IDENTITY, TASK)

    assert caught.value.issue.code == "DEPLOYMENT_OPERATOR_UNAVAILABLE"


@pytest.mark.parametrize(
    ("failure", "recoverable"),
    (
        ("exception", True),
        ("wrong-type", False),
        ("exit-unavailable", True),
        ("stderr", False),
        ("empty", False),
        ("invalid-json", False),
        ("duplicate-key", False),
        ("noncanonical", False),
        ("oversized", False),
    ),
)
def test_operator_client_fails_closed_on_untrusted_execution_results(
    failure: str,
    recoverable: bool,
) -> None:
    def run(_arguments: tuple[str, ...]) -> object:
        if failure == "exception":
            raise OSError("/private/operator token=secret")
        if failure == "wrong-type":
            return object()
        if failure == "exit-unavailable":
            return OperatorExecution(69, b"", b"")
        if failure == "stderr":
            return OperatorExecution(0, b"{}\n", b"private stderr")
        if failure == "empty":
            return OperatorExecution(0, b"", b"")
        if failure == "invalid-json":
            return OperatorExecution(0, b"not-json\n", b"")
        if failure == "duplicate-key":
            return OperatorExecution(
                0,
                b'{"schema_version":"one","schema_version":"two"}\n',
                b"",
            )
        if failure == "noncanonical":
            return OperatorExecution(0, b"{ }\n", b"")
        return OperatorExecution(
            0,
            b"x" * (backend_module._MAX_OPERATOR_RECEIPT_BYTES + 1),
            b"",
        )

    client = SudoOperatorClient(runner=run)  # type: ignore[arg-type]
    with pytest.raises(DeploymentError) as caught:
        client.stage(PLATFORM, IDENTITY, TASK)

    assert caught.value.issue.code == "DEPLOYMENT_OPERATOR_UNAVAILABLE"
    assert caught.value.issue.recoverable is recoverable
    assert "private" not in str(caught.value)
    assert "secret" not in str(caught.value)


def test_operator_client_accepts_stopped_services_and_exact_mutation_receipts() -> None:
    calls: list[tuple[str, ...]] = []

    def run(arguments: tuple[str, ...]) -> OperatorExecution:
        calls.append(arguments)
        operation = arguments[0]
        if operation == "status":
            unit, deployment_identity, task = arguments[1:]
            return OperatorExecution(
                0,
                _receipt(
                    schema_version=OPERATOR_RECEIPT_SCHEMA,
                    operation="status",
                    state="stopped",
                    task_identity=task,
                    deployment_identity=deployment_identity,
                    unit=unit,
                ),
                b"",
            )
        component, deployment_identity, task = arguments[1:]
        return OperatorExecution(
            0,
            _receipt(
                schema_version=OPERATOR_RECEIPT_SCHEMA,
                operation=operation,
                state={"activate": "activated", "rollback": "rolled-back"}[operation],
                task_identity=task,
                deployment_identity=deployment_identity,
                component=component,
            ),
            b"",
        )

    client = SudoOperatorClient(runner=run)
    client.activate(PLATFORM, IDENTITY, TASK)
    client.rollback(PLATFORM, IDENTITY, TASK)
    assert client.status_service("helixweave-api.service", IDENTITY, TASK) == (
        ServiceObservation("stopped")
    )
    assert calls == [
        ("activate", PLATFORM, IDENTITY, TASK),
        ("rollback", PLATFORM, IDENTITY, TASK),
        ("status", "helixweave-api.service", IDENTITY, TASK),
    ]


def test_operator_service_observation_rejects_invalid_unit_and_identity_evidence() -> (
    None
):
    client = SudoOperatorClient(
        runner=lambda _arguments: OperatorExecution(0, b"{}\n", b"")
    )
    with pytest.raises(DeploymentError) as invalid_unit:
        client.status_service("caller-selected.service", IDENTITY, TASK)
    assert invalid_unit.value.issue.code == "DEPLOYMENT_SERVICE_OBSERVATION_INVALID"

    def running_receipt(service_identity: object) -> OperatorExecution:
        return OperatorExecution(
            0,
            _receipt(
                schema_version=OPERATOR_RECEIPT_SCHEMA,
                operation="status",
                state="running",
                task_identity=TASK,
                deployment_identity=IDENTITY,
                unit="helixweave-worker.service",
                service_identity=service_identity,
            ),
            b"",
        )

    for evidence in ({}, _service(deployment_identity=OTHER_IDENTITY).to_dict()):
        malformed = SudoOperatorClient(
            runner=lambda _arguments, evidence=evidence: running_receipt(evidence)
        )
        with pytest.raises(DeploymentError) as caught:
            malformed.status_service("helixweave-worker.service", IDENTITY, TASK)
        assert caught.value.issue.code == "DEPLOYMENT_OPERATOR_UNAVAILABLE"


def test_operator_observation_uses_the_fixed_state_bound_grammar() -> None:
    observation = OperatorObservation.create(
        state_identity=IDENTITY,
        active={
            "platform": OTHER_IDENTITY,
            "encode-runtime": None,
            "bulk-rnaseq-runtime": None,
        },
        database_schema_identity=OTHER_IDENTITY,
        database_schema_heads=("schema-v1",),
        services={unit: None for unit in SERVICE_UNITS},
    )
    calls: list[tuple[str, ...]] = []

    def run(arguments: tuple[str, ...]) -> OperatorExecution:
        calls.append(arguments)
        return OperatorExecution(
            0,
            _receipt(
                schema_version=OPERATOR_RECEIPT_SCHEMA,
                operation="observe",
                state="observed",
                task_identity=TASK,
                deployment_identity=IDENTITY,
                observation=observation.to_dict(),
            ),
            b"",
        )

    result = SudoOperatorClient(runner=run).observe(IDENTITY, TASK)

    assert result == observation
    assert calls == [("observe", IDENTITY, TASK)]


def test_operator_observation_preserves_headless_database_identity() -> None:
    observation = OperatorObservation.create(
        state_identity=IDENTITY,
        active={
            component: None
            for component in (
                "platform",
                "encode-runtime",
                "bulk-rnaseq-runtime",
            )
        },
        database_schema_identity=OTHER_IDENTITY,
        database_schema_heads=(),
        services={unit: None for unit in SERVICE_UNITS},
    )

    assert observation.database_schema_identity == OTHER_IDENTITY
    assert observation.database_schema_heads == ()


def test_operator_verification_uses_the_fixed_state_bound_grammar() -> None:
    verification = _verification()
    calls: list[tuple[str, ...]] = []

    def run(arguments: tuple[str, ...]) -> OperatorExecution:
        calls.append(arguments)
        return OperatorExecution(
            0,
            _receipt(
                schema_version=OPERATOR_RECEIPT_SCHEMA,
                operation="verify",
                state="verified",
                task_identity=TASK,
                deployment_identity=IDENTITY,
                verification=verification.to_dict(),
            ),
            b"",
        )

    result = SudoOperatorClient(runner=run).verify(IDENTITY, TASK)

    assert result == verification
    assert calls == [("verify", IDENTITY, TASK)]


def test_operator_runner_uses_only_fixed_sudo_n_argv_and_environment(
    monkeypatch,
) -> None:
    observed: dict[str, object] = {}

    def run(argv, **kwargs):
        observed["argv"] = argv
        observed.update(kwargs)
        return subprocess.CompletedProcess(argv, 0, stdout=b"{}\n", stderr=b"")

    monkeypatch.setattr(backend_module.subprocess, "run", run)

    result = backend_module._run_operator(("stage", PLATFORM, IDENTITY, TASK))

    assert result == OperatorExecution(0, b"{}\n", b"")
    assert observed["argv"] == (
        "/usr/bin/sudo",
        "-n",
        "/usr/libexec/helixweave-operator",
        "stage",
        PLATFORM,
        IDENTITY,
        TASK,
    )
    assert observed["stdin"] == subprocess.DEVNULL
    assert observed["shell"] is False
    assert observed["timeout"] == 2_700
    assert observed["cwd"] == "/"
    assert observed["env"] == {
        "LANG": "C.UTF-8",
        "LC_ALL": "C.UTF-8",
        "PATH": "/usr/sbin:/usr/bin:/sbin:/bin",
    }


@pytest.mark.parametrize(
    "arguments",
    (
        ("activate", "bulk-rnaseq-runtime", IDENTITY, TASK),
        ("rollback", "bulk-rnaseq-runtime", IDENTITY, TASK),
        ("verify", IDENTITY, TASK),
    ),
)
def test_bulk_materialization_and_verification_have_a_larger_outer_timeout(
    monkeypatch,
    arguments: tuple[str, ...],
) -> None:
    observed: dict[str, object] = {}

    def run(argv, **kwargs):
        observed["argv"] = argv
        observed.update(kwargs)
        return subprocess.CompletedProcess(argv, 0, stdout=b"{}\n", stderr=b"")

    monkeypatch.setattr(backend_module.subprocess, "run", run)

    backend_module._run_operator(arguments)

    assert observed["timeout"] == 15_000


def test_supported_composition_separates_state_reader_group_from_release_owner(
    tmp_path: Path,
    monkeypatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    monkeypatch.setattr(
        backend_module.DeploymentLayout,
        "supported",
        classmethod(lambda cls: layout),
    )
    monkeypatch.setattr(
        backend_module.grp,
        "getgrnam",
        lambda name: SimpleNamespace(
            gr_gid=1234 if name == "helixweave-operators" else 2345
        ),
    )

    backend = ProductionCommandBackend.supported()

    assert backend.manager.ownership.uid == 0
    assert backend.manager.ownership.gid == 0
    assert backend.manager.states.reader_gid == 1234
    assert backend.manager.states.service_gid == 2345
    assert backend.manager.states.verify_environment_content is False
    assert isinstance(
        backend.manager.contract_resolver,
        DeferredNativeContractResolver,
    )
    assert isinstance(
        backend.manager.schema_observer,
        DeferredDatabaseSchemaObserver,
    )


def test_ingress_publication_is_flat_atomic_and_reusable(tmp_path: Path) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manifest, payload = manifest_for(PLATFORM)
    source = write_bundle(tmp_path / "source.tar", manifest, payload)
    component_directory = layout.ingress / PLATFORM
    component_directory.mkdir(parents=True)
    component_directory.chmod(0o2730)
    publisher = IngressPublisher(
        layout,
        directory_uid=os.getuid(),
        directory_gid=os.getgid(),
        file_uid=os.getuid(),
    )

    published = publisher.publish(source, manifest, TASK)
    same = publisher.publish(source, manifest, TASK)

    assert published == component_directory / f"{manifest.identity}.tar"
    assert same == published
    observed = published.lstat()
    assert stat.S_ISREG(observed.st_mode)
    assert observed.st_nlink == 1
    assert stat.S_IMODE(observed.st_mode) == 0o440
    assert not tuple(component_directory.glob("*.partial"))


@dataclass
class _Ingress:
    calls: list[tuple[Path, object, str]]

    def publish(self, path, manifest, task):
        self.calls.append((path, manifest, task))
        return path


class _Operator:
    def __init__(self, manager, bundle: Path) -> None:
        self.manager = manager
        self.bundle = bundle
        self.calls: list[tuple[str, ...]] = []

    def stage(self, component: str, identity: str, task: str) -> None:
        self.calls.append(("stage", component, identity, task))
        self.manager.stage(
            self.bundle,
            expected_owner_uid=os.getuid(),
            expected_owner_gid=os.getgid(),
        )

    def activate(self, component: str, identity: str, task: str) -> None:
        self.calls.append(("activate", component, identity, task))
        self.manager.activate(component, expected_staged_identity=identity)

    def rollback(self, component: str, identity: str, task: str) -> None:
        self.calls.append(("rollback", component, identity, task))
        self.manager.rollback(component, expected_previous_identity=identity)

    def status_service(
        self,
        unit: str,
        deployment_identity: str,
        task: str,
    ) -> ServiceObservation:
        self.calls.append(("status", unit, deployment_identity, task))
        if unit in {"helixweave-api.service", "helixweave-worker.service"}:
            return ServiceObservation("running", IDENTITY)
        return ServiceObservation("stopped")

    def observe(self, state_identity: str, task: str) -> OperatorObservation:
        status = self.manager.status()
        self.calls.append(("observe", state_identity, task))
        return OperatorObservation.create(
            state_identity=status.state.identity,
            active={
                component: status.state.components[component].active
                for component in (
                    "platform",
                    "encode-runtime",
                    "bulk-rnaseq-runtime",
                )
            },
            database_schema_identity=OTHER_IDENTITY,
            database_schema_heads=("schema-v1",),
            services={
                unit: (
                    IDENTITY
                    if unit in {"helixweave-api.service", "helixweave-worker.service"}
                    else None
                )
                for unit in SERVICE_UNITS
            },
        )

    def verify(
        self,
        state_identity: str,
        task: str,
    ) -> DeploymentActionReceipt:
        self.calls.append(("verify", state_identity, task))
        status = self.manager.status()
        frontend_identity = backend_module._frontend_identity(status)
        assert frontend_identity is not None
        return _verification(frontend_identity=frontend_identity)


def test_upgrade_composes_ingress_operator_and_root_state_reread(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    manifest, payload = manifest_for(PLATFORM)
    bundle = write_bundle(tmp_path / "platform.tar", manifest, payload)
    ingress = _Ingress([])
    operator = _Operator(manager, bundle)
    backend = ProductionCommandBackend(
        layout=layout,
        manager=manager,
        operator=operator,
        ingress=ingress,
        task_factory=lambda: TASK,
    )

    result = backend.upgrade(PLATFORM, bundle)

    stage_task = backend_module._operation_task(TASK, "stage")
    activate_task = backend_module._operation_task(TASK, "activate")

    assert result.value == {
        "schema_version": OPERATION_RESULT_SCHEMA,
        "operation": "upgrade",
        "state": "activated",
        "component": PLATFORM,
        "deployment_identity": manifest.identity,
        "state_identity": manager.status().state.identity,
    }
    assert stage_task != activate_task
    assert ingress.calls == [(bundle, manifest, stage_task)]
    assert operator.calls == [
        ("stage", PLATFORM, manifest.identity, stage_task),
        ("activate", PLATFORM, manifest.identity, activate_task),
    ]


def test_install_stages_exactly_the_inspected_and_published_bundle(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    manifest, payload = manifest_for(PLATFORM)
    bundle = write_bundle(tmp_path / "platform.tar", manifest, payload)
    ingress = _Ingress([])
    operator = _Operator(manager, bundle)
    backend = ProductionCommandBackend(
        layout=layout,
        manager=manager,
        operator=operator,
        ingress=ingress,
        task_factory=lambda: TASK,
    )

    result = backend.install(PLATFORM, bundle)

    assert result.value == {
        "schema_version": OPERATION_RESULT_SCHEMA,
        "operation": "install",
        "state": "staged",
        "component": PLATFORM,
        "deployment_identity": manifest.identity,
        "state_identity": manager.status().state.identity,
    }
    assert ingress.calls == [(bundle, manifest, TASK)]
    assert operator.calls == [("stage", PLATFORM, manifest.identity, TASK)]
    assert manager.status().state.components[PLATFORM].active is None


def test_rollback_requires_and_activates_the_exact_previous_identity(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    previous, previous_payload = manifest_for(PLATFORM, version="1.0.0")
    active, active_payload = manifest_for(PLATFORM, version="2.0.0")
    previous_bundle = write_bundle(
        tmp_path / "platform-previous.tar",
        previous,
        previous_payload,
    )
    active_bundle = write_bundle(
        tmp_path / "platform-active.tar",
        active,
        active_payload,
    )
    manager.stage(previous_bundle)
    manager.activate(PLATFORM, expected_staged_identity=previous.identity)
    manager.stage(active_bundle)
    manager.activate(PLATFORM, expected_staged_identity=active.identity)
    operator = _Operator(manager, active_bundle)
    backend = ProductionCommandBackend(
        layout=layout,
        manager=manager,
        operator=operator,
        ingress=_Ingress([]),
        task_factory=lambda: TASK,
    )

    result = backend.rollback(PLATFORM, previous.identity)

    assert result.value == {
        "schema_version": OPERATION_RESULT_SCHEMA,
        "operation": "rollback",
        "state": "rolled-back",
        "component": PLATFORM,
        "deployment_identity": previous.identity,
        "state_identity": manager.status().state.identity,
    }
    assert operator.calls == [("rollback", PLATFORM, previous.identity, TASK)]
    assert manager.status().state.components[PLATFORM].active == previous.identity
    with pytest.raises(DeploymentError) as mismatch:
        backend.rollback(PLATFORM, OTHER_IDENTITY)
    assert mismatch.value.issue.code == "DEPLOYMENT_PREVIOUS_IDENTITY_MISMATCH"
    assert operator.calls == [("rollback", PLATFORM, previous.identity, TASK)]


def test_status_uses_root_observation_when_the_local_database_is_unreadable(
    tmp_path: Path,
    monkeypatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    manifest, payload = manifest_for(PLATFORM)
    bundle = write_bundle(tmp_path / "platform.tar", manifest, payload)
    manager.stage(bundle)
    manager.activate(PLATFORM, expected_staged_identity=manifest.identity)

    def unavailable(_state):
        raise fail(
            "DEPLOYMENT_SCHEMA_OBSERVATION_FAILED",
            "Database schema observation failed.",
        )

    monkeypatch.setattr(manager, "observe_database_schema", unavailable)
    operator = _Operator(manager, bundle)
    backend = ProductionCommandBackend(
        layout=layout,
        manager=manager,
        operator=operator,
        ingress=_Ingress([]),
        task_factory=lambda: TASK,
    )

    result = backend.status()

    assert result.value["schema_version"] == STATUS_RESULT_SCHEMA
    assert result.value["database_schema_identity"] == OTHER_IDENTITY
    assert result.value["database_schema_reason_code"] == "DATABASE_READY"


def test_status_fails_closed_instead_of_zeroing_unknown_operator_counts(
    tmp_path: Path,
    monkeypatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    manifest, payload = manifest_for(PLATFORM)
    bundle = write_bundle(tmp_path / "platform.tar", manifest, payload)
    manager.stage(bundle)
    manager.activate(PLATFORM, expected_staged_identity=manifest.identity)
    operator = _Operator(manager, bundle)

    def unavailable(*_args):
        raise RuntimeError("unavailable")

    monkeypatch.setattr(operator, "observe", unavailable)
    backend = ProductionCommandBackend(
        layout=layout,
        manager=manager,
        operator=operator,
        ingress=_Ingress([]),
        task_factory=lambda: TASK,
    )

    with pytest.raises(DeploymentError) as captured:
        backend.status()

    assert captured.value.issue.code == "OPERATOR_OBSERVATION_UNAVAILABLE"


def test_status_marks_operator_pending_or_recovery_as_interrupted(
    tmp_path: Path,
    monkeypatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    manifest, payload = manifest_for(PLATFORM)
    bundle = write_bundle(tmp_path / "platform.tar", manifest, payload)
    manager.stage(bundle)
    manager.activate(PLATFORM, expected_staged_identity=manifest.identity)
    operator = _Operator(manager, bundle)
    original_observe = operator.observe

    def pending(state_identity: str, task: str) -> OperatorObservation:
        observed = original_observe(state_identity, task)
        return OperatorObservation.create(
            state_identity=observed.state_identity,
            active=observed.active,
            database_schema_identity=observed.database_schema_identity,
            database_schema_heads=observed.database_schema_heads,
            services=observed.services,
            operator_boundary=observed.operator_boundary,
            operator_pending_count=1,
        )

    monkeypatch.setattr(operator, "observe", pending)
    backend = ProductionCommandBackend(
        layout=layout,
        manager=manager,
        operator=operator,
        ingress=_Ingress([]),
        task_factory=lambda: TASK,
    )

    result = backend.status()

    assert result.value["operator_pending_transaction_count"] == 1
    assert result.value["interrupted"] is True


def test_storage_verification_never_crosses_native_or_database_boundaries(
    tmp_path: Path,
    monkeypatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    manifest, payload = manifest_for(PLATFORM)
    bundle = write_bundle(tmp_path / "platform.tar", manifest, payload)
    manager.stage(bundle)
    manager.activate(PLATFORM, expected_staged_identity=manifest.identity)

    def forbidden(*_args, **_kwargs):
        raise RuntimeError("private runtime boundary was crossed")

    monkeypatch.setattr(manager.contract_resolver, "resolve", forbidden)
    monkeypatch.setattr(manager.schema_observer, "observe", forbidden)

    verified = manager.verify_storage()

    assert verified.state.components[PLATFORM].active == manifest.identity


def test_verify_and_doctor_consume_the_same_root_schema_observation(
    tmp_path: Path,
    monkeypatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    manifest, payload = manifest_for(PLATFORM)
    bundle = write_bundle(tmp_path / "platform.tar", manifest, payload)
    manager.stage(bundle)
    manager.activate(PLATFORM, expected_staged_identity=manifest.identity)

    def unavailable(_state):
        raise fail(
            "DEPLOYMENT_SCHEMA_OBSERVATION_FAILED",
            "Database schema observation failed.",
        )

    monkeypatch.setattr(manager, "observe_database_schema", unavailable)
    backend = ProductionCommandBackend(
        layout=layout,
        manager=manager,
        operator=_Operator(manager, bundle),
        ingress=_Ingress([]),
        task_factory=lambda: TASK,
    )

    verification = backend.verify()
    report = backend.doctor()

    assert verification.value["verified"] is True
    assert verification.value["database_schema_identity"] == OTHER_IDENTITY
    database = next(
        item for item in report.value["checks"] if item["check_id"] == "database"
    )
    assert database == {
        "check_id": "database",
        "category": "database",
        "state": "pass",
        "reason_code": "DATABASE_READY",
        "evidence_identity": OTHER_IDENTITY,
    }


def test_verify_rejects_native_evidence_nullability_that_does_not_match_state(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    manifest, payload = manifest_for(PLATFORM)
    bundle = write_bundle(tmp_path / "platform.tar", manifest, payload)
    manager.stage(bundle)
    manager.activate(PLATFORM, expected_staged_identity=manifest.identity)
    operator = _Operator(manager, bundle)
    operator.verify = lambda _state, _task: _verification(  # type: ignore[method-assign]
        native_platform_identity=None
    )
    backend = ProductionCommandBackend(
        layout=layout,
        manager=manager,
        operator=operator,
        ingress=_Ingress([]),
        task_factory=lambda: TASK,
    )

    result = backend.verify()

    assert result.value["verified"] is False


def test_verify_is_static_while_doctor_owns_live_service_and_reference_health(
    tmp_path: Path,
    monkeypatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    manifest, payload = manifest_for(PLATFORM)
    bundle = write_bundle(tmp_path / "platform.tar", manifest, payload)
    manager.stage(bundle)
    manager.activate(PLATFORM, expected_staged_identity=manifest.identity)
    operator = _Operator(manager, bundle)
    frontend_identity = backend_module._frontend_identity(manager.status())
    assert frontend_identity is not None

    monkeypatch.setattr(
        operator,
        "verify",
        lambda _state, _task: _verification(
            frontend_identity=frontend_identity,
            readiness_overrides={
                "references": ReadinessCheck("not-ready", "REFERENCE_NOT_READY"),
                "redis": ReadinessCheck("unavailable", "REDIS_UNAVAILABLE"),
                "worker": ReadinessCheck("unavailable", "SERVICE_STOPPED"),
            },
        ),
    )
    monkeypatch.setattr(
        operator,
        "status_service",
        lambda _unit, _deployment, _task: ServiceObservation("stopped"),
    )
    backend = ProductionCommandBackend(
        layout=layout,
        manager=manager,
        operator=operator,
        ingress=_Ingress([]),
        task_factory=lambda: TASK,
    )

    verification = backend.verify()
    report = backend.doctor()

    assert verification.exit_code == 0
    assert verification.value["verified"] is True
    assert all(
        service["state"] != "running"
        for service in verification.value["deployment"]["services"].values()
    )
    assert report.value["ready"] is False
    checks = {item["check_id"]: item for item in report.value["checks"]}
    assert checks["api"]["state"] == "fail"
    assert checks["redis"]["reason_code"] == "REDIS_UNAVAILABLE"
    assert checks["references"]["state"] == "warning"
    assert checks["references"]["reason_code"] == "REFERENCES_INCOMPLETE"


def test_verify_distinguishes_operator_unavailability_from_invalid_content(
    tmp_path: Path,
    monkeypatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    manifest, payload = manifest_for(PLATFORM)
    bundle = write_bundle(tmp_path / "platform.tar", manifest, payload)
    manager.stage(bundle)
    manager.activate(PLATFORM, expected_staged_identity=manifest.identity)
    operator = _Operator(manager, bundle)

    def unavailable(*_args):
        raise RuntimeError("private operator failure")

    monkeypatch.setattr(operator, "verify", unavailable)
    backend = ProductionCommandBackend(
        layout=layout,
        manager=manager,
        operator=operator,
        ingress=_Ingress([]),
        task_factory=lambda: TASK,
    )

    with pytest.raises(DeploymentError) as captured:
        backend.verify()

    assert captured.value.issue.code == "DEPLOYMENT_OPERATOR_UNAVAILABLE"
    assert captured.value.issue.recoverable is True


def test_doctor_reports_operator_unavailable_instead_of_invalid_content(
    tmp_path: Path,
    monkeypatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    manifest, payload = manifest_for(PLATFORM)
    bundle = write_bundle(tmp_path / "platform.tar", manifest, payload)
    manager.stage(bundle)
    manager.activate(PLATFORM, expected_staged_identity=manifest.identity)
    operator = _Operator(manager, bundle)

    def unavailable(*_args):
        raise RuntimeError("private operator failure")

    monkeypatch.setattr(operator, "observe", unavailable)
    backend = ProductionCommandBackend(
        layout=layout,
        manager=manager,
        operator=operator,
        ingress=_Ingress([]),
        task_factory=lambda: TASK,
    )

    report = backend.doctor()

    checks = {item["check_id"]: item for item in report.value["checks"]}
    assert report.exit_code == 69
    assert checks["deployment-state"]["reason_code"] == "OPERATOR_UNAVAILABLE"
    assert checks["configuration"]["reason_code"] == "OPERATOR_UNAVAILABLE"
    assert checks["database"]["reason_code"] == "OPERATOR_UNAVAILABLE"
    assert checks["frontend"]["reason_code"] == "OPERATOR_UNAVAILABLE"
    assert checks["references"]["reason_code"] == "OPERATOR_UNAVAILABLE"
