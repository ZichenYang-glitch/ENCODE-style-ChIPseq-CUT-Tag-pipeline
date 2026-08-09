from __future__ import annotations

from dataclasses import dataclass
import importlib.util
import json
import os
from pathlib import Path
import runpy
import shutil
import stat
import subprocess
import sys
from types import SimpleNamespace

import pytest

from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.operator import (
    OperatorOutcome,
    OperatorRequest,
    SERVICE_UNITS,
    ServiceIdentity,
    SocketWitness,
    bundle_ingress_path,
    execute_request,
    parse_request,
)


ROOT = Path(__file__).resolve().parents[2]
TEMPLATES = ROOT / "src" / "encode_pipeline" / "deployment" / "templates"
BOOTSTRAP_SCRIPT = ROOT / "scripts" / "bootstrap_helixweave_operator.py"
IDENTITY = f"sha256-{'a' * 64}"
SERVICE_IDENTITY = f"sha256-{'b' * 64}"
TASK_IDENTITY = f"task-{'c' * 32}"


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
                SERVICE_IDENTITY,
            ),
            "status",
            None,
            "helixweave-worker.service",
            SERVICE_IDENTITY,
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
        ("status", "helixweave-api.service", IDENTITY, TASK_IDENTITY),
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

    assert backend.bundle_path == (
        layout.ingress / "platform" / IDENTITY / "bundle.tar"
    )
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


def test_status_must_return_the_exact_requested_service_identity() -> None:
    service = _worker_service_identity()
    backend = RecordingBackend(OperatorOutcome("running", service=service))

    with pytest.raises(DeploymentError) as caught:
        execute_request(
            (
                "status",
                "helixweave-worker.service",
                IDENTITY,
                TASK_IDENTITY,
                SERVICE_IDENTITY,
            ),
            backend=backend,
        )

    assert caught.value.issue.code == "OPERATOR_RESULT_INVALID"


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
    api_socket = SocketWitness(name="api-http", device=1, inode=2)

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


def test_only_supported_units_are_exposed() -> None:
    assert SERVICE_UNITS == (
        "helixweave-api.service",
        "helixweave-worker.service",
        "helixweave-redis.service",
        "helixweave-docker-rootless.service",
    )
    assert "docker.service" not in SERVICE_UNITS


def test_templates_encode_one_bounded_hybrid_topology() -> None:
    api = (TEMPLATES / "helixweave-api.service.in").read_text()
    worker = (TEMPLATES / "helixweave-worker.service.in").read_text()
    redis = (TEMPLATES / "helixweave-redis.service").read_text()
    docker = (TEMPLATES / "helixweave-docker-rootless.service").read_text()
    target = (TEMPLATES / "helixweave.target").read_text()

    for unit in (api, worker, redis, docker):
        assert "User=helixweave" in unit
        assert "ExecStart=/bin/sh" not in unit
        assert "ExecStart=/bin/bash" not in unit
        assert "docker compose" not in unit.lower()
    assert "uvicorn" in api
    assert "encode_pipeline.deployment.frontend:create_app" in api
    assert "encode_pipeline.api.main:create_app" not in api
    assert "static frontend" in api
    assert "npm" not in api
    assert "vite" not in api.lower()
    assert "encode_pipeline.workers.cli" in worker
    assert "RestrictNamespaces=" not in worker
    assert "port 0" in (TEMPLATES / "redis.conf").read_text()
    assert "dockerd-rootless.sh" in docker
    assert "Type=simple" in docker
    assert "Delegate=yes" in docker
    assert "unix:///run/helixweave/docker/docker.sock" in docker
    assert "--data-root=/var/lib/helixweave/docker-rootless" in docker
    assert "/var/run/docker.sock" not in docker
    assert "Wants=helixweave-docker-rootless.service" in target
    assert not any(
        path.name == "helixweave-frontend.service" for path in TEMPLATES.iterdir()
    )


def test_sudoers_and_helper_templates_do_not_expose_a_shell_or_environment() -> None:
    policy = (TEMPLATES / "helixweave-operator.sudoers").read_text()
    helper = (TEMPLATES / "helixweave-operator").read_text()
    gate_cleanup = (TEMPLATES / "helixweave-gate-cleanup").read_text()
    gate_launcher = (TEMPLATES / "helixweave-gate-cleanup-launcher").read_text()
    launcher = (TEMPLATES / "helixweave-operator-launcher").read_text()

    assert policy.count("/usr/libexec/helixweave-operator") == 8
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
    assert "encode_pipeline" not in helper
    assert "os.execve(" in launcher
    assert "/opt/helixweave/operator/current/bin/helixweave-operator" in launcher
    assert "encode_pipeline" not in launcher
    assert gate_cleanup.startswith("#!/usr/bin/python3 -I\n")
    assert "subprocess" not in gate_cleanup
    assert "encode_pipeline" not in gate_cleanup
    assert "GATE_CLEANUP_BACKEND_UNAVAILABLE" in gate_cleanup
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


def test_gate_cleanup_helper_consumes_every_bound_identity(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    namespace = runpy.run_path(str(TEMPLATES / "helixweave-gate-cleanup"))
    main = namespace["main"]
    monkeypatch.setattr(namespace["os"], "geteuid", lambda: 0)
    main.__globals__["_installed_identities"] = lambda: (
        SERVICE_IDENTITY,
        IDENTITY,
        SERVICE_IDENTITY,
    )
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

    assert main(values) == 69
    assert json.loads(capsys.readouterr().err)["issue"]["code"] == (
        "GATE_CLEANUP_BACKEND_UNAVAILABLE"
    )
    changed = (*values[:6], IDENTITY, *values[7:])
    assert main(changed) == 65
    assert json.loads(capsys.readouterr().err)["issue"]["code"] == (
        "GATE_CLEANUP_EXECUTOR_MISMATCH"
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


def test_deployed_operator_helper_matches_the_redacted_receipt_contract() -> None:
    script = TEMPLATES / "helixweave-operator"
    code = """
import json
import runpy
import sys

values = runpy.run_path(sys.argv[1])
values["main"].__globals__["os"].geteuid = lambda: 0
raise SystemExit(values["main"]((
    "stage", "platform", "sha256-" + "a" * 64, "task-" + "c" * 32
)))
"""
    completed = subprocess.run(
        (sys.executable, "-I", "-S", "-c", code, str(script)),
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 69
    receipt = json.loads(completed.stderr)
    assert receipt == {
        "schema_version": "helixweave-operator-receipt-v1",
        "status": "error",
        "issue": {
            "code": "OPERATOR_OPERATION_FAILED",
            "message": "Operator action failed.",
            "recoverable": True,
        },
    }


def test_tmpfiles_separates_immutable_data_and_external_references() -> None:
    content = (TEMPLATES / "helixweave.tmpfiles.conf").read_text()

    assert "/var/lib/helixweave 0755 root root" in content
    assert "/run/helixweave 0755 root root" in content
    assert "/opt/helixweave/releases/platform 0555 root root" in content
    assert "/opt/helixweave/runtimes/bulk-rnaseq 0555 root root" in content
    assert "/var/lib/helixweave/database/backups 0700 helixweave helixweave" in content
    assert "/var/lib/helixweave/workspaces 0750 helixweave helixweave" in content
    assert "/operator/ingress/platform 2730 root helixweave-operators" in content
    assert "reference" not in content.lower()

    for unit_name in ("helixweave-api.service.in", "helixweave-worker.service.in"):
        unit = (TEMPLATES / unit_name).read_text()
        assert "ReadWritePaths=/var/lib/helixweave /run/helixweave" not in unit
        assert "InaccessiblePaths=/var/lib/helixweave/operator" in unit


def test_platform_environment_uses_existing_runtime_contract_names() -> None:
    content = (TEMPLATES / "helixweave-platform.env.in").read_text()

    assert "ENCODE_PIPELINE_DATABASE_URL=" in content
    assert (
        "ENCODE_PIPELINE_REDIS_URL=unix:///run/helixweave/redis/redis.sock" in content
    )
    assert "ENCODE_PIPELINE_REFERENCE_PROFILE_CONFIG=" in content
    assert "HELIXWEAVE_BULK_RNASEQ_RUNTIME_ROOT=@BULK_RNASEQ_RUNTIME_ROOT@" in content
    assert (
        "ENCODE_PIPELINE_MANAGED_DOCKER_SOCKET=/run/helixweave/docker/docker.sock"
        in content
    )
    assert "/var/run/docker.sock" not in content


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
    assert (closure / "closure.json").is_file()
    assert len(validator_calls) == 1
    assert not list(tmp_path.rglob("*.tmp"))
    assert not list(tmp_path.rglob(".partial-*"))


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

    monkeypatch.setattr(backend, "_switch_closure", original_switch)
    backend.apply(operation="update", invoking_user="labadmin")
    second = os.readlink(current)
    assert second != first
    assert (current.parent / second / "closure.json").is_file()


def test_bootstrap_update_refuses_a_stable_boundary_abi_change(
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

    with pytest.raises(bootstrap.BootstrapFailure) as caught:
        backend.apply(operation="update", invoking_user="labadmin")

    assert caught.value.code == "BOOTSTRAP_BOUNDARY_UPDATE_REQUIRED"
    assert os.readlink(current) == original
    with pytest.raises(bootstrap.BootstrapFailure) as reinstall:
        backend.apply(operation="install", invoking_user="labadmin")
    assert reinstall.value.code == "BOOTSTRAP_ALREADY_INSTALLED"
