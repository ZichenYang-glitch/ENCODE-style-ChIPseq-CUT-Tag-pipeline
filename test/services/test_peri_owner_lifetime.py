"""Native lifetime follows the owner launcher, not a finished request thread."""

from __future__ import annotations

from pathlib import Path
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time

import pytest

from encode_pipeline.persistence import (
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.persistence.authentication import (
    SqlAlchemyAuthenticationRepository,
)
from encode_pipeline.persistence.authoring import SqlAlchemyAuthoringRepository
from encode_pipeline.persistence.data_registry import SqlAlchemyDataRegistryRepository
from encode_pipeline.platform.authoring import AuthoringActor
from encode_pipeline.services.authentication import (
    Argon2idParameters,
    PasswordManager,
    digest_session_token,
)
from encode_pipeline.services.authentication_service import (
    AccountAdministrationService,
    AuthenticationService,
)
from encode_pipeline.services.authoring import AuthoringService
from encode_pipeline.services.authoring_bridge import (
    BridgeDisconnected,
    BridgeTurnRequest,
)
from encode_pipeline.services.data_registry import DataRegistryService
from encode_pipeline.services.defaults import (
    create_default_validation_service,
    create_default_workflow_registry,
)
from encode_pipeline.services.peri_runtime.identity import binary_matches
from encode_pipeline.services.peri_runtime.isolation import native_launcher
from encode_pipeline.services.peri_runtime.limits import (
    cgroup_path_is_qualified,
    cgroup_root_for_pid,
    stop_scope,
)
from encode_pipeline.services.peri_runtime.owner import NativeAuthoringBridge
from encode_pipeline.services.peri_runtime.socket_server import (
    OwnerServer,
    RemoteAuthoringBridge,
)
from encode_pipeline.services.validated_inputs import ValidatedInputService

WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


def _binary() -> Path:
    raw = os.environ.get("HELIXWEAVE_PERI_BINARY", "")
    binary = Path(raw) if raw else None
    if binary is None or not binary_matches(binary):
        pytest.fail(
            "HELIXWEAVE_PERI_BINARY must point at the pinned Peri binary",
            pytrace=False,
        )
    return binary


def _service(root: Path, bridge: RemoteAuthoringBridge):
    database_url = f"sqlite:///{root / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    auth_repository = SqlAlchemyAuthenticationRepository(create_session_factory(engine))
    passwords = PasswordManager(
        Argon2idParameters(memory_cost_kib=19456, time_cost=2, parallelism=1)
    )
    AccountAdministrationService(
        repository=auth_repository,
        password_manager=passwords,
    ).bootstrap_initial_administrator("owner", "correct horse battery staple")
    login = AuthenticationService(
        repository=auth_repository,
        password_manager=passwords,
    ).login("owner", "correct horse battery staple", client_identity="127.0.0.1")
    actor = AuthoringActor(
        user_id=login.principal.user_id,
        session_digest=digest_session_token(login.secrets.session_token),
    )
    registry = create_default_workflow_registry()
    data = DataRegistryService(
        repository=SqlAlchemyDataRegistryRepository(create_session_factory(engine))
    )
    project = data.create_project("Peri owner lifetime")
    repository = SqlAlchemyAuthoringRepository(create_session_factory(engine))
    service = AuthoringService(
        repository=repository,
        registry=registry,
        data_registry=data,
        validated_inputs=ValidatedInputService(
            registry=registry,
            validation_service=create_default_validation_service(registry=registry),
            build_identity_provider=None,
            repository=None,
        ),
        snapshots=None,
        bridge=bridge,
    )
    return engine, actor, project, repository, data, service


def _question(value: int) -> dict[str, object]:
    return {
        "intent": "suggest",
        "field_tokens": ["config.threads"],
        "preferences": [{"field_token": "config.threads", "value": value}],
    }


def _wait_until(predicate, timeout: float = 5) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if predicate():
            return
        time.sleep(0.05)
    assert predicate()


def _pid_gone(pid: int) -> bool:
    return not Path(f"/proc/{pid}").exists()


@pytest.mark.platform_real_execution
def test_socket_open_survives_handler_exit_and_completes_a_turn(tmp_path: Path) -> None:
    binary = _binary()
    root = Path(tempfile.mkdtemp(prefix="hw-c4-life-", dir="/tmp"))
    socket_path = root / "owner.sock"
    state_root = root / "state"
    engine = None
    server = None
    local = None
    try:
        remote = RemoteAuthoringBridge(socket_path)
        engine, actor, project, repository, data, service = _service(root, remote)
        registry = create_default_workflow_registry()
        local = NativeAuthoringBridge(
            binary=binary,
            state_root=state_root,
            repository=repository,
        )
        owner_service = AuthoringService(
            repository=repository,
            registry=registry,
            data_registry=data,
            validated_inputs=ValidatedInputService(
                registry=registry,
                validation_service=create_default_validation_service(registry=registry),
                build_identity_provider=None,
                repository=None,
            ),
            snapshots=None,
            bridge=local,
        )
        local.bind_tool_handler(owner_service.read_tool)
        server = OwnerServer(local, repository, socket_path)
        server.start()
        created = service.create_session(
            actor,
            project_id=project.project_id,
            workflow_id=WORKFLOW_ID,
        )
        handlers: list[threading.Thread] = []
        original_process_request = server._server.process_request  # noqa: SLF001

        def record_handler(request, client_address) -> None:
            thread = threading.Thread(
                target=server._server.process_request_thread,  # noqa: SLF001
                args=(request, client_address),
                daemon=True,
            )
            handlers.append(thread)
            thread.start()

        server._server.process_request = record_handler  # noqa: SLF001
        del original_process_request
        opened = service.open_connection(actor, created.sid)
        assert opened.connection_state == "ready"
        assert handlers
        for handler in handlers:
            handler.join(timeout=2)
            assert not handler.is_alive()
        creator = native_launcher().native_id()
        assert creator is not None
        assert Path(f"/proc/{os.getpid()}/task/{creator}").exists()
        assert creator not in {handler.native_id for handler in handlers}
        diag = local.diagnostics(created.sid)
        leader = int(diag["leader_pid"])
        scope_path = Path(str(diag["scope_cgroup"]))
        time.sleep(0.5)
        assert not _pid_gone(leader)
        native_root = cgroup_root_for_pid(int(diag["native_pid"]))
        assert native_root is not None and cgroup_path_is_qualified(native_root)
        assert native_root != cgroup_root_for_pid(os.getpid())

        turn = service.start_turn(
            actor,
            created.sid,
            client_turn_id="turn-1",
            expected_revision=0,
            question=_question(4),
        )
        assert turn.status == "completed", local.closure_notes
        assert turn.result_kind == "proposal"
        assert turn.proposal is not None
        assert turn.proposal["changes"][0]["value"] == 4
        after = local.diagnostics(created.sid)
        assert after["handler_calls"] >= 2
        assert [item["ok"] for item in after["tool_payloads"]] == [True, True]
        service.close_connection(actor, created.sid)
        _wait_until(lambda: _pid_gone(leader))
        _wait_until(lambda: not scope_path.exists())

        other = service.create_session(
            actor,
            project_id=project.project_id,
            workflow_id=WORKFLOW_ID,
        )
        service.open_connection(actor, other.sid)
        second = local.diagnostics(other.sid)
        with pytest.raises(BridgeDisconnected):
            remote.execute_turn(
                BridgeTurnRequest(
                    sid=other.sid,
                    turn_id="atr_" + "ab" * 16,
                    generation=int(second["generation"]) - 1,
                    instance_id=str(second["instance"]),
                    owner_epoch=remote.owner_epoch(),
                    capability_token="stale",
                    question=_question(2),
                    issue_codes=(),
                    schema_digest="a" * 64,
                    projection_digest="b" * 64,
                    base_revision=0,
                )
            )
        assert local.diagnostics(other.sid)["handler_calls"] == second["handler_calls"]
        closed = remote.drop_authorizations((actor.session_digest,))
        assert closed >= 1
        _wait_until(lambda: _pid_gone(int(second["leader_pid"])))
        revoked = service.start_turn(
            actor,
            other.sid,
            client_turn_id="turn-revoked",
            expected_revision=0,
            question=_question(2),
        )
        assert revoked.status == "failed"
        assert revoked.proposal_id is None

        replacement = service.create_session(
            actor,
            project_id=project.project_id,
            workflow_id=WORKFLOW_ID,
        )
        service.open_connection(actor, replacement.sid)
        live = local.diagnostics(replacement.sid)
        stop_scope(str(live["scope_unit"]))
        for pid in (live.get("leader_pid"), live.get("native_pid")):
            if isinstance(pid, int):
                try:
                    os.kill(pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
        disconnected = service.start_turn(
            actor,
            replacement.sid,
            client_turn_id="turn-disconnected",
            expected_revision=0,
            question=_question(2),
        )
        assert disconnected.status == "failed"
        assert disconnected.proposal_id is None
    finally:
        if server is not None:
            server.shutdown()
        if local is not None:
            local.shutdown()
        if engine is not None:
            engine.dispose()
        shutil.rmtree(root, ignore_errors=True)


def _scope_line(text: str) -> dict[str, str]:
    lines = [
        line for line in text.splitlines() if line.startswith("peri-authoring-scope ")
    ]
    assert lines
    match = re.search(
        r"unit=(\S+) pid=(\d+) path=(\S+) memory_max=(\S+) "
        r"swap_max=(\S+) pids_max=(\S+) cpu_max=(\d+ \d+) owner_outside=(\S+)",
        lines[-1],
    )
    assert match is not None
    return {
        "unit": match.group(1),
        "pid": match.group(2),
        "path": match.group(3),
        "memory_max": match.group(4),
        "swap_max": match.group(5),
        "pids_max": match.group(6),
        "cpu_max": match.group(7),
        "owner_outside": match.group(8),
    }


@pytest.mark.platform_real_execution
def test_owner_process_death_reaps_the_native_and_scope() -> None:
    binary = _binary()
    root = Path(tempfile.mkdtemp(prefix="hw-c4-death-", dir="/tmp"))
    socket_path = root / "owner.sock"
    log_path = root / "owner.log"
    proc = None
    engine = None
    try:
        remote = RemoteAuthoringBridge(socket_path)
        engine, actor, project, _repository, _data, service = _service(root, remote)
        env = os.environ.copy()
        env["HELIXWEAVE_PERI_BINARY"] = str(binary)
        env["HELIXWEAVE_PERI_OWNER_SOCKET"] = str(socket_path)
        env["ENCODE_PIPELINE_DATABASE_URL"] = f"sqlite:///{root / 'platform.db'}"
        log = log_path.open("w", encoding="utf-8")
        proc = subprocess.Popen(
            [
                sys.executable,
                "-m",
                "encode_pipeline.services.peri_runtime.socket_server",
            ],
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
        log.close()
        _wait_until(lambda: socket_path.exists() and proc.poll() is None, timeout=20)
        created = service.create_session(
            actor,
            project_id=project.project_id,
            workflow_id=WORKFLOW_ID,
        )
        opened = service.open_connection(actor, created.sid)
        assert opened.connection_state == "ready"
        scope = _scope_line(log_path.read_text(encoding="utf-8"))
        assert scope["memory_max"] == "536870912"
        assert scope["swap_max"] == "0"
        assert scope["pids_max"] == "128"
        assert scope["cpu_max"] == "100000 100000"
        assert scope["owner_outside"] == "True"
        native_pid = int(scope["pid"])
        assert not _pid_gone(native_pid)
        os.kill(proc.pid, signal.SIGKILL)
        proc.wait(timeout=5)
        _wait_until(lambda: _pid_gone(native_pid))
        _wait_until(lambda: not Path(scope["path"]).exists())
    finally:
        if proc is not None and proc.poll() is None:
            proc.kill()
            proc.wait(timeout=5)
        if engine is not None:
            engine.dispose()
        shutil.rmtree(root, ignore_errors=True)
