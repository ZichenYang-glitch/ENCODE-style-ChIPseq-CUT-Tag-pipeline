"""Fixed Peri bridge: admission, two read tools, and old-generation rejection."""

from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import signal
import tempfile

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
from encode_pipeline.services.data_registry import DataRegistryService
from encode_pipeline.services.defaults import (
    create_default_validation_service,
    create_default_workflow_registry,
)
from encode_pipeline.services.peri_runtime.identity import binary_matches
from encode_pipeline.services.peri_runtime.limits import (
    cgroup_path_is_qualified,
    cgroup_root_for_pid,
    stop_scope,
)
from encode_pipeline.services.peri_runtime.owner import NativeAuthoringBridge
from encode_pipeline.services.validated_inputs import ValidatedInputService


NOW = datetime(2026, 11, 10, 12, 0, tzinfo=timezone.utc)


def _provider_contains_tool_result(requests: object) -> bool:
    if isinstance(requests, str):
        return '"ok": true' in requests or '"ok":true' in requests
    if isinstance(requests, dict):
        return any(_provider_contains_tool_result(item) for item in requests.values())
    if isinstance(requests, list):
        return any(_provider_contains_tool_result(item) for item in requests)
    return False


def _pinned_binary() -> Path:
    raw = os.environ.get("HELIXWEAVE_PERI_BINARY", "")
    binary = Path(raw) if raw else None
    if binary is None or not binary_matches(binary):
        pytest.fail(
            "HELIXWEAVE_PERI_BINARY must point at the pinned Peri binary",
            pytrace=False,
        )
    return binary


def test_unverified_binary_stays_unavailable(tmp_path: Path) -> None:
    binary = tmp_path / "peri"
    binary.write_bytes(b"not-the-fixed-runtime")
    bridge = NativeAuthoringBridge(binary=binary, state_root=tmp_path / "state")

    assert bridge.status().available is False


@pytest.mark.platform_real_execution
def test_fixed_native_reads_two_tools_and_rejects_an_old_generation(
    tmp_path: Path,
) -> None:
    fixed_binary = _pinned_binary()
    database_url = f"sqlite:///{tmp_path / 'authoring.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)

    def clock() -> datetime:
        return NOW

    auth_repository = SqlAlchemyAuthenticationRepository(create_session_factory(engine))
    passwords = PasswordManager(
        Argon2idParameters(memory_cost_kib=19456, time_cost=2, parallelism=1)
    )
    AccountAdministrationService(
        repository=auth_repository,
        password_manager=passwords,
        now_factory=clock,
    ).bootstrap_initial_administrator("owner", "correct horse battery staple")
    login = AuthenticationService(
        repository=auth_repository,
        password_manager=passwords,
        now_factory=clock,
    ).login("owner", "correct horse battery staple", client_identity="127.0.0.1")
    actor = AuthoringActor(
        user_id=login.principal.user_id,
        session_digest=digest_session_token(login.secrets.session_token),
    )
    registry = create_default_workflow_registry()
    data = DataRegistryService(
        repository=SqlAlchemyDataRegistryRepository(create_session_factory(engine)),
        now_factory=clock,
    )
    project = data.create_project("Peri authoring")
    authoring_repository = SqlAlchemyAuthoringRepository(
        create_session_factory(engine),
        clock=clock,
    )
    state_root = Path(tempfile.mkdtemp(prefix="hw-a1-native-"))
    bridge = NativeAuthoringBridge(
        binary=fixed_binary,
        state_root=state_root,
        repository=authoring_repository,
    )
    service = AuthoringService(
        repository=authoring_repository,
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
        clock=clock,
    )
    bridge.bind_tool_handler(service.read_tool)
    try:
        assert bridge.status().available is True
        created = service.create_session(
            actor,
            project_id=project.project_id,
            workflow_id="encode-style-chipseq-cuttag-atac-mnase",
        )
        opened = service.open_connection(actor, created.sid)
        first = bridge.diagnostics(created.sid)
        assert opened.connection_state == "ready", bridge.closure_notes
        assert opened.requires_new_connection is False
        assert first["handler_calls"] == 0
        assert first["approvals"]
        assert not any(item["allowed"] for item in first["approvals"])
        assert "bwrap" in " ".join(first["argv"])
        turn = service.start_turn(
            actor,
            created.sid,
            client_turn_id="turn-1",
            expected_revision=0,
            question={
                "intent": "suggest",
                "field_tokens": ["config.threads"],
                "preferences": [{"field_token": "config.threads", "value": 4}],
            },
        )
        assert turn.status == "completed", (turn, bridge.closure_notes)
        assert turn.result_kind == "proposal"
        assert turn.proposal is not None
        assert turn.proposal["changes"][0]["field_token"] == "config.threads"
        assert turn.proposal["changes"][0]["value"] == 4
        after = bridge.diagnostics(created.sid)
        assert after["tool_arrivals"] >= 2
        assert after["handler_calls"] >= 2
        assert "tool_handler_failed" not in after["denied_callbacks"]
        payloads = after["tool_payloads"]
        assert [item["ok"] for item in payloads] == [True, True]
        assert any(
            field.get("token") == "config.threads"
            for field in payloads[0]["data"]["fields"]
        )
        assert any(
            field.get("token") == "config.threads"
            for field in payloads[1]["data"]["fields"]
        )
        provider_bodies = json.dumps(after["provider_requests"])
        assert "config.threads" in provider_bodies
        assert _provider_contains_tool_result(after["provider_requests"])
        native_root = cgroup_root_for_pid(int(after["native_pid"]))
        assert native_root is not None and cgroup_path_is_qualified(native_root)
        assert native_root != cgroup_root_for_pid(os.getpid())
        assert bridge.sweep() == 0
        first_pid = after["pid"]

        other = service.create_session(
            actor,
            project_id=project.project_id,
            workflow_id="encode-style-chipseq-cuttag-atac-mnase",
        )
        service.open_connection(actor, other.sid)
        second = bridge.diagnostics(other.sid)
        assert second["pid"] != first_pid
        assert second["instance"] != after["instance"]
        from encode_pipeline.services.authoring_bridge import (
            BridgeDisconnected,
            BridgeTurnRequest,
        )

        with pytest.raises(BridgeDisconnected):
            bridge.execute_turn(
                BridgeTurnRequest(
                    sid=other.sid,
                    turn_id="atr_" + "ab" * 16,
                    generation=int(second["generation"]) - 1,
                    instance_id=str(second["instance"]),
                    owner_epoch=bridge.owner_epoch(),
                    capability_token="stale",
                    question={
                        "intent": "suggest",
                        "field_tokens": [],
                        "preferences": [],
                    },
                    issue_codes=(),
                    schema_digest="a" * 64,
                    projection_digest="b" * 64,
                    base_revision=0,
                )
            )
        assert bridge.diagnostics(other.sid)["handler_calls"] == second["handler_calls"]
        assert "Maximum CPU cores" not in provider_bodies
        assert project.project_id not in provider_bodies
        original_read = authoring_repository.tool_read_state

        def deny_egress(**kwargs):
            del kwargs
            return "AUTH_SESSION_ENDED", None, None

        authoring_repository.tool_read_state = deny_egress
        denied = service.start_turn(
            actor,
            created.sid,
            client_turn_id="turn-denied",
            expected_revision=0,
            question={
                "intent": "suggest",
                "field_tokens": ["config.threads"],
                "preferences": [{"field_token": "config.threads", "value": 2}],
            },
        )
        authoring_repository.tool_read_state = original_read
        assert denied.status == "failed"
        assert denied.proposal_id is None
        denied_diag = bridge.closure_notes.get(
            created.sid, bridge.diagnostics(created.sid)
        )
        assert denied_diag["handler_calls"] == after["handler_calls"]
        assert "authorization_revoked" in denied_diag["stub_failures"]
        stop_scope(str(second["scope_unit"]))
        for pid in (second.get("leader_pid"), second.get("native_pid")):
            if not isinstance(pid, int):
                continue
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        disconnected = service.start_turn(
            actor,
            other.sid,
            client_turn_id="turn-disconnected",
            expected_revision=0,
            question={
                "intent": "suggest",
                "field_tokens": ["config.threads"],
                "preferences": [{"field_token": "config.threads", "value": 2}],
            },
        )
        assert disconnected.status == "failed", bridge.closure_notes
        assert disconnected.proposal_id is None
        observed = bridge.closure_notes.get(other.sid, bridge.diagnostics(other.sid))
        assert observed["handler_calls"] == 0
    finally:
        bridge.shutdown()
        shutil.rmtree(state_root, ignore_errors=True)
        engine.dispose()
