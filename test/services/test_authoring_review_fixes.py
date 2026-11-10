"""Lowest-layer checks for the reviewed authoring contract gaps."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path
from types import SimpleNamespace
import json
import socket
import threading

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
from encode_pipeline.persistence.models import AuthSessionRow
from encode_pipeline.platform.authentication import SessionRevocationReason
from encode_pipeline.persistence.authoring import TurnCompletion
from encode_pipeline.platform.authoring import (
    ADMISSION_FAILED,
    AUTHENTICATION_REQUIRED,
    CONNECTION_LOST,
    PROPOSAL_INVALIDATED,
    AuthoringActor,
    AuthoringRejected,
    draft_digest,
    empty_draft_content,
    proposal_digest,
)
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
    BridgeAdmission,
    BridgeDisconnected,
    BridgeSessionRequest,
    BridgeStatus,
    BridgeTurnRequest,
)
from encode_pipeline.services.data_registry import DataRegistryService
from encode_pipeline.services.defaults import create_default_workflow_registry
from encode_pipeline.services.peri_runtime import owner as owner_module
from encode_pipeline.services.peri_runtime.owner import NativeAuthoringBridge
from encode_pipeline.services.peri_runtime.session import NativeRejected
from encode_pipeline.services.peri_runtime.socket_server import (
    OwnerServer,
    RemoteAuthoringBridge,
)


NOW = datetime(2026, 11, 10, 12, 0, tzinfo=timezone.utc)
WORKFLOW = "encode-style-chipseq-cuttag-atac-mnase"
FAST = Argon2idParameters(memory_cost_kib=19456, time_cost=2, parallelism=1)


class Clock:
    def __init__(self) -> None:
        self.now = NOW

    def __call__(self) -> datetime:
        return self.now


class ScriptedBridge:
    def __init__(self) -> None:
        self.epoch = "owner-epoch-0001"
        self.closes: list[tuple[int, str | None]] = []
        self.answer = "{}"
        self._admissions = 0
        self.tokens: dict[tuple[str, int], str] = {}

    def status(self) -> BridgeStatus:
        return BridgeStatus(available=True)

    def owner_epoch(self) -> str:
        return self.epoch

    def admit(self, request) -> BridgeAdmission:
        self.tokens[(request.sid, request.generation)] = request.capability_token
        self._admissions += 1
        return BridgeAdmission(
            instance_id=f"instance{self._admissions:02d}",
            state="ready",
        )

    def execute_turn(self, request) -> str:
        if self.answer == "proposal":
            body = {
                "base_revision": request.base_revision,
                "changes": [{"field_token": "config.threads", "value": 4}],
                "kind": "proposal",
                "reason_codes": [],
                "schema_digest": request.schema_digest,
            }
        else:
            body = {
                "base_revision": request.base_revision,
                "changes": [],
                "kind": "explanation",
                "reason_codes": ["ENCODE_OPTIONS_INVALID"],
                "schema_digest": request.schema_digest,
            }
        return json.dumps(body, separators=(",", ":"), sort_keys=True)

    def close(self, request) -> None:
        self.closes.append((request.generation, request.instance_id))


class RecordingBridge:
    def __init__(self, epoch: str) -> None:
        self.epoch = epoch
        self.sessions: dict[str, dict[str, object]] = {}
        self.turns: list[str] = []
        self.closes: list[int] = []

    def status(self) -> BridgeStatus:
        return BridgeStatus(available=True)

    def owner_epoch(self) -> str:
        return self.epoch

    def admit(self, request) -> BridgeAdmission:
        if request.owner_epoch != self.epoch:
            raise BridgeDisconnected
        instance = f"inst{len(self.sessions) + 1:06d}"
        self.sessions[request.sid] = {
            "generation": request.generation,
            "instance": instance,
        }
        return BridgeAdmission(instance_id=instance, state="ready")

    def execute_turn(self, request) -> str:
        self.turns.append(request.owner_epoch)
        current = self.sessions.get(request.sid)
        if (
            request.owner_epoch != self.epoch
            or current is None
            or current["generation"] != request.generation
        ):
            raise BridgeDisconnected
        return (
            '{"kind":"explanation","changes":[],"reason_codes":[],"schema_digest":"%s","base_revision":%s}'
            % (
                request.schema_digest,
                request.base_revision,
            )
        )

    def close(self, request) -> None:
        current = self.sessions.get(request.sid)
        if current is None or current["generation"] != request.generation:
            return
        if request.instance_id not in {None, current["instance"]}:
            return
        self.closes.append(request.generation)
        self.sessions.pop(request.sid, None)

    def drop_authorizations(self, digests) -> int:
        del digests
        return 0

    def sweep(self) -> int:
        return 0

    def shutdown(self) -> None:
        self.sessions.clear()


@pytest.fixture
def clock() -> Clock:
    return Clock()


@pytest.fixture
def stack(tmp_path, clock):
    database_url = f"sqlite:///{tmp_path / 'authoring.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    auth_repository = SqlAlchemyAuthenticationRepository(create_session_factory(engine))
    passwords = PasswordManager(FAST)
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
    data = DataRegistryService(
        repository=SqlAlchemyDataRegistryRepository(create_session_factory(engine)),
        now_factory=clock,
    )
    project = data.create_project("Review")
    repository = SqlAlchemyAuthoringRepository(
        create_session_factory(engine),
        clock=clock,
    )
    registry = create_default_workflow_registry()
    bridge = ScriptedBridge()
    service = AuthoringService(
        repository=repository,
        registry=registry,
        data_registry=data,
        validated_inputs=SimpleNamespace(),
        snapshots=SimpleNamespace(),
        bridge=bridge,
        clock=clock,
    )
    yield SimpleNamespace(
        actor=actor,
        auth=auth_repository,
        bridge=bridge,
        clock=clock,
        database_url=database_url,
        engine=engine,
        project_id=project.project_id,
        repository=repository,
        service=service,
    )
    engine.dispose()


def _extend_login(repository, actor, until: datetime) -> None:
    with repository._session_factory.begin() as session:
        row = session.get(AuthSessionRow, actor.session_digest)
        assert row is not None
        row.expires_at = until


def _ready(stack) -> tuple[str, str]:
    created = stack.service.create_session(
        stack.actor,
        project_id=stack.project_id,
        workflow_id=WORKFLOW,
    )
    opened = stack.service.open_connection(stack.actor, created.sid)
    assert opened.connection_state == "ready"
    record = stack.repository.get_session(stack.actor, created.sid)
    token = stack.bridge.tokens[(created.sid, record.generation)]
    return created.sid, token


def test_read_tools_return_projected_fields(stack) -> None:
    sid, token = _ready(stack)
    record = stack.repository.get_session(stack.actor, sid)
    turn, _replayed = stack.repository.start_turn(
        stack.actor,
        sid,
        client_turn_id="read-1",
        request_digest="ab" * 32,
        question={
            "intent": "suggest",
            "field_tokens": ["config.threads"],
            "preferences": [],
        },
        expected_revision=record.revision,
        schema_digest=record.schema_digest,
        projection_digest=record.projection_digest,
        owner_epoch=record.bridge_owner_epoch,
    )
    schema = stack.service.read_tool(
        sid=sid,
        turn_id=turn.turn_id,
        generation=record.generation,
        instance_id=record.instance_id,
        owner_epoch=record.bridge_owner_epoch,
        capability_token=token,
        name="get_authoring_schema",
    )
    draft = stack.service.read_tool(
        sid=sid,
        turn_id=turn.turn_id,
        generation=record.generation,
        instance_id=record.instance_id,
        owner_epoch=record.bridge_owner_epoch,
        capability_token=token,
        name="get_draft_summary",
    )

    assert schema["ok"] is True
    assert {field["token"] for field in schema["data"]["fields"]} == {
        "config.threads",
        "options.cores",
    }
    assert "description" not in json.dumps(schema)
    assert draft["ok"] is True
    assert draft["data"]["revision"] == record.revision
    assert {field["token"] for field in draft["data"]["fields"]} == {
        "config.threads",
        "options.cores",
    }


def test_expired_and_revoked_reads_return_no_projection(stack) -> None:
    sid, token = _ready(stack)
    record = stack.repository.get_session(stack.actor, sid)
    turn, _replayed = stack.repository.start_turn(
        stack.actor,
        sid,
        client_turn_id="read-2",
        request_digest="cd" * 32,
        question={
            "intent": "suggest",
            "field_tokens": ["config.threads"],
            "preferences": [],
        },
        expected_revision=record.revision,
        schema_digest=record.schema_digest,
        projection_digest=record.projection_digest,
        owner_epoch=record.bridge_owner_epoch,
    )
    _extend_login(
        stack.repository, stack.actor, record.draft_expires_at + timedelta(days=1)
    )
    stack.clock.now = record.draft_expires_at
    expired = stack.service.read_tool(
        sid=sid,
        turn_id=turn.turn_id,
        generation=record.generation,
        instance_id=record.instance_id,
        owner_epoch=record.bridge_owner_epoch,
        capability_token=token,
        name="get_draft_summary",
    )
    assert expired["ok"] is False
    assert expired["data"] is None
    assert "config.threads" not in json.dumps(expired)

    stack.clock.now = NOW
    sid, token = _ready(stack)
    record = stack.repository.get_session(stack.actor, sid)
    turn, _replayed = stack.repository.start_turn(
        stack.actor,
        sid,
        client_turn_id="read-3",
        request_digest="ef" * 32,
        question={
            "intent": "suggest",
            "field_tokens": ["options.cores"],
            "preferences": [],
        },
        expected_revision=record.revision,
        schema_digest=record.schema_digest,
        projection_digest=record.projection_digest,
        owner_epoch=record.bridge_owner_epoch,
    )
    assert stack.auth.revoke_session(
        stack.actor.session_digest,
        stack.clock.now,
        SessionRevocationReason.LOGOUT,
    )
    denied = stack.service.read_tool(
        sid=sid,
        turn_id=turn.turn_id,
        generation=record.generation,
        instance_id=record.instance_id,
        owner_epoch=record.bridge_owner_epoch,
        capability_token=token,
        name="get_authoring_schema",
    )
    assert denied["ok"] is False
    assert denied["data"] is None


def test_explanation_persists_approved_reason_codes(stack) -> None:
    created = stack.service.create_session(
        stack.actor,
        project_id=stack.project_id,
        workflow_id=WORKFLOW,
    )
    schema_digest, projection_digest, _fields = stack.service._identity(WORKFLOW)
    stack.repository.store_validation(
        stack.actor,
        created.sid,
        expected_revision=0,
        schema_digest=schema_digest,
        projection_digest=projection_digest,
        snapshot_id=None,
        issue_codes=("ENCODE_OPTIONS_INVALID",),
        issues=(),
    )
    stack.service.open_connection(stack.actor, created.sid)
    view = stack.service.start_turn(
        stack.actor,
        created.sid,
        client_turn_id="explain-1",
        expected_revision=0,
        question={
            "intent": "explain_issues",
            "field_tokens": ["options.cores"],
            "preferences": [],
        },
    )

    assert view.status == "completed"
    assert view.result_kind == "explanation"
    assert view.proposal_id is None
    assert view.proposal is None
    assert view.adoptable is False
    assert view.explanation == {"reason_codes": ["ENCODE_OPTIONS_INVALID"]}
    stored = stack.repository.peek_turn(created.sid, view.turn_id)
    assert stored is not None
    assert stored.proposal == {"reason_codes": ["ENCODE_OPTIONS_INVALID"]}
    assert stored.decision is None


def test_schema_change_blocks_completion_adoption_and_cached_snapshot(stack) -> None:
    sid, _token = _ready(stack)
    stack.bridge.answer = "proposal"
    turn = stack.service.start_turn(
        stack.actor,
        sid,
        client_turn_id="suggest-1",
        expected_revision=0,
        question={
            "intent": "suggest",
            "field_tokens": ["config.threads"],
            "preferences": [{"field_token": "config.threads", "value": 4}],
        },
    )
    assert turn.adoptable is True
    assert turn.proposal_digest
    original = stack.service._identity

    def changed(workflow_id: str):
        schema_digest, projection_digest, fields = original(workflow_id)
        return "c" * 64, projection_digest, fields

    stack.service._identity = changed
    visible = stack.service.get_turn(stack.actor, sid, turn.turn_id)
    assert visible.adoptable is False
    assert visible.reason_code == PROPOSAL_INVALIDATED
    with pytest.raises(AuthoringRejected) as captured:
        stack.service.adopt(
            stack.actor,
            sid,
            turn.proposal_id,
            expected_revision=0,
            proposal_digest_value=turn.proposal_digest,
        )
    assert captured.value.reason_code == PROPOSAL_INVALIDATED
    stack.service._identity = original
    adopted = stack.service.adopt(
        stack.actor,
        sid,
        turn.proposal_id,
        expected_revision=0,
        proposal_digest_value=turn.proposal_digest,
    )
    stack.service._identity = changed
    replayed = stack.service.adopt(
        stack.actor,
        sid,
        turn.proposal_id,
        expected_revision=0,
        proposal_digest_value=turn.proposal_digest,
    )
    assert adopted.replayed is False
    assert replayed.replayed is True
    assert replayed.receipt["adopted_at"] == adopted.receipt["adopted_at"]


def test_schema_only_save_drops_the_old_snapshot_cache(stack) -> None:
    class Schema:
        schema_version = "1.0.0"
        schema_dialect = "https://json-schema.org/draft/2020-12/schema"
        marker = "before"

        def to_dict(self):
            return {"marker": self.marker, "schema_version": self.schema_version}

    schema = Schema()
    adapter = SimpleNamespace(
        metadata=SimpleNamespace(version="adapter-1"),
        schema=lambda: schema,
    )
    snapshot = SimpleNamespace(
        snapshot_id="vsnap_" + "ab" * 16,
        workflow_id="plain-workflow",
        schema_version="1.0.0",
        schema_dialect=schema.schema_dialect,
        adapter_version="adapter-1",
        expires_at=NOW + timedelta(days=1),
    )
    stack.service._registry = SimpleNamespace(get=lambda workflow_id: adapter)
    stack.service._snapshots = SimpleNamespace(
        get_validated_input_snapshot=lambda snapshot_id: snapshot
    )
    created = stack.repository.create_session(
        stack.actor,
        project_id=stack.project_id,
        workflow_id="plain-workflow",
        schema_digest=stack.service._identity("plain-workflow")[0],
        projection_digest=None,
    )
    stack.repository.store_validation(
        stack.actor,
        created.sid,
        expected_revision=0,
        schema_digest=created.schema_digest,
        projection_digest=None,
        snapshot_id=snapshot.snapshot_id,
        issue_codes=(),
        issues=(),
    )
    schema.marker = "after"
    saved = stack.service.save_draft(
        stack.actor,
        created.sid,
        expected_revision=0,
        content=empty_draft_content(),
    )
    record = stack.repository.get_session(stack.actor, created.sid)
    assert saved.changed is False
    assert record.revision == 0
    assert record.last_snapshot_id is None
    assert record.validated_schema_digest is None
    assert stack.service._cached_validation(record) is None
    assert snapshot.schema_version == "1.0.0"
    assert snapshot.adapter_version == "adapter-1"

    stack.repository.store_validation(
        stack.actor,
        created.sid,
        expected_revision=0,
        schema_digest=record.schema_digest,
        projection_digest=None,
        snapshot_id=snapshot.snapshot_id,
        issue_codes=(),
        issues=(),
    )
    adapter.metadata.version = "adapter-2"
    current = stack.repository.get_session(stack.actor, created.sid)
    assert stack.service._cached_validation(current) is None
    assert snapshot.adapter_version == "adapter-1"


def test_late_close_from_another_connection_keeps_the_new_native(stack) -> None:
    sid, _token = _ready(stack)
    before = stack.repository.get_session(stack.actor, sid)
    other_engine = create_database_engine(stack.database_url)
    other_repository = SqlAlchemyAuthoringRepository(
        create_session_factory(other_engine),
        clock=stack.clock,
    )
    other = AuthoringService(
        repository=other_repository,
        registry=stack.service._registry,
        data_registry=stack.service._data,
        validated_inputs=SimpleNamespace(),
        snapshots=SimpleNamespace(),
        bridge=stack.bridge,
        clock=stack.clock,
    )
    other.close_connection(stack.actor, sid)
    reopened = other.open_connection(stack.actor, sid)
    late = stack.repository.close_connection(
        stack.actor,
        sid,
        expected_generation=before.generation,
        expected_instance_id=before.instance_id,
        expected_owner_epoch=before.bridge_owner_epoch,
    )
    current = other_repository.get_session(stack.actor, sid)

    assert reopened.connection_state == "ready"
    assert late.applied is False
    assert current.connection_state == "ready"
    assert current.generation != before.generation
    assert current.instance_id != before.instance_id
    assert stack.bridge.closes == [(before.generation, before.instance_id)]
    other_engine.dispose()


def test_one_owner_serves_two_callers_and_fences_before_listen(stack, tmp_path) -> None:
    old = stack.repository.create_session(
        stack.actor,
        project_id=stack.project_id,
        workflow_id=WORKFLOW,
        schema_digest="a" * 64,
        projection_digest="b" * 64,
    )
    claim = stack.repository.open_connection(
        stack.actor,
        old.sid,
        owner_epoch="previous-owner",
        projection_digest="b" * 64,
    )
    stack.repository.mark_connection_ready(
        stack.actor,
        old.sid,
        generation=claim.record.generation,
        instance_id="oldinst01",
        owner_epoch="previous-owner",
    )
    recording = RecordingBridge("current-owner")
    socket_path = tmp_path / "owner.sock"
    seen = {}

    def before_listen() -> None:
        seen["socket"] = socket_path.exists()
        seen["state"] = stack.repository.get_session(
            stack.actor, old.sid
        ).connection_state

    server = OwnerServer(
        recording,
        stack.repository,
        socket_path,
        before_listen=before_listen,
    )
    server.start()
    try:
        assert seen == {"socket": False, "state": "closed"}
        assert server.retired == 1
        remote = RemoteAuthoringBridge(socket_path)
        with pytest.raises(BridgeDisconnected):
            remote.execute_turn(
                BridgeTurnRequest(
                    sid=old.sid,
                    turn_id="atr_" + "ab" * 16,
                    generation=claim.record.generation,
                    instance_id="oldinst01",
                    owner_epoch="previous-owner",
                    capability_token="stale",
                    question={},
                    issue_codes=(),
                    schema_digest="a" * 64,
                    projection_digest="b" * 64,
                    base_revision=0,
                )
            )
        assert recording.turns == []
        left = AuthoringService(
            repository=stack.repository,
            registry=stack.service._registry,
            data_registry=stack.service._data,
            validated_inputs=SimpleNamespace(),
            snapshots=SimpleNamespace(),
            bridge=RemoteAuthoringBridge(socket_path),
            clock=stack.clock,
        )
        right_engine = create_database_engine(stack.database_url)
        right_repository = SqlAlchemyAuthoringRepository(
            create_session_factory(right_engine),
            clock=stack.clock,
        )
        right = AuthoringService(
            repository=right_repository,
            registry=stack.service._registry,
            data_registry=stack.service._data,
            validated_inputs=SimpleNamespace(),
            snapshots=SimpleNamespace(),
            bridge=RemoteAuthoringBridge(socket_path),
            clock=stack.clock,
        )
        first = left.create_session(
            stack.actor, project_id=stack.project_id, workflow_id=WORKFLOW
        )
        second = right.create_session(
            stack.actor, project_id=stack.project_id, workflow_id=WORKFLOW
        )
        assert left.open_connection(stack.actor, first.sid).connection_state == "ready"
        assert (
            right.open_connection(stack.actor, second.sid).connection_state == "ready"
        )
        assert (
            left._bridge.owner_epoch() == right._bridge.owner_epoch() == "current-owner"
        )
        assert set(recording.sessions) == {first.sid, second.sid}
        right_engine.dispose()
    finally:
        server.shutdown()


class _IdleNative:
    def __init__(self, binary, directory, *, service_generation, tool_handler) -> None:
        del binary, tool_handler
        directory.mkdir(parents=True, exist_ok=True)
        self.directory = directory
        self.service_generation = service_generation
        self.instance = "instance01"
        self.state = "blocked"
        self.handler_calls = 0
        self.tool_arrivals = 0
        self.capability_token = ""
        self.authorization_check = None
        self.user_id = ""
        self.auth_session_digest = ""
        self.bound_sid = ""
        self.errors: list[str] = []
        self.approvals: list[object] = []
        self.denied_callbacks: list[str] = []
        self.stub = None

    def initialize(self) -> None:
        self.state = "ready"

    def run_authoring_turn(self, text: str) -> str:
        allowed = self.authorization_check is not None and self.authorization_check()
        if not allowed:
            raise NativeRejected("EGRESS_DENIED")
        return text

    def close(self) -> None:
        self.state = "closed"


def test_another_api_worker_runs_the_same_connection(stack, tmp_path, monkeypatch):
    monkeypatch.setattr(owner_module, "binary_matches", lambda _path: True)
    monkeypatch.setattr(owner_module, "isolation_paths_exist", lambda: True)
    monkeypatch.setattr(owner_module, "NativeSession", _IdleNative)
    bridge = NativeAuthoringBridge(
        binary=tmp_path / "peri",
        state_root=tmp_path / "state",
        repository=stack.repository,
    )
    common = {
        "registry": stack.service._registry,
        "data_registry": stack.service._data,
        "validated_inputs": SimpleNamespace(),
        "snapshots": SimpleNamespace(),
        "bridge": bridge,
        "clock": stack.clock,
    }
    opener = AuthoringService(repository=stack.repository, **common)
    other_engine = create_database_engine(stack.database_url)
    other_repository = SqlAlchemyAuthoringRepository(
        create_session_factory(other_engine),
        clock=stack.clock,
    )
    caller = AuthoringService(repository=other_repository, **common)
    replacement_engine = create_database_engine(stack.database_url)
    replacement_repository = SqlAlchemyAuthoringRepository(
        create_session_factory(replacement_engine),
        clock=stack.clock,
    )
    replacement = AuthoringService(repository=replacement_repository, **common)
    try:
        created = opener.create_session(
            stack.actor,
            project_id=stack.project_id,
            workflow_id=WORKFLOW,
        )
        opened = opener.open_connection(stack.actor, created.sid)
        assert opened.connection_state == "ready"
        assert caller._bridge is opener._bridge
        turn = caller.start_turn(
            stack.actor,
            created.sid,
            client_turn_id="from-other-worker",
            expected_revision=0,
            question={
                "intent": "suggest",
                "field_tokens": ["config.threads"],
                "preferences": [{"field_token": "config.threads", "value": 4}],
            },
        )
        assert turn.status == "completed", turn
        assert turn.proposal_digest
        assert turn.adoptable is True
        current = stack.repository.get_session(stack.actor, created.sid)
        assert current.connection_state == "ready"
        again = replacement.start_turn(
            stack.actor,
            created.sid,
            client_turn_id="from-replacement-worker",
            expected_revision=0,
            question={
                "intent": "suggest",
                "field_tokens": ["options.cores"],
                "preferences": [{"field_token": "options.cores", "value": 2}],
            },
        )
        assert again.status == "completed"
        stale = BridgeTurnRequest(
            sid=created.sid,
            turn_id="atr_" + "ab" * 16,
            generation=current.generation - 1,
            instance_id=current.instance_id or "",
            owner_epoch=current.bridge_owner_epoch or "",
            capability_token="stale-token",
            question={},
            issue_codes=(),
            schema_digest=current.schema_digest,
            projection_digest=current.projection_digest or "",
            base_revision=0,
        )
        with pytest.raises(BridgeDisconnected):
            bridge.execute_turn(stale)
        assert (
            stack.repository.get_session(stack.actor, created.sid).connection_state
            == "ready"
        )
        probe, _replayed = stack.repository.start_turn(
            stack.actor,
            created.sid,
            client_turn_id="egress-probe",
            request_digest="ab" * 32,
            question={
                "intent": "suggest",
                "field_tokens": ["config.threads"],
                "preferences": [],
            },
            expected_revision=current.revision,
            schema_digest=current.schema_digest,
            projection_digest=current.projection_digest or "",
            owner_epoch=current.bridge_owner_epoch or "",
        )
        live = BridgeTurnRequest(
            sid=created.sid,
            turn_id=probe.turn_id,
            generation=current.generation,
            instance_id=current.instance_id or "",
            owner_epoch=current.bridge_owner_epoch or "",
            capability_token="",
            question={},
            issue_codes=(),
            schema_digest=current.schema_digest,
            projection_digest=current.projection_digest or "",
            base_revision=current.revision,
        )
        assert bridge._egress_allowed(live) is True
        bridge.suppress_forwarding(
            BridgeSessionRequest(
                sid=created.sid,
                generation=current.generation + 9,
                owner_epoch=current.bridge_owner_epoch or "",
                capability_token="",
                schema_digest=current.schema_digest,
                projection_digest=current.projection_digest or "",
                instance_id=current.instance_id,
            )
        )
        bridge.suppress_forwarding(
            BridgeSessionRequest(
                sid=created.sid,
                generation=current.generation,
                owner_epoch=current.bridge_owner_epoch or "",
                capability_token="",
                schema_digest=current.schema_digest,
                projection_digest=current.projection_digest or "",
                instance_id="other-instance",
            )
        )
        assert bridge._egress_allowed(live) is True
        bridge.bind_tool_handler(opener.read_tool)
        native = bridge._sessions[created.sid]
        native.turn_context = live
        schema = bridge._dispatch("get_authoring_schema", {}, native)
        summary = bridge._dispatch("get_draft_summary", {}, native)
        assert schema["ok"] is True and summary["ok"] is True
        assert {field["token"] for field in schema["data"]["fields"]} == {
            "config.threads",
            "options.cores",
        }
        assert summary["data"]["projection_digest"] == current.projection_digest
        denied_tool = opener.read_tool(
            sid=created.sid,
            turn_id=probe.turn_id,
            generation=current.generation,
            instance_id=current.instance_id or "",
            owner_epoch=current.bridge_owner_epoch or "",
            capability_token="",
            name="get_authoring_schema",
        )
        assert denied_tool["ok"] is False
        native.turn_context = None
        assert bridge._dispatch("get_authoring_schema", {}, native)["ok"] is False
        native.turn_context = live
        original_run = native.run_authoring_turn
        native.run_authoring_turn = lambda _text: ""
        with pytest.raises(BridgeDisconnected):
            bridge.execute_turn(live)
        native.run_authoring_turn = original_run
        assert (
            stack.repository.get_session(stack.actor, created.sid).connection_state
            == "ready"
        )
        bridge.suppress_forwarding(
            BridgeSessionRequest(
                sid=created.sid,
                generation=current.generation,
                owner_epoch=current.bridge_owner_epoch or "",
                capability_token="",
                schema_digest=current.schema_digest,
                projection_digest=current.projection_digest or "",
                instance_id=current.instance_id,
            )
        )
        assert bridge._egress_allowed(live) is False
        with pytest.raises(BridgeDisconnected):
            bridge.execute_turn(live)
        native.forwarding = True
        stack.auth.revoke_session(
            stack.actor.session_digest,
            stack.clock.now,
            SessionRevocationReason.LOGOUT,
        )
        assert bridge._egress_allowed(live) is False
        with pytest.raises(AuthoringRejected) as denied:
            replacement.start_turn(
                stack.actor,
                created.sid,
                client_turn_id="after-logout",
                expected_revision=0,
                question={
                    "intent": "suggest",
                    "field_tokens": ["config.threads"],
                    "preferences": [{"field_token": "config.threads", "value": 1}],
                },
            )
        assert denied.value.reason_code == AUTHENTICATION_REQUIRED
    finally:
        bridge.shutdown()
        other_engine.dispose()
        replacement_engine.dispose()


def test_close_covers_admission_that_finishes_on_the_same_generation(stack):
    created = stack.service.create_session(
        stack.actor,
        project_id=stack.project_id,
        workflow_id=WORKFLOW,
    )
    current = stack.repository.get_session(stack.actor, created.sid)
    claim = stack.repository.open_connection(
        stack.actor,
        created.sid,
        owner_epoch=stack.bridge.epoch,
        projection_digest=current.projection_digest or "",
    )
    assert claim.record.instance_id is None
    armed = {"done": False}
    original = stack.repository.get_session

    def wrapped(actor, sid):
        record = original(actor, sid)
        if (
            sid == created.sid
            and not armed["done"]
            and record.connection_state == "blocked"
        ):
            armed["done"] = True
            stack.repository.mark_connection_ready(
                actor,
                sid,
                generation=claim.record.generation,
                instance_id="became-ready",
                owner_epoch=stack.bridge.epoch,
            )
        return record

    stack.repository.get_session = wrapped
    try:
        closed = stack.service.close_connection(stack.actor, created.sid)
    finally:
        del stack.repository.get_session
    assert closed.connection_state == "closed"
    assert stack.bridge.closes == [(claim.record.generation, "became-ready")]


def test_close_before_admission_blocks_ready_and_old_close_stays_late(stack):
    created = stack.service.create_session(
        stack.actor,
        project_id=stack.project_id,
        workflow_id=WORKFLOW,
    )
    current = stack.repository.get_session(stack.actor, created.sid)
    claim = stack.repository.open_connection(
        stack.actor,
        created.sid,
        owner_epoch=stack.bridge.epoch,
        projection_digest=current.projection_digest or "",
    )
    closed = stack.service.close_connection(stack.actor, created.sid)
    assert closed.connection_state == "closed"
    with pytest.raises(AuthoringRejected) as blocked:
        stack.repository.mark_connection_ready(
            stack.actor,
            created.sid,
            generation=claim.record.generation,
            instance_id="too-late",
            owner_epoch=stack.bridge.epoch,
        )
    assert blocked.value.reason_code == ADMISSION_FAILED
    reopened = stack.service.open_connection(stack.actor, created.sid)
    assert reopened.connection_state == "ready"
    late = stack.repository.close_connection(
        stack.actor,
        created.sid,
        expected_generation=claim.record.generation,
        expected_instance_id=None,
        expected_owner_epoch=stack.bridge.epoch,
    )
    assert late.applied is False
    assert (
        stack.repository.get_session(stack.actor, created.sid).connection_state
        == "ready"
    )


def _generating(stack):
    created = stack.service.create_session(
        stack.actor,
        project_id=stack.project_id,
        workflow_id=WORKFLOW,
    )
    current = stack.repository.get_session(stack.actor, created.sid)
    claim = stack.repository.open_connection(
        stack.actor,
        created.sid,
        owner_epoch=stack.bridge.epoch,
        projection_digest=current.projection_digest or "",
    )
    stack.repository.mark_connection_ready(
        stack.actor,
        created.sid,
        generation=claim.record.generation,
        instance_id="instance01",
        owner_epoch=stack.bridge.epoch,
    )
    turn, _replayed = stack.repository.start_turn(
        stack.actor,
        created.sid,
        client_turn_id="inflight",
        request_digest="cd" * 32,
        question={
            "intent": "suggest",
            "field_tokens": ["config.threads"],
            "preferences": [],
        },
        expected_revision=current.revision,
        schema_digest=current.schema_digest,
        projection_digest=current.projection_digest or "",
        owner_epoch=stack.bridge.epoch,
    )
    payload = {
        "method": "turn",
        "sid": created.sid,
        "generation": claim.record.generation,
        "instance_id": "instance01",
        "owner_epoch": stack.bridge.epoch,
        "turn_id": turn.turn_id,
        "schema_digest": current.schema_digest,
        "projection_digest": current.projection_digest or "",
    }
    return created.sid, turn.turn_id, payload


class _BlockingBridge:
    def __init__(self, epoch: str) -> None:
        self.epoch = epoch
        self.proceed = threading.Event()
        self.aborted = False
        self.closed: list[int] = []

    def status(self) -> BridgeStatus:
        return BridgeStatus(available=True)

    def owner_epoch(self) -> str:
        return self.epoch

    def execute_turn(self, request) -> str:
        del request
        self.proceed.wait(2)
        if self.aborted:
            raise BridgeDisconnected
        return "answer"

    def close(self, request) -> None:
        self.aborted = True
        self.closed.append(request.generation)
        self.proceed.set()


def test_caller_loss_during_turn_fences_without_the_api_worker(stack):
    sid, turn_id, payload = _generating(stack)
    bridge = _BlockingBridge(stack.bridge.epoch)
    server = OwnerServer(bridge, stack.repository, Path("/tmp/unused-owner.sock"))
    client, server_sock = socket.socketpair()
    client.close()
    server.finish_request(
        server_sock,
        (json.dumps(payload) + "\n").encode(),
        server_sock.sendall,
    )
    server_sock.close()
    record = stack.repository.get_session(stack.actor, sid)
    stored = stack.repository.peek_turn(sid, turn_id)
    assert record.connection_state == "closed"
    assert stored is not None and stored.status == "failed"
    assert bridge.closed == [payload["generation"]]
    assert stored.proposal_id is None
    _assert_late_completion_rejected(stack, sid, turn_id)


def test_completed_response_then_socket_close_keeps_the_connection(stack):
    sid, turn_id, payload = _generating(stack)
    bridge = _BlockingBridge(stack.bridge.epoch)
    bridge.proceed.set()
    server = OwnerServer(bridge, stack.repository, Path("/tmp/unused-owner-2.sock"))
    client, server_sock = socket.socketpair()
    box: list[bytes] = []
    server.finish_request(
        server_sock,
        (json.dumps(payload) + "\n").encode(),
        box.append,
    )
    client.close()
    server_sock.close()
    record = stack.repository.get_session(stack.actor, sid)
    stored = stack.repository.peek_turn(sid, turn_id)
    assert box and b'"ok": true' in box[0]
    assert record.connection_state == "ready"
    assert stored is not None and stored.status == "generating"
    assert bridge.closed == []


def test_broken_response_write_fences_before_completion(stack):
    sid, turn_id, payload = _generating(stack)
    bridge = _BlockingBridge(stack.bridge.epoch)
    bridge.proceed.set()
    order: list[str] = []
    bridge.suppress_forwarding = lambda _request: order.append("stop")
    original_close = bridge.close
    original_fail = stack.repository.fail_connection

    def close(request) -> None:
        order.append("close")
        original_close(request)

    def fail(*args, **kwargs):
        order.append("fence")
        return original_fail(*args, **kwargs)

    bridge.close = close
    stack.repository.fail_connection = fail
    server = OwnerServer(bridge, stack.repository, Path("/tmp/unused-owner-3.sock"))
    _client, server_sock = socket.socketpair()
    try:

        def fail_write(_body: bytes) -> None:
            raise BrokenPipeError

        server.finish_request(
            server_sock,
            (json.dumps(payload) + "\n").encode(),
            fail_write,
        )
    finally:
        stack.repository.fail_connection = original_fail
        server_sock.close()
    record = stack.repository.get_session(stack.actor, sid)
    stored = stack.repository.peek_turn(sid, turn_id)
    assert order == ["stop", "fence", "close"]
    assert record.connection_state == "closed"
    assert stored is not None and stored.status == "failed"
    assert stored.proposal_id is None
    _assert_late_completion_rejected(stack, sid, turn_id)


def test_control_timeout_fences_before_the_turn_finishes(stack):
    sid, turn_id, payload = _generating(stack)
    bridge = _BlockingBridge(stack.bridge.epoch)
    server = OwnerServer(
        bridge,
        stack.repository,
        Path("/tmp/unused-owner-timeout.sock"),
        control_timeout=0,
    )
    client, server_sock = socket.socketpair()
    server.finish_request(
        server_sock,
        (json.dumps(payload) + "\n").encode(),
        server_sock.sendall,
    )
    client.close()
    server_sock.close()
    record = stack.repository.get_session(stack.actor, sid)
    stored = stack.repository.peek_turn(sid, turn_id)
    assert record.connection_state == "closed"
    assert record.connection_reason_code == CONNECTION_LOST
    assert stored is not None and stored.status == "failed"
    assert stored.proposal_id is None
    assert bridge.closed == [payload["generation"]]
    _assert_late_completion_rejected(stack, sid, turn_id)


def test_completion_and_disconnect_use_commit_order(stack):
    sid, turn_id, _payload = _generating(stack)
    turn = stack.repository.peek_turn(sid, turn_id)
    assert turn is not None
    release = threading.Event()
    holding = threading.Event()

    def hold() -> None:
        holding.set()
        assert release.wait(2)

    stack.repository._before_commit = hold
    try:
        worker = threading.Thread(
            target=lambda: stack.repository.fail_connection(
                sid,
                generation=turn.generation,
                owner_epoch=turn.bridge_owner_epoch,
                reason_code=CONNECTION_LOST,
            )
        )
        worker.start()
        assert holding.wait(2)
        outcome: dict[str, object] = {}

        def complete() -> None:
            try:
                outcome["turn"] = stack.repository.complete_turn(
                    stack.actor,
                    turn_id,
                    _proposal(turn),
                    live_schema_digest=turn.schema_digest,
                    live_projection_digest=turn.projection_digest,
                )
            except AuthoringRejected as exc:
                outcome["error"] = exc.reason_code

        other = threading.Thread(target=complete)
        other.start()
        release.set()
        worker.join(2)
        other.join(2)
        stored = stack.repository.peek_turn(sid, turn_id)
        assert stored is not None
        assert stored.status == "failed"
        assert stored.proposal_id is None
        assert outcome["error"] == AUTHENTICATION_REQUIRED
    finally:
        stack.repository._before_commit = None

    sid, turn_id, _payload = _generating(stack)
    turn = stack.repository.peek_turn(sid, turn_id)
    assert turn is not None
    entered = threading.Event()
    release = threading.Event()

    def hold_completion() -> None:
        entered.set()
        assert release.wait(2)

    stack.repository._before_commit = hold_completion
    try:
        completed: dict[str, object] = {}

        def complete() -> None:
            completed["turn"] = stack.repository.complete_turn(
                stack.actor,
                turn_id,
                _proposal(turn),
                live_schema_digest=turn.schema_digest,
                live_projection_digest=turn.projection_digest,
            )

        worker = threading.Thread(target=complete)
        worker.start()
        assert entered.wait(2)
        fenced = threading.Event()

        def fence() -> None:
            stack.repository.fail_connection(
                sid,
                generation=turn.generation,
                owner_epoch=turn.bridge_owner_epoch,
                reason_code=CONNECTION_LOST,
            )
            fenced.set()

        other = threading.Thread(target=fence)
        other.start()
        release.set()
        worker.join(2)
        assert fenced.wait(2)
        other.join(2)
        stored = stack.repository.peek_turn(sid, turn_id)
        assert stored is not None
        assert stored.status == "completed"
        assert stored.decision == "invalidated"
        assert stored.proposal_id
        with pytest.raises(AuthoringRejected) as denied:
            _adopt(stack, stored)
        assert denied.value.reason_code == PROPOSAL_INVALIDATED
    finally:
        stack.repository._before_commit = None

    sid, turn_id, _payload = _generating(stack)
    turn = stack.repository.peek_turn(sid, turn_id)
    assert turn is not None
    completed_turn = stack.repository.complete_turn(
        stack.actor,
        turn_id,
        _proposal(turn),
        live_schema_digest=turn.schema_digest,
        live_projection_digest=turn.projection_digest,
    )
    receipt = _adopt(stack, completed_turn)
    stack.repository.fail_connection(
        sid,
        generation=turn.generation,
        owner_epoch=turn.bridge_owner_epoch,
        reason_code=CONNECTION_LOST,
    )
    replay = _adopt(stack, completed_turn)
    assert replay.replayed is True
    assert replay.turn.adopt_request_digest == receipt.turn.adopt_request_digest
    assert replay.turn.decision == "adopted"


def _proposal(turn):
    digest = proposal_digest(
        sid=turn.sid,
        base_revision=turn.base_revision,
        schema_digest=turn.schema_digest,
        projection_digest=turn.projection_digest,
        changes=[{"field_token": "options.cores", "value": 4}],
        reason_codes=[],
    )
    return TurnCompletion(
        kind="proposal",
        result_digest=digest,
        proposal_id=f"prop_{turn.turn_id[-32:]}",
        proposal={
            "changes": [{"field_token": "options.cores", "value": 4}],
            "reason_codes": [],
        },
        proposal_digest=digest,
    )


def _adopt(stack, turn):
    content = empty_draft_content()
    return stack.repository.adopt(
        stack.actor,
        turn.sid,
        turn.proposal_id,
        expected_revision=turn.base_revision,
        proposal_digest_value=turn.proposal_digest,
        request_digest="d" * 64,
        content=content,
        content_digest=draft_digest(content),
        changed=False,
        live_schema_digest=turn.schema_digest,
        live_projection_digest=turn.projection_digest,
    )


def _assert_late_completion_rejected(stack, sid: str, turn_id: str) -> None:
    turn = stack.repository.peek_turn(sid, turn_id)
    assert turn is not None
    with pytest.raises(AuthoringRejected) as captured:
        stack.repository.complete_turn(
            stack.actor,
            turn_id,
            _proposal(turn),
            live_schema_digest=turn.schema_digest,
            live_projection_digest=turn.projection_digest,
        )
    assert captured.value.reason_code == AUTHENTICATION_REQUIRED
