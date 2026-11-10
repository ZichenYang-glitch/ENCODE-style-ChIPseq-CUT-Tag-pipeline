"""Lowest-layer checks for the reviewed authoring contract gaps."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from types import SimpleNamespace
import json

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
from encode_pipeline.platform.authoring import (
    PROPOSAL_INVALIDATED,
    AuthoringActor,
    AuthoringRejected,
    empty_draft_content,
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
    BridgeStatus,
    BridgeTurnRequest,
)
from encode_pipeline.services.data_registry import DataRegistryService
from encode_pipeline.services.defaults import create_default_workflow_registry
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

    def status(self) -> BridgeStatus:
        return BridgeStatus(available=True)

    def owner_epoch(self) -> str:
        return self.epoch

    def admit(self, request) -> BridgeAdmission:
        del request
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
    token = stack.service._capability_tokens[(created.sid, record.generation)]
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
