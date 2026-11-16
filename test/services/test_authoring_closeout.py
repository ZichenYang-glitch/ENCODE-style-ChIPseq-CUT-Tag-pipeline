"""Service paths the full suite left uncovered: expiry, admission, and replay."""

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
from encode_pipeline.persistence.authoring import (
    AdoptionReceiptPage,
    SqlAlchemyAuthoringRepository,
)
from encode_pipeline.persistence.models import AuthSessionRow
from encode_pipeline.platform.authoring import (
    ADMISSION_FAILED,
    BRIDGE_UNAVAILABLE,
    CAPABILITY_CHANGED,
    CAPABILITY_UNAVAILABLE,
    CONNECTION_LOST,
    DRAFT_EXPIRED,
    INVALID_REQUEST,
    MODEL_OUTPUT_REJECTED,
    NOT_FOUND,
    PROPOSAL_EXPIRED,
    REQUEST_TOO_LARGE,
    WORKFLOW_UNAVAILABLE,
    AuthoringActor,
    AuthoringRejected,
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
)
from encode_pipeline.services.data_registry import DataRegistryService
from encode_pipeline.services.defaults import create_default_workflow_registry

NOW = datetime(2026, 11, 10, 12, 0, tzinfo=timezone.utc)
WORKFLOW = "encode-style-chipseq-cuttag-atac-mnase"
QUESTION = {
    "intent": "suggest",
    "field_tokens": ["config.threads"],
    "preferences": [{"field_token": "config.threads", "value": 4}],
}


class Clock:
    def __init__(self) -> None:
        self.now = NOW

    def __call__(self) -> datetime:
        return self.now


class ScriptedBridge:
    def __init__(self) -> None:
        self.epoch = "owner-epoch-0001"
        self.turns: list[str] = []
        self.answer = "proposal"
        self.admit_error: BaseException | None = None
        self.admit_state = "ready"
        self.available = True

    def status(self) -> BridgeStatus:
        return BridgeStatus(available=self.available)

    def owner_epoch(self) -> str:
        return self.epoch

    def admit(self, request):
        if self.admit_error is not None:
            raise self.admit_error
        return BridgeAdmission(instance_id="instance01", state=self.admit_state)

    def execute_turn(self, request) -> str:
        self.turns.append(request.turn_id)
        if self.answer == "proposal":
            body = {
                "base_revision": request.base_revision,
                "changes": [{"field_token": "config.threads", "value": 4}],
                "kind": "proposal",
                "reason_codes": [],
                "schema_digest": request.schema_digest,
            }
        elif self.answer == "disconnect":
            raise BridgeDisconnected
        else:
            return self.answer
        return json.dumps(body, separators=(",", ":"), sort_keys=True)

    def close(self, request) -> None:
        return None


@pytest.fixture
def stack(tmp_path):
    clock = Clock()
    database_url = f"sqlite:///{tmp_path / 'authoring.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
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
    data = DataRegistryService(
        repository=__import__(
            "encode_pipeline.persistence.data_registry",
            fromlist=["SqlAlchemyDataRegistryRepository"],
        ).SqlAlchemyDataRegistryRepository(create_session_factory(engine)),
        now_factory=clock,
    )
    project = data.create_project("Closeout")
    repository = SqlAlchemyAuthoringRepository(
        create_session_factory(engine), clock=clock
    )
    bridge = ScriptedBridge()
    service = AuthoringService(
        repository=repository,
        registry=create_default_workflow_registry(),
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
        data=data,
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


def _ready(stack) -> str:
    created = stack.service.create_session(
        stack.actor, project_id=stack.project_id, workflow_id=WORKFLOW
    )
    opened = stack.service.open_connection(stack.actor, created.sid)
    assert opened.connection_state == "ready"
    return created.sid


def _turn(stack, sid: str, client_turn_id: str):
    record = stack.repository.get_session(stack.actor, sid)
    return stack.service.start_turn(
        stack.actor,
        sid,
        client_turn_id=client_turn_id,
        expected_revision=record.revision,
        question=QUESTION,
    )


def test_project_and_session_pages_honor_cursors(stack) -> None:
    second = stack.data.create_project("Second")
    first_page = stack.service.list_projects(stack.actor, cursor=None, limit=1)
    assert len(first_page.projects) == 1
    assert first_page.next_cursor
    rest = stack.service.list_projects(
        stack.actor, cursor=first_page.next_cursor, limit=10
    )
    assert {item[0] for item in first_page.projects + rest.projects} >= {
        stack.project_id,
        second.project_id,
    }
    with pytest.raises(AuthoringRejected) as invalid:
        stack.service.list_projects(stack.actor, cursor="not-a-cursor", limit=10)
    assert invalid.value.reason_code == INVALID_REQUEST
    stack.data.archive_project(second.project_id)
    visible = stack.service.list_projects(stack.actor, cursor=None, limit=10)
    assert second.project_id not in {item[0] for item in visible.projects}
    with pytest.raises(AuthoringRejected) as missing:
        stack.service.create_session(
            stack.actor, project_id=second.project_id, workflow_id=WORKFLOW
        )
    assert missing.value.reason_code == NOT_FOUND
    with pytest.raises(AuthoringRejected) as unknown:
        stack.service.create_session(
            stack.actor, project_id=stack.project_id, workflow_id="missing-workflow"
        )
    assert unknown.value.reason_code == WORKFLOW_UNAVAILABLE

    created = stack.service.create_session(
        stack.actor, project_id=stack.project_id, workflow_id=WORKFLOW
    )
    other = stack.service.create_session(
        stack.actor, project_id=stack.project_id, workflow_id=WORKFLOW
    )
    page, cursor = stack.service.list_sessions(
        stack.actor,
        project_id=stack.project_id,
        workflow_id=WORKFLOW,
        cursor=None,
        limit=1,
    )
    assert len(page) == 1
    assert cursor
    rest_sessions, next_cursor = stack.service.list_sessions(
        stack.actor,
        project_id=stack.project_id,
        workflow_id=WORKFLOW,
        cursor=cursor,
        limit=10,
    )
    assert {item.sid for item in page + rest_sessions} == {created.sid, other.sid}
    assert next_cursor is None
    with pytest.raises(AuthoringRejected) as bad_session_cursor:
        stack.service.list_sessions(
            stack.actor,
            project_id=None,
            workflow_id=None,
            cursor="%%%",
            limit=10,
        )
    assert bad_session_cursor.value.reason_code == INVALID_REQUEST


def test_sample_and_input_pages_skip_empty_and_archived_rows(stack) -> None:
    earlier = NOW
    later = NOW + timedelta(seconds=1)
    samples = (
        SimpleNamespace(
            sample=SimpleNamespace(sample_id="smp_a", created_at=earlier),
            display_name="A",
            revision=SimpleNamespace(sample_revision_id="smpr_a"),
        ),
        SimpleNamespace(
            sample=SimpleNamespace(sample_id="smp_b", created_at=later),
            display_name="B",
            revision=SimpleNamespace(sample_revision_id="smpr_b"),
        ),
    )
    stack.service._data.list_samples = lambda project_id: samples
    page = stack.service.list_samples(
        stack.actor, stack.project_id, cursor=None, limit=1
    )
    assert page.samples[0][0] == "smp_a"
    assert page.next_cursor
    rest = stack.service.list_samples(
        stack.actor, stack.project_id, cursor=page.next_cursor, limit=10
    )
    assert rest.samples[0][0] == "smp_b"

    class Registry:
        def get_active_project_storage_pool_binding(self, project_id):
            return object()

        def list_input_files(self, project_id):
            return [
                SimpleNamespace(
                    archived_at=None,
                    input_file_id="inpf_b",
                    stable_key="b",
                    created_at=later,
                ),
                SimpleNamespace(
                    archived_at=NOW,
                    input_file_id="inpf_old",
                    stable_key="old",
                    created_at=earlier,
                ),
                SimpleNamespace(
                    archived_at=None,
                    input_file_id="inpf_empty",
                    stable_key="empty",
                    created_at=earlier,
                ),
                SimpleNamespace(
                    archived_at=None,
                    input_file_id="inpf_a",
                    stable_key="a",
                    created_at=earlier,
                ),
            ]

        def list_input_file_revisions(self, file_id):
            if file_id == "inpf_empty":
                return []
            return [
                SimpleNamespace(
                    input_file_revision_id=f"inpfr_{file_id}",
                    revision_number=1,
                    size_bytes=4,
                    content_sha256="ab" * 32,
                )
            ]

    stack.service._input_registry = Registry()
    listed = stack.service.list_input_revisions(
        stack.actor, stack.project_id, cursor=None, limit=1
    )
    assert listed.revisions[0][0] == "inpf_a"
    followed = stack.service.list_input_revisions(
        stack.actor, stack.project_id, cursor=listed.next_cursor, limit=10
    )
    assert [item[0] for item in followed.revisions] == ["inpf_b"]

    class Missing:
        def get_active_project_storage_pool_binding(self, project_id):
            raise KeyError(project_id)

        def list_input_files(self, project_id):
            raise AssertionError("missing binding has no files")

    stack.service._input_registry = Missing()
    empty = stack.service.list_input_revisions(
        stack.actor, stack.project_id, cursor=None, limit=10
    )
    assert empty.revisions == ()

    class Archived:
        def get_active_project_storage_pool_binding(self, project_id):
            raise ValueError("StoragePool is archived")

        def list_input_files(self, project_id):
            raise AssertionError("archived pool is not listed")

    stack.service._input_registry = Archived()
    archived = stack.service.list_input_revisions(
        stack.actor, stack.project_id, cursor=None, limit=10
    )
    assert archived.revisions == ()

    class Broken:
        def get_active_project_storage_pool_binding(self, project_id):
            raise ValueError("missing pool")

    stack.service._input_registry = Broken()
    with pytest.raises(ValueError, match="missing pool"):
        stack.service.list_input_revisions(
            stack.actor, stack.project_id, cursor=None, limit=10
        )
    stack.service._input_registry = None
    with pytest.raises(RuntimeError, match="unavailable"):
        stack.service.list_input_revisions(
            stack.actor, stack.project_id, cursor=None, limit=10
        )


def test_replay_does_not_execute_a_second_native_turn(stack) -> None:
    sid = _ready(stack)
    first = _turn(stack, sid, "same-turn")
    second = _turn(stack, sid, "same-turn")
    assert stack.bridge.turns == [first.turn_id]
    assert second.turn_id == first.turn_id
    assert second.proposal_digest == first.proposal_digest
    assert second.proposal_id == first.proposal_id


def test_lost_response_and_rejected_model_output_do_not_stay_generating(stack) -> None:
    sid = _ready(stack)
    stack.bridge.answer = "disconnect"
    lost = _turn(stack, sid, "lost")
    assert lost.status == "failed"
    assert lost.failure_reason == CONNECTION_LOST
    reopened = stack.service.open_connection(stack.actor, sid)
    assert reopened.connection_state == "ready"
    stack.bridge.answer = "{"
    rejected = _turn(stack, sid, "bad-model")
    assert rejected.status == "failed"
    assert rejected.failure_reason == MODEL_OUTPUT_REJECTED


def test_question_shape_and_size_are_rejected_before_a_turn(stack, monkeypatch) -> None:
    sid = _ready(stack)
    record = stack.repository.get_session(stack.actor, sid)
    cases = [
        {"intent": "suggest"},
        {"intent": "other", "field_tokens": [], "preferences": []},
        {
            "intent": "suggest",
            "field_tokens": ["config.threads", "config.threads"],
            "preferences": [],
        },
        {
            "intent": "suggest",
            "field_tokens": ["missing.field"],
            "preferences": [],
        },
        {
            "intent": "suggest",
            "field_tokens": ["config.threads"],
            "preferences": [{"field_token": "other", "value": 1}],
        },
        {
            "intent": "suggest",
            "field_tokens": ["config.threads"],
            "preferences": [{"field_token": "config.threads", "value": "nope"}],
        },
    ]
    for question in cases:
        with pytest.raises(AuthoringRejected) as captured:
            stack.service.start_turn(
                stack.actor,
                sid,
                client_turn_id=f"bad-{len(question)}",
                expected_revision=record.revision,
                question=question,
            )
        assert captured.value.reason_code == INVALID_REQUEST
    monkeypatch.setattr("encode_pipeline.services.authoring.MAX_QUESTION_BYTES", 20)
    with pytest.raises(AuthoringRejected) as oversized:
        _turn(stack, sid, "oversized")
    assert oversized.value.reason_code == REQUEST_TOO_LARGE


def test_admission_failures_abandon_the_claim(stack) -> None:
    created = stack.service.create_session(
        stack.actor, project_id=stack.project_id, workflow_id=WORKFLOW
    )
    stack.bridge.admit_state = "starting"
    with pytest.raises(AuthoringRejected) as starting:
        stack.service.open_connection(stack.actor, created.sid)
    assert starting.value.reason_code == ADMISSION_FAILED
    stack.bridge.admit_state = "ready"
    stack.bridge.admit_error = BridgeDisconnected()
    with pytest.raises(AuthoringRejected) as disconnected:
        stack.service.open_connection(stack.actor, created.sid)
    assert disconnected.value.reason_code == ADMISSION_FAILED
    stack.bridge.admit_error = RuntimeError("boom")
    with pytest.raises(AuthoringRejected) as exploded:
        stack.service.open_connection(stack.actor, created.sid)
    assert exploded.value.reason_code == ADMISSION_FAILED
    stack.bridge.admit_error = None
    stack.bridge.available = False
    with pytest.raises(AuthoringRejected) as unavailable:
        stack.service.open_connection(stack.actor, created.sid)
    assert unavailable.value.reason_code == BRIDGE_UNAVAILABLE
    stack.bridge.available = True
    original = stack.service._identity
    stack.service._identity = lambda workflow_id: ("a" * 64, None, None)
    with pytest.raises(AuthoringRejected) as capability:
        stack.service.open_connection(stack.actor, created.sid)
    assert capability.value.reason_code == CAPABILITY_UNAVAILABLE
    stack.service._identity = original
    opened = stack.service.open_connection(stack.actor, created.sid)
    assert opened.connection_state == "ready"
    closed = stack.service.close_connection(stack.actor, created.sid)
    again = stack.service.close_connection(stack.actor, created.sid)
    assert closed.connection_state == "closed"
    assert again.connection_state == "closed"


def test_proposal_and_draft_ttl_block_adoption_and_new_turns(stack) -> None:
    sid = _ready(stack)
    proposal = _turn(stack, sid, "expiring")
    assert proposal.adoptable is True
    stack.clock.now = NOW + timedelta(minutes=31)
    _extend_login(stack.auth, stack.actor, stack.clock.now + timedelta(days=1))
    visible = stack.service.get_turn(stack.actor, sid, proposal.turn_id)
    assert visible.adoptable is False
    assert visible.reason_code == PROPOSAL_EXPIRED
    with pytest.raises(AuthoringRejected) as expired_proposal:
        stack.service.reject(
            stack.actor,
            sid,
            proposal.proposal_id,
            proposal_digest_value=proposal.proposal_digest,
        )
    assert expired_proposal.value.reason_code == PROPOSAL_EXPIRED

    fresh = _ready(stack)
    live = _turn(stack, fresh, "reject-me")
    rejected = stack.service.reject(
        stack.actor,
        fresh,
        live.proposal_id,
        proposal_digest_value=live.proposal_digest,
    )
    assert rejected.proposal_state == "rejected"
    assert rejected.adoptable is False
    again = stack.service.reject(
        stack.actor,
        fresh,
        live.proposal_id,
        proposal_digest_value=live.proposal_digest,
    )
    assert again.proposal_state == "rejected"
    pending = _ready(stack)
    pending_turn = _turn(stack, pending, "still-pending")
    stack.clock.now = NOW + timedelta(days=8)
    _extend_login(stack.auth, stack.actor, stack.clock.now + timedelta(days=1))
    hidden = stack.service.get_turn(stack.actor, pending, pending_turn.turn_id)
    assert hidden.adoptable is False
    assert hidden.reason_code == DRAFT_EXPIRED
    with pytest.raises(AuthoringRejected) as expired_draft:
        _turn(stack, fresh, "after-expiry")
    assert expired_draft.value.reason_code == DRAFT_EXPIRED
    receipts = stack.service.list_adoption_receipts(
        stack.actor, fresh, cursor=None, limit=10
    )
    assert receipts.draft_state == "expired"


def test_pending_receipt_rows_are_not_presented_as_adoptions(stack) -> None:
    sid = _ready(stack)
    proposal = _turn(stack, sid, "pending-receipt")
    record = stack.repository.get_session(stack.actor, sid)
    turn = stack.repository.peek_turn(sid, proposal.turn_id)

    def pending_page(*_args, **_kwargs):
        return AdoptionReceiptPage(record, (turn,), None)

    stack.repository.list_adoption_receipts = pending_page
    page = stack.service.list_adoption_receipts(stack.actor, sid, cursor=None, limit=10)
    assert page.entries == ()


def test_schema_change_while_connected_drops_the_old_capability(stack) -> None:
    sid = _ready(stack)
    original = stack.service._identity

    def changed(workflow_id: str):
        _schema, _projection, fields = original(workflow_id)
        return "b" * 64, "c" * 64, fields

    stack.service._identity = changed
    opened = stack.service.open_connection(stack.actor, sid)
    assert opened.schema_digest == "b" * 64
    stack.service._identity = original

    def fail_start(*_args, **_kwargs):
        raise AuthoringRejected(CAPABILITY_CHANGED)

    stack.repository.start_turn = fail_start
    with pytest.raises(AuthoringRejected) as changed_capability:
        _turn(stack, sid, "changed")
    assert changed_capability.value.reason_code == CAPABILITY_CHANGED


def test_saving_a_changed_schema_closes_the_live_connection(stack) -> None:
    sid = _ready(stack)
    record = stack.repository.get_session(stack.actor, sid)
    original = stack.service._identity

    def changed(workflow_id: str):
        _schema, _projection, fields = original(workflow_id)
        return "d" * 64, "e" * 64, fields

    stack.service._identity = changed
    saved = stack.service.save_draft(
        stack.actor,
        sid,
        expected_revision=record.revision,
        content=record.content,
    )
    assert saved.schema_digest == "d" * 64
    assert saved.connection_state != "ready"


def test_malformed_model_output_fails_the_turn(stack) -> None:
    sid = _ready(stack)
    record = stack.repository.get_session(stack.actor, sid)
    schema = record.schema_digest
    base = record.revision

    def body(**overrides: object) -> str:
        value: dict[str, object] = {
            "base_revision": base,
            "changes": [],
            "kind": "explanation",
            "reason_codes": [],
            "schema_digest": schema,
        }
        value.update(overrides)
        return json.dumps(value)

    answers = [
        body(unexpected=True),
        body(schema_digest="f" * 64),
        body(base_revision=base + 1),
        body(kind="other"),
        body(changes="no"),
        body(reason_codes=["NOT_AN_ISSUE"]),
        body(kind="proposal", changes=[{"field_token": "config.threads"}]),
        body(
            kind="explanation",
            changes=[{"field_token": "config.threads", "value": 4}],
        ),
    ]
    for index, answer in enumerate(answers):
        stack.bridge.answer = answer
        failed = _turn(stack, sid, f"model-{index}")
        assert failed.status == "failed"
        assert failed.failure_reason == MODEL_OUTPUT_REJECTED
