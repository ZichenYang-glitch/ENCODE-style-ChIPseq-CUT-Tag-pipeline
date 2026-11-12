"""SQLite authoring drafts: CAS, expiry, and auth-session fences."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
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
from encode_pipeline.persistence.authoring import (
    SqlAlchemyAuthoringRepository,
    TurnCompletion,
)
from encode_pipeline.persistence.models import AuthSessionRow
from encode_pipeline.persistence.data_registry import SqlAlchemyDataRegistryRepository
from encode_pipeline.platform.authentication import SessionRevocationReason
from encode_pipeline.platform.authoring import (
    AUTHENTICATION_REQUIRED,
    AuthoringActor,
    AuthoringRejected,
    DRAFT_EXPIRED,
    NOT_FOUND,
    STALE_DRAFT,
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
from encode_pipeline.services.data_registry import DataRegistryService


NOW = datetime(2026, 11, 10, 12, 0, tzinfo=timezone.utc)
FAST_PROFILE = Argon2idParameters(memory_cost_kib=19456, time_cost=2, parallelism=1)


class Clock:
    def __init__(self) -> None:
        self.now = NOW

    def __call__(self) -> datetime:
        return self.now


@pytest.fixture
def clock() -> Clock:
    return Clock()


@pytest.fixture
def database_url(tmp_path) -> str:
    url = f"sqlite:///{tmp_path / 'authoring.db'}"
    upgrade_database(url)
    return url


@pytest.fixture
def auth_repository(database_url):
    engine = create_database_engine(database_url)
    repository = SqlAlchemyAuthenticationRepository(create_session_factory(engine))
    yield repository
    engine.dispose()


@pytest.fixture
def actor(auth_repository, clock) -> AuthoringActor:
    passwords = PasswordManager(FAST_PROFILE)
    admin = AccountAdministrationService(
        repository=auth_repository,
        password_manager=passwords,
        now_factory=clock,
    )
    auth = AuthenticationService(
        repository=auth_repository,
        password_manager=passwords,
        now_factory=clock,
    )
    admin.bootstrap_initial_administrator("owner", "correct horse battery staple")
    login = auth.login(
        "owner", "correct horse battery staple", client_identity="127.0.0.1"
    )
    return AuthoringActor(
        user_id=login.principal.user_id,
        session_digest=digest_session_token(login.secrets.session_token),
    )


@pytest.fixture
def project_id(database_url, clock) -> str:
    engine = create_database_engine(database_url)
    registry = DataRegistryService(
        repository=SqlAlchemyDataRegistryRepository(create_session_factory(engine)),
        now_factory=clock,
    )
    project = registry.create_project("Authoring project")
    engine.dispose()
    return project.project_id


@pytest.fixture
def repository(database_url, clock):
    engine = create_database_engine(database_url)
    repository = SqlAlchemyAuthoringRepository(
        create_session_factory(engine),
        clock=clock,
    )
    yield repository
    engine.dispose()


def test_unchanged_save_does_not_extend_expiry(repository, actor, project_id, clock):
    created = repository.create_session(
        actor,
        project_id=project_id,
        workflow_id="encode-style-chipseq-cuttag-atac-mnase",
        schema_digest="a" * 64,
        projection_digest="b" * 64,
    )
    clock.now = NOW + timedelta(hours=1)
    saved, changed = repository.save_draft(
        actor,
        created.sid,
        expected_revision=0,
        content=empty_draft_content(),
        schema_digest=created.schema_digest,
        projection_digest=created.projection_digest,
    )

    assert changed is False
    assert saved.revision == 0
    assert saved.draft_expires_at == created.draft_expires_at
    assert saved.purge_after == created.purge_after


def test_real_change_extends_expiry_once(repository, actor, project_id, clock):
    created = repository.create_session(
        actor,
        project_id=project_id,
        workflow_id="encode-style-chipseq-cuttag-atac-mnase",
        schema_digest="a" * 64,
        projection_digest="b" * 64,
    )
    clock.now = NOW + timedelta(hours=1)
    content = empty_draft_content()
    content["options"] = {"cores": 2}
    saved, changed = repository.save_draft(
        actor,
        created.sid,
        expected_revision=0,
        content=content,
        schema_digest=created.schema_digest,
        projection_digest=created.projection_digest,
    )

    assert changed is True
    assert saved.revision == 1
    assert saved.draft_expires_at == clock.now + timedelta(days=7)
    assert saved.purge_after == saved.draft_expires_at + timedelta(days=7)


def test_expiry_boundary_rejects_the_equal_instant(
    repository, actor, project_id, clock
):
    created = repository.create_session(
        actor,
        project_id=project_id,
        workflow_id="encode-style-chipseq-cuttag-atac-mnase",
        schema_digest="a" * 64,
        projection_digest="b" * 64,
    )
    _extend_login(repository, actor, created.draft_expires_at + timedelta(days=1))
    clock.now = created.draft_expires_at - timedelta(seconds=1)
    repository.save_draft(
        actor,
        created.sid,
        expected_revision=0,
        content=empty_draft_content(),
        schema_digest=created.schema_digest,
        projection_digest=created.projection_digest,
    )
    clock.now = created.draft_expires_at
    with pytest.raises(AuthoringRejected) as captured:
        repository.save_draft(
            actor,
            created.sid,
            expected_revision=0,
            content=empty_draft_content(),
            schema_digest=created.schema_digest,
            projection_digest=created.projection_digest,
        )
    assert captured.value.reason_code == DRAFT_EXPIRED


def test_two_connections_let_only_the_first_revision_win(
    database_url, actor, project_id, clock
):
    first_engine = create_database_engine(database_url)
    second_engine = create_database_engine(database_url)
    holding = threading.Event()
    release = threading.Event()
    armed = False

    def before_commit() -> None:
        if not armed:
            return
        holding.set()
        assert release.wait(5)

    writer = SqlAlchemyAuthoringRepository(
        create_session_factory(first_engine),
        clock=clock,
        before_commit=before_commit,
    )
    reader = SqlAlchemyAuthoringRepository(
        create_session_factory(second_engine),
        clock=clock,
    )
    created = writer.create_session(
        actor,
        project_id=project_id,
        workflow_id="encode-style-chipseq-cuttag-atac-mnase",
        schema_digest="a" * 64,
        projection_digest="b" * 64,
    )
    armed = True
    first = empty_draft_content()
    first["options"] = {"cores": 2}
    second = empty_draft_content()
    second["options"] = {"cores": 3}
    outcome: list[str] = []

    def save_second() -> None:
        try:
            reader.save_draft(
                actor,
                created.sid,
                expected_revision=0,
                content=second,
                schema_digest=created.schema_digest,
                projection_digest=created.projection_digest,
            )
        except AuthoringRejected as exc:
            outcome.append(exc.reason_code)
        else:
            outcome.append("saved")

    def save_first() -> None:
        writer.save_draft(
            actor,
            created.sid,
            expected_revision=0,
            content=first,
            schema_digest=created.schema_digest,
            projection_digest=created.projection_digest,
        )

    first_thread = threading.Thread(target=save_first)
    second_thread = threading.Thread(target=save_second)
    first_thread.start()
    assert holding.wait(5)
    second_thread.start()
    release.set()
    first_thread.join(5)
    second_thread.join(5)

    stored = reader.get_session(actor, created.sid)
    first_engine.dispose()
    second_engine.dispose()
    assert outcome == [STALE_DRAFT]
    assert stored.revision == 1
    assert stored.draft_digest == draft_digest(first)


def test_logout_before_completion_discards_the_answer(
    repository, auth_repository, actor, project_id, clock
):
    turn = _ready_turn(repository, actor, project_id)
    assert auth_repository.revoke_session(
        actor.session_digest,
        clock.now,
        SessionRevocationReason.LOGOUT,
    )
    with pytest.raises(AuthoringRejected) as captured:
        repository.complete_turn(
            actor,
            turn.turn_id,
            _proposal(turn),
            live_schema_digest=turn.schema_digest,
            live_projection_digest=turn.projection_digest,
        )
    stored = repository.peek_turn(turn.sid, turn.turn_id)

    assert captured.value.reason_code == AUTHENTICATION_REQUIRED
    assert stored is not None
    assert stored.status == "failed"
    assert stored.decision is None
    assert stored.proposal_id is None


def test_completion_before_logout_stays_historical_but_not_adoptable(
    repository, auth_repository, actor, project_id, clock
):
    turn = _ready_turn(repository, actor, project_id)
    completed = repository.complete_turn(
        actor,
        turn.turn_id,
        _proposal(turn),
        live_schema_digest=turn.schema_digest,
        live_projection_digest=turn.projection_digest,
    )
    assert completed.decision == "pending"
    assert auth_repository.revoke_session(
        actor.session_digest,
        clock.now,
        SessionRevocationReason.LOGOUT,
    )
    stored = repository.peek_turn(turn.sid, turn.turn_id)

    assert stored is not None
    assert stored.status == "completed"
    assert stored.decision == "invalidated"
    assert stored.invalidation_reason == "AUTH_SESSION_ENDED"
    with pytest.raises(AuthoringRejected) as captured:
        repository.adopt(
            actor,
            turn.sid,
            completed.proposal_id,
            expected_revision=0,
            proposal_digest_value=completed.proposal_digest,
            request_digest="c" * 64,
            content=empty_draft_content(),
            content_digest=draft_digest(empty_draft_content()),
            changed=False,
            live_schema_digest=turn.schema_digest,
            live_projection_digest=turn.projection_digest,
        )
    assert captured.value.reason_code == AUTHENTICATION_REQUIRED


def test_repeated_adopt_returns_the_original_receipt(
    repository, actor, project_id, clock
):
    turn = _ready_turn(repository, actor, project_id)
    completed = repository.complete_turn(
        actor,
        turn.turn_id,
        _proposal(turn),
        live_schema_digest=turn.schema_digest,
        live_projection_digest=turn.projection_digest,
    )
    content = empty_draft_content()
    content["options"] = {"cores": 4}
    first = repository.adopt(
        actor,
        turn.sid,
        completed.proposal_id,
        expected_revision=0,
        proposal_digest_value=completed.proposal_digest,
        request_digest="d" * 64,
        content=content,
        content_digest=draft_digest(content),
        changed=True,
        live_schema_digest=turn.schema_digest,
        live_projection_digest=turn.projection_digest,
    )
    clock.now = NOW + timedelta(minutes=45)
    second = repository.adopt(
        actor,
        turn.sid,
        completed.proposal_id,
        expected_revision=0,
        proposal_digest_value=completed.proposal_digest,
        request_digest="d" * 64,
        content=None,
        content_digest=None,
        changed=None,
        live_schema_digest="e" * 64,
        live_projection_digest="f" * 64,
    )

    assert first.replayed is False
    assert second.replayed is True
    assert second.turn.adopted_at == first.turn.adopted_at
    assert second.turn.adoption_draft_expires_at == first.turn.adoption_draft_expires_at
    assert second.record.revision == 1


def test_purge_removes_authoring_rows_and_leaves_snapshot_ids_untouched(
    repository, actor, project_id, database_url, clock
):
    created = repository.create_session(
        actor,
        project_id=project_id,
        workflow_id="encode-style-chipseq-cuttag-atac-mnase",
        schema_digest="a" * 64,
        projection_digest="b" * 64,
    )
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        connection.exec_driver_sql(
            """
            INSERT INTO validated_input_snapshots (
                snapshot_id, workflow_id, adapter_version, schema_version,
                schema_dialect, canonical_payload, payload_digest_scheme,
                payload_digest, validation_outcome, validation_issue_codes,
                validated_at, expires_at, build_adapter_version, build_scheme,
                build_logical_entrypoint, build_digest, build_captured_at
            ) VALUES (
                'vsnap_untouched', 'workflow', '1', '1', 'dialect', '{}',
                'scheme', 'aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa',
                'adapter_validation_succeeded', '[]',
                '2026-11-10 12:00:00', '2026-11-17 12:00:00', '1', 'scheme',
                'entrypoint', 'bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb',
                '2026-11-10 12:00:00'
            )
            """
        )
    _extend_login(repository, actor, created.purge_after + timedelta(days=1))
    clock.now = created.purge_after
    with pytest.raises(AuthoringRejected) as captured:
        repository.get_session(actor, created.sid)
    assert captured.value.reason_code == "NOT_FOUND"
    with engine.connect() as connection:
        sessions = connection.exec_driver_sql(
            "SELECT count(*) FROM agent_authoring_sessions"
        ).scalar()
        turns = connection.exec_driver_sql(
            "SELECT count(*) FROM agent_authoring_turns"
        ).scalar()
        snapshots = connection.exec_driver_sql(
            "SELECT count(*) FROM validated_input_snapshots"
        ).scalar()
    engine.dispose()
    assert sessions == 0
    assert turns == 0
    assert snapshots == 1


def _extend_login(repository, actor, until: datetime) -> None:
    session_factory = repository._session_factory
    with session_factory.begin() as session:
        row = session.get(AuthSessionRow, actor.session_digest)
        assert row is not None
        row.expires_at = until


def _ready_turn(repository, actor, project_id):
    created = repository.create_session(
        actor,
        project_id=project_id,
        workflow_id="encode-style-chipseq-cuttag-atac-mnase",
        schema_digest="a" * 64,
        projection_digest="b" * 64,
    )
    claim = repository.open_connection(
        actor,
        created.sid,
        owner_epoch="epoch-1",
        projection_digest="b" * 64,
    )
    repository.mark_connection_ready(
        actor,
        created.sid,
        generation=claim.record.generation,
        instance_id="instance-1",
        owner_epoch="epoch-1",
    )
    turn, _replayed = repository.start_turn(
        actor,
        created.sid,
        client_turn_id="turn-1",
        request_digest="e" * 64,
        question={"intent": "suggest", "field_tokens": [], "preferences": []},
        expected_revision=0,
        schema_digest="a" * 64,
        projection_digest="b" * 64,
        owner_epoch="epoch-1",
    )
    return turn


def _proposal(turn) -> TurnCompletion:
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
        proposal_id="prop_" + "ab" * 16,
        proposal={
            "changes": [{"field_token": "options.cores", "value": 4}],
            "reason_codes": [],
        },
        proposal_digest=digest,
    )


def test_receipt_discovery_is_bounded_and_survives_login_and_expiry(
    repository, actor, project_id, clock, auth_repository
):
    first = _ready_turn(repository, actor, project_id)
    _adopt(repository, actor, first, cores=2)
    current = repository.get_session(actor, first.sid)
    second, _replayed = repository.start_turn(
        actor,
        first.sid,
        client_turn_id="turn-2",
        request_digest="f" * 64,
        question={"intent": "suggest", "field_tokens": [], "preferences": []},
        expected_revision=current.revision,
        schema_digest=first.schema_digest,
        projection_digest=first.projection_digest,
        owner_epoch="epoch-1",
    )
    _adopt(repository, actor, second, cores=4)
    pending, _replayed = repository.start_turn(
        actor,
        first.sid,
        client_turn_id="turn-3",
        request_digest="a" * 64,
        question={"intent": "suggest", "field_tokens": [], "preferences": []},
        expected_revision=repository.get_session(actor, first.sid).revision,
        schema_digest=first.schema_digest,
        projection_digest=first.projection_digest,
        owner_epoch="epoch-1",
    )
    repository.complete_turn(
        actor,
        pending.turn_id,
        _proposal(pending),
        live_schema_digest=pending.schema_digest,
        live_projection_digest=pending.projection_digest,
    )
    before = repository.get_session(actor, first.sid)
    stored = {
        turn_id: _receipt_snapshot(repository.peek_turn(first.sid, turn_id))
        for turn_id in (first.turn_id, second.turn_id, pending.turn_id)
    }
    page = repository.list_adoption_receipts(
        actor,
        first.sid,
        cursor_time=None,
        cursor_turn=None,
        limit=1,
    )
    assert [item.turn_id for item in page.turns] == [max(first.turn_id, second.turn_id)]
    assert page.next_cursor
    rest = repository.list_adoption_receipts(
        actor,
        first.sid,
        cursor_time=_cursor_time(page.next_cursor),
        cursor_turn=_cursor_turn(page.next_cursor),
        limit=1,
    )
    assert [item.turn_id for item in rest.turns] == [min(first.turn_id, second.turn_id)]
    assert rest.next_cursor is None
    assert pending.turn_id not in {item.turn_id for item in page.turns + rest.turns}
    after = repository.get_session(actor, first.sid)
    assert after.revision == before.revision
    assert after.draft_expires_at == before.draft_expires_at
    assert after.purge_after == before.purge_after
    for turn_id, snapshot in stored.items():
        assert _receipt_snapshot(repository.peek_turn(first.sid, turn_id)) == snapshot

    repository.close_connection(
        actor,
        first.sid,
        expected_generation=after.generation,
        expected_instance_id=after.instance_id,
        expected_owner_epoch=after.bridge_owner_epoch,
    )
    auth_repository.revoke_session(
        actor.session_digest,
        clock.now,
        SessionRevocationReason.LOGOUT,
    )
    passwords = PasswordManager(FAST_PROFILE)
    renewed_login = AuthenticationService(
        repository=auth_repository,
        password_manager=passwords,
        now_factory=clock,
    ).login("owner", "correct horse battery staple", client_identity="127.0.0.1")
    renewed = AuthoringActor(
        user_id=actor.user_id,
        session_digest=digest_session_token(renewed_login.secrets.session_token),
    )
    visible = repository.list_adoption_receipts(
        renewed,
        first.sid,
        cursor_time=None,
        cursor_turn=None,
        limit=20,
    )
    assert {item.turn_id for item in visible.turns} == {first.turn_id, second.turn_id}
    assert visible.record.revision == before.revision
    _extend_login(auth_repository, renewed, before.purge_after + timedelta(days=1))

    clock.now = before.draft_expires_at
    expired = repository.list_adoption_receipts(
        renewed,
        first.sid,
        cursor_time=None,
        cursor_turn=None,
        limit=20,
    )
    assert len(expired.turns) == 2
    clock.now = before.purge_after - timedelta(seconds=1)
    assert (
        len(
            repository.list_adoption_receipts(
                renewed,
                first.sid,
                cursor_time=None,
                cursor_turn=None,
                limit=20,
            ).turns
        )
        == 2
    )
    clock.now = before.purge_after
    with pytest.raises(AuthoringRejected) as missing:
        repository.list_adoption_receipts(
            renewed,
            first.sid,
            cursor_time=None,
            cursor_turn=None,
            limit=20,
        )
    assert missing.value.reason_code == NOT_FOUND


def _adopt(repository, actor, turn, *, cores: int):
    completed = repository.complete_turn(
        actor,
        turn.turn_id,
        _numbered_proposal(turn),
        live_schema_digest=turn.schema_digest,
        live_projection_digest=turn.projection_digest,
    )
    content = empty_draft_content()
    content["options"] = {"cores": cores}
    return repository.adopt(
        actor,
        turn.sid,
        completed.proposal_id,
        expected_revision=turn.base_revision,
        proposal_digest_value=completed.proposal_digest,
        request_digest="d" * 64,
        content=content,
        content_digest=draft_digest(content),
        changed=True,
        live_schema_digest=turn.schema_digest,
        live_projection_digest=turn.projection_digest,
    )


def _numbered_proposal(turn) -> TurnCompletion:
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
        proposal_id="prop_" + turn.turn_id[-32:],
        proposal={
            "changes": [{"field_token": "options.cores", "value": 4}],
            "reason_codes": [],
        },
        proposal_digest=digest,
    )


def _receipt_snapshot(turn):
    assert turn is not None
    return (
        turn.status,
        turn.decision,
        turn.proposal_id,
        turn.adopted_at,
        turn.result_revision,
        turn.adoption_changed,
        turn.adoption_draft_expires_at,
        turn.base_revision,
    )


def _cursor_time(value: str):
    from encode_pipeline.persistence.authoring import decode_receipt_cursor

    return decode_receipt_cursor(value)[0]


def _cursor_turn(value: str):
    from encode_pipeline.persistence.authoring import decode_receipt_cursor

    return decode_receipt_cursor(value)[1]
