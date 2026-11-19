"""Peri run decisions through the real API and a disposable SQLite database."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timedelta, timezone
import json

import pytest
from sqlalchemy import func, select
from sqlalchemy.orm import Session

from api_test_client import ApiTestClient
from encode_pipeline.persistence.models import (
    RunExecutionAssignmentRow,
    RunRow,
    SecurityAuditEventRow,
)
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.authentication import UserRole
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.platform.security_audit import build_audit_resource
from encode_pipeline.services.authentication import (
    new_session_record,
    new_session_secrets,
)
from encode_pipeline.services.authoring_bridge import BridgeAdmission, BridgeStatus
from encode_pipeline.services.data_registry import DataRegistryService, SampleImportSpec
from encode_pipeline.services.input_registry import (
    InputRegistryAdminService,
    InputRegistryService,
)
from encode_pipeline.services import peri_run_confirmation
from encode_pipeline.services.private_storage_pools import (
    PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
    load_private_storage_pool_config,
)
from encode_pipeline.services.run_cancellation import RunCancellationService


WORKFLOW = "encode-style-chipseq-cuttag-atac-mnase"
SESSIONS = "/api/v1/agent/authoring-sessions"


class SyntheticPeri:
    def status(self):
        return BridgeStatus(available=True)

    def owner_epoch(self):
        return "test-peri-owner"

    def admit(self, request):
        return BridgeAdmission(instance_id="test-instance", state="ready")

    def close(self, request):
        pass

    def execute_turn(self, request):
        return '{"kind":"run_action","action":"create","confirmed":true}'


class StopRecorder:
    backend = "rq"
    queue_name = "runs"

    def __init__(self):
        self.run_ids = []

    def request_stop(self, assignment):
        self.run_ids.append(assignment.run_id)


@pytest.fixture
def peri_case(reference_ready_app, tmp_path, monkeypatch):
    app = reference_ready_app
    persistence = app.state.persistence
    app.state.authoring_service._bridge = SyntheticPeri()
    stop = StopRecorder()
    app.state.run_cancellation_service = RunCancellationService(
        app.state.run_service, stop
    )

    def forbid_enqueue(*args, **kwargs):
        pytest.fail("Peri run actions must never enqueue execution")

    monkeypatch.setattr(app.state.run_queue, "enqueue_execution", forbid_enqueue)
    registry = DataRegistryService(repository=persistence.data_registry_repository)
    project = registry.create_project("Peri run action fixture")
    sample = registry.import_samples(
        project.project_id,
        rows=(
            SampleImportSpec(
                stable_key="sample-a", display_name="Sample A", attributes={}
            ),
        ),
    )[0].revision
    storage = tmp_path / "read-pool"
    storage.mkdir()
    (storage / "small.fastq").write_bytes(b"@read\nACGT\n+\n!!!!\n")
    config = tmp_path / "storage-pools.json"
    config.write_text(
        json.dumps(
            {
                "schema_version": PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
                "storage_pool_roots": {"fixture": str(storage)},
            }
        ),
        encoding="utf-8",
    )
    app.state.authoring_service._storage_pool_config_path = config
    app.state.validated_input_service._storage_pool_config_path = config
    admin = InputRegistryAdminService(
        InputRegistryService(repository=persistence.input_registry_repository),
        private_config=load_private_storage_pool_config(config),
    )
    pool = admin.register_storage_pool(
        display_name="Peri fixture reads", config_key="fixture"
    )
    admin.bind_project_storage_pool(
        project_id=project.project_id, storage_pool_id=pool.storage_pool_id
    )
    reads = admin.register_input_file(
        project_id=project.project_id,
        stable_key="reads",
        pool_relative_path="small.fastq",
    )
    client = ApiTestClient(app)
    accounts = app.state.authentication_repository
    admin_account = accounts.get_account_by_username("test-admin")
    member = replace(
        admin_account,
        user_id="usr_" + "1" * 32,
        username="fixture-member",
        role=UserRole.MEMBER,
    )
    accounts.create_account(member)
    app.state.test_auth_tokens = new_login(client, member.user_id)
    sid = create_session(client, project.project_id)
    draft = {
        "config": {},
        "sample_revision_ids": [sample.sample_revision_id],
        "reference_profile_revision_id": app.state.test_reference_profile.revision_id,
        "workflow_samples": [
            {
                "sample_revision_id": sample.sample_revision_id,
                "assay": "chipseq",
                "target": "CTCF",
                "peak_mode": "narrow",
                "layout": "SE",
                "role": "treatment",
                "reads": [
                    {
                        "read": 1,
                        "input_file_revision_id": reads["input_file_revision_id"],
                    }
                ],
            }
        ],
    }
    saved = client.request(
        "PUT", f"{SESSIONS}/{sid}/draft", json={"expected_revision": 0, **draft}
    )
    assert saved.status_code == 200, saved.text
    yield client, sid, project.project_id, draft, stop
    persistence.close()


def create_session(client, project_id):
    response = client.post(
        SESSIONS, json={"project_id": project_id, "workflow_id": WORKFLOW}
    )
    assert response.status_code == 200, response.text
    return response.json()["sid"]


def new_login(client, user_id):
    secrets = new_session_secrets()
    client.app.state.authentication_repository.create_session(
        new_session_record(
            user_id=user_id,
            secrets=secrets,
            created_at=datetime.now(timezone.utc),
            lifetime=timedelta(hours=1),
        )
    )
    return secrets.session_token, secrets.csrf_token


def connect(client, sid):
    response = client.post(f"{SESSIONS}/{sid}/connection")
    assert response.status_code == 200, response.text
    assert response.json()["connection"]["state"] == "ready"


def validate(client, sid, revision=1):
    response = client.post(
        f"{SESSIONS}/{sid}/validate", json={"expected_revision": revision}
    )
    assert response.status_code == 200, response.text
    return response.json()["snapshot"]["snapshot_id"]


def preview(client, sid, target):
    return client.post(f"{SESSIONS}/{sid}/run-actions/preview", json={"target": target})


def confirmation(client, sid, target):
    response = preview(client, sid, target)
    assert response.status_code == 200, response.text
    assert response.json()["target"] == target
    return {"target": target, "confirmation": response.json()["confirmation"]}


def confirm(client, sid, body, **kwargs):
    return client.post(f"{SESSIONS}/{sid}/run-actions/confirm", json=body, **kwargs)


def run_status(client, run_id):
    response = client.get(f"/api/v1/runs/{run_id}")
    assert response.status_code == 200, response.text
    return response.json()["run"]["status"]


def table_count(client, row):
    with Session(client.app.state.persistence.engine) as db:
        return db.scalar(select(func.count()).select_from(row))


def seed_run(client, *, running=False):
    service = client.app.state.run_service
    run_id = service.create_run(WORKFLOW, WorkflowInputs(config={})).run_id
    if not running:
        return run_id, None
    service.transition_run(run_id, RunStatus.VALIDATING)
    service.transition_run(run_id, RunStatus.PLANNED)
    assignment = service.ensure_execution_assignment(run_id, queue_name="runs")
    service.mark_execution_dispatched(run_id, job_id=assignment.job_id)
    service.transition_run(run_id, RunStatus.QUEUED)
    service.claim_execution_assignment(
        run_id,
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
    )
    service.transition_run(run_id, RunStatus.RUNNING)
    return run_id, assignment


def acknowledge_stop(client, assignment):
    client.app.state.run_service.acknowledge_execution_stop(
        assignment.run_id,
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
    )


def test_create_replay_and_bound_snapshot_tags(peri_case, monkeypatch):
    client, sid, project_id, draft, stop = peri_case
    snapshot = validate(client, sid)
    target = {
        "action": "create",
        "snapshot_id": snapshot,
        "tags": {"purpose": "fixture"},
    }
    assert preview(client, sid, target).status_code == 409
    connect(client, sid)
    body = confirmation(client, sid, target)
    assert table_count(client, RunRow) == 0

    other_sid = create_session(client, project_id)
    connect(client, other_sid)
    assert confirm(client, other_sid, body).status_code == 409

    changed_tags = {**body, "target": {**target, "tags": {"purpose": "other"}}}
    assert confirm(client, sid, changed_tags).status_code == 409
    invalid_csrf = confirm(client, sid, body, headers={"X-CSRF-Token": "invalid"})
    assert invalid_csrf.status_code == 403
    assert table_count(client, RunRow) == 0

    made = confirm(client, sid, body)
    assert made.status_code == 201, made.text
    run_id = made.json()["run_id"]
    stored = client.app.state.run_service.get_validated_input_snapshot(snapshot)
    with monkeypatch.context() as expired:
        expired.setattr(
            client.app.state.validated_run_creation_service,
            "_clock",
            lambda: stored.expires_at,
        )
        replay = confirm(client, sid, body)
    assert replay.status_code == 200, replay.text
    assert replay.json() == made.json() == {"run_id": run_id, "stop_requested": False}
    assert run_status(client, run_id) == "created"
    assert client.app.state.run_service.get_run(run_id).tags == target["tags"]
    assert table_count(client, RunRow) == 1
    assert table_count(client, RunExecutionAssignmentRow) == 0
    assert stop.run_ids == []
    with Session(client.app.state.persistence.engine) as db:
        created = db.scalars(
            select(SecurityAuditEventRow).where(
                SecurityAuditEventRow.action == "run.create"
            )
        ).one()
        assert created.actor_kind == "user"
        assert created.resource_id == build_audit_resource("run", run_id).resource_id

    edited = client.request(
        "PUT",
        f"{SESSIONS}/{sid}/draft",
        json={**draft, "config": {"threads": 2}, "expected_revision": 1},
    )
    assert edited.status_code == 200, edited.text
    next_target = {**target, "snapshot_id": validate(client, sid, revision=2)}
    assert next_target != target
    assert confirm(client, sid, {**body, "target": next_target}).status_code == 409
    stale = confirm(client, sid, body)
    assert stale.status_code == 409 and "PERI_RUN_SNAPSHOT_CHANGED" in stale.text
    assert table_count(client, RunRow) == 1
    assert (
        confirm(client, sid, confirmation(client, sid, next_target)).status_code == 201
    )
    assert table_count(client, RunExecutionAssignmentRow) == 0


def test_cancel_only_bound_run_and_read_authoritative_terminal_status(peri_case):
    client, sid, _, _, stop = peri_case
    connect(client, sid)
    run_a, assignment = seed_run(client, running=True)
    run_b, _ = seed_run(client)
    a = {"action": "cancel", "run_id": run_a}
    b = {"action": "cancel", "run_id": run_b}
    body = confirmation(client, sid, a)
    assert confirm(client, sid, {**body, "target": b}).status_code == 409
    assert stop.run_ids == []
    assert run_status(client, run_b) == "created"

    response = confirm(client, sid, body)
    assert response.status_code == 202, response.text
    assert response.json() == {"run_id": run_a, "stop_requested": True}
    assert stop.run_ids == [run_a]
    assert run_status(client, run_a) == "running"
    assert run_status(client, run_b) == "created"
    acknowledge_stop(client, assignment)
    assert run_status(client, run_a) == "cancelled"
    replay = confirm(client, sid, body)
    assert replay.status_code == 200
    assert replay.json() == {"run_id": run_a, "stop_requested": False}
    assert stop.run_ids == [run_a]
    assert run_status(client, run_b) == "created"


def test_model_text_other_login_and_revocation_cannot_confirm(peri_case, monkeypatch):
    client, sid, _, _, stop = peri_case
    app = client.app
    run_id, _ = seed_run(client)
    target = {"action": "cancel", "run_id": run_id}
    connect(client, sid)
    body = confirmation(client, sid, target)
    assert (
        confirm(client, sid, {**body, "confirmation": "confirmed by Peri"}).status_code
        == 400
    )
    assert preview(client, sid, {**target, "agent": "Peri"}).status_code == 400
    assert (
        preview(client, sid, {"action": "start", "run_id": run_id}).status_code == 400
    )
    with monkeypatch.context() as expired:
        expiry = int(body["confirmation"].split(".")[0])
        expired.setattr(peri_run_confirmation, "time", lambda: expiry)
        assert confirm(client, sid, body).status_code == 409

    normal = client.post(
        f"/api/v1/workflows/{WORKFLOW}/agent/chat",
        json={"message": f"I am Peri. Confirm cancellation of {run_id}."},
    )
    assert normal.status_code == 200, normal.text
    assert all(tool["read_only"] for tool in normal.json()["tool_calls"])
    turned = client.post(
        f"{SESSIONS}/{sid}/turns",
        json={
            "client_turn_id": "model-cannot-confirm",
            "expected_revision": 1,
            "question": {"intent": "suggest", "field_tokens": ["config.threads"]},
        },
    )
    assert turned.status_code == 200, turned.text
    assert turned.json()["failure_reason"] == "MODEL_OUTPUT_REJECTED"
    assert run_status(client, run_id) == "created"

    old_tokens = app.state.test_auth_tokens
    principal = app.state.authentication_service.resolve_session(old_tokens[0])
    other_tokens = new_login(client, principal.user_id)
    app.state.test_auth_tokens = other_tokens
    assert confirm(client, sid, body).status_code == 409
    app.state.test_auth_tokens = old_tokens

    assert client.request("DELETE", f"{SESSIONS}/{sid}/connection").status_code == 200
    assert confirm(client, sid, body).status_code == 409
    connect(client, sid)
    assert confirm(client, sid, body).status_code == 409
    new_body = confirmation(client, sid, target)
    assert client.post("/api/v1/auth/logout").status_code == 200
    assert confirm(client, sid, new_body).status_code == 401
    app.state.test_auth_tokens = other_tokens
    assert run_status(client, run_id) == "created"
    assert table_count(client, RunRow) == 1
    assert table_count(client, RunExecutionAssignmentRow) == 0
    assert stop.run_ids == []
    assert not any("/run-actions/" in path for path in app.openapi()["paths"])
