"""API identity and validation wiring for restricted authoring sessions."""

from __future__ import annotations

import json

from encode_pipeline.persistence.data_registry import SqlAlchemyDataRegistryRepository
from encode_pipeline.services.authoring_bridge import BridgeAdmission, BridgeStatus
from encode_pipeline.services.data_registry import DataRegistryService, SampleImportSpec
from api_test_client import ApiTestClient


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


def test_authoring_routes_require_a_login_and_hide_other_owners(reference_ready_app):
    app = reference_ready_app
    saved_tokens = app.state.test_auth_tokens
    app.state.test_auth_tokens = None
    with ApiTestClient(app) as unauthenticated:
        denied = unauthenticated.get("/api/v1/agent/authoring-sessions")
    app.state.test_auth_tokens = saved_tokens
    assert denied.status_code == 401, denied.text[:500]

    registry = _registry(app)
    created_project = registry.create_project("Authoring")
    with ApiTestClient(app) as client:
        created = client.post(
            "/api/v1/agent/authoring-sessions",
            json={
                "project_id": created_project.project_id,
                "workflow_id": WORKFLOW_ID,
            },
        )
        assert created.status_code == 200
        body = created.json()
        assert body["draft_revision"] == 0
        assert body["draft_state"] == "active"
        assert body["draft"]["config"] == {}
        assert body["authoring_capability"]["available"] is True
        assert body["connection"]["state"] == "closed"
        assert body["requires_new_connection"] is True
        tokens = {field["token"] for field in body["authoring_capability"]["fields"]}
        assert tokens == {"config.threads", "options.cores"}
        sid = body["sid"]

        unknown = client.request(
            "PUT",
            f"/api/v1/agent/authoring-sessions/{sid}/draft",
            json={"expected_revision": 0, "config": {}, "unexpected": True},
        )
        assert unknown.status_code == 422
        assert unknown.json()["reason_code"] == "INVALID_REQUEST"
        assert "Traceback" not in unknown.text

        opened = client.post(f"/api/v1/agent/authoring-sessions/{sid}/connection")
        assert opened.status_code == 503
        assert opened.json()["reason_code"] == "BRIDGE_UNAVAILABLE"

        validated = client.post(
            f"/api/v1/agent/authoring-sessions/{sid}/validate",
            json={"expected_revision": 0},
        )
        assert validated.status_code == 200
        payload = validated.json()
        assert payload["ok"] is False
        assert payload["snapshot"] is None
        assert payload["issues"][0]["code"] == "SAMPLE_SELECTION_REQUIRED"

        oversized = client.request(
            "PUT",
            f"/api/v1/agent/authoring-sessions/{sid}/draft",
            content=b'{"config":"' + (b"x" * (256 * 1024)) + b'"}',
            headers={"content-type": "application/json"},
        )
        assert oversized.status_code == 413
        assert oversized.json()["reason_code"] == "REQUEST_TOO_LARGE"
        assert "x" * 40 not in oversized.text


def test_validate_uses_the_real_service_for_failure_null_and_managed_rejection(
    reference_ready_app, tmp_path
):
    app = reference_ready_app
    registry = _registry(app)
    project = registry.create_project("Validator")
    imported = registry.import_samples(
        project.project_id,
        rows=(
            SampleImportSpec(
                stable_key="sample-a",
                display_name="Sample A",
                attributes={},
            ),
        ),
    )
    revision_id = imported[0].revision.sample_revision_id
    reference_id = app.state.test_reference_profile.revision_id
    fastq_1 = tmp_path / "S1.R1.fastq.gz"
    fastq_2 = tmp_path / "S1.R2.fastq.gz"
    with ApiTestClient(app) as client:
        encode = client.post(
            "/api/v1/agent/authoring-sessions",
            json={"project_id": project.project_id, "workflow_id": WORKFLOW_ID},
        ).json()
        invalid = client.request(
            "PUT",
            f"/api/v1/agent/authoring-sessions/{encode['sid']}/draft",
            json={
                "expected_revision": 0,
                "config": {"threads": 0},
                "options": {},
                "sample_revision_ids": [revision_id],
                "reference_profile_revision_id": reference_id,
            },
        )
        assert invalid.status_code == 200, invalid.text
        checked = client.post(
            f"/api/v1/agent/authoring-sessions/{encode['sid']}/validate",
            json={"expected_revision": 1},
        )
        assert checked.status_code == 200, checked.text
        assert checked.json()["ok"] is False
        assert checked.json()["snapshot"] is None
        assert checked.json()["issues"][0]["code"] != "SAMPLE_SELECTION_REQUIRED"

        managed = client.request(
            "PUT",
            f"/api/v1/agent/authoring-sessions/{encode['sid']}/draft",
            json={
                "expected_revision": 1,
                "config": {"threads": 1},
                "sample_revision_ids": [revision_id],
                "reference_profile_revision_id": reference_id,
                "input_selections": [
                    {
                        "input_use_key": "reads",
                        "occurrence": 0,
                        "input_file_revision_ids": ["inpfr_" + "a" * 32],
                    }
                ],
            },
        )
        assert managed.status_code == 200, managed.text
        rejected = client.post(
            f"/api/v1/agent/authoring-sessions/{encode['sid']}/validate",
            json={"expected_revision": managed.json()["draft_revision"]},
        )
        assert rejected.json()["ok"] is False
        assert rejected.json()["snapshot"] is None
        assert rejected.json()["issues"][0]["code"] in {
            "INPUT_USE_CAPABILITY_UNAVAILABLE",
            "INPUT_BINDING_SELECTION_INVALID",
            "ENCODE_REFERENCE_BINDING_CONFLICT",
        }

        hitrac = client.post(
            "/api/v1/agent/authoring-sessions",
            json={"project_id": project.project_id, "workflow_id": "hitrac-preprocess"},
        ).json()
        assert hitrac["authoring_capability"]["available"] is False
        saved = client.request(
            "PUT",
            f"/api/v1/agent/authoring-sessions/{hitrac['sid']}/draft",
            json={
                "expected_revision": 0,
                "config": {},
                "samples": [
                    {
                        "sample_id": "S1",
                        "fastq_1": str(fastq_1),
                        "fastq_2": str(fastq_2),
                    }
                ],
                "sample_revision_ids": [revision_id],
            },
        )
        assert saved.status_code == 200, saved.text
        authoring_only = client.post(
            f"/api/v1/agent/authoring-sessions/{hitrac['sid']}/validate",
            json={"expected_revision": saved.json()["draft_revision"]},
        )
        assert authoring_only.status_code == 409, authoring_only.text
        assert authoring_only.json()["snapshot"] is None
        assert (
            authoring_only.json()["issues"][0]["code"] == "REFERENCE_PROFILE_REQUIRED"
        )


def test_public_turn_digest_adopts_rejects_and_replays(reference_ready_app):
    class Bridge:
        def status(self):
            return BridgeStatus(available=True)

        def owner_epoch(self):
            return "api-owner"

        def admit(self, request):
            del request
            return BridgeAdmission(instance_id="instance01", state="ready")

        def execute_turn(self, request):
            value = 4 if request.question["intent"] == "suggest" else None
            if value is None:
                body = {
                    "base_revision": request.base_revision,
                    "changes": [],
                    "kind": "explanation",
                    "reason_codes": [],
                    "schema_digest": request.schema_digest,
                }
            else:
                body = {
                    "base_revision": request.base_revision,
                    "changes": [{"field_token": "config.threads", "value": value}],
                    "kind": "proposal",
                    "reason_codes": [],
                    "schema_digest": request.schema_digest,
                }
            return json.dumps(body, separators=(",", ":"), sort_keys=True)

        def close(self, request):
            del request

    app = reference_ready_app
    app.state.authoring_service._bridge = Bridge()
    project = _registry(app).create_project("Digest")
    with ApiTestClient(app) as client:
        created = client.post(
            "/api/v1/agent/authoring-sessions",
            json={"project_id": project.project_id, "workflow_id": WORKFLOW_ID},
        )
        assert created.status_code == 200, created.text
        sid = created.json()["sid"]
        opened = client.post(f"/api/v1/agent/authoring-sessions/{sid}/connection")
        assert opened.status_code == 200, opened.text
        suggested = client.post(
            f"/api/v1/agent/authoring-sessions/{sid}/turns",
            json={
                "client_turn_id": "suggest-1",
                "expected_revision": 0,
                "question": {
                    "intent": "suggest",
                    "field_tokens": ["config.threads"],
                    "preferences": [{"field_token": "config.threads", "value": 4}],
                },
            },
        )
        assert suggested.status_code == 200, suggested.text
        body = suggested.json()
        digest = body["proposal_digest"]
        proposal_id = body["proposal_id"]
        assert body["adoptable"] is True
        assert digest and "proposal_digest" not in body["proposal"]
        adopted = client.post(
            f"/api/v1/agent/authoring-sessions/{sid}/proposals/{proposal_id}/adopt",
            json={"expected_revision": 0, "proposal_digest": digest},
        )
        assert adopted.status_code == 200, adopted.text
        assert adopted.json()["replayed"] is False
        replayed = client.post(
            f"/api/v1/agent/authoring-sessions/{sid}/proposals/{proposal_id}/adopt",
            json={"expected_revision": 0, "proposal_digest": digest},
        )
        assert replayed.status_code == 200, replayed.text
        assert replayed.json()["replayed"] is True
        assert replayed.json()["adoption_receipt"] == adopted.json()["adoption_receipt"]
        again = client.post(
            f"/api/v1/agent/authoring-sessions/{sid}/turns",
            json={
                "client_turn_id": "suggest-2",
                "expected_revision": adopted.json()["draft_revision"],
                "question": {
                    "intent": "suggest",
                    "field_tokens": ["options.cores"],
                    "preferences": [{"field_token": "options.cores", "value": 2}],
                },
            },
        )
        assert again.status_code == 200, again.text
        rejected = client.post(
            f"/api/v1/agent/authoring-sessions/{sid}/proposals/"
            f"{again.json()['proposal_id']}/reject",
            json={"proposal_digest": again.json()["proposal_digest"]},
        )
        assert rejected.status_code == 200, rejected.text
        assert rejected.json()["proposal_state"] == "rejected"
        assert rejected.json()["adoptable"] is False


def test_owner_rediscovers_adopted_receipts_after_a_new_login(reference_ready_app):
    from datetime import datetime, timedelta, timezone

    from encode_pipeline.platform.authentication import (
        UserAccount,
        UserRole,
        UserStatus,
    )
    from conftest import TEST_AUTH_PASSWORD_HASH, TEST_AUTH_USER_ID

    class Clock:
        def __init__(self) -> None:
            self.now = datetime.now(timezone.utc).replace(microsecond=0) + timedelta(
                minutes=1
            )

        def __call__(self) -> datetime:
            return self.now

    class Bridge:
        def __init__(self) -> None:
            self.admits = 0
            self.turns = 0

        def status(self):
            return BridgeStatus(available=True)

        def owner_epoch(self) -> str:
            return "owner-epoch-0001"

        def admit(self, request):
            del request
            self.admits += 1
            return BridgeAdmission(instance_id="instance01", state="ready")

        def execute_turn(self, request) -> str:
            self.turns += 1
            token = request.question["preferences"][0]["field_token"]
            value = request.question["preferences"][0]["value"]
            body = {
                "base_revision": request.base_revision,
                "changes": [{"field_token": token, "value": value}],
                "kind": "proposal",
                "reason_codes": [],
                "schema_digest": request.schema_digest,
            }
            return json.dumps(body, separators=(",", ":"), sort_keys=True)

        def close(self, request) -> None:
            del request

    app = reference_ready_app
    clock = Clock()
    service = app.state.authoring_service
    service._clock = clock
    service._repository._clock = clock
    app.state.authentication_service._now_factory = clock
    bridge = Bridge()
    service._bridge = bridge
    renewed = ("c" * 43, "d" * 43)
    _open_login(app, TEST_AUTH_USER_ID, renewed)
    project = _registry(app).create_project("Receipts")
    with ApiTestClient(app) as client:
        created = client.post(
            "/api/v1/agent/authoring-sessions",
            json={"project_id": project.project_id, "workflow_id": WORKFLOW_ID},
        )
        assert created.status_code == 200, created.text
        sid = created.json()["sid"]
        opened = client.post(f"/api/v1/agent/authoring-sessions/{sid}/connection")
        assert opened.status_code == 200, opened.text
        suggested = client.post(
            f"/api/v1/agent/authoring-sessions/{sid}/turns",
            json={
                "client_turn_id": "adopted-turn",
                "expected_revision": 0,
                "question": {
                    "intent": "suggest",
                    "field_tokens": ["config.threads"],
                    "preferences": [{"field_token": "config.threads", "value": 4}],
                },
            },
        )
        assert suggested.status_code == 200, suggested.text
        adopted = client.post(
            f"/api/v1/agent/authoring-sessions/{sid}/proposals/"
            f"{suggested.json()['proposal_id']}/adopt",
            json={
                "expected_revision": 0,
                "proposal_digest": suggested.json()["proposal_digest"],
            },
        )
        assert adopted.status_code == 200, adopted.text
        receipt = adopted.json()["adoption_receipt"]
        hidden = client.post(
            f"/api/v1/agent/authoring-sessions/{sid}/turns",
            json={
                "client_turn_id": "hidden-turn",
                "expected_revision": adopted.json()["draft_revision"],
                "question": {
                    "intent": "suggest",
                    "field_tokens": ["options.cores"],
                    "preferences": [{"field_token": "options.cores", "value": 9}],
                },
            },
        )
        assert hidden.status_code == 200, hidden.text
        before = client.get(f"/api/v1/agent/authoring-sessions/{sid}").json()
        admits = bridge.admits
        turns = bridge.turns
        app.state.test_auth_tokens = renewed
        found = client.get(f"/api/v1/agent/authoring-sessions/{sid}/receipts")
        assert found.status_code == 200, found.text
        body = found.json()
        assert body["draft_state"] == "active"
        assert body["next_cursor"] is None
        assert len(body["receipts"]) == 1
        assert set(body["receipts"][0]) == {"turn_id", "adoption_receipt"}
        assert body["receipts"][0]["adoption_receipt"] == receipt
        assert body["receipts"][0]["turn_id"] == suggested.json()["turn_id"]
        assert "hidden-turn" not in found.text
        assert "config.threads" not in found.text
        assert "originating" not in found.text
        session = client.get(f"/api/v1/agent/authoring-sessions/{sid}")
        assert session.status_code == 200, session.text
        assert session.json()["requires_new_connection"] is True
        assert session.json()["connection"]["state"] == "ready"
        assert session.json()["draft_revision"] == before["draft_revision"]
        assert session.json()["expires_at"] == before["expires_at"]
        assert session.json()["purge_after"] == before["purge_after"]
        assert bridge.admits == admits
        assert bridge.turns == turns
    original_tokens = ("a" * 43, "b" * 43)
    app.state.test_auth_tokens = original_tokens
    with ApiTestClient(app) as client:
        replayed = client.post(
            f"/api/v1/agent/authoring-sessions/{sid}/proposals/"
            f"{receipt['proposal_id']}/adopt",
            json={
                "expected_revision": receipt["base_revision"],
                "proposal_digest": suggested.json()["proposal_digest"],
            },
        )
        assert replayed.status_code == 200, replayed.text
        assert replayed.json()["replayed"] is True
        assert replayed.json()["adoption_receipt"] == receipt
        assert replayed.json()["draft_revision"] == before["draft_revision"]
        bad = client.get(
            f"/api/v1/agent/authoring-sessions/{sid}/receipts",
            params={"cursor": "not-a-cursor"},
        )
        assert bad.status_code == 422
        assert bad.json()["reason_code"] == "INVALID_REQUEST"
        missing = client.get(
            "/api/v1/agent/authoring-sessions/aas_" + "ab" * 16 + "/receipts"
        )
        assert missing.status_code == 404
        logged_out = client.post("/api/v1/auth/logout")
        assert logged_out.status_code == 200, logged_out.text
        app.state.test_auth_tokens = renewed
        still = client.get(f"/api/v1/agent/authoring-sessions/{sid}/receipts")
        assert still.status_code == 200, still.text
        assert still.json()["receipts"][0]["adoption_receipt"] == receipt
        assert still.json()["draft_state"] == "active"
        outsider = ("e" * 43, "f" * 43)
        other = UserAccount(
            user_id="usr_" + "ab" * 16,
            username="receipt-other",
            role=UserRole.MEMBER,
            status=UserStatus.ENABLED,
            password_hash=TEST_AUTH_PASSWORD_HASH,
            created_at=clock.now,
            updated_at=clock.now,
            password_changed_at=clock.now,
        )
        app.state.authentication_repository.create_account(other)
        _open_login(app, other.user_id, outsider)
        app.state.test_auth_tokens = outsider
        assert (
            client.get(f"/api/v1/agent/authoring-sessions/{sid}/receipts").status_code
            == 404
        )
        archived_project = _registry(app).create_project("Archived receipts")
        app.state.test_auth_tokens = renewed
        archived = client.post(
            "/api/v1/agent/authoring-sessions",
            json={
                "project_id": archived_project.project_id,
                "workflow_id": WORKFLOW_ID,
            },
        )
        assert archived.status_code == 200, archived.text
        _registry(app).archive_project(archived_project.project_id)
        assert (
            client.get(
                f"/api/v1/agent/authoring-sessions/{archived.json()['sid']}/receipts"
            ).status_code
            == 404
        )

        clock.now = _parse(before["expires_at"]) - timedelta(seconds=1)
        _cover_clock(app, renewed[0], clock.now)
        active = client.get(f"/api/v1/agent/authoring-sessions/{sid}/receipts")
        assert active.status_code == 200, active.text
        assert active.json()["draft_state"] == "active"
        clock.now = _parse(before["expires_at"])
        expired = client.get(f"/api/v1/agent/authoring-sessions/{sid}/receipts")
        assert expired.status_code == 200, expired.text
        assert expired.json()["draft_state"] == "expired"
        assert expired.json()["receipts"][0]["adoption_receipt"] == receipt
        unchanged = client.get(f"/api/v1/agent/authoring-sessions/{sid}")
        assert unchanged.json()["draft"] is None
        assert unchanged.json()["draft_revision"] == before["draft_revision"]
        assert unchanged.json()["expires_at"] == before["expires_at"]
        clock.now = _parse(before["purge_after"]) - timedelta(seconds=1)
        _cover_clock(app, renewed[0], clock.now)
        assert (
            client.get(f"/api/v1/agent/authoring-sessions/{sid}/receipts").status_code
            == 200
        )
        clock.now = _parse(before["purge_after"])
        _cover_clock(app, renewed[0], clock.now)
        purged = client.get(f"/api/v1/agent/authoring-sessions/{sid}/receipts")
        assert purged.status_code == 404


def _parse(value: str):
    from datetime import datetime, timezone

    return datetime.strptime(value, "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)


def _cover_clock(app, token: str, now) -> None:
    from datetime import timedelta

    from encode_pipeline.persistence.models import AuthSessionRow
    from encode_pipeline.services.authentication import digest_session_token

    digest = digest_session_token(token)
    with app.state.authentication_repository._session_factory.begin() as session:
        row = session.get(AuthSessionRow, digest)
        assert row is not None
        row.created_at = now - timedelta(minutes=1)
        row.expires_at = now + timedelta(days=1)


def _open_login(app, user_id: str, tokens: tuple[str, str]) -> None:
    from datetime import datetime, timedelta, timezone

    from encode_pipeline.services.authentication import (
        SessionSecrets,
        new_session_record,
    )

    app.state.authentication_repository.create_session(
        new_session_record(
            user_id=user_id,
            secrets=SessionSecrets._from_generated(*tokens),
            created_at=datetime.now(timezone.utc),
            lifetime=timedelta(days=7),
        )
    )


def _session_factory(app):
    from encode_pipeline.persistence import create_session_factory

    return create_session_factory(app.state.persistence.engine)


def _registry(app) -> DataRegistryService:
    return DataRegistryService(
        repository=SqlAlchemyDataRegistryRepository(_session_factory(app))
    )
