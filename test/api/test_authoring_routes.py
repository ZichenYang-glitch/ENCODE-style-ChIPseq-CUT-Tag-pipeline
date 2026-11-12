"""API identity and validation wiring for restricted authoring sessions."""

from __future__ import annotations

import json

import pytest

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
        assert rejected.json()["issues"][0]["code"] == "WORKFLOW_SAMPLE_PLAN_REQUIRED"

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


def test_managed_sample_plan_validates_to_a_real_snapshot(
    tmp_path, reference_ready_app
):
    import json

    from encode_pipeline.services.input_registry import (
        InputRegistryAdminService,
        InputRegistryService,
    )
    from encode_pipeline.services.private_storage_pools import (
        PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
        load_private_storage_pool_config,
    )

    app = reference_ready_app
    root = (tmp_path / "pool").resolve()
    (root / "reads").mkdir(parents=True)
    fastq = root / "reads" / "s1.fastq.gz"
    fastq.write_bytes(b"@read\nACGT\n+\n!!!!\n")
    config_path = tmp_path / "pools.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
                "storage_pool_roots": {"authoring-inputs": str(root)},
            }
        ),
        encoding="utf-8",
    )
    app.state.authoring_service._storage_pool_config_path = config_path
    app.state.validated_input_service._storage_pool_config_path = config_path
    registry = _registry(app)
    project = registry.create_project("Managed cohort")
    other = registry.create_project("Other cohort")
    imported = registry.import_samples(
        project.project_id,
        rows=(
            SampleImportSpec(
                stable_key="sample-a",
                display_name="Private sample",
                attributes={"note": "must-not-leak"},
            ),
        ),
    )[0]
    foreign = registry.import_samples(
        other.project_id,
        rows=(
            SampleImportSpec(
                stable_key="sample-b",
                display_name="Foreign",
                attributes={},
            ),
        ),
    )[0]
    admin = InputRegistryAdminService(
        InputRegistryService(
            repository=app.state.persistence.input_registry_repository
        ),
        private_config=load_private_storage_pool_config(config_path),
    )
    pool = admin.register_storage_pool(
        display_name="Authoring inputs",
        config_key="authoring-inputs",
    )
    admin.bind_project_storage_pool(
        project_id=project.project_id,
        storage_pool_id=pool.storage_pool_id,
    )
    registered = admin.register_input_file(
        project_id=project.project_id,
        stable_key="sample-a-r1",
        pool_relative_path="reads/s1.fastq.gz",
    )
    plan = {
        "sample_revision_id": imported.revision.sample_revision_id,
        "assay": "chipseq",
        "target": "CTCF",
        "peak_mode": "narrow",
        "layout": "SE",
        "role": "treatment",
        "control_sample_revision_id": None,
        "reads": [
            {
                "read": 1,
                "input_file_revision_id": registered["input_file_revision_id"],
            }
        ],
    }
    with ApiTestClient(app) as client:
        created = client.post(
            "/api/v1/agent/authoring-sessions",
            json={"project_id": project.project_id, "workflow_id": WORKFLOW_ID},
        )
        assert created.status_code == 200, created.text
        sid = created.json()["sid"]
        saved = client.request(
            "PUT",
            f"/api/v1/agent/authoring-sessions/{sid}/draft",
            json={
                "expected_revision": 0,
                "config": {},
                "samples": None,
                "sample_revision_ids": [imported.revision.sample_revision_id],
                "reference_profile_revision_id": (
                    app.state.test_reference_profile.revision_id
                ),
                "workflow_samples": [plan],
            },
        )
        assert saved.status_code == 200, saved.text
        assert saved.json()["draft"]["samples"] is None
        assert "must-not-leak" not in saved.text
        assert str(root) not in saved.text
        validated = client.post(
            f"/api/v1/agent/authoring-sessions/{sid}/validate",
            json={"expected_revision": saved.json()["draft_revision"]},
        )
        assert validated.status_code == 200, validated.text
        body = validated.json()
        assert body["ok"] is True, body
        assert body["snapshot"]["snapshot_id"].startswith("vsnap_")
        assert body["snapshot"]["sample_revision_ids"] == [
            imported.revision.sample_revision_id
        ]
        binding = app.state.persistence.repository.get_validated_input_use_binding(
            body["snapshot"]["snapshot_id"]
        )
        members = [member for use in binding.input_uses for member in use.members]
        assert [member.input_file_revision_id for member in members] == [
            registered["input_file_revision_id"]
        ]
        assert [member.content_sha256 for member in members] == [
            registered["content_sha256"]
        ]
        assert [use.provenance_mode.value for use in binding.input_uses] == [
            "managed_revision_v1"
        ]
        assert registered["input_file_revision_id"] in validated.text
        assert str(root) not in validated.text
        assert "must-not-leak" not in validated.text
        assert "Private sample" not in validated.text
        ambiguous = dict(plan)
        ambiguous["layout"] = "PE"
        rejected_shape = client.request(
            "PUT",
            f"/api/v1/agent/authoring-sessions/{sid}/draft",
            json={
                "expected_revision": saved.json()["draft_revision"],
                "config": {},
                "sample_revision_ids": [imported.revision.sample_revision_id],
                "reference_profile_revision_id": (
                    app.state.test_reference_profile.revision_id
                ),
                "workflow_samples": [ambiguous],
            },
        )
        assert rejected_shape.status_code == 422
        assert str(root) not in rejected_shape.text
        fastq.write_bytes(b"changed")
        stale = client.post(
            f"/api/v1/agent/authoring-sessions/{sid}/validate",
            json={"expected_revision": saved.json()["draft_revision"]},
        )
        assert stale.json()["ok"] is False
        assert stale.json()["snapshot"] is None
        assert stale.json()["issues"][0]["code"] == "WORKFLOW_SAMPLE_INPUT_UNAVAILABLE"
        assert str(root) not in stale.text
        foreign_plan = dict(plan)
        foreign_plan["sample_revision_id"] = foreign.revision.sample_revision_id
        denied = client.request(
            "PUT",
            f"/api/v1/agent/authoring-sessions/{sid}/draft",
            json={
                "expected_revision": saved.json()["draft_revision"],
                "config": {},
                "sample_revision_ids": [foreign.revision.sample_revision_id],
                "reference_profile_revision_id": (
                    app.state.test_reference_profile.revision_id
                ),
                "workflow_samples": [foreign_plan],
            },
        )
        assert denied.status_code == 200, denied.text
        rejected = client.post(
            f"/api/v1/agent/authoring-sessions/{sid}/validate",
            json={"expected_revision": denied.json()["draft_revision"]},
        )
        assert rejected.json()["ok"] is False
        assert rejected.json()["snapshot"] is None
        assert rejected.json()["issues"][0]["code"] == "WORKFLOW_SAMPLE_PLAN_INVALID"
        archived_path = root / "reads" / "archived.fastq.gz"
        archived_path.write_bytes(b"@read\nAAAA\n+\n!!!!\n")
        archived_input = admin.register_input_file(
            project_id=project.project_id,
            stable_key="sample-a-archived",
            pool_relative_path="reads/archived.fastq.gz",
        )
        InputRegistryService(
            repository=app.state.persistence.input_registry_repository
        ).archive_input_file(str(archived_input["input_file_id"]))
        archived_plan = dict(plan)
        archived_plan["reads"] = [
            {
                "read": 1,
                "input_file_revision_id": archived_input["input_file_revision_id"],
            }
        ]
        archived_save = client.request(
            "PUT",
            f"/api/v1/agent/authoring-sessions/{sid}/draft",
            json={
                "expected_revision": denied.json()["draft_revision"],
                "config": {},
                "samples": None,
                "sample_revision_ids": [imported.revision.sample_revision_id],
                "reference_profile_revision_id": (
                    app.state.test_reference_profile.revision_id
                ),
                "workflow_samples": [archived_plan],
            },
        )
        assert archived_save.status_code == 200, archived_save.text
        archived = client.post(
            f"/api/v1/agent/authoring-sessions/{sid}/validate",
            json={"expected_revision": archived_save.json()["draft_revision"]},
        )
        assert archived.json()["ok"] is False
        assert archived.json()["snapshot"] is None
        assert archived.json()["issues"][0]["code"] == "WORKFLOW_SAMPLE_PLAN_INVALID"
        assert str(root) not in archived.text


def test_managed_sample_plan_rejects_replacement_before_snapshot(
    tmp_path, reference_ready_app, monkeypatch
):
    import json

    from encode_pipeline.services.input_file_access import InputFileAccess
    from encode_pipeline.services.input_registry import (
        InputRegistryAdminService,
        InputRegistryService,
    )
    from encode_pipeline.services.private_storage_pools import (
        PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
        load_private_storage_pool_config,
    )

    app = reference_ready_app
    root = (tmp_path / "pool").resolve()
    (root / "reads").mkdir(parents=True)
    fastq = root / "reads" / "s1.fastq.gz"
    fastq.write_bytes(b"@read\nACGT\n+\n!!!!\n")
    replacement = root / "reads" / "unregistered.fastq.gz"
    replacement.write_bytes(b"@different\nTGCA\n+\n!!!!\n")
    config_path = tmp_path / "pools.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
                "storage_pool_roots": {"authoring-inputs": str(root)},
            }
        ),
        encoding="utf-8",
    )
    app.state.authoring_service._storage_pool_config_path = config_path
    app.state.validated_input_service._storage_pool_config_path = config_path
    registry = _registry(app)
    project = registry.create_project("Managed cohort")
    imported = registry.import_samples(
        project.project_id,
        rows=(
            SampleImportSpec(
                stable_key="sample-a",
                display_name="Private sample",
                attributes={},
            ),
        ),
    )[0]
    admin = InputRegistryAdminService(
        InputRegistryService(
            repository=app.state.persistence.input_registry_repository
        ),
        private_config=load_private_storage_pool_config(config_path),
    )
    pool = admin.register_storage_pool(
        display_name="Authoring inputs",
        config_key="authoring-inputs",
    )
    admin.bind_project_storage_pool(
        project_id=project.project_id,
        storage_pool_id=pool.storage_pool_id,
    )
    registered = admin.register_input_file(
        project_id=project.project_id,
        stable_key="sample-a-r1",
        pool_relative_path="reads/s1.fastq.gz",
    )
    original_reverify = InputFileAccess.reverify

    def reverify_then_replace(access, expected):
        observation = original_reverify(access, expected)
        if fastq.is_symlink():
            return observation
        fastq.unlink()
        fastq.symlink_to(replacement)
        return observation

    monkeypatch.setattr(InputFileAccess, "reverify", reverify_then_replace)
    with ApiTestClient(app) as client:
        created = client.post(
            "/api/v1/agent/authoring-sessions",
            json={"project_id": project.project_id, "workflow_id": WORKFLOW_ID},
        )
        sid = created.json()["sid"]
        saved = client.request(
            "PUT",
            f"/api/v1/agent/authoring-sessions/{sid}/draft",
            json={
                "expected_revision": 0,
                "config": {},
                "samples": None,
                "sample_revision_ids": [imported.revision.sample_revision_id],
                "reference_profile_revision_id": (
                    app.state.test_reference_profile.revision_id
                ),
                "workflow_samples": [
                    {
                        "sample_revision_id": imported.revision.sample_revision_id,
                        "assay": "chipseq",
                        "target": "CTCF",
                        "peak_mode": "narrow",
                        "layout": "SE",
                        "role": "treatment",
                        "control_sample_revision_id": None,
                        "reads": [
                            {
                                "read": 1,
                                "input_file_revision_id": registered[
                                    "input_file_revision_id"
                                ],
                            }
                        ],
                    }
                ],
            },
        )
        assert saved.status_code == 200, saved.text
        validated = client.post(
            f"/api/v1/agent/authoring-sessions/{sid}/validate",
            json={"expected_revision": saved.json()["draft_revision"]},
        )
        body = validated.json()
        assert body["ok"] is False, body
        assert body["snapshot"] is None
        assert body["issues"][0]["code"] == "WORKFLOW_SAMPLE_INPUT_UNAVAILABLE"
        assert "unregistered.fastq.gz" not in validated.text
        assert str(root) not in validated.text


def test_public_validation_rejects_path_and_byte_mismatches(
    tmp_path, reference_ready_app
):
    import json

    from encode_pipeline.services.input_registry import (
        InputRegistryAdminService,
        InputRegistryService,
    )
    from encode_pipeline.services.private_storage_pools import (
        PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
        load_private_storage_pool_config,
    )

    app = reference_ready_app
    root = (tmp_path / "pool").resolve()
    (root / "reads").mkdir(parents=True)
    fastq = root / "reads" / "s1.fastq.gz"
    fastq.write_bytes(b"@read\nACGT\n+\n!!!!\n")
    config_path = tmp_path / "pools.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
                "storage_pool_roots": {"authoring-inputs": str(root)},
            }
        ),
        encoding="utf-8",
    )
    app.state.validated_input_service._storage_pool_config_path = config_path
    registry = _registry(app)
    project = registry.create_project("Managed cohort")
    imported = registry.import_samples(
        project.project_id,
        rows=(
            SampleImportSpec(
                stable_key="sample-a",
                display_name="Private sample",
                attributes={},
            ),
        ),
    )[0]
    admin = InputRegistryAdminService(
        InputRegistryService(
            repository=app.state.persistence.input_registry_repository
        ),
        private_config=load_private_storage_pool_config(config_path),
    )
    pool = admin.register_storage_pool(
        display_name="Authoring inputs",
        config_key="authoring-inputs",
    )
    admin.bind_project_storage_pool(
        project_id=project.project_id,
        storage_pool_id=pool.storage_pool_id,
    )
    registered = admin.register_input_file(
        project_id=project.project_id,
        stable_key="sample-a-r1",
        pool_relative_path="reads/s1.fastq.gz",
    )
    other = root / "reads" / "unregistered.fastq.gz"
    other.write_bytes(b"@different\nTGCA\n+\n!!!!\n")

    def post(path):
        return client.post(
            f"/api/v1/workflows/{WORKFLOW_ID}/validate",
            json={
                "config": {},
                "samples": [
                    {
                        "sample": "sample-a",
                        "assay": "chipseq",
                        "target": "CTCF",
                        "layout": "SE",
                        "role": "treatment",
                        "peak_mode": "narrow",
                        "fastq_1": str(path),
                        "fastq_2": "",
                        "control_sample": "",
                    }
                ],
                "project_id": project.project_id,
                "sample_revision_ids": [imported.revision.sample_revision_id],
                "reference_profile_revision_id": (
                    app.state.test_reference_profile.revision_id
                ),
                "input_selections": [
                    {
                        "input_use_key": "sample-a",
                        "occurrence": 1,
                        "input_file_revision_ids": [
                            registered["input_file_revision_id"]
                        ],
                    }
                ],
            },
        )

    with ApiTestClient(app) as client:
        matched = post(fastq)
        assert matched.status_code == 200, matched.text
        assert matched.json()["ok"] is True, matched.text
        snapshot_id = matched.json()["snapshot"]["snapshot_id"]
        binding = app.state.persistence.repository.get_validated_input_use_binding(
            snapshot_id
        )
        members = [member for use in binding.input_uses for member in use.members]
        assert [member.input_file_revision_id for member in members] == [
            registered["input_file_revision_id"]
        ]
        assert [member.content_sha256 for member in members] == [
            registered["content_sha256"]
        ]
        assert str(root) not in matched.text
        mismatched = post(other)
        assert mismatched.status_code == 409, mismatched.text
        assert mismatched.json()["snapshot"] is None
        assert (
            mismatched.json()["issues"][0]["code"] == "INPUT_BINDING_SELECTION_INVALID"
        )
        assert "unregistered.fastq.gz" not in mismatched.text
        fastq.write_bytes(b"@changed\nTTTT\n+\n!!!!\n")
        changed = post(fastq)
        assert changed.status_code == 409, changed.text
        assert changed.json()["snapshot"] is None
        assert changed.json()["issues"][0]["code"] == "INPUT_BINDING_SELECTION_INVALID"
        assert str(root) not in changed.text
        missing = post(root / "missing" / "s1.fastq.gz")
        assert missing.status_code == 409, missing.text
        assert missing.json()["snapshot"] is None
        assert missing.json()["issues"][0]["code"] == "INPUT_BINDING_SELECTION_INVALID"
        assert str(root) not in missing.text
        linked = root / "linked"
        linked.symlink_to(root / "reads", target_is_directory=True)
        via_link = post(linked / "s1.fastq.gz")
        assert via_link.status_code == 409, via_link.text
        assert via_link.json()["snapshot"] is None
        assert via_link.json()["issues"][0]["code"] == "INPUT_BINDING_SELECTION_INVALID"
        assert str(root) not in via_link.text


@pytest.mark.parametrize("scenario", ("rewrite_restore", "replace_restore"))
def test_validator_window_rejects_restored_replacement(
    tmp_path, reference_ready_app, monkeypatch, scenario
):
    import json

    from encode_pipeline.samples import strict
    from encode_pipeline.services.input_registry import (
        InputRegistryAdminService,
        InputRegistryService,
    )
    from encode_pipeline.services.private_storage_pools import (
        PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
        load_private_storage_pool_config,
    )

    app = reference_ready_app
    root = (tmp_path / "pool").resolve()
    (root / "reads").mkdir(parents=True)
    fastq = root / "reads" / "s1.fastq.gz"
    original_bytes = b"@read\nACGT\n+\n!!!!\n"
    fastq.write_bytes(original_bytes)
    config_path = tmp_path / "pools.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
                "storage_pool_roots": {"authoring-inputs": str(root)},
            }
        ),
        encoding="utf-8",
    )
    app.state.validated_input_service._storage_pool_config_path = config_path
    registry = _registry(app)
    project = registry.create_project("Managed cohort")
    imported = registry.import_samples(
        project.project_id,
        rows=(
            SampleImportSpec(
                stable_key="sample-a",
                display_name="Private sample",
                attributes={},
            ),
        ),
    )[0]
    admin = InputRegistryAdminService(
        InputRegistryService(
            repository=app.state.persistence.input_registry_repository
        ),
        private_config=load_private_storage_pool_config(config_path),
    )
    pool = admin.register_storage_pool(
        display_name="Authoring inputs",
        config_key="authoring-inputs",
    )
    admin.bind_project_storage_pool(
        project_id=project.project_id,
        storage_pool_id=pool.storage_pool_id,
    )
    registered = admin.register_input_file(
        project_id=project.project_id,
        stable_key="sample-a-r1",
        pool_relative_path="reads/s1.fastq.gz",
    )
    original_check = strict._check_fastq_exists
    alternate = b"@read\nTGCA\n+\n!!!!\n"

    def rewrite_during_check(path, sample_id, label, **kwargs):
        if path != str(fastq):
            return original_check(path, sample_id, label, **kwargs)
        parked = fastq.with_name("held-original.fastq.gz")
        if scenario == "replace_restore":
            fastq.rename(parked)
        fastq.write_bytes(alternate)
        try:
            return original_check(path, sample_id, label, **kwargs)
        finally:
            if scenario == "replace_restore":
                fastq.unlink()
                parked.rename(fastq)
            else:
                fastq.write_bytes(original_bytes)

    monkeypatch.setattr(strict, "_check_fastq_exists", rewrite_during_check)
    with ApiTestClient(app) as client:
        response = client.post(
            f"/api/v1/workflows/{WORKFLOW_ID}/validate",
            json={
                "config": {},
                "options": {"strict_inputs": True},
                "samples": [
                    {
                        "sample": "sample-a",
                        "assay": "chipseq",
                        "target": "CTCF",
                        "layout": "SE",
                        "role": "treatment",
                        "peak_mode": "narrow",
                        "fastq_1": str(fastq),
                        "fastq_2": "",
                        "control_sample": "",
                    }
                ],
                "project_id": project.project_id,
                "sample_revision_ids": [imported.revision.sample_revision_id],
                "reference_profile_revision_id": (
                    app.state.test_reference_profile.revision_id
                ),
                "input_selections": [
                    {
                        "input_use_key": "sample-a",
                        "occurrence": 1,
                        "input_file_revision_ids": [
                            registered["input_file_revision_id"]
                        ],
                    }
                ],
            },
        )
    assert response.status_code == 409, response.text
    assert response.json()["snapshot"] is None
    assert response.json()["issues"][0]["code"] == "INPUT_BINDING_SELECTION_INVALID"
    assert fastq.read_bytes() == original_bytes
    assert str(root) not in response.text


def test_catalog_ids_validate_to_a_real_snapshot(tmp_path, reference_ready_app):
    import json

    from encode_pipeline.services.input_registry import (
        InputRegistryAdminService,
        InputRegistryService,
    )
    from encode_pipeline.services.private_storage_pools import (
        PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
        load_private_storage_pool_config,
    )

    app = reference_ready_app
    root = (tmp_path / "pool").resolve()
    (root / "reads").mkdir(parents=True)
    first = root / "reads" / "s1.fastq.gz"
    second = root / "reads" / "s2.fastq.gz"
    first.write_bytes(b"@read\nACGT\n+\n!!!!\n")
    second.write_bytes(b"@read\nTGCA\n+\n!!!!\n")
    archived = root / "reads" / "old.fastq.gz"
    archived.write_bytes(b"@read\nAAAA\n+\n!!!!\n")
    config_path = tmp_path / "pools.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
                "storage_pool_roots": {"authoring-inputs": str(root)},
            }
        ),
        encoding="utf-8",
    )
    app.state.authoring_service._storage_pool_config_path = config_path
    app.state.validated_input_service._storage_pool_config_path = config_path
    registry = _registry(app)
    project = registry.create_project("Managed cohort")
    empty = registry.create_project("Empty cohort")
    imported = registry.import_samples(
        project.project_id,
        rows=(
            SampleImportSpec(
                stable_key="sample-a",
                display_name="Private sample",
                attributes={"note": "must-not-leak"},
            ),
        ),
    )[0]
    admin = InputRegistryAdminService(
        InputRegistryService(
            repository=app.state.persistence.input_registry_repository
        ),
        private_config=load_private_storage_pool_config(config_path),
    )
    pool = admin.register_storage_pool(
        display_name="Authoring inputs",
        config_key="authoring-inputs",
    )
    admin.bind_project_storage_pool(
        project_id=project.project_id,
        storage_pool_id=pool.storage_pool_id,
    )
    admin.bind_project_storage_pool(
        project_id=empty.project_id,
        storage_pool_id=pool.storage_pool_id,
    )
    for name, relative in (
        ("sample-a-r1", "reads/s1.fastq.gz"),
        ("sample-a-r2", "reads/s2.fastq.gz"),
        ("sample-a-old", "reads/old.fastq.gz"),
    ):
        admin.register_input_file(
            project_id=project.project_id,
            stable_key=name,
            pool_relative_path=relative,
        )
    old = app.state.persistence.input_registry_repository.list_input_files(
        project.project_id
    )
    archived_file = next(item for item in old if item.stable_key == "sample-a-old")
    InputRegistryService(
        repository=app.state.persistence.input_registry_repository
    ).archive_input_file(archived_file.input_file_id)

    with ApiTestClient(app) as client:
        missing = client.get(
            "/api/v1/projects/prj_0123456789abcdef0123456789abcdef/input-revisions"
        )
        assert missing.status_code == 404, missing.text
        bad_cursor = client.get(
            f"/api/v1/projects/{project.project_id}/input-revisions",
            params={"cursor": "not-a-cursor"},
        )
        assert bad_cursor.status_code == 422, bad_cursor.text
        vacant = client.get(f"/api/v1/projects/{empty.project_id}/input-revisions")
        assert vacant.status_code == 200, vacant.text
        assert vacant.json()["revisions"] == []
        assert vacant.json()["next_cursor"] is None
        first_page = client.get(
            f"/api/v1/projects/{project.project_id}/input-revisions",
            params={"limit": 1},
        )
        assert first_page.status_code == 200, first_page.text
        assert len(first_page.json()["revisions"]) == 1
        assert first_page.json()["next_cursor"]
        assert str(root) not in first_page.text
        assert "reads/" not in first_page.text
        second_page = client.get(
            f"/api/v1/projects/{project.project_id}/input-revisions",
            params={"limit": 1, "cursor": first_page.json()["next_cursor"]},
        )
        assert len(second_page.json()["revisions"]) == 1
        assert second_page.json()["next_cursor"] is None
        listed = first_page.json()["revisions"] + second_page.json()["revisions"]
        assert {item["stable_key"] for item in listed} == {
            "sample-a-r1",
            "sample-a-r2",
        }
        chosen = next(item for item in listed if item["stable_key"] == "sample-a-r1")
        samples = client.get(f"/api/v1/projects/{project.project_id}/samples")
        sample_revision_id = samples.json()["samples"][0]["sample_revision_id"]
        assert sample_revision_id == imported.revision.sample_revision_id
        created = client.post(
            "/api/v1/agent/authoring-sessions",
            json={"project_id": project.project_id, "workflow_id": WORKFLOW_ID},
        )
        sid = created.json()["sid"]
        saved = client.request(
            "PUT",
            f"/api/v1/agent/authoring-sessions/{sid}/draft",
            json={
                "expected_revision": 0,
                "config": {},
                "samples": None,
                "sample_revision_ids": [sample_revision_id],
                "reference_profile_revision_id": (
                    app.state.test_reference_profile.revision_id
                ),
                "workflow_samples": [
                    {
                        "sample_revision_id": sample_revision_id,
                        "assay": "chipseq",
                        "target": "CTCF",
                        "peak_mode": "narrow",
                        "layout": "SE",
                        "role": "treatment",
                        "control_sample_revision_id": None,
                        "reads": [
                            {
                                "read": 1,
                                "input_file_revision_id": chosen[
                                    "input_file_revision_id"
                                ],
                            }
                        ],
                    }
                ],
            },
        )
        assert saved.status_code == 200, saved.text
        assert (
            saved.json()["draft"]["workflow_samples"][0]["reads"][0][
                "input_file_revision_id"
            ]
            == chosen["input_file_revision_id"]
        )
        validated = client.post(
            f"/api/v1/agent/authoring-sessions/{sid}/validate",
            json={"expected_revision": saved.json()["draft_revision"]},
        )
        assert validated.status_code == 200, validated.text
        body = validated.json()
        assert body["ok"] is True, body
        binding = app.state.persistence.repository.get_validated_input_use_binding(
            body["snapshot"]["snapshot_id"]
        )
        members = [member for use in binding.input_uses for member in use.members]
        assert [member.input_file_revision_id for member in members] == [
            chosen["input_file_revision_id"]
        ]
        assert [member.content_sha256 for member in members] == [
            chosen["content_sha256"]
        ]
        assert str(root) not in validated.text
        assert "must-not-leak" not in validated.text


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
