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


def _session_factory(app):
    from encode_pipeline.persistence import create_session_factory

    return create_session_factory(app.state.persistence.engine)


def _registry(app) -> DataRegistryService:
    return DataRegistryService(
        repository=SqlAlchemyDataRegistryRepository(_session_factory(app))
    )
