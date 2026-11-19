"""HTTP adoption of one selected catalog-proposal subset against SQLite."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
import json
import sys
from pathlib import Path
from threading import Barrier, Thread

import pytest
from sqlalchemy import func, select
from sqlalchemy.orm import Session

pytest.importorskip("fastapi")

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from api_test_client import ApiTestClient  # noqa: E402
from conftest import (  # noqa: E402
    TEST_AUTH_CSRF_TOKEN,
    TEST_AUTH_PASSWORD_HASH,
    TEST_AUTH_SESSION_TOKEN,
    TEST_AUTH_USER_ID,
    seed_test_authentication,
)
from encode_pipeline.api.main import create_app  # noqa: E402
from encode_pipeline.persistence.models import (  # noqa: E402
    CatalogAdoptionReceiptRow,
    ExperimentRevisionRow,
    SampleRevisionRow,
    UserAccountRow,
)
from encode_pipeline.platform.authentication import (  # noqa: E402
    UserAccount,
    UserRole,
    UserStatus,
)
from encode_pipeline.platform.catalog_curation_preview import (  # noqa: E402
    M2_PROJECTION_VERSION,
    M2_RULE_VERSION,
    M2_VOCABULARY_VERSION,
)
from encode_pipeline.services.authentication import (  # noqa: E402
    SessionSecrets,
    new_session_record,
)
from encode_pipeline.services.data_registry import DataRegistryService  # noqa: E402


@pytest.fixture
def app(tmp_path):
    application = create_app(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        workspace_root=tmp_path / "workspaces",
    )
    yield application
    application.state.persistence.close()


@pytest.fixture
def seeded_app(app):
    seed_test_authentication(app)
    return app


def _post(client: ApiTestClient, path: str, body: dict, key: str):
    return client.post(path, json=body, headers={"Idempotency-Key": key})


def _known(term: dict) -> dict:
    return {
        "status": "known",
        "term_id": term["term_id"],
        "term_revision_id": term["revision"]["term_revision_id"],
    }


def _pin(term: dict) -> dict:
    return {
        "term_id": term["term_id"],
        "term_revision_id": term["revision"]["term_revision_id"],
    }


def _open_session(app, user_id: str, role: UserRole, session_token: str, csrf: str):
    now = datetime.now(timezone.utc)
    app.state.authentication_repository.create_account(
        UserAccount(
            user_id=user_id,
            username=user_id,
            role=role,
            status=UserStatus.ENABLED,
            password_hash=TEST_AUTH_PASSWORD_HASH,
            created_at=now,
            updated_at=now,
            password_changed_at=now,
        )
    )
    secrets = SessionSecrets._from_generated(session_token, csrf)
    app.state.authentication_repository.create_session(
        new_session_record(
            user_id=user_id,
            secrets=secrets,
            created_at=now,
            lifetime=timedelta(hours=8),
        )
    )


def _as(app, client: ApiTestClient, session_token: str, csrf: str):
    policy = app.state.auth_cookie_policy
    app.state.test_auth_tokens = None

    def request(method: str, path: str, **kwargs):
        headers = dict(kwargs.pop("headers", {}) or {})
        headers["X-CSRF-Token"] = csrf
        cookies = {
            policy.session_cookie.name: session_token,
            policy.csrf_cookie.name: csrf,
        }
        return client.request(method, path, headers=headers, cookies=cookies, **kwargs)

    return request


def test_admin_adopts_only_the_selected_subset(seeded_app) -> None:
    registry = DataRegistryService(
        repository=seeded_app.state.persistence.data_registry_repository
    )
    project = registry.create_project("Adoption project")
    other = registry.create_project("Other project")
    service = seeded_app.state.experiment_catalog_service
    with ApiTestClient(seeded_app) as client:

        def term(term_type: str, name: str, key: str) -> dict:
            response = _post(
                client,
                "/api/v1/catalog/terms",
                {"term_type": term_type, "display_name": name, "aliases": []},
                key,
            )
            assert response.status_code == 200, response.text
            return response.json()

        assay = term("assay", "CUT&Tag", "assay")
        target = term("target", "H3K4me3", "target")
        batch = term("batch", "batch-2026-11", "batch")
        role = term("replicate_role", "biological replicate", "role")
        cell_a = term("cell_type", "Treg", "cell-a")
        cell_b = term("cell_type", "T cell", "cell-b")
        tag = term("tag", "m2-curated", "tag")

        def sample(key: str) -> dict:
            created = registry.import_sample(
                project.project_id,
                stable_key=key,
                display_name=key,
                attributes={"legacy_note": "preserve-me", "organism": "human"},
            )
            response = _post(
                client,
                f"/api/v1/catalog/samples/{created.sample.sample_id}/catalog-revisions",
                {
                    "expected_revision_number": 1,
                    "catalog": {
                        "schema_version": 1,
                        "species": {"status": "unknown"},
                        "tissue": {"status": "unknown"},
                        "cell_type": {"status": "unknown", "values": []},
                        "cell_line": {"status": "not_applicable", "values": []},
                        "treatment": {"status": "not_applicable", "values": []},
                        "tag": {"status": "unknown", "values": []},
                    },
                },
                f"catalog-{key}",
            )
            assert response.status_code == 200, response.text
            return response.json()

        curated = sample("curated")
        untouched = sample("untouched")

        def experiment(
            key: str, sample_revision_id: str, *, known_target: bool
        ) -> dict:
            response = _post(
                client,
                f"/api/v1/catalog/projects/{project.project_id}/experiments",
                {
                    "stable_key": key,
                    "display_name": key,
                    "assay": _known(assay),
                    "target": _known(target) if known_target else {"status": "unknown"},
                    "batch": {"status": "unknown"},
                    "samples": [
                        {
                            "sample_revision_id": sample_revision_id,
                            "replicate_role": _known(role),
                        }
                    ],
                },
                f"experiment-{key}",
            )
            assert response.status_code == 200, response.text
            return response.json()

        first = experiment(
            "exp-batch", curated["sample_revision_id"], known_target=False
        )
        second = experiment(
            "exp-target",
            untouched["sample_revision_id"],
            known_target=False,
        )
        foreign = registry.import_sample(
            other.project_id,
            stable_key="foreign",
            display_name="foreign",
            attributes={"organism": "human"},
        )
        foreign_sample = _post(
            client,
            f"/api/v1/catalog/samples/{foreign.sample.sample_id}/catalog-revisions",
            {
                "expected_revision_number": 1,
                "catalog": {
                    "schema_version": 1,
                    "species": {"status": "unknown"},
                    "tissue": {"status": "unknown"},
                    "cell_type": {"status": "unknown", "values": []},
                    "cell_line": {"status": "not_applicable", "values": []},
                    "treatment": {"status": "not_applicable", "values": []},
                    "tag": {"status": "unknown", "values": []},
                },
            },
            "foreign-catalog",
        )
        assert foreign_sample.status_code == 200, foreign_sample.text

        def citation(*refs: dict) -> dict:
            return {
                "rule_version": M2_RULE_VERSION,
                "vocabulary_version": M2_VOCABULARY_VERSION,
                "reason_codes": ["BATCH_TERM_MATCH"],
                "source_refs": list(refs),
            }

        batch_item = {
            "item_id": "item_1",
            "target_kind": "experiment",
            "target_id": first["experiment_id"],
            "base_revision_id": first["experiment_revision_id"],
            "base_revision_number": first["revision_number"],
            "allowed_fields": ["batch"],
            "before": {"batch": {"status": "unknown"}},
            "after": {"batch": _known(batch)},
            "citations": citation(
                {
                    "kind": "experiment_revision",
                    "ref_id": first["experiment_revision_id"],
                    "evidence_token": "exp-1-batch",
                },
                {
                    "kind": "term_revision",
                    "ref_id": batch["revision"]["term_revision_id"],
                    "term_revision_id": batch["revision"]["term_revision_id"],
                    "evidence_token": "batch-term",
                },
            ),
        }
        target_item = {
            "item_id": "item_2",
            "target_kind": "experiment",
            "target_id": second["experiment_id"],
            "base_revision_id": second["experiment_revision_id"],
            "base_revision_number": second["revision_number"],
            "allowed_fields": ["target"],
            "before": {"target": {"status": "unknown"}},
            "after": {"target": _known(target)},
            "citations": citation(
                {
                    "kind": "experiment_revision",
                    "ref_id": second["experiment_revision_id"],
                    "evidence_token": "exp-2-target",
                },
                {
                    "kind": "term_revision",
                    "ref_id": target["revision"]["term_revision_id"],
                    "term_revision_id": target["revision"]["term_revision_id"],
                    "evidence_token": "target-term",
                },
            ),
        }
        sample_item = {
            "item_id": "item_3",
            "target_kind": "sample_catalog",
            "target_id": curated["sample_id"],
            "base_revision_id": curated["sample_revision_id"],
            "base_revision_number": curated["revision_number"],
            "allowed_fields": ["cell_type", "tag"],
            "before": {
                "cell_type": {"status": "unknown", "values": []},
                "tag": {"status": "unknown", "values": []},
            },
            "after": {
                "cell_type": {
                    "status": "ambiguous",
                    "values": [],
                    "candidates": [_pin(cell_a), _pin(cell_b)],
                },
                "tag": {"status": "known", "values": [_pin(tag)]},
            },
            "citations": {
                "rule_version": M2_RULE_VERSION,
                "vocabulary_version": M2_VOCABULARY_VERSION,
                "reason_codes": ["ALIAS_AMBIGUOUS_CANDIDATES", "TAG_TERM_MATCH"],
                "source_refs": [
                    {
                        "kind": "sample_revision",
                        "ref_id": curated["sample_revision_id"],
                        "evidence_token": "sample-head",
                    },
                    {
                        "kind": "term_revision",
                        "ref_id": cell_a["revision"]["term_revision_id"],
                        "term_revision_id": cell_a["revision"]["term_revision_id"],
                        "evidence_token": "cell-a",
                        "visibility": "admin_only",
                    },
                    {
                        "kind": "term_revision",
                        "ref_id": tag["revision"]["term_revision_id"],
                        "term_revision_id": tag["revision"]["term_revision_id"],
                        "evidence_token": "tag-term",
                    },
                ],
            },
        }
        registered = _post(
            client,
            f"/api/v1/catalog/projects/{project.project_id}/proposals",
            {
                "source": "offline_suggestion",
                "items": [batch_item, target_item, sample_item],
            },
            "register-main",
        )
        assert registered.status_code == 200, registered.text
        proposal = registered.json()
        assert proposal["proposal_state"] == "pending"
        assert proposal["adoptable_by_caller"] is True
        proposal_id = proposal["proposal_id"]
        replayed_register = _post(
            client,
            f"/api/v1/catalog/projects/{project.project_id}/proposals",
            {
                "source": "offline_suggestion",
                "items": [batch_item, target_item, sample_item],
            },
            "register-main",
        )
        assert replayed_register.status_code == 200, replayed_register.text
        assert replayed_register.json()["proposal_id"] == proposal_id
        moved = _post(
            client,
            f"/api/v1/catalog/experiments/{second['experiment_id']}/revisions",
            {
                "expected_revision_number": second["revision_number"],
                "display_name": "exp-target renamed",
                "assay": _known(assay),
                "target": {"status": "unknown"},
                "batch": {"status": "unknown"},
                "samples": [
                    {
                        "sample_revision_id": untouched["sample_revision_id"],
                        "replicate_role": _known(role),
                    }
                ],
            },
            "move-second",
        )
        assert moved.status_code == 200, moved.text
        frozen = _post(
            client,
            f"/api/v1/catalog/projects/{project.project_id}/proposals",
            {
                "source": "offline_suggestion",
                "items": [batch_item, target_item, sample_item],
            },
            "register-main",
        )
        assert (
            frozen.json()["items"][1]["current_revision_number"]
            == second["revision_number"]
        )
        live = client.get(f"/api/v1/catalog/proposals/{proposal_id}")
        assert live.status_code == 200, live.text
        assert (
            live.json()["items"][1]["current_revision_number"]
            == moved.json()["revision_number"]
        )
        _post(
            client,
            f"/api/v1/catalog/experiments/{second['experiment_id']}/revisions",
            {
                "expected_revision_number": moved.json()["revision_number"],
                "display_name": "exp-target",
                "assay": _known(assay),
                "target": {"status": "unknown"},
                "batch": {"status": "unknown"},
                "samples": [
                    {
                        "sample_revision_id": untouched["sample_revision_id"],
                        "replicate_role": _known(role),
                    }
                ],
            },
            "restore-second",
        )
        # The first registration stays frozen. Adopt against the restored head.
        second_now = client.get(
            f"/api/v1/catalog/experiments/{second['experiment_id']}/revisions"
        )
        assert second_now.status_code == 200, second_now.text
        head = second_now.json()["experiment_revisions"][-1]
        target_item = {
            **target_item,
            "base_revision_id": head["experiment_revision_id"],
            "base_revision_number": head["revision_number"],
            "citations": citation(
                {
                    "kind": "experiment_revision",
                    "ref_id": head["experiment_revision_id"],
                    "evidence_token": "exp-2-target",
                },
                {
                    "kind": "term_revision",
                    "ref_id": target["revision"]["term_revision_id"],
                    "term_revision_id": target["revision"]["term_revision_id"],
                    "evidence_token": "target-term",
                },
            ),
        }
        registered = _post(
            client,
            f"/api/v1/catalog/projects/{project.project_id}/proposals",
            {
                "source": "offline_suggestion",
                "items": [batch_item, target_item, sample_item],
            },
            "register-adoptable",
        )
        assert registered.status_code == 200, registered.text
        proposal = registered.json()
        proposal_id = proposal["proposal_id"]

        member_token = "e" * 43
        member_csrf = "f" * 43
        _open_session(
            seeded_app,
            "usr_" + "e" * 32,
            UserRole.MEMBER,
            member_token,
            member_csrf,
        )
        member = _as(seeded_app, client, member_token, member_csrf)
        denied = member(
            "POST",
            f"/api/v1/catalog/proposals/{proposal_id}/preview",
            json={"selected_item_ids": ["item_1"]},
        )
        assert denied.status_code == 403
        assert denied.json()["issues"][0]["code"] == "ADMINISTRATOR_REQUIRED"
        seeded_app.state.test_auth_tokens = (
            TEST_AUTH_SESSION_TOKEN,
            TEST_AUTH_CSRF_TOKEN,
        )
        viewed = member("GET", f"/api/v1/catalog/proposals/{proposal_id}")
        assert viewed.status_code == 200, viewed.text
        sample_view = next(
            item for item in viewed.json()["items"] if item["item_id"] == "item_3"
        )
        assert "candidates" not in sample_view["after"]["cell_type"]
        assert all(
            ref["ref_id"] != cell_a["revision"]["term_revision_id"]
            for ref in sample_view["citations"]["source_refs"]
        )
        seeded_app.state.test_auth_tokens = (
            TEST_AUTH_SESSION_TOKEN,
            TEST_AUTH_CSRF_TOKEN,
        )

        preview = client.post(
            f"/api/v1/catalog/proposals/{proposal_id}/preview",
            json={"selected_item_ids": ["item_3", "item_1"]},
        )
        assert preview.status_code == 200, preview.text
        preview_body = preview.json()
        assert preview_body["selected_item_ids"] == ["item_1", "item_3"]
        assert preview_body["excluded_item_ids"] == ["item_2"]
        assert preview_body["has_conflicts"] is False
        adopt_body = {
            "preview_digest": preview_body["preview_digest"],
            "selected_item_ids": ["item_1", "item_3"],
            "expected_revisions": preview_body["expected_revisions"],
        }
        adopted = _post(
            client,
            f"/api/v1/catalog/proposals/{proposal_id}/adoptions",
            adopt_body,
            "adopt-main",
        )
        assert adopted.status_code == 200, adopted.text
        receipt = adopted.json()
        assert receipt["replayed"] is False
        assert receipt["selected_item_ids"] == ["item_1", "item_3"]
        assert receipt["excluded_item_ids"] == ["item_2"]
        assert {item["item_id"] for item in receipt["applied_targets"]} == {
            "item_1",
            "item_3",
        }
        replay = _post(
            client,
            f"/api/v1/catalog/proposals/{proposal_id}/adoptions",
            adopt_body,
            "adopt-main",
        )
        assert replay.status_code == 200, replay.text
        assert replay.json()["replayed"] is True
        assert replay.json()["receipt_id"] == receipt["receipt_id"]
        assert _revision_count(seeded_app, first["experiment_id"]) == 2
        assert _revision_count(seeded_app, second["experiment_id"]) == 3
        sample_revision_id = next(
            item["result_revision_id"]
            for item in receipt["applied_targets"]
            if item["item_id"] == "item_3"
        )
        assert _legacy_note(seeded_app, sample_revision_id) == "preserve-me"
        listed = client.get(
            f"/api/v1/catalog/projects/{project.project_id}/adoption-receipts"
        )
        assert listed.status_code == 200, listed.text
        assert [item["receipt_id"] for item in listed.json()["receipts"]] == [
            receipt["receipt_id"]
        ]

        second_admin = "c" * 43
        second_csrf = "d" * 43
        _open_session(
            seeded_app,
            "usr_" + "d" * 32,
            UserRole.ADMINISTRATOR,
            second_admin,
            second_csrf,
        )
        other_admin = _as(seeded_app, client, second_admin, second_csrf)
        seeded_app.state.test_auth_tokens = None
        overlap = other_admin(
            "POST",
            f"/api/v1/catalog/proposals/{proposal_id}/adoptions",
            json={
                "preview_digest": preview_body["preview_digest"],
                "selected_item_ids": ["item_2"],
                "expected_revisions": {
                    second["experiment_id"]: head["revision_number"]
                },
            },
            headers={"Idempotency-Key": "adopt-other"},
        )
        seeded_app.state.test_auth_tokens = (
            TEST_AUTH_SESSION_TOKEN,
            TEST_AUTH_CSRF_TOKEN,
        )
        assert overlap.status_code == 409, overlap.text
        assert overlap.json()["issues"][0]["code"] == "catalog_proposal_invalidated"
        assert _receipt_count(seeded_app, proposal_id) == 1
        assert _revision_count(seeded_app, second["experiment_id"]) == 3

        mismatch = _post(
            client,
            f"/api/v1/catalog/proposals/{proposal_id}/adoptions",
            {**adopt_body, "selected_item_ids": ["item_1"]},
            "adopt-main",
        )
        assert mismatch.status_code == 409
        assert mismatch.json()["issues"][0]["code"] == "catalog_idempotency_conflict"
        digest_registered = _post(
            client,
            f"/api/v1/catalog/projects/{project.project_id}/proposals",
            {"source": "manual", "items": [target_item]},
            "register-digest",
        )
        assert digest_registered.status_code == 200, digest_registered.text
        digest_id = digest_registered.json()["proposal_id"]
        digest_preview = client.post(
            f"/api/v1/catalog/proposals/{digest_id}/preview",
            json={"selected_item_ids": ["item_2"]},
        )
        assert digest_preview.status_code == 200, digest_preview.text
        bad_digest = _post(
            client,
            f"/api/v1/catalog/proposals/{digest_id}/adoptions",
            {
                "preview_digest": "ab" * 32,
                "selected_item_ids": ["item_2"],
                "expected_revisions": digest_preview.json()["expected_revisions"],
            },
            "adopt-digest",
        )
        assert bad_digest.status_code == 409
        assert (
            bad_digest.json()["issues"][0]["code"] == "catalog_preview_digest_mismatch"
        )
        assert _receipt_count(seeded_app, digest_id) == 0

        conflict_item = {
            **batch_item,
            "item_id": "again_batch",
        }
        quiet_item = {
            **target_item,
            "item_id": "again_target",
            "target_id": second["experiment_id"],
        }
        conflict_registered = _post(
            client,
            f"/api/v1/catalog/projects/{project.project_id}/proposals",
            {"source": "manual", "items": [conflict_item, quiet_item]},
            "register-conflict",
        )
        assert conflict_registered.status_code == 200, conflict_registered.text
        conflict_id = conflict_registered.json()["proposal_id"]
        bumped = _post(
            client,
            f"/api/v1/catalog/experiments/{first['experiment_id']}/revisions",
            {
                "expected_revision_number": first["revision_number"] + 1,
                "display_name": "exp-batch bumped",
                "assay": _known(assay),
                "target": {"status": "unknown"},
                "batch": _known(batch),
                "samples": [
                    {
                        "sample_revision_id": curated["sample_revision_id"],
                        "replicate_role": _known(role),
                    }
                ],
            },
            "bump-first",
        )
        assert bumped.status_code == 200, bumped.text
        conflict_preview = client.post(
            f"/api/v1/catalog/proposals/{conflict_id}/preview",
            json={"selected_item_ids": ["again_batch", "again_target"]},
        )
        assert conflict_preview.status_code == 200, conflict_preview.text
        assert conflict_preview.json()["has_conflicts"] is True
        conflict_adopt = _post(
            client,
            f"/api/v1/catalog/proposals/{conflict_id}/adoptions",
            {
                "preview_digest": conflict_preview.json()["preview_digest"],
                "selected_item_ids": ["again_batch", "again_target"],
                "expected_revisions": conflict_preview.json()["expected_revisions"],
            },
            "adopt-conflict",
        )
        assert conflict_adopt.status_code == 409, conflict_adopt.text
        assert conflict_adopt.json()["issues"][0]["code"] == "catalog_revision_conflict"
        detail = conflict_adopt.json()["conflicts"][0]
        assert detail["target_id"] == first["experiment_id"]
        assert detail["expected_revision_id"] == first["experiment_revision_id"]
        assert detail["current_revision_id"] == bumped.json()["experiment_revision_id"]
        assert _receipt_count(seeded_app, conflict_id) == 0
        assert _revision_count(seeded_app, second["experiment_id"]) == 3

        valid_sibling = {
            **batch_item,
            "item_id": "valid_sibling",
            "base_revision_id": bumped.json()["experiment_revision_id"],
            "base_revision_number": bumped.json()["revision_number"],
            "before": {"batch": _known(batch)},
            "after": {"batch": _known(batch)},
            "citations": citation(
                {
                    "kind": "experiment_revision",
                    "ref_id": bumped.json()["experiment_revision_id"],
                    "evidence_token": "bumped-head",
                }
            ),
        }
        invalid_item = {
            **target_item,
            "item_id": "bad_source",
            "citations": citation(
                {
                    "kind": "experiment_revision",
                    "ref_id": "expr_" + "9" * 32,
                    "evidence_token": "missing-revision",
                }
            ),
        }
        invalid = _post(
            client,
            f"/api/v1/catalog/projects/{project.project_id}/proposals",
            {"source": "agent", "items": [invalid_item, valid_sibling]},
            "register-invalid",
        )
        assert invalid.status_code == 200, invalid.text
        invalid_id = invalid.json()["proposal_id"]
        invalid_preview = client.post(
            f"/api/v1/catalog/proposals/{invalid_id}/preview",
            json={"selected_item_ids": ["bad_source", "valid_sibling"]},
        )
        assert invalid_preview.status_code == 409
        assert (
            invalid_preview.json()["issues"][0]["code"]
            == "catalog_proposal_invalidated"
        )
        before_revisions = _revision_count(seeded_app, second["experiment_id"])
        before_first = _revision_count(seeded_app, first["experiment_id"])
        from encode_pipeline.platform.catalog_curation_preview import (
            compute_preview_digest,
        )

        invalid_expected = {
            second["experiment_id"]: head["revision_number"],
            first["experiment_id"]: bumped.json()["revision_number"],
        }
        invalid_adopt = _post(
            client,
            f"/api/v1/catalog/proposals/{invalid_id}/adoptions",
            {
                "preview_digest": compute_preview_digest(
                    proposal_id=invalid_id,
                    proposal_digest=invalid.json()["proposal_digest"],
                    selected_item_ids=["bad_source", "valid_sibling"],
                    expected_revisions=invalid_expected,
                ),
                "selected_item_ids": ["bad_source", "valid_sibling"],
                "expected_revisions": invalid_expected,
            },
            "adopt-invalid",
        )
        assert invalid_adopt.status_code == 409
        assert (
            invalid_adopt.json()["issues"][0]["code"] == "catalog_proposal_invalidated"
        )
        assert _revision_count(seeded_app, second["experiment_id"]) == before_revisions
        assert _revision_count(seeded_app, first["experiment_id"]) == before_first
        assert _receipt_count(seeded_app, invalid_id) == 0

        retired = _post(
            client,
            f"/api/v1/catalog/projects/{project.project_id}/proposals",
            {
                "source": "manual",
                "rule_version": "retired-rules/1",
                "vocabulary_version": M2_VOCABULARY_VERSION,
                "projection_version": M2_PROJECTION_VERSION,
                "items": [quiet_item],
            },
            "register-retired",
        )
        assert retired.status_code == 200, retired.text
        retired_preview = client.post(
            f"/api/v1/catalog/proposals/{retired.json()['proposal_id']}/preview",
            json={"selected_item_ids": ["again_target"]},
        )
        assert retired_preview.status_code == 409
        assert (
            retired_preview.json()["issues"][0]["code"]
            == "catalog_proposal_invalidated"
        )

        expiring = _post(
            client,
            f"/api/v1/catalog/projects/{project.project_id}/proposals",
            {"source": "manual", "items": [quiet_item]},
            "register-expiring",
        )
        assert expiring.status_code == 200, expiring.text
        created_at = datetime.fromisoformat(
            expiring.json()["created_at"].replace("Z", "+00:00")
        )
        service._now_factory = lambda: created_at + timedelta(minutes=31)
        expired = client.post(
            f"/api/v1/catalog/proposals/{expiring.json()['proposal_id']}/preview",
            json={"selected_item_ids": ["again_target"]},
        )
        assert expired.status_code == 410
        assert expired.json()["issues"][0]["code"] == "catalog_proposal_expired"
        service._now_factory = lambda: created_at + timedelta(days=8)
        retained = client.get(
            f"/api/v1/catalog/proposals/{expiring.json()['proposal_id']}"
        )
        assert retained.status_code == 404
        service._now_factory = lambda: datetime.now(timezone.utc)

        engine = seeded_app.state.persistence.engine
        with Session(engine) as db:
            account = db.get(UserAccountRow, TEST_AUTH_USER_ID)
            assert account is not None
            account.role = "member"
            db.commit()
        demoted = _post(
            client,
            f"/api/v1/catalog/proposals/{invalid_id}/adoptions",
            {
                "preview_digest": "cd" * 32,
                "selected_item_ids": ["again_target"],
                "expected_revisions": {
                    second["experiment_id"]: head["revision_number"]
                },
            },
            "adopt-demoted",
        )
        assert demoted.status_code == 403
        assert _receipt_count(seeded_app, invalid_id) == 0
        with Session(engine) as db:
            account = db.get(UserAccountRow, TEST_AUTH_USER_ID)
            assert account is not None
            account.role = "administrator"
            db.commit()

        race_first = experiment(
            "race-a", curated["sample_revision_id"], known_target=True
        )
        race_second = experiment(
            "race-b", untouched["sample_revision_id"], known_target=True
        )
        race_items = [
            {
                **batch_item,
                "item_id": "race_1",
                "target_id": race_first["experiment_id"],
                "base_revision_id": race_first["experiment_revision_id"],
                "base_revision_number": race_first["revision_number"],
                "allowed_fields": ["display_name"],
                "before": {"display_name": "race-a"},
                "after": {"display_name": "race-a adopted"},
                "citations": citation(
                    {
                        "kind": "experiment_revision",
                        "ref_id": race_first["experiment_revision_id"],
                        "evidence_token": "race-a",
                    }
                ),
            },
            {
                **batch_item,
                "item_id": "race_2",
                "target_id": race_second["experiment_id"],
                "base_revision_id": race_second["experiment_revision_id"],
                "base_revision_number": race_second["revision_number"],
                "allowed_fields": ["display_name"],
                "before": {"display_name": "race-b"},
                "after": {"display_name": "race-b adopted"},
                "citations": citation(
                    {
                        "kind": "experiment_revision",
                        "ref_id": race_second["experiment_revision_id"],
                        "evidence_token": "race-b",
                    }
                ),
            },
        ]
        race_registered = _post(
            client,
            f"/api/v1/catalog/projects/{project.project_id}/proposals",
            {"source": "offline_suggestion", "items": race_items},
            "register-race",
        )
        assert race_registered.status_code == 200, race_registered.text
        race_id = race_registered.json()["proposal_id"]
        left = client.post(
            f"/api/v1/catalog/proposals/{race_id}/preview",
            json={"selected_item_ids": ["race_1"]},
        )
        right = client.post(
            f"/api/v1/catalog/proposals/{race_id}/preview",
            json={"selected_item_ids": ["race_2"]},
        )
        assert left.status_code == 200, left.text
        assert right.status_code == 200, right.text
        barrier = Barrier(2)
        statuses: list[int] = []

        def adopt(token: str, csrf: str, body: dict, key: str) -> None:
            barrier.wait()
            response = _as(seeded_app, client, token, csrf)(
                "POST",
                f"/api/v1/catalog/proposals/{race_id}/adoptions",
                json=body,
                headers={"Idempotency-Key": key},
            )
            statuses.append(response.status_code)

        seeded_app.state.test_auth_tokens = None
        threads = [
            Thread(
                target=adopt,
                args=(
                    TEST_AUTH_SESSION_TOKEN,
                    TEST_AUTH_CSRF_TOKEN,
                    {
                        "preview_digest": left.json()["preview_digest"],
                        "selected_item_ids": ["race_1"],
                        "expected_revisions": left.json()["expected_revisions"],
                    },
                    "race-left",
                ),
            ),
            Thread(
                target=adopt,
                args=(
                    second_admin,
                    second_csrf,
                    {
                        "preview_digest": right.json()["preview_digest"],
                        "selected_item_ids": ["race_2"],
                        "expected_revisions": right.json()["expected_revisions"],
                    },
                    "race-right",
                ),
            ),
        ]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        seeded_app.state.test_auth_tokens = (
            TEST_AUTH_SESSION_TOKEN,
            TEST_AUTH_CSRF_TOKEN,
        )
        assert sorted(statuses) == [200, 409]
        assert _receipt_count(seeded_app, race_id) == 1
        revised = {
            _revision_count(seeded_app, race_first["experiment_id"]),
            _revision_count(seeded_app, race_second["experiment_id"]),
        }
        assert revised == {1, 2}

        held = experiment("exp-held", curated["sample_revision_id"], known_target=False)

        def held_item(item_id: str) -> dict:
            return {
                "item_id": item_id,
                "target_kind": "experiment",
                "target_id": held["experiment_id"],
                "base_revision_id": held["experiment_revision_id"],
                "base_revision_number": held["revision_number"],
                "allowed_fields": ["display_name"],
                "before": {"display_name": "exp-held"},
                "after": {"display_name": f"exp-held {item_id}"},
                "citations": citation(
                    {
                        "kind": "experiment_revision",
                        "ref_id": held["experiment_revision_id"],
                        "evidence_token": item_id,
                    }
                ),
            }

        def register_held(item_id: str, key: str) -> str:
            response = _post(
                client,
                f"/api/v1/catalog/projects/{project.project_id}/proposals",
                {"source": "offline_suggestion", "items": [held_item(item_id)]},
                key,
            )
            assert response.status_code == 200, response.text
            return response.json()["proposal_id"]

        left_proposal = register_held("held_left", "register-held-left")
        right_proposal = register_held("held_right", "register-held-right")
        left_preview = client.post(
            f"/api/v1/catalog/proposals/{left_proposal}/preview",
            json={"selected_item_ids": ["held_left"]},
        )
        right_preview = other_admin(
            "POST",
            f"/api/v1/catalog/proposals/{right_proposal}/preview",
            json={"selected_item_ids": ["held_right"]},
        )
        assert left_preview.status_code == 200, left_preview.text
        assert right_preview.status_code == 200, right_preview.text
        collision: list[tuple[int, str]] = []

        def collide(request, proposal_id: str, preview, key: str) -> None:
            barrier = collision_barrier
            barrier.wait()
            response = request(
                "POST",
                f"/api/v1/catalog/proposals/{proposal_id}/adoptions",
                json={
                    "preview_digest": preview["preview_digest"],
                    "selected_item_ids": preview["selected_item_ids"],
                    "expected_revisions": preview["expected_revisions"],
                },
                headers={"Idempotency-Key": key},
            )
            code = ""
            if response.status_code != 200:
                code = response.json()["issues"][0]["code"]
            collision.append((response.status_code, code))

        collision_barrier = Barrier(2)
        primary = _as(
            seeded_app,
            client,
            TEST_AUTH_SESSION_TOKEN,
            TEST_AUTH_CSRF_TOKEN,
        )
        threads = [
            Thread(
                target=collide,
                args=(
                    primary,
                    left_proposal,
                    left_preview.json(),
                    "held-left",
                ),
            ),
            Thread(
                target=collide,
                args=(
                    other_admin,
                    right_proposal,
                    right_preview.json(),
                    "held-right",
                ),
            ),
        ]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        seeded_app.state.test_auth_tokens = (
            TEST_AUTH_SESSION_TOKEN,
            TEST_AUTH_CSRF_TOKEN,
        )
        assert sorted(collision) == [
            (200, ""),
            (409, "catalog_revision_conflict"),
        ]
        assert (
            _receipt_count(seeded_app, left_proposal)
            + _receipt_count(seeded_app, right_proposal)
            == 1
        )
        assert _revision_count(seeded_app, held["experiment_id"]) == 2
        original = client.get(
            "/api/v1/catalog/experiments/"
            f"{held['experiment_id']}/revisions/{held['experiment_revision_id']}"
        )
        assert original.status_code == 200, original.text
        assert original.json()["revision_number"] == held["revision_number"]
        receipts = client.get(
            f"/api/v1/catalog/projects/{project.project_id}/adoption-receipts"
        )
        assert receipts.status_code == 200, receipts.text
        winning = [
            item
            for item in receipts.json()["receipts"]
            if item["proposal_id"] in {left_proposal, right_proposal}
        ]
        assert len(winning) == 1
        assert winning[0]["applied_targets"][0]["target_id"] == held["experiment_id"]

        fresh = experiment(
            "exp-retry", curated["sample_revision_id"], known_target=False
        )
        retry_item = {
            **held_item("retry_once"),
            "target_id": fresh["experiment_id"],
            "base_revision_id": fresh["experiment_revision_id"],
            "base_revision_number": fresh["revision_number"],
            "before": {"display_name": "exp-retry"},
            "after": {"display_name": "exp-retry adopted"},
            "citations": citation(
                {
                    "kind": "experiment_revision",
                    "ref_id": fresh["experiment_revision_id"],
                    "evidence_token": "retry-once",
                }
            ),
        }
        retry_registered = _post(
            client,
            f"/api/v1/catalog/projects/{project.project_id}/proposals",
            {"source": "offline_suggestion", "items": [retry_item]},
            "register-retry-fresh",
        )
        assert retry_registered.status_code == 200, retry_registered.text
        retry_id = retry_registered.json()["proposal_id"]
        retry_preview = client.post(
            f"/api/v1/catalog/proposals/{retry_id}/preview",
            json={"selected_item_ids": ["retry_once"]},
        )
        assert retry_preview.status_code == 200, retry_preview.text
        retry_body = {
            "preview_digest": retry_preview.json()["preview_digest"],
            "selected_item_ids": ["retry_once"],
            "expected_revisions": retry_preview.json()["expected_revisions"],
        }
        retry_barrier = Barrier(2)
        retried: list[dict] = []

        def retry_once() -> None:
            retry_barrier.wait()
            response = primary(
                "POST",
                f"/api/v1/catalog/proposals/{retry_id}/adoptions",
                json=retry_body,
                headers={"Idempotency-Key": "retry-same-request"},
            )
            assert response.status_code == 200, response.text
            retried.append(response.json())

        threads = [Thread(target=retry_once), Thread(target=retry_once)]
        for thread in threads:
            thread.start()
        for thread in threads:
            thread.join()
        seeded_app.state.test_auth_tokens = (
            TEST_AUTH_SESSION_TOKEN,
            TEST_AUTH_CSRF_TOKEN,
        )
        assert {item["receipt_id"] for item in retried} == {retried[0]["receipt_id"]}
        assert {item["replayed"] for item in retried} == {False, True}
        assert _receipt_count(seeded_app, retry_id) == 1
        assert _revision_count(seeded_app, fresh["experiment_id"]) == 2
        readable = client.get(
            f"/api/v1/catalog/projects/{project.project_id}/adoption-receipts"
        )
        assert readable.status_code == 200, readable.text
        stored = [
            item
            for item in readable.json()["receipts"]
            if item["proposal_id"] == retry_id
        ]
        assert stored[0]["receipt_id"] == retried[0]["receipt_id"]
        assert stored[0]["replayed"] is False


def _revision_count(app, experiment_id: str) -> int:
    engine = app.state.persistence.engine
    with Session(engine) as db:
        return int(
            db.scalar(
                select(func.count())
                .select_from(ExperimentRevisionRow)
                .where(ExperimentRevisionRow.experiment_id == experiment_id)
            )
            or 0
        )


def _receipt_count(app, proposal_id: str) -> int:
    engine = app.state.persistence.engine
    with Session(engine) as db:
        return int(
            db.scalar(
                select(func.count())
                .select_from(CatalogAdoptionReceiptRow)
                .where(CatalogAdoptionReceiptRow.proposal_id == proposal_id)
            )
            or 0
        )


def _legacy_note(app, sample_revision_id: str) -> str:
    engine = app.state.persistence.engine
    with Session(engine) as db:
        row = db.get(SampleRevisionRow, sample_revision_id)
        assert row is not None
        return json.loads(row.canonical_payload)["attributes"]["legacy_note"]
