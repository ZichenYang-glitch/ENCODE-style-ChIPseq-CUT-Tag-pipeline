"""HTTP coverage for Slice A catalog maintenance, queries, and trace envelopes."""

from __future__ import annotations

import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

pytest.importorskip("fastapi")

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from api_test_client import ApiTestClient  # noqa: E402
from conftest import (  # noqa: E402
    TEST_AUTH_PASSWORD_HASH,
    TEST_AUTH_USER_ID,
    seed_test_authentication,
)
from encode_pipeline.api.main import create_app  # noqa: E402
from encode_pipeline.platform.results import Result  # noqa: E402
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


def _registry(app) -> DataRegistryService:
    return DataRegistryService(
        repository=app.state.persistence.data_registry_repository
    )


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


def _catalog(cell_type: dict) -> dict:
    return {
        "schema_version": 1,
        "species": {"status": "unknown"},
        "tissue": {"status": "unknown"},
        "cell_type": cell_type,
        "cell_line": {"status": "not_applicable", "values": []},
        "treatment": {"status": "not_applicable", "values": []},
        "tag": {"status": "not_applicable", "values": []},
    }


def test_catalog_http_queries_permissions_and_idempotency(seeded_app) -> None:
    registry = _registry(seeded_app)
    project = registry.create_project("Treg project")
    with ApiTestClient(seeded_app) as client:
        tokens = seeded_app.state.test_auth_tokens
        seeded_app.state.test_auth_tokens = None
        anonymous = client.post(
            "/api/v1/catalog/terms",
            json={"term_type": "tag", "display_name": "x"},
            headers={"Idempotency-Key": "anon"},
        )
        seeded_app.state.test_auth_tokens = tokens
        assert anonymous.status_code == 401
        assert anonymous.json()["issues"][0]["code"] == "AUTHENTICATION_REQUIRED"

        def term(term_type: str, name: str, key: str) -> dict:
            response = _post(
                client,
                "/api/v1/catalog/terms",
                {
                    "term_type": term_type,
                    "display_name": name,
                    "aliases": [name] if term_type == "cell_type" else [],
                },
                key,
            )
            assert response.status_code == 200, response.text
            return response.json()

        cuttag = term("assay", "CUT&Tag", "assay-cut")
        atac = term("assay", "ATAC-seq", "assay-atac")
        target = term("target", "H3K27ac", "target")
        batch = term("batch", "batch-1", "batch")
        role = term("replicate_role", "biological replicate", "role")
        treg = term("cell_type", "regulatory T cell", "treg")
        keratinocyte = term("cell_type", "keratinocyte", "skin")
        replay = _post(
            client,
            "/api/v1/catalog/terms",
            {
                "term_type": "cell_type",
                "display_name": "regulatory T cell",
                "aliases": ["regulatory T cell"],
            },
            "treg",
        )
        assert replay.status_code == 200
        assert (
            replay.json()["revision"]["term_revision_id"]
            == (treg["revision"]["term_revision_id"])
        )
        conflict = _post(
            client,
            "/api/v1/catalog/terms",
            {"term_type": "cell_type", "display_name": "other name"},
            "treg",
        )
        assert conflict.status_code == 409
        assert conflict.json()["issues"][0]["code"] == "catalog_idempotency_conflict"

        listed_terms = client.get("/api/v1/catalog/terms")
        assert listed_terms.status_code == 200
        assert "aliases" not in listed_terms.json()["terms"][0]
        admin_term = client.get(f"/api/v1/catalog/terms/{treg['term_id']}")
        assert admin_term.status_code == 200
        assert admin_term.json()["revision"]["aliases"] == ["regulatory T cell"]

        def classified(key: str, cell_type: dict):
            created = registry.import_sample(
                project.project_id,
                stable_key=key,
                display_name=key,
                attributes={"organism": "human"},
            )
            body = {
                "expected_revision_number": 1,
                "display_name": key,
                "catalog": _catalog(cell_type),
            }
            response = _post(
                client,
                f"/api/v1/catalog/samples/{created.sample.sample_id}/catalog-revisions",
                body,
                f"catalog-{key}",
            )
            assert response.status_code == 200, response.text
            return response.json()

        treg_sample = classified(
            "treg",
            {"status": "known", "values": [_pin(treg)]},
        )
        renamed = _post(
            client,
            f"/api/v1/catalog/terms/{treg['term_id']}/revisions",
            {
                "expected_revision_number": 1,
                "change_kind": "editorial",
                "display_name": "Treg",
                "aliases": ["regulatory T cell"],
            },
            "treg-rename",
        )
        assert renamed.status_code == 200, renamed.text
        pinned = client.get(
            "/api/v1/catalog/samples/"
            f"{treg_sample['sample_id']}/revisions/"
            f"{treg_sample['sample_revision_id']}"
        )
        assert pinned.status_code == 200
        assert pinned.json()["cell_type"]["values"][0]["display_name"] == (
            "regulatory T cell"
        )
        assert (
            client.get(f"/api/v1/catalog/terms/{treg['term_id']}").json()["revision"][
                "display_name"
            ]
            == "Treg"
        )
        skin_sample = classified(
            "skin",
            {"status": "known", "values": [_pin(keratinocyte)]},
        )
        unknown_sample = classified(
            "unknown",
            {"status": "unknown", "values": []},
        )
        plain = registry.import_sample(
            project.project_id,
            stable_key="plain",
            display_name="plain",
            attributes={"organism": "human"},
        )
        invalid = _post(
            client,
            f"/api/v1/catalog/samples/{plain.sample.sample_id}/catalog-revisions",
            {
                "expected_revision_number": 1,
                "catalog": {"schema_version": 1, "species": {"status": "unknown"}},
            },
            "invalid-catalog",
        )
        assert invalid.status_code == 422
        assert invalid.json()["issues"][0]["code"] == "catalog_payload_invalid"
        retried = _post(
            client,
            f"/api/v1/catalog/samples/{plain.sample.sample_id}/catalog-revisions",
            {
                "expected_revision_number": 1,
                "catalog": None,
            },
            "invalid-catalog",
        )
        assert retried.status_code == 200
        assert retried.json()["catalog_state"] == "absent"

        def experiment(key: str, assay: dict, sample_revision_id: str) -> dict:
            response = _post(
                client,
                f"/api/v1/catalog/projects/{project.project_id}/experiments",
                {
                    "stable_key": key,
                    "display_name": key,
                    "assay": assay,
                    "target": _known(target),
                    "batch": _known(batch),
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

        cut_treg = experiment(
            "cut-treg", _known(cuttag), treg_sample["sample_revision_id"]
        )
        experiment("atac-skin", _known(atac), skin_sample["sample_revision_id"])
        unknown = experiment(
            "cut-unknown",
            _known(cuttag),
            unknown_sample["sample_revision_id"],
        )
        matched = client.get(
            f"/api/v1/catalog/projects/{project.project_id}/experiments",
            params=[
                ("assay", cuttag["term_id"]),
                ("cell_type", treg["term_id"]),
            ],
        )
        assert matched.status_code == 200
        assert [item["experiment_id"] for item in matched.json()["experiments"]] == [
            cut_treg["experiment_id"]
        ]
        assert "candidates" not in matched.json()["experiments"][0]["assay"]
        ambiguous = experiment(
            "ambiguous-assay",
            {
                "status": "ambiguous",
                "candidates": [_pin(cuttag), _pin(atac)],
            },
            skin_sample["sample_revision_id"],
        )
        admin_experiment = client.get(
            f"/api/v1/catalog/experiments/{ambiguous['experiment_id']}"
            f"/revisions/{ambiguous['experiment_revision_id']}"
        )
        assert admin_experiment.status_code == 200
        assert len(admin_experiment.json()["assay"]["candidates"]) == 2
        ambiguous_list = client.get(
            f"/api/v1/catalog/projects/{project.project_id}/experiments",
            params={"assay_status": "ambiguous"},
        )
        assert [
            item["experiment_id"] for item in ambiguous_list.json()["experiments"]
        ] == [ambiguous["experiment_id"]]
        assert "candidates" not in ambiguous_list.json()["experiments"][0]["assay"]
        unknown_only = client.get(
            f"/api/v1/catalog/projects/{project.project_id}/experiments",
            params={"cell_type_status": "unknown"},
        )
        assert [
            item["experiment_id"] for item in unknown_only.json()["experiments"]
        ] == [unknown["experiment_id"]]
        absent = client.get(
            f"/api/v1/catalog/projects/{project.project_id}/samples",
            params={"catalog_state": "absent"},
        )
        assert [item["sample_id"] for item in absent.json()["samples"]] == [
            plain.sample.sample_id
        ]
        trace = client.get(
            "/api/v1/catalog/projects/"
            f"{project.project_id}/experiment-revisions/"
            f"{cut_treg['experiment_revision_id']}/run-results"
        )
        assert trace.status_code == 200
        assert trace.json()["trace_status"] == "unresolved"
        assert trace.json()["reason"] == "experiment_snapshot_binding_absent"
        assert trace.json()["results"] == []

        stale = _post(
            client,
            f"/api/v1/catalog/experiments/{cut_treg['experiment_id']}/revisions",
            {
                "expected_revision_number": 1,
                "display_name": "cut-treg revised",
                "assay": _known(cuttag),
                "target": _known(target),
                "batch": _known(batch),
                "samples": [
                    {
                        "sample_revision_id": treg_sample["sample_revision_id"],
                        "replicate_role": _known(role),
                    }
                ],
            },
            "revise-a",
        )
        assert stale.status_code == 200
        loser = _post(
            client,
            f"/api/v1/catalog/experiments/{cut_treg['experiment_id']}/revisions",
            {
                "expected_revision_number": 1,
                "display_name": "cut-treg other",
                "assay": _known(cuttag),
                "target": _known(target),
                "batch": _known(batch),
                "samples": [
                    {
                        "sample_revision_id": skin_sample["sample_revision_id"],
                        "replicate_role": _known(role),
                    }
                ],
            },
            "revise-b",
        )
        assert loser.status_code == 409
        assert loser.json()["issues"][0]["code"] == "catalog_revision_conflict"

        from datetime import datetime, timezone

        from sqlalchemy.orm import Session

        from encode_pipeline.persistence.models import UserAccountRow
        from encode_pipeline.platform.authentication import (
            UserAccount,
            UserRole,
            UserStatus,
        )

        # Keep one enabled administrator so setup stays complete, then demote
        # the seeded session. Member reads and writes are rejected.
        now = datetime.now(timezone.utc)
        seeded_app.state.authentication_repository.create_account(
            UserAccount(
                user_id="usr_" + "d" * 32,
                username="other-admin",
                role=UserRole.ADMINISTRATOR,
                status=UserStatus.ENABLED,
                password_hash=TEST_AUTH_PASSWORD_HASH,
                created_at=now,
                updated_at=now,
                password_changed_at=now,
            )
        )
        engine = seeded_app.state.persistence.engine
        with Session(engine) as db:
            account = db.get(UserAccountRow, TEST_AUTH_USER_ID)
            assert account is not None
            account.role = "member"
            db.commit()
        denied = client.get(f"/api/v1/catalog/terms/{treg['term_id']}")
        assert denied.status_code == 403
        assert denied.json()["issues"][0]["code"] == "ADMINISTRATOR_REQUIRED"
        denied_write = _post(
            client,
            "/api/v1/catalog/terms",
            {"term_type": "tag", "display_name": "member-tag"},
            "member-write",
        )
        assert denied_write.status_code == 403
        assert denied_write.json()["issues"][0]["code"] == "ADMINISTRATOR_REQUIRED"
        with Session(engine) as db:
            account = db.get(UserAccountRow, TEST_AUTH_USER_ID)
            assert account is not None
            account.role = "administrator"
            db.commit()


def test_catalog_read_does_not_change_download_authorization(seeded_app) -> None:
    plan = SimpleNamespace(
        iter_bytes=lambda: iter((b"artifact-bytes",)),
        close=lambda: None,
        media_type="application/octet-stream",
        content_disposition='attachment; filename="artifact.bin"',
        size_bytes=15,
    )

    class StubDownloadService:
        def prepare(self, run_id, artifact_id, **kwargs):
            del run_id, artifact_id, kwargs
            return Result.success(plan)

    seeded_app.state.artifact_download_service = StubDownloadService()
    path = (
        "/api/v1/runs/run-1/artifacts/artifact-1/download"
        "?generation=artifactgen-" + "a" * 64 + "&revision=artifactrev-" + "b" * 64
    )
    with ApiTestClient(seeded_app) as client:
        before = client.get(path)
        catalog = client.get("/api/v1/catalog/terms")
        after = client.get(path)
    assert catalog.status_code == 200
    assert before.status_code == after.status_code == 200
    assert before.content == after.content == b"artifact-bytes"
