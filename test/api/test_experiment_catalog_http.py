"""HTTP coverage for Slice A catalog maintenance, queries, and trace envelopes."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest
from pydantic import TypeAdapter

pytest.importorskip("fastapi")

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from api_test_client import ApiTestClient  # noqa: E402
from conftest import (  # noqa: E402
    TEST_AUTH_PASSWORD_HASH,
    TEST_AUTH_USER_ID,
    seed_test_authentication,
)
from encode_pipeline.api.catalog_models import (  # noqa: E402
    CatalogExperimentListResponse,
    CatalogExperimentResponse,
    CatalogExperimentRunListResponse,
    CatalogSampleListResponse,
    CatalogSampleRevisionResponse,
    CatalogSampleRunListResponse,
)
from encode_pipeline.api.main import create_app  # noqa: E402
from encode_pipeline.platform.results import Result  # noqa: E402
from encode_pipeline.services.data_registry import DataRegistryService  # noqa: E402


def _generated_catalog_urls(calls: list[dict[str, object]]) -> list[str]:
    """Call the committed Orval URL builders and return the paths they emit."""

    repository = Path(__file__).resolve().parents[2]
    node = repository / ".local/node-v22.23.1-linux-x64/bin/node"
    if not node.is_file():
        node = Path("node")
    completed = subprocess.run(
        [
            str(node),
            str(repository / "test/api/generated_catalog_url.mjs"),
            json.dumps(calls),
        ],
        cwd=repository,
        capture_output=True,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        detail = completed.stderr.strip() or completed.stdout.strip()
        raise AssertionError(
            detail or f"generated client exited {completed.returncode}"
        )
    urls = json.loads(completed.stdout)
    if not isinstance(urls, list) or len(urls) != len(calls):
        raise AssertionError(completed.stdout)
    return urls


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
        missing_key = client.post(
            "/api/v1/catalog/terms",
            json={"term_type": "tag", "display_name": "x"},
        )
        assert missing_key.status_code == 400
        assert missing_key.json()["issues"][0]["code"] == "API_REQUEST_INVALID"
        schema = seeded_app.openapi()
        write_operations = {
            ("/api/v1/catalog/terms", "post"),
            ("/api/v1/catalog/terms/{term_id}/revisions", "post"),
            ("/api/v1/catalog/terms/{term_id}/archive", "post"),
            ("/api/v1/catalog/projects/{project_id}/experiments", "post"),
            ("/api/v1/catalog/experiments/{experiment_id}/revisions", "post"),
            ("/api/v1/catalog/experiments/{experiment_id}/archive", "post"),
            ("/api/v1/catalog/samples/{sample_id}/catalog-revisions", "post"),
        }
        for path, method in write_operations:
            parameters = schema["paths"][path][method]["parameters"]
            header = next(
                item
                for item in parameters
                if item["name"] == "Idempotency-Key" and item["in"] == "header"
            )
            assert header["required"] is True
        experiment_query = {
            item["name"]
            for item in schema["paths"][
                "/api/v1/catalog/projects/{project_id}/experiments"
            ]["get"]["parameters"]
        }
        assert {"cell_type_term_id", "assay_term_id", "cell_type_status"} <= (
            experiment_query
        )
        assert {"cell_type", "assay"}.isdisjoint(experiment_query)
        for path, name in (
            (
                "/api/v1/catalog/projects/{project_id}/experiments",
                "cell_type_term_id",
            ),
            (
                "/api/v1/catalog/projects/{project_id}/experiments",
                "assay_term_id",
            ),
            (
                "/api/v1/catalog/projects/{project_id}/samples",
                "cell_type_term_id",
            ),
        ):
            parameter = next(
                item
                for item in schema["paths"][path]["get"]["parameters"]
                if item["name"] == name
            )
            assert parameter["schema"]["type"] == "array"
            assert parameter["schema"]["items"] == {"type": "string"}

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
        atac_skin = experiment(
            "atac-skin", _known(atac), skin_sample["sample_revision_id"]
        )
        cut_skin = experiment(
            "cut-skin", _known(cuttag), skin_sample["sample_revision_id"]
        )
        unknown = experiment(
            "cut-unknown",
            _known(cuttag),
            unknown_sample["sample_revision_id"],
        )

        def experiment_ids(response) -> list[str]:
            assert response.status_code == 200, response.text
            return [item["experiment_id"] for item in response.json()["experiments"]]

        matched = client.get(
            f"/api/v1/catalog/projects/{project.project_id}/experiments",
            params=[
                ("assay_term_id", cuttag["term_id"]),
                ("cell_type_term_id", treg["term_id"]),
            ],
        )
        assert matched.json()["ok"] is True
        CatalogExperimentListResponse.model_validate(matched.json())
        CatalogExperimentResponse.model_validate(cut_treg)
        assert experiment_ids(matched) == [cut_treg["experiment_id"]]
        either_cell = client.get(
            f"/api/v1/catalog/projects/{project.project_id}/experiments",
            params=[
                ("cell_type_term_id", treg["term_id"]),
                ("cell_type_term_id", keratinocyte["term_id"]),
                ("assay_term_id", cuttag["term_id"]),
            ],
        )
        assert set(experiment_ids(either_cell)) == {
            cut_treg["experiment_id"],
            cut_skin["experiment_id"],
        }
        (
            both_url,
            single_url,
            empty_url,
            sample_or_url,
            sample_one_url,
            sample_empty_url,
        ) = _generated_catalog_urls(
            [
                {
                    "fn": "getCatalogListExperimentsUrl",
                    "args": [
                        project.project_id,
                        {
                            "cell_type_term_id": [
                                treg["term_id"],
                                keratinocyte["term_id"],
                            ],
                            "assay_term_id": [cuttag["term_id"]],
                        },
                    ],
                },
                {
                    "fn": "getCatalogListExperimentsUrl",
                    "args": [
                        project.project_id,
                        {
                            "cell_type_term_id": [treg["term_id"]],
                            "assay_term_id": [cuttag["term_id"]],
                        },
                    ],
                },
                {
                    "fn": "getCatalogListExperimentsUrl",
                    "args": [project.project_id, {}],
                },
                {
                    "fn": "getCatalogListSamplesUrl",
                    "args": [
                        project.project_id,
                        {
                            "cell_type_term_id": [
                                treg["term_id"],
                                keratinocyte["term_id"],
                            ],
                        },
                    ],
                },
                {
                    "fn": "getCatalogListSamplesUrl",
                    "args": [
                        project.project_id,
                        {"cell_type_term_id": [keratinocyte["term_id"]]},
                    ],
                },
                {
                    "fn": "getCatalogListSamplesUrl",
                    "args": [project.project_id, {}],
                },
            ]
        )
        assert "%2C" not in both_url
        assert both_url.count("cell_type_term_id=") == 2
        assert both_url.count("assay_term_id=") == 1
        assert set(experiment_ids(client.get(both_url))) == {
            cut_treg["experiment_id"],
            cut_skin["experiment_id"],
        }
        assert "%2C" not in single_url
        assert single_url.count("cell_type_term_id=") == 1
        assert experiment_ids(client.get(single_url)) == [cut_treg["experiment_id"]]
        assert "?" not in empty_url
        assert set(experiment_ids(client.get(empty_url))) == {
            cut_treg["experiment_id"],
            atac_skin["experiment_id"],
            cut_skin["experiment_id"],
            unknown["experiment_id"],
        }
        assert "%2C" not in sample_or_url
        assert sample_or_url.count("cell_type_term_id=") == 2
        sample_or = client.get(sample_or_url)
        assert sample_or.status_code == 200, sample_or.text
        assert {item["sample_id"] for item in sample_or.json()["samples"]} == {
            treg_sample["sample_id"],
            skin_sample["sample_id"],
        }
        sample_one = client.get(sample_one_url)
        assert {item["sample_id"] for item in sample_one.json()["samples"]} == {
            skin_sample["sample_id"],
        }
        sample_empty = client.get(sample_empty_url)
        assert "?" not in sample_empty_url
        assert {item["sample_id"] for item in sample_empty.json()["samples"]} == {
            treg_sample["sample_id"],
            skin_sample["sample_id"],
            unknown_sample["sample_id"],
            plain.sample.sample_id,
        }
        unknown_cut = client.get(
            f"/api/v1/catalog/projects/{project.project_id}/experiments",
            params=[
                ("cell_type_status", "unknown"),
                ("assay_term_id", cuttag["term_id"]),
            ],
        )
        assert experiment_ids(unknown_cut) == [unknown["experiment_id"]]
        mixed_dimension = client.get(
            f"/api/v1/catalog/projects/{project.project_id}/experiments",
            params=[
                ("cell_type_status", "unknown"),
                ("cell_type_term_id", treg["term_id"]),
            ],
        )
        assert mixed_dimension.status_code == 400
        assert mixed_dimension.json()["issues"][0]["code"] == "catalog_filter_invalid"
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
        CatalogSampleListResponse.model_validate(absent.json())
        sample_revision = TypeAdapter(CatalogSampleRevisionResponse)
        sample_revision.validate_python(pinned.json())
        sample_revision.validate_python(treg_sample)
        sample_trace = client.get(
            "/api/v1/catalog/projects/"
            f"{project.project_id}/sample-revisions/"
            f"{treg_sample['sample_revision_id']}/run-results"
        )
        assert sample_trace.status_code == 200
        assert sample_trace.json()["ok"] is True
        CatalogSampleRunListResponse.model_validate(sample_trace.json())
        trace = client.get(
            "/api/v1/catalog/projects/"
            f"{project.project_id}/experiment-revisions/"
            f"{cut_treg['experiment_revision_id']}/run-results"
        )
        assert trace.status_code == 200
        assert trace.json()["ok"] is True
        CatalogExperimentRunListResponse.model_validate(trace.json())
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
