"""Generate a pending target proposal, then preview and adopt it.

The model is a synthetic loopback on the ledger. No provider key is used.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest
from sqlalchemy import func, select
from sqlalchemy.orm import Session

pytest.importorskip("fastapi")

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from api_test_client import ApiTestClient  # noqa: E402
from conftest import seed_test_authentication  # noqa: E402
from encode_pipeline.api.main import create_app  # noqa: E402
from encode_pipeline.platform.experiment_catalog import CatalogError  # noqa: E402
from encode_pipeline.persistence.migrations import downgrade_database  # noqa: E402
from encode_pipeline.persistence.models import (  # noqa: E402
    CatalogProposalRow,
    ExperimentRevisionRow,
)
from encode_pipeline.services.catalog_model_client import LedgerCatalogModel  # noqa: E402
from encode_pipeline.services.data_registry import DataRegistryService  # noqa: E402
from encode_pipeline.services.peri_runtime.bound_model import (  # noqa: E402
    BoundModel,
    Budget,
    binding_id,
)
from encode_pipeline.services.peri_runtime.catalog_classify import (  # noqa: E402
    loopback_decision,
)


class Worker:
    def __init__(self) -> None:
        self.binding = {"model": "bound-a", "protocol": "anthropic"}
        self.sends: list[dict] = []

    def __call__(self, request: dict) -> dict:
        if request.get("op") == "binding":
            return dict(self.binding)
        self.sends.append(request)
        decision = loopback_decision(request["clue"], request["candidates"])
        return {"binding": dict(self.binding), "decision": decision}


@pytest.fixture
def seeded_app(tmp_path):
    application = create_app(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        workspace_root=tmp_path / "workspaces",
    )
    seed_test_authentication(application)
    yield application
    application.state.persistence.close()


def _post(client: ApiTestClient, path: str, body: dict, key: str):
    return client.post(path, json=body, headers={"Idempotency-Key": key})


def _known(term: dict) -> dict:
    return {
        "status": "known",
        "term_id": term["term_id"],
        "term_revision_id": term["revision"]["term_revision_id"],
    }


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


def test_generate_preview_and_adopt_without_moving_the_revision_early(
    seeded_app, tmp_path
):
    registry = DataRegistryService(
        repository=seeded_app.state.persistence.data_registry_repository
    )
    project = registry.create_project("Classification project")
    other = registry.create_project("Other classification project")
    worker = Worker()
    model = BoundModel(tmp_path / "catalog-ledger.sqlite", worker)
    BoundModel.initialize(model.ledger)
    model.catalog_grant(
        project_id=project.project_id,
        field="experiment.target",
        binding=binding_id(worker.binding),
        expires=2_000_000_000,
        budget=Budget(4, 4096 * 4, 64, 100, 16, 100, 100, "CNY"),
    )
    service = seeded_app.state.experiment_catalog_service
    assert service._catalog_model is None
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
        role = term("replicate_role", "biological replicate", "role")
        created = registry.import_sample(
            project.project_id,
            stable_key="classified",
            display_name="classified",
            attributes={"organism": "human"},
        )
        sample = _post(
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
            "catalog-sample",
        )
        assert sample.status_code == 200, sample.text
        experiment = _post(
            client,
            f"/api/v1/catalog/projects/{project.project_id}/experiments",
            {
                "stable_key": "classified-exp",
                "display_name": "classified-exp",
                "assay": _known(assay),
                "target": {"status": "unknown"},
                "batch": {"status": "unknown"},
                "samples": [
                    {
                        "sample_revision_id": sample.json()["sample_revision_id"],
                        "replicate_role": _known(role),
                    }
                ],
            },
            "experiment-classified",
        )
        assert experiment.status_code == 200, experiment.text
        created_experiment = experiment.json()
        suggestion_path = (
            "/api/v1/catalog/projects/"
            f"{project.project_id}/experiments/"
            f"{created_experiment['experiment_id']}/classification-suggestions"
        )
        body = {
            "experiment_revision_id": created_experiment["experiment_revision_id"],
            "field": "target",
            "reason_code": "TARGET_TERM_MATCH",
            "clue": "H3K4me3",
        }
        unbound = _post(client, suggestion_path, body, "generate-unbound")
        assert unbound.status_code == 409, unbound.text
        assert unbound.json()["issues"][0]["code"] == "catalog_model_unbound"
        assert worker.sends == []

        service._catalog_model = LedgerCatalogModel(model)
        ordinal = _post(
            client,
            suggestion_path,
            {**body, "clue": "item-1"},
            "generate-ordinal",
        )
        assert ordinal.status_code == 422, ordinal.text
        assert ordinal.json()["issues"][0]["path"] == "clue"
        assert worker.sends == []

        impersonated = _post(
            client,
            f"/api/v1/catalog/projects/{project.project_id}/proposals",
            {
                "source": "peri_bound_model",
                "items": [{"item_id": "target", "target_kind": "experiment"}],
            },
            "offline-model-source",
        )
        assert impersonated.status_code == 400, impersonated.text
        assert impersonated.json()["issues"][0]["code"] == "API_REQUEST_INVALID"
        with pytest.raises(CatalogError) as rejected_source:
            service.register_proposal(
                actor_user_id="usr_" + "a" * 32,
                idempotency_key="offline-service",
                raw_body=b"{}",
                project_id=project.project_id,
                source="peri_bound_model",
                items=[],
                rule_version=None,
                vocabulary_version=None,
                projection_version=None,
            )
        assert rejected_source.value.status_code == 422
        assert rejected_source.value.path == "source"

        generated = _post(client, suggestion_path, body, "generate-target")
        assert generated.status_code == 200, generated.text
        generated_body = generated.json()
        assert generated_body["outcome"] == "proposal"
        assert generated_body["catalog_revision_changed"] is False
        assert generated_body["decision"] == "h3k4me3"
        assert generated_body["projection"]["candidates"] == [
            {"token": "h3k4me3", "label": "H3K4me3"}
        ]
        assert (
            "term_revision"
            not in generated.text.split('"projection"')[1].split('"proposal"')[0]
        )
        proposal = generated_body["proposal"]
        assert proposal["source"] == "peri_bound_model"
        assert proposal["proposal_state"] == "pending"
        engine = seeded_app.state.persistence.engine
        with Session(engine) as db:
            stored = db.get(CatalogProposalRow, proposal["proposal_id"])
            assert stored is not None
            assert (stored.source, stored.model_source) == ("agent", "peri_bound_model")
        with pytest.raises(
            RuntimeError, match="cannot downgrade while peri_bound_model"
        ):
            downgrade_database(str(engine.url), "20261119_39")
        assert proposal["items"][0]["allowed_fields"] == ["target"]
        assert _revision_count(seeded_app, created_experiment["experiment_id"]) == 1
        assert len(worker.sends) == 1
        assert worker.sends[0]["op"] == "classify"
        assert loopback_decision("H3K4me3", worker.sends[0]["candidates"]) == "h3k4me3"

        replayed = _post(client, suggestion_path, body, "generate-target")
        assert replayed.status_code == 200, replayed.text
        assert replayed.json()["proposal"]["proposal_id"] == proposal["proposal_id"]
        assert len(worker.sends) == 1

        conflict = _post(
            client,
            suggestion_path,
            {**body, "clue": "NotThisTarget"},
            "generate-target",
        )
        assert conflict.status_code == 409, conflict.text
        assert conflict.json()["issues"][0]["code"] == "catalog_idempotency_conflict"
        assert len(worker.sends) == 1

        preview = client.post(
            f"/api/v1/catalog/proposals/{proposal['proposal_id']}/preview",
            json={"selected_item_ids": ["target"]},
        )
        assert preview.status_code == 200, preview.text
        preview_body = preview.json()
        assert preview_body["has_conflicts"] is False
        adopted = _post(
            client,
            f"/api/v1/catalog/proposals/{proposal['proposal_id']}/adoptions",
            {
                "preview_digest": preview_body["preview_digest"],
                "selected_item_ids": ["target"],
                "expected_revisions": preview_body["expected_revisions"],
            },
            "adopt-generated",
        )
        assert adopted.status_code == 200, adopted.text
        assert _revision_count(seeded_app, created_experiment["experiment_id"]) == 2
        head = client.get(
            "/api/v1/catalog/experiments/"
            f"{created_experiment['experiment_id']}/revisions/"
            f"{adopted.json()['applied_targets'][0]['result_revision_id']}"
        )
        assert head.status_code == 200, head.text
        assert head.json()["target"]["term_id"] == target["term_id"]

        declined_sample = registry.import_sample(
            project.project_id,
            stable_key="declined",
            display_name="declined",
            attributes={"organism": "human"},
        )
        declined_catalog = _post(
            client,
            f"/api/v1/catalog/samples/{declined_sample.sample.sample_id}/catalog-revisions",
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
            "catalog-declined",
        )
        assert declined_catalog.status_code == 200, declined_catalog.text
        declined = _post(
            client,
            f"/api/v1/catalog/projects/{project.project_id}/experiments",
            {
                "stable_key": "declined-exp",
                "display_name": "declined-exp",
                "assay": _known(assay),
                "target": {"status": "unknown"},
                "batch": {"status": "unknown"},
                "samples": [
                    {
                        "sample_revision_id": declined_catalog.json()[
                            "sample_revision_id"
                        ],
                        "replicate_role": _known(role),
                    }
                ],
            },
            "experiment-declined",
        )
        assert declined.status_code == 200, declined.text
        declined_body = declined.json()
        insufficient = _post(
            client,
            "/api/v1/catalog/projects/"
            f"{project.project_id}/experiments/"
            f"{declined_body['experiment_id']}/classification-suggestions",
            {
                "experiment_revision_id": declined_body["experiment_revision_id"],
                "field": "target",
                "reason_code": "TARGET_TERM_MATCH",
                "clue": "NotThisTarget",
            },
            "generate-insufficient",
        )
        assert insufficient.status_code == 200, insufficient.text
        assert insufficient.json()["outcome"] == "insufficient_evidence"
        assert insufficient.json()["proposal"] is None
        assert insufficient.json()["catalog_revision_changed"] is False
        assert _revision_count(seeded_app, declined_body["experiment_id"]) == 1
        assert len(worker.sends) == 2

        foreign_sample = registry.import_sample(
            other.project_id,
            stable_key="foreign-classified",
            display_name="foreign-classified",
            attributes={"organism": "human"},
        )
        foreign_catalog = _post(
            client,
            f"/api/v1/catalog/samples/{foreign_sample.sample.sample_id}/catalog-revisions",
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
            "catalog-foreign",
        )
        assert foreign_catalog.status_code == 200, foreign_catalog.text
        foreign = _post(
            client,
            f"/api/v1/catalog/projects/{other.project_id}/experiments",
            {
                "stable_key": "foreign-exp",
                "display_name": "foreign-exp",
                "assay": _known(assay),
                "target": {"status": "unknown"},
                "batch": {"status": "unknown"},
                "samples": [
                    {
                        "sample_revision_id": foreign_catalog.json()[
                            "sample_revision_id"
                        ],
                        "replicate_role": _known(role),
                    }
                ],
            },
            "experiment-foreign",
        )
        assert foreign.status_code == 200, foreign.text
        foreign_body = foreign.json()
        mismatched = _post(
            client,
            "/api/v1/catalog/projects/"
            f"{other.project_id}/experiments/"
            f"{foreign_body['experiment_id']}/classification-suggestions",
            {
                "experiment_revision_id": foreign_body["experiment_revision_id"],
                "field": "target",
                "reason_code": "TARGET_TERM_MATCH",
                "clue": "H3K4me3",
            },
            "generate-foreign",
        )
        assert mismatched.status_code == 409, mismatched.text
        assert (
            mismatched.json()["issues"][0]["code"] == "catalog_model_project_mismatch"
        )
        assert len(worker.sends) == 2
