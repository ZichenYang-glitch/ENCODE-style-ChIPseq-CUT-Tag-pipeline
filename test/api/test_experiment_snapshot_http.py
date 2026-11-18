"""HTTP evidence for freezing an experiment revision onto a validation snapshot."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timezone
from hashlib import sha256

import pytest
from sqlalchemy import text

pytest.importorskip("fastapi")

from api_test_client import ApiTestClient
from encode_pipeline.api.catalog_models import CatalogExperimentRunListResponse
from encode_pipeline.platform.runs import RunArtifactRef, RunStatus
from encode_pipeline.services.data_registry import DataRegistryService
from encode_pipeline.services.run_repositories import RunEventDraft


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


def _post(client: ApiTestClient, path: str, body: dict, key: str):
    return client.post(path, json=body, headers={"Idempotency-Key": key})


def _known(term: dict) -> dict:
    return {
        "status": "known",
        "term_id": term["term_id"],
        "term_revision_id": term["revision"]["term_revision_id"],
    }


def _count(app, table: str) -> int:
    with app.state.persistence.engine.connect() as connection:
        return int(
            connection.execute(text(f"SELECT count(*) FROM {table}")).scalar_one()
        )


def _binding(app, snapshot_id: str):
    with app.state.persistence.engine.connect() as connection:
        return (
            connection.execute(
                text(
                    "SELECT experiment_revision_id, payload_digest "
                    "FROM snapshot_experiment_bindings WHERE snapshot_id=:snapshot_id"
                ),
                {"snapshot_id": snapshot_id},
            )
            .mappings()
            .one()
        )


def _sample_row(tmp_path, name: str) -> dict:
    path = tmp_path / f"{name}.fastq.gz"
    path.write_bytes(b"")
    return {
        "sample": name,
        "fastq_1": str(path.resolve()),
        "layout": "SE",
        "assay": "chipseq",
        "target": "CTCF",
        "peak_mode": "narrow",
    }


def _validate(client, app, tmp_path, project_id, revision_ids, **extra):
    body = {
        "config": {},
        "samples": [_sample_row(tmp_path, f"row-{index}") for index, _ in enumerate(revision_ids)],
        "project_id": project_id,
        "sample_revision_ids": list(revision_ids),
        "reference_profile_revision_id": app.state.test_reference_profile.revision_id,
    }
    body.update(extra)
    return client.post(
        f"/api/v1/workflows/{WORKFLOW_ID}/validate",
        json=body,
    )


def test_validate_freezes_reversed_set_and_traces_published_rows(
    tmp_path,
    reference_ready_app,
) -> None:
    app = reference_ready_app
    registry = DataRegistryService(
        repository=app.state.persistence.data_registry_repository
    )
    project = registry.create_project("Treg project")
    other = registry.create_project("Other project")
    with ApiTestClient(app) as client:
        schema = app.openapi()
        assert set(schema["components"]["schemas"]["RunCreateRequest"]["properties"]) == {
            "snapshot_id",
            "tags",
        }
        assert set(
            schema["components"]["schemas"]["AuthoringValidateRequest"]["properties"]
        ) == {"expected_revision"}
        validation_request = schema["components"]["schemas"]["ValidationRequest"]
        assert "experiment_revision_id" in validation_request["properties"]
        assert "experiment_revision_id" not in validation_request.get("required", [])

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
        target = term("target", "H3K27ac", "target")
        batch = term("batch", "batch-1", "batch")
        role = term("replicate_role", "biological replicate", "role")
        imported_samples = []
        for name in ("alpha", "beta"):
            imported_samples.append(
                registry.import_sample(
                    project.project_id,
                    stable_key=name,
                    display_name=name,
                    attributes={"organism": "human"},
                )
            )
        samples = [
            item.revision.sample_revision_id for item in imported_samples
        ]
        other_sample = registry.import_sample(
            other.project_id,
            stable_key="other",
            display_name="other",
            attributes={"organism": "human"},
        ).revision.sample_revision_id

        def experiment(key: str, revision_ids: list[str], project_id: str) -> dict:
            response = _post(
                client,
                f"/api/v1/catalog/projects/{project_id}/experiments",
                {
                    "stable_key": key,
                    "display_name": key,
                    "assay": _known(assay),
                    "target": _known(target),
                    "batch": _known(batch),
                    "samples": [
                        {
                            "sample_revision_id": revision_id,
                            "replicate_role": _known(role),
                        }
                        for revision_id in revision_ids
                    ],
                },
                f"experiment-{key}",
            )
            assert response.status_code == 200, response.text
            return response.json()

        created = experiment("cut-pair", samples, project.project_id)
        decoy = experiment("other-pair", [other_sample], other.project_id)
        experiment_revision_id = created["experiment_revision_id"]
        before = _count(app, "validated_input_snapshots")

        omitted = _validate(client, app, tmp_path, project.project_id, samples)
        assert omitted.status_code == 200, omitted.text
        assert omitted.json()["snapshot"]["experiment_binding"] == {
            "status": "unresolved"
        }
        omitted_snapshot = omitted.json()["snapshot"]["snapshot_id"]
        unbound_run = client.post(
            f"/api/v1/workflows/{WORKFLOW_ID}/runs",
            json={"snapshot_id": omitted_snapshot},
        )
        assert unbound_run.status_code == 201, unbound_run.text
        assert _count(app, "run_experiment_bindings") == 0

        reversed_ids = list(reversed(samples))
        frozen = _validate(
            client,
            app,
            tmp_path,
            project.project_id,
            reversed_ids,
            experiment_revision_id=experiment_revision_id,
        )
        assert frozen.status_code == 200, frozen.text
        binding = frozen.json()["snapshot"]["experiment_binding"]
        assert binding["status"] == "experiment_snapshot"
        assert binding["experiment_revision_id"] == experiment_revision_id
        stored = _binding(app, frozen.json()["snapshot"]["snapshot_id"])
        assert stored["experiment_revision_id"] == experiment_revision_id
        assert stored["payload_digest"] == binding["payload_digest"]

        missing = _validate(
            client,
            app,
            tmp_path,
            project.project_id,
            samples[:1],
            experiment_revision_id=experiment_revision_id,
        )
        assert missing.status_code == 409
        assert missing.json()["issues"][0]["code"] == "catalog_sample_set_mismatch"
        assert missing.json()["issues"][0]["path"] == "experiment_revision_id"
        assert missing.json()["snapshot"] is None
        assert _count(app, "validated_input_snapshots") == before + 2

        newer = registry.revise_sample(
            imported_samples[0].sample.sample_id,
            display_name="alpha revised",
            attributes={"organism": "human", "note": "revised"},
        )
        different_revision = _validate(
            client,
            app,
            tmp_path,
            project.project_id,
            [newer.sample_revision_id, samples[1]],
            experiment_revision_id=experiment_revision_id,
        )
        assert different_revision.status_code == 409
        assert (
            different_revision.json()["issues"][0]["code"]
            == "catalog_sample_set_mismatch"
        )
        assert _count(app, "validated_input_snapshots") == before + 2

        cross_project = _validate(
            client,
            app,
            tmp_path,
            other.project_id,
            [other_sample],
            experiment_revision_id=experiment_revision_id,
        )
        assert cross_project.status_code == 404
        assert cross_project.json()["issues"][0]["code"] == "catalog_not_found"

        unknown = _validate(
            client,
            app,
            tmp_path,
            project.project_id,
            samples,
            experiment_revision_id="expr_" + "0" * 32,
        )
        assert unknown.status_code == 404
        assert unknown.json()["issues"][0]["code"] == "catalog_not_found"
        assert _count(app, "validated_input_snapshots") == before + 2

        archived = _post(
            client,
            f"/api/v1/catalog/experiments/{decoy['experiment_id']}/archive",
            {"expected_revision_number": 1},
            "archive-decoy",
        )
        assert archived.status_code == 200, archived.text
        archived_validation = _validate(
            client,
            app,
            tmp_path,
            other.project_id,
            [other_sample],
            experiment_revision_id=decoy["experiment_revision_id"],
        )
        assert archived_validation.status_code == 409
        assert (
            archived_validation.json()["issues"][0]["code"]
            == "catalog_experiment_archived"
        )
        assert _count(app, "validated_input_snapshots") == before + 2

        snapshot_id = frozen.json()["snapshot"]["snapshot_id"]
        created_run = client.post(
            f"/api/v1/workflows/{WORKFLOW_ID}/runs",
            json={"snapshot_id": snapshot_id, "tags": {"owner": "lab"}},
        )
        assert created_run.status_code == 201, created_run.text
        run_id = created_run.json()["run"]["run_id"]
        with app.state.persistence.engine.connect() as connection:
            copied = (
                connection.execute(
                    text(
                        "SELECT experiment_revision_id, payload_digest "
                        "FROM run_experiment_bindings WHERE run_id=:run_id"
                    ),
                    {"run_id": run_id},
                )
                .mappings()
                .one()
            )
        assert dict(copied) == dict(stored)

        trace_url = (
            f"/api/v1/catalog/projects/{project.project_id}/experiment-revisions/"
            f"{experiment_revision_id}/run-results"
        )
        unpublished = client.get(trace_url)
        assert unpublished.status_code == 200, unpublished.text
        unpublished_body = unpublished.json()
        assert unpublished_body["trace_status"] == "experiment_snapshot"
        assert "reason" not in unpublished_body
        assert unpublished_body["results"] == []
        CatalogExperimentRunListResponse.model_validate(unpublished_body)

        repository = app.state.persistence.repository
        current = repository.get_run(run_id)
        succeeded_at = datetime.now(timezone.utc)
        repository.update_run(
            replace(
                current,
                status=RunStatus.SUCCEEDED,
                updated_at=succeeded_at,
                ended_at=succeeded_at,
            ),
            expected_status=RunStatus.CREATED,
            event=RunEventDraft(
                event_type="status_changed",
                message="Scientific execution succeeded.",
                status=RunStatus.SUCCEEDED,
            ),
        )
        attempt_id = "resultattempt-" + sha256(b"artifacts").hexdigest()
        repository.begin_artifact_result_attempt(
            run_id,
            attempt_id=attempt_id,
            expected_status=RunStatus.SUCCEEDED,
        )
        repository.replace_artifacts(
            run_id,
            (
                RunArtifactRef(
                    artifact_id="artifact-1",
                    run_id=run_id,
                    artifact_type="file",
                    name="peaks.bed",
                    uri=f"run://runs/{run_id}/artifacts/artifact-1",
                    mime_type="text/plain",
                    produced_at=succeeded_at,
                    revision="artifactrev-" + sha256(b"artifact-1").hexdigest(),
                    metadata={"output_type": "peaks"},
                ),
            ),
            attempt_id=attempt_id,
            expected_status=RunStatus.SUCCEEDED,
            event=RunEventDraft(
                event_type="artifacts_indexed",
                message="Workflow artifacts indexed.",
                status=RunStatus.SUCCEEDED,
            ),
        )
        published = client.get(trace_url)
        assert published.status_code == 200, published.text
        published_body = published.json()
        assert published_body["trace_status"] == "experiment_snapshot"
        assert "reason" not in published_body
        assert len(published_body["results"]) == 1
        assert published_body["results"][0]["run_id"] == run_id
        assert published_body["results"][0]["publication_evidence"] == "present"
        assert (
            published_body["results"][0]["artifact_association_kind"]
            == "run_input_scope"
        )
        CatalogExperimentRunListResponse.model_validate(published_body)

        replay_body = {"snapshot_id": snapshot_id, "tags": {"owner": "lab"}}
        with app.state.persistence.engine.begin() as connection:
            connection.execute(
                text("DROP TRIGGER trg_run_experiment_bindings_no_update")
            )
            connection.execute(
                text(
                    "UPDATE run_experiment_bindings "
                    "SET payload_digest=:digest WHERE run_id=:run_id"
                ),
                {"digest": "b" * 64, "run_id": run_id},
            )
        digest_replay = client.post(
            f"/api/v1/workflows/{WORKFLOW_ID}/runs",
            json=replay_body,
        )
        assert digest_replay.status_code == 409
        assert (
            digest_replay.json()["issues"][0]["code"]
            == "VALIDATED_SNAPSHOT_REPLAY_CONFLICT"
        )
        with app.state.persistence.engine.begin() as connection:
            connection.execute(
                text(
                    "UPDATE run_experiment_bindings "
                    "SET experiment_revision_id=:revision_id, "
                    "payload_digest=:digest WHERE run_id=:run_id"
                ),
                {
                    "revision_id": decoy["experiment_revision_id"],
                    "digest": stored["payload_digest"],
                    "run_id": run_id,
                },
            )
        revision_replay = client.post(
            f"/api/v1/workflows/{WORKFLOW_ID}/runs",
            json=replay_body,
        )
        assert revision_replay.status_code == 409
        assert (
            revision_replay.json()["issues"][0]["code"]
            == "VALIDATED_SNAPSHOT_REPLAY_CONFLICT"
        )
