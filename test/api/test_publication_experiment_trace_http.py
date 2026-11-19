"""Reverse trace from one published result back to its copied experiment revision."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timedelta, timezone
from decimal import Decimal
from hashlib import sha256

import pytest
from sqlalchemy import text

pytest.importorskip("fastapi")

from api_test_client import ApiTestClient
from encode_pipeline.api.catalog_models import (
    CatalogPublicationExperimentTraceResponse,
)
from encode_pipeline.api.models import PublishedResultContextResponse
from encode_pipeline.platform.authentication import UserAccount, UserRole, UserStatus
from encode_pipeline.platform.results import Issue, IssueSeverity
from encode_pipeline.platform.runs import (
    RunArtifactRef,
    RunQcMetric,
    RunStatus,
    build_qc_metric_id,
)
from encode_pipeline.services.authentication import SessionSecrets, new_session_record
from encode_pipeline.services.data_registry import DataRegistryService
from encode_pipeline.services.run_repositories import RunEventDraft


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
TRACE_PATH = (
    "/api/v1/catalog/projects/{project_id}/artifact-publications/"
    "{run_id}/{artifact_id}/experiment-trace"
)
CONTEXT_PATH = (
    "/api/v1/agent/projects/{project_id}/artifact-publications/"
    "{run_id}/{artifact_id}/context"
)


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


def _status(app, run_id: str) -> str:
    with app.state.persistence.engine.connect() as connection:
        return connection.execute(
            text("SELECT status FROM runs WHERE run_id = :run_id"),
            {"run_id": run_id},
        ).scalar_one()


def _generation(app, run_id: str, artifact_id: str) -> str:
    with app.state.persistence.engine.connect() as connection:
        return connection.execute(
            text(
                "SELECT artifact_generation FROM artifact_publications "
                "WHERE run_id = :run_id AND artifact_id = :artifact_id"
            ),
            {"run_id": run_id, "artifact_id": artifact_id},
        ).scalar_one()


def _copied_binding(app, run_id: str):
    with app.state.persistence.engine.connect() as connection:
        return (
            connection.execute(
                text(
                    "SELECT experiment_revision_id, payload_digest "
                    "FROM run_experiment_bindings WHERE run_id = :run_id"
                ),
                {"run_id": run_id},
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
        "samples": [
            _sample_row(tmp_path, f"row-{index}")
            for index, _ in enumerate(revision_ids)
        ],
        "project_id": project_id,
        "sample_revision_ids": list(revision_ids),
        "reference_profile_revision_id": app.state.test_reference_profile.revision_id,
    }
    body.update(extra)
    return client.post(f"/api/v1/workflows/{WORKFLOW_ID}/validate", json=body)


def _publish(app, run_id: str, artifact_name: str) -> None:
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
    attempt_id = "resultattempt-" + sha256(artifact_name.encode()).hexdigest()
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
                name=artifact_name,
                uri=f"run://runs/{run_id}/artifacts/artifact-1",
                mime_type="text/plain",
                produced_at=succeeded_at,
                revision="artifactrev-" + sha256(artifact_name.encode()).hexdigest(),
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


def _trace_url(project_id: str, run_id: str, artifact_id: str = "artifact-1") -> str:
    return TRACE_PATH.format(
        project_id=project_id,
        run_id=run_id,
        artifact_id=artifact_id,
    )


def _context_url(project_id: str, run_id: str, artifact_id: str = "artifact-1") -> str:
    return CONTEXT_PATH.format(
        project_id=project_id,
        run_id=run_id,
        artifact_id=artifact_id,
    )


def _seed_member(app) -> None:
    now = datetime.now(timezone.utc)
    account = UserAccount(
        user_id="usr_" + "1" * 32,
        username="trace-member",
        role=UserRole.MEMBER,
        status=UserStatus.ENABLED,
        password_hash=(
            "$argon2id$v=19$m=65536,t=3,p=4$"
            "pGLQVK/BOQLC0oPSA8RQTg$TdiQWUP9gwBAI8iAXUT7oEtjOPqKjJhqyW0JS8ye/Ag"
        ),
        created_at=now,
        updated_at=now,
        password_changed_at=now,
    )
    app.state.authentication_repository.create_account(account)
    secrets = SessionSecrets._from_generated("c" * 43, "d" * 43)
    app.state.authentication_repository.create_session(
        new_session_record(
            user_id=account.user_id,
            secrets=secrets,
            created_at=now,
            lifetime=timedelta(hours=8),
        )
    )


def test_published_result_traces_the_copied_revision_and_reads_context(
    tmp_path,
    reference_ready_app,
) -> None:
    app = reference_ready_app
    registry = DataRegistryService(
        repository=app.state.persistence.data_registry_repository
    )
    project = registry.create_project("Treg project")
    other = registry.create_project("Other project")
    _seed_member(app)
    with ApiTestClient(app) as client:
        schema = app.openapi()
        operation_paths = {
            operation_id: path
            for path, methods in schema["paths"].items()
            for method, operation in methods.items()
            if method in {"get", "post", "put", "patch", "delete"}
            for operation_id in [operation["operationId"]]
        }
        assert (
            operation_paths["catalogGetArtifactPublicationExperimentTrace"]
            == TRACE_PATH
        )
        assert operation_paths["readPublishedResultContext"] == CONTEXT_PATH
        for operation_id in (
            "catalogGetArtifactPublicationExperimentTrace",
            "readPublishedResultContext",
        ):
            path = operation_paths[operation_id]
            assert list(schema["paths"][path]) == ["get"]

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
        imported = [
            registry.import_sample(
                project.project_id,
                stable_key=name,
                display_name=name,
                attributes={"organism": "human"},
            )
            for name in ("alpha", "beta")
        ]
        samples = [item.revision.sample_revision_id for item in imported]
        created = _post(
            client,
            f"/api/v1/catalog/projects/{project.project_id}/experiments",
            {
                "stable_key": "cut-pair",
                "display_name": "cut-pair",
                "assay": _known(assay),
                "target": _known(target),
                "batch": _known(batch),
                "samples": [
                    {
                        "sample_revision_id": revision_id,
                        "replicate_role": _known(role),
                    }
                    for revision_id in samples
                ],
            },
            "experiment-cut-pair",
        )
        assert created.status_code == 200, created.text
        experiment = created.json()
        experiment_revision_id = experiment["experiment_revision_id"]

        omitted = _validate(client, app, tmp_path, project.project_id, samples)
        assert omitted.status_code == 200, omitted.text
        unbound_run = client.post(
            f"/api/v1/workflows/{WORKFLOW_ID}/runs",
            json={"snapshot_id": omitted.json()["snapshot"]["snapshot_id"]},
        )
        assert unbound_run.status_code == 201, unbound_run.text
        unbound_run_id = unbound_run.json()["run"]["run_id"]

        frozen = _validate(
            client,
            app,
            tmp_path,
            project.project_id,
            list(reversed(samples)),
            experiment_revision_id=experiment_revision_id,
        )
        assert frozen.status_code == 200, frozen.text
        copied_digest = frozen.json()["snapshot"]["experiment_binding"][
            "payload_digest"
        ]
        bound_run = client.post(
            f"/api/v1/workflows/{WORKFLOW_ID}/runs",
            json={
                "snapshot_id": frozen.json()["snapshot"]["snapshot_id"],
                "tags": {"owner": "lab"},
            },
        )
        assert bound_run.status_code == 201, bound_run.text
        run_id = bound_run.json()["run"]["run_id"]
        assert dict(_copied_binding(app, run_id)) == {
            "experiment_revision_id": experiment_revision_id,
            "payload_digest": copied_digest,
        }
        _publish(app, unbound_run_id, "unbound-peaks")
        _publish(app, run_id, "bound-peaks")
        generation = _generation(app, run_id, "artifact-1")
        unbound_generation = _generation(app, unbound_run_id, "artifact-1")

        bound = client.get(
            _trace_url(project.project_id, run_id),
            params={"generation": generation},
        )
        assert bound.status_code == 200, bound.text
        bound_body = bound.json()
        CatalogPublicationExperimentTraceResponse.model_validate(bound_body)
        assert bound_body["trace_status"] == "experiment_snapshot"
        assert bound_body["experiment_revision_id"] == experiment_revision_id
        assert bound_body["payload_digest"] == copied_digest
        assert "reason" not in bound_body
        result = bound_body["result"]
        assert result["artifact_association_kind"] == "run_input_scope"
        assert result["publication_evidence"] == "present"
        assert {
            sample["sample_revision_id"]
            for sample in result["run_sample_binding"]["associated_run_samples"]
        } == set(samples)
        assert "sample_id" not in result

        revised = _post(
            client,
            f"/api/v1/catalog/experiments/{experiment['experiment_id']}/revisions",
            {
                "expected_revision_number": 1,
                "display_name": "cut-pair revised",
                "assay": _known(assay),
                "target": _known(target),
                "batch": _known(batch),
                "samples": [
                    {
                        "sample_revision_id": revision_id,
                        "replicate_role": _known(role),
                    }
                    for revision_id in samples
                ],
            },
            "revise-cut-pair",
        )
        assert revised.status_code == 200, revised.text
        revised_id = revised.json()["experiment_revision_id"]
        assert revised_id != experiment_revision_id
        pinned = client.get(
            _trace_url(project.project_id, run_id),
            params={"generation": generation},
        )
        assert pinned.status_code == 200, pinned.text
        assert pinned.json()["experiment_revision_id"] == experiment_revision_id
        assert pinned.json()["payload_digest"] == copied_digest
        old_forward = client.get(
            f"/api/v1/catalog/projects/{project.project_id}/experiment-revisions/"
            f"{experiment_revision_id}/run-results"
        )
        assert old_forward.status_code == 200, old_forward.text
        assert old_forward.json()["trace_status"] == "experiment_snapshot"
        assert [row["run_id"] for row in old_forward.json()["results"]] == [run_id]
        new_forward = client.get(
            f"/api/v1/catalog/projects/{project.project_id}/experiment-revisions/"
            f"{revised_id}/run-results"
        )
        assert new_forward.status_code == 200, new_forward.text
        assert new_forward.json()["trace_status"] == "unresolved"
        assert new_forward.json()["reason"] == "experiment_snapshot_binding_absent"
        assert new_forward.json()["results"] == []

        unbound = client.get(
            _trace_url(project.project_id, unbound_run_id),
            params={"generation": unbound_generation},
        )
        assert unbound.status_code == 200, unbound.text
        unbound_body = unbound.json()
        CatalogPublicationExperimentTraceResponse.model_validate(unbound_body)
        assert unbound_body["trace_status"] == "unresolved"
        assert unbound_body["reason"] == "experiment_snapshot_binding_absent"
        assert "experiment_revision_id" not in unbound_body
        assert unbound_body["result"]["run_id"] == unbound_run_id
        assert unbound_body["result"]["artifact_association_kind"] == "run_input_scope"

        metric_key = "alignment.mapped_fraction"
        artifacts = app.state.persistence.repository.list_artifacts(run_id)
        app.state.run_service.replace_qc_metrics(
            run_id,
            (
                RunQcMetric(
                    metric_id=build_qc_metric_id(metric_key, "run", None, None),
                    run_id=run_id,
                    metric_key=metric_key,
                    display_name="Mapped fraction",
                    value=Decimal("0.5"),
                    unit="fraction",
                    scope="run",
                    source_artifact_id=artifacts[0].artifact_id,
                    produced_at=datetime.now(timezone.utc),
                ),
            ),
            expected_artifacts=artifacts,
        )
        app.state.run_service.add_event(
            run_id,
            "preflight_failed",
            "Reference profile could not be prepared.",
            stage="preflight",
            issue=Issue(
                code="PREFLIGHT_FAILED",
                message="Reference profile could not be prepared.",
                severity=IssueSeverity.ERROR,
            ),
        )
        bindings_before = _count(app, "run_experiment_bindings")
        events_before = _count(app, "run_events")
        status_before = _status(app, run_id)

        app.state.test_auth_tokens = ("c" * 43, "d" * 43)
        member_context = client.get(
            _context_url(project.project_id, run_id),
            params={"generation": generation},
        )
        assert member_context.status_code == 200, member_context.text
        app.state.test_auth_tokens = ("a" * 43, "b" * 43)
        context = client.get(
            _context_url(project.project_id, run_id),
            params={"generation": generation},
        )
        assert context.status_code == 200, context.text
        assert context.json() == member_context.json()
        context_body = context.json()
        PublishedResultContextResponse.model_validate(context_body)
        assert context_body["experiment_trace"]["experiment_revision_id"] == (
            experiment_revision_id
        )
        assert context_body["experiment_trace"]["payload_digest"] == copied_digest
        run_record = client.get(f"/api/v1/runs/{run_id}")
        assert run_record.status_code == 200, run_record.text
        recorded = run_record.json()["run"]
        assert "experiment_binding" not in recorded
        assert "experiment_revision_id" not in recorded
        assert context_body["run_status"] == recorded["status"]
        assert context_body["reference_evidence"] == "present"
        assert context_body["reference_profile"] == recorded["reference_profile"]
        qc = client.get(f"/api/v1/runs/{run_id}/qc-metrics")
        assert qc.status_code == 200, qc.text
        assert context_body["qc_evidence"] == "present"
        assert context_body["qc_metrics"] == qc.json()["qc_metrics"]
        assert context_body["qc_metrics"][0]["metric_key"] == metric_key
        assert context_body["qc_metrics"][0]["value"] == "0.5"
        assert context_body["preflight_failure_evidence"] == "present"
        assert context_body["preflight_failures"] == [
            {
                "event_id": context_body["preflight_failures"][0]["event_id"],
                "event_type": "preflight_failed",
                "message": "Reference profile could not be prepared.",
                "stage": "preflight",
                "code": "PREFLIGHT_FAILED",
            }
        ]

        unbound_context = client.get(
            _context_url(project.project_id, unbound_run_id),
            params={"generation": unbound_generation},
        )
        assert unbound_context.status_code == 200, unbound_context.text
        unbound_context_body = unbound_context.json()
        PublishedResultContextResponse.model_validate(unbound_context_body)
        assert unbound_context_body["experiment_trace"]["reason"] == (
            "experiment_snapshot_binding_absent"
        )
        assert "experiment_revision_id" not in unbound_context_body["experiment_trace"]
        assert unbound_context_body["qc_evidence"] == "absent"
        assert "qc_metrics" not in unbound_context_body
        assert unbound_context_body["preflight_failure_evidence"] == "absent"
        assert unbound_context_body["reference_evidence"] == "present"

        wrong_project = client.get(
            _trace_url(other.project_id, run_id),
            params={"generation": generation},
        )
        wrong_context = client.get(
            _context_url(other.project_id, run_id),
            params={"generation": generation},
        )
        for response in (wrong_project, wrong_context):
            assert response.status_code == 404, response.text
            assert response.json()["issues"][0]["code"] == "catalog_not_found"
            assert "run_status" not in response.text
            assert experiment_revision_id not in response.text
        missing_generation = client.get(
            _trace_url(project.project_id, run_id),
            params={"generation": "artifactgen-" + "0" * 64},
        )
        assert missing_generation.status_code == 404, missing_generation.text
        publication = client.get(
            f"/api/v1/artifact-publications/{run_id}/artifact-1",
            params={"generation": generation},
        )
        assert publication.status_code == 200, publication.text
        assert "experiment_revision_id" not in publication.text

        app.state.test_auth_tokens = None
        anonymous = client.get(
            _context_url(project.project_id, run_id),
            params={"generation": generation},
        )
        assert anonymous.status_code == 401
        assert "AUTHENTICATION_REQUIRED" in anonymous.text
        app.state.test_auth_tokens = ("a" * 43, "b" * 43)

        with app.state.persistence.engine.begin() as connection:
            connection.execute(
                text("DROP TRIGGER trg_run_experiment_bindings_no_update")
            )
            connection.execute(
                text(
                    "UPDATE run_experiment_bindings "
                    "SET payload_digest = :digest WHERE run_id = :run_id"
                ),
                {"digest": "0" * 64, "run_id": run_id},
            )
        inconsistent = client.get(
            _trace_url(project.project_id, run_id),
            params={"generation": generation},
        )
        assert inconsistent.status_code == 200, inconsistent.text
        inconsistent_body = inconsistent.json()
        CatalogPublicationExperimentTraceResponse.model_validate(inconsistent_body)
        assert inconsistent_body["trace_status"] == "unresolved"
        assert inconsistent_body["reason"] == (
            "experiment_revision_evidence_inconsistent"
        )
        assert "experiment_revision_id" not in inconsistent_body
        assert "payload_digest" not in inconsistent_body
        assert inconsistent_body["result"]["run_id"] == run_id
        assert inconsistent_body["result"]["artifact_association_kind"] == (
            "run_input_scope"
        )
        forward_after = client.get(
            f"/api/v1/catalog/projects/{project.project_id}/experiment-revisions/"
            f"{experiment_revision_id}/run-results"
        )
        assert forward_after.status_code == 200, forward_after.text
        assert forward_after.json()["trace_status"] == "experiment_snapshot"
        assert forward_after.json()["results"] == []

        assert _count(app, "run_experiment_bindings") == bindings_before
        assert _count(app, "run_events") == events_before
        assert _status(app, run_id) == status_before
