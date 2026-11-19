from sqlalchemy import text

from api_test_client import ApiTestClient
from encode_pipeline.services.data_registry import DataRegistryService
from test_publication_experiment_trace_http import (
    WORKFLOW_ID,
    _context_url,
    _mapped_fraction,
    _publish,
    _replace_peaks,
    _seed_member,
    _validate,
)


def test_current_generation_with_unbound_persisted_qc(tmp_path, reference_ready_app):
    app = reference_ready_app
    registry = DataRegistryService(
        repository=app.state.persistence.data_registry_repository
    )
    project = registry.create_project("QC evidence fixture")
    sample = registry.import_sample(
        project.project_id,
        stable_key="one",
        display_name="One",
        attributes={},
    )
    with ApiTestClient(app) as client:
        snapshot = _validate(
            client,
            app,
            tmp_path,
            project.project_id,
            [sample.revision.sample_revision_id],
        )
        run = client.post(
            f"/api/v1/workflows/{WORKFLOW_ID}/runs",
            json={"snapshot_id": snapshot.json()["snapshot"]["snapshot_id"]},
        ).json()["run"]["run_id"]
        _publish(app, run, "generation-one")
        service = app.state.run_service
        (artifact,) = service.list_artifacts(run)
        gen = service.get_result_state(run).artifact_generation
        url = _context_url(project.project_id, run)
        _seed_member(app)
        app.state.test_auth_tokens = ("c" * 43, "d" * 43)
        assert (
            client.get(url, params={"generation": gen}).json()["qc_evidence"]
            == "absent"
        )
        service.replace_qc_metrics(
            run,
            (_mapped_fraction(run, artifact.artifact_id, "0.5"),),
            expected_artifacts=(artifact,),
        )
        assert (
            client.get(url, params={"generation": gen}).json()["qc_metrics"][0]["value"]
            == "0.5"
        )
        with app.state.persistence.engine.begin() as db:
            db.execute(
                text(
                    "UPDATE run_result_states SET qc_generation=NULL, "
                    "qc_manifest_digest=NULL, qc_artifact_generation=NULL "
                    "WHERE run_id=:run"
                ),
                {"run": run},
            )
            assert (
                db.execute(
                    text("SELECT count(*) FROM run_qc_metrics WHERE run_id=:run"),
                    {"run": run},
                ).scalar_one()
                == 1
            )
        qc = client.get(f"/api/v1/runs/{run}/qc-metrics")
        context = client.get(url, params={"generation": gen})
        assert (qc.status_code, qc.json()["issues"][0]["code"]) == (
            500,
            "RUN_QC_METRIC_DATA_INVALID",
        )
        assert (context.status_code, context.json().get("code")) == (
            500,
            "RUN_QC_METRIC_DATA_INVALID",
        )
        assert context.json().get("qc_evidence") != "absent"


def test_superseded_generation_stays_unverifiable(tmp_path, reference_ready_app):
    app = reference_ready_app
    registry = DataRegistryService(
        repository=app.state.persistence.data_registry_repository
    )
    project = registry.create_project("QC generation fixture")
    sample = registry.import_sample(
        project.project_id,
        stable_key="one",
        display_name="One",
        attributes={},
    )
    with ApiTestClient(app) as client:
        snapshot = _validate(
            client,
            app,
            tmp_path,
            project.project_id,
            [sample.revision.sample_revision_id],
        )
        run = client.post(
            f"/api/v1/workflows/{WORKFLOW_ID}/runs",
            json={"snapshot_id": snapshot.json()["snapshot"]["snapshot_id"]},
        ).json()["run"]["run_id"]
        _publish(app, run, "generation-one")
        service = app.state.run_service
        (artifact,) = service.list_artifacts(run)
        first_generation = service.get_result_state(run).artifact_generation
        service.replace_qc_metrics(
            run,
            (_mapped_fraction(run, artifact.artifact_id, "0.5"),),
            expected_artifacts=(artifact,),
        )
        _replace_peaks(app, run, "generation-two")
        (replacement,) = service.list_artifacts(run)
        service.replace_qc_metrics(
            run,
            (_mapped_fraction(run, replacement.artifact_id, "0.9"),),
            expected_artifacts=(replacement,),
        )
        second_generation = service.get_result_state(run).artifact_generation
        assert second_generation != first_generation
        _seed_member(app)
        app.state.test_auth_tokens = ("c" * 43, "d" * 43)
        url = _context_url(project.project_id, run)
        superseded = client.get(url, params={"generation": first_generation})
        current = client.get(url, params={"generation": second_generation})
        assert superseded.status_code == 200, superseded.text
        assert current.status_code == 200, current.text
        assert superseded.json()["qc_evidence"] == "unverifiable"
        assert "qc_metrics" not in superseded.json()
        assert current.json()["qc_metrics"][0]["value"] == "0.9"
