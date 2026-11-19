import asyncio

from sqlalchemy import event

from api_test_client import ApiTestClient
from encode_pipeline.services.data_registry import DataRegistryService
from test_publication_experiment_trace_http import (
    WORKFLOW_ID,
    _context_url,
    _publish,
    _seed_member,
    _trace_url,
    _validate,
)


def test_published_result_database_reads_run_off_the_event_loop(
    tmp_path, reference_ready_app
):
    app = reference_ready_app
    registry = DataRegistryService(
        repository=app.state.persistence.data_registry_repository
    )
    project = registry.create_project("Thread check")
    sample = registry.import_sample(
        project.project_id, stable_key="sample", display_name="Sample", attributes={}
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
        _publish(app, run, "published")
        generation = app.state.run_service.get_result_state(run).artifact_generation
        _seed_member(app)
        app.state.test_auth_tokens = ("c" * 43, "d" * 43)

        # Observe actual publication SELECTs without delaying or modifying a request.
        on_loop = []

        def observe(connection, cursor, statement, parameters, context, executemany):
            if "from artifact_publications" in statement.lower():
                try:
                    asyncio.get_running_loop()
                except RuntimeError:
                    on_loop.append(False)
                else:
                    on_loop.append(True)

        engine = app.state.persistence.engine
        event.listen(engine, "before_cursor_execute", observe)
        try:
            for url in (_context_url, _trace_url):
                on_loop.clear()
                response = client.get(
                    url(project.project_id, run), params={"generation": generation}
                )
                print(
                    url.__name__,
                    response.status_code,
                    "publication SELECTs on event loop:",
                    on_loop,
                )
                assert response.status_code == 200
                assert on_loop and not any(on_loop), (url.__name__, on_loop)
        finally:
            event.remove(engine, "before_cursor_execute", observe)
