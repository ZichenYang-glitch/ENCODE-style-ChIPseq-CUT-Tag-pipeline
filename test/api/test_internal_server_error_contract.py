"""Actual 500 envelopes and their operation-specific OpenAPI declarations."""

from __future__ import annotations

import asyncio
from dataclasses import replace
from datetime import datetime, timedelta, timezone
import json
from types import SimpleNamespace

import httpx
from jsonschema import Draft202012Validator
import pytest

from api_test_client import seeded_auth_async_client
from conftest import seed_test_authentication
from encode_pipeline.api.main import create_app
from encode_pipeline.platform.authentication import UserRole
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.services.authentication import (
    new_session_record,
    new_session_secrets,
)
from encode_pipeline.services.browser_tracks import BrowserTrackService
from encode_pipeline.services.run_repositories import RunEventDraft


async def _application_response(app, request: httpx.Request) -> httpx.Response:
    """Read the status and body the ASGI app actually produced."""

    body_parts: list[bytes] = []
    status_code = 500
    headers: list[tuple[bytes, bytes]] = []

    async def receive() -> dict[str, object]:
        if receive.sent:
            return {"type": "http.disconnect"}
        receive.sent = True
        return {"type": "http.request", "body": request.content, "more_body": False}

    receive.sent = False

    async def send(message: dict[str, object]) -> None:
        nonlocal status_code, headers
        if message["type"] == "http.response.start":
            status_code = int(message["status"])
            headers = list(message.get("headers", []))
        elif message["type"] == "http.response.body":
            body_parts.append(message.get("body", b""))

    await app(
        {
            "type": "http",
            "asgi": {"version": "3.0"},
            "http_version": "1.1",
            "method": request.method,
            "headers": [(key.lower(), value) for key, value in request.headers.raw],
            "scheme": request.url.scheme,
            "path": request.url.path,
            "raw_path": request.url.raw_path.split(b"?")[0],
            "query_string": request.url.query,
            "server": (request.url.host, request.url.port or 80),
            "client": ("127.0.0.1", 123),
        },
        receive,
        send,
    )
    return httpx.Response(status_code, headers=headers, content=b"".join(body_parts))


WORKFLOW = "encode-style-chipseq-cuttag-atac-mnase"
RUN = "run-error-contract"
USER = "usr_" + "1" * 32
GENERATION = "artifactgen-" + "0" * 64
REVISION = "artifactrev-" + "0" * 64
PRIVATE = (
    "PRIVATE_EXCEPTION_MARKER",
    "/private/error-contract/input",
    "TOKEN=secret-marker",
)

# Service calls are interrupted before writes, queue, preflight or agent execution.
# The real router, authentication, error handlers and SQLite stay in place.
FAULTS = {
    "login": ("authentication_service", "login"),
    "logout": ("authentication_service", "logout"),
    "session_state": ("authentication_service", "setup_complete"),
    "get_terminal_email_preference": (
        "authentication_service",
        "get_terminal_email_preference",
    ),
    "set_terminal_email_preference": (
        "authentication_service",
        "set_terminal_email_enabled",
    ),
    "list_accounts": ("account_administration_service", "list_accounts"),
    "create_member_account": ("account_administration_service", "create_member"),
    "set_account_status": ("account_administration_service", "set_account_status"),
    "reset_account_password": ("account_administration_service", "reset_password"),
    "revoke_account_sessions": ("account_administration_service", "revoke_sessions"),
    "listWorkflows": ("registry", "list_metadata"),
    "getWorkflow": ("registry", "get"),
    "listCompatibleReferenceProfiles": ("reference_profile_service", "list_enabled"),
    "getWorkflowSchema": ("registry", "get"),
    "validateWorkflow": ("validated_input_service", "validate"),
    "chatWithWorkflowAgent": ("agent_service", "chat"),
    "createRun": ("validated_run_creation_service", "create_run"),
    "listRuns": ("run_service", "list_run_history"),
    "getRun": ("run_service", "get_run"),
    "startRun": ("run_submission_service", "start_run"),
    "cancelRun": ("run_cancellation_service", "cancel_run"),
    "listRunEvents": ("run_service", "list_events"),
    "listRunLogs": ("run_service", "list_logs"),
    "triggerPreflight": ("run_service", "transition_run"),
    "listRunArtifacts": ("run_service", "list_artifacts_page"),
    "getRunArtifact": ("run_service", "get_artifact_at_generation"),
    "downloadRunArtifact": ("artifact_download_service", "prepare"),
    "listArtifactPublications": (
        "artifact_publication_service",
        "list_artifact_publications",
    ),
    "getArtifactPublication": (
        "artifact_publication_service",
        "get_artifact_publication",
    ),
    "listRunQcMetrics": ("run_service", "list_qc_metrics_page"),
    "getRunBrowserSession": ("browser_track_service", "session"),
    "prepareRunBrowserTrack": ("browser_track_service", "prepare"),
    "summarizeRunBrowserTrackRegion": ("browser_track_service", "region_summary"),
    "readRunBrowserTrackBytes": ("browser_track_service", "open_bytes"),
    "headRunBrowserTrackBytes": ("browser_track_service", "open_bytes"),
    "authoringListProjects": ("authoring_service", "list_projects"),
    "authoringListProjectSamples": ("authoring_service", "list_samples"),
    "authoringListProjectInputRevisions": (
        "authoring_service",
        "list_input_revisions",
    ),
    "authoringListSessions": ("authoring_service", "list_sessions"),
    "authoringCreateSession": ("authoring_service", "create_session"),
    "authoringGetSession": ("authoring_service", "get_session"),
    "authoringListReceipts": ("authoring_service", "list_adoption_receipts"),
    "authoringSaveDraft": ("authoring_service", "save_draft"),
    "authoringOpenConnection": ("authoring_service", "open_connection"),
    "authoringCloseConnection": ("authoring_service", "close_connection"),
    "authoringCreateTurn": ("authoring_service", "start_turn"),
    "authoringGetTurn": ("authoring_service", "get_turn"),
    "authoringAdoptProposal": ("authoring_service", "adopt"),
    "authoringRejectProposal": ("authoring_service", "reject"),
    "authoringValidateDraft": ("authoring_service", "validate"),
    "catalogListRegisteredResults": (
        "registered_result_service",
        "list_current",
    ),
    "catalogListRegisteredResultRevisions": (
        "registered_result_service",
        "list_revisions",
    ),
    "catalogGetRegisteredResultRevision": (
        "registered_result_service",
        "get_revision",
    ),
    "catalogDownloadRegisteredResultRevision": (
        "registered_result_service",
        "prepare_download",
    ),
    "catalogRegisterRegisteredResult": (
        "registered_result_service",
        "register",
    ),
    "catalogReviseRegisteredResult": (
        "registered_result_service",
        "revise",
    ),
    "catalogArchiveRegisteredResult": (
        "registered_result_service",
        "archive",
    ),
    "catalogListExperiments": ("experiment_catalog_service", "list_experiments"),
    "catalogListSamples": ("experiment_catalog_service", "list_samples"),
    "catalogListSampleRevisions": (
        "experiment_catalog_service",
        "list_sample_revisions",
    ),
    "catalogGetSampleRevision": ("experiment_catalog_service", "get_sample_revision"),
    "catalogListExperimentRevisions": (
        "experiment_catalog_service",
        "list_experiment_revisions",
    ),
    "catalogGetExperimentRevision": (
        "experiment_catalog_service",
        "get_experiment_revision",
    ),
    "catalogListTerms": ("experiment_catalog_service", "list_terms"),
    "catalogGetTerm": ("experiment_catalog_service", "get_term"),
    "catalogCreateTerm": ("experiment_catalog_service", "create_term"),
    "catalogReviseTerm": ("experiment_catalog_service", "revise_term"),
    "catalogArchiveTerm": ("experiment_catalog_service", "archive_term"),
    "catalogCreateExperiment": ("experiment_catalog_service", "create_experiment"),
    "catalogReviseExperiment": ("experiment_catalog_service", "revise_experiment"),
    "catalogArchiveExperiment": ("experiment_catalog_service", "archive_experiment"),
    "catalogReviseSampleCatalog": (
        "experiment_catalog_service",
        "revise_sample_catalog",
    ),
    "catalogListSampleRevisionRunResults": (
        "experiment_catalog_service",
        "list_sample_run_results",
    ),
    "catalogListExperimentRevisionRunResults": (
        "experiment_catalog_service",
        "list_experiment_run_results",
    ),
    "catalogRegisterProposal": (
        "experiment_catalog_service",
        "register_proposal",
    ),
    "catalogGetProposal": ("experiment_catalog_service", "get_proposal"),
    "catalogPreviewProposal": (
        "experiment_catalog_service",
        "preview_proposal",
    ),
    "catalogAdoptProposal": ("experiment_catalog_service", "adopt_proposal"),
    "catalogListAdoptionReceipts": (
        "experiment_catalog_service",
        "list_adoption_receipts",
    ),
    "catalogGetArtifactPublicationExperimentTrace": (
        "experiment_catalog_service",
        "get_publication_experiment_trace",
    ),
    "readPublishedResultContext": (
        "artifact_publication_service",
        "get_artifact_publication",
    ),
}
SPECIAL_MODELS = {
    "readPublishedResultContext": "PublishedResultContextErrorResponse",
    "createRun": "RunResponse",
    "listRuns": "RunHistoryResponse",
    "listRunArtifacts": "RunArtifactsResponse",
    "getRunArtifact": "RunArtifactDetailResponse",
    "downloadRunArtifact": "RunArtifactDownloadErrorResponse",
    "listArtifactPublications": "ArtifactPublicationListResponse",
    "getArtifactPublication": "ArtifactPublicationDetailResponse",
    "listRunQcMetrics": "RunQcMetricsResponse",
    "getRunBrowserSession": "BrowserErrorResponse",
    "prepareRunBrowserTrack": "BrowserErrorResponse",
    "summarizeRunBrowserTrackRegion": "BrowserErrorResponse",
    "readRunBrowserTrackBytes": "BrowserErrorResponse",
    "headRunBrowserTrackBytes": "BrowserErrorResponse",
    "authoringValidateDraft": "AuthoringValidationResponse",
}


@pytest.fixture
def error_api(tmp_path):
    app = create_app(
        database_url=f"sqlite:///{tmp_path / 'errors.db'}",
        workspace_root=tmp_path / "workspaces",
        project_root=tmp_path,
    )
    admin = seed_test_authentication(app)
    member = replace(admin, user_id=USER, username="error-member", role=UserRole.MEMBER)
    app.state.authentication_repository.create_account(member)
    secrets = new_session_secrets()
    app.state.authentication_repository.create_session(
        new_session_record(
            user_id=USER,
            secrets=secrets,
            created_at=datetime.now(timezone.utc),
            lifetime=timedelta(hours=1),
        )
    )
    admin_tokens = app.state.test_auth_tokens
    member_tokens = (secrets.session_token, secrets.csrf_token)
    app.state.browser_track_service = BrowserTrackService(
        run_service=app.state.run_service,
        download_service=app.state.artifact_download_service,
        workspace_root=app.state.workspace_root,
        reference_profile_service=app.state.reference_profile_service,
    )
    now = datetime.now(timezone.utc)
    record = RunRecord(
        run_id=RUN,
        workflow_id=WORKFLOW,
        inputs={},
        status=RunStatus.CREATED,
        created_at=now,
        updated_at=now,
        started_at=None,
        ended_at=None,
        current_stage=None,
        cancellation_reason=None,
        error=None,
    )
    app.state.persistence.repository.create_run(
        record,
        RunEventDraft(event_type="created", message="Created for contract test."),
    )
    spec = app.openapi()
    operations = {
        op["operationId"]: (method, path, op)
        for path, methods in spec["paths"].items()
        for method, op in methods.items()
        if "operationId" in op
    }
    responses = []

    async def request(operation):
        method, path, _ = operations[operation]
        path = path.format(
            workflow_id=WORKFLOW,
            run_id=RUN,
            user_id=USER,
            artifact_id="artifact-test",
            generation=GENERATION,
            revision=REVISION,
            file_token=REVISION,
            project_id="prj_" + "1" * 32,
            sample_id="smp_" + "1" * 32,
            sample_revision_id="smpr_" + "1" * 32,
            experiment_id="exp_" + "1" * 32,
            experiment_revision_id="expr_" + "1" * 32,
            term_id="term_" + "1" * 32,
            sid="aas_" + "ab" * 16,
            turn_id="atr_" + "cd" * 16,
            proposal_id="prop_" + "ef" * 16,
            result_id="rreg_" + "2" * 32,
            result_revision_id="rrev_" + "3" * 32,
        )
        payloads = {
            "login": {"username": "test-admin", "password": "private-payload-marker"},
            "create_member_account": {
                "username": "new-member",
                "password": "private-payload-marker",
            },
            "reset_account_password": {"password": "private-payload-marker"},
            "set_account_status": {"enabled": False},
            "set_terminal_email_preference": {"terminal_email_enabled": False},
            "validateWorkflow": {"config": {}},
            "chatWithWorkflowAgent": {"message": "private-payload-marker"},
            "createRun": {"snapshot_id": "vsnap_" + "0" * 32},
            "prepareRunBrowserTrack": {
                "generation": GENERATION,
                "revision": REVISION,
            },
            "authoringCreateSession": {
                "project_id": "prj_" + "1" * 32,
                "workflow_id": WORKFLOW,
            },
            "authoringSaveDraft": {"expected_revision": 0, "config": {}},
            "authoringCreateTurn": {
                "client_turn_id": "turn-1",
                "expected_revision": 0,
                "question": {
                    "intent": "suggest",
                    "field_tokens": ["config.threads"],
                    "preferences": [],
                },
            },
            "authoringAdoptProposal": {
                "expected_revision": 0,
                "proposal_digest": "a" * 64,
            },
            "authoringRejectProposal": {"proposal_digest": "b" * 64},
            "authoringValidateDraft": {"expected_revision": 0},
            "catalogRegisterRegisteredResult": {
                "stable_key": "contract-result",
                "display_name": "Contract result",
                "format": "tsv",
                "origin": {
                    "kind": "external",
                    "source": {"status": "unknown"},
                },
                "method": {"status": "unknown"},
                "reference": {"status": "unknown"},
                "pool_relative_path": "external/contract.tsv",
            },
            "catalogReviseRegisteredResult": {"expected_revision_number": 1},
            "catalogArchiveRegisteredResult": {"expected_revision_number": 1},
            "catalogCreateTerm": {"term_type": "tag", "display_name": "tag"},
            "catalogReviseTerm": {
                "expected_revision_number": 1,
                "change_kind": "editorial",
                "display_name": "tag",
            },
            "catalogArchiveTerm": {"expected_revision_number": 1},
            "catalogCreateExperiment": {
                "stable_key": "exp",
                "display_name": "exp",
                "assay": {},
                "target": {},
                "batch": {},
                "samples": [],
            },
            "catalogReviseExperiment": {
                "expected_revision_number": 1,
                "display_name": "exp",
                "assay": {},
                "target": {},
                "batch": {},
                "samples": [],
            },
            "catalogArchiveExperiment": {"expected_revision_number": 1},
            "catalogReviseSampleCatalog": {
                "expected_revision_number": 1,
                "catalog": None,
            },
            "catalogRegisterProposal": {
                "source": "offline_suggestion",
                "items": [{"item_id": "item_1", "target_kind": "experiment"}],
            },
            "catalogPreviewProposal": {"selected_item_ids": ["item_1"]},
            "catalogAdoptProposal": {
                "preview_digest": "a" * 64,
                "selected_item_ids": ["item_1"],
                "expected_revisions": {"exp_" + "1" * 32: 1},
            },
        }
        params = {}
        if operation in {
            "getRunArtifact",
            "downloadRunArtifact",
            "getArtifactPublication",
            "catalogGetArtifactPublicationExperimentTrace",
            "readPublishedResultContext",
        }:
            params["generation"] = GENERATION
        if operation == "downloadRunArtifact":
            params["revision"] = REVISION
        if operation == "summarizeRunBrowserTrackRegion":
            params.update(chrom="chr1", start=0, end=1)
        app.state.test_auth_tokens = (
            member_tokens if "terminal_email_preference" in operation else admin_tokens
        )
        headers = {}
        if operation.startswith("catalog") and method.lower() == "post":
            headers["Idempotency-Key"] = "contract-fault"
        async with seeded_auth_async_client(
            app,
            transport=httpx.ASGITransport(app=app, raise_app_exceptions=False),
            base_url="http://testserver",
        ) as client:
            if method.lower() == "head":
                # httpx drops HEAD bodies. Keep the JSON envelope the app emitted.
                prepared = client.build_request(
                    method,
                    path,
                    json=payloads.get(operation),
                    params=params,
                    headers=headers,
                )
                response = await _application_response(app, prepared)
            else:
                response = await client.request(
                    method,
                    path,
                    json=payloads.get(operation),
                    params=params,
                    headers=headers,
                )
        responses.append(
            {
                "operation": operation,
                "status": response.status_code,
                "body": response.json(),
            }
        )
        return response

    try:
        yield SimpleNamespace(
            app=app,
            spec=spec,
            operations=operations,
            record=record,
            request=lambda op: asyncio.run(request(op)),
        )
    finally:
        (tmp_path / "http-evidence.json").write_text(
            json.dumps(responses, indent=2) + "\n"
        )
        app.state.run_queue.close()
        try:
            assert app.state.persistence.engine.pool.checkedout() == 0
        finally:
            app.state.persistence.close()


def assert_declared_response(api, operation, response, *, explicit=False):
    assert response.status_code == 500
    body = response.json()
    assert body["ok"] is False
    for marker in (*PRIVATE, "private-payload-marker", "Traceback"):
        assert marker not in response.text
    if operation == "readPublishedResultContext":
        assert body["code"] == "INTERNAL_SERVER_ERROR"
        assert set(body) == {"ok", "code", "message"}
    else:
        assert body["issues"][0]["code"] == (
            "RUN_REFERENCE_EVIDENCE_INVALID" if explicit else "INTERNAL_SERVER_ERROR"
        )
        assert body["issues"][0].get("technical_message") is None
        # Assert null/omission and exact envelope keys independently of permissive models.
        if explicit or operation == "createRun":
            assert set(body) == {"ok", "run", "issues"} and body["run"] is None
        elif operation not in SPECIAL_MODELS:
            assert set(body) == {"ok", "workflow_id", "value", "snapshot", "issues"}
            assert body["value"] is None and body["snapshot"] is None
            assert body["workflow_id"] == (
                WORKFLOW if "{workflow_id}" in api.operations[operation][1] else None
            )
    declared = api.operations[operation][2]["responses"].get("500")
    assert declared is not None, f"{operation}: actual 500 has no declaration"
    schema = declared["content"]["application/json"]["schema"]
    if operation == "getRun":
        assert schema["anyOf"] == [
            {"$ref": "#/components/schemas/RunResponse"},
            {"$ref": "#/components/schemas/ValidationResponse"},
        ]
    else:
        assert schema == {
            "$ref": "#/components/schemas/"
            + SPECIAL_MODELS.get(operation, "ValidationResponse")
        }
    Draft202012Validator({**schema, "components": api.spec["components"]}).validate(
        body
    )


@pytest.mark.parametrize("operation", FAULTS)
def test_actual_internal_error_matches_operation_schema(
    error_api, monkeypatch, operation
):
    state_name, method = FAULTS[operation]
    target = getattr(error_api.app.state, state_name)
    calls = []

    def fail(*args, **kwargs):
        calls.append(True)
        # TypeError bypasses the reference-profile route's intentional 503 catch.
        error = (
            TypeError
            if operation == "listCompatibleReferenceProfiles"
            else RuntimeError
        )
        raise error(" | ".join(PRIVATE))

    async def fail_async(*args, **kwargs):
        fail(*args, **kwargs)

    # Registry is immutable/slotted; patch only the invoked method on its class.
    owner = type(target) if state_name == "registry" else target
    monkeypatch.setattr(
        owner, method, fail_async if operation == "chatWithWorkflowAgent" else fail
    )
    response = error_api.request(operation)
    assert calls == [True], response.text
    assert_declared_response(error_api, operation, response)


@pytest.mark.parametrize("operation", ("createRun", "getRun"))
def test_explicit_reference_failure_matches_operation_schema(
    error_api, monkeypatch, operation
):
    if operation == "createRun":
        # Stop at the creation-service seam; no snapshot capture or real run creation.
        monkeypatch.setattr(
            error_api.app.state.validated_run_creation_service,
            "create_run",
            lambda *a, **kw: SimpleNamespace(record=error_api.record, created=True),
        )

    def fail(*args, **kwargs):
        raise ValueError(" | ".join(PRIVATE))

    monkeypatch.setattr(
        error_api.app.state.run_service, "get_run_reference_binding", fail
    )
    assert_declared_response(
        error_api, operation, error_api.request(operation), explicit=True
    )


def test_normal_reads_keep_success_envelopes(error_api):
    for operation in (
        "session_state",
        "list_accounts",
        "get_terminal_email_preference",
        "listWorkflows",
        "getWorkflow",
        "getWorkflowSchema",
        "listCompatibleReferenceProfiles",
        "listRuns",
        "getRun",
        "listRunEvents",
        "listRunLogs",
        "listRunArtifacts",
        "listRunQcMetrics",
        "listArtifactPublications",
    ):
        response = error_api.request(operation)
        assert response.status_code == 200, (operation, response.text)
        body = response.json()
        assert body.get("ok", True) is True
        schema = error_api.operations[operation][2]["responses"]["200"]["content"][
            "application/json"
        ]["schema"]
        Draft202012Validator(
            {**schema, "components": error_api.spec["components"]}
        ).validate(body)


def test_every_exported_operation_has_an_exercised_internal_error(error_api):
    assert set(error_api.operations) == set(FAULTS)
