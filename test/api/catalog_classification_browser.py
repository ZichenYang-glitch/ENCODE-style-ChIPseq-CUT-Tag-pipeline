"""Opt-in wheel/Chromium acceptance with the production Peri owner and Rust worker.

Set HW_TEST_PERI_NATIVE and HW_TEST_PERI_MODEL_WORKER as for peri_bound_model_browser.
Run node frontend/scripts/catalog-classification-browser.mjs with the checkout's
Python (HELIXWEAVE_PYTHON) and built frontend assets.
"""

from contextlib import closing
import json
from pathlib import Path
import subprocess
import time

import httpx
import pytest
from sqlalchemy.engine import URL

from encode_pipeline.api.main import create_app
from encode_pipeline.deployment.frontend import ManifestStaticApplication
from encode_pipeline.frontend_assets import verify_frontend_api_contract
from encode_pipeline.services.authentication_service import AccountAdministrationService
from encode_pipeline.services.data_registry import DataRegistryService
from encode_pipeline.services.catalog_model_client import load_catalog_model
from encode_pipeline.services.peri_runtime.bound_model import Budget, binding_id
from encode_pipeline.services.peri_runtime.catalog_classify import CatalogTransportError
from peri_bound_model_browser import serve_on_loopback, start_bound_owner
from peri_run_actions_browser import _build_wheel_assets
from test_experiment_catalog_http import _known


ADMIN = {"username": "classification-admin", "password": "browser-fixture-admin"}
TARGET = "H3K4me3"


@pytest.fixture(scope="module")
def wheel_assets(repo_root, tmp_path_factory):
    return _build_wheel_assets(repo_root, tmp_path_factory.mktemp("catalog-wheel"))


def seed_catalog(app, client):
    accounts = AccountAdministrationService(
        repository=app.state.authentication_repository
    )
    accounts.bootstrap_initial_administrator(ADMIN["username"], ADMIN["password"])
    client.post("/api/v1/auth/login", json=ADMIN).raise_for_status()
    client.headers["X-CSRF-Token"] = client.cookies[
        app.state.auth_cookie_policy.csrf_cookie.name
    ]
    registry = DataRegistryService(
        repository=app.state.persistence.data_registry_repository
    )
    project_id = registry.create_project("Synthetic classification").project_id
    sample = registry.import_sample(
        project_id, stable_key="fixture", display_name="Synthetic sample", attributes={}
    )

    def post(path, body, key):
        response = client.post(path, json=body, headers={"Idempotency-Key": key})
        response.raise_for_status()
        return response.json()

    def term(kind, label):
        return post(
            "/api/v1/catalog/terms", {"term_type": kind, "display_name": label}, kind
        )

    assay = _known(term("assay", "CUT&Tag"))
    role = _known(term("replicate_role", "biological replicate"))
    target = _known(term("target", TARGET))

    def experiment(key):
        return post(
            f"/api/v1/catalog/projects/{project_id}/experiments",
            {
                "stable_key": key,
                "display_name": f"Synthetic {key} experiment",
                "assay": assay,
                "target": {"status": "unknown"},
                "batch": {"status": "unknown"},
                "samples": [
                    {
                        "sample_revision_id": sample.revision.sample_revision_id,
                        "replicate_role": role,
                    }
                ],
            },
            key,
        )

    return {
        "project": project_id,
        "target": target,
        "label": TARGET,
        "selected": experiment("selected"),
        "other": experiment("other"),
    }


@pytest.mark.parametrize("locale", ["en-US", "zh-CN"])
def test_packaged_catalog_classification(
    locale, tmp_path, repo_root, request, monkeypatch, wheel_assets
):
    database = URL.create("sqlite", database=str(tmp_path / "catalog.sqlite"))
    env, model, _, calls = start_bound_owner(
        tmp_path, request, database, model_name="fixture-target", reply="h3k4me3"
    )
    # The API needs only the socket, even on a host with an operator profile.
    for name in ("LEDGER", "SETTINGS", "WORKER"):
        monkeypatch.delenv(f"HELIXWEAVE_PERI_MODEL_{name}", raising=False)
    monkeypatch.setenv(
        "HELIXWEAVE_PERI_OWNER_SOCKET", env["HELIXWEAVE_PERI_OWNER_SOCKET"]
    )
    app = create_app(database_url=str(database), workspace_root=tmp_path / "workspaces")
    request.addfinalizer(app.state.persistence.close)
    request.addfinalizer(app.state.run_queue.close)
    verify_frontend_api_contract(wheel_assets.manifest, app.openapi())

    @app.get("/__fixture/provider-count", include_in_schema=False)
    def provider_count():
        return {"count": len(calls)}

    app.mount("/", ManifestStaticApplication(wheel_assets), name="frontend")
    with (
        serve_on_loopback(app) as origin,
        httpx.Client(base_url=origin, trust_env=False) as client,
    ):
        seed = seed_catalog(app, client)

        def ask_other(key):
            experiment = seed["other"]
            return client.post(
                f"/api/v1/catalog/projects/{seed['project']}/experiments/"
                f"{experiment['experiment_id']}/classification-suggestions",
                headers={"Idempotency-Key": key},
                json={
                    "field": "target",
                    "reason_code": "TARGET_TERM_MATCH",
                    "experiment_revision_id": experiment["experiment_revision_id"],
                    "clue": TARGET,
                    "vocabulary_term_revision_id": seed["target"]["term_revision_id"],
                },
            )

        unauthorized = ask_other("before-grant")
        assert unauthorized.status_code == 409
        assert unauthorized.json()["issues"][0]["code"] == "catalog_model_unauthorized"
        assert len(calls) == 0
        grant_id = model.catalog_grant(
            project_id=seed["project"],
            field="experiment.target",
            binding=binding_id(model.worker({"op": "binding"})),
            expires=int(time.time()) + 300,
            budget=Budget(
                requests=1,
                input_tokens=8192,
                output_tokens=128,
                minor_units=10,
                output_per_call=64,
                input_rate=100,
                output_rate=100,
                currency="CNY",
            ),
        )
        finished = subprocess.run(
            ["node", "scripts/catalog-classification-browser.mjs", "--fixture-stdin"],
            cwd=Path(repo_root) / "frontend",
            text=True,
            capture_output=True,
            timeout=100,
            input=json.dumps(
                {"origin": origin, "admin": ADMIN, "locale": locale, **seed}
            ),
        )
        assert finished.returncode == 0, finished.stdout + finished.stderr
        model.catalog_revoke(grant_id)
        revoked = ask_other("after-revocation")
        assert revoked.status_code == 409
        assert revoked.json()["issues"][0]["code"] == "catalog_model_revoked"

    assert len(calls) == 1
    method, path, raw = calls[0]
    assert (method, path) == ("POST", "/v1/messages")
    body = json.loads(raw)
    assert body["model"] == "fixture-target"
    assert body["max_tokens"] == 64 and not body.get("tools")
    (message,) = body["messages"]
    assert message["role"] == "user"
    content = message["content"]
    if isinstance(content, list):
        (block,) = content
        assert block["type"] == "text"
        content = block["text"]
    projection = json.loads(content[content.index("{") :])
    assert projection == json.loads(finished.stdout)["projection"]
    assert projection == {
        "field": "experiment.target",
        "reason_code": "TARGET_TERM_MATCH",
        "clue": TARGET,
        "candidates": [{"label": TARGET, "token": "h3k4me3"}],
    }
    wire = json.dumps(body)
    for private in (
        seed["project"],
        seed["target"]["term_revision_id"],
        seed["selected"]["experiment_revision_id"],
        seed["other"]["experiment_id"],
        ADMIN["password"],
    ):
        assert private not in wire
    with closing(model.db()) as ledger:
        assert (
            ledger.execute(
                "SELECT requests FROM catalog_grants WHERE grant_id=?", (grant_id,)
            ).fetchone()[0]
            == 0
        )
        assert [
            row[0]
            for row in ledger.execute(
                "SELECT delivered FROM catalog_reservations WHERE grant_id=?",
                (grant_id,),
            )
        ] == [1]
    print(finished.stdout.strip())


@pytest.mark.parametrize(
    ("reply", "error"),
    [
        ("insufficient_evidence", None),
        ("1", "catalog_model_malformed"),
        pytest.param(None, "catalog_model_timeout", id="provider-timeout"),
    ],
)
def test_native_classification_decisions_and_reservations(
    tmp_path, request, reply, error
):
    database = URL.create("sqlite", database=str(tmp_path / "catalog.sqlite"))
    env, model, _, calls = start_bound_owner(
        tmp_path, request, database, model_name="fixture-target", reply=reply
    )
    transport = load_catalog_model(env)
    project = "prj_" + "1" * 32
    grant = model.catalog_grant(
        project_id=project,
        field="experiment.target",
        binding=binding_id(model.worker({"op": "binding"})),
        expires=int(time.time()) + 300,
        budget=Budget(1, 4096, 64, 10, 64, 100, 100, "CNY"),
    )
    fields = dict(
        project_id=project,
        clue=TARGET,
        candidates=(("h3k4me3", TARGET),),
        idempotency_key="decision",
        projection_digest="a" * 64,
    )
    with pytest.raises(CatalogTransportError) as foreign:
        transport.classify(**(fields | {"project_id": "prj_" + "2" * 32}))
    assert foreign.value.code == "catalog_model_project_mismatch"
    assert calls == []
    if error is None:
        assert transport.classify(**fields) == reply
        assert transport.classify(**fields) == reply
    else:
        with pytest.raises(CatalogTransportError) as rejected:
            transport.classify(**fields)
        assert rejected.value.code == error
        with pytest.raises(CatalogTransportError) as replay:
            transport.classify(**fields)
        assert replay.value.code == "catalog_model_timeout"
    with closing(model.db()) as ledger, ledger:
        balance = ledger.execute(
            "SELECT requests, input_left, output_left, minor_left FROM catalog_grants WHERE grant_id=?",
            (grant,),
        ).fetchone()
        assert tuple(balance) == (0, 0, 0, 9)
        ledger.execute(
            "UPDATE catalog_grants SET expires=? WHERE grant_id=?",
            (int(time.time()) - 1, grant),
        )
        reservation = ledger.execute(
            "SELECT delivered, decision FROM catalog_reservations WHERE grant_id=?",
            (grant,),
        ).fetchone()
        assert tuple(reservation) == (
            (1, reply) if reply == "insufficient_evidence" else (0, None)
        )
    with pytest.raises(CatalogTransportError) as expired:
        transport.classify(**(fields | {"idempotency_key": "after-expiry"}))
    assert expired.value.code == "catalog_model_expired"
    assert len(calls) == 1
