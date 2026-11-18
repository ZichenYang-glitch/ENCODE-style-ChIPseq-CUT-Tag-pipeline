"""Exercise the catalog CLI with authenticated requests over loopback HTTP."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import replace
from datetime import datetime, timedelta, timezone
import json
import socket
from threading import Event, Thread

import pytest

from encode_pipeline.cli import app, catalog


_PREFIX = "/api/v1/catalog"
_PROJECT = "prj_" + "1" * 32
_EXPERIMENT = "exp_" + "1" * 32
_REVISION = "expr_" + "1" * 32
_ASSAY = "trmv_" + "1" * 32
_TOKEN = "t" * 43
_SESSION_ENV = "HELIXWEAVE_CATALOG_SESSION_TOKEN"


def _preview_args(tmp_path, origin, project, rows):
    source = tmp_path / "classifications.tsv"
    source.write_text(
        "experiment_id\texpected_experiment_revision_id\tassay_term_revision_id\n"
        + "".join("\t".join(row) + "\n" for row in rows),
        encoding="utf-8",
    )
    return [
        "catalog",
        "preview",
        "--api-origin",
        origin,
        "--project-id",
        project,
        str(source),
    ]


@contextmanager
def _serve(application):
    uvicorn = pytest.importorskip("uvicorn")
    ready = Event()

    class Server(uvicorn.Server):
        async def startup(self, sockets=None):
            await super().startup(sockets=sockets)
            if self.started:
                ready.set()

    with socket.socket() as listener:
        listener.bind(("127.0.0.1", 0))
        port = listener.getsockname()[1]
        server = Server(uvicorn.Config(application, log_level="error", http="h11"))
        thread = Thread(target=server.run, kwargs={"sockets": [listener]}, daemon=True)
        thread.start()
        try:
            assert ready.wait(timeout=10), "catalog HTTP test server failed to start"
            yield f"http://127.0.0.1:{port}"
        finally:
            server.should_exit = True
            thread.join(timeout=10)
            assert not thread.is_alive(), "catalog HTTP test server failed to stop"


def test_preview_matches_catalog_permissions_and_revision_history(
    tmp_path, monkeypatch, capsys
):
    httpx = pytest.importorskip("httpx")
    pytest.importorskip("fastapi")
    from conftest import (
        TEST_AUTH_CSRF_TOKEN,
        TEST_AUTH_SESSION_TOKEN,
        seed_test_authentication,
    )
    from encode_pipeline.api.main import create_app
    from encode_pipeline.platform.authentication import UserRole
    from encode_pipeline.services.authentication import (
        new_session_record,
        new_session_secrets,
    )
    from encode_pipeline.services.data_registry import DataRegistryService

    monkeypatch.setenv("ENCODE_PIPELINE_AUTH_SECURE_COOKIES", "0")
    monkeypatch.setattr(catalog, "_PAGE_SIZE", 1)
    application = create_app(
        database_url=f"sqlite:///{tmp_path / 'catalog.db'}",
        workspace_root=tmp_path / "workspaces",
    )
    admin = seed_test_authentication(application)
    member = replace(
        admin,
        user_id="usr_" + "1" * 32,
        username="preview-member",
        role=UserRole.MEMBER,
    )
    repository = application.state.authentication_repository
    repository.create_account(member)
    member_secrets = new_session_secrets()
    repository.create_session(
        new_session_record(
            user_id=member.user_id,
            secrets=member_secrets,
            created_at=datetime.now(timezone.utc),
            lifetime=timedelta(minutes=10),
        )
    )
    registry = DataRegistryService(
        repository=application.state.persistence.data_registry_repository
    )
    project = registry.create_project("Preview fixture").project_id
    other_project = registry.create_project("Other project").project_id

    def sample_revision(project_id):
        return registry.import_sample(
            project_id, stable_key="sample", display_name="Sample", attributes={}
        ).revision.sample_revision_id

    sample = sample_revision(project)
    other_sample = sample_revision(other_project)
    requests = []

    async def observed(scope, receive, send):
        if scope["type"] == "http":
            requests.append((scope["method"], scope["path"], scope["query_string"]))
        await application(scope, receive, send)

    with (
        _serve(observed) as origin,
        httpx.Client(base_url=origin, trust_env=False) as client,
    ):
        client.cookies.set("helixweave_session", TEST_AUTH_SESSION_TOKEN)
        client.cookies.set("helixweave_csrf", TEST_AUTH_CSRF_TOKEN)
        client.headers["X-CSRF-Token"] = TEST_AUTH_CSRF_TOKEN

        def post(path, key, payload):
            response = client.post(
                _PREFIX + path, json=payload, headers={"Idempotency-Key": key}
            )
            assert response.status_code == 200, response.text
            return response.json()

        def assay_pin(name):
            term = post(
                "/terms",
                "assay-" + name.lower().replace(" ", "-"),
                {"term_type": "assay", "display_name": name},
            )
            return {
                "status": "known",
                "term_id": term["term_id"],
                "term_revision_id": term["revision"]["term_revision_id"],
            }

        first_assay = assay_pin("First assay")
        second_assay = assay_pin("Second assay")

        def classification(name, sample_id, assay):
            return {
                "display_name": name,
                "assay": assay,
                "target": {"status": "unknown"},
                "batch": {"status": "unknown"},
                "samples": [
                    {
                        "sample_revision_id": sample_id,
                        "replicate_role": {"status": "unknown"},
                    }
                ],
            }

        def experiment(key, project_id=project, sample_id=sample, assay=first_assay):
            return post(
                f"/projects/{project_id}/experiments",
                key,
                {"stable_key": key, **classification(key, sample_id, assay)},
            )

        revised = experiment("revised")
        revised_id = revised["experiment_id"]
        history_path = f"/experiments/{revised_id}/revisions"
        for number, assay in enumerate((first_assay, second_assay), start=1):
            revised = post(
                history_path,
                f"revise-{number}",
                {
                    "expected_revision_number": number,
                    **classification(f"Revision {number + 1}", sample, assay),
                },
            )
            if number == 1:
                old_revision = revised["experiment_revision_id"]
        matched = experiment(
            "matched",
            assay={
                "status": "ambiguous",
                "candidates": [
                    {key: value for key, value in pin.items() if key != "status"}
                    for pin in (first_assay, second_assay)
                ],
            },
        )
        duplicate = experiment("duplicate")
        wrong_revision = experiment("wrong-revision")
        foreign = experiment("foreign", other_project, other_sample)
        proposed = first_assay["term_revision_id"]

        def row(record, expected=None):
            return (
                record["experiment_id"],
                expected or record["experiment_revision_id"],
                proposed,
            )

        args = _preview_args(
            tmp_path,
            origin,
            project,
            [
                row(matched),
                row(duplicate),
                row(revised, old_revision),
                row(foreign),
                row(wrong_revision, foreign["experiment_revision_id"]),
                row(duplicate),
            ],
        )
        detail_path = _PREFIX + history_path + f"/{old_revision}"
        assert client.get(detail_path).status_code == 200
        client.cookies.set("helixweave_session", member_secrets.session_token)
        assert client.get(detail_path).status_code == 403

        requests.clear()
        monkeypatch.delenv(_SESSION_ENV, raising=False)
        assert app.main(args) == 1
        anonymous = capsys.readouterr()
        assert not anonymous.out and "HTTP 401" in anonymous.err
        assert requests == [
            (
                "GET",
                f"{_PREFIX}/projects/{project}/experiments",
                b"revision_view=linked&limit=1",
            )
        ]

        monkeypatch.setenv(_SESSION_ENV, TEST_AUTH_SESSION_TOKEN)
        assert app.main(args) == 0
        administrator = capsys.readouterr()
        requests.clear()
        monkeypatch.setenv(_SESSION_ENV, member_secrets.session_token)
        assert app.main(args) == 0
        member_output = capsys.readouterr()

    assert not administrator.err and not member_output.err
    assert administrator.out == member_output.out
    assert member_secrets.session_token not in member_output.out
    assert TEST_AUTH_SESSION_TOKEN not in administrator.out
    assert "candidates" not in member_output.out
    report = json.loads(member_output.out)
    assert report["project_id"] == project
    rows = report["rows"]
    assert [entry["line"] for entry in rows] == list(range(2, 8))
    assert [entry["status"] for entry in rows] == [
        "matched",
        "duplicate",
        "stale_revision",
        "unmatched",
        "unknown_revision",
        "duplicate",
    ]
    assert rows[0]["current_assay_term_revision_id"] is None
    assert rows[2] == {
        "line": 4,
        "experiment_id": revised_id,
        "expected_experiment_revision_id": old_revision,
        "assay_term_revision_id": proposed,
        "current_experiment_revision_id": revised["experiment_revision_id"],
        "current_assay_term_revision_id": second_assay["term_revision_id"],
        "status": "stale_revision",
    }
    expected_paths = {
        f"{_PREFIX}/projects/{project}/experiments",
        _PREFIX + history_path,
        f"{_PREFIX}/experiments/{wrong_revision['experiment_id']}/revisions",
    }
    assert {method for method, _, _ in requests} == {"GET"}
    assert {path for _, path, _ in requests} == expected_paths
    assert {path for _, path, query in requests if b"cursor=" in query} == {
        f"{_PREFIX}/projects/{project}/experiments",
        _PREFIX + history_path,
    }


@pytest.mark.parametrize(
    ("mode", "diagnostic"),
    [
        ("redirect", "HTTP 302"),
        ("http-error", "HTTP 500"),
        ("repeated-cursor", "repeated catalog cursor"),
    ],
)
def test_preview_http_failure_does_not_emit_partial_report_or_secrets(
    mode, diagnostic, tmp_path, monkeypatch, capsys
):
    pytest.importorskip("fastapi")
    from fastapi import FastAPI
    from fastapi.responses import JSONResponse

    application = FastAPI()
    requests = []

    @application.get(f"{_PREFIX}/projects/{_PROJECT}/experiments")
    def list_experiments(cursor: str | None = None):
        requests.append(cursor)
        if mode == "redirect":
            return JSONResponse(
                {"private": _TOKEN},
                status_code=302,
                headers={"Location": "/redirected"},
            )
        if mode == "http-error":
            return JSONResponse({"private": _TOKEN}, status_code=500)
        return {
            "ok": True,
            "project_id": _PROJECT,
            "revision_view": "linked",
            "experiments": [
                {
                    "experiment_id": _EXPERIMENT,
                    "experiment_revision_id": _REVISION,
                    "assay": {"status": "unknown"},
                }
            ],
            "next_cursor": "repeated",
        }

    @application.get("/redirected")
    def redirected():
        requests.append("redirected")

    monkeypatch.setenv(_SESSION_ENV, _TOKEN)
    with _serve(application) as origin:
        args = _preview_args(
            tmp_path, origin, _PROJECT, [(_EXPERIMENT, _REVISION, _ASSAY)]
        )
        assert app.main(args) == 1
    output = capsys.readouterr()
    assert not output.out and diagnostic in output.err and _TOKEN not in output.err
    assert requests == ([None, "repeated"] if mode == "repeated-cursor" else [None])


@pytest.mark.parametrize(
    ("case", "diagnostic"),
    [
        ("row", "line 2"),
        ("empty", "no classification rows"),
        ("origin", "HTTPS origin"),
        ("lan", "HTTPS origin"),
        ("malformed-loopback", "HTTPS origin"),
        ("token", "session cookie value"),
    ],
)
def test_preview_rejects_invalid_input_without_echoing_it(
    case, diagnostic, tmp_path, monkeypatch, capsys
):
    rows = [(_EXPERIMENT, _REVISION, _ASSAY)]
    origin = "http://127.0.0.1:1"
    monkeypatch.delenv(_SESSION_ENV, raising=False)
    if case == "row":
        rows = [(_TOKEN, _REVISION, _ASSAY)]
    elif case == "empty":
        rows = []
    elif case == "origin":
        origin = f"http://{_TOKEN}@example.test"
    elif case == "lan":
        origin = "http://192.0.2.1"
    elif case == "malformed-loopback":
        origin = "http://[::1]example.test"
    else:
        monkeypatch.setenv(_SESSION_ENV, f"helixweave_session={_TOKEN}")
    args = _preview_args(tmp_path, origin, _PROJECT, rows)
    with pytest.raises(SystemExit, match="2"):
        app.main(args)
    output = capsys.readouterr()
    assert not output.out and diagnostic in output.err and _TOKEN not in output.err


def test_https_uses_the_secure_session_cookie_name():
    assert catalog._api_origin("https://EXAMPLE.test:8443/") == (
        "https://example.test:8443",
        "__Host-helixweave_session",
    )
