"""Temporary catalog and loopback model for the packaged classification page.

The process owns one private SQLite file. It seeds a synthetic administrator,
controlled vocabulary, and two experiments whose target is unknown. A loopback
Anthropic provider counts every model POST and never uses a paid endpoint.

This does not register an offline proposal. Offline registration is not a
model source. The browser script drives Ask Peri, Preview selected, and
Adopt previewed changes when the packaged page exposes them.
"""

from __future__ import annotations

import argparse
from http.cookiejar import CookieJar
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
from pathlib import Path
import socket
import sys
import tempfile
import threading
import time
import urllib.error
import urllib.request

from encode_pipeline.api.main import create_app
from encode_pipeline.deployment.frontend import ManifestStaticApplication
from encode_pipeline.frontend_assets import load_packaged_frontend_assets
from encode_pipeline.services.authentication_service import (
    AccountAdministrationService,
    AuthenticationActor,
)
from encode_pipeline.services.data_registry import DataRegistryService


ADMIN_USERNAME = "classification-admin"
ADMIN_PASSWORD = "classification-browser-admin"
SELECTED_NAME = "Selected target candidate"
EXCLUDED_NAME = "Excluded target candidate"


def _free_port() -> int:
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return int(sock.getsockname()[1])


def _csrf(jar: CookieJar) -> str:
    for cookie in jar:
        if cookie.name == "helixweave_csrf":
            return cookie.value
    raise RuntimeError("login did not set a CSRF cookie")


def _request(
    opener: urllib.request.OpenerDirector,
    method: str,
    url: str,
    body: dict | None = None,
    headers: dict[str, str] | None = None,
) -> dict:
    data = None if body is None else json.dumps(body).encode()
    request = urllib.request.Request(url, data=data, method=method)
    if body is not None:
        request.add_header("Content-Type", "application/json")
    for key, value in (headers or {}).items():
        request.add_header(key, value)
    try:
        with opener.open(request, timeout=30) as response:
            raw = response.read()
            status = response.status
    except urllib.error.HTTPError as error:
        raw = error.read()
        status = error.code
    payload = json.loads(raw.decode()) if raw else {}
    if status >= 400:
        raise RuntimeError(f"{method} {url} -> {status} {payload}")
    if not isinstance(payload, dict):
        raise RuntimeError(f"{method} {url} returned {payload!r}")
    return payload


def _wait_until_ready(opener: urllib.request.OpenerDirector, base: str) -> None:
    deadline = time.monotonic() + 30
    last_error = "API did not start"
    while time.monotonic() < deadline:
        try:
            _request(opener, "GET", f"{base}/api/v1/auth/session")
            return
        except (RuntimeError, urllib.error.URLError, TimeoutError) as error:
            last_error = str(error)
            time.sleep(0.1)
    raise RuntimeError(last_error)


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


class LoopbackProvider:
    """Count model posts. The event stream matches peri_bound_model_browser."""

    def __init__(self, suggestion: dict) -> None:
        self.suggestion = suggestion
        self.calls: list[dict] = []
        parent = self

        class Handler(BaseHTTPRequestHandler):
            def do_GET(self) -> None:  # noqa: N802
                if self.path != "/__fixture/calls":
                    self.send_response(404)
                    self.end_headers()
                    return
                body = json.dumps({"count": len(parent.calls)}).encode()
                self.send_response(200)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

            def do_POST(self) -> None:  # noqa: N802
                length = int(self.headers.get("Content-Length", "0"))
                raw = self.rfile.read(length)
                try:
                    posted = json.loads(raw.decode() or "{}")
                except json.JSONDecodeError:
                    posted = {"unparsed": True}
                parent.calls.append({"path": self.path, "body": posted})
                self.send_response(200)
                self.send_header("Content-Type", "text/event-stream")
                self.end_headers()
                text = json.dumps(parent.suggestion, separators=(",", ":"))
                events = [
                    (
                        "message_start",
                        {"message": {"id": "fixture", "usage": {"input_tokens": 12}}},
                    ),
                    (
                        "content_block_start",
                        {"index": 0, "content_block": {"type": "text", "text": ""}},
                    ),
                    (
                        "content_block_delta",
                        {"index": 0, "delta": {"type": "text_delta", "text": text}},
                    ),
                    ("content_block_stop", {"index": 0}),
                    (
                        "message_delta",
                        {
                            "delta": {"stop_reason": "end_turn"},
                            "usage": {"output_tokens": 32},
                        },
                    ),
                    ("message_stop", {}),
                ]
                for event, data in events:
                    self.wfile.write(
                        f"event: {event}\ndata: {json.dumps(data)}\n\n".encode()
                    )

            def log_message(self, format: str, *args) -> None:
                return

        self._server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        self._thread = threading.Thread(target=self._server.serve_forever, daemon=True)

    @property
    def base(self) -> str:
        port = self._server.server_address[1]
        return f"http://127.0.0.1:{port}"

    def start(self) -> None:
        self._thread.start()

    def close(self) -> None:
        self._server.shutdown()
        self._thread.join(timeout=5)
        self._server.server_close()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="127.0.0.1")
    args = parser.parse_args()
    root = Path(tempfile.mkdtemp(prefix="catalog-classification-"))
    database = root / "catalog.sqlite"
    if not database.is_relative_to(root):
        raise RuntimeError("database escaped the temporary directory")
    app = create_app(
        database_url=f"sqlite:///{database}",
        workspace_root=root / "workspaces",
    )
    app.mount(
        "/",
        ManifestStaticApplication(load_packaged_frontend_assets()),
        name="frontend",
    )
    accounts = AccountAdministrationService(
        repository=app.state.authentication_repository
    )
    accounts.bootstrap_initial_administrator(ADMIN_USERNAME, ADMIN_PASSWORD)
    registry = DataRegistryService(
        repository=app.state.persistence.data_registry_repository
    )
    project = registry.create_project("Treg classification")
    sample = registry.import_sample(
        project.project_id,
        stable_key="treg",
        display_name="treg",
        attributes={"organism": "synthetic"},
    )
    port = _free_port()
    import uvicorn

    server = uvicorn.Server(
        uvicorn.Config(
            app,
            host=args.host,
            port=port,
            log_level="warning",
            access_log=False,
        )
    )
    server.install_signal_handlers = lambda: None
    thread = threading.Thread(target=server.run, daemon=True)
    thread.start()
    base = f"http://{args.host}:{port}"
    provider: LoopbackProvider | None = None
    try:
        jar = CookieJar()
        opener = urllib.request.build_opener(
            urllib.request.ProxyHandler({}),
            urllib.request.HTTPCookieProcessor(jar),
        )
        _wait_until_ready(opener, base)
        _request(
            opener,
            "POST",
            f"{base}/api/v1/auth/login",
            {"username": ADMIN_USERNAME, "password": ADMIN_PASSWORD},
        )

        def post(path: str, body: dict, key: str) -> dict:
            return _request(
                opener,
                "POST",
                f"{base}{path}",
                body,
                {"Idempotency-Key": key, "X-CSRF-Token": _csrf(jar)},
            )

        def term(term_type: str, name: str, key: str) -> dict:
            return post(
                "/api/v1/catalog/terms",
                {"term_type": term_type, "display_name": name},
                key,
            )

        assay = term("assay", "CUT&Tag", "assay-cuttag")
        target = term("target", "H3K4me3", "target-h3k4me3")
        role = term("replicate_role", "treatment", "role-treatment")
        cell_type = term("cell_type", "regulatory T cell", "cell-treg")
        classified = post(
            f"/api/v1/catalog/samples/{sample.sample.sample_id}/catalog-revisions",
            {
                "expected_revision_number": 1,
                "display_name": "treg",
                "catalog": {
                    "schema_version": 1,
                    "species": {"status": "unknown"},
                    "tissue": {"status": "unknown"},
                    "cell_type": {"status": "known", "values": [_pin(cell_type)]},
                    "cell_line": {"status": "not_applicable", "values": []},
                    "treatment": {"status": "not_applicable", "values": []},
                    "tag": {"status": "not_applicable", "values": []},
                },
            },
            "catalog-treg",
        )
        seeded = {
            "post": post,
            "assay": assay,
            "target": target,
            "role": role,
            "cell_type": cell_type,
        }
        sample_revision_id = classified["sample_revision_id"]

        def experiment(key: str, name: str) -> dict:
            created = seeded["post"](
                f"/api/v1/catalog/projects/{project.project_id}/experiments",
                {
                    "stable_key": key,
                    "display_name": name,
                    "assay": _known(assay),
                    "target": {"status": "unknown"},
                    "batch": {"status": "unknown"},
                    "samples": [
                        {
                            "sample_revision_id": sample_revision_id,
                            "replicate_role": _known(role),
                        }
                    ],
                },
                f"experiment-{key}",
            )
            if created.get("revision_number") != 1:
                raise RuntimeError(f"{name} did not start at revision 1")
            return {
                "experiment_id": created["experiment_id"],
                "experiment_revision_id": created["experiment_revision_id"],
                "revision_number": created["revision_number"],
                "display_name": name,
            }

        selected = experiment("selected-target", SELECTED_NAME)
        excluded = experiment("excluded-target", EXCLUDED_NAME)
        suggestion = {
            "outcome": "suggestion",
            "field": "experiment.target",
            "reason_code": "TARGET_TERM_MATCH",
            "term_id": target["term_id"],
            "term_revision_id": target["revision"]["term_revision_id"],
            "display_name": "H3K4me3",
            "items": [
                {
                    "experiment_id": item["experiment_id"],
                    "base_revision_id": item["experiment_revision_id"],
                    "base_revision_number": item["revision_number"],
                    "evidence_token": "vocab-h3k4me3",
                }
                for item in (selected, excluded)
            ],
        }
        provider = LoopbackProvider(suggestion)
        provider.start()
        settings = {
            "config": {
                "active_alias": "fable",
                "profiles": {"fable": {"provider": "fixture"}},
                "providers": [
                    {
                        "id": "fixture",
                        "type": "anthropic",
                        "apiKey": "fixture-only",
                        "models": {"fable": "fixture-classification"},
                        "baseUrl": provider.base,
                    }
                ],
            }
        }
        # Same shape as peri_bound_model_browser. The API process does not
        # load this file: a partial Peri owner configuration fails closed.
        settings_path = root / "peri-settings.json"
        settings_path.write_text(json.dumps(settings))
        if not database.is_file() or not database.is_relative_to(root):
            raise RuntimeError("synthetic database was not created in the temp root")
        manifest = {
            "api_base": base,
            "admin_username": ADMIN_USERNAME,
            "admin_password": ADMIN_PASSWORD,
            "project_id": project.project_id,
            "project_name": "Treg classification",
            "assay_term_id": assay["term_id"],
            "cell_type_term_id": cell_type["term_id"],
            "target_term_id": target["term_id"],
            "target_term_revision_id": target["revision"]["term_revision_id"],
            "target_name": "H3K4me3",
            "selected": selected,
            "excluded": excluded,
            "provider_base": provider.base,
            "provider_calls": f"{provider.base}/__fixture/calls",
            "database_is_temporary": True,
            "offline_proposal_registered": False,
            "worker_started": False,
            "suggestion": suggestion,
        }
    except Exception:
        server.should_exit = True
        if provider is not None:
            provider.close()
        raise
    print(f"CLASSIFICATION_FIXTURE_READY {json.dumps(manifest)}", flush=True)
    try:
        while thread.is_alive():
            thread.join(timeout=3600)
    finally:
        server.should_exit = True
        if provider is not None:
            provider.close()


if __name__ == "__main__":
    main()
