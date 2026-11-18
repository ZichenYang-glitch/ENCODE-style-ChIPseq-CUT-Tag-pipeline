"""Temporary SQLite API for the catalog page browser path.

Slice A is the parent commit. This process seeds one project through the real
catalog HTTP API, then serves that API until it is stopped.
"""

from __future__ import annotations

import argparse
from http.cookiejar import CookieJar
import json
from pathlib import Path
import socket
import tempfile
import threading
import time
import urllib.error
import urllib.request

from encode_pipeline.api.main import create_app
from encode_pipeline.services.authentication_service import (
    AccountAdministrationService,
    AuthenticationActor,
)
from encode_pipeline.services.data_registry import DataRegistryService


ADMIN_USERNAME = "catalog-admin"
ADMIN_PASSWORD = "catalog-browser-admin"
MEMBER_USERNAME = "catalog-member"
MEMBER_PASSWORD = "catalog-browser-member"


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


def _seed(base: str) -> dict[str, object]:
    jar = CookieJar()
    # Ignore HTTP_PROXY. Python's bypass does not treat 127.0.0.1 as loopback
    # for a 127.* no_proxy pattern, so the seed client would hit the proxy.
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

    cuttag = term("assay", "CUT&Tag", "assay-cuttag")
    atac = term("assay", "ATAC-seq", "assay-atac")
    target = term("target", "H3K27me3", "target-h3k27me3")
    batch = term("batch", "batch-1", "batch-1")
    role = term("replicate_role", "treatment", "role-treatment")
    cell_type = term("cell_type", "regulatory T cell", "cell-treg")
    return {
        "opener": opener,
        "post": post,
        "cuttag": cuttag,
        "atac": atac,
        "target": target,
        "batch": batch,
        "role": role,
        "cell_type": cell_type,
    }


def _experiments(seed: dict, project_id: str, sample_revision_id: str) -> dict:
    post = seed["post"]

    def create(key: str, name: str, assay: dict) -> dict:
        return post(
            f"/api/v1/catalog/projects/{project_id}/experiments",
            {
                "stable_key": key,
                "display_name": name,
                "assay": _known(assay),
                "target": _known(seed["target"]),
                "batch": _known(seed["batch"]),
                "samples": [
                    {
                        "sample_revision_id": sample_revision_id,
                        "replicate_role": _known(seed["role"]),
                    }
                ],
            },
            f"experiment-{key}",
        )

    first = create("cut-page-00", "Treg CUT&Tag v1", seed["cuttag"])
    revised = post(
        f"/api/v1/catalog/experiments/{first['experiment_id']}/revisions",
        {
            "expected_revision_number": 1,
            "display_name": "Treg CUT&Tag v2",
            "assay": _known(seed["cuttag"]),
            "target": _known(seed["target"]),
            "batch": _known(seed["batch"]),
            "samples": [
                {
                    "sample_revision_id": sample_revision_id,
                    "replicate_role": _known(seed["role"]),
                }
            ],
        },
        "experiment-cut-page-00-revise",
    )
    for number in range(1, 21):
        create(f"cut-page-{number:02d}", f"CUT&Tag page {number:02d}", seed["cuttag"])
    create("atac-excluded", "ATAC excluded", seed["atac"])
    return {"first": first, "revised": revised}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="127.0.0.1")
    args = parser.parse_args()
    root = Path(tempfile.mkdtemp(prefix="catalog-page-"))
    database = root / "catalog.sqlite"
    app = create_app(
        database_url=f"sqlite:///{database}",
        workspace_root=root / "workspaces",
    )
    accounts = AccountAdministrationService(
        repository=app.state.authentication_repository
    )
    admin = accounts.bootstrap_initial_administrator(ADMIN_USERNAME, ADMIN_PASSWORD)
    accounts.create_member(
        AuthenticationActor.local_operator(),
        MEMBER_USERNAME,
        MEMBER_PASSWORD,
    )
    registry = DataRegistryService(
        repository=app.state.persistence.data_registry_repository
    )
    project = registry.create_project("Treg project")
    sample = registry.import_sample(
        project.project_id,
        stable_key="treg",
        display_name="treg",
        attributes={"organism": "human"},
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
    try:
        seeded = _seed(base)
        classified = seeded["post"](
            f"/api/v1/catalog/samples/{sample.sample.sample_id}/catalog-revisions",
            {
                "expected_revision_number": 1,
                "display_name": "treg",
                "catalog": {
                    "schema_version": 1,
                    "species": {"status": "unknown"},
                    "tissue": {"status": "unknown"},
                    "cell_type": {
                        "status": "known",
                        "values": [_pin(seeded["cell_type"])],
                    },
                    "cell_line": {"status": "not_applicable", "values": []},
                    "treatment": {"status": "not_applicable", "values": []},
                    "tag": {"status": "not_applicable", "values": []},
                },
            },
            "catalog-treg",
        )
        created = _experiments(
            seeded,
            project.project_id,
            classified["sample_revision_id"],
        )
    except Exception:
        server.should_exit = True
        raise
    manifest = {
        "api_base": base,
        "admin_username": ADMIN_USERNAME,
        "admin_password": ADMIN_PASSWORD,
        "member_username": MEMBER_USERNAME,
        "member_password": MEMBER_PASSWORD,
        "admin_user_id": admin.user_id,
        "project_id": project.project_id,
        "project_name": "Treg project",
        "assay_term_id": seeded["cuttag"]["term_id"],
        "assay_name": "CUT&Tag",
        "cell_type_term_id": seeded["cell_type"]["term_id"],
        "cell_type_name": "regulatory T cell",
        "sample_revision_id": classified["sample_revision_id"],
        "experiment_id": created["revised"]["experiment_id"],
        "older_revision_id": created["first"]["experiment_revision_id"],
        "current_name": "Treg CUT&Tag v2",
        "older_name": "Treg CUT&Tag v1",
        "page_two_name": "CUT&Tag page 20",
        "excluded_name": "ATAC excluded",
    }
    print(f"CATALOG_PAGE_READY {json.dumps(manifest)}", flush=True)
    try:
        while thread.is_alive():
            thread.join(timeout=3600)
    finally:
        server.should_exit = True


if __name__ == "__main__":
    main()
