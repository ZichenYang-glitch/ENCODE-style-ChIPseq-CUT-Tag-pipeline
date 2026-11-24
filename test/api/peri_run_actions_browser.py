"""Opt-in browser acceptance against the real app and an isolated Peri fixture."""

from __future__ import annotations

import json
from pathlib import Path
import socket
import subprocess
import sys
from threading import Thread
from time import monotonic, sleep
from zipfile import ZipFile

from sqlalchemy import select
from sqlalchemy.orm import Session
import uvicorn

from encode_pipeline.deployment.frontend import ManifestStaticApplication
from encode_pipeline.frontend_assets import (
    VerifiedFrontendAssets,
    parse_manifest_bytes,
    verify_frontend_api_contract,
)
from encode_pipeline.persistence.models import RunExecutionAssignmentRow, RunRow
from test_peri_run_actions import (
    WORKFLOW,
    acknowledge_stop,
    peri_case,  # noqa: F401 - pytest discovers the shared fixture here
    seed_run,
)


def test_peri_run_actions_from_wheel_in_chromium(peri_case, repo_root, tmp_path):  # noqa: F811
    client, sid, project_id, _, stop = peri_case
    app = client.app
    run_a, assignment = seed_run(client, running=True)
    run_b, _ = seed_run(client)
    frontend = Path(repo_root) / "frontend"
    assets = _build_wheel_assets(repo_root, tmp_path)

    @app.post("/fixture/acknowledge-stop", include_in_schema=False)
    def acknowledge():
        acknowledge_stop(client, assignment)

    verify_frontend_api_contract(assets.manifest, app.openapi())
    app.mount("/", ManifestStaticApplication(assets), name="frontend")

    policy = app.state.auth_cookie_policy
    session_token, csrf_token = app.state.test_auth_tokens
    with socket.socket() as listener:
        listener.bind(("127.0.0.1", 0))
        listener.listen()
        url = f"http://127.0.0.1:{listener.getsockname()[1]}"
        server = uvicorn.Server(
            uvicorn.Config(app, log_level="warning", loop="asyncio")
        )
        thread = Thread(target=server.run, kwargs={"sockets": [listener]})
        thread.start()
        try:
            deadline = monotonic() + 10
            while not server.started and thread.is_alive() and monotonic() < deadline:
                sleep(0.02)
            assert server.started, "Acceptance HTTP server did not start"
            result = subprocess.run(
                ["node", str(frontend / "e2e/peri-run-actions.acceptance.mjs")],
                cwd=frontend,
                input=json.dumps(
                    {
                        "url": url,
                        "workflow": WORKFLOW,
                        "sid": sid,
                        "projectId": project_id,
                        "runA": run_a,
                        "runB": run_b,
                        "session": [policy.session_cookie.name, session_token],
                        "csrf": [policy.csrf_cookie.name, csrf_token],
                    }
                ),
                text=True,
                capture_output=True,
                timeout=65,
                check=False,
            )
            assert result.returncode == 0, result.stdout + result.stderr
            assert stop.run_ids == [run_a]
            with Session(app.state.persistence.engine) as db:
                assert len(db.scalars(select(RunRow)).all()) == 3
                assert db.scalars(select(RunExecutionAssignmentRow.run_id)).all() == [
                    run_a
                ]
        finally:
            server.should_exit = True
            thread.join(timeout=10)
            assert not thread.is_alive(), "Acceptance HTTP server did not stop"


def _build_wheel_assets(repo_root, tmp_path) -> VerifiedFrontendAssets:
    output = tmp_path / "wheel"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "build",
            "--wheel",
            "--no-isolation",
            "--outdir",
            str(output),
        ],
        cwd=repo_root,
        text=True,
        capture_output=True,
        timeout=60,
        check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    (wheel,) = output.glob("*.whl")
    prefix = "encode_pipeline/frontend_assets/"
    with ZipFile(wheel) as package:
        manifest = parse_manifest_bytes(package.read(prefix + "asset-manifest.json"))
        static = prefix + "static/"
        content = {
            item.filename.removeprefix(static): package.read(item)
            for item in package.infolist()
            if item.filename.startswith(static) and not item.is_dir()
        }
    return VerifiedFrontendAssets.create(manifest, content)
