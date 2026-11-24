"""Opt-in Chromium + qualified native owner + real Peri model worker; loopback only."""

from contextlib import closing
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import json
import os
from pathlib import Path
import shlex
import socket
import subprocess
import sys
from tempfile import TemporaryDirectory
from threading import Thread
import time

import uvicorn

from encode_pipeline.deployment.frontend import ManifestStaticApplication
from encode_pipeline.services.authoring_owner_client import load_authoring_bridge
from encode_pipeline.services.peri_runtime.bound_model import (
    BoundModel,
    binding_id,
    configured_model,
)
from peri_run_actions_browser import _build_wheel_assets
from test_peri_run_actions import WORKFLOW, peri_case  # noqa: F401


def test_browser_to_bound_worker(peri_case, repo_root, tmp_path, request):  # noqa: F811
    client, sid, project, _, _ = peri_case
    app = client.app
    seen = []

    class Provider(BaseHTTPRequestHandler):
        def do_POST(self):
            body = json.loads(self.rfile.read(int(self.headers["Content-Length"])))
            seen.append((self.path, body))
            self.send_response(200)
            self.send_header("Content-Type", "text/event-stream")
            self.end_headers()
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
                    {"index": 0, "delta": {"type": "text_delta", "text": "7"}},
                ),
                ("content_block_stop", {"index": 0}),
                (
                    "message_delta",
                    {
                        "delta": {"stop_reason": "end_turn"},
                        "usage": {"output_tokens": 1},
                    },
                ),
                ("message_stop", {}),
            ]
            for event, data in events:
                self.wfile.write(
                    f"event: {event}\ndata: {json.dumps(data)}\n\n".encode()
                )

        def log_message(self, *args):
            pass

    provider = ThreadingHTTPServer(("127.0.0.1", 0), Provider)
    pt = Thread(target=provider.serve_forever)
    pt.start()
    request.addfinalizer(
        lambda: (provider.shutdown(), pt.join(), provider.server_close())
    )
    root = TemporaryDirectory(prefix="hw-bound-", dir="/tmp")
    request.addfinalizer(root.cleanup)
    settings = tmp_path / "settings.json"
    settings.write_text(
        json.dumps(
            {
                "config": {
                    "active_alias": "fable",
                    "profiles": {"fable": {"provider": "fixture"}},
                    "providers": [
                        {
                            "id": "fixture",
                            "type": "anthropic",
                            "apiKey": "fixture-only",
                            "models": {"fable": "fixture-number"},
                            "baseUrl": f"http://127.0.0.1:{provider.server_port}",
                        }
                    ],
                }
            }
        )
    )
    worker = Path(os.environ["HW_TEST_PERI_MODEL_WORKER"]).resolve(strict=True)
    shim = tmp_path / "local-worker"
    shim.write_text(
        "#!/bin/sh\nexec " + shlex.quote(str(worker)) + ' "$@" --loopback-test\n'
    )
    shim.chmod(0o700)
    url = app.state.persistence.engine.url
    assert Path(url.database).is_relative_to(tmp_path)
    env = {
        **os.environ,
        "HELIXWEAVE_PERI_BINARY": os.environ["HW_TEST_PERI_NATIVE"],
        "HELIXWEAVE_PERI_OWNER_SOCKET": str(Path(root.name) / "owner.sock"),
        "ENCODE_PIPELINE_DATABASE_URL": str(url),
        "HELIXWEAVE_PERI_MODEL_LEDGER": str(tmp_path / "ledger.db"),
        "HELIXWEAVE_PERI_MODEL_SETTINGS": str(settings),
        "HELIXWEAVE_PERI_MODEL_WORKER": str(shim),
        "TMPDIR": root.name,
    }
    model = configured_model(env)
    BoundModel.initialize(model.ledger)
    proc = subprocess.Popen(
        [sys.executable, "-m", "encode_pipeline.services.peri_runtime.socket_server"],
        env=env,
    )

    def stop():
        if proc.poll() is None:
            proc.terminate()
            proc.wait(timeout=25)

    request.addfinalizer(stop)
    bridge = load_authoring_bridge(env)
    deadline = time.monotonic() + 12
    while (
        not bridge.status().available
        and proc.poll() is None
        and time.monotonic() < deadline
    ):
        time.sleep(0.05)
    assert bridge.status().available
    app.state.authoring_service._bridge = bridge
    app.mount(
        "/",
        ManifestStaticApplication(_build_wheel_assets(repo_root, tmp_path)),
        name="frontend",
    )
    p = app.state.auth_cookie_policy
    a, b = app.state.test_auth_tokens
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        sock.listen()
        srv = uvicorn.Server(uvicorn.Config(app, log_level="warning", loop="asyncio"))
        thread = Thread(target=srv.run, kwargs={"sockets": [sock]})
        thread.start()
        try:
            end = time.monotonic() + 5
            while not srv.started and time.monotonic() < end:
                time.sleep(0.02)
            assert srv.started
            grant = [
                sys.executable,
                "-m",
                "encode_pipeline.services.peri_runtime.model_grant",
                "grant",
                "--platform-db",
                url.database,
                "--sid",
                sid,
                "--field",
                "config.threads",
                "--binding",
                binding_id(model.worker({"op": "binding"})),
                "--expires",
                str(int(time.time()) + 120),
                "--requests",
                "1",
                "--input-tokens",
                "4096",
                "--output-tokens",
                "128",
                "--output-per-call",
                "128",
                "--currency",
                "CNY",
                "--minor-units",
                "1",
                "--input-rate",
                "100",
                "--output-rate",
                "100",
            ]
            result = subprocess.run(
                ["node", "e2e/peri-bound-model.acceptance.mjs"],
                cwd=Path(repo_root) / "frontend",
                env=env,
                input=json.dumps(
                    {
                        "url": f"http://127.0.0.1:{sock.getsockname()[1]}",
                        "workflow": WORKFLOW,
                        "projectId": project,
                        "grant": grant,
                        "session": [p.session_cookie.name, a],
                        "csrf": [p.csrf_cookie.name, b],
                    }
                ),
                capture_output=True,
                text=True,
                timeout=65,
            )
            assert result.returncode == 0, result.stdout + result.stderr
            print(result.stdout)
        finally:
            srv.should_exit = True
            thread.join(timeout=10)
    assert len(seen) == 1
    path, body = seen[0]
    assert (
        path == "/v1/messages"
        and body["model"] == "fixture-number"
        and body["max_tokens"] == 128
    )
    assert "tools" not in body or not body["tools"]
    assert "config.threads" in json.dumps(body) and sid not in json.dumps(body)
    with closing(model.db()) as db:
        assert tuple(
            db.execute("SELECT requests,minor_left FROM grants").fetchone()
        ) == (0, 0)


def test_worker_does_not_redirect_retry_or_reuse_old_binding(tmp_path):
    calls = []

    class Reject(BaseHTTPRequestHandler):
        def do_POST(self):
            calls.append(self.path)
            self.send_response(307 if self.path.startswith("/redirect/") else 429)
            self.send_header("Location", "/unexpected")
            self.end_headers()

        def log_message(self, *args):
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Reject)
    thread = Thread(target=server.serve_forever)
    thread.start()
    worker = Path(os.environ["HW_TEST_PERI_MODEL_WORKER"]).resolve(strict=True)
    settings = tmp_path / "profile"

    def configure(prefix):
        settings.write_text(
            json.dumps(
                {
                    "config": {
                        "active_alias": "fable",
                        "profiles": {"fable": {"provider": "fixture"}},
                        "providers": [
                            {
                                "id": "fixture",
                                "type": "anthropic",
                                "apiKey": "fixture",
                                "models": {"fable": "number"},
                                "baseUrl": f"http://127.0.0.1:{server.server_port}/{prefix}",
                            }
                        ],
                    }
                }
            )
        )

    def invoke(body):
        return subprocess.run(
            [str(worker), str(settings), "--loopback-test"],
            input=json.dumps(body),
            capture_output=True,
            text=True,
            timeout=6,
        )

    try:
        for i, prefix in enumerate(["redirect", "rate"]):
            configure(prefix)
            binding = json.loads(invoke({"op": "binding"}).stdout)
            arg = {
                "op": "suggest",
                "binding": binding,
                "field": "config.threads",
                "value": 3,
                "output_tokens": 128,
            }
            assert invoke(arg).returncode != 0
            assert len(calls) == i + 1 and calls[-1] == f"/{prefix}/v1/messages"
        configure("changed")
        assert invoke(arg).returncode != 0 and len(calls) == 2
    finally:
        server.shutdown()
        thread.join()
        server.server_close()
