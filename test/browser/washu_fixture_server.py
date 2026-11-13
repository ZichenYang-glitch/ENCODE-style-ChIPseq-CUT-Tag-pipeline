"""Serve one seeded run so a real WashU page can read same-origin tracks."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys

import uvicorn

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "test"))
sys.path.insert(0, str(ROOT / "test" / "api"))

from encode_pipeline.api.main import create_app  # noqa: E402
from encode_pipeline.services.browser_tracks import BrowserTrackService  # noqa: E402
from conftest import seed_test_authentication  # noqa: E402
from test_routes_browser_tracks import (  # noqa: E402
    FIXTURE,
    _bind_reference,
    _create_run,
    _generation,
    _record,
)


class _Client:
    def __init__(self, app) -> None:
        self.app = app


def build_app(database: Path, workspace: Path):
    app = create_app(
        database_url=f"sqlite:///{database}",
        workspace_root=workspace,
    )
    seed_test_authentication(app)
    client = _Client(app)
    run_id = _create_run(client)
    _bind_reference(client, run_id, assembly="hg38")
    signal = _record(
        client,
        run_id,
        "artifact-signal",
        FIXTURE.read_bytes(),
        name="chr1-signal.bw",
        output_type="cpm_bigwig",
        reference_assembly="hg38",
    )
    peak = _record(
        client,
        run_id,
        "artifact-peak",
        b"chr1\t1500\t1700\tpeak-a\t80\t+\n",
        name="sample_peaks.narrowPeak",
        output_type="macs3_peak",
        revision_kind="content",
        reference_assembly="hg38",
    )
    pet = _record(
        client,
        run_id,
        "artifact-pet",
        b"chr1\t1200\t1400\tchr1\t1600\t1800\tread\t211\t+\t+\n",
        name="sample_all.bedpe.gz",
        output_type="hitrac_bedpe_all",
        revision_kind="content",
        reference_assembly="hg38",
    )
    generation = _generation(client, run_id)
    service = BrowserTrackService(
        run_service=app.state.run_service,
        download_service=app.state.artifact_download_service,
        workspace_root=workspace,
        reference_profile_service=app.state.reference_profile_service,
    )
    for artifact in (pet, peak):
        service.prepare(
            run_id,
            artifact.artifact_id,
            generation=generation,
            revision=artifact.revision,
        )
    return app, run_id, signal.artifact_id, pet.artifact_id


def _fail_after_preflight(app):
    """Allow the one-byte preflight, then fail real data and index reads."""

    async def wrapped(scope, receive, send):
        if scope.get("type") != "http":
            await app(scope, receive, send)
            return
        path = scope.get("path") or ""
        if "/browser-tracks/" not in path or "/bytes/" not in path:
            await app(scope, receive, send)
            return
        headers = {
            name.decode("latin1").lower(): value.decode("latin1")
            for name, value in scope.get("headers", ())
        }
        if headers.get("range") == "bytes=0-0":
            await app(scope, receive, send)
            return
        body = (
            b'{"ok":false,"issues":[{"code":"BROWSER_TRACK_UNAVAILABLE",'
            b'"message":"Track bytes are temporarily unavailable.",'
            b'"severity":"error","path":"file_token","source":"runtime",'
            b'"context":{}}]}'
        )
        await send(
            {
                "type": "http.response.start",
                "status": 500,
                "headers": [
                    (b"content-type", b"application/json"),
                    (b"content-length", str(len(body)).encode("ascii")),
                ],
            }
        )
        await send({"type": "http.response.body", "body": body})

    return wrapped


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    args = parser.parse_args()
    app, run_id, signal_id, pet_id = build_app(args.database, args.workspace)
    if os.environ.get("WASHU_PROOF_FAIL_AFTER_PREFLIGHT") == "1":
        app = _fail_after_preflight(app)
    print(
        f"WASHU_FIXTURE run_id={run_id} signal={signal_id} pet={pet_id}",
        flush=True,
    )
    uvicorn.run(app, host="127.0.0.1", port=args.port, log_level="warning")


if __name__ == "__main__":
    main()
