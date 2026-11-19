"""Serve a paged run whose generation can change between candidate pages."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import uvicorn

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "test"))
sys.path.insert(0, str(ROOT / "test" / "api"))

from encode_pipeline.api.main import create_app  # noqa: E402
from encode_pipeline.deployment.frontend import ManifestStaticApplication  # noqa: E402
from encode_pipeline.frontend_assets import load_packaged_frontend_assets  # noqa: E402
from encode_pipeline.services.browser_tracks import (  # noqa: E402
    BrowserTrackService,
    _ARTIFACT_PAGE_SIZE,
)
from conftest import seed_test_authentication  # noqa: E402
from test_routes_browser_tracks import (  # noqa: E402
    _bind_reference,
    _create_run,
    _generation,
    _record,
    _report_artifact,
)


class _Client:
    def __init__(self, app) -> None:
        self.app = app


def _arm_generation_change(service, run_id: str, flag: Path) -> None:
    real_page = service.list_artifacts_page
    bumped = False

    def change_between_pages(
        requested_run_id,
        *,
        expected_generation,
        after=None,
        limit=None,
    ):
        nonlocal bumped
        if (
            not bumped
            and after is not None
            and requested_run_id == run_id
            and flag.is_file()
        ):
            bumped = True
            current = service.list_artifacts(requested_run_id)
            service.replace_artifacts(
                requested_run_id,
                (*current, _report_artifact(requested_run_id, "report-extra")),
            )
        return real_page(
            requested_run_id,
            expected_generation=expected_generation,
            after=after,
            limit=limit,
        )

    service.list_artifacts_page = change_between_pages


def build_app(database: Path, workspace: Path):
    app = create_app(
        database_url=f"sqlite:///{database}",
        workspace_root=workspace,
    )
    seed_test_authentication(app)
    client = _Client(app)
    run_id = _create_run(client)
    _bind_reference(client, run_id, assembly="hg38")
    reports = tuple(
        _report_artifact(run_id, f"report-{index:04d}")
        for index in range(_ARTIFACT_PAGE_SIZE + 1)
    )
    app.state.run_service.replace_artifacts(run_id, reports)
    track_a = _record(
        client,
        run_id,
        "track-a",
        b"chr1\t1500\t1700\tpeak-a\t80\t+\n",
        name="track-a.narrowPeak",
        output_type="macs3_peak",
        revision_kind="content",
        reference_assembly="hg38",
    )
    track_b = _record(
        client,
        run_id,
        "track-b",
        b"chr1\t1200\t1400\tchr1\t1600\t1800\tread\t211\t+\t+\n",
        name="track-b.bedpe.gz",
        output_type="hitrac_bedpe_all",
        revision_kind="content",
        reference_assembly="hg38",
    )
    track_c = _record(
        client,
        run_id,
        "track-c",
        b"chr1\t1800\t1900\tpeak-c\t40\t+\n",
        name="track-c.narrowPeak",
        output_type="macs3_peak",
        revision_kind="content",
        reference_assembly="hg38",
    )
    generation = _generation(client, run_id)
    browser = BrowserTrackService(
        run_service=app.state.run_service,
        download_service=app.state.artifact_download_service,
        workspace_root=workspace,
        reference_profile_service=app.state.reference_profile_service,
    )
    for artifact in (track_a, track_b, track_c):
        browser.prepare(
            run_id,
            artifact.artifact_id,
            generation=generation,
            revision=artifact.revision,
        )
    flag = workspace / "arm-generation-change"
    _arm_generation_change(app.state.run_service, run_id, flag)
    app.mount(
        "/",
        ManifestStaticApplication(load_packaged_frontend_assets()),
        name="frontend",
    )
    return app, run_id, flag


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    args = parser.parse_args()
    _app, run_id, flag = build_app(args.database, args.workspace)
    print(f"WASHU_SELECTION run_id={run_id} arm={flag}", flush=True)
    uvicorn.run(_app, host="127.0.0.1", port=args.port, log_level="warning")


if __name__ == "__main__":
    main()
