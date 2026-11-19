"""One temporary SQLite for catalog selection, WashU, and the saved download.

The selected experiment, the checked track, and the relocated download are the
same run, artifact, generation, and revision. WashU may transform track bytes.
The relocation digest is compared with the original download body.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta, timezone
from hashlib import sha256
import json
from pathlib import Path
import sys
import threading
import time
from uuid import uuid4

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "test"))
sys.path.insert(0, str(ROOT / "test" / "api"))
sys.path.insert(0, str(ROOT / "test" / "browser"))

import uvicorn  # noqa: E402

from catalog_page_server import (  # noqa: E402
    ADMIN_PASSWORD,
    ADMIN_USERNAME,
    _known,
    _seed,
)
from encode_pipeline.api.main import create_app  # noqa: E402
from encode_pipeline.deployment.frontend import ManifestStaticApplication  # noqa: E402
from encode_pipeline.frontend_assets import load_packaged_frontend_assets  # noqa: E402
from encode_pipeline.platform.adapters import WorkflowInputs  # noqa: E402
from encode_pipeline.platform.builds import WorkflowBuildIdentity  # noqa: E402
from encode_pipeline.platform.data_registry import ProjectSampleSelection  # noqa: E402
from encode_pipeline.platform.runs import RunRecord, RunStatus  # noqa: E402
from encode_pipeline.platform.snapshots import (  # noqa: E402
    PAYLOAD_DIGEST_SCHEME,
    VALIDATION_EVIDENCE_OUTCOME,
    ValidatedInputSnapshot,
    build_workflow_inputs_digest,
    canonical_workflow_inputs_json,
)
from encode_pipeline.services.artifact_disk_relocation import (  # noqa: E402
    ArtifactDiskRelocation,
)
from encode_pipeline.services.authentication_service import (  # noqa: E402
    AccountAdministrationService,
    AuthenticationActor,
)
from encode_pipeline.services.data_registry import DataRegistryService  # noqa: E402
from encode_pipeline.services.run_repositories import RunEventDraft  # noqa: E402
from test_routes_browser_tracks import (  # noqa: E402
    WORKFLOW_ID,
    _bind_reference,
    _record,
)


PEAK_NAME = "shared-peak.narrowPeak"
ARTIFACT_ID = "sharedpeak"
SOURCE_CONTENT = (
    b"chr1\t100\t200\tshared-peak\t50\t+\nchr1\t300\t450\tshared-b\t80\t+\n"
)
EXPERIMENT_NAME = "Treg shared peak"


class _Client:
    def __init__(self, app) -> None:
        self.app = app


def _bind_run(
    app, project_id: str, sample_revision_id: str, experiment_revision_id: str
) -> str:
    now = datetime.now(timezone.utc)
    canonical = canonical_workflow_inputs_json(WorkflowInputs(config={}))
    build_identity = WorkflowBuildIdentity(
        workflow_id=WORKFLOW_ID,
        adapter_version="combined-file-v1",
        scheme="sha256-tree-v1",
        logical_entrypoint="workflow/Snakefile",
        digest="a" * 64,
        captured_at=now,
    )
    snapshot = ValidatedInputSnapshot(
        snapshot_id=f"vsnap_{uuid4().hex}",
        workflow_id=WORKFLOW_ID,
        adapter_version="combined-file-v1",
        schema_version="1.0.0",
        schema_dialect="https://json-schema.org/draft/2020-12/schema",
        workflow_build_identity=build_identity,
        canonical_payload=canonical,
        payload_digest_scheme=PAYLOAD_DIGEST_SCHEME,
        payload_digest=build_workflow_inputs_digest(canonical),
        validation_outcome=VALIDATION_EVIDENCE_OUTCOME,
        validation_issue_codes=(),
        validated_at=now,
        expires_at=now + timedelta(hours=1),
    )
    repository = app.state.persistence.repository
    repository.create_validated_input_snapshot(
        snapshot,
        project_sample_selection=ProjectSampleSelection(
            project_id=project_id,
            sample_revision_ids=(sample_revision_id,),
        ),
        experiment_revision_id=experiment_revision_id,
    )
    run_id = str(uuid4())
    record = RunRecord(
        run_id=run_id,
        workflow_id=WORKFLOW_ID,
        inputs=snapshot.to_workflow_inputs().to_dict(),
        status=RunStatus.CREATED,
        created_at=now,
        updated_at=now,
        started_at=None,
        ended_at=None,
        current_stage=None,
        cancellation_reason=None,
        error=None,
        tags={},
    )
    repository.consume_validated_input_snapshot(
        snapshot.snapshot_id,
        workflow_id=WORKFLOW_ID,
        expected_build_identity=build_identity,
        record=record,
        consumed_at=now,
        event=RunEventDraft(
            event_type="status_changed",
            message="Run created.",
            status=RunStatus.CREATED,
            context={
                "previous_status": None,
                "new_status": RunStatus.CREATED.value,
            },
        ),
    )
    run_service = app.state.run_service
    for status in (
        RunStatus.VALIDATING,
        RunStatus.PLANNED,
        RunStatus.QUEUED,
        RunStatus.RUNNING,
        RunStatus.SUCCEEDED,
    ):
        run_service.transition_run(run_id, status)
    return run_id


def _relocate(app, durable: Path, identity: dict[str, str]) -> dict[str, str]:
    relocation = ArtifactDiskRelocation(
        source_root=app.state.workspace_root,
        durable_root=durable,
    )
    sealed = app.state.artifact_download_service.capture_relocation(
        relocation,
        run_id=identity["run_id"],
        artifact_id=identity["artifact_id"],
        expected_generation=identity["generation"],
        expected_revision=identity["revision"],
    )
    journal = sealed
    for _step in range(4):
        if journal.stage == "switched":
            break
        journal = relocation.resume(
            run_id=identity["run_id"],
            artifact_id=identity["artifact_id"],
        )
    if journal.stage != "switched":
        raise RuntimeError(f"relocation stopped at {journal.stage}")
    return {
        "stage": journal.stage,
        "content_sha256": journal.content_sha256,
        "size_bytes": str(journal.size_bytes),
        "generation": journal.generation,
        "revision": journal.revision,
        "run_id": journal.run_id,
        "artifact_id": journal.artifact_id,
        "relative_path": journal.relative_path,
    }


def _watch(
    app, request_path: Path, result_path: Path, identity: dict[str, str], durable: Path
) -> None:
    while not request_path.is_file():
        time.sleep(0.1)
    try:
        result = _relocate(app, durable, identity)
        result["ok"] = "true"
    except Exception as exc:  # noqa: BLE001 - the browser waits on this file
        result = {"ok": "false", "error": f"{type(exc).__name__}: {exc}"}
    result_path.write_text(json.dumps(result), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, required=True)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--durable", type=Path, required=True)
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    args = parser.parse_args()
    args.workspace.mkdir(parents=True, exist_ok=True)
    args.durable.mkdir(parents=True, exist_ok=True)
    app = create_app(
        database_url=f"sqlite:///{args.database}",
        workspace_root=args.workspace,
        relocation_root=args.durable,
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
    accounts.create_member(
        AuthenticationActor.local_operator(),
        "catalog-member",
        "catalog-browser-member",
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
    server = uvicorn.Server(
        uvicorn.Config(
            app,
            host=args.host,
            port=args.port,
            log_level="warning",
            access_log=False,
        )
    )
    server.install_signal_handlers = lambda: None
    thread = threading.Thread(target=server.run, daemon=True)
    thread.start()
    base = f"http://{args.host}:{args.port}"
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
                    "values": [_known_pin(seeded["cell_type"])],
                },
                "cell_line": {"status": "not_applicable", "values": []},
                "treatment": {"status": "not_applicable", "values": []},
                "tag": {"status": "not_applicable", "values": []},
            },
        },
        "catalog-treg",
    )
    sample_revision_id = classified["sample_revision_id"]
    experiment = seeded["post"](
        f"/api/v1/catalog/projects/{project.project_id}/experiments",
        {
            "stable_key": "shared-peak",
            "display_name": EXPERIMENT_NAME,
            "assay": _known(seeded["cuttag"]),
            "target": _known(seeded["target"]),
            "batch": _known(seeded["batch"]),
            "samples": [
                {
                    "sample_revision_id": sample_revision_id,
                    "replicate_role": _known(seeded["role"]),
                }
            ],
        },
        "experiment-shared-peak",
    )
    run_id = _bind_run(
        app,
        project.project_id,
        sample_revision_id,
        experiment["experiment_revision_id"],
    )
    client = _Client(app)
    _bind_reference(client, run_id, assembly="hg38")
    artifact = _record(
        client,
        run_id,
        ARTIFACT_ID,
        SOURCE_CONTENT,
        name=PEAK_NAME,
        output_type="macs3_peak",
        revision_kind="content",
        reference_assembly="hg38",
    )
    generation = app.state.run_service.get_result_state(run_id).artifact_generation
    if not isinstance(generation, str):
        raise RuntimeError("shared peak has no artifact generation")
    identity = {
        "run_id": run_id,
        "artifact_id": artifact.artifact_id,
        "generation": generation,
        "revision": artifact.revision,
    }
    threading.Thread(
        target=_watch,
        args=(app, args.request, args.result, identity, args.durable),
        daemon=True,
    ).start()
    manifest = {
        "api_base": base,
        "admin_username": ADMIN_USERNAME,
        "admin_password": ADMIN_PASSWORD,
        "project_id": project.project_id,
        "assay_term_id": seeded["cuttag"]["term_id"],
        "cell_type_term_id": seeded["cell_type"]["term_id"],
        "experiment_id": experiment["experiment_id"],
        "experiment_revision_id": experiment["experiment_revision_id"],
        "experiment_name": EXPERIMENT_NAME,
        "sample_revision_id": sample_revision_id,
        "run_id": run_id,
        "artifact_id": artifact.artifact_id,
        "generation": generation,
        "revision": artifact.revision,
        "source_sha256": sha256(SOURCE_CONTENT).hexdigest(),
        "source_size": len(SOURCE_CONTENT),
        "download_path": (
            f"/api/v1/runs/{run_id}/artifacts/{artifact.artifact_id}/download"
            f"?generation={generation}&revision={artifact.revision}"
        ),
        "relocate_request": str(args.request),
        "relocate_result": str(args.result),
        "peak_name": PEAK_NAME,
    }
    print(f"COMBINED_FILE_READY {json.dumps(manifest)}", flush=True)
    try:
        while thread.is_alive():
            thread.join(timeout=3600)
    finally:
        server.should_exit = True


def _known_pin(term: dict) -> dict:
    return {
        "term_id": term["term_id"],
        "term_revision_id": term["revision"]["term_revision_id"],
    }


if __name__ == "__main__":
    main()
