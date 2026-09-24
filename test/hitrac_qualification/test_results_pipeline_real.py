"""Real tiny science -> original result bundle -> SQLite/API/download.

Explicit real_execution tier; no queue/worker/browser claim. Local temporary
SQLite and synthetic auth are used; artifacts and QC are never hand-filled.
"""

import asyncio
from collections import Counter
from dataclasses import asdict, replace
from datetime import datetime, timedelta, timezone
import hashlib
import json
import os
from pathlib import Path
import runpy

import httpx
import pytest
from sqlalchemy import text

from encode_pipeline.api.main import create_app
from encode_pipeline.adapters.hitrac_preprocess.admission import sha256_file
from encode_pipeline.adapters.hitrac_preprocess.execution import admit_runtime
from encode_pipeline.adapters.hitrac_preprocess.results import (
    HiTracPreprocessResultsAdapter,
    METRIC_KEYS,
)
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.authentication import UserRole, UserStatus
from encode_pipeline.services.authentication import SessionSecrets, new_session_record
from encode_pipeline.platform.planning import ExecutionPlan, PlanStatus
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.artifact_extraction import ArtifactExtractionService
from encode_pipeline.services.command_builder import CommandBuilder
from encode_pipeline.services.defaults import create_default_process_runner
from encode_pipeline.services.materialization import WorkspaceMaterializer
from encode_pipeline.services.planning import WorkspacePlanner
from test_qualification_real import (
    coordinates as shared_coordinates,
    expected_rows,
    rows,
)

coordinates = shared_coordinates
pytestmark = pytest.mark.real_execution


@pytest.mark.parametrize(
    "sample_ids",
    [("alpha", "beta"), ("user sample 1", "user sample 2"), ("a  b", "a b ")],
    ids=["existing-qc-id-grammar", "h1-space-ids", "consecutive-trailing-spaces"],
)
def test_real_results_snapshot_bundle_http_download(
    coordinates, tmp_path, monkeypatch, sample_ids
):
    source = Path(__file__).resolve().parents[2]
    for name in tuple(os.environ):
        if name.startswith(("ENCODE_PIPELINE_", "HELIXWEAVE_")):
            monkeypatch.delenv(name)
    monkeypatch.setenv("HELIXWEAVE_TERMINAL_EMAIL_ENABLED", "false")
    monkeypatch.setenv("ENCODE_PIPELINE_REDIS_URL", "redis://127.0.0.1:1/0")
    reference_config = tmp_path / "reference-profiles.json"
    reference_config.write_text(
        json.dumps(
            {
                "schema_version": "helixweave-reference-profiles-v1",
                "profiles": {
                    "tiny": {
                        "bindings": {
                            "hitrac-preprocess": {
                                "schema_version": "hitrac-reference-profile-v1",
                                "binding": coordinates["REFERENCE_BINDING"],
                                "sha256": coordinates["REFERENCE_SHA256"],
                            }
                        }
                    }
                },
            }
        )
    )
    reference_config.chmod(0o600)
    monkeypatch.setenv(
        "ENCODE_PIPELINE_REFERENCE_PROFILE_CONFIG", str(reference_config)
    )
    rb = Path(coordinates["RUNTIME_BINDING"])
    runtime = admit_runtime(rb, sha256_file(rb), timeout=120)
    adapter = HiTracPreprocessResultsAdapter(runtime=runtime)
    registry = WorkflowRegistry([adapter])
    assert adapter.execution_availability().execution == "available"
    app = create_app(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        workspace_root=tmp_path / "workspaces",
        project_root=source,
        registry=registry,
    )
    helpers = runpy.run_path(str(source / "test/conftest.py"))
    account = helpers["seed_test_authentication"](app)
    member = replace(
        account, user_id="usr_" + "1" * 32, username="test-member", role=UserRole.MEMBER
    )
    app.state.authentication_repository.create_account(member)
    member_secrets = SessionSecrets._from_generated("c" * 43, "d" * 43)
    app.state.authentication_repository.create_session(
        new_session_record(
            user_id=member.user_id,
            secrets=member_secrets,
            created_at=datetime.now(timezone.utc),
            lifetime=timedelta(hours=8),
        )
    )
    reference = app.state.reference_profile_service.register(
        safe_key="tiny",
        display_name="Tiny synthetic reference",
        organism="synthetic",
        assembly="tiny",
        config_key="tiny",
    )
    app.state.reference_profile_service.enable(
        reference.profile_id, revision_id=reference.revision_id
    )
    directory = Path(coordinates["TINY_INPUTS"]) / "positive/fastq"
    samples = [
        {
            "sample_id": sample_ids[index - 1],
            "fastq_1": str(path),
            "fastq_2": str(path.with_name(path.name.replace("_R1", "_R2"))),
        }
        for index, path in enumerate(sorted(directory.glob("*_R1.fastq.gz")), 1)
    ]
    inputs = WorkflowInputs(
        config={}, samples=samples, options={"threads": 2, "mapq": 17}
    )
    before = {
        p: sha256_file(Path(p))
        for row in samples
        for k, p in row.items()
        if k.startswith("fastq")
    }
    snapshot = app.state.validated_input_service.validate(
        adapter.metadata.workflow_id,
        inputs,
        reference_profile_revision_id=reference.revision_id,
    )
    assert snapshot.is_success, snapshot.issues
    creation = app.state.validated_run_creation_service.create_run(
        adapter.metadata.workflow_id,
        snapshot.value.snapshot_id,
        requested_by_user_id=account.user_id,
    )
    assert creation.created
    record = creation.record
    service = app.state.run_service
    resolver = app.state.reference_profile_resolver
    workspace = app.state.workspace_root / record.run_id
    plan = ExecutionPlan(
        plan_id="h4-real",
        run_id=record.run_id,
        workflow_id=record.workflow_id,
        status=PlanStatus.PENDING,
        inputs_snapshot=record.inputs,
    )
    planned = WorkspacePlanner(
        registry, reference_profile_resolver=resolver
    ).plan_workspace(plan, workspace)
    assert planned.is_success, planned.issues
    built = CommandBuilder(registry, reference_profile_resolver=resolver).build_command(
        planned.value, workspace
    )
    assert built.is_success, built.issues
    materialized = WorkspaceMaterializer().materialize(
        built.value.workspace_plan, workspace
    )
    assert materialized.is_success, materialized.issues
    service.transition_run(record.run_id, RunStatus.VALIDATING)
    service.complete_preflight(
        record.run_id,
        snapshot.value.workflow_build_identity,
        message="Isolated H4 qualification planning completed.",
    )
    service.transition_run(record.run_id, RunStatus.QUEUED)
    service.transition_run(record.run_id, RunStatus.RUNNING)
    runner = create_default_process_runner(
        registry=registry, settings=app.state.worker_settings
    )
    result = runner.run(built.value.command_spec)
    (tmp_path / "runner-result.json").write_text(
        json.dumps(
            {
                "argv": built.value.command_spec.argv,
                "cwd": built.value.command_spec.cwd,
                "result": None if result.value is None else vars(result.value),
                "issues": [i.to_dict() for i in result.issues],
            },
            indent=2,
        )
    )
    assert result.is_success and result.value.exit_code == 0, result
    service.transition_run(record.run_id, RunStatus.SUCCEEDED)
    complete = json.loads((workspace / "hitrac-attempt/complete.json").read_text())
    expected = expected_rows(coordinates["design"])
    wanted = coordinates["design"]["scenarios"]["positive"][
        "conditional_expected_summary"
    ]
    for token, sample in complete["results"]["samples"].items():
        assert sample["all"] == 8 and sample["noBg"] == 5
        assert sample["metrics"] == pytest.approx(wanted, rel=0, abs=1e-12)
        output = workspace / "hitrac-attempt/output" / token
        assert Counter(map(tuple, rows(output / f"{token}_all.bedpe.gz"))) == Counter(
            expected
        )
        assert Counter(
            map(tuple, rows(output / f"{token}_unique.bedpe.gz"))
        ) == Counter(expected[i] for i in (0, 3, 4, 5, 6))
    extractor = ArtifactExtractionService(
        run_service=service,
        registry=registry,
        build_identity_provider=app.state.build_identity_provider,
        workspace_root=app.state.workspace_root,
        reference_profile_resolver=resolver,
    )
    attempt = extractor.begin_attempt(record.run_id)
    extracted = extractor.extract(record.run_id, attempt_id=attempt)
    assert extracted.is_success, extracted.issues
    assert len(extracted.value) == 5
    state = service.get_result_state(record.run_id)
    assert state.artifact_outcome == state.qc_outcome == "succeeded"
    assert state.artifact_generation == state.qc_artifact_generation
    events = service.list_events(record.run_id, limit=1000)
    repeated = extractor.extract(record.run_id, attempt_id=attempt)
    assert repeated.is_success, repeated.issues
    assert service.get_result_state(record.run_id) == state
    assert service.list_events(record.run_id, limit=1000) == events

    async def http_checks():
        observations = {}
        cookie = app.state.auth_cookie_policy.session_cookie.name
        async with app.router.lifespan_context(app):
            transport = httpx.ASGITransport(app=app, raise_app_exceptions=False)
            async with httpx.AsyncClient(
                transport=transport, base_url="http://test"
            ) as anon:
                for suffix in ("artifacts", "qc-metrics"):
                    response = await anon.get(f"/api/v1/runs/{record.run_id}/{suffix}")
                    assert response.status_code == 401
            async with httpx.AsyncClient(
                transport=transport,
                base_url="http://test",
                cookies={cookie: app.state.test_auth_tokens[0]},
            ) as client:
                listed = await client.get(
                    f"/api/v1/runs/{record.run_id}/artifacts", params={"limit": 100}
                )
                qc = await client.get(
                    f"/api/v1/runs/{record.run_id}/qc-metrics", params={"limit": 100}
                )
                assert listed.status_code == qc.status_code == 200
                artifacts, metrics = listed.json(), qc.json()
                assert (
                    len(artifacts["artifacts"]) == 5
                    and len(metrics["qc_metrics"]) == 30
                )
                assert artifacts["artifact_generation"] == state.artifact_generation
                assert metrics["qc_generation"] == state.qc_generation
                for row in samples:
                    table = {
                        m["metric_key"]: m
                        for m in metrics["qc_metrics"]
                        if m["sample_id"] == row["sample_id"]
                    }
                    assert set(table) == set(METRIC_KEYS)
                    assert [
                        float(table[key]["value"]) for key in METRIC_KEYS
                    ] == pytest.approx(wanted, rel=0, abs=1e-12)
                    assert all(
                        table[key]["unit"]
                        == (
                            "count"
                            if index in (0, 1, 3, 9)
                            else "ratio"
                            if index == 10
                            else "fraction"
                        )
                        for index, key in enumerate(METRIC_KEYS)
                    )
                hashes = {}
                for artifact in artifacts["artifacts"]:
                    assert (
                        artifact["name"].endswith(".bedpe.gz")
                        or artifact["name"] == "tracPre_summary.txt"
                    )
                    downloaded = await client.get(
                        f"/api/v1/runs/{record.run_id}/artifacts/{artifact['artifact_id']}/download",
                        params={
                            "generation": artifacts["artifact_generation"],
                            "revision": artifact["revision"],
                        },
                    )
                    assert downloaded.status_code == 200
                    assert "attachment;" in downloaded.headers["content-disposition"]
                    original = next(
                        (workspace / "hitrac-attempt/output").rglob(artifact["name"])
                    )
                    assert downloaded.content == original.read_bytes()
                    hashes[artifact["artifact_id"]] = hashlib.sha256(
                        downloaded.content
                    ).hexdigest()
                for private in ("bam", "fastq", "private-request", "call-plan"):
                    response = await client.get(
                        f"/api/v1/runs/{record.run_id}/artifacts/{private}/download",
                        params={
                            "generation": state.artifact_generation,
                            "revision": artifacts["artifacts"][0]["revision"],
                        },
                    )
                    assert response.status_code == 404
                public = json.dumps({"artifacts": artifacts, "qc": metrics})
                assert str(tmp_path) not in public and str(rb) not in public
                observations.update(
                    artifacts=artifacts, qc=metrics, artifact_sha256=hashes
                )
            # Existing LAN contract allows any enabled member to read results;
            # it does not define owner-only run access. Administrator actions
            # stay forbidden and a disabled member cannot list or download.
            async with httpx.AsyncClient(
                transport=transport, base_url="http://test", cookies={cookie: "c" * 43}
            ) as member_client:
                shared = await member_client.get(
                    f"/api/v1/runs/{record.run_id}/artifacts"
                )
                assert shared.status_code == 200
                admin_only = await member_client.get("/api/v1/auth/accounts")
                assert admin_only.status_code == 403
                app.state.authentication_repository.save_account(
                    replace(member, status=UserStatus.DISABLED)
                )
                for suffix in ("artifacts", "qc-metrics"):
                    denied = await member_client.get(
                        f"/api/v1/runs/{record.run_id}/{suffix}"
                    )
                    assert denied.status_code == 401
                item = artifacts["artifacts"][0]
                denied = await member_client.get(
                    f"/api/v1/runs/{record.run_id}/artifacts/{item['artifact_id']}/download",
                    params={
                        "generation": state.artifact_generation,
                        "revision": item["revision"],
                    },
                )
                assert denied.status_code == 401
            with app.state.persistence.engine.connect() as connection:
                counts = {
                    table: connection.execute(
                        text(f"SELECT count(*) FROM {table} WHERE run_id=:run_id"),
                        {"run_id": record.run_id},
                    ).scalar_one()
                    for table in (
                        "run_artifacts",
                        "run_qc_metrics",
                        "artifact_publications",
                    )
                }
                assert counts == {
                    "run_artifacts": 5,
                    "run_qc_metrics": 30,
                    "artifact_publications": 5,
                }
            observations["database_counts"] = counts
        return observations

    observed = asyncio.run(http_checks())
    assert {p: sha256_file(Path(p)) for p in before} == before
    (tmp_path / "http-observations.json").write_text(json.dumps(observed, indent=2))
    (tmp_path / "result-state.json").write_text(
        json.dumps(asdict(state), indent=2, default=str)
    )
    private = tmp_path / "browser-source.json"
    private.write_text(
        json.dumps(
            {
                "database_url": app.state.database_url,
                "workspace_root": str(app.state.workspace_root),
                "reference_config": str(reference_config),
                "runtime_binding": str(rb),
                "runtime_sha256": sha256_file(rb),
                "run_id": record.run_id,
                "samples": [row["sample_id"] for row in samples],
                "cookie_name": app.state.auth_cookie_policy.session_cookie.name,
                "session_token": app.state.test_auth_tokens[0],
                "artifact_sha256": observed["artifact_sha256"],
                "implementation_sha256": complete["identity"]["sha256"],
            },
            indent=2,
        )
    )
    private.chmod(0o600)
