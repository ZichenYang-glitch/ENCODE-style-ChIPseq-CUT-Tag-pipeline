"""IRENE-39: register one external file and download an exact revision."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from dataclasses import replace
from hashlib import sha256
import asyncio
import json
import os
from pathlib import Path
import socket
import stat
import subprocess
import threading
import time
from types import SimpleNamespace
from unittest.mock import patch

import httpx
import uvicorn

from sqlalchemy import text

from encode_pipeline.api.main import create_app
from encode_pipeline.platform.authentication import (
    SessionRevocationReason,
    UserRole,
)
from encode_pipeline.platform.security_audit import (
    AuditResourceKind,
    build_audit_resource,
)
from encode_pipeline.services.authentication import (
    digest_session_token,
    new_session_record,
    new_session_secrets,
)
from encode_pipeline.services.data_registry import DataRegistryService
from encode_pipeline.services.input_registry import InputRegistryService
from api_test_client import ApiTestClient, seeded_auth_async_client
from conftest import (
    TEST_AUTH_CSRF_TOKEN,
    TEST_AUTH_SESSION_TOKEN,
    seed_test_authentication,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
FIXTURE = REPOSITORY_ROOT / "test/fixtures/irene39/fixture39.tsv"
PUBLIC_KEYS = {
    "project_id",
    "result_id",
    "result_revision_id",
    "revision_number",
    "revision_status",
    "display_name",
    "format",
    "schema_version",
    "origin",
    "method",
    "reference",
    "bytes",
    "created_at",
}
FORBIDDEN = {
    "run_id",
    "stable_key",
    "storage_pool_id",
    "relative_path",
    "absolute_path",
    "config_key",
}


def test_external_cluster_counts_register_rename_and_exact_download(
    tmp_path: Path,
) -> None:
    fixture_bytes = FIXTURE.read_bytes()
    fixture_sha = sha256(fixture_bytes).hexdigest()
    pool_root = (tmp_path / "pool").resolve()
    source = pool_root / "external" / "fixture39.tsv"
    source.parent.mkdir(parents=True)
    source.write_bytes(fixture_bytes)
    config_path = tmp_path / "storage-pools.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": "storage-pool-roots-v1",
                "storage_pool_roots": {"fixture39-pool": str(pool_root)},
            }
        ),
        encoding="utf-8",
    )
    workspace = tmp_path / "workspace"
    app = create_app(
        database_url=f"sqlite:///{tmp_path / 'irene39.db'}",
        workspace_root=workspace,
        project_root=REPOSITORY_ROOT,
    )
    try:
        admin = seed_test_authentication(app)
        app.state.registered_result_service._storage_pool_config_path = config_path
        member = _member(app, admin)
        project = DataRegistryService(
            repository=app.state.persistence.data_registry_repository
        ).create_project("IRENE-39 fixture")
        other = DataRegistryService(
            repository=app.state.persistence.data_registry_repository
        ).create_project("IRENE-39 other")
        pool = InputRegistryService(
            repository=app.state.persistence.input_registry_repository
        )
        storage = pool.create_storage_pool(
            display_name="Fixture 39",
            config_key="fixture39-pool",
        )
        pool.bind_project_storage_pool(project.project_id, storage.storage_pool_id)
        pool.bind_project_storage_pool(other.project_id, storage.storage_pool_id)

        metadata = {
            "stable_key": "synthetic-cluster-counts",
            "display_name": "Synthetic single-cell cluster counts",
            "format": "tsv",
            "origin": {"kind": "external", "source": {"status": "unknown"}},
            "method": {"status": "unknown"},
            "reference": {"status": "unknown"},
            "pool_relative_path": "external/fixture39.tsv",
        }
        with ApiTestClient(app) as client:
            _use(app, member)
            denied = client.post(
                f"/api/v1/catalog/projects/{project.project_id}/registered-results",
                json=metadata,
                headers={"Idempotency-Key": "irene39-member-register"},
            )
            assert denied.status_code == 403, denied.text
            assert "synthetic-A" not in denied.text

            _use_admin(app)
            registered = client.post(
                f"/api/v1/catalog/projects/{project.project_id}/registered-results",
                json=metadata,
                headers={"Idempotency-Key": "irene39-register-1"},
            )
            assert registered.status_code == 200, registered.text
            revision = registered.json()["revision"]
            _assert_public(revision)
            assert revision["display_name"] == metadata["display_name"]
            assert revision["format"] == "tsv"
            assert revision["revision_number"] == 1
            assert revision["revision_status"] == "current"
            assert revision["origin"] == metadata["origin"]
            assert revision["method"] == {"status": "unknown"}
            assert revision["reference"] == {"status": "unknown"}
            assert revision["bytes"] == {
                "size_bytes": 50,
                "content_sha256": fixture_sha,
            }
            assert source.read_bytes() == fixture_bytes
            assert source.name == "fixture39.tsv"
            _assert_runs_unchanged(app)

            replay = client.post(
                f"/api/v1/catalog/projects/{project.project_id}/registered-results",
                json=metadata,
                headers={"Idempotency-Key": "irene39-register-1"},
            )
            assert replay.status_code == 200, replay.text
            assert (
                replay.json()["revision"]["result_revision_id"]
                == (revision["result_revision_id"])
            )
            conflict = client.post(
                f"/api/v1/catalog/projects/{project.project_id}/registered-results",
                json={**metadata, "display_name": "Different name"},
                headers={"Idempotency-Key": "irene39-register-1"},
            )
            assert conflict.status_code == 409, conflict.text
            assert conflict.json()["issues"][0]["code"] == (
                "catalog_idempotency_conflict"
            )

            downloaded = client.get(_content(project.project_id, revision))
            assert downloaded.status_code == 200, downloaded.text
            assert downloaded.content == fixture_bytes
            _assert_download_headers(downloaded)
            counts = _counts(downloaded.content)
            assert counts["synthetic-A"] == "3"
            assert counts["synthetic-B"] == "2"

            renamed = client.post(
                _revisions(project.project_id, revision["result_id"]),
                json={
                    "expected_revision_number": 1,
                    "display_name": ("Synthetic single-cell cluster counts (renamed)"),
                },
                headers={"Idempotency-Key": "irene39-rename-1"},
            )
            assert renamed.status_code == 200, renamed.text
            current = renamed.json()["revision"]
            assert current["result_id"] == revision["result_id"]
            assert current["result_revision_id"] != revision["result_revision_id"]
            assert current["revision_number"] == 2
            assert current["display_name"].endswith("(renamed)")
            assert current["bytes"] == revision["bytes"]
            assert current["origin"] == revision["origin"]
            assert source.read_bytes() == fixture_bytes
            historical = client.get(
                _revision(project.project_id, revision["result_revision_id"])
            )
            assert historical.status_code == 200, historical.text
            assert (
                historical.json()["revision"]["display_name"]
                == (metadata["display_name"])
            )
            assert historical.json()["revision"]["revision_status"] == "superseded"

            stale = client.post(
                _revisions(project.project_id, revision["result_id"]),
                json={
                    "expected_revision_number": 1,
                    "display_name": "stale",
                },
                headers={"Idempotency-Key": "irene39-stale-1"},
            )
            assert stale.status_code == 409, stale.text
            assert stale.json()["issues"][0]["code"] == "catalog_revision_conflict"
            confirmed = client.post(
                _revisions(project.project_id, revision["result_id"]),
                json={
                    "expected_revision_number": 2,
                    "display_name": "Synthetic single-cell cluster counts (confirmed)",
                },
                headers={"Idempotency-Key": "irene39-stale-1"},
            )
            assert confirmed.status_code == 200, confirmed.text
            assert confirmed.json()["revision"]["revision_number"] == 3
            assert confirmed.json()["revision"]["bytes"] == revision["bytes"]
            assert source.read_bytes() == fixture_bytes

            _use(app, member)
            old = client.get(_content(project.project_id, revision))
            assert old.status_code == 200, old.text
            assert old.content == fixture_bytes
            _assert_binary_contract(app)
            _generated_client_download(
                app,
                project_id=project.project_id,
                result_revision_id=revision["result_revision_id"],
                session_token=member.session_token,
                fixture_sha=fixture_sha,
            )
            member_events = _events_for(
                app,
                actor_user_id=member.user_id,
                project_id=project.project_id,
                result_id=revision["result_id"],
                result_revision_id=revision["result_revision_id"],
            )
            assert len(member_events) == 2
            successes = _success_downloads(app)
            blocked = _download_with_unreadable_copy(
                app,
                client,
                project.project_id,
                revision,
            )
            assert blocked.status_code == 409, blocked.text
            assert blocked.json()["issues"][0]["code"] == (
                "registered_result_bytes_unavailable"
            )
            assert "synthetic-A" not in blocked.text
            assert str(pool_root) not in blocked.text
            assert _success_downloads(app) == successes
            _assert_no_response_copies(workspace)
            failed_audit = _download_with_audit_failure(
                app,
                project.project_id,
                revision,
            )
            assert failed_audit.status_code == 500, failed_audit.text
            assert failed_audit.json()["issues"][0]["code"] == "INTERNAL_SERVER_ERROR"
            assert "audit unavailable" not in failed_audit.text
            assert "synthetic-A" not in failed_audit.text
            assert str(pool_root) not in failed_audit.text
            assert "registered-result-responses" not in failed_audit.text
            assert _success_downloads(app) == successes
            _assert_no_response_copies(workspace)
            assert member_events[-1].outcome.value == "succeeded"
            assert member_events[-1].reason_code is None
            assert revision["result_revision_id"] not in (
                member_events[-1].resource.resource_id
            )

            app.state.authentication_repository.revoke_session(
                digest_session_token(member.session_token),
                datetime.now(timezone.utc),
                SessionRevocationReason.LOGOUT,
            )
            revoked = client.get(_content(project.project_id, revision))
            assert revoked.status_code == 401, revoked.text
            assert "synthetic-A" not in revoked.text
            assert (
                _events_for(
                    app,
                    actor_user_id=member.user_id,
                    project_id=project.project_id,
                    result_id=revision["result_id"],
                    result_revision_id=revision["result_revision_id"],
                )
                == member_events
            )

            app.state.test_auth_tokens = None
            anonymous = client.get(_content(project.project_id, revision))
            assert anonymous.status_code == 401, anonymous.text
            assert "synthetic-A" not in anonymous.text

            _use_admin(app)
            wrong = client.get(_content(other.project_id, revision))
            assert wrong.status_code == 404, wrong.text
            assert wrong.json()["issues"][0]["code"] == "catalog_not_found"
            assert "synthetic-A" not in wrong.text
            assert str(pool_root) not in wrong.text
            assert (
                _events_for(
                    app,
                    actor_user_id=member.user_id,
                    project_id=project.project_id,
                    result_id=revision["result_id"],
                    result_revision_id=revision["result_revision_id"],
                )
                == member_events
            )

            seen: list[Path] = []

            def mutate(copy_path: Path) -> None:
                seen.append(copy_path)
                assert stat.S_IMODE(copy_path.stat().st_mode) == 0o600
                assert stat.S_IMODE(copy_path.parent.stat().st_mode) == 0o700
                response_root = workspace / "registered-result-responses"
                assert response_root in copy_path.parents
                source.write_bytes(b"synthetic-A\t9\n")

            app.state.registered_result_service.after_copy_verified = mutate
            try:
                sealed = client.get(_content(project.project_id, revision))
            finally:
                app.state.registered_result_service.after_copy_verified = None
            assert sealed.status_code == 200, sealed.text
            assert sha256(sealed.content).hexdigest() == fixture_sha
            assert sealed.content == fixture_bytes
            assert seen and not seen[0].exists()
            assert not list(
                (workspace / "registered-result-responses").glob(
                    "helixweave-rreg-response-*"
                )
            )

            changed = client.get(_content(project.project_id, revision))
            assert changed.status_code == 409, changed.text
            assert changed.json()["issues"][0]["code"] == (
                "registered_result_bytes_unavailable"
            )
            assert "synthetic-A" not in changed.text
            assert str(source) not in changed.text
            assert (
                _events_for(
                    app,
                    actor_user_id=member.user_id,
                    project_id=project.project_id,
                    result_id=revision["result_id"],
                    result_revision_id=revision["result_revision_id"],
                )
                == member_events
            )

            replay_after_change = client.post(
                f"/api/v1/catalog/projects/{project.project_id}/registered-results",
                json=metadata,
                headers={"Idempotency-Key": "irene39-register-1"},
            )
            assert replay_after_change.status_code == 200, replay_after_change.text
            assert replay_after_change.json()["revision"]["bytes"][
                "content_sha256"
            ] == (fixture_sha)
            _assert_runs_unchanged(app)
            assert FIXTURE.read_bytes() == fixture_bytes
    finally:
        app.state.run_queue.close()
        app.state.persistence.close()


def _member(app, admin):
    account = replace(
        admin,
        user_id="usr_" + "c" * 32,
        username="irene39-member",
        role=UserRole.MEMBER,
    )
    app.state.authentication_repository.create_account(account)
    secrets = new_session_secrets()
    app.state.authentication_repository.create_session(
        new_session_record(
            user_id=account.user_id,
            secrets=secrets,
            created_at=datetime.now(timezone.utc),
            lifetime=timedelta(hours=1),
        )
    )
    return SimpleNamespace(
        user_id=account.user_id,
        session_token=secrets.session_token,
        csrf_token=secrets.csrf_token,
    )


def _use(app, account) -> None:
    app.state.test_auth_tokens = (account.session_token, account.csrf_token)


def _use_admin(app) -> None:
    app.state.test_auth_tokens = (TEST_AUTH_SESSION_TOKEN, TEST_AUTH_CSRF_TOKEN)


def _content(project_id: str, revision: dict) -> str:
    return (
        f"/api/v1/catalog/projects/{project_id}/registered-result-revisions/"
        f"{revision['result_revision_id']}/content"
    )


def _revision(project_id: str, result_revision_id: str) -> str:
    return (
        f"/api/v1/catalog/projects/{project_id}/registered-result-revisions/"
        f"{result_revision_id}"
    )


def _revisions(project_id: str, result_id: str) -> str:
    return (
        f"/api/v1/catalog/projects/{project_id}/registered-results/"
        f"{result_id}/revisions"
    )


def _counts(body: bytes) -> dict[str, str]:
    rows = {}
    for line in body.decode("utf-8").splitlines():
        key, value = line.split("\t")
        rows[key] = value
    return rows


def _assert_download_headers(response) -> None:
    assert response.headers["content-type"].startswith("application/octet-stream")
    assert response.headers["accept-ranges"] == "none"
    assert response.headers["cache-control"] == "no-store"
    assert response.headers["x-content-type-options"] == "nosniff"
    assert (
        response.headers["content-disposition"]
        == 'attachment; filename="registered-result.tsv"'
    )
    assert response.headers["content-length"] == "50"


def _assert_public(revision: dict) -> None:
    assert set(revision) == PUBLIC_KEYS
    assert FORBIDDEN.isdisjoint(revision)
    assert "run_id" not in json.dumps(revision)


def _assert_runs_unchanged(app) -> None:
    with app.state.persistence.engine.connect() as connection:
        assert connection.scalar(text("SELECT count(*) FROM runs")) == 0


def _success_downloads(app) -> tuple:
    return tuple(
        event
        for event in app.state.authentication_repository.list_security_audit_events(
            limit=100
        )
        if event.action.value == "registered_result.download"
        and event.outcome.value == "succeeded"
    )


def _assert_binary_contract(app) -> None:
    operation = app.openapi()["paths"][
        "/api/v1/catalog/projects/{project_id}/registered-result-revisions/"
        "{result_revision_id}/content"
    ]["get"]
    content = operation["responses"]["200"]["content"]
    assert set(content) == {"application/octet-stream"}
    assert content["application/octet-stream"]["schema"] == {
        "type": "string",
        "format": "binary",
    }
    for status in ("400", "404", "409", "500"):
        assert "application/json" in operation["responses"][status]["content"]
    generated = (
        REPOSITORY_ROOT
        / "frontend/src/api/generated/registered-results/registered-results.ts"
    ).read_text(encoding="utf-8")
    download = generated.split(
        "export const catalogDownloadRegisteredResultRevision",
        1,
    )[1]
    assert "blobFetcher<Blob>" in download
    assert "unknown | Blob" not in download
    assert "return fetcher<" not in download


def _generated_client_download(
    app,
    *,
    project_id: str,
    result_revision_id: str,
    session_token: str,
    fixture_sha: str,
) -> None:
    server, thread, port = _start_api(app)
    try:
        cookie_name = app.state.auth_cookie_policy.session_cookie.name
        environment = os.environ.copy()
        environment.update(
            {
                "HELIXWEAVE_REGISTERED_RESULT_BASE": f"http://127.0.0.1:{port}",
                "HELIXWEAVE_REGISTERED_RESULT_PROJECT": project_id,
                "HELIXWEAVE_REGISTERED_RESULT_REVISION": result_revision_id,
                "HELIXWEAVE_REGISTERED_RESULT_COOKIE": f"{cookie_name}={session_token}",
                "HELIXWEAVE_REGISTERED_RESULT_SHA256": fixture_sha,
            }
        )
        completed = subprocess.run(
            [
                "npm",
                "exec",
                "--",
                "vitest",
                "run",
                "--config",
                "vitest.registered-result.config.ts",
                "--reporter=verbose",
            ],
            cwd=REPOSITORY_ROOT / "frontend",
            env=environment,
            capture_output=True,
            text=True,
            check=False,
            timeout=180,
        )
        assert completed.returncode == 0, completed.stdout + "\n" + completed.stderr
    finally:
        server.should_exit = True
        thread.join(timeout=10)


def _start_api(app):
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]

    class _ApiServer(uvicorn.Server):
        def install_signal_handlers(self) -> None:
            return None

    server = _ApiServer(
        uvicorn.Config(
            app,
            host="127.0.0.1",
            port=port,
            log_level="error",
            access_log=False,
            lifespan="off",
        )
    )
    thread = threading.Thread(target=server.run, name="irene39-registered-result")
    thread.start()
    deadline = time.monotonic() + 10
    while not server.started:
        if not thread.is_alive() or time.monotonic() > deadline:
            raise RuntimeError("registered result server did not start")
        time.sleep(0.05)
    return server, thread, port


def _download_with_unreadable_copy(app, client, project_id: str, revision: dict):
    with patch(
        "encode_pipeline.services.registered_results.os.open",
        side_effect=OSError("unreadable"),
    ):
        return client.get(_content(project_id, revision))


def _download_with_audit_failure(app, project_id: str, revision: dict):
    def fail(*args, **kwargs):
        raise RuntimeError("audit unavailable")

    async def once():
        transport = httpx.ASGITransport(app=app, raise_app_exceptions=False)
        async with seeded_auth_async_client(
            app,
            transport=transport,
            base_url="http://testserver",
        ) as client:
            return await client.get(_content(project_id, revision))

    with patch.object(
        app.state.authentication_repository,
        "record_security_audit",
        fail,
    ):
        return asyncio.run(once())


def _assert_no_response_copies(workspace: Path) -> None:
    root = workspace / "registered-result-responses"
    assert not list(root.glob("helixweave-rreg-response-*"))


def _events_for(app, *, actor_user_id, project_id, result_id, result_revision_id):
    resource = build_audit_resource(
        AuditResourceKind.REGISTERED_RESULT,
        project_id,
        result_id,
        result_revision_id,
    )
    return tuple(
        event
        for event in _success_downloads(app)
        if event.actor_user_id == actor_user_id
        and event.resource is not None
        and event.resource.resource_id == resource.resource_id
    )
