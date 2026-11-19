"""Branch coverage for registered-result reads, writes, and failure paths.

The acceptance test stays on the external happy path. These calls exercise
listing, archive, downstream parents, and the closed error branches that the
fast coverage shards otherwise never enter.
"""

from __future__ import annotations

import asyncio
from dataclasses import replace
from datetime import datetime, timezone
import json
from pathlib import Path
import stat
import threading
from types import SimpleNamespace
from unittest.mock import patch

import httpx

import pytest

from encode_pipeline.api.dependencies import AuthApiError
from encode_pipeline.api.main import create_app
from encode_pipeline.api.routes import registered_results as registered_result_routes
from encode_pipeline.api.routes.registered_results import _failure
from encode_pipeline.persistence.models import RegisteredResultRevisionRow
from encode_pipeline.persistence.registered_results import (
    RegisteredResultStoreError,
    RevisionDraft,
    _require_parent_pins,
)
from encode_pipeline.platform.authentication import (
    AuthenticatedPrincipal,
    UserRole,
)
from encode_pipeline.platform.registered_results import (
    ParentEdge,
    parse_metadata,
    request_sha256,
)
from encode_pipeline.services import registered_results as registered_result_service
from encode_pipeline.services.data_registry import DataRegistryService
from encode_pipeline.services.input_file_access import InputFileAccess
from encode_pipeline.services.input_registry import InputRegistryService
from encode_pipeline.services.registered_results import (
    RegisteredResultAPIError,
    _decode_cursor,
    _encode_cursor,
    _zulu,
)
from api_test_client import ApiTestClient, seeded_auth_async_client
from conftest import seed_test_authentication

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
FIXTURE = REPOSITORY_ROOT / "test/fixtures/irene39/fixture39.tsv"


def test_auth_errors_leave_the_registered_result_route() -> None:
    with pytest.raises(AuthApiError):
        _failure(RegisteredResultAPIError(401, "AUTHENTICATION_REQUIRED", "no"))
    with pytest.raises(AuthApiError):
        _failure(RegisteredResultAPIError(403, "ADMINISTRATOR_REQUIRED", "no"))


def test_registered_result_list_archive_downstream_and_failures(
    tmp_path: Path,
) -> None:
    fixture = FIXTURE.read_bytes()
    pool_root = (tmp_path / "pool").resolve()
    source = pool_root / "external" / "known.tsv"
    source.parent.mkdir(parents=True)
    source.write_bytes(fixture)
    config_path = tmp_path / "storage-pools.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": "storage-pool-roots-v1",
                "storage_pool_roots": {"branch-pool": str(pool_root)},
            }
        ),
        encoding="utf-8",
    )
    app = create_app(
        database_url=f"sqlite:///{tmp_path / 'branches.db'}",
        workspace_root=tmp_path / "workspace",
        project_root=REPOSITORY_ROOT,
    )
    try:
        admin = seed_test_authentication(app)
        app.state.registered_result_service._storage_pool_config_path = config_path
        member = replace(
            admin,
            user_id="usr_" + "c" * 32,
            username="branch-member",
            role=UserRole.MEMBER,
        )
        app.state.authentication_repository.create_account(member)
        registry = DataRegistryService(
            repository=app.state.persistence.data_registry_repository
        )
        project = registry.create_project("IRENE-39 branches")
        unbound = registry.create_project("IRENE-39 unbound")
        pools = InputRegistryService(
            repository=app.state.persistence.input_registry_repository
        )
        storage = pools.create_storage_pool(
            display_name="Branch pool",
            config_key="branch-pool",
        )
        pools.bind_project_storage_pool(project.project_id, storage.storage_pool_id)
        known = _metadata("synthetic-known-source", "Known synthetic counts")
        other = _metadata("synthetic-other-source", "Other synthetic counts")
        client = ApiTestClient(app)
        registered = _post(client, _results(project.project_id), known, "branch-known")
        assert registered.status_code == 200, registered.text
        known_revision = registered.json()["revision"]
        second = _post(client, _results(project.project_id), other, "branch-other")
        assert second.status_code == 200, second.text

        listed = client.get(_results(project.project_id), params={"limit": 1})
        assert listed.status_code == 200, listed.text
        assert listed.json()["next_cursor"]
        followed = client.get(
            _results(project.project_id),
            params={"limit": 1, "cursor": listed.json()["next_cursor"]},
        )
        assert followed.status_code == 200, followed.text
        assert followed.json()["next_cursor"] is None
        assert (
            client.get(
                _results(project.project_id),
                params={"cursor": "not-a-cursor"},
            ).status_code
            == 400
        )
        missing_project = "prj_" + "0" * 32
        assert client.get(_results(missing_project)).status_code == 404

        revised = _post(
            client,
            _revisions(project.project_id, known_revision["result_id"]),
            {
                "expected_revision_number": 1,
                "reference": {"status": "not_applicable"},
            },
            "branch-reference",
        )
        assert revised.status_code == 200, revised.text
        assert revised.json()["revision"]["revision_number"] == 2
        page = client.get(
            _revisions(project.project_id, known_revision["result_id"]),
            params={"limit": 1},
        )
        assert page.status_code == 200, page.text
        assert page.json()["next_cursor"]
        rest = client.get(
            _revisions(project.project_id, known_revision["result_id"]),
            params={"limit": 1, "cursor": page.json()["next_cursor"]},
        )
        assert rest.status_code == 200, rest.text
        assert (
            client.get(
                _revisions(project.project_id, "rreg_" + "0" * 32),
            ).status_code
            == 404
        )

        child_body = _metadata(
            "synthetic-downstream",
            "Downstream synthetic counts",
            origin={
                "kind": "downstream",
                "parents": [
                    {
                        "result_revision_id": known_revision["result_revision_id"],
                        "role": "input",
                    }
                ],
            },
        )
        child = _post(client, _results(project.project_id), child_body, "branch-child")
        assert child.status_code == 200, child.text
        child_revision = child.json()["revision"]
        assert child_revision["origin"]["kind"] == "downstream"
        read_back = client.get(
            _revision(project.project_id, child_revision["result_revision_id"])
        )
        assert read_back.status_code == 200, read_back.text
        parent_view = read_back.json()["revision"]["origin"]["parents"][0]
        assert parent_view["revision_status"] == "superseded"
        assert parent_view["result_revision_id"] == known_revision["result_revision_id"]

        missing_parent = _metadata(
            "synthetic-missing-parent",
            "Missing parent",
            origin={
                "kind": "downstream",
                "parents": [
                    {"result_revision_id": "rrev_" + "f" * 32, "role": "input"}
                ],
            },
        )
        assert (
            _post(
                client,
                _results(project.project_id),
                missing_parent,
                "branch-missing-parent",
            ).status_code
            == 400
        )
        assert (
            _post(
                client,
                _results(project.project_id),
                known,
                "branch-duplicate",
            ).status_code
            == 400
        )
        invalid = _metadata("synthetic-invalid", "bad\nname")
        assert (
            _post(
                client,
                _results(project.project_id),
                invalid,
                "branch-invalid",
            ).status_code
            == 400
        )
        assert (
            _post(
                client,
                _results(unbound.project_id),
                other,
                "branch-unbound",
            ).status_code
            == 404
        )
        stale_cursor = _encode_cursor(
            {
                "include_archived": False,
                "kind": "results",
                "project_id": project.project_id,
                "result_id": "rreg_" + "0" * 32,
                "stable_key": "missing",
            }
        )
        assert (
            client.get(
                _results(project.project_id),
                params={"cursor": stale_cursor},
            ).status_code
            == 404
        )
        bad_number = _encode_cursor(
            {
                "kind": "revisions",
                "project_id": project.project_id,
                "result_id": known_revision["result_id"],
                "revision_number": True,
            }
        )
        assert (
            client.get(
                _revisions(project.project_id, known_revision["result_id"]),
                params={"cursor": bad_number},
            ).status_code
            == 400
        )

        _repository_edges(
            app,
            admin=admin,
            member=member,
            project_id=project.project_id,
            result_id=known_revision["result_id"],
            revision_id=known_revision["result_revision_id"],
            known=known,
        )
        archived = _post(
            client,
            _archive(project.project_id, known_revision["result_id"]),
            {"expected_revision_number": 2},
            "branch-archive",
        )
        assert archived.status_code == 200, archived.text
        replayed = _post(
            client,
            _archive(project.project_id, known_revision["result_id"]),
            {"expected_revision_number": 2},
            "branch-archive",
        )
        assert replayed.status_code == 200, replayed.text
        assert (
            _post(
                client,
                _revisions(project.project_id, known_revision["result_id"]),
                {
                    "expected_revision_number": 2,
                    "display_name": "archived",
                },
                "branch-archived-revise",
            ).status_code
            == 409
        )
        visible = client.get(_results(project.project_id))
        archived_list = client.get(
            _results(project.project_id),
            params={"include_archived": True},
        )
        assert visible.status_code == 200, visible.text
        assert archived_list.status_code == 200, archived_list.text
        visible_names = {item["display_name"] for item in visible.json()["results"]}
        archived_names = {
            item["display_name"] for item in archived_list.json()["results"]
        }
        assert "Known synthetic counts" not in visible_names
        assert any(name.startswith("Known synthetic") for name in archived_names)
        assert (
            client.get(_revision(project.project_id, "rrev_" + "0" * 32)).status_code
            == 404
        )

        _direct_failures(
            app,
            admin=admin,
            project_id=project.project_id,
            result_id=known_revision["result_id"],
            revision_id=known_revision["result_revision_id"],
            source=source,
            fixture=fixture,
            client=client,
        )
    finally:
        app.state.run_queue.close()
        app.state.persistence.close()


def test_slow_registered_result_file_io_leaves_the_event_loop_free(
    tmp_path: Path,
) -> None:
    fixture = FIXTURE.read_bytes()
    pool_root = (tmp_path / "pool").resolve()
    source = pool_root / "external" / "known.tsv"
    source.parent.mkdir(parents=True)
    source.write_bytes(fixture)
    config_path = tmp_path / "storage-pools.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": "storage-pool-roots-v1",
                "storage_pool_roots": {"loop-pool": str(pool_root)},
            }
        ),
        encoding="utf-8",
    )
    app = create_app(
        database_url=f"sqlite:///{tmp_path / 'loop.db'}",
        workspace_root=tmp_path / "workspace",
        project_root=REPOSITORY_ROOT,
    )
    try:
        seed_test_authentication(app)
        app.state.registered_result_service._storage_pool_config_path = config_path
        registry = DataRegistryService(
            repository=app.state.persistence.data_registry_repository
        )
        project = registry.create_project("IRENE-39 loop")
        pools = InputRegistryService(
            repository=app.state.persistence.input_registry_repository
        )
        storage = pools.create_storage_pool(
            display_name="Loop pool",
            config_key="loop-pool",
        )
        pools.bind_project_storage_pool(project.project_id, storage.storage_pool_id)
        asyncio.run(_prove_file_io_stays_off_the_loop(app, project.project_id, fixture))
    finally:
        app.state.run_queue.close()
        app.state.persistence.close()


async def _prove_file_io_stays_off_the_loop(
    app, project_id: str, fixture: bytes
) -> None:
    mode = {"name": ""}
    entered = threading.Event()
    release = threading.Event()
    original_observe = InputFileAccess.observe
    original_seal = registered_result_service.seal_private_copy
    real_os = registered_result_routes.os

    def blocking_observe(self, relative_path: str):
        if mode["name"] == "observe":
            _wait_until_another_request_finishes(entered, release)
        return original_observe(self, relative_path)

    def blocking_seal(**kwargs):
        if mode["name"] == "seal":
            _wait_until_another_request_finishes(entered, release)
        return original_seal(**kwargs)

    class _RouteOs:
        def __getattr__(self, name: str):
            return getattr(real_os, name)

        def read(self, descriptor: int, length: int) -> bytes:
            if mode["name"] == "stream":
                _wait_until_another_request_finishes(entered, release)
            return real_os.read(descriptor, length)

    list_url = _results(project_id)
    transport = httpx.ASGITransport(app=app)
    async with seeded_auth_async_client(
        app,
        transport=transport,
        base_url="http://testserver",
    ) as client:

        async def during(name: str, start):
            entered.clear()
            release.clear()
            mode["name"] = name
            task = asyncio.create_task(start())
            assert await asyncio.to_thread(entered.wait, 5)
            listed = await client.get(list_url)
            assert listed.status_code == 200, listed.text
            mode["name"] = ""
            release.set()
            response = await task
            assert response.status_code == 200, response.text
            return response

        with (
            patch.object(InputFileAccess, "observe", blocking_observe),
            patch.object(
                registered_result_service,
                "seal_private_copy",
                blocking_seal,
            ),
            patch.object(registered_result_routes, "os", _RouteOs()),
        ):
            registered = await during(
                "observe",
                lambda: client.post(
                    list_url,
                    json=_metadata("synthetic-loop-source", "Loop synthetic counts"),
                    headers={"Idempotency-Key": "loop-register"},
                ),
            )
            revision = registered.json()["revision"]
            await during(
                "observe",
                lambda: client.post(
                    _revisions(project_id, revision["result_id"]),
                    json={
                        "expected_revision_number": 1,
                        "display_name": "Loop renamed counts",
                    },
                    headers={"Idempotency-Key": "loop-rename"},
                ),
            )
            content = _content(project_id, revision["result_revision_id"])
            sealed = await during("seal", lambda: client.get(content))
            assert sealed.content == fixture
            streamed = await during("stream", lambda: client.get(content))
            assert streamed.content == fixture


def _wait_until_another_request_finishes(
    entered: threading.Event,
    release: threading.Event,
) -> None:
    entered.set()
    if not release.wait(3):
        raise TimeoutError("event loop stayed blocked during file read")


def _repository_edges(
    app,
    *,
    admin,
    member,
    project_id: str,
    result_id: str,
    revision_id: str,
    known: dict,
) -> None:
    repository = app.state.registered_result_service._repository
    record = repository.get_head(result_id)
    pinned = repository.get_revision(revision_id)
    assert record is not None and pinned is not None
    metadata = parse_metadata(
        {
            key: known[key]
            for key in ("display_name", "format", "origin", "method", "reference")
        }
    )
    body = {
        "display_name": metadata.display_name,
        "format": metadata.format,
        "method": metadata.method,
        "origin": metadata.origin,
        "pool_relative_path": known["pool_relative_path"],
        "reference": metadata.reference,
        "stable_key": known["stable_key"],
    }
    replay, replayed = repository.register(
        actor_user_id=admin.user_id,
        idempotency_key="branch-known",
        request_sha256=request_sha256("register", project_id, None, body),
        project_id=project_id,
        stable_key=known["stable_key"],
        draft=_draft(),
    )
    assert replayed is True
    assert replay.result_id == result_id
    with pytest.raises(RegisteredResultStoreError) as missing_project:
        repository.register(
            actor_user_id=admin.user_id,
            idempotency_key="branch-no-project",
            request_sha256="11" * 32,
            project_id="prj_" + "0" * 32,
            stable_key="branch-no-project",
            draft=_draft(),
        )
    assert missing_project.value.code == "not_found"
    with pytest.raises(RegisteredResultStoreError) as bad_number:
        repository.register(
            actor_user_id=admin.user_id,
            idempotency_key="branch-bad-number",
            request_sha256="22" * 32,
            project_id=project_id,
            stable_key="branch-bad-number",
            draft=_draft(revision_number=2),
        )
    assert bad_number.value.code == "payload"
    with pytest.raises(RegisteredResultStoreError) as unknown_user:
        repository.register(
            actor_user_id="usr_" + "9" * 32,
            idempotency_key="branch-unknown-user",
            request_sha256="33" * 32,
            project_id=project_id,
            stable_key="branch-unknown-user",
            draft=_draft(),
        )
    assert unknown_user.value.code == "unauthenticated"
    with pytest.raises(RegisteredResultStoreError) as forbidden:
        repository.register(
            actor_user_id=member.user_id,
            idempotency_key="branch-member",
            request_sha256="44" * 32,
            project_id=project_id,
            stable_key="branch-member",
            draft=_draft(),
        )
    assert forbidden.value.code == "forbidden"
    with pytest.raises(RegisteredResultStoreError) as stale:
        repository.revise(
            actor_user_id=admin.user_id,
            idempotency_key="branch-stale-direct",
            request_sha256="55" * 32,
            project_id=project_id,
            expected_revision_number=1,
            draft=_draft(result_id=result_id, revision_number=9),
        )
    assert stale.value.code == "revision"
    with pytest.raises(RegisteredResultStoreError) as wrong_archive:
        repository.archive(
            actor_user_id=admin.user_id,
            idempotency_key="branch-archive-wrong",
            request_sha256="66" * 32,
            project_id=project_id,
            result_id=result_id,
            expected_revision_number=9,
            archived_at=datetime.now(timezone.utc),
        )
    assert wrong_archive.value.code == "revision"
    assert repository.get_head("rreg_" + "0" * 32) is None
    assert repository.get_revision("rrev_" + "0" * 32) is None
    assert repository.visible_project("prj_" + "0" * 32) is False
    assert repository.result_in_project(project_id, "rreg_" + "0" * 32) is False
    assert (
        repository.current_anchor_matches(
            project_id,
            "rreg_" + "0" * 32,
            include_archived=False,
        )
        is False
    )
    assert repository.revision_anchor_matches(project_id, result_id, 1) is True
    with repository._session_factory() as session:
        with pytest.raises(RegisteredResultStoreError):
            _require_parent_pins(
                session,
                project_id,
                (
                    ParentEdge(
                        result_revision_id=revision_id,
                        role="input",
                        revision_digest="ab" * 32,
                    ),
                ),
            )
        with pytest.raises(RegisteredResultStoreError):
            _require_parent_pins(
                session,
                "prj_" + "0" * 32,
                (
                    ParentEdge(
                        result_revision_id=revision_id,
                        role="input",
                        revision_digest=pinned.revision_digest,
                    ),
                ),
            )
    with pytest.raises(RegisteredResultAPIError):
        app.state.registered_result_service._archive_body(record)


def _direct_failures(
    app,
    *,
    admin,
    project_id: str,
    result_id: str,
    revision_id: str,
    source: Path,
    fixture: bytes,
    client: ApiTestClient,
) -> None:
    service = app.state.registered_result_service
    repository = service._repository
    principal = AuthenticatedPrincipal.from_account(admin)
    with pytest.raises(RegisteredResultAPIError):
        service._limit(0)
    with pytest.raises(RegisteredResultAPIError):
        service._idempotency(" bad")
    with pytest.raises(RegisteredResultAPIError):
        service._revision_number(True)  # type: ignore[arg-type]
    with pytest.raises(RegisteredResultAPIError):
        service._project("project")
    with pytest.raises(RegisteredResultAPIError):
        service._result("result")
    with pytest.raises(RegisteredResultAPIError):
        _decode_cursor("rcur1." + ("a" * 5000))
    with pytest.raises(RegisteredResultAPIError):
        _decode_cursor("rcur1." + "!!!!")
    assert _zulu(datetime(2026, 11, 19, 1, 2, 3)).endswith("Z")

    record = repository.get_revision(revision_id)
    assert record is not None
    with pytest.raises(RegisteredResultAPIError):
        service._stored_metadata(replace(record, canonical_payload="{"))
    with pytest.raises(RegisteredResultAPIError):
        service._require_digest(
            replace(record, revision_digest="cd" * 32),
            service._stored_metadata(record),
        )
    metadata = parse_metadata(
        {
            "display_name": "Child",
            "format": "tsv",
            "origin": {
                "kind": "downstream",
                "parents": [
                    {"result_revision_id": "rrev_" + "e" * 32, "role": "input"}
                ],
            },
            "method": {"status": "unknown"},
            "reference": {"status": "unknown"},
        }
    )
    with pytest.raises(RegisteredResultAPIError):
        service._parent_edges(project_id, metadata)
    profiles = service._reference_profiles

    class _MissingProfile:
        def get_revision(self, revision_id: str) -> None:
            raise KeyError(revision_id)

    class _MismatchedProfile:
        def get_revision(self, revision_id: str):
            return SimpleNamespace(public_identity_sha256="ab" * 32)

    try:
        service._reference_profiles = _MissingProfile()
        with pytest.raises(RegisteredResultAPIError):
            service._check_reference(
                {
                    "kind": "reference_profile_revision",
                    "reference_profile_revision_id": "refpr_" + "a" * 32,
                    "public_identity_sha256": "cd" * 32,
                }
            )
        service._reference_profiles = _MismatchedProfile()
        with pytest.raises(RegisteredResultAPIError):
            service._check_reference(
                {
                    "kind": "reference_profile_revision",
                    "reference_profile_revision_id": "refpr_" + "a" * 32,
                    "public_identity_sha256": "cd" * 32,
                }
            )
        service._check_reference({"kind": "external", "status": "known"})
    finally:
        service._reference_profiles = profiles

    with pytest.raises(RegisteredResultAPIError):
        service._observe(record.storage_pool_id, "missing.tsv")
    with pytest.raises(RegisteredResultAPIError):
        service._pool_config_key("stgp_" + "0" * 32)
    saved_config = service._storage_pool_config_path
    try:
        service._storage_pool_config_path = None
        with pytest.raises(RegisteredResultAPIError):
            service._pool_config()
        service._storage_pool_config_path = source.parent / "missing.json"
        with pytest.raises(RegisteredResultAPIError):
            service._pool_config()
    finally:
        service._storage_pool_config_path = saved_config
    saved_registry = service._input_registry

    class _BrokenRegistry:
        def get_active_project_storage_pool_binding(self, project_id: str) -> None:
            raise ValueError(project_id)

        def get_storage_pool(self, storage_pool_id: str) -> None:
            raise KeyError(storage_pool_id)

    try:
        service._input_registry = _BrokenRegistry()
        with pytest.raises(RegisteredResultAPIError):
            service._active_binding(project_id)
    finally:
        service._input_registry = saved_registry
    with patch.object(Path, "mkdir", side_effect=OSError("unavailable")):
        with pytest.raises(RegisteredResultAPIError):
            service._ensure_response_root()
    saved_auth = service._authentication
    try:
        service._authentication = SimpleNamespace(
            get_account_by_id=lambda user_id: (_ for _ in ()).throw(KeyError(user_id))
        )
        with pytest.raises(RegisteredResultAPIError) as missing_account:
            service.list_current(
                principal,
                project_id,
                cursor=None,
                limit=1,
                include_archived=False,
            )
        assert missing_account.value.status_code == 401
        service._authentication = SimpleNamespace(
            get_account_by_id=lambda user_id: SimpleNamespace(enabled=False)
        )
        with pytest.raises(RegisteredResultAPIError) as disabled:
            service.list_current(
                principal,
                project_id,
                cursor=None,
                limit=1,
                include_archived=False,
            )
        assert disabled.value.status_code == 401
    finally:
        service._authentication = saved_auth

    original = None
    try:
        with repository._session_factory.begin() as session:
            row = session.get(RegisteredResultRevisionRow, revision_id)
            assert row is not None
            original = row.canonical_payload
            row.canonical_payload = "{"
        with pytest.raises(RegisteredResultStoreError) as broken:
            repository.get_revision(revision_id)
        assert broken.value.code == "unavailable"
    finally:
        if original is not None:
            with repository._session_factory.begin() as session:
                row = session.get(RegisteredResultRevisionRow, revision_id)
                assert row is not None
                row.canonical_payload = original

    source.write_bytes(fixture)
    content = _content(project_id, revision_id)

    def bad_stat(_descriptor: int):
        return SimpleNamespace(st_mode=stat.S_IFREG, st_size=0)

    with patch(
        "encode_pipeline.services.registered_results.os.fstat",
        side_effect=bad_stat,
    ):
        failed_stat = client.get(content)
    assert failed_stat.status_code == 409, failed_stat.text

    def boom(_path: Path) -> None:
        raise RuntimeError("hook failed")

    service.after_copy_verified = boom
    try:
        with pytest.raises(RuntimeError, match="hook failed"):
            client.get(content)
    finally:
        service.after_copy_verified = None
    with (
        patch(
            "encode_pipeline.api.routes.registered_results.StreamingResponse",
            side_effect=RuntimeError("response failed"),
        ),
        pytest.raises(RuntimeError, match="response failed"),
    ):
        client.get(content)


def _metadata(stable_key: str, display_name: str, **overrides: object) -> dict:
    payload: dict = {
        "stable_key": stable_key,
        "display_name": display_name,
        "format": "tsv",
        "origin": {"kind": "external", "source": {"status": "unknown"}},
        "method": {"status": "unknown"},
        "reference": {"status": "unknown"},
        "pool_relative_path": "external/known.tsv",
    }
    payload.update(overrides)
    return payload


def _post(client: ApiTestClient, url: str, body: dict, key: str):
    return client.post(url, json=body, headers={"Idempotency-Key": key})


def _results(project_id: str) -> str:
    return f"/api/v1/catalog/projects/{project_id}/registered-results"


def _revisions(project_id: str, result_id: str) -> str:
    return f"{_results(project_id)}/{result_id}/revisions"


def _archive(project_id: str, result_id: str) -> str:
    return f"{_results(project_id)}/{result_id}/archive"


def _revision(project_id: str, result_revision_id: str) -> str:
    return (
        "/api/v1/catalog/projects/"
        f"{project_id}/registered-result-revisions/{result_revision_id}"
    )


def _content(project_id: str, result_revision_id: str) -> str:
    return f"{_revision(project_id, result_revision_id)}/content"


def _draft(**overrides: object) -> RevisionDraft:
    values: dict = {
        "result_id": "rreg_" + "1" * 32,
        "result_revision_id": "rrev_" + "1" * 32,
        "revision_number": 1,
        "canonical_payload": "{}",
        "revision_digest": "ab" * 32,
        "created_at": datetime.now(timezone.utc),
        "storage_pool_id": "stgp_" + "1" * 32,
        "relative_path": "external/known.tsv",
        "size_bytes": 50,
        "content_sha256": "cd" * 32,
        "parents": (),
    }
    values.update(overrides)
    return RevisionDraft(**values)
