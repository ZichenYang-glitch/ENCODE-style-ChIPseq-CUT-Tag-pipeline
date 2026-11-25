"""HTTP contract for sealed track reads."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from hashlib import sha256
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

pytest.importorskip("fastapi")

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from api_test_client import ApiTestClient  # noqa: E402
from conftest import (  # noqa: E402
    TEST_AUTH_PASSWORD_HASH,
    TEST_AUTH_USER_ID,
    seed_test_authentication,
)
from encode_pipeline.api.main import create_app  # noqa: E402
from encode_pipeline.persistence.database import create_session_factory  # noqa: E402
from encode_pipeline.persistence.models import (  # noqa: E402
    ExperimentRevisionRow,
    ExperimentRow,
    RegisteredResultRevisionRow,
    RegisteredResultRow,
    RunExperimentBindingRow,
    RunProjectBindingRow,
    RunRow,
    StoragePoolRow,
)
from encode_pipeline.platform.authentication import (  # noqa: E402
    AuthenticatedPrincipal,
    UserAccount,
    UserRole,
    UserStatus,
)
from encode_pipeline.services.archived_tracks import (  # noqa: E402
    ArchivedTrackService,
    SqlArchiveReadGate,
)
from encode_pipeline.services.authentication import (  # noqa: E402
    SessionSecrets,
    new_session_record,
)
from encode_pipeline.services.authentication_service import (  # noqa: E402
    AuthenticationActor,
)
from encode_pipeline.services.data_registry import DataRegistryService  # noqa: E402


_HEX = "ab" * 32
_GEN = "artifactgen-" + "a" * 64
_REV = "artifactrev-" + "b" * 64
_EXPERIMENT = "exp_" + "a" * 32
_REVISION = "expr_" + "b" * 32
_OTHER_EXPERIMENT = "exp_" + "c" * 32
_OTHER_REVISION = "expr_" + "d" * 32
_RESULT = "rreg_" + "1" * 32
_RESULT_REVISION = "rrev_" + "2" * 32
_OTHER_RESULT = "rreg_" + "3" * 32
_OTHER_RESULT_REVISION = "rrev_" + "4" * 32
_POOL = "stgp_" + "5" * 32
_RUN = "run-archive-v1"
_V1 = "aver_" + "1" * 32
_V2 = "aver_" + "2" * 32
_V_UNBOUND = "aver_" + "3" * 32
_V_REGISTERED = "aver_" + "4" * 32
_V_FOREIGN = "aver_" + "5" * 32
_V_PAIR = "aver_" + "6" * 32
_V_LIST = "aver_" + "7" * 32
_SECRET_PATH = "pool/secret-signal.bw"
_SECRET_SHA = "cd" * 32


class ArchiveStoreError(Exception):
    def __init__(self, code: str) -> None:
        self.code = code
        super().__init__(code)


class _Member:
    def __init__(self, path: Path, payload: bytes) -> None:
        path.write_bytes(payload)
        self.path = path
        self.size_bytes = len(payload)
        self.sha256 = sha256(payload).hexdigest()


class _Reader:
    def __init__(self, payload: bytes) -> None:
        self._payload = payload
        self.reads: list[tuple[int, int]] = []

    def read(self, start: int, end: int) -> bytes:
        self.reads.append((start, end))
        return self._payload[start:end]

    def close(self) -> None:
        return None

    def __enter__(self) -> "_Reader":
        return self

    def __exit__(self, *_args: object) -> bool:
        self.close()
        return False


class _Objects:
    def __init__(self) -> None:
        self.readers: list[_Reader] = []
        self.attempts: list[_Member] = []

    def open_member(self, member: _Member) -> _Reader:
        self.attempts.append(member)
        if not member.path.is_file():
            raise ArchiveStoreError("archive_offline")
        payload = member.path.read_bytes()
        if (
            len(payload) != member.size_bytes
            or sha256(payload).hexdigest() != member.sha256
        ):
            raise ArchiveStoreError("archive_integrity")
        reader = _Reader(payload)
        self.readers.append(reader)
        return reader


class _Repository:
    def __init__(self) -> None:
        self.rows: list[SimpleNamespace] = []

    def list_versions(self, project_id, experiment_id, *, limit, offset):
        matched = [
            row
            for row in self.rows
            if row.project_id == project_id and row.experiment_id == experiment_id
        ]
        matched.sort(key=lambda row: (row.version_number, row.archive_version_id))
        page = matched[offset : offset + limit]
        next_offset = offset + limit if offset + limit < len(matched) else None
        return page, next_offset

    def get_version(self, project_id, experiment_id, archive_version_id):
        for row in self.rows:
            if (
                row.project_id == project_id
                and row.experiment_id == experiment_id
                and row.archive_version_id == archive_version_id
            ):
                return row
        raise ArchiveStoreError("archive_version_not_found")


def _version(
    *,
    archive_version_id: str,
    project_id: str,
    experiment_id: str,
    experiment_revision_id: str,
    version_number: int,
    source: dict,
    members: dict,
    display_name: str = "signal.bw",
    track_type: str = "bigwig",
    assembly: str = "hg38",
    has_index: bool | None = None,
) -> SimpleNamespace:
    return SimpleNamespace(
        archive_version_id=archive_version_id,
        project_id=project_id,
        experiment_id=experiment_id,
        experiment_revision_id=experiment_revision_id,
        track_key="signal",
        version_number=version_number,
        display_name=display_name,
        assembly=assembly,
        track_type=track_type,
        source=source,
        members=members,
        has_index=has_index,
    )


@pytest.fixture
def app(tmp_path):
    application = create_app(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        workspace_root=tmp_path / "workspaces",
    )
    yield application
    application.state.persistence.close()


@pytest.fixture
def seeded_app(app):
    seed_test_authentication(app)
    return app


def _registry(app) -> DataRegistryService:
    return DataRegistryService(repository=app.state.persistence.data_registry_repository)


def _now() -> datetime:
    return datetime.now(timezone.utc)


def _install(app, tmp_path: Path):
    project = _registry(app).create_project("Archive project")
    other = _registry(app).create_project("Other project")
    moment = _now()
    factory = create_session_factory(app.state.persistence.engine)

    def _commit(*rows) -> None:
        with factory() as session:
            session.add_all(rows)
            session.commit()

    _commit(
        ExperimentRow(
            experiment_id=_EXPERIMENT,
            project_id=project.project_id,
            stable_key="histone",
            created_at=moment,
        ),
        ExperimentRow(
            experiment_id=_OTHER_EXPERIMENT,
            project_id=other.project_id,
            stable_key="histone",
            created_at=moment,
        ),
    )
    _commit(
        ExperimentRevisionRow(
            experiment_revision_id=_REVISION,
            experiment_id=_EXPERIMENT,
            project_id=project.project_id,
            revision_number=1,
            canonical_payload="{}",
            payload_digest_scheme="sha256-framed-experiment-revision-payload-v1",
            payload_digest="a" * 64,
            created_at=moment,
        ),
        ExperimentRevisionRow(
            experiment_revision_id=_OTHER_REVISION,
            experiment_id=_OTHER_EXPERIMENT,
            project_id=other.project_id,
            revision_number=1,
            canonical_payload="{}",
            payload_digest_scheme="sha256-framed-experiment-revision-payload-v1",
            payload_digest="b" * 64,
            created_at=moment,
        ),
    )
    _commit(
        RunRow(
            run_id=_RUN,
            workflow_id="chipseq",
            inputs={},
            status="succeeded",
            created_at=moment,
            updated_at=moment,
            tags={},
        )
    )
    _commit(
        RunExperimentBindingRow(
            run_id=_RUN,
            experiment_revision_id=_REVISION,
            payload_digest="c" * 64,
        ),
        RunProjectBindingRow(
            run_id=_RUN,
            project_id=project.project_id,
            binding_mode="bound_v1",
            provenance="resolved",
            workflow_inputs_digest="f" * 64,
            binding_digest_scheme="sha256-framed-project-sample-binding-v1",
            binding_digest="1" * 64,
            created_at=moment,
        ),
        StoragePoolRow(
            storage_pool_id=_POOL,
            config_key="archive-test-pool",
            display_name="Archive test pool",
            created_at=moment,
        ),
        RegisteredResultRow(
            result_id=_RESULT,
            project_id=project.project_id,
            stable_key="registered-signal",
            created_at=moment,
        ),
        RegisteredResultRow(
            result_id=_OTHER_RESULT,
            project_id=other.project_id,
            stable_key="registered-signal",
            created_at=moment,
        ),
    )
    _commit(
        RegisteredResultRevisionRow(
            result_revision_id=_RESULT_REVISION,
            result_id=_RESULT,
            revision_number=1,
            canonical_payload="{}",
            digest_scheme="sha256-framed-registered-result-revision-v1",
            revision_digest="d" * 64,
            created_at=moment,
            storage_pool_id=_POOL,
            relative_path=_SECRET_PATH,
            size_bytes=4,
            content_sha256=_SECRET_SHA,
        ),
        RegisteredResultRevisionRow(
            result_revision_id=_OTHER_RESULT_REVISION,
            result_id=_OTHER_RESULT,
            revision_number=1,
            canonical_payload="{}",
            digest_scheme="sha256-framed-registered-result-revision-v1",
            revision_digest="e" * 64,
            created_at=moment,
            storage_pool_id=_POOL,
            relative_path=_SECRET_PATH,
            size_bytes=4,
            content_sha256=_SECRET_SHA,
        ),
    )
    objects = _Objects()
    repository = _Repository()
    v1 = _Member(tmp_path / "v1.bw", b"0123456789abcdef")
    v2 = _Member(tmp_path / "v2.bw", b"fedcba9876543210")
    registered = _Member(tmp_path / "registered.bw", b"reg!")
    index = _Member(tmp_path / "pairs.tbi", b"idx-window")
    data = _Member(tmp_path / "pairs.bed", b"data-win")
    run_source = {
        "kind": "run",
        "run_id": _RUN,
        "artifact_id": "signal",
        "generation": _GEN,
        "revision": _REV,
        "path": "/var/archive/secret.bw",
        "sha256": _SECRET_SHA,
    }
    repository.rows.extend(
        [
            _version(
                archive_version_id=_V1,
                project_id=project.project_id,
                experiment_id=_EXPERIMENT,
                experiment_revision_id=_REVISION,
                version_number=1,
                source=run_source,
                members={"data": v1},
            ),
            _version(
                archive_version_id=_V2,
                project_id=project.project_id,
                experiment_id=_EXPERIMENT,
                experiment_revision_id=_REVISION,
                version_number=2,
                source=dict(run_source),
                members={"data": v2},
            ),
            _version(
                archive_version_id=_V_UNBOUND,
                project_id=project.project_id,
                experiment_id=_EXPERIMENT,
                experiment_revision_id=_REVISION,
                version_number=3,
                source={
                    **run_source,
                    "run_id": "run-not-bound",
                },
                members={"data": v2},
            ),
            _version(
                archive_version_id=_V_REGISTERED,
                project_id=project.project_id,
                experiment_id=_EXPERIMENT,
                experiment_revision_id=_REVISION,
                version_number=4,
                source={
                    "kind": "registered",
                    "result_id": _RESULT,
                    "result_revision_id": _RESULT_REVISION,
                    "relative_path": _SECRET_PATH,
                },
                members={"data": registered},
                display_name="registered.bw",
            ),
            _version(
                archive_version_id=_V_FOREIGN,
                project_id=project.project_id,
                experiment_id=_EXPERIMENT,
                experiment_revision_id=_REVISION,
                version_number=5,
                source={
                    "kind": "registered",
                    "result_id": _OTHER_RESULT,
                    "result_revision_id": _OTHER_RESULT_REVISION,
                },
                members={"data": registered},
            ),
            _version(
                archive_version_id=_V_PAIR,
                project_id=project.project_id,
                experiment_id=_EXPERIMENT,
                experiment_revision_id=_REVISION,
                version_number=6,
                source={
                    **run_source,
                    "index_artifact_id": "signalIndex",
                    "index_generation": _GEN,
                    "index_revision": _REV,
                },
                members={"data": data, "index": index},
                display_name="pairs.bedpe",
                track_type="longrange",
            ),
            _version(
                archive_version_id=_V_LIST,
                project_id=project.project_id,
                experiment_id=_EXPERIMENT,
                experiment_revision_id=_REVISION,
                version_number=7,
                source=dict(run_source),
                members={},
                has_index=False,
            ),
        ]
    )
    app.state.archived_track_service = ArchivedTrackService(
        repository=repository,
        objects=objects,
        authentication=app.state.authentication_repository,
        access=SqlArchiveReadGate(factory),
    )
    return SimpleNamespace(
        project=project,
        other=other,
        objects=objects,
        repository=repository,
        v1=v1,
        v2=v2,
        index=index,
    )


def _collection(world) -> str:
    project_id = world.project.project_id
    return (
        "/api/v1/catalog/projects/"
        f"{project_id}/experiments/{_EXPERIMENT}/archived-tracks"
    )


def _content(world, archive_version_id: str, index: bool = False) -> str:
    suffix = "/content.tbi" if index else "/content"
    return f"{_collection(world)}/{archive_version_id}{suffix}"


def _member_session(app):
    repository = app.state.authentication_repository
    moment = _now()
    account = UserAccount(
        user_id="usr_" + "1" * 32,
        username="archive-member",
        role=UserRole.MEMBER,
        status=UserStatus.ENABLED,
        password_hash=TEST_AUTH_PASSWORD_HASH,
        created_at=moment,
        updated_at=moment,
        password_changed_at=moment,
    )
    repository.create_account(account)
    token = "c" * 43
    csrf = "d" * 43
    repository.create_session(
        new_session_record(
            user_id=account.user_id,
            secrets=SessionSecrets._from_generated(token, csrf),
            created_at=moment,
            lifetime=timedelta(hours=8),
        )
    )
    return account, (token, csrf)


def _no_store(response) -> None:
    assert response.headers["cache-control"] == "private, no-store"


def _redacted(response, *secrets: str) -> None:
    text = response.text
    for secret in secrets:
        assert secret not in text


def test_unconfigured_store_and_missing_session_fail_closed(seeded_app) -> None:
    project = "prj_" + "9" * 32
    experiment = "exp_" + "9" * 32
    url = (
        "/api/v1/catalog/projects/"
        f"{project}/experiments/{experiment}/archived-tracks"
    )
    with ApiTestClient(seeded_app) as client:
        offline = client.get(url)
        assert offline.status_code == 409
        assert offline.json()["issues"][0]["code"] == "archive_offline"
        _no_store(offline)
        assert "items" not in offline.text
        seeded_app.state.test_auth_tokens = None
        anonymous = client.get(url)
        assert anonymous.status_code == 401
        assert anonymous.json()["issues"][0]["code"] == "AUTHENTICATION_REQUIRED"
        _no_store(anonymous)


def test_setup_required_does_not_serve_bytes(app) -> None:
    url = (
        "/api/v1/catalog/projects/prj_"
        + "9" * 32
        + "/experiments/exp_"
        + "9" * 32
        + "/archived-tracks/"
        + _V1
        + "/content"
    )
    with ApiTestClient(app) as client:
        response = client.get(url)
    assert response.status_code == 503
    assert response.json()["issues"][0]["code"] == "SETUP_REQUIRED"
    _no_store(response)


def test_exact_version_ranges_and_same_size_replacement(seeded_app, tmp_path) -> None:
    world = _install(seeded_app, tmp_path)
    with ApiTestClient(seeded_app) as client:
        listed = client.get(_collection(world))
        assert listed.status_code == 200, listed.text
        _no_store(listed)
        body = listed.json()
        assert [item["archive_version_id"] for item in body["items"]] == [
            _V1,
            _V2,
            _V_REGISTERED,
            _V_PAIR,
            _V_LIST,
        ]
        assert body["next_offset"] is None
        first = body["items"][0]
        assert first["version_number"] == 1
        assert first["display_name"] == "signal.bw"
        assert first["assembly"] == "hg38"
        assert first["track_type"] == "bigwig"
        assert first["has_index"] is False
        listed_only = body["items"][-1]
        assert listed_only["archive_version_id"] == _V_LIST
        assert listed_only["has_index"] is False
        assert first["source"] == {
            "kind": "run",
            "run_id": _RUN,
            "artifact_id": "signal",
            "generation": _GEN,
            "revision": _REV,
        }
        _redacted(listed, "/var/archive", _SECRET_SHA, _SECRET_PATH)
        detail = client.get(f"{_collection(world)}/{_V1}")
        assert detail.status_code == 200
        assert detail.json() == first
        _no_store(detail)
        page = client.get(_collection(world), params={"limit": 1, "offset": 1})
        assert page.status_code == 200
        assert page.json()["items"][0]["archive_version_id"] == _V2
        assert page.json()["next_offset"] == 2
        invalid = client.get(_collection(world), params={"limit": 101})
        assert invalid.status_code == 400
        assert invalid.json()["issues"][0]["code"] == "archive_invalid"
        _no_store(invalid)
        assert "101" not in invalid.text

        v1 = client.get(_content(world, _V1))
        v2 = client.get(_content(world, _V2))
        assert v1.status_code == 200
        assert v2.status_code == 200
        assert v1.content == b"0123456789abcdef"
        assert v2.content == b"fedcba9876543210"
        assert v1.headers["accept-ranges"] == "bytes"
        assert v1.headers["content-length"] == "16"
        _no_store(v1)
        window = client.get(_content(world, _V1), headers={"Range": "bytes=0-3"})
        assert window.status_code == 206
        assert window.content == b"0123"
        assert window.headers["content-range"] == "bytes 0-3/16"
        assert world.objects.readers[-1].reads == [(0, 4)]
        suffix = client.get(_content(world, _V2), headers={"Range": "bytes=-4"})
        assert suffix.status_code == 206
        assert suffix.content == b"3210"
        assert suffix.headers["content-range"] == "bytes 12-15/16"
        missing = client.get(
            _content(world, _V1),
            headers={"Range": "bytes=20-30"},
        )
        assert missing.status_code == 416
        assert missing.content == b""
        assert missing.headers["content-range"] == "bytes */16"
        _no_store(missing)
        multi = client.get(
            _content(world, _V1),
            headers={"Range": "bytes=0-1,2-3"},
        )
        assert multi.status_code == 416
        head = client.request("HEAD", _content(world, _V1))
        assert head.status_code == 200
        assert head.content == b""
        assert head.headers["content-length"] == "16"
        assert world.objects.readers[-1].reads == []
        ranged_head = client.request(
            "HEAD",
            _content(world, _V1),
            headers={"Range": "bytes=0-3"},
        )
        assert ranged_head.status_code == 206
        assert ranged_head.content == b""
        assert ranged_head.headers["content-length"] == "4"

        world.v1.path.write_bytes(b"XXXXXXXXXXXXXXXX")
        replaced = client.get(_content(world, _V1))
        assert replaced.status_code == 409
        assert replaced.json()["issues"][0]["code"] == "archive_integrity"
        assert replaced.content != b"XXXXXXXXXXXXXXXX"
        assert b"XXXXXXXXXXXXXXXX" not in replaced.content
        _no_store(replaced)
        untouched = client.get(_content(world, _V2))
        assert untouched.status_code == 200
        assert untouched.content == b"fedcba9876543210"


def test_native_and_registered_reads_share_session_scope(seeded_app, tmp_path) -> None:
    world = _install(seeded_app, tmp_path)
    member, tokens = _member_session(seeded_app)
    with ApiTestClient(seeded_app) as client:
        seeded_app.state.test_auth_tokens = tokens
        run_bytes = client.get(_content(world, _V1))
        registered = client.get(_content(world, _V_REGISTERED))
        assert run_bytes.status_code == 200
        assert run_bytes.content == b"0123456789abcdef"
        assert registered.status_code == 200
        assert registered.content == b"reg!"
        _redacted(registered, _SECRET_PATH, _SECRET_SHA)
        unbound = client.get(_content(world, _V_UNBOUND))
        foreign = client.get(_content(world, _V_FOREIGN))
        assert unbound.status_code == 404
        assert foreign.status_code == 404
        assert unbound.json()["issues"][0]["code"] == "archive_version_not_found"
        assert foreign.json()["issues"][0]["code"] == "archive_version_not_found"
        assert world.v2 not in world.objects.attempts
        wrong_experiment = client.get(
            "/api/v1/catalog/projects/"
            f"{world.project.project_id}/experiments/{_OTHER_EXPERIMENT}"
            f"/archived-tracks/{_V1}/content"
        )
        assert wrong_experiment.status_code == 404
        _no_store(wrong_experiment)
        seeded_app.state.test_auth_tokens = ("a" * 43, "b" * 43)
        admin = client.get(_content(world, _V_REGISTERED))
        assert admin.status_code == 200

        resolve_session = seeded_app.state.authentication_service.resolve_session
        principal = AuthenticatedPrincipal(
            user_id=member.user_id,
            username="archive-member",
            role=UserRole.MEMBER,
        )
        admin_principal = AuthenticatedPrincipal.from_account(
            seeded_app.state.authentication_repository.get_account_by_id(
                TEST_AUTH_USER_ID
            )
        )
        seeded_app.state.account_administration_service.set_account_status(
            AuthenticationActor.for_principal(admin_principal),
            member.user_id,
            enabled=False,
        )
        seeded_app.state.authentication_service.resolve_session = (  # type: ignore[method-assign]
            lambda _token: principal
        )
        seeded_app.state.test_auth_tokens = tokens
        attempts = len(world.objects.attempts)
        disabled = client.get(_content(world, _V1))
        assert disabled.status_code == 401
        assert disabled.json()["issues"][0]["code"] == "AUTHENTICATION_REQUIRED"
        _no_store(disabled)
        assert len(world.objects.attempts) == attempts

        seeded_app.state.authentication_service.resolve_session = resolve_session
        seeded_app.state.test_auth_tokens = ("a" * 43, "b" * 43)
        _registry(seeded_app).archive_project(world.project.project_id)
        archived = client.get(_collection(world))
        assert archived.status_code == 404
        assert archived.json()["issues"][0]["code"] == "archive_version_not_found"


def test_index_and_data_stay_on_the_requested_version(seeded_app, tmp_path) -> None:
    world = _install(seeded_app, tmp_path)
    with ApiTestClient(seeded_app) as client:
        data = client.get(_content(world, _V_PAIR))
        index = client.get(_content(world, _V_PAIR, index=True))
        assert data.status_code == 200
        assert data.content == b"data-win"
        assert index.status_code == 200
        assert index.content == b"idx-window"
        assert data.content != index.content
        detail = client.get(f"{_collection(world)}/{_V_PAIR}")
        assert detail.json()["has_index"] is True
        assert detail.json()["track_type"] == "longrange"
        assert detail.json()["source"]["index_artifact_id"] == "signalIndex"
        bigwig_index = client.get(_content(world, _V1, index=True))
        assert bigwig_index.status_code == 404
        assert bigwig_index.content != b"0123456789abcdef"
        assert b"0123456789abcdef" not in bigwig_index.content
        _no_store(bigwig_index)
        world.index.path.unlink()
        missing_index = client.get(_content(world, _V_PAIR, index=True))
        assert missing_index.status_code == 409
        assert missing_index.json()["issues"][0]["code"] == "archive_offline"
        assert b"data-win" not in missing_index.content
        data_still = client.get(_content(world, _V_PAIR))
        assert data_still.status_code == 200
        assert data_still.content == b"data-win"
        other = client.get(_content(world, _V2))
        assert other.status_code == 200
        assert other.content == b"fedcba9876543210"
