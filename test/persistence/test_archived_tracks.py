"""Result versions survive concurrent seals and experiment metadata revisions."""

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

import pytest
from sqlalchemy import text
from sqlalchemy.exc import IntegrityError

from encode_pipeline.persistence.archived_tracks import ArchiveRepository
from encode_pipeline.persistence.database import (
    create_database_engine,
    create_session_factory,
)
from encode_pipeline.persistence.migrations import downgrade_database, upgrade_database
from encode_pipeline.persistence.models import (
    ArtifactPublicationRow,
    ExperimentRevisionRow,
    ExperimentRow,
    ProjectRow,
    RegisteredResultRevisionRow,
    RegisteredResultRow,
    RunExperimentBindingRow,
    RunProjectBindingRow,
    RunRow,
    StoragePoolRow,
)
from encode_pipeline.platform.archived_tracks import ArchiveStoreError
from encode_pipeline.platform.experiment_catalog import (
    EXPERIMENT_REVISION_PAYLOAD_DIGEST_SCHEME,
)
from encode_pipeline.platform.registered_results import REVISION_DIGEST_SCHEME
from encode_pipeline.services.archive_objects import ArchiveObjects


@pytest.fixture
def archive_case(tmp_path):
    url = f"sqlite:///{tmp_path / 'experiment.db'}"
    upgrade_database(url)
    engine = create_database_engine(url)
    sessions = create_session_factory(engine)
    repository = ArchiveRepository(sessions)
    now = datetime.now(timezone.utc)
    project, experiment = "prj_" + "a" * 32, "exp_" + "b" * 32
    result, pool = "rreg_" + "c" * 32, "stgp_" + "d" * 32
    archive = tmp_path / "archive"
    archive.mkdir()
    objects = ArchiveObjects(archive)
    versions = []
    with sessions.begin() as db:
        db.add(
            ProjectRow(
                project_id=project, display_name="Project", kind="user", created_at=now
            )
        )
        db.add(
            StoragePoolRow(
                storage_pool_id=pool,
                display_name="Pool",
                config_key="pool",
                created_at=now,
            )
        )
        db.flush()
        db.add(
            ExperimentRow(
                experiment_id=experiment,
                project_id=project,
                stable_key="exp",
                created_at=now,
            )
        )
        db.add(
            RegisteredResultRow(
                result_id=result,
                project_id=project,
                stable_key="signal",
                created_at=now,
            )
        )
        db.flush()
        for i in (1, 2):
            exp_rev = "expr_" + str(i) * 32
            result_rev = "rrev_" + str(i) * 32
            f = tmp_path / "signal.bw"
            f.write_bytes(bytes([i]) * 32)
            with f.open("rb") as handle:
                member = objects.seal_fd(handle.fileno())
            db.add(
                ExperimentRevisionRow(
                    experiment_revision_id=exp_rev,
                    experiment_id=experiment,
                    project_id=project,
                    revision_number=i,
                    canonical_payload="{}",
                    payload_digest_scheme=EXPERIMENT_REVISION_PAYLOAD_DIGEST_SCHEME,
                    payload_digest=str(i) * 64,
                    created_at=now,
                )
            )
            db.add(
                RegisteredResultRevisionRow(
                    result_revision_id=result_rev,
                    result_id=result,
                    revision_number=i,
                    canonical_payload="{}",
                    digest_scheme=REVISION_DIGEST_SCHEME,
                    revision_digest=str(i) * 64,
                    created_at=now,
                    storage_pool_id=pool,
                    relative_path=f"v{i}/signal.bw",
                    size_bytes=member.size_bytes,
                    content_sha256=member.sha256,
                )
            )
            versions.append(
                dict(
                    project_id=project,
                    experiment_id=experiment,
                    experiment_revision_id=exp_rev,
                    track_key="signal",
                    display_name="signal.bw",
                    assembly="hg38",
                    track_type="bigwig",
                    source=dict(
                        kind="registered",
                        result_id=result,
                        result_revision_id=result_rev,
                    ),
                    members={"data": member},
                )
            )
    yield url, sessions, repository, project, experiment, objects, archive, versions
    engine.dispose()


def test_concurrent_publication_historical_selection_and_rollback(
    tmp_path, archive_case
):
    url, sessions, repository, project, experiment, objects, archive, versions = (
        archive_case
    )
    with ThreadPoolExecutor(3) as runner:
        first = list(runner.map(lambda _: repository.publish(**versions[0]), range(3)))
    assert len({v.archive_version_id for v in first}) == 1
    assert first[0].version_number == 1
    second = repository.publish(**versions[1])
    assert second.version_number == 2
    heads, offset = repository.list_versions(project, experiment, limit=1)
    assert [v.archive_version_id for v in heads] == [second.archive_version_id]
    assert offset == 1
    history, offset = repository.list_versions(
        project, experiment, limit=1, offset=offset
    )
    assert [v.archive_version_id for v in history] == [first[0].archive_version_id]
    assert offset is None
    exact = repository.get_version(project, experiment, first[0].archive_version_id)
    assert exact.experiment_revision_id == versions[0]["experiment_revision_id"]
    with objects.open_member(exact.members["data"]) as reader:
        assert reader.read(0, 1) == bytes([1])
    with pytest.raises(ArchiveStoreError, match="archive_version_conflict"):
        repository.publish(**(versions[0] | {"display_name": "different"}))
    with pytest.raises(ArchiveStoreError, match="archive_integrity"):
        repository.publish(**(versions[0] | {"members": versions[1]["members"]}))
    with pytest.raises(ArchiveStoreError, match="archive_version_not_found"):
        repository.get_version(project, "exp_" + "9" * 32, first[0].archive_version_id)
    # A disk restore is a new filesystem identity, even for the same bytes.
    # An operator re-verifies the complete restored member before republishing.
    restored = tmp_path / "restored"
    restored.write_bytes(bytes([1]) * 32)
    restored.chmod(0o440)
    restored.replace(archive / exact.members["data"].sha256)
    with pytest.raises(ArchiveStoreError, match="archive_integrity"):
        objects.open_member(exact.members["data"])
    with (archive / exact.members["data"].sha256).open("rb") as f:
        proof = objects.seal_fd(f.fileno(), exact.members["data"].sha256)
    rebound = repository.publish(**(versions[0] | {"members": {"data": proof}}))
    assert rebound.archive_version_id == exact.archive_version_id
    with objects.open_member(
        repository.get_version(project, experiment, exact.archive_version_id).members[
            "data"
        ]
    ) as reader:
        assert reader.read(0, 1) == bytes([1])
    with pytest.raises(RuntimeError, match="sealed track versions"):
        downgrade_database(url, "20261119_39")
    assert (
        repository.get_version(
            project, experiment, second.archive_version_id
        ).version_number
        == 2
    )
    with sessions.begin() as db:
        db.get(ProjectRow, project).archived_at = datetime.now(timezone.utc)
    with pytest.raises(ArchiveStoreError, match="archive_version_not_found"):
        repository.get_version(project, experiment, exact.archive_version_id)


def test_published_identity_and_bytes_are_immutable_even_in_sql(archive_case):
    _, sessions, repository, _, _, _, _, versions = archive_case
    repository.publish(**versions[0])
    for statement in (
        "UPDATE archived_track_versions SET archive_version_id='other'",
        "UPDATE archived_track_versions SET members=json_set(members, '$.data.sha256', 'wrong')",
        "DELETE FROM archived_track_versions",
    ):
        with pytest.raises(IntegrityError, match="immutable"), sessions.begin() as db:
            db.execute(text(statement))
    with pytest.raises(IntegrityError, match="scope mismatch"), sessions.begin() as db:
        db.execute(
            text("""
          INSERT INTO archived_track_versions
          SELECT 'another-id', project_id, experiment_id, 'wrong-revision', track_key,
                 3, display_name, assembly, track_type, source_identity, source, members
          FROM archived_track_versions
        """)
        )


def test_publish_requires_distinct_members_and_a_consistent_run_binding(archive_case):
    _, sessions, repo, project, _, _, _, versions = archive_case
    data = versions[0]
    src = data["source"]
    with pytest.raises(ArchiveStoreError, match="archive_invalid"):
        repo.publish(
            **(
                data
                | {
                    "track_type": "longrange",
                    "members": {
                        "data": data["members"]["data"],
                        "index": data["members"]["data"],
                    },
                    "source": {
                        **src,
                        "index_result_id": src["result_id"],
                        "index_result_revision_id": src["result_revision_id"],
                    },
                }
            )
        )
    now = datetime.now(timezone.utc)
    run_id = "changing-binding"
    with sessions.begin() as db:
        db.add(
            RunRow(
                run_id=run_id,
                workflow_id="test",
                tags={},
                inputs={},
                status="succeeded",
                created_at=now,
                updated_at=now,
            )
        )
        db.flush()
        db.add(
            RunProjectBindingRow(
                run_id=run_id,
                project_id=project,
                binding_mode="bound_v1",
                provenance="resolved",
                workflow_inputs_digest="0" * 64,
                binding_digest_scheme="sha256-framed-project-sample-binding-v1",
                binding_digest="0" * 64,
                created_at=now,
            )
        )
        db.add(
            RunExperimentBindingRow(
                run_id=run_id,
                experiment_revision_id=data["experiment_revision_id"],
                payload_digest="9" * 64,
            )
        )
    with sessions.begin() as db:
        db.add(
            ArtifactPublicationRow(
                run_id=run_id,
                artifact_id="signal",
                artifact_generation="artifactgen-" + "a" * 64,
                artifact_revision="artifactrev-" + "b" * 64,
                output_type="bigwig",
                published_at=now,
            )
        )
    with pytest.raises(ArchiveStoreError, match="archive_version_not_found"):
        repo.publish(
            **(
                data
                | {
                    "source": dict(
                        kind="run",
                        run_id=run_id,
                        artifact_id="signal",
                        generation="artifactgen-" + "a" * 64,
                        revision="artifactrev-" + "b" * 64,
                    )
                }
            )
        )


def test_invalid_operator_inputs_are_reported_without_serialization_or_type_errors(
    archive_case,
):
    _, _, repo, project, experiment, _, _, versions = archive_case
    for overrides in (
        {"members": None},
        {"source": {1: "a", "kind": "run"}},
        {"track_key": None},
    ):
        with pytest.raises(ArchiveStoreError, match="archive_invalid"):
            repo.publish(**(versions[0] | overrides))
    with pytest.raises(ArchiveStoreError, match="archive_invalid"):
        repo.list_versions(project, experiment, limit=True)
