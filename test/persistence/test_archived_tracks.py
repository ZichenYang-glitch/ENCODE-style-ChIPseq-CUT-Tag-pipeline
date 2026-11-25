"""Result versions survive concurrent seals and experiment metadata revisions."""

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone

import pytest

from encode_pipeline.persistence.archived_tracks import ArchiveRepository
from encode_pipeline.persistence.database import (
    create_database_engine,
    create_session_factory,
)
from encode_pipeline.persistence.migrations import downgrade_database, upgrade_database
from encode_pipeline.persistence.models import (
    ExperimentRevisionRow,
    ExperimentRow,
    ProjectRow,
    RegisteredResultRevisionRow,
    RegisteredResultRow,
    StoragePoolRow,
)
from encode_pipeline.platform.archived_tracks import ArchiveStoreError
from encode_pipeline.platform.experiment_catalog import (
    EXPERIMENT_REVISION_PAYLOAD_DIGEST_SCHEME,
)
from encode_pipeline.platform.registered_results import REVISION_DIGEST_SCHEME
from encode_pipeline.services.archive_objects import ArchiveObjects


def test_concurrent_publication_historical_selection_and_rollback(tmp_path):
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
    engine.dispose()
