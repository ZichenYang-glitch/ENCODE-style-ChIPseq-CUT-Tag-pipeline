"""SQLAlchemy Project/Sample registry repository contract tests."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone

import pytest

from encode_pipeline.persistence import (
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.persistence.data_registry import (
    SqlAlchemyDataRegistryRepository,
)
from encode_pipeline.platform.data_registry import (
    LEGACY_PROJECT_ID,
    SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
    Project,
    ProjectKind,
    ProjectSampleSelection,
    Sample,
    SampleRevision,
    build_sample_revision_payload_digest,
    canonical_sample_revision_payload,
)
from encode_pipeline.services.data_registry_repositories import (
    DataRegistryConflictError,
    InMemoryDataRegistryRepository,
)


NOW = datetime(2026, 7, 26, 8, 0, tzinfo=timezone.utc)
PROJECT_ID = "prj_11111111111111111111111111111111"
SAMPLE_ID = "smp_11111111111111111111111111111111"


@pytest.fixture
def repository(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'registry.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    yield SqlAlchemyDataRegistryRepository(create_session_factory(engine))
    engine.dispose()


def _project(
    project_id: str = PROJECT_ID,
    *,
    display_name: str = "Project One",
    created_at: datetime = NOW,
) -> Project:
    return Project(
        project_id=project_id,
        display_name=display_name,
        kind=ProjectKind.USER,
        created_at=created_at,
    )


def _sample(
    sample_id: str = SAMPLE_ID,
    *,
    project_id: str = PROJECT_ID,
    stable_key: str = "tumor-a",
    created_at: datetime = NOW,
) -> Sample:
    return Sample(
        sample_id=sample_id,
        project_id=project_id,
        stable_key=stable_key,
        created_at=created_at,
    )


def _revision(
    revision_number: int = 1,
    *,
    sample_revision_id: str | None = None,
    sample_id: str = SAMPLE_ID,
    project_id: str = PROJECT_ID,
    display_name: str = "Tumor A",
    attributes: dict[str, object] | None = None,
    created_at: datetime | None = None,
) -> SampleRevision:
    if sample_revision_id is None:
        sample_revision_id = f"smpr_{revision_number:032x}"
    if attributes is None:
        attributes = {"condition": "treated", "replicate": revision_number}
    if created_at is None:
        created_at = NOW + timedelta(minutes=revision_number)
    canonical_payload = canonical_sample_revision_payload(
        display_name=display_name,
        attributes=attributes,
    )
    return SampleRevision(
        sample_revision_id=sample_revision_id,
        project_id=project_id,
        sample_id=sample_id,
        revision_number=revision_number,
        canonical_payload=canonical_payload,
        payload_digest_scheme=SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
        payload_digest=build_sample_revision_payload_digest(canonical_payload),
        created_at=created_at,
    )


def test_repository_creates_lists_and_archives_projects(repository) -> None:
    project = _project()

    assert repository.create_project(project) == project
    assert repository.get_project(project.project_id) == project
    assert repository.list_projects() == (
        repository.get_project(LEGACY_PROJECT_ID),
        project,
    )

    archived = repository.archive_project(
        project.project_id,
        archived_at=NOW + timedelta(hours=1),
    )
    assert archived.archived_at == NOW + timedelta(hours=1)
    assert repository.list_projects() == (repository.get_project(LEGACY_PROJECT_ID),)
    assert repository.list_projects(include_archived=True) == (
        repository.get_project(LEGACY_PROJECT_ID),
        archived,
    )
    with pytest.raises(DataRegistryConflictError, match="already archived"):
        repository.archive_project(
            project.project_id,
            archived_at=NOW + timedelta(hours=2),
        )
    with pytest.raises(ValueError, match="Legacy Project"):
        repository.archive_project(
            LEGACY_PROJECT_ID,
            archived_at=NOW + timedelta(hours=2),
        )


def test_repository_creates_sample_and_appends_immutable_revisions(repository) -> None:
    project = repository.create_project(_project())
    sample = _sample(project_id=project.project_id)
    initial = _revision()

    assert repository.create_sample(sample, initial) == (sample, initial)
    assert repository.get_sample(sample.sample_id) == sample
    assert repository.list_samples(project.project_id) == (sample,)
    assert repository.get_sample_revision(initial.sample_revision_id) == initial
    assert repository.list_sample_revisions(sample.sample_id) == (initial,)

    revised = _revision(2, display_name="Tumor A revised")
    assert (
        repository.append_sample_revision(
            revised,
            expected_previous_revision_number=1,
        )
        == revised
    )
    assert repository.list_sample_revisions(sample.sample_id) == (initial, revised)
    assert repository.get_sample_revision(initial.sample_revision_id) == initial
    assert initial.display_name == "Tumor A"
    assert revised.display_name == "Tumor A revised"


def test_repository_batch_sample_create_is_atomic_on_late_conflict(repository) -> None:
    repository.create_project(_project())
    first_sample = _sample()
    first_revision = _revision()
    second_sample = _sample(
        "smp_22222222222222222222222222222222",
        stable_key=first_sample.stable_key,
    )
    second_revision = _revision(
        sample_revision_id="smpr_22222222222222222222222222222222",
        sample_id=second_sample.sample_id,
    )

    with pytest.raises(DataRegistryConflictError, match="stable_key"):
        repository.create_samples(
            (
                (first_sample, first_revision),
                (second_sample, second_revision),
            )
        )

    assert repository.list_samples(PROJECT_ID) == ()
    with pytest.raises(KeyError):
        repository.get_sample_revision(first_revision.sample_revision_id)


def test_repository_list_order_matches_memory_for_same_time_entities(
    repository,
) -> None:
    memory = InMemoryDataRegistryRepository(
        legacy_created_at=datetime(1970, 1, 1, tzinfo=timezone.utc)
    )
    project = _project()
    for candidate in (memory, repository):
        for created_project in (
            _project(
                "prj_" + "f" * 32,
                display_name="Project F",
            ),
            _project(
                "prj_" + "2" * 32,
                display_name="Project Two",
            ),
            project,
        ):
            candidate.create_project(created_project)
        entries = tuple(
            (
                _sample(
                    sample_id,
                    stable_key=f"sample-{suffix}",
                    created_at=NOW,
                ),
                _revision(
                    sample_revision_id=f"smpr_{suffix * 32}",
                    sample_id=sample_id,
                    created_at=NOW,
                ),
            )
            for sample_id, suffix in (
                ("smp_" + "f" * 32, "f"),
                ("smp_" + "1" * 32, "1"),
            )
        )
        candidate.create_samples(entries)

    expected_project_ids = (
        LEGACY_PROJECT_ID,
        PROJECT_ID,
        "prj_" + "2" * 32,
        "prj_" + "f" * 32,
    )
    assert (
        tuple(project.project_id for project in memory.list_projects())
        == expected_project_ids
    )
    assert (
        tuple(project.project_id for project in repository.list_projects())
        == expected_project_ids
    )
    assert (
        memory.get_project(LEGACY_PROJECT_ID).created_at
        == repository.get_project(LEGACY_PROJECT_ID).created_at
    )

    expected_ids = ("smp_" + "1" * 32, "smp_" + "f" * 32)
    assert tuple(sample.sample_id for sample in memory.list_samples(PROJECT_ID)) == (
        expected_ids
    )
    assert (
        tuple(sample.sample_id for sample in repository.list_samples(PROJECT_ID))
        == expected_ids
    )


def test_repository_resolves_ordered_same_project_digest_pinned_selection(
    repository,
) -> None:
    repository.create_project(_project())
    first = (_sample(), _revision())
    second = (
        _sample(
            "smp_22222222222222222222222222222222",
            stable_key="control-a",
        ),
        _revision(
            sample_revision_id="smpr_22222222222222222222222222222222",
            sample_id="smp_22222222222222222222222222222222",
            display_name="Control A",
        ),
    )
    repository.create_samples((first, second))

    binding = repository.resolve_project_sample_selection(
        ProjectSampleSelection(
            project_id=PROJECT_ID,
            sample_revision_ids=(
                second[1].sample_revision_id,
                first[1].sample_revision_id,
            ),
        ),
        workflow_inputs_digest="a" * 64,
    )

    assert binding.sample_revision_ids == (
        second[1].sample_revision_id,
        first[1].sample_revision_id,
    )
    assert tuple(ref.payload_digest for ref in binding.sample_revisions) == (
        second[1].payload_digest,
        first[1].payload_digest,
    )
    assert binding.workflow_inputs_digest == "a" * 64


def test_repository_revision_compare_and_append_uses_sqlite_write_lock(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'revision-race.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    repository = SqlAlchemyDataRegistryRepository(create_session_factory(engine))
    repository.create_project(_project())
    repository.create_sample(_sample(), _revision())
    revisions = (
        _revision(
            2,
            sample_revision_id="smpr_22222222222222222222222222222222",
            display_name="Concurrent A",
        ),
        _revision(
            2,
            sample_revision_id="smpr_33333333333333333333333333333333",
            display_name="Concurrent B",
        ),
    )

    def append(revision: SampleRevision) -> str:
        try:
            repository.append_sample_revision(
                revision,
                expected_previous_revision_number=1,
            )
        except DataRegistryConflictError:
            return "conflict"
        return "created"

    with ThreadPoolExecutor(max_workers=2) as pool:
        outcomes = tuple(pool.map(append, revisions))

    assert sorted(outcomes) == ["conflict", "created"]
    assert [
        revision.revision_number
        for revision in repository.list_sample_revisions(SAMPLE_ID)
    ] == [1, 2]
    engine.dispose()


def test_repository_round_trips_all_timestamps_as_utc(repository) -> None:
    local_tz = timezone(timedelta(hours=8))
    local_created = datetime(2026, 7, 26, 16, 0, tzinfo=local_tz)
    project = repository.create_project(_project(created_at=local_created))
    sample = _sample(created_at=local_created)
    revision = _revision(created_at=local_created)
    repository.create_sample(sample, revision)

    assert project.created_at == NOW
    assert repository.get_project(PROJECT_ID).created_at == NOW
    assert repository.get_sample(SAMPLE_ID).created_at == NOW
    assert repository.get_sample_revision(revision.sample_revision_id).created_at == NOW
    assert repository.get_project(PROJECT_ID).created_at.tzinfo is timezone.utc


def test_repository_rejects_archived_and_cross_project_sample_bindings(
    repository,
) -> None:
    repository.create_project(_project())
    repository.create_sample(_sample(), _revision())
    other_project_id = "prj_22222222222222222222222222222222"
    repository.create_project(_project(other_project_id, display_name="Project Two"))

    with pytest.raises(ValueError, match="same Project"):
        repository.resolve_project_sample_selection(
            ProjectSampleSelection(
                project_id=other_project_id,
                sample_revision_ids=("smpr_00000000000000000000000000000001",),
            ),
            workflow_inputs_digest="a" * 64,
        )

    repository.archive_project(
        PROJECT_ID,
        archived_at=NOW + timedelta(hours=1),
    )
    with pytest.raises(ValueError, match="archived"):
        repository.resolve_project_sample_selection(
            ProjectSampleSelection(
                project_id=PROJECT_ID,
                sample_revision_ids=("smpr_00000000000000000000000000000001",),
            ),
            workflow_inputs_digest="a" * 64,
        )
    with pytest.raises(ValueError, match="archived"):
        repository.create_sample(
            _sample("smp_33333333333333333333333333333333", stable_key="late"),
            _revision(
                sample_revision_id="smpr_33333333333333333333333333333333",
                sample_id="smp_33333333333333333333333333333333",
            ),
        )
