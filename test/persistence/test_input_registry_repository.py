"""SQLAlchemy and in-memory Input registry repository parity tests."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from itertools import count
import stat

import pytest

from encode_pipeline.persistence import (
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.persistence.data_registry import (
    SqlAlchemyDataRegistryRepository,
)
from encode_pipeline.persistence.input_registry import (
    SqlAlchemyInputRegistryRepository,
)
from encode_pipeline.platform.data_registry import Project, ProjectKind
from encode_pipeline.platform.input_registry import (
    InputFile,
    InputProvenanceMode,
    InputUseBindingPlan,
    PlannedInputUse,
    ProjectStoragePoolBinding,
    StoragePool,
    build_input_file_revision,
)
from encode_pipeline.services.data_registry_repositories import (
    InMemoryDataRegistryRepository,
)
from encode_pipeline.services.input_file_access import FileObservation
from encode_pipeline.services.input_registry import InputRegistryService
from encode_pipeline.services.input_registry_repositories import (
    InMemoryInputRegistryRepository,
    InputRegistryConflictError,
)


NOW = datetime(2026, 7, 26, 9, 0, tzinfo=timezone.utc)
PROJECT_ID = "prj_11111111111111111111111111111111"
OTHER_PROJECT_ID = "prj_22222222222222222222222222222222"


def _ids(prefix: str):
    serial = count(1)
    return lambda: f"{prefix}_{next(serial):032x}"


def _project(
    project_id: str = PROJECT_ID,
    display_name: str = "Project One",
) -> Project:
    return Project(
        project_id=project_id,
        display_name=display_name,
        kind=ProjectKind.USER,
        created_at=NOW,
    )


def _observation(
    *,
    relative_path: str,
    size_bytes: int,
    content_sha256: str,
) -> FileObservation:
    return FileObservation(
        relative_path=relative_path,
        size_bytes=size_bytes,
        content_sha256=content_sha256,
        path_fingerprint=(
            (1, 2),
            (
                1,
                3,
                stat.S_IFREG | 0o640,
                1,
                size_bytes,
                1_000_000,
                1_000_000,
            ),
        ),
    )


def _service(repository) -> InputRegistryService:
    return InputRegistryService(
        repository=repository,
        storage_pool_id_factory=_ids("stgp"),
        input_file_id_factory=_ids("inpf"),
        input_file_revision_id_factory=_ids("inpfr"),
        now_factory=lambda: NOW,
    )


def _repository_candidates(sql_repositories):
    _sql_projects, sql_inputs = sql_repositories
    memory_projects = InMemoryDataRegistryRepository(legacy_created_at=NOW)
    memory_projects.create_project(_project())
    return (
        sql_inputs,
        InMemoryInputRegistryRepository(project_repository=memory_projects),
    )


def _input_use_plan(
    *,
    revision_id: str | None,
    include_transitional: bool,
    project_id: str = PROJECT_ID,
    closure_contract_version: str = "regular_file_v1",
) -> InputUseBindingPlan:
    input_uses = []
    if revision_id is not None:
        input_uses.append(
            PlannedInputUse(
                key="reads",
                occurrence=0,
                capability_version="reads-v1",
                closure_contract_version=closure_contract_version,
                provenance_mode=InputProvenanceMode.MANAGED_REVISION_V1,
                input_file_revision_ids=(revision_id,),
            )
        )
    if include_transitional:
        input_uses.append(
            PlannedInputUse(
                key="reference",
                occurrence=0,
                capability_version="reference-v1",
                closure_contract_version="index-directory-v1",
                provenance_mode=InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1,
                input_file_revision_ids=(),
            )
        )
    return InputUseBindingPlan(
        project_id=project_id,
        workflow_id="test-workflow",
        adapter_contract_version="adapter-v2",
        input_uses=tuple(input_uses),
    )


@pytest.fixture
def sql_repositories(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'input-registry.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    session_factory = create_session_factory(engine)
    projects = SqlAlchemyDataRegistryRepository(session_factory)
    inputs = SqlAlchemyInputRegistryRepository(session_factory)
    projects.create_project(_project())
    yield projects, inputs
    engine.dispose()


def test_sql_and_memory_repositories_have_registry_contract_parity(
    sql_repositories,
) -> None:
    sql_projects, sql_inputs = sql_repositories
    memory_projects = InMemoryDataRegistryRepository(legacy_created_at=NOW)
    memory_projects.create_project(_project())
    candidates = (
        _service(sql_inputs),
        _service(InMemoryInputRegistryRepository(project_repository=memory_projects)),
    )

    observations = (
        _observation(
            relative_path="reads/a.fastq.gz",
            size_bytes=100,
            content_sha256="a" * 64,
        ),
        _observation(
            relative_path="reads/a-v2.fastq.gz",
            size_bytes=101,
            content_sha256="b" * 64,
        ),
    )
    results = []
    for service in candidates:
        pool = service.create_storage_pool(
            display_name="Ingress",
            config_key="ingress",
        )
        binding = service.bind_project_storage_pool(
            PROJECT_ID,
            pool.storage_pool_id,
        )
        initial = service.register_input_file(
            PROJECT_ID,
            stable_key="reads-a",
            observation=observations[0],
        )
        repeated = service.register_input_file(
            PROJECT_ID,
            stable_key="reads-a",
            observation=observations[0],
        )
        revised = service.register_input_file(
            PROJECT_ID,
            stable_key="reads-a",
            observation=observations[1],
        )
        results.append(
            (
                pool,
                binding,
                initial,
                repeated,
                revised,
                service.list_input_files(PROJECT_ID),
                service.list_input_file_revisions(initial.input_file.input_file_id),
            )
        )

    assert results[0] == results[1]
    assert results[0][2] == results[0][3]
    assert [revision.revision_number for revision in results[0][-1]] == [1, 2]
    assert sql_projects.get_project(PROJECT_ID) == _project()


def test_sql_repository_serializes_compare_and_append_registration(tmp_path) -> None:
    database_url = f"sqlite:///{tmp_path / 'input-revision-race.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    session_factory = create_session_factory(engine)
    project_repository = SqlAlchemyDataRegistryRepository(session_factory)
    project_repository.create_project(_project())
    bootstrap = _service(SqlAlchemyInputRegistryRepository(session_factory))
    pool = bootstrap.create_storage_pool(
        display_name="Ingress",
        config_key="ingress",
    )
    bootstrap.bind_project_storage_pool(PROJECT_ID, pool.storage_pool_id)
    initial = bootstrap.register_input_file(
        PROJECT_ID,
        stable_key="reads-a",
        observation=_observation(
            relative_path="reads/a.fastq.gz",
            size_bytes=100,
            content_sha256="a" * 64,
        ),
    )

    services = tuple(
        InputRegistryService(
            repository=SqlAlchemyInputRegistryRepository(session_factory),
            storage_pool_id_factory=lambda: "stgp_" + "f" * 32,
            input_file_id_factory=lambda: "inpf_" + "f" * 32,
            input_file_revision_id_factory=(
                lambda suffix=suffix: f"inpfr_{suffix * 32}"
            ),
            now_factory=lambda: NOW + timedelta(minutes=1),
        )
        for suffix in ("2", "3")
    )
    observation = _observation(
        relative_path="reads/a-v2.fastq.gz",
        size_bytes=101,
        content_sha256="b" * 64,
    )

    with ThreadPoolExecutor(max_workers=2) as pool_executor:
        registered = tuple(
            pool_executor.map(
                lambda service: service.register_input_file(
                    PROJECT_ID,
                    stable_key="reads-a",
                    observation=observation,
                ),
                services,
            )
        )

    assert registered[0] == registered[1]
    assert registered[0].input_file == initial.input_file
    revisions = bootstrap.list_input_file_revisions(initial.input_file.input_file_id)
    assert [revision.revision_number for revision in revisions] == [1, 2]
    engine.dispose()


def test_sql_repository_round_trips_registry_timestamps_as_utc(
    sql_repositories,
) -> None:
    _projects, repository = sql_repositories
    local_tz = timezone(timedelta(hours=8))
    local_now = datetime(2026, 7, 26, 17, 0, tzinfo=local_tz)
    service = InputRegistryService(
        repository=repository,
        storage_pool_id_factory=_ids("stgp"),
        input_file_id_factory=_ids("inpf"),
        input_file_revision_id_factory=_ids("inpfr"),
        now_factory=lambda: local_now,
    )

    pool = service.create_storage_pool(display_name="Ingress", config_key="ingress")
    binding = service.bind_project_storage_pool(PROJECT_ID, pool.storage_pool_id)
    registered = service.register_input_file(
        PROJECT_ID,
        stable_key="reads-a",
        observation=_observation(
            relative_path="reads/a.fastq.gz",
            size_bytes=100,
            content_sha256="a" * 64,
        ),
    )

    assert pool.created_at == NOW
    assert binding.bound_at == NOW
    assert registered.input_file.created_at == NOW
    assert registered.revision.created_at == NOW
    assert repository.get_storage_pool(pool.storage_pool_id).created_at.tzinfo is (
        timezone.utc
    )


@pytest.mark.parametrize("pool_state", ["missing", "archived"])
def test_sql_and_memory_pure_transitional_plan_does_not_require_active_pool(
    sql_repositories,
    pool_state,
) -> None:
    envelopes = []

    for repository in _repository_candidates(sql_repositories):
        service = _service(repository)
        if pool_state == "archived":
            pool = service.create_storage_pool(
                display_name="Archived ingress",
                config_key="archived-ingress",
            )
            service.bind_project_storage_pool(PROJECT_ID, pool.storage_pool_id)
            service.archive_storage_pool(pool.storage_pool_id)

        envelopes.append(
            repository.resolve_input_use_binding_plan(
                _input_use_plan(
                    revision_id=None,
                    include_transitional=True,
                ),
                project_sample_binding_digest="b" * 64,
                workflow_inputs_digest="c" * 64,
            )
        )

    assert envelopes[0] == envelopes[1]
    assert envelopes[0].contract_mode.value == "declared_input_uses_v1"
    assert envelopes[0].adapter_contract_version == "adapter-v2"
    assert envelopes[0].fully_managed is False
    assert envelopes[0].input_uses[0].provenance_mode is (
        InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1
    )
    assert envelopes[0].input_uses[0].members == ()
    assert envelopes[0].input_uses[0].closure_digest_scheme is None
    assert envelopes[0].input_uses[0].closure_digest is None


@pytest.mark.parametrize(
    "include_transitional",
    [False, True],
    ids=["pure-managed", "mixed"],
)
def test_sql_and_memory_resolve_identical_path_free_managed_evidence(
    sql_repositories,
    include_transitional,
) -> None:
    envelopes = []

    for repository in _repository_candidates(sql_repositories):
        service = _service(repository)
        pool = service.create_storage_pool(
            display_name="Ingress",
            config_key="ingress",
        )
        service.bind_project_storage_pool(PROJECT_ID, pool.storage_pool_id)
        registered = service.register_input_file(
            PROJECT_ID,
            stable_key="reads-a",
            observation=_observation(
                relative_path="private/reads/a.fastq.gz",
                size_bytes=100,
                content_sha256="a" * 64,
            ),
        )

        envelopes.append(
            repository.resolve_input_use_binding_plan(
                _input_use_plan(
                    revision_id=registered.revision.input_file_revision_id,
                    include_transitional=include_transitional,
                ),
                project_sample_binding_digest="b" * 64,
                workflow_inputs_digest="c" * 64,
            )
        )

    assert envelopes[0] == envelopes[1]
    assert envelopes[0].project_sample_binding_digest == "b" * 64
    managed = envelopes[0].input_uses[0]
    assert managed.members[0].logical_member_key == "file"
    assert managed.members[0].input_file_revision_id.startswith("inpfr_")
    assert not hasattr(managed.members[0], "relative_path")
    if include_transitional:
        transitional = envelopes[0].input_uses[1]
        assert transitional.members == ()
        assert transitional.closure_digest is None


@pytest.mark.parametrize(
    ("include_transitional", "pool_state"),
    [
        (False, "missing"),
        (False, "archived"),
        (True, "missing"),
        (True, "archived"),
    ],
    ids=[
        "pure-managed-missing",
        "pure-managed-inactive",
        "mixed-missing",
        "mixed-inactive",
    ],
)
def test_sql_and_memory_managed_uses_require_active_pool(
    sql_repositories,
    include_transitional,
    pool_state,
) -> None:
    for repository in _repository_candidates(sql_repositories):
        service = _service(repository)
        revision_id = "inpfr_" + "9" * 32
        if pool_state == "archived":
            pool = service.create_storage_pool(
                display_name="Archived ingress",
                config_key="archived-ingress",
            )
            service.bind_project_storage_pool(PROJECT_ID, pool.storage_pool_id)
            registered = service.register_input_file(
                PROJECT_ID,
                stable_key="reads-a",
                observation=_observation(
                    relative_path="private/reads/a.fastq.gz",
                    size_bytes=100,
                    content_sha256="a" * 64,
                ),
            )
            revision_id = registered.revision.input_file_revision_id
            service.archive_storage_pool(pool.storage_pool_id)

        error = KeyError if pool_state == "missing" else ValueError
        message = (
            "StoragePool binding"
            if pool_state == "missing"
            else "StoragePool is archived"
        )
        with pytest.raises(error, match=message):
            repository.resolve_input_use_binding_plan(
                _input_use_plan(
                    revision_id=revision_id,
                    include_transitional=include_transitional,
                ),
                project_sample_binding_digest="b" * 64,
                workflow_inputs_digest="c" * 64,
            )


def test_sql_and_memory_reject_cross_project_unknown_and_contract_mismatch(
    sql_repositories,
) -> None:
    sql_projects, sql_inputs = sql_repositories
    memory_projects = InMemoryDataRegistryRepository(legacy_created_at=NOW)
    memory_projects.create_project(_project())
    candidates = (
        (sql_projects, sql_inputs),
        (
            memory_projects,
            InMemoryInputRegistryRepository(project_repository=memory_projects),
        ),
    )

    for projects, repository in candidates:
        projects.create_project(_project(OTHER_PROJECT_ID, "Project Two"))
        service = _service(repository)
        owner_pool = service.create_storage_pool(
            display_name="Owner ingress",
            config_key="owner-ingress",
        )
        service.bind_project_storage_pool(PROJECT_ID, owner_pool.storage_pool_id)
        registered = service.register_input_file(
            PROJECT_ID,
            stable_key="reads-a",
            observation=_observation(
                relative_path="private/reads/a.fastq.gz",
                size_bytes=100,
                content_sha256="a" * 64,
            ),
        )
        other_pool = service.create_storage_pool(
            display_name="Other ingress",
            config_key="other-ingress",
        )
        service.bind_project_storage_pool(
            OTHER_PROJECT_ID,
            other_pool.storage_pool_id,
        )

        with pytest.raises(ValueError, match="not authorized"):
            repository.resolve_input_use_binding_plan(
                _input_use_plan(
                    project_id=OTHER_PROJECT_ID,
                    revision_id=registered.revision.input_file_revision_id,
                    include_transitional=False,
                ),
                project_sample_binding_digest="b" * 64,
                workflow_inputs_digest="c" * 64,
            )
        with pytest.raises(KeyError, match="input_file_revision_id"):
            repository.resolve_input_use_binding_plan(
                _input_use_plan(
                    revision_id="inpfr_" + "9" * 32,
                    include_transitional=False,
                ),
                project_sample_binding_digest="b" * 64,
                workflow_inputs_digest="c" * 64,
            )
        with pytest.raises(ValueError, match="closure contract"):
            repository.resolve_input_use_binding_plan(
                _input_use_plan(
                    revision_id=registered.revision.input_file_revision_id,
                    include_transitional=False,
                    closure_contract_version="directory_index_v1",
                ),
                project_sample_binding_digest="b" * 64,
                workflow_inputs_digest="c" * 64,
            )


def test_sql_create_input_file_rolls_back_after_initial_revision_conflict(
    sql_repositories,
) -> None:
    _projects, repository = sql_repositories
    service = _service(repository)
    pool = service.create_storage_pool(display_name="Ingress", config_key="ingress")
    service.bind_project_storage_pool(PROJECT_ID, pool.storage_pool_id)
    registered = service.register_input_file(
        PROJECT_ID,
        stable_key="reads-a",
        observation=_observation(
            relative_path="reads/a.fastq.gz",
            size_bytes=100,
            content_sha256="a" * 64,
        ),
    )
    new_input_file = InputFile(
        input_file_id="inpf_" + "8" * 32,
        project_id=PROJECT_ID,
        storage_pool_id=pool.storage_pool_id,
        stable_key="reads-b",
        created_at=NOW,
    )
    colliding_revision = build_input_file_revision(
        input_file_revision_id=registered.revision.input_file_revision_id,
        input_file_id=new_input_file.input_file_id,
        project_id=PROJECT_ID,
        storage_pool_id=pool.storage_pool_id,
        revision_number=1,
        relative_path="reads/b.fastq.gz",
        size_bytes=200,
        content_sha256="b" * 64,
        created_at=NOW,
    )

    with pytest.raises(InputRegistryConflictError, match="registry identity"):
        repository.create_input_file(new_input_file, colliding_revision)

    with pytest.raises(KeyError, match="input_file_id"):
        repository.get_input_file(new_input_file.input_file_id)
    assert (
        repository.get_input_file(registered.input_file.input_file_id)
        == registered.input_file
    )
    assert repository.list_input_file_revisions(
        registered.input_file.input_file_id
    ) == (registered.revision,)


def test_sql_append_contract_rejects_scope_order_time_and_identity_conflicts(
    sql_repositories,
) -> None:
    _projects, repository = sql_repositories
    service = _service(repository)
    pool = service.create_storage_pool(display_name="Ingress", config_key="ingress")
    service.bind_project_storage_pool(PROJECT_ID, pool.storage_pool_id)
    registered = service.register_input_file(
        PROJECT_ID,
        stable_key="reads-a",
        observation=_observation(
            relative_path="reads/a.fastq.gz",
            size_bytes=100,
            content_sha256="a" * 64,
        ),
    )

    def revision(
        *,
        revision_id: str,
        input_file_id: str = registered.input_file.input_file_id,
        project_id: str = PROJECT_ID,
        revision_number: int = 2,
        created_at: datetime = NOW + timedelta(minutes=1),
    ):
        return build_input_file_revision(
            input_file_revision_id=revision_id,
            input_file_id=input_file_id,
            project_id=project_id,
            storage_pool_id=pool.storage_pool_id,
            revision_number=revision_number,
            relative_path="reads/a-v2.fastq.gz",
            size_bytes=101,
            content_sha256="b" * 64,
            created_at=created_at,
        )

    with pytest.raises(InputRegistryConflictError, match="concurrently"):
        repository.append_input_file_revision(
            revision(revision_id="inpfr_" + "2" * 32),
            expected_previous_revision_number=2,
        )
    with pytest.raises(InputRegistryConflictError, match="Duplicate"):
        repository.append_input_file_revision(
            revision(revision_id=registered.revision.input_file_revision_id),
            expected_previous_revision_number=1,
        )
    with pytest.raises(ValueError, match="existing InputFile scope"):
        repository.append_input_file_revision(
            revision(
                revision_id="inpfr_" + "3" * 32,
                project_id=OTHER_PROJECT_ID,
            ),
            expected_previous_revision_number=1,
        )
    with pytest.raises(ValueError, match="append exactly once"):
        repository.append_input_file_revision(
            revision(
                revision_id="inpfr_" + "4" * 32,
                revision_number=3,
            ),
            expected_previous_revision_number=1,
        )
    with pytest.raises(ValueError, match="cannot precede history"):
        repository.append_input_file_revision(
            revision(
                revision_id="inpfr_" + "5" * 32,
                created_at=NOW - timedelta(minutes=1),
            ),
            expected_previous_revision_number=1,
        )
    with pytest.raises(KeyError, match="input_file_id"):
        repository.append_input_file_revision(
            revision(
                revision_id="inpfr_" + "6" * 32,
                input_file_id="inpf_" + "9" * 32,
            ),
            expected_previous_revision_number=1,
        )

    assert repository.list_input_file_revisions(
        registered.input_file.input_file_id
    ) == (registered.revision,)


def test_sql_storage_pool_binding_lifecycle_enforces_active_authorization(
    sql_repositories,
) -> None:
    _projects, repository = sql_repositories
    pool = StoragePool(
        storage_pool_id="stgp_" + "1" * 32,
        config_key="ingress-primary",
        display_name="Primary ingress",
        created_at=NOW,
    )

    with pytest.raises(ValueError, match="StoragePool"):
        repository.create_storage_pool(object())
    assert repository.create_storage_pool(pool) == pool
    assert repository.get_storage_pool(pool.storage_pool_id) == pool
    assert repository.list_storage_pools() == (pool,)
    with pytest.raises(InputRegistryConflictError, match="config_key"):
        repository.create_storage_pool(
            StoragePool(
                storage_pool_id="stgp_" + "2" * 32,
                config_key=pool.config_key,
                display_name="Conflicting ingress",
                created_at=NOW,
            )
        )
    with pytest.raises(KeyError, match="storage_pool_id"):
        repository.get_storage_pool("stgp_" + "9" * 32)
    with pytest.raises(KeyError, match="storage_pool_id"):
        repository.archive_storage_pool(
            "stgp_" + "9" * 32,
            archived_at=NOW + timedelta(minutes=1),
        )

    with pytest.raises(ValueError, match="ProjectStoragePoolBinding"):
        repository.create_project_storage_pool_binding(object())
    with pytest.raises(KeyError, match="project_id"):
        repository.create_project_storage_pool_binding(
            ProjectStoragePoolBinding(
                project_id=OTHER_PROJECT_ID,
                storage_pool_id=pool.storage_pool_id,
                bound_at=NOW,
            )
        )
    with pytest.raises(KeyError, match="storage_pool_id"):
        repository.create_project_storage_pool_binding(
            ProjectStoragePoolBinding(
                project_id=PROJECT_ID,
                storage_pool_id="stgp_" + "9" * 32,
                bound_at=NOW,
            )
        )
    with pytest.raises(ValueError, match="cannot precede"):
        repository.create_project_storage_pool_binding(
            ProjectStoragePoolBinding(
                project_id=PROJECT_ID,
                storage_pool_id=pool.storage_pool_id,
                bound_at=NOW - timedelta(microseconds=1),
            )
        )

    binding = ProjectStoragePoolBinding(
        project_id=PROJECT_ID,
        storage_pool_id=pool.storage_pool_id,
        bound_at=NOW,
    )
    assert repository.create_project_storage_pool_binding(binding) == binding
    assert repository.get_project_storage_pool_binding(PROJECT_ID) == binding
    assert repository.get_active_project_storage_pool_binding(PROJECT_ID) == binding
    with pytest.raises(InputRegistryConflictError, match="already bound"):
        repository.create_project_storage_pool_binding(binding)
    with pytest.raises(KeyError, match="StoragePool binding"):
        repository.get_project_storage_pool_binding(OTHER_PROJECT_ID)
    with pytest.raises(KeyError, match="project_id"):
        repository.get_active_project_storage_pool_binding(OTHER_PROJECT_ID)

    archived = repository.archive_storage_pool(
        pool.storage_pool_id,
        archived_at=NOW + timedelta(minutes=1),
    )
    assert repository.list_storage_pools() == ()
    assert repository.list_storage_pools(include_archived=True) == (archived,)
    with pytest.raises(ValueError, match="StoragePool is archived"):
        repository.get_active_project_storage_pool_binding(PROJECT_ID)
    with pytest.raises(InputRegistryConflictError, match="already archived"):
        repository.archive_storage_pool(
            pool.storage_pool_id,
            archived_at=NOW + timedelta(minutes=2),
        )


def test_sql_input_file_lifecycle_rejects_incomplete_or_inactive_evidence(
    sql_repositories,
) -> None:
    _projects, repository = sql_repositories
    service = _service(repository)
    pool = service.create_storage_pool(display_name="Ingress", config_key="ingress")
    service.bind_project_storage_pool(PROJECT_ID, pool.storage_pool_id)
    registered = service.register_input_file(
        PROJECT_ID,
        stable_key="reads-a",
        observation=_observation(
            relative_path="reads/a.fastq.gz",
            size_bytes=100,
            content_sha256="a" * 64,
        ),
    )

    assert (
        repository.get_input_file_by_stable_key(PROJECT_ID, "reads-a")
        == registered.input_file
    )
    assert repository.list_input_files(PROJECT_ID) == (registered.input_file,)
    assert (
        repository.get_input_file_revision(registered.revision.input_file_revision_id)
        == registered.revision
    )
    with pytest.raises(KeyError, match="stable_key"):
        repository.get_input_file_by_stable_key(PROJECT_ID, "unknown")
    with pytest.raises(KeyError, match="project_id"):
        repository.list_input_files(OTHER_PROJECT_ID)
    with pytest.raises(KeyError, match="input_file_id"):
        repository.archive_input_file(
            "inpf_" + "9" * 32,
            archived_at=NOW + timedelta(minutes=1),
        )
    with pytest.raises(KeyError, match="input_file_revision_id"):
        repository.get_input_file_revision("inpfr_" + "9" * 32)
    with pytest.raises(KeyError, match="input_file_id"):
        repository.list_input_file_revisions("inpf_" + "9" * 32)

    candidate = InputFile(
        input_file_id="inpf_" + "8" * 32,
        project_id=PROJECT_ID,
        storage_pool_id=pool.storage_pool_id,
        stable_key="reads-b",
        created_at=NOW + timedelta(minutes=1),
    )
    valid_initial = build_input_file_revision(
        input_file_revision_id="inpfr_" + "8" * 32,
        input_file_id=candidate.input_file_id,
        project_id=candidate.project_id,
        storage_pool_id=candidate.storage_pool_id,
        revision_number=1,
        relative_path="reads/b.fastq.gz",
        size_bytes=200,
        content_sha256="b" * 64,
        created_at=candidate.created_at,
    )
    with pytest.raises(ValueError, match="InputFile"):
        repository.create_input_file(object(), valid_initial)
    with pytest.raises(ValueError, match="InputFileRevision"):
        repository.create_input_file(candidate, object())
    with pytest.raises(ValueError, match="revision 1"):
        repository.create_input_file(
            candidate,
            build_input_file_revision(
                input_file_revision_id="inpfr_" + "7" * 32,
                input_file_id=candidate.input_file_id,
                project_id=candidate.project_id,
                storage_pool_id=candidate.storage_pool_id,
                revision_number=2,
                relative_path="reads/b.fastq.gz",
                size_bytes=200,
                content_sha256="b" * 64,
                created_at=candidate.created_at,
            ),
        )
    with pytest.raises(ValueError, match="cannot precede"):
        repository.create_input_file(
            candidate,
            build_input_file_revision(
                input_file_revision_id="inpfr_" + "7" * 32,
                input_file_id=candidate.input_file_id,
                project_id=candidate.project_id,
                storage_pool_id=candidate.storage_pool_id,
                revision_number=1,
                relative_path="reads/b.fastq.gz",
                size_bytes=200,
                content_sha256="b" * 64,
                created_at=NOW,
            ),
        )
    with pytest.raises(ValueError, match="revision"):
        repository.append_input_file_revision(
            object(),
            expected_previous_revision_number=1,
        )
    for invalid_previous in (True, 0):
        with pytest.raises(ValueError, match="positive"):
            repository.append_input_file_revision(
                valid_initial,
                expected_previous_revision_number=invalid_previous,
            )

    archived = repository.archive_input_file(
        registered.input_file.input_file_id,
        archived_at=NOW + timedelta(minutes=1),
    )
    assert repository.list_input_files(PROJECT_ID) == ()
    assert repository.list_input_files(
        PROJECT_ID,
        include_archived=True,
    ) == (archived,)
    with pytest.raises(ValueError, match="InputFile is archived"):
        repository.append_input_file_revision(
            build_input_file_revision(
                input_file_revision_id="inpfr_" + "6" * 32,
                input_file_id=registered.input_file.input_file_id,
                project_id=PROJECT_ID,
                storage_pool_id=pool.storage_pool_id,
                revision_number=2,
                relative_path="reads/a-v2.fastq.gz",
                size_bytes=101,
                content_sha256="c" * 64,
                created_at=NOW + timedelta(minutes=2),
            ),
            expected_previous_revision_number=1,
        )
    with pytest.raises(InputRegistryConflictError, match="already archived"):
        repository.archive_input_file(
            registered.input_file.input_file_id,
            archived_at=NOW + timedelta(minutes=2),
        )
    assert repository.list_input_file_revisions(
        registered.input_file.input_file_id
    ) == (registered.revision,)
