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
    InputProvenanceMode,
    InputUseBindingPlan,
    PlannedInputUse,
)
from encode_pipeline.services.data_registry_repositories import (
    InMemoryDataRegistryRepository,
)
from encode_pipeline.services.input_file_access import FileObservation
from encode_pipeline.services.input_registry import InputRegistryService
from encode_pipeline.services.input_registry_repositories import (
    InMemoryInputRegistryRepository,
)


NOW = datetime(2026, 7, 26, 9, 0, tzinfo=timezone.utc)
PROJECT_ID = "prj_11111111111111111111111111111111"


def _ids(prefix: str):
    serial = count(1)
    return lambda: f"{prefix}_{next(serial):032x}"


def _project() -> Project:
    return Project(
        project_id=PROJECT_ID,
        display_name="Project One",
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
) -> InputUseBindingPlan:
    input_uses = []
    if revision_id is not None:
        input_uses.append(
            PlannedInputUse(
                key="reads",
                occurrence=0,
                capability_version="reads-v1",
                closure_contract_version="regular_file_v1",
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
        project_id=PROJECT_ID,
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
