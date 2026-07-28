"""Tests for operator-owned StoragePool and InputFile registry orchestration."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from itertools import count
import stat

import pytest

from encode_pipeline.platform.data_registry import (
    LEGACY_PROJECT_ID,
    Project,
    ProjectKind,
)
from encode_pipeline.platform.input_registry import (
    InputProvenanceMode,
    InputUseBindingPlan,
    PlannedInputUse,
)
from encode_pipeline.services.data_registry_repositories import (
    InMemoryDataRegistryRepository,
)
from encode_pipeline.services.input_file_access import FileObservation
from encode_pipeline.services.input_file_access import InputFileAccess
from encode_pipeline.services.input_registry import (
    InputRegistryAdminService,
    InputRegistryService,
)
from encode_pipeline.services.input_registry_repositories import (
    InMemoryInputRegistryRepository,
    InputRegistryConflictError,
)
from encode_pipeline.services.private_storage_pools import PrivateStoragePoolConfig


NOW = datetime(2026, 7, 26, 9, 0, tzinfo=timezone.utc)
PROJECT_ID = "prj_11111111111111111111111111111111"


def _ids(prefix: str):
    serial = count(1)
    return lambda: f"{prefix}_{next(serial):032x}"


def _service() -> tuple[InputRegistryService, InMemoryDataRegistryRepository]:
    projects = InMemoryDataRegistryRepository(legacy_created_at=NOW)
    projects.create_project(
        Project(
            project_id=PROJECT_ID,
            display_name="Project One",
            kind=ProjectKind.USER,
            created_at=NOW,
        )
    )
    repository = InMemoryInputRegistryRepository(project_repository=projects)
    return (
        InputRegistryService(
            repository=repository,
            storage_pool_id_factory=_ids("stgp"),
            input_file_id_factory=_ids("inpf"),
            input_file_revision_id_factory=_ids("inpfr"),
            now_factory=lambda: NOW,
        ),
        projects,
    )


def _observation(
    *,
    relative_path: str = "reads/donor-01.fastq.gz",
    size_bytes: int = 123,
    content_sha256: str = "a" * 64,
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


def test_storage_pool_create_list_archive_never_stores_a_root() -> None:
    service, _projects = _service()

    pool = service.create_storage_pool(
        display_name="  Approved ingress  ",
        config_key="ingress-primary",
    )

    assert pool.storage_pool_id == "stgp_" + "1".zfill(32)
    assert pool.display_name == "Approved ingress"
    assert pool.config_key == "ingress-primary"
    assert pool.created_at == NOW
    assert pool.archived_at is None
    assert not hasattr(pool, "root")
    assert not hasattr(pool, "path")
    assert service.list_storage_pools() == (pool,)

    archived = service.archive_storage_pool(pool.storage_pool_id)

    assert archived.archived_at == NOW
    assert service.list_storage_pools() == ()
    assert service.list_storage_pools(include_archived=True) == (archived,)
    with pytest.raises(InputRegistryConflictError, match="already archived"):
        service.archive_storage_pool(pool.storage_pool_id)


def test_project_binds_once_to_an_active_approved_pool() -> None:
    service, projects = _service()
    pool = service.create_storage_pool(
        display_name="Ingress",
        config_key="ingress",
    )

    binding = service.bind_project_storage_pool(PROJECT_ID, pool.storage_pool_id)

    assert binding.project_id == PROJECT_ID
    assert binding.storage_pool_id == pool.storage_pool_id
    assert binding.bound_at == NOW
    assert service.get_project_storage_pool_binding(PROJECT_ID) == binding
    with pytest.raises(InputRegistryConflictError, match="already bound"):
        service.bind_project_storage_pool(PROJECT_ID, pool.storage_pool_id)
    with pytest.raises(ValueError, match="Legacy"):
        service.bind_project_storage_pool(LEGACY_PROJECT_ID, pool.storage_pool_id)

    archived_project = Project(
        project_id="prj_22222222222222222222222222222222",
        display_name="Archived",
        kind=ProjectKind.USER,
        created_at=NOW,
        archived_at=NOW + timedelta(minutes=1),
    )
    projects.create_project(archived_project)
    with pytest.raises(ValueError, match="archived"):
        service.bind_project_storage_pool(
            archived_project.project_id,
            pool.storage_pool_id,
        )

    other_pool = service.create_storage_pool(
        display_name="Archived ingress",
        config_key="archived-ingress",
    )
    service.archive_storage_pool(other_pool.storage_pool_id)
    unbound = Project(
        project_id="prj_33333333333333333333333333333333",
        display_name="Unbound",
        kind=ProjectKind.USER,
        created_at=NOW,
    )
    projects.create_project(unbound)
    with pytest.raises(ValueError, match="archived"):
        service.bind_project_storage_pool(
            unbound.project_id,
            other_pool.storage_pool_id,
        )


def test_register_is_idempotent_and_changes_append_an_immutable_revision() -> None:
    service, _projects = _service()
    pool = service.create_storage_pool(
        display_name="Ingress",
        config_key="ingress",
    )
    service.bind_project_storage_pool(PROJECT_ID, pool.storage_pool_id)

    created = service.register_input_file(
        PROJECT_ID,
        stable_key="donor-01-r1",
        observation=_observation(),
    )
    repeated = service.register_input_file(
        PROJECT_ID,
        stable_key="donor-01-r1",
        observation=_observation(),
    )
    revised = service.register_input_file(
        PROJECT_ID,
        stable_key="donor-01-r1",
        observation=_observation(
            relative_path="reads/replaced/donor-01.fastq.gz",
            size_bytes=456,
            content_sha256="b" * 64,
        ),
    )

    assert created.input_file.input_file_id == "inpf_" + "1".zfill(32)
    assert created.revision.input_file_revision_id == "inpfr_" + "1".zfill(32)
    assert created.revision.revision_number == 1
    assert repeated == created
    assert revised.input_file == created.input_file
    assert revised.revision.input_file_revision_id == "inpfr_" + "2".zfill(32)
    assert revised.revision.revision_number == 2
    assert revised.revision.relative_path == "reads/replaced/donor-01.fastq.gz"
    assert revised.revision.size_bytes == 456
    assert revised.revision.content_sha256 == "b" * 64
    assert service.list_input_file_revisions(created.input_file.input_file_id) == (
        created.revision,
        revised.revision,
    )
    assert service.list_input_files(PROJECT_ID) == (revised,)


def test_registration_requires_active_project_pool_binding_and_input() -> None:
    service, projects = _service()

    with pytest.raises(KeyError, match="StoragePool binding"):
        service.register_input_file(
            PROJECT_ID,
            stable_key="unbound",
            observation=_observation(),
        )

    pool = service.create_storage_pool(display_name="Ingress", config_key="ingress")
    service.bind_project_storage_pool(PROJECT_ID, pool.storage_pool_id)
    registered = service.register_input_file(
        PROJECT_ID,
        stable_key="bound",
        observation=_observation(),
    )
    service.archive_input_file(registered.input_file.input_file_id)

    assert service.list_input_files(PROJECT_ID) == ()
    assert (
        service.list_input_files(PROJECT_ID, include_archived=True)[
            0
        ].input_file.archived_at
        == NOW
    )
    with pytest.raises(ValueError, match="archived"):
        service.register_input_file(
            PROJECT_ID,
            stable_key="bound",
            observation=_observation(content_sha256="c" * 64),
        )

    projects.archive_project(PROJECT_ID, archived_at=NOW + timedelta(minutes=1))
    with pytest.raises(ValueError, match="archived"):
        service.register_input_file(
            PROJECT_ID,
            stable_key="bound",
            observation=_observation(),
        )


def test_pool_archival_blocks_new_registration_without_invalidating_history() -> None:
    service, _projects = _service()
    pool = service.create_storage_pool(display_name="Ingress", config_key="ingress")
    service.bind_project_storage_pool(PROJECT_ID, pool.storage_pool_id)
    registered = service.register_input_file(
        PROJECT_ID,
        stable_key="registered",
        observation=_observation(),
    )

    service.archive_storage_pool(pool.storage_pool_id)

    assert service.list_input_file_revisions(registered.input_file.input_file_id) == (
        registered.revision,
    )
    with pytest.raises(ValueError, match="archived"):
        service.register_input_file(
            PROJECT_ID,
            stable_key="new",
            observation=_observation(relative_path="reads/new.fastq.gz"),
        )


def test_archived_input_or_pool_cannot_enter_a_new_binding_envelope() -> None:
    service, _projects = _service()
    pool = service.create_storage_pool(display_name="Ingress", config_key="ingress")
    service.bind_project_storage_pool(PROJECT_ID, pool.storage_pool_id)
    registered = service.register_input_file(
        PROJECT_ID,
        stable_key="registered",
        observation=_observation(),
    )
    plan = InputUseBindingPlan(
        project_id=PROJECT_ID,
        workflow_id="test-workflow",
        adapter_contract_version="adapter-v2",
        input_uses=(
            PlannedInputUse(
                key="reads",
                occurrence=0,
                capability_version="reads-v1",
                closure_contract_version="regular_file_v1",
                provenance_mode=InputProvenanceMode.MANAGED_REVISION_V1,
                input_file_revision_ids=(registered.revision.input_file_revision_id,),
            ),
        ),
    )

    service.archive_input_file(registered.input_file.input_file_id)

    with pytest.raises(ValueError, match="InputFile is archived"):
        service.resolve_input_use_binding_plan(
            plan,
            project_sample_binding_digest="b" * 64,
            workflow_inputs_digest="c" * 64,
        )

    other_service, _other_projects = _service()
    other_pool = other_service.create_storage_pool(
        display_name="Ingress",
        config_key="ingress",
    )
    other_service.bind_project_storage_pool(PROJECT_ID, other_pool.storage_pool_id)
    other_registered = other_service.register_input_file(
        PROJECT_ID,
        stable_key="registered",
        observation=_observation(),
    )
    other_plan = InputUseBindingPlan(
        project_id=PROJECT_ID,
        workflow_id="test-workflow",
        adapter_contract_version="adapter-v2",
        input_uses=(
            PlannedInputUse(
                key="reads",
                occurrence=0,
                capability_version="reads-v1",
                closure_contract_version="regular_file_v1",
                provenance_mode=InputProvenanceMode.MANAGED_REVISION_V1,
                input_file_revision_ids=(
                    other_registered.revision.input_file_revision_id,
                ),
            ),
        ),
    )
    other_service.archive_storage_pool(other_pool.storage_pool_id)

    with pytest.raises(ValueError, match="StoragePool is archived"):
        other_service.resolve_input_use_binding_plan(
            other_plan,
            project_sample_binding_digest="b" * 64,
            workflow_inputs_digest="c" * 64,
        )


def test_admin_registration_observes_private_bytes_but_returns_no_path(
    tmp_path,
) -> None:
    service, _projects = _service()
    private_root = tmp_path / "private-inputs"
    private_root.mkdir()
    reads = private_root / "reads"
    reads.mkdir()
    source = reads / "donor.fastq.gz"
    source.write_bytes(b"private-read-bytes")
    admin = InputRegistryAdminService(
        service,
        private_config=PrivateStoragePoolConfig({"ingress": private_root}),
    )
    pool = admin.register_storage_pool(
        display_name="Ingress",
        config_key="ingress",
    )
    admin.bind_project_storage_pool(
        project_id=PROJECT_ID,
        storage_pool_id=pool.storage_pool_id,
    )

    result = admin.register_input_file(
        project_id=PROJECT_ID,
        stable_key="donor-r1",
        pool_relative_path="reads/donor.fastq.gz",
    )

    assert result["input_file_revision_id"] == "inpfr_" + "1".zfill(32)
    assert result["size_bytes"] == len(b"private-read-bytes")
    assert str(private_root) not in repr(result)
    assert "reads/donor.fastq.gz" not in repr(result)


def test_admin_rejects_invalid_stable_key_before_observing_bytes(
    tmp_path,
    monkeypatch,
) -> None:
    service, _projects = _service()
    private_root = tmp_path / "private-inputs"
    private_root.mkdir()
    (private_root / "large.fastq").write_bytes(b"private-read-bytes")
    admin = InputRegistryAdminService(
        service,
        private_config=PrivateStoragePoolConfig({"ingress": private_root}),
    )
    pool = admin.register_storage_pool(display_name="Ingress", config_key="ingress")
    admin.bind_project_storage_pool(
        project_id=PROJECT_ID,
        storage_pool_id=pool.storage_pool_id,
    )

    def unexpected_observe(
        _access: InputFileAccess,
        _relative_path: str,
    ) -> FileObservation:
        raise AssertionError("file bytes must not be observed")

    monkeypatch.setattr(InputFileAccess, "observe", unexpected_observe)

    with pytest.raises(ValueError, match="stable_key"):
        admin.register_input_file(
            project_id=PROJECT_ID,
            stable_key="../invalid",
            pool_relative_path="large.fastq",
        )


def test_admin_rejects_archived_project_before_observing_bytes(
    tmp_path,
    monkeypatch,
) -> None:
    service, projects = _service()
    private_root = tmp_path / "private-inputs"
    private_root.mkdir()
    (private_root / "large.fastq").write_bytes(b"private-read-bytes")
    admin = InputRegistryAdminService(
        service,
        private_config=PrivateStoragePoolConfig({"ingress": private_root}),
    )
    pool = admin.register_storage_pool(display_name="Ingress", config_key="ingress")
    admin.bind_project_storage_pool(
        project_id=PROJECT_ID,
        storage_pool_id=pool.storage_pool_id,
    )
    projects.archive_project(PROJECT_ID, archived_at=NOW + timedelta(minutes=1))

    def unexpected_observe(
        _access: InputFileAccess,
        _relative_path: str,
    ) -> FileObservation:
        raise AssertionError("file bytes must not be observed")

    monkeypatch.setattr(InputFileAccess, "observe", unexpected_observe)

    with pytest.raises(ValueError, match="archived"):
        admin.register_input_file(
            project_id=PROJECT_ID,
            stable_key="donor-r1",
            pool_relative_path="large.fastq",
        )
