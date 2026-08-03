"""SQLAlchemy/in-memory Reference Profile repository parity."""

from __future__ import annotations

from datetime import datetime, timezone

import pytest

from encode_pipeline.persistence import (
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.persistence.reference_profiles import (
    SqlAlchemyReferenceProfileRepository,
)
from encode_pipeline.platform.reference_profiles import (
    AdapterReferenceBindingIdentity,
    ReferenceProfile,
    ReferenceProfileWorkflowBinding,
    build_reference_profile_revision,
)
from encode_pipeline.services.reference_profile_repositories import (
    InMemoryReferenceProfileRepository,
    ReferenceProfileConflictError,
)


NOW = datetime(2026, 8, 3, 9, 0, tzinfo=timezone.utc)
PROFILE_ID = "refp_" + "1" * 32
REVISION_ID = "refpr_" + "2" * 32


def _profile() -> ReferenceProfile:
    return ReferenceProfile(
        profile_id=PROFILE_ID,
        safe_key="grch38",
        created_at=NOW,
    )


def _revision(revision_number: int, revision_id: str):
    return build_reference_profile_revision(
        revision_id=revision_id,
        profile_id=PROFILE_ID,
        revision_number=revision_number,
        display_name=f"GRCh38 r{revision_number}",
        organism="Homo sapiens",
        assembly="GRCh38",
        config_key="grch38-private",
        workflow_bindings=(
            ReferenceProfileWorkflowBinding.from_adapter_identity(
                AdapterReferenceBindingIdentity(
                    workflow_id="bulk-rnaseq",
                    contract_version="bulk-reference-v1",
                    identity_sha256=f"{revision_number}" * 64,
                )
            ),
        ),
        created_at=NOW,
    )


@pytest.fixture
def repositories(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'profiles.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    sql = SqlAlchemyReferenceProfileRepository(create_session_factory(engine))
    yield sql, InMemoryReferenceProfileRepository()
    engine.dispose()


def test_sql_and_memory_reference_profile_lifecycle_parity(repositories) -> None:
    results = []
    for repository in repositories:
        profile, first = repository.create_profile(
            _profile(), _revision(1, REVISION_ID)
        )
        second = repository.append_revision(
            _revision(2, "refpr_" + "3" * 32),
            expected_previous_revision_number=1,
        )
        enabled = repository.set_enabled_revision(PROFILE_ID, second.revision_id)
        results.append(
            (
                profile,
                first,
                second,
                enabled,
                repository.get_profile_by_safe_key("grch38"),
                repository.get_profile_for_revision(second.revision_id),
                repository.list_revisions(PROFILE_ID),
                repository.list_enabled_for_workflow("bulk-rnaseq"),
            )
        )

    assert results[0] == results[1]
    assert results[0][3].enabled_revision_id == "refpr_" + "3" * 32


def test_repository_rejects_duplicate_key_cross_profile_enable_and_revision_race(
    repositories,
) -> None:
    for repository in repositories:
        repository.create_profile(_profile(), _revision(1, REVISION_ID))
        with pytest.raises(ReferenceProfileConflictError):
            repository.create_profile(_profile(), _revision(1, REVISION_ID))
        with pytest.raises(ReferenceProfileConflictError):
            repository.append_revision(
                _revision(2, "refpr_" + "3" * 32),
                expected_previous_revision_number=2,
            )
