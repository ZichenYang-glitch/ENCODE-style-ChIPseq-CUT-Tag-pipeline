"""Tests for the initial Alembic platform schema."""

from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import stat
import tempfile

import pytest
from sqlalchemy import inspect, text
from sqlalchemy.exc import IntegrityError

from encode_pipeline.persistence import migration_admission
from encode_pipeline.persistence import (
    SqlAlchemyRunRepository,
    create_database_engine,
    create_session_factory,
)
from encode_pipeline.persistence.migrations import (
    downgrade_database,
    upgrade_database,
)
from encode_pipeline.platform.data_registry import (
    build_legacy_project_sample_binding,
)
from encode_pipeline.platform.input_registry import InputUseBindingEnvelope
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.platform.snapshots import build_workflow_inputs_digest
from encode_pipeline.services.run_repositories import RunEventDraft


pytestmark = pytest.mark.full_main


LEGACY_PROJECT_ID = "prj_00000000000000000000000000000000"
REV08_CANONICAL_INPUTS = json.dumps(
    {"config": {}, "options": {}, "samples": None},
    sort_keys=True,
    separators=(",", ":"),
)
REV08_INPUTS_DIGEST = build_workflow_inputs_digest(REV08_CANONICAL_INPUTS)
STANDALONE_CANONICAL_INPUTS = json.dumps(
    {"config": {"label": "样本"}, "options": {}, "samples": None},
    ensure_ascii=False,
    sort_keys=True,
    separators=(",", ":"),
)
STANDALONE_NONCANONICAL_INPUTS = (
    '{ "samples": null, "options": {}, "config": {"label": "样本"} }'
)
INPUT_BINDING_DIGEST_SCHEME = "sha256-framed-input-use-binding-envelope-v1"


INPUT_REGISTRY_TABLES = {
    "input_file_revisions",
    "input_files",
    "project_storage_pool_bindings",
    "run_input_bindings",
    "run_input_members",
    "run_input_uses",
    "snapshot_input_bindings",
    "snapshot_input_members",
    "snapshot_input_uses",
    "storage_pools",
}


EXPECTED_TABLES = {
    "alembic_version",
    "artifact_publications",
    "auth_sessions",
    "projects",
    "input_file_revisions",
    "input_files",
    "project_storage_pool_bindings",
    "reference_profile_revisions",
    "reference_profile_workflow_bindings",
    "reference_profiles",
    "run_artifacts",
    "run_events",
    "run_execution_assignments",
    "run_input_bindings",
    "run_input_members",
    "run_input_uses",
    "run_logs",
    "run_project_bindings",
    "run_qc_metrics",
    "run_reference_bindings",
    "run_result_attempts",
    "run_result_states",
    "run_samples",
    "run_workflow_build_identities",
    "runs",
    "sample_revisions",
    "samples",
    "snapshot_project_bindings",
    "snapshot_reference_bindings",
    "snapshot_input_bindings",
    "snapshot_input_members",
    "snapshot_input_uses",
    "snapshot_sample_revisions",
    "security_audit_events",
    "storage_pools",
    "user_accounts",
    "validated_input_snapshots",
}


class _CoverageMigrationTemporaryFiles:
    def __init__(self, root: Path) -> None:
        self._root = root

    def TemporaryDirectory(  # noqa: N802 - mirrors the standard-library API
        self,
        *,
        prefix: str,
        ignore_cleanup_errors: bool,
    ) -> tempfile.TemporaryDirectory[str]:
        if prefix != "helixweave-migration-snapshot-":
            raise AssertionError("unexpected migration snapshot prefix")
        return tempfile.TemporaryDirectory(
            prefix=prefix,
            ignore_cleanup_errors=ignore_cleanup_errors,
            dir=self._root,
        )


@pytest.fixture(autouse=True)
def _route_admitted_migrations_to_controlled_coverage_root(monkeypatch):
    configured_root = os.environ.get("HELIXWEAVE_COVERAGE_MIGRATION_ROOT")
    if configured_root is None:
        yield
        return

    root = Path(configured_root)
    try:
        resolved_root = root.resolve(strict=True)
        root_mode = stat.S_IMODE(resolved_root.stat().st_mode)
        existing_entries = tuple(resolved_root.iterdir())
    except (OSError, RuntimeError):
        pytest.fail("migration coverage root must be an existing private directory")
    if (
        not root.is_absolute()
        or root != resolved_root
        or root.is_symlink()
        or not resolved_root.is_dir()
        or root_mode != 0o700
        or not root.name.startswith("helixweave-coverage-migrations-")
        or existing_entries
    ):
        pytest.fail("migration coverage root must be an empty private directory")

    monkeypatch.setattr(
        migration_admission,
        "tempfile",
        _CoverageMigrationTemporaryFiles(resolved_root),
    )
    yield
    assert tuple(resolved_root.iterdir()) == ()


def test_initial_migration_creates_versioned_run_schema(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"

    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    inspector = inspect(engine)

    assert set(inspector.get_table_names()) == EXPECTED_TABLES
    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260827_15"
        )
        assert connection.scalar(text("PRAGMA foreign_keys")) == 1
        assert connection.scalar(text("PRAGMA journal_mode")) == "wal"

    event_constraints = {
        constraint["name"]
        for constraint in inspector.get_unique_constraints("run_events")
    }
    log_constraints = {
        constraint["name"]
        for constraint in inspector.get_unique_constraints("run_logs")
    }
    assignment_constraints = {
        constraint["name"]
        for constraint in inspector.get_unique_constraints("run_execution_assignments")
    }
    assert "uq_run_events_run_sequence" in event_constraints
    assert "uq_run_logs_run_stream_sequence" in log_constraints
    assert "uq_run_execution_assignments_job_id" in assignment_constraints
    assert {
        constraint["name"]
        for constraint in inspector.get_check_constraints("run_execution_assignments")
    } == {
        "ck_run_execution_assignments_ack_requires_request",
        "ck_run_execution_assignments_claim_requires_dispatch",
        "ck_run_execution_assignments_cleanup_endpoint_format",
        "ck_run_execution_assignments_cleanup_identity_pair",
        "ck_run_execution_assignments_cleanup_scope_format",
        "ck_run_execution_assignments_request_reason_pair",
        "ck_run_execution_assignments_request_requires_claim",
        "ck_run_execution_assignments_requeue_requires_dispatch",
        "ck_run_execution_assignments_requeue_confirm_requires_request",
        "ck_run_execution_assignments_requeue_confirmation_order",
    }
    snapshot_foreign_key = inspector.get_foreign_keys("validated_input_snapshots")[0]
    assert snapshot_foreign_key["referred_table"] == "runs"
    assert snapshot_foreign_key["options"] == {"ondelete": "RESTRICT"}
    assert {
        constraint["name"]
        for constraint in inspector.get_unique_constraints("validated_input_snapshots")
    } == {"uq_validated_input_snapshots_consumed_run_id"}
    assert {
        constraint["name"]
        for constraint in inspector.get_check_constraints("validated_input_snapshots")
    } == {
        "ck_validated_input_snapshots_consumption_pair",
        "ck_validated_input_snapshots_digest_length",
        "ck_validated_input_snapshots_expiry",
        "ck_validated_input_snapshots_success",
    }
    assert inspector.get_foreign_keys("run_events")[0]["options"] == {
        "ondelete": "CASCADE"
    }
    assignment_foreign_key = inspector.get_foreign_keys("run_execution_assignments")[0]
    assert assignment_foreign_key["referred_table"] == "runs"
    assert assignment_foreign_key["options"] == {"ondelete": "CASCADE"}
    assert inspector.get_pk_constraint("run_execution_assignments")[
        "constrained_columns"
    ] == ["run_id"]
    build_foreign_key = inspector.get_foreign_keys("run_workflow_build_identities")[0]
    assert build_foreign_key["referred_table"] == "runs"
    assert build_foreign_key["options"] == {"ondelete": "CASCADE"}
    assert inspector.get_pk_constraint("run_workflow_build_identities")[
        "constrained_columns"
    ] == ["run_id"]
    qc_foreign_key = inspector.get_foreign_keys("run_qc_metrics")[0]
    assert qc_foreign_key["referred_table"] == "runs"
    assert qc_foreign_key["options"] == {"ondelete": "CASCADE"}
    assert {
        constraint["name"]
        for constraint in inspector.get_unique_constraints("run_qc_metrics")
    } == {"uq_run_qc_metrics_run_metric"}
    assert {
        constraint["name"]
        for constraint in inspector.get_check_constraints("run_qc_metrics")
    } == {
        "ck_run_qc_metrics_flag",
        "ck_run_qc_metrics_scope",
        "ck_run_qc_metrics_scope_identifiers",
        "ck_run_qc_metrics_value_text_length",
    }
    assert inspector.get_pk_constraint("run_result_states")["constrained_columns"] == [
        "run_id"
    ]
    assert inspector.get_foreign_keys("run_result_states")[0]["options"] == {
        "ondelete": "CASCADE"
    }
    assert {
        constraint["name"]
        for constraint in inspector.get_check_constraints("run_result_states")
    } == {
        "ck_run_result_states_artifact_attempt",
        "ck_run_result_states_artifact_binding",
        "ck_run_result_states_artifact_outcome",
        "ck_run_result_states_artifact_reason",
        "ck_run_result_states_nonnegative_revisions",
        "ck_run_result_states_qc_attempt",
        "ck_run_result_states_qc_binding",
        "ck_run_result_states_qc_outcome",
        "ck_run_result_states_qc_reason",
    }
    assert inspector.get_pk_constraint("run_result_attempts")[
        "constrained_columns"
    ] == ["attempt_id"]
    assert inspector.get_foreign_keys("run_result_attempts")[0]["options"] == {
        "ondelete": "CASCADE"
    }
    assert {
        constraint["name"]
        for constraint in inspector.get_check_constraints("run_result_attempts")
    } == {"ck_run_result_attempts_binding"}
    assert {
        index["name"] for index in inspector.get_indexes("run_result_attempts")
    } == {"ix_run_result_attempts_run_id"}
    assert {column["name"] for column in inspector.get_columns("run_artifacts")} >= {
        "revision"
    }
    engine.dispose()


def test_initial_migration_can_downgrade_and_reapply(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)

    downgrade_database(database_url, "base")
    downgraded_engine = create_database_engine(database_url)
    assert inspect(downgraded_engine).get_table_names() == ["alembic_version"]
    downgraded_engine.dispose()

    upgrade_database(database_url)
    upgraded_engine = create_database_engine(database_url)
    assert set(inspect(upgraded_engine).get_table_names()) == EXPECTED_TABLES
    upgraded_engine.dispose()


def test_execution_assignment_migration_can_downgrade_independently(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)

    downgrade_database(database_url, "20260711_01")
    downgraded_engine = create_database_engine(database_url)
    assert (
        "run_execution_assignments" not in inspect(downgraded_engine).get_table_names()
    )
    with downgraded_engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260711_01"
        )
    downgraded_engine.dispose()

    upgrade_database(database_url)
    upgraded_engine = create_database_engine(database_url)
    assert "run_execution_assignments" in inspect(upgraded_engine).get_table_names()
    upgraded_engine.dispose()


def test_build_identity_migration_does_not_backfill_legacy_planned_runs(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url, "20260711_02")
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        connection.execute(
            text(
                """
                INSERT INTO runs (
                    run_id, workflow_id, inputs, status, created_at, updated_at,
                    started_at, ended_at, current_stage, cancellation_reason,
                    error, tags
                ) VALUES (
                    :run_id, :workflow_id, :inputs, 'planned', :created_at,
                    :updated_at, NULL, NULL, 'preflight', NULL, NULL, :tags
                )
                """
            ),
            {
                "run_id": "legacy-planned",
                "workflow_id": "encode-style-chipseq-cuttag-atac-mnase",
                "inputs": '{"config":{},"samples":null,"options":{}}',
                "created_at": "2026-07-12 00:00:00",
                "updated_at": "2026-07-12 00:00:00",
                "tags": "{}",
            },
        )
    engine.dispose()

    upgrade_database(database_url)
    upgraded = create_database_engine(database_url)
    with upgraded.connect() as connection:
        assert (
            connection.scalar(
                text("SELECT status FROM runs WHERE run_id = 'legacy-planned'")
            )
            == "planned"
        )
        assert (
            connection.scalar(
                text("SELECT count(*) FROM run_workflow_build_identities")
            )
            == 0
        )
    upgraded.dispose()


def test_project_sample_registry_upgrades_rev08_with_conservative_legacy_bindings(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'project-sample-upgrade.db'}"
    upgrade_database(database_url, "20260717_08")
    engine = create_database_engine(database_url)
    snapshot_id = "vsnap_" + "a" * 32
    unconsumed_snapshot_id = "vsnap_" + "b" * 32
    digest = REV08_INPUTS_DIGEST
    with engine.begin() as connection:
        _insert_rev08_run(connection, "consumed-run", created_at="2026-07-17 01:00:00")
        _insert_rev08_run(
            connection,
            "standalone-run",
            created_at="2026-07-17 02:00:00",
            inputs=STANDALONE_NONCANONICAL_INPUTS,
        )
        _insert_rev08_snapshot(
            connection,
            snapshot_id,
            digest=digest,
            consumed_run_id="consumed-run",
        )
        _insert_rev08_snapshot(
            connection,
            unconsumed_snapshot_id,
            digest=REV08_INPUTS_DIGEST,
            consumed_run_id=None,
        )
        runs_before = (
            connection.execute(text("SELECT * FROM runs ORDER BY run_id"))
            .mappings()
            .all()
        )
        snapshots_before = (
            connection.execute(
                text("SELECT * FROM validated_input_snapshots ORDER BY snapshot_id")
            )
            .mappings()
            .all()
        )
    engine.dispose()

    upgrade_database(database_url)
    upgraded = create_database_engine(database_url)
    with upgraded.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260827_15"
        )
        legacy = (
            connection.execute(
                text(
                    "SELECT project_id, display_name, kind, archived_at "
                    "FROM projects WHERE project_id=:project_id"
                ),
                {"project_id": LEGACY_PROJECT_ID},
            )
            .mappings()
            .one()
        )
        assert dict(legacy) == {
            "project_id": LEGACY_PROJECT_ID,
            "display_name": "Legacy Project",
            "kind": "system",
            "archived_at": None,
        }

        snapshot_bindings = (
            connection.execute(
                text(
                    "SELECT snapshot_id, project_id, binding_mode, provenance, "
                    "workflow_inputs_digest, binding_digest_scheme, binding_digest "
                    "FROM snapshot_project_bindings ORDER BY snapshot_id"
                )
            )
            .mappings()
            .all()
        )
        run_bindings = (
            connection.execute(
                text(
                    "SELECT run_id, project_id, binding_mode, provenance, "
                    "workflow_inputs_digest, binding_digest_scheme, binding_digest "
                    "FROM run_project_bindings ORDER BY run_id"
                )
            )
            .mappings()
            .all()
        )
        assert len(snapshot_bindings) == 2
        assert len(run_bindings) == 2
        assert {row["project_id"] for row in snapshot_bindings + run_bindings} == {
            LEGACY_PROJECT_ID
        }
        assert {row["binding_mode"] for row in snapshot_bindings + run_bindings} == {
            "legacy_v1"
        }
        assert {row["provenance"] for row in snapshot_bindings + run_bindings} == {
            "unresolved"
        }
        assert all(len(row["binding_digest"]) == 64 for row in snapshot_bindings)
        assert all(len(row["binding_digest"]) == 64 for row in run_bindings)
        consumed_snapshot = next(
            row for row in snapshot_bindings if row["snapshot_id"] == snapshot_id
        )
        consumed_run = next(
            row for row in run_bindings if row["run_id"] == "consumed-run"
        )
        expected_binding = build_legacy_project_sample_binding(digest)
        assert (
            consumed_run["binding_digest"]
            == consumed_snapshot["binding_digest"]
            == expected_binding.digest
        )
        assert (
            consumed_run["binding_digest_scheme"]
            == consumed_snapshot["binding_digest_scheme"]
            == expected_binding.digest_scheme
        )
        assert (
            consumed_run["workflow_inputs_digest"]
            == consumed_snapshot["workflow_inputs_digest"]
            == digest
        )

        assert connection.scalar(text("SELECT count(*) FROM samples")) == 0
        assert connection.scalar(text("SELECT count(*) FROM sample_revisions")) == 0
        assert (
            connection.scalar(text("SELECT count(*) FROM snapshot_sample_revisions"))
            == 0
        )
        assert connection.scalar(text("SELECT count(*) FROM run_samples")) == 0
        runs_after = (
            connection.execute(text("SELECT * FROM runs ORDER BY run_id"))
            .mappings()
            .all()
        )
        assert all(row["requested_by_user_id"] is None for row in runs_after)
        assert [
            {key: value for key, value in row.items() if key != "requested_by_user_id"}
            for row in runs_after
        ] == runs_before
        assert (
            connection.execute(
                text("SELECT * FROM validated_input_snapshots ORDER BY snapshot_id")
            )
            .mappings()
            .all()
            == snapshots_before
        )
    repository = SqlAlchemyRunRepository(create_session_factory(upgraded))
    assert repository.get_validated_input_binding(snapshot_id) == expected_binding
    assert repository.get_run_data_binding("consumed-run") == expected_binding
    standalone = repository.get_run_data_binding("standalone-run")
    assert standalone.project_id == LEGACY_PROJECT_ID
    assert standalone.provenance.value == "unresolved"
    assert standalone.sample_revisions == ()
    assert (
        standalone.workflow_inputs_digest
        == sha256(STANDALONE_CANONICAL_INPUTS.encode()).hexdigest()
    )

    unconsumed = repository.get_validated_input_snapshot(unconsumed_snapshot_id)
    consumed_at = datetime(2026, 7, 17, 1, 30, tzinfo=timezone.utc)
    post_upgrade_record = RunRecord(
        run_id="post-upgrade-run",
        workflow_id=unconsumed.workflow_id,
        inputs=unconsumed.to_workflow_inputs().to_dict(),
        status=RunStatus.CREATED,
        created_at=consumed_at,
        updated_at=consumed_at,
        started_at=None,
        ended_at=None,
        current_stage=None,
        cancellation_reason=None,
        error=None,
        tags={},
    )
    creation = repository.consume_validated_input_snapshot(
        unconsumed_snapshot_id,
        workflow_id=unconsumed.workflow_id,
        expected_build_identity=unconsumed.workflow_build_identity,
        record=post_upgrade_record,
        consumed_at=consumed_at,
        event=RunEventDraft(
            event_type="status_changed",
            message="Run created.",
            status=RunStatus.CREATED,
            context={"previous_status": None, "new_status": "created"},
        ),
    )
    assert creation.created is True
    assert repository.get_run_data_binding(
        post_upgrade_record.run_id
    ) == repository.get_validated_input_binding(unconsumed_snapshot_id)
    upgraded.dispose()


def test_project_sample_registry_schema_enforces_evidence_relationships(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'project-sample-schema.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    inspector = inspect(engine)

    foreign_keys = {
        table: {
            foreign_key["name"]: foreign_key["options"].get("ondelete")
            for foreign_key in inspector.get_foreign_keys(table)
        }
        for table in (
            "samples",
            "sample_revisions",
            "snapshot_project_bindings",
            "snapshot_sample_revisions",
            "run_project_bindings",
            "run_samples",
        )
    }
    assert foreign_keys == {
        "samples": {"fk_samples_project": "RESTRICT"},
        "sample_revisions": {"fk_sample_revisions_sample": "RESTRICT"},
        "snapshot_project_bindings": {
            "fk_snapshot_project_bindings_project": "RESTRICT",
            "fk_snapshot_project_bindings_snapshot": "RESTRICT",
        },
        "snapshot_sample_revisions": {
            "fk_snapshot_sample_revisions_binding": "RESTRICT",
            "fk_snapshot_sample_revisions_revision": "RESTRICT",
        },
        "run_project_bindings": {
            "fk_run_project_bindings_project": "RESTRICT",
            "fk_run_project_bindings_run": "CASCADE",
        },
        "run_samples": {
            "fk_run_samples_binding": "CASCADE",
            "fk_run_samples_revision": "RESTRICT",
        },
    }
    for table, foreign_key_name in (
        ("snapshot_sample_revisions", "fk_snapshot_sample_revisions_revision"),
        ("run_samples", "fk_run_samples_revision"),
    ):
        revision_foreign_key = next(
            foreign_key
            for foreign_key in inspector.get_foreign_keys(table)
            if foreign_key["name"] == foreign_key_name
        )
        assert revision_foreign_key["constrained_columns"] == [
            "project_id",
            "sample_revision_id",
            "payload_digest",
        ]
        assert revision_foreign_key["referred_columns"] == [
            "project_id",
            "sample_revision_id",
            "payload_digest",
        ]
    checks = {
        table: {
            constraint["name"] for constraint in inspector.get_check_constraints(table)
        }
        for table in (
            "projects",
            "samples",
            "sample_revisions",
            "snapshot_project_bindings",
            "snapshot_sample_revisions",
            "run_project_bindings",
            "run_samples",
        )
    }
    assert checks == {
        "projects": {
            "ck_projects_archive_order",
            "ck_projects_display_name",
            "ck_projects_id",
            "ck_projects_kind",
            "ck_projects_legacy_identity",
        },
        "samples": {
            "ck_samples_id",
            "ck_samples_not_legacy",
            "ck_samples_stable_key",
        },
        "sample_revisions": {
            "ck_sample_revisions_digest_length",
            "ck_sample_revisions_digest_scheme",
            "ck_sample_revisions_id",
            "ck_sample_revisions_positive_revision",
        },
        "snapshot_project_bindings": {
            "ck_snapshot_project_bindings_digest_length",
            "ck_snapshot_project_bindings_digest_scheme",
            "ck_snapshot_project_bindings_legacy_project",
            "ck_snapshot_project_bindings_mode",
            "ck_snapshot_project_bindings_provenance",
            "ck_snapshot_project_bindings_workflow_digest_length",
        },
        "snapshot_sample_revisions": {
            "ck_snapshot_sample_revisions_ordinal",
        },
        "run_project_bindings": {
            "ck_run_project_bindings_digest_length",
            "ck_run_project_bindings_digest_scheme",
            "ck_run_project_bindings_legacy_project",
            "ck_run_project_bindings_mode",
            "ck_run_project_bindings_provenance",
            "ck_run_project_bindings_workflow_digest_length",
        },
        "run_samples": {"ck_run_samples_ordinal"},
    }
    assert inspector.get_pk_constraint("snapshot_sample_revisions")[
        "constrained_columns"
    ] == ["snapshot_id", "ordinal"]
    assert inspector.get_pk_constraint("run_samples")["constrained_columns"] == [
        "run_id",
        "ordinal",
    ]
    assert {
        constraint["name"]
        for constraint in inspector.get_unique_constraints("snapshot_sample_revisions")
    } == {"uq_snapshot_sample_revisions_revision"}
    assert {
        constraint["name"]
        for constraint in inspector.get_unique_constraints("run_samples")
    } == {"uq_run_samples_revision"}
    assert {
        constraint["name"]
        for constraint in inspector.get_unique_constraints("sample_revisions")
    } == {
        "uq_sample_revisions_project_revision",
        "uq_sample_revisions_project_revision_digest",
        "uq_sample_revisions_sample_number",
    }
    with engine.begin() as connection:
        with pytest.raises(IntegrityError, match="ck_samples_not_legacy"):
            connection.execute(
                text(
                    "INSERT INTO samples "
                    "(sample_id, project_id, stable_key, created_at) VALUES "
                    "('smp_00000000000000000000000000000000', "
                    ":legacy_project_id, 'forbidden', :created_at)"
                ),
                {
                    "legacy_project_id": LEGACY_PROJECT_ID,
                    "created_at": "2026-07-26 00:00:00",
                },
            )
    for statement in (
        "UPDATE projects SET display_name='Renamed' "
        f"WHERE project_id='{LEGACY_PROJECT_ID}'",
        f"DELETE FROM projects WHERE project_id='{LEGACY_PROJECT_ID}'",
    ):
        with engine.begin() as connection:
            with pytest.raises(
                IntegrityError,
                match="reserved Legacy Project is immutable",
            ):
                connection.execute(text(statement))
    assert {index["name"] for index in inspector.get_indexes("sample_revisions")} >= {
        "ix_sample_revisions_project_created",
        "ix_sample_revisions_sample_revision",
    }
    assert {
        index["name"] for index in inspector.get_indexes("run_project_bindings")
    } >= {"ix_run_project_bindings_project"}
    assert {index["name"] for index in inspector.get_indexes("run_samples")} >= {
        "ix_run_samples_revision"
    }
    assert {
        index["name"] for index in inspector.get_indexes("snapshot_project_bindings")
    } >= {"ix_snapshot_project_bindings_project"}
    assert {
        index["name"] for index in inspector.get_indexes("snapshot_sample_revisions")
    } >= {"ix_snapshot_sample_revisions_revision"}

    with engine.begin() as connection:
        connection.execute(
            text(
                "INSERT INTO projects "
                "(project_id, display_name, kind, created_at, archived_at) VALUES "
                "('prj_11111111111111111111111111111111', 'A', 'user', "
                ":created_at, NULL), "
                "('prj_22222222222222222222222222222222', 'B', 'user', "
                ":created_at, NULL)"
            ),
            {"created_at": "2026-07-26 00:00:00"},
        )
        connection.execute(
            text(
                "INSERT INTO samples "
                "(sample_id, project_id, stable_key, created_at) "
                "VALUES ('smp_11111111111111111111111111111111', "
                "'prj_11111111111111111111111111111111', 'sample-a', "
                ":created_at)"
            ),
            {"created_at": "2026-07-26 00:00:00"},
        )
        connection.execute(
            text(
                "INSERT INTO sample_revisions "
                "(sample_revision_id, project_id, sample_id, revision_number, "
                "canonical_payload, payload_digest_scheme, payload_digest, created_at) "
                "VALUES ('smpr_11111111111111111111111111111111', "
                "'prj_11111111111111111111111111111111', "
                "'smp_11111111111111111111111111111111', 1, '{}', "
                "'sha256-framed-sample-revision-payload-v1', "
                ":digest, :created_at)"
            ),
            {"digest": "1" * 64, "created_at": "2026-07-26 00:00:00"},
        )
        for statement in (
            "UPDATE samples SET stable_key='renamed' "
            "WHERE sample_id='smp_11111111111111111111111111111111'",
            "DELETE FROM samples "
            "WHERE sample_id='smp_11111111111111111111111111111111'",
            "UPDATE sample_revisions SET canonical_payload='tampered' "
            "WHERE sample_revision_id="
            "'smpr_11111111111111111111111111111111'",
            "DELETE FROM sample_revisions WHERE sample_revision_id="
            "'smpr_11111111111111111111111111111111'",
        ):
            with pytest.raises(IntegrityError, match="immutable"):
                connection.execute(text(statement))
        connection.execute(
            text(
                "INSERT INTO runs "
                "(run_id, workflow_id, inputs, status, created_at, updated_at, "
                "started_at, ended_at, current_stage, cancellation_reason, "
                "error, tags) VALUES ('run-cross-project', 'workflow', '{}', "
                "'created', :created_at, :created_at, NULL, NULL, NULL, NULL, "
                "NULL, '{}')"
            ),
            {"created_at": "2026-07-26 00:00:00"},
        )
        connection.execute(
            text(
                "INSERT INTO run_project_bindings "
                "(run_id, project_id, binding_mode, provenance, "
                "workflow_inputs_digest, binding_digest_scheme, "
                "binding_digest, created_at) VALUES "
                "('run-cross-project', "
                "'prj_22222222222222222222222222222222', 'bound_v1', "
                "'resolved', :workflow_digest, "
                "'sha256-framed-project-sample-binding-v1', "
                ":digest, :created_at)"
            ),
            {
                "workflow_digest": "3" * 64,
                "digest": "2" * 64,
                "created_at": "2026-07-26 00:00:00",
            },
        )
        with pytest.raises(IntegrityError, match="FOREIGN KEY constraint failed"):
            connection.execute(
                text(
                    "INSERT INTO run_samples "
                    "(run_id, project_id, sample_revision_id, payload_digest, "
                    "ordinal) VALUES "
                    "('run-cross-project', "
                    "'prj_22222222222222222222222222222222', "
                    "'smpr_11111111111111111111111111111111', :payload_digest, 0)"
                ),
                {"payload_digest": "1" * 64},
            )
    engine.dispose()


def test_project_sample_registry_downgrade_and_reapply_is_deterministic(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'project-sample-reapply.db'}"
    upgrade_database(database_url, "20260717_08")
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        _insert_rev08_run(connection, "legacy-run", created_at="2026-07-17 01:00:00")
    engine.dispose()

    upgrade_database(database_url)
    first = create_database_engine(database_url)
    with first.connect() as connection:
        first_digest = connection.scalar(
            text(
                "SELECT binding_digest FROM run_project_bindings "
                "WHERE run_id='legacy-run'"
            )
        )
    first.dispose()

    downgrade_database(database_url, "20260717_08")
    downgraded = create_database_engine(database_url)
    assert {
        "projects",
        "samples",
        "sample_revisions",
        "snapshot_project_bindings",
        "snapshot_sample_revisions",
        "run_project_bindings",
        "run_samples",
    }.isdisjoint(inspect(downgraded).get_table_names())
    with downgraded.connect() as connection:
        assert connection.scalar(text("SELECT count(*) FROM runs")) == 1
    downgraded.dispose()

    upgrade_database(database_url)
    reapplied = create_database_engine(database_url)
    with reapplied.connect() as connection:
        assert (
            connection.scalar(
                text(
                    "SELECT binding_digest FROM run_project_bindings "
                    "WHERE run_id='legacy-run'"
                )
            )
            == first_digest
        )
    reapplied.dispose()


def test_project_sample_registry_upgrade_rejects_inconsistent_consumed_snapshot(
    tmp_path,
) -> None:
    database_path = tmp_path / "inconsistent-consumed.db"
    database_url = f"sqlite:///{database_path}"
    upgrade_database(database_url, "20260717_08")
    engine = create_database_engine(database_url)
    engine.dispose()

    import sqlite3

    connection = sqlite3.connect(database_path)
    connection.execute("PRAGMA foreign_keys=OFF")
    _insert_rev08_snapshot(
        connection,
        "vsnap_" + "f" * 32,
        digest=REV08_INPUTS_DIGEST,
        consumed_run_id="missing-run",
    )
    connection.commit()
    connection.close()

    with pytest.raises(RuntimeError, match="inconsistent consumed snapshot linkage"):
        upgrade_database(database_url)


@pytest.mark.parametrize("mismatch", ["workflow", "inputs", "payload-digest"])
def test_project_sample_registry_upgrade_rejects_mismatched_consumed_evidence(
    tmp_path,
    mismatch,
) -> None:
    database_url = f"sqlite:///{tmp_path / f'mismatched-{mismatch}.db'}"
    upgrade_database(database_url, "20260717_08")
    engine = create_database_engine(database_url)
    run_inputs = (
        json.dumps(
            {"config": {"threads": 2}, "options": {}, "samples": None},
            sort_keys=True,
            separators=(",", ":"),
        )
        if mismatch == "inputs"
        else REV08_CANONICAL_INPUTS
    )
    snapshot_workflow_id = "other-workflow" if mismatch == "workflow" else "workflow"
    snapshot_digest = "e" * 64 if mismatch == "payload-digest" else REV08_INPUTS_DIGEST
    with engine.begin() as connection:
        _insert_rev08_run(
            connection,
            "consumed-run",
            created_at="2026-07-17 01:00:00",
            inputs=run_inputs,
        )
        _insert_rev08_snapshot(
            connection,
            "vsnap_" + "e" * 32,
            digest=snapshot_digest,
            consumed_run_id="consumed-run",
            workflow_id=snapshot_workflow_id,
        )
    engine.dispose()

    with pytest.raises(
        RuntimeError,
        match=r"inconsistent (consumed )?snapshot evidence",
    ):
        upgrade_database(database_url)


def test_project_sample_registry_upgrade_preflights_before_ddl_and_can_retry(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'preflight-retry.db'}"
    upgrade_database(database_url, "20260717_08")
    engine = create_database_engine(database_url)
    invalid_inputs = '{"config":{"bad":NaN},"options":{},"samples":null}'
    with engine.begin() as connection:
        _insert_rev08_run(
            connection,
            "standalone-invalid",
            created_at="2026-07-17 01:00:00",
            inputs=invalid_inputs,
        )
    engine.dispose()

    with pytest.raises(
        RuntimeError,
        match="invalid legacy run inputs block registry upgrade",
    ):
        upgrade_database(database_url)

    failed = create_database_engine(database_url)
    rev09_tables = {
        "projects",
        "samples",
        "sample_revisions",
        "snapshot_project_bindings",
        "snapshot_sample_revisions",
        "run_project_bindings",
        "run_samples",
    }
    assert rev09_tables.isdisjoint(inspect(failed).get_table_names())
    with failed.begin() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260717_08"
        )
        connection.execute(
            text("UPDATE runs SET inputs=:inputs WHERE run_id='standalone-invalid'"),
            {"inputs": STANDALONE_CANONICAL_INPUTS},
        )
    failed.dispose()

    upgrade_database(database_url)
    retried = create_database_engine(database_url)
    with retried.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260827_15"
        )
        assert (
            connection.scalar(
                text(
                    "SELECT count(*) FROM run_project_bindings "
                    "WHERE run_id='standalone-invalid'"
                )
            )
            == 1
        )
    retried.dispose()


@pytest.mark.parametrize("starting_revision", ["20260717_08", "20260726_09"])
def test_input_registry_upgrade_backfills_only_unresolved_compatibility_envelopes(
    tmp_path,
    starting_revision,
) -> None:
    database_url = f"sqlite:///{tmp_path / f'input-registry-{starting_revision}.db'}"
    upgrade_database(database_url, "20260717_08")
    snapshot_id = "vsnap_" + "c" * 32
    workflow_id = "workflow+" + "v" * 247
    assert len(workflow_id) == 256
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        _insert_rev08_run(
            connection,
            "consumed-input-run",
            created_at="2026-07-17 01:00:00",
            workflow_id=workflow_id,
        )
        _insert_rev08_run(
            connection,
            "standalone-input-run",
            created_at="2026-07-17 02:00:00",
            workflow_id=workflow_id,
            inputs=STANDALONE_NONCANONICAL_INPUTS,
        )
        _insert_rev08_snapshot(
            connection,
            snapshot_id,
            digest=REV08_INPUTS_DIGEST,
            consumed_run_id="consumed-input-run",
            workflow_id=workflow_id,
        )
    engine.dispose()
    if starting_revision == "20260726_09":
        upgrade_database(database_url, starting_revision)

    upgrade_database(database_url)
    upgraded = create_database_engine(database_url)
    with upgraded.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260827_15"
        )
        snapshot_binding = (
            connection.execute(
                text(
                    "SELECT project_id, workflow_id, adapter_contract_version, "
                    "binding_mode, workflow_inputs_digest, "
                    "project_sample_binding_digest, binding_digest_scheme, "
                    "binding_digest "
                    "FROM snapshot_input_bindings WHERE snapshot_id=:snapshot_id"
                ),
                {"snapshot_id": snapshot_id},
            )
            .mappings()
            .one()
        )
        run_bindings = {
            row["run_id"]: row
            for row in connection.execute(
                text(
                    "SELECT run_id, project_id, workflow_id, "
                    "adapter_contract_version, binding_mode, "
                    "workflow_inputs_digest, project_sample_binding_digest, "
                    "binding_digest_scheme, binding_digest "
                    "FROM run_input_bindings ORDER BY run_id"
                )
            )
            .mappings()
            .all()
        }
        consumed_binding = run_bindings["consumed-input-run"]
        assert dict(consumed_binding) == {
            "run_id": "consumed-input-run",
            **dict(snapshot_binding),
        }
        project_sample_binding_digest = connection.scalar(
            text(
                "SELECT binding_digest FROM snapshot_project_bindings "
                "WHERE snapshot_id=:snapshot_id"
            ),
            {"snapshot_id": snapshot_id},
        )
        assert (
            connection.scalar(
                text(
                    "SELECT adapter_version FROM validated_input_snapshots "
                    "WHERE snapshot_id=:snapshot_id"
                ),
                {"snapshot_id": snapshot_id},
            )
            == "1.0.0"
        )
        assert dict(snapshot_binding) == {
            "project_id": LEGACY_PROJECT_ID,
            "workflow_id": workflow_id,
            "adapter_contract_version": None,
            "binding_mode": "compatibility_unresolved_v1",
            "workflow_inputs_digest": REV08_INPUTS_DIGEST,
            "project_sample_binding_digest": project_sample_binding_digest,
            "binding_digest_scheme": INPUT_BINDING_DIGEST_SCHEME,
            "binding_digest": _compatibility_input_binding_digest(
                project_id=LEGACY_PROJECT_ID,
                workflow_id=workflow_id,
                adapter_contract_version=None,
                workflow_inputs_digest=REV08_INPUTS_DIGEST,
                project_sample_binding_digest=project_sample_binding_digest,
            ),
        }
        snapshot_envelope = InputUseBindingEnvelope(
            project_id=snapshot_binding["project_id"],
            project_sample_binding_digest=snapshot_binding[
                "project_sample_binding_digest"
            ],
            workflow_id=snapshot_binding["workflow_id"],
            adapter_contract_version=snapshot_binding["adapter_contract_version"],
            workflow_inputs_digest=snapshot_binding["workflow_inputs_digest"],
            contract_mode=snapshot_binding["binding_mode"],
            input_uses=(),
            digest_scheme=snapshot_binding["binding_digest_scheme"],
            digest=snapshot_binding["binding_digest"],
        )
        assert snapshot_envelope.adapter_contract_version is None
        assert snapshot_envelope.input_uses == ()
        standalone = run_bindings["standalone-input-run"]
        assert standalone["adapter_contract_version"] is None
        assert standalone["binding_mode"] == "compatibility_unresolved_v1"
        assert standalone["binding_digest"] == _compatibility_input_binding_digest(
            project_id=LEGACY_PROJECT_ID,
            workflow_id=workflow_id,
            adapter_contract_version=None,
            workflow_inputs_digest=standalone["workflow_inputs_digest"],
            project_sample_binding_digest=standalone["project_sample_binding_digest"],
        )
        assert (
            InputUseBindingEnvelope(
                project_id=standalone["project_id"],
                project_sample_binding_digest=standalone[
                    "project_sample_binding_digest"
                ],
                workflow_id=standalone["workflow_id"],
                adapter_contract_version=standalone["adapter_contract_version"],
                workflow_inputs_digest=standalone["workflow_inputs_digest"],
                contract_mode=standalone["binding_mode"],
                input_uses=(),
                digest_scheme=standalone["binding_digest_scheme"],
                digest=standalone["binding_digest"],
            ).adapter_contract_version
            is None
        )
        assert connection.scalar(text("SELECT count(*) FROM storage_pools")) == 0
        assert (
            connection.scalar(
                text("SELECT count(*) FROM project_storage_pool_bindings")
            )
            == 0
        )
        assert connection.scalar(text("SELECT count(*) FROM input_files")) == 0
        assert connection.scalar(text("SELECT count(*) FROM input_file_revisions")) == 0
        assert connection.scalar(text("SELECT count(*) FROM snapshot_input_uses")) == 0
        assert (
            connection.scalar(text("SELECT count(*) FROM snapshot_input_members")) == 0
        )
        assert connection.scalar(text("SELECT count(*) FROM run_input_uses")) == 0
        assert connection.scalar(text("SELECT count(*) FROM run_input_members")) == 0
    upgraded.dispose()


def test_input_registry_upgrade_preflights_stage2_evidence_before_ddl_and_can_retry(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'input-registry-preflight.db'}"
    upgrade_database(database_url, "20260717_08")
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        _insert_rev08_run(
            connection,
            "preflight-run",
            created_at="2026-07-17 01:00:00",
        )
    engine.dispose()
    upgrade_database(database_url, "20260726_09")

    rev09 = create_database_engine(database_url)
    with rev09.begin() as connection:
        valid_digest = connection.scalar(
            text(
                "SELECT binding_digest FROM run_project_bindings "
                "WHERE run_id='preflight-run'"
            )
        )
        connection.execute(
            text(
                "UPDATE run_project_bindings SET binding_digest=:digest "
                "WHERE run_id='preflight-run'"
            ),
            {"digest": "0" * 64},
        )
    rev09.dispose()

    with pytest.raises(
        RuntimeError,
        match="inconsistent project/sample binding evidence blocks input registry upgrade",
    ):
        upgrade_database(database_url)

    failed = create_database_engine(database_url)
    assert INPUT_REGISTRY_TABLES.isdisjoint(inspect(failed).get_table_names())
    with failed.begin() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260726_09"
        )
        connection.execute(
            text(
                "UPDATE run_project_bindings SET binding_digest=:digest "
                "WHERE run_id='preflight-run'"
            ),
            {"digest": valid_digest},
        )
    failed.dispose()

    upgrade_database(database_url)
    retried = create_database_engine(database_url)
    assert INPUT_REGISTRY_TABLES <= set(inspect(retried).get_table_names())
    with retried.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260827_15"
        )
    retried.dispose()


def test_input_registry_schema_is_path_private_append_only_and_mode_scoped(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'input-registry-schema.db'}"
    upgrade_database(database_url, "20260717_08")
    snapshot_id = "vsnap_" + "d" * 32
    rev08 = create_database_engine(database_url)
    with rev08.begin() as connection:
        _insert_rev08_run(
            connection,
            "compatibility-input-run",
            created_at="2026-07-17 01:00:00",
        )
        _insert_rev08_snapshot(
            connection,
            snapshot_id,
            digest=REV08_INPUTS_DIGEST,
            consumed_run_id="compatibility-input-run",
        )
    rev08.dispose()
    upgrade_database(database_url)

    engine = create_database_engine(database_url)
    inspector = inspect(engine)
    assert {column["name"] for column in inspector.get_columns("storage_pools")} == {
        "storage_pool_id",
        "config_key",
        "display_name",
        "created_at",
        "archived_at",
    }
    assert {
        column["name"]
        for column in inspector.get_columns("project_storage_pool_bindings")
    } == {"project_id", "storage_pool_id", "bound_at"}
    assert all(
        "root" not in column["name"] and "absolute" not in column["name"]
        for table_name in INPUT_REGISTRY_TABLES
        for column in inspector.get_columns(table_name)
    )

    project_id = "prj_" + "1" * 32
    pool_id = "stgp_" + "2" * 32
    input_file_id = "inpf_" + "3" * 32
    revision_id = "inpfr_" + "4" * 32
    created_at = "2026-07-26 01:00:00"
    with engine.begin() as connection:
        connection.execute(
            text(
                "INSERT INTO projects "
                "(project_id, display_name, kind, created_at, archived_at) "
                "VALUES (:project_id, 'Input Project', 'user', :created_at, NULL)"
            ),
            {"project_id": project_id, "created_at": created_at},
        )
        connection.execute(
            text(
                "INSERT INTO storage_pools "
                "(storage_pool_id, config_key, display_name, created_at, "
                "archived_at) VALUES "
                "(:pool_id, 'approved-primary', 'Approved primary', "
                ":created_at, NULL)"
            ),
            {"pool_id": pool_id, "created_at": created_at},
        )
        connection.execute(
            text(
                "INSERT INTO project_storage_pool_bindings "
                "(project_id, storage_pool_id, bound_at) "
                "VALUES (:project_id, :pool_id, :created_at)"
            ),
            {
                "project_id": project_id,
                "pool_id": pool_id,
                "created_at": created_at,
            },
        )
        connection.execute(
            text(
                "INSERT INTO input_files "
                "(input_file_id, project_id, storage_pool_id, stable_key, "
                "created_at, archived_at) VALUES "
                "(:input_file_id, :project_id, :pool_id, 'reads-r1', "
                ":created_at, NULL)"
            ),
            {
                "input_file_id": input_file_id,
                "project_id": project_id,
                "pool_id": pool_id,
                "created_at": created_at,
            },
        )
        connection.execute(
            text(
                "INSERT INTO input_file_revisions "
                "(input_file_revision_id, input_file_id, project_id, "
                "storage_pool_id, revision_number, relative_path, size_bytes, "
                "content_sha256, digest_scheme, digest, created_at) VALUES "
                "(:revision_id, :input_file_id, :project_id, :pool_id, 1, "
                "'incoming/reads.fastq.gz', 5, :content_sha256, "
                "'sha256-framed-input-file-revision-v1', :digest, :created_at)"
            ),
            {
                "revision_id": revision_id,
                "input_file_id": input_file_id,
                "project_id": project_id,
                "pool_id": pool_id,
                "content_sha256": "a" * 64,
                "digest": "b" * 64,
                "created_at": created_at,
            },
        )
        connection.execute(
            text(
                "INSERT INTO runs "
                "(run_id, workflow_id, inputs, status, created_at, updated_at, "
                "started_at, ended_at, current_stage, cancellation_reason, "
                "error, tags) VALUES "
                "('declared-input-run', 'workflow', '{}', 'created', "
                ":created_at, :created_at, NULL, NULL, NULL, NULL, NULL, '{}')"
            ),
            {"created_at": created_at},
        )
        connection.execute(
            text(
                "INSERT INTO run_result_states "
                "(run_id, artifact_revision, qc_revision) "
                "VALUES ('declared-input-run', 0, 0)"
            )
        )
        connection.execute(
            text(
                "INSERT INTO run_project_bindings "
                "(run_id, project_id, binding_mode, provenance, "
                "workflow_inputs_digest, binding_digest_scheme, "
                "binding_digest, created_at) VALUES "
                "('declared-input-run', :project_id, 'bound_v1', 'resolved', "
                ":workflow_digest, 'sha256-framed-project-sample-binding-v1', "
                ":project_sample_digest, :created_at)"
            ),
            {
                "project_id": project_id,
                "workflow_digest": "c" * 64,
                "project_sample_digest": "d" * 64,
                "created_at": created_at,
            },
        )
        connection.execute(
            text(
                "INSERT INTO run_input_bindings "
                "(run_id, project_id, workflow_id, adapter_contract_version, "
                "binding_mode, workflow_inputs_digest, "
                "project_sample_binding_digest, binding_digest_scheme, "
                "binding_digest, created_at) VALUES "
                "('declared-input-run', :project_id, 'workflow', '2.0.0', "
                "'declared_input_uses_v1', :workflow_digest, "
                ":project_sample_digest, "
                "'sha256-framed-input-use-binding-envelope-v1', "
                ":binding_digest, :created_at)"
            ),
            {
                "project_id": project_id,
                "workflow_digest": "c" * 64,
                "project_sample_digest": "d" * 64,
                "binding_digest": "e" * 64,
                "created_at": created_at,
            },
        )
        connection.execute(
            text(
                "INSERT INTO run_input_uses "
                "(run_id, project_id, ordinal, input_use_key, occurrence, "
                "capability_version, closure_contract_version, provenance_mode, "
                "closure_digest_scheme, closure_digest) VALUES "
                "('declared-input-run', :project_id, 0, 'transitional', 0, "
                "'capability-v1', 'regular_file_v1', "
                "'transitional_unmanaged_v1', NULL, NULL)"
            ),
            {"project_id": project_id},
        )
        connection.execute(
            text(
                "INSERT INTO run_input_uses "
                "(run_id, project_id, ordinal, input_use_key, occurrence, "
                "capability_version, closure_contract_version, provenance_mode, "
                "closure_digest_scheme, closure_digest) VALUES "
                "('declared-input-run', :project_id, 1, 'managed', 0, "
                "'capability-v1', 'regular_file_v1', 'managed_revision_v1', "
                "'sha256-framed-input-closure-v1', :closure_digest)"
            ),
            {"project_id": project_id, "closure_digest": "f" * 64},
        )

    with engine.begin() as connection:
        with pytest.raises(
            IntegrityError,
            match="ck_input_file_revisions_safe_relative_path",
        ):
            connection.execute(
                text(
                    "INSERT INTO input_file_revisions "
                    "(input_file_revision_id, input_file_id, project_id, "
                    "storage_pool_id, revision_number, relative_path, "
                    "size_bytes, content_sha256, digest_scheme, digest, "
                    "created_at) VALUES "
                    "(:revision_id, :input_file_id, :project_id, :pool_id, 2, "
                    "'/srv/private/reads.fastq.gz', 5, :content_sha256, "
                    "'sha256-framed-input-file-revision-v1', :digest, "
                    ":created_at)"
                ),
                {
                    "revision_id": "inpfr_" + "5" * 32,
                    "input_file_id": input_file_id,
                    "project_id": project_id,
                    "pool_id": pool_id,
                    "content_sha256": "a" * 64,
                    "digest": "f" * 64,
                    "created_at": created_at,
                },
            )

    member_parameters = {
        "project_id": project_id,
        "input_file_id": input_file_id,
        "revision_id": revision_id,
        "revision_digest": "b" * 64,
        "content_sha256": "a" * 64,
    }
    with engine.begin() as connection:
        with pytest.raises(
            IntegrityError,
            match="input uses require declared_input_uses_v1 binding",
        ):
            connection.execute(
                text(
                    "INSERT INTO run_input_uses "
                    "(run_id, project_id, ordinal, input_use_key, occurrence, "
                    "capability_version, closure_contract_version, "
                    "provenance_mode, closure_digest_scheme, closure_digest) "
                    "VALUES ('compatibility-input-run', :project_id, 0, "
                    "'synthetic', 0, 'capability-v1', 'regular_file_v1', "
                    "'transitional_unmanaged_v1', NULL, NULL)"
                ),
                {"project_id": LEGACY_PROJECT_ID},
            )
    with engine.begin() as connection:
        with pytest.raises(
            IntegrityError,
            match="input members require managed_revision_v1 provenance",
        ):
            connection.execute(
                text(
                    "INSERT INTO run_input_members "
                    "(run_id, project_id, use_ordinal, member_ordinal, "
                    "logical_member_key, input_file_id, input_file_revision_id, "
                    "revision_digest, size_bytes, content_sha256) VALUES "
                    "('declared-input-run', :project_id, 0, 0, 'reads', "
                    ":input_file_id, :revision_id, :revision_digest, 5, "
                    ":content_sha256)"
                ),
                member_parameters,
            )
    with engine.begin() as connection:
        with pytest.raises(IntegrityError, match="FOREIGN KEY constraint failed"):
            connection.execute(
                text(
                    "INSERT INTO run_input_members "
                    "(run_id, project_id, use_ordinal, member_ordinal, "
                    "logical_member_key, input_file_id, input_file_revision_id, "
                    "revision_digest, size_bytes, content_sha256) VALUES "
                    "('declared-input-run', :project_id, 1, 0, 'reads', "
                    ":input_file_id, :revision_id, :revision_digest, 6, "
                    ":content_sha256)"
                ),
                member_parameters,
            )
    with engine.begin() as connection:
        connection.execute(
            text(
                "INSERT INTO run_input_members "
                "(run_id, project_id, use_ordinal, member_ordinal, "
                "logical_member_key, input_file_id, input_file_revision_id, "
                "revision_digest, size_bytes, content_sha256) VALUES "
                "('declared-input-run', :project_id, 1, 0, 'reads', "
                ":input_file_id, :revision_id, :revision_digest, 5, "
                ":content_sha256)"
            ),
            member_parameters,
        )

    immutable_statements = (
        (
            "UPDATE storage_pools SET display_name='Renamed' "
            "WHERE storage_pool_id=:pool_id",
            "StoragePool permits only one archive transition",
        ),
        (
            "DELETE FROM project_storage_pool_bindings WHERE project_id=:project_id",
            "Project/StoragePool binding is immutable",
        ),
        (
            "UPDATE input_files SET stable_key='renamed' "
            "WHERE input_file_id=:input_file_id",
            "InputFile permits only one archive transition",
        ),
        (
            "UPDATE input_file_revisions SET size_bytes=6 "
            "WHERE input_file_revision_id=:revision_id",
            "InputFileRevision is immutable",
        ),
        (
            "UPDATE run_input_bindings SET binding_digest=:replacement_digest "
            "WHERE run_id='compatibility-input-run'",
            "Run input binding is immutable",
        ),
        (
            "DELETE FROM run_input_bindings WHERE run_id='compatibility-input-run'",
            "Run input binding is immutable",
        ),
        (
            "DELETE FROM snapshot_input_bindings WHERE snapshot_id=:snapshot_id",
            "Snapshot input binding is immutable",
        ),
    )
    statement_parameters = {
        "pool_id": pool_id,
        "project_id": project_id,
        "input_file_id": input_file_id,
        "revision_id": revision_id,
        "replacement_digest": "9" * 64,
        "snapshot_id": snapshot_id,
    }
    for statement, message in immutable_statements:
        with engine.begin() as connection:
            with pytest.raises(IntegrityError, match=message):
                connection.execute(text(statement), statement_parameters)

    with engine.begin() as connection:
        connection.execute(
            text(
                "UPDATE storage_pools SET archived_at='2026-07-26 02:00:00' "
                "WHERE storage_pool_id=:pool_id"
            ),
            {"pool_id": pool_id},
        )
    with engine.begin() as connection:
        with pytest.raises(
            IntegrityError,
            match="StoragePool permits only one archive transition",
        ):
            connection.execute(
                text(
                    "UPDATE storage_pools SET archived_at=NULL "
                    "WHERE storage_pool_id=:pool_id"
                ),
                {"pool_id": pool_id},
            )
    engine.dispose()


def test_input_registry_downgrade_and_reapply_is_deterministic(tmp_path) -> None:
    database_url = f"sqlite:///{tmp_path / 'input-registry-reapply.db'}"
    upgrade_database(database_url, "20260717_08")
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        _insert_rev08_run(
            connection,
            "input-reapply-run",
            created_at="2026-07-17 01:00:00",
        )
    engine.dispose()

    upgrade_database(database_url)
    first = create_database_engine(database_url)
    with first.connect() as connection:
        first_binding = dict(
            connection.execute(
                text(
                    "SELECT project_id, workflow_id, adapter_contract_version, "
                    "binding_mode, workflow_inputs_digest, "
                    "project_sample_binding_digest, binding_digest_scheme, "
                    "binding_digest FROM run_input_bindings "
                    "WHERE run_id='input-reapply-run'"
                )
            )
            .mappings()
            .one()
        )
    first.dispose()

    downgrade_database(database_url, "20260726_09")
    downgraded = create_database_engine(database_url)
    assert INPUT_REGISTRY_TABLES.isdisjoint(inspect(downgraded).get_table_names())
    with downgraded.connect() as connection:
        assert connection.scalar(text("SELECT count(*) FROM runs")) == 1
        assert connection.scalar(text("SELECT count(*) FROM run_project_bindings")) == 1
    downgraded.dispose()

    upgrade_database(database_url)
    reapplied = create_database_engine(database_url)
    with reapplied.connect() as connection:
        assert (
            dict(
                connection.execute(
                    text(
                        "SELECT project_id, workflow_id, "
                        "adapter_contract_version, binding_mode, "
                        "workflow_inputs_digest, project_sample_binding_digest, "
                        "binding_digest_scheme, binding_digest "
                        "FROM run_input_bindings "
                        "WHERE run_id='input-reapply-run'"
                    )
                )
                .mappings()
                .one()
            )
            == first_binding
        )
    reapplied.dispose()


def _compatibility_input_binding_digest(
    *,
    project_id: str,
    workflow_id: str,
    adapter_contract_version: str | None,
    workflow_inputs_digest: str,
    project_sample_binding_digest: str,
) -> str:
    canonical = json.dumps(
        {
            "adapter_contract_version": adapter_contract_version,
            "contract_mode": "compatibility_unresolved_v1",
            "input_uses": [],
            "project_id": project_id,
            "project_sample_binding_digest": project_sample_binding_digest,
            "workflow_id": workflow_id,
            "workflow_inputs_digest": workflow_inputs_digest,
        },
        ensure_ascii=True,
        allow_nan=False,
        separators=(",", ":"),
        sort_keys=True,
    )
    digest = sha256()
    for part in (INPUT_BINDING_DIGEST_SCHEME, canonical):
        encoded = part.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, byteorder="big", signed=False))
        digest.update(encoded)
    return digest.hexdigest()


def _insert_rev08_run(
    connection,
    run_id: str,
    *,
    created_at: str,
    workflow_id: str = "workflow",
    inputs: str = REV08_CANONICAL_INPUTS,
) -> None:
    connection.execute(
        text(
            "INSERT INTO runs "
            "(run_id, workflow_id, inputs, status, created_at, updated_at, "
            "started_at, ended_at, current_stage, cancellation_reason, error, tags) "
            "VALUES (:run_id, :workflow_id, :inputs, 'created', :created_at, "
            ":created_at, NULL, NULL, NULL, NULL, NULL, '{}')"
        ),
        {
            "run_id": run_id,
            "workflow_id": workflow_id,
            "inputs": inputs,
            "created_at": created_at,
        },
    )
    connection.execute(
        text(
            "INSERT INTO run_result_states "
            "(run_id, artifact_revision, qc_revision) VALUES (:run_id, 0, 0)"
        ),
        {"run_id": run_id},
    )


def _insert_rev08_snapshot(
    connection,
    snapshot_id: str,
    *,
    digest: str,
    consumed_run_id: str | None,
    workflow_id: str = "workflow",
    canonical_payload: str = REV08_CANONICAL_INPUTS,
) -> None:
    values = {
        "snapshot_id": snapshot_id,
        "workflow_id": workflow_id,
        "canonical_payload": canonical_payload,
        "digest": digest,
        "consumed_run_id": consumed_run_id,
        "consumed_at": ("2026-07-17 01:30:00" if consumed_run_id is not None else None),
    }
    statement = text(
        "INSERT INTO validated_input_snapshots "
        "(snapshot_id, workflow_id, adapter_version, schema_version, "
        "schema_dialect, canonical_payload, payload_digest_scheme, payload_digest, "
        "validation_outcome, validation_issue_codes, validated_at, expires_at, "
        "build_adapter_version, build_scheme, build_logical_entrypoint, "
        "build_digest, build_captured_at, consumed_run_id, consumed_at) VALUES "
        "(:snapshot_id, :workflow_id, '1.0.0', '1.0.0', "
        "'https://json-schema.org/draft/2020-12/schema', :canonical_payload, "
        "'sha256-framed-workflow-inputs-v1', :digest, "
        "'adapter_validation_succeeded', '[]', '2026-07-17 01:00:00', "
        "'2026-07-18 01:00:00', '1.0.0', 'sha256-tree-v1', "
        "'workflow/Snakefile', :digest, '2026-07-17 01:00:00', "
        ":consumed_run_id, :consumed_at)"
    )
    if connection.__class__.__module__ == "sqlite3":
        connection.execute(str(statement), values)
    else:
        connection.execute(statement, values)


def test_pr123_assignment_upgrades_with_empty_cancellation_evidence(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url, "20260711_02")
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        connection.execute(
            text(
                """
                INSERT INTO runs (
                    run_id, workflow_id, inputs, status, created_at, updated_at,
                    started_at, ended_at, current_stage, cancellation_reason,
                    error, tags
                ) VALUES (
                    'legacy-running', 'workflow', '{}', 'running',
                    :timestamp, :timestamp, :timestamp, NULL, 'execution',
                    NULL, NULL, '{}'
                )
                """
            ),
            {"timestamp": "2026-07-12 00:00:00"},
        )
        connection.execute(
            text(
                """
                INSERT INTO run_execution_assignments (
                    run_id, job_id, backend, queue_name, created_at,
                    dispatched_at, claimed_at
                ) VALUES (
                    'legacy-running', 'legacy-job', 'rq', 'default',
                    :timestamp, :timestamp, :timestamp
                )
                """
            ),
            {"timestamp": "2026-07-12 00:00:00"},
        )
    engine.dispose()

    upgrade_database(database_url)
    upgraded = create_database_engine(database_url)
    with upgraded.connect() as connection:
        row = (
            connection.execute(
                text(
                    """
                SELECT job_id, dispatched_at, claimed_at,
                       cancellation_requested_at, cancellation_reason,
                       cancellation_acknowledged_at
                FROM run_execution_assignments
                WHERE run_id = 'legacy-running'
                """
                )
            )
            .mappings()
            .one()
        )
        assert row["job_id"] == "legacy-job"
        assert row["dispatched_at"] is not None
        assert row["claimed_at"] is not None
        assert row["cancellation_requested_at"] is None
        assert row["cancellation_reason"] is None
        assert row["cancellation_acknowledged_at"] is None
    upgraded.dispose()


def test_cancellation_intent_migration_downgrades_without_losing_assignment(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    current = create_database_engine(database_url)
    with current.begin() as connection:
        connection.execute(
            text(
                """
                INSERT INTO runs (
                    run_id, workflow_id, inputs, status, created_at, updated_at,
                    started_at, ended_at, current_stage, cancellation_reason,
                    error, tags
                ) VALUES (
                    'downgrade-running', 'workflow', '{}', 'running',
                    :timestamp, :timestamp, :timestamp, NULL, 'execution',
                    NULL, NULL, '{}'
                )
                """
            ),
            {"timestamp": "2026-07-12 00:00:00"},
        )
        connection.execute(
            text(
                """
                INSERT INTO run_execution_assignments (
                    run_id, job_id, backend, queue_name, created_at,
                    dispatched_at, claimed_at, cancellation_requested_at,
                    cancellation_reason, cancellation_acknowledged_at
                ) VALUES (
                    'downgrade-running', 'downgrade-job', 'rq', 'default',
                    :timestamp, :timestamp, :timestamp, :timestamp,
                    'User requested cancellation.', :timestamp
                )
                """
            ),
            {"timestamp": "2026-07-12 00:00:00"},
        )
    current.dispose()

    downgrade_database(database_url, "20260712_03")
    engine = create_database_engine(database_url)
    inspector = inspect(engine)
    columns = {
        column["name"] for column in inspector.get_columns("run_execution_assignments")
    }
    assert "cancellation_requested_at" not in columns
    assert "cancellation_reason" not in columns
    assert "cancellation_acknowledged_at" not in columns
    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260712_03"
        )
        assignment = (
            connection.execute(
                text(
                    """
                    SELECT job_id, dispatched_at, claimed_at
                    FROM run_execution_assignments
                    WHERE run_id = 'downgrade-running'
                    """
                )
            )
            .mappings()
            .one()
        )
        assert assignment["job_id"] == "downgrade-job"
        assert assignment["dispatched_at"] is not None
        assert assignment["claimed_at"] is not None
    engine.dispose()


def test_qc_metric_migration_upgrades_current_main_without_changing_existing_rows(
    tmp_path,
):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url, "20260712_04")
    current = create_database_engine(database_url)
    with current.begin() as connection:
        connection.execute(
            text(
                """
                INSERT INTO runs (
                    run_id, workflow_id, inputs, status, created_at, updated_at,
                    started_at, ended_at, current_stage, cancellation_reason,
                    error, tags
                ) VALUES (
                    'existing-run', 'workflow', '{}', 'succeeded',
                    :timestamp, :timestamp, :timestamp, :timestamp,
                    'execution', NULL, NULL, '{}'
                )
                """
            ),
            {"timestamp": "2026-07-12 00:00:00"},
        )
        connection.execute(
            text(
                """
                INSERT INTO run_artifacts (
                    artifact_id, run_id, artifact_type, name, uri, mime_type,
                    produced_at, artifact_metadata
                ) VALUES (
                    'artifact-1', 'existing-run', 'file', 'summary.tsv',
                    'run://runs/existing-run/artifacts/artifact-1',
                    'text/tab-separated-values', :timestamp,
                    '{"output_type":"qc_summary"}'
                )
                """
            ),
            {"timestamp": "2026-07-12 00:00:00"},
        )
    current.dispose()

    upgrade_database(database_url)
    upgraded = create_database_engine(database_url)
    with upgraded.connect() as connection:
        assert connection.scalar(text("SELECT count(*) FROM runs")) == 1
        assert connection.scalar(text("SELECT count(*) FROM run_artifacts")) == 1
        assert connection.scalar(text("SELECT count(*) FROM run_qc_metrics")) == 0
        assert (
            connection.scalar(text("SELECT count(*) FROM validated_input_snapshots"))
            == 0
        )
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260827_15"
        )
    upgraded.dispose()


def test_qc_metric_migration_downgrades_independently(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)

    downgrade_database(database_url, "20260712_04")
    downgraded = create_database_engine(database_url)
    inspector = inspect(downgraded)
    assert "run_qc_metrics" not in inspector.get_table_names()
    assert "run_artifacts" in inspector.get_table_names()
    with downgraded.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260712_04"
        )
    downgraded.dispose()


def test_validated_snapshot_migration_upgrades_current_main_without_changing_runs(
    tmp_path,
):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url, "20260712_05")
    current = create_database_engine(database_url)
    with current.begin() as connection:
        connection.execute(
            text(
                """
                INSERT INTO runs (
                    run_id, workflow_id, inputs, status, created_at, updated_at,
                    started_at, ended_at, current_stage, cancellation_reason,
                    error, tags
                ) VALUES (
                    'existing-run', 'workflow', '{}', 'created',
                    :timestamp, :timestamp, NULL, NULL, NULL, NULL, NULL, '{}'
                )
                """
            ),
            {"timestamp": "2026-07-14 00:00:00"},
        )
    current.dispose()

    upgrade_database(database_url)
    upgraded = create_database_engine(database_url)
    with upgraded.connect() as connection:
        assert connection.scalar(text("SELECT count(*) FROM runs")) == 1
        assert (
            connection.scalar(text("SELECT count(*) FROM validated_input_snapshots"))
            == 0
        )
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260827_15"
        )
    upgraded.dispose()


def test_run_history_index_migration_preserves_rows_and_supports_all_query_shapes(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'run-history-upgrade.db'}"
    upgrade_database(database_url, "20260714_06")
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        connection.execute(
            text(
                "INSERT INTO runs "
                "(run_id, workflow_id, inputs, status, created_at, updated_at, "
                "started_at, ended_at, current_stage, cancellation_reason, "
                "error, tags) VALUES "
                "(:run_id, :workflow_id, '{}', :status, :created_at, :updated_at, "
                "NULL, NULL, NULL, NULL, NULL, '{}')"
            ),
            {
                "run_id": "run-1",
                "workflow_id": "workflow-a",
                "status": "created",
                "created_at": "2026-07-14 08:00:00",
                "updated_at": "2026-07-14 08:00:00",
            },
        )
    engine.dispose()

    upgrade_database(database_url)
    upgraded = create_database_engine(database_url)
    expected_indexes = {
        "ix_runs_created_run_id",
        "ix_runs_workflow_created_run_id",
        "ix_runs_status_created_run_id",
        "ix_runs_workflow_status_created_run_id",
    }
    assert expected_indexes <= {
        index["name"] for index in inspect(upgraded).get_indexes("runs")
    }
    with upgraded.connect() as connection:
        assert connection.scalar(text("SELECT count(*) FROM runs")) == 1
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260827_15"
        )
        plans = {
            "ix_runs_created_run_id": (
                "SELECT run_id FROM runs "
                "ORDER BY created_at DESC, run_id DESC LIMIT 51",
                {},
            ),
            "ix_runs_workflow_created_run_id": (
                "SELECT run_id FROM runs WHERE workflow_id = :workflow_id "
                "ORDER BY created_at DESC, run_id DESC LIMIT 51",
                {"workflow_id": "workflow-a"},
            ),
            "ix_runs_status_created_run_id": (
                "SELECT run_id FROM runs WHERE status = :status "
                "ORDER BY created_at DESC, run_id DESC LIMIT 51",
                {"status": "created"},
            ),
            "ix_runs_workflow_status_created_run_id": (
                "SELECT run_id FROM runs "
                "WHERE workflow_id = :workflow_id AND status = :status "
                "ORDER BY created_at DESC, run_id DESC LIMIT 51",
                {"workflow_id": "workflow-a", "status": "created"},
            ),
        }
        for index_name, (statement, parameters) in plans.items():
            detail = " ".join(
                str(row[3])
                for row in connection.execute(
                    text(f"EXPLAIN QUERY PLAN {statement}"), parameters
                )
            )
            assert index_name in detail
            assert "USE TEMP B-TREE FOR ORDER BY" not in detail
    upgraded.dispose()

    downgrade_database(database_url, "20260714_06")
    downgraded = create_database_engine(database_url)
    assert expected_indexes.isdisjoint(
        {index["name"] for index in inspect(downgraded).get_indexes("runs")}
    )
    with downgraded.connect() as connection:
        assert connection.scalar(text("SELECT count(*) FROM runs")) == 1
    downgraded.dispose()

    upgrade_database(database_url)


def test_result_generation_migration_keeps_legacy_results_explicitly_unbound(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'legacy-results.db'}"
    upgrade_database(database_url, "20260714_07")
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        connection.execute(
            text(
                "INSERT INTO runs "
                "(run_id, workflow_id, inputs, status, created_at, updated_at, "
                "started_at, ended_at, current_stage, cancellation_reason, error, tags) "
                "VALUES ('legacy-run', 'workflow', '{}', 'succeeded', :now, :now, "
                ":now, :now, 'execution', NULL, NULL, '{}')"
            ),
            {"now": "2026-07-17 00:00:00"},
        )
        connection.execute(
            text(
                "INSERT INTO run_artifacts "
                "(artifact_id, run_id, artifact_type, name, uri, mime_type, "
                "produced_at, artifact_metadata) VALUES "
                "('artifact-1', 'legacy-run', 'file', 'summary.tsv', "
                "'run://runs/legacy-run/artifacts/artifact-1', 'text/plain', :now, "
                '\'{"output_type":"qc_summary"}\')'
            ),
            {"now": "2026-07-17 00:00:00"},
        )
        connection.execute(
            text(
                "INSERT INTO run_qc_metrics "
                "(metric_id, run_id, metric_key, display_name, value_text, unit, "
                "scope, sample_id, experiment_id, assay, qc_flag, "
                "source_artifact_id, produced_at) VALUES "
                "(:metric_id, 'legacy-run', 'sequencing.total_reads', "
                "'Total reads', '10', 'count', 'sample', 'S1', NULL, 'rnaseq', "
                "NULL, 'artifact-1', :now)"
            ),
            {
                "metric_id": "qcmetric-" + "a" * 64,
                "now": "2026-07-17 00:00:00",
            },
        )
    engine.dispose()

    upgrade_database(database_url)
    upgraded = create_database_engine(database_url)
    with upgraded.connect() as connection:
        state = (
            connection.execute(
                text(
                    "SELECT artifact_revision, artifact_generation, qc_revision, "
                    "qc_generation FROM run_result_states WHERE run_id='legacy-run'"
                )
            )
            .mappings()
            .one()
        )
        assert dict(state) == {
            "artifact_revision": 0,
            "artifact_generation": None,
            "qc_revision": 0,
            "qc_generation": None,
        }
        assert connection.scalar(text("SELECT count(*) FROM run_result_attempts")) == 0
        assert (
            connection.scalar(
                text(
                    "SELECT revision FROM run_artifacts "
                    "WHERE run_id='legacy-run' AND artifact_id='artifact-1'"
                )
            )
            is None
        )

    repository = SqlAlchemyRunRepository(create_session_factory(upgraded))
    with pytest.raises(ValueError, match="generation is unbound"):
        repository.list_artifacts("legacy-run")
    with pytest.raises(ValueError, match="generation is unbound"):
        repository.list_qc_metrics("legacy-run")
    upgraded.dispose()


def test_v030_supported_prior_schema_upgrade_preserves_complete_product_record(
    tmp_path,
) -> None:
    """Revision 07 is the supported pre-generation schema fixture."""
    database_url = f"sqlite:///{tmp_path / 'v030-upgrade.db'}"
    upgrade_database(database_url, "20260714_07")
    engine = create_database_engine(database_url)
    now = "2026-07-17 12:00:00"
    run_id = "release-upgrade-run"
    workflow_id = "bulk-rnaseq"
    artifact_id = "release-artifact"
    metric_id = "qcmetric-" + "b" * 64
    snapshot_id = "vsnap_" + "c" * 32
    build_digest = "d" * 64
    payload = json.dumps(
        {
            "config": {"standard": {}},
            "options": {},
            "samples": [{"sample": "S1"}],
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    payload_digest = build_workflow_inputs_digest(payload)
    with engine.begin() as connection:
        connection.execute(
            text(
                """
                INSERT INTO runs (
                    run_id, workflow_id, inputs, status, created_at, updated_at,
                    started_at, ended_at, current_stage, cancellation_reason,
                    error, tags
                ) VALUES (
                    :run_id, :workflow_id, :inputs, 'succeeded', :now, :now,
                    :now, :now, 'results', NULL, NULL, :tags
                )
                """
            ),
            {
                "run_id": run_id,
                "workflow_id": workflow_id,
                "inputs": payload,
                "now": now,
                "tags": '{"release_fixture":"v0.3.0"}',
            },
        )
        connection.execute(
            text(
                """
                INSERT INTO run_events (
                    event_id, run_id, sequence, event_type, timestamp, status,
                    stage, message, context, issue
                ) VALUES (
                    'event-1', :run_id, 1, 'status_changed', :now, 'succeeded',
                    'results', 'Results indexed.', :context, NULL
                )
                """
            ),
            {"run_id": run_id, "now": now, "context": '{"attempt_id":"attempt-1"}'},
        )
        connection.execute(
            text(
                """
                INSERT INTO run_logs (
                    chunk_id, run_id, stream_name, sequence, timestamp, lines
                ) VALUES (
                    'chunk-1', :run_id, 'stdout', 1, :now, :lines
                )
                """
            ),
            {"run_id": run_id, "now": now, "lines": '["workflow complete"]'},
        )
        connection.execute(
            text(
                """
                INSERT INTO run_execution_assignments (
                    run_id, job_id, backend, queue_name, created_at,
                    dispatched_at, claimed_at, cancellation_requested_at,
                    cancellation_reason, cancellation_acknowledged_at
                ) VALUES (
                    :run_id, 'release-job', 'rq', 'release-queue', :now,
                    :now, :now, NULL, NULL, NULL
                )
                """
            ),
            {"run_id": run_id, "now": now},
        )
        connection.execute(
            text(
                """
                INSERT INTO run_workflow_build_identities (
                    run_id, workflow_id, adapter_version, scheme,
                    logical_entrypoint, digest, captured_at
                ) VALUES (
                    :run_id, :workflow_id, '0.3.0', 'sha256-tree-v1',
                    'main.nf', :digest, :now
                )
                """
            ),
            {
                "run_id": run_id,
                "workflow_id": workflow_id,
                "digest": build_digest,
                "now": now,
            },
        )
        connection.execute(
            text(
                """
                INSERT INTO validated_input_snapshots (
                    snapshot_id, workflow_id, adapter_version, schema_version,
                    schema_dialect, canonical_payload, payload_digest_scheme,
                    payload_digest, validation_outcome, validation_issue_codes,
                    validated_at, expires_at, build_adapter_version, build_scheme,
                    build_logical_entrypoint, build_digest, build_captured_at,
                    consumed_run_id, consumed_at
                ) VALUES (
                    :snapshot_id, :workflow_id, '0.3.0', '1.0.0',
                    'https://json-schema.org/draft/2020-12/schema',
                    :payload, 'sha256-framed-workflow-inputs-v1', :payload_digest,
                    'adapter_validation_succeeded', :issues, :now, :expires,
                    '0.3.0', 'sha256-tree-v1', 'main.nf', :build_digest, :now,
                    :run_id, :now
                )
                """
            ),
            {
                "snapshot_id": snapshot_id,
                "workflow_id": workflow_id,
                "payload": payload,
                "payload_digest": payload_digest,
                "build_digest": build_digest,
                "issues": "[]",
                "now": now,
                "expires": "2026-07-18 12:00:00",
                "run_id": run_id,
            },
        )
        connection.execute(
            text(
                """
                INSERT INTO run_artifacts (
                    artifact_id, run_id, artifact_type, name, uri, mime_type,
                    produced_at, artifact_metadata
                ) VALUES (
                    :artifact_id, :run_id, 'file', 'summary.tsv',
                    :uri, 'text/tab-separated-values', :now, :metadata
                )
                """
            ),
            {
                "artifact_id": artifact_id,
                "run_id": run_id,
                "uri": f"run://runs/{run_id}/artifacts/{artifact_id}",
                "now": now,
                "metadata": '{"output_type":"qc_summary","sample_id":"S1"}',
            },
        )
        connection.execute(
            text(
                """
                INSERT INTO run_qc_metrics (
                    metric_id, run_id, metric_key, display_name, value_text,
                    unit, scope, sample_id, experiment_id, assay, qc_flag,
                    source_artifact_id, produced_at
                ) VALUES (
                    :metric_id, :run_id, 'sequencing.total_reads', 'Total reads',
                    '10', 'count', 'sample', 'S1', NULL, 'rnaseq', 'pass',
                    :artifact_id, :now
                )
                """
            ),
            {
                "metric_id": metric_id,
                "run_id": run_id,
                "artifact_id": artifact_id,
                "now": now,
            },
        )
    engine.dispose()

    upgrade_database(database_url)
    upgraded = create_database_engine(database_url)
    with upgraded.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260827_15"
        )
        run = (
            connection.execute(
                text(
                    "SELECT workflow_id, status, current_stage, tags "
                    "FROM runs WHERE run_id=:run_id"
                ),
                {"run_id": run_id},
            )
            .mappings()
            .one()
        )
        assert run["workflow_id"] == workflow_id
        assert run["status"] == "succeeded"
        assert run["current_stage"] == "results"
        assert json.loads(run["tags"]) == {"release_fixture": "v0.3.0"}
        event = (
            connection.execute(
                text(
                    "SELECT event_id, sequence, event_type, status, stage, message, "
                    "context, issue FROM run_events WHERE run_id=:run_id"
                ),
                {"run_id": run_id},
            )
            .mappings()
            .one()
        )
        assert dict(event) == {
            "event_id": "event-1",
            "sequence": 1,
            "event_type": "status_changed",
            "status": "succeeded",
            "stage": "results",
            "message": "Results indexed.",
            "context": '{"attempt_id":"attempt-1"}',
            "issue": None,
        }
        log = (
            connection.execute(
                text(
                    "SELECT chunk_id, stream_name, sequence, lines "
                    "FROM run_logs WHERE run_id=:run_id"
                ),
                {"run_id": run_id},
            )
            .mappings()
            .one()
        )
        assert dict(log) == {
            "chunk_id": "chunk-1",
            "stream_name": "stdout",
            "sequence": 1,
            "lines": '["workflow complete"]',
        }
        assignment = (
            connection.execute(
                text(
                    "SELECT job_id, backend, queue_name, dispatched_at, claimed_at, "
                    "cancellation_requested_at, cancellation_reason, "
                    "cancellation_acknowledged_at FROM run_execution_assignments "
                    "WHERE run_id=:run_id"
                ),
                {"run_id": run_id},
            )
            .mappings()
            .one()
        )
        assert assignment["job_id"] == "release-job"
        assert assignment["backend"] == "rq"
        assert assignment["queue_name"] == "release-queue"
        assert assignment["dispatched_at"] is not None
        assert assignment["claimed_at"] is not None
        assert assignment["cancellation_requested_at"] is None
        assert assignment["cancellation_reason"] is None
        assert assignment["cancellation_acknowledged_at"] is None
        build_identity = (
            connection.execute(
                text(
                    "SELECT workflow_id, adapter_version, scheme, logical_entrypoint, "
                    "digest FROM run_workflow_build_identities "
                    "WHERE run_id=:run_id"
                ),
                {"run_id": run_id},
            )
            .mappings()
            .one()
        )
        assert dict(build_identity) == {
            "workflow_id": workflow_id,
            "adapter_version": "0.3.0",
            "scheme": "sha256-tree-v1",
            "logical_entrypoint": "main.nf",
            "digest": build_digest,
        }
        snapshot = (
            connection.execute(
                text(
                    "SELECT workflow_id, adapter_version, schema_version, "
                    "canonical_payload, payload_digest_scheme, payload_digest, "
                    "validation_outcome, validation_issue_codes, "
                    "build_adapter_version, build_scheme, "
                    "build_logical_entrypoint, build_digest, consumed_run_id "
                    "FROM validated_input_snapshots "
                    "WHERE snapshot_id=:snapshot_id"
                ),
                {"snapshot_id": snapshot_id},
            )
            .mappings()
            .one()
        )
        assert snapshot["workflow_id"] == workflow_id
        assert snapshot["adapter_version"] == "0.3.0"
        assert snapshot["schema_version"] == "1.0.0"
        assert json.loads(snapshot["canonical_payload"]) == {
            "config": {"standard": {}},
            "samples": [{"sample": "S1"}],
            "options": {},
        }
        assert snapshot["payload_digest_scheme"] == "sha256-framed-workflow-inputs-v1"
        assert snapshot["payload_digest"] == payload_digest
        assert snapshot["validation_outcome"] == "adapter_validation_succeeded"
        assert json.loads(snapshot["validation_issue_codes"]) == []
        assert snapshot["build_adapter_version"] == "0.3.0"
        assert snapshot["build_scheme"] == "sha256-tree-v1"
        assert snapshot["build_logical_entrypoint"] == "main.nf"
        assert snapshot["build_digest"] == build_digest
        assert snapshot["consumed_run_id"] == run_id
        state = (
            connection.execute(
                text(
                    "SELECT artifact_revision, artifact_generation, qc_revision, "
                    "qc_generation FROM run_result_states WHERE run_id=:run_id"
                ),
                {"run_id": run_id},
            )
            .mappings()
            .one()
        )
        assert dict(state) == {
            "artifact_revision": 0,
            "artifact_generation": None,
            "qc_revision": 0,
            "qc_generation": None,
        }
        artifact = (
            connection.execute(
                text(
                    "SELECT artifact_type, name, uri, mime_type, "
                    "artifact_metadata, revision FROM run_artifacts "
                    "WHERE run_id=:run_id AND artifact_id=:artifact_id"
                ),
                {"run_id": run_id, "artifact_id": artifact_id},
            )
            .mappings()
            .one()
        )
        assert artifact["artifact_type"] == "file"
        assert artifact["name"] == "summary.tsv"
        assert artifact["uri"] == f"run://runs/{run_id}/artifacts/{artifact_id}"
        assert artifact["mime_type"] == "text/tab-separated-values"
        assert json.loads(artifact["artifact_metadata"]) == {
            "output_type": "qc_summary",
            "sample_id": "S1",
        }
        assert artifact["revision"] is None
        metric = (
            connection.execute(
                text(
                    "SELECT metric_key, display_name, value_text, unit, scope, "
                    "sample_id, experiment_id, assay, qc_flag, source_artifact_id "
                    "FROM run_qc_metrics WHERE run_id=:run_id"
                ),
                {"run_id": run_id},
            )
            .mappings()
            .one()
        )
        assert dict(metric) == {
            "metric_key": "sequencing.total_reads",
            "display_name": "Total reads",
            "value_text": "10",
            "unit": "count",
            "scope": "sample",
            "sample_id": "S1",
            "experiment_id": None,
            "assay": "rnaseq",
            "qc_flag": "pass",
            "source_artifact_id": artifact_id,
        }

    repository = SqlAlchemyRunRepository(create_session_factory(upgraded))
    with pytest.raises(ValueError, match="generation is unbound"):
        repository.list_artifacts(run_id)
    with pytest.raises(ValueError, match="generation is unbound"):
        repository.list_qc_metrics(run_id)
    upgraded.dispose()
