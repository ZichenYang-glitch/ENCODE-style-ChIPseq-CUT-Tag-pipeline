"""Parity tests for bounded, summary-only run-history persistence."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

import pytest
from sqlalchemy import event, text

from encode_pipeline.persistence import (
    SqlAlchemyRunRepository,
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.platform.run_history import RunHistoryCursor
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.services.run_repositories import (
    InMemoryRunRepository,
    RunEventDraft,
)


NOW = datetime(2026, 7, 14, 8, 0, tzinfo=timezone.utc)


@pytest.fixture(params=["memory", "sqlite"])
def repository(request, tmp_path):
    if request.param == "memory":
        yield InMemoryRunRepository()
        return
    database_url = f"sqlite:///{tmp_path / 'history.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    yield SqlAlchemyRunRepository(create_session_factory(engine))
    engine.dispose()


def _record(
    run_id: str,
    *,
    created_at: datetime = NOW,
    workflow_id: str = "workflow-a",
    status: RunStatus = RunStatus.SUCCEEDED,
) -> RunRecord:
    ended_at = created_at + timedelta(minutes=1) if status.is_terminal else None
    started_at = (
        created_at + timedelta(seconds=1) if status is not RunStatus.CREATED else None
    )
    updated_at = ended_at or started_at or created_at
    return RunRecord(
        run_id=run_id,
        workflow_id=workflow_id,
        inputs={"private_path": "/never/select/me"},
        status=status,
        created_at=created_at,
        updated_at=updated_at,
        started_at=started_at,
        ended_at=ended_at,
        current_stage="execution" if started_at is not None else None,
        cancellation_reason=None,
        error=None,
        tags={"secret": "not-public"},
    )


def _create(repository, record: RunRecord) -> None:
    repository.create_run(
        record,
        RunEventDraft(
            event_type="status_changed",
            message="Run created.",
            status=record.status,
        ),
    )


def test_summary_pages_are_stable_descending_and_filterable(repository) -> None:
    _create(repository, _record("run-a", workflow_id="workflow-a"))
    _create(repository, _record("run-c", workflow_id="workflow-a"))
    _create(
        repository,
        _record(
            "run-b",
            workflow_id="workflow-b",
            status=RunStatus.FAILED,
        ),
    )
    _create(
        repository,
        _record("run-old", created_at=NOW - timedelta(minutes=1)),
    )

    first = repository.list_run_summaries(limit=2)
    assert [run.run_id for run in first] == ["run-c", "run-b"]
    cursor = RunHistoryCursor(
        created_at=first[-1].created_at,
        run_id=first[-1].run_id,
        workflow_id=None,
        status=None,
    )
    assert [
        run.run_id for run in repository.list_run_summaries(after=cursor, limit=3)
    ] == ["run-a", "run-old"]
    assert [
        run.run_id
        for run in repository.list_run_summaries(
            limit=10,
            workflow_id="workflow-a",
        )
    ] == ["run-c", "run-a", "run-old"]
    assert [
        run.run_id
        for run in repository.list_run_summaries(
            limit=10,
            status=RunStatus.FAILED,
        )
    ] == ["run-b"]


def test_same_time_insertions_obey_keyset_without_repeating_prior_rows(
    repository,
) -> None:
    for run_id in ("run-z", "run-m", "run-c"):
        _create(repository, _record(run_id))
    first = repository.list_run_summaries(limit=2)
    assert [run.run_id for run in first] == ["run-z", "run-m"]
    cursor = RunHistoryCursor(
        created_at=NOW,
        run_id="run-m",
        workflow_id=None,
        status=None,
    )

    _create(repository, _record("run-y"))
    _create(repository, _record("run-a"))

    later = repository.list_run_summaries(after=cursor, limit=10)
    assert [run.run_id for run in later] == ["run-c", "run-a"]
    assert not {run.run_id for run in first} & {run.run_id for run in later}


def test_cursor_must_belong_to_exact_filters(repository) -> None:
    _create(repository, _record("run-a"))
    cursor = RunHistoryCursor(
        created_at=NOW,
        run_id="run-a",
        workflow_id="workflow-a",
        status=None,
    )

    with pytest.raises(ValueError, match="filters"):
        repository.list_run_summaries(after=cursor, limit=2)
    failed_cursor = RunHistoryCursor(
        created_at=NOW,
        run_id="run-a",
        workflow_id="workflow-a",
        status=RunStatus.FAILED,
    )
    with pytest.raises(KeyError):
        repository.list_run_summaries(
            after=failed_cursor,
            limit=2,
            workflow_id="workflow-a",
            status=RunStatus.FAILED,
        )


def test_in_memory_rejects_mapping_key_and_selected_row_corruption() -> None:
    repository = InMemoryRunRepository()
    _create(repository, _record("run-good"))
    repository._runs["wrong-key"] = repository._runs.pop("run-good")

    with pytest.raises(ValueError, match="identity"):
        repository.list_run_summaries(limit=2)

    repository._runs.clear()
    _create(repository, _record("run-good"))
    object.__setattr__(repository._runs["run-good"], "workflow_id", "/private")
    with pytest.raises(ValueError, match="public-safe"):
        repository.list_run_summaries(limit=2)


def test_sqlalchemy_summary_query_selects_only_public_columns_and_uses_limit(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'selected.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    repository = SqlAlchemyRunRepository(create_session_factory(engine))
    _create(repository, _record("run-a"))
    statements: list[str] = []

    def capture(_connection, _cursor, statement, _parameters, _context, _many):
        statements.append(statement)

    event.listen(engine, "before_cursor_execute", capture)
    try:
        repository.list_run_summaries(limit=1)
    finally:
        event.remove(engine, "before_cursor_execute", capture)
        engine.dispose()

    statement = next(value for value in statements if "FROM runs" in value)
    assert "runs.inputs" not in statement
    assert "runs.error" not in statement
    assert "runs.tags" not in statement
    assert " LIMIT " in statement


def test_sqlalchemy_rejects_corrupt_selected_summary_row(tmp_path) -> None:
    database_url = f"sqlite:///{tmp_path / 'corrupt.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    repository = SqlAlchemyRunRepository(create_session_factory(engine))
    _create(repository, _record("run-a"))
    with engine.begin() as connection:
        connection.execute(
            text("UPDATE runs SET current_stage = :stage WHERE run_id = :run_id"),
            {"stage": "/private/stage", "run_id": "run-a"},
        )

    with pytest.raises(ValueError, match="public-safe"):
        repository.list_run_summaries(limit=2)
    engine.dispose()
