"""RunService keyset pagination and safe failure mapping tests."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

import pytest

from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.run_history import decode_run_history_cursor
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.services.run_repositories import (
    InMemoryRunRepository,
    RunEventDraft,
)
from encode_pipeline.services.runs import (
    RunHistoryCursorInvalidError,
    RunHistoryCursorNotFoundError,
    RunHistoryDataInvalidError,
    RunHistoryFilterInvalidError,
    RunService,
)


NOW = datetime(2026, 7, 14, 8, 0, tzinfo=timezone.utc)


def _record(run_id: str, *, workflow_id: str = "workflow-a") -> RunRecord:
    return RunRecord(
        run_id=run_id,
        workflow_id=workflow_id,
        inputs={"secret": "/private/input"},
        status=RunStatus.SUCCEEDED,
        created_at=NOW,
        updated_at=NOW + timedelta(minutes=1),
        started_at=NOW + timedelta(seconds=1),
        ended_at=NOW + timedelta(minutes=1),
        current_stage="execution",
        cancellation_reason=None,
        error=None,
        tags={"secret": "value"},
    )


def _service(*records: RunRecord) -> tuple[RunService, InMemoryRunRepository]:
    repository = InMemoryRunRepository()
    for record in records:
        repository.create_run(
            record,
            RunEventDraft(
                event_type="status_changed",
                message="Run created.",
                status=record.status,
            ),
        )
    return RunService(registry=WorkflowRegistry(), repository=repository), repository


def test_service_returns_bounded_page_and_filter_bound_cursor() -> None:
    service, _ = _service(
        _record("run-a"),
        _record("run-c"),
        _record("run-b", workflow_id="workflow-b"),
    )

    page = service.list_run_history(limit=1, workflow_id="workflow-a")

    assert [run.run_id for run in page.runs] == ["run-c"]
    assert page.next_cursor is not None
    cursor = decode_run_history_cursor(
        page.next_cursor,
        workflow_id="workflow-a",
        status=None,
    )
    assert cursor.run_id == "run-c"
    second = service.list_run_history(
        after=page.next_cursor,
        limit=1,
        workflow_id="workflow-a",
    )
    assert [run.run_id for run in second.runs] == ["run-a"]
    assert second.next_cursor is None


def test_service_rejects_invalid_filters_and_cross_filter_cursor() -> None:
    service, _ = _service(_record("run-a"), _record("run-b"))
    page = service.list_run_history(limit=1, workflow_id="workflow-a")

    with pytest.raises(RunHistoryFilterInvalidError):
        service.list_run_history(limit=101)
    with pytest.raises(RunHistoryFilterInvalidError):
        service.list_run_history(workflow_id="/private")
    with pytest.raises(RunHistoryCursorInvalidError):
        service.list_run_history(after="not-a-cursor")
    with pytest.raises(RunHistoryCursorInvalidError):
        service.list_run_history(
            after=page.next_cursor,
            workflow_id="workflow-b",
        )


def test_service_reports_deleted_cursor_boundary_without_exposing_it() -> None:
    service, repository = _service(_record("run-a"), _record("run-b"))
    first = service.list_run_history(limit=1)
    repository._runs.pop("run-b")

    with pytest.raises(RunHistoryCursorNotFoundError):
        service.list_run_history(after=first.next_cursor, limit=1)


def test_service_validates_limit_plus_one_sentinel_before_truncation() -> None:
    service, repository = _service(_record("run-a"), _record("run-b"))
    object.__setattr__(repository._runs["run-a"], "current_stage", "/private")

    with pytest.raises(RunHistoryDataInvalidError):
        service.list_run_history(limit=1)


class _InvalidSummaryRepository(InMemoryRunRepository):
    def list_run_summaries(self, **_kwargs):
        return [object()]


def test_service_fails_closed_for_invalid_repository_collection() -> None:
    service = RunService(
        registry=WorkflowRegistry(),
        repository=_InvalidSummaryRepository(),
    )

    with pytest.raises(RunHistoryDataInvalidError):
        service.list_run_history()
