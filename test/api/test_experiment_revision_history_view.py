"""History lists follow the same revision_view rules as the experiment list.

Default and explicit ``linked`` show the sample catalog frozen on that
experiment revision. ``current_sample_catalog`` shows the sample head, and
adds linked/current revision ids only when those revisions differ. The
history envelope stays ``experiment_revisions`` plus ``next_cursor``.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Callable

import pytest
from pydantic import TypeAdapter
from sqlalchemy import func, select
from sqlalchemy.orm import Session

pytest.importorskip("fastapi")

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from api_test_client import ApiTestClient  # noqa: E402
from conftest import (  # noqa: E402
    TEST_AUTH_PASSWORD_HASH,
    TEST_AUTH_USER_ID,
    seed_test_authentication,
)
from encode_pipeline.api.catalog_models import (  # noqa: E402
    CatalogExperimentResponse,
    CatalogExperimentRevisionListResponse,
    CatalogExperimentRunListResponse,
    CatalogSampleRevisionResponse,
)
from encode_pipeline.api.main import create_app  # noqa: E402
from encode_pipeline.persistence import experiment_catalog as sql_catalog  # noqa: E402
from encode_pipeline.platform.authentication import (  # noqa: E402
    UserAccount,
    UserRole,
    UserStatus,
)
from encode_pipeline.persistence.models import (  # noqa: E402
    ExperimentRevisionRow,
    ExperimentRevisionSampleRow,
    RunExperimentBindingRow,
    RunRow,
    RunSampleRow,
    SampleRevisionRow,
    SnapshotExperimentBindingRow,
    UserAccountRow,
)
from encode_pipeline.services.data_registry import DataRegistryService  # noqa: E402
from encode_pipeline.services.data_registry_repositories import (  # noqa: E402
    InMemoryDataRegistryRepository,
)
from encode_pipeline.services.experiment_catalog import (  # noqa: E402
    ExperimentCatalogService,
)
from encode_pipeline.services.experiment_catalog_repositories import (  # noqa: E402
    InMemoryExperimentCatalogRepository,
)


class _Clock:
    def __init__(self) -> None:
        self.now = datetime(2026, 11, 18, tzinfo=timezone.utc)

    def __call__(self) -> datetime:
        current = self.now
        self.now += timedelta(seconds=1)
        return current


class _Case:
    def __init__(
        self,
        app: Any,
        client: ApiTestClient,
        seed: dict[str, Any],
        snapshot: Callable[[], object],
    ) -> None:
        self.app = app
        self.client = client
        self.seed = seed
        self.snapshot = snapshot

    def close(self) -> None:
        self.app.state.persistence.close()


def _catalog_block(cell_type: dict[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": 1,
        "species": {"status": "unknown"},
        "tissue": {"status": "unknown"},
        "cell_type": cell_type,
        "cell_line": {"status": "not_applicable", "values": []},
        "treatment": {"status": "not_applicable", "values": []},
        "tag": {"status": "not_applicable", "values": []},
    }


def _pin(term: dict[str, Any]) -> dict[str, str]:
    return {
        "term_id": term["term_id"],
        "term_revision_id": term["revision"]["term_revision_id"],
    }


def _known(term: dict[str, Any]) -> dict[str, str]:
    return {"status": "known", **_pin(term)}


def _multi(term: dict[str, Any]) -> dict[str, Any]:
    return {"status": "known", "values": [_pin(term)]}


def _remember(ids: dict[str, str], term: dict[str, Any], label: str) -> None:
    ids[term["term_id"]] = f"TERM_{label}"
    ids[term["revision"]["term_revision_id"]] = f"TRMV_{label}"


def _has_key(value: object, key: str) -> bool:
    if isinstance(value, dict):
        return key in value or any(_has_key(item, key) for item in value.values())
    if isinstance(value, list):
        return any(_has_key(item, key) for item in value)
    return False


def _normalize(body: dict[str, Any], ids: dict[str, str]) -> dict[str, Any]:
    def walk(value: object) -> object:
        if isinstance(value, str):
            return ids.get(value, value)
        if isinstance(value, list):
            return [walk(item) for item in value]
        if isinstance(value, dict):
            return {
                key: walk(item) for key, item in value.items() if key != "next_cursor"
            }
        return value

    normalized = walk(body)
    assert isinstance(normalized, dict)
    normalized["has_next_cursor"] = body.get("next_cursor") is not None
    return normalized


def _assert_no_raw_ids(value: object) -> None:
    if isinstance(value, str):
        assert not value.startswith(
            ("term_", "trmv_", "exp_", "expr_", "smp_", "smpr_", "prj_")
        )
        return
    if isinstance(value, list):
        for item in value:
            _assert_no_raw_ids(item)
        return
    if isinstance(value, dict):
        for item in value.values():
            _assert_no_raw_ids(item)


def _set_role(app: Any, role: str) -> None:
    with Session(app.state.persistence.engine) as session:
        account = session.get(UserAccountRow, TEST_AUTH_USER_ID)
        assert account is not None
        account.role = role
        session.commit()


def _history(client: ApiTestClient, experiment_id: str, **params: object) -> dict:
    response = client.get(
        f"/api/v1/catalog/experiments/{experiment_id}/revisions",
        params=params or None,
    )
    assert response.status_code == 200, response.text
    body = response.json()
    CatalogExperimentRevisionListResponse.model_validate(body)
    return body


def _runs(client: ApiTestClient, project_id: str, revision_id: str) -> dict:
    response = client.get(
        "/api/v1/catalog/projects/"
        f"{project_id}/experiment-revisions/{revision_id}/run-results"
    )
    assert response.status_code == 200, response.text
    body = response.json()
    CatalogExperimentRunListResponse.model_validate(body)
    return body


def _detail(client: ApiTestClient, experiment_id: str, revision_id: str):
    return client.get(
        f"/api/v1/catalog/experiments/{experiment_id}/revisions/{revision_id}"
    )


def _sample(client: ApiTestClient, sample_id: str, revision_id: str) -> dict:
    response = client.get(
        f"/api/v1/catalog/samples/{sample_id}/revisions/{revision_id}"
    )
    assert response.status_code == 200, response.text
    body = response.json()
    TypeAdapter(CatalogSampleRevisionResponse).validate_python(body)
    return body


class _HttpWrites:
    def __init__(self, client: ApiTestClient) -> None:
        self.client = client
        self._n = 0

    def _post(self, path: str, body: dict) -> dict:
        self._n += 1
        response = self.client.post(
            path,
            json=body,
            headers={"Idempotency-Key": f"view-{self._n}"},
        )
        assert response.status_code == 200, response.text
        return response.json()

    def term(self, term_type: str, name: str) -> dict:
        return self._post(
            "/api/v1/catalog/terms",
            {"term_type": term_type, "display_name": name, "aliases": []},
        )

    def revise_sample(
        self,
        sample_id: str,
        expected: int,
        display_name: str,
        catalog: dict,
    ) -> dict:
        return self._post(
            f"/api/v1/catalog/samples/{sample_id}/catalog-revisions",
            {
                "expected_revision_number": expected,
                "display_name": display_name,
                "catalog": catalog,
            },
        )

    def create_experiment(self, project_id: str, body: dict) -> dict:
        return self._post(
            f"/api/v1/catalog/projects/{project_id}/experiments",
            body,
        )

    def revise_experiment(self, experiment_id: str, body: dict) -> dict:
        return self._post(
            f"/api/v1/catalog/experiments/{experiment_id}/revisions",
            body,
        )


class _ServiceWrites:
    def __init__(self, catalog: ExperimentCatalogService) -> None:
        self.catalog = catalog
        self._n = 0

    def _write(self, body: dict) -> dict[str, Any]:
        self._n += 1
        return {
            "actor_user_id": TEST_AUTH_USER_ID,
            "idempotency_key": f"view-{self._n}",
            "raw_body": json.dumps(body, sort_keys=True).encode(),
        }

    def term(self, term_type: str, name: str) -> dict:
        body = {
            "term_type": term_type,
            "display_name": name,
            "namespace": None,
            "accession": None,
            "aliases": [],
        }
        return self.catalog.create_term(**self._write(body), **body)

    def revise_sample(
        self,
        sample_id: str,
        expected: int,
        display_name: str,
        catalog: dict,
    ) -> dict:
        body = {
            "expected_revision_number": expected,
            "display_name": display_name,
            "catalog": catalog,
        }
        return self.catalog.revise_sample_catalog(
            **self._write(body),
            sample_id=sample_id,
            **body,
        )

    def create_experiment(self, project_id: str, body: dict) -> dict:
        return self.catalog.create_experiment(
            **self._write(body),
            project_id=project_id,
            **body,
        )

    def revise_experiment(self, experiment_id: str, body: dict) -> dict:
        return self.catalog.revise_experiment(
            **self._write(body),
            experiment_id=experiment_id,
            **body,
        )


def _seed(registry: DataRegistryService, writes: Any) -> dict[str, Any]:
    project = registry.create_project("Synthetic view")
    cell_x = writes.term("cell_type", "X")
    cell_y = writes.term("cell_type", "Y")
    assay_a = writes.term("assay", "assay-a")
    assay_b = writes.term("assay", "assay-b")
    target = writes.term("target", "target")
    batch = writes.term("batch", "batch")
    role = writes.term("replicate_role", "role")
    created = registry.import_sample(
        project.project_id,
        stable_key="synthetic",
        display_name="synthetic",
        attributes={"organism": "human"},
    )
    frozen = writes.revise_sample(
        created.sample.sample_id,
        1,
        "synthetic",
        _catalog_block(_multi(cell_x)),
    )
    current = writes.revise_sample(
        created.sample.sample_id,
        2,
        "synthetic",
        _catalog_block(_multi(cell_y)),
    )
    assay = {
        "status": "ambiguous",
        "candidates": [_pin(assay_a), _pin(assay_b)],
    }
    samples = [
        {
            "sample_revision_id": frozen["sample_revision_id"],
            "replicate_role": _known(role),
        }
    ]
    first = writes.create_experiment(
        project.project_id,
        {
            "stable_key": "synthetic-exp",
            "display_name": "frozen-v1",
            "assay": assay,
            "target": _known(target),
            "batch": _known(batch),
            "samples": samples,
        },
    )
    second = writes.revise_experiment(
        first["experiment_id"],
        {
            "expected_revision_number": 1,
            "display_name": "frozen-v2",
            "assay": assay,
            "target": _known(target),
            "batch": _known(batch),
            "samples": samples,
        },
    )
    ids = {
        project.project_id: "PROJECT",
        created.sample.sample_id: "SAMPLE",
        created.revision.sample_revision_id: "SMPR_IMPORT",
        frozen["sample_revision_id"]: "SMPR_X",
        current["sample_revision_id"]: "SMPR_Y",
        first["experiment_id"]: "EXP",
        first["experiment_revision_id"]: "EXPR_1",
        second["experiment_revision_id"]: "EXPR_2",
    }
    _remember(ids, cell_x, "X")
    _remember(ids, cell_y, "Y")
    _remember(ids, assay_a, "ASSAY_A")
    _remember(ids, assay_b, "ASSAY_B")
    _remember(ids, target, "TARGET")
    _remember(ids, batch, "BATCH")
    _remember(ids, role, "ROLE")
    return {
        "project_id": project.project_id,
        "sample_id": created.sample.sample_id,
        "import_revision_id": created.revision.sample_revision_id,
        "frozen_revision_id": frozen["sample_revision_id"],
        "current_revision_id": current["sample_revision_id"],
        "experiment_id": first["experiment_id"],
        "revision_ids": (
            first["experiment_revision_id"],
            second["experiment_revision_id"],
        ),
        "cell_x_id": cell_x["term_id"],
        "cell_y_id": cell_y["term_id"],
        "ids": ids,
    }


def _open_app(tmp_path: Path) -> Any:
    application = create_app(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        workspace_root=tmp_path / "workspaces",
    )
    seed_test_authentication(application)
    return application


def _sqlite_snapshot(app: Any, seed: dict[str, Any]) -> object:
    engine = app.state.persistence.engine
    with Session(engine) as session:
        revisions = session.scalars(
            select(ExperimentRevisionRow).where(
                ExperimentRevisionRow.experiment_id == seed["experiment_id"]
            )
        ).all()
        revision_ids = [row.experiment_revision_id for row in revisions]
        pins = session.scalars(
            select(ExperimentRevisionSampleRow).where(
                ExperimentRevisionSampleRow.experiment_revision_id.in_(revision_ids)
            )
        ).all()
        samples = session.scalars(
            select(SampleRevisionRow).where(
                SampleRevisionRow.sample_id == seed["sample_id"]
            )
        ).all()
        counts = tuple(
            session.scalar(select(func.count()).select_from(model))
            for model in (
                RunRow,
                RunSampleRow,
                RunExperimentBindingRow,
                SnapshotExperimentBindingRow,
            )
        )
    return (
        tuple(
            sorted(
                (
                    row.experiment_revision_id,
                    row.revision_number,
                    row.canonical_payload,
                    row.payload_digest,
                )
                for row in revisions
            )
        ),
        tuple(
            sorted(
                (
                    pin.experiment_revision_id,
                    pin.ordinal,
                    pin.sample_id,
                    pin.sample_revision_id,
                    pin.payload_digest,
                )
                for pin in pins
            )
        ),
        tuple(
            sorted(
                (
                    row.sample_revision_id,
                    row.revision_number,
                    row.canonical_payload,
                    row.payload_digest,
                )
                for row in samples
            )
        ),
        counts,
    )


def _open_sqlite(tmp_path: Path) -> _Case:
    app = _open_app(tmp_path)
    app.state.experiment_catalog_service._now_factory = _Clock()
    client = ApiTestClient(app)
    registry = DataRegistryService(
        repository=app.state.persistence.data_registry_repository
    )
    seed = _seed(registry, _HttpWrites(client))
    return _Case(
        app,
        client,
        seed,
        lambda: _sqlite_snapshot(app, seed),
    )


def _open_memory(tmp_path: Path) -> _Case:
    app = _open_app(tmp_path)
    clock = _Clock()
    registry = InMemoryDataRegistryRepository(legacy_created_at=clock.now)
    data = DataRegistryService(repository=registry, now_factory=clock)
    repository = InMemoryExperimentCatalogRepository(registry)
    catalog = ExperimentCatalogService(repository, now_factory=clock)
    seed = _seed(data, _ServiceWrites(catalog))
    app.state.experiment_catalog_service = catalog

    def snapshot() -> object:
        experiment = repository._experiments[seed["experiment_id"]]
        samples = registry.list_sample_revisions(seed["sample_id"])
        return (
            tuple(
                (
                    revision.experiment_revision_id,
                    revision.revision_number,
                    revision.canonical_payload,
                    tuple(
                        (
                            pin.ordinal,
                            pin.sample_id,
                            pin.sample_revision_id,
                            pin.payload_digest,
                        )
                        for pin in revision.payload.samples
                    ),
                )
                for revision in experiment.revisions
            ),
            tuple(
                (
                    revision.sample_revision_id,
                    revision.revision_number,
                    revision.canonical_payload,
                    revision.payload_digest,
                )
                for revision in samples
            ),
        )

    return _Case(app, ApiTestClient(app), seed, snapshot)


def _cell(sample: dict[str, Any]) -> dict[str, Any]:
    return sample["cell_type"]["values"][0]


def _assert_linked_row(row: dict[str, Any], seed: dict[str, Any]) -> None:
    assert row["display_name"] == f"frozen-v{row['revision_number']}"
    sample = row["samples"][0]
    assert _cell(sample)["display_name"] == "X"
    assert _cell(sample)["term_id"] == seed["cell_x_id"]
    assert sample["sample_revision_id"] == seed["frozen_revision_id"]
    assert "linked_sample_revision_id" not in sample
    assert "current_sample_revision_id" not in sample
    assert row["assay"]["status"] == "ambiguous"
    assert "candidates" not in row["assay"]


def _assert_current_row(row: dict[str, Any], seed: dict[str, Any]) -> None:
    assert row["display_name"] == f"frozen-v{row['revision_number']}"
    sample = row["samples"][0]
    assert _cell(sample)["display_name"] == "Y"
    assert _cell(sample)["term_id"] == seed["cell_y_id"]
    assert sample["sample_revision_id"] == seed["current_revision_id"]
    assert sample["linked_sample_revision_id"] == seed["frozen_revision_id"]
    assert sample["current_sample_revision_id"] == seed["current_revision_id"]
    assert row["assay"]["status"] == "ambiguous"
    assert "candidates" not in row["assay"]


def _without_samples(row: dict[str, Any]) -> dict[str, Any]:
    return {key: value for key, value in row.items() if key != "samples"}


def _exercise(case: _Case) -> dict[str, Any]:
    client = case.client
    seed = case.seed
    experiment_id = seed["experiment_id"]
    stored = case.snapshot()
    runs_before = [
        _runs(client, seed["project_id"], revision_id)
        for revision_id in seed["revision_ids"]
    ]
    for body in runs_before:
        assert body["trace_status"] == "unresolved"
        assert body["reason"] == "experiment_snapshot_binding_absent"
        assert body["results"] == []
    detail_before = _detail(client, experiment_id, seed["revision_ids"][0])
    assert detail_before.status_code == 200, detail_before.text
    CatalogExperimentResponse.model_validate(detail_before.json())
    assert len(detail_before.json()["assay"]["candidates"]) == 2
    frozen_before = _sample(client, seed["sample_id"], seed["frozen_revision_id"])
    current_before = _sample(client, seed["sample_id"], seed["current_revision_id"])
    assert _cell(frozen_before)["display_name"] == "X"
    assert _cell(current_before)["display_name"] == "Y"
    assert frozen_before["revision_number"] == 2
    assert current_before["revision_number"] == 3

    default = _history(client, experiment_id)
    linked = _history(client, experiment_id, revision_view="linked")
    current = _history(client, experiment_id, revision_view="current_sample_catalog")
    assert default == linked
    assert [row["revision_number"] for row in default["experiment_revisions"]] == [
        1,
        2,
    ]
    for row in default["experiment_revisions"]:
        _assert_linked_row(row, seed)
        assert not _has_key(row, "candidates")
    for linked_row, current_row in zip(
        linked["experiment_revisions"],
        current["experiment_revisions"],
        strict=True,
    ):
        _assert_current_row(current_row, seed)
        assert not _has_key(current_row, "candidates")
        assert _without_samples(linked_row) == _without_samples(current_row)

    first = _history(client, experiment_id, limit=1)
    assert first["experiment_revisions"][0]["revision_number"] == 1
    assert first["next_cursor"]
    continued = _history(
        client,
        experiment_id,
        limit=1,
        revision_view="linked",
        cursor=first["next_cursor"],
    )
    assert continued["experiment_revisions"][0]["revision_number"] == 2
    assert continued["next_cursor"] is None
    _assert_linked_row(continued["experiment_revisions"][0], seed)
    mixed = client.get(
        f"/api/v1/catalog/experiments/{experiment_id}/revisions",
        params={
            "limit": 1,
            "revision_view": "current_sample_catalog",
            "cursor": first["next_cursor"],
        },
    )
    assert mixed.status_code == 400
    assert mixed.json()["issues"][0]["code"] == "catalog_cursor_invalid"
    current_page = _history(
        client,
        experiment_id,
        limit=1,
        revision_view="current_sample_catalog",
    )
    assert current_page["experiment_revisions"][0]["revision_number"] == 1
    _assert_current_row(current_page["experiment_revisions"][0], seed)
    current_cursor = current_page["next_cursor"]
    assert current_cursor
    for params in (
        {"limit": 1, "cursor": current_cursor},
        {"limit": 1, "revision_view": "linked", "cursor": current_cursor},
    ):
        rejected = client.get(
            f"/api/v1/catalog/experiments/{experiment_id}/revisions",
            params=params,
        )
        assert rejected.status_code == 400
        assert rejected.json()["issues"][0]["code"] == "catalog_cursor_invalid"
    current_next = _history(
        client,
        experiment_id,
        limit=1,
        revision_view="current_sample_catalog",
        cursor=current_cursor,
    )
    assert current_next["experiment_revisions"][0]["revision_number"] == 2
    _assert_current_row(current_next["experiment_revisions"][0], seed)

    assert case.snapshot() == stored
    runs_after = [
        _runs(client, seed["project_id"], revision_id)
        for revision_id in seed["revision_ids"]
    ]
    assert runs_after == runs_before
    detail_after = _detail(client, experiment_id, seed["revision_ids"][0])
    assert detail_after.status_code == 200
    assert detail_after.json() == detail_before.json()
    assert (
        _sample(client, seed["sample_id"], seed["frozen_revision_id"]) == frozen_before
    )
    assert (
        _sample(client, seed["sample_id"], seed["current_revision_id"])
        == current_before
    )
    if len(stored) == 4:
        assert stored[3] == (0, 0, 0, 0)

    case.app.state.authentication_repository.create_account(
        UserAccount(
            user_id="usr_" + "d" * 32,
            username="other-admin",
            role=UserRole.ADMINISTRATOR,
            status=UserStatus.ENABLED,
            password_hash=TEST_AUTH_PASSWORD_HASH,
            created_at=datetime(2026, 11, 18, tzinfo=timezone.utc),
            updated_at=datetime(2026, 11, 18, tzinfo=timezone.utc),
            password_changed_at=datetime(2026, 11, 18, tzinfo=timezone.utc),
        )
    )
    _set_role(case.app, "member")
    try:
        member = _history(
            client,
            experiment_id,
            revision_view="current_sample_catalog",
        )
        assert member == current
        assert not _has_key(member, "candidates")
        denied = _detail(client, experiment_id, seed["revision_ids"][0])
        assert denied.status_code == 403
        assert denied.json()["issues"][0]["code"] == "ADMINISTRATOR_REQUIRED"
    finally:
        _set_role(case.app, "administrator")

    ids = seed["ids"]
    report = {
        "default": _normalize(default, ids),
        "linked": _normalize(linked, ids),
        "current": _normalize(current, ids),
        "linked_page": _normalize(continued, ids),
        "current_page": _normalize(current_next, ids),
        "member": _normalize(member, ids),
    }
    _assert_no_raw_ids(report)
    return report


def test_experiment_revision_history_http_matches_for_both_repositories(
    tmp_path: Path,
) -> None:
    sqlite_case = _open_sqlite(tmp_path / "sqlite")
    try:
        memory_case = _open_memory(tmp_path / "memory")
        try:
            assert _exercise(sqlite_case) == _exercise(memory_case)
        finally:
            memory_case.close()
    finally:
        sqlite_case.close()


def test_sqlite_current_history_keeps_one_snapshot_during_a_sample_revision(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    case = _open_sqlite(tmp_path)
    try:
        seed = case.seed
        frozen = _sample(case.client, seed["sample_id"], seed["frozen_revision_id"])
        cell = _cell(frozen)
        pin = {key: cell[key] for key in ("term_id", "term_revision_id")}
        revision = {
            "expected_revision_number": 3,
            "display_name": "synthetic",
            "catalog": _catalog_block({"status": "known", "values": [pin]}),
        }
        collect_ids = sql_catalog._experiment_name_ids

        def collect_ids_then_revise(*args):
            term_ids = collect_ids(*args)
            # Commit on a second connection after the reader has selected Y.
            case.app.state.experiment_catalog_service.revise_sample_catalog(
                actor_user_id=TEST_AUTH_USER_ID,
                idempotency_key="concurrent-sample-update",
                raw_body=json.dumps(revision).encode(),
                sample_id=seed["sample_id"],
                **revision,
            )
            return term_ids

        query = {"revision_view": "current_sample_catalog", "limit": 1}
        with monkeypatch.context() as patch:
            patch.setattr(sql_catalog, "_experiment_name_ids", collect_ids_then_revise)
            before = _history(case.client, seed["experiment_id"], **query)
        _assert_current_row(before["experiment_revisions"][0], seed)

        after = _history(case.client, seed["experiment_id"], **query)
        sample = after["experiment_revisions"][0]["samples"][0]
        assert _cell(sample)["term_id"] == seed["cell_x_id"]
        assert sample["sample_revision_id"] == sample["current_sample_revision_id"]
        assert sample["sample_revision_id"] != seed["current_revision_id"]
        assert sample["linked_sample_revision_id"] == seed["frozen_revision_id"]
    finally:
        case.close()
