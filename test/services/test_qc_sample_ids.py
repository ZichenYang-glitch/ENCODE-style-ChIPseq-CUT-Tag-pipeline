"""QC sample labels preserve ASCII spaces without widening other identifiers."""

import asyncio
from dataclasses import replace
from decimal import Decimal
from pathlib import Path
import runpy

import httpx
import pytest

from encode_pipeline.api.main import create_app
from encode_pipeline.api.models import QcMetricResponse
from encode_pipeline.platform.adapters import ExtractedQcMetricCandidate
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.runs import (
    build_qc_metric_id,
    validate_qc_identifier_token,
)
from encode_pipeline.services.qc_summary_indexing import QcSummaryIndexingService
from test_atomic_result_bundle import (
    candidates,
    publish,
    service_for,
    store as shared_store,
)

store = shared_store


VALID = ("ab", "a b", "a  b", "a b ", " a b", "A._-9", "a" + " " * 254)
INVALID = (
    "",
    " ",
    "  ",
    ".",
    "..",
    "a/b",
    "a\\b",
    "a:b",
    "a\tb",
    "a\nb",
    "a\rb",
    "a\x00b",
    "a\u00a0b",
    "汉字",
    "x" * 256,
)


def metric_for(sample_id, **changes):
    _artifacts, metrics = candidates()
    return replace(
        metrics[0],
        sample_id=sample_id,
        metric_id=build_qc_metric_id("pets", "sample", sample_id, None),
        **changes,
    )


def candidate_for(sample_id, **changes):
    return ExtractedQcMetricCandidate(
        metric_key="pets",
        display_name="PETs",
        value=Decimal(8),
        unit="count",
        scope="sample",
        sample_id=sample_id,
        source_artifact_id="summary",
        **changes,
    )


@pytest.mark.parametrize("sample_id", VALID)
def test_sample_label_candidate_source_metadata_and_response_preserve_bytes(sample_id):
    candidate = candidate_for(sample_id)
    QcSummaryIndexingService._validate_candidate(candidate, {"summary"})
    assert QcSummaryIndexingService._validated_source_metadata(
        {"sample_id": sample_id, "assay": "Hi-TrAC", "scope": "sample"}
    ) == {"sample_id": sample_id, "assay": "Hi-TrAC", "scope": "sample"}
    response = QcMetricResponse.from_metric(
        metric_for(sample_id), expected_run_id="run-1"
    )
    assert response.sample_id == sample_id
    assert response.model_dump(mode="json")["sample_id"].encode() == sample_id.encode()


@pytest.mark.parametrize("sample_id", INVALID)
def test_sample_label_invalid_characters_rejected_at_each_shared_boundary(sample_id):
    with pytest.raises(ValueError):
        QcSummaryIndexingService._validate_candidate(
            candidate_for(sample_id), {"summary"}
        )
    with pytest.raises(ValueError):
        QcSummaryIndexingService._validated_source_metadata({"sample_id": sample_id})
    with pytest.raises(ValueError):
        QcMetricResponse.from_metric(metric_for(sample_id), expected_run_id="run-1")


@pytest.mark.parametrize("field", ("experiment_id", "assay"))
def test_other_qc_identifiers_still_reject_spaces(field):
    with pytest.raises(ValueError):
        validate_qc_identifier_token("a b")
    with pytest.raises(ValueError):
        QcSummaryIndexingService._validate_candidate(
            candidate_for("sample", **{field: "a b"}), {"summary"}
        )
    with pytest.raises(ValueError):
        QcSummaryIndexingService._validated_source_metadata({field: "a b"})
    with pytest.raises(ValueError):
        QcMetricResponse.from_metric(
            metric_for("sample", **{field: "a b"}), expected_run_id="run-1"
        )


def test_repositories_keep_sample_label_identity_and_trailing_spaces(store):
    service = service_for(store)
    artifacts, _metrics = candidates()
    metrics = tuple(
        metric_for(sample_id) for sample_id in ("ab", "a b", "a  b", "a b ")
    )
    state = publish(service, artifacts, metrics)
    observed = service.list_qc_metrics("run-1")
    assert {m.sample_id for m in observed} == {"ab", "a b", "a  b", "a b "}
    assert len({m.metric_id for m in observed}) == 4
    assert {m.metric_id: m.sample_id.encode() for m in observed} == {
        m.metric_id: m.sample_id.encode() for m in metrics
    }
    assert state.qc_artifact_generation == state.artifact_generation


@pytest.mark.parametrize("sample_id", INVALID)
def test_repositories_reject_invalid_sample_labels_without_publication(
    store, sample_id
):
    service = service_for(store)
    artifacts, _metrics = candidates()
    with pytest.raises(ValueError):
        publish(service, artifacts, (metric_for(sample_id),))
    assert service.list_qc_metrics("run-1") == ()
    assert service.list_artifacts("run-1") == ()


def test_original_http_qc_projection_preserves_sample_label_bytes(
    tmp_path, monkeypatch
):
    # Real temporary SQLite and authentication; seeded contract metrics, no
    # scientific generation claim. The separate real tiny test covers production.
    import os

    for name in tuple(os.environ):
        if name.startswith(("ENCODE_PIPELINE_", "HELIXWEAVE_")):
            monkeypatch.delenv(name)
    monkeypatch.setenv("HELIXWEAVE_TERMINAL_EMAIL_ENABLED", "false")
    monkeypatch.setenv("ENCODE_PIPELINE_REDIS_URL", "redis://127.0.0.1:1/0")
    root = Path(__file__).resolve().parents[2]
    app = create_app(
        database_url=f"sqlite:///{tmp_path / 'api.sqlite'}",
        workspace_root=tmp_path / "workspace",
        project_root=root,
        registry=WorkflowRegistry(),
    )
    try:
        helpers = runpy.run_path(str(root / "test/conftest.py"))
        helpers["seed_test_authentication"](app)
        service_for(app.state.persistence.repository)
        artifacts, _ = candidates()
        names = ("ab", "a b", "a  b", "a b ")
        metrics = tuple(metric_for(name) for name in names)
        publish(app.state.run_service, artifacts, metrics)

        async def request():
            transport = httpx.ASGITransport(app=app, raise_app_exceptions=False)
            async with httpx.AsyncClient(
                transport=transport,
                base_url="http://test",
                cookies={
                    app.state.auth_cookie_policy.session_cookie.name: app.state.test_auth_tokens[
                        0
                    ]
                },
            ) as client:
                response = await client.get("/api/v1/runs/run-1/qc-metrics")
                assert response.status_code == 200, response.text
                returned = response.json()["qc_metrics"]
                assert {m["sample_id"] for m in returned} == set(names)
                assert {m["metric_id"]: m["sample_id"].encode() for m in returned} == {
                    m.metric_id: m.sample_id.encode() for m in metrics
                }

        asyncio.run(request())
    finally:
        app.state.run_queue.close()
        app.state.persistence.close()
