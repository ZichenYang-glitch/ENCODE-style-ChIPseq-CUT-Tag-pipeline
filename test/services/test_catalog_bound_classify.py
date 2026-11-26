"""Ledger and worker counts for experiment.target classification.

The worker is an in-process callable. Nothing here opens a provider or a key.
"""

from __future__ import annotations

from contextlib import closing
from pathlib import Path
import time

import json

import pytest

from encode_pipeline.services.authoring_bridge import BridgeStatus, BridgeTurnRequest
from encode_pipeline.services.catalog_model_client import load_catalog_model
from encode_pipeline.services.peri_runtime.bound_model import (
    BoundModel,
    Budget,
    binding_id,
)
from encode_pipeline.services.peri_runtime.catalog_classify import (
    CatalogTransportError,
    bind_target_candidates,
    loopback_decision,
)
from encode_pipeline.services.peri_runtime.owner import NativeAuthoringBridge
from encode_pipeline.services.peri_runtime.socket_server import OwnerServer


PROJECT = "prj_" + "ab" * 16
OTHER = "prj_" + "cd" * 16
DIGEST = "12" * 32
OTHER_DIGEST = "34" * 32
CANDIDATES = (("h3k4me3", "H3K4me3"), ("h3k27me3", "H3K27me3"))


class Worker:
    def __init__(self) -> None:
        self.binding = {"model": "bound-a", "protocol": "anthropic"}
        self.sends: list[dict] = []
        self.binding_calls = 0
        self.fail: str | None = None

    def __call__(self, request: dict) -> dict:
        if request.get("op") == "binding":
            self.binding_calls += 1
            return dict(self.binding)
        self.sends.append(request)
        if self.fail == "timeout":
            raise CatalogTransportError("catalog_model_timeout")
        if self.fail == "index":
            return {"binding": dict(self.binding), "decision": "1"}
        if request.get("op") == "suggest":
            return {"binding": dict(self.binding), "value": 7}
        decision = loopback_decision(request["clue"], request["candidates"])
        return {"binding": dict(self.binding), "decision": decision}


def _budget(requests: int) -> Budget:
    return Budget(requests, 4096 * requests, 64, 100, 16, 100, 100, "CNY")


def _balances(model: BoundModel, grant_id: str) -> tuple[int, int, int]:
    with closing(model.db()) as db:
        row = db.execute(
            """SELECT requests,
                      (SELECT count(*) FROM catalog_reservations
                       WHERE grant_id=?),
                      (SELECT coalesce(sum(delivered), 0) FROM catalog_reservations
                       WHERE grant_id=?)
               FROM catalog_grants WHERE grant_id=?""",
            (grant_id, grant_id, grant_id),
        ).fetchone()
    assert row is not None
    return int(row[0]), int(row[1]), int(row[2])


@pytest.fixture
def bound(tmp_path: Path):
    worker = Worker()
    model = BoundModel(tmp_path / "ledger.sqlite", worker)
    BoundModel.initialize(model.ledger)
    return worker, model


def _grant(model: BoundModel, worker: Worker, project: str, requests: int = 2) -> str:
    return model.catalog_grant(
        project_id=project,
        field="experiment.target",
        binding=binding_id(worker.binding),
        expires=int(time.time()) + 300,
        budget=_budget(requests),
    )


def _classify(model: BoundModel, **overrides: object) -> str:
    payload = {
        "project_id": PROJECT,
        "idempotency_key": "classify-one",
        "clue": "H3K27me3",
        "candidates": CANDIDATES,
        "projection_digest": DIGEST,
    }
    payload.update(overrides)
    return model.classify(**payload)


@pytest.mark.parametrize(
    "label",
    ["/home/synthetic/run1", "https://example.invalid", "sample@example.invalid"],
)
def test_unsafe_target_label_can_only_use_a_safe_vocabulary_alias(bound, label):
    worker, model = bound
    grant_id = _grant(model, worker, PROJECT, requests=1)
    term = {"display_name": label, "term_id": "t1", "term_revision_id": "v1"}
    assert bind_target_candidates([term], explicit_revision_id="v1") == ()
    eligible = bind_target_candidates(
        [{**term, "aliases": ["H2A.Z"]}], explicit_revision_id="v1"
    )
    assert [(item.label, item.token) for item in eligible] == [("H2A.Z", "h2az")]
    for token, unsafe in (("h2az", label), ("h2az", "CTCF")):
        with pytest.raises(CatalogTransportError) as refused:
            _classify(model, candidates=((token, unsafe),))
        assert refused.value.code == "catalog_model_malformed"
    assert worker.sends == []
    assert _balances(model, grant_id) == (1, 0, 0)


def test_classify_requires_its_own_grant_and_replays_without_another_send(bound):
    worker, model = bound
    with pytest.raises(CatalogTransportError) as unauthorized:
        _classify(model)
    assert unauthorized.value.code == "catalog_model_unauthorized"
    assert worker.sends == []
    assert worker.binding_calls == 0

    model.grant(
        sid="sid",
        epoch="epoch",
        instance="instance",
        generation=1,
        field="config.threads",
        binding=binding_id(worker.binding),
        expires=int(time.time()) + 300,
        budget=_budget(1),
    )
    with pytest.raises(CatalogTransportError) as still_closed:
        _classify(model)
    assert still_closed.value.code == "catalog_model_unauthorized"
    assert worker.sends == []
    request = BridgeTurnRequest(
        "sid",
        "turn",
        1,
        "instance",
        "epoch",
        "secret",
        {
            "intent": "suggest",
            "field_tokens": ["config.threads"],
            "preferences": [{"field_token": "config.threads", "value": 3}],
        },
        ("secret",),
        "schema",
        "projection",
        1,
    )
    suggested = json.loads(model.suggest(request, lambda: True))
    assert suggested["changes"][0]["value"] == 7
    assert len(worker.sends) == 1
    assert worker.sends[0]["op"] == "suggest"

    grant_id = _grant(model, worker, PROJECT, requests=1)
    assert _classify(model) == "h3k27me3"
    sent = worker.sends[-1]
    assert sent["op"] == "classify"
    assert sent["field"] == "experiment.target"
    assert "term_revision_id" not in sent
    assert [item["token"] for item in sent["candidates"]] == [
        "h3k4me3",
        "h3k27me3",
    ]
    assert _classify(model) == "h3k27me3"
    assert len(worker.sends) == 2
    # suggest and classify each read the binding before and after the call.
    assert worker.binding_calls == 4
    requests_left, reservations, delivered = _balances(model, grant_id)
    assert (requests_left, reservations, delivered) == (0, 1, 1)

    with pytest.raises(CatalogTransportError) as exhausted:
        _classify(model, idempotency_key="classify-two")
    assert exhausted.value.code == "catalog_model_quota_exhausted"
    assert len(worker.sends) == 2


def test_cross_project_revocation_expiry_and_timeout_keep_their_counts(bound):
    worker, model = bound
    grant_id = _grant(model, worker, PROJECT, requests=2)
    with pytest.raises(CatalogTransportError) as mismatch:
        _classify(model, project_id=OTHER, idempotency_key="other-project")
    assert mismatch.value.code == "catalog_model_project_mismatch"
    assert worker.sends == []
    assert worker.binding_calls == 0
    assert _balances(model, grant_id) == (2, 0, 0)

    with closing(model.db()) as db, db:
        db.execute(
            "UPDATE catalog_grants SET expires=? WHERE grant_id=?",
            (int(time.time()) - 5, grant_id),
        )
    with pytest.raises(CatalogTransportError) as expired:
        _classify(model, idempotency_key="after-expiry")
    assert expired.value.code == "catalog_model_expired"
    assert worker.sends == []
    assert worker.binding_calls == 1
    assert _balances(model, grant_id) == (2, 0, 0)

    replacement = _grant(model, worker, PROJECT, requests=2)
    model.catalog_revoke(replacement)
    with pytest.raises(CatalogTransportError) as revoked:
        _classify(model, idempotency_key="after-revoke")
    assert revoked.value.code == "catalog_model_revoked"
    assert worker.sends == []
    assert worker.binding_calls == 2
    assert _balances(model, replacement) == (2, 0, 0)

    live = _grant(model, worker, OTHER, requests=1)
    worker.fail = "timeout"
    with pytest.raises(CatalogTransportError) as timed_out:
        _classify(
            model,
            project_id=OTHER,
            idempotency_key="timed-out",
            projection_digest=OTHER_DIGEST,
        )
    assert timed_out.value.code == "catalog_model_timeout"
    assert len(worker.sends) == 1
    assert _balances(model, live) == (0, 1, 0)
    with pytest.raises(CatalogTransportError) as retried:
        _classify(
            model,
            project_id=OTHER,
            idempotency_key="timed-out",
            projection_digest=OTHER_DIGEST,
        )
    assert retried.value.code == "catalog_model_timeout"
    assert len(worker.sends) == 1
    assert _balances(model, live) == (0, 1, 0)


def test_ordinal_and_index_answers_are_not_vocabulary_choices(bound):
    worker, model = bound
    grant_id = _grant(model, worker, PROJECT, requests=2)
    with pytest.raises(CatalogTransportError) as ordinal:
        _classify(model, clue="item-1", idempotency_key="ordinal")
    assert ordinal.value.code == "catalog_model_malformed"
    assert worker.sends == []
    assert worker.binding_calls == 0
    assert _balances(model, grant_id) == (2, 0, 0)

    assert (
        loopback_decision(
            "1",
            [{"token": "h3k4me3", "label": "H3K4me3"}],
        )
        == "insufficient_evidence"
    )
    assert (
        loopback_decision(
            "H3K27me3",
            [{"token": token, "label": label} for token, label in CANDIDATES],
        )
        == "h3k27me3"
    )

    worker.fail = "index"
    with pytest.raises(CatalogTransportError) as indexed:
        _classify(model, idempotency_key="indexed")
    assert indexed.value.code == "catalog_model_malformed"
    assert len(worker.sends) == 1
    assert _balances(model, grant_id) == (1, 1, 0)
    with pytest.raises(CatalogTransportError) as again:
        _classify(model, idempotency_key="indexed")
    assert again.value.code == "catalog_model_timeout"
    assert len(worker.sends) == 1


def test_changed_projection_conflicts_without_a_second_charge(bound):
    worker, model = bound
    grant_id = _grant(model, worker, PROJECT, requests=2)
    assert _classify(model) == "h3k27me3"
    with pytest.raises(CatalogTransportError) as conflict:
        _classify(model, clue="H3K4me3", projection_digest=OTHER_DIGEST)
    assert conflict.value.code == "catalog_idempotency_conflict"
    assert len(worker.sends) == 1
    assert _balances(model, grant_id) == (1, 1, 1)


def test_owner_socket_forwards_classify_and_api_stays_unbound_by_default(tmp_path):
    assert load_catalog_model({}) is None
    worker = Worker()
    model = BoundModel(tmp_path / "owner-ledger.sqlite", worker)
    BoundModel.initialize(model.ledger)
    _grant(model, worker, PROJECT, requests=1)
    bridge = NativeAuthoringBridge(
        binary=tmp_path / "missing-peri",
        state_root=tmp_path / "owner-state",
        bound_model=model,
    )
    unbound = NativeAuthoringBridge(
        binary=tmp_path / "missing-peri",
        state_root=tmp_path / "empty-state",
    )
    with pytest.raises(CatalogTransportError) as missing:
        unbound.classify_catalog(
            {
                "project_id": PROJECT,
                "idempotency_key": "owner-call",
                "clue": "H3K4me3",
                "candidates": [{"token": "h3k4me3", "label": "H3K4me3"}],
                "projection_digest": DIGEST,
            }
        )
    assert missing.value.code == "catalog_model_unbound"
    decision = bridge.classify_catalog(
        {
            "project_id": PROJECT,
            "idempotency_key": "owner-call",
            "clue": "H3K4me3",
            "candidates": [{"token": "h3k4me3", "label": "H3K4me3"}],
            "projection_digest": DIGEST,
        }
    )
    assert decision == {"decision": "h3k4me3"}
    assert len(worker.sends) == 1

    class SocketBridge:
        def status(self) -> BridgeStatus:
            return BridgeStatus(available=True)

        def owner_epoch(self) -> str:
            return "epoch-1"

        def classify_catalog(self, payload: dict) -> dict[str, str]:
            if payload.get("clue") == "bad":
                raise CatalogTransportError("catalog_model_quota_exhausted")
            return {"decision": "h3k4me3"}

    server = OwnerServer(SocketBridge(), object(), tmp_path / "owner.sock")
    assert server.handle(
        {"method": "classify_catalog", "owner_epoch": "epoch-1", "clue": "H3K4me3"}
    ) == {"ok": True, "decision": "h3k4me3"}
    assert (
        server.handle(
            {"method": "classify_catalog", "owner_epoch": "stale", "clue": "H3K4me3"}
        )["error"]
        == "epoch_rejected"
    )
    assert server.handle(
        {"method": "classify_catalog", "owner_epoch": "epoch-1", "clue": "bad"}
    ) == {"ok": False, "error": "catalog_model_quota_exhausted"}


def test_rust_worker_rejects_indexes_and_names_the_catalog_operation():
    source = (
        Path(__file__).resolve().parents[2]
        / "integrations/peri-bound/helixweave-bound.rs"
    ).read_text(encoding="utf-8")
    assert 'r["op"] == "classify"' in source
    assert "insufficient_evidence" in source
    assert "experiment.target" in source
    assert "is_ascii_digit" in source
    assert "Do not answer with a position or index" in source
