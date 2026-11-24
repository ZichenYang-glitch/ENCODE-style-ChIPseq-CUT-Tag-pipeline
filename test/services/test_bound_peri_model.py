from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
from dataclasses import replace
import sqlite3
import time

import pytest

from encode_pipeline.services.authoring_bridge import (
    BridgeDisconnected,
    BridgeTurnRequest,
)
from encode_pipeline.services.peri_runtime.bound_model import (
    BoundModel,
    Budget,
    binding_id,
    configured_model,
)


class Worker:
    binding = {"model": "bound-a"}

    def __init__(self):
        self.sends = []

    def __call__(self, r):
        if r["op"] == "binding":
            return self.binding
        self.sends.append(r)
        return {"binding": self.binding, "value": 7}


@pytest.fixture
def case(tmp_path):
    w = Worker()
    m = BoundModel(tmp_path / "ledger.sqlite", w)
    m.initialize(m.ledger)
    q = {
        "intent": "suggest",
        "field_tokens": ["config.threads"],
        "preferences": [{"field_token": "config.threads", "value": 3}],
    }
    r = BridgeTurnRequest(
        "sid",
        "turn",
        2,
        "instance",
        "epoch",
        "secret",
        q,
        ("secret",),
        "schema",
        "projection",
        1,
    )

    def grant(budget=Budget(1, 4096, 128, 1, 128, 100, 100, "CNY")):
        m.grant(
            sid=r.sid,
            epoch=r.owner_epoch,
            instance=r.instance_id,
            generation=r.generation,
            field="config.threads",
            binding=binding_id(w.binding),
            expires=int(time.time()) + 300,
            budget=budget,
        )

    return w, m, r, grant


@pytest.mark.parametrize(
    "change",
    [
        "no_grant",
        "field",
        "extra",
        "generation",
        "epoch",
        "sid",
        "binding",
        "authorization",
        "input",
        "output",
        "cost",
    ],
)
def test_refuses_before_send(case, change):
    w, m, r, grant = case
    budget = Budget(1, 4096, 128, 1, 128, 100, 100, "CNY")
    if change == "input":
        budget = replace(budget, input_tokens=4095)
    if change == "output":
        budget = replace(budget, output_tokens=127)
    if change == "cost":
        budget = replace(budget, input_rate=10000)
    if change != "no_grant":
        grant(budget)
    if change in {"sid", "epoch", "generation"}:
        key = {"epoch": "owner_epoch"}.get(change, change)
        r = replace(r, **{key: 3 if change == "generation" else "other"})
    if change == "extra":
        r.question["field_tokens"].append("private")
    if change == "field":
        r.question["field_tokens"][0] = "private"
        r.question["preferences"][0]["field_token"] = "private"
    if change == "binding":
        w.binding = {"model": "bound-b"}
    for turn_id in ("first-denial", "second-denial"):
        with pytest.raises(BridgeDisconnected):
            m.suggest(replace(r, turn_id=turn_id), lambda: change != "authorization")
        assert w.sends == []
        with closing(m.db()) as db:
            balances = db.execute(
                "SELECT requests,input_left,output_left,minor_left FROM grants"
            ).fetchall()
            expected = (
                []
                if change == "no_grant"
                else [
                    (
                        budget.requests,
                        budget.input_tokens,
                        budget.output_tokens,
                        budget.minor_units,
                    )
                ]
            )
            assert [tuple(row) for row in balances] == expected
            assert db.execute("SELECT COUNT(*) FROM reservations").fetchone()[0] == 0


def test_concurrent_reservation_and_adoption(case):
    import json

    w, m, r, grant = case
    grant()

    def attempt(i):
        try:
            return json.loads(m.suggest(replace(r, turn_id=str(i)), lambda: True))
        except BridgeDisconnected:
            return None

    with ThreadPoolExecutor(max_workers=2) as pool:
        results = list(pool.map(attempt, range(2)))
    assert len(w.sends) == 1
    assert w.sends[0] == {
        "op": "suggest",
        "binding": w.binding,
        "field": "config.threads",
        "value": 3,
        "output_tokens": 128,
    }
    turn = next(i for i, result in enumerate(results) if result)
    assert results[turn]["changes"] == [{"field_token": "config.threads", "value": 7}]
    args = (r.sid, r.owner_epoch, r.instance_id, r.generation, str(turn))
    m.check(*args)
    with closing(m.db()) as db, db:
        db.execute("UPDATE grants SET revoked=1")
    with pytest.raises(BridgeDisconnected):
        m.check(*args)


@pytest.mark.parametrize("problem", ["timeout", "rebound", "revoked", "invalid"])
def test_failed_delivery_never_refunds(case, problem):
    w, m, r, grant = case
    grant()

    def worker(request):
        response = w(request)
        if request["op"] == "suggest":
            if problem == "timeout":
                raise BridgeDisconnected
            if problem == "rebound":
                w.binding = {"model": "new"}
            if problem == "invalid":
                response["value"] = True
            if problem == "revoked":
                with closing(m.db()) as db, db:
                    db.execute("UPDATE grants SET revoked=1")
        return response

    m.worker = worker
    with pytest.raises(BridgeDisconnected):
        m.suggest(r, lambda: True)
    with closing(m.db()) as db:
        assert tuple(
            db.execute(
                "SELECT requests,input_left,output_left,minor_left FROM grants"
            ).fetchone()
        ) == (0, 0, 0, 0)
        assert db.execute("SELECT delivered FROM reservations").fetchone()[0] == 0
    assert len(w.sends) == 1


def test_incomplete_config_is_never_silent_fallback():
    assert configured_model({}) is None
    with pytest.raises(ValueError):
        configured_model({"HELIXWEAVE_PERI_MODEL_LEDGER": "/x"})


def test_immutable_scope(case):
    _, _, _, grant = case
    grant()
    with pytest.raises(sqlite3.IntegrityError):
        grant()


def test_malformed_ledger_no_network(case):
    w, model, req, grant = case
    grant()
    with closing(model.db()) as db, db:
        db.execute("UPDATE grants SET input_left='unknown'")
    with pytest.raises(BridgeDisconnected):
        model.suggest(req, lambda: True)
    assert w.sends == []


def test_expired_delivery_cannot_be_adopted(case):
    _, model, req, grant = case
    grant()
    model.suggest(req, lambda: True)
    with closing(model.db()) as db, db:
        db.execute("UPDATE reservations SET adopt_until=0")
    with pytest.raises(BridgeDisconnected):
        model.check(
            req.sid, req.owner_epoch, req.instance_id, req.generation, req.turn_id
        )
