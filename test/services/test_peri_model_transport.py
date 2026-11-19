"""Offline P2 transport: suggestion source, allow-list, timeout, and revocation."""

from __future__ import annotations

import ast
from collections import deque
from collections.abc import Callable
from dataclasses import fields
from pathlib import Path

import pytest

from encode_pipeline.services.peri_runtime.model_transport import (
    CONTRACT_VERSION,
    QUALIFICATION_STUB_ROLE,
    TEST_SOURCE,
    TEST_TRANSPORT_ROLE,
    VENDOR_TRANSPORT_ROLE,
    Delivery,
    ModelResponse,
    Refusal,
    Suggestion,
    SyntheticField,
    SyntheticSuggestionTransport,
    TransportProjection,
    TransportRejected,
    TransportRequest,
    UsageRecord,
    VendorUnselected,
    open_vendor_model_transport,
    suggestion_changes,
)
from encode_pipeline.services.peri_runtime.stub import (
    ModelStub,
    open_qualification_model_transport,
)


ROOT = Path(__file__).resolve().parents[2]
Responder = Callable[[TransportProjection], ModelResponse]


def _field(token: str, value: str) -> SyntheticField:
    return SyntheticField(token, value)


def _request(
    *,
    allowed: tuple[str, ...] = ("syn.alpha",),
    projected: tuple[SyntheticField, ...] = (_field("syn.alpha", "from_request"),),
    call_id: str = "call_aaaaaaaa",
) -> TransportRequest:
    return TransportRequest(call_id=call_id, allowed=allowed, projected=projected)


def _response(
    value: str = "from_model",
    *,
    token: str = "syn.alpha",
    response_id: str = "resp_aaaaaaaa",
    fields: tuple[SyntheticField, ...] | None = None,
) -> ModelResponse:
    body = fields if fields is not None else (_field(token, value),)
    return ModelResponse(response_id, body)


def _transport(
    responder: Responder,
    *,
    authorize: Callable[[], str | None] | None = None,
    clock: Callable[[], float] | None = None,
    timeout_seconds: float = 1,
    cancelled: Callable[[], bool] | None = None,
) -> SyntheticSuggestionTransport:
    return SyntheticSuggestionTransport(
        authorize=authorize or (lambda: "grant_aaaaaaaa"),
        clock=clock or (lambda: 0.0),
        timeout_seconds=timeout_seconds,
        responder=responder,
        cancelled=cancelled,
    )


def _clock(*values: float) -> Callable[[], float]:
    pending = iter(values)
    last = values[-1]

    def clock() -> float:
        nonlocal last
        try:
            last = next(pending)
        except StopIteration:
            pass
        return last

    return clock


def _deliver(
    transport: SyntheticSuggestionTransport, request: TransportRequest
) -> Delivery:
    outcome = transport.deliver(request)
    assert isinstance(outcome, Delivery)
    return outcome


def test_roles_stay_distinct() -> None:
    assert QUALIFICATION_STUB_ROLE == "qualification_stub"
    assert TEST_TRANSPORT_ROLE == "test_transport"
    assert VENDOR_TRANSPORT_ROLE == "vendor_transport"
    assert (
        len({QUALIFICATION_STUB_ROLE, TEST_TRANSPORT_ROLE, VENDOR_TRANSPORT_ROLE}) == 3
    )
    assert SyntheticSuggestionTransport.role == TEST_TRANSPORT_ROLE
    assert SyntheticSuggestionTransport.contract_version == CONTRACT_VERSION
    assert not hasattr(ModelStub, "deliver")
    assert not hasattr(ModelStub, "adopt")
    assert open_qualification_model_transport is not open_vendor_model_transport
    with pytest.raises(VendorUnselected) as caught:
        open_vendor_model_transport()
    assert caught.value.reason == "vendor_unselected"
    assert str(caught.value) == "vendor_unselected"


def test_response_values_become_the_suggestion() -> None:
    queued = deque(
        [
            _response("from_first", response_id="resp_aaaaaaaa"),
            _response("from_second", response_id="resp_bbbbbbbb"),
        ]
    )
    seen: list[TransportProjection] = []

    def responder(projection: TransportProjection) -> ModelResponse:
        seen.append(projection)
        return queued.popleft()

    transport = _transport(responder)
    request = _request()
    first_delivery = _deliver(transport, request)
    second_delivery = _deliver(transport, request)
    first = transport.adopt(first_delivery)
    second = transport.adopt(second_delivery)

    assert first.fields == (_field("syn.alpha", "from_first"),)
    assert second.fields == (_field("syn.alpha", "from_second"),)
    assert first != second
    assert first.response_id == "resp_aaaaaaaa"
    assert second.response_id == "resp_bbbbbbbb"
    assert first.source == TEST_SOURCE
    assert suggestion_changes(first) == (
        {"field_token": "syn.alpha", "value": "from_first"},
    )
    assert suggestion_changes(second) == (
        {"field_token": "syn.alpha", "value": "from_second"},
    )
    assert all(item.fields == request.projected for item in seen)
    assert request.projected[0].value not in {
        first.fields[0].value,
        second.fields[0].value,
    }
    assert "from_first" not in repr(first_delivery)
    with pytest.raises(TransportRejected) as repeated:
        transport.adopt(first_delivery)
    assert repeated.value.reason == "not_delivered"
    assert transport.usage() == UsageRecord(
        attempted=2, responder_calls=2, delivered=2, adopted=2
    )


def test_allowed_fields_bound_projection_and_adoption() -> None:
    seen: list[TransportProjection] = []

    def responder(projection: TransportProjection) -> ModelResponse:
        seen.append(projection)
        return _response(
            response_id="resp_cccccccc",
            fields=(_field("syn.beta", "from_response"),),
        )

    transport = _transport(responder)
    request = _request(
        allowed=("syn.alpha", "syn.beta"),
        projected=(_field("syn.alpha", "from_request"),),
    )
    suggestion = transport.adopt(_deliver(transport, request))
    assert seen[0].fields == request.projected
    assert {item.name for item in fields(TransportProjection)} == {"call_id", "fields"}
    assert suggestion.fields == (_field("syn.beta", "from_response"),)
    assert suggestion_changes(suggestion) == (
        {"field_token": "syn.beta", "value": "from_response"},
    )

    blocked = {"called": False}

    def blocked_responder(projection: TransportProjection) -> ModelResponse:
        del projection
        blocked["called"] = True
        return _response("unused", response_id="resp_dddddddd")

    blocked_transport = _transport(blocked_responder)
    hidden = blocked_transport.deliver(
        _request(projected=(_field("syn.hidden", "hidden_value"),))
    )
    assert isinstance(hidden, Refusal)
    assert hidden.reason == "field_not_allowed"
    assert hidden.delivered is False
    assert hidden.adopted is False
    assert blocked["called"] is False
    assert "hidden_value" not in repr(hidden)

    def mixed_responder(projection: TransportProjection) -> ModelResponse:
        del projection
        return _response(
            response_id="resp_dddddddd",
            fields=(
                _field("syn.alpha", "kept_value"),
                _field("syn.other", "blocked_value"),
            ),
        )

    mixed_transport = _transport(mixed_responder)
    mixed = mixed_transport.deliver(
        _request(
            call_id="call_bbbbbbbb",
            allowed=("syn.alpha",),
            projected=(_field("syn.alpha", "from_request"),),
        )
    )
    assert isinstance(mixed, Refusal)
    assert mixed.reason == "field_not_allowed"
    assert mixed.delivered is False
    assert "kept_value" not in repr(mixed)
    assert "blocked_value" not in repr(mixed)
    assert mixed_transport.usage().responder_calls == 1
    assert mixed_transport.usage().delivered == 0
    assert transport.usage().delivered == 1


def test_timeout_discards_the_response() -> None:
    def responder(projection: TransportProjection) -> ModelResponse:
        del projection
        return _response("late_value", response_id="resp_eeeeeeee")

    transport = _transport(responder, clock=_clock(0.0, 5.0), timeout_seconds=1)
    outcome = transport.deliver(_request())
    assert isinstance(outcome, Refusal)
    assert outcome.reason == "timeout"
    assert outcome.delivered is False
    assert "late_value" not in repr(outcome)
    usage = transport.usage()
    assert usage.responder_calls == 1
    assert usage.delivered == 0
    assert usage.adopted == 0
    forged = Delivery("missing", "call_aaaaaaaa", "resp_eeeeeeee", TEST_SOURCE)
    with pytest.raises(TransportRejected) as caught:
        transport.adopt(forged)
    assert caught.value.reason == "not_delivered"
    assert "late_value" not in str(caught.value)

    on_time = _transport(responder, clock=_clock(0.0, 1.0), timeout_seconds=1)
    assert isinstance(on_time.deliver(_request()), Delivery)


def test_revocation_before_send_does_not_call_the_responder() -> None:
    def responder(projection: TransportProjection) -> ModelResponse:
        del projection
        raise AssertionError("responder ran after revocation")

    transport = _transport(responder, authorize=lambda: None)
    outcome = transport.deliver(_request())
    assert isinstance(outcome, Refusal)
    assert outcome.reason == "authorization_revoked"
    assert transport.usage().responder_calls == 0
    assert transport.usage().delivered == 0


def test_revocation_after_the_response_does_not_deliver() -> None:
    checks = {"count": 0}

    def authorize() -> str | None:
        checks["count"] += 1
        if checks["count"] == 1:
            return "grant_aaaaaaaa"
        return None

    def responder(projection: TransportProjection) -> ModelResponse:
        del projection
        return _response("revoked_value", response_id="resp_ffffffff")

    transport = _transport(responder, authorize=authorize)
    outcome = transport.deliver(_request())
    assert isinstance(outcome, Refusal)
    assert outcome.reason == "authorization_revoked"
    assert checks["count"] == 2
    assert "revoked_value" not in repr(outcome)
    assert transport.usage().responder_calls == 1
    assert transport.usage().delivered == 0
    assert transport.usage().adopted == 0


def test_revocation_before_adopt_drops_the_suggestion() -> None:
    grant = {"period": "grant_aaaaaaaa"}

    def responder(projection: TransportProjection) -> ModelResponse:
        del projection
        return _response("adopt_me", response_id="resp_11111111")

    def authorize() -> str | None:
        return grant["period"]

    transport = _transport(responder, authorize=authorize)
    delivery = _deliver(transport, _request())
    assert "adopt_me" not in repr(delivery)
    grant["period"] = None
    with pytest.raises(TransportRejected) as caught:
        transport.adopt(delivery)
    assert caught.value.reason == "authorization_revoked"
    assert "adopt_me" not in str(caught.value)
    with pytest.raises(TransportRejected) as repeated:
        transport.adopt(delivery)
    assert repeated.value.reason == "not_delivered"
    assert transport.usage().adopted == 0


def test_expired_delivery_cannot_be_adopted() -> None:
    now = {"t": 0.0}

    def clock() -> float:
        return now["t"]

    def responder(projection: TransportProjection) -> ModelResponse:
        del projection
        return _response("stale_value", response_id="resp_33333333")

    transport = _transport(responder, clock=clock, timeout_seconds=1)
    delivery = _deliver(transport, _request())
    now["t"] = 2.0
    with pytest.raises(TransportRejected) as caught:
        transport.adopt(delivery)
    assert caught.value.reason == "timeout"
    assert "stale_value" not in str(caught.value)
    with pytest.raises(TransportRejected) as repeated:
        transport.adopt(delivery)
    assert repeated.value.reason == "not_delivered"
    assert transport.usage().delivered == 1
    assert transport.usage().adopted == 0


def test_adopt_at_the_deadline_succeeds() -> None:
    now = {"t": 0.0}

    def clock() -> float:
        return now["t"]

    def responder(projection: TransportProjection) -> ModelResponse:
        del projection
        return _response("on_time", response_id="resp_44444444")

    transport = _transport(responder, clock=clock, timeout_seconds=1)
    delivery = _deliver(transport, _request())
    now["t"] = 1.0
    suggestion = transport.adopt(delivery)
    assert suggestion.fields == (_field("syn.alpha", "on_time"),)
    assert transport.usage().adopted == 1


def test_reauthorization_adopts_only_a_new_delivery() -> None:
    grant = {"period": "grant_aaaaaaaa"}
    queued = deque(
        [
            _response("old_value", response_id="resp_55555555"),
            _response("new_value", response_id="resp_66666666"),
        ]
    )

    def authorize() -> str | None:
        return grant["period"]

    def responder(projection: TransportProjection) -> ModelResponse:
        del projection
        return queued.popleft()

    transport = _transport(responder, authorize=authorize)
    old = _deliver(transport, _request())
    grant["period"] = "grant_bbbbbbbb"
    with pytest.raises(TransportRejected) as caught:
        transport.adopt(old)
    assert caught.value.reason == "authorization_period_mismatch"
    assert "old_value" not in str(caught.value)
    assert "grant_bbbbbbbb" not in str(caught.value)
    with pytest.raises(TransportRejected) as repeated:
        transport.adopt(old)
    assert repeated.value.reason == "not_delivered"
    assert transport.usage().adopted == 0
    fresh = _deliver(transport, _request(call_id="call_bbbbbbbb"))
    suggestion = transport.adopt(fresh)
    assert suggestion.fields == (_field("syn.alpha", "new_value"),)
    assert transport.usage().adopted == 1
    with pytest.raises(TransportRejected) as still_old:
        transport.adopt(old)
    assert still_old.value.reason == "not_delivered"


def test_cancellation_stops_delivery_and_adopt() -> None:
    def responder(projection: TransportProjection) -> ModelResponse:
        del projection
        raise AssertionError("responder ran after cancel")

    stopped = _transport(responder, cancelled=lambda: True)
    outcome = stopped.deliver(_request())
    assert isinstance(outcome, Refusal)
    assert outcome.reason == "cancelled"
    assert stopped.usage().responder_calls == 0

    gate = {"stop": False}

    def live(projection: TransportProjection) -> ModelResponse:
        del projection
        return _response("cancel_me", response_id="resp_22222222")

    transport = _transport(live, cancelled=lambda: gate["stop"])
    delivery = _deliver(transport, _request())
    gate["stop"] = True
    with pytest.raises(TransportRejected) as caught:
        transport.adopt(delivery)
    assert caught.value.reason == "cancelled"
    assert "cancel_me" not in str(caught.value)
    assert transport.usage().adopted == 0


def test_failed_callbacks_do_not_leak_or_deliver() -> None:
    def explode_auth() -> str | None:
        raise RuntimeError("token=secret-marker")

    def responder(projection: TransportProjection) -> ModelResponse:
        del projection
        return _response("should_not_appear")

    refused_transport = _transport(responder, authorize=explode_auth)
    refused = refused_transport.deliver(_request())
    assert isinstance(refused, Refusal)
    assert refused.reason == "authorization_failed"
    assert "secret-marker" not in repr(refused)
    assert refused_transport.usage().responder_calls == 0

    def explode_responder(projection: TransportProjection) -> ModelResponse:
        del projection
        raise RuntimeError("body=late-secret")

    failed_transport = _transport(explode_responder)
    failed = failed_transport.deliver(_request())
    assert isinstance(failed, Refusal)
    assert failed.reason == "responder_failed"
    assert "late-secret" not in repr(failed)
    assert failed_transport.usage().responder_calls == 1
    assert failed_transport.usage().delivered == 0


def test_vendor_opener_stays_offline(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("OPENAI_API_KEY", "not-a-real-key")
    monkeypatch.setenv("HELIXWEAVE_MODEL_ENDPOINT", "http://127.0.0.1/secret")
    with pytest.raises(VendorUnselected) as caught:
        open_vendor_model_transport(
            endpoint="http://example.invalid",
            api_key="sk-test",
        )
    assert caught.value.reason == "vendor_unselected"
    text = str(caught.value)
    assert "not-a-real-key" not in text
    assert "sk-test" not in text
    assert "example.invalid" not in text
    assert {item.name for item in fields(UsageRecord)} == {
        "attempted",
        "responder_calls",
        "delivered",
        "adopted",
    }


def test_transport_module_does_not_import_authoring_or_network() -> None:
    source = (
        ROOT / "src/encode_pipeline/services/peri_runtime/model_transport.py"
    ).read_text(encoding="utf-8")
    modules: set[str] = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.ImportFrom) and node.module is not None:
            modules.add(node.module)
        elif isinstance(node, ast.Import):
            modules.update(alias.name for alias in node.names)
    forbidden = {
        "encode_pipeline.services.authoring",
        "encode_pipeline.persistence",
        "encode_pipeline.persistence.authoring",
        "encode_pipeline.api",
        "encode_pipeline.services.peri_runtime.stub",
        "socket",
        "http",
        "http.client",
        "httpx",
        "urllib",
        "urllib.request",
        "subprocess",
        "openai",
        "os",
    }
    assert modules.isdisjoint(forbidden)
    runtime = (ROOT / "src/encode_pipeline/services/peri_runtime/runtime.py").read_text(
        encoding="utf-8"
    )
    session = (ROOT / "src/encode_pipeline/services/peri_runtime/session.py").read_text(
        encoding="utf-8"
    )
    assert "open_qualification_model_transport" in runtime
    assert "open_qualification_model_transport" in session
    for name in ("SyntheticSuggestionTransport", "open_vendor_model_transport"):
        assert name not in runtime
        assert name not in session


def test_malformed_payloads_are_refusals() -> None:
    def bad_response(projection: TransportProjection) -> ModelResponse:
        del projection
        return "not-a-response"  # type: ignore[return-value]

    outcome = _transport(bad_response).deliver(_request())
    assert isinstance(outcome, Refusal)
    assert outcome.reason == "malformed_response"

    def unused(projection: TransportProjection) -> ModelResponse:
        del projection
        raise AssertionError("responder ran")

    bad_request = _transport(unused).deliver(
        TransportRequest(
            "call_aaaaaaaa",
            ("syn.alpha",),
            ("not-a-field",),  # type: ignore[arg-type]
        )
    )
    assert isinstance(bad_request, Refusal)
    assert bad_request.reason == "malformed_request"
    assert "not-a-field" not in repr(bad_request)


def test_suggestion_changes_reject_a_hand_built_payload() -> None:
    forged = Suggestion(
        contract_version=CONTRACT_VERSION,
        call_id="call_aaaaaaaa",
        response_id="resp_aaaaaaaa",
        fields=(_field("syn.alpha", "donor/secret"),),
        source=TEST_SOURCE,
    )
    with pytest.raises(TransportRejected) as caught:
        suggestion_changes(forged)
    assert caught.value.reason == "malformed_suggestion"
    assert "donor" not in str(caught.value)
