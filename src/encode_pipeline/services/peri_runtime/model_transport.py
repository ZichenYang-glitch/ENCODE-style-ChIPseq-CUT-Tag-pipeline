"""Offline suggestion transport for the P2 boundary.

Three roles stay separate. The qualification stub remains ``ModelStub`` in
``stub.py`` and only qualifies the pinned Peri binary. This module is the
in-process test transport: synthetic fields, no socket, no credentials, and
no experiment records. ``open_vendor_model_transport`` refuses a real
provider. No vendor, outbound field list, or spend limit is selected here.

``adopt`` releases a suggestion this transport still holds. It does not write
a draft and it does not call the authoring service.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
import math
import re
import threading
from typing import NoReturn
import uuid


CONTRACT_VERSION = "p2-transport/1"
QUALIFICATION_STUB_ROLE = "qualification_stub"
TEST_TRANSPORT_ROLE = "test_transport"
VENDOR_TRANSPORT_ROLE = "vendor_transport"
TEST_SOURCE = "test_transport"
_FIELD_LIMIT = 8

_TOKEN = re.compile(r"syn\.[a-z][a-z0-9_]{0,31}\Z")
_VALUE = re.compile(r"[A-Za-z0-9_.:-]{1,64}\Z")
_CALL_ID = re.compile(r"call_[a-z0-9]{8,32}\Z")
_RESPONSE_ID = re.compile(r"resp_[a-z0-9]{8,32}\Z")

_DELIVER_REASONS = frozenset(
    {
        "malformed_request",
        "field_not_allowed",
        "authorization_revoked",
        "authorization_failed",
        "cancelled",
        "timeout",
        "malformed_clock",
        "responder_failed",
        "malformed_response",
    }
)
_EXCEPTION_REASONS = _DELIVER_REASONS | {
    "unknown_delivery",
    "not_delivered",
    "unknown_source",
    "timeout_not_configured",
    "vendor_unselected",
    "malformed_suggestion",
}

Authorize = Callable[[], bool]
Clock = Callable[[], float]
Cancel = Callable[[], bool]


class TransportRejected(Exception):
    """Stable reason code. The message is that code and nothing else."""

    def __init__(self, reason: str) -> None:
        if reason not in _EXCEPTION_REASONS:
            reason = "malformed_request"
        self.reason = reason
        super().__init__(reason)


class VendorUnselected(TransportRejected):
    """A real provider was requested before one was selected."""

    def __init__(self) -> None:
        super().__init__("vendor_unselected")


@dataclass(frozen=True)
class SyntheticField:
    """One fixture field. The token grammar is not a catalog field name."""

    token: str
    value: str


@dataclass(frozen=True)
class TransportRequest:
    """Allowed tokens plus the projection actually handed to the responder."""

    call_id: str
    allowed: tuple[str, ...]
    projected: tuple[SyntheticField, ...]


@dataclass(frozen=True)
class TransportProjection:
    """The only object the responder is allowed to see."""

    call_id: str
    fields: tuple[SyntheticField, ...]


@dataclass(frozen=True)
class ModelResponse:
    """Scripted responder output. This is not a provider payload."""

    response_id: str
    fields: tuple[SyntheticField, ...]


@dataclass(frozen=True)
class Delivery:
    """Handle for one accepted exchange. Field values stay on the transport."""

    delivery_id: str
    call_id: str
    response_id: str
    source: str


@dataclass(frozen=True)
class Refusal:
    """The exchange produced nothing that can be adopted."""

    reason: str
    delivered: bool = False
    adopted: bool = False

    def __post_init__(self) -> None:
        if self.delivered or self.adopted or self.reason not in _DELIVER_REASONS:
            raise TransportRejected("malformed_request")


@dataclass(frozen=True)
class Suggestion:
    """A suggestion whose fields are the accepted response fields."""

    contract_version: str
    call_id: str
    response_id: str
    fields: tuple[SyntheticField, ...]
    source: str


@dataclass(frozen=True)
class UsageRecord:
    """Call counts for this process. Not a provider quota or a spend limit."""

    attempted: int
    responder_calls: int
    delivered: int
    adopted: int


@dataclass(frozen=True)
class _Pending:
    call_id: str
    response_id: str
    fields: tuple[SyntheticField, ...]


Responder = Callable[[TransportProjection], ModelResponse]


def open_vendor_model_transport(*_args: object, **_kwargs: object) -> NoReturn:
    """Refuse a real provider.

    Arguments are ignored. This function does not read the environment, open
    a socket, or construct a client.
    """

    raise VendorUnselected


def suggestion_changes(suggestion: Suggestion) -> tuple[dict[str, str], ...]:
    """Return the authoring change shape for one adopted test suggestion.

    The next wiring step can hand this tuple to the authoring owner. This
    function does not import or call that service.
    """

    if type(suggestion) is not Suggestion or suggestion.source != TEST_SOURCE:
        raise TransportRejected("unknown_source")
    if suggestion.contract_version != CONTRACT_VERSION:
        raise TransportRejected("malformed_suggestion")
    changes: list[dict[str, str]] = []
    seen: set[str] = set()
    for field in suggestion.fields:
        if (
            type(field) is not SyntheticField
            or not _token_ok(field.token)
            or not _value_ok(field.value)
            or field.token in seen
        ):
            raise TransportRejected("malformed_suggestion")
        seen.add(field.token)
        changes.append({"field_token": field.token, "value": field.value})
    if not changes:
        raise TransportRejected("malformed_suggestion")
    return tuple(changes)


class SyntheticSuggestionTransport:
    """In-process transport. Callbacks must not re-enter the same instance."""

    role = TEST_TRANSPORT_ROLE
    contract_version = CONTRACT_VERSION

    def __init__(
        self,
        *,
        authorize: Authorize,
        clock: Clock,
        timeout_seconds: float,
        responder: Responder,
        cancelled: Cancel | None = None,
    ) -> None:
        if not _timeout_ok(timeout_seconds):
            raise TransportRejected("timeout_not_configured")
        self._authorize = authorize
        self._clock = clock
        self._timeout_seconds = float(timeout_seconds)
        self._responder = responder
        self._cancel = cancelled
        self._lock = threading.Lock()
        self._pending: dict[str, _Pending] = {}
        self._attempted = 0
        self._responder_calls = 0
        self._delivered = 0
        self._adopted = 0

    def deliver(self, request: TransportRequest) -> Delivery | Refusal:
        """Project allowed fields, then keep a response only if it still qualifies."""

        with self._lock:
            self._attempted += 1
        refusal = _request_refusal(request)
        if refusal is not None:
            return refusal
        reason = _gate_reason(self._authorization_reason(), self._cancellation_reason())
        if reason is not None:
            return Refusal(reason)
        try:
            started = self._clock()
        except Exception:
            return Refusal("malformed_clock")
        if not _finite(started):
            return Refusal("malformed_clock")
        try:
            response = self._responder(
                TransportProjection(call_id=request.call_id, fields=request.projected)
            )
        except Exception:
            self._count_responder()
            return Refusal("responder_failed")
        self._count_responder()
        try:
            finished = self._clock()
        except Exception:
            return Refusal("malformed_clock")
        if not _finite(finished) or finished < started:
            return Refusal("malformed_clock")
        if finished - started > self._timeout_seconds:
            return Refusal("timeout")
        reason = _gate_reason(self._authorization_reason(), self._cancellation_reason())
        if reason is not None:
            return Refusal(reason)
        if not isinstance(response, ModelResponse):
            return Refusal("malformed_response")
        response_refusal = _response_refusal(response, request.allowed)
        if response_refusal is not None:
            return response_refusal
        delivery_id = uuid.uuid4().hex
        pending = _Pending(
            call_id=request.call_id,
            response_id=response.response_id,
            fields=tuple(response.fields),
        )
        with self._lock:
            self._pending[delivery_id] = pending
            self._delivered += 1
        return Delivery(
            delivery_id=delivery_id,
            call_id=request.call_id,
            response_id=response.response_id,
            source=TEST_SOURCE,
        )

    def adopt(self, delivery: Delivery) -> Suggestion:
        """Return the stored suggestion, or refuse when the grant is gone."""

        if type(delivery) is not Delivery or delivery.source != TEST_SOURCE:
            raise TransportRejected("unknown_delivery")
        with self._lock:
            pending = self._pending.get(delivery.delivery_id)
        if pending is None:
            raise TransportRejected("not_delivered")
        reason = _gate_reason(self._authorization_reason(), self._cancellation_reason())
        if reason is not None:
            self._drop(delivery.delivery_id)
            raise TransportRejected(reason)
        with self._lock:
            pending = self._pending.pop(delivery.delivery_id, None)
            if pending is None:
                raise TransportRejected("not_delivered")
            self._adopted += 1
        return Suggestion(
            contract_version=CONTRACT_VERSION,
            call_id=pending.call_id,
            response_id=pending.response_id,
            fields=pending.fields,
            source=TEST_SOURCE,
        )

    def usage(self) -> UsageRecord:
        with self._lock:
            return UsageRecord(
                attempted=self._attempted,
                responder_calls=self._responder_calls,
                delivered=self._delivered,
                adopted=self._adopted,
            )

    def _authorization_reason(self) -> str | None:
        try:
            allowed = bool(self._authorize())
        except Exception:
            return "authorization_failed"
        if not allowed:
            return "authorization_revoked"
        return None

    def _cancellation_reason(self) -> str | None:
        if self._cancel is None:
            return None
        try:
            cancelled = bool(self._cancel())
        except Exception:
            return "cancelled"
        if cancelled:
            return "cancelled"
        return None

    def _count_responder(self) -> None:
        with self._lock:
            self._responder_calls += 1

    def _drop(self, delivery_id: str) -> None:
        with self._lock:
            self._pending.pop(delivery_id, None)


def _gate_reason(authorization: str | None, cancellation: str | None) -> str | None:
    if authorization is not None:
        return authorization
    return cancellation


def _timeout_ok(value: object) -> bool:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return False
    return _finite(float(value)) and float(value) > 0


def _finite(value: float) -> bool:
    return math.isfinite(value)


def _token_ok(value: object) -> bool:
    return isinstance(value, str) and _TOKEN.fullmatch(value) is not None


def _value_ok(value: object) -> bool:
    return isinstance(value, str) and _VALUE.fullmatch(value) is not None


def _fields_well_typed(fields: tuple[object, ...]) -> bool:
    return all(type(field) is SyntheticField for field in fields)


def _unique_tokens(fields: tuple[SyntheticField, ...]) -> bool:
    tokens = [field.token for field in fields]
    return len(tokens) == len(set(tokens))


def _request_refusal(request: TransportRequest) -> Refusal | None:
    if (
        type(request) is not TransportRequest
        or not isinstance(request.call_id, str)
        or _CALL_ID.fullmatch(request.call_id) is None
        or not isinstance(request.allowed, tuple)
        or not request.allowed
        or len(request.allowed) > _FIELD_LIMIT
        or len(set(request.allowed)) != len(request.allowed)
        or not isinstance(request.projected, tuple)
        or len(request.projected) > _FIELD_LIMIT
        or not _fields_well_typed(request.projected)
    ):
        return Refusal("malformed_request")
    allowed = set(request.allowed)
    if any(not _token_ok(token) for token in request.allowed):
        return Refusal("malformed_request")
    if not _unique_tokens(request.projected):
        return Refusal("malformed_request")
    for field in request.projected:
        if not _value_ok(field.value):
            return Refusal("malformed_request")
        if not _token_ok(field.token) or field.token not in allowed:
            return Refusal("field_not_allowed")
    return None


def _response_refusal(
    response: ModelResponse, allowed: tuple[str, ...]
) -> Refusal | None:
    if (
        not isinstance(response.response_id, str)
        or _RESPONSE_ID.fullmatch(response.response_id) is None
        or not isinstance(response.fields, tuple)
        or not _fields_well_typed(response.fields)
        or not response.fields
        or len(response.fields) > _FIELD_LIMIT
        or not _unique_tokens(response.fields)
    ):
        return Refusal("malformed_response")
    permitted = set(allowed)
    for field in response.fields:
        if not _token_ok(field.token):
            return Refusal("malformed_response")
        if field.token not in permitted:
            return Refusal("field_not_allowed")
        if not _value_ok(field.value):
            return Refusal("malformed_response")
    return None
