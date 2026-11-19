"""Loopback-only P2 acceptance through the existing authoring HTTP routes."""

from __future__ import annotations

from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
import json
import os
from pathlib import Path
from secrets import token_hex
import struct
from tempfile import TemporaryDirectory
from threading import Event, Lock
from time import monotonic

import pytest

from api_test_client import ApiTestClient
from encode_pipeline.services.peri_runtime.identity import binary_matches
from encode_pipeline.services.peri_runtime.model_transport import (
    ModelResponse,
    Refusal,
    SyntheticField,
    SyntheticSuggestionTransport,
    TransportProjection,
    TransportRejected,
    TransportRequest,
)
from encode_pipeline.services.peri_runtime.owner import NativeAuthoringBridge
from encode_pipeline.services.peri_runtime.policy import (
    PORT,
    SYNTHETIC_API_KEY,
    settings,
)
from encode_pipeline.services.peri_runtime.session import NativeRejected, NativeSession
from encode_pipeline.services.peri_runtime.stub import ModelStub


API = "/api/v1/agent/authoring-sessions"
PROJECTED_FIELD = "syn.threads"
# Requests, input tokens, output tokens, and minor currency units. Fixture only.
RESERVATION = (1, 512, 64, 0)


@dataclass
class Endpoint:
    name: str
    value: str
    stub: ModelStub = field(init=False)
    paused: bool = False
    received: Event = field(default_factory=Event)
    release: Event = field(default_factory=Event)

    @property
    def route(self) -> str:
        return f"/fixture-binding-{self.name}/v1"

    @property
    def key(self) -> str:
        return f"fixture-binding-key-{self.name}"

    def config(self) -> bytes:
        document = settings()
        provider = document["config"]["providers"][0]
        provider["baseUrl"] = f"http://127.0.0.1:{PORT}{self.route}"
        provider["apiKey"] = self.key
        if self.name == "missing":
            # Leave the other configured aliases usable to detect fallback.
            document["config"]["active_alias"] = "missing"
        payload = json.dumps(document).encode()
        return struct.pack(">I", len(payload)) + payload

    def open_stub(self, path: Path) -> ModelStub:
        self.stub = ModelStub(path)
        endpoint = self

        class Handler(self.stub.server.RequestHandlerClass):
            def do_POST(self):
                if self.path != endpoint.route + "/chat/completions":
                    self.send_error(404)
                    return
                if self.headers.get("Authorization") != f"Bearer {endpoint.key}":
                    self.send_error(403)
                    return
                # Reuse the qualified stub after checking this endpoint's route and key.
                self.path = "/v1/chat/completions"
                self.headers.replace_header(
                    "Authorization", f"Bearer {SYNTHETIC_API_KEY}"
                )
                super().do_POST()

            def send_response(self, code, message=None):
                if code == 200 and endpoint.paused:
                    endpoint.received.set()
                    endpoint.release.wait(1.5)  # Inside the relay idle deadline.
                super().send_response(code, message)

        self.stub.server.RequestHandlerClass = Handler
        return self.stub

    def answer_once(self) -> None:
        self.stub.script(self.value)


@dataclass
class LocalGrant:
    allowed: frozenset[str] = frozenset({PROJECTED_FIELD})
    remaining: tuple[int, ...] | None = RESERVATION
    lock: Lock = field(default_factory=Lock)

    def reserve(self, request: TransportRequest) -> None:
        if (
            not set(request.allowed).union(f.token for f in request.projected)
            <= self.allowed
        ):
            raise NativeRejected("outbound_field_not_granted")
        with self.lock:
            if self.remaining is None:
                raise NativeRejected("budget_unproven")
            amounts = tuple(
                left - needed
                for left, needed in zip(self.remaining, RESERVATION, strict=True)
            )
            if min(amounts) < 0:
                raise NativeRejected("budget_exhausted")
            self.remaining = amounts


class LocalPeri:
    """Test operator controls injected Peri settings and its private endpoint."""

    def __init__(self, endpoint: Endpoint) -> None:
        self.activate(endpoint)
        self.session: LocalSession
        self.receipts: list[dict[str, object]] = []
        self.clock: Callable[[], float] = monotonic
        self.next_grant = LocalGrant()

    def activate(self, endpoint: Endpoint) -> None:
        self._binding = (endpoint,)

    @property
    def selected(self) -> Endpoint:
        return self._binding[0]

    def open_session(self, *args, **kwargs):
        binding = self._binding
        (target,) = binding
        self.session = LocalSession(
            *args,
            **kwargs,
            model_transport=target.open_stub,
            selected=lambda: self._binding is binding,
            grant=self.next_grant,
            clock=self.clock,
            receipts=self.receipts,
        )
        self.session.bootstrap_payload = target.config()
        return self.session


class LocalSession(NativeSession):
    """Exercise the real bridge and native prompt with the synthetic P2 grammar."""

    def __init__(
        self,
        *args,
        selected: Callable[[], bool],
        grant: LocalGrant,
        clock: Callable[[], float],
        receipts: list[dict[str, object]],
        **kwargs,
    ) -> None:
        self._selected = selected
        self._grant = grant
        self._clock = clock
        self._receipts = receipts
        self._period = "grant_" + token_hex(8)
        self.preflight_failure: str | None = None
        self.delivery_refusal: str | None = None
        super().__init__(*args, **kwargs)

    def rpc(self, method, params, *, timeout=15):
        if method == "initialize":
            params["clientCapabilities"]["_meta"]["peri.sessionWorkspaceV1"] = True
        return super().rpc(method, params, timeout=timeout)

    def run_authoring_turn(self, final_text: str) -> str:
        try:
            return self._local_turn(final_text)
        except TransportRejected as error:
            raise NativeRejected(str(error)) from None

    def _egress_allowed(self) -> bool:
        return self._selected() and super()._egress_allowed()

    def _local_turn(self, final_text: str) -> str:
        proposal = json.loads(final_text)
        changes = proposal["changes"]
        if len(changes) != 1 or changes[0]["field_token"] != "config.threads":
            raise NativeRejected("outbound_field_not_granted")
        request = TransportRequest(
            call_id="call_" + token_hex(8),
            allowed=(PROJECTED_FIELD,),
            projected=(SyntheticField(PROJECTED_FIELD, str(changes[0]["value"])),),
        )

        def respond(projection: TransportProjection) -> ModelResponse:
            try:
                self._grant.reserve(request)
            except NativeRejected as error:
                self.preflight_failure = str(error)
                raise
            text = self.prompt(
                json.dumps(
                    {
                        "call_id": projection.call_id,
                        "projected": {
                            item.token: item.value for item in projection.fields
                        },
                    },
                    separators=(",", ":"),
                ),
                timeout=3,
            )
            return ModelResponse(
                response_id="resp_" + token_hex(8),
                fields=(SyntheticField(PROJECTED_FIELD, text),),
            )

        transport = SyntheticSuggestionTransport(
            authorize=lambda: self._period if self._egress_allowed() else None,
            clock=self._clock,
            timeout_seconds=3,
            responder=respond,
        )
        delivery = transport.deliver(request)
        if isinstance(delivery, Refusal):
            self.delivery_refusal = delivery.reason
            raise NativeRejected(delivery.reason)
        result = transport.adopt(delivery)
        value = result.fields[0].value
        if not value.isdecimal():
            raise NativeRejected("malformed_response")
        self._receipts.append(
            {
                "call_id": result.call_id,
                "response_id": result.response_id,
                "source": result.source,
                "contract_version": result.contract_version,
                "field_tokens": [item.token for item in result.fields],
            }
        )
        proposal["changes"][0]["value"] = int(value)
        return json.dumps(proposal)

    def display_model(self) -> dict[str, str]:
        metadata = self.rpc("session/metadata", {"sessionId": self.session_id})
        # The remaining fields (title and cwd included) stay in the fixture process.
        assert set(metadata) == {
            "sessionId",
            "title",
            "cwd",
            "permissionMode",
            "modelAlias",
            "modelName",
            "providerName",
            "effort",
        }
        return {
            key: metadata[key] for key in ("modelAlias", "modelName", "providerName")
        }


@pytest.fixture
def local_peri(reference_ready_app, monkeypatch):
    binary = Path(os.environ.get("HELIXWEAVE_PERI_BINARY", ""))
    if not binary_matches(binary):
        pytest.fail("HELIXWEAVE_PERI_BINARY must point at the pinned Peri binary")
    app = reference_ready_app
    service = app.state.authoring_service
    project = service._data.create_project("Synthetic Peri binding")
    selection = LocalPeri(Endpoint("a", "7"))
    with TemporaryDirectory(prefix="hw-p2-local-") as state_root:
        bridge = NativeAuthoringBridge(
            binary=binary,
            state_root=Path(state_root),
            repository=service._repository,
            runtime=selection,
        )
        monkeypatch.setattr(service, "_bridge", bridge)
        bridge.bind_tool_handler(service.read_tool)
        try:
            with ApiTestClient(app) as client:
                created = client.post(
                    API,
                    json={
                        "project_id": project.project_id,
                        "workflow_id": "encode-style-chipseq-cuttag-atac-mnase",
                    },
                )
                assert created.status_code == 200
                yield client, f"{API}/{created.json()['sid']}", selection
        finally:
            bridge.shutdown()


def _open(client, url):
    opened = client.post(url + "/connection")
    assert opened.status_code == 200
    assert opened.json()["connection"]["state"] == "ready"
    return opened.content


def _turn(client, url, turn_id):
    response = client.post(
        url + "/turns",
        json={
            "client_turn_id": turn_id,
            "expected_revision": 0,
            "question": {
                "intent": "suggest",
                "field_tokens": ["config.threads"],
                "preferences": [{"field_token": "config.threads", "value": 4}],
            },
        },
    )
    assert response.status_code == 200
    return response


def _check_projection(record, receipt):
    prompt = json.dumps(
        {"call_id": receipt["call_id"], "projected": {PROJECTED_FIELD: "4"}},
        separators=(",", ":"),
    )
    assert record["body"]["messages"][-1] == {"role": "user", "content": prompt}
    assert set(receipt) == {
        "call_id",
        "response_id",
        "source",
        "contract_version",
        "field_tokens",
    }
    assert receipt["source"] == "test_transport"
    assert receipt["field_tokens"] == [PROJECTED_FIELD]


@pytest.mark.platform_real_execution
def test_native_a_and_b_drive_the_same_authoring_entry(
    local_peri, tmp_path, capfd, caplog
):
    client, url, selection = local_peri
    endpoint_a = selection.selected
    endpoint_b = Endpoint("b", "11")
    surfaces = []
    identity = None
    first_socket = None
    first_period = None
    for endpoint in (endpoint_a, endpoint_b):
        selection.activate(endpoint)
        selection.next_grant = LocalGrant()
        surfaces.append(_open(client, url))
        native = selection.session
        model = native.display_model()
        if identity is None:
            identity, first_socket = model, endpoint.stub.path
            first_period = native._period
        else:
            assert model == identity
            assert endpoint.stub.path != first_socket
            assert native._period != first_period
        qualification_count = len(endpoint.stub.requests)
        endpoint.answer_once()
        turn = _turn(client, url, f"suggest-{endpoint.name}")
        body = turn.json()
        assert body["status"] == "completed"
        [change] = body["proposal"]["changes"]
        assert change["field_token"] == "config.threads"
        assert change["value"] == int(endpoint.value)
        readback = client.get(f"{url}/turns/{body['turn_id']}")
        assert readback.status_code == 200
        assert readback.json()["proposal"] == body["proposal"]
        surfaces += [turn.content, readback.content]
        assert len(endpoint.stub.requests) == qualification_count + 1
        record = endpoint.stub.requests[-1]
        assert record["path"] == endpoint.route + "/chat/completions"
        _check_projection(record, selection.receipts[-1])
        assert native.handler_calls == 0
        assert not endpoint.stub.failures
        logs = list(native.directory.rglob("*.log*"))
        assert native.directory / "stderr.log" in logs
        surfaces.extend(log_path.read_bytes() for log_path in logs)
        assert client.request("DELETE", url + "/connection").status_code == 200

    assert len(selection.receipts) == 2
    surfaces.append(json.dumps(selection.receipts).encode())
    assert (tmp_path / "platform.db").is_file()
    surfaces.extend(path.read_bytes() for path in tmp_path.glob("platform.db*"))
    captured = capfd.readouterr()
    surfaces.extend(
        (captured.out.encode(), captured.err.encode(), caplog.text.encode())
    )
    for endpoint in (endpoint_a, endpoint_b):
        for private in (endpoint.key, endpoint.route):
            assert all(private.encode() not in surface for surface in surfaces)


@pytest.mark.platform_real_execution
def test_native_missing_binding_and_injected_switch_errors(local_peri):
    client, url, selection = local_peri
    endpoint = selection.selected
    _open(client, url)
    native = selection.session
    count = len(endpoint.stub.requests)
    with pytest.raises(NativeRejected, match="^ACP_REMOTE_REJECTION$"):
        native.rpc(
            "session/set_config_option",
            {"sessionId": native.session_id, "configId": "model", "value": "haiku"},
        )
    error = next(
        frame["error"] for frame in reversed(native.events) if "error" in frame
    )
    assert error["code"] == -32603
    assert error["message"] == "Failed to persist config"
    assert len(endpoint.stub.requests) == count
    assert client.request("DELETE", url + "/connection").status_code == 200

    missing = Endpoint("missing", "")
    selection.activate(missing)
    opened = client.post(url + "/connection")
    assert opened.status_code == 503
    assert opened.json()["reason_code"] == "ADMISSION_FAILED"
    assert all(
        marker.encode() not in opened.content for marker in (missing.key, missing.route)
    )
    assert missing.stub.requests == []
    assert selection.receipts == []


@pytest.mark.parametrize("gate", ["rebind", "revoke", "deadline"])
@pytest.mark.platform_real_execution
def test_inflight_reply_cannot_use_the_old_grant(local_peri, gate):
    client, url, selection = local_peri
    endpoint = selection.selected
    budget = selection.next_grant
    if gate == "deadline":
        selection.clock = lambda: 4.0 if endpoint.release.is_set() else 0.0
    _open(client, url)
    native = selection.session
    endpoint.answer_once()
    endpoint.paused = True
    before = len(endpoint.stub.requests)
    with ThreadPoolExecutor(max_workers=2) as pool:
        pending = pool.submit(_turn, client, url, "pending")
        try:
            assert endpoint.received.wait(5)
            if gate == "rebind":
                selection.activate(Endpoint("b", "11"))
                selection.next_grant = LocalGrant()
            elif gate == "revoke":
                closing = pool.submit(client.request, "DELETE", url + "/connection")
                # The repository fences the old generation before closing the native.
                deadline = monotonic() + 1
                while client.get(url).json()["connection"]["state"] != "closed":
                    assert monotonic() < deadline
        finally:
            endpoint.release.set()
        result = pending.result(timeout=8).json()
        if gate == "revoke":
            assert closing.result(timeout=8).status_code == 200
    assert result["status"] in {"failed", "interrupted"}
    assert result["proposal"] is None
    assert result["adoptable"] is False
    assert selection.receipts == []
    assert len(endpoint.stub.requests) == before + 1
    assert client.get(url).json()["draft_revision"] == 0
    assert budget.remaining == (0, 0, 0, 0)
    if gate == "deadline":
        assert native.delivery_refusal == "timeout"
        endpoint.paused = False
        _open(client, url)
        before_retry = len(endpoint.stub.requests)
        retried = _turn(client, url, "second-request").json()
        assert retried["status"] == "failed"
        assert selection.session.preflight_failure == "budget_exhausted"
        assert len(endpoint.stub.requests) == before_retry


@pytest.mark.parametrize(
    ("grant", "reason"),
    [
        ({"allowed": frozenset()}, "outbound_field_not_granted"),
        ({"remaining": None}, "budget_unproven"),
        ({"remaining": (0, 512, 64, 0)}, "budget_exhausted"),
    ],
)
@pytest.mark.platform_real_execution
def test_fields_and_budget_refuse_before_the_business_request(
    local_peri, grant, reason
):
    client, url, selection = local_peri
    selection.next_grant = LocalGrant(**grant)
    _open(client, url)
    endpoint = selection.selected
    before = len(endpoint.stub.requests)
    result = _turn(client, url, "denied").json()
    assert result["status"] == "failed"
    assert result["proposal"] is None
    assert len(endpoint.stub.requests) == before
    assert selection.receipts == []
    assert selection.session.preflight_failure == reason
