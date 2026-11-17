"""P1a boundary: the owner keeps grants, and Peri stays behind the runtime seam."""

from __future__ import annotations

import ast
from pathlib import Path

import pytest

from encode_pipeline.services.authoring_bridge import (
    BridgeDisconnected,
    BridgeSessionRequest,
    BridgeTurnRequest,
    UnavailableAuthoringBridge,
)
from encode_pipeline.services.authoring_owner_client import (
    PERI_OWNER_SOCKET_ENV,
    load_authoring_bridge,
)
from encode_pipeline.services.peri_runtime import owner as owner_module
from encode_pipeline.services.peri_runtime import session as session_module
from encode_pipeline.services.peri_runtime.owner import NativeAuthoringBridge
from encode_pipeline.services.peri_runtime.runtime import PeriAuthoringRuntime
from encode_pipeline.services.peri_runtime.socket_server import RemoteAuthoringBridge
from encode_pipeline.services.peri_runtime.stub import (
    ModelStub,
    open_qualification_model_transport,
)


ROOT = Path(__file__).resolve().parents[2]


class _SubstituteSession:
    """Stand-in session. It records calls and refuses a platform capability token."""

    def __init__(self, directory: Path, service_generation: int) -> None:
        self.directory = directory
        self.service_generation = service_generation
        self.instance = "substitute-instance"
        self.state = "blocked"
        self.handler_calls = 0
        self.tool_arrivals = 0
        self.errors: list[str] = []
        self.approvals: list[dict[str, object]] = []
        self.denied_callbacks: list[str] = []
        self.authorization_check = None
        self.turns: list[str] = []

    def __setattr__(self, name: str, value: object) -> None:
        if name == "capability_token":
            raise AssertionError("runtime received the capability token")
        object.__setattr__(self, name, value)

    def initialize(self) -> None:
        self.state = "ready"

    def run_authoring_turn(self, text: str) -> str:
        allowed = self.authorization_check is not None and self.authorization_check()
        if not allowed:
            raise RuntimeError("EGRESS_DENIED")
        self.turns.append(text)
        return text

    def close(self) -> None:
        self.state = "closed"


class _SubstituteRuntime:
    def __init__(self) -> None:
        self.session: _SubstituteSession | None = None

    def open_session(self, binary, directory, *, service_generation, tool_handler):
        del binary, tool_handler
        directory.mkdir(parents=True, exist_ok=True)
        self.session = _SubstituteSession(directory, service_generation)
        return self.session


def test_api_module_does_not_import_the_peri_transport_or_session() -> None:
    source = (ROOT / "src/encode_pipeline/api/main.py").read_text(encoding="utf-8")
    modules = set()
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.ImportFrom) and node.module is not None:
            modules.add(node.module)
        elif isinstance(node, ast.Import):
            modules.update(alias.name for alias in node.names)
    assert "encode_pipeline.services.peri_runtime.socket_server" not in modules
    assert "encode_pipeline.services.peri_runtime.session" not in modules
    assert "NativeSession" not in source


def test_missing_owner_socket_keeps_manual_authoring_available() -> None:
    bridge = load_authoring_bridge({})
    assert isinstance(bridge, UnavailableAuthoringBridge)
    assert bridge.status().available is False
    assert bridge.status().reason_code == "BRIDGE_UNAVAILABLE"


def test_configured_socket_uses_only_the_private_owner_transport(
    tmp_path: Path,
) -> None:
    bridge = load_authoring_bridge(
        {PERI_OWNER_SOCKET_ENV: str(tmp_path / "missing-owner.sock")}
    )
    assert isinstance(bridge, RemoteAuthoringBridge)
    assert bridge.status().available is False


def test_substitute_runtime_cannot_hold_or_qualify_platform_authority(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(owner_module, "binary_matches", lambda _path: True)
    monkeypatch.setattr(owner_module, "isolation_paths_exist", lambda: True)

    def reject_native(*_args, **_kwargs):
        raise AssertionError("owner constructed the Peri session")

    monkeypatch.setattr(session_module, "NativeSession", reject_native)
    runtime = _SubstituteRuntime()
    bridge = NativeAuthoringBridge(
        binary=tmp_path / "peri",
        state_root=tmp_path / "state",
        repository=type(
            "Repository",
            (),
            {"tool_read_state": staticmethod(lambda **_kwargs: ("ok", None, None))},
        )(),
        runtime=runtime,
    )
    request = BridgeSessionRequest(
        sid="aas_" + "ab" * 16,
        generation=3,
        owner_epoch=bridge.owner_epoch(),
        capability_token="platform-token",
        schema_digest="a" * 64,
        projection_digest="b" * 64,
        auth_session_digest="c" * 64,
        user_id="usr_" + "d" * 32,
    )
    admission = bridge.admit(request)
    assert admission.instance_id == "substitute-instance"
    assert runtime.session is not None
    assert not hasattr(runtime.session, "capability_token")
    assert bridge._grants[request.sid].capability_token == "platform-token"
    answer = bridge.execute_turn(
        BridgeTurnRequest(
            sid=request.sid,
            turn_id="atr_" + "cd" * 16,
            generation=3,
            instance_id=admission.instance_id,
            owner_epoch=request.owner_epoch,
            capability_token="caller-token-is-not-authority",
            question={"field_tokens": ["config.threads"], "preferences": []},
            issue_codes=(),
            schema_digest=request.schema_digest,
            projection_digest=request.projection_digest,
            base_revision=1,
        )
    )
    assert answer
    assert runtime.session.turns
    with pytest.raises(BridgeDisconnected):
        bridge.execute_turn(
            BridgeTurnRequest(
                sid=request.sid,
                turn_id="atr_" + "ef" * 16,
                generation=2,
                instance_id=admission.instance_id,
                owner_epoch=request.owner_epoch,
                capability_token="platform-token",
                question={},
                issue_codes=(),
                schema_digest=request.schema_digest,
                projection_digest=request.projection_digest,
                base_revision=1,
            )
        )
    bridge.close(request)
    assert runtime.session.state == "closed"


def test_qualification_transport_is_only_the_labeled_stub(tmp_path: Path) -> None:
    stub = open_qualification_model_transport(tmp_path / "model.sock")
    try:
        assert isinstance(stub, ModelStub)
        assert type(stub) is ModelStub
    finally:
        stub.close()


def test_admitted_peri_runtime_selects_the_qualification_transport(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen: dict[str, object] = {}

    def fake_native(*_args, **kwargs):
        seen["model_transport"] = kwargs["model_transport"]
        return object()

    monkeypatch.setattr(session_module, "NativeSession", fake_native)
    monkeypatch.setattr(
        "encode_pipeline.services.peri_runtime.runtime.NativeSession",
        fake_native,
    )
    PeriAuthoringRuntime().open_session(
        Path("peri"),
        Path("unused"),
        service_generation=1,
        tool_handler=lambda *_args: {},
    )
    assert seen["model_transport"] is open_qualification_model_transport
