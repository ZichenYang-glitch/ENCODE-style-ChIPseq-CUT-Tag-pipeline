"""Fast-suite regressions for the bound Peri grant CLI.

The Chromium acceptance module is not named test_*.py, so the deterministic
suite never imports the operator CLI or the subprocess transport. These tests
use a local worker and a synthetic draft row. They do not call a provider.
"""

from __future__ import annotations

from contextlib import closing
from dataclasses import replace
import json
import os
from pathlib import Path
import shlex
import sqlite3
import stat
import subprocess
import sys
import time

import pytest

from encode_pipeline.platform.authoring import (
    PROPOSAL_INVALIDATED,
    AuthoringActor,
    AuthoringRejected,
)
from encode_pipeline.services.authoring import AuthoringService
from encode_pipeline.services.authoring_bridge import (
    BridgeDisconnected,
    BridgeStatus,
    BridgeTurnRefused,
    BridgeTurnRequest,
)
from encode_pipeline.services.peri_runtime.bound_model import (
    BoundModel,
    Budget,
    binding_id,
    configured_model,
)
from encode_pipeline.services.peri_runtime.model_grant import main
from encode_pipeline.services.peri_runtime.owner import (
    NativeAuthoringBridge,
    _PlatformGrant,
)
from encode_pipeline.services.peri_runtime.socket_server import (
    OwnerServer,
    RemoteAuthoringBridge,
)

_NAMES = (
    "HELIXWEAVE_PERI_MODEL_LEDGER",
    "HELIXWEAVE_PERI_MODEL_WORKER",
    "HELIXWEAVE_PERI_MODEL_SETTINGS",
)
_PROVEN = Budget(1, 4096, 128, 1, 128, 100, 100, "CNY")
_BINDING = {"model": "fixture-number"}


def _install_worker(root: Path, mode: str) -> tuple[Path, Path]:
    script = root / "worker.py"
    log = root / "calls.jsonl"
    script.write_text(
        "\n".join(
            [
                "import json, os, sys",
                f"mode = {mode!r}",
                f"log = {str(log)!r}",
                "req = json.load(sys.stdin)",
                "with open(log, 'a', encoding='utf-8') as handle:",
                "    handle.write(json.dumps({'env': sorted(os.environ), 'req': req}) + '\\n')",
                "if req.get('op') == 'binding':",
                "    json.dump({'model': 'fixture-number'}, sys.stdout)",
                "elif mode == 'exit':",
                "    raise SystemExit(3)",
                "elif mode == 'text':",
                "    sys.stdout.write('not-json')",
                "elif mode == 'list':",
                "    json.dump([], sys.stdout)",
                "else:",
                "    json.dump(",
                "        {'binding': req['binding'], 'value': 7}, sys.stdout",
                "    )",
                "",
            ]
        ),
        encoding="utf-8",
    )
    launcher = root / "worker"
    launcher.write_text(
        "#!/bin/sh\n"
        "exec "
        f"{shlex.quote(sys.executable)} -I -S {shlex.quote(str(script))}\n",
        encoding="utf-8",
    )
    launcher.chmod(0o700)
    return launcher, log


def _recorded(log: Path) -> list[dict[str, object]]:
    if not log.exists():
        return []
    rows = [
        json.loads(line)
        for line in log.read_text(encoding="utf-8").splitlines()
        if line
    ]
    for row in rows:
        env = row["env"]
        assert "HW_BOUND_MODEL_SENTINEL" not in env
        for name in _NAMES:
            assert name not in env
    return rows


def _suggests(log: Path) -> list[dict[str, object]]:
    return [row["req"] for row in _recorded(log) if row["req"]["op"] == "suggest"]


def _grant_balance(ledger: Path) -> tuple[int, int, int, int, int] | None:
    with closing(sqlite3.connect(ledger)) as db:
        row = db.execute(
            "SELECT requests,input_left,output_left,minor_left,revoked FROM grants"
        ).fetchone()
    return None if row is None else tuple(row)


def _delivered(ledger: Path) -> int | None:
    with closing(sqlite3.connect(ledger)) as db:
        row = db.execute("SELECT delivered FROM reservations").fetchone()
    return None if row is None else row[0]


class Operator:
    def __init__(self, root: Path, mode: str = "ok") -> None:
        self.ledger = root / "ledger.sqlite"
        self.settings = root / "settings.json"
        self.platform = root / "platform.sqlite"
        self.settings.write_text("{}", encoding="utf-8")
        self.worker, self.log = _install_worker(root, mode)

    def apply(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setenv(_NAMES[0], str(self.ledger))
        monkeypatch.setenv(_NAMES[1], str(self.worker))
        monkeypatch.setenv(_NAMES[2], str(self.settings))
        monkeypatch.setenv("HW_BOUND_MODEL_SENTINEL", "do-not-forward")

    def connect(self, *, state: str = "ready", epoch: str = "epoch") -> None:
        with sqlite3.connect(self.platform) as db:
            db.execute(
                """CREATE TABLE agent_authoring_sessions (
                    sid TEXT, instance_id TEXT, generation INTEGER,
                    bridge_owner_epoch TEXT, connection_state TEXT
                )"""
            )
            db.execute(
                "INSERT INTO agent_authoring_sessions VALUES (?,?,?,?,?)",
                ("sid", "instance", 2, epoch, state),
            )

    def env(self) -> dict[str, str]:
        return {
            _NAMES[0]: str(self.ledger),
            _NAMES[1]: str(self.worker),
            _NAMES[2]: str(self.settings),
            "HW_BOUND_MODEL_SENTINEL": "do-not-forward",
        }


def _invoke(monkeypatch: pytest.MonkeyPatch, *args: str) -> None:
    monkeypatch.setattr(sys, "argv", ["model_grant", *args])
    main()


def _turn(epoch: str, field: str = "config.threads") -> BridgeTurnRequest:
    return BridgeTurnRequest(
        sid="sid",
        turn_id="turn",
        generation=2,
        instance_id="instance",
        owner_epoch=epoch,
        capability_token="secret",
        question={
            "intent": "suggest",
            "field_tokens": [field],
            "preferences": [{"field_token": field, "value": 3}],
        },
        issue_codes=(),
        schema_digest="schema",
        projection_digest="projection",
        base_revision=1,
    )


def _grant_args(operator: Operator, binding: str, **overrides: str) -> list[str]:
    args = {
        "--platform-db": str(operator.platform),
        "--sid": "sid",
        "--field": "config.threads",
        "--binding": binding,
        "--expires": str(int(time.time()) + 300),
        "--requests": "1",
        "--input-tokens": "4096",
        "--output-tokens": "128",
        "--minor-units": "1",
        "--output-per-call": "128",
        "--input-rate": "100",
        "--output-rate": "100",
        "--currency": "CNY",
    }
    args.update(overrides)
    flat = ["grant"]
    for key, value in args.items():
        flat.extend((key, value))
    return flat


class LiveSession:
    def __init__(self, directory: Path) -> None:
        directory.mkdir()
        self.instance = "instance"
        self.service_generation = 2
        self.state = "ready"
        self.handler_calls = 0
        self.directory = directory
        self.tool_arrivals = 0
        self.errors: list[str] = []
        self.approvals: list[dict[str, object]] = []
        self.denied_callbacks: list[str] = []
        self.authorization_check = None
        self.closed = False

    def close(self) -> None:
        self.closed = True


class ReadRepository:
    def retire_foreign_epochs(self, _epoch: str) -> int:
        return 0

    def tool_read_state(self, **_kwargs: object) -> tuple[str, None, None]:
        return ("ok", None, None)


def _ready_owner(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys):
    operator = Operator(tmp_path)
    operator.apply(monkeypatch)
    operator.connect()
    _invoke(monkeypatch, "init")
    capsys.readouterr()
    _invoke(monkeypatch, "binding")
    binding = json.loads(capsys.readouterr().out)["binding"]
    model = configured_model(os.environ)
    assert model is not None
    bridge = NativeAuthoringBridge(
        binary=tmp_path / "not-peri",
        state_root=tmp_path / "owner-state",
        repository=ReadRepository(),
        bound_model=model,
    )
    with sqlite3.connect(operator.platform) as db:
        db.execute(
            "UPDATE agent_authoring_sessions SET bridge_owner_epoch=?",
            (bridge._epoch,),
        )
    _invoke(monkeypatch, *_grant_args(operator, binding))
    capsys.readouterr()
    session = LiveSession(tmp_path / "session")
    bridge._sessions["sid"] = session
    bridge._grants["sid"] = _PlatformGrant(
        capability_token="secret",
        user_id="user",
        auth_session_digest="digest",
    )
    return operator, bridge, session


@pytest.mark.parametrize(
    "budget",
    [
        Budget(0, 4096, 128, 1, 128, 100, 100, "CNY"),
        Budget(1, 4096, 128, 1, 128, 100, 100, "cn"),
        Budget(1, 4096, 128, 1, 1025, 100, 100, "CNY"),
    ],
)
def test_unproven_budget_does_not_open_a_grant(tmp_path: Path, budget: Budget) -> None:
    ledger = tmp_path / "ledger.sqlite"
    BoundModel.initialize(ledger)
    model = BoundModel(ledger, lambda _request: _BINDING)
    with pytest.raises(ValueError, match="budget_unproven"):
        model.grant(
            sid="sid",
            instance="instance",
            generation=2,
            epoch="epoch",
            field="config.threads",
            binding="a" * 64,
            expires=10**12,
            budget=budget,
        )
    assert _grant_balance(ledger) is None


@pytest.mark.parametrize(
    ("field", "binding", "expires"),
    [
        ("Nope", "a" * 64, 10**12),
        ("config.threads", "not-a-binding", 10**12),
        ("config.threads", "a" * 64, 0),
    ],
)
def test_invalid_grant_fields_do_not_open_a_grant(
    tmp_path: Path, field: str, binding: str, expires: int
) -> None:
    ledger = tmp_path / "ledger.sqlite"
    BoundModel.initialize(ledger)
    model = BoundModel(ledger, lambda _request: _BINDING)
    with pytest.raises(ValueError, match="invalid_grant"):
        model.grant(
            sid="sid",
            instance="instance",
            generation=2,
            epoch="epoch",
            field=field,
            binding=binding,
            expires=expires,
            budget=_PROVEN,
        )
    assert _grant_balance(ledger) is None


def test_cli_init_is_private_and_refuses_replacement(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    operator = Operator(tmp_path)
    operator.apply(monkeypatch)
    _invoke(monkeypatch, "init")
    capsys.readouterr()
    assert stat.S_IMODE(operator.ledger.stat().st_mode) == 0o600
    inode = operator.ledger.stat().st_ino
    with pytest.raises(FileExistsError):
        _invoke(monkeypatch, "init")
    assert operator.ledger.stat().st_ino == inode
    with closing(sqlite3.connect(operator.ledger)) as db:
        names = [
            row[0]
            for row in db.execute(
                "SELECT name FROM sqlite_master WHERE type='table' ORDER BY name"
            )
        ]
    assert names == ["grants", "reservations"]


def test_module_entry_initializes_the_same_ledger(tmp_path: Path) -> None:
    operator = Operator(tmp_path)
    env = os.environ.copy()
    env.update(operator.env())
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "encode_pipeline.services.peri_runtime.model_grant",
            "init",
        ],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert stat.S_IMODE(operator.ledger.stat().st_mode) == 0o600
    again = subprocess.run(
        [
            sys.executable,
            "-m",
            "encode_pipeline.services.peri_runtime.model_grant",
            "init",
        ],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert again.returncode != 0
    assert operator.ledger.stat().st_size > 0


def test_cli_rejects_missing_and_partial_configuration(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    for name in _NAMES:
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setattr(sys, "argv", ["model_grant", "init"])
    with pytest.raises(SystemExit) as missing:
        main()
    assert missing.value.code == 2
    assert "HELIXWEAVE_PERI_MODEL" in capsys.readouterr().err
    monkeypatch.setenv(_NAMES[0], str(tmp_path / "ledger.sqlite"))
    with pytest.raises(ValueError, match="incomplete_bound_model_configuration"):
        main()
    assert not (tmp_path / "ledger.sqlite").exists()


def test_cli_rejects_binding_draft_field_and_budget(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    operator = Operator(tmp_path)
    operator.apply(monkeypatch)
    operator.connect()
    _invoke(monkeypatch, "init")
    capsys.readouterr()
    _invoke(monkeypatch, "binding")
    issued = json.loads(capsys.readouterr().out)
    assert issued == {
        "model": "fixture-number",
        "binding": binding_id(_BINDING),
    }
    binding = issued["binding"]

    monkeypatch.setattr(sys, "argv", ["model_grant", *_grant_args(operator, "ab" * 32)])
    with pytest.raises(SystemExit) as changed:
        main()
    assert changed.value.code == 2
    assert "binding changed" in capsys.readouterr().err
    assert _grant_balance(operator.ledger) is None

    with sqlite3.connect(operator.platform) as db:
        db.execute("UPDATE agent_authoring_sessions SET connection_state='closed'")
    monkeypatch.setattr(sys, "argv", ["model_grant", *_grant_args(operator, binding)])
    with pytest.raises(SystemExit) as unready:
        main()
    assert unready.value.code == 2
    assert "connect the draft first" in capsys.readouterr().err
    assert _grant_balance(operator.ledger) is None

    with sqlite3.connect(operator.platform) as db:
        db.execute("UPDATE agent_authoring_sessions SET connection_state='ready'")
    monkeypatch.setattr(
        sys,
        "argv",
        ["model_grant", *_grant_args(operator, binding, **{"--field": "Nope"})],
    )
    with pytest.raises(ValueError, match="invalid_grant"):
        main()
    assert _grant_balance(operator.ledger) is None
    monkeypatch.setattr(
        sys,
        "argv",
        ["model_grant", *_grant_args(operator, binding, **{"--requests": "0"})],
    )
    with pytest.raises(ValueError, match="budget_unproven"):
        main()
    assert _grant_balance(operator.ledger) is None
    assert _suggests(operator.log) == []


@pytest.mark.parametrize(
    ("requests", "input_tokens", "output_tokens", "minor_units"),
    [
        (1, 8192, 256, 2),
        (2, 4096, 256, 2),
        (2, 8192, 128, 2),
        (2, 8192, 256, 1),
    ],
    ids=["request-count", "input-tokens", "output-tokens", "cost"],
)
def test_first_exhausted_budget_stops_the_next_request(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    requests: int,
    input_tokens: int,
    output_tokens: int,
    minor_units: int,
) -> None:
    operator = Operator(tmp_path)
    operator.apply(monkeypatch)
    operator.connect()
    _invoke(monkeypatch, "init")
    capsys.readouterr()
    _invoke(monkeypatch, "binding")
    binding = json.loads(capsys.readouterr().out)["binding"]
    limits = {
        "--requests": str(requests),
        "--input-tokens": str(input_tokens),
        "--output-tokens": str(output_tokens),
        "--minor-units": str(minor_units),
    }
    _invoke(monkeypatch, *_grant_args(operator, binding, **limits))
    assert json.loads(capsys.readouterr().out) == {
        "instance": "instance",
        "field": "config.threads",
    }
    assert _grant_balance(operator.ledger) == (
        requests,
        input_tokens,
        output_tokens,
        minor_units,
        0,
    )
    model = configured_model(os.environ)
    assert model is not None
    text = model.suggest(_turn("epoch"), lambda: True)
    body = json.loads(text)
    assert body["changes"] == [{"field_token": "config.threads", "value": 7}]
    assert body["kind"] == "proposal"
    assert "sid" not in text
    remaining = (
        requests - 1,
        input_tokens - 4096,
        output_tokens - 128,
        minor_units - 1,
        0,
    )
    assert _grant_balance(operator.ledger) == remaining
    assert _delivered(operator.ledger) == 1
    with pytest.raises(BridgeTurnRefused):
        model.suggest(replace(_turn("epoch"), turn_id="next-turn"), lambda: True)
    assert _suggests(operator.log) == [
        {
            "op": "suggest",
            "binding": _BINDING,
            "field": "config.threads",
            "value": 3,
            "output_tokens": 128,
        }
    ]
    assert _grant_balance(operator.ledger) == remaining
    with closing(sqlite3.connect(operator.ledger)) as db:
        assert db.execute("SELECT COUNT(*) FROM reservations").fetchone()[0] == 1


def test_authorization_drop_after_reserve_does_not_call_the_provider(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    operator = Operator(tmp_path)
    operator.apply(monkeypatch)
    operator.connect()
    _invoke(monkeypatch, "init")
    capsys.readouterr()
    _invoke(monkeypatch, "binding")
    binding = json.loads(capsys.readouterr().out)["binding"]
    _invoke(monkeypatch, *_grant_args(operator, binding))
    capsys.readouterr()
    model = configured_model(os.environ)
    assert model is not None
    seen = {"n": 0}

    def authorized() -> bool:
        seen["n"] += 1
        return seen["n"] < 3

    with pytest.raises(BridgeDisconnected):
        model.suggest(_turn("epoch"), authorized)
    assert seen["n"] == 3
    assert _suggests(operator.log) == []
    assert _grant_balance(operator.ledger) == (0, 0, 0, 0, 0)
    assert _delivered(operator.ledger) == 0

    with pytest.raises(BridgeTurnRefused):
        model.suggest(replace(_turn("epoch"), turn_id="retry"), lambda: True)
    assert _suggests(operator.log) == []


@pytest.mark.parametrize("mode", ["exit", "text", "list"])
def test_failed_mock_delivery_keeps_the_charge(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mode: str
) -> None:
    operator = Operator(tmp_path, mode)
    operator.apply(monkeypatch)
    _invoke(monkeypatch, "init")
    model = configured_model(os.environ)
    assert model is not None
    binding = binding_id(model.worker({"op": "binding"}))
    model.grant(
        sid="sid",
        instance="instance",
        generation=2,
        epoch="epoch",
        field="config.threads",
        binding=binding,
        expires=int(time.time()) + 300,
        budget=_PROVEN,
    )
    with pytest.raises(BridgeDisconnected):
        model.suggest(_turn("epoch"), lambda: True)
    assert len(_suggests(operator.log)) == 1
    assert _grant_balance(operator.ledger) == (0, 0, 0, 0, 0)
    assert _delivered(operator.ledger) == 0
    with pytest.raises(BridgeTurnRefused):
        model.suggest(replace(_turn("epoch"), turn_id="retry"), lambda: True)
    assert len(_suggests(operator.log)) == 1


def test_missing_worker_fails_before_a_charge(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    ledger = tmp_path / "ledger.sqlite"
    BoundModel.initialize(ledger)
    monkeypatch.setenv(_NAMES[0], str(ledger))
    monkeypatch.setenv(_NAMES[1], str(tmp_path / "absent-worker"))
    monkeypatch.setenv(_NAMES[2], str(tmp_path / "settings.json"))
    model = configured_model(os.environ)
    assert model is not None
    with pytest.raises(BridgeDisconnected):
        model.worker({"op": "binding"})
    assert _grant_balance(ledger) is None


def test_revoke_blocks_delivery_and_leaves_the_unused_request(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    operator = Operator(tmp_path)
    operator.apply(monkeypatch)
    operator.connect()
    _invoke(monkeypatch, "init")
    capsys.readouterr()
    _invoke(monkeypatch, "binding")
    binding = json.loads(capsys.readouterr().out)["binding"]
    _invoke(monkeypatch, *_grant_args(operator, binding))
    capsys.readouterr()
    monkeypatch.setattr(sys, "argv", ["model_grant", "revoke", "--instance", "other"])
    with pytest.raises(SystemExit) as unknown:
        main()
    assert unknown.value.code == 2
    assert "unknown instance" in capsys.readouterr().err
    assert _grant_balance(operator.ledger) == (1, 4096, 128, 1, 0)
    _invoke(monkeypatch, "revoke", "--instance", "instance")
    model = configured_model(os.environ)
    assert model is not None
    with pytest.raises(BridgeDisconnected):
        model.suggest(_turn("epoch"), lambda: True)
    assert _suggests(operator.log) == []
    assert _grant_balance(operator.ledger) == (1, 4096, 128, 1, 1)
    assert _delivered(operator.ledger) is None


def test_owner_returns_the_mock_suggestion_and_rechecks_adoption(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    operator, bridge, session = _ready_owner(tmp_path, monkeypatch, capsys)
    text = bridge.execute_turn(_turn(bridge._epoch))
    assert json.loads(text)["changes"] == [
        {"field_token": "config.threads", "value": 7}
    ]
    assert session.closed is False
    bridge.check_adoption("sid", bridge._epoch, "instance", 2, "turn")
    assert _suggests(operator.log) == [
        {
            "op": "suggest",
            "binding": _BINDING,
            "field": "config.threads",
            "value": 3,
            "output_tokens": 128,
        }
    ]
    assert _grant_balance(operator.ledger) == (0, 0, 0, 0, 0)


@pytest.mark.parametrize("refusal", ["field", "unsupported-intent"])
def test_owner_refuses_one_turn_without_losing_the_session(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    refusal: str,
) -> None:
    operator, bridge, session = _ready_owner(tmp_path, monkeypatch, capsys)
    turn = _turn(bridge._epoch, "private")
    if refusal == "unsupported-intent":
        turn = replace(turn, question={"intent": "explain_issues"})
    with pytest.raises(BridgeTurnRefused):
        bridge.execute_turn(turn)
    assert session.closed is False
    assert _suggests(operator.log) == []
    assert _grant_balance(operator.ledger) == (1, 4096, 128, 1, 0)
    answer = bridge.execute_turn(replace(_turn(bridge._epoch), turn_id="authorized"))
    assert json.loads(answer)["changes"] == [
        {"field_token": "config.threads", "value": 7}
    ]
    assert len(_suggests(operator.log)) == 1


def test_unbound_owner_skips_adoption_recheck(tmp_path: Path) -> None:
    bridge = NativeAuthoringBridge(
        binary=tmp_path / "not-peri",
        state_root=tmp_path / "owner-state",
        bound_model=None,
    )
    assert bridge.check_adoption("sid", bridge._epoch, "instance", 2, "turn") is None


def test_owner_ledger_error_is_a_disconnect(tmp_path: Path) -> None:
    model = BoundModel(tmp_path / "absent.sqlite", lambda _request: _BINDING)
    bridge = NativeAuthoringBridge(
        binary=tmp_path / "not-peri",
        state_root=tmp_path / "owner-state",
        bound_model=model,
    )
    with pytest.raises(BridgeDisconnected):
        bridge.check_adoption("sid", "epoch", "instance", 2, "turn")


@pytest.mark.parametrize("has_check", [True, False], ids=["disconnected", "missing"])
def test_failed_adoption_recheck_is_a_public_rejection(has_check: bool) -> None:
    turn = type(
        "Turn",
        (),
        {
            "decision": "pending",
            "bridge_owner_epoch": "stored-epoch",
            "instance_id": "stored-instance",
            "generation": 4,
            "turn_id": "stored-turn",
        },
    )()
    seen: list[tuple[object, ...]] = []

    class Repository:
        def get_proposal(self, actor: object, sid: str, proposal_id: str):
            del actor, sid, proposal_id
            return object(), turn

        def adopt(self, *_args: object, **_kwargs: object) -> None:
            raise AssertionError("rejected proposal was written")

    class Bridge:
        def check_adoption(
            self, sid: str, epoch: str, instance: str, generation: int, turn_id: str
        ) -> None:
            seen.append((sid, epoch, instance, generation, turn_id))
            raise BridgeDisconnected

    service = AuthoringService(
        repository=Repository(),
        registry=object(),
        data_registry=object(),
        validated_inputs=object(),
        snapshots=object(),
        bridge=Bridge() if has_check else object(),
    )
    with pytest.raises(AuthoringRejected) as rejected:
        service.adopt(
            AuthoringActor(user_id="user", session_digest="digest"),
            "stored-sid",
            "proposal",
            expected_revision=0,
            proposal_digest_value="digest",
        )
    assert rejected.value.reason_code == PROPOSAL_INVALIDATED
    assert rejected.value.message == "The proposal can no longer be adopted."
    assert rejected.value.detail == {}
    expected = ("stored-sid", "stored-epoch", "stored-instance", 4, "stored-turn")
    assert seen == ([expected] if has_check else [])


def test_socket_rechecks_adoption_after_a_grant_is_revoked(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    operator, bridge, _session = _ready_owner(tmp_path, monkeypatch, capsys)
    # Native admission is exercised separately with Chromium and the real worker.
    monkeypatch.setattr(bridge, "status", lambda: BridgeStatus(available=True))
    path = tmp_path / "owner.sock"
    server = OwnerServer(bridge, ReadRepository(), path)
    server.start()
    try:
        remote = RemoteAuthoringBridge(path)
        epoch = remote.owner_epoch()
        turn = _turn(epoch)
        assert json.loads(remote.execute_turn(turn))["changes"] == [
            {"field_token": "config.threads", "value": 7}
        ]
        with pytest.raises(BridgeTurnRefused):
            remote.execute_turn(replace(turn, turn_id="over-budget"))
        remote.check_adoption("sid", epoch, "instance", 2, turn.turn_id)
        with pytest.raises(BridgeDisconnected):
            remote.check_adoption("sid", "stale", "instance", 2, turn.turn_id)
        _invoke(monkeypatch, "revoke", "--instance", "instance")
        with pytest.raises(BridgeDisconnected):
            remote.check_adoption("sid", epoch, "instance", 2, turn.turn_id)
        assert len(_suggests(operator.log)) == 1
    finally:
        server.shutdown()
