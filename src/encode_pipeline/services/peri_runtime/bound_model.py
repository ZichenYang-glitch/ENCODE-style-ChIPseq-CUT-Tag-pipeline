"""Owner-only, opt-in suggestions from a separately bound Peri model worker.

Numeric draft suggestions and catalog target classification share the worker
process and the ledger file. They do not share grants. A ``config.threads``
grant cannot authorize ``experiment.target``, and a catalog grant cannot
authorize a draft turn.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from contextlib import closing
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import re
from secrets import token_hex
import sqlite3
import subprocess
import time

from encode_pipeline.services.authoring_bridge import (
    BridgeDisconnected,
    BridgeTurnRefused,
    BridgeTurnRequest,
)
from encode_pipeline.services.peri_runtime.catalog_classify import (
    CATALOG_FIELD,
    INSUFFICIENT_EVIDENCE,
    MAX_CLASSIFY_CANDIDATES,
    MAX_CLASSIFY_REQUEST_BYTES,
    CatalogTransportError,
    accepted_decision,
    classify_worker_request,
    clue_is_auditable,
    evidence_token,
)

# For the one fixed prompt + bounded integer and field name. The operator must
# attest that this covers the binding's tokenizer/framing and tariff (see SOP).
INPUT_RESERVE = 4096
_MAX = 2**53 - 1
_TOKEN = re.compile(r"[a-z][a-z0-9_.]{0,63}")
_SCHEMA = """
CREATE TABLE grants (
 instance TEXT PRIMARY KEY, sid TEXT NOT NULL, epoch TEXT NOT NULL, generation INTEGER NOT NULL,
 binding TEXT NOT NULL, field TEXT NOT NULL, expires INTEGER NOT NULL, output_cap INTEGER NOT NULL,
 input_rate INTEGER NOT NULL, output_rate INTEGER NOT NULL, requests INTEGER NOT NULL,
 input_left INTEGER NOT NULL, output_left INTEGER NOT NULL, minor_left INTEGER NOT NULL,
 currency TEXT NOT NULL, revoked INTEGER NOT NULL DEFAULT 0);
CREATE TABLE reservations (
 instance TEXT NOT NULL, turn TEXT NOT NULL, model TEXT, input_reserved INTEGER, output_reserved INTEGER, minor_reserved INTEGER, adopt_until REAL, delivered INTEGER NOT NULL DEFAULT 0,
 PRIMARY KEY (instance,turn));
"""


def binding_id(binding: object) -> str:
    return hashlib.sha256(
        json.dumps(binding, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


@dataclass(frozen=True)
class Budget:
    requests: int
    input_tokens: int
    output_tokens: int
    minor_units: int
    output_per_call: int
    input_rate: int
    output_rate: int
    currency: str

    def validate(self) -> None:
        if (
            any(
                type(n) is not int or not 1 <= n <= _MAX
                for k, n in vars(self).items()
                if k != "currency"
            )
            or not re.fullmatch("[A-Z]{3}", self.currency)
            or self.output_per_call > 1024
        ):
            raise ValueError("budget_unproven")


class BoundModel:
    def __init__(self, ledger: Path, worker: Callable[[dict], dict]):
        self.ledger = ledger
        self.worker = worker

    def db(self):
        db = sqlite3.connect(
            self.ledger.resolve().as_uri() + "?mode=rw", uri=True, timeout=2
        )
        db.row_factory = sqlite3.Row
        return db

    @staticmethod
    def initialize(path: Path) -> None:
        # Refuse to replace an existing ledger and thereby refill its balances.
        os.close(os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600))
        with closing(sqlite3.connect(path)) as db:
            db.executescript(_SCHEMA)

    def grant(
        self,
        *,
        sid: str,
        instance: str,
        generation: int,
        epoch: str,
        field: str,
        binding: str,
        expires: int,
        budget: Budget,
    ) -> None:
        budget.validate()
        if (
            not _TOKEN.fullmatch(field)
            or not re.fullmatch("[a-f0-9]{64}", binding)
            or expires <= time.time()
        ):
            raise ValueError("invalid_grant")
        with closing(self.db()) as db, db:
            db.execute(
                "INSERT INTO grants VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,0)",
                (
                    instance,
                    sid,
                    epoch,
                    generation,
                    binding,
                    field,
                    expires,
                    budget.output_per_call,
                    budget.input_rate,
                    budget.output_rate,
                    budget.requests,
                    budget.input_tokens,
                    budget.output_tokens,
                    budget.minor_units,
                    budget.currency,
                ),
            )

    @staticmethod
    def _row(db, sid, epoch, instance, generation, binding):
        row = db.execute(
            "SELECT * FROM grants WHERE instance=? AND sid=? AND epoch=? AND generation=? AND binding=? AND revoked=0 AND expires>?",
            (instance, sid, epoch, generation, binding, int(time.time())),
        ).fetchone()
        if (
            row is None
            or not isinstance(row["currency"], str)
            or not re.fullmatch("[A-Z]{3}", row["currency"])
            or any(
                type(row[key]) is not int or not 0 <= row[key] <= _MAX
                for key in (
                    "requests",
                    "input_left",
                    "output_left",
                    "minor_left",
                    "output_cap",
                    "input_rate",
                    "output_rate",
                    "expires",
                )
            )
            or not (
                1 <= row["output_cap"] <= 1024
                and row["input_rate"]
                and row["output_rate"]
            )
        ):
            raise BridgeDisconnected
        return row

    def check(self, sid, epoch, instance, generation, turn):
        binding = binding_id(self.worker({"op": "binding"}))
        with closing(self.db()) as db:
            self._row(db, sid, epoch, instance, generation, binding)
            if not db.execute(
                "SELECT 1 FROM reservations WHERE instance=? AND turn=? AND delivered=1 AND adopt_until>=?",
                (instance, turn, time.time()),
            ).fetchone():
                raise BridgeDisconnected

    def suggest(self, req: BridgeTurnRequest, authorized: Callable[[], bool]) -> str:
        deadline = time.monotonic() + 35
        q = req.question
        tokens = q.get("field_tokens")
        prefs = q.get("preferences")
        if not authorized():
            raise BridgeDisconnected
        if (
            q.get("intent") != "suggest"
            or not isinstance(tokens, list)
            or len(tokens) != 1
            or not isinstance(tokens[0], str)
            or not _TOKEN.fullmatch(tokens[0])
            or not isinstance(prefs, list)
            or len(prefs) != 1
            or not isinstance(prefs[0], dict)
            or prefs[0].get("field_token") != tokens[0]
            or type(prefs[0].get("value")) is not int
            or not -(2**31) <= prefs[0]["value"] < 2**31
        ):
            raise BridgeTurnRefused
        field = tokens[0]
        binding = self.worker({"op": "binding"})
        identity = binding_id(binding)
        args = (req.sid, req.owner_epoch, req.instance_id, req.generation, identity)
        with closing(self.db()) as db, db:
            db.execute("BEGIN IMMEDIATE")
            row = self._row(db, *args)
            cap = row["output_cap"]
            amount = (
                row["input_rate"] * INPUT_RESERVE + row["output_rate"] * cap + 999999
            ) // 1000000
            if not authorized():
                raise BridgeDisconnected
            if field != row["field"]:
                raise BridgeTurnRefused
            updated = db.execute(
                """UPDATE grants SET requests=requests-1,input_left=input_left-?,
                output_left=output_left-?,minor_left=minor_left-? WHERE instance=? AND requests>=1
                AND input_left>=? AND output_left>=? AND minor_left>=?""",
                (
                    INPUT_RESERVE,
                    cap,
                    amount,
                    req.instance_id,
                    INPUT_RESERVE,
                    cap,
                    amount,
                ),
            )
            if updated.rowcount != 1:
                raise BridgeTurnRefused
            db.execute(
                "INSERT INTO reservations(instance,turn,model,input_reserved,output_reserved,minor_reserved) VALUES(?,?,?,?,?,?)",
                (
                    req.instance_id,
                    req.turn_id,
                    binding["model"],
                    INPUT_RESERVE,
                    cap,
                    amount,
                ),
            )
        # Commit the charge BEFORE I/O. Timeouts/restarts never refund it.
        with closing(self.db()) as db:
            self._row(db, *args)
        if not authorized() or time.monotonic() > deadline:
            raise BridgeDisconnected
        response = self.worker(
            {
                "op": "suggest",
                "binding": binding,
                "field": field,
                "value": prefs[0]["value"],
                "output_tokens": cap,
            }
        )
        v = response.get("value")
        if (
            time.monotonic() > deadline
            or response.get("binding") != binding
            or type(v) is not int
            or not -(2**31) <= v < 2**31
            or not authorized()
            or binding_id(self.worker({"op": "binding"})) != identity
        ):
            raise BridgeDisconnected
        with closing(self.db()) as db, db:
            db.execute("BEGIN IMMEDIATE")
            self._row(db, *args)
            db.execute(
                "UPDATE reservations SET delivered=1, adopt_until=? WHERE instance=? AND turn=?",
                (min(time.time() + 35, row["expires"]), req.instance_id, req.turn_id),
            )
        return json.dumps(
            {
                "base_revision": req.base_revision,
                "schema_digest": req.schema_digest,
                "kind": "proposal",
                "reason_codes": [],
                "changes": [{"field_token": field, "value": v}],
            }
        )

    def _ensure_catalog_schema(self) -> None:
        with closing(self.db()) as db:
            db.executescript(
                """
                CREATE TABLE IF NOT EXISTS catalog_grants (
                 grant_id TEXT PRIMARY KEY, project_id TEXT NOT NULL, field TEXT NOT NULL,
                 binding TEXT NOT NULL, expires INTEGER NOT NULL, output_cap INTEGER NOT NULL,
                 input_rate INTEGER NOT NULL, output_rate INTEGER NOT NULL,
                 requests INTEGER NOT NULL, input_left INTEGER NOT NULL,
                 output_left INTEGER NOT NULL, minor_left INTEGER NOT NULL,
                 currency TEXT NOT NULL, revoked INTEGER NOT NULL DEFAULT 0);
                CREATE TABLE IF NOT EXISTS catalog_reservations (
                 grant_id TEXT NOT NULL, call_key TEXT NOT NULL,
                 projection_digest TEXT NOT NULL, model TEXT,
                 input_reserved INTEGER, output_reserved INTEGER, minor_reserved INTEGER,
                 delivered INTEGER NOT NULL DEFAULT 0, decision TEXT,
                 PRIMARY KEY (grant_id, call_key));
                """
            )

    def catalog_grant(
        self,
        *,
        project_id: str,
        field: str,
        binding: str,
        expires: int,
        budget: Budget,
    ) -> str:
        """Authorize one project field. Draft fields such as config.threads are refused."""

        budget.validate()
        if (
            field != CATALOG_FIELD
            or not re.fullmatch(r"prj_[0-9a-f]{32}", project_id)
            or not re.fullmatch(r"[a-f0-9]{64}", binding)
            or type(expires) is not int
            or expires <= time.time()
        ):
            raise ValueError("invalid_grant")
        grant_id = "catg_" + token_hex(16)
        self._ensure_catalog_schema()
        with closing(self.db()) as db, db:
            db.execute("BEGIN IMMEDIATE")
            active = db.execute(
                """SELECT 1 FROM catalog_grants
                WHERE project_id=? AND field=? AND binding=? AND revoked=0 AND expires>?""",
                (project_id, field, binding, int(time.time())),
            ).fetchone()
            if active is not None:
                raise ValueError("invalid_grant")
            db.execute(
                """INSERT INTO catalog_grants (
                    grant_id, project_id, field, binding, expires, output_cap,
                    input_rate, output_rate, requests, input_left, output_left,
                    minor_left, currency, revoked
                ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,0)""",
                (
                    grant_id,
                    project_id,
                    field,
                    binding,
                    expires,
                    budget.output_per_call,
                    budget.input_rate,
                    budget.output_rate,
                    budget.requests,
                    budget.input_tokens,
                    budget.output_tokens,
                    budget.minor_units,
                    budget.currency,
                ),
            )
        return grant_id

    def catalog_revoke(self, grant_id: str) -> None:
        self._ensure_catalog_schema()
        with closing(self.db()) as db, db:
            if (
                db.execute(
                    "UPDATE catalog_grants SET revoked=1 WHERE grant_id=? AND revoked=0",
                    (grant_id,),
                ).rowcount
                != 1
            ):
                raise ValueError("unknown_grant")

    def classify(
        self,
        *,
        project_id: str,
        idempotency_key: str,
        clue: str,
        candidates: Sequence[tuple[str, str]],
        projection_digest: str,
    ) -> str:
        """Reserve one catalog call, then ask the worker. Timeouts stay charged.

        A stored decision for the same key and projection is returned without
        another worker call. Draft grants are not consulted.
        """

        pairs = tuple(candidates)
        self._validate_classify_inputs(
            project_id=project_id,
            idempotency_key=idempotency_key,
            clue=clue,
            candidates=pairs,
            projection_digest=projection_digest,
        )
        self._ensure_catalog_schema()
        replay = self._catalog_replay(project_id, idempotency_key, projection_digest)
        if replay is not None:
            return replay
        if not self._project_has_catalog_grant(project_id):
            self._raise_without_project_grant(project_id)
        binding = self._read_binding()
        identity = binding_id(binding)
        deadline = time.monotonic() + 35
        reserved = self._reserve_catalog_call(
            project_id=project_id,
            idempotency_key=idempotency_key,
            projection_digest=projection_digest,
            clue=clue,
            candidates=pairs,
            binding=binding,
            identity=identity,
        )
        if isinstance(reserved, str):
            return reserved
        request, grant_id = reserved
        if time.monotonic() > deadline:
            raise CatalogTransportError("catalog_model_timeout")
        response = self._invoke_classify(request)
        accepted = accepted_decision(response.get("decision"), pairs)
        if response.get("binding") != binding or accepted is None:
            raise CatalogTransportError("catalog_model_malformed")
        if time.monotonic() > deadline:
            raise CatalogTransportError("catalog_model_timeout")
        try:
            rebound = binding_id(self._read_binding())
        except CatalogTransportError:
            raise CatalogTransportError("catalog_model_timeout") from None
        if rebound != identity or time.monotonic() > deadline:
            raise CatalogTransportError("catalog_model_timeout")
        self._store_catalog_decision(
            grant_id=grant_id,
            idempotency_key=idempotency_key,
            projection_digest=projection_digest,
            decision=accepted,
        )
        return accepted

    def _validate_classify_inputs(
        self,
        *,
        project_id: str,
        idempotency_key: str,
        clue: str,
        candidates: tuple[tuple[str, str], ...],
        projection_digest: str,
    ) -> None:
        if (
            not re.fullmatch(r"prj_[0-9a-f]{32}", project_id)
            or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,127}", idempotency_key)
            or not re.fullmatch(r"[a-f0-9]{64}", projection_digest)
            or not clue_is_auditable(clue)
            or not 1 <= len(candidates) <= MAX_CLASSIFY_CANDIDATES
        ):
            raise CatalogTransportError("catalog_model_malformed")
        seen: set[str] = set()
        for token, label in candidates:
            if (
                not isinstance(token, str)
                or not isinstance(label, str)
                or evidence_token(label) != token
                or token in seen
            ):
                raise CatalogTransportError("catalog_model_malformed")
            seen.add(token)

    def _catalog_replay(
        self,
        project_id: str,
        idempotency_key: str,
        projection_digest: str,
        *,
        db: sqlite3.Connection | None = None,
    ) -> str | None:
        owned = db is None
        connection = self.db() if owned else db
        assert connection is not None
        try:
            if owned:
                connection.execute("BEGIN IMMEDIATE")
            row = connection.execute(
                """SELECT r.projection_digest, r.delivered, r.decision
                FROM catalog_reservations AS r
                JOIN catalog_grants AS g ON g.grant_id = r.grant_id
                WHERE r.call_key=? AND g.project_id=? AND g.field=?""",
                (idempotency_key, project_id, CATALOG_FIELD),
            ).fetchone()
            if row is None:
                return None
            if row["projection_digest"] != projection_digest:
                raise CatalogTransportError("catalog_idempotency_conflict")
            if (
                row["delivered"] == 1
                and isinstance(row["decision"], str)
                and (
                    row["decision"] == INSUFFICIENT_EVIDENCE
                    or re.fullmatch(r"[a-z][a-z0-9]{0,31}", row["decision"])
                )
            ):
                return str(row["decision"])
            raise CatalogTransportError("catalog_model_timeout")
        finally:
            if owned:
                connection.commit()
                connection.close()

    def _project_has_catalog_grant(self, project_id: str) -> bool:
        with closing(self.db()) as db:
            return (
                db.execute(
                    "SELECT 1 FROM catalog_grants WHERE project_id=? AND field=?",
                    (project_id, CATALOG_FIELD),
                ).fetchone()
                is not None
            )

    def _raise_without_project_grant(self, project_id: str) -> None:
        with closing(self.db()) as db:
            other = db.execute(
                "SELECT 1 FROM catalog_grants WHERE project_id!=? AND field=?",
                (project_id, CATALOG_FIELD),
            ).fetchone()
        if other is not None:
            raise CatalogTransportError("catalog_model_project_mismatch")
        raise CatalogTransportError("catalog_model_unauthorized")

    def _read_binding(self) -> dict:
        try:
            binding = self.worker({"op": "binding"})
        except CatalogTransportError:
            raise
        except (BridgeDisconnected, OSError, sqlite3.Error, ValueError, TypeError):
            raise CatalogTransportError("catalog_model_unbound") from None
        if not isinstance(binding, dict) or not isinstance(binding.get("model"), str):
            raise CatalogTransportError("catalog_model_unbound")
        return binding

    def _reserve_catalog_call(
        self,
        *,
        project_id: str,
        idempotency_key: str,
        projection_digest: str,
        clue: str,
        candidates: tuple[tuple[str, str], ...],
        binding: dict,
        identity: str,
    ) -> tuple[dict[str, object], str] | str:
        with closing(self.db()) as db, db:
            db.execute("BEGIN IMMEDIATE")
            replay = self._catalog_replay(
                project_id,
                idempotency_key,
                projection_digest,
                db=db,
            )
            if replay is not None:
                return replay
            rows = db.execute(
                "SELECT * FROM catalog_grants WHERE project_id=? AND field=? AND binding=?",
                (project_id, CATALOG_FIELD, identity),
            ).fetchall()
            sane = [row for row in rows if _catalog_grant_sane(row)]
            now = int(time.time())
            active = [
                row for row in sane if row["revoked"] == 0 and row["expires"] > now
            ]
            if len(active) > 1:
                raise CatalogTransportError("catalog_model_malformed")
            if not active:
                if any(row["revoked"] == 0 for row in sane):
                    raise CatalogTransportError("catalog_model_expired")
                if sane:
                    raise CatalogTransportError("catalog_model_revoked")
                raise CatalogTransportError("catalog_model_unauthorized")
            row = active[0]
            cap = int(row["output_cap"])
            request = classify_worker_request(
                binding=binding,
                clue=clue,
                candidates=candidates,
                output_tokens=cap,
            )
            if len(json.dumps(request).encode()) > MAX_CLASSIFY_REQUEST_BYTES:
                raise CatalogTransportError("catalog_model_malformed")
            amount = (
                row["input_rate"] * INPUT_RESERVE + row["output_rate"] * cap + 999999
            ) // 1000000
            updated = db.execute(
                """UPDATE catalog_grants SET requests=requests-1, input_left=input_left-?,
                output_left=output_left-?, minor_left=minor_left-?
                WHERE grant_id=? AND revoked=0 AND expires>? AND requests>=1
                AND input_left>=? AND output_left>=? AND minor_left>=?""",
                (
                    INPUT_RESERVE,
                    cap,
                    amount,
                    row["grant_id"],
                    now,
                    INPUT_RESERVE,
                    cap,
                    amount,
                ),
            )
            if updated.rowcount != 1:
                raise CatalogTransportError("catalog_model_quota_exhausted")
            grant_id = str(row["grant_id"])
            db.execute(
                """INSERT INTO catalog_reservations (
                    grant_id, call_key, projection_digest, model,
                    input_reserved, output_reserved, minor_reserved
                ) VALUES (?,?,?,?,?,?,?)""",
                (
                    grant_id,
                    idempotency_key,
                    projection_digest,
                    binding["model"],
                    INPUT_RESERVE,
                    cap,
                    amount,
                ),
            )
        return request, grant_id

    def _store_catalog_decision(
        self,
        *,
        grant_id: str,
        idempotency_key: str,
        projection_digest: str,
        decision: str,
    ) -> None:
        with closing(self.db()) as db, db:
            db.execute("BEGIN IMMEDIATE")
            updated = db.execute(
                """UPDATE catalog_reservations SET delivered=1, decision=?
                WHERE grant_id=? AND call_key=? AND projection_digest=? AND delivered=0""",
                (decision, grant_id, idempotency_key, projection_digest),
            )
            if updated.rowcount == 1:
                return
            stored = db.execute(
                """SELECT delivered, decision FROM catalog_reservations
                WHERE grant_id=? AND call_key=? AND projection_digest=?""",
                (grant_id, idempotency_key, projection_digest),
            ).fetchone()
            if (
                stored is not None
                and stored["delivered"] == 1
                and stored["decision"] == decision
            ):
                return
            raise CatalogTransportError("catalog_model_timeout")

    def _invoke_classify(self, request: dict[str, object]) -> dict:
        try:
            response = self.worker(request)
        except CatalogTransportError as exc:
            if exc.code == "catalog_model_timeout":
                raise
            raise CatalogTransportError("catalog_model_malformed") from None
        except (BridgeDisconnected, OSError, sqlite3.Error, ValueError, TypeError):
            raise CatalogTransportError("catalog_model_timeout") from None
        if not isinstance(response, dict):
            raise CatalogTransportError("catalog_model_malformed")
        return response


def _catalog_grant_sane(row: sqlite3.Row) -> bool:
    if not isinstance(row["currency"], str) or not re.fullmatch(
        "[A-Z]{3}", row["currency"]
    ):
        return False
    if row["field"] != CATALOG_FIELD or row["revoked"] not in {0, 1}:
        return False
    return (
        all(
            type(row[key]) is int and 0 <= row[key] <= _MAX
            for key in (
                "requests",
                "input_left",
                "output_left",
                "minor_left",
                "output_cap",
                "input_rate",
                "output_rate",
                "expires",
            )
        )
        and 1 <= row["output_cap"] <= 1024
        and row["input_rate"] >= 1
        and row["output_rate"] >= 1
    )


def configured_model(env: Mapping[str, str]) -> BoundModel | None:
    names = (
        "HELIXWEAVE_PERI_MODEL_LEDGER",
        "HELIXWEAVE_PERI_MODEL_WORKER",
        "HELIXWEAVE_PERI_MODEL_SETTINGS",
    )
    if not any(k in env for k in names):
        return None
    if not all(Path(env.get(k, "")).is_absolute() for k in names):
        raise ValueError("incomplete_bound_model_configuration")
    ledger, worker, settings = (Path(env[k]) for k in names)

    def invoke(request):
        try:
            result = subprocess.run(
                [str(worker), str(settings)],
                input=json.dumps(request),
                capture_output=True,
                text=True,
                check=True,
                timeout=35,
                env={},
                cwd="/",
            )
            reply = json.loads(result.stdout)
            if isinstance(reply, dict):
                return reply
        except (ValueError, OSError, subprocess.SubprocessError):
            pass
        raise BridgeDisconnected

    return BoundModel(ledger, invoke)
