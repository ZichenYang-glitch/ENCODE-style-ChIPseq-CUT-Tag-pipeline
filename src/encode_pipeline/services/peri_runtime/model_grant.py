"""Local operator CLI. No secret values are accepted as arguments."""

import argparse
from contextlib import closing
import json
import os
from pathlib import Path
import sqlite3

from .bound_model import BoundModel, Budget, binding_id, configured_model


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("init")
    commands.add_parser("binding")
    g = commands.add_parser("grant")
    g.add_argument("--platform-db", required=True, type=Path)
    for name in ("sid", "field", "binding", "currency"):
        g.add_argument("--" + name, required=True)
    for name in (
        "expires",
        "requests",
        "input-tokens",
        "output-tokens",
        "minor-units",
        "output-per-call",
        "input-rate",
        "output-rate",
    ):
        g.add_argument("--" + name, required=True, type=int)
    commands.add_parser("revoke").add_argument("--instance", required=True)
    catalog = commands.add_parser("catalog-grant")
    for name in ("project", "field", "binding", "currency"):
        catalog.add_argument("--" + name, required=True)
    for name in (
        "expires",
        "requests",
        "input-tokens",
        "output-tokens",
        "minor-units",
        "output-per-call",
        "input-rate",
        "output-rate",
    ):
        catalog.add_argument("--" + name, required=True, type=int)
    commands.add_parser("catalog-revoke").add_argument("--grant", required=True)
    args = parser.parse_args()
    model = configured_model(os.environ)
    if model is None:
        parser.error("set all three HELIXWEAVE_PERI_MODEL_* paths")
    if args.command == "init":
        BoundModel.initialize(model.ledger)
    elif args.command == "binding":
        b = model.worker({"op": "binding"})
        print(json.dumps({"model": b["model"], "binding": binding_id(b)}))
    elif args.command == "revoke":
        with closing(model.db()) as db, db:
            if (
                db.execute(
                    "UPDATE grants SET revoked=1 WHERE instance=?", (args.instance,)
                ).rowcount
                != 1
            ):
                parser.error("unknown instance")
    elif args.command == "catalog-revoke":
        try:
            model.catalog_revoke(args.grant)
        except ValueError:
            parser.error("unknown catalog grant")
    elif args.command == "catalog-grant":
        if args.binding != binding_id(model.worker({"op": "binding"})):
            parser.error("binding changed")
        if args.field != "experiment.target":
            parser.error("catalog grants only authorize experiment.target")
        try:
            grant_id = model.catalog_grant(
                project_id=args.project,
                field=args.field,
                binding=args.binding,
                expires=args.expires,
                budget=Budget(
                    args.requests,
                    args.input_tokens,
                    args.output_tokens,
                    args.minor_units,
                    args.output_per_call,
                    args.input_rate,
                    args.output_rate,
                    args.currency,
                ),
            )
        except ValueError:
            parser.error("invalid catalog grant")
        print(
            json.dumps(
                {
                    "grant_id": grant_id,
                    "field": args.field,
                    "project_id": args.project,
                }
            )
        )
    else:
        if args.binding != binding_id(model.worker({"op": "binding"})):
            parser.error("binding changed")
        with closing(
            sqlite3.connect(args.platform_db.resolve().as_uri() + "?mode=ro", uri=True)
        ) as db:
            row = db.execute(
                """SELECT instance_id,generation,bridge_owner_epoch FROM agent_authoring_sessions
                WHERE sid=? AND connection_state='ready'""",
                (args.sid,),
            ).fetchone()
        if not row:
            parser.error("connect the draft first")
        model.grant(
            sid=args.sid,
            instance=row[0],
            generation=row[1],
            epoch=row[2],
            field=args.field,
            binding=args.binding,
            expires=args.expires,
            budget=Budget(
                args.requests,
                args.input_tokens,
                args.output_tokens,
                args.minor_units,
                args.output_per_call,
                args.input_rate,
                args.output_rate,
                args.currency,
            ),
        )
        print(json.dumps({"instance": row[0], "field": args.field}))


if __name__ == "__main__":
    main()
