"""Read-only classification previews against the catalog HTTP API."""

from __future__ import annotations

import argparse
from collections import Counter
from collections.abc import Iterator, Sequence
import csv
from dataclasses import dataclass
from http.client import HTTPException
import ipaddress
import json
import os
from pathlib import Path
import re
import sys
from urllib.error import HTTPError
from urllib.parse import urlencode, urlsplit
from urllib.request import HTTPRedirectHandler, ProxyHandler, Request, build_opener


_SESSION_ENV = "HELIXWEAVE_CATALOG_SESSION_TOKEN"
_PAGE_SIZE = 100
_TSV_COLUMNS = (
    "experiment_id",
    "expected_experiment_revision_id",
    "assay_term_revision_id",
)


@dataclass(frozen=True)
class _Proposal:
    line: int
    experiment_id: str
    expected_experiment_revision_id: str
    assay_term_revision_id: str


class _ReadError(Exception):
    pass


class _NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, *args, **kwargs):
        return None


def _is_id(value: object, prefix: str) -> bool:
    return (
        isinstance(value, str)
        and re.fullmatch(f"{prefix}_[0-9a-f]{{32}}", value) is not None
    )


def _read_proposals(path: Path) -> list[_Proposal]:
    try:
        with path.open(encoding="utf-8-sig", newline="") as source:
            table = csv.reader(source, delimiter="\t", quoting=csv.QUOTE_NONE)
            if tuple(next(table, ())) != _TSV_COLUMNS:
                raise ValueError(
                    "invalid TSV header; expected " + ", ".join(_TSV_COLUMNS)
                )
            proposals = []
            for line, cells in enumerate(table, start=2):
                if len(cells) != len(_TSV_COLUMNS) or not all(
                    _is_id(value, prefix)
                    for value, prefix in zip(cells, ("exp", "expr", "trmv"))
                ):
                    raise ValueError(
                        f"invalid TSV row at line {line}; expected three catalog IDs"
                    )
                proposals.append(_Proposal(line, *cells))
    except (OSError, UnicodeError, csv.Error):
        raise ValueError("could not read a UTF-8 TSV") from None
    if not proposals:
        raise ValueError("TSV has no classification rows")
    return proposals


def _api_origin(value: str) -> tuple[str, str]:
    try:
        parsed = urlsplit(value)
        host = parsed.hostname
        port = parsed.port
        if (
            any(char.isspace() or ord(char) < 32 for char in value)
            or any(char in value for char in "?#\\")
            or not host
            or parsed.username is not None
            or parsed.path not in ("", "/")
            or port == 0
        ):
            raise ValueError
        authority = f"[{host}]" if ":" in host else host
        if port is not None:
            authority = f"{authority}:{port}"
        if parsed.netloc.lower() != authority.lower():
            raise ValueError
        if parsed.scheme == "https":
            cookie_name = "__Host-helixweave_session"
        elif parsed.scheme == "http" and ipaddress.ip_address(host).is_loopback:
            cookie_name = "helixweave_session"
        else:
            raise ValueError
    except ValueError:
        raise ValueError(
            "--api-origin must be an HTTPS origin or a numeric loopback HTTP origin"
        ) from None
    return f"{parsed.scheme}://{authority}", cookie_name


class _CatalogReader:
    def __init__(self, origin: str, cookie_name: str, token: str) -> None:
        self._origin = origin
        self._headers = {"Accept": "application/json"}
        if token:
            self._headers["Cookie"] = f"{cookie_name}={token}"
        self._opener = build_opener(ProxyHandler({}), _NoRedirect())

    def _get(self, path: str, params: dict[str, str]) -> dict:
        try:
            request = Request(
                f"{self._origin}/api/v1/catalog{path}?{urlencode(params)}",
                headers=self._headers,
                method="GET",
            )
            with self._opener.open(request, timeout=10) as response:
                page = json.load(response)
        except HTTPError as error:
            error.close()
            raise _ReadError(f"catalog GET failed (HTTP {error.code})") from None
        except (OSError, HTTPException):
            raise _ReadError("catalog GET failed (transport error)") from None
        except (ValueError, RecursionError):
            raise _ReadError("invalid catalog response") from None
        if not isinstance(page, dict):
            raise _ReadError("invalid catalog response")
        return page

    def _rows(
        self, path: str, field: str, *, project_id: str | None = None
    ) -> Iterator[dict]:
        params = {"revision_view": "linked", "limit": str(_PAGE_SIZE)}
        seen_cursors: set[str] = set()
        while True:
            page = self._get(path, params)
            if project_id is not None and (
                page.get("ok") is not True
                or page.get("project_id") != project_id
                or page.get("revision_view") != "linked"
            ):
                raise _ReadError("invalid catalog project response")
            items = page.get(field)
            if not isinstance(items, list) or not all(
                isinstance(item, dict) for item in items
            ):
                raise _ReadError("invalid catalog list")
            cursor = page.get("next_cursor")
            if cursor is not None:
                if not isinstance(cursor, str) or not cursor or cursor in seen_cursors:
                    raise _ReadError("invalid or repeated catalog cursor")
                seen_cursors.add(cursor)
            yield from items
            if cursor is None:
                return
            params["cursor"] = cursor

    def current(
        self, project_id: str, wanted: set[str]
    ) -> dict[str, tuple[str, str | None]]:
        current = {}
        for item in self._rows(
            f"/projects/{project_id}/experiments", "experiments", project_id=project_id
        ):
            experiment_id = item.get("experiment_id")
            if not _is_id(experiment_id, "exp"):
                raise _ReadError("invalid catalog experiment ID")
            if experiment_id not in wanted:
                continue
            revision = item.get("experiment_revision_id")
            assay = item.get("assay")
            if (
                experiment_id in current
                or not _is_id(revision, "expr")
                or not isinstance(assay, dict)
                or not isinstance(assay.get("status"), str)
            ):
                raise _ReadError("invalid or repeated catalog experiment")
            assay_revision = None
            if assay["status"] == "known":
                assay_revision = assay.get("term_revision_id")
                if not _is_id(assay_revision, "trmv"):
                    raise _ReadError("invalid catalog assay revision ID")
            current[experiment_id] = revision, assay_revision
        return current

    def has_revision(self, experiment_id: str, revision_id: str) -> bool:
        for item in self._rows(
            f"/experiments/{experiment_id}/revisions", "experiment_revisions"
        ):
            revision = item.get("experiment_revision_id")
            if item.get("experiment_id") != experiment_id or not _is_id(
                revision, "expr"
            ):
                raise _ReadError("invalid catalog experiment history")
            if revision == revision_id:
                return True
        return False


def _preview(
    reader: _CatalogReader, project_id: str, proposals: list[_Proposal]
) -> dict:
    counts = Counter(proposal.experiment_id for proposal in proposals)
    current = reader.current(
        project_id, {key for key, count in counts.items() if count == 1}
    )
    rows = []
    for proposal in proposals:
        experiment_id = proposal.experiment_id
        expected = proposal.expected_experiment_revision_id
        row = {
            "line": proposal.line,
            "experiment_id": experiment_id,
            "expected_experiment_revision_id": expected,
            "assay_term_revision_id": proposal.assay_term_revision_id,
        }
        if counts[experiment_id] > 1:
            row["status"] = "duplicate"
        elif experiment_id not in current:
            row["status"] = "unmatched"
        else:
            revision, assay_revision = current[experiment_id]
            row["current_experiment_revision_id"] = revision
            row["current_assay_term_revision_id"] = assay_revision
            if revision == expected:
                row["status"] = "matched"
            elif reader.has_revision(experiment_id, expected):
                row["status"] = "stale_revision"
            else:
                row["status"] = "unknown_revision"
        rows.append(row)
    return {"project_id": project_id, "rows": rows}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="helixweave catalog")
    subcommands = parser.add_subparsers(dest="command", required=True)
    preview = subcommands.add_parser(
        "preview", help="preview one assay classification TSV"
    )
    preview.add_argument("--api-origin", required=True)
    preview.add_argument("--project-id", required=True)
    preview.add_argument("tsv", type=Path)
    args = parser.parse_args(argv)

    try:
        origin, cookie_name = _api_origin(args.api_origin)
        if not _is_id(args.project_id, "prj"):
            raise ValueError("--project-id must be a catalog project ID")
        token = os.environ.get(_SESSION_ENV, "")
        if token and not re.fullmatch(r"[A-Za-z0-9_-]{43}", token):
            raise ValueError(
                f"{_SESSION_ENV} must contain only the session cookie value"
            )
        proposals = _read_proposals(args.tsv)
    except ValueError as error:
        preview.error(str(error))

    try:
        report = _preview(
            _CatalogReader(origin, cookie_name, token), args.project_id, proposals
        )
    except _ReadError as error:
        print(f"catalog preview: {error}", file=sys.stderr)
        return 1
    print(json.dumps(report, indent=2))
    return 0
