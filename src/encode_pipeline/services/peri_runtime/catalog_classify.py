"""Bounded experiment.target classification for the owner-side Peri worker.

The production transport is the same one-call worker as numeric suggestions:
``integrations/peri-bound/helixweave-bound.rs`` invoked by
``configured_model``. This module does not open a provider, read a settings
file, or launch an agent. ``op=classify`` is a separate grant field from
draft tokens such as ``config.threads``.

A synthetic loopback may implement the same decision rule in tests. It is not
a stand-in for the owner worker.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import re


CATALOG_FIELD = "experiment.target"
CATALOG_REASON = "TARGET_TERM_MATCH"
INSUFFICIENT_EVIDENCE = "insufficient_evidence"
MODEL_PROPOSAL_SOURCE = "peri_bound_model"
MAX_CLASSIFY_CANDIDATES = 4
# Worker stdin is capped at 4096 bytes. The provider body stays at 1024 inside
# the Rust worker. 2048 leaves room for the binding fingerprint and four labels.
MAX_CLASSIFY_REQUEST_BYTES = 2048

_CLUE = re.compile(r"^[A-Za-z][A-Za-z0-9_.:-]{1,63}\Z")
_TOKEN = re.compile(r"^[a-z][a-z0-9]{0,31}\Z")
_TARGET_LABEL = re.compile(r"[A-Za-z][A-Za-z0-9]*(?:[ ._-][A-Za-z0-9]+)*\Z")
_ORDINAL_PREFIXES = (
    "number",
    "ordinal",
    "candidate",
    "index",
    "item",
    "num",
    "no",
)
_MESSAGES = {
    "catalog_model_unbound": (
        "The owner has no bound model for catalog classification."
    ),
    "catalog_model_unauthorized": (
        "No active catalog classification grant matches this project and binding."
    ),
    "catalog_model_revoked": "The catalog classification grant was revoked.",
    "catalog_model_expired": "The catalog classification grant has expired.",
    "catalog_model_project_mismatch": (
        "The catalog classification grant belongs to another project."
    ),
    "catalog_model_quota_exhausted": (
        "The catalog classification grant has no remaining calls."
    ),
    "catalog_model_timeout": (
        "The model call did not finish and its reservation was kept."
    ),
    "catalog_model_malformed": (
        "The model response was not a bound candidate or insufficient evidence."
    ),
    "catalog_idempotency_conflict": (
        "The idempotency key was already used for a different request."
    ),
}


class CatalogTransportError(Exception):
    """Stable failure on the catalog model path. The code is public."""

    def __init__(self, code: str) -> None:
        if code not in _MESSAGES:
            code = "catalog_model_malformed"
        self.code = code
        self.message = _MESSAGES[code]
        super().__init__(code)


@dataclass(frozen=True)
class BoundCandidate:
    """One vocabulary row bound locally. The model never sees term revisions."""

    token: str
    label: str
    term_id: str
    term_revision_id: str


def compact_clue(value: str) -> str:
    return "".join(character.lower() for character in value if character.isalnum())


def is_uninformative_clue(value: str) -> bool:
    """True when the clue is only a position, index, or empty ordinal."""

    compact = compact_clue(value)
    if not compact or compact.isdigit():
        return True
    for prefix in _ORDINAL_PREFIXES:
        if compact.startswith(prefix):
            rest = compact[len(prefix) :]
            if rest.isdigit():
                return True
    return False


def clue_is_auditable(value: object) -> bool:
    return (
        isinstance(value, str)
        and _CLUE.fullmatch(value) is not None
        and not is_uninformative_clue(value)
    )


def _projection_label(display_name: str, aliases: Sequence[str]) -> str | None:
    for raw in (display_name, *aliases):
        if isinstance(raw, str) and evidence_token(raw) is not None:
            return raw
    return None


def evidence_token(label: str) -> str | None:
    if len(label) > 48 or _TARGET_LABEL.fullmatch(label) is None:
        return None
    token = compact_clue(label)
    if _TOKEN.fullmatch(token) is None or is_uninformative_clue(token):
        return None
    return token


def bind_target_candidates(
    terms: Sequence[Mapping[str, object]],
    *,
    explicit_revision_id: str | None,
) -> tuple[BoundCandidate, ...]:
    """Map current target terms to unique evidence tokens. Indexes are not tokens."""

    bound: list[BoundCandidate] = []
    seen: dict[str, int] = {}
    for term in terms:
        if term.get("archived") is True:
            continue
        display_name = term.get("display_name")
        aliases = term.get("aliases") or ()
        if not isinstance(display_name, str) or isinstance(aliases, (str, bytes)):
            continue
        if not isinstance(aliases, Sequence):
            continue
        label = _projection_label(display_name, tuple(aliases))
        if label is None:
            continue
        token = evidence_token(label)
        term_id = term.get("term_id")
        revision_id = term.get("term_revision_id")
        if (
            token is None
            or not isinstance(term_id, str)
            or not isinstance(revision_id, str)
        ):
            continue
        seen[token] = seen.get(token, 0) + 1
        bound.append(
            BoundCandidate(
                token=token,
                label=label,
                term_id=term_id,
                term_revision_id=revision_id,
            )
        )
    unique = tuple(item for item in bound if seen[item.token] == 1)
    if explicit_revision_id is not None:
        chosen = tuple(
            item for item in unique if item.term_revision_id == explicit_revision_id
        )
        return chosen
    return unique


def classify_worker_request(
    *,
    binding: Mapping[str, object],
    clue: str,
    candidates: Sequence[tuple[str, str]],
    output_tokens: int,
) -> dict[str, object]:
    """The only object the classify worker is allowed to read."""

    return {
        "op": "classify",
        "binding": dict(binding),
        "field": CATALOG_FIELD,
        "reason": CATALOG_REASON,
        "clue": clue,
        "candidates": [{"token": token, "label": label} for token, label in candidates],
        "output_tokens": output_tokens,
    }


def accepted_decision(
    decision: object, candidates: Sequence[tuple[str, str]]
) -> str | None:
    """Accept insufficient evidence or one projected token. Never an index."""

    if decision == INSUFFICIENT_EVIDENCE:
        return INSUFFICIENT_EVIDENCE
    if (
        not isinstance(decision, str)
        or decision.isdigit()
        or is_uninformative_clue(decision)
    ):
        return None
    tokens = {token for token, _label in candidates}
    if decision in tokens:
        return decision
    return None


def loopback_decision(clue: str, candidates: Sequence[Mapping[str, str]]) -> str:
    """Synthetic provider rule: match a label or token, or decline.

    Position in the candidate list is not evidence. An ordinal clue declines
    even when candidates exist.
    """

    if is_uninformative_clue(clue):
        return INSUFFICIENT_EVIDENCE
    folded = clue.casefold()
    matches: list[str] = []
    for candidate in candidates:
        token = candidate.get("token")
        label = candidate.get("label")
        if not isinstance(token, str) or not isinstance(label, str):
            continue
        if folded in {token.casefold(), label.casefold()}:
            matches.append(token)
    if len(matches) == 1:
        return matches[0]
    return INSUFFICIENT_EVIDENCE
