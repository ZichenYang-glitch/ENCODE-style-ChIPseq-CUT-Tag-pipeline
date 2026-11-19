"""Sign browser confirmations with private, generation-scoped host material."""

from __future__ import annotations

import hmac
import json
from collections.abc import Mapping
from time import time

from encode_pipeline.persistence.authoring import AuthoringSessionRecord
from encode_pipeline.platform.authoring import (
    CONNECTION_NOT_READY,
    AuthoringRejected,
)


_CONFIRMATION_SECONDS = 10 * 60


def issue_confirmation(
    connection: AuthoringSessionRecord, target: Mapping[str, object]
) -> str:
    expires_at = int(time()) + _CONFIRMATION_SECONDS
    return f"{expires_at}.{_signature(connection, target, expires_at)}"


def verified_confirmation_expiry(
    connection: AuthoringSessionRecord,
    target: Mapping[str, object],
    confirmation: str,
) -> int | None:
    expires, signature = confirmation.split(".", 1)
    expires_at = int(expires)
    if time() < expires_at and hmac.compare_digest(
        signature, _signature(connection, target, expires_at)
    ):
        return expires_at
    return None


def _signature(
    connection: AuthoringSessionRecord, target: Mapping[str, object], expires_at: int
) -> str:
    # The capability digest is stored only on the platform. The underlying
    # capability is held by the trusted owner, outside the Peri process.
    if connection.connection_capability_digest is None:
        raise AuthoringRejected(CONNECTION_NOT_READY)
    payload = json.dumps(
        {
            "purpose": "helixweave.peri.browser.run-confirmation.v1",
            "user": connection.owner_user_id,
            "login": connection.originating_auth_session_ref,
            "sid": connection.sid,
            "generation": connection.generation,
            "instance": connection.instance_id,
            "owner_epoch": connection.bridge_owner_epoch,
            "workflow": connection.workflow_id,
            "expires_at": expires_at,
            "target": dict(target),
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    return hmac.digest(
        bytes.fromhex(connection.connection_capability_digest),
        payload.encode("ascii"),
        "sha256",
    ).hex()
