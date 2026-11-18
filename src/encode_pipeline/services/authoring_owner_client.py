"""API-side connection to the single authoring owner.

The API process does not import a Peri session or an ACP payload. It either
keeps the explicit disabled bridge or connects to the private Unix transport.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

from encode_pipeline.services.authoring_bridge import (
    AuthoringBridge,
    UnavailableAuthoringBridge,
)


PERI_OWNER_SOCKET_ENV = "HELIXWEAVE_PERI_OWNER_SOCKET"


def load_authoring_bridge(environ: Mapping[str, str]) -> AuthoringBridge:
    """Return the disabled bridge, or the private owner client."""

    raw_socket = environ.get(PERI_OWNER_SOCKET_ENV, "")
    if not raw_socket:
        return UnavailableAuthoringBridge()
    from encode_pipeline.services.peri_runtime.socket_server import (
        RemoteAuthoringBridge,
    )

    return RemoteAuthoringBridge(Path(raw_socket))
