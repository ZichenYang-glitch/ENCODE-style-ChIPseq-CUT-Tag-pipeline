"""API-side catalog classification client.

The API process does not open Peri settings or a provider key. When the owner
socket is configured, classification goes to that owner, which already runs
the bound worker. Tests can pass a ledger model instead.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

from encode_pipeline.services.authoring_bridge import BridgeDisconnected
from encode_pipeline.services.authoring_owner_client import PERI_OWNER_SOCKET_ENV
from encode_pipeline.services.peri_runtime.bound_model import BoundModel
from encode_pipeline.services.peri_runtime.catalog_classify import (
    CatalogTransportError,
)


class LedgerCatalogModel:
    """In-process ledger and worker. Used by tests and by the owner itself."""

    def __init__(self, model: BoundModel) -> None:
        self._model = model

    def classify(
        self,
        *,
        project_id: str,
        idempotency_key: str,
        clue: str,
        candidates: Sequence[tuple[str, str]],
        projection_digest: str,
    ) -> str:
        return self._model.classify(
            project_id=project_id,
            idempotency_key=idempotency_key,
            clue=clue,
            candidates=candidates,
            projection_digest=projection_digest,
        )


class OwnerCatalogModel:
    """Forward one classification to the live owner socket."""

    def __init__(self, socket_path: Path) -> None:
        self._path = socket_path

    def classify(
        self,
        *,
        project_id: str,
        idempotency_key: str,
        clue: str,
        candidates: Sequence[tuple[str, str]],
        projection_digest: str,
    ) -> str:
        from encode_pipeline.services.peri_runtime.socket_server import (
            RemoteAuthoringBridge,
        )

        bridge = RemoteAuthoringBridge(self._path)
        try:
            return bridge.classify_catalog(
                project_id=project_id,
                idempotency_key=idempotency_key,
                clue=clue,
                candidates=[
                    {"token": token, "label": label} for token, label in candidates
                ],
                projection_digest=projection_digest,
            )
        except BridgeDisconnected:
            raise CatalogTransportError("catalog_model_unbound") from None


def load_catalog_model(environ: Mapping[str, str]) -> OwnerCatalogModel | None:
    """Return the owner client only when the existing owner socket is configured."""

    raw_socket = environ.get(PERI_OWNER_SOCKET_ENV, "")
    if not raw_socket:
        return None
    return OwnerCatalogModel(Path(raw_socket))
