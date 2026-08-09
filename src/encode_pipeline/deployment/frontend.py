"""Same-origin production application for verified precompiled frontend bytes."""

from __future__ import annotations

from collections.abc import Mapping
import os
from pathlib import Path, PurePosixPath
import re
from fastapi import FastAPI
from starlette.types import Receive, Scope, Send

from encode_pipeline.frontend_assets import (
    AssetRecord,
    VerifiedFrontendAssets,
    load_packaged_frontend_assets,
)


_SAFE_ROUTE_SEGMENT = re.compile(r"^[0-9A-Za-z][0-9A-Za-z._~:@+-]{0,255}$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
API_CONTRACT_ENVIRONMENT = "HELIXWEAVE_ACTIVE_API_CONTRACT_SHA256"
_CONTENT_TYPES = {
    ".css": "text/css; charset=utf-8",
    ".gif": "image/gif",
    ".html": "text/html; charset=utf-8",
    ".ico": "image/x-icon",
    ".jpeg": "image/jpeg",
    ".jpg": "image/jpeg",
    ".js": "text/javascript; charset=utf-8",
    ".json": "application/json; charset=utf-8",
    ".png": "image/png",
    ".svg": "image/svg+xml",
    ".webp": "image/webp",
    ".woff": "font/woff",
    ".woff2": "font/woff2",
}
_SECURITY_HEADERS = (
    (b"permissions-policy", b"camera=(), microphone=(), geolocation=()"),
    (b"referrer-policy", b"no-referrer"),
    (b"x-content-type-options", b"nosniff"),
    (b"x-frame-options", b"DENY"),
)


class ManifestStaticApplication:
    """Serve only manifest-enumerated bytes plus bounded SPA navigation routes."""

    def __init__(self, assets: VerifiedFrontendAssets) -> None:
        self.assets = assets
        self._records: Mapping[str, AssetRecord] = {
            item.path: item for item in assets.manifest.files
        }

    async def __call__(self, scope: Scope, receive: Receive, send: Send) -> None:
        if scope["type"] == "websocket":
            await send({"type": "websocket.close", "code": 1008})
            return
        if scope["type"] != "http":
            return
        method = str(scope.get("method", "")).upper()
        if method not in {"GET", "HEAD"}:
            await self._plain(
                send,
                405,
                b"Method not allowed.\n",
                extra_headers=((b"allow", b"GET, HEAD"),),
                head=method == "HEAD",
            )
            return

        requested = _request_asset_path(scope)
        if requested is None:
            await self._plain(send, 404, b"Not found.\n", head=method == "HEAD")
            return
        if requested not in self.assets.content:
            if not _spa_navigation(requested, scope):
                await self._plain(
                    send,
                    404,
                    b"Not found.\n",
                    head=method == "HEAD",
                )
                return
            requested = self.assets.manifest.entrypoint

        content = self.assets.content[requested]
        record = self._records[requested]
        etag = f'"{record.sha256}"'.encode("ascii")
        if _header(scope, b"if-none-match") == etag:
            await send(
                {
                    "type": "http.response.start",
                    "status": 304,
                    "headers": (
                        (b"etag", etag),
                        (b"cache-control", _cache_control(requested)),
                        *_SECURITY_HEADERS,
                    ),
                }
            )
            await send({"type": "http.response.body", "body": b""})
            return

        headers = (
            (b"content-type", _content_type(requested).encode("ascii")),
            (b"content-length", str(len(content)).encode("ascii")),
            (b"cache-control", _cache_control(requested)),
            (b"etag", etag),
            *_SECURITY_HEADERS,
        )
        await send({"type": "http.response.start", "status": 200, "headers": headers})
        await send(
            {
                "type": "http.response.body",
                "body": b"" if method == "HEAD" else content,
            }
        )

    async def _plain(
        self,
        send: Send,
        status_code: int,
        content: bytes,
        *,
        extra_headers: tuple[tuple[bytes, bytes], ...] = (),
        head: bool = False,
    ) -> None:
        headers = (
            (b"content-type", b"text/plain; charset=utf-8"),
            (b"content-length", str(len(content)).encode("ascii")),
            (b"cache-control", b"no-store"),
            *_SECURITY_HEADERS,
            *extra_headers,
        )
        await send(
            {
                "type": "http.response.start",
                "status": status_code,
                "headers": headers,
            }
        )
        await send(
            {
                "type": "http.response.body",
                "body": b"" if head else content,
            }
        )


def _header(scope: Scope, name: bytes) -> bytes | None:
    for observed_name, value in scope.get("headers", ()):
        if observed_name.lower() == name:
            return value
    return None


def _request_asset_path(scope: Scope) -> str | None:
    path = scope.get("path")
    if not isinstance(path, str) or not path.startswith("/"):
        return None
    if any(ord(character) < 32 or ord(character) == 127 for character in path):
        return None
    relative = path[1:]
    if relative == "":
        return "index.html"
    if (
        relative.startswith("/")
        or relative.endswith("/")
        or "//" in relative
        or "\\" in relative
    ):
        return None
    parts = PurePosixPath(relative).parts
    if any(part in {"", ".", ".."} for part in parts):
        return None
    return relative


def _spa_navigation(path: str, scope: Scope) -> bool:
    accept = _header(scope, b"accept")
    if accept is not None and b"text/html" not in accept and b"*/*" not in accept:
        return False
    parts = path.split("/")
    if any(_SAFE_ROUTE_SEGMENT.fullmatch(part) is None for part in parts):
        return False
    if parts == ["workflows"] or parts == ["runs"] or parts == ["artifacts"]:
        return True
    if parts[0] == "workflows":
        return len(parts) == 2 or (len(parts) == 3 and parts[2] == "new-run")
    if parts[0] == "runs":
        return len(parts) == 2
    if parts[0] == "artifacts":
        return len(parts) == 3
    return False


def _content_type(path: str) -> str:
    return _CONTENT_TYPES.get(
        PurePosixPath(path).suffix.lower(), "application/octet-stream"
    )


def _cache_control(path: str) -> bytes:
    if path.startswith("assets/"):
        return b"public, max-age=31536000, immutable"
    if path == "index.html":
        return b"no-cache"
    return b"public, max-age=3600"


def create_app(
    *,
    database_url: str | None = None,
    workspace_root: Path | None = None,
    project_root: Path | None = None,
    verified_assets: VerifiedFrontendAssets | None = None,
    api_contract_sha256: str | None = None,
) -> FastAPI:
    """Create the API and verified frontend at one origin.

    Package and release-bound API compatibility verification happen before the
    API factory can open or migrate its database.  The operator materializes the
    expected digest from the verified active platform release; callers cannot
    derive it by opening persistence first.
    """
    assets = verified_assets or load_packaged_frontend_assets()
    observed_contract = (
        os.environ.get(API_CONTRACT_ENVIRONMENT)
        if api_contract_sha256 is None
        else api_contract_sha256
    )
    if (
        not isinstance(observed_contract, str)
        or _SHA256.fullmatch(observed_contract) is None
        or observed_contract != assets.manifest.api_contract_sha256
    ):
        from encode_pipeline.deployment.errors import fail

        raise fail(
            "FRONTEND_API_INCOMPATIBLE",
            "Frontend and API contracts are not compatible.",
        )

    from encode_pipeline.api.main import create_app as create_api_app

    app = create_api_app(
        database_url=database_url,
        workspace_root=workspace_root,
        project_root=project_root,
    )
    app.state.frontend_asset_identity = assets.manifest.identity
    app.state.frontend_api_contract_sha256 = observed_contract
    app.mount("/", ManifestStaticApplication(assets), name="frontend")
    return app


__all__ = [
    "API_CONTRACT_ENVIRONMENT",
    "ManifestStaticApplication",
    "create_app",
]
