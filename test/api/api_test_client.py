"""Small ASGI test client for API route tests.

Starlette's TestClient currently prefers httpx2, which hangs in this local
Python 3.13 environment. This helper keeps route tests synchronous while using
httpx's ASGI transport directly.
"""

from __future__ import annotations

import asyncio
from types import TracebackType
from typing import Any, Self

import httpx


class ApiTestClient:
    """Synchronous wrapper around httpx.AsyncClient and ASGITransport."""

    def __init__(self, app: Any) -> None:
        self.app = app

    def __enter__(self) -> Self:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        return None

    def get(self, url: str, **kwargs: Any) -> httpx.Response:
        return self.request("GET", url, **kwargs)

    def post(self, url: str, **kwargs: Any) -> httpx.Response:
        return self.request("POST", url, **kwargs)

    def request(self, method: str, url: str, **kwargs: Any) -> httpx.Response:
        return asyncio.run(self._request(method, url, **kwargs))

    async def _request(self, method: str, url: str, **kwargs: Any) -> httpx.Response:
        transport = httpx.ASGITransport(app=self.app)
        async with httpx.AsyncClient(
            transport=transport,
            base_url="http://testserver",
            follow_redirects=True,
        ) as client:
            response = await client.request(method, url, **kwargs)
            await response.aread()
            return response
