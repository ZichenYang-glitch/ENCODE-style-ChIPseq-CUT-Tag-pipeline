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
        tokens = getattr(self.app.state, "test_auth_tokens", None)
        cookies = None
        if tokens is not None:
            policy = self.app.state.auth_cookie_policy
            cookies = httpx.Cookies()
            cookies.set(policy.session_cookie.name, tokens[0])
            cookies.set(policy.csrf_cookie.name, tokens[1])
            headers = dict(kwargs.pop("headers", {}) or {})
            headers.setdefault("X-CSRF-Token", tokens[1])
            kwargs["headers"] = headers
        async with httpx.AsyncClient(
            transport=transport,
            base_url="http://testserver",
            follow_redirects=True,
            cookies=cookies,
        ) as client:
            response = await client.request(method, url, **kwargs)
            await response.aread()
            return response


def seeded_auth_cookies(app: Any) -> httpx.Cookies | None:
    """Return the seeded auth cookie jar for raw AsyncClient constructions."""

    tokens = getattr(app.state, "test_auth_tokens", None)
    if tokens is None:
        return None
    policy = app.state.auth_cookie_policy
    cookies = httpx.Cookies()
    cookies.set(policy.session_cookie.name, tokens[0])
    cookies.set(policy.csrf_cookie.name, tokens[1])
    return cookies


def seeded_auth_async_client(app: Any, **kwargs: Any) -> httpx.AsyncClient:
    """Return an AsyncClient carrying the seeded session cookie jar and CSRF header."""

    tokens = getattr(app.state, "test_auth_tokens", None)
    if tokens is not None:
        kwargs.setdefault("cookies", seeded_auth_cookies(app))
        headers = dict(kwargs.pop("headers", {}) or {})
        headers.setdefault("X-CSRF-Token", tokens[1])
        kwargs["headers"] = headers
    return httpx.AsyncClient(**kwargs)
