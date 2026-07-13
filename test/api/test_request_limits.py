"""Tests for the pre-parser workflow authoring request limit."""

from __future__ import annotations

import asyncio
import json
from collections.abc import Sequence

import pytest

from encode_pipeline.api.request_limits import AuthoringRequestLimitMiddleware
from encode_pipeline.platform.adapters import MAX_AUTHORING_REQUEST_BYTES


def _scope(path: str, headers: Sequence[tuple[bytes, bytes]] = ()) -> dict:
    return {
        "type": "http",
        "asgi": {"version": "3.0"},
        "http_version": "1.1",
        "method": "POST",
        "scheme": "http",
        "path": path,
        "raw_path": path.encode("ascii"),
        "query_string": b"",
        "headers": list(headers),
        "client": ("127.0.0.1", 1),
        "server": ("testserver", 80),
        "root_path": "",
    }


async def _invoke(
    *,
    path: str,
    chunks: list[bytes],
    headers: Sequence[tuple[bytes, bytes]] = (),
) -> tuple[int, list[dict]]:
    downstream_calls = 0
    sent: list[dict] = []
    messages = [
        {
            "type": "http.request",
            "body": chunk,
            "more_body": index < len(chunks) - 1,
        }
        for index, chunk in enumerate(chunks)
    ]
    messages.append({"type": "http.disconnect"})

    async def receive() -> dict:
        return messages.pop(0)

    async def send(message: dict) -> None:
        sent.append(message)

    async def downstream(scope, receive, send) -> None:
        nonlocal downstream_calls
        downstream_calls += 1
        received = bytearray()
        while True:
            message = await receive()
            if message["type"] != "http.request":
                break
            received.extend(message.get("body", b""))
            if not message.get("more_body", False):
                break
        await send({"type": "http.response.start", "status": 204, "headers": []})
        await send({"type": "http.response.body", "body": bytes(received)})

    middleware = AuthoringRequestLimitMiddleware(downstream)
    await middleware(_scope(path, headers), receive, send)
    return downstream_calls, sent


def _response(sent: list[dict]) -> tuple[int, dict]:
    start = next(
        message for message in sent if message["type"] == "http.response.start"
    )
    body = b"".join(
        message.get("body", b"")
        for message in sent
        if message["type"] == "http.response.body"
    )
    return start["status"], json.loads(body)


def _status(sent: list[dict]) -> int:
    return next(
        message["status"]
        for message in sent
        if message["type"] == "http.response.start"
    )


def test_exact_limit_with_multiple_chunks_reaches_downstream() -> None:
    calls, sent = asyncio.run(
        _invoke(
            path="/api/v1/workflows/encode/validate",
            chunks=[b"a" * (MAX_AUTHORING_REQUEST_BYTES // 2)] * 2,
        )
    )

    assert calls == 1
    assert _status(sent) == 204


def test_many_empty_chunks_do_not_change_the_bounded_replayed_body() -> None:
    calls, sent = asyncio.run(
        _invoke(
            path="/api/v1/workflows/encode/validate",
            chunks=([b""] * 100) + [b"{}"],
        )
    )

    assert calls == 1
    assert _response(sent) == (204, {})


@pytest.mark.parametrize("operation", ["validate", "runs"])
def test_actual_bytes_over_limit_are_rejected_without_content_length(
    operation: str,
) -> None:
    calls, sent = asyncio.run(
        _invoke(
            path=f"/api/v1/workflows/encode/{operation}",
            chunks=[b"a" * MAX_AUTHORING_REQUEST_BYTES, b"b"],
        )
    )

    status, body = _response(sent)
    assert calls == 0
    assert status == 413
    assert body["ok"] is False
    assert body["issues"][0]["code"] == "API_REQUEST_TOO_LARGE"
    assert body["issues"][0]["technical_message"] is None
    if operation == "validate":
        assert body["workflow_id"] == "encode"
        assert body["value"] is None
    else:
        assert body["run"] is None


def test_falsely_small_content_length_does_not_bypass_actual_byte_limit() -> None:
    calls, sent = asyncio.run(
        _invoke(
            path="/api/v1/workflows/encode/validate",
            chunks=[b"a" * MAX_AUTHORING_REQUEST_BYTES, b"b"],
            headers=[(b"content-length", b"1")],
        )
    )

    assert calls == 0
    assert _response(sent)[0] == 413


def test_declared_oversize_is_rejected_before_body_is_read() -> None:
    calls, sent = asyncio.run(
        _invoke(
            path="/api/v1/workflows/encode/runs",
            chunks=[b"{}"],
            headers=[
                (b"content-length", str(MAX_AUTHORING_REQUEST_BYTES + 1).encode())
            ],
        )
    )

    assert calls == 0
    assert _response(sent)[0] == 413


def test_non_authoring_route_is_not_limited() -> None:
    calls, sent = asyncio.run(
        _invoke(
            path="/api/v1/runs/run-1/start",
            chunks=[b"a" * (MAX_AUTHORING_REQUEST_BYTES + 1)],
        )
    )

    assert calls == 1
    assert _status(sent) == 204
