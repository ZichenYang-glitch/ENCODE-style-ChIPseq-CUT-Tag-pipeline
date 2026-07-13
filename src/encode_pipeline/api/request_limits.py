"""Bounded pre-parser request handling for workflow authoring endpoints."""

from __future__ import annotations

import re

from starlette.responses import JSONResponse
from starlette.types import ASGIApp, Message, Receive, Scope, Send

from encode_pipeline.api.models import IssueResponse, RunResponse, ValidationResponse
from encode_pipeline.platform.adapters import MAX_AUTHORING_REQUEST_BYTES


_AUTHORING_PATH = re.compile(
    r"^/api/v1/workflows/(?P<workflow_id>[^/]+)/(?P<operation>validate|runs)/?$"
)


class AuthoringRequestLimitMiddleware:
    """Buffer and bound validate/create bodies before JSON parsing."""

    def __init__(
        self,
        app: ASGIApp,
        *,
        max_request_bytes: int = MAX_AUTHORING_REQUEST_BYTES,
    ) -> None:
        if (
            isinstance(max_request_bytes, bool)
            or not isinstance(max_request_bytes, int)
            or max_request_bytes <= 0
        ):
            raise ValueError("max_request_bytes must be a positive integer")
        self._app = app
        self._max_request_bytes = max_request_bytes

    async def __call__(self, scope: Scope, receive: Receive, send: Send) -> None:
        target = self._target(scope)
        if target is None:
            await self._app(scope, receive, send)
            return

        declared_length = self._content_length(scope)
        if declared_length is not None and declared_length > self._max_request_bytes:
            await self._send_too_large(scope, receive, send, target)
            return

        body = bytearray()
        while True:
            message = await receive()
            if message["type"] == "http.disconnect":
                return
            if message["type"] != "http.request":
                continue
            chunk = message.get("body", b"")
            if len(body) + len(chunk) > self._max_request_bytes:
                await self._send_too_large(scope, receive, send, target)
                return
            body.extend(chunk)
            if not message.get("more_body", False):
                break

        replayed = False

        async def replay_receive() -> Message:
            nonlocal replayed
            if not replayed:
                replayed = True
                return {
                    "type": "http.request",
                    "body": bytes(body),
                    "more_body": False,
                }
            return await receive()

        await self._app(scope, replay_receive, send)

    @staticmethod
    def _target(scope: Scope) -> tuple[str, str] | None:
        if scope.get("type") != "http" or scope.get("method") != "POST":
            return None
        match = _AUTHORING_PATH.fullmatch(str(scope.get("path", "")))
        if match is None:
            return None
        return match.group("workflow_id"), match.group("operation")

    @staticmethod
    def _content_length(scope: Scope) -> int | None:
        for name, value in scope.get("headers", ()):
            if name.lower() != b"content-length":
                continue
            try:
                parsed = int(value.decode("ascii"), 10)
            except (UnicodeDecodeError, ValueError):
                return None
            return parsed if parsed >= 0 else None
        return None

    async def _send_too_large(
        self,
        scope: Scope,
        receive: Receive,
        send: Send,
        target: tuple[str, str],
    ) -> None:
        workflow_id, operation = target
        issue = IssueResponse(
            code="API_REQUEST_TOO_LARGE",
            message="Request body exceeds the workflow authoring limit.",
            severity="error",
            path="body",
            source="api",
            technical_message=None,
            hint="Reduce the config or sample rows and try again.",
            context={"max_request_bytes": self._max_request_bytes},
        )
        if operation == "runs":
            content = RunResponse(ok=False, run=None, issues=[issue]).model_dump(
                mode="json"
            )
        else:
            content = ValidationResponse(
                ok=False,
                workflow_id=workflow_id,
                value=None,
                issues=[issue],
            ).model_dump(mode="json")
        response = JSONResponse(
            status_code=413,
            content=content,
            headers={"connection": "close"},
        )
        await response(scope, receive, send)
