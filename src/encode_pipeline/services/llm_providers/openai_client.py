"""OpenAI-backed implementation of the LLMClient protocol.

This module is import-safe without the ``openai`` SDK installed. The SDK is
only imported inside :meth:`OpenAILLMClient.__init__` when no client is
injected, which makes the provider constructible in tests with a fake client.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any

from encode_pipeline.services.llm_client import LLMMessage, LLMResponse, ToolDefinition
from encode_pipeline.services.llm_factory import LLMProviderError

if TYPE_CHECKING:
    from openai import AsyncOpenAI  # type: ignore[import-not-found]


class OpenAILLMClient:
    """LLMClient implementation backed by the OpenAI API.

    Args:
        api_key: OpenAI API key. Never exposed by ``__repr__`` or ``__str__``.
        model: Model identifier passed to ``chat.completions.create``.
            Defaults to ``"gpt-4o-mini"``.
        client: Optional injected client. When provided, the real OpenAI SDK
            is not imported and the injected object is used instead. This is
            used by tests to avoid network calls and SDK dependencies.
    """

    def __init__(
        self,
        *,
        api_key: str,
        model: str = "gpt-4o-mini",
        client: Any | None = None,
    ) -> None:
        self._api_key = api_key
        self._model = model
        if client is not None:
            self._client: AsyncOpenAI = client
        else:
            import openai

            self._client = openai.AsyncOpenAI(api_key=api_key)

    async def complete(
        self,
        messages: list[LLMMessage],
        tools: list[ToolDefinition],  # noqa: ARG002
    ) -> LLMResponse:
        """Request a chat completion from OpenAI.

        Args:
            messages: Conversation history as provider-neutral messages.
            tools: Tool definitions (currently ignored; reserved for future use).

        Returns:
            Provider-neutral response.

        Raises:
            LLMProviderError: If the underlying SDK request fails. The API key
                is redacted from the exception message.
        """
        try:
            raw_response = await self._client.chat.completions.create(
                model=self._model,
                messages=[{"role": message.role, "content": message.content} for message in messages],
            )
        except Exception as exc:
            raise LLMProviderError(_redact(str(exc), self._api_key)) from exc

        content = raw_response.choices[0].message.content or ""
        return LLMResponse(content=_redact(content, self._api_key))

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(model={self._model!r})"

    def __str__(self) -> str:
        return f"OpenAILLMClient(model={self._model})"


def _redact(text: str, secret: str) -> str:
    """Return *text* with *secret* replaced by a fixed mask."""
    if not secret:
        return text
    return text.replace(secret, "***REDACTED***")
