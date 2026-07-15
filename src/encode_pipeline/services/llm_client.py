"""Provider-neutral LLM client abstraction and deterministic mock."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Protocol


@dataclass(frozen=True)
class LLMMessage:
    role: str
    content: str


@dataclass(frozen=True)
class ToolDefinition:
    name: str
    description: str
    parameters: dict[str, Any]


@dataclass(frozen=True)
class ToolCall:
    name: str
    arguments: dict[str, Any]


@dataclass(frozen=True)
class LLMResponse:
    content: str
    tool_calls: tuple[ToolCall, ...] = ()


class LLMClient(Protocol):
    async def complete(
        self,
        messages: list[LLMMessage],
        tools: list[ToolDefinition],
    ) -> LLMResponse: ...


class MockLLMClient:
    """Deterministic LLM client for tests and early integration."""

    def __init__(
        self, response_text: str = "", tool_calls: list[ToolCall] | None = None
    ):
        self._response_text = response_text
        self._tool_calls = tuple(tool_calls or ())

    async def complete(
        self,
        messages: list[LLMMessage],
        tools: list[ToolDefinition],
    ) -> LLMResponse:
        return LLMResponse(content=self._response_text, tool_calls=self._tool_calls)
