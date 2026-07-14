"""Tests for the OpenAI-backed LLMClient implementation."""

from __future__ import annotations

import asyncio
from dataclasses import dataclass
from typing import Any

import pytest

from encode_pipeline.services.llm_client import LLMClient, LLMMessage, LLMResponse
from encode_pipeline.services.llm_factory import LLMProviderError
from encode_pipeline.services.llm_providers.openai_client import OpenAILLMClient


@dataclass
class FakeCompletionChoice:
    message: FakeCompletionMessage


@dataclass
class FakeCompletionMessage:
    content: str | None


@dataclass
class FakeCompletionResponse:
    choices: list[FakeCompletionChoice]


@dataclass
class FakeCompletions:
    response_content: str | None = "fake response"

    async def create(
        self,
        *,
        model: str,
        messages: list[dict[str, Any]],
        tools: list[dict[str, Any]] | None = None,
        **kwargs: Any,
    ) -> FakeCompletionResponse:
        self.last_model = model
        self.last_messages = messages
        self.last_tools = tools
        return FakeCompletionResponse(
            choices=[
                FakeCompletionChoice(
                    message=FakeCompletionMessage(content=self.response_content)
                )
            ]
        )


@dataclass
class FakeChat:
    def __init__(self) -> None:
        self.completions = FakeCompletions()


@dataclass
class FakeOpenAIClient:
    def __init__(self) -> None:
        self.chat = FakeChat()


class RaisingFakeCompletions:
    def __init__(self, exception: Exception) -> None:
        self._exception = exception

    async def create(self, **kwargs: Any) -> FakeCompletionResponse:
        raise self._exception


class RaisingFakeChat:
    def __init__(self, exception: Exception) -> None:
        self.completions = RaisingFakeCompletions(exception)


class RaisingFakeOpenAIClient:
    def __init__(self, exception: Exception) -> None:
        self.chat = RaisingFakeChat(exception)


def _run(client: LLMClient, messages: list[LLMMessage], tools: list) -> LLMResponse:
    return asyncio.run(client.complete(messages, tools))


def test_openai_client_implements_llm_client_protocol():
    fake = FakeOpenAIClient()
    client = OpenAILLMClient(api_key="fake-key", client=fake)
    assert hasattr(client, "complete")
    assert asyncio.iscoroutinefunction(client.complete)
    response = _run(client, [], [])
    assert isinstance(response, LLMResponse)


def test_complete_returns_llm_response_with_content():
    fake = FakeOpenAIClient()
    client = OpenAILLMClient(api_key="fake-key", model="gpt-4o", client=fake)
    response = _run(
        client,
        [LLMMessage(role="user", content="hello")],
        [],
    )
    assert isinstance(response, LLMResponse)
    assert response.content == "fake response"
    assert response.tool_calls == ()


def test_complete_maps_messages_to_openai_format():
    fake = FakeOpenAIClient()
    client = OpenAILLMClient(api_key="fake-key", model="gpt-4o", client=fake)
    _run(
        client,
        [
            LLMMessage(role="system", content="you are a tester"),
            LLMMessage(role="user", content="hello"),
        ],
        [],
    )
    assert fake.chat.completions.last_model == "gpt-4o"
    assert fake.chat.completions.last_messages == [
        {"role": "system", "content": "you are a tester"},
        {"role": "user", "content": "hello"},
    ]


def test_default_model_is_gpt_4o_mini():
    fake = FakeOpenAIClient()
    client = OpenAILLMClient(api_key="fake-key", client=fake)
    _run(client, [LLMMessage(role="user", content="hi")], [])
    assert fake.chat.completions.last_model == "gpt-4o-mini"


def test_repr_does_not_expose_api_key():
    client = OpenAILLMClient(api_key="secret-key-123", client=FakeOpenAIClient())
    representation = repr(client)
    assert "secret-key-123" not in representation
    assert "api_key" not in representation


def test_str_does_not_expose_api_key():
    client = OpenAILLMClient(api_key="secret-key-123", client=FakeOpenAIClient())
    text = str(client)
    assert "secret-key-123" not in text
    assert "api_key" not in text


def test_sdk_exception_is_wrapped_and_key_is_redacted():
    api_key = "secret-key-123"
    exception = RuntimeError(f"request failed with api key {api_key}")
    fake = RaisingFakeOpenAIClient(exception)
    client = OpenAILLMClient(api_key=api_key, client=fake)

    with pytest.raises(LLMProviderError) as exc_info:
        _run(client, [LLMMessage(role="user", content="hi")], [])

    assert api_key not in str(exc_info.value)
    assert api_key not in repr(exc_info.value)
    assert "request failed" in str(exc_info.value)
    assert exc_info.value.__cause__ is exception


def test_provider_echo_response_does_not_leak_key_in_client_repr_or_str():
    api_key = "secret-key-123"
    fake = FakeOpenAIClient()
    fake.chat.completions.response_content = api_key
    client = OpenAILLMClient(api_key=api_key, client=fake)

    response = _run(client, [LLMMessage(role="user", content="hello")], [])

    assert api_key not in response.content
    assert api_key not in repr(client)
    assert api_key not in str(client)
