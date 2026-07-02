"""Tests for the environment-gated LLM client factory."""

from __future__ import annotations

import builtins

import pytest

from encode_pipeline.services.llm_client import MockLLMClient
from encode_pipeline.services.llm_factory import (
    DEFAULT_MOCK_RESPONSE,
    LLMProviderConfigurationError,
    LLMProviderError,
    create_llm_client,
)


def _clear_llm_env(monkeypatch) -> None:
    monkeypatch.delenv("ENCODE_AGENT_LLM_PROVIDER", raising=False)
    monkeypatch.delenv("ENCODE_AGENT_LLM_API_KEY", raising=False)
    monkeypatch.delenv("ENCODE_AGENT_LLM_MODEL", raising=False)


def test_default_no_env_returns_mock_llm_client(monkeypatch):
    _clear_llm_env(monkeypatch)
    client = create_llm_client()
    assert isinstance(client, MockLLMClient)
    assert client._response_text == DEFAULT_MOCK_RESPONSE


def test_provider_mock_returns_mock_llm_client(monkeypatch):
    monkeypatch.setenv("ENCODE_AGENT_LLM_PROVIDER", "mock")
    client = create_llm_client()
    assert isinstance(client, MockLLMClient)
    assert client._response_text == DEFAULT_MOCK_RESPONSE


def test_provider_mock_case_and_whitespace_normalized(monkeypatch):
    monkeypatch.setenv("ENCODE_AGENT_LLM_PROVIDER", "  MOCK  ")
    client = create_llm_client()
    assert isinstance(client, MockLLMClient)


def test_provider_openai_missing_key_raises(monkeypatch):
    monkeypatch.setenv("ENCODE_AGENT_LLM_PROVIDER", "openai")
    monkeypatch.delenv("ENCODE_AGENT_LLM_API_KEY", raising=False)
    with pytest.raises(LLMProviderConfigurationError):
        create_llm_client()


def test_provider_openai_whitespace_only_key_raises(monkeypatch):
    monkeypatch.setenv("ENCODE_AGENT_LLM_PROVIDER", "openai")
    monkeypatch.setenv("ENCODE_AGENT_LLM_API_KEY", "   ")
    with pytest.raises(LLMProviderConfigurationError):
        create_llm_client()


def test_provider_openai_missing_key_error_does_not_leak_key(monkeypatch):
    monkeypatch.setenv("ENCODE_AGENT_LLM_PROVIDER", "openai")
    monkeypatch.delenv("ENCODE_AGENT_LLM_API_KEY", raising=False)
    with pytest.raises(LLMProviderConfigurationError) as exc_info:
        create_llm_client()
    assert "secret" not in str(exc_info.value).lower()


def test_provider_openai_missing_sdk_on_import_raises(monkeypatch):
    monkeypatch.setenv("ENCODE_AGENT_LLM_PROVIDER", "openai")
    monkeypatch.setenv("ENCODE_AGENT_LLM_API_KEY", "fake-key")

    real_import = builtins.__import__

    def _raising_import(name, *args, **kwargs):
        if name == "encode_pipeline.services.llm_providers.openai_client":
            raise ImportError("Simulated missing SDK dependency")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _raising_import)
    with pytest.raises(LLMProviderConfigurationError):
        create_llm_client()


def test_provider_openai_missing_sdk_on_client_construction_raises(monkeypatch):
    """Missing optional SDK during provider client construction must wrap."""
    monkeypatch.setenv("ENCODE_AGENT_LLM_PROVIDER", "openai")
    monkeypatch.setenv("ENCODE_AGENT_LLM_API_KEY", "fake-key")

    real_import = builtins.__import__

    def _delayed_raising_import(name, *args, **kwargs):
        result = real_import(name, *args, **kwargs)
        if name == "openai":
            raise ImportError("Simulated missing SDK dependency at construction")
        return result

    monkeypatch.setattr(builtins, "__import__", _delayed_raising_import)
    with pytest.raises(LLMProviderConfigurationError):
        create_llm_client()


def test_provider_openai_missing_sdk_error_does_not_leak_key(monkeypatch):
    monkeypatch.setenv("ENCODE_AGENT_LLM_PROVIDER", "openai")
    monkeypatch.setenv("ENCODE_AGENT_LLM_API_KEY", "secret-key-123")

    real_import = builtins.__import__

    def _raising_import(name, *args, **kwargs):
        if name == "encode_pipeline.services.llm_providers.openai_client":
            raise ImportError("Simulated missing SDK dependency")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _raising_import)
    with pytest.raises(LLMProviderConfigurationError) as exc_info:
        create_llm_client()
    assert "secret-key-123" not in str(exc_info.value)


def test_provider_openai_no_model_preserves_default(monkeypatch):
    """When no model env is set, the factory should not pass an empty model."""
    monkeypatch.setenv("ENCODE_AGENT_LLM_PROVIDER", "openai")
    monkeypatch.setenv("ENCODE_AGENT_LLM_API_KEY", "fake-key")
    monkeypatch.delenv("ENCODE_AGENT_LLM_MODEL", raising=False)

    real_import = builtins.__import__
    captured: dict = {}

    class _FakeProviderClient:
        def __init__(self, *, api_key: str, model: str = "gpt-4o-mini") -> None:
            captured["api_key"] = api_key
            captured["model"] = model

    def _capturing_import(name, *args, **kwargs):
        if name == "encode_pipeline.services.llm_providers.openai_client":
            module = real_import(name, *args, **kwargs)
            module.OpenAILLMClient = _FakeProviderClient
            return module
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _capturing_import)
    create_llm_client()
    assert captured["model"] == "gpt-4o-mini"
    assert captured["api_key"] == "fake-key"


def test_provider_openai_explicit_model_is_passed(monkeypatch):
    """When a model env is set, the factory should pass it through."""
    monkeypatch.setenv("ENCODE_AGENT_LLM_PROVIDER", "openai")
    monkeypatch.setenv("ENCODE_AGENT_LLM_API_KEY", "fake-key")
    monkeypatch.setenv("ENCODE_AGENT_LLM_MODEL", "custom-model")

    real_import = builtins.__import__
    captured: dict = {}

    class _FakeProviderClient:
        def __init__(self, *, api_key: str, model: str = "gpt-4o-mini") -> None:
            captured["api_key"] = api_key
            captured["model"] = model

    def _capturing_import(name, *args, **kwargs):
        if name == "encode_pipeline.services.llm_providers.openai_client":
            module = real_import(name, *args, **kwargs)
            module.OpenAILLMClient = _FakeProviderClient
            return module
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", _capturing_import)
    create_llm_client()
    assert captured["model"] == "custom-model"
    assert captured["api_key"] == "fake-key"


def test_provider_unknown_raises(monkeypatch):
    monkeypatch.setenv("ENCODE_AGENT_LLM_PROVIDER", "unknown-provider")
    with pytest.raises(LLMProviderConfigurationError):
        create_llm_client()


def test_provider_unknown_error_does_not_leak_key(monkeypatch):
    monkeypatch.setenv("ENCODE_AGENT_LLM_PROVIDER", "unknown-provider")
    monkeypatch.setenv("ENCODE_AGENT_LLM_API_KEY", "secret-key-123")
    with pytest.raises(LLMProviderConfigurationError) as exc_info:
        create_llm_client()
    assert "secret-key-123" not in str(exc_info.value)


def test_llm_provider_configuration_error_is_llm_provider_error():
    assert issubclass(LLMProviderConfigurationError, LLMProviderError)
