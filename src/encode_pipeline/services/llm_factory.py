"""Factory for creating LLM clients based on environment configuration."""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

from encode_pipeline.services.llm_client import MockLLMClient

if TYPE_CHECKING:
    from encode_pipeline.services.llm_client import LLMClient


DEFAULT_MOCK_RESPONSE = (
    "This is a deterministic mock explanation from the workflow agent."
)


class LLMProviderError(Exception):
    """Base error for LLM provider failures."""


class LLMProviderConfigurationError(LLMProviderError):
    """Raised when LLM provider configuration is invalid or missing."""


def create_llm_client() -> "LLMClient":
    """Create an LLM client from environment variables.

    Reads ``ENCODE_AGENT_LLM_PROVIDER``, ``ENCODE_AGENT_LLM_API_KEY``, and
    ``ENCODE_AGENT_LLM_MODEL``. The provider value is normalized by stripping
    whitespace and lowercasing. An unset provider, an empty string, or
    ``"mock"`` returns a deterministic ``MockLLMClient``. ``"openai"`` lazily
    imports ``OpenAILLMClient`` and requires an API key. Any other provider or
    missing configuration raises ``LLMProviderConfigurationError``.

    Returns:
        An instance implementing the ``LLMClient`` protocol.

    Raises:
        LLMProviderConfigurationError: If the provider is unknown or required
            configuration (e.g. the real provider API key) is missing.
    """
    provider = os.environ.get("ENCODE_AGENT_LLM_PROVIDER", "").strip().lower()
    api_key = os.environ.get("ENCODE_AGENT_LLM_API_KEY", "").strip()
    model = os.environ.get("ENCODE_AGENT_LLM_MODEL", "").strip()

    if not provider or provider == "mock":
        return MockLLMClient(response_text=DEFAULT_MOCK_RESPONSE)

    if provider == "openai":
        if not api_key:
            raise LLMProviderConfigurationError(
                "Real provider requires ENCODE_AGENT_LLM_API_KEY"
            )
        try:
            from encode_pipeline.services.llm_providers.openai_client import (
                OpenAILLMClient,
            )

            kwargs = {}
            if model:
                kwargs["model"] = model
            return OpenAILLMClient(api_key=api_key, **kwargs)
        except ImportError as exc:
            raise LLMProviderConfigurationError(
                "Real provider is not available; install the optional dependencies"
            ) from exc

    raise LLMProviderConfigurationError(f"Unknown LLM provider: {provider!r}")
