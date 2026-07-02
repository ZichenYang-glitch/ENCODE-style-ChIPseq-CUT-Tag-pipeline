"""Static import-boundary tests for the PR92 agent model/LLM layer.

These tests verify that the agent API models and LLM client module do not
depend on execution, adapter, validator, or provider-specific modules.
"""

from __future__ import annotations

import ast
import importlib
import inspect
from pathlib import Path

import pytest

# Modules introduced or hardened in PR92.
_AGENT_MODULES = [
    "encode_pipeline.api.models",
    "encode_pipeline.services.llm_client",
]

# Top-level module/package names that PR92 agent/LLM code must not import.
_FORBIDDEN_TOP_LEVEL = {
    # Execution and shell access
    "snakemake",
    "subprocess",
    # Adapter / validation internals (out of scope for PR92)
    "encode_pipeline.adapters",
    "encode_pipeline.config.validator",
    "encode_pipeline.samples",
    # Provider SDKs (deferred to PR97)
    "openai",
    "anthropic",
    "google",
    "boto3",
    "azure",
    "langchain",
    "langgraph",
    "mastra",
    "transformers",
    "torch",
}


def _module_source_path(module_name: str) -> Path:
    module = importlib.import_module(module_name)
    file_path = Path(inspect.getfile(module))
    assert file_path.exists(), f"Could not resolve source for {module_name}"
    return file_path


def _collect_imported_top_levels(source: str) -> set[str]:
    tree = ast.parse(source)
    imported: set[str] = set()

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                imported.add(alias.name.split(".")[0])
        elif isinstance(node, ast.ImportFrom):
            if node.module is not None:
                imported.add(node.module.split(".")[0])

    return imported


def _collect_imported_dotted_names(source: str) -> set[str]:
    tree = ast.parse(source)
    imported: set[str] = set()

    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                imported.add(alias.name)
        elif isinstance(node, ast.ImportFrom):
            if node.module is not None:
                prefix = node.module
                if node.level:
                    # Relative imports are not expected in PR92 public modules.
                    prefix = f"{'.' * node.level}{node.module}"
                imported.add(prefix)

    return imported


@pytest.mark.parametrize("module_name", _AGENT_MODULES)
def test_agent_module_does_not_import_forbidden_top_level(module_name: str):
    source = _module_source_path(module_name).read_text(encoding="utf-8")
    imported = _collect_imported_top_levels(source)
    violations = sorted(imported & _FORBIDDEN_TOP_LEVEL)
    assert not violations, f"{module_name} imports forbidden top-level modules: {violations}"


@pytest.mark.parametrize("module_name", _AGENT_MODULES)
def test_agent_module_does_not_import_forbidden_packages(module_name: str):
    """Catch imports like encode_pipeline.adapters even if only top-level encode_pipeline was imported."""
    source = _module_source_path(module_name).read_text(encoding="utf-8")
    imported = _collect_imported_dotted_names(source)
    violations = sorted(
        name
        for name in imported
        if any(name == forbidden or name.startswith(forbidden + ".") for forbidden in _FORBIDDEN_TOP_LEVEL)
    )
    assert not violations, f"{module_name} imports forbidden packages: {violations}"


@pytest.mark.parametrize("module_name", _AGENT_MODULES)
def test_agent_module_does_not_import_adapters_or_validator_submodules(module_name: str):
    """Explicit guard for encode_pipeline.adapters and encode_pipeline.config.validator."""
    source = _module_source_path(module_name).read_text(encoding="utf-8")
    imported = _collect_imported_dotted_names(source)
    assert "encode_pipeline.adapters" not in imported, f"{module_name} imports encode_pipeline.adapters"
    assert "encode_pipeline.config.validator" not in imported, f"{module_name} imports encode_pipeline.config.validator"
    assert "encode_pipeline.samples" not in imported, f"{module_name} imports encode_pipeline.samples"
