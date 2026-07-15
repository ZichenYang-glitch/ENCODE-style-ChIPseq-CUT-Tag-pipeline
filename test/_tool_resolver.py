"""Shared test helpers for resolving external command-line tools."""

import os
import shutil
import tempfile


_TEST_CACHE = os.path.join(tempfile.gettempdir(), "chipseq-test-cache")
os.makedirs(_TEST_CACHE, exist_ok=True)
os.environ.setdefault("XDG_CACHE_HOME", _TEST_CACHE)


def resolve_tool(name, env_var):
    """Return a tool path from an environment variable or PATH.

    Tests should not bake in workstation-specific Conda paths.  Set the
    environment variable when a tool is outside PATH, for example:

        SNAKEMAKE=/path/to/snakemake python3 -m pytest test/workflow/test_smoke_profiles.py
    """
    env_value = os.environ.get(env_var, "")
    if env_value:
        return env_value
    path_value = shutil.which(name)
    if path_value:
        return path_value
    if name == "snakemake":
        conda_snakemake = os.path.join(
            os.path.expanduser("~"),
            "miniconda3",
            "envs",
            "chipseq",
            "bin",
            "snakemake",
        )
        if os.path.exists(conda_snakemake) and os.access(conda_snakemake, os.X_OK):
            return conda_snakemake
    return name


def require_external_tools(missing, context="real execution"):
    """Skip locally or fail a required CI gate when tools are unavailable."""
    if not missing:
        return

    import pytest

    message = f"{context} requires executable tools: {', '.join(missing)}"
    if os.environ.get("HELIXWEAVE_REQUIRE_REAL_EXECUTION") == "1":
        pytest.fail(message, pytrace=False)
    pytest.skip(message)
