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

        SNAKEMAKE=/path/to/snakemake python3 test/test_stage8_smoke_profiles.py
    """
    env_value = os.environ.get(env_var, "")
    if env_value:
        return env_value
    return shutil.which(name) or name
