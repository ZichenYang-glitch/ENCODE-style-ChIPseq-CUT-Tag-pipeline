"""Explicit server composition for verified Hi-TrAC results and execution."""

from __future__ import annotations

import os
import hashlib
from pathlib import Path

from .adapter import HiTracPreprocessAdapter
from .execution import admit_runtime
from .results import HiTracPreprocessResultsAdapter


# This policy resource is outside the installed Bulk Python-file manifest.
# Bind it here because its pins participate in the shared runner's admission
# decision; changing those pins requires new controlled implementation bytes.
_RUNTIME_LOCK_SHA256 = (
    "392e4a4ff93f2b14b380f16572d74a77dbca699112800bcafe50992d692090da"
)


def load_default_hitrac_adapter(environ=None):
    values = os.environ if environ is None else environ
    path = values.get("HELIXWEAVE_HITRAC_RUNTIME_BINDING")
    expected = values.get("HELIXWEAVE_HITRAC_RUNTIME_SHA256")
    if not path or not expected:
        return HiTracPreprocessResultsAdapter()
    try:
        return HiTracPreprocessResultsAdapter(
            runtime=admit_runtime(Path(path), expected)
        )
    except (OSError, TypeError, ValueError):
        return HiTracPreprocessResultsAdapter()


def local_execution_executable(adapter):
    """Trust only this server-composed adapter's live, admitted Python binding."""
    if not isinstance(adapter, HiTracPreprocessAdapter) or adapter._runtime is None:
        return None
    try:
        lock = (
            Path(__file__).resolve().parents[4]
            / "config/hitrac_preprocess/tools.lock.json"
        )
        if hashlib.sha256(lock.read_bytes()).hexdigest() != _RUNTIME_LOCK_SHA256:
            return None
        adapter._runtime.verify()
    except (OSError, TypeError, ValueError):
        return None
    return str(adapter._runtime.python)
