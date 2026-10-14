#!/usr/bin/env python3
"""Spawned subprocesses must resolve encode_pipeline to this checkout.

The session fixture in test/conftest.py prepends <repo>/src to PYTHONPATH so
that tools spawned by tests (notably snakemake running the Snakefile under a
different interpreter) import this checkout's package rather than whatever an
editable install in a shared environment points at.
"""

import os
from pathlib import Path
import subprocess
import sys


REPO = Path(__file__).resolve().parents[1]


def test_subprocess_env_prepends_checkout_src():
    pythonpath = os.environ.get("PYTHONPATH", "")
    assert pythonpath.split(os.pathsep)[0] == str(REPO / "src")


def test_spawned_python_imports_checkout_package():
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import encode_pipeline; print(encode_pipeline.__file__)",
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == str(
        REPO / "src" / "encode_pipeline" / "__init__.py"
    )
