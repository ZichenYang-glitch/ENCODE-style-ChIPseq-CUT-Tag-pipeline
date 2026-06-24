#!/usr/bin/env python3
"""Pytest shim that runs legacy test_stage*.py scripts under pytest.

These scripts use hand-rolled pass/fail counters and expose one of three
entry-point conventions:

* a ``main()`` function
* a ``__main__`` guard block that calls test functions
* neither (pure pytest module, left to normal discovery)

The shim imports each script and invokes the appropriate entry point. Files
with main() or a guard are listed in conftest.py ``collect_ignore`` so pytest
only runs them through this shim.
"""

import importlib
import os
import runpy
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parent.parent
TEST_DIR = REPO_ROOT / "test"


def _discover_stage_modules():
    """Return Paths for all test_stage*.py modules except this shim."""
    return sorted(
        p
        for p in TEST_DIR.glob("test_stage*.py")
        if p.name != "test_stage_shim.py"
    )


def _has_main(src):
    return "def main(" in src


def _has_main_guard(src):
    return 'if __name__ == "__main__":' in src


def _run_main(module_name):
    """Import the module and call its main() function."""
    module = importlib.import_module(module_name)
    try:
        result = module.main()
    except SystemExit as exc:
        code = exc.code if exc.code is not None else 0
        if code != 0:
            pytest.fail(f"{module_name}.main() exited with code {code}")
        return
    if result is not None and result != 0:
        pytest.fail(f"{module_name}.main() returned {result}")


def _run_guard(path):
    """Execute the module as __main__ so its guard block runs."""
    try:
        runpy.run_path(str(path), run_name="__main__")
    except SystemExit as exc:
        code = exc.code if exc.code is not None else 0
        if code != 0:
            pytest.fail(f"{path.name} __main__ guard exited with code {code}")


def _module_name(path):
    """Return a dotted module name for an import."""
    # The test directory is on sys.path when pytest runs, so top-level imports
    # work for modules like test_stage25_manifest_stress.
    return path.stem


@pytest.mark.parametrize(
    "stage_path",
    _discover_stage_modules(),
    ids=lambda p: p.name,
)
def test_stage_shim(stage_path):
    """Run a legacy test_stage*.py script under pytest."""
    src = stage_path.read_text(encoding="utf-8")
    module_name = _module_name(stage_path)

    if _has_main(src):
        _run_main(module_name)
    elif _has_main_guard(src):
        _run_guard(stage_path)
    else:
        pytest.skip(f"{stage_path.name} has no shim entry point; collected normally")
