"""Tests for the ProcessRunner subprocess boundary."""

import ast
import os
import sys
from pathlib import Path

import pytest

from encode_pipeline.platform.adapters import CommandSpec
from encode_pipeline.services.process_runner import ProcessResult, ProcessRunner


def _make_runner(**kwargs):
    """Return a ProcessRunner that allows the current Python executable."""
    exe = sys.executable
    name = Path(exe).name
    defaults = {"allowed_executables": (exe, name)}
    defaults.update(kwargs)
    return ProcessRunner(**defaults)


# ---------------------------------------------------------------------------
# Constructor validation
# ---------------------------------------------------------------------------


def test_process_runner_default_construction():
    runner = ProcessRunner()
    assert runner._allowed_executables == ("snakemake",)
    assert runner._timeout_seconds == 300.0
    assert runner._max_output_bytes == 10_000_000


def test_process_runner_rejects_empty_allowed_executables():
    with pytest.raises(ValueError, match="non-empty"):
        ProcessRunner(allowed_executables=())


def test_process_runner_rejects_empty_string_in_allowed_executables():
    with pytest.raises(ValueError, match="non-empty strings"):
        ProcessRunner(allowed_executables=("snakemake", ""))


def test_process_runner_rejects_non_tuple_allowed_executables():
    with pytest.raises(ValueError, match="must be a tuple"):
        ProcessRunner(allowed_executables=["snakemake"])


def test_process_runner_rejects_zero_timeout():
    with pytest.raises(ValueError, match="positive"):
        ProcessRunner(timeout_seconds=0)


def test_process_runner_rejects_negative_timeout():
    with pytest.raises(ValueError, match="positive"):
        ProcessRunner(timeout_seconds=-1.0)


def test_process_runner_rejects_bool_timeout():
    with pytest.raises(ValueError, match="must be a number"):
        ProcessRunner(timeout_seconds=True)


def test_process_runner_rejects_zero_max_output_bytes():
    with pytest.raises(ValueError, match="positive"):
        ProcessRunner(max_output_bytes=0)


def test_process_runner_rejects_negative_max_output_bytes():
    with pytest.raises(ValueError, match="positive"):
        ProcessRunner(max_output_bytes=-1)


def test_process_runner_rejects_bool_max_output_bytes():
    with pytest.raises(ValueError, match="must be an int"):
        ProcessRunner(max_output_bytes=True)
