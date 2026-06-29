"""Pin and guard ``encode_pipeline.config`` package-level import behavior.

These tests run in subprocesses so that import-state assertions are not
polluted by pytest's own import cache.
"""

import os
import subprocess
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC = REPO_ROOT / "src"


def _run_subprocess(code):
    env = os.environ.copy()
    env["PYTHONPATH"] = str(SRC)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    proc = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        env=env,
    )
    return proc


# ---------------------------------------------------------------------------
# Plain package import should not eagerly load validator/validate
# ---------------------------------------------------------------------------


def test_import_config_does_not_load_validator():
    code = (
        "import sys\n"
        "import encode_pipeline.config\n"
        "print('encode_pipeline.config.validator' in sys.modules)"
    )
    proc = _run_subprocess(code)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "False"


def test_import_config_does_not_load_validate():
    code = (
        "import sys\n"
        "import encode_pipeline.config\n"
        "print('encode_pipeline.config.validate' in sys.modules)"
    )
    proc = _run_subprocess(code)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "False"


def test_import_config_does_not_load_yaml_loader():
    code = (
        "import sys\n"
        "import encode_pipeline.config\n"
        "print('encode_pipeline.config.yaml_loader' in sys.modules)"
    )
    proc = _run_subprocess(code)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "False"


# ---------------------------------------------------------------------------
# Lazy access loads validate/validator only when needed
# ---------------------------------------------------------------------------


def test_config_attribute_access_validate_config():
    code = (
        "import encode_pipeline.config as cfg\n"
        "print(callable(cfg.validate_config))"
    )
    proc = _run_subprocess(code)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "True"


def test_lazy_access_loads_validate_and_validator_only_when_needed():
    code = (
        "import sys\n"
        "import encode_pipeline.config as cfg\n"
        "before_validate = 'encode_pipeline.config.validate' in sys.modules\n"
        "before_validator = 'encode_pipeline.config.validator' in sys.modules\n"
        "_ = cfg.validate_config\n"
        "after_validate = 'encode_pipeline.config.validate' in sys.modules\n"
        "after_validator = 'encode_pipeline.config.validator' in sys.modules\n"
        "print(before_validate, before_validator, after_validate, after_validator)"
    )
    proc = _run_subprocess(code)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "False False True True"


# ---------------------------------------------------------------------------
# Supported public imports still work
# ---------------------------------------------------------------------------


def test_from_config_import_validate_config():
    code = (
        "from encode_pipeline.config import validate_config\n"
        "print(callable(validate_config))"
    )
    proc = _run_subprocess(code)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "True"


def test_from_config_import_resource_validators():
    code = (
        "from encode_pipeline.config import (\n"
        "    validate_picard_reference_resources,\n"
        "    validate_tss_annotation_resources,\n"
        ")\n"
        "print(callable(validate_picard_reference_resources),\n"
        "      callable(validate_tss_annotation_resources))"
    )
    proc = _run_subprocess(code)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "True True"


def test_config_validate_import_still_works():
    code = (
        "from encode_pipeline.config.validate import validate_config\n"
        "print(callable(validate_config))"
    )
    proc = _run_subprocess(code)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "True"


def test_config_validator_import_still_works():
    code = (
        "from encode_pipeline.config.validator import validate_config\n"
        "print(callable(validate_config))"
    )
    proc = _run_subprocess(code)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "True"


# ---------------------------------------------------------------------------
# Unsupported attribute behavior
# ---------------------------------------------------------------------------


def test_unsupported_attribute_raises_attribute_error():
    code = (
        "import encode_pipeline.config as cfg\n"
        "try:\n"
        "    cfg.not_a_real_name\n"
        "except AttributeError as exc:\n"
        "    print('has no attribute' in str(exc))\n"
        "else:\n"
        "    print('False')"
    )
    proc = _run_subprocess(code)
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "True"
