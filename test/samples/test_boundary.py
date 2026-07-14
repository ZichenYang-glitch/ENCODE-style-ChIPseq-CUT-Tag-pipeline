"""Boundary and import-surface tests for the sample loading extraction."""

import ast
import inspect
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest


SRC_ROOT = Path(__file__).resolve().parents[2] / "src"


def _samples_source_modules():
    """Yield all source modules under encode_pipeline.samples."""
    src = SRC_ROOT / "encode_pipeline" / "samples"
    for path in src.glob("*.py"):
        if path.name.startswith("_"):
            continue
        yield path


def _forbidden_imports_in(path: Path):
    """Return any import statements that pull from config.validator/validate."""
    tree = ast.parse(path.read_text())
    bad = []
    for node in ast.walk(tree):
        if isinstance(node, ast.ImportFrom):
            module = node.module or ""
            if module.startswith("encode_pipeline.config.validator"):
                bad.append((node.lineno, module))
            if module.startswith("encode_pipeline.config.validate"):
                bad.append((node.lineno, module))
        elif isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name.startswith("encode_pipeline.config.validator"):
                    bad.append((node.lineno, alias.name))
                if alias.name.startswith("encode_pipeline.config.validate"):
                    bad.append((node.lineno, alias.name))
    return bad


def _run_fresh_import(code: str):
    """Run *code* in a subprocess with clean PYTHONPATH pointing at src."""
    import os

    env = dict(os.environ)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["PYTHONPATH"] = str(SRC_ROOT)
    result = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(code)],
        env=env,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"Fresh import failed:\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    return result


# ---------------------------------------------------------------------------
# Public surface
# ---------------------------------------------------------------------------

def test_validation_error_is_shared_singleton():
    from encode_pipeline.config.errors import ValidationError as errors_cls
    from encode_pipeline.config.validate import ValidationError as validate_cls
    from encode_pipeline.config.validator import ValidationError as validator_cls
    from encode_pipeline.samples.load import ValidationError as samples_cls

    assert errors_cls is validate_cls is validator_cls is samples_cls


def test_load_and_validate_samples_backward_compatible_reexport():
    from encode_pipeline.config.validator import (
        load_and_validate_samples as validator_fn,
    )
    from encode_pipeline.samples.load import (
        load_and_validate_samples as samples_fn,
    )
    from encode_pipeline.samples.loader import (
        load_and_validate_samples as loader_fn,
    )

    assert validator_fn is samples_fn is loader_fn


def test_validate_replicate_groups_backward_compatible_reexport():
    from encode_pipeline.config.validator import (
        validate_replicate_groups as validator_fn,
    )
    from encode_pipeline.samples.replicates import (
        validate_replicate_groups as replicates_fn,
    )

    assert validator_fn is replicates_fn


def test_validate_sample_row_public_helper():
    from encode_pipeline.samples.validate import validate_sample_row

    row = {
        "sample": "s1",
        "fastq_1": "/dev/null",
        "layout": "SE",
        "assay": "chipseq",
        "target": "H3K4me3",
        "peak_mode": "narrow",
        "genome": "hg38",
        "bowtie2_index": "/dev/null",
    }
    record = validate_sample_row(row)
    assert record.id == "s1"
    assert record.layout == "SE"


# ---------------------------------------------------------------------------
# Fresh subprocess imports (detects circular-import issues pytest can hide)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "code",
    [
        'from encode_pipeline.samples.load import load_and_validate_samples, ValidationError; print("ok")',
        'from encode_pipeline.samples import load_and_validate_samples; print("ok")',
        'from encode_pipeline.samples.validate import validate_sample_row; print("ok")',
        'from encode_pipeline.config.validator import load_and_validate_samples, validate_replicate_groups; print("ok")',
        'from encode_pipeline.config.validate import ValidationError; print("ok")',
    ],
    ids=[
        "samples.load",
        "samples.package",
        "samples.validate",
        "config.validator",
        "config.validate",
    ],
)
def test_fresh_import_in_subprocess(code):
    _run_fresh_import(code)


# ---------------------------------------------------------------------------
# Forbidden import check
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("path", list(_samples_source_modules()), ids=lambda p: p.name)
def test_samples_modules_do_not_import_config_validator_or_validate(path):
    bad = _forbidden_imports_in(path)
    assert not bad, f"Forbidden imports in {path.name}: {bad}"


# ---------------------------------------------------------------------------
# Internal default exception class
# ---------------------------------------------------------------------------

def test_public_functions_default_to_validation_error():
    """Public functions in samples.* must default to ValidationError."""
    from encode_pipeline.config.errors import ValidationError
    from encode_pipeline.samples import loader, replicates, rows, strict

    for module, name in [
        (rows, "validate_and_build_sample"),
        (strict, "validate_strict_inputs"),
        (replicates, "validate_replicate_groups"),
        (loader, "load_and_validate_samples"),
    ]:
        fn = getattr(module, name)
        params = inspect.signature(fn).parameters
        error_param = params.get("error_cls")
        if error_param is not None:
            assert error_param.default is ValidationError, (
                f"{module.__name__}.{name} error_cls default is not ValidationError"
            )
