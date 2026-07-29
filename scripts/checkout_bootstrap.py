#!/usr/bin/env python3
"""Clean, checkout-only launcher for HelixWeave development verification."""

from __future__ import annotations

import argparse
import importlib.util
import os
from pathlib import Path
import runpy
import sys
from types import ModuleType


_COMMANDS = (
    "verify-checkout",
    "installed-artifact",
    "pytest",
    "openapi",
    "validate",
    "local-platform",
)


def _abort(reason_code: str, guidance: str) -> None:
    print(
        f"source provenance check failed [{reason_code}]: {guidance}",
        file=sys.stderr,
    )
    raise SystemExit(2)


def _load_provenance(script_path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location(
        "_helixweave_checkout_source_provenance",
        script_path,
    )
    if spec is None or spec.loader is None:
        _abort(
            "bootstrap_source_invalid",
            "restore the checkout provenance source before retrying",
        )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
    except (ImportError, OSError, RuntimeError, SyntaxError):
        _abort(
            "bootstrap_source_invalid",
            "restore the checkout provenance source before retrying",
        )
    return module


def _parse_arguments() -> tuple[argparse.Namespace, list[str]]:
    parser = argparse.ArgumentParser(
        add_help=False,
        description="Run a fixed HelixWeave checkout verification command.",
    )
    parser.add_argument("--repository-root", type=Path, required=True)
    parser.add_argument("command", choices=_COMMANDS)
    parser.add_argument("arguments", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    return args, args.arguments


def _reject_pytest_plugin_injection(arguments: list[str]) -> None:
    if os.environ.get("PYTEST_ADDOPTS") or os.environ.get("PYTEST_PLUGINS"):
        _abort(
            "pytest_plugin_unsafe",
            "remove pytest environment plugin injection before retrying",
        )
    for argument in arguments:
        if (
            argument == "-p"
            or argument.startswith("-p")
            or argument == "--plugins"
            or argument.startswith("--plugins=")
            or argument == "-c"
            or argument.startswith("-c")
            or argument == "--config-file"
            or argument.startswith("--config-file=")
            or argument == "-o"
            or argument.startswith("-o")
            or argument == "--override-ini"
            or argument.startswith("--override-ini=")
            or argument.startswith("@")
        ):
            _abort(
                "pytest_plugin_unsafe",
                "use the repository pytest config and bootstrap-managed plugins",
            )


def _run_script(path: Path, arguments: list[str]) -> int:
    sys.argv = [str(path), *arguments]
    try:
        runpy.run_path(str(path), run_name="__main__")
    except SystemExit as error:
        if error.code is None:
            return 0
        if isinstance(error.code, int):
            return error.code
        return 1
    return 0


def _run_pytest(repository_root: Path, arguments: list[str]) -> int:
    os.environ["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
    try:
        import pytest
        import pytest_cov.plugin as pytest_cov_plugin
    except ImportError:
        _abort(
            "pytest_plugin_missing",
            "install the locked pytest and pytest-cov development dependencies",
        )
    return int(
        pytest.main(
            [
                "-c",
                str(repository_root / "pyproject.toml"),
                "-p",
                "no:cacheprovider",
                *arguments,
            ],
            plugins=[pytest_cov_plugin],
        )
    )


def main() -> int:
    sys.dont_write_bytecode = True
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    if not sys.flags.isolated or not sys.flags.no_site:
        _abort(
            "bootstrap_startup_unsafe",
            "invoke this checkout launcher with both Python -I and -S",
        )

    args, payload_arguments = _parse_arguments()
    try:
        repository_root = args.repository_root.resolve(strict=True)
        bootstrap_root = Path(__file__).resolve(strict=True).parents[1]
    except (OSError, RuntimeError):
        _abort(
            "repository_root_invalid",
            "provide the expected existing repository root explicitly",
        )
    if repository_root != bootstrap_root:
        _abort(
            "repository_root_invalid",
            "run the bootstrap that belongs to the expected repository root",
        )

    provenance = _load_provenance(bootstrap_root / "scripts" / "source_provenance.py")
    if args.command == "pytest":
        _reject_pytest_plugin_injection(payload_arguments)

    try:
        if args.command == "installed-artifact":
            provenance.bootstrap_installed_artifact()
        else:
            provenance.bootstrap_checkout(repository_root)
    except provenance.SourceProvenanceError as error:
        print(error, file=sys.stderr)
        return 2

    if args.command in {"verify-checkout", "installed-artifact"}:
        return 0
    if args.command == "pytest":
        return _run_pytest(repository_root, payload_arguments)

    targets = {
        "openapi": bootstrap_root / "scripts" / "export_openapi.py",
        "validate": bootstrap_root / "scripts" / "validate_samples.py",
        "local-platform": bootstrap_root / "scripts" / "run_local_platform.py",
    }
    return _run_script(targets[args.command], payload_arguments)


if __name__ == "__main__":
    raise SystemExit(main())
