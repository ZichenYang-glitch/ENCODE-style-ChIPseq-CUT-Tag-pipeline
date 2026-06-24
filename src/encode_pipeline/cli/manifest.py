"""CLI entry point: encode-manifest.

Thin wrapper around scripts.make_manifest that adds structured logging and
writes to stdout when no --output is given.
"""

import argparse
import json
import os
import sys
import tempfile

import yaml

from encode_pipeline.cli._logging import (
    emit,
    header_comment_lines,
    make_context,
)

# Make scripts.make_manifest importable whether the package is installed or not.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))
import scripts.make_manifest as _make_manifest
from encode_pipeline.config.validate import validate_config


def _load_config(config_path, config_json):
    """Load and validate config from YAML or JSON."""
    if config_json:
        return json.loads(config_json)
    if config_path:
        with open(config_path) as fh:
            raw = yaml.safe_load(fh)
        return validate_config(raw)
    raise ValueError("--config or --config-json is required")


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Generate the project-level result manifest (TSV)."
    )
    parser.add_argument("--config", help="Path to config.yaml")
    parser.add_argument("--config-json", help="JSON dump of validated config")
    parser.add_argument("--output", "-o", help="Output TSV path (default: stdout)")
    parser.add_argument("--strict", action="store_true", help="Fail if any row is missing")
    parser.add_argument("--outdir", help="Override config outdir")
    args = parser.parse_args(argv)

    if not args.config and not args.config_json:
        parser.error("--config or --config-json is required")

    cfg = _load_config(args.config, args.config_json)
    if args.outdir:
        cfg["outdir"] = args.outdir

    run_id, version, config_hash = make_context(cfg)
    emit(run_id, version, config_hash)

    if args.output:
        # File mode: delegate directly to make_manifest for byte-identical output.
        argv_for_make = ["make_manifest.py"]
        if args.config_json:
            argv_for_make.extend(["--config-json", args.config_json])
        elif args.config:
            argv_for_make.extend(["--config", args.config])
        argv_for_make.extend(["--output", args.output])
        if args.strict:
            argv_for_make.append("--strict")

        old_argv = sys.argv
        sys.argv = argv_for_make
        try:
            _make_manifest.main()
        finally:
            sys.argv = old_argv
        return

    # Stdout mode: generate to a temp file, prepend structured header, stream.
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", delete=False, encoding="utf-8"
    ) as tmp:
        tmp_path = tmp.name

    try:
        argv_for_make = [
            "make_manifest.py",
            "--config-json",
            json.dumps(cfg),
            "--output",
            tmp_path,
        ]
        if args.strict:
            argv_for_make.append("--strict")

        old_stdout = sys.stdout
        old_argv = sys.argv
        # Redirect make_manifest's summary prints to stderr so they don't
        # contaminate the TSV stream.
        sys.stdout = sys.stderr
        sys.argv = argv_for_make
        try:
            _make_manifest.main()
        finally:
            sys.stdout = old_stdout
            sys.argv = old_argv

        comments = header_comment_lines(run_id, version, config_hash)
        with open(tmp_path, "rb") as fh:
            body = fh.read()
        sys.stdout.buffer.write("\n".join(comments).encode() + b"\n" + body)
    finally:
        os.unlink(tmp_path)


if __name__ == "__main__":
    main()
