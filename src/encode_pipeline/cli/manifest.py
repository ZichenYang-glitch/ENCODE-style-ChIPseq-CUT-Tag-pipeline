"""CLI entry point: encode-manifest.

Package-native manifest generator with structured logging.
Writes to stdout when no --output is given.
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
from encode_pipeline.config.validate import validate_config
from encode_pipeline.manifest import build_manifest_rows, write_manifest_tsv


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
    parser.add_argument(
        "--strict", action="store_true", help="Fail if any row is missing"
    )
    parser.add_argument("--outdir", help="Override config outdir")
    args = parser.parse_args(argv)

    if not args.config and not args.config_json:
        parser.error("--config or --config-json is required")

    cfg = _load_config(args.config, args.config_json)
    if args.outdir:
        cfg["outdir"] = args.outdir

    run_id, version, config_hash = make_context(cfg)
    emit(run_id, version, config_hash)

    # Always build rows via the package API.
    rows, missing_count, na_count = build_manifest_rows(cfg)

    if args.output:
        write_manifest_tsv(rows, args.output)
        print(f"Manifest written: {args.output}")
        print(
            f"  {len(rows)} rows: "
            f"{len(rows) - missing_count - na_count} present, "
            f"{missing_count} missing, "
            f"{na_count} not_applicable"
        )
        if args.strict and missing_count > 0:
            print(
                f"ERROR: --strict mode, {missing_count} missing output(s)",
                file=sys.stderr,
            )
            return 1
        return 0

    # Stdout mode: generate to a temp file, prepend structured header, stream.
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", delete=False, encoding="utf-8"
    ) as tmp:
        tmp_path = tmp.name

    try:
        write_manifest_tsv(rows, tmp_path)
        comments = header_comment_lines(run_id, version, config_hash)
        with open(tmp_path, "rb") as fh:
            body = fh.read()
        sys.stdout.buffer.write("\n".join(comments).encode() + b"\n" + body)
    finally:
        os.unlink(tmp_path)

    if args.strict and missing_count > 0:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
