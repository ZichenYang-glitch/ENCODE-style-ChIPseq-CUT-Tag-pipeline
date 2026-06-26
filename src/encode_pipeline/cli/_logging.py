"""Structured logging helpers for CLI entry points."""

import hashlib
import json
import sys
import uuid

from encode_pipeline import __version__


def make_context(config=None):
    """Return a structured logging context tuple.

    Args:
        config: optional validated config dict to hash.

    Returns:
        (run_id, version, config_hash)
    """
    run_id = uuid.uuid4().hex[:12]
    if config is None:
        config_hash = "na"
    else:
        try:
            config_hash = hashlib.sha256(
                json.dumps(config, sort_keys=True).encode()
            ).hexdigest()[:16]
        except TypeError:
            config_hash = "unhashable"
    return run_id, __version__, config_hash


def emit(run_id, version, config_hash, file=sys.stderr):
    """Emit the standard CLI log line."""
    print(
        f"[{run_id}] encode-pipeline v{version} config_hash={config_hash}",
        file=file,
    )


def header_comment_lines(run_id, version, config_hash):
    """Return manifest header comment lines."""
    return [
        f"# run_id={run_id}",
        f"# pipeline_version={version}",
        f"# config_hash={config_hash}",
    ]
