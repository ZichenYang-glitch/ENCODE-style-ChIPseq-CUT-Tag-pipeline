"""CLI entry point for config/sample validation with structured logging."""

import os
import sys

import yaml

from encode_pipeline.cli._logging import emit, make_context
from encode_pipeline.config.validate import validate_config
from encode_pipeline.cli._validator import main as _validator_main


def main(argv=None):
    # Parse just enough to compute the config hash without duplicating the
    # full validator argument handling.
    if argv is None:
        argv = sys.argv[1:]
    config_path = None
    for i, arg in enumerate(argv):
        if arg == "--config" and i + 1 < len(argv):
            config_path = argv[i + 1]
            break
        if arg.startswith("--config="):
            config_path = arg.split("=", 1)[1]
            break

    cfg = None
    if config_path and os.path.isfile(config_path):
        try:
            with open(config_path) as fh:
                raw = yaml.safe_load(fh)
            cfg = validate_config(raw)
        except Exception:
            pass

    run_id, version, config_hash = make_context(cfg)
    emit(run_id, version, config_hash)

    return _validator_main(argv)


if __name__ == "__main__":
    sys.exit(main())
