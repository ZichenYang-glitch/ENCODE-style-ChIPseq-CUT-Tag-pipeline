#!/usr/bin/env python3
"""Generate a project-level result manifest recording core output existence.

Thin wrapper around encode_pipeline.manifest for backward compatibility with
existing Snakemake rules and standalone usage.
"""

import sys

from encode_pipeline.manifest.make import main


if __name__ == "__main__":
    sys.exit(main())
