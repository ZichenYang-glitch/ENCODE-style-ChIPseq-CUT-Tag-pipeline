#!/usr/bin/env python3
"""Generate a project-level result manifest recording core output existence.

Thin wrapper around encode_pipeline.manifest for backward compatibility with
existing Snakemake rules and standalone usage.
"""

import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
# Pin imports to this checkout: the result_manifest rule runs this script
# under the bare chipseq-python conda env (no editable install), and in a
# shared environment an editable install could otherwise resolve
# encode_pipeline to a different checkout. Insert-only (no require_checkout)
# because that conda env has no HelixWeave distribution to verify against.
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from encode_pipeline.manifest.make import main  # noqa: E402


if __name__ == "__main__":
    sys.exit(main())
