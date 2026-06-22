#!/usr/bin/env bash
# chipseq.sh — DEPRECATED legacy single-sample entrypoint
#
# The canonical entry point is the Snakemake workflow. Run:
#
#   snakemake -s workflow/Snakefile --configfile config/config.yaml
#
# See README.md and docs/ for current usage.
# The historical full script is archived at docs/archive/scripts/chipseq-legacy.sh.

echo "ERROR: scripts/chipseq.sh is deprecated." >&2
echo "The canonical entry point is: snakemake -s workflow/Snakefile --configfile config/config.yaml" >&2
exit 1
