# Stage 60 Implementation Plan: Legacy Script Cleanup

**Date:** 2026-06-22
**Status:** Implemented

## Files changed

- `scripts/chipseq.sh`: replace legacy implementation with a deprecation shim
  that exits non-zero and points users to `workflow/Snakefile`.
- `docs/archive/scripts/chipseq-legacy.sh`: preserve the historical full script
  for reference only.
- `README.md`: mark the legacy script as deprecated and document the Snakemake
  command as the canonical entrypoint.
- `KNOWN_ISSUES.md`: close the legacy script hardening issue and clarify that
  the archive is not a supported entrypoint.
- `test/test_stage60_legacy_script.py`: enforce shim behavior and documentation
  consistency.

## Verification

- `bash -n scripts/chipseq.sh`
- `python3 test/test_stage60_legacy_script.py`
- `python3 test/test_stage28_release_readiness.py`
- `python3 test/test_no_hardcoded_paths.py`
- `git diff --check`

## Non-goals

- No `.smk`, `Snakefile`, config, manifest, or scientific behavior changes
- No Artifact runtime adoption
- No Co-Authored-By
