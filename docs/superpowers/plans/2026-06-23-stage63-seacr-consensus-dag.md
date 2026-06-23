# Stage 63 Implementation Plan: SEACR Consensus DAG

**Date:** 2026-06-23
**Status:** Implemented
**Depends on:** Stage 62

## Files modified

| File | Change |
|------|--------|
| `workflow/rules/consensus.smk` | +4 SEACR rules + helpers |
| `workflow/rules/metadata.smk` | +SEACR consensus experiment list |
| `workflow/rules/targets.smk` | +SEACR consensus targets |
| `test/test_stage63_seacr_consensus_dryrun.py` | 11 tests |

## Non-goals

- No SEACR IDR (never)
- No changes to qc.smk, Snakefile, compute_consensus.py
- No Co-Authored-By
