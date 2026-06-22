# Stage 58 Implementation Plan: Mixed Experiment IDR Validation Relax

**Date:** 2026-06-19
**Status:** Implemented

## Files

| File | Change |
|------|--------|
| `scripts/validate_samples.py` | Stage 5 IDR: skip non-chipseq, require >=1 eligible |
| `test/test_stage58_mixed_idr_validation.py` | 10 mixed-sheet tests |

## Steps

1. Relax Stage 5 IDR block to skip non-chipseq experiments
2. Write 10 tests covering mixed sheets, negative cases, invariants
3. Verify regression tests pass

## Non-goals

- No Snakemake, target, config, manifest, or policy changes
- No Co-Authored-By
