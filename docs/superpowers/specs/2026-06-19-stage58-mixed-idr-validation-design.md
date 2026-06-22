# Stage 58: Mixed Experiment IDR Validation Relax — Design Spec

**Date:** 2026-06-19
**Status:** Implemented
**Author:** YangZiChen-glitch (Kaslana)
**Scope:** P0 — allow mixed sample sheets with IDR validation

## 1. Problem

Legacy Stage 5 IDR validation (`validate_samples.py` lines 1636-1661) iterates
all `exp_treatments` and rejects every experiment with `assay != "chipseq"`.
A mixed sample sheet (e.g. ChIP-seq narrow + CUT&Tag + ATAC narrow) fails with
`"Stage 5 supports chipseq only"` even when the ChIP-seq experiments
themselves are valid.

Stage 55 ATAC IDR already uses the correct permissive pattern: skip non-ATAC
experiments, validate only ATAC narrow experiments, require at least one
eligible experiment.

## 2. Fix

Apply the same permissive pattern to Stage 5 IDR validation:

- Skip experiments with `assay != "chipseq"` silently
- Skip experiments with `peak_mode != "narrow"` (chipseq broad)
- Validate exactly 2 biological replicates for chipseq narrow experiments
- Require at least one eligible ChIP-seq narrow experiment
- Do not weaken validation for eligible experiments

## 3. Files

| File | Change |
|------|--------|
| `scripts/validate_samples.py` | Stage 5 IDR block: skip non-chipseq, require ≥1 eligible experiment (~30 lines) |
| `test/test_stage58_mixed_idr_validation.py` | Create — 10 mixed-sheet tests |
| `docs/superpowers/specs/2026-06-19-stage58-mixed-idr-validation-design.md` | Create |
| `docs/superpowers/plans/2026-06-19-stage58-mixed-idr-validation.md` | Create |

## 4. Tests (M1-M10)

| # | Scenario | Expected |
|---|----------|---------|
| M1 | stage5 + ChIP-seq narrow + CUT&Tag | pass |
| M2 | stage5 + ChIP-seq narrow + ATAC narrow | pass |
| M3 | stage5 + ChIP-seq narrow + MNase | pass |
| M4 | stage5 + only CUT&Tag | fail: no eligible |
| M5 | stage5 + only ChIP-seq broad | fail: no eligible (broad skipped) |
| M6 | stage5 + ChIP-seq narrow × 3 bioreps | fail: exactly 2 required |
| M7 | atac_narrow + ATAC narrow + ChIP-seq narrow | pass |
| M8 | stage5 + atac_narrow + both eligible | pass |
| M9 | stage5 + valid CS-narrow + invalid CS-narrow (3 bioreps) | fail: exactly 2 |
| M10 | stage5 + ChIP-seq narrow + ChIP-seq broad | pass |

## 5. Non-goals

- No Snakemake rule, target, config, manifest, or policy changes
- No relaxation of IDR validation for eligible experiments
- No Co-Authored-By

## 6. Verification

```bash
python3 test/test_stage58_mixed_idr_validation.py
python3 test/test_stage55_*.py
python3 test/test_stage53_*.py
python3 test/test_stage54_consensus.py
python3 test/test_stage57_shell_safety.py
python3 test/test_stage28_release_readiness.py
python3 test/test_no_hardcoded_paths.py
git diff --check
```
