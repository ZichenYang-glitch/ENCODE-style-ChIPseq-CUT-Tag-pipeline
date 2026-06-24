# Phase 1 Finish Implementation Plan

**Date:** 2026-06-24
**Status:** Ready for implementation
**Depends on:** Phase 1 IDR consolidation (PR #39)
**Blocks:** Phase 2
**Design spec:** `docs/superpowers/specs/2026-06-24-phase1-finish-design.md`

## Context

Phase 1 left three categories of duplication:

1. `metadata.smk` IDR list builders (Step A — partially done)
2. ATAC/CUT&Tag/broad IDR rule files (Step B)
3. IDR summary scripts (Step C)

This plan completes all three steps on branch `refactor/phase1-finish`.

## Files

### To create

| File | Purpose |
|------|---------|
| `workflow/rules/idr_reproducibility.smk` | Unified ATAC/CUT&Tag/broad IDR rules |
| `scripts/idr_reproducibility_summary.py` | Unified IDR summary script |

### To modify

| File | Change |
|------|--------|
| `workflow/rules/metadata.smk` | Replace 4 IDR list blocks with shared helpers (Step A, done) |
| `workflow/rules/idr_paths.smk` | Add `_idr_filename_infix(assay, peak_suffix)` helper |
| `workflow/rules/idr_atac.smk` | Deprecated; replaced by `idr_reproducibility.smk` |
| `workflow/rules/idr_cuttag.smk` | Deprecated; replaced by `idr_reproducibility.smk` |
| `workflow/rules/idr_broad.smk` | Deprecated; replaced by `idr_reproducibility.smk` |
| `workflow/Snakefile` | Update includes to use `idr_reproducibility.smk` |
| `scripts/atac_idr_summary.py` | Thin wrapper to unified script |
| `scripts/cuttag_idr_summary.py` | Thin wrapper to unified script |
| `scripts/broad_idr_summary.py` | Thin wrapper to unified script |

### Not modified

- `workflow/rules/idr.smk` (legacy Stage 5 IDR remains separate)
- `scripts/stage5b_summary.py` (legacy Stage 5 summary; incompatible TSV schema)
- `workflow/rules/targets.smk` (target list contract already locked)
- `test/test_idr_paths.py` (contract tests preserved)
- `test/test_targets.py` (target list contract preserved)
- `test/test_dag_snapshots.py` (DAG snapshots preserved)

## Implementation steps

### Step A: Parameterize `metadata.smk` IDR list builders

- [x] Add `_idr_exp_matches(exp, assay, peak_mode)`
- [x] Add `_build_idr_experiment_lists(enabled, assay, peak_mode, exps)`
- [x] Add `_build_idr_pseudorep_lists(experiments)`
- [x] Add `_build_broad_idr_pseudorep_lists(experiments, assays)`
- [x] Replace legacy/ATAC/CUT&Tag/broad blocks with helper calls
- [x] Commit `refactor: parameterize IDR experiment list builders in metadata.smk`
- [x] Verify

```bash
python3 -m pytest test/test_dag_snapshots.py test/test_idr_paths.py test/test_targets.py -v
```

### Step B: Unify ATAC/CUT&Tag/broad IDR rule files

- [x] Add `_idr_filename_infix(assay, peak_suffix)` to `workflow/rules/idr_paths.smk`
- [x] Update `idr_repro_peak_input`, `idr_self_thresh_path`, `idr_pooled_peak_input`, `idr_pooled_thresh_path` to use it
- [x] Create `workflow/rules/idr_reproducibility.smk` with parameterized rules
  - Narrow mode (`*_narrow`): `idr_macs3_biorep_narrow`, `idr_true_replicates_narrow`, `idr_split_pseudoreps_narrow`, `idr_macs3_pseudorep_narrow`, `idr_self_pseudoreps_narrow`, `idr_pooled_pseudoreps_narrow`
  - Broad mode (`*_broad`): `idr_macs3_biorep_broad`, `idr_true_replicates_broad`, `idr_split_pseudoreps_broad`, `idr_macs3_pseudorep_broad`, `idr_self_pseudoreps_broad`, `idr_pooled_pseudoreps_broad`
  - Per-assay summary rules: `idr_summary_atac_narrow`, `idr_summary_cuttag_narrow`, `idr_summary_chipseq_broad`, `idr_summary_cuttag_broad`
  - Suffix mapping (`narrowPeak` / `broadPeak`) is hardcoded per rule because the raw `.txt` outputs cannot carry a `{peak_mode}` wildcard.
  - Rule names change (e.g., `atac_macs3_idr_biorep` → `idr_macs3_biorep_narrow`). DAG snapshots are unaffected because the 8 smoke profiles do not enable modern IDR modes.
- [x] Update `workflow/Snakefile` includes to use `idr_reproducibility.smk`
- [x] Commit `refactor: unify ATAC/CUT&Tag/broad IDR rules into idr_reproducibility.smk`
- [x] Run focused tests

```bash
python3 -m pytest test/test_dag_snapshots.py test/test_idr_paths.py test/test_targets.py -v
```

- [x] Verify DAG snapshot diff is zero
- [x] Remove `workflow/rules/idr_atac.smk`, `idr_cuttag.smk`, `idr_broad.smk`
- [x] Commit `refactor: remove deprecated ATAC/CUT&Tag/broad IDR rule files`
- [x] Run focused tests again

### Step C: Unify ATAC/CUT&Tag/broad IDR summary scripts

- [x] Create `scripts/idr_reproducibility_summary.py` with `--assay`/`--peak-mode` flags
- [x] Convert the three modern summary scripts (`atac_idr_summary.py`, `cuttag_idr_summary.py`, `broad_idr_summary.py`) to thin wrappers
- [x] Update unified rules to call the unified script directly
- [x] Commit `refactor: unify IDR reproducibility summary scripts`
- [x] Run focused tests

```bash
python3 -m pytest test/test_dag_snapshots.py test/test_idr_paths.py test/test_targets.py -v
```

### Final verification

- [x] Run full test suite

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest -q
```

- [x] Run whitespace check

```bash
git diff --check
```

## Key design decisions

1. Keep legacy `idr.smk` separate; unify only the three modern `06_reproducibility/idr/` files.
2. Use wildcard constraints `assay = "chipseq|atac|cuttag"` and `peak_mode = "narrow|broad"`.
3. Preserve exact output paths and filenames.
4. Preserve exact rule shell logic from the modern files.
5. Convert old summary scripts to thin wrappers, not delete them immediately.
6. Add `_idr_filename_infix()` helper to eliminate duplicated assay/peak_suffix dispatch.

## Non-goals

- No changes to legacy `idr.smk` rules or outputs.
- No changes to IDR thresholds, p-values, broad cutoffs, or TSV schema.
- No new dependencies.
- No changes to target list ordering.
- No push to main without explicit approval.

## Verification

```bash
git checkout refactor/phase1-finish
python3 -m pytest test/test_dag_snapshots.py test/test_idr_paths.py test/test_targets.py -v
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest -q
git diff --check
```
