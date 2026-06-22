# Stage 55 Implementation Plan: ATAC Narrow IDR

**Date:** 2026-06-19
**Status:** Ready for implementation
**Depends on:** Stage 53 (reproducibility policy), Stage 54 (consensus engine)
**Blocks:** Stage 56 (CUT&Tag narrow IDR)
**Design spec:** `docs/superpowers/specs/2026-06-19-stage55-atac-narrow-idr-design.md`

## Context

Stage 53 defined the reproducibility policy. Stage 55 adds ATAC narrow IDR
as the first consumer of `reproducibility.idr.*` — a parallel path to the
legacy Stage 5 ChIP-seq narrow IDR, using the Stage 53 `06_reproducibility/idr/`
namespace.

## Architecture

Parallel ATAC IDR path. Legacy `idr.smk` is not modified. ATAC IDR gets its
own rules file (`idr_atac.smk`), summary script (`scripts/atac_idr_summary.py`),
metadata structures (`ATAC_IDR_*`), target expansion (`_atac_idr_targets()`),
and output namespace.

## Files

### To create

| File | Purpose |
|------|---------|
| `workflow/rules/idr_atac.smk` | ATAC IDR rules (~200 lines) |
| `scripts/atac_idr_summary.py` | ATAC IDR summary script (~120 lines) |
| `test/test_stage55_config_validation.py` | Config validation tests |
| `test/test_stage55_stage5_invariant.py` | Legacy Stage 5 invariance |
| `test/test_stage55_atac_idr_dryrun.py` | Dry-run DAG tests |
| `test/test_stage55_atac_idr_summary.py` | ATAC IDR summary script tests |

### To modify

| File | Change |
|------|--------|
| `workflow/Snakefile` | REPRO_ENABLED / ATAC_IDR_ENABLED; pass to load_and_validate_samples; rule all; include idr_atac.smk |
| `workflow/rules/metadata.smk` | ATAC_IDR_* derived structures + pseudorep expansion lists |
| `workflow/rules/targets.smk` | `_atac_idr_targets()`; manifest dependency targets |
| `scripts/validate_samples.py` | validate_config: reorder, atac_idr_enabled gate, stage4b requirement, idr settings preservation; validate_replicate_groups: atac idr check; CLI call |

### Not modified

- `workflow/rules/idr.smk`
- `scripts/stage5b_summary.py`
- `config/config.yaml`
- `workflow/schemas/config.schema.yaml`
- `docs/reproducibility-policy.md`
- Any Stage 53/54 files

## Implementation steps

### Step 1: Update validate_config (validate_samples.py)

- Move reproducibility validation before the idr-settings gate
- Compute `atac_idr_enabled` from validated `reproducibility` dict
- Add stage4b requirement for ATAC IDR
- Gate `_validate_idr_settings` on `stage5 or atac_idr_enabled`
- Return `seed: 42` in the default idr dict

### Step 2: Update validate_replicate_groups + call chain

- Add `reproducibility_idr_atac_narrow` parameter
- Add ATAC narrow IDR eligibility check
- Update `load_and_validate_samples()` signature
- Update `Snakefile` call
- Update CLI call in validate_samples.py

### Step 3: Update Snakefile

- Define `REPRO_ENABLED`, `ATAC_IDR_ENABLED` from `VALIDATED_CONFIG`
- Pass `reproducibility_idr_atac_narrow=ATAC_IDR_ENABLED`
- Add `+ _atac_idr_targets()` to `rule all`
- Add `include: "rules/idr_atac.smk"` after `idr.smk`

### Step 3b: Update targets.smk (manifest dependency)

- Add `+ _atac_idr_targets()` inside `_manifest_dependency_targets()`

### Step 4: Add metadata structures (metadata.smk)

- Add `ATAC_IDR_EXPERIMENTS`, `ATAC_IDR_BIOREP_*` lists
- Add pseudorep expansion lists (`ATAC_IDR_SPLIT_*`, `ATAC_IDR_PR_PEAK_*`, `ATAC_IDR_SELF_*`)

### Step 5: Add target expansion (targets.smk)

- Implement `_atac_idr_targets()` using `ATAC_IDR_*` lists and `06_reproducibility/idr/` and `06_reproducibility/final/` paths

### Step 6: Create idr_atac.smk

- Mirror legacy `idr.smk` structure with `atac_` prefix rules and helpers
- **All helper functions use `_atac_` prefix** to avoid Snakemake shared-namespace
  collisions: `_atac_idr_biorep_peaks_inputs`, `_atac_idr_macs3_args`,
  `_atac_idr_peak_input`, `_atac_split_input`, `_atac_idr_pseudorep_inputs`,
  `_atac_self_thresh_path`
- **Never** define or overwrite legacy names: `_idr_biorep_peaks_inputs`,
  `_idr_macs3_args`, `_idr_peak_input`, `_split_input`,
  `_idr_pseudorep_inputs`, `_self_thresh_path`
- Use `06_reproducibility/idr/` for intermediate outputs
- Use `06_reproducibility/final/` for final outputs
- Final peak: `<exp>.atac.macs3.narrow.replicate_validated.idr.narrowPeak`
- Reuse `IDR_THRESHOLD`, `IDR_RANK`, `IDR_SEED` from Snakefile
- Call `scripts/atac_idr_summary.py` from the summary rule

### Step 7: Create scripts/atac_idr_summary.py

- CLI: `--true-peaks`, `--pooled-peaks`, `--self1-peaks`, `--self2-peaks`,
  `--experiment`, `--assay`, `--caller`, `--peak-mode`, `--bio-rep-a`,
  `--bio-rep-b`, `--final-method`, `--final-output`, `--output-tsv`
- Compute rescue ratio and self-consistency ratio
- Copy true-replicate thresholded IDR to final output
- Write summary TSV with 15 columns:
  experiment, assay, peak_mode, caller, bio_rep_a, bio_rep_b,
  true_peaks_Nt, pooled_peaks_Np, self1_peaks_N1, self2_peaks_N2,
  rescue_ratio, self_consistency_ratio, reproducibility_status,
  final_method, final_output

### Step 8: Write tests

| File | Tests |
|------|-------|
| `test/test_stage55_config_validation.py` | V1-V9c (12 tests) |
| `test/test_stage55_stage5_invariant.py` | I1-I4 (4 tests) |
| `test/test_stage55_atac_idr_dryrun.py` | D1-D5 (5 tests) |
| `test/test_stage55_atac_idr_summary.py` | S1-S4 (4 tests) |

### Step 9: Run regression tests

```bash
python3 test/test_stage5a_stress.py
python3 test/test_stage5b_stress.py
python3 test/test_stage53_config_validation.py
python3 test/test_stage53_stage5_invariant.py
python3 test/test_stage53_pooled_not_validated.py
python3 test/test_stage53_experimental_warnings.py
python3 test/test_stage53_reproducibility_policy_contract.py
python3 test/test_stage53_output_path_templates.py
python3 test/test_stage54_consensus.py
python3 test/test_stage28_release_readiness.py
python3 test/test_no_hardcoded_paths.py
git diff --check
```

## Key design decisions

1. Parallel ATAC IDR path, not generalization of legacy rules
2. Own summary script (`scripts/atac_idr_summary.py`), not modifying `stage5b_summary.py`
3. Permissive mixed-run validation: non-ATAC experiments skipped silently;
   ATAC broad remains rejected by baseline assay validation
4. ATAC IDR requires stage4b=true
5. IDR settings validated when stage5 or ATAC IDR enabled
6. Stage 53 namespace: `06_reproducibility/idr/` + `06_reproducibility/final/`
7. Final peak copied from true-replicate IDR thresholded
8. Legacy `idr.smk` completely untouched

## Non-goals

- No CUT&Tag IDR (Stage 56)
- No broad IDR (Stage 57)
- No SEACR IDR
- No consensus DAG integration
- No artifact runtime adoption
- No modification of legacy `idr.smk`
- No modification of `scripts/stage5b_summary.py`

## Verification

```bash
cd /home/irenadler/workflow/chipseq
python3 test/test_stage55_config_validation.py
python3 test/test_stage55_stage5_invariant.py
python3 test/test_stage55_atac_idr_dryrun.py
python3 test/test_stage55_atac_idr_summary.py
python3 test/test_stage5a_stress.py
python3 test/test_stage5b_stress.py
# also: all Stage 53 + Stage 54 + release tests
git diff --check
```
