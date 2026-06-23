# Stage 64 Implementation Plan: CUT&Tag Narrow IDR Runtime

**Date:** 2026-06-23
**Status:** Implemented 2026-06-23
**Depends on:** Stage 55 (ATAC IDR), Stage 62 (MACS3 consensus)
**Blocks:** CUT&Tag broad IDR (future stage, not scoped)
**Design spec:** `docs/superpowers/specs/2026-06-23-stage64-cuttag-idr-design.md`

## Files

### Create

| File | Purpose |
|------|---------|
| `workflow/rules/idr_cuttag.smk` | 7 CUT&Tag IDR rules + `_cuttag_` helpers |
| `scripts/cuttag_idr_summary.py` | 15-column reproducibility QC summary (mirrors ATAC) |
| `test/test_stage64_cuttag_idr_config_validation.py` | 14 config validation tests |
| `test/test_stage64_cuttag_idr_dryrun.py` | 13 dry-run tests |
| `test/test_stage64_cuttag_idr_summary.py` | 8 summary script tests |

### Modify

| File | Change |
|------|--------|
| `workflow/Snakefile` | `CUTTAG_IDR_ENABLED` (after ATAC), include, rule all, validator flag |
| `workflow/rules/metadata.smk` | `CUTTAG_IDR_*` experiment lists (after ATAC IDR block) |
| `workflow/rules/targets.smk` | `_cuttag_idr_targets()` (after ATAC, before consensus) |
| `workflow/rules/consensus.smk` | Final method/output/target suppression scoped to `CUTTAG_IDR_EXPERIMENTS` |
| `workflow/rules/idr_atac.smk` | Summary-rule wildcard constraint to avoid ambiguity with CUT&Tag's generic summary path |
| `scripts/validate_samples.py` | `reproducibility_idr_cuttag_narrow` parameter + eligibility check |

### NOT modified

- `workflow/rules/idr.smk`
- `scripts/stage5b_summary.py`
- `scripts/atac_idr_summary.py`
- `scripts/compute_consensus.py`
- `workflow/rules/qc.smk`
- `workflow/schemas/config.schema.yaml`

## Implementation Steps

### Step 1: Snakefile config gating
Add `CUTTAG_IDR_ENABLED` after `ATAC_IDR_ENABLED`. Pass to validator.

### Step 2: Config validation
Add `reproducibility_idr_cuttag_narrow` parameter to `validate_replicate_groups()`
and `load_and_validate_samples()`. Add CUT&Tag narrow IDR eligibility check
following ATAC pattern. Extend IDR settings gate.

### Step 3: Metadata lists
Add `CUTTAG_IDR_*` experiment lists in metadata.smk, gated on `CUTTAG_IDR_ENABLED`.

### Step 4: idr_cuttag.smk
Create 7-rule file following `idr_atac.smk` pattern. `_cuttag_idr_macs3_args()`
uses layout-aware format, Tn5 shift, relaxed `-p`, no `-q`.

### Step 5: cuttag_idr_summary.py
Create script mirroring `atac_idr_summary.py` with CUT&Tag text.

### Step 6: Targets
Add `_cuttag_idr_targets()` following `_atac_idr_targets()` pattern.

### Step 7: Consensus interaction
Update `_consensus_final_method()`, `_consensus_final_output()`, and
`_consensus_targets()` to handle CUT&Tag narrow IDR final semantics.

### Step 8: Snakefile wiring
Add include and rule all targets.

### Step 9: Tests
Write 35 tests across 3 test files.

### Step 10: Verify
Run full regression suite (~246 tests).

## Non-goals

- No SEACR IDR
- No broad IDR runtime
- No PE-only gating (SE supported)
- No `consensus.enabled` gating IDR
- No Artifact runtime
- No Co-Authored-By
