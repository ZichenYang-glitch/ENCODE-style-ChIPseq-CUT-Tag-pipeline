# Stage 65 Implementation Plan: Broad IDR Experimental Runtime

**Date:** 2026-06-23
**Status:** Implemented 2026-06-23
**Depends on:** Stage 53/55/62/63/64
**Design spec:** `docs/superpowers/specs/2026-06-23-stage65-broad-idr-experimental-design.md`

## Files

### Create

| File | Purpose |
|------|---------|
| `workflow/rules/idr_broad.smk` | 7 broad IDR rules + `_broad_idr_` helpers (shared, assay-dispatching) |
| `scripts/broad_idr_summary.py` | 15-column broad-peak IDR QC summary |
| `test/test_stage65_broad_idr_config_validation.py` | 14 config validation tests |
| `test/test_stage65_broad_idr_dryrun.py` | 13 dry-run DAG tests |
| `test/test_stage65_broad_idr_summary.py` | 8 summary + output contract tests |

### Modify

| File | Change |
|------|--------|
| `workflow/Snakefile` | `BROAD_CHIPSEQ_IDR_ENABLED`, `BROAD_CUTTAG_IDR_ENABLED`, include, rule all, validator flags |
| `workflow/rules/metadata.smk` | Broad IDR experiment lists (both assays) |
| `workflow/rules/targets.smk` | `_broad_idr_targets()` |
| `workflow/rules/consensus.smk` | Final method/output/target suppression for broad |
| `scripts/validate_samples.py` | Broad IDR eligibility, updated experimental warning, IDR settings gate extension |

### NOT modified

- `workflow/rules/idr.smk`
- `workflow/rules/idr_atac.smk`
- `workflow/rules/idr_cuttag.smk`
- `scripts/stage5b_summary.py`
- `scripts/atac_idr_summary.py`
- `scripts/cuttag_idr_summary.py`
- `scripts/compute_consensus.py`
- `workflow/rules/qc.smk`

## Implementation Steps

### Step 1: Snakefile config gates
Add `BROAD_CHIPSEQ_IDR_ENABLED` and `BROAD_CUTTAG_IDR_ENABLED`.
Pass both to validator.

### Step 2: Config validation
Update warning text. Add eligibility checks for both broad flags.
Extend IDR settings validation gate.

### Step 3: Metadata lists
Add broad IDR experiment lists for both chipseq broad and cuttag broad.

### Step 4: idr_broad.smk
Create 7-rule file with assay dispatch for MACS3 args.
Both use `--broad --broad-cutoff`, layout-aware format, relaxed `-p`, no Tn5 shift.
IDR commands use `--input-file-type broadPeak`.
Summary rules are assay-specific (`broad_idr_chipseq_summary`,
`broad_idr_cuttag_summary`) with `wildcard_constraints` for path
disambiguation.

### Step 5: broad_idr_summary.py
Create script following ATAC/CUT&Tag pattern with broadPeak-aware text.

### Step 6: Targets
Add `_broad_idr_targets()` following established pattern.

### Step 7: Consensus interaction
Update `_consensus_final_method()`, `_consensus_final_output()`,
`_consensus_targets()` for both broad modes.

### Step 8: Snakefile wiring
Add include, rule all targets.

### Step 9: Tests
Write 35 tests across 3 test files.
Include dedicated 17-column broadPeak output contract test.

### Step 10: Verify
Run full regression suite (~281 tests).

## Non-goals

- No SEACR IDR
- No narrow IDR changes
- No legacy `06_idr/` changes
- No Artifact runtime
- No Co-Authored-By
