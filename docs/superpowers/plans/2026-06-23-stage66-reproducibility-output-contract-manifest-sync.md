# Stage 66 Implementation Plan: Reproducibility Output Contract / Manifest Sync

**Date:** 2026-06-23
**Status:** Ready for implementation (pending review)
**Depends on:** Stages 53-65 (reproducibility DAG), Stages 48-50 (contract infrastructure)
**Branch note:** Stage 66 depends on Stage 65 being merged and should be implemented from a fresh branch based on main.
**Design spec:** `docs/superpowers/specs/2026-06-23-stage66-reproducibility-output-contract-manifest-sync-design.md`

## Files

### Modified

| File | Change |
|------|--------|
| `docs/output-contract.md` | Add 24 mode-specific reproducibility output types (consensus/IDR final/consensus final/SEACR); add gating conditions |
| `scripts/make_manifest.py` | Add `_compute_reproducibility_eligibility()` helper + `_build_reproducibility_rows()` builder; register in `main()` |
| `docs/architecture/artifact-inventory.yaml` | Add ~24 reproducibility artifact entries with unique IDs; existing 13 Artifact fields only; no new fields |
| `workflow/rules/targets.smk` | `_manifest_dependency_targets()` adds `_consensus_targets()`, `_cuttag_idr_targets()`, `_broad_idr_targets()` |

### Created

| File | Purpose |
|------|---------|
| `test/test_stage66_reproducibility_manifest.py` | ~12 manifest output type tests |

### NOT modified

- All `workflow/rules/*.smk` (except targets.smk manifest wiring only)
- `workflow/Snakefile`
- `scripts/compute_consensus.py`
- All `scripts/*_idr_summary.py`
- `scripts/validate_samples.py`
- `docs/reproducibility-policy.md`
- `workflow/lib/artifact.py`

## Implementation Steps

### Step 1: output-contract.md
Add new "Reproducibility outputs (Stage 53-65)" section with 4 subsections:
consensus (12 types), IDR final (4 types), consensus final (4 types),
IDR summaries (4 types). Total: 24 new output_type rows using existing
format (`output_type | method | rule | path | status`). Paths use `<exp>`.
Add gating conditions for each group.

### Step 2: artifact-inventory.yaml
Add ~24 artifact entries with unique IDs matching output_type names.
Use existing 13 Artifact fields only. No `category` or `tags` fields.
Use existing `assay_gate` values. Path templates use `<experiment>`.
Do not modify `workflow/lib/artifact.py`.

### Step 3: make_manifest.py
Add `_compute_reproducibility_eligibility()` that recomputes eligibility
from validated config + sample sheet, including `stage4b_enabled`.
If `stage4b_enabled` is false, emit zero `06_reproducibility` rows.
Add `_build_reproducibility_rows()` builder emitting mode-specific output types.
Omit disabled modes entirely (no `not_applicable`).
Register in `main()`.

### Step 4: targets.smk
Add three lines to `_manifest_dependency_targets()`:
`_consensus_targets()`, `_cuttag_idr_targets()`, `_broad_idr_targets()`.

### Step 5: Tests
Create `test/test_stage66_reproducibility_manifest.py` with ~12 tests.
Key assertions: mode-specific output types, no `chipseq_macs3_narrow_idr_final_peak`,
no chipseq narrow `replicate_validated.idr` under `06_reproducibility/`,
stage4b:false emits zero rows, disabled modes omitted not `not_applicable`.

### Step 6: Verify
Run full contract/manifest regression suite (~246 tests).

## Non-goals

- No runtime behavior changes
- No new Snakemake rules
- No ChIP-seq narrow `06_reproducibility/` IDR final entries (consensus catalog outputs are allowed)
- No intermediate file cataloguing
- No Co-Authored-By
