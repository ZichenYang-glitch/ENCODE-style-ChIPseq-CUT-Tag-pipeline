# Stage 66: Reproducibility Output Contract / Manifest Sync — Design Spec

**Date:** 2026-06-23
**Status:** Design approved; implementation pending
**Author:** YangZiChen-glitch (Kaslana)
**Depends on:** Stages 53-65 (all reproducibility DAG), Stages 48-50 (contract infrastructure)
**Branch note:** Stage 66 depends on Stage 65 being merged and should be implemented from a fresh branch based on main.

## 1. Purpose

The `06_reproducibility/` namespace has been built up across Stages 53-65 but
was explicitly deferred from output-contract.md, make_manifest.py, and
artifact-inventory.yaml in every stage spec. Stage 66 closes that gap —
documentation and manifest only, zero runtime changes.

## 2. Scope

Add catalog entries for consensus, IDR final, consensus final, and SEACR
reproducibility outputs to:

- `docs/output-contract.md`
- `scripts/make_manifest.py`
- `docs/architecture/artifact-inventory.yaml`

Also wire `_consensus_targets()`, `_cuttag_idr_targets()`, and
`_broad_idr_targets()` into `_manifest_dependency_targets()` so
`result_manifest.tsv` waits for reproducibility outputs to complete.

No new Snakemake rules. No rule logic changes. No new runtime outputs.
`_manifest_dependency_targets()` wiring only.

## 3. Output Type Catalog

All output types are mode-specific, matching the Stage 48/49/50 contract
requirement for unique `manifest_output_type` / `output_type` identifiers.

Tables use the existing output-contract.md format: `output_type | method |
rule | path | status`. Gating conditions are documented in the Gating
Conditions section of output-contract.md. Paths use `<exp>` (Stage 50
normalizes `<exp>` → `<experiment>`).

### 3.1 Consensus outputs (catalog/report)

Consensus peak + summary are always emitted as catalog/report outputs,
regardless of whether IDR or consensus is primary final.

| output_type | method | rule | path | status |
|-------------|--------|------|------|--------|
| `chipseq_macs3_narrow_consensus_peak` | compute_consensus.py | consensus_compute_narrow | `results/experiments/<exp>/06_reproducibility/consensus/<exp>.chipseq.macs3.narrow.consensus.narrowPeak` | implemented |
| `chipseq_macs3_narrow_consensus_summary` | compute_consensus.py | consensus_compute_narrow | `results/experiments/<exp>/06_reproducibility/consensus/<exp>.chipseq.macs3.narrow.consensus.summary.tsv` | implemented |
| `chipseq_macs3_broad_consensus_peak` | compute_consensus.py | consensus_compute_broad | `results/experiments/<exp>/06_reproducibility/consensus/<exp>.chipseq.macs3.broad.consensus.broadPeak` | implemented |
| `chipseq_macs3_broad_consensus_summary` | compute_consensus.py | consensus_compute_broad | `results/experiments/<exp>/06_reproducibility/consensus/<exp>.chipseq.macs3.broad.consensus.summary.tsv` | implemented |
| `cuttag_macs3_narrow_consensus_peak` | compute_consensus.py | consensus_compute_narrow | `results/experiments/<exp>/06_reproducibility/consensus/<exp>.cuttag.macs3.narrow.consensus.narrowPeak` | implemented |
| `cuttag_macs3_narrow_consensus_summary` | compute_consensus.py | consensus_compute_narrow | `results/experiments/<exp>/06_reproducibility/consensus/<exp>.cuttag.macs3.narrow.consensus.summary.tsv` | implemented |
| `cuttag_macs3_broad_consensus_peak` | compute_consensus.py | consensus_compute_broad | `results/experiments/<exp>/06_reproducibility/consensus/<exp>.cuttag.macs3.broad.consensus.broadPeak` | implemented |
| `cuttag_macs3_broad_consensus_summary` | compute_consensus.py | consensus_compute_broad | `results/experiments/<exp>/06_reproducibility/consensus/<exp>.cuttag.macs3.broad.consensus.summary.tsv` | implemented |
| `atac_macs3_narrow_consensus_peak` | compute_consensus.py | consensus_compute_narrow | `results/experiments/<exp>/06_reproducibility/consensus/<exp>.atac.macs3.narrow.consensus.narrowPeak` | implemented |
| `atac_macs3_narrow_consensus_summary` | compute_consensus.py | consensus_compute_narrow | `results/experiments/<exp>/06_reproducibility/consensus/<exp>.atac.macs3.narrow.consensus.summary.tsv` | implemented |
| `cuttag_seacr_consensus_peak` | compute_consensus.py | consensus_compute_seacr | `results/experiments/<exp>/06_reproducibility/consensus/<exp>.cuttag.seacr.<mode>.consensus.bed` | implemented |
| `cuttag_seacr_consensus_summary` | compute_consensus.py | consensus_compute_seacr | `results/experiments/<exp>/06_reproducibility/consensus/<exp>.cuttag.seacr.<mode>.consensus.summary.tsv` | implemented |

Gating: `reproducibility.enabled: true`, `consensus.enabled: true`, `stage4b: true`, ≥2 biological replicates per experiment, assay/peak_mode match.

### 3.2 IDR final outputs (primary replicate-validated)

| output_type | method | rule | path | status |
|-------------|--------|------|------|--------|
| `cuttag_macs3_narrow_idr_final_peak` | cuttag_idr_summary.py | cuttag_idr_summary | `results/experiments/<exp>/06_reproducibility/final/<exp>.cuttag.macs3.narrow.replicate_validated.idr.narrowPeak` | implemented |
| `atac_macs3_narrow_idr_final_peak` | atac_idr_summary.py | atac_idr_summary | `results/experiments/<exp>/06_reproducibility/final/<exp>.atac.macs3.narrow.replicate_validated.idr.narrowPeak` | implemented |
| `chipseq_macs3_broad_idr_final_peak` | broad_idr_summary.py | broad_idr_chipseq_summary | `results/experiments/<exp>/06_reproducibility/final/<exp>.chipseq.macs3.broad.replicate_validated.idr.broadPeak` | implemented |
| `cuttag_macs3_broad_idr_final_peak` | broad_idr_summary.py | broad_idr_cuttag_summary | `results/experiments/<exp>/06_reproducibility/final/<exp>.cuttag.macs3.broad.replicate_validated.idr.broadPeak` | implemented |

Gating: `reproducibility.enabled: true`, respective IDR config flag enabled, `stage4b: true`, assay/peak_mode match, exactly 2 biological replicates.

### 3.3 Consensus final outputs (primary when no IDR)

| output_type | method | rule | path | status |
|-------------|--------|------|------|--------|
| `chipseq_macs3_broad_consensus_final_peak` | cp | consensus_final | `results/experiments/<exp>/06_reproducibility/final/<exp>.chipseq.macs3.broad.replicate_validated.consensus.broadPeak` | implemented |
| `cuttag_macs3_narrow_consensus_final_peak` | cp | consensus_final | `results/experiments/<exp>/06_reproducibility/final/<exp>.cuttag.macs3.narrow.replicate_validated.consensus.narrowPeak` | implemented |
| `cuttag_macs3_broad_consensus_final_peak` | cp | consensus_final | `results/experiments/<exp>/06_reproducibility/final/<exp>.cuttag.macs3.broad.replicate_validated.consensus.broadPeak` | implemented |
| `cuttag_seacr_consensus_final_peak` | cp | consensus_final_seacr | `results/experiments/<exp>/06_reproducibility/final/<exp>.cuttag.seacr.<mode>.replicate_validated.consensus.bed` | implemented |

Gating: `reproducibility.enabled: true`, `consensus.enabled: true`, `stage4b: true`, respective IDR NOT enabled for that mode, ≥2 biological replicates.

### 3.4 IDR QC summaries

| output_type | method | rule | path | status |
|-------------|--------|------|------|--------|
| `atac_macs3_narrow_idr_summary` | atac_idr_summary.py | atac_idr_summary | `results/experiments/<exp>/06_reproducibility/final/reproducibility_summary.tsv` | implemented |
| `cuttag_macs3_narrow_idr_summary` | cuttag_idr_summary.py | cuttag_idr_summary | `results/experiments/<exp>/06_reproducibility/final/reproducibility_summary.tsv` | implemented |
| `chipseq_macs3_broad_idr_summary` | broad_idr_summary.py | broad_idr_chipseq_summary | `results/experiments/<exp>/06_reproducibility/final/reproducibility_summary.tsv` | implemented |
| `cuttag_macs3_broad_idr_summary` | broad_idr_summary.py | broad_idr_cuttag_summary | `results/experiments/<exp>/06_reproducibility/final/reproducibility_summary.tsv` | implemented |

Gating: same as respective IDR final output (shared rule produces both).

### 3.5 NOT included

No ChIP-seq narrow `06_reproducibility/` IDR final entries. ChIP-seq narrow
consensus/report outputs under `06_reproducibility/consensus/` are allowed and
catalogued above. Legacy ChIP-seq narrow IDR remains represented only by
existing `06_idr/` output types: `idr_conservative`, `idr_optimal`,
`idr_reproducibility_summary`.

Intermediate files excluded from catalog:
- `06_reproducibility/idr/idr_peaks/*` — per-biorep IDR-ready MACS3 peaks
- `06_reproducibility/idr/true_replicates/*_idr.txt` — raw IDR output
- `06_reproducibility/idr/self_pseudoreplicates/*` — self-IDR intermediates
- `06_reproducibility/idr/pooled_pseudoreplicates/*` — pooled-IDR intermediates
- `05_pseudorep/*` — pseudorep BAM splits
- `06_reproducibility/consensus/biorep_peaks/*` — per-biorep consensus peaks

## 4. make_manifest.py Design

### 4.1 Config recomputation

`make_manifest.py` cannot rely on Snakemake globals (`REPRO_ENABLED`,
`BROAD_CHIPSEQ_IDR_ENABLED`, etc.). It must recompute eligibility from the
validated config and sample sheet, matching the same logic used by
`validate_samples.py` and `metadata.smk`.

A new module-level helper `_compute_reproducibility_eligibility(config, samples)`
returns a dict:

```python
{
    "stage4b_enabled": bool,
    "repro_enabled": bool,
    "consensus_enabled": bool,
    "atac_narrow_idr": bool,
    "cuttag_narrow_idr": bool,
    "chipseq_broad_idr": bool,
    "cuttag_broad_idr": bool,
    "seacr_enabled": bool,
    "seacr_mode": str,
    # Experiment-eligible experiment IDs:
    "consensus_experiments": {("assay", "peak_mode"): [exp_ids]},
    "atac_idr_experiments": [exp_ids],
    "cuttag_idr_experiments": [exp_ids],
    "chipseq_broad_idr_experiments": [exp_ids],
    "cuttag_broad_idr_experiments": [exp_ids],
    "seacr_consensus_experiments": [exp_ids],
}
```

If `stage4b_enabled` is false, the builder emits **zero** `06_reproducibility`
rows. All other gating is per-mode per-experiment within the builder.

### 4.2 Builder function

`_build_reproducibility_rows(samples, config)` iterates experiments with ≥2
bioreps and emits:

1. **Consensus peak + summary** for each eligible (assay, peak_mode) pair
   (and SEACR) when `repro_enabled and consensus_enabled`.

2. **IDR final peak + idr_summary** for each eligible IDR experiment when the
   respective IDR flag is enabled. Uses mode-specific output_type.

3. **Consensus final peak** for eligible experiments where IDR is NOT the
   primary method. Uses mode-specific output_type.

4. **Experimental broad IDR rows are omitted entirely** when the experimental
   flag is disabled. No `not_applicable` — just skip the row.

### 4.3 Registration

`_build_reproducibility_rows()` is called from `main()` after
`_build_experiment_rows()` and before `_build_project_rows()`.

## 5. output-contract.md Design

Add a new section "Reproducibility outputs (Stage 53-65)" with subsections:

- Consensus outputs (catalog) — 12 output types
- IDR final outputs (primary) — 4 output types
- Consensus final outputs (primary) — 4 output types
- IDR QC summaries — 4 output types

Each row uses the existing format: `output_type | method | rule | path | status`.
Gating conditions are added to the existing Gating Conditions table.
Paths use `<exp>`.

## 6. artifact-inventory.yaml Design

Add ~24 entries with unique IDs matching output_type names. Each entry uses
the existing 13 Artifact fields — no new fields. Use existing `assay_gate`
values only: `all`, `peak_centric`, `mnase`, `cuttag`, `chipseq`, `atac`,
`idr`. Do not modify `workflow/lib/artifact.py`.

Example (all 13 fields):
```yaml
- id: chipseq_macs3_broad_consensus_peak
  description: ChIP-seq broad consensus peak (catalog/report)
  scope: experiment
  level: pooled_experiment
  assay_gate: peak_centric
  path_template: results/experiments/<experiment>/06_reproducibility/consensus/<experiment>.chipseq.macs3.broad.consensus.broadPeak
  producing_rule: consensus_compute_broad
  tool: compute_consensus.py
  manifest_output_type: chipseq_macs3_broad_consensus_peak
  pipeline_done: false
  rule_all: true
  config_gate: reproducibility.enabled and reproducibility.consensus.enabled
  notes: Consensus catalog output; not primary final unless broad IDR is disabled
```

Legacy `06_idr/` entries are already present and unchanged.

### 6.1 Producing rule names

Consensus rules:
- `consensus_compute_narrow` — narrowPeak consensus
- `consensus_compute_broad` — broadPeak consensus
- `consensus_compute_seacr` — SEACR BED consensus
- `consensus_final` — MACS3 consensus final promotion
- `consensus_final_seacr` — SEACR consensus final promotion

IDR summary/final rules:
- `atac_idr_summary` — ATAC narrow IDR final peak + summary
- `cuttag_idr_summary` — CUT&Tag narrow IDR final peak + summary
- `broad_idr_chipseq_summary` — ChIP-seq broad IDR final peak + summary
- `broad_idr_cuttag_summary` — CUT&Tag broad IDR final peak + summary

### 6.2 Stage4b gating

All `06_reproducibility/` manifest rows require `stage4b: true` because the
DAG outputs depend on experiment/biorep structures. `make_manifest.py`
checks `stage4b_enabled` in `_compute_reproducibility_eligibility()` and
emits zero rows when false.

## 7. targets.smk Change

Only `_manifest_dependency_targets()` is modified. Add three lines:

```python
targets += _consensus_targets()
targets += _cuttag_idr_targets()
targets += _broad_idr_targets()
```

`_atac_idr_targets()` is already present. Legacy `_idr_targets()` is already
present. No new targets, no new rules, no rule logic changes.

## 8. Test Plan

### 8.1 New tests

File: `test/test_stage66_reproducibility_manifest.py` (~12 tests)

| # | Test | Checks |
|---|------|--------|
| 1 | `repro.enabled + consensus.enabled` | Mode-specific consensus output_type rows appear |
| 2 | `atac_narrow: true` + ATAC eligible | `atac_macs3_narrow_idr_final_peak` row |
| 3 | `cuttag_narrow: true` + CUT&Tag eligible | `cuttag_macs3_narrow_idr_final_peak` row |
| 4 | `chipseq_broad_experimental: true` + eligible | `chipseq_macs3_broad_idr_final_peak` row |
| 5 | Broad experimental flag false | No broad IDR final row emitted |
| 6 | No `chipseq_macs3_narrow_idr_final_peak` row ever | Assert this output_type never appears |
| 7 | No `replicate_validated.idr` path for chipseq narrow under `06_reproducibility/` | Assert path pattern absent |
| 8 | SEACR enabled + CUT&Tag PE eligible | `cuttag_seacr_consensus_*` rows |
| 9 | `reproducibility.enabled: false` | Zero `06_reproducibility/` rows |
| 10 | Consensus final present when IDR not enabled for that mode | `*_consensus_final_peak` rows |
| 11 | Disabled mode omitted, not `not_applicable` | Row count check |
| 12 | `stage4b: false` | Zero `06_reproducibility/` rows |

### 8.2 Regression tests

| File | Must pass |
|------|-----------|
| `test/test_stage25_manifest_stress.py` | All 18 pass |
| `test/test_stage49_manifest_artifact_contract.py` | No duplicate manifest_output_type; bidirectional coverage |
| `test/test_stage50_output_contract_dry_run.py` | output_type IDs match artifact inventory IDs |
| `test/test_stage62_consensus_dryrun.py` | 13 tests unchanged |
| `test/test_stage63_seacr_consensus_dryrun.py` | 11 tests unchanged |
| `test/test_stage64_cuttag_idr_dryrun.py` | 12 tests unchanged |
| `test/test_stage65_broad_idr_dryrun.py` | 13 tests unchanged |
| `test/test_stage28_release_readiness.py` | 11 tests unchanged |
| `test/test_no_hardcoded_paths.py` | PASS |
| `git diff --check` | Clean |

## 9. File Inventory

### Modified

| File | Change |
|------|--------|
| `docs/output-contract.md` | Add 24 reproducibility output types |
| `scripts/make_manifest.py` | Add `_compute_reproducibility_eligibility()` + `_build_reproducibility_rows()` |
| `docs/architecture/artifact-inventory.yaml` | Add ~24 reproducibility artifact entries |
| `workflow/rules/targets.smk` | `_manifest_dependency_targets()` adds 3 lines |

### Created

| File | Purpose |
|------|---------|
| `test/test_stage66_reproducibility_manifest.py` | ~12 manifest output type tests |

### NOT modified

- All `workflow/rules/*.smk` (except targets.smk manifest wiring)
- `workflow/Snakefile`
- `scripts/compute_consensus.py`
- All `scripts/*_idr_summary.py`
- `scripts/validate_samples.py`
- `docs/reproducibility-policy.md`
- `workflow/lib/artifact.py`

## 10. Non-goals

- No runtime behavior changes
- No new Snakemake rules
- No new conda environments
- No ChIP-seq narrow `06_reproducibility/` IDR final entries (consensus catalog outputs are allowed)
- No intermediate file cataloguing
- No Co-Authored-By

## 11. Verification

```bash
python3 test/test_stage66_reproducibility_manifest.py
python3 test/test_stage25_manifest_stress.py
python3 test/test_stage49_manifest_artifact_contract.py
python3 test/test_stage50_output_contract_dry_run.py
python3 test/test_stage62_consensus_dryrun.py
python3 test/test_stage28_release_readiness.py
python3 test/test_no_hardcoded_paths.py
git diff --check
```
