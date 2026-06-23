# Stage 62: Consensus DAG Integration — Design Spec

**Date:** 2026-06-23
**Status:** Design approved; implementation pending
**Author:** YangZiChen-glitch (Kaslana)
**Depends on:** Stage 54 (consensus engine), Stage 61 (policy alignment)
**Blocks:** Stage 63 (SEACR consensus)
**SEACR consensus:** Deferred to Stage 63

## 1. Overview

Stage 62 integrates the Stage 54 consensus engine (`scripts/compute_consensus.py`)
into the Snakemake DAG for all five MACS3 peak-calling modes. It adds per-biorep
standard MACS3 peak calling, consensus computation, and final-output promotion
where consensus is the primary reproducibility method.

### 1.1 Modes in scope

| # | Assay | Peak mode | Caller | Consensus final? |
|---|-------|-----------|--------|-----------------|
| 1 | chipseq | narrow | MACS3 | No — legacy stage5 IDR is final |
| 2 | chipseq | broad | MACS3 | Yes — consensus final |
| 3 | cuttag | narrow | MACS3 | Yes — until CUT&Tag IDR implemented |
| 4 | cuttag | broad | MACS3 | Yes — consensus final |
| 5 | atac | narrow | MACS3 | No — policy: consensus under consensus/ only until IDR enabled |

SEACR consensus is explicitly out of scope (Stage 63).

### 1.2 Non-goals

- No changes to `idr.smk`, `idr_atac.smk`, `peaks.smk`, or `replicates.smk`
- No reuse of IDR-relaxed-p-value biorep peaks
- No SEACR consensus
- No manifest/output-contract/artifact inventory updates
- No Artifact runtime adoption
- No new IDR rules
- No Co-Authored-By

## 2. Architecture

### 2.1 Per-biorep standard MACS3 rule

**File:** `workflow/rules/consensus.smk`
**Rule:** `consensus_macs3_biorep_peaks`

One shared rule for all five MACS3 modes. Dispatch via `get_macs3_args()`.

```
Input:
  experiments/{experiment}/02_align/biorep{bio_rep}.final.bam
  experiments/{experiment}/02_align/biorep{bio_rep}.final.bam.bai
  [+ experiments/{experiment}/02_align/{experiment}.pooled.control.final.bam
     if experiment has controls]

Params:
  macs3_args = _consensus_macs3_args(wildcards)
  → resolves first treatment sample's assay, calls get_macs3_args()
  → uses standard -q 0.01, NOT relaxed -p 0.1

Output:
  experiments/{experiment}/06_reproducibility/consensus/
    biorep_peaks/
    {experiment}.biorep{bio_rep}.{assay}.macs3.{peak_mode}_peaks.{suffix}
    suffix = narrowPeak or broadPeak
```

The `_peaks` infix matches MACS3's native naming convention
(`<name>_peaks.narrowPeak`, `<name>_peaks.broadPeak`).

**Control policy:** Same pooled-control policy as IDR biorep input helpers
and `_pooled_peaks_inputs()`: if `exp in POOLED_CONTROL_EXPERIMENTS`, add:
`{OUTDIR}/experiments/{exp}/02_align/{exp}.pooled.control.final.bam`
as `-c <control_bam>` to MACS3. The per-biorep rule does NOT require a
separate per-biorep control — it uses the same pooled control BAM that
legacy MACS3 and IDR rules use.

**MACS3 args wrapper (`_consensus_macs3_args`):**
```python
def _consensus_macs3_args(wildcards):
    experiment = wildcards.experiment
    treatment_ids = TREATMENT_SAMPLES_BY_EXPERIMENT.get(experiment, [])
    if not treatment_ids:
        return ""
    # Use the standard assay-specific args (with -q, not -p)
    class _WC: sample = treatment_ids[0]
    return get_macs3_args(_WC)
```

This calls `get_macs3_args()` which dispatches to the assay-specific function
(chipseq/cuttag/atac) with standard q-value 0.01, Tn5 shift for CUT&Tag/ATAC,
and broad mode for broad peaks.

### 2.2 Consensus compute rule

**Rule:** `consensus_compute`

```
Input: per-biorep peak files for an experiment (gathered by _consensus_peak_inputs)

Shell: python3 scripts/compute_consensus.py \
  --peaks <biorep1.narrowPeak> <biorep2.narrowPeak> \
  --bioreps <br1> <br2> \
  --format narrowPeak|broadPeak \
  --min-replicates {min_replicates} \
  --reciprocal-overlap {reciprocal_overlap} \
  --output <consensus_peak> \
  --summary <consensus_summary> \
  --experiment {experiment} --assay {assay} --caller macs3 --peak-mode {peak_mode} \
  --final-method {method} --final-output {path or ""}

Output:
  experiments/{experiment}/06_reproducibility/consensus/
    {experiment}.{assay}.macs3.{peak_mode}.consensus.{suffix}
  experiments/{experiment}/06_reproducibility/consensus/
    {experiment}.{assay}.macs3.{peak_mode}.consensus.summary.tsv
```

### 2.3 Final consensus rule

**Rule:** `consensus_final`

Runs ONLY for modes where consensus is primary final (chipseq broad, cuttag
narrow, cuttag broad). Does NOT run for chipseq narrow (legacy IDR is final)
or ATAC narrow (policy requires IDR for final/).

```
Input: consensus peak file

Shell: cp {input} {output}

Output:
  experiments/{experiment}/06_reproducibility/final/
    {experiment}.{assay}.macs3.{peak_mode}.replicate_validated.consensus.{suffix}
```

**Gating by mode:**

| Mode | consensus_final runs? | Reason |
|------|----------------------|--------|
| chipseq narrow | No | Legacy stage5 IDR is final |
| chipseq broad | Yes | Consensus is primary final |
| cuttag narrow | Yes | CUT&Tag narrow IDR not yet implemented; policy says consensus final until cuttag_narrow IDR is implemented and enabled |
| cuttag broad | Yes | Consensus is primary final |
| atac narrow | No | Policy: consensus under consensus/ only; final/ waits for ATAC IDR |

For CUT&Tag narrow: when `reproducibility.idr.cuttag_narrow: true` but no
CUT&Tag IDR rules exist, the `consensus_final` rule still runs because no
actual IDR output exists to replace it. When CUT&Tag IDR is implemented in a
future stage, the gating will change to skip `consensus_final` when
`cuttag_narrow: true` and CUT&Tag IDR outputs exist.

### 2.4 Output paths (no duplicates, no collisions)

Per-biorep peak filenames include assay/caller/peak_mode:
```
{experiment}.biorep{bio_rep}.{assay}.macs3.{peak_mode}.{suffix}
```
Example: `exp1.biorep1.chipseq.macs3.narrow.narrowPeak`

This avoids collisions between different assay/caller/mode combinations in the
same experiment.

## 3. Metadata and Eligibility

### 3.1 Experiment lists (metadata.smk)

Gated on `CONSENSUS_ENABLED` and `STAGE4B`. Derived from
`MULTI_BIOREP_EXPERIMENTS`:

```python
CONSENSUS_EXPERIMENTS = {
    ("chipseq", "narrow", "macs3"): [...],
    ("chipseq", "broad", "macs3"): [...],
    ("cuttag", "narrow", "macs3"): [...],
    ("cuttag", "broad", "macs3"): [...],
    ("atac", "narrow", "macs3"): [...],
}
```

Each list contains experiments with ≥2 biological replicates and the matching
(assay, peak_mode) from the first treatment sample. Experiments with 1 biorep
are silently excluded.

Pre-computed expansion lists for `expand()` calls mirror the pattern used by
IDR metadata: `CONSENSUS_BIOREP_EXP_LIST`, `CONSENSUS_BIOREP_LIST` etc.

### 3.2 Config gating (Snakefile)

Reuse existing `REPRO_ENABLED` from Stage 55. Add only `CONSENSUS_ENABLED`:

```python
# Stage 62: consensus gating (reuses REPRO_ENABLED from Stage 55)
CONSENSUS_ENABLED = (
    REPRO_ENABLED
    and REPRODUCIBILITY_CONFIG.get("consensus", {}).get("enabled", True)
)
```

### 3.3 Wildcard validation (consensus.smk)

The per-biorep MACS3 rule validates wildcards to prevent accidental targets:

```python
wildcard_constraints:
    assay = r"chipseq|cuttag|atac"
    peak_mode = r"narrow|broad"
    bio_rep = r"\d+"
```

The consensus compute rule derives `suffix` from `peak_mode`:
- `narrow` → `narrowPeak`
- `broad` → `broadPeak`

### 3.4 final_method/final_output summary semantics

Passed to `compute_consensus.py --final-method` / `--final-output`:

| Mode | Condition | final_method | final_output |
|------|-----------|-------------|-------------|
| chipseq narrow | stage5 true | `idr` | legacy `06_idr/final/conservative.narrowPeak` |
| chipseq narrow | stage5 false | `none` | `""` |
| chipseq broad | any | `consensus` | `06_reproducibility/final/...replicate_validated.consensus.broadPeak` |
| cuttag narrow | any (Stage 62) | `consensus` | `06_reproducibility/final/...replicate_validated.consensus.narrowPeak` |
| cuttag broad | any | `consensus` | `06_reproducibility/final/...replicate_validated.consensus.broadPeak` |
| atac narrow | ATAC IDR enabled | `idr` | `06_reproducibility/final/...replicate_validated.idr.narrowPeak` |
| atac narrow | ATAC IDR disabled | `none` | `""` |

Note: CUT&Tag narrow uses `final_method=consensus` in Stage 62 regardless of
`cuttag_narrow` config value, because CUT&Tag IDR is not yet implemented.
This changes in the future CUT&Tag IDR stage.

## 4. Target Expansion (targets.smk)

`_consensus_targets()` adds to rule all:

1. Consensus peak + summary for every eligible experiment (from `consensus_compute`)
2. Final consensus peak for modes where consensus is primary final (from `consensus_final`)

Per-biorep peak files are NOT added to rule all — they are built as
dependencies of `consensus_compute`. Only consensus peak/summary and final
outputs are rule-all targets.

No manifest dependency wiring. Manifest/output-contract updates are deferred.

## 5. File Inventory

| File | Change |
|------|--------|
| `workflow/rules/consensus.smk` | Create — per-biorep MACS3 + consensus + final rules |
| `workflow/Snakefile` | REPRO_ENABLED, CONSENSUS_ENABLED, include, rule all |
| `workflow/rules/metadata.smk` | CONSENSUS experiment lists + biorep expansion |
| `workflow/rules/targets.smk` | `_consensus_targets()` |
| `test/test_stage62_consensus_dryrun.py` | Create — DAG gating tests |

### Explicitly NOT modified

- `workflow/rules/idr.smk`
- `workflow/rules/idr_atac.smk`
- `workflow/rules/peaks.smk`
- `workflow/rules/replicates.smk`
- `scripts/validate_samples.py`
- `config/config.yaml`
- `docs/reproducibility-policy.md`
- Any manifest/output-contract/artifact files

## 6. Test Plan

**File:** `test/test_stage62_consensus_dryrun.py`

| # | Test | Expected |
|---|------|---------|
| D1 | `reproducibility.enabled: false` | no consensus targets |
| D2 | `reproducibility.enabled: true`, `consensus.enabled: false` | no consensus targets |
| D3 | consensus enabled + 2 chipseq narrow bioreps | per-biorep peaks built as deps, consensus peak/summary target |
| D4 | consensus enabled + 2 chipseq broad bioreps | broadPeak consensus + consensus_final target |
| D5 | consensus enabled + 2 cuttag narrow bioreps | consensus_final target (CUT&Tag IDR not yet implemented) |
| D6 | consensus enabled + 2 cuttag broad bioreps | consensus_final target |
| D7 | consensus enabled + 2 atac narrow bioreps | consensus peak/summary, NO consensus_final target |
| D8 | ChIP-seq narrow + stage5: true | `06_idr/final/` legacy target unchanged, NO `06_reproducibility/final/` consensus |
| D9 | ChIP-seq narrow + stage5: false | consensus generated, NO final target (defer final to IDR when enabled) |
| D10 | 1 biorep experiment | no consensus targets (silently excluded) |
| D11 | No SEACR consensus targets | verify zero SEACR paths under `06_reproducibility/consensus/` |
| D12 | No `06_idr/` peaks used as consensus inputs | verify consensus rule input paths are under `06_reproducibility/consensus/biorep_peaks/` |
| D13 | `idr.smk` and `idr_atac.smk` unchanged | git diff confirms zero changes |
| D14 | CUT&Tag narrow + `cuttag_narrow: true` | consensus_final still exists (CUT&Tag IDR not yet implemented) |
| D15 | CUT&Tag narrow + `cuttag_narrow: false` | consensus_final exists (same as D14 — no runtime IDR to replace it) |
| D16 | All new shell blocks have `set -e -o pipefail` | Stage 57 test passes |

## 7. Regression Tests

```bash
python3 test/test_stage54_consensus.py          # consensus engine
python3 test/test_stage53_reproducibility_policy_contract.py
python3 test/test_stage55_config_validation.py
python3 test/test_stage55_stage5_invariant.py
python3 test/test_stage28_release_readiness.py
python3 test/test_no_hardcoded_paths.py
git diff --check
```

## 8. Deferred to Stage 63+

- SEACR consensus (needs per-biorep SEACR rules)
- Manifest/output-contract update for consensus outputs
- CUT&Tag narrow IDR implementation (changes `consensus_final` gating for cuttag narrow)
