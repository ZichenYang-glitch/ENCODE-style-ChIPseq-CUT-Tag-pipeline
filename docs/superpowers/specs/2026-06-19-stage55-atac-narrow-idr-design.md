# Stage 55: ATAC Narrow IDR — Design Spec

**Date:** 2026-06-19
**Status:** Design approved; implementation pending
**Author:** YangZiChen-glitch (Kaslana)
**Depends on:** Stage 53 (reproducibility policy), Stage 54 (consensus engine)
**Blocks:** Stage 56 (CUT&Tag narrow IDR), Stage 57 (broad experimental IDR)

## 1. Overview

Stage 55 adds ATAC-seq narrowPeak IDR support, controlled by the Stage 53
reproducibility config block. It introduces a parallel ATAC IDR path that
coexists with the legacy Stage 5 ChIP-seq narrow IDR without modifying it.

### 1.1 Motivation

Stage 53 defined the reproducibility policy: ATAC narrow IDR is an established
mode gated behind `reproducibility.idr.atac_narrow: true`. The existing Stage 5
IDR pipeline supports only ChIP-seq narrow. ATAC narrowPeak IDR uses the same
narrowPeak + MACS3 + pseudoreplicate machinery and is well-established in the
field.

### 1.2 Scope

| In scope | Out of scope |
|----------|-------------|
| ATAC narrow IDR (parallel path) | CUT&Tag IDR (Stage 56) |
| `workflow/rules/idr_atac.smk` | Broad IDR (Stage 57) |
| `scripts/atac_idr_summary.py` | SEACR IDR |
| Metadata/target expansion | Consensus DAG integration |
| Config validation | Artifact runtime adoption |
| Dry-run tests | Modifying `idr.smk` |

### 1.3 Architecture Decision

**Parallel ATAC IDR path.** Legacy Stage 5 ChIP-seq IDR is not touched.
ATAC IDR uses its own rules (`idr_atac.smk`), metadata structures
(`ATAC_IDR_*`), target expansion (`_atac_idr_targets()`), summary script
(`scripts/atac_idr_summary.py`), and output namespace
(`06_reproducibility/idr/`).

This minimizes legacy regression risk.

## 2. Config and Validation

### 2.1 Snakefile additions

```python
# Stage 55: reproducibility-driven IDR
REPRODUCIBILITY_CONFIG = VALIDATED_CONFIG.get("reproducibility", {})
REPRO_ENABLED = REPRODUCIBILITY_CONFIG.get("enabled", False)
ATAC_IDR_ENABLED = (
    REPRO_ENABLED
    and REPRODUCIBILITY_CONFIG.get("idr", {}).get("atac_narrow", False)
)
```

Pass into `load_and_validate_samples()`:
```python
SAMPLES = load_and_validate_samples(
    VALIDATED_CONFIG["samples"],
    use_control=USE_CONTROL,
    stage5_enabled=STAGE5,
    reproducibility_idr_atac_narrow=ATAC_IDR_ENABLED,
)
```

`rule all`: append `+ _atac_idr_targets()`.
Include: `include: "rules/idr_atac.smk"` after `idr.smk`.

### 2.2 Config validation (validate_config)

```python
# Validate reproducibility first, then determine atac_idr_enabled
validated["reproducibility"] = _validate_reproducibility(
    config.get("reproducibility", {}), validated
)
repro = validated["reproducibility"]
atac_idr_enabled = (
    repro.get("enabled", False)
    and repro.get("idr", {}).get("atac_narrow", False)
)

# ATAC IDR requires stage4b
if atac_idr_enabled and not validated.get("stage4b", True):
    raise ValidationError(
        "config: reproducibility.idr.atac_narrow=true requires "
        "stage4b=true."
    )

# Validate idr settings when stage5 or ATAC IDR is enabled
if validated["stage5"] or atac_idr_enabled:
    validated["idr"] = _validate_idr_settings(config.get("idr", {}))
else:
    validated["idr"] = {"threshold": 0.05, "rank": "p.value", "seed": 42}
```

### 2.3 Sample validation (validate_replicate_groups)

Updated signatures:
- `load_and_validate_samples(..., reproducibility_idr_atac_narrow=False)`
- `validate_replicate_groups(..., reproducibility_idr_atac_narrow=False)`

```python
if reproducibility_idr_atac_narrow:
    atac_narrow_exps = []
    for exp, rows in exp_treatments.items():
        if len(rows) == 0:
            continue
        first = rows[0]
        if first["assay"] != "atac":
            continue       # skip non-ATAC silently
        if first["peak_mode"] != "narrow":
            continue       # defensive; baseline validation rejects ATAC non-narrow
        bio_reps = sorted({r["biological_replicate"] for r in rows})
        if len(bio_reps) != 2:
            raise ValidationError(
                f"reproducibility.idr.atac_narrow: ATAC narrow "
                f"experiment {exp!r} has {len(bio_reps)} biological "
                f"replicate(s). IDR requires exactly 2."
            )
        atac_narrow_exps.append(exp)

    if not atac_narrow_exps:
        raise ValidationError(
            "reproducibility.idr.atac_narrow is true but no eligible "
            "ATAC narrow experiments with exactly 2 biological "
            "replicates were found."
        )
```

Key behavior: non-ATAC experiments are silently skipped. Only ATAC narrow
experiments are checked for ATAC IDR. ATAC broad remains invalid under the
existing baseline assay policy (`assay=atac` currently supports narrow only).
Mixed sample sheets are supported when `stage5: false`.

### 2.4 Config reference

```yaml
reproducibility:
  enabled: true
  idr:
    atac_narrow: true

idr:           # top-level, validated when atac_narrow or stage5 is true
  threshold: 0.05
  rank: "p.value"
  seed: 42
```

## 3. Metadata and Targets

### 3.1 Metadata (metadata.smk)

After the existing legacy Stage 5 block:

```python
# Stage 55: ATAC narrow IDR derived structures
if ATAC_IDR_ENABLED:
    ATAC_IDR_EXPERIMENTS = []
    ATAC_IDR_BIOREP_EXP_LIST = []
    ATAC_IDR_BIOREP_LIST = []

    for exp in MULTI_BIOREP_EXPERIMENTS:
        bioreps = _bioreps_for(exp, "treatment")
        if len(bioreps) == 2:
            first = SAMPLE_MAP[TREATMENT_SAMPLES_BY_EXPERIMENT[exp][0]]
            if first["assay"] == "atac" and first["peak_mode"] == "narrow":
                ATAC_IDR_EXPERIMENTS.append(exp)
                for br in sorted(bioreps):
                    ATAC_IDR_BIOREP_EXP_LIST.append(exp)
                    ATAC_IDR_BIOREP_LIST.append(br)
else:
    ATAC_IDR_EXPERIMENTS = []
    ATAC_IDR_BIOREP_EXP_LIST = []
    ATAC_IDR_BIOREP_LIST = []
```

Pseudorep split lists mirror the legacy structure with `ATAC_` prefix.

### 3.2 Target expansion (targets.smk)

`_atac_idr_targets()` mirrors `_idr_targets()` but uses `ATAC_IDR_*` lists,
`06_reproducibility/idr/` for intermediates, and `06_reproducibility/final/`
for final outputs.

`_manifest_dependency_targets()` appends `+ _atac_idr_targets()`.

### 3.3 idr.smk

Not modified. Legacy `06_idr/` paths, rule names, and gating preserved.

## 4. Output Paths

### 4.1 Intermediates

Paths use Snakemake wildcards `{experiment}`, `{bio_rep}`, `{source}`, `{pr}`.

```
results/experiments/{experiment}/06_reproducibility/idr/
  idr_peaks/
    {experiment}_atac_biorep{bio_rep}_idr.narrowPeak
    {experiment}_atac_{source}_pr{pr}_idr.narrowPeak
  true_replicates/
    {experiment}_atac_idr.txt
    {experiment}_atac_idr.thresholded.narrowPeak
  self_pseudoreplicates/
    {experiment}_atac_biorep{bio_rep}_idr.txt
    {experiment}_atac_biorep{bio_rep}_idr.thresholded.narrowPeak
  pooled_pseudoreplicates/
    {experiment}_atac_idr.txt
    {experiment}_atac_idr.thresholded.narrowPeak
```

### 4.2 Final outputs

```
results/experiments/{experiment}/06_reproducibility/final/
  {experiment}.atac.macs3.narrow.replicate_validated.idr.narrowPeak
  reproducibility_summary.tsv
```

`final/` is a sibling of `idr/` under `06_reproducibility/`, per Stage 53.

### 4.3 Pseudorep BAMs

```
results/experiments/{experiment}/05_pseudorep/
  {experiment}_atac_{source}.pr{pr}.bam
  {experiment}_atac_{source}.pr{pr}.bam.bai
```

Where `{source}` takes values `biorep{br_label}` and `pooled`.

### 4.4 Final peak semantics

The final `replicate_validated.idr.narrowPeak` is copied from the
true-replicate IDR thresholded narrowPeak. This is the primary validated peak
set — equivalent to the legacy conservative output but using the Stage 53
single final filename. The pooled-pseudoreplicate thresholded IDR is recorded
in the summary TSV for rescue ratio / QC.

### 4.5 Legacy paths (unchanged)

```
results/experiments/{experiment}/06_idr/
  ... (unchanged Stage 5 ChIP-seq narrow outputs)
```

## 5. Rules and Scripts

### 5.1 New rule file: `workflow/rules/idr_atac.smk`

| Rule | Purpose |
|------|---------|
| `atac_macs3_idr_biorep` | IDR-ready MACS3 call on a single ATAC biorep BAM |
| `atac_idr_true_replicates` | IDR between the two ATAC biorep peak sets |
| `atac_split_pseudoreps` | Deterministic pseudorep BAM split (ATAC-specific paths) |
| `atac_macs3_idr_pseudorep` | IDR-ready MACS3 on ATAC pseudorep BAM |
| `atac_idr_self_pseudoreps` | Self-pseudorep IDR per ATAC biorep |
| `atac_idr_pooled_pseudoreps` | Pooled-pseudorep IDR for ATAC |
| `atac_idr_summary` | Reproducibility QC summary + final peak copy |

Rules reuse the same MACS3/IDR shell logic as legacy rules but with ATAC
output paths and the `atac_` name prefix. Conda environments (`macs3.yml`,
`idr.yml`, `samtools.yml`, `python.yml`) are shared.

**Helper function prefixing:** All helper functions in `idr_atac.smk` use the
`_atac_` prefix to avoid Snakemake shared-namespace collisions with legacy
`idr.smk`. Specifically, `idr_atac.smk` defines:

- `_atac_idr_biorep_peaks_inputs`
- `_atac_idr_macs3_args`
- `_atac_idr_peak_input`
- `_atac_split_input`
- `_atac_idr_pseudorep_inputs`
- `_atac_self_thresh_path`

and does NOT define/overwrite any of:

- `_idr_biorep_peaks_inputs`
- `_idr_macs3_args`
- `_idr_peak_input`
- `_split_input`
- `_idr_pseudorep_inputs`
- `_self_thresh_path`

This is required because Snakemake includes share one Python namespace.
Legacy rule lambdas resolve helper names at evaluation time; overwriting
them would silently break Stage 5 behavior.

### 5.2 New summary script: `scripts/atac_idr_summary.py`

Purpose: compute rescue ratio and self-consistency ratio from ATAC IDR
outputs, copy the final validated peak set, and write the reproducibility
summary TSV.

```
python3 scripts/atac_idr_summary.py \
  --true-peaks results/experiments/{exp}/06_reproducibility/idr/true_replicates/{exp}_atac_idr.thresholded.narrowPeak \
  --pooled-peaks results/experiments/{exp}/06_reproducibility/idr/pooled_pseudoreplicates/{exp}_atac_idr.thresholded.narrowPeak \
  --self1-peaks results/experiments/{exp}/06_reproducibility/idr/self_pseudoreplicates/{exp}_atac_biorep{bio_rep_a}_idr.thresholded.narrowPeak \
  --self2-peaks results/experiments/{exp}/06_reproducibility/idr/self_pseudoreplicates/{exp}_atac_biorep{bio_rep_b}_idr.thresholded.narrowPeak \
  --experiment {exp} \
  --assay atac \
  --caller macs3 \
  --peak-mode narrow \
  --bio-rep-a {bio_rep_a} --bio-rep-b {bio_rep_b} \
  --final-method idr \
  --final-output results/experiments/{exp}/06_reproducibility/final/{exp}.atac.macs3.narrow.replicate_validated.idr.narrowPeak \
  --output-tsv results/experiments/{exp}/06_reproducibility/final/reproducibility_summary.tsv
```

Summary TSV columns (15):
```
experiment   assay   peak_mode   caller   bio_rep_a   bio_rep_b
true_peaks_Nt   pooled_peaks_Np   self1_peaks_N1   self2_peaks_N2
rescue_ratio   self_consistency_ratio   reproducibility_status
final_method   final_output
```

Key behavior:
- Final peak file is copied from `--true-peaks`.
- Rescue ratio: `max(Np,Nt) / min(Np,Nt)`.
- Self-consistency ratio: `max(N1,N2) / min(N1,N2)`.
- `reproducibility_status`: `pass` if both ratios < 2 and finite; otherwise `fail`.
- Handles zero denominators gracefully (`NA` or `inf`).

`scripts/stage5b_summary.py` is NOT modified. ATAC gets its own summary script
for clean separation from legacy behavior.

### 5.3 MACS3 args

`_atac_idr_macs3_args()` mirrors `_idr_macs3_args()` — extracts layout and
genome from the first treatment sample, uses
`_tool_param("idr_macs3", "pvalue", 0.1)`.

## 6. Test Plan

### 6.1 Config validation tests

**File:** `test/test_stage55_config_validation.py`

| # | Test | Expected |
|---|------|---------|
| V1 | `atac_narrow: true` + `enabled: true` + 2 ATAC narrow bioreps | passes |
| V2 | `atac_narrow: true` + `enabled: false` | passes config validation; ATAC_IDR_ENABLED resolves to false |
| V3 | `atac_narrow: true` + `enabled: true` + `stage4b: false` | fails: requires stage4b=true |
| V4a | `atac_narrow: true` + ATAC broad only | fails: existing assay policy rejects ATAC broad |
| V4b | `atac_narrow: true` + ChIP-seq broad + eligible ATAC narrow | passes; non-ATAC broad skipped |
| V5 | `atac_narrow: true` + ATAC narrow + 3 bioreps | fails: requires exactly 2 |
| V6 | `atac_narrow: true` + ATAC narrow + 1 biorep | fails: requires exactly 2 |
| V7 | `atac_narrow: true` + no ATAC experiments | fails: no eligible experiments |
| V8 | `stage5: false` + `atac_narrow: true` + mixed ATAC narrow + ChIP-seq narrow | passes; only ATAC checked |
| V9a | `atac_narrow: true` + `stage5: false` → `idr` validated | threshold, rank, seed populated |
| V9b | `atac_narrow: true` + invalid `idr.rank` | fails validation |
| V9c | `atac_narrow: true` + invalid `idr.threshold` | fails validation |

### 6.2 Stage 5 invariance tests

**File:** `test/test_stage55_stage5_invariant.py`

| # | Test | Expected |
|---|------|---------|
| I1 | `stage5: true` + ChIP-seq narrow + 2 bioreps | legacy `06_idr/final/conservative.narrowPeak` target unchanged; `optimal.narrowPeak` target unchanged |
| I2 | `stage5: true` + ATAC experiment present | **fails** under legacy Stage 5 validation (assay must be chipseq); proves legacy behavior unchanged |
| I3 | `stage5: false` + `atac_narrow: true` + mixed ATAC narrow + ChIP-seq narrow | ATAC IDR targets under `06_reproducibility/` appear; ChIP-seq legacy `06_idr/` targets absent; ChIP-seq `06_reproducibility/idr/` targets absent |
| I4 | ATAC IDR targets | never under `06_idr/` |

### 6.3 Dry-run tests

**File:** `test/test_stage55_atac_idr_dryrun.py`

| # | Test | Expected |
|---|------|---------|
| D1 | `atac_narrow: true` + 2 ATAC narrow bioreps | ATAC IDR targets under `06_reproducibility/` |
| D2 | `stage5: false` + `atac_narrow: true` + mixed ATAC narrow + ChIP-seq narrow | ATAC IDR targets appear; ChIP-seq `06_idr/` targets absent; ChIP-seq `06_reproducibility/idr/` targets absent |
| D3 | `atac_narrow: false` | no ATAC IDR targets |
| D4 | `enabled: false` + `atac_narrow: true` | no ATAC IDR targets |
| D5 | ATAC IDR final paths | `06_reproducibility/final/<exp>.atac.macs3.narrow.replicate_validated.idr.narrowPeak`; `06_reproducibility/final/reproducibility_summary.tsv` |

### 6.4 Regression tests (must pass unchanged)

| File | Tests |
|------|-------|
| `test/test_stage55_atac_idr_summary.py` | 4/4 |
| `test/test_stage53_*.py` (6 files) | 50/50 |
| `test/test_stage54_consensus.py` | 31/31 |
| `test/test_stage5a_stress.py` | legacy IDR config pass |
| `test/test_stage5b_stress.py` | legacy pseudorep pass |
| `test/test_stage28_release_readiness.py` | 11/11 |
| `test/test_no_hardcoded_paths.py` | pass |

### 6.5 Deferred

| Item | Reason |
|------|--------|
| ATAC IDR real execution | Requires real ATAC FASTQs + Bowtie2 index |
| CUT&Tag IDR | Stage 56 |
| Broad experimental IDR | Stage 57 |

## 7. Files

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
| `workflow/rules/metadata.smk` | ATAC_IDR_* derived structures |
| `workflow/rules/targets.smk` | `_atac_idr_targets()`; manifest dependency targets |
| `scripts/validate_samples.py` | validate_config reorder + stage4b gate; validate_replicate_groups + atac idr check; CLI call |

### Explicitly NOT modified

- `workflow/rules/idr.smk` (legacy Stage 5 IDR)
- `scripts/stage5b_summary.py` (legacy summary)
- `config/config.yaml`
- `workflow/schemas/config.schema.yaml`
- `docs/reproducibility-policy.md`
- Any Stage 53/54 files

## 8. Non-goals

- No CUT&Tag IDR (Stage 56)
- No broad IDR (Stage 57)
- No SEACR IDR
- No consensus DAG integration
- No artifact runtime adoption
- No modification of legacy `idr.smk`
- No migration of ChIP-seq IDR to `06_reproducibility/`
