# Stage 64: CUT&Tag Narrow IDR Runtime — Design Spec

**Date:** 2026-06-23
**Status:** Implemented 2026-06-23
**Author:** YangZiChen-glitch (Kaslana)
**Depends on:** Stage 55 (ATAC IDR), Stage 61 (policy), Stage 62 (MACS3 consensus), Stage 63 (SEACR consensus)

## 1. Purpose

Implement runtime support for CUT&Tag narrow MACS3 IDR, completing the
Stage 61 policy matrix entry for "supported opt-in IDR." CUT&Tag narrow IDR
becomes final when explicitly enabled via `reproducibility.idr.cuttag_narrow:
true`. Consensus remains as report/fallback under `consensus/`.

This is opt-in only — not production-supported default IDR. Production-supported
IDR modes remain ChIP-seq narrow (legacy stage5) and ATAC narrow. CUT&Tag
narrow is supported opt-in in this stage.

## 2. Architecture

**New file:** `workflow/rules/idr_cuttag.smk` — 7-rule mirror of
`workflow/rules/idr_atac.smk`, with `_cuttag_` helper prefix and `cuttag_`
rule prefix. Self-contained in its own file to:

- Preserve legacy ChIP-seq IDR in `idr.smk`
- Preserve ATAC IDR in `idr_atac.smk`
- Avoid helper namespace collisions with `_atac_` helpers

Also new: `scripts/cuttag_idr_summary.py` — mirrors
`scripts/atac_idr_summary.py` with CUT&Tag-specific CLI text/docstring.

## 3. Rules (7 in idr_cuttag.smk)

| # | Rule | Purpose | Conda |
|---|------|---------|-------|
| 1 | `cuttag_macs3_idr_biorep` | IDR-ready MACS3 on biorep BAM | macs3.yml |
| 2 | `cuttag_idr_true_replicates` | True-replicate IDR (raw + thresholded) | idr.yml |
| 3 | `cuttag_split_pseudoreps` | BAM split via `split_pseudoreps.py` | samtools.yml |
| 4 | `cuttag_macs3_idr_pseudorep` | MACS3 on pseudorep BAMs | macs3.yml |
| 5 | `cuttag_idr_self_pseudoreps` | Self-pseudorep IDR per biorep | idr.yml |
| 6 | `cuttag_idr_pooled_pseudoreps` | Pooled pseudorep IDR | idr.yml |
| 7 | `cuttag_idr_summary` | Reproducibility QC + final peak copy | python.yml |

### 3.1 MACS3 args for IDR-ready calls

`_cuttag_idr_macs3_args()`:

- Layout-aware: `-f BAMPE` for PE, `-f BAM` for SE
- Genome size via `_normalize_genome()`
- Relaxed p-value: `-p <idr_macs3.pvalue>` (default 0.1)
- CUT&Tag narrow Tn5 policy: `--nomodel --shift -100 --extsize 200`
- Includes `idr_macs3.extra_args`
- Never emits `-q`

### 3.2 Output namespace

```
06_reproducibility/idr/
  idr_peaks/{exp}_cuttag_biorep{br}_idr.narrowPeak
  idr_peaks/{exp}_cuttag_{source}_pr{pr}_idr.narrowPeak
  true_replicates/{exp}_cuttag_idr.txt
  true_replicates/{exp}_cuttag_idr.thresholded.narrowPeak
  self_pseudoreplicates/{exp}_cuttag_biorep{br}_idr.txt
  self_pseudoreplicates/{exp}_cuttag_biorep{br}_idr.thresholded.narrowPeak
  pooled_pseudoreplicates/{exp}_cuttag_idr.txt
  pooled_pseudoreplicates/{exp}_cuttag_idr.thresholded.narrowPeak

05_pseudorep/
  {exp}_cuttag_{source}.pr{1,2}.bam{,.bai}

06_reproducibility/final/
  {exp}.cuttag.macs3.narrow.replicate_validated.idr.narrowPeak
  reproducibility_summary.tsv
```

## 4. Config and Eligibility

### 4.1 Config gate

```python
CUTTAG_IDR_ENABLED = (
    REPRO_ENABLED
    and REPRODUCIBILITY_CONFIG.get("idr", {}).get("cuttag_narrow", False)
)
```

Config schema already includes `cuttag_narrow` (Stage 61); no schema change needed.

### 4.2 Eligibility

| Condition | Required |
|-----------|----------|
| `reproducibility.enabled` | true |
| `reproducibility.idr.cuttag_narrow` | true |
| `assay` | cuttag |
| `peak_mode` | narrow |
| Biological replicates | exactly 2 |
| Layout | PE or SE (both supported) |
| `stage4b` | true |

### 4.3 Config validation

`validate_samples.py` additions:

1. `cuttag_idr_enabled` gate check (requires stage4b)
2. Extend IDR settings validation gate: `stage5 or atac_idr_enabled or cuttag_idr_enabled`
3. CUT&Tag narrow IDR eligibility in `validate_replicate_groups()`:
   - Skip non-CUT&Tag silently
   - Skip CUT&Tag broad silently
   - Validate exactly 2 biological replicates per CUT&Tag narrow experiment
   - Fail if no eligible experiment exists
4. New parameter: `reproducibility_idr_cuttag_narrow`

### 4.4 Metadata lists

In `metadata.smk`, gated on `CUTTAG_IDR_ENABLED`:

```
CUTTAG_IDR_EXPERIMENTS, CUTTAG_IDR_BIOREP_EXP_LIST, CUTTAG_IDR_BIOREP_LIST,
CUTTAG_IDR_SPLIT_SOURCE_EXP, CUTTAG_IDR_SPLIT_SOURCE_NAME,
CUTTAG_IDR_PR_PEAK_EXP, CUTTAG_IDR_PR_PEAK_SRC, CUTTAG_IDR_PR_PEAK_PR,
CUTTAG_IDR_SELF_EXP, CUTTAG_IDR_SELF_BR
```

Fallback to empty lists when disabled.

## 5. Consensus Interaction

Three narrow changes to `consensus.smk`:

### 5.1 `_consensus_final_method()`

Add before the existing `return "consensus"`:
```python
if assay == "cuttag" and peak_mode == "narrow":
    if CUTTAG_IDR_ENABLED and experiment in CUTTAG_IDR_EXPERIMENTS:
        return "idr"
    return "consensus"
```

### 5.2 `_consensus_final_output()`

Add before the CUT&Tag narrow consensus case:
```python
if (assay == "cuttag" and peak_mode == "narrow"
    and CUTTAG_IDR_ENABLED
    and experiment in CUTTAG_IDR_EXPERIMENTS):
    return (
        f"{OUTDIR}/experiments/{experiment}/06_reproducibility/final/"
        f"{experiment}.cuttag.macs3.narrow."
        f"replicate_validated.idr.narrowPeak"
    )
```

### 5.3 `_consensus_targets()`

When `CUTTAG_IDR_ENABLED`, exclude experiments in `CUTTAG_IDR_EXPERIMENTS`
from the final consensus target expansion for CUT&Tag narrow. Consensus peak
+ summary remain. All other modes unaffected.

## 6. Summary Script

`scripts/cuttag_idr_summary.py` — mirrors `scripts/atac_idr_summary.py`:

- Same 15-column TSV schema
- Same logic: Nt/Np/N1/N2 counts, rescue/self ratios, pass/fail (ratios < 2 and finite)
- Same `shutil.copyfile` for final peak
- CLI text and docstring say "CUT&Tag" not "ATAC"
- No external dependencies (stdlib: argparse, shutil, sys)

## 7. Snakefile Wiring

```python
# After ATAC_IDR_ENABLED:
CUTTAG_IDR_ENABLED = (
    REPRO_ENABLED
    and REPRODUCIBILITY_CONFIG.get("idr", {}).get("cuttag_narrow", False)
)

# In load_and_validate_samples():
reproducibility_idr_cuttag_narrow=CUTTAG_IDR_ENABLED,

# Include (after idr_atac.smk, before consensus.smk):
include: "rules/idr_cuttag.smk"

# Rule all (after _atac_idr_targets(), before _consensus_targets()):
+ _cuttag_idr_targets()
```

## 8. Final Semantics

| Mode | Condition | Final output | Method |
|------|-----------|-------------|--------|
| cuttag narrow | IDR enabled + eligible | `...replicate_validated.idr.narrowPeak` | idr |
| cuttag narrow | IDR disabled or ineligible | `...replicate_validated.consensus.narrowPeak` | consensus |
| cuttag narrow (any) | — | Consensus peak + summary remain under `consensus/` | consensus (report) |

## 9. File Inventory

### Created
| File | Purpose |
|------|---------|
| `workflow/rules/idr_cuttag.smk` | 7 CUT&Tag IDR rules + `_cuttag_` helpers |
| `scripts/cuttag_idr_summary.py` | Reproducibility QC summary |
| `test/test_stage64_cuttag_idr_config_validation.py` | 14 tests |
| `test/test_stage64_cuttag_idr_dryrun.py` | 13 tests |
| `test/test_stage64_cuttag_idr_summary.py` | 8 tests |
| `docs/superpowers/specs/2026-06-23-stage64-cuttag-idr-design.md` | This spec |
| `docs/superpowers/plans/2026-06-23-stage64-cuttag-idr.md` | Implementation plan |

### Modified
| File | Change |
|------|--------|
| `workflow/Snakefile` | `CUTTAG_IDR_ENABLED`, include, rule all, validator flag |
| `workflow/rules/metadata.smk` | `CUTTAG_IDR_*` experiment lists |
| `workflow/rules/targets.smk` | `_cuttag_idr_targets()` |
| `workflow/rules/consensus.smk` | Final method/output/target suppression |
| `workflow/rules/idr_atac.smk` | Summary-rule wildcard constraint so generic `reproducibility_summary.tsv` is not ambiguous |
| `scripts/validate_samples.py` | CUT&Tag narrow IDR eligibility + validation gate |

### Not touched
- `workflow/rules/idr.smk`
- `scripts/stage5b_summary.py`
- `scripts/atac_idr_summary.py`
- `scripts/compute_consensus.py`
- `workflow/rules/qc.smk`
- `workflow/schemas/config.schema.yaml`
- Any manifest/output-contract/artifact files

## 10. Non-goals

- No SEACR IDR (key, rule, path, output)
- No broad IDR runtime (`chipseq_broad_experimental`, `cuttag_broad_experimental`)
- No MNase IDR
- No PE-only gating (SE supported)
- No `consensus.enabled` accidentally gating IDR
- No Artifact runtime adoption
- No Co-Authored-By

## 11. Verification

```bash
# New tests (35)
python3 test/test_stage64_cuttag_idr_config_validation.py
python3 test/test_stage64_cuttag_idr_dryrun.py
python3 test/test_stage64_cuttag_idr_summary.py

# ATAC IDR unchanged
python3 test/test_stage55_atac_idr_dryrun.py
python3 test/test_stage55_atac_idr_summary.py
python3 test/test_stage55_config_validation.py
python3 test/test_stage55_stage5_invariant.py

# Consensus unchanged
python3 test/test_stage62_consensus_dryrun.py
python3 test/test_stage63_seacr_consensus_dryrun.py

# Validation invariant
python3 test/test_stage53_config_validation.py
python3 test/test_stage53_reproducibility_policy_contract.py
python3 test/test_stage58_mixed_idr_validation.py

# Infrastructure
python3 test/test_stage57_shell_safety.py
python3 test/test_stage28_release_readiness.py
python3 test/test_no_hardcoded_paths.py
python3 test/test_stage54_consensus.py

# Git
git diff --check
```
