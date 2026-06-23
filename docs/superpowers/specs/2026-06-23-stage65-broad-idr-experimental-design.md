# Stage 65: Broad IDR Experimental Runtime — Design Spec

**Date:** 2026-06-23
**Status:** Implemented 2026-06-23
**Author:** YangZiChen-glitch (Kaslana)
**Depends on:** Stage 53/55/62/63/64 (narrow IDR, consensus, SEACR consensus)
**Blocks:** None

## 1. Purpose

Implement experimental opt-in broad-peak IDR for ChIP-seq broad and CUT&Tag
broad. The config keys `chipseq_broad_experimental` and
`cuttag_broad_experimental` already exist in the schema and validation.
Stage 65 adds the DAG rules and uses warning text that says broad IDR is
experimental while consensus remains available as fallback/report.

Broad IDR is **experimental opt-in only.** Production-supported IDR remains:
ChIP-seq narrow (legacy stage5) and ATAC narrow. Supported opt-in: CUT&Tag
narrow. Broad IDR is never gated by `stage5`, never default-enabled.

## 2. IDR broadPeak Output Contract

### 2.1 Verified behavior

The `idr` command (bioconda `idr >=2.0,<3`) supports `--input-file-type
broadPeak`. Per the upstream source (`nboley/idr`):

> "The output format mimics the input file type, with some additional fields."

When invoked with `broadPeak` input and `--idr-threshold`:

| Property | Value |
|----------|-------|
| Command | `idr --samples ... --input-file-type broadPeak --rank p.value --idr-threshold <T> --output-file <out>` |
| Output columns | 17 (for 2 replicates) |
| First 9 columns | Standard broadPeak: chrom, start, end, name, score, strand, signalValue, pValue, qValue |
| Columns 10–17 | localIDR, globalIDR, rep1_start, rep1_end, rep1_signalValue, rep2_start, rep2_end, rep2_signalValue |
| Standard broadPeak? | First 9 columns = valid broadPeak; full file is extended (same as narrowPeak IDR output) |
| Final suffix | `.broadPeak` |

### 2.2 Contract test

A dedicated test creates a synthetic 17-column broadPeak-like thresholded file
and asserts:
- The first 9 columns follow broadPeak format (chrom, start, end, name, score, strand, signalValue, pValue, qValue)
- `count_peaks()` ignores header/track lines and counts correctly
- The final copied output preserves all 17 columns
- Final suffix is `.broadPeak`

## 3. Architecture

**One shared file:** `workflow/rules/idr_broad.smk`

| Element | Convention |
|---------|-----------|
| Helper prefix | `_broad_idr_` |
| Rule prefix | `broad_idr_` |
| Modes | `chipseq_broad_experimental`, `cuttag_broad_experimental` |
| Assay dispatch | MACS3 args via `get_macs3_args()` per assay |

### 3.1 MACS3 args for broad IDR-ready calls

`_broad_idr_macs3_args(wildcards)` dispatches by assay:

| Assay | Args |
|-------|------|
| chipseq | `-f {BAMPE/BAM} -g {genome} -p {pvalue} --broad --broad-cutoff {cutoff} {extra}` |
| cuttag | `-f {BAMPE/BAM} -g {genome} -p {pvalue} --broad --broad-cutoff {cutoff} {extra}` |

Both use:
- Layout-aware format (`BAMPE` for PE, `BAM` for SE)
- Relaxed `-p` from `idr_macs3.pvalue` (default 0.1)
- `--broad --broad-cutoff` from `tool_parameters.macs3.broad_cutoff` (default 0.1)
- `idr_macs3.extra_args`
- **No Tn5 shift** for either assay — broad mode in existing `get_macs3_args_cuttag()` does not apply `--nomodel --shift --extsize` for broad peaks

Verification: `workflow/rules/cuttag.smk:get_macs3_args_cuttag()` emits Tn5 shift only for narrow; broad
uses only `--broad --broad-cutoff`. `chipseq.smk:get_macs3_args_chipseq()` uses standard model for narrow
and `--broad --broad-cutoff` for broad. Both match this design.

### 3.2 Rules (7, mirroring ATAC/CUT&Tag narrow pattern)

| # | Rule | Assay dispatch? | Conda |
|---|------|----------------|-------|
| 1 | `broad_idr_macs3_biorep` | Yes (MACS3 args) | macs3.yml |
| 2 | `broad_idr_true_replicates` | No (`--input-file-type broadPeak`) | idr.yml |
| 3 | `broad_idr_split_pseudoreps` | No | samtools.yml |
| 4 | `broad_idr_macs3_pseudorep` | Yes (MACS3 args) | macs3.yml |
| 5 | `broad_idr_self_pseudoreps` | No | idr.yml |
| 6 | `broad_idr_pooled_pseudoreps` | No | idr.yml |
| 7a | `broad_idr_chipseq_summary` | ChIP-seq broad summary + final peak | python.yml |
| 7b | `broad_idr_cuttag_summary` | CUT&Tag broad summary + final peak | python.yml |

## 4. Config and Validation

### 4.1 Config gates (Snakefile)

```python
BROAD_CHIPSEQ_IDR_ENABLED = (
    REPRO_ENABLED
    and REPRODUCIBILITY_CONFIG.get("idr", {}).get(
        "chipseq_broad_experimental", False)
)
BROAD_CUTTAG_IDR_ENABLED = (
    REPRO_ENABLED
    and REPRODUCIBILITY_CONFIG.get("idr", {}).get(
        "cuttag_broad_experimental", False)
)
BROAD_IDR_ENABLED = BROAD_CHIPSEQ_IDR_ENABLED or BROAD_CUTTAG_IDR_ENABLED
```

### 4.2 Config schema

Existing keys at `config.schema.yaml:328-339` — no change needed. Keys:
`chipseq_broad_experimental` (bool, default false),
`cuttag_broad_experimental` (bool, default false).

### 4.3 Validation warning (updated)

When either flag is true, the validation warning text is:

> "Experimental broad IDR enabled: reproducibility.idr.{flag}=true. Broad
> IDR outputs are experimental. Consensus remains available as
> fallback/report."

### 4.4 Eligibility validation

Added to `validate_replicate_groups()` with new parameters:
`reproducibility_idr_chipseq_broad` and `reproducibility_idr_cuttag_broad`.

| Flag | Assay | Peak mode | Bioreps | Layout |
|------|-------|-----------|---------|--------|
| `chipseq_broad_experimental` | chipseq | broad | exactly 2 | PE or SE |
| `cuttag_broad_experimental` | cuttag | broad | exactly 2 | PE or SE |

- Skip non-matching assays silently
- Skip non-broad peak_mode silently
- Validate exactly 2 biological replicates
- Fail if no eligible experiment exists
- Both require stage4b=true

### 4.5 IDR settings gate

Extended: `stage5 or atac_idr_enabled or cuttag_idr_enabled or broad_chipseq_idr_enabled or broad_cuttag_idr_enabled`

## 5. Output Namespace

```
06_reproducibility/idr/
  idr_peaks/{exp}_broad_{assay}_biorep{br}_idr.broadPeak
  idr_peaks/{exp}_broad_{assay}_{source}_pr{pr}_idr.broadPeak
  true_replicates/{exp}_broad_{assay}_idr.txt
  true_replicates/{exp}_broad_{assay}_idr.thresholded.broadPeak
  self_pseudoreplicates/{exp}_broad_{assay}_biorep{br}_idr.txt
  self_pseudoreplicates/{exp}_broad_{assay}_biorep{br}_idr.thresholded.broadPeak
  pooled_pseudoreplicates/{exp}_broad_{assay}_idr.txt
  pooled_pseudoreplicates/{exp}_broad_{assay}_idr.thresholded.broadPeak

05_pseudorep/
  {exp}_broad_{assay}_{source}.pr{1,2}.bam{,.bai}

06_reproducibility/final/
  {exp}.{assay}.macs3.broad.replicate_validated.idr.broadPeak
  reproducibility_summary.tsv
```

`{assay}` in paths = `chipseq` or `cuttag`, matching the config key prefix.

## 6. Final Semantics

| Mode | Condition | Final | Method |
|------|-----------|-------|--------|
| chipseq broad | broad_chipseq enabled + eligible | IDR final, consensus secondary | idr |
| chipseq broad | not enabled | Consensus final | consensus |
| cuttag broad | broad_cuttag enabled + eligible | IDR final, consensus secondary | idr |
| cuttag broad | not enabled | Consensus final | consensus |

### 6.1 Consensus interaction

Three updates to `consensus.smk`:

1. `_consensus_final_method()`: return `"idr"` for broad experiments in broad IDR lists
2. `_consensus_final_output()`: return IDR final path for broad experiments in broad IDR lists
3. `_consensus_targets()`: exclude experiments in `BROAD_CHIPSEQ_IDR_EXPERIMENTS` or `BROAD_CUTTAG_IDR_EXPERIMENTS` from broad consensus final targets

Consensus peak + summary remain under `consensus/` always.

## 7. Summary Script

**New:** `scripts/broad_idr_summary.py`

Clones the 15-column schema from `atac_idr_summary.py` / `cuttag_idr_summary.py`:

```
experiment, assay, peak_mode, caller, bio_rep_a, bio_rep_b,
true_peaks_Nt, pooled_peaks_Np, self1_peaks_N1, self2_peaks_N2,
rescue_ratio, self_consistency_ratio, reproducibility_status,
final_method, final_output
```

Logic: count non-header non-track lines, compute rescue/self ratios, pass if both < 2 and finite, copy true-rep to final. CLI text says "broad-peak IDR."

### 7.1 Summary path ambiguity fix

The broad summary rules write to `final/reproducibility_summary.tsv` — the same
generic path used by ATAC, CUT&Tag narrow, and legacy IDR summaries. Snakemake
also requires every output in one rule to contain the same wildcard set, so
Stage 65 uses two assay-specific summary rules inside the shared
`idr_broad.smk` file:

```python
rule broad_idr_chipseq_summary:
    wildcard_constraints:
        experiment = wildcard_choices(BROAD_CHIPSEQ_IDR_EXPERIMENTS)

rule broad_idr_cuttag_summary:
    wildcard_constraints:
        experiment = wildcard_choices(BROAD_CUTTAG_IDR_EXPERIMENTS)
```

This ensures the generic summary path is only claimed by one rule per
experiment while keeping final peak paths assay-specific.

## 8. Metadata

Two experiment lists in `metadata.smk`, gated on the respective config flags:

```
BROAD_CHIPSEQ_IDR_EXPERIMENTS, BROAD_CHIPSEQ_IDR_BIOREP_EXP_LIST, ...
BROAD_CUTTAG_IDR_EXPERIMENTS, BROAD_CUTTAG_IDR_BIOREP_EXP_LIST, ...

BROAD_IDR_EXPERIMENTS = BROAD_CHIPSEQ_IDR_EXPERIMENTS + BROAD_CUTTAG_IDR_EXPERIMENTS
BROAD_IDR_EXPERIMENT_ASSAY = assay aligned with BROAD_IDR_EXPERIMENTS
```

Plus expanded expansion lists for pseudorep splits, self-IDR, and pooled-IDR.

## 9. Targets

`_broad_idr_targets()` in `targets.smk`, gated on `BROAD_IDR_ENABLED`.
Expands all 6 target groups (biorep peaks, true-rep IDR, pseudorep peaks,
self-IDR, pooled-IDR, final peak + summary) for each eligible experiment.
Added to `rule all` via `+ _broad_idr_targets()`.

## 10. Test Plan

### New tests (~35)

| File | Tests | Focus |
|------|-------|-------|
| `test_stage65_broad_idr_config_validation.py` | 14 | Both flags, eligibility, stage4b, mixed sheets, IDR settings gate |
| `test_stage65_broad_idr_dryrun.py` | 13 | DAG rules, paths, consensus interaction, no legacy/SEACR IDR |
| `test_stage65_broad_idr_summary.py` | 8 | 17-column contract, counting, 15-column TSV, pass/fail, assay in TSV |

### Output contract test (in summary tests)

Synthetic 17-column broadPeak-like thresholded file:
- `chr1\t100\t200\tpeak_1\t1000\t.\t5.0\t3.0\t2.0\t0.01\t0.005\t100\t200\t5.0\t150\t250\t4.5`
- Assert `count_peaks()` returns 1
- Assert final output keeps all 17 columns
- Assert suffix is `.broadPeak`

### Regression (~246 tests unchanged)

| File | Tests |
|------|-------|
| All Stage 53 config/policy tests | unchanged |
| All Stage 55 ATAC IDR tests | unchanged |
| All Stage 62 consensus tests | unchanged |
| All Stage 63 SEACR consensus tests | unchanged |
| All Stage 64 CUT&Tag IDR tests | unchanged |
| `test_stage57_shell_safety.py` | all pass (new rules have `set -e -o pipefail`) |
| `test_stage28_release_readiness.py` | unchanged |
| `test_no_hardcoded_paths.py` | unchanged |

### Files not touched

- `workflow/rules/idr.smk`
- `workflow/rules/idr_atac.smk`
- `workflow/rules/idr_cuttag.smk`
- `scripts/stage5b_summary.py`
- `scripts/atac_idr_summary.py`
- `scripts/cuttag_idr_summary.py`
- `scripts/compute_consensus.py`
- `workflow/rules/qc.smk`

## 11. File Inventory

### Created

| File | Purpose |
|------|---------|
| `workflow/rules/idr_broad.smk` | 7 broad IDR rules + `_broad_idr_` helpers |
| `scripts/broad_idr_summary.py` | 15-column reproducibility QC summary |
| `test/test_stage65_broad_idr_config_validation.py` | 14 tests |
| `test/test_stage65_broad_idr_dryrun.py` | 13 tests |
| `test/test_stage65_broad_idr_summary.py` | 8 tests + output contract |

### Modified

| File | Change |
|------|--------|
| `workflow/Snakefile` | Two config gates, include, rule all, validator flags |
| `workflow/rules/metadata.smk` | Broad IDR experiment lists |
| `workflow/rules/targets.smk` | `_broad_idr_targets()` |
| `workflow/rules/consensus.smk` | Final method/output/target suppression for broad |
| `scripts/validate_samples.py` | Broad IDR eligibility + updated warning text + IDR settings gate |

## 12. Non-goals

- No SEACR IDR
- No MNase IDR
- No narrow IDR changes
- No legacy `06_idr/` changes
- No Artifact runtime
- No manifest/output-contract
- No Co-Authored-By

## 13. Verification

```bash
python3 test/test_stage65_broad_idr_config_validation.py
python3 test/test_stage65_broad_idr_dryrun.py
python3 test/test_stage65_broad_idr_summary.py
python3 test/test_stage53_*.py
python3 test/test_stage55_*.py
python3 test/test_stage62_consensus_dryrun.py
python3 test/test_stage63_seacr_consensus_dryrun.py
python3 test/test_stage64_cuttag_idr_dryrun.py
python3 test/test_stage57_shell_safety.py
python3 test/test_stage28_release_readiness.py
python3 test/test_no_hardcoded_paths.py
git diff --check
```
