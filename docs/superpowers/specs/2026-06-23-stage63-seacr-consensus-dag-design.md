# Stage 63: SEACR Consensus DAG Integration — Design Spec

**Date:** 2026-06-23
**Status:** Implemented
**Author:** YangZiChen-glitch (Kaslana)
**Depends on:** Stage 62 (MACS3 consensus DAG)
**Blocks:** None

## 1. Overview

Stage 63 adds SEACR consensus to the consensus DAG for CUT&Tag PE experiments.
SEACR IDR is NOT planned — consensus is the primary reproducibility method.

### 1.1 Gating

SEACR consensus runs only when ALL are true:
- `REPRO_ENABLED` / `reproducibility.enabled`
- `CONSENSUS_ENABLED` / `reproducibility.consensus.enabled`
- `SEACR_ENABLED` / `cuttag.seacr.enabled`
- `STAGE4B`
- all treatment samples in the experiment have assay == cuttag and layout == PE
- >= 2 biological replicates

### 1.2 Rules added to consensus.smk

| Rule | Purpose |
|------|---------|
| `consensus_seacr_bedgraph` | Per-biorep bedGraph from biorep BAM |
| `consensus_seacr_biorep_peaks` | Per-biorep SEACR BED from bedGraph |
| `consensus_compute_seacr` | Consensus BED + summary via compute_consensus.py |
| `consensus_final_seacr` | Copy consensus to final/ |

Eligible SEACR consensus experiments are discovered in `workflow/rules/metadata.smk`
and consumed by `_consensus_targets()` in `workflow/rules/targets.smk`.

### 1.3 Output paths

```
06_reproducibility/consensus/biorep_peaks/<exp>.biorep{br}.cuttag.seacr.bedgraph
06_reproducibility/consensus/biorep_peaks/<exp>.biorep{br}.cuttag.seacr.<mode>.bed
06_reproducibility/consensus/<exp>.cuttag.seacr.<mode>.consensus.bed
06_reproducibility/consensus/<exp>.cuttag.seacr.<mode>.consensus.summary.tsv
06_reproducibility/final/<exp>.cuttag.seacr.<mode>.replicate_validated.consensus.bed
```

## 2. Non-goals

- No SEACR IDR key, rule, path, or output
- No changes to qc.smk sample-level SEACR rules
- No changes to compute_consensus.py
- No manifest/output-contract/artifact changes
- No Co-Authored-By
