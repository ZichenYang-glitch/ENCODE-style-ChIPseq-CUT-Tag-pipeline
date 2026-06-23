# Stage 62 Implementation Plan: Consensus DAG Integration

**Date:** 2026-06-23
**Status:** Ready for implementation
**Depends on:** Stage 54 (consensus engine), Stage 61 (policy alignment)
**Blocks:** Stage 63 (SEACR consensus)
**Design spec:** `docs/superpowers/specs/2026-06-23-stage62-consensus-dag-integration-design.md`

## Files

### To create

| File | Purpose |
|------|---------|
| `workflow/rules/consensus.smk` | Per-biorep MACS3 + consensus compute + consensus final rules |
| `test/test_stage62_consensus_dryrun.py` | DAG gating tests |

### To modify

| File | Change |
|------|--------|
| `workflow/Snakefile` | Reuse existing REPRO_ENABLED; add CONSENSUS_ENABLED, `include consensus.smk`, `+ _consensus_targets()` to rule all |
| `workflow/rules/metadata.smk` | CONSENSUS experiment lists + biorep expansion lists |
| `workflow/rules/targets.smk` | `_consensus_targets()` |

### NOT modified

- `workflow/rules/idr.smk`
- `workflow/rules/idr_atac.smk`
- `workflow/rules/peaks.smk`
- `workflow/rules/replicates.smk`
- `scripts/validate_samples.py`
- `config/config.yaml`
- Any manifest/output-contract/artifact files

## Implementation steps

### Step 1: Add config gating (Snakefile)

Add REPRO_ENABLED, CONSENSUS_ENABLED after existing Stage 55 config.

### Step 2: Add consensus experiment lists (metadata.smk)

Derive 5 mode-specific lists from MULTI_BIOREP_EXPERIMENTS, plus biorep
expansion lists for per-biorep peak targets.

### Step 3: Create consensus.smk

Three rules:
1. `consensus_macs3_biorep_peaks` — standard MACS3 on biorep BAMs
2. `consensus_compute` — compute_consensus.py on per-biorep peaks
3. `consensus_final` — copy to final/ for modes where consensus is primary

All helpers use `_consensus_` prefix.

### Step 4: Add target expansion (targets.smk)

`_consensus_targets()` — consensus peak/summary for all 5 modes, final for
broad/cuttag-narrow only (not chipseq narrow, not atac narrow).

### Step 5: Wire rule all (Snakefile)

`+ _consensus_targets()` in rule all input list.

### Step 6: Write tests

13 dry-run/DAG tests covering config gating, per-mode targets, final semantics,
and legacy invariance.

### Step 7: Verify

```bash
python3 test/test_stage62_consensus_dryrun.py
python3 test/test_stage54_consensus.py
python3 test/test_stage57_shell_safety.py
python3 test/test_stage55_config_validation.py
python3 test/test_stage55_stage5_invariant.py
python3 test/test_stage28_release_readiness.py
python3 test/test_no_hardcoded_paths.py
git diff --check
```

## Non-goals

- No SEACR consensus (Stage 63)
- No manifest/output-contract/artifact updates
- No new IDR rules
- No changes to idr.smk or idr_atac.smk
- No Co-Authored-By
