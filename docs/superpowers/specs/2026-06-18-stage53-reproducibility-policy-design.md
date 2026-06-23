# Stage 53: Reproducibility Policy and Output Semantics — Design Spec

**Date:** 2026-06-18
**Status:** Design approved; implementation pending
**Author:** YangZiChen-glitch (Kaslana) with IreneAdler (Claude Code / DeepSeek V4 Pro)

## 1. Overview

Stage 53 defines the reproducibility policy and output semantics for
replicate-validated peak outputs across all six peak-calling modes in the
pipeline. It is a **policy/config semantics stage only** — no consensus engine,
no new IDR rules, no DAG integration.

### 1.1 Motivation

The current pipeline has exactly one reproducibility path: Stage 5 IDR for
ChIP-seq narrow with exactly 2 biological replicates. Every other peak-calling
mode — ChIP-seq broad, CUT&Tag narrow/broad/SEACR, ATAC narrow — produces
per-sample peaks and pooled peaks but has no replicate-validated output of any
kind.

Pooled peaks (`04_peaks/pooled/`) are aggregate-signal peaks from merged BAMs,
not replicate-validated peaks. This must be codified as a hard semantic
contract.

### 1.2 Non-goals (whole Stage 53-58 scope)

1. No Artifact runtime adoption (Stage 51 paused).
2. No MNase caller implementation. MNase consensus is not scientifically
   well-defined for this scope.
3. No breaking existing output paths. `06_idr/` and `04_peaks/` remain
   untouched.
4. For broad modes, never enable IDR by default without explicit opt-in.
   SEACR IDR is not planned in this roadmap unless a future stage defines a
   justified rank scheme.
5. Never label pooled peaks as replicate-validated.
6. Never add 3+ biorep IDR. IDR is scoped to exactly 2 biological replicates.
   3+ replicate IDR requires a separate pairwise-IDR or replicate-selection
   policy and is out of scope for Stage 53-58.
7. Never replace `stage5` with `reproducibility`. They coexist.

## 2. Reproducibility Policy

### 2.1 Core Policy Statements

1. **Pooled peaks are not validated peaks.** `04_peaks/pooled/` outputs are
   MACS3 runs on merged BAMs — a statistical enhancement for signal detection,
   not replication-based validation. They MUST NOT be labeled or treated as
   replicate-validated outputs.

2. **Replicate-validated outputs require per-replicate peak calls.** Whether
   consensus or IDR, inputs come from individual peak calls on each biological
   replicate BAM.

3. **Consensus is the baseline method** for all modes. For a candidate peak
   region, count the number of distinct biological replicates with an
   overlapping peak. Keep the region if `support_count >= min_replicates`.
   This is N-of-M replicate support consensus — not strict majority vote.

4. **IDR is an advanced method** for narrowPeak modes. Rank-based
   irreproducible discovery rate analysis, currently scoped to exactly 2
   biological replicates. Well-established for ChIP-seq TF narrow peaks;
   requires separate validation for ATAC/CUT&Tag.

### 2.2 Per-mode Reproducibility Strategy

| # | Assay | Peak mode | Caller | Suffix | Primary reproducibility | Secondary / report | Notes |
|---|-------|-----------|--------|--------|------------------------|-------------------|-------|
| 1 | chipseq | narrow | MACS3 | `.narrowPeak` | IDR (legacy `stage5`) | Consensus | Production-supported IDR; legacy Stage 5 unchanged |
| 2 | chipseq | broad | MACS3 | `.broadPeak` | IDR when experimental flag enabled; consensus otherwise | Consensus fallback/report | Experimental opt-in only |
| 3 | cuttag | narrow | MACS3 | `.narrowPeak` | Consensus | IDR opt-in (supported) | IDR final when explicitly enabled |
| 4 | cuttag | broad | MACS3 | `.broadPeak` | IDR when experimental flag enabled; consensus otherwise | Consensus fallback/report | Experimental opt-in only |
| 5 | cuttag | — | SEACR | `.bed` | Consensus | No IDR planned | SEACR IDR out of scope; MNase IDR out of scope |
| 6 | atac | narrow | MACS3 | `.narrowPeak` | IDR (when `reproducibility.idr.atac_narrow: true`) | Consensus | Production-supported IDR; uses same narrowPeak machinery as ChIP-seq |

### 2.3 Consensus Strategy Details

- **Exactly 2 bioreps:** consensus means both replicates support the peak
  (intersection, `support_count >= min_replicates` with `min_replicates=2`).
- **3+ bioreps:** peaks supported by at least `min_replicates` distinct
  biological replicates are kept. Default `min_replicates: 2`.
- **Reciprocal overlap** threshold (default 0.5) determines whether a peak in
  one replicate "overlaps" a candidate region.
- **No requirement** for all pairwise intersections.

### 2.4 IDR Strategy Details

- **Scope:** Exactly 2 biological replicates. 3+ replicate IDR is out of scope
  for Stage 53-58.
- **Legacy:** `stage5: true` enables ChIP-seq narrow IDR with existing rules
  and paths. This behavior is untouchable.
- **ATAC narrow:** Uses the same narrowPeak IDR pipeline. Gated by
  `reproducibility.idr.atac_narrow: true`.
- **CUT&Tag narrow:** Opt-in via `reproducibility.idr.cuttag_narrow: true`.
  Not default. Fragment length and SNR differences from ChIP-seq warrant
  caution.
- **Broad:** Experimental only. Requires explicit
  `reproducibility.idr.chipseq_broad_experimental: true` or
  `reproducibility.idr.cuttag_broad_experimental: true`. IDR becomes final
  when explicitly enabled and eligible; consensus remains available as
  fallback/report.
- **SEACR IDR:** Not planned. BED output format and score scheme are not
  directly compatible with IDR's rank assumptions.

## 3. Config Design

### 3.1 New `reproducibility` Block

```yaml
# Stage 53+: Reproducibility policy for replicate-validated peak outputs.
# Added in Stage 53 — does not affect existing stage5 ChIP-seq narrow IDR.
reproducibility:
  # Master switch. When false, no new reproducibility outputs are generated.
  # Legacy stage5 IDR is orthogonal to this switch.
  enabled: false

  # ------------------------------------------------------------------
  # Consensus configuration (Stage 54+ implementation)
  # ------------------------------------------------------------------
  consensus:
    enabled: true
    min_replicates: 2           # Keep peaks with >= N supporting replicates
    reciprocal_overlap: 0.5     # Applied to pairwise overlap checks

  # ------------------------------------------------------------------
  # IDR configuration for non-legacy modes (Stage 55+ implementation)
  # ------------------------------------------------------------------
  idr:
    # null/omitted → inferred from stage5
    chipseq_narrow: null

    # Established narrowPeak modes — off by default
    atac_narrow: false
    cuttag_narrow: false

    # Experimental IDR — must be explicitly opted in.
    # IDR becomes final only when explicitly enabled and eligible.
    chipseq_broad_experimental: false
    cuttag_broad_experimental: false
```

### 3.2 `stage5` and `reproducibility` Interaction

| `stage5` | `reproducibility.enabled` | `reproducibility.idr.chipseq_narrow` | Legacy Stage 5 IDR | New chipseq_narrow flag | Behavior |
|----------|---------------------------|--------------------------------------|--------------------|------------------------|----------|
| false | false | null (omitted) | not run | inferred false | No IDR, no reproducibility |
| true | false | null (omitted) | **runs** | inferred true | Legacy IDR only |
| true | true | null (omitted) | **runs** | inferred true | Legacy IDR + new reproducibility |
| false | true | null (omitted) | not run | inferred false | New reproducibility only, no legacy IDR |
| true | true | true | **runs** | true (explicit) | Both, consistent |
| true | true | false | **runs** | false (explicit) | **WARNING:** contradiction. Legacy IDR still runs. |

**Key invariant:** No setting under `reproducibility` may break, disable,
rename, or move existing Stage 5 outputs. `stage5: true` always keeps existing
ChIP-seq narrow IDR behavior.

### 3.3 Validation Rules (in `validate_config`)

When `reproducibility.enabled: true`:

1. `consensus.min_replicates` must be int ≥ 2. Default 2.
2. `consensus.reciprocal_overlap` must be float in (0, 1]. Default 0.5.
3. `idr.chipseq_narrow` accepts true, false, or null. If explicitly set to
   false while `stage5: true`, emit a config contradiction warning but legacy
   Stage 5 still runs.
4. `idr.atac_narrow`, `idr.cuttag_narrow` must be bool. Default false.
5. `idr.chipseq_broad_experimental`, `idr.cuttag_broad_experimental` must be
   bool. Default false. If true, emit an informational warning summarizing
   experimental status.
6. No validation error for experiments with <2 bioreps. Eligibility is
   determined at target expansion time (Stage 55+), not at config validation.

When `reproducibility.enabled: false`, all sub-keys are ignored for DAG
purposes with no validation errors.

## 4. Output Paths

### 4.1 New Namespace

```
results/experiments/<exp>/06_reproducibility/
  consensus/
    <exp>.<assay>.<caller>.<peak_mode>.consensus.<suffix>
    <exp>.<assay>.<caller>.<peak_mode>.consensus.summary.tsv
  idr/
    ... future non-legacy IDR outputs ...
  final/
    <exp>.<assay>.<caller>.<peak_mode>.replicate_validated.<method>.<suffix>
    reproducibility_summary.tsv
```

### 4.2 Legacy Paths (Unchanged)

```
results/experiments/<exp>/06_idr/
  true_replicates/idr.txt, idr.thresholded.narrowPeak
  self_pseudoreplicates/biorep{1,2}.idr.txt, .thresholded.narrowPeak
  pooled_pseudoreplicates/idr.txt, idr.thresholded.narrowPeak
  final/conservative.narrowPeak, optimal.narrowPeak, reproducibility_summary.tsv
```

### 4.3 Consensus Filename Templates

| Mode | Consensus peak | Summary |
|------|---------------|---------|
| chipseq narrow MACS3 | `<exp>.chipseq.macs3.narrow.consensus.narrowPeak` | `<exp>.chipseq.macs3.narrow.consensus.summary.tsv` |
| chipseq broad MACS3 | `<exp>.chipseq.macs3.broad.consensus.broadPeak` | `<exp>.chipseq.macs3.broad.consensus.summary.tsv` |
| cuttag narrow MACS3 | `<exp>.cuttag.macs3.narrow.consensus.narrowPeak` | `<exp>.cuttag.macs3.narrow.consensus.summary.tsv` |
| cuttag broad MACS3 | `<exp>.cuttag.macs3.broad.consensus.broadPeak` | `<exp>.cuttag.macs3.broad.consensus.summary.tsv` |
| cuttag SEACR | `<exp>.cuttag.seacr.<seacr_mode>.consensus.bed` | `<exp>.cuttag.seacr.<seacr_mode>.consensus.summary.tsv` |
| atac narrow MACS3 | `<exp>.atac.macs3.narrow.consensus.narrowPeak` | `<exp>.atac.macs3.narrow.consensus.summary.tsv` |

### 4.4 Final Validated Filename Templates

| Mode | Primary method | Final output |
|------|---------------|-------------|
| chipseq narrow | IDR (legacy `stage5`) | Use legacy `06_idr/final/conservative.narrowPeak` |
| chipseq broad | Consensus | `<exp>.chipseq.macs3.broad.replicate_validated.consensus.broadPeak` |
| cuttag narrow | IDR (when `reproducibility.idr.cuttag_narrow: true`) | `<exp>.cuttag.macs3.narrow.replicate_validated.idr.narrowPeak` |
| cuttag narrow | Consensus (when CUT&Tag IDR not enabled) | `<exp>.cuttag.macs3.narrow.replicate_validated.consensus.narrowPeak` |
| cuttag broad | Consensus | `<exp>.cuttag.macs3.broad.replicate_validated.consensus.broadPeak` |
| cuttag SEACR | Consensus | `<exp>.cuttag.seacr.<seacr_mode>.replicate_validated.consensus.bed` |
| atac narrow | IDR (when `reproducibility.idr.atac_narrow: true`) | `<exp>.atac.macs3.narrow.replicate_validated.idr.narrowPeak` |
| atac narrow | Consensus (when `reproducibility.idr.atac_narrow: false`) | Consensus available under `consensus/`; no `final/` output until IDR is enabled |

For chipseq narrow: do not symlink or copy legacy IDR outputs into
`06_reproducibility/final/`. The `reproducibility_summary.tsv` records the
primary final source as `../06_idr/final/conservative.narrowPeak`.

### 4.5 Consensus Summary TSV Schema

Each consensus run produces a summary TSV with the following columns:

```
experiment   assay   peak_mode   caller   n_bioreps   min_replicates
reciprocal_overlap   consensus_peak_count   support_distribution
biorep_labels   source_peak_files   final_method   final_output
```

- `support_distribution`: JSON mapping, e.g. `{"2": 15000, "3": 8000}` —
  how many consensus peaks are supported by N replicates.
- `biorep_labels`: comma-separated list of biorep labels participating.
- `source_peak_files`: paths to per-biorep peak files used.
- `final_method`: `consensus` or `idr` — which method is primary final.
- `final_output`: path to the primary final peak file.

### 4.6 Consensus Output Format Caveats

For narrowPeak/broadPeak consensus outputs:
- chrom/start/end/name represent the merged consensus interval.
- score/name are derived deterministically.
- signalValue/pValue/qValue columns use a conservative deterministic
  aggregation: max signal among supporting peaks, or documented placeholder
  values.
- The summary TSV is the authoritative source for support_count,
  supporting_bioreps, source_peak_count, and method parameters.

## 5. Verification and Test Plan

### 5.1 Stage 53 Tests (Policy/Config Only)

| # | Test file | What it tests | Implementation requirement |
|---|-----------|--------------|---------------------------|
| 1 | `test_stage53_config_validation.py` | `reproducibility` block parsing: defaults, bad values, null inference | Requires `validate_config` additions |
| 2 | `test_stage53_stage5_invariant.py` | `stage5: true` not suppressed by any `reproducibility` setting | Requires `validate_config` additions |
| 3 | `test_stage53_pooled_not_validated.py` | Policy doc + output-contract entries classify pooled peaks correctly | Documentation-only: reads `docs/reproducibility-policy.md` |
| 4 | `test_stage53_experimental_warnings.py` | Experimental IDR flags emit informational warnings | Requires `validate_config` additions |
| 5 | `test_stage53_reproducibility_policy_contract.py` | Policy matrix includes all 6 modes | Documentation-only: reads `docs/reproducibility-policy.md` |
| 6 | `test_stage53_output_path_templates.py` | Path templates follow `<exp>.<assay>.<caller>.<peak_mode>.<type>.<suffix>` | Documentation-only: reads `docs/reproducibility-policy.md` |

### 5.2 Tests Deferred to Later Stages

| Stage | Test | Reason for deferral |
|-------|------|---------------------|
| 54 | Consensus summary TSV schema | Requires consensus engine to produce real output |
| 54 | Consensus engine unit tests | Requires `scripts/compute_consensus.py` |
| 54 | N-of-M support distribution correctness | Requires consensus engine |
| 55+ | Dry-run DAG tests for reproducibility nodes | Requires DAG integration |
| 55+ | ATAC IDR target existence in dry-run | Requires ATAC IDR rules |
| 56+ | CUT&Tag IDR gating | Requires CUT&Tag IDR policy layer |
| 58 | Manifest entries for new outputs | Requires final output paths |

## 6. Stage Boundaries

| Stage | Scope | Key deliverables |
|-------|-------|-----------------|
| **53** | Policy + config semantics | `docs/reproducibility-policy.md`; optional `validate_config` additions; doc/template-only tests |
| **54** | Consensus engine | `scripts/compute_consensus.py`; summary writer; unit tests across all 6 modes; NO DAG integration |
| **55** | ATAC narrow IDR | Generalize `idr.smk` factories for ATAC; same narrowPeak pipeline with ATAC gating; DAG integration only for enabled modes |
| **56** | CUT&Tag narrow IDR (opt-in) | Policy layer on top of ATAC IDR path; experimental warnings; off by default |
| **57** | Broad experimental IDR policy | Validation warnings for `chipseq_broad_experimental` and `cuttag_broad_experimental`; IDR becomes final only when explicitly enabled and eligible; consensus remains fallback/report; SEACR remains consensus-only |
| **58** | Manifest / contract / report integration + release hardening | Manifest awareness for new outputs; output contract updates; MultiQC integration if meaningful; regression tests since Stage 53 |

### Stage 53 Deliverables

1. `docs/reproducibility-policy.md` — comprehensive policy document
2. Optional `validate_config` additions for the `reproducibility` block
3. 6 tests (policy/config only, see §5.1)

### Deferred to Stage 54+

- Consensus engine implementation and tests (Stage 54)
- DAG target expansion (Stage 55+)
- Snakemake rules for new IDR modes (Stage 55+)
- Final file writing (Stage 55+)
- Manifest updates (Stage 58)
- MultiQC integration (Stage 58)

## 7. Design Decisions Record

1. **N-of-M, not strict majority.** Consensus uses configurable `min_replicates`
   threshold, not a hard-coded majority rule. This matches the config intent
   and extends naturally from 2 replicates to N replicates.

2. **Consensus always, IDR where enabled.** Consensus is generated for every
   mode with ≥2 biological replicates. IDR is the primary final method where
   domain-validated (ChIP-seq narrow, ATAC narrow), supported opt-in
   (CUT&Tag narrow), or explicitly enabled as experimental broad IDR.

3. **Experimental IDR requires explicit opt-in.** For broad modes, IDR becomes
   final only when the corresponding experimental flag is true and the
   experiment is eligible. Consensus remains available as fallback/report.

4. **SEACR IDR excluded.** SEACR BED format's score scheme is not directly
   compatible with IDR's rank assumptions. No experimental flag is exposed.

5. **`stage5` preserved as legacy.** The existing ChIP-seq narrow IDR switch
   is untouchable. `reproducibility` is an orthogonal new layer.

6. **No symlinks to legacy paths.** `06_reproducibility/final/` does not
   contain symlinks to `06_idr/final/`. The reproducibility summary documents
   the source instead.

7. **3+ biorep IDR deferred.** IDR with ≥3 replicates requires pairwise-IDR
   or replicate-selection policy design, which is out of scope for Stage 53-58.
