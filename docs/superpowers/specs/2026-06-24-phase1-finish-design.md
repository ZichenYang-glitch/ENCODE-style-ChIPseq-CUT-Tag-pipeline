# Phase 1 Finish: Consolidate IDR Metadata, Rules, and Summary Scripts

**Date:** 2026-06-24
**Status:** Design approved; implementation in progress
**Author:** Zichen Yang
**Depends on:** Phase 1 IDR consolidation (PR #39)
**Blocks:** Phase 2

## 1. Overview

Phase 1 (PR #39) extracted shared IDR path/args/label helpers into `workflow/rules/idr_paths.smk`,
consolidated the four IDR target builders in `workflow/rules/targets.smk`, and added a
behavioral contract test suite. Three categories of duplication remain:

1. `metadata.smk` contains four nearly identical blocks that build IDR experiment/pseudorep expansion lists.
2. Four IDR `.smk` rule files (`idr.smk`, `idr_atac.smk`, `idr_cuttag.smk`, `idr_broad.smk`) share the same rule structure.
3. Four IDR summary scripts (`stage5b_summary.py`, `atac_idr_summary.py`, `cuttag_idr_summary.py`, `broad_idr_summary.py`) share the same ratio-computation logic.

This design finishes Phase 1 by removing those duplicates while keeping the public behavior of the pipeline unchanged.

### 1.1 Motivation

- Reduce maintenance surface: a change to IDR pseudorep expansion logic currently requires editing four blocks in `metadata.smk`.
- Reduce rule-file sprawl: ATAC/CUT&Tag/broad IDR rules are structurally identical.
- Reduce script sprawl: summary ratio logic is copy-pasted across four scripts.
- Keep the Phase 1 guarantee: zero behavior change for rule outputs, wildcards, target lists, thresholds, and scientific defaults.

### 1.2 Scope

| In scope | Out of scope |
|----------|-------------|
| Parameterize IDR list builders in `metadata.smk` | Refactoring non-IDR metadata |
| Unify ATAC/CUT&Tag/broad IDR rules into one parameterized file | Merging legacy Stage 5 `idr.smk` into the unified file |
| Extract shared `_idr_filename_infix()` helper | Generalizing MACS3 callpeak shell logic outside IDR |
| Unify ATAC/CUT&Tag/broad summary scripts | Rewriting `stage5b_summary.py` internals or unifying it |
| Thin wrappers for backward-compatible modern script CLI | New summary TSV columns or thresholds |
| DAG snapshot preservation | Real execution tests |

### 1.3 Architecture Decision

#### Step A — `metadata.smk`: shared helper functions, not a dataclass

Rejected a formal `IdrProfile` dataclass because:

- The global variables (`IDR_*`, `ATAC_IDR_*`, `CUTTAG_IDR_*`, `BROAD_IDR_*`) are referenced directly by `targets.smk` and by rule files.
- Introducing a dataclass would require either translating attributes back into flat lists or changing every consumer.
- Shared helper functions preserve the existing flat-list contract with minimal churn.

Selected design:

- `_idr_exp_matches(exp, assay, peak_mode)` — predicate for experiment selection.
- `_build_idr_experiment_lists(enabled, assay, peak_mode, exps)` — returns `(experiments, biorep_exp_list, biorep_list)`.
- `_build_idr_pseudorep_lists(experiments)` — returns the seven narrow-ID pseudorep expansion lists.
- `_build_broad_idr_pseudorep_lists(experiments, assays)` — returns the thirteen broad-ID expansion lists, including the assay dimension.

#### Step B — IDR rule files: unify ATAC/CUT&Tag/broad into `idr_reproducibility.smk`, keep legacy `idr.smk` separate

Rejected a single file containing legacy Stage 5 as well because:

- Legacy Stage 5 uses the `06_idr/` namespace and different filename conventions (`{exp}_biorep{br}_idr_peaks.narrowPeak`, no assay infix).
- Forcing legacy into a `{assay}`/`{peak_mode}` wildcard model would require an extra "mode" wildcard and dynamic output dispatch, increasing complexity and regression risk.
- ATAC/CUT&Tag/broad already share the `06_reproducibility/idr/` namespace and differ only in filename infix and rule prefix.

Selected design:

- Create `workflow/rules/idr_reproducibility.smk` parameterized by `{assay}` and `{peak_mode}` wildcards.
- It handles `atac`/`cuttag` narrow and `chipseq`/`cuttag` broad.
- Legacy `workflow/rules/idr.smk` remains untouched.
- Add `_idr_filename_infix(assay, peak_suffix)` in `idr_paths.smk` to replace the duplicated four-way dispatch blocks.

#### Step C — summary scripts: unified `idr_reproducibility_summary.py` with thin wrappers

Rejected deleting the old scripts immediately because external callers (tests, documentation, old config examples) may still reference them.

Selected design:

- Create `scripts/idr_reproducibility_summary.py` that accepts `--assay` and `--peak-mode`.
- Keep the three modern scripts (`atac_idr_summary.py`, `cuttag_idr_summary.py`, `broad_idr_summary.py`) as thin wrappers that delegate to the unified script.
- Legacy `scripts/stage5b_summary.py` is intentionally **not** unified: its TSV schema (12 columns, `--output-cons`/`--output-opt`) is incompatible with the modern 15-column schema (`--final-output`/`--output-peak`). Unifying it would require rewriting its internals, which is out of scope.
- After verification, the modern wrappers can be removed in a follow-up cleanup stage.

## 2. Data Structures

### 2.1 `metadata.smk` globals (unchanged contract)

All flat list names and semantics are preserved:

- `IDR_EXPERIMENTS`, `IDR_BIOREP_EXP_LIST`, `IDR_BIOREP_LIST`
- `IDR_SPLIT_SOURCE_EXP`, `IDR_SPLIT_SOURCE_NAME`
- `IDR_PR_PEAK_EXP`, `IDR_PR_PEAK_SRC`, `IDR_PR_PEAK_PR`
- `IDR_SELF_EXP`, `IDR_SELF_BR`
- `ATAC_IDR_*` mirrors above
- `CUTTAG_IDR_*` mirrors above
- `BROAD_CHIPSEQ_IDR_*`, `BROAD_CUTTAG_IDR_*`, `BROAD_IDR_*`, `BROAD_IDR_EXPERIMENT_ASSAY`, `BROAD_IDR_*_ASSAY`

### 2.2 `idr_paths.smk` helper additions

```python
def _idr_filename_infix(assay, peak_suffix):
    """Return the filename infix for reproducibility IDR peak files.

    Valid combinations:
      - assay in {"atac", "cuttag"}, peak_suffix == "narrowPeak"
      - assay in {"chipseq", "cuttag"}, peak_suffix == "broadPeak"
    """
    if peak_suffix == "narrowPeak" and assay in ("atac", "cuttag"):
        return f"{assay}_"
    if peak_suffix == "broadPeak" and assay in ("chipseq", "cuttag"):
        return f"broad_{assay}_"
    raise ValueError(...)
```

This replaces identical dispatch blocks in `idr_repro_peak_input`, `idr_self_thresh_path`, `idr_pooled_peak_input`, and `idr_pooled_thresh_path`.

The peak-mode to file-suffix mapping is handled directly inside `idr_reproducibility.smk` by hardcoding the appropriate suffix (`narrowPeak` / `broadPeak`) in each rule. A separate helper is unnecessary because narrow and broad modes are split into distinct rules (see §3.1).

### 2.3 `idr_reproducibility.smk` wildcard constraints

```python
wildcard_constraints:
    assay = r"chipseq|atac|cuttag"
    peak_mode = r"narrow|broad"
```

Rule-level wildcard constraints will further restrict `{assay}` for broad rules to `chipseq|cuttag`.

## 3. Rules

### 3.1 Unified rules in `idr_reproducibility.smk`

Narrow and broad modes are split into separate rules. The raw IDR text output (`*.txt`) does not contain a peak-mode suffix, so mixing both modes in one rule would violate Snakemake's requirement that all output files share the same wildcards. Within each mode, the `{assay}` wildcard parameterizes ATAC vs. CUT\&Tag (narrow) or ChIP-seq vs. CUT\&Tag (broad).

| Rule | Mode | Wildcards | Outputs |
|------|------|-----------|---------|
| `idr_macs3_biorep_narrow` | narrow | `experiment`, `assay`, `bio_rep` | `06_reproducibility/idr/idr_peaks/{experiment}_{assay}_biorep{bio_rep}_idr.narrowPeak` |
| `idr_macs3_biorep_broad` | broad | `experiment`, `assay`, `bio_rep` | `06_reproducibility/idr/idr_peaks/{experiment}_broad_{assay}_biorep{bio_rep}_idr.broadPeak` |
| `idr_true_replicates_narrow` | narrow | `experiment`, `assay` | `true_replicates/{experiment}_{assay}_idr.txt`, `.thresholded.narrowPeak` |
| `idr_true_replicates_broad` | broad | `experiment`, `assay` | `true_replicates/{experiment}_broad_{assay}_idr.txt`, `.thresholded.broadPeak` |
| `idr_split_pseudoreps_narrow` | narrow | `experiment`, `assay`, `source` | `05_pseudorep/{experiment}_{assay}_{source}.pr{1,2}.bam` |
| `idr_split_pseudoreps_broad` | broad | `experiment`, `assay`, `source` | `05_pseudorep/{experiment}_broad_{assay}_{source}.pr{1,2}.bam` |
| `idr_macs3_pseudorep_narrow` | narrow | `experiment`, `assay`, `source`, `pr` | `idr_peaks/..._{source}_pr{pr}_idr.narrowPeak` |
| `idr_macs3_pseudorep_broad` | broad | `experiment`, `assay`, `source`, `pr` | `idr_peaks/..._{source}_pr{pr}_idr.broadPeak` |
| `idr_self_pseudoreps_narrow` | narrow | `experiment`, `assay`, `bio_rep` | `self_pseudoreplicates/...biorep{bio_rep}_idr.txt`, `.thresholded.narrowPeak` |
| `idr_self_pseudoreps_broad` | broad | `experiment`, `assay`, `bio_rep` | `self_pseudoreplicates/...biorep{bio_rep}_idr.txt`, `.thresholded.broadPeak` |
| `idr_pooled_pseudoreps_narrow` | narrow | `experiment`, `assay` | `pooled_pseudoreplicates/...idr.txt`, `.thresholded.narrowPeak` |
| `idr_pooled_pseudoreps_broad` | broad | `experiment`, `assay` | `pooled_pseudoreplicates/...idr.txt`, `.thresholded.broadPeak` |
| `idr_summary_atac_narrow` | narrow | `experiment` | `final/reproducibility_summary.tsv`, `final/{experiment}.atac.macs3.narrow.replicate_validated.idr.narrowPeak` |
| `idr_summary_cuttag_narrow` | narrow | `experiment` | `final/reproducibility_summary.tsv`, `final/{experiment}.cuttag.macs3.narrow.replicate_validated.idr.narrowPeak` |
| `idr_summary_chipseq_broad` | broad | `experiment` | `final/reproducibility_summary.tsv`, `final/{experiment}.chipseq.macs3.broad.replicate_validated.idr.broadPeak` |
| `idr_summary_cuttag_broad` | broad | `experiment` | `final/reproducibility_summary.tsv`, `final/{experiment}.cuttag.macs3.broad.replicate_validated.idr.broadPeak` |

Rule names change compared with the old assay-specific files (e.g. `atac_macs3_idr_biorep` becomes `idr_macs3_biorep_narrow`). The 8 smoke profiles do not enable modern IDR modes, so DAG snapshots are unaffected.

For broad-peak experiments, `idr_summary` is instantiated once per assay (chipseq, cuttag) because the final peak filename embeds the assay and the shared `reproducibility_summary.tsv` output cannot carry an `{assay}` wildcard.

### 3.2 Legacy rules in `idr.smk`

Unchanged.

## 4. Scripts

### 4.1 New unified script

`scripts/idr_reproducibility_summary.py` CLI:

```
--true-peaks
--pooled-peaks
--self1-peaks
--self2-peaks
--experiment
--assay
--caller
--peak-mode
--bio-rep-a
--bio-rep-b
--final-method
--final-output
--output-tsv
--output-peak
```

Behavior identical to existing scripts: compute rescue ratio, self-consistency ratio, reproducibility status, copy final peak, write 15-column TSV.

### 4.2 Thin wrappers

- `scripts/atac_idr_summary.py`
- `scripts/cuttag_idr_summary.py`
- `scripts/broad_idr_summary.py`

Each wrapper calls `idr_reproducibility_summary.py` with the appropriate `--assay`/`--peak-mode` defaults and forwards all other arguments.

`scripts/stage5b_summary.py` is **not** converted to a wrapper and remains unchanged.

## 5. Acceptance Criteria

- `metadata.smk` IDR list globals remain identical in content, order, and type.
- All rule output paths in `idr_atac.smk`, `idr_cuttag.smk`, `idr_broad.smk` are preserved.
- Rule names are preserved OR renamed in a backward-compatible way that does not break target resolution.
- `test/test_dag_snapshots.py` passes without snapshot updates.
- `test/test_idr_paths.py` passes.
- `test/test_targets.py` passes.
- Full pytest suite: 331 passed, 5 skipped.
- `git diff --check` clean.
- IDR thresholds, p-values, broad cutoffs, and summary TSV schema unchanged.

## 6. Risks / Rollback

| Risk | Mitigation |
|------|-----------|
| Snakemake wildcard ambiguity when merging broad/narrow into one file | Strict wildcard_constraints per rule; verify DAG snapshots |
| Different shell behavior between legacy and modern rules (e.g. MACS3 output rename) | Keep the modern rules' exact shell blocks when unifying |
| Summary script wrapper introduces subtle CLI parsing differences | Run existing summary tests against wrappers |
| DAG snapshot changes | Stop immediately if any snapshot diff is non-zero |

Rollback: revert the commits on `refactor/phase1-finish`. Legacy `idr.smk` is untouched, so Stage 5 behavior is protected throughout.
