# Stage 23: Rule All Target Builder Cleanup — Design Spec

**Date:** 2026-05-24
**Status:** design review
**Scope:** Extract `rule all` input target generation into helper functions — no DAG behavior change
**Prerequisite:** Stage 20a + 20b + 22 (FE/ppois BigWig) completed and verified

---

## 1. Goals

1. Reduce `rule all` from ~550 inline lines to a small dispatch block (~10 lines) concatenating ~9 helper function return values.
2. Keep every existing target block semantically identical — same paths, same gating, same expand patterns, same named keyword arguments.
3. Make it possible to understand which targets belong together without reading a 550-line `rule all` block.
4. Do NOT introduce new abstractions, classes, or configuration schema.
5. Do NOT change any rule file, output path, or shell command.

---

## 2. Current State Assessment

The current `rule all` (lines 503–1013) contains ~50 target blocks organized chronologically by feature stage. The structure is flat — all targets are keyword arguments to a single `input` dict. This works correctly but is hard to navigate. The target blocks naturally cluster into semantic groups:

| Logical Group | Current Location | # of Blocks |
| :--- | :--- | :--- |
| Base pipeline (done, BigWig, peaks, multiqc) | Lines 505–522 | 4 |
| Single-sample QC (peak_counts through stage3_summary) | Lines 523–653 | 12 |
| Signal tracks (FE/ppois bdg + bw, single-sample) | Lines 573–606 | 4 |
| Blacklist (BAM + peaks) | Lines 607–648 | 4 |
| SEACR / CUT&Tag fragment | Lines 654–682 | 3 |
| Advanced QC (cross-corr, preseq, Picard) | Lines 683–762 | 7 |
| TSS enrichment | Lines 763–803 | 4 |
| Replicate pooled outputs | Lines 804–915 | 12 |
| IDR | Lines 916–1013 | 12 |

---

## 3. Target Group Design

### 3.1 Group assignment

| Helper Function | Targets | Rationale |
| :--- | :--- | :--- |
| `_base_targets()` | `done`, `bigwig`, `peaks`, `multiqc` | Core pipeline outputs always scheduled |
| `_blacklist_targets()` | `bl_bam`, `bl_bai`, `bl_peaks_narrow`, `bl_peaks_broad` | Single gating dimension (blacklist_filter + BLACKLIST_SAMPLES) |
| `_single_sample_qc_targets()` | `peak_counts`, `lib_complexity`, `nrf_pbc`, `frip`, `qc_summary`, `stage3_summary` | Core QC metrics, gated on summary/frip/nrf/library switches |
| `_signal_targets()` | `signal_fe`, `signal_ppois`, `signal_fe_bw`, `signal_ppois_bw`, `pooled_signal_fe`, `pooled_signal_ppois`, `pooled_signal_fe_bw`, `pooled_signal_ppois_bw` | All signal tracks — bedGraph + BigWig, single-sample + pooled. One group for all signal because they share the `signal_tracks` + `chrom_sizes` gating axis |
| `_cuttag_targets()` | `seacr_bg`, `seacr_peaks`, `cuttag_frag_size` | All CUT&Tag-specific outputs |
| `_advanced_qc_targets()` | `cross_corr_qc`, `cross_corr_pdf`, `cross_corr_summary`, `preseq_out`, `picard_align`, `picard_insert`, `picard_qual`, `picard_ins_hist` | All opt-in QC modules (cross-corr, preseq, Picard) |
| `_tss_targets()` | `tss_bed`, `tss_matrix`, `tss_profile_tsv`, `tss_profile_pdf` | TSS enrichment — single gating axis (tss_enrichment + GTF) |
| `_replicate_targets()` | `pooled_bam`, `pooled_bai`, `biorep_bam`, `biorep_bai`, `pooled_ctrl_bam`, `pooled_ctrl_bai`, `pooled_peaks`, `pooled_qc_summary` | Replicate-aware pooled + QC. Signal tracks separated to `_signal_targets()` |
| `_idr_targets()` | All 12 IDR target blocks | Single gating axis (`STAGE5` + `IDR_EXPERIMENTS`) |

### 3.2 Stage 22 BigWig targets placement

**Decision: `_signal_targets()`**

Rationale:
- Signal tracks (bedGraph + BigWig, single + pooled) share the same gating axis: `qc.signal_tracks` + `chrom_sizes`.
- Putting BigWig targets next to their bedGraph counterparts makes the signal contract self-documenting — a reader sees that `signal_tracks: true` produces bedGraph always, and BigWig additionally when `chrom_sizes` is configured.
- `_replicate_targets()` is about pooled BAMs, peaks, and QC — signal belongs elsewhere.
- `_signal_targets()` includes pooled BigWig targets because they share the same gating helpers (`SIGNAL_BW_EXPERIMENTS`).

### 3.3 Helper function signature

All helpers share the same contract:

```python
def _base_targets():
    """Return a list of target lists for rule all input."""
    return [
        expand(...),  # done
        expand(...),  # bigwig
        expand(...),  # peaks
        [f"{OUTDIR}/multiqc/multiqc_report.html"] if MULTIQC else [],
    ]
```

Each helper:
- Takes no arguments (captures top-level variables by closure — same as current inline code).
- Returns a flat list of lists.
- Each inner element is either an `expand()` result (which is itself a list) or a single-element list (like `[path] if condition else []`).
- May define local helper lists internally (e.g., `narrow_bl_samples`, `broad_bl_samples`) — no new global state.

### 3.4 rule all after refactoring

```python
rule all:
    input:
        sum([
            *_base_targets(),
            *_blacklist_targets(),
            *_single_sample_qc_targets(),
            *_signal_targets(),
            *_cuttag_targets(),
            *_advanced_qc_targets(),
            *_tss_targets(),
            *_replicate_targets(),
            *_idr_targets(),
        ], [])
```

`sum(list_of_lists, [])` flattens the nested lists into a single input list. This preserves the exact same input structure (a flat list of file paths and expand results) that Snakemake's DAG resolver expects.

Alternative considered: `list(itertools.chain(...))` — but that adds an import. `sum(..., [])` uses only builtins and mirrors the flattening pattern used in other Snakemake workflows.

---

## 4. Behavior Preservation Proof

The refactoring is a pure code movement — each target block is moved unchanged from an inline `rule all` keyword argument to a helper function's return list. The proof:

1. **Same expand() calls.** Every `expand(pattern, outdir=..., sample=...)` is copied verbatim.
2. **Same gating conditions.** Every `if QC_CONFIG.get(...) and ... else []` is preserved.
3. **Same variable capture.** Helpers access `OUTDIR`, `TREATMENT_SAMPLE_IDS`, `QC_CONFIG`, `BLACKLIST_SAMPLES`, `SIGNAL_BW_SAMPLE_IDS`, `SIGNAL_BW_EXPERIMENTS`, etc. via closure — same as the current inline code.
4. **Same target ordering.** `sum([...], [])` preserves insertion order (Python 3.7+ guarantees dict order; list concatenation preserves order). The relative order of target groups is: base → blacklist → single-sample QC → signal → CUT&Tag → advanced QC → TSS → replicates → IDR.
5. **Same flattening.** The current `rule all` input is a dict of named keyword arguments. After refactoring, it becomes a flat list. Snakemake's DAG resolution only cares about the actual file targets — the keyword argument names are labels within Snakemake's internal representation but do not affect DAG construction. Verified by: the existing smoke profiles check rule scheduling by name (e.g., `"signal_track_fe_bw" in output`), not by keyword argument names.

**Verification method:** Dry-run DAG is identical before and after. The smoke profiles' scheduling checks (which grep for rule names in `snakemake -n` output) still pass identically.

---

## 5. What Is NOT Changed

- No output paths changed.
- No rule shell commands changed.
- No rule files moved or renamed.
- No sample/config schema changed.
- No Stage 22 BigWig gating semantics changed.
- No new Python classes, decorators, or abstractions.
- No changes to `workflow/rules/*.smk`.
- No changes to `scripts/`.
- No changes to `test/` (except possibly adding a Stage 23 stress test — optional in this slice).

---

## 6. Answer to Design Questions

### Q1: 当前 `rule all` 分成哪些 target groups 最合理？
Nine groups (see Section 3.1): base, blacklist, single-sample QC, signal, CUT&Tag, advanced QC, TSS, replicates, IDR.

### Q2: 每个 target group 包含哪些现有 target blocks？
See Section 3.1 table.

### Q3: 如何证明 refactor 后 DAG 行为不变？
- Pure code movement, zero semantic changes.
- All 7 smoke profiles (with scheduling checks) pass identically.
- All 6 Stage 22 BigWig stress tests pass identically.
- All existing stage-specific stress tests pass identically.
- Dry-run output grep for all rule names confirms same set of scheduled rules.

### Q4: 哪些东西本阶段绝对不改？
See Section 5. The mandatory invariant: `snakemake -n` output before and after must produce identical rule scheduling for all 7 smoke profiles.

### Q5: Stage 22 新增的 FE/ppois BigWig targets 应放入哪个 group？
`_signal_targets()`. See Section 3.2 rationale.

### Q6: 是否需要移动 rule 文件？
No. `workflow/rules/*.smk` files stay where they are. This refactoring only touches the target-generation logic in `workflow/Snakefile`.

---

## 7. Out of Scope

- Moving or renaming rule files
- Changing output paths or rule semantics
- Introducing `dataclass` or other Python structural abstractions
- QC summary Python-ification (Stage 24)
- Manifest script (Stage 24a)
- Target builder beyond `rule all` (e.g., `checkpoint` integration)

---

## 8. Verification Plan

```bash
# 1. Validation
python3 scripts/validate_samples.py --config config/config.yaml
python3 test/test_validation_stress.py

# 2. Smoke profiles (7/7) — THE CRITICAL BEHAVIOR PRESERVATION CHECK
python3 test/test_stage8_smoke_profiles.py

# 3. Stage 22 BigWig stress tests (6/6)
python3 test/test_stage22_bigwig_stress.py

# 4. No-hardcoded-paths
python3 test/test_no_hardcoded_paths.py

# 5. Stage-specific stress tests (if available)
python3 test/test_stage4b_stress.py
python3 test/test_stage5a_stress.py
python3 test/test_stage5b_stress.py
python3 test/test_stage6a_stress.py
python3 test/test_stage6b_stress.py
python3 test/test_stage7a_stress.py
python3 test/test_stage7b_stress.py
python3 test/test_stage12_stress.py
python3 test/test_stage18_19_stress.py
```

---

## 9. Self-Review

- No TBDs or placeholders.
- All 6 design questions answered.
- Group assignment rationale is explicit.
- Behavior preservation proof is mechanical (pure code movement).
- Section 5 explicitly lists everything NOT changed.
- Verification covers all existing test suites.
