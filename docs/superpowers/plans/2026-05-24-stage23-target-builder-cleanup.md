# Stage 23: Rule All Target Builder Cleanup — Implementation Plan

> **For agentic workers:** Implement tasks in order. Each step includes verification.
> Full design: `docs/superpowers/specs/2026-05-24-stage23-target-builder-cleanup-design.md`

**Goal:** Extract `rule all` target generation logic into 9 helper functions. Pure code movement — zero DAG behavior change.

**Architecture:** Each helper returns a flat list of lists; `rule all` flattens via `sum([...], [])`. No new imports, classes, or abstractions.

---

### Task 1: Replace rule all with helpers

**Files:**
- Edit: `workflow/Snakefile` — replace lines 498–1013 (current `rule all` block)

**Strategy:** Write all 9 helpers above `rule all`, then replace `rule all` with the dispatch form. Each helper gets a comment header referencing the original line range.

- [ ] **Step 1: Define `_base_targets()`**

```python
def _base_targets():
    """Core pipeline outputs: done markers, CPM BigWig, peaks, MultiQC."""
    return [
        expand(
            "{outdir}/{sample}/logs/{sample}.pipeline.done",
            outdir=OUTDIR,
            sample=TREATMENT_SAMPLE_IDS,
        ),
        expand(
            "{outdir}/{sample}/03_bigwig/{sample}.CPM.bw",
            outdir=OUTDIR,
            sample=ACTIVE_SAMPLE_IDS,
        ),
        expand(
            "{outdir}/{sample}/04_peaks/{sample}",
            outdir=OUTDIR,
            sample=TREATMENT_SAMPLE_IDS,
        ),
        [f"{OUTDIR}/multiqc/multiqc_report.html"] if MULTIQC else [],
    ]
```

- [ ] **Step 2: Define `_blacklist_targets()`**

Compute `_bl_narrow` and `_bl_broad` as local variables inside the function.

- [ ] **Step 3: Define `_single_sample_qc_targets()`**

`peak_counts`, `lib_complexity`, `nrf_pbc`, `frip`, `qc_summary`, `stage3_summary`.

- [ ] **Step 4: Define `_signal_targets()`**

All single-sample bedGraph + BigWig + pooled bedGraph + pooled BigWig (8 targets).

- [ ] **Step 5: Define `_cuttag_targets()`**

`seacr_bg`, `seacr_peaks`, `cuttag_frag_size`.

- [ ] **Step 6: Define `_advanced_qc_targets()`**

`cross_corr_qc`, `cross_corr_pdf`, `cross_corr_summary`, `preseq_out`, `picard_align`, `picard_insert`, `picard_qual`, `picard_ins_hist`.

- [ ] **Step 7: Define `_tss_targets()`**

`tss_bed`, `tss_matrix`, `tss_profile_tsv`, `tss_profile_pdf`.

- [ ] **Step 8: Define `_replicate_targets()`**

`pooled_bam`, `pooled_bai`, `biorep_bam`, `biorep_bai`, `pooled_ctrl_bam`, `pooled_ctrl_bai`, `pooled_peaks`, `pooled_qc_summary`.

- [ ] **Step 9: Define `_idr_targets()`**

All 12 IDR target blocks.

- [ ] **Step 10: Replace `rule all` with dispatch form**

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

---

### Task 2: Verify behavior preservation

- [ ] **Step 1: Run smoke profiles (THE critical check)**

```bash
python3 test/test_stage8_smoke_profiles.py
```

Expected: 7/7 PASS. Scheduling checks must still find all expected rule names.

- [ ] **Step 2: Run Stage 22 BigWig stress tests**

```bash
python3 test/test_stage22_bigwig_stress.py
```

Expected: 6/6 PASS (BigWig gating unchanged).

- [ ] **Step 3: Run validation**

```bash
python3 scripts/validate_samples.py --config config/config.yaml
python3 test/test_validation_stress.py
```

Expected: OK + 15/15 PASS.

- [ ] **Step 4: Run no-hardcoded-paths**

```bash
python3 test/test_no_hardcoded_paths.py
```

Expected: PASS.

- [ ] **Step 5: Run stage-specific stress tests**

```bash
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

- [ ] **Step 6: Final git status and report**

```bash
git status --short
git diff --stat
```
