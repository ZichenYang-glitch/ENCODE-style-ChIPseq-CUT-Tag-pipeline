# Stage 39: MNase-seq Positioning MVP — Design Spec

**Date:** 2026-06-04
**Status:** implemented (Stage 39 completed; 32/32 MNase stress tests pass)
**Scope:** Add `assay: mnase` for PE MNase-seq nucleosome positioning basics
**Excluded:** DANPOS3, iNPS, sub/di fragment classes, rotational QC, NFR/TSS MNase QC, MNase MultiQC custom content

---

## 1. Goals

1. Add `assay: mnase` as a first-class assay alongside `chipseq`, `cuttag`, and `atac`.
2. Add `peak_mode: nucleosome` — required for MNase, rejected for all other assays.
3. MNase is PE-only. SE MNase raises a validation error.
4. Reuse existing preprocessing (FastQC, Trim Galore, Bowtie2, MAPQ filter, duplicate handling, final.bam) unchanged.
5. Skip MACS3 peaks, peak counts, FRiP, MACS3 FE/ppois signal, per-sample qc_summary, pooled peaks, pooled QC summary, and IDR for MNase samples.
6. Keep library complexity, NRF/PBC, CPM BigWig, cross-correlation, preseq, and Picard metrics for MNase samples where configured.
7. Produce MNase-specific outputs: mono-nucleosome BAM, dyad BigWig, mono occupancy BigWig — at both sample and pooled experiment level.
8. Zero new Conda environments — `alignmentSieve` and `bamCoverage` are in existing `deeptools.yml`.
9. Config: optional `mnase.mono_range` (default `[140, 200]`).

---

## 2. Design Decisions

### 2.1 `peak_mode: nucleosome`

**Decision: add `nucleosome` to the peak_mode allowlist; required for MNase; rejected for non-MNase assays.**

The field name `peak_mode` is legacy — it describes the downstream analysis unit. Adding `nucleosome` is the least-invasive approach:
- No new required sample sheet column.
- No need to rename `peak_mode` (a breaking schema change).
- Cross-row consistency checks (assay/target/genome/peak_mode/layout must match within replicate groups) continue to work automatically.

### 2.2 Fragment class split strategy: Option B

Mono BAM via `alignmentSieve` (140–200 bp). Dyad BigWig via `bamCoverage --MNase --binSize 1` from final.bam. Mono occupancy BigWig via `bamCoverage` from mono BAM.

Sub-nucleosome and di-nucleosome layers deferred to Stage 40.

### 2.3 Config: `mnase:` block

Stage 39 has exactly one key:

```yaml
mnase:
  mono_range: [140, 200]
```

The block is optional. If absent, `mono_range` defaults to `[140, 200]`.

`dyad_range` is not needed: `bamCoverage --MNase` uses its own default (130–200 bp) which is correct and well-established.

### 2.4 Target builder strategy

Keep `TREATMENT_SAMPLE_IDS` as the universal set. Introduce filtered subsets:

- `PEAK_SAMPLE_IDS` = treatment samples where assay != "mnase"
- `MNASE_SAMPLE_IDS` = treatment samples where assay == "mnase"
- `PEAK_MULTI_BIOREP_EXPERIMENTS` = multi-biorep experiments whose first treatment is not MNase
- `MNASE_MULTI_BIOREP_EXPERIMENTS` = multi-biorep experiments whose first treatment is MNase

Peak-centric target expansions use `PEAK_SAMPLE_IDS` / `PEAK_MULTI_BIOREP_EXPERIMENTS`.
MNase-centric target expansions use `MNASE_SAMPLE_IDS` / `MNASE_MULTI_BIOREP_EXPERIMENTS`.

### 2.5 Dispatch functions

- `get_remove_dup`: add MNase branch (PE-only, `auto` → `yes`, no broad-peak exception).
- `get_extend_reads`: add MNase branch (PE → no extension; fragments are the biological unit).
- `get_macs3_args`: add MNase branch that raises a clear error (MNase should never reach MACS3).

---

## 3. New Outputs

### Sample level

| Path | Tool |
| :--- | :--- |
| `<outdir>/<s>/03_fragments/<s>.mono.bam` | `alignmentSieve` |
| `<outdir>/<s>/03_fragments/<s>.mono.bam.bai` | `samtools index` |
| `<outdir>/<s>/04_signal/<s>.dyad.CPM.bw` | `bamCoverage --MNase --binSize 1` |
| `<outdir>/<s>/04_signal/<s>.mono.CPM.bw` | `bamCoverage` from mono BAM |

### Pooled experiment level (>=2 treatment bioreps)

| Path | Tool |
| :--- | :--- |
| `<outdir>/experiments/<e>/03_fragments/<e>.pooled.mono.bam` | `alignmentSieve` from pooled BAM |
| `<outdir>/experiments/<e>/03_fragments/<e>.pooled.mono.bam.bai` | `samtools index` |
| `<outdir>/experiments/<e>/04_signal/<e>.pooled.dyad.CPM.bw` | `bamCoverage --MNase --binSize 1` from pooled BAM |
| `<outdir>/experiments/<e>/04_signal/<e>.pooled.mono.CPM.bw` | `bamCoverage` from pooled mono BAM |

---

## 4. Excluded Outputs for MNase

| Output | Reason |
| :--- | :--- |
| `04_peaks/` | No MACS3 peak calling for MNase |
| `01_qc/peak_counts.tsv` | No peaks |
| `01_qc/frip.tsv` | No peaks |
| `01_qc/qc_summary.tsv` | Depends on peak_counts + frip |
| `03_signal/FE.bdg`, `03_signal/ppois.bdg`, `.bw` | MACS3-derived |
| `experiments/<e>/04_peaks/pooled/` | No pooled MACS3 for MNase |
| `experiments/<e>/01_qc/pooled_qc_summary.tsv` | Takes pooled peaks as input |
| `06_idr/` | No IDR for MNase |

---

## 5. Validation Rules

| Check | Condition | Error |
| :--- | :--- | :--- |
| Assay allowlist | `assay` not in `chipseq, cuttag, atac, mnase` | Invalid assay |
| MNase PE-only | `assay == "mnase"` and `layout != "PE"` | Requires PE |
| peak_mode allowlist | `peak_mode` not in `narrow, broad, nucleosome` | Invalid peak_mode |
| MNase nucleosome | `assay == "mnase"` and `peak_mode != "nucleosome"` | MNase requires nucleosome |
| Non-MNase nucleosome | `assay != "mnase"` and `peak_mode == "nucleosome"` | nucleosome is MNase-only |
| mono_range shape | `mnase.mono_range` present, not `[int, int]` with `a < b` | Invalid mono_range |
| mono_range defaults | `mnase` block absent | `mono_range` defaults to `[140, 200]` |

---

## 6. Files to Modify / Create

| File | Action |
| :--- | :--- |
| `workflow/schemas/samples.schema.yaml` | Modify: assay + `mnase`, peak_mode + `nucleosome` |
| `workflow/schemas/config.schema.yaml` | Modify: add `mnase:` block |
| `scripts/validate_samples.py` | Modify: assay, peak_mode, layout, mnase config validation |
| `workflow/rules/mnase.smk` | **Create**: policy + executable rules |
| `workflow/Snakefile` | Modify: derived lists, dispatch, target builder, include mnase.smk |
| `docs/assay-policy.md` | Modify: add MNase-seq section |
| `README.md` | Modify: title, features, assay support, limitations |
| `test/profiles/mnase_pe_noctrl/` | **Create**: config.yaml, samples.tsv |
| `test/test_stage8_smoke_profiles.py` | Modify: add MNase profile |
| `test/test_validation_stress.py` | Modify: add MNase validation cases |
| `test/test_stage39_mnase_stress.py` | **Create**: target builder + DAG tests |
| `test/test_stage25_manifest_stress.py` | Modify: MNase manifest rows |
| `test/test_no_hardcoded_paths.py` | Existing guard passes with MNase additions |
| `docs/superpowers/specs/2026-06-04-stage39-mnase-mvp-design.md` | **Create** |
| `docs/superpowers/plans/2026-06-04-stage39-mnase-mvp.md` | **Create** |

---

## 7. Test Plan

### 7.1 Validation stress (extend existing)

| # | Case | Expected |
|---|---|---|
| 1 | `assay=mnase, layout=PE, peak_mode=nucleosome` → passes | PASS |
| 2 | `assay=mnase, layout=SE` → error | FAIL |
| 3 | `assay=mnase, peak_mode=narrow` → error | FAIL |
| 4 | `assay=chipseq, peak_mode=nucleosome` → error | FAIL |
| 5 | `assay=cuttag, peak_mode=nucleosome` → error | FAIL |
| 6 | `mnase.mono_range: [100, 180]` in config → passes | PASS |
| 7 | `mnase.mono_range: [200, 100]` → error | FAIL |
| 8 | `mnase.mono_range: [140]` → error | FAIL |

### 7.2 Smoke profile

Add `test/profiles/mnase_pe_noctrl/` with 2 MNase treatment bioreps, no controls. Dry-run verifies:
- MNase targets scheduled
- Peak targets not scheduled

### 7.3 Stage 39 stress

| # | Case | Expected |
|---|---|---|
| 1 | `MNASE_SAMPLE_IDS` excludes non-MNase | PASS |
| 2 | `PEAK_SAMPLE_IDS` excludes MNase | PASS |
| 3 | MNase sample dry-run includes mono.bam, dyad.bw, mono.bw | PASS |
| 4 | MNase sample dry-run excludes 04_peaks/, peak_counts, frip | PASS |
| 5 | 2-biorep MNase experiment schedules pooled MNase outputs | PASS |
| 6 | MNase experiment does NOT schedule pooled peaks or pooled QC summary | PASS |
| 7 | IDR not scheduled for MNase experiments | PASS |
| 8 | `mnase.mono_range` config drives alignmentSieve params | PASS |
| 9 | MNase profile integrates with smoke test harness | PASS |

### 7.4 Regression safety

All existing suites must pass unchanged:
- `test_validation_stress.py`
- `test_stage8_smoke_profiles.py` (plus new MNase profile)
- `test_stage22_bigwig_stress.py`
- `test_stage24_qc_summary_unit.py`
- `test_stage25_manifest_stress.py`
- `test_no_hardcoded_paths.py`
- `test_stage27_public_validation_plan.py`
- `test_stage27b_metadata_ci_plan.py`
- `test_stage27c_ci_workflow.py`
- `test_stage28_release_readiness.py`

---

## 8. Risks

| Risk | Mitigation |
| :--- | :--- |
| Target builder regression | `TREATMENT_SAMPLE_IDS` unchanged; all 8 smoke profiles must pass |
| Dispatch ValueError for MNase | Add MNase branch to all 3 dispatch functions |
| `_bamcoverage_inputs` assay check | MNase is PE-only; Chip-seq SE fragment check never fires for MNase |
| `stage3_qc_summary` includes MNase | Gate on `PEAK_SAMPLE_IDS` |
| Doc overclaiming | README says "PE MNase-seq nucleosome positioning basics" with explicit deferred list |
