# Stage 39: MNase-seq Positioning MVP — Implementation Plan

> Full design: `docs/superpowers/specs/2026-06-04-stage39-mnase-mvp-design.md`

**Goal:** Add `assay: mnase` for PE MNase-seq nucleosome positioning basics. Reuse preprocessing, skip peak-centric outputs, add mono BAM + dyad/mono BigWig.

**Status:** Implemented — Stage 39 MNase-seq MVP complete; 32/32 MNase stress tests pass.

---

### Task 1: Update schema files

- [x] `workflow/schemas/samples.schema.yaml`: add `"mnase"` to `assay.allowed`; add `"nucleosome"` to `peak_mode.allowed`; add MNase assay constraint note
- [x] `workflow/schemas/config.schema.yaml`: add `mnase:` block with `mono_range` property

### Task 2: Update validate_samples.py

- [x] Add `"mnase"` to assay allowlist (line 1068)
- [x] Add MNase PE-only constraint after layout validation
- [x] Add `"nucleosome"` to peak_mode allowlist (line 1077)
- [x] Add MNase <=> nucleosome mutual requirement after ATAC constraint
- [x] Add `_validate_mnase_config()` helper called from `validate_config()`
- [x] Update docstrings and CLI description to mention MNase

### Task 3: Create workflow/rules/mnase.smk

- [x] Policy functions: `get_remove_dup_mnase`, `get_extend_reads_mnase`, `get_macs3_args_mnase`
- [x] `mnase_split_mono` rule: `alignmentSieve` from final.bam → mono.bam
- [x] `mnase_dyad_bigwig` rule: `bamCoverage --MNase --binSize 1` from final.bam
- [x] `mnase_mono_bigwig` rule: `bamCoverage` from mono.bam
- [x] `mnase_pooled_mono` rule: `alignmentSieve` from pooled final.bam
- [x] `mnase_pooled_dyad_bigwig` rule: `bamCoverage --MNase --binSize 1` from pooled final.bam
- [x] `mnase_pooled_mono_bigwig` rule: `bamCoverage` from pooled mono.bam
- [x] MNase config helper: `get_mono_range` reading from config

### Task 4: Modify workflow/Snakefile

- [x] Import mnase policy functions
- [x] Add `MNASE_SAMPLE_IDS`, `PEAK_SAMPLE_IDS` derived lists
- [x] Add `MNASE_MULTI_BIOREP_EXPERIMENTS`, `PEAK_MULTI_BIOREP_EXPERIMENTS`
- [x] Wire MNase into dispatch: `get_remove_dup`, `get_extend_reads`, `get_macs3_args`
- [x] Add `_mnase_targets()` helper: sample + pooled MNase outputs
- [x] Update `rule all` to include `_mnase_targets()`
- [x] Gate peak-centric expansions on PEAK_SAMPLE_IDS/PEAK_MULTI_BIOREP_EXPERIMENTS
- [x] Gate `stage3_qc_summary` on `PEAK_SAMPLE_IDS`
- [x] Include `workflow/rules/mnase.smk`

### Task 5: Update documentation

- [x] `docs/assay-policy.md`: add MNase-seq section
- [x] `README.md`: update title, features, assay support, limitations
- [x] MNase scope boundaries documented in README and `docs/assay-policy.md`

### Task 6: Create test profile

- [x] `test/profiles/mnase_pe_noctrl/config.yaml`: minimal config with MNase
- [x] `test/profiles/mnase_pe_noctrl/samples.tsv`: 2 MNase treatment bioreps, PE, nucleosome

### Task 7: Extend test files

- [x] `test/test_validation_stress.py`: add MNase validation cases (assay, layout, peak_mode, mono_range)
- [x] `test/test_stage8_smoke_profiles.py`: add MNase profile (8th profile)
- [x] `test/test_stage39_mnase_stress.py`: new — target builder + DAG tests (9 cases)
- [x] `test/test_no_hardcoded_paths.py`: existing hardcoded-path guard passes with MNase additions

### Task 8: Run tests and report

- [x] Run validation stress
- [x] Run smoke profiles (all 8)
- [x] Run Stage 39 stress
- [x] Run no-hardcoding guard
- [x] Run regression suites: manifest, smoke profiles, validation stress, no-hardcoding guard, release readiness
- [x] Fix any failures
- [x] Report final results
