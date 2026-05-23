# Stage 25: Minimal Manifest — Implementation Plan

> Full design: `docs/superpowers/specs/2026-05-24-stage25-minimal-manifest-design.md`

**Goal:** Minimal project-level result manifest recording core outputs with existence status. TSV-only. No md5 by default. No QC threshold auto-judgment.

---

### Task 1: Create `scripts/make_manifest.py`

- [ ] Read `config.yaml` via PyYAML (already in python.yml env): `outdir`, `qc.signal_tracks`, `stage4b`, `stage5`, `multiqc`, `genome_resources`
- [ ] Read `samples.tsv` via `csv.DictReader`: treatment samples, experiments, assay, target, genome
- [ ] Build sample rows (deterministic order by `sample_id`):
  - Always: `final_bam`, `final_bai`, `cpm_bigwig`, `macs3_peak`, `qc_summary`
  - `signal_tracks: true`: `macs3_fe_bdg`, `macs3_ppois_bdg`
  - `chrom_sizes` non-empty: `macs3_fe_bw`, `macs3_ppois_bw`
  - Gated-off outputs → `status: not_applicable`
- [ ] Build experiment rows (deterministic order by `experiment_id`, requires `stage4b` + multi-biorep):
  - Always: `pooled_final_bam`, `pooled_final_bai`, `biorep_final_bam`, `biorep_final_bai`, `pooled_macs3_peak`, `pooled_qc_summary`
  - `signal_tracks: true`: `pooled_fe_bdg`, `pooled_ppois_bdg`
  - `chrom_sizes` non-empty: `pooled_fe_bw`, `pooled_ppois_bw`
- [ ] Build IDR rows (requires `stage5: true` + eligible experiments):
  - `idr_conservative`, `idr_optimal`, `idr_reproducibility_summary`
- [ ] Build project-level rows:
  - `stage3_qc_summary`, `multiqc_report` (when enabled)
- [ ] `--md5` flag (off by default): compute hex digest for `present` files
- [ ] `--strict` flag: exit non-zero if any row has `status: missing`
- [ ] Write TSV with `csv.DictWriter`, `lineterminator="\n"`, 10 columns

### Task 2: Add `result_manifest` rule to `workflow/rules/report.smk`

- [ ] Add rule with input depending on `stage3_qc_summary`
- [ ] Shell calls `python3 scripts/make_manifest.py --config config/config.yaml --output {output}`
- [ ] Conda env: `python.yml`
- [ ] Output: `results/multiqc/result_manifest.tsv`

### Task 3: Add manifest target to `workflow/Snakefile`

- [ ] Add to `_base_targets()` helper or create a lightweight `_manifest_targets()` helper
- [ ] Gated on `QC_CONFIG.get("summary", True) and TREATMENT_SAMPLE_IDS`

### Task 4: Create `test/test_stage25_manifest_stress.py`

- [ ] Test 1: Default config → manifest contains expected output types
- [ ] Test 2: `signal_tracks: false` → FE/ppois set to `not_applicable`
- [ ] Test 3: `chrom_sizes` empty → BigWig outputs `not_applicable`
- [ ] Test 4: `stage5: false` → IDR outputs `not_applicable`
- [ ] Test 5: Missing file → `status: missing` (create a config pointing to nonexistent output)
- [ ] Test 6: Deterministic row ordering
- [ ] Test 7: No CRLF in output (`lineterminator="\n"`)
- [ ] Test 8: `--strict` flag exits non-zero on missing
- [ ] Test 9: `--md5` flag populates method field (when file exists and is small enough)

### Task 5: Documentation updates

- [ ] `docs/output-contract.md`: add `result_manifest` as project-level output type
- [ ] `README.md`: note manifest in Developer Notes / Release Checklist
- [ ] `CHANGELOG.md`: Stage 25 entry
- [ ] `KNOWN_ISSUES.md`: mark manifest as completed

### Task 6: Full verification

```bash
python3 scripts/validate_samples.py --config config/config.yaml
python3 test/test_validation_stress.py
python3 test/test_stage8_smoke_profiles.py
python3 test/test_stage22_bigwig_stress.py
python3 test/test_stage24_qc_summary_unit.py
python3 test/test_stage25_manifest_stress.py
python3 test/test_no_hardcoded_paths.py
```

## Files

| File | Action |
|------|--------|
| `scripts/make_manifest.py` | Create |
| `workflow/rules/report.smk` | Edit (add `result_manifest` rule) |
| `workflow/Snakefile` | Edit (add manifest target) |
| `test/test_stage25_manifest_stress.py` | Create |
| `docs/output-contract.md` | Edit |
| `README.md` | Edit |
| `CHANGELOG.md` | Edit |
| `KNOWN_ISSUES.md` | Edit |
