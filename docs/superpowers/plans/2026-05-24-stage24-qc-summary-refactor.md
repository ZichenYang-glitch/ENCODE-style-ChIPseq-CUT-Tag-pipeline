# Stage 24: QC Summary Refactor — Implementation Plan

> Full design: `docs/superpowers/specs/2026-05-24-stage24-qc-summary-refactor-design.md`

**Goal:** Replace shell `printf`/`tail`/`cut` in `qc_summary` and `stage3_qc_summary` with Python scripts using `csv.DictReader`. Preserve 37-column output contract byte-compatibly.

---

### Task 1: Create `scripts/assemble_qc_summary.py`

- [ ] Define `_QC_SUMMARY_COLUMNS` (37 columns in exact current order)
- [ ] Parse metadata args (`--sample`, `--assay`, `--target`, `--genome`, `--layout`, `--peak-mode`, `--use-control`, `--control-type`, `--final-bam`, `--peaks-file`, `--has-blacklist`, `--blacklist`, `--bl-bam`, `--bl-peaks`)
- [ ] Parse input TSVs via `csv.DictReader`: `--peak-counts`, `--frip`, `--library-complexity`, `--nrf-pbc`
- [ ] Extract named columns from each input; missing/empty → `"NA"`
- [ ] Blacklist fields → `"NA"` when `has_blacklist == "no"`
- [ ] Write header + single data row in exact 37-column order
- [ ] Required input file missing → `sys.exit(1)` with clear error

### Task 2: Create `scripts/aggregate_qc_summary.py`

- [ ] Use same `_QC_SUMMARY_COLUMNS` constant
- [ ] Accept positional input files + `--output`
- [ ] Zero inputs → write header-only file
- [ ] Validate all input file headers match `_QC_SUMMARY_COLUMNS` exactly; fail with error on mismatch
- [ ] Concatenate data rows preserving order

### Task 3: Update `workflow/rules/qc.smk`

- [ ] Replace `qc_summary` shell block (lines 696-823) with `python3 scripts/assemble_qc_summary.py` call
- [ ] Replace `stage3_qc_summary` shell block (lines 1009-1060) with `python3 scripts/aggregate_qc_summary.py` call
- [ ] Keep all `input:` and `params:` declarations unchanged
- [ ] Keep output paths unchanged

### Task 4: Create `test/test_stage24_qc_summary_unit.py`

- [ ] Test 1: Complete inputs → correct 37-column TSV
- [ ] Test 2: NA blacklist fields when has_blacklist=no
- [ ] Test 3: NA library complexity when Picard unavailable
- [ ] Test 4: Header-only aggregation when zero inputs
- [ ] Test 5: Header mismatch detection in aggregation
- [ ] Test 6: Missing input file → error exit
- [ ] Test 7: Empty peak_counts (0 peaks) handled
- [ ] Test 8: Byte-compare header with expected 37-column header

### Task 5: Documentation

- [ ] Edit `CHANGELOG.md` — Stage 24 entry under `[Unreleased]`

### Task 6: Full verification

```bash
python3 scripts/validate_samples.py --config config/config.yaml
python3 test/test_validation_stress.py
python3 test/test_stage8_smoke_profiles.py
python3 test/test_stage22_bigwig_stress.py
python3 test/test_no_hardcoded_paths.py
python3 test/test_stage6b_stress.py
python3 test/test_stage12_stress.py
python3 test/test_stage24_qc_summary_unit.py
```
