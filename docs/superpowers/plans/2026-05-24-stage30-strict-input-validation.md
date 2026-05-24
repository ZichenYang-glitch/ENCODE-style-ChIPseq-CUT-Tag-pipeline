# Stage 30: Strict Input Validation — Implementation Plan

> Full design: `docs/superpowers/specs/2026-05-24-stage30-strict-input-validation-design.md`

**Goal:** Add optional `--strict-inputs` validation for FASTQ and Bowtie2 index existence, and MACS3 fragment fallback warnings. No DAG changes, no mandatory behavior.

**Status:** Implemented — all 5 tasks completed; 8 stress tests pass; backward-compatible.

---

### Task 1: Add `--strict-inputs` to validate_samples.py

- [x] Add `--strict-inputs` argument to the CLI parser
- [x] Add `_check_fastq_exists(fastq_path, sample_id, label)` helper
- [x] Add `_check_bowtie2_index(bt2_prefix, sample_id)` helper
  - Check complete `.bt2` or `.bt2l` set; raise `ValidationError` with clear message
- [x] Wire both checks into the sample validation loop, gated on `strict_inputs`
- [x] Pass `strict_inputs` flag through `load_and_validate_samples()`
- [x] Default behavior (no `--strict-inputs`) unchanged

### Task 2: Add MACS3 fragment fallback warning

- [x] In `chipseq.smk`, when fallback is used, emit stderr warning
- [x] No DAG change, no exit code change

### Task 3: Create test/test_stage30_strict_inputs_stress.py

- [x] Test 1: `--strict-inputs` with missing `fastq_1` → validation fails
- [x] Test 2: `--strict-inputs` with PE missing `fastq_2` → validation fails
- [x] Test 3: `--strict-inputs` with missing Bowtie2 index → validation fails
- [x] Test 4: `--strict-inputs` with `.bt2` set present → passes
- [x] Test 5: `--strict-inputs` with `.bt2l` set present → passes
- [x] Test 6: Non-strict mode with placeholder paths → passes
- [x] Test 7: `.fq.gz` existing file → passes
- [x] Test 8: MACS3 fallback warning code present

### Task 4: Verify backward compatibility

- [x] `python3 test/test_validation_stress.py` — 15/15 PASS
- [x] `python3 test/test_stage8_smoke_profiles.py` — 7/7 PASS
- [x] `python3 test/test_no_hardcoded_paths.py` — PASS
- [x] `python3 test/test_stage22_bigwig_stress.py` — 6/6 PASS

### Task 5: Documentation updates

- [x] `README.md` — add `--strict-inputs` to validation section
- [x] `docs/configuration.md` — document strict mode as optional pre-run check
- [x] `RELEASE_CHECKLIST.md` — add `--strict-inputs` as optional real-data step
- [x] `docs/assay-policy.md` — note MACS3 fragment fallback warning

### Files

| File | Action |
|------|--------|
| `scripts/validate_samples.py` | Edit (add --strict-inputs logic) |
| `workflow/rules/chipseq.smk` | Edit (add fallback warning) |
| `test/test_stage30_strict_inputs_stress.py` | Create |
| `README.md` | Edit |
| `docs/configuration.md` | Edit |
| `RELEASE_CHECKLIST.md` | Edit |
| `docs/assay-policy.md` | Edit |
