# Stage 24: QC Summary Refactor — Design Spec

**Date:** 2026-05-24
**Status:** design review
**Scope:** Migrate `qc_summary` and `stage3_qc_summary` shell assembly to Python scripts with header-based parsing
**Excluded:** manifest script, frontend, MultiQC plugin, rule all changes, pooled_qc_summary changes, column reordering

---

## 1. Goals

1. Replace the shell `printf` + `tail/cut` logic in `qc_summary` and `stage3_qc_summary` rules with Python scripts that parse input TSVs by header.
2. Gracefully handle missing optional metric values within input files with `NA` (required input files must exist — Snakemake enforces this).
3. Preserve the exact 37-column output contract (same column names, same column order).
4. Enable unit testing of summary assembly logic (shell printf cannot be unit-tested).
5. No pandas — stdlib `csv` only.

## 2. Current State: Shell-Assembled Summaries

| Rule | Lines in qc.smk | Method |
| :--- | :--- | :--- |
| `qc_summary` | 659–823 | Shell `printf` header + `tail -n +2 \| cut -fN` to extract columns |
| `stage3_qc_summary` | 1000–1060 | Shell `head -n 1` + `tail -n +2` loop to concatenate |

`pooled_experiment_qc_summary` already uses `scripts/pooled_qc_summary.py` — unchanged.

## 3. Output Contract: 37-Column qc_summary.tsv

The exact column order (preserved from current shell header, lines 738–773):

```
sample, assay, target, genome, layout, peak_mode,
use_control, control_type, final_bam, peaks,
blacklist, blacklist_filtered_bam, blacklist_filtered_peaks,
total_reads, reads_in_peaks, frip, peak_count, blacklist_filtered_peak_count,
metrics_source, unpaired_reads_examined, read_pairs_examined,
secondary_or_supplementary_reads, unmapped_reads, unpaired_read_duplicates,
read_pair_duplicates, read_pair_optical_duplicates, percent_duplication,
estimated_library_size, total_reads_examined, duplicate_reads_estimate,
total_fragments, distinct_fragments, one_read_fragments, two_read_fragments,
nrf, pbc1, pbc2
```

This is the 37-column constant shared between `assemble_qc_summary.py` and `aggregate_qc_summary.py`.

## 4. Required vs Optional

**Required (file must exist — Snakemake enforces):**
- `peak_counts.tsv`, `frip.tsv`, `library_complexity.tsv`, `nrf_pbc.tsv`, `final.bam`

**Optional within required files (missing/empty fields → `NA`):**
- Blacklist fields when `has_blacklist == "no"` (preserving existing behavior)
- Library complexity fields when Picard is unavailable (already `NA` in input TSV)
- FRiP when `total_reads == 0`
- NRF/PBC fields when BAM has zero fragments

## 5. Architecture

### `scripts/assemble_qc_summary.py`

Takes sample metadata + paths to 4 input TSVs. Reads each via `csv.DictReader`, extracts named columns, writes single-row 37-column TSV.

### `scripts/aggregate_qc_summary.py`

Takes output path + variable input file paths. Validates header consistency across all inputs. Concatenates data rows. If zero inputs, writes header-only file.

### Shared header constant

Both scripts use the same `_QC_SUMMARY_COLUMNS` list. `aggregate_qc_summary.py` validates that each input file's header matches exactly, failing with a clear error on mismatch.

## 6. Out of Scope

- `scripts/make_manifest.py` (Stage 25)
- Column additions, removals, or reordering
- `pooled_experiment_qc_summary` changes
- Stage 22 BigWig or Stage 23 target builder changes
- pandas, numpy, or other heavy dependencies

## 7. Verification

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
