# Stage 16 Real-Run Report Audit

**Date:** 2026-05-22
**Status:** PASS with small fixes
**Scope:** Stage 15 MultiQC custom-content validation on the external CUT&Tag
real-data run from Stage 13.

## Summary

Stage 16 audited the Stage 15 MultiQC integration against an existing
2-replicate CUT&Tag PE run. No real FASTQ, BAM, BigWig, HTML, or MultiQC data
files are stored in the repository.

## CUT&Tag Report Audit

| Check | Result | Notes |
| :--- | :--- | :--- |
| MultiQC report target resolves | PASS after fix | The target must match the absolute `outdir` configured for the run. |
| Cross-correlation custom section appears | PASS | `multiqc_data.json` contains `Cross-correlation QC` with 2 samples. |
| Cross-correlation summary generated | PASS | Both CUT&Tag samples appear in `cross_correlation_summary.tsv`. |
| NSC/RSC threshold flags | PASS | Both samples are flagged `ok` by project thresholds. |
| Existing report rerun behavior | FIXED | MultiQC now runs with `--force` so reruns overwrite `multiqc_report.html` instead of creating suffixed files. |
| Machine-specific paths in summary | FIXED | `cc_qc_file` now stores the `.cc.qc` basename rather than the input absolute path. |

## Observed CUT&Tag Cross-Correlation Metrics

The external CUT&Tag run produced the following pattern:

| Metric | Observation |
| :--- | :--- |
| `estimated_fragment_length` | 0 for both samples |
| `phantom_peak` | Approximately read length |
| `NSC` | Above 1.05 for both samples |
| `RSC` | Above 0.8 for both samples |
| `quality_flag` | `ok` for both samples |

The fragment-length estimate of 0 is consistent with the observed
phantompeakqualtools behavior on CUT&Tag/Tn5-style data in this run. For
CUT&Tag, interpret cross-correlation as an auxiliary signal-quality metric and
use the CUT&Tag fragment-size QC output for insert-size interpretation.

## Fixes Applied

1. `workflow/rules/report.smk`: added `--force` to the MultiQC command so
   Snakemake reruns produce the declared `multiqc_report.html` output.
2. `scripts/parse_cross_correlation.py`: changed `cc_qc_file` in the summary
   TSV from the full input path to the basename.
3. `test/test_parse_cross_correlation.py`: added a CLI regression test that
   prevents absolute input paths from leaking into `cc_qc_file`.

## Remaining Follow-Ups

- Re-run the ChIP-seq Stage 13 report with the same Stage 15/16 code path when
  convenient. The Stage 13 ChIP-seq workflow completed previously, but this
  audit only rechecked the MultiQC custom section on the CUT&Tag report.
- A native MultiQC plugin remains out of scope. The current custom-content
  section is sufficient for this project unless report packaging needs change.
