# CUT&Tag H3K27me3 — Public Data Execution Report

**Status:** not executed yet

This stub follows the [public data execution report template](../public-data-execution-report-template.md).

## Dataset

| Field | Value |
| :--- | :--- |
| Queue | cuttag_h3k27me3 |
| Accession | GSE145187 |
| Source | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145187 |
| Assay | cuttag |
| Target | H3K27me3 (or H3K4me3) |
| Genome | hg38 |
| Layout | PE |
| Treatment replicates | 2 |
| Control samples | 1 (IgG) |

## Expected Validation

- [ ] `cuttag_fragment_size.tsv` reports sub-200 bp enrichment.
- [ ] SEACR outputs under `04_peaks_seacr/` (when `cuttag.seacr.enabled: true`).
- [ ] Low duplication rate (<20%).
- [ ] No IDR outputs (CUT&Tag IDR not supported in v0.2).

## Fill-in Checklist (copy to report before execution)

- [ ] FASTQ/SRA files identified / downloaded
- [ ] Reference resources prepared (Bowtie2 index, chrom_sizes, optional blacklist)
- [ ] Config and sample sheet written
- [ ] `validate_samples.py --strict-inputs` passes
- [ ] Dry-run resolves without errors
- [ ] Full or downsampled run completes
- [ ] Output audit table filled
- [ ] QC observations recorded
- [ ] Runtime/storage notes recorded
