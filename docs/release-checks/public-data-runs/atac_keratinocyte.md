# ATAC Keratinocyte — Public Data Execution Report

**Status:** not executed yet

This stub follows the [public data execution report template](../public-data-execution-report-template.md).

## Dataset

| Field | Value |
| :--- | :--- |
| Queue | atac_keratinocyte |
| Accession | ENCSR254KDA |
| Source | https://www.encodeproject.org/experiments/ENCSR254KDA/ |
| Assay | atac |
| Target | ATAC |
| Genome | hg38 |
| Layout | PE |
| Treatment replicates | 2 |
| Control samples | 0 |

## Expected Validation

- [ ] ATAC dispatch with Tn5-aware MACS3 parameters visible in logs.
- [ ] Non-zero peak counts in `peak_counts.tsv`.
- [ ] TSS profile shows enrichment peak at TSS (when `qc.tss_enrichment: true`).
- [ ] No IDR outputs (ATAC IDR not supported in v0.2).

## Fill-in Checklist (copy to report before execution)

- [ ] FASTQ files identified / downloaded
- [ ] Reference resources prepared (Bowtie2 index, chrom_sizes, optional blacklist, GTF for TSS)
- [ ] Config and sample sheet written
- [ ] `validate_samples.py --strict-inputs` passes
- [ ] Dry-run resolves without errors
- [ ] Full or downsampled run completes
- [ ] Output audit table filled
- [ ] QC observations recorded
- [ ] Runtime/storage notes recorded
