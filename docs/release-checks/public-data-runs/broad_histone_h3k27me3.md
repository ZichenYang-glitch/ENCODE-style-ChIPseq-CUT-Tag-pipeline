# Broad Histone H3K27me3 — Public Data Execution Report

**Status:** not executed yet

This stub follows the [public data execution report template](../public-data-execution-report-template.md).

## Dataset

| Field | Value |
| :--- | :--- |
| Queue | broad_histone_h3k27me3 |
| Accession | ENCSR000AKB |
| Source | https://www.encodeproject.org/experiments/ENCSR000AKB/ |
| Assay | chipseq |
| Target | H3K27me3 |
| Genome | hg38 |
| Layout | PE |
| Treatment replicates | 2 |
| Control samples | 1 |

## Expected Validation

- [ ] MACS3 broadPeak output, not narrowPeak.
- [ ] `pooled_qc_summary.tsv` reports `inferred_histone_class=broad_like`, `peak_mode_status=ok`.
- [ ] No IDR outputs (broad marks not eligible in v0.2).

## Fill-in Checklist (copy to report before execution)

- [ ] FASTQ files identified / downloaded
- [ ] Reference resources prepared (Bowtie2 index, chrom_sizes, optional blacklist)
- [ ] Config and sample sheet written
- [ ] `validate_samples.py --strict-inputs` passes
- [ ] Dry-run resolves without errors
- [ ] Full or downsampled run completes
- [ ] Output audit table filled
- [ ] QC observations recorded
- [ ] Runtime/storage notes recorded
