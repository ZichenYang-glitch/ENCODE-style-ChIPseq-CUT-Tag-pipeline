# TF ChIP CEBPB — Public Data Execution Report

**Status:** not executed yet

This stub follows the [public data execution report template](../public-data-execution-report-template.md).

## Dataset

| Field | Value |
| :--- | :--- |
| Queue | tf_chip_cebpb |
| Accession | ENCSR000DYI |
| Source | https://www.encodeproject.org/experiments/ENCSR000DYI/ |
| Assay | chipseq |
| Target | CEBPB |
| Genome | hg38 |
| Layout | PE |
| Treatment replicates | 2 |
| Control samples | 1 |

## Expected Validation

- [ ] Control, pooled, and IDR paths all produce outputs.
- [ ] `result_manifest.tsv` records `idr_conservative`, `idr_optimal`, `idr_reproducibility_summary`.
- [ ] IDR rescue ratio > 0.5, self-consistency ratio > 0.5.

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
