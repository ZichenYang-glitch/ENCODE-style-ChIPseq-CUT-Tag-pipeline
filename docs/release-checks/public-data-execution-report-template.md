# Public Data Execution Report — <queue_name>

**Date:** (fill in)
**Status:** (not executed / smoke dry-run / downsampled run / full run)
**Reporter:** (name / initials)

---

## 1. Dataset Metadata

| Field | Value |
| :--- | :--- |
| Queue | (e.g., tf_chip_cebpb) |
| Accession | (ENCODE or GEO accession) |
| Source URL | (link to ENCODE/GEO record) |
| Assay | (chipseq / cuttag / atac) |
| Target | (CEBPB / H3K27me3 / ATAC / ...) |
| Genome | (hg38 / mm10 / ...) |
| Layout | (PE / SE) |
| Treatment replicates | (number of treatment bioreps) |
| Control samples | (number of control samples and type) |

---

## 2. Input Preparation

| Item | Detail |
| :--- | :--- |
| FASTQ file IDs / SRR runs | (ENCODE file accessions, GEO SRR numbers, or local paths) |
| Downsampling strategy | (none / N reads per sample / subsampling tool) |
| Bowtie2 index | (path or generation command) |
| chrom_sizes | (path) |
| Blacklist BED | (path, or "not configured") |
| GTF | (path, or "not configured") |
| Reference FASTA + .fai + .dict | (path, or "not configured") |

---

## 3. Run Configuration

| Item | Detail |
| :--- | :--- |
| Config path | (e.g., config/config.tf_chip_cebpb.yaml) |
| Sample sheet path | (e.g., config/samples.tf_chip_cebpb.tsv) |
| Command | (snakemake -s workflow/Snakefile ...) |
| Conda env | (runner env + rule-specific envs) |
| Hardware | (CPU cores, RAM, disk) |
| Runtime | (wall clock time) |

---

## 4. Validation Checklist

- [ ] `python3 scripts/validate_samples.py --config <config>`
- [ ] `python3 scripts/validate_samples.py --config <config> --strict-inputs`
- [ ] Dry-run smoke profiles pass (7/7)
- [ ] Full or downsampled execution completes without error
- [ ] `result_manifest.tsv` reviewed
- [ ] `qc_summary.tsv` values reviewed
- [ ] BigWig outputs present (when chrom_sizes configured)
- [ ] IDR outputs present and reproducible (when stage5 enabled)

---

## 5. Output Audit Table

| output_type | expected | observed | status | notes |
| :--- | :--- | :--- | :--- | :--- |
| final_bam | present | | | |
| final_bai | present | | | |
| cpm_bigwig | present | | | |
| macs3_peak | present | | | |
| qc_summary | present | | | |
| macs3_fe_bdg | present | | | |
| macs3_ppois_bdg | present | | | |
| macs3_fe_bw | gated | | | |
| macs3_ppois_bw | gated | | | |
| pooled_final_bam | gated | | | |
| pooled_macs3_peak | gated | | | |
| pooled_qc_summary | gated | | | |
| idr_conservative | gated | | | |
| idr_optimal | gated | | | |
| stage3_qc_summary | present | | | |
| multiqc_report | present | | | |
| result_manifest | present | | | |

---

## 6. QC Observations

(Selected qc_summary.tsv metrics of interest, or notes on
cross-correlation, FRiP, library complexity, TSS profiles, etc.)

---

## 7. Runtime / Storage Observations

(Disk usage, memory peaks, I/O notes, unexpected slowness.)

---

## 8. Known Deviations

(Differences from expected behavior, skipped outputs, configuration
variations, OS / filesystem quirks.)

---

## 9. Conclusion

**Outcome:** pass / pass with caveats / fail

**Blocker issues:** (list or "none")

**Follow-up actions:** (list or "none")
