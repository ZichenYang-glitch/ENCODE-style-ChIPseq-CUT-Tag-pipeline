# Output Contract

This document defines the pipeline's output types, their paths, generating methods, and status. It serves as the reference for downstream consumers (manifest generation, QC reporting, and user documentation). See [`README.md`](../README.md#output-structure) for the directory tree overview.

The manifest field schema is implemented as Stage 25 via `scripts/make_manifest.py`,
producing `results/multiqc/result_manifest.tsv` with 10 columns:
`sample_id`, `experiment_id`, `assay`, `target`, `genome`, `output_type`, `method`, `path`, `status`, `qc_flag`.

The earlier `project_id` field from the draft schema was removed before implementation
— project-level rows use empty `sample_id` and `experiment_id` fields.

## Current Output Types

### Single-sample outputs

| output_type | method | rule | path | status |
| :--- | :--- | :--- | :--- | :--- |
| `final_bam` | bowtie2+samtools | (common.smk) | `results/<sample>/02_align/<sample>.final.bam` | stable |
| `final_bai` | samtools index | (common.smk) | `results/<sample>/02_align/<sample>.final.bam.bai` | stable |
| `cpm_bigwig` | bamCoverage CPM | (common.smk) | `results/<sample>/03_bigwig/<sample>.CPM.bw` | stable |
| `macs3_peak` | macs3 callpeak | (peaks.smk) | `results/<sample>/04_peaks/<sample>/` | stable |
| `macs3_fe_bdg` | macs3 bdgcmp FE | signal_track_fe | `results/<sample>/03_signal/<sample>.FE.bdg` | stable |
| `macs3_ppois_bdg` | macs3 bdgcmp ppois | signal_track_ppois | `results/<sample>/03_signal/<sample>.ppois.bdg` | stable |
| `macs3_fe_bw` | bedGraphToBigWig | signal_track_fe_bw | `results/<sample>/03_signal/<sample>.FE.bw` | **implemented (Stage 22)** |
| `macs3_ppois_bw` | bedGraphToBigWig | signal_track_ppois_bw | `results/<sample>/03_signal/<sample>.ppois.bw` | **implemented (Stage 22)** |
| `peak_counts` | wc -l | peak_counts | `results/<sample>/01_qc/<sample>.peak_counts.tsv` | stable |
| `frip` | calc_frip.py | frip | `results/<sample>/01_qc/<sample>.frip.tsv` | stable |
| `library_complexity` | parse_dup_metrics.py | library_complexity | `results/<sample>/01_qc/<sample>.library_complexity.tsv` | stable |
| `nrf_pbc` | calc_nrf_pbc.py | nrf_pbc | `results/<sample>/01_qc/<sample>.nrf_pbc.tsv` | stable |
| `qc_summary` | assemble_qc_summary | qc_summary | `results/<sample>/01_qc/<sample>.qc_summary.tsv` | **Stage 24 refactored** |
| `cross_correlation` | phantompeakqualtools | cross_correlation | `results/<sample>/05_qc/cross_correlation/<sample>.cc.qc` | stable |
| `cross_correlation_plot` | phantompeakqualtools | cross_correlation | `results/<sample>/05_qc/cross_correlation/<sample>.cc.plot.pdf` | stable |
| `preseq` | preseq lc_extrap | preseq_complexity | `results/<sample>/05_qc/preseq/<sample>.preseq.txt` | stable |
| `picard_alignment` | Picard CollectMultipleMetrics | picard_collect_multiple_metrics | `results/<sample>/05_qc/picard/<sample>.alignment_summary_metrics` | stable |
| `picard_insert_size` | Picard CollectMultipleMetrics | picard_collect_multiple_metrics | `results/<sample>/05_qc/picard/<sample>.insert_size_metrics` | stable |
| `picard_quality` | Picard CollectMultipleMetrics | picard_collect_multiple_metrics | `results/<sample>/05_qc/picard/<sample>.quality_distribution_metrics` | stable |
| `tss_matrix` | computeMatrix | tss_enrichment_profile | `results/<sample>/05_qc/tss/<sample>.tss_matrix.gz` | stable |
| `tss_profile_tsv` | plotProfile | tss_enrichment_profile | `results/<sample>/05_qc/tss/<sample>.tss_profile.tsv` | stable |
| `tss_profile_pdf` | plotProfile | tss_enrichment_profile | `results/<sample>/05_qc/tss/<sample>.tss_profile.pdf` | stable |
| `cuttag_fragment_size` | calc_cuttag_fragment_size.py | cuttag_fragment_size | `results/<sample>/01_qc/<sample>.cuttag_fragment_size.tsv` | stable |
| `seacr_bedgraph` | bedtools genomecov | seacr_bedgraph | `results/<sample>/04_peaks_seacr/<sample>.bedgraph` | stable |
| `seacr_peaks` | SEACR_1.3.sh | seacr_call | `results/<sample>/04_peaks_seacr/<sample>/<sample>.seacr.<mode>.bed` | stable |
| `blacklist_filtered_bam` | bedtools intersect | blacklist_filter_bam | `results/<sample>/02_align/<sample>.blacklist_filtered.bam` | stable |
| `pipeline_done` | touch | (common.smk) | `results/<sample>/logs/<sample>.pipeline.done` | stable |

### Experiment-level (pooled) outputs

| output_type | method | rule | path | status |
| :--- | :--- | :--- | :--- | :--- |
| `pooled_final_bam` | samtools merge | (replicates.smk) | `results/experiments/<exp>/02_align/<exp>.pooled.final.bam` | stable |
| `biorep_final_bam` | samtools merge / symlink | (replicates.smk) | `results/experiments/<exp>/02_align/biorep<N>.final.bam` | stable |
| `pooled_macs3_peak` | macs3 callpeak | (replicates.smk) | `results/experiments/<exp>/04_peaks/pooled/<exp>_pooled_peaks/` | stable |
| `pooled_fe_bdg` | macs3 bdgcmp FE | pooled_signal_track_fe | `results/experiments/<exp>/03_signal/<exp>.pooled.FE.bdg` | stable |
| `pooled_ppois_bdg` | macs3 bdgcmp ppois | pooled_signal_track_ppois | `results/experiments/<exp>/03_signal/<exp>.pooled.ppois.bdg` | stable |
| `pooled_fe_bw` | bedGraphToBigWig | pooled_signal_track_fe_bw | `results/experiments/<exp>/03_signal/<exp>.pooled.FE.bw` | **implemented (Stage 22)** |
| `pooled_ppois_bw` | bedGraphToBigWig | pooled_signal_track_ppois_bw | `results/experiments/<exp>/03_signal/<exp>.pooled.ppois.bw` | **implemented (Stage 22)** |
| `pooled_qc_summary` | pooled_qc_summary.py | pooled_experiment_qc_summary | `results/experiments/<exp>/01_qc/<exp>.pooled_qc_summary.tsv` | stable |

### IDR outputs (Stage 5, TF ChIP-seq narrowPeak only)

| output_type | method | rule | path | status |
| :--- | :--- | :--- | :--- | :--- |
| `idr_true_replicate_raw` | idr | idr_true_replicates | `results/experiments/<exp>/06_idr/true_replicates/idr.txt` | stable |
| `idr_true_replicate_thresh` | idr | idr_true_replicates | `results/experiments/<exp>/06_idr/true_replicates/idr.thresholded.narrowPeak` | stable |
| `idr_self_raw` | idr | idr_self_pseudoreps | `results/experiments/<exp>/06_idr/self_pseudoreplicates/biorep<N>.idr.txt` | stable |
| `idr_pooled_raw` | idr | idr_pooled_pseudoreps | `results/experiments/<exp>/06_idr/pooled_pseudoreplicates/idr.txt` | stable |
| `idr_conservative` | idr filter | stage5b_summary | `results/experiments/<exp>/06_idr/final/conservative.narrowPeak` | stable |
| `idr_optimal` | idr filter | stage5b_summary | `results/experiments/<exp>/06_idr/final/optimal.narrowPeak` | stable |
| `idr_reproducibility_summary` | (script) | stage5b_summary | `results/experiments/<exp>/06_idr/final/reproducibility_summary.tsv` | stable |

### Project-level outputs

| output_type | method | rule | path | status |
| :--- | :--- | :--- | :--- | :--- |
| `multiqc_report` | MultiQC | (report.smk) | `results/multiqc/multiqc_report.html` | stable |
| `stage3_qc_summary` | aggregate_qc_summary | stage3_qc_summary | `results/multiqc/stage3_qc_summary.tsv` | **Stage 24 refactored** |
| `result_manifest` | make_manifest.py | result_manifest | `results/multiqc/result_manifest.tsv` | **Stage 25 implemented** |
| `cross_correlation_summary` | cross_correlation_summary.py | (report.smk) | `results/multiqc/cross_correlation_summary.tsv` | stable |
| `tss_bed` | gtf_to_tss_bed.py | tss_bed_from_gtf | `results/reference/<genome>.tss.bed` | stable |

## Manifest Field Schema (Stage 25, implemented)

Generated by `scripts/make_manifest.py` → `results/multiqc/result_manifest.tsv`.

```tsv
sample_id	experiment_id	assay	target	genome	output_type	method	path	status	qc_flag
```

| Field | Type | Description |
| :--- | :--- | :--- |
| `sample_id` | string | Sample identifier (empty for per-experiment and project-level rows) |
| `experiment_id` | string | Experiment grouping key (empty for per-sample and project-level rows) |
| `assay` | string | `chipseq`, `cuttag`, or `atac` |
| `target` | string | Antibody or target name from sample sheet |
| `genome` | string | Genome assembly label from sample sheet |
| `output_type` | string | Controlled vocabulary (see table above) |
| `method` | string | Tool or rule name (e.g., `macs3_bdgcmp+bedGraphToBigWig`) |
| `path` | string | Absolute path resolved from the configured `outdir` |
| `status` | string | `present`, `missing`, `not_applicable` |
| `qc_flag` | string | Always `NA` in v0.2 (QC threshold auto-judgment deferred) |

## Gating Conditions

| output_type | Gated by |
| :--- | :--- |
| `macs3_fe_bdg`, `macs3_ppois_bdg` | `qc.signal_tracks: true` |
| `macs3_fe_bw`, `macs3_ppois_bw` | `qc.signal_tracks: true` + non-empty `genome_resources.<genome>.chrom_sizes` |
| `pooled_fe_bdg`, `pooled_ppois_bdg` | `stage4b: true` + `qc.signal_tracks: true` + multi-biorep experiment |
| `pooled_fe_bw`, `pooled_ppois_bw` | `stage4b: true` + `qc.signal_tracks: true` + multi-biorep experiment + `chrom_sizes` |
| `blacklist_filtered_bam` | `qc.blacklist_filter: true` + non-empty `genome_resources.<genome>.blacklist` |
| `cross_correlation` | `qc.cross_correlation: true` |
| `preseq` | `qc.preseq_complexity: true` |
| `picard_*` | `qc.picard_metrics: true` + `genome_resources.<genome>.reference_fasta` |
| `tss_*` | `qc.tss_enrichment: true` + `genome_resources.<genome>.gtf` |
| `cuttag_fragment_size` | `qc.cuttag_fragment_size: true` + `assay: cuttag` |
| `seacr_*` | `cuttag.seacr.enabled: true` + `assay: cuttag` + `layout: PE` |
| `idr_*` | `stage5: true` + `chipseq` + `narrow` + exactly 2 treatment bioreps |
| `multiqc_report` | `multiqc: true` |

## Known Assumptions

- **Pooled FE/ppois BigWig genome resolution:** the first treatment sample's genome is used for `chrom_sizes` lookup in pooled experiments. Mixed-genome pooled experiments are not validated against in this slice. If a pooled experiment contains samples with different genomes, the BigWig conversion may use the wrong chromosome sizes.
- **MACS3 bdgcmp sorting:** MACS3 `bdgcmp` produces coordinate-sorted bedGraph in practice, but the BigWig rules defensively re-sort via `LC_ALL=C sort -k1,1 -k2,2n` before calling `bedGraphToBigWig`.
- **chrom_sizes path validation:** when `genome_resources.<genome>.chrom_sizes` is non-empty, `validate_samples.py` checks that the file exists. An empty string means "no BigWig for this genome" — this is not a validation error.
