# Stage 27a: Public Data Validation Plan

**Date:** 2026-05-24
**Status:** planned (no data downloaded, no execution yet)

## Validation Goals

1. Verify that the v0.2 pipeline produces correct, interpretable outputs on real public data.
2. Cover all three supported assays (ChIP-seq, CUT&Tag, ATAC-seq) plus TF ChIP IDR.
3. Keep the validation matrix small, fixed, and reproducible — not an exhaustive benchmark.
4. Zero large data files committed to the repository.

## Candidate Datasets

| Queue | Dataset | Accession | Assay / Target | Replicates | Use |
| :--- | :--- | :--- | :--- | :--- | :--- |
| TF ChIP | ENCODE A549 CEBPB | ENCSR000DYI (or equivalent) | ChIP-seq / CEBPB (TF) | 2 treatment + 1 control | control, pooled, IDR |
| Broad histone ChIP | ENCODE K562 H3K27me3 | ENCSR000AKB (or equivalent) | ChIP-seq / H3K27me3 (broad) | 2 treatment + 1 control | broad peak, pooled signal, histone policy |
| ATAC | ENCODE primary keratinocyte ATAC | ENCSR254KDA (or equivalent) | ATAC-seq / — | 2 replicates | ATAC dispatch, TSS/FRiP |
| CUT&Tag | GSE145187 H3K27me3 or H3K4me3 + IgG | GSE145187 | CUT&Tag / H3K27me3 or H3K4me3 | 2 treatment + IgG control | SEACR sidecar, fragment-size QC |

### Dataset Selection Notes

- All candidates are from ENCODE or GEO and are publicly accessible without authentication.
- If primary candidates are unavailable, equivalent alternatives with the same assay/target/replicate structure may be substituted.
- The CUT&Tag candidate may need alternative sourcing if GEO access patterns change.
- Ignores any pre-computed peak calls or processed data in the original submissions — the pipeline re-derives everything from FASTQ or BAM.

## Genome / Build Assumptions

| Queue | Genome | Effective Genome Size | Notes |
| :--- | :--- | :--- | :--- |
| TF ChIP (CEBPB) | hg38 / GRCh38 | 2913022398 | ENCODE standard for A549 |
| Broad histone (H3K27me3) | hg38 / GRCh38 | 2913022398 | ENCODE standard for K562 |
| ATAC | hg38 / GRCh38 | 2913022398 | ENCODE standard for keratinocyte |
| CUT&Tag | hg38 / GRCh38 or mm10 / GRCm38 | depends on dataset | Check GEO metadata |

## Required Reference Resources

Per-genome resources needed (paths configured by the user, not committed):

- Bowtie2 index (`.1.bt2` through `.4.bt2`, plus `.rev.*`)
- Chromosome sizes file (`.chrom.sizes`)
- ENCODE blacklist BED (optional; enables blacklist QC)
- Reference FASTA with `.fai` and `.dict` (optional; enables Picard metrics)
- GTF annotation (optional; enables TSS enrichment QC)

Users prepare these once via [`docs/reference-resources.md`](../reference-resources.md).

## Expected Config / Sample Sheet Shape

Each validation queue gets its own `config/config.<queue>.yaml` and `config/samples.<queue>.tsv` (not committed, or committed only as minimal examples with placeholder paths).

Example for TF ChIP:

```yaml
# config/config.tf_chip_cebpb.yaml
samples: "samples.tf_chip_cebpb.tsv"
outdir: "results_tf_chip_cebpb"
threads: 16
genome_resources:
  hg38:
    effective_genome_size: 2913022398
    chrom_sizes: "/path/to/hg38.chrom.sizes"
    blacklist: "/path/to/hg38.blacklist.bed"
stage4b: true
stage5: true
idr:
  seed: 42
  threshold: 0.05
qc:
  signal_tracks: true
  cross_correlation: true
  picard_metrics: true
  tss_enrichment: true
```

```tsv
# config/samples.tf_chip_cebpb.tsv
sample  fastq_1  fastq_2  layout  assay  target  peak_mode  genome  bowtie2_index  experiment  biological_replicate  role  control_sample
CEBPB_rep1  /data/fastq/CEBPB_rep1_R1.fq.gz  /data/fastq/CEBPB_rep1_R2.fq.gz  PE  chipseq  CEBPB  narrow  hg38  /path/to/bt2/hg38  CEBPB  1  treatment  CEBPB_ctrl
CEBPB_rep2  /data/fastq/CEBPB_rep2_R1.fq.gz  /data/fastq/CEBPB_rep2_R2.fq.gz  PE  chipseq  CEBPB  narrow  hg38  /path/to/bt2/hg38  CEBPB  2  treatment  CEBPB_ctrl
CEBPB_ctrl  /data/fastq/CEBPB_ctrl_R1.fq.gz  /data/fastq/CEBPB_ctrl_R2.fq.gz  PE  chipseq  CEBPB  narrow  hg38  /path/to/bt2/hg38  CEBPB  1  control
```

## Expected Outputs to Audit

For each validation queue, after a successful run, verify that the following outputs exist and contain interpretable values:

### All queues
- [ ] `results_<queue>/<sample>/02_align/<sample>.final.bam` + `.bai`
- [ ] `results_<queue>/<sample>/03_bigwig/<sample>.CPM.bw`
- [ ] `results_<queue>/<sample>/04_peaks/<sample>/` (MACS3 peaks)
- [ ] `results_<queue>/<sample>/01_qc/<sample>.qc_summary.tsv` (37 columns)
- [ ] `results_<queue>/multiqc/stage3_qc_summary.tsv`
- [ ] `results_<queue>/multiqc/result_manifest.tsv` (10 columns)
- [ ] `results_<queue>/multiqc/multiqc_report.html`

### Signal tracks (when `chrom_sizes` configured)
- [ ] `results_<queue>/<sample>/03_signal/<sample>.FE.bdg` + `.FE.bw`
- [ ] `results_<queue>/<sample>/03_signal/<sample>.ppois.bdg` + `.ppois.bw`

### TF ChIP (CEBPB)
- [ ] `results_<queue>/experiments/CEBPB/02_align/CEBPB.pooled.final.bam`
- [ ] `results_<queue>/experiments/CEBPB/06_idr/final/conservative.narrowPeak`
- [ ] `results_<queue>/experiments/CEBPB/06_idr/final/optimal.narrowPeak`
- [ ] `results_<queue>/experiments/CEBPB/06_idr/final/reproducibility_summary.tsv`
- [ ] IDR rescue ratio > 0.5, self-consistency ratio > 0.5

### Broad histone (H3K27me3)
- [ ] MACS3 broadPeak output (not narrowPeak)
- [ ] `pooled_qc_summary.tsv` reports inferred_histone_class=broad_like, peak_mode_status=ok
- [ ] No IDR outputs (broad marks not eligible in v0.2)

### CUT&Tag
- [ ] `cuttag_fragment_size.tsv` reports sub-200 bp enrichment
- [ ] SEACR outputs under `04_peaks_seacr/` (when `cuttag.seacr.enabled: true`)
- [ ] Low duplication rate (<20%)

### ATAC
- [ ] MACS3 Tn5-aware parameters visible in logs
- [ ] `peak_counts.tsv` contains non-zero peak counts
- [ ] TSS profile shows enrichment peak at TSS (when `qc.tss_enrichment: true`)

## Known Runtime / Storage Constraints

- FASTQ files: typically 1-10 GB per sample (compressed). Total per queue: ~5-50 GB.
- Bowtie2 index: ~3-4 GB per genome assembly.
- Output BAM/BW/BDG: ~2-10 GB per sample.
- Total storage per queue: ~20-100 GB (FASTQ + index + outputs).
- Runtime: 2-8 hours per queue on a 16-core workstation, depending on read count.
- All intermediate and output files stay on local disk — nothing uploaded.

## Explicit "No Data Committed to Repo" Policy

- **No FASTQ files** (`.fq`, `.fq.gz`, `.fastq`, `.fastq.gz`) are committed.
- **No BAM/BAI files** are committed.
- **No BigWig/BedGraph files** (`.bw`, `.bdg`) are committed.
- **No MultiQC HTML reports** are committed.
- **No Bowtie2 index files** (`.bt2`, `.rev.*`) are committed.
- Only the following are committed:
  - This plan document.
  - Sample config/sample-sheet templates (with placeholder paths, not real data).
  - QC summary audits (small TSV files recording metrics, not raw data).
  - Release check reports.
- Use [`scripts/prepare_public_validation_inputs.py`](../../scripts/prepare_public_validation_inputs.py) to print the planned dataset inventory (no download).

## Execution Tiers

| Tier | Description | Trigger | Data Size |
| :--- | :--- | :--- | :--- |
| Smoke | Dry-run DAG check only | `snakemake -n` | 0 (no data needed except placeholder FASTQs) |
| Downsampled | Small subset (e.g., 1M reads) to verify file creation and QC plumbing | Manual | ~10 MB FASTQ |
| Functional | Moderate downsample to verify QC metrics are reasonable | Manual | ~500 MB - 2 GB FASTQ |
| Reference | Full dataset run for release record | Manual (tagged release) | Full FASTQ |
