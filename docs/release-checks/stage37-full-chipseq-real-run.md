# Stage 37: Full Histone ChIP-seq Real-Run Report

**Date:** 2026-05-26 to 2026-05-27
**Status:** completed (30-sample real ChIP-seq run)
**Scope:** external paired-end histone ChIP-seq dataset; no FASTQ, BAM, BigWig,
bedGraph, HTML, or other generated data committed to this repository

---

## 1. Dataset and Configuration

| Field | Value |
| :--- | :--- |
| Project | Regulation of T helper cell differentiation ChIP-seq |
| Assay | ChIP-seq |
| Genome label | `mm10` |
| Reference build | GRCm38 Bowtie2 index and genome resources |
| Layout | PE |
| Samples | 30 treatment libraries |
| Design | 5 conditions x 3 histone marks x 2 biological replicates |
| Conditions | Naive, Th1 24h, Th1 72h, Th2 24h, Th2 72h |
| Targets | H3K27ac, H3K4me1, H3K4me3 |
| Controls | none (`use_control: false`) |
| MAPQ | 30 |
| Trimming | enabled |
| Duplicate handling | `remove_dup: auto` |
| Blacklist filtering | enabled |
| Signal tracks | enabled |
| Cross-correlation | enabled |
| Preseq / Picard metrics | enabled |
| TSS enrichment | disabled |
| Stage 4b pooled outputs | enabled |
| Stage 5 IDR | disabled |

Peak-mode policy matched the histone targets:

- H3K27ac: `narrow`
- H3K4me1: `broad`
- H3K4me3: `narrow`

Replicate encoding used two biological replicates per experiment and one
technical replicate per library.

---

## 2. Completion Summary

The full Snakemake run completed for all samples:

| Check | Result |
| :--- | :--- |
| Sample-level `pipeline.done` sentinels | 30 / 30 |
| Project MultiQC report | present |
| Project QC summary | present |
| Project result manifest | present |
| Project cross-correlation summary | present |
| Manifest records marked `present` | 452 |

The run validated the current real-data contract for a no-control histone
ChIP-seq project with sample-level and pooled experiment-level outputs.

---

## 3. Outputs Validated

Sample-level outputs:

- 30 final BAMs and BAIs
- 30 CPM BigWigs
- 30 MACS3 peak directories
- 30 blacklist-filtered peak sets
- 30 FE BigWigs and 30 ppois BigWigs
- 30 FE bedGraphs and 30 ppois bedGraphs
- 30 per-sample QC summaries
- FastQC, Trim Galore, Bowtie2, samtools flagstat/idxstats, Picard, preseq,
  NRF/PBC, FRiP, and cross-correlation outputs

Pooled experiment-level outputs:

- 15 pooled final BAMs and BAIs
- 15 pooled MACS3 peak directories
- 15 pooled QC summaries
- 15 pooled FE BigWigs and 15 pooled ppois BigWigs
- 15 pooled FE bedGraphs and 15 pooled ppois bedGraphs
- 30 biological-replicate BAM symlinks or merged BAMs under the experiment
  output tree

Post-run exploratory outputs were also generated externally for 1 kb global
BigWig correlation and per-target PCA plots. These are analysis artifacts, not
workflow source files.

---

## 4. QC Summary

### 4.1 Bowtie2 alignment

| Target | Mean overall alignment rate | Range |
| :--- | ---: | ---: |
| H3K27ac | 97.8% | 95.8-99.0% |
| H3K4me1 | 97.8% | 95.1-99.0% |
| H3K4me3 | 88.8% | 25.9-98.9% |

One H3K4me3 sample had a low Bowtie2 overall alignment rate:

- `ChIP-seq_MLL4_WT_Naive_H3K4me3_rep1`: 25.86%

Despite this alignment outlier, downstream H3K4me3 signal remained coherent:
FRiP was 0.710, cross-correlation passed, and replicate Pearson correlation
against `Naive_H3K4me3_rep2` was 0.9498.

### 4.2 FRiP and peak counts

| Target | Mean reads after processing | FRiP mean | FRiP range | Mean peak count | Peak-count range |
| :--- | ---: | ---: | ---: | ---: | ---: |
| H3K27ac | 47.5M | 0.132 | 0.085-0.189 | 54,774 | 33,726-70,310 |
| H3K4me1 | 56.2M | 0.360 | 0.292-0.428 | 122,786 | 94,643-140,860 |
| H3K4me3 | 35.5M | 0.658 | 0.554-0.773 | 51,596 | 36,660-62,256 |

Lowest FRiP sample:

- `ChIP-seq_MLL4_WT_Th1_72h_H3K27ac_rep2`: FRiP 0.0847, 33,726 peaks

### 4.3 Replicate signal correlation

Replicate correlation was evaluated from 1 kb global CPM BigWig signal.

| Target | Mean rep1-vs-rep2 Pearson | Range |
| :--- | ---: | ---: |
| H3K27ac | 0.9590 | 0.9408-0.9748 |
| H3K4me1 | 0.9630 | 0.9430-0.9820 |
| H3K4me3 | 0.9862 | 0.9498-0.9976 |

The lowest Pearson pair was still 0.9408
(`Th1_72h_H3K27ac`), supporting good replicate-level signal consistency across
the run.

### 4.4 Cross-correlation flags

Cross-correlation status by target:

| Target | Summary |
| :--- | :--- |
| H3K27ac | 1 ok, 3 low RSC, 6 low NSC + low RSC |
| H3K4me1 | 1 low RSC, 9 low NSC + low RSC |
| H3K4me3 | 10 ok |

The broad H3K4me1 profile and some H3K27ac samples showed low NSC/RSC, but this
did not correspond to poor replicate BigWig correlation. H3K4me3
cross-correlation was consistently strong.

---

## 5. Issues Observed

### 5.1 bedGraph-to-BigWig disk pressure

The signal-track rules generate large FE/ppois bedGraphs and also sort temporary
bedGraphs during BigWig conversion. High parallelism can consume hundreds of GB
of scratch/output space and caused an intermediate failure when disk space was
nearly exhausted.

Follow-up candidates:

- Add a configurable concurrency/resource limit for bedGraph-to-BigWig rules.
- Route sort temporary files to an explicit scratch directory.
- Improve sort stderr logging so disk-full failures are obvious.
- Consider optional cleanup or `temp()` handling for large bedGraph
  intermediates when users only need BigWigs.

### 5.2 MultiQC report filename with title

MultiQC generated a title-derived filename in addition to the expected
`multiqc_report.html`. The run was completed with a `--filename
multiqc_report.html` workaround supplied via `tool_parameters.multiqc.extra_args`.

Follow-up candidate:

- Make `--filename multiqc_report.html` part of the workflow rule default, or
  document the required extra arg when a project title is supplied.

### 5.3 MultiQC sample-name normalization

The report showed multiple naming styles for the same biological library
because different tools emitted FASTQ accessions, mate suffixes, sample names,
and `.final` BAM names. This is a display/normalization issue; the underlying
outputs and manifest were complete.

Follow-up candidate:

- Add project-independent MultiQC sample-name cleanup while preserving enough
  stage identity to avoid accidentally overwriting pre-filter and final BAM QC
  metrics.

---

## 6. Artifact Policy

- The real input data and all generated outputs remain external to the
  repository.
- No FASTQ, BAM, BAI, BigWig, bedGraph, peak file, MultiQC HTML, image, or
  Snakemake metadata artifact is committed.
- This report records validation scope and observed QC summaries only.

---

## 7. Conclusion

The workflow successfully completed a 30-library real histone ChIP-seq run with
sample-level QC, peak calling, signal tracks, pooled replicate outputs, project
MultiQC, and a result manifest. The run supports the current no-control histone
ChIP-seq execution path and highlights follow-ups for disk-heavy signal-track
conversion and MultiQC sample-name/report-name ergonomics.
