# QC Interpretation Guide

This document describes each QC module produced by the pipeline and gives
practical guidance for interpreting the output. QC is **descriptive**, not a
universal pass/fail gate. Metrics depend on assay type, antibody quality,
library complexity, sequencing depth, and biological context.

## FastQC / Trim Galore

**Outputs:** `fastqc/` directory (HTML + zip per FASTQ), `*_trimming_report.txt`

FastQC reports per-base quality, GC content, adapter contamination, k-mer
over-representation, and sequence duplication levels. Trim Galore reports
record adapter trimming and quality filtering statistics.

**Interpretation:**
- **Per-base quality:** trailing bases with Q < 20-25 are expected and handled
  by Trim Galore. Persistent low quality in the first 10-15 bp may indicate
  library construction issues.
- **Adapter content:** after trimming, adapter content should be negligible
  (< 1%). Elevated adapter content post-trimming suggests incomplete adapter
  removal.
- **GC content:** ChIP-seq libraries typically show a narrow GC distribution
  centered around 40-50%. Broad or bimodal distributions may indicate
  contamination or PCR bias.
- **Sequence duplication:** high duplication levels are expected for ChIP-seq,
  particularly for TFs with sharp binding. This is normal — the duplicate
  metrics module provides more quantitative assessment.

## Alignment Metrics

**Outputs:** Bowtie2 log, `*.flagstat.txt`, `*.final.flagstat.txt`, `*.idxstats.txt`

**Interpretation:**
- **Overall alignment rate:** TF ChIP-seq: 60-90% is typical for human/mouse.
  Histone marks often exceed 80%. CUT&Tag typically exceeds 80%. Rates below
  50% warrant investigation of reference genome mismatch, contamination, or
  low library quality.
- **MAPQ filtering:** the pipeline filters to MAPQ >= 30 by default
  (configurable via `mapq`). Final flagstat shows post-filter read counts.
- **idxstats:** per-chromosome mapping. Expect uniform coverage across
  chromosomes proportional to chromosome size. Spikes in mitochondrial DNA
  (chrM) or unplaced contigs may occur and are generally not problematic for
  peak calling.
- **PE concordance:** for PE libraries, check the "properly paired" rate in
  flagstat. Low properly-paired rates (< 50%) suggest insert size issues or
  structural variation.

## Duplicate and Library Complexity

**Outputs:** `*.dup_metrics.txt`, `*.library_complexity.tsv`, `*.nrf_pbc.tsv`,
`*.preseq.txt`

### Duplicate metrics (Picard MarkDuplicates)

Recorded for all treatment samples. The pipeline may use samtools or Picard
for duplicate marking, depending on the runtime environment configured.

**Key fields in `library_complexity.tsv`:**
- `percent_duplication`: fraction of reads marked as duplicates. TF ChIP-seq:
  10-50% is common. Histone ChIP-seq: often < 20%. CUT&Tag: typically
  < 20%. High duplication (> 80%) may indicate over-amplification or
  insufficient starting material.
- `estimated_library_size`: estimated unique molecules. Compare with sequencing
  depth — if library size is much smaller than total reads, the library is
  saturated and additional sequencing would yield few new fragments.

### NRF / PBC (library complexity from read positions)

**Key fields in `nrf_pbc.tsv`:**
- **NRF** (Non-Redundant Fraction): number of distinct genomic positions with
  at least one mapped read divided by total mapped reads. Higher is better.
  ENCODE guidelines: NRF > 0.8 is good; < 0.5 indicates poor complexity.
- **PBC1** (PCR Bottleneck Coefficient 1): N1 / N_distinct, where N1 = number
  of positions with exactly one read. Closer to 1 is better. PBC1 > 0.8 is
  good; < 0.5 needs attention.
- **PBC2** (PCR Bottleneck Coefficient 2): N1 / N2, where N2 = number of
  positions with exactly two reads. PBC2 > 1 is the ENCODE-recommended
  threshold for TF ChIP-seq.

### Preseq (library complexity extrapolation)

Predicts the number of unique molecules expected at deeper sequencing depths.
The extrapolation curve shows whether additional sequencing would yield
diminishing returns.

**Interpretation:**
- A curve that plateaus quickly suggests the library is near saturation.
- A curve that continues to rise indicates additional sequencing would yield
  new unique fragments.
- Complements NRF/PBC: NRF/PBC describe current complexity; preseq predicts
  future yield.

## Peak Quality

**Outputs:** `*.peak_counts.tsv`, `*.frip.tsv`, `*_peaks.narrowPeak` or
`*_peaks.broadPeak`

### Peak Count

**Interpretation:**
- **TF ChIP-seq (narrowPeak):** 1,000-100,000 peaks depending on the factor.
  CTCTF typically yields 30,000-80,000 peaks. Abnormally low (< 500) or high
  (> 200,000) counts may indicate poor IP, inappropriate q-value threshold, or
  reference genome issues.
- **Histone ChIP-seq (broadPeak):** 10,000-100,000 peaks. Broad marks
  (H3K27me3, H3K9me3) span large domains and produce fewer but wider peaks.
- **CUT&Tag:** peak counts are generally comparable to ChIP-seq for the same
  target but may be affected by Tn5 insertion bias and lower background.

### FRiP (Fraction of Reads in Peaks)

Ratio of reads falling within called peaks to total mapped reads.

**Interpretation:**
- **TF ChIP-seq:** ENCODE guidelines consider FRiP > 0.01 (1%) acceptable and
  FRiP > 0.05 (5%) good. High-quality TF ChIP-seq often reaches 10-30%.
- **Histone ChIP-seq:** FRiP is typically much higher (20-50% for H3K27ac,
  > 50% for H3K4me3 broad domains) because histone marks cover more of the
  genome.
- **CUT&Tag:** FRiP is generally higher than ChIP-seq for the same target
  due to lower background.
- Low FRiP (< 1% for TF, < 5% for histone) may indicate poor IP efficiency,
  high background, or inappropriate peak-calling parameters.

### MACS3 Signal Tracks

**Outputs:** `*.FE.bdg` (fold-enrichment), `*.ppois.bdg` (Poisson p-value);
`*.FE.bw`, `*.ppois.bw` when `genome_resources.<genome>.chrom_sizes` is configured

FE tracks show enrichment over local background. Ppois tracks show the
statistical significance of enrichment. Both are bedGraph format (BigWig
conversion is planned). Useful for visualizing signal distribution across the
genome in a browser (IGV, UCSC).

## Cross-Correlation (phantompeakqualtools)

**Outputs:** `*.cc.qc`, `*.cc.plot.pdf`, `cross_correlation_summary.tsv`

The cross-correlation analysis is the ENCODE-standard approach for assessing
ChIP-seq library quality. It computes the Pearson correlation between
positive- and negative-strand reads at varying shift distances.

**Key metrics in `cross_correlation_summary.tsv`:**

- **NSC** (Normalized Strand Coefficient): ratio of the cross-correlation
  value at the estimated fragment length to the background cross-correlation.
  ENCODE guidelines: NSC > 1.05 is acceptable; NSC > 1.1 is good. Higher
  values indicate stronger ChIP enrichment relative to background.
- **RSC** (Relative Strand Correlation): ratio of the fragment-length
  cross-correlation to the phantom-peak (read-length) cross-correlation.
  ENCODE guidelines: RSC > 0.8 is acceptable. Higher values indicate that
  the enrichment signal dominates over the read-length artifact.
- **Estimated fragment length:** the shift distance at which the
  cross-correlation peaks. Should be consistent with expected library insert
  size (typically 100-300 bp for standard ChIP-seq).
- **Phantom peak:** the correlation peak at the read length. A prominent
  phantom peak is expected and reflects the strand-separation of paired reads
  at the read length.

**Quality flags** (auto-assigned):
- `ok`: NSC >= 1.05 and RSC >= 0.8
- `low_nsc`: NSC < 1.05 (weak signal-to-noise)
- `low_rsc`: RSC < 0.8 (phantom peak dominating)
- `low_nsc_low_rsc`: both thresholds failed
- `parse_failed`: `.cc.qc` file could not be parsed

**Important caveats:**
- These thresholds are ENCODE guidelines for TF ChIP-seq, not universal
  pass/fail gates. Histone marks and CUT&Tag have different cross-correlation
  profiles.
- CUT&Tag cross-correlation is less canonical than ChIP-seq. Tn5-based
  fragmentation produces different strand distributions.
- A sample flagged `low_nsc` may still yield biologically meaningful peaks;
  use it as a prioritization signal, not a rejection criterion.

The project-level summary at `results/multiqc/cross_correlation_summary.tsv`
aggregates all samples for easy comparison. When MultiQC is enabled, the same
summary is included in `multiqc_report.html` as a custom
**Cross-correlation QC** section.

## Picard CollectMultipleMetrics

**Outputs:** `*alignment_summary_metrics`, `*insert_size_metrics`,
`*quality_distribution_metrics`, `*insert_size_histogram.pdf`

Requires `genome_resources.<genome>.reference_fasta` with matching `.fai` and
`.dict`. Uses `VALIDATION_STRINGENCY=LENIENT` because filtered PE BAMs may
contain mate-field validation warnings.

**Key metrics:**
- **Alignment summary:** total reads, PF (pass-filter) reads, and alignment
  rate. Cross-reference with Bowtie2 alignment log.
- **Insert size metrics:** median, mean, and standard deviation of library
  insert sizes. PE libraries should show a unimodal distribution. Bimodal
  distributions may indicate adapter dimers, large deletions, or structural
  variation.
- **Quality distribution:** per-base quality score profile by cycle. Useful
  for identifying systematic quality degradation.

## TSS Enrichment-Style Profile

**Outputs:** `*.tss_matrix.gz`, `*.tss_profile.tsv`, `*.tss_profile.pdf`,
and `results/reference/<genome>.tss.bed`

When `qc.tss_enrichment: true`, the workflow derives transcription start sites
from `genome_resources.<genome>.gtf`, computes a deepTools matrix over
`-3000/+3000 bp` around each TSS using the CPM BigWig, and writes a profile
plot and tabular profile.

**Interpretation:**
- A clear signal peak around TSSs is expected for open-chromatin assays
  (ATAC-seq) and some promoter-associated histone marks.
- TF and broad histone ChIP-seq targets may not show a strong aggregate TSS
  enrichment profile; absence of a TSS peak is not automatically a QC failure.
- The GTF, FASTA, and Bowtie2 index must be from the same assembly. Annotation
  mismatch can flatten or shift the profile.

## CUT&Tag-Specific QC

### Fragment-Size Distribution

**Output:** `*.cuttag_fragment_size.tsv`

CUT&Tag uses Tn5 transposase which fragments DNA at nucleosome-adjacent sites.
Fragment-size distributions reflect the underlying chromatin structure:

- **Nucleosome-free regions (sub-200 bp):** enrichment of fragments < 200 bp
  suggests TF or other DNA-binding protein targets.
- **Mono-nucleosome (~200-300 bp):** typical for histone modification targets.
- **Multi-nucleosome (> 400 bp):** longer fragments indicate higher-order
  chromatin domains.
- Sub-150 bp enrichment is a hallmark of successful CUT&Tag.

### SEACR Peak Calling

**Output:** `04_peaks_seacr/<sample>/<sample>.seacr.<mode>.bed`

SEACR is a peak caller designed for low-background CUT&Tag data. It uses the
global background distribution rather than a local background model. The
pipeline runs SEACR as a sidecar alongside MACS3.

**Interpretation:**
- SEACR stringent peaks are typically fewer and more conservative than MACS3.
  Use stringent for high-confidence calls.
- SEACR relaxed peaks include lower-signal regions. Use relaxed when
  sensitivity is prioritized.
- Compare SEACR and MACS3 peak sets to assess peak-caller agreement.

## Replicate-Level Outputs

**Outputs:** `experiments/<exp>/` with pooled BAMs, pooled peaks, pooled
signal tracks, and `pooled_qc_summary.tsv`

When `stage4b: true` and the sample sheet defines experiments with >= 2
biological replicates:

- **Pooled BAM:** all treatment biorep BAMs merged. Deeper effective coverage
  for peak calling.
- **Pooled peaks:** MACS3 called on the pooled BAM. Should recover more peaks
  than individual replicates.
- **Pooled QC summary:** histone target classification, peak-mode
  compatibility check, signal track status.

### TF ChIP-seq IDR

**Outputs:** `experiments/<exp>/06_idr/` with true-replicate IDR,
pseudoreplicate IDR, and final conservative/optimal peak sets

IDR (Irreproducible Discovery Rate) quantifies the consistency of peak calls
between biological replicates. The final peak sets are:

- **Conservative:** high-confidence peaks reproducible across both replicates
  and pseudoreplicates.
- **Optimal:** balances sensitivity and specificity.

The `reproducibility_summary.tsv` records N_peaks, rescue ratio, and
self-consistency ratio for each IDR stage.

## Practical Checklists

### ChIP-seq QC Checklist

1. **Alignment rate** > 60% for human/mouse.
2. **FRiP** > 1% (TF) or > 10% (histone); ideally > 5% / > 20%.
3. **NSC** > 1.05; **RSC** > 0.8 (ENCODE recommendation, not a hard gate).
4. **NRF** > 0.5; **PBC1** > 0.5; **PBC2** > 1.
5. **Peak count** within expected range for the target (1,000-100,000 for TF).
6. **Fragment length** estimate consistent with library insert size.
7. **Insert size distribution** unimodal (PE libraries).
8. **Cross-contamination:** idxstats shows expected chromosomal distribution.
9. **Preseq** curve does not plateau immediately (library not saturated).
10. **Adapters** removed (FastQC + Trim Galore).

### CUT&Tag QC Checklist

1. **Alignment rate** > 80%.
2. **Fragment size** enrichment below 200 bp (sub-nucleosomal).
3. **FRiP** higher than equivalent ChIP-seq target.
4. **SEACR peaks** recovered; compare stringent vs relaxed.
5. **Cross-correlation:** interpret cautiously — CUT&Tag strand distributions
   differ from ChIP-seq. Low NSC may not indicate poor quality.
6. **Low duplication rate** (< 20%) — Tn5 does not amplify, so true PCR
   duplicates should be minimal.
7. **NRF/PBC** still informative but thresholds are less strict than ChIP-seq.
8. **MultiQC** shows expected per-base quality and adapter content.

### ATAC-seq QC Checklist

1. **Alignment rate** high and chromosome distribution expected.
2. **Insert size** shows nucleosome-free and nucleosomal periodicity when PE
   data are available.
3. **TSS profile** has a clear enrichment peak when `qc.tss_enrichment: true`.
4. **FRiP** and peak counts are consistent across biological replicates.
5. **Preseq/NRF/PBC** suggest the library is not exhausted.
6. **Peak mode** is `narrow`; broad ATAC peak calling is not supported.

## Known Caveats

- **All thresholds are context-dependent.** Antibody quality, cell type,
  sequencing depth, and biological system affect every metric. Use thresholds
  as guidelines, not gates.
- **CUT&Tag cross-correlation** is less canonical than ChIP-seq. The
  Tn5-based transposition mechanism produces different strand distributions.
  Low NSC values are common for CUT&Tag even in high-quality libraries.
- **Histone marks** (broad domains) and **TFs** (sharp peaks) behave
  differently. Broad marks naturally have higher FRiP and different
  cross-correlation profiles. The pipeline's pooled QC summary
  classifies targets to flag potential mismatches.
- **No hard-fail thresholds.** The pipeline does not abort on low QC metrics.
  Review the QC summary tables and MultiQC report to make informed judgments.
- **FE/ppois bedGraph and BigWig tracks** are generated. BigWig conversion
  requires configuring `genome_resources.<genome>.chrom_sizes`. When
  `chrom_sizes` is empty, only bedGraph files are produced. Visualize
  bedGraph files directly in IGV or convert with UCSC `bedGraphToBigWig`.
- **Reference genome** must match the alignment index. Mismatched genomes
  produce systematically low alignment and meaningless QC values.
