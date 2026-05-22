# Stage 13 Real-Data Validation

**Date:** 2026-05-22
**Status:** PASS with follow-ups
**Scope:** 2-replicate CUT&Tag PE run + 2-replicate ChIP-seq PE run (data
external to repo)

## CUT&Tag Real Run

| Field | Value |
| :--- | :--- |
| Genome | mm39 / GRCm39 |
| Experiment / target | cuttag_SMARCA5 |
| Layout | PE, 2 biological replicates |
| Control | none |
| SEACR | enabled (stringent) |

### Outputs validated

- `final.bam` / `final.bam.bai`
- CPM BigWig
- MACS3 narrowPeak
- MACS3 FE bedGraph
- MACS3 ppois bedGraph
- SEACR stringent BED
- FRiP, NRF/PBC, library_complexity
- Per-sample qc_summary
- CUT&Tag fragment size QC
- cross_correlation (`.cc.qc` + `.cc.plot.pdf`)
- preseq (`.preseq.txt`)
- Picard CollectMultipleMetrics outputs
- Pooled replicate BAM / peaks / signal / QC
- MultiQC report (captured Picard + preseq)

## ChIP-seq Real Run

| Field | Value |
| :--- | :--- |
| Genome | mm39 / GRCm39 |
| Experiment / target | chipseq_H3K27ac |
| Layout | PE, 2 biological replicates |
| Control | none |

### Outputs validated

- `final.bam` / `final.bam.bai`
- CPM BigWig
- MACS3 narrowPeak
- MACS3 FE bedGraph
- MACS3 ppois bedGraph
- FRiP, NRF/PBC, library_complexity
- Per-sample qc_summary
- cross_correlation (`.cc.qc` + `.cc.plot.pdf`)
- preseq (`.preseq.txt`)
- Picard CollectMultipleMetrics outputs
- Pooled replicate BAM / peaks / signal / QC
- MultiQC report (captured Picard + preseq)

## Stage 12 QC Validation

- **cross_correlation:** phantompeakqualtools `.cc.qc` and `.cc.plot.pdf`
  generated for all treatment samples on both assays.
- **preseq:** `lc_extrap -B` outputs generated and visible in MultiQC.
- **Picard metrics:** CollectMultipleMetrics outputs generated and visible in
  MultiQC (alignment summary, insert size, quality distribution, plus side
  outputs: base distribution by cycle, quality by cycle, quality distribution
  PDF, read length histogram).

## Issues Found and Fixed

### 1. Missing `.dict` for Picard reference

Picard CollectMultipleMetrics requires `reference_fasta` plus a matching
samtools FASTA index (`.fai`) and a Picard sequence dictionary (`.dict`) next
to the FASTA. The initial run failed with a Picard error when `.dict` was
absent.

**Resolution:** Documented the `.fai` + `.dict` requirement in
`docs/configuration.md` and `README.md`. Users must create the dictionary with
`picard CreateSequenceDictionary R=genome.fa O=genome.dict` (or `samtools dict`).

### 2. Picard STRICT validation failed on INVALID_FLAG_MATE_UNMAPPED

Real PE BAMs after MAPQ and flag filtering can contain reads whose mapped mate
has no mate reference name set, triggering Picard's default
`VALIDATION_STRINGENCY=STRICT` to abort.

**Resolution:** Added `VALIDATION_STRINGENCY=LENIENT` to the
`picard_collect_multiple_metrics` rule in `workflow/rules/qc.smk`. These are
QC-only metrics; lenient validation allows metric collection to proceed on
filtered BAMs without affecting the underlying data.

## Follow-ups

1. **MultiQC visibility for phantompeakqualtools `cc.qc`:** cross-correlation
   `.cc.qc` files were generated but MultiQC visibility was not confirmed.
   Verify whether MultiQC recognizes these files, or add a custom parser
   configuration if needed.
2. **PE mate-field cleanup:** evaluate `samtools fixmate` after MAPQ filtering
   if `INVALID_FLAG_MATE_UNMAPPED` warnings cause issues for downstream tools
   beyond Picard.
3. **No real data bundled:** real FASTQ, BAM, BAI, BW, BDG, HTML, and MultiQC
   outputs remain external to the repository and are not committed.

## Code Changes (this stage)

| File | Change |
| :--- | :--- |
| `workflow/rules/qc.smk` | Added `VALIDATION_STRINGENCY=LENIENT` to `picard_collect_multiple_metrics` |
| `docs/configuration.md` | Documented `.fai`/`.dict` requirement and LENIENT rationale |
| `README.md` | Added note on `.fai`/`.dict` and LENIENT for Picard metrics |
| `KNOWN_ISSUES.md` | Added Stage 13 real-data validation summary |
| `docs/release-checks/stage13-real-data-validation.md` | This file |
