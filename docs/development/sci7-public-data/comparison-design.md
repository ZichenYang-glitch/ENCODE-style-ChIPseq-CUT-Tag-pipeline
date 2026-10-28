# SCI-7 single-variable comparison design

## Fixed inputs
Use the same treatment and matched-control FASTQ files, proposed reference
identity, aligner, MAPQ, blacklist mode, peak thresholds, and tool versions for
every arm of a comparison. Record the effective config and tool versions.

## Comparisons
1. Duplicate policy: none/retain vs markdup/exclude vs auto, holding the
   analysis BAM independent variable to the duplicate policy only. Compare
   template counts, duplicate fraction, NRF/PBC, FRiP, peak count/width, and
   peak-set overlap.
2. MACS3 keep-dup: same analysis BAM, run MACS3 as a separate controlled
   experiment with `--keep-dup all` vs the alternative under a managed
   parameter contract; do not bypass the production validator or use
   free-form extra_args. Compare peaks against the same input/control.
3. PE max_insert: same PE library, compare 500 vs the assay default
   (ChIP/ATAC 2000, CUT&Tag 700, MNase 500), holding all other aligner
   parameters fixed. Record newly paired templates, coordinate changes,
   stability of high-confidence short fragments, duplicate fraction, peaks,
   and wall time/RAM.

## Data split
Discovery: ENCSR817LUF (histone), ENCSR844TIU (ATAC), ENCSR054ZMK (TF).
Confirmation must use different accessions and must not be used to tune the
decision rule.

## Decision rules
Freeze a default only if the effect is consistent in the confirmation set,
does not reduce replicate concordance, and does not create a clear
signal/background regression. Peak count alone is not a decision criterion.
Record effect sizes and uncertainty; state small-sample limitations.

## No-conclusion conditions
Mismatched controls, unverified reference identity, tool-version drift,
insufficient depth, SE/PE mismatch, phantom-dominated cross-correlation, or
discovery/confirmation overlap make the comparison inconclusive.
