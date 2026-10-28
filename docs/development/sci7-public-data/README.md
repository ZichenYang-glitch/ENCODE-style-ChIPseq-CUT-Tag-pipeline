# SCI-7 public-data preparation (no downloads)

Status: preparation only. No scientific acceptance, no default freeze.

## Scope (from the repo SCI-7 handoff)

Single-variable comparisons for the duplicate-policy, MACS3 keep-dup, and PE
max_insert decisions on public ENCODE-style data. Candidate groups:
ChIP-seq TF, ChIP-seq histone, CUT&Tag, ATAC-seq, MNase-seq; PE and SE where the
assay supports it; treatment + matched control / input.

## Candidate manifest fields

`/tmp/sci7-work/file-manifest.tsv` records:
accession, assay, target, species, reference assembly/version, layout,
replicate, control pairing, file accession/URL, official size, checksum
(md5/sha256), source, query date, verification status.

Every row is currently `unverified`; no URL/size/checksum is asserted until the
official portal response is captured. No FASTQ/BAM download has started.

## Single-variable design

1. Duplicate policy: same file set, only duplicate handling differs
   (none/retain vs markdup/exclude vs auto). Hold aligner, MAPQ, blacklist,
   peak caller, and thresholds fixed.
2. MACS3 keep-dup: same analysis BAM and control, only `--keep-dup` differs.
3. PE max_insert: same PE library, only the insert cap differs
   (e.g. 500/700/2000), with fragment/template counts and peak sets compared.

Each comparison records: input file identities, tool versions, effective
config, command, output checksums, and a pre-registered decision rule.

## Metrics

- Template/fragment counts by validity class; duplicate fraction; NRF/PBC.
- FRiP, peak count, peak width distribution, peak-set Jaccard/overlap.
- Cross-correlation NSC/RSC/phantom status where QC is enabled.
- Qualitative failure modes: SE/PE mismatch, missing control, insufficient
  depth, variable read length, phantom-dominated correlation.

## Data split

Discovery/verification cohorts are separate accessions. Any accession used to
select a default is excluded from the confirmation cohort. No tuning on the
confirmation cohort.

## Non-conclusions

A comparison is inconclusive if the input is underpowered, controls are
mismatched, the reference/assembly identity is unverified, or the tool version
differs from the pinned closure.


## Actual candidates selected 2026-10-28

- TF ChIP-seq: ENCSR054ZMK (NR2F2, human HepG2, PE) + matched control
  ENCSR142JFU; 4 FASTQ, 10.61 GiB.
- Histone ChIP-seq: ENCSR817LUF (H3K36me3, human brain, SE) + matched control
  ENCSR603VTC; 7 FASTQ, 2.34 GiB.
- ATAC-seq: ENCSR844TIU (mouse Treg, PE, 2 bio replicates); 4 FASTQ, 3.68 GiB.
- MNase-seq: ENCSR000CXP/ENCSR000CXQ are SOLiD csfasta/csqual/BAM hg19 with no
  standard FASTQ; recorded as a gap. GEO/SRA query remains.
- CUT&Tag: ENCODE API search did not return CUT&Tag experiments; recorded as a
  gap. GEO/SRA query remains.

File-level rows, official sizes/md5, source APIs and status are in
`file-manifest.tsv`; query evidence is in `source-log.md` and `metadata/`.

## Confirmation candidates (unverified independence)

- TF: ENCSR785RQR (JUND) — independent-accession candidate; donor/library overlap
  check against ENCSR054ZMK still pending.
- Histone: ENCSR213YKJ (H3K4me3) — independent-accession candidate; donor/library
  overlap check still pending.
- ATAC: ENCSR639PPN — independent-accession candidate; species/donor and
  library overlap check still pending.

These are not acceptance data until the donor, biosample, library, and file
hashes are checked against the discovery set.

## GEO/SRA gap queries (2026-10-28)

- NCBI GEO `CUT&Tag`: 19,003 GDS hits; no file-level candidate selected yet.
  Raw responses: `metadata/geo_CUT_Tag.json`.
- NCBI GEO `MNase-seq`: 3,885 GDS hits; no standard-FASTQ candidate selected yet.
  Raw responses: `metadata/geo_MNase-seq.json`.
- ENCODE CUT&Tag API search returned no experiment; ENCODE MNase experiments are
  SOLiD-only. These remain gaps; do not fabricate file identities.
