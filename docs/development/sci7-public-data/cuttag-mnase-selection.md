# CUT&Tag and MNase-seq file-level selection

Status: official metadata and download-link review completed on 2026-10-28. No
FASTQ/BAM was downloaded, decompressed, aligned or scientifically evaluated.
All sizes/MD5 values are ENA Portal API metadata; no checksum is inferred.

## Selected discovery data

Two human Illumina paired-end datasets were selected. Together they contribute
**7 RUNs, 14 FASTQ files and 4,845,640,469 bytes (4.5129 GiB)**. Their complete
file rows are included in [`file-manifest.tsv`](file-manifest.tsv).

### CUT&Tag: GSE159608 / PRJNA669859

ALL-1 KMT2A-AF4 B-ALL patient-derived CD45+ cells, AutoCUT&Tag, H3K27me3
Rep1/Rep2 with same-patient/batch IgG Rep1/Rep2, Illumina HiSeq 2500 paired-end
2x25 bp. BioSamples are SAMN16477147, SAMN16477146, SAMN16477220 and
SAMN16477219; the source does not provide a finer anonymous donor identifier.

The proposed pairings are SRR12842271 -> SRR12842281 (Rep1 -> IgG Rep1) and
SRR12842272 -> SRR12842282 (Rep2 -> IgG Rep2), based on patient, replicate label
and experiment date. They are **proposed**, not explicit SRA/GEO `paired_with`
relationships. Rep1/Rep2 are submitter-labelled replicates on different dates;
GEO does not explicitly call them independent biological replicates.

The source analysis used UCSC hg19 and Bowtie2 `-X 700`. SCI-7 instead proposes
fresh alignment to the GRCh38 primary assembly in
[`reference-identity.md`](reference-identity.md); this is not a liftOver and is
not an assembly attribute of the FASTQ.

| GEO / SRX / RUN | Role | ENA read_count | Official bytes, both mates |
| --- | --- | ---: | ---: |
| GSM4835654 / SRX9309566 / SRR12842271 | H3K27me3 Rep1 | 941,245 | 38,171,016 |
| GSM4835655 / SRX9309567 / SRR12842272 | H3K27me3 Rep2 | 5,515,652 | 232,319,097 |
| GSM4835664 / SRX9309576 / SRR12842281 | IgG Rep1 | 534,078 | 21,774,402 |
| GSM4835665 / SRX9309577 / SRR12842282 | IgG Rep2 | 3,621,858 | 151,314,421 |
| **Total** | **4 RUNs / 8 FASTQ** |  | **443,578,936** |

Rep1 treatment and IgG are shallow. Alignment rate, valid fragment count,
FRiP/background and replicate consistency are admission checks; if underpowered,
Rep1 is a workflow trial only and the affected comparison is inconclusive. The
small download size is not evidence that these libraries can freeze a default.

### MNase-seq: PRJNA928451 LP_pMSCs_FM

Human umbilical-cord mesenchymal stem cells, NEBNext Ultra II, MNase digestion,
gel purification, Illumina NovaSeq 6000 paired-end. The SRA experiment design
explicitly states **“MNase-Seq analyses were performed on three independent
biological replicates.”** BioSamples SAMN32933787/88/89 are the three replicates.

SRA does not directly state per-mate read length. For every selected RUN,
`base_count / read_count = 102`, consistent with 2x51 bp only if mates have equal
length; actual FASTQ records must confirm this. The source design says hg38,
while SCI-7 proposes a specific GRCh38 primary-assembly FASTA. The narrative
instrument text says “HiSeq-4K” but structured metadata says NovaSeq 6000; the
latter is used and the discrepancy is retained.

| BioSample / SRX / RUN | Replicate | ENA read_count | base_count | Official bytes, both mates |
| --- | --- | ---: | ---: | ---: |
| SAMN32933787 / SRX19193323 / SRR23247654 | LP_pMSCs_FM Rep1 | 23,553,446 | 2,402,451,492 | 1,536,213,955 |
| SAMN32933788 / SRX19193324 / SRR23247653 | LP_pMSCs_FM Rep2 | 23,005,066 | 2,346,516,732 | 1,461,205,788 |
| SAMN32933789 / SRX19193325 / SRR23247651 | LP_pMSCs_FM Rep3 | 22,250,126 | 2,269,512,852 | 1,404,641,790 |
| **Total** | **3 RUNs / 6 FASTQ** |  |  | **4,402,061,533** |

MNase has no CUT&Tag-style IgG control. Each R1/R2 pair is a sequencing mate,
not an experimental control pairing.

## Optional confirmation data

The separate
[`geo-sra-confirmation-file-manifest.tsv`](geo-sra-confirmation-file-manifest.tsv)
contains 5 RUNs / 10 FASTQ / 3,371,454,387 bytes. These are optional candidates,
not an already accepted independent cohort.

- **CUT&Tag GSE145187 / PRJNA606273:** crosslinked K562 H3K4me3 Rep1/Rep2 plus
  shared IgG XLink, 6 FASTQ / 232,017,631 bytes. No RUN is shared with the
  ALL-1 discovery set and the cell source differs. Both treatments share
  SRR11923224 rather than having replicate-specific controls. Rep1 used 120,000
  cells versus 500,000 for Rep2, so these are not identical-condition
  replicates. The IgG title/antibody fields say IgG XLink while its description
  says H3K27me3; this source-metadata inconsistency must remain visible.
- **MNase PRJNA928451 HP_pMSCs_FM Rep1/Rep2:** 4 FASTQ / 3,139,436,756 bytes,
  with the same independent-biological-replicate design statement. HP and LP
  use different BioSamples/RUNs but belong to one study, and donor metadata are
  both male umbilical cord from the same provider. This supports same-study
  condition confirmation only; donor independence is unproven.

Neither optional confirmation set may be used to tune the decision rule. The
CUT&Tag candidate is also shallow and must pass the same depth/no-conclusion
checks.

## Rejected or superseded candidates

- ENCODE ENCSR000CXP/ENCSR000CXQ were re-queried and rejected: their raw data
  are SOLiD csfasta/csqual/BAM/bigWig on hg19, not compatible Illumina FASTQ.
- GSM4308159/SRR11074242 and GSM4308163/SRR11074246 are different CUT&Tag
  targets/conditions, not a treatment-control pair or same-condition replicate
  pair. GSE145187 provides a clearer Rep1/Rep2 + IgG candidate.
- PRJNA833071 Capture-MNase runs are probe/TSS-enriched, not whole-genome
  MNase; whole-genome runs in that project are substantially larger.
- PRJNA1173713/GSM8577139/SRR31015809 is TSS-capture, T2T-based and lacks a
  same-condition replicate, so it is not an independent MNase confirmation.
- GSE36979/GSM907783/SRR452459 has no same-condition biological replicate,
  uses hg18 and has only about 2x18 bp inferred reads.
- PRJNA379378 is single-cell scMNase-seq, not bulk MNase input for this
  workflow. PRJNA928451 XM runs are much larger and were superseded by the FM
  three-replicate set.

These are candidate-specific compatibility decisions. An empty or failed API
search was never generalized into a claim that ENCODE/GEO/SRA lacks an assay.

## Queries and link verification

[`metadata/geo-sra/queries.tsv`](metadata/geo-sra/queries.tsv) records all 46
official queries and their saved responses. Sources include GEO SOFT and GDS
ESearch/ESummary, ENA `read_run` searches, project filereports for
PRJNA669859/PRJNA606273/PRJNA928451, ENA experiment/sample XML, and ENCODE
experiment JSON. Rejected-candidate records are retained for audit.

The complete project filereports provide RUN/BioSample/layout/instrument,
read/base counts and per-file FTP path, MD5 and bytes. FTP host paths were
converted to equivalent HTTPS URLs only.

All **24/24** listed discovery and optional-confirmation FASTQ HTTPS URLs
returned HTTP 200 to HEAD, and every Content-Length exactly matched the official
ENA byte count. See
[`metadata/geo-sra/url_head_checks.tsv`](metadata/geo-sra/url_head_checks.tsv)
and [`metadata/geo-sra/selection_summary.json`](metadata/geo-sra/selection_summary.json).
HEAD validation does not verify content MD5, gzip integrity, read-name
synchrony, actual read length, quality encoding or biological QC; those remain
execution-time checks.
