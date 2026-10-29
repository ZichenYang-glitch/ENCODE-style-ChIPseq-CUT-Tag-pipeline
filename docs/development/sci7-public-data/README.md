# SCI-7 public-data validation preparation

Status: one CUT&Tag Rep2 retain-versus-exclude comparison is recorded in
[`cuttag-duplicate-retain-exclude.md`](cuttag-duplicate-retain-exclude.md).
It is a single treatment/control effect-size check. No scientific acceptance
has completed, and no default is frozen. The other planned assay comparisons
have not run. That comparison did not download another FASTQ or reference.

## Scope

The planned comparisons cover ChIP-seq TF, ChIP-seq histone, CUT&Tag, ATAC-seq
and MNase-seq, with PE and SE where supported and matched controls where the
assay design uses them. Duplicate policy, MACS3 keep-dup and PE `max_insert` are
separate single-variable experiments. The executable matrix, effect sizes,
decision rules and no-conclusion conditions are pre-registered in
[`comparison-design.md`](comparison-design.md).

## Verified discovery file metadata

[`file-manifest.tsv`](file-manifest.tsv) records experiment, assay, target,
species, proposed reference, layout, biological/technical replicate,
library/biosample/donor, treatment background, control pairing, file accession,
download URL, official byte size, official MD5, read length, mate pairing,
source query, query date and verification status.

The first 15 rows are marked `official`: their file accessions, sizes, MD5
values and pairing were read from the saved official ENCODE API responses in
[`metadata/`](metadata/). The 14 GEO/SRA rows are marked `official metadata; URL
HEAD verified`: their values come from ENA project filereports and every HTTPS
HEAD Content-Length matched. These statuses describe metadata/link provenance,
not local file possession, downloaded-content checksums, reference identity, or
biological suitability. `reference_proposed` remains a proposed execution
reference; raw FASTQ records do not themselves establish an official assembly.

| Discovery group | Selected experiment(s) | Layout/replicate facts | Official FASTQ bytes |
| --- | --- | --- | ---: |
| TF ChIP-seq | ENCSR054ZMK (NR2F2, HepG2, FOXA3 RNAi background) + control ENCSR142JFU | 100 bp PE; one biological replicate | 11,396,944,028 |
| Histone ChIP-seq | ENCSR817LUF (H3K36me3) + control ENCSR603VTC | 36 bp SE; treatment files belong to one library and one technical replicate, not five independent replicates | 2,510,695,554 |
| ATAC-seq | ENCSR844TIU | 26 bp PE; two biological replicates; short-read applicability remains an explicit limitation | 3,952,560,652 |
| CUT&Tag | GSE159608 / PRJNA669859 H3K27me3 Rep1/Rep2 + matched IgG Rep1/Rep2 | 2x25 bp PE; submitter-labelled replicates; Rep1 is shallow | 443,578,936 |
| MNase-seq | PRJNA928451 LP_pMSCs_FM Rep1/Rep2/Rep3 | PE, 2x51 bp inferred from official counts; SRA explicitly states three independent biological replicates | 4,402,061,533 |
| **Total** | **29 FASTQ files** |  | **22,705,840,703** |

CUT&Tag/MNase selection, rejected candidates, control-pairing caveats and all 46
official GEO/SRA/ENA/ENCODE queries are documented in
[`cuttag-mnase-selection.md`](cuttag-mnase-selection.md) and
[`metadata/geo-sra/`](metadata/geo-sra/). All 24 discovery plus optional
confirmation ENA links returned HTTP 200 with matching Content-Length. ENCODE
ENCSR000CXP/ENCSR000CXQ were rejected because they are SOLiD rather than
standard Illumina FASTQ; an earlier empty ENCODE API search was not generalized
into a claim that public databases lack these assays.

## Confirmation and reference gates

The ENCODE confirmation candidates TF ENCSR785RQR, histone ENCSR213YKJ and ATAC
ENCSR639PPN have now been checked for donor, biosample, library, source-file and
checksum overlap in [`confirmation-review.md`](confirmation-review.md). TF is
independent and uses background-matched control ENCSR194GTV; histone uses only
rep1 with control ENCSR148TRG; ATAC shares donor ENCDO236IKC with discovery and
therefore supports only library/cell-context confirmation. Optional GEO/SRA
CUT&Tag/MNase candidates are recorded separately with shared-control,
same-study/donor and depth limitations. Different experiment accessions alone
are never treated as proof of independence.

[`confirmation-file-manifest.tsv`](confirmation-file-manifest.tsv) contains 28
ENCODE confirmation FASTQ (92,835,234,501 bytes), and
[`geo-sra-confirmation-file-manifest.tsv`](geo-sra-confirmation-file-manifest.tsv)
contains 10 optional GEO/SRA FASTQ (3,371,454,387 bytes). None is downloaded.

[`reference-identity.md`](reference-identity.md) separately records proposed
FASTA assembly/release, contig naming, FAI/chrom_sizes derivation, blacklist
provenance/version/checksum, any required annotation, and Bowtie2 index build
tool/parameters/input identity. The two small ENCODE blacklist payloads were
downloaded and their official sizes/MD5 verified. Large FASTA/GTF/index payloads
remain explicitly unbuilt and unverified.

## Execution limits and budget

[`resource-budget.md`](resource-budget.md) separates official input sizes from
estimated storage, CPU/RAM and wall time, and defines the minimal trial and
expansion order. Downloads and real-library execution require a later explicit
execution step; this preparation commit is not scientific validation.

## Evidence storage

[`evidence-storage.md`](evidence-storage.md) and
[`evidence-integrity.tsv`](evidence-integrity.tsv) describe the deterministic,
lossless gzip representation used for 17 evidence files, their original-byte
SHA256 values, source queries and decompression verification procedure.
