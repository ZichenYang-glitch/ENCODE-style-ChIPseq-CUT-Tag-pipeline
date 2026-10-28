# SCI-7 confirmation-cohort independence review

Status: official metadata review completed on 2026-10-28. No confirmation FASTQ
has been downloaded or executed. `official` sizes/MD5 values below are source
metadata, not local payload verification.

## Decision summary

| Group | Candidate decision | Independence and design limits |
| --- | --- | --- |
| TF | **Keep ENCSR785RQR (JUND), with control ENCSR194GTV** | No donor, biosample, library, source FASTQ accession or MD5 overlap with discovery ENCSR054ZMK/ENCSR142JFU. JUND/K562 CRISPR-CEBPB differs from NR2F2/HepG2 FOXA3-RNAi in target and background. One biological replicate only; no within-candidate replicate-concordance estimate. |
| Histone | **Keep ENCSR213YKJ, but only biological replicate 1, with control ENCSR148TRG** | No donor/biosample/library/FASTQ/MD5 overlap with discovery ENCSR817LUF/ENCSR603VTC. This is mouse B6CASTF1/J hippocampus H3K4me3, not a same-species or same-target replication of human brain H3K36me3. Rep1 has a same-biosample control; rep2 is excluded from the matched-control main manifest because it lacks a corresponding rep2 control. |
| ATAC | **Keep ENCSR639PPN only for library/cell-context confirmation** | No biosample, library, source FASTQ or MD5 overlap with discovery ENCSR844TIU, and four distinct biosamples/biological replicates are available. However, both experiments share donor **ENCDO236IKC**. Th17 versus regulatory T cell is a different cell context, not donor-independent confirmation. |

If a future review makes donor independence mandatory for ATAC, ENCSR639PPN
must be replaced rather than relabelled. The current pre-registered scope is
therefore explicitly library/cell-context confirmation, with shared-donor status
carried into the final interpretation.

ENCSR500KIT was rejected as the TF control because ordinary K562 does not match
ENCSR785RQR's CRISPR-CEBPB deletion background. ENCSR194GTV was selected instead
because its treatment background and donor match the TF treatment.

## Overlap checks

The comparison covered all original `file_type=fastq, output_type=reads` files
for discovery, candidates and selected controls. For histone it also included
the excluded rep2. The computed evidence is in
[`metadata/confirmation/independence-summary.json`](metadata/confirmation/independence-summary.json).

For all three groups, the following discovery/confirmation intersections are
empty:

- original FASTQ accessions and MD5 values;
- library accessions and biosample accessions;
- `(biosample, biological_replicate_number)` pairs.

Donor intersections are empty for TF and histone. ATAC has exactly the shared
donor ENCDO236IKC described above. Within each confirmation experiment, matched
treatment/control sharing a donor is intentional experimental pairing, not an
independent biological replicate.

Technical replicate libraries and multiple files from one library are not
counted as independent biological replicates. In particular, histone rep1 `1_1`
and `1_2` use different technical-replicate libraries but the same biosample.

## File-level confirmation manifest

[`confirmation-file-manifest.tsv`](confirmation-file-manifest.tsv) contains the
selected 28 FASTQ rows with official URL/size/MD5, read length, mate,
replicate/library/biosample/donor and control relationships. All 28 accessions
and all 28 MD5 values are unique, PE mates cross-reference each other within the
same replicate, and every row was checked against a separately fetched official
file JSON.

| Group | FASTQ files | Official bytes | Content |
| --- | ---: | ---: | --- |
| TF treatment ENCSR785RQR | 2 | 32,807,208,347 | PE 151 bp, one biological replicate |
| TF control ENCSR194GTV | 2 | 41,303,197,828 | Background-matched PE control |
| Histone treatment ENCSR213YKJ rep1 | 8 | 4,111,540,921 | SE 76 bp, two technical-replicate libraries |
| Histone control ENCSR148TRG rep1 | 8 | 3,350,521,776 | Same-biosample matched control |
| ATAC ENCSR639PPN | 8 | 11,262,765,629 | PE, four biological replicates |
| **Total planned, not downloaded** | **28** | **92,835,234,501** | **86.460 GiB** |

The TF confirmation payload alone is much larger than the TF discovery payload;
it is a separate budget/execution decision. A low-depth subset cannot substitute
for full-depth evidence when freezing a duplicate default.

## GEO/SRA optional confirmation candidates

[`geo-sra-confirmation-file-manifest.tsv`](geo-sra-confirmation-file-manifest.tsv)
adds 10 FASTQ / 3,371,454,387 bytes for CUT&Tag and MNase. These remain
**optional** until the limitations in
[`cuttag-mnase-selection.md`](cuttag-mnase-selection.md) are accepted:

- GSE145187 CUT&Tag uses different K562 source material and no RUN shared with
  the ALL-1 discovery set, but both treatments share one IgG, the replicate
  cell counts differ, and the IgG description is internally inconsistent.
- PRJNA928451 HP_pMSCs_FM uses different BioSamples/RUNs from LP discovery, but
  it is the same study and donor metadata do not establish independence. It is
  same-study condition confirmation only.

No RUN accession is shared between the GEO/SRA discovery and optional
confirmation rows. That check does not prove donor independence. Across all
currently registered confirmation candidates (ENCODE plus GEO/SRA), the planned
payload is 38 FASTQ / 96,206,688,888 bytes, not a completed download.

## Official evidence

The saved official responses and query logs under
[`metadata/confirmation/`](metadata/confirmation/) include:

- experiment JSON for ENCSR785RQR, ENCSR194GTV, ENCSR500KIT, ENCSR213YKJ,
  ENCSR148TRG, ENCSR639PPN and every discovery experiment;
- 53 individual FASTQ file JSON objects plus the related
  library/biosample/donor objects;
- [`object-query-plan.tsv`](metadata/confirmation/object-query-plan.tsv) and
  [`query-log.tsv`](metadata/confirmation/query-log.tsv), recording 113
  successful metadata/small-file downloads totalling 2,743,943 bytes.

Experiment queries use
`https://www.encodeproject.org/experiments/<ACCESSION>/?format=json`; file
queries use `https://www.encodeproject.org/files/<ACCESSION>/?format=json`.
The embedded experiment FASTQ entries and separately fetched file objects agree
on accession, size, MD5, read length, mate and status. No sequencing payload was
downloaded during this review.
