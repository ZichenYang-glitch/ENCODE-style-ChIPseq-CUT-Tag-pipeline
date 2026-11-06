# SCI-7 proposed reference identity

Status: sources and small metadata are fixed as of 2026-10-28. Large FASTA/GTF
and Bowtie2 indexes have **not** been downloaded or built. Official compressed
MD5 values below are source metadata, not a claim of local payload verification.
Raw FASTQ records do not have an official assembly attribute inferred from this
plan; these are the proposed references for fresh realignment.

The machine-readable summary is
[`metadata/confirmation/reference-identity.tsv`](metadata/confirmation/reference-identity.tsv).

## 1. Human: GENCODE v44 / GRCh38.p14

- Input: `GRCh38.primary_assembly.genome.fa.gz`, GENCODE Release 44,
  GRCh38.p14 primary assembly (reference chromosomes and scaffolds, without
  patches/haplotypes).
- Official URL:
  `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz`
- Official compressed MD5 from the release `MD5SUMS`:
  `9c3fc2ca260a767530dddb0f26721a6b`; HTTPS HEAD reports 844,691,642 bytes.
- Contigs: UCSC-style `chr1`…`chr22`, `chrX`, `chrY`, `chrM` and scaffolds.
  Bind `organelle_chroms: ["chrM"]` explicitly.
- Payload status: **unverified** until downloaded and hashed.

## 2. Mouse: GENCODE vM25 / GRCm38.p6 (mm10)

- Input: `GRCm38.primary_assembly.genome.fa.gz`, GENCODE Release M25,
  GRCm38.p6 primary assembly.
- Official URL:
  `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/GRCm38.primary_assembly.genome.fa.gz`
- Official compressed MD5 from the release `MD5SUMS`:
  `3bc591be24b77f710b6ba5d41022fc5a`; HTTPS HEAD reports 773,507,376 bytes.
- Contigs: `chr1`…`chr19`, `chrX`, `chrY`, `chrM` and scaffolds.
  Bind `organelle_chroms: ["chrM"]` explicitly.
- Payload status: **unverified** until downloaded and hashed.

The mm10/GRCm38 identity must not be silently replaced with mm39/GRCm39 while
reusing this blacklist or index identity.

## 3. Derived FAI, chrom_sizes and Bowtie2 index

The locked execution closure contains Bowtie2 2.5.5 and samtools 1.23.1
(`workflow/envs/align.lock`, package MD5 `af24eb11ce11479ee06d7dc499da0a28`
and `89624fbd17c069abcbc8b19b96d497a0` respectively). Derive every resource
from the same decompressed FASTA for each species:

```bash
samtools faidx ASSEMBLY.primary_assembly.genome.fa
cut -f1,2 ASSEMBLY.primary_assembly.genome.fa.fai \
  > ASSEMBLY.primary_assembly.genome.chrom.sizes
bowtie2-build --threads 8 \
  ASSEMBLY.primary_assembly.genome.fa \
  ASSEMBLY.primary_assembly.genome
```

`--threads 8` is an execution setting, not a scientific parameter; all other
`bowtie2-build` parameters are defaults. The sample sheet's `bowtie2_index` is
the corresponding prefix, and all six `.bt2` files must come from this one
build. Do not mix a prebuilt index or another release's FASTA/FAI/chrom_sizes.

At execution, record compressed and decompressed FASTA SHA256, FAI/chrom_sizes
SHA256, all six index-file SHA256 values, exact tool versions and commands.
These derived artifacts are currently **unbuilt and unverified**.

`.dict` is not required by the default validation plan because
`qc.picard_metrics: false`. If that QC is enabled, derive a dictionary from the
same FASTA using samtools/Picard rather than borrowing another assembly.

## 4. ENCODE blacklist identity and verified small payloads

Both files belong to ENCODE annotation ENCSR636HFF, ENCODE DAC Exclusion List
Regions (formerly DAC Blacklisted Regions).

| Assembly | Fixed file | Official compressed MD5 | Local small-payload verification |
| --- | --- | --- | --- |
| GRCh38 | ENCFF356LFX, GRCh38 unified/v2 distribution, 910 BED intervals | `393688b4f06c9ce26165d47433dd8c37` | Downloaded 8,211 bytes; byte size and MD5 match the official file JSON |
| mm10 | ENCFF547MET, `mm10.blacklist.bed.gz`, 164 BED intervals | `5d19eeb4197ae10540252c61399f3131` | Downloaded 1,263 bytes; byte size and MD5 match the official file JSON |

Official URLs:

- `https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz`
- `https://www.encodeproject.org/files/ENCFF547MET/@@download/ENCFF547MET.bed.gz`

The verified compressed files and decompressed inspection copies are stored as
[`metadata/confirmation/ENCFF356LFX.bed.gz`](metadata/confirmation/ENCFF356LFX.bed.gz)
and
[`metadata/confirmation/ENCFF547MET.bed.gz`](metadata/confirmation/ENCFF547MET.bed.gz).
The GRCh38 BED contains chr1–22/X/Y and the mm10 BED chr1–19; a blacklist need
not cover every reference contig. Once the actual FASTA exists, SCI-4 admission
must still verify every BED contig and `0 <= start < end <= contig length`.
MD5 equality does not prove coordinate compatibility with an arbitrary FASTA.

## 5. Optional annotation

GTF is not required by the default comparison because `qc.tss_enrichment:
false`. If it is enabled, use the primary-assembly annotation from the same
GENCODE release:

- human `gencode.v44.primary_assembly.annotation.gtf.gz`, official MD5
  `b182a9f3b134b9cc2da566a2b1692557`, HTTPS HEAD size 49,730,393 bytes;
- mouse `gencode.vM25.primary_assembly.annotation.gtf.gz`, official MD5
  `c5125258a0a2c5250ddb4c192abbf4e8`, HTTPS HEAD size 28,556,861 bytes.

Neither GTF was downloaded. The HEAD results are saved in
[`metadata/confirmation/reference-head-checks.json`](metadata/confirmation/reference-head-checks.json). When enabled, verify GTF contig compatibility with
the same FASTA; do not substitute CHR-only or ALL annotation while retaining
the primary-assembly identity.

## 6. Admission checks before execution

1. FASTA, FAI, chrom_sizes and BAM `@SQ` agree on same-name contig lengths;
   FAI offsets/line widths match the actual FASTA bytes.
2. Blacklist BED is non-empty, TAB-delimited, has integer legal coordinates,
   known contigs and no out-of-bounds intervals.
3. `organelle_chroms=["chrM"]` resolves in each reference; ATAC auto-exclude
   and ChIP auto-report therefore use a real contig rather than a guessed
   accession.
4. Reference binding v2 records actual content hashes and assembly identity.
5. Raw FASTQ is realigned to this fixed index. An ENCODE BAM's assembly label
   is only supporting evidence for assembly-family selection and never skips
   local reference verification.

## 7. Official source evidence

Saved source files under [`metadata/confirmation/`](metadata/confirmation/)
include the GENCODE release pages, `MD5SUMS` and README files, ENCODE blacklist
file JSON and annotation ENCSR636HFF JSON. Official upstream locations:

- `https://www.gencodegenes.org/human/release_44.html`
- `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/MD5SUMS`
- `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/_README.TXT`
- `https://www.gencodegenes.org/mouse/release_M25.html`
- `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/MD5SUMS`
- `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/_README.TXT`
- `https://www.encodeproject.org/annotations/ENCSR636HFF/?format=json`

No large reference payload was fetched during this preparation step.
