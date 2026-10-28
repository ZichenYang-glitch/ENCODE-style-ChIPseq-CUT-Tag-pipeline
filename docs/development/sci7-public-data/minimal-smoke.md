# Minimal SCI-7 smoke execution checklist

Status: proposed for a later explicit execution authorization. Nothing in this
checklist has been downloaded or run. This is an execution/storage smoke test,
not biological validation and not a default-freezing experiment.

## 1. Authorized input set

Use only CUT&Tag Rep2 treatment SRR12842272 and proposed matched IgG Rep2
SRR12842282. Do not add Rep1 or any other cohort in this step.

| Role/mate | HTTPS URL | Official bytes | Official MD5 |
| --- | --- | ---: | --- |
| Treatment R1 | `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/072/SRR12842272/SRR12842272_1.fastq.gz` | 113,805,899 | `48b1d2a23ecb08c3e7d3d1456f3ece0f` |
| Treatment R2 | `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/072/SRR12842272/SRR12842272_2.fastq.gz` | 118,513,198 | `c5497efbbb374886c25d75aed5366daa` |
| IgG R1 | `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/082/SRR12842282/SRR12842282_1.fastq.gz` | 74,581,655 | `c3c95099d4c523bf34a4b107a8246c49` |
| IgG R2 | `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR128/082/SRR12842282/SRR12842282_2.fastq.gz` | 76,732,766 | `3e87566c84e3b94e3280a2dae6b28e8f` |
| **Total FASTQ** |  | **383,633,518** |  |

Reference input, if not already present:

- `https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz`
- 844,691,642 bytes; official compressed MD5
  `9c3fc2ca260a767530dddb0f26721a6b`.

New download with FASTA: **1,228,325,160 bytes (1.144 GiB)**. GTF is not
required. Use the already verified local ENCFF356LFX blacklist (8,211 bytes,
MD5 `393688b4f06c9ce26165d47433dd8c37`) rather than downloading it again.

## 2. Preflight and download

1. Record exact checkout SHA, tool versions and `df --output=avail`. Apply the
   initial smoke guardrail in `resource-budget.md`: at least 20 GiB free before
   measurements exist. This is an operational guardrail, not a usage estimate.
2. Download only the five files above into task-owned input/reference
   directories. Verify each byte count and MD5 before any alignment. Stop on a
   mismatch; do not silently substitute another file.
3. Record download commands, local paths and SHA256 in the task evidence
   manifest. Keep the original compressed FASTQ and FASTA.

## 3. Reference build and admission

1. Decompress the FASTA in the task reference directory; record compressed and
   decompressed SHA256 and allocated space.
2. Run samtools 1.23.1 `faidx`, derive chrom_sizes from that FAI, and build a
   Bowtie2 2.5.5 index with default build parameters from the same FASTA.
3. Verify FASTA/FAI/chrom_sizes contigs and lengths, blacklist contigs and
   intervals, and `organelle_chroms=["chrM"]`. Record all six index-file
   SHA256 values and reference/index allocated space.
4. Do not enable GTF-dependent QC and do not download another reference or
   prebuilt index.

## 4. One default workflow run

1. Create samples for H3K27me3 Rep2 and IgG Rep2 as CUT&Tag PE, with the
   proposed treatment/control pairing. Use GRCh38 primary reference and
   blacklist analysis.
2. Keep `duplicate_analysis: auto` (must resolve to retain for both treatment
   and IgG) and Bowtie2 `max_insert: auto` (must resolve to 700). Do not add a
   legacy `--keep-dup` or `-X` override.
3. Validate the effective configuration and create a new outdir/run identity.
   For the planned updated-profile reuse measurement, first verify the locked
   Picard admission; if it is unavailable, stop and report the prerequisite
   instead of silently changing profile.
4. Run only this default arm through alignment, duplicate marking, analysis
   BAM, peaks, coverage and applicable QC. Do not run the 500 max_insert arm
   or the duplicate-exclude arm in the smoke.
5. Verify resolved policies, template/count conservation, expected peaks/QC
   and output checksums. Preserve final/analysis BAMs, metrics, logs, commands
   and effective run records as long-term smoke evidence.

## 5. Isolated MACS3 keep-dup check

Freeze the smoke analysis BAM and matched IgG BAM. In two separate evidence
directories, invoke the same pinned MACS3 directly with identical parameters
except `--keep-dup all` versus `--keep-dup 1`. Save commands, versions, logs and
peak outputs. Do not modify production extra_args, the validator or run
identity. This checks the mechanism on a real library; it does not choose a
biological default.

## 6. Measure and stop

At reference build, Bowtie2/MAPQ, duplicate handling, analysis-BAM derivation,
peaks/QC and the direct MACS3 calls, record `df`, allocated `du`, output sizes
and scratch high-water marks as specified in `resource-budget.md`.

Stop after this smoke. Report actual reference, BAM, temporary and long-term
sizes and any failure. Later duplicate-exclude handoff and max_insert 500
realignment require the next execution authorization and a storage projection
using these measurements. No cleanup of raw inputs, references or evidence is
part of this checklist; only explicitly listed task-owned scratch can be
considered under the applicable authorization.
