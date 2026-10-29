# MNase Rep1 discovery: not started

Status: **no conclusion.** This slice has not been downloaded or run. No
default is frozen. The pre-registered decision rules in
[`comparison-design.md`](comparison-design.md) are unchanged.

CUT&Tag Rep1 retain, Rep1 exclude, and Rep2 `max_insert` 500 are already
recorded on draft PR #197 (`57d8097eb96f431a543ac4d2fba94c513244fa7e`). That
result is inconclusive. This file does not rerun those arms and does not
treat their CI skips as a reason to repeat them.

## Why this is the next slice

[`resource-budget.md`](resource-budget.md) section 6 places MNase Rep1 after
the CUT&Tag expansion, before histone, ATAC, and the remaining discovery
groups. The sample is human, so the already built GRCh38 primary index can
be reused. ATAC ENCSR844TIU cannot: it is mouse and the GRCm38 FASTA
(773,507,376 bytes, official MD5 `3bc591be24b77f710b6ba5d41022fc5a`) is not
local. Histone and TF are later and larger. Confirmation is not this slice.

## Sample

PRJNA928451 LP_pMSCs_FM Rep1. BioSample SAMN32933787, experiment
SRX19193323, RUN SRR23247654. PE, no IgG control. SRA design text says the
three LP replicates are independent biological replicates. This hold is
Rep1 only. Rep2 (SRR23247653) and Rep3 (SRR23247651) stay out.

Read length is inferred: ENA `base_count / read_count = 102`, consistent
with 2x51 bp if the mates are equal. That length is checked from the FASTQ
after an authorized download, not before. Structured metadata says NovaSeq
6000. The source analysis used hg38; this slice would realign to the fixed
GRCh38 primary index, not download that FASTA again.

| Mate | URL | Official bytes | Official MD5 |
| --- | --- | ---: | --- |
| `SRR23247654_1.fastq.gz` | `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR232/054/SRR23247654/SRR23247654_1.fastq.gz` | 724,009,428 | `cc5ab28d826b404d0107a7e3e23364e0` |
| `SRR23247654_2.fastq.gz` | `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR232/054/SRR23247654/SRR23247654_2.fastq.gz` | 812,204,527 | `7c4c8288a2b6ce4e3422c293717335f2` |
| **Total** |  | **1,536,213,955** |  |

Local FASTQ directories used for CUT&Tag do not contain these files.
Free space when this hold was written was about 460 GiB on the shared
filesystem. That reading is not an admission check. Remeasure before any
later download and keep the 20 GiB margin in `resource-budget.md`.

## Arms, after authorization only

Candidate: `duplicate_analysis=exclude`, `max_insert=500`. One variable per
arm. New outdir and identity each time. `auto` is a dry-run and must
resolve to those two values. It is not an alignment.

No executed BAM handoff exists. The retain arm is a full workflow, not a
copy of the exclude BAM.

| Arm | Only change | Alignment |
| --- | --- | --- |
| exclude, max_insert 500 | baseline | new |
| retain, max_insert 500 | `duplicate_analysis` | new full workflow |
| exclude, max_insert 200 | `max_insert` | new Bowtie2 alignment |

200 versus 500 is a sensitivity check on the supported mono/di windows.
500 stays the candidate. 200 is not a proposed default. One Rep1 library
cannot supply between-replicate Jaccard. Missing mono/di coverage limits
scope. Confirmation samples are not used, and this discovery RUN is not
the optional HP confirmation RUNs SRR23247650 and SRR23247648.

Isolated MACS3 `--keep-dup all` versus `1` waits until the exclude/500
analysis BAM exists. It is not a production switch and is not a reason to
download anything else.

## What this authorization would not include

| Blocked item | Official bytes |
| --- | ---: |
| MNase Rep2 + Rep3 | 2,865,847,578 |
| Histone discovery ENCSR817LUF + ENCSR603VTC | 2,510,695,554 |
| ATAC discovery ENCSR844TIU | 3,952,560,652 |
| Mouse GRCm38 FASTA, also required for ATAC | 773,507,376 |
| TF discovery ENCSR054ZMK + ENCSR142JFU | 11,396,944,028 |
| Confirmation, all registered candidates | 96,206,688,888 |

The completed CUT&Tag download, about 60 MB, does not authorize any row
above or the 1,536,213,955-byte Rep1 pair.

## Verdict

No scientific result. The next discovery comparison stays inconclusive
until these two FASTQ files are explicitly authorized, downloaded, checked
against the official MD5 values, and run under the arms above.
