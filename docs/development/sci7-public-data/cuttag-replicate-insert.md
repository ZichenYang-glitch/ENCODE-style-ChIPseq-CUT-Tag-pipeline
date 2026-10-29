# CUT&Tag Rep1 duplicate policy, submitter-replicate peaks, and max_insert

Status: **measured, not accepted.** Rep1 retain versus exclude and a Rep2
`max_insert` 500 versus 700 contrast are recorded. No default is frozen.
Rep1 and Rep2 are submitter labels. GEO does not call them independent
biological replicates. Confirmation data were not used.

The comparison rules were fixed in
[`comparison-design.md`](comparison-design.md) before these arms ran.
Peak count is descriptive. Common-universe FRiP uses epsilon `1e-6`.

## Pairing

The pair is the registered proposal in
[`file-manifest.tsv`](file-manifest.tsv) and
[`cuttag-mnase-selection.md`](cuttag-mnase-selection.md): H3K27me3
SRR12842271 (GSM4835654, SRX9309566, SAMN16477147) with IgG SRR12842281
(GSM4835664, SRX9309576, SAMN16477220). The match is patient, Rep1 label,
and experiment date, not an SRA `paired_with` field.

Official ENA `read_count` is 941,245 for the treatment and 534,078 for the
IgG, against 5,515,652 and 3,621,858 for the already measured Rep2 pair.
Rep1 is shallow. That depth limit is carried into the verdict below.

## Provenance

Both Rep1 arms used worktree
`/home/yangzichen/workflow/ENCODE-sci7-cuttag-ext` at
`28bf02ac11a86a52afe4fce029c3961d32dc9f8d`. `workflow/` and `src/` are
unchanged from `f93bf0c110f8c0d859e959be27448e780a79e4c6`.
`verify-checkout` passed. The imported package was
`ENCODE-sci7-cuttag-ext/src/encode_pipeline/__init__.py`.

| Tool | Identity |
| --- | --- |
| Snakemake | 8.30.0 |
| Bowtie2 | 2.5.5 (`bowtie2-align-s`; x86-64-v3 fallback warning, then the default binary) |
| samtools | 1.23.1 |
| MACS3 | 3.0.4 |
| Picard | 3.5.0 |
| bedtools | 2.31.1 |
| FastQC | 0.12.1 |
| Trim Galore | 2.2.0 |
| bamCoverage | 3.5.6 |

Profile `updated`, policy revision 6, schema 2.0.0, MAPQ 30,
`max_insert` 700, organelle `auto` resolved `exclude` with `chrM`,
MACS3 `-f BAMPE -g hs -q 0.01 --keep-dup all -B`, 8 cores.
The two configs differ only in `outdir` and `duplicate_analysis`.

`auto` was a dry-run on the Rep1 sample sheet. It resolved to
`duplicate_analysis=auto(cuttag:retain)` and `max_insert=auto(cuttag:700)`.
The dry-run did not write a result directory and is not an arm.

## Inputs and reference

Rep1 FASTQ bytes and official MD5 values matched. SHA256:

| File | Bytes | SHA256 |
| --- | ---: | --- |
| SRR12842271_1.fastq.gz | 18,696,899 | `3a0e6d38e8ff15a16b8f642d5684498d6a4c1358408dc491759d5279d378deca` |
| SRR12842271_2.fastq.gz | 19,474,117 | `0891f56c0283391f8cb053c1d0fa79862f69f23f9bb2f4ffd63aa48b30c40f22` |
| SRR12842281_1.fastq.gz | 10,738,383 | `ebee7b8e9a9677a877131e31f3bde2aac6b8a5942a27dc2803311302abf54390` |
| SRR12842281_2.fastq.gz | 11,036,019 | `77448e23b1208d492770de8356cc1dc796367914711b185ed599e4ed625d8f11` |

GRCh38 primary FASTA gzip MD5 `9c3fc2ca260a767530dddb0f26721a6b`, SHA256
`e9a2d5a5cd225293646ae298998ad4ea8c11e8f22729bc2d6f5c3dcfefbf8ef8`.
Decompressed FASTA SHA256
`e49b92b3e4f321bf254c042f25b726d9931c4d74c7523e8b6bb530e63b0cfd4b`.
The decompressed blacklist matches official ENCFF356LFX; BED SHA256
`5d214c58c07f33a4fa1d6bcfe61ff00c666510e8fee5520fc1015c5d3ef16f3b`.

## Identities

| Arm | `duplicate_analysis` | Identity | Exit |
| --- | --- | --- | --- |
| retain | retain | `1d570cb768b91cbc4dad7a5a7c49da6e0e5d996e61ae4aa661f6155c54cb2231` | 0, 48/48 |
| exclude | exclude | `be2673cbb6d8ee2a939bce8acd4e64ded5d626e7f58916c203105a751d752006` | 0, 48/48 |

Retain's analysis-BAM rule evaluated `False == True` and did not add
`--exclude-dup`. Exclude evaluated `True == True` and did. Both Bowtie2
commands used `-X 700`.

## Stage counts

NRF, PBC, and Picard percent duplication are pre-duplicate. They match
because that stage is upstream of the only policy difference.

| Metric | Stage | retain | exclude |
| --- | --- | ---: | ---: |
| Picard percent duplication | pre-duplicate MarkDuplicates | 0.023994 | 0.023994 |
| NRF | pre-duplicate MAPQ BAM | 0.976006 | 0.976006 |
| PBC1 | same stage | 0.975919 | 0.975919 |
| PBC2 | same stage | 41.367912 | 41.367912 |
| Analysis valid templates | analysis BAM | 781,696 | 762,944 |
| Analysis read records | `samtools view -c` denominator used by FRiP | 1,563,392 | 1,525,888 |
| Invalid analysis templates | analysis BAM | 0 | 0 |

The analysis-BAM difference is 18,752 templates and 37,504 reads, two
reads per template. MACS3 on the exclude arm read 762,944 treatment
fragments and 39,691 IgG fragments, then wrote an empty narrowPeak.
The raw narrowPeak is 0 bytes, so the empty blacklist-filtered file is
not a blacklist-only loss.

## Peaks and common-universe signal

| | retain | exclude |
| --- | ---: | ---: |
| Narrow peaks | 1,646 | 0 |
| Own-peak FRiP | 0.020901 | 0.000000 |

Merged-base Jaccard is 0. The common universe is the retain peak set,
863,164 bp. Reference arm is retain.

| | retain | exclude |
| --- | ---: | ---: |
| Reads in the common universe | 32,676 | 31,374 |
| FRiP | 0.020901 | 0.020561 |
| Off-peak fraction | 0.979099 | 0.979439 |

Exclude minus retain is -0.034 percentage points. Relative FRiP change
is -0.016267. `log2((0.020561+1e-6)/(0.020901+1e-6))` is -0.023660.
Off-peak background ratio, exclude over retain, is 1.000347.
Per-base enrichment on the 3,099,750,718 bp `chrom.sizes` total is
76.638 (retain) and 75.367 (exclude).

The pre-registered signal/background trigger is a FRiP drop of more than
10% relative and an off-peak increase of more than 5% relative, on every
independent replicate. This pair does not meet that numeric trigger.
It also does not satisfy the independent-replicate clause.

## Submitter-replicate peak agreement

Within retain, merged-base Jaccard of Rep1 versus the existing Rep2
retain peaks is 0.034583 (intersection 329,080 bp, union 9,515,560 bp).
`bedtools jaccard` returns the same intersection, union, and Jaccard.
1,646 Rep1 peaks touch 0.467801 of the Rep1 set against Rep2; 12,458
Rep2 peaks touch 0.057794 of the Rep2 set against Rep1.

Within exclude, Rep1 has no peaks, so Jaccard against the 2,881 Rep2
exclude peaks is 0.

These are submitter-labelled libraries on different dates, and Rep1's
781,696 valid templates are about 19% of Rep2 retain's 4,144,663. The
between-replicate regression trigger is not applicable. Length-binned
coverage correlation was not computed: the design asks for it when there
are at least two independent biological replicates, and these submitter
labels are not that.

## Rep1 resources

`/usr/bin/time -v` is the Snakemake parent. `du` is the final results
directory.

| Arm | Wall | User CPU (s) | System CPU (s) | Average CPU | Parent max RSS | Final results |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| retain | 11:23.89 | 1277.34 | 66.04 | 196% | 3,463,168 KiB | 516M |
| exclude | 11:10.12 | 1214.30 | 62.65 | 190% | 3,462,656 KiB | 497M |

The two Rep1 arms did not run at the same time. Per-arm scratch high-water
from directory samples was 540,447,337 bytes (retain), 520,444,930 bytes
(exclude), and 3,531,647,853 bytes (`max_insert` 500). Those figures are
not added together. Free space on the shared filesystem moved from about
466G to 460G across the whole sequence and is not attributed to one arm.

## max_insert 500 versus reused 700

The 700 side is the existing Rep2 retain arm, identity
`baf233c5c47b1ab001ad2fcb0c2c06709f830f879025c8d2ee9c4fd3518410b6`.
`workflow/` and `src/` match this checkout. Inputs, MAPQ 30, profile,
blacklist, and `duplicate_analysis: retain` match. The 500 arm is a new
Bowtie2 alignment, identity
`be44e0e663dfd1af055279a17d68a338dc90b133b487dc21edc9ba41551e49c1`,
exit 0, 48/48. Both T2 and IgG2 commands used `-X 500`. The analysis-BAM
rule evaluated `False == True` and did not add `--exclude-dup`.

Parent `/usr/bin/time -v`: wall 13:13.91, user CPU 3152.94 s, system CPU
123.18 s, 412% CPU, max RSS 3,481,088 KiB. Final results directory 3.3G.

Analysis-BAM stage, after the same downstream filters:

| | 500 | 700 |
| --- | ---: | ---: |
| T2 valid templates | 3,281,745 | 4,144,663 |
| T2 invalid templates | 0 | 0 |
| T2 narrow peaks | 10,177 | 12,458 |
| T2 own-peak FRiP | 0.033093 | 0.043496 |

Own-peak FRiP values do not share a universe. Picard percent duplication
is 0.009353 at 500 and 0.008728 at 700 because duplicate marking sees
different alignments. That is not a duplicate-policy effect.

Common-universe signal uses the merged peaks, 10,513,096 bp.
`bedtools jaccard` intersection 4,019,353, union 10,513,096, Jaccard
0.382319, matching the script. Reference for the log2 ratio is the 700
candidate. Reads and FRiP denominators are analysis-BAM read records.

| | 500 | 700 |
| --- | ---: | ---: |
| Reads in the common universe | 335,770 | 416,915 |
| Read-record denominator | 6,563,490 | 8,289,326 |
| FRiP | 0.051157 | 0.050295 |
| Off-peak fraction | 0.948843 | 0.949705 |

500 relative to 700 is +0.0862 percentage points. Relative change of the
500 FRiP versus 700 is +0.017139.
`log2((0.051157+1e-6)/(0.050295+1e-6))` is +0.024516.
Off-peak ratio, 500 over 700, is 0.999092.
The candidate arm is 700. Its FRiP relative to 500 changes by -0.016850
and its off-peak fraction changes by +0.000908. That does not meet the
trigger, which requires a FRiP drop of more than 10% relative and an
off-peak increase of more than 5% relative. Enrichment on the
3,099,750,718 bp reference is 15.843 at 500 and 15.562 at 700.

Primary proper pairs from the aligner BAMs, before MAPQ and duplicate
filters:

| Library | Templates at 500 | Templates at 700 | Only in 500 | Newly in 700 | Incomplete pairs |
| --- | ---: | ---: | ---: | ---: | ---: |
| T2 | 3,518,396 | 4,421,146 | 50 | 902,800 | 0 / 0 |
| IgG2 | 2,193,162 | 2,519,642 | 47 | 326,527 | 0 / 0 |

T2 templates newly present at 700 are 902,324 in the 501–700 bp bin
(859,421 with minimum MAPQ ≥ 30 and 42,903 below 30) and 476 at or below
500 bp. None are longer than 700. IgG2 newly present templates are
326,185 in the 501–700 bp bin (296,087 with minimum MAPQ ≥ 30). The 50
and 47 templates present only at 500 show that raising `-X` is not a
filter of the 700 BAM.

High-confidence short templates are minimum MAPQ ≥ 30 and length ≤ 500
in the 500 arm. A missing 700 record stays in the denominator and is not
counted as a coordinate change.

| Library | Denominator | Coordinate or mate changed | Fraction |
| --- | ---: | ---: | ---: |
| T2 | 3,284,411 | 4,019 | 0.001224 |
| IgG2 | 1,878,855 | 3,802 | 0.002024 |

Both fractions are below the 1% trigger.

`auto` was not run as a third alignment. The Rep1 dry-run already
resolved CUT&Tag `auto` to retain and 700.

## Verdict

Inconclusive for a default.

- Rep1 retain and exclude finished 48/48 with exit 0, distinct identities,
  invalid templates 0, and a read difference of exactly two reads per
  dropped template. Common-universe FRiP does not meet the
  signal/background trigger. The exclude arm calls no peaks on this
  shallow library. That peak count is descriptive.
- Submitter-replicate peak Jaccard inside retain is 0.034583. These
  libraries are not established independent biological replicates, and
  Rep1 has about 19% of Rep2 retain's valid templates. The
  between-replicate regression trigger is not applicable. Length-binned
  coverage correlation was not used as quality evidence.
- The 500 arm was realigned. On the common peak universe the 700 candidate
  does not meet the signal/background trigger. Fewer than 1% of
  high-confidence short templates change coordinates or mate assignment.
  Most templates gained at 700 are 501–700 bp, and many of those have
  MAPQ ≥ 30, so the gain is not only low-MAPQ sequence.
- Confirmation has not been run. One shallow pair plus one deeper pair of
  submitter labels cannot freeze `retain` or `max_insert` 700.

TF ChIP, histone ChIP, ATAC-seq, and MNase-seq discovery were not run.
Their registered FASTQ downloads are not covered by this round's CUT&Tag
authorization. ATAC also still needs an unverified GRCm38 reference.
Confirmation FASTQ, 38 files and 96,206,688,888 official bytes, was not
downloaded.

