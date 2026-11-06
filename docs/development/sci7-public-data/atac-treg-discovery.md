# ATAC Treg discovery: measured, no conclusion

Status: **no conclusion.** ENCSR844TIU was downloaded, the four FASTQ files
matched the official sizes and MD5 values, and the three pre-registered arms
each finished with exit status 0. No default is frozen. The decision rules in
[`comparison-design.md`](comparison-design.md) are unchanged. The two
biological replicates share donor ENCDO236IKC, confirmation has not been
evaluated, and 26 bp reads remain an explicit limitation. Exit status 0 is
not scientific acceptance.

The candidate settings are `duplicate_analysis=exclude` and `max_insert=2000`.
`auto` was a dry-run only. It printed
`duplicate_analysis=auto(atac:exclude)` and `max_insert=auto(atac:2000)`
and was not a third alignment.

## Sample

ENCODE experiment ENCSR844TIU, mouse C57BL/6NJ regulatory T cells, donor
ENCDO236IKC. Two biological replicates, libraries ENCLB190GZI and
ENCLB265DZB. Layout is PE. Read length is 26 bp. No control is used.
The four official FASTQ files total 3,952,560,652 bytes.

| File | Replicate | Official bytes | Official MD5 |
| --- | --- | ---: | --- |
| ENCFF789PHQ (R1) | 1 | 740,537,844 | `ffc5e74e3386b021f5290c6422f1ebf1` |
| ENCFF574DTB (R2) | 1 | 748,607,968 | `8ef971a1e83016970c54b61f9c510469` |
| ENCFF488RVS (R1) | 2 | 1,227,795,683 | `ce613825e138a86f52eed450bc0861b7` |
| ENCFF575LXK (R2) | 2 | 1,235,619,157 | `b65f99c481241ae778b7159e587b8a19` |

The reference is GENCODE Release M25 `GRCm38.primary_assembly.genome.fa.gz`,
773,507,376 bytes, official MD5 `3bc591be24b77f710b6ba5d41022fc5a`, SHA256
`617b10dc7ef90354c3b6af986e45d6d9621242b64ed3a94c9abeac3e45f18c17`. FAI,
`chrom_sizes`, and the six Bowtie2 `.bt2` files were derived from that
decompressed FASTA with samtools 1.23.1 and Bowtie2 2.5.5
(`bowtie2-build --threads 8`, other build parameters left at defaults).
The FASTA has 66 contigs. The local UCSC/genomepy GRCm38 index was not used.
Blacklist ENCFF547MET, 164 intervals, official compressed MD5
`5d19eeb4197ae10540252c61399f3131`, was admitted against those contig
lengths. `organelle_chroms` is `chrM`.

The executed checkout is
`634160e0d7c671db466a22a24bbe2287782b9af7`.

All three real arms used the updated profile, schema 2.0.0, policy revision
6, `organelle_mode` resolved to `exclude`, blacklist mode `analysis`, and
MACS3 `-f BAMPE -g mm -q 0.01 --keep-dup all`. Narrow peak mode does not pass
`--broad`. Picard `MarkDuplicates` used `REMOVE_DUPLICATES=false`. An exclude
arm then passed `--exclude-dup` when building the analysis BAM. A retain arm
does not. `max_insert` is Bowtie2 `-X`.

| Tool | Version |
| --- | --- |
| Snakemake | 8.30.0 |
| Bowtie2 | 2.5.5 |
| samtools | 1.23.1 |
| MACS3 | 3.0.4 |
| Picard | 3.5.0 |
| bedtools | 2.31.1 |
| bedGraphToBigWig | 2.10 |

## Arms

Each arm is a full workflow in its own output directory. There is no BAM
handoff. Only one configured variable changes from the candidate.

| Arm | Only change | Identity | Exit | Wall |
| --- | --- | --- | --- | ---: |
| exclude, max_insert 2000 | candidate | `0dd4b53ff0c214d6187b59c152803dba14538bb241bd9eaca4c559b77b152fb8` | 0 | 2010 s on the resumed run |
| retain, max_insert 2000 | `duplicate_analysis` | `b5801aeba447888c4323b7a50472f3557889a3821a71b3517b6ade74c9a4c9b7` | 0, 78/78 | 8130 s |
| exclude, max_insert 500 | `max_insert` | `ab7f621d5c7e85ec5fdb5bcbc884a2b4d0ee5e3a1d4ae77872e64fe7b789c500` | 0, 78/78 | 5820 s |

The candidate's first attempt stopped while converting a fold-enrichment
bedGraph to bigWig because `bedGraphToBigWig` was not on the run PATH. The
binary is bedGraphToBigWig 2.10. The resumed run kept the same config and
finished the remaining 41 jobs with exit 0. Retain and the 500 arm each
scheduled the full 78 jobs.

Pre-duplicate Picard metrics on the MAPQ BAM match between exclude/2000 and
retain/2000. Replicate 1 `percent_duplication` is 0.099213
(2,891,122 / 29,140,461 read pairs). `estimated_library_size` is 137,446,248.
NRF, PBC1, and PBC2 are 0.900787, 0.934258, and 17.302708. Replicate 2
`percent_duplication` is 0.142224 (6,840,489 / 48,096,556).
`estimated_library_size` is 152,931,350. NRF, PBC1, and PBC2 are 0.857776,
0.902455, and 11.272717. The 500 arm is a different alignment, so its
pre-duplicate duplication is 0.105518 and 0.150185 and is not required to
match.

Analysis-BAM read counts, after duplicate policy, organelle exclusion, and
blacklist filtering:

| Arm | Replicate 1 kept | Replicate 2 kept |
| --- | ---: | ---: |
| exclude, max_insert 2000 | 51,253,190 | 80,814,290 |
| retain, max_insert 2000 | 54,633,286 | 89,435,170 |
| exclude, max_insert 500 | 46,406,658 | 73,451,298 |

On exclude/2000, duplicate flags account for 5,782,244 dropped alignments in
replicate 1 and 13,680,978 in replicate 2. The retain arm drops none for that
reason.

## Duplicate effect

Effects are per biological replicate. The common universe for a replicate is
the merged blacklist-filtered narrow peaks of retain/2000 and exclude/2000
for that replicate. Own-peak FRiP is descriptive.

| | Replicate 1 | Replicate 2 |
| --- | ---: | ---: |
| Retain peaks | 62,660 | 73,296 |
| Exclude peaks | 63,052 | 71,826 |
| Merged-base Jaccard, retain versus exclude | 0.947501 | 0.943664 |
| Intersection bp | 44,221,051 | 52,190,695 |
| Union bp | 46,671,242 | 55,306,463 |
| Retain peaks touched by an exclude peak | 0.976716 | 0.963504 |
| Own-peak FRiP, retain | 0.506257 | 0.511346 |
| Own-peak FRiP, exclude | 0.494856 | 0.494654 |
| Common-universe FRiP, retain | 0.509758 | 0.514357 |
| Common-universe FRiP, exclude | 0.498128 | 0.498830 |
| Relative FRiP change | −0.022815 | −0.030187 |
| Relative off-peak change | 0.023723 | 0.031972 |

The joint trigger requires a relative FRiP change strictly below −0.10 and a
relative off-peak change strictly above 0.05, in every replicate. Neither
inequality holds in either replicate. The joint trigger is false.

Within-arm between-replicate merged-base Jaccard is 0.726794 on retain and
0.734423 on exclude. The exclude value is 0.007629 higher. That is not a
fall of more than 0.10. Of the retain peaks that overlap the other
replicate, exclude loses 0.014880 in replicate 1 (881 / 59,205) and 0.008857
in replicate 2 (521 / 58,824). The loss is not above 0.10 in every replicate.

## max_insert effect

The registered short-template denominator is MAPQ ≥ 30 and length ≤ 500 bp
in the 500 arm. The event is an R1 coordinate, R2 coordinate, or ordered
mate-contig change in the 2000 arm. A 500-arm template absent from the 2000
arm stays in the denominator and is not a coordinate change. The trigger is
a fraction strictly greater than 0.01 in every replicate.

| | Replicate 1 | Replicate 2 |
| --- | ---: | ---: |
| Denominator | 26,570,383 | 44,116,609 |
| Changed | 30,475 | 51,984 |
| Fraction | 0.001147 | 0.001178 |
| Templates only in the 2000 arm | 2,732,126 | 4,238,598 |
| Of those, MAPQ ≥ 30 | 2,592,616 | 4,019,533 |
| Of those, length 501–700 bp | 2,225,950 | 3,490,489 |
| Of those, length > 700 bp | 504,362 | 744,599 |

The short-template trigger is false in both replicates. Most templates gained
at 2000 are 501–700 bp with MAPQ ≥ 30. They are not predominantly low-MAPQ,
and they do not coincide with the signal or peak-loss triggers. The 500 arm
has no templates longer than 500 bp, which is what `-X 500` enforces. Peak
counts on exclude/500 are 68,422 and 79,710. Those counts are descriptive.

Within-arm between-replicate Jaccard on exclude/500 is 0.715587. The
exclude/2000 Jaccard is 0.018835 higher, so the 2000 candidate does not fall
by more than 0.10.

Length-binned coverage correlation uses 100 kb windows on the occupied union
only. Empty windows are not included, and the values are not used as quality
evidence. Pearson correlation between the two replicates is about 0.999 in
every occupied length bin of the 2000 alignment. The 500 alignment has the
same pattern below 500 bp and no occupied windows above 500 bp. The
confirmation comparison of those correlations was not run.

## Verdict

No default is frozen. Exclude does not meet the joint signal rule, the
short-template fraction is about 0.11%, and replicate concordance does not
fall by more than 0.10. The two biological replicates still come from one
donor, the reads are 26 bp, and confirmation has not been evaluated. That
is not enough to freeze `exclude` or `max_insert` 2000.
