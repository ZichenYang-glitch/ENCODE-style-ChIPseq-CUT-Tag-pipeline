# CUT&Tag duplicate retain versus exclude

Status: **single-pair mechanism and effect-size check**. This is not biological
acceptance, not a replicate-concordance result, and it does not freeze a
default. Peak count is descriptive only.

The comparison rules were fixed in
[`comparison-design.md`](comparison-design.md) before these arms ran. The
common-universe FRiP formula, epsilon `1e-6`, and off-peak background ratio
come from that document. Per-base enrichment below was fixed before either arm
finished: analysis-BAM read records per reference base, signal over the merged
peak universe and background over the remaining `chrom.sizes` bases.

## Why the earlier smoke is not the retain arm

The 2026-10-28 smoke at workflow commit `c223dd9` used
`duplicate_analysis: auto`, which resolved to retain. Its identity was
`7e5599eebc879bce37fb2147b148f75996b8637bf5e88d846d84de3afd49cbe9`.
No `src/`, `workflow/`, `scripts/`, or `test/` file changed between that
commit and `#192` head `f93bf0c`. The smoke launch did not pass checkout
provenance, and its process log did not record a Bowtie2 version line.

This pair therefore reran both arms in one environment with explicit
`duplicate_analysis: retain` and `exclude`. The smoke directory was left in
place and was not used as an arm. Upstream BAMs were not copied: the identity
guard covers the whole outdir, and there is no executed handoff that can stage
a MAPQ BAM without either forging a marker or adding a cache. Each arm has its
own outdir and identity.

## Fixed inputs

Both arms use the same samples, reference, blacklist, tools, profile, MAPQ,
organelle setting, `max_insert`, and MACS3 `--keep-dup all`. Official FASTQ
MD5 values were recomputed and matched before the runs.

| Item | Value |
| --- | --- |
| Treatment | SRR12842272 R1/R2, sample `T2`, CUT&Tag H3K27me3, PE |
| Control | SRR12842282 R1/R2, sample `IgG2` |
| Treatment R1 SHA256 | `3080d1b6087d7d16c68b4ba6287ca3c53fcb2859fa93be3adc86f583b1da698e` |
| Treatment R2 SHA256 | `99db128438bcd347f5eb0cfd3b3db716b56bd7cc634294b4055edde185fd4f79` |
| IgG R1 SHA256 | `7e677fef03e6261300084d3dfd507436ad4ad539d15a3d0951d19e035b0fdea5` |
| IgG R2 SHA256 | `c3cd26ed7105a56cd3e99def1e9948bdc1c7152b584eb874636eced7e43f512e` |
| Reference | GENCODE GRCh38 primary assembly already admitted for the smoke; compressed MD5 `9c3fc2ca260a767530dddb0f26721a6b` |
| Blacklist | ENCFF356LFX at the smoke path, unchanged |
| Profile | `updated`, policy revision 6, schema 2.0.0 |
| MAPQ | 30 |
| `max_insert` | 700 on both arms |
| Organelle | config `auto`, resolved `exclude` with `chrM` on both arms |
| MACS3 | `-f BAMPE -g hs -q 0.01 --keep-dup all -B` on both arms |
| Cores | 8 |

`#192` head and the Snakemake source SHA were
`f93bf0c110f8c0d859e959be27448e780a79e4c6`. The worktree also contained this
comparison script and its test; Snakemake does not import them. Provenance
`verify-checkout` passed for the task-local clone of the locked `ci-fast`
environment after installing this checkout. The workflow import was
`ENCODE-sci7-cuttag-dup/src/encode_pipeline/__init__.py`. Snakemake was that
environment's 8.30.0. Scientific binaries came from the existing chipseq and
Picard prefixes, not from a new download.

| Tool | Recorded identity |
| --- | --- |
| Snakemake | 8.30.0 |
| Bowtie2 | 2.5.5. The process printed `Failed to launch x86-64-v3 version, staying with default` and then used `bowtie2-align-s` 2.5.5 |
| samtools | 1.23.1 |
| MACS3 | 3.0.4 |
| Picard | 3.5.0 |
| bedtools | 2.31.1 |
| FastQC | 0.12.1 |
| Trim Galore | 2.2.0 |
| bamCoverage | 3.5.6 |

## Identities

| Arm | `duplicate_analysis` | Identity |
| --- | --- | --- |
| retain | retain for T2 and IgG2 | `baf233c5c47b1ab001ad2fcb0c2c06709f830f879025c8d2ee9c4fd3518410b6` |
| exclude | exclude for T2 and IgG2 | `19a4943328ac9d5d6e28b6a49eb88ba644c85f74e61ed5446ce589357377be60` |

Retain did not pass `--exclude-dup` into analysis-BAM construction
(`False == True`). Exclude did (`True == True`). Both runs finished 48/48
steps with exit 0.

## Stage-labelled counts

NRF, PBC, and Picard percent duplication are computed before duplicate
removal. They match across arms because that stage is upstream of the only
policy difference. They are not an effect of retain versus exclude.

| Metric | Stage and denominator | retain | exclude |
| --- | --- | ---: | ---: |
| Picard percent duplication | pre-duplicate MarkDuplicates metrics; fraction of examined pairs | 0.008728 | 0.008728 |
| NRF | pre-duplicate MAPQ BAM; distinct fragments / total fragments | 0.991272 | 0.991272 |
| PBC1 | same stage | 0.991264 | 0.991264 |
| PBC2 | same stage | 114.365519 | 114.365519 |
| Analysis valid templates | analysis BAM after the duplicate policy, organelle filter, and blacklist | 4,144,663 | 4,108,514 |
| Analysis read records | `samtools view -c` on that analysis BAM, the FRiP denominator | 8,289,326 | 8,217,028 |

The analysis-BAM difference is 36,149 templates and 72,298 reads, exactly two
reads per template. Invalid analysis templates stayed 0 on both arms.

## Peaks

Own-peak FRiP uses each arm's own blacklist-filtered narrowPeak file. Those
two fractions do not share a universe, so they are not compared as an effect
size.

| | retain | exclude |
| --- | ---: | ---: |
| Narrow peaks | 12,458 | 2,881 |
| Mean width (bp) | 720.94 | 658.51 |
| Median width (bp) | 560 | 535 |
| Own-peak FRiP | 0.043496 | 0.012022 |
| Own-peak reads in peaks | 360,552 | 98,783 |

Overlap is the fraction of original peaks with at least one shared base, not
a base Jaccard.

| Direction | Fraction |
| --- | ---: |
| retain peaks touched by an exclude peak | 0.219297 |
| exclude peaks touched by a retain peak | 0.998959 |

Merged-base Jaccard uses bedtools-style merging of overlaps and book-ended
intervals, then intersection bases over union bases.

| | bases |
| --- | ---: |
| retain merged | 8,981,476 |
| exclude merged | 1,897,161 |
| intersection | 1,893,095 |
| union | 8,985,542 |
| Jaccard | 0.210682 |

The common universe is that union. It is only 4,066 bp larger than the merged
retain set, so a FRiP on this universe mostly asks whether exclude's analysis
BAM still covers the retain peaks. That is the pre-registered universe, not a
second criterion chosen after the counts.

## Common-universe signal

Denominator for both FRiP values: analysis-BAM read records. Numerator:
`scripts/calc_frip.py` (`bedtools intersect -u` and `samtools view -c`) against
the merged universe. Reference arm is retain. Epsilon is `1e-6`.

| | retain | exclude |
| --- | ---: | ---: |
| Reads in the common universe | 360,804 | 355,806 |
| FRiP | 0.043526 | 0.043301 |
| Off-peak fraction (`1 - FRiP`) | 0.956474 | 0.956699 |

Exclude minus retain is -0.0225 percentage points. The relative FRiP change is
-0.005169. `log2((0.043301+1e-6)/(0.043526+1e-6))` is -0.007477. The off-peak
background ratio, exclude over retain, is 1.000235.

Per-base enrichment uses the verified reference `chrom.sizes` total,
3,099,750,718 bp. The workflow's analysis-BAM rule did not receive a
`chrom_sizes` path; this denominator is the reference file, not a second
workflow output. Signal density is common-universe reads per universe base.
Background density is the other analysis-BAM reads per remaining reference
base.

| | retain | exclude |
| --- | ---: | ---: |
| Reads per peak base | 0.040154 | 0.039598 |
| Reads per off-peak base | 0.002565 | 0.002543 |
| Enrichment | 15.653 | 15.568 |

The pre-registered regression trigger is a FRiP drop of more than 10% relative
and an off-peak background increase of more than 5% relative, on every
independent replicate. This pair does not meet that trigger. It also cannot
satisfy the replicate-consistent peak or between-replicate clauses: there is
only one treatment/control pair. No default is frozen.

## Resources

`/usr/bin/time -v` measured the Snakemake parent, not the sum of child RSS.
`du` is final directory occupancy, not a scratch high-water mark. A
simultaneous scratch peak was not measured. `df` is the shared filesystem and
includes other work on the host.

| Arm | Wall | User CPU (s) | System CPU (s) | Average CPU | Parent max RSS | Final results |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| retain | 20:22.65 | 4903.68 | 171.68 | 415% | 3,480,576 KiB | 3.7G |
| exclude | 15:18.39 | 3742.11 | 147.35 | 423% | 3,480,576 KiB | 3.7G |

Free space on that filesystem moved from 478G to 467G between the retain start
and the exclude end. That delta is not attributed to this comparison alone.

## Large products

FASTQ, BAM, indexes, and run logs stay in the task directory
`/tmp/sci7-cuttag-dup-validation/` and are not in Git. SHA256:

| Product | SHA256 |
| --- | --- |
| retain T2 analysis BAM | `e77f960ada51b0200f9f4e499e2898bd711c069da7c7d048fe286f71d555b4d1` |
| retain IgG2 analysis BAM | `ca2ee416e9ec3d3e3186aec5df0ba16172ff50dd56c90d553b2d52e1c5008ff8` |
| exclude T2 analysis BAM | `1fe44b24644a4eee993470a21224ede20a9db25ab5e29d73baa683a4ccf06e22` |
| exclude IgG2 analysis BAM | `3b2b011a717b792e37d54ae59adb188679e9d8d60ff96083acbe19fb4ea22fab` |
| retain blacklist-filtered peaks | `3a2ca05f58084787167d92f1932539b373a810139d80d41b9996c450f6ddc7b1` |
| exclude blacklist-filtered peaks | `fb256cf45f9831943d436fcedec56686443fe2d892db0643f88a42978c4ee886` |

The machine-readable table is
[`cuttag-duplicate-retain-exclude.json`](cuttag-duplicate-retain-exclude.json).
Regenerate it with `scripts/compare_cuttag_duplicate_arms.py` from the two QC
summaries, the two blacklist-filtered narrowPeak files, the two treatment
analysis BAMs, and the reference `chrom.sizes`.

## Not measured

Length-resolved fragment coverage, a scratch high-water mark, the sum of child
RSS, and any second biological replicate. Confirmation data were not used.
