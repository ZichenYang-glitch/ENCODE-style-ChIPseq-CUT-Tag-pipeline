# TF NR2F2 discovery: measured, no conclusion

Status: **no conclusion.** ENCSR054ZMK and matched control ENCSR142JFU were
downloaded, the four FASTQ files matched the official sizes and MD5 values,
and the three pre-registered arms each finished with exit status 0. No
default is frozen. The decision rules in
[`comparison-design.md`](comparison-design.md) are unchanged. One biological
replicate cannot freeze a default, and confirmation has not been evaluated.
Exit status 0 is not scientific acceptance.

The candidate settings are `duplicate_analysis=exclude` and `max_insert=2000`.
`auto` was a dry-run only. It printed
`duplicate_analysis=auto(chipseq:exclude)` and `max_insert=auto(chipseq:2000)`
and was not a third alignment.

## Sample

ENCODE experiment ENCSR054ZMK, target NR2F2 in HepG2, control experiment
ENCSR142JFU. One biological replicate. Layout is PE. Read length is 100 bp.
The four official FASTQ files total 11,396,944,028 bytes. Alignment reused
the existing GRCh38 primary assembly index. The compressed FASTA is
844,691,642 bytes, MD5 `9c3fc2ca260a767530dddb0f26721a6b`, and that identity
was rechecked before alignment. The executed checkout is
`0ce4ec1ef46eb445a4ae0c2f5b06fb66b0b32be2`.

| File | Official bytes | Official MD5 |
| --- | ---: | --- |
| ENCFF032VKB (treatment R1) | 2,947,099,459 | `3418ca66ca583e45d9bff87ce8c6893d` |
| ENCFF996ZFG (treatment R2) | 3,025,085,699 | `bd36578c69a007effd8cec8b1db7673d` |
| ENCFF429FLQ (control R1) | 2,668,139,292 | `28424afcf1c19a6b2375132429e93603` |
| ENCFF225SRY (control R2) | 2,756,619,578 | `21938488725ce0dec60f03899009221d` |

All three real arms used the updated profile, schema 2.0.0, policy revision
6, `organelle_mode` resolved to `report`, blacklist mode `analysis`, and
MACS3 `-f BAMPE -g hs -q 0.01 --keep-dup all`. Narrow peak mode does not pass
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

## Arms

Each arm is a full workflow in its own output directory. There is no BAM
handoff. Only one configured variable changes from the candidate.

| Arm | Only change | Identity | Exit | Wall |
| --- | --- | --- | --- | ---: |
| exclude, max_insert 2000 | candidate | `8a71d6c8ac0a30652cb89cafb0965c292486b12f856000023a79f2c51439d5fb` | 0, 44/44 | 8850 s |
| retain, max_insert 2000 | `duplicate_analysis` | `5d4e4177199218d0482ea5a74ba3aac94245cd9789c42ff17d37c741b2ef5930` | 0, 45/45 | 9060 s |
| exclude, max_insert 500 | `max_insert` | `1e7c10d6db4e4f95642341ec5ac670aa0093d5fa06dd48e3b4170768e7ef6d08` | 0, 45/45 | 6810 s |

The candidate DAG scheduled 44 jobs. Its identity file was already on disk, so
`analysis_identity_marker` was not scheduled again. The other two arms
scheduled that marker and finished 45 jobs. All three identity files exist.

Pre-duplicate Picard metrics on the MAPQ BAM match between exclude/2000 and
retain/2000. NR2F2 `percent_duplication` is 0.601046
(7,821,755 / 13,013,573 read pairs). `estimated_library_size` is 5,849,048.
Input `percent_duplication` is 0.392451
(17,678,166 / 45,045,502 read pairs). NRF, PBC1, and PBC2 on NR2F2 are
0.398954, 0.359842, and 1.414671 in both of those arms. The 500 arm is a
different alignment, so its pre-duplicate duplication is 0.601058 and is
not required to match.

Analysis-BAM read counts, after duplicate policy and blacklist filtering:
exclude/2000 has 10,361,138 reads and 5,180,569 valid fragments; retain/2000
has 25,972,520 reads and 12,986,260 valid fragments; exclude/500 has
10,366,302 reads and 5,183,151 valid fragments.

## Duplicate effect

Blacklist-filtered narrow peaks: retain 114,667, exclude/2000 88. Merged-base
Jaccard, retain versus exclude, is 0.001190 (intersection 41,193 bp, union
34,616,904 bp). Every exclude peak base sits inside the retain merged set.
The fraction of retain peaks touched by an exclude peak is 0.000567.

Own-peak FRiP is descriptive. It uses each arm's own peaks: retain 0.106103
(2,755,756 / 25,972,520) and exclude/2000 0.000553 (5,730 / 10,361,138).

The registered signal check uses one common peak universe and each arm's own
analysis-BAM denominator. Common-universe FRiP is 0.106103 on retain and
0.065936 on exclude (683,174 / 10,361,138). Relative FRiP change is
−0.378566. Relative off-peak change is 0.044935. The joint trigger requires
both a relative FRiP change strictly below −0.10 and a relative off-peak
change strictly above 0.05. The FRiP side meets that inequality. The
off-peak side does not. The joint trigger is false.

## max_insert effect

The registered short-template denominator is MAPQ ≥ 30 and length ≤ 500 bp
in the 500 arm. The event is an R1 coordinate, R2 coordinate, or ordered
mate-contig change in the 2000 arm. A 500-arm template absent from the 2000
arm stays in the denominator and is not a coordinate change. The trigger is
a fraction strictly greater than 0.01.

Denominator 13,023,039. Changed 890. Fraction 0.00006834. The trigger is
false. Templates present only in the 2000 arm number 7,862. Of those, 7,577
have MAPQ < 30 and 5,777 have length > 700 bp. That count is not the
registered trigger.

Own-peak calls on exclude/500 are 88 peaks and FRiP 0.000553, the same peak
count and the same own-peak FRiP as exclude/2000. That similarity is
descriptive. It is not a second copy of the short-template test.

## Verdict

No default is frozen. The duplicate contrast changes the delivered peak set
by a large amount, and the common-universe FRiP drop is large, but the
pre-registered joint signal rule does not fire. The short-template fraction
is far below 1%. One library, without confirmation, cannot choose
`exclude` or `max_insert` 2000 as a default.
