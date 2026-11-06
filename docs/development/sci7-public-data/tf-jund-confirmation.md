# TF JUND confirmation: measured, no conclusion

Status: **no conclusion.** ENCSR785RQR and matched control ENCSR194GTV were
downloaded, the four FASTQ files matched the official sizes and MD5 values,
and the three pre-registered arms each finished with exit status 0. No
default is frozen. The decision rules in
[`comparison-design.md`](comparison-design.md) are unchanged. This is a
cross-target parameter check against the NR2F2 discovery result in
[`tf-nr2f2-discovery.md`](tf-nr2f2-discovery.md), not a same-target
replicate. One biological replicate cannot freeze a default. Exit status 0
is not scientific acceptance.

The candidate settings are `duplicate_analysis=exclude` and `max_insert=2000`.
`auto` was a dry-run only. It printed
`duplicate_analysis=auto(chipseq:exclude)` and `max_insert=auto(chipseq:2000)`
and was not a third alignment.

## Sample

ENCODE experiment ENCSR785RQR, target JUND in K562 cells treated with sgRNA
targeting CEBPB. Control experiment ENCSR194GTV is input from that same
CEBPB-deletion background. One biological replicate. Layout is PE. Read
length is 151 bp. The four official FASTQ files total 74,110,406,175 bytes.
There is no donor, biosample, library, or source-file overlap with discovery
ENCSR054ZMK and ENCSR142JFU. The target and cell context are different.

| File | Role | Official bytes | Official MD5 |
| --- | --- | ---: | --- |
| ENCFF861BMZ (R1) | JUND | 16,035,460,793 | `9331f325051d9231d43da09c796fee1a` |
| ENCFF919EHZ (R2) | JUND | 16,771,747,554 | `6d3130010f4ef05f4465ae5382bc8705` |
| ENCFF260LDK (R1) | Input | 20,309,688,431 | `afd4058b441e9caee9b4000096ab596a` |
| ENCFF977HJB (R2) | Input | 20,993,509,397 | `1ed068054107050e7d992103ef6e1ebf` |

Alignment used the GRCh38 primary-assembly index whose compressed FASTA is
844,691,642 bytes, official MD5 `9c3fc2ca260a767530dddb0f26721a6b`, SHA256
`e9a2d5a5cd225293646ae298998ad4ea8c11e8f22729bc2d6f5c3dcfefbf8ef8`. That
identity was rechecked before alignment. Blacklist ENCFF356LFX, compressed
MD5 `393688b4f06c9ce26165d47433dd8c37`, is the same blacklist used for the
NR2F2 discovery arms. The executed checkout is
`4ce8c6ee7384083b2b6a8e191dd8e383f76437fe`. Its workflow tree matches
`634160e0d7c671db466a22a24bbe2287782b9af7`.

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
| exclude, max_insert 2000 | candidate | `c9c3836f8ff55fc5a6bd4ea734b699b3953eb044f83ec716fb20b6eae2462d30` | 0, 45/45 | 80581 s |
| retain, max_insert 2000 | `duplicate_analysis` | `f46f39d699f117a046c8a08bb93833488a9bc39cb88e4015fad45b83a5468577` | 0, 45/45 | 63241 s |
| exclude, max_insert 500 | `max_insert` | `b497037b199f53a3233f9d9dff39e88b822cf2592fddeacef46b93790978aae5` | 0, 45/45 | 50160 s |

Pre-duplicate Picard metrics on the MAPQ BAM match between exclude/2000 and
retain/2000. JUND `percent_duplication` is 0.802631
(124,530,229 / 155,152,480 read pairs). `estimated_library_size` is
30,833,501. NRF, PBC1, and PBC2 are 0.197369, 0.132548, and 1.002091.
Input `percent_duplication` is 0.776677
(162,796,607 / 209,606,485 read pairs). `estimated_library_size` is
47,408,853. The 500 arm is a different alignment, so its JUND
pre-duplicate duplication is 0.802650 (124,576,596 / 155,206,704) and is
not required to match.

Analysis-BAM read counts, after duplicate policy and blacklist filtering:
exclude/2000 has 61,147,296 reads and 30,573,648 valid fragments; retain/2000
has 309,865,724 reads and 154,932,862 valid fragments; exclude/500 has
61,160,940 reads and 30,580,470 valid fragments.

## Duplicate effect

Blacklist-filtered narrow peaks: retain 316,361, exclude/2000 7,430.
Merged-base Jaccard, retain versus exclude, is 0.025191 (intersection
3,581,583 bp, union 142,179,377 bp). The fraction of retain peaks touched by
an exclude peak is 0.021684. The fraction of exclude peaks touched by a
retain peak is 0.999731.

Own-peak FRiP is descriptive. It uses each arm's own peaks: retain 0.215922
(66,906,750 / 309,865,724) and exclude/2000 0.016773
(1,025,642 / 61,147,296).

The registered signal check uses one common peak universe and each arm's own
analysis-BAM denominator. Common-universe FRiP is 0.215925 on retain
(66,907,892 / 309,865,724) and 0.188524 on exclude
(11,527,719 / 61,147,296). Relative FRiP change is −0.126901. The log2 ratio
with epsilon `1e-6` is −0.195781, and the absolute difference is −2.740100
percentage points. Relative off-peak change is 0.034947. The joint trigger
requires both a relative FRiP change strictly below −0.10 and a relative
off-peak change strictly above 0.05. The FRiP side meets that inequality.
The off-peak side does not. The joint trigger is false.

Per-base density on the same universe is descriptive. Retain has 0.470588
reads per peak base and 0.082148 reads per off-peak base. Exclude has
0.081079 and 0.016777. These densities are not the joint trigger.

With one biological replicate, within-arm between-replicate Jaccard and
replicate-consistent peak loss are not applicable.

## max_insert effect

The registered short-template denominator is MAPQ ≥ 30 and length ≤ 500 bp
in the 500 arm. MAPQ is the minimum of the two mates. The records are
primary proper pairs from the aligner sorted BAM. The event is an R1
coordinate, R2 coordinate, or ordered mate-contig change in the 2000 arm.
A 500-arm template absent from the 2000 arm stays in the denominator and is
not a coordinate change. The trigger is a fraction strictly greater than
0.01.

Denominator 155,206,704. Changed 1,829. Fraction 0.00001178. Eighteen
qualifying templates are absent from the 2000 arm; they remain in the
denominator and are not counted as changes. The trigger is false.

Templates present only in the 2000 arm number 27,858. Of those, 26,369 have
MAPQ < 30 and 21,877 have length > 700 bp. High-MAPQ templates in that set
number 1,489, of which 822 are 501–700 bp and 667 are longer than 700 bp.
That count is not the registered trigger, and it does not coincide with a
fired joint signal trigger.

Own-peak calls on exclude/500 are 7,425 peaks and own-peak FRiP 0.016773.
That similarity to exclude/2000 is descriptive. It is not a second copy of
the short-template test.

## Verdict

No default is frozen. Exclude lowers common-universe FRiP by more than 10%
relative, and the delivered peak set is much smaller, but the pre-registered
joint signal rule does not fire. The short-template fraction is far below
1%. The library is one biological replicate of JUND in a CEBPB-deletion
K562 background. It does not repeat NR2F2 in HepG2, and it was not used to
change a threshold. That is not enough to freeze `exclude` or `max_insert`
2000.
