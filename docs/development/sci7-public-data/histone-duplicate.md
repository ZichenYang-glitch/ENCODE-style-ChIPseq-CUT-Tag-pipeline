# Histone H3K36me3 retain versus exclude

Status: **measured, not accepted.** One SE library was run twice, changing
only `duplicate_analysis`. No default is frozen. Confirmation data were not
used. `max_insert` is not applicable to SE.

The decision rules were already fixed in
[`comparison-design.md`](comparison-design.md). Peak count is descriptive.
Common-universe FRiP uses epsilon `1e-6`. The signal trigger is strict: a
relative FRiP drop greater than 10% and a relative off-peak rise greater
than 5%, together.

## Sample

ENCSR817LUF H3K36me3, matched control ENCSR603VTC. Seven FASTQ files,
2,510,695,554 official bytes. The treatment files are one library and one
technical replicate, not five biological replicates. Layout is SE. Read
length is 36 bp. Alignment is the existing GRCh38 primary index. The
workflow SHA is `57d8097eb96f431a543ac4d2fba94c513244fa7e`.

Both arms used the updated profile, policy revision 6, `organelle_mode`
resolved to `report`, and MACS3 `--keep-dup all` with `-f BAM`, `--broad`,
and `--broad-cutoff 0.1`. Picard `MarkDuplicates` used
`REMOVE_DUPLICATES=false`. The exclude arm then passed `--exclude-dup` when
building the analysis BAM. Configured `extend_reads` stayed `auto`.

| Tool | Version |
| --- | --- |
| Snakemake | 8.30.0 |
| Bowtie2 | 2.5.5 |
| samtools | 1.23.1 |
| MACS3 | 3.0.4 |
| Picard | 3.5.0 |
| bedtools | 2.31.1 |

## Arms

| Arm | Identity | Exit | Wall |
| --- | --- | --- | --- |
| retain | `a6746190eea7d85f1bff281ce1f03b8ad8211298ce3edb3feab76dc45ed2cb85` | 0, 41/41 | 1283 s |
| exclude | `b19ea352f9b37978834083a6a61771526bf18421046312fb4763711d66dd5826` | 0, 41/41 | 1309 s |

Sorted and MAPQ30 read counts match across arms. H3K36me3 has 38,213,007
sorted reads and 30,224,373 MAPQ30 reads. Input has 17,641,204 sorted reads
and 14,380,500 MAPQ30 reads. Pre-duplicate Picard `percent_duplication` is
0.112055 for H3K36me3. `estimated_library_size` is empty in that metrics
file, not a biological zero.

Exclude removes 3,353,662 H3K36me3 analysis reads and 213,495 Input analysis
reads. Each difference equals the duplicate-flag count in the retain
analysis BAM. The exclude analysis BAMs contain no duplicate-flagged
records.

## Effect

Blacklist-filtered broad peaks: retain 109,311, exclude 99,676. Merged-base
Jaccard is 0.814328 (intersection 56,576,821 bp, union 69,476,708 bp).
Common-universe FRiP is 0.171103 on retain and 0.167317 on exclude.
Relative FRiP change is −0.022131. Off-peak fraction ratio is 1.004568.
The joint signal trigger is false.

MACS3 model fragment length is 238 bp on retain and 240 bp on exclude, so
the auto treatment BigWig extension differs by 2 bp. Input stays on the
200 bp fallback in both arms. That extension difference is a consequence of
the duplicate filter, not a second configured variable.

## Verdict

No default is frozen. One library cannot support a replicate-concordance
claim, and confirmation has not been evaluated. The candidate exclude
policy does not meet the pre-registered signal-regression trigger on this
dataset.
