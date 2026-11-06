# SCI-7 resource budget

Status: execution planning on 2026-10-28. No large FASTQ/BAM, FASTA or index
has been downloaded. Official byte counts are known from portal metadata/HEAD.
BAM sizes, sort temporary space, index size and downstream output sizes are
**unmeasured**; no multiplier in this document is an admission requirement.

The previous approximately 1.2 TiB figure was only an all-input/all-branch
retention scenario using an unmeasured BAM-to-FASTQ multiplier and treating all
28 configurations as independent full chains. It is **not** a minimum execution
requirement and must not be used to request expansion or stop the minimal
trial.

## 1. Reuse and recomputation from the actual workflow

The run matrix still has 28 logical configurations, but it does not require 28
independent alignments.

| Stage/comparison | Reuse decision | What must run again |
| --- | --- | --- |
| Reference FASTA/FAI/chrom_sizes/blacklist/Bowtie2 index | Reuse read-only across every arm and cohort after identity/admission checks | Rebuild only when the reference identity changes |
| Baseline alignment used by duplicate and max_insert comparisons | Reuse the baseline run's artifacts; do not align once per comparison | Downstream products are produced under each isolated run identity |
| Duplicate policy in **updated** profile | `duplicate_handling` always uses Picard `REMOVE_DUPLICATES=false` (`mark`), so for fixed MAPQ input its `final.bam` and duplicate metrics do not depend on retain/exclude | `prepare_analysis_bam` must run again because `--exclude-dup` changes; all downstream analysis/coverage/peaks/QC must be regenerated |
| Duplicate policy in **legacy** profile | Reuse only the upstream MAPQ BAM/BAI. The legacy `remove_dup=yes` branch physically removes duplicates while retain symlinks the input, so `final.bam` is not equivalent | `duplicate_handling` and every downstream stage must run again |
| PE `max_insert` | No aligned BAM is reusable when `-X` changes | Bowtie2 alignment, MAPQ/final/analysis BAMs and all downstream products must run again |
| MACS3 keep-dup | Reuse exactly the same frozen production analysis BAM and matched control; no new alignment or BAM | Only the isolated direct-MACS3 calls, peaks and logs in a separate evidence directory |
| `auto` resolution | No additional scientific arm when auto resolves to an existing explicit value | Only configuration/dry-run/effective-policy checks and artifact comparisons |

A reusable artifact does **not** permit sharing an outdir or bypassing the
identity guard. Every scientific configuration gets a new outdir and its own
`analysis_identity.json`/`effective_config.json`. For updated duplicate-policy
reuse, first create the new run marker through the production entry, then stage
the verified MAPQ/final BAM under a written handoff manifest containing source
SHA, destination SHA, file SHA256 and the scientific equivalence argument
(fixed MAPQ and unchanged Picard `mark`). Never reuse the other run's marker or
write into an unmarked nonempty outdir. Legacy reuse follows the same procedure
for the MAPQ BAM only.

Thus one PE dataset needs one baseline alignment and one alternative-max_insert
alignment; it does not need a third alignment merely to compare duplicate
policy. Histone SE needs no alternative-max_insert alignment. Across all 28
logical configurations there are at most 18 core alignments (9 discovery and 9
confirmation), plus legacy duplicate-handling/downstream reruns.

## 2. Known input and reference downloads

Do not download all registered inputs at once. The totals below describe the
registered manifest, not a required initial download.

| Scope | FASTQ | Official compressed bytes | GiB |
| --- | ---: | ---: | ---: |
| Minimal smoke: CUT&Tag Rep2 treatment + IgG Rep2 | 4 | 383,633,518 | 0.357 |
| CUT&Tag discovery, whole group | 8 | 443,578,936 | 0.413 |
| Histone discovery, whole group | 7 | 2,510,695,554 | 2.338 |
| ATAC discovery, whole group | 4 | 3,952,560,652 | 3.681 |
| MNase discovery, whole group | 6 | 4,402,061,533 | 4.100 |
| TF discovery treatment + control | 4 | 11,396,944,028 | 10.613 |
| **All discovery, only after sequential expansion** | **29** | **22,705,840,703** | **21.146** |
| Optional CUT&Tag confirmation | 6 | 232,017,631 | 0.216 |
| Optional MNase confirmation | 4 | 3,139,436,756 | 2.924 |
| Histone confirmation rep1 + control | 16 | 7,462,062,697 | 6.949 |
| ATAC confirmation | 8 | 11,262,765,629 | 10.489 |
| **TF confirmation treatment + control** | **4** | **74,110,406,175** | **69.021** |
| All confirmation candidates, not an initial download | 38 | 96,206,688,888 | 89.599 |

For human smoke, the GRCh38 primary FASTA is 844,691,642 bytes and the already
verified ENCFF356LFX blacklist is 8,211 bytes. GTF is not needed. If the FASTA
is not already available, the smoke's newly authorized download would be
**1,228,325,160 bytes (1.144 GiB)** for four FASTQ plus FASTA; the local
blacklist need not be fetched again. The Bowtie2 index is built locally, not
downloaded.

Both species' compressed FASTAs total 1.507 GiB. Decompressed FASTA and index
sizes are unmeasured. Keep the shared reference working set separate from
per-run scratch and reuse it across groups.

## 3. Scenario budgets: known and unmeasured components

For every scenario:

`peak = raw inputs currently retained + shared reference working set + long-term products + current-group scratch + parallel-extra`

Only raw input and compressed reference sizes are currently known. The other
terms remain explicitly unmeasured until the smoke and first group report them.

### 3.1 Minimal smoke

- Known new download: 1.144 GiB including GRCh38 FASTA, or 0.357 GiB if the
  reference already exists locally.
- Computation: one default CUT&Tag Rep2 + IgG Rep2 run (auto resolves to
  retain, max_insert 700), then the two isolated MACS3 keep-dup calls on the
  same frozen analysis BAM. No alternative-max_insert alignment and no full
  biological comparison.
- Long-term products: effective config/identity, final and analysis BAMs
  selected for the smoke evidence, duplicate/stage metrics, peaks, QC, logs,
  commands and checksums. Sizes unmeasured.
- Current scratch: one alignment and its sort/index temporary files, one
  analysis-BAM derivation and small MACS3 outputs. Sizes unmeasured.
- Parallel-extra: none; the smoke is serial.

### 3.2 One discovery group

Run one group at a time. Raw input is the corresponding table value; the
largest discovery group is TF at 10.613 GiB. For a PE group, scratch covers at
most the baseline and alternative-max_insert alignments plus the alternative
duplicate analysis/downstream products, not three complete alignments. Histone
SE has only the baseline alignment. Long-term products accumulate only after
that group's evidence is accepted.

### 3.3 All discovery, sequentially

The eventual raw total is 21.146 GiB, but downloads occur per group. At any
time, scratch is only the active group's scratch; it is not the sum of all 14
logical configurations. Long-term products from completed groups accumulate
and must be measured. After each group, recompute the next group's projection
using observed sizes instead of assuming a BAM/FASTQ ratio.

### 3.4 Confirmation by group, with TF separate

Confirmation starts only after discovery decision rules are fixed. The same
reuse rules apply. The largest single group is TF confirmation at 69.021 GiB;
its raw download, two PE alignments and duplicate-analysis/downstream scratch
must be budgeted as one active group. It is not added to a simultaneously
retained 89.599 GiB of other confirmation scratch. Download TF separately and
pause after download, baseline alignment and alternative-max_insert alignment
to remeasure free space.

### 3.5 Optional all-intermediates-retained scenario

The former approximately 1.2 TiB scenario is retained only as an **unmeasured,
non-required upper-retention illustration** based on the old 1.5–3x multiplier.
It is not reliable enough to be a new estimate and is not a prerequisite for
any stage. If a future reviewer explicitly requires all branch intermediates to
remain, calculate a fresh requirement from measured per-stage sizes.

## 4. Long-term retention and task-owned scratch

Long-term retention includes original FASTQ and reference inputs, official and
local checksums, source manifests/query evidence, effective configs and run
identities, selected final/analysis BAMs, duplicate/stage metrics, final peaks,
QC, coverage, logs, commands, effect-size tables and the frozen BAM identity
used for MACS3 keep-dup. These are not cleanup candidates.

Alternative-arm scratch is limited to task-created MAPQ/final/analysis BAMs,
sort temporary files and redundant intermediate coverage whose final metrics
and checksums have already been preserved. Any future cleanup requires the
existing execution authorization, an explicit file list, verification that the
outputs are task-owned and confirmation that long-term evidence is intact.
No raw file, user result, source evidence or traceability record is deleted to
make the budget fit.

Each additional concurrent alignment adds its own MAPQ/final/sort temporary
space and any simultaneous download. Shared reference/index files are not
duplicated, but concurrent downstream scratch and retained outputs still add
space. Parallel-extra is unmeasured and must use the largest observed active
per-alignment scratch, not zero and not a fixed multiplier.

## 5. Measurement and stop rules

At reference decompression/index build, Bowtie2/MAPQ, duplicate handling,
analysis-BAM derivation, peaks/QC and each alternative arm, record:

- `df --output=avail` before and after the stage;
- allocated bytes (`du -s --block-size=1`) for reference and each outdir;
- individual BAM/index/output sizes and the largest temporary files;
- a high-water scratch observation during alignment/sort, since final `du`
  can miss removed temporary files;
- compressed-input to measured-stage-space ratios, retained by assay/layout.

After the smoke, use the largest observed ratio to project the next group.
After every group, replace the projection with the new maximum observed value.

Stop and reassess rather than silently continuing when:

1. available space is below the known next download plus measured retained and
   active working sets plus a **20 GiB operational safety margin**;
2. before any measurement exists, the smoke cannot start with at least 20 GiB
   free, its scratch reaches 10 GiB, or free space falls below 5 GiB. These
   are initial guardrails, not claims about actual usage;
3. the TF confirmation projection using measured alignment and downstream
   ratios does not fit after the 20 GiB margin; pause before its download or
   current alignment stage instead of invoking the obsolete 1.2 TiB figure;
4. checksum/reference admission fails, or a stage's actual scratch exceeds the
   current approved projection without an explained difference.

These rules neither declare the currently available space insufficient nor
promise that the complete matrix will fit. They make the decision after each
measurement.

## 6. Execution order

1. Run the minimal smoke and measure every stage.
2. For CUT&Tag Rep2, add duplicate exclude using the updated-profile final-BAM
   handoff (no new alignment), then max_insert 500 (one new alignment). Expand
   to the full CUT&Tag group and explicitly assess shallow Rep1.
3. Add MNase Rep1, then histone, ATAC and the remaining discovery groups,
   measuring and recomputing storage before each expansion.
4. Run optional CUT&Tag/MNase confirmation, histone and ATAC confirmation.
   Evaluate the 69.021 GiB TF confirmation last as its own storage gate.
5. Keep all outputs and evidence from completed groups; apply only the
   task-owned scratch policy above under the applicable authorization.

No stage in this preparation document is itself authorization for a large
download or a real-library matrix run.

## Measured minimal-smoke checkpoints (2026-10-28)

These are measured final/point values, not scratch high-water marks:

- Final reference directory: **7.7 GiB** (FASTA, FAI, chrom_sizes, six `.bt2`).
- Final results directory: **3.7 GiB**.
- Downloaded FASTQ + GRCh38 FASTA: **1,228,325,160 bytes**.
- Production run wall time: **19:11.70**; user CPU **4,784.74 s**;
  system CPU **183.61 s**; 431% average CPU.
- `/usr/bin/time` max RSS for the timed Snakemake process:
  **3,483,136 KiB**. This is not the simultaneous total RSS of all children.
- Bowtie2-build `/usr/bin/time` max RSS: **7,305,712 KiB**.
- Simultaneous scratch peak: **not measured**. Start/end `df` samples include
  other host tasks and are not attributable to this run.
- The earlier 1.2 TiB all-input/all-intermediates-retained figure is not a
  startup requirement and is not used here. The measured final occupancy is
  about 12.6 GiB for reference + results + inputs; scratch must be measured
  before larger runs.
