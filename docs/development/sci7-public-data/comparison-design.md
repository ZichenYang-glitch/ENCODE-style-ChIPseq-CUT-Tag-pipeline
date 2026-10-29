# SCI-7 single-variable comparison design

Status: execution preparation, plus the CUT&Tag measurements recorded in
[`cuttag-duplicate-retain-exclude.md`](cuttag-duplicate-retain-exclude.md)
and [`cuttag-replicate-insert.md`](cuttag-replicate-insert.md). Those
measurements do not freeze a default. The decision rules below are
unchanged.

## 1. Frozen inputs and evidence identity

Each comparison uses the file-level rows in
[`file-manifest.tsv`](file-manifest.tsv), the reference identity in
[`reference-identity.md`](reference-identity.md), the same aligner/MAPQ/proper-pair,
organelle, blacklist, peak-mode and threshold settings, and the same pinned tool
versions. Record the exact SCI-6+ workflow SHA, effective configuration, resolved
per-sample policies, commands, input/output checksums, CPU/RAM and wall time.
Use a new outdir for every scientific configuration; never rewrite an existing
run identity.

Full-depth libraries are required for default-changing biological conclusions.
If a trial is downsampled, sample read pairs with one seed and apply the identical
record selection to both mates and to every arm. A downsampled run is a controlled
execution check only: it changes duplicate collisions and peak statistics and
cannot freeze a duplicate or peak-calling default.

Discovery data may be used to investigate and pre-register attribution checks.
Confirmation data are evaluated once, after the rules are fixed, and are never
used to choose a value. Accession separation alone is insufficient; the
independence review in [`confirmation-review.md`](confirmation-review.md) must
exclude shared libraries and identical source files/checksums.

## 2. Experiment matrix

`auto` is a configuration-resolution check, not an additional biological arm.
When `auto` resolves to the same value as an explicit arm, do not repeat the
expensive scientific computation and count it as independent evidence. Compare
public artifacts after excluding run identity, logs and other legitimately
metadata-bearing files.

| Data group | Duplicate-policy arms (only this setting changes) | PE `max_insert` arms (only this setting changes) | `auto` resolution check |
| --- | --- | --- | --- |
| TF ChIP PE | `exclude` vs `retain` | 500 vs 2000 | `auto` must equal `exclude`; max_insert must equal 2000 |
| Histone ChIP SE | `exclude` vs `retain` | Not applicable; no `-X` is rendered | `auto` must equal `exclude`; SE records `max_insert` as not applicable |
| ATAC PE | `exclude` vs `retain` | 500 vs 2000 | `auto` must equal `exclude`; max_insert must equal 2000 |
| CUT&Tag PE | `retain` vs `exclude` | 500 vs 700 | `auto` must equal `retain`; max_insert must equal 700 |
| MNase PE | `exclude` vs `retain` | 200 vs 500 is a sensitivity diagnostic only; 500 remains the candidate default | `auto` must equal `exclude`; max_insert must equal 500 |

For duplicate policy, change only `duplicate_analysis` through the production
configuration/validator. `remove_dup` is tested only as the supported legacy
alias; do not run an alias arm that resolves to the same policy and call it an
independent experiment. Keep MAPQ, blacklist and every other setting fixed.
MNase 200 vs 500 measures the effect on the supported mono/di windows and is not
a proposal to adopt 200.

Run the single-variable comparisons before one final combined regression. The
combined run is an integration check and is not attributed to an individual
variable.

## 3. MACS3 keep-dup is a separate controlled tool experiment

Do not inject `--keep-dup` through `extra_args` and do not weaken or bypass the
production validator. The production workflow remains `--keep-dup all` and the
BAM layer owns duplicate policy.

For this methodological comparison only, freeze one production analysis BAM and
its matched control (when applicable), then invoke the pinned MACS3 executable
directly in an isolated evidence directory:

1. arm `all`: explicit `--keep-dup all`;
2. arm `one`: explicit `--keep-dup 1`, representing the previous implicit
   coordinate-collapse behavior.

Everything else, including `-f BAM/BAMPE`, genome size, peak mode, q/broad
threshold and model/shift settings, is identical. Save the complete commands,
MACS3 version, BAM/checksum identity and logs. This direct-tool experiment is not
a production configuration acceptance test and does not authorize a new user
switch. Repeat for the ordinary per-sample case first; pooled/consensus cases
are integration checks where enough independent replicates exist.

## 4. Metrics and effect sizes

Compute all cross-arm signal metrics on a common evaluation peak universe: merge
the peaks from the two arms of that same dataset, then use that frozen universe
for both arms. Do not calculate one arm's FRiP against its own peaks and compare
it directly with another arm using different peaks.

- **Pairing/template accounting:** concordant pair rate, valid template count,
  invalid/orphan/cross-contig/zero-TLEN counts, and count conservation at every
  stage. For max_insert, separately report newly admitted templates by length,
  MAPQ and coordinate/mate changes.
- **Short-template stability:** among MAPQ >= 30 templates with length <= 500 bp
  in the 500 arm, report the fraction whose R1/R2 coordinates or mate assignment
  changes in the comparison arm.
- **Duplicate/library metrics:** duplicate fraction, NRF/PBC and estimated
  library size at their specified pre-duplicate stages; do not recompute these
  from an already deduplicated BAM.
- **Signal/background:** FRiP on the common peak universe, off-peak normalized
  coverage, peak width and length-resolved fragment coverage. Report
  `log2((FRiP_arm + 1e-6)/(FRiP_reference + 1e-6))`, absolute FRiP
  difference in percentage points, and off-peak background ratio; FRiP is a
  fraction in [0, 1] and the epsilon is fixed at `1e-6` for every arm.
- **Peak agreement:** per-dataset peak-set Jaccard and the fraction of the
  reference arm's replicate-consistent peaks recovered by the other arm. Peak
  count is descriptive only.
- **Replicate consistency:** for assays with at least two independent biological
  replicates, compute the within-arm between-replicate peak Jaccard and
  length-binned coverage correlation. A single replicate has no such evidence.
- **Cost:** wall time, maximum RSS and CPU time under identical cores and input.

Report effect sizes per library and per independent replicate, not only pooled
counts. With the currently small candidate set, uncertainty is shown by the
range across replicates and datasets; do not claim population-level statistical
significance.

## 5. Pre-registered decisions

These are operational review triggers for this validation, not universal ENCODE
quality thresholds. No default is frozen from a single replicate, a single peak
count, or discovery data alone.

A candidate may be frozen only when all applicable hard checks pass, the
independent confirmation set has been evaluated once with these same rules, and
none of the following regression conditions is met:

1. **Hard failure in any arm:** configuration override does not take effect,
   tool failure, count non-conservation, invalid-template leakage, SE `-X`
   rendering, incorrect QC stage, or reference/identity mismatch.
2. **Signal/background regression:** the candidate arm both lowers FRiP by more
   than 10% relative and increases off-peak background by more than 5% relative
   on the common evaluation universe in every available independent replicate
   of a dataset.
3. **Loss of established signal:** more than 10% of the reference arm's
   replicate-consistent peaks are lost in every independent replicate, without
   a documented technical explanation consistent with the intended filter.
4. **Replicate regression:** within-arm between-replicate peak Jaccard falls by
   more than 0.10 absolute, or length-binned coverage correlation shows a clear
   contradictory direction between discovery and confirmation.
5. **For max_insert specifically:** more than 1% of the high-confidence short
   templates defined above change coordinates/mate assignment, or newly
   admitted long templates are predominantly low-MAPQ/repetitive and coincide
   with condition 2 or 3. Larger `-X` is not considered better merely because
   more pairs are called proper.

Meeting a trigger means do not freeze that default; retain all outputs and run a
pre-specified attribution comparison (read length, duplication, reference/repeat
context, and, for CUT&Tag, enrichment class). It does not automatically prove
the alternative default is universally correct.

## 6. No-conclusion conditions

Mark the affected dataset or whole comparison inconclusive when there is a
mismatched control, shared discovery/confirmation library or source checksum,
unverified reference identity, tool-version drift, insufficient depth after
downsampling, SE/PE mismatch, variable-read-length behavior that prevents a
controlled comparison, or phantom-dominated correlation used as quality
evidence. Discordant independent replicates are also inconclusive until their
technical cause is resolved. Missing mono/di coverage in MNase and missing
independent replicates in ATAC/CUT&Tag limit the stated scope rather than being
silently treated as success.

## 7. Evidence package

For every arm retain the manifest subset, reference binding, exact SHA/tool
versions, resolved policy, complete command, run record, JUnit where applicable,
stage counts, metric tables with effect sizes, logs, and output checksums. The
final review must distinguish controlled implementation evidence from biological
support and list every assay/condition that remains unvalidated.
