# SCI-7 batch summary

Status: **no default is frozen.** CUT&Tag, histone H3K36me3, MNase Rep1,
TF NR2F2, and ATAC Treg discovery have each been measured under the
pre-registered rules in [`comparison-design.md`](comparison-design.md).
The TF confirmation pair, JUND ENCSR785RQR with control ENCSR194GTV, has
also been measured and is in
[`tf-jund-confirmation.md`](tf-jund-confirmation.md). It does not freeze a
default. Exit status 0 is not scientific acceptance.

This note records effect sizes, independence, and the confirmation
candidates that remain unevaluated. It does not change a production default.

## What was measured

| Batch | Comparison | Independence | Result |
| --- | --- | --- | --- |
| CUT&Tag H3K27me3 | retain versus exclude; max_insert 500 versus 700 | Rep1 is shallow; Rep1 and Rep2 are submitter labels, not established independent biological replicates | Inconclusive. Do not freeze `retain` or `max_insert` 700. |
| Histone H3K36me3 | retain versus exclude. SE, so `max_insert` does not apply | One human library, one technical replicate | No default. Joint signal trigger is false. |
| MNase Rep1 | exclude versus retain; 200 versus 500 is a diagnostic only | One library. Rep2 and Rep3 were not run | No conclusion. The registered 1% short-template metric was not evaluated. |
| TF NR2F2 | exclude versus retain; max_insert 500 versus 2000 | One biological replicate and its matched input | No conclusion. Joint signal trigger is false. Short-template fraction is far below 1%. |
| ATAC Treg | exclude versus retain; max_insert 500 versus 2000 | Two biological replicates, one donor, ENCDO236IKC. Reads are 26 bp | No conclusion. Neither replicate meets the joint signal trigger. Short-template fractions are about 0.00115 and 0.00118. |
| TF JUND confirmation | exclude versus retain; max_insert 500 versus 2000 | One JUND library in K562 cells with a CEBPB deletion, plus matched input. Different target and cell context from NR2F2 | No conclusion. Joint signal trigger is false. Short-template fraction is 1,829 / 155,206,704. |

Candidate settings that were tested, and not frozen:

- ChIP and ATAC: `duplicate_analysis=exclude`, `max_insert=2000`.
- CUT&Tag: `duplicate_analysis=retain`, `max_insert=700`.
- MNase: `duplicate_analysis=exclude`, `max_insert=500`. The 200 arm is not a candidate default.
- Histone SE: `duplicate_analysis=exclude`. No `-X`.

## Effect sizes that matter for the rules

CUT&Tag and histone details stay in their own notes. The new batches:

- TF, one library. Blacklist-filtered narrow peaks are 114,667 on retain and
  88 on exclude. Common-universe relative FRiP change is −0.378566. Relative
  off-peak change is 0.044935, which is not strictly above 0.05, so the joint
  trigger is false. Short-template fraction is 890 / 13,023,039.
- ATAC, both biological replicates. Merged-base Jaccard of retain versus
  exclude is 0.947501 and 0.943664. Common-universe relative FRiP change is
  −0.022815 and −0.030187. Relative off-peak change is 0.023723 and 0.031972.
  The joint trigger is false in both replicates. Short-template fractions are
  30,475 / 26,570,383 and 51,984 / 44,116,609. Within-arm between-replicate
  Jaccard does not fall by more than 0.10 on either contrast. Replicate-consistent
  retain peaks lost by exclude are 1.5% and 0.9%.

A large TF peak-count change and a small ATAC peak-set change are both
descriptive. Neither freezes a default by itself.

The JUND confirmation library is also one replicate. Blacklist-filtered
narrow peaks are 316,361 on retain and 7,430 on exclude. Common-universe
relative FRiP change is −0.126901. Relative off-peak change is 0.034947,
which is not strictly above 0.05, so the joint trigger is false.
Short-template fraction is 1,829 / 155,206,704.

## Independence

Donor, library, and biosample overlap was reviewed before any confirmation
download. TF discovery and its named confirmation candidate do not share a
donor. Histone confirmation is mouse H3K4me3 and does not repeat human
H3K36me3. ATAC discovery and ENCSR639PPN share donor ENCDO236IKC, so that
pair cannot be called donor-independent. Two ATAC biological replicates in
the discovery experiment are still one donor.

## Confirmation gate

The questions locked by these measurements are the pre-registered ones:
whether `exclude` versus `retain`, and whether the candidate `max_insert`,
regress signal, short templates, or replicate concordance under the published
thresholds. Discovery alone cannot answer them for a default.

Candidates that do not answer those questions were not downloaded:

- Histone ENCSR213YKJ is a different species and a different histone mark
  from the completed H3K36me3 result.
- ATAC ENCSR639PPN shares donor ENCDO236IKC with discovery. It is not
  donor-independent confirmation.
- Optional CUT&Tag and MNase confirmation sets were not required to interpret
  the completed no-conclusion records, and they were not used to tune a
  threshold.

TF ENCSR785RQR with control ENCSR194GTV was the donor-independent
confirmation candidate named for the TF policy question. It is a different
target and cell context from NR2F2. Its four FASTQ files are 74,110,406,175
official bytes. Those files were downloaded, the three arms finished, and
the registered comparison is in
[`tf-jund-confirmation.md`](tf-jund-confirmation.md). The result does not
freeze `exclude` or `max_insert` 2000, and it was not used to move a
threshold.

## Unevaluated metrics

- Confirmation for CUT&Tag, histone, MNase, and ATAC.
- MNase registered short-template fraction on the 500-arm denominator.
  A separate 200-arm diagnostic was reported and is not that metric.
- Length-binned coverage correlation as quality evidence. Occupied-window
  correlations were computed for ATAC and were not treated as acceptance.
  No discovery-versus-confirmation direction was available.
- Donor-independent ATAC replication.
- Same-target replication of NR2F2 or H3K36me3.
- Any claim that 26 bp ATAC reads behave like the 100 bp TF libraries.

## Conclusion

No production default is frozen. The measured batches do not agree with one
another strongly enough, and are not independent enough, to promote
`exclude`, `retain`, or a PE `max_insert` value beyond the assay-specific
candidate already written in the comparison design. Inconclusive remains a
valid result.
