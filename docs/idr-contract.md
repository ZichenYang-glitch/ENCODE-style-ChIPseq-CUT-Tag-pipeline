# IDR Contract (v0.2)

This document defines the TF ChIP-seq IDR (Irreproducible Discovery Rate) behavioral
contract as implemented in Stage 5. It describes actual outputs, eligibility, and
scope boundaries — no aspirational features.

## Master gate

`stage5: true` enables all IDR rules. Requires `stage4b: true` (enforced by validator).

## Eligibility

All of the following must be true:

| Condition | Requirement |
| :--- | :--- |
| `assay` | `chipseq` |
| `peak_mode` | `narrow` |
| Biological replicates | Exactly 2 treatment biological replicates per experiment |

Technical replicates are supported — they are merged into biological-replicate BAMs
by Stage 4b before IDR processing. The IDR eligibility check counts unique
`biological_replicate` values, not sample rows.

## Outputs

All IDR outputs land under `results/experiments/<experiment>/`.

### Stage 5a: True-replicate IDR

| Output | Path | Description |
| :--- | :--- | :--- |
| IDR-ready biorep peaks | `04_peaks/idr/<exp>_biorep<N>_idr_peaks.narrowPeak` | MACS3 per-biorep peak calls with relaxed `-p` threshold (`tool_parameters.idr_macs3.pvalue`, default 0.1) |
| True-replicate IDR raw | `06_idr/true_replicates/idr.txt` | Raw IDR output from `idr --samples` between the two biorep peak sets |
| True-replicate IDR thresholded | `06_idr/true_replicates/idr.thresholded.narrowPeak` | IDR-thresholded narrowPeak at configured `idr.threshold` (default 0.05) |

### Stage 5b: Pseudoreplicate IDR

| Output | Path | Description |
| :--- | :--- | :--- |
| Pseudorep BAMs | `05_pseudorep/<exp>_biorep<N>.pr1.bam`, `.pr2.bam` | Deterministic hash-based pseudoreplicate split of each biorep BAM and pooled BAM |
| Pseudorep MACS3 peaks | `04_peaks/idr/<exp>_<source>_pr<1|2>_idr_peaks.narrowPeak` | MACS3 peak calls on each pseudoreplicate BAM |
| Self-IDR per biorep | `06_idr/self_pseudoreplicates/biorep<N>.idr.txt` | Self-consistency IDR between pseudoreps of the same biorep |
| Pooled-IDR | `06_idr/pooled_pseudoreplicates/idr.txt` | IDR between pooled pseudoreps |

### Stage 5b: Final peak sets

| Output | Path | Description |
| :--- | :--- | :--- |
| Conservative peaks | `06_idr/final/conservative.narrowPeak` | Copied from the true-replicate thresholded IDR output (`idr.thresholded.narrowPeak`) |
| Optimal peaks | `06_idr/final/optimal.narrowPeak` | Copied from the pooled-pseudoreplicate thresholded IDR output |
| Reproducibility summary | `06_idr/final/reproducibility_summary.tsv` | Nt (total peaks), Np (passing peaks), N1/N2 (self-consistent peaks per replicate), rescue ratio, self-consistency ratio, and pass/fail status for each IDR stage |

## Configuration

```yaml
stage5: true
stage4b: true           # required
idr:
  seed: 42              # pseudorep split seed (positive int)
  threshold: 0.05       # passed directly to idr --idr-threshold
  rank: "p.value"       # p.value or signal.value

tool_parameters:
  idr_macs3:
    pvalue: 0.1         # relaxed p-value for IDR-ready peak calls
```

## Manifest coverage

The Stage 25 result manifest records only the final IDR outputs:
- `idr_conservative` (`conservative.narrowPeak`)
- `idr_optimal` (`optimal.narrowPeak`)
- `idr_reproducibility_summary` (`reproducibility_summary.tsv`)

Raw and intermediate IDR outputs (true-replicate raw/thresholded, self-IDR,
pooled-IDR, pseudorep peaks) are documented in
[`docs/output-contract.md`](output-contract.md) but not included in the
minimal manifest.

## Out of scope (deferred beyond v0.2)

- **3+ replicate selection:** Automatic pairwise selection among >=3 replicates is not implemented.
- **Broad/histone IDR:** Histone broad-mark experiments do not produce IDR peak sets.
- **CUT&Tag IDR:** CUT&Tag samples are excluded from IDR.
- **ATAC IDR:** ATAC-seq samples are excluded from IDR.
- **Multi-assay IDR:** No cross-assay IDR.

## Relationship to assay policy

See [`docs/assay-policy.md`](assay-policy.md) for the full assay-specific
duplicate, MACS3, and read-extension policies that affect IDR input quality.
