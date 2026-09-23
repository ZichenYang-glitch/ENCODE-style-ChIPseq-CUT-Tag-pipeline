# Reproducibility Policy: Replicate-Validated Peak Outputs

**Status:** Current implemented policy

`replicate_analysis` controls replicate pooling; `chipseq_idr` controls
ChIP-seq narrow IDR. Both are public configuration keys.

## 1. Purpose

This document defines the reproducibility policy for peak outputs in the
ChIP-seq / CUT&Tag / ATAC-seq pipeline. It establishes:

- What counts as a replicate-validated peak set
- What does not
- Which reproducibility methods apply to which assay/caller/mode combinations
- Output path conventions
- Configuration surface
- How reproducibility interacts with the existing `chipseq_idr` ChIP-seq narrow IDR

## 2. Key Concepts

### 2.1 Pooled Peaks Are Not Validated Peaks

Pooled peaks (`results/experiments/<exp>/04_peaks/pooled/`) are MACS3 runs on
merged BAMs pooled across biological replicates. They aggregate signal for
statistical power but do **not** constitute replication-based validation.

**Rule:** Pooled peaks MUST NOT be labeled, classified, or treated as
replicate-validated peak sets. They are aggregate-signal peaks only.

### 2.2 Consensus

Consensus evaluates per-replicate peak calls and retains regions supported by
at least `min_replicates` distinct biological replicates.

- **Graph and components:** [`compute_consensus.py`](../scripts/compute_consensus.py)
  adds an edge when the overlap covers at least `reciprocal_overlap` (default
  0.5) of **both** source peaks. It then finds connected components.
- **N-of-M support:** Count distinct biological replicate labels in each
  component and retain it when `support_count >= min_replicates` (default 2).
  Several peaks from one biorep count as one supporting replicate.
- **Output interval:** Merge the component from its minimum start to maximum
  end. With two bioreps this is not a literal interval intersection. Transitive
  connections can join source peaks that do not directly overlap, and not
  every base of the merged interval must have N-replicate support.
- This is the existing algorithm, not a promise of equivalence to a particular
  ENCODE overlap implementation. Consensus is the baseline method for eligible
  peak modes with at least two treatment bioreps **when enabled**. Default
  targets also require `reproducibility.enabled`, its consensus child switch
  and `replicate_analysis`; the parent defaults to false.

### 2.3 IDR (Irreproducible Discovery Rate)

IDR compares ranked peak calls between two replicates to estimate an
irreproducible discovery rate.

- **Scope:** Exactly 2 biological replicates. 3+ replicate IDR is not in
  scope for this policy version.
- **Rank metric:** Configurable (default `p.value` from MACS3 narrowPeak
  column 8).
- **Production-supported IDR:** ChIP-seq narrow (`chipseq_idr`) and
  ATAC-seq narrow (`reproducibility.idr.atac_narrow`). IDR is final when
  enabled; consensus is secondary/report.
- **Supported opt-in IDR:** CUT&Tag narrow (`reproducibility.idr.cuttag_narrow`).
  Not default. IDR becomes final when explicitly enabled.
- **Experimental opt-in IDR:** ChIP-seq broad, CUT&Tag broad
  (`reproducibility.idr.chipseq_broad_experimental`,
  `cuttag_broad_experimental`). Not default. IDR becomes final when explicitly
  enabled and eligible; consensus remains available as fallback/report.
  Broad experimental IDR becomes final when explicitly enabled and eligible.
- **Not planned for IDR:** SEACR (BED format score scheme not directly
  compatible with IDR rank assumptions) and MNase.

### 2.4 ChIP-seq narrow IDR (`chipseq_idr`)

The existing `chipseq_idr: true` configuration enables ChIP-seq narrow IDR with
exactly 2 biological replicates. This behavior and its output paths are
preserved exactly as-is. No new configuration may disable, rename, or move
the `chipseq_idr` output paths.

The `reproducibility` block is an **orthogonal** layer. `chipseq_idr` controls the
ChIP-seq narrow IDR path; `reproducibility` controls consensus and expanded IDR modes.

## 3. Reproducibility Strategy Matrix

| # | Assay | Peak mode | Caller | Primary reproducibility | Secondary / report | Notes |
|---|-------|-----------|--------|------------------------|-------------------|-------|
| 1 | chipseq | narrow | MACS3 | IDR (`chipseq_idr`) | Consensus | `chipseq_idr` behavior unchanged |
| 2 | chipseq | broad | MACS3 | IDR when experimental flag enabled; consensus otherwise | Consensus fallback/report | Experimental opt-in only |
| 3 | cuttag | narrow | MACS3 | Consensus | IDR opt-in (supported) | IDR final when explicitly enabled; consensus otherwise |
| 4 | cuttag | broad | MACS3 | IDR when experimental flag enabled; consensus otherwise | Consensus fallback/report | Experimental opt-in only |
| 5 | cuttag | — | SEACR | Consensus | No IDR planned | SEACR IDR is outside the current policy |
| 6 | atac | narrow | MACS3 | IDR (when enabled) | Consensus | ATAC narrow IDR uses same narrowPeak IDR machinery |

## 4. Output Paths

### 4.1 Reproducibility namespace

```
results/experiments/<exp>/06_reproducibility/
  consensus/
    <exp>.<assay>.<caller>.<peak_mode>.consensus.<suffix>
    <exp>.<assay>.<caller>.<peak_mode>.consensus.summary.tsv
  idr/
    idr_peaks/...
    true_replicates/...
    self_pseudoreplicates/...
    pooled_pseudoreplicates/...
  final/
    <exp>.<assay>.<caller>.<peak_mode>.replicate_validated.<method>.<suffix>
    reproducibility_summary.tsv
```

### 4.2 ChIP-seq narrow IDR paths

```
results/experiments/<exp>/06_idr/
  true_replicates/
    idr.txt
    idr.thresholded.narrowPeak
  self_pseudoreplicates/
    biorep1.idr.txt
    biorep1.idr.thresholded.narrowPeak
    biorep2.idr.txt
    biorep2.idr.thresholded.narrowPeak
  pooled_pseudoreplicates/
    idr.txt
    idr.thresholded.narrowPeak
  final/
    conservative.narrowPeak
    optimal.narrowPeak
    reproducibility_summary.tsv
```

### 4.3 Consensus Filename Templates

| Assay | Caller | Peak mode | Consensus peak file | Summary file |
|-------|--------|-----------|--------------------|-------------|
| chipseq | macs3 | narrow | `<exp>.chipseq.macs3.narrow.consensus.narrowPeak` | `<exp>.chipseq.macs3.narrow.consensus.summary.tsv` |
| chipseq | macs3 | broad | `<exp>.chipseq.macs3.broad.consensus.broadPeak` | `<exp>.chipseq.macs3.broad.consensus.summary.tsv` |
| cuttag | macs3 | narrow | `<exp>.cuttag.macs3.narrow.consensus.narrowPeak` | `<exp>.cuttag.macs3.narrow.consensus.summary.tsv` |
| cuttag | macs3 | broad | `<exp>.cuttag.macs3.broad.consensus.broadPeak` | `<exp>.cuttag.macs3.broad.consensus.summary.tsv` |
| cuttag | seacr | `<mode>` | `<exp>.cuttag.seacr.<mode>.consensus.bed` | `<exp>.cuttag.seacr.<mode>.consensus.summary.tsv` |
| atac | macs3 | narrow | `<exp>.atac.macs3.narrow.consensus.narrowPeak` | `<exp>.atac.macs3.narrow.consensus.summary.tsv` |

### 4.4 Final Validated Peak Filename Templates

| Assay | Caller | Peak mode | Primary method | Final output |
|-------|--------|-----------|---------------|-------------|
| chipseq | macs3 | narrow | ChIP-seq narrow IDR | Use `06_idr/final/conservative.narrowPeak` |
| chipseq | macs3 | broad | IDR (when `reproducibility.idr.chipseq_broad_experimental: true`) | `<exp>.chipseq.macs3.broad.replicate_validated.idr.broadPeak` |
| chipseq | macs3 | broad | Consensus (when broad IDR not enabled) | `<exp>.chipseq.macs3.broad.replicate_validated.consensus.broadPeak` |
| cuttag | macs3 | narrow | IDR (when `reproducibility.idr.cuttag_narrow: true`) | `<exp>.cuttag.macs3.narrow.replicate_validated.idr.narrowPeak` |
| cuttag | macs3 | narrow | Consensus (when CUT&Tag IDR not enabled) | `<exp>.cuttag.macs3.narrow.replicate_validated.consensus.narrowPeak` |
| cuttag | macs3 | broad | IDR (when `reproducibility.idr.cuttag_broad_experimental: true`) | `<exp>.cuttag.macs3.broad.replicate_validated.idr.broadPeak` |
| cuttag | macs3 | broad | Consensus (when broad IDR not enabled) | `<exp>.cuttag.macs3.broad.replicate_validated.consensus.broadPeak` |
| cuttag | seacr | `<mode>` | Consensus | `<exp>.cuttag.seacr.<mode>.replicate_validated.consensus.bed` |
| atac | macs3 | narrow | IDR (when `reproducibility.idr.atac_narrow: true`) | `<exp>.atac.macs3.narrow.replicate_validated.idr.narrowPeak` |
| atac | macs3 | narrow | Consensus (when ATAC IDR not enabled) | Consensus under `consensus/`; no `final/` output until IDR enabled |

For chipseq narrow, `06_reproducibility/reproducibility_summary.tsv` records
the primary final source as `../06_idr/final/conservative.narrowPeak`.
No symlinks or copies are created into `06_reproducibility/final/`.

### 4.5 Consensus Summary TSV Schema

Each consensus output is accompanied by a summary TSV:

```
experiment   assay   peak_mode   caller   n_bioreps   min_replicates
reciprocal_overlap   consensus_peak_count   support_distribution
biorep_labels   source_peak_files   final_method   final_output
```

| Column | Type | Description |
|--------|------|-------------|
| `experiment` | string | Experiment identifier |
| `assay` | string | `chipseq`, `cuttag`, or `atac` |
| `peak_mode` | string | `narrow`, `broad`, or SEACR mode name |
| `caller` | string | `macs3` or `seacr` |
| `n_bioreps` | int | Number of biological replicates participating |
| `min_replicates` | int | Minimum supporting replicates threshold |
| `reciprocal_overlap` | float | Reciprocal overlap fraction used |
| `consensus_peak_count` | int | Number of consensus peaks retained |
| `support_distribution` | JSON | Map of support count → peak count, e.g. `{"2": 15000, "3": 8000}` |
| `biorep_labels` | string | Comma-separated biorep labels |
| `source_peak_files` | JSON array | Paths to per-biorep peak files used as input |
| `final_method` | string | `consensus`, `idr`, or `none`, as selected by the mode helper |
| `final_output` | string | Path reported by the mode helper; may be empty for a report-only consensus |

### 4.6 Consensus Output Format Caveats

Consensus narrowPeak and broadPeak files use the standard format suffixes for
compatibility, but their field semantics differ from single-sample MACS3
outputs:

- **chrom / start / end / name:** Represent the merged consensus interval.
  `name` is derived deterministically (e.g., consensus peak index).
- **score:** `int(1000 * support_count / n_bioreps + 0.5)`, where support is
  the number of distinct bioreps in the component. Three of three gives 1000;
  two of three gives 667. This is neither constant nor the maximum signal.
- **signalValue:** Maximum signal among supporting source peaks; it is not a
  newly fitted MACS3 or IDR statistic.
- **pValue / qValue:** Both are `-1` (not assigned), as permitted by the
  [UCSC narrowPeak/broadPeak specification](https://genome.ucsc.edu/FAQ/FAQformat.html#format12).
  The sentinel itself is not a format defect.
- **Summary scope:** The table in §4.5 is one aggregate data row with 13
  columns. `support_distribution` counts retained components by support level;
  `biorep_labels` and `source_peak_files` describe the input set. There are no
  per-peak `support_count`, `supporting_bioreps` or `source_peak_count` columns.
  `final_method`/`final_output` describe the selected final method/path; the
  consensus script does not create a second product at `final_output`.

Implementation: `_build_consensus_row` and `write_summary_tsv` in [`compute_consensus.py`](../scripts/compute_consensus.py).
Existing checks: `test_support_fraction_sets_consensus_score`,
`test_overlap_components_respect_interval_topology`, and
`test_summary_schema_and_values` in
[`test_compute_consensus.py`](../test/scripts/test_compute_consensus.py).

## 5. Config Reference

### 5.1 `reproducibility` block

```yaml
reproducibility:
  enabled: false

  consensus:
    enabled: true
    min_replicates: 2
    reciprocal_overlap: 0.5

  idr:
    chipseq_narrow: null
    atac_narrow: false
    cuttag_narrow: false
    chipseq_broad_experimental: false
    cuttag_broad_experimental: false
```

### 5.2 Interaction with `chipseq_idr`

| `chipseq_idr` | `reproducibility.enabled` | Result |
|----------|--------------------------|--------|
| false | false | No IDR, no reproducibility outputs |
| true | false | ChIP-seq narrow IDR only |
| true | true | ChIP-seq narrow IDR + expanded reproducibility outputs |
| false | true | Expanded reproducibility only (no ChIP-seq narrow IDR) |

**Invariant:** `chipseq_idr: true` always keeps existing ChIP-seq narrow IDR
behavior. No `reproducibility` setting may break, disable, rename, or move
existing `chipseq_idr` outputs.

### 5.3 Validation Rules

When `reproducibility.enabled: true`:

| Field | Constraint | Default |
|-------|-----------|---------|
| `consensus.min_replicates` | int ≥ 2 | 2 |
| `consensus.reciprocal_overlap` | float in (0, 1] | 0.5 |
| `idr.chipseq_narrow` | true / false / null | null (inferred from `chipseq_idr`) |
| `idr.atac_narrow` | bool | false |
| `idr.cuttag_narrow` | bool | false |
| `idr.chipseq_broad_experimental` | bool | false |
| `idr.cuttag_broad_experimental` | bool | false |

- If `idr.chipseq_narrow: false` while `chipseq_idr: true`: emit config
  contradiction warning, but `chipseq_idr` still runs.
- If `*_experimental: true`: the validator emits an informational experimental
  warning. Its claim that consensus remains primary is stale; the current DAG
  and manifest select the eligible IDR artifact as final, as described above.
  The wording mismatch is tracked in [Known Issues](../KNOWN_ISSUES.md#scientific-scope).
- If `reproducibility.enabled: false`: all sub-keys are ignored for DAG
  purposes with no validation errors.

## 6. Current non-goals

1. No MNase reproducibility. Nucleosome calling has different semantics.
2. No breaking, renaming, or moving existing output paths.
3. No default broad IDR. Broad IDR remains experimental and requires explicit
   opt-in. SEACR IDR is outside the current policy.
4. No labeling pooled peaks as replicate-validated.
5. No 3+ biorep IDR. Multi-replicate IDR requires separate pairwise-IDR or
   replicate-selection policy design.
6. No replacing `chipseq_idr` with `reproducibility`.

## 7. Implementation status

| Capability | Current state |
| --- | --- |
| ChIP-seq narrow IDR | Implemented through the `chipseq_idr` path. |
| Replicate consensus | Implemented for eligible ChIP-seq, CUT&Tag, ATAC-seq, and SEACR outputs. |
| ATAC-seq narrow IDR | Implemented when explicitly enabled. |
| CUT&Tag narrow IDR | Implemented as an opt-in policy. |
| ChIP-seq and CUT&Tag broad IDR | Implemented as experimental opt-ins; consensus remains available. |
| Manifest and reporting | Integrated with the maintained output contract and QC summaries. |
