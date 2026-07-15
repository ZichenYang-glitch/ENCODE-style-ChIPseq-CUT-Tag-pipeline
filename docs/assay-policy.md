# Assay Policy

This document defines the current behavioral contract for each supported
assay. It describes implemented behavior, not aspirational features.

## ChIP-seq

**Allowed peak modes:** `narrow`, `broad`

### MACS3 parameters

| Peak mode | Arguments |
| :--- | :--- |
| `narrow` | `-f <fmt> -g <genome> -q <qvalue>` (model-based) |
| `broad` | `-f <fmt> -g <genome> -q <qvalue> --broad --broad-cutoff <cutoff>` |

- `--format` (`-f`): `BAMPE` for PE layout, `BAM` for SE.
- `--gsize` (`-g`): resolved from `genome_resources.<genome>.effective_genome_size` (MACS3 shortcut `hs`/`mm` or explicit integer).
- `--qvalue` (`-q`): from `tool_parameters.macs3.qvalue`, default `0.01`.
- `--broad-cutoff`: from `tool_parameters.macs3.broad_cutoff`, default `0.1`.
- Control (`-c`): resolved from `control_sample` or `control_bam` when `use_control: true`.

### Duplicate policy

| `remove_dup` config | Behavior |
| :--- | :--- |
| `auto` | `yes` for narrow, `no` for broad |
| `yes` | Always remove duplicates |
| `no` | Never remove duplicates |

### Read extension (bamCoverage)

| `extend_reads` | PE | SE |
| :--- | :--- | :--- |
| `auto` | No extension (fragments from pairs) | MACS3 predicted fragment size, fallback 200 bp (warning emitted to stderr when fallback used) |
| `yes` | `--extendReads` (pair inference) | Same as `auto` |
| `no` | No extension | No extension |
| `<int>` | `--extendReads <int>` | `--extendReads <int>` |

### Controls

Supported via `use_control: true` with either `control_sample` (FASTQ row) or `control_bam` (external BAM path). Control samples go through full preprocessing but skip peak calling.

### Reproducibility

The legacy `stage5` path supports ChIP-seq narrow IDR with exactly two
treatment biological replicates. The separate `reproducibility` block provides
replicate consensus and experimental opt-in broad IDR under the maintained
[reproducibility policy](reproducibility-policy.md). See the
[legacy IDR contract](idr-contract.md) for `stage5` details.

---

## CUT&Tag

**Allowed peak modes:** `narrow`, `broad`

### MACS3 parameters

| Peak mode | Arguments |
| :--- | :--- |
| `narrow` | `-f <fmt> -g <genome> -q <qvalue> --nomodel --shift -100 --extsize 200` |
| `broad` | `-f <fmt> -g <genome> -q <qvalue> --broad --broad-cutoff <cutoff>` |

The narrow-peak mode uses Tn5-aware parameters (`--nomodel --shift -100 --extsize 200`). Broad mode does not apply the Tn5 shift (per ENCODE CUT&Tag guidelines).

### Duplicate policy

Identical to ChIP-seq (`auto` → narrow-yes/broad-no). The policy delegate (`get_remove_dup_cuttag`) calls `get_remove_dup_chipseq` directly.

### Read extension (bamCoverage)

Identical to ChIP-seq. The policy delegate (`get_extend_reads_cuttag`) calls `get_extend_reads_chipseq` directly.

### SEACR sidecar

Gated by `cuttag.seacr.enabled: true`. Applies only to CUT&Tag PE treatment samples.
Produces bedGraph and SEACR peak BED outputs under `results/<sample>/04_peaks_seacr/`.
MACS3 peaks are still produced as the canonical peak set.

### Fragment-size QC

Enabled by default (`qc.cuttag_fragment_size: true`). Applies only to samples with `assay: cuttag`. Produces `results/<sample>/01_qc/<sample>.cuttag_fragment_size.tsv`.

### Reproducibility

Replicate consensus is available when `reproducibility` is enabled. MACS3
narrow IDR is a supported opt-in for exactly two biological replicates; broad
IDR is experimental and opt-in. SEACR remains consensus-only. See the
[reproducibility policy](reproducibility-policy.md) for eligibility and final
output selection.

---

## ATAC-seq

**Allowed peak modes:** `narrow` only (`broad` raises a validation error).

### MACS3 parameters

`-f <fmt> -g <genome> -q <qvalue> --nomodel --shift -100 --extsize 200`

Tn5-aware shift/extension matching common ATAC-seq practice.

### Duplicate policy

| `remove_dup` config | Behavior |
| :--- | :--- |
| `auto` | `yes` |
| `yes` | Always remove duplicates |
| `no` | Never remove duplicates |

ATAC-seq always removes duplicates in `auto` mode (no broad-peak exception).

### Read extension (bamCoverage)

| `extend_reads` | PE | SE |
| :--- | :--- | :--- |
| `auto` | No extension (fragments from pairs) | `--extendReads 200` |
| `yes` | `--extendReads` (pair inference) | `--extendReads 200` |
| `no` | No extension | No extension |
| `<int>` | `--extendReads <int>` | `--extendReads <int>` |

### Reproducibility and current limits

Replicate consensus and opt-in narrow IDR are implemented under the
[`reproducibility` policy](reproducibility-policy.md).

- ATAC-specific footprinting
- Nucleosome positioning
- Broad peak mode

---

## Cross-Cutting Behaviors

### Controls

- `use_control: true` enables control resolution. When `false`, `control_sample` and `control_bam` are silently ignored.
- Control rows (`role: control`) are preprocessed but skip peak calling.
- External `control_bam` is not processed by the pipeline — used directly in MACS3 `-c`.

### Replicate model

- `stage4b: true` (default) enables replicate-aware outputs.
- Technical replicates (`technical_replicate`) are merged into biological-replicate BAMs.
- Pooled outputs only for experiments with >=2 unique `biological_replicate` values.
- Pooled controls are produced when a multi-biorep treatment experiment has controls referenced.

### Signal tracks

- `qc.signal_tracks: true` produces MACS3 FE and ppois bedGraph files per treatment sample.
- Pooled FE/ppois bedGraph is produced for multi-biorep experiments.
- FE/ppois BigWig conversion is gated on `genome_resources.<genome>.chrom_sizes` being non-empty.
- CPM BigWig (`bamCoverage`) is always produced for all active samples.

### Summary and manifest

- Per-sample `qc_summary.tsv` (37 columns) is assembled by `scripts/assemble_qc_summary.py`.
- Project-level `stage3_qc_summary.tsv` is aggregated by `scripts/aggregate_qc_summary.py`.
- `result_manifest.tsv` records core output existence with 10-column TSV, using `validate_samples` for DAG-consistent gating.
- MNase samples do not produce `qc_summary.tsv` (no peaks); instead they produce `mnase_qc_summary.tsv`.

---

## MNase-seq

**Allowed peak modes:** `nucleosome` only

### PE-only

MNase-seq requires paired-end layout. SE MNase samples raise a validation error.

### No MACS3 peak calling

MNase samples skip the MACS3 peak calling path entirely. No `04_peaks/`,
MACS3 FE/ppois bedGraph, or MACS3 FE/ppois BigWig outputs are produced.
MNase does produce dyad (`bamCoverage --MNase`) and mono occupancy BigWigs
under `04_signal/`.

### Duplicate policy

| `remove_dup` config | Behavior |
| :--- | :--- |
| `auto` | `yes` |
| `yes` | Always remove duplicates |
| `no` | Never remove duplicates |

MNase is PE-only and always removes duplicates in `auto` mode (no broad-peak exception).

### Read extension (bamCoverage)

PE mode: no read extension. Fragment pairs are the biological unit for MNase-seq.

### Nucleosome-centric outputs

| Output | Rule | Tool | Input |
| :--- | :--- | :--- | :--- |
| `03_fragments/<s>.sub.bam` | `mnase_split_sub` | `alignmentSieve` | `final.bam` |
| `03_fragments/<s>.mono.bam` | `mnase_split_mono` | `alignmentSieve` | `final.bam` |
| `03_fragments/<s>.di.bam` | `mnase_split_di` | `alignmentSieve` | `final.bam` |
| `04_signal/<s>.dyad.CPM.bw` | `mnase_dyad_bigwig` | `bamCoverage --MNase --binSize 1` | `final.bam` |
| `04_signal/<s>.mono.CPM.bw` | `mnase_mono_bigwig` | `bamCoverage` | `mono.bam` |

Fragment ranges default to sub [1,139], mono [140,200], di [300,400].
Configurable via `mnase.fragments` in `config/config.yaml`. The legacy
`mnase.mono_range` key is still accepted but `fragments.mono` takes precedence.

Explicit dyad range for `bamCoverage --MNase` is controlled by
`mnase.dyad_range` (default [130, 200]). This may differ from
`fragments.mono` — e.g. a wider dyad capture window is common.

`tool_parameters.bamcoverage.extra_args` must NOT set `--MNase`,
`--minFragmentLength`, or `--maxFragmentLength` when running MNase dyad
BigWig rules; these flags are workflow-managed.

Per-sample QC summary (`01_qc/<s>.mnase_qc_summary.tsv`) records fragment
stratification metadata, read counts, and caller configuration status.

### Pooled outputs

For MNase experiments with >=2 biological replicates, pooled outputs are produced:
`<e>.pooled.mono.bam`, `<e>.pooled.dyad.CPM.bw`, `<e>.pooled.mono.CPM.bw`.

Pooled BAMs reuse the existing `stage4b` replicate merge logic
(assay-agnostic).

### Controls

Supported for preprocessing (control samples produce `final.bam` and CPM BigWig).
MNase downstream rules do not use controls for differential analysis.

### IDR

Not supported.

### Config

```yaml
mnase:
  fragments:
    sub: [1, 139]
    mono: [140, 200]
    di: [300, 400]
  dyad_range: [130, 200]
  callers:
    danpos3: false
    inps: false
    sem: false
```

The `mnase` block is optional. If absent, all defaults apply.
`fragments.mono` takes precedence over top-level `mono_range`.
`dyad_range` defaults to `[130, 200]` (deepTools `--MNase` default).
`callers` surface is reserved; setting any caller to `true` raises a
validation error because caller execution is not implemented.

### Not implemented

- DANPOS3 nucleosome calling
- iNPS high-resolution nucleosome detection
- SEM nucleosome subtype analysis
- Rotational periodicity QC
- NFR/TSS MNase-specific QC profiles
- Pooled sub-nucleosome and di-nucleosome BAM outputs
