# Assay Policy

This document defines the v0.2 behavioral contract for each supported assay.
It describes actual implemented behavior — no aspirational features.

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
| `auto` | No extension (fragments from pairs) | MACS3 predicted fragment size, fallback 200 bp |
| `yes` | `--extendReads` (pair inference) | Same as `auto` |
| `no` | No extension | No extension |
| `<int>` | `--extendReads <int>` | `--extendReads <int>` |

### Controls

Supported via `use_control: true` with either `control_sample` (FASTQ row) or `control_bam` (external BAM path). Control samples go through full preprocessing but skip peak calling.

### IDR eligibility

Only for ChIP-seq, `peak_mode: narrow`, exactly 2 treatment biological replicates per experiment. See [`docs/idr-contract.md`](idr-contract.md).

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

### IDR

Not supported in v0.2. CUT&Tag IDR is deferred to a future release.

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

### Not implemented in v0.2

- ATAC-specific footprinting
- Nucleosome positioning
- Broad peak mode
- ATAC IDR

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
