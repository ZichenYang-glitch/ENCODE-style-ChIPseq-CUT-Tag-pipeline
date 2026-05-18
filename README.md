# ENCODE-style ChIP-seq and CUT&Tag Pipeline

[![Snakemake](https://img.shields.io/badge/Snakemake-%3E%3D8.0-brightgreen.svg?style=flat-square)](https://snakemake.github.io)
[![Conda](https://img.shields.io/badge/conda-supported-blue.svg?style=flat-square)](https://docs.conda.io/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square)](https://opensource.org/licenses/MIT)

## Overview

A Snakemake-based pipeline suite for ChIP-seq and CUT&Tag data analysis. It
handles single-sample preprocessing as well as multi-replicate experiments with
pooled outputs, single-sample QC, and TF ChIP-seq IDR reproducibility analysis.

Default mode is no-input / no-control. Optional controls can be enabled with
`use_control: true` and supplied as an external control BAM or a FASTQ-based
control sample row. Dependencies are managed with Conda.

## Key Features

- **Shared preprocessing:** FastQC, Trim Galore, Bowtie2 alignment, MAPQ
  filtering, Picard duplicate handling, flagstat, idxstats, BigWig generation
- **ChIP-seq / CUT&Tag assay policies:** assay-aware MACS3 parameters, duplicate
  removal, and read extension; optional CUT&Tag SEACR sidecar peak calls
  (`cuttag.seacr.enabled`, output under `results/<sample>/04_peaks_seacr/`)
- **Optional controls:** external control BAM or FASTQ-based control rows;
  control samples processed through the same pipeline
- **Single-sample QC:** FRiP, peak counts, library complexity, NRF/PBC, MACS3
  FE/ppois signal tracks, CUT&Tag fragment-size QC, per-sample QC summaries,
  and project-level aggregation
- **Replicate-aware outputs:** technical replicate merging into biological-replicate
  BAMs, pooled treatment/control BAMs, pooled MACS3 peak calls, and pooled signal
  tracks for multi-replicate experiments
- **TF ChIP-seq IDR:** true-replicate IDR, pseudoreplicate-based self-IDR and
  pooled-IDR, reproducibility metrics, and final conservative/optimal peak sets
  (narrowPeak, exactly 2 biological replicates)
- **Histone-aware QC:** experiment-level pooled QC summaries with target
  classification and peak-mode compatibility status
- **Configurable tool parameters:** optional `tool_parameters` blocks for common
  knobs (MACS3 q-value, bamCoverage normalization, etc.) with `extra_args`
  escape hatches

## Quick Start

### 1. Install

```bash
git clone https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline.git
cd ENCODE-style-ChIPseq-CUT-Tag-pipeline
conda env create -f workflow/envs/chipseq.yml
conda activate chipseq
```

### 2. Configure samples

Edit `config/samples.tsv` (see [Sample Sheet](#sample-sheet) below).

### 3. Adjust workflow options

Edit `config/config.yaml`. Defaults are reasonable for most ChIP-seq runs:

```yaml
samples: "config/samples.tsv"
outdir: "results"
threads: 8
mapq: 30
binsize: 10
remove_dup: "auto"
trim: true
extend_reads: "auto"
use_control: false
multiqc: true
```

### 4. Validate

```bash
python3 scripts/validate_samples.py --config config/config.yaml
snakemake -s workflow/Snakefile --configfile config/config.yaml -n
```

### 5. Run

```bash
# Dry-run
snakemake -s workflow/Snakefile --configfile config/config.yaml -n --use-conda

# Run with 16 cores
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 16 --use-conda

# Resume after interruption
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 16 --use-conda --rerun-incomplete
```

## Sample Sheet

### Required columns

| Column | Description |
| :--- | :--- |
| `sample` | Unique sample ID. Allowed: `A-Z a-z 0-9 _ . -`. |
| `fastq_1` | R1 FASTQ path. |
| `fastq_2` | R2 FASTQ path. Required for PE; leave empty for SE. |
| `layout` | `PE` or `SE`. |
| `assay` | `chipseq` or `cuttag`. |
| `target` | Antibody or target name (e.g. `H3K27ac`, `CTCF`). |
| `peak_mode` | `narrow` or `broad`. |
| `genome` | Genome label used to look up `genome_resources` (e.g. `hs`, `mm`, `hg38`). |
| `bowtie2_index` | Bowtie2 index prefix path. |

### Optional columns

**Replicate metadata:**

| Column | Default | Description |
| :--- | :--- | :--- |
| `experiment` | `<sample>` | Experiment group identifier. |
| `condition` | `<target>` | Condition label (treatment, genotype, timepoint). |
| `replicate` | `1` | Replicate number within an experiment. |
| `biological_replicate` | `<replicate>` | Biological replicate number. |
| `technical_replicate` | `1` | Technical replicate number within a biological replicate. |

**Controls:**

| Column | Description |
| :--- | :--- |
| `role` | `treatment` (default) or `control`. Peak calling runs only for treatment. |
| `control_sample` | Sample ID of a control row to use as MACS3 input. |
| `control_bam` | Path to an external control BAM. |

Technical replicates within a biological replicate are merged into a
`biorep<N>.final.bam`. Experiments with >=2 biological replicates produce
pooled BAMs and pooled peaks.

A minimal sample sheet with only the required columns is fully supported
(single-sample, no controls, no replicates).

### Example

```tsv
sample	fastq_1	fastq_2	layout	assay	target	peak_mode	genome	bowtie2_index	experiment	biological_replicate
H3K27AC_rep1	/data/ac1_R1.fq.gz	/data/ac1_R2.fq.gz	PE	chipseq	H3K27ac	narrow	hs	/path/to/bt2/GRCh38	H3K27AC	1
H3K27AC_rep2	/data/ac2_R1.fq.gz	/data/ac2_R2.fq.gz	PE	chipseq	H3K27ac	narrow	hs	/path/to/bt2/GRCh38	H3K27AC	2
```

## Configuration

Core config keys with their defaults are listed in [Quick Start](#quick-start)
above. Additional feature-specific blocks are described here.

### Genome resources

```yaml
genome_resources:
  hg38:
    effective_genome_size: 2913022398
    chrom_sizes: "/data/genomes/hg38/hg38.chrom.sizes"  # optional
    blacklist: "/data/genomes/hg38/hg38.blacklist.bed"    # optional
    gtf: ""   # optional
    reference_fasta: ""  # optional
```

`effective_genome_size`: MACS3 shortcut (`hs`, `mm`) or positive integer.
Other fields are optional paths. If non-empty they must exist on disk.

### QC switches

```yaml
qc:
  blacklist_filter: true
  frip: true
  library_complexity: true
  nrf_pbc: true
  signal_tracks: true    # MACS3 FE/ppois bedGraph tracks (per-sample + pooled)
  summary: true          # per-sample QC summary + project-level aggregate
```

### Replicate and IDR features

```yaml
stage4b: true            # replicate-aware pooled BAMs and pooled peaks (default on)
stage5: false            # TF ChIP-seq IDR (default off; requires stage4b: true)

idr:
  seed: 42
  threshold: 0.05        # IDR threshold for thresholded peak set
  rank: "p.value"        # ranking measure: p.value or signal.value
```

### Tool parameters

Optional `tool_parameters` blocks let you tune individual tools without editing
rule files. Absent keys use the built-in defaults. Example:

```yaml
tool_parameters:
  macs3:
    qvalue: 0.01
  idr_macs3:
    pvalue: 0.1          # relaxed p-value for IDR-ready peak calls
  bamcoverage:
    normalize_using: "CPM"
```

See `workflow/schemas/config.schema.yaml` for the full set of configurable keys
and their accepted types.

## Workflow Capabilities

### Preprocessing and peak calling

Every sample (treatment or control) flows through: FastQC → Trim Galore
(or symlink when `trim: false`) → Bowtie2 alignment → samtools sort/index →
MAPQ filter → duplicate handling (Picard or samtools fallback) → `final.bam`.
BigWig tracks are generated with deepTools `bamCoverage` (CPM-normalized by
default). MACS3 peak calling runs on treatment samples, with assay-specific
parameters (TF ChIP-seq model-based, CUT&Tag Tn5-aware `--shift -100`).

### Controls

Controls are disabled by default (`use_control: false`). When enabled,
each treatment row may reference a `control_sample` (another sample row with
`role=control`) or an external `control_bam` path. Control rows are processed
through the full preprocessing pipeline and their `final.bam` is passed as
MACS3 `-c`.

### Single-sample QC

When the `qc` block is enabled, each treatment sample receives:
- **Blacklist filtering** (BAM + peaks) when a blacklist BED is configured
- **FRiP** (Fraction of Reads in Peaks)
- **Library complexity** (Picard duplication-derived metrics)
- **NRF/PBC** (BAM-derived library complexity)
- **MACS3 signal tracks**: fold-enrichment (`FE.bdg`) and p-value
  (`ppois.bdg`) bedGraph from `macs3 bdgcmp`
- Per-sample QC summary TSV and a project-level aggregate at
  `results/multiqc/stage3_qc_summary.tsv`

### Replicates and pooled outputs

When `stage4b: true` and the sample sheet defines experiments with >=2
biological replicates:

- Technical replicates are merged into `biorep<N>.final.bam` files (symlink
  when there is only one technical replicate)
- Pooled treatment BAMs are produced by merging all treatment biorep BAMs
- Pooled control BAMs are produced when controls are referenced (symlink for
  1 control, merge for >1)
- Pooled MACS3 peak calling runs on the pooled treatment BAM
- Pooled MACS3 FE/ppois signal tracks are produced when
  `qc.signal_tracks: true`

Outputs land under `results/experiments/<experiment>/`.

### TF ChIP-seq IDR

Enable with `stage5: true` (requires `stage4b: true`). Supported:
chipseq + narrowPeak + exactly 2 treatment biological replicates.

- **True-replicate IDR:** per-biorep IDR-ready MACS3 calls with relaxed `-p`
  threshold (`tool_parameters.idr_macs3`), then `idr --samples` between the
  two biorep peak sets. Produces raw `idr.txt` and
  `idr.thresholded.narrowPeak`.
- **Pseudoreplicate IDR:** deterministic pseudoreplicate BAM splitting,
  self-IDR per biorep, pooled-IDR, and final `conservative.narrowPeak` /
  `optimal.narrowPeak` peak sets with a `reproducibility_summary.tsv`.

IDR output paths are under `results/experiments/<exp>/06_idr/`.

### Histone pooled QC

For any multi-biorep experiment (independent of IDR), the pipeline produces a
pooled experiment QC summary at
`results/experiments/<exp>/01_qc/<exp>.pooled_qc_summary.tsv`. The summary
classifies the target into broad-like, narrow-like, context-dependent, or
unknown, and reports peak-mode compatibility status (`ok` / `mismatch` /
`unknown`). No hard-fail validation on mismatches.

Context-dependent marks (H3K27ac, H3K4me1, H3K4me2) accept both broad and
narrow peak_mode without warning.

### Reports and MultiQC

When `multiqc: true`, MultiQC aggregates QC artifacts from all active samples
into `results/multiqc/multiqc_report.html`.

## Output Structure

```text
results/
├── <sample>/
│   ├── 00_raw/                          # trimmed (or symlinked) FASTQs
│   ├── 01_qc/                           # QC metrics and summaries
│   │   └── <sample>.qc_summary.tsv
│   ├── 02_align/
│   │   ├── <sample>.sorted.bam(.bai)
│   │   ├── <sample>.final.bam(.bai)     # stable downstream BAM contract
│   │   └── <sample>.blacklist_filtered.bam(.bai)  # when blacklist configured
│   ├── 03_bigwig/<sample>.CPM.bw
│   ├── 03_signal/<sample>.FE.bdg        # when qc.signal_tracks: true
│   ├── 04_peaks/<sample>/               # MACS3 peak directory
│   └── logs/
│
├── experiments/<experiment>/
│   ├── 01_qc/
│   │   └── <experiment>.pooled_qc_summary.tsv
│   ├── 02_align/
│   │   ├── biorep<N>.final.bam(.bai)
│   │   ├── <experiment>.pooled.final.bam(.bai)
│   │   └── <experiment>.pooled.control.final.bam(.bai)
│   ├── 03_signal/
│   │   ├── <experiment>.pooled.FE.bdg
│   │   └── <experiment>.pooled.ppois.bdg
│   ├── 04_peaks/
│   │   ├── pooled/<experiment>_pooled_peaks/
│   │   └── idr/
│   │       ├── <experiment>_biorep<N>_idr_peaks.narrowPeak
│   │       ├── <experiment>_biorep<N>_pr1_idr_peaks.narrowPeak
│   │       └── ...
│   ├── 05_pseudorep/
│   │   └── <experiment>_biorep<N>.pr1.bam(.bai)
│   ├── 06_idr/
│   │   ├── true_replicates/
│   │   ├── self_pseudoreplicates/
│   │   ├── pooled_pseudoreplicates/
│   │   └── final/
│   │       ├── conservative.narrowPeak
│   │       ├── optimal.narrowPeak
│   │       └── reproducibility_summary.tsv
│   └── logs/
│
└── multiqc/
    ├── multiqc_report.html
    └── stage3_qc_summary.tsv
```

## Limitations

- **TF ChIP-seq IDR** (`stage5: true`) requires `chipseq` assay, `narrowPeak`
  mode, and exactly 2 biological replicates. CUT&Tag IDR and 3+ biorep IDR
  are not yet supported.
- **ATAC-seq analysis is not included.** This pipeline covers ChIP-seq and
  CUT&Tag only.
- **Histone broad-peak IDR** is deferred to a future release.
- **MultiQC integration** for experiment-level IDR and pooled QC summaries is
  not yet complete.
- **`bamCoverage` RPGC normalization** requires `--effectiveGenomeSize` and is
  not yet wired up in `tool_parameters`.
- `final.bam` is symlinked to the filtered BAM when `remove_dup` is `no` or
  `auto+broad` — this is correct behavior but worth noting for downstream
  consumers that expect a deduplicated BAM.

See [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md) for the full roadmap and planned
follow-ups.

## Developer Notes

### Specs and schemas

Detailed design specs for each major feature live under
[`docs/superpowers/specs/`](docs/superpowers/specs/). Config and sample sheet
contracts are in [`workflow/schemas/`](workflow/schemas/).

### Legacy single-sample script

`scripts/chipseq.sh` is retained for compatibility. It supports single-sample
PE/SE ChIP-seq and CUT&Tag with command-line flags. For new analyses, the
Snakemake batch workflow (`workflow/Snakefile`) is the recommended entry point.

### Repository layout

- `workflow/Snakefile` — entry point, validation, assay dispatch
- `workflow/rules/` — Snakemake rule files (common, peaks, replicates, IDR, QC, report)
- `workflow/envs/chipseq.yml` — Conda environment
- `workflow/schemas/` — human-readable config and sample sheet contracts
- `config/` — default config and sample sheet
- `scripts/` — validation, QC helpers, and analysis scripts
- `test/` — stress tests for validation, DAG, and helper scripts

### Smoke-test profiles

A suite of 7 test profiles under `test/profiles/` covers the major
ChIP-seq and CUT&Tag dispatch paths (SE / PE / control_sample /
control_bam / SEACR / Stage 5 IDR) via dry-run only.  Run with:

```bash
SNAKEMAKE=/path/to/snakemake python3 test/test_stage8_smoke_profiles.py
```

All temporary files stay under `/tmp` — nothing is written into the
repository.  See `docs/superpowers/specs/2026-05-18-stage8a-test-profiles-smoke-design.md`
for the full design.

### Tiny real execution

One real end-to-end run through preprocessing + signal on synthetic
fixtures (20 kb pseudo-random reference, 1 000 PE reads, Bowtie2 index
built at runtime — zero binary files committed):

```bash
SNAKEMAKE=/path/to/snakemake python3 test/test_stage8b_tiny_execution.py
```

All outputs land under `/tmp`. Exit code 0 = PASS, 1 = FAIL, 2 = SKIP.
MACS3 is intentionally skipped (Stage 8a dry-run covers the DAG).
See `docs/superpowers/specs/2026-05-18-stage8b-tiny-real-execution-design.md`
for the full design.

## License

This project is open-source under the [MIT License](LICENSE).
