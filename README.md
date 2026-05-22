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

This is a **v0.1.0-beta research workflow release**: functional and tested
with smoke and tiny-execution profiles, but not a fully ENCODE-compliant
production pipeline. See [Limitations](#limitations) for known gaps.

## Key Features

- **Shared preprocessing:** FastQC, Trim Galore, Bowtie2 alignment, MAPQ
  filtering, samtools duplicate handling, flagstat, idxstats, BigWig generation
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
micromamba create -f workflow/envs/runner.yml
micromamba activate chipseq-runner
```

The `chipseq-runner` environment is intentionally small: it contains Python,
PyYAML, and Snakemake. Bioinformatics tools are installed as rule-specific
Conda environments when you run with `--use-conda`.
See [docs/environments.md](docs/environments.md) for the full environment
layout, cache behavior, and cleanup commands.

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

> For troubleshooting, smoke-test instructions, and detailed run guidance, see
> [docs/quickstart.md](docs/quickstart.md).

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

### Example

```tsv
sample	fastq_1	fastq_2	layout	assay	target	peak_mode	genome	bowtie2_index	experiment	biological_replicate
H3K27AC_rep1	/data/ac1_R1.fq.gz	/data/ac1_R2.fq.gz	PE	chipseq	H3K27ac	narrow	hs	/path/to/bt2/GRCh38	H3K27AC	1
H3K27AC_rep2	/data/ac2_R1.fq.gz	/data/ac2_R2.fq.gz	PE	chipseq	H3K27ac	narrow	hs	/path/to/bt2/GRCh38	H3K27AC	2
```

Optional replicate and control columns (`experiment`, `condition`, `replicate`,
`biological_replicate`, `technical_replicate`, `role`, `control_sample`,
`control_bam`) are documented in `workflow/schemas/samples.schema.yaml`.

> Full column reference, role/control semantics, additional examples
> (treatment+control, biological replicates), and common pitfalls:
> [docs/sample-sheet.md](docs/sample-sheet.md)

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

> Full configuration reference covering core keys, genome resources, QC
> switches, replicates/IDR, CUT&Tag SEACR, tool parameters, `use_control`,
> `multiqc`, and key dependencies:
> [docs/configuration.md](docs/configuration.md)

## Workflow Capabilities

### Preprocessing and peak calling

Every sample (treatment or control) flows through: FastQC → Trim Galore
(or symlink when `trim: false`) → Bowtie2 alignment → samtools sort/index →
MAPQ filter → duplicate handling (samtools by default; Picard if available in
a custom runtime environment) → `final.bam`.
BigWig tracks are generated with deepTools `bamCoverage` (CPM-normalized by
default). MACS3 peak calling runs on treatment samples, with assay-specific
parameters (TF ChIP-seq model-based, CUT&Tag Tn5-aware `--shift -100`).

### Controls

Controls are disabled by default (`use_control: false`). When enabled,
each treatment row may reference a `control_sample` (another sample row with
`role=control`) or an external `control_bam` path. Control rows are processed
through the full preprocessing pipeline and their `final.bam` is passed as
MACS3 `-c`. See [docs/sample-sheet.md](docs/sample-sheet.md) for role and
control semantics.

### Single-sample QC

When the `qc` block is enabled, each treatment sample receives:
- **Blacklist filtering** (BAM + peaks) when a blacklist BED is configured
- **FRiP** (Fraction of Reads in Peaks)
- **Library complexity** (Picard-derived metrics when available; fallback
  status otherwise)
- **NRF/PBC** (BAM-derived library complexity)
- **MACS3 signal tracks**: fold-enrichment (`FE.bdg`) and p-value
  (`ppois.bdg`) bedGraph from `macs3 bdgcmp`
- **Cross-correlation** (opt-in): phantompeakqualtools NSC/RSC and fragment
  length (`qc.cross_correlation: true`)
- **Preseq library complexity** (opt-in): `preseq lc_extrap` extrapolation
  curve (`qc.preseq_complexity: true`)
- **Picard CollectMultipleMetrics** (opt-in): alignment summary, insert size,
  and quality distribution (`qc.picard_metrics: true`; requires
  `reference_fasta` with matching `.fai` and `.dict`; uses
  `VALIDATION_STRINGENCY=LENIENT` for QC-only metrics on filtered PE BAMs)
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
See [docs/configuration.md](docs/configuration.md) for `stage4b` and `stage5`
gating.

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

### Assay support

- **ChIP-seq** and **CUT&Tag** are supported. ATAC-seq is not included.

### TF ChIP-seq IDR

- Requires `chipseq` assay, `narrowPeak` mode, and **exactly 2 treatment
  biological replicates** per experiment.
- 3+ replicate IDR is not yet supported (automatic pairwise selection
  among 3+ replicates is not implemented).
- CUT&Tag IDR is not yet supported.

### Histone broad marks

- Broad-peak IDR is not yet supported. Histone experiments benefit from
  pooled QC summaries (Stage 6b) but do not produce IDR peak sets.

### QC gaps

- GC bias metrics (Picard CollectGcBiasMetrics) are not yet implemented.
- FE/ppois bedGraph tracks are produced; BigWig conversion for those
  tracks is not yet implemented.
- `bamCoverage` RPGC normalization requires `--effectiveGenomeSize` and
  is not yet wired up in `tool_parameters`.
- MultiQC integration for experiment-level IDR and pooled QC summaries
  is not yet complete.

### Data and reproducibility

- No real public dataset is bundled with the pipeline.
- When duplicate removal is disabled (`remove_dup: "no"` or
  `auto` + `broad` peak mode), `final.bam` is the workflow's final
  post-filter BAM and may be a symlink to the filtered BAM. It should
  not be assumed to be duplicate-removed in all modes.

See [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md) for the full roadmap and planned
follow-ups.

## Developer Notes

### Specs and schemas

Detailed design specs for each major feature live under
[`docs/superpowers/specs/`](docs/superpowers/specs/). Config and sample sheet
contracts are in [`workflow/schemas/`](workflow/schemas/).
Path handling rules are documented in
[`docs/developer/no-hardcoding.md`](docs/developer/no-hardcoding.md).

### Legacy single-sample script

`scripts/chipseq.sh` is retained for compatibility. It supports single-sample
PE/SE ChIP-seq and CUT&Tag with command-line flags. For new analyses, the
Snakemake batch workflow (`workflow/Snakefile`) is the recommended entry point.

### Repository layout

- `workflow/Snakefile` — entry point, validation, assay dispatch
- `workflow/rules/` — Snakemake rule files (common, peaks, replicates, IDR, QC, report)
- `workflow/envs/` — Conda environments (lightweight runner plus rule-specific tool envs; see [docs/environments.md](docs/environments.md))
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
for the full design. For troubleshooting and test execution details, see
[docs/quickstart.md](docs/quickstart.md).

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

### CI

A GitHub Actions workflow runs on every PR and push to `main` / `stage*`.
PR/push CI uses a lightweight `ci-fast` environment (python + pyyaml +
snakemake only) for validation and dry-run checks:

- Validate default config + sample sheet
- Run validation stress tests
- Run Stage 8a dry-run smoke profiles (7 profiles)

A manual `workflow_dispatch` job runs the tiny real-execution harness
using the core `chipseq` environment.  See `.github/workflows/ci.yml`
and the `workflow/envs/ci-fast.yml` minimal environment.
Design: `docs/superpowers/specs/2026-05-18-stage8c-github-actions-ci-design.md`

For local environment layout, first-run behavior, and cleanup commands, see
[docs/environments.md](docs/environments.md).

### Local execution

**Prerequisites:**

```bash
# Create the lightweight runner environment
micromamba create -f workflow/envs/runner.yml
micromamba activate chipseq-runner
```

**Validation and dry-run (fast, no data needed):**

```bash
python3 scripts/validate_samples.py --config config/config.yaml
python3 test/test_validation_stress.py
python3 test/test_stage8_smoke_profiles.py
```

**Tiny real execution (requires the core chipseq env):**

```bash
micromamba create -f workflow/envs/chipseq.yml
micromamba activate chipseq
python3 test/test_stage8b_tiny_execution.py
```

**Full workflow run:**

```bash
# 1. Edit config/config.yaml and config/samples.tsv for your data.
# 2. Dry-run:
snakemake -s workflow/Snakefile --configfile config/config.yaml -n --use-conda

# 3. Execute (replace N with core count):
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores N --use-conda
```

On a workstation, use as many cores as available.  On a shared server,
limit cores and consider `--latency-wait` for NFS filesystems.

### Release checklist

Before tagging a release, run these checks locally:

```bash
# 1. Validation
python3 scripts/validate_samples.py --config config/config.yaml
python3 test/test_validation_stress.py

# 2. Default DAG check
snakemake -s workflow/Snakefile --configfile config/config.yaml -n --quiet

# 3. Dry-run smoke profiles (7 profiles, <30 s)
python3 test/test_stage8_smoke_profiles.py

# 4. Tiny real execution (preprocessing + signal, <60 s)
python3 test/test_stage8b_tiny_execution.py

# 5. Repo hygiene
git status --short --untracked-files=all
# Expect: no results/, .snakemake/, *.fq, *.fq.gz, *.bam, *.bai, *.bw

# 6. CI status
# Check GitHub Actions for green on main / current branch.
```

## License

This project is open-source under the [MIT License](LICENSE).
