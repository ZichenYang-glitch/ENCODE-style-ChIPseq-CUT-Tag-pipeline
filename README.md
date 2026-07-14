# HelixWeave

**Reproducible omics workflows, from inputs to evidence.**

[![Snakemake](https://img.shields.io/badge/Snakemake-%3E%3D8.0-brightgreen.svg?style=flat-square)](https://snakemake.github.io)
[![Conda](https://img.shields.io/badge/conda-supported-blue.svg?style=flat-square)](https://docs.conda.io/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square)](https://opensource.org/licenses/MIT)

## Overview

HelixWeave is a workflow-neutral platform for local and small-team omics
analysis. It provides schema-driven input authoring, validation, durable run
tracking, execution, and evidence review without making the platform depend on
one workflow's rules or file layout.

This repository also ships HelixWeave's first scientific adapter: the existing
ENCODE-style ChIP-seq, CUT&Tag, ATAC-seq, and MNase-seq Snakemake workflow.
The adapter and workflow retain their established identities and scientific
behavior while the surrounding product uses the HelixWeave brand.

## Start the local platform

After installing the locked local environment and frontend dependencies, check
the toolchain and start the supervised stack:

```bash
python scripts/run_local_platform.py --doctor
python scripts/run_local_platform.py
```

Open `http://127.0.0.1:5173`. For the deterministic input-to-results path, use
`python scripts/run_local_platform.py --input-authoring-demo`. See
[the local runtime guide](docs/development/local-platform-runtime.md) for setup,
process ownership, storage, and cleanup details.

## Brand and compatibility

This migration changes the product display layer only. Stable scientific and
runtime identities remain compatible.

| Surface | Identity after this change | Compatibility status |
| :--- | :--- | :--- |
| Display brand | `HelixWeave` | Renamed in product UI, browser metadata, docs, and OpenAPI metadata. |
| Repository slug | `ENCODE-style-ChIPseq-CUT-Tag-pipeline` | Unchanged; renaming the GitHub repository requires separate authorization. |
| Python distribution | `encode-pipeline` | Unchanged. |
| Python import namespace | `encode_pipeline` | Unchanged. |
| CLI names | `encode-validate`, `encode-manifest`, `encode-dag`, `encode-worker` | Unchanged. |
| API title and routes | `HelixWeave API`; `/api/v1` | Display metadata renamed; routes and payload contracts unchanged. |
| Workflow and adapter identity | `encode-style-chipseq-cuttag-atac-mnase`; ENCODE-style adapter | Unchanged. |
| Persistence identity | Existing SQLite tables, Alembic history, environment variables, Redis/RQ job identity, and artifact URIs | Unchanged. |
| Frontend package | `helixweave-frontend` | Renamed only in the private, unpublished frontend workspace. |

## Bundled ENCODE-style workflow

A Snakemake-based pipeline suite for ChIP-seq, CUT&Tag, ATAC-seq, and MNase-seq
data analysis. It handles single-sample preprocessing as well as
multi-replicate experiments with pooled outputs, single-sample QC, and TF
ChIP-seq IDR reproducibility analysis.

Default mode is no-input / no-control. Optional controls can be enabled with
`use_control: true` and supplied as an external control BAM or a FASTQ-based
control sample row. Dependencies are managed with Conda.

The latest published pre-release is **v0.2.0-rc1**. The `main` branch includes
substantial post-rc1 changes, including MNase nucleosome positioning,
artifact inventory and contract test infrastructure, Snakefile extraction,
paths.smk MNase path helpers, and an artifact adoption decision record.
See [Limitations](#limitations) for known gaps.

## Key Features

- **Shared preprocessing:** FastQC, Trim Galore, Bowtie2 alignment, MAPQ
  filtering, samtools duplicate handling, flagstat, idxstats, BigWig generation
- **ChIP-seq / CUT&Tag / ATAC-seq assay policies:** assay-aware MACS3
  parameters, duplicate removal, and read extension; optional CUT&Tag SEACR
  sidecar peak calls
  (`cuttag.seacr.enabled`, output under `results/<sample>/04_peaks_seacr/`)
- **MNase-seq nucleosome positioning:** PE-only MNase-seq support with
  sub/mono/di-nucleosome fragment BAMs, dyad BigWig (`bamCoverage --MNase`),
  mono occupancy BigWig, per-sample MNase QC summary, and pooled outputs for
  multi-replicate experiments. Fragment ranges and dyad window are configurable
  via `mnase.fragments` and `mnase.dyad_range`. Caller config surface reserved;
  nucleosome calling (DANPOS3/iNPS/SEM) deferred to v0.3.
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
micromamba create -n chipseq-runner --file workflow/envs/runner.lock
micromamba activate chipseq-runner
```

Or with conda:

```bash
conda create -n chipseq-runner --file workflow/envs/runner.lock
conda activate chipseq-runner
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
python3 scripts/validate_samples.py --config config/config.yaml --strict-inputs  # optional: validate file existence
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
| `assay` | `chipseq`, `cuttag`, `atac`, or `mnase`. |
| `target` | Antibody or target name (e.g. `H3K27ac`, `CTCF`). |
| `peak_mode` | `narrow`, `broad`, or `nucleosome` (`nucleosome` is MNase-only). |
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
Other fields are optional paths unless the feature that uses them is enabled
(`gtf` for TSS profiles, `reference_fasta` for Picard metrics). If non-empty
they must exist on disk.
See [docs/reference-resources.md](docs/reference-resources.md) for preparing
Bowtie2 indexes, FASTA `.fai`/`.dict`, chrom sizes, blacklists, and annotations.

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
default). MACS3 peak calling runs on peak-centric treatment samples (ChIP-seq,
CUT&Tag, ATAC-seq), with assay-specific parameters (TF ChIP-seq model-based,
CUT&Tag Tn5-aware `--shift -100`). MNase-seq samples skip MACS3 and instead
produce nucleosome-centric outputs (mono-nucleosome BAM, dyad BigWig, mono
occupancy BigWig).

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
  length (`qc.cross_correlation: true`); project-level summary at
  `results/multiqc/cross_correlation_summary.tsv` and a MultiQC custom section
- **Preseq library complexity** (opt-in): `preseq lc_extrap` extrapolation
  curve (`qc.preseq_complexity: true`)
- **Picard CollectMultipleMetrics** (opt-in): alignment summary, insert size,
  and quality distribution (`qc.picard_metrics: true`; requires
  `reference_fasta` with matching `.fai` and `.dict`; uses
  `VALIDATION_STRINGENCY=LENIENT` for QC-only metrics on filtered PE BAMs)
- **TSS enrichment-style profiles** (opt-in): deepTools matrix/profile around
  transcript TSSs (`qc.tss_enrichment: true`; requires
  `genome_resources.<genome>.gtf`)
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
│   ├── 05_qc/                           # opt-in cross-correlation, preseq,
│   │                                    # Picard metrics, and TSS profiles
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

- **ChIP-seq**, **CUT&Tag**, baseline **ATAC-seq**, and **PE MNase-seq** are supported.
- ATAC-seq currently supports `peak_mode: narrow` only. ATAC-specific
  insert-size/TSS interpretation is available as QC output, but no
  ATAC-specific footprinting or nucleosome-positioning module is included.
- MNase-seq is PE-only and supports sub/mono/di-nucleosome fragment BAMs,
  dyad BigWig, mono occupancy BigWig, per-sample QC summary, and pooled outputs.
  Nucleosome calling (DANPOS3/iNPS/SEM) and pooled sub/di fragment outputs are
  deferred to v0.3.

### TF ChIP-seq IDR

- Requires `chipseq` assay, `narrowPeak` mode, and **exactly 2 treatment
  biological replicates** per experiment.
- 3+ replicate IDR is not yet supported (automatic pairwise selection
  among 3+ replicates is not implemented).
- CUT&Tag and ATAC-seq IDR are not yet supported.

### Histone broad marks

- Broad-peak IDR is not yet supported. Histone experiments benefit from
  pooled QC summaries (Stage 6b) but do not produce IDR peak sets.

### QC gaps

- GC bias metrics (Picard CollectGcBiasMetrics) are not yet implemented.
- FE/ppois bedGraph tracks are produced; BigWig conversion for those
  tracks is available when `genome_resources.<genome>.chrom_sizes` is configured.
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

See [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md) for known scientific and operational
follow-ups.

## Developer Notes

### Architecture, roadmap, and contract docs

The maintained product roadmap is in
[`docs/development/workflow-platform-agent-roadmap.md`](docs/development/workflow-platform-agent-roadmap.md),
and the durable system boundaries are summarized in
[`docs/architecture/platform-overview.md`](docs/architecture/platform-overview.md).
[`ROADMAP_v0.2.md`](ROADMAP_v0.2.md) is a superseded scientific-release
snapshot.
The release checklist is in [`RELEASE_CHECKLIST.md`](RELEASE_CHECKLIST.md).
The output contract and manifest schema are in [`docs/output-contract.md`](docs/output-contract.md).
Assay-specific behavioral contracts are in [`docs/assay-policy.md`](docs/assay-policy.md).
The IDR reproducibility contract is in [`docs/idr-contract.md`](docs/idr-contract.md).
Public data validation plans are under [`docs/release-checks/`](docs/release-checks/):
[`stage27-public-data-validation-plan.md`](docs/release-checks/stage27-public-data-validation-plan.md),
[`stage27b-metadata-ci-plan.md`](docs/release-checks/stage27b-metadata-ci-plan.md).

### Schemas and maintained decisions

Config and sample sheet contracts are in
[`workflow/schemas/`](workflow/schemas/). Detailed maintained architecture
decisions are linked from the platform overview; completed implementation plans
remain available through Git history.
Path handling rules are documented in
[`docs/developer/no-hardcoding.md`](docs/developer/no-hardcoding.md).

### Container usage

Docker and Apptainer/SingularityCE container images are available for the
runner environment. See [`docs/container-usage.md`](docs/container-usage.md)
for build, run, bind-mount, and troubleshooting instructions.

### Legacy single-sample script (deprecated)

`scripts/chipseq.sh` is deprecated. It now prints an error message directing
users to the Snakemake workflow. The historical full script is archived at
`docs/archive/scripts/chipseq-legacy.sh` for reference.

The canonical entry point is:

```
snakemake -s workflow/Snakefile --configfile config/config.yaml
```

### QC and MultiQC

See [docs/qc-interpretation.md](docs/qc-interpretation.md) for a
comprehensive guide to interpreting every QC metric the pipeline produces
(FastQC, alignment, library complexity, peak quality, cross-correlation,
Picard metrics, TSS profiles, CUT&Tag-specific QC, ATAC-seq QC, and
replicate-level outputs).

### Repository layout

- `workflow/Snakefile` — entry point, validation, assay dispatch
- `workflow/rules/` — Snakemake rule files (common, metadata, peaks, replicates, IDR, QC, report, targets, paths, and per-assay files)
- `workflow/lib/` — shared Python libraries (artifact dataclass, loader, query helpers)
- `workflow/envs/` — Conda environments (lightweight runner plus rule-specific tool envs; see [docs/environments.md](docs/environments.md))
- `workflow/schemas/` — human-readable config and sample sheet contracts
- `config/` — default config and sample sheet
- `scripts/` — validation, QC helpers, and analysis scripts
- `test/` — behavior and contract tests for validation, DAG, and helper scripts

### Smoke-test profiles

A suite of test profiles under `test/profiles/` covers major ChIP-seq and
CUT&Tag dispatch paths (SE / PE / control sample / control BAM / SEACR / IDR)
via dry-run only. Focused contracts cover ATAC, TSS, MNase, artifact inventory,
manifest behavior, and output paths. Run with:

```bash
SNAKEMAKE=/path/to/snakemake python3 -m pytest test/workflow/test_smoke_profiles.py -v
```

All temporary files stay under `/tmp` — nothing is written into the
repository. For troubleshooting and test execution details, see
[docs/quickstart.md](docs/quickstart.md).

### Real execution

The real-execution tier combines focused samtools contracts with one
end-to-end preprocessing + signal run on synthetic fixtures (20 kb
pseudo-random reference, 1 000 PE reads, Bowtie2 index built at runtime — zero
binary files committed):

```bash
SNAKEMAKE=/path/to/snakemake python3 -m pytest -m real_execution \
  test/real_execution -v
```

All outputs land under pytest-managed temporary directories. Missing external
tools are reported as an explicit skip. MACS3 is intentionally omitted because
the dry-run profiles cover its DAG.
See [the real-execution harness](docs/development/real-execution-harness.md)
for tier ownership and prerequisites.

### CI

A GitHub Actions workflow runs on every PR and push to `main` / `stage*`.
PR/push CI uses the locked `ci-fast` environment and currently covers:

- Validate default config + sample sheet
- Durable execution, persistence, API, and worker contracts
- Config validation and dry-run smoke profiles
- No-hardcoded-paths guard
- BigWig, QC summary, manifest, complete Python coverage, and changed-lines checks
- OpenAPI/generated-client drift, frontend tests, typecheck, and build
- Critical browser execution journeys

A manual `workflow_dispatch` job runs the tiny real-execution harness
using the core `chipseq` environment.  See `.github/workflows/ci.yml`
and the `workflow/envs/ci-fast.lock` environment.

For local environment layout, first-run behavior, and cleanup commands, see
[docs/environments.md](docs/environments.md).

### Local execution

**Prerequisites:**

```bash
# Create the lightweight runner environment
micromamba create -n chipseq-runner --file workflow/envs/runner.lock
micromamba activate chipseq-runner
```

**Validation and dry-run (fast, no data needed):**

```bash
python3 scripts/validate_samples.py --config config/config.yaml
python3 -m pytest test/config/test_validation.py -v
python3 -m pytest test/workflow/test_smoke_profiles.py -v
```

**Tiny real execution (requires the core chipseq env):**

```bash
micromamba create -n chipseq --file workflow/envs/chipseq.lock
micromamba activate chipseq
python3 -m pytest -m real_execution \
  test/real_execution -v
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

Before tagging a release, follow the maintained
[`RELEASE_CHECKLIST.md`](RELEASE_CHECKLIST.md). Tagging and publishing require
explicit authorization.

## License

This project is open-source under the [MIT License](LICENSE).
