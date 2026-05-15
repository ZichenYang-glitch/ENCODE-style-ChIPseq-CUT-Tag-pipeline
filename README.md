# ENCODE-style ChIP-seq and CUT&Tag Pipeline

[![Snakemake](https://img.shields.io/badge/Snakemake-%3E%3D8.0-brightgreen.svg?style=flat-square)](https://snakemake.github.io)
[![Conda](https://img.shields.io/badge/conda-supported-blue.svg?style=flat-square)](https://docs.conda.io/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=flat-square)](https://opensource.org/licenses/MIT)

This repository contains a Snakemake-based pipeline suite for ChIP-seq and
CUT&Tag data analysis. The batch workflow is the primary entry point and is
implemented as modular rules under `workflow/rules/`.

Default mode is no-input / no-control. Optional controls can be enabled with
`use_control: true` and supplied either as an external control BAM or as a
FASTQ-based control sample row.

Dependencies are managed with Conda. This pipeline does not include ATAC-seq
analysis.

## Key Features

- ChIP-seq and CUT&Tag assay policies with shared preprocessing rules.
- Snakemake-native DAG tracking, resume support, and per-step reruns.
- Optional control handling through external BAMs or FASTQ-based control rows.
- Conda-managed workflow environment.
- Legacy single-sample shell script retained for compatibility.

## Repository Structure

- `workflow/Snakefile`: Main Snakemake entry point, sample validation, and assay dispatch.
- `workflow/rules/common.smk`: Shared rules for QC, trimming, alignment, filtering, duplicate handling, QC stats, and BigWig generation.
- `workflow/rules/chipseq.smk`: ChIP-seq policy functions.
- `workflow/rules/cuttag.smk`: CUT&Tag policy functions.
- `workflow/rules/peaks.smk`: MACS3 peak calling rule.
- `workflow/rules/report.smk`: Completion sentinels and MultiQC aggregation.
- `workflow/schemas/config.schema.yaml`: Config schema contract (human-readable).
- `workflow/schemas/samples.schema.yaml`: Sample sheet schema contract (human-readable).
- `workflow/envs/chipseq.yml`: Conda environment for the Snakemake workflow.
- `config/config.yaml`: Workflow configuration.
- `config/samples.tsv`: Sample sheet.
- `scripts/validate_samples.py`: Config and sample sheet validator (importable + CLI).
- `scripts/chipseq.sh`: Legacy single-sample script kept for compatibility.
- `KNOWN_ISSUES.md`: Non-blocking follow-ups and current limitations.

## Quick Start: Snakemake Batch Workflow

### 1. Installation

```bash
git clone https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline.git
cd ENCODE-style-ChIPseq-CUT-Tag-pipeline

conda env create -f workflow/envs/chipseq.yml
conda activate chipseq
```

### 2. Configure Samples

Edit `config/samples.tsv`:

```tsv
sample	fastq_1	fastq_2	layout	assay	target	peak_mode	genome	bowtie2_index	control_bam
H3K27AC	/data/ac_R1.fq.gz	/data/ac_R2.fq.gz	PE	chipseq	H3K27ac	narrow	mm	/path/to/bt2/mm10
CTCF_rep1	/data/ctcf_R1.fq.gz		SE	chipseq	CTCF	narrow	hs	/path/to/bt2/hg38
CUTTAG_H3K4me3	/data/ct_R1.fq.gz	/data/ct_R2.fq.gz	PE	cuttag	H3K4me3	narrow	mm	/path/to/bt2/mm10
```

Unused optional fields can be left blank. For example, the `control_bam`
column is empty above because `use_control` is normally `false`.

Required columns:

| Column | Description |
| :--- | :--- |
| `sample` | Unique sample ID. Allowed characters: letters, numbers, `_`, `.`, `-`. |
| `fastq_1` | R1 FASTQ path. |
| `fastq_2` | R2 FASTQ path. Required for `PE`; leave empty for `SE`. |
| `layout` | `PE` or `SE`. |
| `assay` | `chipseq` or `cuttag`. |
| `target` | Antibody or target name. |
| `peak_mode` | `narrow` or `broad`. |
| `genome` | Genome shortcut or MACS3 genome size, for example `mm`, `mm10`, `hs`, `hg38`. |
| `bowtie2_index` | Bowtie2 index prefix. |
| `control_bam` | Optional external control BAM. Used only when `use_control: true`. |

Optional columns for FASTQ-based controls:

| Column | Default | Description |
| :--- | :--- | :--- |
| `role` | `treatment` | `treatment` or `control`. Peak calling runs only for treatment rows. |
| `control_sample` | empty | Sample ID of a control row to use as MACS3 input control. Used only when `use_control: true`. |

### 3. Configure Workflow Options

Edit `config/config.yaml`:

```yaml
samples: "config/samples.tsv"
outdir: "results"
threads: 8
mapq: 30
binsize: 10
remove_dup: "auto"     # auto, yes, no
trim: true             # true or false
extend_reads: "auto"   # auto, yes, no, or a positive integer
use_control: false     # false means no control_bam and no control_sample
multiqc: true
```

### 4. Validate Config and Samples

Validate configuration and sample sheet before running the workflow:

```bash
# Run the standalone validator.
python3 scripts/validate_samples.py --config config/config.yaml

# Or let the Snakefile validate at parse time — invalid configs will exit
# with a clear error before the DAG is built.
snakemake -s workflow/Snakefile --configfile config/config.yaml -n
```

The validator is also importable by the Snakefile, so validation runs
automatically at parse time. Schema contracts in `workflow/schemas/`
document the expected config and sample sheet format.

### 5. Run

```bash
# Dry-run: build the DAG without executing commands.
snakemake -s workflow/Snakefile --configfile config/config.yaml -n --use-conda

# Run with 16 cores. --use-conda creates/uses workflow/envs/chipseq.yml.
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 16 --use-conda

# Resume incomplete jobs after an interrupted run.
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 16 --use-conda --rerun-incomplete
```

## Control/Input Handling

Controls are disabled by default. With `use_control: false`, both
`control_bam` and `control_sample` are ignored.

To use a precomputed control BAM:

```yaml
use_control: true
```

```tsv
sample	fastq_1	fastq_2	layout	assay	target	peak_mode	genome	bowtie2_index	control_bam
H3K27AC	/data/ac_R1.fq.gz	/data/ac_R2.fq.gz	PE	chipseq	H3K27ac	narrow	mm	/path/to/bt2/mm10	/data/input.final.bam
```

To process an input/control sample from FASTQ, add `role` and
`control_sample` columns:

```tsv
sample	fastq_1	fastq_2	layout	assay	target	peak_mode	genome	bowtie2_index	control_bam	role	control_sample
H3K27AC	/data/ac_R1.fq.gz	/data/ac_R2.fq.gz	PE	chipseq	H3K27ac	narrow	mm	/path/to/bt2/mm10		treatment	Input_rep1
Input_rep1	/data/input_R1.fq.gz	/data/input_R2.fq.gz	PE	chipseq	Input	narrow	mm	/path/to/bt2/mm10		control
```

Control rows are processed through the shared preprocessing rules and produce
`final.bam`. MACS3 peak calling is scheduled only for treatment rows. Do not
set both `control_bam` and `control_sample` on the same treatment row.

## Genome Resources (Stage 2+)

The `genome_resources` config block maps genome labels to effective genome
sizes and optional resource paths. This foundation will be used by later
ENCODE-like QC stages (blacklist filtering, FRiP, library complexity, etc.).

**Important distinction:** `hs` and `mm` are **MACS3 effective genome size
shortcuts**, not genome assemblies. They map to:

| Shortcut | MACS3 -g value | Corresponding assembly | Numeric size |
| :--- | :--- | :--- | :--- |
| `hs` | `hs` | GRCh38 / hg38 | 2,913,022,398 |
| `mm` | `mm` | GRCm38 / mm10 | 2,652,783,500 |

**Explicit numeric entries** are also provided:

| Label | Numeric effective genome size |
| :--- | :--- |
| `hg19` | 2,864,785,220 |
| `hg38` | 2,913,022,398 |
| `mm10` | 2,652,783,500 |
| `mm39` | 2,654,621,783 |

For **mm39 / GRCm39**, use the explicit `mm39` entry (or set a numeric
`effective_genome_size`). Do not pass `mm` for mm39 — `mm` maps to the
mm10 size.

When `genome_resources` is configured, `_normalize_genome()` in the
Snakefile prefers the configured `effective_genome_size` over legacy
mappings. The old fallback (`mm39` → `mm`, `hg38` → `hs`, etc.) only
applies when no genome resource entry exists for the sample's genome.

Example configuration (see `config/config.yaml` for defaults):

```yaml
genome_resources:
  hs:
    effective_genome_size: "hs"
    chrom_sizes: ""
    blacklist: ""
    gtf: ""
    reference_fasta: ""

  hg38:
    effective_genome_size: 2913022398
    chrom_sizes: "/data/genomes/hg38/hg38.chrom.sizes"
    blacklist: "/data/genomes/hg38/hg38.blacklist.bed"
    gtf: ""
    reference_fasta: ""
```

## Workflow Steps

1. FastQC on raw FASTQs.
2. Trim Galore, or raw FASTQ symlinks when `trim: false`.
3. Bowtie2 alignment with read group tags.
4. Samtools indexing and MAPQ filtering.
5. Picard duplicate metrics and optional duplicate removal.
6. Samtools `flagstat` and `idxstats`.
7. deepTools `bamCoverage` CPM BigWig generation.
8. MACS3 peak calling for treatment samples.
9. MultiQC aggregation over the current sample directories.

## Assay-Specific Policy

- ChIP-seq uses standard MACS3 parameters. Broad mode adds `--broad --broad-cutoff 0.1`.
- CUT&Tag narrow mode adds Tn5-aware MACS3 parameters: `--nomodel --shift -100 --extsize 200`.
- CUT&Tag broad mode follows the broad MACS3 policy.
- Stage 1 keeps duplicate removal and read extension behavior aligned with the legacy script.

## Stage 3 QC: Single-Sample Quality Metrics

Stage 3a adds ENCODE-like single-sample QC metrics. All features are
resource-gated — blacklist operations only run when a blacklist BED is
configured for the sample's genome in `genome_resources`.

### Configuration

QC is controlled by an optional `qc` block in `config/config.yaml`:

```yaml
qc:
  blacklist_filter: true   # bedtools intersect -v filtering for BAM + peaks
  frip: true               # Fraction of Reads in Peaks
  summary: true            # Per-sample QC summary TSV + project-level aggregate
```

Missing `qc` block defaults to all switches enabled. Each switch accepts
boolean or string boolean (`true`/`false`).

### QC Metrics

#### Blacklist Filtering

- BAM: `bedtools intersect -v -abam final.bam -b blacklist.bed`
- Peaks: `bedtools intersect -v -a peaks_file -b blacklist.bed`
- Only scheduled when the sample's genome has a non-empty `blacklist` path
  in `genome_resources`. Samples without a configured blacklist skip
  blacklist-specific outputs.
- `final.bam` is **not** replaced — blacklist-filtered BAM is an additional
  output (`{sample}.blacklist_filtered.bam`).

#### FRiP (Fraction of Reads in Peaks)

Stage 3a FRiP is **read-record based** (not fragment-based):
```
FRiP = reads_in_peaks / total_reads
```
- `total_reads`: `samtools view -c` on BAM
- `reads_in_peaks`: `bedtools intersect -u -abam BAM -b peaks | samtools view -c`
- Uses blacklist-filtered BAM and peaks when **both** are available for
  a sample; otherwise falls back to unfiltered files.
- Output: `{sample}.frip.tsv`

#### Peak Counts

- Counts lines in the MACS3 peak file (narrowPeak or broadPeak).
- If blacklist-filtered peaks exist, reports both raw and filtered counts.
- Output: `{sample}.peak_counts.tsv`

#### QC Summary

Per-sample TSV (`{sample}.qc_summary.tsv`) with columns:
`sample`, `assay`, `target`, `genome`, `layout`, `peak_mode`,
`use_control`, `control_type`, `final_bam`, `peaks`, `blacklist`,
`blacklist_filtered_bam`, `blacklist_filtered_peaks`, `total_reads`,
`reads_in_peaks`, `frip`, `peak_count`, `blacklist_filtered_peak_count`

Unavailable metrics are filled with `NA` (not empty string).

A project-level summary aggregates all per-sample TSVs at
`multiqc/stage3_qc_summary.tsv`.

### QC Output Structure

```text
results/<sample>/01_qc/
├── <sample>.peak_counts.tsv
├── <sample>.frip.tsv
├── <sample>.qc_summary.tsv
└── ... (existing QC files)

results/<sample>/02_align/
├── <sample>.blacklist_filtered.bam        # only with blacklist
└── <sample>.blacklist_filtered.bam.bai    # only with blacklist

results/<sample>/04_peaks/<sample>/
└── <sample>_peaks.narrowPeak

results/<sample>/04_peaks/<sample>_blacklist_filtered/
└── <sample>_peaks.blacklist_filtered.narrowPeak  # only with blacklist

results/multiqc/
└── stage3_qc_summary.tsv                  # project-level aggregate
```

## Output Structure

The workflow writes per-sample results under `outdir`:

```text
results/
├── <sample>/
│   ├── 00_raw/
│   │   ├── <sample>_R1_val_1.fq.gz
│   │   └── <sample>_R2_val_2.fq.gz
│   ├── 01_qc/
│   │   ├── trim_galore/
│   │   ├── <sample>.flagstat.txt
│   │   ├── <sample>.final.flagstat.txt
│   │   ├── <sample>.idxstats.txt
│   │   └── <sample>.dup_metrics.txt
│   ├── 02_align/
│   │   ├── <sample>.sorted.bam
│   │   ├── <sample>.mapq30.bam
│   │   ├── <sample>.final.bam
│   │   └── <sample>.final.bam.bai
│   ├── 03_bigwig/
│   │   └── <sample>.CPM.bw
│   ├── 04_peaks/
│   │   └── <sample>/
│   └── logs/
│       ├── <sample>.fastqc.done
│       ├── <sample>.trim.done
│       └── <sample>.pipeline.done
└── multiqc/
    └── multiqc_report.html
```

`final.bam` is the stable downstream BAM contract. It points to either the
duplicate-handled BAM or the filtered BAM, depending on `remove_dup`.

## Legacy Single-Sample Script

`scripts/chipseq.sh` remains available for compatibility, but the Snakemake
workflow is the recommended entry point for new analyses.

```bash
# Paired-end ChIP-seq
scripts/chipseq.sh -s MySample -x /path/to/bt2_index -g mm --pe -r1 R1.fq.gz -r2 R2.fq.gz --tech chip

# Paired-end CUT&Tag
scripts/chipseq.sh -s MySample -x /path/to/bt2_index -g mm --pe -r1 R1.fq.gz -r2 R2.fq.gz --tech cuttag

# Single-end ChIP-seq with no trimming
scripts/chipseq.sh -s MySample -x /path/to/bt2_index -g hs --se -r1 sample.fq.gz --tech chip --no-trim
```

Known legacy-script hardening tasks are tracked in `KNOWN_ISSUES.md`.

## Known Limitations

See `KNOWN_ISSUES.md` for non-blocking follow-ups, including schema validation,
plotFingerprint migration, environment cleanup, and legacy script hardening.

## License

This project is open-source under the [MIT License](LICENSE).
