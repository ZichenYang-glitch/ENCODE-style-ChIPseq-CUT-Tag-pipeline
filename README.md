# ENCODE-style ChIP-seq & CUT&Tag Pipeline

This is an automated, high-performance bioinformatic pipeline designed for both **ChIP-seq** and **CUT&Tag** data analysis. It follows ENCODE best practices and features adaptive logic for different sequencing protocols.

**Default mode: no-input / no-control** — the pipeline works without control/input samples. Optional control BAM support can be enabled in the configuration.

This pipeline does **not** include ATAC-seq analysis.

<br>

## Repository Structure

- `chipseq.yml`: Conda environment recipe with all required tools.
- `scripts/chipseq.sh`: Single-sample master script (QC, alignment, peak calling, BigWig).
- `config/`: Snakemake batch configuration (YAML + sample sheet).
- `workflow/Snakefile`: Snakemake batch pipeline.
- `workflow/envs/chipseq.yml`: Conda environment for the Snakemake workflow.
- `.gitignore`: Configured to exclude large genomic files.

<br>

## Quick Start: Single Sample

### 1. Installation

```bash
git clone https://github.com/YangZiChen-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline.git
cd ENCODE-style-ChIPseq-CUT-Tag-pipeline

# Create and activate environment
conda env create -f chipseq.yml
conda activate chipseq
```

### 2. Usage

```bash
# Paired-end ChIP-seq
chipseq.sh -s MySample -x /path/to/bt2_index -g mm --pe -r1 R1.fq.gz -r2 R2.fq.gz --tech chip

# Paired-end CUT&Tag
chipseq.sh -s MySample -x /path/to/bt2_index -g mm --pe -r1 R1.fq.gz -r2 R2.fq.gz --tech cuttag

# Single-end ChIP-seq with no trimming
chipseq.sh -s MySample -x /path/to/bt2_index -g hs --se -r1 sample.fq.gz --tech chip --no-trim

# With optional control/input BAM
chipseq.sh -s MySample -x /path/to/bt2_index -g mm --pe -r1 R1.fq.gz -r2 R2.fq.gz --control_bam /path/to/input.bam
```

<br>

## Script Parameters

| Option | Description | Default |
| :--- | :--- | :--- |
| `-s \| --sample` | Sample ID | **Required** |
| `-x \| --index` | Bowtie2 index prefix | **Required** |
| `-g \| --genome` | Genome build (mm, hs, etc.) | `mm` |
| `--pe / --se` | Paired-end or Single-end | **Required** |
| `--tech` | Protocol: `chip` or `cuttag` | `chip` |
| `--peakMode` | Peak type: `narrow` or `broad` | `narrow` |
| `-p \| --threads` | CPU threads | `8` |
| `-o \| --outdir` | Output directory | `./chip_out` |
| `--trim / --no-trim` | Enable/disable trimming | trim enabled |
| `--removeDup` | Dedup: `auto`, `yes`, `no` | `auto` |
| `--extendReads` | BigWig extension: `auto`, `yes`, `no`, or bp | `auto` |
| `--control_bam` | Optional input/control BAM | (none) |

<br>

## Batch Run (Snakemake)

### 1. Install the workflow environment

```bash
conda env create -f workflow/envs/chipseq.yml
conda activate chipseq
```

### 2. Edit the sample sheet

Edit `config/samples.tsv` with your samples:

```tsv
sample	fastq_1	fastq_2	layout	assay	target	peak_mode	genome	bowtie2_index	control_bam
H3K27AC	/data/R1.fq.gz	/data/R2.fq.gz	PE	chipseq	H3K27ac	narrow	mm	/path/to/bt2/mm10
CTCF_rep1	/data/ctcf_R1.fq.gz		SE	chipseq	CTCF	narrow	hs	/path/to/bt2/hg38
CUTTAG_H3K4me3	/data/ct_R1.fq.gz	/data/ct_R2.fq.gz	PE	cuttag	H3K4me3	narrow	mm	/path/to/bt2/mm10
```

**Columns:**

| Column | Description | Required |
| :--- | :--- | :--- |
| `sample` | Unique sample ID | Yes |
| `fastq_1` | R1 FASTQ path | Yes |
| `fastq_2` | R2 FASTQ path | PE only |
| `layout` | `PE` or `SE` | Yes |
| `assay` | `chipseq` or `cuttag` | Yes |
| `target` | Antibody/target name | Yes |
| `peak_mode` | `narrow` or `broad` | Yes |
| `genome` | Genome build (mm, hs, etc.) | Yes |
| `bowtie2_index` | Bowtie2 index prefix | Yes |
| `control_bam` | Optional input BAM (ignored when `use_control: false`) | No |

### 3. Edit configuration

Edit `config/config.yaml` to set global parameters:

```yaml
samples: "config/samples.tsv"
outdir: "results"
threads: 8
mapq: 30
use_control: false    # default: no-input / no-control
multiqc: true
```

### 4. Run

```bash
# Dry-run (check DAG without executing)
snakemake -s workflow/Snakefile --configfile config/config.yaml -n --use-conda

# Run with 16 cores (--use-conda auto-creates env from workflow/envs/chipseq.yml)
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 16 --use-conda

# Resume incomplete run
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 16 --use-conda --rerun-incomplete

# Alternatively, activate the conda environment first and omit --use-conda:
conda activate chipseq
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 16
```

### 5. Output structure

```
results/
├── <sample>/
│   ├── 00_raw/        # Trimmed FASTQs
│   ├── 01_qc/         # FastQC, flagstat, dup metrics, fingerprint
│   ├── 02_align/      # Sorted, filtered, dedup BAMs
│   ├── 03_bigwig/     # CPM-normalized BigWig
│   ├── 04_peaks/      # MACS3 peak files (.narrowPeak or .broadPeak)
│   └── logs/          # Tool logs + pipeline.done sentinel
└── multiqc/
    └── multiqc_report.html
```

<br>

## Pipeline Workflow

1. **Quality Control**: `FastQC` generates raw data metrics.
2. **Trimming**: `Trim Galore` removes adapters and low-quality bases.
3. **Alignment**: `Bowtie2` maps reads with Read Group (RG) info.
4. **Filtering**: `Samtools` filters MAPQ >= 30, removes unmapped/secondary reads.
5. **Deduplication**: `Picard` or `Samtools` marks/removes PCR duplicates.
6. **Peak Calling**: `MACS3` with adaptive parameters for ChIP-seq vs CUT&Tag.
7. **Visualization**: `deepTools` generates CPM-normalized BigWig files.
8. **Aggregation**: `MultiQC` aggregates QC reports across all samples.

<br>

## Key Design Decisions

- **No-input by default**: The pipeline works without control/input samples. Set `use_control: true` and provide `control_bam` in the sample sheet for optional control support.
- **No ATAC-seq**: This pipeline is focused on ChIP-seq and CUT&Tag only.
- **CUT&Tag mode**: Uses `--tech cuttag` which applies high-resolution Tn5 insertion site centering (`--nomodel --shift -100 --extsize 200`).
- **Breakpoint resume**: Snakemake tracks completion via `.pipeline.done` sentinel files. Use `--rerun-incomplete` to resume failed runs without re-running completed samples.

<br>

## License

This project is open-source under the MIT License.
