# ENCODE-style ChIP-seq and CUT&Tag Pipeline

This repository contains a Snakemake-based pipeline suite for ChIP-seq and
CUT&Tag data analysis. The batch workflow is the primary entry point and is
implemented as modular rules under `workflow/rules/`.

Default mode is no-input / no-control. Optional controls can be enabled with
`use_control: true` and supplied either as an external control BAM or as a
FASTQ-based control sample row.

This pipeline does not include ATAC-seq analysis.

## Repository Structure

- `workflow/Snakefile`: Main Snakemake entry point, sample validation, and assay dispatch.
- `workflow/rules/common.smk`: Shared rules for QC, trimming, alignment, filtering, duplicate handling, QC stats, and BigWig generation.
- `workflow/rules/chipseq.smk`: ChIP-seq policy functions.
- `workflow/rules/cuttag.smk`: CUT&Tag policy functions.
- `workflow/rules/peaks.smk`: MACS3 peak calling rule.
- `workflow/rules/report.smk`: Completion sentinels and MultiQC aggregation.
- `workflow/envs/chipseq.yml`: Conda environment for the Snakemake workflow.
- `config/config.yaml`: Workflow configuration.
- `config/samples.tsv`: Sample sheet.
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

### 2. Configure samples

Edit `config/samples.tsv`:

```tsv
sample	fastq_1	fastq_2	layout	assay	target	peak_mode	genome	bowtie2_index	control_bam
H3K27AC	/data/ac_R1.fq.gz	/data/ac_R2.fq.gz	PE	chipseq	H3K27ac	narrow	mm	/path/to/bt2/mm10
CTCF_rep1	/data/ctcf_R1.fq.gz		SE	chipseq	CTCF	narrow	hs	/path/to/bt2/hg38
CUTTAG_H3K4me3	/data/ct_R1.fq.gz	/data/ct_R2.fq.gz	PE	cuttag	H3K4me3	narrow	mm	/path/to/bt2/mm10
```

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

### 3. Configure workflow options

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

### 4. Run

```bash
# Dry-run: build the DAG without executing commands.
snakemake -s workflow/Snakefile --configfile config/config.yaml -n --use-conda

# Run with 16 cores. --use-conda creates/uses workflow/envs/chipseq.yml.
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 16 --use-conda

# Resume incomplete jobs after an interrupted run.
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 16 --use-conda --rerun-incomplete
```

## Control/Input Handling

Controls are disabled by default. With `use_control: false`, both `control_bam`
and `control_sample` are ignored.

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

## Output Structure

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
│   │   └── <sample>.final.bam
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

This project is open-source under the MIT License.
