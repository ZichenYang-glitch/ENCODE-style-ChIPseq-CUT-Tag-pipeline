# Sample Sheet Reference

The sample sheet (`config/samples.tsv`) defines which samples to process and how.
It is a tab-separated file with one row per sample. This page is the maintained
reference for its required columns, optional metadata, and examples; the
[README](../README.md) provides the product-level entry point.

## Required columns

| Column | Type | Description |
| :--- | :--- | :--- |
| `sample` | string | Unique sample ID. Must match `[A-Za-z0-9_.-]+`. |
| `fastq_1` | path | R1 FASTQ file. |
| `fastq_2` | path | R2 FASTQ file. Required for PE; leave empty for SE. |
| `layout` | string | `PE` or `SE`. |
| `assay` | string | `chipseq`, `cuttag`, `atac`, or `mnase`. |
| `target` | string | Antibody or target name (e.g. `H3K27ac`, `CTCF`). |
| `peak_mode` | string | `narrow`, `broad`, or `nucleosome`. `nucleosome` is MNase-only; other assays reject it. |
| `genome` | string | Genome label matching a key in `genome_resources` (e.g. `hs`, `mm`, `hg38`). |
| `bowtie2_index` | path | Bowtie2 index prefix (not a directory). |

## Optional columns

### Replicate metadata

| Column | Default | Description |
| :--- | :--- | :--- |
| `experiment` | `<sample>` | Experiment group identifier. Samples sharing an `experiment` are treated as replicates of the same condition. |
| `condition` | `<target>` | Condition label (treatment, genotype, timepoint). |
| `replicate` | `1` | Replicate number within an experiment. |
| `biological_replicate` | `<replicate>` | Biological replicate number. |
| `technical_replicate` | `1` | Technical replicate number within a biological replicate. |

Technical replicates within the same biological replicate are merged into a
single `biorep<N>.final.bam`. A single technical replicate is symlinked (no
duplicate data). Experiments with 2+ biological replicates produce pooled BAMs
and pooled peaks when `stage4b: true` (the default).

All replicate columns are optional. A minimal single-sample sheet with only the
required columns is fully supported.

### Control columns

| Column | Default | Description |
| :--- | :--- | :--- |
| `role` | `treatment` | `treatment` or `control`. Controls skip peak calling. |
| `control_sample` | (none) | Sample ID of a control row to use as MACS3 input (`-c`). |
| `control_bam` | (none) | Path to an external control BAM file. |

`control_sample` and `control_bam` are ignored when `use_control: false` in
`config/config.yaml`. Rows with `role: control` are not scheduled as active
outputs unless controls are enabled.

## Role semantics

- **`role: treatment`** (default): The sample flows through the full
  preprocessing pipeline and MACS3 peak calling runs on it.
- **`role: control`**: The sample is preprocessed (FastQC → align → filter →
  `final.bam`) but peak calling is skipped. Its `final.bam` can be referenced
  by treatment rows via `control_sample`.
- **`control_sample`**: references another row's `sample` ID. That row must have
  `role: control`. The referenced sample's `final.bam` becomes the MACS3 input
  control (`-c`).
- **`control_bam`**: an external BAM path, not processed by the pipeline. Used
  directly as MACS3 input control. The file must exist on disk.

## Examples

### Minimal single-sample

A single SE ChIP-seq sample, no controls, no replicates:

```tsv
sample	fastq_1	fastq_2	layout	assay	target	peak_mode	genome	bowtie2_index
H3K27ac_WT	/data/wt_R1.fq.gz		SE	chipseq	H3K27ac	narrow	hg38	/data/genomes/hg38/GRCh38
```

### Treatment with control sample

Two rows: a treatment sample and an input/control sample from the same
experiment:

```tsv
sample	fastq_1	fastq_2	layout	assay	target	peak_mode	genome	bowtie2_index	role	control_sample
H3K27ac_WT	/data/wt_R1.fq.gz	/data/wt_R2.fq.gz	PE	chipseq	H3K27ac	narrow	hg38	/data/genomes/hg38/GRCh38	treatment	INPUT_WT
INPUT_WT	/data/input_R1.fq.gz	/data/input_R2.fq.gz	PE	chipseq	Input	narrow	hg38	/data/genomes/hg38/GRCh38	control
```

Set `use_control: true` in `config/config.yaml` for the control reference
to take effect. The treatment row sets `control_sample: INPUT_WT` to use
the control row's `final.bam` as MACS3 input.

### Biological replicates

Two biological replicates of the same target, pooled for replicate-aware
analysis:

```tsv
sample	fastq_1	fastq_2	layout	assay	target	peak_mode	genome	bowtie2_index	experiment	biological_replicate
H3K27AC_rep1	/data/ac1_R1.fq.gz	/data/ac1_R2.fq.gz	PE	chipseq	H3K27ac	narrow	hs	/path/to/bt2/GRCh38	H3K27AC	1
H3K27AC_rep2	/data/ac2_R1.fq.gz	/data/ac2_R2.fq.gz	PE	chipseq	H3K27ac	narrow	hs	/path/to/bt2/GRCh38	H3K27AC	2
```

With `stage4b: true` (default), the pipeline produces pooled BAMs and pooled
peaks for the `H3K27AC` experiment. With `stage5: true`, TF ChIP-seq IDR runs
on the two biorep peak sets.

### Baseline ATAC-seq

A paired-end ATAC-seq sample uses the same sample-sheet contract with
`assay: atac` and `peak_mode: narrow`:

```tsv
sample	fastq_1	fastq_2	layout	assay	target	peak_mode	genome	bowtie2_index
ATAC_rep1	/data/atac_R1.fq.gz	/data/atac_R2.fq.gz	PE	atac	ATAC	narrow	hg38	/data/genomes/hg38/GRCh38
```

ATAC currently uses MACS3 narrow-peak calling with a Tn5-aware shift/extension
policy. ATAC-specific footprinting and nucleosome-positioning modules are not
included.

### MNase-seq

A paired-end MNase-seq sample for nucleosome positioning analysis:

```tsv
sample	fastq_1	fastq_2	layout	assay	target	peak_mode	genome	bowtie2_index
MNase_WT	/data/mnase_R1.fq.gz	/data/mnase_R2.fq.gz	PE	mnase	H3	nucleosome	hs	/path/to/bt2/GRCh38
```

MNase-seq is **PE-only** — SE samples are rejected during validation. MNase
requires `peak_mode: nucleosome`. MNase samples skip MACS3 peak calling and
instead produce sub/mono/di-nucleosome BAMs, dyad BigWig, mono occupancy
BigWig, and per-sample MNase QC summary outputs. Nucleosome calling
(DANPOS3/iNPS/SEM) is not implemented.

## Common pitfalls

- **PE sample missing `fastq_2`**: The pipeline expects `fastq_2` to be non-empty
  for `layout: PE`. An empty `fastq_2` with `layout: PE` causes Bowtie2 to fail.
- **Sample ID contains spaces or special characters**: Sample IDs must match
  `[A-Za-z0-9_.-]+`. Spaces, slashes, and other characters cause Snakemake
  wildcard parsing errors.
- **`control_sample` references a treatment row**: `control_sample` must point to
  a row with `role: control`. Referencing a treatment row creates a circular
  dependency or silently produces incorrect results.
- **`control_sample` references itself or forms a cycle**: A control row cannot
  reference itself or another row that references back to it.
- **`control_bam` path set while `use_control` is `false`**: Controls are
  disabled unless `use_control: true`; the external BAM will not be used.
- **`genome` key not in `genome_resources`**: The `genome` column value must
  match a key in `config/config.yaml` under `genome_resources`. A mismatch
  causes MACS3 effective genome size lookup to fail.
- **ATAC configured with `peak_mode: broad`**: ATAC support is intentionally
  narrow-peak only in this release; broad ATAC rows are rejected during
  validation.
