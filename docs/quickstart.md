# Quick Start

This guide walks through installing the pipeline, configuring your data, and
running your first analysis. For a compact overview, see the
[README](../README.md).

## 1. Install

```bash
git clone https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline.git
cd ENCODE-style-ChIPseq-CUT-Tag-pipeline
micromamba create -f workflow/envs/runner.yml
micromamba activate chipseq-runner
```

The `chipseq-runner` environment is intentionally small: Python, PyYAML, and
Snakemake. The bioinformatics tools live in rule-specific environments and are
created automatically when you run Snakemake with `--use-conda`.
For details about every environment file and cache cleanup, see
[docs/environments.md](environments.md).

## 2. Configure samples

Edit `config/samples.tsv` with your FASTQ paths and metadata. The minimum
required columns are: `sample`, `fastq_1`, `fastq_2` (empty for SE), `layout`,
`assay`, `target`, `peak_mode`, `genome`, `bowtie2_index`.

See [docs/sample-sheet.md](sample-sheet.md) for the full column reference,
optional replicate/control columns, role semantics, worked examples, and
common pitfalls.

## 3. Adjust workflow options

Edit `config/config.yaml`. The defaults are reasonable for most ChIP-seq runs:

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

See [docs/configuration.md](configuration.md) for the full configuration
reference covering genome resources, QC switches, replicate/IDR features,
CUT&Tag SEACR, tool parameters, and key dependencies.
See [docs/reference-resources.md](reference-resources.md) for preparing
Bowtie2 indexes, FASTA `.fai`/`.dict`, chrom sizes, blacklists, and annotations.

## 4. Validate

Always validate your config and sample sheet before running:

```bash
python3 scripts/validate_samples.py --config config/config.yaml
```

A non-zero exit code means something is wrong — fix the reported issues
before proceeding.

Then check the DAG with a dry-run:

```bash
snakemake -s workflow/Snakefile --configfile config/config.yaml -n --use-conda
```

## 5. Dry-run

The dry-run resolves the full DAG and reports which rules would run:

```bash
snakemake -s workflow/Snakefile --configfile config/config.yaml -n --use-conda
```

This checks rule connectivity, file path resolution, and sample sheet
parsing without executing any tool. `--use-conda` tells Snakemake to use the
rule-specific tool environments declared under `workflow/envs/`.

## 6. Run the workflow

```bash
# Execute with 16 cores
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 16 --use-conda

# Resume after interruption
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 16 --use-conda --rerun-incomplete
```

- **`--cores N`:** use N local CPU cores. Replace `N` with your available core count.
- **`--jobs`:** use with `--cluster` for distributed execution (cluster/cloud).
- **`--rerun-incomplete`:** re-run jobs that produced incomplete output files.
- **`--latency-wait N`:** wait N seconds for NFS filesystem sync; useful on shared
  storage.

## Smoke and test execution

The pipeline ships with three bundled test harnesses that run on synthetic data
(no real FASTQs needed):

```bash
# Validation stress tests — 15 config/sample checks, <5 seconds
python3 test/test_validation_stress.py

# Dry-run smoke profiles — 7 dispatch paths, <30 seconds
# Requires snakemake on PATH.
python3 test/test_stage8_smoke_profiles.py

# Tiny real execution — preprocessing + signal on synthetic 1k reads, <60 seconds
# Requires the core chipseq Conda environment.
python3 test/test_stage8b_tiny_execution.py
```

All temporary output lands under `/tmp`. Exit code 0 = PASS.

The stress tests and smoke profiles are fast and check config parsing and DAG
connectivity. The tiny execution test runs real tools (Bowtie2, samtools,
deepTools) on synthetic fixtures and confirms the preprocessing+signal path
end to end. MACS3 is intentionally skipped in tiny execution (dry-run smoke
profiles cover the peak-calling DAG).

## Troubleshooting

### Slow conda/micromamba solve

Prefer micromamba for environment creation. If plain `conda env create` hangs
on "Solving environment," try:

```bash
# Use the libmamba solver (faster)
conda install -n base conda-libmamba-solver
conda env create -f workflow/envs/runner.yml --experimental-solver=libmamba
```

The environment files use `nodefaults` so local `defaults` channels are not
silently mixed into solves. Heavy or optional tools are isolated in separate
rule environments so the first install stays small.

### `snakemake: command not found`

- Make sure you activated the runner environment: `micromamba activate chipseq-runner`
- Verify snakemake is installed: `which snakemake`
- If the environment was created but snakemake is missing, remove and recreate
  it with `micromamba env remove -n chipseq-runner` followed by
  `micromamba create -f workflow/envs/runner.yml`.

### FASTQ path not found

- Check that paths in `config/samples.tsv` exist on disk.
- Use absolute paths if relative paths don't resolve from the repository root.
- On WSL2, Windows paths need translation: `/mnt/c/Users/...` instead of
  `C:\Users\...`.

### Bowtie2 index prefix issues

- The `bowtie2_index` value is a **prefix**, not a directory. For an index built
  at `/data/genomes/hg38/GRCh38`, the value should be `/data/genomes/hg38/GRCh38`
  (Bowtie2 appends `.1.bt2`, `.2.bt2`, etc.).
- Verify the index files exist: `ls /data/genomes/hg38/GRCh38*.bt2`
- If only a FASTA is available, build the index first:
  `bowtie2-build /data/genomes/hg38.fa /data/genomes/hg38/GRCh38`
- For a complete reference resource checklist, see
  [docs/reference-resources.md](reference-resources.md).

### Control path not found

- `control_bam` points to an external BAM that must exist on disk.
- `control_sample` must match a `sample` ID in the same sample sheet.
- Controls are disabled unless `use_control: true` in `config/config.yaml`.

### Missing tools for tiny execution

The tiny execution test (`test_stage8b_tiny_execution.py`) requires:
- `bowtie2` and `bowtie2-build` on PATH
- `samtools` on PATH
- `deepTools` (`bamCoverage`) on PATH

These are all included in the core `chipseq` Conda environment. Activate it first:
`micromamba activate chipseq` or `conda activate chipseq`.
