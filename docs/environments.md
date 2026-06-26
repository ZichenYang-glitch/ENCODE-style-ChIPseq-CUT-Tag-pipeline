# Environment Guide

This workflow uses a small runner environment plus Snakemake-managed
rule-specific tool environments. The goal is to make the first install fast and
keep heavy or conflicting tools isolated.

For the shortest setup path, create only the runner environment first:

```bash
micromamba create -n chipseq-runner --file workflow/envs/runner.lock
micromamba activate chipseq-runner
```

Then run Snakemake with `--use-conda`:

```bash
snakemake -s workflow/Snakefile --configfile config/config.yaml -n --use-conda
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 16 --use-conda
```

Snakemake creates the rule environments it needs on first use and reuses them
on later runs.

## Locked Environments

Every `workflow/envs/*.yml` has a matching `workflow/envs/*.lock` file generated
with [conda-lock](https://github.com/conda/conda-lock) in **explicit** format.
The YAML files keep semantic version ranges as the intent document; the
lockfiles capture the exact resolved package URLs so CI and users can reproduce
the same environment with any conda-compatible installer.

We evaluated `conda-lock` against `pixi` on the five most constrained
environments (`chipseq`, `idr`, `macs3`, `deeptools`, `picard`) on `linux-64`.
`conda-lock` resolved all of them, including `idr.yml` with its `python <3.10`
constraint. We render the lockfiles in explicit format so they can be installed
with `micromamba create --file` directly, without requiring `conda-lock` at
install time. For this reason the project uses `conda-lock`.

To install from a lockfile:

```bash
micromamba create -n chipseq --file workflow/envs/chipseq.lock
```

Or with conda:

```bash
conda create -n chipseq --file workflow/envs/chipseq.lock
```

## Updating an Environment

1. Edit the source YAML (`workflow/envs/<name>.yml`) to express the new intent.
2. Regenerate the lockfile:

   ```bash
   conda-lock lock -f workflow/envs/<name>.yml -p linux-64 --kind explicit \
       --lockfile workflow/envs/<name>.lock
   ```

3. Commit both the `.yml` and the `.lock`.

The CI `lock-check` workflow fails if a YAML is modified without a matching
`.lock` file.

| File | Used for |
|------|----------|
| `workflow/envs/runner.yml` | First install and local workflow entry point; lockfile is `runner.lock` |
| `workflow/envs/ci-fast.yml` | GitHub Actions fast checks; lockfile is `ci-fast.lock` |
| `workflow/envs/chipseq.yml` | Core runtime for tiny real-execution harness; lockfile is `chipseq.lock` |
| `workflow/envs/fastqc.yml` | FastQC with pinned OpenJDK 17 |
| `workflow/envs/trim.yml` | Trim Galore and cutadapt |
| `workflow/envs/align.yml` | Bowtie2 alignment and samtools conversion |
| `workflow/envs/samtools.yml` | samtools and bedtools operations |
| `workflow/envs/macs3.yml` | MACS3 peak calling and bedGraph signal rules |
| `workflow/envs/deeptools.yml` | deepTools `bamCoverage`, `computeMatrix`, and `plotProfile` |
| `workflow/envs/multiqc.yml` | MultiQC report aggregation |
| `workflow/envs/idr.yml` | IDR, isolated because Bioconda IDR currently requires `python <3.10` |
| `workflow/envs/seacr.yml` | SEACR, isolated because it brings an R runtime |
| `workflow/envs/python.yml` | Python-only helper rules |

All environment files use `conda-forge`, `bioconda`, and `nodefaults`.

## Why Not One Big Environment?

A single all-in-one environment is convenient in small projects, but it became
fragile here:

- IDR has older Python constraints in Bioconda.
- SEACR brings R.
- MultiQC brings a broad reporting stack.
- FastQC brings Java.
- deepTools and MACS3 bring modern Python scientific dependencies.

Putting all of those into one solve made installation slow and brittle. The
split layout lets each rule use the smallest compatible environment.

## What Happens on First Run?

The first `--use-conda` run may still spend time creating tool environments.
That is expected. The difference is that the first install no longer has to
solve every tool at once.

Typical behavior:

- First setup: install `chipseq-runner`.
- First dry-run or run with `--use-conda`: Snakemake prepares required rule
  environments.
- Later runs: Snakemake reuses the cached environments unless the env YAML
  changes.

By default, Snakemake stores rule environments under `.snakemake/conda/` in the
work directory. The directory is generated output and should not be committed.

## Common Commands

Create the runner:

```bash
micromamba create -n chipseq-runner --file workflow/envs/runner.lock
micromamba activate chipseq-runner
```

Dry-run with rule environments:

```bash
snakemake -s workflow/Snakefile --configfile config/config.yaml -n --use-conda
```

Run with rule environments:

```bash
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores 16 --use-conda
```

Use a custom Snakemake environment cache location:

```bash
snakemake -s workflow/Snakefile --configfile config/config.yaml \
  --cores 16 --use-conda --conda-prefix /path/to/snakemake-conda-envs
```

Remove the local runner environment:

```bash
micromamba env remove -n chipseq-runner
```

Remove project-local Snakemake rule environments:

```bash
rm -rf .snakemake/conda
```

## Tiny Real Execution

The Stage 8b tiny real-execution harness intentionally uses the core
`chipseq` runtime instead of creating all rule-specific environments:

```bash
micromamba create -n chipseq --file workflow/envs/chipseq.lock
micromamba activate chipseq
python3 test/test_stage8b_tiny_execution.py
```

This test covers real preprocessing and signal generation on synthetic data.
It does not run IDR, SEACR, or MultiQC.

## CI Behavior

GitHub Actions uses two paths:

- `fast-checks`: uses `workflow/envs/ci-fast.lock` for config validation and
  dry-run smoke profiles.
- `real-execution`: manual `workflow_dispatch`, uses `workflow/envs/chipseq.lock`
  for the tiny real-execution harness.

Full rule-specific environments are primarily for real user runs with
`--use-conda`.

## Troubleshooting

If `snakemake` is not found, activate the runner:

```bash
micromamba activate chipseq-runner
```

If a rule environment fails after an env YAML change, remove the Snakemake
environment cache and retry:

```bash
rm -rf .snakemake/conda
snakemake -s workflow/Snakefile --configfile config/config.yaml -n --use-conda
```

If Conda channel mixing appears in logs, prefer strict channel priority:

```bash
conda config --set channel_priority strict
```

If downloads are slow, rerun the same command. Conda and micromamba reuse
partial package caches where possible, so a second attempt often resumes much
faster than the first.
