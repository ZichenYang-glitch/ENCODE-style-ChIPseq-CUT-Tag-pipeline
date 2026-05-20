# No-Hardcoding Policy

The workflow should not depend on workstation-specific paths or author-local
environment layout. All analysis inputs and resources must enter through
configuration, sample sheets, Snakemake inputs, or Snakemake params.

## Rules

- FASTQ, BAM, Bowtie2 indexes, blacklist BED files, chrom sizes, GTF files,
  FASTA files, and external controls must come from `config.yaml`,
  `samples.tsv`, `input:`, or `params:`.
- Helper scripts should not infer the project output layout. They should accept
  explicit CLI arguments such as `--input`, `--output`, `--sample`, or
  `--threads`.
- Snakemake rules own path construction for workflow outputs.
- Tests may create temporary files, but runtime fixtures should live under
  temporary directories or generated test profiles, not under author-local
  paths.
- Documentation may show placeholder paths such as `/data/...`, but executable
  workflow logic must not depend on those paths.

## Enforcement

Run the guard before merging path-related changes:

```bash
python3 test/test_no_hardcoded_paths.py
```

The guard checks:

- no workstation-specific absolute paths in runtime files;
- no stale Snakemake `conda:` references to missing env files;
- no rule-level use of the legacy monolithic `chipseq.yml`;
- every workflow env file includes `nodefaults`.

Historical design notes under `docs/superpowers/` are not treated as runtime
contracts.
