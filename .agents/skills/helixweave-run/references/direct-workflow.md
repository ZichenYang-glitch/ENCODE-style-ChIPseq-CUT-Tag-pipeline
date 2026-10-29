# Direct Snakemake run

Read [../../_shared/entries.md](../../_shared/entries.md) before choosing a
command. Maintained path: `docs/quickstart.md`. This note is the operational
route for workflow id `encode-style-chipseq-cuttag-atac-mnase`. It does not
restate the sample-sheet or configuration manuals.

Commands below use `CONFIG` for the config file the user is editing. Run them
in the directory where that file's relative paths resolve. `samples` and
other relative paths are resolved from the current directory, not from the
config file's directory. `--repository-root` stays the checkout that
contains the bootstrap, even when the current directory is a profile folder.

## Sample sheet

Edit the sample sheet named by the config `samples` key. The column contract
is `docs/sample-sheet.md`.

Paired-end and single-end checks:

- `layout: PE` requires non-empty `fastq_2`. `fastq_1` and `fastq_2` must be
  different files after path normalization.
- `layout: SE` requires empty `fastq_2`. A filled R2 is rejected.
- `assay: mnase` is PE-only and uses `peak_mode: nucleosome`. Other assays
  reject `nucleosome`. ATAC-seq in this release is narrow-peak only.
- `control_sample` is another row's `sample` id with `role: control`. The
  pipeline uses that row's `final.bam` as the MACS3 control.
- `control_bam` is an external BAM path. The pipeline uses the file directly
  and does not process it as a sample row.
- `use_control: false` ignores both control fields. Control rows are not
  scheduled as active outputs. See `docs/configuration.md`, section
  `use_control`.
- `bowtie2_index` is a prefix. Bowtie2 appends the `.bt2` or `.bt2l`
  suffixes. See `docs/reference-resources.md`, section "Bowtie2 Index Prefix".
- `genome` must equal a key under `genome_resources`. See
  `docs/reference-resources.md` and the genome-resources section of
  `docs/configuration.md`.

## Validator

Primary command, from the inventory and `docs/quickstart.md`. `ROOT` is the
checkout. The current directory is the one in which `CONFIG`'s relative
paths exist:

```bash
python3 "$ROOT/scripts/validate_samples.py" --config CONFIG
python3 "$ROOT/scripts/validate_samples.py" --config CONFIG --strict-inputs
```

Provenance-checked equivalent. Remaining arguments forward to
`scripts/validate_samples.py`:

```bash
python3 -I -S "$ROOT/scripts/checkout_bootstrap.py" \
  --repository-root "$ROOT" validate --config CONFIG
```

Flags implemented by `src/encode_pipeline/cli/_validator.py`:

| Flag | Required | Behavior |
| --- | --- | --- |
| `--config CONFIG` | yes | One config YAML |
| `--strict-inputs` | no | Also require FASTQ files and a complete Bowtie2 index at the prefix |

Default validation leaves FASTQ and Bowtie2 index existence unchecked so
placeholder dry-runs still pass. `docs/configuration.md`, section "Strict
input validation", is the mode table. Genome-resource paths, and `control_bam`
when `use_control` is true, are checked without `--strict-inputs`. Use
`--strict-inputs` only when the user wants real FASTQ and index checks before
a real run.

The validator prints a sample count, `use_control`, and, when set, the
genome-resource keys. It does not print `analysis_profile`. Read that from
the workflow summary
below. An omitted profile is resolved only by
`src/encode_pipeline/config/validator.py` (`_validate_analysis_profile`)
together with `workflow/schemas/config.schema.yaml`. Run the validator and
the dry-run. Do not choose `legacy` or `updated` from a missing key.

Exit 2 with `source provenance check failed [reason_code]` is not a sample
or config error. The validator did not run. Leave the files unchanged and
follow the provenance line. Exit 1 with `ERROR:` is the config failure that
blocks the dry-run.

## Validator and `snakemake --config`

`scripts/validate_samples.py` validates the config file only. It does not see
later `snakemake --config KEY=VALUE` overrides. Snakemake merges those
overrides with `--configfile` before workflow validation, so the merged value
goes through the same validator, the run identity, and the effective-config
record. To check an override, dry-run the workflow with the same `--config`.
See `docs/configuration.md`, section "CLI overrides and the effective
configuration record".

## Dry-run and real run

Dry-run, from the same directory:

```bash
snakemake -s "$ROOT/workflow/Snakefile" --configfile CONFIG -n --use-conda
```

The dry-run loads the Snakefile, prints the effective summary to stderr, and
runs `_check_analysis_outdir_identity`. The marker rule
`analysis_identity_marker` in `workflow/rules/analysis.smk` is a real-run
output. A dry-run does not execute it and writes neither run record.

Real run, only after the user authorizes execution and the outdir check
below allows it:

```bash
snakemake -s "$ROOT/workflow/Snakefile" --configfile CONFIG --cores N --use-conda
```

`N` is the core count the user chose. `--use-conda` selects rule environments
under `workflow/envs/`. The runner environment is `workflow/envs/runner.lock`.
See `docs/environments.md`.

On startup the workflow prints, to stderr, unless the invocation is
`snakemake --lint`:

```text
HelixWeave effective run configuration:
  analysis_profile=<profile> config_schema_version=<version> scientific_policy_revision=<revision>
```

Record those three values from that line. When the summary shows
`duplicate_analysis=auto(...)` and CUT&Tag resolves to `retain`, say that
CUT&Tag auto-retain is a candidate pending the SCI-7 freeze. Read
`docs/assay-policy.md`.

## Resume

Resume an interrupted run only when the marker matches this profile, this
`scientific_policy_revision`, and this identity, and `effective_config.json`
recomputes to that identity. The user must authorize the resume.

```bash
snakemake -s "$ROOT/workflow/Snakefile" --configfile CONFIG --cores N --use-conda \
  --rerun-incomplete
```

`--rerun-incomplete` reruns jobs that left incomplete outputs. It does not
make a mismatched outdir reusable.

## `profiles/default` and `analysis_profile`

`profiles/default` is a Snakemake scheduler profile. Its `config.yaml` sets
`configfile`, `cores`, `use-conda`, `conda-frontend`, `rerun-incomplete`,
`keep-going`, and `show-failed-logs`.
`snakemake --profile profiles/default` starts a real run with those settings.
Read `profiles/default/config.yaml` before using it, and use it only when the
user authorizes that real run.

`analysis_profile` is the scientific profile, `legacy` or `updated`. It is a
config key, validated by the existing validator, and stored in the run marker.
The scheduler profile does not select it.

## Outdir identity

Records, under the configured outdir:

- `.helixweave/analysis_identity.json`
- `.helixweave/effective_config.json`

The marker stores `marker_schema_version`, `analysis_profile`,
`scientific_policy_revision`, and an identity hash. The workflow computes the
hash. The hash covers the validated config plus parsed samples, and excludes
`outdir`, `threads`, and the `samples` path. Path-valued inputs use canonical
absolute locations. `scientific_policy_revision` comes from
`src/encode_pipeline/config/scientific_policy.py`. Users cannot set it.

Inspect the outdir before a real run. A refusal means: use a fresh outdir.
Leave the existing directory unchanged. Do not delete, move, or rewrite it,
and do not reconstruct a record from the current config.

| What you find | What the workflow does |
| --- | --- |
| New directory, or only `.snakemake` and `.helixweave` | A real run creates the marker, then scientific tasks run |
| Marker matches profile, policy revision, and identity, and the effective-config record recomputes to that identity | Resume is allowed |
| Scientific outputs exist and there is no marker | Refused. Unmarked outdir. Use a fresh outdir |
| Outdated marker format | Refused. Use a fresh outdir |
| Marker owned by another `analysis_profile` | Refused. Use a fresh outdir for this profile |
| Marker written by another `scientific_policy_revision` | Refused. Use a fresh outdir for this workflow version |
| Marker identity differs | Refused. Use a fresh outdir or the profile that created it |
| `.helixweave` or a run record is a symlink, or is not a directory or regular file | Refused. Leave the metadata unrepaired |
| Marker is valid and `effective_config.json` is missing, unreadable, or inconsistent with the marker | Refused. Use a fresh outdir. Leave the record unreconstructed |

Inconsistent effective-config records include a missing required field, an
identity that differs from the marker, another profile, a schema version that
disagrees with the recorded config, another policy revision, or a config
block whose recomputed identity differs from the recorded identity.

A dry-run reports the same refusals and writes neither run record. Messages come from
`_check_analysis_outdir_identity`, `_read_analysis_identity_marker`, and
`_check_effective_config_record` in `workflow/Snakefile`. The human table is
`docs/configuration.md`, section "Run identity and effective configuration".

Changing `outdir` points the same scientific identity at a new directory,
because `outdir` is outside the hash. A scientific override that changes the
identity still needs its own fresh outdir. Check that override with a dry-run.
The standalone validator will not see it.
