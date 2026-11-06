---
name: helixweave-run
description: >-
  Prepare and launch an ENCODE-style ChIP-seq, CUT&Tag, ATAC-seq, or
  MNase-seq run from assay, FASTQ, reference, and sample-sheet inputs, or
  route a HelixWeave platform run through the validated-snapshot lifecycle.
  Use when the user asks to configure, validate, dry-run, start, resume, or
  check status of that workflow. Do not use for QC-only interpretation or
  when the main request is diagnosing an existing failure.
---

# HelixWeave run

Operate the checkout that contains this skill. Read versions from that tree.
This skill routes to the existing validator, Snakefile, and platform API. It
does not reimplement validation, default matrices, or identity hashing.

## Confirm the workflow

Use this skill for workflow id `encode-style-chipseq-cuttag-atac-mnase`
(`src/encode_pipeline/platform/registry.py`).

When the user named another adapter, stop. Read that adapter's registered
schema with `GET /api/v1/workflows/{workflow_id}/schema`, or read its adapter
package. Leave ENCODE sample columns, assays, output paths, and thresholds
unused.

QC-only interpretation and diagnosis of an existing failure belong to the
other HelixWeave skills. Stop when that is the main request.

## Read the command inventory

Read [../_shared/entries.md](../_shared/entries.md) before choosing a command.
That file is the command inventory. When it and the checkout disagree, follow
the checkout. `scripts/checkout_bootstrap.py` has no `snakemake` subcommand.
Launch Snakemake with the `snakemake` executable.

## Choose the path

Follow the path the user already chose. Ask once when they have not chosen.

- Direct Snakemake: sample sheet, config, the existing validator, a dry-run,
  then a real run after explicit authorization. Read
  [references/direct-workflow.md](references/direct-workflow.md).
- Platform: doctor, authoring schema, `POST` validate, create a run from the
  returned `snapshot_id` only, then an explicit start. Read
  [references/platform-run.md](references/platform-run.md).

## Check the ENCODE inputs

On the direct path, confirm the points below and leave the column manual in
the docs.

- `layout: PE` needs a non-empty `fastq_2` that is a different file from
  `fastq_1`. `layout: SE` needs an empty `fastq_2`. MNase-seq is PE-only.
- `control_sample` names a `role: control` row. MACS3 uses that row's
  canonical analysis BAM, not its raw `final.bam`. `control_bam` is an
  external path: `updated` rejects it; `legacy` still has to pass admission.
  Both fields apply only when `use_control` is true. Details are in
  [references/direct-workflow.md](references/direct-workflow.md).
- `bowtie2_index` is an index prefix.
- `genome` matches a key under `genome_resources`.

Read `docs/sample-sheet.md`, `docs/configuration.md` (`use_control` and
genome resources), and `docs/reference-resources.md`.

## Validate

Run the existing validator in the directory where the config's relative
paths resolve. The default command leaves FASTQ existence unchecked, so a
placeholder dry-run still validates. `ROOT` is this checkout:

```bash
python3 "$ROOT/scripts/validate_samples.py" --config CONFIG
python3 -I -S "$ROOT/scripts/checkout_bootstrap.py" \
  --repository-root "$ROOT" validate --config CONFIG
```

Pass `--strict-inputs` only when the user wants FASTQ and Bowtie2 index
existence checks before a real run. Both `-I` and `-S` are required.
`--repository-root` is this checkout even when the current directory is not.
A non-zero exit stops the path. Exit 2 whose text is
`source provenance check failed [reason_code]` means the config was not
validated. Do not edit the sample sheet or config for that line. Hand it to
[helixweave-troubleshoot](../helixweave-troubleshoot/SKILL.md). Exit 1 with
`ERROR:` is a config or sample failure; correct that and validate again
before the dry-run. The validator sees one config file. A later
`snakemake --config` override is checked when the workflow loads. See
[references/direct-workflow.md](references/direct-workflow.md).

## Outdir identity

Before any real execution, inspect the target outdir. Refuse it when
scientific outputs are present and the marker is absent, when the marker
mismatches, or when the effective-config record is broken. Choose a fresh
outdir. Leave the old directory in place: do not delete, move, or rewrite it,
and do not forge or repair a marker.

A dry-run performs the same check. It does not execute scientific rules
and does not create `.helixweave/analysis_identity.json`,
`.helixweave/effective_config.json`, or scientific outputs. Snakemake may
still create `.snakemake` scheduler metadata.
Resume only when the marker matches. The refusal table and the resume flag
are in [references/direct-workflow.md](references/direct-workflow.md).
Implementation: `workflow/Snakefile` (`_check_analysis_outdir_identity`,
`_print_effective_config_summary`). Human table: `docs/configuration.md`,
section "Run identity and effective configuration".

## Environment

A direct Snakemake run uses the runner lock `workflow/envs/runner.lock`
(Python, PyYAML, and Snakemake). Rule tools come from `workflow/envs/` when
the command passes `--use-conda`. Leave the platform `ci-fast` stack
uninstalled for that path.

A platform run uses `workflow/envs/ci-fast.lock` and the local platform
runtime. Read `docs/environments.md` and
`docs/development/local-platform-runtime.md`.

## Authorization

Validate and run the `-n` dry-run while preparing a direct run. Start, resume,
cancel, and download only after the user explicitly authorizes that action.
When the user asks how to run the workflow, stop after the dry-run plan.
On the platform path, stop before create and start unless they asked to
submit a run, and still leave `start` uncalled until they authorize execution.

## Evidence

Record the config path, the validated `analysis_profile` and
`config_schema_version`, the checkout `scientific_policy_revision`, the
outdir, and whether both `.helixweave` run records are present. Take the
profile, schema version, and policy revision from the workflow stderr summary
that begins `HelixWeave effective run configuration:`, or from the config
object that validation accepted. Resolve an omitted profile by running the
validator. The snapshot field `schema_version` is the authoring contract
version and is a different value from `config_schema_version`.

Read `scientific_policy_revision` from
`src/encode_pipeline/config/scientific_policy.py` or from that summary line.
Users cannot set it.

The skill text does not embed a historical commit or a personal absolute
tool path as a constant. A record of a run that actually happened keeps
the source SHA, whether the worktree was clean, and the interpreter and
tool paths and versions used by that command. Do not copy those paths
into a public API response. Platform responses omit workspace paths,
environment values, exception text, and private payloads.

For a platform run also record `workflow_id`, `snapshot_id`, `run_id`, and
any artifact or QC generation returned by the inventory's read routes.

## CUT&Tag duplicates

When `duplicate_analysis: auto` resolves to retain for CUT&Tag, say that this
retain behavior is a candidate pending the SCI-7 real-library data freeze.
Read the CUT&Tag duplicate section of `docs/assay-policy.md`. Leave the
per-assay matrix in that doc and in the workflow.

## Change boundary

Leave scientific defaults, admission, and production code unchanged. A full
workflow rerun and the Protected Bulk Gate are outside this skill.
