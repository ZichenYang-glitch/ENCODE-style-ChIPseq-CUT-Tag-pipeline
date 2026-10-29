# Shared command and contract inventory

Read this before choosing a command. It is the single command inventory for
`helixweave-run`, `helixweave-qc`, and `helixweave-troubleshoot`. Do not
copy it into a skill. If this file and the checkout disagree, follow the
checkout and update the skill only when the user asked to maintain skills.

These skills operate the checkout that contains them. Read versions from
that tree. The skill text does not embed one historical commit, one
personal directory, one past result path, or one absolute tool path as a
universal constant.

## What is universal and what is ENCODE-only

Platform lifecycle, validated snapshots, artifact and QC generation cursors,
and error envelopes are workflow-neutral. Sample columns, assays, output
paths, duplicate and organelle defaults, and QC tables below apply only to
workflow id `encode-style-chipseq-cuttag-atac-mnase`.

For any other id, read that adapter's registered schema and its own docs.
`helixweave --doctor` checks the default registry shape in
`src/encode_pipeline/cli/local_platform.py` (`EXPECTED_WORKFLOW_IDS`).
Current members include `bulk-rnaseq` and `hitrac-preprocess`. Do not apply
ENCODE parameters, paths, or thresholds to them.

Do not change scientific defaults, admission, or production code from a
skill. Do not delete result directories, forge identity markers, or relax
validation. A full workflow rerun and the Protected Bulk Gate are not
default checks. Do not download sequencing data. Start, resume, cancel, or
download only after the user explicitly authorizes that action. Explaining
QC is read-only.

## Surfaces

### Source-checkout bootstrap

```bash
python3 -I -S scripts/checkout_bootstrap.py --repository-root . COMMAND
```

Implemented commands are only:

- `verify-checkout`
- `installed-artifact`
- `pytest`
- `openapi`
- `validate`
- `local-platform`

Both `-I` and `-S` are required. `--repository-root` must be the checkout
that contains that bootstrap. There is no `snakemake` subcommand. `validate`
forwards to `scripts/validate_samples.py`. `local-platform` forwards to
`scripts/run_local_platform.py`. Provenance failures are printed by
`scripts/source_provenance.py` and `scripts/checkout_bootstrap.py`; see
`docs/development/source-provenance.md`.

### Direct ENCODE-style Snakemake workflow

Maintained path: `docs/quickstart.md`.

```bash
python3 scripts/validate_samples.py --config CONFIG
python3 scripts/validate_samples.py --config CONFIG --strict-inputs
snakemake -s workflow/Snakefile --workflow-profile none \
  --configfile CONFIG -n --use-conda
snakemake -s workflow/Snakefile --workflow-profile none \
  --configfile CONFIG --cores N --use-conda
```

Resume an interrupted matching run with `--rerun-incomplete`. Snakemake
8.30 loads `profiles/default` from the current directory, then from the
Snakefile directory, unless `--workflow-profile none` is given. This
checkout's `profiles/default` sets another `configfile`, `use-conda`,
`rerun-incomplete`, and `keep-going`. Pass `--workflow-profile none` when
the user chose a config file. Use `snakemake --profile profiles/default`
only when they asked for that scheduler profile. That profile is not
`analysis_profile`.

`scripts/validate_samples.py` is a provenance-checked wrapper around
`encode_pipeline.cli.validate`. It validates one config file. It does not see
later `snakemake --config KEY=VALUE` overrides. Those overrides are merged
and validated when the workflow loads. A dry-run performs the identity check
and writes neither `analysis_identity.json` nor `effective_config.json`.

The small runner environment is `workflow/envs/runner.lock` (Python, PyYAML,
Snakemake). Rule tools come from `workflow/envs/` via `--use-conda`. That
environment is not the platform `workflow/envs/ci-fast.lock` stack.

### Installed product CLI

`helixweave` and `python -m encode_pipeline` dispatch in
`src/encode_pipeline/cli/app.py`:

- `admin` -> `encode_pipeline.cli.admin`
- `bundle` -> deployment bundle producer
- `install`, `status`, `doctor`, `verify`, `upgrade`, `rollback` ->
  `encode_pipeline.deployment.cli`
- anything else -> local platform supervisor (`--doctor` checks without
  starting services)

Compatibility console scripts: `encode-validate`, `encode-manifest`,
`encode-dag`, `encode-worker`.

From a verified source checkout, the maintained doctor command is:

```bash
python3 -I -S scripts/checkout_bootstrap.py \
  --repository-root . local-platform --doctor
```

`helixweave --doctor` is the installed equivalent. Doctor requirements are
in `docs/development/local-platform-runtime.md`: Python 3.12, Snakemake
8.30.0, Redis 7 or newer, Node from `.nvmrc`, and the locked frontend.
Doctor exits before opening ports. A missing Bulk RNA-seq binding is
`not_configured` or `unavailable` and does not fail the core doctor.

### Platform HTTP lifecycle

Base path `/api/v1`. Member routes require a principal and CSRF. Create a
run only from a server-owned validated snapshot:

- `GET /api/v1/workflows`
- `GET /api/v1/workflows/{workflow_id}`
- `GET /api/v1/workflows/{workflow_id}/schema`
- `GET /api/v1/workflows/{workflow_id}/reference-profiles`
- `POST /api/v1/workflows/{workflow_id}/validate`
- `POST /api/v1/workflows/{workflow_id}/runs` with `snapshot_id`
- `POST /api/v1/runs/{run_id}/preflight`
- `POST /api/v1/runs/{run_id}/start`
- `POST /api/v1/runs/{run_id}/cancel`
- `GET /api/v1/runs/{run_id}`
- `GET /api/v1/runs/{run_id}/events`
- `GET /api/v1/runs/{run_id}/logs`
- `GET /api/v1/runs/{run_id}/artifacts`
- `GET /api/v1/runs/{run_id}/artifacts/{artifact_id}`
- `GET /api/v1/runs/{run_id}/artifacts/{artifact_id}/download`
- `GET /api/v1/runs/{run_id}/qc-metrics`

SQLite is lifecycle truth. Redis/RQ is only the queue. Do not enqueue or
patch a run through a private helper. Status values live in
`encode_pipeline.platform.runs.RunStatus`: `created`, `validating`,
`planned`, `queued`, `running`, `succeeded`, `failed`, `cancelled`.
`succeeded` means execution completed. It is not biological acceptance.
Create leaves a run `created`. Start accepts `planned`. The public
preflight route is what accepts `created` and moves the run to
`validating` (`src/encode_pipeline/api/routes/preflight.py`). A second
preflight returns `PREFLIGHT_ALREADY_TRIGGERED`. Start before `planned`
returns `RUN_NOT_READY`.

Public run and snapshot issue codes are raised in
`src/encode_pipeline/api/routes/runs.py` and
`src/encode_pipeline/api/routes/workflows.py`. Artifact and QC reads use
generation cursors; `RUN_ARTIFACT_GENERATION_CHANGED` and
`RUN_QC_METRIC_GENERATION_CHANGED` mean reload from the first page.

## ENCODE run identity

Implementation: `workflow/Snakefile` functions
`_check_analysis_outdir_identity`, `_read_analysis_identity_marker`, and
`_check_effective_config_record`. Human table:
`docs/configuration.md` section "Run identity and effective configuration".

Records, under the configured outdir:

- `.helixweave/analysis_identity.json`
- `.helixweave/effective_config.json`

The marker stores `marker_schema_version`, `analysis_profile`,
`scientific_policy_revision`, and a hash. The hash covers the validated
config plus parsed samples, excluding only `outdir`, `threads`, and the
`samples` path. `scientific_policy_revision` is internal
(`src/encode_pipeline/config/scientific_policy.py`); users cannot set it.
Input paths are canonical absolute locations.

| State | Behavior |
| --- | --- |
| New directory, or only `.snakemake` and `.helixweave` | Marker is created on a real run, then scientific tasks run |
| Marker matches profile, policy revision, and identity, and the effective-config record recomputes to that identity | Resume is allowed |
| Scientific outputs exist and there is no marker | Refused. Keep the directory. Use a fresh outdir |
| Outdated marker, other profile, other policy revision, or different identity | Refused. Use a fresh outdir |
| Symlink or non-regular metadata | Refused. Do not repair or replace it |
| Marker valid but effective-config record missing, unreadable, or inconsistent | Refused. Do not reconstruct it from the current config |

Dry-run reports the same refusals and writes nothing.

`analysis_profile` is `legacy` or `updated`. Omitted configs resolve
`config_schema_version` to `1.0.0`, whose profile default is `legacy`.
Schema `2.0.0` defaults the profile to `updated`. Do not infer a profile
from a missing key without running the validator. Read
`src/encode_pipeline/config/validator.py` (`_validate_analysis_profile`)
and `workflow/schemas/config.schema.yaml`.

## ENCODE QC contract

Authoritative code: `src/encode_pipeline/qc/contract.py` and
`src/encode_pipeline/adapters/encode_qc.py`. Human contract:
`docs/output-contract.md` section "QC Contract (v2)". Interpretation
guide: `docs/qc-interpretation.md`. Policy prose:
`docs/configuration.md`, `docs/assay-policy.md`.

Distinguish:

- numeric `0`: a measured zero only when the accompanying status says the
  value is defined (`ok`, `reported`, or another status that contract code
  allows to carry a number)
- `NA`: the cell is absent or undefined. Do not coerce it to zero
- `not_applicable`: the metric does not apply (for example SE template
  metrics, or SE `max_insert`)
- `not_enabled`: the producing switch is off
- `missing`: the producer was expected and the input is absent
- `parse_failed`: the source was present but rejected by the strict reader
- `no_data`, `tool_failed`, `unresolved_contigs`, `no_mapped_reads`:
  distinct contract statuses; read the column's allowed set before renaming
  them

`result_manifest.tsv` `qc_flag` is always `NA`. Automatic threshold
conclusions are not implemented. Checklist numbers in
`docs/qc-interpretation.md` are context-dependent guidelines, not admission
gates. More peaks, a green test, or run status `succeeded` is not biological
acceptance.

`duplicate_analysis: auto` retains duplicates for CUT&Tag and excludes them
for ChIP-seq, ATAC-seq, and MNase-seq. The CUT&Tag retain default is still a
candidate pending the SCI-7 real-library data freeze. Say that when the
resolved policy is `retain` because of `auto`. Do not describe an older
inconclusive result as a frozen default.

Per-sample resolved policy is `analysis_policy` JSON (duplicate policy, PE
insert bound, organelle mode and status, blacklist mode and application,
input stages, stage counts). Read it. Do not recompute it. With
`blacklist_mode: analysis`, `workflow/Snakefile` sets the analysis root to
`<outdir>/analysis_filtered` and writes the manifest at
`<analysis-root>/multiqc/result_manifest.tsv`. Raw `final.bam` stays under
the configured outdir. `src/encode_pipeline/manifest/make.py` labels several
of those rows `analysis_filtered_*`. Cite the stored `output_type` and
`path`. Platform `analysis_branch` is artifact metadata, not a rename of
that type.

## Evidence to keep

For a run that actually happens, record the source SHA and whether the
worktree was clean, the config path, validated `analysis_profile` and
`config_schema_version`, checkout `scientific_policy_revision`, outdir,
and whether both run records exist. Also record the interpreter and tool
paths and versions used by that command. For a platform run also record
`workflow_id`, `snapshot_id`, `run_id`, and the artifact or QC generation.

Do not copy those absolute paths into skill instructions or into a public
API response. Platform responses omit workspace paths, environment
values, exception text, and private payloads.
