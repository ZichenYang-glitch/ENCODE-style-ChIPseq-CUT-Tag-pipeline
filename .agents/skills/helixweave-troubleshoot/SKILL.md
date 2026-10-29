---
name: helixweave-troubleshoot
description: >-
  Locate a HelixWeave configuration, schema, provenance, conda, Snakemake
  profile, outdir-identity, scheduler, worker, cancellation, scientific-tool,
  or QC-parse failure and choose the smallest confirming check. Use when a
  validator, doctor, run, worker, or artifact read fails. Do not use for a
  first-time run setup or for explaining QC numbers that already parsed.
---

# HelixWeave troubleshoot

Classify the symptom before proposing a command. Read
[../_shared/entries.md](../_shared/entries.md) and pick the lowest-cost check
that can confirm that class. Command text stays in that inventory. The
confirming command and stop condition for each class are in
[references/symptom-routes.md](references/symptom-routes.md).

`scripts/checkout_bootstrap.py` accepts only `verify-checkout`,
`installed-artifact`, `pytest`, `openapi`, `validate`, and `local-platform`.
Both `-I` and `-S` are required. There is no `snakemake` subcommand.
`--repository-root` is the checkout that contains that bootstrap.

If the captured output already names the class, that output is the
confirmation. Run one inventoried check only when the class is still
ambiguous. Re-use the surface that failed. A full workflow rerun and the
Protected Bulk Gate are not checks.

First-time setup stays on `docs/quickstart.md`. QC values that already
parsed belong to [helixweave-qc](../helixweave-qc/SKILL.md).

Lifecycle, validated snapshots, artifact and QC cursors, and error
envelopes are workflow-neutral. Sample columns, assays, output paths, and
ENCODE QC tables apply only to workflow id
`encode-style-chipseq-cuttag-atac-mnase`. For any other id, read
`GET /api/v1/workflows/{workflow_id}/schema` and that adapter's own errors.

## Config, schema, or alias

`ValidationError` from `encode-validate`, `scripts/validate_samples.py`, or
`POST /api/v1/workflows/{workflow_id}/validate`.
`scripts/validate_samples.py` does not see later `snakemake --config`
overrides; the workflow validates the merged config when it loads.

`remove_dup` and `duplicate_analysis` are the common alias conflict
(`docs/configuration.md`, "Duplicate policy"). `yes` maps to `exclude` and
`no` maps to `retain`. Disagreeing explicit values raise `ValidationError`.

`analysis_profile` is `legacy` or `updated`. An omitted key is not evidence
of either value. An omitted `config_schema_version` resolves to `1.0.0`,
and only schema `2.0.0` defaults the profile to `updated`. Read
`src/encode_pipeline/config/validator.py` (`_validate_analysis_profile`,
`_validate_duplicate_settings`, `validate_execution_support`). Legacy
execution support keeps `control_mode=pooled`, an empty
`duplicate_backend`, `fragment_cpm=off`, and `tss_fragment_profile=false`.
Updated execution support requires `duplicate_backend=picard` and
`control_mode=matched`. Path symptoms for this workflow are the
`docs/quickstart.md` Troubleshooting headings for FASTQ, Bowtie2 index
prefix, and control path.

## Source or import provenance

Exit status 2 and one line,
`source provenance check failed [reason_code]`, from
`scripts/checkout_bootstrap.py`, or `SourceProvenanceError` from
`scripts/source_provenance.py`. Reason codes include
`bootstrap_startup_unsafe`, `bootstrap_source_invalid`,
`repository_root_invalid`, `pytest_plugin_unsafe`,
`pytest_plugin_missing`, `pth_mapping_unsafe`, and
`startup_hook_unsafe`. The failure list is
`docs/development/source-provenance.md`. Retry with this checkout's own
bootstrap under `python3 -I -S`. Leave unrelated environments installed.
Do not point `PYTHONPATH` at another worktree. Do not edit the scientific
config to clear this line.

## Conda, profile, or missing tool

The small runner is `workflow/envs/runner.lock`, `snakemake` on `PATH`, and
`--use-conda`. Rule tools come from `workflow/envs/`. A missing rule
environment is a conda failure for that rule. It is not a reason to
recreate the platform stack.

The platform ci-fast stack is Python 3.12, Snakemake 8.30.0, Redis 7 or
newer, and Node from `.nvmrc`. Doctor requires Python 3.12, Snakemake
8.30.0, Redis 7.x, and Node.js 22 or newer. The check is
`local-platform --doctor` or `helixweave --doctor`
(`src/encode_pipeline/cli/local_platform.py`,
`docs/development/local-platform-runtime.md`). Doctor exits before opening
ports or creating runtime data. A missing Bulk RNA-seq binding reported as
`not_configured` or `unavailable` leaves the core doctor passing.

`profiles/default` is a Snakemake profile. `analysis_profile` is the
scientific config key. The profile file sets `configfile`, `cores`,
`use-conda`, `rerun-incomplete`, `keep-going`, and `show-failed-logs`.

Runner symptoms are the `docs/quickstart.md` headings "Slow micromamba
install" and "`snakemake: command not found`". "Missing tools for tiny
execution" applies to that test only.

## Outdir identity

`WorkflowError` from `workflow/Snakefile` `_check_analysis_outdir_identity`
and the helpers it calls. Quote the class:

- unmarked outdir: `Refusing to write analysis outputs into an unmarked outdir`
- different identity: `Existing analysis outdir has a different analysis identity`
- outdated marker: `Existing analysis identity marker uses an outdated format`
- symlink: `must not be a symlink`
- missing effective config: `Effective configuration record is missing next to a valid run marker`
- fresh outdir: `Use a fresh outdir`

The same helpers refuse a non-regular metadata file, a profile mismatch
(`owned by analysis_profile`), and a different scientific policy revision.
`Use a fresh outdir` means a new directory. Keep the existing directory,
`.snakemake`, `.helixweave/analysis_identity.json`, and
`.helixweave/effective_config.json`. The confirming check for a direct
ENCODE Snakemake invocation is a dry-run. It does not execute scientific
rules and creates neither run record nor scientific outputs. Snakemake may
still create `.snakemake` scheduler metadata.
For a platform run, read the run events. A successful preflight includes
that dry-run.

## Scheduler, worker, lifecycle, cancel, or timeout

Statuses in `encode_pipeline.platform.runs.RunStatus`: `created`,
`validating`, `planned`, `queued`, `running`, `succeeded`, `failed`,
`cancelled`. `succeeded` means execution completed. SQLite is lifecycle
truth. Redis/RQ is the queue. Read `GET /api/v1/runs/{run_id}`,
`/events`, and `/logs`.

Public issue codes in `src/encode_pipeline/api/routes/runs.py`:
`RUN_NOT_READY`, `RUN_START_CONFLICT`, `RUN_QUEUE_UNAVAILABLE`,
`RUN_WORKFLOW_BUILD_IDENTITY_MISSING`, `RUN_WORKFLOW_BUILD_CHANGED`,
`RUN_CANCELLATION_CONFLICT`, `RUN_CANCELLATION_UNAVAILABLE`,
`VALIDATED_SNAPSHOT_NOT_FOUND`, `VALIDATED_SNAPSHOT_EXPIRED`,
`VALIDATED_SNAPSHOT_STALE`, `VALIDATED_SNAPSHOT_REPLAY_CONFLICT`.

When the failure write commits, a worker hard timeout records public
issue `RUN_WORKER_FAILED` with event `reason_code` `WORKER_JOB_TIMEOUT`.
A process-runner timeout records public issue `RUN_EXECUTION_FAILED`
with event `reason_code` `PROCESS_RUNNER_TIMEOUT`. Constraints are
"Submission and execution lifecycle", "Explicit administrator recovery",
and "Running cancellation semantics" in
`docs/development/local-platform-runtime.md`. `helixweave admin run` is
only for the host operator, and only when that recovery section says the
observed state is eligible. Leave unrelated workers running.

## Scientific tool non-zero exit

Read the Snakemake log for that rule. `profiles/default/config.yaml` sets
`show-failed-logs: true`. Re-run that incomplete job with
`--rerun-incomplete` only after the user authorizes a resume of the same
identity. Leave the rest of the matrix in place.

## QC or artifact missing or parse_failed

`missing` means the producer was expected and the file is absent.
`parse_failed` means the strict reader rejected a present file
(`src/encode_pipeline/qc/contract.py`). Leave `NA` as `NA`. A parse
failure stays a parse failure; the workflow is not re-run to change it.
After the file is known intact, interpretation belongs to helixweave-qc.
`RUN_ARTIFACT_GENERATION_CHANGED` and
`RUN_QC_METRIC_GENERATION_CHANGED` mean reload the list from the first
page.

## Protect results

Keep existing outputs. Never delete an outdir, remove `.snakemake`,
replace a marker, or hand-edit `effective_config.json`. Never weaken
`--strict-inputs` or schema admission to make a run pass.

## Production defect

When the evidence shows a deterministic production defect, save the
smallest repro: the command, the public error code, and config keys that
are not private paths. Ask the user to report it. This skill does not
patch production code.
