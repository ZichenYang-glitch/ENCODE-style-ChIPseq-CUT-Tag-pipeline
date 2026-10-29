# Platform run

Read [../../_shared/entries.md](../../_shared/entries.md) before choosing a
command. This note routes a run through the public lifecycle for a workflow
id the user already confirmed. Commands and routes absent from that inventory
stay unused.

The platform stack is `workflow/envs/ci-fast.lock` plus the runtime in
`docs/development/local-platform-runtime.md`. Doctor, schema, validate,
create, start, status, cancel, and logs are that stack. A direct Snakemake
run does not install it.

Member routes require a principal. Methods other than GET, HEAD, and OPTIONS
also require CSRF. Use the authenticated local platform session. Leave
private Python helpers, direct SQLite edits, and Redis enqueues unused.

SQLite is lifecycle truth. Redis/RQ is only the queue. Report status from
`GET /api/v1/runs/{run_id}`. Public responses omit workspace paths,
environment values, exception text, and the private input payload. Quote
public status and issue codes.

`succeeded` means execution completed. It is not biological acceptance.

## Doctor

From a verified source checkout:

```bash
python3 -I -S scripts/checkout_bootstrap.py \
  --repository-root . local-platform --doctor
```

`helixweave --doctor` is the installed equivalent. Requirements and the
meaning of a missing Bulk RNA-seq binding are in
`docs/development/local-platform-runtime.md`. Doctor exits before opening
ports. A missing Bulk RNA-seq binding reported as `not_configured` or
`unavailable` does not fail the core doctor.

## Authoring schema

```http
GET /api/v1/workflows/{workflow_id}/schema
```

Build the validate body from that schema. For any id other than
`encode-style-chipseq-cuttag-atac-mnase`, this ENCODE skill stops; follow
the schema you just read. Enabled reference profiles, when the user is
binding one, come from:

```http
GET /api/v1/workflows/{workflow_id}/reference-profiles
```

## Validate

```http
POST /api/v1/workflows/{workflow_id}/validate
```

The body is a `ValidationRequest`: `config`, `samples`, and `options`.
`samples` is inline rows or a server path, as the schema's input modes
declare. Add `project_id` together with `sample_revision_ids`,
`input_selections`, or `reference_profile_revision_id` only when the user is
using those bindings. The route rejects unknown fields.

A successful response carries `snapshot.snapshot_id`. Keep that id. Also
record `snapshot.schema_version` as the authoring contract version. Record
`analysis_profile` and `config_schema_version` from the accepted config.
Those are different fields. Read `scientific_policy_revision` from
`src/encode_pipeline/config/scientific_policy.py` in this checkout.

When validation returns `REFERENCE_PROFILE_REQUIRED`, list reference profiles
and validate again with the revision the user chose. Expired or stale
snapshot errors mean: validate again and use the new `snapshot_id`.

Creating a run from raw config, or from anything other than that server-owned
`snapshot_id`, bypasses snapshot validation. Do not do it.

## Create

Do this when the user asked to create or submit a platform run. Creation
stores a run. It does not execute it.

```http
POST /api/v1/workflows/{workflow_id}/runs
```

Body: `snapshot_id` copied from the validate response, plus optional `tags`.
The first success is HTTP 201. An identical retry returns that snapshot's
canonical run with HTTP 200. A retry that disagrees with the stored run,
including its tags, is `VALIDATED_SNAPSHOT_REPLAY_CONFLICT`; stop and report
it.

Record `workflow_id`, `snapshot_id`, and `run.run_id`.

## Reach planned

Create leaves the run `created`. `POST /api/v1/runs/{run_id}/start` then
returns 409 `RUN_NOT_READY`: "Run must complete preflight before it can
start."

Call `POST /api/v1/runs/{run_id}/preflight`
(`src/encode_pipeline/api/routes/preflight.py`). It accepts a `created` run,
moves it to `validating`, and runs planning in the background. Poll
`GET /api/v1/runs/{run_id}` until the status is `planned` or a terminal
status. The maintained description of that planning step is
`docs/development/local-platform-runtime.md`, section "Submission and
execution lifecycle".

Leave Redis and private helpers unused. A second preflight on a run that has
left `created` returns 409 `PREFLIGHT_ALREADY_TRIGGERED`.

## Start

Call start only after the user explicitly authorizes execution and the run
status is `planned`.

```http
POST /api/v1/runs/{run_id}/start
```

HTTP 202 means the queue accepted the stable job and SQLite records `queued`.
HTTP 409 `RUN_NOT_READY` means preflight has not reached `planned`. HTTP 409
for a missing or changed workflow build identity means: stop and report;
validate again rather than patching the run. HTTP 503 means the queue or
execution backend is unavailable. Do not enqueue the job yourself.

## Status

```http
GET /api/v1/runs/{run_id}
```

Status values in `encode_pipeline.platform.runs.RunStatus` are `created`,
`validating`, `planned`, `queued`, `running`, `succeeded`, `failed`, and
`cancelled`. Events are `GET /api/v1/runs/{run_id}/events`. Artifact and QC
reads use generation cursors. `RUN_ARTIFACT_GENERATION_CHANGED` and
`RUN_QC_METRIC_GENERATION_CHANGED` mean reload from the first page. Record
the generation when you read those routes.

`succeeded` means the worker finished. Biological acceptance is a separate
decision. QC interpretation is outside this skill.

## Cancel

Cancel only after the user explicitly authorizes it.

```http
POST /api/v1/runs/{run_id}/cancel
```

Read `GET /api/v1/runs/{run_id}` afterward and report that SQLite status. A
running run can return HTTP 202 while the stop is being acknowledged. HTTP
409 means durable state rejected the cancel. HTTP 503 means the queue did
not confirm it.

## Logs

```http
GET /api/v1/runs/{run_id}/logs
```

Query `stream_name` (`stdout` by default; `stderr` is the other process
stream), `after`, and `limit`. Follow `next_cursor`. An operator trace of
the run may keep the source SHA and the tool path that was used. The
summary returned through this API quotes status and public issue codes
only. Leave workspace paths and private payloads out of that summary.

The inventory's artifact download route stays unused unless the user
explicitly authorizes a download.
