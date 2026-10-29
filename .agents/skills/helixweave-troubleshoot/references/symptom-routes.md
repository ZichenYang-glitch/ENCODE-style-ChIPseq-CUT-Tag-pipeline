# Symptom routes

Checkout-root commands are the ones in
[../../_shared/entries.md](../../_shared/entries.md). Run the confirming
command for the class already chosen in
[../SKILL.md](../SKILL.md). Stop at the condition below. One confirmation
is enough.

## Config, schema, or alias

Confirm with the surface that raised `ValidationError`:

```bash
python3 -I -S scripts/checkout_bootstrap.py --repository-root . \
  validate --config CONFIG
```

`encode-validate` and `POST /api/v1/workflows/{workflow_id}/validate` are
the same contract on their own surfaces. Add `--strict-inputs` only when
the symptom is a missing declared input. A key, alias, or
`analysis_profile` error is confirmed without that flag.

Stop when the error names `remove_dup` / `duplicate_analysis`, an
`analysis_profile` constraint, or another config key. Read
[../../../../docs/configuration.md](../../../../docs/configuration.md)
("Duplicate policy") and
[../../../../src/encode_pipeline/config/validator.py](../../../../src/encode_pipeline/config/validator.py).
ENCODE path headings are in
[../../../../docs/quickstart.md](../../../../docs/quickstart.md)
(Troubleshooting). Another workflow id uses its own schema from
`GET /api/v1/workflows/{workflow_id}/schema`.

A provenance line from this command moves the case to the next section.

## Source or import provenance

```bash
python3 -I -S scripts/checkout_bootstrap.py --repository-root . \
  verify-checkout
```

Stop on exit status 2 and one
`source provenance check failed [reason_code]` line, including
`bootstrap_startup_unsafe`, `bootstrap_source_invalid`,
`repository_root_invalid`, `pytest_plugin_unsafe`,
`pytest_plugin_missing`, `pth_mapping_unsafe`, or
`startup_hook_unsafe`. Follow that line. The code list is
[../../../../docs/development/source-provenance.md](../../../../docs/development/source-provenance.md).
Implementation:
[../../../../scripts/checkout_bootstrap.py](../../../../scripts/checkout_bootstrap.py)
and
[../../../../scripts/source_provenance.py](../../../../scripts/source_provenance.py)
(`SourceProvenanceError`). The fix is this checkout's bootstrap under
`python3 -I -S`, after the printed guidance. Unrelated environments stay
installed. Do not point `PYTHONPATH` at another worktree.

## Conda, profile, or missing tool

Direct workflow: the inventoried dry-run, or the log of the invocation
that already failed. Stop when the message is `snakemake: command not found`
or a missing environment for one rule. Runner lock:
[../../../../workflow/envs/runner.lock](../../../../workflow/envs/runner.lock).
Headings: "Slow micromamba install" and "`snakemake: command not found`" in
[../../../../docs/quickstart.md](../../../../docs/quickstart.md).
"Missing tools for tiny execution" is that test's prerequisite. A rule-env
failure stops at that rule. The platform environment stays as it is.

Platform stack:

```bash
python3 -I -S scripts/checkout_bootstrap.py --repository-root . \
  local-platform --doctor
```

Installed equivalent: `helixweave --doctor`. Stop when an environment
check fails, or when doctor passes and the failure is outside this stack.
Failed checks include Python 3.12, Snakemake 8.30.0, Redis 7.x, and
Node.js 22 or newer (pin file `.nvmrc`). Doctor creates no runtime data
and opens no ports. Requirements:
[../../../../src/encode_pipeline/cli/local_platform.py](../../../../src/encode_pipeline/cli/local_platform.py)
and "Run locally" in
[../../../../docs/development/local-platform-runtime.md](../../../../docs/development/local-platform-runtime.md).
`not_configured` or `unavailable` for Bulk RNA-seq is not a core failure.

[../../../../profiles/default/config.yaml](../../../../profiles/default/config.yaml)
is a Snakemake profile. `analysis_profile` is the scientific key.

## Outdir identity

Direct ENCODE Snakemake only. Use the directory where `CONFIG`'s relative
paths resolve, and the Snakefile from the checkout:

```bash
snakemake -s "$ROOT/workflow/Snakefile" --configfile CONFIG -n --use-conda
```

Stop on `WorkflowError` from
[../../../../workflow/Snakefile](../../../../workflow/Snakefile)
(`_check_analysis_outdir_identity` and its helpers), including unmarked
outdir, different identity, outdated marker, symlink, or missing effective
config. The dry-run writes neither `analysis_identity.json` nor
`effective_config.json`. The existing directory and both records stay.
`Use a fresh outdir` is a new directory. Human table:
[../../../../docs/configuration.md](../../../../docs/configuration.md)
("Run identity and effective configuration").

A platform run stops on a preflight event from `GET /api/v1/runs/{run_id}/events`
that quotes an identity `WorkflowError`. Absent that event, this is not
the outdir class. A second Snakemake process is not pointed at that workspace.

## Scheduler, worker, lifecycle, cancel, or timeout

```text
GET /api/v1/runs/{run_id}
GET /api/v1/runs/{run_id}/events
GET /api/v1/runs/{run_id}/logs
```

Start with the run. Read events and logs when the status or public code
is not one of the listed issue codes, or when the public issue is
`RUN_WORKER_FAILED` or `RUN_EXECUTION_FAILED`. Stop when SQLite status plus
one public code identifies the window. Issue strings:
[../../../../src/encode_pipeline/api/routes/runs.py](../../../../src/encode_pipeline/api/routes/runs.py).
Lifecycle, recovery eligibility, and cancellation:
[../../../../docs/development/local-platform-runtime.md](../../../../docs/development/local-platform-runtime.md)
("Submission and execution lifecycle", "Explicit administrator recovery",
"Running cancellation semantics").

`helixweave admin run` is the next step only when the user is the host
operator and that recovery section calls the observed state eligible.
Diagnose is the read. `run fail` and `run requeue` wait for that
eligibility. Redis is not the lifecycle record. Other workers stay up.
Cancel, start, and resume wait for an explicit user authorization.

## Scientific tool non-zero exit

Read the Snakemake log already printed for that rule.
[../../../../profiles/default/config.yaml](../../../../profiles/default/config.yaml)
sets `show-failed-logs: true`.

Stop when the log names the rule and a non-zero tool status. The same
identity may resume with `--rerun-incomplete` only after the user
authorizes that resume. The rest of the matrix stays in place. This is
the direct-workflow command in
[../../_shared/entries.md](../../_shared/entries.md).

## QC or artifact missing or parse_failed

```text
GET /api/v1/runs/{run_id}/artifacts
GET /api/v1/runs/{run_id}/qc-metrics
```

Stop when the row says `missing` (expected producer, file absent) or
`parse_failed` (file present, strict reader rejected it). Contract:
[../../../../src/encode_pipeline/qc/contract.py](../../../../src/encode_pipeline/qc/contract.py)
and
[../../../../docs/output-contract.md](../../../../docs/output-contract.md)
("QC Contract (v2)"). `NA` stays `NA`. The workflow is not re-run to
change a parse. An intact file goes to
[../../helixweave-qc/SKILL.md](../../helixweave-qc/SKILL.md) for
interpretation. `RUN_ARTIFACT_GENERATION_CHANGED`
([../../../../src/encode_pipeline/api/routes/artifacts.py](../../../../src/encode_pipeline/api/routes/artifacts.py))
and `RUN_QC_METRIC_GENERATION_CHANGED`
([../../../../src/encode_pipeline/api/routes/qc_metrics.py](../../../../src/encode_pipeline/api/routes/qc_metrics.py))
stop the current page; reload from the first page.

## Deterministic production defect

Stop when the confirming check contradicts the documented contract for
the same public inputs. Save the command, the public error code, and
config keys that are not private paths. Ask the user to report it.
Production code stays unchanged.
