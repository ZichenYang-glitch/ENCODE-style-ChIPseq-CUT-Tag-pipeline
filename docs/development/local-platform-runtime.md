# Durable Local Platform Runtime

The local workflow-platform API stores run lifecycle state in file-backed
SQLite. This is the single-process, single-user persistence layer for the
current MVP; it is not a queue worker or a real execution backend.

## Run locally

Install the API extra and start the FastAPI factory:

```bash
python -m pip install -e ".[api]"
uvicorn encode_pipeline.api.main:create_app --factory --reload
```

By default, the application uses:

```text
~/.encode-pipeline/platform.db
~/.encode-pipeline/workspaces/
```

The app factory applies the bundled Alembic migrations before accepting
requests. Run records, events, stdout/stderr chunks, and artifact references
survive a normal API process restart.

Set `ENCODE_PIPELINE_DATABASE_URL` to use another file-backed SQLite database:

```bash
export ENCODE_PIPELINE_DATABASE_URL="sqlite:////absolute/path/platform.db"
uvicorn encode_pipeline.api.main:create_app --factory
```

Only file-backed SQLite URLs are supported in this phase. An in-memory SQLite
database is deliberately rejected because it cannot preserve state across the
migration and API connections.

## Restart semantics

At application construction, persisted runs in `validating`, `queued`, or
`running` are marked `failed` with `RUN_INTERRUPTED_BY_API_RESTART`. The API
also writes a `run_recovered_after_restart` event. This is intentionally
conservative: the preflight background task belongs to the prior API process,
so retaining an active status after that process exits would be misleading.

`created` and `planned` runs remain unchanged because they do not imply an
active process. Terminal runs remain unchanged. Re-running startup recovery is
idempotent.

## Current boundary

The API's background preflight remains a short, local dry-run path. It is not
a durable job worker, and this release does not implement real Snakemake
execution, process cancellation, Redis/RQ, authentication, HPC scheduling, or
artifact extraction. Those responsibilities require the later worker and
execution milestones.
