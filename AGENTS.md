# Repository Guidance

## Mission

This repository contains a workflow-neutral, local/small-team omics platform
and its first scientific adapter: the ENCODE-style ChIP-seq, CUT&Tag, ATAC-seq,
and MNase-seq Snakemake workflow.

Keep this file limited to durable defaults. Current priorities and delivery
sequence belong in `docs/development/workflow-platform-agent-roadmap.md`.

Use a single agent and a lightweight, risk-driven workflow by default. Do not
invoke Superpowers, swarms, or other multi-agent ceremony unless the user
explicitly requests it or the work is genuinely large and parallel.

## Architecture Boundaries

- `src/encode_pipeline/platform/` owns workflow-neutral domain contracts.
- `src/encode_pipeline/services/` owns workflow-neutral orchestration.
- `src/encode_pipeline/api/` translates HTTP requests and responses; routes
  should stay thin.
- `src/encode_pipeline/persistence/` owns SQLAlchemy repositories and Alembic
  migrations.
- `src/encode_pipeline/workers/` owns Redis/RQ execution mechanics.
- `src/encode_pipeline/adapters/` owns workflow-specific schemas, validation,
  workspace planning, commands, and artifact/QC extraction.
- `workflow/` owns Snakemake rules and scientific behavior.
- `frontend/` owns workflow/run UX and consumes the generated API client.

Preserve these invariants:

- SQLite is the canonical lifecycle and result-metadata store. Redis/RQ is a
  queue and worker boundary, not a second lifecycle source of truth.
- Platform code must not depend on adapter-private validation payloads,
  workflow rule names, sample columns, or output paths.
- Generic API and frontend code must support adapter-owned schema contracts;
  never hard-code the current ENCODE schema version or limits as universal.
- Validated snapshots, durable job identity, lifecycle transitions, and
  cancellation acknowledgement are contracts. Change them deliberately.
- Scientific behavior stays in the adapter or workflow unless it is truly
  reusable across adapters.

## Change Discipline

Read the relevant implementation and tests before editing. Follow existing
frameworks and seams before adding an abstraction.

- Make the smallest coherent change that solves the request.
- Do not combine unrelated refactors, formatting, or roadmap work.
- Small changes need no standalone plan document.
- For medium changes, keep a short in-session plan.
- Write a focused ADR or design note only for material changes to architecture,
  persistence, public contracts, worker semantics, or cross-repo integration.
- Explicit user instructions override roadmap sequencing. Surface conflicts
  instead of silently widening or redirecting scope.

## Contract Rules

- Never hand-edit `frontend/src/api/generated/`.
- When backend OpenAPI changes, run the frontend `openapi:regenerate` script
  and inspect both the specification and generated-client diff.
- Persistence schema changes require an Alembic revision and upgrade coverage
  from the supported prior schema.
- Treat workflow authoring schemas, artifact/QC models, and public error
  envelopes as versioned contracts.
- Preserve fail-closed path handling and redaction. Do not expose workspace
  paths, environment values, exception text, or private payloads through APIs.

## Testing

Run the minimum tests that provide credible evidence for the changed behavior.
Start at the lowest relevant layer and broaden according to blast radius.

- Add tests for real behavior or contracts, not duplicate coverage at every
  unit/API/browser layer.
- Broaden validation for lifecycle, persistence, worker, migration, OpenAPI,
  generated-client, artifact/QC, shared adapter, or browser-flow changes.
- For scientific workflow changes, run the relevant validator, DAG dry-run,
  snapshot, or tiny execution tests.
- Do not weaken assertions, add skips, or delete coverage merely to make a
  change pass.
- Run the full suites only when the change is cross-cutting, a release gate or
  CI requires them, or the user explicitly requests them.

Common commands, run from the repository root:

```bash
python3 -m pytest <targeted paths> -v
python3 -m ruff check <changed paths>
python3 -m ruff format --check <changed paths>
python3 scripts/validate_samples.py --config config/config.yaml
snakemake -s workflow/Snakefile --configfile config/config.yaml -n

npm --prefix frontend test -- --run
npm --prefix frontend run typecheck
npm --prefix frontend run build
npm --prefix frontend run openapi:regenerate
```

Use `npm --prefix frontend run test:e2e` only for relevant browser/runtime
changes. See `docs/development/local-platform-runtime.md` for the real local
stack and its doctor/demo commands.

## Workspace And Git Safety

- Preserve pre-existing user changes, including changes in files you edit.
- Never reset, overwrite, or delete user work, local data, results, secrets,
  environments, or runtime state without explicit authorization.
- Do not merge, push, rebase, tag, rewrite history, or delete branches unless
  the user explicitly requests that action.
- Prefer an isolated worktree when the root tree is dirty or PR work needs a
  clean base.
- Clean up only processes and temporary files that your own work created.

## Definition Of Done

A change is done when its requested behavior is implemented, relevant tests
and contract checks pass, `git diff --check` is clean, and any services started
for verification are stopped. Report tests that were not run and remaining
risks plainly. Do not claim scientific or end-to-end validation beyond the
evidence actually executed.
