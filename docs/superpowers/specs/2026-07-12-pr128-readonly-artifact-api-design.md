# PR128 Read-only Artifact API Design

## Scope

PR128 exposes artifact references already persisted by PR127. It never reads the
run workspace, reruns extraction, changes lifecycle state, or serves artifact
file contents. SQLite remains the only source of artifact metadata.

The API adds:

- `GET /api/v1/runs/{run_id}/artifacts` with `operation_id=listRunArtifacts`
- `GET /api/v1/runs/{run_id}/artifacts/{artifact_id}` with
  `operation_id=getRunArtifact`

Artifact download, QC parsing, frontend artifact UI, workflow changes, and new
adapters remain out of scope.

## Alternatives considered

1. Load every artifact and paginate in the route. This preserves the current
   repository API but permits unbounded database and process memory use.
2. Use offset pagination. This is bounded but cursors are unstable if a
   post-success retry atomically replaces the set.
3. Use a run-scoped artifact-ID cursor backed by artifact-ID ordering.
   This is the selected design because it is bounded, matches existing
   event/log cursor semantics, and never queries outside the requested run.

## Query contract

The list endpoint accepts `after` and `limit`. `limit` defaults to 50 and is
bounded to 1 through 100. The repository orders artifacts by deterministic
`artifact_id` ascending, so the order remains stable even when a retry replaces
rows and assigns new SQL surrogate IDs. When a page has more results,
`next_cursor` is the final returned artifact ID.

The `after` artifact must exist in the same run. A missing or cross-run cursor
returns HTTP 400 with `RUN_ARTIFACT_CURSOR_NOT_FOUND`. A replacement between
pages can invalidate a cursor and safely produces the same response.

The detail query always filters by both `run_id` and `artifact_id`. A missing or
cross-run artifact returns the same HTTP 404 `RUN_ARTIFACT_NOT_FOUND` response,
so the endpoint does not disclose ownership by another run. A missing run
returns the existing HTTP 404 `RUN_NOT_FOUND` issue. An existing run with no
artifacts returns an empty successful list.

## Response projection and disclosure boundary

FastAPI response models project `RunArtifactRef` without temporary dict
serialization. The response contains the public artifact ID, requested run ID,
`artifact_type=file`, safe basename, opaque run-scoped URI, MIME type,
`produced_at`, explicit safe `relative_path`, `output_type`, and `size_bytes`,
plus controlled metadata.

Controlled metadata contains the existing ENCODE catalog and logical dimensions:
`catalog_id`, `scope`, `level`, `sample_id`, `experiment_id`, `assay`, `target`,
`genome`, `method`, and `qc_flag`. Other persisted keys are omitted.

The projection validates that:

- `relative_path` is canonical POSIX, begins with `results/`, and has no empty,
  dot, traversal, NUL, backslash, absolute, home, or drive-like form;
- IDs and output types use bounded logical-token syntax;
- URI exactly equals the run-scoped opaque URI derived from the requested run
  and artifact ID, never `file://`;
- names contain no separators or control characters;
- MIME uses a bounded `type/subtype` token;
- metadata strings contain no control characters or absolute/path URI form;
- size is a non-negative integer.

An unsafe persisted record fails the complete list/detail response with a
stable, path-free `RUN_ARTIFACT_DATA_INVALID` issue. No raw exception, database
path, workspace path, command, or environment data reaches the client.

## Service and persistence boundary

`RunService` owns `list_artifacts(run_id, after, limit)` and
`get_artifact(run_id, artifact_id)`. `RunRepository` implements both for
InMemory and SQLAlchemy persistence. SQLAlchemy resolves the cursor and detail
with predicates containing both identifiers and orders list rows by
`artifact_id`. No schema migration is needed because the current
`(run_id, artifact_id)` unique constraint supplies the required index.
The repository and service retain the legacy no-limit behavior when called
without `limit`; the HTTP list route always supplies its bounded `limit + 1`,
so the public contract is bounded without truncating internal PR127 callers.

The API performs no writes and has no dependency on adapter, workspace, worker,
or filesystem services.

## Verification

Tests cover list, detail, empty sets, stable ordering, pagination limits and
cursors, unknown run, unknown artifact, cross-run artifact/cursor isolation,
unsafe persisted-data redaction, InMemory and SQLAlchemy behavior, SQLite
reopen, operation IDs, OpenAPI export, Orval regeneration, frontend checks, and
the complete Python suite. The final diff review focuses on authorization-like
run isolation, disclosure, pagination, and generated contract drift.
