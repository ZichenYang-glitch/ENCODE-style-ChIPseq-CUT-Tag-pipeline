# PR135 Adapter-Owned Authoring Contract Design

## Objective and boundary

PR135 makes workflow input discovery a versioned adapter contract and lets the
ENCODE-style adapter consume inline sample rows without changing scientific
validation or run lifecycle behavior. It establishes the backend contract that
PR136 can render and that PR137 can submit:

```text
adapter JSON Schema + input modes + limits
  -> config object + inline sample rows + options
  -> existing scientific validation
  -> durable run input snapshot
  -> deterministic workspace config/samples files
```

The adapter continues to own workflow-specific fields and scientific meaning.
The platform owns workflow-neutral transport limits, schema projection, request
handling, and workspace path containment. This PR does not add the new-run UI,
validated snapshot submission semantics, another adapter, or a new scientific
validator.

## Baseline findings

`WorkflowSchema` currently returns unversioned dictionaries explicitly marked
as hints. The ENCODE adapter rejects `list[dict[str, str]]` samples even though
`WorkflowInputs`, durable run snapshots, worker reconstruction, and the
workspace planner already preserve that shape. The scientific config validator
requires an existing sample-sheet path, and the sample loader is the canonical
source for row, cross-row, replicate, and strict-input validation.

The public validation endpoint currently forwards the adapter-private success
value. Doing that for a temporary inline-row bridge would disclose an absolute
temporary path. The endpoint is also async while calling synchronous scientific
validation. Request models constrain types but do not bound actual body bytes,
row counts, column counts, or cell lengths.

## Considered approaches

### Add a second in-memory scientific validator

Duplicating row and cross-row rules would immediately create two scientific
contracts. It is rejected.

### Refactor the scientific loader around an in-memory parser

This could remove the temporary file, but it changes a mature scientific
validation layer and substantially expands PR135. It is deferred unless the
bounded bridge proves inadequate.

### Render inline rows into a private temporary TSV and reuse the validator

The adapter writes a fixed canonical external-column header in a unique
mode-0700 temporary directory, invokes the existing `validate_config` and
`load_and_validate_samples` functions, retains only normalized in-memory
results, and removes the directory in all outcomes. This is selected. It
preserves one scientific authority and does not put the temporary file in the
run workspace, SQLite, events, issues, or HTTP response.

## Versioned schema contract

`WorkflowSchema` becomes a frozen, defensively copied contract with:

- `schema_version`: `1.0.0`;
- `schema_dialect`: `https://json-schema.org/draft/2020-12/schema`;
- per-surface `coverage` values (`partial` or `complete`);
- adapter-declared `authoring_modes` and wire-level `input_modes`;
- bounded `limits`;
- standard JSON Schema documents for config, sample rows, and options.

Each JSON Schema has the same explicit `$schema`, a stable versioned `$id`, and
JSON-safe content. Config and options schemas have object roots. The sample
schema has an array root whose items use external TSV column names.

The ENCODE contract declares:

| Surface | Coverage | Authoring modes | Input modes |
| --- | --- | --- | --- |
| config | partial | schema_form, yaml | object |
| samples | complete for inline rows | tsv_upload, inline_table | inline_rows, server_path |
| options | complete | schema_form | object |

Config is deliberately `partial` with `additionalProperties: true`: the first
form covers every field used by the deterministic tiny profile, but does not
claim to render every advanced validator setting. It omits `samples`, because
sample rows are a sibling adapter input. Unknown advanced config fields remain
representable for the later shared Form/YAML state. The sample schema enumerates
all 17 canonical columns and uses `additionalProperties: false`; the adapter
therefore rejects unknown inline columns. Options accepts only the boolean
`strict_inputs` field.

`SchemaResponse` exposes a typed `schema` value and removes the misleading
`schema_hints` field. Unknown workflows return `schema: null` in the existing
redacted 404 envelope. The existing frontend compatibility adapter performs the
smallest mapping needed to keep the current workflow detail route working;
PR136 will consume the generated schema operation directly.

## Input priority and scientific validation

Input priority is unambiguous:

1. A sample list means inline mode, even when it is empty or invalid.
2. Inline mode never reads or falls back to `config.samples`.
3. A string sample payload overrides `config.samples` as the existing local
   server-path mode.
4. With no sibling samples payload, the legacy `config.samples` path remains
   available for local/API compatibility.

Empty inline rows, unknown columns, forbidden NUL/tab/newline characters, and
structural bound violations fail closed before scientific validation. Valid
rows are rendered with a fixed full header. The temporary directory remains
alive across config and sample validation and is removed in `finally` cleanup.
The returned internal validated config removes the temporary samples path.

After scientific validation the adapter invokes the existing external-input
policy. Relative FASTQ, index, control BAM, and genome-resource paths therefore
produce the same safe failure during validate and workspace planning rather
than a misleading validate success followed by preflight failure.

`plan_workspace` consumes normalized rows and retains the established renderer:
`config/config.yaml` always contains `samples: config/samples.tsv`, and
`config/samples.tsv` uses the canonical normalized representation. Different
temporary names cannot affect either byte stream.

## Bounded transport and domain inputs

The workflow-neutral hard ceiling and ENCODE-declared limits are:

- request body: 2,097,152 UTF-8 bytes;
- inline sample rows: 1–1,000;
- columns per row: 1–64;
- column name: at most 128 characters;
- cell value: at most 4,096 characters.

Pydantic exposes the row/column/cell constraints in OpenAPI. `WorkflowInputs`
uses the same constants so direct service calls cannot bypass the structural
ceiling. Adapter limits may be no wider than the platform ceiling.

A route-scoped pure-ASGI middleware accumulates at most the limit in one byte
buffer before Pydantic parses POST validate/create requests and rejects an
excess chunk before copying it. It counts actual chunks and never trusts
`Content-Length`; chunked, absent, or falsely small headers cannot bypass the
cap, and empty chunk fragmentation cannot grow a message list. Oversized
validation and run-creation requests return a stable 413
`API_REQUEST_TOO_LARGE` in their respective envelopes before adapter or
repository calls. Other routes are unaffected.

Validate and create become synchronous FastAPI endpoints so temporary I/O,
scientific validation, and SQLite work run in FastAPI's threadpool rather than
on the event loop. Public validation responses do not project adapter-private
values; success is represented by `ok`, issues, and a null value.

## Security and failure semantics

- Temporary paths and content never enter SQLite, events, logs, OpenAPI
  examples, or public responses.
- Inline render and temporary-file failures map to stable adapter issues without
  exception text.
- Inline rows cannot cause fallback reads from a user-supplied server path.
- Existing workspace containment and external path policy remain fail closed.
- Request-limit responses contain no body excerpt, path, environment value, or
  exception string.
- The schema is trusted adapter metadata, copied defensively and checked for
  JSON-safe values; it is not executable input.

## Tests and acceptance

Tests cover Draft 2020-12 meta-schema validity, stable IDs, per-surface
coverage, modes, limits, immutability, and invalid JSON values. The tiny profile
without `config.samples`, its inline row, and its options must pass both the
published schema and the existing adapter validator.

Adapter tests cover precedence, empty/unknown/control-character rows,
scientific failures, external path failures, temporary cleanup, concurrent
validation isolation, and absence of temporary paths in results. Repeated
workspace planning must produce byte-identical files; inline and equivalent
path modes must produce the same canonical sample TSV.

Service/API tests cover create, SQLite reopen, workspace reconstruction and
materialization, every structural limit, safe public projection, and exact/
over-limit multi-chunk bodies with absent or misleading `Content-Length`.
OpenAPI and Orval are mechanically regenerated and checked for drift.

Acceptance also requires the full Python suite, relevant real Redis/RQ/
Snakemake gates, frontend tests/typecheck/build, Playwright regression, Ruff,
format, `git diff --check`, a clean worktree, green GitHub checks, and an
independent read-only merge-gate review.

## Non-goals and residual risks

This PR does not add the PR136 input workbench, PR137 immutable validation
snapshot flow, PR138 browser gate/demo, FASTQ upload, omics-intake runtime
coupling, another adapter, Agent writes, auth, multi-user isolation, HPC,
object storage, or immutable workflow bundles.

The bridge briefly stores sample metadata in an OS-private temporary file. Its
permissions and normal lifetime are bounded, but an uncatchable process death
could leave that private directory for operating-system cleanup; a future
in-memory scientific-loader refactor could remove that exposure. The config
form contract is intentionally partial; advanced YAML fields remain
backend-authoritative until explicitly added to a later schema version.
