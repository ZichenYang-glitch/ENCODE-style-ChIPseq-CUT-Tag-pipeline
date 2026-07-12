# PR127 Adapter Artifact Extraction Design

## Objective

Index the ENCODE-style workflow's existing manifest-visible regular-file outputs
as durable `RunArtifactRef` records after a worker has committed `SUCCEEDED`.
SQLite remains canonical. Extraction is retryable, atomic, path-safe, and does
not affect Snakemake targets, Stage 51 decisions, or the terminal run outcome.

This PR does not add artifact HTTP routes, downloads, frontend browsing, QC
parsing, another adapter, object storage, checksums, or artifact-driven DAG
selection.

## Considered approaches

### Adapter candidates plus a platform safety/persistence service (selected)

The adapter uses its own manifest and catalog knowledge to return neutral
candidates containing a logical output type, a workspace-relative POSIX path,
and bounded logical metadata. `ArtifactExtractionService` independently checks
capability, build identity, lifecycle state, path containment, filesystem type,
symlinks, duplicates, and metadata before atomically replacing the complete
artifact set and recording an event.

This keeps ENCODE config, manifest vocabulary, and catalog conventions out of
platform code while retaining a second trust boundary around adapter output.

### Adapter returns `RunArtifactRef` directly (rejected)

This would let adapters choose run IDs, URIs, deterministic identities, local
path policy, and persistence metadata. It would collapse the platform trust
boundary and make cross-adapter safety inconsistent.

### Platform parses manifest TSV or scans `results/` (rejected)

Parsing the generated TSV would create a second manifest parser. Scanning with
glob patterns would index outputs that are not in the catalog and couple the
platform to ENCODE directory conventions. Both violate the existing adapter
boundary and the PR scope.

## Adapter contract

`WorkflowAdapter` gains `extract_artifacts(workspace)` and advertises it with
the `artifact_extract` capability. The return value is a `Result` containing an
immutable tuple of `ExtractedArtifactCandidate` values:

- `output_type`: stable adapter-owned logical type;
- `relative_path`: workspace-relative POSIX path;
- `mime_type`: optional deterministic media type;
- `metadata`: bounded JSON-compatible logical metadata.

The candidate deliberately contains no run ID, absolute path, URI, database
identity, timestamp, or SQLAlchemy type. Unsupported adapters omit the
capability and return an unsupported result if called directly.

## ENCODE extraction

The ENCODE adapter reads the materialized `config/config.yaml` and
`config/samples.tsv`, overrides manifest inputs to absolute in-workspace paths,
and calls `encode_pipeline.manifest.make.build_manifest_rows` directly. It
loads the existing artifact inventory and resolves every manifest output type
before status filtering. Unknown or ambiguous vocabulary fails closed.

Exact catalog output types are resolved directly. The two controlled dynamic
types `biorep<N>_final_bam` and `biorep<N>_final_bai` use a full-match decimal
placeholder; no general wildcard matcher is introduced.

Only `status == "present"` regular-file catalog outputs become candidates.
`missing` and `not_applicable` rows are ignored. The catalog's two known
directory aggregates, `macs3_peak` and `pooled_macs3_peak`, are recognized and
skipped because PR127 rejects directories and may not glob-expand them into a
new file vocabulary. A future file-level manifest vocabulary is required
before their contents can be indexed. If
`results/multiqc/result_manifest.tsv` exists, the adapter adds the catalog's
`result_manifest` entry as a project-level candidate.

Candidate metadata is restricted to catalog ID/scope/level plus the manifest's
sample ID, experiment ID, assay, target, genome, method, and QC flag. It never
contains the workspace path. No BAM, BigWig, peak, or report content is read.

## Platform validation and identities

`ArtifactExtractionService` receives the registry, `RunService`, workspace
root, and `WorkflowBuildIdentityProvider`. It performs these steps in order:

1. Read the canonical run and require `SUCCEEDED`.
2. Require the adapter capability.
3. Compare the persisted preflight build identity with a fresh capture.
4. Call the adapter and validate the complete candidate collection.
5. Atomically persist the complete reference set and success event.

The build identity source manifest adds
`docs/architecture/artifact-inventory.yaml`. This makes catalog changes part of
the same identity used by preflight, worker claim, execution, and later explicit
artifact retries.

Artifact IDs are `artifact-` plus the complete lowercase hex digest of a
framed SHA-256 over output type and relative path. The resulting 73-character
ID remains below the existing 128-character schema bound. URIs are opaque and
run-scoped: `run://runs/{run_id}/artifacts/{artifact_id}`. `produced_at` is the
canonical run `ended_at`, avoiding retry-time drift. Metadata contains only
controlled candidate fields, `relative_path`, `output_type`, and `size_bytes`.

## Path and filesystem policy

Every candidate must be a non-empty relative POSIX path under `results/`.
Absolute paths, Windows paths, backslashes, empty components, `.`, `..`, NUL,
and home expansion are rejected. Duplicate relative paths are rejected even if
their output types differ.

The service rejects a symlink workspace root and any symlink component using
`lstat`. The resolved target must remain under the run workspace and be a
regular file. Directories, FIFOs, sockets, block devices, and character devices
are rejected. Only stat metadata is read; files are never opened or hashed.

Metadata keys and scalar values are length-bounded. Nested containers, NUL,
absolute-path-like values, and values containing the absolute workspace string
are rejected. Events never include candidates, paths, exception strings,
configuration, commands, or environment details.

## Persistence and idempotency

The repository gains an atomic `replace_artifacts` operation. SQLAlchemy uses
one write transaction to re-read the run, require `SUCCEEDED`, delete the old
set, insert the complete sorted set, and append `artifacts_indexed` with only
`artifact_count`. InMemory performs the same operation under one lock.

If the complete set already matches and an `artifacts_indexed` event exists,
the operation is a no-op. Thus concurrent or repeated identical extraction
does not duplicate rows or events. A changed valid complete set replaces the
old set atomically. Insert/event failure rolls back the full transaction.
Validation and adapter failure occur before replace, so the previous complete
set remains intact.

`artifact_extraction_failed` is a separate best-effort event with status
`SUCCEEDED` and a controlled `reason_code` only. The run remains `SUCCEEDED`.
FAILED and CANCELLED runs are rejected before adapter access and receive no
artifact rows.

## Worker integration

Worker runtime composes and exposes `ArtifactExtractionService`. After
`LocalExecutionService.execute()` returns success and SQLite is canonically
`SUCCEEDED`, `_execute_claimed_run` invokes extraction. Expected extraction
failure is already recorded by the service and the RQ job returns normally.
Unexpected extraction exceptions are converted to the same safe failure event
on a best-effort basis and are never routed through execution-failure handling.

Timeout, cancellation, failed execution, identity mismatch before execution,
and work-horse failure do not call extraction.

## Testing

Tests cover the adapter protocol and unsupported capability; manifest status,
mixed assays, dynamic catalog mappings and directory aggregates; path,
symlink, special-file, metadata, duplicate and deterministic-identity attacks;
InMemory and SQLAlchemy atomic replacement, rollback, concurrency, idempotency
and reopen; worker success and failure isolation; lifecycle exclusion; and
catalog build identity drift. Deterministic temporary files replace scientific
datasets. Existing scientific Snakemake and Stage 49/50/51 tests continue to
protect output vocabulary and target selection.
