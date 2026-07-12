# PR130 Workflow-Neutral QC Summary Indexing Design

## Objective

Index a small, trusted set of numeric QC metrics after a canonical successful
run. The ENCODE adapter owns source selection and scientific field mapping;
the platform owns safe file access, type and identity validation, atomic
persistence, idempotency, events, and error redaction. SQLite remains the only
truth source for lifecycle and indexed QC metadata.

This PR adds no HTTP routes, OpenAPI models, frontend, QC dashboard, threshold
editor, Snakemake target, manifest behavior, Stage 51 decision, new adapter, or
runtime queue.

## Evidence that a stable source exists

The repository already produces three machine-readable, manifest-visible
summary types with fixed writer contracts and scientific regression tests:

- `qc_summary`: one 37-column TSV row per peak-centric sample, written by
  `scripts/assemble_qc_summary.py`;
- `mnase_qc_summary`: one fixed-header TSV row per MNase sample, written by
  `scripts/mnase_qc_summary.py`;
- `pooled_qc_summary`: one fixed-header TSV row per eligible pooled experiment,
  written by `scripts/pooled_qc_summary.py`.

All three are catalogued in `docs/architecture/artifact-inventory.yaml`, are
represented by the existing manifest, and are indexed as `RunArtifactRef`
records by PR127. Missing or not-applicable manifest rows do not produce
artifact references, which gives an explicit optional-source absence semantic.

`stage3_qc_summary.tsv` duplicates per-sample `qc_summary` rows and is not read.
MultiQC HTML is never parsed. MultiQC custom data and other component TSVs are
not scanned in v1 because they do not add a single, uniformly manifest-visible
contract across all four assays.

## Considered approaches

### Persisted artifact sources plus a platform safe-document boundary (selected)

The worker first completes existing artifact indexing. The QC service receives
that complete, successful artifact set. The adapter declares the exact output
types it can parse. The platform selects only matching persisted artifact
references, validates and reads each small regular file inside the run
workspace, then passes immutable in-memory source documents to the adapter.
The adapter parses the trusted TSV contracts and returns neutral metric
candidates. The platform validates the complete candidate set and atomically
replaces durable QC rows.

This approach makes optional absence explicit, avoids directory scanning,
retains a durable source-artifact reference, and prevents an adapter parser
from opening an unchecked source path.

### Let the adapter open workspace files directly (rejected)

A single `extract_qc_metrics(workspace)` call would be simpler, but the adapter
would read a path before the platform could enforce containment, symlink, file
type, and size policy. Post-parse validation cannot undo an unsafe read.

### Parse MultiQC data or HTML as the canonical source (rejected)

HTML parsing is explicitly unsafe and unstable. MultiQC data availability and
module vocabulary depend on configuration and tool versions, while the three
summary TSV writers already have repository-owned contracts. Adding a second
normalization layer would weaken rather than strengthen the trust boundary.

## Adapter-neutral contracts

The optional `QcSummaryExtractingAdapter` protocol is separate from the base
`WorkflowAdapter`, so adapters without the `qc_summary_extract` capability do
not need placeholder implementations. It exposes two operations:

1. `qc_source_output_types()` returns a bounded tuple of adapter-owned artifact
   output types. ENCODE returns exactly `qc_summary`, `mnase_qc_summary`, and
   `pooled_qc_summary`.
2. `extract_qc_metrics(inputs, sources)` parses platform-vetted
   `QcSourceDocument` values and returns an immutable tuple of
   `ExtractedQcMetricCandidate` values.

`QcSourceDocument` contains only the persisted source artifact ID, logical
output type, workspace-relative POSIX path, bounded logical artifact metadata,
and bounded file bytes. It contains no absolute path, file handle, SQLAlchemy
row, run service, registry, or environment state.

`ExtractedQcMetricCandidate` contains:

- `metric_key` and `display_name`;
- a finite `Decimal` value with at most 26 integer and 12 fractional digits;
- a controlled unit (`count`, `fraction`, or `ratio` in v1);
- scope (`run`, `sample`, or `experiment`);
- optional `sample_id`, `experiment_id`, and `assay` constrained by scope;
- optional `qc_flag` (`pass`, `warning`, or `fail`);
- `source_artifact_id`.

There is no arbitrary metric metadata JSON.

Platform validation limits `metric_key` to a lowercase dotted token of at most
128 characters, `display_name` to 1..255 printable characters, identifiers and
assay to 1..255 characters from the canonical ENCODE-safe
`[A-Za-z0-9_.-]` vocabulary, and source/output IDs to the existing artifact
identifier grammar. The shared workflow-neutral token validator permits `-`,
`.`, or `_` in the first position, matching scientific sample and experiment
IDs, while rejecting empty strings, exact `.`/`..`, controls, slash,
backslash, and overlong values. Adapter source types are unique, bounded safe
output tokens.

## ENCODE v1 metric catalog

The mapping is code-owned in `encode_pipeline.adapters.encode_qc`; changes are
therefore included in the existing source-tree build identity.

| Source | Scope | Metric key | Display name | Unit | Validation |
| --- | --- | --- | --- | --- | --- |
| `qc_summary` | sample | `sequencing.total_reads` | Total reads | count | integer, >= 0 |
| `qc_summary` | sample | `peaks.frip` | Fraction of reads in peaks | fraction | 0..1 |
| `qc_summary` | sample | `peaks.count` | Peak count | count | integer, >= 0 |
| `qc_summary` | sample | `library.percent_duplication` | Percent duplication | fraction | 0..1 |
| `qc_summary` | sample | `library.estimated_size` | Estimated library size | count | integer, >= 0 |
| `qc_summary` | sample | `library.nrf` | Non-redundant fraction | fraction | 0..1 |
| `qc_summary` | sample | `library.pbc1` | PCR bottleneck coefficient 1 | fraction | 0..1 |
| `qc_summary` | sample | `library.pbc2` | PCR bottleneck coefficient 2 | ratio | >= 0 |
| `mnase_qc_summary` | sample | `mnase.subnucleosomal_reads` | Sub-nucleosomal reads | count | integer, >= 0 |
| `mnase_qc_summary` | sample | `mnase.mononucleosomal_reads` | Mono-nucleosomal reads | count | integer, >= 0 |
| `mnase_qc_summary` | sample | `mnase.dinucleosomal_reads` | Di-nucleosomal reads | count | integer, >= 0 |
| `pooled_qc_summary` | experiment | `peaks.pooled_count` | Pooled peak count | count | integer, >= 0 |
| `pooled_qc_summary` | experiment | `replicates.biological_count` | Biological replicates | count | positive integer |

`NA` means that an optional metric is absent and is skipped. An unknown token,
empty required identity, malformed number, NaN, Infinity, out-of-range value,
wrong header, extra data row, duplicate source, or duplicate semantic metric
fails the complete extraction. The parser checks row scope IDs and assay
against the persisted artifact metadata.

The pooled peak-count candidate carries `pass` when `peak_mode_status == ok`,
`warning` when it is `mismatch`, and no flag when it is `unknown`. Other v1
metrics have no flag. The indexer does not invent universal assay thresholds;
`docs/qc-interpretation.md` explicitly describes QC as contextual and
descriptive.

## Durable read model

`RunQcMetric` is a workflow-neutral immutable domain value with:

- `metric_id`, `run_id`, `metric_key`, `display_name`, `value`, and `unit`;
- `scope`, `sample_id`, `experiment_id`, and `assay`;
- `qc_flag`;
- `source_artifact_id`;
- `produced_at`.

`metric_id` is `qcmetric-` plus a complete SHA-256 digest over a framed tuple of
`metric_key`, scope, sample ID, and experiment ID. It is deterministic across
retries and independent of local paths or database row IDs. The domain layer
owns the only implementation of this algorithm. Both the indexing service and
repository call it, and the repository recomputes the ID from every metric's
semantic coordinates before any write. A well-formed but semantically
mismatched digest and duplicate semantic coordinates under different digests
therefore fail atomically before persistence.

Alembic revision `20260712_05` creates `run_qc_metrics` with a run foreign key,
unique `(run_id, metric_id)`, stable run index, canonical decimal text, scope
and flag checks, and scope/identifier consistency checks. SQLite `NUMERIC`
affinity is deliberately not used: it converts large decimals through binary
floating point and cannot preserve values above `2^53`. The repository stores
only a canonical non-exponent decimal grammar with no leading zeros, at most 26
integer digits and 12 fractional digits, and reconstructs the domain `Decimal`
on read. Tests cover values above `2^53`, the maximum boundary, reopen, and
idempotent equality.

The source artifact ID is a validated logical reference rather than a foreign
key: artifact replacement deletes and reinserts a complete set, and a database
FK would make an otherwise idempotent artifact retry impossible. The QC
repository verifies the complete expected artifact set and every source in the
same write transaction that replaces QC rows.

No absolute source path, workspace root, database URL, command, environment,
exception string, or arbitrary JSON is stored.

## Platform indexing flow

`QcSummaryIndexingService.index(run_id, artifacts)` performs:

1. Read the canonical run and require `SUCCEEDED` with `ended_at`.
2. Recompute and match the durable workflow build identity.
3. Require the registered adapter capability and optional QC protocol.
4. Validate that `artifacts` exactly matches the current complete persisted
   artifact set produced by successful artifact extraction.
5. Ask the adapter for its exact source output types and select only matching
   persisted artifacts.
6. Validate each relative path under `results/`; reject absolute paths,
   backslashes, NUL, empty/dot/traversal components, duplicate paths, directory
   components that are symlinks, and non-regular sources.
7. Open from the filesystem root with directory descriptors. Every absolute
   workspace component and relative source parent is traversed with `dir_fd`,
   `O_NOFOLLOW`, `O_DIRECTORY`, and `O_CLOEXEC`; the final source is opened with
   `O_NOFOLLOW`, `O_NONBLOCK`, and `O_CLOEXEC`. `fstat` must identify a regular
   file no larger than 1 MiB. Persisted `size_bytes` must be a non-boolean,
   non-negative integer within that limit; the opened descriptor's `st_size`
   must equal it. Bytes are read only from that descriptor, bounded to expected
   size + 1, final content length must equal `size_bytes`, and pre/post `fstat`
   device, inode, mode, size, and mtime must match.
   The implementation never calls `Path.open`, `open`, or `read_bytes` on an
   unchecked source. A race canary replaces a path component during traversal
   and proves that no outside sentinel bytes reach the adapter.
8. Pass immutable source documents to the adapter and validate every returned
   candidate, including source membership, identifiers, decimal precision,
   units, flags, and duplicate semantic identities.
9. Construct sorted `RunQcMetric` values with `produced_at = run.ended_at` and
   atomically replace the complete set. The repository receives the complete
   expected `RunArtifactRef` set, re-reads current artifacts under the same
   SQLite `BEGIN IMMEDIATE` or InMemory lock, compares the sets independent of
   order, and verifies every metric source before changing QC rows.

An empty successful artifact set or a set containing no declared QC source is
a successful empty QC index. A source artifact selected by the adapter that is
missing, unreadable, malformed, unsafe, or no longer identical is a failure,
not an empty result.

## Persistence, events, and retry semantics

`replace_qc_metrics` is added to `RunRepository`, with matching InMemory and
SQLAlchemy behavior. SQLAlchemy re-reads the run, requires `SUCCEEDED`,
re-reads and matches the complete expected artifact set, validates every source
and the complete metric replacement, deletes old rows, inserts new rows, and
adds `qc_metrics_indexed` with only `metric_count` in one write transaction.
The in-memory repository performs the same operation under one lock.

If the current complete set is identical and the latest QC outcome is
`qc_metrics_indexed`, retry is a no-op. A changed valid complete set replaces
the old set atomically. Validation, parsing, insert, or event failure leaves
the previous complete set intact.

Expected and defensive failures append `qc_metrics_indexing_failed` only while
the run remains `SUCCEEDED`; context contains one controlled `reason_code` and
nothing else. The run remains `SUCCEEDED`, and the RQ job is not converted into
an execution failure. Repeating the same latest failure reason is an atomic
no-op, avoiding retry event spam. FAILED and CANCELLED runs are rejected before
source selection or file access and receive no metrics.

Any later non-equivalent `replace_artifacts` transaction atomically deletes the
old QC rows and appends `qc_metrics_invalidated` after `artifacts_indexed`. This
also invalidates a previously confirmed empty QC index. An identical artifact
retry remains a no-op and preserves current QC. Therefore old metrics can never
remain current after their complete source artifact generation changes.

Artifact-set equivalence compares the complete `RunArtifactRef` set by identity
and content, independent of tuple order. `qc_metrics_invalidated` participates
in the latest QC outcome set. Artifact replacement writes invalidation only
when a prior QC outcome or QC row exists, and does not repeat it when the latest
outcome is already invalidated. The first artifact index for a run with no QC
history therefore never claims that QC was invalidated.

## Worker integration and ordering

Worker runtime composes `QcSummaryIndexingService`. After
`LocalExecutionService` has committed `SUCCEEDED`, `_execute_claimed_run` runs
artifact extraction first and passes its successful complete result to QC
indexing. Artifact extraction failure records a safe QC indexing failure rather
than treating an unknown source set as an empty index. QC indexing failure is
contained exactly like artifact indexing failure and does not change the RQ
execution result or canonical run terminal state.

Timeout, cancellation, failed execution, build mismatch before execution, and
work-horse failure never index QC metrics.

## Build identity

The current `sha256-tree-v1` identity already includes all Python files under
`src/encode_pipeline`, all top-level scripts, the workflow tree, and the
artifact inventory. The new parser/catalog module and the three existing TSV
writer contracts are therefore fingerprinted without adding a second version
mechanism. Tests explicitly prove a parser mapping change alters the digest.

## Test strategy

Tests cover the optional adapter protocol and unsupported capability; the
three exact TSV contracts, mixed assays, `NA`, empty source sets, pooled flag
mapping, malformed rows, duplicate metrics, units, decimal precision, NaN and
Infinity; absolute/traversal/symlink/directory/FIFO/device and oversized-file
attacks plus a descriptor traversal race canary; deterministic IDs and absence
of absolute paths; InMemory and
SQLAlchemy atomic replacement, rollback, idempotency, concurrency, and reopen;
persisted-size larger/smaller mismatches and mutation during descriptor reads;
semantic-ID recomputation and canonical leading-punctuation identifiers;
QC-versus-artifact replacement races and artifact-driven invalidation; Alembic
upgrade from `20260712_04`; worker success, empty, failure isolation,
and lifecycle exclusion; and build-identity drift.

The existing full Python suite and real Redis/RQ/Snakemake durable gates remain
required. OpenAPI, Orval, React, package-lock, Snakemake targets, Stage 51, and
manifest output are expected to have zero diff.

## Deferred scope and residual risk

- Read-only QC API and UI are later PRs.
- Metric threshold interpretation remains adapter/catalog work; v1 does not
  infer universal pass/fail from contextual scientific guidance.
- Component-level FastQC, cross-correlation, TSS, and MultiQC module metrics
  need an explicit stable adapter mapping before adoption.
- Source documents are capped at 1 MiB. A future legitimate larger summary
  requires an intentional policy change rather than silent acceptance.
- Without a content checksum, replacement by different bytes with the same
  size and filesystem identity evidence outside the descriptor's observed
  lifetime cannot be distinguished. This PR deliberately performs safe stat
  and bounded reads rather than hashing potentially large workflow artifacts.
- There is no heartbeat/reconciler, immutable workflow bundle, object storage,
  checksum, authentication, multi-tenancy, HPC, or second adapter work here.
