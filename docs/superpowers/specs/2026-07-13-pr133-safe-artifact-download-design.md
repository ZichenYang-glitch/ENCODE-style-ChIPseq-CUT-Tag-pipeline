# PR133 Safe Artifact Download Design

## Scope and intent

PR133 adds a run-scoped download endpoint for artifact references already
persisted by PR127 and exposed by PR128. A request supplies only `run_id` and
the opaque `artifact_id`. The platform resolves the canonical reference from
SQLite, validates it again, opens the referenced file through a descriptor-only
workspace boundary, and streams a bounded response.

The route never accepts a path and never opens the workspace itself. The
feature does not add Range requests, checksums, object storage, authentication,
artifact extraction changes, QC behavior, or scientific workflow changes.

## Alternatives considered

1. Return `FileResponse(workspace / relative_path)`. This delegates efficient
   streaming but passes an untrusted mutable pathname to a later open, leaving
   traversal, symlink, and time-of-check/time-of-use boundaries implicit.
2. Read the complete file into memory after validating it. This simplifies
   response errors but is unsuitable for BAM/BigWig-sized artifacts and violates
   the bounded streaming requirement.
3. Resolve the database reference in a platform service, open every path
   component with `dir_fd`, `O_NOFOLLOW`, and `O_CLOEXEC`, retain the verified
   descriptor chain, and stream only the persisted byte count. This is selected.
   It avoids pathname reopening, bounds memory, detects path-entry replacement,
   and gives descriptor ownership an explicit lifecycle.

## Public API contract

The new operation is:

```text
GET /api/v1/runs/{run_id}/artifacts/{artifact_id}/download
operation_id: downloadRunArtifact
```

Successful responses are binary streams with:

- `Content-Length` equal to the persisted and verified `size_bytes`;
- a validated persisted MIME type, or `application/octet-stream` when absent;
- `Content-Disposition: attachment` with an ASCII-safe fallback filename and
  an RFC 5987 UTF-8 filename;
- `X-Content-Type-Options: nosniff` and `Cache-Control: private, no-store`.

The endpoint ignores no hidden path input because none exists. Range is not
implemented; a Range request receives the same full `200` response without
`Content-Range` or partial semantics.

JSON failures use `RunArtifactDownloadErrorResponse` with `ok=false`, the
requested run/artifact identifiers, and disclosure-safe `IssueResponse` values:

- `404 RUN_NOT_FOUND` for an unknown run;
- `404 RUN_ARTIFACT_NOT_FOUND` for both an unknown artifact and an artifact
  belonging to another run;
- `409 RUN_ARTIFACT_DOWNLOAD_CONFLICT` when the persisted reference exists but
  its file is missing, moved, symlinked, non-regular, size-mismatched, or changed;
- `500 RUN_ARTIFACT_DOWNLOAD_DATA_INVALID` for a corrupted persisted reference;
- `500 RUN_ARTIFACT_DOWNLOAD_UNAVAILABLE` for an infrastructure/runtime failure.

Only a stable allowlisted `reason_code` may appear in conflict/unavailable
context. Absolute paths, database URLs, OS messages, exception text, commands,
environment values, source bytes, and `technical_message` never appear.

## Service boundary

`ArtifactDownloadService` is constructed with `RunService` and the absolute
workspace root. `prepare(run_id, artifact_id)` performs the following before a
response is returned:

1. Load the run and then the artifact through run-scoped `RunService` methods.
2. Validate artifact/run identity, opaque URI, file type, name, MIME type,
   `results/`-relative canonical POSIX path, and strict non-negative persisted
   size. A bool is never accepted as an integer size.
3. Validate the run ID as one safe workspace component.
4. Open `/`, every workspace-root component, the run directory, every relative
   parent, and the final file using `dir_fd`. Directories use
   `O_DIRECTORY|O_NOFOLLOW|O_CLOEXEC`; the final node additionally uses
   `O_NONBLOCK` so a FIFO cannot hang before it is rejected. Cleanup ownership
   is recorded immediately after every successful `open`, before `fstat`,
   fingerprinting, or node construction can fail.
5. Require every parent to be a directory and the final descriptor to be a
   regular file whose `st_size` equals persisted `size_bytes`.
6. Compare every descriptor fingerprint to its current parent-directory entry
   before returning the prepared plan.

The service returns a `Result[ArtifactDownloadPlan]`, not a path. The plan
exposes controlled response metadata and a synchronous byte iterator. Route
code cannot convert a relative path into a local file read.

## Descriptor lifecycle and streaming invariants

The plan owns the complete open descriptor chain. Its iterator:

- verifies each descriptor and parent entry before every read;
- reads at most 64 KiB and never more than the persisted remaining byte count;
- fails if EOF arrives early, if an additional byte exists after the expected
  length, or if device/inode/mode/size/mtime changes;
- rechecks the chain after the final byte;
- closes every descriptor in `finally`, even on read failure or generator close.

Before ownership transfers to the plan, preparation keeps an independent list
of every opened descriptor. Any exception between `open` and plan construction
closes that full list in reverse order.

Directory fingerprints intentionally compare only device, inode, and mode:
directory size and mtime change when an unrelated sibling is created and are
not path-entry identity. The final regular file additionally retains its
indexed size and mtime checks.

The `StreamingResponse` installs both a Starlette `BackgroundTask(plan.close)`
and a narrow response `finally` boundary. `close` is idempotent, so normal
completion, client disconnect, response transmission failure, response
cleanup, or an iterator exception all converge on one descriptor release.

Opening the final descriptor pins the bytes to one inode, while the retained
parent chain detects a path entry replaced after preparation. Mutation detected
after some chunks have been sent aborts the stream; HTTP headers cannot be
rewritten into a JSON Issue after streaming begins. This is an inherent HTTP
streaming boundary, not a reason to buffer large artifacts.

## Path and header policy

The relative path must be a canonical non-empty POSIX string under `results/`.
Absolute paths, Windows drive/UNC forms, `~`, backslashes, NUL, empty/dot/dot-dot
components, non-canonical repeated separators, directories, symlink components,
FIFO, socket, and device nodes fail closed.

The download filename must match the persisted basename and contain no
separator or control character. The header fallback contains only ASCII
letters, digits, dot, underscore, and hyphen and is length-bounded; the UTF-8
value is percent-encoded. MIME has the same parameter-free token grammar used
by artifact extraction. These values can never inject another response header.

## FastAPI and OpenAPI composition

The endpoint is a synchronous `def`, so FastAPI runs blocking SQLite and
descriptor preparation in its maintained threadpool. The synchronous stream
iterator is likewise consumed through Starlette's threadpool boundary.

`create_app` constructs and stores one `ArtifactDownloadService`; a dependency
returns it to the thin route. Unexpected pre-stream failures use an
operation-specific generic exception branch so even defensive 500 responses
retain the declared download error envelope.

OpenAPI declares the success schema as `type: string, format: binary` and the
stable JSON 404/409/500 models. Orval mechanically generates
`downloadRunArtifact`. Its operation-specific mutator is a `blobFetcher` that
shares the existing URL, credentials, and redacted error normalization with the
JSON `fetcher`, but decodes successful responses with `Response.blob()`. This
is one API transport boundary with two response decoders, not a handwritten
endpoint or second client. It also handles artifacts whose legitimate MIME is
`application/json` without treating their content as an API envelope.

## Frontend interaction

`ArtifactBrowser` owns a TanStack Query mutation calling the generated
`downloadRunArtifact` operation. The selected artifact inspector receives only
download state and callbacks. On success, it creates a temporary object URL,
uses a controlled filename, activates a temporary anchor, removes it, and
revokes the URL. Selection identity is retained in mutation variables so a late
response cannot be mislabeled as a different artifact.

The inspector uses the existing `Button` and lucide `Download` icon with a
tooltip and accessible name. Pending disables only the download action and says
“Downloading…”. A failure shows a local redacted status and Retry without
removing the confirmed artifact detail. Success is announced through a polite
live region. There is no preview, Range UI, or fabricated availability state.

## Tests and residual risks

Service tests cover run/artifact isolation, path grammar, every symlink level,
directory/FIFO/device rejection, missing files, larger/smaller size mismatch,
path/inode replacement, in-place mutation during iteration, read errors,
bounded chunking, and idempotent descriptor cleanup. API tests cover binary
bytes, safe headers, unknown/cross-run identity, disclosure-safe 409/500,
SQLite reopen, Range non-support, no lifecycle mutation, and stable OpenAPI.

Fetcher and UI tests cover blob decoding even for JSON MIME, redacted API
errors, loading/success/failure/retry, filename control, object-URL cleanup, and
absence of raw error text. The real Playwright success path downloads the
persisted `result_manifest.tsv`, verifies its suggested filename and exact
bytes, and still exercises mobile cancellation cleanup.

Without a checksum, an attacker able to rewrite bytes in place, preserve the
same inode/size/mtime, and race between checks can evade metadata-only change
detection. Immutable artifact storage or checksums would narrow that risk but
are explicitly deferred. A client disconnect or mutation after headers are sent
produces a truncated transport rather than a JSON Issue; descriptors still
close. The browser's generated fetch operation materializes a `Blob` before
handing it to the native download action, so very large artifacts can pressure
client memory even though the server remains chunk-bounded; a future native or
stream-to-disk client boundary can address this without weakening the endpoint.
Authentication and multi-user authorization also remain later work.
