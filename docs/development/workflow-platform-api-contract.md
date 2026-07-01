# Workflow Platform API Contract

This document defines the validation MVP API contract for the workflow
platform. It is a contract document only. It does not implement FastAPI routes,
frontend code, run lifecycle, workflow execution, or agent-facing tools.

## Context

PR77 through PR83 established the current workflow-platform foundation:

- PR77 defined the workflow platform and agent-ready adapter roadmap.
- PR78 added platform `Result` and `Issue` primitives.
- PR79 added the `WorkflowAdapter` contract.
- PR80 added the ENCODE-style adapter wrapper with validation-only capability.
- PR81 added `WorkflowRegistry`.
- PR82 added `ValidationService`.
- PR83 added default composition factories for the bundled ENCODE-style
  adapter and validation service.

The backend validation chain is now default-startable through
`create_default_validation_service()`. This document defines the HTTP contract
that future FastAPI routes and frontend validation UX should follow.

## Goals

- Provide a stable validation MVP API contract before FastAPI and frontend
  implementation.
- Map existing platform models to JSON request and response shapes:
  `WorkflowMetadata`, `WorkflowCapabilities`, `WorkflowSchema`,
  `WorkflowInputs`, `Result`, and `Issue`.
- Keep future route handlers thin over service-layer objects.
- Keep the API workflow-neutral. Workflow-specific behavior remains behind
  adapters.

## Non-Goals

- No FastAPI implementation.
- No frontend code.
- No authentication or authorization.
- No run lifecycle.
- No execution, logs, or SSE streaming.
- No artifacts, QC, or provenance endpoints.
- No plugin discovery or dynamic loading.
- No agent tool implementation.

## Endpoints

The validation MVP exposes three workflow-neutral endpoints.

### `GET /api/v1/workflows`

Lists available workflows from the default workflow registry.

Successful response:

```json
{
  "ok": true,
  "workflows": [
    {
      "metadata": {
        "workflow_id": "encode-style-chipseq-cuttag-atac-mnase",
        "name": "ENCODE-style ChIP-seq/CUT&Tag/ATAC/MNase",
        "version": "0.3.0",
        "description": "Adapter wrapper for the current ENCODE-style workflow.",
        "engines": ["snakemake"],
        "tags": ["chipseq", "cuttag", "atac", "mnase", "encode-style"]
      },
      "capabilities": {
        "supports": ["validation"]
      }
    }
  ],
  "issues": []
}
```

The `metadata` object maps from `WorkflowMetadata.to_dict()`. The
`capabilities` object maps from `WorkflowCapabilities.to_dict()`.

### `GET /api/v1/workflows/{workflow_id}/schema`

Returns adapter-owned schema hints for one workflow.

Successful response:

```json
{
  "ok": true,
  "workflow_id": "encode-style-chipseq-cuttag-atac-mnase",
  "schema": {
    "config_schema": {
      "schema_kind": "hints",
      "required": ["samples"],
      "properties": {
        "samples": {
          "type": "string",
          "description": "Path to the sample sheet TSV."
        }
      }
    },
    "sample_schema": {
      "schema_kind": "hints",
      "format": "tsv",
      "required_columns": ["sample", "fastq_1"]
    },
    "option_schema": {
      "schema_kind": "hints",
      "properties": {
        "strict_inputs": {
          "type": "boolean",
          "default": false
        }
      }
    }
  },
  "issues": []
}
```

The `schema` object maps from `WorkflowSchema.to_dict()`. Schema payloads are
adapter-owned discovery hints unless an adapter explicitly documents stronger
JSON Schema compatibility. The ENCODE-style schema hints may mention
`config.samples` because they reflect the adapter's internal configuration
shape. API clients do not need to duplicate the sample path inside `config`
when they provide the top-level `samples` field.

### `POST /api/v1/workflows/{workflow_id}/validate`

Validates submitted workflow inputs through `ValidationService`.

Request body:

```json
{
  "config": {
    "use_control": false
  },
  "samples": "samples.tsv",
  "options": {
    "strict_inputs": false
  }
}
```

Successful validation response:

```json
{
  "ok": true,
  "workflow_id": "encode-style-chipseq-cuttag-atac-mnase",
  "value": {
    "config": {},
    "samples": []
  },
  "issues": []
}
```

Failed workflow validation response:

```json
{
  "ok": false,
  "workflow_id": "encode-style-chipseq-cuttag-atac-mnase",
  "value": null,
  "issues": [
    {
      "code": "ENCODE_CONFIG_INVALID",
      "message": "Config is invalid.",
      "severity": "error",
      "path": "config",
      "source": "config",
      "technical_message": "Config is invalid.",
      "hint": null,
      "context": {}
    }
  ]
}
```

The request body maps to `WorkflowInputs`:

- `config`: object. Required.
- `samples`: string path, inline list of row objects, or `null`.
- `options`: object. Optional; routes should treat omission as `{}`.

Top-level `samples` is the workflow-platform input field. For the current
ENCODE-style adapter, a string top-level sample path is copied into
`config["samples"]` before validation. Existing ENCODE-style `config.samples`
remains accepted for compatibility. If both top-level `samples` and
`config.samples` are present, top-level `samples` wins. Inline sample rows are
part of the platform request shape, but they are currently unsupported by the
ENCODE-style adapter and should return `ENCODE_ADAPTER_UNSUPPORTED` until
adapter support is added.

The validation response maps from `Result` and `Issue`, but it is not raw
`Result.to_dict()` alone. It keeps `ok`, `value`, and `issues`, and adds
endpoint context through `workflow_id`.

## Response Envelopes

All API responses use endpoint-specific envelopes. Every response includes:

- `ok`: boolean success marker.
- `issues`: list of structured `Issue` objects.

Each endpoint also includes endpoint context:

- Workflow list responses include `workflows`.
- Schema responses include `workflow_id` and `schema`.
- Validation responses include `workflow_id` and `value`.

Route handlers should not expose Python object representations. They should
serialize platform model objects through their `to_dict()` methods or equivalent
JSON-ready conversions.

## HTTP Status Mapping

HTTP status codes describe transport and API-level outcome. The response body
still uses the same envelope and issue shape.

### `200 OK`

Use `200` for successful workflow listing, successful schema lookup, and
syntactically valid validation requests.

Validation can return `200` with `ok: false` when the submitted workflow inputs
are invalid. In that case, the request was well-formed at the API level and the
workflow adapter returned structured validation issues.

```json
{
  "ok": false,
  "workflow_id": "encode-style-chipseq-cuttag-atac-mnase",
  "value": null,
  "issues": [
    {
      "code": "ENCODE_SAMPLES_INVALID",
      "message": "Sample sheet is invalid.",
      "severity": "error",
      "path": "samples",
      "source": "samples",
      "technical_message": "Sample sheet is invalid.",
      "hint": null,
      "context": {}
    }
  ]
}
```

### `404 Not Found`

Use `404` when `workflow_id` is not registered.

```json
{
  "ok": false,
  "workflow_id": "missing-workflow",
  "value": null,
  "issues": [
    {
      "code": "WORKFLOW_NOT_FOUND",
      "message": "Workflow was not found.",
      "severity": "error",
      "path": "workflow_id",
      "source": "registry",
      "technical_message": null,
      "hint": null,
      "context": {
        "workflow_id": "missing-workflow"
      }
    }
  ]
}
```

For the schema endpoint, the same issue shape is used with `schema: null`
instead of `value: null`.

### `409 Conflict`

Use `409` when the workflow exists but does not support the requested
capability.

```json
{
  "ok": false,
  "workflow_id": "some-workflow",
  "value": null,
  "issues": [
    {
      "code": "WORKFLOW_CAPABILITY_UNSUPPORTED",
      "message": "Workflow does not support validation.",
      "severity": "error",
      "path": "workflow.capabilities",
      "source": "registry",
      "technical_message": null,
      "hint": null,
      "context": {
        "workflow_id": "some-workflow",
        "capability": "validation"
      }
    }
  ]
}
```

### `400 Bad Request`

Use `400` when the JSON request is malformed or cannot be translated into
`WorkflowInputs`.

```json
{
  "ok": false,
  "workflow_id": "encode-style-chipseq-cuttag-atac-mnase",
  "value": null,
  "issues": [
    {
      "code": "API_REQUEST_INVALID",
      "message": "Request body is invalid.",
      "severity": "error",
      "path": "body",
      "source": "api",
      "technical_message": null,
      "hint": "Submit an object with config, samples, and options fields.",
      "context": {}
    }
  ]
}
```

### `500 Internal Server Error`

Use `500` for unexpected server errors.

```json
{
  "ok": false,
  "workflow_id": "encode-style-chipseq-cuttag-atac-mnase",
  "value": null,
  "issues": [
    {
      "code": "INTERNAL_SERVER_ERROR",
      "message": "An internal server error occurred.",
      "severity": "error",
      "path": null,
      "source": "runtime",
      "technical_message": null,
      "hint": null,
      "context": {}
    }
  ]
}
```

User-facing messages must not include tracebacks, local filesystem internals, or
secret values. Detailed diagnostics belong in server logs.

## Issue Semantics For Frontend UX

Frontend code should treat `Issue` as the primary error and warning model.

- `message` is the primary user-facing text.
- `technical_message` is debug detail for expandable or developer-oriented
  views.
- `hint` is optional remediation help.
- `path` should drive field-level UI, such as highlighting `config.samples` or
  `options.strict_inputs`.
- `source` groups issues by component, such as `config`, `samples`, `adapter`,
  `registry`, `api`, or `runtime`.
- `context` is structured support data. It is not prose and must not be treated
  as provenance.
- `severity` controls display priority. Error issues block validation success;
  warning and info issues can be shown without marking the result failed.

Frontend screens should preserve issue codes for filtering, support links,
analytics, and future localization. They should display user-facing messages
first and keep technical details secondary.

## Agent Boundaries

For the validation MVP, agents may:

- List workflows.
- Inspect workflow schemas.
- Submit validation requests.
- Explain returned issues.

Agents must not:

- Run, kill, delete, or mutate workflows.
- Bypass services or adapters.
- Treat validation summaries as provenance.
- Infer execution state from validation responses.

Agent-facing behavior should remain read-only for this API contract. Any future
state-changing operation needs an explicit confirmation boundary and a separate
design.

## Implementation Notes

- FastAPI route construction should start from
  `create_default_validation_service()`.
- Workflow list and schema routes may use
  `create_default_workflow_registry()` or a shared application composition
  object when an application layer exists.
- Route handlers should translate HTTP payloads into `WorkflowInputs`.
- Route handlers should translate service `Result` values into the envelopes
  documented here.
- Route handlers must not import validator, sample-loader, or workflow engine
  internals directly.
- The API contract must remain workflow-neutral. Workflow-specific schema,
  validation, and issue details come from adapters.
- The `value` returned by validation is adapter-private. API consumers should
  not depend on ENCODE-specific success payload internals until a separate
  platform model is introduced.

## Verification Checklist

For the docs-only API contract PR:

- `git diff --name-only origin/main..HEAD` shows only
  `docs/development/workflow-platform-api-contract.md`.
- `git diff --check origin/main..HEAD` is clean.
- Placeholder scans over this document have no matches for unfinished markers.
- Commit messages contain no attribution, tool, or vendor trailers.
- `docs/superpowers/`, `research/`, and the untracked external architecture
  report are not committed.
- No Python code, FastAPI code, frontend code, CI files, lockfiles, Snakefile,
  or workflow rule files are changed.
