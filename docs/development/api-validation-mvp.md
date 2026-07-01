# API Validation MVP

This document describes how to start and use the validation-only FastAPI MVP for
the workflow platform. It covers the three endpoints that PR85 implemented and
how to interpret their responses.

The API is optional. Install it with the `api` extra:

```bash
pip install -e ".[api]"
```

## Startup

Start the server with the FastAPI app factory:

```bash
uvicorn encode_pipeline.api.main:create_app --factory
```

By default this binds to `http://127.0.0.1:8000`. The factory builds a FastAPI
app, wires the default workflow registry and validation service, and registers
the validation MVP routes under `/api/v1`.

## Endpoints

The validation MVP exposes three workflow-neutral endpoints:

- `GET /api/v1/workflows` — list registered workflows with metadata and
capabilities.
- `GET /api/v1/workflows/{workflow_id}/schema` — return adapter-owned schema
hints.
- `POST /api/v1/workflows/{workflow_id}/validate` — validate submitted workflow
inputs.

## curl examples

List workflows:

```bash
curl -s http://127.0.0.1:8000/api/v1/workflows
```

Get schema hints for the bundled ENCODE-style workflow:

```bash
curl -s \
  http://127.0.0.1:8000/api/v1/workflows/encode-style-chipseq-cuttag-atac-mnase/schema
```

Validate a workflow:

```bash
curl -s -X POST \
  http://127.0.0.1:8000/api/v1/workflows/encode-style-chipseq-cuttag-atac-mnase/validate \
  -H "Content-Type: application/json" \
  -d '{"config": {"use_control": false}, "samples": "samples.tsv", "options": {"strict_inputs": false}}'
```

Validation failures still return HTTP 200. The response body uses `ok: false`
and includes structured issues:

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

## Response envelope

All endpoints return endpoint-specific envelopes. Every response contains:

- `ok`: boolean success marker for the API call.
- `issues`: list of structured issue objects.

Validation responses also include `workflow_id` and `value`. Schema responses
include `workflow_id` and `schema_hints`. Workflow list responses include
`workflows`.

The `value` returned by validation is adapter-private. Do not depend on its
internal structure until the platform introduces a shared success model.

## HTTP status semantics

- `200 OK` — successful list, schema lookup, or syntactically valid validation
request. Validation failures are returned as `200` with `ok: false`.
- `400 Bad Request` — malformed JSON or a request body that cannot be translated
into `WorkflowInputs`.
- `404 Not Found` — the requested `workflow_id` is not registered.
- `409 Conflict` — the workflow exists but does not support validation.
- `500 Internal Server Error` — unexpected server error. No tracebacks or
internal paths are exposed in the response.

## Non-goals

This MVP intentionally does not include:

- Workflow execution, run lifecycle, or job management.
- Logs, progress, or SSE streaming.
- Artifacts, QC summaries, or provenance endpoints.
- Authentication or authorization.
- Frontend code or server-side rendering.
- Agent tool backend or dynamic plugin loading.

For the formal API contract, see
[workflow-platform-api-contract.md](workflow-platform-api-contract.md).
