# Agent Layer Boundary Design

> **Status:** design review<br>
> **Scope:** read-only validation assistant for the workflow platform<br>
> **Date:** 2026-07-02

This document defines the first agent layer for the workflow validation platform. The agent is intentionally narrow: a **read-only validation assistant bound to one workflow**. Execution, mutation, filesystem access, command generation, provenance claims, and direct adapter access are all out of scope until the relevant platform layers exist.

The design distills the multi-perspective Agency Agent review and the external `agent-layer-architecture-report.md` analysis. The raw external report remains untracked; only the approved decisions are captured here.

---

## 1. Goal and non-goals

### 1.1 Goal

Add a workflow-scoped, read-only agent that helps users understand validation errors and schema constraints for a selected workflow. It runs behind the existing FastAPI backend and is consumed by the existing React frontend.

### 1.2 Non-goals

The agent layer will **not** implement any of the following in this phase:

| Out of scope | Why |
|---|---|
| Execution or run-aware tools (`run_workflow`, `submit_job`, `kill_job`, `get_run_status`, `fetch_logs`, `fetch_artifacts`) | No run lifecycle exists. |
| Execution-adjacent adapter methods (`preview_dag`, `plan_workspace`, `build_command`) | These imply a run lifecycle that does not exist. |
| Mutation tools (`apply_config_edit`, `update_samples`, `save_workflow`, `delete_workflow`) | Agent is read-only in this phase. |
| Human-in-the-loop or approval gates (`HumanApprovalGate`, pending approvals, apply confirmations) | No HITL infrastructure exists. |
| Cross-session persistence (database `AgentSession` table, long-term memory) | In-memory conversation context is sufficient for the MVP. |
| Filesystem tools (`read_file`, `list_dir`, `write_file`, `glob`, `exec`) | Agent must access platform only through services. |
| Snakemake, shell, or CLI command generation | Agent cannot recommend execution steps. |
| Provenance claims | Validation is static checking, not execution success. |
| RAG over assay policy documents | Defer until the core loop is proven. |
| Multi-workflow comparison | Defer until a general agent endpoint is added. |
| LangGraph, LangChain, and provider-specific agent SDKs (e.g., Mastra, Vercel AI SDK) backend usage | Defer until durable state, HITL, execution-aware workflows, or cross-session memory become real requirements. |

---

## 2. Architecture

The agent layer is a horizontal, read-only interpretation layer over the existing service layer. It does not introduce a new execution path.

### 2.1 Dependency direction

```
Frontend (AgentSidebar, IssuePanel)
    │
    ▼
POST /api/v1/workflows/{workflow_id}/agent/chat
    │
    ▼
AgentRouter ──► AgentService ──► LLMClient (mock → real)
    │                  │
    │                  ▼
    │         ReadOnlyToolRegistry
    │                  │
    │      ┌──────────┴──────────┐
    │      ▼                     ▼
    │ WorkflowInfoService    ValidationService
    │      │                     │
    │      ▼                     ▼
    └► WorkflowRegistry ──► WorkflowAdapter
                              ^
                              │
                    Agent code NEVER touches this directly
```

### 2.2 Layer rules

1. `AgentService` may hold only `WorkflowInfoService`, `ValidationService`, and `LLMClient`. It must never hold a `WorkflowAdapter` reference.
2. `ReadOnlyToolRegistry` rejects any tool registered with `read_only=False`.
3. All agent tools return either a platform `Result[T]` or plain text. They never return command specs, workspace plans, DAG previews, or file handles.
4. `AgentRouter` depends only on `AgentService` via `get_agent_service()`.
5. `LLMClient` is an abstraction; `AgentService` does not import provider SDKs.

---

## 3. API endpoint

### 3.1 First endpoint (workflow-scoped)

```http
POST /api/v1/workflows/{workflow_id}/agent/chat
```

The first agent is bound to the selected workflow. Route-level workflow validation resolves the workflow before the agent runs, matching the existing `/workflows/{id}/schema` and `/workflows/{id}/validate` shape.

### 3.2 Reserved for later (general/cross-workflow)

```http
POST /api/v1/agent/chat
```

A centralized endpoint is reserved for a future general assistant that can answer cross-workflow questions ("Which assay should I pick?"). It will be added only after the workflow-scoped validation assistant is proven.

### 3.3 Request/response envelope

Request body:

```json
{
  "session_id": null,
  "message": "Why am I getting FRiP errors?",
  "context": {
    "current_issues": [{"code": "FRIP_LOW", ...}],
    "current_config": {"assay": "chip_tf", ...},
    "current_schema": {...}
  }
}
```

Response envelope:

```json
{
  "ok": true,
  "session_id": null,
  "message": "FRiP (Fraction of Reads in Peaks) measures...",
  "suggestions": [],
  "tool_calls": [
    {
      "tool_name": "explain_issues",
      "input_summary": {...},
      "output_summary": "...",
      "read_only": true
    }
  ],
  "issues": []
}
```

- `ok`, `issues`, and the envelope style match the existing `ValidationResponse` contract.
- `issues` uses the same `IssueResponse` shape as other endpoints.
- `read_only` is always `true` in this phase.

### 3.4 Streaming extension (future)

After the request/response path is proven, the same endpoint may support SSE:

```http
POST /api/v1/workflows/{workflow_id}/agent/chat?stream=true
```

Event types: `message`, `tool_call`, `suggestion`, `done`, `error`.

---

## 4. Service and tool boundaries

### 4.1 WorkflowInfoService

A thin read-only service is introduced so agent tools never call `WorkflowRegistry.get(...).schema()` or touch `WorkflowAdapter` directly.

```python
class WorkflowInfoService:
    def list_workflows(self) -> list[WorkflowMetadata]: ...
    def get_schema(self, workflow_id: str) -> Result[WorkflowSchema]: ...
    def get_capabilities(self, workflow_id: str) -> Result[WorkflowCapabilities]: ...
```

Both the existing `/workflows/{id}/schema` route and the agent `get_schema` tool should consume this service.

### 4.2 First skeleton tools

The first `ReadOnlyToolRegistry` registers exactly four tools:

| Tool | Calls | Output |
|---|---|---|
| `list_workflows` | `WorkflowInfoService.list_workflows()` | `Result[list[WorkflowMetadata]]` |
| `get_workflow_schema` | `WorkflowInfoService.get_schema()` | `Result[WorkflowSchema]` |
| `validate_workflow_input` | `ValidationService.validate()` | `Result[object]` with `Issue[]` |
| `explain_issues` | LLM-only reasoning wrapper | Natural language text |

`explain_issues` is the only tool that uses the `LLMClient` boundary. In PR93 it is backed by `MockLLMClient`, so no real LLM provider is required. The other three tools are data tools that prove the read-only service boundary.

### 4.3 Deferred tools

- `summarize_validation_state` — pure LLM wrapper; add after `explain_issues` is proven.
- `suggest_edits` — generates proposed config changes; requires output filtering, schema-aware validation, disclaimers, and "validated vs unverified" labeling.
- Assay-policy RAG tools.
- Multi-workflow comparison tools.

---

## 5. LLM client boundary

### 5.1 Abstraction first

The `LLMClient` interface is defined in PR92 and used from PR93 onward. A `MockLLMClient` provides deterministic responses so the endpoint and frontend can be tested without a real provider.

```python
class LLMClient:
    async def complete(self, messages: list[Message], tools: list[Tool]) -> LLMResponse: ...
```

### 5.2 Real provider integration

Real LLM provider integration lands in PR97, gated by environment variables and not enabled by default. Tests must pass without API keys. Provider-specific prompts and API keys are kept out of earlier PRs.

---

## 6. Frontend integration

### 6.1 AgentSidebar

- Render a persistent "Validation Assistant — Read Only" label.
- Show message history with Markdown rendering.
- Show tool-call transparency ("Agent looked up schema…").
- Render suggestion cards with amber warning borders, a disclaimer, and only **Copy** and **View Diff** actions.
- Never render **Apply**, **Save**, **Submit**, or **Run** buttons.
- Never auto-inject agent text into the `ValidationWorkspace` inputs.

### 6.2 IssuePanel

Add an **Ask Agent** button to each issue. Clicking prepopulates the chat input with the issue summary and appends the full structured issue to `AgentContext.current_issues`.

### 6.3 Client boundary

- Add `AgentApiClient` in `frontend/src/api/agentClient.ts`.
- Add matching TypeScript types in `frontend/src/api/types.ts`.
- The client calls the workflow-scoped endpoint.

---

## 7. Security and safety

### 7.1 Read-only enforcement

- `ReadOnlyToolRegistry.register()` rejects `read_only=False`.
- Default composition registers only `read_only=True` tools.
- Tool-name deny-list blocks `run_*`, `submit_*`, `execute_*`, `delete_*`, `update_*`, `apply_*`, `write_*`, `modify_*`.

### 7.2 No execution or command generation

- No execution tools are registered.
- System prompt forbids execution language.
- Output filter rejects or replaces command-like strings and shell/Snakemake patterns before prompts, logs, and responses where needed.

### 7.3 Redaction

- `Issue.context` absolute paths and sensitive config keys (`api_key`, `password`, `token`) are redacted to `<REDACTED>` / `<FILE_PATH>` before LLM prompts and logs.

### 7.4 Suggestion safety

- Every `AgentSuggestion` includes a non-empty disclaimer: "This is a model-generated suggestion. Please verify against your experimental design and ENCODE guidelines."
- Suggestions may optionally be re-validated through `ValidationService` and marked "Validated by platform" or "Unverified suggestion."

### 7.5 Audit logging

Log every request, tool call, and LLM call with request ID, workflow ID, tool names, and latency. Logs must not contain raw config, samples, or file paths.

---

## 8. PR sequence

| PR | Branch | Scope |
|---|---|---|
| PR91 | `pr/agent-boundary-spec` | Docs-only: architecture, tool contract, API contract, frontend UX rules, non-goals, risks. |
| PR92 | `pr/agent-models` | Pydantic `AgentRequest`, `AgentResponse`, `AgentContext`, `AgentSuggestion`, `AgentToolCall`; matching TypeScript types; `LLMClient` protocol/interface. |
| PR93 | `pr/agent-service-skeleton` | `WorkflowInfoService`, `ReadOnlyToolRegistry`, four read-only tools, `AgentService` wired to `MockLLMClient`. |
| PR94 | `pr/agent-router` | FastAPI `POST /api/v1/workflows/{workflow_id}/agent/chat`, dependency injection, app-state wiring. |
| PR95 | `pr/frontend-agent-sidebar` | `AgentApiClient`, chat UI, issue-to-agent references, suggestion cards, read-only guardrails. |
| PR96 | `pr/agent-safety-hardening` | Redaction of paths/secrets, output filtering, audit logging, tool-call telemetry for the agent path. |
| PR97 | `pr/llm-provider-integration` | Optional real LLM provider integration behind feature flag and env gating. |

---

## 9. Test/verification matrix

| PR | Tests |
|---|---|
| PR92 | Pydantic model round-trips; TypeScript types match backend contract; `AgentSuggestion.disclaimer` non-empty default; `AgentToolCall.read_only` defaults to `true`. |
| PR93 | Registry rejects `read_only=False`; all registered tools are read-only; agent modules do not import `platform.adapters`; `validate_workflow_input` routes through `ValidationService`; `get_schema` routes through `WorkflowInfoService`; mock LLM client returns deterministic responses. |
| PR94 | Endpoint returns envelope shape; invalid `workflow_id` returns `WORKFLOW_NOT_FOUND` in `issues`; endpoint works without real LLM provider. |
| PR95 | `AgentSidebar` renders read-only label; suggestion cards have no Apply/Save/Submit/Run buttons; `AgentApiClient` uses workflow-scoped endpoint. |
| PR96 | Redaction of absolute paths and sensitive config keys; output filter rejects or replaces command-like strings and shell/Snakemake patterns; audit log entries omit raw config, samples, and file paths; regression tests run without API keys. |
| PR97 | `LLMClient` interface implemented by mock and real client; real provider gated by env var; tests run without API keys; output filtering and redaction active. |
| Cross-cutting | Import-boundary scan fails if agent layer imports `encode_pipeline.adapters`, `subprocess`, `snakemake`, or references `preview_dag`, `plan_workspace`, `build_command`. |

---

## 10. Risks and mitigations

| Risk | Mitigation |
|---|---|
| Hallucinated fixes | Inject schema into prompts; re-validate suggestions; label validated vs unverified. |
| Prompt injection coerces tool misuse | Allow-list tool names; reject non-read-only tools; log unknown tool requests. |
| File path / secret leakage | Redact `Issue.context` and sensitive config keys before prompts and logs. |
| User over-trusts suggestions | Amber warning cards; disclaimer on every suggestion; no Apply button. |
| LLM suggests execution commands | System prompt forbids; output filter rejects or replaces command-like strings and shell/Snakemake patterns; absolute paths are replaced before prompts/logs/responses. |
| Agent bypasses `ValidationService` | Route all data access through `WorkflowInfoService`/`ValidationService`; import-boundary tests. |
| Framework/provider lock-in | `LLMClient` abstraction; `AgentService` and tool registry remain framework-agnostic. |
| UI implies mutation | Read-only label, restricted button vocabulary, no auto-inject. |

---

## 11. External report handling

`agent-layer-architecture-report.md` remains untracked. Only the approved decisions captured in this repo-owned document are implemented. If the team later wants to preserve additional reasoning from the report, it should be distilled into concise repo-owned docs rather than committing the raw report.

---

## 12. Decision log

| Decision | Rationale |
|---|---|
| Workflow-scoped endpoint for MVP | The first agent is a validation assistant bound to a selected workflow; route-level validation keeps the MVP safer and aligns with existing `/workflows/{id}` routes. |
| Centralized endpoint reserved | A general/cross-workflow assistant may use `POST /api/v1/agent/chat` later. |
| `WorkflowInfoService` introduced | Prevents agent tools from touching `WorkflowAdapter` directly; enables redaction and testing. |
| Four tools in first skeleton | Minimum set to prove read-only boundary and LLM path; defer summarization and edit suggestions. |
| `LLMClient` abstraction early | Every code path exercises the real boundary; real-provider swap is configuration, not refactor. |
| Real LLM split into PR97 | Keeps API keys and provider-specific prompts out of the first provider PR. |
| Frontend types in PR92 | Keeps the API contract synchronized across the stack. |
| IssuePanel "Ask Agent" button in PR95 | Primary forcing function for `explain_issues`; makes the frontend PR valuable. |
| `agent-layer-architecture-report.md` untracked | External input, not repo-owned; distill only approved decisions. |
