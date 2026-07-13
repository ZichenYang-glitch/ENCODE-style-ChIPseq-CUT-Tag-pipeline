import {
  getWorkflowSchema,
  listWorkflows,
  validateWorkflow,
} from './generated/workflows/workflows';
import {
  cancelRun,
  createRun,
  getRun,
  listRunEvents,
  listRunLogs,
  startRun,
} from './generated/runs/runs';
import { triggerPreflight } from './generated/preflight/preflight';
import { chatWithWorkflowAgent } from './generated/agent/agent';
import type {
  AgentResponse as GeneratedAgentResponse,
  AgentSuggestion as GeneratedAgentSuggestion,
  AgentToolCall as GeneratedAgentToolCall,
  IssueResponse as GeneratedIssue,
  RunEventResponse as GeneratedRunEvent,
  RunLogChunkResponse as GeneratedRunLogChunk,
  RunRecordResponse as GeneratedRunRecord,
  ValidationRequest,
} from './generated/models';
import { ApiError } from './fetcher';
import type { AgentApiClient } from './agentClient';
import type { WorkflowApiClient } from './client';
import type { RunApiClient } from './runClient';
import type {
  RunEventResponse,
  RunLogChunkResponse,
  RunRecordResponse,
  RunResponse,
} from './runTypes';
import type {
  AgentResponse,
  AgentSuggestion,
  AgentToolCall,
  Issue,
  Severity,
  WorkflowSchema,
} from './types';

function asRecord(value: unknown): Record<string, unknown> {
  return value !== null && typeof value === 'object' && !Array.isArray(value)
    ? (value as Record<string, unknown>)
    : {};
}

function asStringRecord(value: unknown): Record<string, string> {
  const record = asRecord(value);
  return Object.fromEntries(
    Object.entries(record).filter((entry): entry is [string, string] => {
      return typeof entry[1] === 'string';
    }),
  );
}

function toSeverity(value: string | undefined): Severity {
  return value === 'warning' || value === 'info' ? value : 'error';
}

function toIssue(issue: GeneratedIssue): Issue {
  return {
    code: issue.code,
    message: issue.message,
    severity: toSeverity(issue.severity),
    path: issue.path ?? null,
    source: issue.source ?? null,
    technical_message: issue.technical_message ?? null,
    hint: issue.hint ?? null,
    context: asRecord(issue.context),
  };
}

function issuesFromError(error: ApiError): Issue[] {
  if (error.issues.length > 0) {
    return error.issues.map((issue) => ({
      code: issue.code,
      message: issue.message,
      severity: toSeverity(issue.severity),
      path: issue.path ?? null,
      source: issue.source ?? null,
      technical_message: null,
      hint: issue.hint ?? null,
      context: {},
    }));
  }
  return [
    {
      code: error.code,
      message: error.message,
      severity: 'error',
      path: null,
      source: 'api',
      technical_message: null,
      hint: null,
      context: {},
    },
  ];
}

function mapIssues(issues: GeneratedIssue[] | undefined): Issue[] {
  return (issues ?? []).map(toIssue);
}

function toWorkflowSchema(value: unknown): WorkflowSchema {
  const hints = asRecord(value);
  return {
    config_schema: asRecord(hints.config_schema),
    sample_schema: asRecord(hints.sample_schema),
    option_schema: asRecord(hints.option_schema),
  };
}

function toRunRecord(run: GeneratedRunRecord): RunRecordResponse {
  return {
    run_id: run.run_id,
    workflow_id: run.workflow_id,
    inputs: asRecord(run.inputs),
    status: run.status,
    created_at: run.created_at,
    updated_at: run.updated_at,
    started_at: run.started_at,
    ended_at: run.ended_at,
    current_stage: run.current_stage,
    cancellation_reason: run.cancellation_reason,
    error: run.error ? toIssue(run.error) : null,
    tags: asStringRecord(run.tags),
  };
}

function toRunEvent(event: GeneratedRunEvent): RunEventResponse {
  return {
    event_id: event.event_id,
    run_id: event.run_id,
    sequence: event.sequence,
    event_type: event.event_type,
    timestamp: event.timestamp,
    status: event.status ?? null,
    stage: event.stage ?? null,
    message: event.message,
    context: asRecord(event.context),
    issue: event.issue ? toIssue(event.issue) : null,
  };
}

function toRunLogChunk(chunk: GeneratedRunLogChunk): RunLogChunkResponse {
  return {
    chunk_id: chunk.chunk_id,
    run_id: chunk.run_id,
    stream_name: chunk.stream_name,
    sequence: chunk.sequence,
    timestamp: chunk.timestamp,
    lines: chunk.lines,
  };
}

function toRunResponse(
  response: Awaited<ReturnType<typeof getRun>>,
): RunResponse {
  return {
    ok: response.ok,
    run: response.run ? toRunRecord(response.run) : null,
    issues: mapIssues(response.issues),
  };
}

function toAgentToolCall(toolCall: GeneratedAgentToolCall): AgentToolCall {
  return {
    tool_name: toolCall.tool_name,
    input_summary: asRecord(toolCall.input_summary),
    output_summary: toolCall.output_summary ?? '',
    read_only: true,
  };
}

function toAgentSuggestion(
  suggestion: GeneratedAgentSuggestion,
): AgentSuggestion {
  return {
    type: suggestion.type ?? 'general',
    description: suggestion.description,
    target_path: suggestion.target_path ?? null,
    current_value: suggestion.current_value ?? null,
    proposed_value: suggestion.proposed_value ?? null,
    rationale: suggestion.rationale ?? null,
    disclaimer: suggestion.disclaimer ?? 'Verify this suggestion before use.',
  };
}

function toAgentResponse(response: GeneratedAgentResponse): AgentResponse {
  return {
    ok: response.ok,
    session_id: response.session_id ?? null,
    message: response.message,
    suggestions: (response.suggestions ?? []).map(toAgentSuggestion),
    tool_calls: (response.tool_calls ?? []).map(toAgentToolCall),
    issues: mapIssues(response.issues),
  };
}

export function createGeneratedWorkflowClient(): WorkflowApiClient {
  return {
    async listWorkflows() {
      try {
        const response = await listWorkflows();
        return {
          ok: response.ok,
          workflows: response.workflows.map((workflow) => ({
            metadata: {
              workflow_id: workflow.metadata.workflow_id,
              name: workflow.metadata.name,
              version: workflow.metadata.version,
              description: workflow.metadata.description ?? '',
              engines: workflow.metadata.engines ?? [],
              tags: workflow.metadata.tags ?? [],
            },
            capabilities: {
              supports: workflow.capabilities.supports ?? [],
            },
          })),
          issues: mapIssues(response.issues),
        };
      } catch (error) {
        if (!(error instanceof ApiError)) throw error;
        return { ok: false, workflows: [], issues: issuesFromError(error) };
      }
    },

    async getWorkflowSchema(workflowId) {
      try {
        const response = await getWorkflowSchema(workflowId);
        return {
          ok: response.ok,
          workflow_id: response.workflow_id,
          schema_hints: response.ok
            ? toWorkflowSchema(response.schema)
            : null,
          issues: mapIssues(response.issues),
        };
      } catch (error) {
        if (!(error instanceof ApiError)) throw error;
        return {
          ok: false,
          workflow_id: workflowId,
          schema_hints: null,
          issues: issuesFromError(error),
        };
      }
    },

    async validateWorkflow(workflowId, inputs) {
      try {
        const response = await validateWorkflow(workflowId, {
          config: inputs.config,
          samples: inputs.samples as ValidationRequest['samples'],
          options: inputs.options ?? {},
        });
        return {
          ok: response.ok,
          workflow_id: response.workflow_id ?? workflowId,
          value: response.value ?? null,
          snapshot: response.snapshot ?? null,
          issues: mapIssues(response.issues),
        };
      } catch (error) {
        if (!(error instanceof ApiError)) throw error;
        return {
          ok: false,
          workflow_id: workflowId,
          value: null,
          snapshot: null,
          issues: issuesFromError(error),
        };
      }
    },
  };
}

export function createGeneratedRunClient(): RunApiClient {
  return {
    async createRun(workflowId, request) {
      try {
        return toRunResponse(await createRun(workflowId, request));
      } catch (error) {
        if (!(error instanceof ApiError)) throw error;
        return { ok: false, run: null, issues: issuesFromError(error) };
      }
    },

    async preflightRun(runId) {
      try {
        return toRunResponse(await triggerPreflight(runId));
      } catch (error) {
        if (!(error instanceof ApiError)) throw error;
        return { ok: false, run: null, issues: issuesFromError(error) };
      }
    },

    async getRun(runId) {
      try {
        return toRunResponse(await getRun(runId));
      } catch (error) {
        if (!(error instanceof ApiError)) throw error;
        return { ok: false, run: null, issues: issuesFromError(error) };
      }
    },

    async startRun(runId) {
      try {
        return toRunResponse(await startRun(runId));
      } catch (error) {
        if (!(error instanceof ApiError)) throw error;
        return { ok: false, run: null, issues: issuesFromError(error) };
      }
    },

    async listRunEvents(runId, options = {}) {
      try {
        const response = await listRunEvents(runId, {
          after: options.after,
          limit: options.limit,
        });
        return {
          ok: response.ok,
          run_id: response.run_id,
          events: (response.events ?? []).map(toRunEvent),
          next_cursor: response.next_cursor ?? null,
          issues: mapIssues(response.issues),
        };
      } catch (error) {
        if (!(error instanceof ApiError)) throw error;
        return {
          ok: false,
          run_id: runId,
          events: [],
          next_cursor: null,
          issues: issuesFromError(error),
        };
      }
    },

    async listRunLogs(runId, options = {}) {
      const streamName = options.streamName ?? 'stdout';
      try {
        const response = await listRunLogs(runId, {
          stream_name: streamName,
          after: options.after,
          limit: options.limit,
        });
        return {
          ok: response.ok,
          run_id: response.run_id,
          stream_name: response.stream_name,
          chunks: (response.chunks ?? []).map(toRunLogChunk),
          next_cursor: response.next_cursor ?? null,
          issues: mapIssues(response.issues),
        };
      } catch (error) {
        if (!(error instanceof ApiError)) throw error;
        return {
          ok: false,
          run_id: runId,
          stream_name: streamName,
          chunks: [],
          next_cursor: null,
          issues: issuesFromError(error),
        };
      }
    },

    async cancelRun(runId) {
      try {
        return toRunResponse(await cancelRun(runId));
      } catch (error) {
        if (!(error instanceof ApiError)) throw error;
        return { ok: false, run: null, issues: issuesFromError(error) };
      }
    },
  };
}

export function createGeneratedAgentClient(): AgentApiClient {
  return {
    async chat(workflowId, request) {
      return toAgentResponse(await chatWithWorkflowAgent(workflowId, request));
    },
  };
}
