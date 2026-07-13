import type {
  ValidationRequestConfig,
  ValidationRequestOptions,
  ValidationRequestSamples,
} from './generated/models';

export type Severity = 'error' | 'warning' | 'info';

export interface Issue {
  code: string;
  message: string;
  severity: Severity;
  path: string | null;
  source: string | null;
  technical_message: string | null;
  hint: string | null;
  context: Record<string, unknown>;
}

export interface WorkflowMetadata {
  workflow_id: string;
  name: string;
  version: string;
  description: string;
  engines: string[];
  tags: string[];
}

export interface WorkflowCapabilities {
  supports: string[];
}

export interface WorkflowSummary {
  metadata: WorkflowMetadata;
  capabilities: WorkflowCapabilities;
}

export interface WorkflowSchema {
  config_schema: Record<string, unknown>;
  sample_schema: Record<string, unknown>;
  option_schema: Record<string, unknown>;
}

export interface WorkflowInputs {
  config: ValidationRequestConfig;
  samples: ValidationRequestSamples;
  options?: ValidationRequestOptions;
}

export interface ListWorkflowsResponse {
  ok: boolean;
  workflows: WorkflowSummary[];
  issues: Issue[];
}

export interface GetWorkflowSchemaResponse {
  ok: boolean;
  workflow_id: string;
  schema_hints: WorkflowSchema | null;
  issues: Issue[];
}

export interface ValidateWorkflowResponse {
  ok: boolean;
  workflow_id: string;
  value: unknown;
  issues: Issue[];
}

export interface AgentApiError extends Error {
  status: number;
  statusText: string;
}

export interface AgentToolCall {
  tool_name: string;
  input_summary: Record<string, unknown>;
  output_summary: string;
  read_only: true;
}

export interface AgentSuggestion {
  type: 'config_edit' | 'schema_hint' | 'assay_selection' | 'general';
  description: string;
  target_path: string | null;
  current_value: unknown | null;
  proposed_value: unknown | null;
  rationale: string | null;
  disclaimer: string;
}

export interface AgentContext {
  current_issues: Issue[];
  current_config: Record<string, unknown>;
  current_schema: Record<string, unknown>;
}

export interface AgentRequest {
  session_id: string | null;
  message: string;
  context: AgentContext | null;
}

export interface AgentResponse {
  ok: boolean;
  session_id: string | null;
  message: string;
  suggestions: AgentSuggestion[];
  tool_calls: AgentToolCall[];
  issues: Issue[];
}
