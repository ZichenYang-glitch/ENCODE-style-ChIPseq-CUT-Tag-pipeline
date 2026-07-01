export type Severity = 'error' | 'warning' | 'info';

export interface Issue {
  code: string;
  message: string;
  severity: Severity;
  path: string | null;
  source: string;
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
  config: Record<string, unknown>;
  samples: string | unknown[] | null;
  options?: Record<string, unknown>;
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
