import type {
  GetWorkflowSchemaResponse,
  ListWorkflowsResponse,
  ValidateWorkflowResponse,
  WorkflowInputs,
} from './types';
import {
  stubValidationResponses,
  stubWorkflows,
  stubWorkflowSchemas,
} from '../data/stubWorkflows';

export interface WorkflowApiClient {
  listWorkflows(): Promise<ListWorkflowsResponse>;
  getWorkflowSchema(workflowId: string): Promise<GetWorkflowSchemaResponse>;
  validateWorkflow(
    workflowId: string,
    inputs: WorkflowInputs,
  ): Promise<ValidateWorkflowResponse>;
}

function workflowNotFoundIssue(workflowId: string) {
  return {
    code: 'WORKFLOW_NOT_FOUND',
    message: 'Workflow was not found.',
    severity: 'error' as const,
    path: 'workflow_id',
    source: 'registry',
    technical_message: null,
    hint: null,
    context: { workflow_id: workflowId },
  };
}

export function createStubWorkflowClient(): WorkflowApiClient {
  return {
    async listWorkflows() {
      return {
        ok: true,
        workflows: stubWorkflows,
        issues: [],
      };
    },

    async getWorkflowSchema(workflowId) {
      const schemaHints = stubWorkflowSchemas[workflowId];
      if (!schemaHints) {
        return {
          ok: false,
          workflow_id: workflowId,
          schema_hints: null,
          issues: [workflowNotFoundIssue(workflowId)],
        };
      }
      return {
        ok: true,
        workflow_id: workflowId,
        schema_hints: schemaHints,
        issues: [],
      };
    },

    async validateWorkflow(workflowId, inputs) {
      const responder = stubValidationResponses[workflowId];
      if (!responder) {
        return {
          ok: false,
          workflow_id: workflowId,
          value: null,
          snapshot: null,
          issues: [workflowNotFoundIssue(workflowId)],
        };
      }
      return responder(inputs);
    },
  };
}
