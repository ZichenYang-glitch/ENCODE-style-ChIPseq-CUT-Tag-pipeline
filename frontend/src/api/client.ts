import type {
  GetWorkflowResponse,
  GetWorkflowSchemaResponse,
  ListWorkflowsResponse,
  ValidateWorkflowResponse,
  WorkflowInputs,
  WorkflowProjectSampleSelection,
} from './types';
import {
  stubValidationResponses,
  stubWorkflows,
  stubWorkflowSchemas,
} from '../data/stubWorkflows';

export interface WorkflowApiClient {
  listWorkflows(): Promise<ListWorkflowsResponse>;
  getWorkflow(workflowId: string): Promise<GetWorkflowResponse>;
  getWorkflowSchema(workflowId: string): Promise<GetWorkflowSchemaResponse>;
  validateWorkflow(
    workflowId: string,
    inputs: WorkflowInputs,
    projectSampleSelection?: WorkflowProjectSampleSelection,
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

    async getWorkflow(workflowId) {
      const workflow = stubWorkflows.find(
        (candidate) => candidate.metadata.workflow_id === workflowId,
      );
      if (!workflow) {
        return {
          ok: false,
          workflow_id: workflowId,
          workflow: null,
          issues: [workflowNotFoundIssue(workflowId)],
        };
      }
      return {
        ok: true,
        workflow_id: workflowId,
        workflow,
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

    async validateWorkflow(workflowId, inputs, projectSampleSelection) {
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
      if (projectSampleSelection !== undefined) {
        return {
          ok: false,
          workflow_id: workflowId,
          value: null,
          snapshot: null,
          issues: [
            {
              code: 'DATA_BINDING_SELECTION_UNAVAILABLE',
              message:
                'Project/Sample selection is unavailable in the local stub client.',
              severity: 'error',
              path: 'project_id',
              source: 'registry',
              technical_message: null,
              hint: null,
              context: {},
            },
          ],
        };
      }
      return responder(inputs);
    },
  };
}
