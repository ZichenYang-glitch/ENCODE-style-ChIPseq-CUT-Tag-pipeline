import type { WorkflowSchemaResponse } from '../../api/generated/models';
import {
  createAuthoringSchemaFixture,
} from '../input-workbench/test-fixtures';
import {
  readWorkbenchSchema,
  type WorkbenchSchema,
} from '../input-workbench/schemaContract';
import type {
  AuthoringCapabilityView,
  ConfigurationDraftSummary,
  ProposalView,
} from './model';

export function createDraftEditorContract(): WorkflowSchemaResponse {
  const contract = createAuthoringSchemaFixture();
  return {
    ...contract,
    config_schema: {
      ...contract.config_schema,
      properties: {
        limit: { type: 'integer', title: 'Local limit', default: 1 },
      },
      additionalProperties: false,
    },
    option_schema: {
      ...contract.option_schema,
      properties: {
        note: { type: 'string', title: 'Local note' },
      },
    },
  };
}

export function createDraftEditorSchema(): WorkbenchSchema {
  const parsed = readWorkbenchSchema(createDraftEditorContract());
  if (!parsed.ok) {
    throw new Error(parsed.message);
  }
  return parsed.value;
}

export function createCapabilityFixture(): AuthoringCapabilityView {
  return {
    available: true,
    unavailableReason: null,
    intents: ['suggest', 'explain_issues'],
    fields: [
      {
        token: 'resource.limit',
        label: 'Example limit',
        type: 'integer',
        minimum: 0,
        maximum: 8,
      },
      {
        token: 'resource.enabled',
        label: 'Example switch',
        type: 'boolean',
      },
    ],
    issues: [
      {
        code: 'EXAMPLE_ISSUE',
        label: 'Example needs a sample selection.',
      },
    ],
  };
}

export function createDraftSummary(
  sid: string,
  projectName: string,
  state: ConfigurationDraftSummary['state'] = 'active',
): ConfigurationDraftSummary {
  return {
    sid,
    projectName,
    workflowName: 'Example workflow',
    state,
    updatedLabel: '10 Nov 2026',
    expiresLabel: '17 Nov 2026',
  };
}

export function createPendingProposal(
  overrides: Partial<ProposalView> = {},
): ProposalView {
  return {
    state: 'pending',
    digestAvailable: true,
    adoptable: true,
    changes: [
      {
        label: 'Example limit',
        before: { presence: 'missing' },
        after: { presence: 'present', value: 0 },
        recognized: true,
      },
    ],
    ...overrides,
  };
}
