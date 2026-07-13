import { beforeEach, describe, expect, it, vi } from 'vitest';
import { ApiError } from './fetcher';
import {
  createGeneratedAgentClient,
  createGeneratedRunClient,
  createGeneratedWorkflowClient,
} from './generated-client-adapters';
import {
  getWorkflowSchema,
  listWorkflows,
  validateWorkflow,
} from './generated/workflows/workflows';
import { getRun, startRun } from './generated/runs/runs';
import { triggerPreflight } from './generated/preflight/preflight';
import { chatWithWorkflowAgent } from './generated/agent/agent';

vi.mock('./generated/workflows/workflows', () => ({
  getWorkflowSchema: vi.fn(),
  listWorkflows: vi.fn(),
  validateWorkflow: vi.fn(),
}));
vi.mock('./generated/runs/runs', () => ({
  cancelRun: vi.fn(),
  createRun: vi.fn(),
  getRun: vi.fn(),
  startRun: vi.fn(),
  listRunEvents: vi.fn(),
  listRunLogs: vi.fn(),
}));
vi.mock('./generated/preflight/preflight', () => ({
  triggerPreflight: vi.fn(),
}));
vi.mock('./generated/agent/agent', () => ({
  chatWithWorkflowAgent: vi.fn(),
}));

describe('generated client adapters', () => {
  beforeEach(() => {
    vi.clearAllMocks();
  });

  it('normalizes optional workflow transport fields', async () => {
    vi.mocked(listWorkflows).mockResolvedValue({
      ok: true,
      workflows: [
        {
          metadata: {
            workflow_id: 'encode',
            name: 'ENCODE',
            version: '1.0.0',
          },
          capabilities: {},
        },
      ],
    });

    const response = await createGeneratedWorkflowClient().listWorkflows();

    expect(response).toEqual({
      ok: true,
      workflows: [
        {
          metadata: {
            workflow_id: 'encode',
            name: 'ENCODE',
            version: '1.0.0',
            description: '',
            engines: [],
            tags: [],
          },
          capabilities: { supports: [] },
        },
      ],
      issues: [],
    });
  });

  it('uses generated validation and schema operations', async () => {
    vi.mocked(getWorkflowSchema).mockResolvedValue({
      ok: true,
      workflow_id: 'encode',
      schema: {
        schema_version: '1.0.0',
        schema_dialect: 'https://json-schema.org/draft/2020-12/schema',
        coverage: {
          config: 'partial',
          samples: 'complete',
          options: 'complete',
        },
        authoring_modes: {
          config: ['schema_form', 'yaml'],
          samples: ['tsv_upload', 'inline_table'],
          options: ['schema_form'],
        },
        input_modes: {
          config: ['object'],
          samples: ['inline_rows', 'server_path'],
          options: ['object'],
        },
        limits: {
          max_request_bytes: 2_097_152,
          max_sample_rows: 1_000,
          max_sample_columns: 64,
          max_sample_column_name_length: 128,
          max_sample_cell_length: 4_096,
        },
        config_schema: { type: 'object' },
        sample_schema: { type: 'table' },
        option_schema: {},
      },
      issues: [],
    });
    vi.mocked(validateWorkflow).mockResolvedValue({
      ok: true,
      workflow_id: 'encode',
      value: null,
      issues: [],
    });
    const client = createGeneratedWorkflowClient();

    const schema = await client.getWorkflowSchema('encode');
    const validation = await client.validateWorkflow('encode', {
      config: { genome: 'hg38' },
      samples: 'samples.tsv',
      options: {},
    });

    expect(schema.schema_hints?.config_schema).toEqual({ type: 'object' });
    expect(validation.ok).toBe(true);
    expect(validateWorkflow).toHaveBeenCalledWith('encode', {
      config: { genome: 'hg38' },
      samples: 'samples.tsv',
      options: {},
    });
  });

  it('calls the generated preflight operation', async () => {
    vi.mocked(triggerPreflight).mockResolvedValue({
      ok: true,
      run: null,
      issues: [],
    });

    const response = await createGeneratedRunClient().preflightRun('run-1');

    expect(triggerPreflight).toHaveBeenCalledWith('run-1');
    expect(response).toEqual({ ok: true, run: null, issues: [] });
  });

  it('calls the generated start operation', async () => {
    vi.mocked(startRun).mockResolvedValue({ ok: true, run: null, issues: [] });

    const response = await createGeneratedRunClient().startRun('run-1');

    expect(startRun).toHaveBeenCalledWith('run-1');
    expect(response).toEqual({ ok: true, run: null, issues: [] });
  });

  it('converts safe ApiError details into a run issue envelope', async () => {
    vi.mocked(getRun).mockRejectedValue(
      new ApiError(404, 'RUN_NOT_FOUND', 'Run was not found.', [
        {
          code: 'RUN_NOT_FOUND',
          message: 'Run was not found.',
          severity: 'error',
          path: 'run_id',
          source: 'repository',
          hint: 'Refresh the run list.',
        },
      ]),
    );

    const response = await createGeneratedRunClient().getRun('missing');

    expect(response.ok).toBe(false);
    expect(response.issues[0]).toMatchObject({
      code: 'RUN_NOT_FOUND',
      message: 'Run was not found.',
      path: 'run_id',
      source: 'repository',
      hint: 'Refresh the run list.',
      technical_message: null,
      context: {},
    });
  });

  it('normalizes optional agent response collections', async () => {
    vi.mocked(chatWithWorkflowAgent).mockResolvedValue({
      ok: true,
      message: 'Ready.',
    });

    const response = await createGeneratedAgentClient().chat('encode', {
      session_id: null,
      message: 'Explain this.',
      context: null,
    });

    expect(response).toEqual({
      ok: true,
      session_id: null,
      message: 'Ready.',
      suggestions: [],
      tool_calls: [],
      issues: [],
    });
  });
});
