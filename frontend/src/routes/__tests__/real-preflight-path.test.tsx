import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { createMemoryRouter, RouterProvider } from 'react-router-dom';
import { beforeEach, describe, expect, it, vi } from 'vitest';
import { ClientProvider } from '../../api/client-context';
import { appRoutes } from '../../app/router';
import {
  getWorkflowSchema,
  listWorkflows,
  validateWorkflow,
} from '../../api/generated/workflows/workflows';
import {
  createRun,
  getRun,
  listRunEvents,
  listRunLogs,
} from '../../api/generated/runs/runs';
import { triggerPreflight } from '../../api/generated/preflight/preflight';

vi.mock('../../api/generated/workflows/workflows', () => ({
  getWorkflowSchema: vi.fn(),
  listWorkflows: vi.fn(),
  validateWorkflow: vi.fn(),
}));
vi.mock('../../api/generated/runs/runs', () => ({
  cancelRun: vi.fn(),
  createRun: vi.fn(),
  getRun: vi.fn(),
  listRunEvents: vi.fn(),
  listRunLogs: vi.fn(),
  startRun: vi.fn(),
}));
vi.mock('../../api/generated/preflight/preflight', () => ({
  triggerPreflight: vi.fn(),
}));
vi.mock('../../api/generated/agent/agent', () => ({
  chatWithWorkflowAgent: vi.fn(),
}));

const WORKFLOW_ID = 'encode-style-chipseq-cuttag-atac-mnase';

function runRecord(status: string) {
  return {
    run_id: 'run-real-1',
    workflow_id: WORKFLOW_ID,
    inputs: {
      config: { genome: 'hg38' },
      samples: 'samples.tsv',
      options: {},
    },
    status,
    created_at: '2026-07-11T00:00:00.000Z',
    updated_at: '2026-07-11T00:01:00.000Z',
    started_at: null,
    ended_at: null,
    current_stage: 'preflight',
    cancellation_reason: null,
    error: null,
    tags: {},
  };
}

describe('real preflight product path', () => {
  beforeEach(() => {
    vi.clearAllMocks();
    vi.mocked(listWorkflows).mockResolvedValue({
      ok: true,
      workflows: [
        {
          metadata: {
            workflow_id: WORKFLOW_ID,
            name: 'ENCODE workflow',
            version: '0.3.0',
            description: 'Epigenomics workflow.',
            engines: ['snakemake'],
            tags: ['chip-seq'],
          },
          capabilities: { supports: ['validation', 'workspace_plan'] },
        },
      ],
      issues: [],
    });
    vi.mocked(getWorkflowSchema).mockResolvedValue({
      ok: true,
      workflow_id: WORKFLOW_ID,
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
        sample_schema: { type: 'string' },
        option_schema: { type: 'object' },
      },
      issues: [],
    });
    vi.mocked(validateWorkflow).mockResolvedValue({
      ok: true,
      workflow_id: WORKFLOW_ID,
      value: null,
      snapshot: {
        snapshot_id: 'vsnap_0123456789abcdef0123456789abcdef',
        workflow_id: WORKFLOW_ID,
        schema_version: '1.0.0',
        adapter_version: '0.3.0',
        payload_digest: 'a'.repeat(64),
        validated_at: '2026-07-14T00:00:00.000Z',
        expires_at: '2026-07-14T00:30:00.000Z',
      },
      issues: [],
    });
    vi.mocked(createRun).mockResolvedValue({
      ok: true,
      run: runRecord('created'),
      issues: [],
    });
    vi.mocked(triggerPreflight).mockResolvedValue({
      ok: true,
      run: runRecord('validating'),
      issues: [],
    });
    vi.mocked(getRun).mockResolvedValue({
      ok: true,
      run: runRecord('planned'),
      issues: [],
    });
    vi.mocked(listRunEvents).mockResolvedValue({
      ok: true,
      run_id: 'run-real-1',
      events: [
        {
          event_id: 'evt-1',
          run_id: 'run-real-1',
          sequence: 1,
          event_type: 'preflight_completed',
          timestamp: '2026-07-11T00:01:00.000Z',
          status: 'planned',
          stage: 'preflight',
          message: 'Local preflight completed.',
          context: {},
          issue: null,
        },
      ],
      next_cursor: null,
      issues: [],
    });
    vi.mocked(listRunLogs).mockImplementation(async (_, params = {}) => ({
      ok: true,
      run_id: 'run-real-1',
      stream_name: params.stream_name ?? 'stdout',
      chunks:
        params.stream_name === 'stderr'
          ? []
          : [
              {
                chunk_id: 'log-1',
                run_id: 'run-real-1',
                stream_name: 'stdout',
                sequence: 1,
                timestamp: '2026-07-11T00:01:00.000Z',
                lines: ['Dry-run completed successfully.'],
              },
            ],
      next_cursor: null,
      issues: [],
    }));
  });

  it('uses generated operations for validate, create, preflight, and run detail', async () => {
    let resolvePreflight!: (value: {
      ok: boolean;
      run: ReturnType<typeof runRecord>;
      issues: never[];
    }) => void;
    let preflightFinished = false;
    vi.mocked(triggerPreflight).mockImplementation(
      () =>
        new Promise((resolve) => {
          resolvePreflight = resolve;
        }),
    );
    vi.mocked(getRun).mockImplementation(async () => ({
      ok: true,
      run: runRecord(preflightFinished ? 'planned' : 'created'),
      issues: [],
    }));
    const user = userEvent.setup();
    const queryClient = new QueryClient({
      defaultOptions: { queries: { retry: false } },
    });
    const router = createMemoryRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}`],
    });
    render(
      <QueryClientProvider client={queryClient}>
        <ClientProvider>
          <RouterProvider router={router} />
        </ClientProvider>
      </QueryClientProvider>,
    );

    const samplesInput = await screen.findByLabelText(/Samples \(path string\)/i);
    await user.type(
      samplesInput,
      'samples.tsv',
    );
    await user.click(screen.getByTestId('validate-button'));
    expect(await screen.findByTestId('create-run-button')).toBeEnabled();

    await user.click(screen.getByTestId('create-run-button'));

    await waitFor(() => {
      expect(router.state.location.pathname).toBe('/runs/run-real-1');
    });
    expect(triggerPreflight).toHaveBeenCalledWith('run-real-1');
    expect(await screen.findByTestId('run-status-badge')).toHaveTextContent(
      'created',
    );

    preflightFinished = true;
    resolvePreflight({ ok: true, run: runRecord('validating'), issues: [] });
    await waitFor(() => {
      expect(screen.getByTestId('run-status-badge')).toHaveTextContent('planned');
    });
    expect(screen.getByText('Local preflight completed.')).toBeInTheDocument();
    expect(
      screen.getByText('Dry-run completed successfully.'),
    ).toBeInTheDocument();

    expect(validateWorkflow).toHaveBeenCalledTimes(1);
    expect(createRun).toHaveBeenCalledTimes(1);
    expect(triggerPreflight).toHaveBeenCalledWith('run-real-1');
    expect(getRun).toHaveBeenCalledWith('run-real-1');
    expect(listRunEvents).toHaveBeenCalled();
    expect(vi.mocked(listRunLogs).mock.calls.length).toBeGreaterThanOrEqual(2);
  });

  it('does not restore a stale snapshot when legacy inputs change during validation', async () => {
    let resolveValidation!: (value: Awaited<ReturnType<typeof validateWorkflow>>) => void;
    vi.mocked(validateWorkflow).mockImplementation(
      () =>
        new Promise((resolve) => {
          resolveValidation = resolve;
        }),
    );
    const user = userEvent.setup();
    const queryClient = new QueryClient({
      defaultOptions: { queries: { retry: false } },
    });
    const router = createMemoryRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}`],
    });
    render(
      <QueryClientProvider client={queryClient}>
        <ClientProvider>
          <RouterProvider router={router} />
        </ClientProvider>
      </QueryClientProvider>,
    );

    const samplesInput = await screen.findByLabelText(/Samples \(path string\)/i);
    await user.type(samplesInput, 'samples.tsv');
    await user.click(screen.getByTestId('validate-button'));
    await user.type(samplesInput, '.changed');

    resolveValidation({
      ok: true,
      workflow_id: WORKFLOW_ID,
      value: null,
      snapshot: {
        snapshot_id: 'vsnap_0123456789abcdef0123456789abcdef',
        workflow_id: WORKFLOW_ID,
        schema_version: '1.0.0',
        adapter_version: '0.3.0',
        payload_digest: 'a'.repeat(64),
        validated_at: '2026-07-14T00:00:00.000Z',
        expires_at: '2026-07-14T00:30:00.000Z',
      },
      issues: [],
    });

    await waitFor(() => expect(validateWorkflow).toHaveBeenCalledTimes(1));
    expect(screen.getByTestId('create-run-button')).toBeDisabled();
  });

  it('keeps the durable run URL and offers retry when preflight fails', async () => {
    vi.mocked(triggerPreflight).mockResolvedValue({
      ok: false,
      run: null,
      issues: [
        {
          code: 'PREFLIGHT_FAILED',
          message: 'Preflight did not complete.',
          severity: 'error',
          source: 'preflight',
        },
      ],
    });
    vi.mocked(getRun).mockResolvedValue({
      ok: true,
      run: runRecord('created'),
      issues: [],
    });
    const user = userEvent.setup();
    const queryClient = new QueryClient({
      defaultOptions: { queries: { retry: false } },
    });
    const router = createMemoryRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}`],
    });
    render(
      <QueryClientProvider client={queryClient}>
        <ClientProvider>
          <RouterProvider router={router} />
        </ClientProvider>
      </QueryClientProvider>,
    );

    await user.type(await screen.findByLabelText(/Samples \(path string\)/i), 'samples.tsv');
    await user.click(screen.getByTestId('validate-button'));
    await user.click(await screen.findByTestId('create-run-button'));

    await waitFor(() => expect(router.state.location.pathname).toBe('/runs/run-real-1'));
    expect(await screen.findByRole('alert')).toHaveTextContent(
      'PREFLIGHT_FAILED: Preflight did not complete.',
    );
    expect(screen.getByTestId('run-status-badge')).toHaveTextContent('created');
    expect(screen.getByRole('button', { name: 'Run preflight' })).toBeEnabled();
  });
});
