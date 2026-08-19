import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { createMemoryRouter, RouterProvider } from 'react-router-dom';
import { beforeEach, describe, expect, it, vi } from 'vitest';
import { ClientProvider } from '../../api/client-context';
import { AuthProvider } from '../../app/auth';
import { appRoutes } from '../../app/router';
import { createAuthoringSchemaFixture } from '../../features/input-workbench/test-fixtures';
import {
  getWorkflow,
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

const referenceProfileMocks = vi.hoisted(() => ({
  listCompatibleReferenceProfiles: vi.fn(),
}));

vi.mock('../../api/generated/workflows/workflows', () => ({
  getWorkflow: vi.fn(),
  getWorkflowSchema: vi.fn(),
  listCompatibleReferenceProfiles:
    referenceProfileMocks.listCompatibleReferenceProfiles,
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
const REFERENCE_PROFILE = {
  profile_id: `refp_${'1'.repeat(32)}`,
  revision_id: `refpr_${'1'.repeat(32)}`,
  revision_number: 1,
  display_name: 'Human GRCh38',
  organism: 'Homo sapiens',
  assembly: 'GRCh38',
  identity_sha256: '1'.repeat(64),
};

function runRecord(status: string) {
  return {
    run_id: 'run-real-1',
    workflow_id: WORKFLOW_ID,
    status,
    created_at: '2026-07-11T00:00:00.000Z',
    updated_at: '2026-07-11T00:01:00.000Z',
    started_at: null,
    ended_at: null,
    current_stage: 'preflight',
    cancellation_reason: null,
    error: null,
    tags: {},
    reference_profile: REFERENCE_PROFILE,
  };
}

function validationResponse() {
  return {
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
      project_id: 'prj_00000000000000000000000000000000',
      binding_mode: 'legacy_v1' as const,
      provenance: 'unresolved' as const,
      sample_revision_ids: [],
      binding_digest: 'b'.repeat(64),
      input_binding: {
        mode: 'compatibility_unresolved_v1' as const,
        adapter_contract_version: null,
        digest: 'c'.repeat(64),
        fully_managed: false,
        input_uses: [],
      },
      reference_profile: REFERENCE_PROFILE,
    },
    issues: [],
  };
}

function renderProductPath() {
  const queryClient = new QueryClient({
    defaultOptions: { queries: { retry: false } },
  });
  const router = createMemoryRouter(appRoutes, {
    initialEntries: [`/workflows/${WORKFLOW_ID}/new-run`],
  });
  render(
    <QueryClientProvider client={queryClient}>
      <ClientProvider>
        <AuthProvider>
          <RouterProvider router={router} />
        </AuthProvider>
      </ClientProvider>
    </QueryClientProvider>,
  );
  return router;
}

async function authorValidInputs(user: ReturnType<typeof userEvent.setup>) {
  await user.click(await screen.findByRole('tab', { name: 'Samples' }));
  const contents =
    'sample\tfastq_1\tlayout\nS1\t/data/S1.fastq.gz\tSE\n';
  const file = new File([contents], 'samples.tsv', {
    type: 'text/tab-separated-values',
  });
  Object.defineProperty(file, 'text', {
    value: vi.fn().mockResolvedValue(contents),
  });
  await user.upload(
    screen.getByLabelText('Import samples TSV'),
    file,
  );
  await screen.findAllByDisplayValue('S1');
  await user.click(screen.getByRole('tab', { name: 'Review' }));
}

describe('real preflight product path', () => {
  beforeEach(() => {
    vi.clearAllMocks();
    const workflow = {
      metadata: {
        workflow_id: WORKFLOW_ID,
        name: 'ENCODE workflow',
        version: '0.3.0',
        description: 'Epigenomics workflow.',
        engines: ['snakemake'],
        tags: ['chip-seq'],
      },
      schema_version: '1.0.0',
      capabilities: { supports: ['validation', 'workspace_plan'] },
      upstream_identity: null,
      availability: {
        authoring: 'available' as const,
        execution: 'available' as const,
        reason_code: 'WORKFLOW_EXECUTION_READY' as const,
      },
    };
    vi.mocked(listWorkflows).mockResolvedValue({
      ok: true,
      workflows: [workflow],
      issues: [],
    });
    vi.mocked(getWorkflow).mockResolvedValue({
      ok: true,
      workflow_id: WORKFLOW_ID,
      workflow,
      issues: [],
    });
    vi.mocked(getWorkflowSchema).mockResolvedValue({
      ok: true,
      workflow_id: WORKFLOW_ID,
      schema: createAuthoringSchemaFixture(),
      issues: [],
    });
    referenceProfileMocks.listCompatibleReferenceProfiles.mockResolvedValue({
      ok: true,
      workflow_id: WORKFLOW_ID,
      profiles: [REFERENCE_PROFILE],
      issues: [],
    });
    vi.mocked(validateWorkflow).mockResolvedValue(validationResponse());
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
    const router = renderProductPath();
    await authorValidInputs(user);

    await user.click(screen.getByTestId('validate-draft-button'));
    expect(
      await screen.findByTestId('create-validated-run-button'),
    ).toBeEnabled();
    await user.click(screen.getByTestId('create-validated-run-button'));

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
    expect(screen.getByText('Dry-run completed successfully.')).toBeInTheDocument();
    expect(validateWorkflow).toHaveBeenCalledWith(
      WORKFLOW_ID,
      expect.objectContaining({
        reference_profile_revision_id: REFERENCE_PROFILE.revision_id,
      }),
    );
    expect(createRun).toHaveBeenCalledTimes(1);
    expect(getRun).toHaveBeenCalledWith('run-real-1');
    expect(listRunEvents).toHaveBeenCalled();
    expect(vi.mocked(listRunLogs).mock.calls.length).toBeGreaterThanOrEqual(2);
  });

  it('does not restore a stale snapshot when inputs change during validation', async () => {
    let resolveValidation!: (
      value: Awaited<ReturnType<typeof validateWorkflow>>,
    ) => void;
    vi.mocked(validateWorkflow).mockImplementation(
      () =>
        new Promise((resolve) => {
          resolveValidation = resolve;
        }),
    );
    const user = userEvent.setup();
    renderProductPath();
    await authorValidInputs(user);
    await user.click(screen.getByTestId('validate-draft-button'));
    await user.click(screen.getByRole('tab', { name: 'Samples' }));
    await user.type(screen.getAllByLabelText('Sample 1 sample')[0], '-changed');
    resolveValidation(validationResponse());

    await waitFor(() => expect(validateWorkflow).toHaveBeenCalledTimes(1));
    await user.click(screen.getByRole('tab', { name: 'Review' }));
    expect(screen.getByTestId('create-validated-run-button')).toBeDisabled();
    expect(
      await screen.findByText(/Inputs changed while validation was running/i),
    ).toBeVisible();
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
    const router = renderProductPath();
    await authorValidInputs(user);
    await user.click(screen.getByTestId('validate-draft-button'));
    await user.click(await screen.findByTestId('create-validated-run-button'));

    await waitFor(() =>
      expect(router.state.location.pathname).toBe('/runs/run-real-1'),
    );
    expect(await screen.findByRole('alert')).toHaveTextContent(
      'PREFLIGHT_FAILED: Preflight did not complete.',
    );
    expect(screen.getByTestId('run-status-badge')).toHaveTextContent('created');
    expect(screen.getByRole('button', { name: 'Run preflight' })).toBeEnabled();
  });
});
