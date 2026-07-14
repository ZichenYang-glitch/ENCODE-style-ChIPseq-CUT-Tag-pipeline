import { describe, it, expect, vi } from 'vitest';
import { act, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { appRoutes } from '../../app/router';
import { renderWithRouter } from '../../test/test-utils';
import { createStubRunApiClient } from '../../api/runClient';
import { createStubWorkflowClient } from '../../api/client';
import type { RunApiClient } from '../../api/runClient';
import type { RunResponse } from '../../api/runTypes';
import type { WorkflowApiClient } from '../../api/client';
import { stubWorkflows, stubWorkflowSchemas } from '../../data/stubWorkflows';

const WORKFLOW_ID = 'encode-style-chipseq-cuttag-atac-mnase';

function deferred<T>() {
  let resolve!: (value: T) => void;
  const promise = new Promise<T>((resolvePromise) => {
    resolve = resolvePromise;
  });
  return { promise, resolve };
}

async function createRunActionRaceClient() {
  const baseClient = createStubRunApiClient();
  await baseClient.createRun(WORKFLOW_ID, {
    snapshot_id: 'vsnap_0123456789abcdef0123456789abcdef',
  });
  await baseClient.createRun(WORKFLOW_ID, {
    snapshot_id: 'vsnap_abcdef0123456789abcdef0123456789',
  });

  const delayedRefresh = deferred<RunResponse>();
  const delayedResponse = await baseClient.getRun('stub-run-1');
  let runOneGetCount = 0;
  const runClient: RunApiClient = {
    ...baseClient,
    getRun: vi.fn(async (runId) => {
      if (runId === 'stub-run-1') {
        runOneGetCount += 1;
        if (runOneGetCount === 2) {
          return delayedRefresh.promise;
        }
      }
      return baseClient.getRun(runId);
    }),
  };

  return { runClient, delayedRefresh, delayedResponse };
}

describe('Router', () => {
  it('redirects / to /workflows', async () => {
    const { router } = renderWithRouter(appRoutes, {
      initialEntries: ['/'],
    });

    await waitFor(() => {
      expect(router.state.location.pathname).toBe('/workflows');
    });
    expect(
      await screen.findByRole('heading', { name: /Workflows/i }),
    ).toBeInTheDocument();
  });

  it('renders the workflow catalog at /workflows', async () => {
    renderWithRouter(appRoutes, {
      initialEntries: ['/workflows'],
    });

    expect(
      await screen.findByText(/ENCODE-style ChIP-seq/i),
    ).toBeInTheDocument();
  });

  it('renders workflow detail when navigating to /workflows/:workflowId', async () => {
    renderWithRouter(appRoutes, {
      initialEntries: ['/workflows/encode-style-chipseq-cuttag-atac-mnase'],
    });

    expect(await screen.findByText(/Config schema/i)).toBeInTheDocument();
    expect(screen.getByText(/Sample schema/i)).toBeInTheDocument();
    expect(screen.getByText(/Options schema/i)).toBeInTheDocument();
  });

  it('places Author inputs before a default-closed developer schema disclosure', async () => {
    const user = userEvent.setup();
    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}`],
    });

    const authorInputs = await screen.findByRole('link', {
      name: 'Author inputs',
    });
    const disclosure = screen.getByTestId('developer-schema-details');

    expect(authorInputs).toHaveAttribute(
      'href',
      `/workflows/${WORKFLOW_ID}/new-run`,
    );
    expect(
      authorInputs.compareDocumentPosition(disclosure) &
        Node.DOCUMENT_POSITION_FOLLOWING,
    ).toBeTruthy();
    expect(disclosure).not.toHaveAttribute('open');

    await user.click(screen.getByText('Developer schema'));
    expect(disclosure).toHaveAttribute('open');
    expect(screen.getByText('Config schema')).toBeInTheDocument();
    expect(screen.getByText('Sample schema')).toBeInTheDocument();
    expect(screen.getByText('Options schema')).toBeInTheDocument();
  });

  it('navigates to workflow detail when a catalog item is selected', async () => {
    const user = userEvent.setup();
    const { router } = renderWithRouter(appRoutes, {
      initialEntries: ['/workflows'],
    });

    const workflowButton = await screen.findByText(/ENCODE-style ChIP-seq/i);
    await user.click(workflowButton);

    await waitFor(() => {
      expect(router.state.location.pathname).toBe(
        '/workflows/encode-style-chipseq-cuttag-atac-mnase',
      );
    });
    expect(await screen.findByText(/Config schema/i)).toBeInTheDocument();
  });

  it('resets validation draft when workflowId changes', async () => {
    const user = userEvent.setup();
    const secondWorkflowId = 'second-workflow';
    const workflowClient: WorkflowApiClient = {
      async listWorkflows() {
        return {
          ok: true,
          workflows: [
            ...stubWorkflows,
            {
              metadata: {
                workflow_id: secondWorkflowId,
                name: 'Second workflow',
                version: '0.1.0',
                description: 'Second workflow for route tests.',
                engines: ['snakemake'],
                tags: [],
              },
              capabilities: { supports: ['validation'] },
            },
          ],
          issues: [],
        };
      },
      async getWorkflowSchema(workflowId) {
        return {
          ok: true,
          workflow_id: workflowId,
          schema_hints: stubWorkflowSchemas['encode-style-chipseq-cuttag-atac-mnase'],
          issues: [],
        };
      },
      async validateWorkflow(workflowId, inputs) {
        return {
          ok: true,
          workflow_id: workflowId,
          value: { config: inputs.config, samples: [] },
          snapshot: {
            snapshot_id: 'vsnap_0123456789abcdef0123456789abcdef',
            workflow_id: workflowId,
            schema_version: '1.0.0',
            adapter_version: '0.3.0',
            payload_digest: 'a'.repeat(64),
            validated_at: '2026-07-14T00:00:00.000Z',
            expires_at: '2026-07-14T00:30:00.000Z',
          },
          issues: [],
        };
      },
    };

    const { router } = renderWithRouter(appRoutes, {
      initialEntries: ['/workflows/encode-style-chipseq-cuttag-atac-mnase'],
      clients: { workflowClient },
    });

    const samplesInput = await screen.findByLabelText(/Samples \(path string\)/i);
    await user.type(samplesInput, 'samples.tsv');
    expect(samplesInput).toHaveValue('samples.tsv');

    await user.click(
      await screen.findByRole('button', { name: /Second workflow/i }),
    );

    await waitFor(() => {
      expect(router.state.location.pathname).toBe(`/workflows/${secondWorkflowId}`);
    });

    const resetInput = screen.getByLabelText(/Samples \(path string\)/i);
    expect(resetInput).toHaveValue('');
  });

  it('loads an existing run at /runs/:runId', async () => {
    const runClient = createStubRunApiClient();
    await runClient.createRun(WORKFLOW_ID, {
      snapshot_id: 'vsnap_0123456789abcdef0123456789abcdef',
    });

    renderWithRouter(appRoutes, {
      initialEntries: ['/runs/stub-run-1'],
      clients: { runClient },
    });

    expect(await screen.findByTestId('run-status-badge')).toHaveTextContent(
      'created',
    );
    expect(screen.getByTestId('run-event-feed')).toBeInTheDocument();
  });

  it('ignores a stale manual refresh after navigating to another run', async () => {
    const user = userEvent.setup();
    const { runClient, delayedRefresh, delayedResponse } =
      await createRunActionRaceClient();
    const { router } = renderWithRouter(appRoutes, {
      initialEntries: ['/runs/stub-run-1'],
      clients: { runClient },
    });

    expect(await screen.findByText('stub-run-1')).toBeInTheDocument();
    await user.click(screen.getByTestId('refresh-run-button'));
    await waitFor(() => expect(runClient.getRun).toHaveBeenCalledTimes(2));

    await act(async () => {
      await router.navigate('/runs/stub-run-2');
    });
    expect(await screen.findByText('stub-run-2')).toBeInTheDocument();

    await act(async () => {
      delayedRefresh.resolve(delayedResponse);
      await delayedRefresh.promise;
    });

    expect(screen.getByText('stub-run-2')).toBeInTheDocument();
    expect(screen.queryByText('stub-run-1')).not.toBeInTheDocument();
  });

  it('ignores a stale post-cancel refresh after navigating to another run', async () => {
    const user = userEvent.setup();
    const { runClient, delayedRefresh, delayedResponse } =
      await createRunActionRaceClient();
    const { router } = renderWithRouter(appRoutes, {
      initialEntries: ['/runs/stub-run-1'],
      clients: { runClient },
    });

    expect(await screen.findByText('stub-run-1')).toBeInTheDocument();
    await user.click(screen.getByTestId('cancel-run-button'));
    await waitFor(() => expect(runClient.getRun).toHaveBeenCalledTimes(2));

    await act(async () => {
      await router.navigate('/runs/stub-run-2');
    });
    expect(await screen.findByText('stub-run-2')).toBeInTheDocument();

    await act(async () => {
      delayedRefresh.resolve(delayedResponse);
      await delayedRefresh.promise;
    });

    expect(screen.getByText('stub-run-2')).toBeInTheDocument();
    expect(screen.queryByText('stub-run-1')).not.toBeInTheDocument();
  });

  it('renders a 404 page for unknown routes', async () => {
    renderWithRouter(appRoutes, {
      initialEntries: ['/unknown-route'],
    });

    expect(await screen.findByText(/404/i)).toBeInTheDocument();
    expect(screen.getByRole('link', { name: /Back to workflows/i })).toBeInTheDocument();
  });

  it('renders the app without query client errors', async () => {
    renderWithRouter(appRoutes, {
      initialEntries: ['/workflows'],
    });

    expect(
      await screen.findByRole('heading', { name: 'HelixWeave' }),
    ).toBeInTheDocument();
  });

  it('renders not-found state for unknown workflowId', async () => {
    renderWithRouter(appRoutes, {
      initialEntries: ['/workflows/missing-id'],
    });

    expect(await screen.findByText(/Workflow not found/i)).toBeInTheDocument();
    expect(screen.getByText(/missing-id/i)).toBeInTheDocument();
    expect(screen.queryByLabelText(/Samples \(path string\)/i)).not.toBeInTheDocument();
    expect(screen.queryByTestId('create-run-button')).not.toBeInTheDocument();
    expect(screen.queryByText(/Validation Assistant/i)).not.toBeInTheDocument();
  });

  it('does not report workflow API failures as not found', async () => {
    const workflowClient = createStubWorkflowClient();
    workflowClient.listWorkflows = vi.fn().mockResolvedValue({
      ok: false,
      workflows: [],
      issues: [
        {
          code: 'API_UNAVAILABLE',
          message: 'Workflow service is unavailable.',
          severity: 'error',
          path: null,
          source: 'api',
          technical_message: null,
          hint: null,
          context: {},
        },
      ],
    });

    renderWithRouter(appRoutes, {
      initialEntries: [`/workflows/${WORKFLOW_ID}`],
      clients: { workflowClient },
    });

    expect(
      await screen.findByRole('heading', { name: /Workflow unavailable/i }),
    ).toBeInTheDocument();
    expect(screen.queryByText(/Workflow not found/i)).not.toBeInTheDocument();
  });

  it('renders the 404 back link as a single interactive element', async () => {
    renderWithRouter(appRoutes, {
      initialEntries: ['/unknown-route'],
    });

    const link = await screen.findByRole('link', { name: /Back to workflows/i });
    expect(link).toBeInTheDocument();
    expect(link.tagName.toLowerCase()).toBe('a');
    expect(link.querySelector('button')).toBeNull();
    expect(
      screen.getAllByRole('link', { name: /Back to workflows/i }),
    ).toHaveLength(1);
  });
});
