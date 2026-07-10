import { describe, it, expect } from 'vitest';
import { screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { appRoutes } from '../../app/router';
import { renderWithRouter } from '../../test/test-utils';
import { createStubRunApiClient } from '../../api/runClient';
import type { WorkflowApiClient } from '../../api/client';
import { stubWorkflows, stubWorkflowSchemas } from '../../data/stubWorkflows';

const WORKFLOW_ID = 'encode-style-chipseq-cuttag-atac-mnase';

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
      config: { genome: 'hg38' },
      samples: [{ name: 'sample-1', fastq_r1: 's1_R1.fq.gz' }],
      options: { threads: 4 },
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
      await screen.findByRole('heading', { name: /Workflow Platform/i }),
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

  it('renders the 404 back link as a single interactive element', async () => {
    renderWithRouter(appRoutes, {
      initialEntries: ['/unknown-route'],
    });

    const link = await screen.findByRole('link', { name: /Back to workflows/i });
    expect(link).toBeInTheDocument();
    expect(link.tagName.toLowerCase()).toBe('a');
    expect(link.querySelector('button')).toBeNull();
    expect(screen.getAllByRole('link').length).toBe(1);
  });
});
