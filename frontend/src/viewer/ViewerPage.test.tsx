import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { render, screen } from '@testing-library/react';
import { afterEach, describe, expect, it, vi } from 'vitest';
import { ViewerPage } from './ViewerPage';
import { getRunBrowserSession } from '../api/generated/artifacts/artifacts';

vi.mock('wuepgg', () => ({
  GenomeHub: () => <div data-testid="washu-mock-viewer" />,
}));

vi.mock('../api/generated/artifacts/artifacts', () => ({
  getRunBrowserSession: vi.fn(),
  prepareRunBrowserTrack: vi.fn(),
  summarizeRunBrowserTrackRegion: vi.fn(),
}));

vi.mock('../api/generated/runs/runs', () => ({
  listRuns: vi.fn(async () => ({
    ok: true,
    runs: [],
    next_cursor: null,
    issues: [],
  })),
}));

const sessionMock = vi.mocked(getRunBrowserSession);

function renderPage() {
  const client = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={client}>
      <ViewerPage />
    </QueryClientProvider>,
  );
}

describe('ViewerPage', () => {
  afterEach(() => {
    window.history.replaceState(null, '', '/');
    sessionMock.mockReset();
  });

  it('opens a standalone window with no run on hg38', async () => {
    window.history.replaceState(null, '', '/viewer.html');
    renderPage();
    expect(await screen.findByRole('combobox', { name: 'Assembly' })).toHaveValue('hg38');
    expect(screen.getByLabelText('Genomic region')).toBeInTheDocument();
    expect(screen.queryByText(/Open this page with a run id/i)).not.toBeInTheDocument();
    expect(sessionMock).not.toHaveBeenCalled();
  });

  it('opens the run named in the standalone address and leaves an invalid id unloaded', async () => {
    sessionMock.mockResolvedValue({
      ok: true,
      run_id: 'run-1',
      workflow_id: 'hitrac-preprocess',
      artifact_generation: `artifactgen-${'a'.repeat(64)}`,
      reference: {
        status: 'bound',
        profile_id: 'refp_' + 'ab'.repeat(16),
        revision_id: 'refpr_' + 'cd'.repeat(16),
        identity_sha256: 'e'.repeat(64),
        binding_digest: 'f'.repeat(64),
        assembly_label: 'hg38',
        viewer_scaffold: 'hg38',
      },
      annotation_status: 'absent',
      interaction_status: 'available',
      tracks: [
        {
          artifact_id: 'artifact-signal',
          name: 'signal.bw',
          output_type: 'bigwig',
          revision: `artifactrev-${'a'.repeat(64)}`,
          role: 'signal',
          kind: 'signal',
          wash_u_type: 'bigwig',
          label: 'Signal: A',
          detail: 'Signal track.',
          score_meaning: 'signal',
          source_form: 'display_derivative',
          preparation: 'ready',
          content_path: '/api/v1/runs/run-1/browser-tracks/artifact-signal/bytes/g/r',
          index_path: null,
          provenance: null,
        },
      ],
      unavailable: [],
      draw_record_limit: 100000,
      issues: [],
    });
    window.history.replaceState(null, '', '/viewer.html?run=run-1');
    renderPage();
    expect(await screen.findByTestId('browser-run')).toHaveAttribute('data-run-id', 'run-1');
    expect(screen.getByRole('checkbox', { name: /Signal: A/ })).not.toBeChecked();
  });

  it('keeps a framed document waiting for the host', () => {
    const parent = { postMessage: () => undefined };
    const descriptor = Object.getOwnPropertyDescriptor(window, 'parent');
    Object.defineProperty(window, 'parent', { configurable: true, get: () => parent });
    window.history.replaceState(null, '', '/viewer.html#instance=11111111-1111-4111-8111-111111111111');
    const view = renderPage();
    try {
      expect(screen.getByRole('status')).toHaveTextContent('Waiting for the platform viewer host.');
      expect(screen.queryByRole('combobox', { name: 'Assembly' })).not.toBeInTheDocument();
      expect(sessionMock).not.toHaveBeenCalled();
    } finally {
      view.unmount();
      if (descriptor) Object.defineProperty(window, 'parent', descriptor);
    }
  });

  it('does not load a run id the address bar cannot name', async () => {
    window.history.replaceState(null, '', '/viewer.html?run=../secret');
    renderPage();
    expect(await screen.findByRole('combobox', { name: 'Assembly' })).toHaveValue('hg38');
    expect(sessionMock).not.toHaveBeenCalled();
  });
});
