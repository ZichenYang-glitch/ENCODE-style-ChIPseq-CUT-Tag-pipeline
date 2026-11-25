import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { render, screen, waitFor, within } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import type { ReactNode } from 'react';
import { beforeEach, describe, expect, it, vi } from 'vitest';
import GenomeBrowserPanel from './GenomeBrowserPanel';
import { resetTrackReadGuardForTests } from './trackReadGuard';
import { getRunBrowserSession, summarizeRunBrowserTrackRegion } from '../../api/generated/artifacts/artifacts';
import { listRuns } from '../../api/generated/runs/runs';
import type { BrowserSessionResponse, RunSummaryResponse } from '../../api/generated/models';
import { ApiError } from '../../api/fetcher';

vi.mock('wuepgg', () => ({
  GenomeHub: ({
    genomeName,
    viewRegion,
    tracks,
  }: {
    genomeName: string;
    viewRegion?: string;
    tracks: Array<{ type?: string; name?: string; url?: string; fileObj?: unknown; options?: unknown }>;
  }) => {
    (globalThis as { __hubTracks?: unknown }).__hubTracks = tracks;
    for (const track of tracks) {
      if (track.url) void fetch(track.url);
    }
    return (
      <div data-testid="washu-mock-viewer">
        <span>{genomeName}</span>
        <span>{viewRegion}</span>
        <pre data-testid="washu-track-models">{JSON.stringify(tracks)}</pre>
      </div>
    );
  },
}));

vi.mock('../../api/generated/artifacts/artifacts', () => ({
  getRunBrowserSession: vi.fn(),
  prepareRunBrowserTrack: vi.fn(),
  summarizeRunBrowserTrackRegion: vi.fn(),
}));

vi.mock('../../api/generated/runs/runs', () => ({
  listRuns: vi.fn(),
}));

vi.mock('../../viewer/local/LocalTrackPicker', () => ({
  default: ({
    assembly,
    onTracksChange,
  }: {
    assembly: string;
    onTracksChange: (tracks: Array<{
      id: string;
      label: string;
      assembly: string;
      model: { type: string; name: string; [key: string]: unknown };
    }>) => void;
  }) => (
    <div data-testid="local-track-picker" data-assembly={assembly}>
      <button
        type="button"
        onClick={() =>
          onTracksChange([
            {
              id: 'local-signal',
              label: 'Local signal',
              assembly: 'GRCh38',
              model: {
                type: 'bigwig',
                name: 'Local signal',
                url: '',
                fileObj: new File(['a'], 'signal.bw'),
                metadata: { genome: 'hg38' },
                options: { local: true },
              },
            },
            {
              id: 'local-foreign',
              label: 'Foreign signal',
              assembly: 'GRCh38',
              model: {
                type: 'bigwig',
                name: 'Foreign signal',
                url: 'https://epigenome.wustl.edu/proxy/signal.bw',
              },
            },
          ])
        }
      >
        Use local signal
      </button>
      <button
        type="button"
        onClick={() =>
          onTracksChange([
            {
              id: 'local-mouse',
              label: 'Mouse signal',
              assembly: 'mm10',
              model: { type: 'bigwig', name: 'Mouse signal' },
            },
          ])
        }
      >
        Use mouse signal
      </button>
    </div>
  ),
}));

const sessionMock = vi.mocked(getRunBrowserSession);
const listRunsMock = vi.mocked(listRuns);
const summaryMock = vi.mocked(summarizeRunBrowserTrackRegion);

function runSummary(id: string): RunSummaryResponse {
  return {
    run_id: id,
    workflow_id: 'hitrac-preprocess',
    status: 'succeeded',
    created_at: '2026-01-01T00:00:00Z',
    updated_at: '2026-01-02T00:00:00Z',
    started_at: null,
    ended_at: null,
    current_stage: null,
  };
}

function session(
  runId: string,
  identity: string,
  overrides: Partial<BrowserSessionResponse> = {},
): BrowserSessionResponse {
  return {
    ok: true,
    run_id: runId,
    workflow_id: 'hitrac-preprocess',
    artifact_generation: `artifactgen-${'a'.repeat(64)}`,
    reference: {
      status: 'bound',
      profile_id: 'refp_' + 'ab'.repeat(16),
      revision_id: 'refpr_' + 'cd'.repeat(16),
      identity_sha256: identity,
      binding_digest: 'f'.repeat(64),
      assembly_label: 'hg38',
      viewer_scaffold: 'hg38',
    },
    annotation_status: 'absent',
    interaction_status: 'absent',
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
        content_path: `/api/v1/runs/${runId}/browser-tracks/artifact-signal/bytes/g/r`,
        index_path: null,
        provenance: { suggested_region: 'chr1:200-2400' },
      },
    ],
    unavailable: [],
    draw_record_limit: 100000,
    issues: [],
    ...overrides,
  };
}

function renderPanel(ui: ReactNode) {
  const client = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(<QueryClientProvider client={client}>{ui}</QueryClientProvider>);
}

function byteUrls(fetchMock: { mock: { calls: unknown[][] } }): string[] {
  return fetchMock.mock.calls
    .map((call) => String(call[0]))
    .filter((url) => url.includes('/bytes/'));
}

function models(): Array<{ type?: string; name?: string; url?: string }> {
  return JSON.parse(screen.getByTestId('washu-track-models').textContent ?? '[]') as Array<{
    type?: string;
    name?: string;
    url?: string;
  }>;
}

describe('GenomeBrowserPanel overlay', () => {
  beforeEach(() => {
    resetTrackReadGuardForTests();
    sessionMock.mockReset();
    listRunsMock.mockReset();
    summaryMock.mockReset();
    summaryMock.mockResolvedValue({
      ok: true,
      kind: 'signal',
      records: null,
      truncated: false,
      complete: true,
      detail: 'Signal is not counted as records.',
      issues: [],
    });
    listRunsMock.mockResolvedValue({
      ok: true,
      runs: [runSummary('run-1'), runSummary('run-2')],
      next_cursor: null,
      issues: [],
    });
    vi.stubGlobal(
      'fetch',
      vi.fn(async () => new Response(new Uint8Array([1]), { status: 206 })),
    );
  });

  it('opens with no run on hg38 and draws a local model without calling a run session', async () => {
    renderPanel(<GenomeBrowserPanel />);
    expect(await screen.findByRole('combobox', { name: 'Assembly' })).toHaveValue('hg38');
    expect(
      (screen.getByRole('combobox', { name: 'Assembly' }) as HTMLSelectElement).selectedOptions[0]
        ?.textContent,
    ).toBe('hg38 · GRCh38');
    expect(screen.getByText('当前声明 GRCh38')).toBeInTheDocument();
    expect(screen.getByText('用户指定 GRCh38')).toBeInTheDocument();
    expect(sessionMock).not.toHaveBeenCalled();
    await userEvent.selectOptions(screen.getByRole('combobox', { name: 'Assembly' }), 'mm10');
    expect(screen.getByTestId('local-track-picker')).toHaveAttribute('data-assembly', 'mm10');
    await userEvent.clear(screen.getByLabelText('Genomic region'));
    await userEvent.type(screen.getByLabelText('Genomic region'), 'chr1:1-2000');
    await userEvent.click(screen.getByRole('button', { name: 'Use local signal' }));
    await userEvent.selectOptions(screen.getByRole('combobox', { name: 'Assembly' }), 'hg38');
    await waitFor(() => expect(screen.getByTestId('washu-genome-viewer')).toBeInTheDocument());
    const drawn = models();
    expect(drawn.map((track) => track.name)).toEqual(['Ruler', 'Local signal']);
    const hubTracks = (globalThis as { __hubTracks?: Array<{ name?: string; url?: string; fileObj?: unknown }> })
      .__hubTracks;
    const local = hubTracks?.find((track) => track.name === 'Local signal');
    expect(local?.url).toBe('');
    expect(local?.fileObj).toBeInstanceOf(File);
    expect(drawn.some((track) => String(track.url).includes('epigenome.wustl.edu'))).toBe(false);
    expect(screen.getByTestId('reference-conflict')).toHaveTextContent(/another website/i);
    expect(screen.getByTestId('local-browser-track')).toHaveAttribute('data-source', 'local');
    expect(screen.getByTestId('local-browser-track')).toHaveAttribute('data-assembly', 'GRCh38');
    await userEvent.selectOptions(screen.getByRole('combobox', { name: 'Assembly' }), 'mm10');
    expect(screen.getByText('当前声明 GRCm38')).toBeInTheDocument();
    const kept = screen.getAllByTestId('reference-conflict').find((node) =>
      node.getAttribute('data-assembly') === 'GRCh38' && /not relabeled/i.test(node.textContent ?? ''),
    );
    expect(kept).toBeTruthy();
    expect(screen.queryByTestId('washu-genome-viewer')).not.toBeInTheDocument();
  });

  it('adds an authorized run by name and leaves its tracks unread until one is chosen', async () => {
    const fetchMock = vi.mocked(fetch);
    sessionMock.mockImplementation(async (id: string) =>
      session(id, 'e'.repeat(64)),
    );
    renderPanel(<GenomeBrowserPanel runId="run-1" />);
    const opened = await screen.findByRole('article', { name: 'Run run-1' });
    expect(within(opened).getByRole('checkbox', { name: /Signal: A/ })).not.toBeChecked();
    expect(within(opened).getByTestId('browser-track')).toHaveAttribute('data-source', 'platform');
    expect(within(opened).getByTestId('browser-track')).toHaveAttribute(
      'data-revision',
      `artifactrev-${'a'.repeat(64)}`,
    );
    expect(byteUrls(fetchMock)).toEqual([]);
    expect(sessionMock.mock.calls.some((call) => call[0] === 'run-2')).toBe(false);

    await userEvent.selectOptions(screen.getByRole('combobox', { name: 'Authorized runs' }), 'run-2');
    await userEvent.click(screen.getByRole('button', { name: 'Add run' }));
    const added = await screen.findByRole('article', { name: 'Run run-2' });
    expect(within(added).getByRole('checkbox', { name: /Signal: A/ })).not.toBeChecked();
    expect(byteUrls(fetchMock).some((url) => url.includes('/runs/run-2/'))).toBe(false);
    await userEvent.click(within(added).getByRole('checkbox', { name: /Signal: A/ }));
    await waitFor(() =>
      expect(byteUrls(fetchMock).some((url) => url.includes('/runs/run-2/'))).toBe(true),
    );
    expect(byteUrls(fetchMock).some((url) => url.includes('/runs/run-1/'))).toBe(false);
  });

  it('keeps a rejected run from removing another run that shares its artifact id', async () => {
    const fetchMock = vi.fn(async (input: RequestInfo | URL) => {
      const url = String(input);
      if (url.includes('/runs/run-1/')) return new Response('no', { status: 409 });
      return new Response(new Uint8Array([1]), { status: 206 });
    });
    vi.stubGlobal('fetch', fetchMock);
    sessionMock.mockImplementation(async (id: string) =>
      session(id, id === 'run-1' ? 'e'.repeat(64) : 'f'.repeat(64)),
    );
    renderPanel(<GenomeBrowserPanel runId="run-1" />);
    const opened = await screen.findByRole('article', { name: 'Run run-1' });
    await userEvent.click(within(opened).getByRole('checkbox', { name: /Signal: A/ }));
    await userEvent.selectOptions(screen.getByRole('combobox', { name: 'Authorized runs' }), 'run-2');
    await userEvent.click(screen.getByRole('button', { name: 'Add run' }));
    const added = await screen.findByRole('article', { name: 'Run run-2' });
    await userEvent.click(within(added).getByRole('checkbox', { name: /Signal: A/ }));
    expect(await screen.findByRole('alert')).toHaveTextContent(/Signal: A \(run-1\)/);
    await waitFor(() => expect(screen.getByTestId('washu-genome-viewer')).toBeInTheDocument());
    expect(models().some((track) => track.url?.includes('/runs/run-2/'))).toBe(true);
    expect(models().some((track) => track.url?.includes('/runs/run-1/'))).toBe(false);
    expect(byteUrls(fetchMock).some((url) => url.includes('/runs/run-2/'))).toBe(true);
  });

  it('points out an authoritative binding conflict and still draws the matching run', async () => {
    const fetchMock = vi.mocked(fetch);
    sessionMock.mockImplementation(async (id: string) =>
      session(id, id === 'run-1' ? 'e'.repeat(64) : 'f'.repeat(64)),
    );
    renderPanel(<GenomeBrowserPanel runId="run-1" />);
    const opened = await screen.findByRole('article', { name: 'Run run-1' });
    await userEvent.click(within(opened).getByRole('checkbox', { name: /Signal: A/ }));
    await waitFor(() => expect(screen.getByTestId('washu-genome-viewer')).toBeInTheDocument());
    await userEvent.selectOptions(screen.getByRole('combobox', { name: 'Authorized runs' }), 'run-2');
    await userEvent.click(screen.getByRole('button', { name: 'Add run' }));
    const added = await screen.findByRole('article', { name: 'Run run-2' });
    expect(within(added).getByTestId('reference-note')).toHaveTextContent(
      /different reference revisions/i,
    );
    const addedBox = within(added).getByRole('checkbox', { name: /Signal: A/ });
    expect(addedBox).toBeEnabled();
    await userEvent.click(addedBox);
    await waitFor(() =>
      expect(byteUrls(fetchMock).some((url) => url.includes('/runs/run-2/'))).toBe(true),
    );
    expect(screen.getByTestId('washu-genome-viewer')).toBeInTheDocument();
    expect(models().some((track) => track.url?.includes('/runs/run-2/'))).toBe(true);
    await userEvent.click(screen.getByRole('button', { name: 'Use local signal' }));
    expect(await screen.findByTestId('local-browser-track')).toHaveTextContent('用户指定 GRCh38');
    expect(models().some((track) => track.name === 'Local signal')).toBe(true);
    expect(models().some((track) => track.url?.includes('/runs/run-1/'))).toBe(true);
  });

  it('does not relabel a local assembly that disagrees with the window', async () => {
    sessionMock.mockResolvedValue(session('run-1', 'e'.repeat(64)));
    renderPanel(<GenomeBrowserPanel runId="run-1" />);
    const opened = await screen.findByRole('article', { name: 'Run run-1' });
    await userEvent.click(within(opened).getByRole('checkbox', { name: /Signal: A/ }));
    await userEvent.click(screen.getByRole('button', { name: 'Use mouse signal' }));
    expect(await screen.findByTestId('reference-conflict')).toHaveTextContent(/not relabeled/i);
    await waitFor(() => expect(screen.getByTestId('washu-genome-viewer')).toBeInTheDocument());
    expect(models().some((track) => track.name === 'Mouse signal')).toBe(false);
    expect(models().some((track) => track.url?.includes('/runs/run-1/'))).toBe(true);
  });

  it('keeps the second run and local files when the opened run fails', async () => {
    sessionMock.mockImplementation(async (id: string) => session(id, 'e'.repeat(64)));
    renderPanel(<GenomeBrowserPanel runId="run-1" />);
    const opened = await screen.findByRole('article', { name: 'Run run-1' });
    await userEvent.click(within(opened).getByRole('checkbox', { name: /Signal: A/ }));
    await userEvent.selectOptions(screen.getByRole('combobox', { name: 'Authorized runs' }), 'run-2');
    await userEvent.click(screen.getByRole('button', { name: 'Add run' }));
    const added = await screen.findByRole('article', { name: 'Run run-2' });
    await userEvent.click(within(added).getByRole('checkbox', { name: /Signal: A/ }));
    await userEvent.click(screen.getByRole('button', { name: 'Use local signal' }));
    await waitFor(() => {
      expect(models().some((track) => track.url?.includes('/runs/run-1/'))).toBe(true);
      expect(models().some((track) => track.url?.includes('/runs/run-2/'))).toBe(true);
      expect(models().some((track) => track.name === 'Local signal')).toBe(true);
    });

    sessionMock.mockImplementation(async (id: string) => {
      if (id === 'run-1') {
        throw new ApiError(
          409,
          'BROWSER_TRACK_READ_CONFLICT',
          'Artifact generation changed. Reload the browser view.',
          [
            {
              code: 'BROWSER_TRACK_READ_CONFLICT',
              message: 'Artifact generation changed. Reload the browser view.',
              path: 'generation',
            },
          ],
        );
      }
      return session(id, 'f'.repeat(64));
    });
    await userEvent.click(within(opened).getByRole('button', { name: 'Refresh list' }));
    expect(await screen.findByText(/not read until you choose them/i)).toBeInTheDocument();
    await waitFor(() => {
      const drawn = models();
      expect(drawn.some((track) => track.url?.includes('/runs/run-2/'))).toBe(true);
      expect(drawn.some((track) => track.name === 'Local signal')).toBe(true);
      expect(drawn.some((track) => String(track.url).includes('/runs/run-1/'))).toBe(false);
    });
    expect(screen.getByTestId('local-track-picker')).toBeInTheDocument();
    expect(screen.getByRole('combobox', { name: 'Authorized runs' })).toBeInTheDocument();
    expect(
      within(screen.getByRole('article', { name: 'Run run-1' })).queryByRole('checkbox'),
    ).not.toBeInTheDocument();

    sessionMock.mockImplementation(async (id: string) => {
      if (id === 'run-1') {
        throw new ApiError(403, 'FORBIDDEN', 'This run is no longer authorized.');
      }
      return session(id, 'f'.repeat(64));
    });
    await userEvent.click(screen.getByRole('button', { name: 'Refresh list' }));
    expect(await screen.findByRole('button', { name: 'Retry browser view' })).toBeInTheDocument();
    expect(screen.getByRole('alert')).toHaveTextContent(/no longer authorized/i);
    const drawn = models();
    expect(drawn.some((track) => track.url?.includes('/runs/run-2/'))).toBe(true);
    expect(drawn.some((track) => track.name === 'Local signal')).toBe(true);
    expect(drawn.some((track) => String(track.url).includes('/runs/run-1/'))).toBe(false);
    expect(screen.getByRole('article', { name: 'Run run-2' })).toBeInTheDocument();
  });

  it('leaves the open run drawn when another authorized run is refused', async () => {
    sessionMock.mockImplementation(async (id: string) => {
      if (id === 'run-2') {
        throw new ApiError(503, 'BROWSER_TRACK_UNAVAILABLE', 'Run 2 is not available.');
      }
      return session(id, 'e'.repeat(64));
    });
    renderPanel(<GenomeBrowserPanel runId="run-1" />);
    const opened = await screen.findByRole('article', { name: 'Run run-1' });
    await userEvent.click(within(opened).getByRole('checkbox', { name: /Signal: A/ }));
    await waitFor(() => expect(screen.getByTestId('washu-genome-viewer')).toBeInTheDocument());
    await userEvent.selectOptions(screen.getByRole('combobox', { name: 'Authorized runs' }), 'run-2');
    await userEvent.click(screen.getByRole('button', { name: 'Add run' }));
    expect(await screen.findByText(/Run 2 is not available/)).toBeInTheDocument();
    expect(screen.getByText(/other runs stay available/i)).toBeInTheDocument();
    expect(screen.getByTestId('washu-genome-viewer')).toBeInTheDocument();
    expect(screen.getByRole('button', { name: 'Remove run run-2' })).toBeInTheDocument();
  });
});
