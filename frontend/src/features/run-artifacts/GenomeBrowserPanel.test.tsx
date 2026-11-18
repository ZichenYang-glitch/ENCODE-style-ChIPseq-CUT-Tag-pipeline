import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { render, screen, waitFor } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { beforeEach, describe, expect, it, vi } from 'vitest';
import GenomeBrowserPanel from './GenomeBrowserPanel';
import { resetTrackReadGuardForTests } from './trackReadGuard';
import {
  getRunBrowserSession,
  prepareRunBrowserTrack,
  summarizeRunBrowserTrackRegion,
} from '../../api/generated/artifacts/artifacts';
import type { BrowserSessionResponse } from '../../api/generated/models';
import { ApiError } from '../../api/fetcher';

vi.mock('wuepgg', () => ({
  GenomeHub: ({
    genomeName,
    viewRegion,
    tracks,
  }: {
    genomeName: string;
    viewRegion?: string;
    tracks: Array<{ name?: string; url?: string }>;
  }) => {
    for (const track of tracks) {
      if (track.url) void fetch(track.url);
    }
    return (
      <div data-testid="washu-mock-viewer">
        {genomeName} {viewRegion} {tracks.map((track) => track.name).join(',')}
      </div>
    );
  },
}));

vi.mock('../../api/generated/artifacts/artifacts', () => ({
  getRunBrowserSession: vi.fn(),
  prepareRunBrowserTrack: vi.fn(),
  summarizeRunBrowserTrackRegion: vi.fn(),
}));

const sessionMock = vi.mocked(getRunBrowserSession);
const summaryMock = vi.mocked(summarizeRunBrowserTrackRegion);
vi.mocked(prepareRunBrowserTrack);

function session(
  overrides: Partial<BrowserSessionResponse> = {},
): BrowserSessionResponse {
  return {
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
        artifact_id: 'artifact-pet',
        name: 'sample_all.bedpe.gz',
        output_type: 'hitrac_bedpe_all',
        revision: `artifactrev-${'a'.repeat(64)}`,
        role: 'interaction',
        kind: 'pet_records',
        wash_u_type: 'longrange',
        label: 'PET records: sample_all.bedpe.gz',
        detail: 'Paired-end tags, not loops.',
        score_meaning: 'display_flag',
        source_form: 'display_derivative',
        preparation: 'ready',
        content_path: '/api/v1/runs/run-1/browser-tracks/artifact-pet/bytes/g/r',
        index_path: null,
        provenance: {
          suggested_region: 'chr1:200-2400',
          pair_count: 1,
          anchor_row_count: 2,
        },
      },
    ],
    unavailable: [],
    draw_record_limit: 100000,
    issues: [],
    ...overrides,
  };
}

function rangeHeader(init?: RequestInit): string | null {
  const headers = init?.headers;
  if (!headers || headers instanceof Headers || Array.isArray(headers)) {
    return new Headers(headers).get('range');
  }
  const entry = Object.entries(headers).find(([key]) => key.toLowerCase() === 'range');
  return typeof entry?.[1] === 'string' ? entry[1] : null;
}

function renderPanel() {
  const client = new QueryClient({
    defaultOptions: { queries: { retry: false } },
  });
  return render(
    <QueryClientProvider client={client}>
      <GenomeBrowserPanel runId="run-1" />
    </QueryClientProvider>,
  );
}

describe('GenomeBrowserPanel', () => {
  beforeEach(() => {
    resetTrackReadGuardForTests();
    sessionMock.mockReset();
    summaryMock.mockReset();
    summaryMock.mockResolvedValue({
      ok: true,
      kind: 'pet_records',
      records: 0,
      truncated: false,
      complete: true,
      detail: 'No records in this window. This is not a zero-signal measurement.',
      issues: [],
    });
    vi.stubGlobal(
      'fetch',
      vi.fn(async () => new Response(new Uint8Array([1]), { status: 206 })),
    );
  });

  it('shows an unbound reference instead of guessing an assembly', async () => {
    sessionMock.mockResolvedValue(
      session({
        reference: { status: 'unbound', viewer_scaffold: null },
        interaction_status: 'absent',
        tracks: [],
      }),
    );
    renderPanel();
    expect(
      await screen.findByText(/no reference profile binding/i),
    ).toBeInTheDocument();
    expect(screen.queryByTestId('washu-genome-viewer')).not.toBeInTheDocument();
    expect(screen.getByText(/no interaction artifact/i)).toBeInTheDocument();
  });

  it('keeps unsupported files and empty windows distinct from load failure', async () => {
    sessionMock.mockResolvedValue(
      session({
        unavailable: [
          {
            artifact_id: 'artifact-cool',
            name: 'contacts.mcool',
            output_type: 'mcool',
            reason: 'unsupported_contact_matrix',
            detail: 'Cooler matrices are not loaded.',
          },
        ],
      }),
    );
    renderPanel();
    expect(await screen.findByText(/Unsupported format: contacts.mcool/)).toBeInTheDocument();
    expect(screen.getByRole('checkbox', { name: /PET records/ })).not.toBeChecked();
    await userEvent.click(screen.getByRole('checkbox', { name: /PET records/ }));
    expect(await screen.findByText(/not a zero-signal measurement/)).toBeInTheDocument();
    expect(screen.getByTestId('washu-mock-viewer')).toHaveTextContent('hg38');
    expect(screen.getByTestId('washu-mock-viewer')).toHaveTextContent('chr1:200-2400');
  });

  it('shows a failure when a later byte read fails after preflight', async () => {
    vi.stubGlobal(
      'fetch',
      vi.fn(async (_input: RequestInfo | URL, init?: RequestInit) => {
        const range = new Headers(init?.headers).get('range');
        if (range === 'bytes=0-0') {
          return new Response(new Uint8Array([1]), { status: 206 });
        }
        return new Response('changed', { status: 409 });
      }),
    );
    sessionMock.mockResolvedValue(session());
    renderPanel();
    await userEvent.click(await screen.findByRole('checkbox', { name: /PET records/ }));
    expect(await screen.findByRole('alert')).toHaveTextContent('Track failed to load');
    await waitFor(() => {
      expect(screen.queryByTestId('washu-genome-viewer')).not.toBeInTheDocument();
    });
  });

  it('shows a failure when the index read fails after the content preflight', async () => {
    const value = session();
    value.tracks[0] = {
      ...value.tracks[0],
      index_path: `${value.tracks[0].content_path}.tbi`,
    };
    vi.stubGlobal(
      'fetch',
      vi.fn(async (input: RequestInfo | URL) => {
        const url = typeof input === 'string' ? input : input instanceof URL ? input.toString() : input.url;
        if (url.endsWith('.tbi')) return new Response('changed', { status: 409 });
        return new Response(new Uint8Array([1]), { status: 206 });
      }),
    );
    sessionMock.mockResolvedValue(value);
    renderPanel();
    await userEvent.click(await screen.findByRole('checkbox', { name: /PET records/ }));
    expect(await screen.findByRole('alert')).toHaveTextContent('Track failed to load (HTTP 409)');
    await waitFor(() => {
      expect(screen.queryByTestId('washu-genome-viewer')).not.toBeInTheDocument();
    });
  });

  it('does not mount a track whose byte request failed', async () => {
    vi.stubGlobal(
      'fetch',
      vi.fn(async () => new Response('no', { status: 409 })),
    );
    sessionMock.mockResolvedValue(session());
    renderPanel();
    await userEvent.click(await screen.findByRole('checkbox', { name: /PET records/ }));
    expect(await screen.findByRole('alert')).toHaveTextContent('Track failed to load');
    expect(screen.queryByTestId('washu-genome-viewer')).not.toBeInTheDocument();
  });

  it('removes the injected stylesheet when the panel leaves', async () => {
    sessionMock.mockResolvedValue(session());
    const view = renderPanel();
    await screen.findByText(/PET records: sample_all.bedpe.gz/);
    expect(
      document.head.querySelector('style[data-helixweave-washu="run-1"]'),
    ).not.toBeNull();
    localStorage.setItem('persist:helixweave-browser-run-1', 'stale');
    view.unmount();
    expect(
      document.head.querySelector('style[data-helixweave-washu="run-1"]'),
    ).toBeNull();
    expect(localStorage.getItem('persist:helixweave-browser-run-1')).toBeNull();
  });

  it('explains a file that is too large instead of drawing an empty track', async () => {
    const value = session();
    value.tracks[0] = {
      ...value.tracks[0],
      preparation: 'source_too_large',
      content_path: null,
    };
    sessionMock.mockResolvedValue(value);
    renderPanel();
    expect(await screen.findByText(/too large to index/i)).toBeInTheDocument();
    expect(screen.queryByTestId('washu-genome-viewer')).not.toBeInTheDocument();
    await userEvent.click(screen.getByRole('checkbox'));
    expect(screen.queryByTestId('washu-genome-viewer')).not.toBeInTheDocument();
  });

  it('does not keep a cancelled preflight as a failure after the track is selected again', async () => {
    let started = 0;
    vi.stubGlobal(
      'fetch',
      vi.fn(async (_input: RequestInfo | URL, init?: RequestInit) => {
        if (rangeHeader(init) === 'bytes=0-0' && ++started === 1) {
          return await new Promise<Response>((_resolve, reject) => {
            init?.signal?.addEventListener(
              'abort',
              () => reject(new DOMException('cancelled by component', 'AbortError')),
              { once: true },
            );
          });
        }
        return new Response(new Uint8Array([1]), { status: 206 });
      }),
    );
    sessionMock.mockResolvedValue(session());
    renderPanel();
    const checkbox = await screen.findByRole('checkbox');
    expect(checkbox).not.toBeChecked();
    await userEvent.click(checkbox);
    await waitFor(() => expect(started).toBe(1));
    await userEvent.click(checkbox);
    expect(checkbox).not.toBeChecked();
    await userEvent.click(checkbox);
    expect(checkbox).toBeChecked();
    await waitFor(() => expect(started).toBeGreaterThan(1));
    await waitFor(() => expect(screen.queryByRole('alert')).not.toBeInTheDocument());
    expect(screen.getByTestId('washu-genome-viewer')).toBeInTheDocument();
  });

  it('ignores a superseded preflight failure after a later read succeeds', async () => {
    let started = 0;
    let releaseLate: (response: Response) => void = () => {};
    vi.stubGlobal(
      'fetch',
      vi.fn(async (_input: RequestInfo | URL, init?: RequestInit) => {
        if (rangeHeader(init) === 'bytes=0-0' && ++started === 1) {
          return await new Promise<Response>((resolve) => {
            releaseLate = resolve;
          });
        }
        return new Response(new Uint8Array([1]), { status: 206 });
      }),
    );
    sessionMock.mockResolvedValue(session());
    renderPanel();
    const checkbox = await screen.findByRole('checkbox');
    await userEvent.click(checkbox);
    await waitFor(() => expect(started).toBe(1));
    await userEvent.click(checkbox);
    await userEvent.click(checkbox);
    await waitFor(() => expect(started).toBeGreaterThan(1));
    await waitFor(() => expect(screen.getByTestId('washu-genome-viewer')).toBeInTheDocument());
    releaseLate(new Response('changed', { status: 409 }));
    await waitFor(() => expect(screen.queryByRole('alert')).not.toBeInTheDocument());
    expect(screen.getByTestId('washu-genome-viewer')).toBeInTheDocument();
  });
});

function namedTrack(
  artifactId: string,
  label: string,
  overrides: Partial<BrowserSessionResponse['tracks'][number]> = {},
): BrowserSessionResponse['tracks'][number] {
  const base = session().tracks[0];
  return {
    ...base,
    artifact_id: artifactId,
    label,
    name: `${artifactId}.bw`,
    content_path: `/api/v1/runs/run-1/browser-tracks/${artifactId}/bytes/g/r`,
    index_path: null,
    kind: 'signal',
    role: 'signal',
    preparation: 'ready',
    provenance: null,
    ...overrides,
  };
}

function byteUrls(fetchMock: { mock: { calls: unknown[][] } }): string[] {
  return fetchMock.mock.calls
    .map((call) => String(call[0]))
    .filter((url) => url.includes('/bytes/'));
}

describe('GenomeBrowserPanel track selection', () => {
  beforeEach(() => {
    resetTrackReadGuardForTests();
    sessionMock.mockReset();
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
    vi.stubGlobal(
      'fetch',
      vi.fn(async () => new Response(new Uint8Array([1]), { status: 206 })),
    );
  });

  it('reads only the files chosen after an empty open', async () => {
    const fetchMock = vi.mocked(fetch);
    sessionMock.mockResolvedValue(
      session({
        interaction_status: 'absent',
        tracks: [
          namedTrack('artifact-a', 'Signal: A'),
          namedTrack('artifact-b', 'Signal: B'),
          namedTrack('artifact-c', 'Signal: C'),
        ],
      }),
    );
    renderPanel();
    expect(await screen.findByRole('checkbox', { name: /Signal: A/ })).not.toBeChecked();
    expect(screen.getByRole('checkbox', { name: /Signal: C/ })).toBeInTheDocument();
    expect(byteUrls(fetchMock)).toEqual([]);
    expect(summaryMock).not.toHaveBeenCalled();
    await userEvent.click(screen.getByRole('checkbox', { name: /Signal: A/ }));
    await userEvent.click(screen.getByRole('checkbox', { name: /Signal: B/ }));
    await waitFor(() => expect(byteUrls(fetchMock).length).toBeGreaterThan(0));
    const urls = byteUrls(fetchMock);
    expect(urls.some((url) => url.includes('/artifact-a/'))).toBe(true);
    expect(urls.some((url) => url.includes('/artifact-b/'))).toBe(true);
    expect(urls.some((url) => url.includes('/artifact-c/'))).toBe(false);
  });

  it('keeps the remaining choice across a same-generation refresh and stops after the last removal', async () => {
    const fetchMock = vi.mocked(fetch);
    const tracks = [
      namedTrack('artifact-a', 'Signal: A'),
      namedTrack('artifact-b', 'Signal: B'),
    ];
    sessionMock.mockResolvedValue(session({ interaction_status: 'absent', tracks }));
    renderPanel();
    await userEvent.click(await screen.findByRole('checkbox', { name: /Signal: A/ }));
    await userEvent.click(screen.getByRole('checkbox', { name: /Signal: B/ }));
    await waitFor(() => expect(byteUrls(fetchMock).length).toBeGreaterThan(0));
    await userEvent.click(screen.getByRole('checkbox', { name: /Signal: A/ }));
    const afterRemoval = byteUrls(fetchMock).length;
    sessionMock.mockResolvedValue(
      session({
        interaction_status: 'absent',
        tracks: tracks.map((track) => ({ ...track, detail: 'Refreshed list.' })),
      }),
    );
    await userEvent.click(screen.getByRole('button', { name: 'Refresh list' }));
    expect((await screen.findAllByText('Refreshed list.')).length).toBeGreaterThan(0);
    expect(screen.getByRole('checkbox', { name: /Signal: A/ })).not.toBeChecked();
    expect(screen.getByRole('checkbox', { name: /Signal: B/ })).toBeChecked();
    expect(byteUrls(fetchMock).some((url) => url.includes('/artifact-a/'))).toBe(true);
    const refreshedA = byteUrls(fetchMock).filter((url) => url.includes('/artifact-a/')).length;
    await userEvent.click(screen.getByRole('checkbox', { name: /Signal: B/ }));
    expect(screen.getByRole('checkbox', { name: /Signal: B/ })).not.toBeChecked();
    expect(screen.queryByTestId('washu-genome-viewer')).not.toBeInTheDocument();
    await waitFor(() => expect(byteUrls(fetchMock).filter((url) => url.includes('/artifact-a/')).length).toBe(refreshedA));
    expect(byteUrls(fetchMock).length).toBeGreaterThanOrEqual(afterRemoval);
    expect(
      byteUrls(fetchMock)
        .slice(afterRemoval)
        .some((url) => url.includes('/artifact-a/')),
    ).toBe(false);
  });

  it('does not restore an explicit choice after preparation or carry it to another run', async () => {
    const fetchMock = vi.mocked(fetch);
    const chosen = namedTrack('artifact-a', 'Signal: A');
    const waiting = namedTrack('artifact-b', 'Signal: B', {
      preparation: 'needed',
      content_path: null,
    });
    sessionMock.mockResolvedValue(
      session({ interaction_status: 'absent', tracks: [chosen, waiting] }),
    );
    const client = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    const view = render(
      <QueryClientProvider client={client}>
        <GenomeBrowserPanel runId="run-1" artifactIds={['artifact-a']} />
      </QueryClientProvider>,
    );
    const first = await screen.findByRole('checkbox', { name: /Signal: A/ });
    expect(first).toBeChecked();
    await userEvent.click(first);
    expect(first).not.toBeChecked();
    sessionMock.mockResolvedValue(
      session({
        interaction_status: 'absent',
        tracks: [
          chosen,
          {
            ...waiting,
            preparation: 'ready',
            content_path: '/api/v1/runs/run-1/browser-tracks/artifact-b/bytes/g/r',
          },
        ],
      }),
    );
    await userEvent.click(screen.getByRole('button', { name: 'Refresh list' }));
    await waitFor(() =>
      expect(screen.queryByText(/index is not prepared/i)).not.toBeInTheDocument(),
    );
    expect(screen.getByRole('checkbox', { name: /Signal: A/ })).not.toBeChecked();
    expect(screen.getByRole('checkbox', { name: /Signal: B/ })).not.toBeChecked();
    const beforeSwitch = byteUrls(fetchMock).length;
    sessionMock.mockResolvedValue(
      session({
        run_id: 'run-2',
        interaction_status: 'absent',
        tracks: [
          namedTrack('artifact-a', 'Signal: A', {
            content_path: '/api/v1/runs/run-2/browser-tracks/artifact-a/bytes/g/r',
          }),
        ],
      }),
    );
    view.rerender(
      <QueryClientProvider client={client}>
        <GenomeBrowserPanel runId="run-2" />
      </QueryClientProvider>,
    );
    const next = await screen.findByRole('checkbox', { name: /Signal: A/ });
    expect(next).not.toBeChecked();
    expect(byteUrls(fetchMock).filter((url) => url.includes('/runs/run-2/'))).toEqual([]);
    expect(byteUrls(fetchMock).length).toBe(beforeSwitch);
  });

  it('refreshes after a generation conflict and does not read the new files first', async () => {
    const fetchMock = vi.mocked(fetch);
    const original = session({
      interaction_status: 'absent',
      tracks: [namedTrack('artifact-a', 'Signal: A')],
    });
    sessionMock.mockResolvedValue(original);
    renderPanel();
    await userEvent.click(await screen.findByRole('checkbox', { name: /Signal: A/ }));
    await waitFor(() => expect(byteUrls(fetchMock).length).toBeGreaterThan(0));
    const beforeConflict = byteUrls(fetchMock).length;
    sessionMock.mockRejectedValue(
      new ApiError(
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
      ),
    );
    await userEvent.click(screen.getByRole('button', { name: 'Refresh list' }));
    expect(await screen.findByText(/not read until you choose them/i)).toBeInTheDocument();
    expect(byteUrls(fetchMock).filter((url) => url.includes('/new-gen/'))).toEqual([]);
    sessionMock.mockResolvedValue(
      session({
        interaction_status: 'absent',
        artifact_generation: `artifactgen-${'b'.repeat(64)}`,
        tracks: [
          namedTrack('artifact-a', 'Signal: A', {
            content_path: '/api/v1/runs/run-1/browser-tracks/artifact-a/bytes/new-gen/r',
          }),
        ],
      }),
    );
    await userEvent.click(screen.getByRole('button', { name: 'Refresh list' }));
    const box = await screen.findByRole('checkbox', { name: /Signal: A/ });
    expect(box).not.toBeChecked();
    expect(byteUrls(fetchMock).filter((url) => url.includes('/new-gen/'))).toEqual([]);
    expect(byteUrls(fetchMock).length).toBe(beforeConflict);
    await userEvent.click(box);
    await waitFor(() =>
      expect(byteUrls(fetchMock).some((url) => url.includes('/new-gen/'))).toBe(true),
    );
  });
});
