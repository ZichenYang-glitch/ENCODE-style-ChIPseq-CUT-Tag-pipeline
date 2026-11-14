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
});
