import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { fireEvent, render, screen, waitFor } from '@testing-library/react';
import { useEffect, useRef } from 'react';
import { beforeEach, describe, expect, it, vi } from 'vitest';
import GenomeBrowserPanel from './GenomeBrowserPanel';
import { resetTrackReadGuardForTests } from './trackReadGuard';
import {
  getRunBrowserSession,
  summarizeRunBrowserTrackRegion,
} from '../../api/generated/artifacts/artifacts';
import type { BrowserSessionResponse } from '../../api/generated/models';

const hub = vi.hoisted(() => ({
  mounts: [] as Array<{ region: string; storeId: string }>,
  retired: [] as Array<{ region: string; storeId: string }>,
}));

vi.mock('wuepgg', () => ({
  GenomeHub: ({
    genomeName,
    viewRegion,
    storeConfig,
  }: {
    genomeName: string;
    viewRegion?: string;
    storeConfig?: { storeId?: string };
  }) => {
    const mountedRegion = useRef(viewRegion ?? '');
    const storeId = storeConfig?.storeId ?? '';
    useEffect(() => {
      hub.mounts.push({ region: mountedRegion.current, storeId });
      return () => {
        hub.retired.push({ region: mountedRegion.current, storeId });
      };
    }, [storeId]);
    return (
      <div data-testid="washu-mock-viewer">
        {genomeName} {mountedRegion.current}
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
  listRuns: vi.fn(async () => ({
    ok: true,
    runs: [],
    next_cursor: null,
    issues: [],
  })),
}));

const sessionMock = vi.mocked(getRunBrowserSession);
const summaryMock = vi.mocked(summarizeRunBrowserTrackRegion);

function track(
  artifactId: string,
  label: string,
): BrowserSessionResponse['tracks'][number] {
  return {
    artifact_id: artifactId,
    name: `${artifactId}.bw`,
    output_type: 'bigwig',
    revision: `artifactrev-${'a'.repeat(64)}`,
    role: 'signal',
    kind: 'signal',
    wash_u_type: 'bigwig',
    label,
    detail: 'Signal track.',
    score_meaning: 'display_flag',
    source_form: 'display_derivative',
    preparation: 'ready',
    content_path: `/api/v1/runs/run-1/browser-tracks/${artifactId}/bytes/g/r`,
    index_path: null,
    provenance: { suggested_region: 'chr1:4500-6700' },
  };
}

function session(): BrowserSessionResponse {
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
    interaction_status: 'absent',
    tracks: [track('artifact-a', 'Signal: A'), track('artifact-b', 'Signal: B')],
    unavailable: [],
    draw_record_limit: 100000,
    issues: [],
  };
}

function renderPanel() {
  const client = new QueryClient({
    defaultOptions: { queries: { retry: false } },
  });
  const view = render(
    <QueryClientProvider client={client}>
      <GenomeBrowserPanel runId="run-1" />
    </QueryClientProvider>,
  );
  return {
    ...view,
    rerenderPanel() {
      view.rerender(
        <QueryClientProvider client={client}>
          <GenomeBrowserPanel runId="run-1" />
        </QueryClientProvider>,
      );
    },
  };
}

describe('GenomeBrowserPanel navigation after removal', () => {
  beforeEach(() => {
    hub.mounts.length = 0;
    hub.retired.length = 0;
    resetTrackReadGuardForTests();
    sessionMock.mockReset();
    summaryMock.mockReset();
    sessionMock.mockResolvedValue(session());
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

  it('draws the coordinate typed immediately after a track is removed', async () => {
    const view = renderPanel();
    const box = await screen.findByRole('textbox', { name: 'Genomic region' });
    await waitFor(() => expect(box).toHaveValue('chr1:4500-6700'));
    fireEvent.click(screen.getByRole('checkbox', { name: /Signal: A/ }));
    fireEvent.click(screen.getByRole('checkbox', { name: /Signal: B/ }));
    expect(await screen.findByTestId('washu-mock-viewer')).toHaveTextContent('chr1:4500-6700');
    const mounted = hub.mounts.at(-1);
    expect(mounted?.storeId).toBe('helixweave-browser-run-1');

    fireEvent.click(screen.getByRole('checkbox', { name: /Signal: A/ }));
    expect(screen.getByTestId('washu-mock-viewer')).toHaveTextContent('chr1:4500-6700');
    fireEvent.change(box, { target: { value: 'chr1:5200-6400' } });

    await waitFor(() =>
      expect(screen.getByTestId('washu-mock-viewer')).toHaveTextContent('chr1:5200-6400'),
    );
    expect(box).toHaveValue('chr1:5200-6400');
    const drawn = hub.mounts.at(-1);
    expect(drawn).toEqual({
      region: 'chr1:5200-6400',
      storeId: 'helixweave-browser-run-1',
    });
    expect(hub.retired).toContainEqual({
      region: 'chr1:4500-6700',
      storeId: 'helixweave-browser-run-1',
    });
    expect(hub.mounts.filter((entry) => entry.storeId !== 'helixweave-browser-run-1')).toEqual([]);

    const mounts = hub.mounts.length;
    view.rerenderPanel();
    expect(hub.mounts).toHaveLength(mounts);
    expect(screen.getByTestId('washu-mock-viewer')).toHaveTextContent('chr1:5200-6400');
  });
});
