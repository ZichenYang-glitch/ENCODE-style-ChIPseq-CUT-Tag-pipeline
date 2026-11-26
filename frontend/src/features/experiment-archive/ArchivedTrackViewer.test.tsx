import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { act, render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { runInNewContext } from 'node:vm';
import { setUiLanguage } from '../../lib/uiLanguage';
import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { catalogGetArchivedTrack, catalogListArchivedTracks } from '../../api/archivedTracks';
import { ApiError } from '../../api/fetcher';
import { ArchivedTrackViewer } from './ArchivedTrackViewer';
import {
  ARCHIVE_V1,
  ARCHIVE_V2,
  EXPERIMENT_ID,
  PROJECT_ID,
  sealedTrack,
} from './archiveFixtures';
import type { ArchiveViewerRequest } from './archiveView';

vi.mock('wuepgg', () => ({
  GenomeHub: ({
    genomeName,
    viewRegion,
    tracks,
  }: {
    genomeName: string;
    viewRegion?: string;
    tracks: Array<{ type?: string; url?: string; options?: Record<string, unknown> }>;
  }) => (
    <div
      data-testid="washu-mock-viewer"
      data-genome={genomeName}
      data-region={viewRegion ?? ''}
      data-track-url={tracks.find((track) => track.type !== 'ruler')?.url ?? ''}
      data-longrange-options={JSON.stringify(
        tracks.find((track) => track.type === 'longrange')?.options ?? null,
      )}
    />
  ),
}));

vi.mock('../../api/archivedTracks', () => ({
  catalogGetArchivedTrack: vi.fn(),
  catalogListArchivedTracks: vi.fn(),
}));

const getMock = vi.mocked(catalogGetArchivedTrack);
const listMock = vi.mocked(catalogListArchivedTracks);

const ready: ArchiveViewerRequest = {
  kind: 'ready',
  projectId: PROJECT_ID,
  experimentId: EXPERIMENT_ID,
  archiveVersionId: ARCHIVE_V1,
};

function renderViewer(request: ArchiveViewerRequest = ready) {
  const client = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={client}>
      <ArchivedTrackViewer request={request} />
    </QueryClientProvider>,
  );
}

describe('ArchivedTrackViewer', () => {
  afterEach(() => localStorage.removeItem('helixweave.ui-language'));

  beforeEach(() => {
    setUiLanguage('zh-CN');
    getMock.mockReset();
    listMock.mockReset();
    vi.unstubAllGlobals();
  });

  it('draws the requested version from its same-origin content URL', async () => {
    getMock.mockResolvedValue(sealedTrack(1));
    const seen: Array<{ url: string; range: string | null }> = [];
    vi.stubGlobal(
      'fetch',
      vi.fn(async (input: RequestInfo | URL, init?: RequestInit) => {
        const url = typeof input === 'string' ? input : input instanceof URL ? input.toString() : input.url;
        seen.push({ url, range: new Headers(init?.headers).get('range') });
        return new Response(new Uint8Array([4]), { status: 206 });
      }),
    );
    const user = userEvent.setup();
    renderViewer();
    expect(await screen.findByTestId('archived-track-version-label')).toHaveTextContent('第 1 版');
    expect(screen.getByTestId('archived-track-label')).toHaveTextContent('signal.bw');
    expect(screen.getByTestId('archived-track-source')).toHaveTextContent('运行 run-v1');
    expect(screen.getByTestId('archived-track-assembly')).toHaveTextContent('hg38');
    expect(screen.getByTestId('archived-track-viewer')).toHaveAttribute('data-archive-version-id', ARCHIVE_V1);
    expect(screen.getByTestId('archived-track-viewer')).toHaveAttribute(
      'data-content-path',
      `/api/v1/catalog/projects/${PROJECT_ID}/experiments/${EXPERIMENT_ID}/archived-tracks/${ARCHIVE_V1}/content`,
    );
    const initial = await screen.findByTestId('washu-mock-viewer');
    expect(initial).toHaveAttribute('data-track-url', expect.stringContaining(ARCHIVE_V1));
    expect(initial).toHaveAttribute('data-region', 'chr1:1-10000');
    const locus = screen.getByRole('textbox', { name: '基因组位置' });
    expect(locus).toHaveValue('chr1:1-10000');
    await user.clear(locus);
    await user.type(locus, 'chr1:100-200');
    const hub = await screen.findByTestId('washu-genome-viewer');
    expect(hub).toBeInTheDocument();
    expect(screen.getByTestId('washu-mock-viewer')).toHaveAttribute('data-genome', 'hg38');
    expect(screen.getByTestId('washu-mock-viewer')).toHaveAttribute('data-region', 'chr1:100-200');
    expect(screen.getByTestId('washu-mock-viewer').getAttribute('data-track-url')).toContain(
      `/archived-tracks/${ARCHIVE_V1}/content`,
    );
    expect(screen.getByTestId('washu-mock-viewer').getAttribute('data-track-url')).not.toContain(ARCHIVE_V2);
    expect(seen.map((item) => item.url).join(' ')).toContain(ARCHIVE_V1);
    expect(seen.map((item) => item.url).join(' ')).not.toContain(ARCHIVE_V2);
    expect(document.querySelector('input[type="file"]')).toBeNull();
    expect(screen.queryByTestId('local-track-picker')).not.toBeInTheDocument();
    const nativeViewer = screen.getByTestId('washu-mock-viewer');
    const readsBeforeSwitch = seen.length;
    await user.selectOptions(screen.getByRole('combobox', { name: '语言' }), 'en');
    expect(screen.getByRole('heading', { name: 'Archived track' })).toBeVisible();
    expect(screen.getByTestId('archived-track-version-label')).toHaveTextContent('Version 1');
    expect(screen.getByRole('textbox', { name: 'Genomic region' })).toHaveValue('chr1:100-200');
    expect(screen.getByTestId('washu-mock-viewer')).toBe(nativeViewer);
    expect(seen).toHaveLength(readsBeforeSwitch);
    expect(listMock).not.toHaveBeenCalled();
    expect(getMock).toHaveBeenCalledWith(PROJECT_ID, EXPERIMENT_ID, ARCHIVE_V1, expect.any(Object));
  });

  it.each([
    [403, 'FORBIDDEN', 'archive_unauthorized'],
    [409, 'archive_integrity', 'archive_integrity'],
    [503, 'archive_offline', 'archive_offline'],
  ] as const)('removes the plot when a later WashU range returns %s', async (status, serverCode, expectedCode) => {
    getMock.mockResolvedValue(sealedTrack(1));
    vi.stubGlobal('fetch', vi.fn(async () => new Response(new Uint8Array([4]), { status: 206 })));
    const transport = vi.fn(async () =>
      new Response(
        JSON.stringify({ issues: [{ code: serverCode, message: 'private detail' }] }),
        { status, headers: { 'Content-Type': 'application/json' } },
      ),
    );
    const createUrl = vi.fn((_blob: Blob | MediaSource) => 'blob:worker');
    vi.stubGlobal('URL', class extends URL { static createObjectURL = createUrl; });
    vi.stubGlobal('Worker', class {});
    const receiver = { onmessage: (_event: { data: unknown }) => {}, close: vi.fn() };
    vi.stubGlobal('BroadcastChannel', vi.fn(() => receiver));
    const view = renderViewer();
    const hub = await screen.findByTestId('washu-mock-viewer');
    const url = hub.getAttribute('data-track-url')!;
    URL.createObjectURL(new Blob(['/* WashU worker */'], { type: 'text/javascript;charset=utf-8' }));
    const script = await new Promise<string>((resolve, reject) => {
      const reader = new FileReader();
      reader.onload = () => resolve(String(reader.result));
      reader.onerror = reject;
      reader.readAsText(createUrl.mock.calls[0][0] as Blob);
    });
    const worker: { fetch: typeof fetch } = { fetch: transport };
    runInNewContext(script, {
      self: worker, URL, Request,
      BroadcastChannel: class { postMessage(data: unknown) { receiver.onmessage({ data }); } },
    });
    await act(async () => {
      await worker.fetch(url, { headers: { Range: 'bytes=2048-2111' } });
    });
    const alert = await screen.findByRole('alert');
    expect(alert).toHaveAttribute('data-error-code', expectedCode);
    expect(alert).not.toHaveTextContent('private detail');
    expect(screen.queryByTestId('washu-genome-viewer')).toBeNull();
    if (status === 403) expect(screen.queryByTestId('archived-track-source')).toBeNull();
    await expect(worker.fetch('https://proxy.example/track')).rejects.toThrow();
    expect(transport).toHaveBeenCalledTimes(1);
    expect(transport).toHaveBeenCalledWith(url, expect.objectContaining({ mode: 'same-origin', redirect: 'error' }));
    expect(listMock).not.toHaveBeenCalled();
    view.unmount();
    expect(URL.createObjectURL).toBe(createUrl);
    expect(receiver.close).toHaveBeenCalled();
  });

  it('does not draw when the response names a different version', async () => {
    getMock.mockResolvedValue(sealedTrack(2));
    renderViewer();
    expect(await screen.findByRole('alert')).toHaveAttribute('data-error-code', 'archive_version_conflict');
    expect(screen.queryByTestId('washu-genome-viewer')).not.toBeInTheDocument();
    expect(screen.getByTestId('archived-track-request')).toHaveTextContent(ARCHIVE_V1);
    expect(listMock).not.toHaveBeenCalled();
  });

  it('shows a known offline error without the server detail or another version', async () => {
    getMock.mockRejectedValue(new ApiError(503, 'archive_offline', 'missing /mnt/archive/signal.bw'));
    renderViewer();
    const alert = await screen.findByRole('alert');
    expect(alert).toHaveAttribute('data-error-code', 'archive_offline');
    expect(alert).toHaveTextContent('封存副本离线');
    expect(alert).toHaveTextContent('没有改读另一版');
    await userEvent.setup().selectOptions(screen.getByRole('combobox', { name: '语言' }), 'en');
    expect(alert).toHaveTextContent('The archived copy is offline. No other version was opened.');
    expect(alert).toHaveAttribute('data-error-code', 'archive_offline');
    expect(alert).not.toHaveTextContent('/mnt/archive');
    expect(getMock).toHaveBeenCalledTimes(1);
    expect(screen.queryByTestId('washu-genome-viewer')).not.toBeInTheDocument();
    expect(document.querySelector('input[type="file"]')).toBeNull();
  });

  it('checks the paired index and stops when that read fails', async () => {
    getMock.mockResolvedValue(sealedTrack(1, { has_index: true, track_type: 'longrange' }));
    const seen: Array<{ url: string; headers: RequestInit['headers'] | null }> = [];
    vi.stubGlobal(
      'fetch',
      vi.fn(async (input: RequestInfo | URL, init?: RequestInit) => {
        const url = typeof input === 'string' ? input : input instanceof URL ? input.toString() : input.url;
        seen.push({ url, headers: init?.headers ?? null });
        if (url.endsWith('.tbi')) {
          return new Response(
            JSON.stringify({ issues: [{ code: 'archive_integrity', message: 'sha changed' }] }),
            { status: 409, headers: { 'Content-Type': 'application/json' } },
          );
        }
        return new Response(new Uint8Array([1]), { status: 206 });
      }),
    );
    renderViewer();
    expect(await screen.findByRole('alert')).toHaveAttribute('data-error-code', 'archive_integrity');
    expect(seen.map((item) => item.url)).toEqual([
      expect.stringContaining(`/archived-tracks/${ARCHIVE_V1}/content`),
      expect.stringContaining(`/archived-tracks/${ARCHIVE_V1}/content.tbi`),
    ]);
    expect(seen[0]?.url.endsWith('.tbi')).toBe(false);
    expect(seen.map((item) => item.headers)).toEqual([
      { Range: 'bytes=0-0' },
      { Range: 'bytes=0-1' },
    ]);
    expect(screen.queryByTestId('washu-genome-viewer')).not.toBeInTheDocument();
    expect(screen.getByTestId('archived-track-version-label')).toHaveTextContent('第 1 版');
  });

  it('does not fetch when the address is not an exact archive version', async () => {
    renderViewer({ kind: 'invalid' });
    expect(await screen.findByRole('alert')).toHaveAttribute('data-error-code', 'archive_invalid');
    expect(getMock).not.toHaveBeenCalled();
    expect(listMock).not.toHaveBeenCalled();
    expect(document.querySelector('input[type="file"]')).toBeNull();
  });
});
