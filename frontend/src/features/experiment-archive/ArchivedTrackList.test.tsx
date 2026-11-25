import { QueryClient, QueryClientProvider } from '@tanstack/react-query';
import { render, screen, within } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { beforeEach, describe, expect, it, vi } from 'vitest';
import { catalogGetArchivedTrack, catalogListArchivedTracks } from '../../api/archivedTracks';
import { ArchivedTrackList } from './ArchivedTrackList';
import {
  ARCHIVE_V1,
  ARCHIVE_V2,
  EXPERIMENT_ID,
  PROJECT_ID,
  sealedTrack,
} from './archiveFixtures';
import { archivedViewerHref } from './archiveView';

vi.mock('../../api/archivedTracks', () => ({
  catalogListArchivedTracks: vi.fn(),
  catalogGetArchivedTrack: vi.fn(),
}));

const listMock = vi.mocked(catalogListArchivedTracks);
const getMock = vi.mocked(catalogGetArchivedTrack);

function renderList() {
  const client = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={client}>
      <ArchivedTrackList projectId={PROJECT_ID} experimentId={EXPERIMENT_ID} />
    </QueryClientProvider>,
  );
}

describe('ArchivedTrackList', () => {
  beforeEach(() => {
    listMock.mockReset();
    getMock.mockReset();
  });

  it('keeps the first sealed version selectable when a later version has the same name', async () => {
    listMock.mockResolvedValue({
      items: [sealedTrack(2), sealedTrack(1)],
      next_offset: null,
    });
    renderList();
    const list = await screen.findByRole('region', { name: '已封存轨道' });
    const older = await within(list).findByRole('link', { name: '在 WashU 查看 signal.bw 第 1 版' });
    const newer = within(list).getByRole('link', { name: '在 WashU 查看 signal.bw 第 2 版' });
    expect(older).toHaveTextContent('在 WashU 查看第 1 版');
    expect(newer).toHaveTextContent('在 WashU 查看第 2 版');
    expect(older).toHaveAttribute('href', archivedViewerHref(sealedTrack(1)));
    expect(newer).toHaveAttribute('href', archivedViewerHref(sealedTrack(2)));
    expect(older.getAttribute('href')).not.toContain(ARCHIVE_V2);
    expect(newer.getAttribute('href')).toContain(ARCHIVE_V2);
    const rows = within(list).getAllByTestId('archived-track-version');
    expect(rows[0]).toHaveAttribute('data-version-number', '2');
    expect(rows[0]).toHaveAttribute('data-display-name', 'signal.bw');
    expect(rows[1]).toHaveAttribute('data-archive-version-id', ARCHIVE_V1);
    expect(rows[1]).toHaveAttribute('data-version-number', '1');
    expect(within(rows[1] as HTMLElement).getByTestId('archived-track-version-label')).toHaveTextContent(
      '第 1 版',
    );
    expect(within(rows[1] as HTMLElement).getByTestId('archived-track-source')).toHaveTextContent(
      '运行 run-v1',
    );
    expect(list.querySelector('input[type="file"]')).toBeNull();
    expect(within(list).queryByRole('radio')).not.toBeInTheDocument();
    expect(within(list).queryByRole('checkbox')).not.toBeInTheDocument();
    expect(getMock).not.toHaveBeenCalled();
    expect(listMock).toHaveBeenCalledWith(
      PROJECT_ID,
      EXPERIMENT_ID,
      { limit: 50, offset: 0 },
      expect.any(Object),
    );
  });

  it('requests the next offset from the response instead of a newer version', async () => {
    listMock.mockResolvedValueOnce({ items: [sealedTrack(2)], next_offset: 50 });
    renderList();
    expect(await screen.findByRole('link', { name: '在 WashU 查看 signal.bw 第 2 版' })).toBeInTheDocument();
    listMock.mockResolvedValueOnce({ items: [sealedTrack(1)], next_offset: null });
    const user = userEvent.setup();
    await user.click(screen.getByRole('button', { name: '下一页' }));
    expect(await screen.findByRole('link', { name: '在 WashU 查看 signal.bw 第 1 版' })).toHaveAttribute(
      'href',
      archivedViewerHref(sealedTrack(1)),
    );
    expect(listMock).toHaveBeenLastCalledWith(
      PROJECT_ID,
      EXPERIMENT_ID,
      { limit: 50, offset: 50 },
      expect.any(Object),
    );
  });

  it('keeps the older result reachable across a page filtered out by the API', async () => {
    listMock.mockResolvedValueOnce({ items: [], next_offset: 50 })
      .mockResolvedValueOnce({ items: [sealedTrack(1)], next_offset: null });
    renderList();
    expect(await screen.findByText('这一页没有可读取的封存版本。')).toBeInTheDocument();
    await userEvent.setup().click(screen.getByRole('button', { name: '下一页' }));
    expect(await screen.findByRole('link', { name: '在 WashU 查看 signal.bw 第 1 版' }))
      .toHaveAttribute('href', archivedViewerHref(sealedTrack(1)));
  });
});
