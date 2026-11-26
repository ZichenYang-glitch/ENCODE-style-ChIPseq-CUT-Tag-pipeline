import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest';
import { ApiError, fetcher } from './fetcher';
import {
  catalogGetArchivedTrack,
  catalogListArchivedTracks,
  parseArchivedTrack,
} from './archivedTracks';
import {
  ARCHIVE_V1,
  ARCHIVE_V2,
  EXPERIMENT_ID,
  PROJECT_ID,
  sealedTrack,
} from '../features/experiment-archive/archiveFixtures';

vi.mock('./fetcher', async () => {
  const actual = await vi.importActual<typeof import('./fetcher')>('./fetcher');
  return { ...actual, fetcher: vi.fn() };
});

const fetchMock = vi.mocked(fetcher);

function wire(track: ReturnType<typeof sealedTrack>, extra: Record<string, unknown> = {}) {
  return { ...track, source: { ...track.source }, ...extra };
}

describe('archived track client', () => {
  beforeEach(() => {
    fetchMock.mockReset();
  });
  afterEach(() => {
    vi.unstubAllEnvs();
  });

  it('rejects the latest alias and an off-site API before sending metadata requests', async () => {
    await expect(catalogGetArchivedTrack(PROJECT_ID, EXPERIMENT_ID, 'latest')).rejects.toMatchObject({ code: 'archive_invalid' });
    vi.stubEnv('VITE_API_BASE_URL', 'https://other.example');
    await expect(catalogGetArchivedTrack(PROJECT_ID, EXPERIMENT_ID, ARCHIVE_V1)).rejects.toMatchObject({ code: 'archive_invalid' });
    await expect(catalogListArchivedTracks(PROJECT_ID, EXPERIMENT_ID)).rejects.toMatchObject({ code: 'archive_invalid' });
    expect(fetchMock).not.toHaveBeenCalled();
  });

  it('lists every returned version in order and does not keep only the newest', async () => {
    const older = sealedTrack(1);
    const newer = sealedTrack(2);
    fetchMock.mockResolvedValue({ items: [wire(newer), wire(older)], next_offset: null });
    const page = await catalogListArchivedTracks(PROJECT_ID, EXPERIMENT_ID, { limit: 50, offset: 0 });
    expect(page.items.map((item) => item.archive_version_id)).toEqual([ARCHIVE_V2, ARCHIVE_V1]);
    expect(page.items.map((item) => item.version_number)).toEqual([2, 1]);
    expect(fetchMock).toHaveBeenCalledWith(
      `/api/v1/catalog/projects/${PROJECT_ID}/experiments/${EXPERIMENT_ID}/archived-tracks?limit=50&offset=0`,
      expect.objectContaining({ method: 'GET' }),
    );
  });

  it('reads the requested archive id and rejects a different version in the body', async () => {
    fetchMock.mockResolvedValue(wire(sealedTrack(2)));
    await expect(
      catalogGetArchivedTrack(PROJECT_ID, EXPERIMENT_ID, ARCHIVE_V1),
    ).rejects.toBeInstanceOf(ApiError);
    await expect(
      catalogGetArchivedTrack(PROJECT_ID, EXPERIMENT_ID, ARCHIVE_V1),
    ).rejects.toMatchObject({ code: 'archive_version_conflict' });
    expect(fetchMock).toHaveBeenCalledWith(
      `/api/v1/catalog/projects/${PROJECT_ID}/experiments/${EXPERIMENT_ID}/archived-tracks/${ARCHIVE_V1}`,
      expect.objectContaining({ method: 'GET' }),
    );
    expect(fetchMock.mock.calls.some((call) => String(call[0]).includes(ARCHIVE_V2))).toBe(false);
  });

  it('keeps valid experiment labels and track keys in the researcher’s language', () => {
    const original = sealedTrack(1, { display_name: '信号'.repeat(110), track_key: '实验信号一' });
    expect(parseArchivedTrack(wire(original), { projectId: PROJECT_ID, experimentId: EXPERIMENT_ID }))
      .toMatchObject({ display_name: original.display_name, track_key: original.track_key });
  });

  it('does not request a path-like id or a descriptor that carries a file path', async () => {
    await expect(
      catalogGetArchivedTrack(PROJECT_ID, EXPERIMENT_ID, '../secret'),
    ).rejects.toMatchObject({ code: 'archive_invalid' });
    expect(fetchMock).not.toHaveBeenCalled();
    const tainted = wire(sealedTrack(1), { path: '/archive/signal.bw' });
    expect(parseArchivedTrack(tainted, { projectId: PROJECT_ID, experimentId: EXPERIMENT_ID })).toBeNull();
  });
});
