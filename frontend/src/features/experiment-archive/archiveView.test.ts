import { preflightArchivedBytes } from './archivedBytes';
import { afterEach, describe, expect, it, vi } from 'vitest';
import {
  ARCHIVE_V1,
  ARCHIVE_V2,
  EXPERIMENT_ID,
  PROJECT_ID,
  sealedTrack,
} from './archiveFixtures';
import {
  archivedContentPath,
  archivedHubTracks,
  archivedIndexPath,
  archivedSourceText,
  archivedTrackOpenName,
  archivedViewerHref,
  archiveViewerRequest,
} from './archiveView';

describe('archive viewer address', () => {
  it('names one sealed version in the WashU address', () => {
    const older = sealedTrack(1);
    const href = archivedViewerHref(older);
    const request = archiveViewerRequest(href.slice(href.indexOf('?')));
    expect(archivedTrackOpenName(older, 'zh-CN')).toBe('在 WashU 查看 signal.bw 第 1 版');
    expect(archivedTrackOpenName(older, 'en')).toBe('Open version 1 of signal.bw in WashU');
    expect(request).toEqual({
      kind: 'ready',
      projectId: PROJECT_ID,
      experimentId: EXPERIMENT_ID,
      archiveVersionId: ARCHIVE_V1,
    });
    expect(href).not.toContain(ARCHIVE_V2);
    expect(href).not.toContain('latest');
  });

  it('keeps a partial or duplicated archive address from selecting any version', () => {
    expect(archiveViewerRequest('')).toEqual({ kind: 'absent' });
    expect(archiveViewerRequest(`?run=run-1`)).toEqual({ kind: 'absent' });
    expect(archiveViewerRequest(`?archive=${ARCHIVE_V1}`)).toEqual({ kind: 'invalid' });
    expect(
      archiveViewerRequest(
        `?project=${PROJECT_ID}&project=${PROJECT_ID}&experiment=${EXPERIMENT_ID}&archive=${ARCHIVE_V1}`,
      ),
    ).toEqual({ kind: 'invalid' });
    expect(
      archiveViewerRequest(`?project=${PROJECT_ID}&experiment=${EXPERIMENT_ID}&archive=../secret`),
    ).toEqual({ kind: 'invalid' });
  });

  it('builds same-origin content and index paths from the exact version', () => {
    const indexed = sealedTrack(1, { has_index: true, track_type: 'longrange' });
    expect(archivedContentPath(indexed)).toBe(
      `/api/v1/catalog/projects/${PROJECT_ID}/experiments/${EXPERIMENT_ID}/archived-tracks/${ARCHIVE_V1}/content`,
    );
    expect(archivedIndexPath(indexed)).toBe(`${archivedContentPath(indexed)}.tbi`);
    expect(archivedIndexPath(sealedTrack(2))).toBeNull();
    const [ruler, data] = archivedHubTracks(indexed, 'https://lab.example/content');
    expect(ruler).toMatchObject({ type: 'ruler' });
    expect(data).toMatchObject({
      name: 'signal.bw · v1',
      type: 'longrange',
      url: 'https://lab.example/content',
      options: { clampHeight: true },
    });
    expect(String(data?.url)).not.toContain(ARCHIVE_V2);
  });

  it('shows run and registered identities without a file path', () => {
    expect(archivedSourceText(sealedTrack(1), 'zh-CN')).toBe(
      `运行 run-v1 · 产物 signal · 代次 artifactgen-${'a'.repeat(64)} · 修订 artifactrev-${'b'.repeat(64)}`,
    );
    const registered = sealedTrack(1, {
      source: {
        kind: 'registered',
        result_id: `rreg_${'c'.repeat(32)}`,
        result_revision_id: `rrev_${'d'.repeat(32)}`,
        index_result_id: `rreg_${'e'.repeat(32)}`,
        index_result_revision_id: `rrev_${'f'.repeat(32)}`,
      },
    });
    expect(archivedSourceText(registered, 'zh-CN')).toContain(`登记 rreg_${'c'.repeat(32)}`);
    expect(archivedSourceText(registered, 'zh-CN')).toContain(`索引修订 rrev_${'f'.repeat(32)}`);
    expect(archivedSourceText(registered, 'zh-CN')).not.toContain('/');
  });
});


describe('sealed plain BED reader', () => {
  afterEach(() => {
    vi.unstubAllGlobals();
  });

  function serve(text: string, advertisedSize?: number) {
    const bytes = new TextEncoder().encode(text);
    const total = advertisedSize ?? bytes.length;
    vi.stubGlobal('fetch', vi.fn().mockResolvedValue(new Response(bytes, {
      status: 206, headers: { 'Content-Range': `bytes 0-${bytes.length - 1}/${total}` },
    })));
  }

  const peak = sealedTrack(1, { track_type: 'bed', display_name: 'signal.bed' });
  const url = archivedContentPath(peak);
  const open = () => preflightArchivedBytes(peak, url, null, new AbortController().signal);

  it('gives the native text reader the intervals from the selected version', async () => {
    serve('track db=hg38\nchr7\t100\t150\tpeak1\n');
    const file = await open();
    expect(file).toBeInstanceOf(File);
    const [, data] = archivedHubTracks(peak, url, file);
    expect(data).toMatchObject({ type: 'bed', url: '', fileObj: file, isText: true });
    expect(fetch).toHaveBeenCalledTimes(1);
    expect(fetch).toHaveBeenCalledWith(url, expect.objectContaining({ credentials: 'same-origin' }));
  });

  it('rejects a different declared assembly and incomplete transfers', async () => {
    serve('track db=mm10\nchr7\t100\t150\n');
    await expect(open()).rejects.toMatchObject({ code: 'archive_version_conflict' });
    serve('chr7\t100\t150\n', 16 * 1024 * 1024);
    await expect(open()).rejects.toMatchObject({ code: 'archive_invalid' });
  });
});
