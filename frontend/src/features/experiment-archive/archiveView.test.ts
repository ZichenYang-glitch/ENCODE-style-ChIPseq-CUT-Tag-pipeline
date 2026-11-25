import { describe, expect, it } from 'vitest';
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
    expect(archivedTrackOpenName(older)).toBe('在 WashU 查看 signal.bw 第 1 版');
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
      type: 'longrange',
      url: 'https://lab.example/content',
      options: { clampHeight: true },
    });
    expect(String(data?.url)).not.toContain(ARCHIVE_V2);
  });

  it('shows run and registered identities without a file path', () => {
    expect(archivedSourceText(sealedTrack(1))).toBe(
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
    expect(archivedSourceText(registered)).toContain(`登记 rreg_${'c'.repeat(32)}`);
    expect(archivedSourceText(registered)).toContain(`索引修订 rrev_${'f'.repeat(32)}`);
    expect(archivedSourceText(registered)).not.toContain('/');
  });
});
