import { useEffect, useMemo, useState } from 'react';
import { useQuery } from '@tanstack/react-query';
import { GenomeHub } from 'wuepgg';
import washUStyles from 'wuepgg/style.css?raw';
import { catalogGetArchivedTrack } from '../../api/archivedTracks';
import LocusSearch from '../../viewer/search/LocusSearch';
import { ArchivedTrackFacts } from './ArchivedTrackFacts';
import {
  archivedContentPath,
  archivedHubTracks,
  archivedIndexPath,
  archivedTrackDataProps,
  archivedTrackScaffold,
  archiveErrorCopy,
  archiveErrorFromUnknown,
  isDrawableRegion,
  preflightArchivedBytes,
  sameOriginApiUrl,
  type ArchiveViewerRequest,
} from './archiveView';

export function ArchivedTrackViewer({ request }: { request: ArchiveViewerRequest }) {
  const ready = request.kind === 'ready' ? request : null;
  const identity = ready?.archiveVersionId ?? '';
  const [seenIdentity, setSeenIdentity] = useState(identity);
  const [region, setRegion] = useState('');
  if (seenIdentity !== identity) {
    setSeenIdentity(identity);
    setRegion('');
  }

  const descriptorQuery = useQuery({
    queryKey: ['archived-track', ready?.projectId, ready?.experimentId, ready?.archiveVersionId],
    enabled: ready !== null,
    retry: false,
    queryFn: ({ signal }) => {
      if (ready === null) throw new Error('Archive request is incomplete.');
      return catalogGetArchivedTrack(ready.projectId, ready.experimentId, ready.archiveVersionId, {
        signal,
      });
    },
  });

  const descriptor = descriptorQuery.data ?? null;
  const mismatch =
    descriptor !== null &&
    ready !== null &&
    (descriptor.archive_version_id !== ready.archiveVersionId ||
      descriptor.project_id !== ready.projectId ||
      descriptor.experiment_id !== ready.experimentId);
  const track = descriptor !== null && !mismatch ? descriptor : null;
  const scaffold = track === null ? null : archivedTrackScaffold(track);
  const contentPath = track === null ? null : archivedContentPath(track);
  const indexPath = track === null ? null : archivedIndexPath(track);
  const contentUrl = contentPath === null ? null : sameOriginApiUrl(contentPath);
  const indexUrl = indexPath === null ? null : sameOriginApiUrl(indexPath);
  const bytesReady =
    track !== null && scaffold !== null && contentUrl !== null && (indexPath === null || indexUrl !== null);

  const bytesQuery = useQuery({
    queryKey: ['archived-track-bytes', contentUrl, indexUrl],
    enabled: bytesReady,
    retry: false,
    queryFn: ({ signal }) => {
      if (contentUrl === null || track === null) throw new Error('Archive content URL is missing.');
      return preflightArchivedBytes(track, contentUrl, indexUrl, signal);
    },
  });

  useEffect(() => {
    if (identity === '') return undefined;
    const style = document.createElement('style');
    style.dataset.helixweaveWashu = `archive-${identity}`;
    style.textContent = washUStyles;
    document.head.appendChild(style);
    const storeId = `helixweave-archive-${identity}`;
    return () => {
      style.remove();
      localStorage.removeItem(`persist:${storeId}`);
    };
  }, [identity]);

  const storeConfig = useMemo(
    () => ({ storeId: `helixweave-archive-${identity}`, enablePersistence: false }),
    [identity],
  );
  const tracks = useMemo(
    () => (track !== null && contentUrl !== null ? archivedHubTracks(track, contentUrl, bytesQuery.data) : []),
    [contentUrl, track, bytesQuery.data],
  );

  let error = request.kind === 'invalid' ? archiveErrorCopy('archive_invalid') : null;
  if (descriptorQuery.isError) error = archiveErrorFromUnknown(descriptorQuery.error);
  else if (mismatch) error = archiveErrorCopy('archive_version_conflict');
  else if (track !== null && scaffold === null) {
    error = {
      code: 'archive_invalid',
      message: '这个组装没有可用的 WashU 染色体列表。轨道没有被改标。页面没有改读另一版。',
    };
  } else if (track !== null && contentUrl === null) {
    error = {
      code: 'archive_invalid',
      message: '封存轨道只能从当前站点读取。页面没有改读另一版。',
    };
  } else if (bytesQuery.isError) error = archiveErrorFromUnknown(bytesQuery.error);
  const draw =
    error === null && bytesQuery.isSuccess && scaffold !== null && (region === '' || isDrawableRegion(region));

  return (
    <article
      aria-label="已封存轨道"
      lang="zh-CN"
      className="min-w-0 space-y-3"
      data-testid="archived-track-viewer"
      data-content-path={contentPath ?? ''}
      data-index-path={indexPath ?? ''}
      {...(track !== null ? archivedTrackDataProps(track) : {})}
      {...(ready !== null && track === null
        ? {
            'data-project-id': ready.projectId,
            'data-experiment-id': ready.experimentId,
            'data-archive-version-id': ready.archiveVersionId,
          }
        : {})}
    >
      <h1 className="text-lg font-semibold">已封存轨道</h1>
      {ready !== null && track === null ? (
        <p data-testid="archived-track-request">封存版本 {ready.archiveVersionId}</p>
      ) : null}
      {descriptorQuery.isPending && ready !== null ? <p role="status">正在读取这份封存版本…</p> : null}
      {track !== null ? <ArchivedTrackFacts track={track} /> : null}
      {error !== null ? (
        <p role="alert" data-testid="archived-track-error" data-error-code={error.code}>
          {error.message}
        </p>
      ) : null}
      {bytesQuery.isPending && bytesReady ? <p role="status">正在核对封存副本…</p> : null}
      {track !== null && scaffold !== null && contentUrl !== null && error === null ? (
        <LocusSearch assembly={track.assembly} region={region} onNavigate={setRegion} />
      ) : null}
      {draw && scaffold !== null ? (
        <div className="min-h-[28rem] overflow-hidden" data-testid="washu-genome-viewer">
          <GenomeHub
            key={JSON.stringify([scaffold, contentUrl, region])}
            storeConfig={storeConfig}
            genomeName={scaffold}
            viewRegion={region || undefined}
            tracks={tracks}
            showGenomeNavigator={false}
            showNavBar={false}
            showToolBar={false}
            showDisclosure={false}
          />
        </div>
      ) : null}
    </article>
  );
}
