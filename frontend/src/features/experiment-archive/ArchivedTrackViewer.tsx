import { useEffect, useLayoutEffect, useMemo, useState } from 'react';
import { useQuery } from '@tanstack/react-query';
import { GenomeHub } from 'wuepgg';
import washUStyles from 'wuepgg/style.css?raw';
import { catalogGetArchivedTrack } from '../../api/archivedTracks';
import { sameOriginApiUrl } from '../../api/archiveScope';
import LocusSearch from '../../viewer/search/LocusSearch';
import { UiLanguageSelect, useUiLanguage } from '../../lib/uiLanguage';
import { archiveText } from './archiveText';
import { ArchivedTrackFacts } from './ArchivedTrackFacts';
import { preflightArchivedBytes } from './archivedBytes';
import { observeArchiveWorkerReads } from './archiveWorkerReads';
import {
  archivedContentPath,
  archivedHubTracks,
  archivedIndexPath,
  archivedTrackDataProps,
  archivedTrackScaffold,
  archiveErrorCopy,
  archiveErrorFromUnknown,
  isDrawableRegion,
  type ArchiveErrorCode,
  type ArchiveViewerRequest,
} from './archiveView';

// Every supported scaffold contains chr1. GenomeHub 62.0.2 needs an explicit
// interval to issue its first drawing read when a sealed track has no suggested locus.
const INITIAL_REGION = 'chr1:1-10000';

export function ArchivedTrackViewer({ request }: { request: ArchiveViewerRequest }) {
  const language = useUiLanguage();
  const text = archiveText[language];
  const ready = request.kind === 'ready' ? request : null;
  const identity = ready?.archiveVersionId ?? '';
  const [seenIdentity, setSeenIdentity] = useState(identity);
  const [region, setRegion] = useState(INITIAL_REGION);
  const [readError, setReadError] = useState<ArchiveErrorCode | null>(null);
  const [workerScope, setWorkerScope] = useState<string | null>(null);
  if (seenIdentity !== identity) {
    setSeenIdentity(identity);
    setRegion(INITIAL_REGION);
    setReadError(null);
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

  const descriptor = descriptorQuery.isSuccess ? descriptorQuery.data : null;
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
    gcTime: 0,
    queryFn: ({ signal }) => {
      if (contentUrl === null || track === null) throw new Error('Archive content URL is missing.');
      return preflightArchivedBytes(track, contentUrl, indexUrl, signal);
    },
  });

  useLayoutEffect(() => {
    if (contentUrl === null) return;
    const stop = observeArchiveWorkerReads(contentUrl, indexUrl, (failure) => {
      setReadError((current) => failure.code === 'archive_unauthorized' ? failure.code : current ?? failure.code);
    });
    setWorkerScope(contentUrl);
    return stop;
  }, [contentUrl, indexUrl]);

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

  let error = request.kind === 'invalid' ? archiveErrorCopy('archive_invalid', language) : null;
  if (descriptorQuery.isError) error = archiveErrorFromUnknown(descriptorQuery.error, language);
  else if (mismatch) error = archiveErrorCopy('archive_version_conflict', language);
  else if (track !== null && scaffold === null) {
    error = {
      code: 'archive_invalid',
      message: [text.unsupported, text.noFallback].join(language === 'en' ? ' ' : ''),
    };
  } else if (track !== null && contentUrl === null) {
    error = {
      code: 'archive_invalid',
      message: [text.sameOrigin, text.noFallback].join(language === 'en' ? ' ' : ''),
    };
  } else if (bytesQuery.isError) error = archiveErrorFromUnknown(bytesQuery.error, language);
  else if (readError !== null) error = archiveErrorCopy(readError, language);
  const draw =
    error === null && bytesQuery.isSuccess && scaffold !== null &&
    workerScope === contentUrl && isDrawableRegion(region);
  const visibleTrack = error?.code === 'archive_unauthorized' ? null : track;

  return (
    <article
      aria-label={text.viewer}
      lang={language}
      className="min-w-0 space-y-3"
      data-testid="archived-track-viewer"
      data-content-path={contentPath ?? ''}
      data-index-path={indexPath ?? ''}
      {...(visibleTrack !== null ? archivedTrackDataProps(visibleTrack) : {})}
      {...(ready !== null && visibleTrack === null
        ? {
            'data-project-id': ready.projectId,
            'data-experiment-id': ready.experimentId,
            'data-archive-version-id': ready.archiveVersionId,
          }
        : {})}
    >
      <header className="flex flex-wrap items-center justify-between gap-3">
        <h1 className="text-lg font-semibold">{text.viewer}</h1>
        <UiLanguageSelect />
      </header>
      {ready !== null && visibleTrack === null ? (
        <p data-testid="archived-track-request">{text.identifier} {ready.archiveVersionId}</p>
      ) : null}
      {descriptorQuery.isPending && ready !== null ? <p role="status">{text.loadingTrack}</p> : null}
      {visibleTrack !== null ? <ArchivedTrackFacts track={visibleTrack} language={language} /> : null}
      {error !== null ? (
        <p role="alert" data-testid="archived-track-error" data-error-code={error.code}>
          {error.message}
        </p>
      ) : null}
      {bytesQuery.isPending && bytesReady ? <p role="status">{text.verifying}</p> : null}
      {track !== null && scaffold !== null && contentUrl !== null && error === null ? (
        <LocusSearch assembly={track.assembly} region={region} onNavigate={setRegion} />
      ) : null}
      {draw && scaffold !== null ? (
        <div className="min-h-[28rem] overflow-hidden" data-testid="washu-genome-viewer">
          <GenomeHub
            key={JSON.stringify([scaffold, contentUrl, region])}
            storeConfig={storeConfig}
            genomeName={scaffold}
            viewRegion={region}
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
