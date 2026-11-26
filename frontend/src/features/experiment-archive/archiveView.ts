import type { HTMLAttributes } from 'react';
import type { UiLanguage } from '../../lib/uiLanguage';
import { archiveText } from './archiveText';
import type { WashuTrackConfig } from 'wuepgg';
import { ApiError } from '../../api/fetcher';
import { ARCHIVE_VERSION_ID, EXPERIMENT_ID, PROJECT_ID } from '../../api/archiveScope';
import type { ArchivedTrack } from '../../api/archivedTracks';
import { getCatalogReadArchivedTrackContentUrl, getCatalogReadArchivedTrackIndexUrl } from '../../api/generated/catalog/catalog';
import { viewerScaffold } from '../run-artifacts/browserWindow';

const REGION = /^[A-Za-z0-9_.:-]{1,128}:[0-9]{1,10}-[0-9]{1,10}$/;

export const ARCHIVE_PAGE_LIMIT = 50;

export const ARCHIVE_ERROR_CODES = [
  'archive_version_not_found',
  'archive_version_conflict',
  'archive_offline',
  'archive_integrity',
  'archive_invalid',
  'archive_unauthorized',
] as const;

export type ArchiveErrorCode = (typeof ARCHIVE_ERROR_CODES)[number];

export type ArchiveViewerRequest =
  | { kind: 'absent' }
  | { kind: 'invalid' }
  | {
      kind: 'ready';
      projectId: string;
      experimentId: string;
      archiveVersionId: string;
    };

export function archiveErrorCopy(code: string | null | undefined, language: UiLanguage = 'en') {
  const known = ARCHIVE_ERROR_CODES.find((item) => item === code) ?? 'archive_invalid';
  const copy = archiveText[language];
  return { code: known, message: [copy.errors[known], copy.noFallback].join(language === 'en' ? ' ' : '') };
}

export function archiveErrorFromUnknown(error: unknown, language: UiLanguage) {
  if (error instanceof ApiError) return archiveHttpError(error.status, error.code, language);
  return archiveErrorCopy('archive_invalid', language);
}

export function archiveHttpError(status: number, code: string | null, language: UiLanguage = 'en') {
  if (status === 401 || status === 403) return archiveErrorCopy('archive_unauthorized', language);
  if (status === 404) return archiveErrorCopy(code ?? 'archive_version_not_found', language);
  if (status === 409) return archiveErrorCopy(code ?? 'archive_version_conflict', language);
  return archiveErrorCopy(code, language);
}

export function archiveViewerRequest(search: string): ArchiveViewerRequest {
  const params = new URLSearchParams(search.startsWith('?') ? search.slice(1) : search);
  const names = ['project', 'experiment', 'archive'] as const;
  if (names.every((name) => !params.has(name))) return { kind: 'absent' };
  if (names.some((name) => params.getAll(name).length !== 1)) return { kind: 'invalid' };
  const projectId = params.get('project') ?? '';
  const experimentId = params.get('experiment') ?? '';
  const archiveVersionId = params.get('archive') ?? '';
  if (
    !PROJECT_ID.test(projectId) ||
    !EXPERIMENT_ID.test(experimentId) ||
    !ARCHIVE_VERSION_ID.test(archiveVersionId)
  ) {
    return { kind: 'invalid' };
  }
  return { kind: 'ready', projectId, experimentId, archiveVersionId };
}

export function archivedViewerHref(
  track: Pick<ArchivedTrack, 'project_id' | 'experiment_id' | 'archive_version_id'>,
): string {
  const params = new URLSearchParams({
    project: track.project_id,
    experiment: track.experiment_id,
    archive: track.archive_version_id,
  });
  return `/viewer.html?${params.toString()}`;
}

export function archivedTrackOpenName(
  track: Pick<ArchivedTrack, 'display_name' | 'version_number'>,
  language: UiLanguage,
): string {
  return language === 'en'
    ? `Open version ${track.version_number} in WashU · ${track.display_name}`
    : `在 WashU 查看 ${track.display_name} 第 ${track.version_number} 版`;
}

export function archivedVersionText(versionNumber: number, language: UiLanguage): string {
  return language === 'en' ? `Version ${versionNumber}` : `第 ${versionNumber} 版`;
}

export function archivedSourceText(track: ArchivedTrack, language: UiLanguage): string {
  const labels = archiveText[language];
  if (track.source.kind === 'run') {
    const base = `${labels.run} ${track.source.run_id} · ${labels.artifact} ${track.source.artifact_id} · ${labels.generation} ${track.source.generation} · ${labels.revision} ${track.source.revision}`;
    if (
      track.source.index_artifact_id !== undefined &&
      track.source.index_generation !== undefined &&
      track.source.index_revision !== undefined
    ) {
      return `${base} · ${labels.indexArtifact} ${track.source.index_artifact_id} · ${labels.indexGeneration} ${track.source.index_generation} · ${labels.indexRevision} ${track.source.index_revision}`;
    }
    return base;
  }
  const base = `${labels.registered} ${track.source.result_id} · ${labels.revision} ${track.source.result_revision_id}`;
  if (
    track.source.index_result_id !== undefined &&
    track.source.index_result_revision_id !== undefined
  ) {
    return `${base} · ${labels.indexRegistered} ${track.source.index_result_id} · ${labels.indexRevision} ${track.source.index_result_revision_id}`;
  }
  return base;
}

export function archivedContentPath(track: ArchivedTrack): string {
  return getCatalogReadArchivedTrackContentUrl(track.project_id, track.experiment_id, track.archive_version_id);
}

export function archivedIndexPath(track: ArchivedTrack): string | null {
  if (!track.has_index) return null;
  return getCatalogReadArchivedTrackIndexUrl(track.project_id, track.experiment_id, track.archive_version_id);
}

export function archivedHubTracks(track: ArchivedTrack, contentUrl: string, textFile?: File | null): WashuTrackConfig[] {
  const data: WashuTrackConfig = {
    type: track.track_type,
    name: `${track.display_name} · v${track.version_number}`,
    url: contentUrl,
    waitToUpdate: true,
  };
  if (textFile) Object.assign(data, { fileObj: textFile, url: '', isText: true, textConfig: {} });
  if (track.track_type === 'longrange') data.options = { clampHeight: true };
  return [
    { type: 'ruler', name: 'Ruler', waitToUpdate: true },
    data,
  ];
}

export function archivedTrackScaffold(track: ArchivedTrack): string | null {
  return viewerScaffold(track.assembly);
}

export function isDrawableRegion(region: string): boolean {
  return REGION.test(region);
}

export function archivedTrackDataProps(track: ArchivedTrack): HTMLAttributes<HTMLElement> {
  const props: Record<string, string> = {
    'data-archive-version-id': track.archive_version_id,
    'data-version-number': String(track.version_number),
    'data-display-name': track.display_name,
    'data-track-key': track.track_key,
    'data-experiment-revision-id': track.experiment_revision_id,
    'data-assembly': track.assembly,
    'data-track-type': track.track_type,
    'data-has-index': track.has_index ? 'true' : 'false',
    'data-project-id': track.project_id,
    'data-experiment-id': track.experiment_id,
    'data-source-kind': track.source.kind,
  };
  if (track.source.kind === 'run') {
    props['data-run-id'] = track.source.run_id;
    props['data-artifact-id'] = track.source.artifact_id;
    props['data-generation'] = track.source.generation;
    props['data-revision'] = track.source.revision;
    props['data-index-artifact-id'] = track.source.index_artifact_id ?? '';
    props['data-index-generation'] = track.source.index_generation ?? '';
    props['data-index-revision'] = track.source.index_revision ?? '';
  } else {
    props['data-result-id'] = track.source.result_id;
    props['data-result-revision-id'] = track.source.result_revision_id;
    props['data-index-result-id'] = track.source.index_result_id ?? '';
    props['data-index-result-revision-id'] = track.source.index_result_revision_id ?? '';
  }
  return props as HTMLAttributes<HTMLElement>;
}
